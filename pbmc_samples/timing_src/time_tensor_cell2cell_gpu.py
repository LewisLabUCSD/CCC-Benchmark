#!/usr/bin/env python
# coding: utf-8
import os
import tensorly as tl
tl.set_backend('pytorch')

import pandas as pd
import cell2cell as c2c
from tqdm.auto import tqdm

import sys

# **ARGS**
root = str(sys.argv[1])
output = str(sys.argv[2])
n_samples = int(str(sys.argv[3]))
context_aggregation = str(sys.argv[4])

if (context_aggregation.upper() == 'FALSE') or (context_aggregation.upper() == 'F'):
    context_aggregation = False
    context_name = 'IND'
elif (context_aggregation.upper() == 'TRUE') or (context_aggregation.upper() == 'TRUE'):
    context_aggregation = True
    context_name = 'AGG'
else:
    raise ValueError("Specify whether aggregating or not by contexts using TRUE or FALSE. Used variable: {}".format(context_aggregation))

# **INPUTS**
#data_folder = '/data2/hratch/immune_CCI/covid/covid_atlas/interim/timing_inputs/'
data_folder = root + '/data/'

# **OUTPUTS**
output_folder = output + '/pbmc-outputs/'
if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

# # Data
samples = pd.read_csv(data_folder + '/PBMCs/samples_for_timing.csv', index_col=0)
sample_opts = samples.n_samples.unique().tolist()
assert n_samples in sample_opts, "Please select a number of samples in {}".format(sample_opts)
included_samples = samples.loc[samples['n_samples'] == n_samples]['sample_names'].values[0]
included_samples = included_samples.split('; ')
metadata = pd.read_csv(data_folder + '/PBMCs/metadata_for_timing.csv', index_col=0)

datadict = dict()
metadict = dict()
for s in tqdm(included_samples):
    datadict[s] = pd.read_hdf(data_folder + '/PBMCs/umi_per_sample.h5', key=s)
    metadict[s] = metadata.loc[metadata['sampleID'] == s]


# # Preprocessing

# **Metadata**
meta = metadata[['sampleID', 'CoVID-19 severity']]
meta.columns = ['sample', 'severity']
meta = meta.sort_values(['severity', 'sample'])

sample_meta = meta[['sample', 'severity']].drop_duplicates().reset_index(drop=True)
sample_meta = sample_meta.loc[sample_meta['sample'].isin(included_samples)]

sample_disease = dict()
for idx, row in sample_meta.iterrows():
    sample_disease[row['severity']] = row['severity']


# **Ligand-Receptor Pairs**
lr_pairs = pd.read_csv(data_folder + '/Human-2020-Jin-LR-pairs.csv')

# Change complex annotations
lr_pairs['ligand2'] = lr_pairs.interaction_name_2.apply(lambda x: x.split(' - ')[0].upper())
lr_pairs['receptor2'] = lr_pairs.interaction_name_2.apply(lambda x: x.split(' - ')[1].upper().replace('(', '').replace(')', '').replace('+', '&'))

lr_pairs['c2c_interaction'] = lr_pairs.apply(lambda row: row['ligand2'] + '^' + row['receptor2'], axis=1)

lr_pairs = c2c.preprocessing.ppi.remove_ppi_bidirectionality(lr_pairs, ('ligand2', 'receptor2'))

ppi_functions = dict()
for idx, row in lr_pairs.iterrows():
    ppi_label = row['ligand2'] + '^' + row['receptor2']
    ppi_functions[ppi_label] = row['annotation']


# **RNA-seq**
rnaseq_matrices = []

# Use sample_meta to preserve the new order of samples
if context_aggregation:
    for k, df in tqdm(sample_meta.groupby('severity'), total=len(sample_meta.groupby('severity'))):
        context_samples = df['sample'].values.tolist()
        meta_df = pd.concat([metadict[s] for s in context_samples])
        meta_df.index.name = 'Cell'
        meta_df = meta_df.reset_index()
        cells = meta_df['Cell'].values.tolist()
        data_df = pd.concat([datadict[s] for s in context_samples]).loc[cells,:]
        
        rnaseq_matrices.append(c2c.preprocessing.aggregate_single_cells(rnaseq_data=data_df.fillna(0),
                                                                        metadata=meta_df,
                                                                        barcode_col='Cell',
                                                                        celltype_col='majorType',
                                                                        method='nn_cell_fraction',
                                                                        transposed=True))
    context_labels = [k for k, df in sample_meta.groupby('severity')]
else:
    for k, row in tqdm(sample_meta.set_index('sample').iterrows(), total=len(sample_meta)):
        meta_df = metadict[k].copy()
        meta_df.index.name = 'Cell'
        meta_df = meta_df.reset_index()
        cells = meta_df['Cell'].values.tolist()
        data_df = datadict[k].loc[cells,:]

        rnaseq_matrices.append(c2c.preprocessing.aggregate_single_cells(rnaseq_data=data_df.fillna(0),
                                                                        metadata=meta_df,
                                                                        barcode_col='Cell',
                                                                        celltype_col='majorType',
                                                                        method='nn_cell_fraction',
                                                                        transposed=True))
    context_labels = sample_meta['severity'].tolist()


# # Run Analysis

# ### Tensor Factorization

# **Build 4D-Communication Tensor**
tensor = c2c.tensor.InteractionTensor(rnaseq_matrices=rnaseq_matrices,
                                      ppi_data=lr_pairs,
                                      context_names=context_labels,
                                      how='inner',
                                      complex_sep='&',
                                      interaction_columns=('ligand2', 'receptor2'),
                                      communication_score='expression_mean',
                                      device='cuda:0'
                                     )

# **Metadata for TF-Plot**
if context_aggregation:
    metadata_dicts = [None, None, None, None]
else:
    metadata_dicts = [sample_disease, None, None, None]

meta_tf = c2c.tensor.generate_tensor_metadata(interaction_tensor=tensor,
                                              metadata_dicts=metadata_dicts,
                                              fill_with_order_elements=True
                                             )


# **Elbow Analysis**
fig, error = tensor.elbow_rank_selection(upper_rank=25,
                                         runs=10,
                                         init='random',
                                         automatic_elbow=True,
                                         filename=output_folder + '/GPU-COVID-19-Elbow-{}Samples-{}.png'.format(n_samples, context_name),
                                         random_state=888)

# **Perform tensor factorization**
tensor.compute_tensor_factorization(rank=tensor.rank,
                                    init='random',
                                    random_state=888)

# **Plot factors**
cmaps = ['plasma', 'Dark2_r', 'tab20', 'tab20']
fig, axes = c2c.plotting.tensor_factors_plot(interaction_tensor=tensor,
                                             order_labels=['Samples', 'Ligand-Receptor Pairs', 'Sender Cells', 'Receiver Cells'],
                                             metadata = meta_tf,
                                             sample_col='Element',
                                             group_col='Category',
                                             meta_cmaps=cmaps,
                                             fontsize=14,
                                             filename=output_folder + '/GPU-COVID-19-TensorFactorization-{}Samples-{}.png'.format(n_samples, context_name)
                                            )

# **Export Loadings**
tensor.export_factor_loadings(output_folder + '/GPU-COVID-19-Loadings-{}Samples-{}.xlsx'.format(n_samples, context_name))


