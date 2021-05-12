#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import cell2cell as c2c
from tqdm.auto import tqdm
import os

import sys

# **ARGS**
root = str(sys.argv[1])
output = str(sys.argv[2])

# **INPUTS**
data_folder = root + '/data/'

# **OUTPUTS**
output_folder = output + '/outputs/'
if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

# # Data
directory = os.fsencode(data_folder + '/External/')

data = dict()
metadata = dict()

files = os.listdir(directory)
for file in tqdm(files, total=len(files)):
    filename = os.fsdecode(file)
    if filename.startswith("DGE_"): 
        print(filename)
        basename = os.path.basename(filename)
        sample = basename.split('_')[1]
        data[sample] = pd.read_csv(data_folder + '/External/' + filename, index_col=0)
        m = pd.read_csv(data_folder + '/External/' + filename.replace('DGE_', 'Meta_'), index_col=0)
        metadata[sample] = m
    else:
        continue





# # Preprocessing

# **Metadata**

meta = pd.concat(list(metadata.values()))
meta = meta.sort_values(['severity', 'sample'])
sample_meta = meta[['sample', 'severity']].drop_duplicates().reset_index(drop=True)
sample_meta = sample_meta.loc[sample_meta['sample'].isin(list(data.keys()))]

def meta_disease(x):
    if 'HC' in x:
        return 'Control'
    elif 'M' in x:
        return 'Moderate COVID-19'
    elif 'S' in x:
        return 'Severe COVID-19'
    else:
        return 'NA'
    
sample_disease = dict()
for idx, row in sample_meta.iterrows():
    sample_disease[row['severity']] = meta_disease(row['severity'])

# **Ligand-Receptor Pairs**
lr_pairs = pd.read_csv(data_folder + '/Human-2020-Jin-LR-pairs.csv')

# Change complex annotations
lr_pairs['ligand2'] = lr_pairs.interaction_name_2.apply(lambda x: x.split(' - ')[0].upper())
lr_pairs['receptor2'] = lr_pairs.interaction_name_2.apply(lambda x: x.split(' - ')[1].upper() \
                                                          .replace('(', '').replace(')', '').replace('+', '&'))

lr_pairs['c2c_interaction'] = lr_pairs.apply(lambda row: row['ligand2'] + '^' + row['receptor2'], axis=1)

lr_pairs = c2c.preprocessing.ppi.remove_ppi_bidirectionality(lr_pairs, ('ligand2', 'receptor2'))


ppi_functions = dict()
for idx, row in lr_pairs.iterrows():
    ppi_label = row['ligand2'] + '^' + row['receptor2']
    ppi_functions[ppi_label] = row['annotation']


# **RNA-seq**
rnaseq_matrices = []
# Use sample_meta to preserve the new order of samples
for k, row in tqdm(sample_meta.set_index('sample').iterrows(), total=len(sample_meta)):
    df = data[k]
    rnaseq_matrices.append(c2c.preprocessing.aggregate_single_cells(rnaseq_data=df.fillna(0),
                                                                    metadata=metadata[k],
                                                                    barcode_col='Cell',
                                                                    celltype_col='cell_type',
                                                                    method='nn_cell_fraction',
                                                                    transposed=False))


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
                                      communication_score='expression_mean'
                                     )

# **Metadata for TF-Plot**
meta_tf = c2c.tensor.generate_tensor_metadata(interaction_tensor=tensor,
                                              metadata_dicts=[sample_disease, ppi_functions, None, None],
                                              fill_with_order_elements=True
                                             )


# **Elbow Analysis**
fig, error = tensor.elbow_rank_selection(upper_rank=25,
                                         runs=10,
                                         init='random',
                                         automatic_elbow=False,
                                         filename=output_folder + '/COVID-19-Elbow.png',
                                         random_state=888)

# **Perform tensor factorization**
tensor.compute_tensor_factorization(rank=10,
                                    init='svd',
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
                                             filename=output_folder + '/COVID-19-TensorFactorization.png'
                                            )


# **Export Loadings**
tensor.export_factor_loadings(output_folder + '/COVID-19-Loadings.xlsx')

print("--- %s seconds ---" % (time.time() - start_time))
