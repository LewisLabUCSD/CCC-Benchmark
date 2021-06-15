#! /usr/bin/env Rscript
'Time how long it takes to run cellchat given some samples
Usage:
    04_cellchat_timing.r --number=<int> --group=<bool> --root=<str> --outputs=<str> [--seed=<int>]
    
Options:
    -h --help  Show this screen.
    --number=<int>  number of samples to include
    --group=<bool>  whether to group by context (1) or keep separated by samples (0)
    --root=<str>  folder path to root folder of this timing analysis
    --outputs=<str> output folder path
    --seed=<int>  set a seed for random graph generation
' -> doc

suppressPackageStartupMessages({
    library(CellChat, quietly = T)
    library(patchwork, quietly = T)
    library(RhpcBLASctl, quietly = T)
    library(Matrix, quietly = T)
    library(docopt, quietly = T)
    library(rhdf5, quietly = T)
    library(stringr, quietly = T)
})
options(stringsAsFactors = FALSE)
# expression_data_path = '/data2/hratch/immune_CCI/covid/covid_atlas/interim/umi_for_timing/'
#input_path = '/data2/hratch/immune_CCI/covid/covid_atlas/interim/timing_inputs/'

RhpcBLASctl::blas_set_num_threads(1) # no multithreading

cell_grouper<-'majorType'

opts = docopt(doc)
for (opt in c('number', 'group', 'seed')){
    val = opts[[opt]]
    if (!is.null(val)){
        opts[[opt]] = as.numeric(val)
        }}

number<-opts[['number']]
group<-opts[['group']]
seed<-opts[['seed']]
input_path<-paste0(opts[['root']], "/data/PBMCs/")
output_path<-opts[['outputs']]
set.seed(seed)

# Print info
if (group){
    agg <- "_samples_agg"
}else{
    agg <- "_samples_ind"
}
print(paste0("Results will be saved to: ", paste0(output_path, "cellchat_", number, agg, ".rds")))


samples<-read.csv(paste0(input_path, 'samples_for_timing.csv'))
sample.names<-unlist(strsplit(samples[samples$n_samples == number, 'sample_names'], '; '))

# load LR pairs
# filter for the LR pairs used by Tensor cell2cell
# lr_pairs<-read.csv(paste0(input_data_path,'Tensor-cell2cell-LRpairs.csv'))
# lr_pairs<-lr_pairs$interaction_name
humandb<-CellChatDB.human
# humandb$interaction<-CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in% lr_pairs, ] 
# saveRDS(humandb, paste0(output_results_path, 'humandb.rds'))

# load metadata
md.cell <- read.csv(paste0(input_path,'metadata_for_timing.csv'), row.names = 1)
sample_context_map<-readRDS(paste0(input_path, 'sample_context_map.rds'))
contexts<-unique(sample_context_map)

# # load the UMI counts
# if (!group){
#     counts<-lapply(setNames(sample.names, sample.names), function(sn) t(read.csv(paste0(expression_data_path, sn, '.csv'), row.names = 1)))
# }else{ # group by context
#     by.context<-lapply(setNames(contexts, contexts), function(context) names(sample_context_map[sample.names][sample_context_map[sample.names] == context]))
    
#     group_by_context<-function(context){
#         sns<-by.context[[context]]       
#         counts<-lapply(sns, function(sn) t(read.csv(paste0(expression_data_path, sn, '.csv'), row.names = 1)))    
#         counts<-do.call(cbind, counts)
#         return (counts)
#     }
#     counts<-lapply(setNames(contexts, contexts), function(context) group_by_context(context))                   
                    
# }

# load the UMI counts
read_sample_h5<-function(sn){
    counts<-h5read(paste0(input_path, 'umi_per_sample.h5'), sn)
    count<-counts[[4]]
    colnames(count)<-counts[[2]]
    rownames(count)<-sapply(counts[[1]], function(x) str_replace_all(x, '-', '.')) 
    return(count)
}


if (!group){                             
    counts<-lapply(setNames(sample.names, sample.names), function(sn) read_sample_h5(sn))
                                  
}else{ # group by context
    by.context<-lapply(setNames(contexts, contexts), function(context) names(sample_context_map[sample.names][sample_context_map[sample.names] == context]))
    
    group_by_context<-function(context){
        sns<-by.context[[context]]       
        counts<-lapply(sns, function(sn) read_sample_h5(sn))    
        counts<-do.call(cbind, counts)
        return (counts)
    }
    counts<-lapply(setNames(contexts, contexts), function(context) group_by_context(context))                   
                    
}
                   
# filter for intersection of cell types across samples/contexts - to make comparable with Tensor-cell2cell
# cell.types<-Reduce(intersect, lapply(sample.names, function(sn) unique(md.cell[md.cell$sampleID == sn, cell_grouper])))
cell.types<-Reduce(intersect, lapply(counts, function(df) unique(md.cell[colnames(df), cell_grouper])))
cell.ids<-rownames(md.cell[md.cell[[cell_grouper]] %in% cell.types, ])   
counts<-lapply(setNames(counts, counts), function(df) df[colnames(df) %in% cell.ids,])                      

#' Rank the similarity of the shared signaling pathways based on their joint manifold learning
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param comparison1 a numerical vector giving the datasets for comparison. This should be the same as `comparison` in `computeNetSimilarityPairwise`
#' @param comparison2 a numerical vector with two elements giving the datasets for comparison.
#'
#' If there are more than 2 datasets defined in `comparison1`, `comparison2` can be defined to indicate which two datasets used for computing the distance.
#' e.g., comparison2 = c(1,3) indicates the first and third datasets defined in `comparison1` will be used for comparison.
#' @import ggplot2
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
rankSimilarity_ <- function(object, slot.name = "netP", type = c("functional","structural"), comparison1 = NULL,  
                           comparison2 = c(1,2)) {
  type <- match.arg(type)

  if (is.null(comparison1)) {
    comparison1 <- 1:length(unique(object@meta$datasets))
  }
  comparison.name <- paste(comparison1, collapse = "-")
  cat("Compute the distance of signaling networks between datasets", as.character(comparison1[comparison2]), '\n')
  comparison2.name <- names(methods::slot(object, slot.name))[comparison1[comparison2]]

  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  group <- sub(".*--", "", rownames(Y))
  data1 <- Y[group %in% comparison2.name[1], ]
  data2 <- Y[group %in% comparison2.name[2], ]
  rownames(data1) <- sub("--.*", "", rownames(data1))
  rownames(data2) <- sub("--.*", "", rownames(data2))

  pathway.show = as.character(intersect(rownames(data1), rownames(data2)))
  data1 <- data1[pathway.show, ]
  data2 <- data2[pathway.show, ]
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  dist <- NULL
  for(i in 1:nrow(data1)) dist[i] <- euc.dist(data1[i,],data2[i,])
  df <- data.frame(name = pathway.show, dist = dist, row.names = pathway.show)
  df <- df[order(df$dist), , drop = F]
  df$name <- factor(df$name, levels = as.character(df$name))

  return(df)
}

# create cellchat object for each sample or sample.name
covid.list<-list()
for (sample.name in names(counts)){
    # loop through each sample.name and create a cell type future
    expr<-CellChat::normalizeData(counts[[sample.name]])
    cellchat<-createCellChat(object = as(expr, "dgCMatrix"), meta = md.cell[colnames(expr),], 
                                   group.by = cell_grouper)
    cellchat@DB <- humandb # human organism

    cellchat <- subsetData(cellchat) # subset the expression data of signaling genes, assign to @data.signalling 
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat) # generate @ LR slot used by computeCommunProb
    cellchat <- projectData(cellchat, PPI.human) # shallow sequencing depth
    
    cellchat <- computeCommunProb(cellchat, raw.use = F, type = 'triMean', trim = NULL, seed.use = seed, 
                                 population.size = F) 
    
    # The functional similarity analysis requires the same cell population composition between two datasets.
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    covid.list[[sample.name]]<-cellchat
}

# merge and analyze
cellchat <- mergeCellChat(covid.list, add.names = names(covid.list))
cellchat <- computeNetSimilarityPairwise(cellchat, type = 'structural')
cellchat <- netEmbedding(cellchat, type = 'structural')
cellchat <- netClustering(cellchat, type = 'structural',  do.parallel = F, do.plot = F)
# Manifold learning of the signaling networks for datasets 
pairwise_comparisons<-sapply(as.data.frame(combn(seq(1:length(covid.list)),2)), 
                         function(x) as.numeric(x), simplify = F) # pairwise combination of elements

for (pc in names(pairwise_comparisons)){
    path.dist<-rankSimilarity_(cellchat, type = 'structural', comparison1 = 1:length(covid.list),
                               comparison2 = pairwise_comparisons[[pc]]) 
}

# Save file
saveRDS(cellchat, file = paste0(output_path, "cellchat_", number, agg, ".rds"))
saveRDS(path.dist, file = paste0(output_path, "pwc_", number, agg, ".rds"))
