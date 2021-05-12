
suppressPackageStartupMessages({
    library(CellChat)
    library(patchwork)
    library(RhpcBLASctl)
    library(Matrix)
})
options(stringsAsFactors = FALSE)

RhpcBLASctl::blas_set_num_threads(1) # no multithreading

# ARGS
args <- commandArgs(trailingOnly=TRUE)
data_path = paste0(args[1], "/data/")
output_path = args[2]

# load inputs
data_input_map<-readRDS(paste0(data_path, 'data_input_map_by_sample.rds'))
meta_map<-readRDS(paste0(data_path, 'meta_map_by_sample.rds'))
humandb<-readRDS(paste0(data_path, 'humandb.rds'))

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
for (sample.name in names(meta_map)){
    # loop through each sample.name and create a cell type future
    cellchat <- createCellChat(object = data_input_map[[sample.name]], meta = meta_map[[sample.name]], 
                               group.by = 'cell_type')
    cellchat@DB <- humandb # human organism

    cellchat <- subsetData(cellchat) # subset the expression data of signaling genes, assign to @data.signalling 
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat) # generate @ LR slot used by computeCommunProb
    cellchat <- projectData(cellchat, PPI.human) # shallow sequencing depth
    
    cellchat <- computeCommunProb(cellchat, raw.use = F, type = 'triMean', trim = NULL, seed.use = NULL, 
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

path.dist <- list()
for (pc in names(pairwise_comparisons)){
    print(pc)
    path.dist[[pc]] <- rankSimilarity_(cellchat, type = 'structural', comparison1 = 1:length(covid.list),
                                       comparison2 = pairwise_comparisons[[pc]])
}                      

# Save file
saveRDS(cellchat, file = paste0(output_path, "cellchat_bysample.rds"))
saveRDS(path.dist, file = paste0(output_path, "pwc_bysample.rds"))