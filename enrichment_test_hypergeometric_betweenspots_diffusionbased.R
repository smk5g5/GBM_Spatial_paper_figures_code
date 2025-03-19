library(SingleCellExperiment)
library(scuttle)
library(scran)
library(scater)
library(scRNAseq)
library(Seurat)
library(BiocParallel)
library(yaml)
library(Ckmeans.1d.dp)
library(stringr)
library(dplyr)
library(sctransform)
library(ggplot2)
library(sctransform)
library(future)
library(RColorBrewer)
library(circlize)
library(scCustomize)

library(circlize)
library(ComplexHeatmap)
library(stringr)

library(reshape2)

library(clusterProfiler)
library(org.Hs.eg.db)
library("CellTrek")

set.seed(12345L) 

### Load precomputed datasets###
# cluster_meta  contains cluster metadata
cluster_meta = readRDS('/n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/cluster_distances/cluster_meta.rds')

# cluster_distances is a melted matrix of precomputed euclidean distances between spots for each sample

cluster_distances = readRDS('/n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/cluster_distances/cluster_distances_wip_list.rds')

groupdf = readRDS('groupdf')