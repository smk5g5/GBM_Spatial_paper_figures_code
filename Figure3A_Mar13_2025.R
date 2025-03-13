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
library(stringr)library("CellTagR")

library(reshape2)

library(clusterProfiler)
library(org.Hs.eg.db)
library("CellTrek")

set.seed(12345L) 

date = gsub("2025-","25",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

merged.patient.combined.integ_sub = readRDS('/n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/cluster_distances/merged.patient.combined.integ_sub.rds')

spots_to_cells_celltrek_annotations = readRDS('/n/scratch/users/s/sak4832/Dec27_2024/cluster_distances/spots_to_cells_celltrek_annotations.rds')
