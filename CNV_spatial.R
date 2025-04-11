library(SingleCellExperiment)
library(scuttle)
library(scran)
library(scater)
library(scRNAseq)
library(Seurat)
# library(BiocParallel)
library(yaml)
library(Ckmeans.1d.dp)
library(stringr)
library(dplyr)
library(sctransform)
library(ggplot2)
library(sctransform)
library(future)
library(RColorBrewer)
# library(circlize)
# library(scCustomize)

# library(circlize)
# library(ComplexHeatmap)
# library(stringr)
library(reshape2)

# library(clusterProfiler)
# library(org.Hs.eg.db)
library("CellTrek")

library(tidyr)
library(dplyr)

set.seed(12345L) 


date = gsub("2025-","25",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

seurat_obj_list2 = readRDS('/n/scratch/users/s/sak4832/March12th_2025/Feb2_2025/Dec27_2024/cluster_distances/seurat_obj_list2.rds')