# ============================================================================
# Figure 4A code for cluster adjacency analysis
# Author: [Saad Khan]
# Description: Enrichment dotplot of significant contact-dependent (within-spot) ligand receptor interactions
 # in visium data. Size of dot indicates enrichment score, bold circle around spot indicates significance based
 #  on adjusted p.values. Row Barplots show celltype diversity for the said LR whereas column barplots indicate 
 #  celltype diversity in the cluster. The clusters are ordered from edge rich clusters to core-rich clusters
 #  # ============================================================================
##############################
# Load Required Libraries
##############################
# Core single-cell analysis packages
library(SingleCellExperiment)
library(scuttle)
library(scran)
library(scater)
library(Seurat)
library(sctransform)
library(scCustomize)

# Parallelization
library(BiocParallel)
library(future)

# Visualization
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(circlize)
library(ComplexHeatmap)

# Data manipulation
library(dplyr)
library(stringr)
library(reshape2)

# Clustering and analysis
library(Ckmeans.1d.dp)
library(clusterProfiler)
library(org.Hs.eg.db)

# Configuration and YAML support
library(yaml)
# Set seed for reproducibility
set.seed(12345L)

mysample_map = data.frame(sample_name = c("B172-01","B172-02","B176-01-02","B176-02-01","B176-02-04","B178-02-01","B185-01","B185-02","B185-03","B186-01","B186-02","B186-03","B186-04","ST073021-2","ST073021-4","WU1220-1","WU1220-2","WU1220-3","WU1221-Core2","WU1221-Edge1","WU1227-Core5","WU1227-Edge2"))
mysample_map$image_name = make.names(mysample_map$sample_name)

mysample_map$image_name[mysample_map$sample_name=='ST073021-2'] = 'ST073021.02'
mysample_map$image_name[mysample_map$sample_name=='ST073021-4'] = 'ST073021.04'
mysample_map$image_name[mysample_map$sample_name=='WU1220-1'] = 'WU1220.1'
mysample_map$image_name[mysample_map$sample_name=='WU1220-2'] = 'WU1220.02'
mysample_map$image_name[mysample_map$sample_name=='WU1220-1'] = 'WU1220.01'
mysample_map$image_name[mysample_map$sample_name=='WU1220-3'] = 'WU1220.03'
mysample_map$image_name[mysample_map$sample_name=='WU1227-Core5'] = 'WU1227.Core5'
mysample_map$image_name[mysample_map$sample_name=='WU1227-Edge2'] = 'WU1227.Edge2'
mysample_map$image_name[mysample_map$sample_name=='WU1221-Core2'] = 'WU1221.Core'
mysample_map$image_name[mysample_map$sample_name=='WU1221-Edge1'] = 'WU1221.Edge'

mysample_map$image_name[mysample_map$image_name %in% 'WU1221.Core'] = 'WU1221.Core2'
mysample_map$image_name[mysample_map$image_name %in% 'WU1221.Edge'] = 'WU1221.Edge1'


