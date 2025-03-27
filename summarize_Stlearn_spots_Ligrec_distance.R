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

date = gsub("2025-","25",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

source('important_functions.R')

LR_sender_receiver_spot_df_bysample_V2_files = Sys.glob("LR_sender_receiver_spot_df_bysample_V2.3.dist400.*.rds")

names(LR_sender_receiver_spot_df_bysample_V2_files)= gsub("LR_sender_receiver_spot_df_bysample_V2.3.dist400.|.rds","",basename(LR_sender_receiver_spot_df_bysample_V2_files),perl=TRUE);


#read all stlearn Ligand receptor results at distance of 400px
stlearn_withdistance_sigspots = Sys.glob("/n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/stlearn_cci/stlearn_withdistance/*_significant_barcodes_byLR_dist400.csv")
names(stlearn_withdistance_sigspots) = gsub("TWKI-|TWBK-|_significant_barcodes_byLR_dist400.csv","",basename(stlearn_withdistance_sigspots),perl=TRUE);

neighbourhood_spots = Sys.glob("/n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/stlearn_cci/stlearn_withdistance/*_neighbourhood_spots_byLR_dist400.csv")
names(neighbourhood_spots) = gsub("TWKI-|TWBK-|_neighbourhood_spots_byLR_dist400.csv","",basename(neighbourhood_spots),perl=TRUE);

index <- match(names(neighbourhood_spots),mysample_map$sample_name)
#rename spot names
mynew_names = mysample_map$image_name[index]
names(neighbourhood_spots) = mynew_names
