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


library(ggdendro)
library(cowplot)
library(tidyverse)
library(ggtree) # install with `devtools::install_github("YuLab-SMU/ggtree")` as you need a version newer than what bioconductor serves
library(patchwork) 
library(gridExtra)
library(gtools)

set.seed(12345L) 

date = gsub("2025-","25",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

#intraniche_interactions
#interniche_interactions
source('important_functions.R')
# cluster_combs_all = unique(unlist(lapply(names(cluster_combs_list_10max), function(x) {
#   myvec = c(cluster_combs_list_10max[[x]]$bidirectional_LR_no2fc,cluster_combs_list_10max[[x]]$bidirectional_LR_grt2fold,cluster_combs_list_10max[[x]]$unique_LR)
#   return(myvec)
# })))

merged.patient.combined.integ_sub = readRDS('/n/scratch/users/s/sak4832/March12th_2025/Feb2_2025/Dec27_2024/cluster_distances/merged.patient.combined.integ_sub.rds')

single_cell_ref = readRDS('/n/data1/mgh/neuro/petti/lab/Users/khan.saad/RPCA_clustering_and_spots_Ligrec_summary/scrna_seq_ref_celltrek_orig_latestcellstateannotations.rds')

#From cellchat database selct only diffusion based LRs (secreted-signalling & ECM related)
CellChatDB_df_noncomp2_withannot = read.csv('/n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/Novemer22nd_2024/stlearn_cci/CellChatDB_df_noncomp2_withannot.csv')

diffusion_based = subset(CellChatDB_df_noncomp2_withannot,annotation %in% c('Secreted Signaling','ECM-Receptor'))


for(i in 1:nrow(intraniche_interactions)){
	myLR = intraniche_interactions$LR[i]
	sg = intraniche_interactions$sender_group[i]
	rg = intraniche_interactions$receiver_group[i]
	mysub = subset(cluster_distances_with_LR_comb_grouped,LR==myLR&sender_group==sg&receiver_group==rg)
	mydf1 = data.frame(Sample=str_split_i(mysub$sender_spot2,'_',1),spot_id=str_split_i(mysub$sender_spot2,'_',2))
	mydf2 = data.frame(Sample=str_split_i(mysub$receiver_spot2,'_',1),spot_id=str_split_i(mysub$receiver_spot2,'_',2))
	spot_data = unique(rbind(mydf1,mydf2))
	spot_data_list = split(spot_data,spot_data$Sample)
	celltrek_spots_list = lapply(names(spot_data_list),function(x) {
	celltrek_spots = intersect(spot_data_list[[x]]$spot_id,spots_to_cells_celltrek_annotations[[x]]$broad_celltypes_per_spot$spot_id)
	return(celltrek_spots)
	})
	celltrek_spots_len = sum(unlist(lapply(celltrek_spots_list,length)))
	intraniche_interactions$unique_spot_count[i] = length(unique(union(mysub$sender_spot2,mysub$receiver_spot2)))
	intraniche_interactions$celltrek_spot_count[i] = celltrek_spots_len
}

for(i in 1:nrow(interniche_interactions)){
	myLR = interniche_interactions$LR[i]
	sg = interniche_interactions$sender_group[i]
	rg = interniche_interactions$receiver_group[i]
	mysub = subset(cluster_distances_with_LR_comb_grouped,LR==myLR&sender_group==sg&receiver_group==rg)
	mydf1 = data.frame(Sample=str_split_i(mysub$sender_spot2,'_',1),spot_id=str_split_i(mysub$sender_spot2,'_',2))
	mydf2 = data.frame(Sample=str_split_i(mysub$receiver_spot2,'_',1),spot_id=str_split_i(mysub$receiver_spot2,'_',2))
	spot_data = unique(rbind(mydf1,mydf2))
	spot_data_list = split(spot_data,spot_data$Sample)
	celltrek_spots_list = lapply(names(spot_data_list),function(x) {
	celltrek_spots = intersect(spot_data_list[[x]]$spot_id,spots_to_cells_celltrek_annotations[[x]]$broad_celltypes_per_spot$spot_id)
	return(celltrek_spots)
	})
	celltrek_spots_len = sum(unlist(lapply(celltrek_spots_list,length)))
	interniche_interactions$unique_spot_count[i] = length(unique(union(mysub$sender_spot2,mysub$receiver_spot2)))
	interniche_interactions$celltrek_spot_count[i] = celltrek_spots_len
}
