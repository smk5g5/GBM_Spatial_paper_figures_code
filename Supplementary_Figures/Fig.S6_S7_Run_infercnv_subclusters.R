library(infercnv)
library(ggplot2)
library(patchwork)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(SingleCellExperiment)
library(scuttle)
library(scran)
library(scater)
library(scRNAseq)
library(Seurat)
library(yaml)
library(Ckmeans.1d.dp)
library(stringr)
library(dplyr)
library(sctransform)
library(ggplot2)
library(sctransform)
library(future)
library(RColorBrewer)
library(reshape2)
library(tidyr)
library(dplyr)

set.seed(12345L)

date = gsub("2025-","25",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

gene_order_file = readRDS('gene_order_file.rds')

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  args <- c("--help")
}

create_dir <- function(inp_dir) {

if (!dir.exists(inp_dir)){
  dir.create(inp_dir)
} else {
  print("Dir already exists!")
}
}

prefix <- as.character(args[1])

outdir = paste0(getwd(),'/',prefix,'_infercnv2_clonepubversion_0.01_leiden_patient_subclusters/')

seurat_obj_list2 = readRDS('/n/scratch/users/s/sak4832/CNV_spatial/CNV_score_seurat_obj_list2.rds')

select_samp = grep(prefix,names(seurat_obj_list2),value=T)

if(length(select_samp)>1){
selected_seurat_obj_list = seurat_obj_list2[select_samp]

merged <- merge(x = selected_seurat_obj_list[[1]], y = selected_seurat_obj_list[2:length(selected_seurat_obj_list)], add.cell.ids = select_samp)

mydf = merged@meta.data[c('perspot_CNV_classification')]
mydf_sub = mydf[c('perspot_CNV_classification')]

infercnv_data_inp = as.matrix(merged@assays$Spatial@counts)
}else{
mydf = seurat_obj_list2[select_samp]@meta.data[c('perspot_CNV_classification')]
mydf_sub = mydf[c('perspot_CNV_classification')]
infercnv_data_inp = as.matrix(seurat_obj_list2[select_samp]@assays$Spatial@counts)
}


infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=infercnv_data_inp,
                                            gene_order_file=gene_order_file,
                                            annotations_file=mydf_sub,
                                            ref_group_names=c("Low_CNV"))


window_length = 151
sim_method = "meanvar"
pnorm_bayes = 0.3

options(scipen = 100)
infercnv_obj <- infercnv::run(infercnv_obj,window_length = 151,
    cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
    out_dir = outdir, # output files
    cluster_by_groups = T, # cluster or not
    analysis_mode = "subclusters",sim_method=sim_method,
    denoise = T,leiden_resolution=0.05,
    HMM = T,plot_steps=T,BayesMaxPNormal = pnorm_bayes,
    min_cells_per_gene = 3,
    tumor_subcluster_partition_method = "leiden",
    HMM_type = "i6",
    num_threads = 10,
    plot_probabilities = TRUE,
    #plot_probabilities = F,
    save_rds = T,
    save_final_rds = T,
    #no_plot = T,
    #output_format = NA,
    useRaster = T,
    #up_to_step = 17
  )




