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
library("dplyr")
library("Seurat")
library("viridis")
library("RColorBrewer")
library("ggplot2")

library(ComplexHeatmap)
library(cola)
library(biclust)

set.seed(1234567890L)


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
mysample_map$image_name[mysample_map$sample_name=='WU1221-Core2'] = 'WU1221.Core2'
mysample_map$image_name[mysample_map$sample_name=='WU1221-Edge1'] = 'WU1221.Edge1'



clusterprofiler_input2 = clusterprofiler_input %>% left_join(mysample_map, by = 'sample_name')

############################################################
# Input for Consensus clustering analysis
# Celltrek output of (spatially variable) tumor state clusters 
# We created 5 consensus clusters (or tumor states) common across patients
############################################################

clusterprofiler_input = read.table('/n/data1/mgh/neuro/petti/lab/Users/khan.saad/scratch_backup/GBM_scrna_spatial/GBM_spatial_clusterprofiler.txt',sep = "\t",header = F)
colnames(clusterprofiler_input) = c('rds_file','sample_name','tumor_cluster')


Cluster_profiler_list = list()

for(i in 1:nrow(clusterprofiler_input2)){
  input_file = read.table(clusterprofiler_input2$rds_file[i],header = T,sep=",")
  sample_cluster = clusterprofiler_input2$image_name[i]
  Tumor_cluster = clusterprofiler_input2$tumor_cluster[i]
  Cluster_profiler_list[[paste0(sample_cluster,".",Tumor_cluster)]] =   input_file$x
}

Cluster_profiler_binary_htmap = list_to_matrix(Cluster_profiler_list)
# 
#final reduced gene list for each clustered tumor state
#https://github.com/Petti-Lab/GBM_Spatial_paper_figures_code/blob/main/Figure_2/final_genemodules_reduced_mutualinfo_list.rds
final_genemodules_reduced_mutualinfo_list = readRDS('/n/scratch/users/s/sak4832/final_genemodules_reduced_mutualinfo_list.rds')

module_genes = unique(unlist(final_genemodules_reduced_mutualinfo_list))


Cluster_profiler_binary_htmap_sub = Cluster_profiler_binary_htmap[module_genes,]

# Colors and metadata setup
# https://github.com/Petti-Lab/GBM_Spatial_paper_figures_code/blob/main/Figure_2/color_panel_list.rds
color_panel = readRDS('/n/scratch/users/s/sak4832/color_panel_list.rds')

sample_colors = color_panel$Sample
names(sample_colors) = mysample_map$image_name

sample_ord  = str_split_i(colnames(Cluster_profiler_binary_htmap_sub),'.cc',1)

pat_ord =  str_split_i(colnames(Cluster_profiler_binary_htmap_sub),'\\.',1)

# https://github.com/Petti-Lab/GBM_Spatial_paper_figures_code/blob/main/Figure_2/spatial_data_metadata.txt
spatial_data_metadata  = read.table('/n/scratch/users/s/sak4832/spatial_data_metadata.txt',header=T,sep="\t")

kitmeta = merge(spatial_data_metadata[c('sampe_name','Kit')],mysample_map,by.x ='sampe_name',by.y='sample_name')

mydf = data.frame(sample_name=str_split_i(colnames(Cluster_profiler_binary_htmap_sub),'.cc',1))

kit_ord_df = mydf %>% left_join(kitmeta, by = c('sample_name'='image_name'))

kit_ord = kit_ord_df$Kit

patcols = color_panel$Patient

kit_cols = color_panel$Kit

# ha1 = HeatmapAnnotation(foo1 = runif(10), bar1 = sample(c("f", "m"), 10, replace = TRUE),
#     annotation_legend_param = list(
#         foo1 = list(direction = "horizontal"),
#         bar1 = list(nrow = 1)))


colors = structure(c("white", "red"), names = c('0','1') ) # black, red, green, blue


row_split_list = unlist(lapply(names(final_modules2), function(x) {
	rep(x,length(final_modules2[[x]]))
  }))

names(row_split_list) = unname(unlist(final_modules2))
#



V1_sel = c('DPYSL3','CEP170','KIF1B','NFIX','KMT2A','CREBBP','SMAD4','PBRM1','TRIM28','YEATS2','KMT2C')
V2_sel = c("STMN1","MCM7","MCM3","PCNA","CTCF","SMC4","SMC3","SMC1A","NFIB","CENPV","CENPX","CENPF","CDKN2C","MYH10")
V3_sel = c("VIM","CHI3L1","IGFBP7","SERPINA3","SERPING1","SPP1","APOE","MT2A","S100A16","CD63","CD99","HLA-C","CST3","CTSD","CALR","CTSF","B2M")
V4_sel = c("VEGFA","CD44","SOD2","EGLN3","ERO1A","ADM","PLOD1","PLOD2","AK4","CA9","CAV1","IRS2")
V5_sel = c("SOX4","SOX11","DLL3","OLIG2","PDGFRA","BCAN","VCAN","PTPRZ1","THY1","TTYH1","MARCKSL1")

selected_genes = intersect(rownames(Cluster_profiler_binary_htmap_sub),unique(c(V1_sel,V2_sel,V3_sel,V4_sel,V5_sel)))

gene_indices <- match(selected_genes,rownames(Cluster_profiler_binary_htmap_sub))
# gene_annots = rowAnnotation(selgenes = anno_mark(at = gene_indices, labels = selected_genes,which = "row", side = "right"))
gene_annots = rowAnnotation(
  selgenes = anno_mark(
    at = gene_indices,
    labels = selected_genes,
    which = "row",
    side = "right",
    labels_gp = gpar(fontsize = 25)  # Set font size for the annotation labels
  )
)


# Heatmap(mat, name = "mat", column_km = 3)
# row_names_gp = gpar(fontsize = 10),
# row_names_side = "left", row_names_gp = gpar(fontsize = 4))
#
 ha = HeatmapAnnotation(
Patient = pat_ord,
Sample = sample_ord,
Kit = kit_ord,
col = list(Patient = patcols,
           Sample = sample_colors,
           Kit = kit_cols,simple_anno_size = unit(5, "cm")), height = unit(5, "cm"),annotation_name_side="left",annotation_name_gp= gpar(fontsize = 20),
annotation_legend_param = list(
	Patient = list(direction = "vertical",legend_gp = gpar(fontsize = 8),labels_gp = gpar(fontsize = 25,nrow = 1)),
	Sample = list(direction = "vertical",legend_gp = gpar(fontsize = 8),labels_gp = gpar(fontsize = 25,nrow = 2)),
	Kit = list(direction = "vertical",legend_gp = gpar(fontsize = 8),labels_gp = gpar(fontsize = 25,nrow = 3))),simple_anno_size_adjust = TRUE)


binhtmap = Heatmap(Cluster_profiler_binary_htmap_sub, name = c("binary matrix\nspatially variable\ntumor cluster DEGs"),
                             show_column_dend = FALSE,row_split=row_split_list, column_order = col_ord,cluster_column=F,
                             cluster_column_slices = F,col = colors,rect_gp = gpar(col = "black", lwd = 1),
                             row_names_side = "right",right_annotation = gene_annots,
                             column_names_gp = gpar(fontsize = 12),
                             column_title_gp = gpar(fontsize = 15),
                             column_gap = unit(0.5, "mm"),
                             cluster_rows = T,cluster_row_slices = F,
                             top_annotation = ha,
                             show_column_names = F,
                             # column_split = factor(x = merged.combined_sub@meta.data[[integ]], levels =levels(merged.combined_sub@meta.data[[integ]])),
                             show_row_dend = FALSE,
                             heatmap_legend_param = list(direction = "horizontal",title ="rpca integrated clusters", legend_gp = gpar(fontsize = 8),legend_width = unit(5, "cm"),labels_gp = gpar(fontsize = 8)),
                             row_title_gp = gpar(fontsize = 20),
                             show_row_names = F,
                             column_title_rot = 90,
                             column_names_rot = 90)


ht_opt("heatmap_row_names_gp" = gpar(fontsize = 25))


jpeg('Module_genes_binary_heatmap_spatially_variable_tumor_clusters_nolegend.jpg',width = 70,height = 70,units="cm", res=600)
draw(binhtmap,show_annotation_legend = FALSE, show_heatmap_legend = FALSE)
dev.off()


jpeg('Module_genes_binary_heatmap_spatially_variable_tumor_clusters_withlegend.jpg',width = 70,height = 90,units="cm", res=600)
draw(binhtmap,heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()