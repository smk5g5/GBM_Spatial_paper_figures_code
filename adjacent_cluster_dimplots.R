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

load('/n/scratch/users/s/sak4832/March12th_2025/Feb2_2025/Dec27_2024/cluster_distances/jcps_calc_wip2.rda')

selected_cluster_pair = '3_6'


for(i in names(adjacent_cluster_barcodes)){
	Idents(seurat_obj_list2[[i]]) = 'integrated_snn_res.0.9'
spatial_dm_list = list()
if(selected_cluster_pair %in% names(adjacent_cluster_barcodes[[i]])){
select_bars = union(adjacent_cluster_barcodes[[i]][[selected_cluster_pair]]$barcode1,adjacent_cluster_barcodes[[i]][[selected_cluster_pair]]$barcode2)
	mytst = subset(seurat_obj_list2[[i]],cells=select_bars)
  mydmpl = SpatialDimPlot(mytst, image.alpha = 0.8, group.by= 'integrated_snn_res.0.9', pt.size.factor=1.2, stroke=0, crop=FALSE,cols =final_cluster_colors) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+  guides(fill = guide_legend(override.aes = list(size = 10)))
  spatial_dm_list[[i]] = mydmpl
}else{
	next
}

}

lapply(names(spatial_dm_list),function(x) {
  jpeg(paste0('Spatial.integrated_snn_res.0.9_',x,'_nolegend.jpg'), width = 10, height = 7, units="in", res=600);
  print(spatial_dm_list[[x]]);
  dev.off();
  })