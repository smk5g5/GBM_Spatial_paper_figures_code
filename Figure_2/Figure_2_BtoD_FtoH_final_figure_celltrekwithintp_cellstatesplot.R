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
set.seed(12345L) 

############################################################
# Input file list for CellTrek objects
############################################################

celltrek_df = read.table('celltrek_files.txt',header=F)

date = gsub("2025-","25",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

# Map sample names to sanitized image names used by Seurat/CellTrek

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

# Join file list with sample mapping; load CellTrek objects

celltrek_df_merged = merge(celltrek_df,mysample_map,by.x='V1',by.y='sample_name')

celltrek_obj_list = list()

for(i in 1:nrow(celltrek_df_merged)){
	myj = celltrek_df_merged$image_name[i]
	celltrek_obj_list[[myj]] = readRDS(celltrek_df_merged$V2[i])
}

############################################################
# make_broad_celltypes2: collapse finegrained cell types to broad categories
# - Uses `single_cell_ref$Celltypes_Aug8_2024` present in object
# - Adds `broad_celltypes` to the Seurat object metadata
############################################################

make_broad_celltypes2 <- function(single_cell_ref){
  mydf = data.frame(names(table(single_cell_ref$Celltypes_Aug8_2024)))
  names(mydf) = c('all_celltypes')
  mydf$broad_celltypes = NA
  mydf$broad_celltypes[grep('CD4_',mydf$all_celltypes)] = 'CD4_Tcells'
  mydf$broad_celltypes[grep('CD8_',mydf$all_celltypes)] = 'CD8_Tcells'
  mydf$broad_celltypes[grep('MAIT-like|Proliferating|T_EMRA|T_IFNG_stimulated',mydf$all_celltypes)] = 'Other_Tcells'
  mydf$broad_celltypes[grep('TAM',mydf$all_celltypes)] = 'TAMs'
  mydf$broad_celltypes[grep('Vascular|Pericyte',mydf$all_celltypes)] = 'Vascular'
  mydf$broad_celltypes[grep('Fibroblast',mydf$all_celltypes)] = 'Fibroblast'
  mydf$broad_celltypes[grep('Microglia',mydf$all_celltypes)] = 'Microglia'
  mydf$broad_celltypes[grep('Monocyte',mydf$all_celltypes)] = 'Monocyte'
  mydf$broad_celltypes[grep('Choroid plexus|Ependymal|Ex. Neuron|Inhib. Neuron|Neuron',mydf$all_celltypes)] = 'Neuron'
  mydf$broad_celltypes[grep('DC',mydf$all_celltypes)] = 'DC'
  mydf$broad_celltypes[mydf$all_celltypes=='Oligodendrocyte'] = 'Oligodendrocyte'	
  mydf$broad_celltypes[mydf$all_celltypes=='Oligodendrocyte precursor'] = 'Oligodendrocyte precursor'	
  mydf$broad_celltypes[mydf$all_celltypes=='Committed oligodendrocyte precursor'] = 'Committed oligodendrocyte precursor'	
  mydf$broad_celltypes[mydf$all_celltypes=='Astrocyte'] = 'Astrocyte'	
  mydf$broad_celltypes[mydf$all_celltypes=='B cells'] = 'Bcells'
  mydf$broad_celltypes[mydf$all_celltypes=='NK'] = 'NK/NK-like'
  mydf$broad_celltypes[mydf$all_celltypes=='NK-like'] = 'NK/NK-like'
  mydf$broad_celltypes[mydf$all_celltypes=='M1'] = 'M1'
  mydf$broad_celltypes[mydf$all_celltypes=='M2'] = 'M2'
  mydf$broad_celltypes[mydf$all_celltypes=='M3'] = 'M3'
  mydf$broad_celltypes[mydf$all_celltypes=='M4'] = 'M4'
  mydf$broad_celltypes[mydf$all_celltypes=='M5'] = 'M5'
  mydf2 = data.frame(cellnames=names(single_cell_ref$Celltypes_Aug8_2024),Celltypes_Aug8_2024=unname(single_cell_ref$Celltypes_Aug8_2024))
  
  mydf12 = mydf2 %>% left_join(mydf, by = c("Celltypes_Aug8_2024"="all_celltypes"))

  index = match(Cells(single_cell_ref),mydf12$cellnames)
  single_cell_ref$broad_celltypes = mydf12$broad_celltypes[index]
  return(single_cell_ref)
}

############################################################
# Add broad cell type labels to each object
############################################################

celltrek_obj_list2 = lapply(celltrek_obj_list,function(x) {
return(make_broad_celltypes2(single_cell_ref=x))
})


############################################################
# Make tumor-only subsets (M1-M5) for spatial plots
############################################################

celltrek_obj_list_sub = lapply(celltrek_obj_list2,function(x) {
Idents(x) = 'broad_celltypes'
x_sub = subset(x,idents=paste0('M',1:5))
return(x_sub)
})


# Define colors for broad cell types and tweak a few

cellcolor_names = c("M1","M2","M3","M4","M5","TAMs",
	"Microglia","Monocyte","CD4_Tcells","CD8_Tcells","Other_Tcells","Vascular","Fibroblast","Neuron","Astrocyte","Committed oligodendrocyte precursor","Oligodendrocyte","Oligodendrocyte precursor","DC","Bcells","NK/NK-like")
latest_cell_colors2 = c('#E58606','#6BAED6','#E7BA52','#5254A3','#D66B6F','#70CCB0','#91C1B2','#91C292',
	'#3969AC','#4C8CE6','#82C2E8','#FF4219','#D781B5','#E099CA','gray90','gray50','grey30','gray10',"#ffe119","khaki","khaki2")


mycolor_df = data.frame(celltype=cellcolor_names,cellcolors=latest_cell_colors2)

mycolor_df$cellcolors[mycolor_df$celltype=='Fibroblast'] = 'salmon'
mycolor_df$cellcolors[mycolor_df$celltype=='Neuron'] = 'red4'

latest_cell_colors_final = mycolor_df$cellcolors
names(latest_cell_colors_final) = mycolor_df$celltype

############################################################
# Plot tumor states only (M1-M5) per sample with chosen palette (S1 to S5 here shown as M1 to M5)
############################################################

for(j in c("broad_celltypes")){
spatial_dm_list = list()
for(i in names(celltrek_obj_list_sub)){
  Idents(celltrek_obj_list_sub[[i]]) = j
  # scale_fill_manual(values = cluster_colors_list[[j]], breaks = names(cluster_colors_list[[j]]),labels = names(cluster_colors_list[[j]]) ) 
  # spatial_dm_list[[i]] = mydmpl
  jpeg(sprintf("%s/Spatial.Broad.Crop.cellstates.%s.jpg", outdir,i), width = 11, height = 7, units="in", res=300);
  print(SpatialDimPlot(celltrek_obj_list_sub[[i]], image.alpha = 0.7, group.by= j, pt.size.factor=1, stroke=0.09, repel=TRUE,crop=TRUE,cols =latest_cell_colors_final) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+  guides(fill = guide_legend(override.aes = list(size = 10))))
  dev.off()

}
}

outdir='/n/scratch/users/s/sak4832/final_figures_celltrek/allcells'

############################################################
# Plot all cells with full palette
############################################################

for(j in c("broad_celltypes")){
spatial_dm_list = list()
for(i in names(celltrek_obj_list2)){
  Idents(celltrek_obj_list2[[i]]) = j
  # scale_fill_manual(values = cluster_colors_list[[j]], breaks = names(cluster_colors_list[[j]]),labels = names(cluster_colors_list[[j]]) ) 
  mydmpl = SpatialDimPlot(celltrek_obj_list2[[i]], image.alpha = 0.7, group.by= j, pt.size.factor=0.8, stroke=0.09, repel=TRUE,crop=TRUE,cols =latest_cell_colors_final) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+  guides(fill = guide_legend(override.aes = list(size = 10))) 
  # spatial_dm_list[[i]] = mydmpl
	jpeg(sprintf("%s/Spatial.Broad.Crop.allcells.%s.jpg", outdir,i), width = 13, height = 7, units="in", res=300);
	print(mydmpl)
	dev.off()

}
}

############################################################
# Variant palette: gray glia, white tumor (for contrast maps)
############################################################

mycolor_df2 = mycolor_df
mycolor_df2$cellcolors[mycolor_df2$celltype %in% c("Astrocyte","Committed oligodendrocyte precursor","Oligodendrocyte","Oligodendrocyte precursor")] = 'gray50'
mycolor_df2$cellcolors[mycolor_df2$celltype %in% c(paste0('M',1:5))] = 'white'

latest_cell_colors_final = mycolor_df2$cellcolors
names(latest_cell_colors_final) = mycolor_df2$celltype


outdir='/n/scratch/users/s/sak4832/final_figures_celltrek/allcells_tumorwhite'

for(j in c("broad_celltypes")){
spatial_dm_list = list()
for(i in names(celltrek_obj_list2)){
  Idents(celltrek_obj_list2[[i]]) = j
  # scale_fill_manual(values = cluster_colors_list[[j]], breaks = names(cluster_colors_list[[j]]),labels = names(cluster_colors_list[[j]]) ) 
  mydmpl = SpatialDimPlot(celltrek_obj_list2[[i]], image.alpha = 0.7, group.by= j, pt.size.factor=0.8, stroke=0.09, repel=TRUE,crop=TRUE,cols =latest_cell_colors_final) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+  guides(fill = guide_legend(override.aes = list(size = 10))) 
  # spatial_dm_list[[i]] = mydmpl
	jpeg(sprintf("%s/Spatial.Broad.Crop.allcells.tumorwhite.%s.jpg", outdir,i), width = 13, height = 7, units="in", res=300);
	print(mydmpl)
	dev.off()

}
}

# mylegend <- get_legend(
#   # create some space to the left of the legend
#   spatial_dm_list[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12))
# )
library(cowplot)
spatial_dm_list2 = lapply(spatial_dm_list,function(x) {return(x+ theme(legend.position="none") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))})

# pgrid <-do.call(plot_grid, c(spatial_dm_list2, nrow = 5, ncol = 5,greedy=TRUE))
# p <- plot_grid(pgrid,ncol = 2, rel_widths = c(1, .2))

pgrid <-do.call(plot_grid, c(spatial_dm_list2, nrow = 5, ncol = 5))

ggsave(paste0('Spatial_Dimplots_clusters_celltrekwithinterpolation_cellstatesonly_pt.size.factor.0.8_nostroke_cropfalse_',j,'.',date,'.pdf'),pgrid , width = 45, height = 50, units="in",limitsize = FALSE)
# }

############################################################
# Save individual plots with legends (one per page)
############################################################

outdir='/n/scratch/users/s/sak4832/final_figures_celltrek'


for(i in names(spatial_dm_list)){
print(spatial_dm_list[[i]]);
dev.off();
}


library(gridExtra)

ml <- marrangeGrob(spatial_dm_list, nrow=1, ncol=1)
ggsave(paste0('Spatial_Dimplots_clusters_celltrekwithinterpolation_cellstatesonly_pt.size.factor.0.8_nostroke_withlegend',j,'.',date,'.pdf'), ml, width = 11, height = 8, units="in",limitsize = FALSE)




for(j in c("broad_celltypes")){
spatial_dm_list = list()
for(i in names(celltrek_obj_list2)){
  Idents(celltrek_obj_list2[[i]]) = j
  # scale_fill_manual(values = cluster_colors_list[[j]], breaks = names(cluster_colors_list[[j]]),labels = names(cluster_colors_list[[j]]) ) 
  mydmpl = SpatialDimPlot(celltrek_obj_list2[[i]], image.alpha = 0.8, group.by= j, pt.size.factor=0.8, stroke=0.0, repel=TRUE,crop=TRUE,cols =latest_cell_colors_final) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+  guides(fill = guide_legend(override.aes = list(size = 10))) 
  spatial_dm_list[[i]] = mydmpl
}
}
# mylegend <- get_legend(
#   # create some space to the left of the legend
#   spatial_dm_list[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12))
# )
library(cowplot)
spatial_dm_list2 = lapply(spatial_dm_list,function(x) {return(x+ theme(legend.position="none") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))})

pgrid <-do.call(plot_grid, c(spatial_dm_list2, nrow = 5, ncol = 5,greedy=TRUE))
# p <- plot_grid(pgrid, mylegend, ncol = 2, rel_widths = c(1, .2))

ggsave(paste0('Spatial_Dimplots_clusters_celltrekwithinterpolation_allcelltypes_pt.size.factor.0.8_stroke.0.0_',j,'.',date,'.pdf'),pgrid , width = 45, height = 50, units="in",limitsize = FALSE)
# }


library(gridExtra)

ml <- marrangeGrob(spatial_dm_list, nrow=1, ncol=1)
ggsave(paste0('Spatial_Dimplots_clusters_celltrekwithinterpolation_allcelltypes_pt.size.factor.0.8_stroke.0.0_withlegend_',j,'.',date,'.pdf'), ml, width = 15, height = 15, units="in",limitsize = FALSE)

############################################################
# Cluster composition summaries by Patient and Sample
############################################################

merged.patient.combined.integ_sub$integrated_snn_res.0.9 <- droplevels(merged.patient.combined.integ_sub$integrated_snn_res.0.9)

cluster_meta = merged.patient.combined.integ_sub@meta.data[c('Patient','Sample','integrated_snn_res.0.9')]

cluster_meta$Cluster = paste0('C',cluster_meta$integrated_snn_res.0.9)

# dotplot_items = readRDS('dotplot_items_list_march13.2025.rds')
# group_df = dotplot_items$group_df

mypercentdf_patient = as.data.frame(prop.table(table(cluster_meta$Cluster,cluster_meta$Patient),1)*100)
colnames(mypercentdf_patient) <- c("Cluster","Patient","Percentage")


library(dplyr)

df_reordered <- group_df %>% arrange(group)

mypercentdf_patient$Cluster = as.factor(mypercentdf_patient$Cluster)
levels(mypercentdf_patient$Cluster) = df_reordered$Cluster

pat_cols = color_panel_list$Patient

samp_cols = color_panel_list$Sample

# Stacked bar (100%) of patient composition per cluster (simple version)

integ_clus_vec_celltypes_integrated_snn_res.0.9_celltypes_barplot_bypatient = ggplot(mypercentdf_patient, aes(fill=Patient, y=Percentage, x=Cluster))  + geom_bar(position="fill", stat="identity",alpha = 0.9,lwd = 0.5,color = "black", linewidth = 0.05) + ggtitle(paste0('Patient distribution by Cluster : ','integrated_snn_res.0.9')) + scale_fill_manual(values=color_panel_list$Patient) + labs(x = 'integrated_snn_res.0.9', y = "%age Patient distribution by Cluster") + theme_minimal() + theme(text = element_text(size = 15),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + guides(fill=guide_legend(title="Patients")) + theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0('/n/scratch/users/s/sak4832/cluster_distances/integ_clus_vec_celltypes_integrated_snn_res.0.9_barplot_clusters_by_patient.',date,'.jpg'), integ_clus_vec_celltypes_integrated_snn_res.0.9_celltypes_barplot_bypatient, width = 20, height = 14, units="in",dpi=600)

# Add group rectangular annotations around cluster blocks
mypercentdf_patient_grp = mypercentdf_patient %>% left_join(df_reordered, by = c("Cluster"))

group_annotations <- mypercentdf_patient_grp %>%
  group_by(group) %>%
  summarize(
    xmin = min(as.numeric(Cluster)) - 0.5,  # Left boundary (bar edges)
    xmax = max(as.numeric(Cluster)) + 0.5,  # Right boundary (bar edges)
    ymin = 0,  # Bottom boundary
    ymax = 1   # Top boundary (100% stack)
  )

# Enhanced barplot with group outlines

integ_clus_vec_celltypes_integrated_snn_res.0.9_celltypes_barplot_bypatient <- ggplot(mypercentdf_patient_grp, aes(fill = Patient, y = Percentage, x = Cluster)) +
  # Bars
  geom_bar(position = "fill", stat = "identity", alpha = 0.9, lwd = 0.5, color = "black", linewidth = 0.05) +
  
  # Add rectangles around groups
  geom_rect(
    data = group_annotations,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = group),
    inherit.aes = FALSE, fill = NA, linewidth = 3,  # No fill, just outline
  ) + scale_fill_manual(
    name = "Patient",
    values = c(color_panel_list$Patient)  # Combine Patient colors and group colors
  ) +
  
  # Customize group border colors
  scale_color_manual(
    name = "Group Annotation",
    values =group_color_pal  # Add more groups as needed
  ) +
  
  # Titles and labels
  ggtitle(paste0("Patient distribution by Cluster : ", "integrated_snn_res.0.9")) +
  labs(x = "integrated_snn_res.0.9", y = "%age Patient distribution by Cluster") +
  
  # Theme adjustments
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  
  # Legend adjustments
  guides(fill = guide_legend(title = "Patients"), color = guide_legend(title = "Groups"))

# Save the plot
ggsave(
  paste0("/n/scratch/users/s/sak4832/cluster_distances/integ_clus_vec_celltypes_integrated_snn_res.0.9_barplot_clusters_by_patient.", date, ".jpg"),
  integ_clus_vec_celltypes_integrated_snn_res.0.9_celltypes_barplot_bypatient,
  width = 20, height = 14, units = "in", dpi = 600
)


mypercentdf_Sample = as.data.frame(prop.table(table(cluster_meta$Cluster,cluster_meta$Sample),1)*100)
colnames(mypercentdf_Sample) <- c("Cluster","Sample","Percentage")


library(dplyr)

df_reordered <- group_df %>% arrange(group)

mypercentdf_Sample$Cluster = as.factor(mypercentdf_Sample$Cluster)
levels(mypercentdf_Sample$Cluster) = df_reordered$Cluster

pat_cols = color_panel_list$Sample

samp_cols = color_panel_list$Sample

# integ_clus_vec_celltypes_integrated_snn_res.0.9_celltypes_barplot_bySample = ggplot(mypercentdf_Sample, aes(fill=Sample, y=Percentage, x=Cluster))  + geom_bar(position="fill", stat="identity",alpha = 0.9,lwd = 0.5,color = "black", linewidth = 0.05) + ggtitle(paste0('Sample distribution by Cluster : ','integrated_snn_res.0.9')) + scale_fill_manual(values=color_panel_list$Sample) + labs(x = 'integrated_snn_res.0.9', y = "%age Sample distribution by Cluster") + theme_minimal() + theme(text = element_text(size = 15),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + guides(fill=guide_legend(title="Samples")) + theme(plot.title = element_text(hjust = 0.5))
# ggsave(paste0('/n/scratch/users/s/sak4832/cluster_distances/integ_clus_vec_celltypes_integrated_snn_res.0.9_barplot_clusters_by_Sample.',date,'.jpg'), integ_clus_vec_celltypes_integrated_snn_res.0.9_celltypes_barplot_bySample, width = 20, height = 14, units="in",dpi=600)


mypercentdf_Sample_grp = mypercentdf_Sample %>% left_join(df_reordered, by = c("Cluster"))

group_annotations <- mypercentdf_Sample_grp %>%
  group_by(group) %>%
  summarize(
    xmin = min(as.numeric(Cluster)) - 0.5,  # Left boundary (bar edges)
    xmax = max(as.numeric(Cluster)) + 0.5,  # Right boundary (bar edges)
    ymin = 0,  # Bottom boundary
    ymax = 1   # Top boundary (100% stack)
  )


library(ggplot2)
library(dplyr)

# Barplot
library(ggplot2)
library(dplyr)

############################################################
# Composition by Sample (instead of Patient)
############################################################

sample_cols = color_panel_list$Sample
names(sample_cols) = levels(mypercentdf_Sample_grp$Sample)

integ_clus_vec_celltypes_integrated_snn_res.0.9_celltypes_barplot_bySample <- ggplot(mypercentdf_Sample_grp, aes(fill = Sample, y = Percentage, x = Cluster)) +
  # Bars
  geom_bar(position = "fill", stat = "identity", alpha = 0.9, lwd = 0.5, color = "black", linewidth = 0.05) +
  
  # Add rectangles around groups
  geom_rect(
    data = group_annotations,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = group),
    inherit.aes = FALSE, fill = NA, linewidth = 3,  # No fill, just outline
  ) + scale_fill_manual(
    name = "Sample",
    values = c(sample_cols)  # Combine Sample colors and group colors
  ) +
  
  # Customize group border colors
  scale_color_manual(
    name = "Group Annotation",
    values =group_color_pal  # Add more groups as needed
  ) +
  
  # Titles and labels
  ggtitle(paste0("Sample distribution by Cluster : ", "integrated_snn_res.0.9")) +
  labs(x = "integrated_snn_res.0.9", y = "%age Sample distribution by Cluster") +
  
  # Theme adjustments
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  
  # Legend adjustments
  guides(fill = guide_legend(title = "Samples"), color = guide_legend(title = "Groups"))

# Save the plot
ggsave(
  paste0("/n/scratch/users/s/sak4832/cluster_distances/integ_clus_vec_celltypes_integrated_snn_res.0.9_barplot_clusters_by_Sample.", date, ".jpg"),
  integ_clus_vec_celltypes_integrated_snn_res.0.9_celltypes_barplot_bySample,
  width = 20, height = 14, units = "in", dpi = 600
)

