library(colorspace)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(tidyverse)
library(sctransform)
library(ggplot2)
library(sctransform)
library(future)
library(RColorBrewer)
library(reshape2)
library(gtools)
library("CellTrek")
library(tidyr)
library(stringr)
library(dplyr)
library(Seurat)
library(yaml)
library(Ckmeans.1d.dp)


str_split_i <- function(string, pattern, i) {
  stopifnot(is.numeric(i), i == floor(i))  # check_number_whole

  if (i > 0) {
    # Return the i-th piece from the start
    out <- stringr::str_split(string, pattern, simplify = TRUE)
    return(out[, i])
  } else if (i < 0) {
    # Return the i-th piece from the end
    i <- abs(i)
    pieces <- stringr::str_split(string, pattern)
    return(purrr::map_chr(pieces, function(x) {
      n <- length(x)
      if (i > n) NA_character_ else x[[n + 1 - i]]
    }))
  } else {
    stop("i must not be 0.")
  }
}


dotplot_items = readRDS('dotplot_items_list_march13.2025.rds')

CNV_matrix_files = Sys.glob("/scratch/s/szk57/CNV_Spatial_WIP/CNV_spatial/CNV_values_matrix_t2*rds")

names(CNV_matrix_files) = str_split_i(basename(CNV_matrix_files),'\\.',2)

seurat_obj_list2 = readRDS('/scratch/s/szk57/CNV_Spatial_WIP/CNV_spatial/CNV_score_seurat_obj_list2.rds')

seurat_obj_list2_meta = lapply(names(seurat_obj_list2),function(x) {
	mymeta = seurat_obj_list2[[x]]@meta.data[c('perspot_CNV_classification','integrated_snn_res.0.9','Patient','Sample')]
	mysamp = unique(mymeta$Sample)
	mymeta$Barcodes = rownames(mymeta) 
	mymeta$Seurat_barcode = paste0(mysamp,'_',rownames(mymeta))
	mymeta$rpca_clusters = paste0('C',mymeta$`integrated_snn_res.0.9`)
	mymeta = mymeta %>% left_join(dotplot_items$group_df, by = c("rpca_clusters"="Cluster"))
	return(mymeta)
	})

seurat_obj_list2_meta_comb = do.call(rbind,seurat_obj_list2_meta)
seurat_obj_list2_meta_bypat = split(seurat_obj_list2_meta_comb,seurat_obj_list2_meta_comb$Patient)

CNV_matrices = lapply(names(CNV_matrix_files),function(x) {
	myCNVmat = readRDS(CNV_matrix_files[[x]])
	return(myCNVmat)
})

names(CNV_matrices) = names(CNV_matrix_files)

row_meta_files = Sys.glob("/scratch/s/szk57/CNV_Spatial_WIP/CNV_spatial/row_meta*rds")

names(row_meta_files) = str_split_i(basename(row_meta_files),'\\.',2)

# .B186.rds

region_df_files = Sys.glob("/scratch/s/szk57/CNV_Spatial_WIP/CNV_spatial/region_df*rds")
names(region_df_files) = str_split_i(basename(region_df_files),'\\.',2)


region_df_list = lapply(names(region_df_files),function(x) {
	myregion_df = readRDS(region_df_files[[x]])
	return(myregion_df)
})

names(region_df_list) = names(region_df_files)



row_meta_list = lapply(names(row_meta_files),function(x) {
	row_metadf = readRDS(row_meta_files[[x]])
	return(row_metadf)
})

names(row_meta_list) = names(row_meta_files)

seurat_obj_list2_meta_bypat2 = lapply(names(seurat_obj_list2_meta_bypat),function(x) {
	all_meta = merge(seurat_obj_list2_meta_bypat[[x]],row_meta_list[[x]],by=c('Seurat_barcode','Sample'))
	return(all_meta)
	})

names(seurat_obj_list2_meta_bypat2) = names(seurat_obj_list2_meta_bypat)

GBM_meta = read.table('/scratch/s/szk57/GBM_metainfo.txt',header=T,sep="\t")
GBM_meta_sub = GBM_meta[grep('^$',GBM_meta$Patient,invert=T),]
GBM_meta_sub$Patient2 = gsub('TWBK-|TWBK-|TWKI-','',GBM_meta_sub$Patient,perl=TRUE);

GBM_meta_df = GBM_meta_sub
rownames(GBM_meta_df) <- GBM_meta_df$Patient2
GBM_meta_df$Patient2 <- NULL
GBM_meta_df$Patient <- NULL
colnames(GBM_meta_df) = paste0(colnames(GBM_meta_df),'_')

GBM_meta_df[] <- lapply(GBM_meta_df, as.factor)

# Ensure all factor levels are included in the binary matrix
df_bin <- model.matrix(~ . -1, data = GBM_meta_df,
                       contrasts.arg = lapply(GBM_meta_df, contrasts, contrasts = FALSE))
df_bin_mat = as.matrix(df_bin)
Heatmap(matrix = df_bin_mat,cluster_rows = T,cluster_columns = F,column_split = NULL)


seurat_obj_list2_meta_bypat2_2 = lapply(names(seurat_obj_list2_meta_bypat2),function(x) {
	mydf = seurat_obj_list2_meta_bypat2[[x]] %>% left_join(GBM_meta_sub, by = c("Patient"='Patient2'))
	return(mydf)
	})

names(seurat_obj_list2_meta_bypat2_2) = names(seurat_obj_list2_meta_bypat2)


# lapply(names(CNV_matrices),function(x) {CNV_matrices[[x]][1:5,1:5]})

color_panel_list = readRDS('/scratch/s/szk57/color_panel_list.rds')
names(color_panel_list$Sample) = names(seurat_obj_list2)

# setdiff(unlist(lapply(names(seurat_obj_list2_meta_bypat2),function(x) {
# 	unique(seurat_obj_list2_meta_bypat[[x]]$Sample)
# 	})),names(color_panel_list$Sample))

cnv_state_colors <- c(
  `-1` = "purple4",   # deep blue — complete loss
  `-0.5` = "blue1",   # light blue — 0.5x
  `1` = "gray88",   # light gray — neutral
  `1.5` = 'hotpink',	 # pink — 1.5x
  `2` =  "red3",  	 # red — 2x
  `3` = "firebrick4"    # deep red — 3x
)

EGFR_amplification_cols = c('violetred4','grey88','grey40')
names(EGFR_amplification_cols) = c('Yes','No','UNKNOWN')

TERT_mutation_cols = c('darkgreen','greenyellow','grey88','grey40')
names(TERT_mutation_cols) = c('C228T','C250T','WT','UNKNOWN')

MGMT_methylation_cols = c('khaki4','grey88','grey40')
names(MGMT_methylation_cols) = c('Meth.','Unmeth.','UNKNOWN')

seurat_obj_list2_meta_bypat2_2 = lapply(names(seurat_obj_list2_meta_bypat2),function(x) {
	mydf = seurat_obj_list2_meta_bypat2[[x]] %>% left_join(GBM_meta_sub, by = c("Patient"='Patient2'))
	return(mydf)
	})

names(seurat_obj_list2_meta_bypat2_2) = names(seurat_obj_list2_meta_bypat2)

lapply(seurat_obj_list2_meta_bypat2_2,head,n=2)

myCNV_cols = c('gray88',brewer.pal(n=9,name='Reds')[1],brewer.pal(n=9,name='Reds')[5],brewer.pal(n=9,name='Reds')[9])
names(myCNV_cols) = c('NO_CNV','Low_CNV',"intermediate_CNV","high_CNV")
library(RColorBrewer)


final_cluster_colors = readRDS('/scratch/s/szk57/final_cluster_colors.rds')
myclus_names = paste0('C',names(final_cluster_colors))
names(final_cluster_colors) = myclus_names
# dotplot_items$group_color_palette




# infercnv_cnv_mapping_fn <- function(state_id) {
#  return(case_when(
#     state_id==1 ~ -1,
#     state_id==2 ~ -0.5,
#     state_id==3 ~ 1,
#     state_id==4 ~ 1.5,
#     state_id==5 ~ 2,
#     state_id==6 ~ 3,
#   ))
# }



counter = 0
hm_list = lapply(names(seurat_obj_list2_meta_bypat2_2),function(x) {
	counter=counter+1
	if(counter==1){
		my_show_column_names = TRUE
		my_show_heatmap_legend = TRUE
	}else{
		my_show_column_names = FALSE
		my_show_heatmap_legend = FALSE
	}
	CNV_values_matrix_t2 = CNV_matrices[[x]]
	subcluster_df_seurat2 = seurat_obj_list2_meta_bypat2_2[[x]]
	region_df = region_df_list[[x]]
	all(rownames(CNV_values_matrix_t2) %in% subcluster_df_seurat2$Seurat_barcode)
	index <- match(rownames(CNV_values_matrix_t2), subcluster_df_seurat2$Seurat_barcode)
	row_meta <- subcluster_df_seurat2[index, ]
	row_annot <- rowAnnotation(
	Patient = row_meta$Patient,
	TERT_mutation= row_meta$TERT,
	MGMT_methylation = row_meta$MGMT,
	EGFR_amplification = row_meta$EGFR.Amp,
	Sample = row_meta$Sample,
	CNV_classification=row_meta$perspot_CNV_classification,
	Domains = row_meta$group,
	Niches = row_meta$rpca_clusters,
	# CNV_metaclusters = row_meta$metacluster_id,
	# CNV_subclusters = row_meta$subcluster_id,
	col = list(
	Patient = color_panel_list$Patient,
	TERT_mutation= TERT_mutation_cols,
	MGMT_methylation = MGMT_methylation_cols,
	EGFR_amplification = EGFR_amplification_cols,
	Sample = color_panel_list$Sample,
	CNV_classification= myCNV_cols,
	Domains = dotplot_items$group_color_palette,
	Niches = final_cluster_colors
	# CNV_metaclusters = metacluster_cols
	# CNV_subclusters = subcluster_cols
	),    annotation_name_gp = gpar(fontsize = 1),  # Effectively hides annotation titles
    show_annotation_name = FALSE,show_legend = TRUE,
	)

	CNV_values_matrix_t2_2 <- as.matrix(apply(CNV_values_matrix_t2, 2, as.character))

	myrow_split <- factor(row_meta$metacluster_id, levels = mixedsort(unique(row_meta$metacluster_id)))
	mycol_split = factor(region_df$Chr, levels = mixedsort(unique(region_df$Chr)))
	
	hm = Heatmap(CNV_values_matrix_t2_2,
	                           cluster_columns = F,
	                           show_column_dend = FALSE,
	                           cluster_column_slices = F,
	                           row_names_gp = gpar(fontsize = 10),
	                           column_names_gp = gpar(fontsize = 12),
	                           column_title_gp = gpar(fontsize = 15),
	                           column_gap = unit(0.1, "mm"),
                               heatmap_width = unit(20, "cm"),
							   heatmap_height = unit(10, "cm"),
	                           cluster_rows = F,
	                           left_annotation = row_annot,
	                           # top_annotation = col_annot,
	                           show_column_names = my_show_column_names,
	                           column_split =mycol_split ,row_split= myrow_split,
	                           cluster_row_slices = TRUE,
	                           show_row_dend = FALSE,
	                           col = cnv_state_colors,
	                           row_title_gp = gpar(fontsize = 15),
	                           show_row_names = F,
	                           row_names_side = "right",
	                           column_title_rot = 90,
	                           column_names_rot = 90,
	                           raster_by_magick=FALSE,raster_device="agg_png",
	                           use_raster=FALSE,
	                           heatmap_legend_param = list(
	    at = names(cnv_state_colors),
	    labels = c("0x ", "0.5x", "1x", "1.5x", "2x", "3x"),
	    title = "CNV state",
	    legend_gp = gpar(fontsize = 8)
	  ))

	return(hm)

	})

names(hm_list) = names(seurat_obj_list2_meta_bypat2_2)

patchwork_list = lapply(hm_list,function(x) {grid.grabExpr(draw(x,show_annotation_legend = FALSE))})

combined <- wrap_plots(patchwork_list, ncol = 1,nrow = 9)

jpeg("CNV_values_matrix_all_patients.jpg", width = 20, height = 35, units = "in", res = 300)
print(combined)  # or plot(combined)
dev.off()


##############################################################################################################
##############################################################################################################
##############################################################################################################
layer_fun = function(j, i, x, y, w, h, fill) {
target_idx <- which(fill != cnv_state_colors[["1"]])
grid.points(x[target_idx], y[target_idx],
            pch = 16, size = unit(1.8, "mm"),
            gp = gpar(col = fill[target_idx]))
}

dir.create("CNV_heatmap_jpegs", showWarnings = FALSE)

for (i in seq_along(seurat_obj_list2_meta_bypat2_2)) {
  x <- names(seurat_obj_list2_meta_bypat2_2)[i]

  CNV_values_matrix_t2 <- CNV_matrices[[x]]
  subcluster_df_seurat2 <- seurat_obj_list2_meta_bypat2_2[[x]]
  region_df <- region_df_list[[x]]
  index <- match(rownames(CNV_values_matrix_t2), subcluster_df_seurat2$Seurat_barcode)
  row_meta <- subcluster_df_seurat2[index, ]

  row_annot <- rowAnnotation(
    Patient = row_meta$Patient,
    TERT_mutation = row_meta$TERT,
    MGMT_methylation = row_meta$MGMT,
    EGFR_amplification = row_meta$EGFR.Amp,
    Sample = row_meta$Sample,
    CNV_classification = row_meta$perspot_CNV_classification,
    Domains = row_meta$group,
    Niches = row_meta$rpca_clusters,
    col = list(
      Patient = color_panel_list$Patient,
      TERT_mutation = TERT_mutation_cols,
      MGMT_methylation = MGMT_methylation_cols,
      EGFR_amplification = EGFR_amplification_cols,
      Sample = color_panel_list$Sample,
      CNV_classification = myCNV_cols,
      Domains = dotplot_items$group_color_palette,
      Niches = final_cluster_colors
    ),
    annotation_name_gp = gpar(fontsize = 6),
    show_annotation_name = FALSE
  )

  CNV_values_matrix_t2_2 <- as.matrix(apply(CNV_values_matrix_t2, 2, as.character))
  myrow_split <- factor(row_meta$metacluster_id, levels = mixedsort(unique(row_meta$metacluster_id)))
  mycol_split <- factor(region_df$Chr, levels = mixedsort(unique(region_df$Chr)))

  legend_flag <- (i == 1)

  hm <- Heatmap(
    CNV_values_matrix_t2_2,
    name = "CNV state",
    col = cnv_state_colors,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    cluster_column_slices = FALSE,
    cluster_row_slices = TRUE,
    column_split = mycol_split,
    row_split = myrow_split,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 8),
    column_title_gp = gpar(fontsize = 10),
    row_title_gp = gpar(fontsize = 8),
    left_annotation = row_annot,
    heatmap_width = unit(14, "in"),
    heatmap_height = unit(8, "in"),
    raster_device = "agg_png",
    raster_by_magick = FALSE,
    use_raster = FALSE,
    heatmap_legend_param = list(
      at = names(cnv_state_colors),
      labels = c("0x", "0.5x", "1x", "1.5x", "2x", "3x"),
      title = "CNV state",
      legend_gp = gpar(fontsize = 7)
    ),
	layer_fun = function(j, i, x, y, w, h, fill) {
	target_idx <- which(fill != cnv_state_colors[["1"]])
	grid.points(x[target_idx], y[target_idx],
	            pch = 16, size = unit(1.8, "mm"),
	            gp = gpar(col = fill[target_idx]))
	}
  )

  jpeg(
    file = paste0("CNV_heatmap_jpegs/", x, ".jpg"),
    width = 15, height =10 , res = 300, units = "in"
  )
  draw(hm, show_annotation_legend = legend_flag)
  dev.off()
}

##############################################################################################################
#####################Unannotated#########################################################################################
##############################################################################################################

cell_fun = function(j, i, x, y, width, height, fill) {
  grid.points(x = x, y = y, pch = 16, size = unit(2.5, "mm"), gp = gpar(col = fill))
}


for (i in seq_along(seurat_obj_list2_meta_bypat2_2)) {
  x <- names(seurat_obj_list2_meta_bypat2_2)[i]

  CNV_values_matrix_t2 <- CNV_matrices[[x]]
  subcluster_df_seurat2 <- seurat_obj_list2_meta_bypat2_2[[x]]
  region_df <- region_df_list[[x]]
  index <- match(rownames(CNV_values_matrix_t2), subcluster_df_seurat2$Seurat_barcode)
  row_meta <- subcluster_df_seurat2[index, ]

  row_annot <- rowAnnotation(
    Patient = row_meta$Patient,
    TERT_mutation = row_meta$TERT,
    MGMT_methylation = row_meta$MGMT,
    EGFR_amplification = row_meta$EGFR.Amp,
    Sample = row_meta$Sample,
    CNV_classification = row_meta$perspot_CNV_classification,
    Domains = row_meta$group,
    Niches = row_meta$rpca_clusters,
    col = list(
      Patient = color_panel_list$Patient,
      TERT_mutation = TERT_mutation_cols,
      MGMT_methylation = MGMT_methylation_cols,
      EGFR_amplification = EGFR_amplification_cols,
      Sample = color_panel_list$Sample,
      CNV_classification = myCNV_cols,
      Domains = dotplot_items$group_color_palette,
      Niches = final_cluster_colors
    ),
    annotation_name_gp = gpar(fontsize = 6),
    show_annotation_name = FALSE,show_legend=FALSE
  )

  CNV_values_matrix_t2_2 <- as.matrix(apply(CNV_values_matrix_t2, 2, as.character))
  myrow_split <- factor(row_meta$metacluster_id, levels = mixedsort(unique(row_meta$metacluster_id)))
  mycol_split <- factor(region_df$Chr, levels = mixedsort(unique(region_df$Chr)))

  legend_flag <- (i == 1)

  hm <- Heatmap(
    CNV_values_matrix_t2_2,
    name = "CNV state",
    col = cnv_state_colors,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    cluster_column_slices = FALSE,
    cluster_row_slices = FALSE,
    column_split = mycol_split,
    row_split = myrow_split,
    show_row_dend = FALSE,
    column_title = NULL,
	row_title=NULL,
    show_column_dend = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 0),
    column_title_gp = gpar(fontsize = 0),
    row_title_gp = gpar(fontsize = 0),
    left_annotation = row_annot,
    heatmap_width = unit(14, "in"),
    heatmap_height = unit(8, "in"),
    raster_device = "agg_png",
    raster_by_magick = FALSE,
    use_raster = FALSE,show_heatmap_legend = FALSE,
    heatmap_legend_param = list(
      at = names(cnv_state_colors),
      labels = c("0x", "0.5x", "1x", "1.5x", "2x", "3x"),
      title = "CNV state",
      legend_gp = gpar(fontsize = 7)
    ),
	layer_fun = function(j, i, x, y, w, h, fill) {
	target_idx <- which(fill != cnv_state_colors[["1"]])
	grid.points(x[target_idx], y[target_idx],
	            pch = 16, size = unit(1.8, "mm"),
	            gp = gpar(col = fill[target_idx]))
	}
  )

  jpeg(
    file = paste0("CNV_heatmap_jpegs/", x, "_unannotated.jpg"),
    width = 15, height =10, res = 600, units = "in"
  )
  draw(hm, show_annotation_legend = FALSE,show_heatmap_legend = FALSE)
  dev.off()
}

# docker container run -it --name o2r_4.3.1 0b8bab3fedec bash
