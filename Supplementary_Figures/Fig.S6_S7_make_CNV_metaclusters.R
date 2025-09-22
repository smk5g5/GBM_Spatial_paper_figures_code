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
library(reshape2)
library(gtools)
library("CellTrek")

library(tidyr)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(GenomicRanges)
library(colorspace)
library(circlize)
library(infercnv)
library(doParallel)
library(foreach)
library(EnrichedHeatmap)
library(patchwork)

set.seed(1234)


date = gsub("2025-","25",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

seurat_obj_list2 = readRDS('/n/scratch/users/s/sak4832/CNV_spatial/CNV_score_seurat_obj_list2.rds')

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  args <- c("--help")
}

sample_id = as.character(args[1])

infercnv_cnv_mapping_fn <- function(state_id) {
 return(case_when(
    state_id==1 ~ -1,
    state_id==2 ~ -0.5,
    state_id==3 ~ 1,
    state_id==4 ~ 1.5,
    state_id==5 ~ 2,
    state_id==6 ~ 3,
  ))
}

infercnv_dir = str_interp("/n/scratch/users/s/sak4832/CNV_spatial/${sample_id}_infercnv2_clonepubversion_0.01_leiden")
infercnv_path = str_interp("${infercnv_dir}/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.3.pred_cnv_regions.dat")
mcmc_path = str_interp("${infercnv_dir}/BayesNetOutput.HMMi6.leiden.hmm_mode-subclusters/MCMC_inferCNV_obj.rds")

mcmc = readRDS(mcmc_path)

CNV_sim_path = str_interp("/n/scratch/users/s/sak4832/CNV_spatial/${sample_id}_infercnv2_clonepubversion_0.01_leiden/CNV_pearson_similarity_by_spot_${sample_id}.rds")

CNV_similarity = readRDS(CNV_sim_path)

test_df <- map_dfr(
  mcmc@tumor_subclusters$subclusters$malignant,
  ~ as.data.frame(.x),
  .id = "subcluster_id"  # this will add the list name as a new column
)

test_df$`.x` = NULL

test_df$barcode_names = rownames(test_df)

test_df$grp = str_split_i(test_df$barcode_names,'_',2)
test_df$sample_barcode = str_split_i(test_df$barcode_names,'_',1)
test_df_lst = split(test_df,test_df$grp)


mysamps = grep(sample_id,names(seurat_obj_list2),value=T)

seurat_obj_list2_sub = seurat_obj_list2[mysamps]

mymerged_df_spl = lapply(names(seurat_obj_list2_sub), function(x) {
mydf = data.frame(sample_barcode=Cells(seurat_obj_list2_sub[[x]]),Sample=x)
return(mydf)
})

names(mymerged_df_spl) = names(seurat_obj_list2_sub)

mymerged_df = do.call(rbind,mymerged_df_spl)

sample_grp_map <- lapply(split(mymerged_df, mymerged_df$Sample), function(df_sample) {
  sample_barcodes <- unique(df_sample$sample_barcode)

  # For each grp, check if all sample_barcodes are in that grp
  matching_grps <- sapply(split(test_df, test_df$grp), function(df_test) {
    all(sample_barcodes %in% df_test$sample_barcode)
  })

  matched_grp_names <- names(matching_grps)[matching_grps]

  if (length(matched_grp_names) == 1) {
    return(matched_grp_names)
  } else {
    message("No match or multiple matches for Sample: ", unique(df_sample$Sample))
    print(sample_barcodes)
    return(NA)
  }
})


# Convert to data.frame
sample_grp_map_df <- data.frame(
  Sample = names(sample_grp_map),
  grp = unlist(sample_grp_map)
)


subcluster_df_seurat = do.call(rbind,lapply(names(mymerged_df_spl),function(x) {
	mygrp = sample_grp_map[[x]]
	my_merged_pat = merge(test_df_lst[[mygrp]],mymerged_df_spl[[x]],by='sample_barcode')
	return(my_merged_pat)}
	))

compute_subcluster_cnv_similarity <- function(spot_clusters,spot_sim) {
  # # Step 1: Map CNV states to numeric values (fast, vectorized)
  # state_to_value <- c(`1` = -1, `2` = -0.5, `3` = 1, `4` = 1.5, `5` = 2, `6` = 3)
  # cnv_vector <- as.character(cnv_state_matrix)
  # cnv_mapped_vector <- state_to_value[cnv_vector]
  # cnv_mapped <- matrix(
  #   as.numeric(cnv_mapped_vector),
  #   nrow = nrow(cnv_state_matrix),
  #   dimnames = dimnames(cnv_state_matrix)
  # )

  # # Step 2: Compute spot × spot Pearson similarity
  # spot_sim <- cor(t(cnv_mapped), method = "pearson", use = "pairwise.complete.obs")

  # spot_sim = readRDS('/n/scratch/users/s/sak4832/CNV_spatial/B172_infercnv2_clonepubversion_0.01_leiden/CNV_pearson_similarity_by_spot_B172.rds')

  # Step 3: Assign subcluster per spot
  spot_to_cluster <- spot_clusters$subcluster_id
  names(spot_to_cluster) <- spot_clusters$barcode_names

  # Filter to matching barcodes only
  spot_to_cluster <- spot_to_cluster[rownames(spot_sim)]
  unique_clusters <- unique(spot_to_cluster)

  # Step 4: Compute mean spot similarity between subclusters
  subcluster_sim <- matrix(NA, length(unique_clusters), length(unique_clusters),
                           dimnames = list(unique_clusters, unique_clusters))

  for (i in seq_along(unique_clusters)) {
    for (j in seq_along(unique_clusters)) {
      spots_i <- names(spot_to_cluster[spot_to_cluster == unique_clusters[i]])
      spots_j <- names(spot_to_cluster[spot_to_cluster == unique_clusters[j]])

      values <- spot_sim[spots_i, spots_j, drop = FALSE]
      subcluster_sim[i, j] <- mean(values, na.rm = TRUE)
    }
  }

  return(subcluster_sim)
}

# ------------------------------------------------------------------------------
# Function: cluster_subclusters_with_min_spots
# Purpose : Cluster subclusters using a similarity matrix while ensuring:
#           - no more than `max_clusters` total clusters
#           - each final cluster contains at least `min_size` spots
#           - sequential metacluster naming (e.g., meta_1, meta_2, ..., meta_k)
#
# Args:
#   subcluster_sim     : Symmetric subcluster × subcluster similarity matrix
#   spot_clusters      : Data frame with columns:
#                          - subcluster_id
#                          - barcode_names (spot IDs)
#   max_clusters       : Maximum number of metaclusters allowed
#   min_size           : Minimum number of spots per metacluster
#
# Returns:
#   A data frame with:
#     - subcluster_id
#     - metacluster_id (merged & relabeled cluster ID)
# ------------------------------------------------------------------------------

cluster_subclusters_with_min_spots <- function(subcluster_sim, spot_clusters, max_clusters = 20, min_size = 10) {
  library(dplyr)

  # Step 1: Hierarchical clustering on 1 - similarity
  dist_mat <- as.dist(1 - subcluster_sim)
  hc <- hclust(dist_mat, method = "ward.D2")
  meta_clusters <- cutree(hc, k = max_clusters)

  # Initial cluster assignment
  meta_clusters_df <- data.frame(
    subcluster_id = names(meta_clusters),
    metacluster_id = paste0("meta_", unname(meta_clusters)),
    stringsAsFactors = FALSE
  )

  # Map subclusters to their spots
  spot_map <- spot_clusters %>%
    select(subcluster_id, barcode_names)

  # Step 2: Compute spot counts per metacluster
  repeat {
    # Join subcluster-to-metacluster and map to spots
    merged <- left_join(spot_map, meta_clusters_df, by = "subcluster_id")

    # Count spots per metacluster
    cluster_sizes <- merged %>%
      group_by(metacluster_id) %>%
      summarise(size = n(), .groups = "drop")

    # Identify underpopulated metaclusters
    small_clusters <- cluster_sizes %>%
      filter(size < min_size) %>%
      arrange(size)

    if (nrow(small_clusters) == 0) break  # All clusters valid

    # Merge smallest invalid cluster into nearest large cluster
    sc <- small_clusters$metacluster_id[1]
    subcl_in_sc <- meta_clusters_df$subcluster_id[meta_clusters_df$metacluster_id == sc]

    candidates <- setdiff(rownames(subcluster_sim), subcl_in_sc)
    merge_target <- NULL

    for (s in subcl_in_sc) {
      nearest <- names(sort(1 - subcluster_sim[s, candidates]))[1]
      merge_target <- meta_clusters_df$metacluster_id[meta_clusters_df$subcluster_id == nearest]
      if (!is.null(merge_target)) break
    }

    meta_clusters_df$metacluster_id[meta_clusters_df$metacluster_id == sc] <- merge_target
  }

  # Step 3: Relabel metaclusters sequentially (meta_1, meta_2, ...)
  unique_metas <- sort(unique(meta_clusters_df$metacluster_id))
  meta_relabel_map <- setNames(paste0("meta_", seq_along(unique_metas)), unique_metas)
  meta_clusters_df$metacluster_id <- meta_relabel_map[meta_clusters_df$metacluster_id]

  return(meta_clusters_df)
}

# subcluster_sim1 = compute_subcluster_cnv_similarity(spot_clusters=test_df_merged,spot_sim=CNV_similarity)
# meta_clusters_df <- cluster_subclusters_with_min_spots(subcluster_sim=subcluster_sim1, spot_clusters=test_df_merged, max_clusters = 10, min_size = 10)


subcluster_sim1 = compute_subcluster_cnv_similarity(spot_clusters=subcluster_df_seurat,spot_sim=CNV_similarity)
meta_clusters_df <- cluster_subclusters_with_min_spots(subcluster_sim=subcluster_sim1, spot_clusters=subcluster_df_seurat, max_clusters = 20, min_size = 10)

subcluster_df_seurat2 <- subcluster_df_seurat %>%
    left_join(meta_clusters_df, by = c("subcluster_id"))

CNV_per_spot = function(i){

cnv_reg = as.character(mcmc@cell_gene[[i]]$cnv_regions)

all_tumor_spots = colnames(mcmc@expr.data)[mcmc@observation_grouped_cell_indices$malignant]

mystate = mcmc@cell_gene[[i]]$State
myspots = intersect(all_tumor_spots,colnames(mcmc@expr.data)[mcmc@cell_gene[[i]]$Cells])

# observation_grouped_cell_indices = mcmc@reference_grouped_cell_indices$normal
observation_grouped_cell_indices = mcmc@observation_grouped_cell_indices$malignant
# myspots = intersect(colnames(mcmc@expr.data)[observation_grouped_cell_indices],myspots)

other_spots = setdiff(colnames(mcmc@expr.data)[observation_grouped_cell_indices],myspots)

other_matrix =  matrix(rep(3 ,length(other_spots)), nrow = 1)
rownames(other_matrix) = cnv_reg
colnames(other_matrix) = other_spots


my_matrix = matrix(rep(mcmc@cell_gene[[i]]$State ,length(myspots)), nrow = 1)
rownames(my_matrix) = cnv_reg
colnames(my_matrix) = myspots

comb_mat = cbind(other_matrix,my_matrix)

return(comb_mat)
}


HMM_matrix_CNV_states_per_spot <- function(mcmc_obj) {

my_CNV_scores <- lapply(1:length(mcmc_obj@cell_gene), CNV_per_spot)
reference_order <- colnames(my_CNV_scores[[1]])

my_CNV_scores_ordered <- lapply(my_CNV_scores, function(mat) {
  mat[, reference_order, drop = FALSE]  # drop = FALSE preserves matrix structure
})

my_CNV_scores_ordered_comb = do.call(rbind,my_CNV_scores_ordered)

return(my_CNV_scores_ordered_comb)
}

HMM_matrix_CNV_values_per_spot <- function(mcmc_obj) {
  # Mapping from HMM state to numeric CNV value
  state_to_value <- c(`1` = -1, `2` = -0.5, `3` = 1, `4` = 1.5, `5` = 2, `6` = 3)

  # Step 1: Build the CNV state matrix per spot
  my_CNV_scores <- lapply(1:length(mcmc_obj@cell_gene), CNV_per_spot)
  reference_order <- colnames(my_CNV_scores[[1]])

  my_CNV_scores_ordered <- lapply(my_CNV_scores, function(mat) {
    mat[, reference_order, drop = FALSE]
  })

  # Step 2: Combine into full matrix
  CNV_state_matrix <- do.call(rbind, my_CNV_scores_ordered)

  # Step 3: Convert HMM state matrix to numeric CNV values
  CNV_value_matrix <- matrix(
    data = state_to_value[as.character(CNV_state_matrix)],
    nrow = nrow(CNV_state_matrix),
    ncol = ncol(CNV_state_matrix),
    dimnames = dimnames(CNV_state_matrix)
  )

  return(CNV_value_matrix)
}


CNV_state_matrix = HMM_matrix_CNV_states_per_spot(mcmc_obj=mcmc)

CNV_values_matrix = HMM_matrix_CNV_values_per_spot(mcmc_obj=mcmc)

#CNV_state_matrix_t = t(CNV_state_matrix)
CNV_values_matrix_t = t(CNV_values_matrix)

subcluster_df_seurat2$Seurat_barcode = paste0(subcluster_df_seurat2$Sample,'_',subcluster_df_seurat2$sample_barcode)

# Ensure CNV_values_matrix_t has barcode_names as rownames
barcodes <- rownames(CNV_values_matrix_t)

# Create named vector: names = barcode_names, values = Seurat_barcode
barcode_to_seurat <- setNames(subcluster_df_seurat2$Seurat_barcode, subcluster_df_seurat2$barcode_names)

# Map new rownames
new_rownames <- barcode_to_seurat[barcodes]

# Assign if all barcodes have a match
if (all(!is.na(new_rownames))) {
  rownames(CNV_values_matrix_t) <- new_rownames
} else {
  warning("Some barcodes could not be matched to Seurat_barcode.")
}

region_df = data.frame(Chr=str_split_i(mixedsort(colnames(CNV_values_matrix_t)),'-',1),region=mixedsort(colnames(CNV_values_matrix_t)))

CNV_values_matrix_t2 <- CNV_values_matrix_t[, mixedsort(colnames(CNV_values_matrix_t))]

library(ComplexHeatmap)

color_panel_list = readRDS('/n/scratch/users/s/sak4832/color_panel_list.rds')

index = match(rownames(CNV_values_matrix_t2),subcluster_df_seurat2$Seurat_barcode)
metacluster_ord = subcluster_df_seurat2$metacluster_id[index]
subcluster_ord = subcluster_df_seurat2$subcluster_id[index]
Sample_ord = subcluster_df_seurat2$Sample[index]

patcols = color_panel_list$Patient[grep(sample_id, names(color_panel_list$Patient), value = TRUE)]

metacluster_cols <- generate_similar_colors(patcols, n = length(unique(metacluster_ord)))
names(metacluster_cols) = mixedsort(unique(metacluster_ord))

subcluster_cols <- generate_similar_colors(patcols, n = length(unique(subcluster_ord)))
names(subcluster_cols) = mixedsort(unique(subcluster_ord))

sample_cols = color_panel_list$Sample[grep(sample_id,names(color_panel_list$Sample))]

library(circlize)


cnv_state_colors <- c(
  `-1` = "#313695",   # deep blue — complete loss
  `-0.5` = "#74add1",   # light blue — 0.5x
  `1` = "#f0f0f0",   # light gray — neutral
  `1.5` = "#e41a1c",   # red — 1.5x
  `2` = "#67000d"    # deep red — 2x
)

chromosome_cols = DiscretePalette_scCustomize(num_colors = length(unique(region_df$Chr)), palette = "polychrome")
names(chromosome_cols) = unique(region_df$Chr)

all(rownames(CNV_values_matrix_t2) %in% subcluster_df_seurat2$Seurat_barcode)
index <- match(rownames(CNV_values_matrix_t2), subcluster_df_seurat2$Seurat_barcode)

row_meta <- subcluster_df_seurat2[index, ]

row_annot <- rowAnnotation(
  Sample = row_meta$Sample,
  CNV_metaclusters = row_meta$metacluster_id,
  CNV_subclusters = row_meta$subcluster_id,
  col = list(
    Sample = sample_cols,
    CNV_metaclusters = metacluster_cols,
    CNV_subclusters = subcluster_cols
  )
)

myrow_split <- factor(row_meta$metacluster_id, levels = mixedsort(unique(row_meta$metacluster_id)))
mycol_split = factor(region_df$Chr, levels = mixedsort(unique(region_df$Chr)))

hm = Heatmap(CNV_values_matrix_t2, name = "CNV HMM matrix",
                           cluster_columns = F,
                           show_column_dend = FALSE,
                           cluster_column_slices = F,
                           row_names_gp = gpar(fontsize = 10),
                           column_names_gp = gpar(fontsize = 12),
                           column_title_gp = gpar(fontsize = 15),
                           column_gap = unit(0.5, "mm"),
                           cluster_rows = F,
                           left_annotation = row_annot,
                           # top_annotation = col_annot,
                           show_column_names = F,
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
                           use_raster=TRUE,
                           heatmap_legend_param = list(
    at = names(cnv_state_colors),
    labels = c("0x (del)", "0.5x", "1x", "1.5x", "2x"),
    title = "CNV state",
    legend_gp = gpar(fontsize = 8)
  ))

jpeg(paste0('CNV_values_matrix_t2_test_',sample_id,'.jpg'),width = 50,height = 50,units="cm", res=600)
draw(hm,show_annotation_legend = TRUE, show_heatmap_legend = TRUE)
dev.off()

saveRDS(region_df,paste0('region_df.',sample_id,'.rds'))
saveRDS(row_meta,paste0('row_meta.',sample_id,'.rds'))
saveRDS(CNV_values_matrix_t2,paste0('CNV_values_matrix_t2.',sample_id,'.rds'))
