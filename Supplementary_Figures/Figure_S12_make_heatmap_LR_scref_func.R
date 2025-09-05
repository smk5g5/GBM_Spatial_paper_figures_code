#######################################################################
# Function: make_heatmap_LR_scref
#
#   Generate a heatmap visualization of ligand-receptor (LR) communication 
#   probabilities for selected groups and interactions, along with a binary 
#   reference matrix heatmap for comparison.

# Description : Fig. S12. Ligand receptor interactions in single-cell RNA-seq (reference data) recapitulate findings in Xenium and Visium datasets.
# (A) LR interactions where Tumor cells act as senders and other cell types as receivers, grouped by broad categories indicated in the top annotation.
# Each row corresponds to an LR pair, and each column represents a sender–receiver combination as inferred from the single-cell reference, which includes
# spiked-in normal tissue data to improve representation of non-tumor cell types. Cell color in the main heatmap indicates the strength of communication (communication probability),
# and asterisks denote the statistical significance of the interaction. The binary heatmap on the far right of each panel indicates whether each LR interaction was also observed in
# the Visium and Xenium datasets (defined as present in at least three patients per sender–receiver group for Xenium, and present in at least three patients in Visium as shown in 
# Fig. 4, S10, and S11) and in 10x single-cell RNA-seq (defined as having a communication probability greater than the overall mean). (B) LR interactions where non-tumor cell types
# act as senders and tumor cells as receivers. (C) LR interactions between Myeloid and Lymphoid cells. These results highlight consistent intercellular communication patterns across 
# single-cell and spatial transcriptomic platforms and provide insights into potential biological mechanisms.

# Returns:
#   - A list containing:
#       1. heatmap_matrix: ordered scaled probability matrix
#       2. ranked_LRs_bygroup: data frame of selected LRs
#       3. htmap_obj: Heatmap object
#######################################################################


make_heatmap_LR_scref <- function(selected_groups,outfile,mygroup_order,selected_interactions,interaction_list) {
  
  # if(length(selected_interactions)==0){
  # mysel_df =subset(combined_df, interaction_name %in% unique(ranked_LRs$interaction_name) & sender_grp_to_receiver_grp %in% selected_groups) %>% data.frame()
  # }else{
  #   mysel_df =subset(combined_df, interaction_name %in% selected_interactions & sender_grp_to_receiver_grp %in% selected_groups) %>% data.frame()
  # }
  ############################################################
  # Filter input dataframe for selected LR interactions
  ############################################################


  combined_df2 = subset(combined_df,interaction_name %in% selected_interactions)

  # Ligand-receptors (LRs) in selected groups
  ranked_LRs_abovemeanprob_bygroup <- subset(combined_df2, 
                                             sender_grp_to_receiver_grp %in% selected_groups) %>%
    distinct(interaction_name) %>%
    data.frame()
  
  ranked_LRs_bygroup <- subset(combined_df2, 
                               sender_grp_to_receiver_grp %in% selected_groups) %>%
    distinct(interaction_name) %>%
    data.frame()
  
  # interaction_list[['10X-single cell RNAseq']] = as.character(unique(ranked_LRs_abovemeanprob_bygroup$interaction_name))

  # Convert interaction list to binary matrix

  interaction_list_matrix = list_to_matrix(interaction_list)
  
  # Union of all interaction names (from both assays)
  #This is not necessary and a placeholder
  #initially we used only LRs which were present
  # above mean communication probability in 10x single cell reference

  all_interactions = union(as.character(unique(ranked_LRs_abovemeanprob_bygroup$interaction_name)),as.character(unique(ranked_LRs_bygroup$interaction_name)))
  
  # Subset combined_df2 with selected groups and interactions

  mysel_df =subset(combined_df2, interaction_name %in% all_interactions & sender_grp_to_receiver_grp %in% selected_groups) %>% data.frame()
  
  # Create scaled probability and p-value matrices

  scaled_prob_mat <- mysel_df %>%
    select(interaction_name, source.target, scaled_prob) %>%
    pivot_wider(names_from = source.target, values_from = scaled_prob)  %>%
    column_to_rownames("interaction_name") %>%
    as.matrix()
  
  pval_mat <- mysel_df %>%
    select(interaction_name, source.target, pval) %>%
    pivot_wider(names_from = source.target, values_from = pval) %>%
    column_to_rownames("interaction_name") %>%
    as.matrix()
  
  common_interactions <- intersect(rownames(interaction_list_matrix), rownames(scaled_prob_mat))
  
  # Reorder interaction_list_matrix rows to match scaled_prob_mat
  interaction_list_matrix_ordered <- interaction_list_matrix[rownames(scaled_prob_mat)[rownames(scaled_prob_mat) %in% common_interactions], , drop = FALSE]
  
  colors = structure(c("gray99", '#556B2F'), names = c('0','1') ) # black, red, green, blue
  
  # Significance asterisk matrix (based on p-values)

  asterisk_mat <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat),
                         dimnames = dimnames(pval_mat))
  asterisk_mat[!is.na(pval_mat) & pval_mat < 0.001] <- "*"

  # Column grouping and ordering

  col_split <- mysel_df %>%
    distinct(source.target, sender_grp_to_receiver_grp) %>%
    filter(source.target %in% colnames(scaled_prob_mat)) %>%
    arrange(match(source.target, colnames(scaled_prob_mat)))  # Match column order
  
  col_split$split_var = factor(col_split$sender_grp_to_receiver_grp,levels=mygroup_order)
  
  col_split_sorted <- col_split %>%
    group_by(split_var) %>%
    arrange(mixedorder(source.target), .by_group = TRUE) %>%
    ungroup()
  
  scaled_prob_mat_ord = scaled_prob_mat[,col_split_sorted$source.target]
  asterisk_mat_ord = asterisk_mat[,col_split_sorted$source.target]
  
  # Replace missing values with 0

  scaled_prob_mat_ord[is.na(scaled_prob_mat_ord)] <- 0

    # --- Plot the heatmap ---
    nr <- nrow(scaled_prob_mat_ord)
    nc <- ncol(scaled_prob_mat_ord)

  # Adjust rownames (replace "_" with ":")

    mysel_rownames = gsub('_',':',rownames(scaled_prob_mat_ord))

    rownames(scaled_prob_mat_ord) = mysel_rownames
    rownames(asterisk_mat_ord) = mysel_rownames
    rownames(interaction_list_matrix_ordered) = mysel_rownames

    # Create a tibble with height and width determined by row/col logic
    
    dim_tibble <- tibble(
    height = case_when(
      nr <= 10 ~ 12,
      nr > 10 & nr <= 20 ~ 20,
      nr > 20 & nr <= 40 ~ 34,
      nr > 40 & nr <= 60 ~ 48,
      nr > 60 & nr <= 70 ~ 55,
      nr > 70 & nr <= 90 ~ 59,
      nr > 90 & nr <= 100 ~ 65,
    ),
    width = case_when(
      nc <= 15 ~ 30,
      nc > 15 & nc <= 25 ~ 45,
      nc > 25 & nc <= 30 ~ 55,
      nc > 30 & nc <= 40 ~ 65,
      nc > 40 & nc <= 60 ~ 75,
      nc > 60 & nc <= 80 ~ 85,
      nc > 80 & nc <= 95 ~ 95,
      nc > 95 & nc <= 105 ~ 100,
    )
    )

    # Access values
    heatmap_height <- dim_tibble$height
    heatmap_width <- dim_tibble$width

  ht_opt$TITLE_PADDING = unit(c(8.5, 8.5), "points")

  # Heatmap 1: Scaled ligand-receptor probabilities

  myht = Heatmap(scaled_prob_mat_ord,
                 name = "scaled ligand-receptor communication probability",
                 col = col_fun,
                 column_split = col_split_sorted$split_var,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.rect(x = x, y = y, width = width, height = height, 
                             gp = gpar(col = "black", fill = NA, lwd = 0.5))
                   if (!is.na(scaled_prob_mat_ord[i, j])) {
                     grid.text(asterisk_mat_ord[i, j], x = x, y = y, gp = gpar(fontsize = 20))
                   }
                 },show_row_names = FALSE,
                 cluster_rows = TRUE,row_title_gp = gpar(fontsize = 10),column_title_gp = gpar(fill = celltype_palette[mygroup_order],fontsize = 20,col = "white"),
                 cluster_columns = FALSE,column_names_gp = gpar(fontsize = 15),
                 row_names_gp = gpar(fontsize = 15),
                 column_names_rot = 90)
  
  # Heatmap 2: Binary interaction reference matrix

  binhtmap2 = Heatmap(interaction_list_matrix_ordered, name = c("binary matrix\n technique LR threshold"),
                      show_column_dend = FALSE,cluster_columns = TRUE,cluster_rows = TRUE,col = colors,rect_gp = gpar(col = "black", lwd = 1),
                      row_names_gp = gpar(fontsize = 15),
                      column_names_gp = gpar(fontsize = 15),
                      column_title_gp = gpar(fontsize = 12),
                      column_gap = unit(0.5, "mm"),
                      show_column_names = T,width = unit(3, "cm"),
                      # column_split = factor(x = merged.combined_sub@meta.data[[integ]], levels =levels(merged.combined_sub@meta.data[[integ]])),
                      show_row_dend = FALSE,
                      heatmap_legend_param = list(title ="binary matrix LR by technique", legend_gp = gpar(fontsize = 15),legend_width = unit(5, "cm"),labels_gp = gpar(fontsize = 12)),
                      row_title_gp = gpar(fontsize = 20),
                      show_row_names = T,
                      row_names_side = "right",
                      column_title_rot = 90,
                      column_names_rot = 90)

  # Combine heatmaps

  ht_list = myht+binhtmap2
  mysize =calc_ht_size(draw(ht_list),unit = "cm")
  
  print(paste0('heatmap height : ',heatmap_height))
  print(paste0('heatmap width : ',heatmap_width))
  
  jpeg(outfile, width = heatmap_width, height = heatmap_height, units = "cm", res = 600)
  print(ht_list)  # or plot(combined)
  dev.off()

  jpeg(paste0('Nolegends_',outfile), width = heatmap_width, height = heatmap_height, units = "cm", res = 600)
  draw(ht_list, show_heatmap_legend = FALSE, show_annotation_legend = FALSE)
  dev.off()

  # Return objects
  return(list(heatmap_matrix=scaled_prob_mat_ord,ranked_LRs_abovemeanprob_bygroup=ranked_LRs_abovemeanprob_bygroup,htmap_obj=myht))
}
