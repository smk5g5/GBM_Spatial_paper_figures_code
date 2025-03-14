sender_celltype_htmap_receiver_celltype_htmap <- function(sc,rc,interaction_df,mylrs) {

  myLR_spotlevel = subset(cluster_distances_with_LR_comb_grouped,sender_group==sc&receiver_group==rc)
  myLR_spotlevel_sub = myLR_spotlevel %>% subset(LR %in% mylrs)

	ligand_df = unique(myLR_spotlevel_sub[c('barcode1','sender_spot2','sender_group','LR','Sample')])
	receptor_df = unique(myLR_spotlevel_sub[c('barcode2','receiver_spot2','receiver_group','LR','Sample')])

	rownames(receptor_df) = NULL
	rownames(ligand_df) = NULL

	ligand_df_list = split(ligand_df,ligand_df$Sample)
	receptor_df_list = split(receptor_df,receptor_df$Sample)

	ligand_df_celltrekinfo = lapply(names(ligand_df_list), function(x) {

	   mydf_merged = merge(ligand_df_list[[x]],spots_to_cells_celltrek[[x]]$broad_celltypes_per_spot,by.x='barcode1',by.y='spot_id')

	  mydf_merged$Glia = rowSums(mydf_merged[intersect(c("Astrocyte", "Committed oligodendrocyte precursor", "Oligodendrocyte", "Oligodendrocyte precursor"), colnames(mydf_merged))])
	  mydf_merged$Lymphoid = rowSums(mydf_merged[intersect(c("Bcells", "CD4_Tcells", "CD8_Tcells", "NK/NK-like", "Other_Tcells"), colnames(mydf_merged))])
	  mydf_merged$Myeloid = rowSums(mydf_merged[intersect(c("DC", "Microglia", "Monocyte", "TAMs"), colnames(mydf_merged))])
	  mydf_merged$Vascular = rowSums(mydf_merged[intersect(c("Fibroblast", "Vascular"), colnames(mydf_merged))])

	  print(head(mydf_merged))

	    sender_cluster_celltypecount <- mydf_merged %>%
	    group_by(LR) %>%
	    summarize(
	      M1 = sum(M1),
	      M2 = sum(M2),
	      M3 = sum(M3),
	      M4 = sum(M4),
	      M5 = sum(M5),
	      Neuron = sum(Neuron),
	      Glia = sum(Glia),
	      Lymphoid = sum(Lymphoid),
	      Myeloid = sum(Myeloid),
	      Vascular = sum(Vascular)
	    )

	  })

	receptor_df_celltrekinfo = lapply(names(receptor_df_list), function(x) {

	mydf_merged = merge(receptor_df_list[[x]],spots_to_cells_celltrek[[x]]$broad_celltypes_per_spot,by.x='barcode2',by.y='spot_id')

	mydf_merged$Glia = rowSums(mydf_merged[intersect(c("Astrocyte", "Committed oligodendrocyte precursor", "Oligodendrocyte", "Oligodendrocyte precursor"), colnames(mydf_merged))])
	mydf_merged$Lymphoid = rowSums(mydf_merged[intersect(c("Bcells", "CD4_Tcells", "CD8_Tcells", "NK/NK-like", "Other_Tcells"), colnames(mydf_merged))])
	mydf_merged$Myeloid = rowSums(mydf_merged[intersect(c("DC", "Microglia", "Monocyte", "TAMs"), colnames(mydf_merged))])
	mydf_merged$Vascular = rowSums(mydf_merged[intersect(c("Fibroblast", "Vascular"), colnames(mydf_merged))])

	print(head(mydf_merged))

	receiver_cluster_celltypecount <- mydf_merged %>%
	group_by(LR) %>%
	summarize(
	  M1 = sum(M1),
	  M2 = sum(M2),
	  M3 = sum(M3),
	  M4 = sum(M4),
	  M5 = sum(M5),
	  Neuron = sum(Neuron),
	  Glia = sum(Glia),
	  Lymphoid = sum(Lymphoid),
	  Myeloid = sum(Myeloid),
	  Vascular = sum(Vascular)
	)

	})


  ligand_df_celltrekinfo2 = do.call(rbind,ligand_df_celltrekinfo)  %>% group_by(LR)  %>% summarize(
      M1 = sum(M1),
      M2 = sum(M2),
      M3 = sum(M3),
      M4 = sum(M4),
      M5 = sum(M5),
      Neuron = sum(Neuron),
      Glia = sum(Glia),
      Lymphoid = sum(Lymphoid),
      Myeloid = sum(Myeloid),
      Vascular = sum(Vascular)
    )

  receptor_df_celltrekinfo2 = do.call(rbind,receptor_df_celltrekinfo)  %>% group_by(LR)  %>% summarize(
      M1 = sum(M1),
      M2 = sum(M2),
      M3 = sum(M3),
      M4 = sum(M4),
      M5 = sum(M5),
      Neuron = sum(Neuron),
      Glia = sum(Glia),
      Lymphoid = sum(Lymphoid),
      Myeloid = sum(Myeloid),
      Vascular = sum(Vascular)
    )


	# Ensure both ligand and receptor data frames have the same LR pairs
	common_LRs <- intersect(ligand_df_celltrekinfo2$LR, receptor_df_celltrekinfo2$LR)

	ligand_df_celltrekinfo2 <- ligand_df_celltrekinfo2 %>% filter(LR %in% common_LRs)
	receptor_df_celltrekinfo2 <- receptor_df_celltrekinfo2 %>% filter(LR %in% common_LRs)

	myrows = ligand_df_celltrekinfo2$LR
	ligand_df_celltrekinfo2$LR = NULL
	ligand_df_celltrekinfo2_mat = as.matrix(ligand_df_celltrekinfo2)
	rownames(ligand_df_celltrekinfo2_mat) = myrows

	myrows = receptor_df_celltrekinfo2$LR
	receptor_df_celltrekinfo2$LR = NULL
	receptor_df_celltrekinfo2_mat = as.matrix(receptor_df_celltrekinfo2)
	rownames(receptor_df_celltrekinfo2_mat) = myrows

	receptor_df_celltrekinfo2_mat_allvals = rbind(receptor_df_celltrekinfo2_mat,NA_matrix_ret(myrows=setdiff(mylrs,rownames(receptor_df_celltrekinfo2_mat)),mycols=colnames(receptor_df_celltrekinfo2)))

	ligand_df_celltrekinfo2_mat_allvals = rbind(ligand_df_celltrekinfo2_mat,NA_matrix_ret(myrows=setdiff(mylrs,rownames(ligand_df_celltrekinfo2_mat)),mycols=colnames(ligand_df_celltrekinfo2_mat)))
return(list(ligand_df_celltrekinfo=ligand_df_celltrekinfo2_mat_allvals,receptor_df_celltrekinfo=receptor_df_celltrekinfo2_mat_allvals))
}




NA_matrix_ret <- function(myrows,mycols) {
NA_matrix <- matrix(NA, nrow = length(myrows), ncol = length(mycols))
colnames(NA_matrix) <- mycols
rownames(NA_matrix) = myrows
return(NA_matrix)
}
