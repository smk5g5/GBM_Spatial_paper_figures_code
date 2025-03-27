
###Sample map visium to visium image####
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

mysample_map$image_name[mysample_map$sample_name %in% 'WU1221-Core2'] = 'WU1221.Core2'
mysample_map$image_name[mysample_map$sample_name %in% 'WU1221-Edge1'] = 'WU1221.Edge1'
##########################################

###############Classify spots into senders/receivers###########################
sigspots_senders_to_receiverspot_list <- function(sigspots_file, neighbourhood_spots_file, seurat_obj) {
  # Load the neighbourhood spots data
  neighbourhood_spots_df <- read.table(neighbourhood_spots_file, sep = " ", header = TRUE)

  # Create list to store sender spots and their neighborhoods
  sender_spot_neighborhood_list <- list()
  for (i in 1:nrow(neighbourhood_spots_df)) {
    sender_spot <- neighbourhood_spots_df$X[i]
    sender_spot_neighborhood_list[[sender_spot]] <- str_split(neighbourhood_spots_df$neighbour_bcs[i], ',')[[1]]
  }

  # Load significant spots data
  sigspots_df <- read.csv(sigspots_file, sep = ",", header = TRUE)
  sigspots_df_list <- split(sigspots_df, sigspots_df$LR_interaction)

  # Initialize list to store results
  LR_sender_receiver_spot_list <- list()

  # Iterate over each ligand-receptor interaction
  for (LR in names(sigspots_df_list)) {
    myligand <- str_split(LR, '_')[[1]][1]
    myrec <- str_split(LR, '_')[[1]][2]

    # Identify barcodes with expression for ligand and receptor
	# For version 1        receptor_barcodes <- names(which(seurat_obj@assays$Spatial@counts[myrec, ] > 0)) #to further filter neighborhood barcodes to only include higher expression barcodes for receptors

    #version 3 recepter expression > 2 and ligand expression > 2

    receptor_barcodes <- names(which(seurat_obj@assays$Spatial@counts[myrec, ] >= 1)) #to further filter neighborhood barcodes to only include higher expression barcodes for receptors
    ligand_barcodes <- names(which(seurat_obj@assays$Spatial@counts[myligand, ] >= 1))

    # Filter significant sender spots based on ligand expression
    mysig_sender_spots <- intersect(sigspots_df_list[[LR]]$Significant_barcodes, ligand_barcodes)

	# Version 3 is only use spots for a particular LR if they are in Significant_barcodes consider those which are 
    # Find receiver barcodes in the neighborhood of sender spots where receptor is expressed
    #version 3 recepter expression > 2 and ligand expression > 2


	# Version 2 is only use spots for a particular LR if they are in Significant_barcodes consider those which are 
    # Find receiver barcodes in the neighborhood of sender spots where receptor is expressed

    #Version 1 and version 3
    # receptor_barcodes_in_neighbourhood <- lapply(sender_spot_neighborhood_list[mysig_sender_spots], function(x) {
    #   intersect(receptor_barcodes, x)
    # })
    # Remove empty elements

    # Version 2 use receptors spots only if they are in mysig_sender_spots
    receptor_barcodes_in_neighbourhood <- lapply(sender_spot_neighborhood_list[mysig_sender_spots], function(x) {
      abc = intersect(receptor_barcodes, x)
      return(intersect(abc,mysig_sender_spots))
    })

    receptor_barcodes_in_neighbourhood <- Filter(function(x) length(x) > 0, receptor_barcodes_in_neighbourhood)

    # Skip if no valid sender spots or all neighborhoods are empty
    if (length(mysig_sender_spots) == 0 || length(receptor_barcodes_in_neighbourhood) == 0) {
      next
    }

    # Create data frames for each sender-receiver pair
    Sender_spot_receiver_spot_df_list <- lapply(names(receptor_barcodes_in_neighbourhood), function(x) { 
      mydf <- data.frame(receiver_spot = receptor_barcodes_in_neighbourhood[[x]])
      mydf$sender_spot <- x
      mydf$LR <- LR
      return(mydf[c('sender_spot', 'receiver_spot', 'LR')])
    })

    # Combine data frames into one for each LR interaction
    Sender_spot_receiver_spot_df <- do.call(rbind, Sender_spot_receiver_spot_df_list)
    LR_sender_receiver_spot_list[[LR]] <- Sender_spot_receiver_spot_df
  }

  # Combine all LR interaction data frames into a single data frame
  LR_sender_receiver_spot_df <- do.call(rbind, LR_sender_receiver_spot_list)

  # Add sample metadata to LR_sender_receiver_spot_df
  mysample_seurat_meta <- seurat_obj@meta.data[c('integrated_snn_res.0.9')]
  mysample_seurat_meta$barcode <- rownames(mysample_seurat_meta)

  mysample_seurat_meta$integrated_snn_res.0.9 = droplevels(mysample_seurat_meta$integrated_snn_res.0.9)

  LR_sender_receiver_spot_df$Sample <- unique(seurat_obj$Sample)
  LR_sender_receiver_spot_df$Patient <- unique(seurat_obj$Patient)

  # Join metadata and rename columns correctly
  LR_sender_receiver_spot_df2 <- LR_sender_receiver_spot_df %>%
    left_join(mysample_seurat_meta, by = c("sender_spot" = "barcode")) %>%
    rename(sender_cluster = integrated_snn_res.0.9)
  
  LR_sender_receiver_spot_df3 <- LR_sender_receiver_spot_df2 %>%
    left_join(mysample_seurat_meta, by = c("receiver_spot" = "barcode")) %>%
    rename(receiver_cluster = integrated_snn_res.0.9)

  # Return the final data frame
  return(LR_sender_receiver_spot_df3)
}

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


##### Modify single cell reference celltype annotations to more broad annotations######
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


############designate sender spots and receiver spots from significant spots############
# sender_spot_neighborhood_list[mysig_sender_spots]

sigspots_senders_to_receiverspot_list <- function(sigspots_file, neighbourhood_spots_file, seurat_obj) {
  # Load the neighbourhood spots data
  neighbourhood_spots_df <- read.table(neighbourhood_spots_file, sep = " ", header = TRUE)

  # Create list to store sender spots and their neighborhoods
  sender_spot_neighborhood_list <- list()
  for (i in 1:nrow(neighbourhood_spots_df)) {
    sender_spot <- neighbourhood_spots_df$X[i]
    sender_spot_neighborhood_list[[sender_spot]] <- str_split(neighbourhood_spots_df$neighbour_bcs[i], ',')[[1]]
  }

  # Load significant spots data
  sigspots_df <- read.csv(sigspots_file, sep = ",", header = TRUE)
  sigspots_df_list <- split(sigspots_df, sigspots_df$LR_interaction)

  # Initialize list to store results
  LR_sender_receiver_spot_list <- list()

  # Iterate over each ligand-receptor interaction
  for (LR in names(sigspots_df_list)) {
    myligand <- str_split(LR, '_')[[1]][1]
    myrec <- str_split(LR, '_')[[1]][2]

    # Identify barcodes with expression for ligand and receptor
	# For version 1        receptor_barcodes <- names(which(seurat_obj@assays$Spatial@counts[myrec, ] > 0)) #to further filter neighborhood barcodes to only include higher expression barcodes for receptors

    #version 3 recepter expression > 2 and ligand expression > 2

    receptor_barcodes <- names(which(seurat_obj@assays$Spatial@counts[myrec, ] >= 2)) #to further filter neighborhood barcodes to only include higher expression barcodes for receptors
    ligand_barcodes <- names(which(seurat_obj@assays$Spatial@counts[myligand, ] >= 2))

    # Filter significant sender spots based on ligand expression
    mysig_sender_spots <- intersect(sigspots_df_list[[LR]]$Significant_barcodes, ligand_barcodes)

	# Version 3 is only use spots for a particular LR if they are in Significant_barcodes consider those which are 
    # Find receiver barcodes in the neighborhood of sender spots where receptor is expressed
    #version 3 recepter expression > 2 and ligand expression > 2


	# Version 2 is only use spots for a particular LR if they are in Significant_barcodes consider those which are 
    # Find receiver barcodes in the neighborhood of sender spots where receptor is expressed

    #Version 1 and version 3
    # receptor_barcodes_in_neighbourhood <- lapply(sender_spot_neighborhood_list[mysig_sender_spots], function(x) {
    #   intersect(receptor_barcodes, x)
    # })
    # Remove empty elements

    # Version 2 use receptors spots only if they are in mysig_sender_spots
    receptor_barcodes_in_neighbourhood <- lapply(sender_spot_neighborhood_list[mysig_sender_spots], function(x) {
      abc = intersect(receptor_barcodes, x)
      return(intersect(abc,mysig_sender_spots))
    })

    receptor_barcodes_in_neighbourhood <- Filter(function(x) length(x) > 0, receptor_barcodes_in_neighbourhood)

    # Skip if no valid sender spots or all neighborhoods are empty
    if (length(mysig_sender_spots) == 0 || length(receptor_barcodes_in_neighbourhood) == 0) {
      next
    }

    # Create data frames for each sender-receiver pair
    Sender_spot_receiver_spot_df_list <- lapply(names(receptor_barcodes_in_neighbourhood), function(x) { 
      mydf <- data.frame(receiver_spot = receptor_barcodes_in_neighbourhood[[x]])
      mydf$sender_spot <- x
      mydf$LR <- LR
      return(mydf[c('sender_spot', 'receiver_spot', 'LR')])
    })

    # Combine data frames into one for each LR interaction
    Sender_spot_receiver_spot_df <- do.call(rbind, Sender_spot_receiver_spot_df_list)
    LR_sender_receiver_spot_list[[LR]] <- Sender_spot_receiver_spot_df
  }

  # Combine all LR interaction data frames into a single data frame
  LR_sender_receiver_spot_df <- do.call(rbind, LR_sender_receiver_spot_list)

  # Add sample metadata to LR_sender_receiver_spot_df
  mysample_seurat_meta <- seurat_obj@meta.data[c('integrated_snn_res.0.9')]
  mysample_seurat_meta$barcode <- rownames(mysample_seurat_meta)

  mysample_seurat_meta$integrated_snn_res.0.9 = droplevels(mysample_seurat_meta$integrated_snn_res.0.9)

  LR_sender_receiver_spot_df$Sample <- unique(seurat_obj$Sample)
  LR_sender_receiver_spot_df$Patient <- unique(seurat_obj$Patient)

  # Join metadata and rename columns correctly
  LR_sender_receiver_spot_df2 <- LR_sender_receiver_spot_df %>%
    left_join(mysample_seurat_meta, by = c("sender_spot" = "barcode")) %>%
    rename(sender_cluster = integrated_snn_res.0.9)
  
  LR_sender_receiver_spot_df3 <- LR_sender_receiver_spot_df2 %>%
    left_join(mysample_seurat_meta, by = c("receiver_spot" = "barcode")) %>%
    rename(receiver_cluster = integrated_snn_res.0.9)

  # Return the final data frame
  return(LR_sender_receiver_spot_df3)
}
################################################################################################





