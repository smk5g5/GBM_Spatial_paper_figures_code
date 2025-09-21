##############################
# Load Required Libraries
##############################
library("dplyr")              # Data manipulation (filter, mutate, join, etc.)
library("Seurat")             # Single-cell analysis framework
library("viridis")            # Color palettes for plots
library("RColorBrewer")       # Additional color palettes
library("ggplot2")            # General plotting system
library(SingleCellExperiment) # S4 class for single-cell data
library(scran)                # Single-cell RNA-seq analysis tools
library(rcartocolor)          # CARTO color palettes
library(stringr)              # String manipulation functions
library(scCustomize)          # Helper functions to customize Seurat plots
library(gtools)               # Miscellaneous tools (e.g. smart ordering)
library(rjson)                # JSON reading and parsing
library(reshape2)

##############################

date = gsub("2025-","25",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

##############################
# Pre-computed File Lists
##############################
# External shell commands to generate file lists:
# find $PWD -name 'scalefactors_json.json' > /n/scratch/users/s/sak4832/scalefactors_json_inp.txt
# find $PWD -name 'tissue_positions_list.csv' > /n/scratch/users/s/sak4832/tissue_positions_list.csv


scalefactors_json_inp = read.table('/n/scratch/users/s/sak4832/scalefactors_json_inp.txt',header=F,sep=" ")

tissue_list_positions_inp = read.table('/n/scratch/users/s/sak4832/tissue_list_positions_inp.txt',header=F,sep=" ")

##############################
# Extract Sample Names from Paths
##############################
# Derive sample names by stripping path and removing prefixes "TWBK-" or "TWKI-"


tissue_list_positions_inp$sample_name = gsub("TWBK-|TWKI-","",str_split_i(tissue_list_positions_inp$V1,'/',-3),perl=TRUE);

scalefactors_json_inp$sample_name = gsub("TWBK-|TWKI-","",str_split_i(scalefactors_json_inp$V1,'/',-3),perl=TRUE);


##############################
# Define Sample → Image Name Mapping
##############################

# Correct mappings for specific samples (ensure uniqueness / match Visium naming conventions)

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

##############################
# Merge File Lists with Sample Map
##############################
scalefactors_json_inp2 = merge(scalefactors_json_inp,mysample_map,by='sample_name')

tissue_list_positions_inp2 = merge(tissue_list_positions_inp,mysample_map,by='sample_name')

##############################
# Load Tissue Positions Data
##############################

tissue_list_positions_list = list()

for(i in 1:nrow(tissue_list_positions_inp2)){
	sample_name = tissue_list_positions_inp2$image_name[i]
	tissue_list_positions_list[[sample_name]] = read.csv(tissue_list_positions_inp2$V1[i],header=F)
	colnames(tissue_list_positions_list[[sample_name]]) = c("barcode","in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres")
}

##############################
# Load Scale Factors JSON Data
##############################

scalefactors_json_list = list()

for(i in 1:nrow(scalefactors_json_inp2)){
	sample_name = scalefactors_json_inp2$image_name[i]
	scalefactors_json_list[[sample_name]] = fromJSON(paste(readLines(scalefactors_json_inp2$V1[i]), collapse=""))
}

##############################
# Compute Image Coordinates (row/col)
##############################

for(i in names(tissue_list_positions_list)){
	tissue_list_positions_list[[i]]$imagerow = tissue_list_positions_list[[i]]$pxl_row_in_fullres * 65/scalefactors_json_list[[i]][['spot_diameter_fullres']]
	tissue_list_positions_list[[i]]$imagecol = tissue_list_positions_list[[i]]$pxl_col_in_fullres * 65/scalefactors_json_list[[i]][['spot_diameter_fullres']]
}

##############################
# Set Row Names and Remove Barcode Column
##############################

for(i in names(tissue_list_positions_list)){
rownames(tissue_list_positions_list[[i]]) = tissue_list_positions_list[[i]]$barcode
tissue_list_positions_list[[i]]$barcode = NULL
}

##############################
# Compute Distance Matrices for Each Sample
##############################

distance_matrices_by_sample = lapply(names(tissue_list_positions_list),function(x) {
		dist = as.matrix(dist(tissue_list_positions_list[[x]][,c('imagerow','imagecol')]))
    # Pairwise Euclidean distances between spots in image coordinates
		spot.dist = min(dist[upper.tri(dist,diag = FALSE)])     # Compute minimum inter-spot distance (could be used for normalization)
		return(dist)
	})

names(distance_matrices_by_sample) = names(tissue_list_positions_list)

# Load Seurat Object List
seurat_obj_list2 = readRDS('/n/scratch/users/s/sak4832/seurat_obj_list2.rds')


##############################
# Annotate Distance Matrices with RPCA integrated Cluster Information
##############################


meltdf_distance_clusterinfo2 = lapply(names(distance_matrices_by_sample),function(x) {
	meltdf = melt(distance_matrices_by_sample[[x]])
	colnames(meltdf) = c('barcode1','barcode2','distance')   # Melt the distance matrix into long format (barcode1, barcode2, distance)

	meltdf_sub = subset(meltdf,barcode1!=barcode2)   # Remove self-distances (where barcode1 == barcode2)

  # Extract cluster assignments from Seurat object metadata
  # using the "integrated_snn_res.0.9" clustering resolution
	integ_clus_df = seurat_obj_list2[[x]]@meta.data['integrated_snn_res.0.9']
	integ_clus_df$barcode = rownames(integ_clus_df)

	# print(head(integ_clus_df))
	# Keep only those rows where barcode1 exists in Seurat metadata

	meltdf_sub2 = subset(meltdf_sub,meltdf_sub$barcode1 %in% integ_clus_df$barcode)
	# Further filter to keep only rows where barcode2 also exists in metadata

	meltdf_sub22 = subset(meltdf_sub2,meltdf_sub2$barcode2 %in% integ_clus_df$barcode)
	
	# Join cluster identity of barcode1 into the melted distance dataframe

	result2 = meltdf_sub22 %>% left_join(integ_clus_df, by = c("barcode1"="barcode"))

	# print(head(result2))
	# Store cluster identity of barcode1 under "Cluster_1"

	result2$Cluster_1 = result2$`integrated_snn_res.0.9`
	result2$`integrated_snn_res.0.9` = NULL
	# Join cluster identity of barcode2 into the dataframe
	result2 = result2 %>% left_join(integ_clus_df, by = c("barcode2"="barcode"))

	result2$Cluster_2 = result2$`integrated_snn_res.0.9`
	result2$`integrated_snn_res.0.9` = NULL

	# print(head(result2))
	return(result2)   # Return the annotated dataframe for this sample
})

names(meltdf_distance_clusterinfo2) = names(distance_matrices_by_sample)

##############################
# Identify Adjacent Cluster Barcodes
##############################
adjacent_cluster_barcodes = lapply(names(meltdf_distance_clusterinfo2), function(x) {
  
  # Get all unique barcode–barcode–distance rows for this sample
  # (optionally could restrict to Cluster_1 != Cluster_2, but here we keep all)
  mydist = unique(meltdf_distance_clusterinfo2[[x]])
  
  # Extract unique cluster IDs from both cluster columns
  c1 = droplevels(c(unique(mydist$Cluster_1)))
  c2 = droplevels(c(unique(mydist$Cluster_2)))
  
  # Create a unified sorted cluster list
  myclusters = mixedsort(unique(as.character(c1), as.character(c2)))
  
  # Generate all pairwise cluster combinations (inter-cluster)
  cluster_combn = apply(combn(myclusters, 2), 2, paste, collapse = '_')
  
  # Add within-cluster pairs (e.g. Cluster1_Cluster1) alongside inter-cluster pairs
  mycluster_combn2 = union(
    unlist(lapply(myclusters, function(x) { paste0(x, '_', x) })),
    cluster_combn
  )
  
  # For each cluster pair combination, subset the distance dataframe
  unique_cluster_pairdf = lapply(mycluster_combn2, function(x) {
    cluster1 = str_split_i(x, '_', 1)
    cluster2 = str_split_i(x, '_', 2)
    print(cluster1)
    print(cluster2)
    return(subset(mydist, Cluster_1 == cluster1 & Cluster_2 == cluster2))
  })
  
  # Name each subsetted dataframe by its cluster combination (inter-cluster only)
  names(unique_cluster_pairdf) = cluster_combn
  
  # Collapse the list of per-cluster pair dfs into a single dataframe
  unique_cluster_pairdf2 = do.call(rbind, unique_cluster_pairdf)
  rownames(unique_cluster_pairdf2) = NULL
  
  # Keep only pairs of barcodes whose distance is ≤ 102 (adjacent spots)
  adjacent_clusters_only = subset(unique_cluster_pairdf2, distance <= 102)
  
  # Group by cluster pair to prepare per-pair adjacency lists
  min_dist_by_sample = adjacent_clusters_only %>%
    group_by(Cluster_1, Cluster_2)
  
  # Split grouped data into a list, one element per cluster pair
  min_dist_by_sample_list = group_split(min_dist_by_sample)
  
  # Create clean names (e.g. "2_3") for each cluster pair list element
  mycluster_combn = lapply(min_dist_by_sample_list, function(x) {
    myc1 = unique(x$Cluster_1)
    myc2 = unique(x$Cluster_2)
    return(paste0(myc1, '_', myc2))
  })
  
  names(min_dist_by_sample_list) = unlist(mycluster_combn)
  
  # Return the adjacency list for this sample
  return(min_dist_by_sample_list)
})

# Assign sample names to the adjacency list
names(adjacent_cluster_barcodes) = names(meltdf_distance_clusterinfo2)

##############################
# Count Adjacent Cluster Pairs by Sample
##############################
adjacent_cluster_barcodes_count = lapply(names(meltdf_distance_clusterinfo2), function(x) {
  
  # Get all unique barcode–barcode–distance rows for this sample
  # (optionally could restrict to Cluster_1 != Cluster_2, but here we keep all)
  mydist = unique(meltdf_distance_clusterinfo2[[x]])
  
  # Extract cluster IDs observed in both columns
  c1 = droplevels(c(unique(mydist$Cluster_1)))
  c2 = droplevels(c(unique(mydist$Cluster_2)))
  
  # Create a unified sorted cluster list
  myclusters = mixedsort(unique(as.character(c1), as.character(c2)))
  
  # Generate all inter-cluster pair combinations
  cluster_combn = apply(combn(myclusters, 2), 2, paste, collapse = '_')
  
  # Add within-cluster combinations (e.g. "Cluster1_Cluster1")
  mycluster_combn2 = union(
    unlist(lapply(myclusters, function(x) { paste0(x, '_', x) })),
    cluster_combn
  )
  
  # For each cluster pair combination, subset rows from the distance table
  unique_cluster_pairdf = lapply(mycluster_combn2, function(x) {
    cluster1 = str_split_i(x, '_', 1)
    cluster2 = str_split_i(x, '_', 2)
    print(cluster1)  # Debug print (can be removed later)
    print(cluster2)  # Debug print (can be removed later)
    return(subset(mydist, Cluster_1 == cluster1 & Cluster_2 == cluster2))
  })
  
  # Name each subsetted dataframe by its cluster combination
  names(unique_cluster_pairdf) = cluster_combn
  
  # Merge all per-cluster-pair dataframes into one
  unique_cluster_pairdf2 = do.call(rbind, unique_cluster_pairdf)
  rownames(unique_cluster_pairdf2) = NULL
  
  # Keep only "adjacent" spot pairs (≤ 102 pixels apart)
  adjacent_clusters_only = subset(unique_cluster_pairdf2, distance <= 102)
  
  # Count number of adjacent barcode pairs per cluster combination
  min_dist_by_sample_count = adjacent_clusters_only %>%
    group_by(Cluster_1, Cluster_2) %>%
    summarize(count = n())
  
  # Return adjacency counts dataframe for this sample
  return(min_dist_by_sample_count)
})

# Assign sample names to adjacency count list
names(adjacent_cluster_barcodes_count) = names(meltdf_distance_clusterinfo2)

##############################
# Add Sample Names to Adjacency Count Dataframes
##############################
adjacent_cluster_barcodes_count2 = lapply(names(adjacent_cluster_barcodes_count),function(x) {
	adjacent_cluster_barcodes_count[[x]]$sample_name = x   # Add a column "sample_name" to each per-sample adjacency count dataframe

	return(adjacent_cluster_barcodes_count[[x]])   # Return the updated dataframe
	})

names(adjacent_cluster_barcodes_count2) = names(adjacent_cluster_barcodes_count)

##############################
# Remove Within-Cluster Adjacencies and Add Sample Names
##############################

# From adjacency counts, remove rows where Cluster_1 == Cluster_2
adjacent_cluster_barcodes_count2_removesameclusters = lapply(
  names(adjacent_cluster_barcodes_count), 
  function(x) {
    return(subset(adjacent_cluster_barcodes_count[[x]], Cluster_1 != Cluster_2))
  }
)

names(adjacent_cluster_barcodes_count2_removesameclusters) = names(adjacent_cluster_barcodes_count)

# Add a "sample_name" column to each filtered dataframe
adjacent_cluster_barcodes_count2_removesameclusters2 = lapply(
  names(adjacent_cluster_barcodes_count2_removesameclusters),
  function(x) {
    adjacent_cluster_barcodes_count2_removesameclusters[[x]]$sample_name = x
    return(adjacent_cluster_barcodes_count2_removesameclusters[[x]])
  }
)

names(adjacent_cluster_barcodes_count2_removesameclusters2) = names(adjacent_cluster_barcodes_count2_removesameclusters)

adjacent_cluster_barcodes_count_merged = do.call(rbind,adjacent_cluster_barcodes_count2_removesameclusters2)
