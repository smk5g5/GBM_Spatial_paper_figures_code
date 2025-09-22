# ============================================================================
# Figure 3G code for cluster adjacency analysis
# Author: [Saad Khan]
# Description:  Spatial proximity map of niches,
 # where edge width is proportional to Weighted Mean Adjacency 
 # and edge color is proportional to the fraction of samples in which each niche pair is adjacent.
# ============================================================================
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

saveRDS(adjacent_cluster_barcodes,'adjacent_cluster_barcodes.rds')

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

##############################
# Identify Per-Spot Neighbours (within 102 px)
##############################
per_spot_neighbours = lapply(names(meltdf_distance_clusterinfo2), function(x) {
  
  # Filter for spot–spot pairs where distance ≤ 102 pixels
  adjacent_spots_only = subset(meltdf_distance_clusterinfo2[[x]], distance <= 102)
  
  # Split by the "barcode1" column
  # -> creates a list where each element contains all neighbours of one spot
  adjacent_spot_list = split(adjacent_spots_only, adjacent_spots_only$barcode1)
  
  # Extract only the "barcode2" column (neighbour barcodes) for each spot
  neighbour_barcodes = lapply(
    names(adjacent_spot_list),
    function(x) {
      return(adjacent_spot_list[[x]]$barcode2)
    }
  )
  
  # Name each neighbour list element by the corresponding spot barcode
  names(neighbour_barcodes) = names(adjacent_spot_list)
  
  # Return the per-spot neighbour list for this sample
  return(neighbour_barcodes)
})

# Assign sample names to the top-level list
names(per_spot_neighbours) = names(meltdf_distance_clusterinfo2)

##############################
# Compute Per-Cluster Probability per Sample
##############################
per_cluster_probability = lapply(names(seurat_obj_list2), function(x) {
  
  # Extract cluster assignments (resolution 0.9) from Seurat metadata
  integ_clus_df = seurat_obj_list2[[x]]@meta.data['integrated_snn_res.0.9']
  
  # Add barcodes as an explicit column
  integ_clus_df$barcode = rownames(integ_clus_df)
  
  #Count number of spots/cells per cluster
  cluster_count_by_sample = integ_clus_df %>%
    group_by(integrated_snn_res.0.9) %>%
    summarize(cluster_count = n())
  
  # Step 4: Compute cluster probability = (# cells in cluster) / (total cells in sample)
  cluster_count_by_sample$per_cluster_probability =
    cluster_count_by_sample$cluster_count / nrow(integ_clus_df)
  
  # Return dataframe with cluster counts and probabilities
  return(cluster_count_by_sample)
})

# Assign sample names as list element names
names(per_cluster_probability) = names(seurat_obj_list2)

##############################
# Build Cluster-Annotated Spot Adjacency Networks
##############################
adjacency_by_sample = lapply(names(per_spot_neighbours), function(x) {
  
  # Construct adjacency pairs from per-spot neighbours
  # For each spot (Source), create rows mapping to its neighbour barcodes (Target)
  adjacency_pairs <- do.call(
    rbind,
    lapply(names(per_spot_neighbours[[x]]), function(j) {
      data.frame(Source = j, Target = per_spot_neighbours[[x]][[j]])
    })
  )
  
  #Ensure adjacency is treated as undirected
  # (so Source=1, Target=2 is the same as Source=2, Target=1)
  adjacency_pairs <- adjacency_pairs %>%
    rowwise() %>%
    mutate(Pair = paste(sort(c(Source, Target)), collapse = "_")) %>%
    ungroup()
  
  # Remove duplicate undirected pairs
  adjacency_pairs_unique = adjacency_pairs[!duplicated(adjacency_pairs$Pair), ]
  
  #Extract clustering info for this sample
  integ_clus_df = seurat_obj_list2[[x]]@meta.data['integrated_snn_res.0.9']
  integ_clus_df$barcode = rownames(integ_clus_df)
  colnames(integ_clus_df) = c('rpca_cluster', 'barcode')
  
  #Annotate adjacency pairs with cluster identity of Source barcode
  adjacency_pairs_unique2 = adjacency_pairs_unique %>%
    left_join(integ_clus_df, by = c("Source" = "barcode"))
  adjacency_pairs_unique2$Cluster1 = adjacency_pairs_unique2$rpca_cluster
  adjacency_pairs_unique2$rpca_cluster = NULL
  
  # Annotate adjacency pairs with cluster identity of Target barcode
  adjacency_pairs_unique22 = adjacency_pairs_unique2 %>%
    left_join(integ_clus_df, by = c("Target" = "barcode"))
  adjacency_pairs_unique22$Cluster2 = adjacency_pairs_unique22$rpca_cluster
  adjacency_pairs_unique22$rpca_cluster = NULL
  
 
  # Return adjacency list with cluster annotations
  return(adjacency_pairs_unique22)
})

# Assign sample names to adjacency network list
names(adjacency_by_sample) = names(per_spot_neighbours)


##############################
# Generate All Cluster Pair Combinations by Sample
##############################
cluster_combn_by_samples = lapply(names(adjacency_by_sample), function(x) {
  
  # Extract the set of unique clusters observed in adjacency pairs
  c1 = unique(adjacency_by_sample[[x]]$Cluster1)
  c2 = unique(adjacency_by_sample[[x]]$Cluster2)
  
  # Combine both sets, convert to character, and sort consistently
  myclusters = mixedsort(unique(as.character(c1), as.character(c2)))
  
  # Generate all inter-cluster combinations (ClusterA_ClusterB, ClusterB_ClusterC, etc.)
  cluster_combn = apply(combn(myclusters, 2), 2, paste, collapse = '_')
  
  #Add within-cluster combinations (e.g., "Cluster1_Cluster1")
  mycluster_combn2 = union(
    unlist(lapply(myclusters, function(x) { paste0(x, '_', x) })),
    cluster_combn
  )
  
  # Return the full set of cluster combinations, sorted
  return(mixedsort(mycluster_combn2))
})

#Assign sample names to the cluster combination lists
names(cluster_combn_by_samples) = names(adjacency_by_sample)



##############################
# Count Adjacency Edges per Cluster Pair (by Sample)
##############################
adjacency_by_sample_count_by_cluster_pairs = lapply(names(adjacency_by_sample), function(x) {
  
  #Retrieve the full set of possible cluster pair combinations for this sample
  mycluster_combn2 = cluster_combn_by_samples[[x]]
  
  #Retrieve the adjacency dataframe for this sample
  myadj_df = adjacency_by_sample[[x]]
  
  #For each cluster pair, subset adjacency edges that match that pair
  unique_cluster_pairdf = lapply(mycluster_combn2, function(x) {
    cluster1 = str_split_i(x, '_', 1)
    cluster2 = str_split_i(x, '_', 2)
    
    # Debug prints (can be removed later)
    print(cluster1)
    print(cluster2)
    
    # Subset adjacency dataframe for this specific cluster pair
    mydf_sub = subset(myadj_df, Cluster1 == cluster1 & Cluster2 == cluster2)
    print(mydf_sub)  # Debug preview
    
    return(mydf_sub)
  })
  
  #Assign names (cluster pairs) to the list of subsets
  names(unique_cluster_pairdf) = mycluster_combn2
  
  #Combine all per-pair subsets into one dataframe
  unique_cluster_pairdf2 = do.call(rbind, unique_cluster_pairdf)
  rownames(unique_cluster_pairdf2) = NULL
  
  # Count number of adjacency edges per cluster pair
  #  since adjacency_by_sample is already neighbour-based)
  min_dist_by_sample_count = unique_cluster_pairdf2 %>%
    group_by(Cluster1, Cluster2) %>%
    summarize(count = n())
  
  #Return per-sample summary
  return(min_dist_by_sample_count)
})

#Preserve sample names in the output list
names(adjacency_by_sample_count_by_cluster_pairs) = names(adjacency_by_sample)


##############################
# Summarize Total Distinct Adjacency Pairs by Sample
##############################
# For each sample, count the number of distinct undirected adjacency pairs
total_adjacent_distinct_pairs_by_sample = lapply(
  names(adjacency_by_sample),
  function(x) {
    return(n_distinct(adjacency_by_sample[[x]]$Pair))
  }
)

# Assign sample names to the counts
names(total_adjacent_distinct_pairs_by_sample) = names(adjacency_by_sample)


##############################
# Count Adjacency Edges per Cluster Pair (with Sample Labels)
##############################
adjacency_by_sample_count = lapply(names(adjacency_by_sample), function(x) {
  
  # Step 1: Count number of adjacency edges for each Cluster1–Cluster2 pair
  mycnt = adjacency_by_sample[[x]] %>%
    group_by(Cluster1, Cluster2) %>%
    summarize(count = n())
  
  # Step 2: Add sample name column
  mycnt$Sample = x
  
  # Return per-sample adjacency counts
  return(mycnt)
})

# Preserve sample names as list element names
names(adjacency_by_sample_count) = names(adjacency_by_sample)


##############################
# Merge All Per-Sample Counts into One Table
##############################
# Bind all per-sample dataframes into a single dataframe for comparison
adjacency_by_sample_count_merged = do.call(rbind, adjacency_by_sample_count)

##############################
# Compute Adjacency Cluster Probability per Sample
##############################
adjacency_cluster_probability = lapply(names(adjacency_by_sample_count_by_cluster_pairs), function(x) {
  
  # Step 1: Normalize adjacency counts by total # of distinct adjacency pairs
  adjacency_by_sample_count_by_cluster_pairs[[x]]$adjacency_cluster_prob =
    adjacency_by_sample_count_by_cluster_pairs[[x]]$count /
    total_adjacent_distinct_pairs_by_sample[[x]]
  
  # Return annotated dataframe
  return(adjacency_by_sample_count_by_cluster_pairs[[x]])
})

# Step 2: Preserve sample names in the output list
names(adjacency_cluster_probability) = names(adjacency_by_sample_count_by_cluster_pairs)


##############################
# Compute Joint Cluster Probability Score (JCPS)
##############################
adjacency_cluster_probability_jcps = lapply(names(adjacency_cluster_probability), function(x) {
  
  # Join with per-cluster probabilities for Cluster1
  cluster_1_prob = adjacency_cluster_probability[[x]] %>%
    left_join(per_cluster_probability[[x]], by = c("Cluster1" = "integrated_snn_res.0.9"))
  
  # Store cluster1 probability explicitly
  cluster_1_prob$cluster1prob = cluster_1_prob$per_cluster_probability
  cluster_1_prob$per_cluster_probability = NULL
  
  # Join with per-cluster probabilities for Cluster2
  cluster_2_prob = cluster_1_prob %>%
    left_join(per_cluster_probability[[x]], by = c("Cluster2" = "integrated_snn_res.0.9"))
  
  # Store cluster2 probability explicitly
  cluster_2_prob$cluster2prob = cluster_2_prob$per_cluster_probability
  cluster_2_prob$per_cluster_probability = NULL
  
  # Compute JCPS (Joint Cluster Probability Score)
  # JCPS = observed adjacency probability / (expected probability from cluster sizes)
  cluster_2_prob$JCPS =
    cluster_2_prob$adjacency_cluster_prob /
    (cluster_2_prob$cluster1prob * cluster_2_prob$cluster2prob)
  
  # Return dataframe with JCPS
  return(cluster_2_prob)
})

#Assign sample names
names(adjacency_cluster_probability_jcps) = names(adjacency_cluster_probability)

##############################
# Define Global Cluster Combinations
##############################
#Extract all unique clusters from metadata and sort them
myclusters_all = mixedsort(as.character(unique(cluster_meta$integrated_snn_res.0.9)))

# Generate all inter-cluster pairwise combinations (ClusterA_ClusterB)

cluster_combn = apply(combn(myclusters_all,2),2,paste,collapse='_')

##############################
# Add Cluster Pair IDs and Sample Labels to JCPS Data
##############################

  #Create explicit cluster pair combination label for each row

adjacency_cluster_probability_jcps_mod = lapply(names(adjacency_cluster_probability_jcps),function(x) {
      adjacency_cluster_probability_jcps[[x]]$cluster_combn = paste0(adjacency_cluster_probability_jcps[[x]]$Cluster1,'_',adjacency_cluster_probability_jcps[[x]]$Cluster2)
      adjacency_cluster_probability_jcps[[x]]$Sample = x   #Add sample identifier
      return(adjacency_cluster_probability_jcps[[x]])   # Return modified dataframe
  })

names(adjacency_cluster_probability_jcps_mod) = names(adjacency_cluster_probability_jcps)


##############################
# Combine JCPS Data Across Samples
##############################
# Merge all per-sample JCPS results into one dataframe
adjacency_cluster_probability_jcps_mod_comb = do.call(rbind, adjacency_cluster_probability_jcps_mod)

# Split the combined dataframe by cluster pair combination
adjacency_cluster_probability_jcps_mod_comb_spl =
  split(adjacency_cluster_probability_jcps_mod_comb,
        adjacency_cluster_probability_jcps_mod_comb$cluster_combn)


##############################
# Compute Weighted JCPS per Cluster Combination
##############################
adjacency_cluster_probability_jcps_mod = lapply(
  names(adjacency_cluster_probability_jcps_mod_comb_spl),
  function(x) {
    mydf = adjacency_cluster_probability_jcps_mod_comb_spl[[x]]
    
    # Count in how many samples this cluster combination is observed
    sample_proportion =
      adjacency_cluster_probability_jcps_mod_comb_spl[[x]] %>%
      group_by(cluster_combn) %>%
      summarize(comb_by_sample = n())
    
    #Relative sample proportion
    # Ratio of observed samples to unobserved samples (out of 22 total)
    relative_sample_prop =
      sample_proportion$comb_by_sample / (22 - sample_proportion$comb_by_sample)
    
    # Proportion of samples where this cluster combination is present
    my_sample_prop = sample_proportion$comb_by_sample / 22
    
    # Weighted JCPS
    # Weighted average JCPS across samples, weighted by adjacency count
    weighted_jcps =
      (sum(mydf$JCPS * mydf$count) / sum(mydf$count))
    
    # Relative weighted JCPS
    relative_wt_jcps = weighted_jcps * relative_sample_prop
    
    #Create compact dataframe for this cluster combination
    mydf_mod = data.frame(
      Cluster_combs = x,
      Relative_JCP_score = relative_wt_jcps,
      sample_present_percent = my_sample_prop * 100
    )
    
    return(mydf_mod)
  }
)

#Assign names to the list
names(adjacency_cluster_probability_jcps_mod) =
  names(adjacency_cluster_probability_jcps_mod_comb_spl)


##############################
# Collapse Weighted JCPS Results into Final Edge List
##############################
weighted_jcps_df = do.call(rbind, adjacency_cluster_probability_jcps_mod)

#Split cluster combination string into Cluster1 and Cluster2
weighted_jcps_df$Cluster1 = str_split_i(weighted_jcps_df$Cluster_combs, '_', 1)
weighted_jcps_df$Cluster2 = str_split_i(weighted_jcps_df$Cluster_combs, '_', 2)

# #Rename columns to clean, final format
# colnames(weighted_jcps_df) =
#   c('Cluster_combs', 'Relative_JCP_score', 'sample_present_percent',
#     'Cluster1', 'Cluster2')

# #Reorder columns for edge list format
# Edge_df_reord = weighted_jcps_df[c('Cluster1', 'Cluster2',
#                                    'Relative_JCP_score', 'sample_present_percent')]

# colnames(Edge_df_reord) =
#   c('source', 'target', 'Weighted_JCP_score', 'sample_present_percent')


##############################
# Export Edge List as CSV for cytoscape input
##############################
write.csv(Edge_df_reord,
          file = "Edge_data_relative_JCPS.csv",
          quote = FALSE, row.names = FALSE)

# https://github.com/Petti-Lab/GBM_Spatial_paper_figures_code/blob/main/Figure_3/spots_to_cells_celltrek_annotations.rds
##############################
# Prepare Node information for Cytoscape###
spots_to_cells_celltrek_annotations = readRDS('/n/scratch/users/s/sak4832/spots_to_cells_celltrek_annotations.rds')

##############################
# Function: melt_by_samples_sub_mod
# Purpose:
#   - Collapse fine-grained cell type annotations per spot into broader categories
#   - Reshape spot-level data into long format
#   - Link to Seurat cluster metadata
#   - Summarize counts of cell types by cluster
##############################

melt_by_samples_sub_mod <- function(sample_name,
                                    spots_to_cells_celltrek_annotations,
                                    type_name,
                                    seurat_obj_list2,
                                    integ_clus) {
  
  # Extract spot-to-cell annotations for the given sample and type
  all_celltypes_per_spot2 =
    spots_to_cells_celltrek_annotations[[sample_name]][[type_name]]
  
  # Define broad cell type groupings
  Glia = c('Astrocyte','Committed oligodendrocyte precursor',
           'Oligodendrocyte','Oligodendrocyte precursor')
  Lymphoid = c('Bcells','CD4_Tcells','CD8_Tcells',
               'NK/NK-like','Other_Tcells')
  Myeloid = c('DC','Microglia','Monocyte','TAMs')
  Vascular = c('Vascular','Fibroblast')
  
  # Collapse fine-grained types into broad categories
  all_celltypes_per_spot2$Glia =
    rowSums(all_celltypes_per_spot2[intersect(Glia, colnames(all_celltypes_per_spot2))])
  
  all_celltypes_per_spot2$Lymphoid =
    rowSums(all_celltypes_per_spot2[intersect(Lymphoid, colnames(all_celltypes_per_spot2))])
  
  all_celltypes_per_spot2$Myeloid =
    rowSums(all_celltypes_per_spot2[intersect(Myeloid, colnames(all_celltypes_per_spot2))])
  
  all_celltypes_per_spot2$Vascular =
    rowSums(all_celltypes_per_spot2[intersect(Vascular, colnames(all_celltypes_per_spot2))])
  
  # Select key columns (broad categories + neurons + M1–M5 + total count)
  sel_cols = c('spot_id','Glia','Lymphoid','Myeloid','Neuron','Vascular',
               paste0('M',1:5),'total_cell_count')
  
  all_celltypes_per_spot2_sub = all_celltypes_per_spot2[sel_cols]
  
  # Reshape to long format (spot_id, variable, value)
  all_celltypes_per_spot2_melt = melt(all_celltypes_per_spot2_sub, id='spot_id')
  
  # Get Seurat metadata for the same sample (selected cluster resolution)
  seurat_obj_samp = seurat_obj_list2[[sample_name]]
  seurat_obj_meta = seurat_obj_samp@meta.data[c(integ_clus)]
  seurat_obj_meta$spot_id = rownames(seurat_obj_meta)
  
  # Filter melted data to include only spots present in Seurat metadata
  all_celltypes_per_spot2_melt_sub =
    subset(all_celltypes_per_spot2_melt, spot_id %in% seurat_obj_meta$spot_id)
  
  # Join cell type data with cluster metadata
  melt_df = all_celltypes_per_spot2_melt_sub %>%
    left_join(seurat_obj_meta, by = c("spot_id"))
  
  # Summarize counts of each cell type variable per cluster
  celltype_count_rpca_clusters = melt_df %>%
    group_by(variable, !!sym(integ_clus)) %>%
    summarize(count = sum(value))
  
  # Remove total_cell_count from results
  celltype_count_rpca_clusters_sub =
    subset(celltype_count_rpca_clusters, variable != 'total_cell_count')
  
  # Clean factor levels and add sample name
  celltype_count_rpca_clusters_sub$variable =
    droplevels(celltype_count_rpca_clusters_sub$variable)
  
  celltype_count_rpca_clusters_sub$sample_name = sample_name
  
  # Return tidy dataframe
  return(celltype_count_rpca_clusters_sub)
}

##############################
# Function: get_celltype_by_integclus_mod
# Purpose:
#   - Run melt_by_samples_sub_mod across multiple cluster resolutions (integ_clus_vec)
#   - Combine results across all samples per resolution
#   - Summarize cell type counts by cluster for each resolution
##############################

get_celltype_by_integclus_mod <- function(spots_to_cells_celltrek_annotations,
                                          type_name,
                                          integ_clus_vec,
                                          seurat_obj_list2 = seurat_obj_list2) {
  
  integ_clus_vec_celltypes = list()
  
  # Loop over all cluster resolutions and all samples
  for (integs in integ_clus_vec) {
    for (sample_name in names(spots_to_cells_celltrek_annotations)) {
      integ_clus_vec_celltypes[[integs]][[sample_name]] =
        melt_by_samples_sub_mod(
          sample_name = sample_name,
          spots_to_cells_celltrek_annotations = spots_to_cells_celltrek_annotations,
          type_name = type_name,
          seurat_obj_list2 = seurat_obj_list2,
          integ_clus = integs
        )
    }
  }
  
  # Combine per-sample results into a single dataframe for each cluster resolution
  integ_clus_vec_celltypes_combined =
    lapply(names(integ_clus_vec_celltypes), function(x) {
      mytemp = do.call(rbind, integ_clus_vec_celltypes[[x]])
      rownames(mytemp) = NULL
      return(mytemp)
    })
  
  names(integ_clus_vec_celltypes_combined) = names(integ_clus_vec_celltypes)
  
  # Summarize total cell type counts by cluster for each resolution
  integ_clus_vec_celltypes_combined_sumbyclus =
    lapply(names(integ_clus_vec_celltypes_combined), function(x) {
      celltype_count_rpca_clusters =
        integ_clus_vec_celltypes_combined[[x]] %>%
        group_by(variable, !!sym(x)) %>%
        summarize(celltype_count_by_cluster = sum(count))
    })
  
  names(integ_clus_vec_celltypes_combined_sumbyclus) =
    names(integ_clus_vec_celltypes_combined)
  
  return(integ_clus_vec_celltypes_combined_sumbyclus)
}

# Run cell type aggregation across cluster resolutions using CellTrek annotations and Seurat objects
integ_clus_vec_celltypes_integrated_snn_res.0.9 = get_celltype_by_integclus_mod(spots_to_cells_celltrek_annotations=spots_to_cells_celltrek_annotations,type_name='broad_celltypes_per_spot',integ_clus_vec=integ_clus_vec,seurat_obj_list=seurat_obj_list2)

Celltype_composition = integ_clus_vec_celltypes_integrated_snn_res.0.9$integrated_snn_res.0.9

# Reshape cell type composition to wide format: clusters as rows, cell types as columns, values = counts
Celltype_composition_dcast = Celltype_composition %>% dcast(integrated_snn_res.0.9 ~ variable,value.var='celltype_count_by_cluster')

# Build node attribute table: extract unique (sample_name, cluster) pairs across all Seurat objects
Node_attr = lapply(names(seurat_obj_list2),function(x) {
    integ_clus_df = seurat_obj_list2[[x]]@meta.data['integrated_snn_res.0.9']
    integ_clus_df$barcode = rownames(integ_clus_df)
    integ_clus_df$sample_name = x

    integ_clus_df$barcode = NULL

    rownames(integ_clus_df) = NULL

    return(unique(integ_clus_df[c('sample_name','integrated_snn_res.0.9')]))
    }
)

names(Node_attr) =  names(seurat_obj_list2) 
Node_attr_allsamp <- do.call(rbind, Node_attr)

rownames(Node_attr_allsamp) = NULL


Node_attr_allsamp <- do.call(rbind, Node_attr)

rownames(Node_attr_allsamp) = NULL

node_size_by_sample =  Node_attr_allsamp %>% group_by(integrated_snn_res.0.9) %>% summarize(Sample_count=n())

node_size_by_sample$Cluster = paste0('C',node_size_by_sample$`integrated_snn_res.0.9`)

# node_size_by_sample2 = node_size_by_sample[c('Cluster','Sample_count')]

Node_table = merge(node_size_by_sample,Celltype_composition_dcast,by='integrated_snn_res.0.9')


Edge_df = weighted_jcps_df
Edge_df$Cluster1 = paste0('C',Edge_df$Cluster1)
Edge_df$Cluster2 = paste0('C',Edge_df$Cluster2)

# Reorder and rename node/edge tables, then export as CSVs for network analysis

Edge_df_reord = Edge_df[c('Cluster1','Cluster2','Weighted_JCP_score')]
colnames(Edge_df_reord) = c('source','target','Weighted_JCP_score')

Node_table_reord = Node_table[c("Cluster","Sample_count","Glia","Lymphoid","Myeloid","Neuron","Vascular","M1","M2","M3","M4","M5")]

colnames(Node_table_reord) = c(c("id","Sample_count","Glia","Lymphoid","Myeloid","Neuron","Vascular","M1","M2","M3","M4","M5"))

write.csv(Node_table_reord,file = "Node_table_majorcelltypes.csv", quote = F,row.names=F)

write.csv(Edge_df_reord,file = "Edge_data_JCPS.csv",quote = F,row.names=F)

library(RCy3)

cytoscapePing ()
cytoscapeVersionInfo ()

# Load node/edge tables into Cytoscape, build adjacency network, and reshape edge list into a weighted matrix

nodes = read.csv('Node_table_majorcelltypes.csv')
edges = read.csv('Edge_data_JCPS.csv')

### Visualize in Cytoscape ########
createNetworkFromDataFrames(nodes,edges, title="Cluster Adjacency", collection="Joint conditional probability score")

Edge_df_reord_dcast = Edge_df_reord %>% dcast(source ~ target,value.var='Weighted_JCP_score')

Edge_df_reord_dcast[is.na(Edge_df_reord_dcast)] = 0

rownames(Edge_df_reord_dcast) = Edge_df_reord_dcast$source

Edge_df_reord_dcast$source = NULL
