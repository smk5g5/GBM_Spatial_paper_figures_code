# ============================================================================
# Figure 4A code for cluster adjacency analysis
# Author: [Saad Khan]
# Description: Enrichment dotplot of significant contact-dependent (within-spot) ligand receptor interactions
 # in visium data. Size of dot indicates enrichment score, bold circle around spot indicates significance based
 #  on adjusted p.values. Row Barplots show celltype diversity for the said LR whereas column barplots indicate 
 #  celltype diversity in the cluster. The clusters are ordered from edge rich clusters to core-rich clusters
 #  # ============================================================================
##############################
# Load Required Libraries
##############################
# Core single-cell analysis packages
library(SingleCellExperiment)
library(scuttle)
library(scran)
library(scater)
library(Seurat)
library(sctransform)
library(scCustomize)

# Parallelization
library(BiocParallel)
library(future)

# Visualization
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(circlize)
library(ComplexHeatmap)

# Data manipulation
library(dplyr)
library(stringr)
library(reshape2)

# Clustering and analysis
library(Ckmeans.1d.dp)
library(clusterProfiler)
library(org.Hs.eg.db)

# Configuration and YAML support
library(yaml)
# Set seed for reproducibility
set.seed(12345L)

# Create a sample mapping dataframe with sample names

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

mysample_map$image_name[mysample_map$image_name %in% 'WU1221.Core'] = 'WU1221.Core2'
mysample_map$image_name[mysample_map$image_name %in% 'WU1221.Edge'] = 'WU1221.Edge1'

# Load annotated CellChat interaction database (used for stlearn)

CellChatDB_df_noncomp2_withannot = read.csv('/n/scratch/users/s/sak4832/CellChatDB_df_noncomp2_withannot.csv')

# Subset only cell-cell contact based interactions

contact_based = subset(CellChatDB_df_noncomp2_withannot,annotation %in% c('Cell-Cell Contact'))

# Find all CSV files with ligand-receptor info from CellChat within-spot mode
Lig_rec_csvs = Sys.glob("/n/scratch/users/s/sak4832/myhoangdb_cci_withinspotmode_cellchatv2/*.csv")

# Parse ligand, receptor, and sample name information from filenames
myligrec_df = data.frame(filename = Lig_rec_csvs,
	sample_name=str_split_i(str_split_i(basename(Lig_rec_csvs), "_", 3),"\\.",1),
	Ligand=str_split_i(basename(Lig_rec_csvs), "_", 1),
	Receptor=str_split_i(basename(Lig_rec_csvs), "_", 2))

# Remove entries where ligand and receptor are the same
myligrec_df_sub = subset(myligrec_df,myligrec_df$Ligand!=myligrec_df$Receptor)
myligrec_df_list = split(myligrec_df_sub,myligrec_df_sub$sample_name)

# Function to read ligand-receptor CSVs and return a combined dataframe with significance and spot info

read_lig_rec_data <- function(myligrec_df_list,sample_name){
myligrec_csv_df_list = list()

	mysamp = mysample_map$image_name[mysample_map$sample_name==sample_name]

for(i in 1:nrow(myligrec_df_list[[sample_name]])){

	mylig = myligrec_df_list[[sample_name]]$Ligand[i]
	myrec = myligrec_df_list[[sample_name]]$Receptor[i]
	mydf = read.csv(myligrec_df_list[[sample_name]]$filename[i],header=T)
	mydf$significant = 0
	mydf$significant[mydf$p_adjs<=0.05] = 1
	mydf_sub = mydf[c('Barcodes','significant')]
	mydf_sub$spot_id = paste0(mysamp,'_',mydf_sub$Barcodes)
	mydf_sub$Lig_Rec = 	paste0(mylig,"_",myrec)
	myligrec_csv_df_list[[paste0(mylig,"_",myrec)]] = mydf_sub
}

myligrec_csv_df_combined <- do.call(rbind, myligrec_csv_df_list)
rownames(myligrec_csv_df_combined) = NULL
return(myligrec_csv_df_combined)
}

# Initialize list to store binary matrices for significant ligand-receptor interactions (per sample)

ligrec_prep_sub_list = list()
# For each sample, generate a wide-format binary matrix of significant ligand-receptor pairs

for(i in names(myligrec_df_list)){
	ligrec_df = read_lig_rec_data(myligrec_df_list,sample_name=i)
	ligrec_prep = ligrec_df %>% dcast(spot_id ~ Lig_Rec,value.var='significant')
	rownames(ligrec_prep) = ligrec_prep$spot_id
	ligrec_prep$spot_id = NULL
	ligrec_prep_sub_list[[i]] = ligrec_prep_sub = ligrec_prep[(rowSums(ligrec_prep)>1),]
}

# Load integrated Seurat object containing metadata (e.g., clusters, patient IDs)

merged.patient.combined.integ_sub = readRDS('/n/scratch/users/s/sak4832/August_27th_2025/Jul17th_2025/June2_2025/April23_2025/March12th_2025/Feb2_2025/Dec27_2024/cluster_distances/merged.patient.combined.integ_sub.rds')
	# /n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/cluster_distances/merged.patient.combined.integ_sub.rds')

# Extract relevant metadata from integrated Seurat object

myrpca_df = merged.patient.combined.integ_sub@meta.data[c('Sample','Patient','integrated_snn_res.0.9')]

# Split metadata by sample

myrpca_df_list = split(myrpca_df,myrpca_df$Sample)

# Add barcode information to each sample-specific metadata dataframe

myrpca_df_list2 =  lapply(names(myrpca_df_list), function(x) {
	myrpca_df_list[[x]]$Barcode = rownames(myrpca_df_list[[x]])
	# str_split_i(rownames(myrpca_df_list[[x]]),'_',2)
	return(myrpca_df_list[[x]])
	})

names(myrpca_df_list2) = names(myrpca_df_list)

# Rename elements in ligand-receptor list to use standardized image names
names(ligrec_prep_sub_list) = mysample_map$image_name

# Rename elements in metadata list to match standardized image names
names(myrpca_df_list2) = mysample_map$image_name

# For each sample, merge the ligand-receptor matrix with sample metadata by barcode
ligrec_prep_sub_list_rpca_comb = lapply(names(ligrec_prep_sub_list), function(x) {
  return(merge(ligrec_prep_sub_list[[x]], myrpca_df_list2[[x]], by.x = 0, by.y = 'Barcode'))
})

# Assign sample names back to the merged list
names(ligrec_prep_sub_list_rpca_comb) = names(ligrec_prep_sub_list)

# For each sample, create a composite ID from Sample, Patient, and cluster, then remove redundant columns

ligrec_prep_sub_list_rpca_comb =  lapply(names(ligrec_prep_sub_list_rpca_comb), function(x) {
ligrec_prep_sub_list_rpca_comb[[x]]$id = paste0(ligrec_prep_sub_list_rpca_comb[[x]]$Sample,'_',ligrec_prep_sub_list_rpca_comb[[x]]$Patient,'_',ligrec_prep_sub_list_rpca_comb[[x]]$`integrated_snn_res.0.9`)
ligrec_prep_sub_list_rpca_comb[[x]]$Sample=NULL
ligrec_prep_sub_list_rpca_comb[[x]]$Patient=NULL
ligrec_prep_sub_list_rpca_comb[[x]]$integrated_snn_res.0.9=NULL
return(ligrec_prep_sub_list_rpca_comb[[x]])
})

names(ligrec_prep_sub_list_rpca_comb) = names(ligrec_prep_sub_list)

# ligrec_prep_sub_list_rpca_comb_summary =  lapply(names(ligrec_prep_sub_list_rpca_comb), function(x) {
# 	return(ligrec_prep_sub_list_rpca_comb[[x]] %>% group_by(id)  %>% summarize(across(where(is.numeric), sum)))
# })

# names(ligrec_prep_sub_list_rpca_comb_summary) = names(ligrec_prep_sub_list_rpca_comb)

# For each sample, aggregate ligand-receptor values by 'id' (which encodes Sample_Patient_Cluster)
ligrec_prep_sub_list_rpca_comb_summary = lapply(
  names(ligrec_prep_sub_list_rpca_comb), 
  function(x) {
    return(
      ligrec_prep_sub_list_rpca_comb[[x]] %>%
        group_by(id) %>%
        summarize(across(where(is.numeric), sum))
    )
  }
)
names(ligrec_prep_sub_list_rpca_comb_summary) = names(ligrec_prep_sub_list_rpca_comb)

# For each sample, extract the cluster ID from 'id', convert to factor, and aggregate values by cluster
ligrec_prep_sub_list_rpca_comb_clusinfo = lapply(
  names(ligrec_prep_sub_list_rpca_comb),
  function(x) {
    # Extract cluster from 'id' (3rd field in Sample_Patient_Cluster)
    ligrec_prep_sub_list_rpca_comb[[x]]$`integrated_snn_res.0.9` =
      as.factor(str_split_i(ligrec_prep_sub_list_rpca_comb[[x]]$id, '_', 3))
    
    # Drop the composite id column after extracting cluster info
    ligrec_prep_sub_list_rpca_comb[[x]]$id = NULL

    # Group by cluster and summarize ligand-receptor interactions per cluster
    return(
      ligrec_prep_sub_list_rpca_comb[[x]] %>%
        group_by(integrated_snn_res.0.9) %>%
        summarize(across(where(is.numeric), sum))
    )
  }
)
names(ligrec_prep_sub_list_rpca_comb_clusinfo) = names(ligrec_prep_sub_list_rpca_comb)

# Create a binary matrix indicating the presence of each ligand-receptor pair across samples

LR_list_mat = list_to_matrix(lapply(lapply(ligrec_prep_sub_list_rpca_comb_clusinfo,colnames),setdiff,'integrated_snn_res.0.9'))

# Function to extract and summarize one LR pair across clusters in all samples where it's present
read_lig_rec_data <- function(ligrec_prep_sub_list_rpca_comb_clusinfo,LR_list_mat,LR_name){
  # Identify samples where the given LR pair is present

sel_samples = names(which((LR_list_mat[LR_name,]>0)==T))
  # Subset the full cluster-wise list to only those samples

ligrec_prep_sub_list_rpca_comb_clusinfo_sub = ligrec_prep_sub_list_rpca_comb_clusinfo[sel_samples]
  # Extract only the cluster label and the selected LR column from each sample

ligrec_prep_sub_list_rpca_comb_clusinfo_sub2 = lapply(ligrec_prep_sub_list_rpca_comb_clusinfo_sub, function(x) {return(x[c('integrated_snn_res.0.9',LR_name)])})
  # Combine all samples into one dataframe

ligrec_prep_sub_list_rpca_comb_clusinfo_sub2_combined <- do.call(rbind, ligrec_prep_sub_list_rpca_comb_clusinfo_sub2)
  # Aggregate by cluster and sum the LR activity

return(ligrec_prep_sub_list_rpca_comb_clusinfo_sub2_combined %>% group_by(integrated_snn_res.0.9)  %>% summarize(across(where(is.numeric), sum)))
}

# Apply read_lig_rec_data() to each LR pair to compute cluster-level summaries across samples
summarized_LR_int_by_rpcacluster = lapply(
  rownames(LR_list_mat),
  function(x) {
    read_lig_rec_data(
      ligrec_prep_sub_list_rpca_comb_clusinfo = ligrec_prep_sub_list_rpca_comb_clusinfo,
      LR_list_mat = LR_list_mat,
      LR_name = x
    )
  }
)

summarized_LR_int_by_rpcacluster_melt = lapply(summarized_LR_int_by_rpcacluster,function(x) {melt(x)})


summarized_LR_int_by_rpcacluster_melt_comb = do.call(rbind, summarized_LR_int_by_rpcacluster_melt)

summarized_LR_int_by_rpcacluster_melt_comb$Ligand = str_split_i(summarized_LR_int_by_rpcacluster_melt_comb$variable,'_',1)
summarized_LR_int_by_rpcacluster_melt_comb$Receptor = str_split_i(summarized_LR_int_by_rpcacluster_melt_comb$variable,'_',2)


summarized_LR_int_by_rpcacluster_melt_comb_sub = subset(summarized_LR_int_by_rpcacluster_melt_comb,Ligand %in% unique(contact_based$ligand) & Receptor %in% unique(contact_based$receptor))


# Remove collagens integrins Laminins and Osteopontin from this analysis
summarized_LR_int_by_rpcacluster_melt_comb_sub_noitgs_collagens_integs = summarized_LR_int_by_rpcacluster_melt_comb_sub[grep('SPP1|ITG|LAM|COL',summarized_LR_int_by_rpcacluster_melt_comb_sub$variable,invert=T),]

# Pivot long-format LR interaction data into wide format: clusters as rows, LR pairs as columns
summarized_LR_int_by_rpcacluster_melt_comb_dcast = 
  summarized_LR_int_by_rpcacluster_melt_comb_sub_noitgs_collagens_integs %>%
  dcast(integrated_snn_res.0.9 ~ variable, value.var = 'value')


# Replace NA values (missing LR interactions) with 0 to indicate absence
summarized_LR_int_by_rpcacluster_melt_comb_dcast[is.na(summarized_LR_int_by_rpcacluster_melt_comb_dcast)] = 0

# Count number of spots per integrated cluster
rpca_cluster_spot_count = myrpca_df %>%
  group_by(integrated_snn_res.0.9) %>%
  summarise(n())

# Rename columns for clarity: cluster ID and spot count
colnames(rpca_cluster_spot_count) = c('integrated_snn_res.0.9', 'spot_count')


# Perform hypergeometric enrichment test for a specific LR pair in a given cluster
hypergeometric_test_LR <- function(
  summarized_LR_int_by_rpcacluster_melt_comb_dcast,
  rpca_cluster_spot_count,
  sel_LR,
  Cluster
) {
  # Get the number of significant LR interactions (k) in the specified cluster
  mydf = subset(summarized_LR_int_by_rpcacluster_melt_comb_dcast[c('integrated_snn_res.0.9', sel_LR)],
                integrated_snn_res.0.9 == Cluster)
  k = mydf[[sel_LR]]

  # Total number of spots across all clusters
  N = sum(rpca_cluster_spot_count$spot_count)

  # Total number of spots with the selected LR interaction across all clusters
  K = sum(summarized_LR_int_by_rpcacluster_melt_comb_dcast[, c(sel_LR)])

  # Number of spots in the specified cluster
  n = subset(rpca_cluster_spot_count, integrated_snn_res.0.9 == Cluster)$spot_count

  # Enrichment score: observed vs expected proportion
  enrichment_score = (k / K) / (n / N)

  # Hypergeometric p-value for overrepresentation
  pval = phyper(k, K, N - K, n, lower.tail = FALSE)

  # Assemble results in a dataframe
  final_df = data.frame(
    LR = sel_LR,
    enrichment_score = enrichment_score,
    pval = pval,
    Cluster = Cluster,
    K = K,
    k = k
  )

  return(final_df)
}

# Initialize an empty list to store hypergeometric test results for each LR × cluster combination
final_df_list = list()

# Loop over all ligand-receptor pairs (columns, excluding cluster label column)
for (mylr in colnames(summarized_LR_int_by_rpcacluster_melt_comb_dcast[-1])) {
  
  # Loop over all unique clusters
  for (myclus in unique(rpca_cluster_spot_count$`integrated_snn_res.0.9`)) {
    
    # Perform hypergeometric test and store the result in a list
    final_df_list[[paste0(mylr, '_', myclus)]] = hypergeometric_test_LR(
      summarized_LR_int_by_rpcacluster_melt_comb_dcast,
      rpca_cluster_spot_count,
      sel_LR = mylr,
      Cluster = myclus
    )
  }
}

# Apply Bonferroni correction to p-values across all LR-cluster enrichment tests
final_df_list_comb$padj = p.adjust(final_df_list_comb$pval, method = 'bonferroni', n = nrow(final_df_list_comb))

# Select unique rows of interest: one entry per LR-cluster pair with relevant stats
final_df_list_comb_selbyclus = unique(
  final_df_list_comb[c('LR', 'Cluster', 'enrichment_score', 'K', 'k', 'padj')]
)

# Remove rows with NA values to ensure clean downstream analysis
final_df_list_comb_selbyclus_nonan = final_df_list_comb_selbyclus[
  complete.cases(final_df_list_comb_selbyclus), 
]


# Melt each sample's cluster-level LR summary into long format and annotate with sample metadata
ligrec_prep_sub_list_rpca_comb_clusinfo_sampleadd = lapply(
  names(ligrec_prep_sub_list_rpca_comb_clusinfo),
  function(x) {
    # Melt wide-format LR matrix into long format for one sample
    mymlt = melt(ligrec_prep_sub_list_rpca_comb_clusinfo[[x]])
    
    # Add Sample and Patient metadata
    mymlt$Sample = x
    mymlt$Patient = str_split_i(x, '\\.', 1)  # Extract patient name from sample
    
    # Extract Ligand and Receptor from variable name
    mymlt$Ligand = str_split_i(mymlt$variable, '_', 1)
    mymlt$Receptor = str_split_i(mymlt$variable, '_', 2)
    
    # Filter to non-zero values and contact-based ligand-receptor pairs
    mymlt_sub = subset(mymlt, value > 0) %>%
      subset(Ligand %in% unique(contact_based$ligand) &
             Receptor %in% unique(contact_based$receptor))
    
    return(mymlt_sub)
  }
)

# Combine all samples into a single long-format dataframe
ligrec_prep_sub_list_rpca_comb_clusinfo_sampleadd_comb = do.call(
  rbind,
  ligrec_prep_sub_list_rpca_comb_clusinfo_sampleadd
)


# Group LR interactions by ligand, receptor, and cluster, then summarize sample/patient diversity
ligrec_prep_sub_list_rpca_comb_clusinfo_sampleadd_comb_filtered = 
  ligrec_prep_sub_list_rpca_comb_clusinfo_sampleadd_comb %>%
  
  group_by(Ligand, Receptor, integrated_snn_res.0.9) %>%
  summarise(
    distinct_patients = n_distinct(Patient),                     # Number of unique patients per LR-cluster
    distinct_samples = n_distinct(Sample),                       # Number of unique samples per LR-cluster
    B186_count = n_distinct(Sample[Patient == "B186"])           # Number of B186 samples supporting the LR
  ) %>%
  ungroup() %>%
  
  # Keep interactions seen in ≥2 patients and >2 samples, and not dominated by B186
  filter(
    distinct_patients >= 2 & 
    distinct_samples > 2 & 
    B186_count < distinct_samples
  )


# Create combined ligand-receptor pair ID for easier reference and plotting
ligrec_prep_sub_list_rpca_comb_clusinfo_sampleadd_comb_filtered$LR = paste0(
  ligrec_prep_sub_list_rpca_comb_clusinfo_sampleadd_comb_filtered$Ligand, "_",
  ligrec_prep_sub_list_rpca_comb_clusinfo_sampleadd_comb_filtered$Receptor
)

final_df_list_comb_selbyclus_nonan_filtered = subset(final_df_list_comb_selbyclus_nonan,LR %in% unique(ligrec_prep_sub_list_rpca_comb_clusinfo_sampleadd_comb_filtered$LR))

# Join filtered enrichment results with sample and patient metadata for enriched LR–cluster pairs
final_df_list_comb_selbyclus_nonan_filtered_samp_pat = final_df_list_comb_selbyclus_nonan_filtered %>%
  left_join(
    ligrec_prep_sub_list_rpca_comb_clusinfo_sampleadd_comb_filtered,
    by = c("LR" = "LR", "Cluster" = "integrated_snn_res.0.9")
  )


 # Remove rows where join failed (i.e., LR-cluster combinations that did not meet robustness filters)
final_df_list_comb_selbyclus_nonan_filtered_samp_pat = 
  final_df_list_comb_selbyclus_nonan_filtered_samp_pat[
    complete.cases(final_df_list_comb_selbyclus_nonan_filtered_samp_pat),
  ]

# Initialize column to annotate statistical significance
final_df_list_comb_selbyclus_nonan_filtered_samp_pat$significant = NA

# Annotate entries with Bonferroni-adjusted p-value < 0.05 as significant
final_df_list_comb_selbyclus_nonan_filtered_samp_pat$significant[
  final_df_list_comb_selbyclus_nonan_filtered_samp_pat$padj < 0.05
] = 'Significant'

# Annotate others as not significant
final_df_list_comb_selbyclus_nonan_filtered_samp_pat$significant[
  final_df_list_comb_selbyclus_nonan_filtered_samp_pat$padj > 0.05
] = 'Not-Significant'

# Add -log10(padj) column for visualization (e.g., volcano/dot plots)
# Use +1 to avoid -Inf in case padj == 0
final_df_list_comb_selbyclus_nonan_filtered_samp_pat$`-log10padj` =
  -log10(final_df_list_comb_selbyclus_nonan_filtered_samp_pat$padj + 1)


# Replace padj == 0 with a very small value to avoid -Inf in log-scale plots
final_df_list_comb_selbyclus_nonan_filtered_samp_pat = 
  final_df_list_comb_selbyclus_nonan_filtered_samp_pat %>%
  mutate(padj = ifelse(padj == 0, 1e-300, padj))

# Identify LR pairs that are significant in at least one cluster
all_sigLRs = setdiff(
  unique(final_df_list_comb_selbyclus_nonan_filtered_samp_pat$LR),
  (
    final_df_list_comb_selbyclus_nonan_filtered_samp_pat %>%
      group_by(LR) %>%
      summarise(Significant_count = sum(significant == "Significant")) %>%
      filter(Significant_count == 0)
  )$LR
)

# Identify clusters with at least one significant LR pair
all_sigclusters = setdiff(
  unique(final_df_list_comb_selbyclus_nonan_filtered_samp_pat$Cluster),
  (
    final_df_list_comb_selbyclus_nonan_filtered_samp_pat %>%
      group_by(Cluster) %>%
      summarise(Significant_count = sum(significant == "Significant")) %>%
      filter(Significant_count == 0)
  )$Cluster
)

# Subset to only LR–cluster combinations where both LR and cluster were significant at least once
final_df_list_comb_selbyclus_nonan_filtered_samp_pat_allsigLRs_allsigClus = subset(
  final_df_list_comb_selbyclus_nonan_filtered_samp_pat,
  LR %in% all_sigLRs & Cluster %in% all_sigclusters
)

# Load per-spot annotations (e.g., CellTrek-derived spatial mappings)
spots_to_cells_celltrek_annotations = readRDS(
  '/n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/cluster_distances/spots_to_cells_celltrek_annotations.rds'
)

# Export the full filtered and annotated enrichment table to CSV
write.csv(
  final_df_list_comb_selbyclus_nonan_filtered_samp_pat,
  file = 'Contact_dependent_interactions_standard_thresholds.csv',
  row.names = FALSE,
  quote = FALSE
)

