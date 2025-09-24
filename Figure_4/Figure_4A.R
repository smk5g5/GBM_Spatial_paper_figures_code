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
  '/n/scratch/users/s/sak4832/spots_to_cells_celltrek_annotations.rds'
)

# Export the full filtered and annotated enrichment table to CSV
write.csv(
  final_df_list_comb_selbyclus_nonan_filtered_samp_pat,
  file = 'Contact_dependent_interactions_standard_thresholds.csv',
  row.names = FALSE,
  quote = FALSE
)

ligrec_prep_sub_list_rpca_comb_clusinfo_celltrek = lapply(names(ligrec_prep_sub_list), function(x) {

    # Get LR interaction matrix for current sample
    ligrec_prep = ligrec_prep_sub_list[[x]]
    ligrec_prep$spot_id <- rownames(ligrec_prep)

    # Melt matrix into long format: one row per spot × LR pair
    ligrec_long <- melt(ligrec_prep, id.vars = "spot_id", 
                        variable.name = "Lig_Rec", value.name = "significant")

    # Keep only rows with significant interactions (non-zero)
    ligrec_long <- ligrec_long[ligrec_long$significant != 0, ]

    # Extract Ligand and Receptor names from LR column
    ligrec_long$Ligand   = str_split_i(ligrec_long$Lig_Rec, '_', 1)
    ligrec_long$Receptor = str_split_i(ligrec_long$Lig_Rec, '_', 2)

    # Subset to only contact-based LR pairs (biologically meaningful interactions)
    ligrec_long_sub = subset(ligrec_long,
                             Ligand %in% unique(contact_based$ligand) &
                             Receptor %in% unique(contact_based$receptor))

    # Add barcode column to RPCA metadata
    myrpca_df_list2[[x]]$barcode2 = rownames(myrpca_df_list2[[x]])

    # Identify overlapping barcodes between RPCA and LR data
    select_bars = intersect(myrpca_df_list2[[x]]$Barcode, unique(ligrec_long_sub$spot_id))
    ligrec_long_sub_byrpcasub = subset(ligrec_long_sub, spot_id %in% select_bars)

    # Also ensure those spots are present in CellTrek annotations
    select_bars2 = intersect(
        unique(spots_to_cells_celltrek_annotations[[x]][['broad_celltypes_per_spot']]$spot_id),
        unique(ligrec_long_sub_byrpcasub$spot_id)
    )
    ligrec_long_sub_byrpcasub = subset(ligrec_long_sub_byrpcasub, spot_id %in% select_bars2)

    # Merge LR data with RPCA metadata and CellTrek cell type annotations
    ligrec_long_rpca_celltrek = ligrec_long_sub_byrpcasub %>%
        left_join(myrpca_df_list2[[x]], by = c("spot_id" = "Barcode")) %>%
        left_join(spots_to_cells_celltrek_annotations[[x]][['broad_celltypes_per_spot']], by = "spot_id")

    # Collapse CellTrek annotations into broader categories
    ligrec_long_rpca_celltrek$Glia = rowSums(ligrec_long_rpca_celltrek[
        intersect(c("Astrocyte", "Committed oligodendrocyte precursor", 
                    "Oligodendrocyte", "Oligodendrocyte precursor"),
                  colnames(ligrec_long_rpca_celltrek))
    ])

    ligrec_long_rpca_celltrek$Lymphoid = rowSums(ligrec_long_rpca_celltrek[
        intersect(c("Bcells", "CD4_Tcells", "CD8_Tcells", 
                    "NK/NK-like", "Other_Tcells"),
                  colnames(ligrec_long_rpca_celltrek))
    ])

    ligrec_long_rpca_celltrek$Myeloid = rowSums(ligrec_long_rpca_celltrek[
        intersect(c("DC", "Microglia", "Monocyte", "TAMs"),
                  colnames(ligrec_long_rpca_celltrek))
    ])

    ligrec_long_rpca_celltrek$Vascular = rowSums(ligrec_long_rpca_celltrek[
        intersect(c("Fibroblast", "Vascular"),
                  colnames(ligrec_long_rpca_celltrek))
    ])

    return(ligrec_long_rpca_celltrek)
})

# Assign sample names to the resulting list
names(ligrec_prep_sub_list_rpca_comb_clusinfo_celltrek) = names(ligrec_prep_sub_list)


# Function to compute total cell type counts grouped by any variable (e.g., Lig_Rec or Cluster)
get_celltype_freq_mat_by_group <- function(celltype_data_bysamples, group_var) {
  
  # For each sample, summarize cell type composition by the selected grouping variable
  celltype_freq_list = lapply(celltype_data_bysamples, function(x) {
    
    x %>%
      group_by(!!sym(group_var)) %>%  # Dynamically group by e.g., 'Lig_Rec'
      summarize(
        M1 = sum(M1),
        M2 = sum(M2),
        M3 = sum(M3),
        M4 = sum(M4),
        M5 = sum(M5),
        Neuron   = sum(Neuron),
        Glia     = sum(Glia),
        Lymphoid = sum(Lymphoid),
        Myeloid  = sum(Myeloid),
        Vascular = sum(Vascular)
      )
  })

  # Combine the summarized cell type tables across all samples
  celltype_freq_bygrp = do.call(rbind, celltype_freq_list)
  return(celltype_freq_bygrp)
}

# Apply the function to your LR-annotated spatial data, grouped by ligand-receptor pairs
celltype_freq_by_LR = get_celltype_freq_mat_by_group(
  celltype_data_bysamples = ligrec_prep_sub_list_rpca_comb_clusinfo_celltrek,
  group_var = 'Lig_Rec'
)


# Aggregate cell type frequencies per ligand-receptor pair across all samples
celltype_freq_by_LR2 = celltype_freq_by_LR %>%
  group_by(Lig_Rec) %>%
  summarize(
    M1        = sum(M1),
    M2        = sum(M2),
    M3        = sum(M3),
    M4        = sum(M4),
    M5        = sum(M5),
    Neuron    = sum(Neuron),
    Glia      = sum(Glia),
    Lymphoid  = sum(Lymphoid),
    Myeloid   = sum(Myeloid),
    Vascular  = sum(Vascular)
  )

# Filter enriched LR interactions with at least 10 supporting spots in the cluster (k) 
# and 100 across all clusters (K); significance filtering is optional here
final_df_list_comb_Klargerthan100_and_klargerthan10_sig = subset(
  final_df_list_comb_selbyclus, 
  k >= 10 & K >= 100
  # & padj <= 0.05  # ← can be uncommented to filter for statistical significance
)

# For each cluster, select top 5 enriched LR interactions based on enrichment score
final_df_list_comb_selbyclus_sig_Kgrt100_top10byclus = 
  final_df_list_comb_Klargerthan100_and_klargerthan10_sig %>%
  group_by(Cluster) %>%
  slice_max(order_by = enrichment_score, n = 5)

# Add a formatted cluster label (e.g., C1, C2, ...) for use in plotting or reporting
final_df_list_comb_selbyclus_sig_Kgrt100_top10byclus$cluster_id = 
  paste0('C', final_df_list_comb_selbyclus_sig_Kgrt100_top10byclus$Cluster)

# clustering_matrix <- final_df_list_comb_selbyclus_sig_Kgrt100_top10byclus %>%
#   ungroup() %>%  # Ungroup the data frame
#   dplyr::select(cluster_id, LR, `-log10padj`) %>%  # Explicit namespace for select
#   tidyr::pivot_wider(names_from = cluster_id, values_from = `-log10padj`, values_fill = 0) %>%
#   column_to_rownames(var = "LR") %>%
#   as.matrix()

# # Perform hierarchical clustering
# row_dendrogram <- hclust(dist(clustering_matrix))
# col_dendrogram <- hclust(dist(t(clustering_matrix)))

# ordered_rows <- rownames(clustering_matrix)[row_dendrogram$order]
# ordered_cols <- colnames(clustering_matrix)[col_dendrogram$order]

# Subset LR-cluster combinations that passed both LR- and cluster-level significance filters
final_df_list_comb_selbyclus_nonan_filtered_samp_pat_allsigLRs_allsigClus = subset(
  final_df_list_comb_selbyclus_nonan_filtered_samp_pat,
  LR %in% all_sigLRs & Cluster %in% all_sigclusters
)

# Remove ECM-related interactions: integrins (ITG), collagens (COL), laminins (LAM), SPP1
final_df_list_comb_selbyclus_nonan_filtered_samp_pat_allsigLRs_allsigClus_nocolnoint = 
  final_df_list_comb_selbyclus_nonan_filtered_samp_pat_allsigLRs_allsigClus[
    grep('SPP1|ITG|LAM|COL',
         final_df_list_comb_selbyclus_nonan_filtered_samp_pat_allsigLRs_allsigClus$LR,
         invert = TRUE),
  ]

# Identify LR pairs with at least one significant instance
select_LRs = final_df_list_comb_selbyclus_nonan_filtered_samp_pat %>%
  group_by(LR) %>%
  summarise(Significant_count = sum(significant == "Significant")) %>%
  filter(Significant_count != 0)

# Select top 5 LR pairs per cluster based on highest occurrence (`k`), breaking ties by lower `padj`
top_LRs_byclus = final_df_list_comb_selbyclus_nonan_filtered_samp_pat_allsigLRs_allsigClus_nocolnoint %>%
  group_by(Cluster) %>%
  arrange(desc(k), padj) %>%
  slice_max(order_by = k, n = 5)

# Subset only the LR-cluster combinations corresponding to top LR pairs for heatmap visualization
forhtmap = subset(
  final_df_list_comb_selbyclus_nonan_filtered_samp_pat_allsigLRs_allsigClus_nocolnoint,
  LR %in% unique(top_LRs_byclus$LR)
)

# Define named color palette for cell types for use in barplots or heatmaps
htmap_barplot_colors = c(
  "#E58606",  # M1
  "#6BAED6",  # M2
  "#E7BA52",  # M3
  "#5254A3",  # M4
  "#D66B6F",  # M5
  "red4",     # Neuron
  "gray50",   # glia
  "#4C8CE6",  # Lymphoid
  "#70CCB0",  # Myeloid
  "#FF4219"   # Vascular
)
names(htmap_barplot_colors) = c(
  "M1", "M2", "M3", "M4", "M5", 
  "Neuron", "glia", "Lymphoid", "Myeloid", "Vascular"
)

# Define custom color gradient function for dot plot using colorRamp2 (from circlize)
col_fun_dotplot <- colorRamp2(
  c(-0.3, -0.1, -0.01, 0),
  c("#023FA5", "#BEC1D4", "#C6909A", "#8E063B")
)

# Custom function to scale dot size based on a numeric input (e.g., LR occurrence); ignores zeros
dot_size_fun <- function(x) {
  ifelse(x > 0, sqrt(x) * 1.5, NA)  # Scale dots, hide for x <= 0
}

# Create a standalone legend object for the color scale (e.g., in dot plots or heatmaps)
lgd = Legend(
  col_fun = col_fun_dotplot,     # Color mapping function you defined earlier
  title = "-log10(padj+1)"       # Legend title (adjusted to avoid -Inf from padj = 0)
)

# Load cluster annotation or metadata (e.g., cluster labels, groupings, or annotations)
group_df2 = readRDS('groupdf2.rds')

# Subset the dataframe to include only clusters present in `col_order`
selgroup_df2 = group_df2[
  group_df2$integrated_snn_res.0.9 %in% col_order,
]

# Define locale for alphanumeric string sorting (e.g., C1, C10, C11... not C1, C11, C2...)
locale <- list(locale = "en_US", numeric = TRUE)

# Reorder 'group' factor levels to impose biologically meaningful ordering
selgroup_df2$group = factor(selgroup_df2$group, levels = c("edge-rich", "non-specific", "transition", "core-rich"))

# Sort by group and then numerically-aware cluster ID using stringi's locale-aware sorting
selgroup_df2_ord = selgroup_df2 %>%
  arrange(group, stri_rank(integrated_snn_res.0.9, opts_collator = locale))

# Extract final group ordering vector (parallel to cluster ordering)
group_order = as.character(selgroup_df2_ord$group)

# Extract color palette for group annotations (defined earlier in your dotplot settings)
mygroup_cols = dotplot_items$group_color_palette

# Extract final cluster order for columns in heatmap or facets in plot
col_order = as.character(selgroup_df2_ord$integrated_snn_res.0.9)

# Generates a ComplexHeatmap dot plot of ligand-receptor enrichment
# for contact-dependent LRs across clusters or samples, with dot size, 
# color, significance, and cell-type barplot annotations.
make_complexheatmap_signalling_withinspot_enrichment <- function(interaction_df,interactionid_var,out_prefix){
dotplot_matrix <- interaction_df %>%
  group_by(LR, !!sym(interactionid_var)) %>%
  summarise(enrichment_score = mean(enrichment_score, na.rm = TRUE), .groups = "drop") %>%  # mean duplicates
  pivot_wider(names_from = !!sym(interactionid_var), values_from = enrichment_score, values_fill = 0) %>%
  column_to_rownames("LR") %>%
  as.matrix()

dotplot_color_matrix <- interaction_df %>%
group_by(LR, !!sym(interactionid_var)) %>%
summarise(`-log10padj` = mean(`-log10padj`, na.rm = TRUE), .groups = "drop") %>%  # Use mean to avoid duplicates
distinct(LR, !!sym(interactionid_var), .keep_all = TRUE) %>%
pivot_wider(names_from = !!sym(interactionid_var), values_from = `-log10padj`, values_fill = 0) %>%
column_to_rownames("LR") %>%
as.matrix()


celltrek_spot_count_matrix <- interaction_df %>%
group_by(LR, !!sym(interactionid_var)) %>%
summarise(`-log10padj` = mean(`-log10padj`, na.rm = TRUE), .groups = "drop") %>%  # Use mean to avoid duplicates
distinct(LR, !!sym(interactionid_var), .keep_all = TRUE) %>%
pivot_wider(names_from = !!sym(interactionid_var), values_from = `-log10padj`, values_fill = 0) %>%
column_to_rownames("LR") %>%
as.matrix()

significance_matrix <- interaction_df %>%
  dplyr::select(LR, !!sym(interactionid_var), significant) %>%
  pivot_wider(names_from = !!sym(interactionid_var), values_from = significant, values_fill = "Non-significant") %>%
  column_to_rownames("LR") %>%
  as.matrix()


  row_dendrogram <- hclust(dist(dotplot_matrix))

  # hypergeometric_test_LR_comb_df_withpatientcountinfo2 = hypergeometric_test_LR_comb_df_withpatientcountinfo  %>% arrange(desc(celltrek_spot_count),padj)  %>% mutate(padj = ifelse(padj == 0, 1e-300, padj))
  LR_order = rownames(dotplot_matrix)[row_dendrogram$order]

significance_matrix_reord = significance_matrix[LR_order, col_order, drop = FALSE]

dotplot_matrix_reord = dotplot_matrix[LR_order,col_order, drop = FALSE]
dotplot_color_matrix_reord = dotplot_color_matrix[LR_order,col_order, drop = FALSE]

mymat = matrix(0, nrow = nrow(dotplot_matrix_reord), ncol = ncol(dotplot_matrix_reord))  # Empty background
rownames(mymat) = rownames(dotplot_matrix_reord)
colnames(mymat) = colnames(dotplot_matrix_reord)


celltype_freq_by_LR2_sel = subset(celltype_freq_by_LR2, Lig_Rec %in% LR_order)

celltype_freq_by_LR2_sel = as.data.frame(celltype_freq_by_LR2_sel)

rownames(celltype_freq_by_LR2_sel) = celltype_freq_by_LR2_sel$Lig_Rec

celltype_freq_by_LR2_sel$Lig_Rec = NULL

myannot_mat = as.matrix(celltype_freq_by_LR2_sel)

rownames(myannot_mat) = rownames(celltype_freq_by_LR2_sel)

myannot_mat_ord = myannot_mat[LR_order,]

celltype_freq_by_cluster_sel = subset(celltype_freq_by_cluster2, integrated_snn_res.0.9 %in% col_order)

celltype_freq_by_cluster_sel = as.data.frame(celltype_freq_by_cluster_sel)

rownames(celltype_freq_by_cluster_sel) = celltype_freq_by_cluster_sel$`integrated_snn_res.0.9`

celltype_freq_by_cluster_sel$`integrated_snn_res.0.9` = NULL

myannot_mat_cluster = as.matrix(celltype_freq_by_cluster_sel)

myannot_mat_cluster_ord = myannot_mat_cluster[col_order,]
# column_ha = HeatmapAnnotation(foo1 = runif(10), bar1 = anno_barplot(runif(10)))
# row_ha = rowAnnotation(foo2 = runif(10), bar2 = anno_barplot(runif(10)))
# Heatmap(mat, name = "mat", top_annotation = column_ha, right_annotation = row_ha)


myanno_bar_LR = rowAnnotation(Spot_composition_by_LR = anno_barplot(myannot_mat_ord, gp = gpar(fill = htmap_barplot_colors),bar_width = 1, width = unit(6, "cm"))) 

myanno_bar_cluster = HeatmapAnnotation(Spot_composition_by_Cluster = anno_barplot(myannot_mat_cluster_ord, gp = gpar(fill = htmap_barplot_colors),bar_width = 1, width = unit(6, "cm"))) 


# myrow_order = clustered_row_order
dotplot_ht <- ComplexHeatmap::Heatmap(
  mymat,  # Force row order to match dotplot_matrix
  cell_fun = function(j, i, x, y, width, height, fill) {
  grid.rect(x, y, width, height, gp = gpar(col = "black", lty = "dotted"))  
  sig_value <- significance_matrix_reord[i, j]
  significance_border <- ifelse(sig_value == "Significant", 3.5, 0)
  if (dotplot_matrix_reord[i, j] > 0) {
    grid.circle(
      x = x, y = y, 
      r = unit(dot_size_fun(dotplot_matrix_reord[i, j]), "mm"),  
      gp = gpar(fill = col_fun_dotplot(dotplot_color_matrix_reord[i, j]), 
                col = "black",lwd = significance_border)  
    )
  }
},col = c("white"),  
  rect_gp = gpar(col = "black", lty = "dotted"),  
  row_order = LR_order,  # Ensure same row order
  column_order = col_order,
  row_names_side = "right",
  column_names_rot = 45,
  column_names_side = "bottom",
  row_names_gp = grid::gpar(fontsize = 12),
  column_names_gp = grid::gpar(fontsize = 12),
  show_row_dend = TRUE,
  show_heatmap_legend = FALSE,  width = unit(15, "cm"),left_annotation=myanno_bar_LR,
  top_annotation = myanno_bar_cluster)

pdf(file =   paste0('dotplot_ht.',out_prefix,'.pdf'),width =20,height  = 30)
draw(dotplot_ht, heatmap_legend_list = lgd)
dev.off()

}


make_complexheatmap_signalling_withinspot_enrichment(
interaction_df=forhtmap,
interactionid_var='Cluster',
out_prefix='within_spot_enrichment')

