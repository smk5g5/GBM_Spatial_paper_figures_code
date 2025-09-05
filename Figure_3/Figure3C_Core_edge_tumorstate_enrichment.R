# ============================================================================
# Figure 3C code for enrichment comparison of spatial tumor states between domains
# Author: [Saad Khan]   |   Date: 2025-09-05
# Description: Fisher exact test for enrichment of Tumor states between domains. 
# Dot color represents the bonferroni corrected p-value of enrichments
# and dot-size is  the fold change of proportion of tumor states (w.r.t to tumor cells) between domains.
# ============================================================================
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
library(CellTagR)
library(reshape2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(CellTrek)
library(ggdendro)
library(cowplot)
library(tidyverse)
library(ggtree) 
# install with devtools::install_github("YuLab-SMU/ggtree") (newer than bioconductor version)
library(patchwork)
library(gridExtra)
library(gtools)

set.seed(12345L) 

date = gsub("2025-","25",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

######################################################################################
# Load preprocessed data
#output from celltrek spots_to_cells_celltrek : spots_to_cells_celltrek_annotations
#rpca integrated visium spatial object (all samples) : merged.patient.combined.integ_sub
############################################################


spots_to_cells_celltrek_annotations = readRDS('/n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/cluster_distances/spots_to_cells_celltrek_annotations.rds')

merged.patient.combined.integ_sub = readRDS('/n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/cluster_distances/merged.patient.combined.integ_sub.rds')

############################################################
# Extract relevant metadata from integrated Seurat object
############################################################
rpca_meta = merged.patient.combined.integ_sub@meta.data[c("Sample","Patient","integrated_snn_res.0.1","integrated_snn_res.0.3","integrated_snn_res.0.5","integrated_snn_res.0.7","integrated_snn_res.0.9","integrated_snn_res.1.2")]

rpca_meta$cellname = rownames(rpca_meta)

rpca_meta$cellname2 = str_split_i(rpca_meta$cellname,'_',2)

rpca_meta$sample = str_split_i(rpca_meta$cellname,'_',1)

rpca_meta_list = split(rpca_meta,rpca_meta$sample)

# use https://github.com/Petti-Lab/GBM_Spatial_paper_figures_code/blob/main/Figure_3A/dotplot_items_list_march13.2025.rds
#groupdf = dotplot_items$group_df
groupdf = readRDS('groupdf2.rds')

############################################################
# Merge spot-to-cell annotations with metadata and group info
############################################################

merge_by_group_celltrekspots = lapply(names(rpca_meta_list), function(x) {
	spots_to_celltrek_to_rpca = merge(spots_to_cells_celltrek_annotations[[x]]$broad_celltypes_per_spot,rpca_meta_list[[x]],by.x='spot_id',by.y='cellname2') %>% left_join(groupdf,by = c("integrated_snn_res.0.9"))
	spots_to_celltrek_to_rpca_summary = spots_to_celltrek_to_rpca %>% group_by(group) %>% summarize(
	  M1 = sum(M1),
	  M2 = sum(M2),
	  M3 = sum(M3),
	  M4 = sum(M4),
	  M5 = sum(M5))
	return(spots_to_celltrek_to_rpca_summary)
})

names(merge_by_group_celltrekspots) = names(rpca_meta_list)

############################################################
# Annotate patient IDs for each dataset
############################################################
merge_by_group_celltrekspots_pat = lapply(names(merge_by_group_celltrekspots), function(x) {
	merge_by_group_celltrekspots[[x]]$Patient = str_split_i(x,'\\.',1)
	return(merge_by_group_celltrekspots[[x]])
	})

names(merge_by_group_celltrekspots_pat) = names(merge_by_group_celltrekspots)


############################################################
# Run Fisher exact test for each tumor state and group
############################################################

fisher_exact_test_for_enrichment <- function(df) {
	states <- c("M1", "M2","M3","M4","M5")
	groups <- df$group

	results <- list()

	for (state in states) {
	  for (g1 in groups) {
	    for (g2 in groups) {
	      if (g1 != g2) {
	        # Build contingency
	        x1 <- df %>% filter(group == g1) %>% pull(!!sym(state))
	        x2 <- df %>% filter(group == g2) %>% pull(!!sym(state))
	        
	other1 <- df %>%
	  dplyr::filter(group == g1) %>%
	  dplyr::select(-group, -!!sym(state)) %>%
	  rowSums()

	other2 <- df %>%
	  dplyr::filter(group == g2) %>%
	  dplyr::select(-group, -!!sym(state)) %>%
	  rowSums()
	        
	        tbl <- matrix(c(x1, other1, x2, other2), nrow = 2, byrow = TRUE)
	        rownames(tbl) <- c(g1, g2)
	        colnames(tbl) <- c(state, "other")
	        
	        # fisher_p <- fisher.test(tbl)$p.value
	        fisher_p <- fisher.test(tbl , alternative = "greater")$p.value
	        
	        results[[paste(state, g1, "vs", g2)]] <- list(
	          contingency = tbl,
	          p_value = fisher_p
	        )
	      }
	    }
	  }
	}

	mydf_fisher = data.frame(Tumor_states = str_split_i(names(results),' ',1),Comparisons=sapply(str_split(names(results), ' '), function(x) paste(x[-1], collapse = ' ')),p_value=unlist(lapply(names(results), function(x) {return(results[[x]]$p_value)})))
	mydf_fisher$padj = p.adjust(mydf_fisher$p_value,method='bonferroni',n=nrow(mydf_fisher))
	return(mydf_fisher)
}

# Build fold_change dataframe

fold_change_df_enrichment <- function(df) {

	fold_change_df <- data.frame()

	groups <- unique(df$group)
	states <- c("M1", "M2", "M3", "M4", "M5")

	for (state in states) {
	  for (i in seq_along(groups)) {
	    for (j in seq_along(groups)) {
	      if (i != j) {
	        g1 <- groups[i]
	        g2 <- groups[j]
	        
	        # Compute total counts for state
	        state_g1 <- sum(df %>% dplyr::filter(group == g1) %>% pull(!!sym(state)))
	        state_g2 <- sum(df %>% dplyr::filter(group == g2) %>% pull(!!sym(state)))
	        
	        # Compute total counts across ALL states per group
	        total_g1 <- sum(df %>% dplyr::filter(group == g1) %>% dplyr::select(starts_with("M")) %>% as.matrix())
	        total_g2 <- sum(df %>% dplyr::filter(group == g2) %>% dplyr::select(starts_with("M")) %>% as.matrix())
	        
	        # Proportions
	        prop_g1 <- state_g1 / total_g1
	        prop_g2 <- state_g2 / total_g2
	        
	        # Fold Change (you can take log2 later if preferred)
	        fc <- prop_g1 / prop_g2
	        
	        fold_change_df <- rbind(fold_change_df, data.frame(
	          Tumor_state = state,
	          group1 = g1,
	          group2 = g2,
	          proportion_group1 = prop_g1,
	          proportion_group2 = prop_g2,
	          fold_change = fc
	        ))
	      }
	    }
	  }
	}
return(fold_change_df)
}


df_by_patient = do.call(rbind,merge_by_group_celltrekspots_pat) %>% group_by(group,Patient) %>% summarize(
	  M1 = sum(M1),
	  M2 = sum(M2),
	  M3 = sum(M3),
	  M4 = sum(M4),
	  M5 = sum(M5))  %>%  ungroup()  %>%  group_by(Patient)  %>%  group_split()


pat_names = unlist(lapply(df_by_patient, function(x) {
		return(unique(x$Patient))
	}))


df_by_patient_mod = lapply(df_by_patient, function(x) {
		x$Patient = NULL
		return(x)
	})

names(df_by_patient_mod) = pat_names

df_by_patient_long = do.call(rbind,df_by_patient)

df_by_patient_mod_list = lapply(names(df_by_patient_mod),function(x) {
	mymt = melt(df_by_patient_mod[[x]])
	colnames(mymt) = c('group','state','value')
	mymt$Patient = x
	return(mymt)
	})

df_by_patient_long = do.call(rbind,df_by_patient_mod_list)

mydf = do.call(rbind,merge_by_group_celltrekspots) %>% group_by(group) %>% summarize(
	  M1 = sum(M1),
	  M2 = sum(M2),
	  M3 = sum(M3),
	  M4 = sum(M4),
	  M5 = sum(M5))

#run fisher_exact_test_for_enrichment across patients
overl_all_fisher_exact_test_results = fisher_exact_test_for_enrichment(df=mydf)

fold_change_df = fold_change_df_enrichment(df=mydf)

fisher_exact_test_results_bypat = lapply(names(df_by_patient_mod),function(x) { fisher_exact_test_for_enrichment(df=df_by_patient_mod[[x]]) })

names(fisher_exact_test_results_bypat) = names(df_by_patient_mod)


overl_all_fisher_exact_test_results <- overl_all_fisher_exact_test_results %>%
  mutate(sig_label = case_when(
    padj < 0.001 ~ "***",
    padj < 0.01 ~ "**",
    padj < 0.05 ~ "*",
    TRUE ~ ""
  ))

merged_df <- overl_all_fisher_exact_test_results %>%
  separate(Comparisons, into = c("group1", "vs", "group2"), sep = " ") %>%
  dplyr::select(-vs) %>%
  left_join(fold_change_df, by = c("Tumor_states" = "Tumor_state", "group1", "group2"))

merged_df$Tumor_state = gsub("M","S",merged_df$Tumor_states,perl=TRUE);

merged_df$group1_latest = gsub("transition","Interface",merged_df$group1,perl=TRUE);
merged_df$group2_latest = gsub("transition","Interface",merged_df$group2,perl=TRUE);

library(gtools)

max_log_val <- max(log10(merged_df$padj + 1), na.rm = TRUE)

# c("#023FA5","#BEC1D4","#C6909A","#8E063B")

# pdf('Compare_cell_states_v5.pdf',width = 11,height = 10)

merged_df <- merged_df %>%
  mutate(significant = ifelse(padj < 0.05, TRUE, FALSE))

merged_df <- merged_df %>%
  mutate(Comparison = paste(group1_latest, "vs", group2_latest))

merged_df <- merged_df %>%
  mutate(Comparison = factor(Comparison, levels = mixedsort(unique(Comparison))))

custom_order <- c("edge-rich vs non-specific","edge-rich vs Interface","edge-rich vs core-rich","non-specific vs edge-rich","non-specific vs Interface","non-specific vs core-rich","Interface vs edge-rich","Interface vs non-specific", "Interface vs core-rich","core-rich vs edge-rich","core-rich vs non-specific","core-rich vs Interface")

merged_df <- merged_df %>%
  mutate(Comparison = factor(Comparison, levels = custom_order))


############################################################
#Visualize or make the heatmap in figure 3C
############################################################

p <- ggplot(merged_df, aes(x = Comparison, y = Tumor_state)) +
  
  # OUTLINE layer with significance legend!
  geom_point(
    aes(size = fold_change, color = significant),
    shape = 21,
    stroke = 1.5,
    fill = NA
  ) +
  
  # FILL layer (inside dot color)
  geom_point(
    aes(size = fold_change, fill = log10(padj + 1)),
    shape = 21,
    stroke = 0
  ) +
  
  # Color scale for inside fill
  scale_fill_gradientn(
    colors = colorRampPalette(rev(c("#023FA5", "#BEC1D4", "#C6909A", "#8E063B")))(10),
    limits = c(0, max_log_val),
    name = "log10(padj + 1)"
  ) +
  
  # Legend for significance (border)
  scale_color_manual(
    values = c("TRUE" = "black", "FALSE" = "transparent"),
    name = "Significant (padj < 0.05)",
    labels = c("FALSE" = "ns", "TRUE" = "yes")
  ) +
  
  scale_size_continuous(range = c(3, 20), name = "fold_change") +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(linetype = "dotted", size = 1.5, color = "black"),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Dotplot: Fold Change vs Adjusted p-value + Significance Border",
       x = "Group Comparisons", y = "Tumor States")

ggsave(paste0('fisher_exact_test_for_enrichment_', date, '.jpg'),
       p,
      width = 12,height = 10, units = "in", dpi = 600)
