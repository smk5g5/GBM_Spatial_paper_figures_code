### Spatial transcriptomics GBM Paper (visium) figures code

![2025_CGM_GBM_spatial_V1](https://github.com/user-attachments/assets/b1fea9dc-291d-42cd-8451-abaacff07b8b)

## Overview
This is the repo for the manuscript "Mapping the spatial architecture of glioblastoma from core to edge delineates niche-specific tumor cell states and intercellular interactions"

## Repo structure
# Figure 2
 - Figure 2 B-D & F-H (**Figure_2_BtoD_FtoH_final_figure_celltrekwithintp_cellstatesplot.R**) : 
	 - (B) Cell mapping results for representative core tumor sections showing dense 		tumor cellularity, vascularization, and pockets of immune infiltration.
	 -  (C) Cell mapping results for representative edge sections showing minimal tumor content and abundant neurons or glia. 
	 - (D) Cell mapping results for representative transitional 5-ALA+ sections showing dense tumor cellularity,  vascularization, and pockets of immune infiltration. 
	 - (F) Spatial tumor cell states in core tumor slices showing abundant Hypoxic (S4) and Mesenchymal (S3) tumor cells.
	 - (G) Spatial tumor states in edge samples, showing enrichment for Proliferating (S2) and Stem-like (S5) states in regions of high neuronal/glial content. 
	 - (H) Spatial tumor cell states in transitional 5-7866ALA+ tumor slices showing a gradient from core Mesenchymal/Hypoxic to edge (Proliferative/Stem-like) states.
- color_panel_list.rds (color panel list)
- final_genemodules_reduced_mutualinfo_list.rds (final gene list after filtering by mutual information)
- spatial_data_metadata.txt (spatial metadata used in the heatmap)
- Figure 2 E (**Figure_2E_binary_heatmap_final_modules.R**) : 
	- (E) Consensus clustering of marker genes for spatial tumor cell states derived from mapped tumor cells.

# Figure 3
 - Figure 3A (**Figure_3A_celltype_enrichment_visium_doc.R**) : Integrates spatial & single-cell data, performs enrichment and visualization of tumor microenvironment cell states. 
 - Figure 3C (**Figure3C_Core_edge_tumorstate_enrichment.R**) code for enrichment comparison of spatial tumor states between domains 
	 - Fisher exact test for enrichment of Tumor states between domains. Dot color represents the bonferroni corrected p-value of enrichments and dot-size is  the fold change of proportion of tumor states (w.r.t to tumor cells) between domains.
 - Figure 3G (**Figure3G_spot_distances_calculate.R**) Spatial proximity map of niches, where edge width is proportional to Weighted Mean Adjacency and edge color is proportional to the fraction of samples in which each niche pair is adjacent.
 - adjacent_cluster_barcodes.rds : barcodes for each cluster which are adjacent to each other (i.e. at a minimum distance between 2 visium spots which is 100 microns).
 - dotplot_items_list_march13.2025.rds : dotplot items which contain cluster group assignments i.e. edge/core/intermediate/Non-specific
 - final_cluster_colors.rds (final cluster color palette)
 - spots_to_cells_celltrek_annotations.rds (celltrek cells to spot mapping) list of cells (from the scrnaseq reference for each spot)
 - Clusteradjacency.cys (Cytoscape object for generating figure 3G)

# Figure 4
   - Figure 4A (**Figure_4A.R**) : Enrichment dotplot of significant contact-dependent (within-spot) ligand receptor interactions in visium data. Size of dot indicates enrichment score, bold circle around spot indicates significance based on adjusted p.values. Row Barplots show celltype diversity for the said LR whereas column barplots indicate  celltype diversity in the cluster. The clusters are ordered from edge rich clusters to core-rich clusters.

# Supplementary Figures
   - Figure S12 (** Figure_S12_make_heatmap_LR_scref_func.R**) 
	   - Ligand receptor interactions in single-cell RNA-seq (reference data) recapitulate findings in Xenium and Visium datasets.
		   - (A) LR interactions where Tumor cells act as senders and other cell types as receivers, grouped by broad categories indicated in the top annotation. Each row corresponds to an LR pair, and each column represents a sender–receiver combination as inferred from the single-cell reference, which includes spiked-in normal tissue data to improve representation of non-tumor cell types. Cell color in the main heatmap indicates the strength of communication (communication probability), and asterisks denote the statistical significance of the interaction. The binary heatmap on the far right of each panel indicates whether each LR interaction was also observed in the Visium and Xenium datasets (defined as present in at least three patients per sender–receiver group for Xenium, and present in at least three patients in Visium as shown in 
Fig. 4, S10, and S11) and in 10x single-cell RNA-seq (defined as having a communication probability greater than the overall mean). 
		- (B) LR interactions where non-tumor cell types act as senders and tumor cells as receivers. 
		- (C) LR interactions between Myeloid and Lymphoid cells. These results highlight consistent intercellular communication patterns across single-cell and spatial transcriptomic platforms and provide insights into potential biological mechanisms.

- Figure S6 & S7
	-  **Fig.S6_S7_Run_infercnv_subclusters.R** : Run Infercnv on spots against normal brain samples for each of the 9 patients
	- **Fig.S6_make_infercnv_matrix.R** :  Make infercnv matrix from the mcmc object of states and calculate a CNV score for each spot
	- **Fig.S6_S7_make_CNV_metaclusters.R** : Converge the CNV subclusters as identified for each patient into metaclusters based on matrix similarity
	- **Fig.S6_S7_CNV_values_matrix_t2_heatmaps.R** : Make complexheatmap with all different metadata for each patient
