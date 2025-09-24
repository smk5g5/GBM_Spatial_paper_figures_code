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
