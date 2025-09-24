### Spatial transcriptomics GBM Paper (visium) figures code

![2025_CGM_GBM_spatial_V1](https://github.com/user-attachments/assets/b1fea9dc-291d-42cd-8451-abaacff07b8b)

## Overview
This is the repo for the manuscript "Mapping the spatial architecture of glioblastoma from core to edge delineates niche-specific tumor cell states and intercellular interactions"

## Repo structure
# Figure 2
	 - Figure 2 B-D & F-H (Figure_2E_binary_heatmap_final_modules.R) : 
		 - (B) Cell mapping results for representative core tumor sections showing dense 		tumor cellularity, vascularization, and pockets of immune infiltration.
		 -  (C) Cell mapping results for representative edge sections showing minimal tumor content and abundant neurons or glia. 
		 - (D) Cell mapping results for representative transitional 5-ALA+ sections showing dense tumor cellularity,  vascularization, and pockets of immune infiltration. 
		 - (F) Spatial tumor cell states in core tumor slices showing abundant Hypoxic (S4) and Mesenchymal (S3) tumor cells.
		 - (G) Spatial tumor states in edge samples, showing enrichment for Proliferating (S2) and Stem-like (S5) states in regions of high neuronal/glial content. 
		 - (H) Spatial tumor cell states in transitional 5-7866ALA+ tumor slices showing a gradient from core Mesenchymal/Hypoxic to edge (Proliferative/Stem-like) states.
	- color_panel_list.rds (color panel list)
	- final_genemodules_reduced_mutualinfo_list.rds (final gene list after filtering by mutual information)
	- spatial_data_metadata.txt (spatial metadata used in the heatmap)
	- Figure 2 E (Figure_2E_binary_heatmap_final_modules.R) : 
		- (E) Consensus clustering of marker genes for spatial tumor cell states derived from mapped tumor cells.
