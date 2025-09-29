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
	 - (H) Spatial tumor cell states in transitional 5-ALA+ tumor slices showing a gradient from core Mesenchymal/Hypoxic to edge (Proliferative/Stem-like) states.
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


# Figures 5-7 (Xenium data analysis Xenium_Code.R)
-	**Figure 5: Xenium-based Analysis of GBM Tumor Samples**.  
	-	(A) UMAP visualization of all cells from primary cohort samples. 
	-	(B) Spatial plots of GBM026I and GBM030 with cell types colored as per Fig. 5A. 
	-	(C) Hierarchical clustering of tumor clusters from both cohorts (primary and secondary) based upon expression of Visium-based tumor signatures. 
	-	(D) Hierarchical clustering of BANKSY clusters from both cohorts (primary and secondary) based upon normalized values of cell type frequency. 
	-	(E) Spatial plots of GBM026I and GBM030 with niches colored. 
	-	(F) Dot plot of enrichment scores for each cell type per respective niche. The enrichment score was calculated by dividing the cell type frequency per niche by the cell type frequency within the entire data set and signified by the size of each circle. Significance was calculated using the hypergeometric test and denoted by a black border around each respective circle. 
	-	(G) Node-edge plot visualizing the average distance between niches across all samples (in both cohorts). Distances between niches were calculated on a per sample basis and then averaged across samples if applicable. The number of samples with a given niche-niche distance is signified by the circled number. In (B) and (E), specific regions of interest (ROIs) are highlighted by dashed white boxes with the respective ROI matched by the marked number in the top left or bottom right corner.
	
-	**Figure 6: Ligand-Receptor Analysis of Xenium Samples**
	-	(A) Interaction network illustrating unique ligand receptor interactions between general cell types. Only unique ligand receptors detected in more than 3 samples were visualized.The color, opacity, and thickness of arrows are proportional to the number of unique LR interactions connecting any two cell types.
	-	(B and C) Heatmaps of the (B) unique and (C) total number of ligand receptor interactions between general cell types.
	-	(D) Ligand-receptor interactions of tumor cells (Hypoxic: S2 or S4, Stem-like: S1 or S5) as senders and myeloid, neuronal/glial, and vascular cells as receivers. Ligand-receptor pairs unique to a given sender-receiver pair were italicized and labeled with an asterisk.
	-	(E) Interaction network of cell type sender-receiver pairs of the ADM:CALCRL ligand-receptor pair. The color and density of arrows were proportional to the number of samples the ADM:CALCRL interaction was detected as significant for a given sender-receiver pair.
	-	(F) Heatmap of ADM and CALCRL expression per cell type per niche.
	-	(G) Expression of ADM and CALCRL based upon localization in tumor region. Data derived from IVY-GAP.
	
-	**Figure 7: CD8+GZMK+ T cells preferentially localize with LYVE1+CD163+  myeloid cells in vascular niches**
	-	(A) Top 6 cell types enriched within 20 micron radius of CD8+GZMK+ T cells. Enrichment score was calculated by dividing the frequency of a target cell type within a given radius of a root state (CD8+GZMK+ T cells) by the frequency of the target cell type within the sample. Each dot represents a patient sample.
	-	(B) Distance to the nearest target myeloid cell of each CD8+.GZMK+ T cell separated by localization within a vascular or non-vascular niche. Each dot represents a single CD8+.GZMK+ T cell, with all CD8+.GZMK+ T cells across all samples visualized.
	-	(C) Distances from Fig. 6B averaged on a per sample basis. Each dot represents an individual sample. Significance was calculated using the Wilcoxon matched-pairs signed rank test.
	-	(D) Frequency of CD8+.GZMK+ T cells in the vascular niche given proximity to the target myeloid cell (either within 30 μm or further than 30 μm).
	-	(E) Frequency of CD8+GZMK+ T cells within 30 μm of target myeloid cell given presence of T cell in vascular or non-vascular niche.
	-	(F) Correlation of GZMK and LYVE1 expression based upon localization in vascular or non-vascular niche. Data derived from IVY-Gap.

# Supplementary Figures
   - Figure S12 (**Figure_S12_make_heatmap_LR_scref_func.R**) 
	   - Ligand receptor interactions in single-cell RNA-seq (reference data) recapitulate findings in Xenium and Visium datasets.
		   - (A) LR interactions where Tumor cells act as senders and other cell types as receivers, grouped by broad categories indicated in the top annotation. Each row corresponds to an LR pair, and each column represents a sender–receiver combination as inferred from the single-cell reference, which includes spiked-in normal tissue data to improve representation of non-tumor cell types. Cell color in the main heatmap indicates the strength of communication (communication probability), and asterisks denote the statistical significance of the interaction. The binary heatmap on the far right of each panel indicates whether each LR interaction was also observed in the Visium and Xenium datasets (defined as present in at least three patients per sender–receiver group for Xenium, and present in at least three patients in Visium as shown in 
Fig. 4, S10, and S11) and in 10x single-cell RNA-seq (defined as having a communication probability greater than the overall mean). 
		- (B) LR interactions where non-tumor cell types act as senders and tumor cells as receivers. 
		- (C) LR interactions between Myeloid and Lymphoid cells. These results highlight consistent intercellular communication patterns across single-cell and spatial transcriptomic platforms and provide insights into potential biological mechanisms.

- Figure S6 (Fig. S6. **CNV Spot classification**.Spots were classified into Low/Medium/High malignancy based on the frequency and magnitude of copy number alterations. High malignancy spots (i.e. Spots with high CNV scores) usually reside in the Core-rich domains/Niches whereas Low malignancy spots (i.e. Spots with low CNV scores) are in Edge-rich domains as shown in the spatial plots. **Spot malignancy** proportions align well with cell type proportions across clusters, with regions enriched in tumor cells exhibiting correspondingly higher CNV scores as shown in the bar plots.) & S7 (CNV subclusters along with **spot malignancy/CNV score** are shown in greater detail by chromosome for each patient in a **heatmap**)
	-  **Fig.S6_S7_Run_infercnv_subclusters.R** : Run Infercnv on spots against normal brain samples for each of the 9 patients
	- **Fig.S6_make_infercnv_matrix.R** :  Make infercnv matrix from the mcmc object of states and calculate a CNV score for each spot
	- **Fig.S6_S7_make_CNV_metaclusters.R** : Converge the CNV subclusters as identified for each patient into metaclusters based on matrix similarity
	- **Fig.S6_S7_CNV_values_matrix_t2_heatmaps.R** : Make complexheatmap with all different metadata for each patient
