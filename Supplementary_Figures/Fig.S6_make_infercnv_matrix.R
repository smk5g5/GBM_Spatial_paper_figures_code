library(tidyverse)
library(ComplexHeatmap)
library(GenomicRanges)
library(colorspace)
library(circlize)
library(optparse)
library(infercnv)
library(doParallel)
library(foreach)
library(EnrichedHeatmap)
library(patchwork)

set.seed(1234)


# my_matrix <- matrix(c(5, 10, 15), nrow = 1)
# # Set column and row names
# colnames(my_matrix) <- c("A", "B", "C")
# rownames(my_matrix) <- c("Sample1")

# # Print the matrix
# print(my_matrix)

case_when(
    state_id==1 ~ -1,
    state_id==2 ~ -0.5,
    state_id==3 ~ 1,
    state_id==4 ~ 1.5,
    state_id==5 ~ 2,
    state_id==6 ~ 3,
  ))

sg_out <- case_when(
  sg == 'border/transitory' ~ 'border_or_transitory',
  .default = as.character(sg)
)


# cnv_regions = ((unlist(lapply(mcmc@cell_gene, function(x) {
#   length(x$State)})))

# cnv_regions = (unlist(lapply(mcmc@cell_gene, function(x) {
#   x$cnv_regions})))



CNV_per_spot = function(i){

cnv_reg = as.character(mcmc@cell_gene[[i]]$cnv_regions)

all_tumor_spots = colnames(mcmc@expr.data)[mcmc@observation_grouped_cell_indices$malignant]

mystate = mcmc@cell_gene[[i]]$State
myspots = intersect(all_tumor_spots,colnames(mcmc@expr.data)[mcmc@cell_gene[[i]]$Cells])

# observation_grouped_cell_indices = mcmc@reference_grouped_cell_indices$normal
observation_grouped_cell_indices = mcmc@observation_grouped_cell_indices$malignant
# myspots = intersect(colnames(mcmc@expr.data)[observation_grouped_cell_indices],myspots)

other_spots = setdiff(colnames(mcmc@expr.data)[observation_grouped_cell_indices],myspots)

other_matrix =  matrix(rep(3 ,length(other_spots)), nrow = 1)
rownames(other_matrix) = cnv_reg
colnames(other_matrix) = other_spots


my_matrix = matrix(rep(mcmc@cell_gene[[i]]$State ,length(myspots)), nrow = 1)
rownames(my_matrix) = cnv_reg
colnames(my_matrix) = myspots

comb_mat = cbind(other_matrix,my_matrix)

return(comb_mat)
}


# my_CNV_scores = list()

# for(i in 1:length(mcmc@cell_gene)){
# my_CNV_scores[[i]] = CNV_per_spot(i)
# }

Immunefeatures <- c("PTPRC","CD2","CD3D","CD3E","CD3G","CD5","CD7","CD79A",'MS4A1',"CD19")
Mast_cells <- c("GATA2", "MS4A2", "KIT")
Myeloid_cells <- c("FCGR3A", "CD68", "MARCO", "LYZ")
B_lymphocytes <- c("IGHA2", "IGHG3", "IGHM", "CD79A")
NK_cells <- c("KLRD1", "NCAM1", "GNLY", "NKG7")
T_cells <- c("TRAC", "CD3G", "CD3E", "CD3D")

all_immune_feats = unique(c(Immunefeatures,Mast_cells,Myeloid_cells,B_lymphocytes,NK_cells,T_cells))

for(i in names(seurat_obj_list2)){
	seurat_obj_list2[[i]]@meta.data$NormalScore = apply(seurat_obj_list2[[i]]@assays$Spatial[rownames(seurat_obj_list2[[i]]@assays$Spatial) %in% all_immune_feats, ], 2, mean)
	seurat_obj_list2[[i]]@meta.data$NormalScore_sct = apply(seurat_obj_list2[[i]]@assays$SCT[rownames(seurat_obj_list2[[i]]@assays$SCT) %in% all_immune_feats, ], 2, mean)
}


HMM_output_infercnv = read.table('/n/scratch/users/s/sak4832/CNV_spatial/HMM_output_infercnv.txt')

HMM_output_infercnv$Patient = gsub("_infercnv2_clonepubversion_0.01_leiden","",HMM_output_infercnv$V1,perl=TRUE);

HMM_CNV_score_bypatients = list()


for(i in 1:nrow(HMM_output_infercnv)){
mcmc = readRDS(HMM_output_infercnv$V2[i])
mypat = HMM_output_infercnv$Patient[i]
my_CNV_scores <- lapply(1:length(mcmc@cell_gene), CNV_per_spot)
reference_order <- colnames(my_CNV_scores[[1]])

my_CNV_scores_ordered <- lapply(my_CNV_scores, function(mat) {
  mat[, reference_order, drop = FALSE]  # drop = FALSE preserves matrix structure
})

my_CNV_scores_ordered_comb = do.call(rbind,my_CNV_scores_ordered)

my_CNV_scores_ordered_comb_min3 = t(abs(my_CNV_scores_ordered_comb-3))

HMM_CNV_score_bypatients[[mypat]][['CNV_output']] = my_CNV_scores_ordered_comb

HMM_CNV_score_bypatients[[mypat]][['CNV_score']] = my_CNV_scores_ordered_comb_min3
}


normalize_min_max <- function(x) {
    return((x- min(x)) /(max(x)-min(x)))
}
for (i in names(HMM_CNV_score_bypatients)) {
  
  # Compute per-spot CNV score as row sum
  HMM_CNV_score_bypatients[[i]][['perspot_CNV_score']] <- rowSums(HMM_CNV_score_bypatients[[i]][['CNV_score']])
  
  cnv_scores <- HMM_CNV_score_bypatients[[i]][['perspot_CNV_score']]
  
  median_CNV_score <- median(cnv_scores)
  SD_CNV_score <- sd(cnv_scores)
  
  # Classify spots
  HMM_CNV_score_bypatients[[i]][['perspot_CNV_classification']] <- ifelse(
    cnv_scores == 0, "NO_CNV",
    ifelse(cnv_scores < (median_CNV_score - 1 * SD_CNV_score), "Low_CNV",
      ifelse(cnv_scores > (median_CNV_score + 1 * SD_CNV_score), "high_CNV", "intermediate_CNV")
    )
  )
  
  # Scale CNV scores
  HMM_CNV_score_bypatients[[i]][['perspot_CNV_score_scaled']] <- normalize_min_max(cnv_scores)
}

# all(Cells(seurat_obj_list2[['B172.01']]) %in% mydf_list[[1]]$barcodes_seurat)
# sd(normalize_min_max(HMM_CNV_score_bypatients[[i]][['perspot_CNV_score']]))
for (i in names(HMM_CNV_score_bypatients)) {

  mydf <- data.frame(
    barcodes_infercnv = names(HMM_CNV_score_bypatients[[i]][['perspot_CNV_classification']]),
    CNV_class = unname(HMM_CNV_score_bypatients[[i]][['perspot_CNV_classification']])
  )

  mydf2 <- data.frame(
    barcodes_infercnv = names(HMM_CNV_score_bypatients[[i]][['perspot_CNV_score']]),
    perspot_CNV_score = unname(HMM_CNV_score_bypatients[[i]][['perspot_CNV_score']])
  )

  mydf3 <- data.frame(
    barcodes_infercnv = names(HMM_CNV_score_bypatients[[i]][['perspot_CNV_score_scaled']]),
    perspot_CNV_score_scaled = unname(HMM_CNV_score_bypatients[[i]][['perspot_CNV_score_scaled']])
  )

  mydf_combined <- Reduce(function(...) merge(..., by = "barcodes_infercnv"),
                          list(mydf, mydf2, mydf3))

  if(i!='B178'){

  split_bars <- do.call(rbind, strsplit(mydf_combined$barcodes_infercnv, "_"))
  mydf_combined$barcodes_seurat <- split_bars[, 1]
  mydf_combined$barcodes_grp <- split_bars[, 2]

  mydf_list <- split(mydf_combined, mydf_combined$barcodes_grp)

  select_samp <- names(seurat_obj_list2)[grep(i, names(seurat_obj_list2))]

  mytst <- lapply(mydf_list, function(x) {
    mybars <- x$barcodes_seurat
    sapply(select_samp, function(y) {
      all(mybars %in% Cells(seurat_obj_list2[[y]]))
    })
  })

for(myspl in names(mydf_list)){
	mysamp = names(which(mytst[[myspl]]==T))
	index = match(Cells(seurat_obj_list2[[mysamp]]),mydf_list[[myspl]]$barcodes_seurat)
	seurat_obj_list2[[mysamp]]$perspot_CNV_classification = mydf_list[[myspl]]$CNV_class[index]
	seurat_obj_list2[[mysamp]]$perspot_CNV_score = mydf_list[[myspl]]$perspot_CNV_score[index]
	seurat_obj_list2[[mysamp]]$perspot_CNV_score_scaled = mydf_list[[myspl]]$perspot_CNV_score_scaled[index]
}

}else{
	index = match(Cells(seurat_obj_list2[[mysamp]]),mydf_combined$barcodes_seurat)
	seurat_obj_list2[[mysamp]]$perspot_CNV_classification = mydf_combined$CNV_class[index]
	seurat_obj_list2[[mysamp]]$perspot_CNV_score = mydf_combined$perspot_CNV_score[index]
	seurat_obj_list2[[mysamp]]$perspot_CNV_score_scaled = mydf_combined$perspot_CNV_score_scaled[index]
}
}



# for(i in names(seurat_obj_list2)){
# 	seurat_obj_list2[[i]]@meta.data$NormalScore = apply(seurat_obj_list2[[i]]@assays$Spatial[rownames(seurat_obj_list2[[i]]@assays$Spatial) %in% all_immune_feats, ], 2, mean)
# 	seurat_obj_list2[[i]]@meta.data$NormalScore_sct = apply(seurat_obj_list2[[i]]@assays$SCT[rownames(seurat_obj_list2[[i]]@assays$SCT) %in% all_immune_feats, ], 2, mean)
# }

  # plots[[i]] <- SpatialFeaturePlot(object = GBM_spatial_sub, features = i) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +labs(title = paste0('metabolic flux score ',i))

saveRDS(seurat_obj_list2,paste0('CNV_score_seurat_obj_list2.rds'))

lapply(names(seurat_obj_list2),

library(gridExtra)

myCNV_cols = c('gray88',brewer.pal(n=9,name='Reds')[1],brewer.pal(n=9,name='Reds')[5],brewer.pal(n=9,name='Reds')[9])
names(myCNV_cols) = c('NO_CNV','Low_CNV',"intermediate_CNV","high_CNV")
library(RColorBrewer)

supplementary_figures_clustering.R:  mydmpl = SpatialDimPlot(seurat_obj_list2[[i]], image.alpha = 0.5, group.by= mypath, pt.size.factor=1, stroke=0.0, crop=T) + scale_fill_manual(values = mypath_cols) + ggtitle(paste0(i,' Pathologist_annotation : ',mypath)) + theme(plot.title = element_text(hjust = 0.5)) +  guides(fill = guide_legend(override.aes = list(size = 15)))

pdf(file = 'Spatial_Dimplots_CNV_scores.pdf',width = 20,height = 20)
plots=list();
# paste0('Spatial_Dimplots_clusters_allsamples_rpcaclusters_integbypatient_onepageperplot_filteredfinal_withlegend',j,'.',date,'.pdf'),
for (i in names(seurat_obj_list2)){
  plots[[i]][['perspot_CNV_score']] <- SpatialFeaturePlot(object = seurat_obj_list2[[i]], features = 'perspot_CNV_score') + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +labs(title = paste0('perspot_CNV_score ',i))
  plots[[i]][['perspot_CNV_score_scaled']] <- SpatialFeaturePlot(object = seurat_obj_list2[[i]], features = 'perspot_CNV_score_scaled') + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +labs(title = paste0('perspot_CNV_score_scaled ',i))
  plots[[i]][['Immune_score']] <- SpatialFeaturePlot(object = seurat_obj_list2[[i]], features = 'NormalScore') + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +labs(title = paste0('Immune_score ',i))
  plots[[i]][['perspot_CNV_classification']] = SpatialDimPlot(seurat_obj_list2[[i]], image.alpha = 0.8, group.by= 'perspot_CNV_classification', pt.size.factor=1.2, stroke=0, crop=TRUE,cols =myCNV_cols) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+  guides(fill = guide_legend(override.aes = list(size = 10)))

ml <- marrangeGrob(plots[[i]], nrow=2, ncol=2)
print(ml)
}
dev.off()

# # plot_grid(p1, p2, p1, ncol = 3)
# ml <- marrangeGrob(spatial_dm_list, nrow=1, ncol=1)
# ggsave(paste0('Spatial_Dimplots_clusters_allsamples_rpcaclusters_integbypatient_onepageperplot_filteredfinal_withlegend',j,'.',date,'.pdf'), ml, width = 15, height = 15, units="in",limitsize = FALSE)
# lapply(names(seurat_obj_list2),function(x) {
# grep('CNV',colnames(seurat_obj_list2[[x]]@meta.data),ignore.case=T,value=T)
# })

