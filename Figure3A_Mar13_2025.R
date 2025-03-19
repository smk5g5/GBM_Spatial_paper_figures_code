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
library(stringr)library("CellTagR")

library(reshape2)

library(clusterProfiler)
library(org.Hs.eg.db)
library("CellTrek")


library(ggdendro)
library(cowplot)
library(tidyverse)
library(ggtree) # install with `devtools::install_github("YuLab-SMU/ggtree")` as you need a version newer than what bioconductor serves
library(patchwork) 
library(gridExtra)
library(gtools)

set.seed(12345L) 

date = gsub("2025-","25",Sys.Date(),perl=TRUE);
date = gsub("-","",date);


merged.patient.combined.integ_sub = readRDS('/n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/cluster_distances/merged.patient.combined.integ_sub.rds')

# spots_to_cells_celltrek_annotations = readRDS('/n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/cluster_distances/spots_to_cells_celltrek_annotations.rds')

single_cell_ref = readRDS('/n/scratch/users/s/sak4832/scrna_seq_ref_celltrek_orig_latestcellstateannotations.rds')

make_broad_celltypes2 <- function(single_cell_ref){
  mydf = data.frame(names(table(single_cell_ref$Celltypes_Aug8_2024)))
  names(mydf) = c('all_celltypes')
  mydf$broad_celltypes = NA
  mydf$broad_celltypes[grep('CD4_',mydf$all_celltypes)] = 'CD4_Tcells'
  mydf$broad_celltypes[grep('CD8_',mydf$all_celltypes)] = 'CD8_Tcells'
  mydf$broad_celltypes[grep('MAIT-like|Proliferating|T_EMRA|T_IFNG_stimulated',mydf$all_celltypes)] = 'Other_Tcells'
  mydf$broad_celltypes[grep('TAM',mydf$all_celltypes)] = 'TAMs'
  mydf$broad_celltypes[grep('Vascular|Pericyte',mydf$all_celltypes)] = 'Vascular'
  mydf$broad_celltypes[grep('Fibroblast',mydf$all_celltypes)] = 'Fibroblast'
  mydf$broad_celltypes[grep('Microglia',mydf$all_celltypes)] = 'Microglia'
  mydf$broad_celltypes[grep('Monocyte',mydf$all_celltypes)] = 'Monocyte'
  mydf$broad_celltypes[grep('Choroid plexus|Ependymal|Ex. Neuron|Inhib. Neuron|Neuron',mydf$all_celltypes)] = 'Neuron'
  mydf$broad_celltypes[grep('DC',mydf$all_celltypes)] = 'DC'
  mydf$broad_celltypes[mydf$all_celltypes=='Oligodendrocyte'] = 'Oligodendrocyte'	
  mydf$broad_celltypes[mydf$all_celltypes=='Oligodendrocyte precursor'] = 'Oligodendrocyte precursor'	
  mydf$broad_celltypes[mydf$all_celltypes=='Committed oligodendrocyte precursor'] = 'Committed oligodendrocyte precursor'	
  mydf$broad_celltypes[mydf$all_celltypes=='Astrocyte'] = 'Astrocyte'	
  mydf$broad_celltypes[mydf$all_celltypes=='B cells'] = 'Bcells'
  mydf$broad_celltypes[mydf$all_celltypes=='NK'] = 'NK/NK-like'
  mydf$broad_celltypes[mydf$all_celltypes=='NK-like'] = 'NK/NK-like'
  mydf$broad_celltypes[mydf$all_celltypes=='M1'] = 'M1'
  mydf$broad_celltypes[mydf$all_celltypes=='M2'] = 'M2'
  mydf$broad_celltypes[mydf$all_celltypes=='M3'] = 'M3'
  mydf$broad_celltypes[mydf$all_celltypes=='M4'] = 'M4'
  mydf$broad_celltypes[mydf$all_celltypes=='M5'] = 'M5'
  mydf2 = data.frame(cellnames=names(single_cell_ref$Celltypes_Aug8_2024),Celltypes_Aug8_2024=unname(single_cell_ref$Celltypes_Aug8_2024))
  
  mydf12 = mydf2 %>% left_join(mydf, by = c("Celltypes_Aug8_2024"="all_celltypes"))

  index = match(Cells(single_cell_ref),mydf12$cellnames)
  single_cell_ref$broad_celltypes = mydf12$broad_celltypes[index]
  return(single_cell_ref)
}

single_cell_ref2 = make_broad_celltypes2(single_cell_ref)

#read all individual seurat sample files here
spatial_inp = read.table('/n/data1/mgh/neuro/petti/lab/Users/khan.saad/RPCA_clustering_and_spots_Ligrec_summary/spatial_input.txt')


seurat_obj_list = list()

for(i in 1:nrow(spatial_inp)){
  visium_obj = readRDS(spatial_inp$V2[i])
  sample_name = names(visium_obj@images)[1]
  seurat_obj_list[[sample_name]] = visium_obj
}

rpca_meta = merged.patient.combined.integ_sub@meta.data[c("Sample","Patient","integrated_snn_res.0.1","integrated_snn_res.0.3","integrated_snn_res.0.5","integrated_snn_res.0.7","integrated_snn_res.0.9","integrated_snn_res.1.2")]

rpca_meta$cellname = rownames(rpca_meta)

rpca_meta$cellname2 = str_split_i(rpca_meta$cellname,'_',2)

rpca_meta$sample = str_split_i(rpca_meta$cellname,'_',1)

rpca_meta_list = split(rpca_meta,rpca_meta$sample)

for(i in names(seurat_obj_list)){
  rpca_meta_list_sub = rpca_meta_list[[i]]
  index = match(Cells(seurat_obj_list[[i]]),rpca_meta_list_sub$cellname2)
  # c("Sample","Patient","integrated_snn_res.0.1","integrated_snn_res.0.3","integrated_snn_res.0.5","integrated_snn_res.0.7","integrated_snn_res.0.9","integrated_snn_res.1.2","Location","Kit")
  seurat_obj_list[[i]]$Patient = rpca_meta_list_sub$Patient[index]
  seurat_obj_list[[i]]$Sample = rpca_meta_list_sub$Sample[index]
  seurat_obj_list[[i]]$integrated_snn_res.0.1 = rpca_meta_list_sub$integrated_snn_res.0.1[index]
  seurat_obj_list[[i]]$integrated_snn_res.0.3 = rpca_meta_list_sub$integrated_snn_res.0.3[index]
  seurat_obj_list[[i]]$integrated_snn_res.0.5 = rpca_meta_list_sub$integrated_snn_res.0.5[index]
  seurat_obj_list[[i]]$integrated_snn_res.0.7 = rpca_meta_list_sub$integrated_snn_res.0.7[index]
  seurat_obj_list[[i]]$integrated_snn_res.0.9 = rpca_meta_list_sub$integrated_snn_res.0.9[index]
  seurat_obj_list[[i]]$integrated_snn_res.1.2 = rpca_meta_list_sub$integrated_snn_res.1.2[index]
}


celltrek_no_interpolation_celltypes_inp = read.table('/n/data1/mgh/neuro/petti/lab/Users/khan.saad/RPCA_clustering_and_spots_Ligrec_summary/celltrek_no_interpolation_celltypes_inp2.txt')

celltrek_no_interpolation_celltypes_inp$sample_name = celltrek_no_interpolation_celltypes_inp$V1
 # gsub("_split_list_by_spot.rds","",celltrek_no_interpolation_celltypes_inp$V1,perl=TRUE);

celltrek_no_interpolation_celltypes_inp$sample_name2 = gsub("-",".",celltrek_no_interpolation_celltypes_inp$sample_name,perl=TRUE);

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

index <- match(celltrek_no_interpolation_celltypes_inp$sample_name,mysample_map$sample_name)

celltrek_no_interpolation_celltypes_inp$sample_map = mysample_map$image_name[index]
celltrek_no_interpolation_celltypes_inp$sample_name_final = celltrek_no_interpolation_celltypes_inp$sample_map
celltrek_no_interpolation_celltypes_inp$sample_name_final[celltrek_no_interpolation_celltypes_inp$sample_map %in% 'WU1221.Core'] = 'WU1221.Core2'
celltrek_no_interpolation_celltypes_inp$sample_name_final[celltrek_no_interpolation_celltypes_inp$sample_map %in% 'WU1221.Edge'] = 'WU1221.Edge1'


get_celltrek_celltypes_from_no_interpolation_reverse_engineer <- function(sample_name,spot_to_cell_list_rds,single_cell_ref,seurat_obj){
  split_list_by_spot = readRDS(spot_to_cell_list_rds)
  formatted_names_scref = make.names(Cells(single_cell_ref))
  scref_meta = single_cell_ref@meta.data[c('Celltypes_Aug8_2024','broad_celltypes')]
  scref_meta$cellnames = make.names(rownames(scref_meta))
  
split_list_by_spot2 = lapply(split_list_by_spot, function(x) {
    x %>% left_join(scref_meta, by = c("Var1"="cellnames"))
  })

split_list_by_spot3 = lapply(names(split_list_by_spot2), function(x) {
  split_list_by_spot2[[x]]$spot_name = x
  return(split_list_by_spot2[[x]])
})

split_list_by_spot_df  <- do.call(rbind, split_list_by_spot3)

broad_celltype_count = split_list_by_spot_df %>% group_by(spot_name) %>% dplyr::count(broad_celltypes)

all_celltype_count = split_list_by_spot_df %>% group_by(spot_name) %>% dplyr::count(Celltypes_Aug8_2024)

broad_celltype_mat = broad_celltype_count %>% dcast(spot_name ~ broad_celltypes)
all_celltype_mat = all_celltype_count %>% dcast(spot_name ~ Celltypes_Aug8_2024)

broad_celltype_mat[is.na(broad_celltype_mat)] = 0

all_celltype_mat[is.na(all_celltype_mat)] = 0

mysel_broad = colnames(broad_celltype_mat)[-1]

mysel_all = colnames(all_celltype_mat)[-1]

broad_celltype_mat$total_cell_count = rowSums(broad_celltype_mat[mysel_broad])

all_celltype_mat$total_cell_count = rowSums(all_celltype_mat[mysel_all])

spot_id_df = data.frame(spot_id=Cells(seurat_obj),spot_name= make.names(Cells(seurat_obj)))

print(head(spot_id_df))

print(head(broad_celltype_mat))

broad_celltype_mat2 = merge(spot_id_df,broad_celltype_mat,by='spot_name')
all_celltype_mat2 = merge(spot_id_df,all_celltype_mat,by='spot_name')
broad_celltype_mat2$spot_name = NULL
all_celltype_mat2$spot_name = NULL

return(list('broad_celltypes_per_spot' = broad_celltype_mat2,'all_celltypes_per_spot'=all_celltype_mat2))
}


spots_to_cells_celltrek_annotations = list()

for(i in 1:nrow(celltrek_no_interpolation_celltypes_inp)){
  mysamp = celltrek_no_interpolation_celltypes_inp$sample_name_final[i]
  spot_to_cell_list_rds = celltrek_no_interpolation_celltypes_inp$V2[i]
  spots_to_cells_celltrek_annotations[[mysamp]] = get_celltrek_celltypes_from_no_interpolation_reverse_engineer(sample_name=mysamp,spot_to_cell_list_rds=spot_to_cell_list_rds,single_cell_ref=single_cell_ref2,seurat_obj = seurat_obj_list[[mysamp]])
}

integ_clus_vec = 'integrated_snn_res.0.9'

merged.patient.combined.integ_sub$integrated_snn_res.0.9 = droplevels(merged.patient.combined.integ_sub$integrated_snn_res.0.9)

sel_clus2 = levels(merged.patient.combined.integ_sub$integrated_snn_res.0.9)

seurat_obj_list2 =  lapply(names(seurat_obj_list), function(x) {
	Idents(seurat_obj_list[[x]]) = integ_clus_vec
	DefaultAssay(seurat_obj_list[[x]]) = 'Spatial'
	seurat_sub = subset(seurat_obj_list[[x]],idents=intersect(unique(Idents(seurat_obj_list[[x]])),sel_clus2))
	print(x)
	print(length(Cells(seurat_sub)))
	print(length(Cells(seurat_obj_list[[x]])))
	return(seurat_sub)
	})

names(seurat_obj_list2) = names(seurat_obj_list)

final_cluster_colors = readRDS('/n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/cluster_distances/final_cluster_colors.rds')


integ_clus_vec = c("integrated_snn_res.0.9")

noofspots_by_cluster <- function(sample_name,spots_to_cells_celltrek_annotations,type_name,seurat_obj_list2,integ_clus){
  all_celltypes_per_spot2 = spots_to_cells_celltrek_annotations[[sample_name]][[type_name]]
  # all_celltypes_per_spot2_melt = melt(all_celltypes_per_spot2,id='spot_id')

  seurat_obj_samp = seurat_obj_list2[[sample_name]]
  seurat_obj_meta = seurat_obj_samp@meta.data[c(integ_clus)]
  seurat_obj_meta$spot_id = rownames(seurat_obj_meta)

  all_celltypes_per_spot2_melt_sub = subset(all_celltypes_per_spot2,spot_id %in% seurat_obj_meta$spot_id)

  melt_df = all_celltypes_per_spot2_melt_sub %>% left_join(seurat_obj_meta, by = c("spot_id"))

  noofspots_by_cluster2 = melt_df %>% group_by(!!sym(integ_clus)) %>% summarize(spot_count=n())

  return(noofspots_by_cluster2)
}

get_spot_count_bycluster <- function(spots_to_cells_celltrek_annotations,type_name,integ_clus_vec,seurat_obj_list2=seurat_obj_list2){

integ_clus_vec_spot_count = list()

for(integs in integ_clus_vec){
for(sample_name in names(spots_to_cells_celltrek_annotations)){
integ_clus_vec_spot_count[[integs]][[sample_name]] = noofspots_by_cluster(sample_name=sample_name,spots_to_cells_celltrek_annotations,type_name=type_name,seurat_obj_list2=seurat_obj_list2,integ_clus=integs)
}
}


integ_clus_vec_spot_count_combined = lapply(names(integ_clus_vec_spot_count), function(x) {
  mytemp = do.call(rbind,integ_clus_vec_spot_count[[x]])
  rownames(mytemp) = NULL
  return(mytemp)
})

names(integ_clus_vec_spot_count_combined) = names(integ_clus_vec_spot_count)

integ_clus_vec_spot_count_combined_sumbyclus = lapply(names(integ_clus_vec_spot_count_combined), function(x) {
  spot_count_rpca_clusters = integ_clus_vec_spot_count_combined[[x]] %>% group_by(!!sym(x)) %>% summarize(spot_count_bycluster=sum(spot_count))
})

names(integ_clus_vec_spot_count_combined_sumbyclus) = names(integ_clus_vec_spot_count_combined)

return(integ_clus_vec_spot_count_combined_sumbyclus)
}

melt_by_samples_sub <- function(sample_name,spots_to_cells_celltrek_annotations,type_name,seurat_obj_list2,integ_clus){
  
  all_celltypes_per_spot2 = spots_to_cells_celltrek_annotations[[sample_name]][[type_name]]
  all_celltypes_per_spot2_melt = melt(all_celltypes_per_spot2,id='spot_id')

  seurat_obj_samp = seurat_obj_list2[[sample_name]]
  seurat_obj_meta = seurat_obj_samp@meta.data[c(integ_clus)]
  seurat_obj_meta$spot_id = rownames(seurat_obj_meta)

all_celltypes_per_spot2_melt_sub = subset(all_celltypes_per_spot2_melt,spot_id %in% seurat_obj_meta$spot_id)

  melt_df = all_celltypes_per_spot2_melt_sub %>% left_join(seurat_obj_meta, by = c("spot_id"))


  celltype_count_rpca_clusters = melt_df %>% group_by(variable,!!sym(integ_clus)) %>% summarize(count=sum(value))
  celltype_count_rpca_clusters_sub = subset(celltype_count_rpca_clusters,variable != 'total_cell_count')
  celltype_count_rpca_clusters_sub$variable = droplevels(celltype_count_rpca_clusters_sub$variable)
  celltype_count_rpca_clusters_sub$sample_name = sample_name
  return(celltype_count_rpca_clusters_sub)
}

get_celltype_by_integclus <- function(spots_to_cells_celltrek_annotations,type_name,integ_clus_vec,seurat_obj_list2=seurat_obj_list2){

integ_clus_vec_celltypes = list()

for(integs in integ_clus_vec){
for(sample_name in names(spots_to_cells_celltrek_annotations)){
integ_clus_vec_celltypes[[integs]][[sample_name]] = melt_by_samples_sub(sample_name=sample_name,spots_to_cells_celltrek_annotations,type_name=type_name,seurat_obj_list2=seurat_obj_list2,integ_clus=integs)
}
}


integ_clus_vec_celltypes_combined = lapply(names(integ_clus_vec_celltypes), function(x) {
  mytemp = do.call(rbind,integ_clus_vec_celltypes[[x]])
  rownames(mytemp) = NULL
  return(mytemp)
})

names(integ_clus_vec_celltypes_combined) = names(integ_clus_vec_celltypes)

integ_clus_vec_celltypes_combined_sumbyclus = lapply(names(integ_clus_vec_celltypes_combined), function(x) {
  celltype_count_rpca_clusters = integ_clus_vec_celltypes_combined[[x]] %>% group_by(variable,!!sym(x)) %>% summarize(celltype_count_by_cluster=sum(count))
})

names(integ_clus_vec_celltypes_combined_sumbyclus) = names(integ_clus_vec_celltypes_combined)

return(integ_clus_vec_celltypes_combined_sumbyclus)
}

get_tumor_fraction_by_cluster <- function(spots_to_cells_celltrek_annotations,type_name,integ_clus_vec,seurat_obj_list2=seurat_obj_list2,sample_name=sample_name){

mysel_celltrek = spots_to_cells_celltrek_annotations[[sample_name]][[type_name]]

print(head(mysel_celltrek))


all_cell_cols = intersect(colnames(mysel_celltrek),c("Astrocyte","Bcells","CD4_Tcells","CD8_Tcells","Committed oligodendrocyte precursor","DC","Fibroblast","M1","M2","M3","M4","M5","Microglia","Monocyte","Neuron","NK/NK-like","Oligodendrocyte","Oligodendrocyte precursor","Other_Tcells","TAMs","Vascular"))

tumor_cols = intersect(colnames(mysel_celltrek),c("M1","M2","M3","M4","M5"))


print(head(mysel_celltrek))

print(rowSums(mysel_celltrek[tumor_cols]))

print(rowSums(mysel_celltrek[all_cell_cols]))

mysel_celltrek$tumor_fraction = rowSums(mysel_celltrek[tumor_cols])/rowSums(mysel_celltrek[all_cell_cols])

mysel_celltrek_sub = mysel_celltrek[c('spot_id','tumor_fraction')]
  # all_celltypes_per_spot2 = spots_to_cells_celltrek_annotations[[sample_name]][[type_name]]
  # all_celltypes_per_spot2_melt = melt(all_celltypes_per_spot2,id='spot_id')

  seurat_obj_samp = seurat_obj_list2[[sample_name]]
  seurat_obj_meta = seurat_obj_samp@meta.data[c(integ_clus_vec)]
  seurat_obj_meta$spot_id = rownames(seurat_obj_meta)

all_celltypes_per_spot2_melt_sub = subset(mysel_celltrek_sub,spot_id %in% seurat_obj_meta$spot_id)

  mytumor_frac_by_sample = all_celltypes_per_spot2_melt_sub %>% left_join(seurat_obj_meta, by = c("spot_id"))
  mytumor_frac_by_sample$spot_name = paste0(sample_name,'_',mytumor_frac_by_sample$spot_id)

return(mytumor_frac_by_sample[c('spot_name','integrated_snn_res.0.9','tumor_fraction')])
}

get_enrichment_score_for_eachcelltype  <- function(celltype_proportion_byspotcountcluster,celltype,spot_count_bycluster){

	celltype_proportion_byspotcountcluster_sub = subset(celltype_proportion_byspotcountcluster,variable==celltype)

	total_spot_count = sum(spot_count_bycluster$spot_count_bycluster)

	total_celltype_count = sum(celltype_proportion_byspotcountcluster_sub$celltype_count_by_cluster)

	for(i in 1:nrow(celltype_proportion_byspotcountcluster_sub)){

		cluster_name = celltype_proportion_byspotcountcluster_sub$integrated_snn_res.0.9[i]
		celltype_count_by_cluster = celltype_proportion_byspotcountcluster_sub$celltype_count_by_cluster[i]
		spot_count_bycluster2 = celltype_proportion_byspotcountcluster_sub$spot_count_bycluster[i]
		print(spot_count_bycluster2)
		print(celltype_count_by_cluster)
		print(total_spot_count)
		print(total_celltype_count)
		es = (celltype_count_by_cluster/spot_count_bycluster2)/(total_celltype_count/total_spot_count)
		print(es)
		celltype_proportion_byspotcountcluster_sub$enrichment_score_by_celltype[i] = es
	}

return(celltype_proportion_byspotcountcluster_sub)
}

Z_test_One_Proportion <- function(total_celltype_count_bycluster,integ_res,cluster_var,total_celltype_count,celltype){

	#this code compares for each bar of the stacked barplot the observed vs expected Z test for each cell state
	# against the others. e.g. 

	#  variable integrated_snn_res.0.9 celltype_count_by_cluster
	#   <fct>    <fct>                                      <dbl>
	# 1 M1       0                                            292
	# 2 M2       0                                            729
	# 3 M3       0                                           1372
	# 4 M4       0                                            521
	# 5 M5       0                                           2233
	#for cluster 0

	expected_proportion <- total_celltype_count$prop[total_celltype_count$variable %in% celltype]
	print(expected_proportion)

	total_celltype_count_bycluster_sub = total_celltype_count_bycluster[total_celltype_count_bycluster[[integ_res]]==cluster_var,]

	celltype_count = total_celltype_count_bycluster_sub$total_cell_count[total_celltype_count_bycluster_sub$variable==celltype]

	print(celltype_count)

	total_celltype_count_percluster = sum(total_celltype_count_bycluster_sub$total_cell_count)
	print(total_celltype_count_percluster)

	# Perform one-tailed test for proportion
	prop_test <- prop.test(x = celltype_count, n = total_celltype_count_percluster, 
	                   p = expected_proportion, alternative = "greater")

	# View the result
	mydf = data.frame('celltype'=celltype,'ztestpval'=prop_test$p.value,'cluster'=cluster_var)
	return(mydf)
}



# integ_clus_vec_celltypes_integrated_snn_res.0.9 = get_celltype_by_integclus(spots_to_cells_celltrek_annotations=spots_to_cells_celltrek_annotations,type_name='broad_celltypes_per_spot',integ_clus_vec=integ_clus_vec,seurat_obj_list=seurat_obj_list2)


tumor_fraction_list = list()

for(sample_name in names(seurat_obj_list2)){
tumor_fraction_list[[sample_name]] =  get_tumor_fraction_by_cluster(spots_to_cells_celltrek_annotations=spots_to_cells_celltrek_annotations,type_name='broad_celltypes_per_spot',integ_clus_vec=integ_clus_vec,seurat_obj_list=seurat_obj_list2,sample_name=sample_name)
}

tumor_fraction_df = do.call(rbind, tumor_fraction_list)
rownames(tumor_fraction_df) = NULL

integ_clus_vec_spot_count_integrated_snn_res.0.9 = get_spot_count_bycluster(spots_to_cells_celltrek_annotations=spots_to_cells_celltrek_annotations,type_name='broad_celltypes_per_spot',integ_clus_vec=integ_clus_vec,seurat_obj_list=seurat_obj_list2)

integ_clus_vec_celltypes_integrated_snn_res.0.9 = get_celltype_by_integclus(spots_to_cells_celltrek_annotations=spots_to_cells_celltrek_annotations,type_name='broad_celltypes_per_spot',integ_clus_vec=integ_clus_vec,seurat_obj_list=seurat_obj_list2)

integ_clus_vec_celltypes_integrated_snn_res.0.9_celltypes = integ_clus_vec_celltypes_integrated_snn_res.0.9$integrated_snn_res.0.9

integ_clus_vec_celltypes_integrated_snn_res.0.9_celltypes$integrated_snn_res.0.9 = droplevels(integ_clus_vec_celltypes_integrated_snn_res.0.9_celltypes$integrated_snn_res.0.9)

celltype_proportion_byspotcountcluster = integ_clus_vec_celltypes_integrated_snn_res.0.9_celltypes %>% left_join(integ_clus_vec_spot_count_integrated_snn_res.0.9$integrated_snn_res.0.9, by = c("integrated_snn_res.0.9"))


celltype_score_enrichment_list = list()

for(celltype in levels(celltype_proportion_byspotcountcluster$variable)){
celltype_score_enrichment_list[[celltype]] = get_enrichment_score_for_eachcelltype(celltype_proportion_byspotcountcluster=celltype_proportion_byspotcountcluster,celltype=celltype,spot_count_bycluster=integ_clus_vec_spot_count_integrated_snn_res.0.9$integrated_snn_res.0.9)
}

celltype_score_enrichment_df = do.call(rbind, celltype_score_enrichment_list)

rownames(celltype_score_enrichment_df)  = NULL

total_celltype_count = celltype_proportion_byspotcountcluster %>% group_by(variable) %>% summarise(total_cell_count = sum(celltype_count_by_cluster))

allcellcount = sum(total_celltype_count$total_cell_count)

total_celltype_count_bycluster = celltype_proportion_byspotcountcluster %>% group_by(variable,integrated_snn_res.0.9) %>% summarise(total_cell_count = sum(celltype_count_by_cluster))

total_celltype_count$prop = total_celltype_count$total_cell_count/sum(total_celltype_count$total_cell_count)


cellcolor_names = c("M1","M2","M3","M4","M5","TAMs",
	"Microglia","Monocyte","CD4_Tcells","CD8_Tcells","Other_Tcells","Vascular","Fibroblast","Neuron","Astrocyte","Committed oligodendrocyte precursor","Oligodendrocyte","Oligodendrocyte precursor","DC","Bcells","NK/NK-like")
latest_cell_colors2 = c('#E58606','#6BAED6','#E7BA52','#5254A3','#D66B6F','#70CCB0','#91C1B2','#91C292',
	'#3969AC','#4C8CE6','#82C2E8','#FF4219','#D781B5','#E099CA','gray90','gray50','grey30','gray10',"#ffe119","khaki","khaki2")

names(latest_cell_colors2) = cellcolor_names

mycolor_df = data.frame(celltype=cellcolor_names,cellcolors=latest_cell_colors2)

mycolor_df$cellcolors[mycolor_df$celltype=='Fibroblast'] = 'salmon'
mycolor_df$cellcolors[mycolor_df$celltype=='Neuron'] = 'red4'



mydf_list = list()

# cellstate_vec = paste0('M',1:5)

for(myclus in unique(celltype_proportion_byspotcountcluster$integrated_snn_res.0.9)){
	for(mycelltype in cellcolor_names){
		myother_celltype = setdiff(cellcolor_names,mycelltype)
		mydf_list[[paste0(myclus,'_',mycelltype)]] = Z_test_One_Proportion(total_celltype_count_bycluster=total_celltype_count_bycluster,integ_res='integrated_snn_res.0.9',cluster_var=myclus,total_celltype_count=total_celltype_count,celltype=mycelltype)
	}
}

mydf_list_comb = do.call(rbind,mydf_list)
rownames(mydf_list_comb) = NULL

mydf_list_comb_spl = split(mydf_list_comb,mydf_list_comb$celltype)


mydf_list_comb_spl2 = lapply(mydf_list_comb_spl, function(x) {
	mypvals = data.frame(p.adj = p.adjust(x$ztestpval,method = 'bonferroni',n=length(x$ztestpval)))
	mydf_mod = cbind(x,mypvals)
	return(mydf_mod)
})


mydf_list_comb_spl_comb = do.call(rbind,mydf_list_comb_spl2)

rownames(mydf_list_comb_spl_comb) = NULL

mydf_list_comb_spl_comb2 = mydf_list_comb_spl_comb[c('celltype','cluster','p.adj')]

colnames(mydf_list_comb_spl_comb2) = c('variable','integrated_snn_res.0.9','p.adj')

celltype_score_enrichment_df2 = merge(celltype_score_enrichment_df, mydf_list_comb_spl_comb2, by=c('variable','integrated_snn_res.0.9'))

celltype_score_enrichment_df2$Significance = 'Not Significant'

celltype_score_enrichment_df2$Significance[celltype_score_enrichment_df2$p.adj<=0.05] = 'Significant'



group_color_palette = c("#d7191c","#fdae61","#2b83ba","#abdda4")
names(group_color_palette) = c('core-rich','transition','edge-rich','non-specific')

# dotplot_items$group_color_palette
# celltype_score_enrichment_df2_grouped = dotplot_items$celltype_score_enrichment_df2_grouped
# mycolor_df = dotplot_items$mycolor_df


mycolor_df$celltype[mycolor_df$celltype=="Oligodendrocyte precursor"] = "OPC"
mycolor_df$celltype[mycolor_df$celltype=="Committed oligodendrocyte precursor"] = "C-OPC"

mycolor_df$celltype[mycolor_df$celltype=="Oligodendrocyte"] = "Oligo"

group_df = data.frame(cluster=paste0('C',c(0,1,2,13,23,3,6,8,9,12,16,17,18,20,22,19,21,4,7,5,10,15)),group=c(rep('edge-rich',5),rep('core-rich',10),rep('transition',4),rep('non-specific',3)))

colnames(group_df) = c('Cluster','group')

group_df$group = factor(group_df$group, levels=c('core-rich','transition','edge-rich','non-specific'))

group_color_palette = c('#d7191c','#fdae61','#2b83ba','#abdda4')
names(group_color_palette) = levels(group_df$group)

celltype_score_enrichment_df2$Cluster = paste0('C',celltype_score_enrichment_df2$integrated_snn_res.0.9)

celltype_score_enrichment_df2_grouped = celltype_score_enrichment_df2 %>% left_join(group_df, by = c("Cluster"))

tumor_fraction_df$Cluster = paste0('C',tumor_fraction_df$integrated_snn_res.0.9)

tumor_fraction_groupeddf = tumor_fraction_df %>%  group_by(Cluster) %>% summarize(avgtumor_fraction=mean(tumor_fraction))

#Latest naming convention celltypes/cellstates
#change factor levels according to tumor fraction in each cluster higher to lower
# celltype_score_enrichment_df2$celltype <- factor(celltype_score_enrichment_df2$celltype, levels = c("M1","M2","M3","M4","M5","CD4_Tcells","CD8_Tcells","Other_Tcells","Bcells","NK/NK-like","DC","Microglia","Monocyte","TAMs","Astrocyte","C-OPC","Oligo","OPC","Neuron","Vascular","Fibroblast"))
celltype_score_enrichment_df2_grouped$celltype = as.character(celltype_score_enrichment_df2_grouped$variable)
celltype_score_enrichment_df2_grouped$celltype[celltype_score_enrichment_df2_grouped$variable=='M1'] = 'S1'
celltype_score_enrichment_df2_grouped$celltype[celltype_score_enrichment_df2_grouped$variable=='M2'] = 'S2'
celltype_score_enrichment_df2_grouped$celltype[celltype_score_enrichment_df2_grouped$variable=='M3'] = 'S3'
celltype_score_enrichment_df2_grouped$celltype[celltype_score_enrichment_df2_grouped$variable=='M4'] = 'S4'
celltype_score_enrichment_df2_grouped$celltype[celltype_score_enrichment_df2_grouped$variable=='M5'] = 'S5'


celltype_score_enrichment_df2_grouped$celltype[celltype_score_enrichment_df2_grouped$variable=="Oligodendrocyte precursor"] = "OPC"
celltype_score_enrichment_df2_grouped$celltype[celltype_score_enrichment_df2_grouped$variable=="Committed oligodendrocyte precursor"] = "C-OPC"
celltype_score_enrichment_df2_grouped$celltype[celltype_score_enrichment_df2_grouped$variable=="Oligodendrocyte"] = "Oligo"


mycolor_df$celltype[mycolor_df$celltype=='M1'] = 'S1'
mycolor_df$celltype[mycolor_df$celltype=='M2'] = 'S2'
mycolor_df$celltype[mycolor_df$celltype=='M3'] = 'S3'
mycolor_df$celltype[mycolor_df$celltype=='M4'] = 'S4'
mycolor_df$celltype[mycolor_df$celltype=='M5'] = 'S5'

celltype_score_enrichment_df2_grouped$Cluster  = factor(celltype_score_enrichment_df2_grouped$Cluster,levels=c('C3','C20','C8','C17','C22','C16','C9','C18','C12','C6','C7','C21','C19','C4','C1','C13','C23','C0','C2','C5','C10','C15'))

celltype_score_enrichment_df2_grouped$celltype = factor(celltype_score_enrichment_df2_grouped$celltype,levels= c("S1","S2","S3","S4","S5","CD4_Tcells","CD8_Tcells","Other_Tcells","Bcells","NK/NK-like","DC","Microglia","Monocyte","TAMs","Astrocyte","C-OPC","Oligo","OPC","Neuron","Vascular","Fibroblast"))

#Pathologist annotations of core/edge/transition by Spot
mycluster_df = readRDS('mycluster_df.rds')
colnames(mycluster_df) = c('Cluster','region')
mycluster_df$Cluster = paste0('C',mycluster_df$Cluster)

#Pathologist annotations colors
group_color_palette2 = c('#d7191c','#fdae61','#2b83ba')
names(group_color_palette2) = c('Core','Transition','Edge')

group_color_palette2['Transition'] = '#4A0080'

mypercentdf_patient = as.data.frame(prop.table(table(mycluster_df$Cluster,mycluster_df$region),1)*100)
colnames(mypercentdf_patient) <- c("Cluster","region","Percentage")
mypercentdf_patient$Cluster =  factor(mypercentdf_patient$Cluster,levels=c('C3','C20','C8','C17','C22','C16','C9','C18','C12','C6','C7','C21','C19','C4','C1','C13','C23','C0','C2','C5','C10','C15'))



celltype_score_enrichment_df2_grouped <- celltype_score_enrichment_df2_grouped  %>% mutate(
    Cluster = factor(Cluster, levels = tumor_fraction_groupeddf$Cluster)  # Use the same levels as the dotplot
  )

tumor_fraction_grouped$Cluster  = factor(tumor_fraction_grouped$Cluster,levels=c('C3','C20','C8','C17','C22','C16','C9','C18','C12','C6','C7','C21','C19','C4','C1','C13','C23','C0','C2','C5','C10','C15'))
# tumor_fraction_grouped$group.x = NULL
# tumor_fraction_grouped$group = tumor_fraction_grouped$group.y
# tumor_fraction_grouped$group.y  = NULL

# dotplot_items_list = list(
# mypercentdf_patient=mypercentdf_patient,
# group_color_palette2=group_color_palette2,
# tumor_fraction_grouped=tumor_fraction_grouped,
# celltype_score_enrichment_df2_grouped=celltype_score_enrichment_df2_grouped,
# group_color_palette=group_color_palette,
# mycolor_df=mycolor_df,
# group_df=group_df)

# saveRDS(dotplot_items_list,'dotplot_items_list_march13.2025.rds')

#plotting requires R version 4.3.3 (2023-10-31)
#and the following packages
# library(ggdendro)
# library(cowplot)
# library(tidyverse)
# library(ggtree) # install with `devtools::install_github("YuLab-SMU/ggtree")` as you need a version newer than what bioconductor serves
# library(patchwork) 
# library(gridExtra)
# library(gtools)
# library(ggnewscale)
# # Ensure proper ordering of Cluster and celltype
# library(patchwork)


dotplot_items2 = readRDS('dotplot_items_list_march13.2025.rds')
names(dotplot_items2)
mypercentdf_patient = dotplot_items2$mypercentdf_patient
tumor_fraction_grouped = dotplot_items2$tumor_fraction_grouped
group_color_palette = dotplot_items2$group_color_palette
group_color_palette2 = dotplot_items2$group_color_palette2
celltype_score_enrichment_df2_grouped = dotplot_items2$celltype_score_enrichment_df2_grouped
mycolor_df = dotplot_items2$mycolor_df
group_df = dotplot_items2$group_df


text_size <- 20  # Set text size
legend_text_size <- 16  # Set legend text size for all legends
legend_title_size <- 18  # Set legend title size for all legends

barplot <- ggplot(mypercentdf_patient, aes(fill = region, y = Percentage, x = Cluster)) +
  geom_bar(
    position = "fill", stat = "identity",
    alpha = 0.9, color = "black", linewidth = 0.05
  ) +
  scale_fill_manual(values = group_color_palette2) +
  labs(x = 'Cluster', y = "Pathologist annotations") +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.y = element_text(size = text_size),  # Uniform y-axis text size
    axis.title.y = element_text(size = text_size),  # Uniform y-axis title size
    legend.text = element_text(size = legend_text_size),  # Uniform legend text size
    legend.title = element_text(size = legend_title_size),  # Uniform legend title size
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(fill = guide_legend(title = "Pathologist annotations"))


# cluster_names = paste0('C',names(final_cluster_colors))
# names(final_cluster_colors) = cluster_names


library(ggnewscale)
# Ensure proper ordering of Cluster and celltype
library(patchwork)

boxplot <- ggplot(tumor_fraction_grouped, aes(x = Cluster, y = tumor_fraction, fill=group)) +
  geom_boxplot(alpha=0.6) + scale_fill_manual(values = group_color_palette,  # Define a palette for group
                     name = "Clusters") + theme_minimal() + theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = text_size),  # Uniform text size for x-axis
    axis.title.x = element_text(size = text_size),  # Uniform x-axis title size
    axis.title.y = element_text(size = text_size),  # Uniform y-axis title size
    axis.text.y = element_text(size = text_size),  # Uniform y-axis text size
    legend.position = "none"  # Hide legends (shared legends will appear at the top)
  ) + labs(x = "Cluster", y = "Tumor Fraction")

# Create the dotplot without x-axis labels
dotplot <- celltype_score_enrichment_df2_grouped %>%
  ggplot(aes(x = Cluster, 
             y = celltype, 
             fill = celltype,  # Using fill for coloring points
             size = enrichment_score_by_celltype)) + 
  
  # Dotplot layer
  geom_point(aes(stroke = ifelse(Significance == "Significant", 3, 0), 
                 color = Significance),  # Keep color mapped to Significance to show legend
             shape = 21) +  # Shape that supports fill and stroke
  
  # Set minimum and maximum point size
  scale_size_continuous(range = c(2, 22), 
                        name = "Enrichment Score") + 
  
  # Manually set the colors for the significance legend
  scale_color_manual(values = c("Significant" = "black", "Non-significant" = NA),
                     name = "Significance", 
                     labels = c("Significant", "Non-significant"),
                     guide = guide_legend(
                       title = "Significance", 
                       override.aes = list(size = 8, stroke = 2),  # Larger legend points
                       title.theme = element_text(size = legend_title_size,face="bold"), 
                       label.theme = element_text(size = legend_text_size,face="bold")
                     )) + 
  
  # Vertical dotted lines separating clusters
  geom_vline(xintercept = seq(1.5, length(unique(celltype_score_enrichment_df2_grouped$Cluster)) - 0.5, by = 1), 
             linetype = "dotted", 
             color = "black", 
             linewidth = 0.5) + 
  
  # Custom theme with uniform text sizes
  cowplot::theme_cowplot() + 
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.y = element_text(size = text_size,face="bold"),  # Uniform y-axis text size
    axis.title.y = element_text(size = text_size,face="bold"),  # Uniform y-axis title size
    legend.text = element_text(size = legend_text_size),  # Uniform legend text size
    legend.title = element_text(size = legend_title_size),  # Uniform legend title size
    legend.position = "top"
  ) +  
  
  # Custom color scale for celltype fills
  scale_fill_manual(values = setNames(mycolor_df$cellcolors, mycolor_df$celltype),
                    guide = guide_legend(
                      title = "Celltypes", override.aes = list(size = 8, stroke = 1.5),  # Larger legend points
                      title.theme = element_text(size = legend_title_size,face="bold"), 
                      label.theme = element_text(size = legend_text_size,face="bold")
                    )) +
  
  # Add a new scale for ClusterGroup
  ggnewscale::new_scale_fill() +  # Reset fill aesthetic for the next geom
  
  # Add top bar for cluster groupings
  geom_tile(data = celltype_score_enrichment_df2_grouped %>%
              distinct(Cluster, group) %>%
              as.data.frame(),  # Ensure it's a valid data frame
            mapping = aes(x = Cluster, y = -0.5, fill = group), 
            inherit.aes = FALSE, 
            height = 1) +  # Adjust height to make it visually distinct
  
  # Customize bar colors for ClusterGroup
  scale_fill_manual(values = group_color_palette,  # Define a palette for ClusterGroup
                    name = "Cluster Group", 
                    guide = guide_legend(
                      title = "Cluster Group",
                      override.aes = list(size = 8),  # Larger legend points
                      title.theme = element_text(size = legend_title_size), 
                      label.theme = element_text(size = legend_text_size)
                    )) + 
  
  # Horizontal dotted lines separating categories
  geom_hline(yintercept = seq(1.5, length(unique(celltype_score_enrichment_df2_grouped$celltype)) - 0.5, by = 1), 
             linetype = "dotted", 
             color = "black", 
             linewidth = 0.5) +  

  labs(y = "Celltype Enrichment Score")

dotplot = dotplot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
barplot = barplot+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
# barplot = integ_clus_vec_celltypes_integrated_snn_res.0.9_celltypes_barplot_bypatient
# Combine the plots with patchwork
combined_plot <- (dotplot / barplot) / boxplot +  plot_layout(ncol = 1, heights = c(4,1, 0.5))  # Boxplot takes less space than the dotplot

# Display the combined plot
print(combined_plot)


# Save the combined plot
ggsave(paste0('combined_plot_with_uniform_legend_sizes.sortedbyavgtumorfractionwithingroups.rowscolsclustered.', date, '.jpg'), 
       combined_plot, 
       width = 30, height = 22.5, units = "in", dpi = 300)
