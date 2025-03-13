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

set.seed(12345L) 

date = gsub("2025-","25",Sys.Date(),perl=TRUE);
date = gsub("-","",date);


 = readRDS('/n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/cluster_distances/merged.patient.combined.integ_sub.rds')

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






