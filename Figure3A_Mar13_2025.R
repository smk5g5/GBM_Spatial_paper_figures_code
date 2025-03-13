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

merged.patient.combined.integ_sub = readRDS('/n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/cluster_distances/merged.patient.combined.integ_sub.rds')

spots_to_cells_celltrek_annotations = readRDS('/n/scratch/users/s/sak4832/Feb2_2025/Dec27_2024/cluster_distances/spots_to_cells_celltrek_annotations.rds')

single_cell_ref = readRDS('/n/scratch/users/s/sak4832/scrna_seq_ref_celltrek_orig_latestcellstateannotations.rds')

make_broad_celltypes2 <- function(single_cell_ref){
  mydf = data.frame(names(table(single_cell_ref$Celltypes_Aug8_2025)))
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
  mydf2 = data.frame(cellnames=names(single_cell_ref$Celltypes_Aug8_2025),Celltypes_Aug8_2025=unname(single_cell_ref$Celltypes_Aug8_2025))
  
  mydf12 = mydf2 %>% left_join(mydf, by = c("Celltypes_Aug8_2025"="all_celltypes"))

  index = match(Cells(single_cell_ref),mydf12$cellnames)
  single_cell_ref$broad_celltypes = mydf12$broad_celltypes[index]
  return(single_cell_ref)
}