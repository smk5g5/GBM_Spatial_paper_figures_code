library(SingleCellExperiment)
library(scuttle)
library(scran)
library(scater)
library(scRNAseq)
library(Seurat)
# library(BiocParallel)
library(yaml)
library(Ckmeans.1d.dp)
library(stringr)
library(dplyr)
library(sctransform)
library(ggplot2)
library(sctransform)
library(future)
library(RColorBrewer)
# library(circlize)
# library(scCustomize)

# library(circlize)
# library(ComplexHeatmap)
# library(stringr)
library(reshape2)

# library(clusterProfiler)
# library(org.Hs.eg.db)
library("CellTrek")

library(tidyr)
library(dplyr)

set.seed(12345L) 


date = gsub("2025-","25",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

scrna_seq_ref_celltrek = readRDS('/n/scratch/users/s/sak4832/scrna_seq_ref_celltrek_orig_latestcellstateannotations.rds')

Idents(scrna_seq_ref_celltrek) = 'Celltypes_Aug8_2024'
scrna_seq_ref_celltrek_tumor_only = subset(scrna_seq_ref_celltrek,idents=c(paste0('M',1:5)))


neftel_gbm = read.table('/n/data1/mgh/neuro/petti/lab/Users/khan.saad/key.gene.lists/neftel_gbm_signatures.txt',header=T,sep="\t")
neftel_gbm_list = lapply(as.list(neftel_gbm),function(x) {grep("^$",x,invert = T,value = T)})

final_genemodules_reduced_mutualinfo = readRDS('final_genemodules_reduced_mutualinfo_list.rds')


assign_tumor_states <- function(scRNAseq_ref,gene_list,state_pub){

	for(i in names(gene_list)){
	scRNAseq_ref <- AddModuleScore(
	object = scRNAseq_ref,
	features = list(gene_list[[i]]),
	name = i
	)
	}

module_score_names = paste0(names(gene_list),1)
names(module_score_names) = names(gene_list)
sel_df = scRNAseq_ref@meta.data[module_score_names]

mydf = data.frame(cellnames=character(nrow(sel_df)),tumor_state=character(nrow(sel_df)), stringsAsFactors = FALSE)

for(i in 1:nrow(sel_df)){
  cellname = rownames(sel_df[i,])
  max_module = names(module_score_names[module_score_names==names(which.max(sel_df[i,]))])
  mydf$cellnames[i] = cellname
  mydf$tumor_state[i] = max_module
}

tumor_state_celllist =  lapply(df_list, function(x) {x$cellnames})
index <- match(Cells(samp_celltrek),mydf$cellnames)

if(state_pub=='neftel'){
	scRNAseq_ref$neftel_states <- mydf$tumor_state[index]
	}else{
		scRNAseq_ref$tumor_state <- mydf$tumor_state[index]
	}
return(scRNAseq_ref)
}


scrna_seq_ref_celltrek_tumor_only_nef = assign_tumor_states(scRNAseq_ref=scrna_seq_ref_celltrek_tumor_only,gene_list=neftel_gbm_list,state_pub='neftel')

make_heatmap_by_tumor_states <- function(scRNAseq_ref,tumor_ident){

}



