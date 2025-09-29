##########
#Figure 5A
##########

stepped.final.colors <- c(#Tumor
  "#6BAED6","#A61275","#E7BA52","#5254A3","#D66B6F","#FF9408",
  #Myeloid
  "#70CC71","#70CCB0","#70CC90","#95CC70","#70CACC","#B4CCBD","#91C292","#91C1B2",
  #T cells
  "#3969AC","#4C8CE6","#38ABA5","#82C2E8","#2A4E80",
  #Endothelial
  "#B6F52C","#2CF545","#62F52C","#DEF7A8","#2CF589","#9DF57E","#F5EF2C",
  #Vascular
  "#FF19A3","#FF4219","#E31717","#E619FF","#E619FF","#FF6F6F","#D781B5",
  #Neuronal
  "#A13FE0","#4516DC","#8341E8","#BD98E0","#99A2E0","#E099CA","#BAB1E0","#D9B1E0","#B1B8E0","#313861")


outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/General_Figures"

Idents(Aggregated_GBM) <- "Celltype_X"
subsample <- Aggregated_GBM[, sample(colnames(Aggregated_GBM), size =100000, replace=F)]
Idents(subsample) <- "Celltype_X"

pdf(sprintf("%s/UMAP.subsample.Cohort.2.Celltype.Meta.X.pdf", outdir), width = 20, height = 12);
scatter <- DimPlot(subsample, label=FALSE, cols=(stepped.final.colors), raster=FALSE)
print(scatter);
dev.off();

##############
#Figure 5B, 5E
##############

#For overall sample ImageDimPlots w/ Celltype colored

outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/UMAP"

sample <- unique(Aggregated_GBM@meta.data$orig.ident)

Idents(Aggregated_GBM) <- "orig.ident"

objects <- list()
for (i in 1:length(unique(Aggregated_GBM@meta.data$orig.ident))){
  objects[[i]] <- subset(Aggregated_GBM, idents=sample[i])
}

colors <- as.data.frame(stepped.final.colors)
colors$celltype <- levels(Aggregated_GBM@meta.data$Celltype_X)

for (i in 1:length(sample)){
  table <- table(objects[[i]]@meta.data$Celltype_X)
  table <- as.data.frame(table)
  present_cells <- table[table$Freq > 0, ]
  present_cells$celltype <- present_cells$Var1
  
  subset_colors <- merge(present_cells, colors, by = "celltype", sort = FALSE)
  
  subset_colors <- subset_colors[match(present_cells$celltype, subset_colors$celltype), ]
  subset_colors <- as.vector(subset_colors$stepped.final.colors)
  
  Idents(objects[[i]]) <- "Celltype_X"
  
  jpeg(sprintf("%s/ImageDimPlot.%s.jpg", outdir,sample[i]), width = 20, height = 12, units="in", res=300);
  scatter <- ImageDimPlot(objects[[i]],cols=(subset_colors))
  print(scatter);
  dev.off();
  
  pdf(sprintf("%s/ImageDimPlot.%s.pdf", outdir,sample[i]), width = 20, height = 12);
  scatter <- ImageDimPlot(objects[[i]],cols=(subset_colors))
  print(scatter);
  dev.off();
  
}

#For overall sample ImageDimPlots w/ Banksy niche colored

Idents(Aggregated_GBM) <- "BANKSY_X"

niche.colors <- c("#E7BA52","#A61275","#5254A3","#6BAED6",
                  "#673147","#FF8C00","#2CF545","#A13FE0")

outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/UMAP"

sample <- unique(Aggregated_GBM@meta.data$orig.ident)

Idents(Aggregated_GBM) <- "orig.ident"

objects <- list()
for (i in 1:length(unique(Aggregated_GBM@meta.data$orig.ident))){
  objects[[i]] <- subset(Aggregated_GBM, idents=sample[i])
}

colors <- as.data.frame(niche.colors)
colors$celltype <- levels(Aggregated_GBM@meta.data$BANKSY_X)

for (i in 1:length(sample)){
  table <- table(objects[[i]]@meta.data$BANKSY_X)
  table <- as.data.frame(table)
  present_cells <- table[table$Freq > 0, ]
  present_cells$celltype <- present_cells$Var1
  
  subset_colors <- merge(present_cells, colors, by = "celltype", sort = FALSE)
  
  subset_colors <- subset_colors[match(present_cells$celltype, subset_colors$celltype), ]
  subset_colors <- as.vector(subset_colors$niche.colors)
  
  Idents(objects[[i]]) <- "BANKSY_X"
  
  jpeg(sprintf("%s/ImageDimPlot.BANKSY.%s.jpg", outdir,sample[i]), width = 20, height = 12, units="in", res=300);
  scatter <- ImageDimPlot(objects[[i]],cols=(subset_colors))
  print(scatter);
  dev.off();
  
  pdf(sprintf("%s/ImageDimPlot.BANKSY.%s.pdf", outdir,sample[i]), width = 20, height = 12);
  scatter <- ImageDimPlot(objects[[i]],cols=(subset_colors))
  print(scatter);
  dev.off();
  
}

#For sample specific FOVs w/ celltype colored

stepped.final.colors <- c(#Tumor
  "#6BAED6","#A61275","#E7BA52","#5254A3","#D66B6F","#FF9408",
  #Myeloid
  "#70CC71","#70CCB0","#70CC90","#95CC70","#70CACC","#B4CCBD","#91C292","#91C1B2",
  #T cells
  "#3969AC","#4C8CE6","#38ABA5","#82C2E8","#2A4E80",
  #Endothelial
  "#B6F52C","#2CF545","#62F52C","#DEF7A8","#2CF589","#9DF57E","#F5EF2C",
  #Vascular
  "#FF19A3","#FF4219","#E31717","#E619FF","#E619FF","#FF6F6F","#D781B5",
  #Neuronal
  "#A13FE0","#4516DC","#8341E8","#BD98E0","#99A2E0","#E099CA","#BAB1E0","#D9B1E0","#B1B8E0","#313861")

outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/UMAP/ImageDimPlot"

sample <- unique(Aggregated_GBM@meta.data$orig.ident)

Idents(Aggregated_GBM) <- "orig.ident"

objects <- list()
for (i in 4:4){
  objects[[i]] <- subset(Aggregated_GBM, idents=sample[i])
}

colors <- as.data.frame(stepped.final.colors)
colors$celltype <- levels(Aggregated_GBM@meta.data$Celltype_X)

for (i in 4:4){
  table <- table(objects[[i]]@meta.data$Celltype_X)
  table <- as.data.frame(table)
  present_cells <- table[table$Freq > 0, ]
  present_cells$celltype <- present_cells$Var1
  
  subset_colors <- merge(present_cells, colors, by = "celltype", sort = FALSE)
  
  subset_colors <- subset_colors[match(present_cells$celltype, subset_colors$celltype), ]
  subset_colors <- as.vector(subset_colors$stepped.final.colors)
  
  Idents(objects[[i]]) <- "Celltype_X"
  
  cropped.coords <- Crop(objects[[i]][["fov.4"]], x = c(3200, 4200), y = c(2800, 4800), coords = "plot")
  objects[[i]][["zoom"]] <- cropped.coords
  
  jpeg(sprintf("%s/ImageDimPlot.%s.subset.jpg", outdir,sample[i]), width = 20, height = 12, units="in", res=300);
  scatter <- ImageDimPlot(objects[[i]],cols=(subset_colors),fov="zoom", size=1)
  print(scatter);
  dev.off();
  
  pdf(sprintf("%s/ImageDimPlot.%s.subset.pdf", outdir,sample[i]), width = 20, height = 12);
  scatter <- ImageDimPlot(objects[[i]],cols=(subset_colors),fov="zoom", size=1)
  print(scatter);
  dev.off();
  
  cropped.coords <- Crop(objects[[i]][["fov.4"]], x = c(500, 1500), y = c(1800, 3800), coords = "plot")
  objects[[i]][["zoom"]] <- cropped.coords
  
  jpeg(sprintf("%s/ImageDimPlot.%s.subset.2.jpg", outdir,sample[i]), width = 20, height = 12, units="in", res=300);
  scatter <- ImageDimPlot(objects[[i]],cols=(subset_colors),fov="zoom", size=1)
  print(scatter);
  dev.off();
  
  pdf(sprintf("%s/ImageDimPlot.%s.subset.2.pdf", outdir,sample[i]), width = 20, height = 12);
  scatter <- ImageDimPlot(objects[[i]],cols=(subset_colors),fov="zoom", size=1)
  print(scatter);
  dev.off();
  
}


objects <- list()
for (i in 2:2){
  objects[[i]] <- subset(Aggregated_GBM, idents=sample[i])
}

colors <- as.data.frame(stepped.final.colors)
colors$celltype <- levels(Aggregated_GBM@meta.data$Celltype_X)

for (i in 2:2){
  table <- table(objects[[i]]@meta.data$Celltype_X)
  table <- as.data.frame(table)
  present_cells <- table[table$Freq > 0, ]
  present_cells$celltype <- present_cells$Var1
  
  subset_colors <- merge(present_cells, colors, by = "celltype", sort = FALSE)
  
  subset_colors <- subset_colors[match(present_cells$celltype, subset_colors$celltype), ]
  subset_colors <- as.vector(subset_colors$stepped.final.colors)
  
  Idents(objects[[i]]) <- "Celltype_X"
  
  cropped.coords <- Crop(objects[[i]][["fov.2"]], x = c(2500, 4500), y = c(4500, 6500), coords = "plot")
  objects[[i]][["zoom"]] <- cropped.coords
  
  jpeg(sprintf("%s/ImageDimPlot.%s.subset.jpg", outdir,sample[i]), width = 20, height = 12, units="in", res=300);
  scatter <- ImageDimPlot(objects[[i]],cols=(subset_colors),fov="zoom", size=1)
  print(scatter);
  dev.off();
  
  pdf(sprintf("%s/ImageDimPlot.%s.subset.pdf", outdir,sample[i]), width = 20, height = 12);
  scatter <- ImageDimPlot(objects[[i]],cols=(subset_colors),fov="zoom", size=1)
  print(scatter);
  dev.off();
  
}

#For sample specific FOVs w/ BANKSY niche colored

Idents(Aggregated_GBM) <- "orig.ident"

niche.colors <- c("#E7BA52","#A61275","#5254A3","#6BAED6",
                  "#673147","#FF8C00","#2CF545","#A13FE0")

objects <- list()
for (i in 2:2){
  objects[[i]] <- subset(Aggregated_GBM, idents=sample[i])
}

colors <- as.data.frame(niche.colors)
colors$celltype <- levels(Aggregated_GBM@meta.data$BANKSY_X)

for (i in 2:2){
  table <- table(objects[[i]]@meta.data$BANKSY_X)
  table <- as.data.frame(table)
  present_cells <- table[table$Freq > 0, ]
  present_cells$celltype <- present_cells$Var1
  
  subset_colors <- merge(present_cells, colors, by = "celltype", sort = FALSE)
  
  subset_colors <- subset_colors[match(present_cells$celltype, subset_colors$celltype), ]
  subset_colors <- as.vector(subset_colors$niche.colors)
  
  Idents(objects[[i]]) <- "BANKSY_X"
  
  cropped.coords <- Crop(objects[[i]][["fov.2"]], x = c(1500, 6500), y = c(3000, 6000), coords = "plot")
  objects[[i]][["zoom"]] <- cropped.coords
  
  jpeg(sprintf("%s/ImageDimPlot.BANKSY.%s.subset.jpg", outdir,sample[i]), width = 20, height = 12, units="in", res=300);
  scatter <- ImageDimPlot(objects[[i]],cols=(subset_colors),fov="zoom",size=1)
  print(scatter);
  dev.off();
  
  pdf(sprintf("%s/ImageDimPlot.BANKSY.%s.subset.pdf", outdir,sample[i]), width = 20, height = 12);
  scatter <- ImageDimPlot(objects[[i]],cols=(subset_colors),fov="zoom",size=1)
  print(scatter);
  dev.off();
  
}

Idents(Aggregated_GBM) <- "orig.ident"

objects <- list()
for (i in 4:4){
  objects[[i]] <- subset(Aggregated_GBM, idents=sample[i])
}

colors <- as.data.frame(niche.colors)
colors$celltype <- levels(Aggregated_GBM@meta.data$BANKSY_X)

for (i in 4:4){
  table <- table(objects[[i]]@meta.data$BANKSY_X)
  table <- as.data.frame(table)
  present_cells <- table[table$Freq > 0, ]
  present_cells$celltype <- present_cells$Var1
  
  subset_colors <- merge(present_cells, colors, by = "celltype", sort = FALSE)
  
  subset_colors <- subset_colors[match(present_cells$celltype, subset_colors$celltype), ]
  subset_colors <- as.vector(subset_colors$niche.colors)
  
  Idents(objects[[i]]) <- "BANKSY_X"
  
  cropped.coords <- Crop(objects[[i]][["fov.4"]], x = c(500, 5500), y = c(1800, 4800), coords = "plot")
  objects[[i]][["zoom"]] <- cropped.coords
  
  jpeg(sprintf("%s/ImageDimPlot.BANKSY.%s.subset.jpg", outdir,sample[i]), width = 20, height = 12, units="in", res=300);
  scatter <- ImageDimPlot(objects[[i]],cols=(subset_colors),fov="zoom",size=1)
  print(scatter);
  dev.off();
  
  pdf(sprintf("%s/ImageDimPlot.BANKSY.%s.subset.pdf", outdir,sample[i]), width = 20, height = 12);
  scatter <- ImageDimPlot(objects[[i]],cols=(subset_colors),fov="zoom",size=1)
  print(scatter);
  dev.off();
  
}

##########
#Figure 5C
##########

outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/Tumor"
empty_df <- read.csv(sprintf("%s/Avg_Tumor_New_Module_Scores.csv",outdir), row.names=1)

matrix <- as.matrix(empty_df)

outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/General_Figures"

pdf(sprintf("%s/Heatmap.New_Module.Score.by.Sample.pdf", outdir), width = 24, height = 12);
scatter <- pheatmap(matrix, display_numbers=TRUE, color=colorRampPalette(c("blue", "white", "red"))(20), breaks=seq(-1,1,by=0.1))
print(scatter);
dev.off();

##########
#Figure 5D
##########

table <- list()

outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/GBM025A/New"
table[[1]] <- read.csv(sprintf("%s/GBM025A_BANKSY_meta_table_percent.csv",outdir), row.names=1)

outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/GBM026I/New"
table[[2]] <- read.csv(sprintf("%s/GBM026I_BANKSY_meta_table_percent.csv",outdir), row.names=1)

outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/GBM026M/New"
table[[3]] <- read.csv(sprintf("%s/GBM026M_BANKSY_meta_table_percent.csv",outdir), row.names=1)

outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/GBM030/New"
table[[4]] <- read.csv(sprintf("%s/GBM030_BANKSY_meta_table_percent.csv",outdir), row.names=1)

outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/GBM034/New"
table[[5]] <- read.csv(sprintf("%s/GBM034_BANKSY_meta_table_percent.csv",outdir), row.names=1)

tumor.names <- c("GBM024D","GBM024E","GBM024I","GBM024H")
for (i in 1:4){
  outdir <- sprintf("/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Both_Cohort/Cohort_1/BANKSY/%s/New",tumor.names[i])
  table[[i+5]] <- read.csv(sprintf("%s/%s_BANKSY_meta_table_percent.csv",outdir,tumor.names[i]),row.names=1)
}

df <- list()

for (i in 1:length(table)){
  table[[i]]$row <- rownames(table[[i]])
}

final <- table[[1]]

for (i in 2:length(table)){
  final <- merge(final,table[[i]],by='row',all=T) 
}

rownames(final) <- final$row
final <- final[,-1]
final[is.na(final)] <- 0

scale <- scale(final)

#Figures

paletteLength <- 50
myColor <- colorRampPalette(c("blue","white","red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(scale), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scale)/paletteLength, max(scale), length.out=floor(paletteLength/2)))

pdf(sprintf("%s/Heatmap_of_Banksy_clusters.meta.cleaned.pdf", outdir), width = 16, height = 12); 
plot_bar <- pheatmap(scale, color=myColor, breaks=myBreaks,cutree_cols = 6)
print(plot_bar);
dev.off();

##########
#Figure 5F
##########
#Note Mes. was changed to immune infiltrated (BANKSY niches)

outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/BANKSY/Enrichment_Pearson"

#Base file

cohort_2 <- table(Aggregated_GBM@meta.data$Celltype_X, Aggregated_GBM@meta.data$BANKSY_X)
cohort_1 <- table(GBM024.final@meta.data$Celltype_X, GBM024.final@meta.data$BANKSY_X)

cohort_2 <- as.data.frame.matrix(cohort_2)
cohort_1 <- as.data.frame.matrix(cohort_1)

cohort_1_aligned <- cohort_2
cohort_1_aligned[] <- NA  # Set all to NA initially
cohort_1_aligned[rownames(cohort_1), colnames(cohort_1)] <- cohort_1

result <- cohort_2 + replace(cohort_1_aligned, is.na(cohort_1_aligned), 0)

write.csv(result, sprintf("%s/BANKSY.Celltype.Distribution.Total.X.csv",outdir))
result <- read.csv(sprintf("%s/BANKSY.Celltype.Distribution.Total.X.csv",outdir))
result <- result[,-1]

#Hypergeometric test across all 

hypergeometric_test <- function(N, K, n, k) {
  p_value <- 1 - phyper(k - 1, K, N - K, n)
  return(p_value)
}

p_values <- (
  matrix(ncol = ncol(result), nrow = nrow(result))  # Initialize with a matrix of NAs
)
colnames(p_values) <- colnames(result)
rownames(p_values) <- rownames(result)

for (i in 1:length(rownames(result))){
  for (j in 1:length(colnames(result))){
    #Total cells
    N=sum(result)
    #Total number of target cell
    K=sum(result[i,])
    #Number of cells in niche
    n=sum(result[,j])
    #Number of target cells in niche
    k=result[i,j]
    
    p_val <- hypergeometric_test(N,K,n,k)
    p_values[i,j] <- p_val
  }
}

p_values <- as.data.frame.matrix(p_values)
write.csv(p_values,sprintf("%s/Enrichment.Hypergeometric.csv",outdir))

#Figure for OR

outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/BANKSY/Enrichment_Pearson"

result <- read.csv(sprintf("%s/BANKSY.Celltype.Distribution.Total.X.csv",outdir),row.names = 1)
p_values <- read.csv(sprintf("%s/Enrichment.Hypergeometric.csv",outdir),row.names=1)

colnames(result) <- c("Hypoxic.1","Mes.","Mixed.1","Proneural","Hypoxic.2","Mixed.2","Vascular","Neuronal")
colnames(p_values) <- c("Hypoxic.1","Mes.","Mixed.1","Proneural","Hypoxic.2","Mixed.2","Vascular","Neuronal")
rownames(result)[1:6] <- c("M1","M2","S3","S4","S5","S6")
rownames(p_values)[1:6] <- c("M1","M2","S3","S4","S5","S6")

baseline <- rowSums(result)/sum(result)
df_normalized <- sweep(result, 2, colSums(result), FUN = "/")
enrichment_table <- df_normalized/baseline

enrichment_table$CellType <- rownames(enrichment_table)
p_values$CellType <- rownames(p_values)

significance_threshold <- 0.05
significant <- p_values
significant[, -ncol(significant)] <- significant[, -ncol(significant)] < significance_threshold

# Melt the enrichment data for plotting
enrichment_melted <- melt(enrichment_table, id.vars = "CellType", variable.name = "Niche", value.name = "Enrichment")

# Melt the significance data to match the enrichment table
significant_melted <- melt(significant, id.vars = "CellType", variable.name = "Niche", value.name = "Significant")

# Merge enrichment and significance data based on CellType and Niche
enrichment_melted$Significant <- significant_melted$Significant

levels <- levels(Aggregated_GBM@meta.data$Celltype_X)
enrichment_melted$CellType <- factor(enrichment_melted$CellType,
                                     levels = levels)
niche_levels <- c("Mes.","Hypoxic.1","Hypoxic.2","Proneural","Mixed.1","Mixed.2","Vascular","Neuronal")
enrichment_melted$Niche <- factor(enrichment_melted$Niche,
                                  levels = niche_levels)

pdf(sprintf("%s/Niche.Celltype.Enrichment.X.new.pdf", outdir), width = 12, height = 16); 
plot <- ggplot(enrichment_melted, aes(x = Niche, y = CellType)) +
  geom_point(aes(size = Enrichment, fill = CellType, stroke = ifelse(Significant == TRUE, 1, 0)), shape = 21) +  # Conditional border and color
  scale_size(range = c(2, 10), name = "Enrichment Score") +  # Adjust dot size based on enrichment
  scale_fill_manual(values = stepped.final.colors, name = "Cell Type") +  # Use the custom color vector
  labs(title = "Enrichment Dot Plot with Conditional Borders", x = "Niche", y = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot)
dev.off()

##########
#Figure 5G
##########

#Try to make edge node map
outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/BANKSY/BANKSY_Distributions"

set.seed(123)

aggregated.list <- list()

tumor.names <- c("GBM025A","GBM026I","GBM026M","GBM030","GBM034","GBM024D","GBM024E","GBM024H","GBM024I")
cell_names <- c(unique(Aggregated_GBM$BANKSY_Meta))

outdir <- "C:/Users/azwan/Box/Dunn Research/Projects/GBM_Core_Edge/Spatial/Xenium/Xenium_Analysis/Seurat_files/Aggregated_Second_Set/Final/BANKSY/BANKSY_Distributions"

for (a in 1:length(tumor.names)){
  final <- data.frame()
  file_path <- sprintf("%s/%s_Aggregated_Minimum_Distances.csv",outdir,tumor.names[a])
  
  if (file.exists(file_path)) {
    try({
      table <- read.csv(file_path)
      table <- table[,-1]
      long_table <- table %>%
        pivot_longer(
          cols = starts_with("Min_Distance"),  # Columns 2-4 with distances
          names_to = "Target",        # New column name for identities
          values_to = "AvgMinDistance"  # New column name for distances
        )
      long_table$sample <- tumor.names[a]
      aggregated.list <- rbind(aggregated.list, long_table)
    }, silent = FALSE) # Silent to suppress error messages
    
  }
}

aggregated.list <- as.data.frame(aggregated.list)
aggregated.list$Target <- str_remove(aggregated.list$Target, "^Min_Distance_to_")

# Summarize the data
summary_aggregated.list <- aggregated.list %>%
  group_by(BANKSY_X, Target) %>%
  summarise(
    MeanDistance = mean(AvgMinDistance, na.rm = TRUE),    # Calculate the mean distance
    MedianDistance = median(AvgMinDistance, na.rm = TRUE),# Calculate the median distance
    SD = sd(AvgMinDistance, na.rm = TRUE),                # Calculate the standard deviation
    SampleCount = n_distinct(sample)                # Count the number of distinct samples
  ) %>%
  ungroup()  # Ungroup the dataframe after summarizing

write.csv(summary_aggregated.list, sprintf("%s/Aggregated_Minimum_Distance.csv",outdir))
write.csv(aggregated.list, sprintf("%s/Individual_Total_Minimum_Distance.csv",outdir))

new <- summary_aggregated.list %>%
  mutate(pair_id = pmap_chr(list(BANKSY_X, Target), 
                            ~ paste(sort(c(..1, ..2)), collapse = "-")))

# For each pair (BANKSY_Meta, Target), keep the row with the lowest MedianDistance
filtered_data <- new %>%
  group_by(pair_id) %>%
  filter(MedianDistance == min(MedianDistance)) %>%
  ungroup()


# Create a graph object from the summarized dataframe
g <- graph_from_data_frame(filtered_data, directed = TRUE)

# Add edge attributes (e.g., MeanDistance, SD)
E(g)$mean_distance <- filtered_data$MeanDistance
E(g)$median_distance <- filtered_data$MedianDistance
E(g)$sd_distance <- filtered_data$SD
E(g)$sample_count <- filtered_data$SampleCount

E(g)$weight <- filtered_data$MeanDistance

layout <- layout_with_kk(g, weights = E(g)$weight)

pdf(sprintf("%s/Node.Edge.Plot.v4.X.pdf", outdir), width = 18, height = 12);
plot <- ggraph(g, layout = layout) + 
  geom_edge_link(aes(color = mean_distance, width = 0.5 ), 
                 alpha = 0.7)+
  geom_node_point(aes(size = 2), color = "grey") +
  scale_edge_color_gradient(low = "#EBDA3D", high = "#A121EB") +  # Adjust color gradient
  geom_node_text(aes(label = name), repel = TRUE) +
  geom_edge_link(aes(label = sample_count), alpha = 0.5) +
  
  theme_void()
print(plot)
dev.off()

##########
#Figure 6F
##########
#Old BANKSY naming convention used 
#X2E=Hypoxic.1
#X3E=Immune Infiltrated
#X1.4E=Mixed.1
#X1.5E=Stem-like
#X2.4=Hypoxic.2
#X3.5E=Mixed.2

outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/CellChat_X/ADM_CALCRL"

Idents(Aggregated_GBM) <- "BANKSY_X"

niche_names <- c("Vascular","X2E","X2.4")

niche <- subset(Aggregated_GBM, idents=niche_names)

#X2=M2
celltypes <- c("X2","Endothelial.1","Endothelial.2","Endothelial.4","Endothelial.5","Endothelial.KIT+","Proliferating.Endothelial",
               "aSMA.Mural.2","Pericyte.Mural.1","Proliferating.Mural","OPC","Fibroblast.2","Hypoxic.Myeloid.1")

Idents(niche) <- "Celltype_X"

filter <- subset(niche, idents=celltypes)

adm_expression <- FetchData(filter, vars = c("ADM", "Celltype_X", "BANKSY_X"))

# Step 2: Aggregate ADM expression by cell type and niche location
adm_agg <- adm_expression %>%
  group_by(Celltype_X, BANKSY_X) %>%
  summarise(mean_expression = mean(ADM, na.rm = TRUE))

# Step 3: Reshape data to wide format for heatmap
adm_matrix <- adm_agg %>%
  pivot_wider(names_from = BANKSY_X, values_from = mean_expression) %>%
  column_to_rownames("Celltype_X") %>%
  as.matrix()

# Step 4: Plot heatmap
jpeg(sprintf("%s/ADM_by_Celltype_BANKSY.jpg", outdir), width = 8, height = 12, units="in", res=300);
gg1 <- pheatmap(adm_matrix, 
                cluster_rows = TRUE, 
                cluster_cols = TRUE, 
                main = "ADM Expression by Cell Type and Niche Location",
                color=rev(plasma(100)))
print(gg1)
dev.off()

# Step 4: Plot heatmap
pdf(sprintf("%s/ADM_by_Celltype_BANKSY.pdf", outdir), width = 6, height = 12);
gg1 <- pheatmap(adm_matrix, 
                cluster_rows = TRUE, 
                cluster_cols = TRUE, 
                main = "ADM Expression by Cell Type and Niche Location",
                color=rev(plasma(100)))
print(gg1)
dev.off()

calcrl_expression <- FetchData(filter, vars = c("CALCRL", "Celltype_X", "BANKSY_X"))

# Step 2: Aggregate calcrl expression by cell type and niche location
calcrl_agg <- calcrl_expression %>%
  group_by(Celltype_X, BANKSY_X) %>%
  summarise(mean_expression = mean(CALCRL, na.rm = TRUE))

# Step 3: Reshape data to wide format for heatmap
calcrl_matrix <- calcrl_agg %>%
  pivot_wider(names_from = BANKSY_X, values_from = mean_expression) %>%
  column_to_rownames("Celltype_X") %>%
  as.matrix()

# Step 4: Plot heatmap
jpeg(sprintf("%s/calcrl_by_Celltype_BANKSY.jpg", outdir), width = 8, height = 12, units="in", res=300);
gg1 <- pheatmap(calcrl_matrix, 
                cluster_rows = TRUE, 
                cluster_cols = TRUE, 
                main = "calcrl Expression by Cell Type and Niche Location",
                color=rev(plasma(100)))
print(gg1)
dev.off()

# Step 4: Plot heatmap
pdf(sprintf("%s/calcrl_by_Celltype_BANKSY.pdf", outdir), width = 6, height = 12);
gg1 <- pheatmap(calcrl_matrix, 
                cluster_rows = TRUE, 
                cluster_cols = TRUE, 
                main = "calcrl Expression by Cell Type and Niche Location",
                color=rev(plasma(100)))
print(gg1)
dev.off()

##########
#Figure 7A
##########

# Function to calculate cell distribution within a distance range
calculate_distribution <- function(data, Celltype_Sample_Label_Final_Meta, min_distance, max_distance) {
  # Filter the data for the specified cell type
  target_cells <- data %>%
    filter(Celltype_Sample_Label_Final_Meta == !!Celltype_Sample_Label_Final_Meta)
  
  result <- data.frame()
  
  # Vectorize the calculation of distances
  for (i in 1:nrow(target_cells)) {
    x1 <- target_cells$x[i]
    y1 <- target_cells$y[i]
    
    # Calculate distances to all other cells
    distances <- sqrt((data$x - x1)^2 + (data$y - y1)^2)
    
    # Filter cells within the specified distance range
    within_range <- data %>%
      filter(distances >= min_distance & distances < max_distance)
    
    # Calculate the distribution of cell types within the range
    distribution <- within_range %>%
      group_by(Celltype_Sample_Label_Final_Meta) %>%
      summarise(count = n()) 
    
    distribution$cell_index <- i
    result <- bind_rows(result, distribution)
  }
  
  return(result)
}


outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/Distance_Matrices/Enrichment_Zonation"

tumor.names <- c("GBM025A","GBM026I","GBM026M","GBM030","GBM034")

Idents(Aggregated_GBM) <- "orig.ident"

for (a in 5:length(tumor.names)){
  
  object <- subset_opt(Aggregated_GBM, idents=tumor.names[a])
  
  coords <- GetTissueCoordinates(object)
  rownames(coords) <- coords$cell
  meta <- object@meta.data
  subset_meta <- subset(meta, select="Celltype_Sample_Label_Final_Meta")
  
  final_coords <- cbind(coords, subset_meta)
  final_coords <- subset(final_coords, select=-cell)
  
  cell_names <- c(unique(final_coords$Celltype_Sample_Label_Final_Meta))
  
  # Define distance ranges
  distance_ranges <- seq(0.5, 100.5, by = 20)
  
  # Calculate overall frequency of each cell type 
  
  overall <- table(object@meta.data$Celltype_Sample_Label_Final_Meta)
  overall <- as.data.frame(overall)
  colnames(overall) <-c("Celltype_Sample_Label_Final_Meta","Freq")
  overall$overall_frequency <- overall$Freq/sum(overall$Freq)
  
  # Calculate distributions for each distance range and average them
  
  for (b in 1:length(cell_names)){
    average_distributions <- data.frame()
    
    for (i in 1:(length(distance_ranges) - 1)) {
      min_distance <- distance_ranges[i]
      max_distance <- distance_ranges[i + 1]
      
      distribution <- calculate_distribution(final_coords, cell_names[b], min_distance, max_distance)
      
      averaged_distribution <- distribution %>%
        group_by(Celltype_Sample_Label_Final_Meta) %>%
        summarise(total_count = sum(count))
      
      averaged_distribution$distance_range <- paste(min_distance, "-", max_distance, "microns")
      average_distributions <- rbind(average_distributions, averaged_distribution)
    }
    
    # View the average distributions
    print(average_distributions)
    
    figure <- list()
    figure_new <- list()
    names <- unique(average_distributions$distance_range)
    
    final <- data.frame()
    
    for (i in 1:length(names)){
      figure[[i]] <- subset(average_distributions, distance_range==names[i])
      figure[[i]] <- as.data.frame(figure[[i]])
      figure[[i]]$average_percentage <- (figure[[i]]$total_count)/sum(figure[[i]]$total_count)
      figure_new[[i]] <- merge(figure[[i]], overall, by="Celltype_Sample_Label_Final_Meta",all.x=TRUE)
      figure_new[[i]]$enrichment <- figure_new[[i]]$average_percentage/figure_new[[i]]$overall_frequency
      final <- rbind(final, figure_new[[i]])
    }
    
    jpeg(sprintf("%s/%s_%s_Minimum_Distances_BarPlot.Enrichment.jpg", outdir, tumor.names[a],cell_names[b]), width = 12, height = 6, units="in", res=300);
    scatter <-ggplot(final, aes(fill=distance_range, y=enrichment, x=Celltype_Sample_Label_Final_Meta)) + 
      geom_bar(position="dodge", stat="identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(scatter);
    dev.off();
    
    write.csv(final, sprintf("%s/%s/%s_%s_Minimum_Distances_BarPlot.Enrichment.csv",outdir,"CSV_Files",tumor.names[a],cell_names[b]))
    
    print(sprintf("Finished_%s_%s",tumor.names[a],cell_names[b]))
  }
  print(sprintf("Finished_%s",tumor.names[a]))
}

outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Sequential/Final/Distance_Matrices/Enrichment_Zonation"

tumor.names <- c("GBM025A","GBM026I","GBM026M","GBM030","GBM034")
cell_names <- c(unique(Aggregated_GBM$Celltype_Sample_Label_Final_Meta))

final <- data.frame()

for (a in 1:length(tumor.names)){
  for (b in 1:length(cell_names)){
    file_path <- sprintf("%s/%s/%s_%s_Minimum_Distances_BarPlot.Enrichment.csv",outdir,"CSV_Files",tumor.names[a],cell_names[b])
    
    if (file.exists(file_path)) {
      try({
        table <- read.csv(file_path)
        table <- table[,-1]
        table$Root.Celltype <- cell_names[b]
        table$sample <- tumor.names[a]
        final <- rbind(final, table)
      }, silent = TRUE) # Silent to suppress error messages
    }
  }
}

outdir <- "/n/data1/mgh/neuro/petti/lab/Users/Anthony.Wang/GBM_Spatial/GBM_Xenium/Xenium_Second_Run/Aggregated_Both_Cohort/Cohort_1/Distance_Matrices/Enrichment_Zonation"

tumor.names <- c("GBM024D","GBM024E","GBM024H","GBM024I")
cell_names <- c(unique(Aggregated_GBM$Celltype_Sample_Label_Final_Meta))

for (a in 1:length(tumor.names)){
  for (b in 1:length(cell_names)){
    file_path <- sprintf("%s/%s/%s_%s_Minimum_Distances_BarPlot_Enrichment_Python.csv",outdir,"CSV_Files",tumor.names[a],cell_names[b])
    
    if (file.exists(file_path)) {
      try({
        table <- read.csv(file_path)
        names(table)[1:2] <- c("Celltype_Sample_Label_Final_Meta","total_count")
        table$Root.Celltype <- cell_names[b]
        table$sample <- tumor.names[a]
        final <- rbind(final, table)
      }, silent = FALSE) # Silent to suppress error messages
    }
  }
}

subset_final <- final %>% filter(distance_range != "80.5 - 100.5 microns")

#Remove outliers

# Function to remove outliers within each group
remove_outliers_grouped <- function(df, column) {
  df %>% group_by(distance_range, Celltype_Sample_Label_Final_Meta) %>%
    filter(enrichment >= (quantile(enrichment, 0.25) - 1.5 * IQR(enrichment)) & 
             enrichment <= (quantile(enrichment, 0.75) + 1.5 * IQR(enrichment))) %>%
    ungroup()
}

subset_final_filtered <- subset_final %>% filter(!(sample == "GBM024D" & enrichment == "503.2"))


# Create a plot for T cells
for (i in 6:6) {
  # Filter data for the current root cell type
  plot_data <- subset_final_filtered %>% filter(Root.Celltype == root_cell_types[i])
  
  # Remove outliers within each group
  plot_data <- remove_outliers_grouped(plot_data, "enrichment")
  plot_data <- plot_data %>% filter(Celltype_Sample_Label_Final_Meta != "Unidentifiable")
  
  
  # Filter for the closest distance range ("0-50"), calculate average enrichment score, and order cell types
  ordered_cell_types <- plot_data %>%
    filter(distance_range == "0.5 - 20.5 microns") %>%
    group_by(Celltype_Sample_Label_Final_Meta) %>%
    summarize(average_enrichment = mean(enrichment, na.rm = TRUE), .groups = 'drop') %>%
    arrange(desc(average_enrichment))
  
  # Print the ordered cell types
  print(ordered_cell_types)
  
  plot_data$Celltype_Sample_Label_Final_Meta <- factor(plot_data$Celltype_Sample_Label_Final_Meta,
                                                       levels = ordered_cell_types$Celltype_Sample_Label_Final_Meta)
  
  plot_data <- plot_data %>% filter(Celltype_Sample_Label_Final_Meta != "NA")
  
  first_5_categories <- head((ordered_cell_types$Celltype_Sample_Label_Final_Meta), 6)
  
  subset_plot <- plot_data %>% filter(Celltype_Sample_Label_Final_Meta %in% first_5_categories)
  
  # Create the plot
  pdf(sprintf("%s/%s/%s/%s_First.5.Remove.Outliers.pdf", outdir, "Average_Plots","Box_Plots_No_Outliers_PDF", root_cell_types[i]), width = 12, height = 6)
  
  p <- ggplot(subset_plot, aes(x = Celltype_Sample_Label_Final_Meta, y = enrichment, fill = distance_range)) +
    geom_boxplot(outlier.shape=NA) +  # Create boxplot
    geom_jitter(position = position_dodge(width = 0.75), alpha = 0.5, size = 2)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red")+scale_y_continuous(
      trans = 'log10',  # Use log scale if needed
      breaks = c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100),  # Specify the y-axis breaks
      labels = c("1","2","3","4","5","6","7","8","9","10","20","30","40","50","60","70","80","90","100"))+ 
    scale_fill_manual(values=c("#e6ecfa","#c0d0f3","#9bb4ec","#82a2e8"))
  
  # Print the plot
  print(p)
  dev.off()
}


