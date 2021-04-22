library(Seurat)
library(hdf5r)
library(sctransform)
library(ggplot2)
library(gridExtra)
setwd("/home/rstudio/scRNAseq_COVID/")

####set seed
set.seed(2020)

## download count matrices from GEO accession GSE145926, read in and store as Seurat objects
nomsf <- dir("GSE145926_RAW/")[grep(".h5", dir("GSE145926_RAW/"))]
lung_obj <- list()
for(x in nomsf[1:12]){
  print(x)
  olung <- Read10X_h5(filename = paste0("GSE145926_RAW/", x))
  nom <- gsub("_filtered_feature_bc_matrix.h5", "", x)
  olung <- CreateSeuratObject(counts = olung, project = nom, min.cells = 3, min.features = 200)
  print(olung)
  lung_obj[[nom]] <- olung
}
for (i in names(lung_obj)) {
  lung_obj[[i]] <- RenameCells(lung_obj[[i]],
                               add.cell.id = i)
}
# merge all into single object and store
merged_combined <- purrr::reduce(lung_obj,merge,do.normalize = FALSE)
save(merged_combined, file="data/All_Data_combined_Raw.RData")

#### load data 
load("data/All_Data_combined_Raw.RData") # 23,742 genes x 90,696 cells
##load associated metada  
load("Source_data/metadata_scrnaseq.RData")

# Rename samples for clarity
rownames(metadata) <- metadata$geo_accession
old_sample_name <- levels(Idents(merged_combined))
sele <- NULL; for(x in old_sample_name) sele <- c(sele, grep(strsplit(x,split="_")[[1]][1], rownames(metadata)))
metadata <- metadata[sele,] 
tipus <- as.factor(metadata$patient.group.ch1); levels(tipus) <- c("healthy", "mild", "severe", "severe")
new_sample_name <- paste(unlist(strsplit(old_sample_name,"_"))[c(1:12)*2], tipus, sep="_")
levels(merged_combined@active.ident) <- new_sample_name


#####Filtering process
#####
# vcalculate percent pf mitochondria genes
##housekeeping genes
hkgenes <- read.table("Source_data/tirosh_house_keeping.txt", skip = 2)
hkgenes <- as.vector(hkgenes$V1)
merged_combined[["percent.mt"]] <- PercentageFeatureSet(merged_combined, pattern = "^(MT|mt)-")


# diagnostic plots to help with filtering
p <- ggplot(merged_combined@meta.data, aes(nCount_RNA, nFeature_RNA, color=percent.mt))
p <- p + geom_point(size=0.2)
p <- p + geom_hline(yintercept=200, color='red')
p <- p + geom_hline(yintercept=6000, color='red')


hkgenes.found <- which(toupper(rownames(merged_combined@assays$RNA@data)) %in% hkgenes)
# Add_number_of_house_keeping_genes
n.expressed.hkgenes <- Matrix::colSums(merged_combined@assays$RNA@data[hkgenes.found, ] > 0)
merged_combined <- AddMetaData(object = merged_combined, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")
merged_combined <- subset(merged_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA > 1000 & n.exp.hkgenes > 55)



# scTransform pipeline ----------------------------------------------------

##normalize, scale and find variables in one command
merged_combined <- SCTransform(merged_combined, vars.to.regress = c("percent.mt","orig.ident"), 
                               return.only.var.gene=F,verbose = TRUE)

####PCA 
merged_combined<-RunPCA(merged_combined,verbose = T,features = )

# visualize loadings for PC1 and 2
p<-VizDimLoadings(merged_combined, dims = 1:2, reduction = "pca") # PC1: many interferon, TNFSF10 ; PC2: CD3 (T cell) markers


# UMAP
merged_combined <- RunUMAP(merged_combined, dims = 1:50) # increase threads to run faster?

# Clustering
merged_combined <- FindNeighbors(merged_combined, dims = 1:50)
merged_combined <- FindClusters(merged_combined, verbose = FALSE)
DimPlot(merged_combined, label = TRUE) + NoLegend()

### Figure 5A: identify clusters using cell type markers from reference publication
##### Supplementary Figure 9A: UMAP representation of the dataset
DimPlot(object = merged_combined, reduction = 'umap', group.by = 'seurat_clusters', label=TRUE, label.size = 8) + NoLegend()

##### Figure 5B: Plot markers per cluster
features <- c("TPPP3","KRT18","CD68","FCGR3B","CD1C","CLEC9A","LILRA4","TPSB2","CD3D","KLRD1","MS4A1","IGHG4",
              "FCN1","SPP1","FABP4")
DotPlot(merged_combined, features=features) + NoLegend()

##### Cluster identity based on the above features:
macrophages_1 <- c(32,30,29,26,20,13,12,7,6,1,2,3)
macrophages_2 <- c(15)
macrophages_3 <- c(8)
macrophages_4 <- c(22,16,14,10,4,0)
Ts <- c(5,11,18)
epithelial <- c(23,31,17)
others <- c(33,25)
Bs <- c(27,21,24)

merged_combined@meta.data$metaclusters[merged_combined$seurat_clusters %in% macrophages_1] <- "Macrophages FCN1hi"
merged_combined@meta.data$metaclusters[merged_combined$seurat_clusters %in% macrophages_2] <- "Macrophages FCN1loSPP1+"
merged_combined@meta.data$metaclusters[merged_combined$seurat_clusters %in% macrophages_3] <- "Macrophages FCN1+SPP1+"
merged_combined@meta.data$metaclusters[merged_combined$seurat_clusters %in% macrophages_4] <- "Macrophages FABP4+"
merged_combined@meta.data$metaclusters[merged_combined$seurat_clusters == 19] <- "mDC"
merged_combined@meta.data$metaclusters[merged_combined$seurat_clusters == 28] <- "pDC"
merged_combined@meta.data$metaclusters[merged_combined$seurat_clusters %in% Ts] <- "T cells"
merged_combined@meta.data$metaclusters[merged_combined$seurat_clusters == 9] <- "NK cells"
merged_combined@meta.data$metaclusters[merged_combined$seurat_clusters == 23] <- "Epithelial Ciliated cells"
merged_combined@meta.data$metaclusters[merged_combined$seurat_clusters %in% c(31,17)] <- "Epithelial Secretory cells"
merged_combined@meta.data$metaclusters[merged_combined$seurat_clusters %in% others] <- "Others"
merged_combined@meta.data$metaclusters[merged_combined$seurat_clusters == 27] <- "B cells"
merged_combined@meta.data$metaclusters[merged_combined$seurat_clusters %in% c(21,24)] <- "Plasma cells"

# look for cell type and relevant cell markers
# FeaturePlot(merged_combined, features = c("CD80", "CD86", "CTLA4", "IL6","MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ","CD8A"), ncol=2)

# Figure 6 A
FeaturePlot(merged_combined, features = c("CD80", "CD86"), ncol=2)

# Figure 6 B
FeaturePlot(merged_combined, features = c("IL6"), ncol=1)

# Figure 6 C
FeaturePlot(merged_combined, features = c("IL2RA", "CTLA4","CD8A","FOXP3"), ncol=2)

# CD80/86 axis analysis ---------------------------------------------------

status <- tipus
sampleId <- as.factor(merged_combined@meta.data$orig.ident)
clusterAssign <- Idents(merged_combined)
sampleNames <- levels(sampleId)
Tcell_cluster <- "5" # active CD4+ T cell cluster
marker_data <- FetchData(merged_combined, c("IL6", "CD80", "CD86"))

# Number of cells per individual
TCD4_active <- IL6 <- CD80 <- CD86 <- rep(NA, length(sampleNames))
for(i in 1:length(sampleNames)){
  ix <- which(sampleId == sampleNames[i])
  # number of active T cells per sample
  clusterIds <- NULL
  clusterIds <- which(clusterAssign == Tcell_cluster)
  TCD4_active[i] <- length(intersect(ix, clusterIds))
  # number of CD80 and CD86 expressing cells per sample
  CD86[i] <- length(intersect(ix, which(marker_data[,3] > 0)))
  CD80[i] <- length(intersect(ix, which(marker_data[,2] > 0)))
  # number of IL6 producing cells per sample
  IL6[i] <- length(intersect(ix, which(marker_data[,1] > 0)))
}


# Association analysis between the elements of the CD80/86 axis in BALF:
#       CD80/86+ APCs, active CD4+ T cells, and IL6 producing cells
summary(lm(log(datcell$IL6)~log(datcell$TCD4_active)+log(datcell$nTOT)))   
lmA <- lm(log(datcell$IL6)~log(datcell$TCD4_active)+log(datcell$nTOT))
summary(lm(log(datcell$CD80)~log(datcell$TCD4_active)+log(datcell$nTOT)))   
lmB <- lm(log(datcell$CD80)~log(datcell$TCD4_active)+log(datcell$nTOT))
summary(lm(log(datcell$CD86)~log(datcell$TCD4_active)+log(datcell$nTOT)))  
lmC <- lm(log(datcell$CD86)~log(datcell$TCD4_active)+log(datcell$nTOT))
summary(lm(log(datcell$IL6)~log(datcell$CD80)+log(datcell$nTOT))) 
lmD <- lm(log(datcell$IL6)~log(datcell$CD80)+log(datcell$nTOT))
summary(lm(log(datcell$IL6)~log(datcell$CD86)+log(datcell$nTOT))) 


# Figure 7 plots ----------------------------------------------------------

axisData <- data.frame(sampleNames, TCD4_active, CD86, CD80, IL6, status)

## Figure 7A
fill = c("steelblue", "coral", "red")
p8A <- ggplot(axisData, aes(x = CD86 , y = TCD4_active, size=TCD4_active, fill=status)) + geom_point(shape=21) + scale_size_area(max_size = 10)
p8A <- p8A + ggtitle("") + labs(x = "#CD86 cells", y = "#Active T CD4+ cells") + scale_fill_manual(values = fill)
lm_fitA <- lm(TCD4_active ~ CD86, data=datcell)
p8A <- p8A + geom_abline(color="black", slope = lm_fitA$coefficients[2], intercept = lm_fitA$coefficients[1])
p8A


## Figure 7B
fill = c("steelblue", "coral", "red")
p8B <- ggplot(axisData, aes(x = CD80 , y = TCD4_active, size=TCD4_active, fill=status)) + geom_point(shape=21) + scale_size_area(max_size = 10)
p8B <- p8B + ggtitle("") + labs(x = "#CD80 cells",  y = "#Active T CD4+ cells" ) + scale_fill_manual(values = fill)
lm_fitB <- lm(TCD4_active ~ CD80, data=datcell)
p8B <- p8B + geom_abline(color="black", slope = lm_fitB$coefficients[2], intercept = lm_fitB$coefficients[1])
p8B


## Figure 7C
fill = c("steelblue", "coral", "red")
p8C <- ggplot(axisData, aes( x = TCD4_active, y = IL6, size=IL6, fill=status)) + geom_point(shape=21) + scale_size_area(max_size = 10)
p8C <- p8C + ggtitle("") + labs(x = "#Active T CD4+ cells", y = "#IL6+ cells") + scale_fill_manual(values = fill)
lm_fitC <- lm(IL6 ~ TCD4_active, data=datcell)
p8C <- p8C + geom_abline(color="black", slope = lm_fitC$coefficients[2], intercept = lm_fitC$coefficients[1])
p8C

## Figure 7D
fill = c("steelblue", "coral", "red")
p8D <- ggplot(axisData, aes(x = CD86, y = IL6, size=IL6, fill=status)) + geom_point(shape=21) + scale_size_area(max_size = 10)
p8D <- p8D + ggtitle("") + labs(x = "#CD86 cells", y = "#IL6+ cells") + scale_fill_manual(values = fill)
lm_fitD <- lm(IL6 ~ CD86, data=datcell)
p8D <- p8D + geom_abline(color="black", slope = lm_fitD$coefficients[2], intercept = lm_fitD$coefficients[1])
p8D

## Figure 7E
fill = c("steelblue", "coral", "red")
p8E <- ggplot(axisData, aes(x = CD80 , y = IL6, size=IL6, fill=status)) + geom_point(shape=21) + scale_size_area(max_size = 10)
p8E <- p8E + ggtitle("") + labs(x = "#CD80 cells", y = "#IL6+ cells") + scale_fill_manual(values = fill)
lm_fitE <- lm(IL6 ~ CD80, data=datcell)
p8E <- p8E + geom_abline(color="black", slope = lm_fitE$coefficients[2], intercept = lm_fitE$coefficients[1])
p8E
