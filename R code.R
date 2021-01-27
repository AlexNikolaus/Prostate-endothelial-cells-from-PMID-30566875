library(Seurat)
library(tidyverse)
library(SCINA)

Prostate1 <- Read10X('D:/Endothelial scRNA-Seq project/PMID 30566875 prostate/GSE117403_D17')

#Create Seurat object
Prostate1 <- CreateSeuratObject(Prostate1, project = 'Prostate1', min.cells = 3, min.features = 3)

#Compute mt
Prostate1[['percent.mito']] <- PercentageFeatureSet(Prostate1, pattern = '^MT-')
VlnPlot(Prostate1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

#Subset & standard workflow
Prostate1 <- subset(Prostate1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 40)
Prostate1 <- NormalizeData(Prostate1, normalization.method = "LogNormalize", scale.factor = 10000)
Prostate1 <- FindVariableFeatures(Prostate1, selection.method = "vst", nfeatures = 2000)
Prostate1 <- ScaleData(Prostate1)
Prostate1 <- RunPCA(Prostate1, features = VariableFeatures(object = Prostate1))
Prostate1 <- FindNeighbors(Prostate1, dims = 1:10)
Prostate1 <- FindClusters(Prostate1, resolution = 0.1)
Prostate1 <- RunUMAP(Prostate1, dims = 1:10)
DimPlot(Prostate1, reduction = "umap", label = TRUE)
VlnPlot(Prostate1, features = c('CLDN5', 'VWF', 'FLT1', 'PECAM1'), ncol = 2)
Prostate1 <- subset(Prostate1, idents = c('4'))

#Create marker list
Endothelial_Markers <- list(c('CLDN5', 'VWF', 'FLT1', 'PECAM1'))
names(Endothelial_Markers) <- 'Endothelial cells'

#Run SCINA to label endothelial cells
SCINA_results <- SCINA(Prostate1@assays$RNA@data,
                       Endothelial_Markers,
                       max_iter = 2000, 
                       convergence_n = 100, 
                       convergence_rate = 0.999, 
                       sensitivity_cutoff = 0.9, 
                       rm_overlap=FALSE, 
                       allow_unknown=TRUE)
Prostate1$cell_labels <- SCINA_results$cell_labels
DimPlot(Prostate1,reduction = "umap", pt.size = 1, label = TRUE, group.by = 'cell_labels')

#Subset & write object
Prostate1 <- subset(Prostate1, cell_labels == 'Endothelial cells')
Prostate1$organ <- 'Prostate'
write_rds(Prostate1, 'Prostate1.rds')
