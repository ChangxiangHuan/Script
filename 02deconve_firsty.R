library(Seurat)
library(ggplot2)
library(stringr)
library(harmony)
library(patchwork)
library(dplyr)
library(spacexr)
library(tidyverse)
.libPaths("/data1/huanchangxiang/R_library/")
options(scipen = 100)
datasets <- c("a","d","f","h")


################################ Reference #####################################
scRNA <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/merge/merge_raw.rds")
Idents(scRNA) <- scRNA$celltype
scRNA$barcode <- rownames(scRNA@meta.data)
# matrix for reference
sc_counts <- as.matrix(scRNA[['RNA']]@counts)
# UMI for reference
sc_nUMI <- colSums(sc_counts)
# Celltype
cellType <- data.frame(barcode=scRNA$barcode,celltype=scRNA$celltype)
cell_types <- cellType$celltype 
names(cell_types) <- cellType$barcode
cell_types <- as.factor(cell_types)
cell_types
reference <- Reference(sc_counts,cell_types,sc_nUMI)

### 空转处理
for (data in datasets) {
  spatial <- Load10X_Spatial(paste0("/data1/huanchangxiang/datasets/spatial/hcc/HCC_spatial_",data,"/outs/"))
  spatial@images$slice1@coordinates[,2:5] <- lapply(spatial@images$slice1@coordinates[,2:5],as.numeric)
  ########################### RCTD for any samples ###############################
  coords <- Seurat::GetTissueCoordinates(spatial)
  colnames(coords) <- c('x','y')
  
  query <- SpatialRNA(coords = coords,
                      counts = spatial@assays$Spatial@counts,
                      nUMI = structure(spatial$nCount_Spatial,names=colnames(spatial)))
  
  RCTD <- create.RCTD(spatialRNA = query,reference = reference,max_cores = 50)
  RCTD <- run.RCTD(RCTD,doublet_mode = 'full')
  result <- data.frame(RCTD@results$weights)
  spatial <- AddMetaData(spatial,metadata = result)
  rownames(spatial@meta.data) <- paste0(rownames(spatial@meta.data),"_",data)
  saveRDS(spatial,file = paste0("/data1/huanchangxiang/datasets/李莹雪论文/rds/",data,"_celltype.rds"))
}
