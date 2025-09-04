library(Seurat)
library(tidyverse)
library(ggplot2)
library(tidyverse)
library(harmony)
sc <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/merge/merge_raw.rds")
Idents(sc) <- sc$celltype_cellchat
sc <- subset(sc,idents = "LSEC",invert=T)


Endo <- subset(sc,idents = "Endothelial_cells")
DefaultAssay(Endo) <- "RNA"
mat <- as.matrix(GetAssayData(Endo))
# harmony_exclude gene
mit_gene <- read.table('/data1/huanchangxiang/Mitochondria.txt',header = F)
heat_gene <- read.table('/data1/huanchangxiang/Heat_shock_protein.txt',header = F)
Rib_gene <- read.table('/data1/huanchangxiang/Ribosome.txt',header = F)
Dis_gene <- read.table('/data1/huanchangxiang/Dissociation.txt',header = F)
exclu_gene <- rbind(mit_gene,heat_gene,Rib_gene,Dis_gene)
rows_to_delete <- which(rownames(mat) %in% exclu_gene$V1)
nrow(mat)
mat <- mat[-rows_to_delete,]
nrow(mat)

############################ 第一次聚类-过滤 ################################## 
Endo <- CreateSeuratObject(mat) %>%
  AddMetaData(Endo@meta.data) %>%
  NormalizeData(.,normalization.method="LogNormalize",scale.factor=1e4) %>% 
  FindVariableFeatures(.) %>%
  ScaleData(.) %>%
  RunPCA(.,npcs = 30,features=VariableFeatures(object=.)) %>% 
  RunHarmony(c("patient")) %>%
  RunUMAP( reduction="harmony", dims=1:20,min.dist = 0.6) %>%           
  FindNeighbors(reduction="harmony", dims=1:20) %>%             
  FindClusters(resolution=0.4)

Endo_second <- subset(Endo,idents = c("2","7","10"),invert=T)

Endo <- Endo_second %>%
  NormalizeData(.,normalization.method="LogNormalize",scale.factor=1e4) %>% 
  FindVariableFeatures(.) %>%
  ScaleData(.) %>%
  RunPCA(.,npcs = 30,features=VariableFeatures(object=.)) %>% 
  RunHarmony(c("patient")) %>%
  RunUMAP( reduction="harmony", dims=1:20,min.dist = 0.6) %>%           
  FindNeighbors(reduction="harmony", dims=1:20) %>%             
  FindClusters(resolution=0.2)

saveRDS(Endo,file = "/data1/huanchangxiang/datasets/李莹雪论文/rds/Endo_second.rds")

Endo@meta.data$celltype_second <- Endo@meta.data$seurat_clusters
Endo@meta.data$celltype_second <- as.character(Endo@meta.data$celltype_second)
Endo@meta.data[Endo$celltype_second %in% c('1'),"celltype_second"] <- "ECs1"
Endo@meta.data[Endo$celltype_second %in% c('2'),"celltype_second"] <- "ECs2"
Endo@meta.data[Endo$celltype_second %in% c('3'),"celltype_second"] <- "ECs3"
Endo@meta.data[Endo$celltype_second %in% c('4'),"celltype_second"] <- "ECs4"
Endo@meta.data[Endo$celltype_second %in% c('5'),"celltype_second"] <- "ECs5"
Endo@meta.data[Endo$celltype_second %in% c('6'),"celltype_second"] <- "ECs6"
Endo@meta.data[Endo$celltype_second %in% c('0'),"celltype_second"] <- "ECs7"
DimPlot(Endo,label = T,group.by = "celltype_second")
saveRDS(Endo,file = "/data1/huanchangxiang/datasets/李莹雪论文/rds/Endo_second.rds")


sc@meta.data[rownames(Endo@meta.data),"celltype_cellchat"] <- Endo@meta.data[rownames(Endo@meta.data),"celltype_second"]

saveRDS(sc,file = "/data1/huanchangxiang/datasets/李莹雪论文/rds/merge/merge_raw_filter.rds")
