library(Seurat)
library(ggplot2)
library(tidyverse)



marker_genes <- c("TF","COL3A1","DNASE1L3","VIM","ALB","CD3E","MS4A1","CD68")
samples <- c("a","d","f","h")
sp_list <- lapply(samples, function(sample){
  sp <- readRDS(paste0("/data1/huanchangxiang/datasets/李莹雪论文/rds/",sample,".rds"))
  colnames(sp@assays$Spatial@counts) <- paste0(colnames(sp@assays$Spatial@counts),"_",sample)
  colnames(sp@assays$Spatial@data) <- paste0(colnames(sp@assays$Spatial@data),"_",sample)
  rownames(sp@meta.data) <- gsub("-1.*", "-1", rownames(sp@meta.data))
  rownames(sp@meta.data) <- paste0(rownames(sp@meta.data),"_",sample)
  sp@meta.data$barcode   <- rownames(sp@meta.data)
  rownames(sp@images$slice1@coordinates) <- paste0(rownames(sp@images$slice1@coordinates),
                                                   "_",sample)
  return(sp)
})
names(sp_list) <- samples
for (i in samples) {
  
  sp <- sp_list[[i]]
  for (gene in marker_genes) {
    
    p1 <- SpatialFeaturePlot(sp,features = gene)
    ggsave(plot = p1,filename = paste0("sup_Figure2_",i,"_",gene,".pdf"),
           width = 8,height = 7,
           path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Supplementary/Figure2/")
    
  }
  
}
 



