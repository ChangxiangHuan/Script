library(Seurat)
library(ggplot2)

sc <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/merge/merge_raw.rds")


Idents(sc) <- sc$NorT
p1 <- DimPlot(sc)
ggsave(plot = p1,filename = "Figure1_UAMP_NorT.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Supplementary/Figure1/",
       width = 7,height = 6.4)

Idents(sc) <- sc$sample
p1 <- DimPlot(sc)
p1
ggsave(plot = p1,filename = "Figure1_UAMP_Patient.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Supplementary/Figure1/",
       width = 7,height = 6.4)

sample <- c("a","d","f","h")
for (x in sample) {
  sp <- readRDS(paste0("/data1/huanchangxiang/datasets/李莹雪论文/rds/",x,".rds"))
  p1 <- SpatialFeaturePlot(sp,features = "Cancer_cells")
  p2 <- SpatialFeaturePlot(sp,features = "Hepatocyte")
  ggsave(plot = p1,filename = paste0(x,"_spatial_Cancer_cells.pdf"),
         path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Supplementary/Figure1/",
         width = 7,height = 6.4)
  ggsave(plot = p2,filename = paste0(x,"_spatial_Hepatocyte.pdf"),
         path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Supplementary/Figure1/",
         width = 7,height = 6.4)
}

