library(Seurat)
library(ggplot2)
library(tidyverse)
library(ggpubr)
sp <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sp.rds")


data <- sp@meta.data %>% filter(area == "Tumor core") %>%
        mutate(sample_area = paste0(sample,"_",area)) %>%
        select(sample_area,myCAFs) %>% na.omit()
data$sample_area <- factor(data$sample_area, 
                           levels = c("a_Tumor core", "d_Tumor core", "f_Tumor core", "h_Tumor core"))
my_comparisons <- list(
  c("a_Tumor core", "d_Tumor core"),
  c("a_Tumor core", "f_Tumor core"),
  c("a_Tumor core", "h_Tumor core"),
  c("d_Tumor core", "f_Tumor core"),
  c("d_Tumor core", "h_Tumor core"),
  c("f_Tumor core", "h_Tumor core")
)

p1 <- ggplot(data, aes(x = sample_area, y = myCAFs)) +
  geom_boxplot(aes(fill = sample_area), width = 0.6, outlier.shape = NA) +  # 隐藏离群点
  geom_jitter(width = 0.1, size = 0.5, alpha = 0.3) +  # 添加数据点
  stat_compare_means(
    comparisons = my_comparisons, 
    method = "wilcox.test",        # 使用t检验（若数据非正态/方差不齐改用"wilcox.test"）
    p.adjust.method = "bonferroni",  # 多重检验校正方法（可选：BH/fdr等）
    label = "p.signif",       # 显示星号(*)而非具体p值
    symnum.args = list(       # 自定义星号阈值
      cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
      symbols = c("***", "**", "*", "ns")
    ),
    step.increase = 0.08,     
    tip.length = 0.01,label.y = 0.4
  ) +
  scale_fill_brewer(palette = "Set2") +  # 配色方案
  labs(x = "Sample Area", y = "myCAFs Expression") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.line = element_line(colour = "black",size=0.8),
        axis.text = element_text(size=20))  
ggsave(plot = p1,filename = "sup_Figure5B_myCAFs_boxplot.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Supplementary/Figure5/",
       width = 8,height = 8)




Idents(sp) <- sp$sample
sp_list <- SplitObject(sp,split.by = "ident")
sample <- c("a","d","f","h")
for (i in sample) {
  Idents(sp_list[[i]]) <- "seurat_clusters"
  p1 <- VlnPlot(sp_list[[i]],features = "POSTN")
  ggsave(plot = p1,filename = paste0("sup_Vlnplot_",i,".pdf"),
         path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Supplementary/Figure5/",
         width = 9,height = 7.8)
}


sp <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sp.rds")
sp_clusters <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sp_clusters.rds")
sp$seurat_clusters <- sp_clusters$seurat_clusters
for (i in sample) {
  data <- sp@meta.data %>%
    mutate(sample_clusters = paste0(sample,"_",seurat_clusters)) %>%
    dplyr::filter(sample==i) %>%
    select(myCAFs,sample_clusters) %>%
    group_by(sample_clusters) %>%
    summarise(mean_value = mean(myCAFs,na.rm=T))
  p1 <- ggplot(data,aes(x=sample_clusters,y = mean_value)) +
    geom_bar(stat = "identity")
  
  ggsave(plot = p1,filename = paste0("sup_barplot","_myCAFs",i,".pdf"),
         path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Supplementary/Figure5/",
         width = 9,height = 7.8)
  
}
