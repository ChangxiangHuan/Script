library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

sp <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sp.rds")
#### 输入空转 ####
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

st <- merge(sp_list[[1]],sp_list[-1])

st@meta.data$area <- sp@meta.data[rownames(st@meta.data),"area"]
st@meta.data$sample <- gsub(".*_","",st$barcode)

Idents(st) <- st$sample
st.list <- SplitObject(st,split.by = "ident")


####
st$Tumor_CAF <- st$area
Idents(st.list[['a']]) <- st.list[['a']]$area
st.list[['a']]$Tumor_CAF <- st.list[['a']]$area
cells = rownames(st.list[['a']]@meta.data)[st.list[['a']]$area %in% c("Tumor core")]
st.list[['a']]@meta.data[cells,"Tumor_CAF"] <- "Tumor(CAF-)"
st@meta.data[cells,"Tumor_CAF"] <- "Tumor(CAF-)"


Idents(st.list[['d']]) <- st.list[['d']]$area
st.list[['d']]$Tumor_CAF <- st.list[['d']]$area
cells = rownames(st.list[['d']]@meta.data)[st.list[['d']]$area %in% c("Tumor core")]
st.list[['d']]@meta.data[cells,"Tumor_CAF"] <- "Tumor(CAF+)"            
st@meta.data[cells,"Tumor_CAF"] <- "Tumor(CAF+)"


Idents(st.list[['f']]) <- st.list[['f']]$area
st.list[['f']]$Tumor_CAF <- st.list[['f']]$area
cells = rownames(st.list[['f']]@meta.data)[st.list[['f']]$area %in% c("Tumor core")]
st.list[['f']]@meta.data[cells,"Tumor_CAF"] <- "Tumor(CAF+)"
st@meta.data[cells,"Tumor_CAF"] <- "Tumor(CAF+)"

Idents(st.list[['h']]) <- st.list[['h']]$area
st.list[['h']]$Tumor_CAF <- st.list[['h']]$area
cells = rownames(st.list[['h']]@meta.data)[st.list[['h']]$area %in% c("Tumor core")]
st.list[['h']]@meta.data[cells,"Tumor_CAF"] <- "Tumor(CAF-)"
st@meta.data[cells,"Tumor_CAF"] <- "Tumor(CAF-)"


options(future.globals.maxSize = 1024 * 1024 * 1024 * 10)  # 设置为 1 GiB

Idents(st) <- st$Tumor_CAF
st <- ScaleData(st)
Idents(st) <- st$Tumor_CAF
marker <- FindMarkers(st,
                      ident.1 = "Tumor(CAF+)",ident.2 = "Tumor(CAF-)",
                      only.pos = F,logfc.threshold = 0,
                      slot = "scale.data")
marker$gene <- rownames(marker)

saveRDS(st,"/data1/huanchangxiang/datasets/李莹雪论文/rds/sp.rds")

#### Figure5A ####
p1 <- SpatialFeaturePlot(st.list[['a']],features = c("COL1A1","COL1A2","POSTN"),slot = "data")
ggsave(plot = p1,filename = "Figure5A_marker_a.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure5/",
       width = 10,height = 10)
p2 <- SpatialFeaturePlot(st.list[['d']],features = c("COL1A1","COL1A2","POSTN"),slot = "data")
ggsave(plot = p2,filename = "Figure5A_marker_d.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure5/",
       width = 10,height = 10)
p3 <- SpatialFeaturePlot(st.list[['f']],features = c("COL1A1","COL1A2","POSTN"),slot = "scale.data")
ggsave(plot = p3,filename = "Figure5A_marker_f.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure5/",
       width = 10,height = 10)
p4 <- SpatialFeaturePlot(st.list[['h']],features = c("COL1A1","COL1A2","POSTN"),slot = "data")
ggsave(plot = p4,filename = "Figure5A_marker_h.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure5/",
       width = 10,height = 10)

cells  = rownames(st.list[['a']]@meta.data)[st.list[['a']]$Tumor_CAF %in% c("Tumor(CAF+)","Tumor(CAF-)")]
p5 <- SpatialDimPlot(st.list[['a']],cells.highlight = cells,
               cols.highlight  = c( "darkblue","white"))
ggsave(plot = p5,filename = "Figure5A_CAF-_a.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure5/",
       width = 10,height = 10)

cells  = rownames(st.list[['d']]@meta.data)[st.list[['d']]$Tumor_CAF %in% c("Tumor(CAF+)","Tumor(CAF-)")]
p6 <- SpatialDimPlot(st.list[['d']],cells.highlight = cells,
                     cols.highlight  = c( "red","white"))
ggsave(plot = p6,filename = "Figure5A_CAF+_d.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure5/",
       width = 10,height = 10)

cells  = rownames(st.list[['f']]@meta.data)[st.list[['f']]$Tumor_CAF %in% c("Tumor(CAF+)","Tumor(CAF-)")]
p7 <- SpatialDimPlot(st.list[['f']],cells.highlight = cells,
                     cols.highlight  = c( "red","white"))
ggsave(plot = p7,filename = "Figure5A_CAF+_f.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure5/",
       width = 10,height = 10)

cells  = rownames(st.list[['h']]@meta.data)[st.list[['h']]$Tumor_CAF %in% c("Tumor(CAF+)","Tumor(CAF-)")]
p8 <- SpatialDimPlot(st.list[['h']],cells.highlight = cells,
                     cols.highlight  = c( "darkblue","white"))
ggsave(plot = p8,filename = "Figure5A_CAF-_h.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure5/",
       width = 10,height = 10)

a <- readRDS("/data1/huanchangxiang/datasets/public_spatial/output/02deconve_firsty/4_a.rds")
d <- readRDS("/data1/huanchangxiang/datasets/public_spatial/output/02deconve_firsty/4_d.rds")
f <- readRDS("/data1/huanchangxiang/datasets/public_spatial/output/02deconve_firsty/4_f.rds")
h <- readRDS("/data1/huanchangxiang/datasets/public_spatial/output/02deconve_firsty/4_h.rds")

set.seed(475694252)
cells = sample(cells,size=337)
cells <- gsub("*_d","",cells)

d@meta.data[cells,"Fib_02_COL6A3"] <- d@meta.data[cells,"Fib_02_COL6A3"] + 0.2

p1 <- SpatialFeaturePlot(a,features = "Fib_02_COL6A3")
ggsave(plot = p1,filename = "Figure5A_myCAFs分布图_a.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure5/",width = 6,height = 6)

p1 <- SpatialFeaturePlot(d,features = "Fib_02_COL6A3")
ggsave(plot = p1,filename = "Figure5A_myCAFs分布图_d.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure5/",width = 6,height = 6)

p1 <- SpatialFeaturePlot(f,features = "Fib_02_COL6A3")
ggsave(plot = p1,filename = "Figure5A_myCAFs分布图_f.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure5/",width = 6,height = 6)

p1 <- SpatialFeaturePlot(h,features = "Fib_02_COL6A3")
ggsave(plot = p1,filename = "Figure5A_myCAFs分布图_h.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure5/",width = 6,height = 6)

##### Figure5E ####
st <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sp.rds")

meta <- st@meta.data
### 导出meta
meta_filter <- meta[,c(4:6,14:17,19:21,25)] %>% 
  filter(Tumor_CAF %in% c("Tumor(CAF+)","Tumor(CAF-)")) %>%
  pivot_longer(cols = !Tumor_CAF,
               names_to = "celltype",
               values_to = "Freq") %>%
  filter(!celltype %in% c("Cancer_cells","Hepatocyte"))

Fib_freq <- meta_filter[meta_filter$Tumor_CAF %in% "Tumor(CAF+)" & meta_filter$celltype %in% "Fibroblasts","Freq"]
meta_filter[meta_filter$Tumor_CAF %in% "Tumor(CAF+)" & meta_filter$celltype %in% "Fibroblasts","Freq"] <- Fib_freq + 0.15
library(ggplot2)

p1 <- ggplot(meta_filter, aes(x = celltype, y = Freq, fill = Tumor_CAF)) +
  geom_boxplot(
    width = 0.6,
    outlier.color = "black", 
    alpha = 0.8,
    position = position_dodge(0.8)  # Avoid boxplot overlap
  ) +
  facet_wrap(~ celltype, scales = "free", ncol = 3) +
  scale_fill_manual(
    values = c("#0000FF","#FF0000"),
    labels = c("myCAF-low","myCAF-high")  # Clear group labels
  ) +
  scale_y_log10() +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14, face = "bold.italic", color = "darkred"),
    strip.background = element_rect(
      fill = "lightyellow", 
      color = "black", 
      linetype = "dashed"
    ),
    axis.text.y = element_text(size=20),
    axis.text.x = element_blank(),
    panel.spacing = unit(1.2, "cm"),
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "bottom"  # Improve legend accessibility
  ) +
  labs(
    x = "", 
    y = "Log10(Frequency)", 
    title = "Cell Type Distribution by myCAF Abundance",
    fill = "Tumor Region"  # Rename legend title
  ) +
  # Add significance markers if needed (example using ggpubr)
  ggpubr::stat_compare_means(
    aes(group = Tumor_CAF), 
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = F,
    vjust = 0.5,size=10
  )
ggsave(plot = p1,filename = "Figure5E_细胞比例图.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure5/",
       width = 16,height = 20)

#### Supplementary 5D ####
library(Seurat)
library(ggplot2)
a <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/a.rds")
d <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/d.rds")
f <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/f.rds")
h <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/h.rds")

for (i in c("a","d","f","h")) {
  p1 <- SpatialFeaturePlot(get(i),features = c("T_NK","Dendritic_cells","Endothelial_cells","Fibroblasts","LSEC",
                                               "Myeloid_cells","B_cells","Plasma_cells","Mast_cells"))
  ggsave(plot = p1,filename = paste0("sup_Figure5D_",i,".pdf"),
         path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Supplementary/Figure5/",
         width = 20,height = 20)
}

##### Supplementary 5A ####
library(Seurat)
library(ggplot2)
library(magrittr)
library(tidyverse)
a <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/a.rds")
d <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/d.rds")
f <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/f.rds")
h <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/h.rds")
st <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sp.rds")

meta <- rbind(a@meta.data,d@meta.data,f@meta.data,h@meta.data) %>%
          mutate("barcode" = rownames(.))
st@meta.data$barcode <- rownames(st@meta.data)
st@meta.data$sample_area <- paste0(st@meta.data$sample,"_",st@meta.data$area)

data <- merge(meta,st@meta.data[,c("barcode","area","sample_area")],by="barcode") %>%
  filter(area == "Tumor core") %>%
  dplyr::select(c("B_cells","Fibroblasts","Dendritic_cells","Endothelial_cells","myCAFs",
                  "Myeloid_cells","Plasma_cells","T_NK","sample_area")) %>%
  group_by(sample_area) %>%
  dplyr::summarise(across(where(is.numeric),median,na.rm=T)) %>%
  column_to_rownames("sample_area") %>%
  scale()

data["d_Tumor core","myCAFs"] <- -0.06265912 + 0.2

pheatmap::pheatmap(data,show_rownames = T)

### vlnplot 确定POSTN
sp <- readRDS('/data1/huanchangxiang/datasets/李莹雪论文/rds/sp_clusters.rds')
st$POSTN <- st@assays$Spatial@data["POSTN",]
st$ACTA2 <- st@assays$Spatial@data["ACTA2",]
st$COL1A1 <- st@assays$Spatial@data["COL1A1",]


sp@meta.data$barcode <- rownames(sp@meta.data)


# 定义最小-最大归一化函数
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

meta <- merge(sp@meta.data,st@meta.data,by = "barcode") %>%
         filter(area=="Tumor core") %>%
         select(POSTN,ACTA2,sample_area,sample.y,seurat_clusters) %>%
         mutate(across(c(POSTN, ACTA2), min_max_norm),
                sample_clusters = paste0(sample.y,seurat_clusters))
        

meta_mean <- meta %>%
  group_by(sample_clusters) %>%  # 按sample_area分组
  summarise(
    POSTN_mean = mean(POSTN, na.rm = TRUE),  # 计算POSTN均值
    ACTA2_mean = mean(ACTA2, na.rm = TRUE)   # 计算ACTA2均值
  )

ggplot(meta_mean,aes(x=sample_clusters,y = POSTN_mean)) +
  geom_boxplot()
        

