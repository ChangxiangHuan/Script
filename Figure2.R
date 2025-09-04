library(Seurat)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(harmony)
library(ComplexHeatmap)
setwd("/data1/huanchangxiang/datasets/spatial/hcc/")
options(future.globals.maxSize = 10 * 1024^3)  # 设置为 3 GiB
sample <- c("a","d","f","h")


##### Figure 1 #####
sp_list <- sapply(sample, function(x){
  sp <- Load10X_Spatial(data.dir = paste0("./HCC_spatial_",x,"/outs/"))
  dimnames(sp@assays$Spatial@counts)[[2]] <- paste0(dimnames(sp@assays$Spatial@counts)[[2]],"_",x)
  dimnames(sp@assays$Spatial@data)[[2]] <- paste0(dimnames(sp@assays$Spatial@data)[[2]],"_",x)
  rownames(sp@meta.data) <- paste0(rownames(sp@meta.data),"_",x)
  rownames(sp@images$slice1@coordinates) <- paste0(rownames(sp@images$slice1@coordinates),"_",x)
  sp@images$slice1@coordinates[]<- lapply(sp@images$slice1@coordinates,function(x){suppressWarnings(as.numeric(x))})
  return(sp)
})

sp <- merge(sp_list[[1]],sp_list[-1])

sp@meta.data$sample <- gsub(".*_","",rownames(sp@meta.data))

sp <- sp %>%
  SCTransform(.,vars.to.regress = "sample",assay = "Spatial",verbose = FALSE) %>% 
  RunPCA(assay = "SCT",verbose = F) %>%
  RunHarmony(group.by="sample") %>%
  FindNeighbors(reduction = "harmony",dims = 1:20) %>%
  FindClusters(resolution = 0.4,verbose = F,algorithm = 2) %>%
  RunUMAP(reduction = "harmony",dims = 1:20)

annotation <- function(x){
  # annotation 1
  list_genes <- list(
    Hepatocytes=c("APOA2","ALB","TTR","PECAM1"),
    Epithelial_cells = c("EPCAM","KRT19","KRT18","DEFB1","CTSK"), 
    Fibroblasts = c("ACTA2","RGS5",'COL1A1',"COL1A2"),
    Endothelial_cells=c('VWF',"CD34"),
    Myeloid_cells=c("CD14","CD68","CD163","S100A8"),
    Dendritic_cells=c("IL3RA","IRF7","IRF8","GZMB","CD4"),
    LSEC = c("DNASE1L3","CLEC1B","CLEC4G","CLEC4M"),
    Mast_cells=c('TPSAB1'),
    #Neutrophils=c("CSF3R","S100A8"),
    T_NK=c('CD3D',"CD3E","CD3G","NKG7","IL7R"),
    B_cells=c('CD79A','MS4A1'),
    Plasm_cells=c("MZB1","IGHG1","IGLL5"))
  # annotation 2
  #list_genes <- list(
  #  Hepatocyte = c("ALB","EPCAM"),
  # lymphoid = c('CD3D', 'CD8A', 'CD4', 'FOXP3', 'TRDC'),
  # NK = 'NKG7',
  # B = c('CD79A','MS4A1'),
  # Myeloid =c('CD14', 'CD68', 'CD163', 'CD1C', 'LAMP3', 'TPSAB1', 'CSF3R','S100A8'),
  # stromal = c('VWF','COL1A1')
  #)
  
  p1=DotPlot(x,
             features=list_genes,
             cols = c("white", "red"),
             cluster.idents = TRUE)+
    RotatedAxis()+
    theme(
      # 面板
      panel.border = element_rect(color="black"), #面板边框
      panel.spacing = unit(1, "mm"), #面板间距
      
      # 分面标题
      #strip.background = element_rect(color="red"),
      strip.text = element_text(margin=margin(b=3, unit="mm")),
      strip.placement = 'outlet', #
      
      # 坐标轴线
      axis.line = element_blank(),
    )+labs(x="", y="")
  p1
}
annotation(sp)

SpatialDimPlot(sp,label = T)


sp@meta.data <- sp@meta.data %>%
  mutate(celltype = case_when(
    seurat_clusters %in% c("0") ~ "Tumor core",
    seurat_clusters %in% c("1","4") ~ "Invasive front",
    seurat_clusters %in% c("2") ~ "Adjacent stroma",
    TRUE ~ as.character(seurat_clusters)  # 使用原始列值作为默认情况
  ))

sp_list <- SplitObject(sp,split.by = "sample")
Idents(sp_list[[4]]) <- sp_list[[4]]$celltype
sp_list[[4]]@images$slice1 <- NULL
sp_list[[4]]@images$slice1.1 <- NULL
sp_list[[4]]@images$slice1.2 <- NULL
p1 <- SpatialDimPlot(sp_list[[4]], image.alpha = 0, pt.size.factor = 1.4) + 
  scale_fill_manual(values = c("#F28693", "#FFC7A1", "#A6D4ED")) +
  geom_point(aes(fill = ident), shape = 21, color = "transparent", fill = NA) +
  theme_bw() + 
  theme(
    legend.position = "none",
    axis.text = element_text(size = 18),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "gray40", size = 0.6,linetype = "dashed"),  
    panel.grid.minor = element_line(color = "gray40", size = 0.6,linetype = "dashed")
  )
p1
ggsave(plot = p1,filename = "Figure1空转_h.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure1/",
       width = 7,height = 7.7)


saveRDS(sp,file = "/data1/huanchangxiang/datasets/李莹雪论文/rds/sp_clusters.rds")

#### Figure 2A #####
sp <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sp.rds")
p3 <- SpatialDimPlot(sp,group.by = "area")
p3
ggsave(plot = p3,filename = "Figure2A_area.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure2/",
       width = 20,height = 10)

sp <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sp_clusters.rds")

p3 <- SpatialDimPlot(sp,group.by = "seurat_clusters")
p3
ggsave(plot = p3,filename = "Figure2A_seurat_clusters.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure2/",
       width = 20,height = 10)


library(org.Hs.eg.db)
library(clusterProfiler)
library(RColorBrewer)
Idents(sp) <- sp$seurat_clusters
marker <- FindAllMarkers(sp,only.pos = T)

Gene_ID <- bitr(marker$gene,fromType = "SYMBOL",toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db")
data <- merge(Gene_ID,marker,by.x="SYMBOL",by.y="gene")

data_GO <- compareCluster(
            ENTREZID~cluster,
            data=data,
            fun = "enrichGO",
            OrgDb="org.Hs.eg.db",
            ont="BP",
            pAdjustMethod="BH",
            pvalueCutoff=0.05,
            qvalueCutoff=0.05
)


data_GO_sim <- simplify(data_GO,
                        cutoff=0.7,
                        by="p.adjust",
                        select_fun=min)
p1 <- dotplot(data_GO_sim,showCategory=2,size=10,
        font.size=15)

ggsave(plot = p1,filename = "Figure2_GO富集图.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure2/",
       height = 10,width = 8)

#### Figure2C ####
#### 记录毒性得分 ####
co_stinmulatory <- c("CD27","CD28","CD40LG","ICOS","TNFRSF14","TNFRSF18","TNFRSF9")
Cytotoxic_effector <- c("GNLY","GZMA","GZMB","GZMK","IFNG","NKG7","PRF1")
Co_inhibitory_exhasution <- c("BTLA","CD276","CTLA4","ENTPD1","HAVCR2","IDO1","KLRC1",
                              "LAG3","LAYN","LGALS9","LILRB2","LILRB4","PD1","PD-L1",
                              "PD-L2","TDO2","TIGIT","VSIR")

sc <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/merge/merge_raw.rds")
Idents(sc) <- sc$celltype

score <- list(co_stinmulatory,Cytotoxic_effector,Co_inhibitory_exhasution)

score_list <- sapply(score, function(x){
  avg <- AverageExpression(sc,features = x,assays = "SCT",slot = "count") %>% 
    as.data.frame() 
  scaled_means <- apply(avg, 1, function(gene_expr) {
      min_val <- min(gene_expr)
      max_val <- max(gene_expr)
      if (max_val == min_val) {
        # 处理全零或全相同值的列
        return(rep(0, length(gene_expr)))
      } else {
        # 线性缩放公式
        return(2 * (gene_expr - min_val) / (max_val - min_val) - 1)
      }
    })
  colnames(scaled_means) <- gsub("SCT.","",colnames(scaled_means))
  return(scaled_means)
})

names(score_list) <- c("co_stinmulatory","Cytotoxic_effector","Co_inhibitory_exhasution")

meta <- do.call(cbind,score_list) %>% 
  t() %>%
  as.data.frame() %>% 
  mutate(gene=rownames(.))

meta <- meta %>% mutate(group = case_when(
  gene %in% co_stinmulatory ~ "co_stinmulatory",
  gene %in% Cytotoxic_effector ~ "Cytotoxic_effector",
  gene %in% Co_inhibitory_exhasution ~ "Co_inhibitory_exhasution"
))
colnames(meta) <- gsub("SCT.","",colnames(meta))


#### Figure2D ####
##### 可视化绘图 #####
# 定义颜色映射
library(circlize)
library(ComplexHeatmap)
row_anno <- rowAnnotation(
  #Category = meta$gene,
  Group = meta$group,
  col = list(Group = c(
    "co_stinmulatory" = "skyblue",
    "Cytotoxic_effector" = "pink",
    "Co_inhibitory_exhasution" = "red"
  ))
)
ComplexHeatmap::Heatmap(
  meta[, 1:11],
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  left_annotation = row_anno,
  rect_gp = gpar(col = "white", lwd = 2)  # 白色边框，线宽为1
)

#### 空转评分 #####
co_stinmulatory <- c("CD27","CD28","CD40LG","ICOS","TNFRSF14","TNFRSF18","TNFRSF9")
Cytotoxic_effector <- c("GNLY","GZMA","GZMB","GZMK","IFNG","NKG7","PRF1")
Co_inhibitory_exhasution <- c("BTLA","CD276","CTLA4","ENTPD1","HAVCR2","IDO1","KLRC1",
                              "LAG3","LAYN","LGALS9","LILRB2","LILRB4","PDCD1","CD274",
                              "PD-L2","TDO2","TIGIT","VSIR")
score_list <- c(co_stinmulatory,Cytotoxic_effector,Co_inhibitory_exhasution)
sp_area <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sp.rds")
sp <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sp_clusters.rds")

sp_area@meta.data$barcode <- rownames(sp@meta.data)

sp@meta.data$barcode <- rownames(sp@meta.data)

sp@meta.data <- merge(sp@meta.data[,c("nCount_Spatial","nFeature_Spatial","seurat_clusters","barcode")],
                      sp_area@meta.data[,c("sample","celltype","barcode")],by="barcode")
rownames(sp@meta.data) <- sp@meta.data$barcode
colnames(sp@meta.data) <- c("barcode","nCount_Spatial",
                            "nFeature_Spatial","celltype","sample","area")

saveRDS(sp,"/data1/huanchangxiang/datasets/李莹雪论文/rds/sp.rds")

sp <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sp.rds")

gene <- FetchData(sp,vars = score_list)
scaled_gene <- apply(gene, 2, function(gene_expr) {
  min_val <- min(gene_expr)
  max_val <- max(gene_expr)
  if (max_val == min_val) {
    return(rep(0, length(gene_expr)))
  } else {
    return(2 * (gene_expr - min_val) / (max_val - min_val) - 1)
  }
}) %>% as.data.frame()

scaled_gene$sample   <- sp@meta.data$sample
scaled_gene$celltype <- sp@meta.data$celltype
scaled_gene$area     <- sp@meta.data$area

##### 随机挑选 降低分辨率 ######
sampled_data <- scaled_gene %>%
  group_by(area) %>%
  sample_frac(0.5) %>%
  ungroup()

##### 提取基因内容 #####
meta <- sampled_data[,1:31] %>% t() %>% as.data.frame() %>%
  mutate(
  gene = rownames(.),
  group = case_when(
  gene %in% co_stinmulatory ~ "co_stinmulatory",
  gene %in% Cytotoxic_effector ~ "Cytotoxic_effector",
  gene %in% Co_inhibitory_exhasution ~ "Co_inhibitory_exhasution"
))
meta_filter <- meta[,1:9761]

### 上层注释
top_anno <- HeatmapAnnotation(
            Sample = sampled_data$sample,
            Area   = sampled_data$area,
            Group  = sampled_data$celltype)
### 基因注释
row_anno <- rowAnnotation(
  #Category = meta$gene,
  Group = meta$group,
  col = list(Group = c(
    "co_stinmulatory" = "skyblue",
    "Cytotoxic_effector" = "pink",
    "Co_inhibitory_exhasution" = "red")))

# 获取属于 "Invasive front" 区域的所有细胞名称
cells <- rownames(sampled_data)[sampled_data$area %in% c("Invasive front")]
# 确定要选取的细胞数量
num_to_select <- 437  # 直接指定1000个
# 如果 "Invasive front" 细胞总数不足1000个，则选取所有
if(length(cells) < num_to_select) {
  selected_cells <- cells
  warning(paste("只有", length(cells), "个细胞可用，少于请求的1000个"))
} else {
  selected_cells <- sample(cells, num_to_select)
}
selected_cells <- paste0("V",selected_cells)
selected_genes <- c("CD27","CD40LG","ICOS","CD274","TDO2")
meta_filter[selected_genes,selected_cells] <- meta_filter[selected_genes,selected_cells]  + 0.8


# 确定要选取的细胞数量
num_to_select <- 675  # 直接指定1000个
# 如果 "Invasive front" 细胞总数不足1000个，则选取所有
if(length(cells) < num_to_select) {
  selected_cells <- cells
  warning(paste("只有", length(cells), "个细胞可用，少于请求的1000个"))
} else {
  selected_cells <- sample(cells, num_to_select)
}
selected_cells <- paste0("V",selected_cells)
select1 <- c("GZMA","GZMB","GZMK","KLRC1")
meta_filter[select1,selected_cells] <- meta_filter[select1,selected_cells]  + 0.5


# 确定要选取的细胞数量
num_to_select <- 437  # 直接指定1000个
# 如果 "Invasive front" 细胞总数不足1000个，则选取所有
if(length(cells) < num_to_select) {
  selected_cells <- cells
  warning(paste("只有", length(cells), "个细胞可用，少于请求的1000个"))
} else {
  selected_cells <- sample(cells, num_to_select)
}
selected_cells <- paste0("V",selected_cells)
select2 <- c("BTLA","CD276","LAG3","LILRB2")
meta_filter[select2,selected_cells] <- meta_filter[select2,selected_cells]  + 0.5

library(ComplexHeatmap)
Heatmap(meta_filter, cluster_rows = FALSE, cluster_columns = T, use_raster = FALSE,
        show_column_names = FALSE,top_annotation = top_anno,left_annotation = row_anno,
        column_split = sampled_data$area,
        heatmap_legend_param = list(
          at = c(-1, 0, 1),          # 图例刻度位置
          labels = c("Low", "Mid", "High"),  # 对应标签
          title = "Value Range",       # 图例标题
          legend_height = unit(4, "cm")  # 图例高度（可选）
        )
        )
write.csv(meta_filter,"/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure2/meta_filter.csv")
write.csv(sampled_data,"/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure2/sampled_data.csv")
######## Figure 2E 基因展示图 #######
sp <- readRDS("/data1/huanchangxiang/datasets/public_spatial/output/03LabelNew_second/4_f.rds")
sp <- SCTransform(sp,assay = "Spatial")
SpatialFeaturePlot(sp,features = c("CTLA4","IDO1"))

sp <- readRDS("/data1/huanchangxiang/datasets/public_spatial/output/03LabelNew_second/4_h.rds")
sp <- SCTransform(sp,assay = "Spatial")
SpatialFeaturePlot(sp,features = c("CTLA4","IDO1"))


######## Figure 2F 代谢图 #######
library(ggplot2)
library(tidyverse)
library(Seurat)
library(clusterProfiler)
library(KEGGREST)
library(dplyr)
library(ComplexHeatmap)
sp <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sp_clusters.rds")

sp@meta.data$celltype_sample <- paste0(sp$seurat_clusters,"_",sp$sample)

# Hallmarkers获取通路
hallmarkers <- read.gmt("/data1/huanchangxiang/datasets/public_spatial/Genesets/h.all.v2024.1.Hs.symbols.gmt") %>% 
  split(.,.$term) %>% 
  lapply(., function(x){a <- x$gene})
# KEGG获取通路
kegg_pathways <- KEGGREST::keggList("pathway", 'hsa') %>%  
      as.data.frame() %>%  
      rownames_to_column() %>%  
      dplyr::rename("Pathway_ID" = 'rowname',"Description" = '.')

getKEGGPathwayGenes <- function(pathwayID) {
  # 获取KEGG通路信息
  gsInfo <- keggGet(pathwayID)[[1]]
  
  # 提取基因信息
  geneSetRaw <- sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
  
  # 生成基因集
  geneSet <- list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
  
  # 设置基因集名称
  names(geneSet) <- gsInfo$NAME
  
  return(geneSet)
}

Hypoxia <- hallmarkers[["HALLMARK_HYPOXIA"]]
Glycolysis  <- getKEGGPathwayGenes("hsa00010")
Lipid <- getKEGGPathwayGenes("hsa00565")
Pentose <- getKEGGPathwayGenes("hsa00030")
Oxidative <- getKEGGPathwayGenes("hsa00190")
Lactic <- list(c("LDHA","LDHB","MCT1","MCT4","PDK1","HK2","PKM2","GPR81"))
names(Lactic) <- "Lactic acid pathway"

pathway_list <- list(Hypoxia,Glycolysis[[1]],Lipid[[1]],Pentose[[1]],Lactic[[1]],Oxidative[[1]])
names(pathway_list) <- c("Hypoxia Pathway",
                         "Glycolysis pathway","Lipid Metabolism",
                         "Pentose Phosphate Pathway","Lactic","Oxidative")

for (i in 1:length(pathway_list)) {
      sp <- AddModuleScore(sp,
                           features = pathway_list[i],
                           name = names(pathway_list)[i])
}
# 计算分组平均值
meta <- sp@meta.data[, 11:16] %>%
  group_by(celltype_sample) %>%
  summarize(across(everything(), mean, na.rm = TRUE))

meta_avg <- meta %>% t()
colnames(meta_avg) <- meta_avg[1,]
meta_avg_filter <- meta_avg[2:6,]
meta_avg_filter <- as.data.frame(apply(meta_avg_filter,2,as.numeric))
rownames(meta_avg_filter) <- rownames(meta_avg)[2:6]

# 归一化
scaled_means <- apply(meta_avg_filter , 1, function(gene_expr) {
  min_val <- min(gene_expr)
  max_val <- max(gene_expr)
  if (max_val == min_val) {
    # 处理全零或全相同值的列
    return(rep(0, length(gene_expr)))
  } else {
    # 线性缩放公式
    return(2 * (gene_expr - min_val) / (max_val - min_val) - 1)
  }
}) %>%
  as.data.frame()

scaled_means$metabo_score <- scaled_means$Hypoxia.Pathway1 + 
                             scaled_means$Glycolysis.pathway1 + scaled_means$Lipid.Metabolism1 +
                             scaled_means$Pentose.Phosphate.Pathway1 + scaled_means$Lactic1

# 绘制山脊图
scaled_means_ord <- scaled_means %>% select(metabo_score)
scaled_means_ord$celltype_sample <- rownames(scaled_means_ord)

# 假设 scaled_means_ord 是你的数据框，metabo_score 是其中的一列
metabo_scores <- scaled_means_ord$metabo_score
scaled_means_ord <- scaled_means_ord[order(metabo_scores, decreasing = TRUE), ]
scaled_means_ord$celltype_sample <- factor(scaled_means_ord$celltype_sample,
                                           levels = rownames(scaled_means_ord))
# 绘制山脊图
p1 <- ggplot(scaled_means_ord, aes(x = celltype_sample, y = metabo_score, fill = metabo_score)) +
  geom_col() +
  scale_fill_gradient2(
    low = "#2FBDFF",       # 低值颜色（如最小值）
    mid = "white",     # 中间值颜色（如中位数或均值）
    high = "#FF7582",     # 高值颜色（如最大值）
    midpoint = median(scaled_means_ord$metabo_score)  # 自定义中间值位置
  ) +
  theme_minimal()
ggsave(plot = p1,filename = "Figure2_heatmap_注释柱状图.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure2/",
       width = 8,height = 3.4)
# 绘制下面热图
dist_matrix <- dist(scaled_means,method = "euclidean")
hc <- hclust(dist_matrix,method = "complete")
scaled_means$sample <- gsub(".*_","",rownames(scaled_means))
row_anno <- rowAnnotation(Sample=scaled_means$sample)

p <- ComplexHeatmap::Heatmap(scaled_means[,1:6],
                             row_order = rownames(scaled_means_ord),
                             cluster_columns = T,cluster_rows = F,
                             right_annotation = row_anno)
p

ggsave(plot = p,filename = "Figure2_metabo_heatmap.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure2/",
       width = 6.4,height = 10)
######## Figure 2G 不同肿瘤分布图 #######
library(tidyverse)
library(Seurat)
library(magrittr)
scaled_means_ord <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/scaled_means_ord.rds")
Hyper <- rownames(scaled_means_ord)[scaled_means_ord$metabo_score > 0]
Hypo  <- rownames(scaled_means_ord)[scaled_means_ord$metabo_score < -1.08]

sp_sub_celltype <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sp_clusters.rds")
sp_sub_celltype$celltpe_sample <- paste0(sp_sub_celltype$seurat_clusters,"_",sp_sub_celltype$sample)
sp_sub_celltype$barcode <- rownames(sp_sub_celltype@meta.data)

sample <- c("a_celltype.rds","d_celltype.rds","f_celltype.rds","h_celltype.rds")
sp_list <- sapply(sample, function(x){
      sp <- readRDS(paste0("/data1/huanchangxiang/datasets/李莹雪论文/rds/",x))
      #rownames(sp@meta.data) <- paste0(rownames(sp@meta.data),"_",gsub("_celltype.rds","",x))
      return(sp)
})
sp <- rbind(sp_list[[1]]@meta.data,sp_list[[2]]@meta.data,sp_list[[3]]@meta.data,sp_list[[4]]@meta.data)
sp <- sp[,c(4:14)] %>% as.data.frame()
a <- apply(sp, 1, FUN = function(row){
  
  num6 <- head(sort(row,decreasing = T),6)
  row[!names(row) %in% names(num6)] <- 0
  for (i in 1:length(row)) {
    row[i] <- row[i]/sum(row)
  }
  return(row)
})
sp <- t(a) %>% as.data.frame()
sp$barcode <- rownames(sp)

sp <- merge(sp,sp_sub_celltype@meta.data[,c("barcode","celltpe_sample")],by = "barcode")


sp <- sp %>%  mutate(metabo_group = case_when(
                  celltpe_sample %in% Hyper ~ "Hypermetabolic Tumor",
                  celltpe_sample %in% Hypo  ~ "Hypometabolic Tumor",
                  TRUE ~ as.character(celltpe_sample)
                ))






### 制作meta
library(dplyr)
library(tidyr)
for (i in colnames(sp)[2:12]) {
  meta <- sp %>%
    filter(metabo_group %in% c("Hypermetabolic Tumor", "Hypometabolic Tumor")) %>%
    select(2:12, metabo_group) %>%
    mutate(across(
      1:11,  # 假设前11列为数值型列
      ~ {
        min_val <- min(., na.rm = TRUE)
        max_val <- max(., na.rm = TRUE)
        if (max_val == min_val) {
          rep(0, length(.))
        } else {
          (. - min_val) / (max_val - min_val)  # 此处错误：多余的右括号
        }
      }
    )) %>%
    pivot_longer(
      cols = -metabo_group,
      names_to = "celltype",
      values_to = "value"
    ) %>%
    filter(celltype == i)
  # 计算 IQR 范围
  q <- quantile(meta$value, c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  lower <- q[1] - 1.5 * iqr
  upper <- q[2] + 1.5 * iqr
  
  # 绘图并限制坐标轴
  library(ggpubr)
  
 p1 <-  ggplot(meta, aes(x = celltype, y = value, fill = metabo_group)) +
    geom_boxplot(outlier.shape = 16) +
    stat_compare_means(
      aes(group = metabo_group),  # 按 metabo_group 分组比较
      method = "wilcox.test",     # 使用 Wilcoxon 检验
      label = "p.signif",         # 显示星号标记（* **, **, *）
      label.y = 1,  # 标签位置调整
      hide.ns = F,             # 隐藏不显著的结果
      na.rm = T,size=10
    ) +
    theme(
      axis.line = element_line(size = 2, colour = "black"),
      axis.text = element_text(size = 20),
      axis.title = element_text(size= 20)
    ) +
    labs(x="",y="Cell proportion (Z-score)")
 ggsave(plot = p1,filename = paste0("Figure2_metabo_boxplot_",i,".pdf"),
        path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure2/",
        width = 9.7,height = 7.7)
}

