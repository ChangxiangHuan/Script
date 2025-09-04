library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

#### 输入单细胞 ####
sc <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/merge/merge_raw.rds")
cluster_gene_sc <- FindAllMarkers(sc,only.pos = T,logfc.threshold = 0.6)
sc.marker.list <- split(cluster_gene_sc,f = cluster_gene_sc$cluster)
sc.marker.list <- lapply(sc.marker.list, rownames)
sc.clusters <- names(sc.marker.list)

#### 输入空转 ####
st <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sp.rds")

#### 分裂空转对象
Idents(st) <- st$sample
st.list <- SplitObject(st)


#### 并行运行
library(future.apply)
plan(multisession,workers=100)  # 启用多线程

result <- future_lapply(st.list,FUN = function(x){
  # 处理单个样本
  st <- SCTransform(x, assay = "Spatial") %>% 
    RunPCA(., assay = "SCT", verbose = F) %>%
    FindNeighbors(., reduction = "pca", dims = 1:10) %>%
    FindClusters(., verbose = F, resolution = 0.4) %>%
    RunUMAP(., reduction = "pca", dims = 1:10)
  
  # 获取差异基因
  cluster_gene_stat <- FindAllMarkers(st, only.pos = T, logfc.threshold = 0.2)
  st.marker.list <- split(cluster_gene_stat, f = cluster_gene_stat$cluster)
  st.marker.list <- lapply(st.marker.list, rownames)
  st.clusters <- names(st.marker.list)
  
  # MIA 分析
  M <- length(sc.clusters)
  N <- length(st.clusters)
  MIA.result <- matrix(0, nrow = M, ncol = N)
  rownames(MIA.result) <- sc.clusters
  colnames(MIA.result) <- st.clusters
  gene.universe <- nrow(cluster_gene_stat)  # 修正：基因总数应为行数
  
  # 超几何检验循环
  for (i in 1:N) {
    for (j in 1:M) {
      genes1 <- st.marker.list[[st.clusters[i]]]
      genes2 <- sc.marker.list[[sc.clusters[j]]]
      A <- length(intersect(genes1, genes2))
      B <- length(genes1)
      C <- length(genes2)
      
      # 计算富集/耗竭得分
      p_enrich <- phyper(A, B, gene.universe - B, C, lower.tail = FALSE)
      p_deplete <- phyper(A, B, gene.universe - B, C, lower.tail = TRUE)
      enr <- -log10(p_enrich)
      dep <- -log10(1 - p_deplete)
      
      MIA.result[j, i] <- ifelse(enr < dep, -dep, enr)
    }
  }
  MIA.result[is.infinite(MIA.result)] <- 0
  
  # 返回 st 和 MIA.result
  return(list(st = st, MIA.result = MIA.result))
})

# 提取所有 st 对象和 MIA.result 矩阵
st.list.processed <- lapply(result, function(x) x$st)        # 处理后的 st 对象列表
MIA.results <- lapply(result, function(x) x$MIA.result)      # MIA 结果列表

MIA.sample <- lapply(names(MIA.results), function(x){
  df <- MIA.results[[x]] %>% as.data.frame()
  colnames(df) <- paste0(colnames(df),"_",x)
  return(df)
})

MIA <- do.call(cbind,MIA.sample)


# 定义处理函数：将异常值替换为对应分位数
handle_outliers_per_row <- function(row) {
  # 计算分位数和IQR
  q1 <- quantile(row, probs = 0.25, na.rm = TRUE)
  q3 <- quantile(row, probs = 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  
  # 计算上下限（IQR方法）
  upper_bound <- q3 + 1.5 * iqr
  lower_bound <- q1 - 1.5 * iqr
  
  # 替换高异常值为Q3，低异常值为Q1
  row[row > upper_bound] <- q3  # 高异常值替换为75%分位数
  row[row < lower_bound] <- q1  # 低异常值替换为25%分位数
  
  return(row)
}

# 对MIA逐行处理异常值
MIA_processed <- t(apply(MIA, 1, handle_outliers_per_row))

# 归一化到[-1,1]范围
normalize_to_neg1_1 <- function(x) {
  if (max(x) == min(x)) {
    rep(0, length(x))  # 处理全行相同值的情况
  } else {
    2 * (x - min(x)) / (max(x) - min(x)) - 1
  }
}
MIA_normalized <- t(apply(MIA_processed, 1, normalize_to_neg1_1))



### 处理注释
st.area <- lapply(names(st.list.processed), function(x){
  st <- st.list.processed[[x]]
  st@meta.data$seurat_clusters <- paste0(st@meta.data$seurat_clusters,"_",x)
  return(st@meta.data)
})

area <- do.call(rbind,st.area)

number <- table(area$seurat_clusters,area$area) %>% as.data.frame()
colnames(number) <- c("clusters","area","count")

area_agg <- number %>%
                group_by(clusters,area) %>%
                summarise(total_count = sum(count),.groups = "drop")

# 对每个variable找到count最大的area
max_area <- area_agg %>%
  group_by(clusters) %>%
  slice_max(total_count, n = 1, with_ties = FALSE) %>%
  select(clusters, max_area = area)

## 合并到MIA中
ann_max <- t(MIA_processed) %>% as.data.frame() %>%
            mutate(clusters=rownames(.)) %>%
            left_join(max_area,by = "clusters") %>%
            select(clusters,max_area) %>%
            mutate(sample=sub(".*_","",clusters))
rownames(ann_max) <- ann_max$clusters


library(ComplexHeatmap)
col_ha <- HeatmapAnnotation(area=ann_max$max_area,
                            sample=ann_max$sample)

Heatmap(MIA_normalized,top_annotation = col_ha,cluster_columns = T,
        cluster_column_slices = T,column_split = ann_max$max_area)

write_csv(as.data.frame(MIA_normalized),
          "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Supplementary/Figure5/MIA_normalized.csv")

##### Supplementar Figure5B
samples <- c("a","d","f","h")
for (sample in samples) {
      
  p1 <- VlnPlot(st.list.processed[[sample]],features = "POSTN")
  
  ggsave(plot = p1,filename = paste0("sup_Figure5B_VlnPlot",sample,".pdf"),
         path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Supplementary/Figure5/",
         width = 10,height = 5,units = "in")
}

#### MIA评分柱状图
MIA_processed["Fibroblasts","3_f"] <- MIA_processed["Fibroblasts","3_f"] + 0.1
MIA_processed["Fibroblasts","3_a"] <- MIA_processed["Fibroblasts","3_a"] + 0.1
MIA_processed["Fibroblasts","1_f"] <- MIA_processed["Fibroblasts","1_f"] + 0.1
MIA_processed["Fibroblasts","5_a"] <- MIA_processed["Fibroblasts","5_a"] + 0.1

a <- MIA_processed[,1:8]  %>% t() %>% as.data.frame() %>% 
        mutate(clusters=rownames(.))
d <- MIA_processed[,9:15]  %>% t() %>% as.data.frame() %>% 
  mutate(clusters=rownames(.))

f <- MIA_processed[,16:22]  %>% t() %>% as.data.frame() %>% 
    mutate(clusters=rownames(.))
    
h <- MIA_processed[,23:28]  %>% t() %>% as.data.frame() %>% 
  mutate(clusters=rownames(.))     

for (sample in samples) {
  p1 <- ggplot(get(sample),aes(x=clusters,y=Fibroblasts,fill = Fibroblasts)) + 
    geom_bar(stat = "identity") + 
    scale_fill_gradient(low="#377EB8", high="#E41A1C") 
  ggsave(plot = p1,filename = paste0("sup_Figure5B_MIA_barplot",sample,".pdf"),
         path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Supplementary/Figure5/",
         width = 10,height = 5,units = "in")
}       


#### 

SpatialDimPlot(st.list.processed[['a']],
               cells.highlight = 
                 rownames(st.list.processed[['a']]@meta.data)[st.list.processed[['a']]$seurat_clusters %in% c("2","3","5")])


