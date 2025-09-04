library(Seurat)
library(tidyverse)
library(magrittr)
library(sctransform)
library(harmony)
sc <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/merge/merge_raw.rds")
Idents(sc) <- sc$celltype

Fib <- subset(sc,idents = "Fibroblasts")

Fib <-  SCTransform(Fib)
Fib <- RunPCA(Fib,npcs = 30)
Fib <- RunHarmony(Fib,"patient")
Fib <- RunUMAP(Fib,dims = 1:10) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.2)

DimPlot(Fib,label = T)

DotPlot(Fib,features = c("ALB","FAP"))

marker <- FindAllMarkers(Fib,only.pos = T,logfc.threshold = 1,min.pct = 0.45)

annotation <- function(x){
  # annotation 1
  list_genes <- list(
    iCAFs = c("DPT",'CXCL9','CXCL10',"CXCL11","EFEMP1","PDPN"), 
    vCAFs = c("ACTA2","NOTCH3","DSTN","ADIRF"),
    ECM = c("COL6A3","VIM","COL5A1","FAP","FN1","POSTN","MMP11"),
    lipid_Process = c("APOC3","APOC1","APOA2","APOC2","FABP1"),
    apCAFs = c("HLA-DPB1","CXCL12","CCL21","HLA-DRA","IGFBP3"),
    dCAFS = c("MKI67","STMN1","TOP2A","CCNB1","TYMS","HMGB2","PCNA","DEK"),
    CAF_HSC = c("MRC1","SLC9A9","PTPRB","STAB2","SEMA6A","GUCY1B1","CCL2","CX3CL1",
                "MMP2","MMP9","TGFB1","PDGFRB","THY1"),
    CAF_VSMC = c("RGS5", "NDUFA4L2", "MYH11","CNN1"),
    CAF_Port = c("PDGFRA", "MMP23B", "COL1A1", "PRELP"),
    quiescent_HSC = c("DCN","NGFR","LRAT","RBP1","ECM1")
    
  )
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
             cluster.idents = TRUE,col.min = 0,col.max = 1)+
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

Fib_filter <- subset(Fib,idents=c("4","5","6"),invert=T)

Fib_filter_reclust <-  SCTransform(Fib_filter) %>%
        RunPCA(.,npcs = 30) %>%
        RunHarmony(.,"patient") %>%
        RunUMAP(.,dims = 1:10) %>%
        FindNeighbors(dims = 1:10) %>%
        FindClusters(resolution = 0.2)

Fib_filter_reclust_filter <- subset(Fib_filter_reclust,idents="6",invert=T)

Idents(Fib_filter_reclust_filter) <- Fib_filter_reclust_filter$SCT_snn_res.0.2
annotation(Fib_filter_reclust_filter)
saveRDS(Fib_filter_reclust_filter,file = "/data1/huanchangxiang/datasets/李莹雪论文/rds/Fib_second_clusters.rds")

########## Figure3A ##########
if (!dir.exists("/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure3/")) {
  dir.create("/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure3/")
}
Fib_filter_reclust_filter <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/Fib_second_clusters.rds")
DimPlot(Fib_filter_reclust_filter,label = T)

########## Figure3B ##########
VlnPlot(Fib_filter_reclust_filter,features = c("ACTA2","POSTN","COL1A1","FAP"))


########## Figure3C ##########
DotPlot(Fib_filter_reclust_filter,
        features = c("CRABP2","PLA2G2A",
                     "HLA-DRA","CD74",
                     "APOD","S100A4","FBLN1","CXCL12","PDGFB",
                     "COL10A1","POSTN",
                     "ACTA2","MCAM"))

########## Figure3D ##########
library(clusterProfiler)
library(magrittr)
Fib_filter_reclust_filter <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/Fib_second_clusters.rds")
hallmarkers <- read.gmt("/data1/huanchangxiang/datasets/public_spatial/Genesets/h.all.v2024.1.Hs.symbols.gmt") %>% 
  split(.,.$term) %>% 
  lapply(., function(x){a <- x$gene})

for (i in 1:length(hallmarkers)) {
  Fib_filter_reclust_filter <- AddModuleScore(Fib_filter_reclust_filter,
                                              features = hallmarkers[i],
                                              name = names(hallmarkers)[i])
}

meta <- Fib_filter_reclust_filter@meta.data[,13:62] %>% as.data.frame()
colnames(meta) <- gsub("HALLMARK_","",colnames(meta))
colnames(meta) <- gsub("1","",colnames(meta))

scaled_means <- lapply(meta, function(gene_expr) {
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
  as.data.frame() %>%
  setNames(colnames(meta)) %>%  # 保持列名
  `rownames<-`(rownames(meta))  # 保持行名



library(ComplexHeatmap)
library(circlize)
group <- Fib_filter_reclust_filter@meta.data[,12:13] %>% select("SCT_snn_res.0.2")
colnames(group) <- "Cluster"
column_order <- order(group$SCT_snn_res.0.2)
Heatmap(t(scaled_means),cluster_rows = F,cluster_columns = F,show_column_names = F,
        top_annotation = HeatmapAnnotation(df=group),
        col = colorRamp2(c(-0.8, 0, 0.8), c("blue", "white", "red")),
        column_order = column_order,
        column_split = group$Cluster)


###################  Figure3F #########################
#RO/e指数绘图
Idents(data) = data$NorT
data=RenameIdents(data,"T"="Tumor","N"="Adjacent")
data$Tissue=data@active.ident
meta=data@meta.data
meta$Annotation=meta$seurat_clusters

if(TRUE){  
  divMatrix <- function(m1, m2){
    dim_m1 <- dim(m1)
    dim_m2 <- dim(m2)
    if( sum(dim_m1 == dim_m2) == 2 ){
      div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
      row.names(div.result) <- row.names(m1)
      colnames(div.result) <- colnames(m1)
      for(i in 1:dim_m1[1]){
        for(j in 1:dim_m1[2]){
          div.result[i,j] <- m1[i,j] / m2[i,j]
        }
      }
      return(div.result)
    }
    else{
      warning("The dimensions of m1 and m2 are different")
    }
  }
  ROIE <- function(crosstab){
    rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
    rowsum.matrix[,1] <- rowSums(crosstab)
    colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
    colsum.matrix[1,] <- colSums(crosstab)
    allsum <- sum(crosstab)
    roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
    row.names(roie) <- row.names(crosstab)
    colnames(roie) <- colnames(crosstab)
    return(roie)
  }
}
plot_df = meta %>%   
  filter(Tissue %in% c("Adjacent","Tumor"))
plot_df$Annotation = as.character(plot_df$Annotation)
plot_df$Tissue = as.character(plot_df$Tissue)
summary <- table(plot_df[,c('Annotation','Tissue')])
roe <- as.data.frame(ROIE(summary))
roe <- roe[,c("Adjacent","Tumor")]

## 设置颜色
require(ComplexHeatmap)
require(circlize)
col_fun <- colorRamp2(c(0, 1, 1.5, 2),c("#fffde7", "#ffe0b2","#ff9800", "#e65100"))
pdf("/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure3/Figure3F_Roe.pdf",
    width = 9,height = 9)
p2 <- Heatmap(  roe,  
                col = col_fun,  
                cluster_rows = T,  
                cluster_columns = F,  
                clustering_distance_rows = "euclidean",  
                clustering_method_rows = "ward.D",  
                column_names_gp = grid::gpar(fontsize = 6),  
                row_names_gp = grid::gpar(fontsize = 6),  
                cell_fun = function(j, i, x, y, width, height, fill) {    
                  grid.text(sprintf("%.1f", roe[i, j]), x, y, gp = gpar(fontsize = 6))  },  
                width = ncol(roe) * unit(0.3, "inch"),  
                height = nrow(roe) * unit(0.15, "inch"),  
                name = "Ro/e"
)
draw(p2)
dev.off()

#####################  Remodeling of the fibroblast cell  during  cancer ##################### 
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
list_genes <- c("DPT",'CXCL9','CXCL10',"CXCL11","EFEMP1","PDPN","ACTA2","NOTCH3","DSTN","ADIRF",
                "COL6A3","VIM","COL5A1","FAP","FN1","POSTN","MMP11","APOC3","APOC1","APOA2","APOC2","FABP1",
                "DCN","NGFR","LRAT","RBP1","ECM1")
mat <- data.frame(data@assays$RNA@counts[list_genes,]) %>% 
  t() %>% 
  as.data.frame() %>%
  mutate(NorT = Fib$NorT, celltype_second = data$seurat_clusters) %>%
  group_by(celltype_second, NorT) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  mutate(cell_NorT = paste0(celltype_second, "_",NorT)) %>%
  column_to_rownames(var = "cell_NorT")
anno <- mat[,c("celltype_second","NorT")]
mat_filter <- mat %>% select(-celltype_second,-NorT) %>% scale() %>% t()

###  设置颜色
col_fun <- colorRamp2(breaks = c(-2,0,2),colors = c("#1E90FF", "white", "#DC143C"))
NorT_colors <- c("N" = "skyblue", "T" =  "#DC143C")
celltype_colors <- c(
  "0" = "#1F78B4","1" = "#33A02C","2" = "#FB9A99", "3" = "#E31A1C","4" = "#CAB2D6",
  "5" = "#B15928" , "6" = "#A6CEE3"
)
### 设置上面注释
top_annotation <- HeatmapAnnotation(df = anno,col = list(
  celltype_second=celltype_colors,NorT=NorT_colors))

pdf("/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure3/Figure3F_function_clusters.pdf",
    width = 6,height = 10)
Heatmap(mat_filter,cluster_rows = F,cluster_columns = F,
        top_annotation = top_annotation,column_split = anno$celltype_second,
        width = ncol(mat_filter)*unit(5, "mm"), 
        height = nrow(mat_filter)*unit(5, "mm"),
        show_column_names = F,
        col = col_fun)
dev.off()


#####################  Pathway activity score ##################### 
library(progeny)
library(ggplot2)
library(dplyr)
library(Seurat)
library(tidyr)
library(tidyverse)
library(ggh4x)
data <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/Fib_second_clusters.rds")
model <- progeny::model_human_full

model_100 <- model %>% group_by(pathway) %>% slice_min(order_by = p.value,n=100)

Idents(data) <- data$seurat_clusters

CellsClusters <- data.frame(Cell=names(Idents(data)),
                            Celltype=as.character(Idents(data)),
                            stringsAsFactors = F)
Fib <- progeny(data,scale = F,organism = "Human",top = 200,perm = 1,return_assay = T)
Fib <- Seurat::ScaleData(Fib,assay = "progeny")
progen_score_df <- as.data.frame(t(GetAssayData(Fib,slot = "scale.data",assay = "progeny"))) %>%
  rownames_to_column("Cell")     %>%
  gather(Pathway,Activity,-Cell) %>%
  inner_join(.,CellsClusters)    %>%
  group_by(Pathway,Celltype)     %>%
  summarise(avg = mean(Activity), std =sd(Activity)) %>%
  filter(Pathway %in% c("NFkB","Hypoxia","JAK-STAT","MAPK","TGFb","WNT","VEGF","PI3K"))

# 确保Celltype为因子并统一水平（避免分面轴不一致）
progen_score_df$Celltype <- factor(progen_score_df$Celltype)
progen_score_df$group_col <- as.character(progen_score_df$Celltype)  # 确保为字符型
col <- c("1" = "#8B0000","0"="#006400","2"="#DAA520","3"="#4B0082","4"="#FF6347","5"="skyblue")
pathway_colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B")
# 绘制分面散点图
p1 <- ggplot(progen_score_df, aes(x = avg, y = Celltype)) +
  geom_segment(
    aes(x = -1, xend = 1, y = Celltype, yend = Celltype, color = group_col),
    linewidth = 1.2
  ) +
  geom_point(aes(color = group_col), size = 5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap2(~ Pathway, ncol = 8, 
              strip = strip_themed(
                background_x = elem_list_rect(fill = pathway_colors),
                text_x = elem_list_text(color = "white")
              )) +
  labs(x = "Pathway Activity Score", y = "Cell Type") +
  theme_minimal() +
  scale_color_manual(values = col) +
  theme(axis.title = element_text(size=10),
        axis.text  = element_text(size=15),
        axis.text.x  = element_text(angle = 270),
        axis.title.x  = element_text(size=15),
        axis.title.y  = element_text(size=15),
        legend.position = "none")

print(p1)
ggsave(plot = p1,filename = "Figure3_pathway_activity_score.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure3/",
       width = 12,height = 5.5)





