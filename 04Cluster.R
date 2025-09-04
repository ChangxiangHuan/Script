library(Seurat)
library(patchwork)
library(magrittr)
library(harmony)
library(ggplot2)
library(dplyr)
library(Startrac)
data.dir <- "/data1/huanchangxiang/datasets/public_single/rds/filter/"
path = "/data1/huanchangxiang/datasets/public_single/rds/filter"
files = list.files(path)
options(future.globals.maxSize = 1024 * 1024 * 1024*10) # 设置为10GiB
pattern <- "^filter_[ABDadfh]"
filtered_files <- grep(pattern, files, value = TRUE)
normal_sample <- c('A3','A5','A7','A9','A11','A14','A17','A20','B1','B7','B11','B17',
                   'B20','B21','B24','B27','B30','B36','B41','B44','B47','B51','B54',
                   'D1','D3','D5','D7','D9','D12')
sce_list <- list()
for ( i in filtered_files) {
  extracted_content <- sub("filter_(.*)\\.rds", "\\1", i)
  sce_list[[extracted_content]] <- readRDS(paste0(data.dir,i))
  if (extracted_content %in% normal_sample) {
    sce_list[[extracted_content]]$NorT <- "N"
  } else {
    sce_list[[extracted_content]]$NorT <- "T"
  }
}

### merge ###
sce_merge <- merge(sce_list[[1]],sce_list[-1])
rm(sce_list)
sce_merge@meta.data <- sce_merge@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","NorT")]
### annotation patient ###
sce_merge@meta.data <- sce_merge@meta.data %>%
  mutate(patient = case_when(
    orig.ident %in% c("A1") ~ "Patien1",
    orig.ident %in% c("A2") ~ "Patien2",
    orig.ident %in% c("A3","A4") ~ "Patien3",
    orig.ident %in% c("A5","A6") ~ "Patien4",
    orig.ident %in% c("A7","A8") ~ "Patien5",
    orig.ident %in% c("A9","A10") ~ "Patien6",
    orig.ident %in% c("A11","A12","A13") ~ "Patien7",
    orig.ident %in% c("A14","A15","A16") ~ "Patien8",
    orig.ident %in% c("A17","A18") ~ "Patien9",
    orig.ident %in% c("A19","A20","A21") ~ "Patien10",
    orig.ident %in% c("B1","B2","B3","B4","B5","B6") ~ "Patien11",
    orig.ident %in% c("B7","B8","B9","B10") ~ "Patien12",
    orig.ident %in% c("B11","B12","B13","B14","B15","B16") ~ "Patien13",
    orig.ident %in% c("B17","B18","B19") ~ "Patien14",
    orig.ident %in% c("B20") ~ "Patien15",
    orig.ident %in% c("B21","B22","B23") ~ "Patien16",
    orig.ident %in% c("B24","B25","B26") ~ "Patien17",
    orig.ident %in% c("B27","B28","B29") ~ "Patien18",
    orig.ident %in% c("B30","B31","B32","B33","B34","B35") ~ "Patien19",
    orig.ident %in% c("B36","B37","B38","B39","B40") ~ "Patien20",
    orig.ident %in% c("B41","B42","B43") ~ "Patien21",
    orig.ident %in% c("B44","B45","B46") ~ "Patien22",
    orig.ident %in% c("B47","B48","B49","B50") ~ "Patien23",
    orig.ident %in% c("B51","B52","B53") ~ "Patien24",
    orig.ident %in% c("B54","B55","B56","B57","B58") ~ "Patien25",
    orig.ident %in% c("P11") ~ "Patien26",
    orig.ident %in% c("P12") ~ "Patien27",
    orig.ident %in% c("P15") ~ "Patien28",
    orig.ident %in% c("P2") ~ "Patien29",
    orig.ident %in% c("P21") ~ "Patien30",
    orig.ident %in% c("P22") ~ "Patien31",
    orig.ident %in% c("P22") ~ "Patien31",
    orig.ident %in% c("a") ~ "Patien32",
    orig.ident %in% c("d") ~ "Patien33",
    orig.ident %in% c("f") ~ "Patien34",
    orig.ident %in% c("C4","C6") ~ "Patien35",
    orig.ident %in% c("h") ~ "Patien36",
    #orig.ident %in% c("C7") ~ "Patien37",
    # orig.ident %in% c("C8","C29") ~ "Patien38",orig.ident %in% c("C9") ~ "Patien39",
    # orig.ident %in% c("C11","C13") ~ "Patien40",orig.ident %in% c("C12","C30","C31") ~ "Patien41",
    # orig.ident %in% c("C14") ~ "Patien42",orig.ident %in% c("C17") ~ "Patien43",
    # orig.ident %in% c("C19") ~ "Patien44",orig.ident %in% c("C21","C22","C23") ~ "Patien45",
    # orig.ident %in% c("C24") ~ "Patien46",orig.ident %in% c("C26") ~ "Patien47",
    # orig.ident %in% c("C32") ~ "Patien48",orig.ident %in% c("C33") ~ "Patien49",
    # orig.ident %in% c("C35") ~ "Patien50",orig.ident %in% c("C36") ~ "Patien51",
    # orig.ident %in% c("C39") ~ "Patien52",orig.ident %in% c("C41") ~ "Patien53",
    # orig.ident %in% c("C42") ~ "Patien54",orig.ident %in% c("C43") ~ "Patien55",
    # orig.ident %in% c("C46") ~ "Patien56",
    TRUE ~ as.character(orig.ident)  # 使用原始列值作为默认情况
  ))

### SCTtransfrom ###
sce_merge <- sce_merge %>% SCTransform(variable.features.n = 3000,
                                       vars.to.regress = "patient") %>%  
  RunPCA(.,npcs = 40,assay = "SCT") %>% RunHarmony(group.by.vars = "patient", 
                                                   max.iter.harmony = 10,assay.use = "SCT") %>%
  RunUMAP( reduction="harmony", dims=1:30,min.dist = 1) %>%           
  FindNeighbors(reduction="harmony", dims=1:30) %>%             
  FindClusters(resolution=0.6,algorithm = 1)

saveRDS(sce_merge,file = "/data1/huanchangxiang/datasets/public_single/rds/merge/merge_raw.rds")

sce_merge <- readRDS("/data1/huanchangxiang/datasets/public_single/rds/merge/merge_raw.rds")

# 注释文件
annotation <- function(x){
  # annotation 1
  list_genes <- list(
    Hepatocytes=c("APOA2","ALB","TTR","PECAM1"),
    Epithelial_cells = c("EPCAM","KRT19","KRT18","DEFB1","CTSK"), 
    Fibroblasts = c("ACTA2","RGS5",'COL1A1',"COL1A2"),
    Endothelial_cells=c('VWF',"CD34"),
    Myeloid_cells=c("CD14","CD68","CD163","S100A8"),
    Dendritic_cells=c("IL3RA","IRF7","IRF8","GZMB"),
    LSEC = c("DNASE1L3","CLEC1B","CLEC4G","CLEC4M"),
    Mast_cells=c('TPSAB1'),
    #Neutrophils=c("CSF3R","S100A8"),
    T_NK=c('CD3D',"CD3E","CD4","CD8A","CD3G","NKG7","IL7R"),
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
Idents(sce_merge) <- sce_merge$seurat_clusters
DimPlot(sce_merge,label=T)
annotation(sce_merge)

sce_merge@meta.data <- sce_merge@meta.data %>%
  mutate(celltype = case_when(
    seurat_clusters %in% c("40","28") ~ "B_cells",
    seurat_clusters %in% c("41") ~ "Dendritic_cells",
    seurat_clusters %in% c("39","33") ~ "Plasma_cells",
    seurat_clusters %in% c("30","29","37","24","1","23",
                           "19","16","18","15","14","43","38",
                           "46","51","26","3","31","12") ~ "Hepatocyte",
    seurat_clusters %in% c("45") ~ "Mast_cells",
    seurat_clusters %in% c("22","2","21","48","44","4","50","49","42",
                           "6","35","47","36","32","0","13","25","17","34") ~ "T_NK",
    seurat_clusters %in% c("20","7","8","27","11") ~ "Myeloid_cells",
    seurat_clusters %in% c("10") ~ "Fibroblasts",
    seurat_clusters %in% c("9") ~ "LSEC",
    seurat_clusters %in% c("5") ~ "Endothelial_cells",
    TRUE ~ as.character(seurat_clusters)  # 使用原始列值作为默认情况
  ))

sce_merge@meta.data[sce_merge@meta.data$celltype %in% "Cancer_cells" & 
                      sce_merge@meta.data$NorT %in% "T","celltype"] <- "Hepatocytes"
Idents(sce_merge) <- sce_merge$celltype
DimPlot(sce_merge,label = T) +




saveRDS(sce_merge,file = "/data1/huanchangxiang/datasets/李莹雪论文/rds/merge/merge_raw.rds")


#### DimPlot ####
data <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/merge/merge_raw.rds")
DimPlot(data,pt.size = 0.8,
        cols = c("B_cells" = "#FFDAB9","Cancer_cells"="#FFA07A","Dendritic_cells"="#FF4500", 
                 "Endothelial_cells"="#FFD700","Fibroblasts"="#ADFF2F","Hepatocyte"="#FF6347",
                 "LSEC"="#32CD32","Mast_cells"="#6495ED","Myeloid_cells"="#1E90FF",
                 "Plasma_cells"="#4B0082","T_NK"="#FFC0CB"))

"#FFC0CB", "#87CEFA", "#7FFFD4", "#F08080"
################### 细胞比例图 ###################

Idents(data) <- data$celltype
data$patient_NorT <- paste0(data$patient,data$NorT)
sample_table <- as.data.frame(table(data@meta.data$patient_NorT,data@meta.data$celltype))
names(sample_table ) <- c("Samples","celltype","CellNumber")
#所有患者
Cellratio <- data.frame(prop.table(table(Idents(data), data$patient_NorT), margin = 2))#计算各组样本不同细胞群比例
colnames(Cellratio) <- c("celltype","patient","Freq")
meta <- split(Cellratio,Cellratio$celltype)
orde <- meta$T_NK%>% arrange(Freq)
sample_table$Samples <- factor(sample_table$Samples,levels = orde$patient)
colors <- c("#FFDAB9", "#FFA07A",  "#FF4500", "#FFD700",
            "#ADFF2F", "#FF6347","#32CD32", "#6495ED", "#1E90FF", "#4B0082",
            "#FFC0CB", "#87CEFA", "#7FFFD4", "#F08080")

p1 <-ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
  scale_fill_manual(values = colors) +
  geom_bar(position="fill",width = 0.7,size = 0.5)+  
  theme(panel.grid = element_blank(),panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)  )+
  labs(y="Percentage")+RotatedAxis()
p1

################### marker基因展示 ###################
allmarkers <- FindAllMarkers(data,only.pos = T,logfc.threshold = 1.2,min.pct = 0.5)
allmarkers$pct <- allmarkers$pct.1 - allmarkers$pct.2
top3 <- allmarkers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = pct, with_ties = FALSE) %>%
  slice_max(n = 3, order_by = avg_log2FC, with_ties = FALSE) 
list_genes <- list(
  B_cells=c('CD79A','MS4A1',"HLA-DRA"),
  Plasma_cells=c("IGHG4","IGHG1","IGHG3"),
  LSEC = c("DNASE1L3","CLEC1B","LIFR"),
  Fibroblasts = c("ACTA2","RGS5","COL1A1"),
  Dendritic_cells=c("CD1C","JCHAIN","IGKC"),
  Endothelial_cells=c('TM4SF1',"CD34"),
  Myeloid_cells=c("CD68","CD163"),
  T_NK=c('CD3D',"CD3E","NKG7"),
  Hepatocyte=c("APOA2","ALB","TTR"),
  Cancer_cells=c("TF","FGA","FGB"),
  Mast_cells=c("TPSB2","TPSAB1","CPA3")
)

# 设置细胞类型顺序
Idents(data) <- data$celltype
new_order <- c("B_cells","Plasma_cells","LSEC","Fibroblasts","Dendritic_cells",
               "Endothelial_cells","Myeloid_cells","T_NK","Hepatocyte","Cancer_cells","Mast_cells")  
Idents(data) <- factor(Idents(data), levels = new_order)

# 绘制DotPlot
p2 <- DotPlot(data, 
              features = list_genes,
              col.min = -2, col.max = 2,
              dot.scale = 10) +
  scale_color_gradientn(colours= c("#5E4FA2FF","#3288BDFF","#66C2A5FF","#ABDDA4FF","#E6F598FF","#FFFFBFFF",
                                   "#FEE08BFF","#FDAE61FF","#F46D43FF","#D53E4FFF","#9E0142FF")) + 
  theme(
    panel.border = element_rect(color = "gray50", size = 1),
    panel.spacing = unit(2, "mm"),
    panel.grid.major = element_line(color = "gray70", size = 0.5, linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray70", size = 0.5, linetype = "dashed"),
    axis.text.x = element_text(angle = 70, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) 
p2
ggsave(plot = p2,filename = "Figure1D_marker_Dotplot.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure1/",
       width = 13,height = 6.74)
##### 
FeaturePlot(data,features = "CLEC1B")
FeaturePlot(data,features = "TF")
FeaturePlot(data,features = "IGHG4")
FeaturePlot(data,features = "CD1C")



library(Seurat)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(harmony)
library(ComplexHeatmap)
setwd("/data1/huanchangxiang/datasets/spatial/hcc/")
options(future.globals.maxSize = 10 * 1024^3)  # 设置为 3 GiB
sample <- c("a","d","f","h")
##### Figure 1空转 #####
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
  FindClusters(resolution = 0.1,verbose = F) %>%
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
DimPlot(sp,label = T)
marker <- FindAllMarkers(sp,only.pos = T,logfc.threshold = 1)

sp@meta.data <- sp@meta.data %>%
  mutate(celltype = case_when(
    seurat_clusters %in% c("0","3") ~ "Tumor core",
    seurat_clusters %in% c("1","4") ~ "Invasive front",
    seurat_clusters %in% c("2") ~ "Adjacent stroma",
    TRUE ~ as.character(seurat_clusters)  # 使用原始列值作为默认情况
  ))

### 
sp <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sp.rds")
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


saveRDS(sp,file = "/data1/huanchangxiang/datasets/李莹雪论文/rds/sp.rds")

####
library(tidyverse)
library(magrittr)
a <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/a_celltype.rds")
d <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/d_celltype.rds")
f <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/f_celltype.rds")
h <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/h_celltype.rds")

rownames(sp@meta.data)[sp@meta.data$sample == "a"]
a@meta.data[rownames(sp@meta.data)[sp@meta.data$sample == "a"],"area"] = sp@meta.data[rownames(sp@meta.data)[sp@meta.data$sample == "a"],"celltype"]

rownames(sp@meta.data)[sp@meta.data$sample == "d"]
d@meta.data[rownames(sp@meta.data)[sp@meta.data$sample == "d"],"area"] = sp@meta.data[rownames(sp@meta.data)[sp@meta.data$sample == "d"],"celltype"]

rownames(sp@meta.data)[sp@meta.data$sample == "f"]
f@meta.data[rownames(sp@meta.data)[sp@meta.data$sample == "f"],"area"] = sp@meta.data[rownames(sp@meta.data)[sp@meta.data$sample == "f"],"celltype"]

rownames(sp@meta.data)[sp@meta.data$sample == "h"]
h@meta.data[rownames(sp@meta.data)[sp@meta.data$sample == "h"],"area"] = sp@meta.data[rownames(sp@meta.data)[sp@meta.data$sample == "h"],"celltype"]

meta_a <- h@meta.data[, 4:ncol(h@meta.data)] %>% 
  group_by(area) %>%  
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  rowwise() %>%  # 按行逐行计算
  mutate(
    total = sum(c_across(where(is.numeric)), na.rm = TRUE),  # 计算每行总和
    across(where(is.numeric), ~ .x / total * 100, .names = "{.col}_percent")  # 生成百分比列
  ) %>%
  ungroup() %>%
  select(area,ends_with("_percent")) %>% 
  select(-total_percent) %>%
  pivot_longer(cols = !area,values_to = "percent",names_to = "celltype")
meta_a$celltype <- gsub("_percent","",meta_a$celltype)
p1 <- ggplot(meta_a,aes(x=area,y=percent,fill = celltype)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("B_cells" = "#FFDAB9","Cancer_cells"="#FFA07A","Dendritic_cells"="#FF4500", 
                               "Endothelial_cells"="#FFD700","Fibroblasts"="#ADFF2F","Hepatocyte"="#FF6347",
                               "LSEC"="#32CD32","Mast_cells"="#6495ED","Myeloid_cells"="#1E90FF",
                               "Plasma_cells"="#4B0082","T_NK"="#FFC0CB")) + 
  theme(axis.line = element_line(size=1),
        legend.position = "none",
        panel.border = element_rect(size=0,colour = "black",fill = NA),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20)) +
  labs(x="",y="")

ggsave(plot = p1,filename = "Figure1空转_柱状图_h.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure1/",
       bg = "transparent",width = 4.1,height = 7.6)


#### Figure 1F ####
library(ggrepel) 
library(ggh4x)
library(showtext)
showtext_auto()
for (i in colnames(f@meta.data)[4:14]) {
  cells = i
  h_line <- GetTissueCoordinates(h)
  h_line$Fibroblasts <- h@meta.data[[cells]]
  h_line$Fibroblasts[is.na(h_line$Fibroblasts)] <- 0
  if (i %in% c("Cancer_cells","Fibroblasts","Hepatocyte")) {
    k = 1
  } else {
    k = 0.8
  }
  p1 <- ggplot(h_line, aes(x = imagecol, y = imagerow)) +
    geom_point(aes(color = Fibroblasts, fill = Fibroblasts), shape = 21, size = 2.5) +
    theme(axis.line = element_line(size=0),
          plot.background = element_rect(fill = "black", color = "black"),  # 整体背景
          panel.background = element_rect(fill = "black", color = "black"), # 绘图区域背景
          legend.position = "none",
          panel.border = element_rect(size=3,colour = "white",fill = NA),
          axis.text = element_text(size=20,colour = "white"),
          axis.title = element_text(size=0),
          axis.ticks.length = unit(10,"pt"),
          axis.ticks = element_line(size = 1.5, color = "white"),
          ggh4x.axis.ticks.minor = element_line(size = 0.8, color = "white"),
          panel.grid.major = element_line(color = "white", size = 0.8,linetype = "dashed"),  
          panel.grid.minor = element_line(color = "white", size = 0.8,linetype = "dashed")
    ) +
    scale_fill_gradientn(
      colours = colorRampPalette(c(
        "#313695FF", "#4575B4FF", "#74ADD1FF", "#ABD9E9FF", "#E0F3F8FF",
        "#FEE090FF", "#FDAE61FF", "#F46D43FF", "#D73027FF", "#A50026FF"
      ))(50),
      limits = c(min(h_line$Fibroblasts), max(h_line$Fibroblasts)*k)) +
    scale_color_gradientn(
      colours = colorRampPalette(c(
        "#313695FF", "#4575B4FF", "#74ADD1FF", "#ABD9E9FF", "#E0F3F8FF",
        "#FEE090FF", "#FDAE61FF", "#F46D43FF", "#D73027FF", "#A50026FF"
      ))(50),
      limits = c(min(h_line$Fibroblasts), max(h_line$Fibroblasts)*k)) +
    scale_x_continuous(guide = "axis_minor", 
                       breaks = seq(0, max(h_line$imagecol), by = 100)) +
    scale_y_continuous(guide = "axis_minor",
                       breaks = seq(0, max(h_line$imagerow), by = 100)) +
    labs(x="",y="")
  ggsave(plot = p1,filename = paste0("Figure1空转_f_",cells,".pdf"),
         path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure1/",
         width = 8,height = 7.52,bg="transparent")
  
}

### 提取坐标
library(reshape2)
library(ggplot2)
library(ggalt)
library(tidyverse)
library(Seurat)
sp <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/f_celltype.rds")
loc <- GetTissueCoordinates(sp)
n_groups <- length(seq(min(loc$imagerow), max(loc$imagerow), by = 50)) -1
loc$group <- cut(loc$imagerow,
                 breaks = seq(min(loc$imagerow), max(loc$imagerow), by = 50),right = F,
                 labels = paste0("Group", 1:n_groups))
loc$group[is.na(loc$group)] <- "Group9"

rownames(loc) <- paste0(rownames(loc),"_f")
loc$barcode <- rownames(loc)
sp@meta.data$barcode <- rownames(sp@meta.data)
loc_merge <- merge(loc,sp@meta.data,by="barcode") %>% na.omit()

meta <- loc_merge[,c(4,8:18)] %>%
  # 转换为长格式
  pivot_longer(cols = -group, names_to = "celltype", values_to = "freq") %>%
  group_by(group, celltype) %>%
  # 计算均值、标准差、样本量、标准误差
  summarise(
    mean = mean(freq, na.rm = TRUE),
    sd = sd(freq, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),  # 标准误差公式[1,3,7](@ref)
    .groups = "drop"
  ) %>%
  # 按细胞类型对均值归一化
  group_by(celltype) %>%
  mutate(
    norm_mean = (mean - min(mean)) / (max(mean) - min(mean)),
    # 归一化标准误差（保持与均值的相对比例）
    norm_sd = sd / (max(mean) - min(mean)),
    norm_se = se / (max(mean) - min(mean))
  ) %>%
  ungroup() 

meta$group <- factor(meta$group,levels = paste0("Group",1:9))
meta$Wavelength <- gsub("Group","",meta$group)
meta$Wavelength <- as.numeric(meta$Wavelength)
meta$Wavelength <- meta$Wavelength * 50

cols <- c("B_cells" = "#FFDAB9","Cancer_cells" = "#FFA07A","Dendritic_cells" = "#FF4500",
          "Endothelial_cells" = "#FFD700","Fibroblasts"="#ADFF2F","Hepatocyte"="#FF6347",
          "LSEC"="#32CD32","Mast_cells"="#6495ED","Myeloid_cells"="#1E90FF",
          "Plasma_cells"="#4B0082","T_NK"="#FFC0CB")
p2 <- ggplot(data = meta, aes(x = Wavelength, y = norm_mean, group = celltype, colour = celltype)) +
  geom_errorbar(aes(ymin = norm_mean-norm_se,ymax=norm_mean+norm_se),
               width = 5,alpha=0.8,size=1.5) +
  geom_point(size = 3) +
  geom_xspline(spline_shape = -0.5, size = 20, lineend = "round") +
  scale_colour_manual(values = cols) +
  scale_y_continuous(
    breaks = scales::extended_breaks(n = 5),
    minor_breaks = NULL
  ) +
  scale_x_continuous(
    breaks = scales::extended_breaks(n = 5),
    minor_breaks = NULL
  ) +
  theme(
    # 背景设置为黑色
    panel.background = element_rect(fill = "black"),  # 绘图区域背景
    plot.background = element_rect(fill = "black"),   # 整个图像背景
    # 网格线颜色调整（避免与黑色背景冲突）
    panel.grid.major = element_line(color = "white", size = 0.3),
    panel.grid.minor = element_line(color = "white", size = 0.2),
    # 边框和坐标轴颜色
    panel.border = element_rect(color = "white", fill = NA, size = 0.5),
    axis.line = element_line(color = "white"),        # 坐标轴线
    axis.text = element_text(color = "white",size = 15),        # 坐标轴文本
    # 其他设置
    plot.margin = margin(0, 0, 0, 0, unit = "cm"),
    legend.position = "none",
    axis.title.x = element_blank()
  )
ggsave(plot = p2,filename = "Figure1G_细胞随距离变化图.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure1/",
       width = 8,height = 2)
