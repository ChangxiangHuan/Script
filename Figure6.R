#### Figure6C
library(Seurat)
library(ComplexHeatmap)
library(tidyverse)
library(ggplot2)

sc <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/merge/merge_raw.rds")
Idents(sc) <- sc$celltype_cellchat
gene <- c("FAP","MMP1","TNC","POSTN","LOXL1","EGF","VEGFA","PDGFB","TGFB1","CXCL6","SNAI2")
avg <- AverageExpression(sc,assays = "SCT",slot = "data") %>% as.data.frame()
avg_gene <- avg[gene,] %>% na.omit()
colnames(avg_gene) <- gsub("SCT.","",colnames(avg_gene))

min_max_nor <- function(x){(x-min(x)) / (max(x)-min(x))}

nor_avg_gene <- t(apply(avg_gene, 1, min_max_nor))

#### 生成注释
afsplit = data.frame(gene = gene,
                     type = c(rep("ECM remodel",5),rep("Growth factor",4),rep("Immune response",1),
                              rep("Anti-apoptosis",1))
                     )



#### 绘制热图
pdf("/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure6/Figure6C_heatmap.pdf",
    width = 10,height = 20)
ht <- ComplexHeatmap::Heatmap(nor_avg_gene,
                        cluster_columns = T,cluster_rows = F,
                        show_column_dend = F,show_row_dend = F,
                        col=circlize::colorRamp2(c(0,0.5,1),c("blue","#f7f7f7","red")),
                        border_gp = gpar(col = "black",lwd=2),  ### 热图边框与粗细
                        rect_gp   = gpar(col = "white",lwd=2),
                        row_split = afsplit$type,
                        left_annotation = rowAnnotation(fun = afsplit$type,
                                                        col=list(fun=c("ECM remodel" = "red",
                                                                        "Growth factor" = "pink",
                                                                        "Immune response" = "orange",
                                                                        "Anti-apoptosis" = "skyblue"
                                                                        )),
                                                        show_annotation_name=F
                                                        )
                        )
draw(ht)
dev.off()


library(maftools)
library(tidyverse)
library(data.table)
load("/data1/huanchangxiang/datasets/TCGA/TCGA-LIHC_SNP.Rdata")
maf <- data
dim(maf)
maf <- read.maf(maf)
#计算TMB
coad.tmb<- tmb(maf,captureSize=38,logScale=T)
coad.tmb$Tumor_Sample_Barcode <- substr(coad.tmb$Tumor_Sample_Barcode,1,16)
coad.tmb$Tumor_Sample_Barcode <- make.names(coad.tmb$Tumor_Sample_Barcode)
gene <- read.table("/data1/huanchangxiang/datasets/TCGA/datasets/LICH_gene.txt",
                   row.names=1,header = T)
sample <- coad.tmb$Tumor_Sample_Barcode
# 检查sample中的缺失列名
missing_cols <- sample[!sample %in% colnames(gene)]
coad.tmb <- coad.tmb[!coad.tmb$Tumor_Sample_Barcode %in% missing_cols]
gene_tumor <- gene[,coad.tmb$Tumor_Sample_Barcode]
gene_tumor_POSTN <- gene_tumor["POSTN",] %>% t() %>% as.data.frame()
# 重构barcode
colnames(coad.tmb)[1] <- "barcode"
gene_tumor_POSTN$barcode <- rownames(gene_tumor_POSTN)
# merge
data <- merge(gene_tumor_POSTN,coad.tmb,by = "barcode")
library(ggstatsplot)
ggstatsplot::ggscatterstats(data=data,
                            x=POSTN,y=total_perMB_log,
                            type = "noparametric",
                            marginal.type = "density",
                            xlab = "POSTN",ylab = "TMB",
                            bf.messages=FALSE,
                            xfill = "#E69F00",            # X轴颜色
                            yfill = "#56B4E9"             # Y轴颜色
)

### Figure 6D
library(tidyverse)
library(data.table)
library(ggplot2)
library(corpcor)
library(corrplot)

##### 导入数据
gene <- read.table("/data1/huanchangxiang/datasets/TCGA/datasets/LICH_gene.txt",
                   row.names = 1,header = T)
phe  <- fread("/data1/huanchangxiang/datasets/TCGA/datasets/TCGA-LIHC.clinical.tsv.gz")

#### 筛选肿瘤样本
colnames(gene) <- gsub("\\.","-",colnames(gene))
phe_tumor  <- phe[phe$tissue_type.samples == "Tumor"]
sample     <- phe$sample
sample     <- sample[sample %in% colnames(gene)]

#### 选择相关性基因
select_gene <- c("POSTN","KDR","IL10RB","IL10","HAVCR2","CTLA4","CSF1R","CD96","CD274","CD244","CD160","BTLA",
                 "ADORA2A","IDO1","VTCN1","TIGIT","TGFBR1","TGFB1","PDCD1LG2","PDCD1","LGALS9","LAG3",
                 "KIR2DL3","KIR2DL1")
gene_tumor <- gene[select_gene,sample]

##### 用于与POSTN分析的相关性基因
cor_gene <- c("KDR","IL10RB","IL10","HAVCR2","CTLA4","CSF1R","CD96","CD274","CD244","CD160","BTLA",
              "ADORA2A","IDO1","VTCN1","TIGIT","TGFBR1","TGFB1","PDCD1LG2","PDCD1","LGALS9","LAG3",
              "KIR2DL3","KIR2DL1")

cor_list <- lapply(cor_gene, function(i){
  gene_tumor_cor <- gene_tumor[c("POSTN",i),] %>% t() %>% as.data.frame()
  cor_result <- cor.test(gene_tumor_cor[,"POSTN"],gene_tumor_cor[,i])
  r <- cor_result$estimate
})
names(cor_list) <- cor_gene

result <- do.call(cbind,cor_list) %>% as.data.frame()
result$group <- "POSTN"
##### 绘制雷达图
library(ggradar)
library(scales)

ggradar(result,
       font.radar="sans",
       values.radar=c(-0.6,0,0.6),
       grid.min=-0.6,
       grid.mid=0,
       grid.max=0.6,
       group.line.width=1,
       group.point.size=5,
       background.circle.colour='white',#雷达图背景颜色
       gridline.min.linetype="longdash",#中间线条类型，#twodash,longdash,dotted,dashed,solid,blank
       gridline.mid.linetype="longdash",
       gridline.max.linetype="longdash",
       gridline.mid.colour='black',#最里圈网格线颜色
       gridline.min.colour='black',#中间网格线颜色
       gridline.max.colour='black',#最外圈网格线颜色
       axis.line.colour="black",#坐标轴颜色
       axis.label.size=5,#轴标签大小
       legend.position='right',#图例位置
       legend.title='Group',#图例名称
       legend.text.size=10#图例标签大小
      
)+
  theme(plot.title=element_text(size=20,hjust=0.5),
        legend.title=element_text(size=15)
  )

#### Figure6I
library(ggsignif)
library(rstatix)
library(tidyverse)
library(data.table)
##### 导入数据
gene <- read.table("/data1/huanchangxiang/datasets/TCGA/datasets/LICH_gene.txt",
                   row.names = 1,header = T)
phe  <- fread("/data1/huanchangxiang/datasets/TCGA/datasets/TCGA-LIHC.clinical.tsv.gz")
#### 筛选肿瘤样本
colnames(gene) <- gsub("\\.","-",colnames(gene))
phe_tumor  <- phe[phe$sample_type.samples == "Primary Tumor"]
sample     <- phe$sample
sample     <- sample[sample %in% colnames(gene)]

#### 肿瘤基因表达矩阵
gene_tumor <- gene[,sample]

nrow(gene_tumor)
gene_tumor = gene_tumor[rowSums(gene_tumor)>30,]
nrow(gene_tumor)

### 筛选sample
gene_tumor = gene_tumor[apply(gene_tumor, 1, function(x) sum(x > 0) > 0.3*ncol(gene_tumor)), ]

##### 增加Symbol列
gene_tumor <- gene_tumor %>%
      mutate(Symbol = rownames(gene_tumor)) %>%
      select(Symbol,everything())
##### 为满足100M要求，均分为3份
# 假设gene_tumor为待分割的数据框，首列为Symbol（需保留在所有子集中）
cols <- colnames(gene_tumor)[-1]  # 排除第1列Symbol
split_size <- ceiling(length(cols) / 3)  # 计算每份列数

# 分割为3份（保留Symbol列）
list_of_dfs <- lapply(0:2, function(i) {
  start <- i * split_size + 1
  end <- min((i + 1) * split_size, length(cols))
  gene_tumor %>% select(Symbol, cols[start:end])
})


write.table(list_of_dfs[[1]],file = "/data1/huanchangxiang/datasets/TCGA/datasets/gene_tumor_1.txt",
            sep="\t",row.names = F,col.names = T,quote = F)
write.table(list_of_dfs[[2]],file = "/data1/huanchangxiang/datasets/TCGA/datasets/gene_tumor_2.txt",
            sep="\t",row.names = F,col.names = T,quote = F)
write.table(list_of_dfs[[3]],file = "/data1/huanchangxiang/datasets/TCGA/datasets/gene_tumor_3.txt",
            sep="\t",row.names = F,col.names = T,quote = F)


####
immunCellAI1 <- read.table("/data1/huanchangxiang/datasets/TCGA/datasets/ImmuuCellAI_1.txt")
immunCellAI2 <- read.table("/data1/huanchangxiang/datasets/TCGA/datasets/ImmuuCellAI_2.txt")
immunCellAI3 <- read.table("/data1/huanchangxiang/datasets/TCGA/datasets/ImmuuCellAI_3.txt")

immunCellAI <- rbind(immunCellAI1,immunCellAI2,immunCellAI3) %>% 
                mutate(sample = rownames(.))

gene_tumor_POSTN <- gene_tumor["POSTN",] %>% t() %>% as.data.frame() %>%
                    mutate(sample = rownames(.))

data <- merge(immunCellAI,gene_tumor_POSTN,by = "sample")
data$Response <- as.factor(data$Response)

# 执行t检验
stat_test <- data %>% 
  t_test(POSTN ~ response)  # 根据实际分组列名调整

data <- data %>%  mutate(
  response = case_when(
    Response %in% "0" ~ "Resistant",
    Response %in% "1" ~ "Sensitive"
  )
)
ggplot(data,aes(x=response,y = POSTN,fill = response)) +
  geom_boxplot() +
  geom_signif(
    comparisons = list(c("Resistant", "Sensitive")),  # 替换为实际分组名
    map_signif_level = TRUE,   # 显示星号而非p值（***:p<0.001, **:p<0.01等）
    test = "t.test",          # 指定检验方法
    y_position = max(data$POSTN) * 1.1,  # 标记位置高于数据最大值
    tip_length = 0.01,        # 竖线长度
    size = 0.5,                  # 线宽
    textsize = 5               # 星号大小
  ) +
  scale_fill_manual(values = c("Resistant" = "#E41A1C", "Sensitive" = "#377EB8")) +
  theme(axis.line = element_line(colour = "black",size = 1.2),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.position = "none",
        axis.title.y = element_text(size = 20)
        ) +
  labs(x=NULL,y="POSTN(log2(FPKM+1))")
