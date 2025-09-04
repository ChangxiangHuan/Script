library(monocle)
library(ggplot2)
library(Seurat)
library(patchwork)
library(tidyverse)
Fib <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/Fib_second_clusters.rds")
DefaultAssay(Fib) <- "RNA"
expr_matrix <- GetAssayData(Fib,assay = "RNA",slot = "counts")
# 获取分组信息
p_data <- Fib@meta.data
f_data <- data.frame(gene_short_name=row.names(Fib),row.names = row.names(Fib))

### 构建CDS对象
pd <- new("AnnotatedDataFrame",data = p_data)
fd <- new("AnnotatedDataFrame",data = f_data)

cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily = negbinomial.size()
)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


Fib <- FindVariableFeatures(Fib,nfeatures = 3000)
expressed_genes <- VariableFeatures(Fib)

diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~seurat_clusters",cores = 10)
head(diff)

deg <- subset(diff,qval < 0.01)
deg <- deg[order(deg$qval,decreasing = F),]

##### 轨迹基因可视化 ####
ordergene <- rownames(deg)
cds <- setOrderingFilter(cds,ordergene)

plot_ordering_genes(cds)

cds <- reduceDimension(cds,max_components = 2, reduction_method = "DDRTree")

cds <- orderCells(cds)
cds <- orderCells(cds,root_state = 4)

p1 <- plot_cell_trajectory(cds,
                           color_by = "Pseudotime",
                           size=1,
                           show_backbone = T)
p1
ggsave(plot = p1,filename = "monocle_Fib_Pseudotime.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure4/",
       width = 10,height = 10)
p2 <- plot_cell_trajectory(cds,
                           color_by = "seurat_clusters",
                           size=1,
                           show_backbone = T)
p2
ggsave(plot = p2,filename = "monocle_Fib_seurat_clusters.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure4/",
       width = 10,height = 10)

######### Figure4B #########
library(ggpubr)
df <- pData(cds)
p3 <- ggplot(df,aes(Pseudotime,colour = seurat_clusters,fill = seurat_clusters)) + 
  geom_density(bw=0.5,size=1,alpha=0.5) + 
  theme_classic2()
ggsave(plot = p3,filename = "monocle_Fib_pseudotime_denstiy.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure4/",
       width = 10,height = 10)
######### Figure4C #########
s.genes <- c("COL1A2","COL1A1","POSTN","GJ")
plot_genes_in_pseudotime(cds[s.genes,],color_by = "seurat_clusters")

######### Figure4C #########
library(ClusterGVis)
diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)")
num_cells_expressed <- rowSums(exprs(cds) > 0)
diff_test_res$num_cells_expressed <- num_cells_expressed[row.names(diff_test_res)]

diff_test_res <- diff_test_res %>% filter(num_cells_expressed > 10,
                                          use_for_ordering == "TRUE")  # 过滤掉表达细胞数小于 10 的基因
genelist <- row.names(diff_test_res)
p4 <- plot_pseudotime_heatmap2(cds[genelist,],
                              cores = 1,
                              num_clusters = 2,
                              show_rownames = F,
                              return_heatmap = F)
genes <- c("COL1A1","COL1A2","ITGA1","ITGA11","POSTN",  ### module1
           "HES4","SRM","GJA4","YBX1","ARTN","CCN1"
           )
visCluster(object=p4,plot.type="heatmap",
           show_row_dend = F,markGenes = genes)                             


meta <- data.frame(p4$long.res) %>% select(cluster,gene) %>% distinct()
write_csv(meta,file = "/data1/huanchangxiang/datasets/李莹雪论文/rds/monocle2_module_gene.csv")

######### Figure4D #########
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(stringr)
library(Seurat)

genes1 <- read.csv("/data1/huanchangxiang/datasets/李莹雪论文/rds/monocle2_module_gene.csv") %>%
  filter(cluster==1)

Gene_ID <- bitr(genes1$gene,fromType = "SYMBOL",toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db")
enrich.go <- enrichGO(gene = Gene_ID$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "fdr",
                      minGSSize = 10,
                      maxGSSize = 500,
                      qvalueCutoff = 0.05)
dotplot(enrich.go)

write.csv(enrich.go@result,file = "/data1/huanchangxiang/datasets/李莹雪论文/rds/GO_module1.csv")

### module2
genes2 <- read.csv("/data1/huanchangxiang/datasets/李莹雪论文/rds/monocle2_module_gene.csv") %>%
  filter(cluster==2)

Gene_ID <- bitr(genes2$gene,fromType = "SYMBOL",toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db")
enrich.go <- enrichGO(gene = Gene_ID$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "fdr",
                      minGSSize = 10,
                      maxGSSize = 500,
                      qvalueCutoff = 0.05)
dotplot(enrich.go)

write.csv(enrich.go@result,file = "/data1/huanchangxiang/datasets/李莹雪论文/rds/GO_module1.csv")


######### Figure4E #########
library(GENIE3)
Fib <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/Fib_second_clusters.rds")
myCAFs <- subset(Fib,idents="1")

counts <- GetAssayData(myCAFs,assay = "RNA",slot = "counts") %>% as.matrix()

regulators <- read.table("/data1/huanchangxiang/datasets/public_single/rds/pyscenic/hs_hgnc_tfs.txt")
regulator <- regulators$V1
rm(regulators)

# 检查哪些调控基因不在表达矩阵中
missing_genes <- regulator[!regulator %in% rownames(counts)]

# 如果有缺失的基因，从调控基因列表中移除
if(length(missing_genes) > 0) {
  cat("以下调控基因未在表达矩阵中找到，将被移除：", paste(missing_genes, collapse = ", "), "\n")
  regulator <- regulator[regulator %in% rownames(counts)]
} else {
  cat("所有调控基因都已存在于表达矩阵中。\n")
}

# 现在 regulator 中只包含存在于表达矩阵中的基因


weightMatrix <- GENIE3(counts,
                       regulators = regulator,  # 使用更新后的调控基因列表
                       targets = NULL,
                       treeMethod = "RF",
                       K = "sqrt",
                       nTrees = 50,
                       nCores = 10,
                       returnMatrix = TRUE,
                       verbose = TRUE)

linkList_top30 <- getLinkList(weightMatrix,
                               threshold = 0.05)
linkList_top30 <- linkList_top30[,1:2]
write.table(linkList_top30,"/data1/huanchangxiang/datasets/李莹雪论文/rds/TF_30.txt",
            quote=F,row.names = F,sep = "\t")
