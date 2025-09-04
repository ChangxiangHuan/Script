library(Seurat)
library(SeuratObject)
library(ggplot2)
library(gridExtra)
library(pals)
library(spatialLIBD)
library(ExperimentHub)
library(harmony)
library(ggplot2)
library(patchwork)
library(infercnv)
library(future)
plan("sequential")
options(scipen = 100)
sc_hepa <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sc_hepa_filter.rds")
sc_plasma <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/Endo_second.rds")

#####
meta_df <- as.data.frame(sc_hepa@meta.data)
sampled_cells <- lapply(
  split(meta_df, meta_df$celltype),  # 按 celltype 拆分数据框
  function(cluster_df) {
    n_cells <- nrow(cluster_df)
    if (n_cells > 5000) {
      # 计算80%的样本量（四舍五入取整）
      sample_size <- round(n_cells * 0.5)  
      # 随机抽取指定数量的细胞名[1,2](@ref)
      sample(rownames(cluster_df), size = sample_size)  
    } else {
      rownames(cluster_df)  # 直接返回全部细胞名
    }
  }
)
final_cells <- unlist(sampled_cells, use.names = FALSE)
sc_hepa_filter <- subset(sc_hepa,cells=final_cells)
rm(sc_hepa)
###########################
epiMat <- as.matrix(GetAssayData(sc_hepa_filter, slot = 'counts',assay='RNA'))
refMat <- as.matrix(GetAssayData(sc_plasma, slot = 'counts',assay='RNA'))
###########################
gene <- intersect(rownames(refMat), rownames(epiMat))
epiMat <- epiMat[gene,]
refMat <- refMat[gene,]
identical(nrow(refMat),nrow(epiMat))
###########################
dat=cbind(epiMat,refMat)

groupinfo = data.frame(v1=colnames(dat),
                       v2=c(sc_hepa_filter$celltype,rep("ref",ncol(refMat))))
colnames(groupinfo) <- c("barcode","group")
write.table(groupinfo,file = "/data1/huanchangxiang/datasets/public_single/rds/infercnv/groupFile.txt",
            sep = '\t',quote = F,col.names = F,row.names = F)
setwd("/data1/huanchangxiang/datasets/public_single/rds/infercnv")
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = dat,
                                     annotations_file = 'groupFile.txt',
                                     delim = '\t',
                                     gene_order_file = 'hg38_gencode_v27.txt',
                                     ref_group_names = c("ref"))
outdir <- "/data1/huanchangxiang/datasets/public_single/rds/infercnv/"
infercnv_over <- infercnv::run(infercnv_obj,
                               cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir = outdir, # output files
                               cluster_by_groups = T, # cluster or not
                               analysis_mode = "subclusters",
                               denoise = F,
                               HMM = F,
                               #tumor_subcluster_partition_method = "random_trees",
                               #HMM_type = "i6",
                               #BayesMaxPNormal = 0,
                               num_threads = 10,
                               write_expr_matrix = T,
                               output_format = "pdf")

