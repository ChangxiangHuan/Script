library(Seurat)
library(ggplot2)
library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(GSVA)
library(GSEABase)
setwd("/data1/huanchangxiang/datasets/TCGA/")
##### 
exp_cells <- read.csv("./output/theta_state.csv",row.names = 1) 
exp_cells <- read.table("/data1/huanchangxiang/datasets/TCGA/datasets/LICH_gene.txt",
                        header = T,row.names = 1)  %>%t() %>% as.data.frame()
rownames(exp_cells) <- gsub("\\.","-",rownames(exp_cells))
#### 导入临床数据
phe       <- fread("./datasets/TCGA-LIHC.clinical.tsv.gz")
os        <- fread("./datasets/TCGA-LIHC.survival.tsv.gz")
####筛选Tumor组织
phe_tumor <- phe[phe$tissue_type.samples == "Tumor"]
tumor_id <- phe_tumor$sample %>% unique()
if (any(tumor_id %in% rownames(exp_cells))) {
  exp_cells_tumor <- exp_cells[tumor_id,] %>% na.omit() %>% t()
  os_tumor        <- os[os$sample %in% tumor_id]
}

# 初始化结果存储
results_df <- data.frame(
  gene_set_size = integer(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)
for (i in 3:3) {
  gene_sets <- GeneSetCollection(list(GeneSet(marker$gene[1:i],
                                              setName = "my_gene_set")))
  gsvapar <- GSVA::gsvaParam(exprData = exp_cells_tumor,geneSets = gene_sets,maxDiff = T)
  gsva_es <- gsva(gsvapar)
  gsva <- t(gsva_es) %>% as.data.frame()
  colnames(gsva) <- "Mac_02_APOE"
  Mac <- "Mac_02_APOE"
  gsva$sample <- rownames(gsva)
  a <- merge(os_tumor,gsva,by="sample")
  value <- surv_cutpoint(a,time = "OS.time",
                         event = "OS",variables = "Mac_02_APOE")
  cut_off<-as.numeric(value[["cutpoint"]][1,1])
  surv_tmb.cat <- surv_categorize(value)
  cox_model <- coxph(Surv(OS.time, OS) ~ Mac_02_APOE, data = surv_tmb.cat)
  summary_cox <- summary(cox_model)
  # 提取P值
  p_value <- signif(summary_cox$coefficients[1, "Pr(>|z|)"], 3)  # 保留3位有效数字
  
  results_df <- rbind(results_df, data.frame(
    gene_set_size = i,
    p_value = p_value
  ))
  
  cat("Gene set size:", i, "| P =", p_value, "\n")
}

fit1 <- survfit(Surv(OS.time,OS)~Mac_02_APOE,data=surv_tmb.cat)

#### 绘图 
p1 <- ggsurvplot(fit1,
                 pval = TRUE, 
                 conf.int = TRUE,
                 risk.table = FALSE, # Add risk table
                 risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "hv", # Specify median survival
                 ggtheme = theme_bw(),
                 p.adjust.method = "bonferroni",
                 palette = c("#EE8386", "#93C695","orange","pink"),
                 xlab=NULL,ylab=NULL)

p1










# Mac_01_SPP1分组
exp_cells_tumor_group <- exp_cells_tumor %>%
  mutate(mac_group = ifelse(get(Mac) > Mac_cutoff,paste0(Mac,"_high"),paste0(Mac,"_low")),
         os_group  = paste0(mac_group)) %>% 
  tibble::rownames_to_column(var = "sample") %>%
  select(mac_group,os_group,sample)
surv_data <- merge(os_tumor,exp_cells_tumor_group,
                   by.x="sample",by.y="sample")
surv_data$OS.time <- surv_data$OS.time/30
fit1 <- survfit(Surv(OS.time,OS)~os_group,data=surv_data)
fit1

p2 <- ggsurvplot(fit1,
                 pval = TRUE, 
                 conf.int = TRUE,
                 risk.table = FALSE, # Add risk table
                 risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "hv", # Specify median survival
                 ggtheme = theme_bw(),
                 p.adjust.method = "bonferroni",
                 palette = c("#990000", "#000099","orange","pink"),
                 xlab=NULL,ylab=NULL)

p2


# HSC分组
exp_cells_tumor_group <- exp_cells_tumor %>%
  mutate(HSC_group = ifelse(get(HSC) > HSC_cutoff,paste0(HSC,"_high"),paste0(HSC,"_low")),
         os_group  = paste0(HSC_group)) %>% 
  tibble::rownames_to_column(var = "sample") %>%
  select(HSC_group,os_group,sample)
surv_data <- merge(os_tumor,exp_cells_tumor_group,
                   by.x="sample",by.y="sample")
surv_data$OS.time <- surv_data$OS.time/30
fit1 <- survfit(Surv(OS.time,OS)~os_group,data=surv_data)
fit1

p3 <- ggsurvplot(fit1,
                 pval = TRUE, 
                 conf.int = TRUE,
                 risk.table = FALSE, # Add risk table
                 risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "hv", # Specify median survival
                 ggtheme = theme_bw(),
                 p.adjust.method = "bonferroni",
                 palette = c("#990000", "#000099","orange","pink"),
                 xlab=NULL,ylab=NULL)

p3


# LSEC分组
exp_cells_tumor_group <- exp_cells_tumor %>%
  mutate(LSEC_group = ifelse(get(LSEC) > LSEC_cutoff,paste0(LSEC,"_high"),paste0(LSEC,"_low")),
         os_group  = paste0(LSEC_group)) %>% 
  tibble::rownames_to_column(var = "sample") %>%
  select(LSEC_group,os_group,sample)
surv_data <- merge(os_tumor,exp_cells_tumor_group,
                   by.x="sample",by.y="sample")
surv_data$OS.time <- surv_data$OS.time/30
fit1 <- survfit(Surv(OS.time,OS)~os_group,data=surv_data)
fit1

p4 <- ggsurvplot(fit1,
                 pval = TRUE, 
                 conf.int = TRUE,
                 risk.table = FALSE, # Add risk table
                 risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "hv", # Specify median survival
                 ggtheme = theme_bw(),
                 p.adjust.method = "bonferroni",
                 palette = c("#990000", "#000099","orange","pink"),
                 xlab=NULL,ylab=NULL)

p4



# 适合多组比较
fit2 <- survdiff(Surv(OS.time,OS)~os_group,data=surv_data)
fit2


