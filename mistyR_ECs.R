library(mistyR)
library(CARD)
library(Seurat)
library(mistyR)
library(distances)
library(future)
library(tidyverse) 
library(recipes)
library(parallel)
plan(multisession, workers = min(detectCores(), 40))
samples <- c("a","d","f","h")
for (data in samples) {
  sp <- readRDS(paste0("/data1/huanchangxiang/datasets/李莹雪论文/rds/",data,".rds"))
  rownames(sp@meta.data) <- gsub("_a","",rownames(sp@meta.data))
  sp@meta.data[,"ECs3"] <- NULL
  colnames(sp@meta.data)[colnames(sp@meta.data) %in% "ECs4"] <- "ECs3"
  colnames(sp@meta.data)[colnames(sp@meta.data) %in% "ECs5"] <- "ECs4"
  colnames(sp@meta.data)[colnames(sp@meta.data) %in% "ECs6"] <- "ECs5"
  colnames(sp@meta.data)[colnames(sp@meta.data) %in% "ECs7"] <- "ECs6"
  num = 4:21
  composition <- sp@meta.data[,num]
  calculate_cells_proportions <- function(df) { 
    # 对每一行(患者)进行计算 
    proportions <- t(apply(df, 1, function(row) { 
      # 计算该行的总和 
      row_sum <- sum(row) 
      # 计算每个指标占总和的比例 
      return(row / row_sum) })) 
    return(proportions) }
  composition <- calculate_cells_proportions(composition) %>% as.data.frame()
  colnames(composition) <- gsub(" ", "_", colnames(composition), fixed = TRUE)
  composition[is.na(composition)] <- 0
  
  geometry <- GetTissueCoordinates(sp,cols=c("row","col"),scale=NULL)
  
  geom_dist <- as.matrix(distances(geometry))
  dist_nn <- apply(geom_dist, 1, function(x) (sort(x)[2]))
  paraview_radius <- ceiling(mean(dist_nn+ sd(dist_nn)))
  
  ### create view
  GBM_views <- create_initial_view(data.frame(composition))
  
  
  
  run_misty(GBM_views,"/data1/huanchangxiang/datasets/李莹雪论文/rds/mistyR/")
  
  ### 下游result分析
  misty_result <- collect_results("/data1/huanchangxiang/datasets/李莹雪论文/rds/mistyR/")
  
  # misty_result %>%
  #    plot_interaction_heatmap(view="intra",clean = T,cutoff = 0)
  #  
  # result <- misty_result$importances.aggregated %>%
  #    filter(view=="intra",Predictor == "myCAFs") %>%
  #    arrange(-Importance)
  # 
  # 
  # misty_result %>% 
  #   plot_interaction_heatmap(view = "intra", clean = F, 
  #                            trim = 0.05, trim.measure = "multi.R2", cutoff = 0.5)
  # 
  # SpatialFeaturePlot(sp,features = c("myCAFs","LSEC","Endothelial_cells"))
  
  
  result <- misty_result$importances.aggregated %>%
    filter(view=="intra",Predictor == "myCAFs",Target %in% c("ECs1","ECs2","ECs3","ECs4","ECs5","ECs6")) %>%
    arrange(-Importance) %>% 
    select(Target,view,Importance) %>% na.omit() %>%
    pivot_wider(.,names_from = Target,values_from = Importance)
  colnames(result)[1] <- "group"
  
  
  ##### 雷达图可视化 ####
  library(ggradar)
  ggradar(result, 
          grid.min = -1, grid.mid = 1, grid.max = 2,  # 设置网格范围
          group.colours = c( "red"), 
          legend.position = "bottom")
  ggsave(filename = paste0("Figure6A_",data,"ECs.pdf"),
         path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure6/",
         width = 10,height = 10)
}