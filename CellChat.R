library(ggplot2)
library(Seurat)
library(CellChat)
library(harmony)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 1024 * 1024 * 1024 * 20)
sce <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/merge/merge_raw.rds")
Idents(sce) <- sce$NorT
sce_T <- subset(sce,idents="T")
sce_N <- subset(sce,idents="N")
for (data in c("sce_T","sce_N")) {
sce <- get(data)
### staring from Seurat ###
data.input <- sce[["SCT"]]@data ###normalized data
Idents(sce) <- sce$celltype_cellchat
labels <- Idents(sce)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels


### Create a CellChat obj ###
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

### set the L-R interaction database ###
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB) # use Secreted Signaling
# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 20) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

###  probability and infer cellular communication network ###
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")

#### 过滤低细胞数的细胞类型
cellchat <- filterCommunication(cellchat, min.cells = 20)

#### 按照同路统计
# df.net <- subsetCommunication(cellchat, signaling = c("PDGF", "TGFb","CCL"))
# 
# cellchat <- computeCommunProbPathway(cellchat)

#### 计算聚合的细胞间通信网络
### Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat,file = paste0("/data1/huanchangxiang/datasets/李莹雪论文/rds/",data,".rds"))
}

###### 图1
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength")

####### 图2
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#######
pathways.show <- c("PDGF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

#### 画图-myCAFs-other cells
netVisual_bubble(cellchat,sources.use = "myCAFs",targets.use = c("Cancer_cells","Fibroblasts","Mast_cells","LSEC","T_NK","Plasma_cells",
                                                                 "Endothelial_cells","Myeloid_cells","B_cells","Hepatocyte","Dendritic_cells"))

saveRDS(cellchat,"/data1/huanchangxiang/datasets/李莹雪论文/rds/cellchat.rds")


### KS plot
library(CellChat)
library(RColorBrewer)
library(ggplot2)
sce_T <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sce_T.rds")
sce_N <- readRDS("/data1/huanchangxiang/datasets/李莹雪论文/rds/sce_N.rds")

sce_T_net <- subsetCommunication(sce_T)
sce_N_net <- subsetCommunication(sce_N)

net1 <- subset(sce_T_net,source == "myCAFs")
net2 <- subset(sce_N_net,source == "myCAFs")

net1$cell_inter <- paste0(net1$source,"->",net1$target)
net2$cell_inter <- paste0(net2$source,"->",net2$target)

net1$group <- "Tumor"
net2$group <- "Normal"

df.net <- rbind(net1,net2) %>% filter(prob > 0.01)

#分娩图的形式做就ok了
p1 <- ggplot(df.net,aes(x=cell_inter,y=interaction_name)) +
  geom_point(aes(size=prob,color=prob)) +
  geom_point(shape=21,aes(size=prob))+
  facet_wrap(~group)+
  scale_color_gradientn('Communication\nProbability', 
                        colors=colorRampPalette(rev(brewer.pal(9, "PRGn")))(100)) +
  theme_bw() +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size=10),
        axis.text.y = element_text(size=8, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

ggsave(plot = p1,filename = "Figure6B_TvsN.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure6/",
       width = 10,height = 10)

### 只关注myCAFs
net3 <- subset(sce_T_net,target=="myCAFs")
net3$cell_inter <- paste0(net3$source,"->",net3$target)

net3$group <- "receiver"
net1$group <- "sender"

df.net2 <- rbind(net3,net1) %>% filter(!cell_inter %in% c("myCAFs->myCAFs"))

p2 <- ggplot(df.net2,aes(x=cell_inter,y=interaction_name)) +
  geom_point(aes(size=prob,color=prob)) +
  geom_point(shape=21,aes(size=prob))+
  facet_wrap(~group,scales = "free") +
  scale_color_gradientn('Communication\nProbability', 
                        colors=colorRampPalette(rev(brewer.pal(9, "PRGn")))(100)) +
  theme_bw() +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size=10),
        axis.text.y = element_text(size=8, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
ggsave(plot = p2,filename = "Figure6B_SendervsReceiver.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure6/",
       width = 10,height = 10)

library(RColorBrewer)
library(ggh4x)
# 创建颜色映射数据框（确保每个 interaction_name 对应唯一的 annotation）
color_data <- net1 %>%
  distinct(interaction_name, annotation) %>%
  arrange(interaction_name) 

annotation_levels <- unique(color_data$annotation)
color_palette <- scales::brewer_pal(palette = "Set1")(length(annotation_levels))
names(color_palette) <- annotation_levels  
label_colors <- color_palette[color_data$annotation]



p3 <- ggplot(net1, aes(x = cell_inter, y = interaction_name)) +
  geom_point(aes(size = prob, color = prob)) +
  geom_point(shape = 21, aes(size = prob)) +
  scale_color_gradientn(
    'Communication\nProbability',
    colors = colorRampPalette(rev(brewer.pal(9, "PRGn")))(100)
  ) +
  scale_y_discrete(guide = guide_axis_color(colour = label_colors)) +  # 关键修改
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")
  )
p3

ggsave(plot = p3,filename = "Figure6B_annotation_myCAFs.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure6/",
       width = 10,height = 10)

net_ITGA <- net1[grepl("ITGA",net1$receptor,ignore.case = T),] %>%
  select(cell_inter,receptor) %>%
  group_by(cell_inter,receptor) %>%
  summarise(value=n())

p4 <- ggplot(net_ITGA,aes(x=cell_inter,y=value,fill = receptor)) +
  geom_bar(stat = "identity",
           position = "stack") +
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5),
        axis.title.x = element_blank()) +
  labs(y = NULL,x=NULL) +
  coord_flip() +
  theme_classic()
ggsave(plot = p4,filename = "Figure6C_receptors_type.pdf",
       path = "/data1/huanchangxiang/datasets/李莹雪论文/Figure/Figure6/",
       width = 10,height = 10)
