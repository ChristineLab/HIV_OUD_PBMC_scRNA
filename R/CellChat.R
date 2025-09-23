suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(tidyverse)
  library(CellChat)
  library(ggalluvial)
  library(NMF)
})

#Load the Seurat object
sc <- readRDS("/data/combined_named_final_ct.rds")
iterations <- rownames(table(sc$celltype))
Idents(sc) <- "Condition"
table(sc@active.ident)

sc$celltype <- factor(sc$celltype, 
                            levels = c("CD4_naive", "CD4_TEM", "CD4_TCM", "CD8_naive", "CD8_TEM", "CD8_TCM",
                                       "NK", "CTL", "GD_T", "CD8_T", "CD8_prolif", "CD14_mono", "CD16_mono",
                                       "B", "Inter_B", "pDC", "plasma", "platelet"))

sc$celltype_label <- factor(sc$celltype, 
                                  levels = c("CD4_naive", "CD4_TEM", "CD4_TCM", "CD8_naive", "CD8_TEM", "CD8_TCM",
                                             "NK", "CTL", "GD_T", "CD8_T", "CD8_prolif", "CD14_mono", "CD16_mono",
                                             "B", "Inter_B", "pDC", "plasma", "platelet"),
                                  labels = c("CD4 naive", "CD4 TEM", "CD4 TCM", "CD8 naive", "CD8 TEM", "CD8 TCM",
                                             "NK", "Cytotoxic T", "GD T", "CD8 T", "CD8 prolif", "CD14 monocytes", "CD16 monocytes",
                                             "B Cells", "Intermediate B", "pDC", "Plasma", "Platelets"))

sc$Condition <- factor(sc$Condition,
                             levels = c("VS-_OU-", "VS-_OU+", "VS+_OU-",
                                        "VS+_OU-_B","VS+_OU+_B",
                                        "HIV-_OU+", "HIV-_OU-"),
                             labels = c("VS-OUD-", "VS-OUD+", "VS+OUD-",
                                        "VS+OUD- B","VS+OUD+ B",
                                        "HIV-OUD+", "HIV-OUD-")) 
rds_path <- "rds/CellChat/"

#Running CellChat
#VS-OUD-
Idents(combined) <- "Condition"
VSnegOUneg <- subset(combined,idents="VS-OUD-")
Idents(VSnegOUneg ) <- "celltype_label"
drop <-c("GD T","CD8 prolif","CD8 T","pDC", "Plasma", "Platelets")
VSnegOUneg <- subset(VSnegOUneg,idents=drop,invert=T)
cellchat_VSnegOUneg <- createCellChat(object = VSnegOUneg, group.by = "celltype_label")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat_VSnegOUneg@DB <- CellChatDB.use
cellchat_VSnegOUneg <- subsetData(cellchat_VSnegOUneg)
cellchat_VSnegOUneg <- identifyOverExpressedGenes(cellchat_VSnegOUneg)
cellchat_VSnegOUneg <- identifyOverExpressedInteractions(cellchat_VSnegOUneg)
unique(cellchat_VSnegOUneg@idents)
cellchat_VSnegOUneg@idents = droplevels(cellchat_VSnegOUneg@idents, exclude = setdiff(levels(cellchat_VSnegOUneg@idents),unique(cellchat_VSnegOUneg@idents)))
cellchat_VSnegOUneg <- computeCommunProb(cellchat_VSnegOUneg)
cellchat_VSnegOUneg <- filterCommunication(cellchat_VSnegOUneg, min.cells = 10)
cellchat_VSnegOUneg <- computeCommunProbPathway(cellchat_VSnegOUneg)
cellchat_VSnegOUneg <- aggregateNet(cellchat_VSnegOUneg)

#VS+OUD-
Idents(combined) <- "Condition"
VSposOUneg <- subset(combined,idents="VS+OUD-")
Idents(VSposOUneg ) <- "celltype_label"
drop <-c("GD T","CD8 prolif","CD8 T","pDC", "Plasma", "Platelets")
VSposOUneg <- subset(VSposOUneg,idents=drop,invert=T)
cellchat_VSposOUneg <- createCellChat(object = VSposOUneg, group.by = "celltype_label")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat_VSposOUneg@DB <- CellChatDB.use
cellchat_VSposOUneg <- subsetData(cellchat_VSposOUneg)
cellchat_VSposOUneg <- identifyOverExpressedGenes(cellchat_VSposOUneg)
cellchat_VSposOUneg <- identifyOverExpressedInteractions(cellchat_VSposOUneg)
unique(cellchat_VSposOUneg@idents)
cellchat_VSposOUneg@idents = droplevels(cellchat_VSposOUneg@idents, exclude = setdiff(levels(cellchat_VSposOUneg@idents),unique(cellchat_VSposOUneg@idents)))
cellchat_VSposOUneg <- computeCommunProb(cellchat_VSposOUneg)
cellchat_VSposOUneg <- filterCommunication(cellchat_VSposOUneg, min.cells = 10)
cellchat_VSposOUneg <- computeCommunProbPathway(cellchat_VSposOUneg)
cellchat_VSposOUneg <- aggregateNet(cellchat_VSposOUneg)

#VS-OU+
Idents(combined) <- "Condition"
VSnegOUpos <- subset(combined,idents="VS-OUD+")
Idents(VSnegOUpos ) <- "celltype_label"
drop <-c("GD T","CD8 prolif","CD8 T","pDC", "Plasma", "Platelets")
VSnegOUpos <- subset(VSnegOUpos,idents=drop,invert=T)
cellchat_VSnegOUpos <- createCellChat(object = VSnegOUpos, group.by = "celltype_label")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat_VSnegOUpos@DB <- CellChatDB.use
cellchat_VSnegOUpos <- subsetData(cellchat_VSnegOUpos)
cellchat_VSnegOUpos <- identifyOverExpressedGenes(cellchat_VSnegOUpos)
cellchat_VSnegOUpos <- identifyOverExpressedInteractions(cellchat_VSnegOUpos)
unique(cellchat_VSnegOUpos@idents)
cellchat_VSnegOUpos@idents = droplevels(cellchat_VSnegOUpos@idents, exclude = setdiff(levels(cellchat_VSnegOUpos@idents),unique(cellchat_VSnegOUpos@idents)))
cellchat_VSnegOUpos <- computeCommunProb(cellchat_VSnegOUpos)
cellchat_VSnegOUpos <- filterCommunication(cellchat_VSnegOUpos, min.cells = 10)
cellchat_VSnegOUpos <- computeCommunProbPathway(cellchat_VSnegOUpos)
cellchat_VSnegOUpos <- aggregateNet(cellchat_VSnegOUpos)

#VS+OUD-B
Idents(combined) <- "Condition"
VSposOUnegB <- subset(combined,idents="VS+OUD- B")
Idents(VSposOUnegB) <- "celltype_label"
drop <-c("CD8 prolif","CD8 T","pDC", "Plasma", "Platelets")
VSposOUnegB <- subset(VSposOUnegB, idents=drop,invert=T)
cellchat_VSposOUnegB <- createCellChat(object = VSposOUnegB, group.by = "celltype_label")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat_VSposOUnegB@DB <- CellChatDB.use
cellchat_VSposOUnegB <- subsetData(cellchat_VSposOUnegB)
cellchat_VSposOUnegB <- identifyOverExpressedGenes(cellchat_VSposOUnegB)
cellchat_VSposOUnegB <- identifyOverExpressedInteractions(cellchat_VSposOUnegB)
unique(cellchat_VSposOUnegB@idents)
cellchat_VSposOUnegB@idents = droplevels(cellchat_VSposOUnegB@idents, exclude = setdiff(levels(cellchat_VSposOUnegB@idents),unique(cellchat_VSposOUnegB@idents)))
cellchat_VSposOUnegB <- computeCommunProb(cellchat_VSposOUnegB)
cellchat_VSposOUnegB <- filterCommunication(cellchat_VSposOUnegB, min.cells = 10)
cellchat_VSposOUnegB <- computeCommunProbPathway(cellchat_VSposOUnegB)
cellchat_VSposOUnegB <- aggregateNet(cellchat_VSposOUnegB)

#VS+OUD+B
Idents(combined) <- "Condition"
VSposOUposB <- subset(combined,idents="VS+OUD+ B")
Idents(VSposOUposB) <- "celltype_label"
drop <-c("CD8 prolif","CD8 T","pDC", "Plasma", "Platelets")
VSposOUposB <- subset(VSposOUposB, idents=drop,invert=T)
cellchat_VSposOUposB <- createCellChat(object = VSposOUposB, group.by = "celltype_label")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat_VSposOUposB@DB <- CellChatDB.use
cellchat_VSposOUposB <- subsetData(cellchat_VSposOUposB)
cellchat_VSposOUposB <- identifyOverExpressedGenes(cellchat_VSposOUposB)
cellchat_VSposOUposB <- identifyOverExpressedInteractions(cellchat_VSposOUposB)
unique(cellchat_VSposOUposB@idents)
cellchat_VSposOUposB@idents = droplevels(cellchat_VSposOUposB@idents, exclude = setdiff(levels(cellchat_VSposOUposB@idents),unique(cellchat_VSposOUposB@idents)))
cellchat_VSposOUposB <- computeCommunProb(cellchat_VSposOUposB)
cellchat_VSposOUposB<- filterCommunication(cellchat_VSposOUposB, min.cells = 10)
cellchat_VSposOUposB <- computeCommunProbPathway(cellchat_VSposOUposB)
cellchat_VSposOUposB <- aggregateNet(cellchat_VSposOUposB)

#VISUALIZATIONS USING RIVER PLOTS
# Figure 1g
# ###################################
#Plotting river plots for VS-OUD-
selectK(cellchat_VSnegOUneg, pattern = "outgoing")
nPatterns <- 3
cellchat_VSnegOUneg <- identifyCommunicationPatterns(cellchat_VSnegOUneg, pattern = "outgoing", k = nPatterns)
saveRDS(cellchat_VSnegOUneg, file = paste0(rds_path, "cellchat_VSnegOUneg.rds"))

c1 <- netAnalysis_river(cellchat_VSnegOUneg, pattern = "outgoing", 
                        font.size = 6)
c1$layers[[1]]$aes_params$alpha <- 1    
c1$layers[[1]]$geom_params$width <- 1.5  
c1

#Plotting river plots for VS+OUD-
selectK(cellchat_VSposOUneg, pattern = "outgoing")
nPatterns <- 3
cellchat_VSposOUneg <- identifyCommunicationPatterns(cellchat_VSposOUneg, pattern = "outgoing", k = nPatterns)
saveRDS(cellchat_VSposOUneg, file = paste0(rds_path, "cellchat_VSposOUneg.rds"))

c2 <- netAnalysis_river(cellchat_VSposOUneg, pattern = "outgoing", 
                        font.size = 6)
c2$layers[[1]]$aes_params$alpha <- 1    
c2$layers[[1]]$geom_params$width <- 1.5  

c2

#Plotting river plots for VS-OUD+
selectK(cellchat_VSnegOUpos, pattern = "outgoing")
nPatterns <- 3
cellchat_VSnegOUpos <- identifyCommunicationPatterns(cellchat_VSnegOUpos, pattern = "outgoing", k = nPatterns)
saveRDS(cellchat_VSnegOUpos, file = paste0(rds_path, "cellchat_VSnegOUpos.rds"))

c3 <- netAnalysis_river(cellchat_VSnegOUpos, pattern = "outgoing",
                        font.size = 6)
c3$layers[[1]]$aes_params$alpha <- 1    
c3$layers[[1]]$geom_params$width <- 1.5  

c3

#Compare VS-OU- and VS+OU- 
object.list <- list(VSnegOUneg = cellchat_VSnegOUneg, VSposOUneg = cellchat_VSposOUneg)
names(object.list) <- c("VS-OUD-", "VS+OUD-")
cellchat_VS <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
compname <- "cellchat_VS"
comp <- eval(as.symbol(compname))
colors <- c(col_VSNegONeg, col_VSPosONeg)

#VISUALIZATIONS USING RANKNET
# Figure 2e
rankNet(comp, mode = "comparison", stacked = T, do.stat = TRUE, font.size = 10, color.use = colors)+
  theme(text = element_text(face="bold"))+
  theme(aspect.ratio = 5)+
  theme(legend.position = "bottom")

#NETWORK PLOTS for IFN-II Pathway
# Figure 2f
pathway<- "IFN-II"
compname <- "cellchat_VS"
comp <- eval(as.symbol(compname))
colors <- c(col_VSNegONeg, col_VSPosONeg)

pathways.show <- pathway 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, 
                      layout = "circle", edge.weight.max = weight.max[1], 
                      edge.width.max = 10, vertex.label.cex = 1.25,
                      vertex.label.font = 2,
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

#VISUALIZATIONS USING GENE EXPRESSION PLOTS FOR SELPLG
# Figure 2i
pathway<- "SELPLG"
pathways.show <- pathway
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      vertex.label.cex = 1.5,
                      signaling.name = paste(pathways.show, names(object.list)[i]))+
    theme(text = element_text(face="bold"))
}

#VISUALIZATIONS USING GENE EXPRESSION PLOTS FOR SELPLG
# Figure 3j
pathway<- "SELPLG"
pathways.show <- pathway
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), 
                           attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
plotGeneExpression(cellchat_VStoVSOU, signaling = pathway, split.by = "datasets", 
                   colors.ggplot = T, color.use = colors)+ 
  theme(text = element_text(face="bold"),
        legend.position = "right")

#Compare VS-OUD- and VS-OUD+
# Figure 3j
object.list <- list(VSnegOUneg = cellchat_VSnegOUneg, VSnegOUpos = cellchat_VSnegOUpos)
names(object.list)
names(object.list) <- c("VS-OUD-", "VS-OUD+")
cellchat_VStoVSOU <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

#VISUALIZATIONS USING RANKNET
compname <- "cellchat_VStoVSOU"
comp <- eval(as.symbol(compname))
colors <- c(col_VSNegONeg, col_VSNegOPos)

rankNet(comp, mode = "comparison", stacked = T, do.stat = TRUE, font.size = 10, color.use = colors)+
  theme(text = element_text(face="bold"))+
  theme(aspect.ratio = 5)+
  theme(legend.position = "bottom")

#VISUALIZATIONS USING GENE EXPRESSION PLOTS FOR CD99
# Figure 3k
pathway<- "CD99"
pathways.show <- pathway
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
plotGeneExpression(cellchat_VStoVSOU, signaling = pathway, split.by = "datasets", 
                   colors.ggplot = T, color.use = colors)+ 
  theme(text = element_text(face="bold"),
        legend.position = "right")

#VISUALIZATIONS USING NETWORK PLOTS FOR CD99
# Figure 3l
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show,
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]),
                      vertex.label.cex = 1.5,
                      vertex.label.font = 2)+
    theme(text = element_text(face="bold"))
}

#VISUALIZATIONS USING CHORD PLOTS
# Figure 6d
levels <- c("VSnegOUneg", "VSnegOUpos")
object.list <- list(VSnegOUneg = cellchat_VSnegOUneg, VSnegOUpos = cellchat_VSnegOUpos)
names(object.list)
names(object.list) <- c("VS-OUD-", "VS-OUD+")

colors <- c(col_VSNegONeg, col_VSNegOPos)
exp <- cellchat_VSnegOUpos

cellchat_VStoVSOU <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
pos.dataset <- "VS-OUD+"
features.name <- pos.dataset

cellchat_VStoVSOU <- identifyOverExpressedGenes(cellchat_VStoVSOU, 
                                                group.dataset = "datasets", 
                                                pos.dataset = pos.dataset, 
                                                features.name = features.name,
                                                only.pos = FALSE, thresh.pc = 0,
                                                thresh.fc = 0.15, thresh.p = 0.05)

net <- netMappingDEG(cellchat_VStoVSOU, features.name = features.name)

# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat_VStoVSOU,
                              net = net, 
                              datasets = "VS-OUD+",
                              ligand.logFC = 0.2, receptor.logFC = NULL)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat_VStoVSOU, 
                                net = net,
                                datasets = "VS-OUD-",
                                ligand.logFC = -0.2, receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat_VStoVSOU)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_VStoVSOU)

# Chord diagram of all pathways between cell types up and down reg 
plot <- netVisual_chord_gene(object.list[[2]], sources.use = c(9,10), targets.use = c(1:8,11,12), 
                             slot.name = 'net', net = net.up, lab.cex = 3.5, small.gap = 3, 
                             show.legend = TRUE)
title(paste0("Up-regulated signaling in ", names(object.list)[2]), line = -1, adj = 1)

plot1 <- netVisual_chord_gene(object.list[[2]], sources.use = c(9,10), targets.use = c(1:8,11,12),
                              slot.name = 'net', net = net.down, lab.cex = 3.5, small.gap = 3, 
                              show.legend = FALSE)
title(paste0("Down-regulated signaling in ", names(object.list)[2]), line = -1, adj = 1)

# VISUALIZATION USING NETWORK PLOTS IN GRN
#  Figure 6e
pathway<- "GRN"
pathways.show <- pathway
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))+theme(text = element_text(face="bold"))
}

#Compare VS+OUD-_B and VS-OUD+
# Figure 5f
compname <- "cellchat_VStoVSOUB"
comp <- eval(as.symbol(compname))
object.list <- list(VSnegOUpos = cellchat_VSnegOUpos, VSposOUnegB = cellchat_VSposOUnegB)
names(object.list)
names(object.list) <- c("VS-OUD+", "VS+OUD- B")
colors <- c(col_VSNegOPos, col_VSPosONegB)

#VISUALIZATIONS USING RANKNET
rankNet(comp, mode = "comparison", stacked = T, do.stat = TRUE, 
        font.size = 10, color.use = colors)+
  theme(text = element_text(face="bold"))+
  theme(aspect.ratio = 5)+
  theme(legend.position = "bottom")

#Compare VS+OUD+_B and VS-OUD+
# Figure 5f
compname <- "cellchat_VStoOUB"
comp <- eval(as.symbol(compname))
levels <- c("VSnegOUpos", "VSposOUposB")
object.list <- list(VSnegOUpos = cellchat_VSnegOUpos, VSposOUposB = cellchat_VSposOUposB)
names(object.list)
names(object.list) <- c("VS-OUD+", "VS+OUD+ B")
colors <- c(col_VSNegOPos, col_VSPosOPosB)

#VISUALIZATIONS USING RANKNET
rankNet(comp, mode = "comparison", stacked = T, do.stat = TRUE, font.size = 10, color.use = colors)+
  theme(text = element_text(face="bold"))+
  theme(aspect.ratio = 5)+
  theme(legend.position = "bottom")
