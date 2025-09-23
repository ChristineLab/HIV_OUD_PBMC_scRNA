suppressPackageStartupMessages({
  library(circlize)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(tidyverse)
  library(ggplot2)
  library(Seurat)
})

#Load Data
combined <- readRDS("/data/combined_named_final_ct.rds")

CELLTYPES <- c("CD4_naive","CD4_TEM","CD4_TCM","CD8_naive","CD8_TEM","CD8_TCM","NK","CTL","CD14_mono","CD16_mono","B","Inter_B")
col_fun<-colorRamp2(c(-6,-2,-1, -.5, 0,.5,1,1.5,2), c("#313695", "#4575B4","#74ADD1","#ABD9E9","#FFFFBF", "#FDAE61", "#F46D43", "#D73027","#A50026"))

#DE genes overlap with leading edges 

#FIG 2h HALLMARK INTERFERON GAMMA RESPONSE CD4 NAIVE VS+OU- v VS-OU-
# genes <- c("HLA-DRB1","TNFSF10","EPSTI1","XAF1","IFI44L","IFI44","CCL5","GZMA","ISG15","TRIM14","PARP12","SOCS1","TAP1","IFI35","DDX60","UBE2L6","PFKP","MX1","OAS2","DDX58","PIM1","STAT1","SPPL2A","SAMD9L")
# 
# genes_fig <- c("HLA-DRB1", "TNFSF10", "XAF1", "IFI44L", "ISG15", "SOCS1", "TAP1", "UBE2L6", "PIM1", "STAT1", "SAMD9L")
# 
# path <- "HALLMARK INTERFERON GAMMA RESPONSE"
# conds <- c("VS-_OU-","VS+_OU-")
# controls <- c("VS-_OU-")
# conditions <- c("VS+_OU-")
# cells <- "CD4_naive"
# colnames <- c("OU-","OU-")
# colwidth <- 5
# heights <- 10
# fig <- "2h"

# #FIG 2i GO IMMUNE PROCESS CD4 NAIVE VS+OU- v VS-OU- using this instead
# genes <- c("HLA-DRB1","HLA-DQB1","HLA-DPA1","LGALS1","TIGIT","HLA-DPB1","CD8B","WIPF2","IGHA1","CST7","HLA-DRA","SIT1","CCL5","ISG15","GBP1","SPN","PYCARD","BTN2A1","SOCS1","DDX60","LDB1","GBP5","NFKBIZ","LGALS9","CD8A","DDX58","TNFRSF1B","IKBKE","PIM1","STAT1","APOBEC3G","TOB2","PARP9","SPPL2A","BTN3A2","CD59","CD84","PRDM1","UNC13D","ITGB1","PSMA5","BAG6","KLHL6","PSMC4","JUN","ADCY7","CCDC88B","CTR9","FBXW11","IRF3","MYO1G","GATA3","SH2D1A","TAB2","HEXIM1","ITGAL","ZFP36L1","TFRC","ITGA4","HERC5","PSMB8","GPR171","RFTN1","TYROBP","MYC","RASAL3","CD74","PSMD3","PSMB9","STK11","NPLOC4","CTSC","PTPN6","SOCS3","HLA-A","CUEDC2","PSMD8","PSMB3","ID2","TRAC")
# # 
# genes_fig <- c("PSMB3", "LGALS1", "HLA-DRA", "HLA-DRB1", "HLA-DQB1", "SOCS1", "STAT1", "JUN", "IRF3", "HLA-A", "PYCARD", "HLA-DPB1", "ZFP36L1", "PSMA5", "TRAC", "GBP1", "TIGIT", "PSMB8", "NFKBIZ", "HEXIM1")
# #  
# path <- "GO IMMUNE PROCESS"
# conds <- c("VS-_OU-","VS+_OU-")
# controls <- c("VS-_OU-")
# conditions <- c("VS+_OU-")
# cells <- "CD4_naive"
# colnames <- c("OU-","OU-")
# colwidth <- 5
# heights <- 10
# fig <- "2i"

# #SUPPL FIG 2 GO RESPONSE TO CYTOKINE CD4 NAIVE VS+OU- v VS-OU-
# genes <- c("HLA-DRB1","HLA-DQB1","HLA-DPA1","ANXA2","HLA-DPB1","XAF1","MAP3K5","HLA-DRA","CCL5","TCIRG1","ISG15","GBP1","CARD16", "PYCARD","SOCS1","IFI35","GBP5","TRIM56","PFKP","TALDO1","MX1","LGALS9","OAS2","DDX58","TNFRSF1B", "SELPLG","IKBKE","PIM1","ATIC","STAT1","CALCOCO2","PARP9","SPPL2A","IL32","IFI6","YBX3","OTULIN","ACTN4", "NFYB","ITGB1","PSMA5","TRADD","IL2RB","PSMC4","CAMK2G","CTR9","FBXW11","IRF2","IRF3","GATA3","TAB2","GBP4", "ZFP36L1","TFRC","CPNE1","MAT2A","ITGA4","PTPN4","PSMB8","IL10RA","MYC","CD74","ZYX","PSMD3","TANK","PSMB9", "PTPN6","MAPKAPK2","SOCS3","HLA-A","PSMD8","PSMB3")
# 
# genes_fig <- c("HLA-DRB1", "HLA-DQB1", "HLA-DP1", "ANXA2", "HLA-DPB1", "XAF1", "HLA-DRA", "ISG15", "GBP1", "CARD16", "PYCARD", "SOCS1", "TALDO1", "PIM1", "STAT1", "PSMA5", "IRF3", "GBP4", "ZFP36L1", "PSMB8", "IL10RA", "PSMB9", "HLA-A", "PSMB3")
# 
# genes_final <- genes[genes%in%genes_fig]
# path <- "GO RESPONSE TO CYTOKINE"
# conds <- c("VS-_OU-", "VS+_OU-")
# controls <- c("VS-_OU-")
# conditions <- c("VS+_OU-")
# cells <- "CD4_naive"
# colnames <- c("OU-", "OU-")
# colwidth <- 5
# heights <- 10
# # fig <- "S4b"
# 
# # #SUPPL FIG 2 GO RESPONSE TO CYTOKINE CD8 NAIVE VS+OU- v VS-OU-
# genes <- c("HLA-DRB1","HLA-DQB1","HLA-DPA1","ANXA2","HLA-DPB1","XAF1","MAP3K5","HLA-DRA","CCL5","TCIRG1","ISG15","GBP1","CARD16", "PYCARD","SOCS1","IFI35","GBP5","TRIM56","PFKP","TALDO1","MX1","LGALS9","OAS2","DDX58","TNFRSF1B", "SELPLG","IKBKE","PIM1","ATIC","STAT1","CALCOCO2","PARP9","SPPL2A","IL32","IFI6","YBX3","OTULIN","ACTN4", "NFYB","ITGB1","PSMA5","TRADD","IL2RB","PSMC4","CAMK2G","CTR9","FBXW11","IRF2","IRF3","GATA3","TAB2","GBP4", "ZFP36L1","TFRC","CPNE1","MAT2A","ITGA4","PTPN4","PSMB8","IL10RA","MYC","CD74","ZYX","PSMD3","TANK","PSMB9", "PTPN6","MAPKAPK2","SOCS3","HLA-A","PSMD8","PSMB3")
#  
# genes_fig <- c("CISH", "HLA-DRB1","TNFRS1B","GBP1","RIPK2","GSK3A","XAF1","PSMA5","PFKP","BAD","STIP1","HLA-DPB1","PIM1","IFNGR1", "IL2RB", "SOCS3", "PSMC4","DDX58","TRADD","PYCARD","TP53","PHB","SELPLG","PSMB10","PARP9","GSTO1","SOCS1","LDLRAP1", "ADAR","PSMD13","STAT1","JAK3", "PTPN7","STAT2","HDGF","CEBPB","PSMC1","TNIP2","IFI16","P4HB","SMAD4", "MX1","HLA-DP1","IFI35","BST2")
# 
# genes_final <- genes[genes%in%genes_fig]
# path <- "GO RESPONSE TO CYTOKINE"
# conds <- c("VS-_OU-", "VS+_OU-")
# controls <- c("VS-_OU-")
# conditions <- c("VS+_OU-")
# cells <- "CD8_naive"
# colnames <- c("OU-", "OU-")
# colwidth <- 5
# heights <- 10
# fig <- "S4c"

# #FIG3 HALLMARK TNFA CD4 TCM VS+OU- v VS-OU-
# genes <- c("IL7R","TNF","TNFSF9","AREG","NFKBIA","GADD45B","TGIF1","JUNB","PLK2","YRDC","PFKFB3","DUSP1","GPR183","SOCS3","FOS","IER2","ZC3H12A","CEBPB","PNRC1","KDM6B","EGR1","NR4A1","RHOB","IFIH1","CDKN1A","CD83","DDX58","SAT1","NAMPT","ATF3","DUSP5","STAT5A","MYC","DUSP2","PDLIM5","TAP1","IFIT2","SGK1","CCL5")
# 
# genes_fig <- c("SOCS3", "RHOB", "SAT1", "NAMPT", "DUSP2")
# path <- "HALLMARK TNFA SIGNALING VIA NFKB"
# cells <- "CD4_TCM"
# conds <- c("VS-_OU-","VS-_OU+","VS+_OU-")
# colnames <- c("OU-","OU+","OU-")
# colwidth <- 7
# heights <- 5
# fig <- "3f_1"
# controls <- c("VS-_OU+")
# conditions <- c("VS+_OU-")

# #FIG3 HALLMARK TNFA CD8 TCM VS-OU+ v VS+ OU-
genes <- c("IL7R","TNF","TNFSF9","AREG","NFKBIA","GADD45B","TGIF1","JUNB","PLK2","YRDC","PFKFB3","DUSP1","GPR183","SOCS3","FOS","IER2","ZC3H12A","CEBPB","PNRC1","KDM6B","EGR1","NR4A1","RHOB","IFIH1","CDKN1A","CD83","DDX58","SAT1","NAMPT","ATF3","DUSP5","STAT5A","MYC","DUSP2","PDLIM5","TAP1","IFIT2","SGK1","CCL5")
path <- "HALLMARK TNFA SIGNALING VIA NFKB"
cells <- "CD8_TCM"
conds <- c("VS-_OU-","VS-_OU+","VS+_OU-")
colnames <- c("OU-","OU+","OU-")
colwidth <- 5
heights <- 3
fig <- "3"
controls <- c("VS-_OU+")
conditions <- c("VS+_OU-")

#DEG ONLY!
DE_genes <- c()

for (control in controls){
  for (cond in conditions){
    for(ct in cells){
      if (control != cond) {
        tryCatch({
          DEG <- read.csv(paste0("/results/DE/", ct, "_", cond, "_vs_", controls,  "_SIG.csv"), row.names = 1)
          
          DE_genes <- c(DE_genes,DEG$name)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }}}}

DE_genes <- unique(DE_genes)
genes <- Reduce(intersect,list(genes, DE_genes))

genes_final <- genes[genes%in%genes_fig]

combined <- SetIdent(combined,value="Condition")
sub_cond <- conds
temp <- subset(combined, idents=sub_cond)

temp <-SetIdent(temp, value="celltype")
celltype <-cells
temp <- subset(temp, idents=celltype)

genes_final <- genes_final[genes_final %in% rownames(temp@assays$RNA$counts)]

temp <- SetIdent(temp, value="Condition")
temp <- FindVariableFeatures(temp)
temp <- ScaleData(temp,vars.to.regress="nCount_RNA",features=unique(c(temp@assays$RNA@var.features,genes)))
heatcolor <- rev(brewer.pal(11,"RdYlBu"))
what <- AverageExpression(temp,features=genes_final,verbose=T,slot="scale.data")$RNA
min <- quantile(unlist(what),0.1)
mid <- quantile(unlist(what),0.5)
max <-quantile(unlist(what),0.9)
what <- as.matrix(what)
colnames(what) <- conds

p <- Heatmap(what,
             width = unit(colwidth, "cm"),
             name=celltype,
             col=heatcolor,
             clustering_method_columns = "ward.D2",
             cluster_rows=FALSE,
             cluster_columns=FALSE,
             column_split=do.call("rbind",strsplit(colnames(what),split="_"))[,1],
             show_column_names=T,
             show_row_names=T,
             rect_gp = gpar(col = "black", lwd = 2), 
             show_column_dend = FALSE, 
             show_row_dend = FALSE, 
             column_labels = colnames, 
             column_names_rot = 0, 
             column_names_centered = T, heatmap_legend_param = list(title="",direction = "horizontal"),
             column_title_gp = gpar(fontface="bold"), 
             column_names_gp = gpar(fontface="bold"), 
             row_names_gp = gpar(fontface="bold"))

ht_list <- p
p <- draw(p, column_title = paste0(path,"\n",cells),column_title_gp = gpar(fontface = "bold"), heatmap_legend_side ="bottom")


pdf(paste0(de_res_path,"/GO/",cells,"_fgsea_le_heatmap_fig",fig,".pdf"),width = 8, height = heights)
p
dev.off()
