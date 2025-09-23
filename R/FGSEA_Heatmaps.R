suppressPackageStartupMessages({
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(ggplot2)
library(Seurat)
})

#Load Paths
res_path <- "/results/DE/fgsea/fgsea_output/"
fig_path <- "/figures/DE/fgsea/fgsea_output/"

#GROUPED BY CONDITION
#Load Data
CELLTYPES <- c("CD4_naive","CD4_TEM","CD4_TCM","CD8_naive","CD8_TEM","CD8_TCM","NK","CTL","CD14_mono","CD16_mono","B","Inter_B")
colnames <- CELLTYPES

#CONDTIIONS
#set1
# controls <- c("VS+_OU-")
# conditions <- c("VS-_OU-")
controls <- c("VSPos_OUNeg")
conditions <- c("VSNeg_OUNeg")
comp <- c("VS-OU- v VS+OU-")

#set2
# controls <- c("VS-_OU-")
# conditions <- c("VS-_OU+")
conditions <- "VSNeg_OUPos"
controls <- "VSNeg_OUNeg"
comp <- "VS-OU+ v VS-OU-"

#set3
# controls <- c("VS+_OU-")
# conditions <- c("VS-_OU+")
controls <- "VSPos_OUNeg"
conditions <- "VSNeg_OUPos"
comp <- "VS-OU+ v VS+OU-"

#set4
# controls <- c("VS-_OU+")
# conditions <- c("VS+_OU-_B")
controls <- "VSNeg_OUPos"
conditions <- "VSPos_OUNeg_B"
comp <- "VS+OU-_B v VS-OU+"

#set5
# controls <- c("VS-_OU+")
# conditions <- c("VS+_OU+_B")
controls <- "VSNeg_OUPos"
conditions <- "VSPos_OUPos_B"
comp <- "VS+OU+_B v VS-OU+"

#set5
# controls <- c("HIV-_OU-")
# conditions <- c("VS+_OU-")
controls <- "HIVNeg_OUNeg"
conditions <- "VSPos_OUNeg"
comp <- "VS+OU- v HIV-OU-"

gs <- character()
for (control in controls){
  for (cond in conditions){
    if (control != cond) {
      for (ct in CELLTYPES){
        tryCatch({
          t <- read.csv(paste0(res_path,ct,"_",cond,"_vs_",control,"_fgseasig_n5.csv"))       
          gs <- c(gs,t$pathway)
          
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
      }}}}

gs <- unique(gs)
gs <- na.omit(gs)
df <- data.frame(row.names=gs)

for (control in controls){
  for (cond in conditions){
    if (control != cond) {
      for (ct in CELLTYPES){
        tryCatch({
          t <- read.csv(paste0(res_path,ct,"_",cond,"_vs_",control,"_fgseaRes.csv"), sep ="\t")
          t2 <- t[ , c("pathway", "NES")]
          rownames(t2) <- t2$pathway  
          t2$pathway <- NULL
          
          df <- merge(df, t2, by="row.names", all.x=T)
          
          colnames(df)[ncol(df)] <- paste0(control,"_",ct,"_",cond)
          rownames(df) <- df$Row.names 
          df$Row.names <- NULL
          df[is.na(df)] <- 0      
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }}}}

#shorten rownames
row.names(df) <- gsub("_", " ", row.names(df), fixed=TRUE)
row.names(df) <- strtrim(row.names(df), width=60)

#change rowname capitalization
row.names(df) <- str_to_sentence(row.names(df))

#set split
split <- rep(1:length(CELLTYPES), each = 1)
ha <- HeatmapAnnotation(foo = anno_empty(border = F, height = unit(0.1, "cm")))
col_fun <- colorRamp2(c(-6,-2,-1, -.5, 0,.5,1,2,6), c("#313695", "#4575B4","#74ADD1","#ABD9E9","#FFFFBF", "#FDAE61", "#F46D43", "#D73027","#A50026"))
df <- as.matrix(df)

#make heatmap
ht <- Heatmap(df,
              name = "NES",
              col = col_fun,
              column_names_rot=45,
              column_names_centered = TRUE,
              width = ncol(df)*unit(8, "mm"), 
              height = nrow(df)*unit(3, "mm"),
              row_names_gp = gpar(fontsize = 8, fontface = "bold"),
              column_names_side = "top",
              column_names_gp = gpar(fontsize = 8, fontface = "bold"),
              cluster_columns=F,
              cluster_rows=T,
              show_row_dend = F,
              show_column_dend = F,
              column_title = comp,
              column_labels = colnames,
              # column_split = split,
              top_annotation = ha,
              column_title_gp = gpar(fontsize=5, fontface="bold"),
              heatmap_legend_param = list(direction = "horizontal",  
                                          labels_gp = gpar(fontsize =8, fontface="bold"), 
                                          legend_width = unit(4, "cm"),
                                          title_gp = gpar(fontsize = 8, fontface="bold")),
              row_names_max_width = max_text_width( rownames(df), 
                                                    gp = gpar(fontsize = 1, face='bold')))

ht <- ComplexHeatmap::draw(ht, heatmap_legend_side="bottom")

pdf(paste0(fig_path,comp,"_n5.pdf"), width = 10, height = 20)
print(ht)
dev.off()
