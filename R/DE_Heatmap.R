suppressPackageStartupMessages({
  library(circlize)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(tidyverse)
  library(ggplot2)
  library(Seurat)
})

#Set paths
de_res_path <- "/results/DE_MAST_RE/"
de_rds_path <- "/rds/DE_MAST_RE/"
de_fig_path <- "/figures/DE_MAST_RE/"

col_fun <-colorRamp2(c(-2,-1.5,-1, -.5, 0,.5,1,1.5,2), c("#313695", "#4575B4","#74ADD1","#ABD9E9","#FFFFBF", "#FDAE61", "#F46D43", "#D73027","#A50026"))
celltypes <- c("CD4_naive","CD4_TEM","CD4_TCM","CD8_naive","CD8_TEM","CD8_TCM","NK","CTL","CD14_mono","CD16_mono","B","Inter_B")

#CONDITIONS
for (CTS in celltypes){
  controls <- c("VS-_OU+", "VS-_OU-" )
  conditions <- c("VS-_OU-","VS+_OU-")
  direction <- c("UP","DOWN")
  
  gs <- character()
  
  for (control in controls){
    for (cond in conditions){
      for (dir in direction){
        if (control != cond) {
          tryCatch({
            #selecting top genes 
            data1 <- read.csv(paste0(de_res_path, ct, "_", cond,"_vs_",control,"_SIG_",dir,"_DEG.csv"))
            #data2 <- head(data1, 30)
            data1 <- data1[order(abs(data1$avg_log2FC), decreasing = TRUE),]
            data2 <- head(data1, 10)
            data3 <- data2$name
            gs <-c(gs,data3)
          }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }}}}
  
  genes <- unique(gs)
  
  #sest up complexheatmap dataframe 
  gs <- genes
  df <-data.frame(row.names=gs)
  
  for (control in controls){
    for (cond in conditions){
      if (control != cond) {
        tryCatch({
          
          res <- rep(0,length(gs))
          ressig <- rep(0,length(gs))
          data <- read.csv(paste0(de_res_path, ct, "_", cond,"_vs_",control,".csv"))
          rownames(data) <- data$X
          res[gs %in% data$X] <- data[gs[gs %in% data$X],"avg_log2FC"]
          df <- cbind(df,res)
          colnames(df)[ncol(df)] <- paste0(control," v\n ",cond)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }}}
  
  #set heatmap formatting 
  ht_opt(legend_border = "black",heatmap_border = TRUE)
  colwidth<- 4.5
  rowwidth<- 10
  docheight <-6
  docwidth <- 4
  
  ha <- HeatmapAnnotation(foo = anno_empty(border = F, height = unit(0.25, "cm")))
  ht <- Heatmap(df,width = unit(colwidth, "cm"), height = unit(rowwidth, "cm"), col = col_fun,row_names_gp = gpar(fontsize = 6, face="bold"),column_names_gp = gpar(fontsize=6, face="bold"), cluster_columns = F, cluster_rows = F, column_names_rot=0, column_names_centered = T, column_names_side = "top",name = "log2FC", column_title = CTS, column_title_gp = gpar(fontsize= 8, face="bold"), top_annotation = ha, heatmap_legend_param = list(direction = "horizontal"))
  ht <- draw(ht, heatmap_legend_side = "bottom")
  
  print(ht)
}
