suppressPackageStartupMessages({
library(fgsea)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(magrittr)
library(data.table)
library(ggpubr)
})

#Set file paths
de_res_path <- "results/DE/fgsea/input_csv/"
fgsea_res_path <- "/results/DE/")

#Load data
combined <- readRDS("/data/combined_named_final_ct.rds")

#Load gmt
pathways <- gmtPathways("/data/hall_reactome_pid_go_kegg_st.gmt")

#Comparison Set 1
#controls <- c("VS+_OU-")
controls <- c("VSPos_OUNeg")
#conditions <- c("VS-_OU-")
conditions <- c("VSNeg_OUNeg")

#Comparison Set 2
#controls <- c("VS+_OU-", "VS-_OU-")
controls <- c("VSPos_OUNeg", "VSNeg_OUNeg")
#conditions <- c("VS-_OU+")
conditions <- c("VSNeg_OUPos")

#Comparison Set 3
#controls <- c("VS-_OU+")
controls <- c("VSNeg_OUPos")
#conditions <- c("VS+_OU-_B", "VS+_OU+_B")
conditions <- c("VSPos_OUNeg_B", "VSPos_OUPos_B")
# 
#Comparison Set 4
#controls <- c("HIV-_OU-")
controls <- c("HIVNeg_OUNeg")
#conditions <- c("VS+_OU-")
conditions <- c("VSPos_OUNeg")

celltypes <- c("CD4_naive", "CD4_TEM", "CD4_TCM", "CD8_naive", "CD8_TEM", "CD8_TCM", "NK","CTL","CD14_mono","CD16_mono",
               "B","Inter_B",  "plasma", "CD4_prolif", "CD8_T")

for (control in controls){
  for(ct in celltypes){
    for (cond in conditions){
      if (control != cond) {
        tryCatch({
          #rnk
          response <- read.csv(paste0(de_res_path, ct, "_", cond,"_vs_",control,".csv"), row.names = 1)
          response <- response[order(response$avg_log2FC),]
          b <- cbind(rownames(response),response$avg_log2FC)
          write.table(b,paste0(fgsea_res_path, "/fgsea/fgsea_output/", ct, "_", cond,"_vs_",control,"_GSEA.rnk"),
                      quote=F, col.names=F, sep="\t", row.names=F)
          
          ranks <- b
          ranks <- setNames(as.numeric(ranks[,2]), ranks[,1])
          
          #run fgsea
          fgseaRes <- fgsea(pathways, ranks, eps=0,minSize=5, maxSize=500 )
          
          #save
          fwrite(fgseaRes[order(pval), ], paste0(fgsea_res_path, "/fgsea/fgsea_output/", ct, "_", cond,"_vs_",control, "_fgseaRes.csv"),
                 sep="\t", sep2=c("", " ", ""))
          
          fgseaRes$leadingEdge <- NULL
          
          #plot
          topPathwaysUp <- fgseaRes[ES > 0][head(order(padj), n=10), pathway]
          topPathwaysDown <- fgseaRes[ES < 0][head(order(padj), n=10), pathway]
          topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
          
          
          png(paste0(de_res_path, "/fgsea/fgsea_output/", ct, "_", cond,"_vs_",control,"_top10paths.png"),width = 1500, height = 480)
          plotGseaTable(pathways[topPathways], ranks, fgseaRes,
                        gseaParam=0.5)
          dev.off()
          
          #writecsv_split
          PathwaysUp <- fgseaRes[ES > 0][order(padj)]
          PathwaysUp <- PathwaysUp[padj < 0.05]
          PathwaysDown <- fgseaRes[ES < 0][order(padj)]
          PathwaysDown <- PathwaysDown[padj < 0.05]
          
          write.csv(PathwaysUp, paste0(fgsea_res_path, "/fgsea/fgsea_output/", ct, "_", cond,"_vs_",control,"_fgseaRes_UP.csv"))
          write.csv(PathwaysDown, paste0(fgsea_res_path, "/fgsea/fgsea_output/", ct, "_", cond,"_vs_",control,"_fgseaRes_DOWN.csv"))
          
          #write enrichmentmap input
          fgsea_results <- fgseaRes[order(fgseaRes$padj, decreasing=FALSE),]
          
          GSEA_UP <- subset(fgsea_results, NES > 0)
          GSEA_DN <- subset(fgsea_results, NES < 0)
          
          GSEA_UP_trimmed <- GSEA_UP[,c("pathway","pval","padj")]
          colnames(GSEA_UP_trimmed) <- c("Description", "p.Val", "FDR")
          GSEA_UP_trimmed$GO.ID <- GSEA_UP_trimmed$Description
          GSEA_UP_trimmed$Phenotype = "+1"
          GSEA_UP_trimmed <- GSEA_UP_trimmed[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype")]
          head(GSEA_UP_trimmed)
          write.table(GSEA_UP_trimmed, paste0(fgsea_res_path, "/fgsea/fgsea_output/", ct, "_", cond,"_vs_",control,"_GSEA_UP_for_enrichment.txt"), sep = "\t", quote = F, row.names = F)
          
          GSEA_DN_trimmed <- GSEA_DN[,c("pathway","pval","padj")]
          colnames(GSEA_DN_trimmed) <- c("Description", "p.Val", "FDR")
          GSEA_DN_trimmed$GO.ID <- GSEA_DN_trimmed$Description
          GSEA_DN_trimmed$Phenotype = "-1"
          GSEA_DN_trimmed <- GSEA_DN_trimmed[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype")]
          head(GSEA_DN_trimmed)
          write.table(GSEA_DN_trimmed, paste0(fgsea_res_path, "/fgsea/fgsea_output/", ct, "_", cond,"_vs_",control,"_GSEA_DN_for_enrichment.txt"), sep = "\t", quote = F, row.names = F)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
  }
}

#Saving top 10 significant pathways for each celltype 
numba <- 5

for (control in controls){
  for(ct in celltypes){
    for (cond in conditions){
      if (control != cond) {
        tryCatch({
          #sig only 
          gsea <- read.csv(paste0(fgsea_res_path, "/fgsea/fgsea_output/",ct,"_",cond,"_vs_", control, "_fgseaRes.csv"), sep =  "\t")
          data2 <- gsea[gsea$padj < 0.05,]
          data2 <- data2[, c('pathway', 'NES')]
          data2 <- data2[order(data2$NES),]
          topPathwaysUp <- tail(data2,n=numba)
          topPathwaysDown <- head(data2,n=numba)
          topPathways <- rbind(topPathwaysUp, topPathwaysDown)
          topPathways <- unique(topPathways)
          topPathways <- topPathways[order(topPathways$NES),]
          write.csv(topPathways, paste0(fgsea_res_path, "/fgsea/fgsea_output/",ct,"_",cond,"_vs_", control, "_fgseasig_n5.csv"))
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
  }
}
