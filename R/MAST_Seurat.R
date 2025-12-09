#Differential Expression analysis using the MAST Seurat with Age and Sex added as co-variates
suppressPackageStartupMessages({
library(Seurat)
library(ggplot2)
library(tidyverse)
})

#Load the Seurat object
sc <- readRDS("/data/combined_named_final_ct.rds")
iterations <- rownames(table(sc$celltype))

de_res_path <- "results/DE_Seurat_MAST/"
de_rds_path <- "rds/DE_Seurat_MAST/"
conditions <- #Enter condition groups, eg: "VS-OU+" 
controls <- #Enter condition groups, eg: "VS-OU-" 

for (control in controls){
  for(i in iterations){
    for (cond in conditions){
      if (control != cond) {
        tryCatch({
          ct_object <- subset(sc, idents=i)
          ct_object <- SetIdent(ct_object, value="Condition")
          de_ct <- FindMarkers(ct_object, ident.1=cond, ident.2=control,
                               logfc.threshold=0, min.pct=0.10, pseudocount.use = 0.001, 
                               test.use = "MAST", verbose = TRUE, slot = "data",
                               latent.vars = c("Age","Sex"))
  
          saveRDS(de_ct, paste0(de_rds_path,i, "_", cond, "_vs_", control, ".rds"))
          write.csv(de_ct, file = paste0(de_res_path,i, "_", cond,  "_vs_", control, ".csv"),
                    quote=F,
                    row.names=T)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
  }
}


