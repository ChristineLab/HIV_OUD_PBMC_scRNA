# Code for differential expression analysis using MAST with Age, Sex as co-variates and patient as a random effect
suppressPackageStartupMessages({
library(ggplot2)
library(Seurat)
library(reshape2)
library(MAST)
library(tidyverse)
})

#Define paths for results
de_res_path <- "/results/DE_MAST_RE/"

#Loading Seurat object
combined <- readRDS("/data/combined_named_final_ct.rds")

#Extracting expression data 
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
counts <- combined@assays$RNA$data

#Load the meta_data
meta_data <- combined@meta.data
meta_data <- meta_data %>%
  rownames_to_column("wellKey")

#Changing the + and - in the condition to p and m
meta_data$Condition <- gsub("-", "Neg", meta_data$Condition)
meta_data$Condition <- gsub("\\+", "Pos", meta_data$Condition)
meta_data$Condition <- factor(meta_data$Condition)

#Creating singlecell experiment object 
sc <- SingleCellExperiment(assays = list(counts = counts),
                           colData = meta_data)

#Get list of celltypes
iterations <- rownames(table(combined$celltype))
FCTHRESHOLD <- 0.25

#Filter conditions of interest 
control <-  #Enter condition group, eg: "VSNegOUPos"
condition <- #Enter Control group, eg: "VSNegOUNeg"

#Filtering by conditions of interest
sc_filtered <- sc[, sc$Condition %in% c(condition, control)]

#Prepare data for MAST
sca <- FromMatrix(exprsArray = as.matrix(assay(sc_filtered, "counts")),
                  cData = as.data.frame(colData(sc_filtered)))


#Fitting model using MAST
de_results <- lapply(levels(sc_filtered$celltype), function(celltype){
  #Subset SCA for specific cell type
  sca_sub <- sca[, sca$celltype == celltype]
  
  #Relevel conditions 
  cond <- factor(colData(sca_sub)$Condition)
  cond <- relevel(cond, control)
  colData(sca_sub)$Condition <- cond
  
  #Fit the model
  zlm_fit <- zlm(~ Condition + Age + Sex, sca = sca_sub)

  #Extract DE Results
  de_table <- summary(zlm_fit, doLRT = paste0("Condition", condition))$datatable
  
  fcHurdle <- merge(de_table[contrast==paste0("Condition", condition) & component=='H',.(primerid, `Pr(>Chisq)`)], 
                    de_table[contrast==paste0("Condition", condition) & component=='logFC',
                             .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  write.csv(fcHurdle, file = paste0(de_res_path, celltype, "_", condition, "_vs_", control, ".csv"), row.names = FALSE)
  print(paste0("Finished running DE for:", celltype))
})
