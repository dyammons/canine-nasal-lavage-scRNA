#!/usr/bin/Rscript

#load custom functions & packages
source("../../../k9_atlas_scRNA/analysis/scripts/customFunctions_Seuratv5.R")

######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin preprocessing   ######## <<<<<<<<<<<<<<
######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#set output param
outName <- "nasal_lavage_n8_2000_log"

#load in preprocessed data
seu.obj <- readRDS("../output/s3/nasal_lavage_n8.rds") 

#generate viln plots using harmony clusters
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_integrated.harmony", outName = outName,
          outDir = "../output/viln/allCells/", returnViln = T 
         )

#use singleR to ID cells
singleR(seu.obj = seu.obj, clusters = "clusterID_integrated.harmony", reduction = "umap.integrated.harmony", 
        outDir = "../output/singleR/", outName = outName)
