#!/usr/bin/Rscript

#load custom functions & packages
source("../../../k9_atlas_scRNA/analysis/scripts/customFunctions_Seuratv5.R")

######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin preprocessing   ######## <<<<<<<<<<<<<<
######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
outName <- "nasal_lavage_n8"

load10x(din = "../input/", dout = "../output/s1/", outName = outName, testQC = FALSE,
                    nFeature_RNA_high = 9500, nFeature_RNA_low = 200, 
                    nCount_RNA_high = 150000, nCount_RNA_low = 100, 
                    percent.mt_high = 12.5, mt_pattern = "^MT-",
                    nfeatures = 2000,
                    removeDubs = TRUE, removeRBC_pal = FALSE, 
                    pal_feats = NULL, isolatePalRBC = FALSE,
                    featPlots = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                  "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                  "CD4", "MS4A1", "PPBP","HBM")
                   )

#integrate the data using all of the four Seurat v5 integration methods
seu.obj <- integrateData(din = "../output/s1/", dout = "../output/s2/", outName = paste0(outName, "_log"), 
                         runAllMethods = TRUE, indReClus = F)


#use clustree to identify clustering parameters that appear most appropriate
# clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = outName, 
#             test_dims = c(50,45,40), algorithm = 3, prefix = "integrated_snn_res.")

#complete data visualization
for (x in list("integrated.cca", "integrated.harmony", "integrated.joint", "integrated.rcpa")) {
    seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = x, 
                            final.dims = 45, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.1, n.neighbors = 10,
                            prefix = "RNA_snn_res.", assay = "RNA", reduction = x,
                            saveRDS = F, return_obj = T, returnFeats = T,
                            features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                            "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                            "CD4", "MS4A1", "PPBP","HBM")
                          )
}
saveRDS(seu.obj, paste0("../output/s3/", outName,".rds")

###################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end preprocessing   ######## <<<<<<<<<<<<<<
###################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
