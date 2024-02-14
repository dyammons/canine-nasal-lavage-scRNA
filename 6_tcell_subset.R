#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")

######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   customize variables   ######## <<<<<<<<<<<<<<
######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#set output param
outName <- "240202_lav_tcell_n8_2000_log_cfam"
reduction <- "umap.integrated.harmony"
clusMain <- "clusterID_integrated.harmony"
contrast <- c("Post", "Pre")

######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin preprocessing   ######## <<<<<<<<<<<<<<
######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in the proprocessed data & subset on pop of interest
seu.obj <- readRDS("../output/s3/nasal_lavage_n8_log_canFam_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_cfam.csv", groupBy = "clusterID_integrated.harmony", metaAdd = "majorID")
seu.obj$clusterID  <- seu.obj$clusterID_integrated.harmony 

seu.obj <- subset(seu.obj, subset = majorID == "tcell")

#integrate the data using all of the four Seurat v5 integration methods
seu.obj <- integrateData(seu.obj = seu.obj, dout = "../output/s2/", outName = outName, 
                         runAllMethods = TRUE, normalization.method = "LogNormalize", indReClus = T)

#complete data visualization
for (x in list("integrated.cca", "integrated.harmony", "integrated.joint", "integrated.rcpa")) {
    seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = paste0(outName, "_", x), 
                           final.dims = 30, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.3, n.neighbors = 30,
                           prefix = "RNA_snn_res.", assay = "RNA", reduction = x,
                           saveRDS = F, return_obj = T, returnFeats = T,
                           features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
                          )
}

saveRDS(seu.obj, paste0("../output/s3/", outName,"_S3.rds"))

#low quality cells identified, load object back in, remove cells, reintegrate
seu.obj <- readRDS(paste0("../output/s3/", outName,"_S3.rds"))


#inspect inciting clusters
pi <- DimPlot(seu.obj, 
                reduction = reduction, 
                group.by = "clusterID",
#                 cols = colArray$newCol,
                pt.size = 0.25,
                label = TRUE,
                label.box = TRUE
) + NoLegend()
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black", smallAxes = F) 
ggsave(paste0("../output/", outName, "/", outName, "_clusID_orig_UMAP.png"), width = 7, height = 7)


Idents(seu.obj) <- "clusterID_integrated.harmony"
seu.obj.sub <- subset(seu.obj, invert = T, idents = 5)
table(seu.obj.sub$clusterID_integrated.harmony)
seu.obj.sub$clusterID_sub_orig  <- seu.obj.sub$clusterID_integrated.harmony 

#integrate the data using all of the four Seurat v5 integration methods -- NOTE: Warning regarding too few cells!!!
seu.obj <- integrateData(seu.obj = seu.obj.sub, dout = "../output/s2/", outName = paste0(outName,"_clean"), 
                         runAllMethods = TRUE, normalization.method = "LogNormalize", indReClus = T)

#complete data visualization
for (x in list("integrated.cca", "integrated.harmony", "integrated.joint", "integrated.rcpa")) {
    seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = paste0(outName, "_clean_", x), 
                           final.dims = 30, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.3, n.neighbors = 30,
                           prefix = "RNA_snn_res.", assay = "RNA", reduction = x,
                           saveRDS = F, return_obj = T, returnFeats = T,
                           features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
                          )
}

saveRDS(seu.obj, paste0("../output/s3/",outName,"_clean_S3.rds"))

################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin analysis   ######## <<<<<<<<<<<<<<
################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

seu.obj <- readRDS(paste0("../output/s3/", outName,"_clean_S3.rds"))
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "colz")

#generate viln plots using harmony clusters
vilnPlots(seu.obj = seu.obj, groupBy = clusMain, outName = outName,
          outDir = paste0("../output/viln/", outName, "/"), returnViln = T 
         )

### Export data for interactive cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "../output/cb_input/", 
                markers = paste0("../output/viln/", outName, "/", outName, "_", clusMain,"_gene_list.csv"),
                reduction = reduction,  
                colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", 
                                "majorID", clusMain, "name", "cellSource"), 
                skipEXPR = F, test = F,
                feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                          "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                          "CD4", "MS4A1", "PPBP", "HBM")
)    

#use singleR to ID cells
singleR(seu.obj = seu.obj, clusters = clusMain, reduction = reduction, 
        outDir = "../output/singleR/", outName = outName)

### Plot initial cluster UMAP
pi <- DimPlot(seu.obj, 
                reduction = reduction, 
                group.by = clusMain,
#                 cols = colArray$newCol,
                pt.size = 0.25,
                label = TRUE,
                label.box = TRUE
) + NoLegend()
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black", smallAxes = F) 
ggsave(paste0("../output/", outName, "/", outName, "_rawUMAP.png"), width = 7, height = 7)

### UMAP by sample -- if unequal sample size downsample by cellSource
Idents(seu.obj) <- "orig.ident"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj$orig.ident)))
pi <- DimPlot(seu.obj.ds, 
                reduction = reduction, 
                group.by = "name",
                cols = levels(seu.obj.ds$colz),
                pt.size = 0.25,
                label = FALSE,
                shuffle = TRUE
)
p <- formatUMAP(pi) + labs(colour="Cell source:") + theme(legend.position = "top", 
                                                            legend.direction = "horizontal",
                                                            legend.title=element_text(size=12)
                                                            ) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))
ggsave(paste0("../output/", outName, "/", outName, "_umap_bySample.png"), width =7, height = 7)

### Stacked bar graph by clusMain
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "name", groupBy = "name", clusters = clusMain) + scale_x_discrete(expand = c(0, 0)) +
scale_fill_manual(labels = levels(seu.obj$name), 
                    values = levels(seu.obj$colz)) + 
    theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12))
ggsave(paste0("../output/", outName, "/", outName, "_stackedBar_cluster.png"), width =7, height = 5)

