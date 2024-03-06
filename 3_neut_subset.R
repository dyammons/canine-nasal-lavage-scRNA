#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")

#set output param
outName <- "240202_lav_neut_n8_2000_log_cfam"
reduction <- "umap.integrated.harmony"
clusMain <- "clusterID_integrated.harmony"
contrast <- c("Post", "Pre")

######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin preprocessing   ######## <<<<<<<<<<<<<<
######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in the proprocessed data & subset on pop of interest
seu.obj <- readRDS("../output/s3/nasal_lavage_n8_log_canFam_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_cfam.csv", groupBy = "clusterID_integrated.harmony", metaAdd = "majorID")

seu.obj <- subset(seu.obj, subset = majorID == "neut")

#integrate the data using all of the four Seurat v5 integration methods
seu.obj <- integrateData(seu.obj = seu.obj, dout = "../output/s2/", outName = outName, 
                         runAllMethods = TRUE, normalization.method = "LogNormalize", indReClus = T)

#complete data visualization
for (x in list("integrated.cca", "integrated.harmony", "integrated.joint", "integrated.rcpa")) {
    seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = paste0(outName, "_", x), 
                           final.dims = 30, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.1, n.neighbors = 10,
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


#inspect clusters
pi <- DimPlot(seu.obj, 
                reduction = reduction, 
                group.by = "clusterID_integrated.harmony",
#                 cols = colArray$newCol,
                pt.size = 0.25,
                label = TRUE,
                label.box = TRUE
) + NoLegend()
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black", smallAxes = F) 
ggsave(paste0("../output/", outName, "/", outName, "_clusID_orig_UMAP.png"), width = 7, height = 7)


#remove clusters 8& 9 
Idents(seu.obj) <- "clusterID_integrated.harmony"
seu.obj.sub <- subset(seu.obj, invert = T, idents = c(8,9))
table(seu.obj.sub$clusterID_integrated.harmony)

#integrate the data using all of the four Seurat v5 integration methods
seu.obj <- integrateData(seu.obj = seu.obj.sub, dout = "../output/s2/", outName = paste0(outName,"_clean"), 
                         runAllMethods = TRUE, normalization.method = "LogNormalize", indReClus = T)

#complete data visualization
for (x in list("integrated.cca", "integrated.harmony", "integrated.joint", "integrated.rcpa")) {
    seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = paste0(outName, "_clean_", x), 
                           final.dims = 30, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.1, n.neighbors = 10,
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

#set output param
outName <- "240202_lav_neut_n8_2000_log_cfam_clean"
reduction <- "umap.integrated.harmony"
clusMain <- "clusterID_integrated.harmony"
contrast <- c("Post", "Pre")

#load in data and metadata
seu.obj <- readRDS(paste0("../output/s3/",outName,"_S3.rds"))
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "cellSource")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "cellSource2")
seu.obj$cellSource2 <- factor(seu.obj$cellSource2, levels = c("d0","d14","d90"))


### Check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, reduction = reduction, nrow = 1, ncol = 3, features = features, 
                    color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste0("../output/", outName, "/", outName, "_QC_feats.png"), width = 9, height = 3)


### Generate dot plots using vilnPlots of harmony clusters
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


### Use singleR to ID cells
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
                                                            ) + guides(colour = guide_legend(nrow = 2, override.aes = list(size = 4)))
ggsave(paste0("../output/", outName, "/", outName, "_umap_bySample.png"), width = 7, height = 7.5)


### Stacked bar graph by clusMain
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "name", groupBy = "name", clusters = clusMain) + scale_x_discrete(expand = c(0, 0)) +
scale_fill_manual(labels = levels(seu.obj$name), 
                    values = levels(seu.obj$colz)) + 
    theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12))
ggsave(paste0("../output/", outName, "/", outName, "_stackedBar_cluster.png"), width =7, height = 5)


### Frequency plots to run stats - cluster
freqy <- freqPlots(seu.obj, method = 1, nrow = 2, 
                   comp = "cellSource", groupBy = clusMain, legTitle = "Cell source", refVal = "name",
                   namez = "name", 
                   colz = "colz"
) + NoLegend()
ggsave(paste0("../output/", outName, "/",outName, "_freqPlots_cluster.png"), width = 6, height = 6)


#run paired analysis -- TO DO: switch to RM ANOVA
df <- freqy$data %>% separate(name, into = c(NA, "dog"))

lapply(unique(df$clusterID_integrated.harmony), function(x){
    df.sub <- filter(df, clusterID_integrated.harmony == x)
    wilcox.test(freq ~ cellSource, data = df.sub)
})


### Complete pseudobulk DGE by all cells
createPB(seu.obj = seu.obj, groupBy = "allCells", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
         clusters = NULL, outDir = paste0("../output/", outName, "/pseudoBulk/")
)

#add dog metadata to account for paired nature of the data
df <- read.csv(paste0("../output/", outName, "/pseudoBulk/allCells_deg_metaData.csv"), row.names = 1)
df$dog <- unlist(lapply(df$sampleID, function(x){strsplit(x, "_")[[1]][2]}))
df$tp <- c(1,1,1,1,3,3,2,2)
write.csv(df, paste0("../output/", outName, "/pseudoBulk/allCells_deg_metaData.csv"))

pseudoDEG(metaPWD = paste0("../output/", outName, "/pseudoBulk/allCells_deg_metaData.csv"),
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("../output/", outName, "/pseudoBulk/"), 
          outName = outName, 
          paired = T, pairBy = "dog", test.use = "Wald", strict_lfc = F,
          idents.1_NAME = contrast[1], idents.2_NAME = contrast[2],
          inDir = paste0("../output/", outName, "/pseudoBulk/"), title = "All cells", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)



### Complete GSEA using the linDEG results
df <- read.csv(paste0("../output/", outName, "/pseudoBulk/allCells/", outName, "_cluster_allCells_all_genes.csv")) %>% arrange(padj)

upGenes <- df %>% filter(log2FoldChange > 0) %>% pull(gene)
dwnGenes <- df %>% filter(log2FoldChange < 0) %>% pull(gene)
p <- plotGSEA(geneList = upGenes, geneListDwn = dwnGenes, category = "C5", subcategory = "GO:BP", 
              upCol = "red", dwnCol = "blue", size = 4)

minVal <- -15
maxVal <- 5
pi <- p + scale_x_continuous(limits = c(minVal, maxVal), name = "Signed log10(padj)")
ggsave(paste("../output/", outName, "/", outName, "_allCells_gsea.png", sep = ""), width = 7, height = 4)


top.degs <- read.csv(paste0("../output/", outName, "/pseudoBulk/allCells/", outName, "_cluster_allCells_all_genes.csv")) %>% filter(log2FoldChange > 0) %>% arrange(padj) %>% pull(gene)

### Plot top DEGS
features <- top.degs[c(1,2,4:6,8,9,11)]
set.seed(12)
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(seu.obj, downsample = min(table(seu.obj$cellSource)))
p <- FeaturePlot(seu.obj.sub,features = features, 
                 reduction = reduction, pt.size = 0.01, 
                 split.by = "cellSource", order = T, 
                 by.col = F) + labs(x = "UMAP1", y = "UMAP2") & theme(axis.text = element_blank(),
                                                                      axis.title.y.right = element_text(size = 16),
                                                                      axis.ticks = element_blank(),
                                                                      axis.title = element_blank(),
                                                                      axis.line = element_blank(),
                                                                      plot.title = element_text(size=16),
                                                                      title = element_blank(),
                                                                      plot.margin = unit(c(1, 0, 0, 0), "pt")
                                                                      ) & scale_color_gradient(breaks = pretty_breaks(n = 3), limits = c(NA, NA), low = "lightgrey", high = "darkblue")
ggsave(paste0("../output/", outName, "/", outName, "_split_degs.png"), width = 16, height = 4)



########## BONUS ##########

### Investigate if d14 or d90 is driving DEGs
Idents(seu.obj) <- "name"
seu.obj.sub <- subset(seu.obj, idents = c("pre_1", "pre_2", "post_1", "post_2"))
table(seu.obj.sub$name)

seu.obj.sub$d90 <- "d90"
seu.obj.sub$d90 <- as.factor(seu.obj$d90)
### Complete pseudobulk DGE by all cells
createPB(seu.obj = seu.obj.sub, groupBy = "d90", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
         clusters = NULL, outDir = paste0("../output/", outName, "/pseudoBulk/")
)

#add dog metadata to account for paired nature of the data
df <- read.csv(paste0("../output/", outName, "/pseudoBulk/d90_deg_metaData.csv"), row.names = 1)
df$dog <- unlist(lapply(df$sampleID, function(x){strsplit(x, "_")[[1]][2]}))
write.csv(df, paste0("../output/", outName, "/pseudoBulk/d90_deg_metaData.csv"))

pseudoDEG(metaPWD = paste0("../output/", outName, "/pseudoBulk/d90_deg_metaData.csv"),
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("../output/", outName, "/pseudoBulk/"), 
          outName = outName, 
          paired = T, pairBy = "dog", test.use = "Wald", strict_lfc = F,
          idents.1_NAME = contrast[1], idents.2_NAME = contrast[2],
          inDir = paste0("../output/", outName, "/pseudoBulk/"), title = "d90 vs d0", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)

top.degs <- read.csv(paste0("../output/", outName, "/pseudoBulk/d90/", outName, "_cluster_d90_all_genes.csv")) %>% filter(log2FoldChange < 0) %>% arrange(padj) %>% pull(gene)

### Plot top DEGS
features <- top.degs[1:8]
set.seed(12)
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(seu.obj, downsample = min(table(seu.obj$cellSource)))
p <- FeaturePlot(seu.obj.sub,features = features, 
                 reduction = reduction, pt.size = 0.01, 
                 split.by = "cellSource", order = T, 
                 by.col = F) + labs(x = "UMAP1", y = "UMAP2") & theme(axis.text = element_blank(),
                                                                      axis.title.y.right = element_text(size = 16),
                                                                      axis.ticks = element_blank(),
                                                                      axis.title = element_blank(),
                                                                      axis.line = element_blank(),
                                                                      plot.title = element_text(size=16),
                                                                      title = element_blank(),
                                                                      plot.margin = unit(c(1, 0, 0, 0), "pt")
                                                                      ) & scale_color_gradient(breaks = pretty_breaks(n = 3), limits = c(NA, NA), low = "lightgrey", high = "darkblue")
ggsave(paste0("../output/", outName, "/", outName, "_split_d90_degs.png"), width = 16, height = 4)



Idents(seu.obj) <- "name"
seu.obj.sub <- subset(seu.obj, idents = c("pre_3", "pre_4", "post_3", "post_4"))
table(seu.obj.sub$name)

seu.obj.sub$d14 <- "d14"
seu.obj.sub$d14 <- as.factor(seu.obj$d14)
### Complete pseudobulk DGE by all cells
createPB(seu.obj = seu.obj.sub, groupBy = "d14", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
         clusters = NULL, outDir = paste0("../output/", outName, "/pseudoBulk/")
)

#add dog metadata to account for paired nature of the data
df <- read.csv(paste0("../output/", outName, "/pseudoBulk/d14_deg_metaData.csv"), row.names = 1)
df$dog <- unlist(lapply(df$sampleID, function(x){strsplit(x, "_")[[1]][2]}))
write.csv(df, paste0("../output/", outName, "/pseudoBulk/d14_deg_metaData.csv"))

pseudoDEG(metaPWD = paste0("../output/", outName, "/pseudoBulk/d14_deg_metaData.csv"),
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("../output/", outName, "/pseudoBulk/"), 
          outName = outName, 
          paired = T, pairBy = "dog", test.use = "Wald", strict_lfc = F,
          idents.1_NAME = contrast[1], idents.2_NAME = contrast[2],
          inDir = paste0("../output/", outName, "/pseudoBulk/"), title = "d14 vs d0", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)


top.degs <- read.csv(paste0("../output/", outName, "/pseudoBulk/d14/", outName, "_cluster_d14_all_genes.csv")) %>% filter(log2FoldChange > 0) %>% arrange(padj) %>% pull(gene)

### Plot top DEGS
features <- top.degs[c(2:6,8,9)]
set.seed(12)
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(seu.obj, downsample = min(table(seu.obj$cellSource)))
p <- FeaturePlot(seu.obj.sub,features = features, 
                 reduction = reduction, pt.size = 0.01, 
                 split.by = "cellSource", order = T, 
                 by.col = F) + labs(x = "UMAP1", y = "UMAP2") & theme(axis.text = element_blank(),
                                                                      axis.title.y.right = element_text(size = 16),
                                                                      axis.ticks = element_blank(),
                                                                      axis.title = element_blank(),
                                                                      axis.line = element_blank(),
                                                                      plot.title = element_text(size=16),
                                                                      title = element_blank(),
                                                                      plot.margin = unit(c(1, 0, 0, 0), "pt")
                                                                      ) & scale_color_gradient(breaks = pretty_breaks(n = 3), limits = c(NA, NA), low = "lightgrey", high = "darkblue")
ggsave(paste0("../output/", outName, "/", outName, "_split_d14_degs.png"), width = 16, height = 4)


### DGE intersection analysis
top.degs.d90 <- read.csv(paste0("../output/", outName, "/pseudoBulk/d90/", outName, "_cluster_d90_all_genes.csv")) %>% filter(log2FoldChange > 0) %>% arrange(padj) %>% pull(gene)
top.degs.d14 <- read.csv(paste0("../output/", outName, "/pseudoBulk/d14/", outName, "_cluster_d14_all_genes.csv")) %>% filter(log2FoldChange > 0) %>% arrange(padj) %>% pull(gene)
intersect(top.degs.d90, top.degs.d14)
# [1] "GBP1"

top.degs.d90 <- read.csv(paste0("../output/", outName, "/pseudoBulk/d90/", outName, "_cluster_d90_all_genes.csv")) %>% filter(log2FoldChange < 0) %>% arrange(padj) %>% pull(gene)
top.degs.d14 <- read.csv(paste0("../output/", outName, "/pseudoBulk/d14/", outName, "_cluster_d14_all_genes.csv")) %>% filter(log2FoldChange < 0) %>% arrange(padj) %>% pull(gene)
intersect(top.degs.d90, top.degs.d14)
# character(0)

top.degs.all <- read.csv(paste0("../output/", outName, "/pseudoBulk/allCells/", outName, "_cluster_allCells_all_genes.csv")) %>% filter(log2FoldChange > 0) %>% arrange(padj) %>% pull(gene)
top.degs.d90 <- read.csv(paste0("../output/", outName, "/pseudoBulk/d90/", outName, "_cluster_d90_all_genes.csv")) %>% filter(log2FoldChange > 0) %>% arrange(padj) %>% pull(gene)
intersect(top.degs.d90, top.degs.all)
# [1] "CCR1" "GBP1"


top.degs.all <- read.csv(paste0("../output/", outName, "/pseudoBulk/allCells/", outName, "_cluster_allCells_all_genes.csv")) %>% filter(log2FoldChange > 0) %>% arrange(padj) %>% pull(gene)
top.degs.d14 <- read.csv(paste0("../output/", outName, "/pseudoBulk/d14/", outName, "_cluster_d14_all_genes.csv")) %>% filter(log2FoldChange > 0) %>% arrange(padj) %>% pull(gene)
intersect(top.degs.d90, top.degs.d14)
# [1] "GBP1"

################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end analysis   ######## <<<<<<<<<<<<<<
################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
