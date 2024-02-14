#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")
library(circlize)
library(scProportionTest)

################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
###  BEGIN allCells analysis  ### <<<<<<<<<<<<<<
################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

### Set output params
# outName <- "nasal_lavage_n8_2000_log"
outName <- "nasal_lavage_n8_2000_log_cfam"
reduction <- "umap.integrated.harmony"
clusMain <- "clusterID_integrated.harmony"
contrast <- c("Post", "Pre")

### Load in preprocessed data & metadata
# seu.obj <- readRDS("../output/s3/nasal_lavage_n8.rds") #load in GSD (CanFam4 aligned) data
# seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID.csv", groupBy = "clusterID_integrated.harmony", metaAdd = "majorID")
seu.obj <- readRDS("../output/s3/nasal_lavage_n8_log_canFam_S3.rds") #load CanFam3.1 aligned data
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_cfam.csv", groupBy = "clusterID_integrated.harmony", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "cellSource")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "cellSource2")
seu.obj$cellSource2 <- factor(seu.obj$cellSource2, levels = c("d0","d14","d90"))
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj$clusterID_integrated.harmony <- factor(seu.obj$clusterID_integrated.harmony, levels = 0:max(as.numeric(names(table(seu.obj$clusterID_integrated.harmony)))))
seu.obj$clusterID  <- seu.obj$clusterID_integrated.harmony 
seu.obj <- convertTOclusID(seu.obj = seu.obj, metaSlot = "majorID", newMetaName = "major_clusID")

#set colors - run after determining majorIDs -- MANUAL
# colArray <- read.csv("./metaData/allCells_ID.csv")
colArray <- read.csv("./metaData/allCells_ID_cfam.csv")
# colArray <- colArray %>% arrange(majorID) %>% mutate(newCol = gg_color_hue(nrow(colArray)*3)[ c( rep(FALSE, 2), TRUE ) ] ) %>% arrange(clusterID_integrated.harmony)
# write.csv(colArray,"./metaData/allCells_ID.csv", row.names = F)
colArray <- colArray[ ,-c(3)]
colArray.sub <- na.omit(colArray[colArray$majCol == "yes",])

### Check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, reduction = reduction, nrow = 1, ncol = 3, features = features, 
                    color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste0("../output/", outName, "/", outName, "_QC_feats.png"), width = 9, height = 3)

### Generate viln plots using harmony clusters
vilnPlots(seu.obj = seu.obj, groupBy = clusMain, outName = outName,
          outDir = paste0("../output/viln/", outName, "/"), returnViln = T 
         )

### Generate viln plots using majorID
vilnPlots(seu.obj = seu.obj, groupBy = "major_clusID", outName = outName,
          outDir = paste0("../output/viln/", outName, "/"), returnViln = T 
         )

### Generate dot plots using vilnPlots resuts of majorID
pi <- autoDot(seu.integrated.obj = seu.obj, inFile = paste0("../output/viln/", outName, "/", outName, "_major_clusID_gene_list.csv"), groupBy = "major_clusID",
                     MIN_LOGFOLD_CHANGE = 0.5, MIN_PCT_CELLS_EXPR_GENE = 0.1,
                    filterTerm = "ENSCAFG"
                    ) + theme(legend.box="vertical") + scale_fill_manual(values = colArray.sub$newCol)
ggsave(paste0("../output/", outName, "/", outName, "_autodot_major_clusID.png"), width = 5, height = 10)

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
                cols = colArray$newCol,
                pt.size = 0.25,
                label = TRUE,
                label.box = TRUE
) + NoLegend()
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black", smallAxes = F) 
ggsave(paste0("../output/", outName, "/", outName, "_rawUMAP.png"), width = 7, height = 7)


### Plot initial cluster UMAP
pi <- DimPlot(seu.obj, 
                reduction = reduction, 
                group.by = "major_clusID",
                cols = colArray.sub$newCol,
                pt.size = 0.25,
                label = TRUE,
                label.box = TRUE
) + NoLegend()
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black", smallAxes = F) 
ggsave(paste0("../output/", outName, "/", outName, "_majorClusUMAP.png"), width = 7, height = 7)


### Fig supp: umap by major ID
pi <- DimPlot(seu.obj, 
              reduction = reduction, 
              group.by = "majorID",
              cols = colArray.sub$newCol,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
)
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_majorUMAP.png"), width = 7, height = 7)


### Key feature plots
features <- c("PTPRC","CD3E","CTSW", 
                "DRA","CSF3R","S100A12", 
                "CD68","FLT3","FCER1A", 
                "GPNMB","VEGFB","CD34",
                "COL1A2","MS4A1","TOP2A")
p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, title.size = 14, features = features, order = F, legJust = "top", reduction = reduction) 
ggsave(paste0("../output/", outName, "/", outName, "_featPlots.png"), width = 9, height = 15)


### Key dot plot features -- this is best with majorID loaded in
p <- majorDot(seu.obj = seu.obj, groupBy = "majorID",
              features = c("ANPEP", "DLA-DRA", "FLT3", "IGHM", "JCHAIN",
                           "MS4A1", "S100A12", "SERPINA1", "CD4", "IL7R", 
                           "CD52", "CCL5", "GZMB", "KLRB1", "CSTB", "IL5RA", 
                           "IL17RB", "GATA3", "TOP2A", "CENPF", "CD34", "CD109")
) + theme(axis.title = element_blank(),
          axis.text = element_text(size = 12))
ggsave(paste0("../output/", outName, "/", outName, "_majorDot.png"), width = 8, height = 6)


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


### Stacked bar graph by majorID -- preferred once variable is set
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "name", groupBy = "name", clusters = "majorID") + scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(labels = levels(seu.obj$name),
                    values = levels(seu.obj$colz)) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12))
ggsave(paste0("../output/", outName, "/", outName, "_stackedBar_major.png"), width =7, height = 5)


### Frequency plots to run stats - cluster
freqy <- freqPlots(seu.obj, method = 1, nrow = 3, 
                   comp = "cellSource", groupBy = clusMain, legTitle = "Cell source", refVal = "name",
                   namez = "name", 
                   colz = "colz"
)
ggsave(paste0("../output/", outName, "/",outName, "_freqPlots_cluster.png"), width = 12, height = 8)


### Frequency plots to run stats - major with 3 levels
freqy <- freqPlots(seu.obj, method = 1, nrow = 1, 
                   comp = "cellSource2", groupBy = "majorID", legTitle = "Cell source", refVal = "name",
                   namez = "name", 
                   colz = "colz"
)
ggsave(paste0("../output/", outName, "/",outName, "_freqPlots_major2.png"), width = 8, height = 4)


### Frequency plots to run stats - major
freqy <- freqPlots(seu.obj, method = 1, nrow = 1, 
                   comp = "cellSource", groupBy = "majorID", legTitle = "Cell source", refVal = "name",
                   namez = "name", 
                   colz = "colz"
)
ggsave(paste0("../output/", outName, "/",outName, "_freqPlots_major.png"), width = 8, height = 4)

#run paired analysis
df <- freqy$data %>% separate(name, into = c(NA, "dog"))

lapply(levels(df$majorID), function(x){
    df.sub <- filter(df, majorID == x)
    wilcox.test(freq ~ cellSource, data = df.sub, paired = TRUE)
})


### Evaluate cell frequency by cluster using monte carlo permutation
log2FD_threshold <- 0.58
Idents(seu.obj) <- "name"
set.seed(12)
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$name)))
prop_test <- sc_utils(seu.obj.sub)
prop_test <- permutation_test( prop_test, cluster_identity = "majorID", sample_1 = "Pre", sample_2 = "Post", sample_identity = "cellSource" )

p <- permutation_plot(prop_test)  + theme(axis.title.y = element_blank(),
                                          legend.position = "top") + guides(colour = guide_legend("", nrow = 2, 
                                                                                                  byrow = TRUE)) + coord_flip()
res.df <- prop_test@results$permutation
res.df <- res.df %>% mutate(Significance = as.factor(ifelse(obs_log2FD < -log2FD_threshold & FDR < 0.01,"Down",
                                                            ifelse(obs_log2FD > log2FD_threshold & FDR < 0.01,"Up","n.s.")))
                           ) %>% arrange(obs_log2FD)

res.df$clusters <- factor(res.df$clusters, levels = c(res.df$clusters))
p <- ggplot(res.df, aes(x = clusters, y = obs_log2FD)) + 
geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, 
                    color = Significance)) + theme_bw() + geom_hline(yintercept = log2FD_threshold, 
                                                                     lty = 2) + geom_hline(yintercept = -log2FD_threshold, 
                                                                                           lty = 2) + 
geom_hline(yintercept = 0) + scale_color_manual(values = c("blue", "red","grey")
                                               ) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 14),
                                                         axis.title.x = element_blank(),
                                                         legend.position = "top",
                                                         plot.margin = margin(c(3,3,0,24))
                                                        ) + ylab("abundance change (log2FC)")

ggsave(paste("../output/", outName, "/",outName, "_propTest.png", sep = ""), width = 3.5, height = 2, scale = 2 )


### Use lenitent Wilcoxon rank sum test
seu.obj$allCells <- "All cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
linDEG(seu.obj = seu.obj, threshold = 1, thresLine = F, groupBy = "allCells", comparison = "cellSource",contrast= c("Post","Pre"),
       outDir = paste0("../output/", outName, "/"), outName = outName, 
       pValCutoff = 0.01, logfc.threshold = 0.58, saveGeneList = T, addLabs = ""
      )


### Complete GSEA using the linDEG results
df <- read.csv("../output/nasal_lavage_n8_2000_log_cfam/pseudoBulk/allCells/nasal_lavage_n8_2000_log_cfam_cluster_allCells_all_genes.csv") %>% arrange(padj)
upGenes <- df %>% filter(log2FoldChange > 0) %>% pull(gene)
dwnGenes <- df %>% filter(log2FoldChange < 0) %>% pull(gene)
p <- plotGSEA(geneList = upGenes, geneListDwn = dwnGenes, category = "C5", subcategory = "GO:BP", 
              upCol = "blue", dwnCol = "red", size = 3)

minVal <- -15
maxVal <- 25
pi <- p + scale_x_continuous(limits = c(minVal, maxVal), name = "Signed log10(padj)") + 
    theme(axis.title=element_text(size = 16)) + 
    annotate("segment", x = -0.1, 
             y = 17, 
             xend = minVal, 
             yend = 17, 
             lineend = "round", linejoin = "bevel", linetype ="solid", colour = "blue",
             size = 1, arrow = arrow(length = unit(0.1, "inches"))
            ) + 
    annotate(geom = "text", x = (minVal-0.1*1.5)/2-0.1*1.5,
             y = 18,
             label = "Repressed",
             hjust = 0.5,
             vjust = 1.5,
             size = 5) +
    annotate("segment", x = 0.1,
             y = 17,
             xend = maxVal,
             yend = 17,
             lineend = "round", linejoin = "bevel", linetype ="solid", colour = "red",
             size = 1, arrow = arrow(length = unit(0.1, "inches"))
            ) + 
    annotate(geom = "text", x = (maxVal-0.1*1.5)/2+0.1*1.5, 
             y = 18,
             label = "Induced",
             hjust = 0.5,
             vjust = 1.5,
             size = 5)
ggsave(paste("../output/", outName, "/", outName, "_allCells_gsea.png", sep = ""), width = 10, height = 7)


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

top.degs <- read.csv("../output/nasal_lavage_n8_2000_log_cfam/pseudoBulk/allCells/nasal_lavage_n8_2000_log_cfam_cluster_allCells_all_genes.csv") %>% filter(log2FoldChange > 0) %>% arrange(padj) %>% pull(gene)

### Plot top DEGS
features <- top.degs[1:6]
features <- top.degs[7:12]
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
ggsave(paste0("../output/", outName, "/", outName, "_split_degs.png"), width = 12, height = 4)


### Evalute the DGE results for overlap with nanostring data
all.degs <- read.csv("../output/nasal_lavage_n8_2000_log_cfam/pseudoBulk/allCells/nasal_lavage_n8_2000_log_cfam_cluster_allCells_all_genes.csv")
nanoGenes <- read.csv("./metaData/nanoString_genes.csv", header = F)

all.degs$gene[all.degs$gene %in% nanoGenes$V1] %>% length()
#155

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
          inDir = paste0("../output/", outName, "/pseudoBulk/"), title = "All cells", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)

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
          inDir = paste0("../output/", outName, "/pseudoBulk/"), title = "All cells", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)



### Complete linDEG in each major group
linDEG(seu.obj = seu.obj, groupBy = "majorID", comparison = "cellSource", contrast= c("Post","Pre"),
       outDir = paste0("../output/", outName, "/linDEG/"), outName = outName, 
       pValCutoff = 0.01, logfc.threshold = 0.58, saveGeneList = T, addLabs = ""
)

### Complete pseudobulk DGE by each majorID cell type
createPB(seu.obj = seu.obj, groupBy = "majorID", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
         clusters = NULL, outDir = paste0("../output/", outName, "/pseudoBulk/")
)

#add metadata to properly pair the analysis
df <- read.csv(paste0("../output/", outName, "/pseudoBulk/majorID_deg_metaData.csv"), row.names = 1)
df$dog <- unlist(lapply(df$sampleID, function(x){strsplit(x, "_")[[1]][2]}))
write.csv(df, paste0("../output/", outName, "/pseudoBulk/majorID_deg_metaData.csv"))

pseudoDEG(metaPWD = paste0("../output/", outName, "/pseudoBulk/majorID_deg_metaData.csv"),
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("../output/", outName, "/pseudoBulk/"), 
          outName = outName, 
          idents.1_NAME = contrast[1], idents.2_NAME = contrast[2],
          inDir = paste0("../output/", outName, "/pseudoBulk/"), title = "All cells", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)



### Continue to investigate if d14 or d90 is driving DEGs -- this approach is not currently working
Idents(seu.obj) <- "name"
seu.obj.sub <- subset(seu.obj, idents = c("pre_1", "pre_2", "post_1", "post_2"))
table(seu.obj.sub$name)

### Complete pseudobulk DGE by all cells
createPB(seu.obj = seu.obj.sub, groupBy = "majorID", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
         clusters = NULL, outDir = paste0("../output/", outName, "/pseudoBulk/")
)

#add dog metadata to account for paired nature of the data
df <- read.csv(paste0("../output/", outName, "/pseudoBulk/majorID_deg_metaData.csv"), row.names = 1)
df$dog <- unlist(lapply(df$sampleID, function(x){strsplit(x, "_")[[1]][2]}))
write.csv(df, paste0("../output/", outName, "/pseudoBulk/majorID_deg_metaData.csv"))

pseudoDEG(metaPWD = paste0("../output/", outName, "/pseudoBulk/majorID_deg_metaData.csv"),
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("../output/", outName, "/pseudoBulk/"), 
          outName = paste0(outName, "_d90"), 
          paired = T, pairBy = "dog", test.use = "Wald", strict_lfc = F,
          idents.1_NAME = contrast[1], idents.2_NAME = contrast[2],
          inDir = paste0("../output/", outName, "/pseudoBulk/"), title = "All cells", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)

Idents(seu.obj) <- "name"
seu.obj.sub <- subset(seu.obj, idents = c("pre_3", "pre_4", "post_3", "post_4"))
table(seu.obj.sub$name)

seu.obj.sub$d14 <- "d14"
seu.obj.sub$d14 <- as.factor(seu.obj.sub$d14)
### Complete pseudobulk DGE by all cells
createPB(seu.obj = seu.obj.sub, groupBy = "majorID", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
         clusters = NULL, outDir = paste0("../output/", outName, "/pseudoBulk/")
)

#add dog metadata to account for paired nature of the data
df <- read.csv(paste0("../output/", outName, "/pseudoBulk/majorID_deg_metaData.csv"), row.names = 1)
df$dog <- unlist(lapply(df$sampleID, function(x){strsplit(x, "_")[[1]][2]}))
write.csv(df, paste0("../output/", outName, "/pseudoBulk/majorID_deg_metaData.csv"))

pseudoDEG(metaPWD = paste0("../output/", outName, "/pseudoBulk/majorID_deg_metaData.csv"),
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("../output/", outName, "/pseudoBulk/"), 
          outName = paste0(outName, "_d14"), 
          paired = T, pairBy = "dog", test.use = "Wald", strict_lfc = F,
          idents.1_NAME = contrast[1], idents.2_NAME = contrast[2],
          inDir = paste0("../output/", outName, "/pseudoBulk/"), title = "All cells", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)


### heatmap of dge results by major cell types
files <- lapply(levels(seu.obj$majorID), function(x){paste0("../output/", outName, "/pseudoBulk/", x, "/", outName, "_cluster_", x, "_all_genes.csv")})
files <- lapply(levels(seu.obj$majorID), function(x){paste0("../output/", outName, "/pseudoBulk/", x, "/", outName, "_d14_cluster_", x, "_all_genes.csv")})
files <- lapply(levels(seu.obj$majorID), function(x){paste0("../output/", outName, "/pseudoBulk/", x, "/", outName, "_d90_cluster_", x, "_all_genes.csv")})
files <- files[1:4] #run for d90 b/c too few B cells

df.list <- lapply(files, read.csv, header = T)

cnts_mat <- do.call(rbind, df.list)  %>% mutate(direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>% group_by(gs_base,direction) %>% summarize(nRow = n()) %>% pivot_wider(names_from = gs_base, values_from = nRow) %>% as.matrix() %>% t()
colnames(cnts_mat) <- cnts_mat[1,]
cnts_mat <- cnts_mat[-c(1),]
class(cnts_mat) <- "numeric"
cnts_mat[is.na(cnts_mat)] <- 0

#order by number of total # of DEGs
orderList <- rev(rownames(cnts_mat)[order(rowSums(cnts_mat))])
cnts_mat <- cnts_mat[match(orderList, rownames(cnts_mat)),]        

png(file = paste0("../output/", outName, "/", outName, "_deg_heat.png"), width=1500, height=2000, res=400)
par(mfcol=c(1,1))         
ht <- Heatmap(cnts_mat,#name = "mat", #col = col_fun,
              name = "# of DEGs",
              cluster_rows = F,
              row_title = "Cell type",
              col=colorRamp2(c(0,max(cnts_mat)), colors = c("white","red")),
              cluster_columns = F,
              column_title = "# of DEGs",
              show_column_names = TRUE,
              column_title_side = "top",
              column_names_rot = 0,
              column_names_centered = TRUE,
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                                          labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm")),
              cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(sprintf("%.0f", as.matrix(cnts_mat)[i, j]), x, y, gp = gpar(fontsize = 14, col = "black"))
              })
draw(ht, padding = unit(c(2, 12, 2, 5), "mm"),show_heatmap_legend = FALSE)
dev.off()

