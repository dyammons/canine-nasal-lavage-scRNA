#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")
reduction <- "umap.integrated.harmony"
outName <- "misc"

### Plotting for Steve RE: IFN-activated neutrophils email chain
# Done on 03.06.2024
# Upadated 03.07.2024

#load in the full object
seu.obj <- readRDS(file = "../output/s3/240202_lav_neut_n8_2000_log_cfam_clean_S3.rds")


modulez <- list(
    "4T1-M" = c("IL1A", "TNF", "CXCL1", "SAA3", "CCL3", "CXCL3", "LY6E", "CCL4", "IL23A", "CXCL2"),
    #missing: CXCL1, SAA3, CXCL3, CXCL2
    "4T1-P" = c("WFDC17", "IFITM1", "RETNIG", "S100A11", "S100A9", "IFITM2", "LSP1", "JUN", "FOS", "SOC3"),
    #missing: WFDC17, IFITM1, RETNIG, S100A9, IFITM2, JUN, SOC3
    "IFN-hu" = c("IFIT3", "IFIT1", "IFIT2", "STAT1", "ISG15", "MX1", "IFI6", "RSAD2", "IFIT5", "LY6E")
    #missing: MX1
)

names(modulez) <- paste0(names(modulez),"_SIG")
seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

features <- names(modulez)
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

fig_supp <- majorDot(seu.obj = seu.obj, groupBy = "clusterID_integrated.harmony",
                     features = features
                    ) + theme(legend.position = "bottom",
                              axis.title.y = element_blank(),
                              plot.margin = margin(7, 7, 0, 50, "pt"),
                              legend.box="vertical", 
                              legend.margin=margin()
                             ) + scale_y_discrete(position = "right") + guides(size = guide_legend(nrow = 2, byrow = F, title = 'Percent\nenriched')) + guides(color = guide_colorbar(title = 'Scaled\nenrichment score'))  
ggsave(paste("../output/", outName, "/", outName, "_modScores.png", sep = ""), width = 4.5, height = 7)


features <- c("4T1-M_SIG","4T1-P_SIG","IFN-hu_SIG")
p <- FeaturePlot(seu.obj, features = features, pt.size = 0.000000001, order = F, by.col = F, ncol = 3, reduction = reduction
                ) + labs(x = "UMAP1", y = "UMAP2") & theme(axis.text = element_blank(),
                                                           axis.title.y.right = element_text(size = 16),
                                                           axis.ticks = element_blank(),
                                                           axis.title = element_blank(),
                                                           axis.line = element_blank(),
                                                           plot.title = element_text(size=16),
                                                           title = element_blank(),
                                                           plot.margin = unit(c(0, 0, 0, 0), "cm")
                                                          ) & scale_color_gradient(breaks = pretty_breaks(n = 3), 
                                                                                   limits = c(NA, NA), low = "lightgrey", 
                                                                                   high = "darkblue")  & scale_colour_viridis(option="magma", name='Expression')
ggsave(paste("../output/", outName, "/",outName, "_UMAP_enrichment.png", sep = ""), width = 9, height = 3)

features <- c("LY6E")
p <- p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 1, features = features, 
                    color = "black", order = F, pt.size = 0.0000001, title.size = 16, reduction = reduction)
ggsave(paste0("../output/", outName, "/", outName, "_LY6E.png"), width = 3, height = 3)

