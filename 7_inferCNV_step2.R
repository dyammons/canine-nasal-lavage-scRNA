#!/usr/bin/Rscript

# ANALYSIS NOTE: this script extracts the inferCNV output to complete custom 
# plotting of the data.

outName <- "inferCNV_Ind"
reduction <- "umap.integrated.harmony"

#load in the all cells dataset and add metadata
seu.obj <- readRDS("../output/s3/nasal_lavage_n8_log_canFam_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_cfam.csv", 
                    groupBy = "clusterID_integrated.harmony", metaAdd = "majorID")

#subset on all but neuts to be able to load in the data using the inferCNV function
seu.obj <- subset(seu.obj, invert = T, subset = majorID == "neut")

meta.list <- lapply(dirNames, function(x){
    seu.obj <- subset(seu.obj, subset = orig.ident == x)
    seu.obj$rowname <- rownames(seu.obj@meta.data)
    seu.obj <- add_to_seurat(seurat_obj = seu.obj, infercnv_output_path = paste0("../output/inferCNV_Ind/", x))
    return(seu.obj@meta.data)
})

meta.df <- do.call(bind_rows, meta.list)

seu.obj <- readRDS("../output/s3/nasal_lavage_n8_log_canFam_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_cfam.csv", 
                    groupBy = "clusterID_integrated.harmony", metaAdd = "majorID")

seu.obj@meta.data <- seu.obj@meta.data %>%
    rownames_to_column() %>% 
    left_join(meta.df[ , ! colnames(meta.df) %in% colnames(seu.obj@meta.data)], by = "rowname") %>% 
    column_to_rownames()


#plot the results
pi <- DimPlot(seu.obj, 
              reduction = reduction, 
              group.by = "infercnv_subcluster",
              pt.size = 0.1,
              split.by = "orig.ident",
              ncol = 3,
              label = F,
              label.box = F,
              shuffle = TRUE
) 
p <- formatUMAP(pi)
ggsave(paste0("../output/", outName, "/", outName, "_test_quant.png"), width = 10.5, height = 7)



### Obsolete appraoch b/c HMM extraction is a better, more statistically sound approach

#load in the data
dirNames <- list.dirs(path = "../output/inferCNV_Ind/", full.names = F, recursive = F)
files <- paste0("../output/inferCNV_Ind/", dirNames, "/infercnv.observations.txt")
df.list.obs <- lapply(files, read.table, sep = " ")
files <- paste0( "../output/inferCNV_Ind/", dirNames, "/infercnv.references.txt")
df.list.all <- c(df.list.obs, lapply(files, read.table, sep = " "))
cnvStatus <- unlist(lapply(df.list.all, function(df){df <- colSums(abs(1-df))}))
names(cnvStatus) <- gsub("[.]", "-", names(cnvStatus))
df <- as.data.frame(cnvStatus)

#clean data -- split into quartiles
df <- df %>% 
    mutate(
        qunatile = case_when(
            cnvStatus < quantile(cnvStatus)[2] ~ "q1 (low CNV)",
            cnvStatus > quantile(cnvStatus)[2] & cnvStatus < quantile(cnvStatus)[3] ~ "q2",
            cnvStatus > quantile(cnvStatus)[3] & cnvStatus < quantile(cnvStatus)[4] ~ "q3",
            cnvStatus > quantile(cnvStatus)[4] ~ "q4 (high CNV)"),
        samID = substr(rownames(df),1,nchar(rownames(df))-19)
    ) %>% 
    group_by(samID) %>% 
    mutate(
        qunatile_grpd = case_when(
            cnvStatus < quantile(cnvStatus)[2] ~ "q1 (low CNV)",
            cnvStatus > quantile(cnvStatus)[2] & cnvStatus < quantile(cnvStatus)[3] ~ "q2",
            cnvStatus > quantile(cnvStatus)[3] & cnvStatus < quantile(cnvStatus)[4] ~ "q3",
            cnvStatus > quantile(cnvStatus)[4] ~ "q4 (high CNV)"),
        cnvCall = ifelse(qunatile_grpd == "q4 (high CNV)", "aneuploid", "diploid")
    )

#prep the metadata
new.metaData <- as.data.frame(df[ , c("cnvStatus", "cnvCall", "qunatile_grpd")])
new.metaData$rowname <- names(cnvStatus)
seu.obj@meta.data <- seu.obj@meta.data %>% 
    rownames_to_column() %>% 
    left_join(new.metaData, by = "rowname") %>% 
    column_to_rownames()

#plot the data
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 1, features = "cnvStatus", showAxis = F, smallAxis=F, 
                 reduction = reduction, color = "black", order = T, pt.size = 0.0000001, 
                 title.size = 14, titles = "CNV Status", noLegend = T)
ggsave(paste0("../output/", outName, "/", outName, "_inferCNV_status.png"), width = 6, height = 6)

cellz <- names(seu.obj$qunatile_grpd[!is.na(seu.obj$qunatile_grpd)])
seu.obj.sub <- subset(seu.obj, cells = cellz)
pi <- DimPlot(seu.obj.sub, 
              reduction = reduction, 
              group.by = "qunatile_grpd",
              pt.size = 0.1,
              split.by = "orig.ident",
              ncol = 3,
              label = F,
              label.box = F,
              shuffle = TRUE
) 
p <- formatUMAP(pi)
ggsave(paste0("../output/", outName, "/", outName, "_test_quant.png"), width = 10.5, height = 7)
