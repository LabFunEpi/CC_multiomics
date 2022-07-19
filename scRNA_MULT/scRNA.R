
# Get Cell Ranger Summary
setwd("/cellranger_output/scRNA_cr_6.0.1")
samplenames <- c("PM001", "M4666", "C8", "M3399", "C3", "C1", "M5167", "C4", "072", "984", "B36", "B23", "B39", "B33", "454", "B31", "789", "B34", "B38", "570", "120", "B32", "B30", "890", "A6", "A25", "A42", "A16", "108", "A5", "A12", "A40", "A30", "A26")

cr_summary <- lapply(paste0(samplenames, "/outs/metrics_summary.csv"), function(x){read.csv(x) %>% t() %>% t() %>% data.frame()})
names(cr_summary) <- samplenames
cr_summary <- bind_rows(cr_summary, .id = "sample")

summ1 <- cr_summary %>% select(sample, Estimated.Number.of.Cells, Number.of.Reads, Fraction.Reads.in.Cells, Mean.Reads.per.Cell, Median.Genes.per.Cell, Total.Genes.Detected, Median.UMI.Counts.per.Cell, Sequencing.Saturation)
write.table(summ1, file = "summary.table.1", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

summ2 <- cr_summary %>% select(sample, Reads.Mapped.Confidently.to.Genome, Reads.Mapped.Confidently.to.Intergenic.Regions, Reads.Mapped.Confidently.to.Intronic.Regions, Reads.Mapped.Confidently.to.Exonic.Regions, Reads.Mapped.Confidently.to.Transcriptome, Reads.Mapped.Antisense.to.Gene)
write.table(summ2, file = "summary.table.2", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#################################################################################################

# Download and extract Haniffa data
# https://www.covid19cellatlas.org/index.patient.html
# cd /HaniffaData/
# wget https://covid19.cog.sanger.ac.uk/submissions/release1/haniffa21.processed.h5ad
# Raw
# https://drive.google.com/drive/folders/1KHcGqL_UbB6UyG7i-wpD8_0_ZhBldQbK 

library(anndata)
ad_raw <- read_h5ad("/HaniffaData/covid_raw_data_rna.h5ad")
ad_pro <- read_h5ad("/HaniffaData/haniffa21.processed.h5ad")

haniffa_raw <- CreateSeuratObject(counts = t(ad_raw$X[rownames(ad_pro), ]), assay = 'RNA', project = 'haniffa')
haniffa_raw@meta.data <- ad_pro$obs
SaveH5Seurat(haniffa_raw, "/HaniffaData/haniffa21.raw.h5seurat", overwrite = FALSE)

haniffa_raw <- LoadH5Seurat("/HaniffaData/haniffa21.raw.h5seurat")
haniffa_healthy <- subset(x = haniffa_raw, subset = Status == "Healthy")

haniffa_healthy$sample <- haniffa_healthy$patient_id
haniffa_healthy$scrublet.prediction <- "False"
haniffa_healthy$percent.mt <- haniffa_healthy$pct_counts_mt
ha_objs <- SplitObject(haniffa_healthy, split.by = "sample")
ha_objs <- lapply(ha_objs, function(obj){obj %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100) %>% FindNeighbors() %>% FindClusters()})
saveRDS(ha_objs, "/HaniffaData/haniffa21_objs.rds")

###################################################################################################

haniffa_raw <- LoadH5Seurat("/HaniffaData/haniffa21.raw.h5seurat")
haniffa_all <- haniffa_raw
haniffa_all$status <- haniffa_all$Status_on_day_collection_summary
haniffa_all$sample <- haniffa_all$patient_id
nm_map <- read.table("/CC_multiomics/scRNA_MULT/nm_map.txt", sep = "\t") %>% set_colnames(c("nm", "new"))
haniffa_all$predicted.celltype.nm <- (data.frame(nm = haniffa_all$full_clustering) %>% left_join(nm_map))$new
haniffa_all@meta.data <- haniffa_all@meta.data %>% select(Sex, predicted.celltype.nm, Age_interval, sample, status, full_clustering)
haniffa_all <- DietSeurat(haniffa_all, assays = "RNA", features = intersect(rownames(haniffa), rownames(GEX)))
haniffa_all <- haniffa_all %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000)
saveRDS(haniffa_all, "/HaniffaData/HaniffaAllDiet.rds")

haniffa <- readRDS("/HaniffaData/HaniffaAllDiet.rds")
haniffa <- subset(haniffa, subset = status == "Healthy")
library(celldex)
library(SingleR)
monaco <- celldex::MonacoImmuneData()
haniffa.se <- as.SingleCellExperiment(haniffa)
pred.haniffa.monaco <- SingleR(test = haniffa.se, ref = monaco, labels = monaco$label.fine, BPPARAM=MulticoreParam(8))
fine_map <- data.frame(colData(monaco)) %>% select(-label.ont) %>% distinct() %>% 
    mutate(new.fine = case_when(
        label.fine == "MAIT cells" ~ "MAIT", 
        (label.fine == "Vd2 gd T cells" | label.fine == "Non-Vd2 gd T cells") ~ "gdT",
        label.fine == "T regulatory cells" ~ "Treg", 
        label.fine == "Classical monocytes" ~ "CD14 Mono",
        label.fine == "Non classical monocytes" ~ "CD16 Mono",
        label.fine == "Intermediate monocytes" ~ "Int Mono",
        label.fine == "Plasmacytoid dendritic cells" ~ "pDC",
        label.fine == "Myeloid dendritic cells" ~ "cDC",
        label.fine == "Plasmablasts" ~ "Plasmablast",
        label.fine == "Progenitor cells" ~ "HSPC",
        label.fine == "Natural killer cells" ~ "NK",
        TRUE ~ label.main)) %>%
    mutate(new.fine = factor(case_when(
        new.fine == "CD8+ T cells" ~ "CD8 T",
        new.fine == "CD4+ T cells" ~ "CD4 T",
        new.fine == "B cells" ~ "B",
        TRUE ~ new.fine), levels = names(colors.mn))) %>%
    rename(labels = label.fine)
pred.haniffa.monaco <- data.frame(pred.haniffa.monaco) %>% left_join(fine_map)
saveRDS(pred.haniffa.monaco, "/HaniffaData/haniffa.healthy.SingleR.rds")

###################################################################################################

# # Process all data
setwd("/FINAL/scRNA")
cr_out_dir <- "/cellranger_output/scRNA_cr_6.0.1/"
samplenames <- c("PM001", "M4666", "C8", "M3399", "C3", "C1", "M5167", "C4", "072", "984", "B36", "B23", "B39", "B33", "454", "B31", "789", "B34", "B38", "570", "120", "B32", "B30", "890", "A6", "A25", "A42", "A16", "108", "A5", "A12", "A40", "A30", "A26")

print(paste0(Sys.time(), ":: Creating scRNA seurat objects ... "))

# Cell Ranger to Seurat
make_seurat_obj <- function(samplepath){
    counts <- Read10X_h5(paste0(samplepath, "/outs/filtered_feature_bc_matrix.h5"))
    obj <- CreateSeuratObject(counts = counts, assay = 'RNA', project = '10x_RNA')
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    obj <- obj %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100) %>% FindNeighbors() %>% FindClusters()
    saveRDS(obj, paste0(samplepath, "/outs/rna.rds"))
    return(obj)
}
objs <- lapply(paste0(cr_out_dir, samplenames), make_seurat_obj)

# Scrublet
do_scrub <- function(samplepath){
    # system(paste0("gunzip ", samplepath, "/outs/filtered_feature_bc_matrix/features.tsv.gz"))
    # system(paste0("gunzip ", samplepath, "/outs/filtered_feature_bc_matrix/matrix.mtx.gz"))
    system(paste0("python /CC_multiomics/scRNA_MULT/scrub.py ", samplepath, "/outs/filtered_feature_bc_matrix"))
}
invisible(lapply(paste0(cr_out_dir, samplenames), do_scrub))
add_scrub_to_obj <- function(samplepath){
    scrublet <- read.table(paste0(samplepath, "/outs/filtered_feature_bc_matrix/scrublet.tsv"), sep = "\t", col.names = c('scrublet.observed', 'scrublet.simulated', 'scrublet.prediction'))
    obj <- readRDS(paste0(samplepath, "/outs/rna.rds"))
    rownames(scrublet) <- colnames(obj)
    obj <- AddMetaData(obj, metadata = scrublet)
    saveRDS(obj, paste0(samplepath, "/outs/rna.rds"))
    return(obj)
}
objs <- lapply(paste0(cr_out_dir, samplenames), add_scrub_to_obj)

objs <- lapply(paste0(cr_out_dir, samplenames, "/outs/rna.rds"), readRDS)
names(objs) <- samplenames

# Add sample names
add.samplenames <- function(obj, samplename){
    obj$sample <- samplename
    return(obj)
}
mayo_objs <- mapply(add.samplenames, objs, names(objs))
mayo_objs <- lapply(mayo_objs, function(x){RenameCells(x, add.cell.id = x$sample[[1]])})

saveRDS(mayo_objs, "mayo_objs.rds")

############################## End of preprocess ######################################

meta <- read.table("/CC_multiomics/scRNA_MULT/scRNA.metadata.tsv", sep = "\t") %>% set_colnames(c("status", "sample", "age", "WHO", "CRS", "Sex"))
meta <- meta %>% mutate(sample = factor(sample, levels = (meta %>% arrange(age))$sample)) %>% arrange(age) %>% na_if(-1) %>% replace_na(list(age = "  ", WHO = "  ", CRS = "  "))

print(paste0(Sys.time(), ":: Integrating (Seurat v4 batch-correction) ... "))

objs <- mayo_objs

# Combine, integrate and process. 
objs <- lapply(objs, FUN = function(x){DefaultAssay(x) <- "RNA"; return(x)})
features <- SelectIntegrationFeatures(object.list = objs)
objs <- lapply(X = objs, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
# Choosing M3399, C4, B23, B30, A6, A30 as reference
ref_indices <- sort(match(c("M3399", "C4", "B23", "B30", "A6", "A30"), names(objs)))
anchors <- FindIntegrationAnchors(object.list = objs, reference = ref_indices, reduction = "rpca", anchor.features = features, dims = 1:50)
saveRDS(features, "features.rds")
saveRDS(anchors, "anchors.rds")
combined <- IntegrateData(anchorset = anchors, features.to.integrate = c(features, cc.genes$s.genes, cc.genes$g2m.genes), dims = 1:50)
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 50, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:50)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:50)
combined <- FindClusters(combined, resolution = 0.5)
combined$status <- (data.frame(sample = combined$sample) %>% left_join(meta))$status
combined$age <- (data.frame(sample = combined$sample) %>% left_join(meta))$age
combined$WHO <- (data.frame(sample = combined$sample) %>% left_join(meta))$WHO
combined$CRS <- (data.frame(sample = combined$sample) %>% left_join(meta))$CRS
combined$Sex <- (data.frame(sample = combined$sample) %>% left_join(meta))$Sex
saveRDS(combined, "combined.rds")

print(paste0(Sys.time(), ":: Processing integrated QC filtered object ... "))

# combined <- readRDS("combined.rds")
filtered <- subset(combined, subset = percent.mt < 50 & scrublet.prediction == "False" & nFeature_RNA > 200)
filtered <- subset(filtered, features = rownames(filtered@assays$RNA@counts)[rowSums(filtered@assays$RNA@counts > 0) >= 3])
filtered <- DietSeurat(filtered, scale.data = FALSE, assays = c("integrated", "RNA"))

DefaultAssay(filtered) <- "integrated"
filtered <- ScaleData(filtered, verbose = FALSE)
filtered <- RunPCA(filtered, npcs = 50, verbose = FALSE)
filtered <- RunUMAP(filtered, reduction = "pca", dims = 1:50)
filtered <- FindNeighbors(filtered, reduction = "pca", dims = 1:50)
filtered <- FindClusters(filtered, resolution = 0.5)
filtered$status <- (data.frame(sample = filtered$sample) %>% left_join(meta))$status
filtered$age <- (data.frame(sample = filtered$sample) %>% left_join(meta))$age
filtered$WHO <- (data.frame(sample = filtered$sample) %>% left_join(meta))$WHO
filtered$CRS <- (data.frame(sample = filtered$sample) %>% left_join(meta))$CRS
filtered$Sex <- (data.frame(sample = filtered$sample) %>% left_join(meta))$Sex

saveRDS(filtered, "filtered.rds")

#############################################################################

print(paste0(Sys.time(), ":: Cell-type identification using SingleR ... "))
# GEX <- readRDS("filtered.rds")

library(celldex)
monaco <- celldex::MonacoImmuneData()
colors.monaco <- c(`CD4 T` = "#81a66e", `Treg` = "#32584e", `CD8 T` = "#f5df87", `gdT` = "#b4c29f", `MAIT` = "#fdbd92", `NK` = "#546046", `B` = "#7c6289", `Plasmablast` = "#5991ba", `CD14 Mono` = "#d18c5a", `CD16 Mono` = "#f18e8c", `Int Mono` = "#8e402b", `cDC` = "#65a84c", `pDC` = "#b5d75b", `HSPC` = "#ffcbe7", `Neutrophils` = "#a2bcc2", `Basophils` = "#66a5a4")
fine_map <- data.frame(colData(monaco)) %>% select(-label.ont) %>% distinct() %>% 
    mutate(new.fine = case_when(
        label.fine == "MAIT cells" ~ "MAIT", 
        (label.fine == "Vd2 gd T cells" | label.fine == "Non-Vd2 gd T cells") ~ "gdT",
        label.fine == "T regulatory cells" ~ "Treg", 
        label.fine == "Classical monocytes" ~ "CD14 Mono",
        label.fine == "Non classical monocytes" ~ "CD16 Mono",
        label.fine == "Intermediate monocytes" ~ "Int Mono",
        label.fine == "Plasmacytoid dendritic cells" ~ "pDC",
        label.fine == "Myeloid dendritic cells" ~ "cDC",
        label.fine == "Plasmablasts" ~ "Plasmablast",
        label.fine == "Progenitor cells" ~ "HSPC",
        label.fine == "Natural killer cells" ~ "NK",
        TRUE ~ label.main)) %>%
    mutate(new.fine = factor(case_when(
        new.fine == "CD8+ T cells" ~ "CD8 T",
        new.fine == "CD4+ T cells" ~ "CD4 T",
        new.fine == "B cells" ~ "B",
        TRUE ~ new.fine), levels = names(colors.monaco))) %>%
    rename(labels = label.fine)

diet_GEX <- DietSeurat(GEX, assays = "RNA")
GEX.se <- as.SingleCellExperiment(diet_GEX)
library(SingleR)
pred.GEX.monaco <- SingleR(test = GEX.se, ref = monaco, labels = monaco$label.fine, BPPARAM=MulticoreParam(8))
pred.GEX.monaco <- data.frame(pred.GEX.monaco) %>% left_join(fine_map)
saveRDS(pred.GEX.monaco, "SingleR.monaco.fine.rds")


###############################################################################
print(paste0(Sys.time(), ":: Create MAYOPLUSHANIFFA ... "))

haniffa <- readRDS("/HaniffaData/HaniffaAllDiet.rds")
haniffa <- subset(haniffa, subset = status == "Healthy")
pred.haniffa.monaco <- readRDS("/HaniffaData/haniffa.healthy.SingleR.rds")
haniffa$monaco <- pred.haniffa.monaco$new.fine
haniffa$group <- case_when(haniffa$Age_interval %in% c("(20, 29]", "(30, 39]", "(40, 49]") ~ "CTRLY", TRUE ~ "CTRLO")
haniffa$mutation <- case_when(haniffa$Age_interval %in% c("(20, 29]", "(30, 39]", "(40, 49]") ~ "CTRLY", TRUE ~ "CTRLO")

GEX <- readRDS("filtered.rds")
mutmap <- data.frame(
    sample = c("A5", "A6", "A12", "A16", "A25", "A26", "A30", "A40", "A42", "108", "890", "M3399", "M4666", "M5167", "PM001", "C2", "C3", "C4", "C6", "C8", "L010-1", "L010-2", "L027", "C7"), 
    group = c("TET2C", "CHIPC", "TET2C", "CHIPC", "TET2C", "CHIPC", "TET2C", "CHIPC", "TET2C", "TET2C", "CHIPC", "TET2", "TET2", "CHIP", "TET2", "CHIP", "TET2", "TET2", "CHIP", "TET2", "TET2", "TET2", "TET2", "CHIP"),
    mutation = c("2xTET2-C", "DNMT3A-C", "TET2-C", "ASXL1-C", "TET2-C", "DNMT3A-C", "TET2-C", "2xDNMT3A-C", "TET2-C", "2xTET2-C", "EZH2-C", "2xTET2", "TET2", "ASXL1", "TET2", "DNMT3A", "TET2", "TET2/DNMT3A", "SF3B1", "2xTET2", "2xTET2", "2xTET2", "2xTET2", "ASXL1")
)
GEX$group <- (data.frame(sample = GEX$sample) %>% left_join(mutmap) %>% replace_na(list(group = "COVID")))$group
GEX$mutation <- (data.frame(sample = GEX$sample) %>% left_join(mutmap) %>% replace_na(list(mutation = "COVID")))$mutation
pred.GEX.monaco <- readRDS("SingleR.monaco.fine.rds")
GEX$monaco <- pred.GEX.monaco$new.fine

haniffa@meta.data <- haniffa@meta.data %>% select(Sex, nFeature_RNA, nCount_RNA, monaco, Age_interval, sample, group, mutation)
GEX@meta.data <- GEX@meta.data %>% select(Sex, nFeature_RNA, nCount_RNA, monaco, age, sample, group, mutation)
haniffa <- DietSeurat(haniffa, assays = "RNA", features = intersect(rownames(haniffa), rownames(GEX)))
GEX <- DietSeurat(GEX, assays = "RNA", features = intersect(rownames(haniffa), rownames(GEX)))

MAYOPLUSHANIFFA <- merge(GEX, y = haniffa)
MAYOPLUSHANIFFA <- MAYOPLUSHANIFFA %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000)

saveRDS(MAYOPLUSHANIFFA, "MAYOPLUSHANIFFA.rds")

###############################################

# Next, run the computation time intensive DGEA code (separate scripts - degs_auto_1.R, degs_auto_2.R)

###############################################
















