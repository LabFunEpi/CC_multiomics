# Get Cell Ranger Summary
setwd("/cellranger_output/scMULT_v2/")
samplenames <- c("C2", "C5", "C8", "C3", "C1", "C4", "B36", "B23", "B39", "B31", "B34", "B30", "A6", "A25", "A16", "A5", "A12", "A40")

cr_summary <- lapply(paste0(samplenames, "/outs/summary.csv"), function(x){read.csv(x) %>% t() %>% t() %>% data.frame()})
names(cr_summary) <- samplenames
cr_summary <- bind_rows(cr_summary)

write.table(cr_summary, file = "summary.table", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

#####################################################################

setwd("/FINAL/scMULT")
datadir <- "/cellranger_output/scMULT_v2/"
samplenames <- c("C2", "C5", "C8", "C3", "C1", "C4", "B36", "B23", "B39", "B31", "B34", "B30", "A6", "A25", "A16", "A5", "A12", "A40")

make_mult_obj <- function(sampledir){
    h5 <- paste0(sampledir, "/outs/filtered_feature_bc_matrix.h5")
    frag.file <- paste0(sampledir, "/outs/atac_fragments.tsv.gz")
    metadata_csv <- paste0(sampledir, "/outs/per_barcode_metrics.csv")
    
    inputdata.10x <- Read10X_h5(h5)
    rna_counts <- inputdata.10x$`Gene Expression`
    atac_counts <- inputdata.10x$Peaks
    metadata <- read.csv(
		file = metadata_csv,
		header = TRUE,
		row.names = 1
	)
	sobj <- CreateSeuratObject(counts = rna_counts, meta.data = metadata)
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
    
    grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
    grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac_counts <- atac_counts[as.vector(grange.use), ]

    chrom_assay <- CreateChromatinAssay(
        counts = atac_counts,
        sep = c(":", "-"),
        genome = 'hg38',
        fragments = frag.file,
        min.cells = 10,
        annotation = annotations
    )
    sobj[["ATAC"]] <- chrom_assay
    DefaultAssay(sobj) <- "ATAC"
	sobj <- NucleosomeSignal(object = sobj)
	sobj <- TSSEnrichment(object = sobj, fast = FALSE)
	
    # RNA analysis
    DefaultAssay(sobj) <- "RNA"
    sobj <- sobj %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100) %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

    # ATAC analysis
    # We exclude the first dimension as this is typically correlated with sequencing depth
    DefaultAssay(sobj) <- "ATAC"
    sobj <- RunTFIDF(sobj)
    sobj <- FindTopFeatures(sobj, min.cutoff = 'q0')
    sobj <- RunSVD(sobj)
    sobj <- RunUMAP(sobj, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

    sobj <- FindMultiModalNeighbors(sobj, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
    sobj <- RunUMAP(sobj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    sobj <- FindClusters(sobj, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
    
    saveRDS(sobj, paste0(sampledir, "/outs/seurat_object.rds"))
    return(sobj)
}

print(paste0(Sys.time(), ":: Creating MULT object list ... "))

# Load annotations for hg38 genome
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

objs <- lapply(paste0(datadir, samplenames), make_mult_obj)
names(objs) <- samplenames

objs <- lapply(paste0(datadir, samplenames, "/outs/seurat_object.rds"), readRDS)
names(objs) <- samplenames

# Add sample names
add.samplenames <- function(obj, samplename){
    obj$sample <- samplename
    return(obj)
}
objs <- mapply(add.samplenames, objs, names(objs))
objs <- lapply(objs, function(x){RenameCells(x, add.cell.id = x$sample[[1]])})
saveRDS(objs, "mayo_objs.rds")


print(paste0(Sys.time(), ":: Integrating RNA part ... "))

# objs <- readRDS("mayo_objs.rds")
objs <- lapply(objs, function(x){DefaultAssay(x) <- "RNA"; return(x)})
features <- SelectIntegrationFeatures(object.list = objs)
objs <- lapply(X = objs, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = objs, anchor.features = features, dims = 1:50)
saveRDS(features, "rna_features.rds")
saveRDS(anchors, "rna_anchors.rds")
rna_integrated <- IntegrateData(anchorset = anchors, new.assay.name = "integrated_RNA", features.to.integrate = c(features, cc.genes$s.genes, cc.genes$g2m.genes), dims = 1:50)
DefaultAssay(rna_integrated) <- "integrated_RNA"
rna_integrated <- ScaleData(rna_integrated) %>% RunPCA(npcs = 50) %>% RunUMAP(reduction = "pca", dims = 1:50, reduction.name = "umap.irna", reduction.key = "irnaUMAP_")

p1 <- ggplot(data.frame(rna_integrated[["umap.irna"]][[]], label = factor(rna_integrated$sample)), aes(x=irnaUMAP_1, y=irnaUMAP_2, color=label)) + 
    ggrastr::rasterise(geom_point(size=0.2, stroke=0.1, shape=16), dpi = 400) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot()
pdf(file='umap-rna-integrated.pdf', width=8, height=8)
p1
dev.off()
saveRDS(rna_integrated, "rna_integrated.rds")



print(paste0(Sys.time(), ":: Creating ATAC object list ... "))

peakslist <- lapply(paste0(datadir, samplenames), FUN = function(x){
    peaks <- read.table(paste0(x, "/outs/atac_peaks.bed"), col.names = c("chr", "start", "end"));
    return(peaks)
})
allpeaks <- makeGRangesFromDataFrame(bind_rows(peakslist))
combined.peaks <- GenomicRanges::reduce(allpeaks)
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks <- combined.peaks[as.vector(seqnames(combined.peaks) %in% standardChromosomes(combined.peaks))]

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

atacobjlist <- lapply(paste0(datadir, samplenames), FUN = function(x){
    md <- read.table(
        file = paste0(x, "/outs/per_barcode_metrics.csv"),
        stringsAsFactors = FALSE,
        sep = ",",
        header = TRUE,
        row.names = 1
    )[-1, ]
    md <- md %>% filter(is_cell == 1)
    fragments <- CreateFragmentObject(paste0(x, "/outs/atac_fragments.tsv.gz"), cells = rownames(md));
    counts <- FeatureMatrix(fragments = fragments, features = combined.peaks, cells = rownames(md));
    assay <- CreateChromatinAssay(counts, fragments = fragments)
    obj <- CreateSeuratObject(assay, assay = "mATAC")
    return(obj)
})
names(atacobjlist) <- samplenames
add.samplenames <- function(obj, samplename){
    obj$sample <- samplename
    return(obj)
}
atacobjlist <- mapply(add.samplenames, atacobjlist, names(atacobjlist))
atacobjlist <- lapply(atacobjlist, function(x){
    x <- x %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 20) %>% RunSVD();
    return(x)
})
atacobjlist <- lapply(atacobjlist, function(x){RenameCells(x, add.cell.id = x$sample[[1]])})
saveRDS(atacobjlist, "atac_objs.rds")


print(paste0(Sys.time(), ":: Merging ATAC objects ... "))

merged_atac <- merge(
x = atacobjlist[[1]],
y = atacobjlist[2:length(atacobjlist)]
)
DefaultAssay(merged_atac) <- "mATAC"
merged_atac <- RunTFIDF(merged_atac)
merged_atac <- FindTopFeatures(merged_atac, min.cutoff = 20)
merged_atac <- RunSVD(merged_atac)
merged_atac <- RunUMAP(merged_atac, reduction = 'lsi', dims = 2:50, reduction.name = "umap.matac", reduction.key = "matacUMAP_")

p1 <- ggplot(data.frame(merged_atac[["umap.matac"]][[]], label = factor(merged_atac$sample)), aes(x=matacUMAP_1, y=matacUMAP_2, color=label)) + 
ggrastr::rasterise(geom_point(size=0.2, stroke=0.1, shape=16), dpi = 400) +
guides(color = guide_legend(override.aes = list(size=5))) +
theme_cowplot()
pdf(file='umap-atac-merged.pdf', width=8, height=8)
p1
dev.off()
saveRDS(merged_atac, "merged_atac.rds")

print(paste0(Sys.time(), ":: Integrating ATAC part ... "))

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = atacobjlist,
  anchor.features = rownames(atacobjlist[[1]]),
  reduction = "rlsi",
  dims = 2:50
)
saveRDS(integration.anchors, "atac_anchors.rds")

# integrate LSI embeddings
atac_integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = merged_atac[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:50,
  k.weight = 10
)

# create a new UMAP using the integrated embeddings
atac_integrated <- RunUMAP(atac_integrated, reduction = "integrated_lsi", dims = 2:50, reduction.name = "umap.iatac", reduction.key = "iatacUMAP_")
p1 <- ggplot(data.frame(atac_integrated[["umap.iatac"]][[]], label = factor(atac_integrated$sample)), aes(x=iatacUMAP_1, y=iatacUMAP_2, color=label)) + 
    ggrastr::rasterise(geom_point(size=0.2, stroke=0.1, shape=16), dpi = 400) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot()
pdf(file='umap-atac-integrated.pdf', width=8, height=8)
p1
dev.off()

saveRDS(atac_integrated, "atac_integrated.rds")

print(paste0(Sys.time(), ":: Putting integrated RNA and ATAC part together and running WNN analysis ... "))

# setwd("/FINAL/scMULT")
# rna_integrated <- readRDS("rna_integrated.rds")
# atac_integrated <- readRDS("atac_integrated.rds")

mult_integrated <- rna_integrated
mult_integrated[["mATAC"]] <- atac_integrated[["mATAC"]]
mult_integrated[["integrated_lsi"]] <- atac_integrated[["integrated_lsi"]]
mult_integrated[["umap.iatac"]] <- atac_integrated[["umap.iatac"]]
mult_integrated <- FindMultiModalNeighbors(mult_integrated, reduction.list = list("pca", "integrated_lsi"), dims.list = list(1:50, 2:50))
mult_integrated <- RunUMAP(mult_integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
mult_integrated <- FindClusters(mult_integrated, graph.name = "wsnn", algorithm = 3)
p1 <- ggplot(data.frame(mult_integrated[["wnn.umap"]][[]], label = factor(Idents(mult_integrated))), aes(x=wnnUMAP_1, y=wnnUMAP_2, color=label)) + 
    ggrastr::rasterise(geom_point(size=0.2, stroke=0.1, shape=16), dpi = 400) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot()
pdf(file='umap-wnn-integrated-clust.pdf', width=8, height=8)
p1
dev.off()

Annotation(mult_integrated[["mATAC"]]) <- annotations
mult_integrated <- NucleosomeSignal(object = mult_integrated, assay = "mATAC")
mult_integrated <- TSSEnrichment(object = mult_integrated, fast = FALSE, assay = "mATAC")

saveRDS(mult_integrated, "mult_integrated.rds")

# mult_integrated <- readRDS("mult_integrated.rds")

bettercolors <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff', '#000000')

plotdata <- data.frame(
    ident = "all", 
    nFeature_RNA = mult_integrated$nFeature_RNA, 
    nCount_RNA = mult_integrated$nCount_RNA,
    percent.mito = mult_integrated$percent.mt,
    nFeature_mATAC = mult_integrated$nFeature_mATAC, 
    nCount_mATAC = mult_integrated$nCount_mATAC,
    TSS.enrichment = mult_integrated$TSS.enrichment,
    nucleosome_signal = mult_integrated$nucleosome_signal,
    sample = mult_integrated$sample
)
plotdata1 <- data.frame(ident = "all", cutsites = rowSums(mult_integrated[['mATAC']]@counts), cells = rowSums(mult_integrated[['mATAC']]@counts > 0))

forlegend <- ggplot(plotdata, aes(x=ident, y=nFeature_RNA)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.05)) + 
        ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, aes(color = sample)), dpi = 400) +
        scale_y_continuous(trans = pseudolog10_trans) +
        scale_color_manual(values = bettercolors) +
        geom_hline(yintercept = 200) +
        ylab("No. of unique genes detected (RNA)") +
        guides(color = guide_legend(override.aes = list(size=5))) +
        theme_cowplot()
mylegend <- get_legend(forlegend)

p1 <- ggplot(plotdata, aes(x=ident, y=nFeature_RNA)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.05)) + 
        ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, show.legend = FALSE, aes(color = sample)), dpi = 400) +
        scale_color_manual(values = bettercolors) +
        scale_y_continuous(trans = pseudolog10_trans, breaks = log_breaks(n = 10, base = 10), labels = label_number(accuracy = 1, big.mark = "")) +
        geom_hline(yintercept = 200) +
        ylab("No. of unique genes detected (RNA)") +
        theme_cowplot()
p2 <- ggplot(plotdata, aes(x=ident, y=nCount_RNA)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.05)) + 
        ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, show.legend = FALSE, aes(color = sample)), dpi = 400) +
        scale_color_manual(values = bettercolors) +
        scale_y_continuous(trans = pseudolog10_trans, breaks = log_breaks(n = 10, base = 10), labels = label_number(accuracy = 1, big.mark = "")) +
        ylab("Total no. of molecules detected (RNA)") +
        theme_cowplot()
p3 <- ggplot(plotdata, aes(x=ident, y=nFeature_mATAC)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.05)) + 
        ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, show.legend = FALSE, aes(color = sample)), dpi = 400) +
        scale_color_manual(values = bettercolors) +
        scale_y_continuous(trans = pseudolog10_trans, breaks = log_breaks(n = 10, base = 10), labels = label_number(accuracy = 1, big.mark = "")) +
        geom_hline(yintercept = 200) +
        ylab("No. of unique peaks detected (ATAC)") +
        theme_cowplot()
p4 <- ggplot(plotdata, aes(x=ident, y=nCount_mATAC)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.05)) + 
        ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, show.legend = FALSE, aes(color = sample)), dpi = 400) +
        scale_color_manual(values = bettercolors) +
        scale_y_continuous(trans = pseudolog10_trans, breaks = log_breaks(n = 10, base = 10), labels = label_number(accuracy = 1, big.mark = "")) +
        ylab("Total no. of cutsites detected (ATAC)") +
        theme_cowplot()
p5 <- ggplot(plotdata, aes(x=ident, y=percent.mito)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.1)) + 
        ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, show.legend = FALSE, aes(color = sample)), dpi = 400) +
        scale_color_manual(values = bettercolors) +
        scale_y_continuous(breaks = seq(0, 100, 10)) +
        geom_hline(yintercept = 50) +
        ylab("Percent reads mapping to mito (RNA)") +
        theme_cowplot()
p6 <- ggplot(plotdata, aes(x=ident, y=TSS.enrichment)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.1)) + 
        ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, show.legend = FALSE, aes(color = sample)), dpi = 400) +
        scale_color_manual(values = bettercolors) +
        coord_cartesian(ylim = c(0, 20)) +
        scale_y_continuous(breaks = seq(0, 20, 1)) +
        geom_hline(yintercept = 1) +
        ylab("TSS enrichment score (ATAC)") +
        theme_cowplot()
p7 <- ggplot(plotdata, aes(x=ident, y=nucleosome_signal)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.1)) + 
        ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, show.legend = FALSE, aes(color = sample)), dpi = 400) +
        scale_color_manual(values = bettercolors) +
        # scale_y_continuous(breaks = seq(0, 100, 10)) +
        ylab("Nucleosome Signal (ATAC)") +
        theme_cowplot()
p8 <- ggplot(plotdata1, aes(x=ident, y=cutsites)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.1)) + 
        ggrastr::rasterise(geom_jitter(size = 0.1, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        scale_y_continuous(trans = pseudolog10_trans, breaks = log_breaks(n = 10, base = 10), labels = label_number(accuracy = 1, big.mark = "")) +
        ylab("Total no. of cutsites (ATAC)") +
        theme_cowplot()
p9 <- ggplot(plotdata1, aes(x=ident, y=cells)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.1)) + 
        ggrastr::rasterise(geom_jitter(size = 0.1, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        geom_hline(yintercept = 50) +
        scale_y_continuous(trans = pseudolog10_trans, breaks = log_breaks(n = 10, base = 10), labels = label_number(accuracy = 1, big.mark = "")) +
        ylab("Total no. of cells with at least one cutsite (ATAC)") +
        theme_cowplot()

layout <- "
ABCDI
EFGHJ
"
pdf(file='qc.pdf', width=16, height=12)
plot(p1 + p2 + p3 + p4 + p5 + as_ggplot(mylegend) + p6 + p7 + p8 + p9 + plot_layout(design = layout) & theme(axis.title.x=element_blank(), axis.text.x=element_blank()))
dev.off()


print(paste0(Sys.time(), ":: Applying QC filters and running WNN on filtered object ... "))
# mult_integrated <- readRDS("mult_integrated.rds")

filtered <- subset(mult_integrated, subset = percent.mt < 50 & nFeature_RNA > 200 & nFeature_mATAC > 200 & TSS.enrichment > 1)
filtered <- RunUMAP(filtered, reduction = "pca", dims = 1:50, reduction.name = "umap.irna", reduction.key = "irnaUMAP_")
filtered <- RunUMAP(filtered, reduction = "integrated_lsi", dims = 2:50, reduction.name = "umap.iatac", reduction.key = "iatacUMAP_")
filtered <- FindMultiModalNeighbors(filtered, reduction.list = list("pca", "integrated_lsi"), dims.list = list(1:50, 2:50))
filtered <- RunUMAP(filtered, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
filtered$sample_alt <- filtered$sample
filtered$sample <- case_when(filtered$sample_alt == "C1" ~ "PM001", filtered$sample_alt == "C5" ~ "M4666", TRUE ~ filtered$sample_alt)

p1 <- ggplot(data.frame(filtered[["umap.irna"]][[]], label = factor(filtered$sample)), aes(x=irnaUMAP_1, y=irnaUMAP_2, color=label)) + 
    ggrastr::rasterise(geom_point(size=0.2, stroke=0.1, shape=16), dpi = 400) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() + theme(legend.position = "none")
p2 <- ggplot(data.frame(filtered[["umap.iatac"]][[]], label = factor(filtered$sample)), aes(x=iatacUMAP_1, y=iatacUMAP_2, color=label)) + 
    ggrastr::rasterise(geom_point(size=0.2, stroke=0.1, shape=16), dpi = 400) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() + theme(legend.position = "none")
p3 <- ggplot(data.frame(filtered[["wnn.umap"]][[]], label = factor(filtered$sample)), aes(x=wnnUMAP_1, y=wnnUMAP_2, color=label)) + 
    ggrastr::rasterise(geom_point(size=0.2, stroke=0.1, shape=16), dpi = 400) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() + theme(legend.position = "none")

pdf(file='umaps-filtered.pdf', width=15, height=15)
plot(p1 + p2 + p3 + plot_spacer())
dev.off()

saveRDS(filtered, "filtered.rds")

# filtered <- readRDS("filtered.rds")

print(paste0(Sys.time(), ":: Cell-type identification using scRNA-seq data ... "))

scRNAref <- readRDS("/FINAL/scRNA/filtered.rds")
mutmap <- data.frame(
    sample = c("A5", "A6", "A12", "A16", "A25", "A26", "A30", "A40", "A42", "108", "890", "M3399", "M4666", "M5167", "PM001", "C2", "C3", "C4", "C6", "C8", "L010-1", "L010-2", "L027", "C7"), 
    mutation = c("2xTET2-C", "DNMT3A-C", "TET2-C", "ASXL1-C", "TET2-C", "DNMT3A-C", "TET2-C", "2xDNMT3A-C", "TET2-C", "2xTET2-C", "EZH2-C", "2xTET2", "TET2", "ASXL1", "TET2", "DNMT3A", "TET2", "TET2/DNMT3A", "SF3B1", "2xTET2", "2xTET2", "2xTET2", "2xTET2", "ASXL1"),
    mutation1 = c("TET2-C", "CHIP-C", "TET2-C", "CHIP-C", "TET2-C", "CHIP-C", "TET2-C", "CHIP-C", "TET2-C", "TET2-C", "CHIP-C", "TET2", "TET2", "CHIP", "TET2", "CHIP", "TET2", "TET2", "CHIP", "TET2", "TET2", "TET2", "TET2", "CHIP")
)
scRNAref$mutation <- (data.frame(sample = scRNAref$sample) %>% left_join(mutmap) %>% replace_na(list(mutation = "NOCHIP")))$mutation
scRNAref$mutation1 <- (data.frame(sample = scRNAref$sample) %>% left_join(mutmap) %>% replace_na(list(mutation1 = "NOCHIP")))$mutation1
scRNAref$mutation1 <- case_when((scRNAref$status == "COVID" & scRNAref$mutation1 == "NOCHIP") ~ "COVID", TRUE ~ scRNAref$mutation1)
pred.monaco <- readRDS("/FINAL/scRNA/SingleR.monaco.fine.rds")
scRNAref$monaco <- pred.monaco$new.fine

reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")

scRNAref <- SCTransform(
  object = scRNAref,
  assay = "RNA",
  new.assay.name = "refAssay",
  method = 'glmGamPoi',
  ncells = 2000,
  n_genes = 2000,
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)
filtered <- SCTransform(
  object = filtered,
  assay = "RNA",
  new.assay.name = "scRNArefAssay",
  residual.features = rownames(x = scRNAref),
  reference.SCT.model = scRNAref[["refAssay"]]@SCTModel.list$refmodel,
  method = 'glmGamPoi',
  ncells = 2000,
  n_genes = 2000,
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)
anchors <- FindTransferAnchors(
  reference = scRNAref,
  query = filtered,
  k.filter = NA,
  reference.neighbors = NULL,
  reference.assay = "refAssay",
  query.assay = "scRNArefAssay",
  reference.reduction = NULL,
  normalization.method = "SCT",
  features = intersect(rownames(x = scRNAref), VariableFeatures(object = filtered)),
  dims = 1:50,
  n.trees = 20,
  mapping.score.k = 100
)
refdata <- lapply(X = c("monaco"), function(x) {
  scRNAref[[x, drop = TRUE]]
})
names(x = refdata) <- c("scRNA.monaco")
filtered <- TransferData(
  reference = scRNAref,
  query = filtered,
  dims = 1:50,
  anchorset = anchors,
  refdata = refdata,
  n.trees = 20,
  store.weights = TRUE
)

saveRDS(filtered, "filtered.rds")

#########

meta <- read.table("/CC_multiomics/scRNA_MULT/scRNA.metadata.tsv", sep = "\t") %>% set_colnames(c("status", "sample", "age", "WHO", "CRS", "Sex"))
meta <- meta %>% mutate(sample = factor(sample, levels = (meta %>% arrange(age))$sample)) %>% arrange(age) %>% na_if(-1) %>% replace_na(list(age = "  ", WHO = "  ", CRS = "  "))

filtered$status <- (data.frame(sample = filtered$sample) %>% left_join(meta))$status
filtered$age <- (data.frame(sample = filtered$sample) %>% left_join(meta))$age
filtered$WHO <- (data.frame(sample = filtered$sample) %>% left_join(meta))$WHO
filtered$CRS <- (data.frame(sample = filtered$sample) %>% left_join(meta))$CRS
filtered$Sex <- (data.frame(sample = filtered$sample) %>% left_join(meta))$Sex

Annotation(filtered[["mATAC"]]) <- annotations

saveRDS(filtered, "filtered.rds")

###############################################################################

print(paste0(Sys.time(), ":: Adding more information to the object ... "))

mutmap <- data.frame(
    sample = c("A5", "A6", "A12", "A16", "A25", "A26", "A30", "A40", "A42", "108", "890", "M3399", "M4666", "M5167", "PM001", "C2", "C3", "C4", "C6", "C8"), 
    mutation = c("TET2", "DNMT3A", "TET2", "ASXL1", "TET2", "DNMT3A", "TET2", "DNMT3A", "TET2", "TET2", "EZH2", "TET2", "TET2", "ASXL1", "TET2", "DNMT3A", "TET2", "TET2/DNMT3A", "SF3B1", "TET2"))
filtered$mutation <- (filtered@meta.data %>% left_join(mutmap) %>% replace_na(list(mutation = "NONE")))$mutation

print(paste0(Sys.time(), ":: Adding more information to the object: 1) Adding NucleosomeSignal and TSSEnrichment ... "))

filtered <- NucleosomeSignal(object = filtered, assay = "mATAC")
filtered <- TSSEnrichment(object = filtered, fast = FALSE, assay = "mATAC")

print(paste0(Sys.time(), ":: Adding more information to the object: 2) Adding Gene Activity Scores ... "))

gene.activities <- GeneActivity(filtered, assay = "mATAC")
filtered[['GA']] <- CreateAssayObject(counts = gene.activities)
filtered <- NormalizeData(
  object = filtered,
  assay = 'GA',
  normalization.method = 'LogNormalize',
  scale.factor = median(filtered$nCount_GA)
)

print(paste0(Sys.time(), ":: Adding more information to the object: 3) Adding Motifs ... "))

library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
keepBSgenomeSequences <- function(genome, seqnames)
{
    stopifnot(all(seqnames %in% seqnames(genome)))
    genome@user_seqnames <- setNames(seqnames, seqnames)
    genome@seqinfo <- genome@seqinfo[seqnames]
    genome
}
sequences_to_keep <- paste0("chr", c(1:22, "X", "Y"))
genome <- keepBSgenomeSequences(genome, sequences_to_keep)

DefaultAssay(filtered) <- "mATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
filtered <- AddMotifs(filtered, genome = genome, pfm = pwm_set)

write.table(data.frame(peak = rownames(filtered)) %>% separate(peak, into = c("chr", "start", "end"), sep = "-"), 
    file = "MULT_peaks.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

saveRDS(filtered, "filtered.rds")

############################################################################################

print(paste0(Sys.time(), ":: Adding more information to the object: 4) Adding ChromVAR motif scores ... "))

library(BiocParallel)
register(SerialParam())

filtered <- RunChromVAR(
  object = filtered,
  genome = genome,
  assay = "mATAC"
)

saveRDS(filtered, "filtered.rds")

############################################################################################

print(paste0(Sys.time(), ":: Adding more information to the object: 5) Adding ChromVAR ReMap2022 scores ... "))

library(BiocParallel)
register(SerialParam())
MULT_peaks <- data.frame(Peak = rownames(filtered)) %>% separate(col = "Peak", sep = "-", into = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()

remap_mo <- read.table("/public_data/remap2022/remap2022_monocyte_nr_macs2_hg38_v1_0.bed") %>% 
    select(V1:V4) %>% 
    set_colnames(c("chr", "start", "end", "name")) %>%
    separate(col = "name", sep = ":", into = c("TF", "biotype"))
TFs <- unique(remap_mo$TF)
TF_features <- lapply(TFs, function(cur_TF){
    cur_peaks <- makeGRangesFromDataFrame(remap_mo %>% filter(TF == cur_TF))
    cur_peaks <- cur_peaks[as.vector(seqnames(cur_peaks) %in% standardChromosomes(cur_peaks))]
    ovp_features <- suppressWarnings(GRangesToString(MULT_peaks[case_when(is.na(findOverlaps(MULT_peaks, cur_peaks, select = "first", ignore.strand=TRUE)) ~ FALSE, TRUE ~ TRUE)]))
    return(ovp_features)
})
names(TF_features) <- paste0("Mo.", TFs)
remap_mo <- NULL
filtered <- AddChromatinModule(filtered, features = TF_features, genome = genome, assay = "mATAC")

remap_mac <- read.table("/public_data/remap2022/remap2022_macrophage_nr_macs2_hg38_v1_0.bed") %>% 
    select(V1:V4) %>% 
    set_colnames(c("chr", "start", "end", "name")) %>%
    separate(col = "name", sep = ":", into = c("TF", "biotype"))
TFs <- unique(remap_mac$TF)
TF_features <- lapply(TFs, function(cur_TF){
    cur_peaks <- makeGRangesFromDataFrame(remap_mac %>% filter(TF == cur_TF))
    cur_peaks <- cur_peaks[as.vector(seqnames(cur_peaks) %in% standardChromosomes(cur_peaks))]
    ovp_features <- suppressWarnings(GRangesToString(MULT_peaks[case_when(is.na(findOverlaps(MULT_peaks, cur_peaks, select = "first", ignore.strand=TRUE)) ~ FALSE, TRUE ~ TRUE)]))
    return(ovp_features)
})
names(TF_features) <- paste0("Mac.", TFs)
remap_mac <- NULL
filtered <- AddChromatinModule(filtered, features = TF_features, genome = genome, assay = "mATAC")

saveRDS(filtered, "filtered.rds")

###############################################################################
print(paste0(Sys.time(), ":: Create MULTPLUSHANIFFA ... "))

haniffa <- readRDS("/HaniffaData/HaniffaAllDiet.rds")
haniffa <- subset(haniffa, subset = status == "Healthy")
pred.haniffa.monaco <- readRDS("/HaniffaData/haniffa.healthy.SingleR.rds")
haniffa$monaco <- pred.haniffa.monaco$new.fine
haniffa$group <- case_when(haniffa$Age_interval %in% c("(20, 29]", "(30, 39]", "(40, 49]") ~ "CTRLY", TRUE ~ "CTRLO")
haniffa$mutation <- case_when(haniffa$Age_interval %in% c("(20, 29]", "(30, 39]", "(40, 49]") ~ "CTRLY", TRUE ~ "CTRLO")

MULT <- readRDS("filtered.rds")
DefaultAssay(MULT) <- "RNA"
mutmap <- data.frame(
    sample = c("A5", "A6", "A12", "A16", "A25", "A26", "A30", "A40", "A42", "108", "890", "M3399", "M4666", "M5167", "PM001", "C2", "C3", "C4", "C6", "C8", "L010-1", "L010-2", "L027", "C7"), 
    group = c("TET2C", "CHIPC", "TET2C", "CHIPC", "TET2C", "CHIPC", "TET2C", "CHIPC", "TET2C", "TET2C", "CHIPC", "TET2", "TET2", "CHIP", "TET2", "CHIP", "TET2", "TET2", "CHIP", "TET2", "TET2", "TET2", "TET2", "CHIP"),
    mutation = c("2xTET2-C", "DNMT3A-C", "TET2-C", "ASXL1-C", "TET2-C", "DNMT3A-C", "TET2-C", "2xDNMT3A-C", "TET2-C", "2xTET2-C", "EZH2-C", "2xTET2", "TET2", "ASXL1", "TET2", "DNMT3A", "TET2", "TET2/DNMT3A", "SF3B1", "2xTET2", "2xTET2", "2xTET2", "2xTET2", "ASXL1")
)
MULT$group <- (data.frame(sample = MULT$sample) %>% left_join(mutmap) %>% replace_na(list(group = "COVID")))$group
MULT$mutation <- (data.frame(sample = MULT$sample) %>% left_join(mutmap) %>% replace_na(list(mutation = "COVID")))$mutation
MULT$monaco <- MULT$predicted.scRNA.monaco

haniffa@meta.data <- haniffa@meta.data %>% select(Sex, nFeature_RNA, nCount_RNA, monaco, Age_interval, sample, group, mutation)
MULT@meta.data <- MULT@meta.data %>% select(Sex, nFeature_RNA, nCount_RNA, monaco, age, sample, group, mutation)
haniffa <- DietSeurat(haniffa, assays = "RNA", features = intersect(rownames(haniffa), rownames(MULT)))
MULT <- DietSeurat(MULT, assays = "RNA", features = intersect(rownames(haniffa), rownames(MULT)))

MULTPLUSHANIFFA <- merge(MULT, y = haniffa)
MULTPLUSHANIFFA <- MULTPLUSHANIFFA %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000)

saveRDS(MULTPLUSHANIFFA, "MULTPLUSHANIFFA.rds")



