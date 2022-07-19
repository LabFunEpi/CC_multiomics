setwd("/FINAL/scMULT")
filtered <- readRDS("/FINAL/scMULT/filtered.rds")

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
genome_lengths <- seqlengths(genome)
genome.df <- data.frame("chr" = names(genome_lengths), "length" = genome_lengths)

mutmap <- data.frame(
    sample = c("A5", "A6", "A12", "A16", "A25", "A26", "A30", "A40", "A42", "108", "890", "M3399", "M4666", "M5167", "PM001", "C2", "C3", "C4", "C6", "C8", "L010-1", "L010-2", "L027", "C7"), 
    mutation = c("2xTET2-C", "DNMT3A-C", "TET2-C", "ASXL1-C", "TET2-C", "DNMT3A-C", "TET2-C", "2xDNMT3A-C", "TET2-C", "2xTET2-C", "EZH2-C", "2xTET2", "TET2", "ASXL1", "TET2", "DNMT3A", "TET2", "TET2/DNMT3A", "SF3B1", "2xTET2", "2xTET2", "2xTET2", "2xTET2", "ASXL1"),
    mutation1 = c("TET2-C", "CHIP-C", "TET2-C", "CHIP-C", "TET2-C", "CHIP-C", "TET2-C", "CHIP-C", "TET2-C", "TET2-C", "CHIP-C", "TET2", "TET2", "CHIP", "TET2", "CHIP", "TET2", "TET2", "CHIP", "TET2", "TET2", "TET2", "TET2", "CHIP")
)
filtered$mutation1 <- (data.frame(sample = filtered$sample) %>% left_join(mutmap) %>% replace_na(list(mutation1 = "COVID")))$mutation1

statuses <- c("COVID", "TET2-C")
celltypes <- c("CD14 Mono")
celltype.status.levels <- apply(expand.grid(celltypes, statuses) %>% mutate(Var1 = factor(Var1, levels = celltypes)) %>% arrange(Var1), 1, paste, collapse=".")
filtered$celltype.status <- paste(filtered$predicted.scRNA.monaco, filtered$mutation1, sep = ".")

conns.list <- lapply(celltype.status.levels, function(x){
    MULT_sub <- subset(filtered, subset = celltype.status == x)
    DefaultAssay(MULT_sub) <- "mATAC"
    cds <- as.cell_data_set(MULT_sub)
    cicero.obj <- make_cicero_cds(cds, reduced_coordinates = reducedDims(cds)$WNN.UMAP)
    conns <- run_cicero(cicero.obj, genomic_coords = genome.df, sample_num = 100)
    return(conns)
})
names(conns.list) <- celltype.status.levels

saveRDS(conns.list, "conns.list.CD14Mono_COVID_TET2C.rds")



















