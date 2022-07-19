setwd("/FINAL/scRNA")


mutmap <- data.frame(
    sample = c("A5", "A6", "A12", "A16", "A25", "A26", "A30", "A40", "A42", "108", "890", "M3399", "M4666", "M5167", "PM001", "C2", "C3", "C4", "C6", "C8", "L010-1", "L010-2", "L027", "C7"), 
    mutation = c("TET2-C", "DNMT3A-C", "TET2-C", "ASXL1-C", "TET2-C", "DNMT3A-C", "TET2-C", "DNMT3A-C", "TET2-C", "TET2-C", "EZH2-C", "TET2", "TET2", "ASXL1", "TET2", "DNMT3A", "TET2", "TET2", "SF3B1", "TET2", "TET2", "TET2", "TET2", "ASXL1"),
    mutation1 = c("TET2-C", "CHIP-C", "TET2-C", "CHIP-C", "TET2-C", "CHIP-C", "TET2-C", "CHIP-C", "TET2-C", "TET2-C", "CHIP-C", "TET2", "TET2", "CHIP", "TET2", "CHIP", "TET2", "TET2", "CHIP", "TET2", "TET2", "TET2", "TET2", "CHIP")
)

GEX <- readRDS("filtered.rds")
GEX$mutation <- (data.frame(sample = GEX$sample) %>% left_join(mutmap) %>% replace_na(list(mutation = "COVID")))$mutation
GEX$mutation1 <- (data.frame(sample = GEX$sample) %>% left_join(mutmap) %>% replace_na(list(mutation1 = "COVID")))$mutation1

pred.GEX.monaco <- readRDS("SingleR.monaco.fine.rds")
GEX$monaco <- pred.GEX.monaco$new.fine

GEX$celltype.status <- paste(GEX$monaco, GEX$mutation, sep = ".")
Idents(GEX) <- "celltype.status"
DefaultAssay(GEX) <- "RNA"


print("1) Comparing DNMT3A-C vs COVID ...")

statuses <- c("DNMT3A-C", "COVID")

types <- unique(GEX$monaco)
samples <- GEX@meta.data %>% filter(mutation %in% statuses) %>% select(sample) %>% unlist() %>% unname() %>% unique()

degs <- lapply(
    types, 
    FUN = function(celltype){
        curridents <- unique(Idents(GEX))
        if (!(paste0(celltype, ".", statuses[[1]]) %in% curridents) | !(paste0(celltype, ".", statuses[[2]]) %in% curridents)) return(data.frame())
        if (length(WhichCells(GEX, idents = paste0(celltype, ".", statuses[[1]]))) < 3 | length(WhichCells(GEX, idents = paste0(celltype, ".", statuses[[2]]))) < 3) return(data.frame())
        temp1 <- FindMarkers(GEX, ident.1 = paste0(celltype, ".", statuses[[1]]), ident.2 = paste0(celltype, ".", statuses[[2]]), min.pct = 0) %>% rownames_to_column("gene")
        if(is.na(temp1)){ return(data.frame()) } else { return(temp1) }
    }
)
names(degs) <- types
degs <- bind_rows(degs, .id = "celltype")
write.table(degs, file = paste0("degs_", statuses[[1]], "_", statuses[[2]], ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

degs_robust <- lapply(
    samples, 
    FUN = function(sam){
        print(paste0("\t", sam))
        oneout <- subset(GEX, subset = sample != sam)
        curridents <- unique(Idents(oneout))
        ctsp_degs_list <- lapply(
            types,
            FUN = function(celltype){
                print(paste0("\t\t", celltype))
                if (!(paste0(celltype, ".", statuses[[1]]) %in% curridents) | !(paste0(celltype, ".", statuses[[2]]) %in% curridents)) return(c())
                if (length(WhichCells(oneout, idents = paste0(celltype, ".", statuses[[1]]))) < 3 | length(WhichCells(oneout, idents = paste0(celltype, ".", statuses[[2]]))) < 3) return(c())
                temp1 <- FindMarkers(oneout, ident.1 = paste0(celltype, ".", statuses[[1]]), ident.2 = paste0(celltype, ".", statuses[[2]]), min.pct = 0) %>% rownames_to_column("gene")
                if(is.na(temp1)){ return(c()) } else { return((temp1 %>% filter(p_val_adj < 0.05))$gene) }
            }
        )
        names(ctsp_degs_list) <- types
        return(ctsp_degs_list)
    }
)
names(degs_robust) <- samples
degs_robust <- lapply(degs_robust %>% transpose(), function(x){Reduce(intersect, x)})

degsr <- bind_rows(lapply(degs_robust, function(gene){data.frame(gene)}), .id = "celltype") %>% left_join(degs)
write.table(degsr, file = paste0("degsr_", statuses[[1]], "_", statuses[[2]], ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#####################################

print("2) Comparing ASXL1-C vs COVID ...")

statuses <- c("ASXL1-C", "COVID")

types <- unique(GEX$monaco)
samples <- GEX@meta.data %>% filter(mutation %in% c("COVID")) %>% select(sample) %>% unlist() %>% unname() %>% unique()

degs <- lapply(
    types, 
    FUN = function(celltype){
        curridents <- unique(Idents(GEX))
        if (!(paste0(celltype, ".", statuses[[1]]) %in% curridents) | !(paste0(celltype, ".", statuses[[2]]) %in% curridents)) return(data.frame())
        if (length(WhichCells(GEX, idents = paste0(celltype, ".", statuses[[1]]))) < 3 | length(WhichCells(GEX, idents = paste0(celltype, ".", statuses[[2]]))) < 3) return(data.frame())
        temp1 <- FindMarkers(GEX, ident.1 = paste0(celltype, ".", statuses[[1]]), ident.2 = paste0(celltype, ".", statuses[[2]]), min.pct = 0) %>% rownames_to_column("gene")
        if(is.na(temp1)){ return(data.frame()) } else { return(temp1) }
    }
)
names(degs) <- types
degs <- bind_rows(degs, .id = "celltype")
write.table(degs, file = paste0("degs_", statuses[[1]], "_", statuses[[2]], ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

degs_robust <- lapply(
    samples, 
    FUN = function(sam){
        print(paste0("\t", sam))
        oneout <- subset(GEX, subset = sample != sam)
        curridents <- unique(Idents(oneout))
        ctsp_degs_list <- lapply(
            types,
            FUN = function(celltype){
                print(paste0("\t\t", celltype))
                if (!(paste0(celltype, ".", statuses[[1]]) %in% curridents) | !(paste0(celltype, ".", statuses[[2]]) %in% curridents)) return(c())
                if (length(WhichCells(oneout, idents = paste0(celltype, ".", statuses[[1]]))) < 3 | length(WhichCells(oneout, idents = paste0(celltype, ".", statuses[[2]]))) < 3) return(c())
                temp1 <- FindMarkers(oneout, ident.1 = paste0(celltype, ".", statuses[[1]]), ident.2 = paste0(celltype, ".", statuses[[2]]), min.pct = 0) %>% rownames_to_column("gene")
                if(is.na(temp1)){ return(c()) } else { return((temp1 %>% filter(p_val_adj < 0.05))$gene) }
            }
        )
        names(ctsp_degs_list) <- types
        return(ctsp_degs_list)
    }
)
names(degs_robust) <- samples
degs_robust <- lapply(degs_robust %>% transpose(), function(x){Reduce(intersect, x)})

degsr <- bind_rows(lapply(degs_robust, function(gene){data.frame(gene)}), .id = "celltype") %>% left_join(degs)
write.table(degsr, file = paste0("degsr_", statuses[[1]], "_", statuses[[2]], ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#####################################

print("3) Comparing EZH2-C vs COVID ...")

statuses <- c("EZH2-C", "COVID")

types <- unique(GEX$monaco)
samples <- GEX@meta.data %>% filter(mutation %in% c("COVID")) %>% select(sample) %>% unlist() %>% unname() %>% unique()

degs <- lapply(
    types, 
    FUN = function(celltype){
        curridents <- unique(Idents(GEX))
        if (!(paste0(celltype, ".", statuses[[1]]) %in% curridents) | !(paste0(celltype, ".", statuses[[2]]) %in% curridents)) return(data.frame())
        if (length(WhichCells(GEX, idents = paste0(celltype, ".", statuses[[1]]))) < 3 | length(WhichCells(GEX, idents = paste0(celltype, ".", statuses[[2]]))) < 3) return(data.frame())
        temp1 <- FindMarkers(GEX, ident.1 = paste0(celltype, ".", statuses[[1]]), ident.2 = paste0(celltype, ".", statuses[[2]]), min.pct = 0) %>% rownames_to_column("gene")
        if(is.na(temp1)){ return(data.frame()) } else { return(temp1) }
    }
)
names(degs) <- types
degs <- bind_rows(degs, .id = "celltype")
write.table(degs, file = paste0("degs_", statuses[[1]], "_", statuses[[2]], ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

degs_robust <- lapply(
    samples, 
    FUN = function(sam){
        print(paste0("\t", sam))
        oneout <- subset(GEX, subset = sample != sam)
        curridents <- unique(Idents(oneout))
        ctsp_degs_list <- lapply(
            types,
            FUN = function(celltype){
                print(paste0("\t\t", celltype))
                if (!(paste0(celltype, ".", statuses[[1]]) %in% curridents) | !(paste0(celltype, ".", statuses[[2]]) %in% curridents)) return(c())
                if (length(WhichCells(oneout, idents = paste0(celltype, ".", statuses[[1]]))) < 3 | length(WhichCells(oneout, idents = paste0(celltype, ".", statuses[[2]]))) < 3) return(c())
                temp1 <- FindMarkers(oneout, ident.1 = paste0(celltype, ".", statuses[[1]]), ident.2 = paste0(celltype, ".", statuses[[2]]), min.pct = 0) %>% rownames_to_column("gene")
                if(is.na(temp1)){ return(c()) } else { return((temp1 %>% filter(p_val_adj < 0.05))$gene) }
            }
        )
        names(ctsp_degs_list) <- types
        return(ctsp_degs_list)
    }
)
names(degs_robust) <- samples
degs_robust <- lapply(degs_robust %>% transpose(), function(x){Reduce(intersect, x)})

degsr <- bind_rows(lapply(degs_robust, function(gene){data.frame(gene)}), .id = "celltype") %>% left_join(degs)
write.table(degsr, file = paste0("degsr_", statuses[[1]], "_", statuses[[2]], ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

############################################

GEX$celltype.status <- paste(GEX$monaco, GEX$status, sep = ".")
Idents(GEX) <- "celltype.status"
DefaultAssay(GEX) <- "RNA"

print("4) Comparing COVIDCHIP vs COVID ...")

statuses <- c("COVIDCHIP", "COVID")

types <- unique(GEX$monaco)
samples <- GEX@meta.data %>% filter(status %in% statuses) %>% select(sample) %>% unlist() %>% unname() %>% unique()

degs <- lapply(
    types, 
    FUN = function(celltype){
        curridents <- unique(Idents(GEX))
        if (!(paste0(celltype, ".", statuses[[1]]) %in% curridents) | !(paste0(celltype, ".", statuses[[2]]) %in% curridents)) return(data.frame())
        if (length(WhichCells(GEX, idents = paste0(celltype, ".", statuses[[1]]))) < 3 | length(WhichCells(GEX, idents = paste0(celltype, ".", statuses[[2]]))) < 3) return(data.frame())
        temp1 <- FindMarkers(GEX, ident.1 = paste0(celltype, ".", statuses[[1]]), ident.2 = paste0(celltype, ".", statuses[[2]]), min.pct = 0) %>% rownames_to_column("gene")
        if(is.na(temp1)){ return(data.frame()) } else { return(temp1) }
    }
)
names(degs) <- types
degs <- bind_rows(degs, .id = "celltype")
write.table(degs, file = paste0("degs_", statuses[[1]], "_", statuses[[2]], ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

degs_robust <- lapply(
    samples, 
    FUN = function(sam){
        print(paste0("\t", sam))
        oneout <- subset(GEX, subset = sample != sam)
        curridents <- unique(Idents(oneout))
        ctsp_degs_list <- lapply(
            types,
            FUN = function(celltype){
                print(paste0("\t\t", celltype))
                if (!(paste0(celltype, ".", statuses[[1]]) %in% curridents) | !(paste0(celltype, ".", statuses[[2]]) %in% curridents)) return(c())
                if (length(WhichCells(oneout, idents = paste0(celltype, ".", statuses[[1]]))) < 3 | length(WhichCells(oneout, idents = paste0(celltype, ".", statuses[[2]]))) < 3) return(c())
                temp1 <- FindMarkers(oneout, ident.1 = paste0(celltype, ".", statuses[[1]]), ident.2 = paste0(celltype, ".", statuses[[2]]), min.pct = 0) %>% rownames_to_column("gene")
                if(is.na(temp1)){ return(c()) } else { return((temp1 %>% filter(p_val_adj < 0.05))$gene) }
            }
        )
        names(ctsp_degs_list) <- types
        return(ctsp_degs_list)
    }
)
names(degs_robust) <- samples
degs_robust <- lapply(degs_robust %>% transpose(), function(x){Reduce(intersect, x)})

degsr <- bind_rows(lapply(degs_robust, function(gene){data.frame(gene)}), .id = "celltype") %>% left_join(degs)
write.table(degsr, file = paste0("degsr_", statuses[[1]], "_", statuses[[2]], ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
