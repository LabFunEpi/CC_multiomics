setwd("/FINAL/scRNA")

statuses <- c("TET2C", "TET2")

MAYOPLUSHANIFFA <- readRDS("MAYOPLUSHANIFFA.rds")

MAYOPLUSHANIFFA$celltype.status <- paste(MAYOPLUSHANIFFA$monaco, MAYOPLUSHANIFFA$group, sep = ".")
Idents(MAYOPLUSHANIFFA) <- "celltype.status"
DefaultAssay(MAYOPLUSHANIFFA) <- "RNA"

types <- unique(MAYOPLUSHANIFFA$monaco)
samples <- MAYOPLUSHANIFFA@meta.data %>% filter(group %in% statuses) %>% select(sample) %>% unlist() %>% unname() %>% unique()

degs <- lapply(
    types, 
    FUN = function(celltype){
        curridents <- unique(Idents(MAYOPLUSHANIFFA))
        if (!(paste0(celltype, ".", statuses[[1]]) %in% curridents) | !(paste0(celltype, ".", statuses[[2]]) %in% curridents)) return(data.frame())
        temp1 <- FindMarkers(MAYOPLUSHANIFFA, ident.1 = paste0(celltype, ".", statuses[[1]]), ident.2 = paste0(celltype, ".", statuses[[2]]), min.pct = 0) %>% rownames_to_column("gene")
        if(is.na(temp1)){ return(data.frame()) } else { return(temp1) }
    }
)
names(degs) <- types
degs <- bind_rows(degs, .id = "celltype")
write.table(degs, file = paste0("degs_", statuses[[1]], "_", statuses[[2]], ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


degs_robust <- lapply(
    samples, 
    FUN = function(sam){
        print(sam)
        oneout <- subset(MAYOPLUSHANIFFA, subset = sample != sam)
        curridents <- unique(Idents(oneout))
        ctsp_degs_list <- lapply(
            types,
            FUN = function(celltype){
                print(celltype)
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
