
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

setwd("/FINAL")
colors.mn <- c(`CD4 T` = "#99CCCC", `Treg` = "#E794EA", `CD8 T` = "#A9E0AB", `gdT` = "#CE5B5B", `MAIT` = "#EDA6A6", `NK` = "#E2E8A6", `B` = "#D68645", `Plasmablast` = "#48B758", `CD14 Mono` = "#B4AEE5", `CD16 Mono` = "#F9BD95", `Int Mono` = "#B81254", `cDC` = "#67A5E5", `pDC` = "#D8E1A7", `HSPC` = "#C683ED", `Neutrophils` = "#3B69B7", `Basophils` = "black")

mutmap <- data.frame(
    sample = c("A5", "A6", "A12", "A16", "A25", "A26", "A30", "A40", "A42", "108", "890", "M3399", "M4666", "M5167", "PM001", "C2", "C3", "C4", "C6", "C8", "L010-1", "L010-2", "L027", "C7"), 
    mutation = c("2xTET2-C", "DNMT3A-C", "TET2-C", "ASXL1-C", "TET2-C", "DNMT3A-C", "TET2-C", "2xDNMT3A-C", "TET2-C", "2xTET2-C", "EZH2-C", "2xTET2", "TET2", "ASXL1", "TET2", "DNMT3A", "TET2", "TET2/DNMT3A", "SF3B1", "2xTET2", "2xTET2", "2xTET2", "2xTET2", "ASXL1"),
    mutation1 = c("TET2-C", "CHIP-C", "TET2-C", "CHIP-C", "TET2-C", "CHIP-C", "TET2-C", "CHIP-C", "TET2-C", "TET2-C", "CHIP-C", "TET2", "TET2", "CHIP", "TET2", "CHIP", "TET2", "TET2", "CHIP", "TET2", "TET2", "TET2", "TET2", "CHIP")
)

GEX <- readRDS("/FINAL/scRNA/filtered.rds")
GEX$mutation <- (data.frame(sample = GEX$sample) %>% left_join(mutmap) %>% replace_na(list(mutation = "COVID")))$mutation
GEX$mutation1 <- (data.frame(sample = GEX$sample) %>% left_join(mutmap) %>% replace_na(list(mutation1 = "COVID")))$mutation1
GEX$mutation2 <- recode(GEX$mutation1, !!!c(`COVID` = "COVID (no CHIP)", `CHIP` = "Non-TET2 CHIP", `TET2` = "TET2 CHIP", `CHIP-C` = "Non-TET2 CHIP + COVID", `TET2-C` = "TET2 CHIP + COVID"))
GEX$mutation3 <- recode(GEX$mutation1, !!!c(`COVID` = "COVID_NoCHIP", `CHIP` = "NonTET2CHIP", `TET2` = "TET2CHIP", `CHIP-C` = "COVID_NonTET2CHIP", `TET2-C` = "COVID_TET2CHIP"))

pred.GEX.monaco <- readRDS("/FINAL/scRNA/SingleR.monaco.fine.rds")
GEX$monaco <- pred.GEX.monaco$new.fine
level_key <- as.character(1:44)
names(level_key) <- c("C2", "M4666", "C8", "M3399", "C3", "PM001", "M5167", "C4", "072", "984", "B36", "B23", "B1", "B39", "B33", "454", "B31", "O23", "789", "B34", "B38", "B2", "570", "120", "B32", "B30", "O31", "O22", "B7", "890", "A6", "A25", "A42", "A16", "O21", "108", "A5", "A12", "A40", "B07", "A30", "A26", "O24", "O11"
)
GEX$manusample <- recode_factor(GEX$sample, !!!level_key)

MULT <- readRDS("/FINAL/scMULT/filtered.rds")
MULT$mutation <- (data.frame(sample = MULT$sample) %>% left_join(mutmap) %>% replace_na(list(mutation = "COVID")))$mutation
MULT$mutation1 <- (data.frame(sample = MULT$sample) %>% left_join(mutmap) %>% replace_na(list(mutation1 = "COVID")))$mutation1
MULT$mutation2 <- recode(MULT$mutation1, !!!c(`COVID` = "COVID (no CHIP)", `CHIP` = "Non-TET2 CHIP", `TET2` = "TET2 CHIP", `CHIP-C` = "Non-TET2 CHIP + COVID", `TET2-C` = "TET2 CHIP + COVID"))

meta <- read.table("/CC_multiomics/scRNA_MULT/scRNA.metadata.tsv", sep = "\t") %>% set_colnames(c("status", "sample", "age", "WHO", "CRS", "Sex"))
meta <- meta %>% mutate(sample = factor(sample, levels = (meta %>% arrange(age))$sample)) %>% arrange(age) %>% na_if(-1) %>% mutate(age = as.character(age), WHO = as.character(WHO), CRS = as.character(CRS)) %>% replace_na(list(age = "  ", WHO = "  ", CRS = "  "))

haniffa <- readRDS("/HaniffaData/HaniffaAllDiet.rds")
haniffa <- subset(haniffa, subset = status == "Healthy")
pred.haniffa.monaco <- readRDS("/HaniffaData/haniffa.healthy.SingleR.rds")
haniffa$monaco <- pred.haniffa.monaco$new.fine
haniffa$mutation2 <- case_when(haniffa$Age_interval %in% c("(20, 29]", "(30, 39]", "(40, 49]") ~ "Healthy (< 50)", TRUE ~ "Healthy (> 50)")


########################  UMAPS  ##########################


four_umaps <- function(rnaobj, multobj){
    p1 <- ggplot(data.frame(rnaobj[["umap"]][[]], label = factor(rnaobj$monaco)), aes(x=UMAP_1, y=UMAP_2, color=label)) + 
        ggrastr::rasterise(geom_point(size=0.3, stroke=0.1, shape=16), dpi = 400) +
        scale_color_manual(values = colors.mn) +
        guides(color = guide_legend(override.aes = list(size=5))) +
        theme_cowplot() + theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
    p2 <- table(rnaobj$monaco) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn)))) %>%
        ggplot(aes(fill=Var1, y=Freq, x="Blah")) +
        geom_bar(position="fill", stat="identity", width = 0.6) + 
        scale_fill_manual(values = colors.mn) +
        coord_flip() + theme_cowplot() + 
        theme(aspect.ratio=0.1, legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
    p3 <- ggplot(data.frame(multobj[["umap.irna"]][[]], label = factor(multobj$predicted.scRNA.monaco)), aes(x=-irnaUMAP_1, y=-irnaUMAP_2, color=label)) + 
        ggrastr::rasterise(geom_point(size=0.5, stroke=0.1, shape=16), dpi = 400) +
        scale_color_manual(values = colors.mn) +
        guides(color = guide_legend(override.aes = list(size=5))) +
        theme_cowplot() + theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
    p4 <- table(multobj$predicted.scRNA.monaco) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn)))) %>%
        ggplot(aes(fill=Var1, y=Freq, x="Blah")) +
        geom_bar(position="fill", stat="identity", width = 0.6) + 
        scale_fill_manual(values = colors.mn) +
        coord_flip() + theme_cowplot() + 
        theme(aspect.ratio=0.1, legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
    p5 <- ggplot(data.frame(multobj[["umap.iatac"]][[]], label = factor(multobj$predicted.scRNA.monaco)), aes(x=iatacUMAP_1, y=iatacUMAP_2, color=label)) + 
        ggrastr::rasterise(geom_point(size=0.5, stroke=0.1, shape=16), dpi = 400) +
        scale_color_manual(values = colors.mn) +
        guides(color = guide_legend(override.aes = list(size=5))) +
        theme_cowplot() + theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
    p6 <- ggplot(data.frame(multobj[["wnn.umap"]][[]], label = factor(multobj$predicted.scRNA.monaco)), aes(x=wnnUMAP_1, y=wnnUMAP_2, color=label)) + 
        ggrastr::rasterise(geom_point(size=0.5, stroke=0.1, shape=16), dpi = 400) +
        scale_color_manual(values = colors.mn) +
        guides(color = guide_legend(override.aes = list(size=5))) +
        theme_cowplot() + theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
    return((p1 / p2) | (p3 / p4) | (p5 / p4) | (p6 / p4))
}

pdf(file='umaps.pdf', width=20, height=6)
plot(four_umaps(GEX, MULT))
dev.off()

celltype.table1 <- data.frame(celltype = GEX$monaco, sample = GEX$sample) %>% table() %>% data.frame() %>%
    rename(count = Freq) %>%
    mutate(celltype = factor(celltype, levels=names(colors.mn))) %>%
    group_by(sample) %>%
    left_join(GEX@meta.data %>% select(sample, mutation1) %>% distinct()) %>%
    filter(mutation1 %in% c("COVID", "TET2-C"))

celltype.table2 <- celltype.table1 %>% group_by(celltype, mutation1) %>% summarize(count = sum(count))

celltype.table3 <- celltype.table1 %>% group_by(celltype, mutation1) %>% summarize(count = mean(count))
celltype.table3 <- celltype.table3 %>% left_join(celltype.table3 %>% group_by(mutation1) %>% summarize(total = sum(count))) %>% mutate(freq = count/total) %>% filter(celltype %in% c("CD14 Mono", "CD16 Mono", "Int Mono")) %>% ungroup()

celltype.table3 %>% select(-c(total, freq)) %>% pivot_wider(names_from = "mutation1", values_from = "count") %>% column_to_rownames("celltype")
chisq.test(celltype.table3 %>% select(-c(total, freq)) %>% pivot_wider(names_from = "mutation1", values_from = "count") %>% column_to_rownames("celltype"))


pdf(file='umap.pdf', width=10, height=8)
p1 <- ggplot(data.frame(GEX[["umap"]][[]], label = factor(GEX$monaco)), aes(x=UMAP_1, y=UMAP_2, color=label)) + 
    ggrastr::rasterise(geom_point(size=0.2, stroke=0.1, shape=16), dpi = 400) +
    scale_color_manual(values = colors.mn) +
    annotate("text", label = paste0("n = ", nrow(GEX[["umap"]][[]])), x = -2, y = 12, size = 8) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() + theme(legend.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
p2 <- table(GEX$monaco) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn)))) %>%
    ggplot(aes(fill=Var1, y=Freq, x="Blah")) +
    geom_bar(position="fill", stat="identity", width = 0.3) + 
    scale_fill_manual(values = colors.mn) +
    coord_flip() + theme_cowplot() + 
    theme(aspect.ratio=0.1, legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
p1 / p2
dev.off()

pdf(file='umap-others.pdf', width=10, height=8)
p1 <- ggplot(data.frame(GEX[["umap"]][[]], label = factor(GEX$manusample)), aes(x=UMAP_1, y=UMAP_2, color=label)) + 
    ggrastr::rasterise(geom_point(size=0.2, stroke=0.1, shape=16), dpi = 400) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() + theme(legend.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
p1
p1 <- ggplot(data.frame(GEX[["umap"]][[]], label = factor(GEX$status)), aes(x=UMAP_1, y=UMAP_2, color=label)) + 
    ggrastr::rasterise(geom_point(size=0.2, stroke=0.1, shape=16), dpi = 400) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() + theme(legend.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
p1
dev.off()


##################### Celltype proportions ################################


levels.status <- c("TET2CHIP", "NonTET2CHIP", "COVID_NoCHIP", "COVID_TET2CHIP", "COVID_NonTET2CHIP")
celltype.table1 <- data.frame(celltype = GEX$monaco, sample = GEX$sample) %>% table() %>% data.frame() %>%
    rename(count = Freq) %>%
    mutate(celltype = factor(celltype, levels=names(colors.mn))) %>%
    mutate(sample = factor(sample, levels = (meta %>% arrange(age))$sample)) %>%
    left_join(data.frame(sample = GEX$sample, status = factor(GEX$mutation3, levels = levels.status)) %>% distinct()) %>%
    group_by(sample) %>%
    mutate(freq = count / sum(count)) %>%
    left_join(GEX@meta.data %>% select(sample, manusample) %>% distinct())
meta1 <- meta %>% left_join(GEX$sample %>% table() %>% data.frame() %>% set_colnames(c("sample", "total"))) %>% filter(!(status %in% c("HEALTHY", "CMML")))

celltype.table2 <- data.frame(celltype = GEX$monaco, status = GEX$mutation3) %>% table() %>% data.frame() %>%
    rename(count = Freq) %>%
    mutate(celltype = factor(celltype, levels=names(colors.mn)))

p1 <- ggplot(celltype.table1) + 
    geom_bar(aes(fill=celltype, y=count, x=manusample), position=position_fill(reverse = TRUE), stat="identity", color = "black") +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = colors.mn) +
    facet_grid(rows = vars(status), scales = "free", space = "free") +
    coord_flip() + 
    theme_cowplot() + 
    theme(legend.position = "None", strip.background = element_blank(), strip.text = element_blank())
p3 <- ggplot(celltype.table2) + 
    geom_bar(aes(fill=celltype, y=count, x=factor(status, levels = rev(levels.status))), position=position_fill(reverse = TRUE), stat="identity", color = "black") +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = colors.mn) +
    coord_flip() + 
    theme_cowplot() + 
    theme(legend.position = "bottom", axis.title.y = element_blank(), axis.title.x = element_blank())

pdf(file='celltype-proportions-allsamples.pdf', width=18, height=10)
layout <- "
AAAAAAAAAAAAAA
AAAAAAAAAAAAAA
AAAAAAAAAAAAAA
AAAAAAAAAAAAAA
AAAAAAAAAAAAAA
AAAAAAAAAAAAAA
AAAAAAAAAAAAAA
BBBBBBBBBBBBBB
BBBBBBBBBBBBBB
"
plot(p1 + p3 + plot_layout(design = layout))
dev.off()


levels.status <- c("Healthy (< 50)", "Healthy (> 50)", "Non-TET2 CHIP", "TET2 CHIP", "COVID (no CHIP)", "Non-TET2 CHIP + COVID", "TET2 CHIP + COVID")
celltype.table1 <- data.frame(celltype = GEX$monaco, status = GEX$mutation2) %>% table() %>% data.frame() %>%
    rename(count = Freq) %>%
    mutate(celltype = factor(celltype, levels=names(colors.mn))) %>%
    bind_rows(data.frame(celltype = haniffa$monaco, status = haniffa$mutation2) %>% table() %>% data.frame() %>%
    rename(count = Freq) %>%
    mutate(celltype = factor(celltype, levels=names(colors.mn))))
celltype.table2 <- data.frame(celltype = MULT$predicted.scRNA.monaco, status = MULT$mutation2) %>% table() %>% data.frame() %>%
    rename(count = Freq) %>%
    mutate(celltype = factor(celltype, levels=names(colors.mn))) %>%
    bind_rows(data.frame(celltype = haniffa$monaco, status = haniffa$mutation2) %>% table() %>% data.frame() %>%
    rename(count = Freq) %>%
    mutate(celltype = factor(celltype, levels=names(colors.mn))))
celltype.sum <- celltype.table1 %>% group_by(status) %>% summarise(count = sum(count))
p1 <- ggplot(celltype.table1) + 
    geom_bar(aes(fill=celltype, x=count, y=factor(status, levels = rev(levels.status))), position=position_fill(reverse = TRUE), stat="identity", color = "black") +
    geom_text(data = celltype.sum, mapping = aes(x = 1.1, y = factor(status, levels = rev(levels.status)), label = count), hjust = 1) +
    geom_segment(aes(y=-Inf,yend=-Inf,x=0,xend=1)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_x_continuous(expand = expansion(mult = c(0, .08)), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_fill_manual(values = colors.mn) +
    xlab("Cell type proportion") +
    theme_cowplot() + 
    theme(legend.position = "none", legend.title=element_blank(), axis.title.y = element_blank(), axis.line.x=element_blank())
# , axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
celltype.sum <- celltype.table2 %>% group_by(status) %>% summarise(count = sum(count))
p2 <- ggplot(celltype.table2) + 
    geom_bar(aes(fill=celltype, x=count, y=factor(status, levels = rev(levels.status))), position=position_fill(reverse = TRUE), stat="identity", color = "black") +
    geom_text(data = celltype.sum, mapping = aes(x = 1.1, y = factor(status, levels = rev(levels.status)), label = count), hjust = 1) +
    geom_segment(aes(y=-Inf,yend=-Inf,x=0,xend=1)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_x_continuous(expand = expansion(mult = c(0, .08)), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_fill_manual(values = colors.mn) +
    xlab("Cell type proportion") +
    theme_cowplot() + 
    theme(legend.position = "none", legend.title=element_blank(), axis.title.y = element_blank(), axis.line.x=element_blank())

pdf(file='celltype-proportions-final.pdf', width=8, height=3)
p1
p2
dev.off()

celltype.table1 <- data.frame(celltype = GEX$monaco, patient = GEX$sample) %>% table() %>% data.frame() %>%
    rename(count = Freq) %>%
    pivot_wider(names_from = "celltype", values_from = "count") %>%
    left_join(data.frame(patient = GEX$sample, status = GEX$mutation3) %>% distinct()) %>% 
    relocate(status, .after = "patient") %>%
    arrange(status)

write.table(celltype.table1, file = "celltypes.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



celltype.table1 <- data.frame(celltype = GEX$monaco, sample = GEX$sample) %>% table() %>% data.frame() %>%
    rename(count = Freq) %>%
    mutate(celltype = factor(celltype, levels=names(colors.mn))) %>%
    group_by(sample) %>%
    left_join(GEX@meta.data %>% select(sample, mutation2) %>% distinct()) %>%
    filter(mutation2 %in% c("COVID (no CHIP)", "Non-TET2 CHIP + COVID", "TET2 CHIP + COVID")) %>%
    mutate(mutation2 = factor(mutation2, levels = c("TET2 CHIP + COVID", "Non-TET2 CHIP + COVID", "COVID (no CHIP)")))
celltype.table3 <- celltype.table1 %>% group_by(celltype, mutation2) %>% summarize(count = mean(count))
celltype.table3 <- celltype.table3 %>% left_join(celltype.table3 %>% group_by(mutation2) %>% summarize(total = sum(count))) %>% mutate(freq = count/total) %>% filter(celltype %in% c("CD14 Mono", "CD16 Mono", "Int Mono")) %>% ungroup()
celltype.table3 <- celltype.table3 %>% left_join(celltype.table3 %>% group_by(mutation2) %>% summarize(mo_total = sum(count))) %>% mutate(mo_freq = count/mo_total) %>% ungroup() %>% mutate(mutation2 = factor(mutation2, levels = rev(c("COVID (no CHIP)", "Non-TET2 CHIP + COVID", "TET2 CHIP + COVID"))))

celltype.table1 <- data.frame(celltype = MULT$predicted.scRNA.monaco, sample = MULT$sample) %>% table() %>% data.frame() %>%
    rename(count = Freq) %>%
    mutate(celltype = factor(celltype, levels=names(colors.mn))) %>%
    group_by(sample) %>%
    left_join(MULT@meta.data %>% select(sample, mutation2) %>% distinct()) %>%
    filter(mutation2 %in% c("COVID (no CHIP)", "Non-TET2 CHIP + COVID", "TET2 CHIP + COVID"))
celltype.table4 <- celltype.table1 %>% group_by(celltype, mutation2) %>% summarize(count = mean(count))
celltype.table4 <- celltype.table4 %>% left_join(celltype.table4 %>% group_by(mutation2) %>% summarize(total = sum(count))) %>% mutate(freq = count/total) %>% filter(celltype %in% c("CD14 Mono", "CD16 Mono", "Int Mono")) %>% ungroup()
celltype.table4 <- celltype.table4 %>% left_join(celltype.table4 %>% group_by(mutation2) %>% summarize(mo_total = sum(count))) %>% mutate(mo_freq = count/mo_total) %>% ungroup() %>% mutate(mutation2 = factor(mutation2, levels = rev(c("COVID (no CHIP)", "Non-TET2 CHIP + COVID", "TET2 CHIP + COVID"))))

p3 <- ggplot(celltype.table3, aes(y = mutation2, x = freq, fill = celltype)) +
    geom_col(position = position_stack(), color = "black") + 
    scale_y_discrete(expand = c(0,0)) +
    geom_text(aes(label = sprintf("%0.3f", round(freq, digits = 3))), position = position_stack(vjust = 0.5), angle = 270) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 0.25), position = "top") +
    scale_fill_manual(values = colors.mn[c("CD14 Mono", "CD16 Mono", "Int Mono")]) + 
    xlab("Proportion of monocyte sub types in all cells") +
    ggtitle("scRNA-seq")
p4 <- ggplot(celltype.table4, aes(y = mutation2, x = freq, fill = celltype)) +
    geom_col(position = position_stack(), color = "black") + 
    scale_y_discrete(expand = c(0,0)) +
    geom_text(aes(label = sprintf("%0.3f", round(freq, digits = 3))), position = position_stack(vjust = 0.5), angle = 270) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 0.25), position = "top") +
    scale_fill_manual(values = colors.mn[c("CD14 Mono", "CD16 Mono", "Int Mono")]) + 
    xlab("Proportion of monocyte sub types in all cells") +
    ggtitle("scMultiome")
p5 <- ggplot(celltype.table3, aes(y = mutation2, x = count, fill = celltype)) +
    geom_bar(position = position_fill(), color = "black", stat="identity") + 
    scale_y_discrete(expand = c(0,0)) +
    geom_text(aes(label = sprintf("%0.3f", round(mo_freq, digits = 3))), position = position_fill(vjust = 0.5), angle = 270) +
    scale_x_continuous(expand = c(0,0), position = "top") +
    scale_fill_manual(values = colors.mn[c("CD14 Mono", "CD16 Mono", "Int Mono")]) + 
    xlab("Proportion of monocyte sub types in monocytes") +
    ggtitle("scRNA-seq")
p6 <- ggplot(celltype.table4, aes(y = mutation2, x = count, fill = celltype)) +
    geom_bar(position = position_fill(), color = "black", stat="identity") + 
    scale_y_discrete(expand = c(0,0)) +
    geom_text(aes(label = sprintf("%0.3f", round(mo_freq, digits = 3))), position = position_fill(vjust = 0.5), angle = 270) +
    scale_x_continuous(expand = c(0,0), position = "top") +
    scale_fill_manual(values = colors.mn[c("CD14 Mono", "CD16 Mono", "Int Mono")]) + 
    xlab("Proportion of monocyte sub types in monocytes") +
    ggtitle("scMultiome")


pdf(file='monocyte-proportions-final.pdf', width=10, height=10)
(p3 / p4 / p5 / p6) & theme_cowplot() + theme(legend.position="none", axis.title.y = element_blank())
dev.off()

########################  DGEA  #######################
# degs_auto_*.R
degsr <- read.table(file = "/FINAL/scRNA/degsr_TET2C_COVID.tsv", sep = "\t", skip = 1) %>% set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"))

# pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
# goi <- ConvertMotifID(MULT, id = names(pwm_set))
goi <- c("SP1", "EGR1", "NEAT1", "EGR3", "KLF4", "TET2", "MALAT1", "S100A8", "S100A9", "ELF2", rownames(GEX)[str_starts(rownames(GEX), "CEBP")])

temp1 <- degsr %>% 
    filter(p_val_adj < 0.05) %>% 
    mutate(celltype = factor(celltype, levels = names(colors.mn))) %>% 
    mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_adj_pval = -log10(p_val_adj))
pdf(file = "degsr_TET2-C_COVID.pdf", width = 8, height = 8)
p1 <- ggplot(temp1, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
    geom_point(aes(color = celltype), size = 1.5, stroke=0.1, shape = 16) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    annotate("text", label = "Adj. P value = 0.05", x = 3.5, y = -log10(0.05)+5, size = 3) +
    geom_label_repel(data = temp1 %>% mutate(gene = case_when(gene %in% goi ~ gene, TRUE ~ "")), aes(label=gene, segment.color=celltype), arrow = arrow(length = unit(0.010, "npc")), size = 3, max.overlaps = Inf, min.segment.length = 0, box.padding = 1, point.padding = 0.1, show.legend = FALSE, fontface = "italic") +
    coord_cartesian(xlim = c(-4, 4)) +
    scale_color_manual(values = colors.mn, aesthetics = c("color", "segment.color")) +
    labs(x = "log2(Fold change of average expression)", y = "-log10(Adjusted P value)") +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.title=element_blank())
p2 <- table((temp1 %>% filter(avg_log2FC < 0))$celltype) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn))), direction = "DOWN") %>%
    bind_rows(table((temp1 %>% filter(avg_log2FC > 0))$celltype) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn))), direction = "UP")) %>%
    mutate(direction = factor(direction, levels = c("UP", "DOWN"))) %>% 
    ggplot(aes(fill=Var1, y=Freq, x=direction)) +
    geom_bar(stat="identity", width = 0.6) + 
    scale_fill_manual(values = colors.mn) + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip() +
    theme_cowplot() +
    theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), plot.margin = margin(0.5, 0.1, 0.5, 0.1, "cm"))
plot(p1 + p2 + plot_layout(heights = c(15, 1)) + plot_annotation(theme = theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))))
dev.off()

degs <- read.table(file = "/FINAL/scRNA/degs_TET2C_COVID.tsv", sep = "\t", skip = 1) %>% set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% filter(p_val_adj < 0.05)
degsr <- read.table(file = "/FINAL/scRNA/degsr_TET2C_COVID.tsv", sep = "\t", skip = 1) %>% set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% filter(p_val_adj < 0.05)

write.table(degs %>% select(celltype, gene, avg_log2FC, p_val_adj), file = "degs_forIPA.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(degsr %>% select(celltype, gene, avg_log2FC, p_val_adj), file = "degsr_forIPA.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(degsr, file = "degsr_forSuppTable.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

degsr1 <- read.table(file = "/FINAL/scRNA/degsr_TET2C_COVID.tsv", sep = "\t", skip = 1) %>% 
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% 
    filter(p_val_adj < 0.05) %>% 
    mutate(celltype = factor(celltype, levels = names(colors.mn))) %>% 
    mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_adj_pval = -log10(p_val_adj))
degsr2 <- read.table(file = "/FINAL/scRNA/degsr_DNMT3A-C_COVID.tsv", sep = "\t", skip = 1) %>% 
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% 
    filter(p_val_adj < 0.05) %>% 
    mutate(celltype = factor(celltype, levels = names(colors.mn))) %>% 
    mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_adj_pval = -log10(p_val_adj))
degsr3 <- read.table(file = "/FINAL/scRNA/degsr_ASXL1-C_COVID.tsv", sep = "\t", skip = 1) %>% 
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% 
    filter(p_val_adj < 0.05) %>% 
    mutate(celltype = factor(celltype, levels = names(colors.mn))) %>% 
    mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_adj_pval = -log10(p_val_adj))
degsr4 <- read.table(file = "/FINAL/scRNA/degsr_EZH2-C_COVID.tsv", sep = "\t", skip = 1) %>% 
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% 
    filter(p_val_adj < 0.05) %>% 
    mutate(celltype = factor(celltype, levels = names(colors.mn))) %>% 
    mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_adj_pval = -log10(p_val_adj))
degsr5 <- read.table(file = "/FINAL/scRNA/degsr_COVIDCHIP_COVID.tsv", sep = "\t", skip = 1) %>% 
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% 
    filter(p_val_adj < 0.05) %>% 
    mutate(celltype = factor(celltype, levels = names(colors.mn))) %>% 
    mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_adj_pval = -log10(p_val_adj))

p1 <- ggplot(degsr1, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
    geom_point(aes(color = celltype), size = 1.5, stroke=0.1, shape = 16) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    annotate("text", label = "P value = 0.05", x = 2, y = -log10(0.05)+5, size = 4) +
    geom_label_repel(data = degsr1 %>% mutate(gene = case_when(gene %in% goi ~ gene, TRUE ~ "")), aes(label=gene, segment.color=celltype), arrow = arrow(length = unit(0.010, "npc")), size = 3, max.overlaps = Inf, min.segment.length = 0, box.padding = 1, point.padding = 0.1, show.legend = FALSE, fontface = "italic") +
    coord_cartesian(xlim = c(-3, 3)) +
    scale_color_manual(values = colors.mn, aesthetics = c("color", "segment.color")) +
    labs(x = "log2(fold change of average expression)", y = "-log10(adjusted p value)") +
    guides(color = guide_legend(override.aes = list(size=5))) + ggtitle("COVID vs COVID+TET2") +
    theme_cowplot() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0, "cm"))
p2 <- table((degsr1 %>% filter(avg_log2FC < 0))$celltype) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn))), direction = "DOWN") %>%
    bind_rows(table((degsr1 %>% filter(avg_log2FC > 0))$celltype) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn))), direction = "UP")) %>%
    mutate(direction = factor(direction, levels = c("UP", "DOWN"))) %>% 
    ggplot(aes(fill=Var1, y=Freq, x=direction)) +
    geom_bar(stat="identity", width = 0.6) + 
    scale_fill_manual(values = colors.mn) + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip() +
    theme_cowplot() +
    theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
p3 <- ggplot(degsr2, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
    geom_point(aes(color = celltype), size = 1.5, stroke=0.1, shape = 16) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    annotate("text", label = "P value = 0.05", x = 2, y = -log10(0.05)+5, size = 4) +
    geom_label_repel(data = degsr2 %>% mutate(gene = case_when(gene %in% goi ~ gene, TRUE ~ "")), aes(label=gene, segment.color=celltype), arrow = arrow(length = unit(0.010, "npc")), size = 3, max.overlaps = Inf, min.segment.length = 0, box.padding = 1, point.padding = 0.1, show.legend = FALSE, fontface = "italic") +
    coord_cartesian(xlim = c(-3, 3)) +
    scale_color_manual(values = colors.mn, aesthetics = c("color", "segment.color")) +
    labs(x = "log2(fold change of average expression)", y = "-log10(adjusted p value)") +
    guides(color = guide_legend(override.aes = list(size=5))) + ggtitle("COVID vs COVID+DNMT3A") +
    theme_cowplot() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0, "cm"))
p4 <- table((degsr2 %>% filter(avg_log2FC < 0))$celltype) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn))), direction = "DOWN") %>%
    bind_rows(table((degsr2 %>% filter(avg_log2FC > 0))$celltype) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn))), direction = "UP")) %>%
    mutate(direction = factor(direction, levels = c("UP", "DOWN"))) %>% 
    ggplot(aes(fill=Var1, y=Freq, x=direction)) +
    geom_bar(stat="identity", width = 0.6) + 
    scale_fill_manual(values = colors.mn) + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip() +
    theme_cowplot() +
    theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
p5 <- ggplot(degsr3, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
    geom_point(aes(color = celltype), size = 1.5, stroke=0.1, shape = 16) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    annotate("text", label = "P value = 0.05", x = 2, y = -log10(0.05)+5, size = 4) +
    geom_label_repel(data = degsr3 %>% mutate(gene = case_when(gene %in% goi ~ gene, TRUE ~ "")), aes(label=gene, segment.color=celltype), arrow = arrow(length = unit(0.010, "npc")), size = 3, max.overlaps = Inf, min.segment.length = 0, box.padding = 1, point.padding = 0.1, show.legend = FALSE, fontface = "italic") +
    coord_cartesian(xlim = c(-3, 3)) +
    scale_color_manual(values = colors.mn, aesthetics = c("color", "segment.color")) +
    labs(x = "log2(fold change of average expression)", y = "-log10(adjusted p value)") +
    guides(color = guide_legend(override.aes = list(size=5))) + ggtitle("COVID vs COVID+ASXL1") +
    theme_cowplot() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0, "cm"))
p6 <- table((degsr3 %>% filter(avg_log2FC < 0))$celltype) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn))), direction = "DOWN") %>%
    bind_rows(table((degsr3 %>% filter(avg_log2FC > 0))$celltype) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn))), direction = "UP")) %>%
    mutate(direction = factor(direction, levels = c("UP", "DOWN"))) %>% 
    ggplot(aes(fill=Var1, y=Freq, x=direction)) +
    geom_bar(stat="identity", width = 0.6) + 
    scale_fill_manual(values = colors.mn) + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip() +
    theme_cowplot() +
    theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
p7 <- ggplot(degsr4, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
    geom_point(aes(color = celltype), size = 1.5, stroke=0.1, shape = 16) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    annotate("text", label = "P value = 0.05", x = 2, y = -log10(0.05)+5, size = 4) +
    geom_label_repel(data = degsr4 %>% mutate(gene = case_when(gene %in% goi ~ gene, TRUE ~ "")), aes(label=gene, segment.color=celltype), arrow = arrow(length = unit(0.010, "npc")), size = 3, max.overlaps = Inf, min.segment.length = 0, box.padding = 1, point.padding = 0.1, show.legend = FALSE, fontface = "italic") +
    coord_cartesian(xlim = c(-3, 3)) +
    scale_color_manual(values = colors.mn, aesthetics = c("color", "segment.color")) +
    labs(x = "log2(fold change of average expression)", y = "-log10(adjusted p value)") +
    guides(color = guide_legend(override.aes = list(size=5))) + ggtitle("COVID vs COVID+EZH2") +
    theme_cowplot() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0, "cm"))
p8 <- table((degsr4 %>% filter(avg_log2FC < 0))$celltype) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn))), direction = "DOWN") %>%
    bind_rows(table((degsr4 %>% filter(avg_log2FC > 0))$celltype) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn))), direction = "UP")) %>%
    mutate(direction = factor(direction, levels = c("UP", "DOWN"))) %>% 
    ggplot(aes(fill=Var1, y=Freq, x=direction)) +
    geom_bar(stat="identity", width = 0.6) + 
    scale_fill_manual(values = colors.mn) + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip() +
    theme_cowplot() +
    theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
p9 <- ggplot(degsr5, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
    geom_point(aes(color = celltype), size = 1.5, stroke=0.1, shape = 16) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    annotate("text", label = "P value = 0.05", x = 2, y = -log10(0.05)+5, size = 4) +
    geom_label_repel(data = degsr5 %>% mutate(gene = case_when(gene %in% goi ~ gene, TRUE ~ "")), aes(label=gene, segment.color=celltype), arrow = arrow(length = unit(0.010, "npc")), size = 3, max.overlaps = Inf, min.segment.length = 0, box.padding = 1, point.padding = 0.1, show.legend = FALSE, fontface = "italic") +
    coord_cartesian(xlim = c(-3, 3)) +
    scale_color_manual(values = colors.mn, aesthetics = c("color", "segment.color")) +
    labs(x = "log2(fold change of average expression)", y = "-log10(adjusted p value)") +
    guides(color = guide_legend(override.aes = list(size=5))) + ggtitle("COVID vs COVID+CHIP") +
    theme_cowplot() +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))
p10 <- table((degsr5 %>% filter(avg_log2FC < 0))$celltype) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn))), direction = "DOWN") %>%
    bind_rows(table((degsr5 %>% filter(avg_log2FC > 0))$celltype) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn))), direction = "UP")) %>%
    mutate(direction = factor(direction, levels = c("UP", "DOWN"))) %>% 
    ggplot(aes(fill=Var1, y=Freq, x=direction)) +
    geom_bar(stat="identity", width = 0.6) + 
    scale_fill_manual(values = colors.mn) + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip() +
    theme_cowplot() +
    theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())

layout <- '
ACEGI
BDFHJ
'
pdf(file = "degsr_ALL.pdf", width = 30, height = 8)
plot(p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + plot_layout(heights = c(15, 1, 15, 1, 15, 1, 15, 1, 15, 1), design = layout))
dev.off()



######################## Markers ##############################

Idents(GEX) <- "monaco"
all.markers.seurat <- FindAllMarkers(object = GEX)
saveRDS(all.markers.seurat, "all.markers.seurat.rds")

all.markers.seurat <- readRDS("all.markers.seurat.rds")
all.markers.seurat <- all.markers.seurat %>% filter(!(gene %in% str_subset(all.markers.seurat$gene, "^RP") | gene %in% str_subset(all.markers.seurat$gene, "^MT")))

GEX_CVCVCHIP <- subset(GEX, subset = mutation1 %in% c("COVID", "TET2-C"))
# idents_sub <- c("CD4 T", "CD8 T", "MAIT", "NK", "B", "CD14 Mono", "CD16 Mono", "Int Mono", "Neutrophils")
# top.markers.sub <- all.markers.seurat %>% filter(cluster %in% idents_sub) %>% group_by(cluster) %>% slice_head(n = 10)
# top.markers.sub.unique <- top.markers.sub %>% group_by(gene) %>% mutate(n = n()) %>% filter(n == 1)
# top.markers.sub.multi <- top.markers.sub %>% group_by(gene) %>% mutate(n = n()) %>% filter(n > 1)
top.markers.sub <- all.markers.seurat %>% group_by(cluster) %>% slice_head(n = 10)

# source("tempDotPlot.R")
library(scico)
pdf(file='markers.pdf', width=35, height=6)
GEX$monaco <- factor(GEX$monaco, levels = rev(names(colors.mn)))
Idents(GEX) <- "monaco"
plotmarkers <- c(unique(top.markers.sub$gene))
p1 <- DotPlot(GEX, assay = "RNA", features = plotmarkers, cluster.idents = FALSE) +
    scale_color_scico(palette = 'tokyo', direction = -1, limits = c(-2.5, 2.5)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face="italic", size=20), axis.title = element_blank(), axis.text.y = element_text(size=20))
p1
dev.off()

GEX_COVID <- subset(GEX, subset = mutation1 %in% c("COVID"))
GEX_TET2C <- subset(GEX, subset = mutation1 %in% c("TET2-C"))
p1 <- DotPlot(GEX_COVID, assay = "RNA", features = c("TET2", "DNMT3A", "EZH2", "ASXL1", "SF3B1"), scale.max = 60, scale.min = 0) +
    scale_color_scico(palette = 'vik', direction = -1, limits = c(-2.5, 2.5), begin = 0, end = 0.5) +
    scale_radius(breaks = c(0, 20, 40, 60)) +
    guides(color = guide_colorbar(order = 1, title = 'Scaled AvgExpr'), radius = guide_legend(order = 2)) +
    ggtitle("COVID") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face="italic"), axis.title = element_blank())
p2 <- DotPlot(GEX_TET2C, assay = "RNA", features = c("TET2", "DNMT3A", "EZH2", "ASXL1", "SF3B1"), scale.max = 60, scale.min = 0) +
    scale_color_scico(palette = 'vik', direction = 1, limits = c(-2.5, 2.5), begin = 0.5, end = 1) +
    scale_radius(breaks = c(0, 20, 40, 60)) +
    guides(color = guide_colorbar(order = 1, title = 'Scaled AvgExpr'), radius = guide_legend(order = 2)) +
    ggtitle("TET2+COVID") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face="italic"), axis.title = element_blank())
pdf(file='CHIP_Genes.pdf', width=12, height=5)
p1 | p2
dev.off()

GEX_CVCVCHIP <- subset(GEX, subset = mutation1 %in% c("COVID", "TET2-C"))
genes <- c("CD14", "FCGR3A", "CEACAM8", "LTF", "NCAM1", "HLA-DRA", "CD4", "CCR2", "FCGR1A", "CD33", "ANPEP", "CSF1R", "CSF1")
pdf(file='mymarkers.pdf', width=10, height=10)
GEX_CVCVCHIP$monaco <- factor(GEX_CVCVCHIP$monaco, levels = rev(names(colors.mn)))
Idents(GEX_CVCVCHIP) <- "monaco"
DotPlot(GEX_CVCVCHIP, assay = "RNA", features = genes, cluster.idents = FALSE, cols = c("blue", "red"), split.by = "mutation1") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_blank())
dev.off()


####################### Gene Expression #######################

MAYOPLUSHANIFFA <- readRDS("/FINAL/scRNA/MAYOPLUSHANIFFA.rds")
MAYOPLUSHANIFFA$celltype.status <- paste(MAYOPLUSHANIFFA$monaco, MAYOPLUSHANIFFA$group, sep = ".")
MAYOPLUSHANIFFA_avg <- AverageExpression(MAYOPLUSHANIFFA, assays = c("RNA"), group.by = "celltype.status", return.seurat = TRUE)
genes <- c("TET2", "MALAT1", "NEAT1", "EGR1", "S100A8", "S100A9")
# genes <- c("TNF", "LTA", "IL6", "IL10", "CCL2", "IL1B", "CSF2")
celltypes <- c("CD14 Mono", "Int Mono", "CD16 Mono", "cDC")

Idents(MAYOPLUSHANIFFA) <- "celltype.status"
degs <- lapply(
    celltypes, 
    FUN = function(celltype){
        temp1 <- FindMarkers(MAYOPLUSHANIFFA, features = genes, min.cells.group = 0, min.cells.feature = 0, logfc.threshold = 0, ident.1 = paste0(celltype, ".TET2C"), ident.2 = paste0(celltype, ".COVID"), min.pct = 0) %>% rownames_to_column("gene")
    }
)
names(degs) <- celltypes
degs <- bind_rows(degs, .id = "celltype")
brackets <- degs %>% mutate(sig = case_when(p_val_adj <= 0.0001 ~ "****", p_val_adj <= 0.001 ~ "***", p_val_adj <= 0.01 ~ "**", p_val_adj <= 0.05 ~ "*", TRUE ~ "ns"))

plotdata1 <- data.frame(t(MAYOPLUSHANIFFA[["RNA"]]@data[genes,]), check.names = FALSE) %>% 
    mutate(status = MAYOPLUSHANIFFA$group, celltype = factor(MAYOPLUSHANIFFA$monaco, levels = names(colors.mn)), sample = MAYOPLUSHANIFFA$sample) %>%
    filter(status %in% c("CTRLY", "CTRLO", "TET2", "COVID", "TET2C")) %>%
    filter(celltype %in% c("CD14 Mono", "Int Mono", "CD16 Mono", "cDC")) %>%
    # filter(celltype %in% c("CD4 T", "CD8 T", "NK", "B", "CD14 Mono", "Int Mono", "CD16 Mono", "cDC")) %>%
    pivot_longer(!c(status, celltype, sample), names_to = "gene", values_to = "expr") %>%
    mutate(status = factor(recode_factor(status, !!!c(`CTRLY` = "Healthy(<50)", `CTRLO` = "Healthy(>50)", `TET2` = "COVID-/TET2MT", `COVID` = "COVID+/CHIP-", `TET2C` = "COVID+/TET2MT")), levels = c("Healthy(<50)", "Healthy(>50)", "COVID-/TET2MT", "COVID+/CHIP-", "COVID+/TET2MT")))
mu1 <- data.frame(MAYOPLUSHANIFFA_avg[["RNA"]]@data, check.names = FALSE) %>%
    rownames_to_column("gene") %>% 
    pivot_longer(!gene, names_to = c("celltype", "status"), values_to = "expr", names_sep = "\\.") %>%
    pivot_wider(names_from = "gene", values_from = "expr") %>%
    select(c("celltype", "status", all_of(genes))) %>%
    filter(status %in% c("CTRLY", "CTRLO", "TET2", "COVID", "TET2C")) %>%
    filter(celltype %in% c("CD14 Mono", "Int Mono", "CD16 Mono", "cDC")) %>%
    # filter(celltype %in% c("CD4 T", "CD8 T", "NK", "B", "CD14 Mono", "Int Mono", "CD16 Mono", "cDC")) %>%
    pivot_longer(!c(status, celltype), names_to = "gene", values_to = "expr") %>%
    # mutate(celltype = factor(celltype, levels = c("CD14 Mono", "Int Mono", "CD16 Mono", "cDC"))) %>%
    mutate(celltype = factor(celltype, levels = c("CD4 T", "CD8 T", "NK", "B", "CD14 Mono", "Int Mono", "CD16 Mono", "cDC"))) %>%
    mutate(status = factor(recode_factor(status, !!!c(`CTRLY` = "Healthy(<50)", `CTRLO` = "Healthy(>50)", `TET2` = "COVID-/TET2MT", `COVID` = "COVID+/CHIP-", `TET2C` = "COVID+/TET2MT")), levels = c("Healthy(<50)", "Healthy(>50)", "COVID-/TET2MT", "COVID+/CHIP-", "COVID+/TET2MT")))

generate_expr_plots <- function(curgene){return(
ggplot(data = plotdata1 %>% filter(gene == curgene), aes(x = status, y = expr)) +
    geom_violin(aes(fill = celltype)) +
    # ggrastr::rasterise(geom_jitter(aes(color = sample), size=0.2, stroke=0.1, shape=16), dpi = 400) +
    geom_point(data = mu1 %>% filter(gene == curgene), mapping = aes(x = status, y = expr), size=1.5, color="black") +
    scale_fill_manual(values = colors.mn) +
    geom_bracket(data = brackets %>% filter(gene == curgene), mapping = aes(label = sig), xmin = "COVID+/CHIP-", xmax = "COVID+/TET2MT", y.position = plotdata1 %>% filter(gene == curgene) %>% select(expr) %>% max()) +
    # stat_compare_means(aes(group=celltype, label = ..p.signif..), comparisons = list(c("COVID+/CHIP-", "COVID+/TET2MT"))) +
    scale_y_continuous(expand = expansion(mult = c(0.01, .09))) +
    facet_grid(cols = vars(celltype), switch="both") +
    ggtitle(curgene) +
    ylab("LogNormalized Expression") +
    theme_cowplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(), plot.title = element_text(face = "italic"))
)}

pdf(file='GeneExpr.pdf', width=20, height=12)
wrap_plots(lapply(genes, generate_expr_plots))
dev.off()

### MULT ##

MULTPLUSHANIFFA <- readRDS("/FINAL/scMULT/MULTPLUSHANIFFA.rds")
MULTPLUSHANIFFA$celltype.status <- paste(MULTPLUSHANIFFA$monaco, MULTPLUSHANIFFA$group, sep = ".")
MULTPLUSHANIFFA_avg <- AverageExpression(MULTPLUSHANIFFA, assays = c("RNA"), group.by = "celltype.status", return.seurat = TRUE)
genes <- c("TET2", "MALAT1", "NEAT1", "EGR1", "S100A8", "S100A9")
# genes <- c("TNF", "LTA", "IL6", "IL10", "CCL2", "IL1B", "CSF2")
celltypes <- c("CD14 Mono", "Int Mono", "CD16 Mono", "cDC")

Idents(MULTPLUSHANIFFA) <- "celltype.status"
degs <- lapply(
    celltypes, 
    FUN = function(celltype){
        temp1 <- FindMarkers(MULTPLUSHANIFFA, features = genes, min.cells.group = 0, min.cells.feature = 0, logfc.threshold = 0, ident.1 = paste0(celltype, ".TET2C"), ident.2 = paste0(celltype, ".COVID"), min.pct = 0) %>% rownames_to_column("gene")
    }
)
names(degs) <- celltypes
degs <- bind_rows(degs, .id = "celltype")
brackets <- degs %>% mutate(sig = case_when(p_val_adj <= 0.0001 ~ "****", p_val_adj <= 0.001 ~ "***", p_val_adj <= 0.01 ~ "**", p_val_adj <= 0.05 ~ "*", TRUE ~ "ns"))

plotdata1 <- data.frame(t(MULTPLUSHANIFFA[["RNA"]]@data[genes,]), check.names = FALSE) %>% 
    mutate(status = MULTPLUSHANIFFA$group, celltype = factor(MULTPLUSHANIFFA$monaco, levels = names(colors.mn)), sample = MULTPLUSHANIFFA$sample) %>%
    filter(status %in% c("CTRLY", "CTRLO", "TET2", "COVID", "TET2C")) %>%
    filter(celltype %in% c("CD14 Mono", "Int Mono", "CD16 Mono", "cDC")) %>%
    # filter(celltype %in% c("CD4 T", "CD8 T", "NK", "B", "CD14 Mono", "Int Mono", "CD16 Mono", "cDC")) %>%
    pivot_longer(!c(status, celltype, sample), names_to = "gene", values_to = "expr") %>%
    mutate(status = factor(recode_factor(status, !!!c(`CTRLY` = "Healthy(<50)", `CTRLO` = "Healthy(>50)", `TET2` = "COVID-/TET2MT", `COVID` = "COVID+/CHIP-", `TET2C` = "COVID+/TET2MT")), levels = c("Healthy(<50)", "Healthy(>50)", "COVID-/TET2MT", "COVID+/CHIP-", "COVID+/TET2MT")))
mu1 <- data.frame(MULTPLUSHANIFFA_avg[["RNA"]]@data, check.names = FALSE) %>%
    rownames_to_column("gene") %>% 
    pivot_longer(!gene, names_to = c("celltype", "status"), values_to = "expr", names_sep = "\\.") %>%
    pivot_wider(names_from = "gene", values_from = "expr") %>%
    select(c("celltype", "status", all_of(genes))) %>%
    filter(status %in% c("CTRLY", "CTRLO", "TET2", "COVID", "TET2C")) %>%
    filter(celltype %in% c("CD14 Mono", "Int Mono", "CD16 Mono", "cDC")) %>%
    # filter(celltype %in% c("CD4 T", "CD8 T", "NK", "B", "CD14 Mono", "Int Mono", "CD16 Mono", "cDC")) %>%
    pivot_longer(!c(status, celltype), names_to = "gene", values_to = "expr") %>%
    # mutate(celltype = factor(celltype, levels = c("CD14 Mono", "Int Mono", "CD16 Mono", "cDC"))) %>%
    mutate(celltype = factor(celltype, levels = c("CD4 T", "CD8 T", "NK", "B", "CD14 Mono", "Int Mono", "CD16 Mono", "cDC"))) %>%
    mutate(status = factor(recode_factor(status, !!!c(`CTRLY` = "Healthy(<50)", `CTRLO` = "Healthy(>50)", `TET2` = "COVID-/TET2MT", `COVID` = "COVID+/CHIP-", `TET2C` = "COVID+/TET2MT")), levels = c("Healthy(<50)", "Healthy(>50)", "COVID-/TET2MT", "COVID+/CHIP-", "COVID+/TET2MT")))

generate_expr_plots <- function(curgene){return(
ggplot(data = plotdata1 %>% filter(gene == curgene), aes(x = status, y = expr)) +
    geom_violin(aes(fill = celltype)) +
    # ggrastr::rasterise(geom_jitter(aes(color = sample), size=0.2, stroke=0.1, shape=16), dpi = 400) +
    geom_point(data = mu1 %>% filter(gene == curgene), mapping = aes(x = status, y = expr), size=1.5, color="black") +
    scale_fill_manual(values = colors.mn) +
    geom_bracket(data = brackets %>% filter(gene == curgene), mapping = aes(label = sig), xmin = "COVID+/CHIP-", xmax = "COVID+/TET2MT", y.position = plotdata1 %>% filter(gene == curgene) %>% select(expr) %>% max()) +
    # stat_compare_means(aes(group=celltype, label = ..p.signif..), comparisons = list(c("COVID+/CHIP-", "COVID+/TET2MT"))) +
    scale_y_continuous(expand = expansion(mult = c(0.01, .09))) +
    facet_grid(cols = vars(celltype), switch="both") +
    ggtitle(curgene) +
    ylab("LogNormalized Expression") +
    theme_cowplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(), plot.title = element_text(face = "italic"))
)}

pdf(file='GeneExpr_MULT.pdf', width=20, height=12)
wrap_plots(lapply(genes, generate_expr_plots))
dev.off()


###########

genes <- c("S100A8", "S100A9")
celltypes <- c("B", "cDC")

Idents(MAYOPLUSHANIFFA) <- "celltype.status"
degs <- lapply(
    celltypes, 
    FUN = function(celltype){
        temp1 <- FindMarkers(MAYOPLUSHANIFFA, features = genes, min.cells.group = 0, min.cells.feature = 0, logfc.threshold = 0, ident.1 = paste0(celltype, ".TET2C"), ident.2 = paste0(celltype, ".COVID"), min.pct = 0) %>% rownames_to_column("gene")
    }
)
names(degs) <- celltypes
degs <- bind_rows(degs, .id = "celltype")
brackets <- degs %>% mutate(sig = case_when(p_val_adj <= 0.0001 ~ "****", p_val_adj <= 0.001 ~ "***", p_val_adj <= 0.01 ~ "**", p_val_adj <= 0.05 ~ "*", TRUE ~ "ns"))

plotdata1 <- data.frame(t(MAYOPLUSHANIFFA[["RNA"]]@data[genes,]), check.names = FALSE) %>% 
    mutate(status = MAYOPLUSHANIFFA$group, celltype = factor(MAYOPLUSHANIFFA$monaco, levels = names(colors.mn)), sample = MAYOPLUSHANIFFA$sample) %>%
    filter(status %in% c("CTRLY", "CTRLO", "TET2", "COVID", "TET2C")) %>%
    filter(celltype %in% c("B", "cDC")) %>%
    pivot_longer(!c(status, celltype, sample), names_to = "gene", values_to = "expr") %>%
    mutate(status = factor(recode_factor(status, !!!c(`CTRLY` = "Healthy(<50)", `CTRLO` = "Healthy(>50)", `TET2` = "COVID-/TET2MT", `COVID` = "COVID+/CHIP-", `TET2C` = "COVID+/TET2MT")), levels = c("Healthy(<50)", "Healthy(>50)", "COVID-/TET2MT", "COVID+/CHIP-", "COVID+/TET2MT")))
mu1 <- data.frame(MAYOPLUSHANIFFA_avg[["RNA"]]@data, check.names = FALSE) %>%
    rownames_to_column("gene") %>% 
    pivot_longer(!gene, names_to = c("celltype", "status"), values_to = "expr", names_sep = "\\.") %>%
    pivot_wider(names_from = "gene", values_from = "expr") %>%
    select(c("celltype", "status", all_of(genes))) %>%
    filter(status %in% c("CTRLY", "CTRLO", "TET2", "COVID", "TET2C")) %>%
    filter(celltype %in% c("B", "cDC")) %>%
    pivot_longer(!c(status, celltype), names_to = "gene", values_to = "expr") %>%
    mutate(celltype = factor(celltype, levels = c("B", "cDC"))) %>%
    mutate(status = factor(recode_factor(status, !!!c(`CTRLY` = "Healthy(<50)", `CTRLO` = "Healthy(>50)", `TET2` = "COVID-/TET2MT", `COVID` = "COVID+/CHIP-", `TET2C` = "COVID+/TET2MT")), levels = c("Healthy(<50)", "Healthy(>50)", "COVID-/TET2MT", "COVID+/CHIP-", "COVID+/TET2MT")))

generate_expr_plots <- function(curgene){return(
ggplot(data = plotdata1 %>% filter(gene == curgene), aes(x = status, y = expr)) +
    geom_violin(aes(fill = celltype)) +
    # ggrastr::rasterise(geom_jitter(aes(color = sample), size=0.2, stroke=0.1, shape=16), dpi = 400) +
    geom_point(data = mu1 %>% filter(gene == curgene), mapping = aes(x = status, y = expr), size=1.5, color="black") +
    scale_fill_manual(values = colors.mn) +
    geom_bracket(data = brackets %>% filter(gene == curgene), mapping = aes(label = sig), xmin = "COVID+/CHIP-", xmax = "COVID+/TET2MT", y.position = plotdata1 %>% filter(gene == curgene) %>% select(expr) %>% max()) +
    # stat_compare_means(aes(group=celltype, label = ..p.signif..), comparisons = list(c("COVID+/CHIP-", "COVID+/TET2MT"))) +
    scale_y_continuous(expand = expansion(mult = c(0.01, .09))) +
    facet_grid(cols = vars(celltype), switch="both") +
    ggtitle(curgene) +
    ylab("LogNormalized Expression") +
    theme_cowplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(), plot.title = element_text(face = "italic"))
)}

pdf(file='GeneExpr1.pdf', width=8, height=5)
wrap_plots(lapply(genes, generate_expr_plots))
dev.off()

## MULT ###

genes <- c("S100A8", "S100A9")
celltypes <- c("B", "cDC")

Idents(MULTPLUSHANIFFA) <- "celltype.status"
degs <- lapply(
    celltypes, 
    FUN = function(celltype){
        temp1 <- FindMarkers(MULTPLUSHANIFFA, features = genes, min.cells.group = 0, min.cells.feature = 0, logfc.threshold = 0, ident.1 = paste0(celltype, ".TET2C"), ident.2 = paste0(celltype, ".COVID"), min.pct = 0) %>% rownames_to_column("gene")
    }
)
names(degs) <- celltypes
degs <- bind_rows(degs, .id = "celltype")
brackets <- degs %>% mutate(sig = case_when(p_val_adj <= 0.0001 ~ "****", p_val_adj <= 0.001 ~ "***", p_val_adj <= 0.01 ~ "**", p_val_adj <= 0.05 ~ "*", TRUE ~ "ns"))

plotdata1 <- data.frame(t(MULTPLUSHANIFFA[["RNA"]]@data[genes,]), check.names = FALSE) %>% 
    mutate(status = MULTPLUSHANIFFA$group, celltype = factor(MULTPLUSHANIFFA$monaco, levels = names(colors.mn)), sample = MULTPLUSHANIFFA$sample) %>%
    filter(status %in% c("CTRLY", "CTRLO", "TET2", "COVID", "TET2C")) %>%
    filter(celltype %in% c("B", "cDC")) %>%
    pivot_longer(!c(status, celltype, sample), names_to = "gene", values_to = "expr") %>%
    mutate(status = factor(recode_factor(status, !!!c(`CTRLY` = "Healthy(<50)", `CTRLO` = "Healthy(>50)", `TET2` = "COVID-/TET2MT", `COVID` = "COVID+/CHIP-", `TET2C` = "COVID+/TET2MT")), levels = c("Healthy(<50)", "Healthy(>50)", "COVID-/TET2MT", "COVID+/CHIP-", "COVID+/TET2MT")))
mu1 <- data.frame(MULTPLUSHANIFFA_avg[["RNA"]]@data, check.names = FALSE) %>%
    rownames_to_column("gene") %>% 
    pivot_longer(!gene, names_to = c("celltype", "status"), values_to = "expr", names_sep = "\\.") %>%
    pivot_wider(names_from = "gene", values_from = "expr") %>%
    select(c("celltype", "status", all_of(genes))) %>%
    filter(status %in% c("CTRLY", "CTRLO", "TET2", "COVID", "TET2C")) %>%
    filter(celltype %in% c("B", "cDC")) %>%
    pivot_longer(!c(status, celltype), names_to = "gene", values_to = "expr") %>%
    mutate(celltype = factor(celltype, levels = c("B", "cDC"))) %>%
    mutate(status = factor(recode_factor(status, !!!c(`CTRLY` = "Healthy(<50)", `CTRLO` = "Healthy(>50)", `TET2` = "COVID-/TET2MT", `COVID` = "COVID+/CHIP-", `TET2C` = "COVID+/TET2MT")), levels = c("Healthy(<50)", "Healthy(>50)", "COVID-/TET2MT", "COVID+/CHIP-", "COVID+/TET2MT")))

generate_expr_plots <- function(curgene){return(
ggplot(data = plotdata1 %>% filter(gene == curgene), aes(x = status, y = expr)) +
    geom_violin(aes(fill = celltype)) +
    # ggrastr::rasterise(geom_jitter(aes(color = sample), size=0.2, stroke=0.1, shape=16), dpi = 400) +
    geom_point(data = mu1 %>% filter(gene == curgene), mapping = aes(x = status, y = expr), size=1.5, color="black") +
    scale_fill_manual(values = colors.mn) +
    geom_bracket(data = brackets %>% filter(gene == curgene), mapping = aes(label = sig), xmin = "COVID+/CHIP-", xmax = "COVID+/TET2MT", y.position = plotdata1 %>% filter(gene == curgene) %>% select(expr) %>% max()) +
    # stat_compare_means(aes(group=celltype, label = ..p.signif..), comparisons = list(c("COVID+/CHIP-", "COVID+/TET2MT"))) +
    scale_y_continuous(expand = expansion(mult = c(0.01, .09))) +
    facet_grid(cols = vars(celltype), switch="both") +
    ggtitle(curgene) +
    ylab("LogNormalized Expression") +
    theme_cowplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(), plot.title = element_text(face = "italic"))
)}

pdf(file='GeneExpr1_MULT.pdf', width=8, height=5)
wrap_plots(lapply(genes, generate_expr_plots))
dev.off()


###########

genes <- c("HLA-DRA", "S100A8")
celltypes <- c("CD14 Mono", "CD16 Mono", "Int Mono")

Idents(MAYOPLUSHANIFFA) <- "celltype.status"
degs <- lapply(
    celltypes, 
    FUN = function(celltype){
        temp1 <- FindMarkers(MAYOPLUSHANIFFA, features = genes, min.cells.group = 0, min.cells.feature = 0, logfc.threshold = 0, ident.1 = paste0(celltype, ".TET2C"), ident.2 = paste0(celltype, ".TET2"), min.pct = 0) %>% rownames_to_column("gene")
    }
)
names(degs) <- celltypes
degs <- bind_rows(degs, .id = "celltype")
degs1 <- lapply(
    celltypes, 
    FUN = function(celltype){
        temp1 <- FindMarkers(MAYOPLUSHANIFFA, features = genes, min.cells.group = 0, min.cells.feature = 0, logfc.threshold = 0, ident.1 = paste0(celltype, ".COVID"), ident.2 = paste0(celltype, ".TET2"), min.pct = 0) %>% rownames_to_column("gene")
    }
)
names(degs1) <- celltypes
degs1 <- bind_rows(degs1, .id = "celltype")
brackets <- degs %>% mutate(sig = case_when(p_val_adj <= 0.0001 ~ "****", p_val_adj <= 0.001 ~ "***", p_val_adj <= 0.01 ~ "**", p_val_adj <= 0.05 ~ "*", TRUE ~ "ns"))
brackets1 <- degs1 %>% mutate(sig = case_when(p_val_adj <= 0.0001 ~ "****", p_val_adj <= 0.001 ~ "***", p_val_adj <= 0.01 ~ "**", p_val_adj <= 0.05 ~ "*", TRUE ~ "ns"))

plotdata1 <- data.frame(t(MAYOPLUSHANIFFA[["RNA"]]@data[genes,]), check.names = FALSE) %>% 
    mutate(status = MAYOPLUSHANIFFA$group, celltype = factor(MAYOPLUSHANIFFA$monaco, levels = names(colors.mn)), sample = MAYOPLUSHANIFFA$sample) %>%
    filter(status %in% c("CTRLY", "CTRLO", "TET2", "COVID", "TET2C")) %>%
    filter(celltype %in% c("CD14 Mono", "CD16 Mono", "Int Mono")) %>%
    pivot_longer(!c(status, celltype, sample), names_to = "gene", values_to = "expr") %>%
    mutate(status = factor(recode_factor(status, !!!c(`CTRLY` = "Healthy(<50)", `CTRLO` = "Healthy(>50)", `TET2` = "COVID-/TET2MT", `COVID` = "COVID+/CHIP-", `TET2C` = "COVID+/TET2MT")), levels = c("Healthy(<50)", "Healthy(>50)", "COVID-/TET2MT", "COVID+/CHIP-", "COVID+/TET2MT")))
mu1 <- data.frame(MAYOPLUSHANIFFA_avg[["RNA"]]@data, check.names = FALSE) %>%
    rownames_to_column("gene") %>% 
    pivot_longer(!gene, names_to = c("celltype", "status"), values_to = "expr", names_sep = "\\.") %>%
    pivot_wider(names_from = "gene", values_from = "expr") %>%
    select(c("celltype", "status", all_of(genes))) %>%
    filter(status %in% c("CTRLY", "CTRLO", "TET2", "COVID", "TET2C")) %>%
    filter(celltype %in% c("CD14 Mono", "CD16 Mono", "Int Mono")) %>%
    pivot_longer(!c(status, celltype), names_to = "gene", values_to = "expr") %>%
    mutate(celltype = factor(celltype, levels = c("CD14 Mono", "CD16 Mono", "Int Mono"))) %>%
    mutate(status = factor(recode_factor(status, !!!c(`CTRLY` = "Healthy(<50)", `CTRLO` = "Healthy(>50)", `TET2` = "COVID-/TET2MT", `COVID` = "COVID+/CHIP-", `TET2C` = "COVID+/TET2MT")), levels = c("Healthy(<50)", "Healthy(>50)", "COVID-/TET2MT", "COVID+/CHIP-", "COVID+/TET2MT")))

generate_expr_plots <- function(curgene){return(
ggplot(data = plotdata1 %>% filter(gene == curgene), aes(x = status, y = expr)) +
    geom_violin(aes(fill = celltype)) +
    # ggrastr::rasterise(geom_jitter(aes(color = sample), size=0.2, stroke=0.1, shape=16), dpi = 400) +
    geom_point(data = mu1 %>% filter(gene == curgene), mapping = aes(x = status, y = expr), size=1.5, color="black") +
    scale_fill_manual(values = colors.mn) +
    geom_bracket(data = brackets %>% filter(gene == curgene), mapping = aes(label = sig), xmin = "COVID-/TET2MT", xmax = "COVID+/TET2MT", y.position = (plotdata1 %>% filter(gene == curgene) %>% select(expr) %>% max()) + 0.7) +
    geom_bracket(data = brackets1 %>% filter(gene == curgene), mapping = aes(label = sig), xmin = "COVID-/TET2MT", xmax = "COVID+/CHIP-", y.position = (plotdata1 %>% filter(gene == curgene) %>% select(expr) %>% max()) + 0.1) +
    # stat_compare_means(aes(group=celltype, label = ..p.signif..), comparisons = list(c("COVID+/CHIP-", "COVID+/TET2MT"))) +
    scale_y_continuous(expand = expansion(mult = c(0.01, .09))) +
    facet_grid(cols = vars(celltype), switch="both") +
    ggtitle(curgene) +
    ylab("LogNormalized Expression") +
    theme_cowplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(), plot.title = element_text(face = "italic"))
)}

pdf(file='GeneExpr2.pdf', width=10, height=5)
wrap_plots(lapply(genes, generate_expr_plots))
dev.off()

######################### Cutsite Dist ################################

celltypes <- c("CD4 T", "CD8 T", "NK", "B", "CD14 Mono", "CD16 Mono", "Int Mono", "cDC")
temp1 <- data.frame(cutsites = colSums(MULT[["mATAC"]]@counts), celltype = MULT$predicted.scRNA.monaco, status = MULT$mutation1) %>%
    filter(status %in% c("COVID", "TET2-C")) %>%
    mutate(status = factor(recode_factor(status, !!!c(`COVID` = "COVID+/CHIP-", `TET2-C` = "COVID+/TET2MT")), levels = c("COVID+/CHIP-", "COVID+/TET2MT")), celltype = factor(celltype, levels = names(colors.mn)))
temp2 <- data.frame(cutsites = colSums(MULT[["mATAC"]]@data), celltype = MULT$predicted.scRNA.monaco, status = MULT$mutation1) %>%
    filter(status %in% c("COVID", "TET2-C")) %>%
    mutate(status = factor(recode_factor(status, !!!c(`COVID` = "COVID+/CHIP-", `TET2-C` = "COVID+/TET2MT")), levels = c("COVID+/CHIP-", "COVID+/TET2MT")), celltype = factor(celltype, levels = names(colors.mn))) %>%
    filter(celltype %in% celltypes)
temp3 <- temp2 %>% group_by(celltype, status) %>% summarise(n = n())
ord <- temp2 %>% group_by(celltype, status) %>% 
    summarise(mean_cutsites = mean(cutsites)) %>% 
    pivot_wider(names_from = "status", values_from = mean_cutsites) %>% 
    mutate(avg_diff = `COVID+/TET2MT`-`COVID+/CHIP-`) %>% 
    mutate(avg_log2FC = log2(`COVID+/TET2MT`/`COVID+/CHIP-`)) %>% 
    arrange(avg_log2FC) %>% select(celltype) %>% unlist() %>% unname()

library(ggbreak)
pdf(file='cutsite-dist.pdf', width=8, height=5)
ggplot(temp2 %>% mutate(celltype = factor(celltype, levels = ord)), aes(x = status, y = cutsites, color = celltype)) +
    geom_violin() +
    geom_jitter(size=0.3, stroke=0.1, shape=16) +
    stat_summary(fun=mean, geom="point", size=0.7, color="black") +
    stat_compare_means(data = temp2 %>% mutate(celltype = factor(celltype, levels = ord)), aes(group=status, label = ..p.signif..), label.y = 155000, label.x = 1.5) +
    geom_text(data = temp3, mapping = aes(x = status, label = paste0(n)), y = -4000, angle = 25, color = "black") +
    facet_grid(cols = vars(celltype), scales="free") +
    scale_color_manual(values = colors.mn) +
    ylab("Number of cutsites") +
    # scale_y_continuous(trans='log10', labels = comma) +
    # coord_cartesian(ylim = c(0, 50000)) +
    # scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3)) +
    scale_y_break(c(50000, 150000), ticklabels=c(150000, 160000), expand = expansion(mult = c(0.2, 0))) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 25, vjust = 1, hjust=1), strip.text.x = element_text(size = 10), axis.title.x = element_blank())
dev.off()

###################### ReMap2022 Analysis (previous approach) ###########################

remap_mo <- read.table("/public_data/remap2022/remap2022_monocyte_nr_macs2_hg38_v1_0.bed") %>% 
    select(V1:V4) %>% 
    set_colnames(c("chr", "start", "end", "name")) %>%
    separate(col = "name", sep = ":", into = c("TF", "biotype"))
remap_mac <- read.table("/public_data/remap2022/remap2022_macrophage_nr_macs2_hg38_v1_0.bed") %>% 
    select(V1:V4) %>% 
    set_colnames(c("chr", "start", "end", "name")) %>%
    separate(col = "name", sep = ":", into = c("TF", "biotype"))

feats <- c(paste0("Mo.", unique(remap_mo$TF)), paste0("Mac.", unique(remap_mac$TF)))
temp <- MULT@meta.data %>% 
    select(c("mutation1", "predicted.scRNA.monaco", all_of(feats))) %>%
    filter(predicted.scRNA.monaco %in% c("CD4 T", "CD8 T", "B", "CD14 Mono", "CD16 Mono", "Int Mono", "cDC")) %>%
    filter(mutation1 %in% c("COVID", "TET2-C")) %>%
    pivot_longer(!c(mutation1, predicted.scRNA.monaco), names_to = "TF", values_to = "score") %>%
    separate(col = "TF", sep = "\\.", into = c("peakset", "TF")) %>% 
    group_by(mutation1, predicted.scRNA.monaco, peakset, TF) %>% 
    summarise(score = mean(score)) %>% 
    mutate(peakset = factor(peakset, levels = c("Mo", "Mac"))) %>% 
    mutate(TF = factor(TF, levels = c("CDK8", "CEBPB", "CREBBP", "SREBP2", "STAT3", "INTS13", "IRF1", "SMC1", "SMC1A", "MAF", "ATF2", "CHD4", "EPAS1", "HIF1A", "PPARG", "RXR", "CTCF", "NR3C1", "SPI1", "EGR1"))) %>%
    mutate(predicted.scRNA.monaco = factor(predicted.scRNA.monaco, levels = names(colors.mn))) %>%
    mutate(predicted.scRNA.monaco = factor(recode_factor(predicted.scRNA.monaco, !!!c(`CD14 Mono` = "CD14 M", `CD16 Mono` = "CD16 M", `Int Mono` = "Int M", `cDC` = "cDC")), levels = c("CD4 T", "CD8 T", "B", "CD14 M", "CD16 M", "Int M", "cDC"))) %>%
    mutate(mutation1 = factor(case_when(mutation1 == "TET2-C" ~ "COVID+/TET2MT", TRUE ~ "COVID+/CHIP-"), levels = c("COVID+/TET2MT", "COVID+/CHIP-")))

library(scico)
pdf(file='MoMac_TFs.pdf', width=10, height=7)
p1 <- ggplot(temp %>% filter(peakset == "Mo"), aes(x = TF, y = mutation1)) +
    facet_grid(rows = vars(predicted.scRNA.monaco), switch="y") +
    geom_tile(aes(fill = score), color = "black", size=0.5) +
    # scale_fill_gradient2(mid = "lightgray", low = "blue", high = "red", limits = c(-5, 15)) +
    scale_fill_scico(palette = 'roma', direction = -1, limits = c(-15, 15)) +
    scale_y_discrete(position = "right") +
    theme_cowplot() +
    ggtitle("Monocyte specific peaks from ReMap2022") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "italic"), axis.text.y = element_blank(), axis.title = element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), strip.background = element_blank(), strip.text.y = element_text(size = 12), plot.title = element_text(size = 12), axis.ticks.length=unit(0, "cm"))
p2 <- ggplot(temp %>% filter(peakset == "Mac"), aes(x = TF, y = mutation1)) +
    facet_grid(rows = vars(predicted.scRNA.monaco), switch="y") +
    geom_tile(aes(fill = score), color = "black", size=0.5) +
    # scale_fill_gradient2(mid = "lightgray", low = "blue", high = "red", limits = c(-5, 15)) +
    scale_fill_scico(palette = 'roma', direction = -1, limits = c(-15, 15)) +
    scale_y_discrete(position = "right") +
    theme_cowplot() +
    ggtitle("Macrophage specific peaks from ReMap2022") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "italic"), axis.title = element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), strip.background = element_blank(), strip.text.y = element_blank(), plot.title = element_text(size = 12), axis.ticks.length=unit(0, "cm"))
design = "
AAAAAABBBB
"
plot(p1 + p2 + plot_layout(design = design))
dev.off()

###################### ReMap2022 Analysis (current approach) ###########################

remap_my <- read.table("/public_data/remap2022/remap2022_myeloid_nr_macs2_hg38_v1_0.bed") %>% 
    select(V1:V4) %>% 
    set_colnames(c("chr", "start", "end", "name")) %>%
    separate(col = "name", sep = ":", into = c("TF", "biotype"))

DefaultAssay(MULT) <- "mATAC"
library(BiocParallel)
register(SerialParam())
MULT_peaks <- data.frame(Peak = rownames(MULT)) %>% separate(col = "Peak", sep = "-", into = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()

TFs <- unique(remap_my$TF)
TF_features <- lapply(TFs, function(cur_TF){
    cur_peaks <- makeGRangesFromDataFrame(remap_my %>% filter(TF == cur_TF))
    cur_peaks <- keepStandardChromosomes(cur_peaks, pruning.mode="tidy")
    ovp_features <- suppressWarnings(GRangesToString(MULT_peaks[case_when(is.na(findOverlaps(MULT_peaks, cur_peaks, select = "first", ignore.strand=TRUE)) ~ FALSE, TRUE ~ TRUE)]))
    return(ovp_features)
})
names(TF_features) <- TFs
features <- TF_features

library(fastmatch)
object <- MULT
assay <- "mATAC"
# first find index of each feature
feat.idx <- sapply(X = features, FUN = fmatch, rownames(x = object[[assay]]))
j <- sapply(X = seq_along(along.with = features), FUN = function(x) {
    rep(x = x, length(x = features[[x]]))
})

# construct sparse matrix with features
mat <- sparseMatrix(
    i = unlist(x = feat.idx, use.names = FALSE),
    j = unlist(x = j, use.names = FALSE),
    x = 1,
    dims = c(nrow(x = object[[assay]]), length(x = features))
)
rownames(x = mat) <- rownames(x = object[[assay]])
colnames(x = mat) <- names(x = features)

# run chromVAR
cv <- RunChromVAR(
    object = object[[assay]],
    motif.matrix = mat,
    genome = genome
)
MULT[["remap"]] <- cv

########################  Differentially enriched ReMap protein bindings  ##################################

MULT$celltype.status <- paste(MULT$predicted.scRNA.monaco, MULT$mutation1, sep = ".")
Idents(MULT) <- "celltype.status"
DefaultAssay(MULT) <- 'remap'

types <- unique(MULT$predicted.scRNA.monaco)

dars <- lapply(
    types, 
    FUN = function(celltype){
        curridents <- unique(Idents(MULT))
        if (!(paste0(celltype, ".TET2-C") %in% curridents) | !(paste0(celltype, ".COVID") %in% curridents)) return(data.frame())
        if((sum(MULT$celltype.status == paste0(celltype, ".TET2-C")) < 3) | (sum(MULT$celltype.status == paste0(celltype, ".COVID")) < 3)) return(data.frame())
        temp1 <- FindMarkers(MULT, ident.1 = paste0(celltype, ".TET2-C"), ident.2 = paste0(celltype, ".COVID"), mean.fxn = rowMeans, fc.name = "avg_diff") %>% rownames_to_column("motif")
        if(is.na(temp1)){ return(data.frame()) } else { return(temp1) }
    }
)
names(dars) <- types
dars <- bind_rows(dars, .id = "celltype") 
write.table(dars, file = "dars_TET2C_COVID.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

DefaultAssay(MULT) <- "remap"
samples <- MULT@meta.data %>% filter(mutation1 %in% c("COVID", "TET2-C")) %>% select(sample) %>% unlist() %>% unname() %>% unique()
dars_robust <- lapply(
    samples, 
    FUN = function(sam){
        print(sam)
        oneout <- subset(MULT, subset = sample != sam)
        curridents <- unique(Idents(oneout))
        ctsp_dars_list <- lapply(
            types,
            FUN = function(celltype){
                print(celltype)
                if (!(paste0(celltype, ".TET2-C") %in% curridents) | !(paste0(celltype, ".COVID") %in% curridents)) return(c())
                if (length(WhichCells(oneout, idents = paste0(celltype, ".TET2-C"))) < 3 | length(WhichCells(oneout, idents = paste0(celltype, ".COVID"))) < 3) return(c())
                temp1 <- FindMarkers(oneout, ident.1 = paste0(celltype, ".TET2-C"), ident.2 = paste0(celltype, ".COVID"), mean.fxn = rowMeans, fc.name = "avg_diff") %>% rownames_to_column("motif")
                if(nrow(temp1) == 0){ return(c()) } else if (is.na(temp1)) { return(c()) } else { return((temp1 %>% filter(p_val_adj < 0.05))$motif) }
            }
        )
        names(ctsp_dars_list) <- types
        return(ctsp_dars_list)
    }
)
names(dars_robust) <- samples
dars_robust <- lapply(dars_robust %>% transpose(), function(x){Reduce(intersect, x)})
darsr <- bind_rows(lapply(dars_robust, function(motif){data.frame(motif)}), .id = "celltype") %>% left_join(dars)
write.table(darsr, file = "darsr_TET2C_COVID.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

dars <- read.table(file = "dars_TET2C_COVID.tsv", sep = "\t", skip = 1) %>% set_colnames(c("celltype", "motif", "p_val", "avg_diff", "pct.1", "pct.2", "p_val_adj"))
darsr <- read.table(file = "darsr_TET2C_COVID.tsv", sep = "\t", skip = 1) %>% set_colnames(c("celltype", "motif", "p_val", "avg_diff", "pct.1", "pct.2", "p_val_adj"))

temp1 <- dars %>% 
    filter(p_val_adj < 0.05) %>%
    mutate(celltype = factor(celltype, levels = names(colors.mn))) %>% 
    # mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_p_val_adj = -log10(p_val_adj))
# forlabels <- temp1 %>% mutate(TF = case_when((TF %in% c("EGR1")) ~ TF, TRUE ~ ""))
pdf(file = "dars.pdf", width = 8, height = 8)
p1 <- ggplot(temp1, aes(x = avg_diff, y = neg_log10_p_val_adj)) +
    geom_point(aes(color = celltype), size = 1.5, stroke=0.1, shape = 16) +
    # geom_text_repel(data = temp1 %>% filter(neg_log10_p_val_adj > 77), aes(label = TF), size = 3, max.overlaps = 15, fontface = "italic") +
    geom_text_repel(data = temp1 %>% filter(celltype == "CD14 Mono"), aes(label = motif), size = 3, max.overlaps = Inf, fontface = "italic") +
    # geom_text_repel(data = forlabels, aes(label = TF), size = 3, max.overlaps = Inf, fontface = "italic") +
    # coord_cartesian(xlim = c(-5, 5)) +
    labs(x = "Average difference in chromVAR score", y = "-log10(Adjusted P value)") +
    scale_color_manual(values = colors.mn) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.title = element_blank())
plot(p1 + plot_annotation(theme = theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))))
dev.off()

temp1 <- darsr %>% 
    filter(p_val_adj < 0.05) %>%
    mutate(celltype = factor(celltype, levels = names(colors.mn))) %>% 
    # mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_p_val_adj = -log10(p_val_adj))
# forlabels <- temp1 %>% mutate(TF = case_when((TF %in% c("EGR1")) ~ TF, TRUE ~ ""))
pdf(file = "darsr.pdf", width = 8, height = 8)
p1 <- ggplot(temp1, aes(x = avg_diff, y = neg_log10_p_val_adj)) +
    geom_point(aes(color = celltype), size = 5, stroke=0.1, shape = 16) +
    # geom_text_repel(data = temp1 %>% filter(neg_log10_p_val_adj > 77), aes(label = TF), size = 3, max.overlaps = 15, fontface = "italic") +
    geom_text_repel(data = temp1, aes(label = motif), size = 4, max.overlaps = Inf, fontface = "italic") +
    # geom_text_repel(data = forlabels, aes(label = TF), size = 3, max.overlaps = Inf, fontface = "italic") +
    # coord_cartesian(xlim = c(-5, 5)) +
    labs(x = "Average difference in chromVAR score", y = "-log10(Adjusted P value)") +
    scale_color_manual(values = colors.mn[intersect(names(colors.mn), temp1$celltype)]) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.title = element_blank())
    
plot(p1 + plot_annotation(theme = theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))))
dev.off()

########################  Differentially enriched motifs ##################################

MULT$celltype.status <- paste(MULT$predicted.scRNA.monaco, MULT$mutation1, sep = ".")
Idents(MULT) <- "celltype.status"
DefaultAssay(MULT) <- 'chromvar'

types <- unique(MULT$predicted.scRNA.monaco)

dams <- lapply(
    types, 
    FUN = function(celltype){
        curridents <- unique(Idents(MULT))
        if (!(paste0(celltype, ".TET2-C") %in% curridents) | !(paste0(celltype, ".COVID") %in% curridents)) return(data.frame())
        if((sum(MULT$celltype.status == paste0(celltype, ".TET2-C")) < 3) | (sum(MULT$celltype.status == paste0(celltype, ".COVID")) < 3)) return(data.frame())
        temp1 <- FindMarkers(MULT, ident.1 = paste0(celltype, ".TET2-C"), ident.2 = paste0(celltype, ".COVID"), mean.fxn = rowMeans, fc.name = "avg_diff") %>% rownames_to_column("motif")
        if(is.na(temp1)){ return(data.frame()) } else { return(temp1) }
    }
)
names(dams) <- types
dams <- bind_rows(dams, .id = "celltype") 
DefaultAssay(MULT) <- "mATAC"
dams <- dams %>% mutate(TF = ConvertMotifID(MULT, id = dams$motif))
write.table(dams, file = "dams_TET2C_COVID.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

DefaultAssay(MULT) <- "chromvar"
samples <- MULT@meta.data %>% filter(mutation1 %in% c("COVID", "TET2-C")) %>% select(sample) %>% unlist() %>% unname() %>% unique()
dams_robust <- lapply(
    samples, 
    FUN = function(sam){
        print(sam)
        oneout <- subset(MULT, subset = sample != sam)
        curridents <- unique(Idents(oneout))
        ctsp_dams_list <- lapply(
            types,
            FUN = function(celltype){
                print(celltype)
                if (!(paste0(celltype, ".TET2-C") %in% curridents) | !(paste0(celltype, ".COVID") %in% curridents)) return(c())
                if (length(WhichCells(oneout, idents = paste0(celltype, ".TET2-C"))) < 3 | length(WhichCells(oneout, idents = paste0(celltype, ".COVID"))) < 3) return(c())
                temp1 <- FindMarkers(oneout, ident.1 = paste0(celltype, ".TET2-C"), ident.2 = paste0(celltype, ".COVID"), mean.fxn = rowMeans, fc.name = "avg_diff") %>% rownames_to_column("motif")
                if(is.na(temp1)){ return(c()) } else { return((temp1 %>% filter(p_val_adj < 0.05))$motif) }
            }
        )
        names(ctsp_dams_list) <- types
        return(ctsp_dams_list)
    }
)
names(dams_robust) <- samples
dams_robust <- lapply(dams_robust %>% transpose(), function(x){Reduce(intersect, x)})
damsr <- bind_rows(lapply(dams_robust, function(motif){data.frame(motif)}), .id = "celltype") %>% left_join(dams)
write.table(damsr, file = "damsr_TET2C_COVID.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

dams <- read.table(file = "dams_TET2C_COVID.tsv", sep = "\t", skip = 1) %>% set_colnames(c("celltype", "motif", "p_val", "avg_diff", "pct.1", "pct.2", "p_val_adj", "TF"))
damsr <- read.table(file = "damsr_TET2C_COVID.tsv", sep = "\t", skip = 1) %>% set_colnames(c("celltype", "motif", "p_val", "avg_diff", "pct.1", "pct.2", "p_val_adj", "TF"))

temp1 <- damsr %>% 
    filter(p_val_adj < 0.05) %>%
    mutate(celltype = factor(celltype, levels = names(colors.mn))) %>% 
    # mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_p_val_adj = -log10(p_val_adj))
# forlabels <- temp1 %>% mutate(TF = case_when((TF %in% c("EGR1")) ~ TF, TRUE ~ ""))
pdf(file = "damsr.pdf", width = 8, height = 8)
p1 <- ggplot(temp1, aes(x = avg_diff, y = neg_log10_p_val_adj)) +
    geom_point(aes(color = celltype), size = 1.5, stroke=0.1, shape = 16) +
    # geom_text_repel(data = temp1 %>% filter(neg_log10_p_val_adj > 77), aes(label = TF), size = 3, max.overlaps = 15, fontface = "italic") +
    geom_text_repel(data = temp1 %>% filter(celltype == "CD14 Mono"), aes(label = TF), size = 3, max.overlaps = Inf, fontface = "italic") +
    # geom_text_repel(data = forlabels, aes(label = TF), size = 3, max.overlaps = Inf, fontface = "italic") +
    coord_cartesian(xlim = c(-5, 5)) +
    labs(x = "Average difference in chromVAR score", y = "-log10(Adjusted P value)") +
    scale_color_manual(values = colors.mn) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.title = element_blank())
    
plot(p1 + plot_annotation(theme = theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))))
dev.off()

########################  Differentially accessible peaks ##################################

MULT$celltype.status <- paste(MULT$predicted.scRNA.monaco, MULT$mutation1, sep = ".")
Idents(MULT) <- "celltype.status"
DefaultAssay(MULT) <- 'mATAC'

types <- unique(MULT$predicted.scRNA.monaco)

daps <- lapply(
    types, 
    FUN = function(celltype){
        curridents <- unique(Idents(MULT))
        if (!(paste0(celltype, ".TET2-C") %in% curridents) | !(paste0(celltype, ".COVID") %in% curridents)) return(data.frame())
        if((sum(MULT$celltype.status == paste0(celltype, ".TET2-C")) < 3) | (sum(MULT$celltype.status == paste0(celltype, ".COVID")) < 3)) return(data.frame())
        temp1 <- FindMarkers(MULT, ident.1 = paste0(celltype, ".TET2-C"), ident.2 = paste0(celltype, ".COVID"), min.pct = 0, test.use = 'LR', latent.vars = 'atac_peak_region_fragments') %>% rownames_to_column("peak")
        if(is.na(temp1)){ return(data.frame()) } else { return(temp1) }
    }
)
names(daps) <- types
daps <- bind_rows(daps, .id = "celltype")
write.table(daps, file = "daps.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

daps <- read.table(file = "daps.tsv", sep = "\t", skip = 1) %>% set_colnames(c("celltype", "peak", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"))
temp1 <- daps %>% 
    filter(p_val < 0.05) %>% 
    mutate(celltype = factor(celltype, levels = names(colors.mn))) %>% 
    mutate(neg_log10_pval = -log10(p_val)) %>%
    filter(celltype %in% c("CD14 Mono") & avg_log2FC < 0)
pdf(file = "daps.pdf", width = 8, height = 8)
p1 <- ggplot(temp1, aes(x = avg_log2FC, y = neg_log10_pval)) +
    geom_point(aes(color = celltype), size = 0.5, stroke=0.1, shape = 16) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    annotate("text", label = "P value = 0.05", x = 1.5, y = -log10(0.05)+0.1, size = 3) +
    # geom_text_repel(data = temp1 %>% filter(neg_log10_p_val_adj > 77), aes(label = TF), size = 3, max.overlaps = 15, fontface = "italic") +
    # geom_text_repel(data = temp1 %>% filter(celltype == "CD14 Mono"), aes(label = TF), size = 3, max.overlaps = Inf, fontface = "italic") +
    # coord_cartesian(xlim = c(-5, 5)) +
    # labs(x = "Average difference in chromVAR score", y = "-log10(Adjusted P value)") +
    scale_color_manual(values = colors.mn) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot()
plot(p1)
dev.off()

temp1 <- daps %>% 
    filter(p_val < 0.05) %>% 
    mutate(celltype = factor(celltype, levels = names(colors.mn))) %>% 
    mutate(neg_log10_pval = -log10(p_val)) %>% 
    mutate(direction = case_when(avg_log2FC > 0 ~ "UP", TRUE ~ "DOWN"))

pdf(file = "daps_props.pdf", width = 3.75, height = 4)
# %>% filter(!(celltype %in% c("CD4 T", "gdT", "NK", "Basophils")))
p1 <- ggplot(temp1 %>% select(celltype, direction) %>% table() %>% data.frame() %>% mutate(celltype = factor(celltype, levels = rev(names(colors.mn)))) %>% filter(!(celltype %in% c("Basophils")))) +
    geom_bar(aes(fill=direction, y=Freq, x=celltype), position=position_fill(reverse = TRUE), stat="identity") +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c("darkblue", "darkred")) +
    coord_flip() + 
    theme_cowplot()
plot(p1)
dev.off()

############### CoveragePlots ######################

# sort -k 1,1 -k2,2n /public_data/EPICarrayCpGs.bed > /public_data/EPICarrayCpGs_sorted_hg38.bed
EPICarrayCpGs <- read.table("/public_data/EPICarrayCpGs_sorted_hg38.bed") %>% set_colnames(c("chr", "start", "end")) %>% makeGRangesFromDataFrame()

DefaultAssay(MULT) <- "mATAC"
Idents(MULT) <- factor(MULT$predicted.scRNA.monaco, levels = names(colors.mn))
MULT <- LinkPeaks(MULT, peak.assay = "mATAC", expression.assay = "RNA", genes.use = c("EGR1", "MALAT1", "NEAT1"))
MULT_COVID <- subset(MULT, subset = (mutation1 == "COVID"))
Idents(MULT_COVID) <- factor(MULT_COVID$predicted.scRNA.monaco, levels = names(colors.mn))
MULT_COVID <- LinkPeaks(MULT_COVID, peak.assay = "mATAC", expression.assay = "RNA", genes.use = c("EGR1", "MALAT1", "NEAT1"))

myCoveragePlot <- function(obj, reg, gene, ids){
    cov_plot <- CoveragePlot(obj, region = reg, features = gene, annotation = FALSE, peaks = FALSE, links = FALSE, idents = ids, scale.factor = 1e7, ymax = 180) & scale_fill_manual(values = colors.mn)
    expr_plot <- ExpressionPlot(obj, features = gene, assay = "RNA", idents = ids) & scale_fill_manual(values = colors.mn)
    gene_plot <- AnnotationPlot(obj, region = reg)
    peak_plot <- PeakPlot(obj, region = reg)
    CpG_plot <- PeakPlot(obj, region = reg, peaks = EPICarrayCpGs, color = "blue")
    link_plot <- LinkPlot(obj, region = reg)
    bw1 <- BigwigTrack(region = StringToGRanges(reg), bigwig = "/downloads/GSM4042683_EGR1-MONO_hg38.bw.bw", y_label = "EGR1 (MONO)", downsample.rate = 1)
    bw2 <- BigwigTrack(region = StringToGRanges(reg), bigwig = "/downloads/GSM4042681_EGR1-MACRO-D4_hg38.bw.bw", y_label = "EGR1 (MACRO)", downsample.rate = 1)
    bw3 <- BigwigTrack(region = StringToGRanges(reg), bigwig = "/downloads/GSM4042682_EGR1-MACRO-D4-CUTRUN_hg38.bw.bw", y_label = "EGR1 (MACRO_CR)", downsample.rate = 1)
    bw4 <- BigwigTrack(region = StringToGRanges(reg), bigwig = "/downloads/GSM3686934_MO_PU.1-ChIP_donorO_hg38.bw.bw", y_label = "SPI1 (MONO)", downsample.rate = 1)
    bw5 <- BigwigTrack(region = StringToGRanges(reg), bigwig = "/downloads/GSM3686953_MAC7d_PU.1-ChIP_donorO_hg38.bw.bw", y_label = "SPI1 (MACRO)", downsample.rate = 1)
    return(CombineTracks(plotlist = list(cov_plot, peak_plot, CpG_plot, gene_plot, link_plot, bw1, bw2, bw3, bw4, bw5), expression.plot = expr_plot, heights = c(10, 1, 1, 2, 3, 3, 3, 3, 3, 3), widths = c(10, 1)))
}

idents_sub <- c("CD4 T", "CD8 T", "NK", "CD14 Mono", "CD16 Mono", "Int Mono")
pdf("coverage_around_genes.pdf", width=15, height=18)
# p1 <- myCoveragePlot(MULT, "chr5-138365479-138569303", "EGR1", idents_sub)
# # p2 <- myCoveragePlot(MULT_COVID, "chr5-138365479-138569303", "EGR1", idents_sub)
# plot(p1)
p1 <- myCoveragePlot(MULT, "chr11-65397688-65606431", "MALAT1", idents_sub)
# p2 <- myCoveragePlot(MULT_COVID, "chr11-65397688-65606431", "MALAT1", idents_sub)
plot(p1)
# p1 <- myCoveragePlot(MULT, "chr11-65323096-65526513", "NEAT1", idents_sub)
# # p2 <- myCoveragePlot(MULT_COVID, "chr11-65323096-65526513", "NEAT1", idents_sub)
# plot(p1)
dev.off()

DefaultAssay(MULT) <- "mATAC"
MULT_COVID <- subset(MULT, subset = (mutation1 == "COVID"))
Idents(MULT_COVID) <- factor(MULT_COVID$predicted.scRNA.monaco, levels = names(colors.mn))
MULT_TET2C <- subset(MULT, subset = (mutation1 == "TET2-C"))
Idents(MULT_TET2C) <- factor(MULT_TET2C$predicted.scRNA.monaco, levels = names(colors.mn))
MULT_COVID <- LinkPeaks(MULT_COVID, peak.assay = "mATAC", expression.assay = "RNA", genes.use = c("EGR1", "MALAT1", "NEAT1"), method = "pearson", min.cells = 50, score_cutoff = 0.05)
MULT_TET2C <- LinkPeaks(MULT_TET2C, peak.assay = "mATAC", expression.assay = "RNA", genes.use = c("EGR1", "MALAT1", "NEAT1"), method = "pearson", min.cells = 5, score_cutoff = 0.025)

idents_sub <- c("CD4 T", "CD8 T", "NK", "CD14 Mono", "CD16 Mono", "Int Mono")
pdf("coverage_around_genes.pdf", width=30, height=10)
p1 <- CoveragePlot(MULT_COVID, region = "EGR1", features = "EGR1", extend.upstream = 100000, extend.downstream = 100000, expression.assay = "RNA", idents = idents_sub, scale.factor = 1e7, ymax = 100, window = 150) & scale_fill_manual(values = colors.mn)
p2 <- CoveragePlot(MULT_TET2C, region = "EGR1", features = "EGR1", extend.upstream = 100000, extend.downstream = 100000, expression.assay = "RNA", idents = idents_sub, scale.factor = 1e7, ymax = 100, window = 150) & scale_fill_manual(values = colors.mn)
plot(p1 | p2)
p1 <- CoveragePlot(MULT_COVID, region = "MALAT1", features = "MALAT1", extend.upstream = 100000, extend.downstream = 100000, expression.assay = "RNA", idents = idents_sub, scale.factor = 1e7, ymax = 100, window = 150) & scale_fill_manual(values = colors.mn)
p2 <- CoveragePlot(MULT_TET2C, region = "MALAT1", features = "MALAT1", extend.upstream = 100000, extend.downstream = 100000, expression.assay = "RNA", idents = idents_sub, scale.factor = 1e7, ymax = 100, window = 150) & scale_fill_manual(values = colors.mn)
plot(p1 | p2)
p1 <- CoveragePlot(MULT_COVID, region = "NEAT1", features = "NEAT1", extend.upstream = 100000, extend.downstream = 100000, expression.assay = "RNA", idents = idents_sub, scale.factor = 1e7, ymax = 100, window = 150) & scale_fill_manual(values = colors.mn)
p2 <- CoveragePlot(MULT_TET2C, region = "NEAT1", features = "NEAT1", extend.upstream = 100000, extend.downstream = 100000, expression.assay = "RNA", idents = idents_sub, scale.factor = 1e7, ymax = 100, window = 150) & scale_fill_manual(values = colors.mn)
plot(p1 | p2)
dev.off()


############################### Motif enrichment in differentially accessible peaks ########################################

temp1 <- daps %>% 
    filter(p_val < 0.05) %>% 
    mutate(celltype = factor(celltype, levels = names(colors.mn))) %>% 
    mutate(neg_log10_pval = -log10(p_val)) %>% 
    # filter(avg_log2FC < 0) 
    filter(celltype %in% c("CD14 Mono") & avg_log2FC < 0) 
enriched.motifs <- FindMotifs(MULT, features = temp1 %>% select(peak) %>% unlist() %>% unname(), background = NULL) %>%
    filter(pvalue < 0.05) %>%
    mutate(neg_log10_pval = -log10(pvalue))
head(enriched.motifs, n = 20)
write.table(enriched.motifs, file = "enriched.motifs.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
enriched.motifs <- read.table("enriched.motifs.tsv", skip = 1) %>% set_colnames(c("motif", "observed", "background", "percent.observed", "percent.background", "fold.enrichment", "pvalue", "motif.name", "neg_log10_pval"))

pdf(file = "CD14Mono_lostpeaks_enriched_motifs_1.pdf", width = 10, height = 10)
ggplot(enriched.motifs, aes(x = fold.enrichment, y = neg_log10_pval)) +
    geom_point(size=1.5, stroke=0.1, shape=16) +
    geom_text_repel(data = enriched.motifs %>% slice_head(n = 50), aes(label = motif.name), size = 5, fontface = "italic") +
    labs(x = "Fold enrichment", y = "-log10(P value)") +
    theme_cowplot()
dev.off()

pdf(file = "CD14Mono_lostpeaks_enriched_motifs_2.pdf", width = 10, height = 10)
MotifPlot(MULT, motifs = enriched.motifs$motif[1:20]) + theme(strip.text = element_text(face = "italic"))
dev.off()

############## Expression DotPlots ##############################

Idents(GEX) <- "monaco"
# genes <- setdiff(c(ConvertMotifID(MULT, id = enriched.motifs$motif[1:20])), c("SP9"))
genes <- c("CDK8", "CEBPB", "CREBBP", "SREBP2", "STAT3", "INTS13", "IRF1", "SMC1", "SMC1A", "MAF", "ATF2", "CHD4", "EPAS1", "HIF1A", "PPARG", "RXR", "CTCF", "NR3C1", "SPI1", "EGR1")
idents_sub <- c("CD4 T", "CD8 T", "NK", "B", "CD14 Mono", "CD16 Mono", "Int Mono", "cDC")
GEX_COVID <- subset(GEX, subset = mutation1 %in% c("COVID"))
GEX_TET2C <- subset(GEX, subset = mutation1 %in% c("TET2-C"))
p1 <- DotPlot(GEX_COVID, assay = "RNA", features = genes, idents = idents_sub, cluster.idents = FALSE, scale.max = 60, scale.min = 0) +
    scale_color_scico(palette = 'vik', direction = -1, limits = c(-2.5, 2.5), begin = 0, end = 0.5) +
    scale_radius(breaks = c(20, 40, 60)) +
    guides(color = guide_colorbar(order = 1, title = 'Scaled AvgExpr'), radius = guide_legend(order = 2)) +
    ggtitle("COVID") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face="italic", size=20), axis.text.y = element_text(size=20), axis.title = element_blank())
p2 <- DotPlot(GEX_TET2C, assay = "RNA", features = genes, idents = idents_sub, cluster.idents = FALSE, scale.max = 60, scale.min = 0) +
    scale_color_scico(palette = 'vik', direction = 1, limits = c(-2.5, 2.5), begin = 0.5, end = 1) +
    scale_radius(breaks = c(20, 40, 60)) +
    guides(color = guide_colorbar(order = 1, title = 'Scaled AvgExpr'), radius = guide_legend(order = 2)) +
    ggtitle("TET2+COVID") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face="italic", size=20), axis.text.y = element_text(size=20), axis.title = element_blank())
pdf(file='TFs_dotplot_GEX.pdf', width=26, height=4)
p1 | p2
dev.off()


####################### Cicero connections by CpG sites ################


# Run cicero_script_1.R
# cd /FINAL/scMULT/
# sort -k 1,1 -k2,2n /public_data/EPICarrayCpGs.bed > /public_data/EPICarrayCpGs_sorted_hg38.bed
# bedtools intersect -a MULT_peaks.bed -b /public_data/CpG_Islands_hg38.bed -u > MULT_peaks_CpG_islands.bed
# bedtools intersect -a MULT_peaks.bed -b /public_data/EPICarrayCpGs_sorted_hg38.bed -u > MULT_peaks_CpG.bed
# bedtools intersect -a MULT_peaks.bed -b /public_data/EPICarrayCpGs_sorted_hg38.bed -wa > MULT_peaks_CpG_wa.bed
# bedtools intersect -a /public_data/EPICarrayCpGs_sorted_hg38.bed -b MULT_peaks.bed -u > CpG_MULT_peaks.bed

MULT_peaks_CpG_wa <- read.table(file = "/FINAL/scMULT/MULT_peaks_CpG_wa.bed", sep = "\t", skip = 1) %>% 
    set_colnames(c("chr", "start", "end")) %>% unite(everything(), col = "peak", sep = "-", remove = FALSE)

CpG_counts <- MULT_peaks_CpG_wa %>% count(peak) %>% separate(peak, into = c("chr", "start", "end"), sep = "-", remove = FALSE) %>% mutate(len = as.integer(end)-as.integer(start))

conns.list <- readRDS("/FINAL/scMULT/conns.list.CD14Mono_COVID_TET2C.rds")

COVID_conns <- conns.list[["CD14 Mono.COVID"]] %>% drop_na(coaccess) %>% filter(coaccess > 0)
TET2C_conns <- conns.list[["CD14 Mono.TET2-C"]] %>% drop_na(coaccess) %>% filter(coaccess > 0)
conns_per_peak <- COVID_conns %>% count(Peak1) %>% mutate(status = "COVID") %>%
    bind_rows(TET2C_conns %>% count(Peak1) %>% mutate(status = "TET2-C")) %>%
    mutate(status = factor(recode_factor(status, !!!c(`COVID` = "COVID+/CHIP-", `TET2-C` = "COVID+/TET2MT")), levels = c("COVID+/CHIP-", "COVID+/TET2MT")))

temp <- CpG_counts %>% select(peak) %>% unlist() %>% unname()
conns_per_peak_CpG <- COVID_conns %>% filter(Peak1 %in% temp | Peak2 %in% temp) %>% count(Peak1) %>% mutate(status = "COVID") %>%
    bind_rows(TET2C_conns %>% filter(Peak1 %in% temp | Peak2 %in% temp) %>% count(Peak1) %>% mutate(status = "TET2-C")) %>%
    mutate(status = factor(recode_factor(status, !!!c(`COVID` = "COVID+/CHIP-", `TET2-C` = "COVID+/TET2MT")), levels = c("COVID+/CHIP-", "COVID+/TET2MT")))

temp <- CpG_counts %>% select(peak) %>% unlist() %>% unname()
conns_per_peak_CpG_not <- COVID_conns %>% filter(!(Peak1 %in% temp | Peak2 %in% temp)) %>% count(Peak1) %>% mutate(status = "COVID") %>%
    bind_rows(TET2C_conns %>% filter(!(Peak1 %in% temp | Peak2 %in% temp)) %>% count(Peak1) %>% mutate(status = "TET2-C")) %>%
    mutate(status = factor(recode_factor(status, !!!c(`COVID` = "COVID+/CHIP-", `TET2-C` = "COVID+/TET2MT")), levels = c("COVID+/CHIP-", "COVID+/TET2MT")))

pdf("connhisto-CpG-TET2.pdf", width=6, height=5)
p1 <- ggplot(conns_per_peak, aes(x=status, y=n)) + 
    geom_violin(aes(fill = status), show.legend = FALSE, alpha = 0.6) +
    stat_summary(fun=median, geom="point", size=1.5, color="black") +
    stat_compare_means(aes(group=status, label = ..p.signif..), label.y = 100, label.x = 1.5) +
    scale_fill_manual(values = c("#56B4E9", "#009E73")) +
    ylab("Connections per peak") +
    coord_cartesian(ylim = c(0, 100)) +
    theme_cowplot() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 25, vjust = 1, hjust=1))
p2 <- ggplot(conns_per_peak_CpG, aes(x=status, y=n)) + 
    geom_violin(aes(fill = status), show.legend = FALSE, alpha = 0.6) +
    stat_summary(fun=median, geom="point", size=1.5, color="black") +
    stat_compare_means(aes(group=status, label = ..p.signif..), label.y = 100, label.x = 1.5) +
    scale_fill_manual(values = c("#56B4E9", "#009E73")) +
    ylab("Connections with CpG site") +
    coord_cartesian(ylim = c(0, 100)) +
    theme_cowplot() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 25, vjust = 1, hjust=1))
p3 <- ggplot(conns_per_peak_CpG_not, aes(x=status, y=n)) + 
    geom_violin(aes(fill = status), show.legend = FALSE, alpha = 0.6) +
    stat_summary(fun=median, geom="point", size=1.5, color="black") +
    stat_compare_means(aes(group=status, label = ..p.signif..), label.y = 100, label.x = 1.5) +
    scale_fill_manual(values = c("#56B4E9", "#009E73")) +
    ylab("Connections without CpG site") +
    coord_cartesian(ylim = c(0, 100)) +
    theme_cowplot() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 25, vjust = 1, hjust=1))
plot(p1 | p2 | p3)
dev.off()

# library(coin)
# median_test(n~status, data = conns_per_peak)
# median_test(n~status, data = conns_per_peak_CpG)
# median_test(n~status, data = conns_per_peak_CpG_not)

# disc_ks_test((conns_per_peak %>% filter(status == "COVID"))$n, ecdf((conns_per_peak %>% filter(status == "TET2-C"))$n))
# disc_ks_test((conns_per_peak_CpG_not %>% filter(status == "COVID"))$n, ecdf((conns_per_peak_CpG_not %>% filter(status == "TET2-C"))$n))

########################## Comparison of DEGs of various comparisons ####################

degsr_TET2C_COVID <- read.table(file = "/FINAL/scRNA/degsr_TET2C_COVID.tsv", sep = "\t", skip = 1) %>% 
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>%
    filter(p_val_adj < 0.05) %>%
    mutate(direction = ifelse(avg_log2FC > 0, "UP", "DOWN")) %>%
    mutate(celltype = factor(celltype, levels = names(colors.mn)))
degsr_TET2_CTRLY <- read.table(file = "/FINAL/scRNA/degsr_TET2_CTRLY.tsv", sep = "\t", skip = 1) %>% 
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>%
    filter(p_val_adj < 0.05) %>%
    mutate(direction = ifelse(avg_log2FC > 0, "UP", "DOWN")) %>%
    mutate(celltype = factor(celltype, levels = names(colors.mn)))


degsr_TET2C_COVID <- degsr_TET2C_COVID %>% 
    filter(p_val_adj < 0.05) %>% 
    mutate(direction = ifelse(avg_log2FC > 0, "UP", "DOWN")) %>%
    mutate(celltype = factor(celltype, levels = names(colors.mn))) %>% 
    mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_adj_pval = -log10(p_val_adj))
degsr_TET2_CTRLY <- degsr_TET2_CTRLY %>% 
    filter(p_val_adj < 0.05) %>% 
    mutate(direction = ifelse(avg_log2FC > 0, "UP", "DOWN")) %>%
    mutate(celltype = factor(celltype, levels = names(colors.mn))) %>% 
    mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_adj_pval = -log10(p_val_adj))

temp <- degsr_TET2C_COVID %>% 
    left_join(degsr_TET2_CTRLY %>% select(celltype, gene, direction, avg_log2FC) %>% rename(CHIPCTRL = direction, avg_log2FC_COVneg = avg_log2FC)) %>% 
    drop_na(CHIPCTRL) %>% 
    rename(avg_log2FC_COVpos = avg_log2FC) %>%
    mutate(avg_FC_COVpos = 2^avg_log2FC_COVpos, avg_FC_COVneg = 2^avg_log2FC_COVneg) %>%
    filter(celltype == "CD14 Mono")
temp1 <- degsr_TET2C_COVID %>% 
    anti_join(degsr_TET2_CTRLY %>% select(celltype, gene, direction, avg_log2FC) %>% rename(CHIPCTRL = direction, avg_log2FC_COVneg = avg_log2FC)) %>% 
    rename(avg_log2FC_COVpos = avg_log2FC) %>%
    mutate(avg_FC_COVpos = 2^avg_log2FC_COVpos) %>%
    filter(celltype == "CD14 Mono")
temp2 <- degsr_TET2_CTRLY %>% 
    anti_join(degsr_TET2C_COVID %>% select(celltype, gene, direction, avg_log2FC) %>% rename(TET2CCOVID = direction, avg_log2FC_COVneg = avg_log2FC)) %>% 
    rename(avg_log2FC_COVpos = avg_log2FC) %>%
    mutate(avg_FC_COVpos = 2^avg_log2FC_COVpos) %>%
    filter(celltype == "CD14 Mono")

pdf(file = "diamond.pdf", width = 5, height = 9)
forlabels1 <- temp %>% filter(direction == CHIPCTRL) %>% mutate(gene = case_when((gene %in% c("MALAT1", "NEAT1", "EGR1")) ~ gene, TRUE ~ ""))
forlabels2 <- temp %>% filter(direction == CHIPCTRL) %>% mutate(gene = case_when((neg_log10_adj_pval > 100 & !(gene %in% c("MALAT1", "NEAT1", "EGR1"))) ~ gene, TRUE ~ ""))
# avg_log2FC_COVpos > 0 & 
p1 <- ggplot(temp %>% filter(direction == CHIPCTRL), aes(x = avg_log2FC_COVpos, y = neg_log10_adj_pval)) +
    geom_point(size=2, stroke=0.1, shape=16, color = colors.mn["CD14 Mono"]) +
    labs(x = expression(COVID^"+"~": log2(Avg.expr."~italic("TET2")^"MT"~"/ Avg.expr."~italic("TET2")^"WT"~")"), y = "-log10(Adjusted P value)") +
    geom_label_repel(data = forlabels1, aes(label=gene), segment.colour="black", arrow = arrow(length = unit(0.020, "npc")), size = 4, segment.size=0.3, max.overlaps = Inf, min.segment.length = 0, box.padding = 1, fontface = "italic") +
    geom_text_repel(data = forlabels2, aes(label=gene), size = 2.5, max.overlaps = 15, fontface = "italic") +
    ggtitle(expression("Same direction as COVID"^"-")) +
    theme_cowplot() + theme(legend.position = "none", axis.title = element_text(size = 12))
forlabels <- temp %>% filter(direction != CHIPCTRL) %>% mutate(gene = case_when((avg_log2FC_COVpos > 0.5 | neg_log10_adj_pval > 100) ~ gene, TRUE ~ ""))
p2 <- ggplot(temp %>% filter(direction != CHIPCTRL), aes(x = avg_log2FC_COVpos, y = neg_log10_adj_pval)) +
    geom_point(size=2, stroke=0.1, shape=16, color = colors.mn["CD14 Mono"]) +
    labs(x = expression(COVID^"+"~": log2(Avg.expr."~italic("TET2")^"MT"~"/ Avg.expr."~italic("TET2")^"WT"~")"), y = "-log10(Adjusted P value)") +
    geom_text_repel(data = forlabels, aes(label=gene), size = 2.5, max.overlaps = Inf, fontface = "italic") +
    ggtitle(expression("Opposing direction as COVID"^"-")) +
    theme_cowplot() + theme(legend.position = "none", axis.title = element_text(size = 12))
p1 / p2
dev.off()

TET2C_COVID_CD14Mono <- degsr_TET2C_COVID %>% filter(celltype == "CD14 Mono") %>% select(gene, direction) %>% unite("gene_dir", gene:direction, sep = "_") %>% unlist() %>% unname()
TET2_CTRLY_CD14Mono <- degsr_TET2_CTRLY %>% filter(celltype == "CD14 Mono") %>% select(gene, direction) %>% unite("gene_dir", gene:direction, sep = "_") %>% unlist() %>% unname()

length(intersect(TET2_CTRLY_CD14Mono, TET2C_COVID_CD14Mono))
length(setdiff(TET2_CTRLY_CD14Mono, TET2C_COVID_CD14Mono))
length(setdiff(TET2C_COVID_CD14Mono, TET2_CTRLY_CD14Mono))

library(eulerr)
venn <- euler(c("A" = length(TET2_CTRLY_CD14Mono), "B" = length(TET2C_COVID_CD14Mono), "A&B" = length(intersect(TET2_CTRLY_CD14Mono, TET2C_COVID_CD14Mono))))
pdf(file = "diamond_venn.pdf", width = 6, height = 6)
plot(venn, quantities = TRUE, fill = "transparent")
dev.off()

##################### Figure 1 - heatmaps #######################

data1 <- read.table("/CC_multiomics/scRNA_MULT/data1.txt", sep = "\t") %>% 
    set_colnames(c("ID", "CHIP", "MRN", "Age", "DOB", "Sex", "Gene1", "Gene2", "Gene3", "Gene4", "HB", "WBC", "ABC", "AMC", "ALC", "PLT", "CRP", "Ferritin", "Oxygen", "CRS", "AKI", "ALI", "ARDS", "MODS", "Thrombosis")) %>%
    mutate(TET2 = ("TET2" == Gene1 | "TET2" == Gene2 | "TET2" == Gene3 | "TET2" == Gene4)) %>%
    mutate(DNMT3A = ("DNMT3A" == Gene1 | "DNMT3A" == Gene2 | "DNMT3A" == Gene3 | "DNMT3A" == Gene4)) %>%
    mutate(SF3B1 = ("SF3B1" == Gene1 | "SF3B1" == Gene2 | "SF3B1" == Gene3 | "SF3B1" == Gene4)) %>%
    mutate(ASXL1 = ("ASXL1" == Gene1 | "ASXL1" == Gene2 | "ASXL1" == Gene3 | "ASXL1" == Gene4)) %>%
    mutate(MPL = ("MPL" == Gene1 | "MPL" == Gene2 | "MPL" == Gene3 | "MPL" == Gene4)) %>%
    mutate(BRCC3 = ("BRCC3" == Gene1 | "BRCC3" == Gene2 | "BRCC3" == Gene3 | "BRCC3" == Gene4)) %>%
    mutate(CBL = ("CBL" == Gene1 | "CBL" == Gene2 | "CBL" == Gene3 | "CBL" == Gene4)) %>%
    mutate(CREBBP = ("CREBBP" == Gene1 | "CREBBP" == Gene2 | "CREBBP" == Gene3 | "CREBBP" == Gene4)) %>%
    mutate(EZH2 = ("EZH2" == Gene1 | "EZH2" == Gene2 | "EZH2" == Gene3 | "EZH2" == Gene4)) %>%
    mutate(IDH2 = ("IDH2" == Gene1 | "IDH2" == Gene2 | "IDH2" == Gene3 | "IDH2" == Gene4)) %>%
    mutate(JAK2 = ("JAK2" == Gene1 | "JAK2" == Gene2 | "JAK2" == Gene3 | "JAK2" == Gene4)) %>%
    mutate(KIT = ("KIT" == Gene1 | "KIT" == Gene2 | "KIT" == Gene3 | "KIT" == Gene4)) %>%
    mutate(KRAS = ("KRAS" == Gene1 | "KRAS" == Gene2 | "KRAS" == Gene3 | "KRAS" == Gene4)) %>%
    mutate(PPM1D = ("PPM1D" == Gene1 | "PPM1D" == Gene2 | "PPM1D" == Gene3 | "PPM1D" == Gene4)) %>%
    mutate(SETD2 = ("SETD2" == Gene1 | "SETD2" == Gene2 | "SETD2" == Gene3 | "SETD2" == Gene4)) %>%
    mutate(STAT3 = ("STAT3" == Gene1 | "STAT3" == Gene2 | "STAT3" == Gene3 | "STAT3" == Gene4)) %>%
    mutate(SUZ12 = ("SUZ12" == Gene1 | "SUZ12" == Gene2 | "SUZ12" == Gene3 | "SUZ12" == Gene4)) %>%
    mutate(TP53 = ("TP53" == Gene1 | "TP53" == Gene2 | "TP53" == Gene3 | "TP53" == Gene4)) %>%
    filter(MPL == FALSE)
data2 <- read.table("/CC_multiomics/scRNA_MULT/data2.txt", sep = "\t") %>% select(V1:V14) %>%
    set_colnames(c("TNF", "IL6", "IFNb", "IL10", "MCP1", "IL1b", "IFNg", "MIP1a", "GMCSF", "IL2Ra", "IFNa", "IL18", "RunDate", "ID"))
data3 <- read.table("/CC_multiomics/scRNA_MULT/data3.txt", sep = "\t") %>% select(V1:V4) %>%
    set_colnames(c("ID", "CRS.MAX", "CRS.LBL", "CRS.WHO"))
data_in <- data1 %>% left_join(data2) %>% left_join(data3)

library(scico)
factorcols <- c("CHIP", "Sex", "Oxygen", "CRS", "CRS.WHO", "ASXL1", "BRCC3", "CBL", "CREBBP", "DNMT3A", "EZH2", "IDH2", "JAK2", "KIT", "KRAS", "MPL", "PPM1D", "SETD2", "SF3B1", "STAT3", "SUZ12", "TET2", "TP53", "AKI", "ALI", "ARDS", "MODS", "Thrombosis")
hm_data_in <- data_in %>% 
    mutate(ID = factor(ID)) %>%
    mutate(Sex = factor(Sex, levels = c("Female", "Male"))) %>%
    mutate_at(factorcols, as.factor) %>%
    mutate(CRS.WHO = factor(CRS.WHO, levels = 1:8)) %>%
    mutate(CRS.MAX = factor(CRS.MAX, levels = 1:5)) %>%
    mutate(CRS.LBL = case_when(CRS.LBL %in% c("Death", "Life") ~ "Life-threatening/fatal", TRUE ~ CRS.LBL)) %>%
    mutate(CRS.LBL = factor(CRS.LBL, levels = rev(c("Life-threatening/fatal", "Severe", "Moderate", "Mild")))) %>%
    group_by(CHIP, CRS.WHO) %>% 
    arrange(desc(CHIP), desc(CRS.WHO))

hm_row <- function(variable, colorset, show.legend, italicize, aspect){
    p <- ggplot(hm_data, aes(x = ID, y = variable)) +
        geom_tile(aes_string(fill = variable), color = "white") +
        scale_fill_manual(values = colorset, na.value = "grey80") +
        scale_x_discrete(expand = c(0, 0), position = "top")+
        scale_y_discrete(expand = c(0, 0)) +
        theme_cowplot() +
        theme(aspect.ratio=aspect, axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=5), axis.ticks = element_blank(), axis.ticks.length = unit(0, "pt"), axis.line = element_blank(), plot.margin = margin(0, 0, 0, 0, unit = "pt"), legend.key.size = unit(5,"points"), legend.title = element_text(size=6), legend.text = element_text(size=5))
    if(!show.legend){
        p <- p + theme(legend.position = "none")
    }
    if(italicize){
        p <- p + theme(axis.text.y = element_text(face = "italic"))
    }
    return(p)
}

hm_row_cont <- function(variable, palette, direction, show.legend, aspect){
    p <- ggplot(hm_data, aes(x = ID, y = variable)) +
        geom_tile(aes_string(fill = variable), color = "white") +
        scale_fill_scico(palette = palette, direction = direction, na.value = "grey80", limits = c(0, 1)) +
        scale_x_discrete(expand = c(0, 0), position = "top")+
        scale_y_discrete(expand = c(0, 0)) +
        theme_cowplot() +
        theme(aspect.ratio=aspect, axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=5), axis.ticks = element_blank(), axis.ticks.length = unit(0, "pt"), axis.line = element_blank(), plot.margin = margin(0, 0, 0, 0, unit = "pt"), legend.key.size = unit(5,"points"), legend.title = element_text(size=6), legend.text = element_text(size=5))
    if(!show.legend){
        return(p + theme(legend.position = "none"))
    }
    return(p)
}

#########

hm_data <- hm_data_in %>% mutate(ID = factor(ID, levels = hm_data_in$ID))
hm_data <- hm_data %>% filter(CHIP == "yes")
hm_data <- hm_data %>% ungroup() %>% mutate(across(c("CRP", "Ferritin", "TNF", "IL6", "IFNb", "IL10", "MCP1", "IL1b", "GMCSF", "HB", "WBC", "ABC", "AMC", "ALC", "PLT"), function(x) {(x-min(x, na.rm = TRUE))/diff(range(x, na.rm = TRUE))}))

WHO_cols <- scico(8, palette = 'lajolla')
names(WHO_cols) <- 1:8
allcolors <- list(WHO_cols, c("white", "#009E73"), c("#332288", "#DDCC77"))
allshowlegends <- c(TRUE, FALSE, TRUE)
allitalics <- c(FALSE, FALSE, FALSE)
plots1 <- mapply(hm_row, c("CRS.WHO", "Oxygen", "Sex"), allcolors, allshowlegends, allitalics, MoreArgs = list(aspect = 0.04), SIMPLIFY=FALSE)

allpalettes <- c("imola", "imola", "acton", "acton", "acton", "acton", "acton", "acton", "acton", "tokyo", "tokyo", "tokyo", "tokyo", "tokyo", "tokyo")
alldirections <- c(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)
allshowlegends <- c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
plots2 <- mapply(hm_row_cont, c("CRP", "Ferritin", "TNF", "IL6", "IFNb", "IL10", "MCP1", "IL1b", "GMCSF", "HB", "WBC", "ABC", "AMC", "ALC", "PLT"), allpalettes, alldirections, allshowlegends, MoreArgs = list(aspect = 0.04), SIMPLIFY=FALSE)

allcolors <- list(c("white", "#44AA99"), c("white", "#44AA99"), c("white", "#44AA99"), c("white", "#44AA99"))
allshowlegends <- rep(FALSE, 4)
allitalics <- rep(FALSE, 4)
plots3 <- mapply(hm_row, c("AKI", "ARDS", "MODS", "Thrombosis"), allcolors, allshowlegends, allitalics, MoreArgs = list(aspect = 0.04), SIMPLIFY=FALSE)

allcolors <- list(c("white", "#56B4E9"), c("white", "#56B4E9"), c("white", "#56B4E9"), c("white", "#56B4E9"), c("white", "#56B4E9"), c("white", "#56B4E9"), c("white", "#56B4E9"), c("white", "#56B4E9"), c("white", "#56B4E9"), c("white", "#56B4E9"), c("white", "#56B4E9"), c("white", "#56B4E9"), c("white", "#56B4E9"), c("white", "#56B4E9"), c("white", "#56B4E9"), c("white", "#56B4E9"), c("white", "#56B4E9"))
allshowlegends <- rep(FALSE, 17)
allitalics <- rep(TRUE, 17)
plots4 <- mapply(hm_row, c("DNMT3A", "TET2", "SF3B1", "ASXL1", "BRCC3", "CBL", "CREBBP", "EZH2", "IDH2", "JAK2", "KIT", "KRAS", "PPM1D", "SETD2", "STAT3", "SUZ12", "TP53"), allcolors, allshowlegends, allitalics, MoreArgs = list(aspect = 0.04), SIMPLIFY=FALSE)

pdf(file = "hm.pdf", width = 3, height = 5)
wrap_plots(c(plots1, plots4, plots2, plots3), ncol = 1) + plot_layout(guides = "collect")
dev.off()

#########

hm_data <- hm_data_in %>% mutate(ID = factor(ID, levels = hm_data_in$ID))
hm_data <- hm_data %>% filter(CHIP == "no")
hm_data <- hm_data %>% ungroup() %>% mutate(across(c("CRP", "Ferritin", "TNF", "IL6", "IFNb", "IL10", "MCP1", "IL1b", "GMCSF", "HB", "WBC", "ABC", "AMC", "ALC", "PLT"), function(x) {(x-min(x, na.rm = TRUE))/diff(range(x, na.rm = TRUE))}))

WHO_cols <- scico(8, palette = 'lajolla')
names(WHO_cols) <- 1:8
allcolors <- list(WHO_cols, c("white", "#009E73"), c("#332288", "#DDCC77"))
allshowlegends <- c(TRUE, FALSE, TRUE)
allitalics <- c(FALSE, FALSE, FALSE)
plots1 <- mapply(hm_row, c("CRS.WHO", "Oxygen", "Sex"), allcolors, allshowlegends, allitalics, MoreArgs = list(aspect = 0.025), SIMPLIFY=FALSE)

allpalettes <- c("imola", "imola", "acton", "acton", "acton", "acton", "acton", "acton", "acton", "tokyo", "tokyo", "tokyo", "tokyo", "tokyo", "tokyo")
alldirections <- c(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)
allshowlegends <- c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
plots2 <- mapply(hm_row_cont, c("CRP", "Ferritin", "TNF", "IL6", "IFNb", "IL10", "MCP1", "IL1b", "GMCSF", "HB", "WBC", "ABC", "AMC", "ALC", "PLT"), allpalettes, alldirections, allshowlegends, MoreArgs = list(aspect = 0.025), SIMPLIFY=FALSE)

allcolors <- list(c("white", "#44AA99"), c("white", "#44AA99"), c("white", "#44AA99"), c("white", "#44AA99"))
allshowlegends <- rep(FALSE, 4)
allitalics <- rep(FALSE, 4)
plots3 <- mapply(hm_row, c("AKI", "ARDS", "MODS", "Thrombosis"), allcolors, allshowlegends, allitalics, MoreArgs = list(aspect = 0.025), SIMPLIFY=FALSE)

pdf(file = "hm_supp.pdf", width = 4, height = 2)
wrap_plots(c(plots1, plots2, plots3), ncol = 1) + plot_layout(guides = "collect")
dev.off()

################# VAF Violins ############################

vafs <- read.table("/CC_multiomics/scRNA_MULT/vaf.txt", sep = "\t") %>% 
    set_colnames(c("gene", "vaf")) %>%
    filter(gene %in% c("DNMT3A", "TET2", "SF3B1", "ASXL1", "TP53")) %>%
    mutate(gene = factor(gene, levels = c("DNMT3A", "TET2", "SF3B1", "ASXL1", "TP53")))

temp <- vafs %>% group_by(gene) %>% summarise(n = n())

pdf(file='vafviolins.pdf', width=5.5, height=3.5)
ggplot(vafs, aes(x = gene, y = vaf)) +
    geom_violin(scale = "width") +
    geom_boxplot(width=0.1, outlier.shape = NA) +
    geom_text(data = temp, mapping = aes(x = gene, label = paste0("n = ", n)), y = 0.55, color = "black", size = 5) +
    scale_y_continuous(limits = c(0, 0.55)) +
    ylab("VAF") +
    theme_cowplot() + theme(axis.title.x = element_blank())
dev.off()

library(ggokabeito)
temp <- hm_data_in %>% filter(!(CRS.LBL %in% c(NA)))
pdf(file='patientprop.pdf', width=5.5, height=3.5)
ggplot(temp, aes(x = CHIP, fill = CRS.LBL)) +
    geom_bar(position = "fill") +
    scale_fill_okabe_ito(order = c(4, 2, 3, 1)) +
    ylab("Proportion of patients") +
    theme_cowplot()
dev.off()

###########################  TET2 ChIP-seq (signal) from Solary Lab vs. EGR1 ChIP-seq (peaks) in Monocytes from Alessandro Lab ################

# cd /FINAL/
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
# CrossMap.py bed hg19ToHg38.over.chain.gz 1083_PEAKS_GAINING_EGR1_AT_DAY3-5.bed 1083_PEAKS_GAINING_EGR1_AT_DAY3-5_hg38.bed

# awk 'BEGIN {OFS="\t"}; {if ($2-500 < 1) {start = 1} else {start = $2-500}; print $1,start,$3+500,$4;}' 1083_PEAKS_GAINING_EGR1_AT_DAY3-5_hg38.bed > EGR1_Mono_hg38.bed

cd /FINAL/
Mono=/FINAL/EGR1_Mono_hg38.bed
TET2=/FINAL/81224_TET2_blood_mono_hg38.bigwig

computeMatrix reference-point -S ${TET2} \
                              -R ${Mono} \
                              --referencePoint center \
                              -b 1000 -a 1000 -p 16 \
                              --missingDataAsZero \
                              --outFileNameMatrix TET2inMono \
                              --outFileName TET2inMono.tab.gz \
                              --sortRegions descend \
                              --sortUsing mean
plotHeatmap -m TET2inMono.tab.gz --colorList 'white,blue' --sortRegions keep -out TET2inEGR1gainpeaks.pdf --heatmapHeight 12 --heatmapWidth 8

########################## Comparing differentially expressed genes in our study with differentially expressed genes in Alessandro study #########################

EGR1_targets <- read.table("/FINAL/aaz8836_table_s4.txt", sep = "\t", header = TRUE) %>%
    separate(col = "Gene", into = c("GeneID", NA), sep = "\\.")
genes_map <- read.table("/cellranger_output/scRNA_cr_6.0.1/C8/outs/filtered_feature_bc_matrix/features.tsv") %>%
    select(V1, V2) %>% 
    set_colnames(c("GeneID", "gene"))

Mono <- read.table("/FINAL/EGR1_Mono_hg38.bed",col.names = c("chr", "start", "end", "peakname"))
Mono %<>% select(chr, start, end) %>% unite(everything(), sep = "-", col = "peak") %>% unlist() %>% unname()
Mono_targets <- ClosestFeature(MULT, regions = Mono)
Mono_targets %<>% filter(distance < 5000)

degsr <- read.table(file = "/FINAL/scRNA/degsr_TET2C_COVID.tsv", sep = "\t", skip = 1) %>% set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"))

theme_set(theme_cowplot())
temp1 <- degsr %>% 
    filter(p_val_adj < 0.05) %>% 
    mutate(celltype = factor(celltype, levels = names(colors.mn))) %>% 
    filter(celltype %in% c("CD14 Mono")) %>%
    mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_adj_pval = -log10(p_val_adj))
pdf(file = "GenesNearEGR1Peaks.pdf", width = 8, height = 8)
p1 <- ggplot(temp1 %>% filter(gene %in% Mono_targets$gene_name), aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
    geom_point(aes(color = celltype), size = 1.5, stroke=0.1, shape = 16, show.legend = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_label_repel(aes(label=gene), size = 3, max.overlaps = Inf, show.legend = FALSE, fontface = "italic") +
    coord_cartesian(xlim = c(-1.5, 1.5)) +
    scale_color_manual(values = colors.mn, aesthetics = c("color", "segment.color")) +
    labs(x = "log2(Fold change of average expression)", y = "-log10(Adjusted P value)") +
    guides(color = guide_legend(override.aes = list(size=5))) + ggtitle("Genes near EGR1 peaks in monocytes")
plot(p1)
dev.off()

