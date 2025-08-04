setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
source(file.path(general_scripts_folder, "RNA", "create_seurat_from_cr_h5.R"))
source(file.path(general_scripts_folder, "genome_wide", "gtf_filter_gene_names.R"))
library(ggplot2)
library(Seurat)
library(qs)
library(dplyr)
library(reshape2)

cr_path <- file.path("..", "preprocessing", "01_cellranger")
objects_path <- file.path(project_folder, "objects", "R", "seurat")
dir.create(objects_path, showWarnings = FALSE, recursive = TRUE)
output_path <- file.path(project_folder, "output")
qc_path <- file.path(output_path, "qc")
dir.create(qc_path, showWarnings = FALSE, recursive = TRUE)

metadata_path <- file.path(project_folder, "metadata", "metadata.csv")
metadata <- read.csv(metadata_path, comment.char = "#", row.names = 2)

gtf_path <- file.path(genomes_folder, "CR_reference", "refdata-gex-GRCm39-2024-A", "genes", "genes.gtf")
gtf <- gtf_filter_gene_names(gtf_path)
pc_genes <- unique(setdiff(gtf %>% filter(gene_type == "protein_coding") %>% pull(gene_name), ""))

so_list <- lapply(rownames(metadata), function(sample_barcode) {
    sample_name <- metadata[sample_barcode, "Sample.name"]
    sample_condition <- metadata[sample_barcode, "Condition"]

    h5_path <- file.path(cr_path, sample_barcode, "outs", "raw_feature_bc_matrix.h5")

    temp_so <- create_so_default(
        h5_path = h5_path,
        sample_name = sample_name
    )

    temp_so$condition <- sample_condition

    return(temp_so)
})

so <- JoinLayers(merge(
    x = so_list[[1]],
    y = so_list[2:length(so_list)],
    add.cell.ids = metadata$Sample.name
))
so$sample_name <- factor(so$sample_name)
so$condition <- factor(so$condition)

mt_genes <- sort(grep("^MT-", rownames(so), value = FALSE, ignore.case = TRUE))
rp_genes <- sort(c(
    grep("^RP[SL]\\d+", rownames(so), value = FALSE, ignore.case = TRUE),
    grep("^MRP[SL]\\d+", rownames(so), value = FALSE, ignore.case = TRUE)
))

so$percent_mt <- PercentageFeatureSet(so, features = rownames(so)[mt_genes])
so$percent_rp <- PercentageFeatureSet(so, features = rownames(so)[rp_genes])
so <- subset(so, features = -c(mt_genes, rp_genes))

mtd_rna <- c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_rp")

pdf(file.path(qc_path, "unfiltered_violin.pdf"), width = 15)
print(VlnPlot(so, group.by = "sample_name", features = mtd_rna, ncol = 4, pt.size = 0, log = FALSE))
print(VlnPlot(so, group.by = "sample_name", features = mtd_rna, ncol = 4, pt.size = 0, log = TRUE))
dev.off()

qsave(so, file.path(objects_path, "raw_aggregate.qs"), nthreads = ncores)
# so <- qread(file.path(objects_path, "raw_aggregate.qs"), nthreads = ncores)

filtered_so <- subset(so, nFeature_RNA > 1000 & nFeature_RNA < 5000 & nCount_RNA > 1000 & nCount_RNA < 1.5e4 & percent_mt < 0.2 & percent_rp < 1.5)

pdf(file.path(qc_path, "filtered_violin.pdf"), width = 15)
print(VlnPlot(filtered_so, group.by = "sample_name", features = mtd_rna, ncol = 4, pt.size = 0, log = FALSE))
print(VlnPlot(filtered_so, group.by = "sample_name", features = mtd_rna, ncol = 4, pt.size = 0, log = TRUE))

percent_pc_counts_df <- data.frame(sample_name = unique(filtered_so$sample_name))
total_pc_counts <- rep(0, nrow(percent_pc_counts_df))
total_counts <- rep(0, nrow(percent_pc_counts_df))

for (i in seq_len(nrow(percent_pc_counts_df))) {
    sample_val <- percent_pc_counts_df$sample_name[i]
    sample_so <- subset(filtered_so, sample_name == sample_val)

    total_counts[i] <- sum(sample_so$nCount_RNA)
    total_pc_counts[i] <- sum(sample_so@assays$RNA@layers$counts[rownames(sample_so) %in% pc_genes, ])
}

# percent_pc_counts_df$total_counts <- total_counts
# percent_pc_counts_df$total_pc_counts <- total_pc_counts
percent_pc_counts_df$percent_pc_counts <- total_pc_counts / total_counts
percent_pc_counts_df$percent_npc_counts <- 1 - percent_pc_counts_df$percent_pc_counts
percent_pc_counts_df <-  melt(percent_pc_counts_df)
# reorder the levels of variable factor
percent_pc_counts_df$variable <- factor(percent_pc_counts_df$variable, levels = c("percent_npc_counts", "percent_pc_counts"))

# stacked barplot with percent pc npc for each sample
ggplot(percent_pc_counts_df, aes(x = sample_name, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Percent protein coding and non-protein coding counts") +
    xlab("Sample name") +
    ylab("Percent") +
    scale_fill_manual(values = c("percent_pc_counts" = "black", "percent_npc_counts" = "gray"))
dev.off()

pdf(file.path(qc_path, "filtered_scatter.pdf"))

for (x in seq_len(length(mtd_rna) - 1)) {
    for (y in (x+1):length(mtd_rna)) {
        mtd1 <- mtd_rna[x]
        mtd2 <- mtd_rna[y]
        df <- filtered_so@meta.data[, c(mtd1, mtd2, "sample_name")]
        print(ggplot(df, aes_string(x = mtd1, y = mtd2)) + geom_point() + theme_classic() + ggtitle("aggregate"))
    }
}

for (sample_val in unique(filtered_so$sample_name)) {
    sample_so <- subset(filtered_so, sample_name == sample_val)
    for (x in seq_len(length(mtd_rna) - 1)) {
        for (y in (x+1):length(mtd_rna)) {
            mtd1 <- mtd_rna[x]
            mtd2 <- mtd_rna[y]
            df <- sample_so@meta.data[, c(mtd1, mtd2, "sample_name")]
            print(ggplot(df, aes_string(x = mtd1, y = mtd2)) + geom_point() + theme_classic() + ggtitle(sample_val))
        }
    }
}
dev.off()

qsave(filtered_so, file.path(objects_path, "aggregate.qs"), nthreads = ncores)
