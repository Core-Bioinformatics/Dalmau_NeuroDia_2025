setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(qs)
library(ggplot2)
library(rhdf5)
library(dplyr)
library(AnnotationDbi)
library(EnsDb.Mmusculus.v79)
library(confintr)
library(ComplexHeatmap)
library(patchwork)
library(reshape2)
library(ggfittext)

object_path <- file.path(project_folder, "objects", "R")
figure_path <- file.path(project_folder, "output", "figures", "aggregate")
shiny_path <- file.path(project_folder, "shiny_apps", "clustassess")

if (!dir.exists(figure_path)) {
    dir.create(figure_path, recursive = TRUE)
}

aggr_so <- qread(file.path(object_path, "seurat", "aggregate_sct.qs"), nthreads = ncores) 

default_umap_ggplot <- function(umap_df,
                                mtd_value,
                                mtd_name = "",
                                colour_palette = FALSE,
                                facet_plot = FALSE,
                                is_small_dts = FALSE,
                                show_legend = TRUE,
                                label_size = 0) {
    is_continuous <- length(colour_palette) == 1 && is.logical(colour_palette) && !colour_palette
    umap_df <- as.data.frame(umap_df)
    colnames(umap_df) <- c("UMAP1", "UMAP2")
    min_x <- min(umap_df$UMAP1)
    max_x <- max(umap_df$UMAP1)
    min_y <- min(umap_df$UMAP2)
    max_y <- max(umap_df$UMAP2)

    pt_size <- ifelse(is_small_dts, 1, 0.5)

    umap_df$colour_val <- mtd_value

    if (is_continuous) {
        umap_df <- umap_df %>% arrange(.data$colour_val)
    } else {
        if (is.character(umap_df$colour_val)) {
            unique_vals <- unique(umap_df$colour_val)
            unique_vals <- stringr::str_sort(unique_vals, numeric = TRUE)
        } else {
            unique_vals <- levels(umap_df$colour_val)
        }
    }

    if (!is_continuous && is.numeric(facet_plot) && facet_plot > 0) {
        return(wrap_plots(
            lapply(unique_vals, function(unique_val) {
                temp_val <- as.character(umap_df$colour_val)
                temp_val[temp_val != unique_val] <- NA
                new_colour_palette <- colour_palette[unique_val]
                default_umap_ggplot(
                    umap_df[, 1:2],
                    temp_val,
                    mtd_name = paste0(mtd_name, " - ", unique_val),
                    colour_palette = new_colour_palette,
                    facet_plot = FALSE,
                    show_legend = FALSE,
                    is_small_dts = is_small_dts
                )
            }),
            ncol = facet_plot
        ))
    }

    if (!is_continuous) {
        umap_df <- rbind(
            umap_df[is.na(umap_df$colour_val), ],
            umap_df[!is.na(umap_df$colour_val), ]
        )
    }


    gplot_obj <- ggplot(umap_df, aes(x = .data$UMAP1, y = .data$UMAP2, colour = .data$colour_val)) +
        geom_point(size = pt_size) +
        theme_classic() +
        scale_x_continuous(breaks = c(min_x + (max_x - min_x) / 7)) +
        scale_y_continuous(breaks = c(min_y + (max_y - min_y) / 7)) +
        guides(
            x = guide_axis(cap = "upper"),
            y = guide_axis(cap = "upper")
        ) +
        ggtitle(mtd_name) +
        theme(
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(arrow = arrow(type = "closed", length = unit(0.3, "cm"))),
            axis.title.x = element_text(hjust = 0.03),
            axis.title.y = element_text(hjust = 0.03),
            legend.title = element_text(size = 0),
            plot.title = element_text(size = 25),
            axis.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            aspect.ratio = 1,
            legend.position = "bottom",
            legend.direction = "horizontal"
        )

    if (!show_legend) {
        gplot_obj <- gplot_obj + theme(legend.position = "none")
    }

    if (is_continuous) {
        return(gplot_obj + scale_color_viridis_c() + theme(
            legend.key.width = unit(3, "cm")
        ))
    }

    gplot_obj <- gplot_obj +
        scale_colour_manual(values = colour_palette, na.value = "gray85") +
        guides(
            colour = guide_legend(
            override.aes = list(
                size = pt_size * 10,
                shape = 15
            ),
            byrow = TRUE
        ))

    if (label_size == 0) {
        return(gplot_obj)
    }

    medians_values <- t(sapply(seq_along(unique_vals), function(i) {
        mask <- which(umap_df$colour_val == unique_vals[i])
        emb <- umap_df[mask, 1:2]
        if (is.null(emb) || nrow(emb) == 1) {
            return(emb)
        }

        gm_res <- Gmedian::Gmedian(emb)
        return(gm_res)
    }))
    medians_values <- data.frame(medians_values)
    colnames(medians_values) <- c("median_x", "median_y")
    medians_values$labels <- unique_vals

    return(gplot_obj + ggrepel::geom_text_repel(
        data = medians_values,
        ggplot2::aes(x = .data$median_x, y = .data$median_y, label = .data$labels),
        size = label_size,
        color = "white",
        bg.color = "black",
        bg.r = 0.15,
        nudge_x = 0.15,
        nudge_y = 0.15
    ))
}

# sample renaming
sample_renaming <- list(
    "240911_1" = "WT NCD 1",
    "240911_2" = "WT NCD 2",
    "240911_3" = "WT NCD 3",
    "240911_4" = "Sucnr1KO NCD 1",
    "240911_5" = "Sucnr1KO NCD 2",
    "240911_6" = "Sucnr1KO NCD 3",
    "240911_7" = "WT HFD 1",
    "240911_8" = "WT HFD 2",
    "240911_9" = "WT HFD 3",
    "240911_10" = "Sucnr1KO HFD 1",
    "240911_11" = "Sucnr1KO HFD 2",
    "240911_12" = "Sucnr1KO HFD 3"
)

new_sample_name <- as.character(aggr_so$sample_name)
for (i in names(sample_renaming)) {
    new_sample_name[new_sample_name == i] <- sample_renaming[[i]]
}
aggr_so$sample_name <- factor(new_sample_name, levels = unlist(sample_renaming))

###### VIOLIN PLOTS - QC: nfeature ncount mt rp ######
Idents(aggr_so) <- "sample_name"
pal <- h5read(file.path(shiny_path, "aggregate", "stability.h5"), paste0("/colors/", nlevels(aggr_so$sample_name)))
names(pal) <- levels(aggr_so$sample_name)

pdf(file.path(figure_path, "qc_violin.pdf"), width = 25, height = 10)
VlnPlot(aggr_so, features = c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_rp"), ncol = 4, pt.size = 0, cols = pal) 
dev.off()



###### UMAPs - QC: nfeature ncount mt rp cellcycle ######
for (mtd_names in c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_rp")) {
    pdf(file.path(figure_path, paste0("umap_", mtd_names, ".pdf")), width = 10, height = 10)
    print(default_umap_ggplot(
        aggr_so@reductions$umap@cell.embeddings,
        aggr_so[[mtd_names]][,1],
        mtd_name = mtd_names,
        colour_palette = FALSE
    ))
    dev.off()
}

pal <- h5read(file.path(shiny_path, "aggregate", "stability.h5"), paste0("/colors/", 3))
names(pal) <- c("G1", "S", "G2M")

pdf(file.path(figure_path, "umap_cell_cycle.pdf"), width = 10, height = 10)
print(default_umap_ggplot(
    aggr_so@reductions$umap@cell.embeddings,
    aggr_so$Phase,
    mtd_name = "Cell cycle phase",
    colour_palette = pal
))
dev.off()

pdf(file.path(figure_path, "umap_cell_cycle_facet.pdf"), width = 30, height = 10)
print(default_umap_ggplot(
    aggr_so@reductions$umap@cell.embeddings,
    aggr_so$Phase,
    mtd_name = "Cell cycle phase",
    colour_palette = pal,
    facet_plot = 3
))
dev.off()

###### UMAPs - distributions of experimental groups, cell types, clusters ######
# experimental groups
levels(aggr_so$genotype) <- c("Sucnr1KO global", "WT")
aggr_so$diet_genotype <- paste(as.character(aggr_so$genotype), as.character(aggr_so$diet))

pal <- c(
    "WT NCD" = "#B8D6B6",
    "Sucnr1KO global NCD" = "#095256",
    "WT HFD" = "#D6A99A",
    "Sucnr1KO global HFD" = "#880D1E"
)
aggr_so$diet_genotype <- factor(aggr_so$diet_genotype, levels = names(pal))

pdf(file.path(figure_path, "umap_experimental_groups.pdf"), width = 10, height = 10)
print(default_umap_ggplot(
    aggr_so@reductions$umap@cell.embeddings,
    aggr_so$diet_genotype,
    mtd_name = "Experimental groups",
    colour_palette = pal
))
dev.off()

pdf(file.path(figure_path, "umap_experimental_groups_facet.pdf"), width = 40, height = 10)
print(default_umap_ggplot(
    aggr_so@reductions$umap@cell.embeddings,
    aggr_so$diet_genotype,
    mtd_name = "Experimental groups",
    colour_palette = pal,
    facet_plot = 4
))
dev.off()
aggr_so$diet_genotype <- factor(aggr_so$diet_genotype, levels = rev(names(pal)))

# cell types
levels(aggr_so$aggregate_annotation) <- c(
    "Astrocytes",
    "Microglia",
    "Multicellular niche",
    "Neurons",
    "Oligodendrocytes",
    "OPSc"
)
pal <- c(
    "OPSc" = "#BD7293",
    "Oligodendrocytes" = "#C9AC81",
    "Microglia" = "#EB573B",
    "Astrocytes" = "#88C879",
    "Multicellular niche" = "#706EBE",
    "Neurons" = "#7DB2C4"
)

pdf(file.path(figure_path, "umap_cell_types.pdf"), width = 10, height = 10)
print(default_umap_ggplot(
    aggr_so@reductions$umap@cell.embeddings,
    aggr_so$aggregate_annotation,
    mtd_name = "Cell types",
    colour_palette = pal
))
dev.off()

pdf(file.path(figure_path, "umap_cell_types_facet.pdf"), width = 30, height = 20)
print(default_umap_ggplot(
    aggr_so@reductions$umap@cell.embeddings,
    aggr_so$aggregate_annotation,
    mtd_name = "Cell types",
    colour_palette = pal,
    facet_plot = 3
))
dev.off()

# clusters
for (k in c(20, 21)) {
    pal <- h5read(file.path(shiny_path, "aggregate", "stability.h5"), paste0("/colors/", k))
    names(pal) <- as.character(seq_len(k))
    cl <- aggr_so[[paste0("stable_clusters_", k)]][, 1]
    cl <- factor(cl, levels = seq_len(k))

    pdf(file.path(figure_path, paste0("umap_clusters_", k, ".pdf")), width = 10, height = 10)
    print(default_umap_ggplot(
        aggr_so@reductions$umap@cell.embeddings,
        cl,
        mtd_name = paste0("Stable ", k, " Clusters"),
        colour_palette = pal,
        label_size = 10
    ))
    dev.off()

    pdf(file.path(figure_path, paste0("umap_clusters_facet_", k, ".pdf")), width = 30, height = 70)
    print(default_umap_ggplot(
        aggr_so@reductions$umap@cell.embeddings,
        cl,
        mtd_name = paste0("Stable ", k, " Clusters"),
        colour_palette = pal,
        facet_plot = 3
    ))
    dev.off()

}

# genotype
pal <- c(
    "WT" = "#2690b4",
    "Sucnr1KO global" = "#d22187"
)

cl_vec <- aggr_so$genotype
cl_vec[aggr_so$diet != "NCD"] <- NA

pdf(file.path(figure_path, "umap_genotype_on_ncd.pdf"), width = 10, height = 10)
default_umap_ggplot(
    aggr_so@reductions$umap@cell.embeddings,
    cl_vec,
    mtd_name = "NCD",
    colour_palette = pal
)
dev.off()

cl_vec <- aggr_so$genotype
cl_vec[aggr_so$diet == "NCD"] <- NA

pdf(file.path(figure_path, "umap_genotype_on_hfd.pdf"), width = 10, height = 10)
default_umap_ggplot(
    aggr_so@reductions$umap@cell.embeddings,
    cl_vec,
    mtd_name = "HFD",
    colour_palette = pal
)
dev.off()


###### GENE UMAP ######
gene_list <- list(
    "Astrocytes" = list(
        genes = c("Slc1a2", "Aldoc", "Gja1", "Sparcl1", "Sparc", "Sox9", "Gfap", "Gpc5"),
        thresh = 0,
        not_expr = -1
    ),
    "Microglia" = list(
        genes = c("P2ry12", "Cx3cr1", "Siglech", "Hexb", "Sall1", "Trem2", "Aif1", "Ctss", "Tmem119"),
        thresh = 0,
        not_expr = -1
    ),
    "Oligodendrocytes" = list(
        genes = c("Mbp", "Bcas1", "Plp1", "Cnp", "Mog", "Mag"),
        thresh = 0,
        not_expr = -1
    ),
    "OPSc" = list(
        genes = c("Pdgfra", "Cspg4", "Olig2"),
        thresh = 0,
        not_expr = -1
    ),
    "Fibroblasts/ Endothelial Cells/Tanycytes" = list(
        genes = c("Vim", "Flt1"),
        thresh = 0,
        not_expr = -1
    ),
    "Tanycytes" = list(
        genes = c("Rax", "Vim", "Col23a1"),
        thresh = 0,
        not_expr = -1
    ),
    "Endothelial" = list(
        genes = c("Ocln", "Abcb1a", "Lef1", "Flt1", "Col4a6"),
        thresh = 0,
        not_expr = -1
    ),
    "Fibroblasts" = list(
        genes = c("Col1a1", "Vim", "Pdgfrb", "Vtn"),
        thresh = 0,
        not_expr = -1
    ),
    "Epithelial" = list(
        genes = c("Ccdc170", "Col23a1", "Nek5", "Col4a6", "Cfap52", "Dnah6", "Ttc21a", "Pdzrn3", "Ak9"),
        thresh = 0,
        not_expr = -1
    ),
    "Pericytes" = list(
        genes = c("Tbx18"),
        thresh = 0,
        not_expr = -1
    ),
    "Neurons" = list(
        genes = c("Slc17a6", "Slc32a1", "Rbfox3", "Dcx"),
        thresh = 0.6,
        not_expr = 2
    ),
    "Multicellular niche" = list(
        genes = c("Ccdc170", "Col23a1", "Col4a6", "Vim", "Flt1", "Rax"),
        thresh = 0,
        not_expr = -1
    )
)
expr_mat <- GetAssayData(aggr_so, layer = "data")

for (i in names(gene_list)) {
    print(setdiff(gene_list[[i]]$genes, rownames(expr_mat)))
    gene_list[[i]]$genes <- intersect(gene_list[[i]]$genes, rownames(expr_mat))
}

color_values <- function(n) {
    grDevices::colorRampPalette(c("grey85", paletteer::paletteer_d("RColorBrewer::OrRd")))(n)
}
entire_gene_list <- unlist(lapply(gene_list, function(x) x$genes))
expr_mat <- expr_mat[entire_gene_list, , drop = FALSE]
gc()

gene_dir <- file.path(figure_path, "gene_umap")
if (!dir.exists(gene_dir)) {
    dir.create(gene_dir, recursive = TRUE)
}

for (g in entire_gene_list) {
    pdf(file.path(gene_dir, paste0(g, ".pdf")), width = 10, height = 10)
    print(default_umap_ggplot(
        aggr_so@reductions$umap@cell.embeddings,
        expr_mat[g, ],
        mtd_name = g,
    ) + scale_colour_gradientn(colours = color_values(100)))
    dev.off()
}

###### AGGREGATE GENE UMAP ######
gene_dir <- file.path(figure_path, "aggregate_gene_umap")
if (!dir.exists(gene_dir)) {
    dir.create(gene_dir, recursive = TRUE)
}

for (ctype in names(gene_list)) {
    g <- gene_list[[ctype]]$genes
    thresh <- gene_list[[ctype]]$thresh
    not_expr <- gene_list[[ctype]]$not_expr
    if (not_expr == -1) {
        not_expr <- length(g) - 1
    }

    ctype <- stringr::str_replace_all(ctype, "/", ",")
    avg_expr <- colMeans(expr_mat[g, , drop = FALSE])
    print(sum(avg_expr > thresh))
    mask <- colSums(expr_mat[g, , drop = FALSE] > thresh) >= (length(g) - not_expr)
    print(sum(mask))
    avg_expr[!mask] <- 0 

    pdf(file.path(gene_dir, paste0(ctype, ".pdf")), width = 10, height = 10)
    print(default_umap_ggplot(
        aggr_so@reductions$umap@cell.embeddings,
        avg_expr,
        mtd_name = paste0(ctype, "\n", paste(g[seq_len(min(6, length(g)))], collapse = ", ")),
    ) + scale_colour_gradientn(colours = color_values(100)) +
    labs(title = ctype, subtitle = paste(g[seq_len(min(6, length(g)))], collapse = ", ")) +
    theme(
        plot.title = element_text(size = 25, face = "bold"),
        plot.subtitle = element_text(size = 20, face = "italic"),
    )
    )
    dev.off()
}


###### STACKED BAR PLOTS - clusters + annotation, clusters + experimental groups ######
# annotation + exp groups
pal <- list(
    "exp_groups" = list(
        c(
            "WT NCD" = "#B8D6B6",
            "Sucnr1KO global NCD" = "#095256",
            "WT HFD" = "#D6A99A",
            "Sucnr1KO global HFD" = "#880D1E"
        ),
        "diet_genotype"
    ),
    "phase" = list(
        c(
            "G1" = "#b96e69",
            "S" = "#a4d251",
            "G2M" = "#8096d3"
        ),
        "Phase"
    )
)


for (target_mtd in c("aggregate_annotation", "stable_clusters_21")) {
    for (group_mtd in names(pal)) {
        current_df <- aggr_so@meta.data[, c(target_mtd, pal[[group_mtd]][[2]])]
        cont_table <- table(current_df[, 1], current_df[, 2])
        val_df <- melt(cont_table)
        colnames(val_df) <- c("annotation", "group_mtd", "ncells")
        val_df$annotation <- factor(val_df$annotation, levels = rownames(cont_table))
        val_df$group_mtd <- factor(val_df$group_mtd, levels = rev(names(pal[[group_mtd]][[1]])))

        cont_table <- cont_table / rowSums(cont_table) * 100
        val_df$percentage <- melt(cont_table)$value


        pdf(file.path(figure_path, paste0("barplot_aggregate_", target_mtd, "_by_", pal[[group_mtd]][[2]], ".pdf")), width = 5 + length(unique(current_df[,1])) * 0.5, height = 13)
        print(ggplot(val_df, aes(x = .data$annotation, y = .data$percentage, fill = .data$group_mtd)) +
            geom_bar(stat = "identity", position = "fill") +
            scale_fill_manual(values = pal[[group_mtd]][[1]]) +
            geom_fit_text(aes(label = .data$ncells), position = position_fill(vjust = 0.5), contrast = TRUE, angle = 90, size = 20, show.legend = FALSE) +
            theme_classic() +
            theme(
                axis.text = element_text(size = 20),
                legend.text = element_text(size = 15),
                title = element_text(size = 25),
                axis.text.x = element_text(angle = 45, hjust = 1),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.title = element_blank()
            ) +
            guides(
                fill = guide_legend(reverse = TRUE)
            )
        )

        dev.off()

        cont_table <- table(current_df[, 1], current_df[, 2])
        for (i in seq_len(nrow(cont_table))) {
            cont_table[i, ] <- cont_table[i, ] / sum(cont_table[i, ]) * 100
        }

        chi2_mat <- matrix(NA, nrow = nrow(cont_table), ncol = nrow(cont_table))
        rownames(chi2_mat) <- rownames(cont_table)
        colnames(chi2_mat) <- rownames(cont_table)
        for (i in seq(from = 1, to = nrow(chi2_mat) - 1)) {
            for (j in seq(from = i+1, to = nrow(chi2_mat))) {
                temp_mat <- t(cont_table[c(i, j), ])
                chi2_mat[i, j] <- chisq.test(temp_mat)$p.value
                chi2_mat[j, i] <- chi2_mat[i, j]
            }
        }

        col_fun <- circlize::colorRamp2(c(0, 0.05, 1), c("red", "white", "blue"))
        pdf(file.path(figure_path, paste0("heatmap_chi2_aggregate_", target_mtd, "_vs_", pal[[group_mtd]][[2]], ".pdf")), width = 5 + 0.5 * nrow(chi2_mat), height = 5 + 0.5 * nrow(chi2_mat))
        print(
            Heatmap(
                chi2_mat,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                col = col_fun,
                cell_fun = function(j, i, x, y, width, height, fill) {
                    if (is.na(chi2_mat[i,j]) || chi2_mat[i, j] > 0.05) {
                        return(grid.text("", x, y, gp = gpar(fontsize = 20)))
                    }

                    if (chi2_mat[i, j] > 0.01) {
                        return(grid.text("*", x, y, gp = gpar(fontsize = 20)))
                    }

                    if (chi2_mat[i, j] > 0.001) {
                        return(grid.text("**", x, y, gp = gpar(fontsize = 20)))
                    }
                    if (chi2_mat[i, j] > 0.0001) {
                        return(grid.text("***", x, y, gp = gpar(fontsize = 20)))
                    }

                    return(grid.text("****", x, y, gp = gpar(fontsize = 20)))
                },
            )
        )
        dev.off()
    }
}

# stacked barplot cluster by annotation
pal <- c(
    "OPSc" = "#BD7293",
    "Oligodendrocytes" = "#C9AC81",
    "Microglia" = "#EB573B",
    "Astrocytes" = "#88C879",
    "Multicellular niche" = "#706EBE",
    "Neurons" = "#7DB2C4"
)
cont_table <- table(aggr_so$stable_clusters_21, aggr_so$aggregate_annotation)
cont_table <- cont_table / rowSums(cont_table)
val_df <- melt(cont_table)
colnames(val_df) <- c("Clusters", "Cell types", "Fraction of cells")
val_df$Clusters <- factor(val_df$Clusters, levels = as.character(seq_len(nlevels(aggr_so$stable_clusters_21))))

pdf(file.path(figure_path, "barplot_group_cluster_by_annotation.pdf"), width = 10, height = 10)
print(ggplot(val_df, aes(x = .data$Clusters, y = .data$"Fraction of cells", fill = .data$"Cell types")) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = pal) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        title = element_text(size = 25),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # legend.position = "bottom",
        # legend.direction = "horizontal",
        legend.title = element_blank()
    )
)
dev.off()

cont_table <- table(aggr_so[, aggr_so$diet == "HFD"]$genotype, aggr_so[, aggr_so$diet == "HFD"]$aggregate_annotation)
cont_table <- cont_table / rowSums(cont_table)
val_df <- melt(cont_table)
colnames(val_df) <- c("Genotype", "Cell types", "Fraction of cells")

pdf(file.path(figure_path, "barplot_group_hfd_genotype_by_annotation.pdf"), width = 4, height = 10)
print(ggplot(val_df, aes(x = .data$Genotype, y = .data$"Fraction of cells", fill = .data$"Cell types")) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = pal) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        title = element_text(size = 25),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.title = element_blank()
    )
)
dev.off()

cont_table <- table(aggr_so[, aggr_so$diet == "NCD"]$genotype, aggr_so[, aggr_so$diet == "NCD"]$aggregate_annotation)
cont_table <- cont_table / rowSums(cont_table)
val_df <- melt(cont_table)
colnames(val_df) <- c("Genotype", "Cell types", "Fraction of cells")

pdf(file.path(figure_path, "barplot_group_ncd_genotype_by_annotation.pdf"), width = 4, height = 10)
print(ggplot(val_df, aes(x = .data$Genotype, y = .data$"Fraction of cells", fill = .data$"Cell types")) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = pal) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        title = element_text(size = 25),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.title = element_blank()
    )
)
dev.off()

###### Heatmap sample statistics ######
count_mat <- GetAssayData(aggr_so, layer = "counts")
all_genes <- rownames(count_mat)
gene_conv <- select(EnsDb.Mmusculus.v79, keys = all_genes, columns = c("GENENAME", "GENEBIOTYPE"), keytype = "GENENAME")
gene_conv$GENEBIOTYPE[gene_conv$GENEBIOTYPE == "protein_coding"] <- "PC"
gene_conv$GENEBIOTYPE[gene_conv$GENEBIOTYPE != "PC"] <- "NPC"
unique_groups <- unique(gene_conv$GENEBIOTYPE)
gene_groups <- lapply(unique_groups, function(x) {
    gene_conv[gene_conv$GENEBIOTYPE == x, "GENENAME"]
})
names(gene_groups) <- unique_groups

head(aggr_so@meta.data)
aggr_so$sample_name

categ_mtds <- c("sample_name", "diet_genotype", "aggregate_annotation", "stable_clusters_21")
for (mtd in categ_mtds) {
    if (is.factor(aggr_so[[mtd]])) {
        unique_values <- levels(aggr_so[[mtd]][, 1])
    } else {
        unique_values <- unique(aggr_so[[mtd]][, 1])
    }
    unique_values <- stringr::str_sort(unique_values, numeric = TRUE)

    psdb_count_mat <- matrix(0, nrow = length(unique_groups), ncol = length(unique_values))
    colnames(psdb_count_mat) <- unique_values
    rownames(psdb_count_mat) <- unique_groups
    for (i in seq_along(unique_groups)) {
        for (j in seq_along(unique_values)) {
            psdb_count_mat[i, j] <- sum(count_mat[gene_groups[[i]], aggr_so[[mtd]][, 1] == unique_values[j], drop = FALSE])
        }
    }

    chi2_mat <- matrix(NA, nrow = length(unique_values), ncol = length(unique_values))
    cramer_v_mat <- matrix(NA, nrow = length(unique_values), ncol = length(unique_values))
    colnames(chi2_mat) <- unique_values
    rownames(chi2_mat) <- unique_values
    colnames(cramer_v_mat) <- unique_values
    rownames(cramer_v_mat) <- unique_values

    for (i in seq(from = 1, to = length(unique_values) - 1)) {
        for (j in seq(from = i+1, to = length(unique_values))) {
            temp_mat <- psdb_count_mat[, c(i, j)]

            cramer_v_mat[i, j] <- cramersv(temp_mat)
            cramer_v_mat[j, i] <- cramer_v_mat[i, j]

            temp_mat[, 1] <- temp_mat[, 1] / sum(temp_mat[, 1]) * 100
            temp_mat[, 2] <- temp_mat[, 2] / sum(temp_mat[, 2]) * 100

            chi2_mat[i, j] <- chisq.test(temp_mat)$p.value
            chi2_mat[j, i] <- chi2_mat[i, j]
        }
    }

    col_fun <- circlize::colorRamp2(c(0, 0.05, 1), c("red", "white", "blue"))
    pdf(file.path(figure_path, paste0("heatmap_chi2_on_reads_aggregate_", mtd, ".pdf")), width = 5 + 0.5 * nrow(chi2_mat), height = 5 + 0.5 * nrow(chi2_mat))
    print(Heatmap(
        chi2_mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (is.na(chi2_mat[i,j]) || chi2_mat[i, j] > 0.05) {
                return(grid.text("", x, y, gp = gpar(fontsize = 20)))
            }

            if (chi2_mat[i, j] > 0.01) {
                return(grid.text("*", x, y, gp = gpar(fontsize = 20)))
            }

            if (chi2_mat[i, j] > 0.001) {
                return(grid.text("**", x, y, gp = gpar(fontsize = 20)))
            }
            if (chi2_mat[i, j] > 0.0001) {
                return(grid.text("***", x, y, gp = gpar(fontsize = 20)))
            }

            return(grid.text("****", x, y, gp = gpar(fontsize = 20)))
        },
    ))
    dev.off()
}

###### Bubbleplot ######
bubbleplot_fn <- function(so, genes, thresh = 1) {
    exp_mat <- GetAssayData(so, layer = "data")[genes, ]
    target_mtd <- droplevels(Idents(so))
    exp_mat_aggr <- matrix(0, nrow = length(genes), ncol = nlevels(target_mtd))
    perc_expr <- matrix(0, nrow = length(genes), ncol = nlevels(target_mtd))
    rownames(exp_mat_aggr) <- genes
    rownames(perc_expr) <- genes
    colnames(exp_mat_aggr) <- levels(target_mtd)
    colnames(perc_expr) <- colnames(exp_mat_aggr)
    for (i in seq_len(nlevels(target_mtd))) {
        for (j in seq_len(nrow(exp_mat))) {
            exp_mat_aggr[j, i] <- mean(exp_mat[j, target_mtd == colnames(exp_mat_aggr)[i]])
            perc_expr[j, i] <- sum(exp_mat[j, target_mtd == colnames(exp_mat_aggr)[i]] > 0) / sum(target_mtd == colnames(exp_mat_aggr)[i])
        }
    }

    exp_mat_aggr <- t(scale(t(exp_mat_aggr)))
    exp_mat_aggr[exp_mat_aggr > thresh] <- thresh
    exp_mat_aggr[exp_mat_aggr < -thresh] <- -thresh

    df1 <- melt(exp_mat_aggr)
    colnames(df1) <- c("gene", "diet_genotype", "scaled expression")
    df2 <- melt(perc_expr)
    df1$percentage_expressed <- df2$value * 100

    ggplot(df1, aes(x = .data$diet_genotype, y = .data$gene, colour = .data$`scaled expression`, size = .data$percentage_expressed)) +
        geom_point() +
        scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        labs(x = "", y = "") +
        theme_classic() +
        theme(
            axis.text = element_text(size = 20),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 18),
        ) +
        scale_size_continuous(range = c(1, 10))
}


genes <- c("Fkbp5", "Cst3", "Pid1", "Man1a", "Picalm", "B4galt1", "Srgap2", "Ivns1abp", "Runx1", "Tgfbr1", "Lrp1b", "Dlg2")
genes <- intersect(genes, rownames(aggr_so))
msk <- aggr_so$stable_clusters_21 == "19"

Idents(aggr_so) <- "diet_genotype"



pdf(file.path(figure_path, "bubbleplot_genes.pdf"), width = 10, height = 7)
bubbleplot_fn(aggr_so[, msk], genes) 
dev.off()

# set 2
genes <- c("Il1a", "Il1b", "Il4", "Il6", "Il12b", "Il27", "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5", "Cxcl9", "Cxcl10", "Cxcl12", "Cxcl16", "Tnfa", "Ccl4", "Ccl5", "Cd68", "Nfkb1")
genes <- intersect(genes, rownames(aggr_so))
msk <- aggr_so$stable_clusters_21 == "19"
pdf(file.path(figure_path, "bubbleplot_microglia_markers.pdf"), width = 10, height = 7)
bubbleplot_fn(aggr_so[, msk], genes) 

genes <- c(genes, c("Siglech", "P2ry12", "Fcrl1", "Cx3Cr1"))
genes <- intersect(genes, rownames(aggr_so))

bubbleplot_fn(aggr_so[, msk], genes) 
dev.off()

