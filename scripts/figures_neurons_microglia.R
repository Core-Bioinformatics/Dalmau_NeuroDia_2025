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
figure_path <- file.path(project_folder, "output", "figures", "neurons_microglia")
shiny_path <- file.path(project_folder, "shiny_apps", "clustassess")

if (!dir.exists(figure_path)) {
    dir.create(figure_path, recursive = TRUE)
}

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


aggr_so <- qread(file.path(object_path, "seurat", "neurons_microglia_sct.qs"), nthreads = ncores) 

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

# clusters
for (k in c(26)) {
    pal <- h5read(file.path(shiny_path, "neurons_microglia", "stability.h5"), paste0("/colors/", k))
    names(pal) <- as.character(seq_len(k))
    cl <- aggr_so[[paste0("stable_clusters_", k)]][, 1]
    cl <- factor(cl, levels = seq_len(k))


    pdf(file.path(figure_path, paste0("umap_clusters_", k, ".pdf")), width = 10, height = 10)
    # for this plot, I want to add dotted ellipses around groups defined by neuron aggregation
    target_groups <- levels(aggr_so$neurons_microglia_annotation)
    ellipse_coords <- list()
    for (tgroup in target_groups) {
        if (tgroup %in% as.character(1:26)) {
            next
        }
        mask <- aggr_so$neurons_microglia_annotation == tgroup
        subumap <- aggr_so[, mask]@reductions$umap@cell.embeddings

        # define the ellipse that would fit the points from subumap
        ellipse_coords[[tgroup]] <- car::dataEllipse(subumap, levels = 0.99)
    }

    # draw the UMAP with ellipses
    ggobj <- default_umap_ggplot(
        aggr_so@reductions$umap@cell.embeddings,
        cl,
        mtd_name = paste0("Stable ", k, " Clusters"),
        colour_palette = pal,
        label_size = 3
    )
    print(ggobj)

    for (tgroup in target_groups) {
        if (tgroup %in% as.character(1:26)) {
            next
        }
        if (tgroup %in% c(
            "GABAergic (Inhibitory Neurons)",
            "Glutaminergic (Exhibitory Neurons)"
        )) {
            print("skip")
            next
        }
        coords <- ellipse_coords[[tgroup]]
        ggobj <- ggobj + geom_polygon(
            data = coords,
            aes(x = .data$x, y = .data$y),
            fill = NA,
            color = "black",
            linetype = "dotted",
            linewidth = 1
        )
    }
    print(ggobj)

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


for (target_mtd in c("stable_clusters_26")) {
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

##### Cell type violin plot #####
subso <- aggr_so[, aggr_so$stable_clusters_26 == "18"]
expr_mat <- GetAssayData(subso, layer = "data")
gene_cell_types <- rjson::fromJSON(file = file.path(project_folder, "metadata", "annotations.json"))

final_df <- NULL
for (ctype_name in names(gene_cell_types$neurons_microglia[[1]]$mapping)) {
    genes <- gene_cell_types$neurons_microglia[[1]]$mapping[[ctype_name]]$genes
    # thresh <- gene_cell_types$neurons_microglia[[1]]$mapping[[ctype_name]]$threshold
    thresh <- 0
    # relax <- gene_cell_types$neurons_microglia[[1]]$mapping[[ctype_name]]$relaxation
    relax <- 0

    submat <- expr_mat[genes, , drop = FALSE]
    mask <- colSums(submat > thresh) >= (length(genes) - relax)

    if (sum(mask) == 0) {
        next
    }

    temp_df <- data.frame(
        average_gene_expr = colMeans(submat[, mask, drop = FALSE])
    )
    temp_df$cell_type <- ctype_name

    if (is.null(final_df)) {
        final_df <- temp_df
    } else {
        final_df <- rbind(final_df, temp_df)
    }
}

ncells <- final_df %>%
    group_by(cell_type) %>%
    summarise(ncells = n())

pdf(file.path(figure_path, "violin_cell_type_on_stable_18.pdf"), width = 15, height = 10)
ggplot(final_df, aes(x = cell_type, y = average_gene_expr)) +
    geom_violin(aes(fill = cell_type), alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75)) +
    # scale_fill_manual(values = pal) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    # put the number of cells above the violin
    geom_text(
        data = ncells,
        aes(x = cell_type, y = 0.5, label = paste0("n=", ncells)),
        size = 5,
        vjust = -0.5
    ) +
    ggtitle("Cell type expression")
dev.off()

##### Gene expression violin plot #####
otp_path <- file.path(figure_path, "violin_gene_expression")
if (!dir.exists(otp_path)) {
    dir.create(otp_path, recursive = TRUE)
}

genes <- c("Lepr", "Ghsr", "Thra", "Thrb", "Glp1r", "Nr3c1", "Insr", "Gcg", "Acaca", "Fasn", "Pparg", "Srefb1", "Cd36", "Fabp4", "Scarb1", "Ldlr", "Pnpla2", "Lipe", "Cpt1a", "Gck", "Pck1", "Gys2", "Pygl", "Slc2a2", "Slc2a1", "Drd1", "Drd2", "Drd3", "Drd4", "Drd5")
genes <- intersect(genes, rownames(aggr_so))
pdf(file.path(otp_path, "violin_gene_expression_horomone_lipid_glucose_dopamine_markers.pdf"), width = 30, height = 10)
for (g in genes) {
    print(
        Seurat::VlnPlot(
            aggr_so,
            features = g,
            group.by = "stable_clusters_26",
            split.by = "diet_genotype",
            log = FALSE,
            pt.size = 0
        )
    )
}
dev.off()

genes <- c("Ntrk2", "Snx32", "Fgf1", "Cpt1c", "Pde10a", "Snca")
genes <- intersect(genes, rownames(aggr_so))
pdf(file.path(otp_path, "violin_gene_expression_cluster_18_markers.pdf"), width = 30, height = 10)
for (g in genes) {
    print(
        Seurat::VlnPlot(
            aggr_so,
            features = g,
            group.by = "stable_clusters_26",
            split.by = "diet_genotype",
            log = FALSE,
            pt.size = 0
        )
    )
}
dev.off()

genes <- c("Gphn", "Dlg4", "Slc32a1", "Slc17a7", "Psd95", "Vgat", "Vglut")
genes <- intersect(genes, rownames(aggr_so))
pdf(file.path(otp_path, "violin_gene_expression_custom_from_luca.pdf"), width = 30, height = 10)
for (g in genes) {
    print(
        Seurat::VlnPlot(
            aggr_so,
            features = g,
            group.by = "stable_clusters_26",
            split.by = "diet_genotype",
            log = FALSE,
            pt.size = 0
        )
    )
}
dev.off()

##### Bubbleplot #####
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

# first set
genes <- c("Lepr", "Ghsr", "Thra", "Thrb", "Glp1r", "Nr3c1", "Insr", "Gcg", "Acaca", "Fasn", "Pparg", "Srefb1", "Cd36", "Fabp4", "Scarb1", "Ldlr", "Pnpla2", "Lipe", "Cpt1a", "Gck", "Pck1", "Gys2", "Pygl", "Slc2a2", "Slc2a1", "Drd1", "Drd2", "Drd3", "Drd4", "Drd5")
genes <- intersect(genes, rownames(aggr_so))
Idents(aggr_so) <- aggr_so$diet_genotype

pdf(file.path(figure_path, "bubbleplot_horomone_lipid_glucose_dopamine_markers.pdf"), width = 10, height = 12)
for (i in sort(as.numeric(levels(aggr_so$stable_clusters_26)))) {
    print(
        bubbleplot_fn(
            aggr_so[, aggr_so$stable_clusters_26 == i],
            genes
        ) + ggtitle(paste0("Stable cluster ", i))
    )
}
print(
    bubbleplot_fn(
        aggr_so,
        genes
    ) + ggtitle("All stable clusters")
)
dev.off()

# second set
genes <- c("Ntrk2", "Snx32", "Fgf1", "Cpt1c", "Pde10a", "Snca")
genes <- intersect(genes, rownames(aggr_so))
pdf(file.path(figure_path, "bubbleplot_cluster_18_markers.pdf"), width = 8, height = 6)
print(
    bubbleplot_fn(
        aggr_so[, aggr_so$stable_clusters_26 == "18"],
        genes
    )
)
genes <- c(genes, c("Lepr", "Drd1", "Drd2", "Drd3", "Drd4", "Drd5"))
genes <- intersect(genes, rownames(aggr_so))

print(
    bubbleplot_fn(
        aggr_so[, aggr_so$stable_clusters_26 == "18"],
        genes
    )
)
dev.off()

# third set
genes <- c("Snca", "Pcbd2", "Prdm15", "Leng9", "Dalrd3", "Snx32", "Nit1", "Pabpc4", "Lrrc3b", "Cnot3", "Ddx39b", "Ring1", "Nisch", "Snrpf", "Ascc2", "Slc25a3", "Ints3", "Snhg9", "Ilk", "Socs2")
genes <- intersect(genes, rownames(aggr_so))
pdf(file.path(figure_path, "bubbleplot_cluster_12_markers.pdf"), width = 8, height = 9)
print(
    bubbleplot_fn(
        aggr_so[, aggr_so$stable_clusters_26 == "12"],
        genes
    )
)
dev.off()


#### AGGR GENE UMAP #####

gene_list <- list(
    "Glutaminergic (Exhibitory Neurons)" = list(
        genes = c("Slc17a6"),
        thresh = 0.5,
        not_expr = -1
    ),
    "GABAergic (Inhibitory Neurons)" = list(
        genes = c("Slc32a1"),
        thresh = 0.5,
        not_expr = -1
    ),
    "Pomc Neurons" = list(
        genes = c("Pomc"),
        thresh = 0.7,
        not_expr = -1
    ),
    "Dopaminergic Neurons" = list(
        genes = "Th",
        thresh = 0.7,
        not_expr = -1
    ),
    "Suprachiasmatic nucleous (SCN) Neurons" = list(
        genes = c("Rgs16"),
        thresh = 1,
        not_expr = -1
    ),
    "Histaminergic Neurons" = list(
        genes = c("Hdc"),
        thresh = 1.5,
        not_expr = -1
    ),
    "MCH Neurons" = list(
        genes = c("Pmch"),
        thresh = 2,
        not_expr = -1
    ),
    "Cholinergic Neurons" = list(
        genes = c("Chat"),
        thresh = 0.6,
        not_expr = -1
    ),
    "MAO-B+ Neurons" = list(
        genes = c("Maob"),
        thresh = 0.7,
        not_expr = -1
    ),
    "Lepr Neurons" = list(
        genes = c("Lepr"),
        thresh = 1.5,
        not_expr = -1
    ),
    "Subthalamic Neurons" = list(
        genes = c("Pitx2"),
        thresh = 0.6,
        not_expr = -1
    ),
    "SF1 Neurons" = list(
        genes = c("Nr5a1"),
        thresh = 0.5,
        not_expr = -1
    ),
    "AgRP/Npy/Sst Neurons" = list(
        genes = c("Npy", "Sst", "Agrp"),
        thresh = 1,
        not_expr = 1
    ),
    "Glp1r Neurons" = list(
        genes = c("Glp1r"),
        thresh = 0.6,
        not_expr = -1
    ),
    "Parvalbumim Neurons" = list(
        genes = c("Pvalb"),
        thresh = 0.6,
        not_expr = -1
    ),
    "Hipocretin Neurons" = list(
        genes = c("Hcrt"),
        thresh = 1.75,
        not_expr = -1
    ),
    "Microglia" = list(
        genes = c("Cx3cr1", "P2ry12", "Tmem119", "Siglech", "Hexb"),
        thresh = 0,
        not_expr = 2
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


## special case : microglia + gabaergic + glutaminergic
avg_exprs <- lapply(c("GABAergic (Inhibitory Neurons)", "Glutaminergic (Exhibitory Neurons)", "Microglia"), function(ctype) {
    g <- gene_list[[ctype]]$genes
    thresh <- gene_list[[ctype]]$thresh
    not_expr <- gene_list[[ctype]]$not_expr
    if (not_expr == -1) {
        not_expr <- length(g) - 1
    }

    avg_expr <- colMeans(expr_mat[g, , drop = FALSE])
    mask <- colSums(expr_mat[g, , drop = FALSE] > thresh) >= (length(g) - not_expr)
    avg_expr[!mask] <- 0 
    return(avg_expr)
})
names(avg_exprs) <- c("GABAergic (Inhibitory Neurons)", "Glutaminergic (Exhibitory Neurons)", "Microglia")

colour_scales <- list(
    "GABAergic (Inhibitory Neurons)" = paletteer::paletteer_d("RColorBrewer::Blues", n = 9),
    "Glutaminergic (Exhibitory Neurons)" = paletteer::paletteer_d("RColorBrewer::Reds", n = 9),
    "Microglia" = paletteer::paletteer_d("RColorBrewer::Greens", n = 9)
)
umap_df <- data.frame(aggr_so@reductions$umap@cell.embeddings)
min_x <- min(umap_df$umap_1)
max_x <- max(umap_df$umap_1)
min_y <- min(umap_df$umap_2)
max_y <- max(umap_df$umap_2)

gpl <- ggplot() + geom_point(data = umap_df, aes(x = umap_1, y = umap_2), colour = "grey90", size = 0.5) + ggnewscale::new_scale_colour() 
for (ctype in names(avg_exprs)) {
    temp_df <- umap_df
    temp_df$expr <- avg_exprs[[ctype]]
    temp_df <- temp_df %>% dplyr::filter(!is.na(expr)) %>% dplyr::filter(expr > 0) %>% arrange(expr)
    gpl <- gpl + geom_point(
        data = temp_df,
        aes(x = umap_1, y = umap_2, colour = expr),
        size = 1
    ) + scale_colour_gradientn(
        colours = colour_scales[[ctype]],
        limits = c(0, 2),
        na.value = "grey85",
        name = ctype,
        guide = guide_colourbar(title.position = "top")
    ) + labs(title = ctype) 
    if (ctype != names(avg_exprs)[3]) {
        gpl <- gpl + ggnewscale::new_scale_colour()
    }
}

gpl <- gpl +
    theme_classic() +
    scale_x_continuous(breaks = c(min_x + (max_x - min_x) / 7)) +
    scale_y_continuous(breaks = c(min_y + (max_y - min_y) / 7)) +
    guides(
        x = guide_axis(cap = "upper"),
        y = guide_axis(cap = "upper")
    ) +
    labs(
        x = "UMAP1",
        y = "UMAP2"
    ) +
    ggtitle("") +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(arrow = arrow(type = "closed", length = unit(0.3, "cm"))),
        axis.title.x = element_text(hjust = 0.03),
        axis.title.y = element_text(hjust = 0.03),
        plot.title = element_text(size = 25),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 10),
        aspect.ratio = 1,
        legend.position = "bottom",
        legend.direction = "horizontal"
    ) 

pdf(file.path(gene_dir, "microglia_gaba_glutaminergic.pdf"), width = 10, height = 10)
print(gpl)
dev.off()
