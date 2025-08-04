setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(ClustAssess)
library(qs)
library(rhdf5)

objects_folder <- file.path(project_folder, "objects", "R")
so_folder <- file.path(objects_folder, "seurat")
ca_folder <- file.path(objects_folder, "clustassess")

mtd_folder <- file.path(project_folder, "metadata")
shiny_folder <- file.path(project_folder, "shiny_apps", "clustassess")
if (!dir.exists(shiny_folder)) {
    dir.create(shiny_folder, recursive = TRUE)
}

mtd_configs <- read.csv(file.path(mtd_folder, "configurations.csv"), header = TRUE, comment.char = "#")

get_abbv <- function(ftype, fsize) {
    words <- strsplit(ftype, "_")[[1]]
    words <- sapply(words, function(word) stringr::str_to_lower(substr(word, 1, 1)))

    return(paste0(paste(words, collapse = ""), "_", fsize))
}

for (i in seq_len(nrow(mtd_configs))) {
    id <- mtd_configs$id[i]
    ftype <- mtd_configs$ftype[i]

    if (is.null(ftype) || ftype == "") {
        next
    }

    fsize <- as.character(mtd_configs$fsize[i])
    clmethod <- mtd_configs$clmethod[i]
    stable_k <- strsplit(as.character(mtd_configs$stable_k[i]), ";")[[1]]
    chosen_k <- as.character(mtd_configs$chosen_k[i])
    print(paste(id, ftype, fsize, clmethod, paste(stable_k, collapse = ";"), chosen_k, sep = "----"))

    so <- qread(file.path(so_folder, paste0(id, "_sct.qs")), nthreads = ncores)
    DefaultAssay(so) <- "SCT"

    ca_obj_path <- file.path(ca_folder, paste0(id, ".qs"))
    ca_shiny_path <- file.path(shiny_folder, id, "stability.h5")

    ca_obj_exists <- file.exists(ca_obj_path)
    ca_shiny_exists <- file.exists(ca_shiny_path)
    if (!ca_obj_exists && !ca_shiny_exists) {
        print(paste0("ClustAssess object not found for ", id, ". Skipping."))
        next
    }

    mtd_names <- colnames(so@meta.data)
    for (target_mtd_name in mtd_names[grepl("^ecc_|^stable_clusters_", mtd_names)]) {
        so@meta.data[[target_mtd_name]] <- NULL
    }

    if (ca_obj_exists) {
        ca_obj <- qread(file.path(ca_folder, paste0(id, ".qs")), nthreads = ncores)
        ca_obj <- ca_obj[[ftype]][[fsize]]
        pca_emb <- ca_obj$pca
        umap_emb <- ca_obj$umap
        ca_obj <- ca_obj$clustering_stability$split_by_k[[clmethod]]

        for (k_val in stable_k) {
            so@meta.data[[paste0("stable_clusters_", k_val)]] <- factor(ca_obj[[k_val]]$partitions[[1]]$mb)
            so@meta.data[[paste0("ecc_", k_val)]] <- ca_obj[[k_val]]$ecc
        }

        if (!is.null(chosen_k) && chosen_k %in% stable_k) {
            so@meta.data[[paste0(id, "_stable_clusters")]] <- factor(ca_obj[[chosen_k]]$partitions[[1]]$mb)
        }
    } else {
        pca_emb <- h5read(ca_shiny_path, paste0("/", ftype, "/", fsize, "/pca"))
        umap_emb <- h5read(ca_shiny_path, paste0("/", ftype, "/", fsize, "/umap"))

        for (k_val in stable_k) {
            so@meta.data[[paste0("stable_clusters_", k_val)]] <- factor(h5read(ca_shiny_path, paste0("/", ftype, "/", fsize, "/clustering_stability/split_by_k/mbs/", clmethod, "/", k_val)))
            ecc_val <- h5read(ca_shiny_path, paste0("/", ftype, "/", fsize, "/clustering_stability/split_by_k/ecc/", sprintf("%06d", as.integer(k_val)), ";", clmethod))
            ecc_order <- h5read(ca_shiny_path, paste0("/", ftype, "/", fsize, "/clustering_stability/split_by_k/ecc_order/", sprintf("%06d", as.integer(k_val)), ";", clmethod))
            so@meta.data[[paste0("ecc_", k_val)]] <- ecc_val[ecc_order]
        }

        if (!is.null(chosen_k) && chosen_k %in% stable_k) {
            so@meta.data[[paste0(id, "_stable_clusters")]] <- factor(h5read(ca_shiny_path, paste0("/", ftype, "/", fsize, "/clustering_stability/split_by_k/mbs/", clmethod, "/", chosen_k)))
        }
    }

    if (length(so@reductions) == 0) {
        so <- RunPCA(so, npcs = ncol(pca_emb))
        so <- RunUMAP(so, reduction = "pca", dims = seq_len(ncol(pca_emb)))
    }

    old_cln <- colnames(so@reductions$pca@cell.embeddings)[seq_len(ncol(pca_emb))]
    so@reductions$pca@cell.embeddings <- pca_emb
    colnames(so@reductions$pca@cell.embeddings) <- old_cln
    rownames(so@reductions$pca@cell.embeddings) <- colnames(so)

    old_cln <- colnames(so@reductions$umap@cell.embeddings)
    so@reductions$umap@cell.embeddings <- umap_emb
    colnames(so@reductions$umap@cell.embeddings) <- old_cln
    rownames(so@reductions$umap@cell.embeddings) <- colnames(so)

    # qsave(so, file.path(so_folder, paste0(id, "_", get_abbv(ftype, fsize), "_sct.qs")), nthreads = ncores)
    qsave(so, file.path(so_folder, paste0(id, "_sct.qs")), nthreads = ncores)
}

FeaturePlot(so, features = "ecc_26") + scale_colour_viridis_c()
DimPlot(so, group.by = "stable_clusters_26", label = TRUE)
