setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(ClustAssess)
library(qs)

metadata_folder <- file.path(project_folder, "metadata")
ca_app_folder <- file.path(project_folder, "shiny_apps", "clustassess")

anns <- rjson::fromJSON(file = file.path(metadata_folder, "annotations.json"))

for (object_name in names(anns)) {
    print(object_name)
    seurat_path <- file.path(project_folder, "objects", "R", "seurat", paste0(object_name, "_sct.qs"))

    so <- qread(seurat_path, nthreads = ncores)
    current_mtd <- so@meta.data

    for (ann in anns[[object_name]]) {
        print(ann)
        name_existing <- ann$name_existing
        name_new <- ann$name_new
        voting_scheme_filtering <- ann$voting_scheme_filtering
        if (is.null(voting_scheme_filtering)) {
            voting_scheme_filtering <- 0
        }
        mapping <- ann$mapping

        new_mtd <- as.character(current_mtd[[name_existing]])
        names(new_mtd) <- colnames(so)

        if (voting_scheme_filtering == 0) {
            for (new_categ in names(mapping)) {
                new_mtd[current_mtd[[name_existing]] %in% as.character(mapping[[new_categ]])] <- new_categ
            }
        } else {
            expr_matrix <- GetAssayData(so, layer = "data")
            for (new_categ in names(mapping)) {
                current_mp <- mapping[[new_categ]]
                cells_from_mtd <- current_mtd[[name_existing]] %in% as.character(current_mp$mtd)
                gene_expr <- expr_matrix[current_mp$genes, cells_from_mtd, drop = FALSE]
                gene_expr <- colSums(gene_expr > current_mp$threshold)
                remaining_cells <- names(gene_expr)[gene_expr >= (length(current_mp$genes) - current_mp$relaxation)]

                if (length(remaining_cells) > 0) {
                    new_mtd[remaining_cells] <- new_categ
                }
            }
        }

        current_mtd[[name_new]] <- factor(new_mtd)
    }

    so@meta.data <- current_mtd
    qsave(so, seurat_path, nthreads = ncores)

    new_names <- sapply(anns[[object_name]], function(ann) ann$name_new)
    ClustAssess::add_metadata(file.path(ca_app_folder, object_name), current_mtd[ , new_names, drop = FALSE])
}

