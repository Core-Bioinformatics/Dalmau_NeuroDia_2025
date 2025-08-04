library(Seurat)
library(dplyr)
library(ggplot2)
library(ClustAssess)
library(qs)

if (interactive()) {
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
    setwd("scripts")
}

list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)

objects_folder <- file.path(project_folder, "objects", "R")
ca_folder <- file.path(objects_folder, "clustassess")
mtd_folder <- file.path(project_folder, "metadata")
if (!dir.exists(ca_folder)) {
    dir.create(ca_folder, recursive = TRUE)
}

mtd_configs <- read.csv(file.path(mtd_folder, "configurations.csv"), header = TRUE, comment.char = "#")
nreps <- 100
neigh_seq <- seq(from = 5, to = 50, by = 5)
res_seq <- seq(from = 0.1, to = 2, by = 0.1)
assay_name <- "SCT"

for (i in seq_len(nrow(mtd_configs))) {
    id <- mtd_configs$id[i]
    so_path <- file.path(objects_folder, "seurat", paste0(id, "_sct.qs"))
    ca_path <- file.path(ca_folder, paste0(id, ".qs"))

    if (file.exists(ca_path) && file.size(ca_path) > 0) {
        next
    }

    print(id)

    so <- qread(file.path(objects_folder, "seurat", paste0(id, "_sct.qs")), nthreads = ncores)
    DefaultAssay(so) <- assay_name

    # initialise the parallel context
    RhpcBLASctl::blas_set_num_threads(1)
    my_cluster <- parallel::makeCluster(
        ncores,
        type = "PSOCK"
    )
    doParallel::registerDoParallel(cl = my_cluster)

    # extract the matrix and the features
    expr_matrix <- so@assays[[assay_name]]@scale.data
    features <- dimnames(so@assays[[assay_name]])[[1]]
    var_features <- so@assays[[assay_name]]@var.features
    max_ngenes <- length(var_features)

    most_abundant_genes <- rownames(expr_matrix)[order(Matrix::rowSums(expr_matrix), decreasing = TRUE)]

    gene_list <- list(
        "Most_Abundant" = most_abundant_genes[seq_len(max_ngenes)],
        "Highly_Variable" = so@assays[[assay_name]]@var.features[seq_len(max_ngenes)]
    )

    fstep <- mtd_configs$fstep[i]
    steps_list <- list(
        "Most_Abundant" = seq(from = 500, by = fstep, to = max_ngenes),
        "Highly_Variable" = seq(from = 500, by = fstep, to = max_ngenes)
    )

    rm(so)
    gc()

    test_automm <- automatic_stability_assessment(
        expression_matrix = expr_matrix,
        n_repetitions = nreps,
        temp_file = "clustassess_temp.rds",
        save_temp = FALSE,
        n_neigh_sequence = neigh_seq,
        resolution_sequence = res_seq,
        features_sets = gene_list,
        steps = steps_list,
        n_top_configs = 2,
        umap_arguments = list(
            min_dist = 0.3,
            n_neighbors = 30,
            metric = "cosine"
        ),
        verbose = TRUE
    )

    parallel::stopCluster(cl = my_cluster)
    qsave(test_automm, ca_path, nthreads = ncores)

    so <- qread(so_path, nthreads = ncores)
    DefaultAssay(so) <- assay_name

    write_shiny_app(
        object = so,
        assay_name = assay_name,
        clustassess_object = test_automm,
        shiny_app_title = id,
        project_folder = file.path(project_folder, "shiny_apps", "clustassess", id),
        prompt_feature_choice = FALSE
    )
}
