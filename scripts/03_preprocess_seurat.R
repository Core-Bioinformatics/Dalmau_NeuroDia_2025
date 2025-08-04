setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(qs)
library(foreach)

objects_folder <- file.path(project_folder, "objects", "R", "seurat")
metadata_folder <- file.path(project_folder, "metadata")
config_csv <- read.csv(file.path(metadata_folder, "configurations.csv"), comment.char = "#", row.names = 1)

for (cname in rownames(config_csv)) {
    print(cname)
    input_path <- file.path(objects_folder, paste0(cname, ".qs"))
    output_path <- file.path(objects_folder, paste0(cname, "_sct.qs"))

    if (file.exists(output_path) && file.size(output_path) > 0) {
        next
    }

    so <- qread(input_path, nthreads = ncores)
    DefaultAssay(so) <- "RNA"

    n_median_genes <- median(so$nFeature_RNA)
    n_median_genes <- floor(n_median_genes / 500) * 500

    so <- SCTransform(so, return.only.var.genes = FALSE, variable.features.n = n_median_genes, verbose = TRUE)
    DefaultAssay(so) <- "SCT"
    so <- CellCycleScoring(so, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
    RhpcBLASctl::blas_set_num_threads(ncores)
    so <- RunPCA(so, npcs = 30, approx = FALSE, verbose = TRUE)
    so <- RunUMAP(so, reduction = "pca", dims = 1:30, verbose = TRUE, umap.method = "uwot")
    qsave(so, output_path, nthreads = ncores)
}
