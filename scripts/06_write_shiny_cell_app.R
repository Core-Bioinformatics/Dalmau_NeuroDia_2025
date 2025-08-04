setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(qs)
library(ShinyCell)

objects_folder <- file.path(project_folder, "objects", "R", "seurat")
mtd_folder <- file.path(project_folder, "metadata")
shiny_folder <- file.path(project_folder, "shiny_apps", "shiny_cell")
if (!dir.exists(shiny_folder)) {
    dir.create(shiny_folder, recursive = TRUE)
}

get_abbv <- function(ftype, fsize) {
    words <- strsplit(ftype, "_")[[1]]
    words <- sapply(words, function(word) stringr::str_to_lower(substr(word, 1, 1)))

    return(paste0(paste(words, collapse = ""), "_", fsize))
}

mtd_configs <- read.csv(file.path(mtd_folder, "configurations.csv"), header = TRUE, comment.char = "#")

for (i in seq_len(nrow(mtd_configs))) {
    # id <- paste0(mtd_configs$id[i], "_", get_abbv(mtd_configs$ftype[i], mtd_configs$fsize[i]))
    id <- mtd_configs$id[i]
    print(id)
    so <- qread(file.path(objects_folder, paste0(id, "_sct.qs")), nthreads = ncores)
    DefaultAssay(so) <- "SCT"

    # make sure stable clusters are factors, if present
    for (mtd in colnames(so@meta.data)) {
        if (startsWith(mtd, "stable_")) {
            so@meta.data[[mtd]] <- factor(so@meta.data[[mtd]])
        }
    }

    sc_conf <- createConfig(so)
    makeShinyApp(
        so,
        sc_conf,
        gene.mapping = TRUE,
        gex.assay = "SCT",
        shiny.dir = file.path(shiny_folder, id),
        shiny.title = paste0("RNA: ", id)
    )
}
