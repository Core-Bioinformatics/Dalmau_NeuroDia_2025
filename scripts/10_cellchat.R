setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(qs)
library(CellChat)
library(patchwork)

# we follow this
# https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
options(stringsAsFactors = FALSE)

seurat_path <- file.path(project_folder, "objects", "R", "seurat")
cc_path <- file.path(project_folder, "objects", "R", "cellchat")
if (!dir.exists(cc_path)) {
    dir.create(cc_path, recursive = TRUE)
}

cc_output_path <- file.path(project_folder, "output", "cellchat")
if (!dir.exists(cc_output_path)) {
    dir.create(cc_output_path, recursive = TRUE)
}

cc_shiny_path <- file.path(project_folder, "shiny_apps", "cellchat")
if (!dir.exists(cc_shiny_path)) {
    dir.create(cc_shiny_path, recursive = TRUE)
}

so <- qread(file.path(seurat_path, "neurons_microglia_sct.qs"), nthreads = ncores)
good_groups <- setdiff(unique(so$neurons_microglia_annotation), 1:26)
broad_groups <- c("GABAergic (Inhibitory Neurons)", "Glutaminergic (Exhibitory Neurons)", "Microglia")
focused_groups <- c(setdiff(good_groups, broad_groups), "Microglia")
so <- subset(so, subset = neurons_microglia_annotation %in% good_groups)
so_types <- list(
  "broad" = subset(so, subset = neurons_microglia_annotation %in% broad_groups),
  "focused" = subset(so, subset = neurons_microglia_annotation %in% focused_groups)
)

app_text <- paste0("
library(qs)
library(CellChat)

cellchat <- qread('cellchat.qs')
runCellChatApp(cellchat)
")

for (stype in names(so_types)) {

  so_list <- list(
    "neurons_microglia" = so_types[[stype]],
    "neurons_microglia_wt_ncd" = subset(so_types[[stype]], genotype == "wild-type" & diet == "NCD"),
    "neurons_microglia_wt_hfd" = subset(so_types[[stype]], genotype == "wild-type" & diet == "HFD"),
    "neurons_microglia_ko_ncd" = subset(so_types[[stype]], genotype == "Sucnr1KO" & diet == "NCD"),
    "neurons_microglia_ko_hfd" = subset(so_types[[stype]], genotype == "Sucnr1KO" & diet == "HFD")
  )

  if (!dir.exists(file.path(cc_path, stype))) {
    dir.create(file.path(cc_path, stype), recursive = TRUE)
  }

  for (grp_name in names(so_list)) {
    so_list[[grp_name]]$neurons_microglia_annotation <- droplevels(so_list[[grp_name]]$neurons_microglia_annotation)

    cc_obj <- file.path(cc_path, stype, paste0(grp_name, ".qs"))
    if (!file.exists(cc_obj)) {
      cellChat <- createCellChat(object = so_list[[grp_name]], group.by = "neurons_microglia_annotation", assay = "SCT")


      Idents(so_list[[grp_name]]) <- "neurons_microglia_annotation"
      labels <- Idents(so_list[[grp_name]])
      meta <- data.frame(labels = labels, row.names = names(labels))

      cellChat <- addMeta(cellChat, meta = meta)
      cellChat <- setIdent(cellChat, ident.use = "labels")
      levels(cellChat@idents)
      groupSize <- as.numeric(table(cellChat@idents))


      CellChatDB <- CellChatDB.mouse
      showDatabaseCategory(CellChatDB)
      cellChat@DB <- subsetDB(CellChatDB)
      showDatabaseCategory(cellChat@DB)

      # Preprocessing the expression data for cell-cell communication analysis
      cellChat <- subsetData(cellChat)
      cellChat <- identifyOverExpressedGenes(cellChat)
      cellChat <- identifyOverExpressedInteractions(cellChat)

      # Compute the communication probability and infer cellular communication network
      cellChat <- computeCommunProb(cellChat, type = "triMean")
      cellChat <- filterCommunication(cellChat, min.cells = 10)

      # Infer the cell-cell communication at a signaling pathway level
      cellChat <- computeCommunProbPathway(cellChat)

      # Calculate the aggregated cell-cell communication network
      cellChat <- aggregateNet(cellChat)
      cellChat <- addReduction(cellChat, seu.obj = so_list[[grp_name]], dr.use = "umap")
      qs::qsave(cellChat, file = cc_obj, nthreads = ncores)
    } else {
      cellChat <- qread(cc_obj, nthreads = ncores)
    }

    app_path <- file.path(cc_shiny_path, stype, grp_name)
    if (!dir.exists(app_path)) {
      dir.create(app_path, recursive = TRUE)
    }
    qsave(cellChat, file.path(app_path, "cellchat.qs"), nthreads = ncores)
    write(app_text, file.path(app_path, "app.R"))
  }

  # groupSize <- as.numeric(table(cellChat@idents))
  # par(mfrow = c(1,2), xpd=TRUE)
  # netVisual_circle(cellChat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  # netVisual_circle(cellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

  # mat <- cellChat@net$weight
  # par(mfrow = c(ceiling(nlevels(Idents(so)) / 3),3), xpd=TRUE)
  # for (i in 1:nrow(mat)) {
  #   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  #   mat2[i, ] <- mat[i, ]
  #   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  # }

  # # Visualization of cell-cell communication network
  # print(sort(cellChat@netP$pathways))
  # pathways.show <- c("PTPRM")
  # vertex.receiver = seq(1,4) # a numeric vector. 
  # netVisual_aggregate(cellChat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  # # Circle plot
  # par(mfrow=c(1,1))
  # netVisual_aggregate(cellChat, signaling = pathways.show, layout = "circle")

  # # Chord diagram
  # par(mfrow=c(1,1))
  # netVisual_aggregate(cellChat, signaling = pathways.show, layout = "chord")

  # # Heatmap
  # par(mfrow=c(1,1))
  # netVisual_heatmap(cellChat, signaling = pathways.show, color.heatmap = "Reds")

  # runCellChatApp(cellChat)
}
