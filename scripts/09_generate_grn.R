setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
#library(magrittr)
library(tidyr)
library(data.table)
library(biomaRt)
library(org.Mm.eg.db)
library(Seurat)

object_path <- file.path(project_folder, "objects", "R")
figure_path <- file.path(project_folder, "output", "figures", "neurons_microglia", "grn", "pdf")
if (!dir.exists(figure_path)) {
    dir.create(figure_path, recursive = TRUE)
}

gene_list_temp <- list(
    "neurons" = c("Slc17a6", "Slc32a1", "Nr5a1", "Agrp", "Npy", "Sst", "Pomc", "Glp1r", "Pvalb", "Th", "Rgs16", "Hdc", "Hcrt", "Pmch", "Chat", "Maob", "Lepr", "Pitx2", "Pde10a", "Ntrk2"),
    "microglia" = c("Cx3cr1", "P2ry12", "Tmem119", "Siglech", "Hexb"),
    "hypocretin neurons" = c("Hcrt"),
    "hypocretin degs" = c("Ntrk2", "Pde10a", "Fgf1", "Cpt1c", "Snca", "Snx32"),
    "agrp npy sst neurons" = c("Agrp", "Npy", "Sst"),
    "pomc neurons" = c("Pomc"),
    "lepr neurons" = c("Lepr"),
    "sf1 neurons" = c("Nr5a1"),
    "suprachiasmatic neurons" = c("Rgs16"),
    "parvalbumin neurons" = c("Pvalb")
)

gene_list <- list()
gene_list[["microglia"]] <- gene_list_temp[["microglia"]]
gene_list[["neurons"]] <- unique(unlist(gene_list_temp[setdiff(names(gene_list_temp), "microglia")]))


gene_mapping <- mapIds(org.Mm.eg.db, keys = unique(unlist(gene_list)), column = "ENSEMBL", keytype = "SYMBOL")
gene_mapping_list <- lapply(gene_list, function(x) { gene_mapping[x]} )
names(gene_mapping_list) <- names(gene_list)


mapping <- mapIds(org.Mm.eg.db, keys = unlist(gene_mapping_list), column = "SYMBOL", keytype = "ENSEMBL")

# this I just took from bulkAnalyseR
get_link_list_rename <- function(weightMat, plotConnections){
  GENIE3::getLinkList(weightMat, plotConnections) %>%
    dplyr::mutate(from = as.character(.data$regulatoryGene), 
                  to = as.character(.data$targetGene), 
                  value = .data$weight, 
                  regulatoryGene = NULL, 
                  targetGene = NULL,
                  weight = NULL)
}

plot_GRN_multiomics <- function(weightMat, plotConnections){
  
  edges <- get_link_list_rename(weightMat, plotConnections)
  edges$id = paste0(pmin(edges$from,edges$to),pmax(edges$from,edges$to))
  edges = edges[order(-edges$value),]
  edges = edges[!duplicated(edges$id),]
  edges = edges[,1:3]
#   edges$color = '#B32F61'
  edges$color = 'black'

  # add nodes for each end of the edges and colour them by modality
  nodes <- tibble::tibble(
    id = c(edges$to, edges$from),
    label = c(edges$to, edges$from),
    value=1,
    font.color='black',
    font.size=50,
    font.bold=T
  ) %>%
    dplyr::distinct(.data$id, .keep_all = TRUE)

  nodes$value = 100
  
  ## combine the original edges and nodes with the pathway edges and nodes
  edges = edges[edges$from != edges$to,]
  nodes = nodes[nodes$id %in% c(edges$from,edges$to),]
  edges = edges[edges$to %in% nodes$id,]
  edges = edges[edges$from %in% nodes$id,]
  nodes = nodes[nodes$id %in% c(edges$from,edges$to),]
  # add line breaks for long strings
  nodes$label = gsub('_','\n',nodes$label)
  nodes$label = gsub('/','/\n',nodes$label)
  nodes = nodes[nodes$id!='',]

  
  # built the network using the edges and nodes you've defined
  # the physics bit depends how busy your network is, there's good documentation for it
  network <- visNetwork::visNetwork(nodes, edges,width=1080,height=1920) %>%

    visNetwork::visPhysics(solver = "forceAtlas2Based",
                           forceAtlas2Based = list(gravitationalConstant = -150,springConstant=0.1,centralGravity=0.005)) %>%
    visNetwork::visLegend()%>%
    visNetwork::visLayout(randomSeed = 2023)
  nodes$color = '#8bc0e0'
  return(list('nodes'=nodes,'edges'=edges,'network'=network))
}

plot_GRN_igraph <- function(weightMat, plotConnections){
  # get the normal edges and nodes from GENIE3
  edges <- get_link_list_rename(weightMat, plotConnections)
  # this bit removes duplicate edges where just the start and end are reversed
  edges$id = paste0(pmin(edges$from,edges$to),pmax(edges$from,edges$to))
  edges = edges[order(-edges$value),]
  edges = edges[!duplicated(edges$id),]
  edges = edges[,1:3]
  
  # add nodes for each end of the edges and colour them by modality
  nodes <- tibble::tibble(
    id = c(edges$to, edges$from),
    label = c(edges$to, edges$from),
    value=1,
    font.color='black',
    font.size=50,
    font.bold=T
  ) %>%
    dplyr::distinct(.data$id, .keep_all = TRUE)
  
  nodes$value = 100

  edges = edges[edges$from != edges$to,]
  nodes = nodes[nodes$id %in% c(edges$from,edges$to),]
  edges = edges[edges$to %in% nodes$id,]
  edges = edges[edges$from %in% nodes$id,]
  nodes = nodes[nodes$id %in% c(edges$from,edges$to),]
  # add line breaks for long strings
  nodes$label = gsub('_','\n',nodes$label)
  nodes$label = gsub('/','/\n',nodes$label)
  nodes = nodes[nodes$id!='',]
  edges$weight <- edges$value

  
  g <- igraph::graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)
  
  g <- igraph::set_vertex_attr(g, "color", value="#8bc0e0")
  
  return(g)
}


so <- qs::qread(file.path(object_path, "seurat", "neurons_microglia_sct.qs"), nthreads = ncores)
expr_matrix <- GetAssayData(so, slot = "data", assay = "SCT")

psdb_matrix <- AverageExpression(so, group.by = c("stable_clusters_26", "condition"), slot = "data", assays = "SCT")$SCT
psdb_list <- list(
  "neurons" = psdb_matrix[names(gene_mapping_list[['neurons']]), colnames(psdb_matrix)[!startsWith(colnames(psdb_matrix), "g22")]],
  "microglia" = psdb_matrix[names(gene_mapping_list[['microglia']]), colnames(psdb_matrix)[startsWith(colnames(psdb_matrix), "g22")]],
  "all" = psdb_matrix[c(names(gene_mapping_list[['neurons']]), names(gene_mapping_list[["microglia"]])), colnames(psdb_matrix)]
)

edges = c(50, 75, 100, 150)

for (n_edges in edges) {
    for (condition in names(psdb_list)) {


      set.seed(2023)
      res <- GENIE3::GENIE3(as.matrix(psdb_list[[condition]]), targets = rownames(psdb_list[[condition]]),
            nCores = 16)

      g <- plot_GRN_igraph(res, n_edges)
      if (any(igraph::E(g)$weight > 0)) {
        min_weight <- 0.5
        max_weight <- 7
        # rescale edge weights to a range
        igraph::E(g)$weight <- scales::rescale(igraph::E(g)$weight, to = c(min_weight, max_weight))
        pdf(file.path(figure_path, paste(condition, n_edges, '.pdf', sep='')), width = 10, height = 10)
        plot(g,
            edge.width = igraph::E(g)$weight,
            edge.color = "black",
            vertex.size = 5,
            vertex.label.cex = 1.0,
            vertex.label.dist = 1.0,
            vertex.label.degree = -pi/2,
            vertex.color = "#8bc0e0",
            vertex.frame.color = "black",
            vertex.label.color = "black",
            layout = igraph::layout_with_fr,
            main = paste(condition, n_edges, sep = " - ")
          )
        dev.off()
      }

      for (unq_cond in unique(as.character(so$condition))) {
        used_mat <- as.matrix(psdb_list[[condition]][, endsWith(colnames(psdb_list[[condition]]), unq_cond)])
        res <- GENIE3::GENIE3(used_mat, targets = rownames(used_mat), nCores = 16)
        g <- plot_GRN_igraph(res, n_edges)
        if (all(igraph::E(g)$weight == 0)) {
          next
        }
        min_weight <- 0.5
        max_weight <- 7
        # rescale edge weights to a range
        igraph::E(g)$weight <- scales::rescale(igraph::E(g)$weight, to = c(min_weight, max_weight))
        pdf(file.path(figure_path, paste(condition, unq_cond, n_edges, '.pdf', sep='')), width = 10, height = 10)
        plot(g,
            edge.width = igraph::E(g)$weight, 
            edge.color = "black",
            vertex.size = 5,
            vertex.label.cex = 1.0,
            vertex.label.dist = 1.0,
            vertex.label.degree = -pi/2,
            vertex.color = "#8bc0e0",
            vertex.frame.color = "black",
            vertex.label.color = "black",

            layout = igraph::layout_with_fr,
            main = paste(condition, unq_cond, n_edges, sep = " - ")
          )
        dev.off()
      }
    }
}
