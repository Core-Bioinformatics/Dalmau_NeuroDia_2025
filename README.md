# Metabolic challenges regulates hypothalamic SUCNR1 expression
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22976793.svg)](https://doi.org/10.5281/zenodo.22976793)

GEO accession number: GSE299826 


Scripts Summary
- `01_cellranger.sh` - Runs Cell Ranger count for each sample and moves the results into the preprocessing directory.
- `02_create_seurat_from_cr.R` - Creates and merges Seurat objects from Cell Ranger output, generates QC plots, and filters cells.
- `03_preprocess_seurat.R` - Applies SCTransform, cell-cycle scoring, PCA, and UMAP to each configured Seurat object.
- `04_run_clustassess.R` - Runs ClustAssess stability analyses across feature sets, neighbour counts, and clustering resolutions.
- `05_update_seurat_mtd.R` - Adds selected stable-cluster metadata and embeddings from ClustAssess results to Seurat objects.
- `06_write_shiny_cell_app.R` - Generates ShinyCell Shiny applications for interactive exploration of the processed Seurat objects.
- `07_annotate_object.R` - Applies configured cell-type annotations to Seurat metadata and updates the ClustAssess applications.
- `08_subset_by_metadata.R` - Creates and preprocesses Seurat object subsets based on configured metadata values.
- `09_generate_grn.R` - Generates and plots gene-regulatory networks for neuron and microglia marker gene sets.
- `10_cellchat.R` - Infers cell-cell communication networks with CellChat and creates interactive CellChat applications.
- `figures_aggregate.R` - Generates figures for the aggregate dataset, including UMAPs and expression summaries.
- `figures_neurons_microglia.R` - Generates figures for the neurons-plus-microglia dataset, including UMAPs and expression summaries.
