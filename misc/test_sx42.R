library(dplyr)
library(scir)

# Import sample-sheet. ----
metadata <- readr::read_tsv("/omics/groups/OE0538/internal/projects/sexomics/metadata_scirocket/sx42b_samplesheet.tsv")

# Read STARSolo output as a combined cell_data_set. ----
sx42 <- scir::import_samples_monocle(
    folder = "/omics/groups/OE0538/internal/projects/sexomics/runJob/sx42b/alignment/",
    samples  = metadata$sample_name,
    gtf = "/omics/groups/OE0538/internal/projects/sharedData/GRCm39/annotation/gencode.vM31.basic.annotation.gtf"
)

# Run default monocle3 workflow. ----
sx42 <- sx42 %>%
    monocle3::preprocess_cds(method = "PCA") %>%
    monocle3::align_cds(alignment_group = "sample", preprocess_method = "PCA") %>%
    monocle3::reduce_dimension(., reduction_method = "UMAP", cores = 20) %>%
    monocle3::cluster_cells(., cluster_method = "leiden")

# Add additional metadata.
meta_samples <- readr::read_csv("misc/sx42_metadata.csv")
sx42$group <- merge(SummarizedExperiment::colData(sx42), meta_samples)$Group

## Save object. ----
monocle3::save_monocle_objects(cds=sx42, directory_path='/omics/groups/OE0538/internal/projects/sexomics/runJob/sx42b/rdata/sx42_cds')

# Plot UMAP ----
plot_umap <- monocle3::plot_cells(
    cds = sx42,
    reduction_method = "UMAP",
    label_cell_groups = FALSE,
    color_cells_by = "group",
    cell_size = .2,
    rasterize = TRUE,
    cell_stroke = .2,
    show_trajectory_graph = FALSE,
) +
    ggplot2::labs(x = "UMAP - Dim. 1", y = "UMAP - Dim. 2") +
    ggplot2::scale_color_manual(values = c('grey50', '#E5193A', '#0AB066', '#526AB0', '#FE9639')) +
    scir::theme_ggplot()

plot_umap + ggplot2::facet_grid(~group)
