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

# Add additional metadata.
meta_samples <- readr::read_csv("sx42_metadata.csv")
sx42$group <- merge(SummarizedExperiment::colData(sx42), meta_samples)$Group

sx42_nonorm <- sx42 %>%
    monocle3::preprocess_cds(method = "PCA") %>%
    monocle3::reduce_dimension(., reduction_method = "UMAP", cores = 20) %>%
    monocle3::cluster_cells(., cluster_method = "leiden")

sx42_normalized <- sx42 %>%
    monocle3::preprocess_cds(method = "PCA") %>%
    monocle3::align_cds(alignment_group = "sample", preprocess_method = "PCA") %>%
    monocle3::reduce_dimension(., reduction_method = "UMAP", cores = 20) %>%
    monocle3::cluster_cells(., cluster_method = "leiden")


## Save object. ----
monocle3::save_monocle_objects(cds=sx42_nonorm, directory_path='/omics/groups/OE0538/internal/projects/sexomics/runJob/sx42b/rdata/sx42_nonorm')
monocle3::save_monocle_objects(cds=sx42_normalized, directory_path='/omics/groups/OE0538/internal/projects/sexomics/runJob/sx42b/rdata/sx42_normalized')

## Load object. ----

sx42_nonorm <- monocle3::load_monocle_objects('/omics/groups/OE0538/internal/projects/sexomics/runJob/sx42b/rdata/sx42_nonorm')

# Plot UMAP ----
plot_umap <- monocle3::plot_cells(
    cds = sx42_nonorm,
    reduction_method = "UMAP",
    label_cell_groups = FALSE,
    color_cells_by = "group",
    cell_size = .4,
    rasterize = FALSE,
    cell_stroke = .1, alpha = .1,
    show_trajectory_graph = FALSE,
) +
    ggplot2::labs(x = "UMAP - Dim. 1", y = "UMAP - Dim. 2") +
    ggplot2::scale_color_manual(values = c('grey50', '#E5193A', '#0AB066', '#526AB0', '#FE9639')) +
    scir::theme_ggplot()

ggplot2::ggsave(filename = '~/test/umap_sx42.pdf', plot = plot_umap + ggplot2::facet_wrap(~group), units="px", width = 3000, height = 1500, device='pdf', dpi=200)
