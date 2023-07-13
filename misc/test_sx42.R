# Themes ----
theme_job <- ggplot2::theme(
    text = ggplot2::element_text(family = "Roboto"),
    axis.text = ggplot2::element_text(size = 9),
    axis.title.x = ggtext::element_markdown(size = 12, face = "bold"),
    axis.title.y = ggtext::element_markdown(size = 12, face = "bold"),
    legend.title = ggplot2::element_text(face = "bold"),
    legend.text = ggplot2::element_text(size = 10),
    legend.position = "bottom",
    legend.key.size = ggplot2::unit(0.1, "cm"),
    legend.key.width = ggplot2::unit(0.1, "cm"),
    legend.key.height = ggplot2::unit(0.1, "cm"),
    legend.key = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0),
    strip.background = ggplot2::element_rect(fill = "grey50"),
    strip.text = ggplot2::element_text(color = "white"),
    panel.grid.major.y = ggplot2::element_line(color = "grey75", linetype = "dotted", linewidth = ggplot2::rel(.75)),
    panel.grid.major.x = ggplot2::element_line(color = "grey75", linetype = "dotted", linewidth = ggplot2::rel(.75)),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(colour = "white"),
    plot.margin = ggplot2::unit(c(10, 5, 5, 5), "mm"),
    axis.line = ggplot2::element_line(colour = "black", linewidth = ggplot2::rel(1))
)


# Import sample-sheet.
metadata <- readr::read_tsv("/omics/groups/OE0538/internal/projects/sexomics/metadata_scirocket/sx42b_samplesheet.tsv")

sx42 <- import_samples_monocle(
    folder = "/omics/groups/OE0538/internal/projects/sexomics/runJob/sx42b/alignment/",
    samples  = metadata$sample_name,
    gtf = "/omics/groups/OE0538/internal/projects/sharedData/GRCm39/annotation/gencode.vM31.basic.annotation.gtf"
)

# Retrieve the protein-coding genes.
genes <- base::subset(SummarizedExperiment::rowData(x), gene_type == 'protein_coding')
subset(SummarizedExperiment::rowData(x), gene_short_name == "Xist")

sx42 <- sx42 %>%
    monocle3::preprocess_cds(method = "PCA") %>%
    monocle3::align_cds(alignment_group = "sample") %>%
    monocle3::reduce_dimension(., reduction_method = "UMAP") %>%
    monocle3::cluster_cells(., cluster_method = "leiden")

# Add additional metadata.
meta_samples <- readr::read_csv("misc/sx42_metadata.csv")
sx42$group <- merge(SummarizedExperiment::colData(sx42), meta_samples)$Group

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
    theme_job

plot_umap + ggplot2::facet_grid(~group)


monocle3::plot_cells(sx42, color_cells_by="partition", group_cells_by="partition")


