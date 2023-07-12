

x <- scir::import_samples_monocle(
    folder = "/Users/jobbie/DKFZ/odomLab/sexomics/runJob/sx42/alignment/",
    samples  = c("sx39_1", "sx39_2"),
    gtf = "/Users/jobbie/Downloads/gencode.vM32.basic.annotation.gtf.gz"
)

# Retrieve the prootein coding genes.
genes <- base::subset(rowData(x), gene_type == 'protein_coding')

sx42 <- x %>%
    monocle3::preprocess_cds(method = "PCA", use_genes = rownames(genes)) %>%
    monocle3::align_cds(alignment_group = "sample") %>%
    monocle3::reduce_dimension(., reduction_method = "UMAP") %>%
    monocle3::cluster_cells(., cluster_method = "louvain")

monocle3::plot_cells(sx42, color_cells_by = "sample", cell_size = 1, cell_stroke = .2)
