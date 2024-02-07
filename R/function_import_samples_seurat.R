#' @title Import STARSolo results into Seurat.
#' @description Import one or multiple sci-RNA-seqv3 STARSolo (filtered) matrices into
#' a combined Seurat object with correct gene-annotations.
#'
#' @param folder (character): Folder containing the results of STAR / STARSolo (GeneFull_Ex50pAS/filtered).
#' @param samples (character): Samples to import (will search for files).
#' @param gtf (character): Path to GTF file used during alignment and for annotation purposes.
#' @param metadata (tibble): Additional metadata to add to cell_data_set. Should contain a sample_name column.
#'
#' @examples
#' 1 + 1
#'
#' @return (Seurat): Seurat containing all combined samples and annotations without further processing.
#'
#' @importFrom dplyr %>%
#' @export
import_samples_seurat <- function(folder, samples, gtf, metadata = NULL) {
    # Input validation ----
    checkmate::assertAccess(folder, access = "r")
    checkmate::assertCharacter(samples)
    checkmate::assertFileExists(gtf, access = "r")
    checkmate::assertTibble(metadata, null.ok = TRUE)

    if (!is.null(metadata) & is.null(metadata$sample_name)) stop("metadata should contain a column: sample_name")

    # Retrieve  the required files. ----

    futile.logger::flog.info("import_samples_seurat: Retrieving files.")

    files <- tibble::tibble(
        path = list.files(folder, pattern = "barcodes_converted.tsv$|features.tsv$|matrix.mtx$", recursive = TRUE, full.names = TRUE)
    )

    files <- files %>%
        # Only use the GeneFull_Ex50pAS/Filtered data.
        dplyr::filter(base::grepl("GeneFull_Ex50pAS/filtered/", path)) %>%
        # Retrieve features from files.
        dplyr::mutate(
            # Retrieve the sample name from the path.
            sample = base::sub("_Solo.out.*", "", base::sub("/$|^/", "", base::sub(folder, "", path))),
            # Remove everything after the last underscore.
            sample = base::sub("_[^_]*$", "", sample),
            file = base::basename(path)
        ) %>%
        # Filter the specified samples.
        dplyr::filter(sample %in% samples)

    # Retrieve the GTF file. ----

    futile.logger::flog.info("import_samples_seurat: Retrieving and processing GTF file.")

    data_gtf <- rtracklayer::import(gtf, format = "gff") %>%
        tibble::as_tibble() %>%
        dplyr::filter(
            type == "gene",
            !base::is.na(gene_id)
        )

    # Subset on classes-of-interest.
    if ("gene_type" %in% colnames(data_gtf)){
        data_gtf <- data_gtf %>% dplyr::filter(
            base::grepl("protein_coding|lncRNA|IG_.*_gene", gene_type)
        )
    }

    # Mouse-specific filtering.
    if ("mgi_id" %in% colnames(data_gtf)){
        # Remove predicted genes.
        data_gtf <- data_gtf %>%
            dplyr::filter(
                !base::grepl("^Gm[0-9]", gene_name) | gene_name == 'Gm2a'
            ) %>%
            # Remove RIKEN lncRNA genes.
            dplyr::filter(
                ! (base::grepl("Rik$", gene_name) & gene_type == 'lncRNA')
            ) %>%
            # Remove genes without a gene-name.
            dplyr::filter(!grepl("ENSMUS", gene_name))
    }

    # Clean-up.
    data_gtf <- data_gtf %>%
        dplyr::distinct(
            seqnames, start, end, strand, gene_id, gene_type, gene_name
        ) %>%
        # Change names in order not to conflict with GRanges.
        dplyr::select(
            gene_chr = seqnames, gene_start = start, gene_end = end, gene_strand = strand, gene_id, gene_type, gene_name
        ) %>%
        dplyr::filter(!duplicated(gene_name))

    # Generate a seurat object per sample ----

    futile.logger::flog.info("import_samples_seurat: Generating Seurat object for %s samples", dplyr::n_distinct(files$sample))

    p <- progressr::progressor(along = unique(files$sample))
    seurat_samples <- future.apply::future_lapply(unique(files$sample), function(current_sample) {
        # Read STARSolo counts.
        mtx <- Seurat::ReadMtx(
            mtx = files %>% dplyr::filter(sample == current_sample, file == "matrix.mtx") %>% dplyr::pull(path),
            features = files %>% dplyr::filter(sample == current_sample, file == "features.tsv") %>% dplyr::pull(path),
            cells = files %>% dplyr::filter(sample == current_sample, file == "barcodes_converted.tsv") %>% dplyr::pull(path),
        )

        # Create Seurat object of single sample.
        x <- Seurat::CreateSeuratObject(counts = mtx, min.cells = 0, min.features = 50)

        # Subset on selected genes.
        x <- x[base::rownames(x) %in% data_gtf$gene_name, ]

        # Add optional metadata.
        if (!is.null(metadata)) {
            x[[]] <- cbind(x[[]], metadata %>% dplyr::filter(sample_name == current_sample))
        }

        # Visualize progress bar.
        p(base::sprintf("%s", current_sample))

        return(x)
    }, future.seed = TRUE)

    # Merge multiple samples. ----
    seurat_combined <- merge(seurat_samples[[1]], y = seurat_samples[2:length(seurat_samples)])
    seurat_combined <- SeuratObject::JoinLayers(seurat_combined,  assay = "RNA")

    # Add GTF information.
    data_gtf_overlap <- data.frame('gene_name' = rownames(seurat_combined)) %>%
        dplyr::left_join(data_gtf, by = 'gene_name')

    seurat_combined[["RNA"]][[]] <- data_gtf_overlap

    # Return combined samples.
    return(seurat_combined)
}
