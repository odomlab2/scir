#' @title Import STARSolo results into monocle3.
#' @description Import one or multiple sci-RNA-seqv3 STARSolo (filtered) matrices into
#' a combined cds object with correct gene-annotations.
#'
#' @param folder (character): Folder containing the results of STAR / STARSolo (GeneFull_Ex50pAS/filtered).
#' @param samples (character): Samples to import (will search for files).
#' @param gtf (character): Path to GTF file used during alignment and for annotation purposes.
#' @param metadata (tibble): Additional metadata to add to cell_data_set. Should contain a sample_name column.
#'
#' @examples
#' 1 + 1
#'
#' @return (cell_data_set): Cell_data_set containing all combined samples and annotations without further processing.
#'
#' @importFrom dplyr %>%
#' @export
import_samples_monocle <- function(folder, samples, gtf, metadata = NULL) {
    # Input validation ----
    checkmate::assertAccess(folder, access = "r")
    checkmate::assertCharacter(samples)
    checkmate::assertFileExists(gtf, access = "r")
    checkmate::assertTibble(metadata, null.ok = TRUE)

    if (!is.null(metadata) & is.null(metadata$sample_name)) stop("metadata should contain a column: sample_name")

    # Retrieve  the required files. ----

    futile.logger::flog.info("import_samples_monocle: Retrieving files.")

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

    futile.logger::flog.info("import_samples_monocle: Retrieving and processing GTF file.")

    data_gtf <- rtracklayer::import(gtf, format = "gff") %>%
        tibble::as_tibble() %>%
        dplyr::filter(
            type == "gene",
            !base::is.na(gene_id)
        ) %>%
        # Subset on classes-of-interest.
        dplyr::filter(
            base::grepl("protein_coding|lncRNA|IG_.*_gene", gene_type)
        ) %>%
        # Remove predicted genes.
        dplyr::filter(
            !base::grepl("^Gm[0-9]", gene_name) | gene_name == 'Gm2a'
        ) %>%
        dplyr::distinct(
            seqnames, start, end, strand, gene_id, gene_type, gene_name
        ) %>%
        # Change names in order not to conflict with GRanges.
        dplyr::select(
            gene_chr = seqnames, gene_start = start, gene_end = end, gene_strand = strand, gene_id, gene_type, gene_name
        )


    # Generate a cell_data_set per sample ----

    futile.logger::flog.info("import_samples_monocle: Generating cell_data_set object for %s samples", dplyr::n_distinct(files$sample))

    p <- progressr::progressor(along = unique(files$sample))
    cds_samples <- future.apply::future_lapply(unique(files$sample), function(current_sample) {
        # Generate initial cell_data_set.
        cds <- monocle3::load_mm_data(
            mat_path = files %>% dplyr::filter(sample == current_sample, file == "matrix.mtx") %>% dplyr::pull(path),
            feature_anno_path = files %>% dplyr::filter(sample == current_sample, file == "features.tsv") %>% dplyr::pull(path),
            cell_anno_path = files %>% dplyr::filter(sample == current_sample, file == "barcodes_converted.tsv") %>% dplyr::pull(path),
            feature_metadata_column_names = c("gene_short_name", "type"),
            # Use all cells based on the STARSolo threshold.
            umi_cutoff = 0
        )

        # Subset on genes.
        cds <- cds[SummarizedExperiment::rowData(cds)$gene_short_name %in% data_gtf$gene_name, ]

        # Add additional metadata from the GTF.
        monocle3::fData(cds) <- monocle3::fData(cds) %>%
            tibble::as_tibble(rownames = "gene_id") %>%
            dplyr::left_join(data_gtf, by = "gene_id") %>%
            tibble::column_to_rownames("gene_id") %>%
            S4Vectors::DataFrame()

        # Add the sample_name as metadata.
        SummarizedExperiment::colData(cds)$sample_name <- base::factor(current_sample)

        # Visualize progress bar.
        p(base::sprintf("%s", current_sample))

        return(cds)
    }, future.seed = TRUE)

    # Add the names to the list.
    base::names(cds_samples) <- base::vapply(cds_samples, FUN = function(x) {
        base::levels(monocle3::pData(x)$sample_name)
    }, "")

    # Combine samples.
    futile.logger::flog.info("import_samples_monocle: Combining into a single cell_data_set.")
    cds_combined <- monocle3::combine_cds(cds_samples, keep_all_genes = FALSE)
    SummarizedExperiment::colData(cds_combined)$sample_name <- NULL

    # Add metadata.
    if (!is.null(metadata)) {
        futile.logger::flog.info("import_samples_monocle: Adding metadata.")

        # Add the metadata to the cds object.
        monocle3::pData(cds_combined) <- monocle3::pData(cds_combined) %>%
            tibble::as_tibble(rownames = "cell_id") %>%
            dplyr::left_join(metadata, by = c("sample" = "sample_name")) %>%
            tibble::column_to_rownames("cell_id") %>%
            S4Vectors::DataFrame()
    }

    # Add the gene Id to fData as first column.
    monocle3::fData(cds_combined) <- S4Vectors::cbind(data.frame(id = base::gsub("\\..*", "", S4Vectors::rownames(cds_combined))), monocle3::fData(cds_combined))

    # Add log10 umi.
    cds_combined$log.n.umi <- log10(cds_combined$n.umi)

    # Return combined samples.
    return(cds_combined)
}
