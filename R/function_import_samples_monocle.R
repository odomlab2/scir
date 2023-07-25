#' @title Import STARSolo results into monocle3.
#' @description Import one or multiple sci-RNA-seqv3 STARSolo (filtered) matrices into
#' a combined cds object with correct gene-annotations.
#'
#' @param folder (character): Folder containing the results of STAR / STARSolo (GeneFull/filtered).
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

    if(!is.null(metadata) & is.null(metadata$sample_name)) stop("metadata should contain a column: sample_name")

    # Retrieve  the required files. ----

    futile.logger::flog.info("import_samples_monocle: Retrieving files.")

    files <- tibble::tibble(
        path = list.files(folder, pattern = "barcodes_converted.tsv$|features.tsv$|matrix.mtx$", recursive = TRUE, full.names = TRUE)
    )

    files <- files %>%
        # Only use the GeneFull/Filtered data.
        dplyr::filter(base::grepl("GeneFull/filtered/", .data$path)) %>%
        # Retrieve features from files.
        dplyr::mutate(
            # Retrieve the sample name from the path.
            sample = base::sub("_Solo.out.*", "", base::sub("/$|^/", "", base::sub(folder, "", .data$path))),
            # Remove everything after the last underscore.
            sample = base::sub("_[^_]*$", "", .data$sample),
            file = base::basename(.data$path)
        ) %>%
        # Filter the specified samples.
        dplyr::filter(.data$sample %in% samples)

    # Retrieve the GTF file. ----

    futile.logger::flog.info("import_samples_monocle: Retrieving and processing GTF file.")

    data_gtf <- rtracklayer::import(gtf, format = "gff") %>%
        tibble::as_tibble() %>%
        dplyr::filter(
            .data$type == "gene",
            !base::is.na(.data$gene_id)
        ) %>%
        dplyr::distinct(
            .data$seqnames, .data$start, .data$end, .data$width, .data$strand, .data$gene_id, .data$gene_type, .data$gene_name
        )


    # Generate a cell_data_set per sample ----

    futile.logger::flog.info("import_samples_monocle: Generating cell_data_set object for %s samples", dplyr::n_distinct(files$sample))

    p <- progressr::progressor(along = unique(files$sample))
    cds_samples <- future.apply::future_lapply(unique(files$sample), function(current_sample) {

        # Generate initial cell_data_set.
        cds <- monocle3::load_mm_data(
            mat_path = files %>% dplyr::filter(sample == current_sample, file == "matrix.mtx") %>% dplyr::pull(.data$path),
            feature_anno_path = files %>% dplyr::filter(sample == current_sample, file == "features.tsv") %>% dplyr::pull(.data$path),
            cell_anno_path = files %>% dplyr::filter(sample == current_sample, file == "barcodes_converted.tsv") %>% dplyr::pull(.data$path),
            feature_metadata_column_names = c("gene_short_name", "type"),
            # Use all cells based on the STARSolo threshold.
            umi_cutoff = 0
        )

        # Add additional metadata from the GTF.
        SummarizedExperiment::rowData(cds) <- SummarizedExperiment::rowData(cds) %>%
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
        base::levels(SummarizedExperiment::colData(x)$sample_name)
    }, "")

    # Combine samples.
    futile.logger::flog.info("import_samples_monocle: Combining into a single cell_data_set.")
    cds_combined <- monocle3::combine_cds(cds_samples, keep_all_genes = FALSE)
    SummarizedExperiment::colData(cds_combined)$sample_name <- NULL

    # Add metadata.
    if(!is.null(metadata)){
        futile.logger::flog.info("import_samples_monocle: Adding metadata.")

        # Add the metadata to the cds object.
        monocle3::pData(cds) <- S4Vectors::merge(monocle3::pData(cds), metadata, by.x = 'sample', by.y = 'sample_name', all.x = TRUE)
        base::rownames(monocle3::pData(cds)) <- base::sprintf("%s_%s", monocle3::pData(cds)$cell, monocle3::pData(cds)$sample)
    }

    # Return combined samples.
    return(cds_combined)
}
