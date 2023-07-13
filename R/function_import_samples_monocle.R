#' @title import_samples_monocle
#' @description Import one or multiple sci-RNA-seqv3 samples into a monocle3 object.
#'
#' @param folder (character): Folder containing the results of STAR / STARSolo.
#' @param samples (character): Samples to import (will search for files).
#' @param gtf (character): Path to GTF file used during alignment and for annotation purposes.
#'
#' @examples
#' 1 + 1
#'
#' @return (cell_data_set): Cell_data_set containing all combined samples and annotations without further processing.
#'
#' @importFrom dplyr %>%
#' @export
import_samples_monocle <- function(folder, samples, gtf) {
    # Input validation ----
    checkmate::assertAccess(folder, access = "r")
    checkmate::assertCharacter(samples)
    checkmate::assertFileExists(gtf, access = "r")

    # Retrieve  the required files. ----

    futile.logger::flog.info("import_samples_monocle - Retrieving files.")

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

    futile.logger::flog.info("import_samples_monocle - Retrieving and processing GTF file.")

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

    cds_samples <- pbapply::pblapply(unique(files$sample), function(current_sample) {
        futile.logger::flog.info("\t%s", current_sample)

        tryCatch({

            # Generate initial cell_data_set.
            cds <- monocle3::load_mm_data(
                mat_path = files %>% dplyr::filter(sample == current_sample, file == "matrix.mtx") %>% dplyr::pull(.data$path),
                feature_anno_path = files %>% dplyr::filter(sample == current_sample, file == "features.tsv") %>% dplyr::pull(.data$path),
                cell_anno_path = files %>% dplyr::filter(sample == current_sample, file == "barcodes_converted.tsv") %>% dplyr::pull(.data$path),
                feature_metadata_column_names = c("gene_short_name", "type")
            )


            # Add additional metadata from the GTF.
            SummarizedExperiment::rowData(cds) <- SummarizedExperiment::rowData(cds) %>%
                tibble::as_tibble(rownames = "gene_id") %>%
                dplyr::left_join(data_gtf, by = "gene_id") %>%
                tibble::column_to_rownames("gene_id") %>%
                S4Vectors::DataFrame()

            # Add the sample_name as metadata.
            SummarizedExperiment::colData(cds)$sample_name <- base::factor(current_sample)

            return(cds)

        },
        error=function(cond) {
            return(NULL)
        })

    }, cl = 10)

    # Remove the NULL objects.
    cds_samples <- plyr::compact(cds_samples)

    # Add the names to the list.
    base::names(cds_samples) <- base::vapply(cds_samples, FUN = function(x) {
        base::levels(SummarizedExperiment::colData(x)$sample_name)
    }, "")

    # Combine samples.
    cds_combined <- monocle3::combine_cds(cds_samples)
    SummarizedExperiment::colData(cds_combined)$sample_name <- NULL

    # Return combined samples.
    return(cds_combined)
}
