#' @title Generates knee-plot for one or multiple samples based on STARSolo.
#' @description Imports the No. of UMI per barcode and generates a knee-plot.
#' The selected threshol of STARSolo will be visualized, an alternative threshold
#' can be specified to visualize this alternative compared to current threshold.
#'
#' @param folder (character): Folder containing the results of STAR / STARSolo.
#' @param samples (character): Samples to import (will search for files).
#' @param umi_threshold (numeric): Custom threshold to visualize next to STARSolo threshold.
#'
#' @examples
#' 1 + 1
#'
#' @return (list): List of ggplot2 objects.
#'
#' @importFrom dplyr %>%
#' @import future
#' @import progressr
#' @import ggplot2
#' @export
plot_kneeplot <- function(folder, samples, umi_threshold = NULL) {

    # Input validation ----
    checkmate::assertAccess(folder, access = "r")
    checkmate::assertCharacter(samples)
    checkmate::assertNumeric(umi_threshold, null.ok = TRUE)

    futile.logger::flog.info("plot_kneeplot: Retrieving required files.")

    files <- tibble::tibble(
        path = list.files(folder, pattern = "UMIperCellSorted.txt$|matrix.mtx$", recursive = TRUE, full.names = TRUE)
    )

    files <- files %>%
        # Only use the Filtered data.
        dplyr::filter(base::grepl("GeneFull/filtered/", .data$path) | base::grepl("UMIperCellSorted.txt", .data$path)) %>%
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

    futile.logger::flog.info(glue::glue("plot_kneeplot: Generating for {dplyr::n_distinct(files$sample)} samples"))

    # For each sample, generate knee-plot by reading the filtered UMI counts and raw matrix (to determine STARSolo threshold.)
    p <- progressr::progressor(along = unique(files$sample))
    knee_samples <- future.apply::future_lapply(unique(files$sample), function(current_sample) {
        futile.logger::flog.debug(glue::glue("plot_kneeplot: Working on {current_sample}"))

        # Retrieve the required paths.
        path_matrix <- files %>%
            dplyr::filter(sample == current_sample, file == "matrix.mtx") %>%
            dplyr::pull(.data$path)
        path_umi <- files %>%
            dplyr::filter(sample == current_sample, file == "UMIperCellSorted.txt") %>%
            dplyr::pull(.data$path)

        # Determine STARSolo cutoff. ----
        star_cutoff <- Matrix::readMM(file = path_matrix) %>%
            Matrix::colSums(.) %>%
            base::min(.)

        # Retrieve total UMI count per cell. ----
        data_umi <- readr::read_tsv(path_umi, col_names = "n_umi", col_types = "i") %>%
            # Add UMI rank.
            dplyr::mutate(rank_umi = dplyr::min_rank(- n_umi))

        # Generate knee-plot ----
        plot_knee <- ggplot2::ggplot(data_umi, ggplot2::aes(x = rank_umi, y = n_umi)) +
            ggplot2::geom_line(lwd = 1) +
            ggplot2::scale_x_log10() +
            ggplot2::scale_y_log10() +
            ggplot2::annotation_logticks(sides = "lb", scaled = TRUE) +
            ggplot2::labs(x = "Barcodes<br><sub>(Ranked on no. of UMI)</sub>", y = "No. of UMI") +
            ggplot2::geom_hline(yintercept = star_cutoff, lwd = .5, lty = 11, color = "red") +
            ggplot2::annotate("text", x = 100, y = star_cutoff * .75, size = 3, label = glue::glue("Threshold: {star_cutoff}"), color = "red") +
            ggplot2::annotate("text", x = Inf, y = Inf, size = 3, vjust = 1, hjust = 1,
                              label = glue::glue("Sample: {current_sample}\nTotal barcodes: {dim(data_umi)[1]}\nAfter filtering: {sum(data_umi$n_umi >= star_cutoff)}"))

        if (!is.null(umi_threshold)) {
            plot_knee <- plot_knee +
                ggplot2::geom_hline(yintercept = umi_threshold, lwd = .5, lty = 11, color = "darkblue") +
                ggplot2::annotate("text", x = 100, y = umi_threshold * .75, size = 3, label = glue::glue("Threshold (manual): {umi_threshold}"), color = "darkblue") +
                ggplot2::annotate("text", x = Inf, y = Inf, size = 3, vjust = 1, hjust = 1,
                                  label = glue::glue("Sample: {current_sample}\nTotal barcodes: {dim(data_umi)[1]}\nAfter filtering: {sum(data_umi$n_umi >= star_cutoff)}\nAfter filtering (Custom): {sum(data_umi$n_umi >= umi_threshold)}"))
        }

        # Add the ggplot theme.
        plot_knee <- plot_knee + scir::theme_ggplot()

        # Visualize progress bar.
        p(base::sprintf("%s", current_sample))

        # Return plot.
        return(plot_knee)
    }, future.seed = TRUE)

    return(knee_samples)
}
