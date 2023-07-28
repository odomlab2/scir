#' @title Generates an overview of common scRNA QC-features.
#' @description This generates figures for the following QC-features:
#' - Relative mitochondrial expression.
#' - Number of expressed genes.
#' - Number of UMI.
#'
#' Within the expression plot, the mean and mean + 3 * sd are indicated.
#' @param cds (character): Folder containing the results of STAR / STARSolo.
#' @param limits_mt (numeric): Limits for the y-axis of the mitochondrial expression plot.
#' @param limits_expr (numeric): Limits for the y-axis of the number of expressed genes and UMI plot.
#'
#' @examples
#' 1 + 1
#'
#' @return (ggplot2): ggplot2 object.
#'
#' @importFrom dplyr %>%
#' @import ggplot2
#' @import patchwork
#' @export
plot_features <- function(cds, limits_mt = c(0, .02), limits_expr = c(0, 10000)) {
  # Input validation. ----
  checkmate::checkClass(cds, classes = "cell_data_set")

  # Determine rel. exprs. of MT-genes over total expression. ----

  futile.logger::flog.info("plot_features: Calculating rel. Mt content.")

  # Define mitochondrial genes.
  genes_mt <- base::which(monocle3::fData(cds)$gene_chr == "chrM")

  total_exprs_mt <- Matrix::colSums(monocle3::normalized_counts(cds[genes_mt], norm_method = "log"))
  total_exprs <- Matrix::colSums(monocle3::normalized_counts(cds, norm_method = "log"))
  monocle3::pData(cds)$rel_mt <- total_exprs_mt / total_exprs

  ## Visualize mitochondrial gene expression. ----
  futile.logger::flog.info("plot_features: Generating Mt-plot")
  plot_mt <- monocle3::pData(cds) %>%
    tibble::as_tibble() %>%
    ggplot2::ggplot(., ggplot2::aes(x = "MT", y = rel_mt, fill = "MT")) +
    ggplot2::scale_y_continuous(limits = limits_mt, breaks = c(0, .01, 1, by = .01), labels = scales::percent_format()) +
    ggbeeswarm::geom_quasirandom(width = .4, shape = 21) +
    ggplot2::labs(x = "<sub>Features per cell</sub>", y = "Rel. mitochondrial expression<br><sub>(normalized)</sub>") +
    ggplot2::scale_fill_manual(values = c("#FEF2DD"), guide = "none") +
    scir::theme_ggplot()

  ## Visualize total UMI and no. of genes expressed. ----
  futile.logger::flog.info("plot_features: Generating expression + UMI plot")

  plot_expr <- monocle3::pData(cds) %>%
    tibble::as_tibble() %>%
    dplyr::select("No. expressed genes" = num_genes_expressed, "No. UMI" = n.umi) %>%
    tidyr::gather(key = "variable", value = "value") %>%
    ggplot2::ggplot(., ggplot2::aes(x = variable, y = value, fill = variable, group = variable)) +
    ggplot2::scale_y_continuous(labels = scales::comma_format(), limits = c(0, 10000)) +
    ggplot2::labs(x = "<sub>Features per cell</sub>", y = "Frequency") +
    ggplot2::scale_fill_manual(values = c("#F3BFB3", "#5EA9BE"), guide = "none") +
    scir::theme_ggplot() +
    ggbeeswarm::geom_quasirandom(width = .4, shape = 21) +
    # Plot mean + 3* sd(x).
    ggplot2::stat_summary(
      fun.data = function(x) {
        y <- base::mean(x) + 3 * stats::sd(x)
        base::data.frame(y = y[1], ymin = y[1], ymax = y[1])
      },
      geom = "errorbar", color = "red", width = .75, size = 1.25
    ) +
    ggplot2::stat_summary(
      fun.data = function(x) {
        y <- base::mean(x) + 3 * stats::sd(x)
        base::data.frame(y = y[1])
      },
      colour = "red", size = 3, geom = "text",
      ggplot2::aes(label = ggplot2::after_stat(round(y))),
      position = ggplot2::position_nudge(y = .2, x = -.5)
    ) +
    # Plot mean.
    ggplot2::stat_summary(
      fun.data = function(x) {
        y <- base::mean(x)
        base::data.frame(y = y[1], ymin = y[1], ymax = y[1])
      },
      geom = "errorbar", color = "black", width = .75, size = 1.25
    ) +
    ggplot2::stat_summary(
      fun = "mean", colour = "black", size = 3,
      geom = "text",
      ggplot2::aes(label = ggplot2::after_stat(round(y))),
      position = ggplot2::position_nudge(y = .2, x = -.5)
    )

  # Return combined plot. ----
  plot_mt + plot_expr + patchwork::plot_layout(widths = c(1, 2)) + patchwork::plot_annotation(tag_levels = "a")
}
