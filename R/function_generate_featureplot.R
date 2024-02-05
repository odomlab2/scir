#' @title Generates an overview of common scRNA QC-features.
#' @description This generates figures for the following QC-features:
#' - Relative mitochondrial expression.
#' - Number of expressed genes.
#' - Number of UMI.
#'
#' Within the expression plot, the mean and mean + 3 * sd are indicated.
#' @param seurat (Seurat): Seurat object.
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
plot_features <- function(seurat) {
    # Input validation. ----
    checkmate::checkClass(seurat, classes = "Seurat")

    # Determine rel. exprs. of MT-genes over total expression. ----

    futile.logger::flog.info("plot_features: Calculating rel. Mt content.")

    # Define mitochondrial genes.
    genes_mt <- base::rownames(seurat[seurat[['RNA']][[]]$gene_chr == 'chrM'])
    seurat[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat, features = genes_mt) / 100

    ## Visualize mitochondrial gene expression. ----
    futile.logger::flog.info("plot_features: Generating Mt-plot")
    plot_mt <- tibble::as_tibble(seurat@meta.data) %>%
        ggplot2::ggplot(., ggplot2::aes(x = "MT", y = percent.mt, fill = "MT")) +
        ggplot2::scale_y_continuous(labels = scales::percent_format(), breaks = scales::pretty_breaks(), trans = scales::pseudo_log_trans()) +
        ggbeeswarm::geom_quasirandom(width = .4, shape = 21, size = .5) +
        ggplot2::labs(x = NULL, y = "Rel. MT expression") +
        ggplot2::scale_fill_manual(values = c("#FEF2DD"), guide = "none") +
        scir::theme_ggplot()

    ## Visualize total UMI and no. of genes expressed. ----
    futile.logger::flog.info("plot_features: Generating expression + UMI plot")

    plot_expr <- tibble::as_tibble(seurat@meta.data) %>%
        dplyr::select("No. expressed genes" = nFeature_RNA, "No. UMI" = nCount_RNA) %>%
        tidyr::gather(key = "variable", value = "value") %>%
        ggplot2::ggplot(., ggplot2::aes(x = variable, y = value, fill = variable, group = variable)) +
        ggplot2::scale_y_continuous(labels = scales::comma_format()) +
        ggplot2::labs(x = "<sub>Features per cell</sub>", y = "Frequency") +
        ggplot2::scale_fill_manual(values = c("#F3BFB3", "#5EA9BE"), guide = "none") +
        scir::theme_ggplot() +
        ggbeeswarm::geom_quasirandom(width = .4, shape = 21, size = .5) +
        # Plot mean + 3* sd(x).
        ggplot2::stat_summary(
            fun.data = function(x) {
                y <- base::mean(x) + 3 * stats::sd(x)
                base::data.frame(y = y[1], ymin = y[1], ymax = y[1])
            },
            geom = "errorbar", color = "red", width = .75, lwd = 1.25
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
            geom = "errorbar", color = "black", width = .75, lwd = 1.25
        ) +
        ggplot2::stat_summary(
            fun = "mean", colour = "black", size = 3,
            geom = "text",
            ggplot2::aes(label = ggplot2::after_stat(round(y))),
            position = ggplot2::position_nudge(y = .2, x = -.5)
        )

    # Return combined plot. ----
    patchwork::wrap_plots(plot_mt, plot_expr, ncol = 2, widths = c(1, 3)) + patchwork::plot_annotation(tag_levels = "a")
}
