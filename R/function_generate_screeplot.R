#' @title Generates scree-plot for monocle3 object.
#'
#' @param cds (cell_data_set): cell_data_set with preprocess_cds performed.
#' @param pc_total (integer): Total principal components to plot.
#' @param pc_chosen (integer): Principal component to draw a horizontal line.
#' @param pc_threshold (numeric): Custom threshold for min. variance explained (vertical line).
#'
#' @examples
#' 1 + 1
#'
#' @return (ggplot2): ggplot2 object.
#'
#' @import ggplot2
#' @export
plot_screeplot <- function(cds, pc_total = 75, pc_chosen = NULL, pc_threshold = 0.1) {
  # Input validation ----
  checkmate::assertClass(cds, "cell_data_set")
  checkmate::assertNumber(pc_total, null.ok = FALSE)
  checkmate::assertNumber(pc_chosen, null.ok = TRUE)
  checkmate::assertNumeric(pc_threshold, null.ok = TRUE)

  # Check
  cds_pc <- base::length(cds@reduce_dim_aux[["PCA"]][["model"]][["prop_var_expl"]])
  if (cds_pc < pc_total) stop(glue::glue("Requesting too many PCs to plot, there are {cds_pc} PCs available."))

  # Plot.
  plot <- data.frame(
    PC = 1:cds_pc,
    var = cds@reduce_dim_aux[["PCA"]][["model"]][["prop_var_expl"]]
  ) %>%
    dplyr::filter(PC <= pc_total) %>%
    ggplot2::ggplot(., ggplot2::aes(x = PC, y = var))

  if (!is.null(pc_chosen)) plot <- plot + ggplot2::geom_vline(xintercept = 20, lty = 11, color = "red")
  if (!is.null(pc_threshold)) plot <- plot + ggplot2::geom_hline(yintercept = pc_threshold, lty = 11, color = "grey50")

  plot <- plot +
    ggplot2::geom_bar(stat = "identity", width = .2, fill = "black", lwd = 0, na.rm = FALSE) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(.01, .01)), breaks = c(1, seq(5, pc_total, by = 5))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .1)), labels = scales::percent_format(), breaks = c(pc_threshold, seq(0, 1, by = .1))) +
    ggplot2::labs(x = "Principal components", y = "Variance explained") +
    scir::theme_ggplot()

  # Return plot.
  return(plot)
}
