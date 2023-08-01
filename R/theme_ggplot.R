#' @title Common theme.
#' @description Common theme for ggplot2.
#'
#' @return (list): ggplot theme list.
#' @export
theme_ggplot <- function() {
  ggplot2::theme(
    text = ggplot2::element_text(),
    axis.text = ggtext::element_markdown(size = 8),
    axis.title.x = ggtext::element_markdown(size = 10, face = "bold"),
    axis.title.y = ggtext::element_markdown(size = 10, face = "bold"),
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
}
