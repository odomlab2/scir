#' @title Generate Sankey cluster-resolution plot.
#' @description Generates a Sankey plot to visualize the size and identity of clusters across different resolutions.
#'
#' @param seurat (Seurat): Seurat object.
#' @param cluster_prefix (str): Prefix of the cluster column in the Seurat object.
#' @param min_transition (numeric): Min. number of cells in cluster -> cluster transition to include in the plot.
#'
#' @examples
#' \dontrun{
#' cluster_sankey(seurat, cluster_prefix = "SCT_snn_res.", min_transition = 1000)
#' }
#'
#' @return (plotly): Plotly Sankey plot.
#'
#' @importFrom dplyr %>%
#' @export
cluster_sankey <- function(seurat, cluster_prefix, min_transition = 500){

    # Input validation ----
    checkmate::assertClass(seurat, "Seurat")
    checkmate::assertCharacter(cluster_prefix)
    checkmate::assertNumber(min_transition, null.ok = FALSE)

    # Check if the cluster column exists.
    if (!any(grepl(cluster_prefix, colnames(seurat@metadata)))) stop(glue::glue("Cluster column {cluster_prefix} not found in Seurat object."))

    # Retrieve the clustering columns.
    data_clusters <- data_clusters %>%
        dplyr::select(starts_with(cluster_prefix)) %>%
        dplyr::mutate(cell = rownames(.)) %>%
        tidyr::gather(key = "resolution", value = "cluster", -cell) %>%
        dplyr::mutate(resolution = gsub(cluster_prefix, "", resolution))

    # Determine the transition of each cell between clusters for plotly.
    data_clusters <- data_clusters %>%
        dplyr::arrange(cell, resolution) %>%
        dplyr::mutate(
            cluster_next = dplyr::lead(cluster, default = "NA"),
            resolution_next = dplyr::lead(resolution, default = "NA"),
            source = sprintf("Cluster %s\n(%s)", cluster, resolution),
            target = sprintf("Cluster %s\n(%s)", cluster_next, resolution_next),
            label = sprintf("%s -> %s", resolution, resolution_next)
        ) %>%
        # Only keep transitions to higher resolution.
        dplyr::filter(resolution < resolution_next) %>%
        dplyr::group_by(source, target, label) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::ungroup() %>%
        dplyr::filter(n >= min_transition)

    # Generate the Sankey plot.
    p <- plotly::plot_ly(
        data = data_clusters,
        type = "sankey",
        orientation = "h",
        node = list(
            pad = 5,
            thickness = 3,
            line = list( color = "black", width = .5 ),
            label = data_clusters$source %>% unique() %>% sort(),
            color = '#2197FA'
        ),
        link = list(
            source = match(data_clusters$source, data_clusters$source %>% unique() %>% sort()) -1,
            target = match(data_clusters$target, data_clusters$source %>% unique() %>% sort())-1,
            value = data_clusters$n
        ),
        arrangement = "snap"
    ) %>%
        plotly::layout(
            font = list(size = 8)
        )

    return(p)
}
