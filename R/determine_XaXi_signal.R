#' @title Determine Xa/Xi status from haplotyping data.
#' @description Determines cell-wise Xa/Xi status and Xa/Xi status of chromosome X genes.
#'
#' @param seurat (Seurat): Seurat object.
#' @param path_counts (str): Path to the H1/H2/UA counts table.
#' @param min_counts_cells (numeric): Min. number of H1 or H2 counts needed to determine Xa/Xi status for each cell.
#' @param ratioXa (numeric): Threshold for Xa/Xi ratio to call H1 or H2 as Xa/Xi.
#'
#' @examples
#' 1 + 1
#'
#' @return (list): cell_XaXi (status of Xa/Xi per cell), gene_XaXi (Xa/Xi status per gene; on all cells)
#'
#' @importFrom dplyr %>%
#' @export
determine_XaXi_signal <- function(seurat, path_counts, min_counts_cells = 2, ratioXa = 0.1){
    # Input validation ----
    checkmate::assertClass(seurat, "Seurat")
    checkmate::assertFile(path_counts, access = 'r')
    checkmate::assertNumber(min_counts_cells, null.ok = FALSE)
    checkmate::assertNumber(ratioXa, null.ok = FALSE)

    # Import H1/H2/UA counts ----
    futile.logger::flog.info("determine_XiXa_signal: Importing H1/H2/UA counts.")

    counts <- data.table::fread(file = path_counts, sep = "\t")

    # Subset on cells.
    counts <- counts[counts$cell %in% colnames(seurat), ]

    # Subset on chromosome X genes.
    genes_X <- seurat[["RNA"]][[]] %>% dplyr::filter(gene_chr == "chrX") %>% dplyr::distinct(gene_name, gene_id)

    # Add gene-name to counts.
    counts <- counts %>%
        dplyr::inner_join(
            y = genes_X,
            by = c('gene' = 'gene_id')
        )

    # Determine cell-wise Xa/Xi status ----

    futile.logger::flog.info("determine_XiXa_signal: Determining Xa/Xi per cell (n = %s); filtering cells with <%s H1 or H2 counts", ncol(seurat), min_counts_cells)

    # Use summed H1/H2 counts over all chrX genes (except Xist/Tsix)
    cell_XaXi <- counts[! counts$gene_name %in% c('Xist', 'Tsix'), ] %>%
        dplyr::group_by(cell) %>%
        dplyr::summarise(
            total_H1 = sum(H1, na.rm = T),
            total_H2 = sum(H2, na.rm = T),
            total_counts_cell = total_H1 + total_H2,
            ratio = sum(H1, na.rm = T) / total_counts_cell,
        ) %>%

        # Determine if H1 or H2 has the most amount of reads (ratio > pre-defined ratio)
        dplyr::mutate(
            Xa_cell = dplyr::case_when(
                ratio <= ratioXa ~ "H2",
                ratio >= 1-ratioXa ~ "H1",
                .default = "UA"
            ),
            Xi_cell = dplyr::case_when(
                ratio <= ratioXa ~ "H1",
                ratio >= 1-ratioXa ~ "H2",
                .default = "UA"
            ),
        ) %>%
        dplyr::filter(total_counts_cell >= min_counts_cells) %>%
        dplyr::ungroup()

    # Determine gene-wise Xa/Xi status ----

    futile.logger::flog.info("determine_XiXa_signal: Determining Xa/Xi ratio of genes (over all remaining cells: %s)", nrow(cell_XaXi))

    gene_XaXi <- counts %>%
        dplyr::left_join(cell_XaXi, by = 'cell') %>%
        dplyr::filter(Xi_cell != "UA" & Xa_cell != "UA") %>%
        dplyr::group_by(gene_name) %>%
        dplyr::summarise(
            Xa_gene = sum(dplyr::case_when(
                Xa_cell == "H1" ~ H1,
                Xa_cell == "H2" ~ H2,
                .default = NA
            ), na.rm = T),
            Xi_gene = sum(dplyr::case_when(
                Xi_cell == "H1" ~ H1,
                Xi_cell == "H2" ~ H2,
                .default = NA
            ), na.rm = T)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            total = Xa_gene + Xi_gene,
            ratio = Xi_gene / total,
            ratio = dplyr::if_else(is.nan(ratio), 0, ratio)
        )

    # Return list ----
    return(list('cell_XaXi' = cell_XaXi, 'gene_XaXi' = gene_XaXi))

}
