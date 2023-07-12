#' @title import_samples_monocle
#' @description Import one or multiple sci-RNA-seqv3 samples into a monocle3 object.
#'
#' @param samples (character): Path(s) to STARSolo output folders (filtered).
#' @param gtf (character): Path to GTF file used during alignment and for annotation purposes.
#'
#' @examples
#' 1 + 1
#'
#' @return (cell_data_set): Samples and annotations in monocle3 object.
#'
#' @importFrom dplyr %>%
#' @export
import_samples_monocle <- function(samples, gtf) {

    # Input validation ----
    checkmate::assertFileExists(samples, access = "r")
    checkmate::assertFileExists(gtf, access = "r")

    # Import ----

    monocle3::new_cell_data_set(
        expression_data = "",
        cell_metadata = "",
        gene_metadata = ""
    ) %>%
        dplyr::mutate()

}
