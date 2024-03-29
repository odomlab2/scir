% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/determine_XaXi_signal.R
\name{determine_XaXi_signal}
\alias{determine_XaXi_signal}
\title{Determine Xa/Xi status from haplotyping data.}
\usage{
determine_XaXi_signal(seurat, path_counts, min_counts_cells = 2, ratioXa = 0.1)
}
\arguments{
\item{seurat}{(Seurat): Seurat object.}

\item{path_counts}{(str): Path to the H1/H2/UA counts table.}

\item{min_counts_cells}{(numeric): Min. number of H1 or H2 counts needed to determine Xa/Xi status for each cell.}

\item{ratioXa}{(numeric): Threshold for Xa/Xi ratio to call H1 or H2 as Xa/Xi.}
}
\value{
(list): cell_XaXi (status of Xa/Xi per cell), gene_XaXi (Xa/Xi status per gene; on all cells)
}
\description{
Determines cell-wise Xa/Xi status and Xa/Xi status of chromosome X genes.
}
\examples{
1 + 1

}
