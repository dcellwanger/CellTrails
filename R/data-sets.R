#' Example single-cell expression data
#'
#' This dataset contains simulated transcript expression profiles of
#' 25 genes in expressed 100 cells. Simulation was performed using
#' using the Negative Binomial Distribution. Distribution parameters
#' for each feature were sampled from a Gamma distribution.
#' The resulting expression matrix is log2-scaled and was stored in
#' in an object of class 'SingleCellExperiment' (assay \code{logcounts}).
#' The sample metainformation contains the underlying (discrete) simulated
#' \code{age} of the cells.
#' @format An object of class \code{SingleCellExperiment}
#' @usage data(exSCE)
"exSCE"
