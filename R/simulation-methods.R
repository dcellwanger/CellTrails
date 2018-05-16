#' Simulation of RNA-Seq expression data
#'
#' Simple simulation of RNA-Seq expression data estimating counts based
#' on the negative binomial distribution
#' @param n_features Number of genes
#' @param n_samples Number of samples
#' @param prefix_sample Prefix of sample name
#' @param seed Seed of pseudo-random number generator; used for random
#' generation of data
#' @return A numeric matrix with genes in rows and samples in columns
#' @details RNA-Seq counts are generated using the Negative Binomial
#' Distribution. Distribution parameters for each feature are sampled
#' from a Gamma distribution. The
#' resulting expression matrix is log2-scaled.
#' @seealso \code{NegBinomial} and \code{GammaDist}
#' @examples
#' # Matrix with 100 genes and 50 cells
#' dat <- simulate_exprs(n_features=100, n_samples=50)
#' @docType methods
#' @export
#' @author Daniel C. Ellwanger
simulate_exprs <- function(n_features, n_samples,
                           prefix_sample="", seed=1101) {
  set.seed(seed)
  s_means <- 2^rgamma(n_features, shape=3.5, rate=1.2)
  s_vars <- 2^rgamma(n_features, shape=12, rate=1)
  s_counts <- suppressWarnings(
    rnbinom(n_features * n_samples, mu=s_means,
            size=s_means^2/(s_vars-s_means)))
  s_counts[is.na(s_counts)] <- 0
  s_counts <- log2(s_counts + 1)
  s_counts <- matrix(s_counts, nrow=n_features, ncol=n_samples)
  rownames(s_counts) <- paste0("feature_", seq_len(n_features))
  colnames(s_counts) <- paste0(prefix_sample, "sample_", seq_len(n_samples))
  s_counts
}

#' Simulation of example data for trajectory reconstruction
#'
#' This method serves to generate the example data used to
#' demonstrate the usage of individual functions of this package
# #' @param n_features Number of features
# #' @param n_samples Number of samples
# #' @param n_cond Number of conditions
#' @return A \code{SingleCellExperiment} object
#' @details RNA-Seq counts are generated using the Negative Binomial
#' Distribution. Distribution parameters for each feature are sampled from a
#' Gamma distribution. The resulting expression matrix is then log2-scaled.
#' The assay data consist of 25 features and 100 samples.
#' @seealso \code{SingleCellExperiment} \code{simulate_exprs}
#' @examples
#' # Generate example data
#' exDat()
#' @docType methods
#' @import SingleCellExperiment SingleCellExperiment
#' @export
#' @author Daniel C. Ellwanger
exDat <- function() { #n_features=25, n_samples=10, n_cond=10
  expr <- lapply(seq_len(10), function(i)
    simulate_exprs(seed = 1101 + i,
                   n_features=25,
                   n_samples=10,
                   prefix_sample=paste0("C", i, "_")))
  expr <- do.call(cbind, expr)
  meta <- rep(c("2.Mid", "2.Mid", "3.Late", "1.Early", NA,
                "3.Late", "3.Late", NA, NA, "2.Mid"), each = 10)
  SingleCellExperiment(assays = list(logcounts=expr),
                       colData=data.frame(age=meta))
  #ExpressionSet(expr,
  #              phenoData = new("AnnotatedDataFrame",
  #                              data.frame(age=meta,
  #                                         row.names=colnames(expr))))
}
