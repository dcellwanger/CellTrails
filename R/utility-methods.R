#' Pretty string from array
#'
#' Generates short representation of long
#' character vectors
#' @param x A character vector
#' @return String
#' @details Returns pretty print of character vector
#' @keywords internal
#' @author Daniel C. Ellwanger
.prettyString <- function(x, mmax=4) {
  l <- length(x)
  mmax <- max(4, mmax)
  op <- options("useFancyQuotes")
  options(useFancyQuotes=FALSE)
  x <- dQuote(x)
  options(useFancyQuotes = op)
  if(l == 0) {
    return("none")
  } else if (l < mmax) {
    return(paste0(paste0(x, collapse = " "), " (", l, ")"))
  } else {
    return(paste0(x[1], " ", x[2], " ... ", x[l], " (", l, ")"))
  }
}

#' Checks if feature name exists
#'
#' @param x A \code{\link{CellTrailsSet}} object
#' @param feature_name Name of feature
#' @import Biobase
#' @keywords internal
#' @author Daniel C. Ellwanger
.checkFeatureNameExists <- function(x, feature_name) {
  if(!feature_name %in% featureNames(x)) {
    stop("Feature '", feature_name, "' is not contained in data set.")
  }
}

#' Compute differential expression
#'
#' Computes P-value and fold-change between two expression vectors.
#' @param x Feature expression in condition x
#' @param y Feature expression in condition y
#' @param lod Gene limit of detection
#' @return A list containing the following components:
#' @return \item{\code{p.value}}{P-value}
#' @return \item{\code{fold}}{Fold-change}
#' @details For censored data a Peto-Peto test is performed, for non-censored data a Wilcoxon rank sum test.
#' If limit of detection is provided, expectation maximization is used to compute the
#' fold-change.
#' @importFrom EnvStats enormCensored twoSampleLinearRankTestCensored
#' @keywords internal
#' @author Daniel C. Ellwanger
.diffExpr <- function(x, y, lod=NULL, alternative="two.sided") {
  f.mean <- function(x, lod) {
    x.censored <- x == 0

    res <- list()
    if(sum(x.censored) == 0 | length(unique(x)) < 3) {
      res$mean <- mean(x)
      res$sd <- sd(x)
      res$CI <- f.CI(length(x), mean(x), sd(x))
    } else { #expectation maximization
      x[x.censored] <- lod
      x.params <- enormCensored(x=x, censored=x.censored,
                                censoring.side="left", seed=110101,
                                method="mle", ci=TRUE)
      res$mean <- x.params$parameters[1]
      res$sd <- x.params$parameters[2]
      res$CI <- x.params$interval$limits
    }
    res
  }

  f.CI <- function(n, mean, sd) {
    error <- qnorm(0.975) * sd / sqrt(n)
    c(mean - error, mean + error)
  }

  x.censored <- x == 0
  y.censored <- y == 0

  n.x <- sum(!x.censored)
  n.y <- sum(!y.censored)

  pval <- 1
  if(length(x) == 0 | length(y) == 0) {
    pval <- 1
  } else if(n.x + n.y == 0) { #all values are censored
    pval <- 1
  } else {
    if((n.x + n.y) == (length(x) + length(y))) { #no censored values
      pval <- suppressWarnings(wilcox.test(x, y, alternative=alternative, exact=FALSE)$p.value)
    } else {
      pval <- twoSampleLinearRankTestCensored(x=x, x.censored=x.censored,
                                              y=y, y.censored=y.censored,
                                              test="peto-peto",
                                              censoring.side="left",
                                              variance="permutation", alternative=alternative)$p.value
    }
  }
  res <- list()
  res$p.value <- ifelse(is.na(pval), 1, pval)
  res$fold <- mean(x) - mean(y)

  #expectation maximization
  if(!is.null(lod)) {
    res$fold <- NA
    fold <- f.mean(x, lod = lod)$mean - f.mean(y, lod = lod)$mean
    names(fold) <- NULL
    res$fold <- fold
  }
  res
}

#' Radial basis function
#'
#' Computes radial basis function
#' @param x A numeric value, vector or matrix
#' @param sigma Radius of the kernel
#' @return Transformed values
#' @details Also known as Heat kernel or Gaussian kernel
#' @keywords internal
#' @author Daniel C. Ellwanger
.rbf <- function(d, radius) {
  exp((-1) * (d * d) / (2 * radius * radius))
}

#' Reverse inverse hyperbolic sine
#'
#' Computes inverse hyperbolic sine
#' @param x A numeric value, vector or matrix
#' @return Transformed values
#' @details In contrast to sqrt and log, ihs is defined for negative numbers.
#' @keywords internal
#' @author Daniel C. Ellwanger
.ihs <- function(x) {
  log(x + sqrt(x ^ 2 + 1))
}

#' Color palette
#'
#' Generates equally spaced hues around the color wheel
#' @param n Number of colors to be generated
#' @return Colors codes
#' @importFrom grDevices hcl
#' @keywords internal
#' @author Daniel C. Ellwanger
.color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Colors for vector
#'
#' Generates colors for vector
#' @param x A numeric vector
#' @return Color codes
#' @importFrom viridis viridis
#' @keywords internal
#' @author Daniel C. Ellwanger
.color_ramp <- function(x, range=3, colPal=NULL, min.val=1e-10, min.val.col=NA, breaks=100) {
  if(is.null(colPal)) {
    colPal <- viridis
  }

  x[x < min.val] <- NA

  lower <- which(x < mean(x, na.rm=TRUE))
  upper <- which(x >= mean(x, na.rm=TRUE))

  fence <- mean(x, na.rm=TRUE) - range * sd(x[lower])
  x[lower][x[lower] < fence] <- fence

  fence <- mean(x, na.rm=TRUE) + range * sd(x[upper])
  x[upper][x[upper] > fence] <- fence

  m1 <- as.numeric(cut(x[lower],
                       breaks=floor(breaks / 2),
                       include.lowest=TRUE)) + 1
  m2 <- as.numeric(cut(x[upper],
                       breaks=floor(breaks / 2),
                       include.lowest=TRUE)) + max(m1, na.rm=TRUE)

  x[lower] <- m1
  x[upper] <- m2
  x[is.na(x)] <- 1
  c(min.val.col, colPal(max(x) - 1))[x]
}

#' Rescale vector
#'
#' Rescales vector to [ymin, ymax]
#' @param x X-values
#' @param ymin New min
#' @param ymax New max
#' @return Rescaled vector
#' @keywords internal
#' @author Daniel C. Ellwanger
.rescale <- function(x, ymin, ymax) {
  if(length(x) == 1) {
    mean(c(ymin, ymax))
  } else {
    (ymax - ymin) / (max(x, na.rm=TRUE) -  min(x, na.rm=TRUE)) * (x - min(x, na.rm=TRUE)) + ymin
  }
}

#' Denoises expression
#'
#' Blocks factors in expression matrix
#' @param x A numeric expression vector
#' @param design A numeric matrix describing the factors that should be blocked
#' @return A numeric denoised expression vector
#' @keywords internal
#' @author Daniel C. Ellwanger
.denoiseExpression <- function(x, design) {
  isZero <- x == 0 #perform linear fit w/o zeros
  y <- x[!isZero]
  design <- design[!isZero, ]
  fit <- lm(y ~ design + 0)
  x[!isZero] <- residuals(fit) + coefficients(fit)[1]

  # Replacing non-zero values with the scaled residuals
  xvar <- var(y)
  rvar <- var(x)
  x * sqrt(xvar/rvar)
}

#' Simulation of RNA-Seq expression data
#'
#' Simple simulation of RNA-Seq expression data estimating counts based
#' on the negative binomial distribution
#' @param n_features Number of genes
#' @param n_samples Number of samples
#' @param prefix_sample Prefix of sample name
#' @param seed Seed of pseudo-random number generator; used for random generation of data
#' @return A numeric matrix with genes in rows and samples in columns
#' @details RNA-Seq counts are generated using the Negative Binomial Distribution.
#' Distribution parameters for each feature are sampled from a Gamma distribution. The
#' resulting expression matrix is log2-scaled.
#' @seealso \code{\link[stats]{NegBinomial}} and \code{\link[stats]{GammaDist}}
#' @examples
#' # Matrix with 100 genes and 50 cells
#' dat <- simulate_exprs(n_features=100, n_samples=50)
#' @docType methods
#' @export
#' @author Daniel C. Ellwanger
simulate_exprs <- function(n_features, n_samples, prefix_sample="", seed=1101) {
  set.seed(seed)
  s_means <- 2^rgamma(n_features, shape=3.5, rate=1.2)
  s_vars <- 2^rgamma(n_features, shape=12, rate=1)
  s_counts <- suppressWarnings(
    rnbinom(n_features * n_samples, mu=s_means, size=s_means^2/(s_vars-s_means)))
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
#' @return An \code{\link[Biobase]{ExpressionSet}}
#' @details RNA-Seq counts are generated using the Negative Binomial Distribution.
#' Distribution parameters for each feature are sampled from a Gamma distribution. The
#' resulting expression matrix is then log2-scaled. The expression data consist of
#' 25 features and 100 samples
#' @seealso \code{\link[Biobase]{ExpressionSet}} and \code{\link[CellTrails]{simulate_exprs}}
#' @examples
#' # Generate example data
#' exDat
#' @docType methods
#' @export
#' @author Daniel C. Ellwanger
exDat <- function() {
  expr <- lapply(1:10, function(i)
    simulate_exprs(seed = 1101 + i,
                   n_features=25,
                   n_samples=10,
                   prefix_sample=paste0("C", i, "_")))
  expr <- do.call(cbind, expr)
  meta <- rep(c("2.Mid", "2.Mid", "3.Late", "1.Early", NA,
                "3.Late", "3.Late", NA, NA, "2.Mid"), each = 10)
  ExpressionSet(expr,
                phenoData = new("AnnotatedDataFrame",
                                data.frame(age=meta, row.names=colnames(expr))))
}

#' Enrichment test
#'
#' Statistical enrichment analysis using either a Hypergeometric
#' or Fisher's test
#' @param sample_true Number of hits in sample
#' @param sample_size Size of sample
#' @param pop_true Number of hits in population
#' @param pop_size Size of population
#' @param method Statistical method that should be used
#' @return A list containing the following components:
#' @return \item{\code{p.value}}{P-value of the test}
#' @return \item{\code{odds.ratio}}{Odds ratio}
#' @return \item{\code{conf.int}}{Confidence interval for the odds ratio (only shown with method="fisher")}
#' @return \item{\code{method}}{Used statistical test}
#' @details Hypergeometric or one-tailed Fisher exact test is useful for enrichment analyses.
#' For example, one needs to estimate which features are enriched among
#' a set of instances sampled from a population.
#' @seealso \code{\link[stats]{Hypergeometric}} and \code{\link[stats]{fisher.test}}
#' @examples
#' # Population has 13 of total 52 instances positive for a given feature
#' # Sample has 1 of total 5 instances positive for a given feature
#' # Test for significance of enrichment in sample
#' enrichment.test(sample_true=1, sample_size=5, pop_true=13, pop_size=52, method="fisher")
#' @docType methods
#' @export
#' @author Daniel C. Ellwanger
enrichment.test <- function(sample_true, sample_size,
                            pop_true, pop_size, method=c("fisher", "hyper")) {

  method <- toupper(method[1])
  res <- list()
  if(method == "HYPER") {
    p <- phyper(sample_true-1,
                sample_size,
                pop_size-sample_size,
                pop_true,
                lower.tail=FALSE)
    oA <- prod(sample_true, pop_size-pop_true+sample_true)
    oB <- prod(sample_size-sample_true, pop_true-sample_true)
    or <- oA / oB

    res$p.value <- p
    res$odds.ratio <- p
    res$method <- "Hypergeometric test for enrichment"
  } else if (method == "FISHER") {
    m <- matrix(c(sample_true,
                  pop_true-sample_true,
                  sample_size-sample_true,
                  pop_size-pop_true-sample_size+sample_true), 2, 2)
      ft <- fisher.test(m, alternative='greater')

      res$p.value <- ft$p.value
      res$odds.ratio <- ft$estimate
      res$conf.int <- ft$conf.int
      res$method <- "Fisher's exact test for enrichment"
  } else {
    stop('Unknown method select. Please, choose from c("hyper", "fisher").')
  }
  res
}
