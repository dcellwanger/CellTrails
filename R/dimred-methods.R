#' DEF: Spectral embedding of samples
#'
#' For details see \code\link{embedSamples}
#' @param nbins Cubic B-spline discretization is used to compute fuzzy mutual
#' information between pairs of samples; \code{nbins} defines the number of intervals
#' used for discretization of expression data. (default: 10)
#' @param sigma The mutual information matrix is transformed to an adjacency matrix
#' using a heat kernel; \code{sigma} defines the radius of heat kernel (in quantiles;
#' default: \code{sigma} = 0.75 which is the third quantile of the mutual information matrix).
#' @importFrom splines bs
#' @keywords internal
#' @author Daniel C. Ellwanger
.embedSamples_def <- function(x, nbins=10, sigma=.75, design=NULL) {
  f.W <- function(M, nbins) {
    m <- nrow(M)
    n <- ncol(M)
    W <- array(NA, dim = c(nbins, n, m)) #Weight matrix, <bin, genes, cells>

    # Cubic B-splines
    for(i in seq_len(m)) {
      x <- M[i, ]
      k <- seq(from=min(x, na.rm=TRUE), to=max(x, na.rm=TRUE), length.out=nbins - 2)
      W[1:nbins, , i] <- t(bs(x, knots=k[-c(nbins - 2)],
                              Boundary.knots=range(k))) #cubic, degree = 3
    }

    # Avoid NAs
    W[is.na(W)] <- 0
    W
  }

  f.MI <- function(W.x, W.y, pmar.x, pmar.y, n) {
    pjoint <- tcrossprod(W.x, W.y) / n #joint probability
    denom <- tcrossprod(pmar.x, pmar.y) #denominator
    weight <- tcrossprod(1 - pmar.x, 1 - pmar.y) #weights
    sum(weight * pjoint * log2(pjoint / denom), na.rm=TRUE)
  }

  f.I <- function(W) {
    m <- dim(W)[3]
    n <- dim(W)[2]
    pmar <- apply(W, c(1L, 3L), function(x) { # Marginal probabilities
      sum(x, na.rm = TRUE) / (length(x) - sum(is.na(x)))})

    I <- matrix(ncol=m, nrow=m)
    for(i in seq_len(m)) {
      for(j in i:m) {
        I[i, j] <- I[j, i] <- f.MI(W[,,i], W[,,j], pmar[,i], pmar[,j], n)
      }
    }
    I
  }

  f.mi2dist <- function(I) {
    denom <- 1/sqrt(diag(I))
    I <- sweep(I, MARGIN = 2L, STATS = denom, FUN = "*", check.margin = FALSE)
    I <- sweep(I, MARGIN = 1L, STATS = denom, FUN = "*", check.margin = FALSE)
    diag(I) <- 1
    I[I > 1] <- 1 #upper bound
    d <- 1 - I
    sqrt(2 * d)
  }

  #Pre-flight check
  ze <- apply(x[x@useFeature, ], 1L, function(x){sum(x > 0)})
  if(sum(ze == 0) > 0) {
    warning(ze, " samples do not encode trajectory information (i.e. have zero entropy).
            You may want to select a proper set of trajectory features (see '?trajectoryFeatures').")
  }
  ze <- ze > 0

  #Expression matrix
  #M <- t(x[x@useFeature, ]) #No scaling required
  M <- x[x@useFeature, ] #filter by trajectory features
  M <- M[ze, ] #filter by informative features
  if(nrow(M) < 2) {
    stop("Cannot compute embedding, because less than two features were selected. Please, increase
         the number of trajectory features.")
  }

  if(!is.null(design)) { #block uninteresting factors
    message("Blocking nuisance factors ...")
    M <- apply(M, 1L, .denoiseExpression, design)
  } else {
    M <- t(M)
  }

  message("Computing adjacency matrix ...")
  W <- f.W(M=M, nbins=nbins) #Compute weight matrix
  I <- f.I(W) #Compute mutual information
  Dmat <- f.mi2dist(I) #Convert to distance matrix

  #Adjacency matrix
  radius <- ifelse(is.null(sigma), quantile(Dmat, .75), sigma)
  adjmat <- .rbf(Dmat, radius=radius) #is similarity matrix

  #Result
  res <- list()
  res$I <- I
  res$Dmat <- Dmat
  res$adjmat <- adjmat

  #Spectral embedding
  message("Computing spectral embedding ...")
  L <- adjmat
  diag(L) <- diag(L) + apply(L, 1L, sum) / ncol(L) #signless graph Laplacian
  edecomp <- base::eigen(L, symmetric=TRUE)
  res$values <- edecomp$values
  res$values <- .ihs(res$values) #sqrt(edecomp$values)
  res$X <- sweep(edecomp$vectors, MARGIN=2L, STATS=res$values, FUN="*")

  #Set result to obj
  CellTrails::latentSpace(x) <- res$X
  CellTrails::eigenvalues(x) <- res$values

  x
}

#' DEF: Determine number of informative latent dimensions
#'
#' For details see \code\link{findSpectrum}
#' @importFrom utils head
#' @keywords internal
#' @author Daniel C. Ellwanger
.findSpectrum_def <- function(x, frac=100) {
  D <- CellTrails::eigenvalues(x)

  cs <- cumsum(diff(D))
  f <- ifelse(frac <= 1, length(D) * frac, frac)
  h <- head(cs, f)
  fit <- .linear_fit(x.in=seq_along(h), y.in=h)

  new("CellTrailsSpectrum",
      frac=frac,
      n=min(which(diff(which(h > fit$y)) > 1)) + 1,
      cs=cs,
      fit=fit)
}

#' DEF: Truncate eigenbasis
#'
#' For details see \code\link{reduceDimensions}
#' @keywords internal
#' @author Daniel C. Ellwanger
.reduceDimensions_def <- function(x, s) {
  CellTrails::latentSpace(x) <- CellTrails::latentSpace(x)[, seq_len(s@n)]
  CellTrails::eigenvalues(x) <- CellTrails::eigenvalues(x)[seq_len(s@n)]
  x
}

#' t-Distributed Stochastic Neighbor Embedding
#'
#' Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding
#' @param x An \code{\link{CellTrailsSet}} object
#' @param dims Output dimensionality
#' @param perplexity Perplexity parameter (default: 30)
#' @param theta Speed/accuracy trade-off (increase for less accuracy),
#' set to 0.0 for exact tSNE (default: .5)
#' @param max_iter Number of iterations (default: 1000)
#' @param seed Starting value for pseudorandom number generator. Setting
#' seed makes result reproducible (default: 1101)
#' @return A list with the following components:
#' \describe{
#'   \item{\code{Y}}{Matrix containing the new representations for the objects}
#'   \item{\code{perplexity}}{See above}
#'   \item{\code{seed}}{See above}
#' }
#' @importFrom Rtsne Rtsne
#' @keywords internal
#' @author Daniel C. Ellwanger
.bhtsne <- function(x, dims=2, perplexity=30, theta=.5, max_iter=1000, seed=1101){
  if(!is.null(seed)) {
    set.seed(seed)
  }
  result <- list()
  result$Y <- Rtsne(CellTrails::latentSpace(x), dims=dims, pca=FALSE,
                    perplexity=perplexity, theta=theta, max_iter=max_iter)$Y
  result$perplexity <- perplexity
  result$seed <- seed
  result
}

#' DEF: PCA
#'
#' For details see \code\link{pca}
#' @keywords internal
#' @author Daniel C. Ellwanger
.pca <- function(x, do_scaling=TRUE, design=NULL) {
  X <- t(x[x@useFeature, ])

  # Remove features with 0 variance
  f <- apply(X, 2L, var) > 0
  X <- X[, f]

  if(!is.null(design)) {
    message("Blocking nuisance factors ...")
    X <- apply(X, 2L, .denoiseExpression, design)
  }

  message("Performing PCA ...")

  pca_result <- stats::prcomp(X, scale.=do_scaling)

  # Scale matrix
  #if(do_scaling) {
  #  X <- scale(X)
  #}
  # Compute covariance matrix
  #covmat <- cov(X)
  # Eigenvalue decomposition
  #edecomp <- base::eigen(covmat)

  # Result
  res <- list()
  #res$princomp <- X %*% edecomp$vectors
  #colnames(res$princomp) <- paste("PC", 1:ncol(X))
  #res$variance <- edecomp$values #/ sum(edecomp$values)
  #res$sdev <- sqrt(res$variance)
  #res$loadings <- sweep(edecomp$vectors, MARGIN = 2,
  #                      STATS = sqrt(edecomp$values), FUN = "*")

  res$princomp <- pca_result$x
  res$variance <- pca_result$sdev^2
  res$sdev <- pca_result$sdev
  res$loadings <- pca_result$rotation
  res
}
