#' DEF: Imputation of drop-outs
#'
#' For details see \code{\link{imputeDropouts}}
#' @importFrom igraph neighbors
#' @import Biobase
#' @keywords internal
#' @author Daniel C. Ellwanger
.imputeDropouts_def <- function(x, feature_name) {
  X <- x[feature_name, x@useSample]
  X[X == 0] <- NA
  g <- x@trajectory$traj
  rep.total <- sum(is.na(X))

  f.impute <- function(g, y) {
    nds <- which(is.na(y))
    for(nd in nds) {
      nghbr <- as.vector(neighbors(g, V(g)[nd]))
      y[nd] <- mean(y[nghbr])
    }
    y[is.na(y)] <- 0
    y
  }

  message("Performing imputation for ", length(feature_name),  " features ...")
  X.imp <- X
  rep.f <- round(nrow(X) * .1)
  for(i in seq_len(nrow(X))) {
    X.imp[i, ] <- f.impute(g = g, y = X[i, ])
    if(i %% rep.f == 0) {
      rep.ana <- (rep.total - sum(is.na(X.imp)))
      rep.douts <- round((rep.ana - sum(X.imp == 0, na.rm = TRUE)) / rep.ana * 100, 1)
      message(round(rep.ana / rep.total * 100), "% non-detects analyzed (",
              rep.douts, "% drop-outs) ...")
    }
  }
  message("100% non-detects analyzed (",
          round(sum(X.imp[is.na(X)] != 0) / rep.total * 100, 1),
          "% drop-outs).")
  X.imp[is.na(X.imp)] <- 0
  exprs(x)[feature_name, x@useSample] <- X.imp
  x
}
