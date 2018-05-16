#' DEF: Filter feaures by Detection Level
#'
#' For details see \code{filterFeaturesByPOD}
#' @param y An expression vector
#' @param threshold numeric cutoff value
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.filterTrajFeaturesByDL_def <- function(y, threshold) {
  pod <- apply(y, 1L, function(i){sum(i > 0)})

  if(threshold >= 1) {
    threshold <- threshold / max(pod)
    pod <- pod / max(pod)
  }
  f <- pod > threshold

  #Diagnostic plot
  x.steps <- sort(unique(pod))
  y.ecdf <- ecdf(pod)(x.steps)
  dat <- data.frame(X = x.steps, Y = y.ecdf,
                    COLOR = c("Rejected",
                              "Accepted")[(x.steps > threshold) + 1])
  gp <- ggplot(dat, aes_string(x = "X", y = "Y")) +
        geom_step(alpha = .5) +
        geom_point(aes_string(color = "COLOR")) +
        labs(colour = "Filter") +
        xlab("Detection level") + ylab("Features (cumulative fraction)") +
        theme(axis.line = element_line(colour = "black"))
  print(gp)

  #Update attribute
  #trajFeatureNames(x) <- names(which(f))
  #x
  names(which(f))
}

#' DEF: Filter features by coefficient of variation
#'
#' For details see \code{filterFeaturesByCOV}
#' @param y An expression vector
#' @param threshold numeric cutoff value
#' @param design Model matrix
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.filterTrajFeaturesByCOV_def <- function(y, threshold, design=NULL) {
  if(!is.null(design)) {
    message("Blocking nuisance factors ...")
    y <- t(apply(y, 1L, .denoiseExpression, design))
    y[y < 0] <- 0
  }

  fcov <- apply(y, 1L, function(i){sd(i) / mean(i)})
  f <- fcov > threshold

  #Diagnostic plot
  x.steps <- sort(unique(fcov))
  y.ecdf <- ecdf(fcov)(x.steps)

  dat <- data.frame(X = x.steps, Y = y.ecdf,
                    COLOR = c("Rejected",
                              "Accepted")[(x.steps > threshold) + 1])
  gp <- ggplot(dat, aes_string(x = "X", y = "Y")) +
    geom_step(alpha = .5) +
    geom_point(aes_string(color = "COLOR")) +
    labs(colour = "Filter") +
    xlab("Coefficient of variation") + ylab("Features (cumulative fraction)") +
    theme(axis.line = element_line(colour = "black"))
  print(gp)

  #Update attributes
  #trajectoryFeatures(x) <- names(which(f))
  #x
  names(which(f))
}

#' DEF: Filter features by index of dispersion / fano factor
#'
#' For details see \code{filterFeaturesByFF}
#' @param y An expression vector
#' @param z A cutoff z-score
#' @param min_expr Minimum expression level
#' @param design Model matrix
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.filterTrajFeaturesByFF_def <- function(y, z, min_expr=0, design=NULL) {
  if(!is.null(design)) {
    message("Blocking nuisance factors ...")
    y <- t(apply(y, 1L, .denoiseExpression, design))
    y[y < 0] <- 0
  }

  stat.mean <- apply(y, 1L, mean)
  stat.disp <- apply(y, 1L, function(x) var(x) / mean(x))
  stat.bin <- cut(stat.mean, breaks = 20)
  stat.z <- stat.disp
  for(i in levels(stat.bin)) {
    stat.z[stat.bin == i] <- scale(stat.z[stat.bin == i])
  }
  stat.z[is.na(stat.z)] <- 0

  f1 <- stat.z > z
  f2 <- stat.mean > min_expr
  f3 <- f1 & f2

  f <- rep(NA, length(stat.z))
  f[f3] <- "Accepted"
  f[!f3] <- "Rejected"

  #Diagnostic plot
  dat <- data.frame(mean = stat.mean, dispersion = stat.disp, f = f)
  gp <- ggplot() +
        geom_point(data = dat, aes_string(x = "mean", y = "dispersion",
                                          color = "f")) +
        theme(axis.line = element_line(colour = "black")) +
        labs(colour = "Filter") +
        xlab("Mean") + ylab("Fano factor")
  print(gp)

  #Update attribute
  #trajectoryFeatures(x) <- names(which(f3))
  #x
  names(which(f3))
}

# #' DEF: Filter samples by reference gene
# #'
# #' Performs filtering of cells using reference gene information.
# #' @param x An \code{SingleCellExperiment} object
# #' @param refgene Symbol name(s) of reference gene(s)
# #' @param fence Fence to find outlier
# #' @return An \code{SingleCellExperiment} object
# #' @seealso \code{SingleCellExperiment}, \code{quantile}, \code{IQR}
# #' @details Identifies outliers (dead cells and multiplets) based on robust
# #' statistics on expression intensities of reference genes. It filters cells
# #' having a reference gene expression of
# #' \eqn{x > Q_{.25} - f * IQR} or \eqn{x < Q_{.75} + f * IQR},
# #' where \eqn{Q} then quantile, \eqn{IQR} the interquartile function on the
# #' reference gene distribution and \eqn{f} is the fence parameter.
# #' If multiple
# #' reference genes are provided their geometric mean is used. As a rule of
# #' thumb, \code{fence} should be
# #' between 3 (extreme values) and 1.5 (outliers). Further, cells
# #' expressing all reference genes are filtered.
# #' @import Biobase
# #' @keywords internal
# #' @author Daniel C. Ellwanger
#.filterSamplesByReference_def <- function(x, refgene,
#                                          fence=1.5, design=NULL) {
#   .featureNameExists(x, refgene)
#   edat <- .exprs(x[refgene, ])
#
#   if(!is.null(design)) {
#     edat <- apply(edat, 1L, .denoiseExpression, design)
#   }
#
#   # Filter by robust statistics of reference gene expression
#   ref <- apply(edat, 2L, function(x){mean(x, na.rm=TRUE)})
#   lo <- quantile(ref, .25, na.rm=TRUE) - fence * IQR(ref, na.rm=TRUE)
#   hi <- quantile(ref, .75, na.rm=TRUE) + fence * IQR(ref, na.rm=TRUE)
#   f <- ref > lo & ref < hi
#   f[is.na(f)] <- FALSE
#   #x <- x[, f]
#
#   # Filter cells not expressing the reference genes
#   #edat <- .exprs(x[refgene, ])
#   f <- apply(edat, 2L, function(x){all(x[refgene] > 0)}) & f
#   x[, f]
# }
