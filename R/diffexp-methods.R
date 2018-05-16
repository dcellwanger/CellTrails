#' DEF: Differential expression between states
#'
#' @importFrom dtw dtw asymmetric
#' @keywords internal
#' @author Daniel C. Ellwanger
.diffExprState_def <- function(x, state1, state2,
                               feature_name,
                               alternative="two.sided", lod=NULL) {
  state1 <- toupper(state1)
  state2 <- toupper(state2)

  exp1 <- x[feature_name, x$STATE == state1]
  exp2 <- x[feature_name, x$STATE == state2]
  .diffExpr(x = exp2, y = exp1, alternative = alternative, lod = lod)
}

#' DEF: Differential expression between trails
#' @param ptime_1 Pseudotime trail 1
#' @param ptime_2 Pseudotime trail 2
#' @param feature_expr Expression vector
#' @param sts Samples vector
#' @param n Number of fitted values
#' @param score Score type
#' @param dynamicFit_k Freedom param for dynamic fit
#' @param df Freedom param for splines fit
#' @keywords internal
#' @author Daniel C. Ellwanger
.contrastExprTrail_def <- function(ptime_1, ptime_2, feature_expr, sts,
                                   trail_names, n=250, score,
                                   dynamicFit_k=5, df=10) {
  f.rmsd <- function(q, r, m) { #root-mean-squared deviation
    sqrt(mean((q$y - r$y)^2))
  }

  f.td <- function(q, r, m) { #total deviation
    sum(abs((q$y - r$y)))
  }

  f.abc <- function(q, r, m) { #area between curves
    integrate(f = function(x) {
      abs(predict(q$mod, x)$y - predict(r$mod, x)$y)
    }, lower = 0, upper = m)$value
  }

  f.cor <- function(q, r, m) { #correlation
    cor(q$y, r$y, method = "pearson")
  }

  f <- switch(score,
              "rmsd" = f.rmsd,
              "td" = f.td,
              "abc" = f.abc,
              "cor" = f.cor)

  #Fetch trail info
  samples_1 <- !is.na(ptime_1)
  ptime_1 <- ptime_1[samples_1]
  expr_1 <- feature_expr[samples_1]
  states_1 <- sts[samples_1]
  samples_2 <- !is.na(ptime_2)
  ptime_2 <- ptime_2[samples_2]
  expr_2 <- feature_expr[samples_2]
  states_2 <- sts[samples_2]

  #Fit dynamics
  if(max(ptime_1) < max(ptime_2)) {
    reference <- .fitDynamic_def(x=ptime_1, y=expr_1,
                                 z=states_1, k=dynamicFit_k) #, x.out = ptime_1
    query <- .fitDynamic_def(x=ptime_2, y=expr_2,
                             z=states_2, k=dynamicFit_k) #, x.out = ptime_1
    reference_name <- trail_names[1]
    query_name <- trail_names[2]
    ptime <- ptime_2
    query$mod <- query$gam
  } else {
    reference <- .fitDynamic_def(x=ptime_2, y=expr_2,
                                 z=states_2, k=dynamicFit_k) #, x.out = ptime_2
    query <- .fitDynamic_def(x=ptime_1, y=expr_1,
                             z=states_1, k=dynamicFit_k) #, x.out = ptime_2
    reference_name <- trail_names[2]
    query_name <- trail_names[1]
    ptime <- ptime_1
    reference$mod <- reference$gam
  }

  if(length(n) == 1) {
    n <- seq(0, max(ptime), length.out=n)
  }

  if(var(reference$y) == 0 | var(query$y) == 0){ #any gene is constant
    result <- list()

    if(var(reference$y) == 0) {
      if(max(ptime_1) < max(ptime_2)) {
        query <- .fitDynamic_def(x=ptime_2, y=expr_2,
                                 z=states_2, k=dynamicFit_k, n.out=length(n))
      } else {
        query <- .fitDynamic_def(x=ptime_1, y=expr_1,
                                 z=states_1, k=dynamicFit_k, n.out=length(n))
      }
      reference$y <- rep(mean(reference$y), length(query$y))
      reference$x <- query$x
      result[[score]] <- f(q=query, r=reference, m=max(ptime))
    } else if(var(query$y) == 0) {
      if(max(ptime_1) < max(ptime_2)) {
        reference <- .fitDynamic_def(x=ptime_1, y=expr_1,
                                     z=states_1, k=dynamicFit_k,
                                     n.out=length(n))
      } else {
        reference <- .fitDynamic_def(x=ptime_2, y=expr_2,
                                     z=states_2, k=dynamicFit_k,
                                     n.out=length(n))
      }
      query$y <- rep(mean(query$y), length(reference$y))
      query$x <- reference$x
      result[[score]] <- f(q=query, r=reference, m=max(ptime))
    } else {
      reference$y <- rep(mean(reference$y), length(n))
      reference$x <- n
      query$y <- rep(mean(query$y), length(n))
      query$x <- n
      result[[score]] <- f(q=query, r=reference, m=max(ptime))
    }
    return(result)
  }

  #Run DTW
  result <- list()

  # First and last elements of the query are anchored
  # at the boundaries of the reference (= global alignment)
  alignment <- dtw((query$y) / max(query$y),
                        (reference$y) / max(reference$y),
                        keep.internals=TRUE,
                        step.pattern=asymmetric,
                        dist.method="Euclidean",
                        open.end=FALSE,
                        open.begin=FALSE)
  result$alignment <- alignment

  df <- min(length(unique(query$y[alignment$index1])) - 1,
            length(unique(reference$y[alignment$index2])) - 1)
  mod <- smooth.spline(query$x[alignment$index1],
                       query$y[alignment$index1], df=df) #, tol = 1e-7
  query_aligned <- predict(mod, n)
  query_aligned$y[query_aligned$y < 0] <- 0
  query_aligned$mod <- mod

  mod <- smooth.spline(query$x[alignment$index1],
                       reference$y[alignment$index2], df=df) #, tol = 1e-7
  reference_aligned <- predict(mod, n)
  reference_aligned$y[reference_aligned$y < 0] <- 0
  reference_aligned$mod <- mod

  #4. Return aligned dynamics and difference score
  result$input <- list()
  result$input[[query_name]] <- query
  result$input[[reference_name]] <- reference
  result[[query_name]] <- query_aligned
  result[[reference_name]] <- reference_aligned

  result[[score]] <- f(q=query_aligned, r=reference_aligned, m=max(ptime))
  result
}

# #' DEF: Finds marker genes
# #'
# #' Computes P-value and fold-change between two expression vectors.
# #' @param x A \code{CellTrailsMaps} object
# #' @param k Subgroups for which markers should be determined, numeric vector
# #' @param B Number of Monte Carlo samples
# #' @param alternative Alternative hypothesis for P-value simulation
# #' (less, greater, two.sided)
# #' @return A list containing the following components:
# #' @return \item{\code{diffexp}}{Differential expression score}
# #' @return \item{\code{spec}}{Specificity score}
# #' @return \item{\code{p.value}}{Monte Carlo test P-value}
# #' @import Biobase
# #' @keywords internal
# #' @author Daniel C. Ellwanger
# .find_marker_def <- function(x, k=NULL, B=1000, alternative="greater") {
#
#   X <- exprs(x)
#   cl <- factor(x@phenoData$STATE)
#
#   #Signal to noise statistic (PMID:11807556)
#   f.snr_statistic <- function(x, y) {
#     #(mean(x) - mean(y)) / (sd(x) + sd(y)) #(sd(x) + sd(y))
#     (mean(x) - mean(c(x,y))) / (sd(x)/sqrt(length(x)))
#   }
#
#   #scaled signal to noise
#   f.snr <- function(x, sel) {
#     #snr <- NA
#     # if(EM) {
#     #   lod <- min(x[x > 0])
#     #   snr <- f.mean(x[sel], lod)$mean - f.mean(x[!sel], lod)$mean
#     # } else {
#        snr <- mean(x[sel]) - mean(x[!sel])
#     # }
#     #snr <- f.snr_statistic(x[sel], x[!sel]) #(mean(x[sel]) -
#     #                       mean(x[!sel])) / (sd(x[sel]) + sd(x[!sel]))
#
#     fac <- 1
#     x <- sort(x, decreasing=FALSE)
#
#     if(is.nan(snr)) {
#       snr <- NA
#     } else if(snr < 0) {
#       l <- x[seq_len(sum(sel))]
#       u <- x[(length(x) - sum(!sel) + 1):length(x)]
#       #if(EM) {
#       #  fac <- f.mean(l, lod)$mean - f.mean(u, lod)$mean
#       #} else {
#       #  fac <- abs(mean(u) - mean(l))
#       #}
#
#       fac <- f.snr_statistic(u, l) #abs(mean(u) - mean(l)) / (sd(u) + sd(l))
#     } else if (snr > 0) {
#       u <- x[seq_len(sum(!sel))]
#       l <- x[(length(x) - sum(sel) + 1):length(x)]
#       #if(EM) {
#       #  fac <- f.mean(u, lod)$mean - f.mean(l, lod)$mean
#       #} else {
#       #  fac <- abs(mean(u) - mean(l))
#       #}
#       fac <- f.snr_statistic(l, u) #abs(mean(u) - mean(l)) / (sd(u) + sd(l))
#     }
#     snr / fac
#
#     #if(snr < 0) {
#     #  snr <- -(snr / (mean(sort(x, decreasing = F)[seq_len(sum(sel))]) -
#     #         mean(sort(x, decreasing = T)[seq_len(sum(!sel))])))
#     #} else if (snr > 0) {
#     #  snr <- snr / (mean(sort(x, decreasing = T)[seq_len(sum(sel))]) -
#     #         mean(sort(x, decreasing = F)[seq_len(sum(!sel))]))
#     #}
#   }
#
#   #result object
#   res <- list()
#   if(is.null(k)) {
#     k <- levels(cl) #sort(unique(cl))
#   }
#
#   #foreach cluster
#   #cl.avgexpr <- aggregate(X, by = list(cl = cl), mean)[, 2:(ncol(X) + 1)]
#   for(i in k) {
#     sel <- (cl == i)
#     snr <- apply(X, 2L, f.snr, sel = sel)
#
#     #pvalue via Monte Carlo test (Hope, 1968) with B replicates
#     p <- matrix(NA, nrow = ncol(X), ncol = B)
#     rownames(p) <- colnames(X)
#     for(b in seq(B)) {
#       set.seed(10010 + b)
#       p[,b] <- apply(X, 2, f.snr, sel = sample(sel))
#     }
#     pval <- rep(NA, ncol(X))
#     if(alternative == "two.sided") {
#       p.upper <- apply(p, 2L, function(x){x >= snr}) #geq
#       p.lower <- apply(p, 2L, function(x){x <= snr}) #leq
#       pval.upper <- rowSums(p.upper) / B
#       pval.lower <- rowSums(p.lower) / B
#       pval[snr >= 0] <- pval.upper[snr >= 0]
#       pval[snr < 0] <- pval.lower[snr < 0]
#     } else if(alternative == "greater") {
#       p.upper <- apply(p, 2L, function(x){x >= snr}) #geq
#       pval <- rowSums(p.upper) / B
#     } else if(alternative == "less") {
#       p.lower <- apply(p, 2L, function(x){x <= snr}) #leq
#       pval <- rowSums(p.lower) / B
#     } else {
#       stop("")
#     }
#     i <- as.character(i)
#     res[[i]] <- list()
#     res[[i]]$diffexp <- snr
#     res[[i]]$p.value <- pval
#   }
#
#   #snr comparison to other clusters
#   cl.snr <- do.call(rbind, lapply(res, function(x){x$diffexp}))
#   for(i in k) {
#     comp <- apply(cl.snr, 2,
#                   function(x){2 * sum(x[k == i] > x) / (length(x) - 1) - 1})
#     i <- as.character(i)
#     res[[i]]$spec <- comp
#   }
#
#   # Create table
#   lapply(res, function(x){cbind(DiffExp = x$diffexp,
#                                 Spec = x$spec, P.value = x$p.value)})
# }
