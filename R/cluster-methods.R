#' DEF: Find states
#'
#' For details see \code{findStates}
#' @keywords internal
#' @import Biobase
#' @importFrom graphics par plot points abline axis box text
#' @importFrom maptree clip.clust draw.clust
#' @importFrom cba order.optimal
#' @importFrom dendextend rotate
#' @author Daniel C. Ellwanger
.findStates_def <- function(x, link.method="ward.D2", min.size,
                            max.pval=1e-4, min.fc=2, min.g=5, show.plots=FALSE,
                            reverse=FALSE, verbose=FALSE) {

  X <- t(x[x@useFeature, ])
  d <- stats::dist(CellTrails::latentSpace(x))
  ordi <- CellTrails::latentSpace(x)
  if(is.null(rownames(ordi))){
    rownames(ordi) <- sampleNames(x)
  }

  res <- list()

  # Computes average linkage distance between clusters
  f.average_linkage_dist <- function(a, b) {
    sum(as.matrix(dist(a, b))) / (length(a) * length(b))
  }

  n <- attributes(d)$Size

  #find maximal fragmentation
  hcl <- stats::hclust(d = d, method = link.method)
  k <- 0
  ms <- rep(NA, n)
  ms[1] <- n
  for(i in (seq_len(n-2+1)+1)) {
    clsize <- table(stats::cutree(hcl, k = i))
    ms[i] <- min(clsize)
    if(ms[i] > min.size) {
      k <- i
    } else {
      break
    }
  }
  message("Initialized ", k, " clusters with a minimum size of ",
          min.size, " samples each.")

  #quit function if no fragmentation for given params was found
  if(k == 0) {
    stop("No maximal defragmentation found. Choose smaller min_size.")
  }

  #merge cluster based on differential expression
  rownames(X) <- seq_len(n)

  clclip <- clip.clust(hcl, X, k = k) #prune dendrogram
  #plot(clclip)
  #cophenetic distance between cluster
  cldist <- as.matrix(cophenetic(clclip))
  diag(cldist) <- NA
  #if(verbose) {
  #  message("Initialized", dim(cldist)[1], "clusters.")
  #}

  message("Performing post-hoc test ...")
  change <- TRUE
  cl <- rep(NA, n)
  for(i in seq(k)) {
    cl[as.numeric(clclip$membership[[i]])] <- i
  }
  mrgtree <- clclip$merge

  merge.log <- list()
  log.index <- 1
  merge.log2 <- seq(nrow(cldist))
  while(change) {
    change <- FALSE ###NEW
    clorder <- as.numeric(names(sort(table(cl))))
    for(i in clorder) {
      if(verbose) {
        #print(table(cl)) ###BUGFIX
      }

      cld <- cldist[i, ]
      j <- which(cld == min(cld, na.rm = TRUE))

      if(length(j) > 1) { #multiple direct neighbors; no join possible
        next
      } else { #check differential expression
        pvals <- rep(NA, ncol(X))
        fc <- rep(NA, ncol(X))
        for(l in seq(ncol(X))) {
          de <- .diffExpr(X[cl == i, l], X[cl == j, l])
          pvals[l] <- de$p.value
          fc[l] <- median(X[cl == i, l]) - median(X[cl == j, l])
          #suppressWarnings(pvals[l] <- wilcox.test(X[cl == i, l], X[cl == j, l])$p.value)
        }
        pvals[is.na(pvals)] <- 1
        fc[is.na(fc)] <- 0
        pvals <- p.adjust(pvals, method = "fdr")
        #cat(sum(pvals < max.pval), " ", sum(abs(fc) > min.fc), "\n")

        if(sum(pvals < max.pval & abs(fc) > min.fc) < min.g) { #no difference between clusters => join
          cl[cl == i] <- j
          cldist[i, ] <- cldist[, i] <- NA

          if(verbose) {
            cat("Merged ", i, " in ", j, "\n")
          }

          pos <- which(apply(mrgtree, 1, function(x){sum(x %in% c(-i, -j)) == 2}))

          f <- apply(mrgtree, 1, function(x){pos %in% x})
          mrgtree[f, ] <- c(-j, mrgtree[f, ][mrgtree[f, ] != pos])
          mrgtree[pos, ] <- c(NA, NA)

          change <- length(unique(cl)) > 1 #continue only if > 2 clusters remain
          merge.log[[log.index]] <- c(i, j) #merged i in j
          log.index <- log.index + 1
          break;
        }
      }
    }
  }

  #update merge tree
  for(i in seq(nrow(mrgtree))) {
    z <- mrgtree[i, ]
    if(is.na(z[1])) {
      next
    } else {
      if(z[1] > 0) {
        z[1] <- z[1] - sum(is.na(mrgtree[seq_len(z[1]), 1]))
      }
      if(z[2] > 0) {
        z[2] <- z[2] - sum(is.na(mrgtree[seq_len(z[2]), 1]))
      }
      mrgtree[i, ] <- z
    }
  }
  mrgtree <- mrgtree[!is.na(mrgtree[,1]), , drop=FALSE]

  hcd.new.order <- hcl$order
  hcd <- stats::as.dendrogram(hcl)

  #  dmat <- matrix(ncol = ncl, nrow = ncl)
  if(nrow(mrgtree) > 2) { #more than 2 clusters left
    #calculate optimal ordering of clusters
    lbls <- sort(unique(cl))
    ncl <- length(lbls)
    #bari <- matrix(NA, ncol = ifelse(is.null(ordi), ncol(X), ncol(ordi)), nrow = ncl)
    dmat <- matrix(ncol = ncl, nrow = ncl)
    cl <- (-cl)
    for(i in seq(ncl)) {
      cl[cl == -lbls[i]] <- i
      mrgtree[mrgtree == -lbls[i]] <- -i
      #if(is.null(ordi)) {
      #   bari[i, ] <- med(X[cl == i, ], method = "Spatial")$median  #apply(X[cl == i, ], 2L, mean, na.rm = T)
      # } else {
      #   bari[i, ] <- med(ordi[cl == i, ], method = "Spatial")$median #apply(ordi[cl == i, ], 2L, mean, na.rm = T)
      # }

      if(link.method == "average") {
        for(j in seq(ncl)) {
          if(is.null(ordi)) {
            dmat[i, j] <- dmat[j, i] <- f.average_linkage_dist(a=X[cl == i,],
                                                               b=X[cl == j,])
          } else {
            dmat[i, j] <- dmat[j, i] <- f.average_linkage_dist(a=ordi[cl == i, ],
                                                               b=ordi[cl == j, ])
          }
        }
      }
    }

    if(link.method != "average") {
      #amat <- as.matrix(dist(f.adjust_cluster_center(ordi, cl)))
      barycenters <- as.matrix(dist(aggregate(ordi, list(cl), mean)[, -1]))
      #g <- igraph::graph_from_adjacency_matrix(amat, mode = "undirected", weighted = T)
      #res$mst <- igraph::mst(g)
      #dmat <- distances(igraph::mst(g))
      #res$mcentres <- t(sapply(seq(max(cl)), function(i){med(ordi[cl == i, ], method = "Spatial")$median}))
      dmat <- as.matrix(stats::dist(barycenters))
    }

    #compute minimal error between neighbors
    op <- order.optimal(as.dist(dmat), mrgtree)
    #res$dmat <- dmat
    #res$mcentres <- res$mcentres[op$order, ]

    #swap branches in dendrogram
    hcd.new.order <- c()
    for(i in op$order) {
      cells <- which(cl == i)
      cells <- cells[order(hcl$order[cells])]
      hcd.new.order <- c(hcd.new.order, cells)
    }
    hcd <- rotate(hcd, order = as.character(hcd.new.order)) #rownames(ordi)[hcd.new.order]
  }

  # Bugfixing
  #res$mrgtree <- mrgtree

  # Rename clusters (start with 1)
  cl.new <- cl
  for(i in seq(max(cl))) {
    cl.new[cl == unique(cl[hcd.new.order])[i]] <- i
  }
  cl <- cl.new

  #If reverse order is required
  if(reverse) {
    hcd.new.order <- rev(hcd.new.order)
    hcd <- rotate(hcd, order = as.character(hcd.new.order))
    if(length(merge.log) > 0) {
      tmp <- max(unlist(merge.log))
      merge.log <- lapply(merge.log, function(x){tmp - x + 1})
    }
    tmp <- abs(mrgtree[mrgtree < 0])
    tmp <- max(tmp) - tmp + 1
    mrgtree[mrgtree < 0] <- tmp * (-1)
    cl <- max(cl) - cl + 1

    #res$mcentres <- res$mcentres[rev(seq(max(nrow(mcentres)))), ]
  }

  f.find_cluster_plots <- function(i) {
    if(i == 1) {
      par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.65, 0))
      plot(seq_len(k), ms[seq_len(k)], ylab = expression("<" * italic(S) * ">"),
           type = "b",
           xlab = "Number of clusters", col = "blue", axes = FALSE,
           ylim = c(0, max(ms[seq_len(k)])))
      #f.grid(x.grid = seq_len(k))
      points(seq_len(k), ms[seq_len(k)], type = "b", col = "blue")
      abline(h=min.size, lty = 2)
      axis(1, at = seq_len(k)); axis(2); box(bty = "L")
    } else if(i == 2) {
      par(mar = c(0.5, 3, 0.5, 0.5), mgp = c(2, 0.65, 0))
      plot(as.hclust(hcl), main = "", xlab = "", sub = "",
           labels = rep("", length(hcl$order)))
      rect.hclust(as.hclust(hcd), k)
    } else if(i == 3) {
      par(mar = c(0.5, 3, 0.5, 0.5), mgp = c(2, 0.65, 0))
      if(k > 2) {
        draw.clust(clip.clust(as.hclust(hcl), X, k=k), size=1.5,
                   pch=NA, cex=.1)
        #axis(2, font = 2, las = 2)
      }
    } else if(i == 4) {
      par(mar = c(0, 1, 0.5, 0.5), mgp = c(2, 0.65, 0))
      n <- length(merge.log)
      plot(0, 0, type = "n", ylim = c(-2, n), xlim = c(0, 10), axes=FALSE,
           xlab = "", ylab = "")
      text(1, n - 0.5, "Log", adj = c(0,0), font = 2)
      if(n > 0) {
        for(i in seq(n)) {
          text(1, n - i, paste(i, ". Merged: ", merge.log[[i]][1], " & ",
                               merge.log[[i]][2], " => ",
                               merge.log[[i]][2], sep = ""), adj = c(0,0))
        }
      }
      text(1, -1, paste("Number of clusters found:", max(cl)), adj = c(0,0))
    }
  }

  if(show.plots) {
    par(mfrow = c(2,2), mar = c(4, 4, 0.5, 0.5), mgp = c(2, 0.65, 0))
    for(i in seq_len(4)){
      f.find_cluster_plots(i)
    }
  }

  # Result
  res$cl <- cl
  res$order <- hcd.new.order
  res$log <- merge.log
  res$hcd <- hcd
  res

  # Update object
  states <- paste0("S", cl)
  s <- unique(states)
  o <- order(as.numeric(substring(s, 2)))
  phenoData(x)$STATE <- factor(states, levels = s[o])

  message("Found ", length(s), " states.")

  #factor(paste0("S", cl), levels = paste0("S", seq(max(cl))))
  x
}
