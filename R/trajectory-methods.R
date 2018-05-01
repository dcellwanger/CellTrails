#' DEF: Connect states
#'
#' For details see \code{connectStates}
#' @importFrom igraph components induced.subgraph graph_from_adjacency_matrix
#' @keywords internal
#' @author Daniel C. Ellwanger
.connectStates_def <- function(x, l=10, sigma=1){

  f.letters <- function(x) {
    n <- ceiling(log(max(x), 26))

    if(n < 2) {
      return(letters[x])
    } else {
      l <- list()
      for(i in seq(n)) {
        l[[i]] <- letters
      }
      return(apply(expand.grid(l), 1, paste0, collapse = "")[x])
    }
  }

  f.getStates <- function(cells.clust) {
    tmp <- lapply(cells.clust,
                  function(k){d <- D[k, ]; o <- order(d)[seq(l)+1]; d[o]})
    tmp <- unlist(tmp)
    names(tmp) <- cl[as.numeric(names(tmp))]
    tmp
  }

  f.DFS <- function(v, A, visited, parent) {
    visited[v] <- visited[v] + 1
    for(i in seq_len(ncol(A))) {
      if(visited[i] == 0 & A[v, i] > 0) {
        visited <- f.DFS(v = i, A = A, visited = visited, parent = v)
      } else if(visited[i] == 1 & A[v, i] > 0 & i != parent) {
        visited[i] <- Inf
      }
    }
    visited
  }

  dmat <- stats::dist(CellTrails::latentSpace(x))
  cl <- CellTrails::states(x)
  sts <- levels(cl)

  lnn <- list()
  D <- as.matrix(dmat)
  n <- length(sts) #levels(cl)
  #lttrs <- f.letters(seq(n))
  result <- list()

  # Determine nearest neighbors
  lnn <- lapply(sts, function(i) {
    cells.clust <- which(cl == i)
    f.getStates(cells.clust)
    })
  names(lnn) <- sts

  # Filter outliers
  result$lnn <- lnn
  fltr <- exp(median(log(unlist(lnn))) + sigma * mad(log(unlist(lnn))))
  #cat(exp(median(log(unlist(lnn))) + sigma * mad(log(unlist(lnn)))), "\n")
  lnnf <- lapply(lnn, function(x) x[x < fltr])
  #lnnf <- lapply(tmp$lnn, function(x){(x[x < quantile(x, quant)])})
  #lnnf <- lapply(tmp$lnn, function(x){(x[x < median(x) + sigma * mad(x)])})
  ncnt <- lapply(lnnf, function(x){table(names(x))})
  ncnt <- lapply(ncnt, function(x) x / sum(x))

  #print(ncnt) #DEBUG

  result$ncnt <- ncnt

  ### Compute maximum overlap spanning forest (DFS to find circles)
  #dt <- data.frame(do.call(rbind, strsplit(names(unlist(ncnt)), "\\.")), unlist(ncnt), stringsAsFactors = F)
  dt <- data.frame(FROM = rep(names(ncnt), lapply(ncnt, length)),
                   TO = unlist(lapply(ncnt, names)),
                   CNT = as.numeric(unlist(ncnt)), stringsAsFactors=FALSE)

  dt <- dt[!dt[,1] == dt[,2], ]
  dt <- dt[order(dt[,3], decreasing=TRUE), ]

  result$dt <- dt

  amat <- matrix(0, ncol = n, nrow = n)
  colnames(amat) <- rownames(amat) <- sts #levels(cl) #lttrs[seq_len(max(cl))]
  for(i in seq_len(nrow(dt))) {
    amat.new <- amat
    amat.new[dt[i,1], dt[i,2]] <- amat.new[dt[i,2], dt[i,1]] <- 1
    s <- 0
    for(k in seq(n)) {
      s <- s + sum(f.DFS(k, A = amat.new, rep(0, nrow(amat.new)), parent = -1))
    }
    if(!is.infinite(s)) {
      amat <- amat.new
    #  cat(dt[i,1], dt[i,2], "\n")
    }
    #    else {
    #      print(dt[i,seq_len(2)])
    #    }
    if(sum(amat) / 2 == ncol(amat) - 1) {
      i <- l + 1
      break
    }
  }
  #colnames(amat) <- rownames(amat) <- seq(ncol(amat))
  g <- graph_from_adjacency_matrix(amat, mode = "undirected")

  result$amat <- amat

  # Update object
  # Check for components
  maxIfForest <- list()
  cmps <- igraph::components(g)
  for(i in seq(cmps$no)) {
    maxIfForest[[i]] <- induced.subgraph(g, V(g)[cmps$membership == i])
  }
  #result$maxIfForest <- maxIfForest
  #result

  x@spanForest <- maxIfForest
  if(length(x@spanForest) == 1) {
    x@useSample <- rep(TRUE, ncol(x))
  }
  x
}

#' DEF: Select component from trajectory graph
#'
#' For details see \code{selectTrajectory}
#' @importFrom igraph V
#' @keywords internal
#' @author Daniel C. Ellwanger
.selectTrajectory_def <- function(x, component) {
  x@spanForest <- list(x@spanForest[[component]])
  vnames <- names(V(x@spanForest[[1]]))

  #Filter CellTrailsSet object
  x@useSample <- states(x) %in% vnames

  #fltr <- states(x) %in% vnames
  #x@assayData <- x@assayData[fltr, ]
  #x@phenoData <- x@phenoData[fltr, ]
  #x@eigenspace <- x@eigenspace[fltr, ]
  x
}

#' DEF: Trajectory fitting
#'
#' For details see \code{selectTrajectory}
#' @import Biobase
#' @importFrom igraph degree V<-
#' @keywords internal
#' @author Daniel C. Ellwanger
.fitTrajectory_def <- function(x) { #use.comp = F
  orth <- .project_ortho(x)
  orth_graph <- .connect_ortho(orth)

  #Check if projection needs to be expanded
  if(.needsToBeExpanded(x, orth_graph)) {
    orth_graph <- .deleteMedianCentres(orth_graph)
  } else {
    #Refine and remove mediancentres
    vizmap <- .generate_ordination(g = orth_graph,
                                   error = orth$error,
                                   only.ordi = TRUE)
    orth_graph <- .connect_ordi(vizmap$ordi)
  }

  #Update x
  dg <- degree(orth_graph)
  hps <- dg == 1 #trail heads
  bps <- dg > 2  #bifurcations

  blaze.type <- rep(NA, length(dg))
  blaze.type[hps] <- "H"
  blaze.type[bps] <- "B"
  blaze.type <- factor(blaze.type, levels = c("H", "B", "U"))

  blaze.id <- rep(NA, length(dg))
  blaze.id[hps] <- paste0("H", seq(sum(hps)))
  blaze.id[bps] <- paste0("B", seq(sum(bps)))

  V(orth_graph)$sampleName <- sampleNames(x)[x@useSample]

  res <- list()
  res$error <- orth$error
  res$traj <- orth_graph
  res$blaze <- data.frame(type=blaze.type, id=blaze.id,
                          shape=rep("ellipse", length(blaze.type)),
                          stringsAsFactors=FALSE)
  x@trajectory <- res
  x
}

#' AUX: Orthogonal projection
#'
#' Orthogonally projects samples onto trajectory
#' @param x An \code{CellTrailsSet} object
#' @return A list containing the following components:
#' @return \item{\code{Y}}{Projection}
#' @return \item{\code{X.c}}{Mediancentres}
#' @return \item{\code{lambda}}{Projection lambda (position: inside/outside of edge)}
#' @return \item{\code{adjmat}}{Adjacency matrix with sample count per edge}
#' @return \item{\code{error}}{Projection error per sample}
#' @return \item{\code{edgeId}}{ID of each edge}
#' @return \item{\code{edgeId2Edge}}{Mapping of edge to edgeId}
#' @return \item{\code{edge}}{For each sample its assigned edge}
#' @importFrom igraph V neighbors
#' @importFrom utils combn
#' @importFrom depth med
#' @keywords internal
#' @author Daniel C. Ellwanger
.project_ortho <- function(x) {
  cl <- droplevels(states(x)[x@useSample])
  g <- stateTrajectoryGraph(x)[[1]]
  X <- CellTrails::latentSpace(x)[x@useSample, ]
  X.c <- t(vapply(levels(cl), function(i){med(X[cl == i, ],
                                              method = "Spatial")$median},
                  rep(0, dim(X)[2])))
  mx <- max(as.numeric(substr(rownames(X.c), 2, nrow(X.c))))
  cl <- as.character(cl) #avoid automatic factor to number conversion

  #if(use.comp) {
  #  tmp.d <- as.matrix(dist(X, gng$codes))
  #  cl <- apply(tmp.d, 1, which.min)
  #}

  # Orthogonal projection of X on line through P1 and P2
  # I.e. vector projection of P1 on X with lambda = component and P2 = <0...0>
  f.projectXtoLine <- function(X, P1, P2) {
    result <- list()
    u <- P2 - P1
    lambda <- ((X - P1) %*% u) / (u %*% u)

    result$lambda <- lambda
    result$Y <- P1 + lambda %*% u
    result
  }

  Y <- matrix(NA, ncol = ncol(X), nrow = nrow(X))
  lambda <- rep(NA, ncol(X))
  error <- rep(NA, nrow(X))
  edgeId <- rep(NA, nrow(X))

  nCodes <- dim(X.c)[1]
  op_adjmat <- matrix(0, ncol=nCodes, nrow=nCodes,
                      dimnames=list(rownames(X.c), rownames(X.c)))

  for(i in seq_len(nrow(X))) {
    instance <- X[i, ]
    n <- names(neighbors(g, V(g)[cl[i]]))

    min_dist <- Inf      #Euclidean distance (fitting error)
    min_op <- list()     #projection information
    min_ids <- c(NA, NA) #ids of mediancentres of states
    best_dist <- Inf
    best_op <- list()
    best_ids <- c(NA, NA)
    for(j in n) { #project on closest straight line
      op <- f.projectXtoLine(instance, X.c[cl[i], ], X.c[j, ])
      op_dist <- sum((X[i,] - op$Y)^2)

      if(op_dist < best_dist & (op$lambda >= 0 & op$lambda <= 1)) {
        best_op <- op
        best_dist <- op_dist
        best_ids <- c(cl[i], j)
      }

      if(op_dist < min_dist) {
        min_op <- op
        min_dist <- op_dist
        min_ids <- c(cl[i], j)
      }
    }
    # Store projection information
    if(length(best_op) == 0) {
      Y[i, ] <- min_op$Y
      lambda[i] <- min_op$lambda
      up <- op_adjmat[min_ids[1], min_ids[2]] + 1
      op_adjmat[min_ids[1], min_ids[2]] <- op_adjmat[min_ids[2], min_ids[1]] <- up
      error[i] <- min_dist
      min_ids_num <- as.numeric(substr(min_ids, 2, nrow(X.c)))
      edgeId[i] <- (min(min_ids_num) - 1) * mx - (min(min_ids_num) * (min(min_ids_num) + 1) / 2) + max(min_ids_num)
    } else {
      Y[i, ] <- best_op$Y
      lambda[i] <- best_op$lambda
      up <- op_adjmat[best_ids[1], best_ids[2]] + 1
      op_adjmat[best_ids[1], best_ids[2]] <- op_adjmat[best_ids[2], best_ids[1]] <- up
      error[i] <- best_dist
      best_ids_num <- as.numeric(substr(best_ids, 2, nrow(X.c)))
      edgeId[i] <- (min(best_ids_num) - 1) * mx - (min(best_ids_num) * (min(best_ids_num) + 1) / 2) + max(best_ids_num)
    }
  }

  # Mapping: EdgeId to edge
  res <- list()
  res$Y <- Y
  res$X.c <- X.c
  res$lambda <- lambda
  res$adjmat <- op_adjmat
  res$error <- error
  res$edgeId <- as.character(edgeId)
  res$edgeId2Edge <- t(combn(rownames(X.c), 2)) #t(combn(as.numeric(substr(rownames(X.c), 2, nrow(X.c))), 2)) #t(combn(seq(mx), 2))
  rownames(res$edgeId2Edge) <- apply(res$edgeId2Edge, 1, function(x){
    x <- as.numeric(substr(x, 2, nrow(X.c)));
    (min(x) - 1) * mx - (min(x) * (min(x) + 1) / 2) + max(x)})

  res$edge <- res$edgeId2Edge[res$edgeId, ]
  res
}

#' AUX: Connect orthogonal projected samples
#'
#' Connects samples based on their projected position on the trajectory.
#' @param Xorth Result from the orthogonal projection
#' @return An \code{igraph} object
#' @importFrom igraph graph_from_adjacency_matrix
#' @keywords internal
#' @author Daniel C. Ellwanger
.connect_ortho <- function(Xorth) {
  eids <- unique(Xorth$edgeId)
  xcodes <- Xorth$X.c
  dname <- c(seq(nrow(Xorth$Y)), rownames(xcodes))
  amat <- matrix(0, nrow = nrow(Xorth$Y) + nrow(xcodes),
                 ncol = nrow(Xorth$Y) + nrow(xcodes),
                 dimnames = list(dname, dname)) #adjacency between cells

  # Link all cells embedded on line between two coding vectors
  for(i in eids) {
    cells <- which(Xorth$edgeId == i)
    fromTo <- Xorth$edgeId2Edge[i, ]
    cells.in <- cells[Xorth$lambda[cells] >= 0]
    cells.out <- cells[Xorth$lambda[cells] < 0]

    # 1. Link cells located between two coding vectors
    # 1.1 Calculate distance to one of both coding vectors
    d1 <- apply(Xorth$Y[cells.in, , drop = FALSE], 1L,
                function(x){sum((x - xcodes[fromTo[1], ])^2)})
    o1 <- cells.in[order(d1, decreasing = FALSE)]
    d2 <- apply(Xorth$Y[cells.in, , drop = FALSE], 1L,
                function(x){sum((x - xcodes[fromTo[2], ])^2)})
    o2 <- cells.in[order(d2, decreasing = FALSE)]

    # 1.2 Connect cells
    for(j in seq(length(o1) - 1)){
      j1 <- as.character(o1[j])
      j2 <- as.character(o1[j + 1])
      amat[j1, j2] <- amat[j2, j2] <- 1
    }

    # 1.3 Connect cells to coding vector
    amat[as.character(o1[1]), fromTo[1]] <- amat[fromTo[1], as.character(o1[1])] <- 1
    amat[as.character(o2[1]), fromTo[2]] <- amat[fromTo[2], as.character(o2[1])] <- 1

    #amat[o1[1], nrow(Xorth$Y) + fromTo[1]] <- amat[nrow(Xorth$Y) + fromTo[1], o1[1]] <- 1
    #amat[o2[1], nrow(Xorth$Y) + fromTo[2]] <- amat[nrow(Xorth$Y) + fromTo[2], o2[1]] <- 1

    # 2. Link cells located outside the two coding vectors
    d1 <- apply(Xorth$Y[cells.out, , drop = FALSE], 1L,
                function(x){sum((x - xcodes[fromTo[1], ])^2)})
    d2 <- apply(Xorth$Y[cells.out, , drop = FALSE], 1L,
                function(x){sum((x - xcodes[fromTo[2], ])^2)})
    cells.out.1 <- as.vector(cells.out[d1 < d2])
    cells.out.2 <- as.vector(cells.out[d1 > d2])

    # 2.1 Link cells located outside of vector 1
    if(length(cells.out.1) > 1) {
      d1 <- apply(Xorth$Y[cells.out.1, with = TRUE], 1L,
                  function(x){sum((x - xcodes[fromTo[1], ])^2)})
      o1 <- cells.out.1[order(d1, decreasing = FALSE)]
      for(j in seq(length(o1) - 1)){
        j1 <- as.character(o1[j])
        j2 <- as.character(o1[j + 1])
        amat[j1, j2] <- amat[j2, j1] <- 1
      }
      # Link cells to coding vectors
      amat[as.character(o1[1]), fromTo[1]] <- amat[fromTo[1], as.character(o1[1])] <- 1
    } else if(length(cells.out.1) == 1) {
      # Link cell to coding vector
      amat[as.character(cells.out.1), fromTo[1]] <- amat[fromTo[1], as.character(cells.out.1)] <- 1
    }

    # 2.2 Link cells located outside of vector 2
    if(length(cells.out.2) > 1) {
      d2 <- apply(Xorth$Y[cells.out.2, ], 1, function(x){sum((x - xcodes[fromTo[2], ])^2)})
      o2 <- cells.out.2[order(d2, decreasing = FALSE)]
      for(j in seq(length(o2) - 1)){
        j1 <- as.character(o2[j])
        j2 <- as.character(o2[j + 1])
        amat[j1, j2] <- amat[j2, j1] <- 1
      }
      # Link cells to coding vectors
      amat[as.character(o2[1]), fromTo[2]] <- amat[fromTo[2], as.character(o2[1])] <- 1
    } else if(length(cells.out.2) == 1) {
      # Link cell to coding vector
      amat[as.character(cells.out.2), fromTo[2]] <- amat[fromTo[2], as.character(cells.out.2)] <- 1
    }
  }
  w <- matrix(1e-10, nrow=nrow(amat), ncol=ncol(amat))
  w[seq(nrow(Xorth$Y)), seq(nrow(Xorth$Y))] <- as.matrix(dist(Xorth$Y))
  amat <- amat * w
  graph_from_adjacency_matrix(amat, mode="undirected", weighted=TRUE)
}

#' AUX: Orthogonal projection ordination
#'
#' Generates 2D ordination from orthogonal projection
#' @param g Graph from orthogonal projection, igraph object
#' @param error Error from orthogonal projection; numeric vector
#' @param only.ordi Compute only ordination?
#' @return A list containing the following components:
#' @return \item{\code{ordi}}{Ordination}
#' @return \item{\code{ordi.jitter}}{Jittered ordination}
#' @return \item{\code{backbone}}{Backbone of trajectory}
#' @importFrom igraph distances degree
#' @keywords internal
#' @author Daniel C. Ellwanger
.generate_ordination <- function(g, error, factor=7, rev=FALSE,
                                 only.ordi=FALSE) {
  # Offset points
  # x = Target; x.n = neighbors; size = offset width; side = (-1, 1) is left or right, respectively
  f.offset <- function(x, x.n, size, side = 1) {
    #1. Find line through neighbors
    slope <- (x.n[1, 2] - x.n[2, 2]) / (x.n[1, 1] - x.n[2, 1])
    intercept <- x.n[1, 2] - slope * x.n[1, 1]

    #2. Find orthogonal line through x
    slope2 <- -1 / slope
    intercept2 <- x[2] - slope2 * x[1]

    #3. Find any vector orthogonal to line
    v <- NA
    if(side == 1) {
      v <- x[1] + 1
      v <- c(v, v * slope2 + intercept2)  #(e.g. with x1 = 1, x2 according to slope/intercept eq.)
    } else {
      v <- x[1] - 1
      v <- c(v, v * slope2 + intercept2)  #(e.g. with x1 = 1, x2 according to slope/intercept eq.)
    }
    v <- x - v #use sum of vector equation to identify orthogonal vector

    #3. Generate new vector of given size
    # [https://math.stackexchange.com/questions/897723/how-to-resize-a-vector-to-a-specific-length]
    v <- size / sqrt(sum(v^2)) * v #use vector norm to scale it
    v + x #sum of vectors (scaled vector + target point) give new point
  }

  #Result
  res <- list()
  res$ordi <- NA
  res$ordi.jitter <- NA
  res$backbone <- NA

  #Find ordination
  d <- igraph::distances(g)
  rng <- sort(which(d == max(d), arr.ind=TRUE)[1, ])
  if(rev) {
    rng <- rev(rng)
  }
  ordi <- igraph::distances(g)[, rng]
  ordi <- apply(ordi, 2, .rescale, ymin = 0, ymax = 1)

  #Remove mediancentres
  ordi <- ordi[seq(length(error)), ]

  if(!only.ordi) {
    # Add error
    error <- error / max(error) / factor
    d <- as.matrix(dist(ordi))
    eordi <- matrix(NA, nrow = nrow(ordi), ncol = ncol(ordi))
    for(i in seq(nrow(ordi))) {
      d.tmp <- round(d[i, ], 5) #avoid float error
      o <- order(d.tmp, decreasing = FALSE)
      o <- o[2:length(o)]
      nhbr <- c(o[1])
      cnt <- 2
      while(length(nhbr) < 2) {
        curr <- o[cnt]
        if(sum(abs(ordi[nhbr, ] - ordi[curr, ])) > 1e-5) {
          nhbr <- c(nhbr, curr)
        }
        cnt <- cnt + 1
      }
      #d.tmp <- d.tmp[d.tmp > 0]
      #nhbr <- order(d.tmp, decreasing = F)[seq_len(2)]
      eordi[i, ] <- f.offset(ordi[i, ], ordi[nhbr, ], size=error[i], side=1)
    }

    #Reverse ordination
    #if(rev) {
    #  #ordi <- apply(-ordi, 2, .rescale, ymin = 0, ymax = 1)
    #  #eordi <- apply(-eordi, 2, .rescale, ymin = 0, ymax = 1)
    #  ordi[,1] <- max(ordi[,1]) - ordi[,1]
    #  ordi[,2] <- max(ordi[,2]) - ordi[,2]
    #  eordi[,1] <- max(eordi[,1]) - eordi[,1]
    #  eordi[,2] <- max(eordi[,2]) - eordi[,2]
    #}

    res$ordi.jitter <- eordi
    terms <- which(degree(g) == 1)
    #backbone <- get.all.shortest.paths(g, from = terms[1], to = terms[-1])$res
    #res$backbone <- lapply(backbone, function(x) {x <- x[x %in% seq(length(error))]; ordi[x, ]}) #Remove mediancentres
    #res$backbone <- unique(do.call(rbind, Y$backbone))
  }

  # Result
  res$ordi <- ordi
  res
}

#' AUX Check if projection needs to be expanded
#'
#' Simple graphs having only bifurcations along the backbone get
#' simplified by a 2D projection. If a fork with more than two successors
#' or a fork not located along the backbone, but along a side branch, the
#' 2D simplification would collate the states and ignore the bifurcation.
#' Therefore, the simplification step has to be skipped. This methods
#' checks whether a siplification is possible or not.
#' @param ctset An object of class \code{CellTrailsSet}
#' @param g An object of class \code{igraph}
#' @return A logical value
#' @details IF: any node has a degree > 3, return \code{TRUE}; OTHERWISE:
#' Check all sorthest paths from start node to any leaf in trajectory
#' tree. IF any path passes a node with degree > 3, which is not located
#' on the backbone, return \code{TRUE}; OTHERWISE: return \code{FALSE}.
#' @importFrom igraph distances get.shortest.paths degree
#' @keywords internal
#' @author Daniel C. Ellwanger
.needsToBeExpanded <- function(x, g) {
  D <- igraph::distances(g)
  r <- which(D == max(D), arr.ind=TRUE)
  r.s <- as.character(states(x)[which(x@useSample)[r]]) #backbone in trajectory graph
  tg <- stateTrajectoryGraph(x)[[1]]
  bb <- names(get.shortest.paths(tg, from = r.s[1], to = r.s[2])$vpath[[1]])
  dg <- degree(tg)
  if(any(dg > 3)) {
    return(TRUE)
  } else {
    l <- names(which(dg == 1))
    l <- l[!l %in% r.s]

    doExpand <- vapply(l, FUN = function(z) {
      n <- names(get.shortest.paths(tg, from = r.s[1], to = z)$vpath[[1]])
      any(any(dg[n[!n %in% bb]] > 3))
    }, TRUE)
    return(any(doExpand))
  }
}

#' AUX Remove median centres
#'
#' Deletes median centres from trajectory graph
#' if graph has not been simplified/has been expanded.
#' @param g An object of class \code{igraph}
#' @return An unpdated object of class \code{igraph}
#' @importFrom igraph distances degree V neighbors add.edges delete.vertices
#' @keywords internal
#' @author Daniel C. Ellwanger
.deleteMedianCentres <- function(g) {
  mcentres <- which(substr(names(V(g)), 1, 1) == "S") #name = 'Sx'
  D <- igraph::distances(g)
  dg <- degree(g, v = V(g)[mcentres])
  edgs <- rep(NA, (sum(dg) - length(dg)) * 2)
  w <- rep(NA, sum(dg) - length(dg))
  cnt <- 0
  for(i in seq_along(mcentres)) {
    mc <- mcentres[i]
    n <- names(neighbors(g, V(g)[mc]))
    for(j in 2:length(n)) {
      cnt <- cnt + 1
      edgs[cnt + cnt - 1] <- n[1]
      edgs[cnt + cnt] <- n[j]
      w[cnt] <- D[n[1], n[j]]
    }
  }
  g <- add.edges(g, edges = edgs, attr = list(weight = w))
  g <- delete.vertices(g, V(g)[mcentres])
  g
}

#' AUX Generate graph based on ordination
#'
#' Generates 2D ordination from orthogonal projection
#' @param ordi Ordination
#' @return An \code{igraph} object
#' @importFrom igraph graph_from_adjacency_matrix
#' @keywords internal
#' @author Daniel C. Ellwanger
.connect_ordi <- function(ordi) { #, geodmat
  # 1. Identify cells on backbone (slope = -1; intercept = 1)
  #o <- order(viz$ordi[,1])
  #v <- viz$ordi[c(o[1], rev(o)[1]),  ]
  #slope <- (v[2, 2] - v[1, 2]) / (v[2, 1] - v[1, 1]) #-1
  #intercept <- v[2, 2] - slope * v[2, 1] #1

  f.isOnLine <- function(x, slope, intercept) {
    round(slope * x[1] + intercept, 5) == round(x[2], 5)
  }

  onBackbone <- apply(ordi, 1L, f.isOnLine, slope = -1, intercept = 1)

  # 2. Identify side-branches (all of them have same slope = 1 => identify intercepts)
  v <- which(!onBackbone)
  branches <- factor(apply(ordi[v, ], 1L, function(x){round(x[2] - x[1], 5)}))

  # 3. Create adjacency matrix
  adjmat <- matrix(0, nrow = nrow(ordi), ncol = nrow(ordi))
  dmat <- as.matrix(stats::dist(ordi)) #geodmat[seq(nrow(adjmat)), seq(nrow(adjmat))] #

  #backbone
  cells <- which(onBackbone)
  start <- cells[which(dmat[cells, cells] == max(dmat[cells, cells]),
                       arr.ind=TRUE)[1, 1]]
  traj <- cells[order(dmat[start, cells])]
  for(i in seq(length(traj) - 1)) {
    adjmat[traj[i], traj[i + 1]] <- adjmat[traj[i + 1], traj[i]] <- 1
  }

  #branches
  for(i in levels(branches)){
    #connect branch cells
    cells <- v[branches == i]
    start <- cells[which(dmat[cells, cells, drop=FALSE] == max(dmat[cells, cells, drop=FALSE]),
                         arr.ind=TRUE)[1, 1]]
    traj <- cells[order(dmat[start, cells])]
    for(i in seq(length(traj) - 1)) {
      adjmat[traj[i], traj[i + 1]] <- adjmat[traj[i + 1], traj[i]] <- 1
    }
    #connect to backbone
    anker <- which(dmat[cells, onBackbone, drop=FALSE] == min(dmat[cells, onBackbone, drop=FALSE]),
                   arr.ind=TRUE)[1, ]
    anker <- c(cells[anker[1]], which(onBackbone)[anker[2]])
    adjmat[anker[1], anker[2]] <- adjmat[anker[2], anker[1]] <- 1
  }

  adjmat <- (dmat + min(min(dmat[dmat > 0]), 1e-10)) * adjmat
  graph_from_adjacency_matrix(adjmat, mode="undirected", weighted=TRUE)
}
