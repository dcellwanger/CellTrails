#' Adjusting the trajectory graph layout
#'
#' Adjusts the edge distances in the trajectory graph
#' such that edge distances correlates to edge weights;
#' pseudotime is stored in edge weights
#' @param x A \code{\link{CellTrailsSpectrum}} object
#' @return An adjusted layout
#' @importFrom igraph distances mst get.edgelist graph_from_adjacency_matrix
#' @importFrom igraph get.shortest.paths components delete.vertices
#' @keywords internal
#' @author Daniel C. Ellwanger
.adjustLayoutByPtime <- function(x) {
  g <- x@trajectory$traj
  l <- as.matrix(trajectoryLayout(x))
  # Rescale
  yrange <- which.max(apply(l, 2L, function(x) diff(range(x))))
  yrange <- range(l[, yrange])
  if(min(yrange) < 0) {
    yrange <- yrange + abs(min(yrange))
  }
  l <- apply(l, 2L, .rescale, ymin=yrange[1], ymax=yrange[2])

  l_new <- l
  D <- igraph::distances(g)

  # 1. Find internal nodes and leafs
  inodes <- which(degree(g) > 2)
  lnodes <- which(degree(g) == 1)

  if(length(inodes) == 0) { #no internal nodes
    inodes <- lnodes[1] #designate one leaf as internal node
    lnodes <- lnodes[2]
  }

  # 2. Find pairs of internal nodes
  if(length(inodes) > 1) {
    D_inodes <- D[inodes, inodes, drop=FALSE]
    ginodes <- igraph::mst(graph_from_adjacency_matrix(D_inodes, weighted=TRUE, mode="undirected"))
    ge <- get.edgelist(ginodes)
    for(i in seq_len(nrow(ge))) {
      pth <- as.vector(get.shortest.paths(g, from=ge[i, 1], to=ge[i, 2])$vpath[[1]])
      for(j in seq_len(length(pth) - 1)) {
        npair <- pth[j:(j+1)]
        d <- max(dist(l_new[npair, ])[1], 1e-7) #avoiding error due to overlapping nodes
        eps <- D[npair[1], npair[2]] / d
        va <- l_new[npair[1], ]
        vb <- l_new[npair[2], ]
        vc <- va - vb
        vb_new <- va + eps * vc #in one line: va + (eps * (va - vb)) / sum((va - vb)^2)
        l_new[npair[2], ] <- vb_new
      }
      delta <- (vb_new - vb)
      # in last step, push subgraph
      g_tmp <- delete.vertices(g, v=pth[length(pth)])
      comp <- igraph::components(g_tmp)$membership
      compl <- as.numeric(names(comp[comp != comp[as.character(pth[1])]]))
      l_new[compl,] <- sweep(l_new[compl, ], MARGIN=2, STATS=delta, FUN = "+")
    }
  }

  # 3. Foreach internal node find its leafs
  #D <- igraph::distances(g)
  lnode2inode <- cbind(Leaf=lnodes, Internal=inodes[apply(D[lnodes, inodes, drop=FALSE], 1, which.min)])

  # 3. Run from internal node to leafs and update distances (relative)
  for(i in seq_along(lnodes)) {
    from_to <- lnode2inode[i, ]
    pth <- as.vector(get.shortest.paths(g, from=from_to["Internal"], to=from_to["Leaf"])$vpath[[1]])
    for(j in seq_len(length(pth) - 1)) {
      npair <- pth[j:(j+1)]
      d <- max(dist(l_new[npair, ])[1], 1e-7) #avoiding error due to overlapping nodes
      eps <- D[npair[1], npair[2]] / d
      va <- l_new[npair[1], ]
      vb <- l_new[npair[2], ]
      vc <- va - vb
      vb_new <- va + eps * vc
      l_new[npair[2], ] <- vb_new
    }
  }
  if(diff(range(l[,1])) == 0 | diff(range(l[,2])) == 0) {
    #nothing to do; no layout existent
  } else {
    l_new[,1] <- .rescale(x=-l_new[,1], ymin=min(l[,1]), ymax=max(l[,1]))
    l_new[,2] <- .rescale(x=-l_new[,2], ymin=min(l[,2]), ymax=max(l[,2]))
  }
  data.frame(l_new)
}
