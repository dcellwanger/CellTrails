#' DEF: Visualizing states sizes
#'
#' Plots barplot of number of samples per state
#' @param x States
#' @return A \code{ggplot} object
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.plotStateSize_def <- function(x) {
  count <- data.frame(table(x))
  colnames(count) <- c("State", "Freq")
  ggplot(count, aes_string(x="State", y="Freq", fill="State")) +
    geom_bar(stat="identity") +
    labs(x="State", y="Sample count") +
    guides(fill=FALSE) +
    theme(axis.line=element_line(colour="black"))
}

#' DEF: Visualizing the spectrum definition
#'
#' Plots Scree plot with eigengaps
#' @param x A \code{list}
#' @return A \code{ggplot} object
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.plotSpectrum_def <- function(x) {
  if(x$frac <= 1) {
    x$frac <- x$frac * length(x$cs)
  }
  Y <- x$cs[seq(x$frac)]
  gp <- ggplot()
  gp <- gp + aes(x=seq(x$frac), y=Y) +
    geom_point(color=c(rep("red", x$n - 1), rep("gray40", x$frac - x$n + 1)))
  gp <- gp + geom_line(mapping=aes(x=x$fit$x, y=x$fit$y, lty = 'Fit'),
                       color = "blue")
  gp <- gp + xlab("Spectrum")
  gp <- gp + ylab("Total eigengap")
  brks <- sort(c(pretty(seq(x$frac)), x$n - 1))
  gp <- gp + scale_x_continuous(breaks=brks, labels=brks + 1)
  col <- rep("black", length(brks))
  col[brks == x$n - 1] <- "red"
  gp <- gp + theme(axis.line=element_line(colour = "black"),
                   axis.text.x=element_text(color = col),
                   axis.ticks.x=element_line(color = col),
                   legend.title=element_blank())
  #gp <- gp + geom_vline(xintercept = 8)
  gp
}

#' DEF: Violine plots of genes
#'
#' Plot gene expression per state
#' @param x A expression vector
#' @param sts A states assignment vector
#' @param label Label for legend
#' @return A \code{ggplot} object
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.plotStateExpression_def <- function(x, sts, label) {

  df <- data.frame(State=sts, Expression=x)
  brks <- pretty(df$Expression)
  lbl <- brks
  lbl[lbl == 0] <- "nd"

  ggplot(df, aes_string(x = "State", y = "Expression", fill = "State")) +
    geom_jitter(shape = 16, position = position_jitter(0.2), alpha = .5) +
    geom_violin(trim = TRUE, alpha = .8) +
    ylab(label) +
    scale_y_continuous(breaks = brks, labels = lbl) +
    theme(axis.line = element_line(colour = "black"))
}

#' DEF: Visualizing the manifold
#'
#' Plots approximation of sample embedding in latent space in two dimensions
#' @param X Ordination of data points
#' @param y Values to visualize
#' @param name Label for the figure legend
#' @param g The trajectory graph
#' @param weights States along the trajectory
#' @param axis_label Label prefix for each axis
#' @param colors Either the color ramp is colorized or black/white
#' @param type Type of visualization
#' @param samples_only Does not show the whole surface, but colorizes
#' only the samples
#' @param setND Indicates if a '0' value should be set to 'ND' in the
#' figure legend
#' @param alpha Color param
#' @return A \code{ggplot} object
#' @import ggplot2
#' @importFrom igraph V
#' @keywords internal
#' @author Daniel C. Ellwanger
.plotManifold_def <- function(X, y, name, g=NULL, weights=NULL, axis_label="",
                              colors=c("pretty", "bw"),
                              type=c("raw", "surface.fit", "surface.se"),
                              samples_only=FALSE, setND=FALSE,
                              alpha=.1) {

  if(type=="surface.se") {
    colors <- "bw"
  }
  cRamp <- list("bw"=colorRampPalette(c("white", "black")),
                "pretty"=.prettyColorRamp)[[colors[1]]]

  f.ggplayout <- function(gp, X, nas, name) {
    gp <- gp + aes(x = X[nas, 1], y = X[nas, 2]) +
      geom_point(pch = 1, col = "gray40") +
      labs(fill=name)
    gp
  }

  colnames(X) <- c("D1", "D2")

  show_backbone <- !is.null(g)
  if(show_backbone) { #connect samples by trajectory graph
    sname <- V(g)$sampleName
    edgelist <- get.edgelist(g)
    edgelist <- data.frame(X1=V(g)[edgelist[,1]]$sampleName,
                           X2=V(g)[edgelist[,2]]$sampleName,
                           stringsAsFactors=FALSE)
    edges <- data.frame(X[edgelist[,1], ],
                        X[edgelist[,2], ],
                        row.names = seq(nrow(edgelist)))
    colnames(edges) <- c("X1", "Y1", "X2", "Y2")
  }

  # Init plot
  gp <- ggplot()
  type <- tolower(type)[1]

  if(!is.numeric(y)) { #always raw for non-numeric types
    type <- "raw"
  }

  # Plot types
  if(type == "raw") {
    if(show_backbone) {
      gp <- gp + geom_segment(aes_string(x="X1", y="Y1", xend="X2", yend="Y2"),
                              data=edges, size=0.5, colour="darkred")
    }
    if(is.numeric(y)) {
      nas <- rep(FALSE, length(y))
      #limits <- range(y, na.rm=TRUE)
      #values <- pretty(limits)

      if(setND) { #is expression data
        nas <- y == 0
        #limits[1] <- min(y[!nas])
        #lbls <- gsub(lbls, pattern = "^0", replacement = "nd")
      }
      #gp <- f.ggplayout(gp, X, nas, name) +
      #  geom_point(mapping=aes(x=X[!nas,1], y=X[!nas,2], color=y[!nas])) +
      #  scale_color_gradientn(colours=c(.prettyColorRamp(5, grayStart=FALSE)),
      #                        #label=lbls,
      #                        limits=limits, values=seq(0,1,by=0.2))
      brks <- pretty(range(y, na.rm=TRUE, finite=TRUE), 10)
      lbls <- paste(brks[seq_len(length(brks)-1)],
                    brks[seq_len(length(brks)-1)+1], sep = " - ")
      lbls <- gsub(pattern="^0 -", replacement = "nd -", x = lbls)

      equalSpace <- cut(y, breaks=brks, include.lowest=TRUE, labels=lbls)
      breaks <- rev(levels(equalSpace))
      df <- data.frame(X, equalSpace)

      gp <- f.ggplayout(gp, X, nas, name) +
        geom_point(data=df,
                   mapping=aes_string(x="D1", y="D2", fill="equalSpace"),
                   pch=21) +
        scale_fill_manual(name=name, breaks=breaks, labels=breaks,
                            values=cRamp(length(brks) - 1),
                            drop=FALSE)
    } else {
      nas <- is.na(y)
      gp <- f.ggplayout(gp, X, nas, name) +
        geom_point(mapping=aes(x=X[!nas, 1], y=X[!nas, 2],
                               fill=factor(y[!nas])), pch=21)
    }
  } else if(type == "surface.fit") {
    if(samples_only) { #only trajectory
      fit <- .fit_surface(X, y, npoints=1, weights=weights, knots=10)
      gp <- .plotManifold_def(X=X, y=fit$fit$fitted.values, name=name, g=g,
                              weights=weights,
                              axis_label=axis_label, type="raw",
                              samples_only=samples_only, setND=setND)
    } else { #whole 300x300 grid
      fit <- .fit_surface(X, y, npoints=300, weights=weights, knots=10)
      df <- fit$grid
      brks <- pretty(range(df$x3, na.rm=TRUE, finite=TRUE), 10)
      lbls <- paste(brks[seq_len(length(brks)-1)],
                    brks[seq_len(length(brks)-1)+1], sep = " - ")
      lbls <- gsub(pattern="^0 -", replacement = "nd -", x = lbls)
      df$equalSpace <- cut(df$x3, breaks=brks,
                           include.lowest=TRUE, labels=lbls)
      breaks <- rev(levels(df$equalSpace))
      gp <- gp + geom_tile(data=df,
                           aes_string(x="x1", y="x2", fill="equalSpace")) +
        geom_contour(data=df, aes_string(x="x1", y="x2", z="x3"),
                     color="gray40", alpha=0.5, lty=2, lwd=.5) +
        scale_fill_manual(values=cRamp(length(brks) - 1),
                          name=name, breaks=breaks, labels=breaks)
      if(show_backbone) {
        gp <- gp + geom_segment(aes_string(x="X1", y="Y1",
                                           xend="X2", yend="Y2"),
                                data=edges, size=0.5, colour="darkred")
      }
      gp <- gp + geom_point(data=X, aes_string(x="D1", y="D2"),
                            colour="red", size=1.5, alpha=alpha) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0))
    }
  } else if(type == "surface.se") {
    if(samples_only) { #only trajectory
      fit <- .fit_surface(X, y, npoints=1, weights=weights, knots=10)
      y <- predict(fit$fit, type="response", se.fit=TRUE)$se
      gp <- .plotManifold_def(X=X, y=y, name=paste0(name,"(SE)"), g=g,
                              weights=weights,
                              axis_label=axis_label, type="raw",
                              colors="bw",
                              samples_only=samples_only, setND=setND)
    } else {  #whole 300x300 grid
      fit <- .fit_surface(X, y, npoints=300, weights=weights, knots=10)
      df <- fit$grid
      df$x3 <- fit$fit$se
      brks <- pretty(range(df$x3, na.rm=TRUE, finite=TRUE), 10)
      lbls <- paste(brks[seq_len(length(brks)-1)],
                    brks[seq_len(length(brks)-1)+1], sep = " - ")
      df$equalSpace <- cut(df$x3, breaks=brks,
                           include.lowest=TRUE, labels=lbls)
      breaks <- rev(levels(df$equalSpace))

      gp <- gp + geom_tile(data=df,
                           aes_string(x="x1", y="x2", fill="equalSpace")) +
        geom_contour(data=df, aes_string(x="x1", y="x2", z="x3"),
                     color="gray40", alpha=0.5, lty=2, lwd=.5) +
        scale_fill_manual(values=cRamp(length(brks) - 1),
                          name=paste0(name, "(SE)"),
                          breaks=breaks, labels=breaks)
      if(show_backbone) {
        gp <- gp + geom_segment(aes_string(x="X1", y="Y1",
                                           xend="X2", yend="Y2"),
                                data=edges, size=0.5, colour="darkred")
      }
      gp <- gp + geom_point(data=X, aes_string(x="D1", y="D2"),
                            colour="red", size=1.5, alpha=alpha) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0))
    }
  } else {
    stop("Unkown plot type selected.")
  }

  gp <- gp + theme(axis.line=element_line(colour="black"),
                   panel.border=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank()) +
    xlab(paste0(axis_label, " 1")) +
    ylab(paste0(axis_label, " 2"))
  gp
}

#' DEF: Visualizing the trajectory graph
#'
#' Visualizes the trajectory spanning all states
#' of a component
#' @param g State trajectory graph (igraph object)
#' @param y Values; numeric or factor vector
#' @param all_sts All state assignments
#' @param label Label for plot legend
#' @param setND If values of '0' are getting set to 'ND' in the
#' figure legend
#' @param point_size Point size paramenter
#' @param label_offset Label offset parameter
#' @param seed Makes result reproducible
#' @return A \code{ggplot} object
#' @importFrom igraph V layout.fruchterman.reingold get.edgelist
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.plotStateTrajectory_def <- function(g, y, all_sts, name, setND=FALSE,
                                     point_size=3, label_offset=2, seed=1101) {
  # Vertex names
  component_sts <- names(V(g))
  o <- order(as.numeric(substring(component_sts, 2)))
  vnames <- factor(component_sts, levels = component_sts[o])

  # Layout
  set.seed(seed)
  #Y <- data.frame(layout.fruchterman.reingold(g))
  X <- layout.fruchterman.reingold(g)
  X <- matrix(apply(X, 2L, .rescale, ymin = 0, ymax = 1), ncol = 2)
  X <- data.frame(X)

  colnames(X) <- c("X1", "X2")
  rownames(X) <- vnames
  edgelist <- get.edgelist(g)
  edges <- data.frame(X[edgelist[,1],], X[edgelist[,2],])
  colnames(edges) <- c("X1", "Y1", "X2", "Y2")
  gp <- ggplot() +
    geom_segment(aes_string(x="X1", y="Y1", xend="X2", yend="Y2"),
                 data=edges, size=0.5, colour="gray")
  lblfactor <- 1

  # Graph visualization
  y <- y[all_sts %in% vnames]
  all_sts <- all_sts[all_sts %in% vnames]
  if(is.numeric(y)) {
    dat <- aggregate(y, list(all_sts), mean)
    dat <- dat[match(vnames, dat[, 1]), 2]

    nas <- rep(FALSE, length(dat))
    limits <- range(dat, na.rm=TRUE)
    #lbls <- pretty(limits)

    if(setND) { #is expression data
      nas <- dat == 0
      limits[1] <- min(dat[!nas])
      #lbls <- gsub(lbls, pattern = "^0", replacement = "nd")
    }
    gp <- gp + aes(x=X[nas, 1], y=X[nas, 2]) +
      geom_point(pch=1, col="gray40", size=point_size) +
      geom_point(aes(X[!nas, 1], X[!nas, 2],
                     color=dat[!nas]), size=point_size) +
      scale_color_gradientn(colours=.prettyColorRamp(5, grayStart=FALSE),
                            #label=lbls,
                            limits=limits, values=seq(0,1,by=0.2)) +
      labs(colour=name)
  } else {
    y <- factor(y)
    df <- data.frame(X[as.character(all_sts), ], STATE=all_sts, X3=y,
                     K=factor(1))
    l <- list()
    for(i in vnames) {
      l[[i]] <- ggplot(df[df$STATE == i, ], aes_string(x="K",
                                                       fill="X3")) +
        geom_bar(width = 1, color = "black") +
        scale_fill_discrete(drop=FALSE) +
        scale_x_discrete(drop=FALSE) +
        coord_polar("y") +
        #geom_col(position = "fill", alpha = 0.75, colour = "white") +
        #coord_polar("y", start = 1) +  #theta =
        guides(fill=FALSE) +
        theme_void()
    }

    l2 <- list()
    #d1 <- abs(diff(range(Y[,1]))) / (nrow(Y))
    #d2 <- abs(diff(range(Y[,2]))) / (nrow(Y))
    #d <- min(d1, d2)
    d <- 0.05 * point_size

    for(i in seq_along(l)) {
      ac <- annotation_custom(grob=ggplotGrob(l[[i]]),
                              xmin=X[i, 1] - d, xmax=X[i, 1] + d,
                              ymin=X[i, 2] - d, ymax=X[i, 2] + d)
      l2[[i]] <- ac
    }
    gp <- ggplot() + geom_segment(data=edges, aes_string(x="X1",
                                                         y="Y1",
                                                         xend="X2",
                                                         yend="Y2"),
                                  size=0.5, colour="gray")
    #ggplot(data = Y, aes(X1, X2)) + l2 + geom_point()
    df.f <- df[1, ]
    gp <- gp + geom_point(data=df, aes_string(x="X1",
                                              y="X2",
                                              colour="X3")) + l2
    gp <- gp + guides(col=guide_legend(override.aes=list(shape=15, size=5)))
    #gp <- gp + geom_col(data = df, aes(x=0, y=0, fill = X3)) + l2
    gp <- gp + labs(colour=name)
  }

  #lbls
  gp.xrange <- ggplot_build(gp)$layout$panel_ranges[[1]]$x.range
  gp.yrange <- ggplot_build(gp)$layout$panel_ranges[[1]]$y.range
  rng_x <- abs(gp.xrange[1] - gp.xrange[2])
  rng_y <- abs(gp.yrange[1] - gp.yrange[2])

  gp <- gp + geom_text(aes_string(x="X1", y="X2", label="vnames"),
                       hjust=0, vjust=0, data=X,
                       nudge_x= rng_x * 0.01 * label_offset,
                       nudge_y = rng_y * 0.01 * label_offset) +
    xlim(gp.xrange[1] - rng_x * 0.01 * label_offset,
                  gp.xrange[2] + rng_x * 0.01 * label_offset) +
    ylim(gp.yrange[1] - rng_y * 0.01 * label_offset,
                  gp.yrange[2] + rng_y * 0.01 * label_offset) +
    theme(axis.line=element_line(colour="black"),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank()) +
    xlab(paste0("State trajectory graph 1")) +
    ylab(paste0("State trajectory graph 2"))
  gp
}

#' DEF: Visualizing the trajectory fit
#'
#' Shows the residuals along the trajectory backbone
#' @param x The fitting residuals
#' @param g The trajectory graph
#' @param sts The states along the trajectory
#' @param factor The jitter intensity correlates to the projection error
#' @param rev Should the trajectory shown in reverse order?
#' @return A \code{ggplot} object
#' @importFrom igraph get.edgelist
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.plot_trajectoryFit <- function(x, g, sts, factor=7, rev=FALSE) {
  Y <- .generate_ordination(g, x, factor=factor, rev=rev)
  Y$ordi <- Y$ordi * 100
  Y$ordi.jitter <- Y$ordi.jitter * 100
  edgelist <- get.edgelist(g)
  edges <- data.frame(Y$ordi[edgelist[,1],], Y$ordi[edgelist[,2],],
                      row.names=seq(nrow(edgelist)))
  colnames(edges) <- c("X1", "Y1", "X2", "Y2")

  nas <- is.na(sts)
  gp <- ggplot() +
    aes(x = Y$ordi.jitter[nas, 1], y = Y$ordi.jitter[nas, 2]) +
    geom_point(pch=1, col="gray40") +
    labs(colour="State") +
    geom_point(mapping=aes(x=Y$ordi.jitter[!nas, 1],
                           y=Y$ordi.jitter[!nas, 2],
                           color=factor(sts[!nas]))) +
    geom_segment(aes_string(x="X1", y="Y1", xend="X2", yend="Y2"),
                 data=edges, size = 0.5, colour="black") +
    xlab("Pseudotime (%)") +
    ylab("Pseudotime (%)") +
    theme(axis.line=element_line(colour="black"))
  gp
}

#' DEF: Visualizing trailblazing
#'
#' Shows trajectory graph and highlights landmarks
#' @param X The trajectory layout
#' @param g The trajectory graph
#' @param ltype Vector with landmark types
#' @param lid Vector with landmark ids
#' @return A \code{ggplot} object
#' @importFrom igraph V get.edgelist
#' @importFrom ggrepel geom_label_repel
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.plotTrailblazing_def <- function(X, g, ltype, lid) {
  sname <- V(g)$sampleName
  edgelist <- get.edgelist(g)
  edgelist <- data.frame(X1=V(g)[edgelist[,1]]$sampleName,
                         X2=V(g)[edgelist[,2]]$sampleName,
                         stringsAsFactors=FALSE)
  edges <- data.frame(X[edgelist[,1], ],
                      X[edgelist[,2], ],
                      row.names = seq(nrow(edgelist)))
  colnames(edges) <- c("X1", "Y1", "X2", "Y2")

  nnas <- !is.na(ltype)
  tlayout <- data.frame(D1=X[,1], D2=X[,2],
                        label=lid,
                        Type=ltype)
  tlayout <- tlayout[nnas, ]

  gp <- ggplot()
  gp <- gp + geom_segment(data=edges,
                          aes_string(x="X1", y="Y1", xend="X2", yend="Y2"),
                          size=0.5, colour="darkred") +
    geom_label_repel(data=tlayout,
                     mapping=aes_string(x="D1", y="D2",
                                        label="label", fill="Type"),
                     fontface='bold', color='white',
                     #box.padding = unit(0.35, "lines"),
                     #point.padding = unit(0.5, "lines"),
                     #col = "white", size = 2,
                     segment.color="gray40", size=2,
                     min.segment.length=unit(0, 'lines')) +
    geom_point(data=tlayout,
               mapping=aes_string(x="D1", y="D2", color="Type")) +
    theme(axis.line=element_line(colour="black"),
          panel.border=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    xlab('CellTrails 1') + ylab('CellTrails 2')
  #gp <- gp + scale_x_continuous(expand = c(0,0)) +
  #scale_y_continuous(expand = c(0,0))
  gp
}

#' DEF: Visualize trail on map
#'
#' Highlights a trail on map
#' @param X The trajectory layout
#' @param g The trajectory graph
#' @param ptime The trail pseudotime
#' @param name Name of trail
#' @keywords internal
#' @importFrom igraph V get.edgelist induced_subgraph
#' @import ggplot2
#' @author Daniel C. Ellwanger
.plotTrail_def <- function(X, g, ptime, name) {
  sname <- V(g)$sampleName
  edgelist <- get.edgelist(g)
  edgelist <- data.frame(X1=V(g)[edgelist[,1]]$sampleName,
                         X2=V(g)[edgelist[,2]]$sampleName,
                         stringsAsFactors=FALSE)
  edges <- data.frame(X[edgelist[,1], ],
                      X[edgelist[,2], ],
                      row.names=seq(nrow(edgelist)))
  colnames(edges) <- c("X1", "Y1", "X2", "Y2")

  gp <- ggplot() +
    geom_segment(aes_string(x="X1", y="Y1", xend="X2", yend="Y2"),
                 data=edges, size=0.5, colour="gray40")

  X$ptime <- ptime[,1]
  pth <- rownames(ptime)[!is.na(ptime)[,1]]
  ptime <- ptime[,1] / max(ptime[,1], na.rm=TRUE) * 100
  ptime <- ptime[!is.na(ptime)]

  g <- induced_subgraph(g, V(g)[V(g)$sampleName %in% pth])
  sname <- V(g)$sampleName
  edgelist <- get.edgelist(g)
  edgelist <- data.frame(X1=V(g)[edgelist[,1]]$sampleName,
                         X2=V(g)[edgelist[,2]]$sampleName,
                         stringsAsFactors=FALSE)
  edges <- data.frame(X[edgelist[,1], seq_len(2)],
                      X[edgelist[,2], seq_len(2)],
                      row.names = seq(nrow(edgelist)))

  #edgelist <- get.edgelist(g)
  #edges <- data.frame(tlayout[edgelist[,1],], tlayout[edgelist[,2],],
  #                    row.names = seq(nrow(edgelist)))
  colnames(edges) <- c("U1", "V1", "U2", "V2")
  startEnd <- data.frame(X[pth[order(ptime)[c(1, length(ptime))]], ],
                         label = c("Start", "End"))

  gp <- gp + geom_segment(aes_string(x="U1", y="V1", xend="U2", yend="V2"),
                          data=edges, size=0.5, color="black") +
    geom_point(data=X[pth, ],
               aes_string(x="D1", y="D2", color=ptime), size=0.75) +
    geom_label(data=startEnd, mapping=aes_string(x="D1", y="D2", label="label"),
               fontface='bold',
               #box.padding = unit(0.35, "lines"),
               #point.padding = unit(0.5, "lines"),
               size=2) + #col = "white", size = 2,
    theme(axis.line=element_line(colour="black"),
          panel.border=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    xlab('CellTrails 1') + ylab('CellTrails 2') +
    labs(colour = "Time (%)") +
    labs(subtitle = paste0("Trail: ", name))
  gp
}

#' DEF: Visualize expression dynamic
#'
#' Shows feature expression as a function of pseudotime
#' @param x Pseudotime
#' @param Y Expression matrix (samples in columns)
#' @param weights Weight (states)
#' @param k Dimensions of bases used to represent the smooth term
#' @keywords internal
#' @import ggplot2
#' @importFrom reshape2 melt
#' @author Daniel C. Ellwanger
.plotDynamic <- function(x, Y, feature_name, trail_name, weights, k=5) {
  gp <- ggplot()

  if(nrow(Y) == 1) { #show one feature
    dat <- data.frame(X=x[!is.na(x)] / max(x, na.rm=TRUE) * 100,
                      Y=Y[1, !is.na(x)],
                      STATES=weights[!is.na(x)])
    fit <- .fitDynamic_def(x = dat[,1], y = dat[,2], z = dat[,3], k = k)
    fit.dat <- data.frame(FX = fit$x / max(fit$x) * 100, FY = fit$y)

    gp <- gp + geom_point(data = dat, aes_string(x="X", y="Y", color="STATES"))
    gp <- gp + labs(colour = "State")
    gp <- gp + geom_line(data = fit.dat, aes_string(x="FX", y="FY"), lwd = .75)
    gp <- gp + theme(axis.line = element_line(colour = "black"))
    gp <- gp + xlab('Pseudotime (%)') + ylab(feature_name)

    brks <- pretty(dat$Y)
    lbl <- brks
    lbl[lbl == 0] <- "nd"
    gp <- gp + scale_y_continuous(breaks = brks, labels = lbl)
  } else { #show multiple features
    ys <- matrix(nrow = length(x[!is.na(x)]), ncol = length(feature_name))
    for(i in seq_along(feature_name)) {
      ys[, i] <- .fitDynamic_def(x = x[!is.na(x)],
                              y = Y[i, !is.na(x)],
                              z = weights[!is.na(x)],
                              k = k)$y
    }
    colnames(ys) <- feature_name
    dat <- data.frame(ptime=sort(x[!is.na(x)] / max(x, na.rm=TRUE) * 100), ys)
    dat <- melt(dat, id = "ptime", variable.name = "Feature")

    gp <- ggplot(data=dat,
                 aes_string(x="ptime", y="value", colour="Feature"))
    gp <- gp + geom_line(lwd = .75)
    gp <- gp + theme(axis.line = element_line(colour = "black"))
    gp <- gp + xlab('Pseudotime (%)') + ylab("Expression")

    brks <- pretty(dat$value)
    lbl <- brks
    lbl[lbl == 0] <- "nd"
    gp <- gp + scale_y_continuous(breaks=brks, labels=lbl)
  }
  gp <- gp + labs(subtitle = paste0("Trail: ", trail_name))
  gp
}
