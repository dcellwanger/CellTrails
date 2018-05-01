#' DEF: Visualizing the spectrum definition
#'
#' For details see \code{plot.CellTrailsSpectrum}
#' @param x A \code{CellTrailsSpectrum} object
#' @return A \code{ggplot} object
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.plot_spectrum <- function(x) {
  if(x@frac <= 1) {
    x@frac <- x@frac * length(x@cs)
  }
  Y <- x@cs[seq(x@frac)]
  gp <- ggplot()
  gp <- gp + aes(x=seq(x@frac), y=Y) +
    geom_point(color=c(rep("red", x@n - 1), rep("gray40", x@frac - x@n + 1)))
  gp <- gp + geom_line(mapping=aes(x=x@fit$x, y=x@fit$y, lty = 'Fit'),
                       color = "blue")
  gp <- gp + xlab("Spectrum")
  gp <- gp + ylab("Total eigengap")
  brks <- sort(c(pretty(seq(x@frac)), x@n - 1))
  gp <- gp + scale_x_continuous(breaks=brks, labels=brks + 1)
  col <- rep("black", length(brks))
  col[brks == x@n - 1] <- "red"
  gp <- gp + theme(axis.line=element_line(colour = "black"),
                   axis.text.x=element_text(color = col),
                   axis.ticks.x=element_line(color = col),
                   legend.title=element_blank())
  #gp <- gp + geom_vline(xintercept = 8)
  gp
}

#' DEF: Visualizing states sizes
#'
#' For details see \code{plot.CellTrailsSet}
#' @param x A \code{CellTrailsSet} object
#' @return A \code{ggplot} object
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.plot_stateSize <- function(x) {
  count <- data.frame(table(states(x)))
  ggplot(count, aes_string(x = "Var1", y = "Freq", fill = "Var1")) +
    geom_bar(stat = "identity") +
    labs(x = "State", y = "Sample count") +
    guides(fill = FALSE) +
    theme(axis.line = element_line(colour = "black"))
}

#' DEF: Violine plots of genes
#'
#' For details see \code{plot.CellTrailsSet}
#' @param x A \code{CellTrailsSet} object
#' @param feature_name Name of feature
#' @return A \code{ggplot} object
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.plot_stateExpression <- function(x, feature_name) {
  .checkFeatureNameExists(x, feature_name)

  df <- data.frame(State=x$STATE, Expression=x[feature_name, ])
  brks <- pretty(df$Expression)
  lbl <- brks
  lbl[lbl == 0] <- "nd"

  ggplot(df, aes_string(x = "State", y = "Expression", fill = "State")) +
    geom_jitter(shape = 16, position = position_jitter(0.2), alpha = .5) +
    geom_violin(trim = TRUE, alpha = .8) +
    ylab(feature_name) +
    scale_y_continuous(breaks = brks, labels = lbl) +
    theme(axis.line = element_line(colour = "black"))
}

#' DEF: Visualizing the lower-dimensional manifold
#'
#' For details see \code{plot.CellTrailsSet}
#' @param x A \code{CellTrailsSet} object
#' @param pheno_type Type of sample meta information; must match entity
#' in \code{varLabels}.
#' @param feature_name Name of feature; must match entity in \code{featureNames}.
#' @param seed Seed for tSNE computation
#' @param perplexity Perplexity parameter of tSNE calculation
#' @param viz A \code{ggplot} object containing the result of a previous
#' \code{plot_latentSpace} call
#' @return A \code{ggplot} object
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.plot_latentSpace <- function(x, pheno_type=NULL,
                            feature_name=NULL,
                            seed=1101, perplexity=30, viz=NULL, ...) {
  #Pre-flight checks
  if(is.null(pheno_type) && is.null(feature_name)) {
    stop("Wrong parameter setting: define either a pheno_type or a ",
         "feature_name, please.")
  }
  if(!is.null(pheno_type) && !is.null(feature_name)) {
    stop("Wrong parameter setting: define either a pheno_type or a ",
         "feature_name, please.")
  }

  if(is.null(viz)) {
    message("Computing 2D visualization of manifold ...")
    Y <- .bhtsne(x, seed=seed, perplexity=perplexity)$Y

  } else {
    Y <- viz$tsne
    #Y <- ggplot_build(viz)$data[[2]][, c("x", "y")]
  }

  dat <- NULL
  gp <- ggplot()
  gp$tsne <- Y

  if(!is.null(pheno_type)) {
    pheno_type <- ifelse(toupper(pheno_type) == "STATE", "STATE", pheno_type)
    if(!pheno_type %in% varLabels(x)) {
      stop("Variable label '", pheno_type, "' not found. ",
           "Please, check correct spelling. ",
           "To see all available pheno_types call 'varLabels(object)'.")
    }
    dat <- x[[pheno_type]]
    if(is.numeric(dat)) {
      nas <- dat == 0
      #dat[nas] <- NA
      gp <- gp + aes(x = Y[nas, 1], y = Y[nas, 2])
      gp <- gp + geom_point(pch = 1, col = "gray40")
      gp <- gp + xlab("CellTrails tSNE 1")
      gp <- gp + ylab("CellTrails tSNE 2")
      gp <- gp + labs(colour = pheno_type)

      gp <- gp + geom_point(mapping=aes(x=Y[!nas,1], y=Y[!nas,2],
                                          color=dat[!nas]))
      gp <- gp + scale_color_gradientn(colours=viridis(3),
                                       limits=c(1e-10, max(dat)))
    } else {
      nas <- is.na(dat)
      gp <- gp + aes(x = Y[nas, 1], y = Y[nas, 2])
      gp <- gp + geom_point(pch = 1, col = "gray40")
      gp <- gp + xlab("CellTrails tSNE 1")
      gp <- gp + ylab("CellTrails tSNE 2")
      gp <- gp + labs(colour = pheno_type)
      gp <- gp + geom_point(mapping=aes(x=Y[!nas, 1], y=Y[!nas, 2],
                                        color=factor(dat[!nas])))
    }
  } else if(!is.null(feature_name)) {
    .checkFeatureNameExists(x, feature_name)

    dat <- x[feature_name, ]
    nas <- dat == 0
    #dat[nas] <- NA
    gp <- gp + aes(x = Y[nas, 1], y = Y[nas, 2])
    gp <- gp + geom_point(pch = 1, col = "gray40")
    gp <- gp + xlab("CellTrails tSNE 1")
    gp <- gp + ylab("CellTrails tSNE 2")
    gp <- gp + labs(colour = feature_name)

    gp <- gp + geom_point(mapping=aes(x=Y[!nas,1],
                                      y=Y[!nas,2],
                                      color=dat[!nas]))
    gp <- gp + scale_color_gradientn(colours=viridis(3),
                                     limits=c(1e-10, max(dat)))
  }
  gp <- gp + theme(axis.line=element_line(colour="black"),
                   panel.border=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())
  gp
}

#' DEF: Visualizing the trajectory graph
#'
#' For details see \code{plot.CellTrailsSet}
#' @param x A \code{CellTrailsSet} object
#' @param feature_name Name of the feature (optional)
#' @param seed Makes result reproducible
#' @param component Component of the trajectory graph
#' @return A \code{ggplot} object
#' @importFrom igraph V layout.fruchterman.reingold get.edgelist
#' @import Biobase
#' @importFrom viridis viridis
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.plot_stateTrajectoryGraph <- function(x, feature_name=NULL,
                                       pheno_type=NULL, component=NULL,
                                       point_size=3, label_offset=2,
                                       seed=1101) {
  #Pre-flight checks
  if(is.null(component) && length(stateTrajectoryGraph(x)) > 1) {
    stop("Please, specify the component.")
  }
  if(is.null(component)) {
    component <- 1
  }
  if(component > length(stateTrajectoryGraph(x))) {
    stop("Unknown component selected. Please, make sure that the ",
         "correct component number was selected.")
  }
  if(is.null(pheno_type) && is.null(feature_name)) {
    stop("Wrong parameter setting: define either a pheno_type or a ",
         "feature_name, please.")
  }
  if(!is.null(pheno_type) && !is.null(feature_name)) {
    stop("Wrong parameter setting: define either a pheno_type or a ",
         "feature_name, please.")
  }

  g <- stateTrajectoryGraph(x)[[component]]
  sts <- names(V(g))
  o <- order(as.numeric(substring(sts, 2)))
  vnames <- factor(sts, levels = sts[o])

  set.seed(seed)
  #Y <- data.frame(layout.fruchterman.reingold(g))
  Y <- layout.fruchterman.reingold(g)
  Y <- matrix(apply(Y, 2L, .rescale, ymin = 0, ymax = 1), ncol = 2)
  Y <- data.frame(Y)

  colnames(Y) <- c("X1", "X2")
  rownames(Y) <- vnames
  edgelist <- get.edgelist(g)
  edges <- data.frame(Y[edgelist[,1],], Y[edgelist[,2],])
  colnames(edges) <- c("X1", "Y1", "X2", "Y2")
  gp <- ggplot() + geom_segment(aes_string(x="X1",
                                           y="Y1",
                                           xend="X2",
                                           yend="Y2"),
                                data=edges,
                                size = 0.5,
                                colour="gray")
  lblfactor <- 1

  if(!is.null(pheno_type)) {
    pheno_type <- ifelse(toupper(pheno_type) == "STATE", "STATE", pheno_type)
    if(!pheno_type %in% varLabels(x)) {
      stop("Variable label '", pheno_type, "' not found. Please, check correct ",
           "spelling. ", "To see all available pheno_types call ",
           "'varLabels(object)'.")
    }
    pdata <- pData(x)[, c("STATE", pheno_type)]

    impute <- TRUE #set always true
    if(impute) { #perform imputation of missing labels
      nas <- which(is.na(pdata[,2]))
      nnas <- which(!is.na(pdata[,2]))
      D <- as.matrix(dist(CellTrails::latentSpace(x)))[nas, nnas]
      nn <- apply(D, 1L, which.min)
      pdata[nas, 2] <- pdata[nnas[nn], 2]
    }
    pdata <- pdata[pdata$STATE %in% vnames, ]

    if(is.numeric(pdata[, 2])) {
      dat <- aggregate(pdata[, 2], list(pdata[, 1]), mean)
      tmax <- max(dat[, 2])
      dat <- dat[match(vnames, dat[, 1]), 2]
      nas <- dat == 0
      gp <- gp + aes(x = Y[nas, 1], y = Y[nas, 2])
      gp <- gp + geom_point(pch = 1, col = "gray40", size = point_size)
      gp <- gp + geom_point(aes(Y[!nas, 1], Y[!nas, 2], color=dat[!nas]),
                            size=point_size)
      gp <- gp + scale_color_gradientn(colours = viridis(3),
                                       limits=c(1e-10, max(tmax)))
      gp <- gp + labs(colour = pheno_type)
    } else {
      pdata[, 2] <- factor(pdata[, 2])

#      if(length(levels(pdata[,2])) == length(levels(pdata[,1]))) {
#        #plot(Y, col = rainbow(12)[vnames])
#        gp <- gp + geom_point(aes(X1, X2, color = vnames), data=Y, size = point_size)
#        gp <- gp + labs(colour = pheno_type)
#      } else {
          df <- data.frame(Y[as.character(pdata$STATE), ], STATE=pdata$STATE,
                           X3=pdata[,2], K=factor(1))

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
                                    xmin=Y[i, 1] - d, xmax=Y[i, 1] + d,
                                    ymin=Y[i, 2] - d, ymax=Y[i, 2] + d)
            l2[[i]] <- ac
          }
          gp <- ggplot() + geom_segment(data=edges, aes_string(x="X1",
                                                               y="Y1",
                                                               xend="X2",
                                                               yend="Y2"),
                                        size = 0.5, colour="gray")
          #ggplot(data = Y, aes(X1, X2)) + l2 + geom_point()
          df.f <- df[1, ]
          gp <- gp + geom_point(data = df, aes_string(x="X1",
                                                      y="X2",
                                                      colour="X3")) + l2
          gp <- gp + guides(col = guide_legend(override.aes = list(shape=15,
                                                                   size=5)))
          #gp <- gp + geom_col(data = df, aes(x=0, y=0, fill = X3)) + l2
          gp <- gp + labs(colour = pheno_type)
#      }
    }
  } else if(!is.null(feature_name)) {
    .checkFeatureNameExists(x, feature_name)

    dat <- aggregate(x[feature_name, ], list(states(x)), mean)
    tmax <- max(dat[,2])
    dat <- dat[match(vnames, dat[,1]), 2]
    nas <- dat == 0
    gp <- gp + aes(x = Y[nas, 1], y = Y[nas, 2])
    gp <- gp + geom_point(pch = 1, col = "gray40", size = point_size)
    gp <- gp + geom_point(aes(Y[!nas, 1],
                              Y[!nas, 2],
                              color=dat[!nas]), size=point_size)
    gp <- gp + scale_color_gradientn(colours=viridis(3),
                                     limits=c(1e-10, max(tmax)))
    gp <- gp + labs(colour = feature_name)
  }
  #lbls
  gp.xrange <- ggplot_build(gp)$layout$panel_ranges[[1]]$x.range
  gp.yrange <- ggplot_build(gp)$layout$panel_ranges[[1]]$y.range
  gp <- gp + geom_text(aes_string(x="X1", y="X2", label="vnames"),
                       hjust=0, vjust=0, data=Y,
                       nudge_x=abs(gp.xrange[1] - gp.xrange[2]) * 0.01 * label_offset,
                       nudge_y = abs(gp.yrange[1] - gp.yrange[2]) * 0.01 * label_offset)
  #gp <- gp + geom_text_repel(aes(X1, X2, label=vnames), data=Y, segment.size = 0)
  gp <- gp + xlim(gp.xrange[1] - abs(gp.xrange[1] - gp.xrange[2]) * 0.01 * label_offset,
                  gp.xrange[2] + abs(gp.xrange[1] - gp.xrange[2]) * 0.01 * label_offset)
  gp <- gp + ylim(gp.yrange[1] - abs(gp.yrange[1] - gp.yrange[2]) * 0.01 * label_offset,
                  gp.yrange[2] + abs(gp.yrange[1] - gp.yrange[2]) * 0.01 * label_offset)
  gp <- gp + theme(axis.line=element_line(colour="black"),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())

  gp <- gp + xlab(paste0("State trajectory graph 1"))
  gp <- gp + ylab(paste0("State trajectory graph 2"))
  gp
}

#' DEF: Visualizing the trajectory fit
#'
#' For details see \code{plot.CellTrailsSet}
#' @param x A \code{CellTrailsSet} object
#' @param factor The jitter intensity correlates to the projection error
#' @param rev Should the trajectory shown in reverse order?
#' @return A \code{ggplot} object
#' @importFrom igraph get.edgelist
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.plot_trajectoryFit <- function(x, factor=7, rev=FALSE) {
  Y <- .generate_ordination(x@trajectory$traj, x@trajectory$error,
                            factor=factor, rev=rev)
  Y$ordi <- Y$ordi * 100
  Y$ordi.jitter <- Y$ordi.jitter * 100
  edgelist <- get.edgelist(x@trajectory$traj)
  #edgelist <- apply(edgelist, 2, as.numeric)
  #rownames(Y$ordi) <- NULL
  edges <- data.frame(Y$ordi[edgelist[,1],], Y$ordi[edgelist[,2],],
                      row.names=seq(nrow(edgelist)))
  colnames(edges) <- c("X1", "Y1", "X2", "Y2")
  gp <- ggplot()

  dat <- states(x)[x@useSample]
  nas <- is.na(dat)
  gp <- gp + aes(x = Y$ordi.jitter[nas, 1], y = Y$ordi.jitter[nas, 2])
  gp <- gp + geom_point(pch = 1, col = "gray40")
  gp <- gp + labs(colour = "State")
  gp <- gp + geom_point(mapping = aes(x=Y$ordi.jitter[!nas, 1],
                                      y=Y$ordi.jitter[!nas, 2],
                                      color=factor(dat[!nas])))
  gp <- gp + geom_segment(aes_string(x="X1", y="Y1", xend="X2", yend="Y2"),
                                      data=edges, size = 0.5, colour="black")
  gp <- gp + xlab("Pseudotime (%)")
  gp <- gp + ylab("Pseudotime (%)")
  gp <- gp + theme(axis.line=element_line(colour="black"))
  gp
}

#' DEF: Visualizing CellTrails maps
#'
#' For details see \code{plot.CellTrailsSet}
#' @param x A \code{CellTrailsSet} object
#' @param feature_name Name of feature; must match entity in
#' \code{featureNames}.
#' @param pheno_type Type of sample meta information; must match entity
#' in \code{varLabels}.
#' @param npoints Points of the grids (along the x- and y-axis)
#' @param w Weights for GAM
#' @param knots Number of knots for GAM
#' @return A \code{ggplot} object
#' @importFrom igraph V get.edgelist
#' @import Biobase
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.plot_map <- function(x, feature_name=NULL, pheno_type=NULL, npoints=300,
                      weight=TRUE, knots=10, only_backbone=FALSE,
                      map_type=c("full", "backbone", "se")) {
  #Pre-flight check
  if(is.null(pheno_type) && is.null(feature_name)) {
    stop("Wrong parameter setting: define either a pheno_type or a ",
         "feature_name, please.")
  }
  if(!is.null(pheno_type) && !is.null(feature_name)) {
    stop("Wrong parameter setting: define either a pheno_type or a ",
         "feature_name, please.")
  }
  if(is.null(trajectoryLayout(x))) {
    stop("No layout defined. Please, set a layout to your CellTrailsSet object
         first (use 'trajectoryLayout(object) <- value').")
  }
  map_type <- toupper(map_type[1])

  tlayout <- trajectoryLayout(x)
  g <- x@trajectory$traj
  sname <- V(g)$sampleName
  edgelist <- get.edgelist(g)
  edgelist <- data.frame(X1=V(g)[edgelist[,1]]$sampleName,
                         X2=V(g)[edgelist[,2]]$sampleName,
                         stringsAsFactors=FALSE)
  edges <- data.frame(tlayout[edgelist[,1], ],
                      tlayout[edgelist[,2], ],
                      row.names = seq(nrow(edgelist)))
  colnames(edges) <- c("X1", "Y1", "X2", "Y2")

  gp <- ggplot()

  if(!is.null(feature_name)){
    .checkFeatureNameExists(x, feature_name)

    if(map_type == "SINGLE") {
      #Simplified fit on small grid
      fit <- .fit_surface(x, feature_name=feature_name, npoints=1,
                          weight=weight, knots=knots)
      fval <- fit$fit$fitted.values
      brks <- pretty(range(fval, na.rm=TRUE, finite=TRUE), 10)
      #lbls <- paste(c("nd", brks[2:(length(brks) - 1)]), brks[2:(length(brks))], sep = " - ")
      lbls <- paste(brks[seq_len(length(brks)-1)],
                    brks[seq_len(length(brks)-1)+1], sep = " - ")
      lbls <- gsub(pattern="^0 -", replacement = "nd -", x = lbls)

      equalSpace <- cut(fval, breaks=brks, include.lowest=TRUE, labels=lbls)
      breaks <- rev(levels(equalSpace))
      #nbrks <- length(brks)
      #cols <- colorRampPalette(c("gray95", viridis(3)[2:3]))(nbrks - 1)
      #cols <- cols[cut(fval, breaks = brks, include.lowest = T)]
      #cols <- adjustcolor(cols, alpha.f = .5)

      gp <- ggplot() + geom_segment(aes_string(x="X1",
                                               y="Y1",
                                               xend="X2",
                                               yend="Y2"),
                                    data=edges, size = 0.5, colour = "darkred")
      gp <- gp + geom_point(data=tlayout, mapping=aes_string(x="D1", y="D2",
                                                             col="equalSpace"))
      gp <- gp + scale_colour_manual(values=
                                      colorRampPalette(
                                        c("gray95", viridis(3)[2:3]))(length(brks) - 1),
                                      name=feature_name, breaks=breaks, labels=breaks)
      gp <- gp + theme_bw()
      gp <- gp + theme(axis.line=element_line(colour="black"),
                       panel.border=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())
    } else if(map_type == "SE") {
      # Fit on whole grid
      fit <- .fit_surface(x, feature_name=feature_name,
                          npoints=npoints, weight=weight, knots=knots)
      df <- fit$grid
      df$x3 <- fit$fit$se
      brks <- pretty(range(df$x3, na.rm=TRUE, finite=TRUE), 10)
      #lbls <- paste(c("nd", brks[2:(length(brks) - 1)]), brks[2:(length(brks))], sep = " - ")
      lbls <- paste(brks[seq_len(length(brks)-1)],
                    brks[seq_len(length(brks)-1)+1], sep = " - ")
      #lbls <- gsub(pattern="^0 -", replacement = "nd -", x = lbls)
      df$equalSpace <- cut(df$x3, breaks=brks, include.lowest=TRUE, labels=lbls)
      breaks <- rev(levels(df$equalSpace))

      #df <- data.frame(expand.grid(fit$grid$x, fit$grid$y), matrix(fit$grid$z, ncol = 1))
      #colnames(df) <- c("x1", "x2", "x3")
      #library(directlabels)
      #lbls <- paste0(rep("(", length(brks) - 1), lbls, rep("]", length(brks) - 1))
      breaks <- rev(levels(df$equalSpace))

      gp <- ggplot() + geom_tile(data=df, aes_string(x="x1",
                                                     y="x2",
                                                     fill="equalSpace"))
      gp <- gp + geom_contour(data=df, aes_string(x="x1", y="x2", z="x3"),
                              color="gray40", alpha=0.5, lty=2, lwd=.5)
      gp <- gp + scale_fill_manual(values =
                                     colorRampPalette(
                                       c("white", "black"))(length(brks) - 1),
                                   name=paste0(feature_name, "(SE)"),
                                   breaks=breaks,
                                   labels=breaks)
      gp <- gp + geom_segment(aes_string(x="X1", y="Y1", xend="X2", yend="Y2"),
                              data=edges,
                              size=0.5,
                              colour="darkred")
      gp <- gp + theme(axis.line=element_line(colour="black"),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       panel.border=element_blank(),
                       panel.background=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())
      gp <- gp + geom_point(data=tlayout, aes_string(x="D1", y="D2"),
                            colour="red", size=1.5, alpha=.1)
    } else if (map_type == "FULL") {
      # Fit on whole grid
      fit <- .fit_surface(x, feature_name=feature_name,
                          npoints=npoints, weight=weight, knots=knots)
      df <- fit$grid
      brks <- pretty(range(df$x3, na.rm=TRUE, finite=TRUE), 10)
      lbls <- paste(c("nd", brks[2:(length(brks) - 1)]), brks[2:(length(brks))],
                    sep = " - ")
      lbls <- paste(brks[seq_len(length(brks)-1)],
                    brks[seq_len(length(brks)-1)+1], sep = " - ")
      lbls <- gsub(pattern="^0 -", replacement = "nd -", x = lbls)
      df$equalSpace <- cut(df$x3, breaks=brks, include.lowest=TRUE, labels=lbls)
      breaks <- rev(levels(df$equalSpace))

      #df <- data.frame(expand.grid(fit$grid$x, fit$grid$y), matrix(fit$grid$z, ncol = 1))
      #colnames(df) <- c("x1", "x2", "x3")
      #library(directlabels)
      #lbls <- paste0(rep("(", length(brks) - 1), lbls, rep("]", length(brks) - 1))
      breaks <- rev(levels(df$equalSpace))

      gp <- ggplot() + geom_tile(data=df, aes_string(x="x1", y="x2", fill="equalSpace"))
      gp <- gp + geom_contour(data=df, aes_string(x="x1", y="x2", z="x3"),
                              color="gray40", alpha=0.5, lty=2, lwd=.5)
      gp <- gp + scale_fill_manual(values =
                                     colorRampPalette(c("gray95", viridis(3)[2:3]))(length(brks) - 1),
                                     name=feature_name, breaks = breaks, labels = breaks)
      gp <- gp + geom_segment(aes_string(x="X1", y="Y1", xend="X2", yend="Y2"), data=edges,
                              size=0.5, colour="darkred")
      gp <- gp + theme(axis.line=element_line(colour = "black"),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       panel.border=element_blank(),
                       panel.background=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())
      gp <- gp + geom_point(data=tlayout, aes_string(x="D1", y="D2"),
                                          colour="red", size=1.5, alpha=.1)
    } else {
      stop("Unknown map type selected. Please, choose one from ",
           "{full, backbone, se}.")
    }
  } else if(!is.null(pheno_type)) {
    pheno_type <- ifelse(toupper(pheno_type) == "STATE", "STATE", pheno_type)
    if(!pheno_type %in% varLabels(x)) {
      stop("Variable label '", pheno_type, "' not found. Please, check ",
           "correct spelling. ",
           "To see all available pheno_types call 'varLabels(object)'.")
    }

    df <- pData(x)[x@useSample, pheno_type]
    nnas <- !is.na(df)
    df <- df[nnas]
    gp <- gp + geom_segment(aes_string(x="X1", y="Y1", xend="X2", yend="Y2"),
                            data=edges, size=0.5, colour="gray40")
    gp <- gp + labs(colour = pheno_type)
    gp <- gp + geom_point(data=tlayout[nnas, ],
                          mapping=aes_string(x="D1", y="D2", color=factor(df)))
    gp <- gp + theme(axis.line=element_line(colour="black"),
                     panel.border=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
  }
  gp <- gp + xlab('CellTrails 1') + ylab('CellTrails 2')
  gp <- gp + scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0))
  gp
}

#' DEF: Visualizing trailblazing
#'
#' For details see \code{plot.CellTrailsSet}
#' @param x A \code{CellTrailsSet} object
#' @return A \code{ggplot} object
#' @importFrom igraph V get.edgelist
#' @importFrom ggrepel geom_label_repel
#' @import ggplot2
#' @keywords internal
#' @author Daniel C. Ellwanger
.plot_trailblazing <- function(x) {
  tlayout <- trajectoryLayout(x)
  g <- x@trajectory$traj
  sname <- V(g)$sampleName
  edgelist <- get.edgelist(g)
  edgelist <- data.frame(X1=V(g)[edgelist[,1]]$sampleName,
                         X2=V(g)[edgelist[,2]]$sampleName,
                         stringsAsFactors=FALSE)
  edges <- data.frame(tlayout[edgelist[,1], ],
                      tlayout[edgelist[,2], ],
                      row.names = seq(nrow(edgelist)))
  colnames(edges) <- c("X1", "Y1", "X2", "Y2")

  nnas <- !is.na(x@trajectory$blaze$type)
  #Blaze <- x@trajectory$blaze$type[nnas]

  tlayout <- data.frame(D1=tlayout[,1], D2=tlayout[,2],
                        label=x@trajectory$blaze$id,
                        blaze=x@trajectory$blaze$type)
  tlayout <- tlayout[nnas, ]

  gp <- ggplot()
  gp <- gp + geom_segment(data=edges, aes_string(x="X1",
                                                 y="Y1",
                                                 xend="X2",
                                                 yend="Y2"),
                          size=0.5, colour="darkred")
  gp <- gp + geom_label_repel(data=tlayout,
                              mapping=aes_string(x="D1", y="D2",
                              label="label", fill="blaze"),
                              fontface='bold', color='white',
                              #box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines"),
                              segment.color="gray40", size=2,
                              min.segment.length=unit(0, 'lines')) #col = "white", size = 2,
  gp <- gp + geom_point(data=tlayout, mapping=aes_string(x="D1",
                                                         y="D2",
                                                         color="blaze"))
  gp <- gp + theme(axis.line=element_line(colour="black"),
                   panel.border=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())
  gp <- gp + xlab('CellTrails 1') + ylab('CellTrails 2')
  #gp <- gp + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
  gp
}

#' DEF: Visualize trail on map
#'
#' Highlights extracted trail on map
#' @param x A \code{CellTrailsTrail} object
#' @keywords internal
#' @importFrom igraph V get.edgelist induced.subgraph
#' @import ggplot2
#' @author Daniel C. Ellwanger
.plot_trail_on_map <- function(x, trail_name) {
  #Pre-flight check
  if(!trail_name %in% trailNames(x)) {
    stop("Trail with name '", trail_name, "' does not exists.")
  }

  tlayout <- trajectoryLayout(x)
  g <- x@trajectory$traj
  sname <- V(g)$sampleName
  edgelist <- get.edgelist(g)
  edgelist <- data.frame(X1=V(g)[edgelist[,1]]$sampleName,
                         X2=V(g)[edgelist[,2]]$sampleName,
                         stringsAsFactors=FALSE)
  edges <- data.frame(tlayout[edgelist[,1], ],
                      tlayout[edgelist[,2], ],
                      row.names=seq(nrow(edgelist)))
  colnames(edges) <- c("X1", "Y1", "X2", "Y2")

  nnas <- !is.na(x@trajectory$blaze$type)
  Blaze <- x@trajectory$blaze$type[nnas]

  gp <- ggplot()
  gp <- gp + geom_segment(aes_string(x="X1", y="Y1", xend="X2", yend="Y2"),
                          data=edges, size=0.5, colour="gray40")

  ptime <- x@trails[[trail_name]]$ptime
  ptime <- ptime / max(ptime) * 100
  pth <- x@trails[[trail_name]]$samples #x@trails[[trail_name]]$path

  g <- induced.subgraph(x@trajectory$traj, x@trails[[trail_name]]$path)
  sname <- V(g)$sampleName
  edgelist <- get.edgelist(g)
  edgelist <- data.frame(X1=V(g)[edgelist[,1]]$sampleName,
                         X2=V(g)[edgelist[,2]]$sampleName,
                         stringsAsFactors=FALSE)
  edges <- data.frame(tlayout[edgelist[,1], ],
                      tlayout[edgelist[,2], ],
                      row.names = seq(nrow(edgelist)))

  #edgelist <- get.edgelist(g)
  #edges <- data.frame(tlayout[edgelist[,1],], tlayout[edgelist[,2],], row.names = seq(nrow(edgelist)))
  colnames(edges) <- c("U1", "V1", "U2", "V2")
  startEnd <- data.frame(tlayout[sampleNames(x)[pth[c(1, length(pth))]], ],
                         label = c("Start", "End"))

  gp <- gp + geom_segment(aes_string(x="U1", y="V1", xend="U2", yend="V2"),
                          data=edges, size=0.5, color="black")
  gp <- gp + geom_point(data=tlayout[sampleNames(x)[pth], ],
                        aes_string(x="D1", y="D2", color=ptime), size=0.75)
  gp <- gp + geom_label(data=startEnd,
                              mapping=aes_string(x="D1",
                                                 y="D2",
                                                 label="label"),
                              fontface='bold', #box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines"),
                              size=2) #col = "white", size = 2,
  gp <- gp + theme(axis.line=element_line(colour="black"),
                   panel.border=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())
  gp <- gp + xlab('CellTrails 1') + ylab('CellTrails 2')
  gp <- gp + labs(colour = "Time (%)")
  gp <- gp + labs(subtitle = paste0("Trail: ", trail_name))
  gp
}

#' DEF: Visualize expression dynamic
#'
#' Shows feature expression as a function of pseudotime
#' @param x A \code{CellTrailsTrail} object
#' @param feature_name The symbol name of the feature
#' @param trail_name The name of the trail
#' @param k Dimensions of bases used to represent the smooth term
#' @keywords internal
#' @import ggplot2
#' @importFrom reshape2 melt
#' @author Daniel C. Ellwanger
.plot_dynamic <- function(x, feature_name, trail_name, k=5) {
  #Pre-flight check
  for(fn in feature_name) {
    .checkFeatureNameExists(x, fn)
  }
  if(!trail_name %in% trailNames(x)) {
    stop("Trail with name '", trail_name, "' does not exists.")
  }
  if(length(trail_name) > 1) {
    trail_name <- trail_name[1]
    warning("Provided more than one trail_name. Dynamic is only shown for ",
            trail_name, ".")
  }

  trail <- x@trails[[trail_name]]
  gp <- ggplot()

  if(length(feature_name) == 1) { #show one feature
    dat <- data.frame(X=trail$ptime / max(trail$ptime) * 100,
                      Y=x[feature_name, trail$samples],
                      STATES=states(x)[trail$samples])
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
    ys <- matrix(nrow = length(trail$ptime), ncol = length(feature_name))
    for(i in seq_along(feature_name)) {
      ys[, i] <- .fitDynamic_def(x = trail$ptime,
                              y = x[feature_name[i], trail$samples],
                              z = states(x)[trail$samples],
                              k = k)$y
    }
    colnames(ys) <- feature_name
    dat <- data.frame(ptime = trail$ptime / max(trail$ptime) * 100, ys)
    dat <- melt(dat, id = "ptime", variable.name = "Feature")

    gp <- ggplot(data=dat,
                 aes_string(x="ptime", y="value", colour="Feature"))
    gp <- gp + geom_line(lwd = .75)
    gp <- gp + theme(axis.line = element_line(colour = "black"))
    gp <- gp + xlab('Pseudotime (%)') + ylab("Expression")

    brks <- pretty(dat$value)
    lbl <- brks
    lbl[lbl == 0] <- "nd"
    gp <- gp + scale_y_continuous(breaks = brks, labels = lbl)
  }
  gp <- gp + labs(subtitle = paste0("Trail: ", trail_name))
  gp
}
