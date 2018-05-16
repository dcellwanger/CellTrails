#' Linear fit
#'
#' Perform simple linear fit
#' @param x.in Predictor variable, numeric vector
#' @param y.in Response variable, numeric vector
#' @param x.out Predictor variable for linear function, numeric vector
#' @return A list containing the following components:
#' @return \item{\code{p.value}}{Anova P-value}
#' \item{\code{r2}}{Adjusted R-squared}
#' \item{\code{x.out}}{Predictor variable for linear function, numeric vector}
#' \item{\code{y.out}}{Response from linear function}
#' @keywords internal
#' @author Daniel C. Ellwanger
.linear_fit <- function(x.in, y.in, x.out=NULL) {
  fit <- lm(y.in ~ x.in)
  p <- anova(fit)$`Pr(>F)`[1]
  r2 <- summary(fit)$adj.r.squared

  if(is.null(x.out)) {
    x.out <- x.in
  }
  y.out <- x.out * fit$coefficients[2] + fit$coefficients[1]

  res <- list()
  res$p.value <- p
  res$r.squared <- r2
  res$x <- x.out
  res$y <- y.out
  res
}

#' Fitting the map surface
#'
#' Fits the smooth surface of CellTrails maps
#' @param X Ordination; numerical matrix
#' @param y Expression values; numerical vector
#' @param npoints Number of points along x and y axis
#' @param weights Cluster states (defines weights based on fraction
#' of non-detects per state)
#' @param knots Number of knots in spline function
#' @return A list containing the following components:
#' \item{\code{fit}}{GAM fit object}
#' \item{\code{grid}}{Fitted values}
#' @author Daniel C. Ellwanger
#' @importFrom mgcv gam
#' @keywords internal
.fit_surface <- function(X, y, npoints=300, weights=NULL,
                         knots=10, rescale=FALSE) {

  # Remove NAs
  #nnas <- complete.cases(X) & !is.na(y)
  #if (!all(nnas)) {
  #  X <- X[nnas, , drop = FALSE]
  #  y <- y[nnas]
  #  w <- w[nnas]
  #}

  x1r <- range(X[, 1])
  x2r <- range(X[, 2])
  x1 <- X[, 1] #.rescale(X[, 1], 0, 1)
  x2 <- X[, 2] #.rescale(X[, 2], 0, 1)
  if(rescale) {
    x1 <- .rescale(x1, 0, 1)
    x2 <- .rescale(x2, 0, 1)
  }
  eq <- formula(paste0("y ~ s(x1, x2, k = ", knots,
                       ", bs = \"tp\", fx = FALSE)"))

  # dynamic nd weight
  w <- NULL
  if(!is.null(weights)) {
    w <- aggregate(y, list(weights),
                   function(x) {sum(x == 0) / length(x)})
    w <- w[match(weights, w[,1]), 2]
    w[y > 0] <- 1 #set measured value weight always to 1
  }

  ## Fit
  fit <- gam(eq, family="gaussian", weights=w, select=TRUE,
             method="REML", gamma=1)
  x1_proj <- seq(min(x1), max(x1), len = npoints)
  x2_proj <- seq(min(x2), max(x2), len = npoints)
  newd <- expand.grid(x1 = x1_proj, x2 = x2_proj)
  surf <- predict(fit, type = "response", newdata = as.data.frame(newd),
                  se=TRUE, block.size=npoints^2)
  fit$se <- surf$se.fit
  surf <- surf$fit
  surf[surf < 0] <- 0
  surf[surf > max(fit$fitted.values)] <- max(fit$fitted.values)
  fit$fitted.values[fit$fitted.values < 0] <- 0

  result <- list()
  result$grid <- data.frame(x1 = .rescale(newd[,1], x1r[1], x1r[2]),
                            x2 = .rescale(newd[,2], x2r[1], x2r[2]),
                            x3 = surf)
  #list(x = x1_proj, y = x2_proj, z = matrix(surf, nrow = npoints))
  result$fit <- fit
  result
}

#' Fitting dynamics
#'
#' Fits expression as a function of pseudotime using generalized
#' additive models
#' @param x Pseudotime values
#' @param y Expression values
#' @param z Sample states
#' @param k Dimensions of the bases used to represent the smooth term
#' @param n.out Number of predicted expression values within pseudotime range
#' @param x.out Predicted expression values for given set of pseudotime values
#' @return A list containing the following components:
#' \item{\code{x}}{Pseudotime}
#' \item{\code{y}}{Fitted values}
#' \item{\code{mod}}{GAM}
#' @importFrom mgcv gam
#' @keywords internal
#' @author Daniel C. Ellwanger
.fitDynamic_def <- function(x, y, z, k=5, n.out=NULL, x.out=NULL) {
  x.ptime <- x
  y.expr <- y
  z.states <- z

  if(var(y.expr) == 0) {
    #return(list(x = NA, y = NA, mod = NA))
    if(is.null(x.out)) {
      x.out <- x.ptime
    }
    return(list(x=x.out, y=rep(mean(y.expr), length(x.out)), mod=NA))
  }

  #w <- nd.weight
  #if(length(w) < length(x.ptime)) {
  #  w <- ifelse(y.expr == 0, nd.weight, 1)
  #}

  # dynamic nd weight
  w <- aggregate(y.expr, list(z.states), function(x) {sum(x == 0) / length(x)})
  w <- w[match(z.states, w[,1]), 2]
  w[y.expr > 0] <- 1 #set measured value weight always to 1

  y <- y.expr
  x <- x.ptime

  eq <- formula(paste0("y ~ te(x, k = ", k, ", bs = \"tp\", fx=TRUE)"))
  mod <- gam(eq, family="gaussian", weights=w, select=TRUE, method="REML",
             gamma=1)

  y.pred <- NA
  if(!is.null(n.out)) {
    s <- seq(from=min(x.ptime), to=max(x.ptime), length.out=n.out)
    y.pred <- predict(mod, type = "response", newdata = data.frame(x = s))
    x.ptime <- s
  } else if (!is.null(x.out)) {
    y.pred <- predict(mod, type = "response", newdata = data.frame(x = x.out))
    x.ptime <- x.out
  } else {
    y.pred <- predict(mod, type = "response")
  }

  y.pred[y.pred < 0] <- 0
  o <- order(x.ptime)
  list(x = x.ptime[o], y = y.pred[o], mod = mod)
}
