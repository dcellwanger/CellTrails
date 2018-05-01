#' @include AllClasses.R
NULL

###############################################################################
### CellTrailsSet
###############################################################################
# #' CellTrailsSet constructor
# #'
# #' Creates new CellTrailsSet object, which is the container for a
# #' CellTrails analysis
# #' @usage CellTrails(object)
# #' @param object An \code{ExpressionSet} object or a numeric matrix with
# #' expression values (rows = features, columns = samples).
# #' @examples
# #' data(gga_E15_utricle)
# #' # Initialize from ExpressionSet
# #' ctset <- CellTrails(gga_E15_utricle)
# #' # Initialize from numeric matrix
# #' dat <- exprs(gga_E15_utricle)
# #' ctset <- CellTrails(dat)
# #' @docType methods
# #' @rdname CellTrails-method
# #' @aliases CellTrails,ExpressionSet-method
# #' @aliases CellTrails,matrix-method
# #' @return A new \code{CellTrailsSet} object
# #' @import Biobase
# #' @export
# setGeneric("CellTrails", function(object) standardGeneric("CellTrails"))
# setMethod("CellTrails", "ExpressionSet", function(object){
#   as.CellTrailsSet(object)
# })
#
# setMethod("CellTrails", "matrix", function(object){
#   CellTrails(ExpressionSet(object))})

#' Coerce CellTrailsSet into ExpressionSet
#'
#' Function to coerce a \code{CellTrailsSet} object into
#' an \code{ExpressionSet} object.
#' @param object A \code{CellTrailsSet} object
#' @return An \code{ExpressionSet} object
#' @details The resulting \code{ExpressionSet} object
#' contains additional phenoData information from the CellTrails
#' analysis (e.g., the determined states, the pseudotime
#' along each trail).
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @examples
#' # Generate example data
#' dat <- exDat()
#' dat
#' ctset <- as.CellTrailsSet(dat)
#'
#' # From CellTrailsSet to ExpressionSet
#' dat2 <- as.ExpressionSet(ctset)
#' dat2
#' @seealso ExpressionSet
#' @docType methods
#' @aliases as.ExpressionSet,CellTrailsSet-method
#' @import Biobase
#' @export
setGeneric("as.ExpressionSet", function(object)
  standardGeneric("as.ExpressionSet"))
setMethod("as.ExpressionSet", "CellTrailsSet", function(object){
  pd <- pData(object)
  colnames(pd) <- gsub(pattern="^STATE$", x=colnames(pd),
                       replacement="CellTrails.STATE")
  for(n in trailNames(object)) {
    pt <- trailPseudotime(object)[[n]]
    cn <- paste0("CellTrails", ".", n)
    pd[[cn]] <- NA
    pd[trailSamples(object)[[n]], cn] <- pt
  }

  new("ExpressionSet",
      exprs = exprs(object),
      phenoData = new("AnnotatedDataFrame", pd),
      featureData = featureData(object),
      experimentData = experimentData(object),
      annotation = annotation(object),
      protocolData = protocolData(object))
  })

#' Coerce CellTrailsSet into SingleCellExperiment
#'
#' Function to coerce a \code{CellTrailsSet} object into
#' a \code{SingleCellExperiment} object.
#' @param object A \code{CellTrailsSet} object
#' @return A \code{SingleCellExperiment} object
#' @details The resulting \code{SingleCellExperiment} object
#' contains additional \code{colData} information from the CellTrails
#' analysis (e.g., the determined states, the pseudotime
#' along each trail).
#' @seealso \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' @examples
#' # Generate example data
#' dat <- exDat()
#' dat
#' ctset <- as.CellTrailsSet(dat)
#'
#' # From CellTrailsSet to ExpressionSet
#' dat2 <- as.ExpressionSet(ctset)
#' dat2
#' @docType methods
#' @aliases as.SingleCellExperiment,CellTrailsSet-method
#' @import Biobase
#' @importFrom S4Vectors DataFrame SimpleList
#' @export
setGeneric("as.SingleCellExperiment", function(object)
  standardGeneric("as.SingleCellExperiment"))
setMethod("as.SingleCellExperiment", "CellTrailsSet", function(object){
  #m  <- matrix(ncol=2, nrow=ncol(object),
  #             dimnames=list(sampleNames(object), c("X1", "X2")))
  #tl <- as.matrix(trajectoryLayout(object))
  #m[rownames(tl), ] <- tl

  pd <- pData(object)
  colnames(pd) <- gsub(pattern="^STATE$", x=colnames(pd),
                       replacement="CellTrails.STATE")
  for(n in trailNames(object)) {
    pt <- trailPseudotime(object)[[n]]
    cn <- paste0("CellTrails", ".", n)
    pd[[cn]] <- NA
    pd[trailSamples(object)[[n]], cn] <- pt
  }

  SingleCellExperiment(
    assay = list(logcounts = exprs(object)),
    colData = DataFrame(pd),
    rowData = DataFrame(fData(object)),
    reducedDims = SimpleList(
      CellTrails.SE = latentSpace(object)))
      #CellTrails.TL = m))
})

#' Coerce into CellTrailsSet
#'
#' Function to coerce an \code{ExpressionSet} object, numerical
#' matrix or \code{SingleCellExperiment} object into a \code{CellTrailsSet}
#' object.
#' @param object An \code{ExpressionSet} object or numerical matrix
#' with samples in rows and features in columns
#' @return A \code{CellTrailsSet} object
#' @details We recommend to perform a CellTrails analysis on
#' log-normalized expression data without spike-in controls. To properly coerce
#' a \code{SingleCellExperiment} into a \code{CellTrailsSet} a \code{logcounts}
#' assay needs to be defined (see examples).
#' @seealso\ code{\link[CellTrails]{CellTrailsSet}}
#' \code{\link[Biobase]{ExpressionSet}}
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' @examples
#' # Generate example data
#' dat_num <- simulate_exprs(10, 100)
#' dat_eset <- exDat()
#' dat_list <- list(logcounts=dat_num)
#' dat_sce <- SingleCellExperiment::SingleCellExperiment(assay=dat_list)
#'
#' # From numerical matrix to CellTrailsSet
#' ctset1 <- as.CellTrailsSet(dat_num)
#' ctset1
#'
#' # From ExpressionSet to CellTrailsSet
#' ctset2 <- as.CellTrailsSet(dat_eset)
#' ctset2
#'
#' # From SingleCellExperiment to CellTrailsSet
#' ctset3 <- as.CellTrailsSet(dat_sce)
#' @seealso ExpressionSet
#' @docType methods
#' @import Biobase
#' @importFrom SingleCellExperiment SingleCellExperiment colData rowData logcounts
#' @export
setGeneric("as.CellTrailsSet", function(object)
  standardGeneric("as.CellTrailsSet"))

#' @rdname as.CellTrailsSet
#' @aliases as.CellTrailsSet,ExpressionSet-method
setMethod("as.CellTrailsSet", "ExpressionSet", function(object){
  new("CellTrailsSet",
    assayData = object@assayData,
    phenoData = object@phenoData,
    featureData = object@featureData,
    experimentData = object@experimentData,
    annotation = object@annotation,
    protocolData = object@protocolData,
    useFeature = rep(TRUE, nrow(exprs(object))),
    eigen = list(),
    spanForest = NULL,
    trajectory = NULL,
    layout = NULL,
    trails = list())
})

#' @rdname as.CellTrailsSet
#' @aliases as.CellTrailsSet,matrix-method
setMethod("as.CellTrailsSet", "matrix", function(object){
  as.CellTrailsSet(ExpressionSet(object))
})

#' @rdname as.CellTrailsSet
#' @aliases as.CellTrailsSet,SingleCellExperiment-method
setMethod("as.CellTrailsSet", "SingleCellExperiment", function(object){
  eset <- ExpressionSet(logcounts(object))
  ncdat <- ncol(colData(object))
  nrdat <- ncol(rowData(object))
  if(ncdat > 0) { #phenoData exists
    df <- data.frame(colData(object),
                     row.names=sampleNames(eset))
    phenoData(eset) <- new("AnnotatedDataFrame", df)
  }
  if(nrdat > 0) { #featureData exists
    df <- data.frame(rowData(object),
                     row.names=featureNames(eset))
    featureData(eset) <- new("AnnotatedDataFrame", df)
  }
  as.CellTrailsSet(eset)
})
