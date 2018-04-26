#' @include AllClasses.R
NULL

###############################################################################
### CellTrailsSet
###############################################################################
# #' CellTrailsSet constructor
# #'
# #' Creates new CellTrailsSet object, which is the container for a CellTrails analysis
# #' @usage CellTrails(object)
# #' @param object An \code{\link[Biobase]{ExpressionSet}} object or a numeric matrix with
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
#' Function to coerce a \code{\link[CellTrails]{CellTrailsSet}} object into
#' an \code{\link[Biobase]{ExpressionSet}} object.
#' @param object A \code{\link[CellTrails]{CellTrailsSet}} object
#' @return An \code{\link[Biobase]{ExpressionSet}} object
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
setGeneric("as.ExpressionSet", function(object) standardGeneric("as.ExpressionSet"))
setMethod("as.ExpressionSet", "CellTrailsSet", function(object){
  new("ExpressionSet",
      exprs = exprs(object),
      phenoData = phenoData(object),
      featureData = featureData(object),
      experimentData = experimentData(object),
      annotation = annotation(object),
      protocolData = protocolData(object))
  })

#' Coerce into CellTrailsSet
#'
#' Function to coerce an \code{\link[Biobase]{ExpressionSet}} object or numerical
#' matrix into a \code{\link[CellTrails]{CellTrailsSet}} object.
#' @param object An \code{\link[Biobase]{ExpressionSet}} object or numerical matrix
#' with samples in rows and features in columns
#' @return A \code{\link[CellTrails]{CellTrailsSet}} object
#' @examples
#' # Generate example data
#' dat_num <- simulate_exprs(10, 100)
#' dat_eset <- exDat()
#'
#' # From numerical matrix to CellTrailsSet
#' ctset1 <- as.CellTrailsSet(dat_num)
#' ctset1
#'
#' # From ExpressionSet to CellTrailsSet
#' ctset2 <- as.CellTrailsSet(dat_eset)
#' ctset2
#' @seealso ExpressionSet
#' @docType methods
#' @import Biobase
#' @export
setGeneric("as.CellTrailsSet", function(object) standardGeneric("as.CellTrailsSet"))

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
