setOldClass(c("igraph"), prototype=structure(list(), class="igraph"))
setClassUnion("matrixOrNULL", members=c("matrix", "NULL"))
setClassUnion("numericOrNULL", members=c("numeric", "NULL"))
setClassUnion("logicalOrNULL", members=c("logical", "NULL"))
setClassUnion("igraphOrNULL", members=c("igraph", "NULL"))
setClassUnion("listOrNULL", members=c("list", "NULL"))
setClassUnion("dfOrNULL", members=c("data.frame", "NULL"))

#' @import stats methods BiocStyle
NULL

###############################################################################
### CellTrailsSet
###############################################################################
#' Class to reconstruct and visualize trajectories from biological samples
#'
#' Container for high-throughput assays and experimental metadata, as well as,
#' trajectory data. This class is derived from \code{\link[Biobase]{ExpressionSet}}.
# #' @section Slots:
# #'  \emph{Inherited}:
# #'  \describe{
# #'    \item{\code{assayData}}{See \code{\link[Biobase]{ExpressionSet}}.}
# #'    \item{\code{phenoData}}{See \code{\link[Biobase]{ExpressionSet}}.}
# #'    \item{\code{featureData}}{See \code{\link[Biobase]{ExpressionSet}}.}
# #'    \item{\code{experimentData}}{See \code{\link[Biobase]{ExpressionSet}}.}
# #'    \item{\code{annotation}}{See \code{\link[Biobase]{ExpressionSet}}.}
# #'    \item{\code{protocolData}}{See \code{\link[Biobase]{ExpressionSet}}.}
# #'  }
# #'  \emph{Additional}:
# #'  \describe{
# #'    \item{\code{useFeature}}{Object of class \code{\link[base]{logical}}, containing indices of selected trajectory features.}
# #'    \item{\code{useSample}}{Object of class \code{\link[base]{logical}}, containing indices of selected trajectory samples.}
# #'    \item{\code{eigen}}{Object of class \code{\link[base]{list}}, containing the latent space and eigenvalues from the spectral embedding.}
# #'    \item{\code{spanForest}}{Object of class \code{\link[base]{list}}, containing the \code{\link[igraph]{igraph}} components of the maximum interface spanning graph.}
# #'    \item{\code{trajectory}}{Object of class \code{\link[base]{list}}, containing the trajectory fit.}
## '    \item{\code{layout}}{Object of class \code{\link[base]{data.frame}},
# #'    containing the geographic coordinates of the CellTrails maps.}
# #'    \item{\code{trails}}{Object of class \code{\link[base]{list}}, containing trails information.}
# #'  }
# #' @slot .__classVersion__ See \code{\link[Biobase]{ExpressionSet}}
#' @slot assayData See \code{\link[Biobase]{ExpressionSet}}
#' @slot phenoData See \code{\link[Biobase]{ExpressionSet}}
#' @slot featureData See \code{\link[Biobase]{ExpressionSet}}
#' @slot experimentData See \code{\link[Biobase]{ExpressionSet}}
#' @slot annotation See \code{\link[Biobase]{ExpressionSet}}
#' @slot protocolData See \code{\link[Biobase]{ExpressionSet}}
#' @slot useFeature Object of class \code{\link[base]{logical}}, containing indices of selected trajectory features
#' @slot useSample Object of class \code{\link[base]{logical}}, containing indices of selected trajectory samples
#' @slot eigen Object of class \code{\link[base]{list}}, containing the latent space and eigenvalues from the spectral embedding.
#' @slot spanForest Object of class \code{\link[base]{list}}, containing the \code{\link[igraph]{igraph}} components of the maximum interface spanning graph
#' @slot trajectory Object of class \code{\link[base]{list}}, containing the trajectory fit
#' @slot layout Object of class \code{\link[base]{data.frame}}, containing the geographic coordinates of the CellTrails maps
#' @slot trails Object of class \code{\link[base]{list}}, containing trails information
#' @details
#'  \strong{Class-specific methods}
#'  \describe{
#'    \item{\code{show(object)}}{Informatively display object contents.}
#'    \item{\code{object[(index)}}{Conducts subsetting of expression data component.}
#'    \item{\code{trajectoryFeatures(object)}, \code{trajectoryFeatures(object)<-}}{Get/Set features for trajectory reconstruction.}
#'    \item{\code{latentSpace(object)}}{Get latent space from dimension dimension reduction.}
#'    \item{\code{eigenvalues(object)}}{Get eigenvalues from dimension reduction.}
#'    \item{\code{states(object)}}{Get assigned state per sample.}
#'    \item{\code{stateTrajectoryGraph(object)}}{Get components of the trajectory graph.}
#'    \item{\code{trajectoryStates(object)}}{Get states contained in selected trajectory graph component.}
#'    \item{\code{trajectorySamples(object)}}{Get samples contained in selected trajectory graph component.}
#'    \item{\code{trajectoryFit(object)}}{Get trajectory fit information.}
#'    \item{\code{trajectoryLayout(object)}}{Get trajectory layout information.}
#'  }
#' @importClassesFrom Biobase eSet ExpressionSet
#' @import Biobase
#' @name CellTrailsSet
#' @rdname CellTrailsSet-class
#' @aliases CellTrailsSet-class
# #' @exportClass CellTrailsSet
#' @export
#' @seealso \code{\link[Biobase]{ExpressionSet}}, \code{\link[Biobase]{eSet}}
#' @author Daniel C. Ellwanger
setClass("CellTrailsSet",
         contains = "ExpressionSet",
         slots = c(useFeature = "logicalOrNULL",
                   useSample = "logicalOrNULL",
                   eigen = "listOrNULL",
                   spanForest = "listOrNULL",
                   trajectory = "listOrNULL",
                   layout = "dfOrNULL",
                   trails = "listOrNULL"
         )
)

# setMethod("initialize", "CellTrailsSet", function(.Object, ...) {
#   print("BLUBB")
#   #.Object <- new("CellTrailsSet") #callNextMethod()
#   if(class(assayData) == "ExpressionSet") {
#     .Object@assayData = assayData@assayData
#     .Object@phenoData = assayData@phenoData
#     .Object@featureData = assayData@featureData
#     .Object@experimentData = assayData@experimentData
#     .Object@annotation = assayData@annotation
#     .Object@protocolData = assayData@protocolData
#     .Object@useFeature = rep(TRUE, nrow(exprs(assayData)))
#     .Object@eigen = list()
#     .Object@spanForest = NULL
#     .Object@trajectory = NULL
#     .Object@layout = NULL
#     .Object@trails = list() }
#   .Object
# })


# CellTrailsSet <- function(...) {
#   if(class(assayData) == "ExpressionSet") {
#     new("CellTrailsSet",
#         assayData = exprs(assayData),
#         phenoData = assayData@phenoData,
#         featureData = assayData@featureData,
#         experimentData = assayData@experimentData,
#         annotation = assayData@annotation,
#         protocolData = assayData@protocolData,
#         useFeature = rep(TRUE, nrow(exprs(assayData))),
#         eigen = list(),
#         spanForest = NULL,
#         trajectory = NULL,
#         layout = NULL,
#         trails = list()
#     )
#   } else if(class(assayData) == "matrix") {
#     CellTrailsSet(assayData = ExpressionSet(assayData))
#   } else {
#     new("CellTrailsSet", ...)
#   }
# }

#setMethod("$", "CellTrailsSet", function(x, FUN){colnames(x@phenoData)})

###############################################################################
### CellTrailsSpectrum
###############################################################################
#' Class to determine the dimensions of the manifold
#'
#' Container for spectral embedding related information; not
#' to be initialized by user directly.
# #' @section Slots:
# #'  \describe{
# #'    \item{\code{frac}}{Object of class \code{\link[base]{numeric}} containing the number of eigengaps used to perform linear fit}
# #'    \item{\code{n}}{Object of class \code{\link[base]{numeric}} containing the number of relevant dimensions}
# #'    \item{\code{cs}}{Object of class \code{\link[base]{numeric}} contaning the cumulative sum of eigengaps}
# #'    \item{\code{fit}}{Object of class \code{\link[base]{list}} containing the results of the linear fit}
# #'  }
#' @slot frac Object of class \code{\link[base]{numeric}} containing the number of eigengaps used to perform linear fit
#' @slot n Object of class \code{\link[base]{numeric}} containing the number of relevant dimensions
#' @slot cs Object of class \code{\link[base]{numeric}} contaning the cumulative sum of eigengaps
#' @slot fit Object of class \code{\link[base]{list}} containing the results of the linear fit
#' @name CellTrailsSpectrum
#' @rdname CellTrailsSpectrum-class
#' @aliases CellTrailsSpectrum-class
# #' @exportClass CellTrailsSpectrum
#' @export
#' @author Daniel C. Ellwanger
setClass("CellTrailsSpectrum",
         slots = c(frac = "numeric",
                   n = "numeric",
                   cs = "numeric",
                   fit = "list")
)
