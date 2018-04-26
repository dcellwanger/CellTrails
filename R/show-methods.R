#' @include AllClasses.R
NULL

###############################################################################
### CellTrailsSet
###############################################################################
#' Shows content of a CellTrailsSet object
#' @param object A \code{\link[CellTrails]{CellTrailsSet}} object
#' @return \code{show} returns an invisible \code{NULL}
#' @examples
#' # Generate example data
#' ctset <- as.CellTrailsSet(exDat())
#'
#' # Show content
#' show(ctset)
#' @importFrom igraph vcount ecount
#' @import Biobase
#' @author Daniel C. Ellwanger
setMethod("show", "CellTrailsSet", function(object){
  d <- dims(object)
  tfit <- "none"
  if(!is.null(object@trajectory)) {
    tps <- table(object@trajectory$blaze$type)
    tfit <- paste0("MSE=",
                   format(mean(object@trajectory$error), scientific=TRUE, digits=2),
                   " #Branches=", tps[1], " #Terminals=", tps[2])
  }

  out <- paste0("[[ CellTrailsSet ]] \n",
         "assayData: ", d[1], " features, ", d[2], " samples\n",
         "  element names: exprs \n",
         "phenoData: \n",
         "  sampleNames: ", .prettyString(sampleNames(object)), "\n",
         "  varLabels: ", .prettyString(varLabels(object)), "\n",
         "  varMetadata: labelDescription\n",
         "featureData: \n",
         "  featureNames: ", .prettyString(featureNames(object)), "\n",
         "  fvarLabels: ", .prettyString(fvarLabels(object)), "\n",
         "  fvarMetadata: labelDescription\n",
         "mapData: \n",
         "  trajectoryFeatures: ", .prettyString(trajectoryFeatures(object)), "\n",
         "  latentSpace: ", ifelse(is.null(CellTrails::latentSpace(object)), "none", paste0(nrow(CellTrails::latentSpace(object)), " samples, ", ncol(CellTrails::latentSpace(object)), " components")), "\n",
         "  stateTrajectoryGraph: ", ifelse(is.null(stateTrajectoryGraph(object)), "none", paste0("[Component(#Vertices,Edges)]: ", paste0(seq(length(stateTrajectoryGraph(object))), rep("(", length(stateTrajectoryGraph(object))), unlist(lapply(stateTrajectoryGraph(object), function(x){ paste(vcount(x), ecount(x), sep = ',')})), rep(")", length(stateTrajectoryGraph(object))), collapse = " "), "")), "\n",
         "  trajectoryStates: ", .prettyString(trajectoryStates(object)), "\n",
         "  trajectorySamples: ", .prettyString(trajectorySamples(object)), "\n",
         "  trajectoryFit: ", tfit, "\n",
         "  trajectoryLayout: ", ifelse(is.null(trajectoryLayout(object)), "none", "available"), "\n",
         "trailData: \n",
         "  trailNames: ", .prettyString(trailNames(object)))
         #"experiment: \n",
         #"  protocolData, experimentData")
  cat(out)})
