# #' @include AllClasses.R
#NULL

###############################################################################
### CellTrails
###############################################################################
#' Shows relevant content of a SingleCellExperiment object for a CellTrails
#' analysis
#' @param object A \code{SingleCellExperiment} object
#' @return \code{showTrajInfo} returns an invisible \code{NULL}
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Show content
#' showTrajInfo(dat)
#' @importFrom igraph vcount ecount
#' @importFrom SummarizedExperiment rowData
#' @export
#' @docType methods
#' @aliases showTrajInfo,SingleCellExperiment-method
#' @author Daniel C. Ellwanger
setGeneric("showTrajInfo", function(object) standardGeneric("showTrajInfo"))
setMethod("showTrajInfo", "SingleCellExperiment", function(object){
  d <- dim(.exprs(object))
  lmarks <- "none"
  if(!is.null(.trajGraph(object))) {
    tps <- table(.trajLandmark(object, type="type"))
    lmarks <- paste0(" #Branches=", tps["B"],
                     " #Terminals=", tps["H"],
                     " #User=", tps["U"])
  }
  trajs <- ifelse(is.null(.spanForest(object)),
                  "none",
                  paste0("[Component(#Vertices,Edges)]: ",
                         paste0(seq(length(.spanForest(object))),
                                rep("(", length(.spanForest(object))),
                                unlist(lapply(.spanForest(object),
                                              function(x){paste(vcount(x),
                                                                ecount(x),
                                                                sep = ',')})),
                                rep(")", length(.spanForest(object))),
                                collapse = " "), ""))

  mse <- ifelse(is.null(trajResiduals(object)), "NA",
                format(mean(trajResiduals(object), na.rm=TRUE),
                       scientific=TRUE, digits=2))

  latspec <- ifelse(is.null(CellTrails::latentSpace(object)),
                    "none",
                    paste0(nrow(CellTrails::latentSpace(object)),
                           " samples, ",
                           ncol(CellTrails::latentSpace(object)),
                           " dimensions"))

  out <- paste0("[[ CellTrails ]] \n",
         "logcounts: ", d[1], " features, ", d[2], " samples\n",
         "Pheno data: \n",
         "  sampleNames: ", .prettyString(sampleNames(object)), "\n",
         "  phenoNames: ", .prettyString(phenoNames(object)), "\n",
         "Feature data: \n",
         "  featureNames: ", .prettyString(featureNames(object)), "\n",
         "  rowData: ", .prettyString(colnames(rowData(object))), "\n",
         "Trajectory data: \n",
         "  trajFeatureNames: ", .prettyString(trajFeatureNames(object)), "\n",
         "  latentSpace: ", latspec, "\n",
         "Trajectories: ", trajs, "\n",
         #"  states: ", .prettyString(trajectoryStates(object)), "\n",
         "  trajSampleNames: ", .prettyString(trajSampleNames(object)), "\n",
         "  trajResiduals: MSE=", mse, "\n",
         "  landmarks: ", lmarks, "\n",
         "  trajLayout: ", ifelse(is.null(trajLayout(object)),
                                  "none", "available"), "\n",
         "Trail data: \n",
         "  trailNames: ", .prettyString(trailNames(object)))
  cat(out)})
