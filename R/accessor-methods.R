#' @include AllClasses.R
NULL

###############################################################################
# SingleCellExperiment
###############################################################################
###############################################################################
# Internal
###############################################################################
#' SET trajectory features indicator
#'
#' Sets indicator if feature should be used for trajectory
#' reconstruction.
#' @param object An object of class \code{SingleCellExperiment}
#' @param value A logical vector
#' @return An updated object of class \code{SingleCellExperiment}
#' @docType methods
#' @aliases .useFeature<-,SingleCellExperiment-method
#' @keywords internal
#' @author Daniel C. Ellwanger
setGeneric(".useFeature<-", function(object, value)
  standardGeneric(".useFeature<-"))
setMethod(".useFeature<-", "SingleCellExperiment", function(object, value){
  object@int_elementMetadata$CellTrails.isSelected <- value
  object})

#' GET trajectory features indicator
#'
#' Indicates if feature should be used for trajectory
#' reconstruction. Spike-in controls are removed.
#' @param object An object of class \code{SingleCellExperiment}
#' @return A \code{logical} vector
#' @import SingleCellExperiment
#' @docType methods
#' @aliases .useFeature,SingleCellExperiment-method
#' @keywords internal
#' @author Daniel C. Ellwanger
setGeneric(".useFeature", function(object)
  standardGeneric(".useFeature"))
setMethod(".useFeature", "SingleCellExperiment", function(object){
  spike <- isSpike(object)
  if(is.null(spike)){
    spike <- rep(FALSE, nrow(object))
  }
  uF <- object@int_elementMetadata$CellTrails.isSelected
  if(is.null(uF)){
    uF <- rep(TRUE, nrow(object))
  }
  uF & !spike})

#' SET trajectory samples indicator
#'
#' Sets indicator if sample was used for trajectory reconstruction.
#' @param object An object of class \code{SingleCellExperiment}
#' @param value A logical vector
#' @return An updated object of class \code{SingleCellExperiment}
#' @docType methods
#' @aliases useSample<-,SingleCellExperiment-method
#' @keywords internal
#' @author Daniel C. Ellwanger
setGeneric(".useSample<-", function(object, value)
  standardGeneric(".useSample<-"))
setMethod(".useSample<-", "SingleCellExperiment", function(object, value){
  object@int_colData$CellTrails.isSelected <- value
  object})

#' GET trajectory samples indicator
#'
#' Indicates if sample was used for trajectory reconstruction.
#' @param object An object of class \code{SingleCellExperiment}
#' @return A \code{logical} vector
#' @docType methods
#' @aliases .useSample,SingleCellExperiment-method
#' @keywords internal
#' @author Daniel C. Ellwanger
setGeneric(".useSample", function(object)
  standardGeneric(".useSample"))
setMethod(".useSample", "SingleCellExperiment", function(object){
  uS <- object@int_colData$CellTrails.isSelected
  if(is.null(uS)){
    uS <- rep(TRUE, ncol(object))
  }
  uS})

#' GET expression matrix
#'
#' Retrieve numeric matrix of expression values for processing in
#' CellTrails. This wrapper function ensures that all functions in the
#' package receive the proper assay from the
#' \code{SingleCellExperiment} object.
#' @param object An object of class \code{SingleCellExperiment}
#' @return A \code{numeric} matrix
#' @docType methods
#' @aliases .exprs,SingleCellExperiment-method
#' @keywords internal
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @author Daniel C. Ellwanger
setGeneric(".exprs", function(object)
  standardGeneric(".exprs"))
setMethod(".exprs", "SingleCellExperiment", function(object){
  assay(object, "logcounts")})

#' SET state trajectory graph
#'
#' Sets graph object spanning all states (spanning forest)
#' to \code{SingleCellExperiment} object
#' @param object A \code{SingleCellExperiment} object
#' @param value A \code{list} object with an \code{igraph} object
#' per component of the spanning forest
#' @return An updated object of class \code{SingleCellExperiment}
#' @docType methods
#' @aliases .spanForest<-,SingleCellExperiment-method
#' @keywords internal
#' @author Daniel C. Ellwanger
setGeneric(".spanForest<-", function(object, value)
  standardGeneric(".spanForest<-"))
setMethod(".spanForest<-", "SingleCellExperiment", function(object, value){
  object@int_metadata$CellTrails$spanForest <- value
  object})

#' GET state trajectory graph
#'
#' Returns graph object spanning all states (spanning forest)
#' @param object A \code{SingleCellExperiment} object
#' @return A \code{list} object with an \code{igraph} object
#' per component of the spanning forest
#' @docType methods
#' @aliases .spanForest,SingleCellExperiment-method
#' @keywords internal
#' @author Daniel C. Ellwanger
setGeneric(".spanForest", function(object)
  standardGeneric(".spanForest"))
setMethod(".spanForest", "SingleCellExperiment", function(object){
  object@int_metadata$CellTrails$spanForest
})

#' SET trajectory graph
#'
#' Stores trajectory graph in a \code{SingleCellExperiment} object
#' @param object A \code{SingleCellExperiment} object
#' @param value A \code{igraph} object
#' @return An updated object of class \code{SingleCellExperiment}
#' @docType methods
#' @aliases .trajGraph<-,SingleCellExperiment-method
#' @keywords internal
#' @author Daniel C. Ellwanger
setGeneric(".trajGraph<-", function(object, value)
  standardGeneric(".trajGraph<-"))
setMethod(".trajGraph<-", "SingleCellExperiment", function(object, value){
  object@int_metadata$CellTrails$trajGraph <- value
  object})

#' GET trajectory graph
#'
#' Returns trajectory graph from a \code{SingleCellExperiment} object
#' @param object A \code{SingleCellExperiment} object
#' @return A \code{igraph} object
#' @docType methods
#' @aliases .trajGraph,SingleCellExperiment-method
#' @keywords internal
#' @author Daniel C. Ellwanger
setGeneric(".trajGraph", function(object)
  standardGeneric(".trajGraph"))
setMethod(".trajGraph", "SingleCellExperiment", function(object){
  object@int_metadata$CellTrails$trajGraph
})

#' SET trajectory fitting residuals
#'
#' Stores trajectory fitting residuals in \code{SingleCellExperiment} object
#' @param object A \code{SingleCellExperiment} object
#' @param value A \code{numeric} vector
#' @return An updated object of class \code{SingleCellExperiment}
#' @docType methods
#' @aliases .trajResiduals<-,SingleCellExperiment-method
#' @keywords internal
#' @author Daniel C. Ellwanger
setGeneric(".trajResiduals<-", function(object, value)
  standardGeneric(".trajResiduals<-"))
setMethod(".trajResiduals<-", "SingleCellExperiment", function(object, value){
  trsid <- rep(NA, ncol(object))
  trsid[.useSample(object)] <- value
  object@int_colData$CellTrails.residuals <- trsid
  object})

#' GET trajectory fitting residuals
#'
#' Returns trajectory fitting residuals from \code{SingleCellExperiment} object
#' @param object A \code{SingleCellExperiment} object
#' @return A \code{numeric} vector
#' @details The trajectory fitting deviation is defined as the
#' vector rejection from a sample in the latent space to the trajectory
#' backbone. The trajectory backbone is defined by a tree spanning all
#' relevant states. Samples get orthogonally projected onto straight lines
#' connecting related states. This function quantifies the distance between
#' the actual positon of a sample in the latent space and its projectd position
#' on the trajectory backbone. In other words, the higher the distance, the
#' higher its deviation (residual) from the trajectory fit. This function
#' returns all residuals for each projected sample. Residuals of samples which
#' were exluded for trajectory reconstruction are \code{NA}.
#' @seealso \code{fitTrajectory} \code{trajSampleNames}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' trajResiduals(exSCE)[seq_len(5)]
#' @docType methods
#' @aliases trajResiduals,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trajResiduals", function(object)
  standardGeneric("trajResiduals"))
setMethod("trajResiduals", "SingleCellExperiment", function(object){
  object@int_colData$CellTrails.residuals
})

# #' SET tSNE parameters
# #'
# #' Stores tSNE parameters in \code{SingleCellExperiment} object
# #' @param object A \code{SingleCellExperiment} object
# #' @param value A \code{numeric} vector ('seed', 'perplexity')
# #' @return An updated object of class \code{SingleCellExperiment}
# #' @docType methods
# #' @aliases .tsneParams<-,SingleCellExperiment-method
# #' @keywords internal
# #' @author Daniel C. Ellwanger
# setGeneric(".tsneParams<-", function(object, value)
#   standardGeneric(".tsneParams<-"))
# setMethod(".tsneParams<-", "SingleCellExperiment", function(object, value){
#   if(is.null(.tsneParams(object))) {
#     names(value) <- c("seed", "perplexity")
#   }
#   object@int_metadata$CellTrails$tsne_params <- value
#   object})

# #' GET tSNE parameters
# #'
# #' Returns tSNE parameters from \code{SingleCellExperiment} object
# #' @param object A \code{SingleCellExperiment} object
# #' @return A \code{numeric} vector
# #' @docType methods
# #' @aliases .tsneParams,SingleCellExperiment-method
# #' @keywords internal
# #' @author Daniel C. Ellwanger
# setGeneric(".tsneParams", function(object)
#   standardGeneric(".tsneParams"))
# setMethod(".tsneParams", "SingleCellExperiment", function(object){
#   object@int_metadata$CellTrails$tsne_params
# })

#' SET trajectory landmark annotation
#'
#' Stores information on trajectory landmarks
#' in a \code{SingleCellExperiment} object
#' @param object A \code{SingleCellExperiment} object
#' @param value A vector of any type
#' @param type A character in {"type", "id", "shape"}
#' @return An updated object of class \code{SingleCellExperiment}
#' @docType methods
#' @aliases .trajLandmark<-,SingleCellExperiment-method
#' @keywords internal
#' @author Daniel C. Ellwanger
setGeneric(".trajLandmark<-", function(object,
                                       type=c("type", "id", "shape"), value)
  standardGeneric(".trajLandmark<-"))
setMethod(".trajLandmark<-", "SingleCellExperiment", function(object,
                                                              type, value){
  lndmrk <- rep(NA, ncol(object))
  if(is.factor(value)) {
    lndmrk <- factor(lndmrk, levels=levels(value))
  }
  lndmrk[.useSample(object)] <- value
  type <- .capitalize(type)
  type <- paste0("CellTrails.lndmrk", type)
  object@int_colData[, type] <- lndmrk
  object})

#' GET trajectory landmark annotation
#'
#' Returns trajectory landmark information from
#' \code{SingleCellExperiment} object
#' @param object A \code{SingleCellExperiment} object
#' @param type A character in {"type", "id", "shape"}
#' @return A vector of any type
#' @docType methods
#' @aliases .trajLandmark,SingleCellExperiment-method
#' @keywords internal
#' @author Daniel C. Ellwanger
setGeneric(".trajLandmark", function(object, type=c("type", "id", "shape"))
  standardGeneric(".trajLandmark"))
setMethod(".trajLandmark", "SingleCellExperiment", function(object, type){
  type <- .capitalize(type)
  type <- paste0("CellTrails.lndmrk", type)
  if(type %in% colnames(object@int_colData)) {
    object@int_colData[, type]
  } else {
    NULL
  }
})

#' GET phenotype values
#'
#' Returns phenotype values from
#' \code{SingleCellExperiment} object
#' @param object A \code{SingleCellExperiment} object
#' @param name Name of phenotype
#' @return A vector of any type
#' @details Wrapper for colDat(object)[, name] which also accesses
#' internal metadata (e.g., landmarks).
#' @docType methods
#' @aliases .pheno,SingleCellExperiment-method
#' @keywords internal
#' @author Daniel C. Ellwanger
setGeneric(".pheno", function(object, name)
  standardGeneric(".pheno"))
setMethod(".pheno", "SingleCellExperiment", function(object, name){
  lname <- tolower(name)
  cname <- paste0("CellTrails.", name)
  if(lname == "state") {
    d <- states(object)
  } else if(lname == "landmark") {
    d <- landmarks(object)
  } else if(name %in% colnames(colData(object))){
      d <- colData(object)[, name]
  } else if(cname %in% colnames(colData(object))){
    d <- colData(object)[, cname]
  } else {
    d <- NULL
  }
  d})

###############################################################################
# Exported
###############################################################################
#' GET feature names
#'
#' Retrieve feature names from a \code{SingleCellExperiment} object
#' @param object An object of class \code{SingleCellExperiment}
#' @return A \code{character} vector
#' @details Wrapper for \code{rownames(object)}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' featureNames(exSCE)
#' @seealso \code{SingleCellExperiment}
#' @importFrom Biobase featureNames
#' @docType methods
#' @aliases featureNames,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
#setGeneric("featureNames", function(object)
#  standardGeneric("featureNames"))
setMethod("featureNames", "SingleCellExperiment", function(object){
  rownames(object)})

#' GET sample names
#'
#' Retrieve sample names from a \code{SingleCellExperiment} object
#' @param object An object of class \code{SingleCellExperiment}
#' @return A \code{character} vector
#' @details Wrapper for \code{colnames(object)}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' sampleNames(exSCE)[seq_len(5)]
#' @seealso \code{SingleCellExperiment}
#' @importFrom Biobase sampleNames
#' @docType methods
#' @aliases sampleNames,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
#setGeneric("sampleNames", function(object)
#  standardGeneric("sampleNames"))
setMethod("sampleNames", "SingleCellExperiment", function(object){
  colnames(object)})

#' GET phenotype names
#'
#' Retrieve phenotype names from a \code{SingleCellExperiment} object
#' @param object An object of class \code{SingleCellExperiment}
#' @return A \code{character} vector
#' @details Wrapper for \code{colnames(colData(object))}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' phenoNames(exSCE)
#' @seealso \code{SingleCellExperiment}
#' @docType methods
#' @aliases phenoNames,SingleCellExperiment-method
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @export
#' @author Daniel C. Ellwanger
setGeneric("phenoNames", function(object)
  standardGeneric("phenoNames"))
setMethod("phenoNames", "SingleCellExperiment", function(object){
  #nsm <- colnames(colData(object))
  nsm <- gsub(x=colnames(colData(object)), "CellTrails.", "")
  nl <- length(landmarks(object))
  if(nl > 0) {
    c(nsm, "landmark")
  } else {
    nsm
  }})

#' SET states
#'
#' Sets states to a \code{SingleCellExperiment} object
#' @param object An object of class \code{SingleCellExperiment}
#' @param value A numeric, character or factor vector
#' @return An updated object of class \code{SingleCellExperiment}
#' @details State information is added to a
#' \code{SingleCellExperiment} object via \code{colData}. If the
#' vector containing the cluster assignments is numeric, the prefix
#' "S" is added and the vector is converted to type factor.
#' @seealso colData
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Assign clusters
#' cl <- kmeans(logcounts(exSCE), centers=10)$cluster
#' states(exSCE) <- cl
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData<-
#' @docType methods
#' @aliases states<-,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("states<-", function(object, value)
  standardGeneric("states<-"))
setMethod("states<-", "SingleCellExperiment", function(object, value){
  if(is.numeric(value)) {
    s <- unique(value)
    o <- order(s)
    s <- paste0("S", s)
    value <- paste0("S", value)
    value <- factor(value, levels=s[o])
  }
  if(!is.factor(value)) {
    value <- factor(value)
  }
  colData(object)$CellTrails.state <- value
  object})

#' GET states
#'
#' Retrieve computed states from a \code{SingleCellExperiment} object
#' @param object An object of class \code{SingleCellExperiment}
#' @return A \code{factor} vector
#' @details State information is extracted from \code{colData};
#' factor levels are alphanumerically ordered by ID.
#' @examples
#' # Example data
#' data(exSCE)
#'
#' states(exSCE)[seq_len(5)]
#' @seealso \code{SingleCellExperiment} \code{findStates}
#' @docType methods
#' @aliases states,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("states", function(object)
  standardGeneric("states"))
setMethod("states", "SingleCellExperiment", function(object){
  colData(object)$CellTrails.state})

#' SET trajectory features by name
#'
#' Function to set trajectory features by name
#' @param object An object of class \code{SingleCellExperiment}
#' @param value A character vector
#' @return An updated object of class \code{SingleCellExperiment}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Set trajectory features
#' trajFeatureNames(exSCE) <- rownames(exSCE)[seq_len(5)]
#' @docType methods
#' @aliases trajFeatureNames<-,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trajFeatureNames<-", function(object, value)
  standardGeneric("trajFeatureNames<-"))
setMethod("trajFeatureNames<-", "SingleCellExperiment", function(object,
                                                                 value){
  .featureNameExists(x=object, feature_name=value)
  .useFeature(object) <- rownames(object) %in% value
  object})

#' GET trajectory feature names
#'
#' Retrieve names of features that were selected for trajectory reconstruction
#' from a \code{SingleCellExperiment} object.
#' @param object An object of class \code{SingleCellExperiment}
#' @return An object of class \code{character}
#' @details Features can be selected prior to trajectory inference.
#' This method retrieves the user-defined features from a
#' \code{SingleCellExperiment} object. The return value is a character
#' vector containing the feature names.
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Get trajectory features
#' trajFeatureNames(exSCE)[seq_len(5)]
#' @docType methods
#' @aliases trajFeatureNames,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trajFeatureNames", function(object)
  standardGeneric("trajFeatureNames"))
setMethod("trajFeatureNames", "SingleCellExperiment", function(object){
  rownames(object)[.useFeature(object)]})

#' GET trajectory sample names
#'
#' Retrieve names of samples that were aligned onto the trajectory
#' from a \code{SingleCellExperiment} object.
#' @param object An object of class \code{SingleCellExperiment}
#' @return An object of class \code{character}
#' @details A trajectory graph can be initially a forest. Trajectory fitting
#' is performed on one component. This function returns the names of the
#' samples which are member of the selected component.
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Get trajectory samples
#' trajSampleNames(exSCE)[seq_len(5)]
#' @docType methods
#' @aliases trajSampleNames,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trajSampleNames", function(object)
  standardGeneric("trajSampleNames"))
setMethod("trajSampleNames", "SingleCellExperiment", function(object){
  colnames(object)[.useSample(object)]})

#' SET latent space
#'
#' Set CellTrails' latent space to a \code{SingleCellExperiment} object.
#' @param object A \code{SingleCellExperiment} object
#' @param value A numeric matrix with rows are samples and columns are
#' components
#' @return An updated object of class \code{SingleCellExperiment}
#' @details Rows need to be samples and columns to be \emph{d} components
#' (spanning the lower-dimensional latent space).
#' @seealso \code{SingleCellExperiment} \code{reducedDim}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Set latent space
#' latentSpace(exSCE) <- pca(exSCE)$components[, seq_len(10)]
#' @import SingleCellExperiment
#' @docType methods
#' @aliases latentSpace<-,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("latentSpace<-", function(object, value)
  standardGeneric("latentSpace<-"))
setMethod("latentSpace<-", "SingleCellExperiment", function(object, value){
  reducedDim(object, type="CellTrails") <- value
  message("Calculating approximation of CellTrails manifold ",
          "for 2D visualization...")
  tsne_res <- .bhtsne(value)
  X <- tsne_res$Y
  if(!is.null(X)){
    manifold2D(object) <- list(CellTrails.tSNE=X)
    message("Used tSNE perplexity: ", tsne_res$perplexity)
  }
  object})

#' GET CellTrails' latent space
#'
#' Retrieve computed latent space from a \code{SingleCellExperiment} object.
#' @param object A \code{SingleCellExperiment} object
#' @return An object of class \code{matrix}
#' @details Returns the latent space set for a CellTrails analysis. The
#' resulting matrix is numeric. Rows are samples and columns are \emph{d}
#' components. It is a wrapper for \code{reducedDim} to ensure
#' that the proper matrix is received from a \code{SingleCellExperiment}
#' object.
#' @seealso \code{SingleCellExperiment} \code{reducedDim}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Get latent space
#' latentSpace(exSCE)[seq_len(5), ]
#' @docType methods
#' @aliases latentSpace,SingleCellExperiment-method
#' @import SingleCellExperiment
#' @export
#' @author Daniel C. Ellwanger
setGeneric("latentSpace", function(object) standardGeneric("latentSpace"))
setMethod("latentSpace", "SingleCellExperiment", function(object){
  reducedDim(object, type="CellTrails")})

#' SET user-defined landmarks
#'
#' Set user-defined landmarks to a \code{SingleCellExperiment} object.
#' @param object A \code{SingleCellExperiment} object
#' @param value A character vector with sample names
#' @return An updated \code{SingleCellExperiment} object
#' @details Landmarks can be defined on the trajectory and can be necessary to
#' extract individual trails from a trajectory.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the trajectory has not been reconstructed yet.
#' @seealso \code{SingleCellExperiment} \code{fitTrajectory}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Set landmarks
#' userLandmarks(exSCE) <- colnames(exSCE)[5:7]
#' @docType methods
#' @aliases userLandmarks<-,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("userLandmarks<-", function(object, value)
  standardGeneric("userLandmarks<-"))
setMethod("userLandmarks<-", "SingleCellExperiment", function(object, value){
  #Pre-flight check
  if(is.null(.trajGraph(object))) {
    stop("A trajectory has not been computed yet. ",
         "Please, fit the trajectory first.")
  }
  .sampleNameExists(object, value)

  #Run
  #Delete existing user landmarks
  ltypes <- .trajLandmark(object, type="type")
  f <- which(ltypes == "U")
  if(length(f) > 0) {
    .trajLandmark(object[, f], type="id") <- as.factor(NA)
    .trajLandmark(object[, f], type="type") <- as.factor(NA)
    .trajLandmark(object[, f], type="shape") <- "ellipse"
  }
  #Add new user landmarks
  ltypes <- .trajLandmark(object, type="type")
  h_or_b <- which(!is.na(ltypes)) #trail heads and branches
  value <- setdiff(value, colnames(object)[h_or_b]) #keep H or B
  uids <- paste0("U", seq_along(value))
  .trajLandmark(object[, value], type="id") <- uids
  .trajLandmark(object[, value], type="type") <- as.factor("U")
  .trajLandmark(object[, value], type="shape") <- "rectangle"
  object
  })

#' GET user-defined landmarks
#'
#' Gets user-defined landmarks from a \code{SingleCellExperiment} object.
#' @param object A \code{SingleCellExperiment} object
#' @return A character vector with sample names
#' @details Landmarks can be defined on the trajectory by the user
#' with \code{userLandmarks}. Landmarks can be used to extract single
#' trails from a trajectory.
#' @seealso \code{SingleCellExperiment}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Get landmarks
#' userLandmarks(exSCE)
#' @docType methods
#' @aliases userLandmarks,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("userLandmarks", function(object) standardGeneric("userLandmarks"))
setMethod("userLandmarks", "SingleCellExperiment", function(object){
  f <- which(.trajLandmark(object, type="type") == "U")
  ids <- .trajLandmark(object[, f], type="id")
  sNames <- colnames(object[, f])
  names(sNames) <- ids
  sNames})

#' GET landmarks
#'
#' Gets landmarks from a \code{SingleCellExperiment} object.
#' @param object A \code{SingleCellExperiment} object
#' @return A character vector with sample names
#' @details Trail branches (B) and heads (H) are automatically assigned;
#' landmarks can also be defined on the trajectory by the user (U).
#' Landmarks can be used to extract single trails from a trajectory.
#' @seealso \code{userLandmarks}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Get landmarks
#' landmarks(exSCE)[seq_len(5)]
#' @docType methods
#' @aliases landmarks,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("landmarks", function(object) standardGeneric("landmarks"))
setMethod("landmarks", "SingleCellExperiment", function(object){
  .trajLandmark(object, type="id")})

#' ADD trail
#'
#' Function to define a single trail on the trajectory.
#' @param sce An object of class \code{SingleCellExperiment}
#' @param name Name of trail
#' @param from Start landmark
#' @param to End landmark
#' @return An updated object of class \code{SingleCellExperiment}
#' @details A trajectory can be composed of multiple single trails
#' (e.g., developmental progression from a common start towards
#' distinct terminal phenotypes). Start and endpoints of trails can
#' be identified visually using the plot function \code{plotMap}.
#' Here, start (=from) and end (=to) IDs
#' of landmarks are starting with the character "B"
#' (for branching points), "H" (for trail heads, i.e. terminal nodes),
#' and "U" for user-defined landmarks.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the trajectory has not been fitted yet. Please,
#' call \code{fitTrajectory} first. Further, an error is thrown if the
#' provided start or end ID is unknown. A warning is
#' shown if a trail with the same name already exists and gets
#' re-defined.
#' @seealso \code{fitTrajectory} \code{landmarks} \code{plotMap}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Add trail
#' exSCE <- addTrail(exSCE, "H1", "H2", "Tr3")
#' trailNames(exSCE)
#' phenoNames(exSCE)
#' @docType methods
#' @aliases addTrail,SingleCellExperiment-method
#' @importFrom igraph get.shortest.paths distances
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData<-
#' @export
#' @author Daniel C. Ellwanger
setGeneric("addTrail", function(sce, from, to, name)
  standardGeneric("addTrail"))
setMethod("addTrail", "SingleCellExperiment", function(sce, from, to, name){
  from <- toupper(from)
  to <- toupper(to)

  #Pre-flight check
  if(is.null(.trajGraph(sce))) {
    stop("No trajectory information found. Please, compute trajectory first ",
         "(see function 'fitTrajectory').")
  }
  if(!from %in% .trajLandmark(sce, type="id")) {
    stop("Start ID not found.")
  }
  if(!to %in% .trajLandmark(sce, type="id")) {
    stop("End ID not found.")
  }

  ftID <- match(c(from, to), .trajLandmark(sce, type="id")[.useSample(sce)])
  p <- as.vector(get.shortest.paths(.trajGraph(sce),
                                    from = ftID[1],
                                    to = ftID[2])$vpath[[1]])

  smpls <- trajSampleNames(sce)[p]
  ptime <- igraph::distances(.trajGraph(sce), v=ftID[1], to=p)[1, ]
  nm <- paste0("CellTrails.", name)

  if(name %in% trailNames(sce)) { #replace existing trail definition
    warning("A trail with this name already exists and gets replaced.")
    colData(sce)[, nm] <- as.numeric(rep(NA, ncol(sce)))
    colData(sce[, smpls])[, nm] <- ptime
  } else { #new trail definition
    #tc <- sce@int_metadata$CellTrails$trail_cnt
    tnms <- sce@int_metadata$CellTrails$trail_names
    #if(is.null(tnms)) {
      #tc <- 1
      #nm <- name #"tr1"
      #df <- data.frame(cnt=nm, name=name, stringsAsFactors=FALSE)
    #  tnms <- c(name)
    #} else {
      #tc <- tc + 1
      #nm <- name #paste0("tr", tc)
      #df <- data.frame(cnt=nm, name=name, stringsAsFactors=FALSE)
      #df <- rbind(sce@int_metadata$CellTrails$trail_cnt2name, df)
    #  tnms <- c(tnms, name)
    #}
    #sce@int_metadata$CellTrails$trail_cnt <- tc
    #sce@int_metadata$CellTrails$trail_cnt2name <- df
    sce@int_metadata$CellTrails$trail_names <- c(tnms, name)

    colData(sce)[, nm] <- as.numeric(NA) #rep(NA, ncol(sce))
    colData(sce[, smpls])[, nm] <- ptime
  }
  sce})

#' REMOVE trail
#'
#' Removes trail from a \code{SingleCellExperiment} object.
#' @param sce An object of class \code{SingleCellExperiment}
#' @param name Name of trail
#' @return An updated object of class \code{SingleCellExperiment}
#' @details
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the trail name is unknown. All stored trail
#' names can be shown using function \code{trailNames}.
#' @seealso \code{trailNames} \code{addTrail}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Remove trail
#' trailNames(exSCE)
#' exSCE <- removeTrail(exSCE, "Tr1")
#' trailNames(exSCE)
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData<-
#' @docType methods
#' @aliases removeTrail,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("removeTrail", function(sce, name) standardGeneric("removeTrail"))
setMethod("removeTrail", "SingleCellExperiment", function(sce, name){
  #Pre-flight check
  if(!name %in% trailNames(sce)) {
    stop("Could not find a trail with name '", name, "'. Please, check the ",
         "proper spelling of the trail name (case-sensitivity). Valid trail ",
         "names can be received via function 'trailNames'.")
  }
  #f <- which(sce@int_metadata$CellTrails$trail_cnt2name$name == name)
  #nm <- sce@int_metadata$CellTrails$trail_cnt2name$cnt[f]
  #df <- sce@int_metadata$CellTrails$trail_cnt2name[-f, ]
  f <- which(trailNames(sce) == name)
  tnms <- trailNames(sce)[-f]
  sce@int_metadata$CellTrails$trail_names <- tnms
  #sce@int_metadata$CellTrails$trail_cnt2name <- df
  colData(sce)[, paste0("CellTrails.", name)] <- NULL
  sce})

#' GET trail names
#'
#' Function to extract trail names from \code{SingleCellExperiment}
#' object.
#' @param object An object of class \code{SingleCellExperiment}
#' @return A \code{character} vector
#' @seealso \code{addTrail}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' trailNames(exSCE)
#' @docType methods
#' @aliases trailNames,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trailNames", function(object) standardGeneric("trailNames"))
setMethod("trailNames", "SingleCellExperiment", function(object){
  object@int_metadata$CellTrails$trail_names})

#' SET trail names
#'
#' Enables to rename trails stored in a \code{SingleCellExperiment}
#' object.
#' @usage trailNames(object) <- value
#' @param object An object of class \code{SingleCellExperiment}
#' @param value A character vector with the trail names
#' @return An updated object of class \code{SingleCellExperiment}
#' @details
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the number of names does not correspond to the number
#' of trails stored in the object. Further, trail names are required
#' to be unique.
#' @seealso \code{addTrail}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' trailNames(exSCE)
#' trailNames(exSCE) <- c("ABC", "DEF")
#' trailNames(exSCE)
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData<-
#' @docType methods
#' @aliases trailNames<-,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trailNames<-", function(object, value)
  standardGeneric("trailNames<-"))
setMethod("trailNames<-", "SingleCellExperiment", function(object, value) {
  #Pre-flight check
  if(length(trailNames(object)) != length(value)) {
    stop("Number of provided names (", length(value), ") does not correspond ",
         "to number of defined trails (", length(trailNames(object)), ").")
  } else if(any(table(value) > 1)){
    stop("Trail names are required to be unique. Please, choose distinct ",
         "trail names.")
  }
  f <- paste0("CellTrails.", trailNames(object))
  f <- match(f, colnames(colData(object)))
  colnames(colData(object))[f] <- paste0("CellTrails.", value)
  object@int_metadata$CellTrails$trail_names <- value
  object})

#' GET trails
#'
#' Function to extract trail pseudotimes from a
#' \code{SingleCellExperiment} object.
#' @param object An object of class \code{SingleCellExperiment}
#' @return A DataFrame with \code{numeric} columns
#' @seealso \code{addTrail}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' trails(exSCE)
#' @docType methods
#' @aliases trails,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trails", function(object) standardGeneric("trails"))
setMethod("trails", "SingleCellExperiment", function(object){
  if(is.null(trailNames(object))) {
    NULL
  } else {
    f <- paste0("CellTrails.", trailNames(object))
    df <- colData(object)[, f, drop=FALSE]
    colnames(df) <- gsub(x=colnames(df), "CellTrails.", "")
    df
  }})

#' SET 2D manifold representation
#'
#' Stores 2D manifold representation in \code{SingleCellExperiment} object
#' @param object A \code{SingleCellExperiment} object
#' @param value A \code{numeric} matrix with one column per dimension
#' @return An updated object of class \code{SingleCellExperiment}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' gp <- plotManifold(exSCE, color_by="featureName", name="feature_10",
#'                    recalculate=TRUE)
#' manifold2D(exSCE) <- gp
#' @docType methods
#' @aliases manifold2D<-,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("manifold2D<-", function(object, value)
  standardGeneric("manifold2D<-"))
setMethod("manifold2D<-", "SingleCellExperiment", function(object, value){
  #.tsneParams(object) <- c(value$tsne$seed, value$tsne$perplexity)
  if(is.null(value$CellTrails.tSNE)) {
    stop("Plot does not contain proper layout information.")
  }
  #object@int_metadata$CellTrails$manifold2D <- value$CellTrails.tsne
  reducedDim(object, "CellTrails.tSNE") <- value$CellTrails.tSNE
  object})

#' GET 2D manifold representation
#'
#' Returns 2D manifold representation of latent space from
#' \code{SingleCellExperiment} object
#' @param object A \code{SingleCellExperiment} object
#' @return A \code{numeric} matrix
#' @examples
#' # Example data
#' data(exSCE)
#'
#' manifold2D(exSCE)[seq_len(5), ]
#' @docType methods
#' @aliases manifold2D,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("manifold2D", function(object)
  standardGeneric("manifold2D"))
setMethod("manifold2D", "SingleCellExperiment", function(object){
  #object@int_metadata$CellTrails$manifold2D
  reducedDim(object, "CellTrails.tSNE")
})

#' SET state trajectory layout
#'
#' Stores layout of state trajectory in \code{SingleCellExperiment} object
#' @param object A \code{SingleCellExperiment} object
#' @param value A \code{ggplot} object
#' @return An updated object of class \code{SingleCellExperiment}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' gp <- plotStateTrajectory(exSCE, color_by="featureName",
#'                           name="feature_10", recalculate=TRUE)
#' stateTrajLayout(exSCE) <- gp
#' @docType methods
#' @aliases stateTrajLayout<-,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("stateTrajLayout<-", function(object, value)
  standardGeneric("stateTrajLayout<-"))
setMethod("stateTrajLayout<-", "SingleCellExperiment", function(object, value){
  component <- value$CellTrails.stateTrajComponent
  lyt <- value$CellTrails.stateTrajLayout

  if(is.null(component)) {
    stop("No proper layout information detected.")
  }
  stl <- object@int_metadata$CellTrails$stateTrajLayout
  if(is.null(stl) | length(stl) < length(.spanForest(object))) {
    l <- length(.spanForest(object))
    object@int_metadata$CellTrails$stateTrajLayout <- vector("list", l)
  }
  object@int_metadata$CellTrails$stateTrajLayout[[component]] <- lyt
  object})

#' GET state trajectory layout
#'
#' Gets layout of state trajectory from
#' \code{SingleCellExperiment} object
#' @param object A \code{SingleCellExperiment} object
#' @return A \code{numeric} matrix
#' @examples
#' # Example data
#' data(exSCE)
#'
#' stateTrajLayout(exSCE)
#' @docType methods
#' @aliases stateTrajLayout,SingleCellExperiment-method
#' @keywords internal
#' @author Daniel C. Ellwanger
setGeneric(".stateTrajLayout", function(object, component)
  standardGeneric(".stateTrajLayout"))
setMethod(".stateTrajLayout", "SingleCellExperiment", function(object,
                                                              component){
  if(missing(component)) {
    stop("Please, select a component.")
  }
  l <- length(object@int_metadata$CellTrails$stateTrajLayout)
  if(component > l) {
    NULL
  }
  object@int_metadata$CellTrails$stateTrajLayout[[component]]
})

#' SET trajectory layout
#'
#' Sets layout used for trajectory visualization to a
#' \code{SingleCellExperiment} object.
#' @usage trajLayout(object, adjust) <- value
#' @param object An object of class \code{SingleCellExperiment}
#' @param value A data.frame with x- and y-coordinates for
#' each sample (rows = samples, columns = coordinates)
#' @param adjust Indicates if layout has to be adjusted such that edge lengths
#' correlate to pseudotime (default: TRUE)
#' @return An updated object of class \code{SingleCellExperiment}
#' @details
#' CellTrails implements a module which can incorporate pseudotime information
#' into the the graph layout (activated via parameter \code{adjust}). Here,
#' edge lengths between two nodes (samples)
#' will then correspond to the inferred pseudotime that separates two samples
#' along the trajectory.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the number of rows of the layout does not correspond
#' to the number of trajectory samples or if the number of columns is
#' less than 2, or if the row names do not correspond to \code{sampleNames}.
#' @seealso \code{write.ygraphml} \code{trajSampleNames}
#' @examples
#' # Example data
#' data(exSCE)
#' tl <- trajLayout(exSCE)
#'
#' trajLayout(exSCE) <- tl
#' @docType methods
#' @aliases trajLayout<-,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trajLayout<-", function(object, adjust=TRUE, value)
  standardGeneric("trajLayout<-"))
setMethod("trajLayout<-", "SingleCellExperiment", function(object,
                                                           adjust, value) {
  #Pre-flight check
  if(is.null(.trajGraph(object)) & adjust) {
    stop("Please, fit the trajectory first (see 'fitTrajectory').")
  }
  d <- dim(value)
  s <- sum(.useSample(object))
  if(d[1] != s) {
    stop("Number of rows in layout (m=", d[1], ") does not correspond to ",
         "number of trajectory samples (m=", s, ").")
  } else if(d[2] < 2) {
    stop("Number of columns in layout need to be at least 2 ",
         "(numeric columns of x- and y-coordinates per trajectory sample).")
  } else if(all(!rownames(value) %in% trajSampleNames(object))) {
    stop("Rownames of layout do not correspond to trajectory sample names.")
  }
  colvars <- apply(value[, seq_len(2)], 2L, var)
  if(sum(colvars)  == 0) {
    stop("All data points have same coordiates.")
  } else if(colvars[1] == 0) { #make diagonal for linear trajectory
    value[, 1] <- value[, 2]
  } else if(colvars[2] == 0) { #make diagonal for linear trajectory
    value[, 2] <- value[, 1]
  }

  # Run
  X <- value[, seq_len(2)]

  if(adjust){
    X <- .adjustLayoutByPtime(object, X)
  }
  snames <- rownames(value)

  # Store coordinates
  rownames(object@int_colData) <- colnames(object)
  object@int_colData$CellTrails.lytX1 <- as.numeric(rep(NA, ncol(object)))
  object@int_colData$CellTrails.lytX2 <- as.numeric(rep(NA, ncol(object)))
  object[, snames]@int_colData$CellTrails.lytX1 <- X[, 1]
  object[, snames]@int_colData$CellTrails.lytX2 <- X[, 2]

  # Store metadata
  if(ncol(value) == 2) {
    value$shape <- "ellipse"
  }
  .trajLandmark(object[, snames], type="shape") <- value$shape
  userLandmarks <- NULL
  f <- !value$shape == "ellipse"
  if(any(f)) {
    snames <- snames[f]
    userLandmarks(object) <- snames
  }
  object})

#' GET trajectory layout
#'
#' Returns trajectory layout from
#' \code{SingleCellExperiment} object
#' @param object A \code{SingleCellExperiment} object
#' @return A \code{data.frame}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' trajLayout(exSCE)[seq_len(5), ]
#' @docType methods
#' @aliases trajLayout,SingleCellExperiment-method
#' @importFrom BiocGenerics as.data.frame
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trajLayout", function(object)
  standardGeneric("trajLayout"))
setMethod("trajLayout", "SingleCellExperiment", function(object){
  if(is.null(object@int_colData$CellTrails.lytX1) |
     is.null(object@int_colData$CellTrails.lytX2)) {
    NULL
  } else {
    df <- object@int_colData[, c("CellTrails.lytX1", "CellTrails.lytX2")]
    df <- as.data.frame(df)
    colnames(df) <- c("D1", "D2")
    df
  }
})

#' GET trajectory component states
#'
#' Returns states of trajectory components
#' \code{SingleCellExperiment} object
#' @param object A \code{SingleCellExperiment} object
#' @return A \code{character} vector
#' @examples
#' # Example data
#' data(exSCE)
#'
#' trajComponents(exSCE)
#' @docType methods
#' @aliases trajComponents,SingleCellExperiment-method
#' @importFrom igraph V
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trajComponents", function(object)
  standardGeneric("trajComponents"))
setMethod("trajComponents", "SingleCellExperiment", function(object){
  comps <- .spanForest(object)
  if(is.null(comps)){
    NULL
  } else {
    lapply(comps, function(x) V(x)$name)
  }
})
