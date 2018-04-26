#' @include AllClasses.R
NULL

###############################################################################
### CellTrailsSet
###############################################################################
#' SUBSET assay data
#'
#' Extract parts of expression matrix
#' @param x A \code{\link[CellTrails]{CellTrailsSet}} object
#' @param i Feature indices specifying elements to extract
#' @param j Sample indices specifying elements to extract
#' @param drop If TRUE the result is coerced to the lowest possible dimension (default: "missing").
#' @return numerical matrix
#' See \code{\link[base]{drop}} for further details.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Indexing
#' ctset[1:10, 1:10]
#' @export
#' @author Daniel C. Ellwanger
setMethod("[", "CellTrailsSet", function(x, i, j, drop="missing"){
  exprs(x)[i, j, drop = drop]})

#' GET states
#'
#' Retrieve computed states from a \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @return A \code{factor} vector
#' @details State information is extracted from \code{featureData}; factor levels
#' are alphanumerically ordered by ID.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Get state assignments
#' states(ctset)
#' @docType methods
#' @aliases states,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("states", function(object) standardGeneric("states"))
setMethod("states", "CellTrailsSet", function(object){
    object$STATE})
    #s <- unique(object$STATE)
    #o <- order(as.numeric(substring(s, 2)))
    #factor(object$STATE, levels = s[o])})

#' GET trajectory fit
#'
#' Retrieve trajectory fit information from a \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @return An object of type \code{\link[base]{list}} with the following components
#' \describe{
#'   \item{\code{error}}{numerical vector; the actual distance of each sample from the fitted trajectory}
#'   \item{\code{traj}}{\code{\link[igraph]{igraph}} object; graph connecting samples by
#'   their chronological ordering}
#'   \item{\code{blaze}}{data.frame; indicates landmarks (branching points and end points)
#'   in the trajectory graph}
#' }
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Connect states
#' ctset <- connectStates(ctset, l=20)
#'
#' # Select trajectory
#' ctset <- selectTrajectory(ctset, component=1)
#'
#' # Align samples to trajectory
#' ctset <- fitTrajectory(ctset)
#'
#' # Get trajectory fit
#' trajectoryFit(ctset)
#' @docType methods
#' @aliases trajectoryFit,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trajectoryFit", function(object) standardGeneric("trajectoryFit"))
setMethod("trajectoryFit", "CellTrailsSet", function(object){
  object@trajectory})

#' GET latentSpace
#'
#' Retrieve computed latent space from a \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @return An object of class \code{\link[base]{matrix}}
#' @details Returns the latent space computed by spectral embedding of samples. The
#' corresponding matrix is numeric. Rows are samples and columns are \emph{d}
#' components.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Get latent space
#' latentSpace(ctset)
#' @docType methods
#' @aliases latentSpace,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("latentSpace", function(object) standardGeneric("latentSpace"))
setMethod("latentSpace", "CellTrailsSet", function(object){
    object@eigen$space})

#' SET latent space
#'
#' Set computed latent space to a \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}
#' @param value A numeric matrix
#' @return An updated object of class \code{\link[CellTrails]{CellTrailsSet}}
#' @details Rows need to be samples and columns to be \emph{d} components (spanning the
#' lower-dimensional latent space). Usually, this function is not
#' directly accessed by the user.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Perform PCA and set components to object
#' pca_result <- pca(ctset)
#' latentSpace(ctset) <- pca_result$princomp
#' @docType methods
#' @aliases latentSpace<-,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("latentSpace<-", function(object, value) standardGeneric("latentSpace<-"))
setMethod("latentSpace<-", "CellTrailsSet", function(object, value){
  #if(is.null(rownames(value))) {
  #  rownames(value) <- sampleNames(object)
  #}
  rownames(value) <- NULL
  object@eigen$space <- value
  object})

#' GET eigenvalues
#'
#' Retrieve computed eigenvalues from a \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}
#' @return An object of class \code{\link{numeric}}
#' @details Eigenvalues represent the information content for each eigenvector derived by
#' spectral embedding of samples. This function retrieves the derived eigenvalues from a
#' \code{CellTrailsSet} object. The corresponding vector is numeric.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Get eigenvalues
#' eigenvalues(ctset)
#' @docType methods
#' @aliases eigenvalues,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("eigenvalues", function(object) standardGeneric("eigenvalues"))
setMethod("eigenvalues", "CellTrailsSet", function(object){
  object@eigen$values})

#' SET eigenvalues
#'
#' Sets eigenvalues to a \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}
#' @param value A numeric vector
#' @return An updated object of class \code{\link[CellTrails]{CellTrailsSet}}
#' @details Eigenvalues represent the information content for each eigenvector derived by
#' spectral embedding of samples. Vector needs to be numeric and its length needs to correspond
#' to the number of eigenvectors. Usually, this function is not directly accessed by the user.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Perform PCA and set results to object
#' pca_result <- pca(ctset)
#' latentSpace(ctset) <- pca_result$princomp
#' eigenvalues(ctset) <- pca_result$variance
#' @docType methods
#' @aliases eigenvalues<-,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("eigenvalues<-", function(object, value) standardGeneric("eigenvalues<-"))
setMethod("eigenvalues<-", "CellTrailsSet", function(object, value){
  object@eigen$values <- value
  object})

#' GET state trajectory graph
#'
#' Retrieve computed state trajectory graph from a \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @return An object of class \code{\link{list}}
#' @details A graph is computed by spanning states maximizing the overall interface
#' score (i.e., number of neighboring samples between states). This method retrieves
#' the maximum interface spanning forest information from a \code{CellTrailsSet}
#' object. The return value is a \code{list} with an \code{\link{igraph}} object for each
#' component.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Connect states
#' ctset <- connectStates(ctset, l=20)
#'
#' # Get state trajectory graph
#' stateTrajectoryGraph(ctset)
#' @docType methods
#' @aliases stateTrajectoryGraph,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("stateTrajectoryGraph", function(object) standardGeneric("stateTrajectoryGraph"))
setMethod("stateTrajectoryGraph", "CellTrailsSet", function(object){
  object@spanForest})

#' GET trajectory features
#'
#' Retrieve names of features that were selected for trajectory reconstruction from a
#' \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @return An object of class \code{\link{character}}
#' @details Features can be selected prior to trajectory inference. This method retrieves
#' the user-defined features from a \code{CellTrailsSet} object. The return value
#' is a character vector containing the feature names.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Get trajectory features
#' trajectoryFeatures(ctset)
#' @docType methods
#' @aliases trajectoryFeatures,CellTrailsSet-method
#' @import Biobase
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trajectoryFeatures", function(object) standardGeneric("trajectoryFeatures"))
setMethod("trajectoryFeatures", "CellTrailsSet", function(object){
  featureNames(object)[object@useFeature]})

#' SET trajectory features
#'
#' Sets features that are used for trajectory reconstruction to a \code{CellTrailsSet}
#' object.
#' @usage trajectoryFeatures(object) <- value
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @param value A character vector containing the names of the features.
#' @return An updated object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @details Features can be selected prior to trajectory inference. This step is recommended
#' as it increases the accuracy and minimzes the runtime of the trajectory reconstruction
#' algorithm.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' The method states a warning if feature names cannot be found in the stored assay. Since
#' \code{CellTrailsSet} extends class \code{ExpressionSet}, all feature names stored in a
#' \code{CellTrailsSet} object can be retrieved by the function \code{\link[Biobase]{featureNames}}.
#' Further, this method shows a warning if the feature selection generates samples with zero
#' entropy, that is the selected features were not detected in these samples.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Set trajectory features
#' trajectoryFeatures(ctset) <- Biobase::featureNames(ctset)[1:10]
#' @docType methods
#' @aliases trajectoryFeatures<-,CellTrailsSet-method
#' @import Biobase
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trajectoryFeatures<-", function(object, value) standardGeneric("trajectoryFeatures<-"))
setReplaceMethod("trajectoryFeatures", "CellTrailsSet", function(object, value) {
  object@useFeature <- featureNames(object) %in% value
  #Check feature selection
  if(sum(object@useFeature) < length(value)) {
    warning("Some selected features are not listed in expression data, e.g. ",
            value[!value %in% featureNames(object)][1])
  }
  #Check for zero entropy cells
  ze <- sum(apply(object[object@useFeature, ], 2, sum) == 0)
  if(ze > 0) {
    warning(ze, paste0(" sample(s) contain(s) no trajectory information ",
                       "(i.e. have zero entropy) ",
                       "based on selected features. ",
                       "Please, adjust your feature selection."))
  }
  object})

#' GET trajectory states
#'
#' Retrieve state IDs along a selected trajectory from a \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @return A \code{\link[base]{character}} vector
#' @details A state trajectory graph may contain
#' multiple components composed of a set of states. This method retrieves
#' the names of the states contained in the user-selected state trajectory graph
#' component from a \code{CellTrailsSet} object.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Connect states
#' ctset <- connectStates(ctset, l=20)
#'
#' # Select trajectory
#' ctset <- selectTrajectory(ctset, component=1)
#'
#' # Get trajectory states
#' trajectoryStates(ctset)
#' @docType methods
#' @aliases trajectoryStates,CellTrailsSet-method
#' @importFrom igraph V
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trajectoryStates", function(object) standardGeneric("trajectoryStates"))
setMethod("trajectoryStates", "CellTrailsSet", function(object){
  if(length(stateTrajectoryGraph(object)) == 1) {
    V(stateTrajectoryGraph(object)[[1]])$name
  } else {
    NULL
  }})

#' GET trajectory samples
#'
#' Retrieve names of samples that are contained in the trajectory from a \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @return An object of class \code{\link[base]{character}}
#' @details Samples are selected prior to trajectory inference; a state trajectory graph may contain
#' multiple components. Each component is composed of a set of states. This method retrieves
#' the names of the samples that are members of the states contained in the
#' user-selected state trajectory graph component from a \code{CellTrailsSet} object.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Connect states
#' ctset <- connectStates(ctset, l=20)
#'
#' # Select trajectory
#' ctset <- selectTrajectory(ctset, component=1)
#'
#' # Get trajectory samples
#' trajectorySamples(ctset)
#' @docType methods
#' @aliases trajectorySamples,CellTrailsSet-method
#' @import Biobase
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trajectorySamples", function(object) standardGeneric("trajectorySamples"))
setMethod("trajectorySamples", "CellTrailsSet", function(object){
  sampleNames(object)[object@useSample]})

#' GET trajectory layout
#'
#' Retrieve layout for trajectory visualization from
#' \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @return A numeric matrix
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Connect states
#' ctset <- connectStates(ctset, l=20)
#'
#' # Select trajectory
#' ctset <- selectTrajectory(ctset, component=1)
#'
#' # Align samples to trajectory
#' ctset <- fitTrajectory(ctset)
#'
#' # For illustration purposes a layout for each example dataset
#' # is contained in this package.
#' # The function trajectoryLayout allows to set a trajectory layout
#' # to the object; the parameter adjust='TRUE' adjusts the layout
#' # according to the computed pseudotime
#' trajectoryLayout(ctset, adjust=TRUE) <- trajLayouts$example
#'
#' # Get trajectory layout
#' trajectoryLayout(ctset)
#' @docType methods
#' @aliases trajectoryLayout,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trajectoryLayout", function(object) standardGeneric("trajectoryLayout"))
setMethod("trajectoryLayout", "CellTrailsSet", function(object){
  object@layout})

#' SET trajectory layout
#'
#' Sets layout which is used for trajectory visualization to a \code{CellTrailsSet}
#' object.
#' @usage trajectoryLayout(object, adjust) <- value
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @param value A matrix with x- and y-coordinates for each sample; rows =
#' samples, columns = coordinates.
#' @param adjust If adjustment of the layout is required. (default: TRUE)
#' @return An updated object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @details
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the number of rows of the layout does not correspond to
#' the number of \code{trajectorySamples} or if the number of columns does not
#' equal 2 or if the row names do not correspond to \code{sampleNames}.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Connect states
#' ctset <- connectStates(ctset, l=20)
#'
#' # Select trajectory
#' ctset <- selectTrajectory(ctset, component=1)
#'
#' # Align samples to trajectory
#' ctset <- fitTrajectory(ctset)
#'
#' # For illustration purposes a layout for each example dataset
#' # is contained in this package.
#' # The function trajectoryLayout allows to set a trajectory layout
#' # to the object; the parameter adjust='TRUE' adjusts the layout
#' # according to the computed pseudotime
#' trajectoryLayout(ctset, adjust=TRUE) <- trajLayouts$example
#'
#' @docType methods
#' @aliases trajectoryLayout<-,CellTrailsSet-method
#' @import Biobase
#' @export
#' @seealso \code{\link[CellTrails]{trajectorySamples}}
#' @author Daniel C. Ellwanger
setGeneric("trajectoryLayout<-", function(object, adjust=FALSE, value) standardGeneric("trajectoryLayout<-"))
setReplaceMethod("trajectoryLayout", "CellTrailsSet", function(object, adjust, value) {
  if(is.null(value)) {
    object@layout <- value
  } else {
    #Pre-flight check
    d <- dim(value)
    s <- sum(object@useSample)
    if(d[1] != s) {
      stop("Number of rows in layout (m=", d[1], ") does not correspond to ",
           "number of trajectory samples (m=", s, ").")
    }
    if(d[2] != 2) {
      stop("Number of columns in layout need to be 2 (numeric matrix of x- ",
           "and y-coordinates per trajectory sample).")
    }
    if(all(!rownames(value) %in% sampleNames(object)[object@useSample])) {
      stop("Rownames of layout do not correspond to trajectory sample names.")
    }
    vars <- apply(value, 2L, var)
    if(sum(vars)  == 0) {
      stop("All data points have same coordiates.")
    } else if(vars[1] == 0) { #make diagonal for linear trajectory
      value[, 1] <- value[, 2]
    } else if(vars[2] == 0) {
      value[, 2] <- value[, 1]
    }

    colnames(value) <- c("D1", "D2")
    object@layout <- value #data.frame("D1" = value[,1], "D2" = value[,2])

    if(adjust){
      value_new <- .adjustLayoutByPtime(object)
      object@layout <- value_new #data.frame("D1" = value[,1], "D2" = value[,2])
    }
  }
  #rownames(object@layout) <- seq(nrow(value))
  object})

#' ADD trail
#'
#' Function to define a single trail on the trajectory.
#' @param ctset An object of class \code{\link[CellTrails]{CellTrailsSet}}
#' @param name Name of trail
#' @param from Start node ID
#' @param to End node ID
#' @return An updated object of class \code{\link[CellTrails]{CellTrailsSet}}
#' @details A trajectory can be composed of multiple single trails (e.g., developmental
#' progression from a common start towards distinct terminal phenotypes).
#' Start and endpoints of trails can be identified using the plot function
#' (\code{plot(ctset, type = "trailblazing")}. Here, start (=from) and end (=to) IDs
#' of samples are starting with the character "B"
#' (for branching points) or "H" (for trail heads, i.e. terminal nodes).
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the trajectory has not been fitted yet. Please,
#' call \code{\link{fitTrajectory}} first. Further, an error is thrown if the
#' provided start or end ID is unknown.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Connect states
#' ctset <- connectStates(ctset, l=20)
#'
#' # Select trajectory
#' ctset <- selectTrajectory(ctset, component=1)
#'
#' # Align samples to trajectory
#' ctset <- fitTrajectory(ctset)
#'
#' # For illustration purposes a layout for each example dataset
#' # is contained in this package.
#' # The function trajectoryLayout allows to set a trajectory layout
#' # to the object; the parameter adjust='TRUE' adjusts the layout
#' # according to the computed pseudotime
#' trajectoryLayout(ctset, adjust=TRUE) <- trajLayouts$example
#'
#' # Define trail
#' ctset <- addTrail(ctset, from="H1", to="H3", name="Tr1") # Define trail
#' @docType methods
#' @aliases addTrail,CellTrailsSet-method
#' @importFrom igraph get.shortest.paths
#' @export
#' @author Daniel C. Ellwanger
setGeneric("addTrail", function(ctset, from, to, name) standardGeneric("addTrail"))
setMethod("addTrail", "CellTrailsSet", function(ctset, from, to, name){
  from <- toupper(from)
  to <- toupper(to)

  #Pre-flight check
  if(is.null(ctset@trajectory)) {
    stop("No trajectory information found. Please, compute trajectory first ",
         "(see function 'fitTrajectory').")
  }
  if(!from %in% ctset@trajectory$blaze$id) {
    stop("Start ID not found.")
  }
  if(!to %in% ctset@trajectory$blaze$id) {
    stop("End ID not found.")
  }

  ftID <- match(c(from, to), ctset@trajectory$blaze$id)
  p <- as.vector(get.shortest.paths(ctset@trajectory$traj,
                                    from = ftID[1],
                                    to = ftID[2])$vpath[[1]])

  ctset@trails[[name]] <- list()
  ctset@trails[[name]]$path <- p
  ctset@trails[[name]]$samples <- which(ctset@useSample)[p]
  ctset@trails[[name]]$ptime <- distances(ctset@trajectory$traj, v=ftID[1], to=p)[1, ]
  #ctset@trails[[name]]$ptime <- ctset@trails[[name]]$ptime / max(ctset@trails[[name]]$ptime)
  ctset})

#' REMOVE trail
#'
#' Removes trail from a \code{CellTrailsSet} object.
#' @param ctset An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @param name Name of trail
#' @return An updated object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @details
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the trail name is unknown. All stored trail names can be shown
#' using function \code{\link[CellTrails]{trailNames}}.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Connect states
#' ctset <- connectStates(ctset, l=20)
#'
#' # Select trajectory
#' ctset <- selectTrajectory(ctset, component=1)
#'
#' # Align samples to trajectory
#' ctset <- fitTrajectory(ctset)
#'
#' # For illustration purposes a layout for each example dataset
#' # is contained in this package.
#' # The function trajectoryLayout allows to set a trajectory layout
#' # to the object; the parameter adjust='TRUE' adjusts the layout
#' # according to the computed pseudotime
#' trajectoryLayout(ctset, adjust=TRUE) <- trajLayouts$example
#'
#' # Define trail
#' ctset <- addTrail(ctset, from="H1", to="H3", name="Tr1") # Define trail
#' ctset <- addTrail(ctset, from="H1", to="H4", name="Tr2") # Define trail
#' trailNames(ctset)
#'
#' # Remove trail
#' ctset <- removeTrail(ctset, name="Tr1")
#' trailNames(ctset)
#' @docType methods
#' @aliases removeTrail,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("removeTrail", function(ctset, name) standardGeneric("removeTrail"))
setMethod("removeTrail", "CellTrailsSet", function(ctset, name){
  #Pre-flight check
  if(!name %in% trailNames(ctset)) {
    stop("Could not find a trail with name '", name, "'.")
  }
  ctset@trails[[name]] <- NULL
  ctset})

#' GET trail names
#'
#' Get names of trails from \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @return A \code{character} vector
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Connect states
#' ctset <- connectStates(ctset, l=20)
#'
#' # Select trajectory
#' ctset <- selectTrajectory(ctset, component=1)
#'
#' # Align samples to trajectory
#' ctset <- fitTrajectory(ctset)
#'
#' # For illustration purposes a layout for each example dataset
#' # is contained in this package.
#' # The function trajectoryLayout allows to set a trajectory layout
#' # to the object; the parameter adjust='TRUE' adjusts the layout
#' # according to the computed pseudotime
#' trajectoryLayout(ctset, adjust=TRUE) <- trajLayouts$example
#'
#' # Define trail
#' ctset <- addTrail(ctset, from="H1", to="H3", name="Tr1") # Define trail
#' ctset <- addTrail(ctset, from="H1", to="H4", name="Tr2") # Define trail
#'
#' # Get trail names
#' trailNames(ctset)
#' @docType methods
#' @aliases trailNames,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trailNames", function(object) standardGeneric("trailNames"))
setMethod("trailNames", "CellTrailsSet", function(object){
  names(object@trails)})

#' SET trail names
#'
#' Enables to rename trails stored in a \code{CellTrailsSet}
#' object.
#' @usage trailNames(object) <- value
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}
#' @param value A character vector with the trail names
#' @return An updated object of class \code{\link[CellTrails]{CellTrailsSet}}
#' @details
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the number of names does not correspond to the number of
#' trails stored in the object.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Connect states
#' ctset <- connectStates(ctset, l=20)
#'
#' # Select trajectory
#' ctset <- selectTrajectory(ctset, component=1)
#'
#' # Align samples to trajectory
#' ctset <- fitTrajectory(ctset)
#'
#' # For illustration purposes a layout for each example dataset
#' # is contained in this package.
#' # The function trajectoryLayout allows to set a trajectory layout
#' # to the object; the parameter adjust='TRUE' adjusts the layout
#' # according to the computed pseudotime
#' trajectoryLayout(ctset, adjust=TRUE) <- trajLayouts$example
#'
#' # Define trail
#' ctset <- addTrail(ctset, from="H1", to="H3", name="Tr1") # Define trail
#' ctset <- addTrail(ctset, from="H1", to="H4", name="Tr2") # Define trail
#' trailNames(ctset)
#'
#' # Update trail names
#' trailNames(ctset) <- c("TrA", "TrB")
#' trailNames(ctset)
#' @docType methods
#' @aliases trailNames<-,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trailNames<-", function(object, value) standardGeneric("trailNames<-"))
setReplaceMethod("trailNames", "CellTrailsSet", function(object, value) {
  #Pre-flight check
  if(length(object@trails) != length(value)) {
    stop("Number of provided names (", length(value), ") does not correspond to ",
         "number of defined trails (", length(object@trails), ").")
  }

  names(object@trails) <- value
  object})

#' GET trail samples
#'
#' Get names of samples for each trail from \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @param trail_name A character string indicating the name of the trail (optional)
#' @return An object of class \code{list}; each element is a vector of sample
#' names per trail
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Connect states
#' ctset <- connectStates(ctset, l=20)
#'
#' # Select trajectory
#' ctset <- selectTrajectory(ctset, component=1)
#'
#' # Align samples to trajectory
#' ctset <- fitTrajectory(ctset)
#'
#' # For illustration purposes a layout for each example dataset
#' # is contained in this package.
#' # The function trajectoryLayout allows to set a trajectory layout
#' # to the object; the parameter adjust='TRUE' adjusts the layout
#' # according to the computed pseudotime
#' trajectoryLayout(ctset, adjust=TRUE) <- trajLayouts$example
#'
#' # Define trail
#' ctset <- addTrail(ctset, from="H1", to="H3", name="Tr1") # Define trail
#'
#' # Get trail samples
#' trailSamples(ctset, trail_name="Tr1")
#' @docType methods
#' @aliases trailSamples,CellTrailsSet-method
#' @import Biobase
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trailSamples", function(object, trail_name=NULL) standardGeneric("trailSamples"))
setMethod("trailSamples", "CellTrailsSet", function(object, trail_name=NULL){
  if(is.null(trail_name)) {
    lapply(object@trails, function(x) sampleNames(object)[x$samples])
  } else {
    if(trail_name %in% trailNames(object)) {
      sampleNames(object)[object@trails[[trail_name]]$samples]
    } else {
      stop("Please check name of trail. A trail named ", trail_name, " was not defined.")
    }
  }})

#' GET trail pseudotime
#'
#' Get pseudotime of each trail from \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @param trail_name A character string indicating the name of the trail (optional)
#' @return An object of class \code{list}; each element is a vector of pseudotime
#' values per trail
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Connect states
#' ctset <- connectStates(ctset, l=20)
#'
#' # Select trajectory
#' ctset <- selectTrajectory(ctset, component=1)
#'
#' # Align samples to trajectory
#' ctset <- fitTrajectory(ctset)
#'
#' # For illustration purposes a layout for each example dataset
#' # is contained in this package.
#' # The function trajectoryLayout allows to set a trajectory layout
#' # to the object; the parameter adjust='TRUE' adjusts the layout
#' # according to the computed pseudotime
#' trajectoryLayout(ctset, adjust=TRUE) <- trajLayouts$example
#'
#' # Define trail
#' ctset <- addTrail(ctset, from="H1", to="H3", name="Tr1") # Define trail
#'
#' # Get trail pseudotime per sample
#' trailPseudotime(ctset, trail_name="Tr1")
#' @docType methods
#' @aliases trailPseudotime,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trailPseudotime", function(object, trail_name=NULL) standardGeneric("trailPseudotime"))
setMethod("trailPseudotime", "CellTrailsSet", function(object, trail_name=NULL){
  if(is.null(trail_name)) {
    lapply(object@trails, function(x) x$ptime / max(x$ptime))
  } else {
    if(trail_name %in% trailNames(object)) {
      object@trails[[trail_name]]$ptime / max(object@trails[[trail_name]]$ptime)
    } else {
      stop("Please check name of trail. A trail named ", trail_name, " was not defined.")
    }
  }})

#' GET trail states
#'
#' Get states of each trail from \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @param trail_name A character string indicating the name of the trail (optional)
#' @return An object of class \code{list}; each element is a vector of states
#' factor values per trail
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Connect states
#' ctset <- connectStates(ctset, l=20)
#'
#' # Select trajectory
#' ctset <- selectTrajectory(ctset, component=1)
#'
#' # Align samples to trajectory
#' ctset <- fitTrajectory(ctset)
#'
#' # For illustration purposes a layout for each example dataset
#' # is contained in this package.
#' # The function trajectoryLayout allows to set a trajectory layout
#' # to the object; the parameter adjust='TRUE' adjusts the layout
#' # according to the computed pseudotime
#' trajectoryLayout(ctset, adjust=TRUE) <- trajLayouts$example
#'
#' # Define trail
#' ctset <- addTrail(ctset, from="H1", to="H3", name="Tr1") # Define trail
#'
#' # Get trail state per sample
#' trailStates(ctset, trail_name="Tr1")
#' @docType methods
#' @aliases trailStates,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("trailStates", function(object, trail_name=NULL) standardGeneric("trailStates"))
setMethod("trailStates", "CellTrailsSet", function(object, trail_name=NULL){
  if(is.null(trail_name)) {
    lapply(object@trails, function(x) states(object)[x$samples])
  } else {
    if(trail_name %in% trailNames(object)) {
      states(object)[object@trails[[trail_name]]$samples]
    } else {
      stop("Please check name of trail. A trail named ", trail_name, " was not defined.")
    }
  }})

#' ADD landmark
#'
#' Adds a landmark to a \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}.
#' @param sample_name A character string indicating the name of the sample
#' @return An updated object of class \code{\link[CellTrails]{CellTrailsSet}}
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Connect states
#' ctset <- connectStates(ctset, l=20)
#'
#' # Select trajectory
#' ctset <- selectTrajectory(ctset, component=1)
#'
#' # Align samples to trajectory
#' ctset <- fitTrajectory(ctset)
#'
#' # For illustration purposes a layout for each example dataset
#' # is contained in this package.
#' # The function trajectoryLayout allows to set a trajectory layout
#' # to the object; the parameter adjust='TRUE' adjusts the layout
#' # according to the computed pseudotime
#' trajectoryLayout(ctset, adjust=TRUE) <- trajLayouts$example
#'
#' # Add user-defined landmark
#' ctset <- addLandmark(ctset, sample_name="C3_sample_10")
#'
#' plot(ctset, type="trailblazing")
#' @docType methods
#' @aliases addLandmark,CellTrailsSet-method
#' @import Biobase
#' @export
#' @author Daniel C. Ellwanger
setGeneric("addLandmark", function(object, sample_name) standardGeneric("addLandmark"))
setMethod("addLandmark", "CellTrailsSet", function(object, sample_name){
  sNames <- sampleNames(object)[object@useSample]

  if(!sample_name %in% sNames) {
    stop("Please check sample name. The provided sample seems not to be used ",
         "for trajectory reconstruction.")
  } else {
    idx <- sNames == sample_name
    object@trajectory$blaze$type[idx] <- "U"
    f <- which(object@trajectory$blaze$type == "U")
    object@trajectory$blaze$id[f] <- paste0("U", seq_len(length(f)))
    object@trajectory$blaze$shape[idx] <- "rectangle"
  }
  object
  })

#' REMOVE landmark
#'
#' Removes a landmark to a \code{CellTrailsSet} object.
#' @param object An object of class \code{\link[CellTrails]{CellTrailsSet}}
#' @param sample_name A character string indicating the name of the sample
#' @return An updated object of class \code{\link[CellTrails]{CellTrailsSet}}
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#'
#' # Find spectrum
#' spectr <- findSpectrum(ctset, frac=25)
#'
#' # Reduce dimensionality
#' ctset <- reduceDimensions(ctset, spectr)
#'
#' # Find states
#' ctset <- findStates(ctset, max_pval=1e-3, min_feat=4)
#'
#' # Connect states
#' ctset <- connectStates(ctset, l=20)
#'
#' # Select trajectory
#' ctset <- selectTrajectory(ctset, component=1)
#'
#' # Align samples to trajectory
#' ctset <- fitTrajectory(ctset)
#'
#' # For illustration purposes a layout for each example dataset
#' # is contained in this package.
#' # The function trajectoryLayout allows to set a trajectory layout
#' # to the object; the parameter adjust='TRUE' adjusts the layout
#' # according to the computed pseudotime
#' trajectoryLayout(ctset, adjust=TRUE) <- trajLayouts$example
#'
#' # Add user-defined landmark
#' ctset <- addLandmark(ctset, sample_name="C3_sample_10")
#'
#' # Remove landmark
#' ctset <- removeLandmark(ctset, sample_name="C3_sample_10")
#' plot(ctset, type="trailblazing")
#' @docType methods
#' @aliases removeLandmark,CellTrailsSet-method
#' @import Biobase
#' @export
#' @author Daniel C. Ellwanger
setGeneric("removeLandmark", function(object, sample_name) standardGeneric("removeLandmark"))
setMethod("removeLandmark", "CellTrailsSet", function(object, sample_name){
  sNames <- sampleNames(object)[object@useSample]

  if(!sample_name %in% sNames) {
    stop("Please check sample name. The provided sample seems not to be used ",
         "for trajectory reconstruction.")
  } else {
    idx <- sNames == sample_name
    object@trajectory$blaze$type[idx] <- NA
    f <- which(object@trajectory$blaze$type == "U")
    object@trajectory$blaze$id[f] <- paste0("U", seq_len(length(f)))
    object@trajectory$blaze$shape[idx] <- "ellipse"
  }
  object
  })
