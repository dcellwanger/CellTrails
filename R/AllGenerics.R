#' @include AllClasses.R
NULL

#' Filter features by Detection Level (DL)
#'
#' Filters trajectory features that are detected in a minimum number of samples.
#' @param ctset An \code{\link[CellTrails]{CellTrailsSet}} object
#' @param threshold Minimum number of samples; if value < 1 it is interpreted as
#' fraction, otherwise as absolute sample count
#' @return An updated \code{\link[CellTrails]{CellTrailsSet}} object
#' @seealso \code{\link{trajectoryFeatures}}
#' @details The detection level denotes the fraction of samples in which a feature was detected.
#' For each trajectory feature listed in the
#' CellTrailsSet object the relative number of samples having a feature expression
#' value greater than 0 is counted. Features that are expressed in a fraction of all
#' samples greater than \code{threshold} remain labeled as trajectory feature in
#' the \code{CellTrailsSet} object, otherwise
#' they are not considered for dimensionality reduction, clustering and trajectory reconstruction.
#' If the parameter \code{threshold} fullfills
#' \code{threshold} \eqn{>= 1} it becomes converted to a relative fraction of the
#' total sample count.
#' @examples
#' # Simulate example data
#' dat <- simulate_exprs(n_features=15000, n_samples=100)
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Filter features
#' ctset <- filterFeaturesByDL(ctset, threshold=2)
#'
#' length(trajectoryFeatures(ctset))
#' @docType methods
#' @aliases filterFeaturesByDL,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("filterFeaturesByDL", function(ctset, threshold) standardGeneric("filterFeaturesByDL"))
setMethod("filterFeaturesByDL", "CellTrailsSet", function(ctset, threshold){
  .filterFeaturesByDL_def(x=ctset, threshold=threshold)
})

#' Filter features by Coefficient of Variation (COV)
#'
#' Filters trajectory features that have a minimal coefficient of variation.
#' @param ctset An \code{\link[CellTrails]{CellTrailsSet}} object
#' @param threshold Minimum coefficient of variation; numeric value between 0 and 1
#' @param design A numeric matrix describing the factors that should be blocked
#' @return An updated \code{\link[CellTrails]{CellTrailsSet}} object
#' @seealso \code{\link{trajectoryFeatures}}
#' @details For each trajectory feature \emph{x} listed in the CellTrailsSet object
#' the coefficient of variation is computed by
#' \eqn{CoV(x) = sd(x) / mean(x)}. Features with a CoV(x) greater than \code{threshold}
#' remain labeled as trajectory feature in the \code{CellTrailsSet} object, otherwise
#' they are not considered for dimensionality reduction, clustering and trajectory reconstruction.
#' \cr \cr
#' To account for systematic bias in the expression data (e.g., cell cycle effects), a design matrix can be
#' provided for the learning process. It should list the factors that should be blocked and their values per
#' sample. It is suggested to construct a design matrix with \code{\link[stats]{model.matrix}}.
#' @examples
#' # Simulate example data
#' dat <- simulate_exprs(n_features=15000, n_samples=100)
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Filter
#' ctset <- filterFeaturesByDL(ctset, threshold=2)
#' ctset <- filterFeaturesByCOV(ctset, threshold=0.5)
#'
#' length(trajectoryFeatures(ctset))
#' @docType methods
#' @aliases filterFeaturesByCOV,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("filterFeaturesByCOV", function(ctset, threshold, design=NULL) standardGeneric("filterFeaturesByCOV"))
setMethod("filterFeaturesByCOV", "CellTrailsSet", function(ctset, threshold, design){
  .filterFeaturesByCOV_def(x=ctset, threshold=threshold, design)
})

#' Filter features by Fano Factor
#'
#' Filters trajectory features that exhibit a significantly high fano factor
#' (index of dispersion) by considering average expression levels.
#' @param ctset An \code{\link[CellTrails]{CellTrailsSet}} object
#' @param threshold A Z-score cutoff (default: 1.7)
#' @param min_expr Minimum average expression of feature to be considered
#' @param design A numeric matrix describing the factors that should be blocked
#' for filter procedure (default: 0)
#' @return An updated \code{\link[CellTrails]{CellTrailsSet}} object
#' @seealso \code{\link{trajectoryFeatures}}
#' @details To identify the most variable features an unsupervised strategy
#' that controls for the relationship between a features’s average expression
#' intensity and its expression variability is applied. Features are placed into
#' 20 bins based on their mean expression. For each bin the fano factor
#' (a windowed version of the index of dispersion, IOD = variance / mean) distribution
#' is computed and standardized (\emph{Z}-score(\emph{x}) = x/sd(\emph{x}) - mean(\emph{x})/sd(\emph{x})).
#' Features with a \emph{Z}-score
#' greater than \code{threshold} remain labeled as trajectory feature
#' in the \code{CellTrailsSet} object. The parameter \code{min_expr} defines the minimum
#' average expression level of a feature to be considered for this filter method.
#' \cr \cr
#' To account for systematic bias in the expression data (e.g., cell cycle effects), a design matrix can be
#' provided for the learning process. It should list the factors that should be blocked and their values per
#' sample. It is suggested to construct a design matrix with \code{\link[stats]{model.matrix}}.
#' @examples
#' # Simulate example data
#' dat <- simulate_exprs(n_features=15000, n_samples=100)
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Filter
#' ctset <- filterFeaturesByDL(ctset, threshold=2)
#' ctset <- filterFeaturesByCOV(ctset, threshold=0.5)
#' ctset <- filterFeaturesByFF(ctset, threshold=1.7, min_expr=0)
#'
#' length(trajectoryFeatures(ctset))
#' @docType methods
#' @aliases filterFeaturesByFF,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("filterFeaturesByFF", function(ctset, threshold=1.7, min_expr=0, design=NULL) standardGeneric("filterFeaturesByFF"))
setMethod("filterFeaturesByFF", "CellTrailsSet", function(ctset, threshold, min_expr, design){
  .filterFeaturesByIOD_def(x=ctset, z=threshold, min_expr=min_expr, design)
})

#' Spectral embedding of biological samples
#'
#' Non-linear learning of a data representation that captures the intrinsic geometry of
#' the trajectory. CellTrails uses spectral decomposition of a graph encoding
#' sample-to-sample similarities (based on their mutual information) to identify the
#' chronological ordering of samples (e.g., cells).
#' @param ctset A \code{\link[CellTrails]{CellTrailsSet}} object
#' @param design A numeric matrix describing the factors that should be blocked
#' @return An updated \code{\link[CellTrails]{CellTrailsSet}} object
#' @details Single-cell gene expression measurements comprise high-dimensional data of large
#' volume, i.e. many features (e.g., genes) are measured in many samples (e.g., cells); or
#' more formally, \emph{m} samples can be described by the expression of \emph{n} features
#' (i.e. \emph{n} dimensions). The cells’ expression profiles are shaped by many distinct
#' unobserved biological causes related to each cell's geno- and phenotype, such as developmental
#' age, tissue region of origin, cell cycle stage, as well as extrinsic sources such as status
#' of signaling receptors, and environmental stressors, but also technical noise.
#' In other words, a single dimension, despite just containing gene expression information,
#' represents an underlying combination of multiple dependent and independent, relevant and
#' non-relevant factors, whereat each factors’ individual contribution is non-uniform.
#' To obtain a better resolution and to extract underlying information, CellTrails aims to find
#' a meaningful low-dimensional structure - a manifold - that represents cells mainly by their
#' temporal relation along a biological process.
#' \cr \cr
#' This method assumes that the expression vectors are lying on or near a manifold with dimensionality
#' \emph{d} that is embedded in the \emph{n}-dimensional space. By using spectral embedding CellTrails
#' aims to amplify latent temporal information; it reduces noise (ie. truncates non-relevant
#' dimensions) by transforming the expression matrix into a new dataset while retaining the geometry
#' of the original dataset as much as possible. CellTrails captures overall cell-to-cell relations
#' based on the statistical mutual dependency between any two data vectors. A high temporal dependency
#' between two samples should be represented by their close proximity in the lower-dimensional space.
#' \cr \cr
#' First, the mutual depencency between samples is scored using mutual information. This entropy framework
#' naturally requires discretization of data vectors by an indicator function, which assigns each
#' continuous data point (expression value) to exactly one discrete interval (e.g. low, mid or high).
#' However, measurement points located close to the interval borders may get wrongly assigned due to
#' noise-induced fluctuations. Therefore, CellTrails fuzzifies the indicator function by using a
#' piecewise polynomial function, i.e. the domain of each sample expression vector is divided into
#' contiguous intervals (based on Daub \emph{et al.}, 2004). Second, the computed mutual
#' information matrix, which is left-bounded and composed of bits, is scaled to a generalized correlation
#' coefficient. Third, CellTrails constructs a
#' simple complete graph with \emph{m} nodes, one for each data vector (ie. sample), and weights each edge
#' between two nodes by a heat kernel function applied on the generalzied correlation coefficient. Finally,
#' nonlinear spectral embedding (ie. spectral decomposition of the graph's adjacency matrix) is performed
#' (Belkin & Niyogi, 2003; Sussman \emph{et al.}, 2012) unfolding the manifold.
#' \cr \cr
#' To account for systematic bias in the expression data (e.g., cell cycle effects), a design matrix can be
#' provided for the learning process. It should list the factors that should be blocked and their values per
#' sample. It is suggested to construct a design matrix with \code{\link[stats]{model.matrix}}.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' The method throws a warning if selected trajectory features generate samples with zero entropy (e.g., the
#' samples exclusively contain non-detects, that is all expression values are zero).
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Embed samples
#' ctset <- embedSamples(ctset)
#' ctset
#' @references Daub, C.O., Steuer, R., Selbig, J., and Kloska, S. (2004).
#' Estimating mutual information using B-spline functions -- an improved similarity
#' measure for analysing gene expression data. BMC Bioinformatics 5, 118.
#' @references Belkin, M., and Niyogi, P. (2003). Laplacian eigenmaps for dimensionality
#' reduction and data representation. Neural computation 15, 1373-1396.
#' @references Sussman, D.L., Tang, M., Fishkind, D.E., and Priebe, C.E. (2012).
#' A Consistent Adjacency Spectral Embedding for Stochastic Blockmodel Graphs.
#' J Am Stat Assoc 107, 1119-1128.
#' @seealso \code{\link[stats]{model.matrix}}
#' @docType methods
#' @aliases embedSamples,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("embedSamples", function(ctset, design=NULL) standardGeneric("embedSamples"))
setMethod("embedSamples", "CellTrailsSet", function(ctset, design){
  .embedSamples_def(x=ctset, design=design)
})

#' Determine number of informative latent dimensions
#'
#' Identifies the dimensionality of the latent space
#' @param ctset A \code{\link[CellTrails]{CellTrailsSet}} object
#' @param frac Fraction or number (if \code{frac > 1}) of eigengaps used to perform
#' linear fit. (default: 100)
#' @return An object of class \code{CellTrailsSpectrum}
#' @details Similar to a scree plot, this method generates a simple line segement
#' plot showing the lagged differences between ordered eigenvalues (eigengaps). A
#' linear fit is calucated on a fraction of top ranked values to identify informative
#' eigenvectors.
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
#' plot(spectr)
#' @docType methods
#' @aliases findSpectrum,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("findSpectrum", function(ctset, frac=100) standardGeneric("findSpectrum"))
setMethod("findSpectrum", "CellTrailsSet", function(ctset, frac){
  .findSpectrum_def(x=ctset, frac=frac)
})

#' Dimensionality reduction
#'
#' Truncates eigenbasis based on the eigengap criterion.
#' @param ctset A \code{\link[CellTrails]{CellTrailsSet}} object
#' @param ctspec A \code{\link[CellTrails]{CellTrailsSpectrum}} object
#' @return An updated object of class \code{CellTrailsSet}
#' @details The number of informative dimensions was determined
#' using the function \code{\link{findSpectrum}}. The resulting
#' object is of class \code{CellTrailsSpectrum} and is used
#' in this method to reduce the dimensionality of the latent space
#' (ie. the result of the sample embedding).
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the samples stored in the \code{CellTrailsSet} object were
#' not embedded yet (ie. the \code{CellTrailsSet} object does not contain a
#' latent space matrix object; function call \code{\link[CellTrails]{latentSpace}}
#' returns NULL; see \code{\link[CellTrails]{embedSamples}}).
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
#' ctset
#'
#' dim(latentSpace(ctset))
#' plot(ctset, type="latentSpace", feature_name="feature_1")
#' @docType methods
#' @aliases reduceDimensions,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("reduceDimensions", function(ctset, ctspec) standardGeneric("reduceDimensions"))
setMethod("reduceDimensions", "CellTrailsSet", function(ctset, ctspec){
  #Pre-flight check
  if(is.null(latentSpace(ctset))) {
    stop("Samples haven't been embedded yet. Please, call function 'embedSamples' first.")
  }
  #Run
  .reduceDimensions_def(x=ctset, s=ctspec)
})

#' Principal Component Analysis
#'
#' Performs principal component analysis by spectral decomposition of
#' a covariance or correlation matrix
#' @param ctset An \code{\link{CellTrailsSet}} object
#' @param do_scaling FALSE = covariance matrix, TRUE = correlation matrix
#' @param design A numeric matrix describing the factors that should be blocked
#' @return A \code{list} object containing the following components:
#' @return \item{\code{princomp}}{Principal components}
#' @return \item{\code{variance}}{Variance explained by each component}
#' @return \item{\code{loadings}}{Loading score for each feature}
#' @details The calculation is done by a spectral decomposition of the
#' (scaled) covariance matrix of the trajectory features
#' as defined in the \code{CellTrails} object.
#' Features with zero variance get automatically removed. To account for systematic bias in the
#' expression data (e.g., cell cycle effects), a design matrix can be
#' provided for the learning process. It should list the factors that should be blocked and
#' their values per sample. It is suggested to construct a design matrix with
#' \code{\link[stats]{model.matrix}}.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' # Create container
#' ctset <- as.CellTrailsSet(dat)
#'
#' # Principal component analysis
#' pca_result <- pca(ctset)
#'
#' barplot(pca_result$variance[1:10], ylab="Variance",
#'         names.arg=colnames(pca_result$princomp)[1:10], las=2)
#' plot(pca_result$princomp, xlab="PC1", ylab="PC2")
#' @docType methods
#' @aliases pca,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("pca", function(ctset, do_scaling=TRUE, design=NULL) standardGeneric("pca"))
setMethod("pca", "CellTrailsSet", function(ctset, do_scaling, design){
  .pca(x=ctset, do_scaling=do_scaling, design=design)
})

#' Identify trajectory states
#'
#' Determines states using hierarchical spectral clustering with a \emph{post-hoc} test.
#' @param ctset A \code{\link[CellTrails]{CellTrailsSet}} object
#' @param min_size The initial cluster dedrogram is cut at an height such that the minimum cluster
#' size is at least \code{min_size}; if \code{min_size} < 1 than the fraction of total samples
#' is used, otherwise it is used as absoulte count (default: 0.01).
#' @param min_feat Minimum number of differentially expressed features between siblings. If this
#' number is not reached, two neighboring clusters (siblings) in the pruned dendrogram get joined. (default: 5)
#' @param max_pval Maximum \emph{P}-value for differential expression computation. (default: 1e-4)
#' @param min_fc Mimimum fold-change for differential expression computation. (default: 2)
#' @return An unpdated object of class \code{CellTrailsSet}
#' @details To identify cellular subpopulations, CellTrails performs hierarchical clustering via minimization
#' of a square error criterion (Ward, 1963) in the lower-dimensional space. To determine the cardinality
#' of the clustering, CellTrails conducts an unsupervised \emph{post-hoc} analysis. Here, it is assumed that
#' differential expression of assayed features determines distinct cellular stages. First, Celltrails
#' identifies the maximal fragmentation of the data space, i.e. the lowest cutting height in the clustering
#' dendrogram that ensured that the resulting clusters contained at least a certain fraction of samples.
#' Then, processing from this height towards the root, CellTrails iteratively joins siblings if they did
#' not have at least a certain number of differentially expressed features. Statistical significance is tested
#' by means of a two-sample non-parametric linear rank test accounting for censored values (Peto & Peto, 1972).
#' The null hypothesis is rejected using the Benjamini-Hochberg (Benjamini & Hochberg, 1995) procedure for
#' a given significance level.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the samples stored in the \code{CellTrailsSet} object were
#' not embedded yet (ie. the \code{CellTrailsSet} object does not contain a
#' latent space matrix object; \code{latentSpace(object)}is \code{NULL}, see \code{\link[CellTrails]{embedSamples}}).
#' A warning is shown if the dimensionality of the latent space equals the number of samples. If
#' the user intends to use the non-truncated latent space this message can be ignored, otherwise it is
#' suggested to reduce the dimensionality (see \code{\link[CellTrails]{reduceDimensions}}).
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
#' ctset
#'
#' plot(ctset, type="stateSize")
#' @references Ward, J.H. (1963). Hierarchical Grouping to Optimize an Objective Function.
#' Journal of the American Statistical Association, 58, 236-244.
#' @references Peto, R., and Peto, J. (1972).
#' Asymptotically Efficient Rank Invariant Test Procedures (with Discussion).
#' Journal of the Royal Statistical Society of London, Series A 135, 185–206.
#' @references Benjamini, Y., and Hochberg, Y. (1995).
#' Controlling the false discovery rate: a practical and powerful approach to multiple
#' testing. Journal of the Royal Statistical Society Series B 57, 289–300.
#' @docType methods
#' @aliases findStates,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("findStates", function(ctset, min_size=0.01, min_feat=5, max_pval=1e-4,
                                  min_fc=2) standardGeneric("findStates"))
setMethod("findStates", "CellTrailsSet", function(ctset, min_size, min_feat, max_pval, min_fc){
  #Pre-flight check
  if(is.null(latentSpace(ctset))) {
    stop("Samples were not embedded yet, i.e., latent space information is missing. Please,
         embed the samples first (see '?embedSamples').")
  }
  if(nrow(latentSpace(ctset)) == ncol(latentSpace(ctset))) {
    warning("Size of spectrum equals number of samples. Please, make sure that
            dimensionality was reduced (see '?findSpectrum' and '?reduceDimensions').")
  }
  #Run
  if(min_size < 1) {
    min_size <- round(dim(ctset)[2] * min_size)
  }
  .findStates_def(x=ctset, link.method="ward.D2", min.size=min_size,
                 max.pval=max_pval, min.fc=min_fc, min.g=min_feat,
                 show.plots=FALSE, reverse=FALSE, verbose=FALSE)
})

# #' Differential feature expression between trajectory states
# #'
# #' Calculates differential feature expression statistics (\emph{P}-value and fold-change) between
# #' two given states.
# #' @param ctset A \code{\link[CellTrails]{CellTrailsSet}} object
# #' @param state1 Name of first state
# #' @param state2 Name of second state
# #' @param feature_name Name of feature
# #' @param alternative A character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
# #' @return A list containing the following components:
# #' \describe{
# #'   \item{\code{p.value}}{The \emph{P}-value for the test}
# #'   \item{\code{fold}}{The fold-change from state1 to state2}
# #' }
# #' @details Statistical significance is tested by means of a two-sample non-parametric linear rank test
# #' accounting for censored values (Peto & Peto, 1972).
# #' \cr \cr
# #' \emph{Diagnostic messages}
# #' \cr \cr
# #' An error is thrown if the feature name is not known/listed on the assay. Since \code{CellTrailsSet}
# #' extends class \code{ExpressionSet}, all feature names stored in a \code{CellTrailsSet} object
# #' can be retrieved by the function \code{\link[Biobase]{featureNames}}.
# #' @examples
# #' @references Peto, R., and J. Peto. (1972).
# #' Asymptotically Efficient Rank Invariant Test Procedures (with Discussion).
# #' Journal of the Royal Statistical Society of London, Series A 135, 185–206.
# #' @docType methods
# #' @aliases diffExprState,CellTrailsSet-method
# #' @export
# #' @author Daniel C. Ellwanger
# setGeneric("diffExprState", function(ctset, state1, state2, feature_name,
#                                     alternative = "two.sided") standardGeneric("diffExprState"))
# setMethod("diffExprState", "CellTrailsSet", function(ctset, state1, state2, feature_name,
#                                                      alternative){
#   #Pre-flight check
#   .checkFeatureNameExists(ctset, feature_name)
#   #Run
#   .calcDiffExpr_def(ctset, state1=state1, state2=state2,
#                     feature_name=feature_name, alternative=alternative)
# })

#' Connect trajectory states
#'
#' Connects states using maximum interface scoring. For each state an interface score is
#' defined by the relative distribution of states in its local neighborhood. A filter is
#' applied to remove outliers (ie. false positive neighbors). States are spanned by
#' maximizing the total interface score.
#' @param ctset A \code{\link[CellTrails]{CellTrailsSet}} object
#' @param l Neighborhood size (default: 10)
#' @return An updated \code{\link[CellTrails]{CellTrailsSet}} object
#' @details CellTrails assumes that the arrangement of samples in the computed
#' lower-dimensional latent space constitutes a trajectory. Therefore, CellTrails
#' aims to place single samples along a maximum parsimony tree, which resembles a
#' branching developmental continuum. Distances between samples in the latent space
#' are computed using the Euclidean distance. \cr \cr
#' To avoid overfitting and to facilitate
#' the accurate identification of bifurcations, CellTrails simplifies the problem. Analogous
#' to the idea of a ‘broken-stick regression’, CellTrails groups the data and perform linear
#' fits to separate trajectory segments, which are determined by the branching
#' chronology of states. This leaves the optimization problem of finding the minimum number
#' of associations between states while maximizing the total parsimony, which in theory can be solved
#' by any minimum spanning tree algorithm. CellTrails adapts this concept by assuming that adjacent
#' states should be located nearby and therefore share a relative high number of neighboring cells. \cr \cr
#' Each state defines a submatrix of samples that is composed of a distinct set of data vectors,
#' i.e. each state is a distinct set of samples represented in the lower-dimensional space.
#' For each state CellTrails identifies the \emph{l}-nearest neighbors to each state's data vector and
#' takes note of their state memberships and distances. This results in two vectors
#' of length \emph{l} times the state size (i.e., a vector with memberships and a vector with distances).
#' \cr \cr CellTrails removes spurious neighbors (outliers), whose distance to a state is greater than or equal to
#' \deqn{e^{median(log(D)) + MAD(log(D))}}
#' where D is a matrix containing all collected {l}-nearest neighbor sample distances to any
#' state in the latent space. \cr \cr
#' For each state CellTrails calculates the relative frequency on how often a state occurs in the neighborhood
#' of a given state, which is refered to as the interface cardinality scores. \cr \cr
#' CellTrails implements a greedy algorithm to find the tree maximizing the total interface cardinality score,
#' similar to a minimum spanning tree algorithm (Kruskal, 1956). In a nutshell, all interface cardinality
#' scores are organized in a sorted linked list, and a graph with no edges, but k nodes (one for each state)
#' is initialized. During each iteration the highest score is selected, removed from the list and its corresponding
#' edge (connecting two states), if it is not introducing a cycle or is already existent, is added to the graph.
#' The algorithm terminates if the size of the graph is \emph{k}-1 (with \emph{k} equals number of states) or the
#' list is empty. A cycle is determined if nodes were revisited while traversing the graph using depth-first search.
#' Its construction has a relaxed requirement (number of edges < number of nodes) compared to a tree
#' (number of edges = number of nodes - 1), which may result in a graph (forest) having multiple tree components,
#' i.e. several trajectories or isolated nodes.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the states have not been defined yet; function \code{\link[CellTrails]{findStates}}
#' needs to be called first.
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
#' ctset
#'
#' plot(ctset, type = "stateTrajectoryGraph",
#'     feature_name="feature_1", component=1)
#' plot(ctset, type="stateTrajectoryGraph", pheno_type="age",
#'     component=1, point_size=2)
#' @references Kruskal, J.B. (1956). On the shortest spanning subtree of a graph and the
#' traveling salesman problem. Proc Amer Math Soc 7, 48-50.
#' @docType methods
#' @aliases connectStates,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("connectStates", function(ctset, l=10) standardGeneric("connectStates"))
setMethod("connectStates", "CellTrailsSet", function(ctset, l){
  #Pre-flight checks
  if(is.null(states(ctset))) {
    stop("States have not been defined yet. Please, call function 'findStates' first.")
  }
  #Run
  .connectStates_def(ctset, l=l)
})

#' Select component from trajectory graph
#'
#' Retains a single component of a trajectory graph.
#' @param ctset A \code{\link[CellTrails]{CellTrailsSet}} object
#' @param component Number of component to be selected
#' @return An updated \code{\link[CellTrails]{CellTrailsSet}} object
#' @details The construction of a trajectory graph may result in a forest having
#' multiple tree components, which may represent individual trajectories or isolated nodes.
#' This method should be used to extract a single component from the graph. A component is
#' identified by its (integer) number.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the states have not been connected yet; function \code{\link[CellTrails]{connectStates}}
#' needs to be called first. An error is thrown if an unknown component (number) is selected.
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
#' stateTrajectoryGraph(ctset)
#' trajectorySamples(ctset)
#' @docType methods
#' @aliases selectTrajectory,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("selectTrajectory", function(ctset, component) standardGeneric("selectTrajectory"))
setMethod("selectTrajectory", "CellTrailsSet", function(ctset, component){
  #Pre-flight check
  if(is.null(stateTrajectoryGraph(ctset))) {
    stop("Trajectory tree has not been computed yet. Please, call function
         'connectStates' first.")
  }
  if(component > length(stateTrajectoryGraph(ctset))) {
    stop("Component ", component, " is not contained in this 'CellTrailsSet' object.
         Please, make sure the right component number was selected.")
  }
  #Run
  .selectTrajectory_def(ctset, component)
})

#' Align samples to trajectory
#'
#' Orthogonal projection of each sample to the trajectory backbone.
#' @param ctset A \code{\link[CellTrails]{CellTrailsSet}} object
#' @return An updated \code{\link[CellTrails]{CellTrailsSet}} object
#' @details The previously selected component (with \emph{k} states) defines
#' the trajectory backbone. With this function CellTrails embeds the trajectory
#' structure in the latent space by computing \emph{k}-1 straight lines passing
#' through \emph{k} mediancentres (Bedall & Zimmermann, 1979) of adjacent states.
#' Then, a fitting function is learned. Each sample is projected to its most proximal
#' straight line passing through the mediancentre of its assigned state. Here,
#' whenever possible, projections on line segments \emph{between} two mediancentres
#' are preferred. Residuals (fitting deviations) are given by the Euclidean distance
#' between the sample's location and the straight line. Finally, a weighted acyclic
#' trajectory graph can be constructed based on each sample’s position along its straight
#' line. In addition, data vectors are connected to mediancentres to enable the proper
#' determination of branching points. Each edge is weighted by the distance between each node
#' (sample) after orthogonal projection.
#' \cr \cr
#' Of note, the fitting function implies potential side branches in the trajectory graph; those could
#' be caused due to technical variance or encompass samples that were statistically indistinguishable
#' from the main trajectory given the selected genes used for trajectory reconstruction.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if an trajectory graph component was not computed or selected yet; functions \code{\link{connectStates}}
#' and \code{\link{selectTrajectory}} need to be run first.
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
#' plot(ctset, type="trajectoryFit")
#' @references Bedall, F.K., and Zimmermann, H. (1979).
#' Algorithm AS143. The mediancentre. Appl Statist 28, 325-328.
#' @docType methods
#' @aliases fitTrajectory,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("fitTrajectory", function(ctset) standardGeneric("fitTrajectory"))
setMethod("fitTrajectory", "CellTrailsSet", function(ctset){
  #Pre-flight check
  if(length(ctset@spanForest) != 1) {
    stop("None or multiple components in trajectory graph detected. Please, connect states
      (function 'connectStates') and select a single component (function 'selectTrajectory') first.")
  }
  #Run
  .fitTrajectory_def(ctset)
})

# #' Imputation of drop-outs
# #'
# #' Detects and interpolates missing values (drop-outs) of each gene
# #' by evaluating each sample's neighbor expression along the trajectory.
# #' @param ctset A \code{\link[CellTrails]{CellTrailsSet}} object
# #' @param feature_name Names of features which should be imputed (optional)
# #' @return An updated \code{\link[CellTrails]{CellTrailsSet}} object
# #' @details 'Drop-out' events denote measurements failing to produce a signal due to intrinsic
# #' cellular conditions, such as a low target RNA concentration or the
# #' stochastic nature of gene expression (stochastic bursts), or due to technical noise
# #' (e.g., primer dimerization). Drop-outs are false negative expression signals,
# #' i.e. can be interpreted as missing data. Those values can be recovered
# #' computationally increasing the accuracy of downstream analyses, such as
# #' differential expression analysis. CellTrails implements an imputation method
# #' based on a simple assumption: a missing value (non-detect) which is enclosed by
# #' two actual measurements having an expression level > 0 along the trajectory is assumed to
# #' denote a biological unreasonable feature/gene expression fluctuation and is therefore
# #' classified as drop-out; those missing values get linearly approximated. \cr \cr
# #' For example: Let's assume we have given the chronologically
# #' ordered expression subsequence for a given gene \eqn{X}: \eqn{(x_1, x_2, x_3)}.
# #' If \eqn{x_2 = 0} and \eqn{x_1 > 0} and \eqn{x_3 > 0} then
# #' \eqn{x_2} would become \eqn{x_2 = 0.5(x_1 + x_3)}, otherwise \eqn{x_2 = 0}. \cr \cr
# #' Please note that this computation is of reasonable runtime for big data sets.
# #' The parameter \code{feature_name} can be used to define the set of features for which
# #' the imputation should be performed, otherwise all features are used by default.
# #' \cr \cr
# #' \emph{Diagnostic messages}
# #' \cr \cr
# #' An error is thrown if the trajectory has not been computed yet; function
# #' \code{\link[CellTrails]{fitTrajectory}} needs to be called first.
# #' @docType methods
# #' @aliases imputeDropouts,CellTrailsSet-method
# #' @importFrom Biobase featureNames
# #' @export
# #' @author Daniel C. Ellwanger
# setGeneric("imputeDropouts", function(ctset, feature_name=featureNames(ctset)) standardGeneric("imputeDropouts"))
# setMethod("imputeDropouts", "CellTrailsSet", function(ctset, feature_name){
#   #Pre-flight check
#   if(is.null(ctset@trajectory$traj)) {
#     stop("Trajectory fit not found. Please, call function 'fitTrajectory' first.")
#   }
#   cnt <- sum(!feature_name %in% featureNames(ctset))
#   if(cnt > 0) {
#     stop(cnt, " feature name(s) could not be found in this object.")
#   }
#
#   #Run
#   .imputeDropouts_def(ctset, feature_name = feature_name)
# })

#' Export trajectory graph
#'
#' Writes graphml file containing the trajectory graph's structure.
#' @param object A \code{\link[CellTrails]{CellTrailsSet}} object
#' @param file A character string naming a file
#' @param feature_name A character string specifying by which feature expression nodes should be colorized
#' @param pheno_type A character string specifying by which phenotype label nodes should be colorized
#' @param label Defines the label name (optional). Can be either set to the
#' sample name with \code{NAME} or to a phenotype label name.
#' @return \code{write.ygraphml} returns an invisible \code{NULL}
#' @details To visualize the trajectory graph, a proper graph layout has to be computed.
#' Ideally, edges should not cross and nodes should not overlap. CellTrails enables the
#' export and import of the trajectory graph structure using the graphml file format.
#' This file format can be interpreted by most third-party graph analysis applications,
#' allowing the user to subject the trajectory graph to a wide range of (tree)
#' layout algorithms. In particular, its format has additional ygraph attributes best suited to be
#' used with the Graph Visualization Software 'yEd' which is freely available from yWorks GmbH
#' (http://www.yworks.com/products/yed) for all major platforms. \cr\cr
#' The colors of the nodes can be defined by the parameters \code{feature_name} and \code{pheno_type}.
#' Please note that the trajectory landmarks can be colorized via \code{pheno_type = 'landmark'}.
#' If a layout is already present in the provided \code{CellTrailsSet} object, the
#' samples' coordinates will be listed in the graphml file.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the trajectory has not been computed yet; function
#' \code{\link[CellTrails]{fitTrajectory}} needs to be called first.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' \dontrun{
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
#' # Export trajectory graph structure to graphml
#' # color nodes by gene expression (e.g, feature_10)
#' write.ygraphml(ctset, file="yourFilePath",
#'               feature_name="feature_10")
#' # color nodes by metadata (e.g., simulated age of sample) and
#' # label by computed state
#' write.ygraphml(ctset, file="yourFilePath",
#'               pheno_type="age", label="state")}
#' @seealso \code{\link[CellTrails]{read.ygraphml}}
#' @docType methods
#' @aliases write.ygraphml,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("write.ygraphml", function(object, file, feature_name=NULL, pheno_type=NULL, label=NULL) standardGeneric("write.ygraphml"))
setMethod("write.ygraphml", "CellTrailsSet", function(object, file, feature_name, pheno_type, label){
  #Pre-flight check
  if(is.null(object@trajectory$traj)) {
    stop("Trajectory fitting information not found. Please, call function 'fitTrajectory' first.")
  }
  #Run
  .write_ygraphml_def(x=object, file=file, feature_name=feature_name, pheno_type=pheno_type, label=label)
})

#' Import trajectory graph layout
#'
#' Reads ygraphml file containing the trajectory graph's layout
#' @param ctset A \code{\link[CellTrails]{CellTrailsSet}} object
#' @param file A character string naming a file
#' @param adjust Indicates if layout has to be adjusted such that edge lengths
#' correlated to pseudotime. (default: TRUE)
#' @return An updated \code{\link[CellTrails]{CellTrailsSet}} object
#' @details To visualize the trajectory graph, a proper graph layout has to be computed.
#' Ideally, edges should not cross and nodes should not overlap. CellTrails enables the
#' export and import of the trajectory graph structure using the graphml file format.
#' This file format can be interpreted by most third-party graph analysis applications,
#' allowing the user to subject the trajectory graph to a wide range of
#' layout algorithms. Please note that the graphml file needs to contain layout information
#' ("<y:Geometry x=... y=... >" entries) as provided by the 'ygraphml' file definition
#' used by the Graph Visualization Software 'yEd' (freely available from yWorks GmbH,
#' http://www.yworks.com/products/yed). \cr \cr
#' CellTrails implements a module which can incorporate pseudotime information into the the graph
#' layout (activated via parameter \code{adjust}). Here, edge lengths between two nodes (samples)
#' will then correspond to the inferred pseudotime that separates two samples along the trajectory.
#' @examples
#' # Generate example data
#' dat <- exDat()
#'
#' \dontrun{
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
#' ctset <- fitTrajectory(ctset) # Align samples to trajectory
#'
#' # Then: export trajectory graph structure
#' # for layout computation in yEd (see vignette)
#' # and reimport to the R environment
#' write.ygraphml(ctset, file="yourFilePath")
#' ctset <- read.ygraphml(ctset, file="yourFilePath")
#' }
#' @seealso \code{\link[CellTrails]{write.ygraphml}}
#' @docType methods
#' @aliases read.ygraphml,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("read.ygraphml", function(ctset, file, adjust=TRUE) standardGeneric("read.ygraphml"))
setMethod("read.ygraphml", "CellTrailsSet", function(ctset, file, adjust){
  .read_ygraphml_def(x=ctset, file=file, format="graphml", adjust=adjust)
})

#' Fit expression dynamic
#'
#' Fits feature expression as a function of pseudotime along a defined trail.
#' @param ctset A \code{\link[CellTrails]{CellTrailsSet}} object
#' @param feature_name Name of feature
#' @param trail_name Name of trail
#' @return An object of type \code{\link[base]{list}} with the following components
#' \describe{
#'   \item{\code{pseudotime}}{The pseudotime along the trail}
#'   \item{\code{expression}}{The fitted expression values for each value of pseudotime}
#'   \item{\code{gam}}{A object of class \code{\link{gamObject}}}
#' }
#' @details A trail is an induced subgraph of the trajectory graph. A trajectory graph is
#' composed of samples (nodes) that are connected (by weighted edges) if they are chronologically
#' related. A trail has to be defined by the user using \code{\link{addTrail}}.
#' A pseudotime vector is extracted by computing the geodesic distance for each sample
#' from the trail's start node. To infer the expression level of a feature as a function of
#' pseudotime, CellTrails used generalized additive models with a single smoothing term
#' with four basis dimensions. Here, for each feature CellTrails introduces prior weights
#' for each observation to lower the confounding effect of drop-outs to the
#' maximum-likelihood-based fitting process as follows. Each non-detect of feature
#' \emph{j} in state \emph{h} is weighted by the relative fraction of non-detects of
#' feature \emph{j} in state \emph{h}; detected values are always assigned weight = 1.
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
#' # Fit dynamic
#' fit <- fitDynamic(ctset, feature_name="feature_3", trail_name="Tr1")
#'
#' summary(fit)
#' @docType methods
#' @aliases fitDynamic,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("fitDynamic", function(ctset, feature_name, trail_name) standardGeneric("fitDynamic"))
setMethod("fitDynamic", "CellTrailsSet", function(ctset, feature_name, trail_name){
  #Pre-flight check
  .checkFeatureNameExists(ctset, feature_name)
  if(!trail_name %in% trailNames(ctset)) {
    stop("Name of trail is not known. Please, select a valid trail name (see function 'trailNames').")
  }
  trail <- ctset@trails[[trail_name]]
  dat <- data.frame(X = trail$ptime,
                    Y = ctset[feature_name, trail$samples],
                    STATES = states(ctset)[trail$samples])
  fit <- .fitDynamic_def(x = dat[,1] / max(dat[,1]), y = dat[,2], z = dat[,3], k=5) #fixed to 5
  names(fit) <- c("pseudotime", "expression", "gam")
  fit
})

#' Differential trail expression analysis
#'
#' Comparison of feature expression dynamic between two trails.
#' @param ctset A \code{\link[CellTrails]{CellTrailsSet}} object
#' @param feature_name Name of feature; can be multiple names
#' @param trail_names Name of trails
#' @param score Score type; one of \{"RMSD", "TAD", "ABC"\}
#' @return Numeric value
#' @details Genes have non-uniform expression rates and each trail has a distinct
#' set of upregulated genes, but also contains unequal numbers of cells. Because
#' pseudotime is based on transcriptional change, its axis may be distorted,
#' leading to stretched or compressed sections of longitudinal expression data that
#' make comparison of trails challenging. To align different trails, despite these
#' differences, CellTrails employs a dynamic programming based algorithm that has
#' long been known in speech recognition, called dynamic time warping
#' (Sakoe and Chiba, 1978). RNA expression rates are modeled analogous to speaking
#' rates (Aach and Church, 2001); the latter accounts for innate non-linear variation
#' in the length of individual phonemes (i.e. states) resulting in stretching and
#' shrinking of word (i.e. trail) segments. This allows the computation of inter-trail
#' alignment warps of individual expression time series that are similar but
#' locally out of phase. \cr \cr
#' A trail is defined as a chronologically ordered sequence of samples representing
#' feature expression snapshots at distinct points in pseudotime. Expression values
#' of given features per trail are fitted and smoothed using \code{\link{gam}}, which
#' also enables the projection of expression values for non-observed time points.
#' Univariate pairwise alignments are computed resulting in one warp per feature and
#' per trail set. First, the pseudotime axis is unified by projecting feature
#' expression for a given sequence of reference pseudotime points. Then, the
#' corresponding expression value vectors are normalized by \eqn{x / max(x)}.
#' The cross-distance matrix between both time series denotes the input to the
#' dynamic time warping algorithm. Similar to a (global) pairwise protein sequence
#' alignment, monotonicity (i.e. no time loops) and continuity (i.e. no time leaps)
#' constraints have to be imposed on the warping function to preserve temporal sequence
#' ordering. To find the optimal warp, a recursion rule is applied which selects the
#' local minimum of three moves through a dynamic programming matrix: suppose that
#' query snapshot \emph{g} and reference snapshot \emph{h} have already been aligned,
#' then the alignment of \emph{h}+1 with \emph{g}+1 is a (unit slope) diagonal move,
#' \emph{h} with \emph{g}+1 denotes an expansion by repetition of \emph{h}, and
#' \emph{h}+2 with \emph{g}+1 contracts the query by dropping \emph{h}+1. \cr \cr
#' The overall dissimilarity between two aligned expression time series \emph{x} and \emph{y}
#' of length \emph{n} is estimated by either the root-mean-square deviation
#' \eqn{RMSD(x, y) = \sqrt(\sum(x - y)^2/n)}, the total aboslute deviation \eqn{TAD(x, y) = \sum(|x-y|)},
#' the area between the aligned dynamic curves (\code{ABC}), or Pearson's
#' correlation coefficient (\code{cor}) over all aligned elements.
#' @references Sakoe, H., and Chiba, S. (1978). Dynamic programming algorithm
#' optimization for spoken word recognition. IEEE Transactions on Acoustics,
#' Speech, and Signaling Processing 26, 43-49.
#' @references Aach, J., and Church, G.M. (2001). Aligning gene expression
#' time series with time warping algorithms. Bioinformatics 17, 495-508.
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
#' # Differential expression between trails
#' contrastTrailExpr(ctset, feature_name=c("feature_1", "feature_10"),
#'                  trail_names=c("Tr1", "Tr2"), score="rmsd")
#' @docType methods
#' @aliases contrastTrailExpr,CellTrailsSet-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("contrastTrailExpr", function(ctset, feature_name=featureNames(ctset), trail_names, score="rmsd") standardGeneric("contrastTrailExpr"))
setMethod("contrastTrailExpr", "CellTrailsSet", function(ctset, feature_name, trail_names, score){
  #Pre-flight check
  sapply(feature_name, function(x) .checkFeatureNameExists(ctset, x))
  if(any(!trail_names %in% trailNames(ctset))) {
    stop("Name of trail(s) is not known. Please, select valid trail names (see function 'trailNames').")
  }
  if(length(trail_names) > 2) {
    warning("Provided more than two trail names. Only the first two will be used.")
  }
  if(length(trail_names) < 2) {
    stop("Please, provide two trail names.")
  }
  score <- toupper(score)
  if(!score %in% c("RMSD", "TD", "ABC", "COR")) {
    stop("Score method unknown. Please, select one of {'RMSD', 'TD', 'ABC', 'COR'}.")
  }
  sapply(feature_name, function(x) .diffExprTrail_def(ctset,
                                                      feature_name = x,
                                                      trail_names = trail_names,
                                                      score = score)[[score]])
})

###############################################################################
### Plots
###############################################################################
#' Visualize a CellTrailsSpectrum object
#'
#' Method illustrates the automatic determination of informative latent dimensions.
#' @param x A \code{CellTrailsSpectrum} object
#' @param ... Additional arguments; not normally accessed by user directly
#' @return A \code{\link[ggplot2]{ggplot}} object
#' @details Line segement plot showing the lagged differences between
#' ordered top \emph{frac} eigenvalues (eigengaps) and the linear fit.
#' Values above the linear fit are highlighted and corresponding eigenvalues
#' were considered to be part of the informative spectrum.
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
#' #Plot
#' ggp <- plot(spectr)
#' ggp
#' @seealso \code{\link[ggplot2]{ggplot2-package}}
#' @method plot CellTrailsSpectrum
#' @export
#' @author Daniel C. Ellwanger
plot.CellTrailsSpectrum <- function(x, ...) {
  .plot_spectrum(x)
}

#' Visualize a CellTrailsSet object
#'
#' Visualizes the data contained in a \code{\link[CellTrails]{CellTrailsSet}} object
#' @param x A \code{\link[CellTrails]{CellTrailsSet}} object
#' @param type character value; defines which type of data should be shown
#' @param feature_name character value; name of feature(s) (one or multiple of \code{\link{featureNames}})
#' @param pheno_type character value; name of phenotypical information (one of \code{\link{varLabels}})
#' @param viz an object of class \code{\link[ggplot2]{ggplot}}; can be used to avoid tSNE re-calculation
#' @param perplexity numeric value; perplexity parameter of tSNE (default: 30)
#' @param seed integer value; starting number for random number generator (default: 1101)
#' @param component integer value; component of trajectory graph that should be shown
#' @param point_size numeric value; size of shown points or pie charts (default: 3)
#' @param label_offset numeric value; offset of point label (default: 2)
#' @param trail_name character value; name of trail (one of \code{\link{trailNames}})
#' @param map_type character value; type of map that should be shown (default: "full")
#' @param ... Additional arguments; not normally accessed by user directly
#' @return A \code{\link[ggplot2]{ggplot}} object
#' @references van der Maaten, L.J.P. & Hinton, G.E., 2008. Visualizing High-Dimensional
#' Data Using t-SNE. Journal of Machine Learning Research, 9, pp.2579-2605.
#' @details
#' The available usage options:
#' \itemize{
#'   \item{\code{plot(ctset, type="stateSize")}}
#'   \item{\code{plot(ctset, type="stateExpression", feature_name)}}
#'   \item{\code{plot(ctset, type="latentSpace", feature_name, pheno_type, viz, perplexity=30, seed=1101)}}
#'   \item{\code{plot(ctset, type="stateTrajectoryGraph", feature_name, pheno_type, component, point_size=3, label_offset=2, seed=1101)}}
#'   \item{\code{plot(ctset, type="trajectoryFit")}}
#'   \item{\code{plot(ctset, type="map", feature_name, pheno_type, map_type=c("full", "se", "single"))}}
#'   \item{\code{plot(ctset, type="trailblazing")}}
#'   \item{\code{plot(ctset, type="trail", trail_name)}}
#'   \item{\code{plot(ctset, type="dynamic", feature_name, trail_name)}}
#'   }
#'
#' The available visualization options for parameter \code{type}:
#' \describe{
#'   \item{\code{stateSize}}{Barplot showing the absolute number of samples per state.}
#'   \item{\code{stateExpression}}{Violin plots showing the expression distribution of a
#'   feature per state. Each point displays the feature’s expression value in a single sample.
#'   The feature is defined by \code{feature_name}.}
#'   \item{\code{latentSpace}}{Shows the computed manifold in two dimensions using t-distributed
#'   stochastic neighbor embedding (tSNE; van der Maaten and Hinton 2008). Each point in this plot
#'   represents a sample. Points can be colorized according to feature expression or experimental
#'   metadata. The points' coloration can be defined via the attributes \code{feature_name} or
#'   \code{pheno_type}, respectively. A previously computed visualization can be
#'   reused to avoid recalculation of the tSNE via the parameter \code{viz}. The parameters
#'   \code{seed} and \code{perplexity} are used for the tSNE calculation.}
#'   \item{\code{stateTrajectoryGraph}}{Shows a single tree component of the computed trajectory graph.
#'   Each point in this plot represents a state and can be colorized according to feature
#'   expression (mean expression per state) or experimental metadata (arithmetic mean or
#'   percentage distribution of categorial values). The component is defined by parameter
#'   \code{component}. If the trajectory graph contains only a single component, then this
#'   parameter can be left undefined. The points' coloration can be defined via the
#'   attributes \code{feature_name} or \code{pheno_type}. Missing sample lables are recovered using
#'   nearest neighbor learning.}
#'   \item{\code{trajectoryFit}}{Two-dimensional visualization of the trajectory graph fit.
#'   Illustrates the samples’ geodesic distance from the endpoints of the longest path in the
#'   trajectory graph (= pseudotime). The perpendicular dispersion is proportional to the distance
#'   of a sample from the predicted trajectory (black line) in the latent space (= the residuals of
#'   the trajectory fit).}
#'   \item{\code{map}}{Two-dimensional visualization of the trajectory. The red line represents
#'   the trajectory and individual points denote samples. This plot type can either show the
#'   topography of a given feature’s expression landscape or colorizes
#'   individual samples by a metadata label. The feature is selected with parameter \code{feature_name},
#'   the metadata label with \code{pheno_type}, respectively. To show feature expression, a surface is
#'   fitted using isotropic (i.e. same parameters for both map dimensions)
#'   thin-plate spline smoothing in \code{\link{gam}}. It gives an overview of expression dynamics
#'   along all branches of the trajectory. The parameter \code{map_type} defines if either the full
#'   fitted expression surface should be shown (\code{map_type="full"}) or the standard error
#'   of the surface prediction (\code{map_type="se"}), or the expression values of single samples
#'   only (\code{map_type="single"}).}
#'   \item{\code{trailblazing}}{Visualizes landmarks, namely
#'   branch (B), end points (H), and user-defined nodes (U) on the trajectory map. This plot is helpful
#'   to identify and extract individual trails located between two landmarks.}
#'   \item{\code{trail}}{Highlights an individual trail (as selected via the parameter
#'   \code{trail_name}) on the trajectory map.}
#'   \item{\code{dynamic}}{Shows the expression of a feature as a function of pseudotime along a
#'   given trail. Feature name(s) are defined by parameter \code{feature_name}, the trail is defined by
#'   parameter \code{trail_name}. Points represent single samples colorized by state and the line is the
#'   fitted dynamic. If multiple features are selected, then only their fitted dynamics are shown.
#'   The dynamic is fitted using a \code{\link{gam}} with a single
#'   smoothing term with four basis dimension and prior weights as defined in function
#'   \code{\link{fitDynamic}}.}
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
#' # Plot state sizes
#' plot(ctset, type="stateSize")
#'
#' # Plot gene expression per state
#' plot(ctset, type="stateExpression", feature_name="feature_1")
#'
#' # Plot manifold in 2D
#' ggp <- plot(ctset, type="latentSpace", feature_name="feature_1") #gene expression
#' ggp
#' plot(ctset, type="latentSpace", pheno_type="age", viz=ggp) #metadata
#'
#' # Plot maximum interface tree
#' plot(ctset, type="stateTrajectoryGraph", feature_name="feature_1") #gene expression
#' plot(ctset, type="stateTrajectoryGraph", pheno_type="age", point_size=2) #metadata
#'
#' # Plot trajectory fit residuals
#' plot(ctset, type="trajectoryFit")
#'
#' # Plot CellTrails maps
#' plot(ctset, type="map", feature_name="feature_10") #gene expression
#' plot(ctset, type="map", feature_name="feature_10", map_type="se") #standard error
#' plot(ctset, type="map", feature_name="feature_10", map_type="single") #gene expression
#' plot(ctset, type="map", pheno_type="age") #metadata
#'
#' # Plot landmarks on map
#' plot(ctset, type="trailblazing")
#'
#' # Highlight individual trails on map
#' plot(ctset, type="trail", trail_name="Tr1")
#'
#' # Plot expression dynamics
#' plot(ctset, type="dynamic", feature_name="feature_3", trail_name="Tr1")
#' plot(ctset, type="dynamic", feature_name=c("feature_1", "feature_10"), trail_name="Tr2")
#' @seealso \code{\link[ggplot2]{ggplot2-package}}
#' @method plot CellTrailsSet
#' @export
#' @author Daniel C. Ellwanger
plot.CellTrailsSet <- function(x, type, feature_name=NULL, pheno_type=NULL, viz=NULL,
                               perplexity=30, seed=1101,
                               component=NULL, point_size=3, label_offset=2,
                               trail_name=NULL, map_type=c("full", "se", "single"),
                               ...) {
  switch(toupper(type),
         STATESIZE = .plot_stateSize(x),
         STATEEXPRESSION = .plot_stateExpression(x, feature_name=feature_name),
         LATENTSPACE = .plot_latentSpace(x, pheno_type=pheno_type,
                                         feature_name=feature_name,
                                         seed=seed,
                                         perplexity=perplexity,
                                         viz=viz, ...),
         STATETRAJECTORYGRAPH = .plot_stateTrajectoryGraph(x,
                                                           feature_name=feature_name,
                                                           pheno_type=pheno_type,
                                                           component=component,
                                                           point_size=point_size,
                                                           label_offset=label_offset,
                                                           seed=seed),
         TRAJECTORYFIT = .plot_trajectoryFit(x, ...),
         MAP = .plot_map(x, feature_name=feature_name,
                         pheno_type=pheno_type,
                         map_type=map_type),
         TRAILBLAZING = .plot_trailblazing(x),
         TRAIL = .plot_trail_on_map(x, trail_name=trail_name),
         DYNAMIC = .plot_dynamic(x, feature_name=feature_name, trail_name=trail_name),
         stop("Plot type '", type, "' is unknown for class 'CellTrailsSet'.")
  )
}
