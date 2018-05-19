#' @include AllClasses.R
NULL

#' Filter trajectory features by Detection Level (DL)
#'
#' Filters trajectory features that are detected in a minimum number of
#' samples.
#' @param sce An \code{SingleCellExperiment} object
#' @param threshold Minimum number of samples; if value < 1 it is interpreted
#' as fraction, otherwise as absolute sample count
#' @param show_plot Indicates if plot should be shown (default: TRUE)
#' @return A \code{character} vector
#' @details The detection level denotes the fraction of samples in which a
#' feature was detected. For each trajectory feature listed in the
#' CellTrailsSet object the relative number of samples having a feature
#' expression value greater than 0 is counted. Features that are expressed in
#' a fraction of all samples greater than \code{threshold} remain labeled as
#' trajectory feature as listed in the \code{SingleCellExperiment} object,
#' otherwise they may be not considered for dimensionality reduction,
#' clustering and trajectory reconstruction. If the parameter \code{threshold}
#' fullfills \code{threshold} \eqn{>= 1} it becomes converted to a relative
#' fraction of the total sample count. Please note that spike-in controls
#' are ignored and are not listed as trajectory features.
#' @seealso \code{trajFeatureNames} \code{isSpike}
#' @examples
#' # Example data
#' set.seed(1101)
#' dat <- simulate_exprs(n_features=15000, n_samples=100)
#'
#' # Create container
#' alist <- list(logcounts=dat)
#' sce <- SingleCellExperiment(assays=alist)
#'
#' # Filter features
#' tfeat <- filterTrajFeaturesByDL(sce, threshold=2)
#' head(tfeat)
#'
#' # Set trajectory features to object
#' trajFeatureNames(sce) <- tfeat
#'
#' # Number of features
#' length(trajFeatureNames(sce)) #filtered
#' nrow(sce) #total
#' @docType methods
#' @aliases filterTrajFeaturesByDL,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("filterTrajFeaturesByDL", function(sce, threshold, show_plot=TRUE)
  standardGeneric("filterTrajFeaturesByDL"))
setMethod("filterTrajFeaturesByDL", "SingleCellExperiment",
          function(sce, threshold, show_plot){
   y <- .exprs(sce[.useFeature(sce), ])
  .filterTrajFeaturesByDL_def(y=y, threshold=threshold, show_plot=show_plot)
})

#' Filter features by Coefficient of Variation (COV)
#'
#' Filters trajectory features by their coefficient of variation.
#' @param sce An \code{SingleCellExperiment} object
#' @param threshold Minimum coefficient of variation;
#' numeric value between 0 and 1
#' @param design A numeric matrix describing the factors that should be blocked
#' @param show_plot Indicates if plot should be shown (default: TRUE)
#' @return A \code{character} vector
#' @details For each trajectory feature \emph{x} listed in the
#' \code{SingleCellExperiment} object the coefficient of variation is
#' computed by \eqn{CoV(x) = sd(x) / mean(x)}. Features with a CoV(x) greater
#' than \code{threshold} remain labeled as trajectory feature in the
#' \code{SingleCellExperiment} object, otherwise they are not considered
#' for dimensionality reduction, clustering and trajectory reconstruction.
#' Please note that spike-in controls are ignored
#' and are not listed as trajectory features.
#' \cr \cr
#' To account for systematic bias in the expression data
#' (e.g., cell cycle effects), a design matrix can be provided for the
#' learning process. It should list the factors that should be blocked and
#' their values per sample. It is suggested to construct a design
#' matrix with \code{model.matrix}.
#' @seealso \code{trajFeatureNames} \code{isSpike} \code{model.matrix}
#' @examples
#' # Simulate example data
#' set.seed(1101)
#' dat <- simulate_exprs(n_features=15000, n_samples=100)
#'
#' # Create container
#' alist <- list(logcounts=dat)
#' sce <- SingleCellExperiment(assays=alist)
#'
#' # Filter incrementally
#' trajFeatureNames(sce) <- filterTrajFeaturesByDL(sce, threshold=2)
#' trajFeatureNames(sce) <- filterTrajFeaturesByCOV(sce, threshold=0.5)
#'
#' # Number of features
#' length(trajFeatureNames(sce)) #filtered
#' nrow(sce) #total
#' @docType methods
#' @aliases filterTrajFeaturesByCOV,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("filterTrajFeaturesByCOV", function(sce, threshold,
                                               design=NULL, show_plot=TRUE)
  standardGeneric("filterTrajFeaturesByCOV"))
setMethod("filterTrajFeaturesByCOV", "SingleCellExperiment",
          function(sce, threshold, design, show_plot){
  y <- .exprs(sce[.useFeature(sce), ])
  .filterTrajFeaturesByCOV_def(y=y, threshold=threshold, design,
                               show_plot=show_plot)
})

#' Filter features by Fano Factor
#'
#' Filters trajectory features that exhibit a significantly high fano factor
#' (index of dispersion) by considering average expression levels.
#' @param sce An \code{SingleCellExperiment} object
#' @param threshold A Z-score cutoff (default: 1.7)
#' @param min_expr Minimum average expression of feature to be considered
#' @param design A numeric matrix describing the factors that should be blocked
#' @param show_plot Indicates if plot should be shown (default: TRUE)
#' @return A \code{character} vector
#' @details To identify the most variable features an unsupervised strategy
#' that controls for the relationship between a features’s average expression
#' intensity and its expression variability is applied. Features are placed
#' into 20 bins based on their mean expression. For each bin the fano factor
#' (a windowed version of the index of dispersion, IOD = variance / mean)
#' distribution is computed and standardized
#' (\emph{Z}-score(\emph{x}) = x/sd(\emph{x}) - mean(\emph{x})/sd(\emph{x})).
#' Features with a \emph{Z}-score
#' greater than \code{threshold} remain labeled as trajectory feature
#' in the \code{SingleCellExperiment} object. The parameter \code{min_expr}
#' defines the minimum average expression level of a feature to be
#' considered for this filter method. Please note that spike-in controls are
#' ignored and are not listed as trajectory features.
#' \cr \cr
#' To account for systematic bias in the expression data
#' (e.g., cell cycle effects), a design matrix can be provided for the
#' learning process. It should list the factors that should be blocked and
#' their values per sample. It is suggested to construct a design matrix
#' with \code{model.matrix}.
#' @seealso \code{trajFeatureNames} \code{isSpike} \code{model.matrix}
#' @examples
#' # Simulate example data
#' set.seed(1101)
#' dat <- simulate_exprs(n_features=15000, n_samples=100)
#'
#' # Create container
#' alist <- list(logcounts=dat)
#' sce <- SingleCellExperiment(assays=alist)
#'
#' # Filter incrementally
#' trajFeatureNames(sce) <- filterTrajFeaturesByDL(sce, threshold=2)
#' trajFeatureNames(sce) <- filterTrajFeaturesByCOV(sce, threshold=0.5)
#' trajFeatureNames(sce) <- filterTrajFeaturesByFF(sce, threshold=1.7)
#'
#' # Number of features
#' length(trajFeatureNames(sce)) #filtered
#' nrow(sce) #total
#' @docType methods
#' @aliases filterTrajFeaturesByFF,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("filterTrajFeaturesByFF",
           function(sce, threshold=1.7, min_expr=0,
                    design=NULL, show_plot=TRUE)
  standardGeneric("filterTrajFeaturesByFF"))
setMethod("filterTrajFeaturesByFF", "SingleCellExperiment",
          function(sce, threshold, min_expr, design, show_plot){
  y <- .exprs(sce[.useFeature(sce), ])
  .filterTrajFeaturesByFF_def(y=y, z=threshold,
                              min_expr=min_expr, design, show_plot=show_plot)
})

#' Spectral embedding of biological samples
#'
#' Non-linear learning of a data representation that captures the
#' intrinsic geometry of the trajectory. This function performs spectral
#' decomposition of a graph encoding conditional entropy-based
#' sample-to-sample similarities.
#' @param x A \code{SingleCellExperiment} object or a numeric matrix with
#' samples in columns and features in rows
#' @param design A numeric matrix describing the factors that should be blocked
#' @return A list containing the following components:
#'   \item{\code{eigenvectors}}{Ordered components of latent space}
#'   \item{\code{eigenvalues}}{Information content of latent components}
#' @details Single-cell gene expression measurements comprise high-dimensional
#' data of large volume, i.e. many features (e.g., genes) are measured in many
#' samples (e.g., cells); or more formally, \emph{m} samples can be described
#' by the expression of \emph{n} features (i.e., \emph{n} dimensions). The
#' cells’ expression profiles are shaped by many distinct unobserved biological
#' causes related to each cell's geno- and phenotype, such as developmental
#' age, tissue region of origin, cell cycle stage, as well as extrinsic sources
#' such as status of signaling receptors, and environmental stressors, but also
#' technical noise. In other words, a single dimension, despite just containing
#' gene expression information, represents an underlying combination of multiple
#' dependent and independent, relevant and non-relevant factors, whereat each
#' factors’ individual contribution is non-uniform. To obtain a better
#' resolution and to extract underlying information, CellTrails aims to find a
#' meaningful low-dimensional structure - a manifold - that represents cells
#' mainly by their temporal relation along a biological process.
#' \cr \cr
#' This method assumes that the expression vectors are lying on or near a
#' manifold with dimensionality \emph{d} that is embedded in the
#' \emph{n}-dimensional space. By using spectral embedding CellTrails aims to
#' amplify latent temporal information; it reduces noise (ie. truncates
#' non-relevant dimensions) by transforming the expression matrix into a new
#' dataset while retaining the geometry of the original dataset as much as
#' possible.CellTrails captures overall cell-to-cell relations based on the
#' statistical mutual dependency between any two data vectors. A high
#' dependency between two samples should be represented by their close
#' proximity in the lower-dimensional space.
#' \cr \cr
#' First, the mutual depencency between samples is scored using mutual
#' information. This entropy framework naturally requires discretization
#' of data vectors by an indicator function, which assigns each continuous
#' data point (expression value) to exactly one discrete interval (e.g. low,
#' mid or high). However, measurement points located close to the interval
#' borders may get wrongly assigned due to noise-induced fluctuations.
#' Therefore, CellTrails fuzzifies the indicator function by using a piecewise
#' polynomial function, i.e. the domain of each sample expression vector is
#' divided into contiguous intervals (based on Daub \emph{et al.}, 2004).
#' Second, the computed mutual information matrix, which is left-bounded and
#' composed of bits, is scaled to a generalized correlation coefficient. Third,
#' CellTrails constructs a simple complete graph with \emph{m} nodes, one for
#' each data vector (ie. sample), and weights each edge between two nodes by a
#' heat kernel function applied on the generalzied correlation coefficient.
#' Finally, nonlinear spectral embedding (ie. spectral decomposition of the
#' graph's adjacency matrix) is performed
#' (Belkin & Niyogi, 2003; Sussman \emph{et al.}, 2012) unfolding the manifold.
#' Please note that this methods only uses the set of defined trajectory
#' features in a \code{SingleCellExperiment} object; spike-in controls are
#' ignored and are not listed as trajectory features.
#' \cr \cr
#' To account for systematic bias in the expression data
#' (e.g., cell cycle effects), a design matrix can be
#' provided for the learning process. It should list the factors that should be
#' blocked and their values per sample. It is suggested to construct a
#' design matrix with \code{model.matrix}.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' The method throws an error if expression matrix contains samples
#' with zero entropy (e.g., the samples exclusively contain non-detects, that
#' is all expression values are zero).
#' @seealso \code{SingleCellExperiment} \code{trajectoryFeatureNames}
#' \code{model.matrix}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Embed samples
#' res <- embedSamples(exSCE)
#' @references Daub, C.O., Steuer, R., Selbig, J., and Kloska, S. (2004).
#' Estimating mutual information using B-spline functions -- an improved
#' similarity measure for analysing gene expression data.
#' BMC Bioinformatics 5, 118.
#' @references Belkin, M., and Niyogi, P. (2003). Laplacian eigenmaps for
#' dimensionality reduction and data representation. Neural computation 15,
#' 1373-1396.
#' @references Sussman, D.L., Tang, M., Fishkind, D.E., and Priebe, C.E.
#' (2012). A Consistent Adjacency Spectral Embedding for Stochastic Blockmodel
#' Graphs. J Am Stat Assoc 107, 1119-1128.
#' @docType methods
#' @rdname embedSamples
#' @aliases embedSamples,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("embedSamples", function(x, design=NULL)
  standardGeneric("embedSamples"))
setMethod("embedSamples", "SingleCellExperiment", function(x, design){
  M <- .exprs(x[.useFeature(x), ]) #select trajectory features
  .embedSamples_def(x=M, design=design)
})

#' @rdname embedSamples
#' @aliases embedSamples,matrix-method
setMethod("embedSamples", "matrix", function(x, design){
  .embedSamples_def(x=x, design=design)
})

#' Determine number of informative latent dimensions
#'
#' Identifies the dimensionality of the latent space
#' @param x A numeric vector with eigenvalues
#' @param frac Fraction or number (if \code{frac > 1}) of eigengaps
#' used to perform linear fit. (default: 100)
#' @return A \code{numeric} vector with indices of relevant dimensions
#' @details Similar to a scree plot, this method generates a simple line
#' segement plot showing the lagged differences between ordered eigenvalues
#' (eigengaps). A linear fit is calucated on a fraction of top ranked values
#' to identify informative eigenvectors.
#' @seealso \code{pca} \code{embedSamples}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Embedding
#' res <- embedSamples(exSCE)
#'
#' # Find spectrum
#' d <- findSpectrum(res$eigenvalues, frac=30)
#' d
#' @docType methods
#' @aliases findSpectrum,numeric-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("findSpectrum", function(x, frac=100)
  standardGeneric("findSpectrum"))
setMethod("findSpectrum", "numeric", function(x, frac){
  .findSpectrum_def(D=x, frac=frac)
})

#' Principal Component Analysis
#'
#' Performs principal component analysis by spectral decomposition of
#' a covariance or correlation matrix
#' @param sce \code{SingleCellExperiment} object
#' @param do_scaling FALSE = covariance matrix, TRUE = correlation matrix
#' @param design A numeric matrix describing the factors that should be blocked
#' @return A \code{list} object containing the following components:
#'   \item{\code{components}}{Principal components}
#'   \item{\code{eigenvalues}}{Variance per component}
#'   \item{\code{variance}}{Fraction of variance explained by each component}
#'   \item{\code{loadings}}{Loading score for each feature}
#' @details The calculation is done by a spectral decomposition of the
#' (scaled) covariance matrix of the trajectory features
#' as defined in the \code{SingleCellExperiment} object.
#' Features with zero variance get automatically removed.
#' Please note that this methods only uses the set of defined trajectory
#' features in a \code{SingleCellExperiment} object; spike-in controls are
#' ignored and are not listed as trajectory features.
#' \cr \cr
#' To account for systematic bias in the expression data
#' (e.g., cell cycle effects), a
#' design matrix can be provided for the learning process. It should list
#' the factors that should be blocked and
#' their values per sample. It is suggested to construct a design matrix with
#' \code{model.matrix}.
#' @seealso \code{SingleCellExperiment} \code{model.matrix}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Principal component analysis
#' res <- pca(exSCE)
#'
#' # Find relevant number of principal components
#' d <- findSpectrum(res$eigenvalues, frac=20)
#'
#' barplot(res$variance[d] * 100, ylab="Variance (%)",
#'         names.arg=colnames(res$components)[d], las=2)
#' plot(res$component, xlab="PC1", ylab="PC2")
#' @docType methods
#' @aliases pca,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("pca", function(sce, do_scaling=TRUE, design=NULL)
  standardGeneric("pca"))
setMethod("pca", "SingleCellExperiment", function(sce, do_scaling, design){
  M <- .exprs(sce[.useFeature(sce), ])
  if(nrow(M) == nrow(sce) & nrow(M) > 1000) {
    warning("Please note that trajectory features weren't filtered. Thus, ",
            "pca will be performed on all features, which ",
            "may result in lower accuracy and longer computation time.")
  }
  .pca_def(M=M, do_scaling=do_scaling, design=design)
})

#' Identify trajectory states
#'
#' Determines states using hierarchical spectral clustering with a
#' \emph{post-hoc} test.
#' @param sce A \code{SingleCellExperiment} object
#' @param min_size The initial cluster dedrogram is cut at an height such that
#' the minimum cluster size is at least \code{min_size};
#' if \code{min_size} < 1 than the fraction of total samples is used,
#' otherwise it is used as absoulte count (default: 0.01).
#' @param min_feat Minimum number of differentially expressed features between
#' siblings. If this number is not reached, two neighboring clusters (siblings)
#' in the pruned dendrogram get joined. (default: 5)
#' @param max_pval Maximum \emph{P}-value for differential expression
#' computation. (default: 1e-4)
#' @param min_fc Mimimum fold-change for differential expression
#' computation (default: 2)
#' @return A \code{factor} vector
#' @details To identify cellular subpopulations, CellTrails performs
#' hierarchical clustering via minimization of a square error criterion
#' (Ward, 1963) in the lower-dimensional space. To determine the cardinality
#' of the clustering, CellTrails conducts an unsupervised \emph{post-hoc}
#' analysis. Here, it is assumed that differential expression of assayed
#' features determines distinct cellular stages. First, Celltrails identifies
#' the maximal fragmentation of the data space, i.e. the lowest cutting height
#' in the clustering dendrogram that ensured that the resulting clusters
#' contained at least a certain fraction of samples. Then, processing from
#' this height towards the root, CellTrails iteratively joins siblings if
#' they did not have at least a certain number of differentially expressed
#' features. Statistical significance is tested by means of a two-sample
#' non-parametric linear rank test accounting for censored values
#' (Peto & Peto, 1972). The null hypothesis is rejected using the
#' Benjamini-Hochberg (Benjamini & Hochberg, 1995) procedure for
#' a given significance level. \cr
#' Since this methods performs pairwise comparisons, the fold change threshold
#' value is valid in both directions: higher and lower
#' expressed than \code{min_fc}. Thus, input values < 0 are interpreted as a
#' fold-change of 0. For example, \code{min_fc=2} checks for features
#' that are 2-fold differentially expressed in two given states (e.g., S1, S2).
#' Thus, a feature can be either 2-fold higher expressed in state S1 or two-fold
#' lower expressed in state S2 to be validated as differentially expressed. \cr
#' Please note that this methods only uses the set of defined trajectory
#' features in a \code{SingleCellExperiment} object; spike-in controls are
#' ignored and are not listed as trajectory features.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the samples stored in the \code{SingleCellExperiment}
#' object were not embedded yet (ie. the \code{SingleCellExperiment} object
#' does not contain a latent space matrix object; \code{latentSpace(object)}is
#' \code{NULL}).
#' @seealso \code{latentSpace} \code{trajectoryFeatureNames}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Find states
#' cl <- findStates(exSCE, min_feat=2)
#' head(cl)
#' @references Ward, J.H. (1963). Hierarchical Grouping to Optimize
#' an Objective Function. Journal of the American Statistical
#' Association, 58, 236-244.
#' @references Peto, R., and Peto, J. (1972).
#' Asymptotically Efficient Rank Invariant Test Procedures (with Discussion).
#' Journal of the Royal Statistical Society of London, Series A 135, 185–206.
#' @references Benjamini, Y., and Hochberg, Y. (1995).
#' Controlling the false discovery rate: a practical and powerful
#' approach to multiple testing. Journal of the Royal Statistical
#' Society Series B 57, 289–300.
#' @docType methods
#' @aliases findStates,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("findStates", function(sce, min_size=0.01, min_feat=5,
                                  max_pval=1e-4, min_fc=2)
  standardGeneric("findStates"))
setMethod("findStates", "SingleCellExperiment", function(sce, min_size,
                                                         min_feat, max_pval,
                                                         min_fc){
  #Pre-flight check
  if(is.null(latentSpace(sce))) {
    stop("Samples were not embedded yet, i.e., latent space information ",
         "is missing. Please, define the latent space first.")
  }

  #Run
  if(min_size < 1) {
    min_size <- round(ncol(sce) * min_size)
  }
  X <- t(.exprs(sce[.useFeature(sce), ]))
  ordi <- CellTrails::latentSpace(sce)
  .findStates_def(X=X, ordi=ordi, link.method="ward.D2", min.size=min_size,
                 max.pval=max_pval, min.fc=min_fc, min.g=min_feat,
                 reverse=FALSE, verbose=FALSE)
})

# #' Differential feature expression between trajectory states
# #'
# #' Calculates differential feature expression statistics
# #' (\emph{P}-value and fold-change) between
# #' two given states.
# #' @param ctset A \code{CellTrailsSet} object
# #' @param state1 Name of first state
# #' @param state2 Name of second state
# #' @param feature_name Name of feature
# #' @param alternative A character string specifying the alternative
# #' hypothesis, must be one of "two.sided" (default), "greater" or "less".
# #' @return A list containing the following components:
# #' \describe{
# #'   \item{\code{p.value}}{The \emph{P}-value for the test}
# #'   \item{\code{fold}}{The fold-change from state1 to state2}
# #' }
# #' @details Statistical significance is tested by means of a
# #' two-sample non-parametric linear rank test
# #' accounting for censored values (Peto & Peto, 1972).
# #' \cr \cr
# #' \emph{Diagnostic messages}
# #' \cr \cr
# #' An error is thrown if the feature name is not known/listed on the
# #' assay. Since \code{CellTrailsSet}
# #' extends class \code{ExpressionSet}, all feature names stored in
# #' a \code{CellTrailsSet} object
# #' can be retrieved by the function \code{featureNames}.
# #' @examples
# #' @references Peto, R., and J. Peto. (1972).
# #' Asymptotically Efficient Rank Invariant Test Procedures (with Discussion).
# #' Journal of the Royal Statistical Society of London, Series A 135, 185–206.
# #' @docType methods
# #' @aliases diffExprState,CellTrailsSet-method
# #' @export
# #' @author Daniel C. Ellwanger
# #'setGeneric("diffExprState", function(ctset, state1, state2, feature_name,
# #'                                     alternative = "two.sided")
# #'  standardGeneric("diffExprState"))
# #'setMethod("diffExprState", "CellTrailsSet",
# #'         function(ctset, state1, state2, feature_name, alternative){
#   #Pre-flight check
#   .featureNameExists(ctset, feature_name)
#   #Run
#   .calcDiffExpr_def(ctset, state1=state1, state2=state2,
#                     feature_name=feature_name, alternative=alternative)
# })

#' Connect trajectory states
#'
#' Connects states using maximum interface scoring. For each state an
#' interface score is defined by the relative distribution of states in
#' its local l-neighborhood. A filter is applied to remove outliers
#' (ie. false positive neighbors). States are spanned by
#' maximizing the total interface score.
#' @param sce A \code{SingleCellExperiment} object
#' @param l Neighborhood size (default: 10)
#' @return An updated \code{SingleCellExperiment} object
#' @details CellTrails assumes that the arrangement of samples
#' in the computed lower-dimensional latent space constitutes a trajectory.
#' Therefore, CellTrails aims to place single samples along a maximum parsimony
#' tree, which resembles a branching developmental continuum. Distances between
#' samples in the latent space are computed using the Euclidean distance.
#' \cr \cr
#' To avoid overfitting and to facilitate the accurate identification of
#' bifurcations, CellTrails simplifies the problem. Analogous to the idea of
#' a ‘broken-stick regression’, CellTrails groups the data and perform linear
#' fits to separate trajectory segments, which are determined by the branching
#' chronology of states. This leaves the optimization problem of finding the
#' minimum number of associations between states while maximizing the total
#' parsimony, which in theory can be solved by any minimum spanning tree
#' algorithm. CellTrails adapts this concept by assuming that adjacent states
#' should be located nearby and therefore share a relative high number of
#' neighboring cells.
#' \cr \cr
#' Each state defines a submatrix of samples that is composed of a distinct
#' set of data vectors, i.e., each state is a distinct set of samples
#' represented in the lower-dimensional space. For each state CellTrails
#' identifies the \emph{l}-nearest neighbors to each state's data
#' vector and takes note of their state memberships and distances.
#' This results in two vectors of length \emph{l} times the state size
#' (i.e., a vector with memberships and a vector with distances).
#' \cr \cr
#' CellTrails removes spurious neighbors (outliers),
#' whose distance to a state is greater than or equal to
#' \deqn{e^{median(log(D)) + MAD(log(D))}}
#' where D is a matrix containing all collected {l}-nearest neighbor sample
#' distances to any state in the latent space.
#' \cr \cr
#' For each state CellTrails calculates the relative frequency on
#' how often a state occurs in the neighborhood
#' of a given state, which is refered to as the interface cardinality scores.
#' \cr \cr
#' CellTrails implements a greedy algorithm to find the tree maximizing
#' the total interface cardinality score,
#' similar to a minimum spanning tree algorithm (Kruskal, 1956).
#' In a nutshell, all interface cardinality
#' scores are organized in a sorted linked list, and a graph
#' with no edges, but k nodes (one for each state)
#' is initialized. During each iteration the highest score is
#' selected, removed from the list and its corresponding
#' edge (connecting two states), if it is not introducing a cycle or is
#' already existent, is added to the graph.
#' The algorithm terminates if the size of the graph is
#' \emph{k}-1 (with \emph{k} equals number of states) or the
#' list is empty. A cycle is determined if nodes were revisited
#' while traversing the graph using depth-first search.
#' Its construction has a relaxed requirement (number of edges <
#' number of nodes) compared to a tree
#' (number of edges = number of nodes - 1), which may result in a
#' graph (forest) having multiple tree components,
#' i.e. several trajectories or isolated nodes.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the states have not been defined yet;
#' function \code{findStates}
#' needs to be called first.
#' @seealso \code{findStates} \code{states}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Connect states
#' exSCE <- connectStates(exSCE, l=30)
#' @references Kruskal, J.B. (1956). On the shortest
#' spanning subtree of a graph and the traveling salesman problem.
#' Proc Amer Math Soc 7, 48-50.
#' @docType methods
#' @aliases connectStates,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("connectStates", function(sce, l=10)
  standardGeneric("connectStates"))
setMethod("connectStates", "SingleCellExperiment", function(sce, l){
  #Pre-flight checks
  if(is.null(states(sce))) {
    stop("States have not been defined yet. Please, ",
         "cluster samples first.")
  }
  #Run
  dmat <- stats::dist(CellTrails::latentSpace(sce))
  cl <- CellTrails::states(sce)
  .spanForest(sce) <- .connectStates_def(dmat=dmat, cl=cl, l=l)
  sce
})

#' Select component from trajectory graph
#'
#' Retains a single component of a trajectory graph.
#' @param sce A \code{SingleCellExperiment} object
#' @param component Number of component to be selected
#' @return An updated \code{SingleCellExperiment} object
#' @details The construction of a trajectory graph may result in a forest
#' having multiple tree components, which may represent individual
#' trajectories or isolated nodes. This method should be used to extract a
#' single component from the graph. A component is
#' identified by its (integer) number.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the states have not been connected yet;
#' function \code{connectStates} needs to be called first. An
#' error is thrown if an unknown component (number) is selected.
#' @seealso \code{connectStates}
#' @seealso \code{findStates} \code{states}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Select trajectory
#' exSCE <- selectTrajectory(exSCE, component=1)
#' @docType methods
#' @aliases selectTrajectory,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("selectTrajectory", function(sce, component)
  standardGeneric("selectTrajectory"))
setMethod("selectTrajectory", "SingleCellExperiment", function(sce, component){
  #Pre-flight check
  if(is.null(.spanForest(sce))) {
    stop("Trajectory tree has not been computed yet. Please, call function
         'connectStates' first.")
  }
  if(component > length(.spanForest(sce))) {
    stop("Component ", component, " is not contained in this object.",
         "Please, make sure the right component number was selected.")
  }
  #Run
  .spanForest(sce) <- .spanForest(sce)[component]
  vnames <- names(V(.spanForest(sce)[[1]]))
  .useSample(sce) <- states(sce) %in% vnames
  sce
})

#' Align samples to trajectory
#'
#' Orthogonal projection of each sample to the trajectory backbone.
#' @param sce A \code{SingleCellExperiment} object
#' @return An updated \code{SingleCellExperiment} object
#' @details The previously selected component (with \emph{k} states) defines
#' the trajectory backbone. With this function CellTrails embeds the trajectory
#' structure in the latent space by computing \emph{k}-1 straight lines passing
#' through \emph{k} mediancentres (Bedall & Zimmermann, 1979) of adjacent
#' states. Then, a fitting function is learned. Each sample is projected to
#' its most proximal straight line passing through the mediancentre of its
#' assigned state. Here, whenever possible, projections on line segments
#' \emph{between} two mediancentres are preferred. Residuals
#' (fitting deviations) are given by the Euclidean distance between the
#' sample's location and the straight line. Finally, a weighted acyclic
#' trajectory graph can be constructed based on each sample’s position along
#' its straight line. In addition, data vectors are connected to mediancentres
#' to enable the proper determination of branching points. Each edge is
#' weighted by the distance between each node
#' (sample) after orthogonal projection.
#' \cr \cr
#' Of note, the fitting function implies potential side branches in the
#' trajectory graph; those could be caused due to technical variance or
#' encompass samples that were statistically indistinguishable from the main
#' trajectory given the selected genes used for trajectory reconstruction.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if an trajectory graph component was not
#' computed or selected yet; functions \code{connectStates}
#' and \code{selectTrajectory} need to be run first.
#' @seealso \code{connectStates} \code{selectTrajectory}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Align samples to trajectory
#' exSCE <- fitTrajectory(exSCE)
#' @references Bedall, F.K., and Zimmermann, H. (1979).
#' Algorithm AS143. The mediancentre. Appl Statist 28, 325-328.
#' @docType methods
#' @aliases fitTrajectory,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("fitTrajectory", function(sce)
  standardGeneric("fitTrajectory"))
setMethod("fitTrajectory", "SingleCellExperiment", function(sce){
  #Pre-flight check
  if(length(.spanForest(sce)) != 1) {
    stop("None or multiple components in trajectory graph detected. ",
         "Please, connect states",
         "(function 'connectStates') and select a single component ",
         "(function 'selectTrajectory') first.")
  }
  #Run
  cl <- droplevels(states(sce)[.useSample(sce)])
  g <- .spanForest(sce)[[1]]
  X <- CellTrails::latentSpace(sce)[.useSample(sce), ]
  fit <- .fitTrajectory_def(cl=cl, g=g,
                            X=X, snames=colnames(sce)[.useSample(sce)])
  .trajLandmark(sce, type="type") <- fit$lndmrk$type
  .trajLandmark(sce, type="id") <- fit$lndmrk$id
  .trajLandmark(sce, type="shape") <- fit$lndmrk$shape
  .trajResiduals(sce) <- fit$error
  .trajGraph(sce) <- fit$traj
  sce
})

#' Export trajectory graph
#'
#' Writes graphml file containing the trajectory graph's structure.
#' @param sce A \code{SingleCellExperiment} object
#' @param file Character string naming a file
#' @param color_by Indicates if nodes are colorized by a feature expression
#' ('featureName') or phenotype label ('phenoName')
#' @param name A character string specifying the featureName or phenoName
#' @param node_label Defines the node label name (optional). Can be either set
#' to the samples' states ('state') or the samples' names ('name').
#' @return \code{write.ygraphml} returns an invisible \code{NULL}
#' @details To visualize the trajectory graph, a proper graph layout has
#' to be computed. Ideally, edges should not cross and nodes should not
#' overlap (i.e., a planar embedding of the graph). CellTrails enables the
#' export and import of the trajectory
#' graph structure using the graphml file format. This file format can be
#' interpreted by most third-party graph analysis applications,
#' allowing the user to subject the trajectory graph to a wide range of (tree)
#' layout algorithms. In particular, its format has additional ygraph
#' attributes best suited to be used with the Graph Visualization Software
#' 'yEd' which is freely available from yWorks GmbH
#' (http://www.yworks.com/products/yed) for all major platforms.
#' \cr\cr
#' The colors of the nodes can be defined by the parameters
#' \code{color_by} and \code{name}.
#' Please note that the trajectory landmarks are indicated by setting
#' \code{color_by='phenoName'} and {name='landmark'}. States can be indicated
#' by \code{color_by='phenoName'} and {name='state'}.
#' \cr\cr
#' If a layout is already present in the provided \code{CellTrailsSet}
#' object, the samples' coordinates will be listed in the graphml file.
#' \cr \cr
#' \emph{Diagnostic messages}
#' \cr \cr
#' An error is thrown if the trajectory has not been computed yet; function
#' \code{fitTrajectory} needs to be called first. Feature names and phenotype
#' names get checked and will throw an error if not contained in the dataset.
#' Please note, the parameter \code{name} is case-sensitive.
#' @seealso \code{fitTrajectory} \code{featureNames} \code{phenoNames}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' \dontrun{
#' # Export trajectory graph structure to graphml
#' # Color nodes by gene expression (e.g, feature_10)
#' write.ygraphml(sce, file="yourFilePath", color_by="featureName",
#'               name="feature_10")
#'
#' # Color nodes by metadata (e.g., state) and
#' # label nodes by the (simulated) age of each sample
#' write.ygraphml(sce, file="yourFilePath", color_by="phenoName",
#'               name="state", node_label="age")
#'
#' # Color and label nodes by landmark type and id
#' write.ygraphml(sce, file="yourFilePath", color_by="phenoName",
#'               name="landmark", node_label="landmark")}
#' @docType methods
#' @aliases write.ygraphml,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("write.ygraphml", function(sce, file,
                                      color_by=c("phenoName", "featureName"),
                                      name, node_label="state")
  standardGeneric("write.ygraphml"))
setMethod("write.ygraphml", "SingleCellExperiment", function(sce, file,
                                                             color_by,
                                                             name, node_label){
  #Pre-flight check
  if(is.null(.trajGraph(sce))) {
    stop("Trajectory fitting information not found. Please, call function ",
         "'fitTrajectory' first.")
  }
  node_label <- tolower(node_label)
  node_label <- node_label[1]
  if(length(node_label) > 1) {
    warning("Provided more than one node_label. First label was selected.")
    node_label <- node_label[1]
  }
  .phenoNameExists(sce, pheno_name=node_label)

  color_by <- color_by[1]
  col_params <- .validatePlotParams(x=sce, color_by=color_by, name=name)
  if(col_params$color_by == "phenoname" & col_params$name == "landmark") {
    col_params$values <- as.character(col_params$values)
    col_params$values <- substr(col_params$values, 1, 1)
    col_params$values <- factor(col_params$values)
  }
  lbl_values <- as.character(.pheno(sce, node_label))

  #Run
  g <- .trajGraph(sce) #graph
  col_values <- col_params$values[.useSample(sce)] #filter by traj samples
  lbl_values <- lbl_values[.useSample(sce)]
  X <- trajLayout(sce)
  rownames(X) <- colnames(sce)
  shapes <- .trajLandmark(sce, type="shape")[.useSample(sce)]
  .write_ygraphml_def(g=g, X=X, shapes=shapes, file=file,
                      col_values=col_values, lbl_values=lbl_values)
})

#' Fit expression dynamic
#'
#' Fits feature expression as a function of pseudotime along a defined trail.
#' @param sce A \code{SingleCellExperiment} object
#' @param feature_name Name of feature
#' @param trail_name Name of trail
#' @return An object of type \code{list} with the following components
#' \describe{
#'   \item{\code{pseudotime}}{The pseudotime along the trail}
#'   \item{\code{expression}}{The fitted expression values for
#'   each value of pseudotime}
#'   \item{\code{gam}}{A object of class \code{gamObject}}
#' }
#' @details A trail is an induced subgraph of the trajectory graph. A
#' trajectory graph is composed of samples (nodes) that are connected
#' (by weighted edges) if they are chronologically related. A trail has to be
#' defined by the user using \code{addTrail}. A pseudotime vector is extracted
#' by computing the geodesic distance for each sample from the trail's start
#' node. To infer the expression level of a feature as a function of
#' pseudotime, CellTrails used generalized additive models with a single
#' smoothing term with four basis dimensions. Here, for each feature CellTrails
#' introduces prior weights for each observation to lower the confounding
#' effect of drop-outs to the maximum-likelihood-based fitting process as
#' follows. Each non-detect of feature
#' \emph{j} in state \emph{h} is weighted by the relative fraction of
#' non-detects of feature \emph{j} in state \emph{h}; detected values are
#' always assigned weight = 1.
#' @seealso \code{addTrail} \code{gamObject}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Fit dynamic
#' fit <- fitDynamic(exSCE, feature_name="feature_3", trail_name="Tr1")
#'
#' summary(fit)
#' @docType methods
#' @aliases fitDynamic,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("fitDynamic", function(sce, feature_name,
                                  trail_name)
  standardGeneric("fitDynamic"))
setMethod("fitDynamic", "SingleCellExperiment", function(sce, feature_name,
                                                         trail_name){
  #Pre-flight check
  .featureNameExists(sce, feature_name)
  .trailNameExists(sce, trail_name)
  trail <- trails(sce)[[trail_name]]
  trail_samples <- !is.na(trail)
  trail_ptime <- trail[trail_samples]
  fit <- .fitDynamic_def(x=trail_ptime / max(trail_ptime),
                         y=.exprs(sce[feature_name, trail_samples])[1, ],
                         z=states(sce)[trail_samples],
                         k=5) #fixed to 5
  names(fit) <- c("pseudotime", "expression", "gam")
  fit
})

#' Differential trail expression analysis
#'
#' Comparison of feature expression dynamic between two trails.
#' @param sce A \code{SingleCellExperiment} object
#' @param feature_names Name of feature; can be multiple names
#' @param trail_names Name of trails
#' @param score Score type; one of \{"rmsd", "tad", "abc", "cor"\}
#' @return Numeric value
#' @details Genes have non-uniform expression rates and each trail
#' has a distinct set of upregulated genes, but also contains unequal
#' numbers of cells. Because pseudotime is based on transcriptional change,
#' its axis may be distorted, leading to stretched or compressed sections of
#' longitudinal expression data that make comparison of trails challenging.
#' To align different trails, despite these differences, CellTrails employs a
#' dynamic programming based algorithm that has long been known in speech
#' recognition, called dynamic time warping (Sakoe and Chiba, 1978). RNA
#' expression rates are modeled analogous to speaking rates
#' (Aach and Church, 2001);  the latter accounts for innate non-linear
#' variation in the length of individual phonemes (i.e., states) resulting in
#' stretching and shrinking of word (i.e., trail) segments. This allows the
#' computation of inter-trail alignment warps of individual expression time
#' series that are similar but locally out of phase.
#' \cr \cr
#' Univariate pairwise alignments are
#' computed resulting in one warp per feature and per trail set. Similar to a
#' (global) pairwise protein sequence alignment, monotonicity
#' (i.e., no time loops) and continuity (i.e., no time leaps) constraints have
#' to be imposed on the warping function to preserve temporal sequence ordering.
#' To find the optimal warp, a recursion rule is applied which selects the
#' local minimum of three moves through a dynamic programming matrix:
#' suppose that query snapshot \emph{g} and reference snapshot \emph{h}
#' have already been aligned, then the alignment of \emph{h}+1 with
#' \emph{g}+1 is a (unit slope) diagonal move, \emph{h} with
#' \emph{g}+1 denotes an expansion by repetition of \emph{h},
#' and \emph{h}+2 with \emph{g}+1 contracts the query by dropping \emph{h}+1.
#' \cr \cr
#' The overall dissimilarity between two aligned expression time series
#' \emph{x} and \emph{y}
#' of length \emph{n} is estimated by either the root-mean-square deviation
#' \eqn{RMSD(x, y) = \sqrt(\sum(x - y)^2/n)}, the total aboslute deviation
#' \eqn{TAD(x, y) = \sum(|x-y|)},
#' the area between the aligned dynamic curves (\code{ABC}), or Pearson's
#' correlation coefficient (\code{cor}) over all aligned elements.
#' @references Sakoe, H., and Chiba, S. (1978). Dynamic programming algorithm
#' optimization for spoken word recognition. IEEE Transactions on Acoustics,
#' Speech, and Signaling Processing 26, 43-49.
#' @references Aach, J., and Church, G.M. (2001). Aligning gene expression
#' time series with time warping algorithms. Bioinformatics 17, 495-508.
#' @seealso \code{dtw}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Differential expression between trails
#' contrastTrailExpr(exSCE, feature_name=c("feature_1", "feature_10"),
#'                  trail_names=c("Tr1", "Tr2"), score="rmsd")
#' @docType methods
#' @aliases contrastTrailExpr,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("contrastTrailExpr", function(sce,
                                         feature_names=featureNames(sce),
                                         trail_names, score="rmsd")
  standardGeneric("contrastTrailExpr"))
setMethod("contrastTrailExpr", "SingleCellExperiment",
          function(sce, feature_names, trail_names, score){
  #Pre-flight check
  .featureNameExists(sce, feature_names)
  .trailNameExists(sce, trail_names)
  if(length(trail_names) > 2) {
    warning("Provided more than two trail names. Only the first ",
            "two will be used.")
    trail_names <- trail_names[seq_len(2)]
  }
  if(length(trail_names) < 2) {
    stop("Please, provide two trail names.")
  }
  score <- tolower(score)
  if(!score %in% c("rmsd", "td", "abc", "cor")) {
    stop("Score method unknown. Please, select one of {'rmsd', ",
         "'td', 'abc', 'cor'}.")
  }

  trs <- trails(sce)[, trail_names]
  feature_expr <- .exprs(sce[feature_names, ])

  apply(feature_expr, 1L,
         function(x) .contrastExprTrail_def(ptime_1=trs[,1],
                                            ptime_2=trs[,2],
                                            feature_expr=x,
                                            sts=states(sce),
                                            trail_names=trail_names,
                                            score=score)[[score]])
})

###############################################################################
### Plots
###############################################################################
#' Visualize the number of samples per state
#'
#' Shows barplot of state size distribution
#' @param sce A \code{SingleCellExperiment} object
#' @return A \code{ggplot} object
#' @details Barplot showing the absolute number of samples per state.
#' @seealso \code{ggplot} \code{states}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' plotStateSize(exSCE)
#' @docType methods
#' @aliases plotStateSize,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("plotStateSize", function(sce)
  standardGeneric("plotStateSize"))
setMethod("plotStateSize", "SingleCellExperiment", function(sce){
  if(is.null(states(sce))) {
    stop("No state definition found. Please, assign states first.")
  }
  .plotStateSize_def(states(sce))
})

#' Visualize feature expression distribution per state
#'
#' Violin plots showing the expression distribution of a feature per state.
#' @param sce A \code{SingleCellExperiment} object
#' @param feature_name The name of the feature to be visualized
#' @return A \code{ggplot} object
#' @details Each data point displays the feature’s expression value in a
#' single sample. A violine plot shows the density (mirrored on the y-axis) of
#' the expression distribution per sample.
#' @seealso \code{ggplot} \code{states}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' plotStateExpression(exSCE, feature_name="feature_1")
#' @docType methods
#' @aliases plotStateExpression,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("plotStateExpression", function(sce, feature_name)
  standardGeneric("plotStateExpression"))
setMethod("plotStateExpression", "SingleCellExperiment",
          function(sce, feature_name){
  if(is.null(states(sce))) {
    stop("No state definition found. Please, assign states first.")
  }
  .featureNameExists(sce, feature_name)
  .plotStateExpression_def(x=.exprs(sce[feature_name, ])[1, ], sts=states(sce),
                           label=feature_name)
})

#' Visualize the learned manifold
#'
#' Method visualizes an approximation of the manifold in the latent space
#' in two dimensions.
#' @param sce A \code{SingleCellExperiment} object
#' @param color_by Indicates if nodes are colorized by a feature expression
#' ('featureName') or phenotype label ('phenoName')
#' @param name A character string specifying the featureName or phenoName
#' @param seed Seed for tSNE computation (default: 1101)
#' @param perplexity Perplexity parameter for tSNE computation (default: 30)
#' @param only_plot Indicates if only plot should be shown or the tSNE result
#' should be returned (default: FALSE)
#' @return A \code{ggplot} object or a list with the tSNE results and the plot
#' @details Visualizes the learned lower-dimensional manifold in two dimensions
#' using an approximation obtained by Barnes-Hut implementation of
#' t-Distributed Stochastic Neighbor Embedding
#' (tSNE; van der Maaten and Hinton 2008). Each point in this plot represents
#' a sample. Points can be colorized according
#' to feature expression or experimental metadata. The points' coloration can
#' be defined via the attributes \code{feature_name} and \code{name},
#' respectively. A previously computed tSNE visualization will be reused
#' (if the parameter settings are identical). The parameters \code{tsne_seed}
#' and \code{tsne_perplexity} are used for the tSNE calculation.
#' @references van der Maaten, L.J.P. & Hinton, G.E., 2008. Visualizing
#' High-Dimensional Data Using t-SNE. Journal of Machine Learning Research,
#' 9, pp.2579-2605.
#' @seealso \code{Rtsne} \code{latentSpace}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' plotManifold(exSCE, color_by="featureName", name="feature_10", only_plot=TRUE)
#' plotManifold(exSCE, color_by="phenoName", name="age", only_plot=TRUE)
#' @docType methods
#' @aliases plotManifold,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("plotManifold", function(sce,
                                    color_by=c("phenoName", "featureName"),
                                    name, perplexity=30, seed=1101,
                                    only_plot=FALSE)
  standardGeneric("plotManifold"))
setMethod("plotManifold", "SingleCellExperiment",
          function(sce, color_by=c("phenoName", "featureName"),
                   name, perplexity, seed, only_plot){
  #Fetch plot color data
  dat <- .validatePlotParams(sce, color_by=color_by, name=name)
  #Fetch tSNE data
  if(is.null(CellTrails::latentSpace(sce))) {
    stop("Please, define CellTrails' latent space first (see function",
         "latentSpace<-).")
  }
  tsne_params <- c("seed"=seed, "perplexity"=perplexity)
  if(is.null(latentSpaceSNE(sce)) |
     any(!tsne_params == .tsneParams(sce))) { #recalc
    message("Calculating 2D approximation of CellTrails manifold...")
    X <- .bhtsne(latentSpace(sce), perplexity=perplexity, seed=seed)$Y
  } else {
    X <- latentSpaceSNE(sce)
  }
  gp <- .plotManifold_def(X=X, y=dat$values, name=dat$name,
                          axis_label = "CellTrails tSNE",
                          type="raw", setND=(dat$color_by=="featurename"))
  print(gp)
  if(!only_plot) {
    list(tsne=list(X=X, seed=seed, perplexity=perplexity), plot=gp)
  }
})

#' Visualize state trajectory graph
#'
#' Method visualizes the state-to-state relations delineating the
#' trajectory backbone.
#' @param sce A \code{SingleCellExperiment} object
#' @param color_by Indicates if nodes are colorized by a feature expression
#' ('featureName') or phenotype label ('phenoName')
#' @param name A character string specifying the featureName or phenoName
#' @param component Component of trajectory graph that should
#' be shown (integer value)
#' @param point_size Adjusts the point size of the data points shown
#' @param label_offset Adjusts the offset of the data point labels
#' @param seed Seed for layout computation (default: 1101)
#' @return A \code{ggplot} object
#' @details Shows a single tree component of the computed trajectory graph.
#' Each point in this plot represents a state and can be colorized
#' according to feature expression (mean expression per state) or experimental
#' metadata (arithmetic mean or percentage distribution of categorial values).
#' The component is defined by parameter \code{component}. If the trajectory
#' graph contains only a single component, then this parameter can be left
#' undefined. The points' coloration can be defined via the attributes
#' \code{feature_name} or \code{pheno_type}. Missing sample lables are
#' recovered using nearest neighbor learning.
#' @seealso \code{connectStates}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' plotStateTrajectory(exSCE, color_by="phenoName", name="age", component=1)
#' plotStateTrajectory(exSCE, color_by="featureName", name="feature_1",
#'                     component=1)
#' @docType methods
#' @aliases plotStateTrajectory,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("plotStateTrajectory",
           function(sce, color_by=c("phenoName", "featureName"),
                    name, seed=1101, component=NULL, point_size=3,
                    label_offset=2)
             standardGeneric("plotStateTrajectory"))
setMethod("plotStateTrajectory", "SingleCellExperiment",
          function(sce, color_by=c("phenoName", "featureName"),
                   name, seed, component, point_size, label_offset){
  #Pre-flight checks
  if(is.null(.spanForest(sce))) {
    stop("Please, compute the state trajectory graph first.")
  }
  if(is.null(component) & length(.spanForest(sce)) > 1) {
    stop("Please, specify the component.")
  }
  if(is.null(component)) {
    component <- 1
  }
  component <- abs(component)
  if(component > length(.spanForest(sce))) {
    stop("Unknown component selected. Please, make sure that the ",
         "correct component number was selected.")
  }
  # Plotting data
  dat <- .validatePlotParams(sce, color_by=color_by, name=name)
  g <- .spanForest(sce)[[component]]
  y <- .nn_impute(y=dat$values, D=as.matrix(stats::dist(latentSpace(sce))))
  .plotStateTrajectory_def(g=g, y=y, name=dat$name, all_sts=states(sce),
                           point_size=point_size, label_offset=label_offset,
                           seed=seed, setND=(dat$color_by == "featurename"))
})

#' Visualize trajectory fit residuals
#'
#' Method visualizes the fitting residuals along the
#' trajectory backbone.
#' @param sce A \code{SingleCellExperiment} object
#' @return A \code{ggplot} object
#' @details Shows a single tree component of the computed trajectory graph.
#' Each point in this plot represents a state and can be colorized
#' according to feature expression (mean expression per state) or experimental
#' metadata (arithmetic mean or percentage distribution of categorial values).
#' The component is defined by parameter \code{component}. If the trajectory
#' graph contains only a single component, then this parameter can be left
#' undefined. The points' coloration can be defined via the attributes
#' \code{feature_name} or \code{pheno_type}. Missing sample lables are
#' recovered using nearest neighbor learning.
#' @seealso \code{fitTrajectory} \code{trajResiduals}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' plotTrajectoryFit(exSCE)
#' @docType methods
#' @aliases plotTrajectoryFit,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("plotTrajectoryFit", function(sce)
             standardGeneric("plotTrajectoryFit"))
setMethod("plotTrajectoryFit", "SingleCellExperiment", function(sce){
  #Pre-flight checks
  if(is.null(trajResiduals(sce))) {
    stop("Please, compute the trajectory fit first.")
  }
  # Plotting data
  g <- .trajGraph(sce)
  x <- trajResiduals(sce)[.useSample(sce)]
  sts <- states(sce)[.useSample(sce)]
  .plot_trajectoryFit(x=x, g=g, sts=sts, factor=7, rev=FALSE)
})

#' Visualize expression maps
#'
#' Method visualizes topographical expression maps
#' in two dimensions.
#' @param sce A \code{SingleCellExperiment} object
#' @param color_by Indicates if nodes are colorized by a feature expression
#' @param name A character string specifying the featureName or phenoName
#' @param type Type of map; one of {"raw","surface.fit","surface.se"}
#' @param samples_only If only individual samples should be colorized rather
#' than the whole surface (default: FALSE)
#' @return A \code{ggplot} object
#' @details Two-dimensional visualization of the trajectory. The red line
#' representsthe trajectory and individual points denote samples. This plot
#' type can either show thetopography of a given feature’s expression landscape
#' or colorizes individual samples by a metadata label. The feature is selected
#' by setting the parameter \code{color_type} and the respecitve \code{name}.
#' To show feature expression, a surface is fitted using isotropic (i.e., same
#' parameters for both map dimensions) thin-plate spline smoothing in
#' \code{gam}. It gives an overview of expression dynamics along all
#' branches of the trajectory. The parameter \code{type} defines if either the
#' raw/original expression data shoud be shown, the full fitted expression
#' surface should be shown (\code{type="surface.fit"}) or the standard error
#' of the surface prediction (\code{type="surface.se"}), or the expression
#' values of single samples only (\code{type="surface.fit"}
#' and \code{only_samples=TRUE}).
#' \cr\cr
#' To show all landmarks on the map, please use the parameters
#' \code{color_by="phenoName"} and \code{name="landmark"}.
#' @seealso \code{gam}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Plot landmarks
#' plotMap(exSCE, color_by="phenoName", name="landmark")
#'
#' # Plot phenotype
#' plotMap(exSCE, color_by="phenoName", name="age")
#'
#' # Plot feature expression map
#' plotMap(exSCE, color_by="featureName", name="feature_10", type="surface.fit")
#' plotMap(exSCE, color_by="featureName", name="feature_10", type="surface.fit",
#'         samples_only=TRUE)
#'
#' #Plot surface fit standard errors
#' plotMap(exSCE, color_by="featureName", name="feature_10", type="surface.se")
#' @docType methods
#' @aliases plotMap,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("plotMap", function(sce,
                               color_by=c("phenoName", "featureName"),
                               name,
                               type=c("surface.fit", "surface.se", "raw"),
                               samples_only=FALSE)
  standardGeneric("plotMap"))
setMethod("plotMap", "SingleCellExperiment", function(sce, color_by, name,
                                                      type, samples_only){
  #Pre-flight test
  if(is.null(trajLayout(sce))) {
    stop("No graph layout detected. Please, set a trajectory ",
         "graph layout first.")
  }
  type <- tolower(type)[1]
  if(!type %in% c("surface.fit", "surface.se", "raw")) {
    stop("Unknown plot type selected.")
  }

  #Fetch plot color data
  dat <- .validatePlotParams(sce, color_by=color_by, name=name)
  if(dat$color_by=="featurename") {
    weights <- states(sce[, .useSample(sce)])
  } else {
    weights <- NULL
  }
  X <- trajLayout(sce)[.useSample(sce), ]
  rownames(X) <- trajSampleNames(sce)

  if(dat$color_by=="phenoname" & dat$name=="landmark") {
    .plotTrailblazing_def(X=X, g=.trajGraph(sce),
                          ltype=.trajLandmark(sce, "type")[.useSample(sce)],
                          lid=.trajLandmark(sce, "id")[.useSample(sce)])
  } else {
    .plotManifold_def(X=X,
                      g=.trajGraph(sce),
                      weights=weights,
                      y=dat$values[.useSample(sce)],
                      name=dat$name, type=type,
                      axis_label = "CellTrails",
                      samples_only=samples_only,
                      setND=(dat$color_by=="featurename"))
  }
})

#' Visualize single trails
#'
#' Method highlights a single trail on the trajectory map
#' @param sce A \code{SingleCellExperiment} object
#' @param name Name of the trail
#' @return A \code{ggplot} object
#' @details A trail can be defined with the function \code{addTrail} between
#' two landmarks. User-defined landmarks can be set with the function
#' \code{userLandmarks}. This function visualizes the start and endpoints, and
#' the pseudotime of a defined trail along the trajectory. The trail
#' pseudotimes can be directly accessed via the \code{trails}.
#' \cr\cr
#' An error is thrown if the \code{trail_name} is unknown. The function is
#' case-sensitiv. All available trails can be listed by \code{trailNames}.
#' @seealso \code{addTrail} \code{userLandmarks} \code{trailNames}
#' \code{trails}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Plot trail
#' plotTrail(exSCE, name="Tr1")
#' @docType methods
#' @aliases plotTrail,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("plotTrail", function(sce, name)
  standardGeneric("plotTrail"))
setMethod("plotTrail", "SingleCellExperiment", function(sce, name){
  #Pre-flight test
  if(is.null(trajLayout(sce))) {
    stop("No graph layout detected. Please, set a trajectory ",
         "graph layout first.")
  }
  .trailNameExists(sce, name)

  #Fetch data
  X <- trajLayout(sce)[.useSample(sce), ]
  rownames(X) <- trajSampleNames(sce)
  ptime <- trails(sce)[.useSample(sce), name, drop=FALSE]
  .plotTrail_def(X=X, g=.trajGraph(sce), ptime=ptime, name=name)
})

#' Visualize dynamics
#'
#' Shows dynamics of one or multiple features along a given trail
#' @param sce A \code{SingleCellExperiment} object
#' @param feature_name Name of one or multiple features
#' @param trail_name Name of trail
#' @return A \code{ggplot} object
#' @details
#' An error is thrown if the \code{trail_name} or \code{feature_name} are
#' unknown. The function is case-sensitiv. All available trails can be
#' listed by \code{trailNames}, all features with \code{featureNames}.
#' @seealso \code{addTrail} \code{trailNames} \code{featureNames}
#' @examples
#' # Example data
#' data(exSCE)
#'
#' # Plot dynamic of feature_10
#' plotDynamic(exSCE, trail_name="Tr1", feature_name="feature_1")
#' # Plot dynamic of feature_1 and feature_10
#' plotDynamic(exSCE, trail_name="Tr1",
#'             feature_name=c("feature_1", "feature_10"))
#' @docType methods
#' @aliases plotDynamic,SingleCellExperiment-method
#' @export
#' @author Daniel C. Ellwanger
setGeneric("plotDynamic", function(sce, feature_name, trail_name)
  standardGeneric("plotDynamic"))
setMethod("plotDynamic", "SingleCellExperiment",function(sce,
                                                         feature_name,
                                                         trail_name){
  #Pre-flight test
  .featureNameExists(sce, feature_name)
  .trailNameExists(sce, trail_name)
  if(length(trail_name) > 1) {
    trail_name <- trail_name[1]
    warning("Provided more than one trail_name. Dynamic is only shown for ",
            trail_name, ".")
  }

  #Fetch data
  x <- trails(sce[, .useSample(sce)])[[trail_name]]
  Y <- .exprs(sce[feature_name, .useSample(sce)])
  weights <- states(sce[, .useSample(sce)])
  .plotDynamic(x=x, Y=Y, feature_name=feature_name,
               trail_name=trail_name, weights=weights, k=5)
})
