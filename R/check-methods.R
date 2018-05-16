#' Checks if feature exists
#'
#' @param x A \code{SingleCellExperiment} object
#' @param feature_name Name of feature
#' @return logical value
#' @keywords internal
#' @author Daniel C. Ellwanger
.featureNameExists <- function(x, feature_name) {
  test <- feature_name %in% rownames(x)
  if(!all(test)) {
    stop("Feature(s) ' ", paste0(feature_name[!test], " "),
         "' is/are not contained in the data set. ",
         "Please, check for correct spelling (e.g., case sensitivity). ",
         "To see all available feature names call 'featureNames'.")
  }
}

#' Checks if phenotype exists
#'
#' @param x A \code{SingleCellExperiment} object
#' @param pheno_name Name of phenotype
#' @return logical value
#' @keywords internal
#' @author Daniel C. Ellwanger
.phenoNameExists <- function(x, pheno_name) {
  test <- pheno_name %in% phenoNames(x)
  if(!all(test)) {
    stop("Phenotype(s) ' ", paste0(pheno_name[!test], " "),
         "' is/are not contained in the data set. ",
         "Please, check for correct spelling (e.g., case sensitivity). ",
         "To see all available phenotypes call 'phenoNames'.")
  }
}

#' Checks if sample exists
#'
#' @param x A \code{SingleCellExperiment} object
#' @param sample_name Name of sample
#' @return logical value
#' @keywords internal
#' @author Daniel C. Ellwanger
.sampleNameExists <- function(x, sample_name) {
  test <- sample_name %in% sampleNames(x)
  if(!all(test)) {
    stop("Sample(s) ' ", paste0(sample_name[!test], " "),
         "' is/are not contained in the data set. ",
         "Please, check for correct spelling (e.g., case sensitivity). ",
         "To see all available sample names call 'sampleNames'.")
  }
}

#' Checks if trail exists
#'
#' @param x A \code{SingleCellExperiment} object
#' @param trail_name Name of trail
#' @return logical value
#' @keywords internal
#' @author Daniel C. Ellwanger
.trailNameExists <- function(x, trail_name) {
  test <- trail_name %in% trailNames(x)
  if(!all(test)) {
    stop("Trail(s) ' ", paste0(trail_name[!test], " "),
         "' is/are not contained in the data set. ",
         "Please, check for correct spelling (e.g., case sensitivity). ",
         "To see all available trail names call 'trailNames'.")
  }
}

#' Validates parameter settings for plot methods
#'
#' @param x A \code{SingleCellExperiment} object
#' @param color_by Color_by parameter
#' @param name Name parameter
#' @details Parameter \code{color_by} needs to be in
#' {'featureName', 'phenoName'}. Parameter \code{name} needs
#' to be existent in rownames(object) or in
#' {colnames(colData(object)), state, landmark}.
#' If parameters were validated true, this function automatically
#' extracts and returns the values. Character vectors get converted
#' to factors.
#' @return A \code{list} with the parameters and values
#' @keywords internal
#' @author Daniel C. Ellwanger
.validatePlotParams <- function(x, color_by, name) {
  if(length(color_by) != 1) {
    stop("Please, choose 'featureName' or 'phenoName' as a color ",
         "selection criteria (i.e., set parameter 'color_by' properly).")
  }
  # Correct for case sensitivity
  color_by <- tolower(color_by)

  if(!color_by %in% c("featurename", "phenoname")) {
    stop("Unknown selection for parameter 'color_by'. Please, choose either ",
         "'featureName' or 'phenoName'.")
  }
  # Check name
  if(length(name) > 1) {
    warning("Provided more than one name. First name was selected.")
    name <- name[1]
  }

  if(color_by == "featurename") {
    .featureNameExists(x, feature_name=name)
    values <- .exprs(x[name, ])[1, ]
  } else if(color_by == "phenoname") {
    .phenoNameExists(x, pheno_name=name)
    values <- .pheno(x, name)
  }
  if(is.character(values)) {
    values <- factor(values)
  }
  list(color_by=color_by, name=name, values=values)
}
