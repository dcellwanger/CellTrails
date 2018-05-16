#' Reads trajectory graph layout
#'
#' Reads ygraphml file containing the trajectory graph's layout
#' @param file A character string naming a file
#' @return An \code{data.frame} with coordinates of data points and
#' visualization metadata
#' @details To visualize the trajectory graph, a proper graph layout has
#' to be computed. Ideally, edges should not cross and nodes should not
#' overlap. CellTrails enables the export and import of the trajectory
#' graph structure using the graphml file format. This file format can be
#' interpreted by most third-party graph analysis applications, allowing
#' the user to subject the trajectory graph to a wide range of layout
#' algorithms. Please note that the graphml file needs to contain layout
#' information ("<y:Geometry x=... y=... >" entries) as provided by the
#' 'ygraphml' file definition used by the Graph Visualization Software
#' 'yEd' (freely available from yWorks GmbH,
#' http://www.yworks.com/products/yed).
#' @seealso \code{write.ygraphml}
#' @examples
#' # Generate example data
#' sce <- exDat()
#'
#' fn <- system.file("exdata", "exDat.graphml", package="CellTrails")
#' tl <- read.ygraphml(fn)
#' @docType methods
#' @export
#' @author Daniel C. Ellwanger
read.ygraphml <- function(file) {

  ldat <- readLines(file)
  #Sample names
  tmp <- regmatches(ldat, regexpr("@@@.+@@@", ldat))
  snames <- gsub(pattern = "@", replacement = "", x = tmp)

  if(length(snames) == 0){
    stop("Malformed graphml file: No sample name information found: each <node ...> entity
         should have a 'sampleName' attribute.")
  }

  # Geometry
  tmp <- regmatches(ldat, regexpr("<y:Geometry.*>", ldat))
  if(length(tmp) == 0) {
    stop("Malformed graphml file: No geometry information found. Please, make sure that each node
         has x- and y-coordinates: each <node ...> entity should have <y:Geometry ...> attribute
         in the graphml file.")
  }

  x.tmp <- regmatches(tmp, regexpr('x=\"-?[0-9]+\\.?[0-9]*\"', tmp))
  y.tmp <- regmatches(tmp, regexpr('y=\"-?[0-9]+\\.?[0-9]*\"', tmp))
  x1 <- gsub(pattern = "[^-?[:digit:]+\\.?-?[:digit:]+]", replacement = "",
             x = x.tmp, perl = TRUE)
  x2 <- gsub(pattern = "[^-?[:digit:]+\\.?-?[:digit:]+]", replacement = "",
             x = y.tmp, perl = TRUE)

  Y <- data.frame(x1 = as.numeric(x1), x2 = -as.numeric(x2))
  row.names(Y) <- snames

  # Shape by <y:Shape type="ellipse"/>
  tmp <- regmatches(ldat, regexpr("<y:Shape type.*>", ldat))
  if(length(tmp) > 0) {
    shape.tmp <- regmatches(tmp, regexpr('\"[a-zA-Z]+\"', tmp))
    Y$shape <- gsub("\"", replacement = "", x = shape.tmp)
  }
  Y
}
