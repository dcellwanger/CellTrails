#' DEF: Read trajectory graph
#'
#' For details see \code{\link[igraph]{read.ygraphml}}
#' @importFrom utils read.table
#' @keywords internal
#' @author Daniel C. Ellwanger
.read_ygraphml_def <- function(x, file, format="graphml", adjust=TRUE) {
  format <- toupper(format)

  if(format == "GRAPHML") {
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

    x.tmp <- regmatches(tmp, regexpr('x=\"-?[0-9]+\\.?[0-9]*\"', tmp)) #regexpr('x=\"-?[0-9]+\\.?-?[0-9]*\"', tmp)
    y.tmp <- regmatches(tmp, regexpr('y=\"-?[0-9]+\\.?[0-9]*\"', tmp))
    x1 <- gsub(pattern = "[^-?[:digit:]+\\.?-?[:digit:]+]", replacement = "", x = x.tmp, perl = TRUE)
    x2 <- gsub(pattern = "[^-?[:digit:]+\\.?-?[:digit:]+]", replacement = "", x = y.tmp, perl = TRUE)

    Y <- data.frame(D1 = as.numeric(x1), D2 = -as.numeric(x2))
    row.names(Y) <- snames

    # Shape by <y:Shape type="ellipse"/>
    tmp <- regmatches(ldat, regexpr("<y:Shape type.*>", ldat))
    if(length(tmp) > 0) {
      shape.tmp <- regmatches(tmp, regexpr('\"[a-zA-Z]+\"', tmp))
      shape <- gsub("\"", replacement = "", x = shape.tmp)
      f <- !shape == "ellipse"
      x@trajectory$blaze$type[f] <- "U"
      x@trajectory$blaze$id[f] <- paste0("U", seq_len(sum(!shape == "ellipse")))
      x@trajectory$blaze$shape[f] <- shape[f]
    }

    #lines <- grep("y:Geometry", ldat)
    #Y <- matrix(NA, ncol = 2, nrow = length(lines))
    #for(i in seq(lines)){
    #  spl <- strsplit(ldat[lines[i]], split = "=")[[1]]
    #  y1 <- as.numeric(strsplit(strsplit(spl[4], " ")[[1]][1], "\"")[[1]][2])
    #  y2 <- as.numeric(strsplit(strsplit(spl[5], " ")[[1]][1], "\"")[[1]][2])
    #  Y[i, 1] <- y1
    #  Y[i, 2] <- -y2
    #}
  } else if(format == "PLAIN") {
    Y <- read.table(file, row.names = 1, colClasses = c("numeric", "numeric"), sep = "\t")
    colnames(Y) <- c("D1", "D2")
  } else {
    stop("Unknown format ", format, ".")
  }
  #tlayout <- apply(Y, 2, .rescale,  ymin = 0, ymax = 1)
  #data.frame(tlayout)
  #trajectoryLayout(x) <- tlayout
  trajectoryLayout(x, adjust=adjust) <- Y[V(x@trajectory$traj)$sampleName, ]
  x
  #Y
}
