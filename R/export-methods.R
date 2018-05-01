#' DEF: Export trajectory graph
#'
#' For details see \code{write.ygraphml}
#' @importFrom igraph V E vcount ends
#' @importFrom grDevices colorRampPalette
#' @keywords internal
#' @author Daniel C. Ellwanger
.write_ygraphml_def <- function(x, file, pheno_type=NULL,
                                feature_name=NULL, label) {
  #write.graph(x@trajectory$traj, file = file, format = c("graphml"))

  if(!is.null(pheno_type) & !is.null(feature_name)) {
    stop("Both parameters, pheno_type and feature_name, were set. ",
         "Please, choose either a feature name or a pheno type.")
  }

  #Get graph
  g <- x@trajectory$traj

  #Sample name arg
  snames <- V(g)$sampleName

  if(is.null(pheno_type) & is.null(feature_name)) {
    pheno_type <- "STATE"
  }

  cols <- sts <- rep(NA, vcount(g))

  if(!is.null(pheno_type)) {
    if(toupper(pheno_type) == "LANDMARK") {
      #Color and states arg
      sts <- x@trajectory$blaze$type
      cols <- .color_hue(length(levels(sts)))
      cols <- cols[sts]
      #cols[is.na(cols)] <- "#FFFFFF"
      sts <- x@trajectory$blaze$id
      sts[is.na(sts)] <- ""
    } else {
      pheno_type <- ifelse(toupper(pheno_type) == "STATE", "STATE",
                           pheno_type)
      if(!pheno_type %in% varLabels(x)) {
        stop("Variable label '", pheno_type, "' not found. ",
             "Please, check correct spelling. ",
             "To see all available pheno_types call 'varLabels(object)'.")
      }
      #Color and states arg
      sts <- pData(x)[snames, pheno_type]
      cols <- .color_hue(length(levels(sts)))
      cols <- cols[sts]
      #cols[is.na(cols)] <- "#FFFFFF"
      sts <- as.character(sts)
      sts[is.na(sts)] <- ""
    }
  } else {
    .checkFeatureNameExists(x, feature_name)
    fit <- .fit_surface(x, feature_name=feature_name, npoints=1)
    fval <- fit$fit$fitted.values
    brks <- pretty(range(fval, na.rm=TRUE, finite=TRUE), 10)
    equalSpace <- cut(fval, breaks=brks, include.lowest=TRUE)
    cols <- colorRampPalette(
      c("gray95", viridis(3)[2:3]))(length(brks) - 1)[equalSpace]
    sts[seq_along(sts)] <- ""
  }

  if(!is.null(label)) {
    if(toupper(label) == "NAME") {
      sts <- snames
    } else if(toupper(label) == "STATE") {
      sts <- states(x)
    } else {
      if(!label %in% varLabels(x)) {
        stop("Variable label '", pheno_type, "' not found.",
             "Please, check correct spelling for parameter label. ",
             "To see all available pheno_types call 'varLabels(object)'.")
      } else {
        sts <- as.character(sts)
        sts[is.na(sts)] <- ""
      }
    }
  } else {
    label <- snames
  }

  #Layout info
  tlayout <- trajectoryLayout(x)[snames, ]

  #Head
  write(c('<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
               paste('<graphml xmlns="http://graphml.graphdrawing.org/xmlns"',
                     'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"',
                     'xmlns:y="http://www.yworks.com/xml/graphml"',
                     'xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns',
                     'http://www.yworks.com/xml/schema/graphml/1.1/ygraphml.xsd',
                     'http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">',
                     collapse = " "),
               "  <!-- Created by CellTrails -->",
               "  <!-- Compatible w/ yEd 3.14 -->",
               '  <key for="port" id="d0" yfiles.type="portgraphics"/>',
               '  <key for="port" id="d1" yfiles.type="portgeometry"/>',
               '  <key for="port" id="d2" yfiles.type="portuserdata"/>',
               '  <key for="node" id="d3" attr.name="sampleName" attr.type="string" />',
               '  <key for="node" id="d4" yfiles.type="nodegraphics"/>',
               '  <key for="edge" id="d5" attr.name="weight" attr.type="double" />',
               '  <key for="edge" id="d6" yfiles.type="edgegraphics"/>',
               '  <key for="graphml" id="d7" yfiles.type="resources"/>',
               '  <key for="graph" id="d8" attr.name="Description" attr.type="string"/>',
               '  <graph edgedefault="undirected" id="G">',
               '  <data key="d8"/>'),
             file = file, append = FALSE)

  #Add nodes
  for(i in seq_along(V(g))) {
    nx <- ny <- 0
    if(!is.null(tlayout)) {
      nx <- tlayout[snames[i], 1]
      ny <- -tlayout[snames[i], 2]
    }

    lw <- lx <- ly <- 0
    s <- nchar(as.character(sts[i])) #as.numeric(substring(sts[i], 2))
    if(s < 10) {
      lw <- 18.05078125
      lx <- 5.974609375
      ly <- 5.93359375
    } else if(s < 100) {
      lw <- 25.638671875
      lx <- 2.1806640625
      ly <- 5.93359375
    } else {
      lw <- 33.2265625
      lx <- -1.61328125
      ly <- 5.93359375
    }

    nodeshape <- x@trajectory$blaze$shape[i]
    nodecol <- ifelse(is.na(cols[i]),
                      '          <y:Fill hasColor="false" transparent="false"/>',
                      paste0('          <y:Fill color="', cols[i], '" transparent="false"/>'))

    write(c(
      paste0('    <node id="n', V(g)[i], '">'),
      paste0('      <data key="d3"><![CDATA[@@@', snames[i], '@@@]]></data>'), #tagged with @@@
      '      <data key="d4">',
      '        <y:ShapeNode>',
      paste0('          <y:Geometry height="30.0" width="30.0" x="', nx, '" y="', ny, '"/>'),
      nodecol,
      '          <y:BorderStyle color="#000000" type="line" width="1.0"/>',
      paste0('          <y:NodeLabel alignment="center" autoSizePolicy="content" ',
             'fontFamily="Dialog" fontSize="12" fontStyle="plain" hasBackgroundColor="false" ',
             'hasLineColor="false" height="18.1328125" modelName="internal" modelPosition="c" ',
             'textColor="#000000" visible="true" width="', lw, '" x="', lx, '" y="', ly, '">',
             sts[i], '</y:NodeLabel>'),
      paste0('          <y:Shape type="', nodeshape, '"/>'),
      '        </y:ShapeNode>',
      '      </data>',
      '    </node>'),
    file = file, append = TRUE)
  }

  #Add edges
  elist <- ends(g, E(g))
  w <- E(g)$weight
  for(i in seq_len(nrow(elist))) {
    write(c(paste0('    <edge id="e', i, '" source="n', elist[i, 1], '" target="n', elist[i, 2], '">'),
            paste0('      <data key="d5">', w[i] ,'</data>'),
            '      <data key="d6">',
            '        <y:PolyLineEdge>',
            '          <y:Path sx="0.0" sy="0.0" tx="0.0" ty="0.0"/>',
            '          <y:LineStyle color="#000000" type="line" width="1.0"/>',
            '          <y:Arrows source="none" target="none"/>',
            '          <y:BendStyle smoothed="false"/>',
            '        </y:PolyLineEdge>',
            '      </data>',
            '    </edge>'
            ),
          file = file, append = TRUE)

  }

  #Tail
  write(c('  </graph>',
               '  <data key="d7">',
               '    <y:Resources/>',
               '  </data>',
               '</graphml>'),
             file = file, append = TRUE)
}
