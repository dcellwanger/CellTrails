test_plotMap <- function() {
  data(exSCE)
  dat <- SingleCellExperiment(assay=list(logcounts=logcounts(exSCE)))
  RUnit::checkException(plotMap(dat, #Empty
                                color_by="featureName",
                                name=rownames(dat)[1],
                                type="surface.fit",
                                samples_only=FALSE))
  RUnit::checkException(plotMap(exSCE,
                                color_by="A", #wrong input
                                name=rownames(exSCE)[1],
                                type="surface.fit",
                                samples_only=FALSE))
  RUnit::checkException(plotMap(exSCE,
                                color_by="featureName",
                                name="A", #wrong input
                                type="surface.fit",
                                samples_only=FALSE))
  RUnit::checkException(plotMap(exSCE,
                                color_by="featureName",
                                name=rownames(exSCE)[1],
                                type="A", #wrong input
                                samples_only=FALSE))
  ggp <- plotMap(exSCE,
                 color_by="featureName",
                 name=rownames(exSCE)[1],
                 type="surface.fit",
                 samples_only=TRUE)
  RUnit::checkTrue(length(ggp$layers) > 1)
  ggp <- plotMap(exSCE,
                 color_by="featureName",
                 name=rownames(exSCE)[1],
                 type="surface.se",
                 samples_only=TRUE)
  RUnit::checkTrue(length(ggp$layers) > 1)
  ggp <- plotMap(exSCE,
                 color_by="featureName",
                 name=rownames(exSCE)[1],
                 type="surface.fit",
                 samples_only=TRUE)
  ggp <- plotMap(exSCE,
                 color_by="featureName",
                 name=rownames(exSCE)[1],
                 type="raw",
                 samples_only=FALSE) #ignore (always 'samples_only' for 'raw')
  RUnit::checkTrue(length(ggp$layers) > 1)
  ggp <- plotMap(exSCE,
                 color_by="phenoName",
                 name=phenoNames(exSCE)[1],
                 type="surface.fit", #ignore (always 'raw' for discrete input)
                 samples_only=FALSE)
  RUnit::checkTrue(length(ggp$layers) > 1)
}
