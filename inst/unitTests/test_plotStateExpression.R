test_plotStateExpression <- function() {
  data(exSCE)
  dat <- SingleCellExperiment(assay=list(logcounts=logcounts(exSCE)))
  RUnit::checkException(plotStateExpression(dat, #Empty
                                            feature_name=rownames(dat)))
  RUnit::checkException(plotStateExpression(exSCE,
                                            feature_name="A")) #wrong input
  ggp <- plotStateExpression(exSCE,
                             feature_name=rownames(exSCE)[1])
  RUnit::checkTrue(length(ggp$layers) > 1)
}
