test_plotStateSize <- function() {
  data(exSCE)
  dat <- SingleCellExperiment(assay=list(logcounts=logcounts(exSCE)))
  RUnit::checkException(plotStateSize(dat)) #Empty
  ggp <- plotStateSize(exSCE)
  RUnit::checkTrue(length(ggp$layers) > 0)
}
