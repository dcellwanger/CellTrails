test_plotTrajectoryFit <- function() {
  data(exSCE)
  dat <- SingleCellExperiment(assay=list(logcounts=logcounts(exSCE)))
  RUnit::checkException(plotTrajectoryFit(dat)) #Empty
  ggp <- plotTrajectoryFit(exSCE)
  RUnit::checkTrue(length(ggp$layers) > 1)
}
