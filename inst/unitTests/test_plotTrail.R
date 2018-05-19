test_plotTrail <- function(){
  data(exSCE)
  dat <- SingleCellExperiment(assay=list(logcounts=logcounts(exSCE)))
  RUnit::checkException(plotTrail(dat, name="Tr1")) #Empty
  RUnit::checkException(plotTrail(exSCE, name="A")) #wrong input
  ggp <- plotTrail(exSCE, name=trailNames(exSCE)[1])
  RUnit::checkTrue(length(ggp$layers) > 1)
}
