test_plotManifold <- function() {
  data(exSCE)
  dat <- SingleCellExperiment(assay=list(logcounts=logcounts(exSCE)))
  RUnit::checkException(plotManifold(dat, #Empty
                                     color_by="featureName",
                                     name=rownames(dat)[1]))
  RUnit::checkException(plotManifold(exSCE,
                                     color_by="A", #wrong input
                                     name=rownames(exSCE)[1]))
  RUnit::checkException(plotManifold(exSCE,
                                     color_by="featureName",
                                     name="A")) #wrong input
  gp <- plotManifold(exSCE,
                      color_by="featureName",
                      name=rownames(exSCE)[1])
  RUnit::checkTrue(length(gp$layers) > 1)
  RUnit::checkException(manifold2D(exSCE) <- list())
  RUnit::checkException(manifold2D(exSCE) <- gp) #no recalculation
  gp <- plotManifold(exSCE,
                      color_by="featureName",
                      name=rownames(exSCE)[1],
                      recalculate=TRUE) #recalc
  mold <- manifold2D(exSCE)
  manifold2D(exSCE) <- gp
  mnew <- manifold2D(exSCE)
  RUnit::checkTrue(abs(sum(abs(mold) - abs(mnew))) > 0)
}
