test_plotManifold <- function() {
  dat <- CellTrails:::.exDat()
  RUnit::checkException(plotManifold(dat, #Empty
                                     color_by="featureName",
                                     name=rownames(dat)[1],
                                     only_plot=TRUE))

  data(exSCE)
  RUnit::checkException(plotManifold(exSCE,
                                     color_by="A", #wrong input
                                     name=rownames(exSCE)[1],
                                     only_plot=TRUE))
  RUnit::checkException(plotManifold(exSCE,
                                     color_by="featureName",
                                     name="A", #wrong input
                                     only_plot=TRUE))
  res <- plotManifold(exSCE,
                      color_by="featureName",
                      name=rownames(exSCE)[1],
                      only_plot=FALSE)

  RUnit::checkTrue(length(res$plot$layers) > 1)
  latentSpaceSNE(exSCE) <- res
  RUnit::checkEquals(latentSpaceSNE(exSCE), res$tsne$X)

  res <- plotManifold(exSCE,
                      color_by="featureName",
                      name=rownames(exSCE)[1],
                      only_plot=FALSE, seed=101)
  RUnit::checkTrue(!all.equal(latentSpaceSNE(exSCE), res$tsne$X) == TRUE)
}
