test_plotStateTrajectory <- function() {
  dat <- CellTrails:::.exDat()
  RUnit::checkException(plotStateTrajectory(dat, #Empty
                                            color_by="featureName",
                                            name=rownames(dat)[1],
                                            component=1))

  data(exSCE)
  RUnit::checkException(plotStateTrajectory(exSCE,
                                            color_by="A", #wrong input
                                            name=rownames(exSCE)[1],
                                            component=1))
  RUnit::checkException(plotStateTrajectory(exSCE,
                                            color_by="featureName",
                                            name="A", #wrong input
                                            component=1))
  RUnit::checkException(plotStateTrajectory(exSCE,
                                            color_by="featureName",
                                            name=rownames(exSCE)[1],
                                            component=2)) #wrong input
  ggp <- plotStateTrajectory(exSCE,
                             color_by="featureName",
                             name=rownames(exSCE)[1],
                             component=1)
  l1 <- length(ggp$layers)
  RUnit::checkTrue(l1 > 1) #points
  ggp <- plotStateTrajectory(exSCE,
                             color_by="phenoName",
                             name=phenoNames(exSCE)[1],
                             component=1)
  l2 <- length(ggp$layers)
  RUnit::checkTrue(l2 > l1) #pie charts
}
