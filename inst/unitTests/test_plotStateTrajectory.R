test_plotStateTrajectory <- function() {
  data(exSCE)
  dat <- SingleCellExperiment(assay=list(logcounts=logcounts(exSCE)))
  RUnit::checkException(plotStateTrajectory(dat, #Empty
                                            color_by="featureName",
                                            name=rownames(dat)[1],
                                            component=1))
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
                                            component=10)) #wrong input
  exSCE <- connectStates(exSCE, l=10) #has two graph components
  gp <- plotStateTrajectory(exSCE,
                            color_by="featureName",
                            name=rownames(exSCE)[1],
                            component=1)
  RUnit::checkException(stateTrajLayout(exSCE) <- gp)
  RUnit::checkException(stateTrajLayout(exSCE) <- list())
  l1 <- length(gp$layers)
  RUnit::checkTrue(l1 > 1) #points
  gp <- plotStateTrajectory(exSCE,
                            color_by="phenoName",
                            name=phenoNames(exSCE)[1],
                            component=1)
  l2 <- length(gp$layers)
  RUnit::checkTrue(l2 > l1) #pie charts
  gp <- plotStateTrajectory(exSCE,
                            color_by="phenoName",
                            name=phenoNames(exSCE)[1],
                            component=1, recalculate=TRUE)
  stateTrajLayout(exSCE) <- gp

  #Second component
  gp <- plotStateTrajectory(exSCE,
                            color_by="phenoName",
                            name=phenoNames(exSCE)[1],
                            component=2)
  RUnit::checkException(stateTrajLayout(exSCE) <- gp)

  #Select second component
  exSCE <- selectTrajectory(exSCE, component=2)
  gp <- plotStateTrajectory(exSCE,
                            color_by="phenoName",
                            name=phenoNames(exSCE)[1])
}
