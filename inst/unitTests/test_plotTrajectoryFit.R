test_plotTrajectoryFit <- function() {
  dat <- CellTrails:::.exDat()
  RUnit::checkException(plotTrajectoryFit(dat)) #Empty

  data(exSCE)
  ggp <- plotTrajectoryFit(exSCE)
  RUnit::checkTrue(length(ggp$layers) > 1)
}
