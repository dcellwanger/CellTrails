test_plotStateSize <- function() {
  dat <- CellTrails:::.exDat()
  RUnit::checkException(plotStateSize(dat)) #Empty

  data(exSCE)
  ggp <- plotStateSize(exSCE)
  RUnit::checkTrue(length(ggp$layers) > 0)
}
