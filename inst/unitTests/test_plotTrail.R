test_plotTrail <- function(){
  dat <- CellTrails:::.exDat()
  RUnit::checkException(plotTrail(dat, name="Tr1")) #Empty

  data(exSCE)
  RUnit::checkException(plotTrail(exSCE, name="A")) #wrong input
  ggp <- plotTrail(exSCE, name=trailNames(exSCE)[1])
  RUnit::checkTrue(length(ggp$layers) > 1)
}
