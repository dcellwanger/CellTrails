test_plotStateExpression <- function() {
  dat <- CellTrails:::.exDat()
  RUnit::checkException(plotStateExpression(dat, #Empty
                                            feature_name=rownames(dat)))

  data(exSCE)
  RUnit::checkException(plotStateExpression(exSCE,
                                            feature_name="A")) #wrong input
  ggp <- plotStateExpression(exSCE,
                             feature_name=rownames(exSCE)[1])
  RUnit::checkTrue(length(ggp$layers) > 1)
}
