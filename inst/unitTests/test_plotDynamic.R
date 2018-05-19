test_plotDynamic <- function() {
  data(exSCE)
  dat <- SingleCellExperiment(assay=list(logcounts=logcounts(exSCE)))
  RUnit::checkException(plotDynamic(dat,#empty
                                    feature_name=rownames(dat),
                                    trail_name=trailNames(dat)[1]))
  RUnit::checkException(plotDynamic(exSCE,
                                    feature_name="A", #wrong input
                                    trail_name=trailNames(exSCE)[1]))
  RUnit::checkException(plotDynamic(exSCE,
                                    feature_name=c(rownames(exSCE)[1], "A"), #!
                                    trail_name=trailNames(exSCE)[1]))
  RUnit::checkException(plotDynamic(exSCE,
                                    feature_name=rownames(exSCE)[1],
                                    trail_name=c("A", "feature_1"))) #wrong input
  obs <- tryCatch(plotDynamic(exSCE,
                              feature_name=rownames(exSCE)[1],
                              trail_name=trailNames(exSCE)), #multiple trails
                  warning=conditionMessage)
  RUnit::checkTrue(grepl("Provided more than one", obs))
  ggp <- plotDynamic(exSCE,
                     feature_name=rownames(exSCE)[1],
                     trail_name=trailNames(exSCE)[1])
  RUnit::checkTrue(length(ggp$layers) > 1)
}
