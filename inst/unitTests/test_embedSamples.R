test_embedSamples <- function(){
  data(exSCE)
  dat <- SingleCellExperiment(assay=list(logcounts=logcounts(exSCE)))

  #TEST
  se <- embedSamples(dat)
  RUnit::checkEquals(names(se), c("components", "eigenvalues"))
  RUnit::checkEquals(dim(se$components), rep(ncol(dat), 2))
  RUnit::checkEquals(length(se$eigenvalues), ncol(dat))
  RUnit::checkTrue(any(!is.na(se$components)))
  RUnit::checkTrue(any(!is.na(se$eigenvalues)))
}
