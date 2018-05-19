test_connectStates <- function(){
  data(exSCE)
  dat <- SingleCellExperiment(assay=list(logcounts=logcounts(exSCE)))

  se <- embedSamples(dat)
  d <- findSpectrum(se$eigenvalues, frac=30)
  latentSpace(dat) <- se$components[, d]
  states(dat) <- findStates(dat,
                            min_size=0.01,
                            min_feat=2,
                            max_pval=1e-4,
                            min_fc=2)

  #TEST
  RUnit::checkException(connectStates(exSCE, l=-1)) #wrong l
  RUnit::checkException(connectStates(exSCE, l=0)) #wrong l
  dat <- connectStates(dat, l=1) #maximum
  RUnit::checkEqualsNumeric(length(trajComponents(dat)), 5)
  dat <- connectStates(dat, l=30) #minimum
  RUnit::checkEqualsNumeric(length(trajComponents(dat)), 1)
}
