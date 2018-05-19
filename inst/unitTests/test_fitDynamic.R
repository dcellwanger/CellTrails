test_fitDynamic <- function(){
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
  dat <- connectStates(dat, l=30)
  dat <- fitTrajectory(dat)
  dat <- addTrail(dat, from = "H1", to = "H2", name = "Tr1")
  dat <- addTrail(dat, from = "H2", to = "H4", name = "Tr2")

  # TEST
  RUnit::checkException(fitDynamic(dat, feature_name="feature_1",
                                   trail_name="AAA"))
  fit <- fitDynamic(dat, feature_name="feature_1",
                    trail_name="Tr1")
  RUnit::checkEquals(range(fit$pseudotime), c(0, 1))
  RUnit::checkEquals(length(fit$pseudotime), length(fit$expression))
  RUnit::checkEquals(class(fit$gam)[1], "gam")
}
