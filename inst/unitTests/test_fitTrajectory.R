test_fitTrajectory <- function(){
  dat <- CellTrails:::.exDat()
  se <- embedSamples(dat)
  d <- findSpectrum(se$eigenvalues, frac=30)
  latentSpace(dat) <- se$components[, d]
  states(dat) <- findStates(dat,
                            min_size=0.01,
                            min_feat=2,
                            max_pval=1e-4,
                            min_fc=2)
  dat <- connectStates(dat, l=30)

  #TEST
  dat <- fitTrajectory(dat)
  res <- trajResiduals(dat)
  RUnit::checkEqualsNumeric(length(res), ncol(dat))
  RUnit::checkTrue(!any(res < 0))
}
