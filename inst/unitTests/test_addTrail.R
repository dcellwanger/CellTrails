test_addTrail <- function(){
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
  dat <- fitTrajectory(dat)

  #TEST
  dat <- addTrail(dat, from = "H1", to = "H2", name = "Tr1")
  RUnit::checkEquals(nrow(trails(dat)), ncol(dat))
  RUnit::checkEquals(ncol(trails(dat)), 1)
  dat <- addTrail(dat, from = "H2", to = "H4", name = "Tr2")
  RUnit::checkEquals(ncol(trails(dat)), 2)
  RUnit::checkEquals(colnames(trails(dat)), c("Tr1", "Tr2"))
}
