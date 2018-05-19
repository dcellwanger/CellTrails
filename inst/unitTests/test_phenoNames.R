test_phenoNames <- function(){
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
  dat <- addTrail(dat, from = "H1", to = "H2", name = "Tr1")
  dat <- addTrail(dat, from = "H2", to = "H4", name = "Tr2")

  # PHENO NAMES
  n <- c(sub("CellTrails\\.", "", colnames(colData(dat))), "landmark")
  RUnit::checkEquals(intersect(phenoNames(dat), n),
                     union(phenoNames(dat), n))
}
