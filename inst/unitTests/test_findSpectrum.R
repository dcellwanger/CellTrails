test_findSpectrum <- function(){
  dat <- CellTrails:::.exDat()
  se <- embedSamples(dat)

  #TEST
  d <- findSpectrum(se$eigenvalues, frac=30)
  RUnit::checkEquals(d, seq_len(10))
  latentSpace(dat) <- se$components[, d]
  RUnit::checkException(latentSpace(dat) <- se$components[seq_along(2), d])
}
