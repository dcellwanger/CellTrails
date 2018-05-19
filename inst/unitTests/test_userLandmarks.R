test_userLandmarks <- function(){
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
  ul <- c("A", "B")
  RUnit::checkException(userLandmarks(dat) <- ul) #wrong name
  ul <- colnames(dat)[seq_len(2) + 60] #one is U, the other is H
  userLandmarks(dat) <- ul
  RUnit::checkTrue(userLandmarks(dat) == ul[1]) #only one
  ul <- colnames(dat)[seq_len(2) + 62]
  userLandmarks(dat) <- ul
  RUnit::checkEquals(userLandmarks(dat), ul, checkNames=FALSE) #now 2 x U
  RUnit::checkEquals(names(userLandmarks(dat)), c("U1", "U2"))
}
