test_workflow_logic <- function(){
  dat <- CellTrails:::.exDat()
  # FIND STATES
  RUnit::checkException(findStates(dat))
  # CONNECT STATES
  RUnit::checkException(connectStates(dat))
  # FIT TRAJECTORY
  RUnit::checkException(fitTrajectory(dat))
  # ADD USER LANDMARK
  RUnit::checkException(userLandmarks(dat) <- colnames(dat)[2:5])
  # ADD TRAIL
  RUnit::checkException(addTrail(dat, from="H1", to="H2", name="Tr1"))
  # REMOVE TRAIL
  RUnit::checkException(removeTrail(dat, "Tr1"))
  # FIT DYNAMIC
  RUnit::checkException(fitDynamic(dat, feature_name="feature_1",
                                   trail_name="Tr1"))
  # COMPARE DYNAMICS
  RUnit::checkException(contrastTrailExpr(dat, feature_name="feature_1",
                                          trail_names=c("Tr1", "Tr2")))
}
