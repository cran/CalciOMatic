transientConvexPart <-
function(transient,
                                t=1,
                                tOn=1
                                ) {
  ## Function transientConvexPart
  ##
  ## Smoothes the fluorescence transients and determine the convex/concave part
  ## (after the peak). The function returns the index after which the transient
  ## becomes convex/concave.
  ##
  ## t: a vector of time values (in s)
  ## tOn: the stimulus onset time (in s)
  ## transient: a vector of either Calcium or Fluorescence values
  
  transient_smooth <- smooth.spline(t,transient)
  transient_smooth_dt2 <- predict(transient_smooth,deriv=2)
  
  transient_mean0 <- mean(transient[t<tOn])
  transient_min <- min(transient[t>=tOn])
  transient_max <- max(transient[t>=tOn])
  
  if(abs(transient_min-transient_mean0) <= abs(transient_max-transient_mean0)) {
    ## Case Convex
    idxMinMax <- which(t>=tOn & transient==max(transient))[1]
    tMinMax <- t[idxMinMax]
    idx_result <- which(t>=tMinMax & transient_smooth_dt2$y>=0)[1]
  } else {
    ## Case Concave
    idxMinMax <- which(t>=tOn & transient==min(transient))[1]
    tMinMax <- t[idxMinMax]
    idx_result <- which(t>=tMinMax & transient_smooth_dt2$y<=0)[1]
  }
  if(is.na(idx_result)) idx_result <- which(t>=tOn)[1]

  return(idx_result)
}

