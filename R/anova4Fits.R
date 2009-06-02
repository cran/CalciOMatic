`anova4Fits` <-
function(Fit_1, Fit_2) {
  ## Perform an ANalyse Of VAriance to determine the best fit among 2
  ## Useful to check whether a monoexponential or biexponential fit,
  ## preformed by nls, best predicts experimental data
  
  if(inherits(Fit_1,"nls") & inherits(Fit_2,"nls")) {
    result <- anova(Fit_1,Fit_2)
    print(result)
    SumSquares <- result[["Res.Sum Sq"]]
    ind <- which(SumSquares==min(SumSquares))
    Fit <- list(Fit_1,Fit_2)
    cat(sprintf("\nExperimental data were best explained with the %s.\n",attr(Fit[[ind]],"Name")))
  } else {
    cat(sprintf("\nFit_1 or Fit_2 do not correspond to fits performed with nls.\n"))
    ind <- NULL
  }
  return(ind)
}

