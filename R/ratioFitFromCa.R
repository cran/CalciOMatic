ratioFitFromCa <- function(Ca,
                           t,
                           tOn,
                           type="mono",
                           ig=NULL,
                           Plot=FALSE,
                           Fit=TRUE,
                           AfterPeak=FALSE
                           ) {
  ## Function ratioFitFromCa
  ## 
  ## Fit the fluorescence signals present in the data frame 'df',
  ## with the classical ratiometric method. The resulting calcium
  ## signal is fitted with a mono- or bi- exponential decay,
  ## depending on the value of 'type'.
  ## 
  ## The fit uses the covariance matrix of the calcium transients
  ## 
  ## Ca: a vector of [Ca] vs. time
  ## t: a vector time values at which the calcium was recorded
  ## tOn: the time of the calcium jump (in s)
  ## type: a character string indicating the type of calcium decay to fit
  ##       ("mono"exponential or "bi"exponential)
  ## Plot: a logical value; set to TRUE to plot the raw data,
  ##                        the initial guess and the fitted signals
  ## Fit: a logical value; set to TRUE to perform the fit;
  ##                       set to FALSE to estimate an initial guess only
  ## AfterPeak: a logical value; set to TRUE to fit only the convex part of the decay
  ##                             set to FALSE to fit the whole signal
  ## Trace: a logical value; set to TRUE to print the residual sum-of-squares
  ##                         at each iteration of the nls fitting process
  ## WarnOnly: a logical value; set to TRUE to go on even if the fit produced an error
  ## 
  ## The function returns an "nls"/"calcium_ratio_fit" object with the following attributes:
  ##   - Name: "mono-" or "bi-"exponential ratiometric fit
  ##   - Time: a copy of the attribute "Time" of the dataframe df
  ##   - RawData: the raw calcium transient resulting from the ratiometric transformation
  ##   - RawDataFrame: copy of argument df
  ##   - FitFunction: the function used in the nls formula
  
  ## Plot "experimental" data
  if(Plot==TRUE) {X11(); plot(t,Ca,'l')}
  
  ## Define a "mysubset" attribute to the time vector
  attr(t,"mysubset") <- NULL
  
  ## Optionally: keep only the convex/concave part of the adu transients (after the peak)
  idxOn <- which(t>=tOn)[1]
  if(AfterPeak==TRUE & is.logical(AfterPeak)) {
    idxConvex <- transientConvexPart(t=t, tOn=tOn, transient=Ca)
    if(idxConvex > idxOn+1) attr(t,"mysubset") <- idxOn:(idxConvex-1)
  } else if(AfterPeak >= 1 & is.numeric(AfterPeak)) {
    attr(t,"mysubset") <- idxOn:(idxOn+AfterPeak-1)
  }
  
  ## Define the subset of time vector to use for the initial guess
  mysubset <- 1:length(t)
  try(mysubset <- mysubset[-attr(t,"mysubset")],silent=TRUE)
  
  ## Initial guess for the 3 (resp. 5) parameters of the calcium transient
  ## Only the convex part of the signal is taken into account, if relevent
  if(is.null(ig)) {
    ig <- igRatio(Ca=Ca[mysubset], t=t[mysubset], tOn=tOn, type=type)
  }
  
  result <- ig
  
  ## Plot the initial transient over the experimental data
  if(Plot==TRUE) {
    Initial_Transient <- caMonoBiExpFromIG(t=t, tOn=tOn, ig=ig)
    lines(t[-mysubset], Initial_Transient[-mysubset], col='blue', type="p", pch=20)
    lines(t[mysubset], Initial_Transient[mysubset], col='blue')
  }

  if(Fit==TRUE) {
    ## Define the fit formula, which depends on the decay type ("mono", or "bi")
    ## ratio_fit_fct <- mkFunction4RatioFit_v2(type=type)
    
    ## Define the string to evaluate to perform the fit
    ## arg_list <- formals(ratio_fit_fct)
    ## arg_names <- names(arg_list)
    ## arg_2_keep <- 1:(6+2*(type=="bi"))
    ## L <- length(arg_2_keep)
    
    ## Define the superior triangular matrix
    ## to use for the least-squares regression
    cov <- attr(Ca,"cov")
    cov <- cov[mysubset,mysubset]
    vocChol <- t(solve(chol(cov)))
    ## A <- chol(cov)
    ## B <- t(solve(A))
    
    ## Put the variables of ig in the environment of the current function
    for(i in 1:length(ig)) assign(names(ig)[i],ig[[i]])

    ratio_fit_fct <- mkFunction4RatioFit(type=type)
    
    rssFct <- function(p) {
      if(type == "mono") {
        res <- sum((vocChol %*% (Ca[mysubset] - ratio_fit_fct(t[mysubset],
                                                              tOn,
                                                              p[1],
                                                              p[2],
                                                              p[3])
                                 )
                    )^2
                   )
      } else if(type == "bi") {
        res <- sum((vocChol %*% (Ca[mysubset] - ratio_fit_fct(t[mysubset],
                                                              tOn,
                                                              p[1],
                                                              p[2],
                                                              p[3],
                                                              p[4],
                                                              p[5])
                                 )
                    )^2
                   )
      }
      return(res)
    }
    
    keepGoing <- TRUE
    startValue <- as.vector(unlist(ig))
    nbFits <- 1
    
    while(keepGoing == TRUE && nbFits <= 10) {
      calcium_ratio_fit <- optim(startValue, rssFct, method="BFGS", hessian=TRUE)
      keepGoing <- calcium_ratio_fit$convergence > 0
      startValue <- calcium_ratio_fit$par
      nbFits <- nbFits+1
    }
    
    ## Plot the predicted data over the experimental data
    if(Plot==TRUE) lines(t[mysubset], predict(calcium_ratio_fit), col='red', lwd=2)
    
    attr(calcium_ratio_fit,"Name") <- sprintf("%sexponential ratiometric fit",type)
    attr(calcium_ratio_fit,"type") <- type
    attr(calcium_ratio_fit,"Time") <- t
    attr(calcium_ratio_fit,"tOn") <- tOn
    attr(calcium_ratio_fit,"FitFunction") <- ratio_fit_fct
    attr(calcium_ratio_fit,"RawData") <- Ca
    attr(calcium_ratio_fit,"RSSFunction") <- rssFct
    attr(calcium_ratio_fit,"Subset") <- mysubset
    attr(calcium_ratio_fit,"VocChol") <- vocChol
    
    class(calcium_ratio_fit) <- c("ratio_fit")
    
    result <- calcium_ratio_fit
  }
  
  return(result)
}
