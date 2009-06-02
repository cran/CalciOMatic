`igDirect` <-
function(adu_B_340,
                     adu_340,
                     adu_B_380,
                     adu_380,
                     ig_ratio,
                     t,
                     tOn=1,
                     subset=1:length(t),
                     R_min=0.136,
                     R_max=2.701,
                     K_eff=3.637,
                     K_d=0.583,
                     B_T=100,
                     T_340=0.015,
                     T_380=0.006,
                     P=400,
                     P_B=400
                     ) {
  ## Function igDirect
  ## Returns a named list of class "initial_guess" containing initial guesses (IG)
  ## for the remaining parameters of the fluorescence signals F_340 and F_380.
  ## The components of the named list are:
  ##    * log_phi: log of phi
  ##    * log_S_B_340: log of the background fluorescence intensity at 340 nm
  ##    * log_S_B_380: log of the background fluorescence intensity at 380 nm
  ## 
  ## adu_B_340: vector of background measurements at 340 nm
  ## adu_340: vector of "transient" measurements at 340 nm
  ## adu_B_380: vector of background measurements at 380 nm
  ## adu_380: vector of transient measurements at 380 nm
  ## ig_ratio: named list returned by igRatio
  ## t: a vector with the time at which the [Ca] was measured
  ## tOn: the stimulus onset time
  ## subset: a vector of indices to use for the initial guess estimation
  ## R_min: R_min parameter from calibration experiments (dimensionless),
  ##        representing the minimum ratiometric measurement observable
  ##        [see Eq. 21, p10]
  ## R_max: R_max parameter from calibration experiments (dimensionless),
  ##        representing the maximum ratiometric measurement observable
  ##        [see Eq. 22, p10]
  ## K_eff: Effective fura dissociation constant from calibration experiments (in muM),
  ##      [see Eq. 25, p10]
  ## K_d: fura dissociation constant from calibration experiments (in muM),
  ##      [see Sec. 3.2]
  ## B_T: the total fura concentration (in muM)
  ## T_340: exposure time at 340 nm (in s)
  ## T_380: exposure time at 380 nm (in s)
  ## P: number of pixels used for the data binning
  
  ## Estimate the background fluorescence at each wavelength
  S_B_340 <- mean(adu_B_340)/T_340/P_B ##                                 == S_B_340 ==
  S_B_380 <- mean(adu_B_380)/T_380/P_B ##                                 == S_B_380 ==
  
  ## Work differently following the mode of t, adu_340 and adu_380
  if(is.vector(t)) {
    ## Consider a subset of the fluorescence counts at each wavelength
    t <- t[subset]
    adu_340 <- adu_340[subset]
    adu_380 <- adu_380[subset]

    ## Define the data and formula for the linear fit
    Y <- c(adu_340/T_340/P - S_B_340, adu_380/T_380/P - S_B_380)

    Ca_pred <- caMonoBiExpFromIG(t=t, tOn=tOn, ig=ig_ratio)
    X <- c(B_T / (K_d+Ca_pred) * c(R_min*K_eff+R_max*Ca_pred,K_eff+Ca_pred))

  } else if(is.matrix(t)) {
    ## Consider a subset of the fluorescence counts at each wavelength
    t_2 <- matrix(NA,nrow=dim(t)[1],ncol=length(subset))
    t_2[,] <- t[,subset]
    adu_340_2 <- matrix(NA,nrow=dim(t)[1],ncol=length(subset))
    adu_340_2[,] <- adu_340[,subset]
    adu_380_2 <- matrix(NA,nrow=dim(t)[1],ncol=length(subset))
    adu_380_2[,] <- adu_380[,subset]
    B_T_2 <- matrix(NA,nrow=dim(t)[1],ncol=length(subset))
    B_T_2[,] <- B_T[,subset]
    
    ntrans <- dim(t)[1]
    
    ## Define the data and formula for the linear fit
    Y <- c()
    X <- c()
    
    for(k in 1:ntrans) {
      nig <- length(ig_ratio)/ntrans
      Y <- c(Y, adu_340_2[k,]/T_340/P - S_B_340, adu_380_2[k,]/T_380/P - S_B_380)
      Ca_pred <- caMonoBiExpFromIG(t=t_2[k,],
                                   tOn=tOn[k],
                                   ig=ig_ratio[((k-1)*nig+1):(k*nig)])
      X <- c(X, B_T_2[k,] / (K_d+Ca_pred) * c(R_min*K_eff+R_max*Ca_pred,K_eff+Ca_pred))
    }
  }
  
  ## Perform a linear fit with a 0 intercept (with -1, or 0+) to determine phi  
  myfit <- lm(Y ~ X - 1)
  phi <- as.numeric(coef(myfit)[1])
    
  result <- c(ig_ratio,
              log_phi=log(phi),
              log_S_B_340=log(S_B_340),
              log_S_B_380=log(S_B_380)
              )
  
  class(result) <- "initial_guess"
  
  return(result)
}

