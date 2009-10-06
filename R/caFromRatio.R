caFromRatio <-
function(adu_B_340,
                        adu_340,
                        adu_B_380,
                        adu_380,
                        T_340=0.015,
                        T_380=0.006,
                        P,
                        P_B,
                        R_min=0.136,
                        R_max=2.701,
                        K_eff=3.637,
                        Plot=FALSE
                        ) {
  ## Function caFromRatio
  ## Returns a vector of estimated [Ca] values using the "classical"
  ## ratiometric method, from the 4 vectors of adu
  ## The calcium transient comes with a variance/covariance matrix as attribute
  ##
  ## adu_B_340: vector of background adu counts at 340 nm
  ## adu_340: vector of adu counts at 340 nm
  ## adu_B_380: vector of background adu counts at 380 nm
  ## adu_380: vector of adu counts at 380 nm
  ## T_340: exposure time at 340 nm (in s)
  ## T_380: exposure time at 380 nm (in s)
  ## R_min: R_min parameter from calibration experiments (dimensionless), Eq. 21, p10
  ## R_max: R_max parameter from calibration experiments (dimensionless), Eq. 22, p10
  ## K_eff: Effective fura dissociation constant in muM (from calibration experiments),
  ##      Eq. 25, p10
  
  R_min_0 <- R_min
  R_max_0 <- R_max
  K_eff_0 <- K_eff
  
  if(!is.null(R_min_0$mean)) R_min_0 <- R_min_0$mean
  if(!is.null(R_max_0$mean)) R_max_0 <- R_max_0$mean
  if(!is.null(K_eff_0$mean)) K_eff_0 <- K_eff_0$mean
  
  ## Ratiometric transformation
  if(length(adu_B_340) == length(adu_340)) {
    adu_340_B <- adu_340 - adu_B_340*P/P_B
    adu_380_B <- adu_380 - adu_B_380*P/P_B
  } else {
    adu_340_B <- adu_340 - mean(adu_B_340)*P/P_B
    adu_380_B <- adu_380 - mean(adu_B_380)*P/P_B
  }
  
  ## Add 0.5 for continuity correction
  ## adu_340_B <- adu_340_B + 0.5
  ## adu_380_B <- adu_380_B + 0.5
  R <- adu_340_B * T_380 / (adu_380_B * T_340)
  Ca <- K_eff_0 * (R-R_min_0) / (R_max_0-R)
  
  ## Partial derivatives of Ca
  dCa_dadu340 <- K_eff_0 / adu_340_B * R * (R_max_0-R_min_0) / (R_max_0-R)^2
  dCa_dadu340B <- -dCa_dadu340
  
  dCa_dadu380 <- -K_eff_0 / adu_380_B * R * (R_max_0-R_min_0) / (R_max_0-R)^2
  dCa_dadu380B <- -dCa_dadu380
  
  dCa_dRmin <- -Ca / (R-R_min_0)
  dCa_dRmax <- -Ca / (R_max_0-R)
  dCa_dKeff <- Ca / K_eff_0
  
  ## Build the variance-covariance matrix
  
  ## Start with the Ca2+ variance (diagonal terms)
  var_Ca_diag <- dCa_dadu340^2 * adu_340 + dCa_dadu380^2 * adu_380
  
  if(length(adu_B_340) == length(adu_340)) {
    var_Ca_diag <- var_Ca_diag +
                   dCa_dadu340B^2 * adu_B_340 +
                   dCa_dadu380B^2 * adu_B_380
  }
  
  ## Continue with the covariance terms (diagonal + non-diagonal terms)
  cov_Ca_non_diag <- 0
  
  if(length(adu_B_340) != length(adu_340)) {
    cov_Ca_non_diag <- cov_Ca_non_diag +
                       dCa_dadu340B %o% dCa_dadu340B * mean(adu_B_340) +
                       dCa_dadu380B %o% dCa_dadu380B * mean(adu_B_380)
  }
  
  if(!is.null(R_min$USE_se)) {
    if(R_min$USE_se == TRUE) {
      ## cov_Ca_non_diag <- cov_Ca_non_diag + dCa_dRmin %o% dCa_dRmin * R_min$se^2
    }
  }
  
  if(!is.null(R_max$USE_se)) {
    if(R_max$USE_se == TRUE) {
      ## cov_Ca_non_diag <- cov_Ca_non_diag + dCa_dRmax %o% dCa_dRmax * R_max$se^2
    }
  }
  
  if(!is.null(K_eff$USE_se)) {
    if(K_eff$USE_se == TRUE) {
      ## cov_Ca_non_diag <- cov_Ca_non_diag + dCa_dKeff %o% dCa_dKeff * K_eff$se^2
    }
  }
  
  ## Build the entire matrix
  cov_Ca <- diag(var_Ca_diag) + cov_Ca_non_diag
  
  attr(Ca,"cov") <- cov_Ca
  
  if(Plot==TRUE) {
    layout(matrix(c(1,2,3,3),nrow=4,ncol=1))
    plot(adu_340,type="l")
    plot(adu_380,type="l")
    plot(Ca,type="l")
  }
  
  return(Ca)
}

