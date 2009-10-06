mkFluo4DirectFit <-
function(Ca,
                             phi,
                             S_B_340,
                             S_B_380,
                             nb_B,
                             R_min,
                             R_max,
                             K_eff,
                             K_d,
                             B_T,
                             T_340,
                             T_380,
                             P,
                             P_B,
                             SQRT=TRUE
                             ) {
  ## Function mkFluo4DirectFit
  ##
  ## Define predicted fluorescence signals that are useful in the functions:
  ##  - mkCa_MonoBiExp_4_DirectFit
  ##  - mkCa_MonoBiExp_4_DirectFit_allvar
  ##
  ## The arguments are:
  ##  - Ca: the calcium transient
  ##  - phi, S_B_340 and S_B_380: the experiment-specific parameters (see RatioSimulExp for details)
  ##  - nbBackM, R_min, R_max, K_eff, K_d, B_T, T_340, T_380 and P: see RatioSimulExp for details
  ##  - SQRT: a logical value (TRUE or FALSE, or 1 or 0) indicating whether
  ##          the square root of the fluorescence signal is considered (or the signal itself).
  ##          Setting SQRT to TRUE results in a stabilized variance, independent on the signal value.
  
  if(is.list(R_min)) R_min <- R_min$value
  if(is.list(R_max)) R_max <- R_max$value
  if(is.list(K_eff)) K_eff <- K_eff$value
  if(is.list(K_d)) K_d <- K_d$value
  
  ## Define Background and Signal fluorescences at both wavelengths
  B_340 <- rep(fluo(Ca=rep(0,nb_B), R_min=R_min, R_max=R_max, K_eff=K_eff, K_d=K_d,
                    B_T=0, phi=phi, S_B=S_B_340, T_stim=T_340, P=P, P_B=P_B))
  
  F_340 <- rep(fluo(Ca=Ca, R_min=R_min, R_max=R_max, K_eff=K_eff, K_d=K_d,
                    B_T=B_T, phi=phi, S_B=S_B_340, T_stim=T_340, P=P, P_B=P_B))
  
  B_380 <- rep(fluo(Ca=rep(0,nb_B), R_min=1, R_max=1, K_eff=K_eff, K_d=K_d,
                    B_T=0, phi=phi, S_B=S_B_380, T_stim=T_380, P=P, P_B=P_B))
  
  F_380 <- rep(fluo(Ca=Ca, R_min=1, R_max=1, K_eff=K_eff, K_d=K_d,
                    B_T=B_T, phi=phi, S_B=S_B_380, T_stim=T_380, P=P, P_B=P_B))
  
  ## Concatenate all fluorescence signals
  result <- c(B_340,F_340,B_380,F_380)^(1/(1+SQRT))

  return(result)
}

