`mkFunction4DirectFit` <-
function(type="mono",
                                 nb_B=5,
                                 transients=1,
                                 alphamethod=TRUE,
                                 SQRT=TRUE
                                 ) {
  ## Function mkFunction4DirectFit
  ## 
  ## Returns a function prediciting the time courses of F_340, F_380 and the
  ## corresponding background measurements, for an exponential time course of [Ca].
  ##
  ## If type == "mono", the function depends on the following 13 variables:
  ##  - log_Ca0, log_dCa, log_tau, log_phi, log_S_B_340, log_S_B_380,
  ##    R_min, R_max, K_eff, K_d, log_B_T, dB_T and log_tau_B_T
  ##    (see the Ca_MonoExp_fct, B_T_MonoExp_fct and mkFluo_4_DirectFit functions)
  ##
  ## If type == "bi", the function depends on the following 15 variables:
  ##  - log_Ca0, log_dCa, log_tau, mu, log_dtau, log_phi, log_S_B_340, log_S_B_380,
  ##    R_min, R_max, K_eff, K_d, log_B_T, dB_T and log_tau_B_T
  ##    (see the Ca_BiExp_fct, B_T_MonoExp_fct and mkFluo_4_DirectFit functions)
  ##
  ## The arguments are:
  ## nb_B, T_340, T_380 and P: see RatioSimulExp for details
  ##
  ## SQRT: a logical value (TRUE or FALSE, or 1 or 0) indicating whether
  ##       the square root of the fluorescence signal is considered (or the signal itself).
  ##       Setting SQRT to TRUE results in a stabilized variance, independent on the signal value.
  
  ## ------------------------
  ## Create an empty function
  result <- function() NULL

  
  ## ------------------------
  ## Define its arguments ...
  S4F <- "formals(result) <- list(t=NULL, tOn=NULL,"
  if(alphamethod==TRUE) {
    S4F <- c(S4F,"adu_B_340=NULL, adu_340=NULL,")
    S4F <- c(S4F,"adu_B_380=NULL, adu_380=NULL,")
  }
  S4F <- c(S4F,"T_340=NULL, T_380=NULL, P=NULL, P_B=NULL,")
  for(k in transients) {
    if(type=="mono") {
      S4F <- c(S4F, sprintf("log_Ca0_%d=NULL, log_dCa_%d=NULL, log_tau_%d=NULL,",k,k,k))
    } else if(type=="bi") {
      S4F <- c(S4F, sprintf("log_Ca0_%d=NULL, log_dCa_%d=NULL, log_tau_%d=NULL, mu_%d=NULL, log_dtau_%d=NULL,",k,k,k,k,k))
    }
  }
  S4F <- c(S4F,"log_phi=NULL, log_S_B_340=NULL, log_S_B_380=NULL,")
  S4F <- c(S4F,"log_R_min=NULL, log_R_max=NULL, log_K_eff=NULL, log_K_d=NULL,")
  if(alphamethod==TRUE) {
    S4F <- c(S4F,"alpha=NULL);")
  } else {
    S4F <- c(S4F,"B_T=NULL);")
  }  
  eval(parse(text=S4F))
  

  ## ----------------
  ## ... and its body
  S4B <- "body(result) <- expression({"
  
  ## Get estimated Ca transient
  S4B <- c(S4B, "Ca <- c();")
  for(k in transients) {
    if(type=="mono") {
      S4B <- c(S4B, sprintf("Ca <- c(Ca, caMonoExp(t=t['%d',],tOn=tOn['%d'],Ca0=exp(log_Ca0_%d),dCa=exp(log_dCa_%d),tau=exp(log_tau_%d)));",k,k,k,k,k))
    } else if(type=="bi") {
      S4B <- c(S4B, sprintf("Ca <- c(Ca, caBiExp(t=t['%d',],tOn=tOn['%d'],Ca0=exp(log_Ca0_%d),dCa=exp(log_dCa_%d),tau=exp(log_tau_%d),fact=exp(mu_%d)/(1+exp(mu_%d)),dtau=exp(log_dtau_%d)));",k,k,k,k,k,k,k,k))
    }
  }
  
  ## Get estimated B_T transient
  if(alphamethod==TRUE) {
    S4B <- c(S4B, "B_T <- c();")
    for(k in transients) {
      S4B <- c(S4B, sprintf("B_T <- c(B_T, adu_340['%d',]/(T_340*P) - adu_B_340/(T_340*P_B) + alpha * (adu_380['%d',]/(T_380*P) - adu_B_380/(T_380*P_B)));",k,k))
    }
  } else {
    S4B <- c(S4B, "B_T <- B_T;")
  }
    
  ## Get estimated fluorescence signals
  S4B <- c(S4B, "result2 <- mkFluo4DirectFit(Ca=Ca, phi=exp(log_phi), S_B_340=exp(log_S_B_340), S_B_380=exp(log_S_B_380), nb_B=nb_B, R_min=exp(log_R_min), R_max=exp(log_R_max), K_eff=exp(log_K_eff), K_d=exp(log_K_d), B_T=B_T, T_340=T_340, T_380=T_380, P=P, P_B=P_B, SQRT=SQRT);")
  
  ## Return result2
  S4B <- c(S4B, "return(result2)})")
  
  ## print(S4B)
  
  eval(parse(text=S4B))
  
  return(result)
}

