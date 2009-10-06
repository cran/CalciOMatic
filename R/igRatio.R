igRatio <-
function(Ca,
                    t,
                    tOn=1,
                    type="mono"
                    ) {
  ## Function igRatio
  ## Returns a named list of class "initial_guess", containing initial guesses (IG)
  ## for the three/five scalar components of the mono- or bi- exponential calcium decay:
  ##     * log_Ca0: IG for the log of baseline [Ca]
  ##     * log_dCa: IG for the log of the jump amplitude of [Ca]
  ##     * log_tau: IG for the log of the time constant of the monoexp. decay
  ##            OR: IG for the log of the fast time constant of the biexp. decay
  ##  - only for a biexponential decay :
  ##     * mu: IG for the value so that fact_dCa = exp(mu)/(1+exp(mu)),
  ##           'fact' being the weight of the fast time constant of the biexp. decay
  ##     * log_dtau: IG for the log of the delta tau value (tau_s = tau_f + dtau,
  ##                 tau_s being the slow time constant of the biexp. decay)
  ##
  ## Ca: a vector of [Ca] values
  ## t: a vector of times at which the [Ca] was measured
  ## tOn: the stimulus onset time (in s)
  ## type : a character string, named "mono" or "bi", indicating which type of decay is considered
  ##
  ## The IG for Ca0 is obtained by averaging Ca prior to tOn
  ## The IG for dCa is obtained by subtracting the IG for Ca0 form the max of Ca
  ## The IG for tau is obtained from a linear regression on a rescaled and time offset version of Ca
  ## The logarithm of each IG is then stored into the "initial_guess" list
  ##
  ## ... Explain how the IG are obtained for the biexponential decay ...
  
  t_backup <- t
  Ca_backup <- Ca

  ## Estimates the baseline of the calcium transient :                  == Ca0 ==
  ## Add the case where the time vectors contains only values higher than tOn ! 
  Ca0 <- mean(Ca[t<tOn])
  if(Ca0 < 0) {
    stop(sprintf("The baseline calcium concentration is negative : %f !",Ca0))
  }
  
  ## Estimates the jump of the calcium transient :                      == dCa ==
  Ca <- Ca[t>=tOn]
  dCa <- max(Ca)-Ca0
  
  ## Normalizes the calcium transient (decay phase only)
  CaNorm <- (Ca-Ca0)/dCa
  t <- t[t>=tOn]-tOn
  t <- t-t[1]
  
  ## Takes the logarithm of the calcium decay
  N <- length(t)
  logCaNorm <- log(CaNorm)
  t <- t[is.finite(logCaNorm)]
  
  if(length(t) > 0) {
    
    if(type=="mono") {
      ## -------------------------------
      ## Case of a MONOexponential decay
      ## -------------------------------
      
      ## Adjust 1 or 2 linear fits to the logarithm of the normalized
      ## exponential, in order to estimate its slope and then :         == tau ==
      logCaNorm <- logCaNorm[is.finite(logCaNorm)]
      
      ind_limit <- length(t)
      logFit <- lm(logCaNorm[1:ind_limit] ~ t[1:ind_limit])
      tau <- as.numeric(-1/coef(logFit)[2])
      
      if(t[ind_limit] > tau) {
        ind_limit <- which(t>tau)[1]
        logFit <- lm(logCaNorm[1:ind_limit] ~ t[1:ind_limit])
        tau <- as.numeric(-1/coef(logFit)[2])
      }
      
      ## Modify the value of dCa if the first latency in t is not tOn
      ## (valid also if the first latency in t is tOn)
      ## WARNING : depends on the way the functions of the fit are coded
      dCa <- dCa * exp((t[1])/tau)
      
      result <- list(log_Ca0=log(Ca0),
                     log_dCa=log(dCa),
                     log_tau=log(tau))
      
    } else {
      ## -----------------------------
      ## Case of a BIexponential decay
      ## -----------------------------
      
      CaNorm <- CaNorm[is.finite(logCaNorm)]
      logCaNorm <- logCaNorm[is.finite(logCaNorm)]
      
      ## First, fit the biexponential decay with a monoexponential model,
      ## which gives a very rough estimate of the slow time constant :              == tau+dtau ==    
      ind_limit <- length(t)
      logFit <- lm(logCaNorm[1:ind_limit] ~ t[1:ind_limit])
      tau_Start <- as.numeric(-1/coef(logFit)[2])

      tau_slow <- tau_Start
      
      ## Second, fit susbets (of decreasing size) of the biexponential decay with
      ## a monoexponential model. When the fitted time constant reaches a minimum,
      ## this value is considered as the fast time constant of the decay :          == tau ==
      if(t[ind_limit] > min(tau_Start)) {
        ind_limit <- which(t > min(tau_Start))[1]
        logFit <- lm(logCaNorm[1:ind_limit] ~ t[1:ind_limit])
        tau_Start <- c(tau_Start,as.numeric(-1/coef(logFit)[2]))
      }
      
      while(t[ind_limit] > min(tau_Start) & zapsmall(diff(tau_Start))[length(diff(tau_Start))] != 0) {
        ind_limit <- which(t > min(tau_Start))[1]
        logFit <- lm(logCaNorm[1:ind_limit] ~ t[1:ind_limit])
        tau_Start <- c(tau_Start,as.numeric(-1/coef(logFit)[2]))
      }

      tau <- min(tau_Start)
      dtau <- tau_slow-tau

      if(dtau <= 0) dtau <- tau_slow/10
      
      ## Third, the weight of the fact time constant is then fitted with a
      ## linear fit, in order to minimize the residual sum-of-squares :             == fact ==
      BiFormula <- exp(-t/tau)-exp(-t/(tau+dtau))
      BiData <- CaNorm-exp(-t/(tau+dtau))
      BiFit <- lm(BiData ~ BiFormula-1)
      fact <- as.numeric(coef(BiFit)[1])
      if(is.na(fact)) fact <- 0.5
      if(fact<0 | fact > 1) fact <- 0.5
      
      ## In order to get better initial guesses for tau, fact and dtau,
      ## a non-linear fit is performed with optim
      useoptim <- TRUE
      if(useoptim==TRUE) {
        func_fit_optim <- function(x) {
          tau <- x[1]
          fact <- x[2]
          dtau <- x[3]
          y <- caBiExp(t_backup,tOn,Ca0,dCa,tau,fact,dtau)
          Res <- sum((Ca_backup-y)^2)
          return(Res)
        }
        
        Optim_CaNorm <- optim(c(tau,fact,dtau),func_fit_optim,method="Nelder-Mead",
                              control=list(abstol=1e-12,reltol=1e-8))
        
        tau <- Optim_CaNorm$par[1]
        dtau <- Optim_CaNorm$par[3]
        fact <- Optim_CaNorm$par[2]
        if(fact<0 | fact > 1) {
          fact <- 0.5
          func_fit_optim2 <- function(x) {
            tau <- x[1]
            dtau <- x[2]
            y <- caBiExp(t_backup,tOn,Ca0,dCa,tau,fact,dtau)
            Res <- sum((Ca_backup-y)^2)
            return(Res)
          }
          Optim_CaNorm2 <- optim(c(tau,dtau),func_fit_optim2,method="Nelder-Mead",
                                 control=list(abstol=1e-12,reltol=1e-8))
          tau <- Optim_CaNorm2$par[1]
          dtau <- Optim_CaNorm2$par[2]
        }
      }
            
      ## Security test : is tau_f or dtau < 0 ?
      if(dtau < 0 & tau < 0) {
        dtau <- 1
        tau <- 1
      } else {
        if(tau < 0) tau <- dtau
        if(dtau < 0) dtau <- tau/10
      }
      
      ## Modify the value of dCa if the first latency in t is not tOn
      ## (valid also if the first latency in t is tOn)
      ## WARNING : depends on the way the functions of the fit are coded
      dCa <- dCa * (fact * exp(t[1]/tau) + (1-fact) * exp(t[1]/(tau+dtau)))

      result <- list(log_Ca0=log(Ca0),
                     log_dCa=log(dCa),
                     log_tau=log(tau),
                     mu=log(fact/(1-fact)),
                     log_dtau=log(dtau)
                     )
    }
  } else { ## if length(t) == 0, ie if there is some kind of problem in the signal...
    if(type=="mono") {
      result <- list(log_Ca0=log(Ca0), log_dCa=0, log_tau=0)
    } else {
      result <- list(log_Ca0=log(Ca0), log_dCa=0, log_tau=0, mu=0, log_dtau=0)
    }
  }
  
  class(result) <- "initial_guess"
  
  return(result)
}

