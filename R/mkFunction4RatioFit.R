`mkFunction4RatioFit` <-
function(type="mono") {
  ## Function mkFunction4RatioFit
  ##
  ## Returns a function prediciting the exponential time course of [Ca].
  ##
  ## If type == "mono", the function depends on the following 5 variables:
  ##  - t, tOn, log_Ca0, log_dCa and log_tau (see the Ca_MonoExp_fct function)
  ##
  ## If type ~= "mono", the function depends on the following 7 variables:
  ##  - t, tOn, log_Ca0, log_dCa, mu, log_tau_f and log_dtau (see the Ca_BiExp_fct function)
  
  if(type=="mono") {
    result <- function(t,
                       tOn,
                       log_Ca0,
                       log_dCa,
                       log_tau
                       ) {
      caMonoExp(t=t, tOn=tOn,
                Ca0=exp(log_Ca0), dCa=exp(log_dCa), tau=exp(log_tau))
    }
  } else {
    result <- function(t,
                       tOn,
                       log_Ca0,
                       log_dCa,
                       log_tau,
                       mu,
                       log_dtau
                       ) {
      caBiExp(t=t, tOn=tOn,
              Ca0=exp(log_Ca0), dCa=exp(log_dCa), tau=exp(log_tau),
              fact=exp(mu)/(1+exp(mu)), dtau=exp(log_dtau))
    }
  }

  return(result)
}

