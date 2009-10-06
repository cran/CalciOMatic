caMonoExp <-
function(t=1,
                      tOn=1,
                      Ca0=0.05,
                      dCa=0.1,
                      tau=3
                      ) {
  ## Function caMonoExp
  ## Returns a vector of [Ca] vs time values where a [Ca] jump occurs at tOn,
  ## followed by a monoexponential return to baseline
  ##
  ## t: a vector of time values
  ## tOn: time of the concentration jump (in s)
  ## Ca0: baseline [Ca] (in muM)
  ## dCa: change in [Ca] (in muM) occurring at tOn
  ## tau: time constant (in s) of the monoexponential return to baseline

  ## Create the monoexponential calcium transient
  result <- Ca0 + dCa * ifelse(t>=tOn, exp(-(t-tOn)/tau), 0)
  
  attr(result,"Time") <- t
  attr(result,"tOn") <- tOn  
  
  return(result)
}

