caBiExp <-
function(t=1,
                    tOn=1,
                    Ca0=0.05,
                    dCa=0.1,
                    tau=3,
                    fact=1,
                    dtau=2
                    ) {
  ## Function caBiExp
  ## Returns a vector of [Ca] vs time values where a [Ca] jump occurs at tOn,
  ## followed by a biexponential return to baseline value
  ##
  ## t: a vector of time values
  ## tOn: time of the concentration jump (in s)
  ## Ca0: baseline [Ca] in muM
  ## dCa: change in [Ca] (in muM) occurring at tOn
  ## tau: fast time constant (in s) of the biexponential return to baseline
  ## fact: factor (between 0 and 1) applied to dCa to determine the weigth of the fast exponential decay
  ## dtau: delta time constant (in s) of the biexponential return to baseline
  
  ## Create the biexponential calcium transient
  result <- Ca0 + dCa * ifelse(t>=tOn, fact*exp(-(t-tOn)/tau) + (1-fact)*exp(-(t-tOn)/(tau+dtau)), 0)

  attr(result,"Time") <- t
  attr(result,"tOn") <- tOn
  
  return(result)
}

