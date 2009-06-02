`caMonoBiExpFromIG` <-
function(t=1,
                              tOn=1,
                              ig=NULL
                              ) {
  ## Function caMonoBiExpFromIG
  ## Returns a vector of [Ca] vs time values where a [Ca] jump occurs at tOn
  ## followed by a mono- or bi- exponential return to baseline value,
  ## depending on the length of ig
  ##
  ## t: a vector of time values
  ## tOn: time of the concentration jump (in s)
  ## ig: an "initial_guess" object, giving the parameters of the decay
  ##
  ## The class of the returned calcium transient is "calcium_transient"
  
  if(!inherits(ig,"initial_guess") & !inherits(ig,"list")) {
    stop("The third argument of the 'caMonoBiExpFromIG' function must be of class 'initial_guess'")
  }

  if(length(ig)<3 | length(ig)>5) {
    stop("The third argument of the 'caMonoBiExpFromIG' function must contain either 3 or 5 elements")
  }

  N <- names(ig)

  if(is.finite(pmatch("log_Ca0",N)) &
     is.finite(pmatch("log_dCa",N)) &
     is.finite(pmatch("log_tau",N)) &
     is.na(pmatch("mu",N)) &
     is.na(pmatch("log_dtau",N))) {
    result <- caMonoExp(t=t,
                        tOn=tOn,
                        Ca0=exp(ig[[pmatch("log_Ca0",N)]]),
                        dCa=exp(ig[[pmatch("log_dCa",N)]]),
                        tau=exp(ig[[pmatch("log_tau",N)]])
                        )
  } else if(is.finite(pmatch("log_Ca0",N)) &
            is.finite(pmatch("log_dCa",N)) &
            is.finite(pmatch("log_tau",N)) &
            is.finite(pmatch("mu",N)) &
            is.finite(pmatch("log_dtau",N))) {
    result <- caBiExp(t=t,
                      tOn=tOn,
                      Ca0=exp(ig[[pmatch("log_Ca0",N)]]),
                      dCa=exp(ig[[pmatch("log_dCa",N)]]),
                      tau=exp(ig[[pmatch("log_tau",N)]]),
                      fact=exp(ig[[pmatch("mu",N)]])/(1+exp(ig[[pmatch("mu",N)]])),
                      dtau=exp(ig[[pmatch("log_dtau",N)]])
                      )
  } else {
    print(ig)
    stop("The third argument of the 'caMonoBiExpFromIG' function does not have the right format")
  }
  
  return(result)
}

