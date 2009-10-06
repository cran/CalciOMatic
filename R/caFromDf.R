caFromDf <-
function(df,
                     numTransient=1,
                     Plot=FALSE) {
  ## Function caFromDf
  ## Returns a vector of estimated [Ca] values using the "classical"
  ## ratiometric method, from a data frame containing the 4 adu vectors,
  ## the stimulation durations and the calibration parameters
  ## 
  ## df: a simulated data frame returned by RatioSimulExp,
  ##     or an experimental data frame returned by RatioPhysioExp
  ##     [for details about the structure of the data frame, see Ratiometric_Experiments.R]
  ## transient: the number of the transient to use
  ## Plot: a logical value indicating to plot (or not) the ratiometric-obtained calcium transient
  
  ## Background measurements
  lambda_vec <- unique(with(df,lambda))
  
  adu_B_340 <- with(df,adu[is.na(Time) & lambda==lambda_vec[1]])
  adu_B_380 <- with(df,adu[is.na(Time) & lambda==lambda_vec[2]])
  
  adu_340 <- with(df,adu[is.finite(Time) & lambda==lambda_vec[1] & transient==numTransient])
  adu_380 <- with(df,adu[is.finite(Time) & lambda==lambda_vec[2] & transient==numTransient])
  
  t <- with(df,Time[is.finite(Time) & lambda==lambda_vec[1] & transient==numTransient])
  tOn <- attr(df,"tOn")[numTransient]
  
  T_stim <- attr(df,"T_stim")
  T_340 <- T_stim[1]
  T_380 <- T_stim[2]
  P <- attr(df,"P")
  P_B <- attr(df,"P_B")
  R_min <- attr(df,"R_min")
  R_max <- attr(df,"R_max")
  K_eff <- attr(df,"K_eff")
  
  Ca <- caFromRatio(adu_B_340=adu_B_340,
                    adu_340=adu_340,
                    adu_B_380=adu_B_380,
                    adu_380=adu_380,
                    T_340=T_340,
                    T_380=T_380,
                    P=P,
                    P_B=P_B,
                    R_min=R_min,
                    R_max=R_max,
                    K_eff=K_eff,
                    Plot=Plot
                    )
  
  attr(Ca,"Time") <- t
  attr(Ca,"tOn") <- tOn
  
  if(Plot==TRUE) plot(t,Ca,type="l")
  
  return(Ca)
}

