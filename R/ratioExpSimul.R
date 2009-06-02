`ratioExpSimul` <-
function(nb_B=5,
                          Ca,
                          R_min=0.136,
                          R_max=2.701,
                          K_eff=3.637,
                          K_d=0.583,
                          B_T=100,
                          phi=1.25,
                          S_B_340=10,
                          S_B_380=10,
                          T_340=0.015,
                          T_380=0.006,
                          P=400,
                          P_B=400,
                          ntransients=1,
                          G=1,
                          s_ro=0
                          ) {
  ## Function ratioExpSimul
  ## Simulates one or more ratiometric experiment(s) with a given [Ca] transient
  ## (usually a mono- or bi-exponential [Ca] return to baseline)
  ##
  ## nb_B: the number of background measurements performed at each wavelength
  ## Ca: ideal (unnoised) calcium transient, previously defined (in muM)
  ##     with the ca_monoexp_fct or ca_biexp_fct functions
  ## R_min: R_min parameter from calibration experiments (dimensionless),
  ##        representing the minimum ratiometric measurement observable
  ##        [see Eq. 21, p10]
  ## R_max: R_max parameter from calibration experiments (dimensionless),
  ##        representing the maximum ratiometric measurement observable
  ##        [see Eq. 22, p10]
  ## K_eff: Effective fura dissociation constant from calibration experiments (in muM),
  ##      [see Eq. 25, p10]
  ## K_d: fura dissociation constant in muM (from calibration experiments),
  ##      [see Sec. 3.2]
  ## B_T: the total fura concentration (in muM), single number or a vector (same size as Ca)
  ## phi: dimensionless experiment specific parameter, single number or rng function
  ## S_B_340: background (+ dark current) fluorescence intensity at 340 nm
  ##          in count/pixel/sec, single number or rng function
  ## S_B_380: background (+ dark current) fluorescence intensity at 380 nm
  ##          in count/pixel/sec, single number or rng function
  ## T_340: exposure time at 340 nm (in s)
  ## T_380: exposure time at 380 nm (in s)
  ## P: number of pixels used for the fluorescence data binning
  ## P_B: number of pixels used for the background data binning
  ## ntransients: the number of transients to simulate
  ##
  ## A data frame of class "fluo_rawdata" is returned with the following variables:
  ##  - adu: the counts read out of the "camera"
  ##  - Time: the time (in s) at which the adu is read
  ##  - lambda: the excitation wavelength (in nm), a factor
  ##  - transient: the transient number
  ##
  ## The returned data frame has also the following attributes:
  ##  - nb_B: copy of arg bn_B
  ##  - tOn: the time (in s) at which the depolarisation is set
  ##  - T_stim: a vector containing the stimulation duration (in s) at both wavelengths
  ##  - R_min: copy of arg R_min
  ##  - R_max: copy of arg R_max
  ##  - K_eff: copy of arg K_eff
  ##  - K_d: copy of arg K_d
  ##  - P: copy of arg P
  ##  - P_B: copy of arg P_B
  ##  - B_T: the total fura concentration in the pipette (in muM)
  ##  - (opt) Ca_transient: the calcium transient used in the "experiments"
  ##  - (opt) B_340: S_B_340*T_340*P
  ##  - (opt) B_380: S_B_380*T_380*P
  ##  - (opt) F_340: the F_340 time course of the "experiments"
  ##  - (opt) F_380: the F_380 time course of the "experiments"
  ##  - (opt) S_B_340: copy of arg S_B_340
  ##  - (opt) S_B_380: copy of arg S_B_380
  ##  - (opt) phi: phi value used in each "experiment"
  
  ## Ideal fluorescence at each wavelength
  adu <- mkFluo4DirectFit(Ca=rep(Ca,ntransients),
                          phi=phi,
                          S_B_340=S_B_340,
                          S_B_380=S_B_380,
                          nb_B=nb_B,
                          R_min=R_min,
                          R_max=R_max,
                          K_eff=K_eff,
                          K_d=K_d,
                          B_T=B_T,
                          T_340=T_340,
                          T_380=T_380,
                          P=P,
                          P_B=P_B,
                          SQRT=FALSE
                          )
  
  ## Add a poissonian noise and multiply the adu by the gain factor G
  adu <- G*rpois(length(adu),adu)
  
  ## Define the different fields of the data frame
  t <- attr(Ca,"Time")
  tOn <- rep(attr(Ca,"tOn"),ntransients) ## unused except as a data frame attribute  
  Time <- rep(c(rep(NA,nb_B), rep(t,ntransients)),2)

  sample_size <- length(adu)/2
  lambda <- c(rep(340,sample_size),rep(380,sample_size))
  
  transient <- rep(c(rep(0,nb_B),
                     as.vector(matrix(rep(1:ntransients,length(t)),
                                      nrow=length(t), byrow=TRUE))),2)
 
  ## Create the data frame
  result <- data.frame(adu=adu,
                       Time=Time,
                       lambda=lambda,
                       transient=transient
                       )
  
  ## Define the class of the dataframe
  class(result) <- "fluo_rawdata"

  ## Strictly necessary attributes
  attr(result,"nb_B") <- nb_B
  attr(result,"tOn") <- tOn
  attr(result,"R_min") <- R_min
  attr(result,"R_max") <- R_max
  attr(result,"K_eff") <- K_eff
  attr(result,"K_d") <- K_d
  attr(result,"B_T") <- B_T
  attr(result,"T_stim") <- c(T_340,T_380)
  attr(result,"P") <- P
  attr(result,"P_B") <- P_B
  attr(result,"G") <- G
  attr(result,"s_ro") <- s_ro

  return(result)
}

