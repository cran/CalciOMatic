ratioExpPhysio <-
function(dataset="inVitro",
                           expe=1,
                           stim=1,
                           idxOn=10,
                           R_min=0.136,
                           R_max=2.701,
                           K_eff=3.637,
                           K_d=0.583,
                           G=0.146,
                           s_ro=16.4,
                           alphamethod=TRUE
                           ) {
  ## Function ratioExpPhysio
  ## Retrieves the data from a ratiometric experiment,
  ## and converts it as a data frame of class "fluo_rawdata"
  ##
  ## dataset: a string with the name of the variable containing the data of ratiometric measurements
  ##          WARNING : the dataset must have a given structure, detailed elsewhere.
  ##          In a future version, the dataset will be tested to be of class "ratio_dataset".
  ##          
  ## expe: the number of the experiment(s) (field "Exp" of the dataset)
  ## stim: the number of the stimulation performed in the chosen experiment
  ## idxOn: the index of the latency at which the laser is set on.
  ##
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
  ## G: the gain factor applied to the number of photons count by the camera.
  ## alphamethod: a logical value, telling whether (or not) we apply the alpha method.
  ##              This could be used for instance to take time-dependant variations of B_T into account.
  ##
  ## A data frame of class "fluo_rawdata" is returned with the following variables:
  ##  - adu: the counts read out of the "camera"
  ##  - Time: the time (in s) at which the adu is read
  ##  - lambda: the excitation wavelength (in nm), a factor
  ##
  ## The returned data frame has also the following attributes:
  ##  - tOn: the time (in s)at which the depolarisation is set
  ##  - T_stim: the stimulation duration (in s) at both wavelengths
  ##  - R_min : copy of arg R_min
  ##  - R_max : copy of arg R_max
  ##  - K_eff : copy of arg K_eff
  ##  - K_d : copy of arg K_d
  ##  - P : the number of pixels used for data binning (size of the Region Of Interest)
  ##  - P_B : the number of pixels used for background data binning
  ##  - G : the gain factor applied to the number of photons count by the camera
  ##  - (opt) B_T: the total fura concentration (in muM)
  ##  - (opt) nbBackM: the number of background measurements performed at each wavelength

  ## Get the whole experimental DataSet
  data <- eval(parse(text=sprintf("%s$Exp%02d",dataset,expe)))
  if(is.null(data)) data <- eval(parse(text=sprintf("%s[[%02d]]",dataset,expe)))
  
  ## Retrieve background measurements from the experimental DataSet
  ## Test if the user used the same number of background measuements at each wavelength  
  nb_B <- min(length(data$adu340Background),length(data$adu380Background))
  adu340_B <- data$adu340Background[1:nb_B]
  adu380_B <- data$adu380Background[1:nb_B]
  
  ## Retrieve fluorescence measurements from the experimental DataSet
  adu340 <- rbind(data$stim1$adu340, data$stim2$adu340, data$stim3$adu340)
  adu380 <- rbind(data$stim1$adu380, data$stim2$adu380, data$stim3$adu380)
  
  ## Retrieve the time sample from the experimental DataSet
  t <- rbind(data$stim1$time, data$stim2$time, data$stim3$time)
  
  ## Set the stimulus onset
  tOn <- t[,idxOn]
  
  ## Retrieve the time exposures at each wavelength,
  ## the number of pixels of the Region Of Interest,
  ## and the intracellular Fura concentration
  T_340 <- data$exposureTime340
  T_380 <- data$exposureTime380
  P <- data$P
  P_B <- data$PBackground
  B_T <- data$furaPipette

  ## Apply the alpha-method on the three transients at 340 and 380 nm to estimate alpha
  ## The lacking fluorescence values at 360 nm are estimated from a spline-smoothed curve of adu360
  ## We consider the time corresponding to the baseline and the last point of each transient
  ## The alpha-method is applied on these 3*(10+1) points to estimate alpha and beta
  ## The value of alpha is then used to calculate the probable B_T time course for the transient of interest
  if(alphamethod==TRUE) {
    with360 <- TRUE

    if(with360 == TRUE) {
      ## Define the signal to smooth (360 nm)
      data_360 <- data$adu360/(data$exposureTime360*data$P360) -
        data$adu360Background/(data$exposureTime360*data$P360Background)
      
      ## We assume a Poissonian noise -> use sqrt
      data_360_ss <- smooth.spline(data$time360,data_360)
      
      ## Interpolate the fluorescence at 360 nm at times corresponding to the
      ## baselines and the last points of the transients recorded at 340 and 380 nm
      time1 <- t[1,c(1:(idxOn-1),dim(t)[2])]
      time2 <- t[2,c(1:(idxOn-1),dim(t)[2])]
      time3 <- t[3,c(1:(idxOn-1),dim(t)[2])]
      data_360_t <- (predict(data_360_ss,c(time1, time2, time3))$y)
      
      ## Define the signal to fit (340 nm)
      data_340_1 <- adu340[1,c(1:(idxOn-1),dim(t)[2])]/(T_340*P) - adu340_B/(T_340*P_B)
      data_340_2 <- adu340[2,c(1:(idxOn-1),dim(t)[2])]/(T_340*P) - adu340_B/(T_340*P_B)
      data_340_3 <- adu340[3,c(1:(idxOn-1),dim(t)[2])]/(T_340*P) - adu340_B/(T_340*P_B)
      data_340_t <- c(data_340_1, data_340_2, data_340_3)
      
      ## Define the signal to fit (380 nm)
      data_380_1 <- adu380[1,c(1:(idxOn-1),dim(t)[2])]/(T_380*P) - adu380_B/(T_380*P_B)
      data_380_2 <- adu380[2,c(1:(idxOn-1),dim(t)[2])]/(T_380*P) - adu380_B/(T_380*P_B)
      data_380_3 <- adu380[3,c(1:(idxOn-1),dim(t)[2])]/(T_380*P) - adu380_B/(T_380*P_B)    
      data_380_t <- c(data_380_1, data_380_2, data_380_3)
      
      ## Fit alpha and beta
      alphamethod_fit <- nls(data_340_t ~ -alpha*data_380_t+beta*data_360_t,
                             start=list(alpha=0,beta=1))
      alpha <- summary(alphamethod_fit)$parameters[1,1]
      beta <- summary(alphamethod_fit)$parameters[2,1]
      
      alpha_init <- alpha
      
      ## -----------------------------------------------------------------
      ## MONTE-CARLO to get a reliable estimation of ^alpha and sd(^alpha)
      ## -----------------------------------------------------------------
      
      ## Smooth the data      
      data_340_1[-length(data_340_1)] <- mean(data_340_1[-length(data_340_1)])
      data_340_2[-length(data_340_2)] <- mean(data_340_2[-length(data_340_2)])
      data_340_3[-length(data_340_3)] <- mean(data_340_3[-length(data_340_3)])
      data_340_t <- c(data_340_1, data_340_2, data_340_3)

      data_380_1[-length(data_380_1)] <- mean(data_380_1[-length(data_380_1)])
      data_380_2[-length(data_380_2)] <- mean(data_380_2[-length(data_380_2)])
      data_380_3[-length(data_380_3)] <- mean(data_380_3[-length(data_380_3)])
      data_380_t <- c(data_380_1, data_380_2, data_380_3)

      data_360_t <- (data_340_t+alpha*data_380_t)/beta
      
      ## - Draw new data from Poisson distribution (take care of G !)
      ## - Fit these new data
      ## - Take the mean and sd of the vector of alpha values obtained from the fit
      MCtrials <- 1000
      
      alpha_vec <- vector("numeric",MCtrials)
      
      for(k in 1:MCtrials) {
        data_340_t_pois <- G*rpois(length(data_340_t),data_340_t/G)
        data_380_t_pois <- G*rpois(length(data_380_t),data_380_t/G)
        data_360_t_pois <- G*rpois(length(data_360_t),data_360_t/G)
        
        alphamethod_fit <- nls(data_340_t_pois ~ -alpha*data_380_t_pois+beta*data_360_t_pois,
                               start=list(alpha=alpha,beta=beta))
        
        alpha_vec[k] <- summary(alphamethod_fit)$parameters[1,1]
      }
      
      alpha <- mean(alpha_vec)
      alpha <- list(value = alpha,
                    mean = alpha,
                    se = sd(alpha_vec),
                    USE_se = TRUE)
    } else {      
      ## In case there are no data recorded at 360 nm...   
      ## Fit each transient separately, or times before tON
      ## and assign alpha the mean of the different values obtained
      alpha_list <- lapply(1:dim(adu340)[1],
                           function(i) {
                             data_340_i <- adu340[i,1:(idxOn-1)]/(T_340*P) - adu340_B/(T_340*P_B)
                             data_380_i <- adu380[i,1:(idxOn-1)]/(T_380*P) - adu380_B/(T_380*P_B)
                             fit_i <- lm(data_340_i ~ data_380_i)
                             val <- summary(fit_i)$coefficients[2,1]
                             return(val)
                           }
                           )
      
      alpha <- mean(unlist(alpha_list))
      alpha <- list(value = alpha,
                    mean = alpha,
                    se = 0,
                    USE_se = !TRUE)
    }
  }
  
  ## Define the dataframe containing all relevant information for future use for fitting  
  adu <- adu340_B
  for(k in stim) adu <- c(adu, adu340[k,])
  adu <- c(adu,adu380_B)
  for(k in stim) adu <- c(adu, adu380[k,])
  
  Time <- rep(NA,nb_B)
  for(k in stim) Time <- c(Time, t[k,])
  Time <- rep(Time,2)

  sample_size <- length(Time)/2
  lambda <- c(rep(340,sample_size),rep(380,sample_size))

  transient <- rep(c(rep(0,nb_B),
                     as.vector(matrix(rep(stim,dim(t)[2]),
                                      nrow=dim(t)[2], byrow=TRUE))),2)
    
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
  
  if(alphamethod == TRUE) {
    attr(result,"alpha") <- alpha
  }
  
  return(result)
}

