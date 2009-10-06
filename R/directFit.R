directFit <- function(df,
                      transients=c(1,2,3),
                      SQRT=TRUE,
                      ratio=NULL,
                      type="mono",
                      Plot=FALSE,
                      Fit=TRUE,
                      AfterPeak=FALSE,
                      Trace=FALSE,
                      WarnOnly=TRUE
                      ) {
  ## Function directFit
  ## 
  ## Fit the fluorescence signals present in the data frame 'df',
  ## with a direct method. The fit is based on the data generation model
  ## described in the paper ... and on a priori about the "mono-" or "bi-"
  ## exponential decay of the calcium concentration after the light onset.
  ## The CCD acquisition of the fluorescence images being modeled as a Poissonian
  ## process, a way to stabilize the signal variance is to use the square root
  ## transformation. This allows to use classical regression techniques.
  ## 
  ## df: a data frame of class "fluo_rawdata" containing all informative elements
  ##     of the experiment. Its structure is defined in the function RatioExp
  ## transients: the transients (of the data frame) which are to be fitted
  ## SQRT: a logical variable; set to TRUE to use the square root transformation
  ## ratio: EITHER  a list of class "initial_guess", giving IG of calcium variables
  ##                (Ca0, dCa, tau, eventually mu and dtau) from the ratiometric transformation
  ##          OR    a "calcium_fit_ratio" object, whose estimated parameters can be retrieved
  ## type: a character string; indicates the type of calcium decay to fit
  ##       ("mono"exponential or "bi"exponential)
  ## Plot: a logical value; set to TRUE to plot the raw data,
  ##                        the initial guess and the fitted signals
  ## Fit: a logical value; set to TRUE to perform the fit;
  ##                       set to FALSE to estimate an initial guess only
  ## AfterPeak: a logical value; set to TRUE to fit only the convex part of the decay
  ##                             set to FALSE to fit the whole signal
  ## Trace: a logical value; set to TRUE to print the residual sum-of-squares
  ##                         at each iteration of the nls fitting process
  ## WarnOnly: a logical value; set to TRUE to go on even if the fit produced an error
  ## 
  ## The function returns an "nls"/"calcium_ratio_fit" object with the following attributes:
  ##   - Name: "mono-" or "bi-"exponential ratiometric fit
  ##   - Time: a copy of the attribute "Time" of the dataframe df
  ##   - RawData: the raw calcium transient resulting from the ratiometric transformation
  ##   - RawDataFrame: copy of argument df
  ##   - FitFunction: the function used in the nls formula

  ## Required package to smooth adu_340 and adu_380
  ## (needed for the estimation of alpha)
  require(cobs)
  
  ## Extract fluorescence signals from the data frame
  ## adu <- with(df,adu)
  lambda_vec <- unique(df$lambda)
  adu_B_340 <- with(df,adu[is.na(Time) & lambda==lambda_vec[1]])
  adu_B_380 <- with(df,adu[is.na(Time) & lambda==lambda_vec[2]])
  nb_B <- length(adu_B_340)
  df_transients <- unique(with(df,transient))
  if(missing(transients)) transients <- df_transients
  transients <- sort(intersect(df_transients,transients))
  
  ## Construct the t, adu_340 and adu_380 matrices
  rows_names <- c()
  for(k in transients) {
    rows_names <- c(rows_names,sprintf("%d",k))
  }
  
  Ncol <- length(with(df,Time[is.finite(Time) & lambda==lambda_vec[1] & transient==transients[1]]))
  
  t <- matrix(NA, nrow=length(transients), ncol=Ncol, dimnames=list(rows_names,NULL))
  adu_340 <- matrix(NA, nrow=length(transients), ncol=Ncol, dimnames=list(rows_names,NULL))
  adu_380 <- matrix(NA, nrow=length(transients), ncol=Ncol, dimnames=list(rows_names,NULL))
  
  for(k in 1:length(transients)) {
    t[k,] <- with(df,Time[is.finite(Time) & lambda==lambda_vec[1] & transient==transients[k]])
    adu_340[k,] <- with(df,adu[is.finite(Time) & lambda==lambda_vec[1] & transient==transients[k]])
    adu_380[k,] <- with(df,adu[is.finite(Time) & lambda==lambda_vec[2] & transient==transients[k]])
  }
  
  ## Retrieve experimental parameters from the data frame
  tOn <- attr(df,"tOn")[transients]
  names(tOn) <- transients
  T_stim <- attr(df,"T_stim")
  T_340 <- T_stim[1]
  T_380 <- T_stim[2]
  P <- attr(df,"P")
  P_B <- attr(df,"P_B")
  G <- attr(df,"G")
  s_ro <- attr(df,"s_ro")
  
  ## Calibration parameters
  R_min <- attr(df,"R_min")
  if(is.null(R_min$mean)) R_min$mean <- R_min$value
  if(is.null(R_min$USE_se)) R_min$USE_se <- FALSE
  log_R_min <- log(R_min$mean)
  
  R_max <- attr(df,"R_max")
  if(is.null(R_max$mean)) R_max$mean <- R_max$value
  if(is.null(R_max$USE_se)) R_max$USE_se <- FALSE
  log_R_max <- log(R_max$mean)
  
  K_eff <- attr(df,"K_eff")
  if(is.null(K_eff$mean)) K_eff$mean <- K_eff$value
  if(is.null(K_eff$USE_se)) K_eff$USE_se <- FALSE
  log_K_eff <- log(K_eff$mean)

  K_d <- attr(df,"K_d")
  if(is.null(K_d$mean)) K_d$mean <- K_d$value
  if(is.null(K_d$USE_se)) K_d$USE_se <- FALSE
  log_K_d <- log(K_d$mean)
  
  ## B_T - alpha
  if(is.null(alpha <- attr(df,"alpha"))) {
    B_T <- matrix(attr(df,"B_T"), nrow=nrow(adu_340),
                  ncol=ncol(adu_340), dimnames=dimnames(adu_340))
    alpha <- list(value=NA, mean=NA, se=NA, USE_se=FALSE)
  } else {
    B_T <- adu_340/(T_340*P) - adu_B_340/(T_340*P_B) +
           alpha$value * (adu_380/(T_380*P) - adu_B_380/(T_380*P_B))
  }

  alpha_backup <- alpha
  
  ## ---------------------------------------------------------
  ## Get an initial guess for the calcium transient parameters
  
  ## 1) Check that the number of transients corresponds to the
  ##    structure of the input "ratio" argument
  if(length(transients) > 1) {
    if(!inherits(ratio,"ratio_fit_list")) {
      ratio <- NULL
    } else {
      match_transient_ratio <- c()
      for(j in transients) {
        match_transient_ratio <- c(match_transient_ratio, is.null(ratio[[j]]))
      }
      if(max(match_transient_ratio) == 1) {
        ratio <- NULL
      }
    }
  }
  
  if(length(transients) == 1) {
    if(!inherits(ratio,"ratio_fit") |
       !inherits(ratio,"initial_guess") |
       (inherits(ratio,"ratio_fit_list") & transients > length(ratio))) {
      ratio <- NULL
    }
  }
  
  ## 2) Define the appropriate list of initial guesses
  ## names of the parameters
  if(type == "mono") {
    param_names <- c("log_Ca0","log_dCa","log_tau")
  } else if(type =="bi") {
    param_names <- c("log_Ca0","log_dCa","log_tau","mu","log_dtau")
  }
  
  if(inherits(ratio,"initial_guess")) { ## | inherits(ratio,"list")) {
    ig_ratio <- ratio
  } else if(inherits(ratio,"ratio_fit")) {
    ig_ratio <- as.list(ratio$par)
    for(j in 1:length(ig_ratio)) {
      ## names(ig_ratio)[j] <- sprintf("%s_%d",names(ig_ratio)[j],transients)
      names(ig_ratio)[j] <- sprintf("%s_%d",param_names[j],transients)
    }
  } else if(inherits(ratio,"ratio_fit_list")) {
    ig_ratio <- c()
    for(k in transients) {
      ig_ratio_k <- as.list(ratio[[k]]$par)
      for(j in 1:length(ig_ratio_k)) {
        ## names(ig_ratio_k)[j] <- sprintf("%s_%d",names(ig_ratio_k)[j],k)
        names(ig_ratio_k)[j] <- sprintf("%s_%d",param_names[j],k)
      }
      ig_ratio <- c(ig_ratio, ig_ratio_k)
    }
  } else {
    ig_ratio_list <- ratioFitFromDf(df,
                                    transients=transients,
                                    type=type,
                                    Plot=Plot,
                                    Fit=Fit,
                                    AfterPeak=AfterPeak
                                    )
    
    if(inherits(ig_ratio_list,"ratio_fit_list")) {
      ig_ratio <- c()
      for(k in transients) {
        ig_ratio_k <- as.list(ig_ratio_list[[k]]$par)
        for(j in 1:length(ig_ratio_k)) {
          ## names(ig_ratio_k)[j] <- sprintf("%s_%d",names(ig_ratio_k)[j],k)
          names(ig_ratio_k)[j] <- sprintf("%s_%d",param_names[j],k)
        }
        ig_ratio <- c(ig_ratio, ig_ratio_k)
      }
    } else if(inherits(ig_ratio_list,"ratio_fit")) {
      ig_ratio <- as.list(ig_ratio_list$par)
      for(j in 1:length(ig_ratio)) {
        ## names(ig_ratio)[j] <- sprintf("%s_%d",names(ig_ratio)[j],transients)
        names(ig_ratio)[j] <- sprintf("%s_%d",param_names[j],transients)
      }
    }
  }
  class(ig_ratio) <- "initial_guess"
    
  ## -------------------------------------------------------
  ## Define a "mysubset" attribute to the time vector/matrix
  attr(t,"mysubset") <- NULL
  
  ## Optionally: keep only the convex/concave part of the adu transients (after the peak)
  if(AfterPeak >= 1) {
    idxConvex <- 0
    for(k in 1:length(transients)) {
      idxOn <- which(t[k,] >= tOn[k])[1]
      if(is.logical(AfterPeak)) {
        idxConvex_340 <- transientConvexPart(t=t[k,],tOn=tOn[k],transient=adu_340[k,])
        idxConvex_380 <- transientConvexPart(t=t[k,],tOn=tOn[k],transient=adu_380[k,])
        idxConvex_k <- max(c(idxConvex_340,idxConvex_380))
      } else if(is.numeric(AfterPeak)) {
        idxConvex_k <- idxOn+AfterPeak-1
      }
      idxConvex <- max(idxConvex_k,idxConvex)
    }
    if(idxConvex > idxOn) attr(t,"mysubset") <- idxOn:(idxConvex-1)
  }
  
  ## Define the subset of time vector to use for the initial guess
  mysubset <- 1:dim(t)[2]
  try(mysubset <- mysubset[-attr(t,"mysubset")],silent=TRUE)

  ## ---------------------------------------------------------------
  ## Get an initial guess for the remaining parameters of the model:
  ## log_S_B_340, log_S_B_380 and log_phi
  ig <- igDirect(adu_B_340=adu_B_340,
                 adu_340=adu_340,
                 adu_B_380=adu_B_380,
                 adu_380=adu_380,
                 ig_ratio=ig_ratio,
                 t=t,
                 tOn=tOn,
                 subset=mysubset,
                 R_min=R_min$mean,
                 R_max=R_max$mean,
                 K_eff=K_eff$mean,
                 K_d=K_d$mean,
                 B_T=B_T,
                 T_340=T_340,
                 T_380=T_380,
                 P=P,
                 P_B=P_B
                 )
  
  ## Optionnally add initial guesses for the calibration parameters
  ## (log_R_min, log_R_max, log_K_eff and log_K_d) and alpha
  if(R_min$USE_se) ig <- c(ig,log_R_min=log(R_min$mean))
  if(R_max$USE_se) ig <- c(ig,log_R_max=log(R_max$mean))
  if(K_eff$USE_se) ig <- c(ig,log_K_eff=log(K_eff$mean))
  if(K_d$USE_se) ig <- c(ig,log_K_d=log(K_d$mean))
  if(alpha$USE_se) ig <- c(ig,alpha=alpha$mean)
  
  ## Put the variables of ig in the environment of the current function
  for(i in 1:length(ig)) assign(names(ig)[i],ig[[i]])
  
  ## ----------------------------------------------------
  ## Define the function to fit with the "nls" algorithm:
  ## This function depends both on the type of calcium decay,
  ## and of the user choice to consider (or not) the variability
  ## of the parameters supplied by the calibration experiments
  ## (R_min,R_max,K_eff and K_d) and that of the alpha coefficient
  
  ## (1) Define a function prediciting the whole adu transients
  ##     (background and fluorescence at 340 and 380 nm concatenated);
  ##     => function of 13 (resp. 15) variables depending on the type of decay
  function_4_direct_fit <- mkFunction4DirectFit(type=type,
                                                nb_B=nb_B,
                                                transients=transients,
                                                alphamethod=alpha_backup$USE_se,
                                                SQRT=SQRT
                                                )
  
  ## (2) Define the final 'adu' vector to fit and the formula to use
  ##     Define also a vector of weights and the subset
  
  ## Initial adu to fit
  ## /|\  SQUARE ROOT TRANSFORMATION
  adu <- c(adu_B_340, as.vector(t(adu_340)),
           adu_B_380, as.vector(t(adu_380))
           )
  if(SQRT == TRUE) adu <- sqrt(adu + G*s_ro^2)
  adu_backup <- adu
  
  ## Initial time
  t2 <- rep(c(with(df,Time[!is.finite(Time) & lambda==340 & transient==0]),
              as.vector(t(t))
              ),2)
  
  ## Initial weights: If the square root transformation is used,
  ##                  the noise level is sigma=sqrt(G)/2 (photon),
  ##                  where G is the gain factor
  weights <- rep(4/G,length(adu))
  
  ## Initial function for the fit
  function_4_direct_fit_2 <- function() NULL
  
  ## ... Formals
  formals(function_4_direct_fit_2) <- formals(function_4_direct_fit)
  Formals <- names(formals(function_4_direct_fit))
  L <- length(Formals)
  
  ## ... Body
  string4Fit <- "c(function_4_direct_fit("
  for(i in 1:(L-1)) string4Fit <- c(string4Fit,Formals[i],",")
  string4Fit <- c(string4Fit,Formals[L],")")
  
  ## Adding constraints on calibration parameters and alpha
  
  ## R_min
  if(R_min$USE_se) {
    adu <- c(adu,R_min$mean)
    t2 <- c(t2,NA)
    weights <- c(weights,1/(R_min$se^2))
    string4Fit <- c(string4Fit,",","exp(log_R_min)")
  }
  
  ## R_max
  if(R_max$USE_se) {
    adu <- c(adu,R_max$mean)
    t2 <- c(t2,NA)
    weights <- c(weights,1/(R_max$se^2))
    string4Fit <- c(string4Fit,",","exp(log_R_max)")
  }
  
  ## K_eff
  if(K_eff$USE_se) {
    adu <- c(adu,K_eff$mean)
    t2 <- c(t2,NA)
    weights <- c(weights,1/(K_eff$se^2))
    string4Fit <- c(string4Fit,",","exp(log_K_eff)")
  }
  
  ## K_d
  if(K_d$USE_se) {
    adu <- c(adu,K_d$mean)
    t2 <- c(t2,NA)
    weights <- c(weights,1/(K_d$se^2))
    string4Fit <- c(string4Fit,",","exp(log_K_d)")
  }
  
  if(alpha_backup$USE_se) {
    adu <- c(adu,alpha_backup$mean)
    t2 <- c(t2,NA)
    weights <- c(weights,1/(alpha_backup$se^2))
    string4Fit <- c(string4Fit,",","alpha")
  }
  
  string4Fit <- c(string4Fit,")")
  body(function_4_direct_fit_2) <- parse(text=string4Fit)
  
  ## Define a "mysubset" vector for the background and fluorescence data concatenated
  mysubset_B <- c(1:length(adu_B_340))
  mysubset_340 <- mysubset_B
  for(k in 1:dim(adu_340)[1]) {
    mysubset_340 <- c(mysubset_340, length(mysubset_B) + mysubset + (k-1)*dim(adu_340)[2])
  }
  mysubset_380 <- mysubset_340 + (length(mysubset_B) + dim(adu_340)[1]*dim(adu_340)[2])
  mysubset <- c(mysubset_340, mysubset_380)
  
  ## Define the subset of indices to use for the fit
  if(length(weights) > length(adu_backup)) mysubset <- c(mysubset,(length(adu_backup)+1):length(weights))
  
  ## Plot the raw data
  if(Plot==TRUE) {plot(adu,type="l")}
  
  ## Smooth the adu_340 and adu_380 transients to introduce
  ## adu_340+alpha*adu_380 in the fitting process
  if(is.finite(alpha_backup$value)) {
    for(i in 1:dim(adu_340)[1]) {
      t_fit <- 1:(idxOn-1)
      y_fit <- adu_340[i,1:(idxOn-1)]
      LinFit <- lm(y_fit ~ t_fit)
      adu_340[i,1:(idxOn-1)] <- rep(mean(adu_340[i,1:(idxOn-1)]),idxOn-1)
      adu_340[i,idxConvex:ncol(t)] <- cobs(idxConvex:ncol(t),
                                           adu_340[i,idxConvex:ncol(t)],
                                           nknots=10,
                                           constraint="none")$fitted

      t_fit <- 1:(idxOn-1)
      y_fit <- adu_380[i,1:(idxOn-1)]
      LinFit <- lm(y_fit ~ t_fit)
      adu_380[i,1:(idxOn-1)] <- rep(mean(adu_380[i,1:(idxOn-1)]),idxOn-1)
      adu_380[i,idxConvex:ncol(t)] <- cobs(idxConvex:ncol(t),
                                           adu_380[i,idxConvex:ncol(t)],
                                           nknots=10,
                                           constraint="none")$fitted
    }
  }
  
  if(alpha_backup$USE_se == TRUE) {
    B_T <- adu_340/(T_340*P) - mean(adu_B_340)/(T_340*P_B) +
           alpha * (adu_380/(T_380*P) - mean(adu_B_380)/(T_380*P_B))
  }
    
  ## -----------------------------------------------------
  ## Define the function for the initial guess and the fit
  myfunction <- "function_4_direct_fit_2("
  for(i in 1:(L-1)) myfunction <- c(myfunction,Formals[i],",")
  myfunction <- c(myfunction,Formals[L],")")
  
  ## Plot the adu corresponding to the parameters initial guess
  if(Plot==TRUE) {
    Initial_Transient <- c()
    string4IT <- c("Initial_Transient <- ",myfunction)
    eval(parse(text=string4IT))
    lines((1:length(adu))[-mysubset],Initial_Transient[-mysubset],col="blue",type="p",pch=20)
    lines((1:length(adu))[mysubset],Initial_Transient[mysubset],col="blue",lwd=1)
  }
  
  ## Perform the fit and plot the predicted data
  if(Fit==TRUE) {
    ## Define the string to evaluate to perform the fit
    if(length(mysubset) <= length(adu)) {
      String_4_Fit <- c("direct_fit <- nls(adu[mysubset] ~ ",myfunction)
      String_4_Fit <- c(String_4_Fit,"[mysubset],start=ig,weights=weights[mysubset],")
      String_4_Fit <- c(String_4_Fit,"trace=Trace,control=list(warnOnly=WarnOnly))")
    }
    
    ## Perform the fit    
    eval(parse(text=String_4_Fit))
    
    if(Plot==TRUE) lines((1:length(adu))[mysubset],predict(direct_fit),col="red",lwd=2)
    
    attr(direct_fit,"Name") <- sprintf("%sexponential direct fit",type)
    attr(direct_fit,"Time") <- t2
    attr(direct_fit,"RawData") <- adu
    attr(direct_fit,"RawDataFrame") <- df
    attr(direct_fit,"FitFunction") <- myfunction
    attr(direct_fit,"Subset") <- mysubset
    
    class(direct_fit) <- c("direct_fit","nls")
    
    result <- direct_fit
  }
  
  return(result)
}
