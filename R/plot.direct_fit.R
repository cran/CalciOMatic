`plot.direct_fit` <-
function(x,
                            y=NULL,
                            numTransient=1,
                            items=1:7,
                            col="black",
                            col2="darkgray",
                            main="Fluorescence transients: Direct fit",
                            xlabs=c("Time (s)",
                                    "Time (s)",
                                    "Time (s)",
                                    "Time (s)",
                                    "Lag",
                                    "Theoretical quant.",
                                    "Residuals"),
                            ylabs=c(expression(sqrt(adu[340])),
                                    expression(res[340]),
                                    expression(sqrt(adu[380])),
                                    expression(res[380]),
                                    "ACF",
                                    "Sample quant.",
                                    "Counts"),
                            labs=c(expression(A[1]),
                                   expression(A[2]),
                                   expression(B[1]),
                                   expression(B[2]),
                                   "C","D","E"),
                            ylas=1,
                            ask=FALSE,
                            ...
                            ) {
  ## Function plot.direct_fit
  ## Plot fluorescence transients obtained with a ratiometric marker
  ## (raw data and data fitted with the "direct" method)
  ## as well as different analyses arising from the direct fit
  ##
  ## x: an object of class "direct_fit" (see the directFit function for details)
  ## y: NULL. argument not used
  ## items: a vector of indices between 1 and 5, indicating which transients to plot
  ##        (1 for the fluorescence transient at 340 nm (raw data + fitted data),
  ##         2 for the fit residuals at 340 nm,
  ##         3 for the fluorescence transient at 380 nm (raw data + fitted data),
  ##         4 for the fit residuals at 380 nm,
  ##         5 for the auto-correlation function of the residuals
  ##         6 for the quantile-quantile plot of the residuals
  ##         7 for the residuals histogram)
  ## col:  the color used for the raw data and fit analysis
  ## col2: the color used for the fitted data and useful lines
  ## main: a character string; defines the title of the figure
  ## xlabs: a vector of character strings; defines the label to apply to the 'x' axis at the bottom of each plot
  ## ylabs: a vector of character strings; defines the labels to apply to the 'y' axis on the left of each plot
  ## labs: a vector of character strings; defines the label of each panel (useful for publications)
  ## ylas: an integer, defines the orientation of the yticks
  ## oma: a 4-element vector for the outer margins of the figure
  ## mar: a 4-element vector for the current plot margins
  ## ask: a logical value; set to FALSE to plot all transients on the same figure
  ##                       set to TRUE to plot each transient on a different figure
  ## ...: other parameters that can be set to the plot function

  ## -----------------------------------------------
  ## Test the characteristics of the input arguments
  direct_fit <- x

  if(!inherits(direct_fit,"direct_fit")) {
    stop("use only with \"direct_fit\" objects")
  }

  ## Manage items to plot
  items <- sort(unique(items[items>=1 & items<=7]))
  if(length(items) == 0) {
    stop("There is nothing to plot, since length(items) < 1")
  }
  
  ## Manage labels to add to the x-axis
  if(length(xlabs) < length(items)) {
    xlabs <- c("Time (s)","Time (s)","Time (s)","Time (s)",
               "Lag","Theoretical quant.","Residuals")
  } else {
    xlabs0 <- rep("",7)
    if(!missing(xlabs)) {
      for(i in 1:length(items)) xlabs0[items[i]] <- xlabs[i]
    } else {
      for(i in 1:length(items)) xlabs0[items[i]] <- xlabs[items[i]]
    }
    xlabs <- xlabs0
  }
  
  ## Manage labels to add to the y-axis
  if(length(ylabs) < length(items)) {
    ylabs <- c(expression(sqrt(adu[340])),
               expression(res[340]),
               expression(sqrt(adu[380])),
               expression(res[380]),
               "ACF",
               "Sample quant.",
               "Counts")
  } else {
    ylabs0 <- rep("",7)
    if(!missing(ylabs)) {
      for(i in 1:length(items)) ylabs0[items[i]] <- ylabs[i]
    } else {
      for(i in 1:length(items)) ylabs0[items[i]] <- ylabs[items[i]]
    }
    ylabs <- ylabs0
  }
  
  ## Manage labels to add at the top of the y-axis
  if(length(labs) < length(items)) {
    labs <- c(expression(A[1]),expression(A[2]),
              expression(B[1]),expression(B[2]),
              "C","D","E")
  } else {
    labs0 <- rep("",7)
    if(!missing(labs)) {
      for(i in 1:length(items)) labs0[items[i]] <- labs[i]
    } else {
      for(i in 1:length(items)) labs0[items[i]] <- labs[items[i]]
    }
    labs <- labs0
  }
  
  ## Which transients have been fitted
  fitFunction <- as.character(direct_fit$call$formula)[3]
  Nchar <- nchar(fitFunction)
  TF <- which(sapply(1:(Nchar-7), function(i) substr(fitFunction,i,i+7) == "log_Ca0_"))
  transients <- sapply(1:length(TF), function(i) as.numeric(substr(fitFunction,TF[i]+8,TF[i]+8)))
  transients_backup <- transients
  transients <- intersect(transients,numTransient)
  
  ## ----------------------------------------------------------
  ## Plot the figures corresponding to the different transients
  if(length(transients) > 1) {
    for(k in transients) {
      ## Call the plot.direct function with a single transient...
      plot(x=direct_fit,
           y=NULL,
           numTransient=k,
           items=items,
           col=col,
           col2=col2,
           main=main,
           xlabs=xlabs,
           ylabs=ylabs,
           labs=labs,
           ylas=ylas,
           ask=ask,
           ...
           )
      
      if(k < transients[length(transients)]) {
        ## ... and open a new device (identical to the current one) for the next transient
        L <- dev.cur()
        N <- attr(L,"names")
        if(pmatch("X11",N)) N <- "X11"
        eval(parse(text=sprintf("%s()",N)))
      }
    }
  } else if(length(transients) == 1) {
    ## ----
    ## MAIN
    ## ----
    
    ## ---------------------------------------------------------
    ## Plot the different plots on a different figure (ask=TRUE)
    ## or on the same one, with a predefined layout (ask=FALSE)
    idx1234 <- items[items==1 | items==2 | items==3 | items==4]
    L1234 <- length(idx1234)
    idx567 <- items[items==5 | items==6 | items==7]
    L567 <- length(idx567)
    nrow_items <- items-2*floor(items/2)+1
    
    square <- ifelse(length(items)==1,TRUE,FALSE)
    if(ask==FALSE) {
      if(L1234==0) {
        layout <- matrix(1:length(items),nrow=length(items),ncol=1)
      } else if(L1234==1) {
        if(L567==0) {
          layout <- 1
        } else if(L567==1) {
          layout <- matrix(c(1,2),nrow=2,ncol=1,byrow=TRUE)
        } else if(L567==2) {
          layout <- matrix(c(1,1,2,3),nrow=2,ncol=2,byrow=TRUE)
        } else if(L567==3) {
          square <- TRUE
          layout <- matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=TRUE)
        }
      } else if(L1234==2) {
        if(L567==0) {
          layout <- matrix(c(rep(1,nrow_items[1]),rep(2,nrow_items[2])),ncol=1)
        } else if(L567==1) {
          layout <- matrix(c(rep(1,nrow_items[1]),rep(2,nrow_items[2]),rep(3,2)),ncol=1)
        } else if(L567==2) {
          layout <- matrix(c(rep(1,nrow_items[1]),rep(2,nrow_items[2]),rep(3,2),
                             rep(1,nrow_items[1]),rep(2,nrow_items[2]),rep(4,2)),ncol=2)
        } else if(L567==3) {
          layout <- matrix(c(rep(1,nrow_items[1]),rep(2,nrow_items[2]),rep(3,2),
                             rep(1,nrow_items[1]),rep(2,nrow_items[2]),rep(4,2),
                             rep(1,nrow_items[1]),rep(2,nrow_items[2]),rep(5,2)),ncol=3)
        }
      } else if(L1234==3) {
        if(L567==0) {
          layout <- matrix(c(rep(1,nrow_items[1]),rep(2,nrow_items[2]),rep(3,nrow_items[3])),ncol=1)
        } else if(L567==1) {
          layout <- matrix(c(rep(1,nrow_items[1]),rep(2,nrow_items[2]),rep(3,nrow_items[3]),rep(4,2)),ncol=1)
        } else if(L567==2) {
          layout <- matrix(c(rep(1,nrow_items[1]),rep(2,nrow_items[2]),rep(3,nrow_items[3]),rep(4,2),
                             rep(1,nrow_items[1]),rep(2,nrow_items[2]),rep(3,nrow_items[3]),rep(5,2)),ncol=2)
        } else if(L567==3) {
          layout <- matrix(c(rep(1,nrow_items[1]),rep(2,nrow_items[2]),rep(3,nrow_items[3]),rep(4,2),
                             rep(1,nrow_items[1]),rep(2,nrow_items[2]),rep(3,nrow_items[3]),rep(5,2),
                             rep(1,nrow_items[1]),rep(2,nrow_items[2]),rep(3,nrow_items[3]),rep(6,2)),ncol=3)
        }
      } else if(L1234==4) {
        if(L567==0) {
          layout <- matrix(c(rep(1,2),rep(2,1),rep(3,2),rep(4,1)),ncol=1)
        } else if(L567==1) {
          layout <- matrix(c(rep(1,2),rep(2,1),rep(3,2),rep(4,1),rep(5,3)),ncol=1)
        } else if(L567==2) {
          layout <- matrix(c(rep(1,2),rep(2,1),rep(3,2),rep(4,1),rep(5,3),
                             rep(1,2),rep(2,1),rep(3,2),rep(4,1),rep(6,3)),ncol=2)
        } else if(L567==3) {
          layout <- matrix(c(rep(1,2),rep(2,1),rep(3,2),rep(4,1),rep(5,3),
                             rep(1,2),rep(2,1),rep(3,2),rep(4,1),rep(6,3),
                             rep(1,2),rep(2,1),rep(3,2),rep(4,1),rep(7,3)),ncol=3)
        }
      }
    } else {
      layout <- 1
    }
    layout(layout)
    
    ## Retrieve relevant data to plot, from the calcium_ratio_fit object
    ## Retrieve relevant data to plot, from the fit object
    Name <- attr(direct_fit,"Name")
    Time <- attr(direct_fit,"Time")                  ## plots 1&3
    adu_raw <- attr(direct_fit,"RawData")            ## plots 1&3
    adu_fit <- predict(direct_fit)                   ## plots 1&3
    subset <- attr(direct_fit,"Subset")              ## plots 1&2, 3&4
    adu_res <- residuals(direct_fit)                 ## plots 2&4
    rawDF <- attr(direct_fit,"RawDataFrame")
    nb_B <- attr(rawDF,"nb_B")
    paramNames <- names(coef(direct_fit))

    Time <- Time-min(Time,na.rm=TRUE)
    
    ## How many calibration parameters have been fitted ?
    nb_Rmin <- 0
    nb_Rmax <- 0
    nb_Keff <- 0
    nb_Kd <- 0
    nb_alpha <- 0
    ## nb_B_T <- 0

    for(i in paramNames) {
      nb_Rmin <- nb_Rmin + pmatch("log_R_min",i,nomatch=0)
      nb_Rmax <- nb_Rmax + pmatch("log_R_max",i,nomatch=0)
      nb_Keff <- nb_Keff + pmatch("log_K_eff",i,nomatch=0)
      nb_Kd <- nb_Kd + pmatch("log_K_d",i,nomatch=0)
      nb_alpha <- nb_alpha + pmatch("alpha",i,nomatch=0)
      ## nb_B_T <- nb_B_T + pmatch("B_T_",i,nomatch=0)
    }
    nb_Calib <- nb_Rmin + nb_Rmax + nb_Keff + nb_Kd + nb_alpha ## + nb_B_T
    
    ## Create the Time vector for the fit and residuals
    L_subset_backup <- length(subset)
    if(nb_Calib > 0) {
      subset <- subset[-seq((L_subset_backup-nb_Calib+1),L_subset_backup)]
    }
    subset <- matrix(subset,nrow=2,byrow=TRUE)
    subset <- subset[,seq(nb_B+1,ncol(subset))]
    L_tot <- ncol(subset)
    L <- length(transients_backup)
    idx <- which(transients_backup==transients)
    subset_340 <- subset[1,((idx-1)*L_tot/L+1):(idx*L_tot/L)]
    subset_380 <- subset[2,((idx-1)*L_tot/L+1):(idx*L_tot/L)]

    Time_subset_340 <- Time[subset_340]
    Time_subset_380 <- Time[subset_380]

    ## Create the fit vector
    L_fit_backup <- length(adu_fit)
    if(nb_Calib > 0) {
      adu_fit <- adu_fit[-seq((L_fit_backup-nb_Calib+1),L_fit_backup)]
    }
    adu_fit <- matrix(adu_fit,nrow=2,byrow=TRUE)
    adu_fit <- adu_fit[,seq(nb_B+1,ncol(adu_fit))]
    L_tot <- ncol(adu_fit)
    adu_fit_340 <- adu_fit[1,((idx-1)*L_tot/L+1):(idx*L_tot/L)]
    adu_fit_380 <- adu_fit[2,((idx-1)*L_tot/L+1):(idx*L_tot/L)]
    
    ## Create the residuals vector
    if(nb_Calib > 0) {
      adu_res <- adu_res[-seq((L_fit_backup-nb_Calib+1),L_fit_backup)]
    }
    ## adu_res_calib <- NULL
    ## if(nb_Calib > 0) {
    ##   adu_res_calib <- adu_res[seq((L_fit_backup-nb_Calib+1),L_fit_backup)]
    ## }
    adu_res <- matrix(adu_res,nrow=2,byrow=TRUE)
    adu_res <- adu_res[,seq(nb_B+1,ncol(adu_res))]
    L_tot <- ncol(adu_res)
    adu_res_340 <- adu_res[1,((idx-1)*L_tot/L+1):(idx*L_tot/L)]
    adu_res_380 <- adu_res[2,((idx-1)*L_tot/L+1):(idx*L_tot/L)]   
    
    ## Create the Time, raw adu_340 and adu_380 vectors
    Time_340 <- Time[subset_340[1]:subset_340[length(subset_340)]]
    adu_340 <- adu_raw[subset_340[1]:subset_340[length(subset_340)]]
    adu_380 <- adu_raw[subset_380[1]:subset_380[length(subset_380)]]
                          
    ## Fit quality plots
    adu_res <- adu_res_340
    ACF <- acf(adu_res, lag.max=25, plot=FALSE)
    QQNORM <- qqnorm(adu_res, plot.it=FALSE)
    
    ## Define the x and y, x2 and y2 arrays for the plot
    x_list <- vector("list",7)
    x_list[[1]] <- Time_340
    x_list[[2]] <- Time_subset_340
    x_list[[3]] <- Time_340
    x_list[[4]] <- Time_subset_340
    ## x_list[[6]] <- QQNORM$x
    ## x_list[[7]] <- adu_res_340
    x_list[[7]] <- adu_res
    
    y_list <- vector("list",7)
    y_list[[1]] <- adu_340
    y_list[[2]] <- adu_res_340
    y_list[[3]] <- adu_380
    y_list[[4]] <- adu_res_380
    y_list[[5]] <- ACF
    y_list[[6]] <- adu_res
    
    x2_list <- vector("list",7)
    x2_list[[1]] <- Time_subset_340
    x2_list[[3]] <- Time_subset_340
    
    y2_list <- vector("list",7)
    y2_list[[1]] <- adu_fit_340
    y2_list[[3]] <- adu_fit_380
    ## y2_list[[5]] <- adu_res_340
    y2_list[[5]] <- adu_res
    ## y2_list[[6]] <- adu_res_340
    y2_list[[6]] <- adu_res
    
    ## Outer and Figure Margins
    oma=c(0,0,1,0)
    mar=c(4,6,2,1)

    if(is.null(main) | nchar(main)==0) {
      top <- 0
    } else {
      top <- oma[3]
    }
    oma[3] <- top
    
    par(oma=oma, mar=mar)
    
    if(ask==FALSE) {
      ## Insert all plots in a same figure, following the previously defined layout
      if(L1234<=4 & L567==0) {
        par(oma=c(mar[1],oma[2:4]))
      }
      
      for(i in 1:length(items)) {
        n <- items[i]
        if(L1234 > 0) {
          if(n<max(idx1234)) {
            xlabs[n] <- ""
          }
          if(n<=4) {
            par(mar=c(0,mar[2:4]))
          }
        }
        if(square==TRUE) {
          par(mar=mar)
        } else if(square==FALSE & n>=5) {
          if(L1234>=1) {
            par(mar=c(mar[1:2],mar[3]+mar[1],mar[4]))
          } else {
            par(mar=mar)
          }
        }
        
        n <- items[i]
        t <- ifelse(n>=5,n-2,n-2*(ceiling(n/2)-1))

        if(i != 1) {
          main <- ""
        }
        
        plotCalciOMatic(x=x_list[[n]],
                        y=y_list[[n]],
                        n=t,
                        x2=x2_list[[n]],
                        y2=y2_list[[n]],
                        col=col,
                        col2=col2,
                        main=main,
                        xlab=xlabs[n],
                        ylab=ylabs[n],
                        lab=labs[n],
                        ylas=ylas,
                        oma=oma,
                        mar=mar,
                        ask=ask,
                        ...
                        )
        
        ## Back to the standard margins for the plot
        par(mar=mar)
      }
    } else {   
      ## Plot the first transient according to the first element of 'items',
      ## and then plot what the user wants: the user still remains the king !
      par(mar=mar)
      pos <- 0
      L <- length(items)
      action <- "next"
      while(action != "stop") {
        if(action == "next") {
          pos <- pos - L*(pos==L) + 1
        } else if(action == "previous") {
          pos <- pos + L*(pos==1) - 1
        } else if(action=="none") {
          ## Do nothing
        }
        
        n <- items[pos]
        t <- ifelse(n>=5,n-2,n-2*(ceiling(n/2)-1))
        
        action <- plotCalciOMatic(x=x_list[[n]],
                                  y=y_list[[n]],
                                  n=t,
                                  x2=x2_list[[n]],
                                  y2=y2_list[[n]],
                                  col=col,
                                  col2=col2,
                                  main=main,
                                  xlab=xlabs[n],
                                  ylab=ylabs[n],
                                  lab=labs[n],
                                  ylas=ylas,
                                  oma=oma,
                                  mar=mar,
                                  ask=ask,
                                  ...
                                  )
      }
    }
  } else if(length(transients)==0) {
    warning("The requested transient is not the present in the 'direct_fit' object")
  }
}

