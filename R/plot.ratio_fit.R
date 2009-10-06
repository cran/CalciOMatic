plot.ratio_fit <- function(x,
                           y=NULL,
                           items=1:5,
                           col="black",
                           col2="darkgray",
                           main="Intracellular calcium transient: Ratiometric fit",
                           xlabs=c("Time (s)",
                                   "Time (s)",
                                   "Lag",
                                   "Theoretical quant.",
                                   "Residuals"),
                           ylabs=c(expression(paste("[",Ca^{2+phantom()}, "] (",mu,"M)")),
                                   expression(res[Ca]),
                                   "ACF",
                                   "Sample quant.",
                                   "Counts"),
                           labs=c(expression(A[1]),
                                  expression(A[2]),
                                  "B","C","D"),
                           ylas=1,
                           ask=FALSE,
                           ...
                           ) {
  ## Function plot.ratio_fit
  ## Plot an intracellular calcium transient, with raw and fitted data,
  ## as well as different analyses arising from the ratiometric fit
  ##
  ## x: an object of class "calcium_ratio_fit" (see the Ratio_Fit_fct function for details)
  ## y: NULL. argument not used
  ## items: a vector of indices between 1 and 5, indicating items transients to plot
  ##        (1 for the intracellular calcium transient (raw data + fitted data),
  ##         2 for the fit residuals,
  ##         3 for the auto-correlation function of the residuals
  ##         4 for the quantile-quantile plot of the residuals
  ##         5 for the residuals histogram)
  ## col: the color used for the raw data and fit analysis
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
  ratio_fit <- x
  
  if(!inherits(ratio_fit,"ratio_fit")) {
    stop("use only with \"ratio_fit\" objects")
  }

  ## Manage items to plot
  items <- sort(unique(items[items>=1 & items<=5])) 
  if(length(items) < 1) {
    stop("There is nothing to plot, since length(items) < 1")
  }
  
  ## Manage labels to add to the x-axis
  if(length(xlabs) < length(items)) {
    xlabs <- c("Time (s)","Time (s)","Lag","Theoretical quant.","Residuals")
  } else {
    xlabs0 <- rep("",5)
    if(!missing(xlabs)) {
      for(i in 1:length(items)) xlabs0[items[i]] <- xlabs[i]
    } else {
      for(i in 1:length(items)) xlabs0[items[i]] <- xlabs[items[i]]
    }
    xlabs <- xlabs0
  }
    
  ## Manage labels to add to the y-axis
  if(length(ylabs) < length(items)) {
    ylabs <- c(expression(paste("[",Ca^{2+phantom()}, "] (",mu,"M)")),
               expression(res[Ca]),
               "ACF",
               "Sample quant.",
               "Counts")
  } else {
    ylabs0 <- rep("",5)
    if(!missing(ylabs)) {
      for(i in 1:length(items)) ylabs0[items[i]] <- ylabs[i]
    } else {
      for(i in 1:length(items)) ylabs0[items[i]] <- ylabs[items[i]]
    }
    ylabs <- ylabs0
  }
  
  ## Manage labels to add at the top of the panels
  if(length(labs) < length(items)) {
    labs <- c(expression(A[1]),expression(A[2]),"B","C","D")
  } else {
    labs0 <- rep("",5)
    if(!missing(labs)) {
      for(i in 1:length(items)) labs0[items[i]] <- labs[i]
    } else {
      for(i in 1:length(items)) labs0[items[i]] <- labs[items[i]]
    }
    labs <- labs0
  }

  
  ## ---------------------------------------------------------
  ## Plot the different plots on a different figure (ask=TRUE)
  ## or on the same one, with a predefined layout (ask=FALSE) 
  R1and2 <- FALSE
  if(ask==FALSE & length(items) > 1){
    if(items[1]==1 && items[2]==2) {
      R1and2 <- TRUE
      if(length(items)==2) {
        layout(matrix(c(1,1,2),nrow=3,ncol=1))
      } else if(length(items) == 3) {
        layout(matrix(c(1,1,2,3,3),nrow=5,ncol=1))
      } else if(length(items)==4) {
        layout(matrix(c(1,1,2,3,3,1,1,2,4,4),nrow=5,ncol=2))
      } else if(length(items)==5) {
        layout(matrix(c(1,1,2,3,3,1,1,2,4,4,1,1,2,5,5),nrow=5,ncol=3))
      }
    } else if(length(items)<=3) {
      layout(matrix(1:length(items),nrow=length(items),ncol=1))
    } else if(length(items)==4) {
      layout(matrix(c(1,3,2,4),nrow=2,ncol=2))
    }
  } else {
    layout(1)
  }
  
  ## Retrieve relevant data to plot, from the ratio_fit object
  Name <- attr(ratio_fit,"Name")
  Time <- attr(ratio_fit,"Time")                  ## plot 1
  tOn <- attr(ratio_fit,"tOn")                    ## plot 1
  tOn <- tOn-min(Time)                            ## plot 1
  Time <- Time-min(Time)                          ## plot 1
  Ca_raw <- attr(ratio_fit,"RawData")             ## plot 1
  Fit_fct <- attr(ratio_fit,"FitFunction")        ## plot 1
  subset <- attr(ratio_fit,"Subset")              ## plots 1&2
  if(is.null(subset)) subset <- 1:length(Time)
  Param <- ratio_fit$par                          ## plot 1
  if(length(Param) == 3) {                        ## plot 1
    Ca_fit <- Fit_fct(Time[subset], tOn, Param[1], Param[2], Param[3])
  } else if(length(Param) == 5) {
    Ca_fit <- Fit_fct(Time[subset], tOn, Param[1], Param[2], Param[3], Param[4], Param[5])
  }
  Ca_res <- Ca_raw[subset] - Ca_fit               ## plot 2
  ACF <- acf(Ca_res,lag.max=25,plot=FALSE)        ## plot 3
  QQNORM <- qqnorm(Ca_res,plot.it=FALSE)          ## plot 4

  ## Define the x and y, x2 and y2 arrays for the plot
  x_list <- vector("list", 5)
  x_list[[1]] <- Time
  x_list[[2]] <- Time[subset]
  x_list[[4]] <- QQNORM$x
  x_list[[5]] <- Ca_res
  
  y_list <- vector("list", 5)
  y_list[[1]] <- Ca_raw
  y_list[[2]] <- Ca_res
  y_list[[3]] <- ACF
  y_list[[4]] <- QQNORM$y

  x2_list <- vector("list", 5)
  x2_list[[1]] <- Time[subset]

  y2_list <- vector("list", 5)
  y2_list[[1]] <- Ca_fit
  y2_list[[3]] <- Ca_res
  y2_list[[4]] <- Ca_res
    
  ## Outer and Figure Margins
  oma <- c(0,0,1,0)
  mar <- c(4,7,2,1)

  if(is.null(main)) {
    top <- 0
  } else if(main=="") {
    top <- 0
  } else {
    top <- oma[3]
  }
  oma[3] <- top
  
  par(oma=oma, mar=mar)
  
  if(ask==FALSE) {
    ## Insert all plots in a same figure, following the previously defined layout
    for(i in 1:length(items)) {
      n <- items[i]
      if(R1and2==TRUE & n<=2 & max(items)>=3) {
        par(mar=c(0,mar[2:4]))
        if(n==1) {
          xlabs[n] <- ""
        }
      }

      if(n>=3 & min(items)<3) {
        par(mar=c(mar[1:2],mar[3]+mar[1],mar[4]))
      }

      plotCalciOMatic(x=x_list[[n]],
                      y=y_list[[n]],
                      n=n,
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
    ## (1) Change the xlabs[1], so that a x-axis is plotted
    if(items[1]==1 & xlabs[1]=="") {
      if(items[2]==2) {
        xlabs[1] <- xlabs[2]
      }
      else {
        xlabs[1] <- "Time (s)"
      }
    }
    
    ## (2) Plot the first item according to the first element of 'items',
    ## and then plot what the user wants: the user still remains the king !
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
      action <- plotCalciOMatic(x=x_list[[n]],
                                y=y_list[[n]],
                                n=n,
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
}
