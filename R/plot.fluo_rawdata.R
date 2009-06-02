`plot.fluo_rawdata` <-
function(x,
                              y=NULL,
                              numTransient=1,
                              items=1:3,
                              col="black",
                              main="Ratiometric experiment: raw data",
                              xlabs=c("","","Time (s)"),
                              ylabs=c(expression(paste(adu[340]," (photons)")),
                                      expression(paste(adu[380]," (photons)")),
                                      expression(paste("[",Ca^{2+phantom()}, "] (",mu,"M)"))),
                              labs=c(expression(A[1]),expression(A[2]),"B"),
                              ylas=1,
                              ask=FALSE,
                              ...
                              ) {
  ## Function plot.fluo_rawdata
  ## Plot different transients arising from ratiometric experiments
  ## 
  ## x: an object of class "fluo_rawdata" (see the ratioExpSimul function for details)
  ## y: NULL. argument not used
  ## numTransient: a vector of indices indicating which transients to consider
  ## items: a vector of indices between 1 and 3, indicating which signals to plot
  ##        (1 for the fluorescence transient at 340 nm,
  ##         2 for the fluorescence transient at 380 nm,
  ##         3 for the corresponding calcium transient)
  ## col: a character string; defines the color used in the plots
  ## main: a character string; defines the title of the figure
  ## xlabs: a character string; defines the label to apply to the 'x' axis at the bottom of each plot
  ## ylabs: a vector of character strings; defines the labels to apply to the 'y' axes
  ## labs: a vector of character strings; defines the label of each panel (useful for publications)
  ## ylas: an integer, defines the orientation of the yticks
  ## oma: a 4-element vector for the outer margins of the figure
  ## mar: a 4-element vector for the current plot margins
  ## ask: a logical value; set to FALSE to plot all transients on the same figure
  ##                       set to TRUE to plot each transient on a different figure
  ## ...: other parameters that can be set to the plot function
  
  ## -----------------------------------------------
  ## Test the characteristics of the input arguments
  df <- x
  
  if(!inherits(df,"fluo_rawdata")) {
    stop("use plot.fluo_rawdata only with objects of class 'fluo_rawdata'")
  }
  
  ## Manage items to plot
  items <- as.integer(unique(items[items>=1 & items<=3]))
  if(length(items) < 1) {
    stop("There is nothing to plot, since length(items) < 1")
  }
  
  ## Manage labels to add to the x-axis
  if(length(xlabs) < length(items)) {
    xlabs <- c("","","Time (s)")
  } else {
    xlabs0 <- rep("",3)
    if(!missing(xlabs)) {
      for(i in 1:length(items)) xlabs0[items[i]] <- xlabs[i]
    } else {
      for(i in 1:length(items)) xlabs0[items[i]] <- xlabs[items[i]]
    }
    xlabs <- xlabs0
  }

  ## Manage labels to add to the y-axis
  if(length(ylabs) < length(items)) {
    ylabs <- c(expression(paste(adu[340]," (photons)")),
               expression(paste(adu[380]," (photons)")),
               expression(paste("[",Ca^{2+phantom()}, "] (",mu,"M)")))
  } else {
    ylabs0 <- rep("",3)
    if(!missing(ylabs)) {
      for(i in 1:length(items)) ylabs0[items[i]] <- ylabs[i]
    } else {
      for(i in 1:length(items)) ylabs0[items[i]] <- ylabs[items[i]]
    }
    ylabs <- ylabs0
  }
  
  ## Manage labels to add at the top of the panels
  if(length(labs) < length(items)) {
    labs <- c(expression(A[1]),expression(A[2]),"B")
  } else {
    labs0 <- rep("",3)
    if(!missing(labs)) {
      for(i in 1:length(items)) labs0[items[i]] <- labs[i]
    } else {
      for(i in 1:length(items)) labs0[items[i]] <- labs[items[i]]
    }
    labs <- labs0
  }
  
  ## Manage the different transients to work with
  transients <- unique(with(df,transient[transient != 0]))
  if(missing(numTransient)) numTransient <- transients
  transients <- sort(intersect(transients,numTransient))
  
  ## ----------------------------------------------------------
  ## Plot the figures corresponding to the different transients
  if(length(transients) > 1) {
    for(k in transients) {
      ## Call the plot.fluo_rawdata function with a single transient...
      plot(x=df,
           y=NULL,
           numTransient=k,
           items=items,
           col=col,
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
    ## or on the same one with a predefined layout (ask=FALSE)
    if(ask==FALSE & length(items)>1) {
      layout(matrix(1:length(items),nrow=length(items),ncol=1))
    } else {
      layout(1)
    }

    ## Outer and Figure Margins
    oma <- c(4,0,1,0)
    mar <- c(0,7,2,0)
        
    bottom <- ifelse(length(which(xlabs==""))==length(xlabs),0,oma[1])
    if(is.null(main)) {
      top <- 0
    } else if(main == "") {
      top <- 0
    } else {
      top <- oma[3]
    }
    oma <- c(bottom,oma[2],top,oma[4])
    
    par(oma=oma, mar=mar)
    
    
    ## ------------------------------------------------------------
    ## Retrieve relevant data to plot, from the df object
    lambda_vec <- unique(df$lambda)
    Time <- unique(with(df,Time[is.finite(Time) &
                                lambda==lambda_vec[1] &
                                transient==transients]))
    
    ## Build the adu matrix, containing the fluorescence and calcium signals
    L <- length(Time)
    adu <- matrix(data=NA, ncol=length(lambda_vec)+1, nrow=L)
    for(i in 1:length(lambda_vec)) {
      adu[,i] <- with(df,adu[is.finite(Time) &
                             lambda==lambda_vec[i] &
                             transient==transients])
    }
    adu[,length(lambda_vec)+1] <- caFromDf(df,
                                           numTransient=transients,
                                           Plot=FALSE)
    
    if(ask==FALSE) {
      ## Plot all transients on the same figure
      for(i in 1:length(items)) {
        n <- items[i]
        plotCalciOMatic(x=Time,
                        y=adu[,n],
                        n=1,
                        col=col,
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
    } else {
      ## Plot the first item according to the first element of 'items',
      ## and then plot what the user wants: the user still remains the king !
      pos <- 0
      L <- length(items)
      action <- "next"
      
      while(action != "stop") {
        if(action == "next") {
          pos <- pos + 1 - L*(pos==L)
        } else if(action == "previous") {
          pos <- pos - 1 + L*(pos==1)
        } else if(action=="none") {
          ## Do nothing
        }
        
        n <- items[pos]
        action <- plotCalciOMatic(x=Time,
                                  y=adu[,n],
                                  n=1,
                                  col=col,
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
}

