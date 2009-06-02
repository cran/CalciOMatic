`plotCalciOMatic` <-
function(x=NULL,
                            y=NULL,
                            n=1,
                            x2=NULL,
                            y2=NULL,
                            col="black",
                            col2="darkgray",
                            main="MyCalciumRatiometricFit",
                            xlab="",
                            ylab="",
                            lab="A",
                            ylas=1,
                            oma=c(4,0,1,0),
                            mar=c(0,7,2,0),
                            ask=FALSE,
                            ...
                            ) {
  ## Function plot.calciomatic
  ## Perform different kinds of plot according to the value of n (an integer)
  ##
  ##  - n=1: plot y vs. x (type lines and color col),
  ##         and superimpose the plot of y2 vs. x2 (color col2)
  ##         Used to plot raw data and fitted data
  ##  - n=2: plot y vs. x (type lines and color col),
  ##         and add a dashed horizontal line at y=0 (color col2)
  ##         Used to plot fit residuals
  ##  - n=3: plot a bar plot, with y of class "acf" (color col),
  ##         and add horizontal lines at y=+/-1.96/sqrt(length(y$acf)) (color col2)
  ##         Used to plot the auto-correlation function of the residuals
  ##  - n=4: plot y vs. x (type points and color col),
  ##         and add a diagonal dashed line (color col2)
  ##         Used to plot quantile-quantile plots of the fit residuals
  ##  - n=5: plot an histogram of x (color col), and add vertical dashed lines
  ##         at -3/-2/-1/0/1/2/3 times the standard deviation of x (color col2)
  ##         Used to plot the histogram of the fit residuals
  ##
  ## col: the color of the main signal
  ## col2: the color of the secondary signal
  ## main: a character string to write at the top of the figure
  ## xlab: a character string to write below the xaxis
  ## ylab: a character string to write at the left of the yaxis
  ## lab: a character string to write at the top-left of the panel
  ## ylas: the orientation of the yticks
  ## oma: a 4-element vector for the outer margins of the figure
  ## mar: a 4-element vector for the current plot margins
  ## ask: a logical value; if set to TRUE, symbols are plotted at the bottom right of the figure
  ##      and a string "action" is returned at the first time the user clicks on the figure.
  ##      This action can be "previous", "next", "save", "stop" or "none". Useful when
  ##      plot.calciomatic is called from plot.fluorawdata, plot.calcium_ratio_fit or plot.direct_fit
  ## ...: other parameters that can be set to the plot function
  
  ## Type of data to plot
  n <- as.integer(unique(n[n>=1 & n<=5]))
  if(length(n) != 1) {
    stop("There first argument of the plot.calciomatic function must be an integer between 1 and 5")
  }
  
  ## ---------------
  ## Plot parameters
  ## They are given by the user, through the ... way
  cex <- 1.6
  cex.axis <- 1.2
  cex.lab <- 1.1
  cex.main <- 1.8
  font <- 2
  font.axis <- 1
  font.lab <- 1
  font.main <- 2
  line.xlab <- max(c(oma[1],mar[1]))-1.5
  line.ylab <- mar[2]-2
  line.lab <- mar[2]-0.5
  line.main <- oma[3]-1.5
  adj.main <- par("adj")
  xlim <- Inf
  ylim <- Inf
  tcl <- par("tcl")
  mgp.x <- par("mgp")
  mgp.y <- par("mgp")
  
  mycall <- match.call()
  mynames <- names(mycall)
  var2check <- c("cex","cex.axis","cex.lab","cex.main",
                 "font","font.axis","font.lab","font.main",
                 "line.xlab","line.ylab","line.lab","line.main",
                 "adj.main",
                 "xlim","ylim","tcl","mgp.x","mgp.y")
  varNominalValues <- list(1.6, 1.2, 1.1, 1.8,
                           2, 1, 1, 2,
                           max(c(oma[1],mar[1]))-1.5, mar[2]-2, mar[2]-0.5, oma[3]-1.5,
                           par("adj"),
                           Inf, Inf, par("tcl"), par("mgp"), par("mgp"))
  for(k in 1:length(var2check)) {
    if(length(which(mynames==var2check[k])) == 1) {
      T <- sprintf("%s <- eval(mycall[['%s']])",var2check[k],var2check[k])
      eval(parse(text=T))
    } else {
      assign(var2check[k],varNominalValues[[k]])
    }
  }
    
  ## ------------
  ## Compute xlim
  if(xlim == Inf) {
    if(n==1) {
      ## Raw data
      xlim <- c(min(x),max(x))
    } else if(n==2) {
      ## Residuals
      xlim <- c(min(x),max(x))
    } else if(n==3) {
      ## Auto-correlation function
      xlim <- c(0,length(y$acf)-1)
    } else if(n==4) {
      ## Quantile-Quantile plot
      xlim <- c(-3,3)
    } else if(n==5) {
      ## Residuals histogram
      delta_x <- max(abs(x))
      x_max <- delta_x+0.1*delta_x
      x_max <- signif(x_max,2)
      xlim <- c(-x_max,x_max)
    }
  }
  
  ## ------------
  ## Compute ylim
  if(ylim == Inf) {
    if(n==1) {
      ## Raw data
      delta_y <- max(y)-min(y)
      ylim <- c(min(y)-0.05*delta_y,max(y)+0.05*delta_y)    
    } else if(n==2) {
      ## Residuals
      delta_y <- max(abs(y))
      y_max <- delta_y+0.1*delta_y
      y_max <- signif(y_max,1)
      ylim <- c(-y_max,y_max)      
    } else if(n==3) {
      ## Auto-correlation function
      ylim <- c(min(0,floor(10*min(y$acf))/10,-1.96/sqrt(length(y2))),1)
    } else if(n==4) {
      ## Quantile-Quantile plot
      delta_y <- max(abs(y))
      y_max <- delta_y+0.1*delta_y
      y_max <- signif(y_max,1)
      ylim <- c(-y_max,y_max)
    } else if(n==5) {
      ## Residuals histogram
      x_mean <- mean(x)
      x_sd <- sd(x)
      x_iqr <- IQR(x)
      
      N1 <- ceiling((max(xlim)-x_mean)/x_sd)
      
      ## Choose the bin size according to Scott's formula (1979)
      ## See "Modern Applied Statistics with S", 4th edition
      N2 <- 3.5*(length(x)^(-1/3))
      N2_vec <- 1/(1:20)
      N2 <- N2_vec[which.min(abs(N2_vec-N2))]
      
      ## Choose the bin size according to Freedman-Diaconis's formula (1979)
      ## See "Modern Applied Statistics with S", 4th edition
      N3 <- 2*x_iqr/x_sd*(length(x)^(-1/3))
      N3_vec <- 1/(1:20)
      N3 <- N3_vec[which.min(abs(N3_vec-N3))]
      
      ## Choose the bin size manually
      N4 <- 1/4
      
      breaks <- x_mean + seq(-N1-0.5*N2, N1+0.5*N2, N2) * x_sd
      breaks <- x_mean + seq(-N1-0.5*N3, N1+0.5*N3, N3) * x_sd
      breaks <- x_mean + seq(-N1-0.5*N4, N1+0.5*N4, N4) * x_sd
      
      H <- hist(x, breaks=breaks, plot=FALSE)
      ylim <- c(0, max(H$counts))
    }
  }
  
  ## ---------
  ## Main plot
  ## Define a string which depends on the value of n
  if(n==1) {
    String_4_Plot <- "plot(x, y, type='l', col=col,"
  } else if(n==2) {
    subset <- 1:length(x)
    subset2 <- c()
    x_diff <- diff(x)
    limit <- which(x_diff > 1.01*min(x_diff))
    if(length(limit) > 0) {
      subset <- 1:limit
      subset2 <- (limit+1):length(x)
    }
    String_4_Plot <- "plot(x[subset], y[subset], type='l', col=col, lwd=par('lwd'),"
  } else if(n==3) {
    String_4_Plot <- "plot(y, col='white',"
  } else if(n==4) {
    y <- sort(y)
    sampleSize <- length(y)
    proba <- ppoints(sampleSize)
    x <- qnorm(proba, 0, 1)
    String_4_Plot <- "plot(x, y, col='white', pch=20, cex=cex.lab,"
  } else if(n==5) {
    String_4_Plot <- "H <- hist(x, breaks=breaks, border=col,"
  }
 
  String_4_Plot <- c(String_4_Plot,
                     "main='',
                      xlab='',
                      xlim=xlim,
                      ylab='',
                      ylim=ylim,
                      axes=FALSE,
                      ann=FALSE)"
                     )
  
  ## Perform the main plot
  eval(parse(text=String_4_Plot))

  ## -----------
  ## Decorations
  if(n==1) {
    ## Superimpose  the fitted data
    if(!is.null(x2) & !is.null(y2)) {
      subset <- 1:length(x2)
      subset2 <- c()
      x2_diff <- diff(x2)
      limit <- which(x2_diff > 1.01*min(x2_diff))
      if(length(limit) > 0) {
        subset <- 1:limit
        subset2 <- (limit+1):length(x2)
      }
      lines(x2[subset],
            y2[subset],
            col=col,
            lwd=3*par("lwd")
            )
      try(lines(x2[subset2],
                y2[subset2],
                col=col,
                lwd=3*par("lwd")),
          silent=TRUE
          )
    }
    
  } else if(n==2) {
    lines(x[subset2], y[subset2], col=col, lwd=par("lwd"))
    ## Add a horizontal dotted line
    abline(h=0,
           col=col,
           lty=2
           )
    
  } else if(n==3) {
    ## Plot the confidence interval
    L <- y$n.used    
    factor <- qt(0.975,L)    
    polygon(c(-5,length(y$acf)+2,length(y$acf)+2,-5),
            c(-1,-1,1,1)*factor/sqrt(L),
            col=col2,
            border=col2,
            lty=0
            )    
    lines(y$lag,
          y$acf,
          type="h",
          col=col)    
    abline(h=0, col=col)
    
  } else if(n==4) {
    ## Add the 95% confidence interval on the quantiles
    CI <- 0.95
    ci_0.95 <- qnorm(qbeta((1-CI)/2, 1:sampleSize, sampleSize-(1:sampleSize)+1), 0, 1)
    ci0.95 <- qnorm(qbeta(1-(1-CI)/2, 1:sampleSize, sampleSize-(1:sampleSize)+1), 0, 1)
    polygon(c(x,rev(x)),
            c(ci_0.95*sd(y),rev(ci0.95*sd(y))),
            col=col2,
            border=col2
            )
    points(x, y, col=col, pch=20)
    ## Add the qqline
    abline(a=0,
           b=sd(y),
           lty=2
           )

  } else if(n==5) {
    ## Add vertical lines corresponding to the mean, mean+/-123*sd, the qqline
    for(k in -3:3) {
      lines(x=c(x_mean+k*x_sd,x_mean+k*x_sd),
            y=c(0,H$counts[which(H$breaks > x_mean+k*x_sd)[1]-1]),
            col=col,
            lty=2
            )
    }
  }

  ## Add automatic X ticks and a X label
  if(xlab != "") {
    ## X-Axis with X-ticks
    par(las=1);
    axis(side=1,
         cex.axis=cex.axis,
         font.axis=font.axis,
         tcl=tcl,
         mgp=mgp.x,
         lwd=par("lwd")
         )
    ## X-label
    par(las=1);
    mtext(xlab,
          side=1,
          line=line.xlab,
          cex=cex.lab,
          font=font.lab
          )
  }
  
  ## Add automatic Y ticks and a Y label
  ## (1) Compute the ytick location, if necessary
  if(n==1) {
    at <- NULL    
  } else if(n==2) {
    yat <- signif(y_max,1)
    ifelse(yat > y_max, yat <- signif(0.9*y_max,1), yat <- yat)
    at <- seq(-yat,yat,yat)
  } else if(n==3) {
    at <- seq(floor(5*min(y$acf))/5,1,0.2)
  } else if(n==4) {
    yat <- signif(y_max,1)
    ifelse(yat > y_max, yat <- signif(0.9*y_max,1), yat <- yat)
    at <- seq(-yat,yat,yat/2)
  } else if(n==5) {
    at <- NULL
  }
  
  ## (2) Add the yticks to the plot and write the ylabel
  ## Y-axis with Y-ticks
  par(las=ylas);
  axis(side=2,
       at=at,
       cex.axis=cex.axis,
       font.axis=font.axis,
       tcl=tcl,
       mgp=mgp.y,
       lwd=par("lwd")
       )
  ## Y-labels
  par(las=0);
  mtext(ylab,
        side=2,
        line=line.ylab,
        cex=cex.lab,
        font=font.lab
        )
  
  ## (3) Add a label at the left-top of the plot
  par(las=1);
  mtext(lab,
        side=2,
        line=line.lab,
        at=max(ylim)+0.05*diff(ylim),
        adj=0,
        cex=cex,
        font=font
        )
  
  ## (4) Add a main title to the figure
  par(las=1);
  title(main=main,
        line=line.main,
        outer=TRUE,
        cex.main=cex.main,
        font.main=font.main,
        adj=adj.main
        )
  
  ## mtext(text=main,
  ##       side=3,
  ##       line=line.main, ## oma[3]-1.5,
  ##       outer=FALSE,
  ##       cex=cex.main,
  ##       font=font.main,
  ##       )
  
  ## (5) Add symbols to go from one plot to another (if ask is set to TRUE)
  if(ask==TRUE) {
    X <- xlim[1] + c(0.91,0.95,0.99)*diff(xlim)
    Y <- ylim[1] - 0.02*diff(ylim)
    lines(X[1], Y, type="p", pch="<", lwd=1, cex=2)
    lines(X[2], Y, type="p", pch=15,  lwd=2, cex=1.5)
    lines(X[3], Y, type="p", pch=">", lwd=1, cex=2)
    
    coord <- unlist(locator(1))
    ACTION <- c("previous","stop","next")
    
    distance <- abs(X-coord[1])
    
    if(coord[1] > min(X)-0.02*diff(xlim) & coord[1] < max(X)+0.02*diff(xlim)
       & coord[2] > min(Y)-0.02*diff(ylim) & coord[2] < max(Y)+0.02*diff(ylim)) {
      action <- ACTION[which(distance == min(distance))]
    } else {
      action <- "none"
    }
    return(action)
  }
}

