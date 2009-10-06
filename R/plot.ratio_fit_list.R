plot.ratio_fit_list <-
function(x,
                                y=NULL,
                                numTransient=1,
                                items=1:5,
                                col="black",
                                col2="darkgray",
                                main="Intracellular calcium transient: Ratiometric fit",
                                xlabs=c("Time (s)",
                                        "Time (s)",
                                        "Lag",
                                        "Theoretical quantiles",
                                        "Residuals"),
                                ylabs=c(expression(paste("[",Ca^{2+phantom()}, "] (",mu,"M)")),
                                        "Residuals",
                                        "Autocorrelation function",
                                        "Sample quantiles",
                                        "Counts"),
                                labs=c(expression(A[1]),
                                       expression(A[2]),
                                       "B","C","D"),
                                ylas=1,
                                ask=FALSE,
                                ...
                                ) {
  ## Function plot.ratio_fit_list
  ## Plot an intracellular calcium transient, with raw and fitted data,
  ## as well as different analyses arising from a ratiometric fit list
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
  ## save: a logical value; set to TRUE to save the plotted figures
  ## file: a character string; defines the name of the file containing the figure(s)
  ## width, height: numbers; give the width and height of the pdf file (when save is set to TRUE)
  ## ...: other parameters that can be set to the plot function
  
  ## -----------------------------------------------
  ## Test the characteristics of the input arguments
  ratio_fit_list <- x
  
  if(!inherits(ratio_fit_list,"ratio_fit_list")) {
    stop("use only with \"ratio_fit_list\" objects")
  }
  
  ## Which transients have been fitted
  transients <- c()
  for(i in 1:length(ratio_fit_list)) {
    if(!is.null(ratio_fit_list[[i]])) {
      transients <- c(transients, i)
    }
  }
  ## transients <- sapply(1:length(ratio_fit_list), function(i) ifelse(is.null(ratio_fit_list[[i]]),NULL,i))
  transients_backup <- transients
  transients <- intersect(transients,numTransient)
 
  ## ----------------------------------------------------------
  ## Plot the figures corresponding to the different transients
  for(k in transients) {
    ## Call the plot.ratio_fit function for single 'calcium_ratio_fit' objects
    plot(x=ratio_fit_list[[k]],
         y=NULL,
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
}

