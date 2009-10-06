ratioFitFromDf <- function(df,
                           transients=1,
                           type="mono",
                           ig=NULL,
                           Plot=FALSE,
                           Fit=TRUE,
                           AfterPeak=FALSE
                           ) {
  ## Function ratioFitFromDf
  ## 
  ## Fit one or several of the fluorescence signals contained
  ## in the data frame 'df', with the ratiometric method.
  ## The resulting calcium signal is fitted with a mono- or
  ## bi- exponential decay, depending on the value of 'type'.
  ## 
  ## df: a data.frame of class "fluo_rawdata" containing all informative elements
  ##     of the experiment. Its structure is defined in the function RatioExp
  ## transient: the number of the transient(s) to consider
  ## type: a character string indicating the type of calcium decay to fit
  ##       ("mono"exponential or "bi"exponential)
  ## Plot: a logical value; set to TRUE to plot the raw data,
  ##                        the initial guess and the fitted signals
  ## Fit: a logical value; set to TRUE to perform the fit;
  ##                       set to FALSE to estimate an initial guess only
  ## AfterPeak: a logical/numerical value;
  ##                * set to TRUE to fit only the convex part of the decay
  ##                * set to FALSE to fit the whole signal
  ##                * if set to a positive numerical value, skip Trace time samples
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
  
  ## Test if the data frame df is of class "fluo_rawdata"
  if(class(df)!="fluo_rawdata") {
    stop("The first argument of the Ratio_Fit_fct function must be of class 'fluo_rawdata'")
  }
  
  ## Retrieve the important components of the input data frame
  ## adu <- df$adu
  ## lambda_vec <- unique(df$lambda)
  df_transients <- unique(df$transient)
  if(missing(transients)) transients <- df_transients
  transients <- sort(intersect(df_transients,transients))
  nb_B <- attr(df,"nb_B")
  
  calcium_ratio_fit <- vector("list",max(transients))
  
  for(k in transients) {
    if(length(which(transients==k)) == 1) {
      ## Get the calcium transient and associated time
      Ca <- caFromDf(df=df, numTransient=k, Plot=FALSE)
      t <- attr(Ca,"Time")
      tOn <- attr(Ca,"tOn")
      
      ## Fit it
      calcium_ratio_fit[[k]] <- ratioFitFromCa(Ca,
                                               t,
                                               tOn,
                                               type=type,
                                               ig=ig,
                                               Plot=Plot,
                                               Fit=Fit,
                                               AfterPeak=AfterPeak
                                               )
      
      attr(calcium_ratio_fit[[k]],"RawDataFrame") <- df
    }
  }
  
  class(calcium_ratio_fit) <- "ratio_fit_list"
  
  ## Special case: length(df_transients) == 1
  if(length(transients) == 1) {
    calcium_ratio_fit <- calcium_ratio_fit[[transients]]
  }
  
  return(calcium_ratio_fit)
}
