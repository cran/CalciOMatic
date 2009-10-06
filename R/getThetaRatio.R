getThetaRatio <-
function(calcium_ratio_fit,
                          ciMode=c("normalApprox","likelihoodRatio"),
                          ...
                          ) {
  ## ##########################################################################
  ## Estimate 95% CI ot the parameters of an arbitrary model
  ## given a [Ca] transient estimated with the ratiometric method.
  ## -------------------------------------------------------------------------
  ## Arguments:
  ## calcium_ratio_fit: an object of class "ratio_fit_optim"
  ## ciMode: should the normal approximation ("normal")
  ##           or the likelihood ratio ("ratio") be used to obtain the CI?
  ## ...: not used
  
  ## Get information from the "calcium_ratio_fit" object
  type <- attr(calcium_ratio_fit,"type")
  t <- attr(calcium_ratio_fit,"Time")
  tOn <- attr(calcium_ratio_fit,"tOn")
  Ca <- attr(calcium_ratio_fit,"RawData")
  rssFct <- attr(calcium_ratio_fit,"RSSFunction")
  mysubset <- attr(calcium_ratio_fit,"Subset")
  vocChol <- attr(calcium_ratio_fit,"VocChol")
  attr(Ca,"cov") <- NULL
  
  est <- calcium_ratio_fit$par
  se <- try(sqrt(diag(solve(calcium_ratio_fit$hessian))), silent=TRUE)
  
  if (ciMode[1] == "normalApprox") {
    if (!inherits(se,"try-error")) {
      ## result <- exp(sapply(seq(along=est), function(idx) est[idx] +qt(0.975,df=length(t))*se[idx]*c(-1,1)))
      result <- sapply(seq(along=est), function(idx) est[idx] +qt(0.975,df=length(t))*se[idx]*c(-1,1))
    } else {
      result <- matrix(NA,nr=2,nc=length(est))
    } ## End of conditional on exists("se")
  } else if(ciMode[1] == "likelihoodRatio") {
    nbPar <- length(est)
    result <- matrix(0,nr=2,nc=nbPar)
    pSeq <- seq(along=est)
    for (pIdx in pSeq) {
      par.m <- est[pIdx]
      profFct <- function(par,x) {
        p <- numeric(nbPar)
        p[pIdx] <- x
        p[pSeq[-pIdx]] <- par
        rssFct(p)
      }
      profRSS <- function(x) optim(par = est[-pIdx],
                                   fn = profFct,
                                   method = "BFGS",
                                   x = x)$value - calcium_ratio_fit$value
      thres <- qchisq(0.95,df=1)
      tFct <- function(x) profRSS(x) - thres      
      par.l <- par.m - 0.5*se[pIdx]
      theSign <- tFct(par.l)*tFct(par.m)
      while(theSign >= 0) {
        par.l <- par.l - se[pIdx]
        theSign <- tFct(par.l)*tFct(par.m)
      }
      par.u <- par.m + 0.5*se[pIdx]
      theSign <- tFct(par.u)*tFct(par.m)
      while(theSign >= 0) {
        par.u <- par.u + se[pIdx]
        theSign <- tFct(par.u)*tFct(par.m)
      }
      ##browser()
      lwr <- uniroot(tFct,c(par.l,par.m))$root
      upr <- uniroot(tFct,c(par.m,par.u))$root
      result[1,pIdx] <- lwr
      result[2,pIdx] <- upr
    } ## End of loop on pIdx
    ## result <- exp(result)
  } ## End of conditional on ciMode[1] == "normalApprox"
  return(result)
}

