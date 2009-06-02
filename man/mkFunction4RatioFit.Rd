\name{mkFunction4RatioFit}
\alias{mkFunction4RatioFit}
\title{Function for Ratiometric Fit}
\description{
  The function \code{mkFunction4RatioFit} returns a function
  prediciting the exponential time course of an intracellular
  concentration. The function, the arguments of which depends on the
  value of \code{type}, can be passed to the \code{nls} formula, for a
  ratiometric fit
}
\usage{
mkFunction4RatioFit(type = "mono")
}
\arguments{
  \item{type}{a character string (either \code{"mono"} or \code{"bi"}),
    specifying which type of exponential decay should be considered}
}
\value{
  A function that has five or seven arguments, depending on the
  value of \code{type}. If \code{type} is set to \code{"mono"}, the five
  arguments are the following: \code{t},\code{tOn}, \code{log_Ca0},
  \code{log_dCa} and \code{log_tau}. If \code{type} is set to
  \code{"bi"}, the output function has two more arguments: \code{mu} and
  \code{log_dtau}. For details about the meaning of each argument, see
  the \code{Ca_MonoExp_fct} or \code{Ca_BiExp_fct} functions
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{caMonoExp}}, \code{\link{caBiExp}},
  \code{\link{ratioFitFromCa}}, \code{\link{ratioFitFromDf}}}
\examples{
## Time parameters
tOn <- 1
Time <- seq(0,30,0.1)
Ca0 <- 0.10
dCa <- 0.25
tau <- 1.5
mu <- 2
dtau <- 10

## Define a calcium biexponential decay with the
## mkCa_MonoBiExp_4_RatioFit function
Ca_biexp_fct <- mkFunction4RatioFit(type = "bi")

Ca_biexp_decay1 <- Ca_biexp_fct(t = Time, tOn = tOn,
                                log_Ca0 = log(Ca0),
                                log_dCa = log(dCa),
                                log_tau = log(tau),
                                mu = mu,
                                log_dtau = log(dtau))

## Define a calcium biexponential decay with the
## caBiExp function
Ca_biexp_decay2 <- caBiExp(t = Time, tOn = tOn, Ca0 = Ca0,
                           dCa = dCa, tau = tau,
                           fact = 1/(1+exp(-mu)), dtau = dtau)

## Check that both decays are similar
plot(Time, Ca_biexp_decay1, type="l", lwd=2)
lines(Time, Ca_biexp_decay2, col="red", lty=2, lwd=2)
}
