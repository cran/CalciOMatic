\name{caMonoBiExpFromIG}
\alias{caMonoBiExpFromIG}
\title{Mono- or Bi- Exponential Time Course of Intracellular Calcium
  Concentration}
\description{
  The \code{caMonoBiExpFromIG} function returns a vector of
  intracellular calcium concentration (Ca). A \code{dCa} jump occurs at
  \code{tOn}, followed by either a monoexponential or a biexponential
  return to the baseline value \code{Ca0}, depending on the fields of
  the input list \code{ig}. The \code{caMonoBiExpFromIG} function is a
  low-level function of the \code{ratioFitFromCa} and \code{directFit}
  functions
}
\usage{caMonoBiExpFromIG(t = 1, tOn = 1, ig = NULL)}
\arguments{
  \item{t}{a vector of time values at which Ca is computed (in s)}
  \item{tOn}{the time of the Ca jump (in s)}
  \item{ig}{an object of class \code{"initial_guess"}, giving the
    parameters of the decay. This object is a list with the follwing
    numerical fields: \code{"log_Ca0"}, \code{"log_dCa"},
    \code{"log_tau"} in the case of a monoexponential decay. In the case
    of a biexponential decay, the two following fields are also
    included: \code{"mu"} and \code{"log_dtau"}}
}
\value{
  A vector containing the Ca values. The vector has the two following
  attributes:
  \item{Time}{a copy of argument \code{t}}
  \item{tOn}{a copy of argument \code{tOn}}
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{caMonoExp}}, \code{\link{caBiExp}}}
\examples{
  ## Parameters of the  calcium transient
  Ca0 <- 0.10
  dCa <- 0.25
  tau <- 1

  ## List of initial guesses
  ig_mono <- list(log_Ca0 = log(Ca0),
                  log_dCa = log(dCa),
                  log_tau = log(tau))
  class(ig_mono) <- "initial_guess"

  ## Build the calcium transient
  Ca <- caMonoBiExpFromIG(t   = seq(0,10,0.1),
                          tOn = 2,
                          ig  = ig_mono)

  ## Plot the calcium transient vs. time
  plot(attr(Ca,"Time"), Ca, type="l")

  ## Add a vertical dashed line at tOn
  abline(v = attr(Ca,"tOn"), lty = 2)
}
