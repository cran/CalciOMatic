\name{caMonoExp}
\alias{caMonoExp}
\title{MonoExponential Time Course of Intracellular Calcium
  Concentration}
\description{
  The function \code{caMonoExp} returns a vector of intracellular
  calcium concentration (\eqn{[Ca^{2+}]}{[Ca^2+]}) vs time values \code{t}. A \code{dCa} jump
  occurs at \code{tOn}, followed by a monoexponential return to baseline
  value \code{Ca0}, with time constant \code{tau}.
}
\usage{
caMonoExp(t = 1, tOn = 1, Ca0 = 0.05, dCa = 0.1, tau = 3)
}
\arguments{
  \item{t}{a vector of time values at which \eqn{[Ca^{2+}]}{[Ca^2+]} is computed (in s)}
  \item{tOn}{the time of the \eqn{[Ca^{2+}]}{[Ca^2+]} jump (in s)}
  \item{Ca0}{the baseline \eqn{[Ca^{2+}]}{[Ca^2+]} (in \eqn{\mu}{mu}M)}
  \item{dCa}{the \eqn{[Ca^{2+}]}{[Ca^2+]} jump occuring at \code{tOn} (in \eqn{\mu}{mu}M)}
  \item{tau}{the time constant of the \eqn{[Ca^{2+}]}{[Ca^2+]} monoexponential return to
    baseline (in s)}
}
\value{
  A vector containing the \eqn{[Ca^{2+}]}{[Ca^2+]} values. The vector has the two following
  attributes:
  \item{Time}{a copy of argument \code{t}}
  \item{tOn}{a copy of argument \code{tOn}}
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{caBiExp}}, \code{\link{caMonoBiExpFromIG}}}
\examples{
  ## Simulate a monoexponential calcium transient
  Ca <- caMonoExp(t    = seq(0,10,0.1),
                  tOn  = 2,
                  Ca0  = 0.25,
                  dCa  = 1,
                  tau  = 2)

  ## Plot the calcium transient vs. time
  plot(attr(Ca,"Time"),Ca,type="l")

  ## Add a vertical dashed line at tOn
  abline(v = attr(Ca,"tOn"), lty = 2)
}
