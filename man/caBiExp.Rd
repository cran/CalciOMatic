\name{caBiExp}
\alias{caBiExp}
\title{BiExponential Time Course of Intracellular Calcium Concentration}
\description{
  The function \code{caBiExp} returns a vector of intracellular
  calcium concentration (\eqn{[Ca^{2+}]}{[Ca^2+]}) vs time values \code{t}. A \code{dCa} jump
  occurs at \code{tOn}, followed by a biexponential return to baseline
  value \code{Ca0}. The fast time constant (\code{tau}) has a weight
  \code{fact} (between 0 and 1), the slow time constant
  (\code{tau}+\code{dtau}) has a weight \code{1-fact}.
}
\usage{
caBiExp(t = 1, tOn = 1, Ca0 = 0.05, dCa = 0.1, tau = 3, fact = 1, dtau = 2)
}
\arguments{
  \item{t}{a vector of time values at which \eqn{[Ca^{2+}]}{[Ca^2+]} is computed (in s)}
  \item{tOn}{the time of the \eqn{[Ca^{2+}]}{[Ca^2+]} jump (in s)}
  \item{Ca0}{the baseline \eqn{[Ca^{2+}]}{[Ca^2+]} (in \eqn{\mu}{mu}M)}
  \item{dCa}{the \eqn{[Ca^{2+}]}{[Ca^2+]} jump occuring at \code{tOn} (in \eqn{\mu}{mu}M)}
  \item{tau}{the fast time constant of the \eqn{[Ca^{2+}]}{[Ca^2+]} biexponential return to
    baseline (in s)}
  \item{fact}{the relative weight of the fast time constant of the
    biexponential decay (a real number between 0 and 1). The relative
    weight of the slow time constant is given by \code{1-fact}}
  \item{dtau}{added to \code{tau}, defines the slow time constant of the
    \eqn{[Ca^{2+}]}{[Ca^2+]} biexponential return to baseline (in s)}
}
\value{
  A vector containing the \eqn{[Ca^{2+}]}{[Ca^2+]} values. The vector has the two following
  attributes:
  \item{Time}{a copy of argument \code{t}}
  \item{tOn}{a copy of argument \code{tOn}}
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{caMonoExp}}, \code{\link{caMonoBiExpFromIG}}}
\examples{
  ## Simulate a biexponential calcium transient
  Ca <- caBiExp(t    = seq(0,20,0.1),
                tOn  = 2,
                Ca0  = 0.25,
                dCa  = 1,
                tau  = 2,
                fact = 0.7,
                dtau = 8)

  ## Plot the calcium transient vs. time
  plot(attr(Ca,"Time"), Ca, type="l")

  ## Add a vertical dashed line at tOn
  abline(v = attr(Ca,"tOn"), lty = 2)
}
