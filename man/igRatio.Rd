\name{igRatio}
\alias{igRatio}
\title{Provide an Initial Guess For a Calcium Concentration Ratiometric
  Fit}
\description{
  The function \code{IG_Ratio} provides an initial guess for the
  parameters of an intracellular calcium concentration transient
  obtained after a ratiometric transformation. The transients considered
  here are either mono- or biexponential. Parameters of a \eqn{F_T}{F_T}
  transient (instead of a \eqn{[Ca^{2+}]}{Ca^2+} transient) can also be
  estimated.
}
\usage{
igRatio(Ca, t, tOn = 1, type = "mono")
}
\arguments{
  \item{Ca}{a vector of \eqn{[Ca^{2+}]}{Ca^2+} values (in muM)}
  \item{t}{a vector of time values at which Ca is computed (in s)}
  \item{tOn}{the time of the Ca jump (in s)}
  \item{type}{a character string (either \code{"mono"} or \code{"bi"}),
    indicating which type of exponential decay to consider}
}
\details{
  This function provides initial guesses for three or five parameters,
  depending on the type of exponential return to baseline. If
  \code{type} is set to \code{"mono"}, three parameters are guessed:
  \itemize{
  \item{\code{log_Ca0}}{is obtained by averaging the signal prior to
    \code{tOn} (and taking its logarithm)}
  \item{\code{log_dCa}}{is obtained by substracting the IG for
    \code{Ca0} from the maximum value of the signal (and taking its
    logarithm)}
  \item{\code{log_tau}}{is obtained from a linear regression on a
    rescaled (by \code{dCa}) and time offset (by \code{Ca0}) version of
    the original signal (and taking its logarithm)} \newline \newline
  }
If \code{type} is set to \code{"bi"}, two more parameters (\code{mu} and
\code{log_dtau}) are guessed. For that purpose, the slow time constant
of the signal is first guessed: successive linear fits of the end part
of the signal log-normalized signal (of increasing lengths, from T to
Tend, with decreasing T) are performed, until the fitted time constant
reaches T. Then, \eqn{\tau_s}{tau_s} is set to T and its relative weight
(from which \code{mu} arises) is deduced from the fit
intercept. Considering the original slope of the signal then leads to
the fast time constant of the biexponential decay (\code{tau}), thus, to
\code{dtau}.
}
\value{
  A named list of class \code{"initial_guess"}, containing initial
  guesses (IG) for the three/five scalar following components of the
  mono- or bi- exponential calcium decay:
  \item{log_Ca0}{IG for the logarithm of the signal baseline (which is
    always positive in the case of a \eqn{[Ca^{2+}]}{Ca^2+} or
    \eqn{[F_T]}{F_T} concentration)}
  \item{log_dCa}{IG for the logarithm of the signal jump, if \code{F_T}
    is set to \code{TRUE} (indeed, the \eqn{[Ca^{2+}]}{Ca^2+}
    concentration jump is always positive)}
  \item{dCa}{IG for the signal jump, if \code{F_T} is set to \code{TRUE}
    (indeed, the \eqn{[F_T]}{F_T} concentration jump can be eiter
    positive or negative)}
  \item{log_tau}{IG for the logarithm of the time constant of the
    monexponential decay (if \code{type} is set to \code{"mono"}) or the
    fast time constant of the biexponential decay (if \code{type} is set
    to \code{"bi"})}
  \item{mu}{IG for the real number (between -\code{Inf} and +\code{Inf})
    defining the relative weight of the fast and slow time constants of
    the biexponential decay (if \code{type} is set to
    \code{"bi"}). The weight of the fast time constant is given by
    \eqn{\frac{exp(\mu)}{1+exp(\mu)}}{exp(mu)/(1+exp(mu))}}
  \item{log_dtau}{IG for the logarithm of the d\eqn{\tau}{tau} defining
    the slow time constant of the biexponential decay (if \code{type} is
    set to \code{"bi"}). This slow time constant is given by
    \eqn{\tau_s=\tau+d\tau}{tau_s=tau+dtau}}
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{igDirect}},
 \code{\link{caMonoBiExpFromIG}}}
\examples{
## Parameters of the monoexponential calcium transient
tOn <- 1
Time <- seq(0,12,length.out=160)
Ca0 <- 0.10
dCa <- 0.25
tau <- 1.5

## Calibration parameters
R_min <- list(value=0.136, mean=0.136, se= 0.00363, USE_se = FALSE)
R_max <- list(value=2.701, mean=2.701, se= 0.151,   USE_se = FALSE)
K_eff <- list(value=3.637, mean=3.637, se= 0.729,   USE_se = FALSE)
K_d   <- list(value=0.583, mean=0.583, se= 0.123,   USE_se = FALSE)

## Experiment-specific parameters
nb_B    <- 1
B_T     <- 100.0
T_340   <- 0.015
T_380   <- 0.006
P       <- 200
P_B     <- 200
phi     <- 2
S_B_340 <- 30
S_B_380 <- 80

## Create a monoexponential calcium decay
Ca <- caMonoExp(t=Time,
                tOn=tOn,
                Ca0=Ca0,
                dCa=dCa,
                tau=tau
                )

## Simulate the corresponding ratiometric experiment
df <- ratioExpSimul(nb_B    = nb_B,
                    Ca      = Ca,
                    R_min   = R_min,
                    R_max   = R_max,
                    K_eff   = K_eff,
                    K_d     = K_d,
                    B_T     = B_T,
                    phi     = phi,
                    S_B_340 = S_B_340,
                    S_B_380 = S_B_380,
                    T_340   = T_340,
                    T_380   = T_380,
                    P       = P,
                    P_B     = P_B,
                    ntransients = 1,
                    G       = 1,
                    s_ro    = 0)

## Get the noisy calcium transient from the data frame
Ca_noisy <- caFromDf(df, numTransient=1, Plot=FALSE)

## Find an Initial Guess for the calcium transient parameters
ig_mono <- igRatio(Ca=Ca_noisy, t=Time, tOn=tOn, type="mono")

## Plot the simulated noisy calcium transient over the original
## calcium transient
plot(Time, Ca_noisy, type="l")
lines(Time, Ca, col="blue")

## Add the calcium transient corresponding to the initial guess
lines(Time, caMonoBiExpFromIG(t=Time, tOn=tOn, ig=ig_mono), lwd=2, col="red")

## Add the corresponding legend
legend("topright",c("Ideal","Noisy","Initial Guess"),
       col=c("blue","black","red"),lwd=c(1,1,2))
}
