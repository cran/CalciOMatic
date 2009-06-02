\name{mkFunction4DirectFit}
\alias{mkFunction4DirectFit}
\title{Function for Direct Fit}
\description{
  The function \code{mkFunction4DirectFit} returns a function
  prediciting the exponential time course of fluorescence transients at
  two wavelengths (including background fluorescence). The function, the
  arguments of which depends on the value of \code{type}, can be passed
  to the \code{nls} formula, for a direct fit
}
\usage{
mkFunction4DirectFit(type = "mono", nb_B = 5, transients = 1,
                     alphamethod = TRUE, SQRT = TRUE)
}
\arguments{
  \item{type}{a character string (either \code{"mono"} or \code{"bi"}),
    specifying which type of exponential decay should be considered}
  \item{nb_B}{the number of background measurements performed at each
    wavelength}
  \item{transients}{a vector of integers, specifying the numbers of the
    transients to fit (determines the names of the calcium dynamics
    parameters (\code{log_Ca0_1}, \code{log_Ca0_2}, ...))}
  \item{alphamethod}{a logical value. Set to \code{TRUE} (respectively
    \code{FALSE}) to include \code{alpha} (respectively \code{B_T}) in
    the output function formals}
  \item{SQRT}{a logical value. Set to TRUE (the default value) to return
    the square root of the fluorescence signals}
}
\value{
  A function, the arguments of which depend on the value of \code{type},
  \code{alphamethod} and \code{transients}. The list below describes all
  possible arguments:
  \item{t}{the times at which the fluorescence values are expected (in s)}
  \item{tOn}{the time of the fluorescence jump (in s)}
  \item{adu_B_340}{the background fluorescence measurement(s) at 340 nm}
  \item{adu_340}{the fluorescence transient(s) at 340 nm}
  \item{adu_B_380}{the background fluorescence measurement(s) at 380 nm}
  \item{adu_380}{the fluorescence transient(s) at 380 nm}
  \item{T_340}{the exposure time at 340 nm}
  \item{T_380}{the exposure time at 380 nm}
  \item{P}{the number of pixels of the ROI}
  \item{P_B}{the number of pixels of the background region}
  \item{log_Ca0_1}{the logarithm of the \eqn{[Ca^{2+}]}{Ca^2+} baseline}
  \item{log_dCa_1}{the logarithm of the \eqn{[Ca^{2+}]}{Ca^2+} jump}
  \item{log_tau_1}{the logarithm of the \eqn{[Ca^{2+}]}{Ca^2+} time constant}
  \item{mu_1}{the real number (between -\code{Inf} and +\code{Inf})
    defining the relative weight of the fast and slow time constants of
    the \eqn{[Ca^{2+}]}{Ca^2+} biexponential decay (if \code{type} is set to
    \code{"bi"}). The weight of the fast time constant is given by
    \eqn{\frac{exp(\mu)}{1+exp(\mu)}}{exp(mu)/(1+exp(mu))}}
  \item{log_dtau}{the logarithm of the d\eqn{\tau}{tau} defining
    the slow time constant of the \eqn{[Ca^{2+}]}{Ca^2+} biexponential
    decay (if \code{type} is set to \code{TRUE"bi"}). This slow time
    constant is given by \eqn{\tau_s=\tau+d\tau}{tau_s=tau+dtau}}
  \item{log_phi}{the logarithm of the experiment-specific amplitude
    coefficient}
  \item{log_S_B_340}{the logarithm of the background fluorescence at 340
    nm}
  \item{log_S_B_380}{the logarithm of the background fluorescence at 380
    nm}
  \item{log_R_min}{the logarithm of the minimum ratiometric measurement
    observable}
  \item{log_R_max}{the logarithm of the maximum ratiometric measurement
    observable}
  \item{log_K_eff}{the logarithm of the effective Fura dissociation
    constant in the cell (in muM)}
  \item{log_K_d}{the logarithm of the Fura dissociation constant (in muM)}
  \item{alpha}{the isocoefficient. Present only if \code{alphamethod} is
    set to \code{TRUE}. In this case, \code{B_T} is replaced with
    \eqn{adu_{340}/(T_{340} \cdot P) - adu_{B,340}/(T_{340} \cdot P_B) +
      alpha \cdot (adu_{380}/(T_{380} \cdot P) - adu_{B,380}/(T_{380}
      \cdot P_B))}{adu_340/(T_340*P) - adu_B_340/(T_340*P_B) +
      alpha * (adu_380/(T_380*P) - adu_B_380/(T_380*P_B))} in the
    fluorescence model}
  \item{B_T}{the total concentration of the dye in the cell. Present
    only if \code{alphamethod} is set to \code{FALSE}}
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{mkFluo4DirectFit}}, \code{\link{directFit}}}
\examples{
## Parameters of the biexponential calcium transient
tOn  <- 1
Time <- seq(0,12,length.out=160)
Ca0  <- 0.10
dCa  <- 0.25
tau  <- 1.5
mu <- 0
dtau <- 10

## Calibrated parameters
R_min <- 0.136
R_max <- 2.701
K_eff <- 3.637
K_d   <- 0.583

## Experiment-specific parameters
B_T     <- 100.0
T_340   <- 0.015
T_380   <- 0.006
P       <- 200
P_B     <- 200
phi     <- 2
S_B_340 <- 30
S_B_380 <- 80

## Define a function for fluorescence transients based on
## a monoexponential calcium concentration time course
Fluo_bi_fct <- mkFunction4DirectFit(type = "bi",
                                    nb_B = 5,
                                    transients = 1,
                                    alphamethod = FALSE,
                                    SQRT = TRUE)

## Create the fluorescence transients
Time <- matrix(Time, nrow=1, dimnames=list("1",NULL))
names(tOn) <- "1"
Fluo_bi <- Fluo_bi_fct(t = Time, tOn = tOn,
                       T_340 = T_340, T_380 = T_380,
                       P = P, P_B = P_B,
                       log_Ca0_1 = log(Ca0),
                       log_dCa_1 = log(dCa),
                       log_tau_1 = log(tau),
                       mu_1 = mu,
                       log_dtau_1 = log(dtau),
                       log_phi = log(phi),
                       log_S_B_340 = log(S_B_340),
                       log_S_B_380 = log(S_B_380),
                       log_R_min = log(R_min), log_R_max = log(R_max),
                       log_K_eff = log(K_eff), log_K_d = log(K_d),
                       B_T = B_T)
}
