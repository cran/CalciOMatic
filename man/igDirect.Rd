\name{igDirect}
\alias{igDirect}
\title{Provide an Initial Guess For the Experiment-Specific Parameters
  of a Direct Fit}
\description{
  The function \code{igDirect} provides an initial guess for the
  experiment-specific parameters of fluorescence transients obtained
  with a ratiometric dye (the background fluorescence \code{log_S_B_340}
  and \code{log_S_B_380}, as well as the scaling coefficient
  \code{log_phi})
}
\usage{
igDirect(adu_B_340, adu_340, adu_B_380, adu_380,
         ig_ratio, t, tOn = 1, subset = 1:length(t),
         R_min = 0.136, R_max = 2.701, K_eff = 3.637, K_d = 0.583,
         B_T = 100, T_340 = 0.015, T_380 = 0.006, P = 400, P_B = 400)
}
\arguments{
  \item{adu_B_340}{the background fluorescence at 340 nm}
  \item{adu_340}{the fluorescence transient at 340 nm}
  \item{adu_B_380}{the background fluorescence at 380 nm}
  \item{adu_380}{the fluorescence transient at 380 nm}
  \item{ig_ratio}{the initial guess list for the parameters of the
  \eqn{[Ca^{2+}]}{Ca^2+} transient, returned by the \code{IG_Ratio}
  function}
  \item{t}{a vector of time values at which the fluorescence values were
  obtained (in s)}
  \item{tOn}{the time of the fluorescence jump (in s)}
  \item{subset}{a vector of time indices to consider (generally the
  whole fluorescence signals)}
  \item{R_min}{the minimum fluorescence ratio between the measurements
    at 340 and 380 nm. This parameter is obtained from calibration
    experiments}
  \item{R_max}{the maximum fluorescence ratio between the measurements
    at 340 and 380 nm. This parameter is obtained from calibration
    experiments}
  \item{K_eff}{the effective dissociation constant of the dye in the
    cell (in \eqn{\mu M}{muM}). This parameter is obtained from
    calibration experiments}
  \item{K_d}{the dissociation constant of the dye in the cell (in
    \eqn{\mu M}{muM}). This parameter is obtained from calibration
    experiments}
  \item{B_T}{the total concentration of the dye in the cell (in \eqn{\mu
      M}{muM})}
  \item{T_340}{the exposure time at 340 nm (in s)}
  \item{T_380}{the exposure time at 380 nm (in s)}
  \item{P}{the number of pixels od the Region Of Interest (ROI)}
  \item{P_B}{the number of pixels od the Beckground Region}
}
\details{
  The intial guesses for \code{log_S_B_340} and \code{log_S_B_380} are
  obtained by averaging the \code{adu_B_340} and \code{adu_B_380}
  signals, divided by \code{T_340*P_B} and \code{T_380*P_B}
  respectively, and by taking their logarithm.

  The initial guess for \code{log_phi} is obtained by applying a linear
  (zero-intercept) regression between the following signals, and by
  taking the logarithm of the slope:
  \deqn{%
    \code{c} \left(\frac{adu_{340}}{T_{340} \cdot P} - S_{B,340},
                   \frac{adu_{380}}{T_{380} \cdot P} - S_{B,380}
		   \right)
  }{%
    c(adu_340/T_340/P - S_B_340, adu_380/T_380/P - S_B_380)
  }

  \deqn{%
    \frac{[B_T]}{K_f+[Ca^{2+}]_{ratio}} \cdot \code{c}
    \left(R_{min} \cdot K_{eff} + R_{max} \cdot [Ca^{2+}]_{ratio},
          K_{eff} + [Ca^{2+}]_{ratio} \right)
  }{
    B_T / (K_d+Ca_ratio) * c(R_min*K_eff+R_max*Ca_ratio,K_eff+Ca_ratio)
  }

  In these formulas, \eqn{[Ca^{2+}]_{ratio}}{Ca_ratio} refers to the
  calcium concentration transient estimated with the initial guess
  parameters listed in the \code{ig_ratio} argument
}
\value{
  A named list of class \code{"initial_guess"}, containing initial
  guesses (IG) for the logarithms of the three experiment-specific
  parameters: The background fluorescences at 340 and 380 nm
  (\code{log_S_B_340} and \code{log_S_B_380} respectively) and the
  amplitude coefficient \code{log_phi}
}
\references{see the \code{\link{fluo}} documentation for details
    about the data generation model}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{igRatio}},
  \code{\link{fluo}}}
\examples{
## Parameters of the monoexponential calcium transient
tOn  <- 1
Time <- seq(0,12,length.out=160)
Ca0  <- 0.10
dCa  <- 0.25
tau  <- 1.5

## Create a monoexponential calcium decay
Ca <- caMonoExp(t = Time, tOn = tOn,
                Ca0 = Ca0, dCa = dCa, tau = tau)

## Calibrated parameters
R_min <- list(value=0.136, mean=0.136, se= 0.00363, USE_se = TRUE)
R_max <- list(value=2.701, mean=2.701, se= 0.151,   USE_se = TRUE)
K_eff <- list(value=3.637, mean=3.637, se= 0.729,   USE_se = TRUE)
K_d   <- list(value=0.583, mean=0.583, se= 0.123,   USE_se = TRUE)

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

## Define Background and Transient fluorescence
## signals at 340 and 380 nm
adu_B_340 <- fluo(Ca=rep(0,nb_B),
                  R_min=R_min$value, R_max=R_max$value,
                  K_eff=K_eff$value, K_d=K_d$value,
                  B_T=0, phi=phi, S_B=S_B_340,
                  T_stim=T_340, P=P, P_B=P_B)

adu_340   <- fluo(Ca=Ca,
                  R_min=R_min$value, R_max=R_max$value,
                  K_eff=K_eff$value, K_d=K_d$value,
                  B_T=B_T, phi=phi, S_B=S_B_340,
                  T_stim=T_340, P=P, P_B=P_B)
  
adu_B_380 <- fluo(Ca=rep(0,nb_B),
                  R_min=1, R_max=1,
                  K_eff=K_eff$value, K_d=K_d$value,
                  B_T=0, phi=phi, S_B=S_B_380,
                  T_stim=T_380, P=P, P_B=P_B)

adu_380   <- fluo(Ca=Ca,
                  R_min=1, R_max=1,
                  K_eff=K_eff$value, K_d=K_d$value,
                  B_T=B_T, phi=phi, S_B=S_B_380,
                  T_stim=T_380, P=P, P_B=P_B)

## Add Poissonian noise to these signals
adu_B_340 <- rpois(length(adu_B_340), adu_B_340)
adu_340   <- rpois(length(adu_340), adu_340)
adu_B_380 <- rpois(length(adu_B_380), adu_B_380)
adu_380   <- rpois(length(adu_380), adu_380)

## Extract the noisy calcium transient
## (from the ratiometric transformation)
Ca_noisy <- caFromRatio(adu_B_340, adu_340,
                        adu_B_380, adu_380,
                        T_340, T_380,
                        P, P_B,
                        R_min, R_max, K_eff,
                        Plot = FALSE)

## Perform a ratiometric fit to determine
## the calcium dynamics parameters
ratio_fit <- ratioFitFromCa(Ca_noisy, t=Time, tOn, type="mono")

## List the fitted parameters and create
## the corresponding calcium transient
ig_mono <- as.list(coefficients(ratio_fit))
class(ig_mono) <- "initial_guess"

## Perform an Initial Guess for the Experiment-Specific Parameters
ig_direct <- igDirect(adu_B_340 = adu_B_340,
                      adu_340 = adu_340,
                      adu_B_380 = adu_B_380,
                      adu_380 = adu_380,
                      ig_ratio = ig_mono,
                      t = Time, tOn = tOn, subset = 1:length(Time),
                      R_min = R_min$value, R_max = R_max$value,
                      K_eff = K_eff$value, K_d = K_d$value,
                      B_T = 100, T_340 = T_340, T_380 = T_380, P = P, P_B = P_B)

## Compare the initial guess with the known values of the parameters
print(exp(as.vector(unlist(ig_direct))))
print(c(Ca0, dCa, tau, phi, S_B_340, S_B_380))
}
