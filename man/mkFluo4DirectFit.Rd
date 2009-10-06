\name{mkFluo4DirectFit}
\alias{mkFluo4DirectFit}
\title{Define Predicted Fluorescence Signals}
\description{
  The function \code{mkFluo4DirectFit} defines a single signal made of
  four fluorescence signals (the background and fluorescence transients
  at both wavelengths), for use with within the
  \code{mkFunction4DirectFit} function
}
\usage{
mkFluo4DirectFit(Ca, phi, S_B_340, S_B_380, nb_B,
                 R_min, R_max, K_eff, K_d, B_T,
                 T_340, T_380, P, P_B, SQRT = TRUE)
}
\arguments{
  \item{Ca}{the time course of the intracellular calcium concentration}
  \item{phi}{the scaling experiment-specific parameter}
  \item{S_B_340}{the background fluorescence at 340 nm}
  \item{S_B_380}{the background fluorescence at 380 nm}
  \item{nb_B}{the number of background measurements}
  \item{R_min}{the minimum fluorescence ratio between the measurements
    at 340 and 380 nm. This parameter is obtained from calibration experiments}
  \item{R_max}{the maximum fluorescence ratio between the measurements
    at 340 and 380 nm. This parameter is obtained from calibration experiments}
  \item{K_eff}{the effective dissociation constant of the dye in the cell (in
    \eqn{\mu}{mu}M). This parameter is obtained from calibration experiments}
  \item{K_d}{the dissociation constant of the dye in the cell (in
    \eqn{\mu}{mu}M). This parameter is obtained from calibration experiments}
  \item{B_T}{the total concentration of the dye inside the cell (in \eqn{\mu}{mu}M)}
  \item{T_340}{the exposure time at 340 nm}
  \item{T_380}{the exposure time at 380 nm}
  \item{P}{the number of pixels of the ROI}
  \item{P_B}{the number of pixels of the background region}
  \item{SQRT}{a logical value. Set to \code{TRUE} to apply the square
    root transformation to the fluorescence signals}
}
\value{
  a vector containing, in this order: the background fluorescence (in
  count) at 340 nm, the fluorescence transient at 340 nm, the background
  fluorescence at 380 nm and the fluorescence transient at 380 nm. If the
  \code{SQRT} argument is set to \code{TRUE}, the square root of the
  whole signal is returned
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{mkFunction4DirectFit}},
  \code{\link{directFit}}}
\examples{
## Parameters of the biexponential calcium transient
tOn <- 1
Time <- seq(0,12,length.out=160)
Ca0 <- 0.10
dCa <- 0.25
tau <- 1.5
mu <- 0
dtau <- 10

## Calibration parameters
R_min <- 0.136
R_max <- 2.701
K_eff <- 3.637
K_d   <- 0.583

## Experiment-specific parameters
nb_B <- 1
B_T <- 100.0
T_340 <- 0.015
T_380 <- 0.006
P <- 200
P_B <- 200
phi <- 2
S_B_340 <- 30
S_B_380 <- 80

## Create a biexponential calcium decay
Ca_Bi <- caBiExp(t=Time,
                 tOn=tOn,
                 Ca0=Ca0,
                 dCa=dCa,
                 tau=tau,
                 fact=1/(1+exp(-mu)),
                 dtau=dtau)

## Define the whole original fluorescence vector
Fluo_bi <- mkFluo4DirectFit(Ca = Ca_Bi, phi,
                            S_B_340, S_B_380, nb_B,
                            R_min, R_max, K_eff, K_d,
                            B_T, T_340, T_380, P, P_B,
                            SQRT = TRUE)
}
