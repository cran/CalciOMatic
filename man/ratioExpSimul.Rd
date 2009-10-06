\name{ratioExpSimul}
\alias{ratioExpSimul}
\title{Simulate Ratiometric Experiment}
\description{
  The function \code{ratioExpSimul} simulates the results of one
  ratiometric experiment, i.e, the photon counts obtained at both
  wavelengths (340 and 380 nm), knowing the time course of the
  intracellular calcium concentration. The photon counts are described
  as the realization of a Poissonian process.
}
\usage{
ratioExpSimul(nb_B = 5, Ca,
              R_min = 0.136, R_max = 2.701, K_eff = 3.637, K_d = 0.583,
              B_T = 100, phi = 1.25, S_B_340 = 10, S_B_380 = 10,
              T_340 = 0.015, T_380 = 0.006, P = 400, P_B = 400,
              ntransients = 1, G = 1, s_ro = 0, noise = TRUE)
}
\arguments{
  \item{nb_B}{the number of background measurements to simulate before
    the fluorescence transients}
  \item{Ca}{the ideal calcium transient from which fluorescence signales
    arise}
  \item{R_min}{the minimum fluorescence ratio between the measurements
    at 340 and 380 nm. This parameter is obtained from calibration
    experiments}
  \item{R_max}{the maximum fluorescence ratio between the measurements
    at 340 and 380 nm. This parameter is obtained from calibration
    experiments}
  \item{K_eff}{the effective dissociation constant of the dye in the
    cell (in \eqn{\mu}{mu}M). This parameter is obtained from
    calibration experiments}
  \item{K_d}{the dissociation constant of the dye in the cell (in
    \eqn{\mu}{mu}M). This parameter is obtained from calibration
    experiments}
  \item{B_T}{the total concentration of the dye in the cell (in
    \eqn{\mu}{mu}M)}
  \item{phi}{the scaling experiment-specific parameters}
  \item{S_B_340}{the background fluorescence at 340 nm}
  \item{S_B_380}{the background fluorescence at 380 nm}
  \item{T_340}{the exposure time at 340 nm}
  \item{T_380}{the exposure time at 380 nm}
  \item{P}{the number of pixels of the ROI}
  \item{P_B}{the number of pixels of the background region}
  \item{ntransients}{a vector of integers (above or equal to 1)
    specifying the indices of the transients to simulate}
  \item{G}{the gain of the CCD camera}
  \item{s_ro}{the standard deviation of the read-out process of the
    camera}
  \item{noise}{a logical. see details below}
}
\details{
  The way fluorescence values arise from intracellular calcium
  concentration values is described in the \code{\link{fluo}}
  function. Recording fluorescence with a CCD camera noises the photon
  counts, which can be described as the realization of a Poissonian
  process, the parameter of which is the fluorescence value
  itself. Ratiometric experiments are simulated thus simulated by drawing
  Poissonian samples from ideal fluorescence transients. These noisy
  data are then multiplied bu the gain \code{G} of the CCD camera, and
  the standard deviation of the read-out noise (\code{s_ro}) is finally
  added. In case \code{noise} is set to \code{FALSE}, fluorescence data
  arising from the fluorescence model are only digitized (rounded
  towards the nearest integer), but not drawn from a Poisson
  distribution.
}
\value{
  An object of class \code{"fluo_rawdata"}, which is a data frame with
  four columns:
  \tabular{ll}{
    \code{adu} \tab the photon counts (or Analog-to-Digital Units) at both
    wavelengths, \cr
    \code{ } \tab including background fluorescence \cr
    \code{ } \tab   \cr
    \code{Time} \tab the times at which each value in \code{adu} was
    recorded. \cr
    \code{ } \tab For the background fluorescence, \code{Time} is set to NA \cr
    \code{ } \tab   \cr
    \code{lambda} \tab the wavelength at which each value in \code{adu} was
    recorded (a factor) \cr
    \code{transient} \tab the number of the fluorescence transient in
    the input data (can be 1, 2 or 3 \cr
    \tab for transient signals, and 0 for
  background measurements) \cr
  }  
  Data appear in this order : (1) the background fluorescence at 340 nm,
  (2) the fluorescence transient(s) at 340 nm, (3) the background
  fluorescence at 380 nm, (4) the fluorescence transient(s) at 380 nm.
  The object has also the following attributes:
  \item{tOn}{the time at which the stimulation is applied (in s)}
  \item{T_stim}{a vector containing the exposure time at 340 nm and 380
    nm}
  \item{R_min}{a copy of arg \code{R_min}}
  \item{R_max}{a copy of arg \code{R_max}}
  \item{K_eff}{a copy of arg \code{K_eff}}
  \item{K_d}{a copy of arg \code{K_d}}
  \item{P}{a copy of arg \code{P}}
  \item{P_B}{a copy of arg \code{P_B}}
  \item{B_T}{a copy of arg \code{B_T}}
  \item{nb_B}{a copy of arg \code{nb_B}}
  \item{G}{a copy of arg \code{G}}
  \item{s_ro}{a copy of arg \code{s_ro}}
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{ratioExpPhysio}}}
\examples{
## Parameters of the monoexponential calcium transient
tOn  <- 1
Time <- seq(0,12,length.out=160)
Ca0  <- 0.10
dCa  <- 0.25
tau  <- 1.5

## Calibrated parameters
R_min <- list(value=0.136, mean=0.136, se=0.00363, USE_se=FALSE)
R_max <- list(value=2.701, mean=2.701, se=0.151,   USE_se=FALSE)
K_eff <- list(value=3.637, mean=3.637, se=0.729,   USE_se=FALSE)
K_d   <- list(value=0.583, mean=0.583, se=0.123,   USE_se=FALSE)

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
Ca <- caMonoExp(t = Time, tOn = tOn,
                Ca0 = Ca0, dCa = dCa, tau = tau)

## Simulate the corresponding ratiometric experiment
simulData <- ratioExpSimul(nb_B    = nb_B,
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

## Plot the raw data
plot(simulData, numTransient=1)
}