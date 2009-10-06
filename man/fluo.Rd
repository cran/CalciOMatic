\name{fluo}
\alias{fluo}
\title{Convert Intracellular Calcium Concentration into Fluorescence
  Values}
\description{
  The function \code{fluo} converts an intracellular calcium
  concentration to a photon count, depending on the values of the
  calibration parameters (\code{R_min}, \code{R_max}, \code{K_eff} and
  \code{K_d}) and the experiment-specific parameters (\code{B_T},
  \code{phi}, \code{S_B}, \code{T_stim}, \code{P} and \code{P_B})
}
\usage{
fluo(Ca = 1, R_min = 0.136, R_max = 2.701, K_eff = 3.637, K_d = 0.583,
     B_T = 100, phi = 1.25, S_B = 10, T_stim = 0.015, P = 400, P_B = 400)
}
\arguments{
  \item{Ca}{the intracellular calcium concentration (in \eqn{\mu}{mu}M)}
  \item{R_min}{the minimum fluorescence ratio between the measurements
    at 340 and 380 nm. This parameter is obtained from calibration
    experiments}
  \item{R_max}{the maximum fluorescence ratio between the measurements
    at 340 and 380 nm. This parameter is obtained from calibration
    experiments}
  \item{K_eff}{the effective dissociation constant of the dye in the cell (in
    \eqn{\mu}{mu}M). This parameter is obtained from calibration experiments}
  \item{K_d}{the dissociation constant of the dye in the cell (in
    \eqn{\mu}{mu}M). This parameter is obtained from calibration experiments}
  \item{B_T}{the total dye concentration in the cell (in \eqn{\mu}{mu}M)}
  \item{phi}{a dimensionless scaling experiment-specific parameter}
  \item{S_B}{the background (+ dark current) fluorescence intensity (in
    count/pixel/sec)}
  \item{T_stim}{the exposure time (in s)}
  \item{P}{the number of pixels of the Region Of Interest (ROI)}
  \item{P_B}{the number of pixels of the Background Region}
}
\details{
  The calcium imaging technique makes use of the ability of a
  fluorescent dye (e.g. Fura) to bind with calcium ions presents inside
  a neural cell. Briefly, the photons emitted by a neural tissue are
  recorded by a CCD camera, following the illumination of the tissue at
  a relevant wavelength (corresponding to the excitation properties of
  the free and/or bound dye). The amount of photons emitted depends on
  the intracellular calcium concentration (with which the dye is bound),
  and, in the case of a ratiometric dye, a algebraic relationship links
  both variables. The latter is given by:
\deqn{%
F_{340} = \left(\frac{B_T \cdot \phi}{K_d+[Ca^{2+}]} \cdot \left(R_{min} \cdot K_{eff}+R_{max} \cdot [Ca^{2+}]\right)+S_{B,340}\right) \cdot T_{stim,340} \cdot P
}{%
F_340 = (B_T*phi/(K_d+[Ca^{2+}]) * (R_min*K_eff+R_max*[Ca^{2+}]) + S_B,340) * T_stim,340 * P
}
\deqn{%
F_{380} = \left(\frac{B_T \cdot \phi}{K_d+[Ca^2+]} \cdot \left(K_{eff}+[Ca^2+]\right)+S_{B,380}\right) \cdot T_{stim,380} \cdot P
}{%
F_380 = (B_T*phi/(K_d+[Ca^2+]) * (K_eff+[Ca^2+]) + S_B,380) * T_stim,380 * P
}
The function \code{fluo} determines photon counts according to one of
these two equations, depending on the values of \code{R_min},
\code{R_max} and \code{B_T}
}
\value{
  An object of class \code{"fluo_transient"}, which is a vector
  containing the fluorescence values calculated as described above. The
  object has several attributes, which are:
  \item{Ca}{a copy of argument \code{Ca}}
  \item{R_min}{a copy of argument \code{R_min}}
  \item{R_max}{a copy of argument \code{R_max}}
  \item{K_eff}{a copy of argument \code{K_eff}}
  \item{K_d}{a copy of argument \code{K_d}}
  \item{B_T}{a copy of argument \code{B_T}}
  \item{T_stim}{a copy of argument \code{T_stim}}
  \item{P}{a copy of argument \code{P}}
  \item{S_B}{a copy of argument \code{S_B}}
  \item{phi}{a copy of argument \code{phi}}
}
\references{Joucla S, Pippow A, Kloppenburg P and Pouzat C (2009)
  Quantitative estimation of calcium dynamics from ratiometric
  measurements: a direct, non-ratioing, method, Journal of
  Neurophysiology, in revision}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{caBiExp}},
  \code{\link{caMonoExp}}
}
\examples{
## Parameters of the monoexponential calcium transient
tOn  <- 1
Time <- seq(0,12,length.out=160)
Ca0  <- 0.10
dCa  <- 0.25
tau  <- 1.5

## Calibrated parameters
R_min <- 0.136
R_max <- 2.701
K_eff <- 3.637
K_d <- 0.583

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

## Define Background and Signal fluorescences at 340 and 380 nm
B_340 <- fluo(Ca=rep(0,nb_B), R_min=R_min, R_max=R_max, K_eff=K_eff, K_d=K_d,
              B_T=0, phi=phi, S_B=S_B_340, T_stim=T_340, P=P, P_B=P_B)

F_340 <- fluo(Ca=Ca, R_min=R_min, R_max=R_max, K_eff=K_eff, K_d=K_d,
              B_T=B_T, phi=phi, S_B=S_B_340, T_stim=T_340, P=P, P_B=P_B)
  
B_380 <- fluo(Ca=rep(0,nb_B), R_min=1, R_max=1, K_eff=K_eff, K_d=K_d,
              B_T=0, phi=phi, S_B=S_B_380, T_stim=T_380, P=P, P_B=P_B)

F_380 <- fluo(Ca=Ca, R_min=1, R_max=1, K_eff=K_eff, K_d=K_d,
              B_T=B_T, phi=phi, S_B=S_B_380, T_stim=T_380, P=P, P_B=P_B)

## Plot the fluorescence transients at 340 and 380 nm
layout(matrix(1:2,nrow=2))
plot(Time, F_340, type="l", bty="n")
plot(Time, F_380, type="l", bty="n")
}
