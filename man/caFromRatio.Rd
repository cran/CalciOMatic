\name{caFromRatio}
\alias{caFromRatio}
\title{Get Calcium Concentration From Fluorescence Signals, Using the
  Ratiometric Transformation}
\description{
  The function \code{caFromRatio} applies the ratiometric transformation
  to vectors of fluorescence (including background fluorescence) and
  returns the corresponding intracellular calcium concentration.
}
\usage{
caFromRatio(adu_B_340, adu_340,
            adu_B_380, adu_380,
            T_340 = 0.015, T_380 = 0.006,
            P, P_B,
            R_min = 0.136, R_max = 2.701, K_eff = 3.637,
            Plot = FALSE)
}
\arguments{
  \item{adu_B_340}{a vector of background fluorescence values (photon
    counts) recorded at 340 nm}
  \item{adu_340}{a vector of fluorescence values recorded at 340 nm}
  \item{adu_B_380}{a vector of background fluorescence values recorded
    at 380 nm}
  \item{adu_380}{a vector of fluorescence values recorded at 380 nm}
  \item{T_340}{the exposure time at 340 nm (in s)}
  \item{T_380}{the exposure time at 380 nm (in s)}
  \item{P}{the number of pixels od the Region Of Interest (ROI)}
  \item{P_B}{the number of pixels od the Beckground Region}
  \item{R_min}{the minimum fluorescence ratio between the measurements
    at 340 and 380 nm. This parameter is obtained from calibration experiments}
  \item{R_max}{the maximum fluorescence ratio between the measurements
    at 340 and 380 nm. This parameter is obtained from calibration experiments}
  \item{K_eff}{the effective dissociation constant of the dye in the cell (in
    muM). This parameter is obtained from calibration experiments}
  \item{Plot}{a logical value: Set to TRUE to plot the calcium transient
    deduced from the ratiometric transformation}
}
\details{
  The calcium imaging technique makes use of the ability of a
  fluorescent dye (e.g. Fura) to bind with calcium ions present inside
  a neural cell. Briefly, photons emitted by the calcium-free and
  calcium-bound forms of the dye are recorded by a CCD camera, following
  the illumination of the tissue by a laser at relevant wavelengths
  (corresponding to the excitation maxima of the free and bound forms of
  the dye). In the case of a ratiometric dye, an algebraic relationship
  links the intracellular calcium concentration and the photon counts at
  both wavelnegths (340 and 380 nm, in the case of Fura-2). It is thus
  possible to retrieve the intracellular calcium concentration from the
  ratio of the photon counts recorded at these two wavelenghts (after
  substraction of the background fluorescence): This is the ratiometric
  transformation. The ratio \eqn{R}{R} is defined as:
\deqn{%
R = \frac{\frac{1}{P} \cdot adu_{340} - \frac{1}{P_B} \cdot
  adu_{B,340}}{\frac{1}{P} \cdot adu_{380} - \frac{1}{P_B} \cdot adu_{B,380}} \cdot \frac{T_{380}}{T_{340}} = %
\frac{R_{min} \cdot K_{eff} + R_{max} \cdot [Ca^{2+}]}{K_{eff} + [Ca^{2+}]}.
}{%
R = (1/P*adu_340-1/P_B*adu_B_340)/(1/P*adu_380-1/P_B*adu_B_380)*T_380/T_340 = (R_min*K_eff+R_max*[Ca^2+])/(K_eff+[Ca^2+]).
}
Then, the intracellular calcium concentration is given by:
\deqn{%
[Ca^{2+}] = K_{eff} \cdot \frac{R-R_{min}}{R_{max}-R}
}{%
Ca^2+ = K_eff * (R-R_min) / (R_max-R)
}
}
\value{
  A vector of intracellular calcium concentrations calculated with the
  ratiometric transformation described above.
}
\references{the CalciOMatic manuscript submitted to the Journal of Neuropysiology}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\examples{
## (0) "Experimental" parameters

## Parameters of the monoexponential calcium transient
tOn <- 1
Time <- seq(0,10,0.1)
Ca0 <- 0.10
dCa <- 0.25
tau <- 1.5

## Calibration parameters
R_min <- list(value=0.136, mean=0.136, se= 0.00363, USE_se = TRUE)
R_max <- list(value=2.701, mean=2.701, se= 0.151,   USE_se = TRUE)
K_eff <- list(value=3.637, mean=3.637, se= 0.729,   USE_se = TRUE)
K_d   <- list(value=0.583, mean=0.583, se= 0.123,   USE_se = TRUE)

## Experiment-specific parameters
nb_B    <- 5
B_T     <- 100.0
T_340   <- 0.015
T_380   <- 0.006
P       <- 1000
P_B     <- 1000
phi     <- 1.25
S_B_340 <- 100/P/T_340
S_B_380 <- 100/P/T_380

## (1) Create a monoexponential calcium decay
Ca_Mono <- caMonoExp(t = Time, tOn = tOn,
                     Ca0 = Ca0, dCa = dCa, tau = tau)

## (2) Create the background and transient fluorescence signals
adu_B_340 <- rep(fluo(Ca=rep(0,nb_B),
                      R_min=R_min$value, R_max=R_max$value,
                      K_eff=K_eff$value, K_d=K_d$value,
                      B_T=0, phi=phi, S_B=S_B_340, T_stim=T_340, P=P, P_B=P_B))
  
adu_340 <- rep(fluo(Ca=Ca_Mono,
                    R_min=R_min$value, R_max=R_max$value,
                    K_eff=K_eff$value, K_d=K_d$value,
                    B_T=B_T, phi=phi, S_B=S_B_340, T_stim=T_340, P=P, P_B=P_B))

adu_B_380 <- rep(fluo(Ca=rep(0,nb_B),
                      R_min=1, R_max=1, K_eff=K_eff$value, K_d=K_d$value,
                      B_T=0, phi=phi, S_B=S_B_380, T_stim=T_380, P=P, P_B=P_B))
  
adu_380 <- rep(fluo(Ca=Ca_Mono,
                    R_min=1, R_max=1, K_eff=K_eff$value, K_d=K_d$value,
                    B_T=B_T, phi=phi, S_B=S_B_380, T_stim=T_380, P=P, P_B=P_B))

## (3) Get the noisy calcium transient from the ratiometric transformation
Ca <- caFromRatio(adu_B_340, adu_340,
                  adu_B_380, adu_380,
                  T_340 = 0.015, T_380 = 0.006,
                  P, P_B,
                  R_min = R_min, R_max = R_max, K_eff = K_eff,
                  Plot = TRUE)

## (4) Superimpose the original calcium transient
lines(Ca_Mono, lty=2, col="red")
}
