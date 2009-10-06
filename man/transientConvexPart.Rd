\name{transientConvexPart}
\alias{transientConvexPart}
\title{Select the Convex or Concave Part of a Transient}
\description{
  The functon \code{transientConvexPart} extracts the indices of a
  given transient where the signal is monotonically convex or concave,
  after a local peak (maximum or minimum) at the beginning of the
  transient
}
\usage{
transientConvexPart(transient, t = 1, tOn = 1)
}
\arguments{
  \item{transient}{the vector to work on}
  \item{t}{a vector of time values at which \code{transient} has been
    obtained (in s)}
  \item{tOn}{the time of the transient jump (in s)}
}
\details{
  The function \code{transientConvexPart} is designed to work on
  transients of the following form: First, prior to \code{tOn}, a
  baseline; Then, at \code{tOn}, a sharp (positive or negative) jump,
  which leads to a global maximum or minimum; Finally, a monotonic
  return to baseline. Real \eqn{[Ca^{2+}]}{Ca^2+} or fluorescence
  transients, on which this function is applied, are generally of this
  form. The function smoothes the input transient, finds the time (after
  the peak) at which the second derivative changes sign, and returns its
  index
}
\value{
  An integer, which is the index of the transient (after the peak) at
  which the second derivative changes sign, and returns its index
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\examples{
## Parameters of the monoexponential calcium transient
tOn <- 1
Time <- seq(0,12,length.out=160)
Ca0 <- 0.10
dCa <- 0.25
tau <- 1.5

## Calibration parameters
R_min <- 0.136
R_max <- 2.701
K_eff <- 3.637
K_d   <- 0.583

## Experiment-specific parameters
nb_B <- 5
B_T <- 100.0
T_340 <- 0.015
T_380 <- 0.006
P <- 200
P_B <- 200
phi <- 20
S_B_340 <- 300
S_B_380 <- 800

## Create a monoexponential calcium decay
Ca_Mono <- caMonoExp(t=Time,
                     tOn=tOn,
                     Ca0=Ca0,
                     dCa=dCa,
                     tau=tau)

## Simulate the corresponding ratiometric experiment
df_Mono <- ratioExpSimul(nb_B = nb_B, Ca = Ca_Mono,
                         R_min = R_min, R_max = R_max,
                         K_eff = K_eff, K_d = K_d,
                         B_T = B_T, phi = phi, P = P, P_B = P_B,
                         ntransients = 1,
                         S_B_340 = S_B_340, S_B_380 = S_B_380,
                         T_340 = T_340, T_380 = T_380, G = 1, s_ro = 0)

## Get the fluorescence transients at 340 and 380 nm, respectively
t <- with(df_Mono,Time[!is.na(Time) & lambda==340])
adu_340 <- with(df_Mono,adu[!is.na(Time) & lambda==340])
adu_380 <- with(df_Mono,adu[!is.na(Time) & lambda==380])

## Calculate the indices of convex/concave starts at both wavelengths
idx_340 <- transientConvexPart(t = t, tOn = tOn, transient = adu_340)
idx_380 <- transientConvexPart(t = t, tOn = tOn, transient = adu_380)

## Plot both transients, with a specific color for the
## portions of interest
layout(matrix(c(1,2),ncol=1))

plot(t[c(1:idx_340)], adu_340[c(1:idx_340)], type="l",
     xlim = c(Time[1],Time[length(Time)]))
lines(t[c(idx_340:length(adu_340))],
      adu_340[c(idx_340:length(adu_340))], col="blue")

plot(t[c(1:idx_380)], adu_380[c(1:idx_380)], type="l",
     xlim = c(Time[1], Time[length(Time)]))
lines(t[c(idx_380:length(adu_380))],
      adu_380[c(idx_380:length(adu_380))], col="red")
}
