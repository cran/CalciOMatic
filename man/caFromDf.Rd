\name{caFromDf}
\alias{caFromDf}
\title{Get Calcium Concentration From a Fluorescence Data Frame, Using
  the Ratiometric Transformation}
\description{
  The function \code{caFromDf} applies the ratiometric transformation to
  vectors of fluorescence (including background fluorescence) contained
  in a data frame and returns the corresponding intracellular calcium
  concentration. The structure of the data frame is defined as in the
  \code{ratioExpsimul} function.
}
\usage{caFromDf(df, numTransient = 1, Plot = FALSE)}
\arguments{
  \item{df}{a data frame of class \code{"fluo_rawdata"} containing all
    relevant information (fluorescence transients, background
    fluorescence, calibration parameters and exposure times). The
    structure of the input data frame must be defined as in the
    \code{ratioExpsimul} function.}
  \item{numTransient}{an integer: The index of the transient to analyse
    in the input data frame \code{df}.}
  \item{Plot}{a logical value: Set to TRUE to plot the calcium transient
    deduced from the ratiometric transformation.}
}
\details{see the help of the \code{\link{caFromRatio}} function.}
\value{
  A vector of intracellular calcium concentrations calculated with the
  ratiometric transformation.
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{ratioExpSimul}}, \code{\link{caFromRatio}}}
\examples{
## (0) 'Experimental' parameters

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

## (2) Simulate the corresponding ratiometric experiment
df_Mono <- ratioExpSimul(nb_B    = nb_B,
                         Ca      = Ca_Mono,
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

## (3) Get the noisy calcium transient from the data frame
Ca_noisy <- caFromDf(df           = df_Mono,
                     numTransient = 1,
                     Plot         = FALSE)


## (4) Plot the simulated noisy calcium transient
##     over the ideal calcium transient
## plot(attr(Ca_noisy,"Time"), Ca_noisy, type = "l", col = "blue")
## lines(Time, Ca_Mono, col="red", lwd = 2)
## abline(v = attr(Ca_noisy,"tOn"), lty = 2)
}
