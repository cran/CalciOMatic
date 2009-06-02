\name{CalciOMatic-package}
\alias{CalciOMatic-package}
\alias{CalciOMatic}
\docType{package}
\title{Automatic Calcium Imaging Analysis}
\description{Simulate and analyse calcium imaging data obtained with
  ratiometric dyes. The package provides tools to fit parametric models
  of calcium dynamics on experimental data. Two methods are available:
  the classical 'ratiometric' method and a new 'direct' method, which
  does not imply any data ratioing and fits directly the fluorescence
  transients recorded at two excitation wavelengths. The latter method
  allows for the construction of meaningful confidence intervals on the
  calcium dynamics parameters
}
\details{
\tabular{ll}{
  Package: \tab CalciOMatic\cr
  Type: \tab Package\cr
  Version: \tab 1.1-2\cr
  Date: \tab 2009-06-02\cr
  License: \tab GPL (>= 2)\cr
  Depends: \tab cobs\cr
}
}
\author{Sebastien Joucla, Christophe Pouzat
  
Maintainer: Sebastien Joucla <sebastien.joucla@parisdescartes.fr>
}
\references{Joucla S, Pippow A, Kloppenburg P and Pouzat C (2009)
  Quantitative estimation of calcium dynamics from ratiometric
  measurements: A direct, non-ratioing, approach. J Neurophyiol,
  submitted}
\keyword{package}
\examples{
## Load the data set from cockroach olfactory interneurons
data(inVitro)

## Define the calibrated parameters of the calcium indicator (Fura-2)
R_min <- list(value=0.136, mean=0.136, se= 0.00363, USE_se = TRUE)
R_max <- list(value=2.701, mean=2.701, se= 0.151,   USE_se = TRUE)
K_eff <- list(value=3.637, mean=3.637, se= 0.729,   USE_se = TRUE)
K_d   <- list(value=0.583, mean=0.583, se= 0.123,   USE_se = TRUE)

## Characteristics of the CCD camera, obtained from 'calibration' experiments
G    <- 0.146
s_ro <- 16.4

## Create the data.frame containing the physiological data:
physioData <- ratioExpPhysio(dataset="inVitro",
                             expe=2, stim=2,
                             idxOn=10,
                             R_min=R_min, R_max=R_max,
                             K_eff=K_eff, K_d=K_d,
                             G=0.146, s_ro=16.4,
                             alphamethod=TRUE)

## Fit the physiological data with the direct method:
physioDirectFit <- directFit(physioData,
                             transients=2,
                             SQRT=TRUE,
                             type="mono",
                             AfterPeak=14)

## Plot the raw and fitted data, as well as plots of goodness of fit
plot(physioDirectFit, numTransient=2, items=1:6)
}
