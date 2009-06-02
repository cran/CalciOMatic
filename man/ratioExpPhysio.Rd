\name{ratioExpPhysio}
\alias{ratioExpPhysio}
\title{Gathers the Results of a Ratiometric Experiment in an Object of
  Class "fluo\_rawdata"}
\description{
  The function \code{ratioExpPhysio} gathers the results of a single
  ratiometric experiment with 1 or more fluorescence transients in a
  data frame of class \code{"fluo_rawdata"}, usable by the following
  functions: \code{ratioFitFromDf}, \code{directFit} and
  \code{plot.fluo_rawdata}
}
\usage{
ratioExpPhysio(dataset = "inVitro", expe = 1, stim = 1, idxOn = 10,
               R_min = 0.136, R_max = 2.701, K_eff = 3.637, K_d = 0.583,
               G = 0.146, s_ro = 16.4, alphamethod = TRUE)
}
\arguments{
  \item{dataset}{a character string. The name of the variable containing
    results of ratiometric experiments. The minimal structure of this
    variable is detailed in \code{\link{inVitro}}}
  \item{expe}{the number of the experiment to consider (field "Exp.." of
    the dataset)}
  \item{stim}{a vector of integers specifying the number of the
    stimulations to consider (field "stim.")}
  \item{idxOn}{the index of the time at which the light is set on}
  \item{R_min}{the minimum fluorescence ratio between the measurements
    at 340 and 380 nm. This parameter is obtained from calibration
    experiments}
  \item{R_max}{the maximum fluorescence ratio between the measurements
    at 340 and 380 nm. This parameter is obtained from calibration
    experiments}
  \item{K_eff}{the effective dissociation constant of the dye in the
    cell (in muM). This parameter is obtained from calibration
    experiments}
  \item{K_d}{the dissociation constant of the dye in the cell (in
    muM). This parameter is obtained from calibration experiments}
  \item{G}{the gain of the CCD camera}
  \item{s_ro}{the standard deviation of the read-out process of the
    camera}
  \item{alphamethod}{a logical value. If set to \code{TRUE}, the
    fluorescence measurements at 360 nm (contained in the dataset) are
    used to estimate the isocoefficient \code{alpha}}
}
\details{
  Details about the estimation of the isocoefficient \code{alpha} with
  the fluorescence measurements at 340, 380 and 360 nm are given in
  Joucla et al. (2009, J Neurophysiol) (see Methods and Appendix C)

  The variable entitled \code{"dataset"}, which contains all experiment
  informations, should be a list with fields named \code{"Exp01"},
  \code{"Exp02"}, etc. Each of this field should be a list with (at
  least) the following fields. The information contained in these fields
  are retrieved and put at the right place in the output data frame:
  \tabular{ll}{
    \code{stim1} \tab (eventually \code{stim2}, \code{stim3}, etc.)\cr
    \code{ } \tab   \cr
    \code{adu340Background} \tab a vector of background fluorescence
      recorded at 340 nm\cr
    \code{ } \tab   \cr
    \code{adu380Background} \tab a vector of background fluorescence
      recorded at 380 nm\cr
    \code{ } \tab   \cr
    \code{P} \tab the number of pixels used for the data binning of the raw
    image, \cr
    \code{ } \tab for fluorescence transient\cr
    \code{ } \tab   \cr
    \code{PBackground} \tab the number of pixels used for the data binning of
      the raw image,\cr
    \code{ } \tab for background fluorescence\cr
    \code{ } \tab   \cr
    \code{furaPipette} \tab the total Fura concentration in the cell (in muM)\cr
    \code{ } \tab   \cr
    \code{exposureTime340} \tab the exposure time at 340 nm (in s)\cr
    \code{ } \tab   \cr
    \code{exposureTime380} \tab the exposure time at 380 nm (in s)
  }
  Each field of \code{"stim1"} should be a list with at least the
  following fields:
  \tabular{ll}{  
    \code{time} \tab the times at which the fluorescence transient was
    acquired \cr
    \code{ } \tab   \cr
    \code{adu340} \tab the fluorescence transient obtained at 340 nm \cr
    \code{ } \tab   \cr
    \code{adu380} \tab the fluorescence transient obtained at 380 nm \cr
  }
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
  \item{tOn}{the time of the light onset (in s)}
  \item{T_stim}{a vector containing the exposure time at 340 nm and 380
    nm}
  \item{R_min}{a copy of arg \code{R_min}}
  \item{R_max}{a copy of arg \code{R_max}}
  \item{K_eff}{a copy of arg \code{K_eff}}
  \item{K_d}{a copy of arg \code{K_d}}
  \item{P}{the number of pixels used for data binning of the fluorescence
    transients}
  \item{P_B}{the number of pixels used for data binning of the
    fluorescence transients}
  \item{B_T}{the total Fura concentration in the cell (in muM)}
  \item{nb_B}{the number of background measurements performed at each
    wavelength}
  \item{alpha}{an estimation of the isocoefficient (only if
  \code{alphamethod} is set to \code{TRUE})}
  \item{G}{a copy of rag \code{G}}
  \item{s_ro}{a copy of rag \code{s_ro}}
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{ratioExpSimul}}}
\examples{
## Load the data from cockroach olfactory interneurons
data(inVitro)

## Calibrated parameters
R_min <- list(value=0.136, mean=0.136, se= 0.00363, USE_se = TRUE)
R_max <- list(value=2.701, mean=2.701, se= 0.151,   USE_se = TRUE)
K_eff <- list(value=3.637, mean=3.637, se= 0.729,   USE_se = TRUE)
K_d   <- list(value=0.583, mean=0.583, se= 0.123,   USE_se = TRUE)

## Create the data frame containing the physiological data
## (experiment #2, stimulation #2)
## G and s_ro are the respectively the gain of the CCD camera
## and the standard deviation of its read-out process
physioData <- ratioExpPhysio(dataset="inVitro",
                             expe=2, stim=2,
                             idxOn=10,
                             R_min=R_min, R_max=R_max,
                             K_eff=K_eff, K_d=K_d,
                             G=0.146, s_ro=16.4,
                             alphamethod=TRUE)

## Plot the raw data
plot(physioData)
}
