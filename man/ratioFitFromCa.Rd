\name{ratioFitFromCa}
\alias{ratioFitFromCa}
\title{Perform a Fit on an Intracellular Calcium Concentration Vector}
\description{
  The function \code{ratioFitFromCa} performs a fit on an intracellular
  calcium concentration transient. The transient is fitted with a mono-
  or a biexponential decay, depending on the value of \code{type}
}
\usage{
ratioFitFromCa(Ca, t, tOn, type = "mono", ig = NULL,
               Plot = FALSE, Fit = TRUE, AfterPeak = FALSE,
               Trace = FALSE, WarnOnly = TRUE)
}
\arguments{
  \item{Ca}{a vector of calcium concentration (in muM)}
  \item{t}{a vector of latencies at which the calcium concentration was
    obtained (is s)}
  \item{tOn}{the time of the calcium concentration jump (in s)}
  \item{type}{a character string (either \code{"mono"} or \code{"bi"}),
    specifying the type of calcium exponential decay to consider}
  \item{ig}{an object of class \code{"initial_guess"}, giving values of
    the calcium dynamics parameters to initiate the fitting process with
    using \code{nls}. This is a list with the following fields:
    (\code{"log_Ca0"}, \code{"log_dCa"}, \code{"log_tau"} for a
    monoexponential decay, eventually \code{"mu"} and \code{"log_dtau"}
    for a biexponential decay). If \code{ig} is not an object of class
    \code{"initial_guess"}, initial guesses are estimated using the
    \code{"igRatio"} function}
  \item{Plot}{a logical value. Set to \code{TRUE} to plot the original
    signals, the initial guess and the fit results}
  \item{Fit}{a logical value. Set to \code{TRUE} to perform the fit, or
    to \code{FALSE} to compute an initial guess only}
  \item{AfterPeak}{a logical or numerical value. Set to \code{FALSE} to perform the
    fit on the whole fluorescence transients, to \code{TRUE} to consider
    only the part before the fluorescence jump and the convex part after the
    fluorescence peak, for both signals, or to an integer to skip a
    given number of samples after the fluorescence jump}
  \item{Trace}{a logical value. Set to \code{TRUE} to print results for
    successive steps of the optimization algorithm}
  \item{WarnOnly}{a logical value. Set to \code{TRUE} to go on even if
    the fit produced an error}
}
\details{
  This function can be used to fit any signal made with a first part
  corresponding to a baseline signal, then an (almost) instantaneous
  rise, and a mono- or bi-exponential return to baseline. Since the
  logarithm of all parameters (\code{Ca0}, \code{dCa} and \code{tau})
  are fitted, the baseline and jumpmustbe positive values
}
\value{
  An object that inherits from both \code{"nls"} and \code{"ratio_fit"}
  classes. The object has the following attributes:
  \item{"Name"}{a character string telling which type of fit has been performed}
  \item{"Time"}{the whole time vector, which includes
    \code{NA}s/\code{NaN}s for the background fluorescence signals and
    the supplementary calibration parameters, when relevant}
  \item{"RawData"}{the \eqn{[Ca^{2+}]}{Ca^2+} signal deduced from the
  ratiometric transformation. This signal, which is the one passed to
  the \code{nls} formula, has two attributes: \code{"var"} is the vector
  of variances estimated from the error propagation method, and
  \code{"Time"} is the vector of latencies at which fluorescence
  measurements were performed}
  \item{"RawDataFrame"}{a copy of the input data frame}
  \item{"FitFunction"}{the function passed to the \code{nls} formula}
  \item{"Subset"}{the indices of the \code{Time} vector used for the fit}
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{transientConvexPart}}, 
  \code{\link{caFromDf}}, 
  \code{\link{igRatio}},
  \code{\link{ratioFitFromDf}}
}
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
                             alphamethod=FALSE)

## Retrieve the calcium concentration from the data frame
Ca_noisy <- caFromDf(df           = physioData,
                     numTransient = 2,
                     Plot         = FALSE)

## Perform a ratiometric fit
physioRatioFit <- ratioFitFromCa(Ca  = Ca_noisy,
                                 t   = attr(Ca_noisy,"Time"),
                                 tOn = attr(Ca_noisy, "tOn"),
                                 type = "mono",
                                 AfterPeak = 14)
}
