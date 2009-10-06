\name{directFit}
\alias{directFit}
\title{Perform a Direct Fit of Fluorescence Signals Obtained with a
  Ratiometric Dye}
\description{
  The function \code{directFit} performs a direct fit on fluorescence
  signals obtained with a ratiometric dye. The calcium dynamics are
  fitted with either a mono- or a biexponential decay, depending on the
  value of \code{type}.
}
\usage{
directFit(df, transients = c(1, 2, 3), SQRT = TRUE, ratio = NULL,
          type = "mono", Plot = FALSE, Fit = TRUE,
          AfterPeak = FALSE, Trace = FALSE, WarnOnly = TRUE)
}
\arguments{
  \item{df}{a data frame of class \code{"fluo_rawdata"} containing all
    relevant information (fluorescence transients, background
    fluorescence, calibration parameters and exposure times). The
    structure of the input data frame must be the same as the one
    defined in \code{\link{ratioExpSimul}}}
  \item{transients}{a vector of integers giving the numbers of the
    transients to fit simultaneously}
  \item{SQRT}{a logical value. Set to TRUE (default) to fit the square
    root of the fluorescence signals}
  \item{ratio}{an object of class \code{"initial_guess"} or
    \code{"calcium_fit_ratio"}. If \code{ratio} is a list of class
    \code{"initial_guess"}, it gives intial guesses (IGs) for the
    calcium dynamics parameters (\code{log_Ca0}, \code{log_dCa},
    \code{log_tau} for a monoexponential decay, eventually \code{mu} and
    \code{log_dtau} for a biexponential decay). If \code{ratio} is a
    \code{"ratio_fit"} object, the fitted values of its parameters are
    used as initial guesses. If \code{ratio} is none of the above
    objects, a ratiometric fit is performed on \code{df} to find initial
    guesses for the calcium dynamics parameters}
  \item{type}{a character string (either \code{"mono"} or \code{"bi"}),
    specifying the type of calcium exponential decay to consider}
  \item{Plot}{a logical value. Set to \code{TRUE} to plot the original
    signals, the initial guess and the fit results}
  \item{Fit}{a logical value. Set to \code{TRUE} to perform the fit, or
    to \code{FALSE} to compute an initial guess only}
  \item{AfterPeak}{a logical or numerical value. Set to \code{FALSE} to perform the
    fit on the whole fluorescence transients, to \code{TRUE} to consider
    only the part before the fluorescence jump and the convex part after the
    fluorescence peak (for both signals), or to an integer to skip a
    given number of samples after the fluorescence jump}
  \item{Trace}{a logical value. Set to \code{TRUE} to print results for
    successive steps of the optimization algorithm}
  \item{WarnOnly}{a logical value. Set to \code{TRUE} to go on even if
    the fit produced an error}
}
\details{
  The fit is performed using the \code{nls} function, which determines
  the nonlinear (weighted) least-squares estimates of the parameters of
  a nonlinear model. The algorithm is set to the default
  Gauss-Newton.
  
  The initial guesses for the experiment-specific parameters are
  calculated with the \code{\link{igDirect}} function. If the \code{"USE_se"}
  field of the calibration parameters and \code{alpha} are set to
  \code{TRUE}, the initial guesses for these parameters are given by
  their experimental mean value.
  
  The quality of the direct fit is based on the probabilistic properties
  of the fluorescence signals, which are described as realizations of
  Poisson processes. For values of parameter above about 10, the
  Poissonian distribution can be approximated by a Gaussian distribution
  with variance equal to the mean. Applying the square root
  transformation (by setting the \code{SQRT} logical argument to
  \code{TRUE}) to such a process leads to a stabilization of the
  variance, which becomes equal to \eqn{\frac{1}{4}}{1/4}. Thus, one is
  brought back to a standard nonlinear regression setting. Moreover, in
  this situation, it is possible to account for the limited precision
  with which the calibration parameters are known. For this purpose, these
  parameters are also fitted, and a weight of
  \eqn{\frac{1}{\sigma_{exp}^2}}{1/sigma_exp^2} is applied to each of
  them (\eqn{\frac{1}{\sigma_{exp}}}{1/sigma_exp} referring to the
  experimental standard error to the mean (sem)). The whole signal to fit is thus the following:
  \deqn{%
    \left(B_{340}, F_{340}, B_{380}, F_{380}, R_{min}, R_{max}, K_{eff}, K_d
    \right),
  }{c(B_340, F_340, B_380, F_380, R_min, R_max, K_eff, K_d),
  }
  with the following weights:
  \deqn{%
    \left(4, 4, 4, 4, \frac{1}{\sigma_{Rmin}^2},
    \frac{1}{\sigma_{Rmax}^2}, \frac{1}{\sigma_{Keff}^2},
    \frac{1}{\sigma_{Kd}^2}
    \right)
  }{c(4, 4, 4, 4, 1/sigma_Rmin^2, 1/sigma_Rmax^2, 1/sigma_Keff^2,
    1/sigma_Kd^2)
  }
}
\value{
  An object that inherits from both \code{"nls"} and
  \code{"direct_fit"} classes. The object has the following attributes:
  \item{"Name"}{a character string telling which type of fit has been performed}
  \item{"Time"}{the whole time vector, which includes
    \code{NA}s/\code{NaN}s for the background fluorescence signals and
    the supplementary calibration parameters, when relevant}
  \item{"RawData"}{the raw signal, which is created by the concatenation
  of the background fluorescence at 340 nm, the fluorescence transient
  at 340 nm, the background fluorescence at 380 nm, the fluorescence
  transient at 380 nm, and, when relevant, the mean values of the
  selected calibration parameters. This signal is the one passed to the
  \code{nls} formula}
  \item{"RawDataFrame"}{a copy of the input data frame}
  \item{"FitFunction"}{the function passed to the \code{nls} formula}
  \item{"Subset"}{the indices of the \code{Time} vector used for the fit}
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{transientConvexPart}}, 
  \code{\link{mkFluo4DirectFit}}, 
  \code{\link{igDirect}},
  \code{\link{ratioExpSimul}},
  \code{\link{ratioExpPhysio}}
}
\examples{
## Direct Fit On Simulated Data
## (parameters set to the value in Table 2 from
## Joucla et al. (2009) (Journal of Neurophysiology))
## ==================================================

## Parameters of the monoexponential calcium transient
tOn  <- 1
Time <- seq(0,12,length.out=160)
Ca0  <- 0.10
dCa  <- 0.25
tau  <- 1.5

## Calibrated parameters
R_min <- list(value=0.136, mean=0.136, se=0.00363, USE_se=TRUE)
R_max <- list(value=2.701, mean=2.701, se=0.151,   USE_se=TRUE)
K_eff <- list(value=3.637, mean=3.637, se=0.729,   USE_se=TRUE)
K_d   <- list(value=0.583, mean=0.583, se=0.123,   USE_se=TRUE)

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
df <- ratioExpSimul(nb_B    = nb_B,
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

## Perform a monoexponential and a biexpoential ratiometric fit
direct_fit <- directFit(df = df,
                        transients = 1,
                        SQRT = TRUE,
                        ratio = NULL,
                        type = "mono")

## Plot the raw and fitted data as well as goodness of fit tests
plot(direct_fit, numTransient=1, items=1:6)

## Direct Fit On Physiological Data
## (reproduces Fig. 6 of Joucla et al. (2009))
## ===========================================

## Load the data from cockroach olfactory interneurons
data(inVitro)

## Calibrated parameters
R_min <- list(value=0.136, mean=0.136, se=0.00363, USE_se=TRUE)
R_max <- list(value=2.701, mean=2.701, se=0.151,   USE_se=TRUE)
K_eff <- list(value=3.637, mean=3.637, se=0.729,   USE_se=TRUE)
K_d   <- list(value=0.583, mean=0.583, se=0.123,   USE_se=TRUE)

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

## Fit the physiological data with the direct method
## (skip 14 samples of the signal after tOn)
physiolDirectFit <- directFit(physioData,
                              transients=2,
                              SQRT=TRUE,
                              type="mono",
                              AfterPeak=14)

## Plot the raw and fitted data as well as goodness of fit tests
plot(physiolDirectFit, numTransient=2, items=1:6)
}
