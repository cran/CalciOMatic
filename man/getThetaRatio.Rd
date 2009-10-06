\name{getThetaRatio}
\alias{getThetaRatio}
\title{
  Compute confidence intervals for parameters of calcium dynamics model
  fitted with the ratiometric approach
}
\description{
  The function \code{getThetaRatio} estimates confidence intervals for
  parameters a calcium dynamic model fitted to data estimated with the
  ratiometric transformation. The way these CIs are computed depends on
  the \code{ciMode} parameter.
}
\usage{
getThetaRatio(calcium_ratio_fit,
              ciMode = c("normalApprox", "likelihoodRatio"), ...)
}
\arguments{
  \item{calcium_ratio_fit}{an object of class \code{"ratio_fit"}}
  \item{ciMode}{should the normal approximation (\code{"normal"}) or the
    likelihood ratio (\code{"ratio"}) be used to obtain the CI?}
  \item{\dots}{not used}
}
\details{
  The \code{ciMode} argument specifies which approach to use to estimate
  the CIs. If set to \code{"normal"}, the quadratic approximation of the
  log-likelihood applies, and the \eqn{95\%}{95\%} CIs are given as
  \eqn{t(0.975,dof) \times se(p)}{t(0.975,dof)*se(p)}, where \eqn{t}{t}
  is the Student quantile function, \eqn{dof}{dof} is the number of
  degrees of freedom, \eqn{se(p)}{se(p)} is the standard error
  associated to the estimation of parameter \eqn{p}{p} (given by the
  inverse of the square root of the diagonal of the hessian matrix
  returned by \code{"optim"}). If \code{ciMode} is set to
  \code{"likelihoodRatio"}, we make use of the likelihood ratio
  statistics (Davison, 2003).
}
\value{
  A matrix with 2 rows and \code{N} columns, corresponding to the number
  of parameters of the calcium dynamics model. Each column gives the
  lower and upper bound of the 95\% confidence interval for each
  parameter.
}
\references{
  Davison AC (2003), Statistical Models, Cambridge University Press
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\examples{
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

## Compute the confidence interval
## using the likelihood ratio statistics
CI <- getThetaRatio(physioRatioFit,
                    ciMode = "likelihoodRatio")

print(CI)
}