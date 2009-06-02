\name{anova4Fits}
\alias{anova4Fits}
\title{Perform an ANalysis Of VAriance between two fit objects}
\description{
  The function \code{Anova_4_Fits} performs an ANOVA between two objects
  inheriting from the \code{"nls"} class, in order to determine which
  one is best adapted to fit the raw data
}
\usage{
anova4Fits(Fit_1, Fit_2)
}
\arguments{
  \item{Fit_1}{the first \code{"nls"} object to compare}
  \item{Fit_2}{the second \code{"nls"} object to compare}
}
\details{
  The sums of the square residuals of both models are compared, the
  least of both tells which model is the most appropriate to fit the raw
  data
}
\value{
  An integer (1 or 2) indicating which model is best appropriate to fit
  the raw data
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{directFit}}}
\examples{
## Parameters of the biexponential calcium transient
tOn  <- 1
Time <- seq(0,30,0.1)
Ca0  <- 0.10
dCa  <- 0.25
tau  <- 1.5
mu   <- 0
dtau <- 10

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

## Create a biexponential calcium decay
Ca_Bi <- caBiExp(t = Time, tOn = tOn,
                 Ca0 = Ca0, dCa = dCa, tau = tau,
                 fact=1/(1+exp(-mu)), dtau = dtau)

## Simulate the corresponding ratiometric experiment
df_Bi <- ratioExpSimul(nb_B    = nb_B,
                       Ca      = Ca_Bi,
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
direct_fit_mono <- directFit(df = df_Bi,
                             transients = 1,
                             SQRT = TRUE,
                             ratio = NULL,
                             type = "mono")
direct_fit_bi   <- directFit(df = df_Bi,
                             transients = 1,
                             SQRT = TRUE,
                             ratio = NULL,
                             type = "bi")

## Test which model ('mono' or 'bi') bests predicts the 'experimental' data
idx <- anova4Fits(Fit_1 = direct_fit_mono, Fit_2 = direct_fit_bi)
print(idx)
}
