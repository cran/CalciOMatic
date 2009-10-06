\name{plotCalciOMatic}
\alias{plotCalciOMatic}
\title{Low-Level Plot Function for CalciOMatic}
\description{
  The low-level function \code{plotCalciOMatic} performs different
  kinds of predefined plots, depending on the value of \code{n}. It is
  generally called by plot methods dedicated to several types of objects
  (\code{"fluo_rawdata"}, \code{"ratio_fit"}, \code{"ratio_fit_list"} or
  \code{"direct_fit"})
}
\usage{
plotCalciOMatic(x = NULL, y = NULL, n = 1, x2 = NULL, y2 = NULL,
                col = "black", col2 = "darkgray",
                main = "MyCalciumRatiometricFit",
                xlab = "", ylab = "", lab = "A", ylas = 1,
                oma = c(4, 0, 1, 0), mar = c(0, 7, 2, 0),
                ask = FALSE, ...)
}
\arguments{
  \item{x}{the abscissa of the main signal to plot (see details below)}
  \item{y}{the main signal to plot (see details below)}
  \item{n}{an integer, between 1 and 5, telling which kind of plot is
    used (see details below)}
  \item{x2}{the abscissa of the secondary signal to plot (see details
    below)}
  \item{y2}{the secondary signal to plot (see details below)}
  \item{col}{the color of the main signal (either an integer or a
    character string)}
  \item{col2}{the color of the secondary signal (either an integer or a
    character string)}
  \item{main}{a character string specifying the main title of the plot}
  \item{xlab}{a character string specifying the label of the x-axis}
  \item{ylab}{a character string specifying the label of the y-axis}
  \item{lab}{a character string specifying the label at the top-left
    of the plot (generally a letter, a letter with an integer, or an
    expression)}
  \item{ylas}{an integer specifying the orientation of the
    yticks. Possible values are 0 or 3 (vertical), or 1 or 2
    (horizontal)}
  \item{oma}{a vector of length 4 specifying the outer margin of the figure}
  \item{mar}{a vector of length 4 specifying the margin of the panel to add}
  \item{ask}{a logical value. Set to \code{FALSE} to draw all plots on
    the same figure. If set to \code{TRUE}, a single device will be
    opened, in which the first plot will be drawn. Three symbols (left
    arrow, black square and right arrow) will be added at the bottom
    right of the plot, for interactions with the user (see details
    below)}
  \item{\dots}{one or more of the following plot parameters:
    \code{"cex"}, \code{"cex.axis"}, \code{"cex.lab"},
    \code{"cex.main"}, \code{"font"}, \code{"font.axis"},
    \code{"font.lab"}, \code{"font.main"}, \code{"line.xlab"},
    \code{"line.ylab"}, \code{"line.lab"}, \code{"line.main"},
    \code{"adj.main"}, \code{"xlim"}, \code{"ylim"}, \code{"tcl"},
    \code{"mgp.x"}, \code{"mgp.y"}}
}
\details{
  The third argument of the \code{plot.calciomatic} function is an
  integer \code{n} comprised between 1 and 5. The plot drawn directly
  depends on the value of \code{n}, as follows:
  \itemize{
  \item{\code{n=1}}{plot \code{y} vs. \code{x} (type \code{lines} and
                    color \code{col}), and superimpose the plot of
		    \code{y2} vs. \code{x2} (color \code{col2}). Used to
		    plot raw data and fitted data vs. time values}
  \item{\code{n=2}}{plot \code{y} vs. \code{x} (type \code{lines} and
                    color \code{col}), and add a dashed horizontal line
		    at y=0 (color \code{col2}). Used to plot fit
		    residuals}
  \item{\code{n=3}}{plot a bar plot, with \code{y} of class \code{"acf"}
                    (color \code{col}), and add a \eqn{95\%}{95\%}
		    confidence band between 
                    \eqn{
		      y= - 1.96 / \sqrt{\code{length(y[["acf"]])}}
	            }{
                      y= - 1.96 /sqrt(length(y$acf))
	            }
		    and
                    \eqn{
		      y= + 1.96 / \sqrt{\code{length(y[["acf"]])}}
	            }{
                      y= + 1.96 /sqrt(length(y$acf))
	            }
		    (color \code{col2}). Used to plot the
		    auto-correlation function of the residuals}
  \item{\code{n=4}}{plot \code{y} vs. \code{x} (type \code{points} and
                    color \code{col}), add a diagonal dashed line
		    (color \code{col}) and a \eqn{95\%}{95\%}
		    confidence band (color \code{col2}). Used to plot
		    quantile-quantile plots of the fit residuals}
  \item{\code{n=5}}{plot an histogram of \code{x} (color \code{col}),
                    and add vertical dashed lines at -3/-2/-1/0/1/2/3
		    times the standard deviation of \code{x} (color
		    \code{col2}). Used to plot the histogram of the fit
		    residuals}
  }
}
\value{
  When \code{ask} is set to its default \code{FALSE} value, a plot is
  drawn, and nothing else is returned by the function. If \code{ask} is
  set to \code{TRUE}, the action to execute is returned, as a character
  string
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{plot.fluo_rawdata}},
  \code{\link{plot.ratio_fit}},
  \code{\link{plot.ratio_fit_list}},
  \code{\link{plot.direct_fit}}
}
\examples{
## Plot simulated fluorescence raw data. A simpler way would be to
## use the high-level plot.fluo_rawdata function, so, this example
## is just for the sake of the form

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

## Extract relevant data from the data frame
Time <- with(df,Time[!is.na(Time) & lambda==340])
adu_340 <- with(df,adu[!is.na(Time) & lambda==340])
adu_380 <- with(df,adu[!is.na(Time) & lambda==380])

## Plot the fluorescence data in the same figure
par(oma = c(4, 0, 1, 0), mar = c(0, 7, 2, 0))
layout(matrix(c(1,2),ncol=1))

plotCalciOMatic(x = Time, y = adu_340, n = 1, xlab = "",
                ylab = expression(paste(adu[340], "(photons)")),
                lab = "A", main = "Fluorescence Raw Data")

plotCalciOMatic(x = Time, y = adu_380, n = 1, xlab = "Time (s)",
                ylab = expression(paste(adu[380], " (photons)")),
                lab = "B", main = "")
}
