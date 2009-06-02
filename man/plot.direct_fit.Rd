\name{plot.direct_fit}
\alias{plot.direct_fit}
\title{Plot Function for Objects of Class "direct\_fit"}
\description{
  The function \code{plot.direct_fit} performs different
  kinds of predefined plots for objects of class \code{direct_fit}
}
\usage{
\method{plot}{direct_fit}(x, y = NULL, numTransient = 1, items = 1:7,
     col = "black", col2 = "darkgray",
     main = "Fluorescence transients: Direct fit",
     xlabs = c("Time (s)", "Time (s)", "Time (s)", "Time (s)",
               "Lag", "Theoretical quant.", "Residuals"),
     ylabs = c(expression(sqrt(adu[340])),
               expression(res[340]),
               expression(sqrt(adu[380])),
               expression(res[380]),
               "ACF", "Sample quant.", "Counts"),
     labs = c(expression(A[1]), expression(A[2]),
              expression(B[1]), expression(B[2]),
              "C", "D", "E"),
     ylas = 1, ask = FALSE, ...)
}
\arguments{
  \item{x}{a data frame of class \code{direct_fit}, as
    returned by the \code{Direct_Fit_fct} function}
  \item{y}{argument not used (NULL by default)}
  \item{numTransient}{a vector of integers specigying which of the
    fitted transients whould be plotted}
  \item{items}{a vector of integers (between 1 and 7), telling which
    plots to draw}
  \item{col}{the color of the main signals to plot (either an integer or a
    character string)}
  \item{col2}{the color of the secondary signals to plot (either an integer or a
    character string)}
  \item{main}{a character string, the main title of the figure}
  \item{xlabs}{a vector of character strings to add to the x-axes}
  \item{ylabs}{a vector of character strings to add to the y-axes}
  \item{labs}{a vector of character strings to add to the top left of
    each panel}
  \item{ylas}{an integer specifying the orientation of the
    yticks. Possible values are 0 or 3 (vertical), or 1 or 2
    (horizontal)}
  \item{ask}{a logical value (\code{FALSE} by default). If set to
    \code{TRUE}, three or four symbols (among the followings: left
    arrow, red circle, black square and right arrow) are added at the
    bottom right of the plot, on which the user is allowed to click to
    perform any action. In that case, a character string is returned,
    telling what the user asked for. This option is useful when
    \code{plot.calciomatic} is called by a high-level plot method}
  \item{\dots}{one or more of the following plot parameters:
    \code{"cex"}, \code{"cex.axis"}, \code{"cex.lab"},
    \code{"cex.main"}, \code{"font"}, \code{"font.axis"},
    \code{"font.lab"}, \code{"font.main"}, \code{"line.xlab"},
    \code{"line.ylab"}, \code{"line.lab"}, \code{"line.main"},
    \code{"adj.main"}, \code{"xlim"}, \code{"ylim"}, \code{"tcl"},
    \code{"mgp.x"}, \code{"mgp.y"}}
}
\details{
  If the user does not want to draw all plots in the same figure, (s)he
  can set to logical \code{ask} value to \code{TRUE}. In that case,
  three symbols (left arrow, black square and right arrow) will be added
  at the bottom right of the plot, allowing user interactions. By
  clicking on the arrowhead oriented left (resp. right), the user will
  draw the previous (resp. next) plot (within \code{items}). By clicking
  on the black square, the user will stay on the current plot and none
  of the symbols will be available anymore
}
\value{
  This plot function does not return anything else that the plotted figures
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{plotCalciOMatic}},
  \code{\link{plot.fluo_rawdata}},
  \code{\link{plot.direct_fit}},
  \code{\link{plot.ratio_fit_list}}
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
                             alphamethod=TRUE)

## Perform a direct fit
physioDirectFit <- directFit(physioData,
                             transients=2,
                             SQRT=TRUE,
                             type="mono",
                             AfterPeak=14)

## Plot the raw and fitted data as well as goodness of fit tests
## on the same figure
plot(x=physioDirectFit,
     numTransient=2,
     items=1:6)

## Plot the raw and fitted data as well as goodness of fit tests
## on the same figure
## plot(x=physioDirectFit,
##      numTransient=2,
##      items=1:6,
##      ask=TRUE)
}
