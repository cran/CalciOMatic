\name{plot.ratio_fit_list}
\alias{plot.ratio_fit_list}
\title{Plot Method for Objects of Class "ratio\_fit\_list"}
\description{
  The function \code{plot.ratio_fit_list} performs different kinds of
  predefined plots for objects of class \code{ratio_fit_list}
}
\usage{
\method{plot}{ratio_fit_list}(x, y = NULL, numTransient = 1,
     items = 1:5, col = "black", col2 = "darkgray",
     main = "Intracellular calcium transient: Ratiometric fit",
     xlabs = c("Time (s)", "Time (s)", "Lag",
               "Theoretical quantiles", "Residuals"),
     ylabs = c(expression(paste("[", Ca^{2 + phantom()},
               "] (", mu, "M)")), "Residuals",
               "Autocorrelation function",
               "Sample quantiles", "Counts"),
     labs = c(expression(A[1]), expression(A[2]),
              "B", "C", "D"),
     ylas = 1, ask = FALSE, ...)
}
\arguments{
  \item{x}{an object of class \code{ratio_fit_list}, as
    returned by the \code{ratioFitFromDf} function}
  \item{y}{argument not used (NULL by default)}
  \item{numTransient}{a vector of integers specifying which of the
    fitted transients whould be plotted}
  \item{items}{a vector of integers (between 1 and 5), telling which
    panels to draw}
  \item{col}{the color of the main signals to plot (either an integer or a
    character string)}
  \item{col2}{the color of the secondary signals to plot (either an
    integer or a character string)}
  \item{main}{a character string specifying the main title of the
    figure}
  \item{xlabs}{a vector of character strings specifying the label of the
    x-axes}
  \item{ylabs}{a vector of character strings specifying the label of the
    x-axes}
  \item{labs}{a vector of character strings specifying the label at the
    top-left of each panel (generally a letter, a letter with an integer, or an
    expression)}
  \item{ylas}{an integer specifying the orientation of the
    yticks. Possible values are 0 or 3 (vertical), or 1 or 2
    (horizontal)}
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
  If the user does not want to draw all plots in the same figure, (s)he
  can set to logical \code{ask} value to \code{TRUE}. In that case,
  three symbols (left arrow, black square and right arrow) will be added
  at the bottom right of the plot, allowing user interactions. By
  clicking on the arrowhead oriented left (resp. right), the user will
  draw the previous (resp. next) plot (within \code{items}). By clicking
  on the black square, the user will stay on the current plot and none
  of the symbols will be available anymore.
  If \code{ask} is set to \code{FALSE}, \eqn{N}{N} =
  \code{length(numTransient)} figures will be opened simultaneously
  within a device similar to the one in which the first figure was
  plotted
}
\value{
  This plot method does not return anything else that the plotted
  figures
}
\author{Sebastien Joucla \email{sebastien.joucla@parisdescartes.fr}}
\seealso{\code{\link{plotCalciOMatic}},
  \code{\link{plot.fluo_rawdata}},
  \code{\link{plot.direct_fit}},
  \code{\link{plot.ratio_fit}}
}
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

## Perform a ratiometric fit
physioRatioFit <- ratioFitFromDf(df = physioData,
                                 transients = 2,
                                 AfterPeak = 14)

## Print the class of 'physioRatioFit'
print(class(physioRatioFit))

## Create the data frame containing the physiological data
## (experiment #2, stimulations #2 and #3)
## G and s_ro are the respectively the gain of the CCD camera
## and the standard deviation of its read-out process
physioData <- ratioExpPhysio(dataset="inVitro",
                             expe=2, stim=c(2,3),
                             idxOn=10,
                             R_min=R_min, R_max=R_max,
                             K_eff=K_eff, K_d=K_d,
                             G=0.146, s_ro=16.4,
                             alphamethod=FALSE)

## Perform a ratiometric fit
physioRatioFit <- ratioFitFromDf(df = physioData,
                                 transients = c(2,3),
                                 AfterPeak = 14)

## Print the class of 'physioRatioFit'
print(class(physioRatioFit))

## Plot the results on the same figure
## plot(x = physioRatioFit, numTransient = c(2,3), items=1:4)

## Plot the results on separate figures and wait for a user action
## Click on the black square to make the second figure appear
## plot(x = physioRatioFit, numTransient = c(2,3), items=1:4, ask=TRUE)
}
