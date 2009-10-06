\name{inVitro}
\alias{inVitro}
\docType{data}
\title{Calcium Transients and Fura 2 Loading Curves}
\description{
  Data from 20 \emph{in vitro} experiments performed on cockroaches
  (\emph{Periplaneta americana}) antennal lobe neurons. Each data set
  contains a list of experimental results made of the fura 2 loading
  curve (measured at 360 nm) and 3 (or 2) calcium transients (measured
  at 340 and 380 nm).
}
\usage{data(inVitro)}
\format{
  The data are stored as \code{integer}. \code{inVitro} contains 20
  sublists. Each sublist corresponds to one experiment and contains:
  \tabular{ll}{
    \code{time360} \tab A vector of times (in s) at which fluorescence
    measurements \cr
    \tab at the isosbestic wavelength (360 nm) were taken\cr
    \code{adu360} \tab A vector of fluorescence measurements at 360
      nm\cr
    \code{adu360Background} \tab Background fluorescence measurements at 360
      nm\cr
    \code{P360} \tab Number of on-chip binned pixels used for the
      measurements at 360 nm\cr
    \code{P360Background} \tab Number of on-chip binned pixels used for the
    measurements of background \cr
    \tab fluorescence at 360 nm. If different from \code{P360} then the
    pixels used for \cr
    \tab background measurement were different than the one
    used for \cr
    \tab loading curve measurements\cr
    \code{onChipBinning360} \tab Number of physical pixels binned
      on-chip for the 360 nm measurements\cr
    \code{exposureTime360} \tab Exposure time used at 360 nm (in s)\cr
    \code{adu340Background} \tab Background fluorescence measurements at 340
      nm\cr
    \code{adu380Background} \tab Background fluorescence measurements at 380
      nm\cr
    \code{furaPipette} \tab The total fura concentration in the
      pipette (in \eqn{\mu}{micro}M)\cr
    \code{exposureTime340} \tab Exposure time used at 340 nm (in s)\cr
    \code{exposureTime380} \tab Exposure time used at 380 nm (in s)\cr
    \code{P} \tab Number of on-chip binned pixels\cr
    \code{PBackground} \tab Number of on-chip binned pixels used for the
    measurements of background \cr
    \tab fluorescence at 340 and 380 nm. If different
    from \code{P} then the pixels used for \cr
    \tab background measurement
    were different than the one used for \cr
    \tab transient measurements\cr
    \code{onChipBinning} \tab Number of physical pixels binned
      on-chip for the 340 and 380 nm measurements\cr    
    \code{stim1, stim2 and stim3} \tab List with calcium
    transient data (see below)\cr
    }

  Each element \code{stim1}, \code{stim2} and when a third transient was
  measured, \code{stim3}, are lists with the following components:
  \tabular{ll}{
    \code{time} \tab A vector of times (in s) at which fluorescence
      measurements were taken\cr
    \code{adu340} \tab A vector of fluorescence measurements at 340
      nm\cr
    \code{adu380} \tab A vector of fluorescence measurements at 380
    nm\cr
    }
}
\details{
  Break-in time is 0.

  \code{inVitro} contains data recorded in vitro (from cultured cells)
  from antennal lobe neurons for a total of 20 experiments. Experiments
  5 and 16 in the \code{inVitro} data set have only 2 transients instead
  of 3.
}
\source{
  Andreas Pippow \email{andreas.pippow@uni-koeln.de} did the experiments
}
\examples{
data(inVitro)
}
\keyword{datasets}
