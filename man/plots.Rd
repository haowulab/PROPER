\name{plots}
\alias{plotPower}
\alias{plotPowerTD}
\alias{plotPowerFD}
\alias{plotFDR}
\alias{plotFDcost}
\alias{plotPowerAlpha}
\alias{plotAll}

\title{
  Plot the stratified curves for power-related quantities
}
\description{
  These are a group of functions to generate plot to visualize the
  stratified power metrics, including stratified power, FDR, true
  discoveries, false discoveries, false discovery cost, and type I
  error. There will be a line for each sample size. X-axis are strata
  and y-axis are the metrics. Error bars will be plotted if requested.
}

\usage{
plotPower(powerOutput, cols = 1:ncol(powerOutput$FD),
    lty = 1:ncol(powerOutput$power), main = "", ylab = "Power",
    leg = TRUE, error.bar=TRUE)

plotFDR(powerOutput, cols = 1:ncol(powerOutput$FDR), lty=1:ncol(powerOutput$FDR),
    main = "", ylab="FDR", leg = TRUE, error.bar=TRUE)

plotPowerTD(powerOutput, cols = 1:ncol(powerOutput$TD),
    lty=1:ncol(powerOutput$TD), main="", ylab = "# True Discoveries",
    leg = TRUE, error.bar=TRUE)

plotPowerFD(powerOutput, cols = 1:ncol(powerOutput$FD), lty=1:ncol(powerOutput$FD),
    main = "", ylab = "# False Discoverie", leg = TRUE, error.bar=TRUE)

plotFDcost(powerOutput, cols = 1:ncol(powerOutput$FD), lty=1:ncol(powerOutput$FD),
    main="", ylab = "False discovery cost", leg = TRUE, error.bar=TRUE) 

plotPowerAlpha(powerOutput, cols = 1:ncol(powerOutput$alpha), lty=1:ncol(powerOutput$alpha),
    main = "", ylab = "False positive rate", leg = TRUE, error.bar=TRUE)

plotAll(powerOutput)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{powerOutput}{Result object from "comparePower" function.}
  \item{cols, lty}{Colors and line types of the TD curves.}
  \item{main}{Figure title.}
  \item{ylab}{Y-axis label.}
  \item{leg}{Indicator for having figure legend or not.}
  \item{error.bar}{Indicator for whether plotting error bars in the figures.}
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\seealso{
  comparePower
}

\examples{
\dontrun{
simOptions = RNAseq.SimOptions.2grp()
## run a few simulations
simRes = runSims(Nreps=c(3,5,7), sim.opts=simOptions, nsims=5,
                 DEmethod="edgeR")

## using FDR 0.1 to call DE, then look at power curves and summary
powers = comparePower(simRes)
## plot 
par(mfrow=c(2,3))
plotPower(powers)
plotPowerTD(powers)
plotFDR(powers)
plotFDcost(powers)
plotPowerAlpha(powers)

## in one figure
par(mfrow=c(2,3))
plotAll(powers)
}

}
