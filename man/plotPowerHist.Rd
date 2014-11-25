\name{plotPowerHist}
\alias{plotPowerHist}

\title{
  Plot the histogram of power 
}
\description{
  Using histogram-alike graph to visualize the distribution of all genes
  and DE genes in all user specified strata. 
}
\usage{
plotPowerHist(powerOutput, simResult, main = "Histogram of power", return = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{powerOutput}{Result object from "comparePower" function.}
  \item{simResult}{Result object from "runSims" function.}
  \item{main}{Figure caption.}
  \item{return}{Return a matrix for average number of genes and true DE
    genes in each strata.}
}

\value{
  A matrix of two rows. First row is for number of genes, and the second
  row is for number of DE genes. Columns are strata.
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\examples{
\dontrun{
simOptions = RNAseq.SimOptions.2grp()
## run a few simulations
simRes = runSims(Nreps=c(3,5,7), sim.opts=simOptions, nsims=5,
                 DEmethod="edgeR")

## using FDR 0.1 to call DE, then look at power curves and summary
powers = comparePower(simRes)
plotPowerHist(powers, simRes)
}
}
