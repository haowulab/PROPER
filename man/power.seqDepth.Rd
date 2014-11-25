\name{power.seqDepth}
\alias{power.seqDepth}

\title{
  Estimate the marginal power under different sequencing depth
}

\description{
  The function helps study the power under different sequencing depth
  and sample sizes. It estimates the marginal powers based on existing
  simulation results, given new sequencing depth.
}

\usage{
power.seqDepth(simResult, powerOutput, depth.factor = c(0.2, 0.5, 1, 2, 5, 10))
}


\arguments{
  \item{simResult}{Result object from "runSims" function.}
  \item{powerOutput}{Result object from "comparePower" function.}
  \item{depth.factor}{A vector of numbers specifying the *relative* sequencing depth,
    comparing to the depth used in the simulation. 1 means using the
    same number of total reads as the simulation. }
}

\details{
  The powers under different sequencing depth and sample sizes provides
  important guidence in experimental design. Under the same total number
  of sequence reads, investigator can choose to use more replicates and
  shallower coverage for each, or less replicates and deeper coverage.

  This function provides estimated marginal power holding all
  experimental variables fixed (biological variation, effect sizes,
  sample sizes, etc.) except the sequencing depth. Changing sequencing
  depth will only alter the marginal distribution of average
  counts. Since the stratified power (by average counts) won't change,
  those numbers are used in estimating the powers under different
  depth. This approaches allows skipping new simulations, which saves
  computation. 
}

\value{A matrix for marginal powers. Each row is for a 
  sequencing depth, each columns is for a sample size.}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\seealso{
  comparePower, summary.power
}

\examples{
\dontrun{
simOptions = RNAseq.SimOptions.2grp()
simRes = runSims(Nreps=c(3,5,7), sim.opts=simOptions, nsims=5,
                 DEmethod="edgeR")
powers = comparePower(simRes)
power.seqDepth(simRes, powers)
}
}
