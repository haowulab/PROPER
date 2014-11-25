\name{summaryPower}
\alias{summaryPower}

\title{
  Generate a summary table for power analysis result
}
\description{
  Takes a result object from "comparePower" function and print out a
  table to summarize important power-related quantities. The results are
  marginalized, meaning that they are averaged quantities over all
  strata and simulations. This provides a quick view of the marginal
  results for power-related quantities under all sample sizes. 
}

\usage{
summaryPower(powerOutput)
}

\arguments{
  \item{powerOutput}{The result object from "comparePower" function.}
}

\value{
  A matrix with multiple rows, each for a different sample size. Columns
  include sample size, specified nomial type I control value (for FDR or
  p-values), actual FDR, marginal power, averaged number of true and
  false discoveries, and false discovery costs. 
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\seealso{
  runSims, comparePower  
}

\examples{
\dontrun{
simOptions = RNAseq.SimOptions.2grp()
simRes = runSims(Nreps=c(3,5,7), sim.opts=simOptions, nsims=5,
                 DEmethod="edgeR")
powers = comparePower(simRes)
summaryPower(powers)
}
}

