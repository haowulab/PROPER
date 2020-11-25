\name{runSims}
\alias{runSims}

\title{
  Run a number of RNA-seq simulations and DE detections
}
\description{
  This is the "wrapper" function for running RNA-seq DE detection
  simulation. It runs simulations under different sample sizes 
  (replicates in each group) for a certain numbers. In each
  simulation, the RNA-seq data are generated, and then DE detection
  (using user specified method/software) is performed. The return object
  contains DE test results (test statistics, p-values, FDRs, etc.) from
  all simulations. 
}

\usage{
runSims(Nreps = c(3, 5, 7, 10), Nreps2, nsims = 100, sim.opts,
        DEmethod = c("edgeR", "DSS", "DESeq2"), verbose =TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Nreps}{Sample sizes one wants to perform simulation on. This is
    a vector, each element specifies the number of biological replicates
    in each group. Default value is c(3, 5, 7, 10), which means we want
    to perform simulation when there are 3, 5, 7, or 10 replicates in
    each group.}
  \item{Nreps2}{Sample sizes for the second treatment group. If this is
	missing, it'll take the same value as Nreps and then two groups will
	be assumed to have the same sample size. This parameter allows two
	treatment groups have different sample sizes. 	When specified, it
	must be a vector of the same length as Nreps.}
  \item{nsims}{Number of simulations.}
  \item{sim.opts}{An object for simulation option. This should be the
    return object from "RNAseq.SimOptions.2grp" function.}
  \item{DEmethod}{String to specify the DE detection method to be
    used. Available options are "edgeR", "DSS", and "DESeq2". }
  \item{verbose}{Logical value to indicate whether to output some
    messages (progress report).}
}

\details{
  This is the main simulation function in the packge. After specifying
  the simulation parameters (from "RNAseq.SimOptions.2grp" function), one
  wants to evaluate the power vs different sample sizes under that simulation
  setting. This function simulates the count data and performs
  statistical tests for DE detection. It only stores  and returns the DE
  test results (test statistics, p-values, FDRs, etc.) but doesn't make
  inferences. The inferences will be conducted in the "comparePower"
  function. The advantage is that for one simulation setting, 
  the simulation only need to be run once. The inferences using
  different critical values and type I error controls can then be drawn
  from the same results. This greatly save the computation because the
  simulation part is the most computationally intensive.

  This function can be slow, depends on the setting (number of genes,
  replicates, simulations, etc). For the default (50000 genes, 100
  simulations, for 3, 5, 7, or 10 replicates), it takes about an hour to
  run on a single core of i7 2.7GHz CPU. But again, this only need to be
  run once for a particular simulation setting. 
}

\value{
  A list with following fields:
  \item{pvalue, fdrs}{3D array for p-values and FDR from each
    simulation. The dimension of the array is ngenes * N * nsims. Here N
    is length(Nreps), of the number of different sample sizes settings.}
  \item{xbar}{3D array for average read counts for genes. Dimension is
    the same as pvalue/fdr. This will be used in "comparePower" function
    for stratified power quantities. }
  \item{DEid}{A list of length nsims. Each contains the id of DE genes.}
  \item{lfcs}{A list of length nsims. Each contains the log fold changes
    of DE genes.}
  \item{Nreps}{The input Nreps.}
  \item{sim.opts}{The input sim.opts.}
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\seealso{
RNAseq.SimOptions.2grp, simRNAseq
}

\examples{
simOptions = RNAseq.SimOptions.2grp()
## using 3 different sample sizes, run 2 simulations, using edgeR
simRes = runSims(Nreps=c(3,5,7), sim.opts=simOptions, nsims=2,
                 DEmethod="edgeR")
names(simRes)

\dontrun{
## using 5 different sample sizes, run 100 simulations, using edgeR.
## This will be slow.
simRes = runSims(Nreps=c(3,5,7,10,15), sim.opts=simOptions, nsims=100,
                 DEmethod="edgeR")

## for different sample sizes in two groups
simRes = runSims(Nreps=c(3,5,7,10), Nreps2=c(5,7,9,11),
                 sim.opts=simOptions, nsims=100, DEmethod="edgeR")

}
}
