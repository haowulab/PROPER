\name{comparePower}
\alias{comparePower}
\title{
  Compute the power-related quantities from simulation results
}
\description{
  This function take the simulation output from "runSims" function and
  compute a variety of power-related quantities.
}
\usage{
comparePower(simOutput, alpha.type = c("fdr", "pval"), alpha.nominal=0.1,
             stratify.by = c("expr", "dispersion"), strata,
             filter.by = c("none","expr"), strata.filtered = 1,
             target.by=c("lfc", "effectsize"), delta=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simOutput}{The result from "runSims" function.}
  \item{alpha.type}{A string to represent the way to call DE genes. 
    Available options are "fdr" and "pval", for calling DE genes based
    on FDR or p-values. Default is "fdr".}
  \item{alpha.nominal}{The nomial value for call DE genes. Recommend
    values are 0.1 when using FDR, and 0.01 for using p-values.}
  \item{stratify.by}{A string to represent the way to stratify
    genes. Available options are "expr" and "dispersion", for
    stratifying gene by average expression levels or over-dispersion.} 
  \item{strata}{The strata used for stratification.}
  \item{filter.by}{A string to represent the way to filter
    genes. This is used in conjunction with strata.filtered for gene
    filtering. Available options are "none" and "expr". "none" stands
    for no filtering, thus all genes will be considered. "expr" stands
    for filtering based on average expression levels.}
  \item{strata.filtered}{The strata to be filtered out in computing
    power-related quantities. Genes fall into these strata will be
    excluded when computing power-related quantities. See "Details" for
    more description of gene filtering.}
  \item{target.by}{A string to specify the method to define
	"biologically important" DE genes. Available options are (1) "lfc": interesting
	genes are defined by absolute log fold changes.  (2) "effectsize":
	interesting genes are defined by absolute log fold changes divided
	by the square root of dispersion.}
  \item{delta}{A threshold used for defining "biologically
    important" genes. Genes with absolute log fold changes (when
	target.by is "lfc") or effect sizes (when target.by is "effectsize") 
	greater than this value are deemed DE in power calculations. See
	"Details" for more description.}
}

\details{
  This is the main function to compute various power-related quantities,
  under stratification and filtering of all genes.

  Gene stratification: we advocate to compute and visualize the powers
  at different stratification of genes. Because of the characteristics
  of RNA-seq data such as average expression level (counts), powers for
  calling DE genes are different even when the magnitude of changes are
  the same. The stratified results will provide a more comprehensive
  power assessment and better guide the investigators in experimental
  designs and analysis strategies.

  Gene filtering: sometimes it is advisible to filter out some genes
  (such as the ones with very low count) before DE detection. The
  filtering option here provides an opportunity to compare the powers
  before and after filtering.

  Define biologically interesting genes: we advocate to compute powers
  for "biologically interesting genes", because there might be many true
  DE genes with very low changes which are difficult to detect. In this
  sense, we are only interested in genes with adequate changes between
  groups (defined as "biologically interesting"). We provide two options
  to define biologically interesting genes: by absolute values of log fold
  changes or effect sizes (absolute values of log fold changes divided
  by the square root of dispersions). Genes with these quantities over
  a threshold are deemed intersting, and the power calculation are based
  on these genes.
}

\value{
  A list with following fields:
  \item{TD, FD, alpha, FDR, power}{3D array representing the number of true
    discoveries, false discoveries, type I error rate, FDR, and power for all simulation
    settings. The dimension of the arrays are nstrata * N * nsims. Here
    nstrata is number of specified strata. N is number of different
    sample sizes settings, and nsims is number of simulations.}
  \item{alpha.marginal, FDR.marginal, power.marginal}{Matrix
    representing the marginal type I error rate, FDR, and power for all
    simulation settings. The dimension of the matrices are N *
    nsims. These are the marginalized (over all strata) values of the
    stratified quantities.}
  \item{stratify.by}{The input stratify.by.}
  \item{strata}{The input strata.}
  \item{Nreps}{Sample sizes one wants to perform simulation on. This is
    taken from the simulation options.}
  \item{target.by}{The input method to define "biologically important"
	DE genes.} 
  \item{delta}{The input delta for biologically important genes.}
  
}

\author{
  Hao Wu <hao.wu@emory.edu>
}


\seealso{
RNAseq.SimOptions.2grp, simRNAseq, runSims
}

\examples{
\dontrun{
simOptions = RNAseq.SimOptions.2grp()
## run a few simulations
simRes = runSims(Nreps=c(3,5,7), sim.opts=simOptions, nsims=5,
                 DEmethod="edgeR")

## using FDR 0.1 to call DE, then look at power curves and summary
powers = comparePower(simRes)
summaryPower(powers)
par(mfrow=c(2,2))
plotPower(powers)
plotPowerTD(powers)
plotFDR(powers)
plotFDcost(powers)

## filter out the genes with low counts (<10) and redo power calculation
## Marginal powers are significantly higher.
powers = comparePower(simRes, filter.by="expr", strata.filtered=1)
summaryPower(powers)
par(mfrow=c(2,2))
plotPower(powers)
plotPowerTD(powers)
plotFDR(powers)
plotFDcost(powers)

## Provide higher threshold for log fold change to define true DE.
## This will result in higher power.
powers2 = comparePower(simRes, delta=2)
summaryPower(powers2)
par(mfrow=c(2,2))
plotPower(powers2)
plotPowerTD(powers2)
plotFDR(powers2)
plotFDcost(powers2)

## use effect size to define biologically interesting genes
powers3 = comparePower(simRes, filter.by="expr", strata.filtered=1,
                      target.by="effectsize", delta=1)
summaryPower(powers3)
par(mfrow=c(2,2))
plotPower(powers3)
plotPowerTD(powers3)
plotFDR(powers3)
plotFDcost(powers3)
}
}
