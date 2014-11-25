\name{RNAseq.SimOptions.2grp}
\alias{RNAseq.SimOptions.2grp}
\title{
Set up options for simulating RNA-seq data in two-group comparison.
}
\description{
This function takes user provided options for for simulating
RNA-seq data, and return a list. The result of  this function will be
the input for "runSims" and "simRNAseq" function. 
}
\usage{
RNAseq.SimOptions.2grp(ngenes, seqDepth, lBaselineExpr, lOD, p.DE, lfc, sim.seed)
}

\arguments{
  \item{ngenes}{Number of genes in the simulation.}
  \item{lBaselineExpr}{log baseline expression for each gene. This can be: (1) a
    constant; (2) a vector with length ngenes; (3) a function that
    takes an integer n, and generate a vector of length n; (4) a string
    specifying the name of existing datasets, from which the mean
    expressions will be sampled. Details for this option is provide in
	"Details" section.}
  \item{lOD}{log over-dispersion for each gene. Available options are
    the same as for lBaselineExpr.}
  \item{seqDepth}{Sequencing depth, in terms of total read counts. This
    will be ignored if lBaselineExpr is specified.}
  \item{p.DE}{Percentage of genes being differentially expressed
  (DE). By default it's 5\%.}
  \item{lfc}{log-fold change for DE genes. This can be: (1) a
    constant; (2) a vector with length being number of DE genes; (3) a function that
    takes an integer n, and generate a vector of length n.
	If the input is a vector and the length is not the number of DE
	genes, it will be sampled with replacement to generate log-fold change.}
  \item{sim.seed}{Simulation seed.}
}

\value{
  A list with following fields:
  \item{ngenes}{An integer for number of genes. }
  \item{p.DE}{Percentage of DE genes.}
  \item{lBaselineExpr}{A vector of length ngenes for log baseline expression.}
  \item{lOD}{A vector of length ngenes for log over-dispersion.}
  \item{lfc}{A vector of length (ngenes*p.DE) for log fold change of the
    DE genes.}
  \item{sim.seed}{The specified simulation seed.}
  \item{design}{A string representing the experimental design. From this
    function it is '2grp', standing for two-group comparison.}
}

\details{
  The simulation of RNA-seq data requires a lot of parameters. This
  function provides users an interface to specify the simulation
  parameters. The result from this function will be used for simulating
  RNA-seq count data. By default, the simulation parameters are similar to
  that from Cheung data (for unrelated individuals, with large
  biological variance).

  The baseline expression levels and log over-dispersions can be sampled
  from real data. There are parameters estimated from several real
  datasets distributed with the package. Available string options for
  "lBaselineExpr" and "lOD" include: (1) "cheung": parameters from
  Cheung data, which measures unrelated individuals, so the dispersions
  are large; (2) "gilad": from Gilad data which are for Human liver
  sample comparisons between male and female. This dataset has moderate
  dispersions; (3) "bottomly": from Bottmly data which are from
  comparing two strains of inbred mice. The dispersions are small. (4)
  "maqc": from MAQC data which are technical replicates. There are
  no biological variation from the replicates because the data are
  technical replicates. The dispersions from this dataset is very
  small.

  The effect sizes (log fold changes of the DE genes) are arbitrarily
  specified. It is possible to estimate those from real data. We provide
  a simple example in the package vignette for doing so. 
}


\examples{
## default
simOptions=RNAseq.SimOptions.2grp()
summary(simOptions)

## specify some parameters: generate baseline expression and
## dispersion from Bottom data, and specify a function for
## alternative log fold changes.
fun.lfc=function(x) rnorm(x, mean=0, sd=1.5)
simOptions=RNAseq.SimOptions.2grp(ngenes=30000, lBaselineExpr="bottomly",
       lOD="bottomly", p.DE=0.05, lfc=fun.lfc)
summary(simOptions)
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\seealso{
simRNAseq,runSims
}
