\name{simRNAseq}
\alias{simRNAseq}
\title{
  Simulate RNA-seq count data
}
\description{
  This function generate simulated RNA-seq count data for a two-group
  comparison design. It takes an object for simulation option (return object
  from 'RNAseq.SimOptions.2grp' function, and sample sizes (replicates
  in each condition, then generate a matrix of counts.  
}
\usage{
simRNAseq(simOptions, n1, n2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simOptions}{An object for simulation option. This should be the return object
    from 'RNAseq.SimOptions.2grp' function.}
  \item{n1, n2}{Sample size in two treatment groups.}
}

\details{
  The count data are generated based a negative binomial (NB)
  distribution. Parameters for NB are provided in the input object. 
}

\value{
  A list with following fields:
  \item{counts}{A matrix of dimension ngenes x (n1+n2). Each row is for
  a gene, each column is for a sample. Columns 1 to n1 are for the first
  condition. The rest columns are for the other condition.}
  \item{designs}{A vector/matrix representing the experimental designs. In a
    two-group comparison, it's a simple 0/1 vector.}
  \item{DEid}{The ID (row index) for true DE genes.}
  \item{simOptions}{The simulation option object. This is exactly the
    same as the input simOptions.}
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\seealso{
RNAseq.SimOptions.2grp, runSims
}

\examples{
simOptions = RNAseq.SimOptions.2grp()
data = simRNAseq(simOptions, n1=3, n2=3)
names(data)
}

