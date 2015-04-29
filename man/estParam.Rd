\name{estParam}
\alias{estParam}
\title{
  estimate simulation parameters from a count table
}
\description{
  This function estimate simulation parameters from a count table, which
  can be used to set simulation scenarios. The parameters estimated
  include baseline expression, biological variation in the form of
  dispersion paramter.
}
\usage{
estParam(X,type = c(1,2))
}
\arguments{
  \item{X}{An ExpressionSet or a matrix, with entries being RNA sequencing 
    counts}
  \item{type}{\code{type=1} estimates with minimal assumption based on sample
    means and variances for each genes. \code{type=2} uses \code{edgeR}
    package.}
}

\details{
 This is a function that allows a user to establish simulation basis from his/her own data
  set. Type 1 uses simple average count and sample variance to estimate
 over dispersion, after normalizing by library size for each
 sample. Estimated dispersion is bounded at minimum 0.001. Type 2 uses
 \code{edgeR}'s \code{AveLogCPM} as the estiamte for baseline expression
 and \code{tagwiseDispersion} as the estimate for dispersion. 
}

\value{
  A list with following fields:
  \item{seqDepth}{Estimated sequencing depth for each sample}
  \item{lmean}{A vector of baseline expression rate per gene in log
  scale, with length \code{nrow(X)}}
  \item{lOD}{A vector of log dispersion for each gene, with length \code{nrow(X)}}
}

\author{
  Jean Wu <zhijin_wu@brown.edu>
}


\examples{
\dontrun{
lmu0=rnorm(20000,5,2)
lOD0=rnorm(20000,-4,1)
lmu=exp(rnorm(20000*10,lmu0,exp(lOD0)))
X=matrix(rpois(20000*10,lmu),20000,10)
param=estParam(X)
}
}
