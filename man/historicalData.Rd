\name{param}
\alias{param}
\alias{cheung}
\alias{gilad}
\alias{bottomly}
\alias{maqc}
\alias{GE.human}
\alias{pbmc}

\docType{data}

\title{
  Some existing RNA-seq and gene expression microarray data
}

\description{
  These are several datasets distributed with the package. They are
  simulation parameters estimated from existing RNA-seq data. 
  The original raw count data were obtained from reCount website.
}

\usage{
data(cheung)
data(gilad)
data(bottomly)
data(maqc)
data(GE.human)
data(pbmc)
}

\details{
  The data include:
  \itemize{
	\item cheung: parameters from Cheung data. They are measurements
	from unrelated individuals, so the dispersions are large.
	\item gilad: parameters from Gilad data. They are for Human liver
	sample comparisons between male and female. This dataset has
	moderate dispersions.
	\item bottomly: parameters from Bottomly data. They compare two
	strains of inbred mice and the within group dispersions are small.
	\item maqc: MAQC data which are technical replicates. There are
	no biological variation from the replicates, so the dispersions are
	close to 0.
	\item GE.human: a vector of aggregated counts from all samples in
	Cheung data. This is used for generating marginal expression
	distribution when provide sequnencing depth.
	\item pbmc: data from gene expression microarray data. This is for
	an example of using historical data to estimate effect sizes.
  }
}

\author{
  Hao Wu <hao.wu@emory.edu>
}


