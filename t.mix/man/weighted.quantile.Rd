\name{weighted.quantile}
\alias{weighted.quantile}

\title{
Weighted quantile
}
\description{
Computes weighted quantile
}
\usage{
weighted.quantile(x, w=rep(1/length(x), times=length(x)), probs=0.5)
}

\arguments{
  \item{x}{Data (numeric vector of length N)}
  \item{w}{Weights (numeric vector of length N, entries non-negative)}
  \item{probs}{Which quantile(-s) shall be computed (numeric vector, entries within (0, 1))}
}

\value{
\item{q}{Weighted sample quantile(-s)}
}

\author{
Lisa Allmesberger
Fabian Bergs
Stefan Immler
Michael Kässmann
}

\keyword{Weighted quantile}
