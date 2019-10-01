\name{visualization.t.mixture}
\alias{visualization.t.mixture}

\title{
Visualization of the results from a t mixture distribution
}
\description{
Visualization of the results from a t mixture distribution
}
%usage too large
%\usage{
%visualization.t.mixture(n, j, pi, mu, s, df, ecoff=NULL, kernel.density=TRUE, hist=TRUE, new.plot=hist)
%}
\arguments{
  \item{n}{Frequencies (integer vector of length J)
}
\item{j}{Possible values (numeric vector of length J)
}
  \item{pi}{Mixture proportions (numeric vector of length K > 0; sum(pi) = 1, entries non-negative)}
  \item{mu}{Component medians (numeric vector of lentgh K)}
  \item{s}{Component sigma's (numeric vector of lentgh K; entries positive)}
  \item{df}{Component degrees of freedom (numeric vector of lentgh K; entries positive)}
  \item{ecoff}{Should ECOFF be plotted as well? (booelan scalar)}
  \item{kernel.density}{Should the kernel density be plotted?}
  \item{hist}{Plots a histogram}
  \item{new.plot}{Which type should be plotted}
  \item{pi_resistant}{Proportion of resistant observations (numerical scalar within (0, 1))}
}

\value{
None (only plots)
}

\author{
Lisa Allmesberger
Fabian Bergs
Stefan Immler
Michael Kässmann
}
\seealso{
dt.mixture, ECOFF.t.mixture
}
\keyword{visualization}


