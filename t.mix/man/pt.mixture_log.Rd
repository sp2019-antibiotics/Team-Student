\name{pt.mixture_log}
\alias{pt.mixture_log}
\title{
CdF of a mixture of t distributions.
}
\description{
CdF of a mixture of t distributions.
}
\usage{
pt.mixture_log(y, pi, mu, s, df, log.p=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Value(-s) at which the CDF shall be evaluated (numeric vector)}
  \item{pi}{Mixture proportions (numeric vector of length K > 0; sum(pi)=1, entries non-negativ!)}
  \item{mu}{Component medians (numeric vector of length K)}
  \item{s}{Component sigma's (numeric vector of length K; entries positive!)}
  \item{df}{Component degrees of freedom (numeric vector of length K; entries positive!)}
  \item{log.p}{if log=TRUE it plots a logarithm plot}
}

\value{
\item{F(y)}{Result}
}

\author{
Lisa Allmesberger
Fabian Bergs
Stefan Immler
Michael Kässmann
}

\seealso{
dt.mixture, qt.mixture, rt.mixture
}


