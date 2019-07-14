\name{dt.mixture_log}
\alias{dt.mixture_log}
\title{
Distribution of a Log Mixture Function
}
\description{
Log PDF of a mixture t distributions
}
\usage{
dt.mixture_log(y, pi, mu, s, df, log=FALSE)
}
\arguments{
  \item{y}{Value(-s) at which the PDF shall be evaluated (numeric vector)}
  \item{pi}{Mixture proportions (numeric vector of length K > 0; sum(pi) = 1, entries non-negative)}
  \item{mu}{Component medians (numeric vector of lentgh K)}
  \item{s}{Component sigma's (numeric vector of lentgh K; entries positive)}
  \item{df}{Component degrees of freedom (numeric vector of lentgh K; entries positive)}
  \item{log}{if log=TRUE it plots a logarithm plot}
}

\value{
\item{f(y)}{Result}
}


\author{
Lisa Allmesberger
Fabian Bergs
Stefan Immler
Michael Kässmann
}

\seealso{
qt.mixture, dt.mixture, pt.mixture, rt.mixture
}

\keyword{dt.mixture}
