\name{qt.mixture}
\alias{qt.mixture}
\title{
Quantile function
}
\description{
Quantile function of a mixture of t distributions
}
\usage{
qt.mixture(p, pi, mu, s, df, tol=10*.Machine$double.eps)
}
\arguments{
  \item{p}{Value at which the Quantile function shall be evaluated (numeric scalar within (0, 1))}
  \item{pi}{Mixture proportions (numeric vector of length K > 0; sum(pi) = 1, entries non-negative)}
  \item{mu}{Component medians (numeric vector of lentgh K)}
  \item{s}{Component sigma's (numeric vector of lentgh K; entries positive)}
  \item{df}{Component degrees of freedom (numeric vector of lentgh K; entries positive)}
  \item{tol}{Numerical tolerance of result (numeric scalar, positive, "small")}
}

\value{
\item{Q(y)}{Result}
}

\author{
Lisa Allmesberger
Fabian Bergs
Stefan Immler
Michael Kässmann
}

\seealso{
pt.mixture, bisection.method, weighted.quantile
}
\keyword{Quantile}
