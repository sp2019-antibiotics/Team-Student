\name{rt.mixture}
\alias{rt.mixture}
\title{
Simulates N observations from a mixture of t distributions.
}
\description{
Simulates N observations from a mixture of t distributions.
}
\usage{
rt.mixture(N, pi, mu, s, df, round=FALSE)
}

\arguments{
  \item{N}{Number of realizations (integer scalar, positive)}
  \item{pi}{Mixture proportions (numeric vector of length K > 0; sum(pi) = 1, entries non-negative!)}
  \item{mu}{Component medians (numeric vector of length K)}
  \item{s}{compnent sigma's (numeric vector length K; entries positive)}
  \item{df}{Also called nu; component degrees of freedom (numeric vector of length K; entries positive)}
  \item{round}{Should number be rounded (boolean scalar)}
}

\value{
\item{y}{N accordingly distributed realizations (numeric vector of length N)}
}

\author{
Lisa Allmesberger
Fabian Bergs
Stefan Immler
Michael Kässmann
}

\seealso{
qt.mixture, dt.mixture, pt.mixture
}

