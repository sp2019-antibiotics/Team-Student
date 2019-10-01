\name{AIC.t.mixture}
\alias{AIC.t.mixture}
\title{
Calculate the AIC
}
\description{
Calculate the AIC
}
\usage{
AIC.t.mixture(n, j, pi, mu, s, df)
}
\arguments{
   \item{n}{Frequencies (integer vector of length J; non-negative entries)}
 \item{j}{Possible values (numeric vector of length J > 0)}
  \item{pi}{Mixture proportions (numeric vector of length K > 0; sum(pi)=1, entries non-negativ!)}
  \item{mu}{Component medians (numeric vector of length K)}
  \item{s}{Component sigma's (numeric vector of length K; entries positive!)}
  \item{df}{Component degrees of freedom (numeric vector of length K; entries positive!)}
}

\author{
Lisa Allmesberger
Fabian Bergs
Stefan Immler
Michael Kässmann
}
