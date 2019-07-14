\name{BIC_t_mixture}
\alias{BIC_t_mixture}
\title{
Calculate the BIC
}
\description{
Calculate the BIC
}
\usage{
BIC_t_mixture(n, j, pi, mu, s, df)
}
%- maybe also 'usage' for other objects documented here.
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
