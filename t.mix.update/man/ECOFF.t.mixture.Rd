\name{ECOFF.t.mixture}
\alias{ECOFF.t.mixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes ECOFF for given parameter values
}
\description{
computes ECOFF for given paramter values.
Returns group index (integer scalar within [1,k], ECOFF (numeric scalar))
}
\usage{
ECOFF.t.mixture(par=NULL, pi=NULL, mu=NULL, s=NULL, df=NULL,
                            quantile=0.01, pi_cutoff = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{Matrix of parameters (nmeric 4xK - matrix); may also be delivered separately.}
  \item{pi}{Mixture proportions (numeric vector of length K > 0; sum(pi) = 1, entries non-negative!)}
  \item{mu}{Component medians (numeric vector of length K)}
  \item{s}{Component sigma's (numeric vector of length K; entries positive.)}
  \item{df}{Component degrees of freedom (numeric vector of length K; entries positive)}
  \item{quantile}{Cut-off value (numeric scalar within (0,1))}
  \item{pi_cutoff}{Lower bound for groupsize of 'wild type' (numeric scalar within (0, 1))}
}

\value{
result: group index (integer scalar within [1, K], ECOFF (numeric scalar))
}

\author{
Lisa Allmesberger
Fabian Bergs
Stefan Immler
Michael Kässmann
}
\seealso{
qt.mixture
}
