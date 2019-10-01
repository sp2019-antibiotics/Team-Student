\name{BIC.t.mixture}
\alias{BIC.t.mixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the BIC
}
\description{
Calculate the BIC
}
\usage{
BIC.t.mixture(n, j, pi, mu, s, df)
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




\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}

