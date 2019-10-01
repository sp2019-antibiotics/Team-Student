\name{pt.mixture}
\alias{pt.mixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CdF of a mixture of t distributions.
}
\description{
CdF of a mixture of t distributions.
}
\usage{
pt.mixture(y, pi, mu, s, df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Value(-s) at which the CDF shall be evaluated (numeric vector)}
  \item{pi}{Mixture proportions (numeric vector of length K > 0; sum(pi)=1, entries non-negativ!)}
  \item{mu}{Component medians (numeric vector of length K)}
  \item{s}{Component sigma's (numeric vector of length K; entries positive!)}
  \item{df}{Component degrees of freedom (numeric vector of length K; entries positive!)}
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





\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}

