\name{Qmax}
\alias{Qmax}
\title{
Quantile Max
}
\description{
Quantile Max
}
\usage{
Qmax(par, tau, j, n, a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{Matrix of parameters (nmeric 4xK - matrix); may also be delivered separately.}
  \item{tau}{A value.}
  \item{j}{Possible values (numeric vector of length J > 0)}
  \item{n}{Frequencies (integer vector of length J; non-negative entries)}
  \item{a}{a}
  \item{b}{a}
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
function (par, tau, j, n, a, b)
{
    return(-(Q2(par, tau, j, n) + Q3(par[2], a, b) + Q4(par[3])))
  }
}

