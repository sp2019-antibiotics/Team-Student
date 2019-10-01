\name{Q3}
\alias{Q3}
\title{
Quantile3
}
\description{
Quantile 3
}
\usage{
Q3(s, a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{Logarithm of standard deviation}
  \item{a}{Parameters of invers Gamma distribution}
  \item{b}{Parameters of invers Gamma distribution}
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
function (s, a, b)
{
    return(dinvgamma(exp(s)^2, shape = a, scale = b, log = TRUE))
  }
}

