\name{bisection.method}
\alias{bisection.method}

\title{
Bisection method
}

\description{
Solves equation f(x) = y numerically for x
}

\usage{
bisection.method(y, f, limits, tol=10*.Machine$double.eps, ...)
}

\arguments{
  \item{y}{Value that f should attain (numeric scalar)}
  \item{f}{Function to be inverted (function; should be strictly monotonous)}
  \item{limits}{Search interval (vector of length 2; must contain solution!)}
  \item{tol}{Value}
  \item{...}{Argument to be passed to methods.}
}

\value{
\item{x}{Such that f(x) = y}
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

