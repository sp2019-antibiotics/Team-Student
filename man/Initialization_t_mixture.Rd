\name{Initialization.t.mixture}
\alias{Initialization.t.mixture}
\title{
Starting Values for EM algorithm
}
\description{
Computes starting values (mixture proportions) for EM algorithm
}
\usage{
Initialization.t.mixture(n, j, K, memb.exp=2, draw=FALSE)
}
\arguments{
\item{n}{Frequencies (integer vector of length J)
}
\item{j}{Possible values (numeric vector of length J)
}
\item{K}{Number of components (numeric scalar)
}
\item{memb.exp}{Clustering parameter (numeric scalar)
}
\item{draw}{Should clustering be drawn? (Boolean scalar)
}
}

\value{
\item{tau}{Class membership values for all classes and all possible values (numeric JxK - matrix)}
}

\author{
Lisa Allmesberger
Fabian Bergs
Stefan Immler
Michael Kässmann
}

\keyword{initialization}




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

