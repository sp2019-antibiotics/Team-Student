\name{dt.mixture}
\alias{dt.mixture}
\title{
Dt Mixture Function
}
\description{
PDF of a mixture of t distributions
}
\usage{
dt.mixture(y, pi, mu, s, df)
}
\arguments{
  \item{y}{Value(-s) at which the PDF shall be evaluated (numeric vector)}
  \item{pi}{Mixture proportions (numeric vector of length K > 0; sum(pi) = 1, entries non-negative)}
  \item{mu}{Component medians (numeric vector of lentgh K)}
  \item{s}{Component sigma's (numeric vector of lentgh K; entries positive)}
  \item{df}{Component degrees of freedom (numeric vector of lentgh K; entries positive)}
}

\value{
\item{f(y)}{Result}
}

\author{
Lisa Allmesberger
Fabian Bergs
Stefan Immler
Michael Kässmann
}

\seealso{
qt.mixture, pt.mixture, rt.mixture
}

\keyword{dt.mixture}




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

