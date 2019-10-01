\name{Q2}
\alias{Q2}
\title{
Quantile2
}
\description{
Quantile 2
}
\usage{
Q2(par, tau, j, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{Matrix of parameters (nmeric 4xK - matrix); may also be delivered separately.}
  \item{tau}{A value}
   \item{j}{Possible values (numeric vector of length J > 0)}
   \item{n}{Frequencies (integer vector of length J; non-negative entries)}
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
function (par, tau, j, n)
{
    upper_log = pt((j + 1/2 - par[1])/exp(par[2]), exp(par[3]),
        log.p = TRUE)
    lower_log = pt((j - 1/2 - par[1])/exp(par[2]), exp(par[3]),
        log.p = TRUE)
    return(sum((upper_log + VGAM::log1mexp(upper_log - lower_log)) *
        tau * n))
  }
}

