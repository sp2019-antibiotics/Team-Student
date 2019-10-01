\name{Q}
\alias{Q}
\title{
Quantile
}
\description{
Calculates the quantile
}
\usage{
Q(par, j, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{Matrix of parameters (nmeric 4xK - matrix); may also be delivered separately.}
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
function (par, j, n)
{
    K = dim(par)[1]
    uppers_log = pt((rep(j, each = K) + 1/2 - par[, 2])/exp(par[,
        3]), exp(par[, 4]), log.p = TRUE)
    lowers_log = pt((rep(j, each = K) - 1/2 - par[, 2])/exp(par[,
        3]), exp(par[, 4]), log.p = TRUE)
    numerators_log <- log(par[, 1]) + matrix(uppers_log + VGAM::log1mexp(uppers_log -
        lowers_log), nrow = K)
    denominators_log = matrixStats::colLogSumExps(numerators_log)
    tau <- exp(t(t(numerators_log) - denominators_log))
    return(sum(n * tau * log(par[, 1])) + Q2(par[, 2:4], tau,
        j, n))
  }
}

