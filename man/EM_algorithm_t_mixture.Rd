\name{EM.algorithm.t.mixture}
\alias{EM.algorithm.t.mixture}
\title{
EM algorithm for t-mixture
}
\description{
EM algorithm for t-mixture
}
%usage too large
%\usage{
%EM.algorithm.t.mixture(n,j,tau_hat,alpha,beta,maxiter=1000,optim.method=c("Nelder%-Mead","BFGS"),eps=10^-3)
%}
\arguments{
\item{n}{Frequencies (integer vector of length J)
}
\item{j}{Possible values (numeric vector of length J)
}
\item{tau_hat}{estimated component memberhsips from Initialization (numeric JxK - Matrix)
}
\item{alpha}{Hyperparameters for MAP estimation (numeric scalar)
}
\item{beta}{Hyperparameters for MAP estimation (numeric scalar)
}
\item{maxiter}{Maximum number of iterations (integer scalar)
}
\item{eps}{Convergence criterion
}
\item{optim.method}{Optimization algorithm to be used in M-Step}
}


\value{
\item{par}{Estimated parameter values (numeric 4xK - matrix)}
}

\author{
Lisa Allmesberger
Fabian Bergs
Stefan Immler
Michael Kässmann
}




