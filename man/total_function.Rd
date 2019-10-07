\name{total.function}
\alias{total.function}
\title{
Total Function
}
\description{
The main function for the analysis
}

\usage{
total.function(n, j, K, atoms=NULL, draw=FALSE, Ecoff.quantile=0.01, pi_cutoff=0.2, alpha=0.05, beta=NULL, memb.exp=2, maxiter=1000, eps=10^-3, optim.method=c("Nelder-Mead", "BFGS"))
}
\arguments{
 \item{n}{Frequencies (integer vector of length J; non-negative entries)}
 \item{j}{Possible values (numeric vector of length J > 0)}
 \item{K}{Number of components (numeric scalar)}
 \item{atoms}{Values marking 'resistant' observations (numeric vector of length < J; elements must also be in j)}
 \item{draw}{Should results be visualized? (boolean scalar)}
 \item{Ecoff.quantile}{Which quantile should be used for Ecoff? (numeric scalar within (0, 1))}
 \item{pi_cutoff}{Lower bound for group size of 'wild type' (numeric sclara within (0, 1))}
 \item{alpha}{Hyperparameters for MAP estimation (numeric scalar)}
 \item{beta}{Hyperparameters for MAP estimation (numeric scalar)}
 \item{memb.exp}{Clustering parameter (numeric scalar > 1)}
 \item{maxiter}{Maximum number of iterationsn (integer scalar)}
 \item{eps}{Convergence criterion; relative change in likelihood smaller than eps}
 \item{optim.method}{optimization method to be used (one of "Nelder-Mead" (default), "BFGS"; see documentation of optim())}
}

\author{
Lisa Allmesberger
Fabian Bergs
Stefan Immler
Michael Kässmann
}

\value{
An object of class `\code{list}' with elements:
  \item{parameters}{Kx4 matrix containing the calculated parameters.}
  \item{ECOFF}{1x2 matrix containing the group index and ECOFF.}
  \item{AIC}{AIC value.}
  \item{BIC}{BIC value.}
  \item{log_likelihood}{1xM matrix containing the log likelihood value for each EM iteration. (M number of EM iterations)}
  \item{EM_iterations}{Number of EM iterations.}
}

\keyword{total funcion}



