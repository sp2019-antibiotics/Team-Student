\name{model.select.t.mixture}
\alias{model.select.t.mixture}

\title{
Select model of a t-mixture distribution}
\description{
Selecting a model
}
\arguments{
  \item{n}{Frequencies (integer vector of length J; non-negative entries)}
 \item{j}{Possible values (numeric vector of length J > 0)}
 \item{Kmax}{Number of maximum components (numeric scalar)}
  \item{atoms}{Values marking 'resistant' observations (numeric vector of length < J; elements must also be in j)}
  \item{draw}{Should results be visualized? (boolean scalar)}
 \item{Ecoff.quantile}{Which quantile should be used for Ecoff? (numeric scalar within (0, 1))}
  \item{pi_cutoff}{Lower bound for group size of 'wild type' (numeric sclara within (0, 1))}
 \item{alpha}{Hyperparameters for MAP estimation (numeric scalar)}
 \item{beta}{Hyperparameters for MAP estimation (numeric scalar)}
  \item{memb.exp}{Clustering parameter (numeric scalar > 1)}
 \item{maxiter}{Maximum number of iterationsn (integer scalar)}
 \item{eps}{Convergence criterion}
  \item{optim.method}{a}
}

\author{
Lisa Allmesberger
Fabian Bergs
Stefan Immler
Michael Kässmann
}
