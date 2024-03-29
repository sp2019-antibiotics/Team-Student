\name{ungroup.data}
\alias{ungroup.data}
\title{
Ungroup data
}

\description{
A function that ungroups data (if N>700, the resulting vector is not the actual sample, but a close representation of it with length 500)
}

\usage{
ungroup.data(n, j)
}

\arguments{
  \item{n}{possible values (numeric vector of length J > 0)}
  \item{j}{frequencies (integer vector of length J; non-negative entries)}
}

\value{
\item{y}{ungrouped data vector (length: min(500, N))}
}

\author{
Lisa Allmesberger
Fabian Bergs
Stefan Immler
Michael Kässmann
}



