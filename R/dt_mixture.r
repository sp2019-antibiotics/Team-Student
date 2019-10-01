####################### dt.mixture(): ################################################

#PDF of a mixture of t distributions

#Input:
#y: value(-s) at which the PDF shall be evaluated (numeric vector)
#pi: mixture proportions (numeric vector of length K > 0; sum(pi) = 1, entries non-negative)
#mu: component medians (numeric vector of lentgh K)
#s: component sigma's (numeric vector of lentgh K; entries positive)
#df: component degrees of freedom (numeric vector of lentgh K; entries positive)

#Output:
#result: f(y)

#arguments are to be checked
dt.mixture <- function(y, pi, mu, s, df)
{
  stopifnot(is.numeric(y))

  stopifnot(is.numeric(pi))
  stopifnot(is.vector(pi))

  stopifnot(is.numeric(mu))
  stopifnot(is.vector(mu))

  stopifnot(is.numeric(s))
  stopifnot(is.vector(s))

  stopifnot(is.numeric(df))
  stopifnot(is.vector(df))

  K = length(pi)
  stopifnot(K > 0)
  stopifnot(length(mu)==K & length(s)==K & length(df)==K)

  stopifnot(all(pi>=0))

  stopifnot(all(s>0))
  stopifnot(all(df>0))

  by_components = sapply(1:K, function(i) {
    (1/s[i]) * dt((y - mu[i])/s[i], df[i])})

  result = as.numeric(by_components %*% pi)
  return(result)
}
