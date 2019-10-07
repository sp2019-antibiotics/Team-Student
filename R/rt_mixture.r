######################### rt.mixture() #############################

#Simulates N observations from a mixture of t distributions

#Input:
#N: Number of realizations (integer scalar, positive)
#pi: mixture proportions (numeric vector of length K > 0; sum(pi) = 1, entries non-negative)
#mu: component medians (numeric vector of lentgh K)
#s: component sigma's (numeric vector of lentgh K; entries positive)
#df: component degrees of freedom (numeric vector of lentgh K; entries positive)
#round: Should number be rounded (booelan scalar)

#Output:
#y: N accordingly distributed realizations (numeric vector of length N)

rt.mixture<-function(N, pi, mu, s, df, round=FALSE)
{
  #prepare/check arguments
  stopifnot(is.numeric(N))
  stopifnot(length(N)==1)
  N <- as.integer(N)
  K = length(pi)
  stopifnot(K>=1)

  stopifnot(is.numeric(pi))
  stopifnot(is.numeric(mu))
  stopifnot(is.numeric(s))
  stopifnot(is.numeric(df))
  stopifnot(is.logical(round))


  stopifnot(is.vector(pi))
  stopifnot(is.vector(mu))
  stopifnot(is.vector(s))
  stopifnot(is.vector(df))
  stopifnot(length(round)==1)

  K = length(pi)
  stopifnot(K > 0)
  stopifnot(length(mu)==K & length(s)==K & length(df)==K)

  stopifnot(all(s>0))
  stopifnot(all(df>0))
  stopifnot(all(pi>=0))

  pi = pi/sum(pi)
  param = list("pi"=pi, "mu"=mu, "s"=s, "df"=df)
  lengths = sapply(param, length)
  K = lengths[1]
  if(any(lengths != K))
  {
    warning("Parameter vectors have different lengths. The number of classes is taken to be the minimal length. Other parameters are deleted.")
    K = min(lengths)
    stopifnot(K>0)
    param = sapply(param, function(x){x[1:K]})
  }

  #generate class memberships
  z<-sample(1:K, size=N, prob=param$pi, replace = TRUE)

  #count class sizes
  freq = sapply(1:K, function(i){sum(z==i)})
  cum_freq = cumsum(freq)

  #generate observations
  y<-numeric(N)
  for(i in 1:K)
  {
    y[(cum_freq[i]-freq[i]+1):(cum_freq[i])] =
      rt(freq[i], df[i]) * s[i] + mu[i]
  }

  #round (if requested)
  if(round==TRUE)
  {
    y<-round(y, 0)
  }

  #create random permutation of values
  y = y[sample(1:N, N, replace = FALSE)]

  return(y)
}

