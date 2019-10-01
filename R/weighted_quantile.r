######################## weighted.quantile(): ##################################

#computes weighted quantile

#Input:
#x: data (numeric vector of length N)
#w: weights (numeric vector of length N, entries non-negative)
#probs: which quantile(-s) shall be computed (numeric vector, entries within (0, 1))

#Output:
#q: weighted sample quantile(-s)

weighted.quantile <- function(x, w=rep(1/length(x), times=length(x)), probs=0.5)
{
  #check arguments
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(w))

  N=length(x)
  stopifnot(length(w)==N)

  stopifnot(all(w>=0))
  stopifnot(any(w>0))
  stopifnot(all(0 < probs & probs < 1))

  #prepare arguments
  w = w/sum(w)

  order <- order(x)
  x = x[order]
  w = w[order]

  p = length(probs)
  probs = sort(probs)

  #prepare loop variables
  q = numeric(p)
  W_sum = 0
  i=0
  #calculate quantile
  for(j in 1:p)
  {
    while(W_sum < probs[j])
    {
      i = i+1
      W_sum = W_sum + w[i]
    }

    q[j] = x[i]
  }

  return(q)
}
