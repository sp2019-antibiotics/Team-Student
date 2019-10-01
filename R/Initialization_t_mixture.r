######################## Initialization.t.mixture(): #################

#computes starting values (mixture proportions) for EM algorithm

#Input:
#j: possible values (numeric vector of length J)
#n: frequencies (integer vector of length J)
#K: number of components (numeric scalar)
#memb.exp: clustering parameter (numeric scalar)
#draw: Should clustering be drawn? (Boolean scalar):

#Output:
#tau: class membership values for all classes and all possible values (numeric JxK - matrix)

#Uses:
#weighted.quantile()

#Used package(-s):
#cluster

Initialization.t.mixture<-function(n, j, K, memb.exp=2, draw=FALSE)
{
  #check arguments
  stopifnot(is.numeric(n))
  stopifnot(is.vector(n))

  n = as.integer(n)
  stopifnot(all(n>=0))
  N = sum(n)
  stopifnot(N>0)

  stopifnot(is.numeric(j))
  stopifnot(is.vector(j))
  stopifnot(length(j) == length(n))

  j = j[n>0]
  n = n[n>0]

  J = length(j)
  stopifnot(J>0)

  ord = order(j)
  j = j[ord]
  n = n[ord]

  stopifnot(is.numeric(K))
  stopifnot(length(K)==1)
  K <- as.integer(K)
  stopifnot(K >= 1)

  stopifnot(is.numeric(memb.exp))
  stopifnot(length(memb.exp)==1)
  stopifnot(memb.exp > 1)

  stopifnot(is.logical(draw))

  if(N > 700)
  {#draw "sample"
    y = weighted.quantile(j, n, probs=(1:499)/500)
    n_sample = table(y)
    j_sample = as.numeric(names(n_sample))
  }
  else #ungroup data (necessary for clustering)
  {
    y = numeric(N)
    a = 1
    b = 0
    for(k in 1:J)
    {
      b = b + n[k]
      y[a:b] = j[k]
      a = b+1
    }
  }

  #perform  fuzzy clustering
  clustering = cluster::fanny(y, K, memb.exp=memb.exp, metric = "euclidean", cluster.only = TRUE)
  clustervector = clustering$clustering
  tau_cluster = clustering$membership

  #draw clustering (if requested)
  if(draw==TRUE)
  {
    draw.tau = cbind(y, tau_cluster, clustervector)
    draw.tau = draw.tau[order(y),]
    freq = table(clustervector)

    print("Best Hard Clustering and Membership functions")
    plot(x=y[clustervector==1], y=rep(1, times=freq[1]),
         col=1, xlim=c(min(y), max(y)), ylim=c(0, 1.1), main="Clustering", xlab="y", ylab="Membership")
    lines(draw.tau[,1], draw.tau[,2], col=1)
    if(K > 1)
    {
      for(i in 2:K)
      {
        points(x=y[clustervector==i], y=rep(1, times=freq[i]), col=i)
        lines(draw.tau[,1], draw.tau[,i+1], col=i)
      }
    }

  }

  if(N<=700)
  {
    tau = tau_cluster[cumsum(n),]
    return(tau)
  }
  else #search closest point that got in random sample
  {
    tau = matrix(0, ncol=K, nrow=J)

    l=1
    while(j[l] <= j_sample[1])
    {
      tau[l,] = tau_cluster[1,]
      l = l+1
    }

    current = 1
    while(j[l] <= tail(j_sample, 1))
    {
      while(j_sample[current] < j[l])
        current = current + 1

      if(j_sample[current] - j[l] < j[l] - j_sample[current-1])
      {
        tau[l,] = tau_cluster[current,]
      }
      else
      {
        tau[l,] = tau_cluster[current-1,]
      }
      l = l+1

      if(l > J)
        break
    }
    length(j)
    while(l <= J)
    {
      tau[l,] = as.vector(tail(tau_cluster, 1))
      l = l+1
    }
    return(tau)
  }
}
