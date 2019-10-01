model.select.t.mixture <- function(n, j, Kmax, criterion=c("BIC", "AIC"), atoms=NULL, draw=FALSE, Ecoff.quantile=0.01, pi_cutoff=0.2, alpha=0.05, beta=NULL, memb.exp=2, maxiter=1000, eps=10^-3, optim.method=c("Nelder-Mead", "BFGS"))
{
  #check arguments
  criterion = match.arg(criterion)

  stopifnot(is.numeric(Kmax))
  stopifnot(length(Kmax)==1)
  Kmax <- as.integer(Kmax)
  stopifnot(Kmax > 0)

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

  Ntilde = sum(n[!(j %in% atoms)])
  if(Kmax >= Ntilde/2)
  {
    warning("number of clusters needs to be <N/2, but Kmax >= N/2. Hence, Kmax is set to N/2 - 1.")
    Kmax = Ntilde/2 - 1
  }


  stopifnot(is.numeric(Kmax))
  stopifnot(length(Kmax)==1)
  k <- as.integer(Kmax)
  stopifnot(Kmax>=1)

  stopifnot(is.logical(draw))

  stopifnot(is.numeric(Ecoff.quantile))
  stopifnot(length(Ecoff.quantile) == 1)
  stopifnot(0<Ecoff.quantile & Ecoff.quantile<1)

  stopifnot(is.numeric(pi_cutoff))
  stopifnot(length(pi_cutoff) == 1)
  stopifnot(0<=pi_cutoff & pi_cutoff<1)

  stopifnot(is.numeric(alpha))
  stopifnot(length(alpha) == 1)
  stopifnot(0<alpha)

  stopifnot(is.numeric(alpha))
  stopifnot(length(alpha) == 1)
  stopifnot(0<alpha)



  stopifnot(is.numeric(memb.exp))
  stopifnot(length(memb.exp) == 1)
  stopifnot(1<memb.exp)

  stopifnot(is.numeric(eps))
  stopifnot(length(eps) == 1)
  stopifnot(0<eps)

  stopifnot(is.numeric(maxiter))
  stopifnot(length(maxiter) == 1)
  stopifnot(0<maxiter)


  optim.method = match.arg(optim.method)

  print("Current number of components:")
  optim.method = match.arg(optim.method)
  result = total.function(n, j, 1, atoms, draw, Ecoff.quantile, pi_cutoff, alpha, beta, memb.exp, maxiter, eps, optim.method)
  BIC_max  = result$BIC
  AIC_max = result$AIC

  if(Kmax == 1)
  {
    return(result)
  }

  if(criterion=="BIC")
  {
    for(k in 2:Kmax)
    {
      print(k)
      result_new = total.function(n, j, k, atoms, draw, Ecoff.quantile, pi_cutoff, alpha, beta, memb.exp, maxiter, eps, optim.method)
      BIC_new = result_new$BIC
      if(BIC_new > BIC_max)
      {
        result = result_new
        BIC_max = BIC_new
      }
    }
  }
  else
  {
    for(k in 2:Kmax)
    {
      print(k)
      result_new = total.function(n, j, k, atoms, draw, Ecoff.quantile, pi_cutoff, alpha, beta, memb.exp, maxiter, eps, optim.method)
      AIC_new = result_new$AIC

      if(AIC_new > AIC_max)
      {
        result = result_new
        AIC_max = AIC_new
      }
    }
  }

  return(result)
}
