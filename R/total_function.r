####################### total.function(): ########################

#Input:
#j: possible values (numeric vector of length J > 0)
#n: frequencies (integer vector of length J; non-negative entries)
#K: number of components (numeric scalar)
#atoms: values marking 'resistant' observations (numeric vector of length < J; elements must also be in j)
#draw: Should results be visualized? (boolean scalar)
#Ecoff.quantile: Which quantile should be used for Ecoff? (numeric scalar within (0, 1))
#pi_cutoff: lower bound for group size of 'wild type' (numeric sclara within (0, 1))
#memb.exp: clustering parameter (numeric scalar > 1)
#alpha, beta: Hyperparameters for MAP estimation (numeric scalar)
#maxiter: maximum number of iterationsn (integer scalar)
#eps: convergence <=> relative change in likelihood smaller than eps
#optim.method: optimization method to be used (one of "Nelder-Mead" (default), "BFGS"; see documentation of optim())

#Output:
#par: estimated parameter values (numeric 4xK - Matrix)
#Ecoff: calculated Ecoff (numeric scalar)
#AIC, BIC: AIC and BIC of fitted model
#log_likelihood: Value of log-likelihood function in each EM iteration
#EM_iterations: number of EM iterations

#Uses:
#Initialization.t.mixture()
#EM.algorithm.e.mixture()
#weigthed.quantile()
#ECOFF.t.mixture()
#visualization.t.mixture()

#packages used:
#invgamma
#cluster

#needs to revised yet (maybe also its description)
total.function<-function(n, j, K, atoms=NULL, draw=FALSE, Ecoff.quantile=0.01, pi_cutoff=0.2, alpha=0.05, beta=NULL, memb.exp=2, maxiter=1000, eps=10^-3, optim.method=c("Nelder-Mead", "BFGS"))
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
  stopifnot(all(atoms %in% j))
  stopifnot(length(j) == length(n))

  j = j[n>0]
  n = n[n>0]

  J = length(j)
  stopifnot(J>0)

  stopifnot(is.numeric(K))
  stopifnot(length(K)==1)
  k <- as.integer(K)
  stopifnot(K>=1)

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

  n_with_resistant = n
  j_with_resistant = j
  pi_resistant = 0

  #remove "resistant" observations
  if(!is.null(atoms))
  {
    stopifnot(is.vector(atoms))
    stopifnot(is.numeric(atoms))

    atoms = atoms[atoms %in% j]
    if(length(atoms) == 0)
    {
      atoms = NULL
      print("No resistant observations in sample.")
    }
    else
    {
      r = sapply(atoms, function(a){which(j==a)})
      N_resistant = sapply(r, function(x){sum(n[x])})
      pi_resistant = N_resistant/N
      N = N - sum(N_resistant)

      n = n[-r]
      j = j[-r]
    }

  }

  #calculate starting values
  tau_hat = Initialization.t.mixture(n=n, j=j, K=K, memb.exp=memb.exp, draw=draw)

  #define beta
  if(is.null(beta))
  {
    beta = (1+alpha)*0.36*(weighted.quantile(x=j, w=n, probs=0.9) - weighted.quantile(x=j, w=n, probs=0.1))/K
  }
  else
  {
    stopifnot(is.numeric(beta))
    stopifnot(length(beta) == 1)
    stopifnot(0<beta)
  }


  #perform EM algorithm
  result_EM = EM.algorithm.t.mixture(n, j, t(tau_hat), alpha, beta, maxiter=maxiter, optim.method=optim.method, eps=eps)
  final.values = result_EM$parameters
  final.values[,1] = final.values[,1]/sum(final.values[,1])

  #sort by mu
  if(K>1)
  {
    final.values = final.values[order(final.values[, 2]),]
  }


  #calculate ECOFF
  ECOFF = ECOFF.t.mixture(par=final.values, quantile=Ecoff.quantile, pi_cutoff=pi_cutoff)
  if(ECOFF[2] > max(j) | ECOFF[2] < min(j))
    warning("Calculated ECOFF is out of range of j")

  #calculate BIC & AIC
  BIC = BIC.t.mixture(n, j, final.values[,1], final.values[,2], final.values[,3], final.values[,4])
  AIC = AIC.t.mixture(n, j, final.values[,1], final.values[,2], final.values[,3], final.values[,4])

  #names
  rownames(final.values) <- paste("group", 1:K, sep="")
  colnames(final.values) <- c("pi", "mu", "s", "df")

  #deal with "resistant" observations
  if(!is.null(atoms))
  {

    atoms.values = matrix(c(pi_resistant, atoms, rep(NA, times=length(atoms)), rep(NA, times=length(atoms))), ncol=4)
    final.with_resistant = rbind(atoms.values, final.values)
    final.with_resistant[,1] = final.with_resistant[, 1]/sum(final.with_resistant[, 1])

    rownames(final.with_resistant)[1:length(atoms)] <- paste("resistant", 1:length(atoms), sep="")
  }

  #draw density (if requested)
  if(draw==TRUE)
  {
    visualization.t.mixture(n_with_resistant, j_with_resistant,
                            final.values[,1], final.values[,2], final.values[,3], final.values[,4],
                            ecoff=ECOFF[2], pi_resistant=pi_resistant)
  }

  #return estimates
  if(!is.null(atoms))
    return(list(parameters = final.with_resistant, ECOFF = ECOFF, AIC = AIC, BIC = BIC,
                log_likelihood = result_EM$log_likelihood, EM_iterations = result_EM$iterations))
  else
    return(list(parameters = final.values, ECOFF = ECOFF, AIC = AIC, BIC = BIC,
                log_likelihood = result_EM$log_likelihood, EM_iterations = result_EM$iterations))
}


