###################### EM.algorithm.t.mixture(): ################################

#Input:
#j: possible values (numeric vector of length J)
#n: frequencies (integer vector of length J)
#tau_hat: estimated component memberhsips from Initialization (numeric JxK - Matrix)
#alpha, beta: Hyperparameters for MAP estimation (numeric scalar)
#maxiter: maximum number of iterations (integer scalar)
#optim.method: optimizationn algorith to be used in M-Step
#eps: convergence criterion

#Output:
#par: estimated parameter values (numeric 4xK - matrix)

#Uses:
#Qii etc.
#weighted.quantile()

#Used package:
#invgamma
#VGAM
#MatrixStats

#argument have to be checked
EM.algorithm.t.mixture <- function(n, j, tau_hat, alpha, beta, maxiter=1000, optim.method=c("Nelder-Mead", "BFGS"), eps=10^-3)
{
  optim.method = match.arg(optim.method)

  K <- dim(tau_hat)[1]
  N <- sum(n)

  #define parameters
  par <- matrix(0,K,4)
  par[,2] <- rep(weighted.quantile(j,n,0.5),times=K)
  par[,3] <- rep(beta/(1+alpha),times=K)
  par[,4] <- rep(10,times=K)

  par[,3:4] <- log(par[,3:4])

  #i=1 step
  #Computation of pi_k, mu_k, s_k, df_k in the i-th iteration
  par[,1] <- apply(t(tau_hat)*n, 2, sum)/N

  for(k in 1:K)
  {
    par[k,2:4] <- optim(par[k,2:4], Qmax, tau=tau_hat[k,], n=n, j=j, a=alpha, b=beta)$par
  }

  #Computation of tau_hat
  uppers_log = pt((rep(j,each=K)+1/2-par[,2])/exp(par[,3]), exp(par[,4]), log.p=TRUE)
  lowers_log = pt((rep(j,each=K)-1/2-par[,2])/exp(par[,3]), exp(par[,4]), log.p=TRUE)

  numerators_log <- log(par[,1]) + matrix(uppers_log + VGAM::log1mexp(uppers_log-lowers_log), nrow=K)
  denominators_log = matrixStats::colLogSumExps(numerators_log)

  tau_hat_log <- t(t(numerators_log)-denominators_log)

  l = numeric(maxiter)
  l[1] = sum(denominators_log*n)

  #tau_hat <- exp(t(t(numerators_log)-denominators_log))

  #i>=2
  i <- 2
  while(i <= maxiter)
  {
    #Computation of pi_k, mu_k, s_k, df_k   in the i-th iteration
    par[,1] <- apply(t(tau_hat)*n, 2, sum)/N

    for(k in 1:K)
    {
      par[k,2:4] <- optim(par[k,2:4], Qmax, tau=tau_hat[k,], n=n, j=j, a=alpha, b=beta)$par
    }

    #Computation of tau_hat


    uppers_log = pt((rep(j,each=K)+1/2-par[,2])/exp(par[,3]), exp(par[,4]), log.p=TRUE)
    lowers_log = pt((rep(j,each=K)-1/2-par[,2])/exp(par[,3]), exp(par[,4]), log.p=TRUE)

    numerators_log <- log(par[,1]) + matrix(uppers_log + VGAM::log1mexp(uppers_log-lowers_log), nrow=K)
    denominators_log = matrixStats::colLogSumExps(numerators_log)

    tau_hat <- exp(t(t(numerators_log)-denominators_log))

    #check convergence
    l[i] = sum(denominators_log*n)
    if(i>5)
    {
      if(abs(l[i] - l[i-5]) < eps)
      {
        break
      }
    }
    i <- i+1
  }

  if(i > maxiter)
    i = i-1

  par[,3:4] = exp(par[,3:4])
  l = l[1:i]

  result = list(parameters = par, log_likelihood = l, iterations = i)
  return(result)
}


