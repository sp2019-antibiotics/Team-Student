######################################
######################################
### 'Main' functions for analysis

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
      #print(c(l, j[l], j_sample[current]))
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
        #print(c(l, j[l], j_sample[current]))
      }
      else
      {
        tau[l,] = tau_cluster[current-1,]
        #print(c(l, j[l], j_sample[current]))
      }
      l = l+1

      if(l > J)
        break
    }
    length(j)
    while(l <= J)
    {
      tau[l,] = as.vector(tail(tau_cluster, 1))
      #print(c(l, j[l], j_sample[current]))
      l = l+1
    }
    return(tau)
  }
}


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
  #print(par[,1])


  for(k in 1:K)
  {

    #print(paste("k =", k))
    par[k,2:4] <- optim(par[k,2:4], Qmax, tau=tau_hat[k,], n=n, j=j, a=alpha, b=beta)$par
  }

  #Computation of tau_hat
  #numerators <- par[,1]*matrix(pt((rep(j,each=K)+1/2-par[,2])/par[,3], par[,4])-
  #                             pt((rep(j,each=K)-1/2-par[,2])/par[,3], par[,4]), nrow=K)

  #numerators <- par[,1]*matrix(pt((rep(j,each=K)+1/2-par[,2])/exp(par[,3], exp(par[,4]))-
  #                              pt((rep(j,each=K)-1/2-par[,2])/exp(par[,3]), exp(par[,4])), nrow=K)
  #denominators <- apply(numerators, 2, sum)

  #tau_hat <- t(t(numerators)/denominators)

  #c(print(length(J)*K/length(par[,2])), print(length(J)*K/length(par[,3])), print(length(J)*K/length(par[,4])))

  uppers_log = pt((rep(j,each=K)+1/2-par[,2])/exp(par[,3]), exp(par[,4]), log.p=TRUE)
  lowers_log = pt((rep(j,each=K)-1/2-par[,2])/exp(par[,3]), exp(par[,4]), log.p=TRUE)

  numerators_log <- log(par[,1]) + matrix(uppers_log + VGAM::log1mexp(uppers_log-lowers_log), nrow=K)

  #denominators_log <- apply(numerators_log, 2, matrixStats::logSumExp(numerators_log))
  denominators_log = matrixStats::colLogSumExps(numerators_log)

  tau_hat_log <- t(t(numerators_log)-denominators_log)

  l = numeric(maxiter)
  l[1] = sum(denominators_log*n)

  #tau_hat <- exp(t(t(numerators_log)-denominators_log))

  #colSums(tau_hat)
  #summary(as.vector(exp(tau_hat_log)))
  #colSums(exp(tau_hat_log))
  #summary(as.vector(tau_hat))

  #i>=2
  i <- 2
  while(i <= maxiter)
  {
    #print(paste("i = ", i))

    #regularized version of l (i-1 step)
    #lR <- sum(t(n*t(log(numerators))))+sum(Q3(par[,3], alpha, beta))

    #Computation of pi_k, mu_k, s_k, df_k   in the i-th iteration
    par[,1] <- apply(t(tau_hat)*n, 2, sum)/N

    for(k in 1:K)
    {
      #print(paste("k =", k))
      par[k,2:4] <- optim(par[k,2:4], Qmax, tau=tau_hat[k,], n=n, j=j, a=alpha, b=beta)$par
    }

    #Computation of tau
    #numerators <- par[,1]*matrix(pt((rep(j,each=K)+1/2-par[,2])/par[,3], par[,4])-
    #                              pt((rep(j,each=K)-1/2-par[,2])/par[,3], par[,4]), nrow=K)
    #numerators <- par[,1]*matrix(pt((rep(j,each=K)+1/2-par[,2])/exp(par[,3]), exp(par[,4]))-
    #                              pt((rep(j,each=K)-1/2-par[,2])/exp(par[,3]), exp(par[,4])), nrow=K)


    #denominators <- apply(numerators, 2, sum)

    #tau_hat <- t(t(numerators)/denominators)

    #c(print(length(J)*K/length(par[,2])), print(length(J)*K/length(par[,3])), print(length(J)*K/length(par[,4])))


    uppers_log = pt((rep(j,each=K)+1/2-par[,2])/exp(par[,3]), exp(par[,4]), log.p=TRUE)
    lowers_log = pt((rep(j,each=K)-1/2-par[,2])/exp(par[,3]), exp(par[,4]), log.p=TRUE)

    numerators_log <- log(par[,1]) + matrix(uppers_log + VGAM::log1mexp(uppers_log-lowers_log), nrow=K)
    #denominators_log <- apply(numerators_log, 2, matrixStats::logSumExp(numerators_log))

    denominators_log = matrixStats::colLogSumExps(numerators_log)

    #tau_hat_log <- t(t(numerators_log)-denominators_log)
    tau_hat <- exp(t(t(numerators_log)-denominators_log))

    #regularized version of l (i step)
    #lRnew <- sum(t(n*t(log(numerators))))+sum(Q3(par[,3], alpha, beta))
    l[i] = sum(denominators_log*n)
    if(i>5)
    {
      #print(abs(l[i] - l[i-5]) < eps)
      if(abs(l[i] - l[i-5]) < eps)
      {
        break
      }
    }


    #convergence criterion
    #if(abs(lR - lRnew) < eps)
    # i <- 10000

    i <- i+1
  }
  par[,3:4] = exp(par[,3:4])
  l = l[1:i]

  result = list(parameters = par, log_likelihood = l, iterations = i)
  return(result)
}



####################### visualization.t.mixture(): ###########################################

#Visualizes results of analysis

#Input:
#y: data (numeric vector)
#pi: mixture proportions (numeric vector of length K > 0; sum(pi) = 1, entries non-negative)
#mu: component medians (numeric vector of lentgh K)
#s: component sigma's (numeric vector of lentgh K; entries positive)
#df: component degrees of freedom (numeric vector of lentgh K; entries positive)
#ECOFF: Should ECOFF be plotted as well? (booelan scalar)

#Output: None (only plots)

#Uses:
#dt.mixture()
#ECOFF.t.mixture()

#hist() may not work well yet
#grouped data instead of y
#arguments are to be checked

visualization.t.mixture <- function(n, j, pi, mu, s, df, ecoff=NULL, kernel.density=TRUE, hist=TRUE, new.plot=hist, pi_resistant=0)
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
  stopifnot(sum(pi)==1)
  stopifnot(all(s>0))
  stopifnot(all(df>0))

  stopifnot(is.logical(kernel.density))
  stopifnot(is.logical(hist))
  stopifnot(is.logical(new.plot))

  stopifnot(is.numeric(pi_resistant))
  stopifnot(length(pi_resistant)==1)
  stopifnot(0<=pi_resistant && pi_resistant<=1)

  #ungroup and plot data
  y = ungroup.data(n, j)

  max.data <- max(y)
  min.data <- min(y)
  hist.step <- (max.data - min.data)

  t = seq(from = min.data, to = max.data, by = 0.1)
  ft = dt.mixture(t, pi, mu, s, df)
  ft = ft/(1+pi_resistant)

  dens = density(y)

  ylim2 = max(ft, dens$y)

  if (hist.step > 50)
    hist.step <- hist.step
  h=hist(
    y,
    breaks = seq(min.data, max.data, l = hist.step),
    main="Density Estimates",
    probability = TRUE,
    xlim = c(min.data, max.data),
    ylim = c(0, ylim2),
    col="grey"
    ,xaxt = 'n'
  )

  if(max(h$density) > ylim2)
  {
    ylim2 = max(h$density)
    h=hist(
      y,
      breaks = seq(min.data, max.data, l = hist.step),
      main="Density Estimates",
      probability = TRUE,
      xlim = c(min.data, max.data),
      ylim = c(0, max(h$density)),
      col="grey"
      ,xaxt = 'n'
    )

  }
  axis(side = 1,
       at = seq(min.data, max.data, (max.data - min.data) / hist.step))
  lines(dens, col = "red", lwd = 2)
  lines(t, ft, col = "blue", lwd = 2)

  if(!is.null(ecoff))
  {
    stopifnot(is.numeric(ecoff))
    stopifnot(length(ecoff)==1)

    if(ecoff < min.data | ecoff < max.data)
      warning("Forwarded ECOFF is out of range of j and will be 'left' or 'right' to the plot.")

    abline(v=ecoff, lty="dashed", lwd=2)
    legend("topleft", col = c("red", "blue", "black"),
           legend = c("Kernel Density", "t-Density", "ECOFF"), lty = c(1, 1, 2), cex = 0.75, lwd = 2)
  }
  else
  {
    legend("topleft", col = c("red", "blue"),
           legend = c("Kernel Density", "t-Density"), lty = 1, cex = 0.75, lwd = 2)
  }
}


####################### ECOFF.t.mixture(): #####################

#computes ECOFF for given parameter values

#Input:
#par: Matrix of parameters (numeric 4xK - matrix); may also be delivered separately:
#pi: mixture proportions (numeric vector of length K > 0; sum(pi) = 1, entries non-negative)
#mu: component medians (numeric vector of lentgh K)
#s: component sigma's (numeric vector of lentgh K; entries positive)
#df: component degrees of freedom (numeric vector of lentgh K; entries positive)
#pi_resistant: proportion of resistant observations (numerical scalar within (0, 1))
#quantile: cut-off value (numeric scalar within (0, 1))
#pi_cutoff: lower bound for groupsize of 'wild type' (numeric scalar within (0, 1))

#Ouput:
#result: group index (integer scalar within [1, K]), ECOFF (numeric scalar)

#Uses:
#qt.mixture()

#arguments are to be checked
ECOFF.t.mixture <- function(par=NULL, pi=NULL, mu=NULL, s=NULL, df=NULL,
                            quantile=0.01, pi_cutoff = 0.2)
{

  stopifnot(is.numeric(quantile))
  stopifnot(is.vector(quantile))
  stopifnot(0<quantile & quantile<1)

  stopifnot(is.numeric(pi_cutoff))
  stopifnot(length(pi_cutoff) == 1)
  stopifnot(0<=pi_cutoff & pi_cutoff<1)

  if(is.null(par))
  {
    if(is.null(pi) | is.null(mu) | is.null(s) | is.null(df))
      stop("No parameters forwarded")

    K = length(pi)
    stopifnot(K>0)
    parList = list(mu, s, df)
    lengths = sapply(parList, length)

    stopifnot(all(lengths == K))

    par = cbind(pi, mu, s, df)
  }

  stopifnot(is.numeric(par))
  stopifnot(all(par[,1]>0))
  stopifnot(sum(par[,1])==1)
  stopifnot(all(par[,3]>0))
  stopifnot(all(par[,4]>0))


  if(dim(par)[1] == 1)
  {
    return(c(1, qt(quantile, par[,4])*par[,3] + par[2]))
  }

  par = par[order(par[,2], decreasing=TRUE),]

  i=0
  sum_pi = 0
  while(sum_pi < pi_cutoff)
  {
    i = i+1
    sum_pi = sum_pi + par[i,1]
  }

  if(i==length(pi))
    warning("There is no non-trivial collection of groups that exceed 'pi_cutoff'.
            Maybe try a lower value for 'pi_cutoff'")

  par_wild = matrix(par[1:i,], ncol=4)
  par_wild[,1] = par_wild[,1]/sum(par_wild[,1])

  result = c(i, qt.mixture(p=quantile, pi=par_wild[,1],
                           mu=par_wild[,2], s=par_wild[,3], df=par_wild[,4]))
  names(result) = c("Group Index", "ECOFF")
  return(result)
}

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

#Output:
#par: estimated parameter values (numeric 4xK - Matrix)
#Ecoff: calculated Ecoff (numeric scalar)

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
  #else #pi_resistant is needed for ECOFF later on
  #{
  #  pi_resistant = 0
  #}

  #calculate starting values
  #starting.values = Starting.Values.Mixtures.t(y, K, starting.estimation, cluster.method, draw=draw)
  tau_hat = Initialization.t.mixture(n=n, j=j, K=K, memb.exp=memb.exp, draw=draw)

  #define beta
  if(is.null(beta))
  {
    #centers = j[apply(tau_hat, 1, which.max)]
    #range_center = range(centers)
    #IQR = weighted.quantile(x=j, w=n, p=0.75) - weighted.quantile(x=j, w=n, p=0.25)
    #beta = (1+alpha)*0.6*(weighted.quantile(x=j, w=n, p=0.75) - weighted.quantile(x=j, w=n, p=0.25))/K
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
  final.values = final.values/sum(final.values)
  #sort by mu
  if(K>1)
  {
    final.values = final.values[order(final.values[, 2]),]
  }


  #calculate ECOFFS
  #ECOFF1 = unname(qt(Ecoff.quantile, final.values[1, 4]) * final.values[1, 3] + final.values[1, 2])
  #ECOFF2 = unname(qt(Ecoff.quantile, final.values[K, 4]) * final.values[K, 3] + final.values[K, 2])

  #ECOFF = ECOFF.t.mixture(par=final.values, pi_resistant=pi_resistant,
  #quantile=Ecoff.quantile, pi_cutoff=pi_cutoff)
  #print(final.values)
  ECOFF = ECOFF.t.mixture(par=final.values, quantile=Ecoff.quantile, pi_cutoff=pi_cutoff)
  if(ECOFF[2] > max(j) | ECOFF[2] < min(j))
    warning("Calculated ECOFF is out of range of j")

  #calculate BIC
  BIC = BIC_t_mixture(n, j, final.values[,1], final.values[,2], final.values[,3], final.values[,4])
  AIC = AIC_t_mixture(n, j, final.values[,1], final.values[,2], final.values[,3], final.values[,4])

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


model.select.t.mixture <- function(n, j, Kmax, atoms=NULL, draw=FALSE, Ecoff.quantile=0.01, pi_cutoff=0.2, alpha=0.05, beta=NULL, memb.exp=2, maxiter=1000, eps=10^-3, optim.method=c("Nelder-Mead", "BFGS"))
{
  #check arguments
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






  #print("New Bacteria culture:")
  print("Current number of components:")
  optim.method = match.arg(optim.method)
  result = total.function(n, j, 1, atoms, draw, Ecoff.quantile, pi_cutoff, alpha, beta, memb.exp, maxiter, eps, optim.method)
  BIC_max  = result$BIC
  #print(c(1, BIC_max))
  if(Kmax == 1)
  {
    return(result)
  }

  for(k in 2:Kmax)
  {
    print(k)
    result_new = total.function(n, j, k, atoms, draw, Ecoff.quantile, pi_cutoff, alpha, beta, memb.exp, maxiter, eps, optim.method)
    BIC_new = result_new$BIC
    #print(c(k, BIC_new))
    if(BIC_new > BIC_max)
    {
      result = result_new
      BIC_max = BIC_new
    }
  }

  return(result)
}



############################################
############################################
## Complementary functions:


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
  #param$pi = param$pi/sum(param$pi)

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
  stopifnot(sum(pi)==1)
  stopifnot(all(s>0))
  stopifnot(all(df>0))

  by_components = sapply(1:K, function(i) {
    (1/s[i]) * dt((y - mu[i])/s[i], df[i])})

  result = as.numeric(by_components %*% pi)
  return(result)
}

dt(1:9, df=2*(1:9))
sapply(1:9, function(i){dt(i, 2*i)})

dt.mixture_log <- function(y, pi, mu, s, df, log=FALSE)
{
  K = length(pi)
  total_matrix_log =  matrix(log(pi) - log(s) + dt((rep(y, each=K) - mu)/s, df, log=TRUE), ncol=K, byrow=TRUE)
  result_log = matrixStats::rowLogSumExps(total_matrix_log)

  if(log==FALSE)
  {
    return(exp(result_log))
  }
  else
  {
    return(result_log)
  }
}



#################### pt.mixture(): ############################################

#CDF of a mixture of t distributions

#Input:
#y: value(-s) at which the CDF shall be evaluated (numeric vector)
#pi: mixture proportions (numeric vector of length K > 0; sum(pi) = 1, entries non-negative)
#mu: component medians (numeric vector of lentgh K)
#s: component sigma's (numeric vector of lentgh K; entries positive)
#df: component degrees of freedom (numeric vector of lentgh K; entries positive)

#Output:
#result: F(y)

#arguments are to be checked
pt.mixture <- function(y, pi, mu, s, df)
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
  stopifnot(sum(pi)==1)
  stopifnot(all(s>0))
  stopifnot(all(df>0))

  by_components = sapply(1:K, function(i) {
    pt((y - mu[i])/s[i], df[i])})

  result = as.numeric(by_components %*% pi)
  return(result)
}


pt.mixture_log <- function(y, pi, mu, s, df, log.p=FALSE)
{
  K = length(pi)
  total_matrix_log =  matrix(log(pi) + pt((rep(y, each=K) - mu)/s, df, log.p=TRUE), byrow=TRUE, ncol=K)
  result_log = matrixStats::rowLogSumExps(total_matrix_log)

  if(log.p==FALSE)
  {
    return(exp(result_log))
  }
  else
  {
    return(result_log)
  }
}


#################### qt.mixture(): ############################################

#Quantile function of a mixture of t distributions

#Input:
#p: value at which the Quantile function shall be evaluated (numeric scalar within (0, 1))
#pi: mixture proportions (numeric vector of length K > 0; sum(pi) = 1, entries non-negative)
#mu: component medians (numeric vector of lentgh K)
#s: component sigma's (numeric vector of lentgh K; entries positive)
#df: component degrees of freedom (numeric vector of lentgh K; entries positive)
#tol: numerical tolerance of result (numeric scalar, positive, "small")

#Output:
#result: Q(y)

#Uses:
#pt.mixture()
#bisection.method()
#weighted.quantile()

#arguments are to be checked
qt.mixture <- function(p, pi, mu, s, df, tol=10*.Machine$double.eps)
{
  stopifnot(is.numeric(p))
  stopifnot(length(p) > 0)
  if(length(p) > 1)
  {
    warning("Only first value of p is used. To apply qt.mixture() to multiple values, use sapply()")
    p = p[1]
  }
  if(p==1)
    return(Inf)
  else if(p==0)
    return(-Inf)
  else if(p>1 | p<0)
    stop("p must be a real number within [0, 1]")

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
  stopifnot(sum(pi)==1)
  stopifnot(all(s>0))
  stopifnot(all(df>0))


  if(length(pi)==1)
    return(mu + s*qt(p, df))

  #Calculate starting values for bisection method
  mu_ges = weighted.quantile(x=mu, w=pi, probs=0.5)
  s_ges = as.vector(s%*%pi * length(pi))
  df_ges = min(df)

  c = abs(qt(p, df_ges))

  if(p<0.5)
  {
    a = mu_ges - c*s_ges
    b = mu_ges
  }
  else if(p>0.5)
  {
    a = mu_ges
    b = mu_ges + c*s_ges
  }
  else #p=0.5
  {
    a = mu_ges - s_ges/2
    b = mu_ges + s_ges/2
  }

  #check if interval actually contains quantile
  fa = pt.mixture(a, pi=pi, mu=mu, s=s, df=df)
  fb = pt.mixture(b, pi=pi, mu=mu, s=s, df=df)

  #widen interval (if necessary)
  while(fa > p)
  {
    a = a - s_ges/2
    fa = pt.mixture(a, pi=pi, mu=mu, s=s, df=df)
  }
  while(fb < p)
  {
    b = b + s_ges/2
    fb = pt.mixture(b, pi=pi, mu=mu, s=s, df=df)
  }

  #solve equation F(x) == p for x by bisection method
  return(bisection.method(y=p, f=pt.mixture, limits=c(a, b),
                          tol=tol, pi=pi, mu=mu, s=s, df=df))
}

qt.mixture_log <- function(p, pi, mu, s, df, tol=100*.Machine$double.eps)
{
  stopifnot(is.numeric(p))
  stopifnot(length(p) > 0)
  if(length(p) > 1)
  {
    warning("Only first value of p is used. To apply qt.mixture() to multiple values, use sapply()")
    p = p[1]
  }
  if(p==1)
    return(Inf)
  else if(p==0)
    return(-Inf)
  else if(p>1 | p<0)
    stop("p must be a real number within [0, 1]")


  #print(tol)
  if(length(pi)==1)
    return(mu + s*qt(p, df))

  p_log = log(p)

  #Calculate starting values for bisection method
  mu_ges = weighted.quantile(x=mu, w=pi, probs=0.5)
  s_ges = as.vector(s%*%pi * length(pi))
  df_ges = min(df)

  c = abs(qt(p, df_ges))

  if(p<0.5)
  {
    a = mu_ges - c*s_ges
    b = mu_ges
  }
  else if(p>0.5)
  {
    a = mu_ges
    b = mu_ges + c*s_ges
  }
  else #p=0.5
  {
    a = mu_ges - s_ges/2
    b = mu_ges + s_ges/2
  }

  #check if interval actually contains quantile & widen interval (if necessary)
  fa_log = pt.mixture_log(a, pi=pi, mu=mu, s=s, df=df, log.p=TRUE)
  fb_log = pt.mixture_log(b, pi=pi, mu=mu, s=s, df=df, log.p=TRUE)

  while(fa_log > p_log)
  {
    #print("Pepi")
    a = a - s_ges/2
    fa_log = pt.mixture_log(a, pi=pi, mu=mu, s=s, df=df, log.p=TRUE)
  }
  while(fb_log < p_log)
  {
    #print("Gustav")
    b = b + s_ges/2
    fb_log = pt.mixture_log(b, pi=pi, mu=mu, s=s, df=df, log.p=TRUE)
  }

  #print("Hansi")
  #print(c(a, b))

  #solve equation F(x) == p for x by bisection method
  return(bisection.method(y=p_log, f=pt.mixture_log, limits=c(a, b),
                          tol=tol, pi=pi, mu=mu, s=s, df=df, log.p=TRUE))
}


#####################################################
#####################################################
##Auxiliary functions (not for public use)

######################## weighted.quantile(): ##################################

#computes weighted quantile

#Input:
#x: data (numeric vector of length N)
#w: weights (numeric vector of length N, entries non-negative)
#z: which quantile(-s) shall be computed (numeric vector, entries within (0, 1))

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
  stopifnot(all(0 < z & z < 1))

  #prepare arguments
  w = w/sum(w)

  order <- order(x)
  x = x[order]
  w = w[order]

  p = length(z)
  z = sort(z)

  #prepare loop variables
  q = numeric(p)
  W_sum = 0
  i=0
  #calculate quantile
  for(j in 1:p)
  {
    while(W_sum < z[j])
    {
      i = i+1
      W_sum = W_sum + w[i]
    }

    q[j] = x[i]
  }

  return(q)
}

############################## bisection.method() ########################

#Solves equation f(x) = y numerically for x

#Input:
#y: value that f should attain (numeric scalar)
#f: function to be inverted (function; should be strictly monotonous)
#limits: search interval (vector of length 2; must contain solution)

#Output.
#x: such that f(x) = y

bisection.method <- function(y, f, limits, tol=10*.Machine$double.eps, ...)
{
  a=min(limits)
  b=max(limits)
  k=0
  while(b-a > tol)
  {
    print(k)
    print(c(a, (a+b)/2, b))

    k = k+1
    fa = f(a, ...)
    fb = f(b, ...)
    fm = f((a+b)/2, ...)

    if(abs(fa - y) <= tol)
      return(a)
    else if(abs(fb - y) <= tol)
      return(b)
    else if(abs(fm - y) <= tol)
      return((a+b)/2)
    else if(fa <= y & y <= fm)
      b = (a+b)/2
    else if(fa >= y & y >= fm)
      b = (a+b)/2
    else if(fm <= y & y <= fb)
      a = (a+b)/2
    else if(fm >= y & y >= fb)
      a = (a+b)/2
    else
      stop("Function 'f' might not be monotonous within search space")

    #print(b-a <= tol)
  }

  return((a+b)/2)
}


################ ungroup.data(): ############################
#function that ungroups data (if N>700, the resulting vector is not the actual sample, but a close representation of it with length 500)


#Input:
#j: possible values (numeric vector of length J > 0)
#n: frequencies (integer vector of length J; non-negative entries)

#Output:
#y: ungrouped data vector (length: min(500, N))

#Uses:
#weighted.quantile()

ungroup.data <- function(n, j)
{
  N = sum(n)
  J=length(j)
  if(N > 500) #take representative observations
  {
    y = weighted.quantile(j, n, probs=(1:499)/500)
    n_sample = table(y)
    j_sample = as.numeric(names(n_sample))
  }
  else #ungroup data
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
  return(y)
}


logLik_t_mixture<-function(n, j, pi, mu, s, df)
{
  K=length(pi)
  uppers_log = pt((rep(j,each=K)+1/2-mu)/s, df, log.p=TRUE)
  lowers_log = pt((rep(j,each=K)-1/2-mu)/s, df, log.p=TRUE)

  numerators_log <- log(pi) + matrix(uppers_log + VGAM::log1mexp(uppers_log-lowers_log), nrow=K)
  denominators_log = matrixStats::colLogSumExps(numerators_log)

  return(sum(denominators_log*n))
}


BIC_t_mixture<-function(n, j, pi, mu, s, df)
{
  N = sum(n)
  K = length(pi)
  #print(c(2*logLik.t.mixture(n, j, pi, mu, s, df), logLik.t.mixture(n, j, pi, mu, s, df)))
  #print(c(log(N), (4*K-1), (4*K-1)*log(N)))

  return(2*logLik_t_mixture(n, j, pi, mu, s, df) - (4*K-1) * log(N))
}


AIC_t_mixture<-function(n, j, pi, mu, s, df)
{
  N = sum(n)
  K = length(pi)
  return(2*logLik_t_mixture(n, j, pi, mu, s, df) - (4*K-1) * 2)
}


#################### Q1k etc. #####################################
Q2 <- function(par, tau, j, n)
{
  upper_log = pt((j+1/2-par[1])/exp(par[2]), exp(par[3]), log.p=TRUE)
  lower_log = pt((j-1/2-par[1])/exp(par[2]), exp(par[3]), log.p=TRUE)
  return(sum((upper_log + VGAM::log1mexp(upper_log - lower_log)) * tau * n))
}


Q3 <- function(s, a, b)
{
  return(dinvgamma(exp(s)^2, shape= a, scale= b, log=TRUE))
}

Q4 <- function(df)
{
  return(ifelse(exp(df) <= 10^12, 0, -Inf))
}


Qmax <- function(par, tau, j, n, a, b)
{
  return(-(Q2(par, tau, j, n) + Q3(par[2], a, b) + Q4(par[3])))
}

Q <- function(par, j, n)
{
  K=dim(par)[1]
  uppers_log = pt((rep(j,each=K)+1/2-par[,2])/exp(par[,3]), exp(par[,4]), log.p=TRUE)
  lowers_log = pt((rep(j,each=K)-1/2-par[,2])/exp(par[,3]), exp(par[,4]), log.p=TRUE)
  numerators_log <- log(par[,1]) + matrix(uppers_log + VGAM::log1mexp(uppers_log-lowers_log), nrow=K)
  denominators_log = matrixStats::colLogSumExps(numerators_log)

  tau <- exp(t(t(numerators_log)-denominators_log))

  return(sum(n * tau * log(par[,1])) + Q2(par[,2:4], tau, j, n))
}


