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
  stopifnot(1 - 10*.Machine$double.eps <= sum(par[,1]) & sum(par[,1]) <= 1 + 10*.Machine$double.eps)
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
