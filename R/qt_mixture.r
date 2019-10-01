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
  stopifnot(1 - 10*.Machine$double.eps <= sum(pi) & sum(pi) <= 1 + 10*.Machine$double.eps)
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
