BIC.t.mixture<-function(n, j, pi, mu, s, df)
{
  N = sum(n)
  K = length(pi)
  return(2*logLik.t.mixture(n, j, pi, mu, s, df) - (4*K-1) * log(N))
}
