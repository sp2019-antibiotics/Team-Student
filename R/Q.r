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


