logLik.t.mixture.old<-function(n, j, pi, mu, s, df)
{
  uppers_log = pt((rep(j,each=K)+1/2-mu)/s, df, log.p=TRUE)
  lowers_log = pt((rep(j,each=K)-1/2-mu)/s, df, log.p=TRUE)

  numerators_log <- log(pi) + matrix(uppers_log + VGAM::log1mexp(uppers_log-lowers_log), nrow=K)
  denominators_log = matrixStats::colLogSumExps(numerators_log)
  return(sum(n*denominators_log))
}

