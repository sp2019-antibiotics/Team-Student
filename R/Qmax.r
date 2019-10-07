Qmax <- function(par, tau, j, n, a, b)
{
  return(-(Q2(par, tau, j, n) + Q3(par[2], a, b) + Q4(par[3])))
}


