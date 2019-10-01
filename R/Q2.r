Q2 <- function(par, tau, j, n)
{
  upper_log = pt((j+1/2-par[1])/exp(par[2]), exp(par[3]), log.p=TRUE)
  lower_log = pt((j-1/2-par[1])/exp(par[2]), exp(par[3]), log.p=TRUE)
  return(sum((upper_log + VGAM::log1mexp(upper_log - lower_log)) * tau * n))
}
