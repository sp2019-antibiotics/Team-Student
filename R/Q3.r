Q3 <- function(s, a, b)
{
  return(invgamma::dinvgamma(exp(s)^2, shape= a, scale= b, log=TRUE))
}
