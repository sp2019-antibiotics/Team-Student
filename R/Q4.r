Q4 <- function(df)
{
  return(ifelse(exp(df) <= 10^12, 0, -Inf))
}

