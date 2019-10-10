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
  }
  else #ungroup data
  {
    y = rep(j, n)
  }
  return(y)
}
