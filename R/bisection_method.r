############################## bisection.method() ########################

#Solves equation f(x) = y numerically for x

#Input:
#y: value that f should attain (numeric scalar)
#f: function to be inverted (function; should be strictly monotonous)
#limits: search interval (vector of length 2; must contain solution)

#Output.
#x: such that f(x) = y

bisection.method <- function(y, f, limits, tol=10*.Machine$double.eps, ...)
{
  a=min(limits)
  b=max(limits)
  k=0
  while(b-a > tol)
  {
    k = k+1
    fa = f(a, ...)
    fb = f(b, ...)
    fm = f((a+b)/2, ...)

    if(abs(fa - y) <= tol)
      return(a)
    else if(abs(fb - y) <= tol)
      return(b)
    else if(abs(fm - y) <= tol)
      return((a+b)/2)
    else if(fa <= y & y <= fm)
      b = (a+b)/2
    else if(fa >= y & y >= fm)
      b = (a+b)/2
    else if(fm <= y & y <= fb)
      a = (a+b)/2
    else if(fm >= y & y >= fb)
      a = (a+b)/2
    else
      stop("Function 'f' might not be monotonous within search space")
  }

  return((a+b)/2)
}
