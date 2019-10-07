visualization.t.mixture <- function(n, j, pi, mu, s, df, ecoff=NULL, kernel.density=TRUE, hist=TRUE, new.plot=hist, pi_resistant=0)
{
  #check arguments
  stopifnot(is.numeric(n))
  stopifnot(is.vector(n))

  n = as.integer(n)
  stopifnot(all(n>=0))
  N = sum(n)
  stopifnot(N>0)

  stopifnot(is.numeric(j))
  stopifnot(is.vector(j))
  stopifnot(length(j) == length(n))

  j = j[n>0]
  n = n[n>0]

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

  stopifnot(is.logical(kernel.density))
  stopifnot(is.logical(hist))
  stopifnot(is.logical(new.plot))

  stopifnot(is.numeric(pi_resistant))
  stopifnot(length(pi_resistant)==1)
  stopifnot(0 <= pi_resistant & pi_resistant <= 1)

  #ungroup and plot data
  y = ungroup.data(n, j)

  max.data <- max(y)
  min.data <- min(y)
  hist.step <- (max.data - min.data)

  t = seq(from = min.data, to = max.data, by = 0.1)
  ft = dt.mixture(t, pi, mu, s, df)
  ft = ft/(1+pi_resistant)

  dens = density(y)

  ylim2 = max(ft, dens$y)

  if (hist.step > 50)
    hist.step <- hist.step
  h=hist(
    y,
    breaks = seq(min.data, max.data, l = hist.step),
    main="Density Estimates",
    probability = TRUE,
    xlim = c(min.data, max.data),
    ylim = c(0, ylim2),
    col="grey"
    ,xaxt = 'n'
  )

  if(max(h$density) > ylim2)
  {
    ylim2 = max(h$density)
    h=hist(
      y,
      breaks = seq(min.data, max.data, l = hist.step),
      main="Density Estimates",
      probability = TRUE,
      xlim = c(min.data, max.data),
      ylim = c(0, max(h$density)),
      col="grey"
      ,xaxt = 'n'
    )

  }
  axis(side = 1,
       at = seq(min.data, max.data, (max.data - min.data) / hist.step))
  lines(dens, col = "red", lwd = 2)
  lines(t, ft, col = "blue", lwd = 2)

  if(!is.null(ecoff))
  {
    stopifnot(is.numeric(ecoff))
    stopifnot(length(ecoff)==1)

    if(ecoff < min.data | ecoff > max.data)
      warning("Forwarded ECOFF is out of range of j and will be 'left' or 'right' to the plot.")

    abline(v=ecoff, lty="dashed", lwd=2)
    legend("topleft", col = c("red", "blue", "black"),
           legend = c("Kernel Density", "t-Density", "ECOFF"), lty = c(1, 1, 2), cex = 0.75, lwd = 2)
  }
  else
  {
    legend("topleft", col = c("red", "blue"),
           legend = c("Kernel Density", "t-Density"), lty = 1, cex = 0.75, lwd = 2)
  }
}
