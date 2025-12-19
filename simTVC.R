# df: ID, time, des

sim_TVC <- function(seed, df, chfun, hstr, theta, beta, npar) {

  # Calculating the cumulative hazard function at all time points
  if (npar == 2) {
    CH <- CH_TVC(
      df    = df,
      beta  = beta,
      npar  = npar,
      ae0   = theta[1],
      be0   = theta[2],
      chfun = chfun,
      hstr = hstr
    )
  }

  if (npar == 3) {
    CH <- CH_TVC(
      df    = df,
      beta  = beta,
      npar  = npar,
      ae0   = theta[1],
      be0   = theta[2],
      ce0   = theta[3],
      chfun = chfun,
      hstr = hstr
    )
  }

  sim <- rep(0,n)
  status <- rep(0,n)

  n <- length(unique(df$ID))

  set.seed(seed)
  u = runif(n)

  if (npar == 2) {
 survs <- SPred_TVC(df = df,
            beta = beta,
            npar = npar,
            ae0   = theta[1],
            be0   = theta[2],
            chfun = chfun,
            hstr = hstr)
  }

  if (npar == 3) {
   survs <- SPred_TVC(df = df,
              beta = beta,
              npar = npar,
              ae0   = theta[1],
              be0   = theta[2],
              ce0   = theta[3],
              chfun = chfun,
              hstr = hstr)
  }


  times <- as.vector(with(df, tapply(time, ID, max)))


ind_sol <- (survs < u)

sim[!ind_sol] <- times[!ind_sol]

  for (i in c(1:n)[ind_sol]) {

    print(i)

  }

  out = list(sim = sim, status = status)

  return(out)
}
