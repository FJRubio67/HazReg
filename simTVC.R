sim_TVC <- function(seed, df, chfun, hstr, theta, beta, npar) {

  set.seed(seed)

  n <- length(unique(df$ID))
  sim    <- numeric(n)
  status <- integer(n)

  ## Survival at last observed time
  if (npar == 2) {
    survs <- SPred_TVC(
      df    = df,
      beta  = beta,
      npar  = npar,
      ae0   = theta[1],
      be0   = theta[2],
      chfun = chfun,
      hstr  = hstr
    )
  } else {
    survs <- SPred_TVC(
      df    = df,
      beta  = beta,
      npar  = npar,
      ae0   = theta[1],
      be0   = theta[2],
      ce0   = theta[3],
      chfun = chfun,
      hstr  = hstr
    )
  }

  ## Maximum follow-up times
  times <- as.vector(with(df, tapply(time, ID, max)))

  ## Uniform draws
  u <- runif(n)

  ## Censoring indicator
  censored <- (survs > u)

  sim[censored]    <- times[censored]
  status[censored] <- 0

  ## Event times
  for (i in which(!censored)) {

    rootfun <- function(t) {
      SPred_TVC_i(
        df    = df,
        i     = i,
        t     = t,
        beta  = beta,
        npar  = npar,
        ae0   = theta[1],
        be0   = theta[2],
        ce0   = if (npar == 3) theta[3] else NULL,
        chfun = chfun,
        hstr  = hstr
      ) - u[i]
    }

    ## Safety check for bracketing
    if (rootfun(0) < 0 || rootfun(times[i]) > 0) {
      stop("Root not bracketed for individual ", i)
    }

    sim[i] <- uniroot(rootfun, interval = c(0, times[i]))$root
    status[i] <- 1
  }

  list(time = sim, status = status)
}
