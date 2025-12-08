# Compute the cumulative hazard function for at several times points for each individual
# Proportional hazards models
# 2-parameter models

compute_CHPH2 <- function(df, beta, ae0, be0, chfun) {

  # Extract design matrix
  Xmat <- as.matrix(df[, grep("^des", names(df))])

  # exp(x beta)
  exp_xb <- exp(Xmat %*% beta)

  # baseline cumulative hazard at each time
  ch0_t <- chfun(df$time, ae0, be0)

  # final cumulative hazard contribution
  CH_t <- as.vector(ch0_t * exp_xb)

  return(CH_t)
}

# Compute the cumulative hazard function for at several times points for each individual
# Proportional hazards models
# 3-parameter models

compute_CHPH3 <- function(df, beta, ae0, be0, ce0, chfun) {

  # Extract design matrix
  Xmat <- as.matrix(df[, grep("^des", names(df))])

  # exp(x beta)
  exp_xb <- exp(Xmat %*% beta)

  # baseline cumulative hazard at each time
  ch0_t <- chfun(df$time, ae0, be0, ce0)

  # final cumulative hazard contribution
  CH_t <- as.vector(ch0_t * exp_xb)

  return(CH_t)
}


# Compute the cumulative hazard function for at several times points for each individual
# AFT models
# 2-parameter models

compute_CHAFT2 <- function(df, beta, ae0, be0, chfun) {

  # Extract design matrix
  Xmat <- as.matrix(df[, grep("^des", names(df))])

  # exp(x beta)
  exp_xb <- exp(Xmat %*% beta)

  # baseline cumulative hazard at each time * exp(x^T beta)
  ch0_t <- chfun(df$time*exp_xb, ae0, be0)

  # final cumulative hazard contribution
  CH_t <- as.vector(ch0_t)

  return(CH_t)
}

# Compute the cumulative hazard function for at several times points for each individual
# AFT models
# 3-parameter models

compute_CHAFT3 <- function(df, beta, ae0, be0, ce0, chfun) {

  # Extract design matrix
  Xmat <- as.matrix(df[, grep("^des", names(df))])

  # exp(x beta)
  exp_xb <- exp(Xmat %*% beta)

  # baseline cumulative hazard at each time * exp(x^T beta)
  ch0_t <- chfun(df$time*exp_xb, ae0, be0, ce0)

  # final cumulative hazard contribution
  CH_t <- as.vector(ch0_t)

  return(CH_t)
}

# df: matrix or data frame containing 'ID' (individual ID), 'time' (contact time points), 'des' (design matrix).
# 'time' must be in increasing order.
# The last element in 'time' must be the time to event (either observed or right-censoring time)
# The design matrix 'des' must contain the covariate values at each of the contact time points in 'time'

HMLE_TVC <- function (init, df, status, hstr = NULL, dist = NULL, des = NULL,
          method = "Nelder-Mead", maxit = 100)
{

  df <- df[order(df$ID, df$time), ]  # ensure sorted

  times <- as.vector(with(df, tapply(time, ID, max)))

  last_rows <- df[ave(df$time, df$ID, FUN = max) == df$time, ]
  des <- as.matrix(last_rows[, grep("^des", names(df))])


  status <- as.vector(as.logical(status))
  times.obs <- times[status]
  des_obs <- des[status, ]


#-------------------------------------------------------------------------------
# PH
#-------------------------------------------------------------------------------
  if (hstr == "PH") {
    if (dist == "PGW") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        beta <- par[4:(3 + p)]
        # Hazard calculations
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hpgw(times.obs, ae0, be0, ce0, log = TRUE) +
          x.beta.obs

      # Cumulative hazard calculations
       df$CH_t <- compute_CHPH3(df, beta, ae0, be0, ce0, chpgw)

       # split by ID
       lst <- split(df, df$ID)

       # compute lagged differences within each ID
       res_list <- lapply(lst, function(d) {
         data.frame(
           ID = d$ID[-1],                # drop the first (no diff)
           time = d$time[-1],            # time of the difference
           diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
         )
       })


       # sum diff_chaz for each ID
       chaz0 <- as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

       # Negative log-likelihood
        val <- -sum(lhaz0) + sum(chaz0)

        return(val)
      }
    }
    if (dist == "EW") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        beta <- par[4:(3 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        lhaz0 <- hew(times.obs, ae0, be0, ce0, log = TRUE) +
          x.beta.obs

        # Cumulative hazard calculations
        df$CH_t <- compute_CHPH3(df, beta, ae0, be0, ce0, chew)

        # split by ID
        lst <- split(df, df$ID)

        # compute lagged differences within each ID
        res_list <- lapply(lst, function(d) {
          data.frame(
            ID = d$ID[-1],                # drop the first (no diff)
            time = d$time[-1],            # time of the difference
            diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
          )
        })


        # sum diff_chaz for each ID
        chaz0 <- as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

        # Negative log-likelihood
        val <- -sum(lhaz0) + sum(chaz0)

        return(val)
      }
    }
    if (dist == "GenGamma") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        beta <- par[4:(3 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        lhaz0 <- hggamma(times.obs, ae0, be0, ce0, log = TRUE) +
          x.beta.obs

        # Cumulative hazard calculations
        df$CH_t <- compute_CHPH3(df, beta, ae0, be0, ce0, chggama)

        # split by ID
        lst <- split(df, df$ID)

        # compute lagged differences within each ID
        res_list <- lapply(lst, function(d) {
          data.frame(
            ID = d$ID[-1],                # drop the first (no diff)
            time = d$time[-1],            # time of the difference
            diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
          )
        })


        # sum diff_chaz for each ID
        chaz0 <- as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

        # Negative log-likelihood
        val <- -sum(lhaz0) + sum(chaz0)

        return(val)
      }
    }
    if (dist == "LogNormal") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        lhaz0 <- hlnorm(times.obs, ae0, be0, log = TRUE) +
          x.beta.obs

        # Cumulative hazard calculations
        df$CH_t <- compute_CHPH2(df, beta, ae0, be0, chlnorm)

        # split by ID
        lst <- split(df, df$ID)

        # compute lagged differences within each ID
        res_list <- lapply(lst, function(d) {
          data.frame(
            ID = d$ID[-1],                # drop the first (no diff)
            time = d$time[-1],            # time of the difference
            diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
          )
        })


        # sum diff_chaz for each ID
        chaz0 <- as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

        # Negative log-likelihood
        val <- -sum(lhaz0) + sum(chaz0)

        return(val)
      }
    }
    if (dist == "LogLogistic") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        lhaz0 <- hllogis(times.obs, ae0, be0, log = TRUE) +
          x.beta.obs

        # Cumulative hazard calculations
        df$CH_t <- compute_CHPH2(df, beta, ae0, be0, chllogis)

        # split by ID
        lst <- split(df, df$ID)

        # compute lagged differences within each ID
        res_list <- lapply(lst, function(d) {
          data.frame(
            ID = d$ID[-1],                # drop the first (no diff)
            time = d$time[-1],            # time of the difference
            diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
          )
        })


        # sum diff_chaz for each ID
        chaz0 <- as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

        # Negative log-likelihood
        val <- -sum(lhaz0) + sum(chaz0)


      }
    }
    if (dist == "Gamma") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        lhaz0 <- hgamma(times.obs, ae0, be0, log = TRUE) +
          x.beta.obs

        # Cumulative hazard calculations
        df$CH_t <- compute_CHPH2(df, beta, ae0, be0, chgamma)

        # split by ID
        lst <- split(df, df$ID)

        # compute lagged differences within each ID
        res_list <- lapply(lst, function(d) {
          data.frame(
            ID = d$ID[-1],                # drop the first (no diff)
            time = d$time[-1],            # time of the difference
            diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
          )
        })


        # sum diff_chaz for each ID
        chaz0 <- as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

        # Negative log-likelihood
        val <- -sum(lhaz0) + sum(chaz0)

        return(val)
      }
    }
    if (dist == "Weibull") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        lhaz0 <- hweibull(times.obs, ae0, be0, log = TRUE) +
          x.beta.obs

        # Cumulative hazard calculations
        df$CH_t <- compute_CHPH2(df, beta, ae0, be0, chweibull)

        # split by ID
        lst <- split(df, df$ID)

        # compute lagged differences within each ID
        res_list <- lapply(lst, function(d) {
          data.frame(
            ID = d$ID[-1],                # drop the first (no diff)
            time = d$time[-1],            # time of the difference
            diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
          )
        })


        # sum diff_chaz for each ID
        chaz0 <- as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

        # Negative log-likelihood
        val <- -sum(lhaz0) + sum(chaz0)

        return(val)
      }
    }
  }
#-------------------------------------------------------------------------------
# AFT
#-------------------------------------------------------------------------------
  if (hstr == "AFT") {
    if (dist == "PGW") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        beta <- par[4:(3 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hpgw(times.obs * exp.x.beta.obs, ae0,
                      be0, ce0, log = TRUE) + x.beta.obs

        # Cumulative hazard calculations
        df$CH_t <- compute_CHAFT3(df, beta, ae0, be0, ce0, chpgw)

        # split by ID
        lst <- split(df, df$ID)

        # compute lagged differences within each ID
        res_list <- lapply(lst, function(d) {
          data.frame(
            ID = d$ID[-1],                # drop the first (no diff)
            time = d$time[-1],            # time of the difference
            diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
          )
        })


        # sum diff_chaz for each ID
        chaz0 <- as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

        # Negative log-likelihood
        val <- -sum(lhaz0) + sum(chaz0)

        return(val)
      }
    }
    if (dist == "EW") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        beta <- par[4:(3 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hew(times.obs * exp.x.beta.obs, ae0,
                     be0, ce0, log = TRUE) + x.beta.obs

        # Cumulative hazard calculations
        df$CH_t <- compute_CHAFT3(df, beta, ae0, be0, ce0, chew)

        # split by ID
        lst <- split(df, df$ID)

        # compute lagged differences within each ID
        res_list <- lapply(lst, function(d) {
          data.frame(
            ID = d$ID[-1],                # drop the first (no diff)
            time = d$time[-1],            # time of the difference
            diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
          )
        })


        # sum diff_chaz for each ID
        chaz0 <- as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

        # Negative log-likelihood
        val <- -sum(lhaz0) + sum(chaz0)

        return(val)
      }
    }
    if (dist == "GenGamma") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        beta <- par[4:(3 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hggamma(times.obs * exp.x.beta.obs,
                         ae0, be0, ce0, log = TRUE) + x.beta.obs

        # Cumulative hazard calculations
        df$CH_t <- compute_CHAFT3(df, beta, ae0, be0, ce0, chggamma)

        # split by ID
        lst <- split(df, df$ID)

        # compute lagged differences within each ID
        res_list <- lapply(lst, function(d) {
          data.frame(
            ID = d$ID[-1],                # drop the first (no diff)
            time = d$time[-1],            # time of the difference
            diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
          )
        })


        # sum diff_chaz for each ID
        chaz0 <- as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

        # Negative log-likelihood
        val <- -sum(lhaz0) + sum(chaz0)

        return(val)
      }
    }
    if (dist == "LogNormal") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hlnorm(times.obs * exp.x.beta.obs, ae0,
                        be0, log = TRUE) + x.beta.obs

        # Cumulative hazard calculations
        df$CH_t <- compute_CHAFT2(df, beta, ae0, be0, chlnorm)

        # split by ID
        lst <- split(df, df$ID)

        # compute lagged differences within each ID
        res_list <- lapply(lst, function(d) {
          data.frame(
            ID = d$ID[-1],                # drop the first (no diff)
            time = d$time[-1],            # time of the difference
            diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
          )
        })


        # sum diff_chaz for each ID
        chaz0 <- as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

        # Negative log-likelihood
        val <- -sum(lhaz0) + sum(chaz0)

        return(val)
      }
    }
    if (dist == "LogLogistic") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hllogis(times.obs * exp.x.beta.obs,
                         ae0, be0, log = TRUE) + x.beta.obs

        # Cumulative hazard calculations
        df$CH_t <- compute_CHAFT2(df, beta, ae0, be0, chllogis)

        # split by ID
        lst <- split(df, df$ID)

        # compute lagged differences within each ID
        res_list <- lapply(lst, function(d) {
          data.frame(
            ID = d$ID[-1],                # drop the first (no diff)
            time = d$time[-1],            # time of the difference
            diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
          )
        })


        # sum diff_chaz for each ID
        chaz0 <- as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

        # Negative log-likelihood
        val <- -sum(lhaz0) + sum(chaz0)

        return(val)
      }
    }
    if (dist == "Gamma") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hgamma(times.obs * exp.x.beta.obs, ae0,
                        be0, log = TRUE) + x.beta.obs

        # Cumulative hazard calculations
        df$CH_t <- compute_CHAFT2(df, beta, ae0, be0, chgamma)

        # split by ID
        lst <- split(df, df$ID)

        # compute lagged differences within each ID
        res_list <- lapply(lst, function(d) {
          data.frame(
            ID = d$ID[-1],                # drop the first (no diff)
            time = d$time[-1],            # time of the difference
            diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
          )
        })


        # sum diff_chaz for each ID
        chaz0 <- as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

        # Negative log-likelihood
        val <- -sum(lhaz0) + sum(chaz0)

        return(val)
      }
    }
    if (dist == "Weibull") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hweibull(times.obs * exp.x.beta.obs,
                          ae0, be0, log = TRUE) + x.beta.obs

        # Cumulative hazard calculations
        df$CH_t <- compute_CHAFT2(df, beta, ae0, be0, chweibull)

        # split by ID
        lst <- split(df, df$ID)

        # compute lagged differences within each ID
        res_list <- lapply(lst, function(d) {
          data.frame(
            ID = d$ID[-1],                # drop the first (no diff)
            time = d$time[-1],            # time of the difference
            diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
          )
        })


        # sum diff_chaz for each ID
        chaz0 <- as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

        # Negative log-likelihood
        val <- -sum(lhaz0) + sum(chaz0)

        return(val)
      }
    }
  }
# Optimisation step
  if (method != "nlminb") {
    OPT <- optim(init, log.lik, control = list(maxit = maxit),
                 method = method)
  }
  if (method == "nlminb") {
    OPT <- nlminb(init, log.lik, control = list(iter.max = maxit))
  }
  OUT <- list(log_lik = log.lik, OPT = OPT)
  return(OUT)
}
