

################################################################################
################################################################################
################################################################################
# Time-varying covariates
################################################################################
################################################################################
################################################################################


#----------------------------------------------------------------------------------------
#' Compute the Cumulative Hazard for a Proportional Hazards or Accelerated Failure Model
#' (2- and 3-parameter baseline)
#'
#' Computes the cumulative hazard \eqn{H(t \mid x(t))} at multiple time points
#' for each individual under a proportional hazards (PH) model or
#' Accelerated Failure Time (AFT) model with a
#' two-parameter or three-parameter parametric baseline hazard.
#'
#' The PH model assumes
#' \deqn{H(t \mid x(t)) = H_0(t; a_0, b_0, c_0)\exp(x(t)^\top \beta).}
#'
#' In the AFT model, event time is rescaled as
#' \deqn{H(t \mid x(t)) = H_0(t \exp(x(t)^\top\beta); a_0, b_0, c_0).}
#'
#' @param df A data frame containing:
#'   \itemize{
#'     \item `time`: numeric vector of time points.
#'     \item Covariate columns named with prefix `"des"` (e.g., `des1`, `des2`, ...),
#'           representing \eqn{x(t)}.
#'   }
#' @param beta Numeric vector of regression coefficients.
#' @param theta Numeric baseline parameters of the cumulative hazard.
#' @param chfun A function computing the baseline cumulative hazard:
#'   `chfun(time, theta[1], theta[2])` or `chfun(time, theta[1], theta[2], theta[3])`.
#' @param hstr Hazard structure ("PH" or "AFT")
#'
#' @return A numeric vector with the cumulative hazard evaluated at each time
#'   point in `df`.
#'
#' @export
CH_TVC <- function(df, beta,
                   theta,
                   chfun, hstr) {

  npar = length(theta)

  # Original order
  df$original_order_idx <- seq_len(nrow(df))

  ## Ensure sorted data
  df <- df[order(df$ID, df$time), ]

  ## Extract design matrix
  Xmat <- as.matrix(df[, grep("^des", names(df))])
  exp_xb <- as.vector(exp(Xmat %*% beta))

  ## Baseline cumulative hazard
  H0 <- function(tt) {
    if (npar == 2){
      out <-    chfun(tt, theta[1], theta[2])
    }
    if (npar == 3){
      out <-    chfun(tt,  theta[1], theta[2], theta[3])
    }
    return(out)
  }

  ## Split by individual
  split_df <- split(seq_len(nrow(df)), df$ID)

  ## Storage
  H_out <- numeric(nrow(df))

  for (idx in split_df) {

    t <- df$time[idx]
    xb <- exp_xb[idx]

    if (hstr == "PH") {

      H0_t <- H0(t)
      dH0  <- diff(c(0, H0_t))
      H_i  <- cumsum(dH0 * xb)

    }

    if (hstr == "AFT") {

      H0_t <- H0(t * xb)
      dH0  <- diff(c(0, H0_t))
      H_i  <- cumsum(dH0)

    }

    H_out[idx] <- H_i
  }

  df$cum_hazard <- H_out

  ## 4. Restore original order and remove the temporary index
  df <- df[order(df$original_order_idx), ]
  df$original_order_idx <- NULL

  return(df)
}


#----------------------------------------------------------------------------------------
#' Compute the Survival Function for a Proportional Hazards or Accelerated Failure Model
#' (2- and 3-parameter baseline)
#'
#' Computes the survival function \eqn{exp(-H(t \mid x(t)))} at the last time point
#' for each individual under a proportional hazards (PH) model or
#' Accelerated Failure Time (AFT) model with a
#' two-parameter or three-parameter parametric baseline hazard.
#'
#' The PH model assumes
#' \deqn{H(t \mid x(t)) = H_0(t; a_0, b_0, c_0)\exp(x(t)^\top \beta).}
#'
#' In the AFT model, event time is rescaled as
#' \deqn{H(t \mid x(t)) = H_0(t \exp(x(t)^\top\beta); a_0, b_0, c_0).}
#'
#' @param df A data frame containing:
#'   \itemize{
#'     \item `time`: numeric vector of time points.
#'     \item Covariate columns named with prefix `"des"` (e.g., `des1`, `des2`, ...),
#'           representing \eqn{x(t)}.
#'   }
#' @param beta Numeric vector of regression coefficients.
#' @param theta Numeric baseline parameters of the cumulative hazard.
#' @param chfun A function computing the baseline cumulative hazard:
#'   `chfun(time, theta[1], theta[2])` or `chfun(time, theta[1], theta[2], theta[3])`.
#' @param hstr Hazard structure ("PH" or "AFT")
#'
#' @return A numeric vector with the survival function evaluated at the last time
#'   point in `df`.
#'
#' @export
SPred_TVC <-
  function(df, beta, theta, chfun, hstr) {
    # Sample size
    n <- max(df$ID)


    # Calculating the cumulative hazard function at all time points

    CH <- CH_TVC(
      df    = df,
      beta  = beta,
      theta = theta,
      chfun = chfun,
      hstr = hstr
    )$cum_hazard

    ## Extract last cumulative hazard per individual
    #    last_idx <- unique(ave(seq_along(CH), df$ID, FUN = max))
    #    H_last   <- CH[last_idx]
    H_last <- tapply(CH, df$ID, function(x) x[length(x)])

    ## Survival at last time
    S_last <- exp(-H_last)

    return(S_last)
  }


#----------------------------------------------------------------------------------------
#' Compute the Individual Survival Function at an Arbitrary Time Point
#' for a Proportional Hazards or Accelerated Failure Model
#' (2- and 3-parameter baseline)
#'
#' Computes the individual-specific survival function
#' \eqn{S_i(t) = \exp\{-H_i(t \mid x_i(t))\}}
#' at a user-specified time point \eqn{t} for a given individual \eqn{i},
#' under a proportional hazards (PH) model or an accelerated failure time (AFT) model
#' with time-varying covariates.
#'
#' The cumulative hazard for individual \eqn{i} is defined as:
#'
#' \deqn{
#' H_i(t \mid x_i(t)) = \int_0^t h_0(u; a_0, b_0, c_0)
#' \exp\{x_i(u)^\top \beta\} \, du
#' }
#'
#' for the PH model, and
#'
#' \deqn{
#' H_i(t \mid x_i(t)) = H_0\!\left(t \exp\{x_i(t)^\top \beta\};
#' a_0, b_0, c_0 \right)
#' }
#'
#' for the AFT model.
#'
#' Time-varying covariates are assumed to be piecewise constant between
#' the observation times provided in `df`. If the evaluation time `t`
#' does not coincide with an observed time point for individual `i`,
#' the covariate values are taken to be those at the most recent time
#' strictly less than `t`.
#'
#' @param df A data frame containing:
#'   \itemize{
#'     \item `ID`: individual identifier.
#'     \item `time`: numeric vector of observation times.
#'     \item Covariate columns named with prefix `"des"` (e.g., `des1`, `des2`, ...),
#'           representing the time-varying covariate process \eqn{x_i(t)}.
#'   }
#' @param i Integer specifying the individual for whom the survival
#'   function is to be evaluated.
#' @param t Numeric value giving the time point at which the survival
#'   function is evaluated. Must lie within the observation window of
#'   individual `i`.
#' @param beta Numeric vector of regression coefficients.
#' @param theta Numeric baseline parameters of the cumulative hazard.
#' @param chfun A function computing the baseline cumulative hazard:
#'   `chfun(time, theta[1], theta[2])` or `chfun(time, theta[1], theta[2], theta[3])`.
#' @param hstr Character string specifying the hazard structure:
#'   `"PH"` for proportional hazards or `"AFT"` for accelerated failure time.
#'
#' @return A numeric scalar giving the survival probability
#'   \eqn{S_i(t)} for individual `i` at time `t`.
#'
#' @details
#' This function relies on `CH_TVC()` to compute the cumulative hazard
#' over the observed time grid augmented with the evaluation time `t`.
#' The cumulative hazard is then extracted at time `t`, and the survival
#' function is obtained as \eqn{\exp\{-H_i(t)\}}.
#'
#' The resulting survival function is continuous in time but may exhibit
#' changes in slope at covariate change points, reflecting the
#' piecewise-constant nature of the time-varying covariates.
#'
#' @seealso \code{\link{CH_TVC}}, \code{\link{SPred_TVC}}
#'
#' @export
SPred_TVC_i <- function(
    df, i, t,
    beta,
    theta,
    chfun,
    hstr
) {

  ## 1. Subset individual data
  dfi <- df[df$ID == i, ]
  dfi <- dfi[order(dfi$time), ]

  if (t < min(dfi$time) || t > max(dfi$time))
    stop("t must lie within the individual's observation window")

  ## 2. Add time t if needed (piecewise-constant covariates)
  if (!any(abs(dfi$time - t) < .Machine$double.eps)) {
    idx <- max(which(dfi$time < t))
    newrow <- dfi[idx, ]
    newrow$time <- t
    dfi <- rbind(dfi, newrow)
    dfi <- dfi[order(dfi$time), ]
  }

  ## 3. Compute cumulative hazard via existing engine
  CH <- CH_TVC(
    df    = dfi,
    beta  = beta,
    theta = theta,
    chfun = chfun,
    hstr  = hstr
  )$cum_hazard

  ## 4.  Extract cumulative hazard at time t
  H_i_t <- CH[which.min(abs(dfi$time - t))]

  ## 5. Survival
  S_i_t <- exp(-H_i_t)

  return(S_i_t)
}


#----------------------------------------------------------------------------------------
#' Simulate Event Times from a Time-Varying Covariates Survival Model
#'
#' Simulates event times from a survival model with time-varying covariates
#' using the probability integral transform. The function supports
#' proportional hazards (PH) and accelerated failure time (AFT) models
#' with parametric baseline hazards.
#'
#' For each individual \eqn{i}, a uniform random variable \eqn{U_i \sim \mathrm{Unif}(0,1)}
#' is generated. If the survival probability at the last observed time
#' \eqn{t_{\max,i}} satisfies
#' \deqn{S_i(t_{\max,i}) > U_i,}
#' the individual is right-censored at \eqn{t_{\max,i}}. Otherwise, the event
#' time \eqn{T_i} is obtained by solving
#' \deqn{S_i(t) = U_i}
#' for \eqn{t \in (0, t_{\max,i}]}, where \eqn{S_i(t)} is the individual-specific
#' survival function accounting for time-varying covariates.
#'
#' The survival function is evaluated via \code{\link{SPred_TVC}} and
#' \code{\link{SPred_TVC_i}}, ensuring consistency with the underlying
#' cumulative hazard specification.
#'
#' @param seed Integer random seed for reproducibility.
#' @param df A data frame containing the longitudinal covariate information with:
#'   \itemize{
#'     \item \code{ID}: individual identifier.
#'     \item \code{time}: observation times for covariate measurements.
#'     \item Covariate columns named with prefix \code{"des"}
#'           (e.g., \code{des1}, \code{des2}, ...).
#'   }
#'   Covariates are assumed to be piecewise constant between observation times.
#' @param chfun A function computing the baseline cumulative hazard,
#'   e.g. \code{chfun(time, ...)}.
#' @param hstr Character string specifying the hazard structure.
#'   Either \code{"PH"} for proportional hazards or \code{"AFT"} for
#'   accelerated failure time models.
#' @param theta Numeric vector of baseline hazard parameters. The length of
#'   \code{theta} determines whether a two- or three-parameter baseline is used.
#' @param beta Numeric vector of regression coefficients corresponding to the
#'   time-varying covariates.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{time}: numeric vector of simulated event or censoring times.
#'     \item \code{status}: event indicator (1 = event, 0 = right-censored).
#'   }
#'
#' @details
#' The simulation assumes that covariates are left-continuous and piecewise
#' constant over time. Event times are generated conditional on the observed
#' covariate history up to the last follow-up time for each individual.
#'
#' This function is suitable for simulation studies involving time-varying
#' covariates under parametric PH or AFT models.
#'
#' @seealso \code{\link{CH_TVC}}, \code{\link{SPred_TVC}}, \code{\link{SPred_TVC_i}}
#'
#' @export
sim_TVC <- function(seed, df, chfun, hstr, theta, beta) {

  set.seed(seed)

  n <- length(unique(df$ID))
  sim    <- numeric(n)
  status <- integer(n)

  ## Survival at last observed time

  survs <- SPred_TVC(
    df    = df,
    beta  = beta,
    theta = theta,
    chfun = chfun,
    hstr  = hstr
  )


  ## Maximum follow-up times
  times <- as.vector(with(df, tapply(time, ID, max)))

  ## Uniform draws
  u <- runif(n)

  ## Censoring indicator
  censored <- (survs > u)

  if(sum(censored) > 0){
    sim[censored]    <- times[censored]
    status[censored] <- 0
  }

  ## Event times
  for (i in which(!censored)) {

    rootfun <- function(t) {
      SPred_TVC_i(
        df    = df,
        i     = i,
        t     = t,
        beta  = beta,
        theta = theta,
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



#----------------------------------------------------------------------------------------
#' Maximum Likelihood Estimation for Parametric Hazard Models with Time-Varying Covariates
#'
#' @description
#' `HMLE_TVC()` fits parametric survival models in the presence of
#' **time-varying covariates**, using maximum likelihood estimation.
#'
#' The function supports:
#'
#' * **Proportional Hazards (PH)** models with time-varying covariates
#' * Fully parametric baseline hazards (2-parameter or 3-parameter)
#'
#' * **Accelerated Failure Time (AFT)** models with time-varying covariates
#' * Fully parametric baseline hazards (2-parameter or 3-parameter)
#'
#' The likelihood is constructed from the cumulative hazard differences across
#' observation intervals for each individual, using a counting-process representation.
#'
#' For each individual, the data must contain several rows:
#' one per time-varying covariate measurement, along with the corresponding time.
#'
#----------------------------------------------------------------------------------------
#' @param init    : initial point for optimisation step
#' under the parameterisation (log(scale), log(shape1), log(shape2), beta) for scale-shape1-shape2 models or
#' (mu, log(scale), beta) for log-location scale models.
#----------------------------------------------------------------------------------------
#' @param df
#' A data frame in **long format**, containing one row per individual per
#' covariate-measurement time. Required columns:
#'
#' * `ID` — individual identifier
#' * `time` — time at which the covariates are measured
#' * `des*` — covariate columns used in the model (e.g., `des1`, `des2`, …)
#'
#' The last row for each ID represents the individual's event/censoring time,
#' even if the event time does not coincide with a measurement time.
#'
#----------------------------------------------------------------------------------------
#' @param status vector of event indicators (1 = event at the final time; 0 = censored)
#'
#----------------------------------------------------------------------------------------
#' @param hstr Hazard structure ("PH" or "AFT")
#----------------------------------------------------------------------------------------
#' @param dist    : distribution for the baseline hazard:
#'                 Power Generalised Weibull ("PGW")
#'                 Generalised Gamma ("GenGamma"))
#'                 Exponentiated Weibull ("EW")
#'                 Weibull ("Weibull")
#'                 Gamma ("Gamma")
#'                 LogNormal ("LogNormal")
#'                 LogLogistic ("LogLogistic")
#'
#----------------------------------------------------------------------------------------
#' @param method
#' Optimisation method for the likelihood.
#' Either `"nlminb"` or a valid `optim()` method.
#'
#' @param maxit
#' Maximum number of optimisation iterations.
#'
#----------------------------------------------------------------------------------------
#' @return
#' A list containing:
#'
#' * The full output from `optim()` or `nlminb()`
#' * The **negative log-likelihood function** used for optimisation
#' * A vector giving, for each ID, the cumulative hazard increments used in the likelihood
#'
#' Returned invisibly where appropriate.
#'
#----------------------------------------------------------------------------------------
#' @details
#'
#' ## Likelihood formulation for PH models
#'
#' For each individual \(i\), let
#' \(t_{i1} < t_{i2} < \cdots < t_{iK_i}\)
#' denote the *observation / measurement times*.
#'
#' The cumulative hazard contribution over interval \((t_{ij-1}, t_{ij})\) is:
#'
#' \deqn{
#' \Delta H_{ij}
#'   = \left[ H_0(t_{ij}) - H_0(t_{ij-1}) \right]
#'     \exp(x_{ij}^\top \beta),
#' }
#'
#' where \(x_{ij}\) is the vector of covariates measured at time \(t_{ij}\).
#'
#' The full log-likelihood is:
#'
#' \deqn{
#' \ell = \sum_i \left(
#'   - \sum_j \Delta H_{ij}
#'   + \delta_i \log \left[
#'        h_0(T_i) \exp(x_{iK}^\top \beta)
#'     \right]
#' \right),
#' }
#'
#' where:
#'
#' * \(\Delta H_{ij}\) comes from cumulative hazard increments
#' * \(T_i = t_{iK}\) is the final event or censoring time
#' * \(\delta_i\) is the event indicator
#' * hazard and cumulative hazard are computed based on `dist`
#'
#' The function internally:
#' 1. Splits the data by ID
#' 2. Computes cumulative hazard at all measurement times
#' 3. Computes increments \(\Delta H_{ij}\) for each ID
#' 4. Constructs the likelihood
#' 5. Optimises over \(\beta\) and baseline parameters
#'
#----------------------------------------------------------------------------------------
#' @section Data structure:
#' The input data frame must contain:
#'
#' * varying number of rows per ID
#' * strictly increasing `time` within each ID
#' * last row containing the event/censoring time
#'
#' Covariates must be named as `des1`, `des2`, etc.
#'
#' @export
HMLE_TVC <-  function (init,
            df,
            status,
            hstr = NULL,
            dist = NULL,
            method = "Nelder-Mead",
            maxit = 100)
  {
    df <- df[order(df$ID, df$time),]  # ensure sorted

    times <- as.vector(with(df, tapply(time, ID, max)))

    last_rows <- df[ave(df$time, df$ID, FUN = max) == df$time,]
    des <- as.matrix(last_rows[, grep("^des", names(df))])


    status <- as.vector(as.logical(status))
    times.obs <- times[status]
    des_obs <- des[status,]


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
          theta = c(ae0,be0,ce0)
          beta <- par[4:(3 + p)]
          # Hazard calculations
          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hpgw(times.obs, ae0, be0, ce0, log = TRUE) +
            x.beta.obs

          # Cumulative hazard calculations
          df$CH_t <- CH_TVC(
            df = df,
            beta = beta,
            theta = theta,
            chfun = chpgw,
            hstr = "PH"
          )$cum_hazard

          # split by ID
          lst <- split(df, df$ID)

          # compute lagged differences within each ID
          res_list <- lapply(lst, function(d) {
            data.frame(
              ID = d$ID[-1],
              # drop the first (no diff)
              time = d$time[-1],
              # time of the difference
              diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
            )
          })


          # sum diff_chaz for each ID
          chaz0 <-
            as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

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
          theta = c(ae0,be0,ce0)
          beta <- par[4:(3 + p)]
          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hew(times.obs, ae0, be0, ce0, log = TRUE) +
            x.beta.obs

          # Cumulative hazard calculations
          df$CH_t <- CH_TVC(
            df = df,
            beta = beta,
            theta = theta,
            chfun = chew,
            hstr = "PH"
          )$cum_hazard

          # split by ID
          lst <- split(df, df$ID)

          # compute lagged differences within each ID
          res_list <- lapply(lst, function(d) {
            data.frame(
              ID = d$ID[-1],
              # drop the first (no diff)
              time = d$time[-1],
              # time of the difference
              diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
            )
          })


          # sum diff_chaz for each ID
          chaz0 <-
            as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

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
          theta = c(ae0,be0,ce0)
          beta <- par[4:(3 + p)]
          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hggamma(times.obs, ae0, be0, ce0, log = TRUE) +
            x.beta.obs

          # Cumulative hazard calculations
          df$CH_t <- CH_TVC(
            df = df,
            beta = beta,
            theta = theta,
            chfun = chggama,
            hstr = "PH"
          )$cum_hazard

          # split by ID
          lst <- split(df, df$ID)

          # compute lagged differences within each ID
          res_list <- lapply(lst, function(d) {
            data.frame(
              ID = d$ID[-1],
              # drop the first (no diff)
              time = d$time[-1],
              # time of the difference
              diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
            )
          })


          # sum diff_chaz for each ID
          chaz0 <-
            as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

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
          theta = c(ae0,be0)
          beta <- par[3:(2 + p)]
          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hlnorm(times.obs, ae0, be0, log = TRUE) +
            x.beta.obs

          # Cumulative hazard calculations
          df$CH_t <- CH_TVC(
            df = df,
            beta = beta,
            ae0 = ae0,
            be0 = be0,
            chfun = chlnorm,
            hstr = "PH"
          )$cum_hazard

          # split by ID
          lst <- split(df, df$ID)

          # compute lagged differences within each ID
          res_list <- lapply(lst, function(d) {
            data.frame(
              ID = d$ID[-1],
              # drop the first (no diff)
              time = d$time[-1],
              # time of the difference
              diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
            )
          })


          # sum diff_chaz for each ID
          chaz0 <-
            as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

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
          theta = c(ae0,be0)
          beta <- par[3:(2 + p)]
          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hllogis(times.obs, ae0, be0, log = TRUE) +
            x.beta.obs

          # Cumulative hazard calculations
          df$CH_t <- CH_TVC(
            df = df,
            beta = beta,
            theta = theta,
            chfun = chllogis,
            hstr = "PH"
          )$cum_hazard

          # split by ID
          lst <- split(df, df$ID)

          # compute lagged differences within each ID
          res_list <- lapply(lst, function(d) {
            data.frame(
              ID = d$ID[-1],
              # drop the first (no diff)
              time = d$time[-1],
              # time of the difference
              diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
            )
          })


          # sum diff_chaz for each ID
          chaz0 <-
            as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

          # Negative log-likelihood
          val <- -sum(lhaz0) + sum(chaz0)


        }
      }
      if (dist == "Gamma") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          theta = c(ae0,be0)
          beta <- par[3:(2 + p)]
          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hgamma(times.obs, ae0, be0, log = TRUE) +
            x.beta.obs

          # Cumulative hazard calculations
          df$CH_t <- CH_TVC(
            df = df,
            beta = beta,
            theta = theta,
            chfun = chgamma,
            hstr = "PH"
          )$cum_hazard

          # split by ID
          lst <- split(df, df$ID)

          # compute lagged differences within each ID
          res_list <- lapply(lst, function(d) {
            data.frame(
              ID = d$ID[-1],
              # drop the first (no diff)
              time = d$time[-1],
              # time of the difference
              diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
            )
          })


          # sum diff_chaz for each ID
          chaz0 <-
            as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

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
          theta = c(ae0,be0)
          beta <- par[3:(2 + p)]
          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hweibull(times.obs, ae0, be0, log = TRUE) +
            x.beta.obs

          # Cumulative hazard calculations
          df$CH_t <- CH_TVC(
            df = df,
            beta = beta,
            theta = theta,
            chfun = chweibull,
            hstr = "PH"
          )$cum_hazard

          # split by ID
          lst <- split(df, df$ID)

          # compute lagged differences within each ID
          res_list <- lapply(lst, function(d) {
            data.frame(
              ID = d$ID[-1],
              # drop the first (no diff)
              time = d$time[-1],
              # time of the difference
              diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
            )
          })


          # sum diff_chaz for each ID
          chaz0 <-
            as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

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
          theta = c(ae0,be0,ce0)
          beta <- par[4:(3 + p)]
          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hpgw(times.obs * exp.x.beta.obs, ae0,
                        be0, ce0, log = TRUE) + x.beta.obs

          # Cumulative hazard calculations
          df$CH_t <- CH_TVC(
            df = df,
            beta = beta,
            theta = theta,
            chfun = chpgw,
            hstr = "AFT"
          )$cum_hazard

          # split by ID
          lst <- split(df, df$ID)

          # compute lagged differences within each ID
          res_list <- lapply(lst, function(d) {
            data.frame(
              ID = d$ID[-1],
              # drop the first (no diff)
              time = d$time[-1],
              # time of the difference
              diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
            )
          })


          # sum diff_chaz for each ID
          chaz0 <-
            as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

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
          theta = c(ae0,be0,ce0)
          beta <- par[4:(3 + p)]
          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hew(times.obs * exp.x.beta.obs, ae0,
                       be0, ce0, log = TRUE) + x.beta.obs

          # Cumulative hazard calculations
          df$CH_t <- CH_TVC(
            df = df,
            beta = beta,
            theta = theta,
            chfun = chew,
            hstr = "AFT"
          )$cum_hazard

          # split by ID
          lst <- split(df, df$ID)

          # compute lagged differences within each ID
          res_list <- lapply(lst, function(d) {
            data.frame(
              ID = d$ID[-1],
              # drop the first (no diff)
              time = d$time[-1],
              # time of the difference
              diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
            )
          })


          # sum diff_chaz for each ID
          chaz0 <-
            as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

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
          theta = c(ae0,be0,ce0)
          beta <- par[4:(3 + p)]
          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hggamma(times.obs * exp.x.beta.obs,
                           ae0, be0, ce0, log = TRUE) + x.beta.obs

          # Cumulative hazard calculations
          df$CH_t <- CH_TVC(
            df = df,
            beta = beta,
            theta = theta,
            chfun = chggamma,
            hstr = "AFT"
          )$cum_hazard

          # split by ID
          lst <- split(df, df$ID)

          # compute lagged differences within each ID
          res_list <- lapply(lst, function(d) {
            data.frame(
              ID = d$ID[-1],
              # drop the first (no diff)
              time = d$time[-1],
              # time of the difference
              diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
            )
          })


          # sum diff_chaz for each ID
          chaz0 <-
            as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

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
          theta = c(ae0,be0)
          beta <- par[3:(2 + p)]
          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hlnorm(times.obs * exp.x.beta.obs, ae0,
                          be0, log = TRUE) + x.beta.obs

          # Cumulative hazard calculations
          df$CH_t <- CH_TVC(
            df = df,
            beta = beta,
            theta = theta,
            chfun = chlnorm,
            hstr = "AFT"
          )$cum_hazard

          # split by ID
          lst <- split(df, df$ID)

          # compute lagged differences within each ID
          res_list <- lapply(lst, function(d) {
            data.frame(
              ID = d$ID[-1],
              # drop the first (no diff)
              time = d$time[-1],
              # time of the difference
              diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
            )
          })


          # sum diff_chaz for each ID
          chaz0 <-
            as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

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
          theta = c(ae0,be0)
          beta <- par[3:(2 + p)]
          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hllogis(times.obs * exp.x.beta.obs,
                           ae0, be0, log = TRUE) + x.beta.obs

          # Cumulative hazard calculations
          df$CH_t <- CH_TVC(
            df = df,
            beta = beta,
            theta = theta,
            chfun = chllogis,
            hstr = "AFT"
          )$cum_hazard

          # split by ID
          lst <- split(df, df$ID)

          # compute lagged differences within each ID
          res_list <- lapply(lst, function(d) {
            data.frame(
              ID = d$ID[-1],
              # drop the first (no diff)
              time = d$time[-1],
              # time of the difference
              diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
            )
          })


          # sum diff_chaz for each ID
          chaz0 <-
            as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

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
          theta = c(ae0,be0)
          beta <- par[3:(2 + p)]
          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hgamma(times.obs * exp.x.beta.obs, ae0,
                          be0, log = TRUE) + x.beta.obs

          # Cumulative hazard calculations
          df$CH_t <- CH_TVC(
            df = df,
            beta = beta,
            theta = theta,
            chfun = chgamma,
            hstr = "AFT"
          )$cum_hazard

          # split by ID
          lst <- split(df, df$ID)

          # compute lagged differences within each ID
          res_list <- lapply(lst, function(d) {
            data.frame(
              ID = d$ID[-1],
              # drop the first (no diff)
              time = d$time[-1],
              # time of the difference
              diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
            )
          })


          # sum diff_chaz for each ID
          chaz0 <-
            as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

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
          theta = c(ae0,be0)
          beta <- par[3:(2 + p)]
          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hweibull(times.obs * exp.x.beta.obs,
                            ae0, be0, log = TRUE) + x.beta.obs

          # Cumulative hazard calculations
          df$CH_t <- CH_TVC(
            df = df,
            beta = beta,
            theta = theta,
            chfun = chweibull,
            hstr = "AFT"
          )$cum_hazard

          # split by ID
          lst <- split(df, df$ID)

          # compute lagged differences within each ID
          res_list <- lapply(lst, function(d) {
            data.frame(
              ID = d$ID[-1],
              # drop the first (no diff)
              time = d$time[-1],
              # time of the difference
              diff_chaz = diff(d$CH_t)    # value[t] - value[t-1]
            )
          })


          # sum diff_chaz for each ID
          chaz0 <-
            as.vector(tapply(df$CH_t, df$ID, sum, na.rm = TRUE))

          # Negative log-likelihood
          val <- -sum(lhaz0) + sum(chaz0)

          return(val)
        }
      }
    }
    # Optimisation step
    if (method != "nlminb") {
      OPT <- optim(init,
                   log.lik,
                   control = list(maxit = maxit),
                   method = method)
    }
    if (method == "nlminb") {
      OPT <- nlminb(init, log.lik, control = list(iter.max = maxit))
    }
    OUT <- list(log_lik = log.lik, OPT = OPT)
    return(OUT)
  }


