

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
#' @param df A data frame in longitudinal format containing:
#'   \itemize{
#'     \item `ID`: an identifier for each patient, potentially across multiple time points. 
#'     \item `time`: numeric vector of time points, strictly monotonically increasing within each ID. 
#'   }
#' @param beta Numeric vector of regression coefficients.
#' @param theta Numeric baseline parameters of the cumulative hazard.
#' @param chfun A function computing the baseline cumulative hazard:
#'   `chfun(time, theta[1], theta[2])` or `chfun(time, theta[1], theta[2], theta[3])`.
#' @param hstr Hazard structure ("PH" or "AFT")
#' @param `ID_nm` - The name of the `ID` variable in `df`. Default "ID".
#' @param `time_nm` - The name of the `time` variable in `df`. Default "time". 
#' @return A numeric vector with the cumulative hazard evaluated at each time
#'   point in `df`.
#'
#' @export
CH_TVC <- function(df, beta,
                   theta,
                   chfun, hstr, ID_nm = "ID", time_nm = "time"){
  npar = length(theta)
  ##Convert name into des* structure.
  colnames(df)[!(colnames(df) %in% c(ID_nm, time_nm) ) ] = 
    sapply(FUN = function(x) paste("des", x, sep = ""), X = as.character(seq_len((ncol(df)-2) ) ) )
  # Original order
  df = cbind(df,original_order_idx = seq_len(nrow(df)))
  print(df)
  ## Ensure sorted data 
  df <- df[order(df[,"original_order_idx"], df[, time_nm]), ]
  ## Extract design matrix and associated linear predictor for every row.
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
  split_df <- split(seq_len(nrow(df)), df[, ID_nm])
  ## Storage
  H_out <- numeric(nrow(df))
  #for every single unique ID:
  for (idx in split_df) {
    #get times
    t <- df[idx, time_nm]
    #get associated factors
    xb <- exp_xb[idx]
    #get delta t's by taking the difference between t and t shifted down by 1 index. 
    dt = c(0, diff(t))
    
    if (hstr == "PH") {

      H0_t <- H0(t)
      dH0 <- H0(dt)
      H_i  <- cumsum((H0_t - dH0) * xb)

    }

    if (hstr == "AFT") {
      #evaluate times and differences between times, then take their sum. 
      H0_t <- H0(t * xb)
      print(H0_t)
      dH0 <- H0(dt * xb)
      print(dH0)
      H_i  <- cumsum(H0_t - dH0)
    }
    H_out[idx] <- H_i
  }
  
  df = cbind(df, cum_hazard = H_out)

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
#' @param df A data frame in longitudinal format containing:
#'   \itemize{
#'     \item `ID`: an identifier for each patient, potentially across multiple time points. 
#'     \item `time`: numeric vector of time points, strictly monotonically increasing within each ID. 
#'   }
#' @param beta Numeric vector of regression coefficients.
#' @param theta Numeric baseline parameters of the cumulative hazard.
#' @param chfun A function computing the baseline cumulative hazard:
#'   `chfun(time, theta[1], theta[2])` or `chfun(time, theta[1], theta[2], theta[3])`.
#' @param hstr Hazard structure ("PH" or "AFT")
#' @param `ID_nm` - The name of the `ID` variable in `df`. Default "ID".
#' @param `time_nm` - The name of the `time` variable in `df`. Default "time". 
#'
#' @return A numeric vector with the survival function evaluated at the last time
#'   point in `df`.
#'
#' @export
SPred_TVC <-
  function(df, beta, theta, chfun, hstr, ID_nm = "ID", time_nm = "time") {
    # Sample size
    n <- max(df[, ID_nm])
    # Calculating the cumulative hazard function at all time points
    CH <- CH_TVC(
      df    = df,
      beta  = beta,
      theta = theta,
      chfun = chfun,
      hstr = hstr, 
      ID_nm = ID_nm, 
      time_nm = time_nm
    )
    #Find an individual's cumulative hazard as the sum of their CHs for each 
    #time split. 
    H_last = aggregate(CH, by = formula(CH$cum_hazard ~ CH[, ID_nm]), FUN = sum)[, 2]

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
#'
#' @param df A data frame in longitudinal format containing:
#'   \itemize{
#'     \item `ID`: an identifier for each patient, potentially across multiple time points. 
#'     \item `time`: numeric vector of time points, strictly monotonically increasing within each ID. 
#'   }
#' @param beta Numeric vector of regression coefficients.
#' @param theta Numeric baseline parameters of the cumulative hazard.
#' @param chfun A function computing the baseline cumulative hazard:
#'   `chfun(time, theta[1], theta[2])` or `chfun(time, theta[1], theta[2], theta[3])`.
#' @param hstr Hazard structure ("PH" for Proportional Hazards or "AFT" for Accelerated Failure Time)
#' @param `ID_nm` - The name of the `ID` variable in `df`. Default "ID".
#' @param `time_nm` - The name of the `time` variable in `df`. Default "time". 
#' @param i Integer specifying the individual for whom the survival
#'   function is to be evaluated.
#' @param t Numeric value giving the time point at which the survival
#'   function is evaluated. Must lie within the observation window of
#'   individual `i`.
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
    hstr, 
    ID_nm = "ID", 
    time_nm = "time"
) {
  ## 1. Subset individual data
  dfi <- df[df[, ID_nm] == i, ]
  dfi <- dfi[order(dfi[, time_nm]), ]
  
  
  ## 2. Add time t if needed (piecewise-constant covariates, even from time 0)
  #if the time lies before the first observation, we copy all the covariates 
  #from the first time point backwards, and use those to construct a new dataframe
  #for time t < min(dfi[, time]). Note under the case with t < min(dfi[, time_nm])
  #there will be a warning, but this is a superfluous warning and can be ignored. 
  
  #If t exceeds the last observation time, then this function will extrapolate
  #(with questionable certainty), via the piecewise-constant assumption. To do 
  #so, we just set the time at the last row to be t. 
  
  if (t >= max(dfi[, time_nm])){
    pre_obs = 0
    dfi[nrow(dfi), time_nm] = t
  }
  
  else if (t <= min(dfi[, time_nm])) {
    pre_obs = 1
    dfi = dfi[1, ]
    dfi[, time_nm] = t
  }
  
  else{
    pre_obs = 0
    idx <- max(which(dfi[, time_nm] < t))
    newrow <- dfi[idx, ]
    newrow[, time_nm] <- t
    dfi <- rbind(dfi, newrow)
    dfi <- dfi[order(dfi[, time_nm]), ]
  }

  ## 3. Compute cumulative hazard via existing engine
  CH <- CH_TVC(
    df    = dfi,
    beta  = beta,
    theta = theta,
    chfun = chfun,
    hstr  = hstr, 
    ID_nm = ID_nm,
    time_nm = time_nm
  )
  print(CH)
  ## 4.  Extract cumulative hazard at time t, as the SUM OF THE CUMULATIVE 
  #HAZARDS UP TO POINT T. 
  assign(x = "H_i_t", value = ifelse(pre_obs, yes = CH$cum_hazard, 
                                     no = sum(CH$cum_hazard[1:which.min(abs(dfi[,time_nm] - t))])))
  

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
#'     \item `ID`: an identifier for each patient, potentially across multiple time points. 
#'     \item `time`: numeric vector of time points, strictly monotonically increasing within each ID.
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
sim_TVC <- function(n = NULL, seed, df, chfun, hstr, theta, beta, ID_nm = "ID", 
                    time_nm = "time"){

  set.seed(seed)
  
  if (is.null(n)) n <- length(unique(df[, ID_nm]))
  sim    <- rep(NA, n)

  ## Maximum follow-up times
  times <- unlist(lapply(X = split(df, ~ df[,ID_nm]), FUN = function(x) max(x[, time_nm])),
                  use.names = F)

  ## Uniform draws
  u <- runif(n)
  
  ## Event times
  for (i in seq_len(n)) {
    ##evaluate the analytical expression for t using uniroot
    rootfun <- function(t) {
      SPred_TVC_i(
        df    = df,
        i     = i,
        t     = t,
        beta  = beta,
        theta = theta,
        chfun = chfun,
        hstr  = hstr,
        ID_nm = ID_nm,
        time_nm = time_nm
      ) - u[i]
    }
    #find t which solves S(t) = u
    sim[i] <- uniroot(rootfun, interval = c(0, times[i]), extendInt = "downX", tol = 1e-8)$root
    }

  return(sim)
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


