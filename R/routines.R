#----------------------------------------------------------------------------------------
#' Power Generalised Weibull (PGW) probability density function.
#' http://rpubs.com/FJRubio/PGW
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @param t       : positive argument
#' @param log: log scale (TRUE or FALSE)
#' @return the value of the PGW probability density function
#' @export

dpgw <- function(t, sigma, nu, gamma, log = FALSE){
  val <- log(nu) - log(gamma) - nu*log(sigma) + (nu-1)*log(t) +
    (1/gamma - 1)*log( 1 + (t/sigma)^nu ) +
    ( 1 - ( 1 + (t/sigma)^nu )^(1/gamma) )
  if(log) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Power Generalised Weibull (PGW) survival function.
#' http://rpubs.com/FJRubio/PGW
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @param t       : positive argument
#' @param log.p: log scale (TRUE or FALSE)
#' @return the value of the PGW survival function
#' @export
spgw <- function(t, sigma, nu, gamma, log.p = FALSE){
  val <- 1 - ( 1 + (t/sigma)^nu )^(1/gamma)
  if(log.p) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Power Generalised Weibull (PGW) hazard function.
#' http://rpubs.com/FJRubio/PGW
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @param t       : positive argument
#' @param log: log scale (TRUE or FALSE)
#' @return the value of the PGW hazard function
#' @export
hpgw <- function(t, sigma, nu, gamma, log = FALSE){
  val <- log(nu) - log(gamma) - nu*log(sigma) + (nu-1)*log(t) +
    (1/gamma - 1)*log( 1 + (t/sigma)^nu )
  if(log) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Power Generalised Weibull (PGW) cumulative hazard function.
#' http://rpubs.com/FJRubio/PGW
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @param t       : positive argument
#' @return the value of the PGW cumulative hazard function
#' @export
chpgw <- function(t, sigma, nu, gamma){
  val <- -1 + ( 1 + (t/sigma)^nu )^(1/gamma)
  return(val)
}

#----------------------------------------------------------------------------------------
#' Power Generalised Weibull (PGW) random number generation.
#' http://rpubs.com/FJRubio/PGW
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @param n       : 	number of observations
#' @return generates random deviates
#' @export
rpgw <- function(n, sigma, nu, gamma){
  p <- runif(n)
  out <- sigma*(  ( 1 - log(1-p) )^gamma - 1 )^(1/nu)
  return(as.vector(out))
}

#----------------------------------------------------------------------------------------
#' Power Generalised Weibull (PGW) quantile function.
#' http://rpubs.com/FJRubio/PGW
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @param p       : 	probability. A value in (0,1)
#' @return the value of the PGW quantile function
#' @export
qpgw <- function(p, sigma, nu, gamma){
  out <- sigma*(  ( 1 - log(1-p) )^gamma - 1 )^(1/nu)
  return(out)
}


#----------------------------------------------------------------------------------------
#' Power Exponentiated Weibull (EW) hazard function.
#' https://rpubs.com/FJRubio/EWD
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @param t       : positive argument
#' @param log: log scale (TRUE or FALSE)
#' @return the value of the EW hazard function
#' @export
hew <- function(t, sigma, nu, gamma, log=FALSE){
  log.pdf <-  log(gamma) + (gamma-1)*pweibull(t,scale=sigma,shape=nu,log.p=TRUE) +
    dweibull(t,scale=sigma,shape=nu,log=TRUE)
  cdf <- exp(gamma*pweibull(t,scale=sigma,shape=nu,log.p=TRUE) )
  log.h <- log.pdf - log(1-cdf)
  ifelse(log, return(log.h), return(exp(log.h)))
}

#----------------------------------------------------------------------------------------
#' Power Exponentiated Weibull (EW) cumulative hazard function.
#' https://rpubs.com/FJRubio/EWD
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @param t       : positive argument
#' @return the value of the EW cumulative hazard function
#' @export
chew <- function(t, sigma, nu, gamma){
  cdf <- exp(gamma*pweibull(t,scale=sigma,shape=nu,log.p=TRUE) )
  return(-log(1-cdf))
}

#----------------------------------------------------------------------------------------
#' Power Exponentiated Weibull (EW) quantile function.
#' https://rpubs.com/FJRubio/EWD
#----------------------------------------------------------------------------------------
#' @param p       : 	probability. A value in (0,1)
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @return the value of the EW quantile function
#' @export
qew <- function(p, sigma, nu, gamma){
  quant <-  qweibull(p^(1/gamma),scale=sigma,shape=nu)
  return(quant)
}

#----------------------------------------------------------------------------------------
#' Weibull (W) hazard function.
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param t       : positive argument
#' @param log: log scale (TRUE or FALSE)
#' @return the value of the Weibull hazard function
#' @export
hweibull <- function(t, sigma, nu, log=FALSE){
  log.pdf <-  dweibull(t,scale=sigma,shape=nu, log = TRUE)
  log.s <- pweibull(t,scale=sigma,shape=nu, lower.tail = FALSE, log.p = TRUE)
  log.h <- log.pdf - log.s
  ifelse(log, return(log.h), return(exp(log.h)))
}


#----------------------------------------------------------------------------------------
#'  Weibull (W) cumulative hazard function.
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param t       : positive argument
#' @return the value of the Weibull cumulative hazard function
#' @export
chweibull <- function(t, sigma, nu){
  H0 <- -pweibull(t,scale=sigma,shape=nu, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}

#----------------------------------------------------------------------------------------
#' Lognormal (LN) hazard function.
#----------------------------------------------------------------------------------------
#' @param mu      : mean parameter in the log scale
#' @param sigma     : scale parameter in the log scale
#' @param t       : positive argument
#' @param log: log scale (TRUE or FALSE)
#' @return the value of the LN hazard function
#' @export

hlnorm <- function(t, mu, sigma, log = FALSE){
  lpdf0 <-  dlnorm(t, mu, sigma, log = TRUE)
  ls0 <- plnorm(t, mu, sigma, lower.tail = FALSE, log.p = TRUE)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Lognormal (LN) cumulative hazard function.
#----------------------------------------------------------------------------------------
#' @param mu      : mean parameter in the log scale
#' @param sigma     : scale parameter in the log scale
#' @param t       : positive argument
#' @return the value of the LN cumulative hazard function
#' @export
chlnorm <- function(t, mu, sigma){
  H0 <- -plnorm(t, mu, sigma, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}

#----------------------------------------------------------------------------------------
#' Log-logistic (LL) hazard function.
#----------------------------------------------------------------------------------------
#' @param mu      : mean parameter in the log scale
#' @param sigma     : scale parameter in the log scale
#' @param t       : positive argument
#' @param log: log scale (TRUE or FALSE)
#' @return the value of the LL hazard function
#' @export

hllogis <- function(t, mu, sigma, log = FALSE){
  lpdf0 <-  dlogis(log(t),mu,sigma, log = TRUE) - log(t)
  ls0 <- plogis(log(t),mu,sigma, lower.tail = FALSE, log.p = TRUE)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Log-logistic (LL) cumulative hazard function.
#----------------------------------------------------------------------------------------
#' @param mu      : mean parameter in the log scale
#' @param sigma     : scale parameter in the log scale
#' @param t       : positive argument
#' @return the value of the LL cumulative hazard function
#' @export
chllogis <- function(t, mu, sigma){
  H0 <- -plogis(log(t),mu,sigma, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}

#----------------------------------------------------------------------------------------
#' Log-logistic (LL) quantile function.
#----------------------------------------------------------------------------------------
#' @param mu      : mean parameter in the log scale
#' @param sigma     : scale parameter in the log scale
#' @param p       : 	probability. A value in (0,1)
#' @return the value of the LL quantile function
#' @export
qllogis <- function(p, mu, sigma){
  qq <- exp(qlogis(p, mu, sigma))
  return(qq)
}

#----------------------------------------------------------------------------------------
#'  Gamma (G) hazard function.
#----------------------------------------------------------------------------------------
#' @param scale     : scale parameter
#' @param shape      : shape parameter
#' @param t       : positive argument
#' @param log: log scale (TRUE or FALSE)
#' @return the value of the Gamma hazard function
#' @export
hgamma <- function(t, shape, scale, log = FALSE){
  lpdf0 <-  dgamma(t, shape = shape, scale = scale, log = TRUE)
  ls0 <- pgamma(t, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#'  Gamma (G) cumulative hazard function.
#----------------------------------------------------------------------------------------
#' @param scale     : scale parameter
#' @param shape      : shape parameter
#' @param t       : positive argument
#' @param log: log scale (TRUE or FALSE)
#' @return the value of the Weibull hazard function
#' @export
chgamma <- function(t, shape, scale){
  H0 <- -pgamma(t, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}

#----------------------------------------------------------------------------------------
#' Generalised Gamma (GG) probability density function.
#' https://rpubs.com/FJRubio/GG
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @param t       : positive argument
#' @param log: log scale (TRUE or FALSE)
#' @return the value of the GG probability density function
#' @export

dggamma <- function(t, sigma, nu, gamma, log = FALSE){
  val <- log(gamma) - nu*log(sigma) - lgamma(nu/gamma) + (nu - 1)*log(t) -
    (t/sigma)^gamma
  if(log) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Generalised Gamma (GG) cumulative distribution function.
#' https://rpubs.com/FJRubio/GG
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @param t       : positive argument
#' @param log.p: log scale (TRUE or FALSE)
#' @return the value of the GG cumulative distribution function
#' @export
pggamma <- function(t, sigma, nu, gamma, log.p = FALSE){
  val <- pgamma( t^gamma, shape = nu/gamma, scale = sigma^gamma, log.p = TRUE)
  if(log.p) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Generalised Gamma (GG) survival function.
#' https://rpubs.com/FJRubio/GG
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @param t       : positive argument
#' @param log.p: log scale (TRUE or FALSE)
#' @return the value of the GG survival function
#' @export
sggamma <- function(t, sigma, nu, gamma, log.p = FALSE){
  val <- pgamma( t^gamma, shape = nu/gamma, scale = sigma^gamma, log.p = TRUE, lower.tail =  FALSE)
  if(log.p) return(val) else return(exp(val))
}


#----------------------------------------------------------------------------------------
#' Generalised Gamma (GG) hazard function.
#' https://rpubs.com/FJRubio/GG
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @param t       : positive argument
#' @param log: log scale (TRUE or FALSE)
#' @return the value of the GG hazard function
#' @export
hggamma <- function(t, sigma, nu, gamma, log = FALSE){
  log.pdf <- dggamma(t, sigma, nu, gamma, log = TRUE)
  log.s <- sggamma(t, sigma, nu, gamma, log.p = TRUE)
  val <- log.pdf - log.s
  if(log) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Generalised Gamma (GG) cumulative hazard function.
#' https://rpubs.com/FJRubio/GG
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @param t       : positive argument
#' @return the value of the GG cumulative hazard function
#' @export
chggamma <- function(t, sigma, nu, gamma){
  val <- -pgamma( t^gamma, shape = nu/gamma, scale = sigma^gamma, log.p = TRUE, lower.tail =  FALSE)
  return(val)
}

#----------------------------------------------------------------------------------------
#' Generalised Gamma (GG) quantile function.
#' https://rpubs.com/FJRubio/GG
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @param p       : 	probability. A value in (0,1)
#' @return the value of the GG quantile function
#' @export
qggamma <- function(p, sigma, nu, gamma){
  out <- qgamma(p, shape = nu/gamma, scale = sigma^gamma)^(1/gamma)
  return(out)
}

#----------------------------------------------------------------------------------------
#' Generalised Gamma (GG) random number generation.
#' https://rpubs.com/FJRubio/GG
#----------------------------------------------------------------------------------------
#' @param sigma     : scale parameter
#' @param nu      : shape parameter
#' @param gamma   : shape parameter
#' @param n       : 	number of observations
#' @return generates random deviates
#' @export
rggamma <- function(n, sigma, nu, gamma){
  p <- runif(n)
  out <- qgamma(p, shape = nu/gamma, scale = sigma^gamma)^(1/gamma)
  return(as.vector(out))
}



#----------------------------------------------------------------------------------------
#' GHMLE function: Hazard Regression Models with a parametric baseline hazard
#----------------------------------------------------------------------------------------
#' @param init    : initial point for optimisation step
#' under the parameterisation (log(scale), log(shape1), log(shape2), alpha, beta) for scale-shape1-shape2 models or
#' (mu, log(scale), alpha, beta) for log-location scale models.
#' @param hstr    : hazard structure:
#'            No covariates ("baseline"),
#'           AFT model with PGW baseline hazard ("AFT"),
#'           PH model with PGW baseline hazard ("PH"),
#'           AH model with PGW baseline hazard ("AH"),
#'           GH model with PGW baseline hazard ("GH")
#'           *GH is not available with Weibull dist
#' @param dist    : distribution for the baseline hazard:
#'                 Power Generalised Weibull ("PGW")
#'                 Generalised Gamma ("GenGamma"))
#'                 Exponentiated Weibull ("EW")
#'                 Weibull ("Weibull")
#'                 Gamma ("Gamma")
#'                 LogNormal ("LogNormal")
#'                 LogLogistic ("LogLogistic")
#' @param method  : "nlminb" or optimisation method to be used in optim (see ?optim)
#' @param maxit   : maximum number of iterations in optim or nlminb
#' @param times   : times to event
#' @param status    : vital status indicators (TRUE or 1 = observed, FALSE  or 0 = censored)
#' @param des     : design matrix for hazard-level effects
#' @param des_t   : design matrix for time-level effects (it is recommended not to use splines here)
#' @return It returns the output from optim or nlminb for the selected model and the negative log likelihood function
#' @export
GHMLE <-
  function(init,
           times,
           status,
           hstr = NULL,
           dist = NULL,
           des = NULL,
           des_t = NULL,
           method = "Nelder-Mead",
           maxit = 100) {
    # Required variables
    times <- as.vector(times)
    status <- as.vector(as.logical(status))
    times.obs <- times[status]
    if (!is.null(des)){
      des <- as.matrix(des)
      des.obs <- des[status, ]
    }
    if (!is.null(des_t)){
      des_t <- as.matrix(des_t)
      des_t.obs <- des_t[status, ]
    }

    #------------------------------------------------------------------------------------
    # baseline models
    #------------------------------------------------------------------------------------

    if (hstr == "baseline") {
      # PGW
      if (dist == "PGW") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])

          lhaz0 <- hpgw(times.obs, ae0, be0, ce0, log = TRUE)

          val <- -sum(lhaz0) + sum(chpgw(times, ae0, be0, ce0))
          return(val)
        }
      }

      # EW
      if (dist == "EW") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])

          lhaz0 <- hew(times.obs, ae0, be0, ce0, log = TRUE)

          val <- -sum(lhaz0) + sum(chew(times, ae0, be0, ce0))
          return(val)
        }
      }

      # GenGamma
      if (dist == "GenGamma") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])

          lhaz0 <- hggamma(times.obs, ae0, be0, ce0, log = TRUE)

          val <- -sum(lhaz0) + sum(chggamma(times, ae0, be0, ce0))
          return(val)
        }
      }

      # LogNormal
      if (dist == "LogNormal") {
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])

          lhaz0 <- hlnorm(times.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) + sum(chlnorm(times, ae0, be0))
          return(val)
        }
      }

      # LogLogistic
      if (dist == "LogLogistic") {
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])

          lhaz0 <- hllogis(times.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) + sum(chllogis(times, ae0, be0))
          return(val)
        }
      }

      # Gamma
      if (dist == "Gamma") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])

          lhaz0 <- hgamma(times.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) + sum(chgamma(times, ae0, be0))
          return(val)
        }
      }

      # Weibull
      if (dist == "Weibull") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])

          lhaz0 <- hweibull(times.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) + sum(chweibull(times, ae0, be0))
          return(val)
        }
      }

    }

    #------------------------------------------------------------------------------------
    # Proportional Hazards models
    #------------------------------------------------------------------------------------

    if (hstr == "PH") {

      # PGW
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
          lhaz0 <- hpgw(times.obs, ae0, be0, ce0, log = TRUE)  + x.beta.obs

          val <-  -sum(lhaz0) +
            sum(chpgw(times, ae0, be0, ce0) * exp.x.beta)
          return(val)
        }
      }

      # EW
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
          lhaz0 <- hew(times.obs, ae0, be0, ce0, log = TRUE)  + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chew(times, ae0, be0, ce0) * exp.x.beta)
          return(val)
        }
      }

      # GenGamma
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
          lhaz0 <- hggamma(times.obs, ae0, be0, ce0, log = TRUE)  + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chggamma(times, ae0, be0, ce0) * exp.x.beta)
          return(val)
        }
      }

      # LogNormal
      if (dist == "LogNormal") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hlnorm(times.obs, ae0, be0, log = TRUE)  + x.beta.obs

          val <- -sum(lhaz0) + sum(chlnorm(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }

      # LogLogistic
      if (dist == "LogLogistic") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hllogis(times.obs, ae0, be0, log = TRUE)  + x.beta.obs

          val <- -sum(lhaz0) + sum(chllogis(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }

      # Gamma
      if (dist == "Gamma") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hgamma(times.obs, ae0, be0, log = TRUE)  + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chgamma(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }

      # Weibull
      if (dist == "Weibull") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hweibull(times.obs, ae0, be0, log = TRUE)  + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chweibull(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }
    }

    #------------------------------------------------------------------------------------
    # Accelerated Failure Time models
    #------------------------------------------------------------------------------------

    if (hstr == "AFT") {

      # PGW
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
          lhaz0 <- hpgw(times.obs * exp.x.beta.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chpgw(times * exp.x.beta, ae0, be0, ce0))
          return(val)
        }
      }

      # EW
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
          lhaz0 <- hew(times.obs * exp.x.beta.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chew(times * exp.x.beta, ae0, be0, ce0))
          return(val)
        }
      }

      # GenGamma
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
          lhaz0 <- hggamma(times.obs * exp.x.beta.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs

          val <-  -sum(lhaz0) +
            sum(chggamma(times * exp.x.beta, ae0, be0, ce0))
          return(val)
        }
      }

      # LogNormal
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
          lhaz0 <- hlnorm(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) + x.beta.obs

          val <-  -sum(lhaz0) +
            sum(chlnorm(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }

      # LogLogistic
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
          lhaz0 <- hllogis(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chllogis(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }

      # Gamma
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
          lhaz0 <- hgamma(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chgamma(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }

      # Weibull
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
          lhaz0 <- hweibull(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chweibull(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }

    }


    #------------------------------------------------------------------------------------
    # Accelerated Hazards models
    #------------------------------------------------------------------------------------

    if (hstr == "AH") {

      # PGW
      if (dist == "PGW") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]
          lhaz0 <- hpgw(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE)

          val <- -sum(lhaz0) +
            sum(chpgw(times * exp.x.alpha, ae0, be0, ce0) / exp.x.alpha)
          return(val)
        }
      }

      # EW
      if (dist == "EW") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]
          lhaz0 <- hew(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE)

          val <- -sum(lhaz0) +
            sum(chew(times * exp.x.alpha, ae0, be0, ce0) / exp.x.alpha)
          return(val)
        }
      }

      # GenGamma
      if (dist == "GenGamma") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]
          lhaz0 <- hggamma(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE)

          val <- -sum(lhaz0) +
            sum(chggamma(times * exp.x.alpha, ae0, be0, ce0) / exp.x.alpha)
          return(val)
        }
      }

      # LogNormal
      if (dist == "LogNormal") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]
          lhaz0 <- hlnorm(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) +
            sum(chlnorm(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }

      # LogLogistic
      if (dist == "LogLogistic") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]
          lhaz0 <- hllogis(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) +
            sum(chllogis(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }

      # Gamma
      if (dist == "Gamma") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]
          lhaz0 <- hgamma(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) +
            sum(chgamma(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }

      # Weibull
      if (dist == "Weibull") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]
          lhaz0 <- hweibull(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) +
            sum(chweibull(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }

    }

    #------------------------------------------------------------------------------------
    # General Hazards models
    #------------------------------------------------------------------------------------

    if (hstr == "GH") {

      # PGW
      if (dist == "PGW") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p0)]
          beta <- par[(4 + p0):(3 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          lhaz0 <- hpgw(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs

          val <-  -sum(lhaz0) +
            sum(chpgw(times * exp.x.alpha, ae0, be0, ce0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      # EW
      if (dist == "EW") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p0)]
          beta <- par[(4 + p0):(3 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          lhaz0 <- hew(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs

          val <-  -sum(lhaz0) +
            sum(chew(times * exp.x.alpha, ae0, be0, ce0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      # GenGamma
      if (dist == "GenGamma") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p0)]
          beta <- par[(4 + p0):(3 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          lhaz0 <- hggamma(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs

          val <-  -sum(lhaz0) +
            sum(chggamma(times * exp.x.alpha, ae0, be0, ce0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      # LogNormal
      if (dist == "LogNormal") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p0)]
          beta <- par[(3 + p0):(2 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          lhaz0 <- hlnorm(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chlnorm(times * exp.x.alpha, ae0, be0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      # LogLogistic
      if (dist == "LogLogistic") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p0)]
          beta <- par[(3 + p0):(2 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          lhaz0 <- hllogis(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chllogis(times * exp.x.alpha, ae0, be0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      # Gamma
      if (dist == "Gamma") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p0)]
          beta <- par[(3 + p0):(2 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          lhaz0 <- hgamma(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chgamma(times * exp.x.alpha, ae0, be0) * exp.x.beta.dif)
          return(sum(val))
        }
      }
    }

    #------------------------------------------------------------------------------------
    # Optimisation step
    #------------------------------------------------------------------------------------

    if (method != "nlminb") {
      OPT <- optim(init,
                   log.lik,
                   control = list(maxit = maxit),
                   method = method)
    }
    if (method == "nlminb") {
      OPT <- nlminb(init, log.lik, control = list(iter.max = maxit))
    }

    #------------------------------------------------------------------------------------
    # Output
    #------------------------------------------------------------------------------------
    OUT <- list(log_lik = log.lik, OPT = OPT)
    return(OUT)
  }


#----------------------------------------------------------------------------------------
#' Function to calculate the normal confidence intervals.
#' The parameters indicated with "index" are transformed to the real line using log().
#----------------------------------------------------------------------------------------
#' @param FUN   : minus log-likelihood function to be used to calculate the confidence intervals
#' @param MLE   : maximum likelihood estimator of the parameters of interest
#' @param level : confidence level
#' @param index : position of the positive parameters under the original parameterisation
#' @return a list containing the upper and lower conf.int limits, the transformed MLE, and std errors
#' @export

Conf_Int <- function(FUN,MLE,level=0.95,index=NULL){
  sd.int <- abs(qnorm(0.5*(1-level)))
  tempf <- function(par){
    par[index] = exp(par[index])
    return(FUN( par ))
  }
  r.MLE <- MLE
  r.MLE[index] <- log(MLE[index])
  HESS <- hessian(tempf,x=r.MLE)
  Fisher.Info <- solve(HESS)
  Sigma <- sqrt(diag(Fisher.Info))
  U<- r.MLE + sd.int*Sigma
  L<- r.MLE - sd.int*Sigma
  C.I <- cbind(L,U,r.MLE, Sigma)
  names.row <- paste0("par", seq_along(1:length(MLE)))
  names.row[index] <- paste0("log.par", seq_along(index))
  rownames(C.I)<- names.row
  colnames(C.I)<- c("Lower","Upper","Transf MLE", "Std. Error")
  return(C.I)
}

#----------------------------------------------------------------------------------------
#' simGH function: Function to simulate times to event from a model with a GH structure
#' for different parametric baseline hazards.
#' Distributions: LogNormal, LogLogistic, GenGamma, Gamma, Weibull, PGW, EW.
#' See: https://github.com/FJRubio67/HazReg
#----------------------------------------------------------------------------------------
#' @param seed  : seed for simulation
#' @param n : sample size (number of individuals)
#' @param theta  :  parameters of the baseline hazard
#' @param beta_h  : regression parameters multiplying the hazard for GH model
#' @param beta_t  : regression parameters multiplying the time scale for GH model
#' @param beta  : regression parameters for AFT, PH, and AH models
#' @param des_h : Design matrix for GH model (hazard scale)
#' @param des_t : Design matrix for GH model (time scale)
#' @param des : Design matrix for AFT, PH, and AH models
#' @param hstr  : hazard structure (AH, AFT, PH, GH)
#' @param baseline  : baseline hazard distribution
#' @return a vector containing the simulated times to event
#' @export
#----------------------------------------------------------------------------------------
simGH <- function(seed, n, des = NULL, des_h = NULL, des_t = NULL, theta,
                  beta_h = NULL, beta_t = NULL, beta = NULL, hstr, baseline){

  if(!is.null(des))   des <- as.matrix(des)
  if(!is.null(des_h)) des_h <- as.matrix(des_h)
  if(!is.null(des_t)) des_t <- as.matrix(des_t)

  # Baseline hazard
  if(baseline == "LogNormal")      quantf <- function(p) qlnorm(p, theta[1], theta[2])
  if(baseline == "LogLogistic")      quantf <- function(p) qllogis(p, theta[1], theta[2])
  if(baseline == "Gamma")       quantf <- function(p) qgamma(p, theta[1], theta[2])
  if(baseline == "Weibull")       quantf <- function(p) qweibull(p, theta[1], theta[2])
  if(baseline == "PGW")     quantf <- function(p) qpgw(p, theta[1], theta[2], theta[3])
  if(baseline == "EW")      quantf <- function(p) qew(p, theta[1], theta[2], theta[3])
  if(baseline == "GenGamma")      quantf <- function(p) qggamma(p, theta[1], theta[2], theta[3])

  # Uniform variates used in the simulation
  set.seed(seed)
  u = runif(n)

  # GH simulation
  if( hstr == "GH" ){
    # Linear predictors
    exp.xbeta_t  <- exp(des_t%*%beta_t)
    exp.dif <- exp(des_t%*%beta_t - des_h%*%beta_h )

    # Simulating the times to event
    p0 <- as.vector(1 - exp(log(1-u)*exp.dif))
    times <- as.vector(quantf(p0)/exp.xbeta_t)
  }

  # PH simulation
  if( hstr == "PH" ){
    # Linear predictors
    exp.xbeta_t  <- 1
    exp.dif <- exp( - des%*%beta )

    # Simulating the times to event
    p0 <- as.vector(1 - exp(log(1-u)*exp.dif))
    times <- as.vector(quantf(p0)/exp.xbeta_t)
  }

  # AFT simulation
  if( hstr == "AFT" ){
    # Linear predictors
    exp.xbeta_t  <- exp(des%*%beta)
    exp.dif <- 1

    # Simulating the times to event
    p0 <- as.vector(1 - exp(log(1-u)*exp.dif))
    times <- as.vector(quantf(p0)/exp.xbeta_t)
  }

  # AH simulation
  if( hstr == "AH" ){
    # Linear predictors
    exp.xbeta_t  <- exp(des%*%beta)
    exp.dif <- exp(des%*%beta)

    # Simulating the times to event
    p0 <- as.vector(1 - exp(log(1-u)*exp.dif))
    times <- as.vector(quantf(p0)/exp.xbeta_t)
  }

  return(as.vector(times))
}




###############################################################################################
###############################################################################################
###############################################################################################
#' Relative (Net) Survival models
###############################################################################################
###############################################################################################
###############################################################################################

########################################################################################################
#' Log likelihood and MLE for the GH excess hazards model.
#' Baseline hazards: Lognormal, Log-logistic, Gamma, PGW, EW, GG
########################################################################################################
#' @param init    : initial point for optimisation step
#' under the parameterisation (log(scale), log(shape1), log(shape2), alpha, beta) for scale-shape1-shape2 models or
#' (mu, log(scale), alpha, beta) for log-location scale models.
#' @param hstr    : hazard structure:
#'            No covariates ("baseline"),
#'           AFT model with PGW baseline hazard ("AFT"),
#'           PH model with PGW baseline hazard ("PH"),
#'           AH model with PGW baseline hazard ("AH"),
#'           GH model with PGW baseline hazard ("GH")
#'           *GH is not available with Weibull dist
#' @param dist    : distribution for the baseline hazard:
#'                 Power Generalised Weibull ("PGW")
#'                 Generalised Gamma ("GenGamma"))
#'                 Exponentiated Weibull ("EW")
#'                 Weibull ("Weibull")
#'                 Gamma ("Gamma")
#'                 LogNormal ("LogNormal")
#'                 LogLogistic ("LogLogistic")
#' @param method  : "nlminb" or optimisation method to be used in optim (see ?optim)
#' @param maxit   : maximum number of iterations in optim or nlminb
#' @param times   : times to event
#' @param status  : vital status indicators (TRUE or 1 = observed, FALSE  or 0 = censored)
#' @param hp      : population hazard (for all individuals)
#' @param des     : design matrix for hazard-level effects
#' @param des_t   : design matrix for time-level effects (it is recommended not to use splines here)
#' @return It returns the output from optim or nlminb for the selected model and the negative log likelihood function
#' @export

GEHMLE <-
  function(init,
           times,
           status,
           hp,
           hstr = NULL,
           dist = NULL,
           des = NULL,
           des_t = NULL,
           method = "Nelder-Mead",
           maxit = 100) {

    # Required variables
    times <- as.vector(times)
    status <- as.vector(as.logical(status))
    times.obs <- times[status]
    if (!is.null(des)){
      des <- as.matrix(des)
      des.obs <- des[status, ]
    }
    if (!is.null(des_t)){
      des_t <- as.matrix(des_t)
      des_t.obs <- des_t[status, ]
    }
    hp <- as.vector(hp)
    hp.obs <- as.vector(hp[status])
    lhp.obs <- log(hp.obs)

    #------------------------------------------------------------------------------------
    # baseline models
    #------------------------------------------------------------------------------------

    if (hstr == "baseline") {
      # PGW
      if (dist == "PGW") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])

          MAT = cbind(lhp.obs , hpgw(times.obs, ae0, be0, ce0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chpgw(times, ae0, be0, ce0))
          return(val)
        }
      }

      # EW
      if (dist == "EW") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])

          MAT <- cbind(lhp.obs , hew(times.obs, ae0, be0, ce0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chew(times, ae0, be0, ce0))
          return(val)
        }
      }

      # GenGamma
      if (dist == "GenGamma") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])

          MAT <- cbind(lhp.obs , hggamma(times.obs, ae0, be0, ce0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chggamma(times, ae0, be0, ce0))
          return(val)
        }
      }

      # LogNormal
      if (dist == "LogNormal") {
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])

          MAT <- cbind(lhp.obs , hlnorm(times.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chlnorm(times, ae0, be0))
          return(val)
        }
      }

      # LogLogistic
      if (dist == "LogLogistic") {
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])

          MAT <- cbind(lhp.obs , hllogis(times.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chllogis(times, ae0, be0))
          return(val)
        }
      }

      # Gamma
      if (dist == "Gamma") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])

          MAT <- cbind(lhp.obs , hgamma(times.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chgamma(times, ae0, be0))
          return(val)
        }
      }

      # Weibull
      if (dist == "Weibull") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])

          MAT <- cbind(lhp.obs , hweibull(times.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chweibull(times, ae0, be0))
          return(val)
        }
      }

    }



    #------------------------------------------------------------------------------------
    # Proportional Hazards models
    #------------------------------------------------------------------------------------

    if (hstr == "PH") {

      # PGW
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

          MAT <- cbind(lhp.obs , hpgw(times.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <-  -sum(lhaz0) +
            sum(chpgw(times, ae0, be0, ce0) * exp.x.beta)
          return(val)
        }
      }

      # EW
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

          MAT <- cbind(lhp.obs , hew(times.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chew(times, ae0, be0, ce0) * exp.x.beta)
          return(val)
        }
      }

      # GenGamma
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

          MAT <- cbind(lhp.obs , hggamma(times.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chggamma(times, ae0, be0, ce0) * exp.x.beta)
          return(val)
        }
      }

      # LogNormal
      if (dist == "LogNormal") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))

          MAT <- cbind(lhp.obs , hlnorm(times.obs, ae0, be0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chlnorm(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }

      # LogLogistic
      if (dist == "LogLogistic") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))

          MAT <- cbind(lhp.obs , hllogis(times.obs, ae0, be0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chllogis(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }

      # Gamma
      if (dist == "Gamma") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))

          MAT <- cbind(lhp.obs , hgamma(times.obs, ae0, be0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chgamma(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }

      # Weibull
      if (dist == "Weibull") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))

          MAT <- cbind(lhp.obs , hweibull(times.obs, ae0, be0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chweibull(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }
    }


    #------------------------------------------------------------------------------------
    # Accelerated Failure Time models
    #------------------------------------------------------------------------------------

    if (hstr == "AFT") {

      # PGW
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

          MAT <- cbind(lhp.obs ,
                         hpgw(times.obs * exp.x.beta.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chpgw(times * exp.x.beta, ae0, be0, ce0))
          return(val)
        }
      }

      # EW
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

          MAT <- cbind(lhp.obs ,
                         hew(times.obs * exp.x.beta.obs, ae0, be0, ce0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chew(times * exp.x.beta, ae0, be0, ce0))
          return(val)
        }
      }

      # GenGamma
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

          MAT <- cbind(lhp.obs ,
                         hggamma(times.obs * exp.x.beta.obs, ae0, be0, ce0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <-  -sum(lhaz0) +
            sum(chggamma(times * exp.x.beta, ae0, be0, ce0))
          return(val)
        }
      }

      # LogNormal
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

          MAT <- cbind(lhp.obs ,
                         hlnorm(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <-  -sum(lhaz0) +
            sum(chlnorm(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }

      # LogLogistic
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

          MAT <- cbind(lhp.obs ,
                         hllogis(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chllogis(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }

      # Gamma
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

          MAT <- cbind(lhp.obs ,
                         hgamma(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chgamma(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }

      # Weibull
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

          MAT <- cbind(lhp.obs ,
                         hweibull(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chweibull(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }

    }

    #------------------------------------------------------------------------------------
    # Accelerated Hazards models
    #------------------------------------------------------------------------------------

    if (hstr == "AH") {

      # PGW
      if (dist == "PGW") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]

          MAT <- cbind(lhp.obs ,
                         hpgw(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chpgw(times * exp.x.alpha, ae0, be0, ce0) / exp.x.alpha)
          return(val)
        }
      }

      # EW
      if (dist == "EW") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]

          MAT <- cbind(lhp.obs ,
                         hew(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chew(times * exp.x.alpha, ae0, be0, ce0) / exp.x.alpha)
          return(val)
        }
      }

      # GenGamma
      if (dist == "GenGamma") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]

          MAT <- cbind(lhp.obs ,
                         hggamma(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chggamma(times * exp.x.alpha, ae0, be0, ce0) / exp.x.alpha)
          return(val)
        }
      }

      # LogNormal
      if (dist == "LogNormal") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]

          MAT <- cbind(lhp.obs ,
                         hlnorm(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chlnorm(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }

      # LogLogistic
      if (dist == "LogLogistic") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]

          MAT <- cbind(lhp.obs ,
                         hllogis(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chllogis(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }

      # Gamma
      if (dist == "Gamma") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]

          MAT <- cbind(lhp.obs ,
                         hgamma(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chgamma(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }

      # Weibull
      if (dist == "Weibull") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]

          MAT <- cbind(lhp.obs ,
                         hweibull(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chweibull(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }

    }

    #------------------------------------------------------------------------------------
    # General Hazards models
    #------------------------------------------------------------------------------------

    if (hstr == "GH") {

      # PGW
      if (dist == "PGW") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p0)]
          beta <- par[(4 + p0):(3 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          MAT <- cbind(lhp.obs ,
                         hpgw(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <-  -sum(lhaz0) +
            sum(chpgw(times * exp.x.alpha, ae0, be0, ce0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      # EW
      if (dist == "EW") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p0)]
          beta <- par[(4 + p0):(3 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          MAT <- cbind(lhp.obs ,
                         hew(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <-  -sum(lhaz0) +
            sum(chew(times * exp.x.alpha, ae0, be0, ce0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      # GenGamma
      if (dist == "GenGamma") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p0)]
          beta <- par[(4 + p0):(3 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          MAT <- cbind(lhp.obs ,
                         hggamma(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <-  -sum(lhaz0) +
            sum(chggamma(times * exp.x.alpha, ae0, be0, ce0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      # LogNormal
      if (dist == "LogNormal") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p0)]
          beta <- par[(3 + p0):(2 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          MAT <- cbind(lhp.obs ,
                         hlnorm(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chlnorm(times * exp.x.alpha, ae0, be0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      # LogLogistic
      if (dist == "LogLogistic") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p0)]
          beta <- par[(3 + p0):(2 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          MAT <- cbind(lhp.obs ,
                         hllogis(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chllogis(times * exp.x.alpha, ae0, be0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      # Gamma
      if (dist == "Gamma") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p0)]
          beta <- par[(3 + p0):(2 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          MAT <- cbind(lhp.obs ,
                         hgamma(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chgamma(times * exp.x.alpha, ae0, be0) * exp.x.beta.dif)
          return(sum(val))
        }
      }
    }

    #------------------------------------------------------------------------------------
    # Optimisation step
    #------------------------------------------------------------------------------------

    if (method != "nlminb") {
      OPT <- optim(init,
                   log.lik,
                   control = list(maxit = maxit),
                   method = method)
    }
    if (method == "nlminb") {
      OPT <- nlminb(init, log.lik, control = list(iter.max = maxit))
    }

    #------------------------------------------------------------------------------------
    # Output
    #------------------------------------------------------------------------------------
    OUT <- list(log_lik = log.lik, OPT = OPT)
    return(OUT)

  }

