GHMLE_TVC <- function (init, times, status, hstr = NULL, dist = NULL, des = NULL, 
          des_t = NULL, method = "Nelder-Mead", maxit = 100) 
{
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  times.obs <- times[status]
  if (!is.null(des)) {
    des <- as.matrix(des)
    des.obs <- des[status, ]
  }
  if (!is.null(des_t)) {
    des_t <- as.matrix(des_t)
    des_t.obs <- des_t[status, ]
  }
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
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hpgw(times.obs, ae0, be0, ce0, log = TRUE) + 
          x.beta.obs
        val <- -sum(lhaz0) + sum(chpgw(times, ae0, be0, 
                                       ce0) * exp.x.beta)
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
        val <- -sum(lhaz0) + sum(chew(times, ae0, be0, 
                                      ce0) * exp.x.beta)
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
        val <- -sum(lhaz0) + sum(chggamma(times, ae0, 
                                          be0, ce0) * exp.x.beta)
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
        val <- -sum(lhaz0) + sum(chlnorm(times, ae0, 
                                         be0) * exp.x.beta)
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
        val <- -sum(lhaz0) + sum(chllogis(times, ae0, 
                                          be0) * exp.x.beta)
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
        lhaz0 <- hgamma(times.obs, ae0, be0, log = TRUE) + 
          x.beta.obs
        val <- -sum(lhaz0) + sum(chgamma(times, ae0, 
                                         be0) * exp.x.beta)
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
        val <- -sum(lhaz0) + sum(chweibull(times, ae0, 
                                           be0) * exp.x.beta)
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
        val <- -sum(lhaz0) + sum(chpgw(times * exp.x.beta, 
                                       ae0, be0, ce0))
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
        val <- -sum(lhaz0) + sum(chew(times * exp.x.beta, 
                                      ae0, be0, ce0))
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
        val <- -sum(lhaz0) + sum(chggamma(times * exp.x.beta, 
                                          ae0, be0, ce0))
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
        val <- -sum(lhaz0) + sum(chlnorm(times * exp.x.beta, 
                                         ae0, be0))
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
        val <- -sum(lhaz0) + sum(chllogis(times * exp.x.beta, 
                                          ae0, be0))
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
        val <- -sum(lhaz0) + sum(chgamma(times * exp.x.beta, 
                                         ae0, be0))
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
        val <- -sum(lhaz0) + sum(chweibull(times * exp.x.beta, 
                                           ae0, be0))
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