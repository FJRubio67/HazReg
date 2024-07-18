## -------------------------------------------------------------------------------------------------------
rm(list=ls())

# Required packages
#library(devtools)
#install_github("FJRubio67/HazReg")
library(HazReg)
library(numDeriv)
library(knitr)
library(survival)
library(mvtnorm)

# Load the simulated data set
# https://github.com/FJRubio67/HazReg
df <- as.data.frame(read.table("dataGH.txt"))
colnames(df)
dim(df)

# Design matrix for hazard level effects
X <- as.matrix(cbind( df$agec, df$sex, df$TTT ))
colnames(X) <- cbind("std age", "sex", "trt")

# Design matrix for time level effects
Xt <- as.matrix(X)

q <- dim(Xt)[2]
p <- dim(X)[2]


# Vital status
status <- as.vector(df$status)

# Survival times
times <- as.vector(df$surv.time) # in years

# Population hazard rates for all individuals
hp <- df$rate

# Histogram of survival times
hist(times, probability = TRUE, breaks = 30)
box()


## -------------------------------------------------------------------------------------------------------
# PGWGH
OPTPGWGH <- GEHMLE(init = rep(0, 3 + ncol(X) + ncol(Xt)), times = times, status = status, hp = hp, 
                   hstr = "GH", dist = "PGW", des = X, des_t = Xt, method = "nlminb", maxit = 10000)

# PGWAFT
OPTPGWAFT <- GEHMLE(init = rep(0, 3 + ncol(X)), times = times, status = status, hp = hp, 
                    hstr = "AFT", dist = "PGW", des = X, method = "nlminb", maxit = 10000)

# PGWPH
OPTPGWPH <- GEHMLE(init = rep(0, 3 + ncol(X)), times = times, status = status, hp = hp, 
                   hstr = "PH", dist = "PGW", des = X, method = "nlminb", maxit = 10000)

# PGWAH
OPTPGWAH <- GEHMLE(init = rep(0, 3 + ncol(X)), times = times, status = status, hp = hp, 
                   hstr = "AH", dist = "PGW", des_t = X, method = "nlminb", maxit = 10000)


# LLGH
OPTLLGH <- GEHMLE(init = rep(0, 2 + ncol(X) + ncol(Xt)), times = times, status = status, hp = hp, 
                 hstr = "GH", dist = "LogLogistic", des = X, des_t = Xt, method = "nlminb", maxit = 10000)

# LLAFT
OPTLLAFT <- GEHMLE(init = rep(0, 2 + ncol(X)), times = times, status = status, hp = hp, 
                  hstr = "AFT", dist = "LogLogistic", des = X, method = "nlminb", maxit = 10000)

# LLPH
OPTLLPH <- GEHMLE(init = rep(0, 2 + ncol(X)), times = times, status = status, hp = hp, 
                 hstr = "PH", dist = "LogLogistic", des = X, method = "nlminb", maxit = 10000)

# LLAH
OPTLLAH <- GEHMLE(init = rep(0, 2 + ncol(X)), times = times, status = status, hp = hp, 
                 hstr = "AH", dist = "LogLogistic", des_t = X, method = "nlminb", maxit = 10000)


# EWGH
OPTEWGH <- GEHMLE(init = rep(0, 3 + ncol(X) + ncol(Xt)), times = times, status = status, hp = hp, 
                 hstr = "GH", dist = "EW", des = X, des_t = Xt, method = "nlminb", maxit = 10000)

# GGGH
OPTGGGH <- GEHMLE(init = rep(0, 3 + ncol(X) + ncol(Xt)), times = times, status = status, hp = hp, 
                 hstr = "GH", dist = "GenGamma", des = X, des_t = Xt, method = "nlminb", maxit = 10000)

# LNGH
OPTLNGH <- GEHMLE(init = rep(0, 2 + ncol(X) + ncol(Xt)), times = times, status = status, hp = hp, 
                 hstr = "GH", dist = "LogNormal", des = X, des_t = Xt, method = "nlminb", maxit = 10000)

# GGH
OPTGGH <- GEHMLE(init = rep(0, 2 + ncol(X) + ncol(Xt)), times = times, status = status, hp = hp, 
               hstr = "GH", dist = "Gamma", des = X, des_t = Xt, method = "nlminb", maxit = 10000)

# MLEs in the original parameterisations
MLEPGWGH <- c(exp(OPTPGWGH$OPT$par[1:3]),OPTPGWGH$OPT$par[-c(1:3)])
MLEEWGH <- c(exp(OPTEWGH$OPT$par[1:3]),OPTEWGH$OPT$par[-c(1:3)])
MLEGGGH <- c(exp(OPTGGGH$OPT$par[1:3]),OPTGGGH$OPT$par[-c(1:3)])
MLEGGH <- c(exp(OPTGGH$OPT$par[1:2]),OPTGGH$OPT$par[-c(1:2)])
MLELNGH <- c(OPTLNGH$OPT$par[1], exp(OPTLNGH$OPT$par[2]), OPTLNGH$OPT$par[-c(1,2)])
MLELLGH <- c(OPTLLGH$OPT$par[1], exp(OPTLLGH$OPT$par[2]), OPTLLGH$OPT$par[-c(1,2)])

MLES <- cbind(MLEPGWGH, MLEEWGH, MLEGGGH, c(MLEGGH[1:2], NA, MLEGGH[-(1:2)]), 
              c(MLELNGH[1:2], NA, MLELNGH[-(1:2)]), c(MLELLGH[1:2], NA, MLELLGH[-(1:2)]))
colnames(MLES) <- c("PGWGH", "EWGH", "GGGH", "GGH", "LNGH", "LLGH")
rownames(MLES) <- c("theta[1]","theta[2]","theta[3]","age_t", "sex_t", "trt_t","age", "sex", "trt")

# MLEs for GH models
kable(MLES, digits = 3)



## -------------------------------------------------------------------------------------------------------
# AIC for models with PGW baseline hazard
AICPGWGH <- 2*OPTPGWGH$OPT$objective + 2*length(OPTPGWGH$OPT$par)
AICPGWAFT <- 2*OPTPGWAFT$OPT$objective + 2*length(OPTPGWAFT$OPT$par)
AICPGWPH <- 2*OPTPGWPH$OPT$objective + 2*length(OPTPGWPH$OPT$par)
AICPGWAH <- 2*OPTPGWAH$OPT$objective + 2*length(OPTPGWAH$OPT$par)

# AICs for models with LL baseline hazard
AICLLGH <- 2*OPTLLGH$OPT$objective + 2*length(OPTLLGH$OPT$par)
AICLLAFT <- 2*OPTLLAFT$OPT$objective + 2*length(OPTLLAFT$OPT$par)
AICLLPH <- 2*OPTLLPH$OPT$objective + 2*length(OPTLLPH$OPT$par)
AICLLAH <- 2*OPTLLAH$OPT$objective + 2*length(OPTLLAH$OPT$par)

# AICs for GH models with GG, EW, LN, and G hazards
AICGGGH <- 2*OPTGGGH$OPT$objective + 2*length(OPTGGGH$OPT$par)
AICEWGH <- 2*OPTEWGH$OPT$objective + 2*length(OPTEWGH$OPT$par)
AICLNGH <- 2*OPTLNGH$OPT$objective + 2*length(OPTLNGH$OPT$par)
AICGGH <- 2*OPTGGH$OPT$objective + 2*length(OPTGGH$OPT$par)



# All AICs
AICs <- c(AICPGWGH, AICPGWAFT, AICPGWPH, AICPGWAH,
          AICLLGH, AICLLAFT, AICLLPH, AICLLAH,
          AICGGGH, AICEWGH, AICLNGH, AICGGH)

round(AICs, digits = 2)

# Best model: EWGH
which.min(AICs)


## -------------------------------------------------------------------------------------------------------
# Fitted baseline hazard functions for GH models
PGWGHhaz <- Vectorize(function(t) hpgw(t, MLEPGWGH[1], MLEPGWGH[2], MLEPGWGH[3]) ) 
EWGHhaz <- Vectorize(function(t) hew(t, MLEEWGH[1], MLEEWGH[2], MLEEWGH[3]) ) 
GGGHhaz <- Vectorize(function(t) hggamma(t, MLEGGGH[1], MLEGGGH[2], MLEGGGH[3]) ) 
GGHhaz <- Vectorize(function(t) hgamma(t, MLEGGH[1], MLEGGH[2]) ) 
LNGHhaz <- Vectorize(function(t) hlnorm(t, MLELNGH[1], MLELNGH[2])) 
LLGHhaz <- Vectorize(function(t) hllogis(t, MLELLGH[1], MLELLGH[2])) 

# Note that the baseline hazards associated to the top models look similar
curve(PGWGHhaz,1e-6, max(times), xlab = "Time (years)", ylab = "Baseline Hazard", main = "",
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, ylim = c(0,0.55), n = 1000)
curve(EWGHhaz,1e-6, max(times), lwd = 2, n = 1000, lty = 2, add = TRUE) 
curve(GGGHhaz,1e-6, max(times), lwd = 2, n = 1000, lty = 3, add = TRUE) 
curve(GGHhaz,1e-6, max(times), lwd = 2, n = 1000, lty = 4, add = TRUE) 
curve(LNGHhaz,1e-6, max(times), lwd = 2, n = 1000, lty = 5, add = TRUE) 
curve(LLGHhaz,1e-6, max(times), lwd = 2, n = 1000, lty = 6, add = TRUE) 
legend("topright", legend = c("PGWGH", "EWGH", "GGGH", "GGH", "LNGH", "LLGH"), lwd = rep(2,6), lty = 1:6)


## -------------------------------------------------------------------------------------------------------
# MLE in the original parameterisation
MLE <- c(exp(OPTEWGH$OPT$par[1:3]),OPTEWGH$OPT$par[-c(1:3)])

round(MLE, digits = 3)

# 95% Confidence intervals under the reparameterisation
CI <- Conf_Int(FUN = OPTEWGH$log_lik, MLE = OPTEWGH$OPT$par, level = 0.95)
rownames(CI) <- c("theta[1]","theta[2]","theta[3]","age_t", "sex_t", "trt_t","age", "sex", "trt")

kable(CI, digits = 3)

# Fitted baseline hazard function
fit_haz <- Vectorize(function(t) hew(t, MLE[1], MLE[2], MLE[3])) 

curve(fit_haz,0.001, max(times), xlab = "Time (years)", ylab = "Baseline Hazard", main = "",
      cex.axis = 1.5, cex.lab = 1.5, lwd = 2, ylim = c(0,0.3), n = 1000)

# Average population net survival function 

  p0 <- dim(Xt)[2]
  p1 <- dim(X)[2]
  theta1 <- MLE[1]; theta2 <- MLE[2]; theta3 <- MLE[3]; alpha <- MLE[4:(3+p0)]; beta <- MLE[(4+p0):(3+p0+p1)]
  x.alpha <- Xt%*%alpha
  x.dif <- X%*%beta - x.alpha

net_surv <- Vectorize(function(t){
  out <- mean( exp( - chew(t*exp(x.alpha), theta1, theta2, theta3)*exp(x.dif)  )  )
  return(out)
})


# Comparison
curve(net_surv,0.001,5, type = "l", col = "black", lwd = 2, lty = 1, ylim = c(0,1),
     xlab = "Time (years)", ylab = "Net Survival", main = "",
     cex.axis = 1.5, cex.lab = 1.5)


# Confidence intervals for the net survival function based on a normal approximation
# at specific time points t0

# Hessian and asymptotic covariance matrix
HESS <- hessian(func = OPTEWGH$log_lik, x = OPTEWGH$OPT$par)
Sigma <- solve(HESS)

# Reparameterised MLE 
r.MLE <- OPTEWGH$OPT$par

# The function to obtain approximate CIs based on Monte Carlo simulations 
# from the asymptotic normal distribution of the MLEs
# t0 : time where the confidence interval will be calculated
# level : confidence level
# n.mc : number of Monte Carlo iterations

conf.int.nsurv <- function(t0, level, n.mc){
  mc <- vector()
  S.par <- function(par){ mean( exp( - chew(t0*exp(Xt%*%par[4:(3+p0)]), par[1], par[2], par[3])*
                                      exp(X%*%par[(4+p0):(3+p0+p1)]-Xt%*%par[4:(3+p0)])  )  )
  }
  
  for(i in 1:n.mc) {
    val <- rmvnorm(1,mean = r.MLE, sigma = Sigma)
    val[1:3] <- exp(val[1:3])
    mc[i] <- S.par(val)
  }
  
  L <- quantile(mc,(1-level)*0.5)
  U <- quantile(mc,(1+level)*0.5)
  
  M <- S.par(MLE)
  
  return(c(L,M,U))
}


# times for CIs calculations
timesCI <- c(1,2,3,4,5)

CIS <- matrix(0, ncol = 4, nrow = length(timesCI))

for(k in 1:length(timesCI)) CIS[k,] <- c(timesCI[k],conf.int.nsurv(timesCI[k],0.95,10000))

colnames(CIS) <- cbind("year","lower","net survival","upper")
print(kable(CIS,digits=4))


