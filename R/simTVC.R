rm(list=ls())

# Required packages
#library(devtools)
#install_github("FJRubio67/SimLT")
library(SimLT)


#------------------------------------------------------------------------
# Baseline design matrix simulation
#------------------------------------------------------------------------
# Seed
seed <- 123
set.seed(seed)

# Sample size
n <- 1000
# Design matrix
Xcf <- simDesMatrix(seed = seed, n = n, admin.cens = "31-12-2015", scale.age = FALSE,
                    site = "colon", sex = "female")
head(Xcf)

# Number of individual contact points

lambda0 <- 2
ncp <- rpois(n = n, lambda = lambda0)
ncps <- sum(ncp)

# IDs
id <- rep(seq_along(ncp), times = ncp)


# Contact points
# ncp: vector of counts, length n
# fup: vector of follow-up times, length n


