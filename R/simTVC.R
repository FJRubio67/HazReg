rm(list=ls())

# Required packages
#library(devtools)
#install_github("FJRubio67/SimLT")
library(SimLT)


#------------------------------------------------------------------------
# Baseline design matrix simulation
#------------------------------------------------------------------------
# Sample size
n <- 1000
# Design matrix
Xcf <- simDesMatrix(seed = 123, n = n, admin.cens = "31-12-2015", scale.age = FALSE,
                    site = "colon", sex = "female")
head(Xcf)

