set.seed(1234)
des1 <- cbind(rnorm(10),rnorm(10))
des2 <- cbind(rnorm(10),rnorm(10))
des3 <- cbind(rnorm(10),rnorm(10))
des4 <- cbind(rnorm(10),rnorm(10))
des5 <- cbind(rnorm(10),rnorm(10))
des6 <- cbind(rnorm(10),rnorm(10))


des_list <- list(des1, des2, des3, des4, des5, des6)
time <- c(0, 1, 2, 3, 4, 5)

n <- nrow(des1)   # number of individuals
p <- ncol(des1)   # number of covariates


df <- do.call(
  rbind,
  lapply(seq_along(time), function(k) {
    data.frame(
      ID   = 1:n,
      time = time[k],
      des1 = des_list[[k]][, 1],
      des2 = des_list[[k]][, 2]
    )
  })
)


beta <- c(0.1, 0.2)
ae0  <- 0
be0  <- 1

npar = 2

SPred_TVC(df = df,
          beta = beta,
          npar = npar,
          ae0 = ae0,
          be0 = be0,
          chfun = chlnorm,
          hstr = "PH")
