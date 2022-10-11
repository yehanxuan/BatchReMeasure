
Estimate_ReMeasure_S1_Pair = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c = 1e-7, B, a0.Ini = NULL, a1.Ini=NULL,
                                      rho.Ini=NULL, beta.Ini = NULL, sigma1.Ini=NULL, sigma2.Ini=NULL) {

  start = proc.time()[1]
  nc1 = nrow(Zc1)
  nt2 = nrow(Zt2)
  nc2 = nrow(Zc2)
  out = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c)
  a0_hat = out$a0
  a0Var = out$a0Var
  a1_hat = out$a1
  beta_hat = out$beta
  sigma1_hat = out$sigma1
  sigma2_hat = out$sigma2
  rho_hat = out$rho
  ztest = a0_hat/sqrt(out$a0Var)
  objVec = out$objVec

  ztest_b = rep(NA, B)

  for (j in 1:B){
    ind_nc1_un = sample(1:(nc1-nc2), nc1-nc2, replace=TRUE)
    ind_nt2 = sample(1:nt2, nt2, replace=TRUE)
    ind_nc2 = sample(1:nc2, nc2, replace=TRUE)


    Zc1_b = Zc1
    Zc1_b[Index, ] = Zc1[Index, ,drop = F][ind_nc2, , drop = F]
    Zc1_b[-Index, ] = Zc1[-Index, , drop = F][ind_nc1_un, , drop = F]
    Yc1_b = as.matrix(Yc1)
    Yc1_b[Index] = Yc1[Index][ind_nc2]
    Yc1_b[-Index] = Yc1[-Index][ind_nc1_un]

    Yt2_b = Yt2[ind_nt2]
    Zt2_b = Zt2[ind_nt2, , drop = F]
    Yc2_b = Yc2[ind_nc2]
    Zc2_b = Zc2[ind_nc2, , drop = F]
    # The problem will occur when all the rows are the same value!
    # The initial estimation of sigma1,sigma2 = 0, not make sense!
    out_b = Estimate_ReMeasure_S1(Zc1_b, Zt2_b, Zc2_b, Yc1_b, Yt2_b, Yc2_b, Index, tol.c)
    ztest_b[j] = (out_b$a0-a0_hat)/sqrt(out_b$a0Var)
    if ( j%%200 == 0) {
      print(j)
    }
  }
  Time = proc.time()[1] - start
  return(list("a0" = a0_hat, "a0Var" = a0Var, "a1" = a1_hat, "beta" = beta_hat,
              "rho" = rho_hat, "sigma1" = sigma1_hat, "sigma2" = sigma2_hat, "Time" = Time,
              "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b))
}


Estimate_ReMeasure_S1_Res = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c = 1e-7, B, a0.Ini = NULL, a1.Ini=NULL,
                                      rho.Ini=NULL, beta.Ini = NULL, sigma1.Ini=NULL, sigma2.Ini=NULL) {
  nc1 = nrow(Zc1)
  nt2 = nrow(Zt2)
  nc2 = nrow(Zc2)
  # B Time of bootstrap
  start = proc.time()[1]
  out = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c)
  a0_hat = out$a0
  a0Var = out$a0Var
  a1_hat = out$a1
  beta_hat = out$beta
  sigma1_hat = out$sigma1
  sigma2_hat = out$sigma2
  rho_hat = out$rho
  ztest = a0_hat/sqrt(out$a0Var)
  objVec = out$objVec
  ## Estimation of noise
  ec1_hat = Yc1 - Zc1%*%beta_hat
  et2_hat = Yt2 - a0_hat - a1_hat - Zt2%*%beta_hat
  ec2_hat = Yc2 - a1_hat - Zc2%*%beta_hat
  ec1_hat_re = ec1_hat[Index, ,drop=F]
  ec1_hat_unre = ec1_hat[-Index, ,drop=F]

  ztest_b = rep(NA,B)
  ## Using bootstrap to estimate the standard error
  for (j in 1:B){
    ind_nc1_un = sample(1:(nc1-nc2), nc1-nc2, replace=TRUE)
    ind_nt2 = sample(1:nt2, nt2, replace=TRUE)
    ind_nc2 = sample(1:nc2, nc2, replace=TRUE)
    Yc1_b = as.matrix(Yc1)
    Yc1_b[-Index,] = Zc1[-Index,]%*%beta_hat + ec1_hat_unre[ind_nc1_un, ,drop=F]
    Yc1_b[Index,] = Zc1[Index,]%*%beta_hat + ec1_hat_re[ind_nc2, ,drop=F]
    Yt2_b = a0_hat + a1_hat + Zt2%*%beta_hat + et2_hat[ind_nt2, ,drop=F]
    Yc2_b = a1_hat + Zc2%*%beta_hat + ec2_hat[ind_nc2, ,drop=F]
    out_b = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1_b, Yt2_b, Yc2_b, Index, tol.c)
    ztest_b[j] = (out_b$a0-a0_hat)/sqrt(out_b$a0Var)
    if ( j%%200 == 0) {
      print(j)
    }
  }
  Time = proc.time()[1] - start
  # ztest_b approximate the centralized distribution of ztest
  # so we can compare the ztest with the quantile of ztest_b
  return(list("a0" = a0_hat, "a0Var" = a0Var, "a1" = a1_hat, "beta" = beta_hat,
              "rho" = rho_hat, "sigma1" = sigma1_hat, "sigma2" = sigma2_hat, "Time" = Time,
              "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b))
}


Estimate_ReMeasure_S1_Wild = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c = 1e-7, B, a0.Ini = NULL, a1.Ini=NULL,
                                         rho.Ini=NULL, beta.Ini = NULL, sigma1.Ini=NULL, sigma2.Ini=NULL) {
  nc1 = nrow(Zc1)
  nt2 = nrow(Zt2)
  nc2 = nrow(Zc2)
  start = proc.time()[1]
  out = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c)
  a0_hat = out$a0
  a0Var = out$a0Var
  a1_hat = out$a1
  beta_hat = out$beta
  sigma1_hat = out$sigma1
  sigma2_hat = out$sigma2
  rho_hat = out$rho
  ztest = a0_hat/sqrt(out$a0Var)
  objVec = out$objVec

  ztest_b = rep(NA,B)
  for (j in 1:B) {
    ec1 = rnorm(nc1, mean = 0, sd = sigma1_hat)
    Yc1_b = ec1 + Zc1%*%beta_hat
    et2 = rnorm(nt2, mean = 0, sd = sigma2_hat)
    Yt2_b = a0_hat + a1_hat + et2 + Zt2%*%beta_hat
    ecInd = rnorm(nc2, mean = 0, sd = sigma1_hat)
    ec2 = (rho_hat * ec1[Index] + sqrt(1 - rho_hat^2) * ecInd) * sigma2_hat/sigma1_hat
    Yc2_b = a1_hat + ec2 + Zc2%*%beta_hat
    out_b = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1_b, Yt2_b, Yc2_b, Index, tol.c)
    ztest_b[j] = (out_b$a0-a0_hat)/sqrt(out_b$a0Var)
    if ( j%%500 == 0) {
      print(j)
    }
  }

  Time = proc.time()[1] - start
  return(list("a0" = a0_hat, "a0Var" = a0Var, "a1" = a1_hat, "beta" = beta_hat,
              "rho" = rho_hat, "sigma1" = sigma1_hat, "sigma2" = sigma2_hat, "Time" = Time,
              "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b))
}


#' Pair bootstrap method to estimate the uncertainty. Sample pairs are drawn
#' with replacement from each group of different batches
#'
#' @param Y the response vector of control and case samples
#' @param X binary vector indicates control/case status
#' @param Z model matrix (sample x variable dimensions)
#' @param ind.r index of remeasured samples
#' @param Y.r the response vector of remeasured sample
#'
#' @return The estimated coefficients, z-test statistics provided by Remeasure method and pair bootstrap.
#' The statistics can be used to estimate the power
#'
#' @examples
#' n = 100; n1 = 40; r1 = 1; r2 = 0.6; a0 = 0.8; a1 = 0.5
#' v1 = r1^2; v2 = 1
#' X =  as.numeric(gl(2, n / 2)) - 1
#' Z <- cbind(rep(1, n), rnorm(n))
#' b <- c(0, -0.5)
#' Et <- rnorm(n, sd = ifelse (X == 0, sqrt(v1), sqrt(v2)))
#' Y <- Z %*% b + cbind(X, X) %*% c(a0, a1) + Et
#' Z.r.a <- Z[1 : (n / 2), ]
#' Et.r.a <- Et[1 : (n / 2)]
#' Y.r.a <- a1 + Z.r.a %*% b +
#'  r2 * sqrt(v2) * Et.r.a/ sqrt(v1) + rnorm(n/2, sd = sqrt( (1 - r2^2) * v2 ) )
#' ind.r <- 1:n1
#' Y.r = Y.r.a[ind.r]
#'
#' boot = batch.ReMeasure.S1.pair(Y, X, Z, ind.r, Y.r)
#' boot$a0
#' boot$a0Var
#' boot$a1
#' boot$ztest
#' boot$ztest_b
#' @export
#'
batch.ReMeasure.S1.pair = function(Y, X, Z, ind.r, Y.r) {
  ind0 <- X == 0
  ind1 <- X == 1
  Yc1 <- Y[ind0]
  Yt2 <- Y[ind1]
  Zc1 <- Z[ind0, , drop = F]
  Zt2 <- Z[ind1, , drop = F]

  Zc2 = Zc1[ind.r, , drop = F]
  Yc2 = Y.r
  Estimate = Estimate_ReMeasure_S1_Pair(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, ind.r, B = 1000)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  betaH = Estimate$beta
  rhoH = Estimate$rho
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  Time = Estimate$Time
  ztest = Estimate$ztest
  ztest_b = Estimate$ztest_b
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "beta" = betaH, "rho" = rhoH,
              "sigma1" = sigma1H, "sigma2" = sigma2H, "Time" = Time, "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b))
}

#' Residual bootstrap method to estimate the uncertainty.
#' The residuals are sampled with replacement and new bootstrap samples are generated
#'
#' @param Y the response vector of control and case samples
#' @param X binary vector indicates control/case status
#' @param Z model matrix (sample x variable dimensions)
#' @param ind.r index of remeasured samples
#' @param Y.r the response vector of remeasured sample
#'
#' @return The estimated coefficients, z-test statistics provided by Remeasure method and residual bootstrap.
#' The statistics can be used to estimate the power
#'
#' @examples
#' n = 100; n1 = 40; r1 = 1; r2 = 0.6; a0 = 0.8; a1 = 0.5
#' v1 = r1^2; v2 = 1
#' X =  as.numeric(gl(2, n / 2)) - 1
#' Z <- cbind(rep(1, n), rnorm(n))
#' b <- c(0, -0.5)
#' Et <- rnorm(n, sd = ifelse (X == 0, sqrt(v1), sqrt(v2)))
#' Y <- Z %*% b + cbind(X, X) %*% c(a0, a1) + Et
#' Z.r.a <- Z[1 : (n / 2), ]
#' Et.r.a <- Et[1 : (n / 2)]
#' Y.r.a <- a1 + Z.r.a %*% b +
#' r2 * sqrt(v2) * Et.r.a/ sqrt(v1) + rnorm(n/2, sd = sqrt( (1 - r2^2) * v2 ) )
#' ind.r <- 1:n1
#' Y.r = Y.r.a[ind.r]
#'
#' boot = batch.ReMeasure.S1.res(Y, X, Z, ind.r, Y.r)
#' boot$a0
#' boot$a0Var
#' boot$a1
#' boot$ztest
#' boot$ztest_b
#' @export
#'
batch.ReMeasure.S1.res = function(Y, X, Z, ind.r, Y.r) {
  ind0 <- X == 0
  ind1 <- X == 1
  Yc1 <- Y[ind0]
  Yt2 <- Y[ind1]
  Zc1 <- Z[ind0, , drop = F]
  Zt2 <- Z[ind1, , drop = F]

  Zc2 = Zc1[ind.r, , drop = F]
  Yc2 = Y.r
  Estimate = Estimate_ReMeasure_S1_Res(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, ind.r, B = 1000)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  betaH = Estimate$beta
  rhoH = Estimate$rho
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  Time = Estimate$Time
  ztest = Estimate$ztest
  ztest_b = Estimate$ztest_b

  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "beta" = betaH, "rho" = rhoH,
              "sigma1" = sigma1H, "sigma2" = sigma2H, "Time" = Time, "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b))
}


#' Wild bootstrap method to estimate the uncertainty.
#' We fix the covariates and generate independent Gaussian noise with the estimated variance
#'
#' @param Y the response vector of control and case samples
#' @param X binary vector indicates control/case status
#' @param Z model matrix (sample x variable dimensions)
#' @param ind.r index of remeasured samples
#' @param Y.r the response vector of remeasured sample
#'
#' @return The estimated coefficients, z-test statistics provided by Remeasure method and wild bootstrap.
#' The statistics can be used to estimate the power
#'
#' @examples
#' n = 100; n1 = 40; r1 = 1; r2 = 0.6; a0 = 0.8; a1 = 0.5
#' v1 = r1^2; v2 = 1
#' X =  as.numeric(gl(2, n / 2)) - 1
#' Z <- cbind(rep(1, n), rnorm(n))
#' b <- c(0, -0.5)
#' Et <- rnorm(n, sd = ifelse (X == 0, sqrt(v1), sqrt(v2)))
#' Y <- Z %*% b + cbind(X, X) %*% c(a0, a1) + Et
#' Z.r.a <- Z[1 : (n / 2), ]
#' Et.r.a <- Et[1 : (n / 2)]
#' Y.r.a <- a1 + Z.r.a %*% b +
#' r2 * sqrt(v2) * Et.r.a/ sqrt(v1) + rnorm(n/2, sd = sqrt( (1 - r2^2) * v2 ) )
#' ind.r <- 1:n1
#' Y.r = Y.r.a[ind.r]
#'
#' boot = batch.ReMeasure.S1.wild(Y, X, Z, ind.r, Y.r)
#' boot$a0
#' boot$a0Var
#' boot$a1
#' boot$ztest
#' boot$ztest_b
#' @export
#'
batch.ReMeasure.S1.wild = function(Y, X, Z, ind.r, Y.r) {
  ind0 <- X == 0
  ind1 <- X == 1
  Yc1 <- Y[ind0]
  Yt2 <- Y[ind1]
  Zc1 <- Z[ind0, , drop = F]
  Zt2 <- Z[ind1, , drop = F]
  Zc2 = Zc1[ind.r, , drop = F]
  Yc2 = Y.r
  Estimate = Estimate_ReMeasure_S1_Wild(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, ind.r, B = 1000)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  betaH = Estimate$beta
  rhoH = Estimate$rho
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  Time = Estimate$Time
  ztest = Estimate$ztest
  ztestb = Estimate$ztest_b
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "beta" = betaH, "rho" = rhoH,
              "sigma1" = sigma1H, "sigma2" = sigma2H, "Time" = Time, "objVec" = objVec, "ztest" = ztest, "ztestb" = ztestb))
}











