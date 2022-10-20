
Variance_a0_OnlyBatch2 = function(Zt2, Zc2, Yt2, Yc2, sigma2H) {
  nt2 = nrow(Zt2)
  nc2 = nrow(Zc2)
  Zc2_c =  t( t(Zc2) - colMeans(Zc2) )
  Yc2_c = Yc2 - mean(Yc2)
  Zt2_c = t( t(Zt2) - colMeans(Zt2) )
  Yt2_c = Yt2 - mean(Yt2)

  Cov =  (t(Zt2_c)%*%Zt2_c + t(Zc2_c)%*%Zc2_c) +  0.001*diag(ncol(Zc2))

  diff = colMeans(Zt2) - colMeans(Zc2)
  c1 = t(rep(1, nt2))/nt2 - t(diff)%*%solve(Cov, t(Zt2_c))
  c2 = -t(rep(1, nc2))/nc2 - t(diff)%*%solve(Cov, t(Zc2_c))
  a0Var = (sigma2H^2)*( c1%*%t(c1) + c2%*%t(c2) )
  return(a0Var)
}


Estimate_OnlyBatch2 = function(Zt2, Zc2, Yt2, Yc2) {
  start = proc.time()[1]
  nt2 = nrow(Zt2)
  nc2 = nrow(Zc2)
  Zc2_c =  t( t(Zc2) - colMeans(Zc2) )
  Yc2_c = Yc2 - mean(Yc2)
  Zt2_c = t( t(Zt2) - colMeans(Zt2) )
  Yt2_c = Yt2 - mean(Yt2)

  Cov =  (t(Zt2_c)%*%Zt2_c + t(Zc2_c)%*%Zc2_c) + 0.001*diag(ncol(Zc2))
  Cor = (t(Zt2_c)%*%Yt2_c + t(Zc2_c)%*%Yc2_c)
  betaH = solve(Cov, Cor)
  a1H = mean(Yc2 - Zc2%*%betaH)
  a0H = mean(Yt2 - Zt2%*%betaH) - mean(Yc2 - Zc2%*%betaH)

  tmp  = sum( (Yt2_c - Zt2_c%*%betaH)^2) + sum( (Yc2_c - Zc2_c%*%betaH)^2)
  sigma2HSq = tmp/(nc2 + nt2)
  sigma2H = sqrt(sigma2HSq)


  a0Var = Variance_a0_OnlyBatch2(Zt2, Zc2, Yt2, Yc2, sigma2H)
  Time = proc.time()[1] - start

  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "beta" = betaH, "sigma2" = sigma2H,
              "Time" = Time))
}



#' Estimate the true effect and adjust for the batch effect only using the samples in the second batch
#'
#' @param Y the response vector for both control and case samples
#' @param X binary vector indicates control/case status.
#' X = 0 represents control, X = 1 represents case
#' @param Z model matrix (sample x variable dimensions) consists of covariates that affect the response
#' @param ind.r index for samples remeasured from control samples in the first batch.
#' The length of the index should be less or equal to the sample size of control samples in the first batch
#' @param Y.r the response vector of remeasured samples.
#' Due to batch effects, it is usually not equal to the responses of corresponding control samples in the first batch
#'
#'
#' @return The estimates of parameters including the true effect,
#' batch effect, variances, correlation and computational time
#'
#' @examples
#'
#' n = 100; n1 = 40 ;r1 = 2; r2 = 0.6;a0 = 0.5; a1 = 0.5
#' v1 = r1^2; v2 = 1
#' X =  as.numeric(gl(2, n / 2)) - 1
#' Z =cbind(rep(1, n), rnorm(n))
#' b =c(0, -0.5)
#' Et =rnorm(n, sd = ifelse (X == 0, sqrt(v1), sqrt(v2)))
#' Y =Z %*% b + cbind(X, X) %*% c(a0, a1) + Et
#' Z.r.a =Z[1 : (n / 2), ]
#' Et.r.a =Et[1 : (n / 2)]
#' Y.r.a =a1 + Z.r.a %*% b +
#' r2 * sqrt(v2) * Et.r.a/ sqrt(v1) + rnorm(n/2, sd = sqrt( (1 - r2^2) * v2 ) )
#' ind.r = 1:n1
#' Y.r = Y.r.a[ind.r]
#' # estimate the parameters
#' Estimate = batch.Batch2.S1(Y, X, Z, ind.r, Y.r)
#' Estimate$a0
#' Estimate$a0Var
#' Estimate$a1
#'
#' @export
#'
batch.Batch2.S1 = function(Y, X, Z, ind.r, Y.r) {

  ind0 = X == 0
  ind1 = X == 1
  Yc1 = Y[ind0]
  Yt2 = Y[ind1]
  Zc1 = Z[ind0, , drop = F]
  Zt2 = Z[ind1, , drop = F]

  Zc2 = Zc1[ind.r, , drop = F]
  Yc2 = Y.r
  Estimate = Estimate_OnlyBatch2(Zt2, Zc2, Yt2, Yc2)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  betaH = Estimate$beta
  rhoH = NULL
  sigma1H = NULL
  sigma2H = Estimate$sigma2
  objVec = NULL
  pv = 2 * stats::pnorm(-abs(a0H / sqrt(a0Var)))
  Time = Estimate$Time
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "beta" = betaH, "rho" = rhoH, "p.value" = pv,
              "sigma1" = sigma1H, "sigma2" = sigma2H, "objVec" = objVec, "Time" = Time))
}











