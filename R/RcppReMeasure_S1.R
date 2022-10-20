
RcppReMeasure_Estimate_S1 = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c = 1e-7, a0.Ini = NULL, a1.Ini=NULL,
                                     rho.Ini=NULL, beta.Ini = NULL, sigma1.Ini=NULL, sigma2.Ini=NULL) {
  nc2 = nrow(Zc2)
  if ( is.null(sigma1.Ini) ) {
    sigma1.Ini = sqrt( mean( (Yc1[Index] - mean(Yc1[Index]) )^2 ) )
  }

  if (sigma1.Ini <= 0.001 ){
    sigma1.Ini = 0.5
  }

  if ( is.null(sigma2.Ini) ) {
    sigma2.Ini = sqrt( mean( (Yc2 - mean(Yc2))^2 ) )
  }

  if (sigma2.Ini <= 0.001) {
    sigma2.Ini = 0.5
  }
  if ( is.null(rho.Ini) ) {
    # All use Control in case rho larger than 1
    rho.Ini = as.numeric( t(  (Yc1[Index] - mean(Yc1[Index])) ) %*%(Yc2 - mean(Yc2))/(nc2*sigma1.Ini*sigma2.Ini))
  }

  if ( rho.Ini > 0.99 ) {
    rho.Ini = 0.95
  } else if (rho.Ini <= -0.99) {
    rho.Ini = -0.95
  }
  if ( is.null(beta.Ini) ) {
    beta.Ini = solve( t(Zc1)%*%Zc1, t(Zc1)%*%Yc1)
  }
  if ( is.null(a1.Ini) ) {
    a1.Ini = mean(Yc2 - Zc2%*%beta.Ini) - rho.Ini*sigma2.Ini/sigma1.Ini*mean(Yc1[Index] - Zc2%*%beta.Ini)
  }
  if (is.null(a0.Ini)) {
    a0.Ini = mean(Yt2 - Zt2%*%beta.Ini) - a1.Ini
  }

  a0H = a0.Ini
  a1H = a1.Ini
  betaH = beta.Ini
  sigma1H = sigma1.Ini
  sigma2H = sigma2.Ini
  rhoH = rho.Ini

  obj_old = Rcpp_Objective_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, a1H, betaH, rhoH, sigma1H, sigma2H, Index)

  objVec = obj_old
  i = 0
  gap = 1e7

  start_S1 = proc.time()[1]
  while ( (i < 100)&&(gap > tol.c) ) {
    i = i + 1
    sigma1H = Update_sigma1_S1(Zc1, Zc2, Yc1, Yc2, a0H, a1H, betaH, rhoH, sigma2H, Index)
    sigma2H = Update_sigma2_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, a1H, betaH, rhoH, sigma1H, Index)
    rhoH = Update_rho_S1(Zc1, Zc2, Yc1, Yc2, a0H, a1H, betaH, sigma1H, sigma2H, Index)
    betaH = Update_beta_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, a1H, rhoH, sigma1H, sigma2H, Index)
    aVec = Update_a_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, a1H, betaH, rhoH, sigma1H, sigma2H, Index)
    a0H = aVec[2]
    a1H = aVec[1]
    obj_new = Rcpp_Objective_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, a1H, betaH, rhoH,
                                sigma1H, sigma2H, Index)
    gap = abs(obj_old - obj_new)
    objVec = c(objVec, obj_new)
    obj_old = obj_new
  }
  Time_S1 = proc.time()[1] - start_S1

  a0Var = Variance_a0(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, sigma1H, sigma2H, rhoH, Index)
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "beta" = betaH, "rho" = rhoH,
              "sigma1" = sigma1H, "sigma2" = sigma2H, "Time" = Time_S1, "objVec" = objVec))
}



#' The Rcpp Wrapper for batch.ReMeasure.S1 function for batch effect correction with remeasured samples.
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
#' @return The estimates of parameters through optimizing a joint model.
#' The parameters include true and batch effect, variances, correlation and time
#'
#' @examples
#' n = 100; n1 = 50; r1 = 1; r2 = 0.6; a0 = 0.5; a1 = 0.5
#' v1 = r1^2; v2 = 1
#' X =  as.numeric(gl(2, n / 2)) - 1
#' Z = cbind(rep(1, n), rnorm(n))
#' b = c(0, -0.5)
#' Et = rnorm(n, sd = ifelse (X == 0, sqrt(v1), sqrt(v2)))
#' Y = Z %*% b + cbind(X, X) %*% c(a0, a1) + Et
#' Z.r.a = Z[1 : (n / 2), ]
#' Et.r.a = Et[1 : (n / 2)]
#' Y.r.a = a1 + Z.r.a %*% b + r2 * sqrt(v2) * Et.r.a/ sqrt(v1) +
#' rnorm(n/2, sd = sqrt( (1 - r2^2) * v2 ) )
#' ind.r = 1:n1
#' Y.r = Y.r.a[ind.r]
#' # estimate the parameters
#' Estimate = Rcpp.batch.ReMeasure.S1(Y, X, Z, ind.r, Y.r)
#' Estimate$a0
#' Estimate$a1
#' Estimate$a0Var
#' Estimate$p.value
#'
#' @export
#'
Rcpp.batch.ReMeasure.S1 = function(Y, X, Z, ind.r, Y.r) {
  ind0 = X == 0
  ind1 = X == 1
  Yc1 = Y[ind0]
  Yt2 = Y[ind1]
  Zc1 = Z[ind0, , drop = F]
  Zt2 = Z[ind1, , drop = F]
  Zc2 = Zc1[ind.r, , drop = F]
  Yc2 = Y.r
  Estimate = RcppReMeasure_Estimate_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, ind.r, tol.c = 1e-7)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  betaH = Estimate$beta
  rhoH = Estimate$rho
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  pv = 2 * stats::pnorm(-abs(a0H / sqrt(a0Var)))
  Time = Estimate$Time
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "beta" = betaH, "rho" = rhoH, "p.value" = pv,
              "sigma1" = sigma1H, "sigma2" = sigma2H, "objVec" = objVec, "Time" = Time))
}









