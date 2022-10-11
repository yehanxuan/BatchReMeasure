
#' Power calculator using the theoretical formula
#'
#' @param nc1 size of control sample in the first batch
#' @param nt2 size of treatment sample in the second batch
#' @param r1 scale effect that determines the variance of the first batch
#' @param r2 correlation parameter
#' @param a0 true effect
#' @param a1 batch effect
#' @param nc2 remeasured size of control sample
#' @param alpha the nominal level
#' @param seedJ the seed
#'
#' @return The power calculated from theoretical formula
#'
#' @examples
#' alpha = 0.1
#' Calculate_Power_S1(nc1 = 50, nt2 = 50, r1=0.5, r2=0.6, a0 = 0.5, a1 = 0.5,
#' nc2 = 20, alpha = 0.1)
#' @export
Calculate_Power_S1 = function(nc1, nt2, r1, r2, a0, a1, nc2, alpha = 0.05,seedJ = 2 ){

  set.seed(seedJ)
  n = nc1 + nt2
  X =  as.numeric(gl(2, n / 2)) - 1
  Z <- cbind(rep(1, n), rnorm(n))
  b <- c(0, -0.5)

  v1 = r1^2
  v2 = 1
  Et <- rnorm(n, sd = ifelse (X == 0, sqrt(v1), sqrt(v2)))
  Y <- Z %*% b + cbind(X, X) %*% c(a0, a1) + Et

  Z.r.a <- Z[1 : (n / 2), ]
  Et.r.a <- Et[1 : (n / 2)]
  Y.r.a <- a1 + Z.r.a %*% b + r2 * sqrt(v2) * Et.r.a/ sqrt(v1) +
    rnorm(n/2, sd = sqrt( (1 - r2^2) * v2 ) )

  ind.r <- 1:nc2
  Y.r = Y.r.a[ind.r]

  ind0 <- X == 0
  ind1 <- X == 1
  Yc1 <- Y[ind0]
  Yt2 <- Y[ind1]
  Zc1 <- Z[ind0, , drop = F]
  Zt2 <- Z[ind1, , drop = F]
  Zc2 = Zc1[ind.r, , drop = F]
  Yc2 = Y.r

  a0Var = Oracle_Variance_a0(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, sqrt(v1), sqrt(v2), r2, ind.r)

  C = a0/sqrt(a0Var)
  C1 = qnorm(alpha/2, lower.tail = TRUE) - C
  C2 = qnorm(alpha/2, lower.tail = TRUE) + C
  Power = pnorm(C1) + pnorm(C2)
  #qnorm(1 - alpha/2) - C
  return(Power)
}



# Calculate_Power_S2 = function(nc1, nt2, r1, r2, r3, a0, a1, a3, nc3, nt3, seedJ = 2){
#   set.seed(seedJ)
#   n = nc1 + nt2
#   X =  as.numeric(gl(2, n / 2)) - 1
#   Z <- cbind(rep(1, n), rnorm(n))
#   b <- c(0, -0.5)
#   v1 = 2/(1+r1)
#   v2 = r1 * v1
#   v3 = r1 * v2
#   Et <- rnorm(n, sd = ifelse (X == 0, sqrt(v1), sqrt(v2)))
#   Y <- Z %*% b + cbind(X, X) %*% c(a0, a1) + Et
#
#   Z.r1.a <- Z[1 : (n / 2), ]
#   Et.r1.a <- Et[1 : (n / 2)]
#   Y.r1.a <- a3 + Z.r1.a %*% b + r2 * sqrt(v3) * Et.r1.a/sqrt(v1) +
#     rnorm(n/2, sd = sqrt( (1 - r2^2) * v3 ) )
#   Z.r2.a <- Z[(n/2+1) : n, , drop = F]
#   Et.r2.a <- Et[(n/2+1) : n]
#   Y.r2.a <- a0 + a3 + Z.r2.a %*% b + r3 * sqrt(v3) * Et.r2.a/sqrt(v2) +
#     rnorm(n/2, sd = sqrt( (1 - r3^2)*v3 ) )
#
#   ind.r1 <- 1:n1
#   Y.r1 = Y.r1.a[ind.r1]
#   ind.r2 <- 1:n2
#   Y.r2 = Y.r2.a[ind.r2]
#
#   ind0 <- X == 0
#   ind1 <- X == 1
#   Yc1 <- Y[ind0]
#   Yt2 <- Y[ind1]
#   Zc1 <- Z[ind0, , drop = F]
#   Zt2 <- Z[ind1, , drop = F]
#   Zc3 = Zc1[ind.r1, , drop = F]
#   Yc3 = Y.r1
#   Zt3 = Zt2[ind.r2, , drop = F]
#   Yt3 = Y.r2
#   a0Var = S2_Oracle_Variance_a0(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, sqrt(v1), sqrt(v2), sqrt(v3),
#                                 r2, r3, ind.r1, ind.r2)
#
#   C = a0/sqrt(a0Var)
#   C1 = qnorm(alpha/2, lower.tail = TRUE) - C
#   C2 = qnorm(alpha/2, lower.tail = TRUE) + C
#   Power = pnorm(C1) + pnorm(C2)
#   return(Power)
# }






