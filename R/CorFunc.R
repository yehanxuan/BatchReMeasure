# Need to consider which the level of X corresponds control/case
#require(lmvar)
#require(MASS)
loglik.r1 <- function (paras, Y, X, Z, ind.r, Y.r) {
	n <- length(Y)
	n1 <- length(ind.r)
	p <- ncol(Z)
	X <- as.numeric(factor(X)) - 1
	XX <- cbind(X, X)
	a0 <- paras[1]
	a1 <- paras[2]
	b <- paras[3 : (p + 2)]
	v1 <- exp(paras[p + 3])
	v2 <- exp(paras[p + 4])
	rho <- (1 - exp(paras[p + 5])) / (1 + exp(paras[p + 5]))

	mu <- Z %*% b + XX %*% c(a0, a1)
	s2 <- c(v1, v2)[X + 1]

	Z.r <- Z[ind.r, ]
	Y.o <- Y[ind.r]
	mu.r <- a1 + Z.r %*% b + rho * sqrt(v2 / v1) * (Y.o - Z.r %*% b)
	s2.r <- rep((1 - rho^2) * v2, n1)

	-sum(stats::dnorm(c(Y, Y.r), mean = c(mu, mu.r), sd = sqrt(c(s2, s2.r)), log = TRUE))

}



#' Generic maximum likelihood estimation. The objective is the same with remeasured method
#' except that the generic optimization methods are used like BFGS or conjugate gradient (CG)
#'
#' @param Y the response vector of control and case samples
#' @param X binary vector indicate the control or case
#' @param Z model matrix (sample x variable dimensions)
#' @param ind.r index of remeasured samples
#' @param Y.r the response vector of remeasured samples
#' @param start the initial value of parameters, default is NULL
#' @param method the optimization method chosen for the optimization, default is BFGS
#'
#' @return The estimated coefficient of MLE
#'
#' @examples
#'
#' # generate the data
#' n = 100; n1 = 10 ;r1 = 1; r2 = 0.6;a0 = 0.8; a1 = 0.5
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
#' # estimate the parameters
#' Estimate = batch.correct.r1(Y, X, Z, ind.r, Y.r)
#' Estimate$a0
#' Estimate$a0Var
#' Estimate$a1
#'
#' @export
#'
batch.correct.r1 <- function (Y, X, Z, ind.r, Y.r, start = NULL, method = 'BFGS') {
	# Z include the intercept
  # X treatment and case
  start_S1 = proc.time()[1]
	call <- match.call()
	if (is.null(start)) {

		Y.o <- Y[ind.r]
		c2 <-  mean(Y.r - Y.o) ## a1
		lm.obj <- lm(Y ~ X + Z - 1)
		lm.coef <- coef(lm.obj)
		c1 <- lm.coef[1] - c2
		c3 <- lm.coef[2 : length(lm.coef)]
		c4 <- c5 <-  log(sigma(lm.obj)^2)   # sigma residual standard error
		c6 <- cor(Y.o - resid(lm.obj)[ind.r], Y.r - resid(lm.obj)[ind.r])
		c6 <- log((1 - c6) / (1 + c6))

		start <- c(c1, c2, c3, c4, c5, c6)  # a0, a1, beta, sigma, rho, here is slightly
		#c5 not c6
	}

	oout <- stats::optim(start, loglik.r1, method = method, hessian = TRUE,
			Y = Y, X = X, Z = Z, ind.r = ind.r, Y.r = Y.r)
	coef <- oout$par
	try.obj <- try({	vcov <- solve(oout$hessian)
				stat <- coef[1] / sqrt(vcov[1, 1])
				pv <- 2 * stats::pnorm(-abs(stat))
	#			pv <- 2 * stats::pt(-abs(coef[1] / sqrt(vcov[1, 1])), df = length(ind.r) + sum(X == 1) - ncol(Z))
			})
    if (inherits(try.obj, 'try-error')) {
		stat <- NA
		pv <- NA
    }

	a0H = as.numeric(coef[1])
	a1H = as.numeric(coef[2])
	betaH = as.vector( coef[3:4])
	sigma1H = exp( coef[5]/2 )
	sigma2H = exp( coef[6]/2 )

	rhoH = (1-exp(coef[7]))/(1 + exp(coef[7]))
	a0Var = as.numeric(vcov[1, 1])
	Time_S1 = proc.time()[1] - start_S1
	#return(list(call = call, stat = stat, p.value = pv, optim.out = oout, optim.st = start ))
	return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "sigma1" = sigma1H,
	            "sigma2" = sigma2H, "rho" = rhoH, "beta" = NULL, "objVec" = NULL, "Time" = Time_S1))
}


loglik.r2 <- function (paras, Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2) {
	n <- length(Y)
	n1 <- length(ind.r1)
	n2 <- length(ind.r2)
	ind0 <- X==0
	ind1 <- X==1
	Yc1 <- Y[ind0]
	Yt2 <- Y[ind1]
	Zc1 <- Z[ind0, , drop = F]
	Zt2 <- Z[ind1, , drop = F]

	p <- ncol(Z)
	X <- as.numeric(factor(X)) - 1
	XX <- cbind(X, X)
	a0 <- paras[1]
	a1 <- paras[2]
	a2 <- paras[3]
	b <- paras[4 : (p + 3)]
	v1 <- exp(paras[p + 4])
	v2 <- exp(paras[p + 5])
	v3 <- exp(paras[p + 6])


	rho1 <- (1 - exp(paras[p + 7])) / (1 + exp(paras[p + 7]))
	## Modify this
	rho2 = (1 - exp(paras[p + 8])) / (1 + exp(paras[p + 8]))
	# rho2 <- sqrt(v1 / v2)

	mu <- Z %*% b + XX %*% c(a0, a1)
	s2 <- c(v1, v2)[X + 1]

	Z.r1 <- Zc1[ind.r1, , drop =F]
	Y.o1 <- Yc1[ind.r1]
 	mu.r1 <- a2 + Z.r1 %*% b + rho1 * sqrt(v3 / v1) * (Y.o1 - Z.r1 %*% b)
	s2.r1 <- rep((1 - rho1^2) * v3, n1)

	Z.r2 <- Zt2[ind.r2, , drop = F]
	Y.o2 = Yt2[ind.r2]
	mu.r2 <- a0 + a2 + Z.r2 %*% b + rho2 * sqrt(v3 / v2) * (Y.o2 - a0 - a1 -  Z.r2 %*% b)
	s2.r2 <- rep((1 - rho2^2) * v3, n2)

	-sum(stats::dnorm(c(Y, Y.r1, Y.r2), mean = c(mu, mu.r1, mu.r2), sd = sqrt(c(s2, s2.r1, s2.r2)), log = TRUE))

}

batch.correct.r2 <- function (Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2, start = NULL, method = 'BFGS') {
	# Z include the intercept
  start_S2 = proc.time()[1]
	call <- match.call()
	if (is.null(start)) {

		Y.o1 <- Y[ind.r1]
		c3 <-  mean(Y.r1 - Y.o1) # a3
		Y.o2 <- Y[ind.r2]
		c23 <-  mean(Y.r2 - Y.o2)
		c2 <- c3 - c23  # a1

		lm.obj <- lm(Y ~ X + Z - 1)

		lm.coef <- coef(lm.obj)
		c1 <- lm.coef[1] - c2  # a0
		c4 <- lm.coef[2 : length(lm.coef)]
		c5 <- c6 <- c7 <- log(sigma(lm.obj)^2)
		c8 <- cor(Y.o1 - resid(lm.obj)[ind.r1], Y.r1 - resid(lm.obj)[ind.r1])
		c8 <- log((1 - c8) / (1 + c8))
		c9 <- cor(Y.o2 - resid(lm.obj)[ind.r2], Y.r2 - resid(lm.obj)[ind.r2])
		c9 <- log((1 - c9) / (1 + c9))

		start <- c(c1, c2, c3, c4, c5, c6, c7, c8, c9)
	}

	p = ncol(Z)
	oout <- stats::optim(start, loglik.r2, method = method, hessian = TRUE,
			Y = Y, X = X, Z = Z, ind.r1 = ind.r1, ind.r2 = ind.r2,  Y.r1 = Y.r1, Y.r2 = Y.r2)
	coef <- oout$par
	try.obj <- try({	vcov <- solve(oout$hessian)
				pv <- 2 * stats::pnorm(-abs(coef[1] / sqrt(vcov[1, 1])))
			})
	if (inherits(try.obj, 'try-error')) {
		pv <- NA
	}
	a0H = as.numeric(coef[1])
	a1H = as.numeric(coef[2])
	a3H = as.numeric(coef[3])
	betaH = as.vector(coef[4 : (p+3) ])
	sigma1H = exp(coef[p+4]/2)
	sigma2H = exp(coef[p+5]/2)
	sigma3H = exp(coef[p+6]/2)
	rho1H = (1 - exp(coef[p+7]) )/(1 + exp(coef[p+7]))
	rho2H = (1 - exp(coef[p+8]) )/(1 + exp(coef[p+8]))
	a0Var = as.numeric(vcov[1, 1])
	Time_S2 = proc.time()[1] - start_S2
	# return(list(call = call, p.value = pv, optim.out = oout, optim.st = start ))
	return(list("a0"=a0H, "a0Var"=a0Var, "a1"=a1H, "a3" = a3H, "sigma1"=sigma1H,
	            "sigma2"=sigma2H, "sigma3"=sigma3H, "rho1"=rho1H, "rho2"=rho2H,
	            "beta" = NULL, "objVec"=NULL, "Time"=Time_S2 ))
}




