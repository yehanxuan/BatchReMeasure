# BatchReMeasure
Correct Batch effect by remeasure samples in completely confounded case-control studies. Currently, we only consider the control and case samples are collected in the first and second batches, respectively, and a subset of control samples are re-measured in the second batch.
We maximize the joint the model to obtain estimates of parameters like the true effects, batch effects. This package can evaluate the power of the developed procedure and provide a power calculation tool. 


# Installation 
The package can be installed using the following code.
```{r}
install.packages(BatchReMeasure)
```
or
```{r}
devtools::install_github('yehanxuan/BatchReMeasure')
library(BatchReMeasure)
```

# Example with simulated data 
The following script estimates the true and batch effects on a simulated data. We have both the R-version remeasure functions and the more computationally efficient C++ wrapper.
```{r}
# Loading packages 
library(RConics)
library(RcppArmadillo)
library(Rcpp)
library(BatchReMeasure)
```
Set sample size n=100 and remeasured sample size n1=40.
Samples consists of 50 control samples and 50 case samples.
Set the std of the first batch as 2 and the std of the second batch as 1; The correlation, true effect, batch effect are set as 0.6, 0.5 and 0.5, respectively 
```{r}
n = 100
n1 = 40 
r1 = 2
r2 = 0.6 
a0 = 0.5
a1 = 0.5
```
Simulate the data 
```{r}
X =  as.numeric(gl(2, n / 2)) - 1
Z <- cbind(rep(1, n), rnorm(n))
b <- c(0, -0.5)
v1 = r1^2
v2 = 1
Et <- rnorm(n, sd = ifelse (X == 0, sqrt(v1), sqrt(v2)))
Y <- Z %*% b + cbind(X, X) %*% c(a0, a1) + Et
Z.r.a <- Z[1 : (n / 2), ]
Et.r.a <- Et[1 : (n / 2)]
Y.r.a <- a1 + Z.r.a %*% b + r2 * sqrt(v2) * Et.r.a/ sqrt(v1) + rnorm(n/2, sd = sqrt( (1 - r2^2) * v2 ) )
ind.r <- 1:n1
Y.r = Y.r.a[ind.r]
```
Since the batch effects and the biological effects are indistinguishable in this example, we conduct batch effect correction with remeasurement
```{r}
Estimate = batch.ReMeasure.S1(Y, X, Z, ind.r, Y.r)
a0H = Estimate$a0
a0Var = Estimate$a0Var
a1H = Estimate$a1
```
The package also includes functions that only consider case and remeasured control samples in the second batch
```{r}
Estimate = batch.Batch2.S1(Y, X, Z, ind.r, Y.r)
```
If we ignore the batch effect, simply run
```{r}
Estimate = batch.Ignore.S1(Y, X, Z, ind.r, Y.r)
```
The bootstrap method is also provided to improve the estimation of uncertainty
```{r}
# The residual bootstrap 
boot = batch.ReMeasure.S1.res(Y, X, Z, ind.r, Y.r)
boot$ztest
boot$ztest_b
```



# References 
Eddelbuettel, Dirk, and Conrad Sanderson (2014). "RcppArmadillo: Accelerating R with high-performance C++ linear algebra." Computational Statistics & Data Analysis 71: 1054-1063.


