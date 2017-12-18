### 1a
set.seed(1)
x <- sort(runif(100)) + 100
y <- .2*(x-100 -.5)+(x-100 -.5)^2 + rnorm(100)*.1
#
plot(x,y)
b <- lm(y~x+I(x^2))
lines(x,fitted(b),lwd=2)
#
X <- model.matrix(b)
beta.hat <- solve(t(X)%*%X, t(X)%*%y)
#
x1 <- x+1000
plot(x,y)
b1 <- lm(y~x1+I(x1^2))
lines(x,fitted(b1),lwd=2)

### 1b
X1 <- model.matrix(b1)
kappa(X)
kappa(X1) # with numerical error matrix is singular => second order coeff not available

### 1c
plot(X[,2],X[,3])
plot(X1[,2],X1[,3])
cor(X[,2],X[,3]) # x & x^2 cols strongly correlated => high condition number
cor(X1[,2],X1[,3]) #with numerical error, x & x^2 cols are linearly dependent => condition number approaches infinity

### 1d
plot(x,y)
b2 <- lm(y~scale(x1)+I(scale(x1^2)))
lines(x,fitted(b2),lwd=2)
X2 <- model.matrix(b2)
kappa(X2)

### 1e
xs <- scale(x)
plot(x,y)
b3 <- lm(y~xs+I(xs^2))
lines(x, fitted(b3), lwd=2)
X3 <- model.matrix(b3)
kappa(X3) #much better
cor(X3[,2],X3[,3]) #much better

### 1f
b4 <- lm(y~poly(x,2))
X4 <- model.matrix(b4)
kappa(X4) #bit bigger than part (e)
cor(X4[,2],X4[,3]) #now super-tiny

### 1g
n <- nrow(x)
XX <- cbind(n^-.5,poly(x,2))
# QR decomp: Q=XX, R=I
# lm can be found using just t(Q)%*%y
t(XX)%*%y # gives same numbers as b4

### 2a
devtools::install_bitbucket("finnlindgren/FLtools")
library(FLtools)
An <- c(10, 20, 50, 100, 200, 500, 1000)
Atime <- c()
for (n in An) {
  loop.max <- max(1, 10000/n)
  A <- matrix(rnorm(n^2), n, n); A <- t(A) %*% A
  ## Use B <- A*1 to make sure R doesn't use any hidden precomputations
  Atime <- rbind(Atime,
                 system.time({
                   for (loop in 1:loop.max) {
                     B <- A*1
                     Matrix::chol(B)
                   }}) / loop.max)
}
Atime
plot(An, Atime[,1], log="", pch=16)
plot(An, Atime[,1], log="xy", pch=16)

### 2b
Atime2 <- c()
for (n in An) {
  loop.max <- max(1, 10000/n)
  A <- covariance.ar1(n, 0.6) #now an AR1 cov. matrix
  Atime2 <- rbind(Atime2,
                 system.time({
                   for (loop in 1:loop.max) {
                     B <- A*1
                     Matrix::chol(B)
                   }}) / loop.max)
}
Atime2
plot(An, Atime2[,1], log="", pch=16)
points(An, Atime[,1], pch=16, col=2)
plot(An, Atime2[,1], log="xy", pch=16)
points(An, Atime[,1], pch=16, col=2)
legend("topleft", c("unstructured", "AR1 covariance"), pch=16, col=c(1,2))

### 2c
Atime3 <- c()
for (n in An) {
  loop.max <- max(1, 10000/n)
  A <- precision.ar1(n, 0.6) #now an AR1 cov. matrix
  Atime3 <- rbind(Atime3,
                  system.time({
                    for (loop in 1:loop.max) {
                      B <- A*1
                      Matrix::chol(B)
                    }}) / loop.max)
}
Atime3
plot(An, Atime2[,1], log="xy", pch=16)
points(An, Atime[,1], pch=16, col=2)
points(An, Atime3[,1], pch=16, col=4)
legend("topleft", c("unstructured", "AR1 covariance", "AR1 precision"), pch=16, col=c(1,2,4))

# now over bigger An range
An2 <- c(10, 1e2, 1e3, 1e4, 1e5, 1e6)
Atime3a <- c()
for (n in An2) {
  loop.max <- max(1, 10000/n)
  A <- precision.ar1(n, 0.6) #now an AR1 cov. matrix
  Atime3a <- rbind(Atime3a,
                  system.time({
                    for (loop in 1:loop.max) {
                      B <- A*1
                      Matrix::chol(B)
                    }}) / loop.max)
}
Atime3a
plot(An2, Atime3a[,1], log="xy", pch=16)

### 2d
Atime4 <- c()
for (n in An) {
  loop.max <- max(1, 10000/n)
  A <- precision.ar1(n, 0.6) #now an AR1 cov. matrix
  Atime4 <- rbind(Atime4,
                   system.time({
                     for (loop in 1:loop.max) {
                       B <- A*1
                       base::chol(B)
                     }}) / loop.max)
}
Atime4
#
Atime4b <- c()
for (n in An) {
  loop.max <- max(1, 10000/n)
  A <- precision.ar1(n, 0.6) #now an AR1 cov. matrix
  Atime4b <- rbind(Atime4b,
                  system.time({
                    for (loop in 1:loop.max) {
                      B <- A*1
                      Matrix::chol(B)
                    }}) / loop.max)
}
Atime4b
#
plot(An, Atime4[,1], log="xy", pch=16)
points(An, Atime4b[,1], pch=16, col=2)
legend("topleft", c("dense", "sparse"), col=c(1,2), pch=16)

### 3a
FLtools::optimisation()

### 4a
set.seed(10)
n <- 100
n.b <- 10
n.beta <- 5
#
X <- cbind(1,matrix(runif(n*n.beta-n),n,n.beta-1))
Z <- matrix(runif(n*n.b),n,n.b)
beta <- rep(1,n.beta)
b <- rnorm(n.b)
y <- X%*%beta + Z%*%b + rnorm(n)

### 4b
# E(Y) = X%*%beta
# V(Y) = E(t(y)%*%y) = t(Z)%*%Z*sigma.b^2 + I*sigma^2
# Y ~ N(E(Y), V(Y))

### 4c
logLik <- function(theta,y,X,Z) {
  n <- length(y)
  beta <- theta[1:ncol(X)]
  theta <- theta[-(1:ncol(X))]
  V <- diag(n)*exp(theta[2]) + Z %*% t(Z)*exp(theta[1])
  R <- chol(V)
  z <- forwardsolve(t(R), y-X %*% beta)
  ll <- -n*log(2*pi)/2 - sum(log(diag(R))) - sum(z*z)/2
  ll
}

### 4d
optim(rep(-1,ncol(X)+2), logLik, y=y, X=X, Z=Z, control=list(fnscale=-1), method="BFGS")
#other methods agree with this maximum, except N-M which finds (presumably) another local max

### 4e

