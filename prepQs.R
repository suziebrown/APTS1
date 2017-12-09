#1
eps<-1
x<-0.1
while (x+eps != x){eps <- eps/2}
eps

#2a
x <- matrix(runif(1e5),1e3,1e2)
sum1 <- function(x) {
    z <- rep(0,1e3)
    for (i in 1:1e3) {
        for (j in 1:1e2) {
            z[i] <- z[i] + x[i,j]
        }
    }
    z
}
sum1a <- function(x) {
    z <- rep(0,1e3)
    for (i in 1:1e3) {
        z[i] <- sum(x[i,])
    }
    z
}
sum2 <- function(x) {
    apply(x, 1, sum)
}
sum3 <- function(x) {
    rowSums(x)
}
microbenchmark(sum1(x),sum1a(x),sum2(x),sum3(x))

#2b
n<-1e5
z <- rnorm(n)
neg1 <- function(n,z){
    zneg <- 0
    j <- 1
    for (i in 1:n) {
        if (z[i]<0) {
            zneg[j] <- z[i]
            j <- j+1
        }
    }
    zneg
}
neg2 <- function(n,z){
    z[z<0]
}
neg2a <- function(z){
    z[z<0]
}
microbenchmark(neg1(n,z), neg2(n,z), neg2a(z))

#3
set.seed(1)
n<-1e3
A <- matrix(runif(n*n),n,n)
x <- runif(n)
t(x)%*%A%*%x
sum(diag(A))
sum(diag(t(A)%*%diag(x)%*%A))

#4a
n <- 1e3
A <- matrix(runif(n*n),n,n)
xx <- runif(n)
y <- A%*%xx
#4b
Ainv <- solve(A)
x <- Ainv%*%y
mean(abs(x-xx))
#4c
x <- solve(A,y)
mean(abs(x-xx))

#5a
empcdf1 <- function(x) { #slow version, for reference
    cdf <-numeric(length(x))
    for (i in 1:length(x)) {
        cdf[i] <- sum(x<x[i])
    }
    cdf/length(x)
}
empcdf <- function(x) { #fast version
    (sort.int(sort.int(x, index.return=TRUE)$ix, index.return = T)$ix - 1)/length(x)
}
#test:
x <- rnorm(10000)
plot(x, empcdf(x))
plot(sort.int(x), sort.int(empcdf(x)), type='s')

#5b
empcdf <- function(x, plot.cdf=TRUE, ...) { 
    cdf <- (sort.int(sort.int(x, index.return=TRUE)$ix, index.return = T)$ix - 1)/length(x)
    if (plot.cdf) {plot(sort.int(x), sort.int(cdf), type='s', xlab='', ylab='', main="Empirical CDF", ...)}
    else {return(cdf)}
}

#7a
rb <- function(x,z) {
    if (length(x) != length(z)){stop("arguments should be vectors of the same length")}
    100*(z-x^2)^2 + (1-x)^2
}
#7b
s1<-seq(-1.5,1.5, 0.01)
s2<-seq(-0.5, 1.5, 0.01)
mat <- outer(s1, s2, FUN=rb)
contour(s1,s2,mat)
#7c
contour(s1,s2,log(mat, base=10))
#7d
rb.grad <- function(x,z) {
    dfx <- -400*x*(z-x^2) -2*(1-x)
    dfz <- 200*(z-x^2)
    c(dfx,dfz)
}
#7e
x<-1.8
z<-4.1
Delta <- 1e-7
diffx <- (rb(x+Delta,z)-rb(x,z))/Delta
diffz <- (rb(x,z+Delta)-rb(x,z))/Delta
c(diffx,diffz)
rb.grad(x,z)
#7f
rb.hess <- function(x,z) {
    h11 <- -400*z + 1200*x^2 + 2
    h12 <- -400*x
    h22 <- 200
    matrix(c(h11,h12,h12,h22),2,2)
}
#7g
x<-1.8
z<-4.1
Delta <- 1e-7
diffx <- (rb.grad(x+Delta,z)-rb.grad(x,z))/Delta
diffz <- (rb.grad(x,z+Delta)-rb.grad(x,z))/Delta
rbind(diffx,diffz)
rb.hess(x,z)
#7h
