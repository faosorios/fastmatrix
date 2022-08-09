pkgname <- "fastmatrix"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('fastmatrix')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Mahalanobis")
### * Mahalanobis

flush(stderr()); flush(stdout())

### Name: Mahalanobis
### Title: Mahalanobis distance
### Aliases: Mahalanobis
### Keywords: multivariate

### ** Examples

x <- cbind(1:6, 1:3)
xbar <- colMeans(x)
S <- matrix(c(1,4,4,1), ncol = 2) # is negative definite
D2 <- mahalanobis(x, center = xbar, S)
all(D2 >= 0) # several distances are negative

## next command produces the following error:
## Covariance matrix is possibly not positive-definite
## Not run: D2 <- Mahalanobis(x, center = xbar, S)



cleanEx()
nameEx("array.mult")
### * array.mult

flush(stderr()); flush(stdout())

### Name: array.mult
### Title: Array multiplication
### Aliases: array.mult
### Keywords: array algebra

### ** Examples

x <- array(0, dim = c(2,3,3)) # 2 x 3 x 3 array
x[,,1] <- c(1,2,2,4,3,6)
x[,,2] <- c(2,4,4,8,6,12)
x[,,3] <- c(3,6,6,12,9,18)

a <- matrix(1, nrow = 2, ncol = 3)
b <- matrix(1, nrow = 3, ncol = 2)

y <- array.mult(a, b, x) # a 2 x 2 x 2 array
y



cleanEx()
nameEx("asSymmetric")
### * asSymmetric

flush(stderr()); flush(stdout())

### Name: asSymmetric
### Title: Force a matrix to be symmetric
### Aliases: asSymmetric
### Keywords: array

### ** Examples

a <- matrix(1:16, ncol = 4)
isSymmetric(a) # FALSE
a <- asSymmetric(a) # copy lower triangle into upper triangle



cleanEx()
nameEx("bezier")
### * bezier

flush(stderr()); flush(stdout())

### Name: bezier
### Title: Computation of Bezier curve
### Aliases: bezier
### Keywords: smooth

### ** Examples

# a tiny example
x <- c(1.0, 0.25, 1.25, 2.5, 4.00, 5.0)
y <- c(0.5, 2.00, 3.75, 4.0, 3.25, 1.0)
plot(x, y, type = "o")
z <- bezier(x, y, ngrid = 50)
lines(z$xgrid, z$ygrid, lwd = 2, lty = 2, col = "red")

# other simple example
x <- c(4,6,4,5,6,7)
y <- 1:6
plot(x, y, type = "o")
z <- bezier(x, y, ngrid = 50)
lines(z$xgrid, z$ygrid, lwd = 2, lty = 2, col = "red")



cleanEx()
nameEx("bracket.prod")
### * bracket.prod

flush(stderr()); flush(stdout())

### Name: bracket.prod
### Title: Bracket product
### Aliases: bracket.prod
### Keywords: array algebra

### ** Examples

x <- array(0, dim = c(2,3,3)) # 2 x 3 x 3 array
x[,,1] <- c(1,2,2,4,3,6)
x[,,2] <- c(2,4,4,8,6,12)
x[,,3] <- c(3,6,6,12,9,18)

a <- matrix(1, nrow = 3, ncol = 2)

y <- bracket.prod(a, x) # a 3 x 3 x 3 array
y



cleanEx()
nameEx("cg")
### * cg

flush(stderr()); flush(stdout())

### Name: cg
### Title: Solve linear systems using the conjugate gradients method
### Aliases: cg
### Keywords: algebra array

### ** Examples

a <- matrix(c(4,3,0,3,4,-1,0,-1,4), ncol = 3)
b <- c(24,30,-24)
z <- cg(a, b)
z # converged in 3 iterations



cleanEx()
nameEx("cholupdate")
### * cholupdate

flush(stderr()); flush(stdout())

### Name: cholupdate
### Title: Rank 1 update to Cholesky factorization
### Aliases: cholupdate
### Keywords: algebra array

### ** Examples

a <- matrix(c(1,1,1,1,2,3,1,3,6), ncol = 3)
r <- chol(a)
x <- c(0,0,1)
b <- a + outer(x,x)
r1 <- cholupdate(r, x)
r1
all(r1 == chol(b)) # TRUE



cleanEx()
nameEx("circulant")
### * circulant

flush(stderr()); flush(stdout())

### Name: circulant
### Title: Form a symmetric circulant matrix
### Aliases: circulant
### Keywords: array algebra

### ** Examples

x <- c(2,3,5,7,11,13)
circulant(x)



cleanEx()
nameEx("comm.info")
### * comm.info

flush(stderr()); flush(stdout())

### Name: comm.info
### Title: Compact information to construct the commutation matrix
### Aliases: comm.info
### Keywords: array

### ** Examples

z <- comm.info(m = 3, n = 2, condensed = FALSE)
z # where are the ones in commutation matrix of order '3,2'?

K32 <- commutation(m = 3, n = 2, matrix = TRUE)
K32 # only recommended if m and n are very small



cleanEx()
nameEx("comm.prod")
### * comm.prod

flush(stderr()); flush(stdout())

### Name: comm.prod
### Title: Matrix multiplication envolving the commutation matrix
### Aliases: comm.prod
### Keywords: array algebra

### ** Examples

K42 <- commutation(m = 4, n = 2, matrix = TRUE)
x <- matrix(1:24, ncol = 3)
y <- K42 %*% x

z <- comm.prod(m = 4, n = 2, x) # K42 is not stored
all(z == y) # matrices y and z are equal!



cleanEx()
nameEx("commutation")
### * commutation

flush(stderr()); flush(stdout())

### Name: commutation
### Title: Commutation matrix
### Aliases: commutation
### Keywords: array algebra

### ** Examples

z <- commutation(m = 100, condensed = TRUE)
object.size(z) # 40.6 Kb of storage

z <- commutation(m = 100, condensed = FALSE)
object.size(z) # 80.7 Kb of storage

K100 <- commutation(m = 100, matrix = TRUE) # time: < 2 secs
object.size(K100) # 400 Mb of storage, do not request this matrix!

# a small example
K32 <- commutation(m = 3, n = 2, matrix = TRUE)
a <- matrix(1:6, ncol = 2)
v <- K32 %*% vec(a)
all(vec(t(a)) == as.vector(v)) # vectors are equal!



cleanEx()
nameEx("corAR1")
### * corAR1

flush(stderr()); flush(stdout())

### Name: corAR1
### Title: AR(1) Correlation Structure
### Aliases: corAR1
### Keywords: array

### ** Examples

R <- corAR1(rho = 0.8, p = 5)



cleanEx()
nameEx("corCS")
### * corCS

flush(stderr()); flush(stdout())

### Name: corCS
### Title: Compound Symmetry Correlation Structure
### Aliases: corCS
### Keywords: array

### ** Examples

R <- corCS(rho = 0.8, p = 5)



cleanEx()
nameEx("cov.MSSD")
### * cov.MSSD

flush(stderr()); flush(stdout())

### Name: cov.MSSD
### Title: Mean Square Successive Difference (MSSD) estimator of the
###   covariance matrix
### Aliases: cov.MSSD
### Keywords: multivariate

### ** Examples

x <- cbind(1:10, c(1:3, 8:5, 8:10))
z0 <- cov(x)
z0
z1 <- cov.MSSD(x)
z1



cleanEx()
nameEx("cov.weighted")
### * cov.weighted

flush(stderr()); flush(stdout())

### Name: cov.weighted
### Title: Weighted covariance matrices
### Aliases: cov.weighted
### Keywords: multivariate

### ** Examples

x <- cbind(1:10, c(1:3, 8:5, 8:10))
z0 <- cov.weighted(x) # all weights are 1
D2 <- Mahalanobis(x, center = z0$mean, cov = z0$cov)
p <- ncol(x)
wts <- (p + 1) / (1 + D2) # nice weights!
z1 <- cov.weighted(x, weights = wts)
z1



cleanEx()
nameEx("dupl.cross")
### * dupl.cross

flush(stderr()); flush(stdout())

### Name: dupl.cross
### Title: Matrix crossproduct envolving the duplication matrix
### Aliases: dupl.cross
### Keywords: array algebra

### ** Examples

D2 <- duplication(n = 2, matrix = TRUE)
D3 <- duplication(n = 3, matrix = TRUE)
x <- matrix(1, nrow = 9, ncol = 4)
y <- t(D3) %*% x %*% D2

z <- dupl.cross(n = 3, k = 2, x) # D2 and D3 are not stored
all(z == y) # matrices y and z are equal!

x <- matrix(1, nrow = 9, ncol = 9)
z <- dupl.cross(n = 3, x = x) # same matrix is used to pre- and post-multiplying x
z # print result



cleanEx()
nameEx("dupl.info")
### * dupl.info

flush(stderr()); flush(stdout())

### Name: dupl.info
### Title: Compact information to construct the duplication matrix
### Aliases: dupl.info
### Keywords: array

### ** Examples

z <- dupl.info(n = 3, condensed = FALSE)
z # where are the ones in duplication of order 3?

D3 <- duplication(n = 3, matrix = TRUE)
D3 # only recommended if n is very small



cleanEx()
nameEx("dupl.prod")
### * dupl.prod

flush(stderr()); flush(stdout())

### Name: dupl.prod
### Title: Matrix multiplication envolving the duplication matrix
### Aliases: dupl.prod
### Keywords: array algebra

### ** Examples

D4 <- duplication(n = 4, matrix = TRUE)
x <- matrix(1, nrow = 16, ncol = 2)
y <- crossprod(D4, x)

z <- dupl.prod(n = 4, x, transposed = TRUE) # D4 is not stored
all(z == y) # matrices y and z are equal!



cleanEx()
nameEx("duplication")
### * duplication

flush(stderr()); flush(stdout())

### Name: duplication
### Title: Duplication matrix
### Aliases: duplication
### Keywords: array algebra

### ** Examples

z <- duplication(n = 100, condensed = TRUE)
object.size(z) # 40.5 Kb of storage

z <- duplication(n = 100, condensed = FALSE)
object.size(z) # 80.6 Kb of storage

D100 <- duplication(n = 100, matrix = TRUE)
object.size(D100) # 202 Mb of storage, do not request this matrix!

# a small example
D3 <- duplication(n = 3, matrix = TRUE)
a <- matrix(c( 1, 2, 3,
               2, 3, 4,
               3, 4, 5), nrow = 3)
upper <- vech(a)
v <- D3 %*% upper
all(vec(a) == as.vector(v)) # vectors are equal!



cleanEx()
nameEx("equilibrate")
### * equilibrate

flush(stderr()); flush(stdout())

### Name: equilibrate
### Title: Equilibration of a rectangular or symmetric matrix
### Aliases: equilibrate
### Keywords: array algebra

### ** Examples

x <- matrix(c(1, 1, 1,
              1, 2, 1,
              1, 3, 1,
              1, 1,-1,
              1, 2,-1,
              1, 3,-1), ncol = 3, byrow = TRUE)
z <- equilibrate(x)
apply(z, 2, function(x) sum(x^2)) # all 1

xx <- crossprod(x)
equilibrate(xx)



cleanEx()
nameEx("frank")
### * frank

flush(stderr()); flush(stdout())

### Name: frank
### Title: Frank matrix
### Aliases: frank
### Keywords: array algebra

### ** Examples

F5 <- frank(n = 5)
det(F5) # equals 1



cleanEx()
nameEx("geomean")
### * geomean

flush(stderr()); flush(stdout())

### Name: geomean
### Title: Geometric mean
### Aliases: geomean
### Keywords: univar

### ** Examples

set.seed(149)
x <- rlnorm(1000)
mean(x)    # 1.68169
median(x)  # 0.99663
geomean(x) # 1.01688



cleanEx()
nameEx("hadamard.prod")
### * hadamard.prod

flush(stderr()); flush(stdout())

### Name: hadamard
### Title: Hadamard product of two matrices
### Aliases: hadamard
### Keywords: array algebra

### ** Examples

x <- matrix(rep(1:10, times = 5), ncol = 5)
y <- matrix(rep(1:5, each = 10), ncol = 5)
z <- hadamard(x, y)
z



cleanEx()
nameEx("harris.test")
### * harris.test

flush(stderr()); flush(stdout())

### Name: harris.test
### Title: Test for variance homogeneity of correlated variables
### Aliases: harris.test
### Keywords: htest

### ** Examples

x <- iris[,1:4]
z <- harris.test(x, test = "robust")
z



cleanEx()
nameEx("helmert")
### * helmert

flush(stderr()); flush(stdout())

### Name: helmert
### Title: Helmert matrix
### Aliases: helmert
### Keywords: array algebra

### ** Examples

n <- 1000
set.seed(149)
x <- rnorm(n)

H <- helmert(n)
object.size(H) # 7.63 Mb of storage
K <- H[2:n,]
z <- c(K %*% x)
sum(z^2) # 933.1736

# same that
(n - 1) * var(x)



cleanEx()
nameEx("is.lower.tri")
### * is.lower.tri

flush(stderr()); flush(stdout())

### Name: is.lower.tri
### Title: Check if a matrix is lower or upper triangular
### Aliases: is.lower.tri is.upper.tri
### Keywords: array

### ** Examples

  x <- matrix(rnorm(10 * 3), ncol = 3)
  R <- chol(crossprod(x))

  is.lower.tri(R)
  is.upper.tri(R)



cleanEx()
nameEx("jacobi")
### * jacobi

flush(stderr()); flush(stdout())

### Name: jacobi
### Title: Solve linear systems using the Jacobi method
### Aliases: jacobi
### Keywords: algebra array

### ** Examples

a <- matrix(c(5,-3,2,-2,9,-1,3,1,-7), ncol = 3)
b <- c(-1,2,3)
start <- c(1,1,1)
z <- jacobi(a, b, start)
z # converged in 15 iterations



cleanEx()
nameEx("kronecker.prod")
### * kronecker.prod

flush(stderr()); flush(stdout())

### Name: kronecker.prod
### Title: Kronecker product on matrices
### Aliases: kronecker.prod
### Keywords: array

### ** Examples

# block diagonal matrix:
a <- diag(1:3)
b <- matrix(1:4, ncol = 2)
kronecker.prod(a, b)

# examples with vectors
ones <- rep(1, 4)
y <- 1:3
kronecker.prod(ones, y) # 12-dimensional vector
kronecker.prod(ones, t(y)) # 3 x 3 matrix



cleanEx()
nameEx("krylov")
### * krylov

flush(stderr()); flush(stdout())

### Name: krylov
### Title: Computes a Krylov matrix
### Aliases: krylov
### Keywords: array

### ** Examples

a <- matrix(c(1, 3, 2, -5, 1, 7, 1, 5, -4), ncol = 3, byrow = TRUE)
b <- c(1, 1, 1)
k <- krylov(a, b, m = 4)
k



cleanEx()
nameEx("kurtosis")
### * kurtosis

flush(stderr()); flush(stdout())

### Name: kurtosis
### Title: Mardia's multivariate skewness and kurtosis coefficients
### Aliases: kurtosis skewness
### Keywords: multivariate

### ** Examples

setosa <- iris[1:50,1:4]
kurtosis(setosa)
skewness(setosa)



cleanEx()
nameEx("ldl")
### * ldl

flush(stderr()); flush(stdout())

### Name: ldl
### Title: The LDL decomposition
### Aliases: ldl
### Keywords: algebra array

### ** Examples

a <- matrix(c(2,-1,0,-1,2,-1,0,-1,1), ncol = 3)
z <- ldl(a)
z # information of LDL factorization

# computing det(a)
prod(z$d) # product of diagonal elements of D

# a non-positive-definite matrix
m <- matrix(c(5,-5,-5,3), ncol = 2)
try(chol(m)) # fails
ldl(m)



cleanEx()
nameEx("lu-methods")
### * lu-methods

flush(stderr()); flush(stdout())

### Name: lu-methods
### Title: Reconstruct the L, U, or X Matrices from an LU object
### Aliases: constructX extractL extractU
### Keywords: algebra array

### ** Examples

a <- matrix(c(10,-3,5,-7,2,-1,0,6,5), ncol = 3)
z <- lu(a)
L <- extractL(z)
L
U <- extractU(z)
U
X <- constructX(z)
all(a == X)



cleanEx()
nameEx("lu")
### * lu

flush(stderr()); flush(stdout())

### Name: lu
### Title: The LU factorization of a square matrix
### Aliases: lu lu.default solve.lu is.lu
### Keywords: algebra array

### ** Examples

a <- matrix(c(3,2,6,17,4,18,10,-2,-12), ncol = 3)
z <- lu(a)
z # information of LU factorization

# computing det(a)
prod(diag(z$lu)) # product of diagonal elements of U

# solve linear equations
b <- matrix(1:6, ncol = 2)
solve(z, b)



cleanEx()
nameEx("lu2inv")
### * lu2inv

flush(stderr()); flush(stdout())

### Name: lu2inv
### Title: Inverse from LU factorization
### Aliases: lu2inv
### Keywords: algebra array

### ** Examples

a <- matrix(c(3,2,6,17,4,18,10,-2,-12), ncol = 3)
z <- lu(a)
a %*% lu2inv(z)



cleanEx()
nameEx("matrix.inner")
### * matrix.inner

flush(stderr()); flush(stdout())

### Name: matrix.inner
### Title: Compute the inner product between two rectangular matrices
### Aliases: matrix.inner
### Keywords: array math

### ** Examples

x <- matrix(c(1, 1, 1,
              1, 2, 1,
              1, 3, 1,
              1, 1,-1,
              1, 2,-1,
              1, 3,-1), ncol = 3, byrow = TRUE)
y <- matrix(1, nrow = 6, ncol = 3)
matrix.inner(x, y)

# must be equal
matrix.norm(x, type = "Frobenius")^2
matrix.inner(x)



cleanEx()
nameEx("matrix.norm")
### * matrix.norm

flush(stderr()); flush(stdout())

### Name: matrix.norm
### Title: Compute the norm of a rectangular matrix
### Aliases: matrix.norm
### Keywords: array math

### ** Examples

# a tiny example
x <- matrix(c(1, 1, 1,
              1, 2, 1,
              1, 3, 1,
              1, 1,-1,
              1, 2,-1,
              1, 3,-1), ncol = 3, byrow = TRUE)
matrix.norm(x, type = "Frobenius")
matrix.norm(x, type = "1")
matrix.norm(x, type = "Inf")

# an example not that small
n <- 1000
x <- .5 * diag(n) + 0.5 * matrix(1, nrow = n, ncol = n)
matrix.norm(x, type = "Frobenius")
matrix.norm(x, type = "1")
matrix.norm(x, type = "Inf")
matrix.norm(x, type = "maximum") # equal to 1



cleanEx()
nameEx("mediancenter")
### * mediancenter

flush(stderr()); flush(stdout())

### Name: mediancenter
### Title: Mediancenter
### Aliases: mediancenter
### Keywords: multivariate

### ** Examples

x <- cbind(1:10, c(1:3, 8:5, 8:10))
z <- mediancenter(x)$median # degenerate solution
xbar <- colMeans(x)
plot(x, xlab = "", ylab = "")
points(x = xbar[1], y = xbar[2], pch = 16, col = "red")
points(x = z[1], y = z[2], pch = 3, col = "blue", lwd = 2)



cleanEx()
nameEx("minkowski")
### * minkowski

flush(stderr()); flush(stdout())

### Name: minkowski
### Title: Computes the p-norm of a vector
### Aliases: minkowski
### Keywords: math

### ** Examples

# a tiny example
x <- rnorm(1000)
minkowski(x, p = 1)
minkowski(x, p = 1.5)
minkowski(x, p = 2)
minkowski(x, p = Inf)

x <- x / minkowski(x)
minkowski(x, p = 2) # equal to 1



cleanEx()
nameEx("moments")
### * moments

flush(stderr()); flush(stdout())

### Name: moments
### Title: Central moments
### Aliases: moments
### Keywords: univar

### ** Examples

set.seed(149)
x <- rnorm(1000)
z <- moments(x)
z



cleanEx()
nameEx("ols")
### * ols

flush(stderr()); flush(stdout())

### Name: ols
### Title: Fit linear regression model
### Aliases: ols
### Keywords: regression

### ** Examples

# tiny example of regression
y <- c(1, 3, 3, 2, 2, 1)
x <- matrix(c(1, 1,
              2, 1,
              3, 1,
              1,-1,
              2,-1,
              3,-1), ncol = 2, byrow = TRUE)
f0 <- ols(y ~ x) # intercept is included by default
f0 # printing results (QR method was used)

f1 <- ols(y ~ x, method = "svd") # using SVD method instead
f1



cleanEx()
nameEx("ols.fit-methods")
### * ols.fit-methods

flush(stderr()); flush(stdout())

### Name: ols.fit-methods
### Title: Fit a Linear Model
### Aliases: ols.fit.cg ols.fit.chol ols.fit.qr ols.fit.svd ols.fit.sweep
### Keywords: regression array

### ** Examples

set.seed(151)
n <- 100
p <- 2
x <- matrix(rnorm(n * p), n, p) # no intercept!
y <- rnorm(n)
z <- ols.fit.chol(x, y)
z



cleanEx()
nameEx("ols.fit")
### * ols.fit

flush(stderr()); flush(stdout())

### Name: ols.fit
### Title: Fitter Functions for Linear Models
### Aliases: ols.fit
### Keywords: regression array

### ** Examples

set.seed(151)
n <- 100
p <- 2
x <- matrix(rnorm(n * p), n, p) # no intercept!
y <- rnorm(n)
fm <- ols.fit(x = x, y = y, method = "chol")
fm



cleanEx()
nameEx("power.method")
### * power.method

flush(stderr()); flush(stdout())

### Name: power.method
### Title: Power method to approximate dominant eigenvalue and eigenvector
### Aliases: power.method
### Keywords: array algebra

### ** Examples

n <- 1000
x <- .5 * diag(n) + 0.5 * matrix(1, nrow = n, ncol = n)

# dominant eigenvalue must be (n + 1) / 2
z <- power.method(x, only.value = TRUE)



cleanEx()
nameEx("ridge")
### * ridge

flush(stderr()); flush(stdout())

### Name: ridge
### Title: Ridge regression
### Aliases: ridge
### Keywords: models

### ** Examples

z <- ridge(GNP.deflator ~ ., data = longley, lambda = 4, method = "grid")
z # ridge regression on a grid over seq(0, 4, length = 200)

z <- ridge(GNP.deflator ~ ., data = longley)
z # ridge parameter selected using GCV (default)



cleanEx()
nameEx("seidel")
### * seidel

flush(stderr()); flush(stdout())

### Name: seidel
### Title: Solve linear systems using the Gauss-Seidel method
### Aliases: seidel
### Keywords: algebra array

### ** Examples

a <- matrix(c(5,-3,2,-2,9,-1,3,1,-7), ncol = 3)
b <- c(-1,2,3)
start <- c(1,1,1)
z <- seidel(a, b, start)
z # converged in 10 iterations



cleanEx()
nameEx("sherman.morrison")
### * sherman.morrison

flush(stderr()); flush(stdout())

### Name: sherman.morrison
### Title: Sherman-Morrison formula
### Aliases: sherman.morrison
### Keywords: array algebra

### ** Examples

n <- 10
ones <- rep(1, n)
a <- 0.5 * diag(n)
z <- sherman.morrison(a, ones, 0.5 * ones)
z



cleanEx()
nameEx("sweep.operator")
### * sweep.operator

flush(stderr()); flush(stdout())

### Name: sweep.operator
### Title: Gauss-Jordan sweep operator for symmetric matrices
### Aliases: sweep.operator
### Keywords: array algebra

### ** Examples

# tiny example of regression, last column contains 'y'
xy <- matrix(c(1, 1, 1, 1,
               1, 2, 1, 3,
               1, 3, 1, 3,
               1, 1,-1, 2,
               1, 2,-1, 2,
               1, 3,-1, 1), ncol = 4, byrow = TRUE)
z <- crossprod(xy)
z <- sweep.operator(z, k = 1:3)
cf <- z[1:3,4] # regression coefficients
RSS <- z[4,4]  # residual sum of squares

# an example not that small
x <- matrix(rnorm(1000 * 100), ncol = 100)
xx <- crossprod(x)
z <- sweep.operator(xx, k = 1)



cleanEx()
nameEx("symm.info")
### * symm.info

flush(stderr()); flush(stdout())

### Name: symm.info
### Title: Compact information to construct the symmetrizer matrix
### Aliases: symm.info
### Keywords: array

### ** Examples

z <- symm.info(n = 3)
z # elements in symmetrizer matrix of order 3

N3 <- symmetrizer(n = 3, matrix = TRUE)
N3 # only recommended if n is very small



cleanEx()
nameEx("symm.prod")
### * symm.prod

flush(stderr()); flush(stdout())

### Name: symm.prod
### Title: Matrix multiplication envolving the symmetrizer matrix
### Aliases: symm.prod
### Keywords: array algebra

### ** Examples

N4 <- symmetrizer(n = 4, matrix = TRUE)
x <- matrix(1:32, ncol = 2)
y <- N4 %*% x

z <- symm.prod(n = 4, x) # N4 is not stored
all(z == y) # matrices y and z are equal!



cleanEx()
nameEx("symmetrizer")
### * symmetrizer

flush(stderr()); flush(stdout())

### Name: symmetrizer
### Title: Symmetrizer matrix
### Aliases: symmetrizer
### Keywords: array algebra

### ** Examples

z <- symmetrizer(n = 100)
object.size(z) # 319 Kb of storage

N100 <- symmetrizer(n = 100, matrix = TRUE) # time: < 2 secs
object.size(N100) # 800 Mb of storage, do not request this matrix!

# a small example
N3 <- symmetrizer(n = 3, matrix = TRUE)
a <- matrix(rep(c(2,4,6), each = 3), ncol = 3)
a
b <- 0.5 * (a + t(a))
b
v <- N3 %*% vec(a)
all(vec(b) == as.vector(v)) # vectors are equal!



cleanEx()
nameEx("vec")
### * vec

flush(stderr()); flush(stdout())

### Name: vec
### Title: Vectorization of a matrix
### Aliases: vec
### Keywords: array

### ** Examples

x <- matrix(rep(1:10, each = 10), ncol = 10)
x
y <- vec(x)
y



cleanEx()
nameEx("vech")
### * vech

flush(stderr()); flush(stdout())

### Name: vech
### Title: Vectorization the lower triangular part of a square matrix
### Aliases: vech
### Keywords: array

### ** Examples

x <- matrix(rep(1:10, each = 10), ncol = 10)
x
y <- vech(x)
y



cleanEx()
nameEx("whitening")
### * whitening

flush(stderr()); flush(stdout())

### Name: whitening
### Title: Whitening transformation
### Aliases: whitening
### Keywords: multivariate

### ** Examples

x <- iris[,1:4]
species <- iris[,5]
pairs(x, col = species) # plot of Iris

# whitened data
z <- whitening(x)
pairs(z, col = species) # plot of



cleanEx()
nameEx("wilson.hilferty")
### * wilson.hilferty

flush(stderr()); flush(stdout())

### Name: wilson.hilferty
### Title: Wilson-Hilferty transformation
### Aliases: wilson.hilferty
### Keywords: multivariate

### ** Examples

x <- iris[,1:4]
z <- wilson.hilferty(x)
par(pty = "s")
qqnorm(z, main = "Transformed distances Q-Q plot")
abline(c(0,1), col = "red", lwd = 2, lty = 2)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
