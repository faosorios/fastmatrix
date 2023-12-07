c ID: mchol_GMW.f, last updated 2023-12-06, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE mchol_dcmp(a, lda, n, d, macheps, info)
      INTEGER          lda, n, info
      DOUBLE PRECISION a(lda,*), d(*)
      DOUBLE PRECISION macheps
c
c     Computes the modified Cholesky factorization of a real symmetric but 
c     not necessarily positive definite matrix 'a'. The factorization has 
c     the form a + e = l %*% d %*% t(l), where 'd' is a diagonal matrix, 'l' 
c     is unitary lower triangular and 'e' is a nonnegative diagonal matrix 
c     that is zero if 'a' es sufficiently positive definite.
c     This routine is based in algorithm of Gill, Murray and Wright (1981). 
c     Practical Optimization. Academic Press.
c
c     parameters:
c     a       (input/output) DOUBLE PRECISION array, dimension (lda, n)
c             on entry, the symmetric matrix a. the leading n-by-n upper
c             triangular part of a contains the upper triangular part of
c             the matrix a. on exit, if info = 0, the factor 'l' from the
c             modified Cholesky factorization is stored on the strictly 
c             lower triangular part of a.
c     lda     (input) INTEGER
c             lda is the leading dimension of the array a. lda >= max(1,n).
c     n       (input) INTEGER
c             the order of the matrix a. n >= 0.
c     d       (output) DOUBLE PRECISION array, dimension (n)
c             contains the diagonal elements of matrix 'd'.
c     macheps (input) DOUBLE PRECISION
c             the smallest positive floating-point number 'macheps' such 
c             as 1 + macheps != 1, normally 2.220446e-16.
c     info    (output) INTEGER
c             = 0: successful exit
c             < 0: if info = -i, the i-th argument had an illegal value
c
c     .. parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER       (ZERO = 0.0d0, ONE = 1.0d0)
c     .. local scalars ..
      INTEGER          i, j, k
      DOUBLE PRECISION accum, beta, delta, gamma, eps, theta, tmp
c     .. intrinsic ..
      INTRINSIC        abs, max, sqrt
c
c     test the input parameters
c
      info = 0
      if (n .LT. 0) then
        info = -3
      else if (lda .LT. max(1, n)) then
        info = -2
      else if (macheps .LT. ZERO) then
        info = -5
      end if
      if (info .NE. 0) return
c
c     quick return if possible
c
      if (n .EQ. 0) return
c
c     computing largest-magnitude diagonal and off-diagonal 
c     elements of the matrix a
c
      gamma = ZERO
      do i = 1, n
        gamma = max(gamma, abs(a(i,i)))
      end do
      eps = ZERO
      do j = 1, n
        do i = j+1, n
          eps = max(eps, abs(a(i,j)))
        end do
      end do 
c
c     computing delta and beta parameters
c
      delta = macheps * max(gamma + eps, ONE)
      beta = max(gamma, eps / sqrt(n * n - ONE))
      beta = sqrt(max(beta, macheps))
c
c     initialize diagonal elements
c
      do i = 1, n
        d(i) = a(i,i)
      end do
c
c     main loop
c
      do j = 1, n
        do k = 1, j-1
          a(j,k) = a(j,k) / d(k)
        end do
        do i = j+1, n
          accum = ZERO
          do k = 1, j-1
            accum = accum + a(i,k) * a(j,k)
          end do
          a(i,j) = a(i,j) - accum
        end do
        theta = ZERO
        if (j .LE. n) then
          do i = j+1, n
            theta = max(theta, a(i,j))
          end do
        end if
        tmp = max(delta, abs(d(j)))
        d(j) = max(tmp, (theta / beta)**2)
        if (j .LT. n) then
          do i = j+1, n
            d(i) = d(i) - a(i,j)**2 / d(j)
          end do
        end if
      end do
c
      return
      END
