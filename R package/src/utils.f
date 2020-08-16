c ID: utils.f, last updated 2020-08-16, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE inner_frobenius(a, lda, b, ldb, n, p, value)
      INTEGER          lda, ldb, n, p
      DOUBLE PRECISION a(lda,*), b(ldb,*), value
c
c     computes the Frobenius inner product
c
c     .. parameters ..
      INTEGER          ONE
      PARAMETER       (ONE = 1)
c     .. BLAS functions ..
      DOUBLE PRECISION ddot
c     .. local scalars ..
      INTEGER          j
      DOUBLE PRECISION temp
c
      value = 0.0d0
c
      if (n .LE. 0) return
      if (p .LE. 0) return
      if (lda .LT. max(1,n)) return
      if (ldb .LT. max(1,n)) return
c
      do j = 1, p
        temp  = ddot(n, a(1,j), ONE, b(1,j), ONE)
        value = value + temp
      end do
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE hadamard(x, y, n, prod)
      INTEGER          n
      DOUBLE PRECISION x(*), y(*), prod(*)
c
c     hadamard returns the element-wise product between vectors 'x' and 'y'
c     using unrolled loops
c
c     parameters:
c     x     (input) DOUBLE PRECISION array, dimension (n)
c           on entry, an n-dimensional vector, unchanged on exit.
c     y     (input) DOUBLE PRECISION array, dimension (n)
c           on entry, an n-dimensional vector, unchanged on exit.
c     n     (input) INTEGER
c           order of vectors x and y, n > 0.
c     prod  (output) DOUBLE PRECISION array, dimension (n)
c           on exit, prod contains the element-wise product between x and y.
c
      INTEGER i, m, mp1
c
      if (n .LE. 0) return
c
      m = mod(n, 5)
      if (m .EQ. 0) goto 20
      do i = 1, m
        prod(i) = x(i) * y(i)
      end do
      if (n .LT. 5) return
c
   20 mp1 = m + 1
      do i = mp1, n, 5
        prod(i) = x(i) * y(i)
        prod(i + 1) = x(i + 1) * y(i + 1)
        prod(i + 2) = x(i + 2) * y(i + 2)
        prod(i + 3) = x(i + 3) * y(i + 3)
        prod(i + 4) = x(i + 4) * y(i + 4)
      end do
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE sweepop(a, lda, p, k, work, reverse, info)
      INTEGER          lda, p, k, info
      DOUBLE PRECISION a(lda,*), work(*)
c
c     Gauss-Jordan sweep operator for symmetric matrices
c
c     parameters:
c     a       (input/output) DOUBLE PRECISION array, dimension (lda, p)
c             a symmetric matrix. on exit if info = 0, a is the sweeped matrix
c     lda     (input) INTEGER
c             lda is the leading dimension of the array a. lda >= max(1,p).
c     p       (input) INTEGER
c             p is the order of the matrix. k > 0.
c     k       (input) INTEGER
c             k element of the diagonal which will be sweeped, k > 0.
c     work    (workspace) DOUBLE PRECISION array, dimension (p)
c             working array.
c     reverse (input) INTEGER
c             reverse sweep operator must be applied?
c     info    (output) INTEGER
c             = 0:  successful exit
c             < 0:  if info = i, the i-th argument had an illegal value
c
      DOUBLE PRECISION alpha
c
      info = 0
      if (p .LE. 0) then
        info = 3
      else if (lda .LT. max(1, p)) then
        info = 2
      else if ((k .LE. 0) .OR. (k .GT. p)) then
        info = 4
      end if
c
      if (info .NE. 0) return
c
      if (a(k,k) .EQ. 0.d0) then
        info = -4
        return
      end if
c
      alpha = 1.0d0 / a(k,k)
      call dcopy(k, a(1,k), 1, work, 1)
      call dcopy(p - k, a(k,k+1), lda, work(k+1), 1)
      call dsyr('U', p, -alpha, work, 1, a, lda)
      if (reverse .EQ. 0) then
        call dscal(p, alpha, work, 1)
      else
        call dscal(p, -alpha, work, 1)
      end if
      call dcopy(k - 1, work, 1, a(1,k), 1)
      call dcopy(p - k, work(k+1), 1, a(k,k+1), lda)
      a(k,k) = -alpha
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE vech(x, ldx, n, y)
      INTEGER          ldx, n
      DOUBLE PRECISION x(ldx,*), y(*)
c
c     vech returns the element-wise product between vectors 'x' and 'y'
c     using unrolled loops
c
c     parameters:
c     x     (input) DOUBLE PRECISION array, dimension (ldx, n)
c           a square matrix x. only the upper half of x need be stored.
c           the lower part of the array x is not referenced.
c     ldx   (input) INTEGER
c           ldx is the leading dimension of the array x. ldx >= max(1,n).
c     n     (input) INTEGER
c           n is the order of the matrix. n > 0.
c     y     (output) DOUBLE PRECISION array, dimension (n * (n + 1) / 2).
c           on exit, y contains the vectorization of the upper part of
c           the square matrix x.
c
      INTEGER i, j, k
c
      if (n .LE. 0) return
      if (ldx .LT. max(1,n)) return
c
      k = 0
      do i = 1, n
        do j = i, n
          k = k + 1
          y(k) = x(i,j)
        end do
      end do
c
      return
      END
