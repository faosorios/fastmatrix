c ID: omo.f, last updated 2022-07-06, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION FUNCTION blinf(a, lda, n, p, x, y)
      INTEGER          lda, n, p
      DOUBLE PRECISION a(lda,*), x(*), y(*)
c
c     this function computes the bilinear form, t(x) %*% A %*% y
c
c     parameters:
c     a       DOUBLE PRECISION array, dimension (lda, p)
c             a rectangular matrix.
c     lda     INTEGER
c             lda is the leading dimension of the array a. lda >= max(1,n).
c     n       INTEGER
c             the number of rows of the matrix a. n > 0.
c     p       INTEGER
c             the number of columns of the matrix a. p > 0.
c     x       DOUBLE PRECISION array, dimension (n)
c             vector of length n.
c     y       DOUBLE PRECISION array, dimension (p)
c             vector of length p.
c
c     .. local scalars ..
      INTEGER          i, j
      DOUBLE PRECISION accum
c     .. parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER       (ZERO = 0.0d0)
c
      blinf = ZERO
c
c     quick return if possible
c
      if (n .LE. 0) return
      if (p .LE. 0) return
      if (lda .LT. max(1, n)) return
c
c     start the operations
c
      accum = ZERO
      do i = 1, n
        do j = 1, p
          accum = accum + a(i,j) * x(i) * y(j)
        end do
      end do
      blinf = accum
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION FUNCTION quadf(a, lda, n, x)
      INTEGER          lda, n
      DOUBLE PRECISION a(lda,*), x(*)
c
c     this function computes the quadratic form, t(x) %*% A %*% x
c
c     parameters:
c     a       DOUBLE PRECISION array, dimension (lda, n)
c             an n-by-n square matrix.
c     lda     INTEGER
c             lda is the leading dimension of the array a. lda >= max(1,n).
c     n       INTEGER
c             the number of rows of the matrix a. n > 0.
c     x       DOUBLE PRECISION array, dimension (n)
c             vector of length n.
c
c     .. local scalars ..
      INTEGER          i, j
      DOUBLE PRECISION accum
c     .. parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER       (ZERO = 0.0d0)
c
      quadf = ZERO
c
c     quick return if possible
c
      if (n .LE. 0) return
      if (lda .LT. max(1, n)) return
c
c     start the operations
c
      accum = ZERO
      do i = 1, n
        do j = 1, n
          accum = accum + a(i,j) * x(i) * x(j)
        end do
      end do
      quadf = accum
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE murrv(y, a, lda, n, p, x, info)
      INTEGER          lda, n, p, info
      DOUBLE PRECISION a(lda,*), x(*), y(*)
c
c     multiplies a real rectangular matrix by a vector, y = a %*% x
c
c     parameters:
c     a       (input) DOUBLE PRECISION array, dimension (lda, p)
c             a rectangular matrix.
c     lda     (input) INTEGER
c             lda is the leading dimension of the array a. lda >= max(1,n).
c     n       (input) INTEGER
c             the number of rows of the matrix a. n > 0.
c     p       (input) INTEGER
c             the number of columns of the matrix a. p > 0.
c     x       (input) DOUBLE PRECISION array, dimension (p)
c             vector of length p.
c     y       (output) DOUBLE PRECISION array, dimension (n)
c             vector of length n.
c     info    (output) INTEGER
c             = 0: successful exit
c             < 0: if info = -i, the i-th argument had an illegal value
c
c     .. local scalars ..
      INTEGER          i, j
      DOUBLE PRECISION accum
c     .. parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER       (ZERO = 0.0d0)
c
c     test the input parameters
c
      info = 0
      if (n .LT. 0) then
        info = -3
      else if (p .LT. 0) then
        info = -4
      else if (lda .LT. max(1, n)) then
        info = -2
      end if
      if (info .NE. 0) return
c
c     quick return if possible
c
      if ((n .EQ. 0) .OR. (p .EQ. 0)) return
c
c     start the operations
c
      do i = 1, n
        accum = ZERO
        do j = 1, p
          accum = accum + a(i,j) * x(j)
        end do
        y(i) = accum
      end do
c
      return
      END
