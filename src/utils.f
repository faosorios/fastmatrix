c ID: utils.f, last updated 2020-08-08, F.Osorio

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
      m = mod(n, 4)
      if (m .EQ. 0) goto 20
      do i = 1, m
        prod(i) = x(i) * y(i)
      end do
      if (n .LT. 4) return
c
   20 mp1 = m + 1
      do i = mp1, n, 4
        prod(i) = x(i) * y(i)
        prod(i + 1) = x(i + 1) * y(i + 1)
        prod(i + 2) = x(i + 2) * y(i + 2)
        prod(i + 3) = x(i + 3) * y(i + 3)
      end do
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
