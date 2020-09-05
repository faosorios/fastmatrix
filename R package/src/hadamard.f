c ID: hadamard.f, last updated 2020-09-03, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE hadamard_prod(x, y, n, prod)
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
c     quick return if possible
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
