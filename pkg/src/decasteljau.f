c ID: decasteljau.f, last updated 2021-11-12, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE decasteljau(x, y, n, s, z)
      INTEGER          n
      DOUBLE PRECISION x(*), y(*), s, z(*)
c
c     returns a point on the Bezier curve using the De Casteljau's method
c
c     parameters:
c     x, y    (input) DOUBLE PRECISION array, dimension (n)
c             coordinates of the points in the scatter plot
c     n       (input) INTEGER
c             the number of observations. n > 0.
c     s       (input) DOUBLE PRECISION
c             value in the interval [0,1]
c     z       (output) DOUBLE PRECISION array, dimension (2)
c             on exit the Bezier curve at point 's'
c
c     .. parameters ..
      DOUBLE PRECISION ONE
      PARAMETER       (ONE = 1.0d0)
c     .. local scalars ..
      DOUBLE PRECISION p, q
      INTEGER          i, k
c     .. local arrays ..
      DOUBLE PRECISION u(n), v(n)
c
c     initializing
c
      p = s
      q = ONE - p
c
c     save input
c
      do i = 1, n
        u(i) = x(i)
        v(i) = y(i)
      end do
c
c     the De Casteljau's recurrence
c
      do k = 2, n
        do i = 1, n - k + 1
          u(i) = q * u(i) + p * u(i+1)
          v(i) = q * v(i) + p * v(i+1)
        end do
      end do
c
c     output
c
      z(1) = u(1)
      z(2) = v(1)
c
      return
      END
