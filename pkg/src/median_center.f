c ID: median_center.f, last updated 10-16-2021, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE median_center(x, ldx, n, p, median, iter, info)
      INTEGER          ldx, n, p, iter, info
      DOUBLE PRECISION x(ldx,*), median(*)
c
c     Computes the mediancenter for a sample of multivariate observations,
c     AS 78: Applied Statistics 23, 1974, 466-470. doi: 10.2307/2347150
c
c     parameters:
c     x       (input) DOUBLE PRECISION array, dimension (ldx, p)
c             a rectangular matrix.
c     ldx     (input) INTEGER
c             ldx is the leading dimension of the array x. ldx >= max(1,n).
c     n       (input) INTEGER
c             the number of rows of the matrix x. n > 0.
c     p       (input) INTEGER
c             the number of columns of the matrix x. p > 0.
c     median  (output) DOUBLE PRECISION array, dimension(p)
c             the coordinates of the mediancenter
c     iter    (output) INTEGER
c             the number of iterations performed, negative if a degenerate
c             solution is found
c     info    (output) INTEGER
c             = 0: successful exit
c             < 0: if info = -i, the i-th argument had an illegal value
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     .. parameters ..
      DOUBLE PRECISION ZERO, ONE, LEPSD, LEPSR, LEPSI
      PARAMETER       (ZERO = 0.0d0, ONE = 1.0d0)
      PARAMETER       (LEPSD = 1.0d-4, LEPSR = 1.0d-5, LEPSI = 1.0d-6)
      INTEGER          ICOUNT, LCOUNT
      PARAMETER       (ICOUNT = 200, LCOUNT = 100)
c     .. local scalars ..
      INTEGER          i, j, k, ii, l, lc, ll
      DOUBLE PRECISION accum, comp, corner, d, dd, diam, delta, epsd,
     *                 epsi, epsr, lambda, slam, u1, u2, xx
c     .. local arrays ..
      DOUBLE PRECISION c(p), z(p)
*     .. intrinsic functions ..
      INTRINSIC        max
c
c     test the input parameters
c
      info = 0
      if (n .LT. 0) then
        info = -3
      else if (p .LT. 0) then
        info = -4
      else if (ldx .LT. max(1, n)) then
        info = -2
      end if
      if (info .NE. 0) return
c
c     quick return if possible
c
      if ((n .EQ. 0) .OR. (p .EQ. 0)) then
        info = 1
        iter = 0
      end if
      if (info .NE. 0) return
c
c     initial settings
c
      ll = 0
      ii = 1
      if (n .EQ. 1) goto 25
c
c     calculate the diameter
c
      diam = 0.0d0
      do i = 2, n
        do j = 1, i-1
          accum = ZERO
          do k = 1, p
            accum = accum + (x(i,k) - x(j,k))**2
          end do
          diam = max(accum, diam)
        end do
      end do
      diam = sqrt(diam)
      epsr = LEPSR * diam
      epsi = LEPSI * diam
      epsd = LEPSD * diam
c
c     initial median centre = the centroid
c
      u1 = ONE / dble(n)
      do j = 1, p
        accum = ZERO
        do i = 1, n
          accum = accum + x(i,j)
        end do
        median(j) = accum * u1
      end do
c
c     main loop
c
      do l = 1, ICOUNT
c
c     direction cosines and resultant
c
        corner = ZERO
        do j = 1, p
          c(j) = ZERO
        end do
        ll = l
        do 11 i = 1, n
          d = ZERO
          do j = 1, p
            d = d + (x(i,j) - median(j))**2
          end do
          dd = sqrt(d)
          if (dd .GT. epsd) goto 9
          corner = corner + ONE
          ii = i
          goto 11
    9     d = ONE / dd
          do j = 1, p
            c(j) = c(j) + (x(i,j) - median(j)) * d
          end do
   11   continue
        d = ZERO
        do j = 1, p
          d = d + c(j)**2
        end do
        d = sqrt(d)
        dd = d
c
c     tests for zero resultant or degenerate solution
c
        if (corner .EQ. ZERO) goto 13
        if (d .LE. corner) goto 25
        d = d - corner
   13   if (d .LE. epsr) goto 24
        dd = ONE / dd
        do j = 1, p
          c(j) = c(j) * dd
        end do
c
c     step by bisection to give zero component at lambda
c
        u1 = ZERO
        u2 = diam
        do lc = 1, LCOUNT
          comp = ZERO
          lambda = 0.5 * (u1 + u2)
          slam = lambda * lambda
          do j = 1, p
            z(j) = median(j) + lambda * c(j)
          end do
          do i = 1, n
            delta = ZERO
            d = slam
            do j = 1, p
              xx = x(i,j)
              d = d - (xx - median(j))**2
              delta = delta + (xx - z(j))**2
            end do
            dd = sqrt(delta)
            if (dd .LT. epsd) goto 21
            comp = comp - (d + delta) / dd
          end do
          if (comp .GT. ZERO) goto 18
          u2 = lambda
          goto 19
   18     u1 = lambda
   19     if ((u2 - u1) .LE. epsi) goto 21
        end do
   21   do j = 1, p
          median(j) = median(j) + c(j) * lambda
        end do
      end do
c
   24 iter = ll
      return
c
   25 iter = -ll
      do j = 1, p
        median(j) = x(ii,j)
      end do
      return
c
      END
