c ID: ccc.f, last updated 2025-06-23, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE rhoc_ustat(x, y, n, p1, p2, p3)
      INTEGER          n
      DOUBLE PRECISION x(*), y(*), p1(*), p2(*), p3(*)
c
c     contructs kernels of the U-statistics
c
c     parameters:
c     x       (input) DOUBLE PRECISION array, dimension (n)
c             observations of 1st measurement instrument
c     y       (input) DOUBLE PRECISION array, dimension (n)
c             observations of 2nd measurement instrument
c     n       (input) INTEGER
c             length of vectors x and y. n > 0.
c     p1      (output) DOUBLE PRECISION array, dimension (n)
c             kernel of the 1st component of U-statistic
c     p2      (output) DOUBLE PRECISION array, dimension (n)
c             kernel of the 2nd component of U-statistic
c     p3      (output) DOUBLE PRECISION array, dimension (n)
c             kernel of the 3rd component of U-statistic
c
c     .. local scalars ..
      INTEGER          i, j
      DOUBLE PRECISION acc1, acc2, acc3, up1, up2, up3
      DOUBLE PRECISION d1, d2, d3, d4, s1, s2, s3, s4
c
c     quick return if possible
c
      if (n .LE. 0) return
c
c     computing kernels of U-statistics
c
      do i = 1, n
        acc1 = 0.d0
        acc2 = 0.d0
        acc3 = 0.d0
        do j = 1, n
          if (i .NE. j) then
            d1 = (x(i) - y(i))**2
            d2 = (x(j) - y(j))**2
            s1 = (x(i) + y(i))**2
            s2 = (x(j) + y(j))**2
            d3 = (x(i) - y(j))**2
            d4 = (x(j) - y(i))**2
            s3 = (x(i) + y(j))**2
            s4 = (x(j) + y(i))**2
            up1 = (d1 + d2) / 2 - (s1 + s2) / 2
            up2 = (x(i)**2 + x(j)**2) + (y(i)**2 + y(j)**2)
            up3 = (d3 - s3) / 2 + (d4 - s4) / 2
            acc1 = acc1 + up1
            acc2 = acc2 + up2
            acc3 = acc3 + up3
          end if
        end do
        p1(i) = acc1 / (n - 1)
        p2(i) = acc2 / (n - 1)
        p3(i) = acc3 / (n - 1)
      end do
c
      return
      END

