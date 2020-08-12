c ID: norms_extra.f, last updated 2020-08-11, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION FUNCTION minkowski(n, x, inc, p)
      INTEGER          n, inc
      DOUBLE PRECISION x(*), p
c
c     computes the p-norm of x using an unrolled loop
c
      DOUBLE PRECISION accum
      INTEGER i, m, mp1, ninc
c
      minkowski = 0.0d0
      accum = 0.0d0
      if ((n .LE. 0) .OR. (inc .LE. 0)) return
      if (inc .EQ. 1) goto 20
c
c     code for increment not equal to 1
c
      ninc = n * inc
      do i = 1, ninc, inc
        accum = accum + dabs(x(i))**p
      end do
      minkowski = accum
      return
c
c     code for increment equal to 1
c
   20 m = mod(n, 8)
      if (m .EQ. 0) goto 21
      do i = 1, m
        accum = accum + dabs(x(i))**p
      end do
      if (n .LT. 8) goto 22
   21 mp1 = m + 1
      do i = mp1, n, 8
        accum = accum + dabs(x(i))**p + dabs(x(i + 1))**p
     *  + dabs(x(i + 2))**p + dabs(x(i + 3))**p + dabs(x(i + 4))**p
     *  + dabs(x(i + 5))**p + dabs(x(i + 6))**p + dabs(x(i + 7))**p
      end do
   22 minkowski = accum
c
      return
      END
