c ID: frank.f, last updated 2022-08-01, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE frank_mat(f, ldf, n, info)
      INTEGER          ldf, n, info
      DOUBLE PRECISION f(ldf,*)
c
c     constructs the n-by-n Frank matrix
c
c     parameters:
c     f       (output) DOUBLE PRECISION array, dimension (ldf, n)
c             a square matrix. on exit if info = 0, f is the Frank matrix
c     ldf     (input) INTEGER
c             ldf is the leading dimension of the array f. ldf >= max(1,n).
c     n       (input) INTEGER
c             n is the order of the matrix. n > 0.
c     info    (output) INTEGER
c             = 0: successful exit
c             < 0: if info = -i, the i-th argument had an illegal value
c
c     .. local scalars ..
      INTEGER i, j
c
c     test the input parameters
c
      info = 0
      if (n .LE. 0) then
        info = -2
      else if (ldf .LT. max(1, n)) then
        info = -3
      end if
      if (info .NE. 0) return
c
c     start the operations
c
      do i = 1, n
        do j = 1, n
          if (i .LE. j) then
            f(i,j) = dble(n - j + 1)
          elseif (i .EQ. j + 1) then
            f(i,j) = dble(n - j)
          else
            f(i,j) = 0.d0
          end if
        end do
      end do
c
      return
      END
