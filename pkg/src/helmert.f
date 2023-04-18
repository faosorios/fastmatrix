c ID: helmert.f, last updated 2021-05-09, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE helmert_mat(x, ldx, n, info)
      INTEGER          ldx, n, info
      DOUBLE PRECISION x(ldx,*)
c
c     helmert_mat contructs the helmert matrix of order 'n'.
c
c     parameters:
c     x       (output) DOUBLE PRECISION array, dimension (ldx, n)
c             the helmert matrix of order 'n'
c     ldx     (input) INTEGER
c             leading dimension of the array x. ldx >= max(1, n).
c     n       (input) INTEGER
c             order of helmert matrix. n > 0.
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
      if (n .LT. 0) then
        info = -3
      else if (ldx .LT. max(1, n)) then
        info = -2
      end if
      if (info .NE. 0) return
c
c     quick return if possible
c
      if (n .EQ. 0) return
c
c     1st row of Helmert matrix
c
      do j = 1, n
        x(1,j) = 1.d0 / sqrt(dble(n))
      end do
c
c     rows 2 to n:
c
      do i = 2,n
        do j = 1,i-1
          x(i,j) = 1.d0 / sqrt(dble(i * (i - 1)))
        enddo
        x(i,i) = -(i - 1) / sqrt(dble(i * (i - 1)))
      enddo
c
      return
      END
