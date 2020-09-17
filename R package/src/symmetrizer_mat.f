c ID: symmetrizer_mat.f, last updated 2020-09-26, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE symmetrizer_mat(x, ldx, n, row, col, val, len, info)
      INTEGER          ldx, n, len, info
      INTEGER          row(*), col(*)
      DOUBLE PRECISION x(ldx,*), val(*)
c
c     symmetrizer_mat contructs the symmetrizer matrix of order 'n' from
c     their compact representation.
c
c     parameters:
c     x       (output) DOUBLE PRECISION array, dimension (ldx, n * n)
c             the symmetrizer matrix of order 'n'
c     ldx     (input) INTEGER
c             leading dimension of the array x. ldx >= max(1, n * n).
c     n       (input) INTEGER
c             order of symmetrizer matrix. n > 0.
c     row     (input) INTEGER array, dimension (len)
c             provide compact information to build x matrix.
c     col     (input) INTEGER array, dimension (len)
c             provide compact information to build x matrix.
c     val     (input) DOUBLE PRECISION array, dimension (len)
c             provide compact information to build x matrix.
c     info    (output) INTEGER
c             = 0: successful exit
c             < 0: if info = -i, the i-th argument had an illegal value
c
c     .. local scalars ..
      INTEGER k
c
c     test the input parameters
c
      info = 0
      if (n .LT. 0) then
        info = -3
      else if (ldx .LT. max(1, n * n)) then
        info = -2
      else if (len .LT. 0) then
        info = -7
      end if
      if (info .NE. 0) return
c
c     quick return if possible
c
      if ((n .EQ. 0) .OR. (len .EQ. 0)) return
c
c     start the operations
c
      do k = 1, len
        x(row(k),col(k)) = val(k)
      end do
c
      return
      END
