c ID: circulant.f, last updated 2022-02-26, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE circulant_mat(x, n, c, ldc, info)
      INTEGER          n, ldc, info
      DOUBLE PRECISION c(ldc,*), x(*)
c
c     contructs the n-by-n circular matrix based on the v vector
c
c     parameters:
c     x       (input) DOUBLE PRECISION array, dimension (n)
c             first column of the circular matrix
c     n       (input) INTEGER
c             length of vector x. n > 0.
c     c       (input/output) DOUBLE PRECISION array, dimension (ldc, n)
c             a square matrix. on exit if info = 0, is the circulant matrix.
c     ldc     (input) INTEGER
c             ldc is the leading dimension of the array c. ldc >= max(1,n).
c     info    (output) INTEGER
c             = 0: successful exit
c             < 0: if info = -i, the i-th argument had an illegal value
c
c     .. local scalars ..
      INTEGER          j
c
c     test the input parameters
c
      info = 0
      if (n .LT. 0) then
        info = -2
      else if (ldc .LT. max(1, n)) then
        info = -4
      end if
      if (info .NE. 0) return
c
c     quick return if possible
c
      if ((n .EQ. 0) .OR. (ldc .EQ. 0)) return
c
c     setup circular matrix c, column by column
c
      call dcopy(n, x(1), 1, c(1,1), 1)
      do j = 2, n
        call dcopy(n - 1, c(2,j-1), 1, c(1,j), 1)
        c(n,j) = c(1,j-1)
      end do
c
      return
      END
