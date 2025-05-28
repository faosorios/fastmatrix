c ID: hankel.f, last updated 2025-05-26, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE hankel_mat(x, y, n, c, ldc, info)
      INTEGER          n, ldc, info
      DOUBLE PRECISION c(ldc,*), x(*), y(*)
c
c     contructs the n-by-n Hankel matrix based on the x vector
c
c     parameters:
c     x       (input) DOUBLE PRECISION array, dimension (n)
c             first column of the Hankel matrix
c     y       (input) DOUBLE PRECISION array, dimension (n)
c             last column of the Hankel matrix, x(n) = y(1)
c     n       (input) INTEGER
c             length of vector x. n > 0.
c     c       (input/output) DOUBLE PRECISION array, dimension (ldc, n)
c             a square matrix. on exit if info = 0, is the Hankel matrix.
c     ldc     (input) INTEGER
c             ldc is the leading dimension of the array c. ldc >= max(1,n).
c     info    (output) INTEGER
c             = 0: successful exit
c             < 0: if info = -i, the i-th argument had an illegal value
c
c     .. local scalars ..
      INTEGER          i, j, k
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
c     setup Hankel matrix c
c
      do i = 1, n
        do j = 1, n
          k = i + j - 1
          if (k .GT. n) then 
            c(i,j) = y(i + j - n)
          else 
            c(i,j) = x(i + j - 1)
          end if
        end do
      end do
c
      return
      END
