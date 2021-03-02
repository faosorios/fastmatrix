c ID: commutation.f, last updated 2020-08-26, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE comm_rows(m, n, row)
      INTEGER m, n
      INTEGER row(*)
c
c     compact information to build a commutation matrix
c
      INTEGER i, mn
      REAL    ratio
c
c     quick return if possible
c
      if ((m .LE. 0) .OR. (n .LE. 0)) return
c
c     start the operations
c
      mn = m * n
      do i = 1, mn
        ratio  = (i - 1) / m
        row(i) = modulo(i - 1, m) * n + floor(ratio) + 1
      end do
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE commutation_mat(x, ldx, m, n, row, info)
      INTEGER ldx, m, n, info
      INTEGER x(ldx,*), row(*)
c
c     commutation_mat contructs the commutation matrix of order 'mn' from
c     their compact representation. Let a an m-by-n matrix, commutation
c     matrix of order 'mn' transform vec(a) into vec(a')
c
c     parameters:
c     x       (output) INTEGER array, dimension (ldx, m * n)
c             the commutation matrix of order 'mn'
c     ldx     (input) INTEGER
c             leading dimension of the array x. ldx >= max(1, m * n).
c     m       (input) INTEGER
c             the number of rows of the matrix a. m > 0.
c     n       (input) INTEGER
c             the number of columns of the matrix a. n > 0.
c     row     (input) INTEGER array, dimension (m * n)
c             provide compact information to build x matrix.
c     info    (output) INTEGER
c             = 0: successful exit
c             < 0: if info = -i, the i-th argument had an illegal value
c
c     .. local scalars ..
      INTEGER j, ncol
c
c     test the input parameters
c
      info = 0
      if (m .LT. 0) then
        info = -3
      else if (n .LT. 0) then
        info = -4
      else if (ldx .LT. max(1, m * n)) then
        info = -2
      end if
      if (info .NE. 0) return
c
c     quick return if possible
c
      if ((m .EQ. 0) .OR. (n .EQ. 0)) return
c
c     start the operations
c
      ncol = m * n
c
      do j = 1, ncol
        x(row(j),j) = 1
      end do
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE comm_left_mult(col, m, n, a, lda, arow, acol, b, ldb,
     *                          info)
      INTEGER          m, n, lda, arow, acol, ldb, info
      INTEGER          col(*)
      DOUBLE PRECISION a(lda,*), b(ldb,*)
c
c     .. local scalars ..
      INTEGER          i, j, nrow, ncol
c
c     test the input parameters
c
      info = 0
      if (m .LT. 0) then
        info = -2
      else if (n .LT. 0) then
        info = -3
      else if (lda .LT. max(1, m * n)) then
        info = -5
      else if (arow .LT. 0) then
        info = -6
      else if (acol .LT. 0) then
        info = -7
      else if (ldb .LT. max(1, m * n)) then
        info = -9
      end if
      if (info .NE. 0) return
c
c     quick return if possible.
c
      if ((m .EQ. 0) .OR. (n .EQ. 0) .OR. (arow .EQ. 0) .OR.
     *    (acol .EQ. 0)) return
c
c     start the operations
c
      nrow = m * n
      ncol = m * n
c
c     check if matrices are compatible
c
      if (ncol .NE. arow) then
        info = 1
        return
      end if
c
      do j = 1, acol
        do i = 1, nrow
          b(i,j) = a(col(i),j)
        end do
      end do
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE comm_right_mult(col, m, n, a, lda, arow, acol, b, ldb,
     *                           info)
      INTEGER          m, n, lda, arow, acol, ldb, info
      INTEGER          col(*)
      DOUBLE PRECISION a(lda,*), b(ldb,*)
c
c     .. local scalars ..
      INTEGER          i, j, nrow, ncol
c
c     test the input parameters
c
      info = 0
      if (m .LT. 0) then
        info = -2
      else if (n .LT. 0) then
        info = -3
      else if (arow .LT. 0) then
        info = -6
      else if (acol .LT. 0) then
        info = -7
      else if (lda .LT. max(1, arow)) then
        info = -5
      else if (ldb .LT. max(1, arow)) then
        info = -9
      end if
      if (info .NE. 0) return
c
c     quick return if possible.
c
      if ((m .EQ. 0) .OR. (n .EQ. 0) .OR. (arow .EQ. 0) .OR.
     *    (acol .EQ. 0)) return
c
c     start the operations
c
      nrow = m * n
      ncol = m * n
c
c     check if matrices are compatible
c
      if (acol .NE. nrow) then
        info = 1
        return
      end if
c
      do j = 1, ncol
        do i = 1, arow
          b(i,j) = a(i, col(j))
        end do
      end do
c
      return
      END
