c ID: array.f, last updated 2020-08-16, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE arraymult(a, lda, r, p, b, ldb, q, s, x, ldx, n,
     *                     y, ldy, info)
      INTEGER          lda, r, p, ldb, q, s, ldx, n, ldy, info
      DOUBLE PRECISION a(lda,*), b(ldb,*), x(ldx,p,*), y(ldy,r,*)
c
c     array multiplication: y = a * x * b, with x a 3D array
c
c     parameters:
c     a       (input) DOUBLE PRECISION array, dimension (lda, p)
c             a rectangular matrix.
c     lda     (input) INTEGER
c             lda is the leading dimension of the array a. lda >= max(1,r).
c     r       (input) INTEGER
c             the number of rows of the matrix a. r > 0.
c     p       (input) INTEGER
c             the number of columns of the matrix a. p > 0.
c     b       (input) DOUBLE PRECISION array, dimension (lda, p)
c             a rectangular matrix.
c     ldb     (input) INTEGER
c             ldb is the leading dimension of the array b. ldb >= max(1,q).
c     q       (input) INTEGER
c             the number of rows of the matrix b. q > 0.
c     s       (input) INTEGER
c             the number of columns of the matrix b. s > 0.
c     x       (input) DOUBLE PRECISION array, dimension (ldx, p, q)
c             an 3D array.
c     ldx     (input) INTEGER
c             ldx is the leading dimension of array x. ldx >= max(1,n).
c     n       (input) INTEGER
c             the number of 'faces' of array x. n > 0.
c     y       (output) DOUBLE PRECISION array, dimension (ldy, r, s)
c             an 3D array.
c     ldy     (input) INTEGER
c             ldy is the leading dimension of array y. ldy >= max(1,n).
c     info    (output) INTEGER
c             = 0: successful exit
c             < 0: if info = -i, the i-th argument had an illegal value
c
c     .. local scalars ..
      INTEGER          i, j, k, l, t
      DOUBLE PRECISION accum, temp
c     .. parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER       (ZERO = 0.0d0)
c
c     test the input parameters
c
      info = 0
      if (r .LT. 0) then
        info = -3
      else if (p .LT. 0) then
        info = -4
      else if (q .LT. 0) then
        info = -7
      else if (s .LT. 0) then
        info = -8
      else if (n .LT. 0) then
        info = -11
      else if (lda .LT. max(1, r)) then
        info = -2
      else if (ldb .LT. max(1, q)) then
        info = -6
      else if (ldx .LT. max(1, n)) then
        info = -10
      else if (ldy .LT. max(1, n)) then
        info = -13
      end if
      if (info .NE. 0) return
c
c     quick return if possible.
c
      if ((r .EQ. 0) .OR. (p .EQ. 0) .OR. (q .EQ. 0) .OR.
     *    (s .EQ. 0) .OR. (n .EQ. 0)) return
c
c     start the operations.
c
      do t = 1, n
        do k = 1, r
          do l = 1, s
            accum = ZERO
            do i = 1, p
              do j = 1, q
                temp  = a(k,i) * x(t,i,j) * b(j,l)
                accum = accum + temp
              end do
            end do
            y(t,k,l) = accum
          end do
        end do
      end do
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE bracketprod(a, lda, m, n, x, ldx, p, q, y, ldy, info)
      INTEGER          lda, m, n, ldx, p, q, ldy, info
      DOUBLE PRECISION a(lda,*), x(ldx,p,*), y(ldy,p,*)
c
c     bracket product: y = [a][x], with x a 3D array
c
c     parameters:
c     a       (input) DOUBLE PRECISION array, dimension (lda, n)
c             a rectangular matrix.
c     lda     (input) INTEGER
c             lda is the leading dimension of the array a. lda >= max(1,m).
c     m       (input) INTEGER
c             the number of rows of the matrix a. m > 0.
c     n       (input) INTEGER
c             the number of columns of the matrix a. n > 0.
c     x       (input) DOUBLE PRECISION array, dimension (ldx, p, q)
c             an 3D array.
c     ldx     (input) INTEGER
c             ldx is the leading dimension of array x. ldx >= max(1,n).
c     p       (input) INTEGER
c             the number of 'rows' of array x. p > 0.
c     q       (input) INTEGER
c             the number of 'columns' of array x. q > 0.
c     y       (output) DOUBLE PRECISION array, dimension (ldy, p, q)
c             an 3D array.
c     ldy     (input) INTEGER
c             ldy is the leading dimension of array y. ldy >= max(1,m).
c     info    (output) INTEGER
c             = 0: successful exit
c             < 0: if info = -i, the i-th argument had an illegal value
c
c     .. local scalars ..
      INTEGER          i, j, k, t
      DOUBLE PRECISION accum
c     .. parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER       (ZERO = 0.0d0)
c
c     test the input parameters
c
      info = 0
      if (m .LT. 0) then
        info = -3
      else if (n .LT. 0) then
        info = -4
      else if (p .LT. 0) then
        info = -7
      else if (q .LT. 0) then
        info = -8
      else if (lda .LT. max(1, m)) then
        info = -2
      else if (ldx .LT. max(1, n)) then
        info = -6
      else if (ldy .LT. max(1, m)) then
        info = -10
      end if
      if (info .NE. 0) return
c
c     quick return if possible.
c
      if ((m .EQ. 0) .OR. (n .EQ. 0) .OR. (p .EQ. 0) .OR.
     *    (q .EQ. 0)) return
c
c     start the operations.
c
      do t = 1, m
        do i = 1, p
          do j = 1, q
            accum = ZERO
            do k = 1, n
              accum = accum + a(t,k) * x(k,i,j)
            end do
            y(t,i,j) = accum
          end do
        end do
      end do
c
      return
      END
