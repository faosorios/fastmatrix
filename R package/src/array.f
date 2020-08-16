c ID: array.f, last updated 2020-08-16, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE arraymult(a, lda, r, p, b, ldb, q, s, x, ldx, n,
     *                     y, ldy, info)
      INTEGER          lda, r, p, ldb, q, s, ldx, n, ldy, info
      DOUBLE PRECISION a(lda,*), b(ldb,*), x(ldx,p,*), y(ldy,r,*)
c
c     .. local scalars ..
      INTEGER          i, j, k, l, t
      DOUBLE PRECISION accum, temp
c     .. parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER       (ZERO = 0.0d0)
c     ..
c     .. executable statements ..
c
c     test the input parameters
c
      info = 0
      if (r .LT. 0) then
        info = 3
      else if (p .LT. 0) then
        info = 4
      else if (q .LT. 0) then
        info = 7
      else if (s .LT. 0) then
        info = 8
      else if (lda .LT. max(1, r)) then
        info = 2
      else if (ldb .LT. max(1, q)) then
        info = 6
      else if (ldx .LT. max(1, n)) then
        info = 10
      else if (ldy .LT. max(1, n)) then
        info = 13
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
c
c     end of arraymult
c
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE bracketprod(a, lda, m, n, x, ldx, p, q, y, ldy, info)
      INTEGER          lda, m, n, ldx, p, q, ldy, info
      DOUBLE PRECISION a(lda,*), x(ldx,p,*), y(ldy,p,*)
c
c     .. local scalars ..
      INTEGER          i, j, k, t
      DOUBLE PRECISION accum
c     .. parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER       (ZERO = 0.0d0)
c     ..
c     .. executable statements ..
c
c     test the input parameters
c
      info = 0
      if (m .LT. 0) then
        info = 3
      else if (n .LT. 0) then
        info = 4
      else if (p .LT. 0) then
        info = 7
      else if (q .LT. 0) then
        info = 8
      else if (lda .LT. max(1, m)) then
        info = 2
      else if (ldx .LT. max(1, n)) then
        info = 6
      else if (ldy .LT. max(1, m)) then
        info = 10
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
c
c     end of bracketprod
c
      END
