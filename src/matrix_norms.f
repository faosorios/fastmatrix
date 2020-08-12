c ID: matrix_norms.f, last updated 2020-08-11, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE maxcol_norm(a, lda, n, p, value)
      INTEGER          lda, n, p
      DOUBLE PRECISION a(lda,*), value
c
c     computes the 1-norm of matrix 'a'
c
c     .. parameters ..
      INTEGER          ONE
      PARAMETER       (ONE = 1)
c     .. BLAS functions ..
      DOUBLE PRECISION dasum
c     .. local scalars ..
      INTEGER          j
      DOUBLE PRECISION temp
c
      value = 0.0d0
c
      if (n .LE. 0) return
      if (p .LE. 0) return
      if (lda .LT. max(1,n)) return
c
      do j = 1, p
        temp  = dasum(n, a(1,j), ONE)
        value = max(temp, value)
      end do
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE maxrow_norm(a, lda, n, p, value)
      INTEGER          lda, n, p
      DOUBLE PRECISION a(lda,*), value
c
c     computes the infinite-norm of matrix 'a'
c
c     .. BLAS functions ..
      DOUBLE PRECISION dasum
c     .. local scalars ..
      INTEGER          i
      DOUBLE PRECISION temp
c
      value = 0.0d0
c
      if (n .LE. 0) return
      if (p .LE. 0) return
      if (lda .LT. max(1,n)) return
c
      do i = 1, n
        temp  = dasum(p, a(i,1), lda)
        value = max(temp, value)
      end do
c
      return
      END
