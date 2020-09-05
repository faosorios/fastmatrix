c ID: inner_frobenius.f, last updated 2020-09-03, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE inner_frobenius(a, lda, b, ldb, n, p, value)
      INTEGER          lda, ldb, n, p
      DOUBLE PRECISION a(lda,*), b(ldb,*), value
c
c     computes the Frobenius inner product between matrices
c
c     parameters:
c     a       (input) DOUBLE PRECISION array, dimension (lda, p)
c             a rectangular matrix.
c     lda     (input) INTEGER
c             lda is the leading dimension of the array a. lda >= max(1,n).
c     b       (input) DOUBLE PRECISION array, dimension (ldb, p)
c             a rectangular matrix.
c     ldb     (input) INTEGER
c             ldb is the leading dimension of the array b. ldb >= max(1,n).
c     n       (input) INTEGER
c             the number of rows of matrices a and b. n > 0.
c     p       (input) INTEGER
c             the number of columns of matrices a and b. p > 0.
c     value   (output) DOUBLE PRECISION
c             Euclidean norm between matrices a and b.
c
c     .. BLAS functions ..
      DOUBLE PRECISION ddot
c     .. local scalars ..
      INTEGER          j
      DOUBLE PRECISION temp
c
      value = 0.0d0
c
c     quick return if possible
c
      if (n .LE. 0) return
      if (p .LE. 0) return
      if (lda .LT. max(1,n)) return
      if (ldb .LT. max(1,n)) return
c
c     compute Frobenius inner product between a and b, column by
c     column, updating partial result
c
      do j = 1, p
        temp  = ddot(n, a(1,j), 1, b(1,j), 1)
        value = value + temp
      end do
c
      return
      END
