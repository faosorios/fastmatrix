c ID: dupl_mat.f, last updated 2020-07-31, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE dupl_mat(a, lda, n, col)
      INTEGER lda, n
      INTEGER a(lda,*), col(*)
c
c     dupl_mat contructs the duplication matrix of order 'n' from their compact
c     representation
c
c     parameters:
c     a     - integer array of dimension n**2-by-n*(n+1)/2.
c           on exit a is the duplication matrix of order n.
c     lda   - integer.
c           on entry, lda specifies the leading dimension of a.
c           unchanged on exit.
c     n     - integer.
c           on entry, n specifies the order of duplication matrix a.
c           unchanged on exit.
c     col   - integer array of dimension n**2
c           col contains integers that allow to construct the duplication matrix
c           unchanged on exit.
c
      INTEGER i, nrow
c
      nrow = n * n
c
      do i = 1, nrow
        a(i,col(i)) = 1
      end do
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE dlmm(col, n, a, lda, arow, acol, b, ldb)
      INTEGER          n, lda, arow, acol, ldb
      INTEGER          col(*)
      DOUBLE PRECISION a(lda,*), b(ldb,*)
c
      INTEGER i, j, nrow, ncol
c
      nrow = n * n
      ncol = n * (n + 1) / 2
c
      if (arow .NE. ncol) then
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
