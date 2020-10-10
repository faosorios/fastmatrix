c ID: cpy_triangle.f, last updated 2020-09-27, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE cpy_upper(a, lda, nrow, ncol, b, ldb)
      INTEGER          lda, nrow, ncol, ldb
      DOUBLE PRECISION a(lda,*), b(ldb,*)
c
c     copy upper triangular part of a into b
c
      INTEGER i, j, p
c
      p = min(nrow, ncol)
      do j = 1, p
        do i = 1, j
          b(i,j) = a(i,j)
        end do
      end do
c
      return
      END
