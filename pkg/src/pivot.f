c ID: pivot.f, last updated 2020-09-14, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE pivot_mat(a, lda, n, pivot)
      INTEGER          lda, n
      INTEGER          pivot(*)
      DOUBLE PRECISION a(lda,*)
c
c     apply column interchanges provided by 'pivot' information
c
      INTEGER j, jp
c
c     quick return if possible
c
      if ((lda .LT. max(1, n)) .OR. (n .LE. 0)) return
c
c     start the operations
c
      do j = n - 1, 1, -1
        jp = pivot(j)
        if (jp .NE. j) then
          call dswap(n, a(1,j), 1, a(1, jp), 1)
        end if
      end do
c
      return
      END
