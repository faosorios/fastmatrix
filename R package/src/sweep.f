c ID: sweep.f, last updated 2020-08-21, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE sweepop(a, lda, p, k, reverse, info)
      INTEGER          lda, p, k, reverse, info
      DOUBLE PRECISION a(lda,*)
c
c     Gauss-Jordan sweep operator for symmetric matrices
c
c     parameters:
c     a       (input/output) DOUBLE PRECISION array, dimension (lda, p)
c             a symmetric matrix. on exit if info = 0, a is the sweeped matrix
c     lda     (input) INTEGER
c             lda is the leading dimension of the array a. lda >= max(1,p).
c     p       (input) INTEGER
c             p is the order of the matrix. k > 0.
c     k       (input) INTEGER
c             k element of the diagonal which will be sweeped, k > 0.
c     reverse (input) INTEGER
c             reverse sweep operator must be applied?
c     info    (output) INTEGER
c             = 0: successful exit
c             < 0: if info = -i, the i-th argument had an illegal value
c             > 0: if info = k,  the element a(k,k) is zero
c
c     .. local scalars ..
      INTEGER          i
      DOUBLE PRECISION div, piv
c
c     test the input parameters
c
      info = 0
      if (p .LE. 0) then
        info = -3
      else if (lda .LT. max(1, p)) then
        info = -2
      else if ((k .LE. 0) .OR. (k .GT. p)) then
        info = -4
      end if
      if (info .NE. 0) return
c
c     quick return if possible
c
      if (a(k,k) .EQ. 0.d0) then
        info = k
        return
      end if
c
      div = a(k,k)
      call dscal(p, 1.0d0 / div, a(k,1), lda)
c
      do i = 1, p
        if (i .NE. k) then
          piv = a(i,k)
          call daxpy(p, -piv, a(k,1), lda, a(i,1), lda)
          if (reverse .EQ. 1) then
            a(i,k) = piv / div
          else
            a(i,k) = -piv / div
          end if
        end if
      end do
      a(k,k) = 1.0d0 / div
c
      return
      END
