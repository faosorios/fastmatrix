c ID: ldl_decomp.f, last updated 2020-09-18, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE ldl_dcmp(a, lda, n, d, info)
      INTEGER          lda, n, info
      DOUBLE PRECISION a(lda,*), d(*)
c
c     Computes the LDL decomposition of a real symmetric matrix 'a'. The
c     factorization has the form a = l %*% d %*% t(l), where 'd' is a diagonal
c     matrix and 'l' is unitary lower triangular.
c
c     parameters:
c     a       (input/output) DOUBLE PRECISION array, dimension (lda, n)
c             on entry, the symmetric matrix a. the leading n-by-n upper
c             triangular part of a contains the upper triangular part of
c             the matrix a. on exit, if info = 0, the factor 'l' from the
c             LDL factorization is stored on the strictly lower triangular
c             part of a.
c     lda     (input) INTEGER
c             lda is the leading dimension of the array a. lda >= max(1,n).
c     n       (input) INTEGER
c             the order of the matrix a. n >= 0
c     d       (output) DOUBLE PRECISION array, dimension (n)
c             contains the diagonal elements of matrix 'd'.
c     info    (output) INTEGER
c             = 0: successful exit
c             < 0: if info = -i, the i-th argument had an illegal value
c
c     .. parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER       (ZERO = 0.0d0)
c     .. local scalars ..
      INTEGER i, j, k
      DOUBLE PRECISION accum
c
c     test the input parameters
c
      info = 0
      if (n .LT. 0) then
        info = -3
      else if (lda .LT. max(1, n)) then
        info = -2
      end if
      if (info .NE. 0) return
c
c     quick return if possible
c
      if (n .EQ. 0) return
c
c     start the operations
c
      do j = 1, n
        accum = ZERO
        do k = 1, j - 1
          accum = accum + d(k) * a(j,k)**2
        end do
        d(j) = a(j,j) - accum
        do i = j + 1, n
          accum = ZERO
          do k = 1, j - 1
            accum = accum + d(k) * a(i,k) * a(j,k)
          end do
          a(i,j) = (a(j,i) - accum) / d(j)
        end do
      end do
c
      return
      END
