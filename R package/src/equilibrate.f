c ID: equilibrate.f, last updated 2020-09-03, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE equilibrate_cols(a, lda, n, p, scales, cond, job, info)
      INTEGER          lda, n, p, job, info
      DOUBLE PRECISION a(lda,*), scales(*), cond
c
c     equilibrate the columns of a rectangular matrix using 2-norm
c
c     parameters:
c     a       (input/output) DOUBLE PRECISION array, dimension (lda, p)
c             a rectangular matrix. on exit if info = 0, and job = 1, is
c             the equilibrated matrix.
c     lda     (input) INTEGER
c             lda is the leading dimension of the array a. lda >= max(1,n).
c     n       (input) INTEGER
c             the number of rows of the matrix a. n > 0.
c     p       (input) INTEGER
c             the number of columns of the matrix a. p > 0.
c     scales  (output) DOUBLE PRECISION array, dimension (p)
c             if info = 0, scales contains the column equilibration factors.
c     cond    (output) DOUBLE PRECISION
c             if info = 0, cond contains the ratio between the smallest and
c             the largest equilibration factors. If cond >= 0.1, it is not
c             worth scaling by 'scales'
c     job     (input) INTEGER
c             columns must be equilibrated (scaled)?
c     info    (output) INTEGER
c             = 0: successful exit
c             < 0: if info = -i, the i-th argument had an illegal value
c             > 0: if info = i,  the i-th column of a is exactly zero
c
c     .. parameters ..
      DOUBLE PRECISION ONE, ZERO
      PARAMETER       (ONE = 1.0d0, ZERO = 0.0d0)
c     .. BLAS functions ..
      DOUBLE PRECISION dnrm2
c     .. LAPACK functions ..
      DOUBLE PRECISION dlamch
c     .. local scalars ..
      INTEGER          j
      DOUBLE PRECISION condmin, condmax, smallnum, bignum
c
c     test the input parameters
c
      info = 0
      if (n .LT. 0) then
        info = -3
      else if (p .LT. 0) then
        info = -4
      else if (lda .LT. max(1, n)) then
        info = -2
      end if
      if (info .NE. 0) return
c
c     quick return if possible
c
      if ((n .EQ. 0) .OR. (p .EQ. 0)) then
         cond = ONE
         return
      end if
c
c     get machine constants
c
      smallnum = dlamch('S')
      bignum = ONE / smallnum
c
c     compute column scale factors
c
      do j = 1, p
        scales(j) = dnrm2(n, a(1,j), 1)
      end do
c
c     find the maximum and minimum scale factors
c
      condmin = bignum
      condmax = ZERO
      do j = 1, p
        condmin = min(condmin, scales(j))
        condmax = max(condmax, scales(j))
      end do
c
      if (condmin .EQ. ZERO) then
c
c     find the first zero scale factor and return an error code
c
         do j = 1, p
           if (scales(j) .EQ. ZERO) then
             info = j
             return
           end if
         end do
      else
c
c     invert the scale factors and equilibrate columns if requested
c
        do j = 1, p
          scales(j) = ONE / min(max(scales(j), smallnum), bignum)
          if (job .EQ. 1) then
            call dscal(n, scales(j), a(1,j), 1)
          end if
        end do
c
c     compute condition.
c
        cond = max(condmin, smallnum) / min(condmax, bignum)
      end if
c
      return
      END
