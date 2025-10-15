c ID: schur_decomp.f, last updated 2025-10-13, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE schur_decomp(a, lda, n, task, wr, wi, v, ldv, work, 
     *                        lwork, dummy, info)
      INTEGER          lda, n, task, ldv, lwork, info
      LOGICAL          dummy(*)
      DOUBLE PRECISION a(lda,*), v(ldv,*), wr(*), wi(*), work(*)
c
c     Computes the Schur decomposition of a real n-by-n nonsymmetric 
c     matrix 'a'. The factorization has the form a = q %*% v %*% t(q), 
c     where 'v' is an upper triangular matrix and 'q' is orthogonal.
c
c     This function is just a simplified wrapper of the subroutine 
c     'DGEES' from LAPACK (next a verbatim copy of )
c
c     parameters:
c     a       (input/output) DOUBLE PRECISION array, dimension (lda, n)
c             on entry, the n-by-n matrix a. On exit, a has been overwritten
c             by its real Schur form upper triangular matrix 't'
c     lda     (input) INTEGER
c             leading dimension of the array a. lda >= max(1,n).
c     n       (input) INTEGER
c             the order of the matrix a. n >= 0.
c     task   (input) INTEGER
c             if task = 0, Schur vectors are not computed
c             if task = 1, Schur vectors are computed
c     wr      (output) DOUBLE PRECISION array, dimension (n)
c             real part of the computed eigenvalues
c     wi      (output) DOUBLE PRECISION array, dimension (n)
c             imaginary part of the computed eigenvalues
c     v       (output) DOUBLE PRECISION array, dimension (ldv, n)
c             if job = 'V', v contains the orthogonal matrix q of 
c             Schur vectors.
c             if job = 'N', v is not referenced
c     ldv     (input) INTEGER
c             leading dimension of the array v. 
c             if job = 'V', ldv >= max(1,n).
c             if job = 'N', ldv >= 1
c     work    (output) DOUBLE PRECISION array, dimension (max(1,lwork))
c             in exit, if info = 0, work(1) contains the optimal lwork.
c     lwork   (input) INTEGER
c             the dimension of the array work. lwork >= max(1,3*n).
c             for good performance, lwork must generally be larger.
c     dummy   (dummy) LOGICAL array, dimension(n)
c             this argument is not referenced
c     info    (output) INTEGER
c             = 0: successful exit
c             otherwise see documentation of DGEES
c
c     .. local scalars ..
      CHARACTER        job, sort
      INTEGER          sdim
c
c     initialization (eigenvalues are not ordered)
c
      sort = 'N'
      sdim = 0
      if (task .EQ. 1) then
        job = 'V'
      else 
        job = 'N'
      end if
c
c     calling DGEES routine
c
      call dgees(job, sort, dummy_select, n, a, lda, sdim, wr, wi,
     *           v, ldv, work, lwork, dummy, info)
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      LOGICAL FUNCTION dummy_select(wr, wi)
      DOUBLE PRECISION wr, wi
c
c     dummy SELECT function called by DGEES
c
      dummy_select = .FALSE.
c
      return
      END
