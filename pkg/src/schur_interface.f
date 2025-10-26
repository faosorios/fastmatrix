c ID: schur_interface.f, last updated 2025-10-13, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE schur_wrapper(a, lda, n, task, re, im, v, ldv, work, 
     *                         lwork, bwork, info)
      INTEGER          lda, n, task, ldv, lwork, bwork, info
      DOUBLE PRECISION a(lda,*), v(ldv,*), re(*), im(*), work(*)
c
c     interface to DGEES routine from LAPACK, using a simplified call
c     with specific values for certain parameters
c
c     .. local scalars ..
      CHARACTER        job, sort
      LOGICAL          dummy
      INTEGER          sdim
c
c     initialization
c
      sort  = 'N'
      sdim  = 0
      dummy = (bwork .NE. 0)
      if (task .EQ. 1) then
        job = 'V'
      else 
        job = 'N'
      end if
c
c     calling DGEES routine (eigenvalues are not ordered)
c
      call dgees(job, sort, dummy_select, n, a, lda, sdim, re, im,
     *           v, ldv, work, lwork, dummy, info)
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      LOGICAL FUNCTION dummy_select()
c
c     SELECT function for unordered result (called by DGEES)
c
      dummy_select = .FALSE.
c
      return
      END
