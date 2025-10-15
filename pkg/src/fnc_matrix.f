c ID: fnc_matrix.f, last updated 2025-10-12, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE fnc_parlett(t, ldt, n, f, ldf)
      INTEGER          ldt, ldf, n
      DOUBLE PRECISION t(ldt,*), f(ldf,*)
c
c     Parlett method for the evaluation of a function of a triangular 
c     matrix
c
c     parameters:
c     t       (input) DOUBLE PRECISION array, dimension (ldt, n)
c             a n-by-n triangular matrix.
c     ldt     (input) INTEGER
c             ldt is the leading dimension of the array t. ldt >= max(1,n).
c     n       (input) INTEGER
c             the number of rows of the matrix t. n > 0.
c     f       (output) DOUBLE PRECISION array, dimension (ldf, n)
c             an square matrix.
c     ldf     (input) INTEGER
c             ldt is the leading dimension of the array f. ldf >= max(1,n).
c
c     .. local scalars ..
      INTEGER          i, j, k, p
      DOUBLE PRECISION s
c
c     quick return if possible
c
      if (n .LE. 0) return
      if (ldt .LT. max(1,n)) return
      if (ldf .LT. max(1,n)) return
c
c     scalar Parlett recurrence
c
      do p = 1, n - 1
        do i = 1, n - p
          j = i + p
          s = t(i,j) * (f(j,j) - f(i,i))
          do k = i + 1, j - 1
            s = s + t(i,k) * f(k,j) - f(i,k) * t(k,j)
          end do
          f(i,j) = s / (t(j,j) - t(i,i))
        end do
      end do
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE sqrt_parlett(t, ldt, n, u, ldu, info)
      INTEGER          ldt, ldu, n, info
      DOUBLE PRECISION t(ldt,*), u(ldu,*)
c
c     use a Parlett recurrence to obtain a primary square root of a 
c     triangular matrix
c
c     parameters:
c     t       (input) DOUBLE PRECISION array, dimension (ldt, n)
c             a n-by-n triangular matrix.
c     ldt     (input) INTEGER
c             ldt is the leading dimension of the array t. ldt >= max(1,n).
c     n       (input) INTEGER
c             the number of rows of the matrix t. n > 0.
c     u       (output) DOUBLE PRECISION array, dimension (ldu, n)
c             an square matrix.
c     ldu     (input) INTEGER
c             ldt is the leading dimension of the array u. ldu >= max(1,n).
c     info    (output) INTEGER
c             = 0: successful exit
c             < 0: if info = -i, the i-th argument had an illegal value
c             > 0: if info =  i, the element of triangular matrix have 
c                  an illegal value or it is zero.
c
c     .. local scalars ..
      INTEGER          i, j, k
      DOUBLE PRECISION accum
c     .. parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER       (ZERO = 0.0d0)
c
c     test the input parameters
c
      info = 0
      if (n .LT. 0) then
        info = -3
      else if (ldt .LT. max(1, n)) then
        info = -2
      else if (ldu .LT. max(1, n)) then
        info = -5
      end if
      if (info .NE. 0) return
c
c     quick return if possible
c
      if (n .EQ. 0) return
c
c     checking diagonal elements
c
      do i = 1, n
        if (t(i,i) .LE. ZERO) then
          info = i
          return
        else 
          u(i,i) = dsqrt(t(i,i))
        end if
      end do
c
c     start the recurrence
c
      do j = 2, n
        do i = j - 1, 1, -1
          accum = ZERO
          do k = i + 1, j - 1
            accum = accum + u(i,k) * u(k,j)
          end do
          u(i,j) = (t(i,j) - accum) / (u(i,i) + u(j,j))
        end do
      end do
c
      return
      END

