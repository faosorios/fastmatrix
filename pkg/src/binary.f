c ID: binary.f, last updated 2025-10-23, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE int2logical(y, x, n)
      INTEGER    n, x(*)
      LOGICAL    y(*)
c
c     convert integers to logical type
c
c     parameters:
c     y       (output) LOGICAL array, dimension (n)
c             a logical vector
c     x       (input) INTEGER array, dimension (n)
c             a vector of integers
c     n       (input) INTEGER
c             dimension of the vector, n > 0.
c
c     .. local scalars ..
      INTEGER i, m
c
c     quick return if possible
c
      if (n .LE. 0) return
c
c     start the operations
c
      m = mod(n, 7)
      do i = 1, m, 7
        y(i)     = (x(i)     .NE. 0)
        y(i + 1) = (x(i + 1) .NE. 0)
        y(i + 2) = (x(i + 2) .NE. 0)
        y(i + 3) = (x(i + 3) .NE. 0)
        y(i + 4) = (x(i + 4) .NE. 0)
        y(i + 5) = (x(i + 5) .NE. 0)
        y(i + 6) = (x(i + 6) .NE. 0)
      end do 
      do i = m + 1, n
        y(i) = (x(i) .NE. 0)
      end do
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     WARNING:
c
c     The following recommendation should be taken into consideration 
c     when invoking the routines defined below
c
c     To guarantee enough space for any binary vector 'x', its dimension 
c     'n' would need to be set to 32 (the routines do not verify this).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE add_binary(z, x, y, n)
      INTEGER    n, x(*), y(*), z(*)
c
c     add two (signed) binary vectors
c
c     parameters:
c     z       (output) INTEGER array, dimension (n)
c             a vector of binary digits
c     y       (input) INTEGER array, dimension (n)
c             a vector of binary digits
c     x       (input) INTEGER array, dimension (n)
c             a vector of binary digits
c     n       (input) INTEGER
c             dimension of the vector, n > 0.
c
c     .. parameters ..
      INTEGER    BASE, ONE
      PARAMETER (BASE = 2, ONE = 1)
c     .. local scalars ..
      LOGICAL    overflow
      INTEGER    i
c
c     quick return if possible
c
      if (n .LE. 0) return
c
c     start the operations
c
      overflow = .FALSE.
      do i = 1, n
        z(i) = x(i) + y(i)
      end do

      do i = 1, n
        if (BASE .LE. z(i)) then
          z(i) = z(i) - BASE
          if (ONE .LT. i) then 
            z(i-1) = z(i-1) + ONE
          else 
            overflow = .TRUE.
          end if
        end if
      end do
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE check_binary(x, n, info)
      INTEGER    info, n, x(*)
c
c     check entries of a binary vector
c     (the only check made is that the entries are all 0 or 1)
c
c     parameters:
c     x       (input) INTEGER array, dimension (n)
c             a vector of binary digits
c     n       (input) INTEGER
c             dimension of the vector, n > 0.
c     info    (output) INTEGER
c             = 0: successful exit
c             > 0: if info = i, the i-th element had an illegal value 
c                  routine returns when the first illegal element is 
c                  found
c
c     .. parameters ..
      INTEGER    BASE, ZERO
      PARAMETER (BASE = 2, ZERO = 0)
c     .. local scalars ..
      INTEGER    i
c
c     quick return if possible
c
      if (n .LE. 0) return
c
c     start the operations
c
      info = 0
      do i = 1, n
        if ((x(i) .LT. ZERO) .OR. (BASE .LE. x(i))) then
          info = i
          return
        end if
      end do
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE two_complement(y, x, n)
      INTEGER    n, x(*), y(*)
c
c     computes the two's complement of a binary vector
c
c     parameters:
c     y       (output) INTEGER array, dimension (n)
c             the two's complemented vector
c     x       (input) INTEGER array, dimension (n)
c             the vector to be complemented
c     n       (input) INTEGER
c             dimension of the vector, n > 0.
c
c     .. parameters ..
      INTEGER    BASE, ONE, ZERO
      PARAMETER (BASE = 2, ONE = 1, ZERO = 0)
c     .. local scalars ..
      INTEGER    i
c     .. local arrays ..
      INTEGER    u(100), v(100)
c
c     quick return if possible
c
      if (n .LE. 0) return
c
c     setting vectors
c
      do i = 1, n
        u(i) = (BASE - ONE) - x(i)
      end do

      do i = 1, n - 1
        v(i) = ZERO
      end do
      v(n) = ONE

      call add_binary(y, u, v, n)
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE int2binary(y, x, n)
      INTEGER    n, x, y(*)
c
c     returns the signed binary representation of an integer
c
c     parameters:
c     y       (output) INTEGER array, dimension (n)
c             the signed binary representation
c     x       (input) INTEGER
c             an integer to be represented
c     n       (input) INTEGER
c             dimension of the vector, n > 0.
c
c     .. parameters ..
      INTEGER    BASE, ZERO
      PARAMETER (BASE = 2, ZERO = 0)
c     .. local scalars ..
      INTEGER    i, p
c
c     quick return if possible
c
      if (n .LE. 0) return
c
c     start the operations
c
      p = ABS(x)
      do i = n, 2, -1
        y(i) = mod(p, BASE)
        p = p / BASE
      end do

      y(1) = ZERO
      if (p .LT. ZERO) then
        call two_complement(y, y, n)
      end if
c
      return
      END
