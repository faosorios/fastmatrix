c ID: floyd_warshall.f, last updated 2025-09-20, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE floyd_init(costs, ldc, paths, ldp, n, BIG)
      INTEGER          ldc, ldp, n
      DOUBLE PRECISION BIG
      INTEGER          paths(ldp,*)
      DOUBLE PRECISION costs(ldc,*)
c
c     initialization of Floyd-Warshall algorithm
c
c     parameters:
c     costs   (input/output) DOUBLE PRECISION array, dimension (ldc, n)
c             an square matrix containing the edge costs.
c     ldc     (input) INTEGER
c             ldc is the leading dimension of costs array. ldc >= max(1,n).
c     paths   (input/output) INTEGER array, dimension (ldp, n)
c             an square matrix of indices.
c     ldp     (input) INTEGER
c             ldp is the leading dimension of paths array. ldp >= max(1,n).
c     n       (input) INTEGER
c             number of nodes (order of costs and paths matrices), n > 0.
c     BIG     (input) DOUBLE PRECISION
c             a big number for representing large costs.
c
c     .. local scalars ..
      INTEGER          i, j
c     .. parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER       (ZERO = 0.0d0)
c
c     quick return if possible
c
      if (n .LE. 0) return
      if (ldc .LT. max(1,n)) return
      if (ldp .LT. max(1,n)) return
c
c     initialization
c
      do i = 1, n
        do j = 1, n
          if (i .EQ. j) then 
            costs(i,j) = ZERO
            paths(i,j) = 0
          else 
            paths(i,j) = i
          end if
        end do
      end do
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE floyd_warshall(costs, ldc, paths, ldp, n)
      INTEGER          ldc, ldp, n
      INTEGER          paths(ldp,*)
      DOUBLE PRECISION costs(ldc,*)
c
c     Floyd-Warshall algorithm for finding shortest paths in a directed 
c     graph
c
c     parameters:
c     costs   (input/output) DOUBLE PRECISION array, dimension (ldc, n)
c             an square matrix containing the edge costs.
c     ldc     (input) INTEGER
c             ldc is the leading dimension of costs array. ldc >= max(1,n).
c     paths   (input/output) INTEGER array, dimension (ldp, n)
c             an square matrix of indices.
c     ldp     (input) INTEGER
c             ldp is the leading dimension of paths array. ldp >= max(1,n).
c     n       (input) INTEGER
c             number of nodes (order of costs and paths matrices), n > 0.
c
c     .. local scalars ..
      INTEGER          i, j, k
c
c     quick return if possible
c
      if (n .LE. 0) return
      if (ldc .LT. max(1,n)) return
      if (ldp .LT. max(1,n)) return
c
c     Floyd-Warshall iterations
c
      do k = 1, n
        do i = 1, n
          do j = 1, n
            if (costs(i,k) + costs(k,j) .LT. costs(i,j)) then 
              costs(i,j) = costs(i,k) + costs(k,j)
              paths(i,j) = paths(k,j)
            end if 
          end do
        end do
      end do
c
      return
      END
