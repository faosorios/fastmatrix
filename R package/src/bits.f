c ID: bits.f, last updated 2020-08-13, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE bitreversal(c, t)
      INTEGER c(0:*), t
c
c     Elster's linear bit reversal algorithm
c
      INTEGER j, k, n, q
c
      n = 2**t
      k = 1
      c(0) = 0
c
      do q = 0, t - 1
        n = n / 2
        do j = 0, k - 1
          c(k + j) = c(j) + n
        end do
        k = 2 * k
      end do
c
      return
      END
