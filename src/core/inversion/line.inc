      SUBROUTINE LINE(A, T, N)

      IMPLICIT NONE
!$OMP DECLARE TARGET
!$ACC ROUTINE SEQ
      REAL A, AM, P9, T
      INTEGER J, K, KK, N
      DIMENSION A(3), T(3, 3)

      DO KK = 1, N
         AM = A(KK)
         DO K = KK, N
            IF (A(K)>AM) AM = A(K)
         ENDDO
         DO K = KK, N
            IF (A(K)==AM) EXIT
         ENDDO
         A(K) = A(KK)
         A(KK) = AM
         DO J = 1, 3
            P9 = T(J, K)
            T(J, K) = T(J, KK)
            T(J, KK) = P9
         ENDDO
      ENDDO
!
!     write(*,*) 'line'
      END SUBROUTINE LINE
