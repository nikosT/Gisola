      SUBROUTINE LUBKSB(A, N, NP, INDX, B)

      IMPLICIT NONE
!$OMP DECLARE TARGET
!$ACC ROUTINE SEQ
      REAL A, B, SUM
      INTEGER I, II, INDX, J, LL, N, NP
      DIMENSION A(NP, NP), INDX(N), B(N)
      II = 0
      DO I = 1, N
         LL = INDX(I)
         SUM = B(LL)
         B(LL) = B(I)
         IF (II/=0) THEN
            DO J = II, I - 1
               SUM = SUM - A(I, J)*B(J)
            ENDDO
         ELSEIF (SUM/=0.) THEN
            II = I
         ENDIF
         B(I) = SUM
      ENDDO
      DO I = N, 1, -1
         SUM = B(I)
         DO J = I + 1, N
            SUM = SUM - A(I, J)*B(J)
         ENDDO
         B(I) = SUM/A(I, I)
      ENDDO
      END SUBROUTINE LUBKSB
