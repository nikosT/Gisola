      SUBROUTINE LUDCMP(A, N, NP, INDX, D)

      USE COMMON_VARS

      IMPLICIT NONE

!$OMP DECLARE TARGET
!$ACC ROUTINE SEQ
      REAL A, AAMAX, D, DUM, SUM, TINY, VV
      INTEGER I, IMAX, INDX, J, K, N, NP
      PARAMETER (TINY = 1.0E-20)
      DIMENSION A(NP, NP), INDX(N), VV(400)
      D = 1.
      DO I = 1, N
         AAMAX = 0.
         DO J = 1, N
            IF (ABS(A(I, J))>AAMAX) AAMAX = ABS(A(I, J))
         ENDDO
         IF (AAMAX==0.) THEN
!$OMP ATOMIC UPDATE
!$ACC ATOMIC UPDATE
            ERROR_ON_DEVICE = IOR(ERROR_ON_DEVICE, 2)
         ENDIF
         VV(I) = 1./AAMAX
      ENDDO
      DO J = 1, N
         DO I = 1, J - 1
            SUM = A(I, J)
            DO K = 1, I - 1
               SUM = SUM - A(I, K)*A(K, J)
            ENDDO
            A(I, J) = SUM
         ENDDO
         AAMAX = 0.
         DO I = J, N
            SUM = A(I, J)
            DO K = 1, J - 1
               SUM = SUM - A(I, K)*A(K, J)
            ENDDO
            A(I, J) = SUM
            DUM = VV(I)*ABS(SUM)
            IF (DUM>=AAMAX) THEN
               IMAX = I
               AAMAX = DUM
            ENDIF
         ENDDO
         IF (J/=IMAX) THEN
            DO K = 1, N
               DUM = A(IMAX, K)
               A(IMAX, K) = A(J, K)
               A(J, K) = DUM
            ENDDO
            D = -D
            VV(IMAX) = VV(J)
         ENDIF
         INDX(J) = IMAX
         IF (A(J, J)==0.) A(J, J) = TINY
         IF (J/=N) THEN
            DUM = 1./A(J, J)
            DO I = J + 1, N
               A(I, J) = A(I, J)*DUM
            ENDDO
         ENDIF
      ENDDO
      END SUBROUTINE LUDCMP
