!********************************************JACOBI  Num Rec *************
      SUBROUTINE JACOBINR(A, N, NP, D, V, NROT)

      USE COMMON_VARS

      IMPLICIT NONE

!$OMP DECLARE TARGET
!$ACC ROUTINE SEQ
      INTEGER N, NP, NROT, NMAX
      REAL A(NP, NP), D(NP), V(NP, NP)
      PARAMETER (NMAX = 500)
      INTEGER I, IP, IQ, J
      REAL*8 C, G, H, S, SM, T, TAU, THETA, TRESH, P(NMAX), Z(NMAX)
      DO IP = 1, N
         DO IQ = 1, N
            V(IP, IQ) = 0.
         ENDDO
         V(IP, IP) = 1.
      ENDDO
      DO IP = 1, N
         P(IP) = A(IP, IP)
         D(IP) = P(IP)
         Z(IP) = 0.
      ENDDO
      NROT = 0
      DO I = 1, 50
         SM = 0.
         DO IP = 1, N - 1
            DO IQ = IP + 1, N
               SM = SM + ABS(A(IP, IQ))
            ENDDO
         ENDDO
         IF (SM==0.) RETURN
         IF (I<4) THEN
            TRESH = 0.2*SM/N**2
         ELSE
            TRESH = 0.
         ENDIF
         DO IP = 1, N - 1
            DO IQ = IP + 1, N
               G = 100.*ABS(A(IP, IQ))
               IF ((I>4) .AND. (ABS(D(IP)) + G==ABS(D(IP))) .AND. (ABS(D(IQ)) + G==ABS(D(IQ)))) THEN
                  A(IP, IQ) = 0.
               ELSEIF (ABS(A(IP, IQ))>TRESH) THEN
                  H = D(IQ) - D(IP)
                  IF (ABS(H) + G==ABS(H)) THEN
                     T = A(IP, IQ)/H
                  ELSE
                     THETA = 0.5*H/A(IP, IQ)
                     T = 1./(ABS(THETA) + SQRT(1. + THETA**2))
                     IF (THETA<0.) T = -T
                  ENDIF
                  C = 1./SQRT(1 + T**2)
                  S = T*C
                  TAU = S/(1. + C)
                  H = T*A(IP, IQ)
                  Z(IP) = Z(IP) - H
                  Z(IQ) = Z(IQ) + H
                  D(IP) = D(IP) - H
                  D(IQ) = D(IQ) + H
                  A(IP, IQ) = 0.
                  DO J = 1, IP - 1
                     G = A(J, IP)
                     H = A(J, IQ)
                     A(J, IP) = G - S*(H + G*TAU)
                     A(J, IQ) = H + S*(G - H*TAU)
                  ENDDO
                  DO J = IP + 1, IQ - 1
                     G = A(IP, J)
                     H = A(J, IQ)
                     A(IP, J) = G - S*(H + G*TAU)
                     A(J, IQ) = H + S*(G - H*TAU)
                  ENDDO
                  DO J = IQ + 1, N
                     G = A(IP, J)
                     H = A(IQ, J)
                     A(IP, J) = G - S*(H + G*TAU)
                     A(IQ, J) = H + S*(G - H*TAU)
                  ENDDO
                  DO J = 1, N
                     G = V(J, IP)
                     H = V(J, IQ)
                     V(J, IP) = G - S*(H + G*TAU)
                     V(J, IQ) = H + S*(G - H*TAU)
                  ENDDO
                  NROT = NROT + 1
               ENDIF
            ENDDO
         ENDDO
         DO IP = 1, N
            P(IP) = P(IP) + Z(IP)
            D(IP) = P(IP)
            Z(IP) = 0.
         ENDDO
      ENDDO
!$OMP ATOMIC UPDATE
!$ACC ATOMIC UPDATE
      ERROR_ON_DEVICE = IOR(ERROR_ON_DEVICE, 1)
      END SUBROUTINE JACOBINR
