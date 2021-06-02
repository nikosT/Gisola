      SUBROUTINE JACOBI(N, ITMAX, EPS1, EPS2, EPS3, A, T, EIGEN, AIK)

      IMPLICIT NONE

!$OMP DECLARE TARGET
!$ACC ROUTINE SEQ
      REAL A, AIK, CSA, EIGEN, EPS1, EPS2, EPS3, HOLDIK, HOLDKI, OFFDSQ, P, Q, S, SIGMA1, SIGMA2, SNA, SPQ, T
      INTEGER I, IP1, ITER, ITMAX, J, K, N, NM1
!
      DIMENSION A(3, 3), T(3, 3), AIK(3), EIGEN(3)

      DO I = 1, N
         DO J = 1, N
            T(I, J) = 0.
         ENDDO
      ENDDO
      NM1 = N - 1
      SIGMA1 = 0.
      OFFDSQ = 0.
      DO I = 1, N
         SIGMA1 = SIGMA1 + A(I, I)**2
         T(I, I) = 1.
         IP1 = I + 1
         IF (I>=N) EXIT
         DO J = IP1, N
            OFFDSQ = OFFDSQ + A(I, J)**2
         ENDDO
      ENDDO
      S = 2.*OFFDSQ + SIGMA1
!
      DO ITER = 1, ITMAX
         DO I = 1, NM1
            IP1 = I + 1
            DO J = IP1, N
               Q = ABS(A(I, I) - A(J, J))
               IF (Q>EPS1) THEN
                  IF (ABS(A(I, J))<=EPS2) GOTO 10
                  P = 2.*A(I, J)*Q/(A(I, I) - A(J, J))
                  SPQ = SQRT(P*P + Q*Q)
                  CSA = SQRT((1. + Q/SPQ)*.5)
                  SNA = P/(2.*CSA*SPQ)
               ELSE
                  CSA = .707106781186547
                  SNA = CSA
               ENDIF
               DO K = 1, N
                  HOLDKI = T(K, I)
                  T(K, I) = HOLDKI*CSA + T(K, J)*SNA
                  T(K, J) = HOLDKI*SNA - T(K, J)*CSA
               ENDDO
               DO K = 1, N
                  IF (K<=J) THEN
                     AIK(K) = A(I, K)
                     A(I, K) = CSA*AIK(K) + SNA*A(K, J)
                     IF (K==J) A(J, K) = SNA*AIK(K) - CSA*A(J, K)
                  ELSE
                     HOLDIK = A(I, K)
                     A(I, K) = CSA*HOLDIK + SNA*A(J, K)
                     A(J, K) = SNA*HOLDIK - CSA*A(J, K)
                  ENDIF
               ENDDO
               AIK(J) = SNA*AIK(I) - CSA*AIK(J)
               DO K = 1, J
                  IF (K>I) THEN
                     A(K, J) = SNA*AIK(K) - CSA*A(K, J)
                  ELSE
                     HOLDKI = A(K, I)
                     A(K, I) = CSA*HOLDKI + SNA*A(K, J)
                     A(K, J) = SNA*HOLDKI - CSA*A(K, J)
                  ENDIF
               ENDDO
 10            A(I, J) = 0.
            ENDDO
         ENDDO
         SIGMA2 = 0.0D+0
         DO I = 1, N
            EIGEN(I) = A(I, I)
            SIGMA2 = SIGMA2 + EIGEN(I)**2
         ENDDO
!       write(*,*) eigen(1), eigen(2), eigen(3)
!        write(*,*) "here"
!        write(*,*) sigma1,sigma2,sigma1/sigma2,eps3
         IF (1.0 - SIGMA1/SIGMA2<EPS3) GOTO 99999
         SIGMA1 = SIGMA2
      ENDDO
!!!!!      WRITE (*, '(//'' No convergence with '',I3,'' iterations''//)') ITMAX
!
!     write(*,*) 'jacobi'
99999 END SUBROUTINE JACOBI
