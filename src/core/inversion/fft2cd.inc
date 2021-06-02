      SUBROUTINE FFT2CD(A, M, IWK)

      IMPLICIT NONE

!     specifications for arguments
      INTEGER M, IWK(*)
      COMPLEX*16 A(*)

!     specifications for local variables
      INTEGER I, ISP, J, JJ, JSP, K, K0, K1, K2, K3, KB, KN, MK, MM, MP, N, N4, N8, N2, LM, NN, JK
      REAL*8 RAD, C1, C2, C3, S1, S2, S3, CK, SK, SQ, A0, A1, A2, A3, B0, B1, B2, B3, TWOPI, TEMP, ZERO, ONE
      REAL*8 Z0(2), Z1(2), Z2(2), Z3(2)
      COMPLEX*16 ZA0, ZA1, ZA2, ZA3, AK2
      EQUIVALENCE (ZA0, Z0(1))
      EQUIVALENCE (ZA1, Z1(1))
      EQUIVALENCE (ZA2, Z2(1))
      EQUIVALENCE (ZA3, Z3(1))
      EQUIVALENCE (A0, Z0(1))
      EQUIVALENCE (B0, Z0(2))
      EQUIVALENCE (A1, Z1(1))
      EQUIVALENCE (B1, Z1(2))
      EQUIVALENCE (A2, Z2(1))
      EQUIVALENCE (B2, Z2(2))
      EQUIVALENCE (A3, Z3(1))
      EQUIVALENCE (B3, Z3(2))
      DATA SQ, SK, CK, TWOPI/.7071068, .3826834, .9238795, 6.283185/
      DATA ZERO/0.0/, ONE/1.0/

      MP = M + 1
      N = 2**M
      IWK(1) = 1
      MM = (M/2)*2
      KN = N + 1

!     initialize work vector
      DO I = 2, MP
         IWK(I) = IWK(I - 1) + IWK(I - 1)
      ENDDO
      RAD = TWOPI/N
      MK = M - 4
      KB = 1
      IF (MM/=M) THEN
         K2 = KN
         K0 = IWK(MM + 1) + KB
         DO
            K2 = K2 - 1
            K0 = K0 - 1
            AK2 = A(K2)
            A(K2) = A(K0) - AK2
            A(K0) = A(K0) + AK2
            IF (K0<=KB) EXIT
         ENDDO
      ENDIF
      C1 = ONE
      S1 = ZERO
      JJ = 0
      K = MM - 1
      J = 4
      IF (K<1) GOTO 500
      GOTO 200
 100  JJ = IWK(J) + JJ
      J = 4
 200  ISP = IWK(K)
      IF (JJ==0) GOTO 400

!     reset trigonometric parameter(s)
      C2 = JJ*ISP*RAD
      C1 = COS(C2)
      S1 = SIN(C2)
 300  C2 = C1*C1 - S1*S1
      S2 = C1*(S1 + S1)
      C3 = C2*C1 - S2*S1
      S3 = C2*S1 + S2*C1
 400  JSP = ISP + KB

!     determine fourier coefficients in groups of 4
      DO I = 1, ISP
         K0 = JSP - I
         K1 = K0 + ISP
         K2 = K1 + ISP
         K3 = K2 + ISP
         ZA0 = A(K0)
         ZA1 = A(K1)
         ZA2 = A(K2)
         ZA3 = A(K3)
         IF (S1/=ZERO) THEN
            TEMP = A1
            A1 = A1*C1 - B1*S1
            B1 = TEMP*S1 + B1*C1
            TEMP = A2
            A2 = A2*C2 - B2*S2
            B2 = TEMP*S2 + B2*C2
            TEMP = A3
            A3 = A3*C3 - B3*S3
            B3 = TEMP*S3 + B3*C3
         ENDIF
         TEMP = A0 + A2
         A2 = A0 - A2
         A0 = TEMP
         TEMP = A1 + A3
         A3 = A1 - A3
         A1 = TEMP
         TEMP = B0 + B2
         B2 = B0 - B2
         B0 = TEMP
         TEMP = B1 + B3
         B3 = B1 - B3
         B1 = TEMP
         A(K0) = DCMPLX(A0 + A1, B0 + B1)
         A(K1) = DCMPLX(A0 - A1, B0 - B1)
         A(K2) = DCMPLX(A2 - B3, B2 + A3)
         A(K3) = DCMPLX(A2 + B3, B2 - A3)
      ENDDO
      IF (K<=1) THEN
         KB = K3 + ISP

!        check for completion of final iteration
         IF (KN>KB) THEN
            IF (J/=1) THEN
               J = J - 1
               C2 = C1
               IF (J/=2) THEN
                  C1 = (C1 - S1)*SQ
                  S1 = (C2 + S1)*SQ
               ELSE
                  C1 = C1*CK + S1*SK
                  S1 = S1*CK - C2*SK
               ENDIF
               GOTO 300
            ELSE
               K = 3
               J = MK
               DO WHILE (IWK(J)<=JJ)
                  JJ = JJ - IWK(J)
                  J = J - 1
                  IF (IWK(J)>JJ) EXIT
                  JJ = JJ - IWK(J)
                  J = J - 1
                  K = K + 2
               ENDDO
               GOTO 100
            ENDIF
         ENDIF
      ELSE
         K = K - 2
         GOTO 200
      ENDIF

!     permute the complex vector in reverse binary order to normal order
 500  IF (M<=1) GOTO 99999
      MP = M + 1
      JJ = 1

!     initialize work vector
      IWK(1) = 1
      DO I = 2, MP
         IWK(I) = IWK(I - 1)*2
      ENDDO
      N4 = IWK(MP - 2)
      IF (M>2) N8 = IWK(MP - 3)
      N2 = IWK(MP - 1)
      LM = N2
      NN = IWK(MP) + 1
      MP = MP - 4

!     determine indices and switch a
      J = 2
 600  JK = JJ + N2
      AK2 = A(J)
      A(J) = A(JK)
      A(JK) = AK2
      J = J + 1
      IF (JJ>N4) THEN
         JJ = JJ - N4
         IF (JJ>N8) THEN
            JJ = JJ - N8
            K = MP
            DO WHILE (IWK(K)<JJ)
               JJ = JJ - IWK(K)
               K = K - 1
            ENDDO
            JJ = IWK(K) + JJ
         ELSE
            JJ = JJ + N8
         ENDIF
      ELSE
         JJ = JJ + N4
      ENDIF
      IF (JJ>J) THEN
         K = NN - J
         JK = NN - JJ
         AK2 = A(J)
         A(J) = A(JJ)
         A(JJ) = AK2
         AK2 = A(K)
         A(K) = A(JK)
         A(JK) = AK2
      ENDIF
      J = J + 1

!     cycle repeated until limiting number of changes is achieved
      IF (J<=LM) GOTO 600

99999 END SUBROUTINE FFT2CD

