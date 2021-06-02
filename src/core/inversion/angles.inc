      SUBROUTINE ANGLES(ANN1, ANN2, DIP, STRIKE, RAKE)
      IMPLICIT NONE
!$OMP DECLARE TARGET
!$ACC ROUTINE SEQ
      REAL AN1, AN2, ANN1, ANN2, CF, CR, DIP, EPS, FAC, RAKE, SDI, SF, STRIKE
      INTEGER J
!
!     procedure ANGLES determines dip, strike and rake angles from
!       nodal plane normals;
!     to obtain both the solutions it has to be called twice:
!       CALL ANGLES (AN1,AN2,D1,S1,R1)
!       CALL ANGLES (AN2,AN1,D2,S2,R2)
!
      DIMENSION ANN1(3), AN1(3), ANN2(3), AN2(3)

      EPS = 0.001
      FAC = 180./3.1415927
      IF (ANN1(3)>0.) THEN
         DO J = 1, 3
            AN2(J) = -ANN2(J)
            AN1(J) = -ANN1(J)
         ENDDO
      ELSE
         DO J = 1, 3
            AN2(J) = ANN2(J)
            AN1(J) = ANN1(J)
         ENDDO
      ENDIF
      IF (AN1(3)< - 0.99) THEN       !!!!!!!!!!!!!!!!
         DIP = ACOS(0.999)
      ELSE
         DIP = ACOS( - AN1(3))
      ENDIF      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (ABS(ABS(AN1(3)) - 1.)<EPS) THEN
         RAKE = 0.
         STRIKE = ATAN2(AN2(2), AN2(1))
         IF (STRIKE<0.) STRIKE = STRIKE + 6.2831853
      ELSE
         SDI = 1./SQRT(1. - AN1(3)*AN1(3))
         STRIKE = ATAN2( - AN1(1)*SDI, AN1(2)*SDI)
         IF (STRIKE<0.) STRIKE = STRIKE + 6.2831853
         IF (ABS(AN1(3))>=EPS) THEN
            CF = COS(STRIKE)
            SF = SIN(STRIKE)
            RAKE = ATAN2((AN2(1)*SF - AN2(2)*CF)/( - AN1(3)), AN2(1)*CF + AN2(2)*SF)
         ELSEIF (ABS(STRIKE)<EPS) THEN
            RAKE = ATAN2( - AN2(3), AN2(1))
         ELSEIF (ABS(ABS(STRIKE) - 1.5707963)<EPS) THEN
            RAKE = ATAN2( - AN2(3), AN2(2))
         ELSE
            CF = COS(STRIKE)
            SF = SIN(STRIKE)
            IF (ABS(CF)>ABS(SF)) THEN
               CR = AN2(1)/CF
            ELSE
               CR = AN2(2)/SF
            ENDIF
            RAKE = ATAN2( - AN2(3), CR)
         ENDIF
      ENDIF
      DIP = DIP*FAC
      STRIKE = STRIKE*FAC
      RAKE = RAKE*FAC
!
!     write(*,*) 'angles'
      END SUBROUTINE ANGLES
