      PROGRAM CORR_KAG
      IMPLICIT NONE
      REAL DIV, ROTSELECT, STAT, TABLE
 
!
      DIMENSION ROTSELECT(100000)
      DIMENSION DIV(200), TABLE(200), STAT(13)
      REAL*8 XMSTRIKE, XMDIP, XMRAKE, ROTANGLE, THETA, PHI
      REAL*8 XASTRIKE, XADIP, XARAKE
 
      CHARACTER(5) BUF1, BUF2, BUF3, BUF4, BUF5, BUF6
 
 
      CALL GETARG(1, BUF1)
      CALL GETARG(2, BUF2)
      CALL GETARG(3, BUF3)
 
      CALL GETARG(4, BUF4)
      CALL GETARG(5, BUF5)
      CALL GETARG(6, BUF6)
 
 
      READ (BUF1, *) XMSTRIKE
      READ (BUF2, *) XMDIP
      READ (BUF3, *) XMRAKE
 
      READ (BUF4, *) XASTRIKE
      READ (BUF5, *) XADIP
      READ (BUF6, *) XARAKE
 
!        write(*,*) xmstrike,xmdip,xmrake
!        write(*,*) xastrike,xadip,xarake
 
 
      IF (XMSTRIKE==XASTRIKE .AND. XMDIP==XADIP .AND. XMRAKE==XARAKE) THEN
         ROTANGLE = 0.0
         WRITE (*, *) ROTANGLE
      ELSE
         CALL MINDCROT(XMSTRIKE, XMDIP, XMRAKE, XASTRIKE, XADIP, XARAKE, ROTANGLE, THETA, PHI)
         WRITE (*, *) ROTANGLE
      ENDIF
 
 
      END PROGRAM CORR_KAG
 
 
 ! After Kagan, Y. Y. (1991). 3-D rotation of double-couple earthquake sources, Geophys. J. Int., 106(3), 709-716.
 
      SUBROUTINE MINDCROT(STRIKE1, DIP1, RAKE1, STRIKE2, DIP2, RAKE2, ROTANGLE, THETA, PHI)
    ! Provides minimum rotation between two DC mechanisms
    ! minimum rotation ROTANGLE along axis given by THETA and PHI
      IMPLICIT NONE
      REAL*8 STRIKE1, DIP1, RAKE1, STRIKE2, DIP2, RAKE2, ROTANGLE, THETA, PHI
      REAL*8, DIMENSION(4) :: Q1, Q2, QDUM, QDUM2, Q
      REAL*8 MAXVAL, DUM
      INTEGER ICODE, I
 
      CALL QUATFPS(Q1, STRIKE1, DIP1, RAKE1)
      CALL QUATFPS(Q2, STRIKE2, DIP2, RAKE2)
      MAXVAL = 180.
      DO I = 1, 4
         CALL F4R1(Q1, Q2, QDUM, I)
         CALL SPHCOOR(QDUM, ROTANGLE, THETA, PHI)
         IF (ROTANGLE<MAXVAL) THEN
            MAXVAL = ROTANGLE
            Q = QDUM
         ENDIF
      ENDDO
      CALL SPHCOOR(Q, ROTANGLE, THETA, PHI)
 
      END SUBROUTINE MINDCROT
 
 
 
      SUBROUTINE QUATFPS(QUAT, DDI, DAI, SAI)
    ! calculates rotation quaternion corresponding to given focal mechanism
    ! input: strike (DD), dip (DA), rake (SA)
    ! output: QUAT
      IMPLICIT NONE
      REAL*8 AN1, AN2, AN3, CDA, CDD, CSA, D2, DA, DAI, DD, DDI, ERR, P1, P2, P3, S1, S2, S3, SA, SAI
      REAL*8 SDA, SDD, SSA, T1, T2, T3, TEMP, U0, U1, U2, U3, UM, V1, V2, V3
      INTEGER IC, ICOD
      REAL*8, PARAMETER :: RAD = 57.295779513082320876798154814105D0
      REAL*8 QUAT(4)
 
      ERR = 1.D-15
      IC = 1
      DD = DDI/RAD
      DA = DAI/RAD
      SA = SAI/RAD
      CDD = DCOS(DD)
      SDD = DSIN(DD)
      CDA = DCOS(DA)
      SDA = DSIN(DA)
      CSA = DCOS(SA)
      SSA = DSIN(SA)
      S1 = CSA*SDD - SSA*CDA*CDD
      S2 = -CSA*CDD - SSA*CDA*SDD
      S3 = -SSA*SDA
      V1 = SDA*CDD
      V2 = SDA*SDD
      V3 = -CDA
 
      AN1 = S2*V3 - V2*S3
      AN2 = V1*S3 - S1*V3
      AN3 = S1*V2 - V1*S2
 
      D2 = 1.D0/DSQRT(2.D0)
      T1 = (V1 + S1)*D2
      T2 = (V2 + S2)*D2
      T3 = (V3 + S3)*D2
      P1 = (V1 - S1)*D2
      P2 = (V2 - S2)*D2
      P3 = (V3 - S3)*D2
 
      U0 = (T1 + P2 + AN3 + 1.D0)/4.D0
      U1 = (T1 - P2 - AN3 + 1.D0)/4.D0
      U2 = ( - T1 + P2 - AN3 + 1.D0)/4.D0
      U3 = ( - T1 - P2 + AN3 + 1.D0)/4.D0
      UM = DMAX1(U0, U1, U2, U3)
      IF (UM/=U0) THEN
         IF (UM==U1) THEN
 
            ICOD = 2*IC
            U1 = DSQRT(U1)
            U2 = (T2 + P1)/(4.D0*U1)
            U3 = (AN1 + T3)/(4.D0*U1)
            U0 = (P3 - AN2)/(4.D0*U1)
            GOTO 100
         ELSEIF (UM==U2) THEN
 
            ICOD = 3*IC
            U2 = DSQRT(U2)
            U1 = (T2 + P1)/(4.D0*U2)
            U0 = (AN1 - T3)/(4.D0*U2)
            U3 = (P3 + AN2)/(4.D0*U2)
            GOTO 100
         ELSEIF (UM==U3) THEN
 
            ICOD = 4*IC
            U3 = DSQRT(U3)
            U0 = (T2 - P1)/(4.D0*U3)
            U1 = (AN1 + T3)/(4.D0*U3)
            U2 = (P3 + AN2)/(4.D0*U3)
            GOTO 100
         ELSE
            WRITE (*, *) 'INTERNAL ERROR'
         ENDIF
      ENDIF
 
      ICOD = 1*IC
      U0 = DSQRT(U0)
      U3 = (T2 - P1)/(4.D0*U0)
      U2 = (AN1 - T3)/(4.D0*U0)
      U1 = (P3 - AN2)/(4.D0*U0)
 
 100  TEMP = U0*U0 + U1*U1 + U2*U2 + U3*U3
      IF (DABS(TEMP - 1.D0)>ERR) WRITE (*, *) 'INTERNAL ERROR'
 
      QUAT(1) = U1
      QUAT(2) = U2
      QUAT(3) = U3
      QUAT(4) = U0
 
      END SUBROUTINE QUATFPS
 
 
      SUBROUTINE SPHCOOR(QUAT, ANGL, THETA, PHI)
    !returns rotation angle (ANGL) of a counterclockwise rotation
    !and spherical coordinates (colatitude THETA and azimuth PHI) of the
    !rotation pole. THETA=0 corresponds to vector pointing down.
      IMPLICIT NONE
      REAL*8 Q4N
      REAL*8 QUAT(4), THETA, COSTH, PI, RAD2DEG, ANGL, PHI
 
      PI = ACOS( - 1.0)
      RAD2DEG = 180./PI
 
      IF (QUAT(4)<0.D0) QUAT(:) = -QUAT(:)
      Q4N = DSQRT(1.D0 - QUAT(4)**2)
      COSTH = 1.D0
      IF (DABS(Q4N)>1.D-10) COSTH = QUAT(3)/Q4N
      IF (DABS(COSTH)>1.D0) COSTH = IDINT(COSTH)
      THETA = ACOS(COSTH)*RAD2DEG
 
      ANGL = 2.D0*ACOS(QUAT(4))*RAD2DEG
 
      PHI = 0.D0
      IF (DABS(QUAT(1))>1.D-10 .OR. DABS(QUAT(2))>1.D-10) PHI = ATAN2(QUAT(2), QUAT(1))*RAD2DEG
      IF (PHI<0.D0) PHI = PHI + 360.D0
 
      END SUBROUTINE SPHCOOR
 
 
      SUBROUTINE QUATP(Q1, Q2, Q3)
    !calculates quaternion product Q3=Q2*Q1
      IMPLICIT NONE
      REAL*8, DIMENSION(4) :: Q1, Q2, Q3
 
      Q3(1) = Q1(4)*Q2(1) + Q1(3)*Q2(2) - Q1(2)*Q2(3) + Q1(1)*Q2(4)
      Q3(2) = -Q1(3)*Q2(1) + Q1(4)*Q2(2) + Q1(1)*Q2(3) + Q1(2)*Q2(4)
      Q3(3) = Q1(2)*Q2(1) - Q1(1)*Q2(2) + Q1(4)*Q2(3) + Q1(3)*Q2(4)
      Q3(4) = -Q1(1)*Q2(1) - Q1(2)*Q2(2) - Q1(3)*Q2(3) + Q1(4)*Q2(4)
 
      END SUBROUTINE QUATP
 
 
      SUBROUTINE QUATD(Q1, Q2, Q3)
    !quaternion division Q3=Q2*Q1^-1
      IMPLICIT NONE
      REAL*8, DIMENSION(4) :: Q1, QC1, Q2, Q3
 
      QC1(1:3) = -Q1(1:3)
      QC1(4) = Q1(4)
      CALL QUATP(QC1, Q2, Q3)
 
      END SUBROUTINE QUATD
 
 
      SUBROUTINE BOXTEST(Q1, Q2, QM, ICODE)
    !if ICODE==0 finds minimal rotation quaternion
    !if ICODE==N finds rotation quaternion Q2=Q1*(i,j,k,1) for N=(1,2,3,4)
      IMPLICIT NONE
      INTEGER ICODE, IXC
      REAL*8 QM
      REAL*8, DIMENSION(4) :: Q1, Q2, QUATT
      REAL*8 QUAT(4, 3)/1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0/
 
      IF (ICODE==0) THEN
         ICODE = 1
         QM = DABS(Q1(1))
         DO IXC = 2, 4
            IF (DABS(Q1(IXC))>QM) THEN
               QM = DABS(Q1(IXC))
               ICODE = IXC
            ENDIF
         ENDDO
      ENDIF
 
      IF (ICODE==4) THEN
         Q2 = Q1
      ELSE
         QUATT(:) = QUAT(:, ICODE)
         CALL QUATP(QUATT, Q1, Q2)
      ENDIF
 
      IF (Q2(4)<0.D0) Q2 = -Q2
 
      QM = Q2(4)
 
      END SUBROUTINE BOXTEST
 
 
      SUBROUTINE F4R1(Q1, Q2, Q, ICODE)
    ! Q=Q2*(Q1*(i,j,k,1))^-1 for N=(1,2,3,4)
    ! if N=0, then it finds it of the minimum
      IMPLICIT NONE
      INTEGER ICODE
      REAL*8 QM
      REAL*8, DIMENSION(4) :: Q, Q1, Q2, QR1
 
      CALL BOXTEST(Q1, QR1, QM, ICODE)
      CALL QUATD(QR1, Q2, Q)
 
      END SUBROUTINE F4R1
