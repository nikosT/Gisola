      SUBROUTINE SILSUB(A, STR1, DIP1, RAKE1, STR2, DIP2, RAKE2, AMT, ADC, AMV1)

      IMPLICIT NONE

      INCLUDE "interface/interface_ang.inc"
      INCLUDE "interface/interface_angles.inc"
      INCLUDE "interface/interface_jacobi.inc"
      INCLUDE "interface/interface_line.inc"

!$OMP DECLARE TARGET
!$ACC ROUTINE SEQ
      REAL A, ACLVD, AD, ADC, AK, AL, AL1, ALM, AM, AMM, AMT, AMV, AMV1, AN1, AN2, ANGDL, ANGDLL, ANGLL, ANGLL1
      REAL ANGLL2, DIP1, DIP2, EN, EN1, EN1MAX, EN1MIN, EPS, EPS1, FAC, FACM, HLP, RAKE1, RAKE2, SQR1, SQR2, SQR3, STR1, STR2, T
      INTEGER I, J, N
!     Author: J. Sileny
 
 
!     COMMON/COMSGN/amm(3,3)
      DIMENSION AMM(3, 3)
                         ! moje zmena
      DIMENSION AM(3, 3), T(3, 3), EN(3), EN1(3)
      DIMENSION AN1(3), AN2(3), AK(3), HLP(3), AD(3)
      DIMENSION AL(3), AL1(3), ALM(3)
      DIMENSION A(6)

!
      FAC = 3.14159/180.
!
!     open(2,file='moment.lst')
 
!     imom=6            ! a( ) enter as parameters
!     do n=1,imom
!     read(111,*) a(n)
!     enddo
 
!     imom=5
!     do n=1,imom
!     read(111,*) a(n)
!     enddo
!     a(6)=0.
 
 
      AM(1, 1) = -1.*A(4) + A(6)
      AM(2, 2) = -1.*A(5) + A(6)
      AM(3, 3) = A(4) + A(5) + A(6)
      AM(1, 2) = A(1)
      AM(1, 3) = A(2)
      AM(2, 3) = -A(3)
 
 
      AM(2, 1) = AM(1, 2)
      AM(3, 1) = AM(1, 3)
      AM(3, 2) = AM(2, 3)
!
!     gisola: divide by 1000 and multiply by 1000 to avoid overflow
!     scaling of the tensor Mij
      AMT = SQRT(.5*((AM(1, 1)/1000)*(AM(1, 1)/1000) + (AM(2, 2)/1000)*(AM(2, 2)/1000) + &
            (AM(3, 3)/1000)*(AM(3, 3)/1000)) + (AM(1, 2)/1000)*(AM(1, 2)/1000) + (AM(1, 3)/1000)*(AM(1, 3)/1000) + (AM(2, 3)/1000) &
            *(AM(2, 3)/1000)) !gisola

      AMT=AMT*1000 !gisola

      FACM = 1./AMT/SQRT(2.)
      DO I = 1, 3
         DO J = 1, 3
            AM(I, J) = AM(I, J)*FACM
            AMM(I, J) = AM(I, J)
         ENDDO
      ENDDO
!      write(*,*) am
!
!     computation of eigenvalues and eigenfunctions, their ordering
      EPS1 = 1.0D-9
      CALL JACOBI(3, 50, EPS1, EPS1, EPS1, AM, T, EN, AK)
      CALL LINE(EN, T, 3)
!
!     create right-handed system from eigenvectors
      HLP(1) = T(2, 1)*T(3, 2) - T(3, 1)*T(1, 2)
      HLP(2) = T(3, 1)*T(1, 2) - T(1, 1)*T(3, 2)
      HLP(3) = T(1, 1)*T(2, 2) - T(2, 1)*T(1, 2)
      IF (HLP(1)*T(1, 3) + HLP(2)*T(2, 3) + HLP(3)*T(3, 3)<0.) THEN
         DO J = 1, 3
            T(J, 3) = -T(J, 3)
         ENDDO
      ENDIF
!
!
!     processing moment tensor
!
!     write(2,*)
!     write(2,*)
!     write(2,*)
!     write(2,*)
!     write(2,'(/'' M0:  '',e10.3)')amt
 
!     write(2,'(/'' MT (normalized):'')')
!     do 2 i=1,3
!   2 write(2,'(8g10.4)')(amm(i,j),j=1,3)
 
!     write(2,'(/'' MT eigenvalues:'')')
!     write(2,'(3e13.6)')(en(j),j=1,3)
 
!     write(2,'(/'' eigenvectors of moment tensor:'')')
!     write(2,'(8g10.4)')(t(j,1),j=1,3)
!     write(2,'(8g10.4)')(t(j,2),j=1,3)
!     write(2,'(8g10.4)')(t(j,3),j=1,3)
 
      AMV = EN(1) + EN(2) + EN(3)
      DO J = 1, 3
         EN1(J) = EN(J) - AMV/3.
      ENDDO
 
      EN1MAX = MAX(ABS(EN1(1)), ABS(EN1(2)), ABS(EN1(3)))
      EN1MIN = MIN(ABS(EN1(1)), ABS(EN1(2)), ABS(EN1(3)))
      IF (ABS(EN1MAX)<1.E-20) THEN
!!!!!         WRITE (*, *) 'problem in silsub'     !!!!!!!!!!!!!!!!!!
!!!!!         WRITE (*, *) 'MT is 100% ISO'     !!!!!!!!!!!!!!!!!!
!!!!!         WRITE (*, *) 's/d/r are undefined'     !!!!!!!!!!!!!!!!!!
         ADC = 0.
         AMV1 = 100.
         RETURN
      ELSE
         EPS = -EN1MIN/ABS(EN1MAX)
      ENDIF
 
!
!       construct vectors lambda_dev, d, l and l1'
      SQR2 = 1./SQRT(2.)
      SQR3 = 1./SQRT(3.)
      SQR1 = 1./SQRT(1.5)
      DO J = 1, 3
         ALM(J) = EN1(1)*T(J, 1) + EN1(2)*T(J, 2) + EN1(3)*T(J, 3)
         AD(J) = (T(J, 1) - T(J, 3))*SQR2
         AL(J) = (T(J, 1) - .5*T(J, 2) - .5*T(J, 3))*SQR1
         AL1(J) = (.5*T(J, 1) + .5*T(J, 2) - T(J, 3))*SQR1
      ENDDO
!     compute angles: lambda_dev^l(min), lambda_dev^d, d^l
      ANGDL = ANG(AD, ALM)
      ANGLL1 = ANG(AL, ALM)
      ANGLL2 = ANG(AL1, ALM)
      ANGDLL = ANG(AL, AD)
      ANGLL = MIN(ANGLL1, ANGLL2)
      ANGLL = ANGLL/ANGDLL
      ANGDL = ANGDL/ANGDLL
!
!     MT decomposition (1)
!     write(2,*)' '
!     write(2,*)' MT decomposition (1)'
      AMV1 = AMV**2/(EN(1)**2 + EN(2)**2 + EN(3)**2)/3.*100.
      IF (AMV>=0.) THEN
!     write(2,'('' V(explosive):'',t16,f5.1,''%'')')amv1
!     write(2,'('' V(implosive):'',t16,f5.1,''%'')')amv1
      ENDIF
      ACLVD = MIN(ABS(EN1(1)), ABS(EN1(2)), ABS(EN1(3)))/MAX(ABS(EN1(1)), ABS(EN1(2)), ABS(EN1(3)))*200.
!     write(2,'('' DC:'',t16,f5.1,''%'')')100.-aclvd
      IF (ANGLL1<ANGLL2) THEN
!     write(2,'('' CLVD(T-axis):'',t16,f5.1,''%'')')aclvd
!     write(2,'('' CLVD(P-axis):'',t16,f5.1,''%'')')aclvd
      ENDIF
 
!     MT decomposition (2)       ! this is routinely used in ISOLA (if we need decom (1), all this section should be 'c' except adc=
!     According eq. 8 of Vavrycuk, JGR,2001.
!     write(2,*)' '
!     write(2,*)' MT decomposition (2)'
      AMV1 = ABS(AMV)/MAX(ABS(EN(1)), ABS(EN(2)), ABS(EN(3)))/3.*100.
      IF (AMV>=0.) THEN
!     write(2,'('' V(explosive):'',t16,f5.1,''%'')')amv1
         AMV1 = AMV1
                  !!! new march2012
      ELSE
!     write(2,'('' V(implosive):'',t16,f5.1,''%'')')amv1
         AMV1 = -1.*AMV1
                      !!! VOL% is positive or negative; new march2012
      ENDIF
      ACLVD = 2.*ABS(EPS)*(100. - ABS(AMV1))
                                         !this is abs(CLVD%) !! new march2012
!c    adc=100.-amv1-aclvd              ! CAUTION even if we want decomp (1) this must stay here uncommented!!!
      ADC = 100. - ABS(AMV1) - ACLVD  !!! DC% always positive by definition; new march2012
 
 
!     write(2,'('' DC:'',t16,f5.1,''%'')')adc
      IF (ANGLL1<ANGLL2) THEN
!     write(2,'('' CLVD(T-axis):'',t16,f5.1,''%'')')aclvd
!     write(2,'('' CLVD(P-axis):'',t16,f5.1,''%'')')aclvd
      ENDIF
 
!     write(2,'(/'' T and P axes:'')')
!     write(2,'(8g10.4)')(t(j,1),j=1,3)
!     write(2,'(8g10.4)')(t(j,3),j=1,3)
 
      DO J = 1, 3
         AN1(J) = (T(J, 1) + T(J, 3))*.7071068
         AN2(J) = (T(J, 1) - T(J, 3))*.7071068
      ENDDO
!     write(2,'(/'' Nodal plane normals:'')')
!     write(2,'(8g10.4)')(an1(j),j=1,3)
!     write(2,'(8g10.4)')(an2(j),j=1,3)
 
      CALL ANGLES(AN1, AN2, DIP1, STR1, RAKE1)
      CALL ANGLES(AN2, AN1, DIP2, STR2, RAKE2)
!     write(2,'(/'' strike, dip and rake angles:'')')
!     write(2,'(8g10.4)')str1,dip1,rake1
!     write(2,'(8g10.4)')str2,dip2,rake2
!     close(2)
 
!cc   CALL RJPLT(EN,T)  (this subroutine, and its subr. 'proj' are DELETED
                                              ! full in sileny.for)
      END SUBROUTINE SILSUB
