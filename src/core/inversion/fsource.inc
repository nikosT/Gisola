! @(#) fsource.F       FSOURCE 1.4      11/17/92 1
!**********************************************************
!       FSOURCE
!
!       Different kind of source function defined in the
!       frequency domain
!       Source function for dislocation (step, ramp, haskell)
!       are normalized so that in far-field the low-frequency
!       level is proportionnal to the seismic moment with a factor
!       equals to: Rad/(4.PI.rho.beta^3) * 1./r
!       where Rad= Radiation coefficient with (possibly)
!       free surface effect
!
!       input:
!               type    ->      see below
!               omega   ->      angular frequency
!               t0,t1   ->      time constant when needed
!               dt      ->      sampling rate
!**********************************************************

      FUNCTION FSOURCE(TYPE, OMEGA, T0, T1, DT)

      USE PARAMETERS

      IMPLICIT NONE

      INTEGER TYPE
      REAL DT, T0, T1
      REAL*8 UUR, UUI, TRISE, TRUPT
      COMPLEX*16 FSOURCE, UU, UEX, UXX, OMEGA, SHX, OMEGAC

! TYPE=0        Source = Dirac en deplacement
      IF (TYPE==0) FSOURCE = 1

! TYPE=1        Source = Ricker en deplacement
      IF (TYPE==1) THEN
         UU = OMEGA*T0
         UU = UU*UU/PI2/PI2
         UU = EXP( - UU)
         UU = OMEGA*OMEGA*UU*DT
         FSOURCE = UU
      ENDIF

! TYPE=2        Source = step en deplacement
!               2 steps possibles (1) real=1/(ai*omega)
!                                 (2) bouchon's
      IF (TYPE==2) THEN
         SHX = EXP(OMEGA*PI*T0/2.)     !Bouchon's
         SHX = 1./(SHX - 1./SHX)
         UU = -AI*T0*PI*SHX
         FSOURCE = UU
      ENDIF

! TYPE=7        Source = step en deplacement
      IF (TYPE==7) THEN
         UU = 1./AI/OMEGA
         FSOURCE = UU
      ENDIF

! TYPE=3        Source = file
      IF (TYPE==3) THEN
         READ (130, *) UUR, UUI
         FSOURCE = DCMPLX(UUR, UUI)
      ENDIF

! TYPE=4        Source = triangle en deplacement
      IF (TYPE==4) THEN
         UU = EXP(AI*OMEGA*T0/4.)
         UU = (UU - 1./UU)/2./AI
         UU = UU/(OMEGA*T0/2.)
         FSOURCE = UU*UU*4.
      ENDIF

! TYPE=5        Source = rampe causale
!               rise time T=t0
      IF (TYPE==5) THEN
         TRISE = T0
         UU = AI*OMEGA*TRISE
         UU = (1. - EXP( - UU))/UU
         FSOURCE = UU/(AI*OMEGA)
      ENDIF

! TYPE=9        Source = Brune
!               corner freq= 1 / t0
      IF (TYPE==9) THEN
         TRISE = T0
         OMEGAC = PI2*(1./TRISE)
         UU = OMEGAC + AI*OMEGA
         FSOURCE = OMEGAC*OMEGAC/(UU*UU)
      ENDIF

! TYPE=6,8      Source = modele d'haskell, trapezoide
!                        1 ere cste de temps rise time: riset
!                        2 eme cste de temps, duree de la rupture
!                        trupt = Length/2/rupt_velocity (Haskell)
      IF ((TYPE==6) .OR. (TYPE==8)) THEN
         TRISE = T0
         TRUPT = T1
         UU = AI*OMEGA*TRISE
         UU = (1. - EXP( - UU))/UU     ! ramp
         UXX = AI*OMEGA*TRUPT/2.       ! finite fault
         UEX = EXP(UXX)
         UXX = (UEX - 1./UEX)/UXX/2.
         FSOURCE = UU*UXX/(AI*OMEGA)
      ENDIF

      END FUNCTION FSOURCE

