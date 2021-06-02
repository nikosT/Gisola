      SUBROUTINE FILTER15(KEYFIL, KEYDIS, SHIFT, F1, F4, TIMSTEP, B)

      USE PARAMETERS

      IMPLICIT NONE

      INCLUDE "interface/interface_timefilters_xapiir.inc"

!$OMP DECLARE TARGET
!$ACC ROUTINE SEQ
      REAL AA, B, B0, B1, BASE, DT, F1, F4, SHIFT, SUM1, SUM2, TIMSTEP, TREND, X2, XAUX
      INTEGER I, ISH, ITIM, KEYDIS, KEYFIL, MAX_NT, NBASE, NTIM, NTREND
 
! Time domain causal filter using XAPIIR
 
      DIMENSION B(NUM_OF_TIME_SAMPLES)
      DIMENSION XAUX(NUM_OF_TIME_SAMPLES), TREND(NUM_OF_TIME_SAMPLES)
!      DIMENSION AA( - 100000:100000)   ! for shift ??
     DIMENSION AA(10000)
! declaration needed for XAPIIR
      CHARACTER*2 TYPE, APROTO
      INTEGER NSAMPS, PASSES, IORD
      REAL*4 TRBNDW, A, FLO, FHI, TS
 
      DT = TIMSTEP
      NTIM = NUM_OF_TIME_SAMPLES
 
!         keyfil= ... band pass (0) or low pass (1) or high pass (2)
!     keydis= ... output is velocity (0) or displacement (1)
!     shift=  ... shift (in seconds) to the right
 
!
!     READING VELOCITY and SHHIFTING 'to the right'
!

      ISH = IFIX(SHIFT/DT)  ! shift is in parameters

      DO I = 1, NTIM
         AA(I + ISH) = B(I)
      ENDDO
 
 
!
!     FILTRATION IN THE TIME DOMAIN  (using XAPIIR by D. Harris 1990)
!
 
      NSAMPS = NTIM
      APROTO = 'BU' ! Butterworth
      TRBNDW = 0.        ! not used for 'BU'
      A = 0.                ! not used for 'BU'
      IORD = 4       ! number of poles (4-5 recommended) the 4 for compatibility with ISOLA GUI Butt amplitude response 4th, power 8th
      TS = DT
      PASSES = 1       ! (1 for forward; 2 for forward and reverse = zero phase)
      MAX_NT = NUM_OF_TIME_SAMPLES
 
      IF (KEYFIL==0) TYPE = 'BP'        ! band pass
      IF (KEYFIL==1) TYPE = 'LP'        ! low pass
      IF (KEYFIL==2) TYPE = 'HP'        ! high pass
!              write(*,*) type
      FLO = F1       !   ignored for LP
      IF (FLO==0.) FLO = 1./(DT*REAL(NUM_OF_TIME_SAMPLES))   ! BandPass cannot start at 0
      FHI = F4
 
      DO ITIM = 1, NTIM
         XAUX(ITIM) = AA(ITIM)
                             ! auxiliary xaux needed due to dimension of aa
      ENDDO
 
      CALL XAPIIR(XAUX, NSAMPS, APROTO, TRBNDW, A, IORD, TYPE, FLO, FHI, TS, PASSES, MAX_NT)
 
 
      DO ITIM = 1, NTIM
         AA(ITIM) = XAUX(ITIM)
      ENDDO
 
 
!
!     LINEAR TREND correction of velocity
!
 
      SUM1 = 0.
      SUM2 = 0.
!     ntrend=ifix(40./dt) from 40 seconds ? (consider shift and travel time)
      NTREND = NTIM  ! in samples; from complete record ??
      DO I = 1, NTREND
         SUM1 = SUM1 + AA(I)
         SUM2 = SUM2 + AA(I)*FLOAT(I)
      ENDDO
      B0 = (2.*(2.*FLOAT(NTREND) + 1)*SUM1 - 6.*SUM2)/(FLOAT(NTREND)*(FLOAT(NTREND) - 1.))
      B1 = (12.*SUM2 - 6.*(FLOAT(NTREND) + 1)*SUM1)/(DT*FLOAT(NTREND)*(FLOAT(NTREND) - 1)*(FLOAT(NTREND) + 1))
      DO I = 1, NTIM
                  ! now applied to all samples
         TREND(I) = B0 + B1*DT*FLOAT(I)
      ENDDO
      DO I = 1, NTIM
!      aa(i)=aa(i)-trend(i) ! trend correction enabled
         AA(I) = AA(I)      ! trend correction disabled
      ENDDO
 
 
!
!   BASELINE correction of velocity
!
 
!         xlim2=10. ! in seconds ??? (consider shift and travel time)
!         nbase=ifix(xlim2/dt)
      NBASE = 50
                ! samples ??
      BASE = 0.
      DO I = 1, NBASE
         BASE = BASE + AA(I)
      ENDDO
      BASE = BASE/FLOAT(NBASE)
      DO I = 1, NTIM
!          aa(i)=aa(i)-base ! baseline correction enabled
         AA(I) = AA(I)        ! baseline correction disabled
      ENDDO
 
 
 
!
!     INTEGRATION  in time domain
!
 
      IF (KEYDIS/=0) THEN
         X2 = 0.       ! from VELOCITY to DISPLACEMENT
         DO I = 1, NTIM
            X2 = X2 + AA(I)*DT
            AA(I) = X2
         ENDDO
      ENDIF
 
!
!     LINEAR TREND CORRECTION of displacement
!
 
      SUM1 = 0.
      SUM2 = 0.
      NTREND = IFIX(40./DT)
                          ! trend from part between P wave, 40s ????????(consider shift and travel time)
      DO I = 1, NTREND
         SUM1 = SUM1 + AA(I)
         SUM2 = SUM2 + AA(I)*FLOAT(I)
      ENDDO
      B0 = (2.*(2.*FLOAT(NTREND) + 1)*SUM1 - 6.*SUM2)/(FLOAT(NTREND)*(FLOAT(NTREND) - 1.))
      B1 = (12.*SUM2 - 6.*(FLOAT(NTREND) + 1)*SUM1)/(DT*FLOAT(NTREND)*(FLOAT(NTREND) - 1)*(FLOAT(NTREND) + 1))
      DO I = 1, NTIM
                  ! now applied to all samples
         TREND(I) = B0 + B1*DT*FLOAT(I)
      ENDDO
      DO I = 1, NTIM
!     aa(i)=aa(i)-trend(i) ! trend correction enabled
         AA(I) = AA(I)     ! trend correction disabled
      ENDDO
 
 
      DO I = 1, NTIM
!         b(i)=aa(i-ish) ! Possible backward time shift
         B(I) = AA(I)
      ENDDO
 
 
      END SUBROUTINE FILTER15
 
      INCLUDE "timefilters.inc"

