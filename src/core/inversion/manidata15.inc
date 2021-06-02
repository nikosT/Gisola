      SUBROUTINE MANIDATA15(KEYGPS, X, RR)

      USE COMMON_VARS
      USE PARAMETERS

      IMPLICIT NONE

      INCLUDE "interface/interface_filter15.inc"

      REAL F1, F4, RR, SHIFT, TIME, X
      INTEGER ICOM, IR, ITIM, KEYDIS, KEYFIL, KEYGPS, NFILE
 
! Manipulating the data. Input data x, output data also x (re-written).
! Power of the data rr.
 
      DIMENSION X(NUM_OF_TIME_SAMPLES, MAX_STATIONS, 3)
      DIMENSION KEYGPS(MAX_STATIONS)
 
!      CHARACTER *5 statname(MAX_STATIONS)
!      character *17 statfil1,statfil2,statfil3  !!g77
 
! *******************************************************************
! *************MANIPULATING OBSERVED DATA****************************
! *******************************************************************
 
      DO IR = 1, NR
         DO ITIM = 1, NTIM
            X(ITIM, IR, 1) = 0.
            X(ITIM, IR, 2) = 0.
            X(ITIM, IR, 3) = 0.
         ENDDO
      ENDDO
 
!
!     READING velocity SEISMOGRAMS TO BE INVERTED (=DATA)
!
 
      DO IR = 1, NR
         NFILE = 1000 + 1*IR
!      statfil1=trim(statname(ir))//'raw.dat'
!      open(nfile,file=statfil1)
         DO ITIM = 1, NTIM
            READ (NFILE, *) TIME, X(ITIM, IR, 1), X(ITIM, IR, 2), X(ITIM, IR, 3)
         ENDDO
      CLOSE(NFILE)
      ENDDO
 
 
 
 
!
!     FILTERING THE DATA - OLD !!!
!
!      do ir=1,nr
!        do icom=1,3
!         call filter(f1,f2,f3,f4,dt,x(1,ir,icom))
!        enddo
!      enddo
 
 
 
 
!
!     FILTERING THE DATA
!     GPS data supposed to be in ***raw as shifted by 50s, LF displacement
!
 
 
      DO IR = 1, NR
         IF (KEYGPS(IR)==0) THEN
                               !!!!!! filter only for SEIS; but suppressed for gps
            F1 = FF1(IR)
                  !not used for LP!!!! 'filter' involves 50s shift, filtr, integration, stabilizations
            F4 = FF4(IR)
                  !not used for HP!!!!
            KEYFIL = 0
               !(BP=0, LP=1, HP=2)   ! LP to avoid double suppression of low freq after pre-filtering in GUI
            KEYDIS = 1
               !(vel=0, displ=1)
            SHIFT = 0.
               ! shift 50 s is impossible if TL<50s !!!!!!!!!
            DO ICOM = 1, 3
               CALL FILTER15(KEYFIL, KEYDIS, SHIFT, F1, F4, DT, X(1, IR, ICOM))
                                                                     ! seis. has the OLD filter
!         call filter(f1,f2,f3,f4,dt,x(1,ir,icom))   ! seis. has the OLD filter
!        call filter_stat(f1,f2,f3,f4,dt,x(1,ir,icom)) ! this is a NEW ('clever') filter to support statics
! CAUTION same option as here (the old filter or new) must be also in elemat2_12c_mix
            ENDDO
         ENDIF !!!!!!!!!!!!!!!!!!!
      ENDDO
 
 
      DO IR = 1, NR       ! output of filtered input data
         NFILE = 2000 + 1*IR
!      statfil2=trim(statname(ir))//'fil.dat'
!      open(nfile,file=statfil2)
         DO ITIM = 1, NTIM
            TIME = FLOAT(ITIM - 1)*DT
            WRITE (NFILE, '(4(1x,e12.6))') TIME, X(ITIM, IR, 1), X(ITIM, IR, 2), X(ITIM, IR, 3)
         ENDDO
!         close(nfile)
      ENDDO
 
 
!
!     COMPUTING 'POWER' OF THE original DATA (the USED stations ONLY !)
!
 
      RR = 0.
      DO ICOM = 1, 3                                    !icom=seismo comp.
         DO IR = 1, NR                                  !  ir=station
            IF (STAT(IR)) THEN
               DO ITIM = 1, NTIM
                         ! ntim or 4000! For ACRE can change ntim here to 4000 and in ISOLA15 to get VR
 
!          rr=rr+x(itim,ir,icom)*x(itim,ir,icom)     ! old (no weights)
                  RR = RR + (X(ITIM, IR, ICOM)*WEIG(IR, ICOM))**2
                                                     ! new (with weights) 9.9. 2015
 
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      RR = RR*DT
 
      END SUBROUTINE MANIDATA15
