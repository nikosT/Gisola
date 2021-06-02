      SUBROUTINE NORM15_POSTVAR(allstatfile, F1, DT)

      USE PARAMETERS

      IMPLICIT NONE
      REAL ANORM1, ANORM2, ANORM3, DUREF, POSTVAR, TIME, VARRED, WEIG, X1, X2, X3
      INTEGER I, ICOM, IR, ITIM, NFILE, NR, NRUSED, NTIM, NUMBD, NUMBDF, NUMBP, NUSE
 
      REAL F1, DT
      CHARACTER *40 allstatfile

      DIMENSION X1(NUM_OF_TIME_SAMPLES, MAX_STATIONS, 3)
      DIMENSION X2(NUM_OF_TIME_SAMPLES, MAX_STATIONS, 3)
      DIMENSION X3(NUM_OF_TIME_SAMPLES, MAX_STATIONS, 3)
      DIMENSION NUSE(MAX_STATIONS), WEIG(3, MAX_STATIONS)
 
!     CHARACTER *12 filename
      CHARACTER*5 STATNAME(MAX_STATIONS)
      CHARACTER*17 STATFIL1, STATFIL2, STATFIL3


 
!      OPEN (99, FILE = 'inpinv.dat')
      OPEN (100, FILE = allstatfile)
      OPEN (101, FILE = 'inv4.dat')
 
 
!      DO I = 1, 3
!         READ (99, *)
!      ENDDO
!      READ (99, *) DT
!      DO I = 1, 9
!         READ (99, *)
!      ENDDO
!      READ (99, *) F1
!                    ! low freq limit of inpinv (not station variable)
 
 
      IR = 1
      DO
         READ (100, *, END = 100) STATNAME(IR), NUSE(IR), WEIG(1, IR), WEIG(2, IR), WEIG(3, IR)
         IR = IR + 1
         IF (IR>MAX_STATIONS) EXIT
      ENDDO
 100  NR = IR - 1
      WRITE (*, *) 'nr=', NR
!     if(nr.gt.MAX_STATIONS) then
!     write(*,*) 'TOO MANY STATIONS !!!!!'
!     stop
!     endif
 
      DO IR = 1, NR
         IF (NUSE(IR)==0) THEN
            WEIG(1, IR) = 0.
                       !!!!!!!! weights REDEFINED to ZERO if station not used
            WEIG(2, IR) = 0.
            WEIG(3, IR) = 0.
         ENDIF
      ENDDO
 
 
      NTIM = NUM_OF_TIME_SAMPLES
!        ntim=4000
!
!     READING SEISMOGRAMS (=DATA)
!
 
 
      DO IR = 1, NR
         !NFILE = 1000 + 1*IR
         NFILE = 2000 + 1*IR
         !STATFIL1 = TRIM(STATNAME(IR))//'fil.dat'
         REWIND(NFILE)
         !OPEN (NFILE, FILE = STATFIL1)
         DO ITIM = 1, NTIM
            READ (NFILE, '(4(1x,e12.6))') TIME, X1(ITIM, IR, 1), X1(ITIM, IR, 2), X1(ITIM, IR, 3)
         ENDDO
         CLOSE (NFILE)
      ENDDO
 
 
      DO IR = 1, NR
         !NFILE = 2000 + 1*IR
         NFILE = 3000 + 1*IR
         !STATFIL2 = TRIM(STATNAME(IR))//'res.dat'
         REWIND(NFILE)
         !OPEN (NFILE, FILE = STATFIL2)
         DO ITIM = 1, NTIM
            READ (NFILE, '(4(1x,e12.6))') TIME, X2(ITIM, IR, 1), X2(ITIM, IR, 2), X2(ITIM, IR, 3)
         ENDDO
         CLOSE (NFILE)
      ENDDO
 
      DO ICOM = 1, 3
         DO IR = 1, NR
            DO ITIM = 1, NTIM
               X3(ITIM, IR, ICOM) = X1(ITIM, IR, ICOM) - X2(ITIM, IR, ICOM)
            ENDDO
         ENDDO
      ENDDO
 
 
      ANORM1 = 0.
      DO ICOM = 1, 3
         DO IR = 1, NR
            DO ITIM = 1, NTIM
!        if(ir.ne.7.or.icom.ne.2)       !temporary only
!     *    		anorm1=anorm1 +x1(itim,ir,icom)*x1(itim,ir,icom)
!     	anorm1=anorm1+x1(itim,ir,icom)*x1(itim,ir,icom)
               ANORM1 = ANORM1 + X1(ITIM, IR, ICOM)*X1(ITIM, IR, ICOM)*WEIG(ICOM, IR)**2.
            ENDDO
         ENDDO
      ENDDO
      ANORM1 = ANORM1*DT      !sum of data squared
 
 
      ANORM2 = 0.
      DO ICOM = 1, 3
         DO IR = 1, NR
            DO ITIM = 1, NTIM
!    	anorm2=anorm2+x2(itim,ir,icom)*x2(itim,ir,icom)
               ANORM2 = ANORM2 + X2(ITIM, IR, ICOM)*X2(ITIM, IR, ICOM)*WEIG(ICOM, IR)**2.
            ENDDO
         ENDDO
      ENDDO
      ANORM2 = ANORM2*DT      !integral of residuals squared
 
 
      ANORM3 = 0.
      DO ICOM = 1, 3
         DO IR = 1, NR
            DO ITIM = 1, NTIM
               ANORM3 = ANORM3 + X3(ITIM, IR, ICOM)*X3(ITIM, IR, ICOM)
            ENDDO
         ENDDO
      ENDDO
      ANORM3 = ANORM3*DT
 
      VARRED = 1. - ANORM2/ANORM1       ! variance reduction WEIGHTED
 
      NRUSED = 0
      DO IR = 1, NR
         NRUSED = NRUSED + NUSE(IR)
                              !number of used stations
      ENDDO
 
      NUMBD = (3*NRUSED*5)
                         ! number of data
!	  (assuming 5 independent samples per comp, only!!!)
!     (assuming, for simplicity, that all 3 comp. are taken at each used station)
!     (assuming just 1 subevent)
      NUMBP = 7
               ! number of parameters ( 1 sub with 5 MT components + 1 position + 1 time)
!                ! if more subevents, then substitute 7 by 7 x number of subevents
!                ! unimportant if numbd >> numbp
      NUMBDF = NUMBD - NUMBP
 
!     postvar=anorm2/float(numbdf)  ! type0 (old; worked well, but not justified)
      DUREF = (1./F1)
!     postvar=anorm2/(3.*float(nrused)*(1.-varred)*duref) ! type2
!     postvar=anorm1/(3.*float(nrused)*duref) ! type3 = type2 = estimate of mean ampl.
      POSTVAR = ANORM2/(3.*FLOAT(NRUSED)*DUREF)
                                              ! type4 (often will be too low)
 
 
      DO IR = 1, NR
         !NFILE = 3000 + 1*IR
         NFILE = 4000 + 1*IR
         STATFIL3 = TRIM(STATNAME(IR))//'syn.dat'
         OPEN (NFILE, FILE = STATFIL3)
         DO ITIM = 1, NTIM
            TIME = FLOAT(ITIM - 1)*DT
            WRITE (NFILE, '(4(1x,e12.6))') TIME, X3(ITIM, IR, 1), X3(ITIM, IR, 2), X3(ITIM, IR, 3)
         ENDDO
         CLOSE (NFILE)
      ENDDO
 
 
      WRITE (*, *)
      WRITE (*, *) 'variance reduction, weighted'
      WRITE (*, *) 'varred=', VARRED
!      write(*,*) 'posterior data variance'
!      write(*,*) 'postvar=', postvar
 
 
      WRITE (101, *) ANORM2, ANORM1, VARRED, POSTVAR, NUMBD, NUMBP, NUMBDF
 
      WRITE (101, *)
      WRITE (101, *) 'integral of residuals squared, weighted'
      WRITE (101, *) ANORM2
 
      WRITE (101, *)
      WRITE (101, *) 'integral of data squared, , weighted'
      WRITE (101, *) ANORM1
 
      WRITE (101, *)
      WRITE (101, *) 'variance reduction, weighted'
      WRITE (101, *) VARRED
 
      WRITE (101, *)
      WRITE (101, *) 'number of data'
      WRITE (101, *) '(assuming 5 indep. samples per comp.)'
      WRITE (101, *) '(this is a problematic estimate)'
      WRITE (101, *) NUMBD
 
      WRITE (101, *)
      WRITE (101, *) 'number of parameters'
      WRITE (101, *) '(assuming 1 sub, 5 MT components,1 position,1 time)'
      WRITE (101, *) NUMBP
 
      WRITE (101, *)
      WRITE (101, *) 'number of degrees of freedom'
      WRITE (101, *) '(this is a problematic estimate)'
      WRITE (101, *) NUMBDF
 
      WRITE (101, *)
      WRITE (101, *) 'posterior data variance'
      WRITE (101, *) '(this is a problematic estimate)'
      WRITE (101, *) POSTVAR
 
 
      END SUBROUTINE NORM15_POSTVAR
 
 
