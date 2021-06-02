
      PROGRAM ISOLA15

      USE COMMON_VARS
      USE PARAMETERS

      IMPLICIT NONE

      INCLUDE "interface/interface_elemat15.inc"
      INCLUDE "interface/interface_oneinv15.inc"
      INCLUDE "interface/interface_silsub.inc"

      REAL ACCUMERR, AINB, AINBS, AINP, AINPS, AINT, AINTS, AMIN, AMOMENT, AMOMENTS, AMOPP, AMORP, AMORR, AMORT, AMOTP, AMOTT
      REAL AMOXX, AMOXY, AMOXZ, AMOYY, AMOYZ, AMOZZ, AOPT, AOPTSUM, ASAVE, AVOL, AZB, AZBS, AZP, AZPS, AZT, AZTS
      REAL CBEST, CBESTALL, CORR, CORRMAX, DCP, DCPERC, DCPERCS, DIP1, DIP1S, DIP2, DIP2S
      REAL EIGMAX, EIGMIN, EIGRAT, EN, ERROR2, FINR, GOLD, ORI, POSTVAR
      REAL RAKE1, RAKE1S, RAKE2, RAKE2S, RINV, ROLD, RRORI, STR1, STR1S, STR2, STR2S, SX, SYN
      REAL TIME, TWOAAAA, TWOCORR, TWODCPER, TWODIP, TWODIP2, TWOMOM, TWOMSFT, TWORAK, TWORAK2, TWOSHFT, TWOSTR, TWOSTR2, TWOVOLPER
      REAL TWOVVVV, VARDAT, VARRED, VOPT, VV, W, XINV, XMO, XMOMMAG, XMOMMAGS
      INTEGER I, IA, IBEGIN, IBEST, ICOM, ID1, ID2, IERR, IM, IOPTSHF, IR, IR1, IR2
      INTEGER IRECALL, IS1, IS2, ISELECT, ISEQ, ISEQM, ISEQUEN, ISH, ISOUR, ISOURMAX, ISUB, ITIM, JM
      INTEGER KEYGPS, KEYINV, N, NFILE, NROT, NUMF1, NUMF2, NUSE

!     POZOR - pokud se zmeni pocet bodu pri vypocet VR tak zmenit take manidata15.inc
!     !!!!!!!!!!!!!!!!!!!!

!     mixing GPS and seis stations
!     ALLSTAT: GPS pure statics (2), gps-grams (1), normal seismograms (0)

!     freq range is station dependent        !!!!!!!!!!!!!!
!     !!!!! use LF for GPS and other for seis

!     Multiple-point source inversion of  COMPLETE waveforms from
!     r e g i o n a l   or    l o c a l  stations.
!     Iterative deconvolution, similar to Kikuchi and Kanamori (1991).
!     Green's functions by discrete wavenumber method of Bouchon (1981).
!     Moment tensor of subevents is found by least-square minimization
!     of misfit between observed and synthetic waveforms, while position
!     and time of subevents is found by maximization of correlation
!     through grid search.
!     Possible modes of the moment tensor (MT) retrieval are as follows:
!     - full MT : DC+CLVD+VOL
!     - deviatoric MT  DC+CLVD; assuming VOL%=0 (=RECOMMENDED OPTION)
!     - constrained double-couple MT: assuming VOL%=CLVD%=0
!     - known and fixed 100% DC (only position and time is retrieved)

!     Author: J. Zahradnik 2003-2011-2015

!     DIMENSION   if you want to increase dim for w(...), change also:
!                 oneinv2, oneinv2_vol, subevnt, subevnt_vol, elemat2, eleonly
!                 adding more stations requires also change of filenumbers in open(...)

      CHARACTER*2 CHR(100)
      CHARACTER*5 STATNAME(MAX_STATIONS)
      CHARACTER*40 STATFIL1, STATFIL2, STATFIL3 !!g77
      CHARACTER*10 CORRFILE

      CHARACTER*40 allstatfile, inpinvfile, grdatfile, crustalfile, sourcefile, statfile, grhesfile, rawdir

      DIMENSION XINV(NUM_OF_TIME_SAMPLES, MAX_STATIONS, 3)
      DIMENSION ORI(NUM_OF_TIME_SAMPLES, MAX_STATIONS, 3)
      DIMENSION SX(NUM_OF_TIME_SAMPLES, MAX_STATIONS, 3)
      DIMENSION SYN(NUM_OF_TIME_SAMPLES, MAX_STATIONS, 3)
      DIMENSION FINR(NUM_OF_TIME_SAMPLES, MAX_STATIONS, 3)
      DIMENSION W(MIN_TIME_SHIFT:MAX_TIME_SHIFT, MAX_STATIONS, 3, 6, MAX_SOURCE_POSITIONS)    !  6 = full moment, 5 = deviatoric only
      DIMENSION ROLD(6, 6, MAX_SOURCE_POSITIONS)
      DIMENSION RINV(6, 6, MAX_SOURCE_POSITIONS)
      DIMENSION CORR(MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION ISH(MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION ASAVE(6, MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION AOPT(6)
      DIMENSION VOPT(6)
      DIMENSION AOPTSUM(6)
      DIMENSION TWOCORR(MAX_SOURCE_POSITIONS, MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION TWOSHFT(MAX_SOURCE_POSITIONS, MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION TWOSTR(MAX_SOURCE_POSITIONS, MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION TWODIP(MAX_SOURCE_POSITIONS, MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION TWORAK(MAX_SOURCE_POSITIONS, MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION TWOSTR2(MAX_SOURCE_POSITIONS, MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION TWODIP2(MAX_SOURCE_POSITIONS, MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION TWORAK2(MAX_SOURCE_POSITIONS, MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION TWODCPER(MAX_SOURCE_POSITIONS, MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION TWOVOLPER(MAX_SOURCE_POSITIONS, MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION TWOMSFT(MAX_SOURCE_POSITIONS, MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION TWOMOM(MAX_SOURCE_POSITIONS, MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION TWOAAAA(MAX_SOURCE_POSITIONS, MAX_NUM_OF_TIME_SHIFTS, 6)
      DIMENSION TWOVVVV(MAX_SOURCE_POSITIONS, MAX_NUM_OF_TIME_SHIFTS, 6)
      DIMENSION IBEST(MAX_SOURCE_POSITIONS)
      DIMENSION CBEST(MAX_SOURCE_POSITIONS)
      DIMENSION XMO(MAX_SOURCE_POSITIONS)
      DIMENSION DCP(MAX_SOURCE_POSITIONS)
      DIMENSION IS1(MAX_SOURCE_POSITIONS)
      DIMENSION ID1(MAX_SOURCE_POSITIONS)
      DIMENSION IR1(MAX_SOURCE_POSITIONS)
      DIMENSION IS2(MAX_SOURCE_POSITIONS)
      DIMENSION ID2(MAX_SOURCE_POSITIONS)
      DIMENSION IR2(MAX_SOURCE_POSITIONS)
      DIMENSION NUSE(MAX_STATIONS)
      DIMENSION KEYGPS(MAX_STATIONS)
      DIMENSION GOLD(6, 6), EN(6), VV(6, 6)

      INTEGER*8 RATE, START_TIME, END_TIME

      DATA CHR/'00', '01', '02', '03', '04', '05', '06', '07', '08', '09', &
               '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', &
               '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', &
               '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', &
               '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', &
               '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', &
               '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', &
               '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', &
               '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', &
               '90', '91', '92', '93', '94', '95', '96', '97', '98', '99'/

      INTEGER NC, NFREQ, NS, IKMAX
      REAL    TL, AW, XL, UCONV, FREF
      integer iargc

      !NAMELIST /INPUT/ NC, NFREQ, TL, AW, NR, NS, XL, IKMAX, UCONV, FREF
      NAMELIST /INPUT/ NFREQ, TL, AW, XL, IKMAX, UCONV, FREF

      WRITE (*, *)
      WRITE (*, *) 'This is ISOLA'
      WRITE (*, *)

      NTIM = NUM_OF_TIME_SAMPLES

      allstatfile="allstat.dat"
      inpinvfile="inpinv.dat"
      grdatfile="grdat.hed"
      grhesfile="gr.hes"
      crustalfile="crustal.dat" ! count nc
      sourcefile="source.dat" ! count ns
      statfile="station.dat" ! count nr
      rawdir="./" ! current dir

      if (iargc() .ge. 0) then
        call getarg(1, allstatfile)
        call getarg(2, inpinvfile)
        call getarg(3, grdatfile)
        call getarg(4, grhesfile)
        call getarg(5, crustalfile) ! count nc
        call getarg(6, sourcefile) ! count ns
        call getarg(7, statfile) ! count nr
        call getarg(8, rawdir)
      endif

      OPEN (892, FILE = allstatfile)                               ! input: stations, weights,...
      OPEN (893, FILE = inpinvfile)                                ! input: control parameters
      OPEN (898, FORM = 'formatted', FILE = grdatfile)             ! input: more config 

      OPEN (894, FILE = 'inv1.dat')                                  ! output: all details
      OPEN (895, FILE = 'inv2.dat')                                  !       : less details
      OPEN (896, FILE = 'inv2c.dat')                                 !       : less details, summed subevents
      OPEN (897, FILE = 'inv3.dat')                                  !       : moment tensor (for GMT)

    !  READ (893, *)
    !  READ (893, *) KEYINV
    !  READ (893, *)
    !  READ (893, *) DT
    !  READ (893, *)
    !  READ (893, *) ISOURMAX
    !  READ (893, *)
    !  READ (893, *)
      READ (893, *) IBEGIN, ISTEP, ILAST
      IFIRST = IBEGIN - ISTEP
    !  READ (893, *)
    !  READ (893, *) ISUBMAX
    !  READ (893, *)
    !  READ (893, *)
    !  READ (893, *) F1, F2, F3, F4          ! not used (read from allstat)
    !  READ (893, *)
    !  READ (893, *) VARDAT
      CLOSE(893)

      READ (898, INPUT)
      CLOSE(898)

      ! read CRUSTAL
      ! crustal format: 2 comment lines + the rest layers
      open (110, form='formatted', file=crustalfile)
      ! count number of crustals
      NC=0
      do while(.true.)
        read(110,*,end=91)
        NC=NC+1
      end do
91 close(110)
      NC=NC-2 ! throw 2 comment lines away

      ! read SOURCE
      open (110, form='formatted', file=sourcefile)
      ! count number of sources
      NS=0
      do while(.true.)
        read(110,*,end=92)
        NS=NS+1
      end do
92 close(110)

      ! read STATION
      open (110, form='formatted', file=statfile)
      ! count number of stations
      NR=0
      do while(.true.)
        read(110,*,end=93)
        NR=NR+1
      end do
93 close(110)

      KEYINV=2 !fixed to deviatoric type
      DT=TL/8192.
      ISOURMAX=NS
      ISUBMAX=1
      VARDAT=2.0e-12


      IF (ISOURMAX > MAX_SOURCE_POSITIONS) THEN
         WRITE (*, *) 'Number of source positions (', ISOURMAX, ') exceeds maximum allowed (', MAX_SOURCE_POSITIONS, ').'
         STOP
      ENDIF

      IF (NS > MAX_SOURCE_POSITIONS) THEN
         WRITE (*, *) 'Number of source positions (', NS, ') exceeds maximum allowed (', MAX_SOURCE_POSITIONS, ').'
         STOP
      ENDIF

      IF (NR > MAX_STATIONS) THEN
         WRITE (*, *) 'Number of stations (', NR, ') exceeds maximum allowed (', MAX_STATIONS, ').'
         STOP
      ENDIF

      IF (NC > MAX_CRUSTAL_LAYERS) THEN
         WRITE (*, *) 'Number of crustal layers (', NC, ') exceeds maximum allowed (', MAX_CRUSTAL_LAYERS, ').'
         STOP
      ENDIF

      CALL SYSTEM_CLOCK(COUNT_RATE = RATE)

      WRITE (*, *) 'Initializing array w'

      CALL SYSTEM_CLOCK(START_TIME)

      W = 0.

      CALL SYSTEM_CLOCK(END_TIME)

      WRITE (*, *) 'Time required: ', REAL(END_TIME - START_TIME)/REAL(RATE)
      WRITE (*, *)

!*********************************************************************
!******************** LOOP OVER SOURCE POSITIONS *********************
!*********************************************************************
      WRITE (*, *) 'Reading input file gr.hes'

      CALL SYSTEM_CLOCK(START_TIME)

      OPEN(898, FORM = 'unformatted', FILE = grhesfile)

      DO ISOUR = 1, ISOURMAX
         CALL ELEMSE(W, NTIM, TL, AW, NFREQ, ISOUR, 898)
      ENDDO

      CLOSE(898)

      CALL SYSTEM_CLOCK(END_TIME)

!$ACC ENTER DATA COPYIN(W) ASYNC

      WRITE (*, *) 'Time required: ', REAL(END_TIME - START_TIME)/REAL(RATE)
      WRITE (*, *)

      AOPT(6) = 0.
      VOPT(6) = 0.
      AOPTSUM(6) = 0.

      DO IR = 1, NR
                                  ! nuse=0 station not used
                                  ! weig=0 component not used
                                  ! keygps not read, but defined below as FIXED
         READ (892, *) STATNAME(IR), NUSE(IR), WEIG(IR, 1), WEIG(IR, 2), WEIG(IR, 3), FF1(IR), FF2(IR), FF3(IR), FF4(IR)

         STAT(IR) = .TRUE.
         IF (NUSE(IR)==0) STAT(IR) = .FALSE.
         NTM(IR) = NUM_OF_TIME_SAMPLES + 1
      ENDDO

      WRITE (*, *) 'number of stations (max ', MAX_STATIONS, '), nr=', NR
      WRITE (*, *) 'stations used in the inversion:'
      DO IR = 1, NR
         IF (STAT(IR)) WRITE (*, *) STATNAME(IR)
      ENDDO
      CLOSE (892)

      DO IR = 1, NR          !!!!!!!!!1.2.2015
         KEYGPS(IR) = 0.     ! seismic stations
      ENDDO

      DO IR = 1, NR
         NUMF1 = 1000 + (IR*1)
         NUMF2 = 2000 + (IR*1)
         STATFIL1 = TRIM(rawdir)//'/'//TRIM(STATNAME(IR))//'raw.dat'
         STATFIL2 = TRIM(STATNAME(IR))//'fil.dat'
         write(*,*) STATFIL1
         write(*,*) STATFIL2
         OPEN (NUMF1, FILE = STATFIL1)
         OPEN (NUMF2, FILE = STATFIL2)
      ENDDO

      NMOM = 6
      IF (KEYINV/=1) NMOM = 5

     ! IF (IBEGIN<MIN_TIME_SHIFT .OR. ILAST>IABS(MIN_TIME_SHIFT)) THEN
      IF (IBEGIN<MIN_TIME_SHIFT .OR. ILAST>IABS(MAX_TIME_SHIFT)) THEN
         WRITE (*, *) 'limit of time shift exceeded; check ifirst, ilast'
         STOP
      ENDIF
      ISEQM = (ILAST - IFIRST)/ISTEP        ! number of tested time shifts
      IF (ISEQM>MAX_NUM_OF_TIME_SHIFTS) THEN
         WRITE (*, *) 'too many shifts requested; check ifirst,ilast,istep'
         STOP
      ENDIF

!*********************************************************************
!********************** MANIPULATE OBSERVED DATA *********************
!*********************************************************************

!      input: read from file in the subroutine
!     output: ori=ori data filtered
!             rrori=data power

      CALL MANIDATA15(KEYGPS, ORI, RRORI)   ! manipulate OBSERVED data

!*********************************************************************
!******                  LOOP OVER SUBEVENTS                   *******
!****** For each one, ALL possible source positions are tested *******
!*********************************************************************

      WRITE (*, *)

      DO ICOM = 1, 3
         DO IR = 1, NR
            DO ITIM = 1, NTIM
               SYN(ITIM, IR, ICOM) = 0.     ! initialization
            ENDDO
         ENDDO
      ENDDO


      ISUB = 0                    ! counting subevents
 200  ISUB = ISUB + 1

      WRITE (*, *)
      WRITE (*, *) 'searching subevent #', ISUB
      WRITE (*, *)

      IF (ISUB==1) THEN
         DO ICOM = 1, 3
            DO IR = 1, NR
               DO ITIM = 1, NTIM
                  XINV(ITIM, IR, ICOM) = ORI(ITIM, IR, ICOM)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO ICOM = 1, 3
            DO IR = 1, NR
               DO ITIM = 1, NTIM
                  XINV(ITIM, IR, ICOM) = FINR(ITIM, IR, ICOM)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!*************************************************************************
!**** FITTING DATA WITH A GIVEN SOURCE POSITION AND VARYING TIME SHIFT ***
!*************************************************************************

      WRITE (*, *) 'Starting calculations for multiple source points'

      CALL SYSTEM_CLOCK(START_TIME)

!$ACC WAIT

!$ACC UPDATE DEVICE(WEIG,FF1,FF2,FF3,FF4,ANDC,A2,B,DT,NTM,NR,NTIM,NMOM,ISUBMAX,IFIRST,ISTEP,ILAST,STAT)
!$ACC PARALLEL LOOP GANG WORKER INDEPENDENT                                                                     &
!$ACC COPYIN(ROLD,RINV,KEYINV,KEYGPS,XINV,ISEQM,RRORI,VARDAT,ISOURMAX)                                                   &
!$ACC COPY(TWOCORR,TWOMSFT,TWOSHFT,TWOAAAA,TWOVVVV,TWOSTR,TWODIP,TWORAK,TWOSTR2,TWODIP2,TWORAK2)                &
!$ACC COPY(TWODCPER,TWOVOLPER,TWOMOM,CBEST,IBEST,IR1,IR2,ID1,ID2,IS1,IS2,DCP,XMO)                               &
!$ACC PRIVATE(ISOUR,ASAVE,CORR,ISH,RAKE1,RAKE2,DIP1,DIP2,STR1,STR2,CORRMAX,IRECALL,IOPTSHF,AVOL,DCPERC,AMOMENT) &
!$ACC FIRSTPRIVATE(AOPT)

!$OMP PARALLEL DO DEFAULT(NONE)                                                                                 &
!$OMP SHARED(W,ROLD,RINV,KEYINV,KEYGPS,XINV,ISEQM,RRORI,DT,NMOM,VARDAT,ISOURMAX)                                         &
!$OMP SHARED(TWOCORR,TWOMSFT,TWOSHFT,TWOAAAA,TWOVVVV,TWOSTR,TWODIP,TWORAK,TWOSTR2,TWODIP2,TWORAK2)              &
!$OMP SHARED(TWODCPER,TWOVOLPER,TWOMOM,CBEST,IBEST,IR1,IR2,ID1,ID2,IS1,IS2,DCP,XMO)                             &
!$OMP PRIVATE(ISOUR,ASAVE,CORR,ISH,RAKE1,RAKE2,DIP1,DIP2,STR1,STR2,CORRMAX,IRECALL,IOPTSHF,AVOL,DCPERC,AMOMENT) &
!$OMP FIRSTPRIVATE(AOPT)

!!!!!!$OMP TARGET UPDATE TO(WEIG,FF1,FF2,FF3,FF4,ANDC,A2,B,DT,NTM,NR,NTIM,NMOM,ISUBMAX,IFIRST,ISTEP,ILAST,STAT)
!!!!!!$OMP TARGET MAP(TO:W,ROLD,RINV,KEYINV,KEYGPS,XINV,ISEQM,RRORI,VARDAT,ISOURMAX)                                          &
!!!!!!$OMP MAP(TOFROM:TWOCORR,TWOMSFT,TWOSHFT,TWOAAAA,TWOVVVV,TWOSTR,TWODIP,TWORAK,TWOSTR2,TWODIP2,TWORAK2)          &
!!!!!!$OMP MAP(TOFROM:TWODCPER,TWOVOLPER,TWOMOM,CBEST,IBEST,IR1,IR2,ID1,ID2,IS1,IS2,DCP,XMO)
!!!!!!$OMP TEAMS DISTRIBUTE PARALLEL DO                                                                              &
!!!!!!$OMP SHARED(WEIG,FF1,FF2,FF3,FF4,ANDC,A2,B,DT,NTM,NR,NTIM,NMOM,ISUBMAX,IFIRST,ISTEP,ILAST,STAT)                &
!!!!!!$OMP SHARED(W,ROLD,RINV,KEYINV,KEYGPS,XINV,ISEQM,RRORI,VARDAT)                                                 &
!!!!!!$OMP SHARED(TWOCORR,TWOMSFT,TWOSHFT,TWOAAAA,TWOVVVV,TWOSTR,TWODIP,TWORAK,TWOSTR2,TWODIP2,TWORAK2)              &
!!!!!!$OMP SHARED(TWODCPER,TWOVOLPER,TWOMOM,CBEST,IBEST,IR1,IR2,ID1,ID2,IS1,IS2,DCP,XMO)                             &
!!!!!!$OMP PRIVATE(ISOUR,ASAVE,CORR,ISH,RAKE1,RAKE2,DIP1,DIP2,STR1,STR2,CORRMAX,IRECALL,IOPTSHF,AVOL,DCPERC,AMOMENT) &
!!!!!!$OMP FIRSTPRIVATE(AOPT)
      DO ISOUR = 1, ISOURMAX
!       !!! ATTENTION: rold comes * 1.e20; rinv comes / 1.e20 !!!
!               input: w=elem seismograms for a given source position
!              output: rold=data matrix
!                      rinv=inverse matrix

         CALL ELEMAT15(KEYGPS, W, ROLD, RINV, ISOUR)

         IF (KEYINV==1) CALL ONEINV15(XINV, W, ROLD, RINV, ASAVE, CORR, ISH, ISOUR)     ! nmom=6
         IF (KEYINV==2) CALL ONEINV15(XINV, W, ROLD, RINV, ASAVE, CORR, ISH, ISOUR)     ! nmom=5

!        action: 'one inversion'; fitting xinv data by a single subevent
!                right hand side formed from xinv data and elemse (green)
!                elemse data for each source used repeatedly with several time shifts
!                1 source position = 1 Green (elemse w); 1 RINV inv. matrix
!                each time shift of w = its own right-hand side and solution
!         input: xinv=data to be inverted for a set of time shifts
!                w =elem seis
!                rold =system matrix
!                rinv =inverse  matrix
!        output: asave=moment tensor coefficients (array for all shifts)
!                corr= correlation (array for all shifts)
!                ish=shifts (array of integer shifts for all shift steps)

         DO I = 1, ISEQM                                  ! saving correlation and shift values
            TWOCORR(ISOUR, I) = CORR(I)
            TWOMSFT(ISOUR, I) = (1. - CORR(I)**2)*RRORI   !!!!!!new march2012
            TWOSHFT(ISOUR, I) = FLOAT(ISH(I))*DT
            DO N = 1, NMOM
               TWOAAAA(ISOUR, I, N) = ASAVE(N, I)
               TWOVVVV(ISOUR, I, N) = SQRT(VARDAT*RINV(N, N, ISOUR)*1.E20)    !  sigma =sqrt(var); corrected for formal 1.e20
               AOPT(N) = ASAVE(N, I)
            ENDDO
            CALL SILSUB(AOPT, STR1, DIP1, RAKE1, STR2, DIP2, RAKE2, AMOMENT, DCPERC, AVOL)

            TWOSTR(ISOUR, I) = STR1
            TWODIP(ISOUR, I) = DIP1
            TWORAK(ISOUR, I) = RAKE1
            TWOSTR2(ISOUR, I) = STR2
            TWODIP2(ISOUR, I) = DIP2
            TWORAK2(ISOUR, I) = RAKE2
            TWODCPER(ISOUR, I) = DCPERC
            TWOVOLPER(ISOUR, I) = AVOL                    !!!!!new march2012
            TWOMOM(ISOUR, I) = AMOMENT
         ENDDO

!******************************************************************************
!*** INSPECTING RESULTS OF THE SHIFT LOOP (SEARCHING SHIFT WITH OPT. CORR.) ***
!***    SAVING THOSE OF THE BEST CORRELATION FOR A GIVEN SOURCE POSITION    ***
!******************************************************************************

         CORRMAX = -100.                   ! 0. changed to -100: 31.12.2013
         DO I = 1, ISEQM
            IF (CORR(I)>CORRMAX) THEN
               IRECALL = I                 ! sequential number of the best shift (1,2,...iseqm)
               IOPTSHF = ISH(I)            ! best value of 'ishift' (in time steps, not in sec)
               CORRMAX = CORR(I)           ! best value of correlation
            ENDIF
         ENDDO
         DO N = 1, 6
            AOPT(N) = ASAVE(N, IRECALL)    ! a's for optimum shift
         ENDDO

         ! saving results for the best correlation
         ! (for each source position)

         CALL SILSUB(AOPT, STR1, DIP1, RAKE1, STR2, DIP2, RAKE2, AMOMENT, DCPERC, AVOL)
         XMO(ISOUR) = AMOMENT
         DCP(ISOUR) = DCPERC               ! zde mozno DOCASNE mit dcp(isour)=avol
         IS1(ISOUR) = IFIX(STR1)
         ID1(ISOUR) = IFIX(DIP1)
         IR1(ISOUR) = IFIX(RAKE1)
         IS2(ISOUR) = IFIX(STR2)
         ID2(ISOUR) = IFIX(DIP2)
         IR2(ISOUR) = IFIX(RAKE2)
         IBEST(ISOUR) = IOPTSHF
         CBEST(ISOUR) = CORRMAX
      ENDDO        ! End of ISOUR loop
!!!!!!$OMP END TEAMS DISTRIBUTE PARALLEL DO
!!!!!!$OMP END TARGET
!!!!!!$OMP END PARALLEL DO

!$ACC UPDATE SELF(ERROR_ON_DEVICE)

      CALL SYSTEM_CLOCK(END_TIME)

      WRITE (*, *) 'Time required: ', REAL(END_TIME - START_TIME)/REAL(RATE)
      WRITE (*, *)

      IF (ERROR_ON_DEVICE /= 0) THEN
         IF (IBITS(ERROR_ON_DEVICE, 0, 1) == 1) THEN
            WRITE (*, *) 'Too many iterations in JACOBINR for at least one thread.'
         ENDIF
         IF (IBITS(ERROR_ON_DEVICE, 1, 1) == 1) THEN
            WRITE (*, *) 'Singular matrix found in LUDCMP for at least one thread.'
         ENDIF
         STOP
      ENDIF

! ==============================================================
      CORRFILE = 'corr'//CHR(ISUB)//'.dat'
      OPEN (899, FILE = CORRFILE)
      WRITE (899, *)
      WRITE (899, *) '2D correlation for isub=', ISUB
      DO ISOUR = 1, ISOURMAX
         DO ISEQ = 1, ISEQM                      !!! new march2012
                                                 !!! new march2012
            WRITE (899, '(1x,i5,2(1x,f9.4),5x,            3(1x,i5),5x,3(1x,i5),2(1x,f7.2),2(1x,e12.6))') ISOUR,                &
                   TWOSHFT(ISOUR, ISEQ), TWOCORR(ISOUR, ISEQ), IFIX(TWOSTR(ISOUR, ISEQ)), &
                   IFIX(TWODIP(ISOUR, ISEQ)), IFIX(TWORAK(ISOUR, ISEQ)),      &
                   IFIX(TWOSTR2(ISOUR, ISEQ)), IFIX(TWODIP2(ISOUR, ISEQ)), IFIX(TWORAK2(ISOUR, ISEQ)), &
                   TWODCPER(ISOUR, ISEQ), TWOVOLPER(ISOUR, ISEQ),&
                   TWOMSFT(ISOUR, ISEQ), TWOMOM(ISOUR, ISEQ)
                                                 !!! new march2012
         ENDDO
      ENDDO
      CLOSE (899)

      WRITE (894, *) 'All trial positions and shifts for subevent#', ISUB
      WRITE (894, *) '(isour = source position,ishift*dt=time shift)'
      WRITE (894, *) 'isour,ishift,corr,moment,DC%,str,dip,rak,str,dip,rak'
      DO ISOUR = 1, ISOURMAX
         WRITE (894, '(2x,i4,2x,i5,2x,f9.4,2x,e15.4,2x,f8.3,2x,i5,2x,i5,2x,i5,2x,i5,2x,i5,2x,i5)') ISOUR, IBEST(ISOUR),           &
                CBEST(ISOUR), XMO(ISOUR), DCP(ISOUR), IS1(ISOUR), ID1(ISOUR), IR1(ISOUR), IS2(ISOUR), ID2(ISOUR), IR2(ISOUR)
      ENDDO

!*********************************************************************
!**************** SEARCHING THE BEST SOURCE POSITION *****************
!*********************************************************************

      CBESTALL = -100.                      ! changed 31.12.2013
      DO ISOUR = 1, ISOURMAX                ! AUTOMATIC
         IF (CBEST(ISOUR)>CBESTALL) THEN
            CBESTALL = CBEST(ISOUR)
            ISELECT = ISOUR                 ! optimum position
         ENDIF
      ENDDO

      IOPTSHF = IBEST(ISELECT)      ! optimum time shift

!
!     A possibility to MANUALLY change the selected source position 'iselect'
!
      WRITE (*, *) 'Trial source positions and shifts for subevent #', ISUB
      WRITE (*, *) '(position #, shift (multiples of dt), correl., DC%)'
      WRITE (*, *) '(strike1, dip1, rake 1   and  strike2, dip2, rake2)'

      DO ISOUR = 1, ISOURMAX
         WRITE (*, '(i4,2x,i4,2x,f7.4,2x,f6.2,2x,6i6)') ISOUR, IBEST(ISOUR), CBEST(ISOUR), DCP(ISOUR), IS1(ISOUR), ID1(ISOUR),     &
                IR1(ISOUR), IS2(ISOUR), ID2(ISOUR), IR2(ISOUR)
      ENDDO
      WRITE (*, *) 'automatic search suggests trial source #', ISELECT

      ISEQUEN = (IOPTSHF - IFIRST)/ISTEP

      DO N = 1, NMOM
         AOPT(N) = TWOAAAA(ISELECT, ISEQUEN, N)
         VOPT(N) = TWOVVVV(ISELECT, ISEQUEN, N)
      ENDDO

      AMOXX = -1.*AOPT(4) + AOPT(6)
      AMOYY = -1.*AOPT(5) + AOPT(6)
      AMOZZ = AOPT(4) + AOPT(5) + AOPT(6)
      AMOXY = AOPT(1)
      AMOXZ = AOPT(2)
      AMOYZ = -1.*AOPT(3)

      AMOTT = AMOXX     !t=theta, p=phi (delta), r=r
      AMOPP = AMOYY
      AMORR = AMOZZ
      AMOTP = -1.*AMOXY
      AMORT = AMOXZ
      AMORP = -1.*AMOYZ

      WRITE (*, *) 'Results for subevent #', ISUB
      WRITE (*, *) 'trial source position #', ISELECT
      WRITE (*, *) 'time shift (multiple of dt)', IOPTSHF
      CALL SILSUB(AOPT, STR1, DIP1, RAKE1, STR2, DIP2, RAKE2, AMOMENT, DCPERC, AVOL)
      WRITE (*, *) 'strike,dip,rake', IFIX(STR1), IFIX(DIP1), IFIX(RAKE1)
      WRITE (*, *) 'strike,dip,rake', IFIX(STR2), IFIX(DIP2), IFIX(RAKE2)

!        ---------start of -----NEW-----------------        Feb 2011

! we are already out of source loop , so we have to calculate w, sx

      OPEN(898, FORM = 'unformatted', FILE = grhesfile)

      ! For each source point, gr_xyz writes:
      !  - 3 records, each containing 6 COMPLEX*16 values (and two record markers, each 4 bytes, added by Fortran) : 3 * (6 * 16 + 8) bytes
      !  - This is repeated for each of the NFREQ frequencies and the NR stations
      CALL FSEEK(898, (ISELECT - 1) * NFREQ * NR * 3 * (6 * 16 + 8), 0)

      CALL ELEMSE(W, NTIM, TL, AW, NFREQ, ISELECT, 898)

      CLOSE(898)

      CALL ELEMAT15(KEYGPS, W, ROLD, RINV, ISELECT)   ! analyzing system matrix for the best position (independent of time shift)
                                                      ! iselect and w correspond to the best position

      DO IM = 1, NMOM
         DO JM = 1, NMOM
            GOLD(IM, JM) = ROLD(IM, JM, ISELECT)      ! rold comes multiplied by 1e20
         ENDDO
      ENDDO

      CALL JACOBINR(GOLD, NMOM, 6, EN, VV, NROT)      ! 6 is dimension, cannot be nmom
      DO IM = 1, NMOM                                 ! converting eigenvalues of GTG into sing. values of G
                                                      ! introducing (a constant) data variance, vardat, here is the same as to divide GTG by vardat
         EN(IM) = SQRT(EN(IM)/VARDAT)/1.E10           ! sqrt to get SING values; corrected for previous formal 1.e20
      ENDDO

      EIGMIN = 1.E30
      EIGMAX = 1.E-30
      DO IM = 1, NMOM
         IF (ABS(EN(IM))<EIGMIN) EIGMIN = ABS(EN(IM))
         IF (ABS(EN(IM))>EIGMAX) EIGMAX = ABS(EN(IM))
      ENDDO
      EIGRAT = EIGMIN/EIGMAX
      WRITE (*, *)
      WRITE (*, *) 'SING. values, inlcuding vardat (min., max, min/max):'
      WRITE (*, *) EIGMIN, EIGMAX, EIGRAT

!        ------------end of -------NEW---------------
      WRITE (*, *) AOPT

      CALL SILSUB(AOPT, STR1, DIP1, RAKE1, STR2, DIP2, RAKE2, AMOMENT, DCPERC, AVOL)

      CALL PL2PT(STR1, DIP1, RAKE1, AZP, AINP, AZT, AINT, AZB, AINB, IERR)               ! JZ Feb26,2011

      WRITE (894, *)
      WRITE (894, *) 'Selected source position for subevent #', ISUB
      WRITE (894, *) 'isour,ishift', ISELECT, IOPTSHF
      WRITE (894, *)
      WRITE (894, *) 'SINGULAR values, incl. vardat (min., max, max/min)'
      WRITE (894, *) EIGMIN, EIGMAX, 1./EIGRAT
      WRITE (894, *)

      WRITE (894, *) 'Inversion result:'
      WRITE (894, *) 'coefficients of elem.seismograms a(1),a(2),...a(6):'
      WRITE (894, '(6(1x,e12.6))') (AOPT(N), N = 1, 6)
      WRITE (894, *) 'and their sigma_a(1), sigma_a(2),... sigma_a(6):'
      WRITE (894, '(6(1x,e12.6))') (VOPT(N), N = 1, 6)

      WRITE (894, *)
      WRITE (894, *) 'moment (Nm):', AMOMENT
      XMOMMAG = (2.0/3.0)*LOG10(AMOMENT) - 6.0333
       ! Hanks&Kanamori(1979) !thimios

      WRITE (894, '(a20,f6.1)') 'moment magnitude:', XMOMMAG
      WRITE (894, '(a20,f6.1)') 'VOL % :', AVOL
      WRITE (894, '(a20,f6.1)') 'DC % :', DCPERC
      WRITE (894, '(a20,f6.1)') 'abs(CLVD) % :', 100. - ABS(AVOL) - DCPERC               ! instead avol should be abs(avol), see silsub.inc
      WRITE (894, *) 'strike,dip,rake:', IFIX(STR1), IFIX(DIP1), IFIX(RAKE1)
      WRITE (894, *) 'strike,dip,rake:', IFIX(STR2), IFIX(DIP2), IFIX(RAKE2)
      WRITE (894, *) 'P-axis azimuth and plunge:', IFIX(AZP), IFIX(AINP)
      WRITE (894, *) 'T-axis azimuth and plunge:', IFIX(AZT), IFIX(AINT)
      WRITE (894, *) 'B-axis azimuth and plunge:', IFIX(AZB), IFIX(AINB)
      WRITE (894, *)

!*********************************************************************
!**** subevent seismo for the best soure position and time shift *****
!*********************************************************************

      CALL SUBEVNT15(AOPT, W, IOPTSHF, SX, ISELECT)    ! seismo sx for optimum position and for optimum time shift ioptshf

!*********************************************************************
!********** SYNTH SEIMOGRAM BEING SUMNED UP FROM SUBEVENTS ***********
!*********************************************************************

      DO ICOM = 1, 3
         DO IR = 1, NR
            DO ITIM = 1, NTIM
               SYN(ITIM, IR, ICOM) = SYN(ITIM, IR, ICOM) + SX(ITIM, IR, ICOM)
            ENDDO
         ENDDO
      ENDDO

      DO ICOM = 1, 3
         DO IR = 1, NR
            DO ITIM = 1, NTIM
               FINR(ITIM, IR, ICOM) = ORI(ITIM, IR, ICOM) - SYN(ITIM, IR, ICOM)
            ENDDO
         ENDDO
      ENDDO

!
!     COMPUTING L2 MISFIT (error2) DIRECTLY = computing L2 norm of final resi
! 
      ERROR2 = 0.
      DO ICOM = 1, 3
         DO IR = 1, NR
            IF (STAT(IR)) THEN
               DO ITIM = 1, NTIM                                               ! if changing ntim to 4000 here, must change in manidata15 to get rrori
                  ERROR2 = ERROR2 + (FINR(ITIM, IR, ICOM)*WEIG(IR, ICOM))**2   ! new 9.9.2015 WEIGHTS included in VR
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      ERROR2 = ERROR2*DT

      POSTVAR = ERROR2/FLOAT(3*NR*NTIM)     ! a posteriori data variance
                                            ! good meaning only if all nr stations and points ntim were used

      ACCUMERR = ERROR2/RRORI               ! normalization by power of ORI data
                                            ! rrori (from manidata.inc) includes weights ! new 9.9. 2015
                                            ! (weights considered = VR refers to used components only)
      VARRED = 1. - ACCUMERR                ! variance reduction
                                            ! accumerr and varred consider only the used components
                                            ! (both in error2 and in rrori, assuming ntim=NUM_OF_TIME_SAMPLES)
                                            ! for other measures of the misfit, see NORM.FOR

      WRITE (894, *) 'After subtraction of subevent #', ISUB
      WRITE (894, *) 'weighted variance reduction (used components only):'
      WRITE (894, *) 'varred=', VARRED
      WRITE (894, *) '======================================='
      WRITE (894, *)

      CALL PL2PT(STR1, DIP1, RAKE1, AZP, AINP, AZT, AINT, AZB, AINB, IERR)

      WRITE (895, '(1x,i5,1x,f6.2,1x,e12.6,12(1x,f6.0),1x,f6.1,1x,e12.4)') ISELECT, IOPTSHF*DT, AMOMENT, STR1, DIP1, RAKE1, STR2,  &
             DIP2, RAKE2, AZP, AINP, AZT, AINT, AZB, AINB, DCPERC, VARRED

      AMIN = MIN(ABS(AMORR), ABS(AMOTT), ABS(AMOPP), ABS(AMORT), ABS(AMORP), ABS(AMOTP))
      IA = IFIX(LOG10(AMIN))
      AMORR = AMORR/10**FLOAT(IA)
      AMOTT = AMOTT/10**FLOAT(IA)
      AMOPP = AMOPP/10**FLOAT(IA)
      AMORT = AMORT/10**FLOAT(IA)
      AMORP = AMORP/10**FLOAT(IA)
      AMOTP = AMOTP/10**FLOAT(IA)
      WRITE (897, 99001) ISELECT, IOPTSHF, AMORR, IA, AMOTT, IA, AMOPP, IA, AMORT, IA, AMORP, IA, AMOTP, IA
99001 FORMAT (2(1x, i5), 3x, 6(f20.4, 'e+', i2.2))

      DO I = 1, NMOM
         AOPTSUM(I) = AOPTSUM(I) + AOPT(I)    ! summary MT
      ENDDO

      CALL SILSUB(AOPTSUM, STR1S, DIP1S, RAKE1S, STR2S, DIP2S, RAKE2S, AMOMENTS, DCPERCS, AVOL)

      CALL PL2PT(STR1S, DIP1S, RAKE1S, AZPS, AINPS, AZTS, AINTS, AZBS, AINBS, IERR)

      XMOMMAGS = (2.0/3.0)*LOG10(AMOMENTS) - 6.0333        ! summary Mw ! Jiri Apr10, 2011
      WRITE (896, '(1x,i5,1x,f6.2,1x,e12.6,1x,f6.2,12(1x,f6.0),1x,f6.1,1x,e12.4)') ISELECT, IOPTSHF*DT, AMOMENTS, XMOMMAGS, STR1S, &
             DIP1S, RAKE1S, STR2S, DIP2S, RAKE2S, AZPS, AINPS, AZTS, AINTS, AZBS, AINBS, DCPERCS, VARRED

      IF (ISUB<ISUBMAX) GOTO 200

      CLOSE (894)
      CLOSE (895)
      CLOSE (896)
      CLOSE (897)

!****************************************************************************
!* SAVING FINAL RESIDUAL SEISMOGRAM (final = after subtract all subevents) **
!****************************************************************************

      DO IR = 1, NR
         NFILE = 3000 + 1*IR
         STATFIL3 = TRIM(STATNAME(IR))//'res.dat'
         OPEN (NFILE, FILE = STATFIL3)
         DO ITIM = 1, NTIM
            TIME = FLOAT(ITIM - 1)*DT
            WRITE (NFILE, '(4(1x,e12.6))') TIME, FINR(ITIM, IR, 1), FINR(ITIM, IR, 2), FINR(ITIM, IR, 3)
         ENDDO
         !CLOSE (NFILE)
      ENDDO

      WRITE (*, *)
      WRITE (*, *) 'The following files were created:'
      WRITE (*, *) 'INV1.DAT, INV2.DAT, INV2c.DAT, INV3.DAT'
      WRITE (*, *) 'CORR01.DAT, ... (for all subevents)'
      WRITE (*, *) '*FIL.DAT, *RES.DAT, ... (for all stations)'

      CALL NORM15_POSTVAR(allstatfile, FF1(1), DT) ! assuming here that the first station contains the lowest freq (F1)
      CALL DSRETC(DIP1, STR1, RAKE1)
      END PROGRAM ISOLA15

! =================================================================

      INCLUDE "manidata15.inc"
      INCLUDE "elemat15.inc"
      INCLUDE "oneinv15.inc"
      INCLUDE "subevnt15.inc"
      INCLUDE "filter15.inc"
      INCLUDE "silsub.inc"
      INCLUDE "jacobi.inc"
      INCLUDE "line.inc"
      INCLUDE "ang.inc"
      INCLUDE "angles.inc"
      INCLUDE "lubksb.inc"
      INCLUDE "ludcmp.inc"
      INCLUDE "jacobinr.inc"
      INCLUDE "pl2pt.inc"
      INCLUDE "elemse.inc"
      INCLUDE "fft2cd.inc"
      INCLUDE "fsource.inc"
      INCLUDE "norm15.inc"
      INCLUDE "dsretc15.inc"

