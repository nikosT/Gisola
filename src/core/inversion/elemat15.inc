      SUBROUTINE ELEMAT15(KEYGPS, W, ROLD, RINV, ISOUR)

      USE COMMON_VARS
      USE PARAMETERS

      IMPLICIT NONE

      INCLUDE "interface/interface_filter15.inc"
      INCLUDE "interface/interface_jacobinr.inc"
      INCLUDE "interface/interface_lubksb.inc"
      INCLUDE "interface/interface_ludcmp.inc"

!$OMP DECLARE TARGET
!$ACC ROUTINE SEQ
      REAL DDD, EN, F1, F4, GOLD, R, RINV, ROLD, SHIFT, VV, W, Y
      INTEGER I, ICOM, INDX, IR, ISOUR, IT, ITIM, J, KEYDIS, KEYFIL, KEYGPS, M, N, NMAX
      INTEGER NROT

! ATTENTION: on output, rold is * 1.e20, rinv is / 1.e20

!   Input: elem se (w) for a fixed point source
!   Output: system matrix, inverted matrix
 
      DIMENSION KEYGPS(MAX_STATIONS)
      DIMENSION W(MIN_TIME_SHIFT:MAX_TIME_SHIFT, MAX_STATIONS, 3, 6, MAX_SOURCE_POSITIONS)
      DIMENSION R(6, 6)
      DIMENSION RINV(6, 6, MAX_SOURCE_POSITIONS)
      DIMENSION ROLD(6, 6, MAX_SOURCE_POSITIONS)
      DIMENSION GOLD(6, 6), EN(6), VV(6, 6)
      DIMENSION INDX(6), Y(6, 6)
 
      DOUBLE PRECISION WNWEIG, WMWEIG
                                     !!!!!! EXTREMELY IMPORTANT for gfort!!!!!!!
      REAL*4 MULTIP

! *******************************************************************
! *********** MANIPULATING ELEMENTARY SEISMOGRAMS (Green's) *********
! *******************************************************************
 
!
!     FILTERING THE ELEM. SEISMO - OLD !!!
!
 
!      do ir=1,nr
!         f1=ff1(ir)
!         f2=ff2(ir)
!         f3=ff3(ir)
!         f4=ff4(ir)
!        do icom=1,3
!          do it=1,6 ! here MUST be 6, not nmom
!          call filter(f1,f2,f3,f4,dt, w(1,ir,icom,it) )
!          enddo
!        enddo
!      enddo
!
 
 
 
 
!
!     FILTERING THE ELEM. SEISMO - NEW
!
!
      DO IR = 1, NR
         IF (KEYGPS(IR)==0) THEN
 
!         f1=0.1 ! pre-filter from raw data preparation (assuming HP)
!         f4=ff4(ir) ! not used in HP
!      keyfil=2 !(BP=0, LP=1, HP=2)
!      keydis=0 !(vel=0, displ=1)   !!!!!!! POZOR imitate pre-filter, no integration
!         shift=0.  ! POZOR imitate pre-filter 0 sec shift
!          do icom=1,3
!          do it=1,6 ! here MUST be 6, not nmom
!           call filter15(keyfil,keydis,shift,f1,f4,dt,w(1,ir,icom,it,isour))        ! seis. has the OLD filter
!             enddo
!          enddo
 
 
            F1 = FF1(IR)
            F4 = FF4(IR)
            KEYFIL = 0
               !(BP=0, LP=1, HP=2)   ! puvodni
            KEYDIS = 1
               !(vel=0, displ=1)
            SHIFT = 0.
                 ! shift 50 s is impossible if TL<50s !!!!!!!!!
            DO ICOM = 1, 3
               DO IT = 1, 6
                ! here MUST be 6, not nmom
                  CALL FILTER15(KEYFIL, KEYDIS, SHIFT, F1, F4, DT, W(1, IR, ICOM, IT, ISOUR))
               ENDDO
            ENDDO
 
         ELSE
 
            F1 = FF1(IR)
                     !!!!! for keygps = 1
            F4 = FF4(IR)
            KEYFIL = 1
               !(BP=0, LP=1, HP=2)
            KEYDIS = 1
               !(vel=0, displ=1)
            SHIFT = 0.
                     ! shift 50 s is impossible if TL<50s !!!!!!!!!
            DO ICOM = 1, 3
               DO IT = 1, 6
                ! here MUST be 6, not nmom
                  CALL FILTER15(KEYFIL, KEYDIS, SHIFT, F1, F4, DT, W(1, IR, ICOM, IT, ISOUR))
               ENDDO
            ENDDO
         ENDIF
 
         IF (KEYGPS(IR)==2) THEN
                               ! for 2 pure statics (after integration in filter CAUTION what filtration is used
            DO IT = 1, 6
               DO ITIM = MIN_TIME_SHIFT, MAX_TIME_SHIFT
!        do itim=1,ntim-1
 
!      w(itim,ir,1,it)=w(NUM_OF_TIME_SAMPLES,ir,1,it)  ! in case of problems can be substitued by some average at the 'end part'
!      w(itim,ir,2,it)=w(NUM_OF_TIME_SAMPLES,ir,2,it)        ! put 'c' if you want to see COMPLETE synth seis in FORWARD  modeling of GPS
!      w(itim,ir,3,it)=w(NUM_OF_TIME_SAMPLES,ir,3,it)        ! use of complete synth in case of keyGPS=2 not tested , probably possible but with lower VR
 
               ENDDO
            ENDDO
         ENDIF
      ENDDO ! over stations
 
 
 
!      do ir=1,nr  !!!only 1st station          ! output of elemse for plotting WITH FILTER
!        do it=1,6 !!! instead of 1,6
!          do itim=1,ntim
!          time=float(itim-1)*dt
!          write(8256,'(4(1x,e12.6))') time,
!     *          w(itim,ir,1,it),w(itim,ir,2,it),w(itim,ir,3,it)
!          enddo
!          write(8256,*) '*'
!        enddo
!      enddo
 
!       open(8256,file='testunc.dat')
!        do ir=1,nr   !!only 1st station         ! output of elemse for plotting WITH FILTER
!          do it=1,6
!            do itim=1,ntim
!           time=float(itim-1)*dt
!           write(8256,'(4(1x,e12.6))') time,
!     *          w(itim,ir,1,it),w(itim,ir,2,it),w(itim,ir,3,it)
!           enddo
!           write(8256,*) '*'
!         enddo
!       enddo
!       close(8256)
 
 
 
 
! ***********************************************************************
! ***********FORMING AND INVERTING SYSTEM MATRIX (from elem seismograms)
! ***********************************************************************

      MULTIP = 1.0E20
!
!     FORMING SYSTEM MATRIX
!
      DO N = 1, NMOM                                 !n,m mom tensor comp
         DO M = 1, NMOM
            R(N, M) = 0.
            DO ICOM = 1, 3
               DO IR = 1, NR
                  IF (STAT(IR)) THEN
                     DO ITIM = 1, NTIM
 
                        WNWEIG = W(ITIM, IR, ICOM, N, ISOUR)*WEIG(IR, ICOM)
                        WMWEIG = W(ITIM, IR, ICOM, M, ISOUR)*WEIG(IR, ICOM)
!                   write(*,*) wNweig, wMweig,wNweig*wMweig
                        R(N, M) = R(N, M) + WNWEIG*WMWEIG*MULTIP
!                                                     1.e20 = FORMAL SCALING
!        this formal scaling 1.e20 is essential to prevent some underflow
!        1.e20 is later compensated by the same multiplication of the
!        right-side of the equation
!        Note that matrix and right-hand side are formally multiplied
!             without affecting w, x, ...!!!
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
            R(N, M) = R(N, M)*DT
         ENDDO
      ENDDO
 
 
!        do i=1,nmom          ! printing sytem matrix
!        write(897,'(6(1x,e12.6))') (r(i,j),j=1,nmom)
!        enddo
!        write(897,*)
 
      DO I = 1, NMOM         ! saving r(n,m) before inversion
         DO J = 1, NMOM      ! (because after inversion r(n,m) is destroyed !)
            ROLD(I, J, ISOUR) = R(I, J)
            GOLD(I, J) = R(I, J)
         ENDDO
      ENDDO

      CALL JACOBINR(GOLD, NMOM, 6, EN, VV, NROT)
!        write(*,*) 'eigenvalues'
!              write(*,*) (en(i),i=1,nmom)
!        write(*,*) 'eigenvectors'
!              do i=1,nmom
!              write(*,*) (vv(i,j),j=1,nmom)
!              enddo
 
!
!     INVERTING THE SYSTEM MATRIX
!
 
! r(i,j)  matrix inversion (from Numerical Recipes) -> rinv(i,j)

      DO I = 1, NMOM
         DO J = 1, NMOM
            Y(I, J) = 0.
         ENDDO
         Y(I, I) = 1.
      ENDDO
      NMAX = 6    !!! dim. of array r
      CALL LUDCMP(R, NMOM, NMAX, INDX, DDD)   ! indx, ddd = output
      DO J = 1, NMOM                     ! former r is DESTROYED !!!
         CALL LUBKSB(R, NMOM, NMAX, INDX, Y(1, J))
      ENDDO
      DO I = 1, NMOM
         DO J = 1, NMOM
            RINV(I, J, ISOUR) = Y(I, J)
         ENDDO
      ENDDO
 
!       write(897,*)
!       do i=1,nmom          ! printing inverse matrix
!      write(897,'(6(1x,e12.6))') (rinv(i,j,isour),j=1,nmom)
!       write(897,'(6(1x,e12.6))') sqrt(rinv(i,i,isour))
!       enddo
!       write(897,*)
!
!       close(897)

      END SUBROUTINE ELEMAT15
