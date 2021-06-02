      SUBROUTINE ONEINV15(X, W, ROLD, RINV, ASAVE, CORR, ISH, ISOUR)

      USE COMMON_VARS
      USE PARAMETERS

      IMPLICIT NONE
!$OMP DECLARE TARGET
!$ACC ROUTINE SEQ
      REAL A, ASAVE, CORR, DELTA, G, PSI, RINV, ROLD, RR, SUM1, W, X, XMOM, XWEIG
      INTEGER I, ICOM, IR, ISEQ, ISEQM, ISH, ISHIFT, ISOUR, IT, ITIM, M, N
 
! Input: single elemse and 'data' (ori or residual)
! Temporal shift of elemse (in loop); finding moment tensor for each step
!finding correlation for each step
! Output: iseqm number of loop steps
! Output: array: ish(1,... iseqm) = shifts
! Output: array: corr(1,... iseqm) = correlation; sqrt(psi)=sqrt(1-delta/rr)
! delta=residual L2 misfit, rr= data squared; both with weights
! Output: array: asave(1,... iseqm) = moment tensor coefficients
 
      DIMENSION A(6), G(6)
      DIMENSION RINV(6, 6, MAX_SOURCE_POSITIONS)
      DIMENSION ROLD(6, 6, MAX_SOURCE_POSITIONS)
      DIMENSION X(NUM_OF_TIME_SAMPLES, MAX_STATIONS, 3)
      DIMENSION W(MIN_TIME_SHIFT:MAX_TIME_SHIFT, MAX_STATIONS, 3, 6, MAX_SOURCE_POSITIONS)
      DIMENSION CORR(MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION ISH(MAX_NUM_OF_TIME_SHIFTS)
      DIMENSION ASAVE(6, MAX_NUM_OF_TIME_SHIFTS)

      DO I = 1, 6
         A(I) = 0.
      ENDDO

!
!     COMPUTING 'POWER' OF THE DATA (observed or residual)
!                  (used stations ONLY !)
 
      RR = 0.
      DO ICOM = 1, 3                                    !icom seismo comp.
         DO IR = 1, NR                                  ! ir =station
            IF (STAT(IR)) THEN
               DO ITIM = 1, NTIM                           ! itim=time
                  XWEIG = X(ITIM, IR, ICOM)*WEIG(IR, ICOM)
                                                         ! icom seismo comp
                  RR = RR + XWEIG*XWEIG
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      RR = RR*DT
 
 
 
! *********************************************************************
! ********NON-LINEAR PART; THE SHIFT LOOP = GRID SEARCH OVER SHIFT ****
! *********************************************************************
 
!
!     STARTING  LOOP on SHIFT
!       (time shift = ishift*dt; shift > 0 ... to the RIGHT)
 
      ISEQM = (ILAST - IFIRST)/ISTEP
      ISHIFT = IFIRST
 
      DO ISEQ = 1, ISEQM
         ISHIFT = ISHIFT + ISTEP
                           ! value of the shift
 
!
!     FORMING RIGHT-HAND SIDE OF THE SYSTEM (dependent on the SHIFT)
!
 
 
         DO IT = 1, NMOM                                  !n moment tensor comp
            G(IT) = 0.                                    !it time
            DO ICOM = 1, 3                                !icom seismo comp.
               DO IR = 1, NR                              !ir station
                  IF (STAT(IR)) THEN
                     DO ITIM = 1, NTIM
                        XWEIG = X(ITIM, IR, ICOM)*WEIG(IR, ICOM)**2.
                        G(IT) = G(IT) + W(ITIM - ISHIFT, IR, ICOM, IT, ISOUR)*XWEIG*1.E20
!                                                  1.e20 = FORMAL SCALING
!                     (compensating the previous scaling of the system matrix)
!               ! this MINUS ishift in w is equivalent to PLUS in x
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
            G(IT) = G(IT)*DT
         ENDDO
 
 
!
!     APPLYING INVERSE MATRIX TO THE RIGHT-HAND SIDE (WITHIN SHIFT LOOP)
!
 
 
         DO N = 1, NMOM                   ! you MUST disable this part  !
            A(N) = 0.                     !                             !
            DO M = 1, NMOM                !      C A U T I O N          !
               A(N) = A(N) + RINV(N, M, ISOUR)*G(M)
                                          !                             !
            ENDDO                         ! when foc mech is prescribed !
         ENDDO
 
         SUM1 = 0.
         DO N = 1, NMOM
            SUM1 = SUM1 + A(N)*G(N)/1.E20
                                     ! g de-scaled by 1.e20
         ENDDO
 
!      sum2=0.
!      do n=1,nmom
!      do m=1,nmom
!      sum2=sum2 + rold(n,m,isour)*a(n)*a(m)/1.e20  ! rold de-scaled
!      enddo
!      enddo
!      write(567,*) 'sum1,sum2=',sum1,sum2
 
         XMOM = 1.
 
!     delta=rr - 2.*sum1 + sum2   !
         DELTA = RR - SUM1        ! the same (since sum1=sum2)
 
         PSI = SUM1/RR      ! psi = correlation
 
 
         IF (PSI<0) THEN                         !!!!! new July 5, 2014
            PSI = 0.
!cc         xmom=-1.*xmom
         ELSE
            PSI = SQRT(PSI)
                         ! new
         ENDIF
 
 
! *******************************************************
!      no DC constraint here (it is only in cnsinv)
! *******************************************************
 
 
         CORR(ISEQ) = PSI
                        ! saving correlation; iseq counts the shifts
         ISH(ISEQ) = ISHIFT
 
         DO N = 1, 6  ! here it must be 6, not nmom
            ASAVE(N, ISEQ) = A(N)*XMOM
                                   ! saving  MOMENT COEFFICIENTS
         ENDDO
 
 
      ENDDO         ! SHIFT LOOP FINISHED

      END SUBROUTINE ONEINV15
