      SUBROUTINE ELEMSE(W, LOCAL_NTIM, TL, AW, NFREQ, ISOUR, FILENUM_HES)

      USE COMMON_VARS
      USE PARAMETERS

      IMPLICIT NONE

      INTEGER LOCAL_NTIM, NFREQ, ISOUR, FILENUM_HES
      REAL W, TL, AW

      DIMENSION W(MIN_TIME_SHIFT:MAX_TIME_SHIFT, MAX_STATIONS, 3, 6, MAX_SOURCE_POSITIONS)

      INTEGER    IWK(NUM_OF_TIME_SAMPLES), JF, ICC, ICS, IFR, IR, IT, ITIM
      REAL       DFREQ, FREQ, CK, CC, PAS, T0, T1, TIME, LOCAL_DT, LOCAL_AW
      COMPLEX*16 UX(NUM_OF_TIME_SAMPLES, MAX_STATIONS, 6)
      COMPLEX*16 UY(NUM_OF_TIME_SAMPLES, MAX_STATIONS, 6)
      COMPLEX*16 UZ(NUM_OF_TIME_SAMPLES, MAX_STATIONS, 6)
      COMPLEX*16 UXF(6), UYF(6), UZF(6), DERIV, FSOURCE, US, OMEGA

      UX = (0., 0.)
      UY = (0., 0.)
      UZ = (0., 0.)

      ICS = 7  ! true step
      T0  = 0.
      T1  = 0.
      ICC = 1  ! output = velocity

      DFREQ = 1./TL
      PAS = TL/FLOAT(LOCAL_NTIM)
      LOCAL_AW = -PI*AW/TL           ! aw re-defined here

      DO JF = 1, NFREQ         ! loop over frequency
         FREQ = FLOAT(JF - 1)/TL
         OMEGA = DCMPLX(PI2*FREQ, LOCAL_AW)
         DERIV = (AI*OMEGA)**ICC
         US = FSOURCE(ICS, OMEGA, T0, T1, PAS)*DERIV

         DO IR = 1, NR         ! loop over stations
            READ (FILENUM_HES) (UXF(IT), IT = 1, 6)
            READ (FILENUM_HES) (UYF(IT), IT = 1, 6)
            READ (FILENUM_HES) (UZF(IT), IT = 1, 6)

            DO IT = 1, 6       ! moment tensor components
               UX(JF, IR, IT) = UXF(IT)*US
               UY(JF, IR, IT) = UYF(IT)*US
               UZ(JF, IR, IT) = UZF(IT)*US
            ENDDO
         ENDDO                 ! end loop over stations
      ENDDO                    ! ending loop over frequency

!++++++++++++
!                ELEMENTARY  SEISMO  CALCUL
!++++++++++++

      DO IR = 1, NR                          ! stations
         DO IT = 1, 6                        ! moment tensor
            DO IFR = LOCAL_NTIM + 2 - NFREQ, LOCAL_NTIM  ! frequency
               UX(IFR, IR, IT) = CONJG(UX(LOCAL_NTIM + 2 - IFR, IR, IT))
               UY(IFR, IR, IT) = CONJG(UY(LOCAL_NTIM + 2 - IFR, IR, IT))
               UZ(IFR, IR, IT) = CONJG(UZ(LOCAL_NTIM + 2 - IFR, IR, IT))
            ENDDO

            CALL FFT2CD(UX(1, IR, IT), NUM_OF_TIME_SAMPLES_LOG2, IWK)
            CALL FFT2CD(UY(1, IR, IT), NUM_OF_TIME_SAMPLES_LOG2, IWK)
            CALL FFT2CD(UZ(1, IR, IT), NUM_OF_TIME_SAMPLES_LOG2, IWK)

            DO ITIM = 1, LOCAL_NTIM                ! time
               CK = FLOAT(ITIM - 1)/LOCAL_NTIM
               CC = EXP( - LOCAL_AW*TL*CK)/TL      ! removing artific. absorption

!!!!!!!!!!       (factor 1./tl has nothing to do with it, must be always)
!      cc=1./tl ! this will suppress exp() but keep division by tl
!!!!!!!!!! in other words - if we experimentate with calcelling this operation
!!!!!!!!!!!!!!    we must at least divite u(...) by TL

               UX(ITIM, IR, IT) = UX(ITIM, IR, IT)*CC  ! this is regular (then in filter_stat we need artif. atten)
               UY(ITIM, IR, IT) = UY(ITIM, IR, IT)*CC
               UZ(ITIM, IR, IT) = UZ(ITIM, IR, IT)*CC
            ENDDO                            ! frequency
         ENDDO                               ! moment tensor
      ENDDO                                  ! stations

!ccc    ELEMENTARY SEISMO  (it=1,2...6)  OUTPUT   !!!  always all 6 !!!

      LOCAL_DT = TL/FLOAT(LOCAL_NTIM)

      DO IR = 1, NR            ! stations
         DO IT = 1, 6          ! seismo for 6 basic moment tensors
            DO ITIM = 1, LOCAL_NTIM  ! time
               TIME = FLOAT(ITIM - 1)*LOCAL_DT
               W(ITIM, IR, 1, IT, ISOUR) = REAL(UX(ITIM, IR, IT))
               W(ITIM, IR, 2, IT, ISOUR) = REAL(UY(ITIM, IR, IT))
               W(ITIM, IR, 3, IT, ISOUR) = REAL(UZ(ITIM, IR, IT))
!!!!!               IF (ITIM < 6) THEN
!!!!!                  WRITE(*, *) TIME, W(ITIM, IR, 1, IT, ISOUR), W(ITIM, IR, 2, IT, ISOUR), W(ITIM, IR, 3, IT, ISOUR)
!!!!!               ENDIF
            ENDDO              ! time
         ENDDO                 ! basic moment tensors
      ENDDO                    ! stations

      END SUBROUTINE ELEMSE

