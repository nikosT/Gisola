INTERFACE

      SUBROUTINE LP(P, Z, RTYPE, DCVALUE, NSECTS, SN, SD)
      IMPLICIT NONE
!$OMP DECLARE TARGET
!$ACC ROUTINE SEQ
      INTEGER NSECTS
      COMPLEX P(*), Z(*)
      CHARACTER*4 RTYPE(*)
      REAL*4 SN(*), SD(*), DCVALUE

      END SUBROUTINE LP

END INTERFACE

