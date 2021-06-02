!      real strk,dip,rake
 
 
!      write(*,*) 'Strike,Dip, rake (deg)=?'
!      read(*,*) strk, dip, rake
 
 
!      call pl2pt(strike,dip,rake,trendp,plungp,trendt,plungt,
!     1trendb,plungb,ierr)
 
!      write(*,*)trendp,plungp,trendt,plungt
 
!
!      stop
!      end
 
 
      SUBROUTINE PL2PT(STRIKE, DIP, RAKE, TRENDP, PLUNGP, TRENDT, PLUNGT, TRENDB, PLUNGB, IERR)
      IMPLICIT NONE
!
!     compute trend and plunge of P, T and B axes
!     from strike, dip and rake of a nodal plane
!
!
!     usage:
!     call pl2pt(strike,dip,rake,trendp,plungp,trendt,plungt,trendb,plungb,ierr)
!
!     arguments:
!     strike         strike angle in degrees of the first nodal plane (INPUT)
!     dip            dip angle in degrees of the first nodal plane (INPUT)
!     rake           rake angle in degrees of the first nodal plane (INPUT)
!     trendp         trend of P axis (OUTPUT)
!     plungp         plunge or P axis (OUTPUT)
!     trendt         trend of T axis (OUTPUT)
!     plungt         plunge or T axis (OUTPUT)
!     trendb         trend of B axis (OUTPUT)
!     plungb         plunge or B axis (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input STRIKE angle out of range
!     2              input DIP angle out of range
!     4              input RAKE angle out of range
!     3              1+2
!     5              1+4
!     7              1+2+4
!     8,9,10,11      internal error
!
!      implicit none
!-------------------------------------------------------------------------------
      INTEGER IO
      REAL AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR, C360, C90,&
           C0, C1, C2, C3
      COMMON /FPSCOM/ AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR,&
                      C360, C90, C0, C1, C2, C3, IO
!-------------------------------------------------------------------------------
      REAL STRIKE, DIP, RAKE, ANX, ANY, ANZ, DX, DY, DZ, PX, PY, PZ, TX, TY, TZ, BX, BY, BZ, TRENDP, PLUNGP, TRENDT, PLUNGT,       &
           TRENDB, PLUNGB
      INTEGER IERR
!
      CALL FPSSET
!
      CALL PL2ND(STRIKE, DIP, RAKE, ANX, ANY, ANZ, DX, DY, DZ, IERR)
      IF (IERR/=0) THEN
         WRITE (IO, '(1x,a,i3)') 'PL2PT: ierr=', IERR
         RETURN
      ENDIF
      CALL ND2PT(DX, DY, DZ, ANX, ANY, ANZ, PX, PY, PZ, TX, TY, TZ, BX, BY, BZ, IERR)
      IF (IERR/=0) THEN
         IERR = 8
         WRITE (IO, '(1x,a,i3)') 'PL2PT: ierr=', IERR
      ENDIF
      CALL CA2AX(PX, PY, PZ, TRENDP, PLUNGP, IERR)
      IF (IERR/=0) THEN
         IERR = 9
         WRITE (IO, '(1x,a,i3)') 'PL2PT: ierr=', IERR
      ENDIF
      CALL CA2AX(TX, TY, TZ, TRENDT, PLUNGT, IERR)
      IF (IERR/=0) THEN
         IERR = 10
         WRITE (IO, '(1x,a,i3)') 'PL2PT: ierr=', IERR
      ENDIF
      CALL CA2AX(BX, BY, BZ, TRENDB, PLUNGB, IERR)
      IF (IERR/=0) THEN
         IERR = 11
         WRITE (IO, '(1x,a,i3)') 'PL2PT: ierr=', IERR
      ENDIF
      END SUBROUTINE PL2PT
 
 
 
      SUBROUTINE FPSSET
      IMPLICIT NONE
!
!     define constants (i.e. input ranges and tolerances) used throughout the
!     package. It is called by every subroutines to setup constants
!
!     usage:
!     call fpsset
!
!     constants in fpscom common block:
!
!     amistr         strike lower limit
!     amastr         strike upper limit
!     amidip         dip lower limit
!     amadip         dip upper limit
!     amirak         rake lower limit
!     amarak         rake upper limit
!     amitre         trend lower limit
!     amatre         trend upper limit
!     amiplu         plunge lower limit
!     amaplu         plunge upper limit
!     orttol         orthogonality tolerance
!     ovrtol         dip overtaking tolerance
!     tentol         moment tensor symmetry tolerance
!     dtor            degree to radians
!     c360            360.
!     c90            90.
!     c0             0.
!     c1             1.
!     c2             2.
!     c3             3.
!     io             error messages file unit
!
!      implicit none
!-------------------------------------------------------------------------------
      INTEGER IO
      REAL AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR, C360, C90,&
           C0, C1, C2, C3
      COMMON /FPSCOM/ AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR,&
                      C360, C90, C0, C1, C2, C3, IO
!-------------------------------------------------------------------------------
      INTEGER IFL
      SAVE IFL
      DATA IFL/0/
      IF (IFL==0) THEN
         AMISTR = -360.
         AMASTR = 360.
         AMIDIP = 0.
         AMADIP = 90.
         AMIRAK = -360.
         AMARAK = 360.
         AMITRE = -360.
         AMATRE = 360.
         AMIPLU = 0.
         AMAPLU = 90.
         ORTTOL = 2.
         OVRTOL = 0.001
         TENTOL = 0.0001
         DTOR = 0.017453292519943296
         C360 = 360.
         C90 = 90.
         C0 = 0.
         C1 = 1.
         C2 = 2.
         C3 = 3.
         IO = 6
         IFL = 1
      ENDIF
      END SUBROUTINE FPSSET
 
 
 
      SUBROUTINE PL2ND(STRIKE, DIP, RAKE, ANX, ANY, ANZ, DX, DY, DZ, IERR)
      IMPLICIT NONE
!
!     compute Cartesian components of outward normal and slip
!     vectors from strike, dip and rake
!
!     usage:
!     call pl2nd(strike,dip,rake,anx,any,anz,dx,dy,dz,ierr)
!
!     arguments:
!     strike         strike angle in degrees (INPUT)
!     dip            dip angle in degrees (INPUT)
!     rake           rake angle in degrees (INPUT)
!     anx,any,anz    components of fault plane outward normal versor in the
!                    Aki-Richards Cartesian coordinate system (OUTPUT)
!     dx,dy,dz       components of slip versor in the Aki-Richards
!                    Cartesian coordinate system (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input STRIKE angle out of range
!     2              input DIP angle out of range
!     4              input RAKE angle out of range
!     3              1+2
!     5              1+4
!     7              1+2+4
!
!      implicit none
!-------------------------------------------------------------------------------
      INTEGER IO
      REAL AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR, C360, C90,&
           C0, C1, C2, C3
      COMMON /FPSCOM/ AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR,&
                      C360, C90, C0, C1, C2, C3, IO
!-------------------------------------------------------------------------------
      REAL ANX, ANY, ANZ, DX, DY, DZ, STRIKE, DIP, RAKE, WSTRIK, WDIP, WRAKE
      INTEGER IERR
!
      CALL FPSSET
!
      ANX = C0
      ANY = C0
      ANZ = C0
      DX = C0
      DY = C0
      DZ = C0
      IERR = 0
      IF (STRIKE<AMISTR .OR. STRIKE>AMASTR) THEN
         WRITE (IO, '(1x,a,g10.4,a)') 'PL2ND: input STRIKE angle ', STRIKE, ' out of range'
         IERR = 1
      ENDIF
      IF (DIP<AMIDIP .OR. DIP>AMADIP) THEN
         IF (DIP<AMADIP .AND. DIP> - OVRTOL) THEN
            DIP = AMIDIP
         ELSEIF (DIP>AMIDIP .AND. DIP - AMADIP<OVRTOL) THEN
            DIP = AMADIP
         ELSE
            WRITE (IO, '(1x,a,g10.4,a)') 'PL2ND: input DIP angle ', DIP, ' out of range'
            IERR = IERR + 2
         ENDIF
      ENDIF
      IF (RAKE<AMIRAK .OR. RAKE>AMARAK) THEN
         WRITE (IO, '(1x,a,g10.4,a)') 'PL2ND: input RAKE angle ', RAKE, ' out of range'
         IERR = IERR + 4
      ENDIF
      IF (IERR/=0) RETURN
      WSTRIK = STRIKE*DTOR
      WDIP = DIP*DTOR
      WRAKE = RAKE*DTOR
!
      ANX = -SIN(WDIP)*SIN(WSTRIK)
      ANY = SIN(WDIP)*COS(WSTRIK)
      ANZ = -COS(WDIP)
      DX = COS(WRAKE)*COS(WSTRIK) + COS(WDIP)*SIN(WRAKE)*SIN(WSTRIK)
      DY = COS(WRAKE)*SIN(WSTRIK) - COS(WDIP)*SIN(WRAKE)*COS(WSTRIK)
      DZ = -SIN(WDIP)*SIN(WRAKE)
      END SUBROUTINE PL2ND
 
 
      SUBROUTINE ND2PT(WANX, WANY, WANZ, WDX, WDY, WDZ, PX, PY, PZ, TX, TY, TZ, BX, BY, BZ, IERR)
      IMPLICIT NONE
!
!     compute Cartesian component of P, T and B axes from outward normal
!     and slip vectors
!
!     usage:
!     call nd2pt(anx,any,anz,dx,dy,dz,px,py,pz,tx,ty,tz,bx,by,bz,ierr)
!
!     arguments:
!     anx,any,anz    components of fault plane outward normal vector in the
!                    Aki-Richards Cartesian coordinate system (INPUT)
!     dx,dy,dz       components of slip vector in the Aki-Richards
!                    Cartesian coordinate system (INPUT)
!     px,py,pz       components of downward P (maximum dilatation) axis versor
!                    in the Aki-Richards Cartesian coordinate system (OUTPUT)
!     tx,ty,tz       components of downward T (maximum tension) axis versor
!                    in the Aki-Richards Cartesian coordinate system (OUTPUT)
!     bx,by,bz       components of downward B (neutral) axis versor in the
!                    Aki-Richards Cartesian coordinate system (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input vectors not perpendicular among each other
!
!      implicit none
!-------------------------------------------------------------------------------
      INTEGER IO
      REAL AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR, C360, C90,&
           C0, C1, C2, C3
      COMMON /FPSCOM/ AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR,&
                      C360, C90, C0, C1, C2, C3, IO
!-------------------------------------------------------------------------------
      REAL WANX, WANY, WANZ, AMN, ANX, ANY, ANZ, WDX, WDY, WDZ, AMD, DX, DY, DZ, ANG, PX, PY, PZ, TX, TY, TZ, BX, BY, BZ, AMP
      INTEGER IERR
!
      CALL FPSSET
!
      IERR = 0
      CALL NORM(WANX, WANY, WANZ, AMN, ANX, ANY, ANZ)
      CALL NORM(WDX, WDY, WDZ, AMD, DX, DY, DZ)
      CALL ANGLE(ANX, ANY, ANZ, DX, DY, DZ, ANG)
      IF (ABS(ANG - C90)>ORTTOL) THEN
         WRITE (IO, '(1x,a,g15.7,a)') 'ND2PT: input vectors not '//'perpendicular, angle=', ANG
         IERR = 1
      ENDIF
      PX = ANX - DX
      PY = ANY - DY
      PZ = ANZ - DZ
      CALL NORM(PX, PY, PZ, AMP, PX, PY, PZ)
      IF (PZ<C0) CALL INVERT(PX, PY, PZ)
      TX = ANX + DX
      TY = ANY + DY
      TZ = ANZ + DZ
      CALL NORM(TX, TY, TZ, AMP, TX, TY, TZ)
      IF (TZ<C0) CALL INVERT(TX, TY, TZ)
      CALL VECPRO(PX, PY, PZ, TX, TY, TZ, BX, BY, BZ)
      IF (BZ<C0) CALL INVERT(BX, BY, BZ)
      END SUBROUTINE ND2PT
 
      SUBROUTINE CA2AX(WAX, WAY, WAZ, TREND, PLUNGE, IERR)
      IMPLICIT NONE
!
!     compute trend and plunge from Cartesian components
!
!     usage:
!     call ca2ax(ax,ay,az,trend,plunge,ierr)
!
!     arguments:
!     ax,ay,az       components of axis direction vector in the Aki-Richards
!                    Cartesian coordinate system (INPUT)
!     trend          clockwise angle from North in degrees (OUTPUT)
!     plunge         inclination angle in degrees (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!
!      implicit none
!-------------------------------------------------------------------------------
      INTEGER IO
      REAL AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR, C360, C90,&
           C0, C1, C2, C3
      COMMON /FPSCOM/ AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR,&
                      C360, C90, C0, C1, C2, C3, IO
!-------------------------------------------------------------------------------
      REAL WAX, WAY, WAZ, WNORM, AX, AY, AZ, TREND, PLUNGE
      INTEGER IERR
!
      CALL FPSSET
!
      IERR = 0
      CALL NORM(WAX, WAY, WAZ, WNORM, AX, AY, AZ)
      IF (AZ<C0) CALL INVERT(AX, AY, AZ)
      IF (AY/=C0 .OR. AX/=C0) THEN
         TREND = ATAN2(AY, AX)/DTOR
      ELSE
         TREND = C0
      ENDIF
      TREND = MOD(TREND + C360, C360)
      PLUNGE = ASIN(AZ)/DTOR
      END SUBROUTINE CA2AX
 
 
      SUBROUTINE NORM(WAX, WAY, WAZ, ANORM, AX, AY, AZ)
      IMPLICIT NONE
!
!     compute euclidean norm and versor components
!
!     usage:
!     call norm(wax,way,waz,anorm,ax,ay,az)
!
!     arguments:
!     wax,way,waz    Cartesian component of input vector (INPUT)
!     anorm          Euclidean norm of input vector (OUTPUT)
!     ax,ay,az       normalized Cartesian component of the vector (OUTPUT)
!
!      implicit none
!-------------------------------------------------------------------------------
      INTEGER IO
      REAL AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR, C360, C90,&
           C0, C1, C2, C3
      COMMON /FPSCOM/ AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR,&
                      C360, C90, C0, C1, C2, C3, IO
!-------------------------------------------------------------------------------
      REAL WAX, WAY, WAZ, ANORM, AX, AY, AZ
!
      CALL FPSSET
!
      ANORM = SQRT(WAX*WAX + WAY*WAY + WAZ*WAZ)
      IF (ANORM==C0) RETURN
      AX = WAX/ANORM
      AY = WAY/ANORM
      AZ = WAZ/ANORM
      END SUBROUTINE NORM
 
 
      SUBROUTINE VECPRO(PX, PY, PZ, TX, TY, TZ, BX, BY, BZ)
      IMPLICIT NONE
!
!     compute vector products of two vectors
!
!     usage:
!     call vecpro(px,py,pz,tx,ty,tz,bx,by,bz)
!
!     arguments:
!
!     px,py,pz       Cartesian component of first vector (INPUT)
!     tx,ty,tz       Cartesian component of second vector (INPUT)
!     bx,by,bz       Cartesian component of vector product (OUTUT)
!
!      implicit none
!-------------------------------------------------------------------------------
      INTEGER IO
      REAL AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR, C360, C90,&
           C0, C1, C2, C3
      COMMON /FPSCOM/ AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR,&
                      C360, C90, C0, C1, C2, C3, IO
!-------------------------------------------------------------------------------
      REAL PX, PY, PZ, TX, TY, TZ, BX, BY, BZ
!
      CALL FPSSET
!
      BX = PY*TZ - PZ*TY
      BY = PZ*TX - PX*TZ
      BZ = PX*TY - PY*TX
      END SUBROUTINE VECPRO
 
 
      SUBROUTINE INVERT(AX, AY, AZ)
      IMPLICIT NONE
!
!     invert vector
!
!     usage:
!     utility routine for internal use only
!
!      implicit none
!-------------------------------------------------------------------------------
      INTEGER IO
      REAL AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR, C360, C90,&
           C0, C1, C2, C3
      COMMON /FPSCOM/ AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR,&
                      C360, C90, C0, C1, C2, C3, IO
!-------------------------------------------------------------------------------
      REAL AX, AY, AZ
!
      CALL FPSSET
!
      AX = -AX
      AY = -AY
      AZ = -AZ
      END SUBROUTINE INVERT
 
      SUBROUTINE ANGLE(WAX, WAY, WAZ, WBX, WBY, WBZ, ANG)
      IMPLICIT NONE
!
!     compute the angle (in degrees) between two vectors
!
!     usage:
!     call angle(wax,way,waz,wbx,wby,wbz,ang)
!
!     arguments:
!     wax,way,waz    Cartesian component of first vector (INPUT)
!     wbx,wby,wbz    Cartesian component of second vector (INPUT)
!     ang            angle between the two vectors in degrees (OUTPUT)
!
!      implicit none
!-------------------------------------------------------------------------------
      INTEGER IO
      REAL AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR, C360, C90,&
           C0, C1, C2, C3
      COMMON /FPSCOM/ AMISTR, AMASTR, AMIDIP, AMADIP, AMIRAK, AMARAK, AMITRE, AMATRE, AMIPLU, AMAPLU, ORTTOL, OVRTOL, TENTOL, DTOR,&
                      C360, C90, C0, C1, C2, C3, IO
!-------------------------------------------------------------------------------
      REAL WAX, WAY, WAZ, WBX, WBY, WBZ, AX, AY, AZ, BX, BY, BZ, ANG, ANORM, BNORM, PROD
!
      CALL FPSSET
!
      CALL NORM(WAX, WAY, WAZ, ANORM, AX, AY, AZ)
      CALL NORM(WBX, WBY, WBZ, BNORM, BX, BY, BZ)
      PROD = AX*BX + AY*BY + AZ*BZ
      ANG = ACOS(MAX( - C1, MIN(C1, PROD)))/DTOR
      END SUBROUTINE ANGLE
