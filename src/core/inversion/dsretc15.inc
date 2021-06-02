!+
      SUBROUTINE DSRETC(strike, dip, rake)
      IMPLICIT NONE
      REAL strike, dip, rake
      REAL ANBTP, ANGS, ANGS2, PTTP
      INTEGER J
!
!	Based on the Aki & Richards convention, gives all
!	  representations of a focal mechanism for input
!	  Dip, Strike, Rake    or
!	  A and N trend and plunge    or
!	  P and T trend and plunge.
!
!	The code is Fortran 77 for a VAX/VMS system
!	Subroutine specific to this program are
!	  FMREPS,PTTPIN,ANTPIN,DSRIN,AN2MOM,V2TRPL,TRPL2V,AN2DSR
!	21 August 1991:  sun version
!	30 December 1993 Incuded Stuart Sipkin's beachball
!	15 March 2002: Added moment tensor input
!-
      LOGICAL AN, PT, DSR, TRUTH, FIRST, MT
      CHARACTER*80 GETSTRING, COMMNT
      REAL*4 MOMTEN(6)
      DIMENSION PTTP(4), ANGS(3), ANGS2(3), ANBTP(6)
!
!	If your compiler complains about form='print', leave it out.
!
      OPEN (UNIT = 2, FILE = 'dsretc.lst', STATUS = 'unknown')
!      OPEN (UNIT = 3, FILE = 'dsr.dat')
!	write(2,'(/,5x,a)') commnt(1:lenc(commnt))
!	WRITE(*,*) 'Can enter P & T or D, S & R or A & N'
      DSR = .FALSE.
 
      AN = .FALSE.
      PT = .FALSE.
      MT = .FALSE.
!	IF (TRUTH('Dip, Strike and Rake?..[Y]')) THEN
      DSR = .TRUE.
!	  CALL PRINTX(' Enter Dip, Strike and Rake (degrees)')
 
 
!      READ (3, *) (ANGS(J), J = 1, 3)

     ANGS(1)=strike
     ANGS(2)=dip
     ANGS(3)=rake
!	ELSE IF (TRUTH('P and T axes trend and plunge?..[Y]')) THEN
!	  PT = .TRUE.
!	  CALL PRINTX
!     .	  (' Enter trend and plunge for P and T (t,p,t,p)')
!	  READ(*,*) (PTTP(J),J=1,4)
!	else if (truth('Moment-tensor input?...[Y]')) then
!	  MT = .true.
!	ELSE
!	  AN = .TRUE.
!	  CALL PRINTX
!     .	  (' Enter trend and plunge for A and N (t,p,t,p)')
!	  READ(*,*) (ANBTP(J),J=1,4)
!	END IF
 
!	WRITE(2,'(1H0)')
 
 
      CALL FMREPS(ANBTP, ANGS, PTTP, ANGS2, AN, PT, DSR, MT, MOMTEN, 2, 6)
      FIRST = .TRUE.
      CALL BBALL(MOMTEN, PTTP(1), PTTP(2), PTTP(3), PTTP(4), 2, FIRST)
!	IF (.NOT.TRUTH('Run some more?...[Y]')) STOP
!	write(2,'(1H1)')
!        first = .true.
!	IF (TRUTH('Same comment?..[Y]')) GO TO 200
!	GO TO 100
 
 
      END SUBROUTINE DSRETC
 
!+
      SUBROUTINE BBALL(G, PAZIM, PPLNG, TAZIM, TPLNG, UNIT, FIRST)
      IMPLICIT NONE
      REAL AOI, DD, DEL, DX, DY, G, PAZIM, PI, POLAR, PPLNG, R0, RADIUS, TAZIM, THETA, TPLNG, X0, XPAZIM, XTAZIM, XTHETA, Y0
      INTEGER I, IX0, IXP, IXT, IY0, IYP, IYT, J
 
! ...... generate printer plot rendition of lower hemisphere
!        equal area projection
!	g has the six elements of the moment tensor, the rest are the
!	  plunge and trends of the P and T axes in degrees. unit is the output
!	  unit.
!	From Stuart Sipkin and Bob Uhrhammer 1993
!-
      DIMENSION G(6)
      INTEGER UNIT
      CHARACTER*1 ACH(39, 72), APLUS, AMINUS, APAXIS, ATAXIS, ABLANK
      LOGICAL FIRST
!
      DATA APLUS, AMINUS, APAXIS, ATAXIS, ABLANK/'#', '-', 'P', 'T', ' '/
      DATA RADIUS/1.41/
!
! ...... construct lower hemisphere fps
!
 
      PI = 3.14159265
 
      R0 = RADIUS
      X0 = R0 + 0.250
      Y0 = R0 + 0.500
      IX0 = 12.*X0
      IY0 = 6.5*Y0
      DO I = 1, 2*IX0
         DO J = 1, 2*IY0
            DX = REAL(I - IX0)/12.
            DY = -REAL(J - IY0)/6.5
            DD = DX*DX + DY*DY
            IF (DD>0.) THEN
               DEL = SQRT(DD)
            ELSE
               DEL = 0.
            ENDIF
            IF ((DX==0.) .AND. (DY==0.)) THEN
               THETA = 0.
            ELSE
               XTHETA = ATAN2(DX, DY)
               THETA = 180.*(XTHETA/PI)
            ENDIF
            IF (DEL>R0) THEN
               ACH(J, I) = ABLANK
               CYCLE
            ENDIF
            IF (DEL>=R0) THEN
               AOI = 90.0
            ELSE
               AOI = 90.*DEL/R0
            ENDIF
            IF (POLAR(G, AOI, THETA, FIRST)>0.) THEN
               ACH(J, I) = APLUS
            ELSE
               ACH(J, I) = AMINUS
            ENDIF
         ENDDO
      ENDDO
!
! ...... add P & T axis
!
      XPAZIM = PI*PAZIM/180.
      IXP = NINT(R0*12.*(90. - PPLNG)*SIN(XPAZIM)/90. + REAL(IX0))
      IYP = NINT( - R0*6.5*(90. - PPLNG)*COS(XPAZIM)/90. + REAL(IY0))
 
      DO I = IXP - 1, IXP + 1
         DO J = IYP - 1, IYP + 1
            ACH(J, I) = ABLANK
         ENDDO
      ENDDO
 
      ACH(IYP, IXP) = APAXIS
 
      XTAZIM = PI*TAZIM/180.
      IXT = NINT(R0*12.*(90. - TPLNG)*SIN(XTAZIM)/90. + REAL(IX0))
      IYT = NINT( - R0*6.5*(90. - TPLNG)*COS(XTAZIM)/90. + REAL(IY0))
 
      DO I = IXT - 1, IXT + 1
         DO J = IYT - 1, IYT + 1
            ACH(J, I) = ABLANK
         ENDDO
      ENDDO
 
      ACH(IYT, IXT) = ATAXIS
!
! ...... add fps plot
!
      DO I = 1, 2*IY0 - 2
         WRITE (UNIT, '(5x,72a1)') (ACH(I, J), J = 1, 2*IX0)
      ENDDO
 
      DO I = 1, 2*IY0 - 2
         WRITE (*, '(5x,72a1)') (ACH(I, J), J = 1, 2*IX0)
      ENDDO
 
      END SUBROUTINE BBALL
 
 
      REAL*4 FUNCTION POLAR(G, AOI, THETA, FIRST)
      IMPLICIT NONE
      REAL AOI, G, PI, THETA, X, XAOI, XTHETA, Y, Z
!
! ...... compute first motion polarity as a function of aoi & theta
!        for a moment tensor for a double-couple solution.
!	Conventions differ slightly from Sipkin.  My moment tensor is derived
!	  from the outer product of two vectors and is hence normalized.  The
!	  order is also different from his, apparently.
 
      DIMENSION G(6)
      REAL MXX, MXY, MXZ, MYY, MYZ, MZZ
      LOGICAL FIRST
 
      PI = 3.14159265
 
!      if(first) then   ! t.2012
      MXX = G(2)
      MXY = -G(6)
      MXZ = G(4)
      MYY = G(3)
      MYZ = -G(5)
      MZZ = G(1)
!        first = .false.
!      endif
 
      XTHETA = PI*THETA/180.
      XAOI = PI*AOI/180.
      X = COS(XTHETA)*SIN(XAOI)
      Y = SIN(XTHETA)*SIN(XAOI)
      Z = COS(XAOI)
!
      POLAR = X*MXX*X + 2*X*MXY*Y + 2*X*MXZ*Z + 2*Y*MYZ*Z + Y*MYY*Y + Z*MZZ*Z
!
      END FUNCTION POLAR
 
 
!+
      SUBROUTINE FMREPS(ANBTP, ANGS, PTTP, ANGS2, AN, PT, DSR, MT, MOMTEN, LUNIT1, LUNIT2)
      IMPLICIT NONE
      REAL ANBTP, ANGS, ANGS2, PI, PTTP, RDEG
      INTEGER I, J
!
!	Input: A and N (trned and plunge), or P and T or dip, strike
!	  and rake (A&R convention)
!	Output: the other representations plus the auxiliary plane.
!	PTTP:  4 parameters, trend and plunge for P and T
!	  P trend and plunge, T trend and plunge
!	ANBTP  6 parameters, t and p for A, N and B respectively.
!	ANGS   3 parameters,  dip, strike and rake for first plane
!	ANGS2  3 parameters, dip strike and rake for auxiliary plane
!	AN, PT, DSR are LOGICAL variables which are true if
!	  A and N, P and T or dip-strike-rake are input
!	MOMTEN  6 parameters:  the moment tensor (unit scalar
!	  magnitude) D & W convention
!	Angles come in and go out in degrees.
!	If LUNIT1 and/or LUNIT2 are positive, the representations
!	  are written on those logical unit numbers.
!	22 July 1985:  Added moment tensor output
!	30 October 1992:  Fixed up problem if lunit2=5.
!	7 March 2002:  Added moment-tensor input
!-
      LOGICAL AN, PT, DSR, MT
      REAL*4 MOMTEN(6)
      INTEGER LUNIT1, LUNIT2
      DIMENSION PTTP(4), ANGS(3), ANGS2(3), ANBTP(6)
      RDEG = 45.0/ATAN(1.0)
      PI = 4.0*ATAN(1.0)
      IF (MT) THEN
         CALL MT_IN(PTTP, PI)
         CALL PTTPIN(PTTP, ANGS, ANGS2, ANBTP, MOMTEN, PI)
      ELSEIF (PT) THEN
         DO J = 1, 4
            PTTP(J) = PTTP(J)/RDEG
         ENDDO
         CALL PTTPIN(PTTP, ANGS, ANGS2, ANBTP, MOMTEN, PI)
      ELSEIF (AN) THEN
         DO J = 1, 4
            ANBTP(J) = ANBTP(J)/RDEG
         ENDDO
         CALL ANTPIN(ANBTP, ANGS, ANGS2, PTTP, MOMTEN, PI)
      ELSEIF (DSR) THEN
         DO J = 1, 3
            ANGS(J) = ANGS(J)/RDEG
         ENDDO
         CALL DSRIN(ANGS, ANBTP, ANGS2, PTTP, MOMTEN, PI)
      ENDIF
      DO I = 1, 3
         ANGS(I) = ANGS(I)*RDEG
         ANGS2(I) = ANGS2(I)*RDEG
         PTTP(I) = PTTP(I)*RDEG
         ANBTP(I) = ANBTP(I)*RDEG
      ENDDO
      ANBTP(4) = ANBTP(4)*RDEG
      ANBTP(5) = ANBTP(5)*RDEG
      ANBTP(6) = ANBTP(6)*RDEG
      PTTP(4) = PTTP(4)*RDEG
      IF (LUNIT1>0) THEN
!	  WRITE (LUNIT1,1) (ANGS(I),I=1,3)
!	  WRITE(LUNIT1,2)(ANGS2(I),I=1,3),'   Auxiliary Plane'
!	  WRITE (LUNIT1,3) (ANBTP(J),J=1,4)
!	  WRITE(LUNIT1,4) (ANBTP(J),J=5,6)
!	  WRITE(LUNIT1,5) PTTP
!	  WRITE(LUNIT1,6) MOMTEN
      ENDIF
      IF (LUNIT2>0) THEN
         IF (LUNIT2==5) THEN
!	    WRITE (*,1) (ANGS(I),I=1,3)
!	    WRITE(*,2)(ANGS2(I),I=1,3),'   Auxiliary Plane'
!	    WRITE (*,3) (ANBTP(J),J=1,4)
!	    WRITE(*,4) (ANBTP(J),J=5,6)
!	    WRITE(*,5) PTTP
!	    WRITE(*,6) MOMTEN
!	    WRITE (LUNIT2,1) (ANGS(I),I=1,3)
!	    WRITE(LUNIT2,2)(ANGS2(I),I=1,3),'   Auxiliary Plane'
!	    WRITE (LUNIT2,3) (ANBTP(J),J=1,4)
!	    WRITE(LUNIT2,4) (ANBTP(J),J=5,6)
!	    WRITE(LUNIT2,5) PTTP
!	    WRITE(LUNIT2,6) MOMTEN
         ENDIF
      ENDIF
      RETURN
!
99001 FORMAT (5X, 'Dip,Strike,Rake ', 3F9.2)
99002 FORMAT (5X, 'Dip,Strike,Rake ', 3F9.2, A)
99003 FORMAT (5X, 'Lower Hem. Trend, Plunge of A,N ', 4F9.2)
99004 FORMAT (5X, 'Lower Hem. Trend & Plunge of B ', 2F9.2)
99005 FORMAT (5X, 'Lower Hem. Trend, Plunge of P,T ', 4F9.2)
99006 FORMAT (5X, 'MRR =', F5.2, '  MTT =', F5.2, '  MPP =', F5.2, '  MRT =', F5.2, '  MRP =', F5.2, '  MTP =', F5.2)
      END SUBROUTINE FMREPS
 
 
      SUBROUTINE PRINTX(LINE)
      IMPLICIT NONE
      INTEGER J, KK, LENC
!+
!	SUBROUTINE PRINTX(LINE)
!  OUTPUTS A MESSAGE TO THE TERMINAL
!  PRINTX STARTS WITH A LINE FEED BUT DOES NOT END WITH A CARRIAGE RETURN
!  THE PRINT HEAD REMAINS AT THE END OF THE MESSAGE
!
!  IF THE MESSAGE LENGTH IS LESS THAN 40,
!	DOTS ARE INSERTED UP TO COL. 39
!	AND A COLON IS PUT IN COL. 40.
!
!  USE FOR CONVERSATIONAL INTERACTION
!			Alan Linde ... April 1980.
!	10 Sugust 1985:  Corrected a minor error for  strings > 40 bytes
!	20 June 1986:  Made it compatible with Fortran 77
!-
      CHARACTER*(*) LINE
      CHARACTER*60 BUF
      CHARACTER*2 COLON
      CHARACTER*1 DOT, DELIM
      DATA DELIM/'$'/, DOT/'.'/, COLON/': '/
      KK = LENC(LINE)   !  length minus right-hand blanks
      IF (LINE(KK:KK)==DELIM) KK = KK - 1
      IF (KK>58) KK = 59
      BUF(1:KK) = LINE(1:KK)
      IF (KK<49) THEN
         DO J = KK + 1, 49
            BUF(J:J) = DOT
         ENDDO
         KK = 49
      ENDIF
      BUF(KK:KK + 1) = COLON
      KK = KK + 2
      WRITE (*, '(1x,A,$)') BUF(1:KK)
      END SUBROUTINE PRINTX
 
 
      FUNCTION LENC(STRING)
      IMPLICIT NONE
      INTEGER J, LENC, LENGTH
!+
!	function lenc(string)
!
!	Returns length of character variable STRING excluding right-hand
!	  most blanks or nulls
!-
      CHARACTER*(*) STRING
      LENGTH = LEN(STRING)      ! total length
      IF (LENGTH==0) THEN
         LENC = 0
         RETURN
      ENDIF
      IF (ICHAR(STRING(LENGTH:LENGTH))==0) STRING(LENGTH:LENGTH) = ' '
      DO J = LENGTH, 1, -1
         LENC = J
         IF (STRING(J:J)/=' ' .AND. ICHAR(STRING(J:J))/=0) RETURN
      ENDDO
      LENC = 0
      END FUNCTION LENC
 
 
!+
 
 
!+
      SUBROUTINE IYESNO(MSG, IANS)
      IMPLICIT NONE
      INTEGER IANS, KK, NCHAR
!
!
!     PURPOSE:
!	     THIS LITTLE SUBROUTINE ASKS A QUESTION AND RETURNS A
!	     RESPONSE TO THAT QUESTION. THE ANSWER TO THE QUESTION
!	     MUST BE EITHER 'Y' FOR YES, 'N' FOR NO, OR NOTHING
!	     (i.e. simply hitting carrage return) FOR THE DEFAULT
!	     REPONSE TO THE QUESTION.
!
!     ON INPUT:
!	    MSG = BYTE STRING CONTAINING THE QUESTION
!
!     ON OUTPUT:
!	    IANS = THE LOGICAL REPONSE TO THE QUESTION (1 or 0)
!     EXTRA FEATURES:
!	    DEFAULT SITUATION IS:
!	    IF LAST 3 CHARACTERS IN 'MSG' ARE
!	  	     [Y]  OR  [N]
!	    THEN 'IANS' = 1   OR   0
!
!	    IF LAST 3 CHARACTERS ARE NOT ONE OF ABOVE PAIRS
!	    THEN 'IANS' = 0
!	    (i.e. default for no supplied default is N)
!	30 JULY 1989:  IF ENTERED CHARACTER IS A BLANK OR A TAB,
!	    TREATS AS A NULL ENTRY.
!       27 July 1993: Did input read through cstring so can have
!         comment lines
!-
      CHARACTER*1 DELIM/'$'/, CHARIN, BLANK/' '/
      CHARACTER*3 TEST, UCY, LCY
      CHARACTER*80 STRING_IN
      CHARACTER*(*) MSG
      DATA UCY/'[Y]'/, LCY/'[y]'/
      KK = LEN(MSG)
      IF (MSG(KK:KK)==DELIM) KK = KK - 1
      TEST = MSG(KK - 2:KK)
      CALL PRINTX(MSG)
      CALL CSTRING(STRING_IN, NCHAR)
      IF ((NCHAR>0) .AND. (STRING_IN(1:1)==BLANK)) NCHAR = 0
      IF (NCHAR/=0) THEN
         CHARIN = STRING_IN(1:1)
         IF (CHARIN==UCY(2:2) .OR. CHARIN==LCY(2:2)) THEN
            IANS = 1
         ELSE
            IANS = 0
         ENDIF
      ELSEIF ((TEST==UCY) .OR. (TEST==LCY)) THEN
         IANS = 1
      ELSE
         IANS = 0
      ENDIF
      END SUBROUTINE IYESNO
 
 
!+
      SUBROUTINE CSTRING(STRING, NSTRING)
      IMPLICIT NONE
      INTEGER LENC, NSTRING, NTAB
!
!	Input a character string with a read(*,'(A)') string
!	If first two characters are /* it will read the next entry
!	Tab is a delimiter.
!	Returns string and nstring, number of characters to tab.
!	string stars with first non-blank character.
!       25 May 2001.  Took out parameter statement for tab.
!-
      LOGICAL MORE
      CHARACTER*1 TAB
      CHARACTER*(*) STRING
!
      TAB = CHAR(9)
      MORE = .TRUE.
      DO WHILE (MORE)
         READ (*, '(A)') STRING
         NSTRING = LENC(STRING)
         MORE = (NSTRING>=2 .AND. STRING(1:2)=='/*')
      ENDDO
      IF (NSTRING>0) THEN
         NTAB = INDEX(STRING(1:NSTRING), TAB)
         IF (NTAB>0) NSTRING = NTAB - 1
      ENDIF
      END SUBROUTINE CSTRING
 
 
!+
      SUBROUTINE DSRIN(ANGS, ANBTP, ANGS2, PTTP, MOMTEN, PI)
      IMPLICIT NONE
      REAL A, ANBTP, ANGS, ANGS2, B, DIP, P, PI, PTTP, RAKE, SR2, STR, T
      INTEGER J
!
!	Calculates other representations of fault planes with
!		dip, strike and rake (A&R convinention) input.  All
!		angles are in radians.
!	22 July 1985:  Added moment tensor output (D&W convention)
!	               normalized to unit scalar moment
!-
      REAL N(3), MOMTEN(6)
      DIMENSION PTTP(4), ANGS(3), ANGS2(3), ANBTP(6), P(3), T(3), A(3), B(3)
      DATA SR2/0.707107/
      RAKE = ANGS(3)
      STR = ANGS(2)
      DIP = ANGS(1)
      A(1) = COS(RAKE)*COS(STR) + SIN(RAKE)*COS(DIP)*SIN(STR)
      A(2) = COS(RAKE)*SIN(STR) - SIN(RAKE)*COS(DIP)*COS(STR)
      A(3) = -SIN(RAKE)*SIN(DIP)
      N(1) = -SIN(STR)*SIN(DIP)
      N(2) = COS(STR)*SIN(DIP)
      N(3) = -COS(DIP)
      CALL V2TRPL(A, ANBTP(1), PI)
      CALL V2TRPL(N, ANBTP(3), PI)
      DO J = 1, 3
         T(J) = SR2*(A(J) + N(J))
         P(J) = SR2*(A(J) - N(J))
      ENDDO
      B(1) = P(2)*T(3) - P(3)*T(2)
      B(2) = P(3)*T(1) - P(1)*T(3)
      B(3) = P(1)*T(2) - P(2)*T(1)
      CALL V2TRPL(P, PTTP(1), PI)
      CALL V2TRPL(T, PTTP(3), PI)
      CALL V2TRPL(B, ANBTP(5), PI)
      CALL AN2DSR(N, A, ANGS2, PI)
      CALL AN2MOM(A, N, MOMTEN)
      END SUBROUTINE DSRIN
 
 
!+
      SUBROUTINE ANTPIN(ANBTP, ANGS, ANGS2, PTTP, MOMTEN, PI)
      IMPLICIT NONE
      REAL A, ANBTP, ANGS, ANGS2, B, P, PI, PTTP, SR2, T
      INTEGER J
!
!	Calculates other representations of fault planes with
!		trend and plunge of A and N as input.  All
!		angles are in radians.
!	22 July 1985:  Added moment tensor output.
!-
      REAL N(3), MOMTEN(6)
      DIMENSION PTTP(4), ANGS(3), ANGS2(3), ANBTP(6), P(3), T(3), A(3), B(3)
      DATA SR2/0.707107/
      CALL TRPL2V(ANBTP(1), A)
      CALL TRPL2V(ANBTP(3), N)
      DO J = 1, 3
         T(J) = SR2*(A(J) + N(J))
         P(J) = SR2*(A(J) - N(J))
      ENDDO
      B(1) = P(2)*T(3) - P(3)*T(2)
      B(2) = P(3)*T(1) - P(1)*T(3)
      B(3) = P(1)*T(2) - P(2)*T(1)
      CALL V2TRPL(P, PTTP(1), PI)
      CALL V2TRPL(T, PTTP(3), PI)
      CALL V2TRPL(B, ANBTP(5), PI)
      CALL AN2DSR(A, N, ANGS, PI)
      CALL AN2DSR(N, A, ANGS2, PI)
      CALL AN2MOM(A, N, MOMTEN)
      END SUBROUTINE ANTPIN
 
 
!+
      SUBROUTINE PTTPIN(PTTP, ANGS, ANGS2, ANBTP, MOMTEN, PI)
      IMPLICIT NONE
      REAL A, ANBTP, ANGS, ANGS2, B, P, PI, PTTP, SR2, T
      INTEGER J
!
!	Calculates other representations of fault planes with
!		trend and plunge of P and T as input.  All
!		angles are in radians.
!	22 July 1985:  Added moment tensor output
!-
      REAL N(3), MOMTEN(6)
      DIMENSION PTTP(4), ANGS(3), ANGS2(3), ANBTP(6), P(3), T(3), A(3), B(3)
      DATA SR2/0.707107/
      CALL TRPL2V(PTTP(1), P)
      CALL TRPL2V(PTTP(3), T)
      DO J = 1, 3
         A(J) = SR2*(P(J) + T(J))
         N(J) = SR2*(T(J) - P(J))
      ENDDO
      B(1) = P(2)*T(3) - P(3)*T(2)
      B(2) = P(3)*T(1) - P(1)*T(3)
      B(3) = P(1)*T(2) - P(2)*T(1)
      CALL V2TRPL(A, ANBTP(1), PI)
      CALL V2TRPL(N, ANBTP(3), PI)
      CALL V2TRPL(B, ANBTP(5), PI)
      CALL AN2DSR(A, N, ANGS, PI)
      CALL AN2DSR(N, A, ANGS2, PI)
      CALL AN2MOM(A, N, MOMTEN)
      END SUBROUTINE PTTPIN
 
 
!+
      SUBROUTINE AN2MOM(A, N, MOMTEN)
      IMPLICIT NONE
      INTEGER J
!
!	Starting with the A and N axis, calculates the elements
!	  of the moment tensor with unit scalar moment.
!	  Convention used is that of Dziewonski & Woodhouse
!	  (JGR 88, 3247-3271, 1983) and Aki & Richards (p 118)
!	24 September 1985: If an element is < 0.000001 (ABS), set to zero
!-
      REAL*4 A(3), N(3), MOMTEN(6)
!	      Moment tensor components:  M(I,j) = A(I)*N(J)+A(J)*N(I)
      MOMTEN(1) = 2.0*A(3)*N(3)         !  MRR = M(3,3)
      MOMTEN(2) = 2.0*A(1)*N(1)         !  MTT = M(1,1)
      MOMTEN(3) = 2.0*A(2)*N(2)         !  MPP = M(2,2)
      MOMTEN(4) = A(1)*N(3) + A(3)*N(1) !  MRT = M(1,3)
      MOMTEN(5) = -A(2)*N(3) - A(3)*N(2)
                                        !  MRP = -M(2,3)
      MOMTEN(6) = -A(2)*N(1) - A(1)*N(2)
                                        !  MTP = -M(2,1)
      DO J = 1, 6
         IF (ABS(MOMTEN(J))<0.000001) MOMTEN(J) = 0.0
      ENDDO
      END SUBROUTINE AN2MOM
 
 
!+
!	SUBROUTINE AN2DSR(A,N,ANGS,PI)
!
!	Calculates dip, strike and rake (ANGS) - A&R convention,
!		from A and N.
!	12 January 2000:  Fixed a divide by zero when angs(1) .eq. 0
!-
      SUBROUTINE AN2DSR(A, N, ANGS, PI)
      IMPLICIT NONE
      REAL A1, A2, PI, RR
      REAL N(3), A(3), ANGS(3)
      IF (N(3)== - 1.0) THEN
         ANGS(2) = ATAN2(A(2), A(1))
         ANGS(1) = 0.0
      ELSE
         ANGS(2) = ATAN2( - N(1), N(2))
         IF (N(3)==0.0) THEN
            ANGS(1) = 0.5*PI
         ELSEIF (ABS(SIN(ANGS(2)))>=0.1) THEN
            ANGS(1) = ATAN2( - N(1)/SIN(ANGS(2)), -N(3))
         ELSE
            ANGS(1) = ATAN2(N(2)/COS(ANGS(2)), -N(3))
         ENDIF
      ENDIF
      A1 = A(1)*COS(ANGS(2)) + A(2)*SIN(ANGS(2))
      IF (ABS(A1)<0.0001) A1 = 0.0
      IF (A(3)==0.0) THEN
         A2 = A(1)*SIN(ANGS(2)) - A(2)*COS(ANGS(2))
         IF (ABS(A2)<0.0001) A2 = 0.0
         IF (ABS(SIN(2*ANGS(2)))>=0.0001) THEN
            ANGS(3) = ATAN2(A2/SIN(2*ANGS(2)), A1)
         ELSEIF (ABS(SIN(ANGS(2)))>=0.0001) THEN
            RR = A(2)/SIN(ANGS(2))
            IF (RR>1) RR = 1.
            IF (RR< - 1.) RR = -1.
            ANGS(3) = ACOS(RR)
         ELSE
            IF (A1>1) A1 = 1.
            IF (A1< - 1.) A1 = -1.
            ANGS(3) = ACOS(A1)
         ENDIF
      ELSEIF (ANGS(1)/=0.0) THEN
         ANGS(3) = ATAN2( - A(3)/SIN(ANGS(1)), A1)
      ELSE
         ANGS(3) = ATAN2( - 1000000.0*A(3), A1)
      ENDIF
      IF (ANGS(1)<0.0) THEN
         ANGS(1) = ANGS(1) + PI
         ANGS(3) = PI - ANGS(3)
         IF (ANGS(3)>PI) ANGS(3) = ANGS(3) - 2*PI
      ENDIF
      IF (ANGS(1)>0.5*PI) THEN
         ANGS(1) = PI - ANGS(1)
         ANGS(2) = ANGS(2) + PI
         ANGS(3) = -ANGS(3)
         IF (ANGS(2)>=2*PI) ANGS(2) = ANGS(2) - 2*PI
      ENDIF
      IF (ANGS(2)<0.0) ANGS(2) = ANGS(2) + 2.0*PI
      END SUBROUTINE AN2DSR
 
 
!+
      SUBROUTINE V2TRPL(XYZ, TRPL, PI)
      IMPLICIT NONE
      REAL C, PI, S, TRPL, XYZ
      INTEGER J
!
!	Transforms from XYZ components of a unit vector to
!	  the trend and plunge for the vector.
!	Trend is the azimuth (clockwise from north looking down)
!	Plunge is the downward dip measured from the horizontal.
!	All angles in radians
!	X is north, Y is east, Z is down
!	If the component of Z is negative (up), the plunge,TRPL(2),
!	  is replaced by its negative and the trend, TRPL(1),
!	  Is changed by PI.
!	The trend is returned between 0 and 2*PI, the plunge
!	  between 0 and PI/2.
!	12 January 2000: If xyz(3) = -1.0, make the trend PI.  Made
!	  consistency in the roundoff -- all are now 0.0001
!-
      DIMENSION XYZ(3), TRPL(2)
      DO J = 1, 3
         IF (ABS(XYZ(J))<=0.0001) XYZ(J) = 0.0
         IF (ABS(ABS(XYZ(J)) - 1.0)<0.0001) XYZ(J) = XYZ(J)/ABS(XYZ(J))
      ENDDO
      IF (ABS(XYZ(3))==1.0) THEN
!
!	plunge is 90 degrees
!
         IF (XYZ(3)<0.0) THEN
            TRPL(1) = PI
         ELSE
            TRPL(1) = 0.0
         ENDIF
         TRPL(2) = 0.5*PI
         RETURN
      ENDIF
      IF (ABS(XYZ(1))>=0.0001) THEN
         TRPL(1) = ATAN2(XYZ(2), XYZ(1))
      ELSEIF (XYZ(2)>0.0) THEN
         TRPL(1) = PI/2.
      ELSEIF (XYZ(2)<0.0) THEN
         TRPL(1) = 3.0*PI/2.0
      ELSE
         TRPL(1) = 0.0
      ENDIF
      C = COS(TRPL(1))
      S = SIN(TRPL(1))
      IF (ABS(C)>=0.1) TRPL(2) = ATAN2(XYZ(3), XYZ(1)/C)
      IF (ABS(C)<0.1) TRPL(2) = ATAN2(XYZ(3), XYZ(2)/S)
      IF (TRPL(2)<0.0) THEN
         TRPL(2) = -TRPL(2)
         TRPL(1) = TRPL(1) - PI
      ENDIF
      IF (TRPL(1)<0.0) TRPL(1) = TRPL(1) + 2.0*PI
      END SUBROUTINE V2TRPL
 
 
!+
!	SUBROUTINE TRPL2V(TRPL,XYZ)
!
!	Transforms to XYZ components of a unit vector from
!		the trend and plunge for the vector.
!	Trend is the azimuth (clockwise from north looking down)
!	Plunge is the downward dip measured from the horizontal.
!	All angles in radians
!	X is north, Y is east, Z is down
!-
      SUBROUTINE TRPL2V(TRPL, XYZ)
      IMPLICIT NONE
      INTEGER J
      REAL TRPL, XYZ
      DIMENSION XYZ(3), TRPL(2)
      XYZ(1) = COS(TRPL(1))*COS(TRPL(2))
      XYZ(2) = SIN(TRPL(1))*COS(TRPL(2))
      XYZ(3) = SIN(TRPL(2))
      DO J = 1, 3
         IF (ABS(XYZ(J))<0.0001) XYZ(J) = 0.0
         IF (ABS(ABS(XYZ(J)) - 1.0)<0.0001) XYZ(J) = XYZ(J)/ABS(XYZ(J))
      ENDDO
      END SUBROUTINE TRPL2V
!+
      SUBROUTINE MT_IN(PTTP, PI)
      IMPLICIT NONE
      REAL A, DEVMOM, EPS, EPS1, EPS2, PI, PTTP, U, W, XYZ
      INTEGER J
!
!	eigenvalues/vectors using EISPACK routines from www.netlib.no
!	Much of code adapted from Jost/Herrmann mteig.f and mtdec.f
!	Uses original EISPACK routines for TRED2 and IMTQL2, not NR
!	Also includes subroutine eig, which calls TRED2 and IMTQL2.
!-
      DIMENSION A(3, 3), U(3, 3), W(3), PTTP(4), XYZ(3)
      REAL*4 MRR, MTT, MPP, MRT, MRP, MTP, ISOTROP
      WRITE (*, *) 'Input MRR MTT MPP MRT MRP MTP (free format)'
      READ (*, *) MRR, MTT, MPP, MRT, MRP, MTP
!      write(2,*) '  Input is moment tensor (Dziewonski convention)'
!     write(2,*) 'MRR MTT MPP MRT MRP MTP'
!      write(2,'(1p6g11.4)') MRR, MTT, MPP, MRT, MRP, MTP
!
!	Convention is X north, Y east, Z down
!
      A(3, 3) = MRR
      A(1, 1) = MTT
      A(2, 2) = MPP
      A(1, 3) = MRT
      A(3, 1) = A(1, 3)
      A(2, 3) = -MRP
      A(3, 2) = A(2, 3)
      A(1, 2) = -MTP
      A(2, 1) = A(1, 2)
      CALL EIG(A, U, W)
!
!	Ordering returned is from smallest to highest (P, B, T for DC)
!
!      write(*,*) ' '
!     write(*,*) 'EIGENVALUES                EIGENVECTORS'
!     do j=1,3
!       write(*,'(1pg11.4,2x,0p5f11.4)') W(j),(U(i,j),i=1,3)
!     end do
      ISOTROP = 0.0
      DO J = 1, 3
         ISOTROP = ISOTROP + W(J)
      ENDDO
      DEVMOM = 0.0
      DO J = 1, 3
         W(J) = W(J) - ISOTROP/3.0
         DEVMOM = DEVMOM + W(J)*W(J)
      ENDDO
      DEVMOM = SQRT(0.5*DEVMOM)
      IF (DEVMOM<0.001*ISOTROP) DEVMOM = 0.0
      WRITE (*, *) ' '
      WRITE (*, '(a,1pg11.4,a,g11.4)') '  Trace of moment tensor = ', ISOTROP, '  Deviatoric tensor moment = ', DEVMOM
      WRITE (2, *) ' '
      WRITE (2, '(a,1pg11.4,a,g11.4)') '  Trace of moment tensor = ', ISOTROP, '  Deviatoric tensor moment = ', DEVMOM
      IF (DEVMOM==0.0) THEN
         WRITE (2, *) 'Exiting because purely isotropic source'
         WRITE (*, *) 'Exiting because purely isotropic source'
         STOP
      ENDIF
      WRITE (2, '(a,1p3g11.4)') '  Deviatoric moment tensor eigenvalues:', (W(J), J = 1, 3)
      WRITE (*, '(a,1p3g11.4)') '  Deviatoric moment tensor eigenvalues:', (W(J), J = 1, 3)
!---- Dziewonski, Chou, Woodhouse,  JGR 1981 2825-2852
!---- eps=0 pure double couple
!---- eps=0.5 pure CLVD
      EPS = ABS(W(2)/AMAX1( - W(1), W(3)))
      EPS1 = EPS*200.0
      EPS2 = 100.0 - EPS1
      WRITE (*, '(a)') '  EPSILON    % OF CLVD     % OF DC'
      WRITE (*, '(f9.4,2f12.4)') EPS, EPS1, EPS2
      WRITE (2, '(a)') '  EPSILON    % OF CLVD     % OF DC'
      WRITE (2, '(f9.4,2f12.4)') EPS, EPS1, EPS2
      WRITE (2, *) ' '
      IF (EPS>=0.25) THEN
         WRITE (2, *) ' Exiting because less than 50% double couple'
         WRITE (*, *) ' Exiting because less than 50% double couple'
         STOP
      ENDIF
!
!	Get trend and plunge for P
!
      DO J = 1, 3
         XYZ(J) = U(J, 1)
      ENDDO
      CALL V2TRPL(XYZ, PTTP(1), PI)
!
!	Get trend and plunge for T
!
      DO J = 1, 3
         XYZ(J) = U(J, 3)
      ENDDO
      CALL V2TRPL(XYZ, PTTP(3), PI)
!     do j=1,4
!        pttp(j) = rdeg*pttp(j)
!      end do
!      write(*,*) '  '
!      write(*,*) '  Trend and Plunge of P and T'
!     write(*,'(4f11.4)') pttp
      END SUBROUTINE MT_IN
      SUBROUTINE EIG(A, U, W)
      IMPLICIT NONE
      REAL A, U, W, WORK
      INTEGER I, IERR, J, N, NP
      DIMENSION A(3, 3), U(3, 3), W(3), WORK(3)
      NP = 3
      DO I = 1, NP
         DO J = 1, 3
            U(I, J) = A(I, J)
         ENDDO
      ENDDO
      N = 3
      CALL TRED2(NP, N, A, W, WORK, U)
      CALL IMTQL2(NP, N, W, WORK, U, IERR)
!
!	This system has P, B, T as a right-hand coordinate system.
!	I prefer P, T, B
!
      DO J = 1, 3
         U(J, 1) = -U(J, 1)
      ENDDO
      END SUBROUTINE EIG
!---------------
      SUBROUTINE TRED2(NM, N, A, D, E, Z)
      IMPLICIT NONE
!
      INTEGER I, J, K, L, N, NM
      REAL A(NM, N), D(N), E(N), Z(NM, N)
      REAL F, G, H, HH, SCALE
!
!     this subroutine is a translation of the algol procedure tred2,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix to a
!     symmetric tridiagonal matrix using and accumulating
!     orthogonal similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        z contains the orthogonal transformation matrix
!          produced in the reduction.
!
!        a and z may coincide.  if distinct, a is unaltered.
!
!     Questions and comments should be directed to Alan K. Cline,
!     Pleasant Valley Software, 8603 Altus Cove, Austin, TX 78759.
!     Electronic mail to cline@cs.utexas.edu.
!
!     this version dated january 1989. (for the IBM 3090vf)
!
!     ------------------------------------------------------------------
!
!?      call xuflow(0)
      DO I = 1, N
         DO J = I, N
            Z(J, I) = A(J, I)
         ENDDO
      ENDDO
!
      DO I = 1, N
         D(I) = A(N, I)
      ENDDO
!
      DO I = N, 2, -1
         L = I - 1
         H = 0.0E0
         SCALE = 0.0E0
         IF (L>=2) THEN
!     .......... scale row (algol tol then not needed) ..........
            DO K = 1, L
               SCALE = SCALE + ABS(D(K))
            ENDDO
!
            IF (SCALE/=0.0E0) THEN
!
               DO K = 1, L
                  D(K) = D(K)/SCALE
                  H = H + D(K)*D(K)
               ENDDO
!
               F = D(L)
               G = -SIGN(SQRT(H), F)
               E(I) = SCALE*G
               H = H - F*G
               D(L) = F - G
!     .......... form a*u ..........
               DO J = 1, L
                  E(J) = 0.0E0
               ENDDO
!
               DO J = 1, L
                  F = D(J)
                  Z(J, I) = F
                  G = E(J) + Z(J, J)*F
!
                  DO K = J + 1, L
                     G = G + Z(K, J)*D(K)
                     E(K) = E(K) + Z(K, J)*F
                  ENDDO
!
                  E(J) = G
               ENDDO
!     .......... form p ..........
               F = 0.0E0
!
               DO J = 1, L
                  E(J) = E(J)/H
                  F = F + E(J)*D(J)
               ENDDO
!
               HH = -F/(H + H)
!     .......... form q ..........
               DO J = 1, L
                  E(J) = E(J) + HH*D(J)
               ENDDO
!     .......... form reduced a ..........
               DO J = 1, L
                  F = -D(J)
                  G = -E(J)
!
                  DO K = J, L
                     Z(K, J) = Z(K, J) + F*E(K) + G*D(K)
                  ENDDO
!
                  D(J) = Z(L, J)
                  Z(I, J) = 0.0E0
               ENDDO
               GOTO 50
            ENDIF
         ENDIF
         E(I) = D(L)
!
!"    ( ignore recrdeps
         DO J = 1, L
            D(J) = Z(L, J)
            Z(I, J) = 0.0E0
            Z(J, I) = 0.0E0
!
         ENDDO
!
 50      D(I) = H
      ENDDO
!     .......... accumulation of transformation matrices ..........
      DO I = 2, N
         L = I - 1
         Z(N, L) = Z(L, L)
         Z(L, L) = 1.0E0
         H = D(I)
         IF (H/=0.0E0) THEN
!
            DO K = 1, L
               D(K) = Z(K, I)/H
            ENDDO
!     ( ignore recrdeps
!     ( prefer vector
            DO J = 1, L
               G = 0.0E0
!
               DO K = 1, L
                  G = G + Z(K, I)*Z(K, J)
               ENDDO
!
               G = -G
!
               DO K = 1, L
                  Z(K, J) = Z(K, J) + G*D(K)
               ENDDO
            ENDDO
         ENDIF
!
         DO K = 1, L
            Z(K, I) = 0.0E0
         ENDDO
!
      ENDDO
!
!"    ( prefer vector
      DO I = 1, N
         D(I) = Z(N, I)
         Z(N, I) = 0.0E0
      ENDDO
!
      Z(N, N) = 1.0E0
      E(1) = 0.0E0
      END SUBROUTINE TRED2
      SUBROUTINE IMTQL2(NM, N, D, E, Z, IERR)
      IMPLICIT NONE
!
      INTEGER I, J, K, L, M, N, NM, IERR
      REAL D(N), E(N), Z(NM, N)
      REAL B, C, F, G, P, R, S, TST1, TST2
!
!     this subroutine is a translation of the algol procedure imtql2,
!     num. math. 12, 377-383(1968) by martin and wilkinson,
!     as modified in num. math. 15, 450(1970) by dubrulle.
!     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the implicit ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1.
!
!        e has been destroyed.
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     Questions and comments should be directed to Alan K. Cline,
!     Pleasant Valley Software, 8603 Altus Cove, Austin, TX 78759.
!     Electronic mail to cline@cs.utexas.edu.
!
!     this version dated january 1989. (for the IBM 3090vf)
!
!     ------------------------------------------------------------------
!
!?      call xuflow(0)
      IERR = 0
      IF (N/=1) THEN
!
         DO I = 2, N
            E(I - 1) = E(I)
         ENDDO
!
         E(N) = 0.0E0
!
         DO L = 1, N
            J = 0
!     .......... look for small sub-diagonal element ..........
 20         DO M = L, N - 1
               TST1 = ABS(D(M)) + ABS(D(M + 1))
               TST2 = TST1 + ABS(E(M))
               IF (TST2==TST1) EXIT
            ENDDO
!
            P = D(L)
            IF (M/=L) THEN
               IF (J==30) GOTO 100
               J = J + 1
!     .......... form shift ..........
               G = (D(L + 1) - P)/(2.0E0*E(L))
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
! *      r = pythag(g,1.0d0)
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
               IF (ABS(G)<=1.0E0) THEN
                  R = SQRT(1.0E0 + G*G)
               ELSE
                  R = G*SQRT(1.0E0 + (1.0E0/G)**2)
               ENDIF
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
               G = D(M) - P + E(L)/(G + SIGN(R, G))
               S = 1.0E0
               C = 1.0E0
               P = 0.0E0
!     .......... for i=m-1 step -1 until l do -- ..........
               DO I = M - 1, L, -1
                  F = S*E(I)
                  B = C*E(I)
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
! *         r = pythag(f,g)
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
                  IF (ABS(F)>=ABS(G)) THEN
                     R = ABS(F)*SQRT(1.0E0 + (G/F)**2)
                  ELSEIF (G/=0.0E0) THEN
                     R = ABS(G)*SQRT((F/G)**2 + 1.0E0)
                  ELSE
                     R = ABS(F)
                  ENDIF
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
                  E(I + 1) = R
                  IF (R==0.0E0) THEN
!     .......... recover from underflow ..........
                     D(I + 1) = D(I + 1) - P
                     E(M) = 0.0E0
                     GOTO 20
                  ENDIF
                  S = F/R
                  C = G/R
                  G = D(I + 1) - P
                  R = (D(I) - G)*S + 2.0E0*C*B
                  P = S*R
                  D(I + 1) = G + P
                  G = C*R - B
!     .......... form vector ..........
                  DO K = 1, N
                     F = Z(K, I + 1)
                     Z(K, I + 1) = S*Z(K, I) + C*F
                     Z(K, I) = C*Z(K, I) - S*F
                  ENDDO
!
               ENDDO
!
               D(L) = D(L) - P
               E(L) = G
               E(M) = 0.0E0
               GOTO 20
            ENDIF
         ENDDO
!     .......... order eigenvalues and eigenvectors ..........
         DO I = 1, N - 1
            K = I
            P = D(I)
!
            DO J = I + 1, N
               IF (D(J)<P) THEN
                  K = J
                  P = D(J)
               ENDIF
            ENDDO
!
            D(K) = D(I)
            D(I) = P
!
            DO J = 1, N
               P = Z(J, I)
               Z(J, I) = Z(J, K)
               Z(J, K) = P
            ENDDO
!
!
         ENDDO
      ENDIF
      GOTO 99999
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 100  IERR = L
99999 END SUBROUTINE IMTQL2
