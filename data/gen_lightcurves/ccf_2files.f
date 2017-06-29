      PROGRAM crossc
*     
*     THIS PROGRAM CALCULATES THE CCF USING NO PADDING.  ALSO,
*     THE NORMALIZATION IS CALCULATED BY USING ONLY THE REAL DATA
*     POINTS WHICH ARE BEING CORRELATED.  THE USER HAS THE OPTION
*     OF SELECTING A MAXIMUM SHIFT TO BE CORRELATED.
*
*     THIS PROGRAM IS THE DRIVER FOR C.M. GASKELL'S SUBROUTINE CCOR,
*     WHICH PERFORMS CROSS-CORRELATIONS ON SERIES OF DATA
*
*     THIS IS A "SPECIAL PERFORMANCE" VERSION OF THE DRIVER. IT
*     ALLOWS USE OF SEPARATE FILES FOR THE CONTINUUM AND LINE
*     LIGHT CURVES, AND ALLOWS INTERPOLATION IN ONLY ONE ARRAY
*     AS WELL AS IN BOTH
*
*     22-JAN-90  BMP
*
*     THIS VERSION WAS WRITTEN BY RUSS WHITE IN AUTUMN 1992, AND
*     ADAPTED TO WORK WITH PETERSON'S PROGRAMS 
*     31-DEC-92 BMP
*
*     increased maxpix from 500 to 2000 to handle NGC 5548 data
*     24-Jan-94 BMP
*
*     experimental version to test variable centroid threshold
*     03-Mar-94 BMP
*
*     corrected centroid calculation
*     25-Jul-95 BMP
*    
      PARAMETER (MAXPIX=5000,MAXCCF=40100)
      COMMON /CCFAR/ CCF(MAXCCF),INDEXA(MAXPIX),INDEXB(MAXPIX)
      DIMENSION CONT(MAXPIX), FLUX(MAXPIX), TIME1(MAXPIX)
      DIMENSION TIME2(MAXPIX), ARRAY(MAXPIX) 
      DIMENSION AOUT(2*MAXPIX), BOUT(2*MAXPIX)
      CHARACTER*80 DUMSTR, NAME, NAME2
c      BYTE SRES
      INTEGER*4 WINDOW, N2, ICODE, FINA, FINB, NPTS, LASTSH, SHI, RANG
      REAL*8 TUNIT, LEVEL, SPAN, DTMIN, SHFT, JOUT(2*MAXPIX)
      REAL*8 FRONT, TAIL,LPEAK, TEMP, VAL,POS,JDEL
      REAL*8 VMIN
      DATA LEVEL/0.80/, WINDOW/30/
c      DATA IYES/'Y'/, JYES/'y'/
      CHARACTER*1 ANSW,SRES

***** PROMPT FOR DATA FILE INFORMATION *****

      WRITE(6,10)
 10   FORMAT(' Cross correlation of data from 1 or 2 seperate '/ 
     &     '  light curve files. This version assumes that the'/
     &     '  data are in 5-column light-curve format and that'/
     &     '  the light curves are non-stationary, i.e., '/
     &     '  mean and variance are recomputed for each lag.')
      WRITE(6,20)
 20   FORMAT(' Enter number of files to be read:',$)
      READ(5,*) NFILE
      IF(NFILE.NE.2) NFILE=1
      IF(NFILE.EQ.1) THEN
         WRITE(6,50)
 50      FORMAT(' Enter name of data file for cross-correlation:',$)
      ELSE
         WRITE(6,51)
 51      FORMAT(' Enter name of first data file for'/ 
     &        '  cross-correlation:',$)
      ENDIF
      READ(5,'(A)') NAME   
      OPEN(UNIT=1,FILE=NAME,STATUS='OLD',IOSTAT=IERR)
      IF(NFILE.EQ.2) THEN
         WRITE(6,70)
 70      FORMAT(' Select CONTINUUM (=1, DEFAULT) or LINE (=2)'/
     &        '  from data file:',$)
         READ(5,*) K1
         IF(K1.NE.2) K1=1
      ENDIF

***** READ DATA FROM DATA FILE *****
 
      NPTS=0
      DO 230 I = 1, MAXPIX
            READ(1,*,IOSTAT=IERR)TIME1(I),CONT(I),CDUM
c,FLUX(I),FDUM
            IF (IERR.LT.0) GOTO 250
cm            WRITE(6,200) TIME1(I),CONT(I),CDUM,
cm     &           FLUX(I),FDUM
            
cm 200        FORMAT(1X,1P,5(1X,E12.4),0P)      
         NPTS=NPTS+1
         IF (IERR.NE.0) GOTO 990
 230  CONTINUE
 250  IF (NPTS.EQ.MAXPIX) THEN
         READ(1,'(A)',IOSTAT=IERR) DUMSTR
         IF (IERR.EQ.0) THEN
            WRITE(6,255) MAXPIX
 255        FORMAT( I3,' observational times have been read '/
     &           '   from the data file.  There is more data in '/
     &           '   the file, but the arrays storing this data are'/
     &           '   full. Should I continue (Y/N)?:',$)
            READ(5,'(A)') ANSW
            IF(ANSW.NE.'Y'.AND.ANSW.NE.'y') THEN
               GOTO 995   
            ENDIF
         ENDIF
      ENDIF
      CLOSE(UNIT=1)
      
      IF(NFILE.EQ.1) THEN
         DO 258 I=1,NPTS
            JOUT(I)=TIME1(I)
 258     CONTINUE
         WRITE(6,260) NPTS
 260     FORMAT(' Number of Observations = ',I4)
      ELSE
         WRITE(6,261)
 261     FORMAT(' Enter name of second data file for cross'/
     &        '  correlation:',$)
         READ(5,'(A)') NAME2
         OPEN(UNIT=1,FILE=NAME2,STATUS='OLD',IOSTAT=IERR)
         WRITE(6,70)
         READ(5,*) K2
         IF(K2.NE.2) K2=1

***** READ DATA FROM SECOND FILE *****

         N2=0
         DO 270 I=1, MAXPIX
            READ(1,*,IOSTAT=IERR) TIME2(I),FVAL,FDUM
            IF (IERR.LT.0) GOTO 275
cm            WRITE(6,200) TIME2(I),CVAL,CDUM,
cm     &           FVAL,FDUM
         N2=N2+1
         IF (IERR.NE.0) GOTO 990
         IF(K2.EQ.2) THEN
            ARRAY(N2)=FVAL
         ELSE
            ARRAY(N2)=CVAL
         ENDIF
 270  CONTINUE
 275  CLOSE(UNIT=1)

***** MERGE DATA FROM FILES ONE AND TWO *****
     
         INDA=1
         INDB=1
         CDUM=0.
         KOUNT=0
 281     KOUNT=KOUNT+1
         JDEL=TIME1(INDA)-TIME2(INDB)
         IF(JDEL.GT.0) THEN
            JOUT(KOUNT)=TIME2(INDB)
            AOUT(KOUNT)=0.
            BOUT(KOUNT)=ARRAY(INDB)
            INDB=INDB+1
            IF(INDB.GT.N2) GO TO 285
            GO TO 281
         ENDIF
         IF(JDEL.EQ.0) THEN
            JOUT(KOUNT)=TIME1(INDA)
            IF(K1.EQ.1) THEN
               AOUT(KOUNT)=CONT(INDA)
            ELSE
               AOUT(KOUNT)=FLUX(INDA)
            ENDIF
            BOUT(KOUNT)=ARRAY(INDB)
            INDA=INDA+1
            INDB=INDB+1
            IF(INDA.GT.NPTS) THEN
               IF(INDB.GT.N2) GO TO 310
               GO TO 300
            ELSE
               IF(INDB.GT.N2) GO TO 285
               GO TO 281
            ENDIF
         ELSE
            JOUT(KOUNT)=TIME1(INDA)
            IF(K1.EQ.1) THEN
               AOUT(KOUNT)=CONT(INDA)
            ELSE
               AOUT(KOUNT)=FLUX(INDA)
            ENDIF
            BOUT(KOUNT)=0.
            INDA=INDA+1
            IF(INDA.GT.NPTS) GOTO 300
            GO TO 281
         ENDIF
 285     KOUNT=KOUNT+1
         JOUT(KOUNT)=TIME1(INDA)
         IF(K1.EQ.1) THEN
            AOUT(KOUNT)=CONT(INDA)
         ELSE
            AOUT(KOUNT)=FLUX(INDA)
         ENDIF
         BOUT(KOUNT)=0.
         INDA=INDA+1
         IF(INDA.GT.NPTS) GOTO 310
         GO TO 285
 300     KOUNT=KOUNT+1
         JOUT(KOUNT)=TIME2(INDB)
         AOUT(KOUNT)=0.
         BOUT(KOUNT)=ARRAY(INDB)
         INDB=INDB+1
         IF(INDB.LE.N2) GOTO 300
 310     INDA=NPTS
         NPTS=KOUNT
         INDB=INDA+N2-KOUNT         
         DO 320 I=1,NPTS
            CONT(I)=AOUT(I)
            FLUX(I)=BOUT(I)
            IF (I.EQ.1 ) PRINT*,'Merged data from both files:'
cm            WRITE(6,200) REAL(JOUT(I)),CONT(I),FLUX(I)
 320        CONTINUE
         WRITE(6,330) INDA,N2,INDB
 330     FORMAT('  No. of nights in first array = ',I4,
     &        '  No. of nights in second array = ',I4,
     &        '  No. of nights overlap = ',I4)
      ENDIF

***** SELECT AN INTERPOLATION UNIT (days per bin) *****

      IF(NFILE.EQ.1)THEN
         SPAN=TIME1(NPTS)-TIME1(1)
      ELSE
         SPAN=JOUT(NPTS)-JOUT(1)
      ENDIF
      DTMIN=(SPAN)/95000
      WRITE(6,335) DTMIN
 335  FORMAT(' Select an interpolating time unit that is',
     &     '  a divisor of 1 and is larger than ',F5.3,': ',$)
      READ(5,*)TUNIT

***** SET FIRST CONTINUUM = FIRST POINT *****

      FINA = 0
      FINB = 0
      DO 342 I = 1, NPTS
         IF(CONT(I).NE.0.) THEN
            FINA=FINA+1
            INDEXA(FINA)=I
         ENDIF
         IF(FLUX(I).NE.0.) THEN
            FINB=FINB+1
            INDEXB(FINB)=I
         ENDIF
 342  CONTINUE

***** DETERMINE MAXIMUM RANGE FOR SHIFT *****

      IF (JOUT(INDEXA(1)).GE.JOUT(INDEXB(1))) THEN
         FRONT = JOUT(INDEXA(1))
      ELSE
         FRONT = JOUT(INDEXB(1))
      ENDIF
      IF (JOUT(INDEXA(FINA)).LE.JOUT(INDEXB(FINB))) THEN
         TAIL = JOUT(INDEXA(FINA))
      ELSE
         TAIL = JOUT(INDEXB(FINB))
      ENDIF
      RANG = INT((TAIL-FRONT)/2)
      WRITE(6,350) RANG
 350  FORMAT (' The default range for calculation is ',I4,' days'/
     & ' Do you wish to change this (Y/N)?: ',$)
      READ(5,'(A)') SRES
      IF (SRES.EQ.'Y'.OR.SRES.EQ.'y') THEN
         WRITE(6,360)
 360     FORMAT(' Enter maximum shift: ',$)
         READ(5,*) SHI
         IF(SHI.LT.RANG) THEN
            LASTSH = INT(SHI/TUNIT)
            ELSE
            LASTSH = INT(RANG/TUNIT)
         ENDIF
         ELSE
         LASTSH=INT(RANG/TUNIT)
      ENDIF

***** SELECT THRESHOLD FOR COMPUTATION OF CENTROID
      WRITE(6,370)
 370  FORMAT(' Enter threshold (fraction of CCF maximum value,'/
     & '  default = 0.8) for centroid calculation: ',$)
      READ(5,*) RCENT
      IF(RCENT.LT.0..OR.RCENT.GE.1.) THEN
         WRITE(6,371)
 371     FORMAT(' *** value out of range. Default value assumed ***')
         RCENT=0.8
         ENDIF

***** SELECT TYPE OF INTERPOLATION *****

      WRITE(6,380)
 380  FORMAT(' Select interpolation scheme:'/
     &     '   0 = both arrays (DEFAULT)'/
     &     '   1 = first array only'/
     &     '   2 = second array only'/
     &     '   Enter option:',$)
      READ(5,*) ICODE
      IF(ICODE.NE.1.AND.ICODE.NE.2) ICODE=0

***** CALL THE CROSS-CORRELATION SUBROUTINE *****
      CALL  CCOR(NPTS,JOUT,LEVEL,WINDOW,POS,
     &     VAL,VMIN,ICODE,TUNIT,LPEAK,FINA,FINB,
     &     CONT,FLUX,LASTSH)

***** PRINT RESULTS OF CROSS-CORRELATION *****

      WRITE(6,385) INT(LASTSH*TUNIT), INT(LASTSH*TUNIT)
 385  FORMAT('Cross-correlation coefficients have been calculated for'/
     &     ' shifts from -',I4,' to ',I4,' days.')
      WRITE(6,390) LPEAK,VAL,POS
      WRITE(43,*)NPTS,LPEAK
 390  FORMAT(' Peak lag = ',F9.3/
     &     ' Maximum correlation coefficient = ',F6.3/
     &     ' First moment peak =',F9.3)

***** WRITE CROSS-CORRELATION FUNCTION TO PLOT FILE *****

      OPEN(UNIT=10,FILE='ccf.dat',STATUS='unknown',IOSTAT=IERR)
      DO 450 I = 1, LASTSH*2+1
         SHFT= (real(I)-(LASTSH+1))*TUNIT
         WRITE(10,400) SHFT, CCF(I)
 400     FORMAT(1X,F11.3,1X,F8.4)
 450  CONTINUE
      CLOSE(UNIT=10)

***** LOCATE PEAK IN CCF
      MAXR=1
      RMAX=CCF(1)
      DO 500 I=2,2*LASTSH+1
         IF(CCF(I).GT.RMAX) THEN
            RMAX=CCF(I)
            MAXR=I
            ENDIF
 500  CONTINUE

***** LOCATE RCENT*MAXIMUM VALUE ON LEFT SIDE OF PEAK

      I=1
      RLIMIT=RCENT*RMAX
 510  IND=MAXR-I
      IF(IND.LT.1) GO TO 520
      IF(CCF(IND).GT.RLIMIT) THEN
         I=I+1
         GO TO 510
         ELSE
         DX=(RLIMIT-CCF(IND))/(CCF(IND+1)-CCF(IND))
         TL=(real(IND-LASTSH-1)+DX)*TUNIT
         INDEXL=IND+1
         GO TO 550
         ENDIF
 520  WRITE(6,530)
 530  FORMAT(' *** Unable to compute error in peak from CCF ***')
      IF(RCENT.GT.0.5) GO TO 605
      STOP

**** LOCATE HALF MAXIMUM VALUE ON RIGHT SIDE OF PEAK

 550  I=1
 560  IND=MAXR+I
      IF(IND.GT.(LASTSH*2+1)) GO TO 520
      IF(CCF(IND).GT.RLIMIT) THEN
         I=I+1
         GO TO 560
         ELSE
         DX=(CCF(IND-1)-RLIMIT)/(CCF(IND-1)-CCF(IND))
         TR=(real(IND-LASTSH-2)+DX)*TUNIT
         INDEXR=IND-1
         ENDIF
      TMID=0.5*(TR+TL)
      HWHM=0.5*(TR-TL)
      FWHM=2.*HWHM
      TEMP=NPTS - 2
      GPERR=(0.75*HWHM)/(1. +RMAX*DSQRT(TEMP))
      WRITE(6,600) FWHM,TMID,GPERR
 600  FORMAT(' Full-width half-maximum (FWHM) = ',F8.3/
     &       ' Midpoint of half-maximum points = ',F8.3/
     &       ' Gaskell-Peterson error in peak  = ', F8.3)

***** COMPUTE CENTROID ABOVE LEVEL RCENT
***** LOCATE THRESHOLD ON LEFT SIDE OF PEAK
      IF(RCENT.EQ.0.5) GO TO 670 
 605  I=1
      RTHRES=RCENT*RMAX
 610  IND=MAXR-I
      IF(IND.LT.1) GO TO 620
      IF(CCF(IND).GT.RTHRES) THEN
         I=I+1
         GO TO 610
         ELSE
         INDEXL=IND+1
         GO TO 650
         ENDIF
 620  WRITE(6,630)
 630  FORMAT(' *** Unable to compute centroid of CCF peak   ***'/
     &       ' *** for this threshold. Try larger threshold ***')
      STOP

**** LOCATE THRESHOLD ON RIGHT SIDE OF PEAK

 650  I=1
 660  IND=MAXR+I
      IF(IND.GT.(LASTSH*2+1)) GO TO 620
      IF(CCF(IND).GT.RTHRES) THEN
         I=I+1
         GO TO 660
         ELSE
         INDEXR=IND-1
         ENDIF
 670  SUMN=0.
      SUMD=0.
      DO 680 I = INDEXL,INDEXR
         TIME=(real(I)-LASTSH-1)*TUNIT
         SUMN=SUMN + CCF(I)*TIME
         SUMD=SUMD + CCF(I)
 680     CONTINUE
      CENTRO=SUMN/SUMD
      WRITE(6,700) RCENT,CENTRO
c      WRITE(22,*)INDA,CENTRO
      WRITE(22,*)NPTS,CENTRO
 700  FORMAT(' Centroid at ',f4.2,' maximum = ',f8.3)
      STOP

***** WARNING OF UNFAMILIAR DATA IN DATA FILE *****

 990  PRINT*,'The data file contains non-standard data.'

 995  STOP
      END
      SUBROUTINE CCOR(NPOINT,JOUT,LEVEL,WINDOW,
     &     POSITI,PEAKVA,MINVAL,ICODE,TUNIT,LPEAK,
     &     FINA,FINB,A,B,LASTSH)
     
*     MARTIN GASKELL'S CROSS-CORRELATION PROGRAM,
*     MODIFIED FOR SEPARATE INTERPOLATIONS IN A AND B

*     NPOINT = NUMBER OF OBSERVATIONS
*     JOUT = REALARRAY[1..NPOINT] OF OBSEVATION DATES        
*     A = REALARRAY[1..NPOINT] OF FIRST THING TO BE CORRELATED
*     B = ..DITTO..FOR SECOND THING TO BE CORRELATED
*     LEVEL [REAL] = THRESHOLD FOR FINDING FIRST MOMENT PEAK
*          (RECOMMEND 0.80). MUST NOT BE >< 1.0
*     WINDOW [INTEGER] = WIDTH OF PEAK IN DAYS ( 30 IS USED )
*     POSITI = FIRST MOMENT PEAK (RETURNED). IF A OF -9999 IS
*          RETURNED, THEN NO PEAK HAS BEEN FOUND.
*     PEAKVA  [REAL]=
*     MINVAL [REAL] =
*     ICODE = CROSS-CORRELATION CODE:
*          0 (CORRELATE IN A AND B AND AVERAGE)
*          1 (CORRELATE IN A ONLY)
*          2 (CORRELATE IN B ONLY)
*     TUNIT = TIME UNIT USED FOR INTERPOLATION 
*
*     increased maxpix from 500 to 2000    BMP 24-Jan-94      
      PARAMETER (MAXCCF=40100,MAXPIX=5000)
      COMMON/CCFAR/CCF(MAXCCF),INDEXA(MAXPIX),INDEXB(MAXPIX)
 
      INTEGER*4  FINA, FINB, FRSTSH, BN, AN, WTA, WTB,
     &     ABEG, BBEG, AEND, BEND, A1, B1, AF1, BF1,
     &     LASTSH, LEFT, MAXSHI, MINLAG, ICODE,
     &     PEAKLA, NPOINT, RANGE, RIGHT, SHIFT, WINDOW

      REAL*8  ACORRE, AXSUM, AYSUM, AXYSUM, BCORRE, BXSUM, BYSUM,
     &     BXYSUM, LEVEL, MEANAI, MEANAS, MEANBI, MEANBS,
     &     MINVAL, PEAKVA, TUNIT, JOUT(NPOINT),
     &     POSITI, THRESH, BOTTOM, TOP, LPEAK 

      DIMENSION A(NPOINT)
      DIMENSION B(NPOINT)

***** CALCULATION OF CCF LOOP *****
     
      MINVAL = +1.1
      PEAKVA = -1.1
      FRSTSH = -1*LASTSH

***** FIND CORRELATION COEFFICIENTS *****

      DO 600 SHIFT = frstsh,lastsh

         AXSUM = 0.0
         AYSUM = 0.0
         AXYSUM = 0.0
         ACORRE = 0.0
         BXSUM = 0.0
         BYSUM = 0.0
         BXYSUM = 0.0
         BCORRE = 0.0

***** CALCULATION OF THE BEGINING AND ENDING OF THE CC SHIFT ***** 

         B1=1
 660     IF((JOUT(INDEXB(B1))-SHIFT*TUNIT).GE.JOUT(INDEXA(1))) THEN
            BBEG = B1
         ELSE
            B1=B1+1
            GOTO 660
         ENDIF
         BF1=FINB
 665     IF((JOUT(INDEXB(BF1))-SHIFT*TUNIT).LE.
     &        JOUT(INDEXA(FINA))) THEN
            BEND = BF1
         ELSE
            BF1=BF1-1
            GOTO 665
         ENDIF

         A1=1
 680     IF((JOUT(INDEXA(A1))+SHIFT*TUNIT).GE.JOUT(INDEXB(1))) THEN
            ABEG = A1
         ELSE
            A1=A1+1
            GOTO 680
         ENDIF
         AF1=FINA
 685     IF((JOUT(INDEXA(AF1))+SHIFT*TUNIT).LE.
     &        JOUT(INDEXB(FINB))) THEN
            AEND = AF1
         ELSE
            AF1=AF1-1
            GOTO 685
         ENDIF

***** DETERMINE THE RELATIVE WEIGHT OF EACH CORRELATION *****

         WTA = 0
         WTB = 0
         WTA = AEND - ABEG + 1
         WTB = BEND - BBEG + 1
         
***** CALCULATE ACORRE  *****

         MEANBI = 0.0
         MEANAS = 0.0
         DO 698 BN = BBEG, BEND
            MEANAS = MEANAS + AINTER(BN,SHIFT,JOUT,A,FINA,NPOINT,TUNIT)
            MEANBI = MEANBI + B(INDEXB(BN))
 698     CONTINUE
         MEANAS = MEANAS/WTB
         MEANBI = MEANBI/WTB

         IF (ICODE.NE.2) THEN
            DO 700 BN = BBEG,BEND
               AXSUM = AXSUM + (AINTER(BN,SHIFT,JOUT,A,FINA,
     &              NPOINT,TUNIT) -MEANAS)**2
               AYSUM = AYSUM + (B(INDEXB(BN))-MEANBI)**2
               AXYSUM = AXYSUM + (AINTER(BN,SHIFT,JOUT,A,FINA,
     &              NPOINT,TUNIT) -MEANAS)*(B(INDEXB(BN))-MEANBI)
 700        CONTINUE
            IF(AXSUM*AYSUM.LE.0.0) GOTO 11
            ACORRE = AXYSUM/SQRT(AXSUM*AYSUM)
         ENDIF

***** CALCULATE BCORRE *****

         MEANAI = 0.0
         MEANBS = 0.0
         DO 798 AN = ABEG, AEND
            MEANBS = MEANBS+BINTER(AN,SHIFT,JOUT,B,FINB,NPOINT,TUNIT)
            MEANAI = MEANAI+A(INDEXA(AN))
 798     CONTINUE
         MEANBS = MEANBS/WTA
         MEANAI = MEANAI/WTA

         IF (ICODE.NE.1) THEN
            DO 800 AN = ABEG,AEND
               BXSUM = BXSUM + (BINTER(AN,SHIFT,JOUT,B,FINB,NPOINT,
     &              TUNIT)-MEANBS)**2
               BYSUM = BYSUM + (A(INDEXA(AN))-MEANAI)**2
               BXYSUM = BXYSUM + (BINTER(AN,SHIFT,JOUT,B,FINB,NPOINT,
     &              TUNIT)-MEANBS)*(A(INDEXA(AN))-MEANAI)
 800        CONTINUE
            IF (BXSUM*BYSUM.LE.0.0) GOTO 11
            BCORRE = BXYSUM/SQRT(BXSUM*BYSUM)
         ENDIF

***** COMPUTE CCF BASED ON INTERPOLATION CODE *****
          
         ADJ=INT(LASTSH)+1
         IF(ICODE.EQ.0) THEN
     
*     CCF IS WEIGHTED AVERAGE OF TWO CORRELATIONS

            CCF(SHIFT+ADJ) = (ACORRE*WTB+BCORRE*WTA)/(WTA+WTB)
            ELSE
            IF(ICODE.EQ.1) THEN
     
*     CCF IS BASED ON INTERPOLATION IN A ONLY
     
               CCF(SHIFT+ADJ) = ACORRE
            ELSE

*     CCF IS BASED ON INTERPOLATION IN B ONLY
     
               CCF(SHIFT+ADJ) = BCORRE
            ENDIF
         ENDIF
 11      CONTINUE

***** FIND PEAK IN CCF *****
     
         IF (CCF(SHIFT+ADJ).GE.PEAKVA) THEN 
            PEAKLA = SHIFT
            PEAKVA = CCF(SHIFT+ADJ)
         ENDIF
         IF (CCF(SHIFT+ADJ).LE.MINVAL) THEN
            MINLAG = SHIFT
            MINVAL = CCF(SHIFT+ADJ)
         ENDIF
 600  CONTINUE

***** DECIDE ON CRITERION FOR PEAK ( LEVEL=0.80, WINDOW=30 ) *****

      MAXSHI = INT(LASTSH)
      TOP = 0.0
      BOTTOM =0.0
      THRESH = LEVEL*(PEAKVA-MINVAL) + MINVAL 
      RANGE = INT(WINDOW/(TUNIT*2))
      IF(PEAKLA-RANGE.LT.-MAXSHI) RANGE = MAXSHI + PEAKLA
      IF(PEAKLA+RANGE.GT.MAXSHI) RANGE = MAXSHI - PEAKLA
      LEFT = PEAKLA - RANGE + ADJ
      RIGHT = PEAKLA + RANGE + ADJ
      IF(CCF(LEFT).GT.THRESH) THRESH = CCF(LEFT)
      IF(CCF(RIGHT).GT.THRESH) THRESH = CCF(RIGHT)

***** FIND FIRST MOMENT PEAK *****

      DO 14, SHIFT = LEFT, RIGHT
         IF(CCF(SHIFT).LE.THRESH) GO TO 15
         TOP = TOP + (CCF(SHIFT)-THRESH)*(SHIFT-ADJ)
         Bottom = BOTTOM + (CCF(SHIFT)-THRESH)
 15      CONTINUE
 14   CONTINUE
      POSITI = -9999

***** RETURN -9999 IF NO PEAK *****
     
      IF (BOTTOM.EQ.0.0) GO TO 10
      POSITI = TOP/BOTTOM
 10   CONTINUE
      LPEAK = real(PEAKLA)*TUNIT
      POSITI = POSITI*TUNIT
      RETURN
      END


      REAL*4 FUNCTION AINTER(BN,SHIFT,JOUT,A,FINA,NPOINT,TUNIT)

      PARAMETER(MAXCCF=40100,MAXPIX=5000)
      COMMON/CCFAR/CCF(MAXCCF),INDEXA(MAXPIX),INDEXB(MAXPIX)
      INTEGER*4 FINA, NPOINT, BN,  SHIFT
      REAL*8 LASTA, NEXTA, LASTTM, NEXTTM, JOUT(NPOINT),TUNIT
      DIMENSION A(NPOINT)
      
      K=1
      SHFTTM = JOUT(INDEXB(BN))-SHIFT*TUNIT
      
 710  IF(SHFTTM.EQ.JOUT(INDEXA(K))) THEN
         AINTER = A(INDEXA(K))
         GOTO 720
      ENDIF

      NEXTTM = JOUT(INDEXA(K+1))
      IF(SHFTTM.LT.NEXTTM) THEN
         LASTTM = JOUT(INDEXA(K))
         LASTA = A(INDEXA(K))
         NEXTA = A(INDEXA(K+1))
         AINTER = LASTA + (NEXTA - LASTA)*(SHFTTM-LASTTM)/
     &        (NEXTTM-LASTTM)
         GOTO 720
      ELSE
         K=K+1
         GOTO 710
      ENDIF
  
 720  RETURN
      END

      REAL*4 FUNCTION BINTER(AN,SHIFT,JOUT,B,FINB,NPOINT,TUNIT)
     
      PARAMETER (MAXPIX=5000, MAXCCF=40100)
      COMMON/CCFAR/CCF(MAXCCF), INDEXA(MAXPIX), INDEXB(MAXPIX)
      INTEGER*4 FINB, NPOINT, AN, SHIFT
      REAL*8 LASTB, NEXTB, LASTTM, NEXTTM, JOUT(NPOINT),TUNIT
      DIMENSION  B(NPOINT)

      K = 1
      SHFTTM = JOUT(INDEXA(AN))+SHIFT*TUNIT

 730  IF (SHFTTM.EQ.JOUT(INDEXB(K))) THEN
         BINTER = B(INDEXB(K))
         GOTO 740
      ENDIF

      NEXTTM = JOUT(INDEXB(K+1))
      IF (SHFTTM.LT.NEXTTM) THEN
         LASTTM = JOUT(INDEXB(K))
         LASTB = B(INDEXB(K))
         NEXTB = B(INDEXB(K+1))
         BINTER = LASTB + (NEXTB-LASTB)*(SHFTTM-LASTTM)/
     &        (NEXTTM-LASTTM)
         GOTO 740
      ELSE
         K=K+1
         GOTO 730
      ENDIF

 740  RETURN
      END

