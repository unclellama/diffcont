!Last-modified: 24 Apr 2011 05:18:59 PM

!----------------------------------------------------------------------

Module MCMC
use TransferFunctions
use ModelParameters
implicit none
INTERFACE
   FUNCTION func(x)
   IMPLICIT NONE
   REAL(kind=8), DIMENSION(:), INTENT(IN) :: x
   REAL(kind=8) :: func
   END FUNCTION func
END INTERFACE

CONTAINS

SUBROUTINE Run_MCMC(idomcmc,mcmcfile)
implicit none
INTEGER(kind=4) :: idomcmc,iseed,iopen,ikeep,T(8)
INTEGER(kind=4) :: i,j,k,ntrymax,ntry,nacc,nwrite,nskip
INTEGER(kind=4),parameter :: nburn = 1000, nmcmc = 10000
REAL(kind=8),allocatable :: dmcmc(:),parold(:),parnew(:),ptry(:)
REAL(kind=8) :: plast,pnew,ratio,ran2,gasdev,facc
CHARACTER(len=*) :: mcmcfile
CHARACTER(len=2) :: vfe_str
EXTERNAL ran2,gasdev

IF(idomcmc .eq. 1)THEN
PRINT*,'STARTING MCMC...'
CALL DATE_AND_TIME(VALUES = T)
iseed = T(1)+70*(T(2)+12*(T(3)+31*(T(5)+23*(T(6)+59*T(7)))))
IF (MOD(iseed,2).EQ.0) iseed = iseed-1
iseed = 0 - iseed
! variable format express
WRITE(vfe_str,'(i2)')ndim
! assume at least 20% acceptance rate
ntrymax = 5 * nmcmc
!----------------------------------------------------------------------
ALLOCATE(dmcmc(ndim),parold(ndim),parnew(ndim),ptry(nvar))
! logarithmic steps for tau and sigma
dmcmc(1) = 0.025D0
dmcmc(2) = 0.15D0
! linear steps for lag, width, and scale
DO i=2,ncurve
  dmcmc((i-1)*3)   = 0.4D0
  dmcmc((i-1)*3+1) = 0.025D0
  dmcmc((i-1)*3+2) = 0.1D0
ENDDO
! linear steps for alpha
IF(convmode .eq. BERET)dmcmc(6) = 0.03D0
!----------------------------------------------------------------------
ntry = 0
nacc = 0
nwrite = 0
iopen = 0
nskip = 0
TRIAL: DO i=1,ntrymax
    IF (i.eq.ntrymax) THEN
        PRINT*,'WARNING: CHAIN ABOUT TO DIE BECAUSE IT HITS NTRYMAX'
    ENDIF
    IF ((ntry.EQ.0).AND.(iopen.EQ.0)) THEN
        parold = psave
        plast  = plikesave
        PRINT*,'START BURN-IN'
    ENDIF
   IF ((ntry.GE.nburn).AND.(iopen.EQ.0)) THEN
        CLOSE(UNIT=25)
        nwrite = 0
        OPEN(UNIT=30,FILE=mcmcfile,FORM='FORMATTED',STATUS='UNKNOWN')
        iopen = 1
        PRINT*,'OPENING OUTPUT MCMC FILE, SET IOPEN TO 1 '
    ENDIF
    IF (nwrite.EQ.nmcmc) THEN
        CLOSE(UNIT=30)
        PRINT*, 'FINISHING WRITING MCMC...QUIT'
        EXIT
    ENDIF
! logarithmic wander in tau and sigma
    parnew(1) = parold(1)*exp(dmcmc(1)*(-1.0D0+2.0D0*ran2(iseed)))
    parnew(2) = parold(2)*exp(dmcmc(2)*(-1.0D0+2.0D0*ran2(iseed)))
! linear wander in lag, amplitude, width and alpha
    DO j=3,ndim
        parnew(j) = parold(j) + dmcmc(j)*gasdev(iseed)
    ENDDO
! keep parameters if ivary == 0
    DO j=1,ndim
      IF (ivary(j).eq.0) parnew(j) = psave(j)
    ENDDO
! copy the relevant variables for func, which only takes varying parameters
    k = 0
    DO j=1,ndim
      IF (ivary(j).eq.1) THEN
        k          = k + 1
        ptry(k) = parnew(j)
      ENDIF
    ENDDO
! evaluate likelihood
    pnew = -func(ptry)
    ntry = ntry + 1
    ikeep= 0
    IF(iposdef .ne. 0) THEN
        ikeep = 0
        nskip = nskip + 1
        PRINT*,'WARNING: SKIP AN ILL-DEFINED COVARIANCE MATRIX'
        IF(nskip .gt. 100)THEN
            PRINT*,'ERROR: TOO MANY STEPS SKIPPED FOR BAD MATRICES'
            IF(iopen .eq. 1)THEN
                CLOSE(30)
            ELSE
                CLOSE(25)
            ENDIF
            EXIT TRIAL
        ENDIF
        CONTINUE
    ENDIF
! Metropolis-Hastings algorithm
    IF (pnew.ge.plast) THEN
      ikeep = 1
    ELSE
      ratio= exp(max(-72.0D0,pnew-plast))
      IF (ran2(iseed).le.ratio) THEN
        ikeep = 1
      ELSE
        CONTINUE
      ENDIF
    ENDIF

    IF (ikeep.EQ.0) THEN
!        CYCLE TRIAL
         CONTINUE
    ELSE
        nacc = nacc + 1
        IF (100*int(float(i)/100.D0).eq.i) THEN
            WRITE(*,'(a32,1x,i6,a20,i6)')'ACCEPTED:',nacc,'AFTER TRYING ',ntry
            IF(iopen .eq. 1)THEN
                CALL FLUSH(30)
            ELSE
                CALL FLUSH(25)
            ENDIF
        ENDIF
        IF(iopen .eq. 0)THEN
            WRITE(25,'(2(1x,i5),4(1x,f13.6),$)')ntry,nacc,pnew,plast,&
            LOG10(parnew(1)),LOG10(abs(parnew(2)))
            WRITE(25,'('//vfe_str//'(1x,f13.6),$)')(parnew(j),j=1,ndim)
            WRITE(25,'(f13.5)')chisave
        ELSE IF(iopen .eq. 1)THEN
            WRITE(30,'(2(1x,i5),4(1x,f13.6),$)')ntry,nacc,pnew,plast,&
            LOG10(parnew(1)),LOG10(abs(parnew(2)))
            WRITE(30,'('//vfe_str//'(1x,f13.6),$)')(parnew(j),j=1,ndim)
            WRITE(30,'(f13.5)')chisave
            nwrite = nwrite + 1
        ENDIF
! update the starting point
        plast = pnew
        parold= parnew
    ENDIF

    IF (ntry.GT.100) THEN
        facc = float(nacc)/float(ntry)
        IF(facc .lt. 0.2D0) THEN
            ntry = 0
            nacc = 0
            PRINT*,'RESETTING MCMC RANGE DOWNWARDS '
            dmcmc = 0.5D0*dmcmc
            IF (iopen.EQ.1) THEN
                CLOSE(UNIT=30)
                iopen = 0
            ELSE
                CLOSE(UNIT=25)
            ENDIF
        ELSE IF (facc.GT.0.5D0) THEN
            ntry = 0
            nacc = 0
            PRINT*,'RESETTING MCMC RANGE UPWARDS '
            dmcmc = 2.0D0*dmcmc
            IF (iopen.EQ.1) THEN
                CLOSE(UNIT=30)
                iopen = 0
            ELSE
                CLOSE(UNIT=25)
            ENDIF
        ENDIF
    ENDIF
ENDDO TRIAL


DEALLOCATE(dmcmc,parold,parnew,ptry)
ENDIF
END SUBROUTINE

!----------------------------------------------------------------------
End Module
!----------------------------------------------------------------------

FUNCTION ran2(idum)
implicit none
INTEGER(kind=4) :: idum
REAL(kind=8) :: ran2
INTEGER(kind=4),PARAMETER :: IM1=2147483563,IM2=2147483399,&
IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
NTAB=32,NDIV=1+IMM1/NTAB
REAL(kind=8),PARAMETER :: AM=1.D0/float(IM1),EPS=1.2D-7,RNMX=1.D0-EPS
INTEGER(kind=4) :: j,k
INTEGER(kind=4), SAVE :: idum2=123456789,iy=0
INTEGER(kind=4),DIMENSION(NTAB),SAVE :: iv = 0

if (idum.le.0) then
  idum=max(-idum,1)
  idum2=idum
  do j=NTAB+8,1,-1
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      if (j.le.NTAB) iv(j)=idum
  enddo
  iy=iv(1)
endif
k=idum/IQ1
idum=IA1*(idum-k*IQ1)-k*IR1
if (idum.lt.0) idum=idum+IM1
k=idum2/IQ2
idum2=IA2*(idum2-k*IQ2)-k*IR2
if (idum2.lt.0) idum2=idum2+IM2
j=1+iy/NDIV
iy=iv(j)-idum2
iv(j)=idum
if(iy.lt.1)iy=iy+IMM1
ran2=min(AM*iy,RNMX)
return
END FUNCTION

!----------------------------------------------------------------------

FUNCTION gasdev(idum)
implicit none
INTEGER(kind=4) :: idum
REAL(kind=8) :: gasdev
INTEGER(kind=4),SAVE :: iset = 0
REAL(kind=8),SAVE :: gset
REAL(kind=8) :: fac,rsq,v1,v2,ran2

if (iset.eq.0) then
    do
        v1=2.*ran2(idum)-1.
        v2=2.*ran2(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.lt.1.D0.and.rsq.gt.0.D0)exit 
    enddo
    fac=sqrt(-2.D0*log(rsq)/rsq)
    gset=v1*fac
    gasdev=v2*fac
    iset=1
else
    gasdev=gset
    iset=0
endif
return
END FUNCTION

!----------------------------------------------------------------------
