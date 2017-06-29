! Last-modified: 24 Apr 2011 06:09:17 PM

!***********************************************************************
! FUN_FUNC
!**********************************************************************
! func = -likelihood = function to be minmized
!***********************************************************************
FUNCTION func(pin)
use LightCurves
use TransferFunctions
use ModelParameters
use BeretMode
implicit none
REAL(kind=8), DIMENSION(:) :: pin
REAL(kind=8),allocatable :: p(:)
REAL(kind=8) :: func
!----------------------------------------------------------------------
REAL(kind=8), allocatable :: lcovar(:,:)
REAL(kind=8) :: vuse,tuse,vusel,tusel,argv,argt,salpha
REAL(kind=8),allocatable :: scale(:),swid(:),slag(:)
REAL(kind=8) :: plikeb,chi2b
INTEGER(kind=4) :: i,j,k
!-----------------------------------------------------------------------
interface
    subroutine Copy_ModelPara(conv_mode,n_curve,var,tau,alpha,lag,wid,scale,p)
    REAL(kind=8) var,tau,alpha
    REAL(kind=8),DIMENSION(:) :: lag,wid,scale,p
    INTEGER(kind=4) :: conv_mode,n_curve
    end subroutine
    subroutine fastlikelin(sigma,tau,slag,swid,scale,salpha,lcovar,plike,chi2)
    REAL(kind=8) :: sigma,tau,variance,plike,chi2
    REAL(kind=8),DIMENSION(:) :: scale,slag,swid
    REAL(kind=8),DIMENSION(:,:) :: lcovar
    REAL(kind=8) :: salpha
    end subroutine
end interface
!-----------------------------------------------------------------------
! initialize the parameters to defaults
! set the linear parameters to work out separate means for each light curve
ALLOCATE(lcovar(nlin,nlin),p(ndim))
ALLOCATE(scale(ncurve),swid(ncurve),slag(ncurve))
!-----------------------------------------------------------------------
! set parameters for whether they are being minimized and
! come from the input vector pin or whether they are being
! held constant and come from the saved vector psave
k = 0
DO i=1,ndim
  IF (ivary(i).eq.1) THEN
    k = k + 1
    p(i) = pin(k)
  ELSE
! psave comes with the common block
    p(i) = psave(i)
  ENDIF
ENDDO
!-----------------------------------------------------------------------
do i=1,ncurve
  slag(i)  = 0.0D0
  swid(i)  = 0.0D0
  scale(i) = 1.0D0
enddo
! unroll these into the variables being passed to fastlikelin
CALL Copy_ModelPara(convmode,ncurve,vuse,tuse,salpha,slag,swid,scale,p)
!-----------------------------------------------------------------------
! the minimizer sometimes wanders into unphysical regimes that can cause
! crashes -- avoid the crash by setting the parameter to something
! reasonable, add a huge penality to the function
IF (tuse.lt.0.0D0) THEN
  print*,'WARNING: NEGATIVE TIME SCALE ',tuse,' SETTING TO &
  ONE-TENTH OF A DAY AND PENALIZING!'
  tuse    = 0.1D0
  penalty = penalty + 1.0D12
ELSE IF (vuse.lt.0.0D0) THEN
  print*,'WARNING: NEGATIVE DISPERSION ',vuse,' SETTING TO &
  ONE-TENTH OF A DAY AND PENALIZING!'
  vuse    = 0.1D0
  penalty = penalty + 1.0D12
ELSE
  penalty = 0.0D0
ENDIF
!-----------------------------------------------------------------------
! compute the likelihood of the model
CALL fastlikelin(vuse,tuse,slag,swid,scale,salpha,lcovar,plikeb,chi2b)
IF (iposdef .ne. 0) THEN 
    print*,'WARNING: MATRIX IS NOT POSITIVE DEFINITE(func)!'
    penalty = penalty + 1.0D12
ENDIF
!-----------------------------------------------------------------------
! if multiple light curves, use the "Gaussian" prior on the tau/sigma values 
! derived from the continuum fits alone

! if just the continuum, use log priors on sigma=vuse and 
! for tau use log(tmed/tau) when tau > tmed = median lightcurve point spacing
!             log(tau/tmed)          < tmed
IF (convmode.ne.SINGLE) THEN
  tusel  = LOG10(tuse)
  vusel  = LOG10(vuse)
  IF (tusel.lt.taucen) THEN
    argt = (tusel-taucen)/tauscatlow
  ELSE
    argt = (tusel-taucen)/tauscathig
  ENDIF
  IF (vusel.lt.sigcen) THEN
    argv = (vusel-sigcen)/sigscatlow
  ELSE
    argv = (vusel-sigcen)/sigscathig
  ENDIF
  plikeb = plikeb - 0.5*(argt*argt+argv*argv)
ELSE
  IF (tuse.ge. tmed) THEN
    plikeb = plikeb - log(tuse/tmed) - log(vuse)
  ELSE IF (tuse.lt. tmed) THEN
    plikeb = plikeb + log(tuse/tmed) - log(vuse)
  ENDIF
ENDIF
!-----------------------------------------------------------------------
! save the chi^2 value into the common block
chisave= chi2b
func   = -plikeb + penalty
!-----------------------------------------------------------------------
IF (chi2b.lt.0.0) THEN
   print*,'WARINING: NEGATIVE CHI^2',chi2b,'REJECT OR STOP!'
   print*,'PARAMETERS SIGMA,TAU:',vuse,tuse
   do j=1,ncurve
      print*,'FUNC: OBJ',j,'SCALE,LAG,WIDTH,ALPHA',&
      scale(j),slag(j),swid(j),salpha
   enddo
   RETURN
ENDIF
!-----------------------------------------------------------------------
RETURN
END
!***********************************************************************
