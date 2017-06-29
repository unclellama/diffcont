!   Last-modified: 13 Jan 2011 02:21:01 AM

!**********************************************************************

module run_amoeba
use LightCurves
use TransferFunctions
use ModelParameters
use BeretMode
implicit none
interface
    function func(pin)
    REAL(kind=8), DIMENSION(:) :: pin
    REAL(kind=8) :: func
    end function
    SUBROUTINE amoeba(p,y,ftol,func,iter)
    INTEGER(kind=4) :: iter
    REAL(kind=8) :: ftol
    REAL(kind=8), DIMENSION(:) :: y
    REAL(kind=8), DIMENSION(:,:) :: p
    REAL(kind=8) :: func
    EXTERNAL func
    END SUBROUTINE
end interface
!----------------------------------------------------------------------
CONTAINS

subroutine init_plike(plike1st)
implicit none
REAL(kind=8) :: plike1st
REAL(kind=8),allocatable :: p1sttry(:)
INTEGER(kind=4) :: i,j,k
penalty = 0.0D0
ALLOCATE(p1sttry(nvar))
k = 0
DO i=1,ndim
    IF (ivary(i).eq.1) THEN
        k = k + 1
        p1sttry(k) = psave(i)
    ENDIF
ENDDO
plike1st = -func(p1sttry)
DEALLOCATE(p1sttry)
end subroutine

subroutine amoeba_minimize(pbest)
implicit none
REAL(kind=8) :: pbest
REAL(kind=8),allocatable :: p(:,:),y(:),ptry(:)
INTEGER(kind=4) :: i,j,k,it,iter
REAL(kind=8) ::  ftol

ALLOCATE(p(nvar+1,nvar),y(nvar+1),ptry(nvar))
j = 0
DO i=1,ndim
    IF (ivary(i).eq.1) THEN
        j      = j+1
        ptry(j)= psave(i)
    ENDIF
ENDDO
! initialize the minimization by generating the vertices
! of the amoeba -- initial point and THEN nvar additional
! ones sequentially making the variables 10% larger
DO i=1,nvar+1
  DO j=1,nvar
    p(i,j) = ptry(j)
  ENDDO
  IF (i.ne.1) THEN
    it       = i-1
    if(ptry(it).ne.0.0D0) THEN
        p(i,it) = 1.1D0*ptry(it)
    ELSE
! assign 1.0 to zero initial value (arbitrary and only 
! for optimization over alpha)
        p(i,it) = 0.02
    ENDIF
  ENDIF
ENDDO
DO i=1,nvar+1
  DO j=1,nvar
    ptry(j) = p(i,j)
  ENDDO
! calculate the minus of likelihood
  y(i) = func(ptry)
ENDDO

! run the minimiziation
ftol   = 0.000001D0
CALL amoeba(p,y,ftol,func,iter)
WRITE(*,'(a32,1x,i4,a11)')'AMOEBA DONE IN:',iter,'ITERATIONS'
! find the best vertex of the amoeba (minimum value of y)
pbest = 1.D32
DO i=1,nvar+1
  IF (y(i).le.pbest) THEN
    pbest = y(i)
    DO j=1,nvar
      ptry(j) = p(i,j)
    ENDDO
  ENDIF
ENDDO
! func is the thing minimized by amoeba, which means
! it is the NEGATIVE of the likelihood
pbest   = -func(ptry)
! now unfold the minimized variables back into the 
! full variable array
j = 0
DO i=1,ndim
  IF (ivary(i).eq.1) THEN
    j          = j + 1
! update psave with the new best model found by amoeba
    psave(i)   = ptry(j)
  ENDIF
ENDDO
end subroutine amoeba_minimize
!----------------------------------------------------------------------
end module run_amoeba

!**********************************************************************
