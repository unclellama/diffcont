! Last-modified: 24 Apr 2011 06:10:25 PM



Module HighRes
use ModelParameters
use LightCurves
use BeretMode
implicit none
interface
    subroutine Chol_GetSol(INFO,mat,vecb,nsize,vecsol)
    integer(kind=4) :: nsize,INFO
    real(kind=8) :: mat(nsize,nsize),vecb(nsize),vecsol(nsize)
    end subroutine
    subroutine LU_GetSol(INFO,mat,vecb,nsize,vecsol)
    integer(kind=4) :: nsize,INFO
    real(kind=8) :: mat(nsize,nsize),vecb(nsize),vecsol(nsize)
    end subroutine
end interface
REAL(kind=8),allocatable :: jwant(:),mwant(:),ewant(:),magvar(:)
INTEGER(kind=4),allocatable :: iwant(:)
INTEGER(kind=4) :: nwant

contains

SUBROUTINE Get_HighResLF(icurve,sigma,tau,slag,swid,scale,salpha)
implicit none
REAL(kind=8) :: sigma,tau,salpha
REAL(kind=8),DIMENSION(:) ::  slag,swid,scale
INTEGER(kind=4) :: icurve,ihires,ialloc,i,j,k
REAL(kind=8) :: trange,tstart,tend,dt

IF(icurve .eq. 1)THEN
    nwant    = 10 * npt
    ALLOCATE(jwant(nwant),mwant(nwant),ewant(nwant),magvar(npt),iwant(nwant),STAT=ialloc)
    IF(ialloc .ne. 0) STOP 'ERROR: ARRAY ALLOCATION FAILED IN Get_HighResLF'
    ! presumably the light curve was sorted
    trange   = jd(npt)-jd(1)
    tstart   = jd(1)  - 0.2*trange
    tend     = jd(npt)+ 0.2*trange
    dt       = (tend-tstart)/float(nwant-1)
!    PRINT*,' DATA  START, END ',jd(1),jd(npt)
!    PRINT*,' MODEL START, END ',tstart,tend
    DO i=1,nwant
        jwant(i) = tstart + float(i-1)*dt
        iwant(i) = 1
    ENDDO
    DO i=1,npt
        magvar(i) = mag(i) - lpar(id(i))
    ENDDO

    PRINT*,'PREDICT CONTINUUM LIGHT CURVE'
    CALL PREDICT(sigma,tau,slag,swid,scale,salpha)
    ihires = 13
    WRITE(*,'(a32,1x,i2)')'WRITING TO FILE UNIT:',ihires

    DO i=1,nwant
        WRITE(ihires,'(3(1x,f13.6))')jwant(i),(mwant(i)+lpar(1))*renorm(1)+bar1(1),ewant(i)*renorm(1)
    ENDDO
    CLOSE(UNIT=ihires)

    DO i=2,ncurve
        DO j=1,nwant
            iwant(j) = i
        ENDDO
        PRINT*,'PREDICT LINE LIGHT CURVE', i
        CALL PREDICT(sigma,tau,slag,swid,scale,salpha)
        ihires = ihires - 1
        WRITE(*,'(a32,1x,i2)')'WRITING TO FILE UNIT:',ihires
        DO j=1,nwant
           WRITE(ihires,'(3(1x,f13.6))')jwant(j),(mwant(j)+lpar(i))*renorm(i)+bar1(i),ewant(j)*renorm(i)
        ENDDO
        CLOSE(UNIT=ihires)
    ENDDO
    DEALLOCATE(jwant,mwant,ewant,magvar,iwant)

ENDIF
END SUBROUTINE

!----------------------------------------------------------------------

SUBROUTINE PREDICT(sigma,tau,slag,swid,scale,salpha)
implicit none
REAL(kind=8),DIMENSION(:) ::  slag,swid,scale
REAL(kind=8) :: sigma,tau,salpha
REAL(kind=8) ::  cpnmatrix(npt,npt)
REAL(kind=8) ::  cplusninvdoty(npt),covar(npt),cplusninvdotcovar(npt)
REAL(kind=8) :: variance
INTEGER(kind=4) :: INFO,i,j,k
interface
    function getcmatrix(id1,id2,jd1,jd2,variance,tau,slag,swid,scale,salpha)
    REAL(kind=8) :: getcmatrix
    INTEGER(kind=4) ::  id1,id2
    REAL(kind=8) ::  jd1,jd2
    REAL(kind=8) ::  variance,tau,salpha
    REAL(kind=8),DIMENSION(:) ::  slag,swid,scale
    end function
end interface

variance = sigma*sqrt(0.5D0*tau)
! build the covariance plus noise matrix of the data
DO i=1,npt
  DO j=i,npt
    cpnmatrix(i,j)   = getcmatrix(id(i),id(j),jd(i),jd(j),&
    variance,tau,slag,swid,scale,salpha)
  ENDDO
  cpnmatrix(i,i) = cpnmatrix(i,i) + emag(i)*emag(i)
enddo
! fill in the matrix
DO i=1,npt
  DO j=1,i-1
    cpnmatrix(i,j) = cpnmatrix(j,i)
  ENDDO
ENDDO
! do the LU decomposition of the matrix
! now we want cpnmatrix^(-1)*mag = x, which is the same as
!    mag = cpnmatrix*x, so we solve this equation for x
!call Chol_GetSol(INFO,cpnmatrix,magvar,npt,cplusninvdoty)
call LU_GetSol(INFO,cpnmatrix,magvar,npt,cplusninvdoty)

! computes the model light curve C(t) (C+N)^(-1) mag
!                          error C(t) (C+N)^(-1) C(t)
! note -- it might be possible to do the error a faster way..
do i=1,nwant
    ! get the minimum variance estimate of flux at jwant(i) and its associated
    ! error
  do j=1,npt
    covar(j)   = getcmatrix(iwant(i),id(j),jwant(i),jd(j),&
    variance,tau,slag,swid,scale,salpha)
  enddo
!  call Chol_GetSol(INFO,cpnmatrix,covar(j),npt,cplusninvdotcovar)
  call LU_GetSol(INFO,cpnmatrix,covar,npt,cplusninvdotcovar)
  mwant(i) = 0.0D0
  ewant(i) = getcmatrix(iwant(i),iwant(i),jwant(i),jwant(i),&
  variance,tau,slag,swid,scale,salpha)
  do j=1,npt
    mwant(i) = mwant(i) + covar(j)*cplusninvdoty(j)
    ewant(i) = ewant(i) - covar(j)*cplusninvdotcovar(j)
  enddo
  ewant(i) = sqrt(ewant(i))
enddo
return
END SUBROUTINE

!----------------------------------------------------------------------
End Module
