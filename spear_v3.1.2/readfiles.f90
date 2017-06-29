!Last-modified: 24 Apr 2011 05:20:01 PM


!**********************************************************************

Module ReadFiles
use TransferFunctions
use LightCurves
use ModelParameters
use BeretMode
implicit none
interface
    subroutine SORT3_LightCurve(n,ra,rb,rc,rd)
    INTEGER(kind=4) :: n
    REAL(kind=8) :: ra(n),rb(n),rc(n)
    INTEGER(kind=4) :: rd(n)
    end subroutine SORT3_LightCurve
    subroutine SORT(n,arr)
    INTEGER(kind=4) :: n
    REAL(kind=8) :: arr(n)
    end subroutine SORT
end interface
CONTAINS
!----------------------------------------------------------------------
subroutine Read_LC_File_A(filename)
implicit none
!----------------------------------------------------------------------
CHARACTER(len=*) :: filename
INTEGER(kind=4),parameter :: NMAX = 2000
INTEGER(kind=4) :: i,j,nline,ialloc
REAL(kind=8) :: tjd,tmag,temag
REAL(kind=8) :: datarray(NMAX,3)
INTEGER(kind=4) :: datid(NMAX)
OPEN(UNIT=10,FILE=filename,FORM='FORMATTED',STATUS='OLD')
! number of light curves in the input (including continuum)
READ(10,*)ncurve
! for those modes the light curve number is predetermined
IF (convmode.EQ.SINGLE) ncurve = 1
IF (convmode.EQ.TOPHAT) ncurve = 2
IF (convmode.EQ.BERET ) ncurve = 2
! let nlin equal to ncurve
nlin = ncurve
npt  = 0
DO i=1,ncurve
    READ(10,*)nline
    DO j=1,nline
        READ(10,*)tjd,tmag,temag
        IF(temag .GE. 1.D-3)THEN
        npt = npt + 1
        IF(npt .GE. NMAX) STOP 'ERROR: NUMBER OF LC POINTS EXCEEDS NMAX'
        datarray(npt,1)   = tjd
        datarray(npt,2)   = tmag
        datarray(npt,3)   = temag
        datid(npt)        = i
        ENDIF
    ENDDO
ENDDO
CLOSE(UNIT=10)
ALLOCATE(jd(npt),mag(npt),emag(npt),id(npt),lmat(npt,nlin),lpar(nlin),STAT=ialloc)
IF(ialloc .ne. 0) STOP 'ERROR: TRY TO ASSIGN A HUGE ARRAY'
jd   = datarray(1:npt,1)
mag  = datarray(1:npt,2)
emag = datarray(1:npt,3)
id   = datid(1:npt)
! sort the data points to be in temporal order (at least to ensure that)
CALL SORT3_LightCurve(npt,jd,mag,emag,id)
! set up Linear Parameter Matrix
do i=1,npt
  do j=1,nlin
    lmat(i,j) = 0.0D0
    IF (j.eq.id(i)) lmat(i,j) = 1.0D0
  enddo
enddo
end subroutine Read_LC_File_A
!----------------------------------------------------------------------
subroutine Read_Prior_File_A(filename)
implicit none
!----------------------------------------------------------------------
CHARACTER(len=*) :: filename
OPEN(UNIT=11,FILE=filename,FORM='FORMATTED',STATUS='OLD')
READ(11,*)taucen,tauscatlow,tauscathig,sigcen,sigscatlow,sigscathig
CLOSE(UNIT=11)
WRITE(*,'(a32,2(1x,f13.6))')'SIGMA AND TAU CENTERED AT:',10.D0**sigcen,10.D0**taucen
! WARNING: STATS.DAT HAS MEDIAN, 15% and 85% VALUES NOT SIGMAS
tauscatlow = taucen-tauscatlow
tauscathig = tauscathig-taucen
IF((tauscatlow.LE.0.0D0).OR.(tauscathig.LE.0.0D0))STOP 'ERROR: &
NEGATIVE WIDTH OF PRIOR FOUND'
end subroutine Read_Prior_File_A
!----------------------------------------------------------------------
subroutine Get_LC_STAT()
implicit none
!----------------------------------------------------------------------
REAL(kind=8),allocatable :: tlist(:),delts(:),bar0(:)
INTEGER(kind=4) :: nlist,ialloc,i,j,ilist,nmed
! work out the mean time spacings -- this is sometimes used
! in the prior on tau -- basically produces a list of the
! time differences between adjacent points, sorts that
! list and saves the median time spacing
nlist = count(id.eq.1)
ilist = 0
ALLOCATE(tlist(nlist),delts(nlist-1),STAT=ialloc)
IF(ialloc .ne. 0) STOP 'ERROR: TRY TO ASSIGN A HUGE ARRAY'
DO i=1,npt
  IF (id(i).EQ.1) THEN
    ilist = ilist + 1
    tlist(ilist) = jd(i)
  ENDIF
ENDDO
DO i=1,nlist-1
  delts(i) = tlist(i+1)-tlist(i)
  IF(delts(i).LT.0.0D0) PRINT*, 'WARNING: LIGHT CURVES NOT IN TEMPORAL ORDER'
ENDDO
CALL sort(nlist-1,delts)
nmed = int(float(nlist)/2.0D0+0.5D0)
! get the median of the continuum gap
tmed = delts(nmed)
! get the span of the continuum
tspan= tlist(nlist)-tlist(1)
! determine and subtract off the averages of each 
! light curve -- not strictly necessary except for BERET mode
ALLOCATE(bar0(ncurve),bar1(ncurve),renorm(ncurve))
bar0 = 0.0D0
bar1 = 0.0D0
DO i=1,npt
  j       = id(i)
  bar0(j) = bar0(j) +  1.0D0/emag(i)**2
  bar1(j) = bar1(j) + mag(i)/emag(i)**2
ENDDO
bar1 = bar1/bar0
DO i=1,ncurve
    renorm(i) = amplitude(bar1(i))
ENDDO
! record the renorm factor for continuum
cnorm = renorm(1)
! subtract the mean
DO i=1,npt
    mag(i) = (mag(i) - bar1(id(i)))/renorm(id(i))
   emag(i) =                emag(i)/renorm(id(i))
ENDDO
DEALLOCATE(bar0,tlist,delts)
end subroutine
!----------------------------------------------------------------------
subroutine Set_BeretMode()
implicit none
!----------------------------------------------------------------------
INTEGER(kind=4) :: i,j
IF (convmode .eq. BERET) THEN
    PRINT*,'SET UP FOR BERET MODE'
    mmagc = bar1(1)
! deviation of emission line flux from the mean
    nem =COUNT(id .eq. 2)
    ALLOCATE(dmage(nem),jde(nem))
    j = 0
    DO i=1,npt
      IF(id(i).EQ.2)THEN
          j        = j+1
          dmage(j) = mag(i)
          jde(j)   = jd(i)
      ENDIF
    ENDDO
! double check of nem
    IF(nem .NE. j)THEN
        PRINT*,'ERROR: NUMBER OF LINES DOESNT MATCH',nem,j
        STOP
    ENDIF
ELSE
    PRINT*,'NO SETUP FOR NON-BERET MODE'
ENDIF
end subroutine Set_BeretMode
!----------------------------------------------------------------------
function amplitude(xx)
implicit none
!----------------------------------------------------------------------
REAL(kind=8) :: xx,x,amplitude
! in case in magnitude system
x = abs(xx)
if ((x <= 100D0) .and. (x >= 0D0)) then
    ! IMPORTANT
    ! the numbers are generally safe in this regime,
    ! more importantly, numbers between 0 and 100 are
    ! possibly in magnitudes, which cannot be rescaled
    ! by constant factors.
    amplitude = 1.0D0
else
    ! we want to rescale light curves (in flux units) 
    ! that have (abs)numerical values outside of [0,100] to 
    ! avoid numerical instabilities in matrix operation.
    amplitude = 1.0D0
    if (x>=1.0D0) then
        do while(x > amplitude)
            amplitude=10.D0*amplitude
        enddo
        amplitude = amplitude/10.D0
    else
        do while(x < amplitude)
            amplitude=amplitude/10.D0
        enddo
    endif
endif
end function
!----------------------------------------------------------------------

End Module ReadFiles

!**********************************************************************
