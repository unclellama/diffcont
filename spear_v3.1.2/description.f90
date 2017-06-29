!Last-modified: 24 Apr 2011 05:05:18 PM

!**********************************************************************

Module LightCurves
IMPLICIT NONE
! LIGHT CURVE DATA
!   jd = Julian Dates
!  mag = magnitude or flux [we use flux in our paper]
! emag = error in magnitude of flux
REAL(kind=8), allocatable :: jd(:),mag(:),emag(:)
!   id = id number of data point -- which light curve it belongs to
INTEGER(kind=4), allocatable ::  id(:)
REAL(kind=8), allocatable :: lmat(:,:)
End Module LightCurves

!**********************************************************************

Module TransferFunctions
IMPLICIT NONE
! TRANSFER FUNCTION TYPE
! options for what those modes are
!  DELTAFUNC    = expecting two curves but use a delta function (slow=shig)
!  TOPHAT       = tophat function
!  SINGLE       = just a single unsmoothed light curve
!  PHOTOECHO    = light curves from two broad bands, one with lines, the other without.
!  DOUBLETOPHAT = tophats for different emission lines or velocity bins
!  BERET        = luminosity dependent transfer function 
INTEGER(kind=4),parameter :: DELTAFUNC    = 0, &
                             TOPHAT       = 1, &
                             SINGLE       = 2, &
                             PHOTOECHO    = 3, &
                             DOUBLETOPHAT = 4, &
                             BERET        = 5
End Module TransferFunctions

!**********************************************************************

Module ModelParameters
IMPLICIT NONE
! VARIABLE CONTROL
! these are used to allow you to hold variables fixed or to
! minimize them: ivary(i) = 1 means to minimize, =0 means
! to not minimize and used the original, saved values
! psave(i)  
! ndim = number of total variables that could be minimized
REAL(kind=8), allocatable :: psave(:),lpar(:)
INTEGER(kind=4), allocatable :: ivary(:)
INTEGER(kind=4) :: nvar,iposdef
INTEGER(kind=4) :: convmode,ncurve,npt,ndim,nlin
! PRIOR ON CONTINUUM MODEL
! this sets up a "Gaussian" prior on the values of tau and sigma -- if
! you are fitting something other than a single light curve, the code
! presently assumes you have done an independent fit of the continuum
! light curve and generated its statistics -- then when fitting the
! multiple light curves it restricts tau/sigma by a prior set by
! their allowed range when fitting the continuum alone
REAL(kind=8) :: taucen,tauscatlow,tauscathig,sigcen,sigscatlow,sigscathig
! tmed   = the median time difference between adjacent light curve
!            points - sometimes used in the prior on tau
REAL(kind=8) :: tmed,tspan,chisave,plikesave,penalty
REAL(kind=8), allocatable :: bar1(:),renorm(:)
! cnorm = renorm(1)
REAL(kind=8) :: cnorm
End Module ModelParameters

!**********************************************************************

Module BeretMode
implicit none
!----------------------------------------------------------------------
! for BERET mode only
!       mmagc  = mean continuum flux
!  dmage(NMAX) = emission light curve subtracted by the emission flux mean
!    jde(NMAX) = time epochs of line light curve measurements
REAL(kind=8) :: mmagc
REAL(kind=8),allocatable :: dmage(:),jde(:)
! number of data points of the first emission line
INTEGER(kind=4) :: nem
End Module BeretMode

!**********************************************************************



