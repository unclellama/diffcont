! Last-modified: 24 Apr 2011 05:14:14 PM

!***********************************************************************
! Description 
!***********************************************************************
! JD2LAG implements the luminosity dependence of lag for the BERET mode
!***********************************************************************

!***********************************************************************
! Changelog 
!***********************************************************************
! 07 Oct 2009:
! Adapted from NR routine (locate.f)
! 07 Oct 2009:
! We know that our jd dates are increasing (jde)
! 15 Oct 2009:
! Fixed "Negative Number to Fractional Power" problem
!***********************************************************************

!***********************************************************************
! SUB_JD2LAG
! given a julian date x, the routine will find out the correspondant 
! luminosity dependent lag
!***********************************************************************
SUBROUTINE jd2lag(lag,tjd,scale,alpha,lag0)
use BeretMode
! output
REAL(kind=8) ::  lag
! inputs (scale(2) and lag0(2))
REAL(kind=8) ::  tjd,scale,alpha,lag0
! the index of tjd in jde array
INTEGER(kind=4) ::  jid
! flux deviation at jid
REAL(kind=8) ::  magejid
! fractional variation of flux
REAL(kind=8) ::  df

INTEGER(kind=4) ::  jl,jm,ju
! bisection search
jl=0
ju=nem+1
DO WHILE(ju-jl.gt.1)
    jm=(ju+jl)/2
    if(tjd.gt.jde(jm))then
      jl=jm
    else
      ju=jm
    endif
ENDDO
jid=jl+1
magejid = dmage(jid)
scale = abs(scale)
! calculate lag (ignore lpar here)
! exact formula should be,
! lag = lag0*((((magejid-lpar2)/scale)+mmagc+lpar1)/mmagc)**alpha)
! but we will use
! lag = lag0*((magejid/(scale*mmagc)+1.)**alpha)
df = magejid/(scale*mmagc)+1.D0
if(df.lt.0.0D0)df=0.D0
lag = lag0*(df**alpha)
RETURN
END
!************************************************************************
