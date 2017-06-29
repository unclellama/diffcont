!   Last-modified: 24 Apr 2011 06:03:33 PM

SUBROUTINE Copy_ModelPara(conv_mode,ncurve,sig,tau,alpha,lag,wid,scale,p)
USE TransferFunctions
implicit none
REAL(kind=8) sig,tau,alpha
REAL(kind=8),DIMENSION(:) :: lag,wid,scale,p
INTEGER(kind=4),intent(IN) :: conv_mode,ncurve
INTEGER(kind=4) :: i
sig     = p(1)
tau     = p(2)
alpha   = 0.0D0
lag(1)  = 0.0D0
wid(1)  = 0.0D0
! scale(1) will be used for PHOTOECHO mode only
IF (conv_mode.ne.PHOTOECHO) scale(1)= 1.0D0

IF (conv_mode.eq.TOPHAT) THEN
  lag(2)  = p(3)
  wid(2)  = abs(p(5))
  scale(2)= abs(p(4))
ELSE IF (conv_mode.eq.BERET ) THEN
  lag(2)  = p(3)
  wid(2)  = abs(p(5))
  scale(2)= abs(p(4))
  alpha   = p(6)
ELSE IF (conv_mode.eq.DOUBLETOPHAT) THEN
  do i=2,ncurve
     lag(i)  = p((i-1)*3)
     wid(i)  = abs(p((i-1)*3+2))
     scale(i)= abs(p((i-1)*3+1))
  enddo
ELSE IF (conv_mode.eq.DELTAFUNC) THEN
   lag(2)  = p(3)
   wid(2)  = 0.0D0
   scale(2)= abs(p(4))
! this is a tophat centered at zero lag
ELSE IF (conv_mode.eq.PHOTOECHO) THEN
   lag(2)  = p(3)
   scale(1)= abs(p(4))
   scale(2)= abs(p(5))
   wid(2)  = 0.0D0
ENDIF
END SUBROUTINE Copy_ModelPara
