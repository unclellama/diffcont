!Last-modified: 24 Apr 2011 04:51:49 PM

SUBROUTINE amoeba(p,y,ftol,func,iter)
IMPLICIT NONE
INTEGER(kind=4), INTENT(OUT) :: iter
REAL(kind=8), INTENT(IN) :: ftol
REAL(kind=8), DIMENSION(:), INTENT(INOUT) :: y
REAL(kind=8), DIMENSION(:,:), INTENT(INOUT) :: p
INTERFACE
   FUNCTION func(x)
   IMPLICIT NONE
   REAL(kind=8), DIMENSION(:), INTENT(IN) :: x
   REAL(kind=8) :: func
   END FUNCTION func
END INTERFACE
INTEGER(kind=4), PARAMETER :: ITMAX=5000
REAL(kind=8), PARAMETER :: TINY=1.0D-10
INTEGER(kind=4) :: ihi,ndim
REAL(kind=8), DIMENSION(size(p,2)) :: psum
call amoeba_private
CONTAINS

SUBROUTINE amoeba_private
IMPLICIT NONE
INTEGER(kind=4) :: i,ilo,inhi
REAL(kind=8) :: rtol,ysave,ytry,ytmp
IF((size(p,2) == size(p,1) - 1) .and. (size(p,2) == size(y) -1))THEN
    ndim = size(y) - 1
ELSE
    STOP 'ERROR: terminated in amoeba for inconsistent arr dimensions'
ENDIF
iter=0
psum(:)=sum(p(:,:),dim=1)
do
   ilo=iminloc(y(:))
   ihi=imaxloc(y(:))
   ytmp=y(ihi)
   y(ihi)=y(ilo)
   inhi=imaxloc(y(:))
   y(ihi)=ytmp
   rtol=2.0D0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
   if (rtol < ftol) then
      call swap_scalar(y(1),y(ilo))
      call swap_vector(p(1,:),p(ilo,:))
      RETURN
   end if
   if (iter >= ITMAX) STOP 'ERROR: ITMAX exceeded in amoeba'
   ytry=amotry(-1.0D0)
   iter=iter+1
   if (ytry <= y(ilo)) then
      ytry=amotry(2.0D0)
      iter=iter+1
   else if (ytry >= y(inhi)) then
      ysave=y(ihi)
      ytry=amotry(0.5D0)
      iter=iter+1
      if (ytry >= ysave) then
         p(:,:)=0.5D0*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
         do i=1,ndim+1
            if (i /= ilo) y(i)=func(p(i,:))
         end do
         iter=iter+ndim
         psum(:)=sum(p(:,:),dim=1)
      end if
   end if
end do
END SUBROUTINE amoeba_private

FUNCTION amotry(fac)
IMPLICIT NONE
REAL(kind=8), INTENT(IN) :: fac
REAL(kind=8) :: amotry
REAL(kind=8) :: fac1,fac2,ytry
REAL(kind=8), DIMENSION(size(p,2)) :: ptry
fac1=(1.0D0-fac)/ndim
fac2=fac1-fac
ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
ytry=func(ptry)
if (ytry < y(ihi)) then
   y(ihi)=ytry
   psum(:)=psum(:)-p(ihi,:)+ptry(:)
   p(ihi,:)=ptry(:)
end if
amotry=ytry
END FUNCTION amotry

FUNCTION imaxloc(arr)
REAL(kind=8), DIMENSION(:), INTENT(IN) :: arr
INTEGER(kind=4) :: imaxloc
INTEGER(kind=4), DIMENSION(1) :: imax
imax=maxloc(arr(:))
imaxloc=imax(1)
END FUNCTION imaxloc
FUNCTION iminloc(arr)
REAL(kind=8), DIMENSION(:), INTENT(IN) :: arr
INTEGER(kind=4) :: iminloc
INTEGER(kind=4), DIMENSION(1) :: imax
imax=minloc(arr(:))
iminloc=imax(1)
END FUNCTION iminloc
SUBROUTINE swap_scalar(a,b)
REAL(kind=8), INTENT(INOUT) :: a,b
REAL(kind=8) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_scalar
SUBROUTINE swap_vector(a,b)
REAL(kind=8), DIMENSION(:), INTENT(INOUT) :: a,b
REAL(kind=8), DIMENSION(SIZE(a)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_vector
END SUBROUTINE amoeba
