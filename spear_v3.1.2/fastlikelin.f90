! Last-modified: 24 Apr 2011 06:08:05 PM

! Description
! fastlikelin computes the likelihood with linear parameters the fast 
! way, using cholesky decompositions

! Changelog
! 03 Nov 2010:
! Change to F90
! 15 Jul 2009:
! Add DOUBLETOPHAT mode (YZ)
! 25 Sep 2009:
! Change from slow,shig to slag,swid (YZ)
! 07 Oct 2009:
! Add BERET mode

SUBROUTINE fastlikelin(sigma,tau,slag,swid,scale,salpha,lcovar,plike,chi2)
use LightCurves
use ModelParameters
use BeretMode
implicit none
REAL(kind=8) :: cmatrix(npt,npt),cpnmatrix(npt,npt)
REAL(kind=8) :: cpninvmatrix(npt,npt),det_cpn,det_tlmat2
REAL(kind=8) :: linpar(nlin)
REAL(kind=8) :: tlmat1(nlin,npt),tlmat3(nlin,npt)
REAL(kind=8) :: tlmat2(nlin,nlin),tlmat2inv(nlin,nlin)
REAL(kind=8) :: sigma,tau,variance,plike,chi2
REAL(kind=8),DIMENSION(:) :: scale,slag,swid
REAL(kind=8),DIMENSION(:,:) :: lcovar
REAL(kind=8) :: salpha
INTEGER(kind=4) :: INFO,i,j,k
interface
    function getcmatrix(id1,id2,jd1,jd2,variance,tau,slag,swid,scale,salpha)
    REAL(kind=8) :: getcmatrix
    INTEGER(kind=4) ::  id1,id2
    REAL(kind=8) ::  jd1,jd2
    REAL(kind=8) ::  variance,tau,salpha
    REAL(kind=8),DIMENSION(:) ::  slag,swid,scale
    end function
    subroutine LU_GetInvDet(INFO,mat,nsize,det,inv)
    integer(kind=4) :: nsize,INFO
    real(kind=8) :: mat(nsize:nsize),inv(nsize:nsize)
    real(kind=8) :: det
    end subroutine
    subroutine Chol_GetInvDet(INFO,mat,nsize,det,inv)
    integer(kind=4) :: nsize,INFO
    real(kind=8) :: mat(nsize:nsize),inv(nsize:nsize)
    real(kind=8) :: det
    end subroutine
    subroutine PrtMat(mat,m,n,mwant,nwant)
    integer(kind=4) :: m,n,i,j,mwant,nwant
    real(kind=8) :: mat(m,n)
    end subroutine
end interface

! NOTE: the variance here actually is the standard deviation of the light curve, 
!       the REAL variance is the square of this quantity.
variance = sigma*sqrt(0.5D0*tau)

! build upper half the covariance plus noise matrix
do i=1,npt
  do j=i,npt
! cmatrix: S
    cmatrix(i,j)   = getcmatrix(id(i),id(j),jd(i),jd(j),variance,tau,slag,swid,scale,salpha)
    cpnmatrix(i,j) = cmatrix(i,j)
  enddo
! cpn    : S+N
  cpnmatrix(i,i) = cpnmatrix(i,i) + emag(i)*emag(i)
enddo

! fill in the rest of the matrix
do i=1,npt
  do j=1,i-1
    cmatrix(i,j)   = cmatrix(j,i)
    cpnmatrix(i,j) = cpnmatrix(j,i)
  enddo
enddo

! -------------------------------------------------------
! if you want cholesky decompostion (more efficient but less stable)
call Chol_GetInvDet(iposdef,cpnmatrix,npt,det_cpn,cpninvmatrix)
! if you want LU decompostion (less efficient but more stable)
!call LU_GetInvDet(iposdef,cpnmatrix,npt,det_cpn,cpninvmatrix)
if (iposdef.ne.0) then
  print*,'WARNING: NON POSITIVE DEFINITE MATRIX(fastlikelin)'
  do j=1,nlin
    print*,'line',j,'sca,lag,wid,alpha',scale(j),slag(j),swid(j),salpha
  enddo
  return
endif

! ---------------------------------------------------------
! work out the optimal linear parameters
! lpar: q 
! lmat: L
!  mag: y
!      q = (L^T invC L)^(-1) L^T invC y
!   lpar = (lmat^T * cpninvmatrix * lmat)^(-1) lmat^T cpninvmatrix mag
!               tlmat1 = L^T invC
! lets work out tlmat1 = lmat^T cpninvmatrix 
CALL DGEMM ( 'T', 'N', nlin, npt, npt,1.0D0, lmat, npt, &
cpninvmatrix, npt, 0.0D0, tlmat1, nlin )
!                 q =          L^T invC y     = tlmat1 y
! now work out lpar = lmat^T cpninvmatrix mag = tlmat1 mag
CALL DGEMV('N',nlin,npt,1.0D0,tlmat1,nlin,mag,1,0.0D0,linpar,1)
!             tlmat2 =        L^T invC L        = tlmat1 L
! now workout tlmat2 = lmag^T cpninvmatrix lmat = tlmat1 lmat
CALL DGEMM ( 'N', 'N',nlin,nlin,npt,1.0D0, tlmat1,nlin, &
lmat,npt, 0.0D0,tlmat2,nlin)
! tlmat2 is of small size (nlin x nlin)
! do the LU decomposition of the matrix
!CALL LU_GetInvDet(INFO,tlmat2,nlin,det_tlmat2,tlmat2inv)
! do the Chol decomposition of the matrix
CALL Chol_GetInvDet(INFO,tlmat2,nlin,det_tlmat2,tlmat2inv)
! normally we should do LU substitute to get the linear parameters lpar
! solve the equation to get lpar/q:
!              q = (L^T invC L)^(-1) L^T invC y  
! but here the size is small and we've got the inverse (and we have to!), so
CALL DGEMV('N',nlin,nlin,1.0D0,tlmat2inv,nlin,linpar,1,0.0D0,lpar,1)

! now if we just wanted to solve for the fit, we could use 
!       chi^2 = y^T (S+N)^(-1) (y - Lq)
!       chi^2 = mag^T cpninvmatrix^(-1) (mag-lmat*lpar)

! for Lcovar, it dates back to the derivation in the Appendix of RH
! but we also need the determinant, which means updating the matrix
!       (S+N)^(-1) - L^T  (S+N)^(-1) (L^T invC L)^(-1)  L^T (S+N)^(-1)
!       cpninvmatrix -> cpninvmatrix - tlmat1 lcovar tlmat1 
! so we actually need to invert tlmat2 (which will also be the 
! covariance matrix of the linear parameters)
! note, however, that this could be done faster if we do not want 
! the covariance matrix by simultaneously
! doing lcovar  tlmat1
!           tlmat2 =        L^T invC L        = tlmat1 L
!           lcovar =       (L^T invC L)^(-1)     
do i=1,nlin
    do j=1,nlin
        lcovar(i,j) = tlmat2inv(i,j)
    enddo
enddo
! now compute lcovar *tlmat1
CALL DGEMM('N','N',nlin,npt,nlin,1.0D0,tlmat2inv,nlin,&
tlmat1,nlin,0.0D0,tlmat3,nlin)
! now use cmatrix for storage
CALL DGEMM('T','N',npt,npt,nlin,1.0D0,tlmat1,nlin,&
tlmat3,nlin,0.0D0,cmatrix,npt)
do i=1,npt
  do j=1,npt
      cpninvmatrix(i,j) = cpninvmatrix(i,j) - cmatrix(i,j)
  enddo
enddo
! determine the chi^2
!        Chi = y^T[invC - invC L (L^T invC L)^(-1) L^T invC] y
chi2  = 0.0D0
do i=1,npt
  do j=1,npt
    chi2   = chi2   + mag(i)*mag(j)*cpninvmatrix(i,j)
  enddo
enddo
! NOTE!!! the determinants were already taken log_e in the subroutines
plike = -0.5D0*chi2-0.5D0*(det_tlmat2+det_cpn)
return
END SUBROUTINE
