! Last-modified: 24 Apr 2011 06:01:16 PM

!-----------------------------------------------------------------------------------
! Solving Ax = b:
! for an input symmetric matrix A, LU decompose it to get the 
! product of A^-1 b
SUBROUTINE LU_GetSol(INFO,mat,vecb,nsize,vecsol)
implicit none
INTEGER(kind=4) :: nsize
REAL(kind=8) :: mat(nsize,nsize),vecb(nsize),vecsol(nsize)
REAL(kind=8) :: matlu(nsize,nsize),b(nsize,1)
INTEGER(kind=4) :: INFO,i
INTEGER(kind=4) :: IPIV(nsize)
! save the input matrix
matlu = mat
b(1:nsize,1) = vecb
CALL DGETRF(nsize, nsize, matlu, nsize, IPIV, INFO)
if(INFO .gt. 0) then
    print*,'WARNING: Matrix is singular'
    return
else if(INFO .lt. 0)then
    print*,'WARNING: Parameter input illegal!'
    stop
endif
CALL DGETRS( 'N', nsize, 1, matlu, nsize, IPIV, b, nsize, INFO )
if (INFO .ne. 0) then
    print*,'ERROR: Matrix cannot be solved!'
    stop
else
    vecsol = b(1:nsize,1) 
    return
endif
END  SUBROUTINE
!-----------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------
! for an input symmetric matrix A, Cholesky decompose it to get the 
! product of A^-1 b
SUBROUTINE Chol_GetSol(INFO,mat,vecb,nsize,vecsol)
implicit none
INTEGER(kind=4) :: nsize
REAL(kind=8) :: mat(nsize,nsize),vecb(nsize),vecsol(nsize)
REAL(kind=8) :: matchol(nsize,nsize),b(nsize,1)
INTEGER(kind=4) :: INFO,i
CHARACTER(len=1) :: UPLO
! save the input matrix
matchol = mat
b(1:nsize,1) = vecb
! cholesky decomposition of the matrix into a Lower triangular one
UPLO = 'L'
call DPOTRF(UPLO, nsize , matchol, nsize, INFO )
if(INFO .gt. 0) then
    print*,'WARNING: Matrix is not positive definite(Chol_GetSol)'
    return
else if(INFO .lt. 0)then
    print*,'WARNING: Parameter input illegal!'
    stop
endif
call DPOTRS(UPLO, nsize, 1, matchol, nsize, b, nsize, INFO )
if(INFO .ne. 0)then
    print*,'WARNING: DPOTRS failed!'
    stop
else
    vecsol = b(1:nsize,1) 
    return
endif
END SUBROUTINE
!-----------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------
! for an input symmetric matrix A, LU decompose it to get the 
! log(determinant) and the inversion.
subroutine LU_GetInvDet(INFO,mat,nsize,det,inv)
implicit none
integer(kind=4) :: nsize
real(kind=8) :: mat(nsize,nsize),inv(nsize,nsize)
!real(kind=8),DIMENSION(:,:) :: mat,inv
real(kind=8) :: det
integer(kind=4) :: INFO,i,flag
integer(kind=4) :: IPIV(nsize),LWORK
real(kind=8),allocatable :: WORK(:)
real(kind=8) :: testwork(1)

inv = mat

CALL DGETRF( nsize, nsize, inv, nsize, IPIV, INFO )
if(INFO .gt. 0) then
    print*,'WARNING: Matrix is singular'
    return
else if(INFO .lt. 0)then
    print*,'WARNING: Parameter input illegal!'
    stop
endif

det = 1.0D0
do i = 1, nsize
  det = det*inv(i,i)
  if ( IPIV(i) .ne. i ) then
      IPIV(IPIV(i)) = IPIV(i)
      det = -det
  endif
enddo

if(det .lt. 0.0D0)then
    print*,'ERROR: Negative Determinant in LU Decomposition!'
else
    det = LOG(det)
endif

LWORK = -1
CALL DGETRI( nsize, inv, nsize, IPIV, testwork, LWORK, INFO )
if (INFO .eq. 0) then
    LWORK = testwork(1)
    allocate(WORK(LWORK))
endif
CALL DGETRI( nsize, inv, nsize, IPIV, WORK, LWORK, INFO )
deallocate(WORK)
if (INFO .ne. 0) then
    print*,'ERROR: Matrix cannot be inverted!'
    stop
else
    return
endif
end  subroutine
!-----------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------
! for an input symmetric matrix A, Cholesky decompose it to get the 
! determinant and the inversion.
subroutine Chol_GetInvDet(INFO,mat,nsize,det,inv)
implicit none
! to enhance efficiency, we implicitly assume that the size of the 
! matrix is strictly nsize, so in the parent program, arrays should 
! be dynamically allocated.
integer(kind=4) :: nsize
real(kind=8) :: mat(nsize,nsize),inv(nsize,nsize)
!real(kind=8),DIMENSION(:,:) :: mat,inv
real(kind=8) :: det
character(len=1) :: UPLO
integer(kind=4) :: INFO

! save the input matrix
inv = mat
! cholesky decomposition of the matrix into a Lower triangular one
UPLO = 'L'
call DPOTRF(UPLO, nsize , inv, nsize, INFO )
if(INFO .gt. 0) then
    print*,'WARNING: Matrix is not positive definite(Chol_GetInvDet)'
    return
else if(INFO .lt. 0)then
    print*,'WARNING: Parameter input illegal!'
    stop
endif
! get determinant of inv(actually the real matrix)
call GetTriDet(inv,nsize,det)
! on exit, inv carries the lower triangle of the inverse of 'matrix'
call DPOTRI(UPLO, nsize, inv, nsize, INFO )
if(INFO .ne. 0)then
    print*,'WARNING: DPOTRI failed!'
    stop
endif
! expand the whole matrix so that the upper triangle is filled with a
! mirrored lower triangle, get the 'real' inv
call ExpMat(UPLO,inv,nsize)
return
end subroutine


SUBROUTINE PrtMat(mat,m,n,mwant,nwant)
implicit none
integer(kind=4) :: m,n,i,j,mwant,nwant
real(kind=8) :: mat(m,n)
print*,'printing matrix...'
do i=1,mwant
    do j=1,nwant
        WRITE(6,'((g14.5),$)')mat(i,j)
    enddo
    write(6,*)achar(14)
enddo
write(6,*)achar(10)
END SUBROUTINE

SUBROUTINE ExpMat(UPLO,mat,m)
implicit none
integer(kind=4) :: m,i,j
real(kind=8) :: mat(m,m)
character(len=1) :: UPLO

if(UPLO .eq. 'L') then
    do i=1,m
        do j=1,m
            if(i.lt.j)then
                mat(i,j) = mat(j,i)
            endif
        enddo
    enddo
elseif(UPLO .eq. 'U') then
    do i=1,m
        do j=1,m
            if(i.gt.j)then
                mat(i,j) = mat(j,i)
            endif
        enddo
    enddo
else
    print*,'UPLO is neither U or L, quit'
    stop
endif
END SUBROUTINE

SUBROUTINE GetTriDet(triangle,m,det)
implicit none
integer(kind=4) :: m,i
real(kind=8) :: triangle(m,m)
real(kind=8) :: det
! logarithmic determinant!
det = 0.0
do i = 1,m
! for pos-def matrices, L have strictly positive diagonal
! entries 
    det = det + LOG(triangle(i,i))
enddo
det = 2.0D0*det
END SUBROUTINE


