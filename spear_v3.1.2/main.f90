!   Last-modified: 09 Feb 2011 04:13:47 PM

!**********************************************************************
PROGRAM MAIN
use LightCurves
use TransferFunctions
use ModelParameters
use BeretMode
use ReadFiles
use run_amoeba
use HighRes
use MCMC
IMPLICIT NONE
!**********************************************************************
CHARACTER(len=100) :: fname,statfile,mcmcfile
! variable format expression
CHARACTER(len=2) :: vfe_str,vfe_lin,vfe_cur
INTEGER(kind=4) :: i,j,k
INTEGER(kind=4) :: ialloc,iamoeba,icurve,idomcmc,ndof
LOGICAL :: ok
REAL(kind=8) :: chibest,pbest
REAL(kind=8),ALLOCATABLE :: bscale(:),bswid(:),bslag(:)
REAL(kind=8) :: sbest,tbest,bsalpha,plike1st
! for benchmarking
INTEGER(kind=4) :: count_0, count_1, count_rate, count_max
!----------------------------------------------------------------------
INTERFACE
    SUBROUTINE Copy_ModelPara(conv_mode,n_curve,var,tau,alpha,lag,wid,scale,p)
    REAL(kind=8) var,tau,alpha
    REAL(kind=8),DIMENSION(:) :: lag,wid,scale,p
    INTEGER(kind=4) :: conv_mode,n_curve
    END SUBROUTINE
END INTERFACE
!----------------------------------------------------------------------
!PRINT*,'REMINDER: if the program ever gets SEGMENTATION FAULT because of'
!PRINT*,'REMINDER: running our of memory for large arrays, please'
!PRINT*,'REMINDER: set the stack size as unlimited: ulimit -s unlimited'
PRINT*,'PLEASE ENTER THE LIGHTCURVE FIlE:'
READ(*,'(A)')fname
INQUIRE(FILE=fname, exist=ok)
IF( .NOT. ok )THEN
    WRITE(*,'(a32,1x,a40)')'ERROR: CANNOT FIND:',fname
    STOP
ELSE
    WRITE(*,'(a32,1x,a40)')'LIGHT CURVE DATA FILE:',fname
ENDIF
!----------------------------------------------------------------------
PRINT*,'*****************************************************'
PRINT*,'PLEASE ENTER THE CONVOLUTION MODE:'
PRINT*,'  DELTAFUNC    = 0                PHOTOECHO    = 3  *'
PRINT*,'  TOPHAT       = 1                DOUBLETOPHAT = 4  *'
PRINT*,'  SINGLE       = 2                BERET        = 5  *'
READ*,convmode
SELECT CASE (convmode)
 CASE(0);      PRINT*,'                    USING MODE: DELTAFUNC'
 CASE(1);      PRINT*,'                    USING MODE: TOPHAT'
 CASE(2);      PRINT*,'                    USING MODE: SINGLE'
 CASE(3);      PRINT*,'                    USING MODE: PHOTOECHO(TEST)'
 CASE(4);      PRINT*,'                    USING MODE: DOUBLETOPHAT'
 CASE(5);      PRINT*,'                    USING MODE: BERET'
 CASE DEFAULT; STOP   '                    ERROR: NO THIS MODE EXISTS'
END SELECT
PRINT*,'*****************************************************'
!----------------------------------------------------------------------
CALL Read_LC_File_A(fname)
! set the number of parameters given the convolution mode
IF (convmode.EQ.DELTAFUNC)    ndim = 4
IF (convmode.EQ.TOPHAT)       ndim = 5
IF (convmode.EQ.BERET )       ndim = 6
IF (convmode.EQ.SINGLE)       ndim = 2
IF (convmode.EQ.PHOTOECHO)    ndim = 5
! ndim of DOUBLETOPHAT mode depends on the # of light curves
IF (convmode.EQ.DOUBLETOPHAT) ndim = 2 + 3*(ncurve-1)
!----------------------------------------------------------------------
WRITE(vfe_str, '(i2)' )ndim
WRITE(vfe_lin, '(i2)' )nlin
WRITE(vfe_cur, '(i2)' )ncurve
!----------------------------------------------------------------------
ALLOCATE(psave(ndim),ivary(ndim),STAT=ialloc)
IF(ialloc .ne. 0) STOP 'ERROR: NO ENOUGH MEMORY'
PRINT*,'PLEASE ENTER THE INITIAL PARAMETERS:'
READ*,(psave(i),i=1,ndim) 
WRITE(*,'('//vfe_str//'(1x,f13.6))')(psave(i),i=1,ndim)
PRINT*,'PLEASE SET THE PARAMETERS TO BE FIXED (FIXED:0 ; VARY: 1):'
READ*,(ivary(i),i=1,ndim)
WRITE(*,'('//vfe_str//'(1x,i13))')(ivary(i),i=1,ndim)
! set the number of parameters that will be varied during minimization
!----------------------------------------------------------------------
nvar = count(ivary .eq. 1)
PRINT*,'*****************************************************'
WRITE(*,'(a32,1x,i4)')'NUMBER OF LIGHT CURVES:',ncurve
WRITE(*,'(a32,1x,i4)')'NUMBER OF LINEAR PARAMETERS:',nlin
WRITE(*,'(a32,1x,i4)')'NUMBER OF GOOD POINTS OBTAINED:',npt
WRITE(*,'(a32,1x,i4)')'NUMBER OF MODEL PARAMETERS:',ndim
WRITE(*,'(a32,1x,i4)')'NUMBER OF VARYING PARAMETERS:',nvar
PRINT*,'*****************************************************'
!----------------------------------------------------------------------
PRINT*,'FIND THE LOCAL OPTIMAL PARAMETERS?(1=YES)'
READ*,iamoeba
IF(iamoeba .eq. 1)THEN 
    WRITE(*,'(a32)')'YES'
ELSE
    WRITE(*,'(a32)')'NO'
ENDIF
!----------------------------------------------------------------------
PRINT*,'GENERATE A LIGHT CURVE?(1=YES)'
READ*,icurve
IF(icurve .eq. 1)THEN 
    WRITE(*,'(a32)')'YES'
ELSE
    WRITE(*,'(a32)')'NO'
ENDIF
!----------------------------------------------------------------------
PRINT*,'DO MCMC CHAIN?(1=YES)'
READ*,idomcmc
IF(idomcmc .eq. 1)THEN 
    WRITE(*,'(a32)')'YES'
ELSE
    WRITE(*,'(a32)')'NO'
ENDIF
IF (idomcmc.eq.1) THEN
    PRINT*,'PLEASE ENTER THE NAME FOR MCMC OUTPUT:'
    READ(*,'(a)')mcmcfile
    WRITE(*,'(a32,1x,a40)')'LIGHT CURVE DATA FILE:',mcmcfile
ENDIF
!----------------------------------------------------------------------
IF (convmode.NE.SINGLE) THEN
    statfile = 'stats.dat'
    INQUIRE(FILE=statfile, EXIST=ok)
    IF( .NOT. ok )THEN
        statfile = '../single/stats.dat'
        INQUIRE(FILE=statfile, EXIST=ok)
        IF( .NOT. ok )THEN
             STOP 'ERROR: CAN NOT FIND stats.dat'
        ENDIF
    ENDIF
    PRINT*,'*****************************************************'
    CALL Read_PRIOR_File_A(statfile)
ENDIF
!----------------------------------------------------------------------
PRINT*,'*****************************************************'
CALL Get_LC_STAT()
WRITE(*,'(a32,1x,f13.6)')'MEDIAN TIME SAMPLING IS:',tmed
DO i=1,ncurve
  WRITE(*,'(a32,1x,i2,1x,a4,f13.6)')'MEAN FLUX/MAG FOR LIGHT CURVE:',i,'IS:',bar1(i)
  WRITE(*,'(a32,1x,i2,1x,a4,f13.6)')'NORMALIZED AMP FOR LIGHT CURVE:',i,'IS:',renorm(i)
ENDDO
PRINT*,'*****************************************************'
!----------------------------------------------------------------------
! for BERET mode
CALL Set_BeretMode()
PRINT*,'*****************************************************'
!----------------------------------------------------------------------
! Calculate the likelihood for the initial parameters (for benchmarking)
CALl SYSTEM_CLOCK(count_0, count_rate, count_max)
CALL init_plike(plike1st)
WRITE(*,'(a32,1x,f13.6)')'LIKELIHOOD OF THE INPUT MODEL:',plike1st
PRINT*,'*****************************************************'
CALL SYSTEM_CLOCK(count_1, count_rate, count_max)
WRITE(*,'(a32,1x,f13.6)')'TIME ELAPSED FOR ABOVE RUN:',float(count_1 - count_0)*1.0D0/float(count_rate);
! if mcmc without amoeba
plikesave = plike1st
!----------------------------------------------------------------------
IF(iamoeba .eq. 1)THEN
    WRITE(*,'(a32,1x,i2,1x,a14,i2)')'ABOUT TO RUN CONVMODE:',convmode,'WITH NPARAMS:',nvar
    CALL amoeba_minimize(pbest)
    chibest   = chisave
    plikesave = pbest
    WRITE(*,'(a32)')'LINEAR PARAMETERS:'
    WRITE(*,'('//vfe_lin//'(1x,f13.6))')(lpar(i),i=1,nlin)
    PRINT*,'******************************************************'
!----------------------------------------------------------------------
! write out and save the this best new model
! output to the stdout
    WRITE(*,'(a32,1x,f13.6)')'LIKELIHOOD OF THE BEST MODEL:',pbest
    WRITE(*,'(a32)')'PARAMETERS OF THE BEST MODEL:'
    WRITE(*,'('//vfe_str//'(1x,f13.6))')(psave(i),i=1,ndim)
    PRINT*,'******************************************************'
! output of best fit model to file 
    ndof = npt - nvar - ncurve
    OPEN(UNIT=12,FILE='best3.new',FORM='FORMATTED',STATUS='UNKNOWN')
    WRITE(12,'(i13,1x,a22)')convmode,'convmode'
    WRITE(12,'(i13,1x,a22)')ndim,'npara'
    WRITE(12,'(f13.6,1x,a22)')pbest,'maxlikelihood'
    WRITE(12,'(f13.6,1x,a22)')chibest,'chi'
    WRITE(12,'(i13,1x,a22)')ndof,'ndof'
    WRITE(12,'('//vfe_cur//'(f13.6),1x,a22)')(renorm(j),j=1,ncurve),'rescaled'
    WRITE(12,'('//vfe_cur//'(f13.6),1x,a22)')(bar1(j),j=1,ncurve),'mean'
    WRITE(12,'('//vfe_lin//'(f13.6),1x,a22)')(lpar(i),i=1,nlin),'linearshifts'
    WRITE(12,'('//vfe_str//'(f13.6),1x,a22)')(psave(j),j=1,ndim),'bestfit'
    CLOSE(unit=12)
ENDIF
!----------------------------------------------------------------------
! copy these out into variables that make more
! "physical" sense
! best estimated sigma and tau
ALLOCATE(bscale(ncurve),bswid(ncurve),bslag(ncurve))
CALL Copy_ModelPara(convmode,ncurve,sbest,tbest,bsalpha,bslag,bswid,bscale,psave)
!----------------------------------------------------------------------
CALL Get_HighResLF(icurve,sbest,tbest,bslag,bswid,bscale,bsalpha)
!----------------------------------------------------------------------
! assume you will always use the psave from amoeba
CALL Run_MCMC(idomcmc,mcmcfile)
!----------------------------------------------------------------------
DEALLOCATE(bscale,bswid,bslag)
DEALLOCATE(psave,ivary)
DEALLOCATE(jd,mag,emag,id,lmat,lpar,bar1,renorm)
IF (convmode .eq. BERET)DEALLOCATE(dmage,jde)
END PROGRAM MAIN

!**********************************************************************
