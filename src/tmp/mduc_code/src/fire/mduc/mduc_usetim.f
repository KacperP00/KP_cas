c=======================================================================
      SUBROUTINE mduc_usetim(mph)
c=======================================================================
c
c     USETIM IS A SPECIAL PURPOSE ROUTINE FOR GENERATING USER
c     DEFINED OUTPUT OR FOR CHANGING THE RUN PARAMETERS
c     ALL USER-SPECIFIED EXPRESSIONS MUST BE WRITTEN IN STANDARD
c     FORTRAN 90 OR IN THE FORTRAN VERSION AVAILABLE ON YOUR MACHINE
c     TO ACTIVATE THE STATEMENTS, PLEASE REMOVE THE 'C' FROM THE
c     FIRST COLUMN
c
c     THIS ROUTINE WILL BE EXECUTED AT THE START OF EVERY TIME-STEP
c
c     See the list of variables for the use of this routine
c
c.....contact cfd_support@avl.com
c-----
c=======================================================================
c-----
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : i_usetim,dt,dto,itst,iter,iters,time,timo,
     x                     nsp,nep,pref,iampro
      USE comm1,    ONLY : p,vol,t,te,ed,pasphi
      USE mduc
c-----
      IMPLICIT NONE
c-----
      INTEGER, INTENT(in) :: mph
c-----
      INTEGER    :: mat,isc,nc
      REAL(prec) :: totVol
c-----
c-----------------------------------------------------------------------
c-----
      mat= 1 
c      write(*,*) iampro, ': in mduc_usetim'
      IF (iter == iters) THEN
c       make sure that time step number is updated for restart
        isc = itst-1
        CALL MDUCupdateTime(timo,dt,isc,mduc_mem)
c-----  set Pavg_old to the right value for restart
        Pavg_old  = 0.0
        totVol    = 0.0
        DO nc=nsp(mat),nep(mat)
          Pavg_old = Pavg_old + (pref(mat)+p(nc))*vol(nc)
          totVol   = totVol + vol(nc)
        END DO
c     reduce across processors
        CALL dglsum(Pavg_old)
        CALL dglsum(totVol)
c     finish normalization
        Pavg_old = Pavg_old/totVol
        CALL MDUCsetPressure(-1,Pavg_old,mduc_mem)
        WRITE (*,*) 'Pavg_old:', Pavg_old
c-----
        RETURN                  ! don't do this before first time step
      END IF
      mat= 1                          ! number of materials
c-----
      IF(I_USETIM == 1) THEN          ! activation of usetim in the GUI (= 1)
c----- add code here for every timestep
        DO nc=nsp(mat),nep(mat)
c----- reset mean mixture fraction sources for next spray loop
          pasphi(nc+iMEVAP1) = 0.0
          IF (mduc_nInj .EQ. 2) THEN
            pasphi(nc+iMEVAP2) = 0.0
          END IF
        END DO
      END IF
c-----
c-----
      END SUBROUTINE mduc_usetim
