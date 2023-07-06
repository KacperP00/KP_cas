c=======================================================================
      SUBROUTINE mduc_useent(mat,mph)
c=======================================================================
c
c     USEEENT IS A SPECIAL PURPOSE ROUTINE TO MODIFY
c     SOURCE TERMS IN ETHALPY EQUATION
c     ALL USER-SPECIFIED EXPRESSIONS MUST BE WRITTEN IN STANDARD
c     FORTRAN 90 OR IN THE FORTRAN VERSION AVAILABLE ON YOUR MACHINE
c     TO ACTIVATE THE STATEMENTS, PLEASE REMOVE THE 'C' FROM THE
c     FIRST COLUMN
c
c     THIS ROUTINE WILL BE EXECUTED AFTER EACH ITERATION
c
c     See the list of variables for the use of this routine
c
c.....contact cfd_support@avl.com
c-----
c-----------------------------------------------------------------------
c-----
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : nsp,nep,dtin,itst,i_useent,ispray,iampro
      USE comm1,    ONLY : su1
      USE cydrnoz,  ONLY : ssmass,ieulag
      USE mduc,     ONLY : H_vap_ref,iHCHEM
      IMPLICIT NONE
c-----
      INTEGER, INTENT(in) :: mat,mph
      INTEGER  :: nc
c-----
c-----------------------------------------------------------------------
c-----
c      IF (I_USEENT == 1 .AND. itst>1) THEN
      IF (I_USEENT == 1 .AND. ISPRAY == 1) THEN
c-----
        DO nc=nsp(mat),nep(mat)
c     add the formation enthalpies to the total enthalpy source
          su1(nc) = su1(nc) + H_vap_ref*ssmass(nc)*dtin*ieulag(4)
c-----
        END DO
      END IF
c-----
      RETURN
c-----
      END SUBROUTINE mduc_useent
