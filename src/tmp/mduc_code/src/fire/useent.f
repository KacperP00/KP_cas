c=======================================================================
      SUBROUTINE useent(mat,mph)
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
      USE comm0, ONLY : i_useent,nsp,nep
      USE comm1, ONLY : te,ed
c-----
      IMPLICIT NONE
c-----
      INTEGER, INTENT(IN) :: mat,mph
c-----
      INTEGER :: nc
c-----------------------------------------------------------------------
c-----
      IF (I_USEENT == 1) THEN
        CALL mduc_useent(mat,mph)
      END IF
c-----
      RETURN
      END SUBROUTINE useent
