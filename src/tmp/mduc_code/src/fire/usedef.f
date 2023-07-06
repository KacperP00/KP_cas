c=======================================================================
      SUBROUTINE usedef(mat,mph)
c=======================================================================
c
c     USEDEF IS A SPECIAL PURPOSE ROUTINE TO MODIFY DEFINITION
c     OF SOME VARIABLES AND CONSTANTS
c     ALL USER-SPECIFIED EXPRESSIONS MUST BE WRITTEN IN STANDARD
c     FORTRAN 90 OR IN THE FORTRAN VERSION AVAILABLE ON YOUR MACHINE
c     TO ACTIVATE THE STATEMENTS, PLEASE REMOVE THE 'C' FROM THE
c     FIRST COLUMN
c
c     THIS ROUTINE WILL BE EXECUTED AT PROGRAM START/RESTART 
c
c     See the list of variables for the use of this routine
c
c.....contact cfd_support@avl.com
c-----
c-----------------------------------------------------------------------
c-----
      USE comm0, ONLY : i_useini,num_actscl,num_passcl
c-----
      IMPLICIT NONE
c-----
      INTEGER, INTENT(IN) :: mat,mph
c-----------------------------------------------------------------------
c-----
      IF (I_USEINI .GT. 0) THEN
        CALL mduc_usedef(mat,mph)
      END IF
c-----
      RETURN
      END SUBROUTINE usedef
