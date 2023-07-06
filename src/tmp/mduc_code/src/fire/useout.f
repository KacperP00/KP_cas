c=======================================================================
      SUBROUTINE useout(mph)
c=======================================================================
c
c     USEOUT IS A SPECIAL PURPOSE ROUTINE FOR GENERATING USER
c     DEFINED OUTPUT OR FOR CHANGING THE RUN PARAMETERS
c     ALL USER-SPECIFIED EXPRESSIONS MUST BE WRITTEN IN STANDARD
c     FORTRAN 90 OR IN THE FORTRAN VERSION AVAILABLE ON YOUR MACHINE
c     TO ACTIVATE THE STATEMENTS, PLEASE REMOVE THE 'C' FROM THE
c     FIRST COLUMN
c
c     THIS ROUTINE WILL BE EXECUTED AT THE END OF EVERY TIME-STEP
c
c     See the list of variables for the use of this routine
c
c.....contact cfd_support@avl.com
c-----
c=======================================================================
c-----
      USE comm0, ONLY : i_useout
c-----
      IMPLICIT NONE
c-----
      INTEGER, INTENT(in) :: mph
c-----
c-----------------------------------------------------------------------
c-----
      IF(I_USEOUT == 1) THEN     ! activation of useout in the GUI (= 1)
        CALL mduc_useout(mph)
      END IF
c-----
      END SUBROUTINE useout
