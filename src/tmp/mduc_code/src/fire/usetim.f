c=======================================================================
      SUBROUTINE usetim(mph)
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
      USE comm0, ONLY : i_usetim
c-----
      IMPLICIT NONE
c-----
      INTEGER, INTENT(in) :: mph
c-----------------------------------------------------------------------
c-----
      IF(I_USETIM == 1) THEN  ! activation of usetim in the GUI (= 1)
c     mduc_usetim is not needed anymore
        CALL mduc_usetim(mph)
      END IF
c-----
      END SUBROUTINE usetim
