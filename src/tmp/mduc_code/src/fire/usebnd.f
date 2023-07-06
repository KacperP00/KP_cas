c=======================================================================
      SUBROUTINE usebnd(k,mph)
c=======================================================================
c
c     USEBND IS A SPECIAL PURPOSE ROUTINE FOR CHANGING THE
c     VALUES AT THE BOUNDARIES
c     ALL USER-SPECIFIED EXPRESSIONS MUST BE WRITTEN IN STANDARD
c     FORTRAN 90 OR IN THE FORTRAN VERSION AVAILABLE ON YOUR MACHINE
c     TO ACTIVATE THE STATEMENTS, PLEASE REMOVE THE 'C' FROM THE
c     FIRST COLUMN
c
c     THIS ROUTINE WILL BE EXECUTED AT THE START OF CALCULATION
c
c     See the list of variables for the use of this routine
c
c.....parameter k       value which shows the boundary region,
c                       see .log file
c
c.....contact cfd_support@avl.com
c-----
c--------------------------------------------------------------------
c-----
      USE comm0, ONLY : i_usebnd
c-----
      IMPLICIT NONE
c-----
      INTEGER, INTENT(in) :: k,mph
c-----
c--------------------------------------------------------------------
c----- activation of usebnd in the GUI
      IF(I_USEBND .GT. 0) THEN
        CALL mduc_usebnd(k,mph)
      END IF
c-----
      END SUBROUTINE usebnd
c--------------------------------------------------------------------
