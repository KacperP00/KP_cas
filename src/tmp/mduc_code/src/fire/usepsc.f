c=======================================================================
      SUBROUTINE usepsc(iconv,mph,kpsc)
c=======================================================================
c
c     USEPSC IS A SPECIAL PURPOSE ROUTINE TO MODIFY SOURCE TERM IN
c     SCALAR EQUATION
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
      USE comm0, ONLY : num_actscl,num_passcl,itst,nsp,nep
      USE comm1, ONLY : te,ed
c-----
      IMPLICIT NONE
c-----
      INTEGER, INTENT(IN) :: iconv,mph,kpsc
      INTEGER :: nc,mat
      mat = 1
c-----------------------------------------------------------------------
c-----
      IF (num_actscl > 0) THEN
        CALL mduc_usepsc(iconv,mph,kpsc)
      END IF
c-----
      IF (num_passcl > 0) THEN
      END IF
c-----
      RETURN
c-----
      END SUBROUTINE usepsc
