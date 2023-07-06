c=======================================================================
      subroutine use_cthchk
c=======================================================================
c
c     description: user function for resetting the reference
c                  temperature for the enthalpy equation
c
c     called from: cthdis
c
c-----
c-----------------------------------------------------------------------
c-----
      USE comm0, ONLY : i_useini
c-----
      IMPLICIT NONE
c-----------------------------------------------------------------------
c----- 
      IF (I_USEINI .GT. 0) THEN
        CALL mduc_cthchk
      END IF
c-----
      END SUBROUTINE use_cthchk
