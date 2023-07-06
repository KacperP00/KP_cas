c=======================================================================
      SUBROUTINE use_cthtmp
c=======================================================================
c
c     description: user function for calculating the temperature 
c                  of gas mixture
c
c     called from: cthdis
c
c-----------------------------------------------------------------------
c-----
      USE comm0, ONLY : i_useini
c-----
      IMPLICIT NONE
c-----------------------------------------------------------------------
c----- 
      IF (I_USEINI .GT. 0) THEN
        CALL mduc_cthtmp
      ELSE
        CALL cthtmp
      END IF
c-----      
      END SUBROUTINE use_cthtmp
c-----
