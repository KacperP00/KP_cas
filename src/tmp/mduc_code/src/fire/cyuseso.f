c=======================================================================
      SUBROUTINE cyuseso(nd)
c=======================================================================
c
c-----------------------------------------------------------------------
c     cyuseso allows for user defined set-up of sources for additional 
c     transport equations. This routine is called inside spray integration
c     loop for each parcel within each subcycling time step. 
c
c     cyuseso is called by: cysourc
c     cyuseso calls       : div
c
c-----
c-----------------------------------------------------------------------
c-----
      USE cydrnoz
c-----
      IMPLICIT NONE
c-----
      INTEGER :: nd
c-----------------------------------------------------------------------
c---- compute evaporated mass from each nozzle for mean mixture fraction 
c---- source terms
      CALL mduc_cyuseso(nd)
c     
      END SUBROUTINE
c=======================================================================
