C=======================================================================
      subroutine mduc_cyuseso(nd)
C=======================================================================
c
C-----------------------------------------------------------------------
C     cyuseso allows for user defined set-up of sources for additional 
C     transport equations. This routine is called inside spray integration
C     loop for each parcel within each subcycling time step. 
C
C     cyuseso is called by: cysourc
C     cyuseso calls       : div
C
c=======================================================================
c=======================================================================
c=======================================================================
c
      USE prec_mod, ONLY : prec
      USE cydrnoz,  ONLY : cdropo,ddropo,sdropo,
     x                     cdropn,ddropn,sdropn,
     x                     ncelln,nozzln
      USE mduc,     ONLY : iMEVAP1, iMEVAP2
      USE comm1,    ONLY : pasphi
c-----
      IMPLICIT NONE
c-----
c
c---- input
      INTEGER :: nd
c---- local
      INTEGER :: nc
      REAL(prec) :: dmo, dmn, ddm
      REAL(prec) :: M_PI = 3.141592653589793
c
c---- compute evaporated mass from each nozzle for mean mixture fraction
c---- source terms
      dmo = cdropo*M_PI*ddropo*sdropo**3/6.d0             ! old mass
      dmn = cdropn(nd)*M_PI*ddropn(nd)*sdropn(nd)**3/6.d0 ! new mass
      ddm  = dmo - dmn
c---- store evaporated mass for each nozzle (kg)
      nc = ncelln(nd)
      IF (nozzln(nd).eq.1) then
        pasphi(nc+iMEVAP1) = pasphi(nc+iMEVAP1) + ddm
      else if (nozzln(nd).eq.2) then
        pasphi(nc+iMEVAP2) = pasphi(nc+iMEVAP2) + ddm
      end if
c
      end subroutine
