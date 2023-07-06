c=======================================================================
      SUBROUTINE usecps(mat,mph)
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
      USE prec_mod, ONLY : prec
      USE comm0,    ONLY : nsp,nep,ncell,nbfac,nreg,itst,i_usecps,iampro
      USE comm1,    ONLY : cpcof,cpcofb,t,ner,nsr,ibc
      USE cthmod,   ONLY : numscl,ymfrac,ymfrab
c      USE mduc,     ONLY : MDUCgetNspecies,indexY
c----- 
      IMPLICIT NONE
c-----
      INTEGER, INTENT(IN) :: mat,mph
      INTEGER :: nc,nb,ir,ib,i,nSpecies
      REAL(prec), ALLOCATABLE, DIMENSION(:) :: h_tmp, cp_tmp
c-----
c      IF (itst < 2) RETURN
      IF (i_usecps == 1) THEN
        write(*,*) iampro, ': in usecps'
      END IF
c-----
c$$$      nSpecies = MDUCgetNspecies()
c$$$      IF(.NOT.ALLOCATED(h_tmp))  ALLOCATE(h_tmp(nSpecies))
c$$$      IF(.NOT.ALLOCATED(cp_tmp)) ALLOCATE(cp_tmp(nSpecies))
c$$$c-----
c$$$      print*, 'calculating cps'
c$$$      DO nc=nsp(mat),nep(mat)
c$$$c----- call the thermo mechanism from MDUC
c$$$        CALL MDUCcomputeThermoData(h_tmp, cp_tmp, t(nc))
c$$$        cpcof(nc) = 0.0
c$$$        DO i=1,numscl
c$$$          IF (indexY(i) .GE. 0) THEN
c$$$            cpcof(nc) = cpcof(nc)
c$$$     x           + ymfrac(nc+(i-1)*ncell)*cp_tmp(indexY(i))
c$$$          END IF
c$$$        END DO
c$$$      END DO
c$$$
c$$$      DO ir=0,nreg
c$$$        IF (ibc(2,ir) == mat) THEN
c$$$          DO ib=nsr(ir),ner(ir)
c$$$c----- call the thermo mechanism from MDUC
c$$$            CALL MDUCcomputeThermoData(h_tmp, cp_tmp, t(nc))
c$$$            cpcofb(nc) = 0.0
c$$$            DO i=1,numscl
c$$$              IF (indexY(i) .GE. 0) THEN
c$$$                cpcofb(nc) = cpcofb(nc)
c$$$     x               + ymfrab(ib+(i-1)*nbfac)*cp_tmp(indexY(i))
c$$$              END IF
c$$$            END DO
c$$$          END DO
c$$$        END IF
c$$$      END DO
c$$$c-----
c$$$      IF(ALLOCATED(h_tmp))  DEALLOCATE(h_tmp)
c$$$      IF(ALLOCATED(cp_tmp)) DEALLOCATE(cp_tmp)
c-----
      END SUBROUTINE usecps
