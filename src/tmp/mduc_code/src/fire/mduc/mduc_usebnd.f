c=======================================================================
      SUBROUTINE mduc_usebnd(k,mph)
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
      USE prec_mod, ONLY : prec
      USE comm0
      USE comm1
      USE cthmod
      USE mduc
c-----
      IMPLICIT NONE
c-----
      INTEGER, INTENT(in) :: k,mph
c-----
      INTEGER    :: ir,i,ib
      REAL(prec) :: usqr,h_chem,wmix
c-----
c--------------------------------------------------------------------
c-----
      IF (I_USEBND == 1) THEN
        ir=k
        IF (ibc(1,ir)==4) THEN ! wall
          DO ib=nsr(ir),ner(ir)
            usqr    = DOT_PRODUCT(ub(:,ib), ub(:,ib))
            usqr    = usqr*iener
c-----      compute the chemical enthalpy
            h_chem = 0.0
            wmix   = 0.0
            DO i=1,numscl
              h_chem = h_chem + ymolfb(ib+(i-1)*nbfac)*enthfo(i) ! kJ/kmol
              wmix   = wmix   + ymolfb(ib+(i-1)*nbfac)*wmolar(i) ! kg/kmol
            END DO
            h_chem = (h_chem/wmix)*1.d3 ! J/kg
c-----      compute the total enthalpy
            hsb(ib)=(tb(ib)-H_Tref)*cpcofb(ib) + 0.5*usqr + h_chem
          END DO
        END IF
c-----
      ELSE IF (I_USEBND == 2) THEN
c-----
      ELSE IF (I_USEBND == 99) THEN
         CALL USE_FORMULA('USEBND', 2, k, mph, 0)
      END IF
c-----
      RETURN
      END SUBROUTINE mduc_usebnd
c-----------------------------------------------------------------------
c-----
c=======================================================================
      SUBROUTINE CVS_USEBND_F()
c=======================================================================
        CHARACTER(LEN=150) :: CVS_USEBND_F_ID
        DATA CVS_USEBND_F_ID
     1/'@(#) $RCSfile: usebnd.f,v $
     1  $Revision: 1.10 $
     1  $Date: 2010/12/30 15:09:34 $
     1  $Author: vanderma $'/
        WRITE(*,*) CVS_USEBND_F_ID
      END SUBROUTINE CVS_USEBND_F
c-----------------------------------------------------------------------
