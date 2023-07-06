c=======================================================================
      subroutine mduc_cthchk
c=======================================================================
c
c     description: user function for overwriting GUI input for 
c                  species transport
c
c-----------------------------------------------------------------------
c
c     author:            Peter Priesching
c     date:              11-04-06
c     project:           Fire v2009
c     project manager:   Peter Priesching
c
c     literature reference:
c
c-----------------------------------------------------------------------
c
c     modifications:
c
c     author:
c     date:
c     purpose:
c     identification:
c
c-----------------------------------------------------------------------
c
c     argument list:
c
c-----------------------------------------------------------------------
c
c     subroutine is called by: cthdis
c
c-----------------------------------------------------------------------
c
c     subroutine calls:
c
c     name: cthucp            purpose: determine specific heat
c
c=======================================================================
c-----modules
      USE comm0, ONLY: tref_h,iampro
c-----
      IMPLICIT NONE
c-----local dimension
      INTEGER :: len1
      REAL    :: rdum
      CHARACTER(len=128) :: strin
      CHARACTER(len=256) :: text1,text2
c-----------------------------------------------------------------------
c-----
      tref_h=298.15
c-----
      IF (iampro < 2) THEN
        rdum=tref_h
        CALL uustrr(rdum,'*',len1,strin)
        write (*,*) 'rdum   = ', rdum
        write (*,*) 'tref_h = ', tref_h
        text1='Reference Temperature for H-Equation re-set!'
        text2='New value is '//strin(1:len1)//' [K].'
        CALL FIO_UUMESS('USE_CTHCHK','I',text1,text2)
      END IF
c-----
      RETURN
      END SUBROUTINE mduc_cthchk
