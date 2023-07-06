c=======================================================================
      SUBROUTINE cyusesu
c=======================================================================
c      Implementation of MDUC via User Routines
c-----------------------------------------------------------------------
c     cyusesu is called in cysetup after reading the spray GUI from
c     routine cyispra. This userfunction allows to set up additional 
c     model activation flags, modify predefined model flags  or read 
c     additional input data from an input file.
c
c     cyusesu is called by subroutine cyscalc
c     cyusesu calls       : -
c-----
c-----------------------------------------------------------------------
c-----
      USE cydrnoz
c-----
      IMPLICIT NONE
c-----
c-----------------------------------------------------------------------
c---- Mixture fraction source terms => cyuseso called by cysourc.f
      i_useso = 1
      RETURN
      END
c=======================================================================
