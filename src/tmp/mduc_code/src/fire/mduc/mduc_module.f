c=======================================================================
      MODULE mduc
c=======================================================================
c
c     This module defines variables required for coupling the FIRE
c     code with the MDUC combustion module.  It also declares the type
c     for MDUC library functions that return a value.
c
c     NOTE: this module must be compiled before the following user
c     routines that use it:
c        useini.f, usepsc.f, usetim.f
c     Put usemduc.f before these in the list of routines compiled
c
c.....contact emdoran@stanford.edu
c-----
c=======================================================================
c-----
      USE prec_mod, ONLY : prec
      IMPLICIT NONE
c-----
c----- declare variables to be called throughout the programme
c-----
c----- main memory storage
      INTEGER(8) :: mduc_mem
c----- input file
      CHARACTER(len=256) :: mduc_infile
c----- global variable for mpi root process
      INTEGER    :: MPI_COMM_MDUC, rootpro
c----- for source terms
      REAL(prec) :: dPdt,T_drop_ref,H_vap_ref, Pavg_old
c----- reference temperature for enthalpy calculation
      REAL(prec) :: H_Tref
c----- flag for enthalpy version
      INTEGER    :: H_method
c----- flag for number of injections
      INTEGER    :: mduc_nInj
c----- flags for outputing source terms
      INTEGER    :: mduc_outDebug,mduc_outSpraySrc,mduc_outSDM
      INTEGER    :: mduc_outSootSrc,mduc_outSootTP
c----- indices for active scalars
      INTEGER    :: mduc_actoff
      INTEGER    :: iZ1MEAN,iZ1VAR,iZ2MEAN,iZ2VAR,iZ0VAR,iHTOT
      INTEGER    :: bZ1MEAN,bZ1VAR,bZ2MEAN,bZ2VAR,bZ0VAR,bHTOT
c----- indices for passive scalars
      INTEGER    :: mduc_pasoff
      INTEGER    :: iCHI0,iCHI1,iCHI2,iLIQV1,iLIQV2,iMEVAP1,iMEVAP2
      INTEGER    :: iZ0VSRC,iZ1VSRC,iZ2VSRC,iZ0VEPS,iZ1VEPS,iZ2VEPS
      INTEGER    :: iHCHEM,iHRATE,iHS,iCPMIX,iDCPDT,iTCP
      INTEGER    :: bCHI0,bCHI1,bCHI2
      INTEGER    :: bZ0VSRC,bZ1VSRC,bZ2VSRC,bZ0VEPS,bZ1VEPS,bZ2VEPS
      INTEGER    :: iZ0VPR,iZ1VPR,iZ2VPR,bZ0VPR,bZ1VPR,bZ2VPR
      INTEGER    :: iZ0VEV,iZ1VEV,iZ2VEV,bZ0VEV,bZ1VEV,bZ2VEV
      INTEGER    :: iZ0VDI,iZ1VDI,iZ2VDI,bZ0VDI,bZ1VDI,bZ2VDI
      INTEGER    :: bHCHEM,bHRATE,bHS,bCPMIX,bDCPDT,bTCP
c----- injection information
      INTEGER    :: mduc_inj, mduc_mix      ! current injection
      REAL(prec) :: Z1_zero, Z2_zero
      REAL(prec) :: mduc_tmain, mduc_tpilot ! delay times
c----- parameters for SDM model
      INTEGER    :: debug_sdm
      INTEGER    :: iFDT0,iFDT1,iFDT2,iXI0,iXI1,iXI2
      INTEGER    :: iALPH0,iALPH1,iALPH2,iSPALD
      INTEGER    :: bFDT0,bFDT1,bFDT2,bXI0,bXI1,bXI2
      INTEGER    :: bALPH0,bALPH1,bALPH2,bSPALD
c----- flamelet index table of species in solver
      INTEGER    :: mduc_all_species
      INTEGER, ALLOCATABLE, DIMENSION(:) :: indexY !, indexYFav
c----- options for scalar dissipation rate computation
      INTEGER, DIMENSION(3) :: chiForm,chiNorm
c----- variables for soot calculations
      INTEGER    :: mduc_soot,sootSrcUpdate
      INTEGER    :: nSootMom,nSootSpecies
      INTEGER    :: spH,spH2,spOH,spH2O,spO2,spC2H2,spA4
      INTEGER    :: iSOOT_MOM,bSOOT_MOM
      INTEGER    :: iSOOT_M00,iSOOT_M10,iSOOT_M01,iSOOT_N0
      INTEGER    :: bSOOT_M00,bSOOT_M10,bSOOT_M01,bSOOT_N0
      INTEGER    :: iSOOT_M20,iSOOT_M02,iSOOT_M11
      INTEGER    :: bSOOT_M20,bSOOT_M02,bSOOT_M11
      INTEGER    :: iSOOT_N,iSOOT_VOL,iSOOT_AREA
      INTEGER    :: iSOOT_TP,iSOOT_SRC
      INTEGER, DIMENSION(7)           :: sootIndx
      REAL(prec), DIMENSION(7)        :: sootMomZ,sootMomChemSrc
      CHARACTER(len=32)               :: sootSrcUpdateDef
      CHARACTER(len=32), DIMENSION(7) :: sootSpecies
c-----
c----- declarations for MDUC functions that return a variable
c-----
      INTEGER, EXTERNAL    :: MDUCgetNspecies,MDUCgetNSpecs
      INTEGER, EXTERNAL    :: MDUCgetSpeciesIndex
      INTEGER, EXTERNAL    :: MDUCgetNdim,MDUCgetGridSize
      INTEGER, EXTERNAL    :: MDUCcalcPDF,MDUCcalcPDF_1D,MDUCcalcPDF_2D
      INTEGER, EXTERNAL    :: MDUCcalcMeanY_1D,MDUCcalcMeanY_2D
      REAL(prec), EXTERNAL :: MDUCgetPressure,MDUCgetZref,MDUCgetChiRef
      REAL(prec), EXTERNAL :: MDUCgetPointT
      REAL(prec), EXTERNAL :: MDUCgetMolWeight,MDUCcalcTfromH
      REAL(prec), EXTERNAL :: MDUCcalcChiLogForm,MDUCcalcChiErfcForm
c----- strings to read from csf file
      CHARACTER(len=32)  :: H_method_name
      CHARACTER(len=512) :: mduc_satData
c-----
      END MODULE mduc
