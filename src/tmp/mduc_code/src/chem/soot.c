#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "soot.h"
#include "sootHMOM.h"
#include "utils.h"
#include "display.h"

void initSoot(chemStruct *chem)
{
  int i,j,k;
  int counter;
  double NuclNbrC2;
  double NuclNbrH;

  switch (chem->soot) {
  case 1: /* HMOM */

    if (chem->nSootPAH == 0) {  
      fprintf(stderr, "ERROR: No PAH species found in input file\n");
      exit(2);
    }
    nPAH = chem->nSootPAH; 
    chem->sPAH = (int*)  calloc(nPAH, sizeof(int));
    PAHNbrC2 = (double*) calloc(nPAH, sizeof(double));
    PAHNbrH = (double*) calloc(nPAH, sizeof(double));
    PAHStick = (double*) calloc(nPAH, sizeof(double));

    /* ----- initialize the reactions ----- */
    chem->nSootRxn     = 7 + nPAH;
    chem->sSoot	       = chem->nRxnTotal;
    chem->nRxnTotal   += chem->nSootRxn;
    chem->sootRxnRates = (double*)  calloc(chem->nSootRxn, sizeof(double));
   
    counter = 0; 
    for (i=0; i<nPAH; i++ ) {
      /* FM format */
      if (strcmp(chem->sootPAHSpec[i], "A2-C10H8") == 0)
      {
        PAHNbrC2[i] = 5.0;
        PAHNbrH[i] = 8.0;
#ifdef FMSTICK
        PAHStick[i] = 1.0e-3; 
#else
        PAHStick[i] = 2.0e-3; 
#endif /* FMSTICK */
        counter++;
      }
      /* Chemkin format */
      if (strcmp(chem->sootPAHSpec[i], "A2XC10H8") == 0)
      {
        PAHNbrC2[i] = 5.0;
        PAHNbrH[i] = 8.0;
#ifdef FMSTICK
        PAHStick[i] = 1.0e-3; 
#else
        PAHStick[i] = 2.0e-3;
#endif /* FMSTICK */
        counter++;
      }
      /* FM format */
      if (strcmp(chem->sootPAHSpec[i], "A2R5-C12H8") == 0)
      {
        PAHNbrC2[i] = 6.0;
        PAHNbrH[i] = 8.0;
#ifdef FMSTICK
        PAHStick[i] = 3.0e-3; //FM
#else
        PAHStick[i] = 4.0e-3;
#endif /* FMSTICK */
        counter++;
      }
      /* Chemkin format */ 
      if (strcmp(chem->sootPAHSpec[i], "A2R5") == 0)
      {
        PAHNbrC2[i] = 6.0;
        PAHNbrH[i] = 8.0;
#ifdef FMSTICKS
        PAHStick[i] = 3.0e-3; //FM
#else
        PAHStick[i] = 4.0e-3;
#endif /* FMSTICK */
        counter++;
      }
      /* FM format */
      if (strcmp(chem->sootPAHSpec[i], "P2-C12H10") == 0)
      {
        PAHNbrC2[i] = 6.0;
        PAHNbrH[i] = 10.0;
        PAHStick[i] = 8.5e-3;
        counter++;
      }
      /* Chemkin format */
      if (strcmp(chem->sootPAHSpec[i], "P2XC12H10") == 0)
      {
        PAHNbrC2[i] = 6.0;
        PAHNbrH[i] = 10.0;
        PAHStick[i] = 8.5e-3;
        counter++;
      }
      /* FM format */
      if (strcmp(chem->sootPAHSpec[i], "A3-C14H10") == 0)
      {
        PAHNbrC2[i] = 7.0;
        PAHNbrH[i] = 10.0;
        PAHStick[i] = 1.5e-2;
        counter++;
      }
      /* Chemkin format */
      if (strcmp(chem->sootPAHSpec[i], "A3XC14H10") == 0)
      {
        PAHNbrC2[i] = 7.0;
        PAHNbrH[i] = 10.0;
        PAHStick[i] = 1.5e-2;
        counter++;
      }
      /* FM format */
      if (strcmp(chem->sootPAHSpec[i], "A3R5-C16H10") == 0)
      {
        PAHNbrC2[i] = 8.0;
        PAHNbrH[i] = 10.0;
        PAHStick[i] = 2.5e-2;
        counter++;
      }   
      /* Chemkin format */
      if (strcmp(chem->sootPAHSpec[i], "A3R5") == 0)
      {
        PAHNbrC2[i] = 8.0;
        PAHNbrH[i] = 10.0;
        PAHStick[i] = 2.5e-2;
        counter++;
      }   
      /* FM format */
      if (strcmp(chem->sootPAHSpec[i], "A4-C16H10") == 0)
      {
        PAHNbrC2[i] = 8.0;
        PAHNbrH[i] = 10.0;
        PAHStick[i] = 2.5e-2;
        counter++;
      }
      /* Chemkin format */
      if (strcmp(chem->sootPAHSpec[i], "A4XC16H10") == 0)
      {
        PAHNbrC2[i] = 8.0;
        PAHNbrH[i] = 10.0;
        PAHStick[i] = 2.5e-2;
        counter++;
      }
      /* FM format */
      if (strcmp(chem->sootPAHSpec[i], "FLTN-C16H10") == 0)
      {
        PAHNbrC2[i] = 8.0;
        PAHNbrH[i] = 10.0;
        PAHStick[i] = 2.5e-2;
        counter++;
      }
      /* Chemkin format */
      if (strcmp(chem->sootPAHSpec[i], "FLTN") == 0)
      {
        PAHNbrC2[i] = 8.0;
        PAHNbrH[i] = 10.0;
        PAHStick[i] = 2.5e-2;
        counter++;
      }
      /* FM format */
      if (strcmp(chem->sootPAHSpec[i], "A4R5-C18H10") == 0)
      {
        PAHNbrC2[i] = 9.0;
        PAHNbrH[i] = 10.0;
        PAHStick[i] = 3.9e-2;
        counter++;
      }
      /* Chemkin format */
      if (strcmp(chem->sootPAHSpec[i], "A4R5") == 0)
      {
        PAHNbrC2[i] = 9.0;
        PAHNbrH[i] = 10.0;
        PAHStick[i] = 3.9e-2;
        counter++;
      }
    }

    if (counter != nPAH) {
      fprintf(stderr, "ERROR: Wrong PAH species in input file\n");
      exit(2);
    }
    
    NuclNbrC2 = 2.0*PAHNbrC2[0];
    NuclNbrH  = 2.0*PAHNbrH[0];

    for (i = 1; i < nPAH; i++)
    {
      if (NuclNbrC2 > 2.0*PAHNbrC2[i]) NuclNbrC2 = 2.0*PAHNbrC2[i];
      if (NuclNbrH  > 2.0*PAHNbrH[i]) NuclNbrH = 2.0*PAHNbrH[i];

    }

    NuclVol = NuclNbrC2;
    NuclSurf = pow(NuclVol, 2.0/3.0); 

    /* set the species for the soot reactions */
    chem->nSootSp      = nPAH + 7;
    chem->iSootSp      = (int*)     calloc(chem->nSootSp, sizeof(int));
    chem->sootRxnSpec  = (char**)   calloc(chem->nSootSp, sizeof(char*));
    chem->sootRxnConc  = (double*) calloc(chem->nSootSp, sizeof(double));
    for (k=0; k<chem->nSootSp; k++) {
      chem->sootRxnSpec[k] = (char*) calloc(25, sizeof(char));
      /* set reaction indices to -1 for error checking later */
      chem->iSootSp[k] = -1;
    }

    /* assign species names */
    for (i=0; i<nPAH; i++){
      chem->sootRxnSpec[i] = chem->sootPAHSpec[i];
    }
    //chem->sootRxnSpec[0] = "A4-C16H10"; /* pyrene */
    chem->sootRxnSpec[nPAH+0] = "OH";
    chem->sootRxnSpec[nPAH+1] = "H2O";
    chem->sootRxnSpec[nPAH+2] = "H";
    chem->sootRxnSpec[nPAH+3] = "H2";
    chem->sootRxnSpec[nPAH+4] = "C2H2";
    chem->sootRxnSpec[nPAH+5] = "O2";
    chem->sootRxnSpec[nPAH+6] = "CO";

    /* get additional locations for soot (rest are set in initChem) */
    //chem->sA4   = getSpeciesIndex(chem->nSpecies, chem->species, "A4-C16H10");
    for (i=0; i<nPAH; i++) {
      chem->sPAH[i] = getSpeciesIndex(chem->nSpecies, chem->species, chem->sootPAHSpec[i]);
    }
    chem->sC2H2 = getSpeciesIndex(chem->nSpecies, chem->species, "C2H2");

    /* look for the relevant species for the soot reactions */
    for (i=0; i<chem->nSootSp; i++) {
      j = getSpeciesIndex(chem->nSpecies, chem->species, chem->sootRxnSpec[i]);
      if (j < 0) {
        fprintf(stderr, "Error: soot model HMOM requires species %s",
                chem->sootRxnSpec[k]);
        fprintf(stderr, " not found in mechanism\n");
      }
      else
        chem->iSootSp[i] = j;
    }

    /* ----- intitalize the moment information ----- */
    if (chem->sootOrder == 1)
      chem->nSootMom = 4;
    else if (chem->sootOrder == 2)
      chem->nSootMom = 7;
    sootHMOMsetOrder(chem->nSootMom);
    /* Clipping value should be higer than tolerance */
    sootHMOMsetTolerance(1.0E-20);
    //sootHMOMsetTolerance(chem->sootTolerance);

    int foundSel;
    chem->sootMomOn    = (int*)   calloc(chem->nSootMom, sizeof(int));
    chem->sootMomNames = (char**) calloc(chem->nSootMom, sizeof(char*));
    chem->sootMomAbbrv = (char**) calloc(chem->nSootMom, sizeof(char*));
    for (i=0; i<chem->nSootMom; i++) {
      chem->sootMomNames[i] = (char*) calloc(32, sizeof(char));
      chem->sootMomAbbrv[i] = (char*) calloc(32, sizeof(char));
    }
    strcpy(chem->sootMomNames[0], "number-density");
    strcpy(chem->sootMomNames[1], "volume");
    strcpy(chem->sootMomNames[2], "surface");
    strcpy(chem->sootMomNames[chem->nSootMom-1], "delta-weight");
    strcpy(chem->sootMomAbbrv[0], "M00");
    strcpy(chem->sootMomAbbrv[1], "M10");
    strcpy(chem->sootMomAbbrv[2], "M01");
    strcpy(chem->sootMomAbbrv[chem->nSootMom-1], "N0");
    if (chem->nSootMom == 7) {
      strcpy(chem->sootMomNames[3], "volume-var");
      strcpy(chem->sootMomNames[4], "vol-surf-covar");
      strcpy(chem->sootMomNames[5], "surface-var");
      strcpy(chem->sootMomAbbrv[3], "M20");
      strcpy(chem->sootMomAbbrv[4], "M11");
      strcpy(chem->sootMomAbbrv[5], "M02");
    }

    if (chem->nSootMomOn < 0) { /* if undefined, set to default */
      for (i=0; i<chem->nSootMom; i++)
        chem->sootMomOn[i] = 1;
      chem->nSootMomOn = chem->nSootMom;
    }
    else {
      for (i=0; i<chem->nSootMomOn; i++) {
        if (!strcmp(chem->sootMomSelect[i], "default") ||
            !strcmp(chem->sootMomSelect[i], "all")) {
          for (j=0; j<chem->nSootMom; j++)
            chem->sootMomOn[j] = 1;
          chem->nSootMomOn = chem->nSootMom;
          break;
        }
        else {
          foundSel = 0;
          for (j=0; j<chem->nSootMom; j++) {
            if (!strcmp(chem->sootMomSelect[i],chem->sootMomNames[j]) ||
                !strcmp(chem->sootMomSelect[i],chem->sootMomAbbrv[j])) {
              chem->sootMomOn[j] = 1;
              foundSel = 1;
            }
          }
          if (!foundSel)
            fprintf(stderr, "Warning: soot moment not recognized: %s\n",
                    chem->sootMomSelect[i]);
        }
      }
    }

    /* number of soot source term contributions */
    chem->nSootSrc = 6;
    chem->sootSrcOn    = (int*)   calloc(chem->nSootSrc, sizeof(int));
    chem->sootSrcNames = (char**) calloc(chem->nSootSrc, sizeof(char*));
    for (i=0; i<chem->nSootSrc; i++) {
      chem->sootSrcOn[i] = 0;
      chem->sootSrcNames[i] = (char*) calloc(32, sizeof(char));
    }
    strcpy(chem->sootSrcNames[0], "nucleation");
    strcpy(chem->sootSrcNames[1], "condensation");
    strcpy(chem->sootSrcNames[2], "coagulation");
    strcpy(chem->sootSrcNames[3], "surface-growth");
    strcpy(chem->sootSrcNames[4], "oxidation");
    strcpy(chem->sootSrcNames[5], "fragmentation");

    if (chem->nSootSrcOn < 0) { /* if undefined, set to default */
      for (i=0; i<chem->nSootSrc; i++)
        chem->sootSrcOn[i] = 1;
      chem->nSootSrcOn = chem->nSootSrc;
    }
    else {
      for (i=0; i<chem->nSootSrcOn; i++) {
        if (!strcmp(chem->sootSrcSelect[i], "default") ||
            !strcmp(chem->sootSrcSelect[i], "all")) {
          for (j=0; j<chem->nSootSrc; j++)
            chem->sootSrcOn[j] = 1;
          chem->nSootSrcOn = chem->nSootSrc;
          break;
        }
        else {
          foundSel = 0;
          for (j=0; j<chem->nSootSrc; j++) {
            if (!strncmp(chem->sootSrcSelect[i],chem->sootSrcNames[j],4)) {
              chem->sootSrcOn[j] = 1;
              foundSel = 1;
            }
          }
          if (!foundSel)
            fprintf(stderr, "Warning: soot source not recognized: %s\n",
                    chem->sootSrcSelect[i]);
        }
      }
    }

    /* write a summary */
    if (MDUC_PRINT_PROC) fprintf(stderr, "===== Soot Initialization Summary =====\n");
    if (MDUC_PRINT_PROC) fprintf(stderr, "  Soot moments : %i\n", chem->nSootMom);
    if (MDUC_PRINT_PROC) fprintf(stderr, "  Active soot moments : %i\n", chem->nSootMomOn);
    for (i=0; i<chem->nSootMom; i++) {
      if (chem->sootMomOn[i])
        if (MDUC_PRINT_PROC) fprintf(stderr, "    %-3s (%-18s) : active\n", chem->sootMomAbbrv[i],
                chem->sootMomNames[i]);
      else
        if (MDUC_PRINT_PROC) fprintf(stderr, "    %-3s (%-18s) : inactive\n", chem->sootMomAbbrv[i],
                chem->sootMomNames[i]);
    }
    if (MDUC_PRINT_PROC) fprintf(stderr, "  Active soot sources : %i\n", chem->nSootSrcOn);
    for (i=0; i<chem->nSootSrc; i++) {
      if (chem->sootSrcOn[i])
        if (MDUC_PRINT_PROC) fprintf(stderr, "    %-17s : active\n", chem->sootSrcNames[i]);
      else
        if (MDUC_PRINT_PROC) fprintf(stderr, "    %-17s : inactive\n", chem->sootSrcNames[i]);
    }
    if (MDUC_PRINT_PROC) fprintf(stderr, "=======================================\n");

    break;
  default:
    if (MDUC_PRINT_PROC) fprintf(stderr, "Warning: soot model not recognized\n");
  } /* end switch chem->soot */
}
/* -------------------------------------------------------------------------- */
void reinitSootMoments(int nMom, double *moments,double tol)
{
  /* reinitialize the soot moments according to the solver tolerance */
  /* ----- first order moments ----- */
  double loc_tol = 1.0E-20;
//  loc_tol = tol; 
  /* M_00: number density [-] */
  moments[0] = fmax(moments[0],2.0*loc_tol);
  /* M_10: volume fraction [m^3] */
  moments[1] = fmax(moments[1],2.0*loc_tol*pow(NuclVol*2.0,1.0+2.0/3.0*0.0));
  /* M_01: surface area density [m^2] */
  moments[2] = fmax(moments[2],2.0*loc_tol*pow(NuclVol*2.0,0.0+2.0/3.0*1.0));

  /* ----- second order moments ----- */
  if (nMom == 7) {
    /* M_20: volume variance [m^6] */
    moments[3] = fmax(moments[3],loc_tol*pow(NuclVol*2.0,2.0+2.0/3.0*0.0));
    /* M_11: volume-surface covariance [m^5] */
    moments[4] = fmax(moments[4],loc_tol*pow(NuclVol*2.0,1.0+2.0/3.0*1.0));
    /* M_02: surface variance [m^4] */
    moments[5] = fmax(moments[5],loc_tol*pow(NuclVol*2.0,0.0+2.0/3.0*2.0));
  }

  /* N_0: weight of the delta function [-] */
  moments[nMom-1] = fmax(moments[nMom-1],loc_tol);
}

void initSootMoments(int nMom, double *moments, double tol)
{
  /* initialize the soot moments to a small number density, with other
     moments determined by M_xy = M_00*40^(x+2/3*y) */

  /* ----- first order moments ----- */
  /* M_00: number density [-] */
  //moments[0] = 1.0e-20;
  ///* M_10: volume fraction [m^3] */
  //moments[1] = 1.0e-20*pow(NuclVol*2.0,1.0+2.0/3.0*0.0);
  ///* M_01: surface area density [m^2] */
  //moments[2] = 1.0e-20*pow(NuclVol*2.0,0.0+2.0/3.0*1.0);

  ///* ----- second order moments ----- */
  //if (nMom == 7) {
  //  /* M_20: volume variance [m^6] */
  //  moments[3] = 1.0e-20*pow(NuclVol*2.0,2.0+2.0/3.0*0.0);
  //  /* M_11: volume-surface covariance [m^5] */
  //  moments[4] = 1.0e-20*pow(NuclVol*2.0,1.0+2.0/3.0*1.0);
  //  /* M_02: surface variance [m^4] */
  //  moments[5] = 1.0e-20*pow(NuclVol*2.0,0.0+2.0/3.0*2.0);
  //}

  ///* N_0: weight of the delta function [-] */
  //moments[nMom-1] = 1.0e-20;
  /* ----- first order moments ----- */
  double loc_tol = 1.0E-20;

  double V0 = NuclVol;
  double S0 = pow(V0,2.0/3.0);
//  double loc_tol=tol;
  /* M_00: number density [-] */
  moments[0] = loc_tol*2.0;
  /* M_10: volume fraction [m^3] */
  moments[1] = loc_tol*2.0*V0;
  /* M_01: surface area density [m^2] */
  moments[2] = loc_tol*2.0*S0;

  /* ----- second order moments ----- */
  if (nMom == 7) {
    /* M_20: volume variance [m^6] */
    moments[3] = loc_tol*pow(NuclVol*2.0,2.0+2.0/3.0*0.0);
    /* M_11: volume-surface covariance [m^5] */
    moments[4] = loc_tol*pow(NuclVol*2.0,1.0+2.0/3.0*1.0);
    /* M_02: surface variance [m^4] */
    moments[5] = loc_tol*pow(NuclVol*2.0,0.0+2.0/3.0*2.0);
  }

  /* N_0: weight of the delta function [-] */
  moments[nMom-1] = loc_tol;
}
/* -------------------------------------------------------------------------- */

void updateSootRxnRates(double *cdot, chemStruct *chem)
{
}
/* -------------------------------------------------------------------------- */
