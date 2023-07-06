#include "mducUDF.h"
#include "udfScalars.h"
#include "udfMemory.h"
#include "udfPremixed.h"

extern udfStruct udf;

void UDFsetMemoryIndex()
{
  int i, k;

  /* ----- set the default values ----- */
  udf.iCHIR[0]   = udf.iCHIR[1]   = udf.iCHIR[2]   = -1;
  udf.iFuelV[0]  = udf.iFuelV[1]  = udf.iFuelV[2]  = -1;
  udf.iZV_SRC[0] = udf.iZV_SRC[1] = udf.iZV_SRC[2] = -1;
  udf.iZV_EPS[0] = udf.iZV_EPS[1] = udf.iZV_EPS[2] = -1;
  udf.iFDT[0]    = udf.iFDT[1]    = udf.iFDT[2]    = -1;
  udf.iHRrate    = -1;

  /* premixed variables */
  udf.iSL        = udf.iST        = udf.iCrif      = -1;
  udf.iRegime    = udf.iProgSrc   = udf.iPscale    = -1;
  udf.iHRprem    = udf.iRIFSrc    = -1;
  udf.iDiffSc    = -1;
  udf.iProgMax   = udf.iProgMin   = -1;

  /* SDM debug variables */
  udf.iSPALD       = udf.iReTURB                         = -1;
  udf.iXI[0]       = udf.iXI[1]       = udf.iXI[2]       = -1;
  udf.iAlphaB[0]   = udf.iAlphaB[1]   = udf.iAlphaB[2]   = -1;
  udf.iZV_PRD[0]   = udf.iZV_PRD[1]   = udf.iZV_PRD[2]   = -1;
  udf.iZV_D2[0]    = udf.iZV_D2[1]    = udf.iZV_D2[2]    = -1;
  udf.iZV_SRC_R[0] = udf.iZV_SRC_R[1] = udf.iZV_SRC_R[2] = -1;
  udf.iZV_EPS_R[0] = udf.iZV_EPS_R[1] = udf.iZV_EPS_R[2] = -1;
  udf.iZV_SRC_L[0] = udf.iZV_SRC_L[1] = udf.iZV_SRC_L[2] = -1;

  /* ----- assign memory as necessary ----- */
  udf.nUDM = 0;
  udf.udmNames = (char**) calloc(1, sizeof(char*));
  k = 0;

  /* memory necessary for fuel mixture fraction (1) */
  if (udf.iZ1 >= 0) {
    udf.nUDM += 6;
    udf.udmNames = (char**) realloc(udf.udmNames, udf.nUDM*sizeof(char*));

    k = UDFassignMem(k, "HR-rate",         &udf.iHRrate,    &udf.udmNames[k]);
    k = UDFassignMem(k, "fuel-volume-1",   &udf.iFuelV[1],  &udf.udmNames[k]);
    k = UDFassignMem(k, "chi-ref-1",       &udf.iCHIR[1],   &udf.udmNames[k]);
    k = UDFassignMem(k, "Z1-variance-src", &udf.iZV_SRC[1], &udf.udmNames[k]);
    k = UDFassignMem(k, "Z1-variance-eps", &udf.iZV_EPS[1], &udf.udmNames[k]);
    k = UDFassignMem(k, "SDM-fdt-1",       &udf.iFDT[1],    &udf.udmNames[k]);
  }

  /* memory necessary for fuel mixture fraction (1) */
  if (udf.iZ2 >= 0) {
    udf.nUDM += 5;
    udf.udmNames = (char**) realloc(udf.udmNames, udf.nUDM*sizeof(char*));

    k = UDFassignMem(k, "fuel-volume-2",   &udf.iFuelV[2],  &udf.udmNames[k]);
    k = UDFassignMem(k, "chi-ref-2",       &udf.iCHIR[2],   &udf.udmNames[k]);
    k = UDFassignMem(k, "Z2-variance-src", &udf.iZV_SRC[2], &udf.udmNames[k]);
    k = UDFassignMem(k, "Z2-variance-eps", &udf.iZV_EPS[2], &udf.udmNames[k]);
    k = UDFassignMem(k, "SDM-fdt-2",       &udf.iFDT[2],    &udf.udmNames[k]);
  }

  /* memory necessary for fuel mixture fraction (0) */
  if (udf.iZ0V >= 0) {
    udf.nUDM += 4;
    udf.udmNames = (char**) realloc(udf.udmNames, udf.nUDM*sizeof(char*));

    k = UDFassignMem(k, "chi-ref-0",       &udf.iCHIR[0],   &udf.udmNames[k]);
    k = UDFassignMem(k, "Z0-variance-src", &udf.iZV_SRC[0], &udf.udmNames[k]);
    k = UDFassignMem(k, "Z0-variance-eps", &udf.iZV_EPS[0], &udf.udmNames[k]);
    k = UDFassignMem(k, "SDM-fdt-0",       &udf.iFDT[0],    &udf.udmNames[k]);
  }

  /* if progress variable is present, assume premixed model is turned on */
  if (udf.iProg >= 0) {
    udf.nUDM += 11;
    udf.udmNames = (char**) realloc(udf.udmNames, udf.nUDM*sizeof(char*));

    k = UDFassignMem(k, "comb-regime",      &udf.iRegime,  &udf.udmNames[k]);
    k = UDFassignMem(k, "prog-src",         &udf.iProgSrc, &udf.udmNames[k]);
    k = UDFassignMem(k, "rif-src",          &udf.iRIFSrc,  &udf.udmNames[k]);
    k = UDFassignMem(k, "pressure-scale",   &udf.iPscale,  &udf.udmNames[k]);
    k = UDFassignMem(k, "lam-flame-speed",  &udf.iSL,      &udf.udmNames[k]);
    k = UDFassignMem(k, "turb-flame-speed", &udf.iST,      &udf.udmNames[k]);
    k = UDFassignMem(k, "prog-rif",         &udf.iCrif,    &udf.udmNames[k]);
    k = UDFassignMem(k, "HR-premixed",      &udf.iHRprem,  &udf.udmNames[k]);
    k = UDFassignMem(k, "progMax",          &udf.iProgMax, &udf.udmNames[k]);
    k = UDFassignMem(k, "progMin",          &udf.iProgMin, &udf.udmNames[k]);
    k = UDFassignMem(k, "Dt-scaling",       &udf.iDiffSc,  &udf.udmNames[k]);
  }

#if DEBUG_SDM
  /* Spalding number and cell turbulent Reynolds number */
  udf.nUDM += 2;
  udf.udmNames = (char**) realloc(udf.udmNames, udf.nUDM*sizeof(char*));

  k = UDFassignMem(k, "spalding-number",  &udf.iSPALD[1],  &udf.udmNames[k]);
  k = UDFassignMem(k, "cell-reynolds-no", &udf.iReTurb[1], &udf.udmNames[k]);

  if (udf.iZ1 >= 0) {
    udf.nUDM += 7;
    udf.udmNames = (char**) realloc(udf.udmNames, udf.nUDM*sizeof(char*));
    
    k = UDFassignMem(k,"SDM-xi-1",          &udf.iXI[1],      &udf.udmNames[k]);
    k = UDFassignMem(k,"SDM-alpha-1",       &udf.iAlphaB[1],  &udf.udmNames[k]);
    k = UDFassignMem(k,"Z1-variance-prod",  &udf.iZV_PRD[1],  &udf.udmNames[k]);
    k = UDFassignMem(k,"Z1-variance-d2",    &udf.iZV_D2[1],   &udf.udmNames[k]);
    k = UDFassignMem(k,"Z1-variance-src-r", &udf.iZV_SRC_R[1],&udf.udmNames[k]);
    k = UDFassignMem(k,"Z1-variance-eps-r", &udf.iZV_EPS_R[1],&udf.udmNames[k]);
    k = UDFassignMem(k,"Z1-variance-src-l", &udf.iZV_SRC_L[1],&udf.udmNames[k]);
  }

  if (udf.iZ2 >= 0) {
    udf.nUDM += 7;
    udf.udmNames = (char**) realloc(udf.udmNames, udf.nUDM*sizeof(char*));

    k = UDFassignMem(k,"SDM-xi-2",          &udf.iXI[2],      &udf.udmNames[k]);
    k = UDFassignMem(k,"SDM-alpha-2",       &udf.iAlphaB[2],  &udf.udmNames[k]);
    k = UDFassignMem(k,"Z2-variance-prod",  &udf.iZV_PRD[2],  &udf.udmNames[k]);
    k = UDFassignMem(k,"Z2-variance-d2",    &udf.iZV_D2[2],   &udf.udmNames[k]);
    k = UDFassignMem(k,"Z2-variance-src-r", &udf.iZV_SRC_R[2],&udf.udmNames[k]);
    k = UDFassignMem(k,"Z2-variance-eps-r", &udf.iZV_EPS_R[2],&udf.udmNames[k]);
    k = UDFassignMem(k,"Z2-variance-src-l", &udf.iZV_SRC_L[2],&udf.udmNames[k]);
  }

  if (udf.iZ0V >= 0) {
    udf.nUDM += 7;
    udf.udmNames = (char**) realloc(udf.udmNames, udf.nUDM*sizeof(char*));

    k = UDFassignMem(k,"SDM-xi-0",          &udf.iXI[0],      &udf.udmNames[k]);
    k = UDFassignMem(k,"SDM-alpha-0",       &udf.iAlphaB[0],  &udf.udmNames[k]);
    k = UDFassignMem(k,"Z0-variance-prod",  &udf.iZV_PRD[0],  &udf.udmNames[k]);
    k = UDFassignMem(k,"Z0-variance-d2",    &udf.iZV_D2[0],   &udf.udmNames[k]);
    k = UDFassignMem(k,"Z0-variance-src-r", &udf.iZV_SRC_R[0],&udf.udmNames[k]);
    k = UDFassignMem(k,"Z0-variance-eps-r", &udf.iZV_EPS_R[0],&udf.udmNames[k]);
    k = UDFassignMem(k,"Z0-variance-src-l", &udf.iZV_SRC_L[0],&udf.udmNames[k]);
  }
#endif DEBUG_SDM

  Message0("  Memory for MDUC library (%i):\n", udf.nUDM);
  for (i=0; i<udf.nUDM; i++)
    Message0("   %3i %32s\n", i, udf.udmNames[i]);
}
/* ------------------------------------------------------------------------- */

int UDFassignMem(int index, char *id, int *var, char **name)
{
  var[0]  = index;
  name[0] = (char*) calloc(32, sizeof(char));
  strcpy(name[0], id);
  return index+1;
}
/* ------------------------------------------------------------------------- */
