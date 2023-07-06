#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mducChem.h"
#include "parser.h"
#include "chemMech.h"
#include "chemTable.h"
#include "chemInit.h"
#include "chemUtils.h"
#include "init.h"
#include "fileIO.h"
#include "chemKinIO.h"
#include "ceaIO.h"
#include "chemTable.h"
#include "solveUtils.h"
#include "utils.h"
#include "display.h"

void printHelp()
{
  printf("\n");
  printf("options :\n");
  printf("  -h                                      # print help\n");
  printf("  -i  <input file>                        "
         "# read info from input file\n");
  printf("  -m  <mech file>                         # chemical mechanism\n");
  printf("  -t  <thermo file>                       "
         "# file with thermodynamic data\n");
  printf("  -tr <trans file>                        "
         "# file with transport data\n");
  printf("  -f  <format>                            # mechanism format\n");
  printf("commands :\n");
  printf("  --info                                  # write mechanism info\n");
  printf("  --convert <type>                        "
         "# convert mechanism format\n");
  printf("  --chemtab-info <file> <type> <dim>      "
         "# print chemtable header\n");
  printf("  --chemtab-lookup <file> <x_i>           "
         "# lookup from chemtable at x_i\n");
  printf("  --chemtab-lookup-expansion <file> <x_i> "
         "# lookup from expansion table at x_i\n");
  printf("  --chemtab-PT\n");
  printf("  --cea\n");
  exit(0);
}
/* -------------------------------------------------------------------------- */

void setDefaultOptions(optStruct *opts)
{
  strcpy(opts->inputFile, "");
  strcpy(opts->mechFile,  "chem.dat");
  strcpy(opts->mechBase,  "chem");
  strcpy(opts->thermFile, "therm.dat");
  strcpy(opts->transFile, "trans.dat");
  strcpy(opts->format,    "chemkin");
  opts->convert = 0;
  opts->info    = 0;

  /* chemtable options */
  strcpy(opts->chemTabType, "premixed");
  opts->chemTable        = 0;
  opts->chemTabInfo      = 0;
  opts->chemTabDim       = 4;
  opts->chemTabLookup    = 0;
  opts->chemTabLookupExp = 0;
  opts->cea = 0;
}
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  int i, j, k;
  FILE *pF = NULL;
  void *tab,*tabExp;

  /* ----- set options ----- */
  optStruct opts; /* struct to store options from command line */
  setDefaultOptions(&opts);
  readOptions(argc, argv, &opts);
  MDUC_PRINT_PROC = 1;

  /* ----- options that don't require input ----- */
  if (opts.cea) {
    fprintf(stderr, "reading CEA\n");
    chemStruct *chem = (chemStruct*) calloc(1, sizeof(chemStruct));
    ioStruct io;
    io.parser = parseFile("mduc.inp");
    rxnStruct *rxn;
    rxn = paramGetReaction("global-reaction",io.parser);
    fprintf(stderr, "global-reaction:\n");
    writeChemKinRxn(stderr, rxn[0]);
    chemInitPeriodicTable(chem);
    CEAreadThermo(opts.thermFile,chem);
    CEAwriteThermo("cea_therm_check.dat",chem);
    return 0;
  }
  if (opts.chemTabInfo) {
    tab = loadChemTable(opts.chemTabFile,opts.chemTabDim,opts.chemTabType);
    printChemTableInfo(tab);
    return 0;
  }
  if (opts.chemTabLookup) {
    /* load the chemtable and display the info */
    tab = loadChemTable(opts.chemTabFile, opts.chemTabDim, "premixed");
    printChemTableInfo(tab);

    /* allocate some memory for display */
    int nV = getChemTabNvar(tab);
    double v0, *vLook = (double*) calloc(nV, sizeof(double));
    char **vNames = (char**) calloc(nV, sizeof(char*));
    for (i=0; i<nV; i++)
      vNames[i] = (char*) calloc(64, sizeof(char));
    getChemTabVarNames(vNames, tab);

    /* look up the interpolated value */
    v0 = lookupChemTable("all", opts.chemTabX, vLook, tab);

    /* display the results */
    fprintf(stderr, "\n Values at at: (");
    for (i=0; i<opts.chemTabDim; i++) {
      fprintf(stderr, "%.5f", opts.chemTabX[i]);
      if (i < opts.chemTabDim-1) fprintf(stderr, ", ");
    }
    fprintf(stderr, ")\n");
    for (i=0; i<nV; i++) {
      fprintf(stderr, "%3s%-25s %15.6e\n", "", vNames[i], vLook[i]);
    }
    return 0;
  }
  if (opts.chemTabLookupExp) {
    /* load the chemtable and display the info */
    tab = loadChemTable(opts.chemTabFile, opts.chemTabDim, "pressure");
    printChemTableInfo(tab);

    /* allocate some memory for display */
    int nV = getChemTabNvar(tab);
    double v0, *vLook = (double*) calloc(nV, sizeof(double));
    char **vNames = (char**) calloc(nV, sizeof(char*));
    for (i=0; i<nV; i++)
      vNames[i] = (char*) calloc(64, sizeof(char));
    getChemTabVarNames(vNames, tab);

    /* look up the interpolated value */
    v0 = lookupChemTable("all", opts.chemTabX, vLook, tab);

    /* display the results */
    fprintf(stderr, "\n Values at at: (");
    for (i=0; i<opts.chemTabDim; i++) {
      fprintf(stderr, "%.5f", opts.chemTabX[i]);
      if (i < opts.chemTabDim-1) fprintf(stderr, ", ");
    }
    fprintf(stderr, ")\n");
    for (i=0; i<nV; i++) {
      fprintf(stderr, "%3s%-25s %15.6e\n", "", vNames[i], vLook[i]);
    }
    return 0;
  }
  if (opts.chemTable) {
    chemTable_fitP(&opts);
    return 0;
  }

  /* ----- check the filenames ----- */
  pF = fopen(opts.mechFile, "r");
  if (!pF) {
    fprintf(stderr, "Error: cannot open mechanism file %s\n", opts.mechFile);
    exit(2);
  }
  else {
    fclose(pF);
  }
  /* thermo file*/
  pF = fopen(opts.thermFile, "r");
  if (!pF) {
    /* if default, check to see if base name is intended */
    if (!strcmp(opts.thermFile, "therm.dat")) {
      strcpy(opts.thermFile, opts.mechBase);
      strcat(opts.thermFile, ".therm");
      pF = fopen(opts.thermFile, "r");
    }
    if (!pF) {
      fprintf(stderr, "Error: cannot open thermo file %s\n", opts.thermFile);
      exit(2);
    }
  }
  else {
    fclose(pF);
  }
  /* transport file*/
  pF = fopen(opts.transFile, "r");
  if (!pF) {
    /* if default, check to see if base name is intended */
    if (!strcmp(opts.transFile, "trans.dat")) {
      strcpy(opts.transFile, opts.mechBase);
      strcat(opts.transFile, ".trans");
      pF = fopen(opts.transFile, "r");
    }
    if (!pF) {
      fprintf(stderr, "Error: cannot open transport file %s\n", opts.transFile);
      exit(2);
    }
  }
  else {
    fclose(pF);
  }
  /* try to determine the mechanism format from the suffix */
  if (!strcmp(opts.mechSuffix, "mech"))
    strcpy(opts.format, "chemkin");
  else if (!strcmp(opts.mechSuffix, "flm"))
    strcpy(opts.format, "flamemaster");

  /* ----- execute the desired option ----- */
  if (opts.info) {
    ioStruct   io;
    io.parser = parseFile(opts.inputFile);
    chemStruct *chem = chemInit(opts.mechBase,io.parser);
  }
  if (opts.convert) {
    convertMechanism(&opts);
  }

  return 0;
}
/* -------------------------------------------------------------------------- */

void convertMechanism(optStruct *opts)
{
  if (!strcmp(opts->format, opts->convFormat)) {
    fprintf(stderr,"Error: mechanism input and conversion formats are equal\n");
    exit(2);
  }

  char buff[256];
  fprintf(stderr, " ===== converting mechanism ===== \n");
  fprintf(stderr, "  mechanism:      %s\n", opts->mechFile);  
  fprintf(stderr, "  input format:   %s\n", opts->format);
  fprintf(stderr, "  output format:  %s\n", opts->convFormat);

  /* ----- read in the mechanism ----- */
  ioStruct io;
  io.parser = parseFile(opts->inputFile);
  strcpy(buff, strtok(opts->mechFile, "."));

  /* initialize the chemistry */
  chemStruct *chem = chemInit(buff,io.parser);
  chem->readTrans = 1;

  /* convert from FlameMaster format */
  if ( !strcmp(opts->convFormat, "C") ) {
    /* ----- write C compiled mechansim ----- */
    writeMechCompiledC(chem);
  }
  else if ( !strcmp(opts->convFormat, "chemkin") ) {
    /* ----- write the mechanism file ----- */
    sprintf(buff, "convert_%s.mech", chem->mechName);
    writeChemKinMech(buff, chem);
    
    /* ----- write the thermo data ----- */
    sprintf(buff, "convert_%s.therm", chem->mechName);
    writeChemKinThermo(buff, chem);
    
    /* ----- write the transport data ----- */
    sprintf(buff, "convert_%s.trans", chem->mechName);
    writeChemKinTrans(buff, chem);
  }
  else if ( !strcmp(opts->format, "flamemaster") ) {

  }
  else {
    fprintf(stderr, "Error: conversion format not recognized\n");
  }
  fprintf(stderr, " ================================\n");
}
/* -------------------------------------------------------------------------- */

int chemTable_fitP(optStruct *opts)
{
  fprintf(stderr, "----- Fitting premixed flamelets for variable pressure\n");
  FILE *dFile;
  int nt, i, j, k, iSp, nP, nTu, nPhi, iV, indx, newP, newT, newPhi;
  int nData = opts->nFiles;
  char fileName[512], tmpSp[64];

  /* least squares fitting info */
  int nDat, nPoly, *ia;
  double *xDat, *yDat, *sig, *pCoeff, **covar, chisq;

  // allocate pointers to internal structs
  mducStruct  *mem    = (mducStruct*)  calloc(1, sizeof(mducStruct));
  ioStruct    *io     = (ioStruct*)    calloc(nData, sizeof(ioStruct));
  driveStruct *driver = (driveStruct*) calloc(nData, sizeof(driveStruct));
  timeStruct  *time   = (timeStruct*)  calloc(nData, sizeof(timeStruct));
  flowStruct  *flow   = (flowStruct*)  calloc(nData, sizeof(flowStruct));
  flameStruct *flame  = (flameStruct*) calloc(nData, sizeof(flameStruct));

  /* initialize chemical mechanism */
  io[0].parser = parseFile(opts->inputFile);
  chemStruct *chem = chemInit("compiled",io[0].parser);

  /* allocate the info for the fitting */
  double *P_raw   = (double*)  calloc(nData, sizeof(double));
  double *T_raw   = (double*)  calloc(nData, sizeof(double));
  double *phi_raw = (double*)  calloc(nData, sizeof(double));
  double *P_range = (double*)  calloc(nData, sizeof(double));
  double *T_unb   = (double*)  calloc(nData, sizeof(double));
  double *H_unb   = (double*)  calloc(nData, sizeof(double));
  double *S_L     = (double*)  calloc(nData, sizeof(double));
  double *t_f     = (double*)  calloc(nData, sizeof(double));
  double **prog   = (double**) calloc(nData, sizeof(double*));
  double *Y_prog0 = (double*)  calloc(chem->nSpecies, sizeof(double));

  /* progress variable species */
  /* TODO: read in from input file */
  int nProg = 4;
  char progSpecies[4][25] = {"CO", "CO2", "H2O", "H2"};

  int *progIndx = (int*) calloc(nProg, sizeof(int));

  /* loop over all the files and read in the flamelets */
  nP = nTu = nPhi = 0;
  for (nt=0; nt<nData; nt++) {
    fprintf(stderr, "  reading flamelet: %s\n", opts->dataFile[nt]);
    /* ----- set default values ----- */
    setDefaults(&io[nt], &driver[nt], &time[nt], &flame[nt], &flow[nt], NULL);

    /* ----- read in flamelet file ----- */
    readFlamelet(opts->dataFile[nt],
                 &io[nt],&driver[nt],&time[nt],&flame[nt],&flow[nt]);
    
    P_raw[nt]   = flow[nt].P;
    T_raw[nt]   = io[nt].oxT;
    phi_raw[nt] = flame[nt].phi;

    newP = newT = newPhi = 1;
    for (i=0; i<nt; i++) {

      if (fabs(P_raw[i]   - P_raw[nt])   < 1.0e-9) newP   = 0;
      if (fabs(T_raw[i]   - T_raw[nt])   < 1.0e-9) newT   = 0;
      if (fabs(phi_raw[i] - phi_raw[nt]) < 1.0e-9) newPhi = 0;
    }

    if (newP)   nP++;
    if (newT)   nTu++;
    if (newPhi) nPhi++;

    P_range[nt] = flow[nt].P;
    T_unb[nt]   = io[nt].oxT;
    S_L[nt]     = flame[nt].S_L;
    t_f[nt]     = flame[nt].thick;

    /* get the unburnt species concentrations */
    /* TODO: generalize, currently relies on mass fractions being first half
             of species in boundary condition */
    for (i=0; i<io[nt].nOx/2; i++) {
      iSp = getSpeciesIndex(chem->nSpecies, chem->species, io[nt].oxName[i]);
      Y_prog0[iSp] = io[nt].oxY[i];
    }

    /* compute the unburnt enthalpy */
    H_unb[nt] = 0.0;
    chemComputeThermoData(chem->h, chem->cp, T_unb[nt], chem);
    for (i=0; i<chem->nSpecies; i++)
      H_unb[nt] += Y_prog0[i]*chem->h[i];

    /* compute the progress variable */
    prog[nt] = (double*) calloc(flame[nt].n1, sizeof(double));
    for (j=0; j<flame[nt].n1; j++) {
      for (i=0; i<io[nt].nInputs; i++) {
        for (k=0; k<nProg; k++) {
          if (!strcmp(progSpecies[k], io[nt].inputSpecies[i]))
            prog[nt][j] += io[nt].initY[i][j];
        } /* end loop over progress variable species */
      } /* end loop over input species */
    } /* end loop over grid points */

  } /* end loop over files */

  fprintf(stderr, "nP:   %i\n", nP);
  fprintf(stderr, "nTu:  %i\n", nTu);
  fprintf(stderr, "nPhi: %i\n", nPhi);

  /* return 0;  */
  /* find out the number of unique P, T, phi */


  int nFit = io[0].nInputs;

/*   double ***vFit = (double***) calloc(nData, sizeof(double**)); */
/*   for (nt=0; nt<nData; nt++) { */
/*     vFit[nt] = (double**) calloc(nFit, sizeof(double*)); */

/*     /\* assign the variables *\/ */
/*     iV = 0; */
/*     vFit[nt][iV] = &flame[nt].x1;     iV++; */
/*     vFit[nt][iV] = &io[nt].miscIn[0]; iV++; */
/*     vFit[nt][iV] = &io[nt].initT;     iV++; */
/*     for (i=0; i<io[nt].nInputs; i++, iV++) { */
/*       vFit[nt][iV] = &io[nt].initY[i]; */
/*     } */
/*   } */

  /* define the number of variables to be fit 
     species, temperature, density, TotalEnthalpy, */
  double **Yi_max  = (double**) calloc(nData, sizeof(double*));
  double **var_max = (double**) calloc(nData, sizeof(double*));
  double *T_max    = (double*)  calloc(nData, sizeof(double));

/*   for (j=0; j<nData; j++) { */
/*     Yi_max[j]  = (double*) calloc(io->nInputs, sizeof(double)); */
/*     var_max[j] = (double*) calloc(io->nMisc,   sizeof(double)); */

/*     for (i=0; i<io[nt].nInputs; i++) Yi_max[j][i]  = -1.0E20; */
/*     for (i=0; i<io[nt].nMisc; i++)   var_max[j][i] = -1.0E20; */
/*   } */

  /* find the maximum value and store it's location */

  /* compute the coefficients for each variable */
/*   for (j=0; j<nData; j++) { */
/*     /\* get the max value of the species *\/ */
/*     for (i=0; i<io->nInputs; i++) { */
/*       Yi_max[j][i] = getMaxValue(flame[nt].n1,flame[nt].n2,0,1,io[nt].initY[i]); */
/*     } */
/*   } */

  /* build data arrays for least-squares fitting */
  nP = nData; /* TODO: for now, only varying pressure */
  nPoly = 3;
  /* nPoly = 5; /\* enthalpy fitting *\/ */
  nDat  = flame[0].n1*nData;
  nDat  = nData;

  int refType = 1;  /* 0: all values, 1: maximum values */

  if (refType == 0) {
    nPoly = 5;
    nDat = flame[0].n1*nData;
  }

  /* allocate memory for polynomial least squares fit */
  xDat   = (double*)  calloc(nDat,  sizeof(double));
  yDat   = (double*)  calloc(nDat,  sizeof(double));
  sig    = (double*)  calloc(nDat,  sizeof(double));
  for (nt=0; nt<nDat; nt++) sig[nt] = 1.0;
  pCoeff = (double*)  calloc(nPoly, sizeof(double));
  covar  = (double**) calloc(nPoly, sizeof(double*));
  ia     = (int*)     calloc(nPoly, sizeof(int));
  for (i=0; i<nPoly; i++) {
    covar[i] = (double*) calloc(nPoly, sizeof(double));
    ia[i] = 1;
  }

  double **presF  = (double**) calloc(nFit, sizeof(double*));
  for (i=0; i<nFit; i++) {
    presF[i] = (double*) calloc(nPoly, sizeof(double));
  }

  double **presS  = (double**) calloc(nFit, sizeof(double*));
  for (i=0; i<nFit; i++) {
    presS[i] = (double*) calloc(2, sizeof(double));
  }

  
  double maxVal, maxRef;

  dFile = fopen("pressure_coeffs.dat", "w");
  fprintf(dFile, "# index<1> variable<2> alpha<3> beta<4>\n");
  
  indx = 1;

  /* fit the flame location */
/*   maxRef = getMaxValue(flame[0].n1,flame[0].n2,0,1,flame[0].x1); */
/*   for (nt=0; nt<nData; nt++) { */
/*     maxVal = getMaxValue(flame[nt].n1,flame[nt].n2,0,1,flame[nt].x1); */
/*     xDat[nt] = log(P_range[nt]/P_range[0]); */
/*     yDat[nt] = log(maxVal/maxRef); */
/*     fprintf(stderr, "%15.6e %15.6e\n", xDat[nt], yDat[nt]); */
/*   } */
/*   lfit(xDat,yDat,sig,nDat,pCoeff,ia,nPoly,covar,&chisq,fitPoly); */

/*   fprintf(dFile, "%5i %40s ", indx, "y"); indx++; */
/*   for (k=0; k<nPoly; k++) */
/*     fprintf(dFile, "%15.6e ", pCoeff[k]); */
/*   fprintf(dFile, "\n"); */

/*   /\* dummy for massflowrate *\/ */
/*   fprintf(dFile, "%5i %40s ", indx, "massflowrate"); indx++; */
/*   for (k=0; k<nPoly; k++) */
/*     fprintf(dFile, "%15.6e ", 0.0); */
/*   fprintf(dFile, "\n"); */

/*   /\* fit the temperature  *\/ */
/*   maxRef = getMaxValue(flame[0].n1,flame[0].n2,0,1,io[0].initT); */
/*   for (nt=0; nt<nData; nt++) { */
/*     maxVal = getMaxValue(flame[nt].n1,flame[nt].n2,0,1,io[nt].initT); */
/*     xDat[nt] = log(P_range[nt]/P_range[0]); */
/*     yDat[nt] = log(maxVal/maxRef); */
/*     fprintf(stderr, "%15.6e %15.6e\n", xDat[nt], yDat[nt]); */
/*   } */
/*   lfit(xDat,yDat,sig,nDat,pCoeff,ia,nPoly,covar,&chisq,fitPoly); */

/*   fprintf(dFile, "%5i %40s ", indx, "temperature"); indx++; */
/*   for (k=0; k<nPoly; k++) */
/*     fprintf(dFile, "%15.6e ", pCoeff[k]); */
/*   fprintf(dFile, "\n"); */

  for (i=0; i<io[0].nInputs; i++) {

    /* fill the data according to all points */
    if (refType == 0) {
      for (nt=0; nt<nData; nt++) {
        for (j=0; j<flame[nt].n1; j++) {
          xDat[nt*flame[nt].n1+j] = log(P_range[nt]/P_range[0]);
          yDat[nt*flame[nt].n1+j] = log(io[nt].initY[i][j]/io[0].initY[i][j]);
        }
      }
    }
    else if (refType == 1) { /* use maximum value */
      maxRef = getMaxValue(flame[0].n1,flame[0].n2,0,1,io[0].initY[i]);
      
      for (nt=0; nt<nData; nt++) {
        maxVal = getMaxValue(flame[nt].n1,flame[nt].n2,0,1,io[nt].initY[i]);
        xDat[nt] = log(P_range[nt]/P_range[0]);
        yDat[nt] = log(maxVal/maxRef);
      }
    }

    
/*     for (nt=0; nt<nData; nt++) { */
/*       for (j=0; j<flame[nt].n1; j++) { */
/*         xDat[nt*flame[nt].n1+j] = H_unb[nt]/H_unb[0]; */
/*         yDat[nt*flame[nt].n1+j] = io[nt].initY[i][j]/io[0].initY[i][j]; */
/*       } */
/*     } */

    /* for (nt=0; nt<nDat; nt++) {
/*       fprintf(stderr, "%15.6e %15.6e\n", xDat[nt], yDat[nt]); */
/*     } */

    lfit(xDat,yDat,sig,nDat,pCoeff,ia,nPoly,covar,&chisq,fitPoly);

    /* save the coefficients */
    for (k=0; k<nPoly; k++)
      presF[i][k] = pCoeff[k];
  }
  
  /* dump the coefficients */
  for (i=0; i<io[0].nInputs; i++, indx++) {
    fprintf(dFile, "%5i %40s ", indx, io[0].inputSpecies[i]);
    for (k=0; k<nPoly; k++)
      fprintf(dFile, "%15.6e ", presF[i][k]);
    fprintf(dFile, "\n");
  }
  fclose(dFile);

  /* get the density coefficients */
  double *rhoCoeff = (double*) calloc(nPoly, sizeof(double));
  for (nt=0; nt<nData; nt++) {
    if (refType == 0) {
      for (j=0; j<flame[nt].n1; j++) {
        /* xDat[nt*flame[nt].n1+j] = H_unb[nt]/H_unb[0]; */
        /* yDat[nt*flame[nt].n1+j] = flame[nt].rho[j]/flame[0].rho[j]; */
        xDat[nt*flame[nt].n1+j] = log(P_range[nt]/P_range[0]);
        yDat[nt*flame[nt].n1+j] = log(flame[nt].rho[j]/flame[0].rho[j]);
      }
    }
    else if (refType == 1) {
      maxRef   = getMaxValue(flame[0].n1,flame[0].n2,0,1,flame[0].rho);
      maxVal   = getMaxValue(flame[nt].n1,flame[nt].n2,0,1,flame[nt].rho);
      xDat[nt] = log(P_range[nt]/P_range[0]);
      yDat[nt] = log(maxVal/maxRef);
    }
  }
  lfit(xDat,yDat,sig,nDat,rhoCoeff,ia,nPoly,covar,&chisq,fitPoly);

  /* find the location shift based on progress variable at location of maximum
     value */
  double vScaled;
  double x;
  double xx;
  double maxValSc, maxProg, maxProg_ref;
  double *rhoCoeff_s = (double*) calloc(2, sizeof(double));

  /* get the reference value */
  maxProg_ref = getMaxValue(flame[0].n1,1,0,1,prog[0])/exp(rhoCoeff[0]);
  for (nt=0; nt<nData; nt++) { /* FIXME: only for refType 1 */
    maxProg = getMaxValue(flame[nt].n1,1,0,1,prog[nt]);

    x  = log(P_range[nt]/P_range[0]);
    xx = 1.0;
    for (k=0; k<nPoly; k++) {
      vScaled += rhoCoeff[k]*xx;
      xx *= x;
    }
    xDat[nt] = x;
    yDat[nt] = log(maxProg/exp(vScaled)/maxProg_ref);
  }
  lfit(xDat,yDat,sig,nDat,rhoCoeff_s,ia,2,covar,&chisq,fitPoly);
  

  fprintf(stderr,"density coefficients: %15.6e %15.6e %15.6e\n",
          rhoCoeff[0],rhoCoeff[1],rhoCoeff[2]);
  fprintf(stderr,"                      %15.6e %15.6e\n",
          rhoCoeff_s[0],rhoCoeff_s[1]);

  dFile = fopen("test_output.dat", "w");
  fprintf(dFile,"# pressure[bar]<1> temperature[K]<2> S_L[m/s]<3> t_f[m]<4> H_unb[J/kg]<5>\n");
  for (nt=0; nt<nData; nt++) {
    fprintf(dFile, "%15.6e %15.6e %15.6e %15.6e %15.6e\n", P_range[nt],
            T_unb[nt], S_L[nt], t_f[nt], H_unb[nt]);
  }
  fclose(dFile);

  /* write out the data */
  int ii;
  double *prog_scaled = (double*) calloc(flame[0].n1, sizeof(double));

  strcpy(fileName, "progress_variable_data.dat");
  dFile = fopen(fileName, "w");
  fprintf(dFile, "# prog_var[-]<1> density[kg/m^3]<2> scaled[]<3>\n");
  for (nt=0; nt<nData; nt++) {
    /* sprintf(fileName, "output_data_%04i", nt); */
    x = log(P_range[nt]/P_range[0]);

    /* compute the scaled progress variable */
    /* FIXME: this doesn't work yet */
    prog_scaled = 0;
    for (j=0; j<flame[nt].n1; j++) {
      for (i=0; i<io[nt].nInputs; i++) {
        for (k=0; k<nProg; k++) {
          if (!strcmp(progSpecies[k], io[nt].inputSpecies[i])) {
            xx = 1.0;
            for (ii=0; ii<nPoly; ii++) {
              prog_scaled[j] += rhoCoeff[ii]*xx;
              xx *= x;
            }
          }
        } /* end loop over progress variable species */
      } /* end loop over input species */
    } /* end loop over grid points */

    xx=1.0;
    vScaled=0.0;
    for (k=0; k<nPoly; k++) {
      vScaled += rhoCoeff[k]*xx;
      /* xx *= H_unb[nt]/H_unb[0]; */
      xx *= log(P_range[nt]/P_range[0]);
    }
    for (j=0; j<flame[nt].n1; j++) {
      /* fprintf(dFile, "%15.6e %15.6e\n", prog[nt][j], flame[nt].rho[j]); */
      fprintf(dFile, "%15.6e %15.6e %15.6e\n", prog[nt][j], flame[nt].rho[j],
              flame[0].rho[j]/exp(vScaled));
    }

    fprintf(dFile, "\n\n");
  } /* end loop over files */
  fclose(dFile);

  return 0;
}
/* -------------------------------------------------------------------------- */

void readArray(int n, int iStart, double **var, char **argv)
{
  int i, j;
  var[0] = (double*) calloc(n, sizeof(double));
  for (i=0, j=iStart; i<n; i++, j++)
    var[0][i] = atof(argv[j]);
}
/* -------------------------------------------------------------------------- */

void readOptions(int argc, char **argv, optStruct *opts)
{
  int i, j, k, nOpt, iOpt, iEnd, nTmp;
  int *optLocs;
  char buffer[256];
  
  /* read through the argument list and get the option locations */
  optLocs    = (int*)  calloc(argc, sizeof(int));
  opts->args = (char*) calloc(argc, sizeof(char));
  nOpt = 0;
  for (i=1; i<argc; i++) {
    if (argv[i][0] == '-') {
      optLocs[nOpt] = i;
      opts->args[nOpt] = argv[i][1];
      nOpt++;
    }
  } /* end loop over arguments */
  opts->nOpts = nOpt;
  
  for (i=0; i<nOpt; i++) {
    iOpt = optLocs[i];
    
    if (i == nOpt-1)
      iEnd = argc-1;
    else
      iEnd = optLocs[i+1]-1;
    
    switch (argv[optLocs[i]][1]) {
    case 'h':
      printHelp();
      break;
    case 'i':
      strcpy(opts->inputFile, argv[iOpt+1]);
      break;
    case 'm':
      strcpy(opts->mechFile, argv[iOpt+1]);
      strcpy(buffer, argv[iOpt+1]);
      strcpy(opts->mechBase, strtok(buffer, "."));
      strcpy(opts->mechSuffix, strtok(NULL, "."));
      break;
    case 't':
      if (argv[optLocs[i]][2] == '\0') { /* thermodynamic file */
        strcpy(opts->thermFile, argv[iOpt+1]);
      }
      else if (argv[optLocs[i]][2] == 'r') { /* transport file */
        strcpy(opts->transFile, argv[iOpt+1]);
      }
      break;
    case 'f':
      strcpy(opts->format, argv[iOpt+1]);
      break;
    case 'o':
      strcpy(opts->output, argv[iOpt+1]);
      break;
    case '-': /* has a longer word */
      if (!strcmp(argv[optLocs[i]], "--info")) {
        opts->info = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--convert")) {
        opts->convert = 1;
        strcpy(opts->convFormat, argv[iOpt+1]);
      }
      if (!strcmp(argv[optLocs[i]], "--chemtab-info")) {
        opts->chemTabInfo = 1;
        strcpy(opts->chemTabFile, argv[iOpt+1]);
        if (argv[iOpt+2][0] != '-') /* read in type */
          strcpy(opts->chemTabType, argv[iOpt+2]);
        if (argv[iOpt+3][0] != '-') /* read in dimension */
          opts->chemTabDim = atoi(argv[iOpt+3]);
        /* fprintf(stderr, "%s %s %s\n", opts->chemTabType */
      }
      if ((!strncmp(argv[optLocs[i]], "--chemtab-lookup", 16))) {
        if (!strcmp(argv[optLocs[i]], "--chemtab-lookup"))
          opts->chemTabLookup = 1;
        else
          opts->chemTabLookupExp = 1;
        strcpy(opts->chemTabFile, argv[iOpt+1]);
        opts->chemTabDim = iEnd - (iOpt+1);
        fprintf(stderr, "chemTabDim: %i (%i %i)\n",opts->chemTabDim,iEnd,iOpt);
        readArray(opts->chemTabDim, iOpt+2, &opts->chemTabX, argv);
      }
      if (!strcmp(argv[optLocs[i]], "--chemtab-PT")) {
        opts->chemTable = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--cea")) {
        opts->cea = 1;
      }
      break;
    } /* end switch */
  } /* end loop over options */
  free(optLocs);
}
/* -------------------------------------------------------------------------- */
