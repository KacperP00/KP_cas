#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "mducTools.h"
#include "accessMDUC.h"
#include "chemMech.h"
#include "chemInit.h"
#include "accessChem.h"
#include "utils.h"
#include "metric.h"
#include "init.h"
#include "solve2D.h"
#include "init.h"
#include "fileIO.h"
#include "vtkIO.h"
#include "parser.h"
#include "destroy.h"
#include "statsUtils.h"
#include "postProcess.h"
#include "errFunc.h"
#include "solveUtils.h"
#include "scalarDiss.h"
#include "polyInt.h"
#include "engineUtils.h"
#include "display.h"

double calcLiquidFuelTemp(double T_L, double h_v, double cp_L, double T_ref,
                          int nSpecies, int nFuel, int *fuelIndex, double *Y,
                          chemStruct *chem)
{
  int    i, count = 0;
  double hv_fuel;  // enthalpy of pure fuel vapour
  double hv_ref;   // enthalpy of composition at reference temperature
  double T_fuel = 300;
  double *cp, *h, cpSum, hSum, deltaT, temp;

  cp = (double*) calloc(nSpecies, sizeof(double));
  h  = (double*) calloc(nSpecies, sizeof(double));

  // compute the enthalpy of the vapour at the reference temperature
  chemComputeThermoData(h, cp, T_ref, chem);
  hv_ref = 0.0;
  for (i=0; i<nFuel; i++)
    hv_ref += Y[i]*h[fuelIndex[i]];

  // compute the enthalpy of the vapour at fuel boundary
  hv_fuel = cp_L*(T_L - T_ref) - h_v + hv_ref;

  // compute the temperature at the pure fuel boundary
  do {
    chemComputeThermoData(h, cp, T_fuel, chem);
    hSum = 0.0;
    cpSum = 0.0;
    for (i=0; i<nFuel; i++) {
      cpSum += cp[fuelIndex[i]]*Y[i];
      hSum  +=  h[fuelIndex[i]]*Y[i];
    }
    deltaT = -(hSum - hv_fuel)/cpSum;
    T_fuel += deltaT;
    
    count++;
    if (++count > 1000) {
      fprintf(stderr,"Error: temperature iteration for h = %g not converged\n",
              hv_fuel);
      exit(2);
    }
  } while (fabs(deltaT/T_fuel) > 1.0e-8);
  
  free(cp);
  free(h);

  return T_fuel;
}
/* -------------------------------------------------------------------------- */

void interpolateSolution(optStruct *opts)
{
  int i, j, k, r, row, col;
  double dx;
  char  fileHead[100], interpFile[100], nInfo[30], *head;

  // allocate pointers to internal structs
  ioStruct    *io        = (ioStruct*)    calloc(1, sizeof(ioStruct));
  driveStruct *driver    = (driveStruct*) calloc(1, sizeof(driveStruct));
  timeStruct  *time      = (timeStruct*)  calloc(1, sizeof(timeStruct));
  flowStruct  *flow      = (flowStruct*)  calloc(1, sizeof(flowStruct));
  flameStruct *flame     = (flameStruct*) calloc(1, sizeof(flameStruct));

  flameStruct flameI;
  flowStruct  flowI;

  /* ----- set default values ----- */
  setDefaults(io, driver, time, flame, flow, NULL);
  setDefaults(io, driver, time, &flameI, &flowI, NULL);

  readFlamelet(opts->mducFile, io, driver, time, flame, flow);

  /* ----- initialize necessary parameters from input file read ----- */
  flame->nVars = io->nInputs + 1;
  flame->sT    = io->nInputs;

  io->nOutputs = io->nInputs;
  io->outputSpecies = (char**) calloc(io->nOutputs, sizeof(char*));
  io->inputIndex    = (int*)   calloc(io->nInputs, sizeof(int));
  io->outputIndex   = (int*)   calloc(io->nOutputs, sizeof(int));
  for (i=0; i<io->nOutputs; ++i) {
    io->outputSpecies[i] = (char*) calloc(25, sizeof(char));
    strcpy(io->outputSpecies[i], io->inputSpecies[i]);
    io->outputIndex[i] = io->inputIndex[i] = i;
  }

  // initialise solution vector
  flame->Y = (double*) calloc(flame->nVars*flame->n1*flame->n2, sizeof(double));
  for (k=0; k<flame->n2; ++k) {
    r   = k*flame->n1;
    row = k*flame->n1*flame->nVars;
    for (j=0; j<flame->n1; ++j) {
      col = j*flame->nVars;
      for (i=0; i<io->nInputs; ++i)
        flame->Y[row+col+io->inputIndex[i]] = io->initY[i][r+j];
      
      flame->Y[row+col+flame->sT] = io->initT[r+j];
    } // end loop over columns
  } // end loop over rows
  
  // allocate the memory for the interpolated flamelet
  flameI.n1    = opts->n1;
  if (flame->nDim == 1)
    flameI.n2  = 1;
  else if (flame->nDim == 2)
    flameI.n2  = opts->n2;
  flameI.nDim  = flame->nDim;
  flameI.nVars = flame->nVars;
  flameI.sT    = flame->sT;

  //flameI.x1  = (double*) calloc(flameI.n1, sizeof(double));
  //flameI.x2  = (double*) calloc(flameI.n2, sizeof(double));
  flameI.Y = (double*) calloc(flameI.n1*flameI.n2*flame->nVars,sizeof(double));
  flameI.H   = (double*) calloc(flameI.n1*flameI.n2, sizeof(double));
  flowI.chi1 = (double*) calloc(flameI.n1*flameI.n2, sizeof(double));
  flowI.chi2 = (double*) calloc(flameI.n1*flameI.n2, sizeof(double));

  /* fill the fine mesh (equispaced grid for now) */
  /* dx = (flame->x1[flame->n1-1]-flame->x1[0])/(flameI.n1-1); */
/*   flameI.x1[0] = flame->x1[0]; */
/*   for (j=1; j < flameI.n1; j++) */
/*     flameI.x1[j] = flameI.x1[j-1] + dx; */
  initGrid(flameI.n1,0,1,opts->nPivot,opts->pivot,opts->gridPerc,
           &flameI.x1,&flameI.dx1);

/*   if (flame->nDim == 2) */
/*     { */
/*       dx = (flame->x2[flame->n2-1]-flame->x2[0])/(flameI.n2-1); */
/*       flameI.x2[0] = flame->x2[0]; */
/*       for (j=1; j < flameI.n2; j++) */
/*         flameI.x2[j] = flameI.x2[j-1] + dx; */
/*     } */
  initGrid(flameI.n2,0,1,opts->nPivot,opts->pivot,opts->gridPerc,
           &flameI.x2,&flameI.dx2);

  /* ----- interpolate solution of all variables ----- */
  if (flame->nDim == 1) {
    interpSolution1D(flame->n1, flame->x1, flameI.n1, flameI.x1, flame->nVars,
                     flame->Y, flameI.Y);
    if (flame->H != NULL)
      interpSolution1D(flame->n1, flame->x1, flameI.n1, flameI.x1, 1,
                       flame->H, flameI.H);
    if (flow->chi1 != NULL)
      interpSolution1D(flame->n1, flame->x1, flameI.n1, flameI.x1, 1,
                       flow->chi1, flowI.chi1);
  }
  else if (flame->nDim == 2) {
    interpSolution2D(flame->n1, flame->n2, flame->x1, flame->x2, flameI.n1,
                     flameI.n2, flameI.x1, flameI.x2, flame->nVars,
                     flame->Y, flameI.Y);
    if (flow->chi1 != NULL)
      interpSolution2D(flame->n1, flame->n2, flame->x1, flame->x2, flameI.n1,
                       flameI.n2, flameI.x1, flameI.x2, 1,
                       flow->chi1, flowI.chi1);
    if (flow->chi1 != NULL)
      interpSolution2D(flame->n1, flame->n2, flame->x1, flame->x2, flameI.n1,
                       flameI.n2, flameI.x1, flameI.x2, 1,
                       flow->chi2, flowI.chi2);
  }

  // make filename for interpolated output
  strcpy(fileHead, opts->mducFile);
  head = strtok(fileHead, ".");
  strcpy(interpFile, "interp_");
  strcat(interpFile, head);
  sprintf(nInfo, "_n%i", flameI.n1);
  strcat(interpFile, nInfo);
  strcat(interpFile, ".rif");

  writeFlamelet(interpFile, io, driver, time, &flameI, &flowI);
}
/* -------------------------------------------------------------------------- */

void convert1D_2D(optStruct *opts)
{
  int i, j, k, r, row, col;

  // allocate pointers to internal structs
  ioStruct    *io        = (ioStruct*)    calloc(1, sizeof(ioStruct));
  driveStruct *driver    = (driveStruct*) calloc(1, sizeof(driveStruct));
  timeStruct  *time      = (timeStruct*)  calloc(1, sizeof(timeStruct));
  flowStruct  *flow      = (flowStruct*)  calloc(1, sizeof(flowStruct));
  flameStruct *flame     = (flameStruct*) calloc(1, sizeof(flameStruct));

  /* ----- set default values ----- */
  setDefaults(io, driver, time, flame, flow, NULL);

  readFlamelet(opts->mducFile, io, driver, time, flame, flow);

  // initialise the second grid dimension
  flame->nDim = 2;
  flame->n2   = opts->n2;
  initGrid(flame->n2,0,1,io->nPivot2,io->gPivot2,io->gPerc2,
           &flame->x2,&flame->dx2);

  io->nOutputs = io->nInputs;
  io->outputSpecies = (char**) calloc(io->nOutputs, sizeof(char*));
  io->outputIndex = (int*) calloc(io->nOutputs, sizeof(int));
  for (i=0; i<io->nOutputs; i++) {
    io->outputSpecies[i] = (char*) calloc(25, sizeof(char));
    strcpy(io->outputSpecies[i], io->inputSpecies[i]);
    io->outputIndex[i] = i;
  }
  
  flame->nVars = io->nInputs+1;
  flame->sT    = io->nInputs;

  // allocate memory
  flame->Y   = (double*) calloc(flame->nVars*flame->n1*flame->n2,sizeof(double));
  flame->H   = (double*) realloc(flame->H,flame->n1*flame->n2*sizeof(double));
  flow->chi1 = (double*) realloc(flow->chi1,flame->n1*flame->n2*sizeof(double));

  if (opts->dir == 1) {
    // init the io stuff
    io->nFuel2 = io->nOx;
    io->fuelT2 = io->oxT;
    io->fuelIndex2 = (int*) calloc(io->nFuel2, sizeof(int));
    io->fuelY2     = (double*) calloc(io->nFuel2, sizeof(double));
    io->fuelName2  = (char**) calloc(io->nFuel2, sizeof(char*));
    for (i=0; i<io->nFuel2; i++) {
      // assign all fuel mass fractions and names
      io->fuelY2[i] = io->oxY[i];
      io->fuelName2[i] = (char*) calloc(25, sizeof(char));
      strcpy(io->fuelName2[i], io->oxName[i]);
    }
      
    // copy the row solution to all rows
    for (k=0; k<flame->n2; k++) {
      r   = k*flame->n1;
      row = k*flame->n1*flame->nVars;
      for (j=0; j < flame->n1; j++) {
        col = j*flame->nVars;
        
        for (i=0; i<io->nInputs; i++)
          flame->Y[row+col+i] = io->initY[i][j];//flame->Y[col+i];
        
        flame->Y[row+col+flame->sT] = io->initT[j];
        
        flame->H[r+j]   = flame->H[j];
        
        flow->chi1[r+j] = flow->chi1[j];
      }
    }
    
  } // end all row init
  else if (opts->dir == 2) {
    // init the io stuff
    io->nFuel2 = io->nFuel1;
    io->fuelT2 = io->fuelT1;
    io->fuelIndex2 = (int*) calloc(io->nFuel2, sizeof(int));
    io->fuelY2     = (double*) calloc(io->nFuel2, sizeof(double));
    io->fuelName2  = (char**) calloc(io->nFuel2, sizeof(char*));
    for (i=0; i<io->nFuel2; i++) {
      // assign all fuel mass fractions and names
      io->fuelY2[i] = io->fuelY1[i];
      io->fuelName2[i] = (char*) calloc(25, sizeof(char));
      strcpy(io->fuelName2[i], io->fuelName1[i]);
    }
    
    for (i=0; i<io->nFuel1; i++) free(io->fuelName1[i]);
    free(io->fuelIndex1);
    free(io->fuelY1);
    free(io->fuelName1);
    
    // init the io stuff
    io->nFuel1 = io->nOx;
    io->fuelT1 = io->oxT;
    io->fuelIndex1 = (int*) calloc(io->nFuel1, sizeof(int));
    io->fuelY1     = (double*) calloc(io->nFuel1, sizeof(double));
    io->fuelName1  = (char**) calloc(io->nFuel1, sizeof(char*));
    for (i=0; i<io->nFuel1; i++) {
      // assign all fuel mass fractions and names
      io->fuelY1[i] = io->oxY[i];
      io->fuelName1[i] = (char*) calloc(25, sizeof(char));
      strcpy(io->fuelName1[i], io->oxName[i]);
    }
    
    // transpose the row solution to the first column
    transposeRowCol(flame->n1, flame->n2, 1, flow->chi1);
    
    // switch the grid
    int    tmpN;
    double *tmpGrid;
    
    tmpN    = flame->n1;
    tmpGrid = flame->x1;
    
    flame->n1 = flame->n2;
    flame->x1 = flame->x2;
    flame->n2 = tmpN;
    flame->x2 = tmpGrid;
    
    for (k=0; k<flame->n2; k++) {
      r   = k*flame->n1;
      row = k*flame->n1*flame->nVars;
      for (j=0; j < flame->n1; j++) {
        col = j*flame->nVars;
        
        for (i=0; i<io->nInputs; i++)
          flame->Y[row+col+i] = io->initY[i][k];
        
        flame->Y[row+col+flame->sT] = io->initT[k];
        
        flame->H[r+j]   = flame->H[r];
        
        flow->chi1[r+j] = flow->chi1[r];
      }
    }
  } // end init cols
  
  writeFlamelet(opts->outFile, io, driver, time, flame, flow);
}
/* ------------------------------------------------------------------------- */

void genEnsightTime(optStruct *opts)
{
  int i, nSteps;
  double *tFile, *var;

  nSteps = readMonitorFile(opts->monitorFile,opts->iT,&tFile,opts->iV,&var);

  FILE *pFile = fopen("ensightTime.dat", "w");
  fprintf(pFile, "%i\n", nSteps);

  for (i=0; i<nSteps; i++) {
    fprintf(pFile, "%15.6e ", tFile[i]);
    if (i%5 == 0)
      fprintf(pFile, "\n");
  }
  fclose(pFile);
  free(tFile); free(var);
}
/* ------------------------------------------------------------------------- */

void mduc2ensight(optStruct *opts)
{
  int i, j, nt, iFile, checkTime = 0, tNext = 0;
  char fileName[200];
  double tCur, *timeVals = NULL;
  FILE *pFile;
  //char  tag[50], buffer[200], *ident;

  // allocate pointers to internal structs
  ioStruct    *io        = (ioStruct*)    calloc(1, sizeof(ioStruct));
  driveStruct *driver    = (driveStruct*) calloc(1, sizeof(driveStruct));
  timeStruct  *time      = (timeStruct*)  calloc(1, sizeof(timeStruct));
  flowStruct  *flow      = (flowStruct*)  calloc(1, sizeof(flowStruct));
  flameStruct *flame     = (flameStruct*) calloc(1, sizeof(flameStruct));

  /* if a current ensight case file is specified to match the frequency,
     read in the time array */
  if (strcmp(opts->ensCase, "none")) {
    checkTime = 1;
    nt = readEnsightCaseTime(opts->ensCase, &timeVals);
  }
  
  fprintf(stderr, "Converting data files to Ensight Gold format:\n");
  iFile = 0;
  for (j=0; j<opts->nFiles; j++) {
    // check to see if this file corresponds to the right time
    if (checkTime) {
      if (tNext >= nt) {
        fprintf(stderr, "No more timesteps. Exiting\n");
        break;
      }
      tCur = readFileTime(opts->dataFile[j]);
      if (fabs(tCur - timeVals[tNext]) < 1e-7) {
        /* ----- set default values and read in the flamelet ----- */
        setDefaults(io, driver, time, flame, flow, NULL);
        readFlamelet(opts->dataFile[j], io, driver, time, flame, flow);
        
        tNext++;
        iFile++;
        fprintf(stderr, "\t--> converting %s ...\n", opts->dataFile[j]);
      }
      else if (tCur > timeVals[tNext]) {
        fprintf(stderr, "Error: file time %15.6e > ensight time %15.6e\n",
                tCur, timeVals[tNext]);
        exit(2);
      }
      else {
        fprintf(stderr, "\t--> skipping   %s ...\n", opts->dataFile[j]);
        continue;
      }
    }
    else {
      fprintf(stderr, "\t--> %s ...\n", opts->dataFile[j]);
      /* ----- set default values and read in the flamelet ----- */
      setDefaults(io, driver, time, flame, flow, NULL);
      readFlamelet(opts->dataFile[j], io, driver, time, flame, flow);
      if (j == 0)
        iFile = 1;
    }
    
    strcpy(io->ensightDir, opts->ensDir);
    if (iFile == 1) {
      initEnsight(io, flow);
      
      writeEnsightGeom(io, driver, time, flame, flow);
      
      writeEnsightCaseHead(io, driver, time, flame, flow);
    }
    
    // update the case file
    updateEnsightCaseTime(time->t, io, time);
    
    // write species data
    for (i=0; i<io->nInputs; i++) {
      sprintf(fileName, "%s/%s/%s.%06i", io->ensightDir, io->inputSpecies[i],
              io->inputSpecies[i], iFile);
      writeEnsightData(fileName, io->inputSpecies[i], flame->n2, flame->n1,
                       1, 0, io->initY[i]);
    }
    
    // write temperature
    if (io->initT != NULL) {
      sprintf(fileName, "%s/%s/%s.%06i", io->ensightDir, "T", "T", iFile);
      writeEnsightData(fileName, "T", flame->n2, flame->n1, 1, 0, io->initT);
    }
    
    // write probability distribution
    if (flow->pdf != NULL) {
      sprintf(fileName, "%s/%s/%s.%06i", io->ensightDir,"PDF","PDF",iFile);
      writeEnsightData(fileName, "PDF", flame->n2,flame->n1,1,0,flow->pdf);
    }
    
    // write scalar dissipation rates
    if (flow->chi1 != NULL) {
      sprintf(fileName, "%s/%s/%s.%06i",io->ensightDir,"CHI1","CHI1",iFile);
      writeEnsightData(fileName, "CHI1",flame->n2,flame->n1,1,0,flow->chi1);
    }
    if (flow->chi12 != NULL) {
      sprintf(fileName, "%s/%s/%s.%06i", io->ensightDir,"CHI12","CHI12",iFile);
      writeEnsightData(fileName, "CHI12",flame->n2,flame->n1,1,0,flow->chi12);
    }
    if (flow->chi2 != NULL) {
      sprintf(fileName, "%s/%s/%s.%06i",io->ensightDir,"CHI2","CHI2",iFile);
      writeEnsightData(fileName, "CHI2",flame->n2,flame->n1,1,0,flow->chi2);
    }

    // write all other variables
    for (i=0; i<io->nMisc; i++) {
      sprintf(fileName, "%s/%s/%s.%06i", io->ensightDir, io->miscNames[i],
              io->miscNames[i], iFile);
      writeEnsightData(fileName,io->miscNames[i],flame->n2,flame->n1,1,0,
                       io->miscIn[i]);
    }

    // destroy the flamelet memory
    destroyFlame(flame);
    destroyFlow(flow);
    destroyIO(io);
  }
}
/* ------------------------------------------------------------------------- */

void writeBinStr(FILE *pFile, int len, char *str)
{
  char *buf = (char*) calloc(len, sizeof(char));
  sprintf(buf, "%s", str);
  fwrite(buf, sizeof(char), len, pFile);
  free(buf);
}
void writeBinInt(FILE *pFile, int len, int *var)
{
  fwrite(var, sizeof(int), len, pFile);
}
void writeBinReal(FILE *pFile, int len, float *var)
{
  fwrite(var, sizeof(float), len, pFile);
}
/* ------------------------------------------------------------------------- */

void convertEnsightParticle(optStruct *opts)
{
  int i, j, n, geomFile = 0, nData, *ids;
  float var, *xcoord, *ycoord, *zcoord;
  FILE *pFile, *binFile;
  char binFileName[200], buffer[80];

  for (n=0; n<opts->nFiles; n++) {
    /* open the regular file */
    pFile = fopen(opts->dataFile[n], "r");
    
    /* create an alternate data name */
    strcpy(binFileName, opts->dataFile[n]);
    strcat(binFileName, ".bin");
    binFile = fopen(binFileName, "wb");
    
    /* determine whether it is a geometry file or variable file */
    fgets(buffer, 80, pFile);
    fgets(buffer, 80, pFile);
    geomFile = 0;
    if (!strcmp(buffer, "particle coordinates\n"))
      geomFile = 1;
    rewind(pFile);
    
    fprintf(stderr, "converting %s --> %s\n", opts->dataFile[n], binFileName);
    
    if (geomFile) {
      /* write the description and header */
      strcpy(buffer, "C Binary");
      fwrite("C Binary", sizeof(char), 80, binFile);
      fgets(buffer, 80, pFile);
      fwrite(buffer, sizeof(char), 80, binFile);
      fgets(buffer, 80, pFile);
      fwrite(buffer, sizeof(char), 80, binFile);
      fscanf(pFile,"%10d\n",&nData);
      fwrite(&nData, sizeof(int), 1, binFile);
      
      /* allocate memory for geometry data */
      ids    = (int*)   calloc(nData, sizeof(int));
      xcoord = (float*) calloc(nData, sizeof(float));
      ycoord = (float*) calloc(nData, sizeof(float));
      zcoord = (float*) calloc(nData, sizeof(float));
      
      /* read in the geometry data */
      for (i=0; i<nData; i++)
        fscanf(pFile,"%i%f%f%f\n",&ids[i],&xcoord[i],&ycoord[i],&zcoord[i]);
      
      /* write the particle ids as a single list */
      fwrite(ids, sizeof(int), nData, binFile);
      
      /* write the coordinates as a list of tuples */
      for (i=0; i<nData; i++) {
        fwrite(&xcoord[i], sizeof(float), 1, binFile);
        fwrite(&ycoord[i], sizeof(float), 1, binFile);
        fwrite(&zcoord[i], sizeof(float), 1, binFile);
      }
    }
    else {
      fwrite("C Binary", sizeof(char), 80, binFile);
      fgets(buffer, 80, pFile);
      fwrite(buffer, sizeof(char), 80, binFile);
      do {
        fscanf(pFile, "%f", &var);
        fwrite(&var, sizeof(float), 1, binFile);
      } while(!feof(pFile));
    }
    
    fclose(pFile); fclose(binFile);
  }
}
/* ------------------------------------------------------------------------- */

void getZfromPDF(optStruct *opts)
{
  int n, i, j, k, row, *nn, n1, n2, pdfFlag, pdfBound[2], **pdfBound2D, iM, nM;
  char fileName[200];

  // allocate pointers to internal structs
  ioStruct    *io        = (ioStruct*)    calloc(1, sizeof(ioStruct));
  driveStruct *driver    = (driveStruct*) calloc(1, sizeof(driveStruct));
  timeStruct  *time      = (timeStruct*)  calloc(1, sizeof(timeStruct));
  flowStruct  *flow      = (flowStruct*)  calloc(1, sizeof(flowStruct));
  flameStruct *flame     = (flameStruct*) calloc(1, sizeof(flameStruct));

  // allocate the arrays
  double *t      = (double*) calloc(opts->nFiles+1, sizeof(double));
  double *Z1     = (double*) calloc(opts->nFiles+1, sizeof(double));
  double *Z2     = (double*) calloc(opts->nFiles+1, sizeof(double));
  double *sumPDF = (double*) calloc(opts->nFiles+1, sizeof(double));
  double **x, *x1, *x2, meanPDF[3], **Y, **pdf, *tMass, *mass, m;

  x  = (double**) calloc(2, sizeof(double*));
  nn = (int*)     calloc(2, sizeof(int));

  // read in the mass from a monitor file
  if (strcmp(opts->monitorFile, "none"))
    nM = readMonitorFile(opts->monitorFile, opts->iT, &tMass, opts->iV, &mass);

  fprintf(stderr, "Getting mean mixture fraction from PDF\n");
  for (n=0; n<opts->nFiles; n++) {
    /* ----- set default values ----- */
    setDefaults(io, driver, time, flame, flow, NULL);
    
    // read in the flow file
    readFlamelet(opts->dataFile[n], io, driver, time, flame, flow);
    readFlowData(opts->dataFile[n], driver, flow);
    
    // store the time
    t[n+1] = time->t;
    fprintf(stderr, "%s: %8.6f\n", opts->dataFile[n], t[n+1]);
    
    if (strcmp(opts->monitorFile, "none")) {
      // find the corresponding index in the mass array
      if (t[n+1] < tMass[0] || t[n+1] > tMass[nM-1] ) {
        fprintf(stderr,"No mass data for t = %8.6f, continuing...\n",t[n+1]);
        continue;
      }
      
      iM = 0;
      do { iM++; } while ( tMass[iM] < t[n+1] && iM < nM-1 );
      iM--;
      
      if (iM < nM-1)
        m = linInterp(t[n+1], tMass[iM], mass[iM], tMass[iM+1], mass[iM+1]);
      else {
        if (tMass[iM] == t[n+1])
          m = mass[iM];
        else
          continue;
      }
    }
    else
      m = 1.0;
    
    // integrate the pdf for this time
    if (flame->nDim == 1) {
      pdfBound[0] = 0;
      pdfBound[1] = flame->n1-1;
      pdfFlag = 1;
      
      // get the pdf normalisation
      sumPDF[n+1] = 0.0;
      for (j=0; j<flame->n1; j++) sumPDF[n+1] += flow->pdf[j];
      
      Z1[n+1] = m*calcLinMeanY(0, 1, -1, -1, flame->n1, flame->x1, flame->x1,
                               flow->pdf, pdfBound, pdfFlag);
      Z2[n+1] = 0.0;
    }
    else if (flame->nDim == 2) {
      x[0]  = x1 = flame->x1;
      x[1]  = x2 = flame->x2;
      nn[0] = n1 = flame->n1;
      nn[1] = n2 = flame->n2;
      if (!strcmp(io->transZ2, "doran") || !strcmp(io->transZ2, "standard"))
        driver->TRANS_Z2 = 1;
      else if (!strcmp(io->transZ2, "hasse"))
        driver->TRANS_Z2 = 2;
      
      // allocate and fill memory to compute means
      // for computing the mean, compute the function Y
      Y          = (double**) calloc(3, sizeof(double*));
      pdf        = (double**) calloc(n2, sizeof(double*));
      pdfBound2D = (int**)    calloc(n2, sizeof(int*));
      for (i=0; i<3; i++) Y[i] = (double*) calloc(n2*n1, sizeof(double));
      
      sumPDF[n+1] = 0.0;
      for (k=0; k<n2; k++) {
        // compute over entire bounds
        pdf[k]        = (double*) calloc(n1, sizeof(double));
        pdfBound2D[k] = (int*) calloc(2, sizeof(int));
        pdfBound2D[k][0] = 0;
        pdfBound2D[k][1] = n1-1;
        
        row = k*n1;
        for (j=0; j<n1; j++) {
          if (driver->TRANS_Z2 == 1) {
            Y[0][row+j] = x1[j] + x2[k]*(1.0 - x1[j]);
            Y[1][row+j] = x1[j];
            Y[2][row+j] = x2[k]*(1.0 - x1[j]);
          }
          else if (driver->TRANS_Z2 == 2) {
            Y[0][row+j] = x1[j];
            Y[1][row+j] = x1[j]*(1.0 - x2[k]);
            Y[2][row+j] = x2[k]*x1[j];
          }
          
          sumPDF[n+1] += flow->pdf[row+j];
        } // end loop over columns
      } // end loop over rows
      
      // copy the pdf to a structure array
      convertArray1D_MD(n2, n1, 1, flow->pdf, &pdf);
      
      // compute the total mean and the mean of each direction
      for (i=0; i<3; i++)
        meanPDF[i] = calcLinMeanY_2D(0, 1, nn, x, Y[i], -1, -1, -1,
                                     pdf, pdfBound2D, -1, driver->TRANS_Z2);
      
      Z1[n+1] = m*meanPDF[1];
      Z2[n+1] = m*meanPDF[2];
      
      for (k=0; k<n2; k++) {
        free(pdfBound2D[k]); free(pdf[k]);
      }
      free(pdf); free(pdfBound2D);
      for (i=0; i<3; i++) free(Y[i]);
      free(Y);
    }
    
    destroyFlame(flame);
    destroyFlow(flow);
    destroyIO(io);
  }

  // write the results to a file
  FILE *pFile = fopen("Zmean_integrate.dat", "w");

  fprintf(pFile, "# time Z1 Z2 sumPDF\n");
  for (n=1; n<opts->nFiles; n++)
    fprintf(pFile, "%15.6e %15.6e %15.6e %15.6e\n",t[n],Z1[n],Z2[n],sumPDF[n]);

  fclose(pFile);

  free(t); free(Z1); free(Z2); free(x); free(nn);
  if (strcmp(opts->monitorFile, "none")) {
    free(mass); free(tMass);
  }
}
/* ------------------------------------------------------------------------- */

void checkGrid(optStruct *opts)
{
  int i;
  double *x, *dx;
  FILE *pFile;

  /* initGrid(n, 0.0, 1.0, xP, pFine, &x, &dx); */
  initGrid(opts->nGrid, 0.0, 1.0, opts->nPivot, opts->pivot, opts->gridPerc,
           &x, &dx);

  if (!strcmp(opts->outFile, "none")) { /* write evaluation grid */
    pFile = fopen("testGrid.dat", "w");
    
    fprintf(pFile, "# x  dx  r\n");
    fprintf(pFile, "%8.6f %8.6f %8.6f\n", x[0], dx[0], 0.0);
    for (i=1; i<opts->nGrid-1; i++)
      fprintf(pFile, "%8.6f %8.6f %8.6f\n", x[i], dx[i], dx[i]/dx[i-1]);
    fprintf(pFile, "%8.6f %8.6f %8.6f\n", x[opts->nGrid-1], 0.0 , 0.0);
  }
  else { /* write in NGA format */
    pFile = fopen(opts->outFile, "w");
    
    fprintf(pFile, "%i\n", opts->nGrid);
    for (i=0; i<opts->nGrid; i++)
      fprintf(pFile, "%8.6f\n", x[i]);
  }
  
  fclose(pFile);
  free(x);
  free(dx);
}
/* ------------------------------------------------------------------------- */

void writeEnsightConst(optStruct *opts, int n, int freq)
{
  int i, j, nt, L1 = 3, LN = 8, checkTime = 0, nEns, tNext;
  double *t, *var, *timeVals = NULL;
  char fileOut[200];
  FILE *pFile;

  /* if a current ensight case file is specified to match the frequency,
     read in the time array */
  if (strcmp(opts->ensCase, "none")) {
    checkTime = 1;
    nEns = readEnsightCaseTime(opts->ensCase, &timeVals);
  }
  
  // read in the monitor file
  nt = readMonitorFile(opts->monitorFile, opts->iT, &t, opts->iV, &var);

  if (!strcmp(opts->outFile, "none")) {
    // generate the output file name
    strcpy(fileOut, opts->monitorFile);
    strcat(fileOut, ".ens");
  }
  else
    strcpy(fileOut, opts->outFile);
  pFile = fopen(fileOut, "w");

  // write just the variable at the appropriate frequency
  tNext = 0;
  for (i=0, j=0; i<nt; i++) {
    if (checkTime) {
      if (tNext >= nEns) {
        fprintf(stderr, "No more timesteps. Exiting\n");
        break;
      }
      if (fabs(t[i] - timeVals[tNext]) < 1.0e-7) {
        tNext++;
        fprintf(pFile, "%12.5e\n", var[i]);
      }
      else if (t[i] > timeVals[tNext]) {
        fprintf(stderr, "Error: file time %15.6e > ensight time %15.6e\n",
                t[i], timeVals[tNext]);
        exit(2);
      }
      else {
        //fprintf(stderr, "\t--> skipping %15.6e ...\n", t[i]);
        continue;
      }
    }
    else {
      // write just the variable at the appropriate frequency
      if ((i+1)%freq == 0) {
        fprintf(pFile, "%12.5e\n", var[i]);
      }
    }
  }

  fclose(pFile);
  free(t);
  free(var);
}
/* ------------------------------------------------------------------------- */

int  readEnsightCaseTime(char *fileName, double **time)
{
  int nt;
  char  tag[50], buffer[200], *ident;
  FILE *pFile = fopen(fileName, "r");

  // find the beginning of the TIME section
  while (!feof(pFile)) {
    fscanf(pFile, "%s", tag); // read in the first string
    if (!strcmp(tag, "TIME"))
      break;
  }
  
  // read in the time info
  fgets(buffer, 100, pFile);    // time set line
  fgets(buffer, 100, pFile);    // time set
  fgets(buffer, 100, pFile);    // number of steps
  ident = strtok(buffer, ":");
  nt = atoi(strtok(NULL,""));   // store number of steps
  fgets(buffer, 100, pFile);    // file start number
  fgets(buffer, 100, pFile);    // file increment
  fgets(buffer, 12, pFile);     // "time values:"

  // read in already existing time steps
  readRealArray(pFile, nt, time);  
      
  fclose(pFile);
  return nt;
}
/* ------------------------------------------------------------------------- */

void cleanMonitorFile(optStruct *opts)
{
  int i, j;
  char header[2][200], fileNew[200];
  FILE *pFile;
  int    nt, ntNew, jL, jR;
  double *time, *var, *t_clean, *v_clean;

  for (i=0; i<opts->nFiles; i++) {
    /* get the header rows for storing */
    pFile = fopen(opts->dataFile[i], "r");
    fgets(header[0], 200, pFile);
    fgets(header[1], 200, pFile);
    fclose(pFile);
    
    /* read in the data from the monitor file */
    nt = readMonitorFile(opts->dataFile[i], opts->iT, &time, opts->iV, &var);
      
    /* run through the file to see if it needs to be cleaned */
    jR = -1;
    for (j=nt-2; j>=0; j--)
      if (time[j] > time[j+1]) {
        jR = j+1;
        break;
      }
    
    if (jR > 0) { /* need to update the file */
      t_clean = (double*) calloc(nt, sizeof(double));
      v_clean = (double*) calloc(nt, sizeof(double));
      
      /* always keep from the end */
      t_clean[nt-1] = time[nt-1];
      v_clean[nt-1] = var[nt-1];
      jR = nt-1;
      ntNew = 0;
      for (j=nt-2; j>=0; j--) {
        if (time[j] < t_clean[jR]) {
          t_clean[jR-1] = time[j];
          v_clean[jR-1] = var[j];
          jR--;
          ntNew++;
        }
      }

      /* write the new file */
      strcpy(fileNew, opts->dataFile[i]);
      strcat(fileNew, ".clean");
      
      pFile = fopen(fileNew, "w");
      fprintf(pFile, "%s", header[0]);
      fprintf(pFile, "%s", header[1]);

      for (j=jR; j<nt; j++)
        fprintf(pFile, "%.9f %.11f\n", t_clean[j], v_clean[j]);
      
      fclose(pFile);
      free(t_clean); free(v_clean);
    } 
    
    free(time); free(var);
  }
}
/* ------------------------------------------------------------------------- */

void combineMonitorFiles(optStruct *opts)
{
  int i, j;
  char header[2][200], fileOut[200], *head, buffer[1000], id[50], **names;
  FILE *pFile;
  int    *nt, ntNew, jL, jR;
  double **time, **var;

  nt    = (int*)     calloc(opts->nFiles, sizeof(int));
  names = (char**)   calloc(opts->nFiles, sizeof(char*));
  time  = (double**) calloc(opts->nFiles, sizeof(double*));
  var   = (double**) calloc(opts->nFiles, sizeof(double*));

  /* get the case id */
  strcpy(buffer, opts->dataFile[0]);
  head = strtok(buffer, ".");
  head = strtok(NULL, ".\n");
  strcpy(id, head);
  fprintf(stderr, "combining files for case %s\n", id);

  /* read in all the monitor files */
  for (i=0; i<opts->nFiles; i++) {
    names[i] = (char*) calloc(50, sizeof(char));
    
    /* get the names from the file */
    strcpy(buffer, opts->dataFile[i]);
    head = strtok(buffer, ".");
    strcpy(names[i], head);
    
    /* read in the data from the monitor file */
    nt[i]=readMonitorFile(opts->dataFile[i],opts->iT,&time[i],opts->iV,&var[i]);
  }

  for (i=1; i<opts->nFiles; i++) {
    if (nt[i] != nt[0]) {
      fprintf(stderr, "number of time steps in file don't match\n");
      exit(2);
    }
    for (j=0; j<nt[i]; j++) {
      if (fabs(time[i][j]-time[0][j]) > 1.0e-6) {
        fprintf(stderr, "time in all files does not match\n");
        exit(2);
      }
    }
  }

  /* write the new file */
  if (!strcmp(opts->outFile, "none"))
    sprintf(fileOut, "monitor.%s", id);
  else
    sprintf(fileOut, "%s.%s", opts->outFile, id);
  pFile = fopen(fileOut, "w");
  
  // write the header 
  fprintf(pFile, "# time[s] ");
  for (i=0; i<opts->nFiles; i++) fprintf(pFile, "%s ", names[i]);
  fprintf(pFile, "\n");

  for (j=0; j<nt[0]; j++) {
    fprintf(pFile, "%15.6e ", time[0][j]);
    for (i=0; i<opts->nFiles; i++)
      fprintf(pFile, "%15.6e ", var[i][j]);
    fprintf(pFile, "\n");
  }

  fclose(pFile);
}
/* ------------------------------------------------------------------------- */

int readMonitorFile(char *fileName, int iT, double **time, int iV, double **var)
{
  int   i, j, len, nHeader, firstChar, nCol;
  double tmp;
  
  // strings to store the tag for each line and comments (if any)
  char  tag[50], comment[1000], *col;
  
  FILE  *pFile = fopen(fileName, "r");
  fpos_t pos;

  if (!pFile)
    mducErr("monitor file not found", "readMonitorFile");

  len = nHeader = 0;
  // get the number of data lines
  while (!feof(pFile)) {
    strcpy(comment, "\n");
    fgets(comment, 1000, pFile);
    // assume header if doesn't start with digit ( or '-' for negative number)
    firstChar = strspn (comment," ");
    if (!isdigit(comment[firstChar]) && comment[firstChar] != '-' && 
        comment[firstChar] != '\n')
      nHeader++;
    else if (strcmp(comment, "\n") && comment[0] != '\"')
      len++;
  }
  /* len -= 2; */
  rewind(pFile);

  // allocate the memory
  time[0] = (double*) calloc(len, sizeof(double));
  var[0]  = (double*) calloc(len, sizeof(double));

  // get the two header lines
  for (i=0; i<nHeader; i++)
    fgets(comment, 1000, pFile);

  // get the number of columns
  fgetpos(pFile, &pos);  // get the of the buffer at the line start
  fgets(comment, 1000, pFile);
  fsetpos(pFile, &pos);
  nCol = 0;
  col = strtok(comment, " ");
  while (col != NULL) {
    nCol++;
    col = strtok(NULL, " ");
  }
  nCol--;
  fprintf(stderr, "%i columns in file\n", nCol);

  // read in the variables from the file
  for (i=0; i<len; i++) {
    fgetpos(pFile, &pos);  // get the of the buffer at the line start
    // read in variables until time index
    for (j=0; j<iT-1; j++) tmp = readReal(pFile);
    time[0][i] = readReal(pFile);
    fgets(comment, 1000, pFile); // read in rest of line
    
    // reset file to line start
    fsetpos(pFile, &pos);
    for (j=0; j<iV-1; j++) tmp = readReal(pFile);
    var[0][i]  = readReal(pFile);
    if (iV < nCol) fgets(comment, 1000, pFile); // read in rest of line
  }
  
  fclose(pFile);

  return len;
}
/* ------------------------------------------------------------------------- */

void checkChi2D(optStruct *opts, void *mem)
{
  int i, j, k, r;

  // declare flamelet structs so that internal functions can be used
  driveStruct *driver;
  ioStruct    *io;
  chemStruct  *chem;
  timeStruct  *time;
  flameStruct *flame;
  flowStruct  *flow;

  getMDUCStructs(&io, &driver, &chem, &time, &flow, &flame, mem);

  int     *n = (int*)     calloc(2, sizeof(int));  
  double **x = (double**) calloc(2, sizeof(double*));
  n[0] = flame->n1;
  n[1] = flame->n2;
  x[0] = flame->x1;
  x[1] = flame->x2;

  double **chi1  = (double**) calloc(flame->n2, sizeof(double*));
  double **chi12 = (double**) calloc(flame->n2, sizeof(double*));
  double **chi2  = (double**) calloc(flame->n2, sizeof(double*));
  for (k=0; k<flame->n2; k++) {
    chi1[k]  = (double*) calloc(flame->n1, sizeof(double));
    chi12[k] = (double*) calloc(flame->n1, sizeof(double));
    chi2[k]  = (double*) calloc(flame->n1, sizeof(double));
  }

  // compute the 2D scalar dissipation rate profile
  calcChi2D(io->chiMode, n, x, flow->ZRef, flow->chiRef, chi1, chi12, chi2,
            driver->TRANS_Z2, driver->CHI_SCALE);
  
  // convert to unstructure format for writing
  convertArrayMD_1D(n[0], n[1], 1, &chi1,  flow->chi1);
  convertArrayMD_1D(n[0], n[1], 1, &chi12, flow->chi12);
  convertArrayMD_1D(n[0], n[1], 1, &chi2,  flow->chi2);

  // write output
  writeFlowData(io->flowData, io, driver, time, flame, flow);

  for (k=0; k<flame->n2; k++) {
    free(chi1[k]);
    free(chi12[k]);
    free(chi2[k]);
  }
  free(chi1);
  free(chi12);
  free(chi2);
}
/* ------------------------------------------------------------------------- */

void chiProfile(optStruct *opts)
{
  int j;
  double *x, *dx;
  double *chiErfc, *chiLog;

  double chiRef = 0.125;
  double Zref = 0.1;

  chiErfc = (double*) calloc(opts->nGrid, sizeof(double));
  chiLog  = (double*) calloc(opts->nGrid, sizeof(double));

  initGrid(opts->nGrid, 0.0, 1.0, opts->nPivot, opts->pivot, opts->gridPerc,
           &x, &dx);

  for (j=0; j<opts->nGrid; j++) {
    chiErfc[j] = chiRef*calcErfcForm(x[j], Zref, 1.0);
    chiLog[j]  = chiRef*calcLogForm(x[j], Zref, 1.0);
  }

  FILE *pFile = fopen("chiErfc.dat", "w");

  fprintf(pFile, "# Z chi\n");
  for (j=0; j<opts->nGrid; j++)
    fprintf(pFile, "%8.6f %8.6f\n", x[j], chiErfc[j]);
  fclose(pFile);

  pFile = fopen("chiLog.dat", "w");
  fprintf(pFile, "# Z chi\n");
  for (j=0; j<opts->nGrid; j++)
    fprintf(pFile, "%8.6f %8.6f\n", x[j], chiLog[j]);
  fclose(pFile);

  free(x); free(dx);
  free(chiErfc); free(chiLog);
}
/* ------------------------------------------------------------------------- */

void postProcessChi(optStruct *opts)
{
  int    i, j, k, r, nt, jM, kM, n1, n2, jj, kk;
  char   fileName[200], prefix[200], buffer[200], transform[20], *head;
  double Z1, Z2, Zmean[3], *x1, *x2, tmp;
  double *Z1ref, *Z2ref;
  double den0, den1, den2;
  double chiRef[4], chiMean[4], *chi0;
  double **chi1_m, **chi2_m;
  double **chi1_F, **chi12_F, **chi2_F, **chi0_F, **theta_F;
  double ***chi1_m_F, ***chi2_m_F;
  double **chiRef_F, **chi1_r_F, **chi2_r_F;
  double *pdf_Z1, *pdf_Z2;
  double *theta, thetaMean, thetaRef;
  double **JSD, JSD_log[4], JSD_erfc[4];
  FILE   *pFile, *meanFile;

  double (**chiForm)(double, double, double);

  // allocate pointers to internal structs
  mducStruct  *mem    = (mducStruct*)  calloc(1, sizeof(mducStruct));
  ioStruct    *io     = (ioStruct*)    calloc(1, sizeof(ioStruct));
  driveStruct *driver = (driveStruct*) calloc(1, sizeof(driveStruct));
  timeStruct  *time   = (timeStruct*)  calloc(1, sizeof(timeStruct));
  flowStruct  *flow   = (flowStruct*)  calloc(1, sizeof(flowStruct));
  flameStruct *flame  = (flameStruct*) calloc(1, sizeof(flameStruct));

  /* assign the methods */
  chiForm = (double (**)(double, double, double))
    calloc(opts->nChiMethod, sizeof(double (*)(double, double, double)));
  for (i=0; i<opts->nChiMethod; i++) {
    if (!strcmp(opts->chiMethod[i], "log"))
      chiForm[i] = &calcLogForm;
    else if (!strcmp(opts->chiMethod[i], "erfc"))
      chiForm[i] = &calcErfcForm;
    else if (!strcmp(opts->chiMethod[i], "quad"))
      chiForm[i] = &calcQuadForm;
    else
      fprintf(stderr, "%s form not a valid option\n", opts->chiMethod[i]);
  }

  /* write header for mean data */
  meanFile = fopen("chiMean.dat", "w");
  fprintf(meanFile, "# time ");
  fprintf(meanFile, "chi1-mean chi12-mean chi2-mean chi0-mean theta-mean ");
  fprintf(meanFile, "chi1-ref chi12-ref chi2-ref chi0-ref theta-ref ");
  for (i=0; i<opts->nChiMethod; i++)
    fprintf(meanFile, "chi1-%s-JSD chi12-%s-JSD chi2-%s-JSD ", 
            opts->chiMethod[i], opts->chiMethod[i], opts->chiMethod[i]);
  fprintf(meanFile, "\n");

  /* set up pointer for methods */
  chi0_F  = (double**) calloc(opts->nChiMethod, sizeof(double*));
  chi1_F  = (double**) calloc(opts->nChiMethod, sizeof(double*));
  chi12_F = (double**) calloc(opts->nChiMethod, sizeof(double*));
  chi2_F  = (double**) calloc(opts->nChiMethod, sizeof(double*));
  theta_F = (double**) calloc(opts->nChiMethod, sizeof(double*));

  /* marginal chi distributions (_F means presumed) */
  chi1_m   = (double**)  calloc(3, sizeof(double*));
  chi2_m   = (double**)  calloc(3, sizeof(double*));
  chi1_m_F = (double***) calloc(opts->nChiMethod, sizeof(double**));
  chi2_m_F = (double***) calloc(opts->nChiMethod, sizeof(double**));
  chiRef_F = (double**)  calloc(opts->nChiMethod, sizeof(double*));
  chi1_r_F = (double**)  calloc(opts->nChiMethod, sizeof(double*));
  chi2_r_F = (double**)  calloc(opts->nChiMethod, sizeof(double*));
  JSD      = (double**)  calloc(opts->nChiMethod, sizeof(double*));

  strcpy(transform, "none");
  fprintf(stderr, "Post-processing scalar dissipation DNS data:\n");
  for (nt=0; nt<opts->nFiles; nt++) {
    /* ----- set default values ----- */
    setDefaults(io, driver, time, flame, flow, NULL);

    // read in the flow file
    driver->PDF_MODE = driver->CHI_MODE == -1;
    readFlamelet(opts->dataFile[nt], io, driver, time, flame, flow);
    readFlowData(opts->dataFile[nt], driver, flow);
    fprintf(stderr,"\t--> %s\n",opts->dataFile[nt]);

    // assign local grid parameters for convenience
    n1 = flame->n1; n2 = flame->n2;
    x1 = flame->x1; x2 = flame->x2;

    /* assign memory for additional data from the dns*/
    chi0      = (double*) calloc(n1*n2, sizeof(double));
    theta     = (double*) calloc(n1*n2, sizeof(double));
    pdf_Z1    = (double*) calloc(n1, sizeof(double));
    pdf_Z2    = (double*) calloc(n2, sizeof(double));
    Z1ref     = (double*) calloc(n2, sizeof(double));
    Z2ref     = (double*) calloc(n1, sizeof(double));
    chi1_m[1] = (double*) calloc(n1, sizeof(double));
    chi1_m[2] = (double*) calloc(n2, sizeof(double));
    chi2_m[1] = (double*) calloc(n1, sizeof(double));
    chi2_m[2] = (double*) calloc(n2, sizeof(double));

    // assign local memory for presumed chi shapes
    for (i=0; i<opts->nChiMethod; i++) {
      chi0_F[i]   = (double*)  calloc(n1*n2, sizeof(double));
      chi1_F[i]   = (double*)  calloc(n1*n2, sizeof(double));
      chi12_F[i]  = (double*)  calloc(n1*n2, sizeof(double));
      chi2_F[i]   = (double*)  calloc(n1*n2, sizeof(double));
      theta_F[i]  = (double*)  calloc(n1*n2, sizeof(double));
      chi1_m_F[i] = (double**) calloc(3, sizeof(double*));
      chi2_m_F[i] = (double**) calloc(3, sizeof(double*));
      chi1_m_F[i][1] = (double*) calloc(n1, sizeof(double));
      chi1_m_F[i][2] = (double*) calloc(n2, sizeof(double));
      chi2_m_F[i][1] = (double*) calloc(n1, sizeof(double));
      chi2_m_F[i][2] = (double*) calloc(n2, sizeof(double));
      chiRef_F[i] = (double*) calloc(3, sizeof(double));
      chi1_r_F[i] = (double*) calloc(n2, sizeof(double));
      chi2_r_F[i] = (double*) calloc(n1, sizeof(double));
      JSD[i]      = (double*) calloc(4, sizeof(double));
    }
      
    /* ----- compute mean quantities  ----- */
    chiMean[0] = chiMean[1] = chiMean[2] = chiMean[3] = thetaMean = 0.0;
    Zmean[0] = Zmean[1] = Zmean[2] = 0.0;
    for (k=0; k<n2; k++) {
      r = k*n1;
      for (j=0; j<n1; j++) {
        if (x1[j] + x2[k] > 1.0)
          flow->chi1[r+j] = flow->chi2[r+j] = flow->chi12[r+j] = 0.0;
        else {
          chi0[r+j]=flow->chi1[r+j]+2*flow->chi12[r+j]+flow->chi2[r+j];
          chiMean[0] +=       chi0[r+j]*flow->pdf[r+j];
          chiMean[1] += flow->chi1[r+j]*flow->pdf[r+j];
          chiMean[2] += flow->chi2[r+j]*flow->pdf[r+j];
          chiMean[3] += flow->chi12[r+j]*flow->pdf[r+j];
          if (flow->chi1[r+j]*flow->chi2[r+j] > 0.0)
            theta[r+j] = acos(flow->chi12[r+j]/
                              sqrt(flow->chi1[r+j]*flow->chi2[r+j]));
          thetaMean += theta[r+j]*flow->pdf[r+j];

          // compute the marginals
          pdf_Z1[j] += flow->pdf[r+j];
          pdf_Z2[k] += flow->pdf[r+j];
        }
        Zmean[0] += (x1[j]+x2[k])*flow->pdf[r+j];
        Zmean[1] += x1[j]*flow->pdf[r+j];
        Zmean[2] += x2[k]*flow->pdf[r+j];
      } /* end loop over columns */
    } /* end loop over rows */

      /* determine which bin Zmean falls into */
    jM = -1;
    do { jM++; } while (0.5*(x1[jM]+x1[jM+1]) < Zmean[1] && jM < n1-1);
    kM = -1;
    do { kM++; } while (0.5*(x2[kM]+x2[kM+1]) < Zmean[2] && kM < n2-1);

    /* ----- compute the conditional means to use as a reference ----- */
    for (k=0; k<n2; k++) {
      r = k*n1;
      for (j=0; j<n1; j++) {
        Z1ref[k] += x1[j]*flow->pdf[r+j]/pdf_Z2[k];
        Z2ref[j] += x2[k]*flow->pdf[r+j]/pdf_Z1[j];
        chi1_m[1][j] += flow->chi1[r+j]*flow->pdf[r+j]/pdf_Z1[j];
        chi1_m[2][k] += flow->chi1[r+j]*flow->pdf[r+j]/pdf_Z2[k];
        chi2_m[1][j] += flow->chi2[r+j]*flow->pdf[r+j]/pdf_Z1[j];
        chi2_m[2][k] += flow->chi2[r+j]*flow->pdf[r+j]/pdf_Z2[k];
      } /* end loop over columns */
    } /* end loop over rows */

      /* ----- compute a reference chi as the mean chi at Z_j mean ----- */
    chiRef[1] = tmp = 0.0;
    for (k=0; k<n2; k++) {
      chiRef[1] += flow->chi1[k*n1+jM]*flow->pdf[k*n1+jM];
      tmp += flow->pdf[k*n1+jM];
    }
    chiRef[1] /= tmp;
    fprintf(stderr, "chiRef[1]_old = %15.6e\n", chiRef[1]);

    chiRef[2] = tmp = 0.0;
    for (j=0; j<n2; j++) {
      chiRef[2] += flow->chi2[kM*n1+j]*flow->pdf[kM*n1+j];
      tmp += flow->pdf[kM*n1+j];
    }
    chiRef[2] /= tmp;

    chiRef[0] = tmp = 0.0;
    for (k=0; k<n2; k++) {
      r = k*n1;
      if (1 - x2[k] > Zmean[0]) {
        jM = n1;
        do {
          jM--;
        } while (1 - 0.5*(x1[jM]+x1[jM-1]) - x2[k] < Zmean[0] && jM > 0);
          
        chiRef[0] += chi0[r+jM]*flow->pdf[r+jM];
        tmp += flow->pdf[r+jM];
      }
    }
    chiRef[0] /= tmp;

    // get the angle at the reference value
    thetaRef = theta[kM*n1+jM];

    /* ----- write mean/error data to file ----- */
    fprintf(meanFile, "%15.6e ", time->t);
    fprintf(meanFile, "%15.6e %15.6e %15.6e %15.6e %15.6e ",
            chiMean[1], chiMean[3], chiMean[2], chiMean[0], thetaMean);
    fprintf(meanFile, "%15.6e %15.6e %15.6e %15.6e %15.6e ",
            chiRef[1], chiRef[3], chiRef[2], chiRef[0], thetaRef);
      
    /* ----- start compute the presumed distributions ----- */
    for (i=0; i<opts->nChiMethod; i++) {
      /* integrate the presumed form with Zmean as reference to mean chi */
      den0 = den1 = den2 = 0.0;
      for (k=0; k<n2; k++) {
        r = k*n1;
        for (j=0; j<n1; j++) {
          den0 += chiForm[i](x1[j]+x2[k],Zmean[1]+Zmean[2],1)*flow->pdf[r+j];
          den1 += chiForm[i](x1[j],Zmean[1],1)*flow->pdf[r+j];
          den2 += chiForm[i](x2[k],Zmean[2],1)*flow->pdf[r+j];
        } /* end loop over columns */
      } /* end loop over rows */

      /* compute the reference chi for the mean presumed distribution */
      chiRef_F[i][0] = chiMean[0]/den0;
      chiRef_F[i][1] = chiMean[1]/den1;
      chiRef_F[i][2] = chiMean[2]/den2;

      /* compute the mean presumed shape */
      for (j=0; j<n1; j++)
        chi1_m_F[i][1][j] = chiRef_F[i][1]*chiForm[i](x1[j],Zmean[1],1);

      for (k=0; k<n2; k++)
        chi2_m_F[i][2][k] = chiRef_F[i][2]*chiForm[i](x2[k],Zmean[2],1);

      /* compute the distribution of each chi in the other direction */
      for (k=0, kk=n2-1; k<n2; k++, kk--) {
        r = k*n1;
        for (j=0, jj=n1-1; j<n1; j++, jj--) {
          chi1_m_F[i][2][k] += chi1_m_F[i][1][j]*flow->pdf[r+j]/pdf_Z2[k];
          chi2_m_F[i][1][j] += chi1_m_F[i][2][k]*flow->pdf[r+j]/pdf_Z1[j];
        } /* end loop over columns */
      } /* end loop over rows */

      /* compute the reference values as a function of the other mixture */
      for (k=0; k<n2; k++) {
        r = k*n1;
        den1 = 0.0;
        for (j=0; j<n1; j++)
          den1+=chiForm[i](x1[j],Z1ref[k],1)*flow->pdf[r+j]/pdf_Z2[k];
        chi1_r_F[i][k] = chi1_m_F[i][2][k]/den1;
        chi1_r_F[i][k] = chi1_m[2][k]/den1;
      } /* end loop over rows */
          
      for (j=0; j<n1; j++) {
        den2 = 0.0;
        for (k=0; k<n2; k++) {
          r = k*n1;
          den2 += chiForm[i](x2[j],Z2ref[j],1)*flow->pdf[r+j]/pdf_Z1[j];
        }
        chi2_r_F[i][j] = chi2_m_F[i][1][j]/den2;
        chi2_r_F[i][j] = chi2_m[1][j]/den2;
      } /* end loop over columns */

      /* compute the presumed forms */
      for (k=0; k<n2; k++) {
        r = k*n1;
        for (j=0; j<n1; j++) {
          Z1 = x1[j];
          Z2 = x2[k];

          if (Z1 + Z2 <= 1.0) {
            chi0_F[i][r+j]  = chiRef_F[i][1]*chiForm[i](Z1+Z2,1-Zmean[1]-Zmean[2],1);
            chi1_F[i][r+j]  = chiRef_F[i][1]*chiForm[i](Z1     ,Zmean[1],1);
            /* chi1_F[i][r+j]  = chi1_r_F[i][k]*chiForm[i](Z1,Z1ref[k],1); */
            chi2_F[i][r+j]  = chiRef_F[i][2]*chiForm[i](Z2     ,Zmean[2],1);
            /* chi2_F[i][r+j]  = chi2_r_F[i][j]*chiForm[i](Z2     ,Zmean[2],1); */
            chi12_F[i][r+j] = 0.5*(chi0_F[i][r+j]
                                   - chi1_F[i][r+j] - chi2_F[i][r+j]);

            if (chi1_F[i][r+j]*chi2_F[i][r+j] > 0.0)
              theta_F[i][r+j] = acos(chi12_F[i][r+j]/
                                     sqrt(chi1_F[i][r+j]*chi2_F[i][r+j]));
          }
        } /* end loop over columns */
      } /* end loop over rows */

      /* compute the divergence to the actual distribution */
      JSD[i][1] = calcJensonShannon(n1*n2, 1, &flow->chi1,  &chi1_F[i]);
      JSD[i][2] = calcJensonShannon(n1*n2, 1, &flow->chi2,  &chi2_F[i]);
      JSD[i][3] = calcJensonShannon(n1*n2, 1, &flow->chi12, &chi12_F[i]);
      for (j=1; j<4; j++) fprintf(meanFile, "%15.6e ", JSD[i][j]);

    } /* end loop over methods */

      /* finish the line in the mean tracking file */
    fprintf(meanFile, "\n");
    fflush(meanFile);
      

      /*     for (k=0; k<n2; k++) */
/*         { */
/*           r = k*n1; */
/*           for (j=0; j<n1; j++) */
/*             { */
/*               chiMeanD1[j] += flow->chi1[r+j]*flow->pdf[r+j]/pdf_Z2[k]; */
/*             } */
/*         } */
      
 
   /*    den1 = 0.0; */
/*       for (k=0; k<n2; k++) */
/*         den1 += calcErfcForm(x2[k],Zmean[2],1)*marg2[k]; */
/*       chiRef[2] = chiMean[2]/den1; */
     

 /*      for (k=0; k<n2; k++) */
/*         { */
/*           r = k*n1; */
/*           /\* unscale the marginal *\/ */
/*           if (k == 0) */
/*             dZ = 0.5*(x2[k+1]-x2[k]); */
/*           else if (k == n2-1) */
/*             dZ = 0.5*(x2[k]-x2[k-1]); */
/*           else */
/*             dZ = 0.5*(x2[k+1] - x2[k-1]); */
              
/*           /\* if (marg2[k] > 1.0e-4) *\/ */
/* /\*             chiMarg1[k] = chiMean[1]*dZ/marg2[k]; *\/ */
/*           chiMarg1[k] = 0.0; */
/*           for (j=0; j<n1; j++) */
/*             { */
/*               if (x1[j]+x2[k] <= 1.0) */
/*                 chiMarg1[k] += chiMeanQ1[j]*flow->pdf[r+j]/pdf_Z2[k]; */
/*             } */
/*           /\* chiMarg1[k] /= chiMean[1]; *\/ */
/*             /\* chiMarg1[k] = chiRef[2]*calcErfcForm(x2[k],Zmean[2],1); *\/ */
/*         } */

/*       for (j=0; j<n1; j++) */
/* /\*         if (marg1[j] > 1.0e-2) *\/ */
/*           chiMarg2[j] = chiMean[2]*pdf_Z1[j]; */

      

/*       /\* compute the reference values *\/ */
/*       for (k=0; k<n2; k++) */
/*         { */
/*           r = k*n1; */
/*           den1 = den2 = 0.0; */
/*           /\* unscale the marginal *\/ */
/*           if (k == 0) */
/*             dZ = 0.5*(x2[k+1]-x2[k]); */
/*           else if (k == n2-1) */
/*             dZ = 0.5*(x2[k]-x2[k-1]); */
/*           else */
/*             dZ = 0.5*(x2[k+1] - x2[k-1]); */

/*           for (j=0; j<n1; j++) */
/*             { */
/*               Z1 = x1[j]; */
/*               Z2 = x2[k]; */

/*               /\* if (Z1 + Z2 <= 1.0) *\/ */
/* /\*                 den1 += calcErfcForm(Z1,Z1ref[k],1)*flow->pdf[r+j]/marg2[k]; *\/ */
/*                 den1 += calcQuadForm(Z1,Z1ref[k],1)*flow->pdf[r+j]/pdf_Z2[k]; */
/*             } */
/*           chiRef1[k] = chiMarg1e[k]/den1; */
/*         } */

/*       for (j=0; j<n1; j++) */
/*         { */
/*           Z1 = x1[j]; */
/*           den2 = 0.0; */
/*           for (k=0; k<n2; k++) */
/*             { */
/*               r = k*n1; */
/*               Z2 = x2[k]; */
/*               if (Z1 + Z2 <= 1.0) */
/*                 den2 += calcErfcForm(Z2,Z2ref[j],1)*flow->pdf[r+j]; */
/*             } */
/*           chiRef2[j] = chiMarg2[j]/den2; */
/*         } */

      /* write some info to a file */
/*       pFile = fopen("chiMarg.dat", "w"); */



/*       fclose(pFile); */

/*       // compute the presumed forms */
/*       for (i=0; i<4; i++) */
/*         JSD_log[i] = JSD_erfc[i] = 0.0; */
/*       for (i=0; i<opts->nChiMethod; i++) */
/*         { */
/*       for (k=0; k<n2; k++) */
/*         { */
/*           r = k*n1; */
/*           for (j=0; j<n1; j++) */
/*             { */
/*               Z1 = x1[j]; */
/*               Z2 = x2[k]; */

/*               if (Z1 + Z2 <= 1.0) */
/*                 { */

/*                   chi0_F[i][r+j]  = chiRef[0]*chiForm[i](1-Z1-Z2,Zmean[0],1); */
/*                   /\* chi1_F[i][r+j]  = chiRef[1]*chiForm[i](Z1     ,Zmean[1],1); *\/ */
/*                   chi1_F[i][r+j]  = chi1_r_F[i][k]*chiForm[i](Z1,Z1ref[k], 1); */
/*                   chi2_F[i][r+j]  = chiRef[2]*chiForm[i](Z2     ,Zmean[2],1); */
/*                   chi12_F[i][r+j] = 0.5*(chi0_F[i][r+j] */
/*                                          - chi1_F[i][r+j] - chi2_F[i][r+j]); */

/*                   if (chi1_F[i][r+j]*chi2_F[i][r+j] > 0.0) */
/*                     theta_F[i][r+j] = acos(chi12_F[i][r+j]/ */
/*                                            sqrt(chi1_F[i][r+j]*chi2_F[i][r+j])); */
                                         
                /*   // calculate log profile */
/*                   /\* if (flow->pdf[r+j] > 0.0) *\/ */
/*                     { */
/*                       chi0Log[r+j] = chiRef[0]*calcLogForm(1-Z1-Z2,Zmean[0],1); */
/*                       chi1Log[r+j] = chiRef1[k]*calcLogForm(Z1,Z1ref[k],1); */
/*                       chi2Log[r+j] = chiRef2[j]*calcLogForm(Z2,Z2ref[j],1); */
/*                     } */
/*                   chi12Log[r+j] = 0.5*(chi0Log[r+j]-chi1Log[r+j]-chi2Log[r+j]); */
/*                   if (chi1Log[r+j]*chi2Log[r+j] > 0.0) */
/*                     thetaLog[r+j] = acos(chi12Log[r+j]/ */
/*                                          sqrt(chi1Log[r+j]*chi2Log[r+j])); */
                  
/*                   // calculate inverfc profile */
/*                   /\* if (flow->pdf[r+j] > 0.0) *\/ */
/*                     { */
/*                       chi0Erfc[r+j] = chiRef[0]*calcErfcForm(1-Z1-Z2,Zmean[0],1); */
/*                       /\* chi1Erfc[r+j] = chiRef[1]*calcQuadForm(Z1,Zmean[1],1); *\/ */
/*                       /\* chi1Erfc[r+j] = chiRef1[k]*calcErfcForm(Z1,Z1ref[k],1); *\/ */
/*                       chi1Erfc[r+j] = chiRef1[k]*calcQuadForm(Z1,Z1ref[k],1); */
/*                       chi2Erfc[r+j] = chiRef2[j]*calcErfcForm(Z2,Z2ref[j],1); */
/*                     } */
/*                   chi12Erfc[r+j]=0.5*(chi0Erfc[r+j]-chi1Erfc[r+j]-chi2Erfc[r+j]); */
/*                   if (chi1Erfc[r+j]*chi2Erfc[r+j] > 0.0) */
/*                     thetaErfc[r+j] = acos(chi12Erfc[r+j]/ */
/*                                           sqrt(chi1Erfc[r+j]*chi2Erfc[r+j])); */
                  
/*                   if (!isnan(chi1Erfc[r+j])) */
/*                     chi1MeanErfc += chi1Erfc[r+j]*flow->pdf[r+j]; */
/*                 } */
/*             } /\* end loop over columns *\/ */
/*         } /\* end loop over rows *\/ */
/*         } /\* end loop over methods *\/ */

/*       fprintf(stderr, "chi1Mean: %15.6e/%15.6e (%8.6f)\n", chi1MeanErfc, */
/*               chiMean[1], fabs(chi1MeanErfc-chiMean[1])/chiMean[1]*100); */
      
    /* ----- write the output file ----- */
    // make file name
    strcpy(buffer, opts->dataFile[nt]);
    head = strtok(buffer, ".");
    strcpy(fileName, head);
    /*       sprintf(buffer, "_%05i.chi", nt); */
    /*       strcat(fileName, buffer); */
    strcat(fileName, ".chi");
    pFile = fopen(fileName, "w");

    // write the header
    fprintf(pFile, "header\n");
    fprintf(pFile, "title  = \"Scalar dissipation rate comparison\"\n");
    fprintf(pFile, "author = \"MDUC tools\"\n");
    fprintf(pFile, "date   = \"%s\"\n\n", io->date);

    fprintf(pFile, "time   = %13.6e [s]\n", time->t);
    fprintf(pFile, "Z1_min = %13.6e\n", flow->phi1min);
    fprintf(pFile, "Z1_max = %13.6e\n", flow->phi1max);
    if (flame->nDim == 2) {
      fprintf(pFile, "Z2_min = %13.6e\n", flow->phi2min);
      fprintf(pFile, "Z2_max = %13.6e\n", flow->phi2max); 
    }
    fprintf(pFile, "\n");

    fprintf(pFile, "chi1Mean  = %15.6e [1/s]\n", chiMean[1]);
    if (flame->nDim == 2) {
      fprintf(pFile, "chi12Mean = %15.6e [1/s]\n", chiMean[3]);
      fprintf(pFile, "chi2Mean  = %15.6e [1/s]\n", chiMean[2]);
    }
    fprintf(pFile, "Z1mean    = %15.6e\n", Zmean[1]);
    if (flame->nDim == 2)
      fprintf(pFile, "Z2mean    = %15.6e\n", Zmean[2]);
    fprintf(pFile, "\n");

    fprintf(pFile, "gridPoints     = %i\n", n1);
    if (flame->nDim == 2) {
      fprintf(pFile, "gridPoints2    = %i\n", n2);
      fprintf(pFile, "transformation = %s\n", transform);
    }
    fprintf(pFile, "\n");

    // write body
    fprintf(pFile, "body\n");
    if (!strcmp(transform, "hasse"))
      fprintf(pFile, "Z1+Z2\n");
    else
      fprintf(pFile, "Z1\n");
    writeArray(pFile, 1, n1, 1, 0, 5, x1);
    if (flame->nDim == 2) {
      if (!strcmp(transform, "doran"))
        fprintf(pFile, "Z2/(1-Z1)\n");
      else if (!strcmp(transform, "hasse"))
        fprintf(pFile, "Z2/(Z1+Z2)\n");
      else
        fprintf(pFile, "Z2\n");
      writeArray(pFile, 1, n2, 1, 0, 5, x2);
    }
    fprintf(pFile, "pdf\n");
    writeArray(pFile, n1, n2, 1, 0, 5, flow->pdf);
    fprintf(pFile, "chi1\n");
    writeArray(pFile, n1, n2, 1, 0, 5, flow->chi1);
    fprintf(pFile, "chi12\n");
    writeArray(pFile, n1, n2, 1, 0, 5, flow->chi12);
    fprintf(pFile, "chi2\n");
    writeArray(pFile, n1, n2, 1, 0, 5, flow->chi2);
    fprintf(pFile, "theta\n");
    writeArray(pFile, n1, n2, 1, 0, 5, theta);

    for (i=0; i<opts->nChiMethod; i++) {
      fprintf(pFile, "chi1-%s\n", opts->chiMethod[i]);
      writeArray(pFile, n1, n2, 1, 0, 5, chi1_F[i]);
      fprintf(pFile, "chi12-%s\n", opts->chiMethod[i]);
      writeArray(pFile, n1, n2, 1, 0, 5, chi12_F[i]);
      fprintf(pFile, "chi2-%s\n", opts->chiMethod[i]);
      writeArray(pFile, n1, n2, 1, 0, 5, chi2_F[i]);
      fprintf(pFile, "theta-%s\n", opts->chiMethod[i]);
      writeArray(pFile, n1, n2, 1, 0, 5, theta_F[i]);
    }

    fprintf(pFile, "trailer\n");
    fclose(pFile);
    /* ----- finish writing file ----- */

    /* ----- write the marginal file ----- */
    strcpy(buffer, opts->dataFile[nt]);
    head = strtok(buffer, ".");
    strcpy(fileName, head);
    /*       sprintf(buffer, "_%05i.chi_marg", nt); */
    /*       strcat(fileName, buffer); */
    strcat(fileName, ".chi_marg");
    pFile = fopen(fileName, "w");

    /* write the header */
    fprintf(pFile, "# Z1/Z2 pdf_marg Zref chiRef ");
    for (i=0; i<opts->nChiMethod; i++)
      fprintf(pFile, "chiRef-%s ", opts->chiMethod[i]);
    fprintf(pFile, "chi1-marg chi2-marg ");
    for (i=0; i<opts->nChiMethod; i++)
      fprintf(pFile, "chi1-marg-%s chi2-marg-%s ", opts->chiMethod[i],
              opts->chiMethod[i]);
    fprintf(pFile, "\n");

    for (j=0; j<n1; j++) {
      fprintf(pFile, "%8.6f %15.6e %15.6e ",x1[j], pdf_Z1[j], Z2ref[j]);
      fprintf(pFile, "%15.6e ", chiRef[1]);
      for (i=0; i<opts->nChiMethod; i++)
        fprintf(pFile, "%15.6e ", chi2_r_F[i][j]);
      fprintf(pFile, "%15.6e %15.6e ", chi1_m[1][j], chi2_m[1][j]);
      for (i=0; i<opts->nChiMethod; i++)
        fprintf(pFile,"%15.6e %15.6e ",chi1_m_F[i][1][j],chi2_m_F[i][1][j]);
      fprintf(pFile, "\n");
    }
    fprintf(pFile, "\n\n");
    for (k=0; k<n2; k++) {
      fprintf(pFile, "%8.6f %15.6e %15.6e ",x2[k], pdf_Z2[k], Z1ref[k]);
      fprintf(pFile, "%15.6e ", chiRef[2]);
      for (i=0; i<opts->nChiMethod; i++)
        fprintf(pFile, "%15.6e ", chi1_r_F[i][k]);
      fprintf(pFile, "%15.6e %15.6e ", chi1_m[2][k], chi2_m[2][k]);
      for (i=0; i<opts->nChiMethod; i++)
        fprintf(pFile,"%15.6e %15.6e ",chi1_m_F[i][2][k],chi2_m_F[i][2][k]);
      fprintf(pFile, "\n");
    }
    fclose(pFile);
    /* ----- finish writing file ----- */

    // free memory allocated in loop
    for (i=0; i<opts->nChiMethod; i++) {
      free(chi0_F[i]); free(chi1_F[i]); free(chi12_F[i]); free(chi2_F[i]);
      free(theta_F[i]);
      free(chi1_m_F[i][1]); free(chi1_m_F[i][2]);
      free(chi2_m_F[i][1]); free(chi2_m_F[i][2]);
      free(chi1_m_F[i]);    free(chi2_m_F[i]);

      free(chiRef_F[i]); free(chi1_r_F[i]); free(chi2_r_F[i]);
      free(JSD[i]);
    }
    free(chi0); free(theta); free(pdf_Z1); free(pdf_Z2);
    free(chi1_m[1]); free(chi1_m[2]); free(chi2_m[1]); free(chi2_m[2]);

    /* free(chi1Log); free(chi12Log); free(chi2Log); free(chi0Log); */
    /* free(chi1Erfc); free(chi12Erfc); free(chi2Erfc); free(chi0Erfc); */
    /* free(theta); free(thetaLog); free(thetaErfc); */
    destroyFlame(flame);
    destroyFlow(flow);
    destroyIO(io);
  }

  fclose(meanFile);
}
/* ------------------------------------------------------------------------- */

void heatReleaseBudget(optStruct *opts)
{
  int i, j, k, n1, n2, r;
  double *HR, *pdf, HRbudget[3], t; //

  // get the grid size from the first file
  n1 = n2 = 1;
  n1 = (int)(readFileScalar(opts->dataFile[0], "gridPoints"));
  n2 = (int)(readFileScalar(opts->dataFile[0], "gridPoints2"));

  // allocate memory
  HR  = (double*) calloc(n1*n2, sizeof(double));
  pdf = (double*) calloc(n1*n2, sizeof(double));

  // open the output file
  FILE *pFile = fopen("HRbudget.dat", "w");

  // write the file header
  fprintf(pFile, "# time   HR_pilot   HR_main   HR_mix\n");

  // loop over all files
  for (i=0; i<opts->nFiles; i++) {
    fprintf(stderr, "\t--> %s ... \n", opts->dataFile[i]);
    // get file time
    t = readFileScalar(opts->dataFile[i], "time");

    // read in the heat release and pdf for this time step
    readVariableArray(opts->dataFile[i], "HeatRelease", &HR);
    readVariableArray(opts->dataFile[i], "pdf", &pdf);

    // convolute to compute the heat release
    // heat release of the pilot is the integral along the Z1 axis
    HRbudget[1] = 0.0;
    for (j=0; j<n1; j++)
      HRbudget[1] += HR[j]*pdf[j];

    // heat release of the main is the integral along Z2 axis
    HRbudget[2] = 0.0;
    for (k=0; k<n1; k++)
      HRbudget[2] += HR[k*n1]*pdf[k*n1];

    // heat release of mixed regions is integral of region not on axis
    HRbudget[0] = 0.0;
    for (k=1; k<n2; k++)
      for (j=1; j<n1; j++)
        HRbudget[0] += HR[k*n1+j]*pdf[k*n1+j];

    // write to the output file
    fprintf(pFile, "%15.6e %15.6e %15.6e %15.6e\n", t, HRbudget[1],
            HRbudget[2], HRbudget[0]);
  }

  fclose(pFile);
  free(HR);
  free(pdf);
}
/* ------------------------------------------------------------------------- */

void transformFlow(optStruct *opts)
{
  // allocate pointers to internal structs
  ioStruct    *io        = (ioStruct*)    calloc(1, sizeof(ioStruct));
  driveStruct *driver    = (driveStruct*) calloc(1, sizeof(driveStruct));
  timeStruct  *time      = (timeStruct*)  calloc(1, sizeof(timeStruct));
  flowStruct  *flow      = (flowStruct*)  calloc(1, sizeof(flowStruct));
  flameStruct *flame     = (flameStruct*) calloc(1, sizeof(flameStruct));

  /* ----- set default values ----- */
  setDefaults(io, driver, time, flame, flow, NULL);

  // read in the flow file
  readFlamelet(opts->mducFile, io, driver, time, flame, flow);
  readFlowData(opts->mducFile, driver, flow);

  MDUCtransformVar(0,driver->TRANS_Z2,flame->n1,flame->n2,flame->x1,flame->x2,
                   &flow->chi1);
  MDUCtransformVar(0,driver->TRANS_Z2,flame->n1,flame->n2,flame->x1,flame->x2,
                   &flow->chi2);
  MDUCtransformVar(0,driver->TRANS_Z2,flame->n1,flame->n2,flame->x1,flame->x2,
                   &flow->pdf);

  writeFlowData("transformed.flo", io, driver, time, flame, flow);
}
/* ------------------------------------------------------------------------- */

void testMetric(optStruct *opts)
{
  int nGrids = 10;
  int i, j, jst;
  int nx;
  double *x, *dx;
  double *f, *f_x, *f_xx;
  double *f_x_a, *f_xx_a;
  double *err_x, *err_xx;
  double L2_x, L2_xx, tmp, L2_x0, L2_xx0;
  double L_inf_x, L_inf_xx, L_inf_x0, L_inf_xx0;
  double PI = 3.1415926535897931;

  flameStruct *flame = (flameStruct*) calloc(1, sizeof(flameStruct));

  flame->nDim  = 1;
  flame->n2    = 1;
  flame->order = opts->metricOrder;

  // output data files
  FILE *pFile = fopen("metricTest.dat", "w");
  fprintf(pFile, "# n L2_x L2_xx\n");
  FILE *dFile = fopen("metricFunc.dat", "w");
  fprintf(dFile, "# x f f_x_a f_xx_a f_x f_xx\n");

  for (i=1; i<=nGrids; i++) {
    // initialise grid and metric operators
    nx = flame->n1 = i*opts->nGrid;
    initGrid(nx,0,2*PI,opts->nPivot,opts->pivot,opts->gridPerc, &flame->x1,
             &flame->dx1);
    initMetric(flame);

    // allocate the memory for this grid
    f      = (double*) calloc(nx, sizeof(double)); 
    f_x    = (double*) calloc(nx, sizeof(double)); 
    f_xx   = (double*) calloc(nx, sizeof(double)); 
    f_x_a  = (double*) calloc(nx, sizeof(double)); 
    f_xx_a = (double*) calloc(nx, sizeof(double)); 
    err_x  = (double*) calloc(nx, sizeof(double)); 
    err_xx = (double*) calloc(nx, sizeof(double)); 

    // compute the analytic function and its derivatives
    for (j=0; j<nx; j++) {
      f[j]      =  sin(flame->x1[j]-PI/2);
      f_x_a[j]  =  cos(flame->x1[j]-PI/2);
      f_xx_a[j] = -sin(flame->x1[j]-PI/2);
      /* f[j]      =  sin(2*flame->x1[j]); */
      /* f_x_a[j]  =  2*cos(2*flame->x1[j]); */
      /* f_xx_a[j] = -4*sin(2*flame->x1[j]); */
    }

    // use the gradient operators
    L2_x = L2_xx = 0.0;
    L_inf_x = L_inf_xx = 0.0;
    for (j=0; j<nx; j++) {
      // compute gradients
      jst = stencil_x(j, flame);
      f_x[j]  = grad_x(jst, 0, j, 0, 1, f, flame);
      f_xx[j] = grad_xx(jst, 0, j, 0, 1, f, flame);

      // compute the residual, L2 and Linfinity norms
      tmp   = f_x[j] - f_x_a[j];
      if (j > 1 && j < nx-2) { // don't include boundary points
        L2_x += tmp*tmp;
        if (fabs(tmp) > L_inf_x) L_inf_x = fabs(tmp);
      }
      err_x[j] = fabs(tmp/f_x_a[j])*100; // percent error

      // for second derivative
      tmp   = f_xx[j] - f_xx_a[j];
      if (j > 1 && j < nx-2) { // don't include boundary points
        L2_xx += tmp*tmp;
        if (fabs(tmp) > L_inf_xx) L_inf_xx = fabs(tmp);
      }
      err_xx[j] = fabs(tmp/f_xx_a[j])*100; // percent error

      // write the function, derivatives and error for this grid to file
      fprintf(dFile, "%8.6f %15.6e %15.6e %15.6e %15.6e %15.6e\n",
              flame->x1[j], f[j], f_x_a[j], f_xx_a[j], f_x[j], f_xx[j]);
    }
    fprintf(dFile, "\n\n");
    L2_x  = sqrt(L2_x);
    L2_xx = sqrt(L2_xx);
    if (i == 1) { // store norms for first grid to normalize later
      L2_x0 = L2_x;
      L2_xx0 = L2_xx;
      L_inf_x0 = L_inf_x;
      L_inf_xx0 = L_inf_xx;
    }

    // write the norm data to file (normalised by initial)
    fprintf(pFile, "%i %15.6e %15.6e %15.6e %15.6e\n", flame->n1,
            L2_x/L2_x0, L2_xx/L2_xx0, L_inf_x/L_inf_x0, L_inf_xx/L_inf_xx0);

    free(f); free(f_x); free(f_xx); free(f_x_a); free(f_xx_a);
    free(flame->x1); free(flame->dx1);
    for (j=0; j<flame->n1; j++) {
      free(flame->grad.x[j]);
      free(flame->grad.xx[j]);
    }
    free(flame->grad.x);
    free(flame->grad.xx);
  } // end loop over grids

  fclose(pFile);
  fclose(dFile);
}
/* ------------------------------------------------------------------------- */

void checkThermo(optStruct *opts, chemStruct *chem, flameStruct *flame)
{
  void *mem;
  // chemStruct chem;
  // initChem(&chem);
  //mem = MDUCinit(opts->mducFile);

  fprintf(stderr, "using mechanism with %i species\n", chem->nSpecies);

  fprintf(stderr, "Initial temperature: %6.1f K\n", flame->Y[flame->sT]);

  int i,j;
  double T  = flame->Y[flame->sT];
  double Tref = 298.15;

  int nT = 100;
  int Tmin, Tmax, dT;
  double *Trange = (double*) calloc(nT, sizeof(double));
  double totMol;
  double *Y   = (double*) calloc(chem->nSpecies, sizeof(double));
  double *X   = (double*) calloc(chem->nSpecies, sizeof(double));
  double *h   = (double*) calloc(chem->nSpecies, sizeof(double));
  double *cp  = (double*) calloc(chem->nSpecies, sizeof(double));
  double *hr  = (double*) calloc(chem->nSpecies, sizeof(double));
  double *cpr = (double*) calloc(chem->nSpecies, sizeof(double));

  Trange[0] = Tmin = Tref;
  Tmax = 2000;
  dT = (Tmax-Tmin)/(nT-1);
  for (i=1; i<nT; i++) Trange[i] = Trange[i-1] + dT;
  
  totMol = 0.0;
  for (i=0; i<chem->nSpecies; i++) {
    Y[i] = flame->Y[i];
    X[i] = Y[i]/chem->W[i];
    totMol += X[i];
  }

  for (i=0; i<chem->nSpecies; i++) X[i] /= totMol;

  /*   for (i=0; i<chem->nSpecies; i++) */
  /*     fprintf(stderr, "%10s %8.6f %8.6f\n", chem->species[i], Y[i], X[i]); */

  /* compute the total enthalpy */  
  chemComputeThermoData(h, cp, T, chem);
  chemComputeThermoData(hr, cpr, Tref, chem);

  double formH = 0.0, cpmix = 0.0, sensH = 0.0, totH = 0.0;
  double wmix = 0.0, cpmix_r = 0.0;
  double cpavg, newH, newT;
  for (i=0; i<chem->nSpecies; i++) {
    wmix    += X[i]*chem->W[i];
    totH    += Y[i]*h[i];
    formH   += Y[i]*hr[i];
    cpmix   += Y[i]*cp[i];
    cpmix_r += Y[i]*cpr[i];
    sensH   += Y[i]*(h[i]-hr[i]);
  }
  cpavg = sensH/(T-Tref);
  newH  = 0.5*(cpmix+cpmix_r)*(T-Tref)+formH;
  //  newH  = cpmix*(T-Tref)+formH;

  fprintf(stderr, "Mixture mol. weight:   %15.3f g/mol\n", wmix);
  fprintf(stderr, "Total enthalpy:        %15.6e J/kg\n", totH);
  fprintf(stderr, "Sensible enthalpy:     %15.6e J/kg\n", sensH);
  fprintf(stderr, "Enthalpy of formation: %15.6e J/kg\n", formH);
  fprintf(stderr, "Mixture specific heat: %15.6e J/(kg-K)\n", cpmix);
    
  fprintf(stderr, "cp*(T-Tref):           %15.6e J/kg\n", cpmix*(T-Tref));
  fprintf(stderr, "Average specific heat: %15.6e J/(kg-K)\n", cpavg);
  fprintf(stderr, "New total enthalpy:    %15.6e J/kg\n", newH);

  FILE *pFile = fopen("therm_compare.dat", "w");

  // write the header
  fprintf(pFile, "T[K] cpmix cpavg error_cp  formH sensH totH newH errorH\n");

  for (j=0; j<nT; j++) {
    chemComputeThermoData(h, cp, Trange[j], chem);
    cpmix = totH = sensH = formH = 0.0;
    for (i=0; i<chem->nSpecies; i++) {
      totH  += Y[i]*h[i];
      formH += Y[i]*hr[i];
      cpmix += Y[i]*cp[i];
      sensH += Y[i]*(h[i]-hr[i]);
    }
    cpavg = sensH/(Trange[j]-Tref);
    newH  = cpmix*(Trange[j]-Tref)+formH;

    fprintf(pFile, "%15.6e %15.6e %15.6e %15.6e ",
            Trange[j], cpmix, cpavg, fabs(cpmix-cpavg)/cpmix*100);
    fprintf(pFile, "%15.6e %15.6e %15.6e %15.6e %15.6e\n",
            formH, sensH, totH, newH, fabs(newH-totH)/totH*100);
  }

  fclose(pFile);


  free(Y); free(X); free(h); free(cp); free(hr); free(cpr); free(Trange);
}
/* ------------------------------------------------------------------------- */

void genFireSpeciesList(optStruct *opts)
{
  int i,j, nSpecies;
  char **species;
  char **speciesH;

  char buffer[100], tag[20];

  FILE *pFile = fopen(opts->mducFile, "r");
  FILE *hFile = fopen(opts->monitorFile, "r");
  FILE *mechFile = fopen("mech_names.out", "w");
  FILE *fireFile = fopen("fire_names.out", "w");

  fscanf(pFile, "%i", &nSpecies);
  for (i=0; i<8; i++)
    fgets(buffer, 100, hFile);
  
  species  = (char**) calloc(nSpecies, sizeof(char*));
  speciesH = (char**) calloc(nSpecies, sizeof(char*));
  for (i=0; i<nSpecies; i++) {
    species[i]  = (char*) calloc(25, sizeof(char));
    speciesH[i] = (char*) calloc(25, sizeof(char));
    fscanf(pFile, "%s", species[i]);
    fscanf(hFile, "%s", speciesH[i]);
    fgets(buffer, 100, hFile);
    fprintf(fireFile, "%s\n", species[i]);
    for (j=0; j<25; j++) species[i][j] = toupper(species[i][j]);
    fprintf(mechFile, "strcpy( names[%s], \"%-25s\" );\n",
            speciesH[i], species[i]);
  }

  fclose(pFile); fclose(fireFile); fclose(mechFile); fclose(hFile);
}
/* ------------------------------------------------------------------------- */

void calcInduction(optStruct *opts)
{
  double PI     = 3.141592653589793;
  engineStruct *engine = (engineStruct*) calloc(1, sizeof(engineStruct));

  engine->RPM    = 1500;
  engine->CAD_0  = 420;
  engine->bore   = 0.0820;
  engine->stroke = 0.0850;
  engine->rod    = 0.1435;
  engine->CR     = 13;

  int i, j, jL, nt, nAv;
  double *time, *CAD, *mass, *Pc, *Pi, *Av, *AvIn, *vol, *CAD_Av, *gamma;

  double *mdot_CFD, *mdot_steady, *mdot_unsteady, *mass_unsteady, **F;
  double *u1;

  double Dp = 0.030;
  double Ai = 0.25*PI*Dp*Dp;
  double L  = 500.0;
  double a, a_2, g, R, P0, rhoi;
  double Ti = 300;

  /* intake speed of sound */
  P0 = 101325;
  g = 1.4;
  R = 287;
  a_2 = g*R*Ti;
  a = pow(a_2, 0.5);
  rhoi = 1.2;

  /* read monitor files */
  nt = readMonitorFile(opts->dataFile[0], opts->iT, &time, opts->iV, &mass);
  free(time);
  nt = readMonitorFile(opts->dataFile[1], opts->iT, &time, opts->iV, &Pc);
  free(time);
  nt  = readMonitorFile(opts->dataFile[2], opts->iT, &time, opts->iV, &Pi);
/*   free(time); */
/*   nt  = readMonitorFile(opts->dataFile[3], opts->iT, &time, opts->iV, &gamma); */
  nAv = readMonitorFile(opts->dataFile[3], opts->iT, &CAD_Av, opts->iV, &AvIn);

  /* convert time values to CAD */
  vol = (double*) calloc(nt, sizeof(double));
  CAD = (double*) calloc(nt, sizeof(double));
  Av  = (double*) calloc(nt, sizeof(double));

  mdot_CFD      = (double*) calloc(nt, sizeof(double));
  mdot_steady   = (double*) calloc(nt, sizeof(double));
  mdot_unsteady = (double*) calloc(nt, sizeof(double));
  mass_unsteady = (double*) calloc(nt, sizeof(double));
  u1            = (double*) calloc(nt, sizeof(double));
  F    = (double**) calloc(nt, sizeof(double*));
  for (i=0; i<nt; i++)
    F[i] = (double*) calloc(3, sizeof(double));

  for (i=0; i<nt; i++)
    {
      CAD[i] = time[i]*6*engine->RPM + engine->CAD_0;
      vol[i] = calcCylVol(CAD[i], engine);

      /* interpolate the valve area profile */
      jL = gridLocate(nAv, CAD_Av, CAD[i]);
      Av[i] = linInterp(CAD[i], CAD_Av[jL], AvIn[jL], CAD_Av[jL+1], AvIn[jL+1]);
    }

  u1[0] = pow(a_2*2.0/(g-1.0)*(pow(P0/Pi[0],(g-1)/g) - 1), 0.5);
  u1[0] = 0.1;
  
  FILE *pFile = fopen("induction_out", "w");

  mass_unsteady[0] = mass[0];

  /* loop over time */
  for (i=1; i<nt; i++)
    {
      /* compute cfd mass flow */
      mdot_CFD[i] = (mass[i] - mass[i-1])/(time[i]-time[i-1]);

      /* compute steady mass flow */
      if (Pi[i] > Pc[i])
        mdot_steady[i] = calcValveSteadyMassFlow(g, rhoi, Pi[i], Pc[i], Av[i]);
      else
        mdot_steady[i] = -calcValveSteadyMassFlow(g, rhoi, Pc[i], Pi[i], Av[i]);
      
      /* compute unsteady mass flow */
/*       mass_unsteady[i] = mass_unsteady[i-1] + (time[i]-time[i-1])*(1.0*Av[i]/(Ai*L)*(Pi[i]*Av[i] - Pc[i]*Av[i] - a*a*Av[i]/vol[i]*(mass_unsteady[i] - mass_unsteady[0]))); */
/*       mass_unsteady[i] = mass_unsteady[i-1] + (time[i]-time[i-1])*(1.0*Av[i]/(Ai*L)*(Pi[i]*Av[i] - Pc[i]*Av[i] - a*a*Av[i]/vol[i]*mdot_unsteady[i]*(time[i]-time[i-1]))); */

      /* ----- most current version ----- */
/*       F[i][0] = Pi[i]*Av[i] - Pc[i]*Av[i]; /\* forcing term *\/ */
/*       F[i][1] = a*a*Av[i]/vol[i]*(mass_unsteady[i] - mass_unsteady[0]); /\* inertial *\/ */
/*       mdot_unsteady[i] = mdot_unsteady[i-1] + (time[i]-time[i-1])*(1.0*Av[i]/(Ai*L)*(F[i][0] - F[i][1])); */

/* /\*       mdot_unsteady[i] = (mass_unsteady[i] - mass_unsteady[i-1])/(time[i]-time[i-1]); *\/ */
/*       mass_unsteady[i] = mass_unsteady[i] + (time[i]-time[i-1])*mdot_unsteady[i]; */
      /* ----- end ----- */

      F[i][0] = (Pi[i] - Pc[i])*Av[i]; /* forcing term */
      F[i][1] = a_2*Av[i]/vol[i]*(mass_unsteady[i] - mass_unsteady[0]); /* inertial */
      u1[i] = u1[i-1] + (time[i]-time[i-1])*1.0/(Ai*1.2*L)*(F[i][0] - F[i][1]);

      mdot_unsteady[i] = 1.2*u1[i]*Av[i];

      mass_unsteady[i] = mass_unsteady[i-1]+(time[i]-time[i-1])*mdot_unsteady[i];


      fprintf(pFile, "%15.6e %15.6e %15.6e %15.6e %15.6e %15.6e ",
              time[i], CAD[i], vol[i], Pc[i], Pi[i], Av[i]);
      fprintf(pFile, "%15.6e %15.6e %15.6e %15.6e %15.6e ", mdot_CFD[i], mdot_steady[i],
              mdot_unsteady[i], mass_unsteady[i], u1[i]);
      for (j=0; j<2; j++)
        fprintf(pFile, "%15.6e  ", F[i][j]);
      fprintf(pFile, "\n");

    }

  /* write the results */
  fclose(pFile);

}
/* ------------------------------------------------------------------------- */

void integrateHeatRelease(optStruct *opts)
{
  double PI     = 3.141592653589793;
  engineStruct *engine = (engineStruct*) calloc(1, sizeof(engineStruct));

  engine->RPM    = 1500;
  engine->CAD_0  = 420;
  engine->bore   = 0.0820;
  engine->stroke = 0.0850;
  engine->rod    = 0.1435;
  engine->CR     = 13;

  int i, j, jL, nt, nAv, nt_dum;
  double *time, *Q_h, *CAD, *mass, *Pc, *Pi, *Av, *AvIn, *vol, *CAD_Av, *gamma, *t_dum;

  double *mdot_CFD, *mdot_steady, *mdot_unsteady, *mass_unsteady, **F;
  double *u1;

  double Dp = 0.030;
  double Ai = 0.25*PI*Dp*Dp;
  double L  = 500.0;
  double a, a_2, g, R, P0, rhoi;
  double Ti = 300;
  double dt, dVdt, dPdt, Q_h_tot;

  /* intake speed of sound */
  P0 = 101325;
  g = 1.4;
  R = 287;
  a_2 = g*R*Ti;
  a = pow(a_2, 0.5);
  rhoi = 1.2;

  /* read monitor files */
  opts->iT = 1;
  opts->iV = 2;

  nt = readMonitorFile(opts->dataFile[0], opts->iT, &time, opts->iV, &Pc);
  //  nt_dum = readMonitorFile(opts->dataFile[0], opts->iT, &t_dum, 10, &gamma);
  nt_dum = readMonitorFile(opts->dataFile[1], opts->iT, &t_dum, opts->iV, &gamma);

  vol = (double*) calloc(nt, sizeof(double));
  CAD = (double*) calloc(nt, sizeof(double));
  Q_h = (double*) calloc(nt, sizeof(double));
  //  gamma = (double*) calloc(nt, sizeof(double));

  for (i=0; i<nt; i++) {
    CAD[i] = time[i]*6*engine->RPM + engine->CAD_0;
    vol[i] = calcCylVol(CAD[i], engine);
  }

  /* write the header */

  FILE *pFile = fopen("heat_release.out", "w");

  fprintf(pFile, "# time[s]<1> CAD<2> Pc[Pa]<3> Q_h[J/s]<4> gamma<5>\n");

  Q_h_tot = 0.0;
  for (i=1; i<nt; i++) {
    dt = (time[i]-time[i-1]);
    dPdt = (Pc[i] - Pc[i-1])/dt;
    dVdt = (vol[i] - vol[i-1])/dt;
    /* if (CAD[i] < 720) */
    /*         gamma[i] = 1.32; */
    /*       else if (CAD[i] >= 720) */
    /*         gamma[i] = 1.27; */
    //      gamma[i] = log(Pc[i]/Pc[i-1])/log(vol[i-1]/vol[i]);
    Q_h[i] = gamma[i]/(gamma[i]-1.0)*Pc[i]*dVdt + 1.0/(gamma[i]-1.0)*vol[i]*dPdt;

    Q_h_tot += fmax(0.0,Q_h[i])*dt;
    fprintf(pFile, "%15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n",
            time[i], CAD[i], Pc[i], Q_h[i], gamma[i], dPdt, dVdt);
  }

  fprintf(stderr, "integrated heat release: %15.6e\n", Q_h_tot);



  /* write the results */
  fclose(pFile);

}
/* ------------------------------------------------------------------------- */

void printPoly(int na, gPoint * a, int nb, gPoint * b)
{
  printf("%11.7f, %11.7f\n", polyIntersect(na, a, nb, b), polyIntersect(na, a, na, a));
}

void testPolyA() {
  { gPoint a[] = {{2,3}, {2,3}, {2,3}, {2,4}, {3,3}, {2,3}, {2,3}},
      b[] = {{1,1}, {1,4}, {4,4}, {4,1}, {1,1}};
      printPoly(7, a, 5, b);}// 1/2, 1/2
  // The redundant vertices above are to provoke errors as good test cases should.
  // It is not necessary to duplicate the first vertex at the end.
  {gPoint  a[] = {{1,7}, {4,7}, {4, 6}, {2,6}, {2, 3}, {4,3}, {4,2}, {1,2}},
      b[] = {{3,1}, {5,1}, {5,4}, {3,4}, {3,5}, {6,5}, {6,0}, {3,0}};
      printPoly(8, a, 8, b);}// 0, 9
  {gPoint  a[] = {{1,1}, {1,2}, {2,1}, {2,2}},
      b[] = {{0,0}, {0,4}, {4,4}, {4,0}}; printPoly(4, a, 4, b);}// 0, 1/2
  {gPoint  a[] = {{0,0}, {3,0}, {3,2}, {1,2}, {1,1}, {2,1}, {2,3}, {0,3}},
      b[] = {{0,0}, {0,4}, {4,4}, {4,0}}; printPoly(8, a, 4, b);}// -9, 11
  {gPoint  a[] = {{0,0}, {1,0}, {0,1}},
      b[] = {{0,0}, {0,1}, {1,1}, {1,0}}; printPoly(3, a, 4, b);}// -1/2, 1/2
  {gPoint  a[] = {{1,3}, {2,3}, {2,0}, {1,0}},
      b[] = {{0,1}, {3,1}, {3,2}, {0,2}}; printPoly(4, a, 4, b);}// -1, 3
  {gPoint  a[] ={{0,0}, {0,2}, {2,2}, {2,0}},
      b[] = {{1, 1}, {3, 1}, {3, 3}, {1, 3}}; printPoly(4, a, 4, b);}// -1, 4
  {gPoint  a[] = {{0,0}, {0,4}, {4,4}, {4,0}},
      b[] = {{1,1}, {1,2}, {2,2}, {2,1}}; printPoly(4, a, 4, b);}// 1, 16
}
/* -------------------------------------------------------------------------- */

void flameletFileAddWmol(optStruct *opts)
{
  int i,j,k,iF,col;
  double *Wmol, *Y;

  // allocate pointers to internal structs
  ioStruct    *io        = (ioStruct*)    calloc(1, sizeof(ioStruct));
  driveStruct *driver    = (driveStruct*) calloc(1, sizeof(driveStruct));
  timeStruct  *time      = (timeStruct*)  calloc(1, sizeof(timeStruct));
  flowStruct  *flow      = (flowStruct*)  calloc(1, sizeof(flowStruct));
  flameStruct *flame     = (flameStruct*) calloc(1, sizeof(flameStruct));

  /* initialize the mechanism */
  io->parser = parseFile(opts->mducFile);
  chemStruct *chem = chemInit("compiled",io->parser);
  printSpeciesList(chem);

  /* loop over all the files specified */
  for (iF=0; iF<opts->nFiles; iF++) {

    fprintf(stderr, "Adding Molecular weight to file: %s\n",opts->dataFile[iF]);

    /* ----- set default values and read in the flamelet ----- */
    setDefaults(io, driver, time, flame, flow, NULL);
    readFlamelet(opts->dataFile[iF], io, driver, time, flame, flow);

    io->inputIndex = (int*) calloc(io->nInputs, sizeof(int));
    for (i=0; i<io->nInputs; ++i) {
      io->inputIndex[i] = getSpeciesIndex(chem->nSpecies, chem->species, 
                                          io->inputSpecies[i]);
      if (io->inputIndex[i] < 0)
        fprintf(stderr, "Warning: species %s not found\n", io->inputSpecies[i]);
    }

    Wmol = (double*) calloc(flame->n1, sizeof(double));
    Y    = (double*) calloc(flame->n1*flame->n2*chem->nSpecies, sizeof(double));

    /* store the species in a standard array format */
    /* TODO: only works for 1D at the moment */
    for (j=0; j<flame->n1; j++) {
      for (i=0; i<io->nInputs; ++i)
        Y[j*chem->nSpecies+io->inputIndex[i]] = io->initY[i][j];
    } /* end loop over columns */

    /* compute the molecular weight */
    for (j=0; j<flame->n1; j++) {
      Wmol[j] = calcMixtureW(chem->nSpecies, &Y[j*chem->nSpecies], chem->W);
    }

    addVarToBody(opts->dataFile[iF],flame->n1,1,Wmol,"WMIX");

    destroyFlame(flame);
    destroyFlow(flow);
    destroyIO(io);
    free(Wmol); free(Y);
  }
}
/* -------------------------------------------------------------------------- */
