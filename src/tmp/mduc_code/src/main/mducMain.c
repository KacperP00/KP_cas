#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "mducMain.h"
#include "mducTools.h"
#include "accessMDUC.h"
#include "postProcess.h"
#include "fileIO.h"
#include "parser.h"
#include "display.h"
#include "steadyFlamelet.h"
#include "vtkIO.h"
#include "iFileIO.h"

void printHelp()
{
  printf("\n");
  printf("options :\n");
  printf("  -h                       # print help\n");
  printf("  -i <input file>          # input file/flamelet file \n");
  printf("  -o <output file>         # output filename\n");
  printf("  -d <data file(s)>        # data files to modify \n");
  printf("  -m <file>                # monitor file\n");
  printf("  -n <int>                 # grid size (various)\n");
  printf("  -f <int>                 # for frequency\n");
  printf("  -g <pivot> ...           # generate grid\n");
  printf("  -p <percentage> ...      # grid point distribution\n");
  printf("commands :\n");
  printf("  --post-process <file>    # post-process specified file\n");
  printf("  --fuel-temp              # compute liquid fuel temp\n");
  printf("  --interpolate <n1> <n2>  # interpolate solution to new grid\n");
  printf("  --init <1> (<2>)         # init 2D field from 2 1D solutions\n");
  printf("  --transform-flow <file>  # transform flow file\n");
  printf("  --convert-1D <dir> <n2>  # convert 1D solution to 2D with n\n");
  printf("  --data2ensight <dir>     # write specified data to ensight <dir>\n");
  printf("  --ens-case <file>        # case file to read time array from\n");
  printf("  --ens-time               # convert time from monitor to ensight\n");
  printf("  --range <int> <int>      # only use a specified range of data\n");
  printf("  --ens-const              # convert monitor file to ensight const\n");
  printf("  --ens-conv               # convert particle file to binary\n");
  printf("  --integrate-pdf\n");
  printf("  --hr-budget\n");
  printf("  --metric\n");
  printf("  --clean\n");
  printf("  --output-species-list\n");
  printf("  --FIRE-species-list\n");
  printf("  --vtk-check\n");
  printf("  --add-wmol\n");
  printf("  --ifile-info\n");
  printf("  --test-parser\n");
  exit(0);
}
/* -------------------------------------------------------------------------- */

void setDefaultOptions(optStruct *opts)
{
  // set default file names
  strcpy(opts->mducFile, "mduc.inp");
  strcpy(opts->outFile, "none");
  strcpy(opts->postProcFile, "none");
  strcpy(opts->monitorFile, "none");
  strcpy(opts->ensCase, "none");

  opts->postProcess  = 1;  // by default, post-process CA.in
  opts->fuelTemp     = 0;
  opts->interp       = 0;
  opts->init2D       = 0;
  opts->convert1D    = 0;
  opts->writeEnsight = 0;
  opts->ensConst     = 0;
  opts->ensConv      = 0;
  opts->ensTime      = 0;
  opts->checkGrid    = 0;
  opts->nGrid        = 101;
  opts->integratePDF = 0;
  opts->checkChi     = 0;
  opts->postChi      = 0;
  opts->tFlow        = 0;
  opts->freq         = 2;
  opts->hrBudg       = 0;
  opts->clean        = 0;
  opts->combine      = 0;
  opts->iT           = 1;
  opts->iV           = 2;
  opts->metric       = 0;
  opts->polyA        = 0;
  opts->multi1D      = 0;
  opts->chiProfile   = 0;
  opts->checkThermo  = 0;
  opts->listSpecies  = 0;
  opts->fireSpecies  = 0;
  opts->checkVTK     = 0;
  opts->addWmol      = 0;
  opts->iFile        = 0;
  opts->testParse    = 0;
}
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  int i, j, k;
  optStruct opts; // struct to store options from command line
  void *workMDUC;

  // declare flamelet structs so that internal functions can be used
  driveStruct *driver;
  ioStruct    *io;
  chemStruct  *chem;
  timeStruct  *time;
  flameStruct *flame;
  flowStruct  *flow;

#ifdef ACCESS_MPI
  // call MPI_Init before handling any other options
  MPI_Init(&argc, &argv);
  MDUCsetPrintProc(MPI_COMM_WORLD);
#else
  MDUCsetPrintProc(NULL);
#endif // ACCESS_MPI

  // set options
  setDefaultOptions(&opts);
  readOptions(argc, argv, &opts);
 
  /* ----- options that don't require initialisation/input file ----- */
  if (opts.init2D){
    void *work1, *work2;
    if (opts.nFiles == 1) {
      printf("initialising flamelet using %s:\n", opts.dataFile[0]);
      work1 = MDUCinit(opts.dataFile[0]);
      MDUCwriteFlamelet("init_solution.rif", work1);
    }
    else if (opts.nFiles == 2) {
      printf("initialising a 2D field from two 1D flamelets:\n");
      printf(" flamelet 1: %s\n", opts.dataFile[0]);
      printf(" flamelet 2: %s\n", opts.dataFile[1]);
      
      // read in the two 1D files to initialise from
      work1 = MDUCinit(opts.dataFile[0]);
      work2 = MDUCinit(opts.dataFile[1]);
      
      // initialise and output the 2D field
      MDUCinit2Dfrom_1D_and_1D(work2, work1);
      printf("writing flamelet to %s\n", "init2D_solution.rif");
      MDUCwriteFlamelet("init_solution.rif", work1);
      
      MDUCdestroy(work1);
      MDUCdestroy(work2);
    }
    return 0;
  }
  else if (opts.convert1D) {
    convert1D_2D(&opts);
    return 0;
  }
  else if (opts.writeEnsight) {
    mduc2ensight(&opts);
    return 0;
  }
  else if (opts.ensConst) {
    writeEnsightConst(&opts, opts.nGrid, opts.freq);
    return 0;
  }
  else if (opts.ensConv) {
    convertEnsightParticle(&opts);
    return 0;
  }
  else if (opts.ensTime) {
    genEnsightTime(&opts);
    return 0;
  }
  else if (opts.interp) {
    interpolateSolution(&opts);
    return 0;
  }
  else if (opts.metric) {
    testMetric(&opts);
    return 0;
  }
  else if (opts.polyA) {
    testPolyA();
    return 0;
  }
  else if (opts.chiProfile) {
    chiProfile(&opts);
    return 0;
  }
  else if (opts.checkGrid) {
    checkGrid(&opts);
    return 0;
  }
  else if (opts.integratePDF) {
    getZfromPDF(&opts);
    return 0;
  }
  else if (opts.hrBudg) {
    heatReleaseBudget(&opts);
    return 0;
  }
  else if (opts.tFlow) {
    transformFlow(&opts);
    return 0;
  }
  else if (opts.clean) {
    cleanMonitorFile(&opts);
    return 0;
  }
  else if (opts.combine) {
    combineMonitorFiles(&opts);
    return 0;
  }
  else if (opts.postChi) {
    postProcessChi(&opts);
    return 0;
  }
  else if (opts.fireSpecies) {
    genFireSpeciesList(&opts);
    return 0;
  }
  else if (opts.checkVTK) {
    readVTK(opts.mducFile);
    return 0;
  }
  else if (opts.addWmol) {
    flameletFileAddWmol(&opts);
    return 0;
  }
  else if (opts.iFile) {
    iFileStruct *iFile;
    iFile = iFileRead(opts.mducFile);
    return 0;
  }
  else if (opts.testParse) {
    parseStruct *parser;
    parser = parseFile(opts.mducFile);
    parsePrint(parser);
    char testString[32];
    double testR = paramGetReal("non-exist-r", 1.1, parser);
    fprintf(stderr, "\n\nnon-exist-r: %f\n", testR);
    testR = paramGetReal("var-2", 1.1, parser);
    fprintf(stderr, "var-2: %f\n", testR);
    paramGetString(testString, "var-1", "default string", parser);
    fprintf(stderr, "string test: %s\n", testString);
    return 0;
  }
  else if (opts.multi1D) {
    void *work, *work1, *work2;
    work  = MDUCinit(opts.dataFile[0]);
    work1 = MDUCinit(opts.dataFile[1]);
    work2 = MDUCinit(opts.dataFile[2]);
    getMDUCStructs(&io, &driver, &chem, &time, &flow, &flame, work);
    postProcessMulti1D(opts.postProcFile,io,driver,work,work1,work2);
    return 0;
  }
  

  // get the data from the input file
  workMDUC = MDUCinit(opts.mducFile);
  getMDUCStructs(&io, &driver, &chem, &time, &flow, &flame, workMDUC);

#ifdef ACCESS_MPI
  // intialize parallel MDUC
  MDUCsetParallel(MPI_COMM_WORLD, workMDUC);
#else
  // TODO: Write a function that intializes variables in case of serial
  flame->j0m = 0;
  flame->jNp = flame->n1;
#endif // ACCESS_MPI

  if (opts.fuelTemp) {
    double Tfuel;
    fprintf(stderr, "Computing fuel vapour temperature...\n");
    Tfuel = calcLiquidFuelTemp(io->fuelT1,io->hVap1,io->cpLiq1,io->fuelTref1,
                               chem->nSpecies,io->nFuel1,io->fuelIndex1,
                               io->fuelY1, chem);
    fprintf(stderr, "fuel temp = %8.4f\n", Tfuel);
    return 0;
  }
  else if (opts.checkThermo) {
    checkThermo(&opts, chem, flame);
    return 0;
  }
  else if (opts.listSpecies) {
    printSpeciesList(chem);
    return 0;
  }
  else if (opts.checkChi) {
    checkChi2D(&opts, workMDUC);
    return 0;
  }

  if (opts.postProcess) {
    if (driver->TIME_MODE == -1) { /* steady flamelet */
      if (flame->nDim < 2) {
        solveSteadyFlamelet1D(workMDUC);
      }
      else if (flame->nDim == 2) {
        solveSteadyFlamelet2D(workMDUC);
      }
      return 0;
    } /* end steady solution */

    if (flame->nDim < 2)
      postProcess1D(opts.postProcFile, io, driver, workMDUC);
    else if (flame->nDim == 2)
      postProcess2D(opts.postProcFile, io, driver, workMDUC);

    return 0;
  }

  return 0;
}
/* -------------------------------------------------------------------------- */

/* functions to retrieve internal structs from solver */
ioStruct     * MDUCgetIOStruct(void *mem)
{ mducStruct *mduc = (mducStruct*) mem; return mduc->io; }

driveStruct  * MDUCgetDriveStruct(void *mem)
{ mducStruct *mduc = (mducStruct*) mem; return mduc->driver; }

timeStruct   * MDUCgetTimeStruct(void *mem)
{ mducStruct *mduc = (mducStruct*) mem; return mduc->time; }

chemStruct   * MDUCgetChemStruct(void *mem)
{ mducStruct *mduc = (mducStruct*) mem; return mduc->chem; }

flameStruct  * MDUCgetFlameStruct(void *mem)
{ mducStruct *mduc = (mducStruct*) mem; return mduc->flame; }

flowStruct   * MDUCgetFlowStruct(void *mem)
{ mducStruct *mduc = (mducStruct*) mem; return mduc->flow; }

engineStruct * MDUCgetEngineStruct(void *mem)
{ mducStruct *mduc = (mducStruct*) mem; return mduc->engine; }

void getMDUCStructs(ioStruct **io, driveStruct **driver, chemStruct **chem,
                    timeStruct **time, flowStruct **flow, flameStruct **flame,
                    void *mem)
{ // dereference all internal flamelet solver structs
  mducStruct *mduc = (mducStruct*) mem;
  io[0]     = mduc->io;
  driver[0] = mduc->driver;
  chem[0]   = mduc->chem;
  time[0]   = mduc->time;
  flow[0]   = mduc->flow;
  flame[0]  = mduc->flame;
}
/* -------------------------------------------------------------------------- */

#ifdef ACCESS_MPI
void getMPIStruct(mpiStruct **mpi, void *mem)
{
  mducStruct *mduc = (mducStruct*) mem;
  mpi[0] = mduc->mpi;
}
/* -------------------------------------------------------------------------- */
#endif // ACCESS_MPI

void readList(int n, int iStart, char ***var, char **argv)
{
  int i, j;
  
  var[0] = (char**) calloc(n, sizeof(char*));
  for (i=0, j=iStart; i<n; i++, j++) {
    var[0][i] = (char*) calloc(strlen(argv[j])+10, sizeof(char));
    strcpy(var[0][i], argv[j]);
  }
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
  int i, j, k, nOpt, iOpt, iEnd;
  int *optLocs;
  
  // read through the argument list and get the option locations
  optLocs    = (int*)  calloc(argc, sizeof(int));
  opts->args = (char*) calloc(argc, sizeof(char));
  nOpt = 0;
  for (i=1; i<argc; i++) {
    if (argv[i][0] == '-') {
      optLocs[nOpt] = i;
      opts->args[nOpt] = argv[i][1];
      nOpt++;
    }
  } // end loop over arguments
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
      strcpy(opts->mducFile, argv[iOpt+1]);
      break;
    case 'I':
      opts->nFiles = iEnd - iOpt;
      readList(opts->nFiles, iOpt+1, &opts->dataFile, argv);
      opts->multi1D = 1;
      break;
    case 'o':
      strcpy(opts->outFile, argv[iOpt+1]);
      break;
    case 'd':
      opts->nFiles = iEnd - iOpt;
      readList(opts->nFiles, iOpt+1, &opts->dataFile, argv);
      break;
    case 'g':
      opts->checkGrid = 1;
      opts->nPivot = iEnd - iOpt;
      readArray(opts->nPivot, iOpt+1, &opts->pivot, argv);
      break;
    case 'p':
      readArray(opts->nPivot, iOpt+1, &opts->gridPerc, argv);
      break;
    case 'n':
      opts->nGrid = atoi(argv[iOpt+1]);
      break;
    case 'f':
      opts->freq = atoi(argv[iOpt+1]);
      break;
    case 'm':
      strcpy(opts->monitorFile, argv[iOpt+1]);
      break;
    case '-': {// has a longer word
      if (!strcmp(argv[optLocs[i]], "--post-process")) {
        opts->postProcess = 1;
        strcpy(opts->postProcFile, argv[iOpt+1]);
      }
      if (!strcmp(argv[optLocs[i]], "--fuel-temp")) {
        opts->fuelTemp = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--interpolate")) {
        opts->interp = 1;
        opts->n1 = atoi(argv[iOpt+1]);
        opts->n2 = atoi(argv[iOpt+2]);
      }
      if (!strcmp(argv[optLocs[i]], "--init")) {
        opts->init2D = 1;
        opts->nFiles = iEnd - iOpt;
        readList(opts->nFiles, iOpt+1, &opts->dataFile, argv);
      }
      if (!strcmp(argv[optLocs[i]], "--convert-1D")) {
        opts->convert1D = 1;
        opts->dir = atoi(argv[iOpt+1]);
        opts->n2  = atoi(argv[iOpt+2]);
      }
      if (!strcmp(argv[optLocs[i]], "--data2ensight")) {
        opts->writeEnsight = 1;
        strcpy(opts->ensDir, argv[iOpt+1]);
      }
      if (!strcmp(argv[optLocs[i]], "--range")) {
        opts->range[0] = atoi(argv[iOpt+1]);
        opts->range[1] = atoi(argv[iOpt+2]);
      }
      if (!strcmp(argv[optLocs[i]], "--integrate-pdf")) {
        opts->integratePDF = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--hr-budget")) {
        opts->hrBudg = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--transform-flow")) {
        opts->tFlow = 1;
        strcpy(opts->mducFile, argv[iOpt+1]);
      }
      if (!strcmp(argv[optLocs[i]], "--check-chi")) {
        opts->checkChi = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--chi-profile")) {
        opts->chiProfile = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--ens-const")) {
        opts->ensConst = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--ens-conv")) {
        opts->ensConv = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--ens-case")) {
        strcpy(opts->ensCase, argv[iOpt+1]);
      }
      if (!strcmp(argv[optLocs[i]], "--ens-time")) {
        opts->ensTime = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--clean")) {
        opts->clean = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--combine")) {
        opts->combine = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--time-index")) {
        opts->iT = atoi(argv[iOpt+1]);
      }
      if (!strcmp(argv[optLocs[i]], "--var-index")) {
        opts->iV = atoi(argv[iOpt+1]);
      }
      if (!strcmp(argv[optLocs[i]], "--post-chi")) {
        opts->postChi = 1;
        opts->nChiMethod = iEnd - iOpt;
        readList(opts->nChiMethod, iOpt+1, &opts->chiMethod, argv);
      }
      if (!strcmp(argv[optLocs[i]], "--metric")) {
        opts->metric = 1;
        opts->metricOrder = atoi(argv[iOpt+1]);
      }
      if (!strcmp(argv[optLocs[i]], "--poly-area")) {
        opts->polyA = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--thermo")) {
        opts->checkThermo = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--output-species-list")) {
        opts->listSpecies = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--FIRE-species-list")) {
        opts->fireSpecies = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--vtk-check")) {
        opts->checkVTK = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--add-wmol")) {
        opts->addWmol = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--ifile-info")) {
        opts->iFile = 1;
      }
      if (!strcmp(argv[optLocs[i]], "--test-parser")) {
        opts->testParse = 1;
      }
      break;
    }
    }
  } // end loop over options
  
  free(optLocs); free(opts->args);
}
/* -------------------------------------------------------------------------- */
