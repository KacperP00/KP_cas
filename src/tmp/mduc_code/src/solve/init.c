/*==========================================================*/
/*  This function is the first routine to be called         */
/*  this will initialise the basic routines and allocate    */
/*  memory to the variables needed later                    */
/*==========================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <dirent.h> 
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include "init.h"
#include "fileIO.h"
#include "ioStruct.h"
#include "transport.h"
#include "chemMech.h"
#include "chemInit.h"
#include "chemUtils.h"
#include "wiebe.h"
#include "accessChem.h"
#include "metric.h"
#include "parser.h"
#include "utils.h"
#include "solveUtils.h"
#include "engine.h"
#include "flameletInit.h"
#include "flamelet1D.h"
#include "flamelet2D.h"
#include "solve1D.h"
#include "solve2D.h"
#include "soot.h"
#include "sootHMOM.h"
#include "display.h"

void setDefaults(ioStruct *io, driveStruct *driver, timeStruct *timeVars,
                 flameStruct *flame, flowStruct *flow, engineStruct *engine)
{
  time_t rawtime;
  struct tm *timeinfo;
  
  // get time for run
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  
  /* ----- set critical parameters to -1 for warning later ----- */
  driver->ENGINE      = 0;
  driver->SDIAL       = 0;  // use CVODE as equation solver
  driver->PRES_MODE   = -1; // read pressure from file
  driver->MIX_MODE    = -1;
  driver->CHI_MODE    = 0;  // compute log form
  driver->CHI_INTERP  = 0;  // no interpolation
  driver->CHI_SCALE   = 0;  // no scaling
  driver->PDF_MODE    = 0;  // compute beta-pdf
  driver->TRANS_Z2    = 1;  // use standard transformation
  driver->PHI_HAT     = 0;  // no filtering of physical bounds
  driver->VISC_MODE   = 0;  // treat viscosity explicitly for soot

  /* ----- set driver defaults ----- */
  driver->TIME_MODE    = 0;       // time variables based on CFD solver
  driver->DIFF_MODE[0] = 1;       // Diffusion in dir1 is on
  driver->DIFF_MODE[1] = 0;       // Cross diff is off
  driver->DIFF_MODE[2] = 1;       // Diffusion in dir1 is on
  driver->CONV_MODE    = 1;       // convection is on
  driver->ENERGY       = 0;       // use temperature for energy equation
  driver->Q_LOSS       = 0;       // adiabatic
  driver->H_DIM        = 0;       // enthalpy is not a dimension
  driver->DIM          = 1;       // set direction to x1
  driver->ABSTOL       = 1.0e-9; // absolute tolerance
  driver->RELTOL       = 1.0e-9;  // relative tolerance
  driver->TABSTOL      = -1.0; // absolute tolerance for temperature
  driver->SABSTOL      = -1.0; // absolute tolerance for soot moments
  driver->QABSTOL      = -1.0; // absolute tolerance for quadrature
  driver->QRELTOL      = -1.0; // relative tolerance for quadrature
  driver->maxFlamelets = 5;       // maximum number of flamelets allowed
  driver->output_tsrc  = 0;       // don't output temp. sources by default
  driver->output_w     = 0;    // don't export reaction rates
  driver->radiation    = 0;    // don't consider radation heat loss
  driver->ignite    = 0;    // don't consider radation heat loss
  driver->thermophoresis = 0;  // dot't consider thermophoresis in soot model
  driver->sootDiff = 0;  // dot't consider diffusion in soot model
  driver->sootDiffusivity = 0.25;  // dot't consider diffusion in soot model
  driver->detailedDiff = 0;  // dot't consider detailed diffusio
  driver->readBody     = 1;       // read body of data files
  driver->DNS          = 0;       // Check if it is DNS

  /* ----- set input/output defaults ----- */
  io->debug        = 0;   // don't write debug comments
  io->readData     = 0;   // start from initialisation
  io->readFlow     = 0;   // don't read in flow data
  io->writeFlow    = 0;   // don't write flow data file (post-processing)
  io->freqType     = 0;
  io->writeFreq    = 1;
  io->printFreq    = 10;
  io->tFreq        = 1.0;
  io->tWrite       = 0.0;
  io->nMonitorY    = 0;   // don't print any species to the screen
  io->nFuel2       = 0;   // 1D flamelet
  io->initScalar   = 1;   // need to initialize scalars
  io->nOutputs     = 0;   // this will ensure default to output all species
  io->nOutputVar   = 0;   // no additional variables specified
  io->nInputs      = 0;
  io->nMonitorZ    = 10;  // monitor 10 grid points
  io->outputHtot   = 0;   // by default, output total enthalpy
  io->outputViscSrc= 1;   // by default, output viscous heat production
  io->outputGPDF   = 1;   // by default, output global pdf
  io->outputHR     = 0;   // by default, don't output hr term in flamelet file
  io->outputChi    = 1;   // by default, output scalar dissipation rate
  io->outputRho    = 0;   // by default, don't output Rho
  io->outputRhoSigma = 0;   // by default, don't output Rho
  io->outputCp     = 1;   // by default, output specific heat
  io->outputSumYi  = 1;
  io->outputMu     = 0;   // by default, output viscosity
  io->outputPDF    = 0;   // don't write pdf in flamelet (should be in flow)
  io->outputQdot   = 0;   // heat losses
  io->outputGamma  = 0;   // ratio of specific heats
  io->outputCFL    = 0;
  io->outputDryNOX = 0;   // by default, don't output dry NOx 6 percent
  io->outputConv   = 0;   // don't output conv. coeff. by default
  io->outputSootSrc= 0;
  io->outputSootMxmyp = 0;
  io->outputConvVel = 0;

  strftime(io->date, 100, "%c", timeinfo); // set date
  strcpy(io->sdialSolver, "cvode");
  strcpy(io->engineType, "HCCI");
  // set all boundary conditions to temperature
  io->bcType[0] = io->bcType[1] = io->bcType[2] = 'T';
  io->nOx = io->nFuel1 = 0;
  io->monitor = NULL;
  io->initSoot = NULL;

  // set all the pointer memory to NULL
  io->nDim = io->n1 = io->n2 = io->inputIndex = io->outputIndex = NULL;
  io->inputSpecies = io->outputSpecies = io->monitorSpecies = NULL;
  io->chemFuel = io->outputVar = NULL;
  io->initT = io->oxY = io->fuelY1 = io->fuelY2 = io->chemFuelY = NULL;
  io->gPDF = NULL;
  io->oxName = io->fuelName1 = io->fuelName2 = io->miscNames = NULL;
  io->oxIndex = io->fuelIndex1 = io->fuelIndex2 = io->monitorIndex = NULL;
  io->initY = io->miscIn = NULL;

  /* ----- set grid defaults ----- */
  io->nPivot1 = io->nPivot2 = 0;
  io->gPivot1 = io->gPivot2 = io->gPerc1 = io->gPerc2 = NULL;

  /* ----- set chi struct defaults -----*/
  io->chiData = NULL;
  io->nFiles = 0;

  /* ----- set default file names ----- */
  strcpy(io->dataDir, "");
  strcpy(io->inputFile,   "mduc.inp");
  strcpy(io->outputData,  "mduc_.rif");
  strcpy(io->monitorFile, "mduc_monitor");
  strcpy(io->errorFile,   "mduc.err");
  strcpy(io->flowData,    "mduc_.flo");
  strcpy(io->numStyle, "time-step"); // output file numbering style is flow time
  strcpy(io->ensightDir, "ensight");
  strcpy(io->postFile, "pdfs/"); // Default file prefix if needed
  strcpy(io->chiInterp, "none");
  strcpy(io->heatLoss, "none");
  strcpy(io->boundFilter, "none");
  io->fileNumPrec = 5;
  
  /* ----- set time defaults ----- */
  timeVars->t_m1    = timeVars->t0 = 0.0;
  timeVars->nt      = 0;
  timeVars->dt      = 1.0e-9;
  timeVars->dt_pre  = 0.0;
  timeVars->nIter   = -1;
  timeVars->t_end   = -1.0;
  timeVars->CFL     = 1.0;
  timeVars->RPM     = -1;
  timeVars->CAD_0   = 0.0;

  /* ----- set flamelet defaults ----- */
  flame->nFlamelets = 1;       // solve at least one flamelet
  flame->n1 = flame->n2 = 1;   // homogeneous solver
  flame->nDim       = 1;       // default to 1D solver
  flame->RGAS       = 8314.34; // universal gas constant J/kg-K
  flame->order      = 2;       // default to second order
  flame->intSrcQ = 0.0;

  /* ----- set engine defaults ----- */
  if (engine)
    defaultEngine(engine);

  // set all the pointer memory to NULL
  flame->sP = -1;
  flame->x1 = flame->x2 = flame->dx1 = flame->dx2 = NULL;
  flame->v1 = flame->v2 = flame->chiA = flame->chiB = flame->chiC = NULL;
  flame->vT = NULL;
  flame->vMxmyp = flame->vMxy = NULL;
  flame->Y = flame->Yp = flame->Yp_m1 = flame->Y_m1 = flame->Y_m2 = NULL;
  flame->Y_pre = flame->sumYi = flame->rhoDotSoot = flame->rhoDotSpray = NULL;
  flame->xi = flame->D = flame->rhoSigma = NULL; 
  flame->convCoeff = flame->chemSrcY = flame->chemSrcY_tmp = NULL;
  flame->rho = flame->cpmix = flame->sumCp = flame->sumCp2 = flame->T = NULL;
  flame->Mw = flame->R = flame->mu = flame->D = flame->gamma = flame->LeZ_LeI = NULL;
  flame->srcH = flame->srcP = flame->srcQ = flame->HRrate = flame->H = NULL;
  flame->chem_src = flame->diff_src = flame->conv_src = flame->visc_src = NULL;
  flame->chem_temp = flame->diff_temp = flame->conv_temp = flame->gPDF = NULL;
  flame->dY.x = flame->dY.y = flame->dY.xx = flame->dY.yy = flame->dY.xy = NULL;
  flame->dCp.x = flame->dCp.y = NULL;
  flame->dCp.xx = flame->dCp.yy = flame->dCp.xy = NULL;
  flame->grad.x = flame->grad.xx = flame->grad.y = flame->grad.yy = NULL;
  flame->nSootMom = 0;
  flame->sootMom = flame->sootMomSrc = flame->sootMomSrcPart = NULL;
  flame->sootMomxmyp = flame->sootMomTrans = NULL;
  flame->sootMomTransConv = flame->sootMomTransDiff = NULL;

  flame->radiation = NULL;

  /* ----- set flow defaults ----- */
  flow->P = flow->P_m1 = flow->P_m2 = 0.0;
  flow->phi1min = flow->phi2min = flow->phi1min_m1 = flow->phi2min_m1 = 0.0;
  flow->phi1max = flow->phi2max = flow->phi1max_m1 = flow->phi2max_m1 = 1.0;
  flow->phi1min_corr = flow->phi2min_corr = 0.0;
  flow->phi1max_corr = flow->phi2max_corr = 1.0;
  flow->ZRef[0]    = flow->ZRef[1]    = flow->ZRef[2]    = -1.0;
  flow->chiRef[0]  = flow->chiRef[1]  = flow->chiRef[2]  = -1.0;
  flow->strain[0]  = flow->strain[1]  = flow->strain[2]  = -1.0;
  flow->phiMean[0] = flow->phiMean[1] = flow->phiMean[2] = -1.0;
  flow->phiVar[0]  = flow->phiVar[1]  = flow->phiVar[2]  = -1.0;
  flow->meanRho = 1.0;
  flow->chi1 = flow->chi2 = flow->chi12 = NULL;
  flow->chi1_m1 = flow->chi2_m1 = flow->chi12_m1 = NULL;
  flow->avgSdot = flow->avgSdot_m1 = flow->pdf = flow->avgH = NULL;
  flow->avgSdot_M = flow->strainCondZ = NULL;
  flow->phiStar = NULL;
  flow->phi = flow->Dt_phi = flow->phiRaw = flow->phiHat = NULL;
  flow->Dt_phiRaw = flow->Dt_phiHat = flow->DtHat_phi = NULL;
  flow->fhat.w = NULL;

  /* ----- set statistics defaults ----- */
  io->mean[0] = io->mean[1] = io->mean[2] = -1;
  io->var[0]  = io->var[1]  = io->var[2]  = io->var[3] = -1;
}
/* -------------------------------------------------------------------------- */

void initIO(ioStruct *io, driveStruct *driver, chemStruct *chem)
{
  int i;
  char *dir, rootDir[200], buffer[200];
  struct stat st;

  /* set up the output directory */
  if (strlen(io->dataDir) > 0) {
    strcpy(rootDir, io->dataDir);
    dir = strtok(rootDir, "/");
    strcpy(buffer, dir);
    if (io->dataDir[0] == '/')
      prefixString("/", buffer);
    while (dir != NULL) {
      if (stat(buffer, &st) != 0)
        mkdir(buffer, 0755);
      dir = strtok(NULL, "/");
      if (dir != NULL) {
        strcat(buffer, "/");
        strcat(buffer, dir);
      }
    }
    strcat(io->dataDir, "/");
  }
  
  /* add directory to output data file names */
  prefixString(io->dataDir, io->outputData);
  if (io->writeFlow) prefixString(io->dataDir, io->flowData);

  /* add directory to data files to read (if absolute path not given) */
  if (io->readData && (io->dataDir[0] != '/'))
    prefixString(io->dataDir, io->inputData);
  if (io->readFlow && (io->dataDir[0] != '/'))
    prefixString(io->dataDir, io->inputFlow);

  /* construct scalars file name */
  strcpy(buffer, io->outputData);
  dir = strtok(buffer, ".");
  strcpy(io->scalars, strcat(dir, "scalars"));

  /* set equation solver */
  if (driver->ENERGY) {
    driver->SDIAL = 1;   /* must use DAE solver */
    io->outputHtot = 1;  /* ensure the solution is output */
  }
  else if (!strcmp(io->sdialSolver, "cvode"))
    driver->SDIAL = 0;
  else if (!strcmp(io->sdialSolver, "ida"))
    driver->SDIAL = 1;
  else if (!strcmp(io->sdialSolver, "rk"))
    driver->SDIAL = -1;

  /* set up the heat loss model */
  if ( !strcmp(io->heatLoss, "woschni") )
    driver->Q_LOSS = 1;

  /* set problem types */
  if ( !strcmp(io->pressureMode, "isobaric") ||
       !strcmp(io->pressureMode, "constant-pressure") )
    driver->PRES_MODE = 0;
  else if ( !strcmp(io->pressureMode, "isochoric") ||
            !strcmp(io->pressureMode, "constant-volume") )
    driver->PRES_MODE = 1;
  else if ( !strcmp(io->pressureMode, "compression") )
    driver->PRES_MODE = 2;
  else if ( !strcmp(io->pressureMode, "implicit") )
    driver->PRES_MODE = 3;
  else if ( !strcmp(io->pressureMode, "explicit") )
    driver->PRES_MODE = 4;
  else if  (!strcmp(io->pressureMode, "from-file"))
    driver->PRES_MODE = -1;

  /* set pdf treatment */
  if (io->pdfMode==PDF_TYPE_BETA)
    driver->PDF_MODE = 0;
  else if (io->pdfMode==PDF_TYPE_FROM_FILE)
    driver->PDF_MODE = -1;

  /* set scalar dissipation functional form */
  if (!strncmp(io->chiMode, "log", 3))
    driver->CHI_MODE = 0;
  else if (!strcmp(io->chiMode, "erfc"))
    driver->CHI_MODE = 1;
  else if (!strcmp(io->chiMode, "equation") || !strcmp(io->chiMode, "solve"))
    driver->CHI_MODE = 2;
  else if (!strcmp(io->chiMode, "from-flamelet"))
    driver->CHI_MODE = 3;
  else if (!strcmp(io->chiMode, "from-file"))
    driver->CHI_MODE = -1;

  /* set chi interpolation treatment */
  if (!strncmp(io->chiInterp, "linear", 3))
    driver->CHI_INTERP = 1;
  else if (!strncmp(io->chiInterp, "least-squares", 5))
    driver->CHI_INTERP = 2;
  else if (!strncmp(io->chiInterp, "laplace", 5))
    driver->CHI_INTERP = 3;
  else if (!strcmp(io->chiInterp, "none"))
    driver->CHI_INTERP = 0;

  /* default to reading from file list if any variable is specified in file */
  if (driver->PRES_MODE < 0 || driver->PDF_MODE < 0 || driver->CHI_MODE < 0)
    driver->TIME_MODE = 1;

  /* set time series mode for post-processing */
  if (!strcmp(io->timeMode, "steady"))
    driver->TIME_MODE = -1;
  else if (!strcmp(io->timeMode, "input"))
    driver->TIME_MODE = 0;
  else if (!strncmp(io->timeMode, "file-list", 5))
    driver->TIME_MODE = 1;
  else if (!strncmp(io->timeMode, "scalar-file", 6))
    driver->TIME_MODE = 2;
  else if (!strcmp(io->timeMode, "engine")) {
    driver->TIME_MODE = 3;
    if (driver->PRES_MODE > 2) {
      if (!strcmp(io->engineType,"HCCI") || !strcmp(io->engineType,"hcci"))
        driver->ENGINE  = 1;
      else if (!strcmp(io->engineType,"SI") || !strcmp(io->engineType,"si"))
        driver->ENGINE  = 2;
      else if (!strcmp(io->engineType,"SACI") || !strcmp(io->engineType,"saci"))
        driver->ENGINE  = 3;
      fprintf(stderr, "driver->ENGINE:    %i\n", driver->ENGINE);
      fprintf(stderr, "driver->PRES_MODE: %i\n", driver->PRES_MODE);
    }
  }

  /* set up the output frequency */
  char freq[32];
  paramGetString(freq,"frequency-type","step",io->parser);
  if (!strcmp(freq, "step"))
    io->freqType = 0;
  else if (!strcmp(freq, "time"))
    io->freqType = 1;
  else if (!strcmp(freq, "cad"))
    io->freqType = 2;
  else if (isdigit(freq[0]))
    io->freqType = atoi(freq);
  else
    fprintf(stderr, "Warning: frequency-type %s not recognized\n",freq);

  /* get the frequency value */
  if (io->freqType == 0)
    io->writeFreq = paramGetInt("write-frequency",1,io->parser);
  else if (io->freqType == 1)
    io->tFreq = paramGetReal("write-frequency",1e-3,io->parser);

  /* determine the viscosity treatment (mainly for soot) */
  if (!strcmp(io->viscUpdate, "implicit") || !strcmp(io->viscUpdate, "1"))
    driver->VISC_MODE = 1;

  /* set the transformation for two mixture fraction formulation */
  if (!strcmp(io->transZ2, "none"))
    driver->TRANS_Z2 = 0;
  else if (!strcmp(io->transZ2, "doran") || !strcmp(io->transZ2, "standard"))
    driver->TRANS_Z2 = 1;
  else if (!strcmp(io->transZ2, "hasse"))
    driver->TRANS_Z2 = 2;
  if (driver->H_DIM != 0)
    driver->TRANS_Z2 = 0;

  /* set the transformation for two mixture fraction formulation */
  if (!strcmp(io->boundFilter, "none"))
    driver->PHI_HAT = 0;
  else if (!strcmp(io->boundFilter, "value"))
    driver->PHI_HAT = 1;
  else if (!strcmp(io->boundFilter, "derivative"))
    driver->PHI_HAT = 2;

  io->oxIndex    = (int*) calloc(io->nOx, sizeof(int));
  io->fuelIndex1 = (int*) calloc(io->nFuel1, sizeof(int));
  
  /* get the indices from the mechanism for the fuel, oxidizer, and outputs */
  for (i=0; i<io->nOx; ++i)
    io->oxIndex[i] = getSpeciesIndex(chem->nSpecies, chem->species, 
                                     io->oxName[i]);
  for (i=0; i<io->nFuel1; ++i)
    io->fuelIndex1[i] = getSpeciesIndex(chem->nSpecies, chem->species,
                                        io->fuelName1[i]);

  if (io->nFuel2 > 0) {
    io->fuelIndex2 = (int*) calloc(io->nFuel2, sizeof(int));
    for (i=0; i<io->nFuel2; ++i)
      io->fuelIndex2[i] = getSpeciesIndex(chem->nSpecies, chem->species, 
                                          io->fuelName2[i]);
  }
  
  /* assign indices for the species to be output */
  /* if no species were specified explicitly in the input file, default all */
  if (io->nOutputs == 0) {
    io->nOutputs = chem->nSpecSolve;
    io->outputSpecies = (char**) calloc(io->nOutputs, sizeof(char*));
    io->outputIndex   = (int*) calloc(io->nOutputs, sizeof(int));
    for (i=0; i<io->nOutputs; ++i) {
      io->outputSpecies[i] = (char*) calloc(25, sizeof(char));
      strcpy(io->outputSpecies[i], chem->species[i]);
      io->outputIndex[i] = i;
    }
  }
  else {
    io->outputIndex = (int*) calloc(io->nOutputs, sizeof(int));
    for (i=0; i<io->nOutputs; ++i)
      io->outputIndex[i] = getSpeciesIndex(chem->nSpecies, chem->species, 
                                           io->outputSpecies[i]);
  }

  /* check whether fuel species still need to be initialized */
  if (chem->nFuel == 0) {
    /* TODO: check whether a specific boundary was specified, for now assume
       it is Z1 boundary */
    chem->nFuel = io->nFuel1;
    chem->fuel  = (char**)  calloc(chem->nFuel, sizeof(char*));
    chem->fuelY = (double*) calloc(chem->nFuel, sizeof(double));
    for (i=0; i<chem->nFuel; i++) {
      chem->fuel[i] = (char*) calloc(32, sizeof(char));
      strcpy(chem->fuel[i], io->fuelName1[i]);
      chem->fuelY[i] = io->fuelY1[i];
    }
    initChemFuel(chem);
  }
  
  /* check if other variables were specified for output */
  for (i=0; i<io->nOutputVar; i++) {
    if (!strcmp(io->outputVar[i], "heat-release"))
      io->outputHR = 1;
    else if (!strcmp(io->outputVar[i], "density") ||
             !strcmp(io->outputVar[i], "rho"))
      io->outputRho = 1;
    else if (!strcmp(io->outputVar[i], "rhoSigma"))
      io->outputRhoSigma = 1;
    else if (!strncmp(io->outputVar[i], "viscosity", 4))
      io->outputMu = 1;
    else if (!strcmp(io->outputVar[i], "dry-nox"))
      io->outputDryNOX = 1;
    else if (!strcmp(io->outputVar[i], "scalar-dissipation") ||
             !strcmp(io->outputVar[i], "chi"))
      io->outputChi = 1;
    else if (!strcmp(io->outputVar[i], "specific-heat") ||
             !strcmp(io->outputVar[i], "cp"))
      io->outputCp = 1;
    else if (!strcmp(io->outputVar[i], "heat-loss") ||
             !strcmp(io->outputVar[i], "qdot"))
      io->outputQdot = 1;
    else if (!strcmp(io->outputVar[i], "gamma"))
      io->outputGamma = 1;
    else if (!strncmp(io->outputVar[i], "convection-coefficient", 16))
      io->outputConv = 1;
  }
  
  /* either open a monitor file or print to stderr */
  if (io->dataDir[0] != '/' && strlen(io->dataDir) > 0)
    prefixString(io->dataDir, io->monitorFile);

  /* assign the indices of species to monitor on the screen */
  /* if no species are specified, default to none */
  if (io->nMonitorY > 0) {
    io->monitorIndex = (int*) calloc(io->nMonitorY, sizeof(int));
    for (i=0; i<io->nMonitorY; ++i)
      io->monitorIndex[i] = getSpeciesIndex(chem->nSpecies, chem->species,
                                            io->monitorSpecies[i]);
  }

  /* if read from a flamelet, get the indices for the species read */
  if (io->readData) {
    io->inputIndex = (int*) calloc(io->nInputs, sizeof(int));
    for (i=0; i<io->nInputs; ++i)
      io->inputIndex[i] = getSpeciesIndex(chem->nSpecies, chem->species, 
                                          io->inputSpecies[i]);
  }
 
  /* set chi struct required for stand alone calculations from chi data*/
  if (driver->CHI_MODE == 3) {
    int counter = 0, i= 0, j=0;
    int *tInd;
    double *tSorted;
    char **fileNames;
    DIR *d;
    struct dirent *dir;
    chiStruct *chiDataOld;
	
    d = opendir(io->chiDataDir);
    if (d)
    {
      // count number of files in directory
      // Note: number of files inculdes current and previous directory
      while((dir = readdir(d)) != NULL)
      {
	if (strcmp(dir->d_name,".") != 0 && strcmp(dir->d_name,"..") != 0) {
	  counter++;
        }
      }
      closedir(d);
    }
    io->nFiles = counter;

    // do memory allocation
    io->chiData = (chiStruct*) calloc(io->nFiles, sizeof(chiStruct));
    chiDataOld = (chiStruct*) calloc(io->nFiles, sizeof(chiStruct));
    fileNames = (char**) calloc(io->nFiles, sizeof(char*));
    tSorted = (double*) calloc(io->nFiles, sizeof(double));
    tInd = (int*) calloc(io->nFiles, sizeof(int));

    d = opendir(io->chiDataDir);
    if (d)
    {

      // get filenames
      while ((dir = readdir(d)) != NULL)
      {
	if (strcmp(dir->d_name,".") != 0 && strcmp(dir->d_name,"..") != 0) {
 	  fileNames[i] = (char*) calloc(256, sizeof(char));
          strcpy(fileNames[i], io->chiDataDir);
	  strcat(fileNames[i], dir->d_name);
	  i++;
        }
      }

      closedir(d);
    } else {
       fprintf(stderr, "ERROR: Reading chi data files failed");
       exit(2);
    }

    // read data from files 
    for(i=0; i<io->nFiles; i++) {
       /* initialize structs */
      io->chiData[i].tChi = 0.0;
      io->chiData[i].chiVal = NULL;

      /* read time and chi data */
      readChi(fileNames[i], io, i);
 
      /* tSorted is required for sorting data*/
      tSorted[i] = io->chiData[i].tChi;
      chiDataOld[i] = io->chiData[i];
    }

    // sort data with respect of time
    qsort(tSorted, io->nFiles, sizeof(double), cmpfunc);
    
    /* geting indiced */
    for (i=0; i<io->nFiles; i++) {
      for (j=0; j<io->nFiles; j++) {
 	if (tSorted[i] == io->chiData[j].tChi)
	  tInd[i] = j;
      }
    }
    
    for (i=0; i<io->nFiles; i++) {
      io->chiData[i] = chiDataOld[tInd[i]];
    } 

    // Gitter interpolation, input gitter nicht vorhanden?
  }
}

/* -------------------------------------------------------------------------- */
/* compare function required for sorting chi time */

int cmpfunc(const void *x, const void *y) {
  double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
}

/*----------------------------------------------------------------------------*/

void initGrid(int n, double xmin, double xmax, int nP, double *pivot,
              double *perc, double **x, double **dx)
{
  int i, j, j0, *nPart, jEnd;
  double dx0, *r, *xEnd;

  /* error checking */
  if ( (pivot && !perc) || (!pivot && perc) ) {
    fprintf(stderr, "Warning: insufficient information about grid pivots\n");
    fprintf(stderr, "         using equidistant grid\n");
    nP = 0;
  }

  /* allocate memory for grid parameters */
  x[0]  = (double*) calloc(n, sizeof(double));
  dx[0] = (double*) calloc(n-1, sizeof(double));
  nPart = (int*)    calloc(nP+1, sizeof(int));
  r     = (double*) calloc(nP+1, sizeof(double));
  xEnd  = (double*) calloc(nP+1, sizeof(double));

  j0 = 0;
  for (i=0; i<nP; i++) {
    nPart[i] = floor(perc[i]*n);
    xEnd[i]  = pivot[i];
    j0 += nPart[i];
  }
  nPart[nP] = n - j0;
  xEnd[nP]  = xmax;

  /* compute the growth rates for each section */
  r[0] = 1.0;
  dx0 = (xEnd[0]-xmin)/(nPart[0]-1);
  for (i=1; i<=nP; i++) {
    r[i] = calcGridGrowth(nPart[i], dx0, xEnd[i]-xEnd[i-1]);
    dx0 *= pow(r[i],nPart[i]);
  }
  
  /* fill the grid */
  x[0][0] = xmin;
  j0 = 1;
  jEnd = nPart[0];
  dx0  = (xEnd[0]-xmin)/(nPart[0]-1);
  for (i=0; i<=nP; i++) {
    if (i==0)
      jEnd = nPart[0];
    else
      jEnd = j0+nPart[i];
    
    for (j=j0; j<jEnd; j++) {
      x[0][j]    = x[0][j-1] + dx0;
      dx[0][j-1] = dx0;
      dx0 *= r[i];
    } /* end loop over mesh section */
    j0 += nPart[i];
    if (i==0) j0--;
  } /* end loop over sections */
  
  free(xEnd); free(r); free(nPart);
}
/* -------------------------------------------------------------------------- */

void initSolver(ioStruct *io, driveStruct *driver, chemStruct *chem,
                flameStruct *flame, flowStruct *flow)
{
  int i, j, k, k_pre, r, row, col, ii, maxDim, nGrid, nEq, iX1;
  int found,f;
  double oppSideY, *initH, H_ref, *energy;
  const double BOLTZMANN=5.67051E-8; // [W/(m^2 K^4)]
  
  /* read in some input parameters */
  driver->PHYS_COORD = paramGetInt("physical-coordinate",0,io->parser);

  /* ----- initialise parameters from chemical mechanism ----- */
  flame->nVars = chem->nSpecSolve+1;
  flame->sT    = chem->nSpecSolve;
  if (driver->CHI_MODE == 2) { /* solve equation for chi */
    flame->sChi   = flame->nVars;
    flame->nVars += 1;
  }
  if (chem->sootMode == 1) { /* solve equations for soot moments */
    flame->sSoot  = flame->nVars;
    flame->nVars += chem->nSootMom;
//    fprintf(stderr, "adding soot moments to solution: %i\n", chem->nSootMom);
  }

  /* assign some indices for specific species */
  flame->sO2   = getSpeciesIndex(chem->nSpecies, chem->species, "O2");
  flame->sN2   = getSpeciesIndex(chem->nSpecies, chem->species, "N2");

  maxDim = IMAX(flame->n1, flame->n2);
  nGrid  = flame->n1*flame->n2;
  nEq    = nGrid*flame->nVars;
  flame->j0 = 0;
  flame->jN = flame->n1-1;
  flame->k0 = 0;
  flame->kN = flame->n2-1;

  /* solve additional equations for engine */
  if (driver->ENGINE) {
    if (driver->PRES_MODE == 3) { /* implicit pressure */
      flame->sP  = nEq;
      nEq       += 1;
    }
    flame->sM  = nEq;
    flame->sMa = nEq; nEq += 1;   /* air mass */
    flame->sMf = nEq; nEq += 1;   /* fuel mass */
    flame->sMe = nEq; nEq += 1;   /* egr mass */
    if (driver->ENGINE > 1) { /* SI */
      flame->sMb = nEq; nEq += 1; /* burned mass */
      flame->sTb = nEq; nEq += 1; /* burned zone temperature */
      flame->sDF = nEq; nEq += 1; /* flame diameter */
    }
  }

  /* ----- allocate solution ----- */
  flame->Y     = (double*) calloc(nEq, sizeof(double));
  flame->Yp    = (double*) calloc(nEq, sizeof(double));
  flame->Y_m1  = (double*) calloc(nEq, sizeof(double));
  flame->Y_m2  = (double*) calloc(nEq, sizeof(double));
  flame->Y_pre = (double*) calloc(nEq, sizeof(double));
  flame->sumYi = (double*) calloc(nEq, sizeof(double));
  flame->Yp_m1 = (double*) calloc(nEq, sizeof(double));
  flame->chemSrcY     = (double*) calloc(nEq, sizeof(double));
  flame->chemSrcY_tmp = (double*) calloc(nEq, sizeof(double));

  if (flame->sP >= 0)
    flame->Y[flame->sP] = flow->P;

  /* ----- allocate local properties ----- */
  flame->chiA        = (double*) calloc(nGrid, sizeof(double));
  flame->rhoDotSoot  = (double*) calloc(nGrid, sizeof(double));
  flame->rhoDotSpray = (double*) calloc(nGrid, sizeof(double));
  flame->xi          = (double*) calloc(nGrid, sizeof(double));
  flame->D           = (double*) calloc(nGrid, sizeof(double));
  flame->rhoSigma    = (double*) calloc(nGrid, sizeof(double));
  flame->v1          = (double*) calloc(nGrid, sizeof(double));
  flame->vT          = (double*) calloc(nGrid, sizeof(double));
  if (chem->soot) {
     flame->vMxmyp          = (double*) calloc(nGrid, sizeof(double));
     flame->vMxy          = (double*) calloc(nGrid, sizeof(double));
  }
  if (flame->rho == NULL)
    flame->rho  = (double*) calloc(nGrid, sizeof(double));

  flame->gamma  = (double*) calloc(nGrid, sizeof(double));
  flame->cpmix  = (double*) calloc(nGrid, sizeof(double));
  flame->Mw     = (double*) calloc(nGrid, sizeof(double));
  flame->R      = (double*) calloc(nGrid, sizeof(double));
  flame->mu     = (double*) calloc(nGrid, sizeof(double));
  if (chem->soot) 
     flame->D   = (double*) calloc(nGrid, sizeof(double));

  if (driver->detailedDiff)
     flame->LeZ_LeI = (double*) calloc(nGrid*chem->nSpecSolve, sizeof(double));

  flame->srcP   = (double*) calloc(nGrid, sizeof(double));
  flame->srcH   = (double*) calloc(nGrid, sizeof(double));
  flame->srcQ   = (double*) calloc(nGrid, sizeof(double));
  flame->sumCp  = (double*) calloc(nGrid, sizeof(double));
  flame->HRrate = (double*) calloc(nGrid, sizeof(double));
  flame->visc_src = (double*) calloc(nGrid, sizeof(double));
  flame->radiation = (double*) calloc(nGrid, sizeof(double));
  flame->gPDF     = (double*) calloc(nGrid, sizeof(double));

  if (driver->ENERGY)
    flame->T = (double*) calloc(nGrid, sizeof(double));

  flame->convCoeff = (double*) calloc(nGrid, sizeof(double));
  if (driver->output_tsrc) {
    flame->conv_src  = (double*) calloc(nGrid, sizeof(double));
    flame->chem_src  = (double*) calloc(nGrid, sizeof(double));
    flame->diff_src  = (double*) calloc(nGrid, sizeof(double));
    flame->conv_temp = (double*) calloc(nGrid, sizeof(double));
    flame->chem_temp = (double*) calloc(nGrid, sizeof(double));
    flame->diff_temp = (double*) calloc(nGrid, sizeof(double));
  }

  if (chem->soot) {
    driver->sootMode = chem->sootMode;
    if (chem->sootMode <= 3) {
       flame->nSootMom = chem->nSootMom;
       flame->sootMom    = (double*) calloc(flame->nSootMom*nGrid, sizeof(double));
       flame->sootMomxmyp = (double*) calloc((flame->nSootMom-1)*nGrid, sizeof(double));
       flame->sootMomSrc = (double*) calloc(flame->nSootMom*nGrid, sizeof(double));
       flame->sootMomSrcPart =
         (double*) calloc(chem->nSootSrc*flame->nSootMom*nGrid, sizeof(double));
       flame->sootMomTrans = (double*) calloc(flame->nSootMom*nGrid, sizeof(double));
       flame->sootMomTransConv = (double*) calloc(flame->nSootMom*nGrid, sizeof(double));
       flame->sootMomTransDiff = (double*) calloc(flame->nSootMom*nGrid, sizeof(double));
    }
    if (chem->sootMode >= 3) {
       flame->sqrtT    		= (double*) calloc(nGrid, sizeof(double));
       flame->T_mu    		= (double*) calloc(nGrid, sizeof(double));
       flame->muSqrtMw_rhoSqrtT = (double*) calloc(nGrid, sizeof(double));
       flame->sgCoeff    	= (double*) calloc(nGrid, sizeof(double));
       flame->oxCoeff  		= (double*) calloc(nGrid, sizeof(double));
       flame->o2Coeff  		= (double*) calloc(nGrid, sizeof(double));
       flame->dimProdRate  	= (double*) calloc(nGrid, sizeof(double));
    }
  }

  // if flamelet file was read, place data in appropriate format
  if (io->readData) {
    // compute dx according to grid read in
    if (flame->n1 > 1) {
      flame->dx1 = (double*) calloc(flame->n1-1, sizeof(double));
      for (j=0; j<flame->n1-1; j++)
        flame->dx1[j] = flame->x1[j+1] - flame->x1[j];
    }
    
    // set the whole solution vector to zero
    for (i=0; i<flame->nVars*flame->n1*flame->n2; ++i)
      flame->Y[i] = 0.0;
    
    // set the whole solution vector to zero
    for (i=0; i<flame->nVars*flame->n1*flame->n2; ++i) {
      flame->chemSrcY[i] = 0.0;
      flame->chemSrcY_tmp[i] = 0.0;
    }
    
    if (driver->ENERGY) {
      if (!flame->H) { // file only had temperature
        // compute enthalpy
        flame->H = (double*) calloc(nGrid, sizeof(double));
        
        double *cp = (double*) calloc(chem->nSpecies, sizeof(double));
        double *h  = (double*) calloc(chem->nSpecies, sizeof(double));
        
        for (k=0; k<flame->n2; k++) {
          r   = k*flame->n1;
          for (j=0; j<flame->n1; j++) {
            chemComputeThermoData(h,cp,io->initT[r+j],chem);
            flame->H[r+j] = 0.0;
            for (i=0; i<chem->nSpecSolve; i++)
              flame->H[r+j] += h[i]*io->initY[i][r+j];
          } // end loop over columns
        } // end loop over rows
        
        free(cp);
        free(h);
      }
      for (i=0; i<nGrid; i++) flame->T[i] = io->initT[i];
      energy = flame->H;
    }
    else
      energy = io->initT;

    if (!flame->H)  // file only had temperature
      flame->H = (double*) calloc(nGrid, sizeof(double));

    // initialise the solution vector based on the species
    // read from the input file
    for (k=0; k<flame->n2; ++k) {
      r   = k*flame->n1;
      row = k*flame->n1*flame->nVars;
      for (j=0; j<flame->n1; ++j) {
        col = j*flame->nVars;
        for (i=0; i<io->nInputs; ++i)
	  if (io->inputIndex[i] >= 0)
#ifdef SUNDIALS_FROM_3_2
	     // TODO: Prohibit negative values if contraints are enforced
	     // Negative initial values mostly likely caused by importing 
	     // text file rather that binary
             flame->Y[row+col+io->inputIndex[i]] = fmax(io->initY[i][r+j],0.0);
#else
             flame->Y[row+col+io->inputIndex[i]] = io->initY[i][r+j];
#endif 
        flame->Y[row+col+flame->sT] = io->initT[r+j];
        //flame->T[row+col] = energy[r+j];//io->initT[r+j];
      } // end loop over columns
    } // end loop over rows

    if (chem->soot) {
      if (chem->sootMode <= 3 && io->initSootMom==0) {
         for (k=0; k<flame->n2; k++) {
           r   = k*flame->n1;
           row = k*flame->n1*flame->nSootMom;
           for (j=0; j<flame->n1; j++) {
             col = j*flame->nSootMom;
             flame->sootMom[row+col+0] = io->initSoot[0][r+j];
             flame->sootMom[row+col+1] = io->initSoot[1][r+j];
             flame->sootMom[row+col+2] = io->initSoot[2][r+j];
             if (flame->nSootMom > 4) {
               flame->sootMom[row+col+3] = io->initSoot[3][r+j];
               flame->sootMom[row+col+4] = io->initSoot[4][r+j];
               flame->sootMom[row+col+5] = io->initSoot[5][r+j];
               flame->sootMom[row+col+6] = io->initSoot[6][r+j];
             }
             else
               flame->sootMom[row+col+3] = io->initSoot[6][r+j];

             /* Make sure that the minimum values are not smaller than the
 	      * the solver tolerance */ 
             //reinitSootMoments(flame->nSootMom,&flame->sootMom[row+col],driver->SABSTOL);
             reinitSootMoments(flame->nSootMom,&flame->sootMom[row+col],1.0E-20);

             if (chem->sootMode == 1) {
               for (i=0; i<flame->nSootMom; i++) {
                 flame->Y[(k*flame->n1+j)*flame->nVars+flame->sSoot+i] = 
                   flame->sootMom[row+col+i];
               } /* end loop over soot moments */
             } /* implicit soot mode */

           } /* end loop over columns */
         } /* end loop over rows */
      } /* soot mode */
      else if (chem->sootMode <= 3) {
        for (k=0; k<flame->n2; k++) {
          for (j=0; j<flame->n1; j++) {
            initSootMoments(chem->nSootMom,
                            &flame->sootMom[(k*flame->n1+j)*flame->nSootMom],
			    1.0E-20);
        //    initSootMoments(chem->nSootMom,
        //                    &flame->sootMom[(k*flame->n1+j)*flame->nSootMom],
	//		    driver->SABSTOL);

            if (chem->sootMode == 1) {
              /* copy initial soot moments into solver variables */
              for (i=0; i<chem->nSootMom; i++)
                flame->Y[(k*flame->n1+j)*flame->nVars+flame->sSoot+i] = 
                  flame->sootMom[(k*flame->n1+j)*flame->nSootMom+i];
            }
          } /* end loop over columns */
        } /* end loop over rows */
      }
      for (k=0; k<flame->n2; k++) {
        r   = k*flame->n1;
        row = k*flame->n1*flame->nSootMom;
        for (j=0; j<flame->n1; j++) {
          col = j*flame->nSootMom;
          /* make sure viscosity is up-to-date */
          if (chem->readTrans)
             flame->mu[r+j] = calcMixViscosity(chem->nSpecies,
                                            flame->Y[j*flame->nVars+flame->sT],
                                            &flame->Y[j*flame->nVars], chem);
        } /* end loop over columns */
      } /* end loop over rows */
    } /* end soot */
    
    flame->nSootMom = chem->nSootMom;

    if (flame->rho == NULL) 
      flame->rho = (double*) calloc(nGrid, sizeof(double));
    
    // allocate memory specific to 2D flamelet
    if (flame->nDim == 2) {
      // compute dx according to grid read in
      flame->dx2 = (double*) calloc(flame->n2-1, sizeof(double));
      for (j=0; j<flame->n2-1; j++)
        flame->dx2[j] = flame->x2[j+1] - flame->x2[j];
      
      flame->chiB    = (double*) calloc(nGrid, sizeof(double));
      flame->chiC    = (double*) calloc(nGrid, sizeof(double));
      flame->v2      = (double*) calloc(nGrid, sizeof(double));
      flame->sumCp2  = (double*) calloc(nGrid, sizeof(double));
      flame->dY.y    = (double*) calloc(flame->nVars*maxDim,sizeof(double));
      flame->dY.yy   = (double*) calloc(flame->nVars*maxDim,sizeof(double));
      flame->dY.xy   = (double*) calloc(flame->nVars*maxDim,sizeof(double));
      flame->dCp.y   = (double*) calloc(maxDim, sizeof(double));
    }
    
    /* initialise the metric quantities */
    if (flame->nDim > 0)
      initMetric(flame);
    
    /* compute the current properties */
    double wmix;
    for (k=0; k<flame->n2; k++) {
      row = k*flame->n1*flame->nVars;
      r   = k*flame->n1;
      for (j=0; j<flame->n1; j++) {
        col  = j*flame->nVars;
        wmix = calcMixtureW(chem->nSpecSolve, &flame->Y[row+col], chem->W);
        flame->Mw[r+j] = wmix;
        flame->rho[r+j] = wmix*flow->P/(flame->RGAS*flame->Y[row+col+flame->sT]);
        if (chem->readTrans)
          flame->mu[r+j]  = calcMixViscosity(chem->nSpecies,
                                             flame->Y[row+col+flame->sT],
                                             &flame->Y[row+col], chem);
      } /* end loop over columns */
    } /* end loop over rows */

    if (chem->soot) {
       if (chem->sootMode >= 3) {
          for (k=0; k<flame->n2; k++) {
            row = k*flame->n1*flame->nVars;
            r   = k*flame->n1;
            for (j=0; j<flame->n1; j++) {
              col  = j*flame->nVars;
              /*TODO: Compute concentration only for PAH */
              /* compute and store concentrations of all species at local point */
              for (i = 0; i < chem->nSpecies; i++)
                chem->c[i] = flame->rho[r+j]*flame->Y[row+col+i]/chem->W[i];
              /* copy the species necessary for soot reactions */
              int iS;
              for (iS=0; iS<chem->nSootSp; iS++)
             chem->sootRxnConc[iS] = chem->c[chem->iSootSp[iS]];

	      double sgCoeff, oxCoeff,o2Coeff, dimProdRate;
              sootHMOM_GasPhase_Coeff(flame->Y[row+col+flame->sT], flame->mu[r+j],
              	                flame->rho[r+j], flame->Mw[r+j], chem->sootRxnConc, 
              			&sgCoeff, &oxCoeff,&o2Coeff, &dimProdRate);

	      flame->sgCoeff[r+j] 	    = sgCoeff;
	      flame->oxCoeff[r+j]  	    = oxCoeff;
	      flame->o2Coeff[r+j] 	    = o2Coeff;
	      flame->dimProdRate[r+j] 	    = dimProdRate;
              flame->sqrtT[r+j] 	    = sqrt(flame->Y[row+col+flame->sT]);
              flame->T_mu[r+j]  	    = flame->Y[row+col+flame->sT]/flame->mu[r+j];
              flame->muSqrtMw_rhoSqrtT[r+j] = flame->mu[r+j]*sqrt(flame->Mw[r+j]/
              			(1000.0*flame->Y[row+col+flame->sT]*flame->rho[r+j]*flame->rho[r+j]));
            } 
          }
       }
    }
    if (driver->output_w) {
      if (io->gPDF != NULL) 
      for (i=0; i<nGrid; i++) flame->gPDF[i] = io->gPDF[i];
      else
         for (i=0; i<nGrid; i++) flame->gPDF[i] = 1.0;
    }
    if (driver->CHI_MODE == 3) {
    for (f=0; f < io->nFiles; f++) {/* loop over files */
      /* allocate memory */ 
      io->chiData[f].chiVal = (double*) calloc(flame->n1, sizeof(double));
      
      /* set first and last value to zero */
      io->chiData[f].chiVal[0] = 0.0;
      io->chiData[f].chiVal[flame->n1-1] = 0.0;
      for (j=1; j < flame->n1-1; j++)
         io->chiData[f].chiVal[j] = io->chiData[f].chiRaw[j];
    }
    }
    if (driver->radiation) {
       for (j=0; j<flame->n1; j++) {
       k   = j*flame->nVars + flame->sT; /* array index temperature */
       /* update radiation source term */
       /* Gas phase radiation (only CO2 / H2O / CH4 / CO considered) */
       // Grosshandler, W. L., RADCAL: A Narrow-Band Model for Radiation Calculations in a Combustion Environment, NIST technical note 1402, 1993.
       /* 0: H2O 1:CO2 2:CH4 3:CO */
       /* compute partial pressures of radiating species */ 
       int is=0;
       for (i = 0; i < chem->nRad; i++) {
         is = chem->iRad[i];
         chem->pRad[i]  = flame->Y[j*flame->nVars+is]/chem->W[is]*flame->Mw[j]*flow->P;
       }
       
       /* compute emission coefficients */
       double T = flame->Y[k];    
       double Tm = 1000.0 / T;
       double T2 = T*T;
       chem->alpha[0] = -0.23093 + Tm * (-1.12390 + Tm * (9.41530    + Tm * 
           	     (-2.99880 + Tm * (0.51382  + Tm * -1.86840E-05))));  

       chem->alpha[1] =   18.741 + Tm * (-121.310 + Tm * (273.500    + Tm *
           	     (-194.050 + Tm * (56.310   + Tm * -5.8169))));

       chem->alpha[2] = 6.6334 - 0.0035686 * T + 1.6682e-08 * T2 +
                2.5611e-10 * T2 * T - 2.6558e-14 * T2 * T2;

       if (T < 750.0) 
           chem->alpha[3] = 4.7869+ T * (-0.06953 + T * (2.95775e-4 + T *
           	     (-4.25732e-7 + T * 2.02894e-10)));      
       else
           chem->alpha[3] = 10.09 + T * (-0.01183 + T * (4.7753e-6  + T *
           	     (-5.87209e-10+ T * -2.5334e-14)));
       double rad = 0.0;
       for (i = 0; i < chem->nRad; i++) 
          rad += chem->alpha[i] * chem->pRad[i] / 1.013E5;

#ifdef FM_VS_MDUC
       rad *= 2.0;
#endif /* FM_VS_MDUC */ 

       double fv=0.0;
       /* Soot radiation */
       if (chem->soot>0) {
          if (chem->sootMode<4) {
             /* compute soot volume fraction in [1/m^3] */
             /* NOTE: for soot mode == 1 we are solving M/rho */
             if (chem->sootMode == 1) 
                   fv = flame->Y[j*flame->nVars+flame->sSoot+1]*flame->rho[j]*SootMolarMass/SootDensity;
             else
                   fv = flame->sootMom[j*chem->nSootMom+1]*SootMolarMass/SootDensity;
          } /* soot mode < 4 */
       } /* soot */
       double alphas = -3.75e5 + 1735.0 * T;
       if (alphas < 0.0) alphas = 0.0;
#ifdef FM_VS_MDUC
       alphas *= 0.33*2.0;
#endif  /* FM_VS_MDUC */

       rad += alphas * fv;
#ifdef FM_VS_MDUC
       flame->radiation[j] = rad * 2.0*BOLTZMANN * (T2*T2);
#else
       flame->radiation[j] = rad * 2.0*BOLTZMANN * (T2*T2-pow(chem->TuRad,4));
#endif /* FM_VS_MDUC */
       } /* grid loop */
    } /* radiation */
    return;
  } // end init from data file

  if (driver->output_w)
      for (i=0; i<nGrid; i++) flame->gPDF[i] = 1.0;

  // initialise the grid in the x1 direction
  initGrid(flame->n1,0,1,io->nPivot1,io->gPivot1,io->gPerc1,
           &flame->x1,&flame->dx1);

  iX1 = (flame->n1-1)*flame->nVars;

  if (io->outputHtot && flame->H == NULL)
    flame->H = (double*) calloc(nGrid, sizeof(double));

  /* if 0D, interpolate the species from the mass fractions based on the mean */
  if (flame->nDim == 0 && flow->phiMean[1] >= 0.0) {

    double *Y_Z0 = (double*) calloc(chem->nSpecies, sizeof(double));
    double *Y_Z1 = (double*) calloc(chem->nSpecies, sizeof(double));

    for (i=0; i<chem->nSpecies; i++)
      Y_Z0[i] = Y_Z1[i]  = 0.0; 

    /* fill in the mass fractions */
    for (i=0; i<io->nOx;    i++) Y_Z0[io->oxIndex[i]]    = io->oxY[i];
    for (i=0; i<io->nFuel1; i++) Y_Z1[io->fuelIndex1[i]] = io->fuelY1[i];
    
    if (flow->phiMean[2] >= 0) {
      /* first interpolate the streams */
      double *Y_N  = (double*) calloc(chem->nSpecies, sizeof(double));
      double *Y_Z2 = (double*) calloc(chem->nSpecies, sizeof(double));
      double x1_N, T_N;

      for (i=0; i<chem->nSpecies; i++)
        Y_Z2[i] = Y_N[i] = 0.0; 
      
      for (i=0; i<io->nFuel2; i++)
        Y_Z2[io->fuelIndex2[i]] = io->fuelY2[i];

      /* intersection of mixing line with x1 axis */
      if (flow->phiMean[2] < 1.0)
        x1_N = flow->phiMean[1]/(1.0-flow->phiMean[2]);
      else 
        x1_N = 0.0;

      /* first mix the oxidizer side with the fuel 1 side */
      for (i=0; i<chem->nSpecies; i++)
        Y_N[i] = Y_Z0[i]*(1-x1_N) + Y_Z1[i]*x1_N;

      T_N = io->oxT*(1-x1_N) + io->fuelT1*x1_N;

      /* interpolate along mixing line to fuel 2 */
      for (i=0; i<chem->nSpecies; i++)
        flame->Y[i] = (1.0-flow->phiMean[2])*Y_N[i] + flow->phiMean[2]*Y_Z2[i];

      flame->Y[flame->sT] = (1.0-flow->phiMean[2])*T_N 
        + flow->phiMean[2]*io->fuelT2;

      free(Y_N); free(Y_Z2);
    }
    else { /* 1D interpolation */
      for (i=0; i<chem->nSpecies; i++)
        flame->Y[i] = Y_Z0[i]*(1-flow->phiMean[1]) + Y_Z1[i]*flow->phiMean[1];

      /* linearly mix the temperature */
      flame->Y[flame->sT] = io->oxT*(1-flow->phiMean[1])
        + io->fuelT1*flow->phiMean[1];
    }

    fprintf(stderr, "initial species info:\n");
    for (i=0; i<chem->nSpecies; i++)
      if (flame->Y[i] > 1.0e-7)
        fprintf(stderr, " %25s: %15.6e\n", chem->species[i], flame->Y[i]);
    
    if (driver->CHI_MODE == 3) {
      for (j=0; j < flame->n1; j++) 
          flame->chiA[j] = io->chiData[0].chiVal[j]; 
    }
    free(Y_Z0); free(Y_Z1);
    return;
  } /* end 0D interpolation */

  // allocate some local memory
  initH  = (double*) calloc(flame->n1, sizeof(double));
  
  // initialise the field from the boundary conditions given or from a given
  // flamelet solution
  
  // init from input file
  // calculate the enthalpy boundary conditions from temperature and species
  // mass fractions

  /* ----- set boundary conditions for oxidizer (left, x1=0) ----- */
  /* check to see if 1D solution should be interpolated */
  double *Y0  = (double*) calloc(chem->nSpecies, sizeof(double));
  double *Y1  = (double*) calloc(chem->nSpecies, sizeof(double));
  double *Yox = (double*) calloc(chem->nSpecies, sizeof(double));
  if (flame->nDim == 1 && flow->phiMean[2] >= 0.0 && io->nFuel2 > 0) {
    for (i=0; i<io->nOx; i++)    Y0[io->oxIndex[i]]    = io->oxY[i];
    for (i=0; i<io->nFuel2; i++) Y1[io->fuelIndex2[i]] = io->fuelY2[i];
    fprintf(stderr, " Initial oxidizer species (interpolated):\n");
    for (i=0; i<chem->nSpecies; i++) {
      Yox[i] = Y0[i]*(1-flow->phiMean[2]) + Y1[i]*flow->phiMean[2];
      if (Yox[i] > 1.0e-7)
        fprintf(stderr, " %25s: %15.6e\n", chem->species[i], Yox[i]);
    }
  }
  else { /* use boundary conditions directly */
    for (i=0; i<io->nOx; i++) Yox[io->oxIndex[i]] = io->oxY[i];
  }

  /* set the species */
  for (i=0; i<chem->nSpecies; i++)
    flame->Y[i] = Yox[i];

  initH[0] = 0.0;
  if (io->bcType[0] == 'T') { /* temperature specified */
    flame->Y[0 + flame->sT] = io->oxT;

    /*  compute enthalpy at boundary */
    chemComputeThermoData(chem->h, chem->cp, io->oxT, chem);
    for (i=0; i<chem->nSpecies; i++) {
      initH[0] += flame->Y[i]*chem->h[i];
    }
  }
  else if (io->bcType[0] == 'H') { /* enthalpy specified */
    initH[0] = io->oxH;
    
    /* compute temperature */
    io->oxT = calcTfromH(chem->nSpecies, chem->nSpecSolve, initH[0],
                         &flame->Y[0], 300, chem);
    flame->Y[0 + flame->sT] = io->oxT;
  }
  else if (io->bcType[0] == 'U' || io->bcType[0] == 'L')
    printf("Error: oxidiser boundary conidition specified incorrectly\n");
  
  if (driver->ENERGY) /* set the enthalpy if using diff-alg. eqs. */
    flame->H[0] = flame->Y[0 + flame->sT] = initH[0];

  /* ----- set boundary conditions for fuel 1 (right, Z1=1) ----- */
  for (i=0; i<io->nFuel1; i++)
    flame->Y[iX1+io->fuelIndex1[i]] = io->fuelY1[i];

  initH[flame->n1-1] = 0.0;
  if (io->bcType[1] == 'T') { // temperature specified
    flame->Y[iX1+flame->sT] = io->fuelT1;
    
    // compute enthalpy
    chemComputeThermoData(chem->h, chem->cp, io->fuelT1, chem);
    for (i=0; i<io->nFuel1; i++)
      initH[flame->n1-1] += io->fuelY1[i]*chem->h[io->fuelIndex1[i]];
  }
  else if (io->bcType[1] == 'H') { // enthalpy specified
    initH[flame->n1-1] = io->fuelH1;
    
    // compute temperature
    io->fuelT1 = calcTfromH(chem->nSpecies, chem->nSpecSolve,
                            initH[flame->n1-1], &flame->Y[iX1], 300, chem);
    flame->Y[iX1+flame->sT] = io->fuelT1;
  }
  else if (io->bcType[1] == 'L') { // liquid fuel specified
    // compute the enthalpy of the vapour at the reference temperature
    chemComputeThermoData(chem->h, chem->cp, io->fuelTref1, chem);
    H_ref = 0.0;
    for (i=0; i<io->nFuel1; i++)
      H_ref += io->fuelY1[i]*chem->h[io->fuelIndex1[i]];
    
    // compute the enthalpy of the vapour at fuel boundary
    initH[flame->n1-1] =   io->cpLiq1*(io->fuelT1 - io->fuelTref1)
                         - io->hVap1 + H_ref;
    
    // compute temperature of vapour at boundary
    io->fuelT1 = calcTfromH(chem->nSpecies, chem->nSpecSolve,
                            initH[flame->n1-1], &flame->Y[iX1], 300, chem);
    
    flame->Y[iX1+flame->sT] = io->fuelT1;
    
    if (MDUC_PRINT_PROC) printf("Fuel 1 specified as liquid, vapour temperature set to: %8.4f K\n",
           io->fuelT1);
  }
  else if (io->bcType[1] == 'D') { // liquid fuel using droplet
    /* double *dropData = calloc(6, sizeof(double)); */
    double hfg, cpG, cpL, Psat, Q, B, Xfs, Yfs, Ts, Hs, Textrap;
    
/*       dropData[0] = io->fuelT1;      // liquid fuel temperature */
/*       dropData[1] = flow->pressure;  // ambient pressure */
/*       dropData[2] = io->oxT;         // ambient temperature */
/*       dropData[3] = io->Yinit[io->sFuel]; */
/*       dropData[4] = molecular weight fuel */
/*       dropData[5] = molecular weight oxidizer */

    // compute the molecular weight of the fuel and oxidizer
    for (i=0, io->Mfuel1=0.0; i<io->nFuel1; i++)
      io->Mfuel1 += io->fuelY1[i]*chem->W[io->fuelIndex1[i]];
    for (i=0, io->Mox=0.0; i<io->nOx; i++)
      io->Mox    += io->oxY[i]*chem->W[io->oxIndex[i]];
    io->oxP = flow->P;
    
    fprintf(stderr, "Mol. weight fuel:     %8.6f [g/mol]\n", io->Mfuel1);
    fprintf(stderr, "Mol. weight oxidizer: %8.6f [g/mol]\n", io->Mox);
    fprintf(stderr, "Oxidizer pressure:    %15.6e [Pa]\n", io->oxP);
    
    // compute the theoretical surface temperature
    Ts = brent(MDUCcalcSpaldingFunc, 183, 538, 1.0e-4, io);
    NIST_getThermo(io->satData, &hfg, &cpG, &cpL, &Psat, Ts);
    Q = hfg + cpL*(Ts - io->fuelT1);
    B = cpG*(io->oxT - Ts)/Q;
    Xfs = Psat/io->oxP;
    Yfs = (B + 0.0)/(1.0 + B);
    
    fprintf(stderr, "Droplet surface properties\n");
    fprintf(stderr, "    temperature:        %8.3f\n", Ts);
    fprintf(stderr, "    fuel mass fraction: %8.3f\n", Yfs);
    fprintf(stderr, "    fuel mole fraction: %8.3f\n", Xfs);
    fprintf(stderr, "Droplet Spalding number :    %8.3f\n", B);
    
    double *YoxS   = (double*) calloc(io->nOx, sizeof(double));
    double *YfuelS = (double*) calloc(io->nFuel1, sizeof(double));
    for (i=0; i<io->nOx; i++)
      YoxS[i] = linInterp(Yfs,flame->x1[0],io->oxY[i],flame->x1[flame->n1-1],0);
    for (i=0; i<io->nFuel1; i++)
      YfuelS[i] = linInterp(Yfs, flame->x1[0],0.0, flame->x1[flame->n1-1],
                            io->fuelY1[i]);
    
    // compute the enthalpy
    chemComputeThermoData(chem->h, chem->cp, Ts, chem);
    for (i=0, Hs=0.0; i<io->nOx; i++)
      Hs += YoxS[i]*chem->h[io->oxIndex[i]];
    for (i=0; i<io->nFuel1; i++)
      Hs += YfuelS[i]*chem->h[io->fuelIndex1[i]];
    
    // extrapolate the enthalpy to the boundary
    double Hslope = (Hs-initH[0])/Yfs;
    initH[flame->n1-1] = initH[0] + Hslope;
    
    Textrap = (Ts-io->oxT)/Yfs + io->oxT;
    fprintf(stderr, "Extrapolated temperature: %8.3f K\n", Textrap);
    
    fprintf(stderr, "Fuel enthalpy extrapolated to: %15.6e\n",
            initH[flame->n1-1]);
    
    io->fuelT1 = calcTfromH(chem->nSpecies, chem->nSpecSolve,
                            initH[flame->n1-1], &flame->Y[iX1], 300, chem);
    flame->Y[iX1+flame->sT] = io->fuelT1;
    
    fprintf(stderr, "Fuel 1 specified as liquid droplet:\n");
    fprintf(stderr, "   setting vapour temperature to %8.4f K\n", io->fuelT1);
  }
  else
    printf("Error: fuel 1 boundary condition specified incorrectly\n");

  if (driver->ENERGY) /* set the enthalpy if using diff-alg. eqs. */
    flame->H[flame->n1-1] = flame->Y[iX1+flame->sT] = initH[flame->n1-1];

  if (driver->H_DIM == 1) { /* set flow bounds to enthalpy, just in case */
    flow->phi1min = initH[0];
    flow->phi1max = initH[flame->n1-1];
  }

  /* Ignite: TODO: make this more general (multiple fuel species)*/
  double W,nuO2,nuH2O,nuCO2;
  int iFuel,iO2,iH2O,iCO2;
  double *X0  = (double*) calloc(chem->nSpecies, sizeof(double));

  if (driver->ignite == 1) {
     nuH2O = 2.0*chem->fuelH;
     nuCO2 = chem->fuelC;
     nuO2  = (nuH2O/2.0+nuCO2-chem->fuelO);
     iFuel = chem->iFuel; 
     /* TODO: Use chem->si */
     iO2   = getSpeciesIndex(chem->nSpecies, chem->species, "O2");
     iH2O  = getSpeciesIndex(chem->nSpecies, chem->species, "H2O");
     iCO2  = getSpeciesIndex(chem->nSpecies, chem->species, "CO2");

  }

  /* ----- linearly interpolate enthalpy and mass fractions ----- */
  for (j=1; j<flame->n1-1; j++) {
    k     = j*flame->nVars;
    k_pre = (j-1)*flame->nVars;


    /* linearly interpolate the enthalpy */
    initH[j] = linInterp(flame->x1[j], flame->x1[0], initH[0],
                         flame->x1[flame->n1-1], initH[flame->n1-1]);
      
    /* interpolate the species mass fractions */
    for (i=0; i<chem->nSpecies; i++)
      flame->Y[k+i] = linInterp(flame->x1[j], flame->x1[0], flame->Y[i],
                                flame->x1[flame->n1-1], flame->Y[iX1+i]);

    W = calcMixtureW(chem->nSpecSolve,&flame->Y[k], chem->W);

    /* Ignite: TODO: make this more general (multiple fuel species)*/
    if (driver->ignite == 1) {

       /* convert to mole fraction */
       for (i=0;i<chem->nSpecies;i++)
	  X0[i] = flame->Y[k+i] * W / chem->W[i];

       if (X0[iO2] >= nuO2 * X0[iFuel]) {// lean 
	  X0[iO2]  -= X0[iFuel]*nuO2;
	  X0[iH2O] += X0[iFuel]*nuH2O;
	  X0[iCO2] += X0[iFuel]*nuCO2;
          X0[iFuel] = 0.0;
       }
       else {// rich
          X0[iFuel] -= X0[iO2] / nuO2;
	  X0[iH2O]  += X0[iO2] / nuO2*nuH2O;
	  X0[iCO2]  += X0[iO2] / nuO2*nuCO2;
	  X0[iO2]   = 0.0;
       }
       /* renormalize mole fractions and recalculate W*/	
       double sum=0.0;
       for (i=0;i<chem->nSpecies;i++)
	  sum+=X0[i];

       W = 0.0;
       for (i=0;i<chem->nSpecies;i++) {
   	  X0[i]/=sum;
	  W += X0[i] * chem->W[i];
       }
       
       /* convert to mass fraction */       
       for (i=0;i<chem->nSpecies;i++)
	  flame->Y[k+i] = X0[i] * chem->W[i]/W;

    }

    /* get the temperature based on the enthalpy */
    if (driver->ENERGY)
      flame->H[j] = flame->Y[k+flame->sT] = initH[j];
    else
      flame->Y[k+flame->sT] = calcTfromH(chem->nSpecies, chem->nSpecSolve,
                                         initH[j], &flame->Y[k],
                                         flame->Y[k_pre+flame->sT], chem);
  } /* end loop over x1 */

  /* ----- set up the 2D field ----- */
  if (flame->nDim == 2) {
    // initialise the grid in the second dimension
    initGrid(flame->n2, 0, 1, io->nPivot2, io->gPivot2, io->gPerc2,
             &flame->x2, &flame->dx2);
    initX2dimension(driver, io, chem, flame, flow, NULL);
  }

  /* initialise the metric quantities */
  initMetric(flame);

  /* initialize soot quantities if necessary */
  if (chem->soot) {
    if (chem->sootMode <= 3) {
      for (k=0; k<flame->n2; k++) {
        for (j=0; j<flame->n1; j++) {
          initSootMoments(chem->nSootMom,
                          &flame->sootMom[(k*flame->n1+j)*flame->nSootMom],
			  1.0E-20);
        //  initSootMoments(chem->nSootMom,
        //                  &flame->sootMom[(k*flame->n1+j)*flame->nSootMom],
	//		  driver->SABSTOL);

          if (chem->sootMode == 1) {
            /* copy initial soot moments into solver variables */
            for (i=0; i<chem->nSootMom; i++)
              flame->Y[(k*flame->n1+j)*flame->nVars+flame->sSoot+i] = 
                flame->sootMom[(k*flame->n1+j)*flame->nSootMom+i];
          }
        } /* end loop over columns */
      } /* end loop over rows */
    } /* soot mode */
  } /* soot */

  /* compute the current properties */
  double wmix;
  for (k=0; k<flame->n2; k++) {
    row = k*flame->n1*flame->nVars;
    r   = k*flame->n1;
    for (j=0; j<flame->n1; j++) {
      col  = j*flame->nVars;
      wmix = calcMixtureW(chem->nSpecSolve, &flame->Y[row+col], chem->W);
      flame->Mw[r+j] = wmix;
      flame->rho[r+j] = wmix*flow->P/(flame->RGAS*flame->Y[row+col+flame->sT]);
      if (chem->readTrans)
        flame->mu[r+j]  = calcMixViscosity(chem->nSpecies,
                                           flame->Y[row+col+flame->sT],
                                           &flame->Y[row+col], chem);

    } /* end loop over columns */
  } /* end loop over rows */

  if (chem->soot) {
     if (chem->sootMode >= 3) {
        for (k=0; k<flame->n2; k++) {
          row = k*flame->n1*flame->nVars;
          r   = k*flame->n1;
          for (j=0; j<flame->n1; j++) {
            col  = j*flame->nVars;
	    double sgCoeff, oxCoeff,o2Coeff, dimProdRate;
            sootHMOM_GasPhase_Coeff(flame->Y[row+col+flame->sT], flame->mu[r+j],
            	                flame->rho[r+j], flame->Mw[r+j], chem->sootRxnConc, 
            			&sgCoeff, &oxCoeff,&o2Coeff, &dimProdRate);

	    flame->sgCoeff[r+j] 	    = sgCoeff;
	    flame->oxCoeff[r+j]  	    = oxCoeff;
	    flame->o2Coeff[r+j] 	    = o2Coeff;
	    flame->dimProdRate[r+j] 	    = dimProdRate;
            flame->sqrtT[r+j] 	    = sqrt(flame->Y[row+col+flame->sT]);
            flame->T_mu[r+j]  	    = flame->Y[row+col+flame->sT]/flame->mu[r+j];
            flame->muSqrtMw_rhoSqrtT[r+j] = flame->mu[r+j]*sqrt(flame->Mw[r+j]/
            			(1000.0*flame->Y[row+col+flame->sT]*flame->rho[r+j]*flame->rho[r+j]));

	  } 
        }
     }
  }
  if (driver->CHI_MODE == 3) {
    for (f=0; f < io->nFiles; f++) {/* loop over files */
      /* allocate memory */ 
      if (io->chiData[f].chiRaw != NULL)
      io->chiData[f].chiVal = (double*) calloc(flame->n1, sizeof(double));
      
      if (io->chiData[f].pdfRaw != NULL)
         io->chiData[f].pdfVal = (double*) calloc(flame->n1, sizeof(double));

      if (io->chiData[f].hSrcRaw != NULL)
         io->chiData[f].hSrcVal = (double*) calloc(flame->n1, sizeof(double));

      if (io->chiData[f].rhoDotRaw != NULL)
         io->chiData[f].rhoDotVal = (double*) calloc(flame->n1, sizeof(double));

      if (io->chiData[f].TDotRaw != NULL)
         io->chiData[f].TDotVal = (double*) calloc(flame->n1, sizeof(double));
      
      for (j=0; j < flame->n1; j++) { /* loop over target grid */
        found = 0;i=0; /* FIXME */
        while (found<1 && i<io->chiData[f].n1-1) { 
          i++;
          if (io->chiData[f].x1[i]>flame->x1[j]) found=1;
        }
        if (found>0 && i>0) {
           if (io->chiData[f].chiRaw != NULL)
	   io->chiData[f].chiVal[j] = (io->chiData[f].chiRaw[i] - io->chiData[f].chiRaw[i-1])
      		     /(io->chiData[f].x1[i] - io->chiData[f].x1[i-1])
      		     *(flame->x1[j] - io->chiData[f].x1[i-1]) + 
      		      io->chiData[f].chiRaw[i-1];
           if (io->chiData[f].pdfRaw != NULL)
	      io->chiData[f].pdfVal[j] = (io->chiData[f].pdfRaw[i] - io->chiData[f].pdfRaw[i-1])
      		                        /(io->chiData[f].x1[i] - io->chiData[f].x1[i-1])
      		                        *(flame->x1[j] - io->chiData[f].x1[i-1]) + 
      		                          io->chiData[f].pdfRaw[i-1];
           if (io->chiData[f].hSrcRaw != NULL)
	      io->chiData[f].hSrcVal[j] = (io->chiData[f].hSrcRaw[i] - io->chiData[f].hSrcRaw[i-1])
      		                        /(io->chiData[f].x1[i] - io->chiData[f].x1[i-1])
      		                        *(flame->x1[j] - io->chiData[f].x1[i-1]) + 
      		                          io->chiData[f].hSrcRaw[i-1];
           if (io->chiData[f].rhoDotRaw != NULL)
	      io->chiData[f].rhoDotVal[j] = (io->chiData[f].rhoDotRaw[i] - io->chiData[f].rhoDotRaw[i-1])
      		                        /(io->chiData[f].x1[i] - io->chiData[f].x1[i-1])
      		                        *(flame->x1[j] - io->chiData[f].x1[i-1]) + 
      		                          io->chiData[f].rhoDotRaw[i-1];
           if (io->chiData[f].TDotRaw != NULL)
	      io->chiData[f].TDotVal[j] = (io->chiData[f].TDotRaw[i] - io->chiData[f].TDotRaw[i-1])
      		                        /(io->chiData[f].x1[i] - io->chiData[f].x1[i-1])
      		                        *(flame->x1[j] - io->chiData[f].x1[i-1]) + 
      		                          io->chiData[f].TDotRaw[i-1];
        }
        else {
           if (io->chiData[f].chiRaw != NULL)
   	      io->chiData[f].chiVal[j] = 0.0;
           if (io->chiData[f].pdfRaw != NULL)
	      io->chiData[f].pdfVal[j] = 0.0;
           if (io->chiData[f].hSrcRaw != NULL)
	      io->chiData[f].hSrcVal[j] = 0.0;
           if (io->chiData[f].rhoDotRaw != NULL)
	      io->chiData[f].rhoDotVal[j] = 0.0;
           if (io->chiData[f].TDotRaw != NULL)
	      io->chiData[f].TDotVal[j] = 0.0;
        }
      }
      /* set first and last value to zero */
      if (io->chiData[f].chiRaw != NULL) {
         io->chiData[f].chiVal[0] = 0.0;
         io->chiData[f].chiVal[flame->n1-1] = 0.0;
      }
    }
  }
  free(Y0); free(Y1); free(Yox); free(initH);free(X0);
}
/* -------------------------------------------------------------------------- */

void initTime(driveStruct *driver, timeStruct *time)
{
  if (driver->TIME_MODE != 0)
    {
      if (driver->TIME_MODE == 1)
        {
          
        }
      else if (driver->TIME_MODE == 2)
        {
          
        }
      else if (driver->TIME_MODE == 3)
        {
          
        }
      //else
      //	printf(stderr, "!!!!! ERROR: TIME_MODE not valid\n");
    } /* end TIME_MODE != 0 */
}
/* -------------------------------------------------------------------------- */

void initSUNDIALS(ioStruct *io, sdialStruct *sdial)
{
  char tmp[128];
  /* ensure memory is null */
  sdial->mem = NULL;
  sdial->Ydata = sdial->YpData = sdial->resData = NULL;
  sdial->YscaleData = sdial->FscaleData = NULL;
  
  paramGetString(sdial->solver,"sundials-solver","cvode",io->parser);
  strcpy(io->sdialSolver,sdial->solver);

  /* store linear system solver */
  paramGetString(sdial->method,"sundials-method","gmres",io->parser);
  if (!strcmp(sdial->method, "dense")) {
    fprintf(stderr, "Warning: using dense solver "
            " (this is inefficient in most cases)\n");
    sdial->iMethod = 0;
  }
  else if (!strncmp(sdial->method, "band", 4))
    sdial->iMethod = 1;
  else if (!strcmp(sdial->method, "gmres"))
    sdial->iMethod = 2;
  else if (!strcmp(sdial->method, "bicgstab"))
    sdial->iMethod = 3;
  else if (!strcmp(sdial->method, "tfqmr"))
    sdial->iMethod = 4;

#ifdef ACCESS_MPI
  /* make sure that we're not trying to use a serial method when parallel */
  if (sdial->iMethod < 2) {
    if (MDUC_PRINT_PROC) {
      fprintf(stderr, "Warning: Sundials method %s not available in parallel: "
              "using GMRES instead\n", sdial->method);
    }
    strcpy(sdial->method, "gmres");
    sdial->iMethod = 2;
  }
#endif /* ACCESS_MPI */

  /* assign preconditioner  */
  paramGetString(sdial->precnd,"sundials-preconditioner","left",io->parser);
  if (!strcmp(sdial->precnd, "left"))
    sdial->iPreCnd = 1;
  else if (!strcmp(sdial->precnd, "right"))
    sdial->iPreCnd = 2;
  else if (!strcmp(sdial->precnd, "both"))
    sdial->iPreCnd = 3;

  /* for 2D problems, deterimine integration scheme to use */
  paramGetString(tmp,"sundials-2d","ADI",io->parser);
  if (!strcmp(tmp,"ADI"))
    sdial->ADI = 1;
  else if (!strncmp(tmp,"implicit",4))
    sdial->ADI = 0;
  else if (MDUC_PRINT_PROC)
    fprintf(stderr, "Warning: 2D integration method %s not defined\n",tmp);

  if (MDUC_PRINT_PROC == 1) {
    fprintf(stdout,"----- Sundials Setup Summary -----\n");
    fprintf(stdout,"  Solver: %s\n", sdial->solver);
    fprintf(stdout,"  Method: %s\n", sdial->method);
    if (sdial->iMethod > 1)
      fprintf(stdout,"  Preconditioner: %s\n", sdial->precnd);
    fprintf(stdout,"----------------------------------\n");
  }
}
/* -------------------------------------------------------------------------- */

void initCVODE(int dim, int nMesh, int nVars, double tInit, double *YInit,
               driveStruct *driver, sdialStruct *cvode)
{
  int flag;  // for checking return values from cvode functions
  long int mu = nVars;
  long int ml = nVars;

  cvode->mem  = NULL;         // set the memory to NULL
  cvode->Y = cvode->Yp = NULL;
  cvode->Yscale = cvode->Fscale = cvode->tmp = NULL;
  
  cvode->nEq = nMesh*nVars;  // set the number of equations in the system
  cvode->row = cvode->col = 0;

  /* DEBUG: add engine components */
  if (driver->ENGINE) {
    cvode->nEq += 3;
    if (driver->ENGINE > 1) /* SI */
      cvode->nEq += 3;
    if (driver->PRES_MODE == 3) /* implicit pressure */
      cvode->nEq += 1;
  }

  /* set the solution vector and initialise from solution YInit
     also assign pointer to first memory location in Y
     NOTE: use compiler directive here later for parallel */
  if (dim == 1)
    cvode->Y = N_VMake_Serial(cvode->nEq, YInit);
  else if (dim == 0 && YInit)
    cvode->Y = N_VMake_Serial(cvode->nEq, YInit);
  else
    cvode->Y = N_VNew_Serial(cvode->nEq);

  cvode->Ydata = NV_DATA_S(cvode->Y);

  /* create the memory object for CVODE:
     The system is stiff due to chemistry, thus choose solver appropriately
       CV_BDF:    linear multistep method (Backward Differentiation Formulas) 
       CV_NEWTON: Newton iteration (modified newton iteration for banded) */
#ifdef SUNDIALS_FROM_4_0
  cvode->mem = CVodeCreate(CV_BDF);
#else
  cvode->mem = CVodeCreate(CV_BDF, CV_NEWTON);
#endif /* SUNDIALS_FROM_4_0 */

  /* set right-hand-side function based on flamelet dimension */
  if (dim == 0)
    cvodeRHS = cvodeRHS_0D;
  else if (dim == 1)
    cvodeRHS = cvodeFlamelet1D;
  else if (dim == 2)
    cvodeRHS = cvodeRHS_2D;

  /* set the general functions */
  SDialGetNumSteps = CVodeGetNumSteps;
  SDialGetLastStep = CVodeGetLastStep;
  SDialSetStopTime = CVodeSetStopTime;
  SDialSetInitStep = CVodeSetInitStep;
  SDialFree        = CVodeFree;

  /* allocate internal memory object:
     initialise time:tInit
     RHS function:   cvodeRHS
     CV_SS:          scalar relative and absolute tolerances
     reltol:         relative tolerance
     abstol:         absolute tolerance
     
     NOTE: Sundials v2.4.0 is now default, support for v2.3.0 available with
           compiler flag SUNDIALS_2_3 */
#ifdef SUNDIALS_2_3
  flag = CVodeMalloc(cvode->mem, cvodeRHS, tInit, cvode->Y,
                     CV_SS, driver->RELTOL, &driver->ABSTOL);

  // set pointer to parameters (structs) needed by the RHS function
  flag = CVodeSetFdata(cvode->mem, cvode);

  // set the solver to return after every internal timestep
  cvode->itask = CV_ONE_STEP_TSTOP;
  cvode->iret  = CV_TSTOP_RETURN;
#else /* Sundials v2.4.0 (default) */
  flag = CVodeInit(cvode->mem, cvodeRHS, tInit, cvode->Y);
  flag = CVodeSStolerances(cvode->mem, driver->RELTOL, driver->ABSTOL);

  // set pointer to parameters (structs) needed by the RHS function
  flag = CVodeSetUserData(cvode->mem, cvode);

  // set the solver to return after every internal timestep
  cvode->itask = CV_ONE_STEP;
  cvode->iret  = CV_TSTOP_RETURN;
#endif /* SUNDIALS_2_3 */

  if (flag != CV_SUCCESS) {
    if (flag == CV_MEM_NULL)
      fprintf(stderr, "CVODE ERROR: CVodeMalloc requires CVodeCreate");
    printf("exiting: unable to allocate with CVodeMalloc\n");
    exit(1);
  }
  
  // defined banded matrix, set number of equations, upper/lower bandwidths
  if (dim == 0) {
#ifdef SUNDIALS_FROM_3_0
    cvode->A = SUNDenseMatrix(cvode->nEq,cvode->nEq);
#else /* SUNDIALS LOWER THAN 3.0 */
    flag = CVDense(cvode->mem, cvode->nEq);
#endif /* SUNDIALS_FROM_3_0 */
    cvode->row = cvode->col = 0;
  }
  else {
    if (cvode->iMethod == 0) { /* dense */
      mducWarn("Warning: using dense solver"
               "(this is inefficient in most cases)", "cvodeInitFlamelet");
#ifdef SUNDIALS_FROM_3_0
      cvode->A = SUNDenseMatrix(cvode->nEq,cvode->nEq);
#else /* SUNDIALS LOWER THAN 3.0 */
      flag = CVDense(cvode->mem,cvode->nEq);
#endif /* SUNDIALS_FROM_3_0 */
      SDialGetNumJacEvals = CVDlsGetNumJacEvals;
    }
    else if (cvode->iMethod == 1) { /* banded */
#ifdef SUNDIALS_FROM_3_0
      cvode->sunIndTypeTmp = nVars;
#ifdef SUNDIALS_FROM_4_0
      cvode->A  = SUNBandMatrix(cvode->nEq,cvode->sunIndTypeTmp,cvode->sunIndTypeTmp);
#else
      cvode->A  = SUNBandMatrix(cvode->nEq,cvode->sunIndTypeTmp,cvode->sunIndTypeTmp,
			        cvode->sunIndTypeTmp+cvode->sunIndTypeTmp);
#endif
      cvode->LS = SUNBandLinearSolver(cvode->Y,cvode->A);
      flag      = CVDlsSetLinearSolver(cvode->mem, cvode->LS, cvode->A);
#else /* SUNDIALS LOWER THAN 3.0 */
      flag = CVBand(cvode->mem,cvode->nEq,nVars,nVars);
#endif /* SUNDIALS_FROM_3_0 */
      SDialGetNumJacEvals = CVDiagGetNumRhsEvals;
    }
    else if (cvode->iMethod == 2) { /* GMRES */
#ifdef SUNDIALS_FROM_3_0
      cvode->LS = SUNSPGMR(cvode->Y,cvode->iPreCnd,0);
      flag = CVSpilsSetLinearSolver(cvode->mem, cvode->LS);
#else /* SUNDIALS LOWER THAN 3.0 */
      flag = CVSpgmr(cvode->mem,cvode->iPreCnd,0);
#endif /* SUNDIALS_FROM_3_0 */
      SDialGetNumJacEvals = CVSpilsGetNumJtimesEvals;
    }
    else if (cvode->iMethod == 3) { /* BICGSTAB */
#ifdef SUNDIALS_FROM_3_0
      cvode->LS = SUNSPBCGS(cvode->Y,cvode->iPreCnd,0);
      flag = CVSpilsSetLinearSolver(cvode->mem,cvode->LS);
#else /* SUNDIALS LOWER THAN 3.0 */
      flag = CVSpbcg(cvode->mem,cvode->iPreCnd,0);
#endif /* SUNDIALS_FROM_3_0 */
      SDialGetNumJacEvals = CVSpilsGetNumJtimesEvals;
    }
    else if (cvode->iMethod == 4) { /* TFQMR */
#ifdef SUNDIALS_FROM_3_0
      cvode->LS = SUNSPTFQMR(cvode->Y,cvode->iPreCnd,0);
      flag =  CVSpilsSetLinearSolver(cvode->mem, cvode->LS);
#else /* SUNDIALS LOWER THAN 3.0 */
      flag = CVSptfqmr(cvode->mem,cvode->iPreCnd,0);
#endif /* SUNDIALS_FROM_3_0 */
      SDialGetNumJacEvals = CVSpilsGetNumJtimesEvals;
    }
    
    if (cvode->iMethod > 1) { /* initialize precondition for Krylov methods */
#ifdef ACCESS_MPI
#ifdef SUNDIALS_FROM_3_0
      cvode->sunIndTypeTmp = nVars;
      flag = CVBBDPrecInit(cvode->mem,cvode->nEqL,cvode->sunIndTypeTmp,cvode->sunIndTypeTmp,
			   cvode->sunIndTypeTmp,cvode->sunIndTypeTmp,0.0,cvodeRHS_L,NULL);
#else /* SUNDIALS LOWER THAN 3.0 */
      flag = CVBBDPrecInit(cvode->mem,cvode->nEqL,nVars,nVars,nVars,nVars,
                           0.0,cvodeRHS_L,NULL);
#endif /*SUNDIALS_FROM_3_0 */
#else /* SERIAL */
      flag = CVBandPrecInit(cvode->mem,cvode->nEq,nVars,nVars);
#endif /* ACCESS_MPI */
    }
  }
}
/* -------------------------------------------------------------------------- */

void initIDA(int dim, int nMesh, int nVars, double tInit, double reltol,
             double abstol, double *YInit, sdialStruct *ida)
{
  int j, flag;  // for checking return values from cvode functions
  // dereference pointers to structs
  driveStruct  *driver = (driveStruct*)  ida->args[0];
  chemStruct   *chem   = (chemStruct*)   ida->args[1];
  flowStruct   *flow   = (flowStruct*)   ida->args[2];
  flameStruct  *flame  = (flameStruct*)  ida->args[3];

  ida->mem = NULL;         // set the memory to NULL
  ida->nEq = nMesh*nVars;  // set the number of equations in the system

  /* set the solution vector and initialise from solution YInit
     also assign pointer to first memory location in Y
     NOTE: use compiler directive here later for parallel */
  if (dim == 1)
    ida->Y    = N_VMake_Serial(ida->nEq, YInit);
  else
    ida->Y    = N_VNew_Serial(ida->nEq);
  ida->Yp = N_VNew_Serial(ida->nEq);

  ida->Ydata  = NV_DATA_S(ida->Y);
  ida->YpData = NV_DATA_S(ida->Yp);

  /* create the memory object for IDA */
  ida->mem = IDACreate();

  /* set right-hand-side function based on flamelet dimension */
  if (dim == 0)
    idaRHS = idaRHS_0D;
  else if (dim == 1)
    idaRHS = idaFlamelet1D;
  else if (dim == 2)
    idaRHS = idaFlamelet2D;

  /* compute the initial derivative */
  double *temp;
  if (dim == 1)
    {
      // ensure temperature is up-to-date
      if (driver->ENERGY)
        {
          updateEnergy1D(driver,chem,flame,flow,0,ida->Ydata,ida->YpData);
          temp = flame->T;
        }
      else
        {
          temp = (double*) calloc(flame->n1, sizeof(double));
          extractVarArray(nVars, nMesh, 1, flame->sT, ida->Ydata, temp);
        }

      // get species derivatives
      flameletY_RHS_1D(ida,driver,chem,flame,flow,temp,ida->Ydata,ida->YpData);

      // energy equation
      if (driver->ENERGY)
        {
          for (j=0; j<flame->n1; j++)
            ida->YpData[j*flame->nVars + flame->sT] = flame->H[j];
        }
      else
        flameletT_RHS_1D(driver,chem,flame,flow,ida->Ydata,ida->YpData);
    }

  /* set the general functions */
  SDialGetNumSteps = IDAGetNumSteps;
  SDialGetLastStep = IDAGetLastStep;
  SDialSetStopTime = IDASetStopTime;
  SDialSetInitStep = IDASetInitStep;
  SDialFree        = IDAFree;

  /* allocate internal memory object:
     initial time: tInit

     NOTE: Sundials v2.4.0 is now default, support for v2.3.0 available with
           compiler flag SUNDIALS_2_3 */
#ifdef SUNDIALS_2_3
  flag = IDAMalloc(ida->mem, idaRHS, tInit, ida->Y, ida->Yp,
                   IDA_SS, reltol, &abstol);

  // set pointer to parameters (structs) needed by the RHS function
  flag = IDASetRdata(ida->mem, ida);

  ida->itask = IDA_ONE_STEP_TSTOP;
  ida->iret  = IDA_TSTOP_RETURN;
#else /* Sundials v2.4.0 (default) */
  flag = IDAInit(ida->mem, idaRHS, tInit, ida->Y, ida->Yp);
  flag = IDASStolerances(ida->mem, reltol, abstol);

  // set pointer to parameters (structs) needed by the RHS function
  flag = IDASetUserData(ida->mem, ida);

  ida->itask = IDA_ONE_STEP;
  ida->iret  = IDA_TSTOP_RETURN;
#endif /* SUNDIALS_2_3 */

  if (flag != IDA_SUCCESS)
    {
      if (flag == IDA_MEM_NULL)
        fprintf(stderr, "IDA ERROR: IDAMalloc requires IDACreate");
      printf("exiting: unable to allocate with IDAMalloc\n");
      exit(1);
    }
  
  // defined banded matrix, set number of equations, upper/lower bandwidths
  if (dim == 0) {
#ifdef SUNDIALS_FROM_3_0
    ida->A  = SUNDenseMatrix(ida->nEq,ida->nEq);
    ida->LS = SUNDenseLinearSolver(ida->Y,ida->A);
    flag    = IDADlsSetLinearSolver(ida->mem,ida->LS,ida->A);
#else
    flag = IDADense(ida->mem, ida->nEq);
#endif /* SUNDIALS_FROM_3_0 */
  }
  else {
#ifdef SUNDIALS_FROM_3_0
    ida->sunIndTypeTmp = flame->nVars;
#ifdef SUNDIALS_FROM_4_0
    ida->A  = SUNBandMatrix(ida->nEq, ida->sunIndTypeTmp,ida->sunIndTypeTmp);
#else
    ida->A  = SUNBandMatrix(ida->nEq, ida->sunIndTypeTmp,ida->sunIndTypeTmp,
			    ida->sunIndTypeTmp+ida->sunIndTypeTmp);
#endif
    ida->LS = SUNBandLinearSolver(ida->Y,ida->A);
    flag    = IDADlsSetLinearSolver(ida->mem,ida->LS,ida->A);
#else
    flag = IDABand(ida->mem, ida->nEq, nVars, nVars);
#endif /* SUNDIALS_FROM_3_0 */
  }
}
/* -------------------------------------------------------------------------- */

void initFlow(driveStruct *driver, flameStruct *flame, flowStruct *flow)
{
  int i,j;
  double PI = 3.1415926535897931;
  /* if the strain rate is set, assume that it will be used as the reference */
  if (flow->strain[0] > 0.0 || flow->strain[1] > 0.0 || flow->strain[2] > 0.0) {
    if (driver->CHI_MODE == 0) {
      fprintf(stderr,"Error: log profile not yet compatible with strain\n");
      exit(2);
    }
    
    flow->chiRef[0] = flow->strain[0]/PI;
    flow->chiRef[1] = flow->strain[1]/PI;
    flow->chiRef[2] = flow->strain[2]/PI;
    
    /* set reference to 1 so that the second term in erfc is not computed */
    flow->ZRef[0] = flow->ZRef[1] = flow->ZRef[2] = 1.0;
  }

  /* allocate memory for scalar dissipation rates if not already done */
  if (flow->chi1 == NULL)
    flow->chi1 = (double*) calloc(flame->n1*flame->n2, sizeof(double));
  if (flame->nDim == 2) {
    if (flow->chi2 == NULL)
      flow->chi2 = (double*) calloc(flame->n1*flame->n2, sizeof(double));
    if (flow->chi12 == NULL)
      flow->chi12 = (double*) calloc(flame->n1*flame->n2, sizeof(double));
  }

  /* allocate memory for scalar dissipation rates if not already done */
  if (flow->chi1_m1 == NULL)
    flow->chi1_m1 = (double*) calloc(flame->n1*flame->n2, sizeof(double));
  if (flame->nDim == 2) {
    if (flow->chi2_m1 == NULL)
      flow->chi2_m1 = (double*) calloc(flame->n1*flame->n2, sizeof(double));
    if (flow->chi12_m1 == NULL)
      flow->chi12_m1 = (double*) calloc(flame->n1*flame->n2, sizeof(double));
  }
  
  if (flow->pdf == NULL)
    flow->pdf = (double*) calloc(flame->n1*flame->n2, sizeof(double));

  /* allocate memory for the conditional mean of enthalpy source term */
  if (driver->H_DIM == 1) {
    flow->avgSdot = (double*) calloc(flame->n2, sizeof(double));
    flow->avgH    = (double*) calloc(flame->n2, sizeof(double));
  }
  else if (driver->H_DIM == 2) {
    flow->avgSdot = (double*) calloc(flame->n1, sizeof(double));
    flow->avgH    = (double*) calloc(flame->n1, sizeof(double));
  }
  if (driver->CHI_MODE == 2) {
    flow->strainCondZ = (double*) calloc(flame->n1*flame->n2, sizeof(double));
    flow->avgSdot_M   = (double*) calloc(flame->n1*flame->n2, sizeof(double));
  }
  
  /* initialize the filter for the bounds */
  double x_tmp, sumW, n_sigma = 3.0;
  flow->fhat.nW = 10;
  flow->phiRaw    = (flowBound**) calloc(3, sizeof(flowBound*));
  flow->phiHat    = (flowBound**) calloc(3, sizeof(flowBound*));
  flow->Dt_phiRaw = (flowBound**) calloc(3, sizeof(flowBound*));
  flow->Dt_phiHat = (flowBound**) calloc(3, sizeof(flowBound*));
  flow->DtHat_phi = (flowBound**) calloc(3, sizeof(flowBound*));
  flow->phiStar   = (flowBound*)  calloc(3, sizeof(flowBound));
  for (i=0; i<3; i++) {
    flow->phiRaw[i]    = (flowBound*) calloc(flow->fhat.nW, sizeof(flowBound));
    flow->phiHat[i]    = (flowBound*) calloc(flow->fhat.nW, sizeof(flowBound));
    flow->Dt_phiRaw[i] = (flowBound*) calloc(flow->fhat.nW, sizeof(flowBound));
    flow->Dt_phiHat[i] = (flowBound*) calloc(flow->fhat.nW, sizeof(flowBound));
    flow->DtHat_phi[i] = (flowBound*) calloc(flow->fhat.nW, sizeof(flowBound));
  }
    
  /* compute the filter weights */
  flow->fhat.w = (double*) calloc(flow->fhat.nW, sizeof(double));
  sumW = 0.0;
  for (j=0; j<flow->fhat.nW; j++) {
    x_tmp  = (double)(j)/n_sigma;
    flow->fhat.w[j] = 0.39894228*exp(-0.5*x_tmp*x_tmp);
    sumW += flow->fhat.w[j];
  }

  /* normalize the filter */
  for (j=0; j<flow->fhat.nW; j++) {
    flow->fhat.w[j] /= sumW;
  }
  
  /* fill the arrays */
  for (j=0; j<flow->fhat.nW; j++) {
    flow->phiRaw[1][j].min = flow->phi1min;
    flow->phiRaw[1][j].max = flow->phi1max;
    flow->phiRaw[1][j].del = flow->phi1max - flow->phi1min;
    flow->Dt_phiRaw[1][j].min = flow->Dt_phiRaw[1][j].max = 0.0;
    flow->Dt_phiRaw[1][j].del = 0.0;
    flow->Dt_phiHat[1][j].min = flow->Dt_phiHat[1][j].max = 0.0;
    flow->DtHat_phi[1][j].min = flow->DtHat_phi[1][j].max = 0.0;
    flow->Dt_phiHat[1][j].del = flow->DtHat_phi[1][j].del = 0.0;
    
    flow->phiRaw[2][j].min = flow->phi2min;
    flow->phiRaw[2][j].max = flow->phi2max;
    flow->phiRaw[2][j].del = flow->phi2max - flow->phi2min;
    flow->Dt_phiRaw[2][j].min = flow->Dt_phiRaw[2][j].max = 0.0;
    flow->Dt_phiRaw[2][j].del = 0.0;
    flow->Dt_phiHat[2][j].min = flow->Dt_phiHat[2][j].max = 0.0;
    flow->DtHat_phi[2][j].min = flow->DtHat_phi[2][j].max = 0.0;
    flow->Dt_phiHat[2][j].del = flow->DtHat_phi[2][j].del = 0.0;
  }
}
/* -------------------------------------------------------------------------- */

void initX2dimension(driveStruct *driver, ioStruct *io, chemStruct *chem,
                     flameStruct *flame, flowStruct *flow, double *Y2)
{
  int i, j, k;
  int r, row, row_pre, col, col_pre;
  int iX1, iX2, iX3, iZ2, jZ2, kZ2, jN;
  double Y_N;
  double Z1, Z2, Z2m, x1_N;
  double *Y_j, **initH, H_ref, *mod_H, *avg_H;

  // allocate local grid variables for readability
  int     n1   = flame->n1;
  int     n2   = flame->n2;
  int     nVar = flame->nVars;
  double *x1   = flame->x1;
  double *x2   = flame->x2;

  int maxDim = IMAX(n1, n2); // max number of points in 1D
  int nGrid  = n1*n2;        // total number of grid points
  int nEq    = nGrid*nVar;   // number of equations to solve

  // allocate local memory
  Y_j    = (double*) calloc(chem->nSpecSolve, sizeof(double));
  initH  = (double**) calloc(n2, sizeof(double*));
  for (i=0; i<n2; ++i) 
    initH[i] = (double*) calloc(n1, sizeof(double));
  
  // reallocate the solution vector memory
  flame->Y      = (double*) realloc(flame->Y,      nEq*sizeof(double));
  flame->Yp     = (double*) realloc(flame->Yp,     nEq*sizeof(double));
  flame->Y_m1   = (double*) realloc(flame->Y_m1,   nEq*sizeof(double));
  flame->Y_m2   = (double*) realloc(flame->Y_m2,   nEq*sizeof(double));
  flame->Y_pre  = (double*) realloc(flame->Y_pre,  nEq*sizeof(double));
  flame->sumYi  = (double*) realloc(flame->sumYi,  nEq*sizeof(double));
  flame->Yp_m1  = (double*) realloc(flame->Yp_m1,  nEq*sizeof(double));
  flame->chiA   = (double*) realloc(flame->chiA,   nGrid*sizeof(double));
  flame->rhoDotSoot = (double*) realloc(flame->rhoDotSoot, nGrid*sizeof(double));
  flame->rhoDotSpray = (double*) realloc(flame->rhoDotSpray, nGrid*sizeof(double));
  flame->xi     = (double*) realloc(flame->xi,     nGrid*sizeof(double));
  flame->D      = (double*) realloc(flame->D,     nGrid*sizeof(double));
  flame->rhoSigma = (double*) realloc(flame->rhoSigma,     nGrid*sizeof(double));
  flame->v1     = (double*) realloc(flame->v1,     nGrid*sizeof(double));
  flame->vT     = (double*) realloc(flame->vT,     nGrid*sizeof(double));
  if (chem->soot) {
     flame->vMxmyp  = (double*) realloc(flame->vMxmyp,     nGrid*sizeof(double));
     flame->vMxy = (double*) realloc(flame->vMxy,     nGrid*sizeof(double));
  }
  flame->rho    = (double*) realloc(flame->rho,    nGrid*sizeof(double));
  flame->cpmix  = (double*) realloc(flame->cpmix,  nGrid*sizeof(double));
  flame->gamma  = (double*) realloc(flame->gamma,  nGrid*sizeof(double));
  flame->mu     = (double*) realloc(flame->mu,     nGrid*sizeof(double));
  flame->D      = (double*) realloc(flame->D,      nGrid*sizeof(double));
  flame->Mw     = (double*) realloc(flame->Mw,     nGrid*sizeof(double));
  flame->HRrate = (double*) realloc(flame->HRrate, nGrid*sizeof(double));
  flame->chemSrcY    =(double*) realloc(flame->chemSrcY,nEq*sizeof(double));
  flame->chemSrcY_tmp=(double*) realloc(flame->chemSrcY_tmp,nEq*sizeof(double));
  flame->convCoeff = (double*) realloc(flame->convCoeff, nGrid*sizeof(double));
  flame->visc_src  = (double*) realloc(flame->visc_src,  nGrid*sizeof(double));
  flame->radiation  = (double*) realloc(flame->radiation,  nGrid*sizeof(double));
  flame->gPDF  = (double*) realloc(flame->gPDF,  nGrid*sizeof(double));
  if (driver->ENERGY)
    flame->T = (double*) realloc(flame->T, nGrid*sizeof(double));
  if (driver->ENERGY || io->outputHtot)
    flame->H = (double*) realloc(flame->H, nGrid*sizeof(double));
  if (driver->output_tsrc) {
    flame->conv_src  = (double*) realloc(flame->conv_src, nGrid*sizeof(double));
    flame->chem_src  = (double*) realloc(flame->chem_src, nGrid*sizeof(double));
    flame->diff_src  = (double*) realloc(flame->diff_src, nGrid*sizeof(double));
    flame->conv_temp = (double*) realloc(flame->conv_temp,nGrid*sizeof(double));
    flame->chem_temp = (double*) realloc(flame->chem_temp,nGrid*sizeof(double));
    flame->diff_temp = (double*) realloc(flame->diff_temp,nGrid*sizeof(double));
  }
  
  // allocate memory specific to 2D flamelet
  flame->chiB    = (double*) calloc(nGrid, sizeof(double));
  flame->chiC    = (double*) calloc(nGrid, sizeof(double));
  flame->v2      = (double*) calloc(nGrid, sizeof(double));
  flame->sumCp2  = (double*) calloc(nGrid, sizeof(double));
  flame->dY.y    = (double*) calloc(nVar*maxDim, sizeof(double));
  flame->dY.yy   = (double*) calloc(nVar*maxDim, sizeof(double));
  flame->dY.xy   = (double*) calloc(nVar*maxDim, sizeof(double));
  flame->dCp.y   = (double*) calloc(maxDim, sizeof(double));
  
  if (driver->H_DIM==2) {
    mod_H = (double*) calloc(n2, sizeof(double));
    avg_H = (double*) calloc(n1, sizeof(double));
  }
  else if (driver->H_DIM==1) {
    mod_H = (double*) calloc(n1, sizeof(double));
    avg_H = (double*) calloc(n2, sizeof(double));
  }
  else if (driver->H_DIM==0)
    mod_H = avg_H = NULL;

  // set all the reallocated memory to zero
  for (k=1; k < n2; k++) {
    row = k*n1*nVar;
    r   = k*n1;
    for (j=0; j < n1; j++) {
      col = j*nVar;
      
      flame->chiA[r+j] = 0.0;
      flame->v1[r+j]   = 0.0;
      
      for (i=0; i < chem->nSpecSolve; i++) {
        flame->Y[row+col+i] = flame->Yp[row+col+i] = 0.0;
        flame->Y_m1[row+col+i]  = flame->Y_m2[row+col+i]  = 0.0;
        flame->Yp_m1[row+col+i] = flame->Y_pre[row+col+i] = 0.0;
        flame->chemSrcY[row+col+i] = 0.0;
        flame->chemSrcY_tmp[row+col+i] = 0.0;
      } // end loop over species
    } // end loop over columns
  } // end loop over rows
  
  // set the indices for boundaries
  iX1 = (n1-1)*nVar;     // start of x1 = 1, x2 = 0
  iX2 = (n2-1)*n1*nVar;  //          x1 = 0, x2 = 1
  iX3 = (n1*n2-1)*nVar;  //          x1 = 1, x2 = 1
  if (driver->TRANS_Z2 < 2 || driver->H_DIM != 0) {
    iZ2 = iX2;
    jZ2 = 0;
    kZ2 = n2 - 1;
  }
  else if (driver->TRANS_Z2 == 2) {
    iZ2 = iX1 + iX2;
    jZ2 = n1 - 1;
    kZ2 = n2 - 1;
  }
  
  /* ----- initialise the solution along the Z2 axis (Z1=0)  ----- */
  // make sure the enthalpy along the Z1 axis is current
  for (j=0; j<n1; j++) {
    k = j*nVar;
    initH[0][j] = 0.0;
    chemComputeThermoData(chem->h, chem->cp, flame->Y[k+flame->sT], chem);
    for (i=0; i<chem->nSpecSolve; i++)
      initH[0][j] += chem->h[i]*flame->Y[k+i];
  }
  
  /* ----- set boundary conditions for fuel 2 (top, Z2=1) ----- */
  for (i=0; i<io->nFuel2; i++)
    flame->Y[iZ2 + io->fuelIndex2[i]] = io->fuelY2[i];
  
  initH[kZ2][jZ2] = 0.0;
  if (io->bcType[2] == 'T') // temperature specified
    {
      flame->Y[iZ2 + flame->sT] = io->fuelT2;

      // compute enthalpy
      chemComputeThermoData(chem->h, chem->cp, io->fuelT2, chem);
      for (i=0; i<io->nFuel2; i++)
        initH[kZ2][jZ2] += io->fuelY2[i]*chem->h[io->fuelIndex2[i]];
    }
  else if (io->bcType[2] == 'H') // enthalpy specified
    {
      initH[kZ2][jZ2] = io->fuelH2;

      // compute temperature
      io->fuelT2 = calcTfromH(chem->nSpecies, chem->nSpecSolve, initH[kZ2][jZ2],
                              &flame->Y[iZ2], 300, chem);
      flame->Y[iZ2 + flame->sT] = io->fuelT2;
    }
  else if (io->bcType[2] == 'L') // liquid fuel specified
    {
      // compute the enthalpy of the vapour at the reference temperature
      chemComputeThermoData(chem->h, chem->cp, io->fuelTref2, chem);
      H_ref = 0.0;
      for (i=0; i<io->nFuel2; i++)
        H_ref += io->fuelY2[i]*chem->h[io->fuelIndex2[i]];

      // compute the enthalpy of the vapour at fuel boundary
      initH[kZ2][jZ2] = io->cpLiq2*(io->fuelT2 - io->fuelTref2)-io->hVap2+H_ref;

      // compute temperature of vapour at boundary
      io->fuelT2 = calcTfromH(chem->nSpecies,chem->nSpecSolve, initH[kZ2][jZ2],
                              &flame->Y[iZ2], 300, chem);
      flame->Y[iZ2 + flame->sT] = io->fuelT2;

      printf("Fuel 2 specified as liquid, vapour temperature set to: %8.4f K\n",
             io->fuelT2);
    }
  else if (io->bcType[2] == 'E')
    {
      // initialising second axis from existing flamelet
      if (driver->H_DIM == 0)
        {
          if (driver->TRANS_Z2 < 2) // standard transform
            {
              // initialise x1 = 0 axis to existing solution
              for (k=1; k<n2; k++)
                {
                  row = k*n1*nVar;

                  // set the temperature
                  flame->Y[row + flame->sT] = Y2[k*nVar + flame->sT];

                  initH[k][0] = 0.0;
                  chemComputeThermoData(chem->h, chem->cp,
                                        flame->Y[row+flame->sT], chem);
                  for (i=0; i<chem->nSpecies; i++)
                    {
                      flame->Y[row + i] = Y2[k*nVar + i];
                      initH[k][0] += chem->h[i]*flame->Y[row + i];
                    }

                } // end loop over rows
            }
          else if (driver->TRANS_Z2 == 2) // Hasse transform
            {
              // initialise x2 = 1 axis to existing solution
              for (j=1; j<n1; j++)
                {
                  col = j*nVar;

                  flame->Y[iX2 + col + flame->sT] = Y2[col + flame->sT];

                  initH[n2-1][j] = 0.0;
                  chemComputeThermoData(chem->h, chem->cp, Y2[col + flame->sT],
                                        chem);
                  for (i=0; i<nVar; i++)
                    {
                      flame->Y[iX2 + col + i] = Y2[col + i];
                      initH[n2-1][j] += chem->h[i]*Y2[col+i];
                    }
                } // end loop over columns
            }
        } // end two-mixture formulation
    }
  else
    printf("Error: fuel 2 boundary condition specified incorrectly\n");

  if (driver->ENERGY)
    flame->Y[iZ2 + flame->sT] = initH[kZ2][jZ2];

  if (driver->H_DIM != 0)
    {
      for (i=0; i<n1; i++)
        {
          avg_H[i] = flow->avgH[i];
        }
      // Setup the min/max of modified enthalpy
      if (driver->H_DIM == 2)
        {
          mod_H[0] = flow->phi2min;
          mod_H[n2-1] = flow->phi2max;
        }
      else
        {
          mod_H[0] = flow->phi1min;
          mod_H[n1-1] = flow->phi1max;
        }
      
      // Go through and fill in avgH if it is not present in flow domain
      for (j=0; j<n1; j++)
        {
          if (avg_H[j]==0.0)
            {
              // Get mod_H for middle row
              k = (n2)/2;
              mod_H[k] = linInterp(x2[k],x2[0],mod_H[0],x2[n2-1],mod_H[n2-1]);
              
              // Get linear interpolated enthalpy
              initH[k][0] = mod_H[k] + avg_H[0];
              initH[k][n1-1] = mod_H[k] + avg_H[n1-1];
              initH[k][j] = linInterp(x1[j], x1[0], initH[k][0],
                                      x1[n1-1], initH[k][n1-1]);
              avg_H[j] = initH[k][j]-mod_H[k];
            }
        }
      
      for (k=0; k<n2; k++)
        {
          // Set the grid in between
          row = k*(n1*nVar); // start of current row
          if (k != 0)
            row_pre = (k-1)*(n1*nVar); // start of previous row
          else
            row_pre = row;
          // Get mod_H
          mod_H[k] = linInterp(x2[k], x2[0], mod_H[0], x2[n2-1], mod_H[n2-1]);
          
          // set the right boundary (Z1 = 1)
          initH[k][n1-1] = mod_H[k] + avg_H[n1-1];
          col = (n1-1)*nVar;
          for (i=0; i<chem->nSpecSolve; i++)
            flame->Y[row+col+i] = flame->Y[iX1+i];
          for (i=0; i<chem->nSpecSolve; i++)
            Y_j[i] = flame->Y[row+col+i];
          
          if (driver->ENERGY)
            flame->Y[row+col+flame->sT] = initH[k][n1-1];
          else // get the temperature based on the enthalpy
            flame->Y[row+col+flame->sT] =
              calcTfromH(chem->nSpecies, chem->nSpecSolve, initH[k][n1-1],
                         Y_j, flame->Y[row_pre+flame->sT], chem);
          
          for (j=0; j<n1-1; j++)
            {
              col = j*nVar; // current row/column position
              initH[k][j] = mod_H[k] + avg_H[j];
              
              // linear interpolate species
              for (i=0; i<chem->nSpecSolve; i++)
                flame->Y[row+col+i] = linInterp(x1[j], x1[0], flame->Y[i],
                                                x1[n1-1], flame->Y[iX1 + i]);
              
              if (driver->ENERGY)
                flame->Y[row+col+flame->sT] = initH[k][j];
              else // get the temperature based on the enthalpy
                flame->Y[row+col+flame->sT] =
                  calcTfromH(chem->nSpecies,chem->nSpecSolve,initH[k][j],
                             &flame->Y[row+col],flame->Y[row_pre + flame->sT],
                             chem);
            } // end loop over columns
        } // end loop over rows
      free(avg_H);
    }
  else if (driver->TRANS_Z2 < 2) // standard transform
    {
      if (driver->TRANS_Z2 == 2)
        {
          initH[n2-1][0] = initH[0][0];
          for (i=0; i<nVar; i++)
            flame->Y[iX2+i] = flame->Y[i];
        }
      /* ----- loop over all the interior points ----- */
      for (k=1; k<n2-1; k++)
        {
          row = k*(n1*nVar); // start of current row
          row_pre = (k-1)*(n1*nVar); // start of previous row
          
          /* ----- set the left boundary ----- */
          if (driver->TRANS_Z2 == 2)
            {
              /* for this transformation, left boundary (Z1 = 0) is a 
                 singularity. Therefore, set all values to (Z1,Z2) = (0,0) */
              initH[k][0] = initH[0][0];
              for (i=0; i<nVar; i++)
                flame->Y[row+i] = flame->Y[i];
            }
          else if (io->bcType[2] != 'E') /* x1 = 0 --> Z1 = 0 */
            { /* only set if initialising from b.c.s */
              // interpolate enthalpy
              initH[k][0] = linInterp(x2[k], x2[0], initH[0][0],
                                      x2[n2-1], initH[n2-1][0]);

              // interpolate species
              for (i=0; i<chem->nSpecSolve; i++)
                flame->Y[row+i] = linInterp(x2[k], x2[0], flame->Y[i],
                                            x2[n2-1], flame->Y[iX2 + i]);
	  
              // set the energy variable
              if (driver->ENERGY)
                flame->Y[row+flame->sT] = initH[k][0];
              else // compute temperature from enthalpy
                flame->Y[row+flame->sT] =
                  calcTfromH(chem->nSpecies, chem->nSpecSolve, initH[k][0],
                             &flame->Y[row], flame->Y[row_pre+flame->sT], chem);
            } // end computation of left boundary
      
          /* ----- loop over interior points in x1 direction ----- */
          for (j=1; j<n1-1; j++)
            {
              col = j*nVar; // current row/column position
              
              // compute the actual Z2 for this point
              if (driver->TRANS_Z2 == 0)
                {
                  Z1 = x1[j];
                  Z2 = x2[k];
                }
              else if (driver->TRANS_Z2 == 1)
                {
                  Z1 = x1[j];
                  Z2 = x2[k]*(1.0-x1[j]);
                }
              else if (driver->TRANS_Z2 == 2)
                {
                  Z1 = x1[j]*(1.0-x2[k]);
                  Z2 = x1[j]*x2[k];
                }
              
              // intersection of the mixing line with the x1 axis
              x1_N = Z1/(1-Z2);
              
              if (Z1 + Z2 > 1.0)
                {
                  /* This should only occur on an untransformed mesh.
                     Set to interpolation between fuel boundaries so that
                     the Cp calc is not crazy */
                  Z2m = 0.5*(1.0+Z2-Z1);
                  initH[k][j] = Z2m*initH[n2-1][0] + (1-Z2m)*initH[0][n1-1];
                  for (i=0; i<chem->nSpecSolve; i++)
                    flame->Y[row+col+i] =       Z2m*flame->Y[iX2+i]
                      + (1-Z2m)*flame->Y[iX1+i];
                }
              else /* standard interpolation along mixing line */
                {
                  // find the location of the mesh point before the actual
                  jN = gridLocate(n1, x1, x1_N);
                  
                  // interpolate enthalpy at x1_N between neighbours
                  Y_N = linInterp(x1_N, x1[jN], initH[0][jN],
                                  x1[jN+1],initH[0][jN+1]);
                  initH[k][j] = Z2*initH[kZ2][jZ2] + (1-Z2)*Y_N;
                  
                  // loop over all species
                  for (i=0; i<chem->nSpecSolve; i++)
                    {
                      // interpolate along Z2=0 boundary to get current value
                      Y_N = linInterp(x1_N, x1[jN], flame->Y[jN*nVar + i],
                                      x1[jN+1], flame->Y[(jN+1)*nVar + i]);
                      
                      // interpolate along mixing line from Z2 = 1 to Z1 = x_N
                      flame->Y[row+col+i] =  Z2*flame->Y[(kZ2*n1+jZ2)*nVar+i]
                        + (1-Z2)*Y_N;
                    }
                }
              
              if (driver->ENERGY)
                flame->Y[row+col+flame->sT] = initH[k][j];
              else // get the temperature based on the enthalpy
                flame->Y[row+col+flame->sT] = 
                  calcTfromH(chem->nSpecies, chem->nSpecSolve, initH[k][j],
                             &flame->Y[row+col],flame->Y[row_pre+col+flame->sT],
                             chem);
            } // end loop over columns
	  
          /* ----- set the right boundary ----- */
          col = (n1-1)*nVar;
          if (driver->TRANS_Z2 == 0)
            { // interpolate between Z1+Z2=0.5 and Z1=1
              initH[k][n1-1] =       0.5*x2[k] *initH[n2-1][0]
                + (1 - 0.5*x2[k])*initH[0][n1-1];
              for (i=0; i<chem->nSpecSolve; i++)
                flame->Y[row+col+i] =       0.5*x2[k]* flame->Y[iX2+i]
                  + (1 - 0.5*x2[k])*flame->Y[iX1+i];
            }
          else if (driver->TRANS_Z2 == 1)
            { // right boundary is singular: set all to value at (1,0)
              initH[k][n1-1] = initH[0][n1-1];
              for (i=0; i<chem->nSpecSolve; i++)
                flame->Y[row+col+i] = flame->Y[iX1+i];
            }
          else if (driver->TRANS_Z2 == 2)
            { // Z1+Z2=1 boundary
              initH[k][n1-1] =     x2[k] *initH[n2-1][n1-1]
                + (1-x2[k])*initH[0][n1-1];
              for (i=0; i<flame->nVars; i++)
                flame->Y[row+col+i] =      x2[k] *flame->Y[iZ2+i]
                  + (1-x2[k])*flame->Y[iX1+i];
            }
          if (driver->ENERGY)
            flame->Y[row+col+flame->sT] = initH[k][n1-1];
          else
            flame->Y[row+col+flame->sT] = 
              calcTfromH(chem->nSpecies, chem->nSpecSolve, initH[k][n1-1],
                         &flame->Y[row+col],flame->Y[row_pre+col+flame->sT],
                         chem);
          
        } /* end loop over rows (x2) */
      
      /* ----- set the top boundary ----- */
      for (j=0; j<n1; j++)
        {
          row_pre = (n2-2)*n1*nVar;
          col = j*nVar;     // current column

          if (driver->TRANS_Z2 == 0 && j != 0)
            {
              // interpolate between Z1+Z2=0.5 and Z1=1
              initH[n2-1][j] = (1 - 0.5*x1[j])*initH[n2-1][0]
                                  + 0.5*x1[j] *initH[0][n1-1];
              for (i=0; i<chem->nSpecSolve; i++)
                flame->Y[iX2+col+i] = (1 - 0.5*x1[j])*flame->Y[iX2+i]
                                        +  0.5*x1[j] *flame->Y[iX1+i];
            }
          else if (driver->TRANS_Z2 == 1 && j != 0)
            { // boundary: Z1+Z2 = 1
              Z2  = 1 - x1[j];  // current Z2 based on Z1+Z2=1              
	  
              // interpolate the enthalpy along the boundary
              initH[n2-1][j] = Z2*initH[n2-1][0] + (1-Z2)*initH[0][n1-1];
	  
              // interpolate species mass fractions
              for (i=0; i<chem->nSpecSolve; i++)
                flame->Y[iX2+col+i] = Z2*flame->Y[iX2+i]+(1-Z2)*flame->Y[iX1+i];
            }
          else if (driver->TRANS_Z2 == 2 && io->bcType[2] != 'E' && j != n1-1)
            { // boundary: Z1 = 0
              initH[n2-1][j] = linInterp(x1[j], x1[0], initH[n2-1][0],
                                         x1[n1-1], initH[n2-1][n1-1]);
              for (i=0; i<chem->nSpecSolve; i++)
                flame->Y[iX2+col+i] = linInterp(x1[j], x1[0], flame->Y[iX2+i],
                                                x1[n1-1], flame->Y[iZ2+i]);
            }

          if (driver->ENERGY)
            flame->Y[iX2+col+flame->sT] = initH[n2-1][j];
          else // get the temperature based on the enthalpy
            flame->Y[iX2+col+flame->sT] =
              calcTfromH(chem->nSpecies,chem->nSpecSolve,initH[n2-1][j],
                         &flame->Y[iX2+col],flame->Y[row_pre+col+flame->sT],
                         chem);
        }
    } // end standard transform
  else if (driver->TRANS_Z2 == 2) // Hasse transform
    {
      // set oxidizer boundary
      for (k=1; k<n2; k++)
        {
          row = k*n1*nVar;
          for(i=0; i<nVar; i++)
            flame->Y[row+i] = flame->Y[i];
          initH[k][0] = initH[0][0];
        } // end loop over rows

      /* ----- loop over all interior points ----- */
      for (j=1; j<n1; j++)
        {
          col     = j*nVar;
          col_pre = (j-1)*nVar;

          // set the top boundary (x2 = 1) if initialising from b.c.s
          if (io->bcType[2] != 'E')
            {
              initH[kZ2][j] = linInterp(x1[j], x1[0], initH[kZ2][0],
                                        x1[n1-1], initH[kZ2][n1-1]);

              for (i=0; i<chem->nSpecSolve; i++)
                flame->Y[iX2+col+i] = linInterp(x1[j], x1[0], flame->Y[iX2+i],
                                                x1[n1-1], flame->Y[iZ2+i]);
	  
              if (driver->ENERGY)
                flame->Y[iX2+col+flame->sT] = initH[kZ2][j];
              else // get the temperature based on the enthalpy
                flame->Y[iX2+col+flame->sT] =
                  calcTfromH(chem->nSpecies, chem->nSpecSolve, initH[kZ2][j],
                             &flame->Y[iX2+col],flame->Y[iX2+col_pre+flame->sT],
                             chem);
            }

          // loop over interior row points
          for (k=1; k<n2-1; k++)
            {
              row = k*n1*nVar;

              // for now, linearly interpolate along constant Z (i.e. Z1+Z2)
              initH[k][j] = linInterp(x2[k], x2[0], initH[0][j],
                                      x2[kZ2], initH[kZ2][j]);

              for (i=0; i<chem->nSpecSolve; i++)
                flame->Y[row+col+i] = linInterp(x2[k], x2[0], flame->Y[col+i],
                                                x2[kZ2], flame->Y[iX2+col+i]);

              if (driver->ENERGY)
                flame->Y[row+col+flame->sT] = initH[k][j];
              else
                flame->Y[row+col+flame->sT] =
                  calcTfromH(chem->nSpecies,chem->nSpecSolve,initH[k][j],
                             &flame->Y[row+col],flame->Y[row+col_pre+flame->sT],
                             chem);

            } // end loop over rows
        } // end loop over columns
    } // end Hasse transform

 for (i=0; i<n2; ++i) free(initH[i]);
 free(initH);
 free(Y_j);
}
/* -------------------------------------------------------------------------- */

double MDUCcalcSpaldingFunc(double Ts, void *args)
{
  ioStruct *io = (ioStruct*) args;

  double BT, Bf, Psat, Q, cpG, cpL, hfg, Xfs, Yfs;

  /* get the thermo tabulated from NIST */
  NIST_getThermo(io->satData, &hfg, &cpG, &cpL, &Psat, Ts);

  /* compute mass fraction at the surface based on saturation pressure*/
  Xfs = fmin(1.0,Psat/io->oxP); // vars[1] == P_inf

  Yfs = io->Mfuel1*Xfs/(io->Mfuel1*Xfs+(1-Xfs)*io->Mox);

  Q = hfg + cpL*(Ts - io->fuelT1);   // vars[0] == T_liq

  BT = cpG*(io->oxT - Ts)/Q;        // vars[2] == T_inf
  if (Yfs < 1.0)
    Bf = (Yfs - 0.0)/(1.0 - Yfs); // vars[3] == Yf_inf
  else
    Bf = 1e9; // big value

  return BT - Bf;
}
/* -------------------------------------------------------------------------- */

