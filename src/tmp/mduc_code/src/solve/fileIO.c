#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <math.h>
#include "fileIO.h"
#include "parser.h"
#include "accessChem.h" // should remove dependence on this
#include "chemMech.h"
#include "utils.h"
#include "errFunc.h"

void createDataFileName(const char *fileName, ioStruct *io, timeStruct *time,
                        char *dataFileName)
{
  int i, j;
  char iStr[10], numString[50]; /* for number string */
  char file[256], *head, *tail;  /* strings to hold beginning and end */
  
  /* read in the leading file name and extension based on the first */
  strcpy(file, fileName);
  head = strtok(file, ".");
  tail = strtok(NULL, " ");
  
  if (!strcmp(io->numStyle, "flow-time")) {
    /* convert the time to a string */
    sprintf(numString, "%0*.*fms", io->fileNumPrec+5, io->fileNumPrec,
            time->t*1e3);
    /* replace the period with an underscore */
    numString[strcspn(numString,".")] = '_';
  }
  else if (!strcmp(io->numStyle, "time-step")) {
    /* conver the time-step to a string */
    sprintf(numString, "%0*i", io->fileNumPrec, time->nt);
  }
  else if ( !strcmp(io->numStyle, "crank-angle") ||
            !strcmp(io->numStyle, "CAD") ) {
    double CAD     = time->t*6*time->RPM + time->CAD_0;
    double CAD_rel = fmod(CAD, 720.0);
    int    cycle   = (int)(CAD)/720;
    /* conver the time-step to a string */
    sprintf(numString, "cyc%02i_%04.1fCAD", cycle, CAD_rel);
    for (i=0; i<strlen(numString); i++)
      if (numString[i] == '.')
        numString[i] = 'p';
  }
  
  /* compile new file name */
  sprintf(dataFileName, "%s%s.%s", head, numString, tail);
}
/* -------------------------------------------------------------------------- */

void steadyFlameletFileName(const char *fileName, ioStruct *io,flowStruct *flow,
                            char *dataFileName)
{
  int i, j;
  char iStr[10], numString[50]; // for number string
  char file[256], *head, *tail;  // strings to hold beginning and end
  
  // read in the leading file name and extension based on the first .
  strcpy(file, fileName);
  head = strtok(file, ".");
  tail = strtok(NULL, " ");
  
  /* fprintf(stderr, "P: %03.0f\n", flow->P*1e-5); */
  /* fprintf(stderr, "P: %05.0f\n", ); */

  sprintf(dataFileName, "%sP%03.0f_%05.0f.%s", head, flow->P*1e-5,
          remainder(flow->P, 1.0e5), tail);
}
/* -------------------------------------------------------------------------- */

void checkMonitor(ioStruct *io)
{
  /* exit if monitor stream already exists */
  if (io->monitor)
    return;

  if (!strcmp(io->monitorFile, "stderr"))
    io->monitor = stderr;
  else
    io->monitor = fopen(io->monitorFile, "w");
}
/* -------------------------------------------------------------------------- */

int  checkFileExists(char *fileName)
{
  int exists = 0;
  FILE *pF = fopen(fileName, "r");
  if (pF) {
    exists = 1;
    fclose(pF);
  }
  return exists;
}
/* -------------------------------------------------------------------------- */

void readInitFile(ioStruct *io, driveStruct *driver, timeStruct *time,
                  flameStruct *flame, flowStruct *flow, sdialStruct *sdial,
                  engineStruct *engine)
{
  int   i, j;
  int   readFlame; // flag for ending the reading of a flamelet

  // strings to store the tag for each line and comments (if any)
  char  tag[50], comment[200], pdfTag[50];
  
  FILE   *pFile;   // input stream
  pFile = fopen(io->inputFile, "r"); // open input file

  if (!pFile) {
    printf("ERROR: cannot open input file %s\n", io->inputFile);
    exit(1);
  }
  
  /* ----- loop until end of file is reached  ----- */
  while (!feof(pFile)) { 
    fscanf(pFile, "%s", tag); // read in the first string
    
    // first check to see if end of input is desired
    if (!strcmp(tag, "ENDINPUT"))
      break;
    // check to see if it is a comment
    else if (tag[0] == '!' || tag[0] == '#')
      fgets(comment, 200, pFile);
    
    /* ----- driver variables ----- */
    else if (!strcmp(tag, "number-of-flamelets"))
      driver->nFlamelets = readInteger(pFile);
    else if (!strcmp(tag, "pressure-mode"))
      readString(pFile, io->pressureMode);
    else if (!strcmp(tag, "mixing-mode"))
      driver->MIX_MODE  = readInteger(pFile);
    else if (!strcmp(tag, "chi-mode"))
      readString(pFile, io->chiMode);
    else if (!strcmp(tag, "chi-interp"))
      readString(pFile, io->chiInterp);
    else if (!strcmp(tag, "chi-scale"))
      driver->CHI_SCALE = readInteger(pFile);
    else if (!strcmp(tag, "pdf-mode")) {
      readString(pFile,pdfTag);
      if (!strcmp(pdfTag, "beta") || !strcmp(pdfTag, "beta-pdf")) 
	io->pdfMode==PDF_TYPE_BETA;
      if (!strcmp(pdfTag, "sml"))
	io->pdfMode==PDF_TYPE_SML;
      if (!strcmp(pdfTag, "gaussian") || !strcmp(pdfTag, "gauss")) 
	io->pdfMode==PDF_TYPE_GAUSS;
      if (!strcmp(pdfTag, "lognormal") || !strcmp(pdfTag, "lognorm")) 
	io->pdfMode==PDF_TYPE_LOG;
      if (!strcmp(pdfTag, "dirichlet") || !strcmp(pdfTag, "diric")) 
	io->pdfMode==PDF_TYPE_DIRICHLET;
      if (!strcmp(pdfTag, "joint-beta") || !strcmp(pdfTag, "jbeta")) 
	io->pdfMode==PDF_TYPE_JBETA;
      if (!strcmp(pdfTag, "independent") || !strcmp(pdfTag, "independent-beta")) 
	io->pdfMode==PDF_TYPE_INDEPENDENT_BETA;
      if (!strcmp(pdfTag, "independent-sml")) 
	io->pdfMode==PDF_TYPE_INDEPENDENT_SML;
      if (!strcmp(pdfTag, "from-file")) 
	io->pdfMode==PDF_TYPE_FROM_FILE;

    }
    else if (!strcmp(tag, "time-mode"))
      readString(pFile, io->timeMode);
    else if (!strcmp(tag, "diffusion-mode") || !strcmp(tag, "diff-mode"))
      readOptArray(pFile, 3, driver->DIFF_MODE);
    else if (!strcmp(tag, "convection-mode"))
      driver->CONV_MODE = readInteger(pFile);
    else if (!strcmp(tag, "enthalpy-dimension"))
      driver->H_DIM     = readInteger(pFile);
    else if (!strcmp(tag, "read-data"))
      io->readData = readInteger(pFile);
    else if (!strcmp(tag, "read-rates"))
      io->readRates = readInteger(pFile);
    else if (!strcmp(tag, "output-temp-src"))
      driver->output_tsrc = readReal(pFile);
    else if (!strcmp(tag, "heat-loss"))
      readString(pFile, io->heatLoss);
    else if (!strcmp(tag, "dns-mode"))
      driver->DNS = readInteger(pFile);
    else if (!strcmp(tag, "bounds-filter"))
      readString(pFile, io->boundFilter);

    /* ----- solver settings ----- */
    else if (!strcmp(tag, "sundials-solver")) {
      readString(pFile, sdial->solver);
      strcpy(io->sdialSolver, sdial->solver);
    }
    else if (!strcmp(tag, "sundials-method"))
      readString(pFile, sdial->method);
    else if (!strcmp(tag, "sundials-preconditioner"))
      readString(pFile, sdial->precnd);
    else if (!strcmp(tag, "spatial-order"))
      flame->order  = readInteger(pFile);
    else if (!strcmp(tag, "energy-equation"))
      driver->ENERGY = readInteger(pFile);
    else if (!strcmp(tag, "temperature-tolerance"))
      driver->TABSTOL = readReal(pFile);
    else if (!strcmp(tag, "soot-tolerance"))
      driver->SABSTOL = readReal(pFile);
    else if (!strcmp(tag, "absolute-tolerance"))
      driver->ABSTOL = readReal(pFile);
    else if (!strcmp(tag, "relative-tolerance"))
      driver->RELTOL = readReal(pFile);
    else if (!strcmp(tag, "quadrature-abs-tolerance"))
      driver->QABSTOL = readReal(pFile);
    else if (!strcmp(tag, "quadrature-rel-tolerance"))
      driver->QRELTOL = readReal(pFile);
    else if (!strcmp(tag, "viscosity-update"))
      readString(pFile, io->viscUpdate);
    else if (!strcmp(tag, "thermophoresis"))
      driver->thermophoresis = readInteger(pFile);
    else if (!strcmp(tag, "soot-diffusion"))
      driver->sootDiff = readInteger(pFile);
    else if (!strcmp(tag, "soot-diffusivity"))
      driver->sootDiffusivity = readReal(pFile);
    else if (!strcmp(tag, "detailed-diffusion"))
      driver->detailedDiff = readInteger(pFile);
    else if (!strcmp(tag, "radiation"))
      driver->radiation = readInteger(pFile);
    else if (!strcmp(tag, "ignite"))
      driver->ignite = readInteger(pFile);
    
    /* ----- possible time variables ----- */
    else if (!strcmp(tag, "start-time"))
      time->t0 = readReal(pFile);
    else if (!strcmp(tag, "time-step-number"))
      time->nt = readInteger(pFile);
    else if (!strcmp(tag, "time-step-size"))
      time->dt = readReal(pFile);
    else if (!strcmp(tag, "final-time"))
      time->t_end = readReal(pFile);
    else if (!strcmp(tag, "time-steps"))
      time->nIter = readInteger(pFile);
    else if (!strcmp(tag, "CFL"))
      time->CFL = readReal(pFile);
    
    /* ----- IO variables ----- */
    else if (!strcmp(tag, "title"))
      readString(pFile, io->title);
    else if (!strcmp(tag, "author"))
      readString(pFile, io->author);
    else if (!strcmp(tag, "date"))
      readString(pFile, io->date);
    else if (!strcmp(tag, "data-directory"))
      readString(pFile, io->dataDir);
    else if (!strcmp(tag, "read-data-file")) {
      readString(pFile, io->inputData);
      io->readData = 1;
    }
    else if (!strcmp(tag, "read-rates-file")) {
      readString(pFile, io->inputRates);
      io->readRates = 1;
    }
    else if (!strcmp(tag, "chi-data-directory")) 
      readString(pFile, io->chiDataDir);
    else if (!strcmp(tag, "read-flow-file")) {
      readString(pFile, io->inputFlow);
      io->readFlow = 1;
    }
    else if (!strcmp(tag, "write-data-file"))
      readString(pFile, io->outputData);
    else if (!strcmp(tag, "flow-data-file")) {
      readString(pFile, io->flowData);
      io->writeFlow = 1;
    }
    else if (!strcmp(tag, "write-frequency"))
      io->writeFreq = readInteger(pFile);
    else if (!strcmp(tag, "print-frequency"))
      io->printFreq = readInteger(pFile);
    else if (!strncmp(tag, "output-var", 10))
      io->nOutputVar = readNames(pFile, &io->outputVar);
    else if (!strcmp(tag, "output-visc-src"))
      io->outputViscSrc = readInteger(pFile);
    else if (!strcmp(tag, "output-global-pdf"))
      io->outputGPDF = readInteger(pFile);
    else if (!strcmp(tag, "output-total-enthalpy"))
      io->outputHtot = readInteger(pFile);
    else if (!strcmp(tag, "output-heat-release"))
      io->outputHR = readInteger(pFile);
    else if (!strcmp(tag, "output-chi"))
      io->outputChi = readInteger(pFile);
    else if (!strcmp(tag, "output-w"))
      driver->output_w = readInteger(pFile);
    else if (!strcmp(tag, "output-rho"))
      io->outputRho = readInteger(pFile);
    else if (!strcmp(tag, "output-rhoSigma"))
      io->outputRhoSigma = readInteger(pFile);
    else if (!strcmp(tag, "output-pdf"))
      io->outputPDF = readInteger(pFile);
    else if (!strcmp(tag, "output-cfl"))
      io->outputCFL = readInteger(pFile);
    else if (!strcmp(tag, "output-dry-nox"))
      io->outputDryNOX = readInteger(pFile);
    else if (!strcmp(tag, "convection-coefficient"))
      io->outputConv = readReal(pFile);
    else if (!strcmp(tag, "post-process-file") || !strcmp(tag, "flow-prefix"))
      readString(pFile, io->postFile);
    else if (!strcmp(tag, "file-number-style"))
      readString(pFile, io->numStyle);
    else if (!strcmp(tag, "file-number-precision"))
      io->fileNumPrec = readInteger(pFile);
    else if (!strcmp(tag, "file-number-start"))
      io->fileMin = readInteger(pFile);
    else if (!strcmp(tag, "file-number-stop"))
      io->fileMax = readInteger(pFile);
    else if (!strcmp(tag, "monitor-file"))
      readString(pFile, io->monitorFile);
    else if (!strcmp(tag, "error-file"))
      readString(pFile, io->errorFile);
    else if (!strcmp(tag, "ensight-directory")) {
      io->writeEnsight = 1;
      readString(pFile, io->ensightDir);
    }

    /* ----- grid drivers ----- */
    else if (!strcmp(tag, "MaxGridPoints"))
      driver->maxNum1 = readInteger(pFile);
    else if (!strcmp(tag, "MaxGridPoints2"))
      driver->maxNum2 = readInteger(pFile);
    else if (!strcmp(tag, "Z2-transformation"))
      readString(pFile, io->transZ2);
    
    /* ----- reference scalars ----- */
    else if (!strcmp(tag, "pressure"))
      flow->P = readReal(pFile);
    else if (!strcmp(tag, "Z1-max")) {
      flow->chiZ1max = readReal(pFile);
    }
    else if (!strcmp(tag, "Z-ref"))
      readOptArray(pFile, 3, flow->ZRef);
    else if (!strcmp(tag, "Z0-ref"))
      flow->ZRef[0] = readReal(pFile);
    else if (!strcmp(tag, "Z1-ref"))
      flow->ZRef[1] = readReal(pFile);
    else if (!strcmp(tag, "Z2-ref"))
      flow->ZRef[2] = readReal(pFile);
    else if (!strcmp(tag, "Z-mean"))
      readOptArray(pFile, 3, flow->phiMean);
    else if (!strcmp(tag, "Z0-mean"))
      flow->phiMean[0] = readReal(pFile);
    else if (!strcmp(tag, "Z1-mean"))
      flow->phiMean[1] = readReal(pFile);
    else if (!strcmp(tag, "Z2-mean"))
      flow->phiMean[2] = readReal(pFile);
    else if (!strcmp(tag, "Z-variance"))
      readOptArray(pFile, 3, flow->phiVar);
    else if (!strncmp(tag, "Z0-variance", 6))
      flow->phiVar[0] = readReal(pFile);
    else if (!strncmp(tag, "Z1-variance", 6))
      flow->phiVar[1] = readReal(pFile);
    else if (!strncmp(tag, "Z2-variance", 6))
      flow->phiVar[2] = readReal(pFile);
    else if (!strcmp(tag, "chi-ref"))
      readOptArray(pFile, 3, flow->chiRef);
    else if (!strcmp(tag, "chi0-ref"))
      flow->chiRef[0] = readReal(pFile);
    else if (!strcmp(tag, "chi1-ref"))
      flow->chiRef[1] = readReal(pFile);
    else if (!strcmp(tag, "chi2-ref"))
      flow->chiRef[2] = readReal(pFile);
    else if (!strcmp(tag, "strain-rate"))
      readOptArray(pFile, 3, flow->strain);
    else if (!strcmp(tag, "strain-rate-0"))
      flow->strain[0] = readReal(pFile);
    else if (!strcmp(tag, "strain-rate-1"))
      flow->strain[1] = readReal(pFile);
    else if (!strcmp(tag, "strain-rate-2"))
      flow->strain[2] = readReal(pFile);
    else if (!strcmp(tag, "Min-enthalpy"))
      io->oxH = readReal(pFile);
    else if (!strcmp(tag, "Max-enthalpy"))
      io->fuelH1 = readReal(pFile);
    
    /* ----- chemistry drivers ----- */
    else if (!strcmp(tag, "soot-output-src"))
      io->outputSootSrc = readInteger(pFile);
    else if (!strcmp(tag, "soot-outout-Mxmyp"))
      io->outputSootMxmyp = readInteger(pFile);
    else if (!strcmp(tag, "output-conv-vel"))
      io->outputConvVel = readInteger(pFile);
    else if (!strcmp(tag, "dump-error-file"))
      io->dumpErrorFile = readInteger(pFile);
      
    /* ----- engine setup ----- */
    else if (!strcmp(tag, "engine-setup-file"))
      readString(pFile, io->engineFile);
    else if (!strcmp(tag, "engine-type")) {
      readString(pFile, engine->type);
      strcpy(io->engineType, engine->type);
    }
    else if (!strcmp(tag, "engine-mass"))
      engine->mass = readReal(pFile);
    else if (!strcmp(tag, "engine-bore"))
      engine->bore = readReal(pFile);
    else if (!strcmp(tag, "engine-stroke"))
      engine->stroke = readReal(pFile);
    else if (!strcmp(tag, "engine-rod-length"))
      engine->rod = readReal(pFile);
    else if (!strcmp(tag, "engine-comp-ratio"))
      engine->CR = readReal(pFile);
    else if (!strcmp(tag, "engine-volume-tdc"))
      engine->V_TDC = readReal(pFile);
    else if (!strcmp(tag, "engine-speed"))
      time->RPM = engine->RPM = readReal(pFile);
    else if (!strcmp(tag, "engine-start-cad"))
      time->CAD_0 = engine->CAD_0 = readReal(pFile);
    else if (!strcmp(tag, "engine-stop-cad"))
      engine->CAD_end = readReal(pFile);
    else if (!strcmp(tag, "engine-step-size"))
      engine->Dt_CAD = readReal(pFile);
    else if (!strcmp(tag, "engine-manifold-pressure"))
      readString(pFile, engine->manifoldFile);
    /* intake valve */
    else if (!strcmp(tag, "intake-valve-open"))
      engine->IVO = readReal(pFile);
    else if (!strcmp(tag, "intake-valve-close"))
      engine->IVC = readReal(pFile);
    else if (!strcmp(tag, "intake-port-diameter"))
      engine->intake.D_port = readReal(pFile);
    else if (!strcmp(tag, "intake-stem-diameter"))
      engine->intake.D_stem = readReal(pFile);
    else if (!strcmp(tag, "intake-head-diameter"))
      engine->intake.D_head = readReal(pFile);
    else if (!strcmp(tag, "intake-seat-diameter"))
      engine->intake.D_seat = readReal(pFile);
    else if (!strcmp(tag, "intake-seat-width"))
      engine->intake.W = readReal(pFile);
    else if (!strcmp(tag, "intake-seat-angle"))
      engine->intake.theta = readReal(pFile);
    else if (!strcmp(tag, "intake-valve-lift"))
      readString(pFile, engine->intake.liftFile);
    else if (!strcmp(tag, "intake-manifold-pressure"))
      readString(pFile, engine->intake.P_file);
    else if (!strcmp(tag, "intake-mass-flow"))
      readString(pFile, engine->intake.mdotFile);
/*     else if (!strcmp(tag, "intake-temperature")) /\* redundant? *\/ */
/*       readString(pFile, engine->intake.P_file); */
    /* exhaust valve */
    else if (!strcmp(tag, "exhaust-valve-open"))
      engine->EVO = readReal(pFile);
    else if (!strcmp(tag, "exhaust-valve-close"))
      engine->EVC = readReal(pFile);
    else if (!strcmp(tag, "exhaust-port-diameter"))
      engine->exhaust.D_port = readReal(pFile);
    else if (!strcmp(tag, "exhaust-stem-diameter"))
      engine->exhaust.D_stem = readReal(pFile);
    else if (!strcmp(tag, "exhaust-head-diameter"))
      engine->exhaust.D_head = readReal(pFile);
    else if (!strcmp(tag, "exhaust-seat-diameter"))
      engine->exhaust.D_seat = readReal(pFile);
    else if (!strcmp(tag, "exhaust-seat-width"))
      engine->exhaust.W = readReal(pFile);
    else if (!strcmp(tag, "exhaust-seat-angle"))
      engine->exhaust.theta = readReal(pFile);
    else if (!strcmp(tag, "exhaust-valve-lift"))
      readString(pFile, engine->exhaust.liftFile);
    else if (!strcmp(tag, "exhaust-manifold-pressure"))
      readString(pFile, engine->exhaust.P_file);
    else if (!strcmp(tag, "exhaust-mass-flow"))
      readString(pFile, engine->exhaust.mdotFile);
    /* injection */
    else if (!strcmp(tag, "injection-profile"))
      engine->nInj = readNames(pFile, &engine->injProf);
    else if (!strcmp(tag, "injection-timing"))
      readRealArraySize(pFile, &engine->inj_t);
    else if (!strcmp(tag, "injection-duration"))
      readRealArraySize(pFile, &engine->injDur);
    /* spark options */
    else if (!strcmp(tag, "engine-spark-model"))
      readString(pFile, engine->sparkModel);
    else if (!strcmp(tag, "engine-spark-advance"))
      engine->sparkAdv = readReal(pFile);
    else if (!strncmp(tag, "engine-spark-product-species", 25))
      engine->nPremSp = readNames(pFile, &engine->premSp);
    else if (!strncmp(tag, "engine-spark-product-frac", 25))
      readRealArraySize(pFile, &engine->Yburnt);
    else if (!strncmp(tag, "engine-spark-diameter", 17))
      engine->sparkD = readReal(pFile);
    /* wiebe */
    else if (!strcmp(tag, "wiebe-a"))
      engine->wiebeHCCI.a = readReal(pFile);
    else if (!strcmp(tag, "wiebe-m"))
      engine->wiebeHCCI.m = readReal(pFile);
    else if (!strcmp(tag, "wiebe-soc"))
      engine->wiebeHCCI.t0 = readReal(pFile);
    else if (!strcmp(tag, "wiebe-alpha"))
      engine->wiebeHCCI.alpha = readReal(pFile);
    else if (!strcmp(tag, "wiebe-beta"))
      engine->wiebeHCCI.beta = readReal(pFile);
    else if (!strcmp(tag, "wiebe-spark-a"))
      engine->wiebeSI.a = readReal(pFile);
    else if (!strcmp(tag, "wiebe-spark-m"))
      engine->wiebeSI.m = readReal(pFile);
    else if (!strcmp(tag, "wiebe-spark-soc"))
      engine->wiebeSI.t0 = readReal(pFile);
    else if (!strcmp(tag, "wiebe-spark-alpha"))
      engine->wiebeSI.alpha = readReal(pFile);
    else if (!strcmp(tag, "wiebe-spark-beta"))
      engine->wiebeSI.beta = readReal(pFile);
    else if (!strcmp(tag, "engine-burn-rate"))
      readString(pFile, engine->burnRateFile);
    else if (!strcmp(tag, "burn-rate-scale"))
      engine->burnRateScale = readReal(pFile);

    /* heat transfer */
    else if (!strcmp(tag, "engine-wall-temp"))
      engine->T_wall = readReal(pFile);
    else if (!strcmp(tag, "woschni-a"))
      engine->woschni.a = readReal(pFile);
    else if (!strcmp(tag, "woschni-b"))
      engine->woschni.b = readReal(pFile);
    else if (!strcmp(tag, "woschni-c"))
      engine->woschni.c = readReal(pFile);
    else if (!strcmp(tag, "woschni-d"))
      engine->woschni.d = readReal(pFile);
    else if (!strcmp(tag, "woschni-c1"))
      engine->woschni.C2 = readReal(pFile);
    else if (!strcmp(tag, "woschni-c2"))
      engine->woschni.C2 = readReal(pFile);
    else if (!strcmp(tag, "woschni-alpha"))
      engine->woschni.alpha = readReal(pFile);
    else if (!strcmp(tag, "woschni-ref-temp"))
      engine->woschni.Tref = readReal(pFile);
    else if (!strcmp(tag, "woschni-ref-pressure"))
      engine->woschni.Pref = readReal(pFile);
    else if (!strcmp(tag, "woschni-pressure-scale"))
      engine->woschni.Pscale = readReal(pFile);
    else if (!strcmp(tag, "hohen-a"))
      engine->hohen.a = readReal(pFile);
    else if (!strcmp(tag, "hohen-b"))
      engine->hohen.b = readReal(pFile);
    else if (!strcmp(tag, "hohen-c"))
      engine->hohen.c = readReal(pFile);
    else if (!strcmp(tag, "hohen-d"))
      engine->hohen.d = readReal(pFile);
    else if (!strcmp(tag, "hohen-c2"))
      engine->hohen.C2 = readReal(pFile);
    else if (!strcmp(tag, "hohen-c3"))
      engine->hohen.C2 = readReal(pFile);
    else if (!strcmp(tag, "hohen-alpha"))
      engine->hohen.alpha = readReal(pFile);
    else if (!strcmp(tag, "hohen-ref-temp"))
      engine->hohen.Tref = readReal(pFile);
    else if (!strcmp(tag, "hohen-ref-pressure"))
      engine->hohen.Pref = readReal(pFile);
    else if (!strcmp(tag, "hohen-pressure-scale"))
      engine->hohen.Pscale = readReal(pFile);
    
    /* ----- flamelet parameters ----- */
    else if (!strcmp(tag, "FLAMELET")) {
      readFlame = 1;
      driver->FLNO = readInteger(pFile);
      if (driver->FLNO > driver->maxFlamelets) {
        printf("\n!!!!!!!!!!     WARNING      !!!!!!!!!!\n");
        printf("Flamelet number (%i) > Max number of Flamelets (%i)\n",
               driver->FLNO, driver->maxFlamelets);
        printf("Flamelet not read\n");
        readFlame = 0;
      }
      else 
        driver->FLNO--; // adjust flamelet number to array numbering
      
      while (readFlame) {
        fscanf(pFile, "%s", tag);
        if (!strcmp(tag, "ENDFLAMELET")) {
          fscanf(pFile, "\n");
          readFlame = 0;
        }
        else if (tag[0] == '!' || tag[0] == '#')
          fgets(comment, 200, pFile);
        else if (!strcmp(tag, "dimensions"))
          flame->nDim = readInteger(pFile);
        else if (!strcmp(tag, "gridPoints")||!strcmp(tag, "gridPoints1"))
          flame->n1 = readInteger(pFile);
        else if (!strcmp(tag, "gridPoints2"))
          flame->n2 = readInteger(pFile);
        else if (!strcmp(tag, "gridPivot1"))
          io->nPivot1 = readRealArraySize(pFile, &io->gPivot1);
        else if (!strcmp(tag, "gridPercent1"))
          readRealArraySize(pFile, &io->gPerc1);
        else if (!strcmp(tag, "gridPivot2"))
          io->nPivot2 = readRealArraySize(pFile, &io->gPivot2);
        else if (!strcmp(tag, "gridPercent2"))
          readRealArraySize(pFile, &io->gPerc2);
        else if (!strcmp(tag, "OxidizerSide"))
          io->bcType[0]=readBC(pFile, &io->nOx, &io->oxName, &io->oxT,
                               &io->oxY, &io->oxH);
        else if (!strcmp(tag, "FuelSide") || !strcmp(tag, "FuelSide1"))
          io->bcType[1]=readFuelBC(pFile,&io->nFuel1,&io->fuelName1,
                                   &io->fuelT1,&io->fuelY1,&io->fuelH1,
                                   &io->hVap1,&io->cpLiq1,&io->fuelTref1,
                                   io->satData);
        else if (!strcmp(tag, "FuelSide2"))
          io->bcType[2]=readFuelBC(pFile,&io->nFuel2,&io->fuelName2,
                                   &io->fuelT2,&io->fuelY2,&io->fuelH2,
                                   &io->hVap2,&io->cpLiq2,&io->fuelTref2,
                                   io->satData);
        else if (!strcmp(tag, "output-species"))
          io->nOutputs = readNames(pFile, &io->outputSpecies);
        else if (!strcmp(tag, "monitorPoints"))
          io->nMonitorZ = readInteger(pFile);
        else if (!strcmp(tag, "monitor-species"))
          io->nMonitorY = readNames(pFile, &io->monitorSpecies);
        else {
          fprintf(stderr, "\nWarning: ");
          fprintf(stderr, "%s not a valid flamelet variable\n", tag);
        }
      } // end while readFlame
    } // end read flamelet
    else {
/*       fprintf(stderr, "\nWarning: "); */
/*       fprintf(stderr, "%s not associated with an input variable\n", tag); */
    }
  } // end loop over file
  
  fclose(pFile);
}
/* -------------------------------------------------------------------------- */

void readFlamelet(char *flameFile, ioStruct *io, driveStruct *driver, 
                  timeStruct *time, flameStruct *flame, flowStruct *flow)
{
  int i, j, iCurSpec;

  /* strings to store the tag for each line and comments (if any) */
  char  tag[50], comment[200];
  
  FILE   *pFile;
  pFile = fopen(flameFile, "r"); /* open input file */
  
  if (!pFile) {
    fprintf(stderr,"Error: could not open flamelet file \"%s\"\n",flameFile);
    exit(2);
  }

  while (!feof(pFile)) {
    /* read in string */
    fscanf(pFile, "%s", tag);

    /* first check to see if end of input is desired */
    if (!strcmp(tag, "ENDINPUT"))
      break;
    /* check to see if it is a comment */
    else if (tag[0] == '!' || tag[0] == '#')
      fgets(comment, 200, pFile);
    
    /* ----- read in header data ----- */
    else if (!strcmp(tag, "title"))
      readString(pFile, io->title);
    else if (!strcmp(tag, "author"))
      readString(pFile, io->author);
    else if (!strcmp(tag, "date"))
      readString(pFile, io->date);
    
    /* ----- read in any scalar parameters ----- */
    else if (!strcmp(tag, "time")) 
      time->t = time->t_m1 = readReal(pFile);
    else if (!strcmp(tag, "pressure"))
      flow->P = flow->P_m1 = readReal(pFile);
    else if (!strcmp(tag, "DPDT"))
      flow->Dt_P = readReal(pFile);
    else if (!strcmp(tag, "Z1_min") || !strcmp(tag, "H1_min"))
      flow->phi1min_m1 = readReal(pFile);
    else if (!strcmp(tag, "Z1_max") || !strcmp(tag, "H1_max"))
      flow->phi1max_m1 = readReal(pFile);
    else if (!strcmp(tag, "Z2_min") || !strcmp(tag, "H2_min"))
      flow->phi2min_m1 = readReal(pFile);
    else if (!strcmp(tag, "Z2_max") || !strcmp(tag, "H2_max"))
      flow->phi2max_m1 = readReal(pFile);
    else if (!strcmp(tag, "Z1_mean") || !strcmp(tag, "H1_mean"))
      io->mean[1] = readReal(pFile);
    else if (!strcmp(tag, "Z2_mean") || !strcmp(tag, "H2_mean"))
      io->mean[2] = readReal(pFile);
    else if (!strcmp(tag, "Z1_variance") || !strcmp(tag, "H1_variance"))
      io->var[1] = readReal(pFile);
    else if (!strcmp(tag, "Z2_variance") || !strcmp(tag, "H2_variance"))
      io->var[2] = readReal(pFile);
    else if (!strcmp(tag, "Z1_Z2_variance") || !strcmp(tag, "Z1_H2_variance")
             || !strcmp(tag, "H1_Z2_variance") )
      io->var[3] = readReal(pFile);
    else if (!strcmp(tag, "FlameLoc"))
      flame->flameZ = readReal(pFile);
    else if (!strcmp(tag, "Tmax"))
      flame->Tmax = readReal(pFile);
    else if (!strcmp(tag, "transformation"))
      readString(pFile, io->transZ2);
    else if (!strcmp(tag, "burningVelocity"))
      flame->S_L = readReal(pFile);
    else if (!strcmp(tag, "FlameThickness"))
      flame->thick = readReal(pFile);
    else if (!strcmp(tag, "fuel-air-equivalence-ratio"))
      flame->phi = readReal(pFile);
    
    /* ----- read in the fuel and oxidizer boundary conditions here ---- */
    else if (!strcmp(tag, "OxidizerSide") || !strcmp(tag, "unburnt")) {
      io->bcType[0] = readBC(pFile, &io->nOx, &io->oxName, &io->oxT, &io->oxY,
                             &io->oxH);
      for (i=0; i<io->nOx; i++) stripHeader(io->oxName[i]);
    }
    else if (!strcmp(tag, "FuelSide") || !strcmp(tag, "FuelSide1")) {
      io->bcType[1] = readFuelBC(pFile,&io->nFuel1,&io->fuelName1,&io->fuelT1,
                                 &io->fuelY1,&io->fuelH1,&io->hVap1,&io->cpLiq1,
                                 &io->fuelTref1,io->satData);
      for (i=0; i<io->nFuel1; i++) stripHeader(io->fuelName1[i]);
    }
    else if (!strcmp(tag, "FuelSide2")) {
      io->bcType[2] = readFuelBC(pFile,&io->nFuel2,&io->fuelName2,&io->fuelT2,
                                 &io->fuelY2,&io->fuelH2,&io->hVap2,&io->cpLiq2,
                                 &io->fuelTref2,io->satData);
      for (i=0; i<io->nFuel2; i++) stripHeader(io->fuelName2[i]);
    }
    else if (!strcmp(tag, "numOfSpecies"))
      io->nInputs = readInteger(pFile);
    else if (!strcmp(tag, "gridPoints") || !strcmp(tag, "gridPoints1"))
      flame->n1 = readInteger(pFile);
    else if (!strcmp(tag, "gridPoints2")) {
      flame->n2 = readInteger(pFile);
      flame->nDim = 2;
    }
      
    /* ----- read in the body of data here ----- */
    else if (!strcmp(tag, "body")) {
      /* perform some memory allocation */
      io->initY = (double**) calloc(io->nInputs, sizeof(double*));
      io->inputSpecies = (char**) calloc(io->nInputs, sizeof(char*));
      for (i=0; i < io->nInputs; i++)
        io->inputSpecies[i] = (char*) calloc(25, sizeof(char));
      /* initialise storage for miscellaneous variables */
      io->miscIn    = (double**) calloc(1, sizeof(double*));
      io->miscNames = (char**) calloc(1, sizeof(char*));
      io->nMisc = 0;
      iCurSpec  = 0;
      io->initSoot = (double**) calloc(7, sizeof(double*));
      io->initSootMom = 1;
      flame->nSootMom = 0;
      
      /* read in all the desired arrays until trailer is reached */
      fscanf(pFile, "%s", tag);
      do {
        if (!strcmp(tag, "Z1") || !strcmp(tag, "y"))
          readRealArray(pFile, flame->n1, &flame->x1);
        else if (!strcmp(tag, "Z1+Z2")) {
          readRealArray(pFile, flame->n1, &flame->x1);
          driver->TRANS_Z2 = 2;
        }
        else if (!strcmp(tag, "H1")) {
          readRealArray(pFile, flame->n1, &flame->x1);
          driver->H_DIM = 1;
        }
        else if (!strcmp(tag, "Z2"))
          readRealArray(pFile, flame->n2, &flame->x2);
        else if (!strcmp(tag, "Z2/(1-Z1)")) {
          readRealArray(pFile, flame->n2, &flame->x2);
          driver->TRANS_Z2 = 1;
        }
        else if (!strcmp(tag, "Z2/(Z1+Z2)")) {
          readRealArray(pFile, flame->n2, &flame->x2);
          driver->TRANS_Z2 = 2;
        }
        else if (!strcmp(tag, "H2")) {
          readRealArray(pFile, flame->n2, &flame->x2);
          driver->H_DIM = 2;
        }
        else if (!strncmp(tag, "massfraction", 12)) {
          stripHeader(tag);
          strcpy(io->inputSpecies[iCurSpec], tag);
          readRealArray(pFile, flame->n1*flame->n2, &io->initY[iCurSpec]);
          iCurSpec++;
        }
/*         else if (!strncmp(tag, "ProdRate", 8)) { */
/*           stripHeader(tag); */
/*           strcpy(io->inputSpecies[iCurSpec], tag); */
/*           readRealArray(pFile, flame->n1*flame->n2, &io->initY[iCurSpec]); */
/*           iCurSpec++; */
/*         } */
        else if (!strcmp(tag, "temperature"))
          readRealArray(pFile, flame->n1*flame->n2, &io->initT);
        else if (!strcmp(tag, "chi") || !strcmp(tag, "chi1"))
          readRealArray(pFile, flame->n1*flame->n2, &flow->chi1);
        else if (!strcmp(tag, "chi2")) {
          if (flame->n2 == 1) {
            fprintf(stderr, "ERROR: chi2 requires n2\n");
            exit(2);
          }
          readRealArray(pFile, flame->n1*flame->n2, &flow->chi2);
        }
        else if (!strcmp(tag, "chi12"))
          readRealArray(pFile, flame->n1*flame->n2, &flow->chi12);
        else if (!strcmp(tag, "pdf"))
          readRealArray(pFile, flame->n1*flame->n2, &flow->pdf);
        else if (!strcmp(tag, "chiA"))
          readRealArray(pFile, flame->n1*flame->n2, &flame->chiA);
        else if (!strcmp(tag, "chiB"))
          readRealArray(pFile, flame->n1*flame->n2, &flame->chiB);
        else if (!strcmp(tag, "chiC"))
          readRealArray(pFile, flame->n1*flame->n2, &flame->chiC);
        else if (!strcmp(tag, "TotalEnthalpy"))
          readRealArray(pFile, flame->n1*flame->n2, &flame->H);
        else if (!strcmp(tag, "Density") || !strcmp(tag, "density"))
          readRealArray(pFile, flame->n1*flame->n2, &flame->rho);
        else if (!strcmp(tag, "soot-number-density_M00")) {
          readRealArray(pFile, flame->n1*flame->n2, &io->initSoot[0]);
          io->initSootMom = 0;
          flame->nSootMom++;
        }
        else if (!strcmp(tag, "soot-volume_M10")) {
            readRealArray(pFile, flame->n1*flame->n2, &io->initSoot[1]);
            flame->nSootMom++;
        }
        else if (!strcmp(tag, "soot-surface_M01")) {
          readRealArray(pFile, flame->n1*flame->n2, &io->initSoot[2]);
          flame->nSootMom++;
        }
        else if (!strcmp(tag, "soot-volume-var_M20")) {
          readRealArray(pFile, flame->n1*flame->n2, &io->initSoot[3]);
          flame->nSootMom++;
        }
        else if (!strcmp(tag, "soot-vol-surf-var_M11")) {
          readRealArray(pFile, flame->n1*flame->n2, &io->initSoot[4]);
          flame->nSootMom++;
        }
        else if (!strcmp(tag, "soot-surface-var_M02")) {
          readRealArray(pFile, flame->n1*flame->n2, &io->initSoot[5]);
          flame->nSootMom++;
        }
        else if (!strcmp(tag, "soot-delta-weight_N0")) {
          readRealArray(pFile, flame->n1*flame->n2, &io->initSoot[6]);
          flame->nSootMom++;
        }
        //else if (!strcmp(tag, "visc-src")) {
	//  /* nothing to be done */
        //}
        //else if (!strcmp(tag, "global-pdf")) {
        //  readRealArray(pFile, flame->n1*flame->n2, &io->gPDF);
        //}
        else {
          io->miscIn = (double**) realloc(io->miscIn,
                                          (io->nMisc+1)*sizeof(double*));
          io->miscIn[io->nMisc] = (double*) calloc(flame->n1*flame->n2,
                                                   sizeof(double));
          io->miscNames = (char**) realloc(io->miscNames,
                                           (io->nMisc+1)*sizeof(char*));
          io->miscNames[io->nMisc] = (char*) calloc(25, sizeof(char));
          strcpy(io->miscNames[io->nMisc], tag);
          readRealArray(pFile,flame->n1*flame->n2,&io->miscIn[io->nMisc]);
          io->nMisc++;
        }
        
        /* scan in next line to check whether body is fully read */
        fscanf(pFile, "%s", tag);
        if (flame->nDim == 2 && !strcmp(tag, "chi1-min")) {
          fclose(pFile);
          return;
        }
      } while (strcmp(tag, "trailer"));
      
    } /* end read body */
  } /* end loop over file */
  
  fclose(pFile);
}
/* -------------------------------------------------------------------------- */

void readChi(char *flameFile, ioStruct *io, int i)
{

  /* strings to store the tag for each line and comments (if any) */
  char  tag[50], comment[200];
  int nGrid1=1, nGrid2=1;
  
  FILE   *pFile;
  pFile = fopen(flameFile, "r"); /* open input file */
  
  if (!pFile) {
    fprintf(stderr,"Error: could not open flamelet file \"%s\"\n",flameFile);
    exit(2);
  }

  while (!feof(pFile)) {
    /* read in string */
    fscanf(pFile, "%s", tag);

    /* first check to see if end of input is desired */
    if (!strcmp(tag, "ENDINPUT"))
      break;
    /* check to see if it is a comment */
    else if (tag[0] == '!' || tag[0] == '#')
      fgets(comment, 200, pFile);
    
    /* ----- read in any scalar parameters ----- */
    else if (!strcmp(tag, "time"))
      io->chiData[i].tChi = readReal(pFile);
    
    else if (!strcmp(tag, "Pressure"))
      io->chiData[i].P = readReal(pFile);
    
    /* ----- read in the fuel and oxidizer boundary conditions here ---- */
    else if (!strcmp(tag, "gridPoints") || !strcmp(tag, "gridPoints1"))
      nGrid1 = readInteger(pFile);
    else if (!strcmp(tag, "gridPoints2")) {
      nGrid2 = readInteger(pFile);
    }

    /* ----- read in the body of data here ----- */
    else if (!strcmp(tag, "body")) {
      /* initialise storage for miscellaneous variables */
      io->miscIn    = (double**) calloc(1, sizeof(double*));
      io->miscNames = (char**) calloc(1, sizeof(char*));
      io->nMisc = 0;

      /* read in all the desired arrays until trailer is reached */
      fscanf(pFile, "%s", tag);
      do {
        if (!strcmp(tag, "chiA")) {
          readRealArray(pFile, nGrid1*nGrid2, &io->chiData[i].chiRaw);
	  io->chiData[i].n1 = nGrid1;
	  io->chiData[i].n2 = nGrid2;
        }
        if (!strcmp(tag, "Z1")) 
          readRealArray(pFile, nGrid1*nGrid2, &io->chiData[i].x1);
        else {
        if (!strcmp(tag, "global-pdf")) {
          readRealArray(pFile, nGrid1*nGrid2, &io->chiData[i].pdfRaw);
        }
        if (!strcmp(tag, "h-src")) {
          readRealArray(pFile, nGrid1*nGrid2, &io->chiData[i].hSrcRaw);
        }
 	// TODO: visc-src, P
          /* skip */ 
          //io->miscIn = (double**) realloc(io->miscIn,
          //                                (io->nMisc+1)*sizeof(double*));
          //io->miscIn[io->nMisc] = (double*) calloc(nGrid1*nGrid2,
          //                                         sizeof(double));
          //io->miscNames = (char**) realloc(io->miscNames,
          //                                 (io->nMisc+1)*sizeof(char*));
          //io->miscNames[io->nMisc] = (char*) calloc(25, sizeof(char));
          //strcpy(io->miscNames[io->nMisc], tag);
          //readRealArray(pFile, nGrid1*nGrid2, &io->miscIn[io->nMisc]);
          //io->nMisc++;
        }
        
        /* scan in next line to check whether body is fully read */
        fscanf(pFile, "%s", tag);
        /*if (flame->nDim == 2 && !strcmp(tag, "chi1-min")) {
          fclose(pFile);
          return;
        }*/
      } while (strcmp(tag, "trailer"));
      
    } /* end read body */
  } /* end loop over file */
  fclose(pFile);
}
/* -------------------------------------------------------------------------- */

double readFlowData(char *fileName, driveStruct *driver, flowStruct *flow)
{
  int n1, n2, nDim;
  // strings to store the tag for each line and comments (if any)
  char  tag[50], comment[200];
  double t_cur = -1.0;

  FILE *pFile;
  pFile = fopen(fileName, "r");

  if (!pFile)
    {
      fprintf(stderr, "Error: flow data file \"%s\" not found\n", fileName);
      exit(2);
    }

  while (!feof(pFile))
    { 
      /* read in string */
      fscanf(pFile, "%s", tag);
      
      /* first check to see if end of input is desired */
      if (!strcmp(tag, "ENDINPUT"))
        break;
      /* check to see if it is a comment */
      else if (tag[0] == '!' || tag[0] == '#')
        fgets(comment, 200, pFile);
      else if (!strcmp(tag, "time"))
        t_cur = readReal(pFile);
      else if (!strcmp(tag, "Z1_min") || !strcmp(tag, "H1_min"))
        flow->phi1min = readReal(pFile);
      else if (!strcmp(tag, "Z1_max") || !strcmp(tag, "H1_max"))
        flow->phi1max = readReal(pFile);
      else if (!strcmp(tag, "Z2_min") || !strcmp(tag, "H2_min"))
        flow->phi2min = readReal(pFile);
      else if (!strcmp(tag, "Z2_max") || !strcmp(tag, "H2_max"))
        flow->phi2max = readReal(pFile);
      else if (!strcmp(tag, "H1_min_corr"))
        flow->phi1min_corr = readReal(pFile);
      else if (!strcmp(tag, "H2_min_corr"))
        flow->phi2min_corr = readReal(pFile);
      else if (!strcmp(tag, "H1_max_corr"))
        flow->phi1max_corr = readReal(pFile);
      else if (!strcmp(tag, "H2_max_corr"))
        flow->phi2max_corr = readReal(pFile);
      else if (!strcmp(tag, "pressure") && driver->PRES_MODE == -1)
        flow->P = readReal(pFile);
      else if (!strcmp(tag, "strain-rate"))
        readOptArray(pFile, 3, flow->strain);
      else if (!strcmp(tag, "strain-rate-0"))
        flow->strain[0] = readReal(pFile);
      else if (!strcmp(tag, "strain-rate-1"))
        flow->strain[1] = readReal(pFile);
      else if (!strcmp(tag, "strain-rate-2"))
        flow->strain[2] = readReal(pFile);
      else if (!strcmp(tag, "ConstantLewisNumbers"))
        {
          readString(pFile, tag);
          if (!strcmp(tag, "True"))
            flow->constLewisNum = 1;
          else
            flow->constLewisNum = 0;
        }
      else if (!strcmp(tag, "gridPoints") || !strcmp(tag, "gridPoints1"))
        {
          n1 = readInteger(pFile);
          n2 = 1;
        }
      else if (!strcmp(tag, "gridPoints2"))
        {
          n2 = readInteger(pFile);
          nDim = 2;
        }
      
      else if (!strcmp(tag, "body"))
        {
          if (driver->readBody == 0)
            {
              fclose(pFile);
              return t_cur;
            }
          
          // read in all the desired arrays until trailer is reached
          do {
            if (!strcmp(tag, "chi1"))
              if (driver->CHI_MODE == -1)
                readRealArray(pFile, n1*n2, &flow->chi1);
              else
                readRealArray(pFile, n1*n2, &flow->chi1_m1);
            else if (!strcmp(tag, "chi2"))
              {
                if (n2 == 1)
                  {
                    fprintf(stderr, "Error: chi2 requires n2\n");
                    exit(2);
                  }
                if (driver->CHI_MODE == -1)
                  readRealArray(pFile, n1*n2, &flow->chi2);
                else
                  readRealArray(pFile, n1*n2, &flow->chi2_m1);
              }
            else if (!strcmp(tag, "chi12"))
              {
                if (n2 == 1)
                  {
                    fprintf(stderr, "Error: chi12 requires n2\n");
                    exit(2);
                  }
                if (driver->CHI_MODE == -1)
                  readRealArray(pFile, n1*n2, &flow->chi12);
                else
                  readRealArray(pFile, n1*n2, &flow->chi12_m1);
              }
            else if (!strcmp(tag, "GlobalPDF") || !strcmp(tag, "pdf12") ||
                     !strcmp(tag, "pdf") || !strcmp(tag, "pdf1") )
              readRealArray(pFile, n1*n2, &flow->pdf);
            else if (!strcmp(tag, "Sdot1"))
              readRealArray(pFile, n1, &flow->avgSdot_m1);
            else if (!strcmp(tag, "Sdot2"))
              readRealArray(pFile, n2, &flow->avgSdot_m1);
            else if (!strcmp(tag, "SdotM"))
              readRealArray(pFile, n2, &flow->avgSdot_M);
            else if (!strcmp(tag, "cond-strain-rate"))
              readRealArray(pFile, n2, &flow->strainCondZ);
            else if (!strcmp(tag, "enthalpy_c_z1"))
              readRealArray(pFile, n1, &flow->avgH);
            else if (!strcmp(tag, "enthalpy_c_z2"))
              readRealArray(pFile, n2, &flow->avgH);
            // scan in next line to check whether body is fully read
            fscanf(pFile, "%s", tag);
          } while (strcmp(tag, "trailer"));
        } // end read body
    } // end loop over file

  fclose(pFile);
  return t_cur;
}
/* -------------------------------------------------------------------------- */

double readFileTime(char *fileName)
{
  double t = -1.0;
  char   tag[50];
  FILE *pFile = fopen(fileName, "r");

  while (!feof(pFile))
    { 
      /* read in string */
      fscanf(pFile, "%s", tag);
      
      if (!strcmp(tag, "time"))
        {
          t = readReal(pFile);
          break;
        }
    }

  if (t < 0.0)
    fprintf(stderr, "No time found in file: %s\n", fileName);

  fclose(pFile);
  return t;
}
/* -------------------------------------------------------------------------- */

double readFileScalar(char *fileName, char *varName)
{
  int    found = 0;
  double var   = -1.0;
  char   tag[50];
  FILE   *pFile = fopen(fileName, "r");

  while (!feof(pFile))
    { 
      /* read in string */
      fscanf(pFile, "%s", tag);
      
      if (!strcmp(tag, varName))
        {
          var = readReal(pFile);
          found = 1;
          break;
        }
    }

  if (!found)
    fprintf(stderr, "%s not found in file: %s\n", varName, fileName);

  fclose(pFile);
  return var;
}
/* -------------------------------------------------------------------------- */

void readVariableArray(char *fileName, char *varName, double **var)
{
  int n1 = 1, n2 = 1, found = 0;
  // strings to store the tag for each line and comments (if any)
  char  tag[50], comment[200];

  FILE *pFile = fopen(fileName, "r");

  if (!pFile)
    {
      fprintf(stderr, "Error: Can't open file %s\n", fileName);
      exit(2);
    }
  while (!feof(pFile))
    { 
      /* read in string */
      fscanf(pFile, "%s", tag);
      
      /* first check to see if end of input is desired */
      if (!strcmp(tag, "ENDINPUT"))
        break;
      /* check to see if it is a comment */
      else if (tag[0] == '!' || tag[0] == '#')
        fgets(comment, 200, pFile);
      else if (!strcmp(tag, "gridPoints") || !strcmp(tag, "gridPoints1"))
        n1 = readInteger(pFile);
      else if (!strcmp(tag, "gridPoints2"))
        n2 = readInteger(pFile);
      else if (!strcmp(tag, varName))
        {
          found = 1;
          readRealArray(pFile, n1*n2, var);
        }

    } // end loop over file

  if (!found)
    fprintf(stderr, "%s not found in file: %s\n", varName, fileName);

  fclose(pFile);
}
/* -------------------------------------------------------------------------- */
void writeRates(char *fileName, chemStruct *chem)
{
  int i;
  FILE   *pFile;                 /* output stream */
  pFile = fopen(fileName, "w"); /* open output file */

  if (!pFile) {
    fprintf(stderr, "Error: Could not open file %s to write flamelet\n",
            fileName);
    exit(2);
  }

  /* ----- write soot model parameters ----- */
  if(chem->nSootPAH > 0){
    
    fprintf(pFile, " ! ------------------------- soot model parameters ---------------------------- ! \n");
    fprintf(pFile, " soot-model %6s\n", chem->sootModel);
    fprintf(pFile, " soot-mode %6i\n", chem->sootMode);
    fprintf(pFile, " soot-moment-order %6i\n", chem->sootOrder);
    
    fprintf(pFile, " soot-active-source* ");
    for (i=0; i<chem->nSootSrc; i++)
      if (chem->sootSrcOn[i])
        fprintf(pFile, "%s ",chem->sootSrcNames[i]);
    fprintf(pFile, "\n");
    
    fprintf(pFile, " soot-active-mom* ");
    for (i=0; i < chem->nSootMomOn; i++)
      if (chem->sootMomOn[i])
        fprintf(pFile, "%s ",chem->sootMomNames[i]);
    fprintf(pFile, "\n");

    fprintf(pFile, " soot-PAH-species* ");
    for (i=0; i<chem->nSootPAH; i++ ) 
      fprintf(pFile,"%s ", chem->sootPAHSpec[i]);
     fprintf(pFile,"\n");
  }
 
  /* ----- write body ----- */
  fprintf(pFile, "body \n");
  
  for (i=0;i<chem->nRxnTotal;i++)
      fprintf(pFile,"%6i%25.15e%25.15e\n",i,chem->w_int[i],chem->w_int_pre[i]);

  fclose(pFile);
}

void readRates(char *fileName, chemStruct *chem)
{
  int i;
  double tmp;
  FILE   *pFile;                 /* output stream */
  char tag[200],tag1[40],tag2[40],tag3[40];
  pFile = fopen(fileName, "r"); /* open output file */

  if (!pFile) {
    fprintf(stderr, "Error: Could not open file %s to read rates\n",
            fileName);
    exit(2);
  }

  while (!feof(pFile))
    {
      fscanf(pFile, "%[^\n]\n", tag);  // read line
      sscanf(tag, "%s %s %s", tag1,tag2,tag3);  // parse values
      if (atoi(tag1) < chem->nRxnTotal) {
         chem->w_int[atoi(tag1)] = atof(tag2);
         chem->w_int_pre[atoi(tag1)] = atof(tag3);
      }
    }

  fclose(pFile);
}

void writeFlamelet(char *fileName, ioStruct *io, driveStruct *driver,
                   timeStruct *time, flameStruct *flame, flowStruct *flow)
{
  int i, j, k, r, row, col, perLine = 5;
  char dim1, dim2;

  /* assign the appropriate co-ordinate names */
  if (driver->H_DIM == 0) {
    dim1 = 'Z';
    dim2 = 'Z';
  }
  else if (driver->H_DIM == 1) {
    dim1 = 'H';
    dim2 = 'Z';
  }
  else if (driver->H_DIM == 2) {
    dim1 = 'Z';
    dim2 = 'H';
  }
  
  FILE   *pFile;                 /* output stream */
  pFile = fopen(fileName, "w"); /* open output file */

  if (!pFile) {
    fprintf(stderr, "Error: Could not open file %s to write flamelet\n",
            fileName);
    exit(2);
  }
  
  /* ----- write header ----- */
  fprintf(pFile, "header\n");
  fprintf(pFile, "title = \"%s\"\n", io->title);
  fprintf(pFile, "mechanism = \"%s\"\n", io->mechanism);
  fprintf(pFile, "author = \"%s\"\n", io->author);
  fprintf(pFile, "date = \"%s\"\n\n", io->date);
  
  fprintf(pFile, "fuel = \"");
  for (i=0; i<io->nFuel1; i++) {
    fprintf(pFile, "%s", io->fuelName1[i]);
    if (i+1 < io->nFuel1)
      fprintf(pFile, "/");
  }
  fprintf(pFile, "\"\n");
  fprintf(pFile, "time = %15.6e [s]\n", time->t);
  fprintf(pFile, "pressure = %15.6e [Pa]\n", flow->P);
  fprintf(pFile, "DPDT     = %15.6e [Pa/s]\n", flow->Dt_P);
  fprintf(pFile, "%c1_min = %15.6e\n", dim1, flow->phi1min);
  fprintf(pFile, "%c1_max = %15.6e\n", dim1, flow->phi1max);
  if (flame->nDim == 2) {
    fprintf(pFile, "%c2_min = %15.6e\n", dim2, flow->phi2min);
    fprintf(pFile, "%c2_max = %15.6e\n", dim2, flow->phi2max); 
  }
  fprintf(pFile, "\n");
  /*fprintf(pFile, "chi = %5.2f [1/s]\n", flow->chi1Ref);
  if (flow->constLewisNum)
    fprintf(pFile, "ConstantLewisNumbers = \"True\"\n");
  else
    fprintf(pFile, "ConstantLewisNumbers = \"False\"\n");
  fprintf(pFile, "FlameLoc = %.8f\n", flame->flameZ);
  fprintf(pFile, "Tmax = %5.2f [K]\n", flame->Tmax);*/
  
  /* print boundary conditions */
  fprintf(pFile, "FuelSide\n");
  fprintf(pFile, "begin\n");
  if (io->bcType[1] == 'T')
    fprintf(pFile, " Temperature = %8.3f [K]\n", io->fuelT1);
  else if (io->bcType[1] == 'H')
    fprintf(pFile, " Enthalpy = %8.4e [J/kg]\n", io->fuelH1);
  else if (io->bcType[1] == 'L') {
    fprintf(pFile, " Temperature  = %8.3f [K]\n", io->fuelT1);
    fprintf(pFile, " Latent-heat = %8.4e [J/kg]\n", io->hVap1);
    fprintf(pFile, " Liquid-cp    = %8.2f [J/(kg-K)]\n", io->cpLiq1);
    fprintf(pFile, " Temp-ref     = %8.3f [K]\n", io->fuelTref1);
  }
  for (i=0; i<io->nFuel1; i++)
    fprintf(pFile," Massfraction-%s = %7.5f\n",io->fuelName1[i],io->fuelY1[i]);
  fprintf(pFile, "end\n\n");
  
  if (flame->nDim == 2) {
    fprintf(pFile, "FuelSide2\n");
    fprintf(pFile, "begin\n");
    if (io->bcType[2] == 'T' || io->bcType[2] == 'E')
      fprintf(pFile, " Temperature = %8.3f [K]\n", io->fuelT2);
    else if (io->bcType[2] == 'H')
      fprintf(pFile, " Enthalpy = %8.4e [J/kg]\n", io->fuelH2);
    else if (io->bcType[2] == 'L') {
      fprintf(pFile, " Temperature  = %8.3f [K]\n", io->fuelT2);
      fprintf(pFile, " Latent-heat = %8.4e [J/kg]\n", io->hVap2);
      fprintf(pFile, " Liquid-cp    = %8.2f [J/(kg-K)]\n", io->cpLiq2);
      fprintf(pFile, " Temp-ref     = %8.3f [K]\n", io->fuelTref2);
    }
    for (i=0; i<io->nFuel2; i++)
      fprintf(pFile, " Massfraction-%s = %7.5f\n", io->fuelName2[i],
              io->fuelY2[i]);
    fprintf(pFile, "end\n\n");
  }
  
  fprintf(pFile, "OxidizerSide\n");
  fprintf(pFile, "begin\n");
  if (io->bcType[1] == 'T')
    fprintf(pFile, " Temperature = %8.3f [K]\n", io->oxT);
  else if (io->bcType[1] == 'H')
    fprintf(pFile, " Enthalpy = %8.4e [J/kg]\n", io->oxH);
  for (i=0; i<io->nOx; i++)
    fprintf(pFile, " Massfraction-%s = %7.5f\n", io->oxName[i], io->oxY[i]);
  fprintf(pFile, "end\n\n");
    
  fprintf(pFile, "numOfSpecies = %i\n", io->nOutputs);
  fprintf(pFile, "gridPoints = %i\n", flame->n1);
  if (flame->nDim == 2)
    fprintf(pFile, "gridPoints2 = %i\n", flame->n2);
  fprintf(pFile, "\n");
  
  /* ----- write body ----- */
  fprintf(pFile, "body\n");
  /* print the mesh */
  fprintf(pFile, "%c1\n", dim1);
  writeArray(pFile, 1, flame->n1, 1, 0, perLine, flame->x1);
  if (flame->nDim == 2) {
    fprintf(pFile, "%c2\n", dim2);
    writeArray(pFile, 1, flame->n2, 1, 0, perLine, flame->x2);
  }
  
  /* print the temperature profile */
  fprintf(pFile, "temperature [K]\n");
  if (driver->ENERGY)
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->T);
  else
    writeArray(pFile, flame->n2, flame->n1, flame->nVars, flame->sT, perLine,
               flame->Y);
  
  /* print the species that are tracked (from input file) */
  for (i=0; i<io->nOutputs; ++i) {
    fprintf(pFile, "massfraction-%s\n", io->outputSpecies[i]);
    writeArray(pFile, flame->n2, flame->n1, flame->nVars, io->outputIndex[i], 
               perLine, flame->Y);
  } /* end print species */
  
  /* print the scalar dissipation rate */
  if (flame->chiA != NULL && io->outputChi) {
    if (driver->CHI_MODE == 2)
      extractVarArray(flame->nVars,flame->n1,1,flame->sChi,flame->Y,
                      flame->chiA);
    fprintf(pFile, "chiA [1/s]\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->chiA);
  }
  if (flame->nDim == 2) {
    if (flame->chiB != NULL && io->outputChi) {
      fprintf(pFile, "chiB [1/s]\n");
      writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->chiB);
    }
    if (flame->chiC != NULL && io->outputChi) {
      fprintf(pFile, "chiC [1/s]\n");
      writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->chiC);
    }
  }
  if (io->outputViscSrc) {
    /* print the viscous heat production */
    fprintf(pFile, "h-src [J/(kg*s)]\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->visc_src);
  }
  if (io->outputConvVel) {
    /* print the viscous heat production */
    fprintf(pFile, "vT\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->vT);
    if (flame->vMxmyp != NULL) {
       fprintf(pFile, "vMxmyp\n");
       writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->vMxmyp);
    }
    if (flame->vMxy != NULL) {
       fprintf(pFile, "vMxy\n");
       writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->vMxy);
    }
  }
  if (io->outputGPDF) {
    /* print the viscous heat production */
    fprintf(pFile, "global-pdf [kg]\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->gPDF);
  }
  if (flame->sgCoeff != NULL) {
    fprintf(pFile, "HMOM-sgCoeff\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->sgCoeff);
  }
  if (flame->oxCoeff != NULL) {
    fprintf(pFile, "HMOM-oxCoeff\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->oxCoeff);
  }
  if (flame->o2Coeff != NULL) {
    fprintf(pFile, "HMOM-o2Coeff\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->o2Coeff);
  }
  if (flame->dimProdRate != NULL) {
    fprintf(pFile, "HMOM-dimProdRate\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->dimProdRate);
  }
  if (flame->sqrtT != NULL) {
    fprintf(pFile, "HMOM-sqrtT\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->sqrtT);
  }
  if (flame->T_mu != NULL) {
    fprintf(pFile, "HMOM-Tmu\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->T_mu);
  }
  if (flame->muSqrtMw_rhoSqrtT != NULL) {
    fprintf(pFile, "HMOM-muSqrtMw_rhoSqrtT\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->muSqrtMw_rhoSqrtT);
  }
  if (io->outputHtot) { 
    /* print the total enthalpy */
    fprintf(pFile, "TotalEnthalpy [J/kg]\n");
    if (driver->ENERGY)
      writeArray(pFile,flame->n2,flame->n1,flame->nVars,flame->sT,perLine,
                 flame->Y);
    else
      writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->H);
  }

  if (io->outputHR) {
    /* print the src term of enthalpy */
    fprintf(pFile, "HeatRelease [J/(kg-s)]\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->HRrate);
  }

  if (driver->output_tsrc) {
    /* print the src term of convection */
    fprintf(pFile, "Convection-Src\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->conv_src);
      
    /* print the src term of chemistry */
    fprintf(pFile, "Chemistry-Src\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->chem_src);
    
    /* print the src term of convection */
    fprintf(pFile, "Diffusion-Src\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->diff_src);
  }

  if (io->outputConv) {
    fprintf(pFile, "convection-coefficient\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->convCoeff);
  }

  if (io->outputRho) { /* print the density */
    fprintf(pFile, "Density [kg/m3]\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->rho);
  }
  if (io->outputRhoSigma) { /* print the density */
    fprintf(pFile, "rhoSigma [kg/m3]\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->rhoSigma);
  }

  if (io->outputMu) { /* viscosity */
    fprintf(pFile, "Viscosity [Pa-s]\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->mu);
  }

  if (io->outputCp) { /* specific heat */
    fprintf(pFile, "Specific-heat [J/(kg-K)]\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->cpmix);
  }
  if (io->outputSumYi) { /* specific heat */
    fprintf(pFile, "Sum Yi\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->sumYi);
  }

  if (io->outputPDF) {
    fprintf(pFile, "pdf\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flow->pdf);
  }

  if (driver->Q_LOSS > 0) {
    fprintf(pFile, "qdot\n");
    writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flame->srcQ);
  }

  /* write soot source terms if present */
  char sootMomName[7][32] = {"number-density","volume","surface","volume-var",
                               "vol-surf-var","surface-var","delta-weight"};
  char sootMomAbbrv[7][8] = {"M00", "M10", "M01", "M20", "M11", "M02", "N0"};
  // FIXME: Corrent units!
  char sootMomUnit[7][16] = {"1/m^3","m^3/m^3","m^2/m^3","m^6/m^3","m^5/m^3",
                                "m^4/m^3","1/m^3"};
  char sootMomSrcUnit[7][16] = {"1/(m^3-s)","m^3/(m^3-s)","m^2/(m^3-s)",
                                "m^6/(m^3-s)","m^5/(m^3-s)","m^4/(m^3-s)",
                                "1/(m^3-s)"};
  char sootSrcNames[6][32] = {"nucleation","condensation","coagulation",
                              "surface-growth","oxidation","fragmentation"};
  int  sI[7] = {0, 1, 2, 3, 4, 5, 6};

  if (driver->sootMode > 0) {
    
    if (driver->sootMode<4) {
    if (flame->nSootMom < 5)
      sI[3] = 6;

    /* ----- moments ----- */
    for (i=0; i<flame->nSootMom; i++) {
      fprintf(pFile,"soot-%s_%s [%s]\n",sootMomName[sI[i]],sootMomAbbrv[sI[i]],
              sootMomUnit[sI[i]]);
      writeArray(pFile,flame->n2,flame->n1,flame->nSootMom,i,perLine,
                 flame->sootMom);
    }

    /* ----- moment source terms ----- */
    for (i=0; i<flame->nSootMom; i++) {
      fprintf(pFile,"soot-%s-src_%s [%s]\n",sootMomName[sI[i]],
              sootMomAbbrv[sI[i]], sootMomSrcUnit[sI[i]]);
      writeArray(pFile,flame->n2,flame->n1,flame->nSootMom,i,perLine,
                 flame->sootMomSrc);
    }
    if (io->outputSootMxmyp) { /* soot source term budget */
      for (i=0; i<flame->nSootMom-1; i++) {
        fprintf(pFile,"soot-Mxmyp_%i \n",i);
        writeArray(pFile,flame->n2,flame->n1,flame->nSootMom-1,i,perLine,
                   flame->sootMomxmyp);
      }
      for (i=0; i<flame->nSootMom; i++) {
        fprintf(pFile,"soot-%s-diff_%s [%s]\n",sootMomName[sI[i]],
                sootMomAbbrv[sI[i]], sootMomSrcUnit[sI[i]]);
        writeArray(pFile,flame->n2,flame->n1,flame->nSootMom,i,perLine,
                   flame->sootMomTransDiff);
      }
      for (i=0; i<flame->nSootMom; i++) {
        fprintf(pFile,"soot-%s-conv_%s [%s]\n",sootMomName[sI[i]],
                sootMomAbbrv[sI[i]], sootMomSrcUnit[sI[i]]);
        writeArray(pFile,flame->n2,flame->n1,flame->nSootMom,i,perLine,
                   flame->sootMomTransConv);
      }
    }
          
    if (io->outputSootSrc) { /* soot source term budget */
      for (j=0; j<6; j++) {
        for (i=0; i<flame->nSootMom-1; i++) {
          fprintf(pFile,"soot-%s-%s_%s [%s]\n",sootMomName[i],
                  sootSrcNames[j],sootMomAbbrv[i],sootMomSrcUnit[i]);
          writeArray(pFile,flame->n2,flame->n1,flame->nSootMom*6,
                     j*flame->nSootMom+i,perLine,flame->sootMomSrcPart);
        }
        fprintf(pFile,"soot-%s-%s_%s [%s]\n",sootMomName[6],sootSrcNames[j],
                sootMomAbbrv[6],sootMomSrcUnit[6]);
        writeArray(pFile,flame->n2,flame->n1,flame->nSootMom*6,
                   j*flame->nSootMom+(flame->nSootMom-1),perLine,
                   flame->sootMomSrcPart);
      }
    }
    for (i=0; i<flame->nSootMom; i++) {
      fprintf(pFile,"soot-%s-trans_%s [%s]\n",sootMomName[sI[i]],
              sootMomAbbrv[sI[i]], sootMomSrcUnit[sI[i]]);
      writeArray(pFile,flame->n2,flame->n1,flame->nSootMom,i,perLine,
                 flame->sootMomTrans);
    }
    }
  }

  /* ----- write trailer ----- */
  fprintf(pFile, "trailer\n");
  
  fclose(pFile);
}
/* -------------------------------------------------------------------------- */

void writeFlowData(char *fileName, ioStruct *io, driveStruct *driver,
                   timeStruct *time, flameStruct *flame, flowStruct *flow)
{
  int i, j, k, r, row, col, perLine = 5;
  char dim1, dim2;
  
  // assign the appropriate co-ordinate names
  if (driver->H_DIM == 0)
    {
      dim1 = 'Z';
      dim2 = 'Z';
    }
  else if (driver->H_DIM == 1)
    {
      dim1 = 'H';
      dim2 = 'Z';
    }
  else if (driver->H_DIM == 2)
    {
      dim1 = 'Z';
      dim2 = 'H';
    }
  
  FILE   *pFile;                 // output stream
  pFile = fopen(fileName, "w"); // open output file
  
  /* ----- write header ----- */
  fprintf(pFile, "header\n");
  fprintf(pFile, "title = \"%s\"\n", io->title);
  fprintf(pFile, "author = \"%s\"\n", io->author);
  fprintf(pFile, "date = \"%s\"\n\n", io->date);
  
  fprintf(pFile, "time = %15.6e [s]\n", time->t);
  fprintf(pFile, "pressure = %15.6e [Pa]\n", flow->P);
  fprintf(pFile, "DPDT     = %15.6e [Pa/s]\n", flow->Dt_P);
  fprintf(pFile, "%c1_min = %15.6e\n", dim1, flow->phi1min);
  fprintf(pFile, "%c1_max = %15.6e\n", dim1, flow->phi1max);
  if (flame->nDim == 2)
    {
      fprintf(pFile, "%c2_min = %15.6e\n", dim2, flow->phi2min);
      fprintf(pFile, "%c2_max = %15.6e\n", dim2, flow->phi2max); 
    }
  fprintf(pFile, "\n");

  fprintf(pFile, "gridPoints = %i\n", flame->n1);
  if (flame->nDim == 2)
    fprintf(pFile, "gridPoints2 = %i\n", flame->n2);
  fprintf(pFile, "\n");
  
  /* ----- write body ----- */
  fprintf(pFile, "body\n");
  // print the mesh
  fprintf(pFile, "%c1\n", dim1);
  writeArray(pFile, 1, flame->n1, 1, 0, perLine, flame->x1);
  if (flame->nDim == 2)
    {
      fprintf(pFile, "%c2\n", dim2);
      writeArray(pFile, 1, flame->n2, 1, 0, perLine, flame->x2);
    }
  
  // print the scalar dissipation rate
  if (flow->chi1 != NULL)
    {
      fprintf(pFile, "chi1 [1/s]\n");
      writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flow->chi1);
    }
  if (flame->nDim == 2)
    {
      if (flow->chi12 != NULL)
        {
          fprintf(pFile, "chi12 [1/s]\n");
          writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flow->chi12);
        }
      if (flow->chi2 != NULL)
        {
          fprintf(pFile, "chi2 [1/s]\n");
          writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flow->chi2);
        }
    }

  // print the pdf
  if (flow->pdf != NULL)
    {
      fprintf(pFile, "pdf \n");
      writeArray(pFile, flame->n2, flame->n1, 1, 0, perLine, flow->pdf);
    }

  if (flow->avgSdot != NULL)
    {
      if (driver->H_DIM == 1)
        {
          fprintf(pFile, "Sdot2 [J/(kg-s)]\n");
          writeArray(pFile, flame->n2, 1, 1, 0, perLine, flow->avgSdot);
        }
      else if (driver->H_DIM == 2)
        {
          fprintf(pFile, "Sdot1 [J/(kg-s)]\n");
          writeArray(pFile, 1, flame->n1, 1, 0, perLine, flow->avgSdot);
        }
    }

  /* ----- write trailer ----- */
  fprintf(pFile, "trailer\n");
  
  fclose(pFile);
}
/* -------------------------------------------------------------------------- */

int readTabFileVar(char *fileName, int iT, double **time, int iV, double **var)
{
  int   i, j, len, nHeader, firstChar, nCol;
  double tmp;
  
  /* strings to store the tag for each line and comments (if any) */
  char  tag[50], comment[1000], *col;
  
  FILE  *pFile = fopen(fileName, "r+");
  fpos_t pos;

  if (!pFile) {
    fprintf(stderr, "Error: monitor file not found: %s\n", fileName);
    exit(2);
  }

  len = nHeader = 0;
  /* get the number of data lines */
  while (!feof(pFile)) {
    strcpy(comment, "\n");
    fgets(comment, 1000, pFile);
    /* assume header if doesn't start with digit (or '-' for negative number) */
    firstChar = strspn (comment," ");
    if (!isdigit(comment[firstChar]) && comment[firstChar] != '-' && 
        comment[firstChar] != '\n')
      nHeader++;
    else if (strcmp(comment, "\n") && comment[0] != '\"')
      len++;
  }
  /* len -= 2; */
  fprintf(stderr, "%i headers found (%i)\n", nHeader, len);
  rewind(pFile);

  /* allocate the memory */
  time[0] = (double*) calloc(len, sizeof(double));
  var[0]  = (double*) calloc(len, sizeof(double));

  /* get the header lines */
  for (i=0; i<nHeader; i++)
    fgets(comment, 1000, pFile);

  /* get the number of columns */
  fgetpos(pFile, &pos);  /* get the of the buffer at the line start */
  fgets(comment, 1000, pFile);
  fsetpos(pFile, &pos);
  nCol = 0;
  col = strtok(comment, " ");
  while (col != NULL) {
    nCol++;
    col = strtok(NULL, " ");
  }
  fprintf(stderr, "columns: %i\n", nCol);

  /* read in the variables from the file */
  for (i=0; i<len; i++) {
    fgetpos(pFile, &pos);  /* get the of the buffer at the line start */
    /* read in variables until time index */
    for (j=0; j<iT-1; j++) tmp = readReal(pFile);
    time[0][i] = readReal(pFile);
    fgets(comment, 1000, pFile); /* read in rest of line */
    
    /* reset file to line start */
    fsetpos(pFile, &pos);
    for (j=0; j<iV-1; j++) tmp = readReal(pFile);
    var[0][i]  = readReal(pFile);
    if (iV < nCol) fgets(comment, 1000, pFile); /* read in rest of line */
  }
  
  fclose(pFile);

  return len;
}
/* ------------------------------------------------------------------------- */

void addVarToBody(char *fileName, int n1, int n2, double *var, char *varName)
{
  int  perLine = 5;
  char tag[100], *trailBuffer;
  fpos_t trail;
  long int fileSize, bodySize, trailSize;
  
  FILE *pFile = fopen(fileName, "r+");

  if (!pFile)
    {
      printf("Can't add variable to %s, file does not exist!\n", fileName);
      return;
    }

  // find the end of the body of file
  while (!feof(pFile))
    { 
      fscanf(pFile, "%s", tag);  // read in string

      if (!strcmp(tag, "trailer"))
        break;
    }

  // get the total length of the header/body in bytes
  bodySize = ftell(pFile) - 7;

  // get the total length of the file
  fseek(pFile, 0, SEEK_END);
  fileSize = ftell(pFile);
  trailSize = fileSize - bodySize; // size of trailer to read in

  // set position to point before trailer
  fseek(pFile, bodySize, SEEK_SET);

  // read in the whole trailer
  trailBuffer = (char*) calloc(trailSize, sizeof(char));
  fread(trailBuffer, 1, trailSize, pFile);
  
  // go back to beginning of trailer
  fseek(pFile, bodySize, SEEK_SET);

  // write the variable name
  fprintf(pFile, "%s\n", varName);

  // write the variable array at this point
  writeArray(pFile, n2, n1, 1, 0, perLine, var);

  // re-write the trailer at the end of the file
  fwrite(trailBuffer, 1, trailSize, pFile);

  free(trailBuffer);
  fclose(pFile);
}
/* -------------------------------------------------------------------------- */

void writeSolution(char *fileName, ioStruct *io, chemStruct *chem,
                   flameStruct *flame)
{
  int i, j, k;
  int row, col;
  double **H;
  
  /* ----- write data file ----- */
  FILE *pFile;
  pFile  = fopen(fileName, "w");
  
  if (flame->nDim == 2)
    fprintf(pFile, "%5s","x2 ");
  fprintf(pFile, "%5s","x1");
  for (i=0; i < io->nOutputs; i++)
    fprintf(pFile, " Y_%s ", io->outputSpecies[i]);
  fprintf(pFile, " Temp[K] Htot[J/kg]\n");
  
  // this should be the default, but this will ensure that the loop below works
  if (flame->nDim == 1)
    flame->n2 = 1;
  
  H  = (double**) calloc(flame->n2, sizeof(double*));
  for (i=0; i<flame->n2; ++i) 
    H[i] = (double*) calloc(flame->n1, sizeof(double));  
  
  /* compute the enthalpy field */
  /* ****** IMPROVEMENT: store this or compute using a function */
  for (k=0; k<flame->n2; k++) {
    row = k*(flame->n1*flame->nVars);
    for (j=0; j<flame->n1; j++) {
      col = j*flame->nVars;
      H[k][j] = 0.0;
      chemComputeThermoData(chem->h,chem->cp,flame->Y[row+col+flame->sT],chem);
      for (i=0; i<chem->nSpecies; i++)
        H[k][j] += chem->h[i]*flame->Y[row+col+i];
    }
  }
  
  for (k=0; k<flame->n2; k++)
    {
      row = k*(flame->n1*flame->nVars); // start of current row
      for (j=0; j<flame->n1; j++)
        {
          col = j*flame->nVars; // current row/column position
          
          fprintf(pFile, "%6.4f ", flame->x1[j]);
          if (flame->nDim == 2)
            fprintf(pFile, "%6.4f ", flame->x2[k]);
          for (i=0; i<io->nOutputs; i++)
            fprintf(pFile, "%15.6e ", flame->Y[row+col+io->outputIndex[i]]);
          fprintf(pFile,"%15.6e %15.6e\n",flame->Y[row+col+flame->sT],H[k][j]);
        } // end loop over x1
      fprintf(pFile, "\n\n"); // spacing between data blocks
    } // end loop over x2
  
  fclose(pFile);
  
  for (i=0; i<flame->n2; ++i) free(H[i]);
  free(H);
}
/* -------------------------------------------------------------------------- */

void writeInitSolution(char *fileName, ioStruct *io, flameStruct *flame)
{
  int i, j, firstY;
  char fileHead[100], gnuPlot[200], *head;
  strcpy(fileHead, fileName);
    
  /* ----- write data file ----- */
  FILE *pFile;
  pFile  = fopen(fileName, "w");
  
  fprintf(pFile, "%5s %10s","x1","temp_[k]");
  for (i=0; i < io->nInputs; i++)
    fprintf(pFile, " Y_%s ", io->inputSpecies[i]);
  fprintf(pFile, "\n");
  for (j=0; j<flame->n1; j++) 
    {
      fprintf(pFile, "%5.3f %10.3f ", flame->x1[j], io->initT[j]);
      for (i=0; i < io->nInputs; i++)
        fprintf(pFile, "%10.3e ", io->initY[i][j]);
      fprintf(pFile, "\n");
    }
  
  fclose(pFile);
  
  /* ----- write gnuplot script ----- */
  // create file name to generate gnuplot script to plot data in this file
  head = strtok(fileHead, ".");
  strcpy(gnuPlot, "plot_");
  strcat(gnuPlot, head);
  strcat(gnuPlot, ".gp");
  
  FILE *gpFile;
  gpFile = fopen(gnuPlot, "w");
  
  // write header info
  fprintf(gpFile, "set terminal X11 1 persist\n");
  fprintf(gpFile, "set title \"initial solution read from flamelet\"\n");
  fprintf(gpFile, "set ylabel \"mass fraction\"\n");
  fprintf(gpFile, "set xlabel \"x1\"\n");
  fprintf(gpFile, "set key out vert center right\n");
  fprintf(gpFile, "\n");
  
  // write plotting info for species
  firstY = 3; // location of first species
  fprintf(gpFile,"plot '%s' using 1:%i title %i, \\\n",fileName,firstY,firstY);
  for (i=firstY+1; i < firstY+io->nInputs; i++)
    {
      fprintf(gpFile, "      '' using 1:%i title %i", i, i);
      if (i == firstY + io->nInputs - 1)
        fprintf(gpFile, "\n");
      else
        fprintf(gpFile, ", \\\n");
    }
  
  // write plotting info for temperature
  fprintf(gpFile, "\n");
  fprintf(gpFile, "set terminal X11 2 persist\n");
  fprintf(gpFile, "set title \"initial temperature profile from flamelet\"\n");
  fprintf(gpFile, "set ylabel \"T [K]\"\n");
  fprintf(gpFile, "set xlabel \"x1\"\n");
  fprintf(gpFile, "set key inside vert bottom right\n");
  fprintf(gpFile, "\n");
  fprintf(gpFile, "plot '%s' using 1:2 title \"Temperature\"\n", fileName);
  
  fclose(gpFile);
}
/* -------------------------------------------------------------------------- */

void writeArray(FILE *pFile, int m, int n, int offset, int index, int newLine,
                double *data)
{
  int j, k, row, col, curLine;
  
  for (k = 0; k < m; ++k)
    {
      curLine = newLine;   // reset the counter for a newline
      row = k*n*offset;    // set the row index
      for (j = 0; j < n; ++j)
        {
          col = j*(offset); // set the column index
          fprintf(pFile, "%15.6e", data[row+col+index]);
          curLine--;
          if (curLine == 0)
            {
              fprintf(pFile, "\n");
              curLine = newLine;
            }
        }
      if (curLine != newLine)
        fprintf(pFile, "\n");
    }
}
/* -------------------------------------------------------------------------- */

void writeScalarsHeader(driveStruct *driver, ioStruct *io)
{
  int i, col;

  FILE *pFile = fopen(io->scalars, "w");
  if (pFile) {
    fprintf(pFile, "# time[s]<1> ");
    col = 2;
    if (driver->TIME_MODE == 3) { /* engine */
      fprintf(pFile, "# CAD[deg]<%i> ", col); col++;
    }
    for (i=0; i<io->nOutputs; ++i, ++col)
      fprintf(pFile, "Y_%s<%i> ", io->outputSpecies[i], col);
    fprintf(pFile, "Y_sum<%i> ", col);        col++;
    fprintf(pFile, "T[K]<%i> ", col);         col++;
    fprintf(pFile, "Pressure[Pa]<%i> ", col); col++;
    if (io->outputHtot) {
      fprintf(pFile, "TotalEnthalpy[J/(kg-K)]<%i> ", col); col++;
    }
    if (io->outputChi) {
      fprintf(pFile, "Chi[1/s]<%i> ", col); col++;
    }
    if (io->outputHR) {
      fprintf(pFile, "Heat-Release[J/(kg-s)]<%i> ", col); col++;
    }
    if (io->outputRho) {
      fprintf(pFile, "Rho[kg/m3]<%i> ", col); col++;
    }
    if (io->outputCp) {
      fprintf(pFile, "Cp[J/(kg-K)]<%i> ", col); col++;
    }
    if (io->outputQdot) {
      fprintf(pFile, "Qdot[J/kg]<%i> ", col);      col++;
      fprintf(pFile, "integralQdot[J]<%i> ", col); col++;
    }
    fprintf(pFile, "1/rho_0*(dP/dt)[J/kg-s]<%i> ", col); col++;
    fprintf(pFile, "dHmin/dt[J/kg-s]<%i> ", col); col++;
    fprintf(pFile, "1/rho_1*(dP/dt)[J/kg-s]<%i> ", col); col++;
    fprintf(pFile, "dHmax/dt[J/kg-s]<%i> ", col); col++;
    fprintf(pFile, "dPhiMin/dt[J/kg-s]<%i> ", col); col++;
    fprintf(pFile, "dPhiMax/dt[J/kg-s]<%i> ", col); col++;
    fprintf(pFile, "d<PhiMin>/dt[J/kg-s]<%i> ", col); col++;
    fprintf(pFile, "d<PhiMax>/dt[J/kg-s]<%i> ", col); col++;
    fprintf(pFile, "<dPhiMin/dt>[J/kg-s]<%i> ", col); col++;
    fprintf(pFile, "<dPhiMax/dt>[J/kg-s]<%i> ", col); col++;
    if (io->outputGamma) {
      fprintf(pFile, "gamma<%i> ", col); col++;
    }
    if (io->outputDryNOX) {
      fprintf(pFile, "NO_dry_6perc<%i> ", col); col++;
    }
    if (driver->DNS) {
      fprintf(pFile, "dPdt<%i> ", col);           col++;
      fprintf(pFile, "d(rho*Cp*T)/dt<%i> ", col); col++;
      fprintf(pFile, "diff1[1/s]<%i> ", col);     col++;
      fprintf(pFile, "diff12[1/s]<%i> ", col);    col++;
      fprintf(pFile, "diff2[1/s]<%i> ", col);     col++;
    }
    fprintf(pFile, "\n");
    fclose(pFile);
  }
  else {
    fprintf(stderr, "ERROR: could not open file %s\n", io->scalars);
    exit(2);
  }
}
/* -------------------------------------------------------------------------- */

int truncateFileTime(char *fileName, double t_end)
{
  int i;
  float t_line;

  FILE *pFile = NULL, *buFile;
  size_t check_read;
  char buffer[5000];
  long sizeData;
  char *data;
  /* TODO: figure a better way of reading lines to avoid large buffer */

  pFile = NULL;
  pFile = fopen(fileName, "r");
  if (!pFile) /* return if file doesn't exist */
    return 0;

  /* ----- make a backup of the current file ----- */
  /* read in all the data from existing file */
  fseek(pFile,0,SEEK_END);
  sizeData = ftell(pFile);
  data = (char*) calloc(sizeData, sizeof(char));
  fseek(pFile,0,SEEK_SET); /* rewind file */
  fread(data,1,sizeData,pFile);

  /* write the data to the backup file */
  sprintf(buffer, "%s~", fileName);
  buFile = fopen(buffer, "w");
  fwrite(data,1,sizeData,buFile);
  fclose(buFile);
  free(data);

  /* ----- locate the end of the data before current time ----- */
  /* make sure file is at start */
  fseek(pFile,0,SEEK_SET);
  t_line = -1;
  while (!feof(pFile))
    {
      sizeData = ftell(pFile);
      fgets(buffer, 5000, pFile);

      if (buffer[0] != '#' && buffer[0] != '\n')
        sscanf(buffer, "%f %*s", &t_line);
      t_line = (double)(t_line);

      if ((t_line >= t_end))
        break;
    } /* end loop over file */

  /* go to the beginning of the file */
  fseek(pFile, 0, SEEK_SET);

  /* read in everything up until the current time step */
  data = (char*) calloc(sizeData, sizeof(char));
  check_read = fread(data,1,sizeData,pFile);

  /* reopen the file and write the previous information */
  freopen(fileName, "w", pFile);
  fwrite(data,1,sizeData,pFile);

  /* free up data and close file */
  free(data);
  fclose(pFile);
  return 1;
}
/* -------------------------------------------------------------------------- */

void stripHeader(char *str)
{
  char *head, *species;
  
  head    = strtok(str, "-"); // read filename till -
  species = strtok(NULL, ""); // get the species name after the dash
  if (species == NULL)
    strcpy(str, head);        // there was no header!
  else
    strcpy(str, species);     // return species name
}
/* -------------------------------------------------------------------------- */

void NGA_printFileInfo(char *file)
{
  int    i, nx, ny, nz, nVars;
  double t, dt;
  char   **varNames;

  // read in the info from the file
  FILE *pFile = fopen(file, "rb");

  // check if file exists
  if (!pFile)
    {
      fprintf(stderr, "Error: file %s could not be opened\n", file);
      exit(2);
    }

  // read in the header info
  NGA_readInfo(pFile, &nx, &ny, &nz, &nVars, &varNames, &t, &dt);

  fclose(pFile);

  // print info to screen
  fprintf(stderr, "\n----- NGA Data File Info -----\n\n");
  fprintf(stderr, "File:      %s\n", file);
  fprintf(stderr, "Grid Size: %i x %i x %i\n", nx, ny, nz);
  fprintf(stderr, "time:      %8.6f\n", t);
  fprintf(stderr, "dt:        %8.6f\n", dt);
  fprintf(stderr, "\n%i Variables: ", nVars);
  for (i=0; i<nVars; i++) fprintf(stderr, "%s ", varNames[i]);
  fprintf(stderr, "\n\n-------------------------------\n");

  for (i=0; i<nVars; i++) free(varNames[i]);
  free(varNames);
}
/* -------------------------------------------------------------------------- */

void NGA_readInfo(FILE *pFile, int *nx, int *ny, int *nz, int *nVars,
                  char ***varNames, double *t, double *dt)
{
  int i;

  // read in grid size data and number of variables
  fread(nx, sizeof(int), 1, pFile);
  fread(ny, sizeof(int), 1, pFile);
  fread(nz, sizeof(int), 1, pFile);
  fread(nVars, sizeof(int), 1, pFile);

  // read in time variables
  fread(dt, sizeof(double), 1, pFile);
  fread(t,  sizeof(double), 1, pFile);

  // read variable names
  varNames[0] = (char**) calloc(*nVars, sizeof(char*));
  for (i=0; i<*nVars; i++)
    {
      varNames[0][i] = (char*) calloc(8, sizeof(char));
      fread(varNames[0][i], sizeof(char), 8, pFile);
      stripString(varNames[0][i]);
    }
}
/* -------------------------------------------------------------------------- */

void NGA_readData(char *file, int *nx, int *ny, int *nz, char *varName,
                  double ****var, double *t, double *dt)
{
  int i,j,k, nVars, iVar;
  long int varStart;
  char **names;

  double dummy;
  
  FILE *pFile = fopen(file, "rb");
  
  if (!pFile)
    {
      fprintf(stderr, "Error: file %s could not be opened\n", file);
      exit(2);
    }

  // read in the header info
  NGA_readInfo(pFile, nx, ny, nz, &nVars, &names, t, dt);

  // make sure the variable name is valid and get its index
  iVar = -1;
  for (i=0; i<nVars; i++)
    if (!strcmp(varName, names[i]))
      iVar = i;

  if (iVar < 0)
    {
      fprintf(stderr, "Error: variable %s does not exist\n", varName);
      exit(2);
    }

  // go to the position of the variable in the file
  fseek(pFile, iVar*(*nx)*(*ny)*(*nz)*sizeof(double), SEEK_CUR);

  // allocate the memory and read the variable
  var[0]  = (double***) calloc(*nz, sizeof(double**));
  for (k=0; k<*nz; k++)
    {
      var[0][k] = (double**) calloc(*ny, sizeof(double*));
      for (j=0; j<*ny; j++)
        {
          var[0][k][j] = (double*) calloc(*nx, sizeof(double));

          for (i=0; i<*nx; i++)
            fread(&var[0][k][j][i], sizeof(double), 1, pFile);
        }  // end loop over y
    } // end loop over z

  for (i=0; i<nVars; i++)
    free(names[i]);
  free(names);

  fclose(pFile);
}
/* -------------------------------------------------------------------------- */

void NGA_readDataAll(char *file, int *nx, int *ny, int *nz, int *nVars,
                     char ***varNames, double *****vars, double *t, double *dt)
{
  int i, j, k, m;
  
  FILE *pFile = fopen(file, "rb");
  
  if (!pFile)
    {
      fprintf(stderr, "Error: file %s could not be opened\n", file);
      exit(2);
    }

  // read in the header info
  NGA_readInfo(pFile, nx, ny, nz, nVars, varNames, t, dt);

  // allocate memory and read data
  vars[0] = (double****) calloc(*nVars, sizeof(double***));
  for (m=0; m<*nVars; m++)
    {
      vars[0][m] = (double***) calloc(*nz, sizeof(double**));
      for (k=0; k<*nz; k++)
        {
          vars[0][m][k] = (double**) calloc(*ny, sizeof(double*));
          for (j=0; j<*ny; j++)
            {
              vars[0][m][k][j] = (double*) calloc(*nx, sizeof(double));
              for (i=0; i<*nx; i++)
                {
                  fread(&vars[0][m][k][j][i], sizeof(double), 1, pFile);
                } // end loop over x
            }  // end loop over y
        } // end loop over z
    } // end loop over variables

  fclose(pFile);
}
/* -------------------------------------------------------------------------- */

int NIST_readThermo(char *fileName, double **hfg, double **cpG, double **cpL,
                    double **Psat, double **T)
{
  int i, j, len;
  double HL, HG, dT;
  int nVars, sT, sPsat, sCpL, sCpG, sHG, sHL;

  char comment[800];
  float rTmp;
  double line[25];

  nVars = 25;
  sT    = 0;
  sPsat = 1;
  sHL   = 5;
  sCpL  = 8;
  sHG   = 17;
  sCpG  = 20;

  FILE *pFile = fopen(fileName, "r");

  if (!pFile)
    {
      fprintf(stderr, "Error: thermo file %s not found in getThermoNIST\n",
              fileName);
      exit(1);
    }

  len = 0;
  // get the number of data lines
  while (!feof(pFile))
    {
      strcpy(comment, "\n");
      fgets(comment, 1000, pFile);
      if (strcmp(comment, "\n"))
        len++;
    }
  len -= 1;
  rewind(pFile);

  /* allocate memory */
  T[0]    = (double*) calloc(len, sizeof(double));
  Psat[0] = (double*) calloc(len, sizeof(double));
  cpG[0]  = (double*) calloc(len, sizeof(double));
  cpL[0]  = (double*) calloc(len, sizeof(double));
  hfg[0]  = (double*) calloc(len, sizeof(double));

  // get the header line
  fgets(comment, 800, pFile);

  for (j=0; j<len; j++)
    {
      /* scan in the line */
      for (i=0; i<25; i++)
        {
          fscanf(pFile, "%f", &rTmp);
          line[i] = (double)(rTmp);
        }

      /* store relevant data */
      T[0][j]    = line[sT];
      Psat[0][j] = line[sPsat]*1.0e6;
      cpG[0][j]  = line[sCpG];
      cpL[0][j]  = line[sCpL];
      hfg[0][j]  = (line[sHG]-line[sHL])*1.0e3;
    }

  fclose(pFile);
  return len;
}
/* ------------------------------------------------------------------------- */

void NIST_getThermo(char *fileName, double *hfg, double *cpG, double *cpL,
                    double *Psat, double T)
{
  int i;
  double HL, HG, dT;
  int nVars, sT, sPsat, sCpL, sCpG, sHG, sHL;

  nVars = 25;
  sT    = 0;
  sPsat = 1;
  sHL   = 5;
  sCpL  = 8;
  sHG   = 17;
  sCpG  = 20;

  FILE *pFile = fopen(fileName, "r");
  char comment[800];
  float rTmp;
  double line[25], line_m1[25];
  fpos_t pos;

  // get the header line
  fgets(comment, 800, pFile);

  // read the first line
  for (i=0; i<nVars; i++)
    {
      fscanf(pFile, "%f", &rTmp);
      line_m1[i] = (double)(rTmp);
    }
  if (T == line_m1[0])
    {
      *Psat  = line_m1[sPsat]*1e6;
      *cpL = line_m1[sCpL];
      *cpG = line_m1[sCpG];
      *hfg   = (line_m1[sHG]-line_m1[sHL])*1e3;
      return;
    }

  while (!feof(pFile))
    {
      fgetpos(pFile, &pos);
      for (i=0; i<25; i++)
        {
          fscanf(pFile, "%f", &rTmp);
          line[i] = (double)(rTmp);
        }

      if (T == line[0])
        {
          *Psat = line[sPsat]*1e6;
          *cpL  = line[sCpL];
          *cpG  = line[sCpG];
          *hfg  = (line[sHG]-line[sHL])*1e3;
          return;
        }
      else if (T < line[0] && T > line_m1[0])
        {
          dT    = line[sT]-line_m1[sT];
          *Psat = line_m1[sPsat] + (line[sPsat]-line_m1[sPsat])/dT;
          *cpL  = line_m1[sCpL]  + (line[sCpL]-line_m1[sCpL])/dT;
          *cpG  = line_m1[sCpG]  + (line[sCpG]-line_m1[sCpG])/dT;
          HG    = line_m1[sHG]   + (line[sHG]-line_m1[sHG])/dT;
          HL    = line_m1[sHL]   + (line[sHL]-line_m1[sHL])/dT;
          *Psat = *Psat*1e6;
          *hfg  = (HG - HL)*1e3;
          return;
        }

      for (i=0; i<25; i++) line_m1[i] = line[i];
    }

  printf("ERROR: temperature range not sufficient for T = %8.2f\n", T);
  fclose(pFile);
}
/* ------------------------------------------------------------------------- */

void initEnsight(ioStruct *io, flowStruct *flow)
{
  int i;
  struct stat st;
  char buffer[100];

  // check to see if the directory exists
  if (stat(io->ensightDir, &st) != 0) mkdir(io->ensightDir, 0755);

  // species to plot
  for (i=0; i<io->nInputs; i++)
    {
      strcpy(buffer, io->ensightDir);
      strcat(buffer, "/");
      strcat(buffer, io->inputSpecies[i]);
      if (stat(buffer, &st) != 0) mkdir(buffer, 0755);
    }

  // temperature
  if (io->initT != NULL)
    {
      strcpy(buffer, io->ensightDir);
      strcat(buffer, "/T");
      if (stat(buffer, &st) != 0) mkdir(buffer, 0755);
    }

  // probability distribution
  if (flow->pdf != NULL)
    {
      strcpy(buffer, io->ensightDir);
      strcat(buffer, "/PDF");
      if (stat(buffer, &st) != 0) mkdir(buffer, 0755);
    }

  // scalar dissipation rates
  if (flow->chi1 != NULL)
    {
      strcpy(buffer, io->ensightDir);
      strcat(buffer, "/CHI1");
      if (stat(buffer, &st) != 0) mkdir(buffer, 0755);
    }
  if (flow->chi12 != NULL)
    {
      strcpy(buffer, io->ensightDir);
      strcat(buffer, "/CHI12");
      if (stat(buffer, &st) != 0) mkdir(buffer, 0755);
    }
  if (flow->chi2 != NULL)
    {
      strcpy(buffer, io->ensightDir);
      strcat(buffer, "/CHI2");
      if (stat(buffer, &st) != 0) mkdir(buffer, 0755);
    }

  // miscellaneous variables
  for (i=0; i<io->nMisc; i++)
    {
      strcpy(buffer, io->ensightDir);
      strcat(buffer, "/");
      strcat(buffer, io->miscNames[i]);
      if (stat(buffer, &st) != 0) mkdir(buffer, 0755);
    }
}
/* -------------------------------------------------------------------------- */

void writeEnsightCaseHead(ioStruct *io, driveStruct *driver, timeStruct *time,
                          flameStruct *flame, flowStruct *flow)
{
  int i;
  char  caseFile[100], buffer[100];

  strcpy(caseFile, io->ensightDir);
  strcat(caseFile, "/mduc.case");

  FILE *pFile = fopen(caseFile, "w");

  fprintf(pFile, "FORMAT\n");
  fprintf(pFile, "type: ensight gold\n");
  fprintf(pFile, "GEOMETRY\n");
  fprintf(pFile, "model: geometry\n");
  fprintf(pFile, "VARIABLE\n");

  for (i=0; i<io->nInputs; i++)
    fprintf(pFile, "scalar per node: 1 %s %s/%s.******\n",
            io->inputSpecies[i], io->inputSpecies[i], io->inputSpecies[i]);

  if (io->initT != NULL)
    fprintf(pFile, "scalar per node: 1 Temperature T/T.******\n");

  if (flow->pdf != NULL)
    fprintf(pFile, "scalar per node: 1 PDF PDF/PDF.******\n");

  if (flow->chi1 != NULL)
    fprintf(pFile, "scalar per node: 1 CHI1 CHI1/CHI1.******\n");
  if (flow->chi12 != NULL)
    fprintf(pFile, "scalar per node: 1 CHI12 CHI12/CHI12.******\n");
  if (flow->chi2 != NULL)
    fprintf(pFile, "scalar per node: 1 CHI2 CHI2/CHI2.******\n");

  for (i=0; i<io->nMisc; i++)
    fprintf(pFile, "scalar per node: 1 %s %s/%s.******\n", io->miscNames[i],
            io->miscNames[i], io->miscNames[i]);
  
  fprintf(pFile, "TIME\n");
  fprintf(pFile, "time set: 1\n");
  fprintf(pFile, "number of steps: %i\n", 0);
  fprintf(pFile, "filename start number: %i\n", 1);
  fprintf(pFile, "filename increment: %i\n", 1);
  fprintf(pFile, "time values: ");
  fprintf(pFile, "\n");

  fclose(pFile);
}
/* -------------------------------------------------------------------------- */

void updateEnsightCaseTime(double newTime, ioStruct *io, timeStruct *time)
{
  int i, nt, iStart, increment;
  long int headSize;
  char  caseFile[100], tag[50], buffer[100], *headBuffer, *ident;
  double *timeVals = NULL;

  strcpy(caseFile, io->ensightDir);
  strcat(caseFile, "/mduc.case");

  FILE *pFile = fopen(caseFile, "r");

  // find the beginning of the TIME section
  while (!feof(pFile))
    { 
      fscanf(pFile, "%s", tag);  // read in string

      if (!strcmp(tag, "TIME"))
        break;
    }

  // get the total length of the header/body in bytes
  headSize = ftell(pFile) - 4;

  // read in the time info
  fgets(buffer, 80, pFile); // time set line
  fgets(buffer, 80, pFile); // time set
  fgets(buffer, 80, pFile); // number of steps
  ident = strtok(buffer, ":");
  nt = atoi(strtok(NULL,""));
  fgets(buffer, 80, pFile); // file start number
  ident = strtok(buffer, ":");
  iStart = atoi(strtok(NULL,""));
  fgets(buffer, 80, pFile); // file increment
  ident = strtok(buffer, ":");
  increment = atoi(strtok(NULL,""));

  // "time values:"
  fgets(buffer, 12, pFile);

  // read in already existing time steps
  readRealArray(pFile, nt, &timeVals);
  
  // go back to the beginning of the file and read in the whole header
  rewind(pFile);
  headBuffer = (char*) calloc(headSize, sizeof(char));
  fread(headBuffer, 1, headSize, pFile);

  freopen(caseFile, "w", pFile);

  // re-write the header
  fwrite(headBuffer, 1, headSize, pFile);

  // update the time vectors
  nt++;
  timeVals = (double*) realloc(timeVals, nt*sizeof(double));
  timeVals[nt-1] = newTime;

  // rewrite the time info
  fprintf(pFile, "TIME\n");
  fprintf(pFile, "time set: 1\n");
  fprintf(pFile, "number of steps: %i\n", nt);
  fprintf(pFile, "filename start number: %i\n", iStart);
  fprintf(pFile, "filename increment: %i\n", increment);
  fprintf(pFile, "time values:\n");
  writeArray(pFile, 1, nt, 1, 0, 5, timeVals);

  fclose(pFile);

  free(timeVals);
  free(headBuffer);
}
/* -------------------------------------------------------------------------- */

void writeEnsightData(char *fileName, char *varName, int m, int n,
                      int offset, int index, double *data)
{
  int  iTmp = 1;
  char buffer[80];

  int j, k, row, col, curLine;
  float *tmp = (float*) calloc(m*n, sizeof(float));
  
  for (k = 0; k < m; ++k)
    {
      row = k*n*offset;    // set the row index
      for (j = 0; j < n; ++j)
        {
          col = j*(offset); // set the column index
          tmp[k*n + j] = (float)(data[row+col+index]);
        }
    }

  FILE *pFile = fopen(fileName, "wb");

  sprintf(buffer, "Variable file: %s", varName);
  fprintf(pFile, "%-80s", buffer);
  fprintf(pFile, "%-80s", "part");
  fwrite(&iTmp, sizeof(int), 1, pFile);
  fprintf(pFile, "%-80s", "block");
  fwrite(tmp, sizeof(float), n*m, pFile);

  fclose(pFile);

  free(tmp);
}
/* -------------------------------------------------------------------------- */

void writeEnsightGeom(ioStruct *io, driveStruct *driver,
                      timeStruct *time, flameStruct *flame, flowStruct *flow)
{
  int   iTmp = 1;
  char  geomFile[100];
  float fTmp = 0.0;

  strcpy(geomFile, io->ensightDir);
  strcat(geomFile, "/geometry");

  float xMin, xMax, yMin, yMax, zMin, zMax;
  xMin = 0.0;
  xMax = 1.0;
  yMin = 0.0;
  if (flame->nDim == 2)
    yMax = 1.0;
  else
    yMax = 0.0;
  zMin = 0.0;
  zMax = 0.0;

  // write geometry
  FILE *pFile = fopen(geomFile, "wb");

  fprintf(pFile, "%-80s", "C BINARY");
  fprintf(pFile, "%-80s", "Ensight Geometry File");
  fprintf(pFile, "%-80s", "Structured Geometry from MDUC");
  fprintf(pFile, "%-80s", "node id assign");
  fprintf(pFile, "%-80s", "element id assign");
  fprintf(pFile, "%-80s", "extents");
  fwrite(&xMin, sizeof(float), 1, pFile);
  fwrite(&xMax, sizeof(float), 1, pFile);
  fwrite(&yMin, sizeof(float), 1, pFile);
  fwrite(&yMax, sizeof(float), 1, pFile);
  fwrite(&zMin, sizeof(float), 1, pFile);
  fwrite(&zMax, sizeof(float), 1, pFile);
  fprintf(pFile, "%-80s", "part");
  fwrite(&iTmp, sizeof(int), 1, pFile);
  fprintf(pFile, "%-80s", "MDUC Grid");
  fprintf(pFile, "%-80s", "block rectilinear");
  fwrite(&flame->n1, sizeof(int), 1, pFile);
  fwrite(&flame->n2, sizeof(int), 1, pFile);
  fwrite(&iTmp, sizeof(int), 1, pFile);
  writeBinaryArrayReal(pFile, flame->n1, flame->x1);
  if (flame->nDim == 2)
    writeBinaryArrayReal(pFile, flame->n2, flame->x2);
  else
    fwrite(&fTmp, sizeof(float), 1, pFile);
  fwrite(&fTmp, sizeof(float), 1, pFile);

  fclose(pFile);
}
/* -------------------------------------------------------------------------- */

void writeBinaryArrayReal(FILE *pFile, int n, double *val)
{
  int i;
  float *tmp = (float*) calloc(n, sizeof(float));
  
  for (i=0; i<n; i++) tmp[i] = (float)(val[i]);
  fwrite(tmp, sizeof(float), n, pFile);

  free(tmp);
}
/* -------------------------------------------------------------------------- */

int binaryFileRead(void *buffer, int size, int count, FILE *pFile)
{
  int i,n;
  unsigned char *tmp;
  char *pos;

  tmp = (unsigned char*) malloc(size);
  pos = (char*) buffer;
   
  for (n=0; n<count; n++) {
    for (i=0; i<size; i++) {
      tmp[i] = fgetc(pFile);
    }
    memcpy(pos,tmp,size);
    pos += size;
  }
   
  free(tmp);
  return 0;
} 
/* ------------------------------------------------------------------------- */

int binaryFileWrite(void *buffer, int size, int count, FILE *pFile)
{
  int i,n;
  unsigned char *tmp;
  char *pos;

  tmp = (unsigned char*) malloc(size);
  pos = (char*) buffer;
   
  for (n=0; n<count; n++) {
    memcpy(tmp,pos,size);
    for (i=0; i<size; i++) {
      fputc(tmp[i], pFile);
    }
    pos += size;
  }
   
  free(tmp);
  return 0;
} 
/* ------------------------------------------------------------------------- */
