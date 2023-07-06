#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "accessMDUC.h"
#include "accessChem.h"
#include "chemInit.h"
#include "chemMech.h"
#include "chemTable.h"
#include "chemUtils.h"
#include "mducStructs.h"
#include "init.h"
#include "chemInit.h"
#include "destroy.h"
#include "fileIO.h"
#include "display.h"
#include "parser.h"
#include "flameletInit.h"
#include "flameletUtils.h"
#include "flamelet1D.h"
#include "flamelet2D.h"
#include "solve1D.h"
#include "solve2D.h"
#include "scalarDiss.h"
#include "metric.h"
#include "solveUtils.h"
#include "engine.h"
#include "engineUtils.h"
#include "utils.h"
#include "statsUtils.h"
#include "errFunc.h"
#include "sootHMOM.h"


/* sets MDUC_PRINT_PROC before MRIF Parallel Init
   when used with CIAO */
void MDUCtempPrintProc(int printproc)
{
  MDUC_PRINT_PROC=0;
  if (printproc==1) {
	  MDUC_PRINT_PROC = 1;
	}
}

void *MDUCinit(char *inputFileName)
{
  int j;
  char outputFileName[256];
  mducStruct *mduc = (mducStruct*) calloc(1, sizeof(mducStruct));
  
  /* allocate pointers to internal structs */
  mduc->io     = (ioStruct*)     calloc(1, sizeof(ioStruct));
  mduc->driver = (driveStruct*)  calloc(1, sizeof(driveStruct));
  mduc->chem   = NULL;
  mduc->time   = (timeStruct*)   calloc(1, sizeof(timeStruct));
  mduc->flow   = (flowStruct*)   calloc(1, sizeof(flowStruct));
  mduc->flame  = (flameStruct*)  calloc(1, sizeof(flameStruct));
  mduc->sdial  = (sdialStruct*)  calloc(1, sizeof(sdialStruct));
  mduc->engine = (engineStruct*) calloc(1, sizeof(engineStruct));
  mduc->rk     = (rkStruct*)     calloc(1, sizeof(rkStruct));
#ifdef ACCESS_MPI
  mduc->mpi    = (mpiStruct*)   calloc(1, sizeof(mpiStruct));
#endif /* ACCESS_MPI */

  /* local pointers */
  ioStruct     *io     = mduc->io;
  driveStruct  *driver = mduc->driver;
  chemStruct   *chem   = mduc->chem;
  timeStruct   *time   = mduc->time;
  flowStruct   *flow   = mduc->flow;
  flameStruct  *flame  = mduc->flame;
  sdialStruct  *sdial  = mduc->sdial;
  engineStruct *engine = mduc->engine;

  /* ----- set default values ----- */
  setDefaults(io, driver, time, flame, flow, engine);
  strcpy(sdial->solver, "cvode");
  strcpy(sdial->method, "gmres");
  strcpy(sdial->precnd, "left");

  /* store the input file name */
  strcpy(mduc->io->inputFile, inputFileName);
  
  /* ----- read input data file ----- */
  /* read into parser */
  io->parser = parseFile(inputFileName);
  /* TODO: should get rid of this function */
  readInitFile(io, driver, time, flame, flow, sdial, engine);

  /* start the monitor */
#ifndef ACCESS_MPI
  checkMonitor(io);
#endif

  /* set time to initial time */
  time->t   = time->t0;
  time->nt0 = time->nt;

  /* read in data if starting from existing solution */
  if (io->readData)
    readFlamelet(io->inputData, io, driver, time, flame, flow);

  if (io->readFlow)
    readFlowData(io->inputFlow, driver, flow);

  if (driver->ENERGY && flame->nDim == 2)
    mducErr("enthalpy energy equation not yet implemented in 2D","MDUCinit");

  /* ----- initialise parameters after reading input file ----- */
  mduc->chem = chemInit("compiled",io->parser); /* variables from mechanism */
  chem = mduc->chem;

  if (driver->output_w)
     if (io->readRates)
        readRates(io->inputRates,chem);

  
  initIO(io, driver, chem);           // input/output parameters
  initTime(driver, time);             // time variables
  initFlow(driver, flame, flow);
  /* initialise the main working memory and solution for the solver */
  initSolver(io, driver, chem, flame, flow);

  /* initialize engine parameters if necessary */
  initEngine(io, chem, engine);
  if (driver->ENGINE) {
    initEngineMass(mduc);
  }

  /* initialise internal SUNDIALS memory for solver */
  /* cross-reference some pointers for SUNDIALS arguments to be passed */
  sdial->args    = (void**) calloc(8, sizeof(void*));
  sdial->args[0] = driver;
  sdial->args[1] = chem;
  sdial->args[2] = flow;
  sdial->args[3] = flame;
  sdial->args[4] = time;
  sdial->args[5] = engine;
  sdial->args[6] = io;
#ifdef ACCESS_MPI
  sdial->args[7] = mduc->mpi;
#endif /* ACCESS_MPI */

  /* call the solver setup */
  initSUNDIALS(io, sdial);
  
  /* write solution file at this step (change output name here if desired) */
  if (!io->readData) {
    if (!driver->ENERGY && io->outputHtot)
      calcTotalEnthalpy(flame->n1,flame->n2,flame->nVars,flame->sT,flame->Y,
                        flame->H,chem);
    
    /* createDataFileName(io->outputData, io, time, outputFileName); */
    /* writeFlamelet(outputFileName, io, driver, time, flame, flow); */
  }
  /* parsePrint(io->parser); */
  if (driver->output_w)
     if (io->readRates && MDUC_PRINT_PROC)
       writeRates("restart_rates",chem);

  return (void*)mduc;
}
/* -------------------------------------------------------------------------- */

void *MDUCinitChem(char *inputFileName)
{
/*   if (inputFileName == NULL) { /\* only works for compiled mechanism *\/ */
    
    mducStruct *mduc = (mducStruct*) calloc(1, sizeof(mducStruct));
    mduc->chem = chemInit("compiled",NULL);
    return (void*)mduc;
/*   } */
}
/* -------------------------------------------------------------------------- */

void MDUCdestroy(void *mem)
{
  /* dereference internal flamelet solver structs */
  mducStruct   *mduc   = (mducStruct*) mem;
  ioStruct     *io     = mduc->io;
  driveStruct  *driver = mduc->driver;
  chemStruct   *chem   = mduc->chem;
  timeStruct   *time   = mduc->time;
  flowStruct   *flow   = mduc->flow;
  flameStruct  *flame  = mduc->flame;
  sdialStruct  *sdial  = mduc->sdial;
  engineStruct *engine = mduc->engine;
  
  destroyFlame(flame);
  destroyFlow(flow);
  destroyChem(chem);
  destroyIO(io);
  destroyEngine(engine);
  destroySDial(driver->SDIAL, sdial);
  free(sdial->args);
  
  free(io); free(driver); free(chem); free(time); free(flow); free(flame);
  free(sdial); free(engine); free(mem); 
}
/* -------------------------------------------------------------------------- */

void MDUCinitFlow(double t, double P, double phi1min, double phi1max,
                  double phi2min, double phi2max, void *mem)
{
  // dereference internal flamelet solver memory
  mducStruct  *mduc   = (mducStruct*) mem;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  timeStruct  *time   = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;
  sdialStruct *sdial  = mduc->sdial;

  time->t_m1 = time->t = t;
  flow->P_m1       = P;
  if (driver->ENGINE && driver->PRES_MODE == 3)
    flame->Y[flame->sP] = flame->Y_m1[flame->sP] = P;
  flow->phi1min_m1 = phi1min;
  flow->phi1max_m1 = phi1max;
  flow->phi2min_m1 = phi2min;
  flow->phi2max_m1 = phi2max;

  int i, j, k, r, row, col;
  for (j=0; j<flow->fhat.nW; j++) {
    flow->phiRaw[1][j].min = flow->phiHat[1][j].min = phi1min; 
    flow->phiRaw[1][j].max = flow->phiHat[1][j].max = phi1max; 
    flow->phiRaw[1][j].del = flow->phiHat[1][j].del = phi1max - phi1min; 
    flow->phiRaw[2][j].min = flow->phiHat[2][j].min = phi2min; 
    flow->phiRaw[2][j].max = flow->phiHat[2][j].max = phi2max; 
    flow->phiRaw[2][j].del = flow->phiHat[2][j].del = phi2max - phi2min;
  }

  /* compute the current properties */
  double wmix;
  for (k=0; k<flame->n2; k++) {
    row = k*flame->n1*flame->nVars;
    r   = k*flame->n1;
    for (j=0; j<flame->n1; j++) {
      col  = j*flame->nVars;
      wmix = calcMixtureW(chem->nSpecSolve, &flame->Y[row+col], chem->W);
      flame->rho[r+j] = wmix*flow->P/(flame->RGAS*flame->Y[row+col+flame->sT]);
    } /* end loop over columns */
  } /* end loop over rows */
}
/* -------------------------------------------------------------------------- */

/* Function that checks whether MDUC was compiled for MPI usage*/
int MDUCisParallel(void)
{
#ifdef ACCESS_MPI
return 1;
#else
return 0;
#endif /*ACCESS_MPI*/
}

#ifdef ACCESS_MPI
void MDUCsetPrintProc(MPI_Comm GLOBAL_COMM)
{
  MDUC_PRINT_PROC = 1;
  if (GLOBAL_COMM == MPI_COMM_NULL) return;
  int myid;
  MPI_Comm_rank(GLOBAL_COMM, &myid);
  if (myid != 0) MDUC_PRINT_PROC = 0;
  MPI_Barrier(GLOBAL_COMM);
}
/* -------------------------------------------------------------------------- */

void MDUCsetParallel(MPI_Comm GLOBAL_COMM, void *mem)
{
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  timeStruct  *time   = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;
  sdialStruct *sdial  = mduc->sdial;
  mpiStruct   *mpi    = mduc->mpi;

  /* duplicate the communicator */
  MPI_Comm_dup(GLOBAL_COMM, &mpi->world);

  /* define the communicator */
  int i,minCol = 1, nLocal, *iKeep, key;
  MPI_Comm_size(mpi->world, &mpi->nWorld);
  MPI_Comm_rank(mpi->world, &mpi->iWorld);
  MPI_Comm_group(mpi->world, &mpi->worldG);
  mpi->rxn = mpi->nWorld-1;
  nLocal = mpi->nWorld;
  if (flame->n1/mpi->nWorld < minCol) {
    /* too few columns per processor, create a local communicator with subset
       of world processors */
    nLocal = flame->n1/minCol;
    iKeep = (int*) calloc(nLocal,sizeof(int));
    for (i=0; i<nLocal; i++)
      iKeep[i] = i;

    MPI_Group_incl(mpi->worldG,nLocal,iKeep,&mpi->localG);
    MPI_Comm_create(mpi->world,mpi->localG,&mpi->comm);
    if (mpi->iWorld < nLocal) {
      MPI_Comm_size(mpi->comm, &mpi->nProc);
      MPI_Comm_rank(mpi->comm, &mpi->myid);
    }
    else {
      mpi->nProc = 0;
      mpi->myid  = -1;
    }
    free(iKeep);
  }
  else {
    /* use global communicator */
    MPI_Comm_create(mpi->world,mpi->worldG,&mpi->comm);
    MPI_Comm_size(mpi->comm, &mpi->nProc);
    MPI_Comm_rank(mpi->comm, &mpi->myid);
  }

  /* set the host to processor 0 by default */
  mpi->host = 0;

  if (mpi->myid != mpi->host)
    MDUC_PRINT_PROC = 0;

  /* if using local processor set, this processor will not solve flamelet, so 
     no need to do further initialization */
  if (mpi->myid < 0) return;

  if (flame->nDim == 1) {
    initParallel1D(mem);
  }
  else if (flame->nDim == 2) {
    if (sdial->ADI)
      initParallel2D_ADI(mem);
    else
      initParallel2D(mem);
  }
}

void MDUCsetParallelmRIF(MPI_Comm GLOBAL_COMM, void *mem, int nLocal_in, int hostLocal)
{
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  timeStruct  *time   = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;
  sdialStruct *sdial  = mduc->sdial;
  mpiStruct   *mpi    = mduc->mpi;

  /* duplicate the communicator */
  MPI_Comm_dup(GLOBAL_COMM, &mpi->world);

  /* define the communicator */
  int i,minCol = 1, *iKeep, key, nLocal;
  MPI_Comm_size(mpi->world, &mpi->nWorld);
  MPI_Comm_rank(mpi->world, &mpi->iWorld);
  MPI_Comm_group(mpi->world, &mpi->worldG);

  if (flame->n1/minCol<nLocal_in) {
	nLocal = flame->n1/minCol;
  }
  else {
   	nLocal = nLocal_in;
  }
 
  int iKeepMe = 0;
  iKeep = (int*) calloc(nLocal,sizeof(int));

  for (i=0; i<nLocal; i++) {
    iKeep[i] = (i+hostLocal) % mpi->nWorld;
//    while (iKeep[i]>=nLocal) {
//      iKeep[i] -= nLocal;
//    }
//    if (mpi->iWorld == hostLocal) {
//       fprintf(stderr, "%i - Keeping : %i\n",hostLocal,iKeep[i]);
//    }
    if (mpi->iWorld == iKeep[i]) iKeepMe = 1;

  }
//  if (iKeepMe) {
//     fprintf(stdout, "Host : %i\n",mpi->iWorld);
//  }
//  if (mpi->iWorld == hostLocal) {
//     for (i=0; i<nLocal; i++) {
//	 fprintf(stdout,"Keeping: %i\n",iKeep[i]);
//     }
//  }
  MPI_Barrier(mpi->world);
  MPI_Group_incl(mpi->worldG,nLocal,iKeep,&mpi->localG);
  MPI_Comm_create(mpi->world,mpi->localG,&mpi->comm);
  //if (mpi->iWorld < nLocal + hostLocal && mpi->iWorld >=hostLocal) {
  if (iKeepMe == 1) {
    MPI_Comm_size(mpi->comm, &mpi->nProc);
    MPI_Comm_rank(mpi->comm, &mpi->myid);
  }
  else {
    mpi->nProc = 0;
    mpi->myid  = -1;
  }
  free(iKeep);

  /* set the host to processor 0 by default */
  mpi->host = 0; // hostLocal;

  if (mpi->myid != mpi->host)
    MDUC_PRINT_PROC = 0;

  /* if using local processor set, this processor will not solve flamelet, so 
     no need to do further initialization */
  if (mpi->myid < 0) return;

  if (flame->nDim == 1) {
    initParallel1D(mem);
  }
  else if (flame->nDim == 2) {
    if (sdial->ADI)
      initParallel2D_ADI(mem);
    else
      initParallel2D(mem);
  }
}
/* -------------------------------------------------------------------------- */

void initParallel1D(void *mem)
{
  int i,j,nUA,myid;
  int nBcomm,bComm[2],bSend[2],bRcv[2],bCount[2];
  int colsPer;
  
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  flameStruct *flame  = mduc->flame;
  mpiStruct   *mpi    = mduc->mpi;
  MPI_Status status;
  myid = mpi->myid;

  /* ----- allocate memory for mpi variables ----- */
  /* column variables */
  mpi->nColLocal  = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->colStart   = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->colEnd     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->cStart     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->cCount     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->jStart     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->jCount     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->cPlus      = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->cCountPlus = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->jPlus      = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->jCountPlus = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->soot0      = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->sootN      = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->soot0xmyp  = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->sootNxmyp  = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->sootSrcP0  = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->sootSrcPN  = (int*) calloc(mpi->nProc, sizeof(int));

  /* ----- assign columns to processors ----- */
  /* determine number of columns per processor */
  colsPer = flame->n1/mpi->nProc;
  for (i=0; i<mpi->nProc; i++)
    mpi->nColLocal[i] = colsPer;
  
  /* assign any remaining columns */
  if (flame->n1%mpi->nProc != 0) {
    nUA = flame->n1 - colsPer*mpi->nProc;
    for (i=0; i<nUA; i++)
      mpi->nColLocal[mpi->nProc-i-1]++;
  }

  /* ----- set column decomposition ----- */
  /* this is done very simply at the moment with basic striping.  The first
     processor is set and and then the rest are set based on marching */
  mpi->colStart[0] = 0;
  mpi->colEnd[0]   = mpi->nColLocal[0] - 1;
  mpi->cStart[0]   = 0;
  mpi->cCount[0]   = mpi->nColLocal[0];
  mpi->jStart[0]   = 0;
  mpi->jCount[0]   = mpi->nColLocal[0]*flame->nVars;
  mpi->cPlus[0]    = 0;
  mpi->jPlus[0]    = 0;
  mpi->soot0[0]    = 0;
  mpi->sootN[0]    = mpi->nColLocal[0]*chem->nSootMom;
  mpi->soot0xmyp[0]    = 0;
  mpi->sootNxmyp[0]    = mpi->nColLocal[0]*(chem->nSootMom-1);
  mpi->sootSrcP0[0]= 0;
  mpi->sootSrcPN[0]= mpi->nColLocal[0]*chem->nSootMom*chem->nSootSrc;
  if (mpi->nColLocal[1] == 1) {
    mpi->cCountPlus[0] = 1; /*3;*/
    mpi->jCountPlus[0] = flame->nVars; /*3*flame->nVars;*/
  }
  else {
    mpi->cCountPlus[0] =  mpi->nColLocal[0] + 1;
    mpi->jCountPlus[0] = (mpi->nColLocal[0] + 1)*flame->nVars;
  }

  for (i=1; i<mpi->nProc; i++) {
    mpi->colStart[i] = mpi->colEnd[i-1]+1;
    mpi->colEnd[i]   = mpi->colStart[i]+mpi->nColLocal[i]-1;
    mpi->cStart[i]   = mpi->colStart[i];
    mpi->cCount[i]   = mpi->nColLocal[i];
    mpi->jStart[i]   = mpi->colStart[i]*flame->nVars;
    mpi->jCount[i]   = mpi->nColLocal[i]*flame->nVars;
    mpi->soot0[i]    = mpi->colStart[i]*chem->nSootMom;
    mpi->sootN[i]    = mpi->nColLocal[i]*chem->nSootMom;
    mpi->soot0xmyp[i]= mpi->colStart[i]*(chem->nSootMom-1);
    mpi->sootNxmyp[i]= mpi->nColLocal[i]*(chem->nSootMom-1);
    mpi->sootSrcP0[i]= mpi->colStart[i]*chem->nSootMom*chem->nSootSrc;
    mpi->sootSrcPN[i]= mpi->nColLocal[i]*chem->nSootMom*chem->nSootSrc;
  }

  /* ----- determine the process neighbours ----- */
  /* assume that the partitioning puts the boundaries on the host and max
     processor */
  if (myid == 0 || myid == mpi->nProc-1)
    mpi->nBcomm = 1;
  else
    mpi->nBcomm = 2;

  mpi->bComm  = (int*) calloc(mpi->nBcomm, sizeof(int));
  mpi->bCount = (int*) calloc(mpi->nBcomm, sizeof(int));
  mpi->bSend  = (int*) calloc(mpi->nBcomm, sizeof(int));
  mpi->bRcv   = (int*) calloc(mpi->nBcomm, sizeof(int));

  mpi->nSRcomm = 2*mpi->nBcomm;
  mpi->srComm  = (int*) calloc(mpi->nSRcomm, sizeof(int));
  mpi->srCount = (int*) calloc(mpi->nSRcomm, sizeof(int));
  mpi->srSender  = (int*) calloc(mpi->nSRcomm, sizeof(int));
  mpi->srRecver  = (int*) calloc(mpi->nSRcomm, sizeof(int));
  mpi->srIndex  = (int*) calloc(mpi->nSRcomm, sizeof(int));
 
  int iSR = 0; 
  
  int ii = 0;
  int nGhost = 1;
  /* set the left boundary comm */
  if (myid != 0) {
    mpi->bComm[ii]  = myid-1;
    if (nGhost > mpi->nColLocal[myid-1])
      fprintf(stderr, "Warning: ghost cells exceeds neighbouring processor: "
              "%i --> %i/%i\n",myid,nGhost,mpi->nColLocal[myid-1]);

    mpi->bCount[ii] = nGhost*flame->nVars;
    mpi->bSend[ii]  = mpi->colStart[myid]*flame->nVars;
    mpi->bRcv[ii]   = (mpi->colStart[myid]-nGhost)*flame->nVars;
    ii++;

    /* receive from left neighbour */
    mpi->srComm[iSR]    = mpi->nProc + (myid)*2-1; /* tag */
    mpi->srCount[iSR]   = nGhost*flame->nVars; /* number of elements */
    mpi->srSender[iSR]  = myid-1; /* sender */
    mpi->srRecver[iSR]  = myid;   /* receiver */
    mpi->srIndex[iSR]   = (mpi->colStart[myid]-nGhost)*flame->nVars; 
    iSR++;
 
    /* send to left neighbour */
    mpi->srComm[iSR]    = mpi->nProc + (myid-1)*2+1; /* tag */
    mpi->srCount[iSR]   = nGhost*flame->nVars; /* number of elements */
    mpi->srSender[iSR]  = myid; /* sender */
    mpi->srRecver[iSR]  = myid-1;   /* receiver */
    mpi->srIndex[iSR]   = (mpi->colStart[myid])*flame->nVars; 
    iSR++;  
  }
  /* set the right boundary */
  if (myid != mpi->nProc - 1) {
    mpi->bComm[ii]  = myid+1;
    if (nGhost > mpi->nColLocal[myid+1])
      fprintf(stderr, "Warning: ghost cells exceeds neighbouring processor: "
              "%i --> %i/%i\n",myid,nGhost,mpi->nColLocal[myid+1]);

    mpi->bCount[ii] = nGhost*flame->nVars;
    mpi->bSend[ii]  = (mpi->colEnd[myid]+1-nGhost)*flame->nVars;
    mpi->bRcv[ii]   = (mpi->colEnd[myid]+1)*flame->nVars;

    /* send to right neighbour */
    mpi->srComm[iSR]    = mpi->nProc + (myid+1)*2-1; /* tag */
    mpi->srCount[iSR]   = nGhost*flame->nVars; /* number of elements */
    mpi->srSender[iSR]  = myid; /* sender */
    mpi->srRecver[iSR]  = myid+1;   /* receiver */
    mpi->srIndex[iSR]   = (mpi->colStart[myid+1]-nGhost)*flame->nVars; 
    iSR++;
 
    /* recv from right neighbour */
    mpi->srComm[iSR]    = mpi->nProc + (myid)*2+1; /* tag */
    mpi->srCount[iSR]   = nGhost*flame->nVars; /* number of elements */
    mpi->srSender[iSR]  = myid+1; /* sender */
    mpi->srRecver[iSR]  = myid;   /* receiver */
    mpi->srIndex[iSR]   = (mpi->colStart[myid+1])*flame->nVars; 
    iSR++;  
  }

  /* make sure that the flamelet know what it has to solve */
  flame->j0 = mpi->colStart[myid];
  flame->jN = mpi->colEnd[myid];
  flame->j0m = IMAX(flame->j0 - nGhost, 0);
  flame->jNp = IMIN(flame->jN + nGhost, flame->n1-1);

  /* ----- print summary ----- */
  /* ----- write parallel info to table ----- */
  int tabWidth = 80;
  if (myid == 0) {
    for (i=0; i<tabWidth; i++) fprintf(stderr, "=");
    fprintf(stderr, "\n");
    fprintf(stderr,"   Parallel initialisation for %i processors\n",
            mpi->nProc);
    if (mpi->nProc < mpi->nWorld)
      fprintf(stderr, "   ATTENTION: subset of %i global processors\n",
              mpi->nWorld);
    fprintf(stderr,"   %i grid points, %i variables\n\n",flame->n1,
            flame->nVars);
    for (i=0; i<tabWidth; i++) fprintf(stderr, "-");
    fprintf(stderr, "\n");
    fprintf(stderr,"%5s | %18s column %10s | %16s comm %9s |\n",
            "","","","","");
    fprintf(stderr, " node | %3s %5s %5s %6s %6s %6s | ", "n", "start", "end",
            "r_j", "l_j", "c_j");
    fprintf(stderr, " %6s %6s %6s %6s    |\n", "", "send", "recv", "count");
    for (i=0; i<tabWidth; i++) fprintf(stderr, "-");
    fprintf(stderr, "\n");
    fprintf(stderr, " %4i | %3i %5i %5i %6i %6i %6i | ", 0,
            mpi->nColLocal[0], mpi->colStart[0], mpi->colEnd[0],
            mpi->jStart[0], mpi->jStart[0]+mpi->jCount[0], mpi->jCount[0]);
    fprintf(stderr, " L %3s %6s %6s %6s     |\n", "-", "-", "-", "-");
    fprintf(stderr, " %4s | %36s | ", "", "");
    fprintf(stderr, " R %3i %6i %6i %6i     |\n",mpi->bComm[0],mpi->bSend[0],
            mpi->bRcv[0], mpi->bCount[0]);
    for (i=1; i<mpi->nProc; i++) {
      nBcomm = 2;
      if (i == mpi->nProc-1)
        nBcomm = 1;
      MPI_Recv(bComm,  nBcomm, MPI_INT, i, i, mpi->comm, &status);
      MPI_Recv(bSend,  nBcomm, MPI_INT, i, i, mpi->comm, &status);
      MPI_Recv(bRcv,   nBcomm, MPI_INT, i, i, mpi->comm, &status);
      MPI_Recv(bCount, nBcomm, MPI_INT, i, i, mpi->comm, &status);
      fprintf(stderr, " %4i | %3i %5i %5i %6i %6i %6i | ", i,
              mpi->nColLocal[i], mpi->colStart[i], mpi->colEnd[i],
              mpi->jStart[i], mpi->jStart[i]+mpi->jCount[i], mpi->jCount[i]);
      fprintf(stderr, " L %3i %6i %6i %6i     |\n", bComm[0],bSend[0],bRcv[0],
              bCount[0]);
      fprintf(stderr, " %4s | %36s | ", "", "");
      if (i == mpi->nProc-1)
        fprintf(stderr, " R %3s %6s %6s %6s     |", "-", "-", "-", "-");
      else
        fprintf(stderr, " R %3i %6i %6i %6i     |",bComm[1],bSend[1],bRcv[1],
                bCount[1]);
      fprintf(stderr, "\n");
    }
    for (i=0; i<tabWidth; i++) fprintf(stderr, "-");
    fprintf(stderr, "\n");
  }
  else {
    MPI_Send(mpi->bComm,  mpi->nBcomm, MPI_INT, 0, mpi->myid, mpi->comm);
    MPI_Send(mpi->bSend,  mpi->nBcomm, MPI_INT, 0, mpi->myid, mpi->comm);
    MPI_Send(mpi->bRcv,   mpi->nBcomm, MPI_INT, 0, mpi->myid, mpi->comm);
    MPI_Send(mpi->bCount, mpi->nBcomm, MPI_INT, 0, mpi->myid, mpi->comm);
  }
}
/* -------------------------------------------------------------------------- */

void initParallel2D(void *mem)
{
  int i,j,nUA,myid;
  int nBcomm,bComm[2],bSend[2],bRcv[2],bCount[2];
  int rowsPer,rowLength,colsPer,colLength,totLength;
  
  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;
  sdialStruct *sdial  = mduc->sdial;
  mpiStruct   *mpi    = mduc->mpi;
  MPI_Status status;
  myid = mpi->myid;

  /* ----- allocate memory for mpi variables ----- */
  /* row variables */
  mpi->nRowLocal  = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->rowStart   = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->rowEnd     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->rStart     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->rCount     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->kStart     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->kCount     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->rPlus      = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->rCountPlus = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->kPlus      = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->kCountPlus = (int*) calloc(mpi->nProc, sizeof(int));

  /* column variables */
  mpi->nColLocal  = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->colStart   = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->colEnd     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->cStart     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->cCount     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->jStart     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->jCount     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->cPlus      = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->cCountPlus = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->jPlus      = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->jCountPlus = (int*) calloc(mpi->nProc, sizeof(int));

  /* ----- assign rows to processors ----- */
  /* determine number of rows per processor */
  rowsPer = flame->n2/mpi->nProc;
  for (i=0; i<mpi->nProc; i++)
    mpi->nRowLocal[i] = rowsPer;
  
  /* assign any remaining rows */
  if (flame->n2%mpi->nProc != 0) {
    nUA = flame->n2 - rowsPer*mpi->nProc;
    for (i=0; i<nUA; i++)
      mpi->nRowLocal[i]++;
  }

  /* ----- assign columns to processors ----- */
  if (sdial->ADI) {
    /* determine number of columns per processor */
    colsPer = flame->n1/mpi->nProc;
    for (i=0; i<mpi->nProc; i++)
      mpi->nColLocal[i] = colsPer;
    
    /* get the number of unallocated columns and start assigning them */
    if (flame->n1%mpi->nProc != 0) {
      nUA = flame->n1 - colsPer*mpi->nProc;
      for (i=0; i<nUA; i++)
        mpi->nColLocal[i]++;
    }
  }
  else {
    /* implicit solve uses striping */
    mpi->nColLocal[myid] = flame->n1;
  }

  totLength = flame->n1*flame->n2*flame->nVars;
  rowLength = flame->n1*flame->nVars;
  colLength = flame->n2*flame->nVars;
    
  /* ----- set row decomposition ----- */
  /* this is done very simply at the moment with basic striping.  The first
     processor is set and and then the rest are set based on marching */
  mpi->rowStart[0] = 0;
  mpi->rowEnd[0]   = mpi->nRowLocal[0] - 1;
  mpi->rStart[0]   = 0;
  mpi->rCount[0]   = mpi->nRowLocal[0]*flame->n1;
  mpi->kStart[0]   = 0;
  mpi->kCount[0]   = mpi->nRowLocal[0]*rowLength;
  mpi->rPlus[0]    = 0;
  mpi->kPlus[0]    = 0;
  if (mpi->nRowLocal[0] == 1) {
    mpi->rCountPlus[0] = 3;
    mpi->kCountPlus[0] = 3*flame->nVars;
  }
  else {
    mpi->rCountPlus[0] = (mpi->nRowLocal[0] + 1)*flame->n1;
    mpi->kCountPlus[0] = (mpi->nRowLocal[0] + 1)*rowLength;
  }

  /* set columns */
  mpi->colStart[0] = 0;
  mpi->colEnd[0]   = mpi->nColLocal[0] - 1;
  mpi->cStart[0]   = 0;
  mpi->cCount[0]   = mpi->nColLocal[0]*flame->n2;
  mpi->jStart[0]   = 0;
  mpi->jCount[0]   = mpi->nColLocal[0]*colLength;
  mpi->cPlus[0]    = 0;
  mpi->jPlus[0]    = 0;
  if (mpi->nColLocal[0] == 1) {
    mpi->cCountPlus[0] = 3*flame->n2;
    mpi->jCountPlus[0] = 3*colLength;
  }
  else {
    mpi->cCountPlus[0] = (mpi->nColLocal[0] + 1)*flame->n2;
    mpi->jCountPlus[0] = (mpi->nColLocal[0] + 1)*colLength;
  }

  for (i=1; i<mpi->nProc; i++) {
    mpi->rowStart[i] = mpi->rowEnd[i-1]+1;
    mpi->rowEnd[i]   = mpi->rowStart[i]+mpi->nRowLocal[i]-1;
    mpi->rStart[i]   = mpi->rowStart[i]*flame->n1;
    mpi->rCount[i]   = mpi->nRowLocal[i]*flame->n1;
    mpi->kStart[i]   = mpi->rowStart[i]*rowLength;
    mpi->kCount[i]   = mpi->nRowLocal[i]*rowLength;
  }
  
  /* ----- determine the process neighbours ----- */
  /* assume that the partitioning puts the boundaries on the host and max
     processor */
  if (myid == 0 || myid == mpi->nProc-1)
    mpi->nBcomm = 1;
  else
    mpi->nBcomm = 2;


  mpi->bComm  = (int*) calloc(mpi->nBcomm, sizeof(int));
  mpi->bCount = (int*) calloc(mpi->nBcomm, sizeof(int));
  mpi->bSend  = (int*) calloc(mpi->nBcomm, sizeof(int));
  mpi->bRcv   = (int*) calloc(mpi->nBcomm, sizeof(int));

  mpi->nSRcomm = 2*mpi->nBcomm;
  mpi->srComm  = (int*) calloc(mpi->nSRcomm, sizeof(int));
  mpi->srCount = (int*) calloc(mpi->nSRcomm, sizeof(int));
  mpi->srSender  = (int*) calloc(mpi->nSRcomm, sizeof(int));
  mpi->srRecver  = (int*) calloc(mpi->nSRcomm, sizeof(int));
  mpi->srIndex  = (int*) calloc(mpi->nSRcomm, sizeof(int));
 
  int iSR = 0; 
  int ii = 0;
  int nGhost = 1;
  /* set the left boundary comm */
  if (myid != 0) {
    mpi->bComm[ii]  = myid-1;
    if (nGhost > mpi->nRowLocal[myid-1])
      fprintf(stderr, "Warning: ghost cells exceed neighbouring processor\n");

    mpi->bCount[ii] = nGhost*rowLength;
    mpi->bSend[ii]  = mpi->rowStart[myid]*rowLength;
    mpi->bRcv[ii]   = (mpi->rowStart[myid]-nGhost)*rowLength;

    ii++;

    /* receive from left neighbour */
    mpi->srComm[iSR]    = mpi->nProc + (myid)*2-1; /* tag */
    mpi->srCount[iSR]   = nGhost*rowLength; /* number of elements */
    mpi->srSender[iSR]  = myid-1; /* sender */
    mpi->srRecver[iSR]  = myid;   /* receiver */
    mpi->srIndex[iSR]   = (mpi->rowStart[myid]-nGhost)*rowLength; 
    iSR++;
 
    /* send to left neighbour */
    mpi->srComm[iSR]    = mpi->nProc + (myid-1)*2+1; /* tag */
    mpi->srCount[iSR]   = nGhost*rowLength; /* number of elements */
    mpi->srSender[iSR]  = myid; /* sender */
    mpi->srRecver[iSR]  = myid-1;   /* receiver */
    mpi->srIndex[iSR]   = (mpi->rowStart[myid])*rowLength; 
    iSR++;  
  }
  /* set the right boundary */
  if (myid != mpi->nProc - 1) {
    mpi->bComm[ii]  = myid+1;
    if (nGhost > mpi->nRowLocal[myid+1])
      fprintf(stderr, "Warning: ghost cells exceeds neighbouring processor\n");

    mpi->bCount[ii] = nGhost*rowLength;
    mpi->bSend[ii]  = (mpi->rowEnd[myid]+1-nGhost)*rowLength;
    mpi->bRcv[ii]   = (mpi->rowEnd[myid]+1)*rowLength;

    /* send to right neighbour */
    mpi->srComm[iSR]    = mpi->nProc + (myid+1)*2-1; /* tag */
    mpi->srCount[iSR]   = nGhost*rowLength; /* number of elements */
    mpi->srSender[iSR]  = myid; /* sender */
    mpi->srRecver[iSR]  = myid+1;   /* receiver */
    mpi->srIndex[iSR]   = (mpi->rowStart[myid+1]-nGhost)*rowLength; 
    iSR++;
 
    /* recv from right neighbour */
    mpi->srComm[iSR]    = mpi->nProc + (myid)*2+1; /* tag */
    mpi->srCount[iSR]   = nGhost*rowLength; /* number of elements */
    mpi->srSender[iSR]  = myid+1; /* sender */
    mpi->srRecver[iSR]  = myid;   /* receiver */
    mpi->srIndex[iSR]   = (mpi->rowStart[myid+1])*rowLength; 
    iSR++;  
  }

  
  
  /* make sure that the flamelet know what it has to solve */
  flame->j0  = 0;
  flame->jN  = flame->n1-1;
  flame->j0m = IMAX(flame->j0 - nGhost, 0);
  flame->jNp = IMIN(flame->jN + nGhost, flame->n1-1);
  flame->k0  = mpi->rowStart[myid];
  flame->kN  = mpi->rowEnd[myid];
  flame->k0m = IMAX(flame->k0 - nGhost, 0);
  flame->kNp = IMIN(flame->kN + nGhost, flame->n2-1);

  /* ----- print summary ----- */
  /* ----- write parallel info to table ----- */
  int tabWidth = 80;
  if (myid == 0) {
    for (i=0; i<tabWidth; i++) fprintf(stderr, "=");
    fprintf(stderr, "\n");
    fprintf(stderr,"   Parallel initialisation for %i processors\n",
            mpi->nProc);
    fprintf(stderr,"   %i grid points, %i variables\n\n",flame->n1,
            flame->nVars);
    for (i=0; i<tabWidth; i++) fprintf(stderr, "-");
    fprintf(stderr, "\n");
    fprintf(stderr,"%5s | %20s row %12s | %16s comm %9s |\n",
            "","","","","");
    fprintf(stderr, " node | %3s %5s %5s %6s %6s %6s | ", "n", "start", "end",
            "r_j", "l_j", "c_j");
    fprintf(stderr, " %6s %6s %6s %6s    |\n", "", "send", "recv", "count");
    for (i=0; i<tabWidth; i++) fprintf(stderr, "-");
    fprintf(stderr, "\n");
    fprintf(stderr, " %4i | %3i %5i %5i %6i %6i %6i | ", 0,
            mpi->nRowLocal[0], mpi->rowStart[0], mpi->rowEnd[0],
            mpi->kStart[0], mpi->kStart[0]+mpi->kCount[0], mpi->kCount[0]);
    fprintf(stderr, " L %3s %6s %6s %6s     |\n", "-", "-", "-", "-");
    fprintf(stderr, " %4s | %36s | ", "", "");
    fprintf(stderr, " R %3i %6i %6i %6i     |\n",mpi->bComm[0],mpi->bSend[0],
            mpi->bRcv[0], mpi->bCount[0]);
    for (i=1; i<mpi->nProc; i++) {
      nBcomm = 2;
      if (i == mpi->nProc-1)
        nBcomm = 1;
      MPI_Recv(bComm,  nBcomm, MPI_INT, i, i, mpi->comm, &status);
      MPI_Recv(bSend,  nBcomm, MPI_INT, i, i, mpi->comm, &status);
      MPI_Recv(bRcv,   nBcomm, MPI_INT, i, i, mpi->comm, &status);
      MPI_Recv(bCount, nBcomm, MPI_INT, i, i, mpi->comm, &status);
      fprintf(stderr, " %4i | %3i %5i %5i %6i %6i %6i | ", i,
              mpi->nRowLocal[i], mpi->rowStart[i], mpi->rowEnd[i],
              mpi->kStart[i], mpi->kStart[i]+mpi->kCount[i], mpi->kCount[i]);
      fprintf(stderr, " L %3i %6i %6i %6i     |\n", bComm[0],bSend[0],bRcv[0],
              bCount[0]);
      fprintf(stderr, " %4s | %36s | ", "", "");
      if (i == mpi->nProc-1)
        fprintf(stderr, " R %3s %6s %6s %6s     |", "-", "-", "-", "-");
      else
        fprintf(stderr, " R %3i %6i %6i %6i     |",bComm[1],bSend[1],bRcv[1],
                bCount[1]);
      fprintf(stderr, "\n");
    }
    for (i=0; i<tabWidth; i++) fprintf(stderr, "-");
    fprintf(stderr, "\n");
  }
  else {
    MPI_Send(mpi->bComm,  mpi->nBcomm, MPI_INT, 0, mpi->myid, mpi->comm);
    MPI_Send(mpi->bSend,  mpi->nBcomm, MPI_INT, 0, mpi->myid, mpi->comm);
    MPI_Send(mpi->bRcv,   mpi->nBcomm, MPI_INT, 0, mpi->myid, mpi->comm);
    MPI_Send(mpi->bCount, mpi->nBcomm, MPI_INT, 0, mpi->myid, mpi->comm);
  }
}
/* -------------------------------------------------------------------------- */

void initParallel2D_ADI(void *mem)
{
  int i, nUA;
  int rowsPer, colsPer;
  int rowLength, colLength, totLength;

  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;
  mpiStruct   *mpi    = mduc->mpi;

  /* ----- allocate memory for mpi variables ----- */
  /* row variables */
  mpi->nRowLocal  = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->rowStart   = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->rowEnd     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->rStart     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->rCount     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->kStart     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->kCount     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->rPlus      = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->rCountPlus = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->kPlus      = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->kCountPlus = (int*) calloc(mpi->nProc, sizeof(int));

  /* column variables */
  mpi->nColLocal  = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->colStart   = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->colEnd     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->cStart     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->cCount     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->jStart     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->jCount     = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->cPlus      = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->cCountPlus = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->jPlus      = (int*) calloc(mpi->nProc, sizeof(int));
  mpi->jCountPlus = (int*) calloc(mpi->nProc, sizeof(int));

  /* define the number of variables in row/column for solution array */
  totLength = flame->n1*flame->n2*flame->nVars;
  rowLength = flame->n1*flame->nVars;
  colLength = flame->n2*flame->nVars;

  /* return proper values if number of processors not specified */
  /* !!! This is redundant, function shouldn't be called in serial */
  if (mpi->nProc == 1) {
    /* set row variables */
    mpi->nRowLocal[0]  = flame->n2;
    mpi->rowStart[0]   = 0;
    mpi->rowEnd[0]     = flame->n2-1;
    mpi->rStart[0]     = 0;
    mpi->rCount[0]     = flame->n1*flame->n2;
    mpi->kStart[0]     = 0;
    mpi->kCount[0]     = flame->n1*flame->n2*flame->nVars;
    mpi->rPlus[0]      = 0;
    mpi->rCountPlus[0] = flame->n1*flame->n2;
    mpi->kPlus[0]      = 0;
    mpi->kCountPlus[0] = flame->n1*flame->n2*flame->nVars;
    
    /* set column variables */
    mpi->nColLocal[0]  = flame->n1;
    mpi->colStart[0]   = 0;
    mpi->colEnd[0]     = flame->n1-1;
    mpi->cStart[0]     = 0;
    mpi->cCount[0]     = flame->n1*flame->n2;
    mpi->jStart[0]     = 0;
    mpi->jCount[0]     = flame->n1*flame->n2*flame->nVars;
    mpi->cPlus[0]      = 0;
    mpi->cCountPlus[0] = flame->n1*flame->n2;
    mpi->jPlus[0]      = 0;
    mpi->jCountPlus[0] = flame->n1*flame->n2*flame->nVars;
    
    return;
  }
  
  /* ----- assign rows to processors ----- */
  /* determine number of rows per processor */
  rowsPer = flame->n2/mpi->nProc;
  for (i=0; i<mpi->nProc; i++)
    mpi->nRowLocal[i] = rowsPer;

  /* get the number of unallocated rows and start assigning them */
  if (flame->n2%mpi->nProc != 0) {
    nUA = flame->n2 - rowsPer*mpi->nProc;
    for (i=0; i<nUA; i++)
      mpi->nRowLocal[i]++;
  }

  /* ----- assign columns to processors ----- */
  /* determine number of columns per processor */
  colsPer = flame->n1/mpi->nProc;
  for (i=0; i<mpi->nProc; i++)
    mpi->nColLocal[i] = colsPer;
  
  /* get the number of unallocated columns and start assigning them */
  if (flame->n1%mpi->nProc != 0) {
    nUA = flame->n1 - colsPer*mpi->nProc;

    /* for 2D mixture fraction (my transformation), solve one less column
       if (flame->nDim == 2 && driver->H_DIM == 0)
       nUA--; */

    for (i=0; i<nUA; i++)
      mpi->nColLocal[i]++;
  }
  else {
    /* for 2D mixture fraction (my transformation), solve one less column
       if (flame->nDim == 2 && driver->H_DIM == 0)
       mpi->nColLocal[mpi->nProc-1]--; */
  }

  /* ----- assign starting and end rows/cols to processors ----- */
  
  /* ----- set host rows/columns ----- */
  mpi->rowStart[0] = mpi->nRowLocal[1];
  mpi->rowEnd[0]   = mpi->rowStart[0] + mpi->nRowLocal[0]-1;
  mpi->colStart[0] = mpi->nColLocal[1];
  mpi->colEnd[0]   = mpi->colStart[0] + mpi->nColLocal[0]-1;
  
  /* NOTE: with multiple processors, the indices/length of host are not needed
     since the data will not be sent anywhere */
  
  /* ----- set 1st node to boundary (save communication) ----- */
  /* set rows */
  mpi->rowStart[1] = 0;
  mpi->rowEnd[1]   = mpi->nRowLocal[1] - 1;
  mpi->rStart[1]   = 0;
  mpi->rCount[1]   = mpi->nRowLocal[1]*flame->n1;
  mpi->kStart[1]   = 0;
  mpi->kCount[1]   = mpi->nRowLocal[1]*rowLength;
  mpi->rPlus[1]    = 0;
  mpi->kPlus[1]    = 0;
  if (mpi->nRowLocal[1] == 1) {
    mpi->rCountPlus[1] = 3*flame->n1;
    mpi->kCountPlus[1] = 3*rowLength;
  }
  else {
    mpi->rCountPlus[1] = (mpi->nRowLocal[1] + 1)*flame->n1;
    mpi->kCountPlus[1] = (mpi->nRowLocal[1] + 1)*rowLength;
  }

  /* set columns */
  mpi->colStart[1] = 0;
  mpi->colEnd[1]   = mpi->nColLocal[1] - 1;
  mpi->cStart[1]   = 0;
  mpi->cCount[1]   = mpi->nColLocal[1]*flame->n2;
  mpi->jStart[1]   = 0;
  mpi->jCount[1]   = mpi->nColLocal[1]*colLength;
  mpi->cPlus[1]    = 0;
  mpi->jPlus[1]    = 0;
  if (mpi->nColLocal[1] == 1) {
    mpi->cCountPlus[1] = 3*flame->n2;
    mpi->jCountPlus[1] = 3*colLength;
  }
  else {
    mpi->cCountPlus[1] = (mpi->nColLocal[1] + 1)*flame->n2;
    mpi->jCountPlus[1] = (mpi->nColLocal[1] + 1)*colLength;
  }
  
  /* assign the rest of the rows/columns to respective processors */
  for (i=2; i<mpi->nProc; i++) {
    /* set start rows/columns */
    if (i == 2) {
      mpi->rowStart[i] = mpi->rowEnd[0] + 1;
      mpi->colStart[i] = mpi->colEnd[0] + 1;
    }
    else {
      mpi->rowStart[i] = mpi->rowEnd[i-1] + 1;
      mpi->colStart[i] = mpi->colEnd[i-1] + 1;
    }

    /* set end rows/columns */
    mpi->rowEnd[i] = mpi->rowStart[i] + mpi->nRowLocal[i] - 1;
    mpi->colEnd[i] = mpi->colStart[i] + mpi->nColLocal[i] - 1;
    
    /* set row indices/lengths */
    mpi->rStart[i] = mpi->rowStart[i]*flame->n1;
    mpi->rCount[i] = mpi->nRowLocal[i]*flame->n1;
    mpi->kStart[i] = mpi->rowStart[i]*rowLength;
    mpi->kCount[i] = mpi->nRowLocal[i]*rowLength;
    
    if (mpi->rowStart[i] == flame->n2 - 1) {
      if (mpi->nRowLocal[i] > 1) {
        fprintf(stderr, "Error: rowStart at boundary with n > 1\n");
        exit(2);
      }
      mpi->rPlus[i] = (mpi->rowStart[i]-2)*flame->n1;
      mpi->kPlus[i] = (mpi->rowStart[i]-2)*rowLength;
    }
    else {
      mpi->rPlus[i] = (mpi->rowStart[i]-1)*flame->n1;
      mpi->kPlus[i] = (mpi->rowStart[i]-1)*rowLength;
    }

    if (mpi->rowEnd[i] == flame->n2 - 1) {
      if (mpi->nRowLocal[i] == 1) {
        mpi->rCountPlus[i] = (mpi->nRowLocal[i]+2)*flame->n1;
        mpi->kCountPlus[i] = (mpi->nRowLocal[i]+2)*rowLength;
      }
      else {
        mpi->rCountPlus[i] = (mpi->nRowLocal[i]+1)*flame->n1;
        mpi->kCountPlus[i] = (mpi->nRowLocal[i]+1)*rowLength;
      }
    }
    else {
      mpi->rCountPlus[i] = (mpi->nRowLocal[i]+2)*flame->n1;
      mpi->kCountPlus[i] = (mpi->nRowLocal[i]+2)*rowLength;
    }

    /* set column indices/lengths */
    mpi->cStart[i] = mpi->colStart[i]*flame->n2;
    mpi->cCount[i] = mpi->nColLocal[i]*flame->n2;
    mpi->jStart[i] = mpi->colStart[i]*colLength;
    mpi->jCount[i] = mpi->nColLocal[i]*colLength;
    
    if (mpi->colStart[i] == flame->n1 - 1) {
      if (mpi->nColLocal[i] > 1) {
        fprintf(stderr, "Error: colStart at boundary with n > 1\n");
        exit(2);
      }
      mpi->cPlus[i] = (mpi->colStart[i]-2)*flame->n2;
      mpi->jPlus[i] = (mpi->colStart[i]-2)*colLength;
    }
    else {
      mpi->cPlus[i] = (mpi->colStart[i]-1)*flame->n2;
      mpi->jPlus[i] = (mpi->colStart[i]-1)*colLength;
    }
    
    if (mpi->colEnd[i] == flame->n1 - 1) {
      if (mpi->nColLocal[i] == 1) {
        mpi->cCountPlus[i] = (mpi->nColLocal[i]+2)*flame->n2;
        mpi->jCountPlus[i] = (mpi->nColLocal[i]+2)*colLength;
      }
      else {
        mpi->cCountPlus[i] = (mpi->nColLocal[i]+1)*flame->n2;
        mpi->jCountPlus[i] = (mpi->nColLocal[i]+1)*colLength;
      }
    }
    else {
      mpi->cCountPlus[i] = (mpi->nColLocal[i]+2)*flame->n2;
      mpi->jCountPlus[i] = (mpi->nColLocal[i]+2)*colLength;
    }
  }

  /* ----- write parallel info to table ----- */
  int tabWidth = 99;
  if (mpi->myid == 0) {
    for (i=0; i<tabWidth; i++) fprintf(stderr, "=");
    fprintf(stderr, "\n");
    fprintf(stderr,"   Parallel initialisation for %i processors\n",
            mpi->nProc);
    fprintf(stderr,"   %i rows, %i columns, %i variables\n\n",
            flame->n2, flame->n1, flame->nVars);
    for (i=0; i<tabWidth; i++) fprintf(stderr, "-");
    fprintf(stderr, "\n");
    fprintf(stderr,"%5s | %19s row %19s | %18s column %17s |\n",
            "", "", "", "", "");
    fprintf(stderr, " node | %3s %5s %5s %6s %6s %6s %6s |", "n", "start",
            "end", "r_i", "l_i", "r_fi", "l_fi");
    fprintf(stderr, " %3s %5s %5s %6s %6s %6s %6s |\n", "n", "start", "end",
            "r_i", "l_i", "r_fi", "l_fi");
    for (i=0; i<tabWidth; i++) fprintf(stderr, "-");
    fprintf(stderr, "\n");
    for (i=0; i<mpi->nProc; i++) {
      fprintf(stderr, "%5i | %3i %5i %5i %6i %6i %6i %6i |", i,
              mpi->nRowLocal[i], mpi->rowStart[i], mpi->rowEnd[i],
              mpi->rStart[i],mpi->rCount[i],mpi->kStart[i],mpi->kCount[i]);
      fprintf(stderr, " %3i %5i %5i %6i %6i %6i %6i |\n", mpi->nColLocal[i],
              mpi->colStart[i], mpi->colEnd[i], mpi->cStart[i],
              mpi->cCount[i], mpi->jStart[i], mpi->jCount[i]);
    }
    for (i=0; i<tabWidth; i++) fprintf(stderr, "-");
    fprintf(stderr, "\n");
  }
}
/* -------------------------------------------------------------------------- */
#else
void MDUCsetPrintProc(void *dum)
{
  MDUC_PRINT_PROC = 1;
}
/* -------------------------------------------------------------------------- */
#endif /* ACCESS_MPI */

void MDUCinit2Dfrom_BC_and_1D(int nGrid2,int nFuel2,char **fuelName,double *Y,
                              double Tfuel2,void *mem)
{
  int i, j, k, r, nGrid;
  char **useFuelName;
  double *useY;

  // dereference the solution of the 1D flamelet
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;
  timeStruct  *time   = mduc->time;

  //broadcast data from previous time step to all processors
#ifdef ACCESS_MPI
  mpiStruct *mpi = mduc->mpi; // dereference mpi info
  MPI_Status stat;  // for receive statements

  MPI_Bcast(&time->t_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&time->t, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&time->dt_pre, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&time->dt_cvode_pre, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&time->nt, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
 
  MPI_Bcast(&flow->P_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow->P_pre, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow->phi1min_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow->phi1max_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&(flow->chi1_m1[0]),flame->n1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow->phi2min_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow->phi2max_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);

  MPI_Bcast(&(flame->Y_m1[0]), flame->n1*flame->nVars, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&(flame->Y_m2[0]), flame->n1*flame->nVars, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  if (driver->SDIAL == 1)
    MPI_Bcast(&(flame->Yp_m1[0]), flame->n1*flame->nVars, MPI_DOUBLE_PRECISION, 0, mpi->comm);
 
#endif

  if (nGrid2 == -1) /* use the same b.c.s as Z1 */
    {
      nGrid2      = flame->n1;
      nFuel2      = io->nFuel1;
      Tfuel2      = io->fuelT1;
      useFuelName = io->fuelName1;
      useY        = io->fuelY1;
      io->nPivot2 = io->nPivot1;
      io->gPivot2 = io->gPivot1;
      io->gPerc2  = io->gPerc1;
    }

  // init the io stuff
  io->nFuel2     = nFuel2;
  io->fuelT2     = Tfuel2;
  io->bcType[2]  = 'T';
  io->fuelY2     = (double*) calloc(io->nFuel2, sizeof(double));
  io->fuelIndex2 = (int*)    calloc(io->nFuel2, sizeof(int));
  io->fuelName2  = (char**)  calloc(io->nFuel2, sizeof(char*));
  for (i=0; i<io->nFuel2; i++)
    {
      io->fuelName2[i] = (char*) calloc(25, sizeof(char));

      // assign all fuel mass fractions and names
      io->fuelY2[i] = useY[i];
      strcpy(io->fuelName2[i], useFuelName[i]);
      io->fuelIndex2[i] = getSpeciesIndex(chem->nSpecies, chem->species,
                                          io->fuelName2[i]);
    }
  
  // initialise the second grid dimension
  flame->nDim = 2;
  flame->n2   = nGrid2;
  initGrid(flame->n2,0,1,io->nPivot2,io->gPivot2,io->gPerc2,
           &flame->x2,&flame->dx2);
  
  /* ----- initialise the 2D field ----- */
  initX2dimension(driver, io, chem, flame, flow, NULL);

  /* ----- reinitialise the metric ----- */
  for (i=0; i<flame->n1; i++)
    {
      free(flame->grad.x[i]); free(flame->grad.xx[i]);
    }
  free(flame->grad.x); free(flame->grad.xx);
  initMetric(flame);
  
  /* ----- reallocate storage for scalar dissipation ----- */
  nGrid = flame->n1*flame->n2;
  flow->chi1     = (double*) realloc(flow->chi1,     nGrid*(sizeof(double)));
  flow->chi12    = (double*) realloc(flow->chi12,    nGrid*(sizeof(double)));
  flow->chi2     = (double*) realloc(flow->chi2,     nGrid*(sizeof(double)));
  flow->chi1_m1  = (double*) realloc(flow->chi1_m1,  nGrid*(sizeof(double)));
  flow->chi12_m1 = (double*) realloc(flow->chi12_m1, nGrid*(sizeof(double)));
  flow->chi2_m1  = (double*) realloc(flow->chi2_m1,  nGrid*(sizeof(double)));

  // set all the newly allocated memory for chi1 to zero
  for (k=1; k<flame->n2; k++)
    {
      r = k*flame->n1;
      for (j=0; j < flame->n1; j++)
        {
          flow->chi1[r+j]    = flow->chi2[r+j]    = flow->chi12[r+j] = 0.0;
          flow->chi1_m1[r+j] = flow->chi2_m1[r+j] = flow->chi12_m1[r+j] = 0.0;
        }
    }
}
/* -------------------------------------------------------------------------- */

void MDUCinit2Dfrom_1D_and_1D(void *mem2, void *mem1)
{
  int i, j, k, r, nGrid;

  // dereference the solution of the 1D flamelet
  mducStruct  *mduc   = (mducStruct*) mem1;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  timeStruct  *time   = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;

  // dereference the solution of the second 1D flamelet
  mducStruct  *mduc2   = (mducStruct*) mem2;
  ioStruct    *io2     = mduc2->io;
  driveStruct *driver2 = mduc2->driver;
  chemStruct  *chem2   = mduc2->chem;
  timeStruct  *time2   = mduc2->time;
  flowStruct  *flow2   = mduc2->flow;
  flameStruct *flame2  = mduc2->flame;

#ifdef ACCESS_MPI

  mpiStruct *mpi = mduc->mpi; // dereference mpi info
  mpiStruct *mpi2 = mduc2->mpi; // dereference mpi info
  MPI_Status stat;  // for receive statements

  /* -----  Broadcast the time from 1D solution to all processors ----- */
  MPI_Bcast(&time->t_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&time->t, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&time2->t_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi2->comm);
  MPI_Bcast(&time2->t, 1, MPI_DOUBLE_PRECISION, 0, mpi2->comm);
  MPI_Bcast(&time->dt_pre, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&time2->dt_pre, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&time->dt_cvode_pre, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&time2->dt_cvode_pre, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&time->nt, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&time2->nt, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
 
  MPI_Bcast(&flow->P_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow2->P_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow->P_pre, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow2->P_pre, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow->phi1min_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow2->phi1min_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow->phi1max_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow2->phi1max_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&(flow->chi1_m1[0]),flame->n1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&(flow2->chi1_m1[0]),flame2->n1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow->phi2min_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow2->phi2min_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow->phi2max_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&flow2->phi2max_m1, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);

  MPI_Bcast(&(flame->Y_m1[0]), flame->n1*flame->nVars, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&(flame2->Y_m1[0]), flame2->n1*flame2->nVars, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&(flame->Y_m2[0]), flame->n1*flame->nVars, MPI_DOUBLE_PRECISION, 0, mpi->comm);
  MPI_Bcast(&(flame2->Y_m2[0]), flame2->n1*flame2->nVars, MPI_DOUBLE_PRECISION, 0, mpi->comm);

  if (driver->SDIAL == 1)
    {
      MPI_Bcast(&(flame->Yp_m1[0]), flame->n1*flame->nVars, MPI_DOUBLE_PRECISION, 0, mpi->comm);
      MPI_Bcast(&(flame2->Yp_m1[0]), flame2->n1*flame2->nVars, MPI_DOUBLE_PRECISION, 0, mpi->comm);
    }

#endif

  // init the io stuff
  io->nFuel2 = io2->nFuel1;
  io->fuelT2 = io2->fuelT1;
  io->bcType[2] = 'E';
  io->fuelY2     = (double*) calloc(io->nFuel2, sizeof(double));
  io->fuelIndex2 = (int*)    calloc(io->nFuel2, sizeof(int));
  io->fuelName2  = (char**)  calloc(io->nFuel2, sizeof(char*));
  for (i=0; i<io->nFuel2; i++)
    {
      // assign all fuel mass fractions and names
      io->fuelY2[i] = io2->fuelY1[i];
      io->fuelName2[i] = (char*) calloc(25, sizeof(char));
      strcpy(io->fuelName2[i], io2->fuelName1[i]);
      io->fuelIndex2[i] = getSpeciesIndex(chem->nSpecies, chem->species, 
                                          io->fuelName2[i]);
    }
  
  // initialise the second grid dimension from the other 1D flamelet
  flame->nDim = 2;
  flame->n2   = flame2->n1;
  flame->x2   = (double*) calloc(flame->n2, sizeof(double));
  flame->dx2  = (double*) calloc(flame->n2-1, sizeof(double));
  for (k=0; k<flame->n2; k++)   flame->x2[k]  = flame2->x1[k];
  for (k=0; k<flame->n2-1; k++) flame->dx2[k] = flame2->dx1[k];

  /* ----- initialise the 2D field ----- */
  initX2dimension(driver, io, chem, flame, flow, flame2->Y);

  /* ----- reinitialise the metric ----- */
  for (i=0; i<flame->n1; i++)
    {
      free(flame->grad.x[i]); free(flame->grad.xx[i]);
    }
  free(flame->grad.x); free(flame->grad.xx);
  initMetric(flame);

  /* ----- reallocate storage for scalar dissipation ----- */
  nGrid = flame->n1*flame->n2;
  flow->chi1     = (double*) realloc(flow->chi1,     nGrid*(sizeof(double)));
  flow->chi12    = (double*) realloc(flow->chi12,    nGrid*(sizeof(double)));
  flow->chi2     = (double*) realloc(flow->chi2,     nGrid*(sizeof(double)));
  flow->chi1_m1  = (double*) realloc(flow->chi1_m1,  nGrid*(sizeof(double)));
  flow->chi12_m1 = (double*) realloc(flow->chi12_m1, nGrid*(sizeof(double)));
  flow->chi2_m1  = (double*) realloc(flow->chi2_m1,  nGrid*(sizeof(double)));

  // set all the newly allocated memory for chi1 to zero
  for (k=1; k<flame->n2; k++)
    {
      r = k*flame->n1;
      for (j=0; j < flame->n1; j++)
        {
          flow->chi1[r+j]    = flow->chi2[r+j]    = flow->chi12[r+j] = 0.0;
          flow->chi1_m1[r+j] = flow->chi2_m1[r+j] = flow->chi12_m1[r+j] = 0.0;
        }
    }
}
/* -------------------------------------------------------------------------- */

void MDUCupdateTime(double t, double dt, int timeStep, void *mem)
{
  timeStruct  *time   = (timeStruct*) ((mducStruct*)mem)->time;
  time->t  = t;
  time->dt = dt;
  time->nt = timeStep;
}

void MDUCupdateTimestep(int timeStep, void *mem)
{
  timeStruct  *time   = (timeStruct*) ((mducStruct*)mem)->time;
  time->nt = timeStep;
}
/* -------------------------------------------------------------------------- */

void MDUCsetAvgSdot(int m, int n, double *Sdot, void *mem)
{
  int i;
  flowStruct  *flow   = (flowStruct*) ((mducStruct*)mem)->flow;
  for (i=0; i<n; i++) {
    if (m == 0)
      flow->avgSdot[i] = Sdot[i];
    else if (m == -1)
      flow->avgSdot_m1[i] = Sdot[i];
  }
}
/* -------------------------------------------------------------------------- */

void MDUCsetAvgH(int n, double *avg_H, void *mem)
{
  int i,j, k ;
  mducStruct  *mduc   = (mducStruct*) mem;
  flowStruct  *flow   = mduc->flow;
  driveStruct *driver = mduc->driver;
  flameStruct *flame  = mduc->flame;
  chemStruct  *chem   = mduc->chem;
  double **initH, *mod_H, *avg_Hl,*Y_j;
  double *x1   = flame->x1;
  double *x2   = flame->x2;
  int     nVar = flame->nVars;
  int     n1   = flame->n1;
  int     n2   = flame->n2;
  int  iX1, row, row_pre, col;

  Y_j    = (double*) calloc(chem->nSpecSolve, sizeof(double));
  initH  = (double**) calloc(n2, sizeof(double*));
  for (i=0; i<n2; ++i) 
    initH[i] = (double*) calloc(n1, sizeof(double));

  if (driver->H_DIM==2)
    {
      avg_Hl   = (double*) calloc(n2, sizeof(double));
      mod_H    = (double*) calloc(n2, sizeof(double));
    }
  else if (driver->H_DIM==1)
    {
      avg_Hl   = (double*) calloc(n1, sizeof(double));
      mod_H    = (double*) calloc(n1, sizeof(double));
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

  // Copy to avg_Hl
  for (j=0; j<n1; j++)
    avg_Hl[j] = avg_H[j];

  // Go through and fill in avgH if it is not present in flow domain
  for (j=0; j<n1; j++)
    {
      if (avg_Hl[j]==0.0)
	{
	  // Get mod_H for middle row
	  k = (n2)/2;
	  mod_H[k] = linInterp(x2[k],x2[0],mod_H[0],x2[n2-1],mod_H[n2-1]);

	  // Get linear interpolated enthalpy
	  initH[k][0] = mod_H[k] + avg_Hl[0];
	  initH[k][n1-1] = mod_H[k] + avg_Hl[n1-1];
	  initH[k][j] = linInterp(x1[j], x1[0], initH[k][0],
				  x1[n1-1], initH[k][n1-1]);
	  avg_Hl[j] = initH[k][j]-mod_H[k];
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
      initH[k][n1-1] = mod_H[k] + avg_Hl[n1-1];
      col = (n1-1)*nVar;
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
	  initH[k][j] = mod_H[k] + avg_Hl[j];
	      
	  if (driver->ENERGY)
	    flame->Y[row+col+flame->sT] = initH[k][j];
	  else // get the temperature based on the enthalpy
	    flame->Y[row+col+flame->sT] =
	      calcTfromH(chem->nSpecies,chem->nSpecSolve,initH[k][j],
                   &flame->Y[row+col],flame->Y[row_pre + flame->sT], chem);
	} // end loop over columns
    } // end loop over rows


  for (i=0; i<n2; ++i) free(initH[i]);
  free(initH);
  free(avg_Hl);
  free(Y_j);
}
/* -------------------------------------------------------------------------- */

void MDUCsetScalarDiss(int n1, int n2, double **chi1, double **chi12,
                       double **chi2, void *mem)
{
  int j, k, flag;
  driveStruct *driver = (driveStruct*) ((mducStruct*)mem)->driver;
  flowStruct  *flow   = (flowStruct*)  ((mducStruct*)mem)->flow;
  flameStruct *flame  = (flameStruct*) ((mducStruct*)mem)->flame;
  sdialStruct *sdial  = (sdialStruct*) ((mducStruct*)mem)->sdial;
  timeStruct  *time   = (timeStruct*)  ((mducStruct*)mem)->time;

  for (k=0; k<n2; k++)
    {
      for (j=0; j<n1; j++)
        {
          if (chi1 != NULL) flow->chi1_m1[k*n1+j] = chi1[k][j];
          if (chi2 != NULL) flow->chi2_m1[k*n1+j] = chi2[k][j];
          if (chi12 != NULL) flow->chi12_m1[k*n1+j]  = chi12[k][j];

          // !!! generalize for 2D later
          if (driver->CHI_MODE == 2)
            {
              // for now, specify a constant strain rate profile for post-proc.
              flow->strainCondZ[j] = flow->strain[1];

              // reinitialise cvode solver
              sdial->Ydata[j*flame->nVars + flame->sChi] = chi1[k][j];
#ifdef SUNDIALS_2_3
              if (driver->SDIAL == 0)
                flag = CVodeReInit(sdial->mem, cvodeRHS,time->t,sdial->Y,CV_SS,
                                   driver->RELTOL, &driver->ABSTOL);
              else if (driver->SDIAL == 1)
                flag = IDAReInit(sdial->mem, idaRHS,time->t,sdial->Y,sdial->Yp,
                                 IDA_SS, driver->RELTOL, &driver->ABSTOL);
#else
              if (driver->SDIAL == 0)
                flag = CVodeReInit(sdial->mem, time->t, sdial->Y);
              else if (driver->SDIAL == 1)
                flag = IDAReInit(sdial->mem, time->t, sdial->Y, sdial->Yp);
#endif /* SUNDIALS_2_3 */
            }
        } // end loop over columns
    } // end loop over rows
}
/* -------------------------------------------------------------------------- */

void MDUCsetPDF(int n1, int n2, double **pdf, void *mem)
{
  int j, k, flag;
  flowStruct  *flow   = (flowStruct*)  ((mducStruct*)mem)->flow;
  flameStruct *flame  = (flameStruct*) ((mducStruct*)mem)->flame;

  if (flow->pdf == NULL || n1 != flame->n1 || n2 != flame->n2) {
    fprintf(stderr, "Warning: pdf memory not correct, returning\n");
    return;
  }

  for (k=0; k<n2; k++) {
    for (j=0; j<n1; j++) {
      flow->pdf[k*n1+j] = pdf[k][j];
    }
  }
}
/* -------------------------------------------------------------------------- */

void MDUCsetSpeciesIndex(int nY, char **species, int **index, void *mem)
{
  int i;
  // error checking
  if (mem == NULL)
    mducErr("species list must be defined in flamelet memory",
            "MDUCsetSpeciesIndex");

  // allocate the memory for storing indices
  index[0] = (int*) calloc(nY, sizeof(int));

  // get information from flamelet memory
  mducStruct  *mduc   = (mducStruct*) mem;

  // assign species indices
  for (i=0; i<nY; i++)
    {
      if (!strcmp(species[i], "dummy")) // if dummy species, set negative
        index[0][i] = -1;
      else
        index[0][i] = getSpeciesIndex(mduc->chem->nSpecies,mduc->chem->species,
                                      species[i]);
    } // end loop over CFD species
}
/* -------------------------------------------------------------------------- */

int  MDUCgetGridSize(int dir, void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  if (dir == 1)
    return mduc->flame->n1;
  else if (dir == 2)
    return mduc->flame->n2;
  return 0;
}
/* -------------------------------------------------------------------------- */

int  MDUCgetHdim(void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  driveStruct *driver = mduc->driver;
  return driver->H_DIM;
}

/* -------------------------------------------------------------------------- */

void MDUCgetGridBound(double *x1min, double *x1max, double *x2min,
                      double *x2max, void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  flowStruct  *flow   = (flowStruct*) ((mducStruct*)mem)->flow;

  if (x1min != NULL) *x1min = flow->phi1min;
  if (x1max != NULL) *x1max = flow->phi1max;
  if (x2min != NULL) *x2min = flow->phi2min;
  if (x2max != NULL) *x2max = flow->phi2max;
}
/* -------------------------------------------------------------------------- */

void MDUCsetGridBound(double x1min, double x1max, double x2min,
                      double x2max, void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  flowStruct  *flow   = (flowStruct*) ((mducStruct*)mem)->flow;

  flow->phi1min = x1min;
  flow->phi1max = x1max;
  flow->phi2min = x2min;
  flow->phi2max = x2max;

}
/* -------------------------------------------------------------------------- */

int  MDUCgetNspecies(void *mem)
{
  mducStruct *mduc = (mducStruct*) mem;
  return mduc->chem->nSpecies;
}
/* -------------------------------------------------------------------------- */

int  MDUCgetNspecs(void *mem)
{
  mducStruct *mduc = (mducStruct*) mem;
  return mduc->chem->nSpecSolve;
}
int  MDUCgetNsootMom(void *mem)
{
  mducStruct *mduc = (mducStruct*) mem;
  return mduc->chem->nSootMom;
}
/* -------------------------------------------------------------------------- */

void MDUCgetSpeciesNames(char ***names, int buf_len, void *mem)
{
  int i,n;
  mducStruct *mduc = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;

  n = chem->nSpecies;

  (*names) = (char**) calloc(n, sizeof(char*));
  for (i=0; i<n; i++)
    (*names)[i] = (char*) calloc(buf_len, sizeof(char));
  for (i=0; i<n; i++)
    {
      strncpy((*names)[i], chem->species[i],buf_len);
    }
}
/* -------------------------------------------------------------------------- */

void MDUCcomputeThermoData(double *h, double *cp, double T, void *mem)
{
  mducStruct *mduc = (mducStruct*) mem;
  chemComputeThermoData(h, cp, T, mduc->chem);
}
/* -------------------------------------------------------------------------- */

int  MDUCgetNdim(void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  return mduc->flame->nDim;
}
/* -------------------------------------------------------------------------- */

double *MDUCgetZrefV(void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  return mduc->flow->ZRef;
}
/* -------------------------------------------------------------------------- */

double *MDUCgetChiRefV(void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  return mduc->flow->chiRef;
}
/* -------------------------------------------------------------------------- */

double MDUCgetZref(int dir, void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  if (dir > 2)
    return -1.0;
  else
    return mduc->flow->ZRef[dir];
}
/* -------------------------------------------------------------------------- */
double MDUCgetChiZ1max(void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  return mduc->flow->chiZ1max;
}
/* -------------------------------------------------------------------------- */

double MDUCgetChiRef(int dir, void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  if (dir > 2)
    return -1.0;
  else
    return mduc->flow->chiRef[dir];
}
/* -------------------------------------------------------------------------- */

void MDUCgetZmoments(int dir, double *mean, double *var, void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  if (dir > 2)
    {
      *mean = -1.0;
      *var  = -1.0;
    }
  else
    {
      *mean = mduc->flow->phiMean[dir];
      *var  = mduc->flow->phiVar[dir];
    }
}
/* -------------------------------------------------------------------------- */

double MDUCgetStrainRate(int dir, void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  if (dir > 2)
    return -1.0;
  else
    return mduc->flow->strain[dir];
}
/* -------------------------------------------------------------------------- */

// this function and the next are really not both necessary for the general case
int MDUCgetGrid(double **grid, void *mem)
{
  flameStruct *flame  = (flameStruct*) ((mducStruct*)mem)->flame;
  
  *grid = (double*) calloc(flame->n1, sizeof(double));
  
  int i;
  for (i=0; i<flame->n1; i++)
    (*grid)[i] = flame->x1[i];
  
  return flame->n1;
}
/* -------------------------------------------------------------------------- */

int * MDUCgetGrid2D(double ***grid, void *mem)
{
  int i;
  int *gridSize;
  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;
  
  // allocate the memory for the grid
  gridSize = (int*) calloc(2, sizeof(int));
  *grid = (double**) calloc(flame->nDim, sizeof(double*));
  (*grid)[0] = (double*) calloc(flame->n1, sizeof(double));
  gridSize[0] = flame->n1;
  gridSize[1] = 1; // set default to single row
  
  for (i=0; i<flame->n1; i++)
    (*grid)[0][i] = flame->x1[i];
  
  if (flame->nDim == 2)
    {
      (*grid)[1] = (double*) calloc(flame->n2, sizeof(double));
      gridSize[1] = flame->n2;
      for (i=0; i<flame->n2; i++)
        (*grid)[1][i] = flame->x2[i];
    }
  
  return gridSize;
}
/* -------------------------------------------------------------------------- */

void MDUCgetZGrid(double ***Z1grid, double ***Z2grid, void *mem)
{
  int j, k;
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;
  
  *Z1grid = (double**) calloc(flame->n2, sizeof(double*));
  *Z2grid = (double**) calloc(flame->n2, sizeof(double*));
  for (k=0; k<flame->n2; ++k)
    {
      (*Z1grid)[k] = (double*) calloc(flame->n1, sizeof(double));
      (*Z2grid)[k] = (double*) calloc(flame->n1, sizeof(double));
    }
  
  for (k=0; k<flame->n2; ++k)
    {
      for (j=0; j<flame->n1; ++j)
        retrieveZ(flame->x1[j], flame->x2[k],&(*Z1grid)[k][j],&(*Z2grid)[k][j]);
    } // end loop over rows
}
/* -------------------------------------------------------------------------- */

double MDUCgetPressure(int n, void *mem)
{
  flowStruct  *flow   = (flowStruct*) ((mducStruct*)mem)->flow;
  if (n == -1)
    return flow->P_m1;
  else if (n == 0)
    return flow->P;
  else
    return flow->P; // return P by default
}
/* -------------------------------------------------------------------------- */

double MDUCgetHbounds(int n, void *mem)
{
  flowStruct  *flow   = (flowStruct*) ((mducStruct*)mem)->flow;
  if (n == -1)
    return flow->phi1min_m1;
  else if (n == 1)
    return flow->phi1max_m1;
  else
    return flow->phi1min_m1; // return phi1min by default
}
/* -------------------------------------------------------------------------- */

double MDUCgetRhobounds(int n, void *mem)
{
  flameStruct  *flame   = (flameStruct*) ((mducStruct*)mem)->flame;
  if (n == -1)
    return flame->rho[0];
  else if (n == 1)
    return flame->rho[flame->n1-1];
  else
    return flame->rho[0]; // return rho0 by default
}
/* -------------------------------------------------------------------------- */

void MDUCsetPressure(int n, double P, void *mem)
{
  flowStruct  *flow   = (flowStruct*) ((mducStruct*)mem)->flow;
  if (n == -1)
    flow->P_m1 = P;
  else if (n == 0)
    flow->P = P;
}
/* -------------------------------------------------------------------------- */

void MDUCsetTempBounds(double T1, double T2, double T3, void *mem)
{
  int k;
  flameStruct *flame = (flameStruct*) ((mducStruct*)mem)->flame;
  if (T1 > 0.0) /* xi = 0 */
    flame->Y[flame->sT] = T1;
  if (T2 > 0.0) /* x1 = 1 */
    for (k=0; k<flame->n2; k++)
      flame->Y[(k*flame->n1+flame->n1-1)*flame->nVars+flame->sT] = T2;
  if (flame->nDim == 2 && T3 > 0.0) { /* x2 = 1 */
    flame->Y[(flame->n2-1)*flame->n1*flame->nVars+flame->sT] = T3;
  }
}
/* -------------------------------------------------------------------------- */

double MDUCgetTime(int n, void *mem)
{
  timeStruct  *time   = (timeStruct*) ((mducStruct*)mem)->time;
  if (n == -1)
    return time->t_m1;
  else if (n == 0)
    return time->t;
  else
    return time->t; // return P by default
}
/* -------------------------------------------------------------------------- */

void MDUCsetDataDir(char *directory, void *mem)
{
  ioStruct *io = (ioStruct*) ((mducStruct*)mem)->io;
  strcpy(io->dataDir, directory);
}
/* -------------------------------------------------------------------------- */

void MDUCsetFilePrefix(char *prefix, void *mem)
{
  ioStruct *io = (ioStruct*) ((mducStruct*)mem)->io;
  char buffer[1024];

  /* construct file names */
  strcpy(io->outputData, prefix);
  strcat(io->outputData, "_.rif");
  strcpy(io->flowData, prefix);
  strcat(io->flowData, "_.flo");
  strcpy(io->scalars, prefix);
  strcat(io->scalars, "_scalars");

  /* add the path */
  if (io->dataDir[0] != '\0') {
    strcat(io->dataDir, "/");
    prefixString(io->dataDir, io->outputData);
    prefixString(io->dataDir, io->flowData);
    prefixString(io->dataDir, io->scalars);
    prefixString(io->dataDir, io->monitorFile);
  }
}
/* -------------------------------------------------------------------------- */

double MDUCgetDPDT(void *mem)
{
  flowStruct  *flow   = (flowStruct*) ((mducStruct*)mem)->flow;
  return flow->Dt_P;
}
/* -------------------------------------------------------------------------- */

void MDUCgetPointY(double x1, double x2, double *Y, void *mem)
{
  int i, j, k, jL, jR, kA, kB, rowA, rowB, colL, colR;

  // dereference flame for grid and solution
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  flameStruct *flame  = mduc->flame;

  jL = jR = kA = kB = 0;
  // find the indices to the left/right of x1
  do { jR++; } while (flame->x1[jR] < x1 && jR < flame->n1);
  jL = jR-1;

  colL = jL*flame->nVars;
  colR = jR*flame->nVars;

  if (flame->nDim == 1)
    {
      for (i=0; i<chem->nSpecSolve; i++)
        Y[i] = linInterp(x1, flame->x1[jL], flame->Y[colL+i], flame->x1[jR],
                         flame->Y[colR+i]);
    }
  if (flame->nDim == 2)
    {
      // find the indices above/below x2
      do { kA++; } while (flame->x2[kA] < x2 && kA < flame->n2);
      kB = kA-1;

      rowA = kA*flame->n1*flame->nVars;
      rowB = kB*flame->n1*flame->nVars;

      for (i=0; i<chem->nSpecSolve; i++)
        Y[i] = biLinInterp(x1, x2, flame->x1[jL], flame->x1[jR], flame->x2[kB],
                           flame->x2[kA], flame->Y[rowB+colL+i],
                           flame->Y[rowB+colR+i], flame->Y[rowA+colL+i],
                           flame->Y[rowA+colR+i]);
    }
}
/* -------------------------------------------------------------------------- */

double MDUCgetPointT(double x1, double x2, void *mem)
{
  int i, j, k, jL, jR, kA, kB, rowA, rowB, colL, colR;

  // dereference flame for grid and solution
  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;
  i = flame->sT;

  jL = jR = kA = kB = 0;
  // find the indices to the left/right of x1
  do { jR++; } while (flame->x1[jR] < x1 && jR < flame->n1);
  jL = jR-1;

  colL = jL*flame->nVars;
  colR = jR*flame->nVars;

  if (flame->nDim == 1)
    return linInterp(x1, flame->x1[jL], flame->Y[colL+i], flame->x1[jR],
                     flame->Y[colR+i]);

  // 2D solution
  // find the indices above/below x2
  do { kA++; } while (flame->x2[kA] < x2 && kA < flame->n2);
  kB = kA-1;

  rowA = kA*flame->n1*flame->nVars;
  rowB = kB*flame->n1*flame->nVars;

  return biLinInterp(x1, x2, flame->x1[jL], flame->x1[jR],
                     flame->x2[kB], flame->x2[kA],
                     flame->Y[rowB+colL+i],flame->Y[rowB+colR+i],
                     flame->Y[rowA+colL+i],flame->Y[rowA+colR+i]);
}
/* -------------------------------------------------------------------------- */

double MDUCgetPointH(double x1, double x2, void *mem)
{
  int i;
  // dereference flame for grid and solution
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  flameStruct *flame  = mduc->flame;

  double Hpoint, Tpoint;
  double *Ypoint = (double*) calloc(chem->nSpecies, sizeof(double));

  Tpoint = MDUCgetPointT(x1, x2, mem);
  MDUCgetPointY(x1, x2, Ypoint, mem);

  Hpoint = 0.0;
  chemComputeThermoData(chem->h, chem->cp, Tpoint, chem);
  for (i=0; i<chem->nSpecies; i++)
    Hpoint += Ypoint[i]*chem->h[i];

  free(Ypoint);

  return Hpoint;
}
/* -------------------------------------------------------------------------- */

void MDUCstoreCvodePre(void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  timeStruct  *time   = mduc->time;
  sdialStruct *sdial  = mduc->sdial;

  SDialGetLastStep(sdial->mem, &time->dt_cvode_pre);
}
/* -------------------------------------------------------------------------- */

void MDUCsetSolution(double t, double P, double ***Y, double **Temp, void *mem)
{
  int i, j, k, row, col, flag;
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  timeStruct  *time   = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;
  sdialStruct *sdial  = mduc->sdial;
  
  /* ----- copy given solution to the flamelet array ----- */
  for (k=0; k<flame->n2; k++) {
    row = k*flame->n1*flame->nVars;
    for (j=0; j<flame->n1; j++) {
      col = j*(flame->nVars);
      // loop over all species
      for (i=0; i<chem->nSpecSolve; i++)
        flame->Y[row+col+i] = Y[k][j][i];
      
      flame->Y[row+col+flame->sT] = Temp[k][j];
    } // end loop over columns
  } // end loop over rows

  // reset the previous pressure and time
  if (P > 0.0)  flow->P_m1 = flow->P =  P;
  if (t >= 0.0) time->t = t;

  // reinitialise if a 1D flamelet (2D done during solve anyway)
  MDUCresetODE(mem);
}
/* -------------------------------------------------------------------------- */

void MDUCresetODE(void *mem)
{
  int flag;
  mducStruct  *mduc   = (mducStruct*) mem;
  driveStruct *driver = mduc->driver;
  timeStruct  *time   = mduc->time;
  flameStruct *flame  = mduc->flame;
  sdialStruct *sdial  = mduc->sdial;
#ifdef ACCESS_MPI
  mpiStruct   *mpi    = mduc->mpi;
#endif

  /* if ADI scheme is used, reset is done in actual solve anyway */
  if (flame->nDim == 2 && sdial->ADI) return;

  if (sdial->mem == NULL || sdial->Ydata == NULL) {
#ifdef ACCESS_MPI
    if (mpi->myid == mpi->host) 
#endif /* ACCESS_MPI */
    	fprintf(stderr, "Warning: no Sundial object to reset\n");
    return;
  }

  int offset = 0;
  int nEq    = sdial->nEq;
#ifdef ACCESS_MPI
  /* TODO: should jStart be used for both 1D and 2D? */
  if (flame->nDim == 1) {
    offset = mpi->jStart[mpi->myid];
    nEq    = mpi->jCount[mpi->myid];
  }
  else if (flame->nDim == 2) {
    offset = mpi->kStart[mpi->myid];
    nEq    = mpi->kCount[mpi->myid];
  }
#endif /* ACCESS_MPI */

  /* make sure solver solution array is up-to-date */
  copyArray(nEq, &flame->Y[offset], sdial->Ydata);
  if (driver->SDIAL == 1)
    copyArray(nEq, &flame->Yp[offset], sdial->YpData);

  /* reinitialise cvode solver */
#ifdef SUNDIALS_2_3
  if (driver->SDIAL == 0)
    flag = CVodeReInit(sdial->mem, cvodeRHS, time->t, sdial->Y, CV_SS,
                       driver->RELTOL, &driver->ABSTOL);
  else if (driver->SDIAL == 1)
    flag = IDAReInit(sdial->mem, idaRHS, time->t, sdial->Y, sdial->Yp,
                     IDA_SS, driver->RELTOL, &driver->ABSTOL);
#else
  if (driver->SDIAL == 0)
    flag = CVodeReInit(sdial->mem, time->t, sdial->Y);
  else if (driver->SDIAL == 1)
    flag = IDAReInit(sdial->mem, time->t, sdial->Y, sdial->Yp);
#endif /* SUNDIALS_2_3 */
}
/* -------------------------------------------------------------------------- */

void MDUCgetSolution(double ***Y, double **Temp, void *mem)
{
  int i, j, k, row, col;
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  flameStruct *flame  = mduc->flame;
  
  /* ----- copy flamelet solution to flow solver solution arrays ----- */
  for (k=0; k<flame->n2; k++)
    {
      row = k*flame->n1*flame->nVars;
      for (j=0; j<flame->n1; j++)
        {
          col = j*flame->nVars;
          // loop over all species
          for (i=0; i<chem->nSpecSolve; i++)
            Y[k][j][i] = flame->Y[row+col+i];
          
          Temp[k][j] = flame->Y[row+col+flame->sT];
        } /* end loop over x1 */
    } /* end loop over x2 */
}
/* -------------------------------------------------------------------------- */
void MDUCgetSootCoeff(double **sqrtT, double **T_mu, double **muSqrtMw_rhoSqrtT,
		      double **sgCoeff, double **oxCoeff, double **o2Coeff,
		      double **dimProdRate, void *mem)
{
  int j, k, row, col;
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;
  chemStruct  *chem   = mduc->chem;

  if (chem->sootMode<3) 
     fprintf(stdout,"Wrong soot mode selected!\n");
  
  /* ----- copy flamelet solution to flow solver solution arrays ----- */
  for (k=0; k<flame->n2; k++)
    {
      row = k*flame->n1*flame->nVars;
      for (j=0; j<flame->n1; j++)
        {
          col = j;
          sqrtT[k][j] = flame->sqrtT[row+col];
          T_mu[k][j] = flame->T_mu[row+col];
          muSqrtMw_rhoSqrtT[k][j] = flame->muSqrtMw_rhoSqrtT[row+col];
          sgCoeff[k][j] = flame->sgCoeff[row+col];
          oxCoeff[k][j] = flame->oxCoeff[row+col];
          o2Coeff[k][j] = flame->o2Coeff[row+col];
          dimProdRate[k][j] = flame->dimProdRate[row+col];
        } /* end loop over x1 */
    } /* end loop over x2 */
}
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void MDUCgetSootPAH(double **yPAH, double **srcPAHpos, double **srcPAHneg,
		      double **rhoDot, void *mem)
{
  int j, k, row, col;
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;
 
  /* ----- copy flamelet solution to flow solver solution arrays ----- */
  for (k=0; k<flame->n2; k++)
    {
      row = k*flame->n1*flame->nVars;
      for (j=0; j<flame->n1; j++)
        {
          col = j*flame->nVars;
          yPAH[k][j]      = flame->yPAH[row+col];
          srcPAHpos[k][j] = flame->srcPAHpos[row+col];
          srcPAHneg[k][j] = flame->srcPAHneg[row+col];
          rhoDot[k][j]    = flame->rhoDotSoot[row+col];
        } /* end loop over x1 */
    } /* end loop over x2 */
}

void MDUCgetVar(char *name, double **var, void *mem)
{
  int i, j, k, r;
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;
  chemStruct  *chem   = mduc->chem;
  double *var_internal;

  if (!strcmp(name, "HRrate"))
    var_internal = flame->HRrate;
  else if (!strcmp(name, "enthalpy"))
    {
      /* Recalculate the enthalpy */
      calcTotalEnthalpy(flame->n1, flame->n2, flame->nVars, flame->sT, flame->Y,
                        flame->H, chem);
      var_internal = flame->H;
    }
  else if (!strcmp(name, "ChemSrc"))
    var_internal = flame->chem_src;
  
  /* ----- copy flamelet solution to flow solver solution arrays ----- */
  for (k=0; k<flame->n2; k++)
    {
      r = k*flame->n1;
      for (j=0; j<flame->n1; j++)
        {
          var[k][j] = var_internal[r+j];
        } /* end loop over columns */
    } /* end loop over rows */
}
/* -------------------------------------------------------------------------- */

void MDUCupdateInert1D(double dt, double P, double phimin, double phimax,
                       double Tox, double Tfuel, double *chi, double **Y,
                       double *Temp, void *mem)
{
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  timeStruct  *time   = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;
  sdialStruct *sdial  = mduc->sdial;

  int i, j, k, k_pre, iX1, flag;
  char outputFileName[256];
  
  double *Htot  = (double*) calloc(flame->n1, sizeof(double));
  double *Y_j   = (double*) calloc(chem->nSpecies, sizeof(double));
  double *hox   = (double*) calloc(chem->nSpecies, sizeof(double));
  double *hfuel = (double*) calloc(chem->nSpecies, sizeof(double));
#ifndef ACCESS_MPI 
  checkMonitor(io);
#endif

  // update time
  time->dt = dt;
  time->t += dt;
  time->nt++;
  
  // assign flow variables to struct
  flow->P       = P;
  flow->phi1min = phimin;
  flow->phi1max = phimax;

  // compute the pressure rate
  flow->Dt_P  = (flow->P - flow->P_m1)/time->dt;
  
  for (j=0; j<flame->n1; j++) flow->chi1[j] = chi[j];
  
  iX1 = (flame->n1-1)*flame->nVars; // location of start of fuel boundary
  
  chemComputeThermoData(hox, chem->cp, Tox, chem);
  chemComputeThermoData(hfuel, chem->cp, Tfuel, chem);
  
  Htot[0]= 0.0; // total enthalpy at oxidiser boundary
  Htot[flame->n1-1] = 0.0; // total enthalpy at fuel boundary
  for (i = 0; i < chem->nSpecies; i++)
    {
      Htot[0] += hox[i]*flame->Y[i];
      Htot[flame->n1-1] += hfuel[i]*flame->Y[iX1 + i];
    }
  
  // set the boundary temperatures
  flame->Y[flame->sT] = Tox;
  flame->Y[iX1 + flame->sT] = Tfuel;
  // update the interior profile
  for (j=1; j < flame->n1-1; j++)
    {
      k = j*flame->nVars;
      k_pre = (j-1)*flame->nVars;
      
      // linearly interpolate the enthalpy
      Htot[j] = linInterp(flame->x1[j], flame->x1[0], Htot[0],
                          flame->x1[flame->n1-1], Htot[flame->n1-1]);
      
      // get the species mass fractions for this mesh point
      for (i=0; i<chem->nSpecies; i++)
        Y_j[i] = flame->Y[k+i];
      
      // iterate for temperature based on the enthalpy and species
      flame->Y[k+flame->sT] = calcTfromH(chem->nSpecies, chem->nSpecSolve,
                                         Htot[j], Y_j,flame->Y[k_pre+flame->sT],
                                         chem);
    } // end loop over columns
  
  // reinitialise CVODE
  copyArray(sdial->nEq, flame->Y, sdial->Ydata);
#ifdef SUNDIALS_2_3
  if (driver->SDIAL == 0)
    flag = CVodeReInit(sdial->mem, cvodeRHS, time->t, sdial->Y, CV_SS,
                       driver->RELTOL, &driver->ABSTOL); 
  else if (driver->SDIAL == 1)
    {
      copyArray(sdial->nEq, flame->Yp_m1, sdial->YpData);
      flag = IDAReInit(sdial->mem, idaRHS, time->t, sdial->Y,
                       sdial->Yp, IDA_SS, driver->RELTOL, &driver->ABSTOL);
    }
#else
  if (driver->SDIAL == 0)
    flag = CVodeReInit(sdial->mem, time->t, sdial->Y); 
  else if (driver->SDIAL == 1)
    {
      copyArray(sdial->nEq, flame->Yp_m1, sdial->YpData);
      flag = IDAReInit(sdial->mem, time->t, sdial->Y, sdial->Yp);
    }
#endif

  // udpate io variables for writing files
  io->oxT = Tox;
  io->fuelT1 = Tfuel;

  if (io->outputHtot)
    calcTotalEnthalpy(flame->n1, flame->n2, flame->nVars, flame->sT, flame->Y,
                      flame->H, chem);
  
  // write solution file at this step (change output file name here if desired)
  if (time->nt%io->writeFreq == 0)
    {
      createDataFileName(io->outputData, io, time, outputFileName);
      writeFlamelet(outputFileName, io, driver, time, flame, flow);
    }
  
  // update normalisation for next time step
  flow->P_m1       = flow->P;
  flow->phi1min_m1 = flow->phi1min;
  flow->phi1max_m1 = flow->phi1max;
  time->dt_pre = dt;
  time->t_m1 = time->t;
  
  /* ----- copy flamelet solution to flow solver solution arrays ----- */
  for (j=0; j<flame->n1; j++)
    {
      k = j*(flame->nVars);
      // loop over all species
      for (i=0; i<chem->nSpecSolve; i++)
        Y[j][i] = flame->Y[k+i];
      
      Temp[j] = flame->Y[k+flame->sT];
    } // end mesh loop
  
  free(Htot);
  free(Y_j);
  free(hox);
  free(hfuel);
}
/* -------------------------------------------------------------------------- */

void MDUCflamelet1D(double dt,double P,double phimin,double phimax,double *chi,
                    int inert,double **Y,double *Temp,void *mem)
{
  int i,j,k;

  /* dereference internal flamelet solver structs */
  mducStruct   *mduc   = (mducStruct*) mem;
  chemStruct   *chem   = mduc->chem;
  flameStruct  *flame  = mduc->flame;

  int myid = 0;
#ifdef ACCESS_MPI
  MPI_Status comm_stat;
  mpiStruct *mpi = mduc->mpi; /* dereference mpi info */
  myid = mpi->myid;
#endif /* ACCESS_MPI */

  /* call the flamelet solver */
  if (myid >= 0)
    solveFlamelet1D(dt,P,phimin,phimax,chi,inert,mem);

  int fileId;
#ifdef ACCESS_MPI
  /* communicate solution to processors that are not solving the flamelet */
/* FIXME!!!
  if (myid == 0) {
    for (i=mpi->nProc; i<mpi->nWorld; i++)
      MPI_Send(flame->Y,flame->n1*flame->nVars,MPI_DOUBLE_PRECISION,
               i,i,mpi->world);
  }
  else {
    if (mpi->iWorld >= mpi->nProc)
      MPI_Recv(flame->Y,flame->n1*flame->nVars,MPI_DOUBLE_PRECISION,
               0,mpi->iWorld,mpi->world,&comm_stat);
  }
*/
  /* not sure why scatter doesn't work, when above clearly does */
  /*
  MPI_Scatter(flame->Y,flame->n1*flame->nVars,MPI_DOUBLE_PRECISION,
              flame->Y,flame->n1*flame->nVars,MPI_DOUBLE_PRECISION,0,mpi->world);
  */
#endif /* ACCESS_MPI */

  /* ----- copy flamelet solution to flow solver solution arrays ----- */
  for (j=0; j<flame->n1; j++) {
    k = j*(flame->nVars);
    /* loop over all species */
    for (i=0; i<chem->nSpecSolve; i++)
      Y[j][i] = flame->Y[k+i];
    
    Temp[j] = flame->Y[k+flame->sT];
  } /* end mesh loop */
}
/* -------------------------------------------------------------------------- */

void MDUCflamelet2D(double dt,double P,double phi1min,double phi1max,
                    double phi2min,double phi2max,double **chi1,double **chi12,
                    double **chi2,int inert,double ***Y,double **Temp,void *mem)
{
  int i,j,k,row,col;

  /* dereference internal flamelet solver structs */
  mducStruct   *mduc   = (mducStruct*) mem;
  chemStruct   *chem   = mduc->chem;
  flameStruct  *flame  = mduc->flame;
  sdialStruct  *sdial  = mduc->sdial;

  /* call the flamelet solver */
  if (sdial->ADI)
    /* TODO: need to update temperature if want to call inert from here */
    MDUCsolve2D(dt,P,phi1min,phi1max,phi2min,phi2max,chi1,chi12,chi2,
                Y,Temp,mem);
  else
    solveFlamelet2D(dt,P,phi1min,phi1max,phi2min,phi2max,chi1,chi12,chi2,
                    inert,mem);

  /* ----- copy flamelet solution to flow solver solution arrays ----- */
  for (k=0; k<flame->n2; k++) {
    row = k*flame->n1*flame->nVars;
    for (j=0; j<flame->n1; j++) {
      col = j*flame->nVars;
      for (i=0; i<chem->nSpecSolve; i++)
        Y[k][j][i] = flame->Y[row+col+i];
      
      Temp[k][j] = flame->Y[row+col+flame->sT];
    } /* end loop over columns */
  } /* end loop over rows */
}
/* -------------------------------------------------------------------------- */

void MDUCsolve1D(double dt, double P, double phimin, double phimax, double *chi,
                 double **Y, double *Temp, void *mem)
{
  /* keep for compatibility with old versions */
  MDUCflamelet1D(dt,P,phimin,phimax,chi,0,Y,Temp,mem);
}
/* -------------------------------------------------------------------------- */

void MDUCupdateInert2D(double dt, double P, double phi1min, double phi1max,
                       double phi2min, double phi2max, double Tox,double Tfuel1,
                       double Tfuel2, double **chi1, double **chi12,
                       double **chi2, double ***Y, double **Temp, void *mem)
{
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  timeStruct  *time   = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;
  
  int i, j, k, r, row, row_pre, col, k_pre, jN;
  int iZ2, jZ2, kZ2, iX1, iX2;
  double Z1, Z2, x1_N, Y_N;
  char outputFileName[256];
  
  /* variables for enthalpy */
  double *hox    = (double*) calloc(chem->nSpecies, sizeof(double));
  double *hfuel1 = (double*) calloc(chem->nSpecies, sizeof(double));
  double *hfuel2 = (double*) calloc(chem->nSpecies, sizeof(double));
  double **initH = (double**) calloc(flame->n2, sizeof(double*));
  for (i=0; i<flame->n2; ++i) 
    initH[i] = (double*) calloc(flame->n1, sizeof(double));
#ifndef ACCESS_MPI
  checkMonitor(io);
#endif

  // update time
  time->dt = dt;
  time->t += dt;
  time->nt++;
  
  // assign flow variables to struct
  flow->P       = P;
  flow->phi1min = phi1min;
  flow->phi1max = phi1max;
  flow->phi2min = phi2min;
  flow->phi2max = phi2max;

  // udpate io variables for writing files
  io->oxT    = Tox;
  io->fuelT1 = Tfuel1;
  io->fuelT2 = Tfuel2;

  // compute the pressure rate
  flow->Dt_P  = (flow->P - flow->P_m1)/time->dt;
  
  /* */
  for (k=0; k<flame->n2; k++)
    {
      r = k*flame->n1;
      for (j=0; j<flame->n1; j++)
        {
          flow->chi1[r+j]  = chi1[k][j];
          flow->chi12[r+j] = chi12[k][j];
          flow->chi2[r+j]  = chi2[k][j];
        }
    }

  /* solvePrestep2D(chi1, chi12, chi2, driver, flame, flow); */
  
  // set the indices for boundaries
  iX1 = (flame->n1-1)*flame->nVars;             // start of x1 = 1, x2 = 0
  iX2 = (flame->n2-1)*(flame->n1*flame->nVars); // x1 = 0, x2 = 1
  if (driver->TRANS_Z2 < 2 || driver->H_DIM != 0)
    {
      iZ2 = iX2;
      jZ2 = 0;
      kZ2 = flame->n2 - 1;
    }
  else if (driver->TRANS_Z2 == 2)
    {
      iZ2 = iX1 + iX2;
      jZ2 = flame->n1 - 1;
      kZ2 = flame->n2 - 1;
    }

  chemComputeThermoData(hox, chem->cp, Tox, chem);
  chemComputeThermoData(hfuel1, chem->cp, Tfuel1, chem);
  chemComputeThermoData(hfuel2, chem->cp, Tfuel2, chem);

  /* ----- initialise enthalpies at boundaries  ----- */
  initH[0][0]           = 0.0; // total enthalpy at oxidiser boundary
  initH[0][flame->n1-1] = 0.0; // total enthalpy at fuel 1 boundary
  initH[flame->n2-1][0] = 0.0; // total enthalpy at fuel 2 boundary
  for (i=0; i < chem->nSpecies; i++)
    {
      initH[0][0]           += hox[i]*flame->Y[i];
      initH[0][flame->n1-1] += hfuel1[i]*flame->Y[iX1+i];
      initH[flame->n2-1][0] += hfuel2[i]*flame->Y[iX2+i];
    }

  // set the boundary temperatures
  flame->Y[flame->sT]       = Tox;
  flame->Y[iX1 + flame->sT] = Tfuel1;
  flame->Y[iX2 + flame->sT] = Tfuel2;

  // update the solution along Z1 (Z2 = 0)
  for (j=1; j < flame->n1-1; j++)
    {
      k = j*flame->nVars;
      k_pre = (j-1)*flame->nVars;
      
      // linearly interpolate the enthalpy
      initH[0][j] = linInterp(flame->x1[j], flame->x1[0], initH[0][0],
                              flame->x1[flame->n1-1], initH[0][flame->n1-1]);
      
      // iterate for temperature based on the enthalpy and species
      flame->Y[k+flame->sT] =
        calcTfromH(chem->nSpecies, chem->nSpecSolve, initH[0][j],
                   &flame->Y[k], flame->Y[k_pre+flame->sT], chem);
    } // end loop over columns

  /* ----- loop over all the interior points ----- */
  for (k=1; k<flame->n2-1; k++)
    {
      row = k*(flame->n1*flame->nVars);         // start of current row
      row_pre = (k-1)*(flame->n1*flame->nVars); // start of previous row
      
      /* ----- set the left boundary ----- */
      if (driver->TRANS_Z2 == 2)
        { /* for this transformation, left boundary (Z1 = 0) is a 
             singularity. Therefore, set all values to (Z1,Z2) = (0,0) */
          initH[k][0] = initH[0][0];
        }
      else 
        { 
          initH[k][0] = linInterp(flame->x2[k], flame->x2[0], initH[0][0],
                                  flame->x2[flame->n2-1],initH[flame->n2-1][0]);
        }

      // get the temperature based on the enthalpy
      flame->Y[row+flame->sT] =
        calcTfromH(chem->nSpecies, chem->nSpecSolve, initH[k][0],
                   &flame->Y[row], flame->Y[row_pre+flame->sT], chem);
      
      /* ----- loop over interior points in x1 direction ----- */
      for (j=1; j<flame->n1-1; j++)
        {
          col = j*flame->nVars; // current row/column position
          
          // compute the actual Z2 for this point
          if (driver->TRANS_Z2 == 0)
            {
              Z1 = flame->x1[j];
              Z2 = flame->x2[k];
            }
          else if (driver->TRANS_Z2 == 1)
            {
              Z1 = flame->x1[j];
              Z2 = flame->x2[k]*(1.0-flame->x1[j]);
            }
          else if (driver->TRANS_Z2 == 2)
            {
              Z1 = flame->x1[j]*(1.0-flame->x2[k]);
              Z2 = flame->x1[j]*flame->x2[k];
            }

          // intersection of the mixing line with the x1 axis
          x1_N = Z1/(1-Z2);
          
          // find the location of the mesh point before the actual
          jN = gridLocate(flame->n1, flame->x1, x1_N);
          
          // interpolate the enthalpy at the value of x1_N between neighbours
          Y_N = linInterp(x1_N, flame->x1[jN], initH[0][jN],
                          flame->x1[jN+1], initH[0][jN+1]);
          initH[k][j] = Z2*initH[kZ2][jZ2] + (1-Z2)*Y_N;

          // get the temperature based on the enthalpy
          flame->Y[row+col+flame->sT] = 
            calcTfromH(chem->nSpecies,chem->nSpecSolve,initH[k][j],
                       &flame->Y[row+col], flame->Y[row_pre+col+flame->sT],
                       chem);
          
        } // end loop over columns
      
      /* ----- set the right boundary ----- */
      col = (flame->n1-1)*flame->nVars;
      if (driver->TRANS_Z2 == 0) // interpolate between Z1+Z2=0.5 and Z1=1
        initH[k][flame->n1-1] = 0.5*flame->x2[k]*initH[flame->n2-1][0]
          + (1 - 0.5*flame->x2[k])*initH[0][flame->n1-1];
      else if (driver->TRANS_Z2 == 1)
        // right boundary is singular: set all to value at (1,0)
        initH[k][flame->n1-1] = initH[0][flame->n1-1];
      else if (driver->TRANS_Z2 == 2) // Z1+Z2=1 boundary
        initH[k][flame->n1-1] = flame->x2[k]*initH[flame->n2-1][flame->n1-1]
          + (1-flame->x2[k])*initH[0][flame->n1-1];

      // get the temperature based on the enthalpy
      flame->Y[row+col+flame->sT] = 
        calcTfromH(chem->nSpecies, chem->nSpecSolve, initH[k][flame->n1-1],
                   &flame->Y[row+col],flame->Y[row_pre+col+flame->sT], chem);

    } // end loop over rows (x2)
  
  /* ----- set the top boundary (Z1+Z2 = 1) ----- */
  for (j=1; j<flame->n1; j++)
    {
      row_pre = (flame->n2-2)*(flame->n1*flame->nVars); // start of prev row 
      col = j*flame->nVars;     // current column

      if (driver->TRANS_Z2 == 0 && j != 0)
        {
          // interpolate between Z1+Z2=0.5 and Z1=1
          initH[flame->n2-1][j] = (1 - 0.5*flame->x1[j])*initH[flame->n2-1][0]
            + 0.5*flame->x1[j] *initH[0][flame->n1-1];
        }
      else if (driver->TRANS_Z2 == 1 && j != 0)
        { // boundary: Z1+Z2 = 1
          Z2  = 1 - flame->x1[j];  // current Z2 based on Z1+Z2=1              
          
          // interpolate the enthalpy along the boundary
          initH[flame->n2-1][j] = Z2*initH[flame->n2-1][0]
            + (1-Z2)*initH[0][flame->n1-1];
          
        }
      else if (driver->TRANS_Z2 == 2 && j != flame->n1-1)
        { // boundary: Z1 = 0
          initH[flame->n2-1][j] =
            linInterp(flame->x1[j], flame->x1[0], initH[flame->n2-1][0],
                      flame->x1[flame->n1-1], initH[flame->n2-1][flame->n1-1]);
        }

      // compute the temperature along the boundary from species and enthalpy
      flame->Y[iX2+col+flame->sT] =
        calcTfromH(chem->nSpecies,chem->nSpecSolve,initH[flame->n2-1][j],
                   &flame->Y[iX2+col],flame->Y[row_pre+col+flame->sT], chem);
    }
  
  if (io->outputHtot)
    calcTotalEnthalpy(flame->n1, flame->n2, flame->nVars, flame->sT, flame->Y,
                      flame->H, chem);

  // write solution file at this step (change output file name here if desired)
  if (time->nt%io->writeFreq == 0)
    {
      createDataFileName(io->outputData, io, time, outputFileName);
      writeFlamelet(outputFileName, io, driver, time, flame, flow);
    }

  /* MDUCwriteFlowData(chi1, chi2, NULL, NULL, mem); */
  
  // update normalisation for next time step
  flow->P_m1       = flow->P;
  flow->P_pre      = flow->P;
  flow->phi1min_m1 = flow->phi1min;
  flow->phi1max_m1 = flow->phi1max;
  flow->phi2min_m1 = flow->phi2min;
  flow->phi2max_m1 = flow->phi2max;
  time->dt_pre = dt;
  time->t_m1 = time->t;
  
  /* ----- copy flamelet solution to flow solver solution arrays ----- */
  for (k=0; k<flame->n2; k++)
    {
      row = k*flame->n1*flame->nVars;
      for (j=0; j<flame->n1; j++)
        {
          col = j*flame->nVars;
          // loop over all species
          for (i=0; i<chem->nSpecSolve; i++)
            Y[k][j][i] = flame->Y[row+col+i];
          
          Temp[k][j] = flame->Y[row+col+flame->sT];
        } // end loop over x1
    } // end loop over x2
  
  for (k=0; k<flame->n2; k++) free(initH[k]);
  free(initH);
  free(hox); free(hfuel1); free(hfuel2);
}
/* -------------------------------------------------------------------------- */

void MDUCupdateY_pre(void *mem)
{
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;

  copyArray(flame->n1*flame->n2*flame->nVars, flame->Y, flame->Y_pre);
}
/* -------------------------------------------------------------------------- */

void MDUCupdateHRrate(double **HR, double **H, void *mem)
{
  int j, k, row, col, r;
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;

  for (k=0; k<flame->n2; k++)
    {
      r = k*flame->n1;
      for (j=0; j<flame->n1; j++)
        {
          flame->HRrate[r+j] = HR[k][j];
          flame->H[r+j]      = H[k][j];
        }
    }
}
/* -------------------------------------------------------------------------- */

void MDUCsolve2D(double dt, double P, double phi1min, double phi1max,
                 double phi2min, double phi2max, double **chi1, double **chi12,
                 double **chi2, double ***Y, double **Temp, void *mem)
{
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  timeStruct  *time   = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;
  sdialStruct *sdial  = mduc->sdial;
  
  int i, j, k, node, tag = 0, flag = 0, subit = 1;
  int localIts, maxFails = 6000;
  int rowStart, rowEnd, colStart, colEnd;
  int r, row, col;                   // indices for row and column
  int squareGrid = 0;                // flag for whether n1 == n2
  int singularBound = 0;             // flag for whether bound is singular
  int index0D, solve0D = 0;          // flag to solve 0D at the bound
  long int nst;                      // number of solver steps
  char outputFileName[256];
  double tStart, tNow, tNext, tEnd;  // various time variables
  double *temp;
  double dtMax, dtNow, dtPre;        // time step sizes
  double *Cp_pre1;                   // storage for explicit mix. specific heat
  flameStruct *flameCol;             // for use when solving columns
  flowStruct  *flowCol;              // for use when solving columns
  sdialStruct *sdial0D;              // for use if singular point is present

  // set local tolerances to those specified by input file
  realtype abstol = driver->ABSTOL;
  realtype reltol = driver->RELTOL;

  /* ----- initialise whether parallel or not ----- */
  int myid, host, rxn;
#ifdef ACCESS_MPI
  mpiStruct *mpi = mduc->mpi; // dereference mpi info
  MPI_Status stat;  // for receive statements

  // get process ID
  MPI_Comm_rank(mpi->comm, &myid);
  mpi->myid = myid; // store just in case

  // assign the indices of the rows/columns solved on this process
  rowStart = mpi->rowStart[myid];
  rowEnd   = mpi->rowEnd[myid];
  colStart = mpi->colStart[myid];
  colEnd   = mpi->colEnd[myid];
  rxn = mpi->nProc-1;

#else
  myid = host = 0;

  // solve all rows/columns
  rowStart = colStart = 0;
  rowEnd   = flame->n2 - 1;
  colEnd   = flame->n1 - 1;
#endif // ACCESS_MPI
  
  // allocate local memory
  flameCol = (flameStruct*) calloc(1, sizeof(flameStruct));
  flowCol  = (flowStruct*)  calloc(1, sizeof(flowStruct));
  Cp_pre1  = (double*) calloc(flame->n1*flame->n2, sizeof(double));

  if (myid == 0)
    checkMonitor(io);
  
  // Copy old solution to pre array
  copyArray(flame->n1*flame->n2*flame->nVars, flame->Y, flame->Y_pre);

  /* if 2D mixture fraction problem, singularity exists at column boundary.
     Use 0D to solve homogeneous reactor for half time step */
  if (driver->H_DIM == 0 && driver->TRANS_Z2 > 0)
    {
      singularBound = 1;

      if (driver->TRANS_Z2 == 1) // standard transform
        {
          index0D = flame->n1-1; // set index of singular boundary

          // check if this node has singular boundary
          if (colEnd == index0D)
            {
              solve0D = 1;
              colEnd--;
            }
        }
      else if (driver->TRANS_Z2 == 2) // Hasse transform
        {
          index0D = 0;  // set index of singular boundary       

          // check if this node has singular boundary
          if (colStart == index0D)
            {
              solve0D = 1;
              colStart++;
            }
        }

      // allocate memory for solving 0D problem
      sdial0D = (sdialStruct*) calloc(1, sizeof(sdialStruct));
      sdial0D->args    = (void**) calloc(6, sizeof(void*));
      sdial0D->args[0] = driver;
      sdial0D->args[1] = chem;
      sdial0D->args[2] = flowCol;
      sdial0D->args[3] = flameCol;
      sdial0D->args[4] = NULL;
      sdial0D->args[5] = NULL;
      if (driver->SDIAL == 0)
        initCVODE(0, 1, flame->nVars, time->t, NULL, driver, sdial0D);
      else if (driver->SDIAL == 1)
        {
          flow->P  = fmax(flow->P, 1.0e5);
          initIDA(0, 1, flame->nVars, time->t, reltol, abstol, NULL, sdial0D);
        }
    }

  if (flame->n1 == flame->n2) // check whether grid is square
    squareGrid = 1;

  // update time step to take
  time->dt = dt;
  tStart   = time->t;
  tNow     = time->t;
  tEnd     = time->t_n = time->t + dt;

  // check to make sure the flow field has been initialised for the flamelet
  if (flow->P_m1 == 0.0)
    mducErr("initial pressure has not been set", "MDUCsolve2D");
  
  // assign flow variables to struct
  flow->P = P;
  flow->phi1min = phi1min;
  flow->phi1max = phi1max;
  flow->phi2min = phi2min;
  flow->phi2max = phi2max;

  if (driver->H_DIM == 0)
    {
      flow->phi2min = phi2min/(1.0 - phi1min);
      if (phi1max < 1.0)
        flow->phi2max = phi2max/(1.0 - phi1max);
      else
        flow->phi2max = 1.0;
    }
 
  // compute pressure rate
  flow->Dt_P = (flow->P - flow->P_m1)/time->dt;
  // get normalisation for variables
  flow->D_phi1 = flow->phi1max - flow->phi1min;
  flow->D_phi2 = flow->phi2max - flow->phi2min;

  if (driver->H_DIM == 1 && driver->DNS == 1)
    {
      flow->Dt_phi1min = (flow->phi1min_corr - flow->phi1min_m1)/time->dt;
      flow->Dt_phi1max = (flow->phi1max_corr - flow->phi1max_m1)/time->dt;
    }
  else
    {
      flow->Dt_phi1min = (flow->phi1min - flow->phi1min_m1)/time->dt;
      flow->Dt_phi1max = (flow->phi1max - flow->phi1max_m1)/time->dt;
    }
  if (driver->H_DIM == 2 && driver->DNS == 1)
    {
      flow->Dt_phi2min = (flow->phi2min_corr - flow->phi2min_m1)/time->dt;
      flow->Dt_phi2max = (flow->phi2max_corr - flow->phi2max_m1)/time->dt;
    }
  else
    {
      flow->Dt_phi2min = (flow->phi2min - flow->phi2min_m1)/time->dt;
      flow->Dt_phi2max = (flow->phi2max - flow->phi2max_m1)/time->dt;
    }

  // Clip
  if (flow->Dt_phi2min<0.0) flow->Dt_phi2min = 0.0;
  if (flow->Dt_phi1min<0.0) flow->Dt_phi1min = 0.0;
  if (flow->Dt_phi2max>0.0) flow->Dt_phi2max = 0.0;
  if (flow->Dt_phi1max>0.0) flow->Dt_phi1max = 0.0;

  // zero the variables for heat release computation
  for (k=0; k<flame->n2; k++)
    {
      row = k*(flame->n1*flame->nVars);
      r   = k*flame->n1;
      for (j=0; j<flame->n1; j++)
        {
          col = j*flame->nVars;
          for (i=0; i<flame->nVars; i++)
            flame->chemSrcY[row+col+i] = flame->chemSrcY_tmp[row+col+i] = 0.0;

          flame->HRrate[r+j] = 0.0;
        }
    }

  if ( (driver->H_DIM != 0) && (flow->avgSdot == NULL) )
    fprintf(stderr, "Error: average source term for enthalpy not set\n");

  /* compute quantities that are functions of domain only
     These should be invariant for purposes of flamelet solver
     NOTE: includes the transformed scalar dissipation rates and convective
     velocities */
  solvePrestep2D(chi1, chi12, chi2, driver, flame, flow);
  
  /* ----- assign pointers for column solve ----- */
  setFlameCol(flame, flameCol);
  setFlowCol(flow, flowCol);

  // initialise SUNDIALS solver for the rows
  if (driver->SDIAL == 0)
    initCVODE(flame->nDim, flame->n1, flame->nVars, time->t, &flame->Y[0],
              driver, sdial);
  else if (driver->SDIAL == 1)
    {
      fprintf(stderr, "Warning: 2D ida not fully tested \n");
      initIDA(flame->nDim, flame->n1, flame->nVars, time->t, reltol, abstol,
              &flame->Y[0], sdial);

      calcCpmix2D(chem, flame, flame->Y, Cp_pre1);

      if (driver->ENERGY)
        {
          // update to most current temperature and enthalpy
          for (k=0; k<flame->n2; k++)
            updateEnergy2D(driver,chem,flame,flow,0.0,0,flame->Y,flame->Yp_m1);
          temp = flame->T;
        }
      else
        { // temperature is part of solution
          temp = (double*) calloc(flame->n1*flame->n2, sizeof(double));
          extractVarArray(flame->nVars,flame->n1,flame->n2,flame->sT,flame->Y,
                          temp);
        }

      // compute the initial Yprime
      for (k=0; k<flame->n2; k++)
        {
          calcExplicit_2D(driver, chem, flame, flow, k, Cp_pre1, flame->Y_m1);
          calcY_RHS_imp(driver,chem,flame,flow,k,&temp[k*flame->n1],flame->Y,
                        flame->Yp_m1);

          if (driver->ENERGY)
            {
              // update RHS energy to enthalpy that was computed in calcY_RHS
              for (j=0; j<flame->n1; j++)
                flame->Yp_m1[(k*flame->n1+j)*flame->nVars+flame->sT] = 
                  flame->srcH[k*flame->n1+j];
            }
          else
            calcT_RHS_imp(driver, chem, flame, flow, k, flame->Y, flame->Yp_m1);
        }
    } // end sundials init for rows

  if (myid == 0 && time->nt == time->nt0) // print monitor header
    {
      fprintf(io->monitor, "%5s %15s %15s %15s %5s %15s %15s %5s %15s\n",
              "step", "time[s]", "P[Pa]", "end-time[s]", "ADI", "dt_max[s]",
              "dt[s]", "r/c", "time[s]");
      fflush(io->monitor);
    }

  /* ----- begin time loop ----- */
  while (fabs(tEnd - time->t) > 1.0e-12)
    {
      // compute the maximum time step on the host
      if (myid == 0)
        {
          /* make sure that nothing has decreased below zero */
          for (k=0; k<flame->n2; k++)
            {
              row = k*flame->n1*flame->nVars;
              for(j=0; j<flame->n1; j++)
                {
                  col = j*flame->nVars;
                  for (i=0; i<flame->nVars; i++)
                    flame->Y[row+col+i] = fmax(flame->Y[row+col+i], 0.0);
                }
            }

          // update the solution at N-1 to N
          // NOTE: should see if can combine with calcCpmix2D
          updateYpre(flame);
     
          // compute the mixture specific heat for explicit gradients
          calcCpmix2D(chem, flame, flame->Y_m1, Cp_pre1);
      
          // determine the maximum time step based on CFL
          dtMax = calcMaxTimeStep(driver, chem, flame, flow, time->CFL, Cp_pre1,
                                  flame->Y_m1);

          if (dtMax < tEnd - tNow) 
            dtNow = dtMax;
          else
            dtNow = tEnd - tNow;

          // write info to monitor
          if (subit == 1)
            fprintf(io->monitor, "%5i %15.6e %15.6e %15.6e ", time->nt, time->t,
                    flow->P, tEnd);
          else
            fprintf(io->monitor, "%54s","");
          
          fprintf(io->monitor, "%5i %15.6e %15.6e %5s ",subit,dtMax,dtNow,"r");
          fflush(io->monitor);
        }

#ifdef ACCESS_MPI
      /* -----  broadcast the time step to all processors ----- */
      MPI_Bcast(&dtNow, 1, MPI_DOUBLE_PRECISION, 0, mpi->comm);

      /* ----- transfer data to nodes ----- */
      // send current solution
      MPI_Scatterv(flame->Y_m1, mpi->kCountPlus, mpi->kPlus,
                   MPI_DOUBLE_PRECISION, &flame->Y_m1[mpi->kPlus[myid]],
                   mpi->kCountPlus[myid], MPI_DOUBLE_PRECISION, 0, mpi->comm);

      // send current solution derivative if using IDA
      if (driver->SDIAL == 1)
        MPI_Scatterv(flame->Yp_m1, mpi->kCountPlus, mpi->kPlus,
                     MPI_DOUBLE_PRECISION, &flame->Yp_m1[mpi->kPlus[myid]],
                     mpi->kCountPlus[myid], MPI_DOUBLE_PRECISION, 0, mpi->comm);

      // send explicit specific heat
      MPI_Scatterv(Cp_pre1, mpi->rCountPlus, mpi->rPlus, MPI_DOUBLE_PRECISION,
                   &Cp_pre1[mpi->rPlus[myid]], mpi->rCountPlus[myid],
                   MPI_DOUBLE_PRECISION, 0, mpi->comm);

      // send source terms required for heat release computation
      if (io->outputHR)
        MPI_Scatterv(flame->chemSrcY, mpi->kCountPlus, mpi->kPlus,
                     MPI_DOUBLE_PRECISION, &flame->chemSrcY[mpi->kPlus[myid]],
                     mpi->kCountPlus[myid], MPI_DOUBLE_PRECISION, 0, mpi->comm);

      // send source terms for temperature equation
      if (driver->output_tsrc)
        {
          MPI_Scatterv(flame->chem_src, mpi->rCount, mpi->rStart,
                       MPI_DOUBLE_PRECISION, &flame->chem_src[mpi->rPlus[myid]],
                       mpi->rCountPlus[myid], MPI_DOUBLE_PRECISION,0,mpi->comm);
          MPI_Scatterv(flame->conv_src, mpi->rCount, mpi->rStart,
                       MPI_DOUBLE_PRECISION, &flame->conv_src[mpi->rPlus[myid]],
                       mpi->rCountPlus[myid], MPI_DOUBLE_PRECISION,0,mpi->comm);
          MPI_Scatterv(flame->diff_src, mpi->rCount, mpi->rStart,
                       MPI_DOUBLE_PRECISION, &flame->diff_src[mpi->rPlus[myid]],
                       mpi->rCountPlus[myid], MPI_DOUBLE_PRECISION,0,mpi->comm);
        }
#endif // ACCESS_MPI

      /* ----- loop over rows for ADI ----- */
      driver->DIM = 1; // set the current solve direction
      for (k=rowStart; k <= rowEnd; k++)
        {
          row = k*(flame->n1*flame->nVars);
          sdial->row = k;
          r = k*flame->n1;
          tNow  = time->t;
          tNext = tNow + 0.5*dtNow;

          // set the stop time for this integration
          flag = SDialSetStopTime(sdial->mem, tNext);

          // compute and store off-direction and cross derivatives (explicit)
          calcExplicit_2D(driver, chem, flame, flow, k, Cp_pre1, flame->Y_m1);

          // re-initialise the SUNDIALS memory for the current row
          copyArray(sdial->nEq, &flame->Y_m1[row], sdial->Ydata);
#ifdef SUNDIALS_2_3
          if (driver->SDIAL == 0)
            flag = CVodeReInit(sdial->mem, cvodeRHS_2D, tNow, sdial->Y, CV_SS,
                               reltol, &abstol);
          else if (driver->SDIAL == 1)
            {
              copyArray(sdial->nEq, &flame->Yp_m1[row], sdial->YpData);
              flag = IDAReInit(sdial->mem, idaRHS_2D, tNow, sdial->Y,
                               sdial->Yp, IDA_SS, reltol, &abstol);
            }
#else
          if (driver->SDIAL == 0)
            flag = CVodeReInit(sdial->mem, tNow, sdial->Y);
          else if (driver->SDIAL == 1)
            {
              copyArray(sdial->nEq, &flame->Yp_m1[row], sdial->YpData);
              flag = IDAReInit(sdial->mem, tNow, sdial->Y, sdial->Yp);
            }
#endif /* SUNDIALS_2_3 */

          // call SUNDIALS to solve the line to time + 0.5*dt
          localIts = 0;
          do {
            if (driver->SDIAL == 0)
              flag = CVode(sdial->mem, tNext, sdial->Y, &tNow, sdial->itask);
            else if (driver->SDIAL == 1)
              flag = IDASolve(sdial->mem, tNext, &tNow, sdial->Y, sdial->Yp,
                              sdial->itask);

            if (flag < 0)
              {
                fprintf(io->monitor,"\nSUNDIALS failed on row %i: flag = %i\n",
                        k,flag);
                fflush(io->monitor);
              }

            // get solver data about last time step
            SDialGetNumSteps(sdial->mem, &nst);
            SDialGetLastStep(sdial->mem, &dtPre);

            // update sources
            if (io->outputHR || driver->output_tsrc) 
              {
                for (j=0; j<flame->n1; j++)
                  {
                    if (io->outputHR)
                      {
                        for (i=0; i<chem->nSpecSolve; i++)
                          flame->chemSrcY[row+j*flame->nVars+i] +=
                            flame->chemSrcY_tmp[row+j*flame->nVars+i]*dtPre;
                      }
                    if (driver->output_tsrc)
                      {
                        flame->chem_src[r+j] += flame->chem_temp[r+j]*dtPre;
                        flame->conv_src[r+j] += flame->conv_temp[r+j]*dtPre;
                        flame->diff_src[r+j] += flame->diff_temp[r+j]*dtPre;
                      }
                  } // end loop over columns for this row
              }

            if (localIts+1 > maxFails)
              {
                fprintf(io->monitor,
                        "\nSUNDIALS: row %i exceeded max iterations\n", k);
                fprintf(io->monitor, "         processor %i\n", myid);
                fflush(io->monitor);
                exit(2);
              }
            localIts++;
          } while (flag != sdial->iret);

          // store solution (from SUNDIALS -> flame)
          copyArray(sdial->nEq, sdial->Ydata, &flame->Y[row]);
          if (driver->SDIAL == 1)
            copyArray(sdial->nEq, sdial->YpData, &flame->Yp[row]);
        } // end loop over rows (x2)
      
      // update the time to n+1/2
      time->t = tNow; // should also be tNext
      if (myid == 0) {
        fprintf(io->monitor, "%15.6e\n", time->t);
        fflush(io->monitor);
      }
      
      /* ----- update solution to time n+1/2 ----- */
#ifdef ACCESS_MPI
      // retrieve solution parts from all nodes
      MPI_Gatherv(&flame->Y[mpi->kStart[myid]], mpi->kCount[myid],
                  MPI_DOUBLE_PRECISION, flame->Y, mpi->kCount, mpi->kStart,
                  MPI_DOUBLE_PRECISION, 0, mpi->comm);

      // retrieve solution derivative if using IDA
      if (driver->SDIAL == 1)
        MPI_Gatherv(&flame->Yp[mpi->kStart[myid]], mpi->kCount[myid],
                    MPI_DOUBLE_PRECISION, flame->Yp, mpi->kCount,
                    mpi->kStart, MPI_DOUBLE_PRECISION, 0, mpi->comm);

      // retrieve source terms for heat release computation
      if (io->outputHR)
        MPI_Gatherv(&flame->chemSrcY[mpi->kStart[myid]], mpi->kCount[myid],
                    MPI_DOUBLE_PRECISION, flame->chemSrcY, mpi->kCount,
                    mpi->kStart, MPI_DOUBLE_PRECISION, 0, mpi->comm);

      // retrieve temperature source terms
      if (driver->output_tsrc)
        {
          MPI_Gatherv(&flame->chem_src[mpi->rStart[myid]],mpi->rCount[myid],
                      MPI_DOUBLE_PRECISION, flame->chem_src, mpi->rCount,
                      mpi->rStart, MPI_DOUBLE_PRECISION, 0, mpi->comm);
          MPI_Gatherv(&flame->conv_src[mpi->rStart[myid]],mpi->rCount[myid],
                      MPI_DOUBLE_PRECISION, flame->conv_src, mpi->rCount,
                      mpi->rStart, MPI_DOUBLE_PRECISION, 0, mpi->comm);
          MPI_Gatherv(&flame->diff_src[mpi->rStart[myid]],mpi->rCount[myid],
                      MPI_DOUBLE_PRECISION, flame->diff_src, mpi->rCount,
                      mpi->rStart, MPI_DOUBLE_PRECISION, 0, mpi->comm);
        }
#endif // ACCESS_MPI
      
      // transpose solution vector such that columns become rows, and vice-versa
      // this allows use of same RHS function
      transposeRowCol(flame->n1, flame->n2, flame->nVars, flame->Y);
      transposeRowCol(flame->n1, flame->n2, 1, flame->chiA);
      transposeRowCol(flame->n1, flame->n2, 1, flame->chiB);
      transposeRowCol(flame->n1, flame->n2, 1, flame->chiC);
      transposeRowCol(flame->n1, flame->n2, 1, flame->v1);
      transposeRowCol(flame->n1, flame->n2, 1, flame->v2);
      transposeRowCol(flame->n1, flame->n2, 1, flame->rho);
      transposeRowCol(flame->n1, flame->n2, flame->nVars, flame->chemSrcY);
      transposeRowCol(flame->n1, flame->n2, 1, flow->chi1);
      transposeRowCol(flame->n1, flame->n2, 1, flow->chi2);
      transposeRowCol(flame->n1, flame->n2, 1, flow->chi1_m1);
      transposeRowCol(flame->n1, flame->n2, 1, flow->chi2_m1);
      transposeRowCol(flame->n1, flame->n2, 1, flow->chi12);
      transposeRowCol(flame->n1, flame->n2, 1, flow->chi12_m1);
      // Transpose driver-Diff
      int tmpv = driver->DIFF_MODE[0];
      driver->DIFF_MODE[0] = driver->DIFF_MODE[2];
      driver->DIFF_MODE[2] = tmpv;

      if (driver->ENERGY)
        {
          transposeRowCol(flame->n1, flame->n2, 1, flame->H);
          transposeRowCol(flame->n1, flame->n2, 1, flame->T);
        }
      if (driver->output_tsrc) 
        {
          transposeRowCol(flame->n1, flame->n2, 1, flame->chem_src);
          transposeRowCol(flame->n1, flame->n2, 1, flame->conv_src);
          transposeRowCol(flame->n1, flame->n2, 1, flame->diff_src);
        }
      if (driver->SDIAL == 1)
        transposeRowCol(flame->n1, flame->n2, flame->nVars, flame->Yp);
      
      if (myid == 0)
        {
          /* make sure that nothing has decreased below zero */
          for (k=0; k<flame->n2; k++)
            {
              row = k*flame->n1*flame->nVars;
              for(j=0; j<flame->n1; j++)
                {
                  col = j*flame->nVars;
                  for (i=0; i<flame->nVars; i++)
                    flame->Y[row+col+i] = fmax(flame->Y[row+col+i], 0.0);
                }
            }

          // store the explicit solution
          updateYpre(flameCol);
      
          // calculate the new mixture specific heat for the interim solution
          calcCpmix2D(chem, flameCol, flameCol->Y_m1, Cp_pre1);
        }

#ifdef ACCESS_MPI
      /* ----- transfer data to nodes ----- */
      // send current solution
      MPI_Scatterv(flame->Y_m1, mpi->jCountPlus, mpi->jPlus,
                   MPI_DOUBLE_PRECISION, &flame->Y_m1[mpi->jPlus[myid]],
                   mpi->jCountPlus[myid], MPI_DOUBLE_PRECISION, 0, mpi->comm);
      
      // send current solution derivative if using IDA
      if (driver->SDIAL == 1)
        MPI_Scatterv(flame->Yp_m1, mpi->jCountPlus, mpi->jPlus,
                     MPI_DOUBLE_PRECISION, &flame->Yp_m1[mpi->jPlus[myid]],
                     mpi->jCountPlus[myid], MPI_DOUBLE_PRECISION, 0, mpi->comm);
      
      // send explicit specific heat
      MPI_Scatterv(Cp_pre1, mpi->cCountPlus, mpi->cPlus, MPI_DOUBLE_PRECISION,
                   &Cp_pre1[mpi->cPlus[myid]], mpi->cCountPlus[myid],
                   MPI_DOUBLE_PRECISION, 0, mpi->comm);

      // send source terms required for heat release computation
      if (io->outputHR)
        MPI_Scatterv(flame->chemSrcY, mpi->jCountPlus, mpi->jPlus,
                     MPI_DOUBLE_PRECISION, &flame->chemSrcY[mpi->jPlus[myid]],
                     mpi->jCountPlus[myid], MPI_DOUBLE_PRECISION, 0, mpi->comm);

      // send source terms for temperature equation
      if (driver->output_tsrc)
        {
          MPI_Scatterv(flame->chem_src, mpi->cCount, mpi->cStart,
                       MPI_DOUBLE_PRECISION, &flame->chem_src[mpi->cPlus[myid]],
                       mpi->cCountPlus[myid], MPI_DOUBLE_PRECISION,0,mpi->comm);
          MPI_Scatterv(flame->conv_src, mpi->cCount, mpi->cStart,
                       MPI_DOUBLE_PRECISION, &flame->conv_src[mpi->cPlus[myid]],
                       mpi->cCountPlus[myid], MPI_DOUBLE_PRECISION,0,mpi->comm);
          MPI_Scatterv(flame->diff_src, mpi->cCount, mpi->cStart,
                       MPI_DOUBLE_PRECISION, &flame->diff_src[mpi->cPlus[myid]],
                       mpi->cCountPlus[myid], MPI_DOUBLE_PRECISION,0,mpi->comm);
        }
#endif // ACCESS_MPI
      
      // create a new CVODE memory if necessary (i.e. n1 != n2)
      if (!squareGrid)
        {
          // dellocate the current CVODE memory
          destroySDial(driver->SDIAL, sdial);

          if (driver->SDIAL == 0)
            initCVODE(flame->nDim,flameCol->n1,flameCol->nVars,time->t,
                      &flameCol->Y[0],driver,sdial);
          else if (driver->SDIAL == 1)
            initIDA(flame->nDim,flameCol->n1,flameCol->nVars,time->t,reltol,
                    abstol,&flameCol->Y[0],sdial);
        }

      /* ----- loop over columns for ADI ----- */
      driver->DIM = 2; // set the current solve direction
      sdial->args[2] = (int*) flowCol;
      sdial->args[3] = (int*) flameCol;
      for (j=colStart; j <= colEnd; j++)
        {
          col = j*(flameCol->n1*flameCol->nVars);
          sdial->row = j;
          r = j*flameCol->n1;
          tNow  = time->t;
          tNext = tNow + 0.5*dtNow;

          // set the stop time for this integration
          flag = SDialSetStopTime(sdial->mem, tNext);
          
          // compute the off-direction and cross-derivatives
          calcExplicit_2D(driver,chem,flameCol,flowCol,j,Cp_pre1,
                          flameCol->Y_m1);

          // re-initialise the CVODE memory for the current column
          copyArray(sdial->nEq, &flameCol->Y_m1[col], sdial->Ydata);
#ifdef SUNDIALS_2_3
          if (driver->SDIAL == 0)
            flag = CVodeReInit(sdial->mem, cvodeRHS_2D, tNow, sdial->Y, CV_SS,
                               reltol, &abstol);
          else if (driver->SDIAL == 1)
            {
              copyArray(sdial->nEq, &flame->Yp_m1[col], sdial->YpData);
              flag = IDAReInit(sdial->mem, idaRHS_2D, tNow, sdial->Y,
                               sdial->Yp, IDA_SS, reltol, &abstol);
            }
#else
          if (driver->SDIAL == 0)
            flag = CVodeReInit(sdial->mem, tNow, sdial->Y);
          else if (driver->SDIAL == 1)
            {
              copyArray(sdial->nEq, &flame->Yp_m1[col], sdial->YpData);
              flag = IDAReInit(sdial->mem, tNow, sdial->Y, sdial->Yp);
            }
#endif /* SUNDIALS_2_3 */

          // call CVODE to solve the line to time + 0.5*dt
          localIts = 0;
          do {
            if (driver->SDIAL == 0)
              flag = CVode(sdial->mem, tNext, sdial->Y, &tNow, sdial->itask);
            else if (driver->SDIAL == 1)
              flag = IDASolve(sdial->mem, tNext, &tNow, sdial->Y, sdial->Yp,
                              sdial->itask);
            if (flag < 0)
              {
                fprintf(io->monitor,"\nSUNDIALS failed: column %i: flag = %i\n",
                        j, flag);
                fflush(io->monitor);
              }
  
            // get solver data about last time step
            SDialGetNumSteps(sdial->mem, &nst);
            SDialGetLastStep(sdial->mem, &dtPre);

            // update sources
            if (io->outputHR || driver->output_tsrc) 
              {
                for (k=0; k<flameCol->n1; k++)
                  {
                    if (io->outputHR)
                      {
                        for (i=0; i<chem->nSpecSolve; i++)
                          flame->chemSrcY[col+k*flame->nVars+i] +=
                            flame->chemSrcY_tmp[col+k*flame->nVars+i]*dtPre;
                      }
                    if (driver->output_tsrc)
                      {
                        flame->chem_src[r+k] += flame->chem_temp[r+k]*dtPre;
                        flame->conv_src[r+k] += flame->conv_temp[r+k]*dtPre;
                        flame->diff_src[r+k] += flame->diff_temp[r+k]*dtPre;
                      }
                  }
              }

            if (localIts+1 > maxFails)
              {
                fprintf(io->monitor,
                        "\nSUNDIALS: column %i exceeded max iterations\n", j);
                fprintf(io->monitor, "         processor %i\n", myid);
                fflush(io->monitor);
                exit(2);
              }
            localIts++;
          } while (flag != sdial->iret);

          // store solution from CVODE memory
          copyArray(sdial->nEq, sdial->Ydata, &flameCol->Y[col]);
          if (driver->SDIAL == 1)
            copyArray(sdial->nEq, sdial->YpData, &flameCol->Yp[col]);
        } // end loop over columns (x1)

      /* ----- solve singular point as homogeneous reactor if necessary ----- */
      if (solve0D)
        {
          col = index0D*(flameCol->n1*flameCol->nVars);
          sdial0D->row = index0D;
          sdial0D->col = 0;
          tNow  = time->t;
          tNext = tNow + 0.5*dtNow;

          // set the stop time for this integration
          flag = SDialSetStopTime(sdial0D->mem, tNext);

          // copy the point solution into the array
          copyArray(sdial0D->nEq, &flameCol->Y_m1[col], sdial0D->Ydata);
#ifdef SUNDIALS_2_3
          if (driver->SDIAL == 0)
            flag = CVodeReInit(sdial0D->mem,cvodeRHS_0D,tNow,sdial0D->Y,CV_SS,
                               reltol, &abstol);
          else if (driver->SDIAL == 1)
            {
              copyArray(sdial0D->nEq,&flameCol->Yp_m1[col],sdial0D->YpData);
              flag = IDAReInit(sdial0D->mem, idaRHS_0D, tNow, sdial0D->Y,
                               sdial0D->Yp, IDA_SS, reltol, &abstol);
            }
#else
          if (driver->SDIAL == 0)
            flag = CVodeReInit(sdial0D->mem, tNow, sdial0D->Y);
          else if (driver->SDIAL == 1)
            {
              copyArray(sdial0D->nEq, &flame->Yp_m1[col], sdial0D->YpData);
              flag = IDAReInit(sdial0D->mem, tNow, sdial0D->Y, sdial0D->Yp);
            }
#endif /* SUNDIALS_2_3 */

          // call CVODE to solve the line to time + 0.5*dt
          do {
            if (driver->SDIAL == 0)
              flag = CVode(sdial0D->mem,tNext,sdial0D->Y,&tNow,sdial->itask);
            else if (driver->SDIAL == 1)
              flag = IDASolve(sdial0D->mem, tNext, &tNow, sdial0D->Y,
                              sdial0D->Yp, sdial0D->itask);
            
            if (flag < 0)
              printf("SUNDIALS failed column %i: flag = %i\n", index0D, flag);

            SDialGetLastStep(sdial->mem, &dtPre);

            // update sources
            if (io->outputHR || driver->output_tsrc) 
              {
                for (k=0; k<flameCol->n1; k++)
                  {
                    if (io->outputHR)
                      for (i=0; i<chem->nSpecSolve; i++)
                        flame->chemSrcY[col+index0D*flame->nVars+i] +=
                          flame->chemSrcY_tmp[col+index0D*flame->nVars+i]*dtPre;
                    if (driver->output_tsrc)
                      {
                        flame->chem_src[col] += flame->chem_temp[col]*dtPre;
                        flame->conv_src[col] += flame->conv_temp[col]*dtPre;
                        flame->diff_src[col] += flame->diff_temp[col]*dtPre;
                      }
                  }
              }
          } while (flag != sdial0D->iret);

          // store solution from CVODE memory
          copyArray(sdial0D->nEq, sdial0D->Ydata, &flameCol->Y[col]);
          if (driver->SDIAL == 1)
            copyArray(sdial0D->nEq,sdial0D->YpData,&flameCol->Yp[col]);
        } // end reactor solve
      
      /* ----- update solution to time n+1 ----- */
#ifdef ACCESS_MPI 
      // retrieve solution parts from all nodes
      MPI_Gatherv(&flame->Y[mpi->jStart[myid]], mpi->jCount[myid],
                  MPI_DOUBLE_PRECISION, flame->Y, mpi->jCount, mpi->jStart,
                  MPI_DOUBLE_PRECISION, 0, mpi->comm);

      // retrieve solution derivative if using IDA
      if (driver->SDIAL == 1)
        MPI_Gatherv(&flame->Yp[mpi->jStart[myid]], mpi->jCount[myid],
                    MPI_DOUBLE_PRECISION, flame->Yp, mpi->jCount,
                    mpi->jStart, MPI_DOUBLE_PRECISION, 0, mpi->comm);

      // retrieve source terms for heat release computation
      if (io->outputHR)
        MPI_Gatherv(&flame->chemSrcY[mpi->jStart[myid]], mpi->jCount[myid],
                    MPI_DOUBLE_PRECISION, flame->chemSrcY, mpi->jCount,
                    mpi->jStart, MPI_DOUBLE_PRECISION, 0, mpi->comm);

      // retrive temperature source terms
      if (driver->output_tsrc)
        {
          MPI_Gatherv(&flame->chem_src[mpi->cStart[myid]],mpi->cCount[myid],
                      MPI_DOUBLE_PRECISION, flame->chem_src, mpi->cCount,
                      mpi->cStart, MPI_DOUBLE_PRECISION, 0, mpi->comm);
          MPI_Gatherv(&flame->conv_src[mpi->cStart[myid]],mpi->cCount[myid],
                      MPI_DOUBLE_PRECISION, flame->conv_src, mpi->cCount,
                      mpi->cStart, MPI_DOUBLE_PRECISION, 0, mpi->comm);
          MPI_Gatherv(&flame->diff_src[mpi->cStart[myid]],mpi->cCount[myid],
                      MPI_DOUBLE_PRECISION, flame->diff_src, mpi->cCount,
                      mpi->cStart, MPI_DOUBLE_PRECISION, 0, mpi->comm);
        }
#endif // ACCESS_MPI

      // transpose back to original orientation
      transposeRowCol(flame->n2, flame->n1, flame->nVars, flame->Y);
      transposeRowCol(flame->n2, flame->n1, 1, flame->chiA);
      transposeRowCol(flame->n2, flame->n1, 1, flame->chiB);
      transposeRowCol(flame->n2, flame->n1, 1, flame->chiC);
      transposeRowCol(flame->n2, flame->n1, 1, flame->v1);
      transposeRowCol(flame->n2, flame->n1, 1, flame->v2);
      transposeRowCol(flame->n2, flame->n1, 1, flame->rho);
      transposeRowCol(flame->n2, flame->n1, flame->nVars, flame->chemSrcY);
      transposeRowCol(flame->n2, flame->n1, 1, flow->chi1);
      transposeRowCol(flame->n2, flame->n1, 1, flow->chi2);
      transposeRowCol(flame->n2, flame->n1, 1, flow->chi1_m1);
      transposeRowCol(flame->n2, flame->n1, 1, flow->chi2_m1);
      transposeRowCol(flame->n2, flame->n1, 1, flow->chi12);
      transposeRowCol(flame->n2, flame->n1, 1, flow->chi12_m1);
      // Transpose driver-Diff
      int tmpv2 = driver->DIFF_MODE[0];
      driver->DIFF_MODE[0] = driver->DIFF_MODE[2];
      driver->DIFF_MODE[2] = tmpv2;
      if (driver->ENERGY)
        {
          transposeRowCol(flame->n2, flame->n1, 1, flame->H);
          transposeRowCol(flame->n2, flame->n1, 1, flame->T);
        }
      if (driver->output_tsrc) 
        {
          transposeRowCol(flame->n2, flame->n1, 1, flame->chem_src);
          transposeRowCol(flame->n2, flame->n1, 1, flame->conv_src);
          transposeRowCol(flame->n2, flame->n1, 1, flame->diff_src);
        }
      if (driver->SDIAL == 1)
        transposeRowCol(flame->n2, flame->n1, flame->nVars, flame->Yp);
      
      // copy the singular point to the entire boundary
      if (singularBound)
        {
          for (k=1; k<flame->n2; k++)
            {
              r   = k*flame->n1;
              row = k*flame->n1*flame->nVars;
              col = index0D*flame->nVars;
          
              for (i=0; i<flame->nVars; i++)
                {
                  flame->Y[row+col+i] = flame->Y[col+i];
                  flame->chemSrcY[row+col+i] = flame->chemSrcY[col+i];
                }

              flame->rho[r+index0D]    = flame->rho[index0D];

             if (driver->output_tsrc)
                {
                  flame->chem_src[r+index0D] = flame->chem_src[index0D];
                  flame->diff_src[r+index0D] = flame->conv_src[index0D];
                  flame->conv_src[r+index0D] = flame->diff_src[index0D];
                }
            } // end loop over columns
        }
      
      // update time to n+1
      time->t = tNow;
      if (myid == 0) {
        fprintf(io->monitor, "%91s %5s %15.6e\n","", "c", time->t);
        fflush(io->monitor);
      }

      // create a new CVODE memory if necessary
      if ( (!squareGrid) && (time->t < tEnd) )
        {
          // dellocate the current SUNDIALS memory
          destroySDial(driver->SDIAL, sdial);

          // initialise the memory for column mesh size
          if (driver->SDIAL == 0)
            initCVODE(flame->nDim,flame->n1,flame->nVars,time->t,&flame->Y[0],
                      driver,sdial);
          else if (driver->SDIAL == 1)
            initIDA(flame->nDim,flame->n1,flame->nVars,time->t,reltol,abstol,
                    &flame->Y[0],sdial);
        }
      sdial->args[2] = (int*) flow;
      sdial->args[3] = (int*) flame;
      subit++;
    } // end time loop
 
#ifdef ACCESS_MPI
  // accumulate all required data for output
  if (io->outputRho)
    MPI_Gatherv(&flame->rho[mpi->rStart[myid]],mpi->rCount[myid],
                MPI_DOUBLE_PRECISION, flame->rho, mpi->rCount, mpi->rStart,
                MPI_DOUBLE_PRECISION, 0, mpi->comm);   
#endif // ACCESS_MPI
 
  // finish averaging of heat release rate over time step
  if (myid == 0 && io->outputHR) MDUCcalcHRrate(mem);

  // update normalisation for next time step
  flow->P_m1       = flow->P;
  flow->P_pre      = flow->P;
  flow->phi1min_m1 = flow->phi1min;
  flow->phi1max_m1 = flow->phi1max;
  flow->phi2min_m1 = flow->phi2min;
  flow->phi2max_m1 = flow->phi2max;
  
  // update time variables
  time->t_m1   = time->t;
  time->dt_pre = dt;
  time->nt++;

  // copy scalar dissipation rate
  copyArray(flame->n1*flame->n2, flow->chi1, flow->chi1_m1);
  copyArray(flame->n1*flame->n2, flow->chi2, flow->chi2_m1);
  copyArray(flame->n1*flame->n2, flow->chi12, flow->chi12_m1);

  if (myid == 0) {
    // compute total enthalpy if output
    if (io->outputHtot)
      calcTotalEnthalpy(flame->n1,flame->n2,flame->nVars,flame->sT,flame->Y,
                        flame->H,chem);

    // write flamelet file
    if (io->freqType == 0 && time->nt%io->writeFreq == 0) {
          createDataFileName(io->outputData, io, time, outputFileName);
          writeFlamelet(outputFileName, io, driver, time, flame, flow);
    }
    else if (io->freqType == 1 && (fmod(time->t,io->tFreq) < 0.5*time->dt ||
                                   time->t > io->tWrite+io->tFreq)) {
      createDataFileName(io->outputData, io, time, outputFileName);
      writeFlamelet(outputFileName, io, driver, time, flame, flow);
      io->tWrite = time->t;
    }
  
    /* ----- copy flamelet solution to flow solver solution arrays ----- */
    for (k=0; k<flame->n2; k++) {
      row = k*flame->n1*flame->nVars;
      for (j=0; j<flame->n1; j++) {
        col = j*flame->nVars;
        // loop over all species
        for (i=0; i<chem->nSpecSolve; i++)
          Y[k][j][i] = flame->Y[row+col+i];
        
        Temp[k][j] = flame->Y[row+col+flame->sT];
      } // end loop over x1
    } // end loop over x2
  }

  // deallocate CVODE memory
  if (sdial->mem != NULL) {
    destroySDial(driver->SDIAL, sdial);
  }
  if (singularBound && sdial0D->mem != NULL) {
    destroySDial(driver->SDIAL, sdial0D);
    free(sdial0D->args);
    free(sdial0D);
  }
  free(Cp_pre1); free(flameCol); free(flowCol);
}
/* -------------------------------------------------------------------------- */

double MDUCcalcResidual(void *mem)
{
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  timeStruct  *time   = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;
  sdialStruct *sdial  = mduc->sdial;

  int i,j,k,row,col;
  double L2norm, tmp;

  L2norm = 0.0;
  for (k=0; k<flame->n2; k++) {
    row = k*flame->n1*flame->nVars;
    for (j=0; j<flame->n1; j++) {
      col = j*flame->nVars;
      for (i=0; i<flame->nVars; i++) {
        tmp = flame->Y[row+col+i] - flame->Y_pre[row+col+i];
        L2norm += tmp*tmp;
      }
    }
  }

  if (L2norm > 0.0) return pow(L2norm, 0.5);
  return 0.0;
}
/* -------------------------------------------------------------------------- */

int  MDUCgetSpeciesIndex(char species[25], void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  return getSpeciesIndex(chem->nSpecies, chem->species, species);
}
/* -------------------------------------------------------------------------- */

double MDUCgetMolWeight(char species[25], void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  if (!strcmp(species, "dummy")) return 1.0;

  int i = getSpeciesIndex(chem->nSpecies, chem->species, species);
  return chem->W[i];
}
/* -------------------------------------------------------------------------- */

void MDUCgetMolWeightAll(double *W, void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  int i;
  for (i=0; i<chem->nSpecies; i++) W[i] = chem->W[i];
}
/* -------------------------------------------------------------------------- */

void  *MDUCloadChemTable(char *chemTabFile, int nDim, char *type)
{
  return (void*)loadChemTable(chemTabFile, nDim, type);
};
/* -------------------------------------------------------------------------- */

double MDUClookupChemTable(char *name, double *x, double *vars, void *tab)
{
  return lookupChemTable(name, x, vars, tab);
};
/* -------------------------------------------------------------------------- */

double MDUClookupChemTableIndex(int *index, double *x, double *vars, void *tab)
{
  return lookupChemTableIndex(index, x, vars, tab);
};
/* -------------------------------------------------------------------------- */

int    MDUCgetChemTabNdim(void *tab)     { return getChemTabNdim(tab); };
int    MDUCgetChemTabNvar(void *tab)     { return getChemTabNvar(tab); };
int    MDUCgetChemTabNspecies(void *tab) { return getChemTabNspecies(tab); };
double MDUCgetChemTabPref(void *tab)     { return getChemTabPref(tab); };
double MDUCgetChemTabTref(void *tab)     { return getChemTabTref(tab); };
void   MDUCgetChemTabVarNames(char **names, void *tab)
{ getChemTabVarNames(names, tab); }
void   MDUCgetChemTabSpecies(char **species, void *tab)
{ getChemTabSpecies(species, tab); }

double MDUCchemTableVarMin(char *name, void *tab)
{ 
  return chemTableVarMin(name, tab);
};
/* -------------------------------------------------------------------------- */

double MDUCchemTableVarMax(char *name, void *tab)
{
  return chemTableVarMax(name, tab);
};
/* -------------------------------------------------------------------------- */

double MDUCchemTableProgMax(double *x, void *tab)
{ 
  double V = chemTableProgMax(x, tab);
  return V;
}
/* -------------------------------------------------------------------------- */

double MDUCchemTableProgMin(double *x, void *tab)
{ 
  double V = chemTableProgMin(x, tab);
  return V;
}
/* -------------------------------------------------------------------------- */

double MDUCchemTableProgScaleP(double P, double *x, void *tab)
{
  double prog_scaled = chemTableProgScaleP(P, x, tab);
  return prog_scaled;
}
/* -------------------------------------------------------------------------- */

void MDUCwriteFlamelet(char *fileName, void *mem)
{
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  timeStruct  *time   = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;

  if (fileName == NULL)
    {
      char outputFileName[256];
      createDataFileName(io->outputData, io, time, outputFileName);
      writeFlamelet(outputFileName, io, driver, time, flame, flow);
    }
  else
    writeFlamelet(fileName, io, driver, time, flame, flow);
}
/* -------------------------------------------------------------------------- */
void MDUCwriteRates(char *fileName, void *mem)
{
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem     = mduc->chem;

  writeRates(fileName, chem);
}
/* -------------------------------------------------------------------------- */

void MDUCprintStream(void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  printStreamProperties(io);
}
/* -------------------------------------------------------------------------- */

void MDUClistSpecies(void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  printSpeciesList(chem);
}
/* -------------------------------------------------------------------------- */

void MDUCsetMonitor(int nY, int nZ, char **monitorSpecies, void *mem)
{
  int i;
  
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  chemStruct  *chem   = mduc->chem;
  
  if (chem->species == NULL)
    {
      fprintf(stderr,"Error: species id must be initialised from mechanism\n");
      return;
    }
  
  // free any memory if it is already set previously
  if (io->monitorIndex)
    {
      for (i=0; i<io->nMonitorY; i++) free(io->monitorSpecies[i]);
      free(io->monitorSpecies);
      free(io->monitorIndex);
    }
  
  // set the number of species and mesh points
  io->nMonitorY = nY;
  io->nMonitorZ = nZ;
  
  // allocate and get index of species to monitor
  io->monitorIndex = (int*) calloc(nY, sizeof(int));
  io->monitorSpecies = (char**) calloc(nY, sizeof(char*));
  for (i=0; i<nY; ++i)
    io->monitorSpecies[i] = (char*) calloc(25, sizeof(char));
  
  for (i=0; i<nY; ++i)
    {
      strcpy(io->monitorSpecies[i], monitorSpecies[i]);
      io->monitorIndex[i] = getSpeciesIndex(chem->nSpecies, chem->species,
                                            monitorSpecies[i]);
    }
}
/* -------------------------------------------------------------------------- */

/* ----- access functions for statistics computations ----- */
int MDUCcalcPDF(double mean1, double var1, double mean2, double var2,double var,
                double **pdf, int **pdfBound, void *mem)
{
  int pdfFlag;

  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;


  if (flame->nDim == 1)
    {
      pdfFlag = MDUCcalcPDF_1D(io->pdfMode, flame->n1, flame->x1, mean1, var1,
                               pdf[0], pdfBound[0]);
    }
  else if (flame->nDim == 2)
    {
      int     *n = (int*)     calloc(2, sizeof(int));
      double **x = (double**) calloc(2, sizeof(double*));

      n[0] = flame->n1;
      n[1] = flame->n2;
      x[0] = flame->x1;
      x[1] = flame->x2;

      pdfFlag = MDUCcalcPDF_2D(io->pdfMode, n, x, mean1, var1, mean2, var2, var,
                               pdf, pdfBound, driver->TRANS_Z2);

      free(n);
      free(x);
    }

  return pdfFlag;
}
/* -------------------------------------------------------------------------- */

int MDUCcalcPDF_1D(int method, int n, double *x, double mean, double var,
                   double *pdf, int *pdfBound)
{
  int pdfFlag;
  if (method==PDF_TYPE_BETA) {
    /*pdfFlag = calcBetaPDF(n, x, mean, var, pdf, pdfBound);*/
    pdfFlag = calcBetaPDFall(n, x, mean, var, pdf, pdfBound);
  }
  else if (method==PDF_TYPE_SML)
    pdfFlag = calcSML_PDF(n, x, mean, var, pdf, pdfBound);
  else if (method==PDF_TYPE_GAUSS)
    pdfFlag = calcGaussPDF(n, x, mean, var, pdf, pdfBound);
  else if (method==PDF_TYPE_LOG)
    pdfFlag = calcLogNormPDF(n, x, mean, var, pdf, pdfBound);
  else
    {
      fprintf(stderr, "calcPDF_1D: method %i not a valid option\n", method);
      exit(2);
    }
  return pdfFlag;
}
/* -------------------------------------------------------------------------- */

int MDUCcalcPDF_2D(int method, int *n, double **x, double mean1, double var1,
                   double mean2, double var2, double var, double **pdf,
                   int **pdfBound, int transform)
{
  int pdfFlag;

  if (method==PDF_TYPE_DIRICHLET)
    pdfFlag = calcDirichletPDF(n,x,mean1,var1,mean2,var2,var,pdf,pdfBound,
                               transform);
  else if (method==PDF_TYPE_JBETA)
    pdfFlag = calcJointBetaPDF(n,x,mean1,var1,mean2,var2,var,pdf,pdfBound,
                               transform);
  else if (method==PDF_TYPE_SML)
    pdfFlag = calcSML_PDF2D(n,x,mean1,var1,mean2,var2,var,pdf,pdfBound,
                            transform);
  else if (method==PDF_TYPE_GAUSS)
    pdfFlag = calcGaussPDF2D(n,x,mean1,var1,mean2,var2,var,pdf,pdfBound,
                             transform);
  else if (method==PDF_TYPE_LOG)
    pdfFlag = calcLogNormPDF2D(n,x,mean1,var1,mean2,var2,var,pdf,pdfBound,
                             transform);
  /*else if (!strcmp(method, "independent") || !strcmp(method, "independent-beta")
           || !strcmp(method, "independent-sml"))*/
  else if (method==PDF_TYPE_INDEPENDENT_BETA
           || method==PDF_TYPE_INDEPENDENT_SML)
    {
      int    j, k;
      int    n1  = n[0];
      int    n2  = n[1];
      double *x1 = x[0];
      double *x2 = x[1];
      double *pdf1 = (double*) calloc(n[0], sizeof(double));
      double *pdf2 = (double*) calloc(n[0], sizeof(double));

      // compute 1D betaPDF in each direction
      if (method==PDF_TYPE_INDEPENDENT_BETA)
        {
          calcBetaPDF(n[0], x[0], mean1, var1, pdf1, pdfBound[0]);
          calcBetaPDF(n[1], x[1], mean2, var2, pdf2, pdfBound[0]);
        }
      else if (method==PDF_TYPE_INDEPENDENT_SML)
        {
          calcSML_PDF(n[0], x[0], mean1, var1, pdf1, pdfBound[0]);
          calcSML_PDF(n[1], x[1], mean2, var2, pdf2, pdfBound[0]);
        }

      // do pointwise multiplication to get a 2D PDF
      for (k=0; k<n2; k++)
        {
          for (j=0; j<n1; j++)
            {
              pdf[k][j] = pdf1[j]*pdf2[k];
            } // end loop over columns

          pdfBound[k][0] = 0;
          pdfBound[k][1] = n1-1;
        } // end loop over rows 

      pdfFlag = -1;

      free(pdf1);
      free(pdf2);
    }
  else
    {
      fprintf(stderr, "calcPDF_2D: method %i not a valid option\n", method);
      exit(2);
    }

  return pdfFlag;
}
/* -------------------------------------------------------------------------- */

int MDUCcalcMeanY(double mean1, double var1, double mean2, double var2,
                  double var, double ***Y, double *Ymean, double **pdf,
                  int **pdfBound, void *mem)
{
  int i, j, k, row, col;
  int pdfFlag, mducY = 0;

  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;

  // if the solution from the flamelet should be used, allocate memory
  if (Y == NULL)
    {
      mducY = 1;
      Y = (double***) calloc(flame->n2, sizeof(double**));
      for (k=0; k<flame->n2; k++)
        {
          row = k*flame->n1*flame->nVars;
          Y[k] = (double**) calloc(flame->n1, sizeof(double*));
          for (j=0; j<flame->n1; j++)
            {
              col = j*flame->nVars;
              Y[k][j] = (double*) calloc(chem->nSpecSolve, sizeof(double));

              for (i=0; i<chem->nSpecSolve; i++)
                Y[k][j][i] = flame->Y[row+col+i];
            } // end loop over columns
        } // end loop over rows
    }

  // compute the mean
  if (flame->nDim == 1)
    {
      pdfFlag = MDUCcalcMeanY_1D(io->pdfMode, flame->n1, flame->x1,
                                 chem->nSpecSolve, mean1, var1, Y[0], Ymean,
                                 pdf[0], pdfBound[0]);
    }
  else if (flame->nDim == 2)
    {
      int     *n = (int*) calloc(2, sizeof(int));
      double **x = (double**) calloc(2, sizeof(double*));

      n[0] = flame->n1;
      n[1] = flame->n2;
      x[0] = flame->x1;
      x[1] = flame->x2;

      pdfFlag = MDUCcalcMeanY_2D(io->pdfMode, n, x, chem->nSpecSolve, mean1,
                                 var1, mean2, var2, var, Y, Ymean, pdf,
                                 pdfBound, driver->TRANS_Z2);

      free(n);
      free(x);
    }

  if (mducY)
    {
      for (k=0; k<flame->n2; k++)
        {
          for (j=0; j<flame->n1; j++)
            free(Y[k][j]);

          free(Y[k]);
        }
      free(Y);
    }

  return pdfFlag;
}
/* -------------------------------------------------------------------------- */

int MDUCcalcMeanY_1D(int method, int n, double *x, int nSpecies, double mean,
                     double var, double **Y, double *Ymean, double *pdf,
                     int *pdfBound)
{
  int i, j, pdfFlag;
  double alpha = 1.0, beta = 1.0, tmp;
  double *Yi;

  Yi = (double*) calloc(n*nSpecies, sizeof(double));

  pdfFlag = MDUCcalcPDF_1D(method, n, x, mean, var, pdf, pdfBound);

  if (pdfFlag == 0)
    {
      // ensure that the mean and variance are properly bounded
      mean  = fmin(1.0, fmax(0.0, mean));
      var   = fmax(0.0, var);
        
      // compute the shape parameters for the beta distribution
      tmp   = mean*(1-mean)/var - 1;
      alpha = mean*tmp;
      beta  = (1-mean)*tmp;
    }

  // put solution into 1D array for mean computation function
  convertArrayMD_1D(nSpecies, n, 1, &Y, Yi);

  // loop over species
/*  for (i=0; i<nSpecies; i++)
    Ymean[i] = calcLinMeanY(i, nSpecies, alpha, beta, n, x, Yi, pdf, pdfBound,
                            pdfFlag);*/
  calcLinAllMeanY(nSpecies, alpha, beta, n, x, Yi, Ymean, pdf, pdfBound,
                            pdfFlag);

  free(Yi);
  return pdfFlag;
}
/* -------------------------------------------------------------------------- */

int MDUCcalcMomentsPDF_2D(int *n, double **x, double **X, double *meanPDF,
                          double *varPDF, double **pdf, int **pdfBound,
                          int pdfFlag, int transform)
{
  int i,j,k;
  double beta[4] = {-1.0, -1.0, -1.0, -1.0};

  int     n1 = n[0];
  int     n2 = n[1];
  double *x1 = x[0];
  double *x2 = x[1];

  // compute the total mean and the mean of each direction
  for (i=0; i<3; i++)
    meanPDF[i] = calcLinMeanY_2D(0, 1, n, x, X[i], beta[1], beta[2],
                                 beta[3], pdf, pdfBound, pdfFlag,
                                 transform);

  // compute all variances
  for (i=0; i<3; i++)
    varPDF[i] = calcVar(n,0,1,meanPDF[i],meanPDF[i],X[i],X[i],pdf);

  // compute the covariance
  varPDF[3] = calcVar(n,0,1,meanPDF[1],meanPDF[2],X[1],X[2],pdf);
  return 0;
}
/* -------------------------------------------------------------------------- */

double MDUCcalcMeanYi(int n, double *x, int iY, double mean, double var,
                      double **Y, double *pdf, int *pdfBound, int pdfFlag)
{
  int i, j;
  double alpha = 1.0, beta = 1.0, tmp;
  double *Yi, Ymean;

  if (pdfFlag == 0)
    {
      // ensure that the mean and variance are properly bounded
      mean  = fmin(1.0, fmax(0.0, mean));
      var   = fmax(0.0, var);
        
      // compute the shape parameters for the beta distribution
      tmp   = mean*(1-mean)/var - 1;
      alpha = mean*tmp;
      beta  = (1-mean)*tmp;
    }

  Yi = (double*) calloc(n, sizeof(double));
  for (j=0; j<n; j++)
    Yi[j] = Y[j][iY];

  Ymean = calcLinMeanY(0, 1, alpha, beta, n, x, Yi, pdf, pdfBound, pdfFlag);

  free(Yi);
  return Ymean;
}
/* -------------------------------------------------------------------------- */

double MDUCcalcMeanVar(char *varName, double mean1, double var1, double mean2,
                       double var2, double var0, double **pdf, int **pdfBound,
                       int pdfFlag, void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;
  driveStruct *driver = mduc->driver;

  int    i, j, index;
  double Ymean, tmp;

  if (!strncmp(varName, "temp", 4))
    index = flame->sT;
  else if (!strncmp(varName, "pres", 4))
    index = flame->sP;
  else
    index = MDUCgetSpeciesIndex(varName, mem);

  if (flame->nDim == 1)
    {
      double alpha = 1.0, beta = 1.0;
      if (pdfFlag == 0)
        {
          // ensure that the mean and variance are properly bounded
          mean1  = fmin(1.0, fmax(0.0, mean1));
          var1   = fmax(0.0, var1);
        
          // compute the shape parameters for the beta distribution
          tmp   = mean1*(1-mean1)/var1 - 1;
          alpha = mean1*tmp;
          beta  = (1-mean1)*tmp;
        }

      Ymean = calcLinMeanY(index, flame->nVars, alpha, beta, flame->n1,
                           flame->x1, flame->Y, pdf[0], pdfBound[0], pdfFlag);
    }
  else if (flame->nDim == 2)
    {
      double mean, S, Q, beta1 = 1.0, beta2 = 1.0, beta3 = 1.0;
      int     *n = (int*)     calloc(2, sizeof(int));
      double **x = (double**) calloc(2, sizeof(double*));
      n[0] = flame->n1; n[1] = flame->n2;
      x[0] = flame->x1; x[1] = flame->x2;

      if (pdfFlag == 0)
        {
          // ensure that the mean and variance are properly bounded
          mean1 = fmin(1.0, fmax(0.0, mean1));
          mean2 = fmin(1.0, fmax(0.0, mean2));
          mean  = mean1 + mean2;
          var1  = fmax(0.0, var1);
          var2  = fmax(0.0, var2);
          var0  = fmax(0.0, var0);
          
          // compute mean and variance parameters
          S = mean1*mean1 + mean2*mean2 + (1-mean)*(1-mean);
          Q = var1 + var2 + var0;
          
          // compute the shape parameters for the beta distribution
          tmp   = (1.0 - S)/Q - 1.0;
          beta1 = mean1*tmp;
          beta2 = mean2*tmp;
          beta3 = (1.0-mean)*tmp;
        }
      
      Ymean = calcLinMeanY_2D(index,flame->nVars,n,x,flame->Y,beta1,beta2,beta3,
                               pdf,pdfBound,pdfFlag,driver->TRANS_Z2);
    }
  else if (flame->nDim == 0)
    {
      Ymean = flame->Y[index];
    }

  return Ymean;
}
/* -------------------------------------------------------------------------- */

double MDUCcalcMeanDensity(double P, double **TIn, double **pdfIn,
                           int **pdfBound, int pdfFlag, void *mem)
{
  int i, j, k, row, col, r;
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;

  double rho, mean, wmix;
  double *pdf, *T;

  if (P <= 0.0) P = flow->P;

  if (pdfIn == NULL)
    {
      pdf = flow->pdf;
      pdfBound = (int**) calloc(flame->n2, sizeof(int*));
      for (k=0; k<flame->n2; k++)
        {
          pdfBound[k] = (int*) calloc(2, sizeof(int));
          pdfBound[k][0] = 0;
          pdfBound[k][1] = flame->n1-1;
        }
    }
  else
    {
      pdf = (double*) calloc(flame->n1*flame->n2, sizeof(double));
      convertArrayMD_1D(flame->n1, flame->n2, 1, &pdfIn, pdf);
    }

  T = (double*) calloc(flame->n1*flame->n2, sizeof(double));
  if (TIn == NULL)
    extractVarArray(flame->nVars,flame->n1,flame->n2,flame->sT,flame->Y,T);
  else
    convertArrayMD_1D(flame->n1, flame->n2, 1, &TIn, T);
  
  mean = 0.0;
  for (k=0; k<flame->n2; k++)
    {
      row = k*flame->n1*flame->nVars;
      r   = k*flame->n1;
      for (j=0; j<flame->n1; j++)
        {
          col   = j*flame->nVars;
          wmix  = calcMixtureW(chem->nSpecSolve, &flame->Y[row+col], chem->W);
          rho   = wmix*P/(flame->RGAS*T[r+j]);
          mean += rho*pdf[r+j];
        } // end loop over columns
    } // end loop over rows

  if (pdfIn == NULL)
    {
      for (k=0; k<flame->n2; k++) free(pdfBound[k]);
      free(pdfBound);
    }
  else
    free(pdf);
  free(T);
  return mean;
}
/* -------------------------------------------------------------------------- */

double MDUCcalcLocalMeanVar(char *varName, double *mean, double *var,
                            double **pdf, int **pdfBound, int pdfFlag,void *mem)
{
  int i,j,k,r;
  mducStruct  *mduc   = (mducStruct*) mem;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;

  double beta1 = 1.0, beta2 = 1.0, beta3 = 1.0, meanTot, S, Q, tmp;
  double *v, vMean;

  if (!strncmp(varName, "heat-release", 8)) {
    v = flame->HRrate;
  }
  else if (!strncmp(varName, "density", 4)) {
    v = flame->rho;
  }
  else if (!strncmp(varName, "specific-heat", 10) || !strcmp(varName, "cp")) {
    v = flame->cpmix;
  }
  else {
    fprintf(stderr, "MDUCcalcLocalMean Warning: %s not recognized\n", varName);
    return -1.0;
  }
  
  if (flame->nDim == 1) {
    if (pdfFlag == 0) {
      /* ensure that the mean and variance are properly bounded */
      mean[0]  = fmin(1.0, fmax(0.0, mean[0]));
      var[0]   = fmax(0.0, var[0]);
      
      /* compute the shape parameters for the beta distribution */
      tmp   = mean[0]*(1-mean[0])/var[0] - 1;
      beta1 = mean[0]*tmp;
      beta2 = (1-mean[0])*tmp;
    }
    
    vMean = calcLinMeanY(0, 1, beta1, beta2, flame->n1, flame->x1,
                         v, *pdf, *pdfBound, pdfFlag);
  }
  else if (flame->nDim == 2) {
    int n[2];
    double **x;
    x = (double**) calloc(2, sizeof(double*));
    n[0] = flame->n1; n[1] = flame->n2;
    x[0] = flame->x1; x[1] = flame->x2;
    
    if (pdfFlag == 0) {
      /* ensure that the mean and variance are properly bounded */
      mean[0] = fmin(1.0, fmax(0.0, mean[0]));
      mean[1] = fmin(1.0, fmax(0.0, mean[1]));
      meanTot = mean[0] + mean[1];
      var[0]  = fmax(0.0, var[0]);
      var[1]  = fmax(0.0, var[1]);
      var[2]  = fmax(0.0, var[2]);
      
      /* compute mean and variance parameters */
      S = mean[0]*mean[0] + mean[1]*mean[1] + (1-meanTot)*(1-meanTot);
      Q = var[0] + var[1] + var[2];
      
      /* compute the shape parameters for the beta distribution */
      tmp   = (1.0 - S)/Q - 1.0;
      beta1 = mean[0]*tmp;
      beta2 = mean[1]*tmp;
      beta3 = (1.0-meanTot)*tmp;
    }

    vMean = calcLinMeanY_2D(0, 1, n, x, v, beta1, beta2, beta3, pdf,
                            pdfBound, pdfFlag, driver->TRANS_Z2);
    free(x);
  }

  return vMean;
}
/* -------------------------------------------------------------------------- */

double MDUCcalcMeanGamma(double **pdfIn, int **pdfBound, int pdfFlag,
                         double *Rmean, void *mem)
{
  int i, j, k, row, col, r;
  mducStruct   *mduc   = (mducStruct*) mem;
  chemStruct   *chem   = mduc->chem;
  flowStruct   *flow   = mduc->flow;
  flameStruct  *flame  = mduc->flame;
  engineStruct *engine = mduc->engine;

  double gamma, cpmix, wmix, R;
  double *pdf, *T;

  if (pdfIn == NULL) {
    pdf = flow->pdf;
    pdfBound = (int**) calloc(flame->n2, sizeof(int*));
    for (k=0; k<flame->n2; k++) {
      pdfBound[k] = (int*) calloc(2, sizeof(int));
      pdfBound[k][0] = 0;
      pdfBound[k][1] = flame->n1-1;
    }
  }
  else {
    pdf = (double*) calloc(flame->n1*flame->n2, sizeof(double));
    convertArrayMD_1D(flame->n1, flame->n2, 1, &pdfIn, pdf);
  }

  gamma = 0.0;
  if (Rmean) Rmean[0] = 0.0;
  for (k=0; k<flame->n2; k++) {
    row = k*flame->n1*flame->nVars;
    r   = k*flame->n1;
    for (j=0; j<flame->n1; j++) {
      col   = j*flame->nVars;
      
      /* compute average molecular weight */
      wmix  = calcMixtureW(chem->nSpecSolve, &flame->Y[row+col], chem->W);
      R     = flame->RGAS/wmix;
      
      if (Rmean) Rmean[0] += R*pdf[r+j];

      /* compute mixture average specific heat */
      chemComputeThermoData(chem->h,chem->cp,flame->Y[row+col+flame->sT],chem);
      cpmix = calcMixtureCp(chem->nSpecSolve, chem->cp, &flame->Y[row+col]);
      gamma += cpmix/(cpmix - R)*pdf[r+j];
    } /* end loop over columns */
  } /* end loop over rows */

  if (engine->spark && engine->ignited) {
    double mu = 0, mb = 0, m;
    for (i=0; i<3; i++) mu += flame->Y[flame->sM+i];
    mb = flame->Y[flame->sMb];
    m = mu + mb;
    gamma = mu/m*gamma + mb/m*engine->cpPrem/(engine->cpPrem - engine->Rprem);
  }

  if (pdfIn == NULL) {

  }
  else
    free(pdf);

  return gamma;
}
/* -------------------------------------------------------------------------- */

void MDUCsetSootMom(int n1, int nSootMom, double **sootMom, void *mem)
{
  int i, j;

  /* dereference internal flamelet solver structs */
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  flameStruct *flame  = mduc->flame;

  for (j=0; j<flame->n1; j++) {
    for (i=0; i<nSootMom; i++) {
      if (sootMom[j][i] > 0.0) 
        flame->sootMom[j*nSootMom+i] = sootMom[j][i];
    } /* end loop over moments */
  } /* end loop over grid */
}
/* -------------------------------------------------------------------------- */
void MDUCsetViscSrc(int n1,double *viscSrc, void *mem)
{
  int j;

  /* dereference internal flamelet solver structs */
  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;

  for (j=1; j<flame->n1-1; j++) {
        flame->visc_src[j] = viscSrc[j];
  } /* end loop over grid */

  /* set bounds to zero */ 
  flame->visc_src[0] = flame->visc_src[flame->n1-1] = 0.0;
}
void MDUCsetRhoDot(int n1,double *rhoDot, void *mem)
{
  int j;

  /* dereference internal flamelet solver structs */
  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;

  for (j=1; j<flame->n1-1; j++) {
        flame->rhoDotSpray[j] = rhoDot[j];
  } /* end loop over grid */

  /* set bounds to zero */ 
  flame->rhoDotSpray[0] = flame->rhoDotSpray[flame->n1-1] = 0.0;
}
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void MDUCsetGlobalPDF(int n1,double *gPDF, void *mem)
{
  int j;

  /* dereference internal flamelet solver structs */
  mducStruct  *mduc   = (mducStruct*) mem;
  driveStruct *driver = mduc->driver;
  flameStruct *flame  = mduc->flame;

  if (driver->output_w) {
     for (j=0; j<flame->n1; j++) {
           flame->gPDF[j] = gPDF[j];
     } /* end loop over grid */
  }
}
/* -------------------------------------------------------------------------- */

int MDUCcalcMeanSootSrc(double mean1, double var1, double mean2, double var2,
                        double var, double *meanSootSrc, void *mem)
{
  int i, j, k, row, col;
  int pdfFlag, mducY = 0;

  /* dereference internal flamelet solver structs */
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;

  /* compute the mean */
  if (flame->nDim == 1) {
    int pdfBound1D[2];
    double *pdf1D = (double*) calloc(flame->n1, sizeof(double));
    
    pdfFlag = MDUCcalcPDF_1D(PDF_TYPE_BETA,flame->n1,flame->x1,mean1,var1,
                             pdf1D,pdfBound1D);

    for (i=0; i<chem->nSootMom; i++) meanSootSrc[i] = 0.0;
    for (j=0; j<flame->n1; j++) {
      for (i=0; i<chem->nSootMom; i++) {
        meanSootSrc[i] += flame->sootMomSrc[j*chem->nSootMom+i]*pdf1D[j];
      } /* end loop over moments */
    } /* end loop over grid */
    free(pdf1D);
  }
  else if (flame->nDim == 2) {
    fprintf(stderr, "Warning: soot not implemented for 2D solution yet\n");
  }

  return pdfFlag;
}
/* -------------------------------------------------------------------------- */
void MDUCgetSootSrc(double **Ssrc, void *mem)
{
  int i, j, k, row, col;
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  flameStruct *flame  = mduc->flame;
  
  /* ----- copy flamelet solution to flow solver solution array ----- */

  for (j=0; j<flame->n1; j++) {
    for (i=0; i<chem->nSootMom; i++) {
      Ssrc[j][i] = flame->sootMomSrc[j*chem->nSootMom+i];
        } /* end loop over moments */
  } /* end loop over x1 */
}
void MDUCgetSootSolution(double **S, void *mem)
{
  int i, j, k, row, col;
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  flameStruct *flame  = mduc->flame;
  
  /* ----- copy flamelet solution to flow solver solution array ----- */

  for (j=0; j<flame->n1; j++) {
    for (i=0; i<chem->nSootMom; i++) {
      S[j][i] = flame->sootMom[j*chem->nSootMom+i];
        } /* end loop over moments */
  } /* end loop over x1 */
}
void MDUCgetRadiation1D(double *S, void *mem)
{
  int i, j, k, row, col;
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;
  driveStruct *driver = mduc->driver;
  
  /* ----- copy flamelet solution to flow solver solution array ----- */
  if (driver->radiation) 
    for (j=0; j<flame->n1; j++) 
      S[j] = flame->radiation[j];
  else
    for (j=0; j<flame->n1; j++) 
        S[j] = 0.0;
 
}
/* -------------------------------------------------------------------------- */

void MDUCsootHMOM(int nMom,double T_SI, double Visc_SI, double Rho_SI, double W_SI,
                  double *c_spec, double *mom_SI, double *src_SI, double *rxnrates,
                  double *rhodot_SI, double *srcP_SI, int *srcActive)
{
  double MolarMass;

  /*
    compute source terms for six phyical/chemical processes,taking into account 7 species
    units : T [K], Visc [Pa*s], Rho [kg/m^3], W[kg/mol], c_spec(7) [kmol/m^3], 
            mom(nMom), src(nMom), rxnrates(7) [kmol/m^3/s], rhodot(7) [kg/m^3/s], 
            srcP(6*nMom), srcActive(6)
    species need to be provided in following order:
    A4 (pyrene), OH, H2O, H, H2, C2H2, O2
    moments in order: M_00: number density (1/m^3); M_10: volume fraction (m^3/m^3);
    M_01: surface area density (m^2/m^3); [optional: M_20: volume variance (m^6/m^3); 
    M_11: volume-surface covariance (m^5/m^3); M_02: surface variance (m^4/m^3)]; 
    N_0: weight of the delta function (1/m^3)
    srcActive needs to be provided like this:
    Nucleation,Condensation,Coagulation,SurfaceGrowth,Oxidation,Fragmentation
   */

  MolarMass = W_SI*1000.0; // convert to g/mol

  sootHMOMsetOrder(nMom);

  sootHMOM(T_SI, Visc_SI, Rho_SI, MolarMass, c_spec, mom_SI, src_SI, rxnrates,
           rhodot_SI, srcP_SI, srcActive,1.0E-20);
}
/* -------------------------------------------------------------------------- */

int MDUCcalcMeanY_2D(int method, int *n, double **x, int nSpecies,
                     double mean1, double var1, double mean2, double var2,
                     double var, double ***Y, double *Ymean, double **pdf,
                     int **pdfBound, int transform)
{
  int    i, j, k, row, col, pdfFlag;
  double mean, beta1, beta2, beta3, S, Q, tmp;
  double *Yi;
  
  pdfFlag = MDUCcalcPDF_2D(method, n, x, mean1, var1, mean2, var2, var, pdf,
                           pdfBound, transform);

  if (pdfFlag == 0)
    {
      // ensure that the mean and variance are properly bounded
      mean1 = fmin(1.0, fmax(0.0, mean1));
      mean2 = fmin(1.0, fmax(0.0, mean2));
      mean  = mean1 + mean2;
      var1  = fmax(0.0, var1);
      var2  = fmax(0.0, var2);
      var   = fmax(0.0, var);
        
      // compute mean and variance parameters
      S = mean1*mean1 + mean2*mean2 + (1-mean)*(1-mean);
      Q = var1 + var2 + var;

      // compute the shape parameters for the beta distribution
      tmp   = (1.0 - S)/Q - 1.0;
      beta1 = mean1*tmp;
      beta2 = mean2*tmp;
      beta3 = (1.0-mean)*tmp;
    }

  Yi = (double*) calloc(n[0]*n[1]*nSpecies, sizeof(double));
  convertArrayMD_1D(nSpecies, n[0], n[1], Y, Yi);

  // loop over species
  for (i=0; i<nSpecies; i++)
    Ymean[i] = calcLinMeanY_2D(i, nSpecies, n, x, Yi, beta1, beta2, beta3,
                               pdf, pdfBound, pdfFlag, transform);

  free(Yi);
  return pdfFlag;
}
/* -------------------------------------------------------------------------- */

void MDUCupdateTempDPDT(double P, double dPdt, double **TIn, void *mem)
{
  int i, j, k, row, col, r;
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;
  timeStruct  *time   = mduc->time;
  double wmix, cpmix, *T;

  if (P < 0.0)   P = flow->P;
  
  T = (double*) calloc(flame->n1*flame->n2, sizeof(double));
  if (TIn == NULL)
    extractVarArray(flame->nVars,flame->n1,flame->n2,flame->sT,flame->Y,T);
  else
    convertArrayMD_1D(flame->n1, flame->n2, 1, &TIn, T);

  /* ----- compute the contribution from the center points ----- */
  for (k=0; k<flame->n2; k++)
    {
      row = k*flame->n1*flame->nVars; 
      r   = k*flame->n1;
      for (j=0; j<flame->n1; j++)
        {
          col = j*flame->nVars;
          
          // get enthalpy and specific heat from chemical mechanism
          chemComputeThermoData(chem->h, chem->cp, T[r+j], chem);
          
          // compute local mixture properties
          wmix  = calcMixtureW(chem->nSpecSolve, &flame->Y[row+col], chem->W);
          flame->rho[r+j] = wmix*P/(flame->RGAS*T[r+j]);
          cpmix = calcMixtureCp(chem->nSpecSolve, chem->cp, &flame->Y[row+col]);
          
          // reset the temperature to the value at t_n
          T[r+j]  = flame->Y[row+col+flame->sT];
          
          // add the pressure source term
          T[r+j] += dPdt*time->dt/(flame->rho[r+j]*cpmix);
        } // end loop over columns
    } // end loop over rows
  
  if (TIn != NULL)
    convertArray1D_MD(flame->n1, flame->n2, 1, T, &TIn);
  free(T);
}
/* -------------------------------------------------------------------------- */

#ifdef ACCESS_CFD
void MDUCcalcAvgSdot(int nCFD1, int nCFD2, int nCFD3, double ***Zmean,
                     double ***Zvar, double ***Hmean, double ***Hvar,
                     double ***rhoMean, void *mem)
{
  int ic, jc, kc, j, k;
  int row, col;
  int nZ, nH;
  double *Z, *H;
  double *pdfZ, *pdfH;
  double wmix, temp, rho;
  int pdfBoundZ[2], pdfBoundH[2];

  double *sumH;
  double **Sdot, *SdotNum, *SdotDen;

  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  timeStruct  *time   = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;

  if (driver->H_DIM == 1)
    {
      nH = flame->n1;
      nZ = flame->n2;
      H = flame->x1;
      Z = flame->x2;
    }
  else if (driver->H_DIM == 2)
    {
      nH = flame->n2;
      nZ = flame->n1;
      H = flame->x2;
      Z = flame->x1;
    }

  // assign local memory
  pdfZ    = (double*) calloc(nZ, sizeof(double));
  pdfH    = (double*) calloc(nH, sizeof(double));
  sumH    = (double*) calloc(nZ, sizeof(double));
  Sdot    = (double**) calloc(nZ, sizeof(double*));
  SdotNum = (double*) calloc(nZ, sizeof(double));
  SdotDen = (double*) calloc(nZ, sizeof(double));

  for (j=0; j<nZ; j++)
    {
      SdotNum[j] = 0.0;
      SdotDen[j] = 0.0;
      Sdot[j]    = (double*) calloc(nH, sizeof(double));
    }
 
  for (kc=0; kc<nCFD3; kc++)
    {
      for (jc=0; jc<nCFD2; jc++)
        {
          for (ic=0; ic<nCFD1; ic++)
            {
              // compute the local PDF values for mixture fraction and enthalpy
              calcBetaPDF(nZ, Z, Zmean[kc][jc][ic], Zvar[kc][jc][ic], pdfZ,
                          pdfBoundZ);
              calcBetaPDF(nH, H, Hmean[kc][jc][ic], Hvar[kc][jc][ic], pdfH,
                          pdfBoundH);

              // compute Sdot for all non-zero PDF values
              for (k=pdfBoundH[0]; k<=pdfBoundH[1]; k++)
                {
                  // NOTE: this will be dependent on the enthalpy direction
                  if (driver->H_DIM == 1)
                    col = k*flame->nVars;
                  else if (driver->H_DIM == 2)
                    row = k*flame->nVars*nZ;
            

                  for (j=pdfBoundZ[0]; j<=pdfBoundZ[1]; j++)
                    {

                      if (driver->H_DIM == 1)
                        row = j*flame->nVars*nH;
                      else if (driver->H_DIM == 2)
                        col = j*flame->nVars;

                      // compute local density and pressure source term
                      temp = flame->Y[row + col + flame->sT];
                      wmix = calcMixtureW(chem->nSpecies, &flame->Y[row+col],
                                          chem->W);
                      rho  = wmix*flow->P/(flame->RGAS*temp);
                      Sdot[k][j] = 1.0/rho*flow->Dt_P;
                    } // end loop over enthalpy
                } // end loop over mixture fraction

              // now integrate over enthalpy and sum
              for (j=pdfBoundZ[0]; j<=pdfBoundZ[1]; j++)
                {

                  sumH[j] = 0.0;
                  for (k=pdfBoundH[0]; k<=pdfBoundH[1]; k++)
                    sumH[j] += Sdot[k][j]*pdfH[k];

                  SdotNum[j] += rhoMean[kc][jc][ic]*sumH[j]*pdfZ[j];
                  SdotDen[j] += rhoMean[kc][jc][ic]*pdfZ[j];
                  // NOTE: pdf denominator should = 1, therefore not included
                }
            } // end loop over all CFD cells
        }
    }

  // compute final volume averaged source term
  for (j=0; j<nZ; j++)
    {
      if (SdotDen[j] > 1.0e-6)
        flow->avgSdot[j] = SdotNum[j]/SdotDen[j];
    }

  for (j=0; j<nZ; j++) free(Sdot[j]);
  free(pdfZ);
  free(pdfH);
  free(sumH);
  free(Sdot);
  free(SdotNum);
  free(SdotDen);
}
/* -------------------------------------------------------------------------- */
#endif // ACCESS_CFD

void MDUCwriteScalars(double curTime, double *mean, double *var, double **pdf,
                      int **pdfBound, int pdfFlag, void *mem)
{
  int i;
  
  // dereference internal flamelet solver structs
  mducStruct   *mduc   = (mducStruct*) mem;
  ioStruct     *io     = mduc->io;
  driveStruct  *driver = mduc->driver;
  chemStruct   *chem   = mduc->chem;
  timeStruct   *time   = mduc->time;
  flowStruct   *flow   = mduc->flow;
  flameStruct  *flame  = mduc->flame;
  engineStruct *engine = mduc->engine;
  
  double beta1 = 1.0, beta2 = 1.0, beta3 = 1.0, meanTot, S, Q, tmp;

  double Ysum = 0.0;
  double *Ymean = (double*) calloc(io->nOutputs, sizeof(double));

  int   fileExists;
  FILE *pFile;

  /* check whether the file needs to be initialized */
  if (io->initScalar) {
    /* check whether file already exists */
    fileExists = truncateFileTime(io->scalars, curTime);
    
    if (!fileExists)
      writeScalarsHeader(driver, io);
    io->initScalar = 0;
  }

  /* append the current timestep to the file */
  pFile = fopen(io->scalars, "a");
  if (flame->nDim == 1) {
    if (pdfFlag == 0) {
      /* ensure that the mean and variance are properly bounded */
      mean[0]  = fmin(1.0, fmax(0.0, mean[0]));
      var[0]   = fmax(0.0, var[0]);
      
      /* compute the shape parameters for the beta distribution */
      tmp   = mean[0]*(1-mean[0])/var[0] - 1;
      beta1 = mean[0]*tmp;
      beta2 = (1-mean[0])*tmp;
    }
    
    fprintf(pFile, "%15.6e ", curTime);
    if (driver->TIME_MODE == 3) /* engine */
      fprintf(pFile, "%15.6e ", time_to_cad(curTime, engine));
    for (i=0; i<io->nOutputs; i++) {
      Ymean[i] = calcLinMeanY(io->outputIndex[i],flame->nVars,beta1,beta2,
                              flame->n1,flame->x1,flame->Y,*pdf,*pdfBound,
                              pdfFlag);
      fprintf(pFile, "%15.6e ", Ymean[i]);
      Ysum += Ymean[i];
    }
    fprintf(pFile, "%15.6e ", Ysum);

    /* temperature */
    if (driver->ENERGY)
      fprintf(pFile, "%15.6e ", 
              calcLinMeanY(0, 1, beta1, beta2, flame->n1, flame->x1,
                           flame->T, *pdf, *pdfBound, pdfFlag));
    else
      fprintf(pFile, "%15.6e ", 
              calcLinMeanY(flame->sT, flame->nVars, beta1, beta2, flame->n1,
                           flame->x1, flame->Y, *pdf, *pdfBound, pdfFlag));
    
    /* pressure */
    if (driver->ENGINE)
      fprintf(pFile, "%15.6e ", flame->Y[flame->sP]);
    else
      fprintf(pFile, "%15.6e ", flow->P);
    
    if (io->outputHtot) {
      if (driver->ENERGY)
        fprintf(pFile, "%15.6e ",
                calcLinMeanY(0, 1, beta1, beta2, flame->n1, flame->x1, flame->T,
                             *pdf, *pdfBound, pdfFlag));
      else
        fprintf(pFile, "%15.6e ", 
                calcLinMeanY(0, 1, beta1, beta2, flame->n1, flame->x1, flame->H,
                             *pdf, *pdfBound, pdfFlag));
    }
    
    if (io->outputChi) {
      fprintf(pFile, "%15.6e ",
              calcLinMeanY(0, 1, beta1, beta2, flame->n1, flame->x1,
                           flow->chi1, *pdf, *pdfBound, pdfFlag));
    }
    if (io->outputHR) {
      fprintf(pFile, "%15.6e ",
              calcLinMeanY(0, 1, beta1, beta2, flame->n1, flame->x1,
                           flame->HRrate, *pdf, *pdfBound, pdfFlag));
    }
    if (io->outputRho) {
      fprintf(pFile, "%15.6e ",
              calcLinMeanY(0, 1, beta1, beta2, flame->n1, flame->x1,
                           flame->rho, *pdf, *pdfBound, pdfFlag));
    }
    if (io->outputCp) {
      fprintf(pFile, "%15.6e ",
              calcLinMeanY(0, 1, beta1, beta2, flame->n1, flame->x1,
                           flame->cpmix, *pdf, *pdfBound, pdfFlag));
    }
    if (io->outputQdot) {
      fprintf(pFile, "%15.6e ",
              calcLinMeanY(0, 1, beta1, beta2, flame->n1, flame->x1,
                           flame->srcQ, *pdf, *pdfBound, pdfFlag));
    }
    fprintf(pFile, "%15.6e %15.6e %15.6e %15.6e ",
            1.0/flame->rho[0]*flow->Dt_P,           flow->Dt_phi1min,
            1.0/flame->rho[flame->n1-1]*flow->Dt_P, flow->Dt_phi1max);
    fprintf(pFile, "%15.6e %15.6e %15.6e %15.6e %15.6e %15.6e ",
            flow->Dt_phiRaw[1][0].min, flow->Dt_phiRaw[1][0].max,
            flow->Dt_phiHat[1][0].min, flow->Dt_phiHat[1][0].max,
            flow->DtHat_phi[1][0].min, flow->DtHat_phi[1][0].max);
    
  }
  else if (flame->nDim == 2) {
    int n[2];
    double **x;
    x = (double**) calloc(2, sizeof(double*));
    x[0] = flame->x1;
    x[1] = flame->x2;
    n[0] = flame->n1;
    n[1] = flame->n2;
    
    if (pdfFlag == 0) {
      /* ensure that the mean and variance are properly bounded */
      mean[0] = fmin(1.0, fmax(0.0, mean[0]));
      mean[1] = fmin(1.0, fmax(0.0, mean[1]));
      meanTot = mean[0] + mean[1];
      var[0]  = fmax(0.0, var[0]);
      var[1]  = fmax(0.0, var[1]);
      var[2]  = fmax(0.0, var[2]);
      
      /* compute mean and variance parameters */
      S = mean[0]*mean[0] + mean[1]*mean[1] + (1-meanTot)*(1-meanTot);
      Q = var[0] + var[1] + var[2];
      
      /* compute the shape parameters for the beta distribution */
      tmp   = (1.0 - S)/Q - 1.0;
      beta1 = mean[0]*tmp;
      beta2 = mean[1]*tmp;
      beta3 = (1.0-meanTot)*tmp;
    }

    fprintf(pFile, "%15.6e ", curTime);
    if (driver->TIME_MODE == 3) /* engine */
      fprintf(pFile, "%15.6e ", time_to_cad(curTime, engine));
    for (i=0; i<io->nOutputs; i++) {
      Ymean[i] = calcLinMeanY_2D(io->outputIndex[i], flame->nVars, n, x,
                                 flame->Y,beta1,beta2,beta3,pdf,pdfBound,
                                 pdfFlag,driver->TRANS_Z2);
      fprintf(pFile, "%15.6e ", Ymean[i]);
      Ysum += Ymean[i];
    }
    fprintf(pFile, "%15.6e ", Ysum);
    fprintf(pFile, "%15.6e ",
            calcLinMeanY_2D(flame->sT, flame->nVars, n , x, flame->Y,
                            beta1, beta1, beta3, pdf, pdfBound, pdfFlag,
                            driver->TRANS_Z2));
    if (io->outputHtot)
      fprintf(pFile, "%15.6e ",
              calcLinMeanY_2D(0, 1, n , x, flame->H, beta1, beta2, beta3, pdf,
                              pdfBound, pdfFlag, driver->TRANS_Z2));
    if (io->outputChi) {
      fprintf(pFile, "%15.6e ",
              calcLinMeanY_2D(0, 1, n , x, flow->chi1, beta1, beta2, beta3, pdf,
                              pdfBound, pdfFlag, driver->TRANS_Z2));
    }  
    if (io->outputHR) {
      fprintf(pFile, "%15.6e ",
              calcLinMeanY_2D(0, 1, n , x, flame->HRrate, beta1, beta2, beta3,
                              pdf, pdfBound, pdfFlag, driver->TRANS_Z2));
    }
    if (io->outputRho) {
      fprintf(pFile, "%15.6e ",
              calcLinMeanY_2D(0, 1, n , x, flame->rho, beta1, beta2, beta3, pdf,
                              pdfBound, pdfFlag, driver->TRANS_Z2));
    }
    if (io->outputCp) {
      fprintf(pFile, "%15.6e ",
              calcLinMeanY_2D(0, 1, n , x, flame->cpmix, beta1, beta2, beta3,
                              pdf, pdfBound, pdfFlag, driver->TRANS_Z2));
    }
    free(x);
  }
  else if (flame->nDim == 0) {
    fprintf(pFile, "%15.6e ", curTime);
    if (driver->TIME_MODE == 3) /* engine */
      fprintf(pFile, "%15.6e ", time_to_cad(curTime, engine));
    for (i=0; i<io->nOutputs; i++) {
      fprintf(pFile, "%15.6e ", flame->Y[io->outputIndex[i]]);
      Ysum += flame->Y[io->outputIndex[i]];
    }
    fprintf(pFile, "%15.6e ", Ysum);
    fprintf(pFile, "%15.6e ", flame->Y[flame->sT]);
    if (driver->ENGINE)
      fprintf(pFile, "%15.6e ", flame->Y[flame->sP]);
    else
      fprintf(pFile, "%15.6e ", flow->P);
    if (io->outputHtot)
      fprintf(pFile, "%15.6e ", flame->H[0]);
    if (io->outputChi)
      fprintf(pFile, "%15.6e ", flow->chi1[0]);
    if (io->outputHR)
      fprintf(pFile, "%15.6e ", flame->HRrate[0]);
    if (io->outputRho)
      fprintf(pFile, "%15.6e ", flame->rho[0]);
    if (io->outputCp)
      fprintf(pFile, "%15.6e ", flame->cpmix[0]);
    if (io->outputQdot) {
      fprintf(pFile, "%15.6e ", flame->srcQ[0]);
      flame->intSrcQ += flame->srcQ[0]*time->dt;
      fprintf(pFile, "%15.6e ", flame->intSrcQ);
    }
    if (io->outputGamma)
      fprintf(pFile, "%15.6e ", flame->gamma[0]);
    if (driver->ENGINE) {
      fprintf(pFile, "%15.6e ", flame->Y[flame->sM]);
      fprintf(pFile, "%15.6e ", flame->Y[flame->sM+1]);
      fprintf(pFile, "%15.6e ", flame->Y[flame->sM+2]);
    }
  }
  
  /* calc dry NOX 6 percent */
  if (io->outputDryNOX) {
    int iH2O = getSpeciesIndex(chem->nSpecies, chem->species, "H2O");
    int iNO  = getSpeciesIndex(chem->nSpecies, chem->species, "NO");
    int iO2  = getSpeciesIndex(chem->nSpecies, chem->species, "O2");
    double NODry6 = 0.0, H2O, XNO_dry, XO2_dry;
    
    /* only do if all the species necessary are present in the mechanism */
    if ((iH2O >= 0) && (iNO >= 0) && (iO2 >= 0)) {
      if (flame->nDim == 1) {
        H2O     = calcLinMeanX(iH2O,flame->nVars,beta1,beta2,flame->n1,
                               flame->x1,flame->Y,*pdf,*pdfBound,
                               pdfFlag,chem->nSpecies,chem->W);
        XNO_dry = calcLinMeanX(iNO,flame->nVars,beta1,beta2,flame->n1,
                               flame->x1,flame->Y,*pdf,*pdfBound,pdfFlag,
                               chem->nSpecies,chem->W)/(1.0 - H2O);
        XO2_dry = calcLinMeanX(iO2,flame->nVars,beta1,beta2,flame->n1,
                               flame->x1,flame->Y,*pdf,*pdfBound,pdfFlag,
                               chem->nSpecies,chem->W)/(1.0 - H2O);
        NODry6 = 1.e6*XNO_dry*(4.76*(3+(1-XO2_dry)*2)/(1-4.76*XO2_dry)-2)/
          (4.76*(3+(1-0.06)*2)/(1-4.76*0.06)-2);
        fprintf(pFile, "%15.6e ", NODry6);
      }
    }
  }
  
  /* if comparing to DNS, compute the balance terms in the temperature eq. */
  if (driver->DNS && flame->nDim == 2) {
    int    j, k, jst, kst, r, row, col;
    double wmix, *Ypre, dx, dy;
    double Cp, Cp_m1, rho, rho_m1;           /* density and specific heat */
    double dTdt, Cp_dTdt;                    /* temperature rate term */
    double Dx_T, *Dy_T, Dxx_T, Dxy_T, Dyy_T; /* spatial derivatives */
    double diff1, diff12, diff2;             /* diffusion terms */
    
    /* make sure we have the correct solution at n-1 */
    if (flame->nDim == 1)
      Ypre = flame->Y_m1;
    else if (flame->nDim == 2)
      Ypre = flame->Y_pre;
    
    /* write the dPdt term */
    fprintf(pFile, "%15.6e ", (flow->P - flow->P_pre)/time->dt);
    
    /* store the gradient in y for computing the cross-derivative */
    Dy_T = (double*) calloc(flame->n1*flame->n2, sizeof(double));
    for (k=0; k<flame->n2; k++) {
      r   = k*flame->n1;
      kst = stencil_y(k, flame);
      for (j=0; j<flame->n1; j++)
        Dy_T[r+j] = grad_y(kst,k,j,flame->sT,flame->nVars,flame->Y,flame);
    }
    
    /* compute the average rho*Cp*dTdt */
    Cp_dTdt = diff1 = diff12 = diff2 = 0.0;
    for (k=0; k<flame->n2; k++) {
      row = k*flame->n1*flame->nVars;
      r   = k*flame->n1;
      kst = stencil_y(k, flame);
      if (k == 0)
        dy = 0.5*(flame->x2[k+1] - flame->x2[k]);
      else if (k == flame->n2-1)
        dy = 0.5*(flame->x2[k] - flame->x2[k-1]);
      else
        dy = 0.5*(flame->x2[k+1] - flame->x2[k-1]);
      
      for (j=0; j<flame->n1; j++) {
        col = j*flame->nVars;
        jst = stencil_x(j, flame);
        if (j == 0)
          dx = 0.5*(flame->x1[j+1] - flame->x1[j]);
        else if (j == flame->n1-1)
          dx = 0.5*(flame->x1[j] - flame->x1[j-1]);
        else
          dx = 0.5*(flame->x1[j+1] - flame->x1[j-1]);
        
        /* compute the average Cp at this point */
        chemComputeThermoData(chem->h,chem->cp,flame->Y[row+col+flame->sT],chem);
        Cp = calcMixtureCp(chem->nSpecSolve,chem->cp,&flame->Y[row+col]);
        chemComputeThermoData(chem->h, chem->cp, Ypre[row+col+flame->sT], chem);
        Cp_m1 = calcMixtureCp(chem->nSpecSolve,chem->cp,&Ypre[row+col]);
        
        /* calculate current and previous density */
        wmix   = calcMixtureW(chem->nSpecSolve,&flame->Y[row+col],chem->W);
        rho    = wmix*flow->P/(flame->RGAS*flame->Y[row+col+flame->sT]);
        wmix   = calcMixtureW(chem->nSpecSolve, &Ypre[row+col], chem->W);
        rho_m1 = wmix*flow->P_pre/(flame->RGAS*Ypre[row+col+flame->sT]);
        
        /* compute the rate of d/dt(rho*cp*T) */
        dTdt     = (rho*Cp*flame->Y[row+col+flame->sT] -
                    rho_m1*Cp_m1*Ypre[row+col+flame->sT])/time->dt;
        Cp_dTdt += dTdt*pdf[k][j];
        
        /* compute the gradients and the diffusion terms */
        Dx_T  = grad_x (jst,k,j,flame->sT,flame->nVars,flame->Y,flame);
        Dxx_T = grad_xx(jst,k,j,flame->sT,flame->nVars,flame->Y,flame);
        Dxy_T = grad_x (jst,k,j,0,1,Dy_T,flame);
        Dyy_T = grad_yy(kst,k,j,flame->sT,flame->nVars,flame->Y,flame);
        
        diff1 += 0.5*flow->chi1[r+j]*Dxx_T*pdf[k][j]*dx*dy;
        diff12 +=    flow->chi12[r+j]*Dxy_T*pdf[k][j]*dx*dy;
        diff2  += 0.5*flow->chi2[r+j]*Dyy_T*pdf[k][j]*dx*dy;
      } /* end loop over columns */
    } /* end loop over rows */
    
    /* write terms to file */
    fprintf(pFile,"%15.6e %15.6e %15.6e %15.6e ",Cp_dTdt,diff1,diff12,diff2);
    free(Dy_T);
  }

  fprintf(pFile, "\n");
  fclose(pFile);
  free(Ymean);
}
/* -------------------------------------------------------------------------- */

double MDUCupdateOxTemp(double H, double guess, void *mem)
{
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  flameStruct *flame  = mduc->flame;
  
  flame->Y[flame->sT] = calcTfromH(chem->nSpecies, chem->nSpecSolve, H,
                                   flame->Y, guess, chem);

  return flame->Y[flame->sT];
}
/* -------------------------------------------------------------------------- */

double MDUCcalcTfromH(int nSolve, double H, double *Y, double guess, void *mem)
{
  int nSpecies = MDUCgetNspecies(mem);
  mducStruct *mduc = (mducStruct*) mem;
  return calcTfromH(nSpecies, nSolve, H, Y, guess, mduc->chem);
}
/* -------------------------------------------------------------------------- */

double MDUCcalcHtot(int iRow, int iCol, double T, void *mem)
{
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  flameStruct *flame  = mduc->flame;

  int i,row, col;
  double H;

  // reference indices for the specified row and column
  row = iRow*flame->nVars*flame->n1;
  col = iCol*flame->nVars;

  // get data from mechanism
  chemComputeThermoData(chem->h, chem->cp, T, chem);
  
  // sum over species to get total enthalpy
  H = 0.0;
  for (i=0; i<chem->nSpecies; i++)
    H += chem->h[i]*flame->Y[row+col+i];

  return H;
}
/* -------------------------------------------------------------------------- */

double MDUCreadFlowData(double **chi1, double **chi12, double **chi2,
                        double **pdf, double *sdot, double timeval, void *mem)
{
  int j, k, r;
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  timeStruct  *time   = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;

  char filestring[200], temp[40];
  int timeint = (int) (timeval);
  double *pChi1, *pChi2, *pChi12, tFile;

  sprintf(temp, "%0*i", io->fileNumPrec, timeint);
  strcpy(filestring, io->postFile);
  strcat(filestring, temp);
  strcat(filestring, ".flo");

  // Read flow file
  tFile = readFlowData(filestring, driver, flow);

  // assign correct time level to assign to based on option
  if (driver->CHI_MODE == -1)
    {
      pChi1  = flow->chi1;
      pChi12 = flow->chi12;
      pChi2  = flow->chi2;
    }
   else
    {
      pChi1  = flow->chi1_m1;
      pChi12 = flow->chi12_m1;
      pChi2  = flow->chi2_m1;
    }

  // copy scalar dissipation and pdf arrays if passed
  for (k=0; k<flame->n2; k++)
    {
      r = k*flame->n1;
      for (j=0; j<flame->n1; j++)
        {
          if (pChi1 != NULL && chi1 != NULL)
            chi1[k][j] = pChi1[r+j];
          if (pChi12 != NULL && chi12 != NULL)
            chi12[k][j] = pChi12[r+j];
          if (pChi2 != NULL && chi2 != NULL)
            chi2[k][j] = pChi2[r+j];
          if (flow->pdf != NULL && pdf != NULL)
            pdf[k][j]  = flow->pdf[r+j];
        } // end loop over columns
    } // end loop over rows

  if (driver->H_DIM == 1)
    {
      if (flow->avgSdot != NULL && sdot != NULL)
        copyArray(flame->n2, flow->avgSdot, sdot);
    }
  else if (driver->H_DIM == 2)
    {
      if (flow->avgSdot != NULL && sdot != NULL)
        copyArray(flame->n1, flow->avgSdot, sdot);
    }

  return tFile;
}
/* -------------------------------------------------------------------------- */

void MDUCwriteFlowData(double **chi1, double **chi2, double **pdf,
                       double *sdot, void *mem)
{
  int j, k, r;
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  timeStruct  *time   = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;

  char outputFileName[256];

  // make sure memory is allocated
  flow->chi1 = (double*) realloc(flow->chi1,flame->n1*flame->n2*sizeof(double));
  flow->chi2 = (double*) realloc(flow->chi2,flame->n1*flame->n2*sizeof(double));
  flow->pdf  = (double*) realloc(flow->pdf, flame->n1*flame->n2*sizeof(double));

  for (k=0; k<flame->n2; k++)
    {
      r = k*flame->n1;
      for (j=0; j<flame->n1; j++)
        {
          if (chi1 != NULL) flow->chi1[r+j] = chi1[k][j];
          if (chi2 != NULL) flow->chi2[r+j] = chi2[k][j];
          if (pdf  != NULL) flow->pdf[r+j]  = pdf[k][j];
        }
    }

  if (sdot != NULL) 
    {
      if (driver->H_DIM == 1)
        {
          for (j=0; j<flame->n2; j++)
            flow->avgSdot[j] = sdot[j];
        }
      if (driver->H_DIM == 2)
        {
          for (j=0; j<flame->n1; j++)
            flow->avgSdot[j] = sdot[j];
        }
    }

  // create data file using same style as flamelet file
  createDataFileName(io->flowData, io, time, outputFileName);

  writeFlowData(outputFileName, io, driver, time, flame, flow);
}
/* -------------------------------------------------------------------------- */

void MDUCwriteFlowDataMem(char *fileName, void *mem)
{
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  timeStruct  *time   = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;

  if (fileName == NULL)
    {
      char outputFileName[256];
      createDataFileName(io->flowData, io, time, outputFileName);
      writeFlowData(outputFileName, io, driver, time, flame, flow);
    }
  else
    writeFlowData(fileName, io, driver, time, flame, flow);
}
/* -------------------------------------------------------------------------- */

void MDUCwriteVariableToFlamelet(char *varName, double **var, void *mem)
{
  // dereference flamelet structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  timeStruct  *time   = mduc->time;
  flameStruct *flame  = mduc->flame;

  int j, k;
  char outputFileName[256];
  double *tmp = (double*) calloc(flame->n1*flame->n2, sizeof(double));

  // convert to 1D array
  for (k=0; k<flame->n2; k++)
    for (j=0; j<flame->n1; j++)
      tmp[k*flame->n1 + j] = var[k][j];

  // re-create current data-file name
  createDataFileName(io->outputData, io, time, outputFileName);

  addVarToBody(io->inputData, flame->n1, flame->n2, tmp, varName);

  free(tmp);
}
/* -------------------------------------------------------------------------- */

void MDUCreadVariable(char *fileName,char *varName,int n1,int n2,double **var)
{
  int j, k;

  double *varTmp = (double*) calloc(n1*n2, sizeof(double));

  readVariableArray(fileName, varName, &varTmp);

  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      var[k][j] = varTmp[k+j*n1];

  free(varTmp);
}
/* -------------------------------------------------------------------------- */

void MDUCtransformVar(int s, int trans, int n1, int n2, double *x1, double *x2,
                      double **var)
{
  int j, k, iZ1, iZ2;
  double Z1, Z2;
  double interp;

  // transform 2D unstructured vector according to transformation rules
  double **tmp = (double**) calloc(n2, sizeof(double*));
  for (k=0; k<n2; k++) tmp[k] = (double*) calloc(n1, sizeof(double));

  // copy original variable to temporary array
  if (s == 0) // unstructured array
    convertArray1D_MD(n1, n2, 1, var[0], &tmp);
  else // structured array
    for (k=0; k<n2; k++) copyArray(n1, var[k], tmp[k]);

  // set values outside of realizable domain to zero
  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      if (x1[j] + x2[k] > 1.0)
        tmp[k][j] = 0.0;
  
  if (trans == 1) // standard
    {
      for (k=0; k<n2; k++)
        {
          for (j=0; j<n1; j++)
            {
              // only need to interpolate in the eta direction
              Z2  = x2[k]*(1.0-x1[j]);
              iZ2 = gridLocate(n2, x2, Z2);

              // compute the interpolated value
              interp=linInterp(Z2,x2[iZ2],tmp[iZ2][j],x2[iZ2+1],tmp[iZ2+1][j]);

              // assign to either appropriate array type
              if (s == 0) // unstructured
                var[0][k*n1+j] = interp;
              else // structured
                var[k][j] = interp;
            }
        }
    }
  else if (trans == 2) // Hasse
    {
      for (k=0; k<n2; k++)
        {
          for (j=0; j<n1; j++)
            {
              // need to do bilinear interpolation
              Z1  = x1[j]*(1.0-x2[k]);
              Z2  = x1[j]*x2[k];
              iZ1 = gridLocate(n1, x1, Z1);
              iZ2 = gridLocate(n2, x2, Z2);

              // compute the interpolated value
              interp = biLinInterp(Z1,Z2,x1[iZ1],x1[iZ1+1],x2[iZ2],x2[iZ2+1],
                                   tmp[iZ2][iZ1],   tmp[iZ2][iZ1+1],
                                   tmp[iZ2+1][iZ1], tmp[iZ2+1][iZ1+1]);

              // assign to either appropriate array type
              if (s == 0) // unstructured
                var[0][k*n1+j] = interp;
              else // structured
                var[k][j] = interp;
            }
        }
    }
  
  for (k=0; k<n2; k++) free(tmp[k]);
  free(tmp);
}
/* -------------------------------------------------------------------------- */

void MDUCfitChi1D(int n1, double xmin, double xmax, double *x1, double *chi1)
{
  int    i, j, k, jDat, jLast, jFirst, nDat, nPoly, *ia;
  double *xDat, *chiDat, *sig, *pCoeff, **covar, chisq;
  double tmp, TINY = 1.0e-6;

  /* find the number of non-zero chi elements (as well as last value) */
  for (nDat=0, jLast=0, j=0; j<n1; j++)
    if (chi1[j] > TINY)
      {
        nDat++;
        jLast = j;
      }

  /* check to see if the grid is at the maximum extents */
  if ((xmax - x1[n1-1]) > 1.0e-3) /* normalised grid, use final point as zero */
    jLast = n1-1;
  else
    jLast++;

  nDat += 2;  // add space for zero end points
  /* allocate data arrays */
  xDat   = (double*) calloc(nDat, sizeof(double));
  chiDat = (double*) calloc(nDat, sizeof(double));
  sig    = (double*) calloc(nDat, sizeof(double));
  for (j=0; j<nDat; j++) sig[j] = 1.0; // no info know about standard deviation

  /* copy the non-zero values into a data array */
  xDat[0] = chiDat[0] = 0.0;
  for (jDat=1, j=1; j<n1; j++)
    if (chi1[j] > TINY)
      {
        xDat[jDat]   = x1[j];
        chiDat[jDat] = chi1[j];
        jDat++;
      }
  xDat[nDat-1]   = x1[jLast];
  chiDat[nDat-1] = 0.0;

  /* base the order of the polynomial on the number of data points available,
     with minimum 4th order and maximum 10th order.  Otherwise, order is half
     the number of data points */
  nPoly = IMAX(5, IMIN(11, nDat/2));

  /* allocate memory for polynomial least squares fit */
  pCoeff = (double*)  calloc(nPoly, sizeof(double));
  covar  = (double**) calloc(nPoly, sizeof(double*));
  ia     = (int*)     calloc(nPoly, sizeof(int));
  for (i=0; i<nPoly; i++)
    {
      covar[i] = (double*) calloc(nPoly, sizeof(double));
      ia[i] = 1;
    }

  /* call the linear least squares regression, assuming polynomial fit */
  lfit(xDat, chiDat, sig, nDat, pCoeff, ia, nPoly, covar, &chisq, fitPoly);

  /* calculate the function using the polynomial fit */
  chi1[0] = chi1[jLast-1] = 0.0;
  for (j=1; j<jLast-1; j++)
    {
      chi1[j] = 0.0;
      for (i=0; i<nPoly; i++)
        {
          tmp = 1.0;
          for (k=0; k<i; k++) tmp *= x1[j];
          chi1[j] += pCoeff[i]*tmp;
        }
      chi1[j] = fmax(chi1[j], 0.0); // ensure no negative values
    }
  
  free(xDat);
  free(chiDat);
  free(sig);
  free(pCoeff);
  free(ia);
  for (i=0; i<nPoly; i++) free(covar[i]);
  free(covar);
}
/* -------------------------------------------------------------------------- */

void MDUClaplaceChi1D(int n1,double xmin,double xmax,double *x1,double *chi1)
{
  // construct the matrix for the 1D
  int i, j, jmax, *indx, *mask, iter, interp = 0;
  double *b, *y, v1 = 0.0, err;

  mask = (int*) calloc(n1, sizeof(int));
  b    = (double*) calloc(n1, sizeof(double));
  y    = (double*) calloc(n1, sizeof(double));

  int norm = 0;
  if ((xmax - x1[n1-1]) > 1.0e-3) norm = 1;

  jmax = n1-1;
  if (!norm)
    {
      for (j=jmax; j>0; j--)
        {
          if (chi1[j] > 0.0)
            {
              jmax = j;
              break;
            }
        }
    }

  // enforce boundary conditions
  chi1[0] = chi1[n1-1] = 0.0;

  // set up the system of equations
  for (j=0; j<n1; j++)
    {
      if (chi1[j] > 0.0 || j == 0 || j == n1-1 || j > jmax )
        {
          b[j] = y[j] = v1 = chi1[j];
          mask[j] = 1;
        }
      else
        {
          b[j] = 0.0;
          y[j] = v1;
          mask[j] = 0;
          interp = 1; // interpolation need for a point
        }
    }
  
  if (!interp)
    {
      free(y);
      free(b);
      free(mask);
      return;
    }

  LapIntStruct *LapArgs = (LapIntStruct*) calloc(1, sizeof(LapIntStruct));
  LapArgs->n1 = n1;
  LapArgs->x1 = x1;
  LapArgs->mask = mask;

  Linbcg(n1, b, y, 1, 1.0e-6, 2*n1, &iter, &err, LapArgs,
         LaplaceInterp_atimes_1D, LaplaceInterp_asolve);

  for (j=0; j<n1; j++) chi1[j] = y[j];
  
  free(mask);
  free(b);
  free(y);
  free(LapArgs);
}
/* -------------------------------------------------------------------------- */

void MDUCinterpChiLap(double **chi1, double **chi12, double **chi2, void *mem)
{
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;

  // construct the matrix for the 1D
  int i, j, k, r, iter, interp[3] = {0,0,0};
  double v1, v2, v3, err;
  int n1  = flame->n1;
  int n2  = flame->n2;
  int nEq = n1*n2;

  int *jmin  = (int*) calloc(n2, sizeof(int));
  int *jmax  = (int*) calloc(n2, sizeof(int));
  int *kmin  = (int*) calloc(n1, sizeof(int));
  int *kmax  = (int*) calloc(n1, sizeof(int));
  int *mask  = (int*) calloc(nEq, sizeof(int));
  int *mask2 = (int*) calloc(nEq, sizeof(int));
  int *mask3 = (int*) calloc(nEq, sizeof(int));
  double *b  = (double*) calloc(nEq, sizeof(double));
  double *y  = (double*) calloc(nEq, sizeof(double));
  double *b2 = (double*) calloc(nEq, sizeof(double));
  double *y2 = (double*) calloc(nEq, sizeof(double));
  double *b3 = (double*) calloc(nEq, sizeof(double));
  double *y3 = (double*) calloc(nEq, sizeof(double));

  double **chi3;
  // compute chi3 from the raw data
  chi3 = (double**) calloc(flame->n2, sizeof(double*));
  for (k=0; k<flame->n2; k++)
    {
      chi3[k] = (double*) calloc(flame->n1, sizeof(double));
      for (j=0; j<flame->n1; j++)
        chi3[k][j] = chi1[k][j] + 2*chi12[k][j] + chi2[k][j];
    }

  // handle normalized grids
  int norm1 = 0;
  if ((flow->phi1max - flame->x1[n1-1]) > 1.0e-3) norm1 = 1;

  for (k=0; k<n2; k++) // set default to edges
    {
      jmin[k] = 0;
      jmax[k] = n1 - 1;
    }
  if (!norm1) // find range of non zero points if non-normalized grid
    {
      for (k=0; k<n2; k++)
        {
          do { jmin[k]++; } while (chi1[k][jmin[k]] <= 0.0 && jmin[k] < n1-1);
          do { jmax[k]--; } while (chi1[k][jmax[k]] <= 0.0 && jmax[k] > 0);
        }
    }

  int norm2 = 0;
  if ((flow->phi2max - flame->x2[n2-1]) > 1.0e-3) norm2 = 1;
  for (j=0; j<n1; j++) // set default to edges
    {
      kmin[j] = 0;
      kmax[j] = n2 - 1;
    }
  if (!norm2) // find range of non zero points if non-normalized grid
    {
      for (j=0; j<n1; j++)
        {
          do { kmin[j]++; } while (chi2[kmin[j]][j] <= 0.0 && kmin[j] < n2-1);
          do { kmax[j]--; } while (chi2[kmax[j]][j] <= 0.0 && kmax[j] > 0);
        }
    }
  
  // set up the system(s) of equations
  v1 = v2 = v3 = 0.0;
  for (k=0; k<n2; k++)
    {
      r = k*n1;
      for (j=0; j<n1; j++)
        {
          // for chi1
          if (chi1[k][j] > 0.0 || j < jmin[k] || j > jmax[k]
              || k < kmin[j] || k > kmax[j])
            {
              b[r+j] = y[r+j] = v1 = chi1[k][j];
              mask[r+j] = 1;
            }
          else
            {
              b[r+j] = 0.0;
              y[r+j] = v1;
              mask[r+j] = 0;
              interp[0] = 1; // interpolation need for a point
            }

          // for chi2
          if (chi2[k][j] > 0.0 || j < jmin[k] || j > jmax[k]
              || k < kmin[j] || k > kmax[j] )
            {
              b2[r+j] = y2[r+j] = v2 = chi2[k][j];
              mask2[r+j] = 1;
            }
          else
            {
              b2[r+j] = 0.0;
              y2[r+j] = v2;
              mask2[r+j] = 0;
              interp[1] = 1; // interpolation need for a point
            }

          // for chi3
          if (chi3[k][j] > 0.0 || j < jmin[k] || j > jmax[k]
              || k < kmin[j] || k > kmax[j] )
            {
              b3[r+j] = y3[r+j] = v3 = chi3[k][j];
              mask3[r+j] = 1;
            }
          else
            {
              b3[r+j] = 0.0;
              y3[r+j] = v3;
              mask3[r+j] = 0;
              interp[2] = 1; // interpolation need for a point
            }
        }
    }

  // set up data needed by Laplace interpolator  
  LapIntStruct *LapArgs = (LapIntStruct*) calloc(1, sizeof(LapIntStruct));
  LapArgs->n1 = n1;
  LapArgs->n2 = n1;
  LapArgs->x1 = flame->x1;
  LapArgs->x2 = flame->x2;
  LapArgs->mask = mask;

  // Laplace interpolate each chi using biconjugate gradient method
  if (interp[0])
    Linbcg(nEq, b, y, 1, 1.0e-6, 2*IMAX(n1,n2), &iter, &err, LapArgs,
           LaplaceInterp_atimes_2D, LaplaceInterp_asolve);

  LapArgs->mask = mask2;
  if (interp[1])
    Linbcg(nEq, b2, y2, 1, 1.0e-6, 2*IMAX(n1,n2), &iter, &err, LapArgs,
           LaplaceInterp_atimes_2D, LaplaceInterp_asolve);

  LapArgs->mask = mask3;
  if (interp[2])
    Linbcg(nEq, b3, y3, 1, 1.0e-6, 2*IMAX(n1,n2), &iter, &err, LapArgs,
           LaplaceInterp_atimes_2D, LaplaceInterp_asolve);

  // copy the final solutions back
  for (k=0; k<n2; k++)
    {
      for (j=0; j<n1; j++)
        {
          if (interp[0])
            chi1[k][j] = y[k*n1+j];
          if (interp[1])
            chi2[k][j] = y2[k*n1+j];
          if (interp[2])
            chi3[k][j] = y3[k*n1+j];

          // compute missing chi12 from interpolated values
          if (!mask[k*n1+j])
            chi12[k][j] = 0.5*(chi3[k][j] - chi1[k][j] - chi2[k][j]);
        }
    }
  
  for (k=0; k<n2; k++) free(chi3[k]);
  free(chi3);
  free(mask); free(mask2); free(mask3);
  free(b); free(b2); free(b3); free(y); free(y2); free(y3);
  free(jmin); free(jmax); free(kmin); free(kmax);
  free(LapArgs);
}
/* -------------------------------------------------------------------------- */

void MDUCinterpChiDNS(double **chi1, double **chi12, double **chi2, void *mem)
{
  // dereference internal flamelet solver structs
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;

  int    n1, n2;
  int    j, k, jL, jR, jMax, kA, kB, norm;
  double xi, eta, chi3Local;
  double **chi3, **chi3H, **chi3Int; // scalar dissipation of "third" scalar
  double **chi2T;                    // transposed arrays

  // compute chi3 from the raw data
  if (flame->nDim > 1 && chi12 != NULL && chi2 != NULL)
    {
      chi3 = (double**) calloc(flame->n2, sizeof(double*));
      for (k=0; k<flame->n2; k++)
        {
          chi3[k] = (double*) calloc(flame->n1, sizeof(double));
          for (j=0; j<flame->n1; j++)
            chi3[k][j] = chi1[k][j] + 2*chi12[k][j] + chi2[k][j];
        }
    }

  /* ----- chi1 interpolation ----- */
  // check to see if x1 grid is at the maximum extents
  norm = 0;
  if ((flow->phi1max - flame->x1[flame->n1-1]) > 1.0e-3) norm = 1;

  // interpolate chi1 along x1 for each row
  for (k=0; k<flame->n2; k++)
    interpDNS_1D(flame->n1, flame->x1, norm, flow->phi1max, chi1[k]);

  if (flame->nDim == 1 || chi2 == NULL) // finished if only 1D flamelet
    return;

  /* ----- chi2 interpolation ----- */
  // compute transposed chi2
  chi2T = (double**) calloc(flame->n1, sizeof(double*));
  for (j=0; j<flame->n1; j++)
    {
      chi2T[j] = (double*) calloc(flame->n2, sizeof(double));
      for (k=0; k<flame->n2; k++)
        chi2T[j][k] = chi2[k][j];
    }

  // check to see if x2 grid is at the maximum extents
  norm = 0;
  if ((flow->phi2max - flame->x2[flame->n2-1]) > 1.0e-3) norm = 1;
  // interpolate chi2 along x2 for each column
  for (j=0; j<flame->n1; j++)
    interpDNS_1D(flame->n2, flame->x2, norm, flow->phi2max, chi2T[j]);

  // retranspose chi2
  for (k=0; k<flame->n2; k++)
    for (j=0; j<flame->n1; j++)
      chi2[k][j] = chi2T[j][k];

  if (driver->H_DIM == 0)
    {
      // allocate and fill chi3 transformed and interpolated arrays
      chi3H   = (double**) calloc(flame->n2, sizeof(double*));
      chi3Int = (double**) calloc(flame->n2, sizeof(double*));
      for (k=0; k<flame->n2; k++)
        {
          chi3H[k]   = (double*) calloc(flame->n1, sizeof(double));
          chi3Int[k] = (double*) calloc(flame->n1, sizeof(double));
          for (j=0; j<flame->n1; j++)
            chi3H[k][j] = chi3Int[k][j] = chi3[k][j];
        }
      // interpolate onto transformed grid with eta = x1+x2
      MDUCtransformVar(1,2,flame->n1,flame->n2,flame->x1,flame->x2,chi3H);
      // make sure that the left hand column is constant
      for (k=0; k<flame->n2; k++)
        chi3H[k][0] = 0.0;

      /* interpolate chi3 along Z1 direction for now.  Really should be
         interpolated normal to lines of constant Z1+Z2 */
      for (k=0; k<flame->n2; k++)
        interpDNS_1D(flame->n1, flame->x1, 0, 1, chi3Int[k]);

      // just compute missing chi12 from interpolated chis for now
      for (k=0; k<flame->n2; k++)
        {
          for (j=0; j<flame->n1; j++)
            {
              if (chi3[k][j] <= 0.0  && (chi1[k][j] > 0.0 || chi2[k][j] > 0.0))
                chi12[k][j] = 0.5*( fmax(0.0,chi3Int[k][j])
                                    - chi2[k][j] - chi1[k][j]);
            }
        }

      // enforce boundary conditions
      for (k=0; k<flame->n2; k++)
        {
          chi1[k][0]  = chi1[k][flame->n1-1]  = 0.0;
          chi2[k][flame->n1-1] = 0.0;
          chi12[k][0] = chi12[k][flame->n1-1] = 0.0;
        }
      for (j=0; j<flame->n1; j++) chi2[0][j] = chi12[0][j] = 0.0;

      // free memory
      for (k=0; k<flame->n2; k++)
        {
          free(chi3H[k]);
          free(chi3Int[k]);
        }
      free(chi3H);
      free(chi3Int);
    } // end two mixture fraction
  else
    {
      /* Interpolation of cross-dissipation of Z-H not yet implemented, for now
         it is usually assumed to be negligible.  Potentially should be simple
         bilinear or Laplace interpolation to fill missing data */
    }

  for (k=0; k<flame->n2; k++) free(chi3[k]);
  free(chi3);
  for (j=0; j<flame->n1; j++) free(chi2T[j]);
  free(chi2T);
}
/* -------------------------------------------------------------------------- */

double MDUCgetHRrate(double volume, void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;
  /* TODO: generalize to all dimensions, maybe pass pdf as well */
  return volume*flame->HRrate[0]; /* J/s */
}
/* -------------------------------------------------------------------------- */

double MDUCgetHeatLoss(void *mem)
{
  mducStruct  *mduc   = (mducStruct*) mem;
  flameStruct *flame  = mduc->flame;
  /* TODO: generalize to all dimensions, maybe pass pdf as well */
  return flame->srcQ[0]; /* J/s */
}
/* -------------------------------------------------------------------------- */

void MDUCcalcHRrate(void *mem) /* TODO: this function does not belong here */
{
  int i, j, k, row, col, r;
  mducStruct  *mduc   = (mducStruct*) mem;
  chemStruct  *chem   = mduc->chem;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;
  timeStruct  *time   = mduc->time;
  double rho, wmix, temp, hr_local;
  

  double *Y_pre;
  if (flame->nDim == 2)
    Y_pre = flame->Y_pre;
  else
    Y_pre = flame->Y_m1;

  /* ----- compute and store the heat release ----- */
  for (k=0; k<flame->n2; k++) {
    row = k*flame->n1*flame->nVars;
    r   = k*flame->n1;
    for (j=0; j<flame->n1; j++) {
      col  = j*flame->nVars;
      // calculate mixture properties from previous time step
      wmix = calcMixtureW(chem->nSpecSolve,&Y_pre[row+col],chem->W);
      temp = Y_pre[row+col+flame->sT];
      rho  = wmix*flow->P_m1/(flame->RGAS*temp);
      chemComputeThermoData(chem->h, chem->cp, temp, chem);
      
      // sum contribution from all species
      hr_local = 0.0;
      for (i=0; i<chem->nSpecSolve; i++)
        hr_local += flame->chemSrcY[row+col+i]*rho*chem->h[i];
      flame->HRrate[r+j] = -hr_local/time->dt;
    }
  }
}
/* -------------------------------------------------------------------------- */

double MDUCcalcChiLogForm(double Z, double Zref, double Zmax)
{
  return calcLogForm(Z, Zref, Zmax);
}
/* -------------------------------------------------------------------------- */

double MDUCcalcChiErfcForm(double Z, double Zref, double Zmax)
{
  return calcErfcForm(Z, Zref, Zmax);
}
/* -------------------------------------------------------------------------- */

double MDUCcalcChiQuadForm(double Z, double Zref, double Zmax)
{
  return calcQuadForm(Z, Zref, Zmax);
}
/* -------------------------------------------------------------------------- */

double MDUCerf(double x) { return erf(x); }

#ifdef ACCESS_CFD
int MDUCcheckInjection(int nCFD1, int nCFD2, int nCFD3, double threshold,
                       double ***Z1mean, double ***Z2mean)
{
  int i, j, k;

  // CFD variables
  double ***Z1, ***Z2; // mean mixture fractions

  // set the pointers appropriately for the dimension of the array
  // NOTE: 1D and 2D passed as address of pointer
  if (nCFD3 == 1 && nCFD2 == 1)
    {
      Z1 = &Z1mean;
      Z2 = &Z2mean;
    }
  else
    {
      Z1 = Z1mean;
      Z2 = Z2mean;
    }

  // loop over all cells to check whether mixture fraction is present
  // look for the second injection
  for (k=0; k<nCFD3; k++)
    {
      for (j=0; j<nCFD2; j++)
        {
          for (i=0; i<nCFD1; i++)
            {
              if (Z2[k][j][i] > threshold)
                return 2;
            }
        }
    }

  // look for the first injection
  for (k=0; k<nCFD3; k++)
    {
      for (j=0; j<nCFD2; j++)
        {
          for (i=0; i<nCFD1; i++)
            {
              if (Z2[k][j][i] > threshold)
                return 1;
            }
        }
    }

  // return no injection if no mixture fraction found
  return 0;
}
/* -------------------------------------------------------------------------- */

void MDUCupdateCFD_1D(int nCFD1, int nCFD2, int nCFD3, int nYFav,
                      int *indexYFav, double ***C_Y, double ***C_T,
                      double ***C_H, double ***C_PHIMEAN, double ***C_PHIVAR,
                      int dim, int n1, double *x1, double **Y, double *T)
{
  // CFD variables
  double ****YFav;   // mean species mass fraction tracked in CFD
  double ***TFav;    // mean temperature
  double ***HFav;    // mean total enthalpy
  double ***phiFav;  // mean co-ordinate
  double ***phiVar;  // co-ordinate variance

  // set the pointers appropriately for the dimension of the array
  // NOTE: 1D and 2D passed as address of pointer
  if (nCFD3 == 1 && nCFD2 == 1)
    {
      YFav   = &C_Y;
      TFav   = &C_T;
      HFav   = &C_H;
      phiFav = &C_PHIMEAN;
      phiVar = &C_PHIVAR;
    }
  else
    {
      YFav   = C_Y;
      TFav   = C_T;
      HFav   = C_H;
      phiFav = C_PHIMEAN;
      phiVar = C_PHIVAR;
    }

  // local variables
  int    i, j, k, m, p, q, nZ1 = 0;
  int    counter, nCount = 500, lineCount = 200;
  double *pdf, *Ymean;
  int    pdfBound[2];
  double sumPDF;

  double tmp;
  double TINY = 1e-6;

  int nSpecies = GetNSpecies(); /* consider passing as argument */
  int nSpecSolve = GetNSpecs();

  // allocate all local memory
  pdf   = (double*) calloc(n1, sizeof(double));
  Ymean = (double*) calloc(nSpecies, sizeof(double));

  counter = 0;
  for (k=0; k<nCFD3; k++)
    {
      for (j=0; j<nCFD2; j++)
        {
          for (i=0; i<nCFD1; i++)
            {
              /* ----- compute mean mixture fraction of all species  -----  */
              if ( (phiFav[k][j][i] > 0.0) && (phiVar[k][j][i] > TINY))
                {
                  // compute the PDF for this cell
                  MDUCcalcPDF_1D(PDF_TYPE_BETA,n1,x1,phiFav[k][j][i],phiVar[k][j][i],
                                  pdf, pdfBound);

                  // error checking on the PDF
                  sumPDF = 0.0;
                  for (p=0; p<n1; p++)
                    {
                      if ( (1.0-pdf[p] < -1.0e-9) || (pdf[p] < -1.0e-9) )
                        {
                          printf("PDF computation out of bounds!\n");
                          printf("Z1mean = %15.6e\n", phiFav[k][j][i]);
                          printf("Z1var  = %15.6e\n", phiVar[k][j][i]);
                          printf("x pdf\n");
                          for (p=0; p<n1; p++)
                            printf("%8.6f %15.6e\n", x1[p], pdf[p]);

                          printf("PDF error!\n");
                          exit(2);
                        }
                      sumPDF += pdf[p];
                    }

                  /* compute mean mass fractions of all species by integrating
                     over the grid with the PDF */
                  for (q=0; q<nSpecies; q++) Ymean[q] = 0.0;
                  for(p=pdfBound[0]; p<=pdfBound[1]; p++)
                    {
                      for (q=0; q<nSpecies; q++)
                        Ymean[q] += Y[p][q]*pdf[p];
                    } /* end loop over columns */

                  nZ1++; /* count the number of cells computed */
                } /* end cell computation for non-zero Z1 */
              else
                {
                  /* no mixture fraction, update based on the oxidizer species,
                     essentially delta function at pure oxidizer */
                  for (q=0; q < nSpecies; q++)
                    Ymean[q] = Y[0][q];
                } // finished computing mean mixture fraction

              /* ----- update the cell properties ----- */
              /* update the mixture fractions tracked by Fluent */
              for (m=0; m<nYFav; m++)
                {
                  if (indexYFav[m] >= 0)
                    YFav[k][j][i][m] = Ymean[indexYFav[m]];
                }
              /* update the cell temperature by iterating from the total
                 enthalpy and mass fractions */
              TFav[k][j][i] = calcTfromH(nSpecies, nSpecSolve, HFav[k][j][i],
                                         Ymean,1000.0, chem);

              counter++;
              if (counter%nCount == 0)
                printf(".");
              if (counter%(lineCount*nCount) == 0)
                printf("\n");
            }
        }
    } // end loop over CFD cells
  printf("\n%i cells with Z1 updated (percent = %6.2f)\n",
         nZ1, (double)(nZ1)/counter*100.0);

  free(pdf);
  free(Ymean);
}
/* -------------------------------------------------------------------------- */
#endif // ACCESS_CFD
