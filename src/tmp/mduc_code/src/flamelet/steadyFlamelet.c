#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "steadyFlamelet.h"
#include "chemMech.h"
#include "accessChem.h"
#include "solve1D.h"
#include "flamelet1D.h"
#include "flamelet2D.h"
#include "fileIO.h"
#include "scalarDiss.h"
#include "utils.h"
#include "solveUtils.h"
#include "display.h"

void solveSteadyFlamelet1D(void *mem)
{
  /* dereference internal flamelet solver structs */
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  timeStruct  *time   = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;
  sdialStruct *sdial  = mduc->sdial;

  int i,j, flag;
  char outFileName[512];

  /* initialize the solver */
  initSteadyFlamelet(flame, sdial);
  writeFlamelet("init_steady_flamelet.rif", io, driver, time, flame, flow);

  /* set up an I/O */
  checkMonitor(io);

  /* ensure that all unsteady terms are zero */
  flow->Dt_P = 0.0;
  for (j=0; j<flame->n1; j++) {
    flame->v1[j] = 0.0;
  }

  /* make sure time is gone */
  time->t   = time->t_m1 = 0.0;
  time->t_n = 1.0;

  /* default to solving between 0 and 1 */
  flow->phi1min = 0.0;
  flow->phi1max = 1.0;
  flow->phi     = flow->phiRaw;
  flow->Dt_phi  = flow->Dt_phiRaw;
  flow->D_phi1  = flow->phiStar[1].del = 1.0;

  /* calculate the scalar dissipation rate */
  /* TODO: generalize the form of the dissipation rate here with options */
  /* calcChi1DcounterflowDiffFlame(flame->n1,flame->x1,0,strain,chi); */
  calcChi1D(io->chiMode, flame->n1, flame->x1, flow->ZRef[1], 1.0,
            flow->chiRef[1], flow->chi1_m1);

  /* make sure to setup the transformed scalar dissipation rate */
  flow->D_phi1_m1 = flow->phi1max_m1 - flow->phi1min_m1;
  flameletPrestep1D(driver, flame, chem, flow, io, time);

  /* call the KINSOL solver */
  fprintf(stderr, "Solving steady flamelet\n");
  flag = 0;
  flag = KINSol(sdial->mem,sdial->Y,KIN_LINESEARCH,sdial->Yscale,sdial->Fscale);

  /* use a basic newton solver */
  /* newt(sdial->Ydata, sdial->nEq, &flag, sdial, newtFlamelet1D, NULL); */

  if (flag < 0)
    fprintf(stderr, "KINSOL Error: %i\n", flag);
  double fNorm;
  KINGetFuncNorm(mem, &fNorm);
  if (flag == KIN_SUCCESS)
    fprintf(stderr, "KINSol successfully converged (%i): Fnorm = %15.6e\n",
            flag, fNorm);
  else 
    fprintf(stderr, "KINSol failed: %i\n", flag);

  steadyFlameletFileName(io->outputData, io, flow, outFileName);
  writeFlamelet(outFileName, io, driver, time, flame, flow);
}
/* -------------------------------------------------------------------------- */

void initSteadyFlamelet(flameStruct *flame, sdialStruct *kin)
{
  int i,j,k,flag;

  /* set the total number of equations to solve */
  kin->nEq = (flame->n1-2)*flame->n2*flame->nVars;
  
  /* allocate the N_Vector space needed */
  kin->Y      = N_VNew_Serial(kin->nEq);
  kin->Yscale = N_VNew_Serial(kin->nEq);
  kin->Fscale = N_VNew_Serial(kin->nEq);

  /* set pointers to data in N_Vector elements */
  kin->Ydata      = NV_DATA_S(kin->Y);
  kin->YscaleData = NV_DATA_S(kin->Yscale);
  kin->FscaleData = NV_DATA_S(kin->Fscale);
  
  if (flame->nDim == 1)
    kinFlamelet = kinFlamelet1D;
  else if (flame->nDim == 2)
    kinFlamelet = kinFlamelet2D;

  /* create and allocate the KINSOL memory object */
  kin->mem = KINCreate();
  /* set pointers to data and function calls */
#if SUNDIAL_2_3
  flag = KINMalloc(kin->mem, kinFlamelet, kin->Y);
  flag = KINSetFdata(kin->mem, kin);
#else
  flag = KINInit(kin->mem, kinFlamelet, kin->Y);
  flag = KINSetUserData(kin->mem, kin);
  /* flag = KINSetNoInitSetup(kin->mem, TRUE); */
  flag = KINSetNumMaxIters(kin->mem, 1000);
#endif /* SUNDIAL_2_3 */
  flag = KINSetPrintLevel(kin->mem, 3);
  flag = KINSetMaxSetupCalls(kin->mem, 30);

  /* set the solution scaling vectors (not at all sure what to do here) */
  for (i=0; i<kin->nEq; i++)
    kin->YscaleData[i] = 1.0;

  /* set the function scaling vectors */
  for (i=0; i<kin->nEq; i++)
    kin->FscaleData[i] = 10.0;

  for (k=0; k<flame->n2; k++) {
    for (j=0; j<flame->n1; j++) {
      kin->YscaleData[(k*flame->n1+j)*flame->nVars+flame->sT] = 0.001;
      kin->FscaleData[(k*flame->n1+j)*flame->nVars+flame->sT] = 0.000001;
    }
  }
  
  /* set the solver method */
  if (kin->iMethod == 0) { /* dense */
    fprintf(stderr, "Warning: using dense solver ");
    fprintf(stderr, " (this is inefficient in most cases)\n");
#ifdef SUNDIALS_FROM_3_0
    fprintf(stderr, "Extend to SUNDIALS_FROM_3_0\n");
#else
    flag = KINDense(kin->mem, kin->nEq);
#endif /* SUNDIALS_FROM_3_0 */
  }
  else if (kin->iMethod == 1) { /* banded */
    fprintf(stderr, "Kinsolve set to banded\n");
    fprintf(stderr, "  number of equations:   %i\n", kin->nEq);
    fprintf(stderr, "  upper/lower bandwidth: %i\n", flame->nVars);
#ifdef SUNDIALS_FROM_3_0
    fprintf(stderr, "Extend to SUNDIALS_FROM_3_0\n");
#else
    flag = KINBand(kin->mem, kin->nEq, flame->nVars, flame->nVars);
#endif /* SUNDIALS_FROM_3_0 */
  }
  else if (kin->iMethod == 2) { /* GMRES */
    fprintf(stderr, "Kinsolve set to GMRES\n");
#ifdef SUNDIALS_FROM_3_0
    fprintf(stderr, "Extend to SUNDIALS_FROM_3_0\n");
#else
    flag = KINSpgmr(kin->mem,0);
#endif /* SUNDIALS_FROM_3_0 */
  }
  else if (kin->iMethod == 3) { /* BICGSTAB */
    fprintf(stderr, "Kinsolve set to BICGSTAB\n");
#ifdef SUNDIALS_FROM_3_0
    fprintf(stderr, "Extend to SUNDIALS_FROM_3_0\n");
#else
    flag = KINSpbcg(kin->mem,0);
#endif /* SUNDIALS_FROM_3_0 */
  }
  else if (kin->iMethod == 4) { /* TFQMR */
    fprintf(stderr, "Kinsolve set to TFQMR\n");
#ifdef SUNDIALS_FROM_3_0
    fprintf(stderr, "Extend to SUNDIALS_FROM_3_0\n");
#else
    flag = KINSptfqmr(kin->mem,0);
#endif /* SUNDIALS_FROM_3_0 */
  }

  /* TODO: should initialize a preconditioner or Jacobian here */
  if (kin->iMethod > 1) { /* initialize precondition for Krylov methods */
    fprintf(stderr, "setting up preconditioner\n");
    flag = KINBBDPrecInit(kin->mem,kin->nEq,flame->nVars,flame->nVars,
                          flame->nVars,flame->nVars,0.0,kinFlamelet1D_Loc,NULL);
  }

  for (i=0; i<kin->nEq; i++)
    kin->Ydata[i] = flame->Y[flame->nVars+i];
}
/* -------------------------------------------------------------------------- */
#ifdef SUNDIALS_FROM_2_7
int kinFlamelet1D_Loc(long int Nlocal, N_Vector u, N_Vector fval, void *data)
#else
#ifdef SUNDIALS_FROM_3_0
int kinFlamelet1D_Loc(sunindextype Nlocal, N_Vector u, N_Vector fval, void *data)
#else /* SUNDIALS LOWER THAN 2.7 */
int kinFlamelet1D_Loc(int Nlocal, N_Vector u, N_Vector fval, void *data)
#endif /* SUNDIALS_FROM_3_0 */
#endif /* SUNDIALS_FROM_2_7 */
{
  return kinFlamelet1D(u, fval, data);
}
/* -------------------------------------------------------------------------- */

int kinFlamlet1D_comm(int Nlocal, N_Vector u, void *data)
{
  return 0;
}
/* -------------------------------------------------------------------------- */

int kinFlamelet1D(N_Vector u, N_Vector fval, void *data)
{
  int i,j,k,iX1;
  realtype *Y, *FVAL;
  double   *temp;

  /* dereference pointers to structs */
  sdialStruct  *kin    = (sdialStruct*)  data;
  driveStruct  *driver = (driveStruct*)  kin->args[0];
  chemStruct   *chem   = (chemStruct*)   kin->args[1];
  flowStruct   *flow   = (flowStruct*)   kin->args[2];
  flameStruct  *flame  = (flameStruct*)  kin->args[3];
  timeStruct   *time   = (timeStruct*)   kin->args[4];
  ioStruct     *io     = (ioStruct*)     kin->args[6];

  /* assign pointers to beginning of KINSOL vector memory */
  Y    = NV_DATA_S(u);
  FVAL = NV_DATA_S(fval);

  double *Ytmp = (double*) calloc(flame->n1*flame->nVars, sizeof(double));
  double *Ftmp = (double*) calloc(flame->n1*flame->nVars, sizeof(double));

  /* set the position of the last grid point */
  iX1 = (flame->n1-1)*flame->nVars;

  /* store the constant boundaries in the temporary array */
  for (i=0; i<flame->nVars; i++) {
    Ytmp[i]     = flame->Y[i];
    Ytmp[iX1+i] = flame->Y[iX1+i];
  }
  /* copy the interior points from the solution vector */
  for (j=1; j<flame->n1-1; j++) {
    for (i=0; i<flame->nVars; i++) {
      Ytmp[j*flame->nVars+i] = Y[(j-1)*flame->nVars+i];
    }
  }

  /* get the flamelet function */
  flameletRHS_1D(0.0, Ytmp, Ftmp, data);

  /* transfer the interior points to the solution array */
  for (i=0; i<kin->nEq; i++)
    FVAL[i] = Ftmp[flame->nVars+i];

  /* update the flamelet solution if writing out iteration solutions */
  /* TODO: this should only be done once at the end of the solution */
  for (i=0; i<flame->n1*flame->nVars; i++)
    flame->Y[i] = Ytmp[i];


/*   double srcT, minSrcT = 1e20, maxSrcT = -1e20, minT = 1e20, maxT = -1e20; */
/*   for (j=0; j<flame->n1; j++) { */
/*     srcT = FVAL[j*flame->nVars+flame->sT]; */
/*     if (srcT < minSrcT) minSrcT = srcT; */
/*     if (srcT > maxSrcT) maxSrcT = srcT; */
/*     if (Y[j*flame->nVars+flame->sT] < minT) minT = Y[j*flame->nVars+flame->sT]; */
/*     if (Y[j*flame->nVars+flame->sT] > maxT) maxT = Y[j*flame->nVars+flame->sT]; */
/*   } */
  /* fprintf(stderr, "min/max srcT: %8.2f/%8.2f (%8.2f/%8.2f)\n", */
/*           minSrcT, maxSrcT, minT, maxT); */
  long int nNonLin;
  KINGetNumNonlinSolvIters(kin->mem, &nNonLin);

  char *flameFileName = (char*) calloc(256, sizeof(char));

  /* write the u and F vectors */
  sprintf(flameFileName, "uF_vectors_it%li.dat", nNonLin);
  FILE *pF = fopen(flameFileName, "w");
  fprintf(pF, "# no u F\n");
  for (i=0; i<kin->nEq; i++)
    fprintf(pF, "%5i %15.6e %15.6e\n", i, Y[i], FVAL[i]);
  fclose(pF);

  /* fprintf(stderr, " boundary temperature: %f\n", flame->Y[iX1+flame->sT]); */
  sprintf(flameFileName, "flamelet_it%li.rif", nNonLin);
  writeFlamelet(flameFileName, io, driver, time, flame, flow);
  free(flameFileName);

  steadyFlameletMonitor(kin->mem);
  free(Ytmp); free(Ftmp);
  return 0;
}
/* -------------------------------------------------------------------------- */

void newtFlamelet1D(int nF, double *a, double *f, void *args)
{
  /* dereference pointers to structs */
  sdialStruct  *newt   = (sdialStruct*)  args;
  flameStruct  *flame  = (flameStruct*)  newt->args[3];

  int i,j,iX1;
  double *Ytmp = (double*) calloc(flame->n1*flame->nVars, sizeof(double));
  double *Ftmp = (double*) calloc(flame->n1*flame->nVars, sizeof(double));

  /* set the position of the last grid point */
  iX1 = (flame->n1-1)*flame->nVars;

  /* store the constant boundaries in the temporary array */
  for (i=0; i<flame->nVars; i++) {
    Ytmp[i]     = flame->Y[i];
    Ytmp[iX1+i] = flame->Y[iX1+i];
  }
  /* copy the interior points from the solution vector */
  for (j=1; j<flame->n1-1; j++) {
    for (i=0; i<flame->nVars; i++) {
      Ytmp[j*flame->nVars+i] = a[(j-1)*flame->nVars+i];
    }
  }

  /* get the flamelet function */
  flameletRHS_1D(0.0, Ytmp, Ftmp, args);

  /* transfer the interior points to the solution array */
  for (i=0; i<nF; i++)
    f[i] = Ftmp[flame->nVars+i];
  
  /* update the flamelet solution if writing out iteration solutions */
  /* TODO: this should only be done once at the end of the solution */
  for (i=0; i<flame->n1*flame->nVars; i++)
    flame->Y[i] = Ytmp[i];

  free(Ytmp); free(Ftmp);
}
/* -------------------------------------------------------------------------- */

void solveSteadyFlamelet2D(void *mem)
{
  /* dereference internal flamelet solver structs */
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  timeStruct  *timeS  = mduc->time;
  flowStruct  *flow   = mduc->flow;
  flameStruct *flame  = mduc->flame;
  sdialStruct *sdial  = mduc->sdial;

  int i,j,k,r,flag;
  char outFileName[512];

  /* initialize the solver */
  initSteadyFlamelet(flame, sdial);

  /* set up an I/O */
  checkMonitor(io);

  /* ensure that all unsteady terms are zero */
  flow->Dt_P = 0.0;
  for (k=0; k<flame->n2; k++) {
    r = k*flame->n1;
    for (j=0; j<flame->n1; j++) {
      flame->v1[r+j] = flame->v2[r+j] = 0.0;
    }
  }

  /* default to solving between 0 and 1 */
  flow->phi1min = flow->phi2min = 0.0;
  flow->phi1max = flow->phi2max = 1.0;
  flow->D_phi1  = flow->phiStar[1].del = 1.0;
  flow->D_phi2  = flow->phiStar[2].del = 1.0;

  /* calculate the scalar dissipation rate */
  /* TODO: generalize the form of the dissipation rate here with options */
  /* calcChi1DcounterflowDiffFlame(flame->n1,flame->x1,0,strain,chi); */
  /* calcChi1D(io->chiMode, flame->n1, flame->x1, flow->ZRef[1], 1.0, */
  /*           flow->chiRef[1], flow->chi1); */

  time_t startRun, endRun;
  double elapsedTime;
  time(&startRun); // start recording run time

  /* call the KINSOL solver */
  fprintf(stderr, "Solving steady 2D flamelet\n");
  /* write the header for the monitor */
  fprintf(stderr, "# iter fEval jacEval norm\n");

  flag = KINSol(sdial->mem, sdial->Y, KIN_NONE, sdial->Yscale,
                sdial->Fscale);
  /* if (flag < 0) */
  double fNorm;
  KINGetFuncNorm(mem, &fNorm);
  if (flag == KIN_SUCCESS)
    fprintf(stderr, "KINSol cuccessfully converged (%i): Fnorm = %15.6e\n",
            flag, fNorm);
  else 
    fprintf(stderr, "KINSol failed: %i\n", flag);

  time(&endRun);
  elapsedTime = difftime(endRun, startRun);
  fprintf(stderr, "Steady flamelet run time: %.2f\n", elapsedTime); 

  steadyFlameletFileName(io->outputData, io, flow, outFileName);
  writeFlamelet(outFileName, io, driver, timeS, flame, flow);
}
/* -------------------------------------------------------------------------- */

int kinFlamelet2D(N_Vector u, N_Vector fval, void *data)
{
  int i,j,k,r,row,col;
  realtype *Y, *FVAL;
  double   *temp;

  /* dereference pointers to structs */
  sdialStruct  *kin    = (sdialStruct*)  data;
  driveStruct  *driver = (driveStruct*)  kin->args[0];
  chemStruct   *chem   = (chemStruct*)   kin->args[1];
  flowStruct   *flow   = (flowStruct*)   kin->args[2];
  flameStruct  *flame  = (flameStruct*)  kin->args[3];
  timeStruct   *time   = (timeStruct*)   kin->args[4];

  /* assign pointers to beginning of KINSOL vector memory */
  Y    = NV_DATA_S(u);
  FVAL = NV_DATA_S(fval);

  /* ----- update quantities before computing derivatives ----- */
  if (driver->ENERGY) {
    /* update to most current temperature and enthalpy */
    for (k=0; k<flame->n2; k++) {
      r   = k*flame->n1;
      row = r*flame->nVars;
      for (j=0; j<flame->n1; j++) {
        col = j*flame->nVars;
        flame->T[r+j] = calcTfromH(chem->nSpecies, chem->nSpecSolve,
                                   Y[row+col+flame->sT], &Y[row+col],1000,chem);
      } /* end loop over columns */
    } /* end loop over rows */
    temp = flame->T;
  }
  else { /* temperature is part of solution */
    temp = (double*) calloc(flame->n1*flame->n2, sizeof(double));
    extractVarArray(flame->nVars, flame->n1, flame->n2, flame->sT, Y, temp);
  }

  /* ----- compute species derivatives ----- */
  flameletY_RHS_2D(driver, chem, flame, flow, temp, Y, FVAL);

  /* ----- compute energy equation ----- */
  if (driver->ENERGY) {
    /* update the energy based on the species current species and enthalpy */
    for (k=0; k<flame->n2; k++) {
      r   = k*flame->n1;
      row = r*flame->nVars;
      for (j=0; j<flame->n1; j++) {
        col = j*flame->nVars;
        chemComputeThermoData(chem->h, chem->cp, flame->T[r+j], chem);
        FVAL[row+col+flame->sT] = 0.0;
        for (i=0; i<chem->nSpecies; i++) {
          FVAL[row+col+flame->sT] += Y[row+col+i]*chem->h[i];
        } /* end loop over species */
      } /* end loop over columns */
    } /* end loop over rows */
  }
  else {
    flameletT_RHS_2D(driver, chem, flame, flow, Y, FVAL);
    free(temp);
  }
  steadyFlameletMonitor(kin->mem);
  return 0;
}
/* -------------------------------------------------------------------------- */

void steadyFlameletMonitor(void *mem)
{
  /* summarize the solver */
#ifdef SUNDIALS_FROM_2_7
  long int lsflag;
#else
  int lsflag;
#endif /* SUNDIALS_FROM_2_7 */
  long int nFeval, nNonLin, nIter, nbacktr, nbcfails, nJacEval = 0;
  double fNorm, step;
  /* KINSpilsGetNumFuncEvals(mem, &nFeval); */
  /* KINSpilsGetNumLinIters(mem, &nIter); */
  KINGetNumNonlinSolvIters(mem, &nNonLin);
  KINGetFuncNorm(mem, &fNorm);
  KINGetStepLength(mem, &step);
  /* KINGetNumJacEvals(sdial->mem, &nJacEval); */
  /* KINGetNumFuncEvals(mem, &nFeval); */
  KINDlsGetNumFuncEvals(mem, &nFeval);
  KINDlsGetNumJacEvals(mem, &nJacEval);
  KINDlsGetLastFlag(mem, &lsflag);
  KINGetNumBetaCondFails(mem, &nbcfails);
  KINGetNumBacktrackOps(mem, &nbacktr);
  fprintf(stderr, "%5li %8li %8li %5i %8li %8li %15.6e %15.6e\n",
          nNonLin, nFeval, nJacEval, lsflag, nbcfails, nbacktr, step, fNorm);
  /* fprintf(stderr, "Function Norm: %15.6e\n", funcNorm); */
  /* fprintf(stderr, "Function evaluations: %li\n", nFeval); */
  /* fprintf(stderr, "Jacobian evaluations: %li\n", nJacEval); */
  /* fprintf(stderr, "Non-linear iterations: %li\n", nIter); */

}
/* -------------------------------------------------------------------------- */
