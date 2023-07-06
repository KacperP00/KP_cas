#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "flameletInit.h"
#include "flamelet2D.h"
#include "accessMDUC.h"
#include "solve2D.h"
#include "chemMech.h"
#include "transport.h"
#include "accessChem.h"
#include "metric.h"
#include "solveUtils.h"
#include "utils.h"
#include "fileIO.h"
#include "display.h"

void solveFlamelet2D(double dt,double P,double phi1min,double phi1max,
                     double phi2min,double phi2max,double **chi1,double **chi12,
                     double **chi2,int inert,void *mem)
{
  int myid = 0;

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
#ifdef ACCESS_MPI
  mpiStruct *mpi = mduc->mpi; /* dereference mpi info */
  myid = mpi->myid;
#endif /* ACCESS_MPI */

  int i,j,k,r,c,row,col,jMaxT;
  char outputFileName[256];
  
  /* set CVODE variables */
  int      flag;        /* flags for CVODE return */
  realtype tEnd, dtPre; /* final time to integrate to and previous time step */
  long int nst;         /* current time step number */

  /* check whether Sundials needs to be initialized */
  if (sdial->mem == NULL) {
    if (driver->SDIAL == 0)
      cvodeInitFlamelet(time->t, flame->Y, sdial);
    else if (driver->SDIAL == 1)
      idaInitFlamelet(time->t, flame->Y, sdial);
  }

  checkMonitor(io);

  /* update time step to take */
  time->dt = dt;
  tEnd = time->t_n = time->t + dt;
  
  /* check to make sure the flow field has been initialised for the flamelet */
  if (flow->P_m1 == 0.0)
    mducErr("initial pressure has not been set","solveFlamelet1D");
  
  /* assign flow variables to struct */
  flow->P = P;
  flow->phi1min = phi1min;
  flow->phi1max = phi1max;
  flow->phi2min = phi2min;
  flow->phi2max = phi2max;

  /* compute pressure rate */
  flow->Dt_P = (flow->P - flow->P_m1)/time->dt;

  /* TODO: check if this makes sense */
  if (driver->H_DIM == 0) {
    flow->phi2min = phi2min/(1.0 - phi1min);
    if (phi1max < 1.0)
      flow->phi2max = phi2max/(1.0 - phi1max);
    else
      flow->phi2max = 1.0;
  }
 
  /* get normalisation for variables */
  flow->D_phi1     = flow->phi1max - flow->phi1min;
  flow->D_phi1_m1  = flow->phi1max_m1 - flow->phi1min_m1;
  flow->D_phi2     = flow->phi2max - flow->phi2min;
  flow->D_phi2_m1  = flow->phi2max_m1 - flow->phi2min_m1;

  if (driver->H_DIM == 1 && driver->DNS == 1) {
    flow->Dt_phi1min = (flow->phi1min_corr - flow->phi1min_m1)/time->dt;
    flow->Dt_phi1max = (flow->phi1max_corr - flow->phi1max_m1)/time->dt;
  }
  else {
    flow->Dt_phi1min = (flow->phi1min - flow->phi1min_m1)/time->dt;
    flow->Dt_phi1max = (flow->phi1max - flow->phi1max_m1)/time->dt;
  }
  if (driver->H_DIM == 2 && driver->DNS == 1) {
    flow->Dt_phi2min = (flow->phi2min_corr - flow->phi2min_m1)/time->dt;
    flow->Dt_phi2max = (flow->phi2max_corr - flow->phi2max_m1)/time->dt;
  }
  else {
    flow->Dt_phi2min = (flow->phi2min - flow->phi2min_m1)/time->dt;
    flow->Dt_phi2max = (flow->phi2max - flow->phi2max_m1)/time->dt;
  }

  /* clip */
  if (flow->Dt_phi2min<0.0) flow->Dt_phi2min = 0.0;
  if (flow->Dt_phi1min<0.0) flow->Dt_phi1min = 0.0;
  if (flow->Dt_phi2max>0.0) flow->Dt_phi2max = 0.0;
  if (flow->Dt_phi1max>0.0) flow->Dt_phi1max = 0.0;

  /* write the monitor header if first time step*/
  if (time->nt == time->nt0 && myid == 0) {
    fprintf(io->monitor, "# nst<1> time[s]<2> dt[s]<3> Tmax[K]<4>\n");
    fflush(io->monitor);
  }

  /* the pre-step does the following:
     - store solution of Y(N) to Y(N-1)
     - computes transformation variables, v1 and chiA
     (these are invariant over the time step)
     - compute initial RHS */
  solvePrestep2D(chi1,chi12,chi2,driver,flame,flow);

  if (inert) { /* ----- update based on new boundary temperature ----- */
    int    k_pre;
    int    iX1  = (flame->n1-1)*flame->nVars;           /* x1 = 1, x2 = 0 */
    int    iX2  = (flame->n2-1)*flame->n1*flame->nVars; /* x1 = 0, x2 = 1 */
    double *hox = (double*) calloc(chem->nSpecies, sizeof(double));
    double *hf1 = (double*) calloc(chem->nSpecies, sizeof(double));
    double *hf2 = (double*) calloc(chem->nSpecies, sizeof(double));

    /* calculate the thermo data at the boundaries */
    chemComputeThermoData(hox,chem->cp,flame->Y[flame->sT],chem);
    chemComputeThermoData(hf1,chem->cp,flame->Y[iX1 + flame->sT],chem);
    chemComputeThermoData(hf2,chem->cp,flame->Y[iX2 + flame->sT],chem);

    int rX2 = (flame->n2-1)*flame->n1;
    flame->H[0] = flame->H[flame->n1-1] = 0.0;
    flame->H[rX2] = 0.0;
    for (i=0; i<chem->nSpecies; i++) {
      flame->H[0]           += hox[i]*flame->Y[i];
      flame->H[flame->n1-1] += hf1[i]*flame->Y[iX1+i];
      flame->H[rX2]         += hf2[i]*flame->Y[iX2+i];
    }
    free(hox); free(hf1); free(hf2);
    
    /* interpolate the updated state from new boundary conditions */
    for (k=flame->j0; k<=flame->jN; k++) {
      
      /* set the left boundary */

      /* set the right boundary */
      

      /* update the interior points */
      for (j=1; j<flame->n1-1; j++) {

        /* TODO: finish inert 2D */

      } /* end loop over columns */
    } /* end loop over rows */
  
    /* reinitialise CVODE */
    MDUCresetODE(mem);
  }
  else { /* ----- call 2D flamelet solver ----- */
    /* set such that it does not integrate past the end of time step */
    flag = SDialSetStopTime(sdial->mem, tEnd);

    /* ----- call SUNDIALS to integrate to next time ----- */
    do {
      if (driver->SDIAL == 0)
        flag = CVode(sdial->mem, tEnd, sdial->Y, &time->t, sdial->itask);
      else if (driver->SDIAL == 1)
        flag = IDASolve(sdial->mem, tEnd, &time->t, sdial->Y, sdial->Yp,
                        sdial->itask);
      
      if (flag < 0)
      fprintf(stderr, "SUNDIALS execution failed: flag = %i\n", flag);

      /* get solver data about last time step */
      SDialGetNumSteps(sdial->mem, &nst);
      SDialGetLastStep(sdial->mem, &dtPre);
      sdial->totalIt += nst;
      
      /* update sources */
      if (io->outputHR || driver->output_tsrc) {
        for (k=flame->k0; k<=flame->kN; k++) {
          row = k*flame->n1*flame->nVars;
          r   = k*flame->n1;
          for (j=flame->j0; j<=flame->jN; j++) {
            col = j*flame->nVars;
            if (io->outputHR)
              for (i=0; i<chem->nSpecSolve; i++)
                flame->chemSrcY[row+col+i] += 
                  flame->chemSrcY_tmp[row+col+i]*dtPre;
            
            if (driver->output_tsrc) {
              flame->chem_src[r+j] += flame->chem_temp[r+j]*dtPre;
              flame->conv_src[r+j] += flame->conv_temp[r+j]*dtPre;
              flame->diff_src[r+j] += flame->diff_temp[r+j]*dtPre;
            }
          } /* end loop over columns */
        } /* end loop over rows */
      }
      
      /* TODO: print monitor */
      if (myid == 0) {
        fprintf(io->monitor, "%5li %15.6e %15.6e\n",nst,time->t,dtPre);
        fflush(io->monitor);
      }
    } while (flag != sdial->iret);
    
    /* finish averaging of heat release rate over time step */
    if (io->outputHR) {
      for (k=0; k<flame->n2; k++)
        for (j=0; j<flame->n1; j++)
          flame->HRrate[j] = 0.0;
      MDUCcalcHRrate(mem);
    }
  } /* end solve */
  
  if (io->outputMu || chem->soot > 0) {
    for (k=flame->k0; k<flame->kN; k++) {
      row = k*flame->n1*flame->nVars;
      r   = k*flame->n1;
      for (j=0; j<flame->n1; j++) {
        col = j*flame->nVars;
        flame->mu[r+j] = calcMixViscosity(chem->nSpecies,
                                          flame->Y[row+col+flame->sT],
                                          &flame->Y[row+col],chem);
      }
    }

    /* TODO: soot mode */
  }

#ifdef ACCESS_MPI
  /* make sure all processors have the up-to-date solution */
/*   MPI_Allgatherv(&flame->Y[mpi->kStart[myid]],mpi->kCount[myid], */
/*                  MPI_DOUBLE_PRECISION,flame->Y,mpi->kCount,mpi->kStart, */
/*                  MPI_DOUBLE_PRECISION,mpi->comm); */
/*   if (driver->SDIAL == 1) */
/*     MPI_Allgatherv(&flame->Yp[mpi->kStart[myid]],mpi->kCount[myid], */
/*                    MPI_DOUBLE_PRECISION,flame->Yp,mpi->kCount,mpi->kStart, */
/*                    MPI_DOUBLE_PRECISION,mpi->comm); */

  MPI_Gatherv(&flame->Y[mpi->kStart[myid]],mpi->kCount[myid],
              MPI_DOUBLE_PRECISION,flame->Y,mpi->kCount,mpi->kStart,
              MPI_DOUBLE_PRECISION,0,mpi->comm);

  if (driver->SDIAL == 1)
    MPI_Gatherv(&flame->Yp[mpi->kStart[myid]],mpi->kCount[myid],
                MPI_DOUBLE_PRECISION,flame->Yp,mpi->kCount,mpi->kStart,
                MPI_DOUBLE_PRECISION,0,mpi->comm);

  /* retrieve variables to be output onto the root node */  
  MPI_Gatherv(&flame->cpmix[mpi->rStart[myid]],mpi->rCount[myid],
              MPI_DOUBLE_PRECISION,flame->cpmix,mpi->rCount,mpi->rStart,
              MPI_DOUBLE_PRECISION,0,mpi->comm);
  MPI_Gatherv(&flame->gamma[mpi->rStart[myid]],mpi->rCount[myid],
              MPI_DOUBLE_PRECISION,flame->gamma,mpi->rCount,mpi->rStart,
              MPI_DOUBLE_PRECISION,0,mpi->comm);
  MPI_Gatherv(&flame->rho[mpi->rStart[myid]],mpi->rCount[myid],
              MPI_DOUBLE_PRECISION,flame->rho,mpi->rCount,mpi->rStart,
              MPI_DOUBLE_PRECISION,0,mpi->comm);
  MPI_Gatherv(&flame->R[mpi->rStart[myid]],mpi->rCount[myid],
              MPI_DOUBLE_PRECISION,flame->R,mpi->rCount,mpi->rStart,
              MPI_DOUBLE_PRECISION,0,mpi->comm);
  MPI_Gatherv(&flame->mu[mpi->rStart[myid]],mpi->rCount[myid],
              MPI_DOUBLE_PRECISION,flame->mu,mpi->rCount,mpi->rStart,
              MPI_DOUBLE_PRECISION,0,mpi->comm);
  /* retrieve source terms for heat release computation */
  if (io->outputHR)
    MPI_Gatherv(&flame->chemSrcY[mpi->kStart[myid]],mpi->kCount[myid],
                MPI_DOUBLE_PRECISION,flame->chemSrcY,mpi->kCount,mpi->kStart,
                MPI_DOUBLE_PRECISION,0,mpi->comm);

  if (chem->soot > 0) {
    MPI_Gatherv(&flame->mu[mpi->rStart[myid]],mpi->rCount[myid],
                MPI_DOUBLE_PRECISION,flame->chem_src,mpi->rCount,mpi->rStart,
                MPI_DOUBLE_PRECISION,0,mpi->comm);
  }
#else
  /* copy the final solution into the flame vector */
  copyArray(flame->n1*flame->n2*flame->nVars,sdial->Ydata,flame->Y);
#endif

  /* update normalisation for next time step */
  flow->P_m1       = flow->P;
  flow->P_pre      = flow->P;
  flow->phi1min_m1 = flow->phi1min;
  flow->phi1max_m1 = flow->phi1max;
  flow->phi2min_m1 = flow->phi2min;
  flow->phi2max_m1 = flow->phi2max;

  /* copy scalar dissipation rates */
  copyArray(flame->n1*flame->n2, flow->chi1, flow->chi1_m1);
  copyArray(flame->n1*flame->n2, flow->chi2, flow->chi2_m1);
  copyArray(flame->n1*flame->n2, flow->chi12, flow->chi12_m1);

  /* update flow time */
  time->t_m1   = time->t;
  time->dt_pre = dt;
  time->dt_cvode_pre = dtPre;
  time->nt++;

  if (driver->ENERGY) {
    /* copy most current enthalpy from solution array for next time step */
    extractVarArray(flame->nVars,flame->n1,flame->n2,flame->sT,
                    flame->Y,flame->H);
    
    /* make sure the temperature is up-to-date */
    for (k=0; k<flame->n2; k++) {
      row = k*flame->n1*flame->nVars;
      r   = k*flame->n1;
      for (j=0; j<flame->n1; j++) {
        col = j*flame->nVars;
        flame->T[r+j] = calcTfromH(chem->nSpecies,chem->nSpecSolve,
                                   flame->Y[row+col+flame->sT],
                                   &flame->Y[row+col],1000,chem);
      }
    }
  }
  else { /* compute the total enthalpy if output to file */
    if (io->outputHtot)
      calcTotalEnthalpy(flame->n1,flame->n2,flame->nVars,flame->sT,flame->Y,
                        flame->H,chem);
  }
  
#ifdef ACCESS_MPI
  /* distribute solution data to all processors for output */
/*   MPI_Allgatherv(&flame->Y[mpi->kStart[myid]],mpi->kCount[myid], */
/*                  MPI_DOUBLE_PRECISION,flame->Y,mpi->kCount,mpi->kStart, */
/*                  MPI_DOUBLE_PRECISION,mpi->comm); */
#endif /* ACCESS_MPI */

  /* write solution file at this step (change file name here if desired) */
  if (myid == 0 && flame->nDim > 0) {
    if (io->freqType == 0 && time->nt%io->writeFreq == 0) {
      createDataFileName(io->outputData, io, time, outputFileName);
      writeFlamelet(outputFileName, io, driver, time, flame, flow);
    }
    else if (io->freqType == 1 && (fmod(time->t,io->tFreq) < 0.5*time->dt ||
                                   time->t > io->tWrite+io->tFreq) ) {
      createDataFileName(io->outputData,io,time,outputFileName);
      writeFlamelet(outputFileName,io,driver,time,flame,flow);
      io->tWrite = time->t;
    }
  }
}
/* -------------------------------------------------------------------------- */

int cvodeFlamelet2D(realtype t, N_Vector u, N_Vector udot, void *data)
{
  realtype *Y, *RHS;
  /* assign pointers to vector memory for solution vector and right-hand-side */
#ifdef ACCESS_MPI
  Y   = NV_DATA_P(u);
  RHS = NV_DATA_P(udot);

  /* communicate data from neighbouring processors */
  bCommFlamelet2D(Y, data);

  sdialStruct *sdial  = (sdialStruct*) data;
  flameStruct *flame  = sdial->args[3];

  /* use the data now copied/communicated in flame->Y for calculating RHS */
  Y = flame->Y;
#else /* SERIAL */
  Y   = NV_DATA_S(u);
  RHS = NV_DATA_S(udot);
#endif /* ACCESS_MPI */

  /* compute the current flamelet RHS */
  flameletRHS_2D(t, Y, RHS, data);
  
  return 0;
}
/* -------------------------------------------------------------------------- */

#ifdef ACCESS_MPI
#ifdef SUNDIALS_FROM_2_7
int cvodeFlamelet2D_local(long int Nlocal, realtype t, N_Vector u, N_Vector udot,
                          void *data)
#else
#ifdef SUNDIALS_FROM_3_0
int cvodeFlamelet2D_local(sunindextype Nlocal, realtype t, N_Vector u, N_Vector udot,
                          void *data)
#else
int cvodeFlamelet2D_local(int Nlocal, realtype t, N_Vector u, N_Vector udot,
                          void *data)
#endif /* SUNDIALS_FROM_3_0 */
#endif /* SUNDIALS_FROM_2_7 */
{
  realtype *Y, *RHS;
  /* assign pointers to beginning of Sundials vector memory */
  Y   = NV_DATA_P(u);
  RHS = NV_DATA_P(udot);

  sdialStruct *sdial  = (sdialStruct*) data;
  flameStruct *flame  = sdial->args[3];
  mpiStruct   *mpi    = sdial->args[7];

  /* copy current solution to flamelet memory for RHS computation */
  copyArray(mpi->kCount[mpi->myid], Y, &flame->Y[mpi->kStart[mpi->myid]]);
  Y = flame->Y;

  /* compute the current flamelet RHS */
  flameletRHS_2D(t, Y, RHS, data);
  return 0;
}
/* -------------------------------------------------------------------------- */

int bCommFlamelet2D(realtype *Y, void *data)
{
  sdialStruct *sdial  = (sdialStruct*) data;
  flameStruct *flame  = sdial->args[3];
  mpiStruct   *mpi    = sdial->args[7];

  int i, myid=mpi->myid;
  MPI_Status  status;
  MPI_Request request[4];

  /* copy current solution to flamelet memory for RHS computation */
  copyArray(mpi->kCount[myid], Y, &flame->Y[mpi->kStart[myid]]);

  /* make non-blocking request for required cells from all boundary processes */
  for (i=0; i<mpi->nBcomm; i++)
    MPI_Irecv(&flame->Y[mpi->bRcv[i]], mpi->bCount[i], MPI_DOUBLE_PRECISION,
              mpi->bComm[i], mpi->bRcv[i], mpi->comm, &request[i]);

  /* send all required cells to boundary processors */
  for (i=0; i<mpi->nBcomm; i++)
    MPI_Send(&flame->Y[mpi->bSend[i]], mpi->bCount[i], MPI_DOUBLE_PRECISION,
             mpi->bComm[i], mpi->bSend[i], mpi->comm);

  /* make sure that all data has been successfully received before continuing */
  for (i=0; i<mpi->nBcomm; i++)
    MPI_Wait(&request[i], &status);

  return 0;
}
/* -------------------------------------------------------------------------- */
#endif /* ACCESS_MPI */

int idaFlamelet2D(realtype t, N_Vector y, N_Vector yp, N_Vector r, void *data)
{
  return 0;
}
/* -------------------------------------------------------------------------- */

void flameletRHS_2D(double t, double *Y, double *RHS, void *args)
{
  int i,j,k,row,col,rowL,r,j0D;
  double *T,*RHS_s;

  /* dereference the necessary structs */
  sdialStruct  *sdial  = (sdialStruct*)  args;
  driveStruct  *driver = (driveStruct*)  sdial->args[0];
  chemStruct   *chem   = (chemStruct*)   sdial->args[1];
  flowStruct   *flow   = (flowStruct*)   sdial->args[2];
  flameStruct  *flame  = (flameStruct*)  sdial->args[3];
  timeStruct   *time   = (timeStruct*)   sdial->args[4];
  engineStruct *engine = (engineStruct*) sdial->args[5];

  /* update the flow variables for current time */
  updateFlow2D(t, sdial->row, driver, time, flame, flow);

  /* set the current temperature from the solution array */
  T = (double*) calloc(flame->n1*flame->n2, sizeof(double));
  extractVarArray(flame->nVars, flame->n1, flame->n2, flame->sT, Y, T);

  /* compute the r.h.s. for the species */
  flameletY_RHS_2D(driver, chem, flame, flow, T, Y, RHS);

  /* compute the r.h.s for the temperature */
  flameletT_RHS_2D(driver, chem, flame, flow, Y, RHS);

  /* copy singularity if necessary */
  if (driver->TRANS_Z2 > 0) { /* !!!!!!! DEBUG */
    /* pick the appropriate column for the current transformation */
    if (driver->TRANS_Z2 == 1)
      j0D = flame->n1-1;
    else if (driver->TRANS_Z2 == 2)
      j0D = 0;
    col = j0D*flame->nVars;

    /* assign memory for singularity RHS and fill */
    RHS_s = (double*) calloc(flame->nVars, sizeof(double));
    for (i=0; i<flame->nVars; i++)
      RHS_s[i] = RHS[col+i];

#ifdef ACCESS_MPI
    mpiStruct   *mpi    = sdial->args[7];
    /* broadcast the first row solution from root */
    MPI_Bcast(RHS_s,flame->nVars,MPI_DOUBLE_PRECISION,0,mpi->comm);
    MPI_Bcast(&flame->cpmix[j0D],1,MPI_DOUBLE_PRECISION,0,mpi->comm);
    MPI_Bcast(&flame->gamma[j0D],1,MPI_DOUBLE_PRECISION,0,mpi->comm);
    MPI_Bcast(&flame->rho[j0D],1,MPI_DOUBLE_PRECISION,0,mpi->comm);
    MPI_Bcast(&flame->Mw[j0D],1,MPI_DOUBLE_PRECISION,0,mpi->comm);
    MPI_Bcast(&flame->R[j0D],1,MPI_DOUBLE_PRECISION,0,mpi->comm);
    MPI_Bcast(&flame->mu[j0D],1,MPI_DOUBLE_PRECISION,0,mpi->comm);
    if (driver->output_tsrc) {
      MPI_Bcast(&flame->chem_src[j0D],1,MPI_DOUBLE_PRECISION,0,mpi->comm);
      MPI_Bcast(&flame->diff_src[j0D],1,MPI_DOUBLE_PRECISION,0,mpi->comm);
      MPI_Bcast(&flame->conv_src[j0D],1,MPI_DOUBLE_PRECISION,0,mpi->comm);
    }
#endif /* ACCESS_MPI */
      
    /* loop over current rows */
    for (k=flame->k0; k<=flame->kN; k++) {
      r    = k*flame->n1;
      row  = k*flame->n1*flame->nVars;
      rowL = (k-flame->k0)*flame->n1*flame->nVars;
      /* fprintf(stderr, "%i %i %i %i\n", k, r, row, rowL); */
      
      for (i=0; i<flame->nVars; i++) {
/*         if (fabs(RHS[rowL+col+i] - RHS_s[i]) > 1.0e-6) */
/*           fprintf(stderr, "%i %i %i: %f --> %f\n",k,i,rowL+col+i, */
/*                   RHS[rowL+col+i],RHS_s[i]); */
        /* TODO: DEBUG, this for some reason doesn't work */
        /* RHS[rowL+col+i]            = RHS_s[i]; */
        flame->chemSrcY[row+col+i] = flame->chemSrcY[col+i];
      }
      flame->cpmix[r+j0D] = flame->cpmix[j0D];
      flame->rho[r+j0D]   = flame->rho[j0D];
      flame->gamma[r+j0D] = flame->gamma[j0D];
      flame->Mw[r+j0D]    = flame->Mw[j0D];
      flame->R[r+j0D]     = flame->R[j0D];
      flame->mu[r+j0D]    = flame->mu[j0D];
      
      if (driver->output_tsrc) {
        flame->chem_src[r+j0D] = flame->chem_src[j0D];
        flame->diff_src[r+j0D] = flame->conv_src[j0D];
        flame->conv_src[r+j0D] = flame->diff_src[j0D];
      }
    } /* end loop over rows */
  } /* end singularity */

  free(T);
  if (driver->TRANS_Z2 > 0) free(RHS_s);
}
/* -------------------------------------------------------------------------- */

void flameletY_RHS_2D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                      flowStruct *flow, double *T, double *Y, double *RHS)
{
  int    i,j,k,jst,kst,r,row,col,rowL,iX1;
  double Dx_Y,Dxx_Y,Dxy_Y,Dyy_Y; /* spatial derivatives */
  double vT1,vT2,vH1,vH2;        /* additional velocities in temperature eq. */
  double src,diff,conv;          /* source, diffusion, and convection terms */

  /* set some local pointers for readability */
  double P     = flow->P;
  double dPdt  = flow->Dt_P;
  int nSpecies = chem->nSpecSolve;
  double *cp   = chem->cp;
  
  /* loop over mesh points to solve the species equations */

  for (k=flame->k0; k<=flame->kN; k++) {
    row  = k*(flame->n1*flame->nVars); /* start of current row */
    rowL = (k-flame->k0)*(flame->n1*flame->nVars);
    r    = k*flame->n1;                /* start of current row (single var.) */
    kst  = stencil_y(k, flame);        /* start of stencil */

    /* precompute the derivative in the y direction to enable cross-derivative
       computation in column loop */
    for (j=flame->j0; j<=flame->jN; j++) {
      col = j*flame->nVars;
      for (i=0; i<nSpecies; i++)
        flame->dY.y[col+i] = grad_y(kst,k,j,i,flame->nVars,Y,flame);
    }

    /* loop over columns */
    for (j=flame->j0; j<=flame->jN; j++) {
      col  = j*flame->nVars;      /* column index */
      jst  = stencil_x(j, flame); /* start of stencil */
      
      /* compute local density and pressure source term */
      flame->Mw[r+j]   = calcMixtureW(nSpecies,&Y[row+col],chem->W);
      flame->rho[r+j]  = flame->Mw[r+j]*P/(flame->RGAS*T[r+j]);
      flame->srcP[r+j] = 1.0/flame->rho[r+j]*dPdt;

      /* add additional terms implicitly to convective velocity if current
         direction is enthalpy */
      vH1 = vH2 = 0.0;
      if (driver->H_DIM != 0) {
        if (driver->DIM == driver->H_DIM && flow->D_phi1 != 0.0)
          vH1 = (flame->srcP[r+j] - flow->avgSdot[j])/flow->D_phi1;
        else if (driver->DIM != driver->H_DIM && flow->D_phi2 != 0.0)
          vH2 = (flame->srcP[r+j] - flow->avgSdot[k])/flow->D_phi2;
      }
      
      /* get enthalpy and specific heat from chemical mechanism */
      chemComputeThermoData(chem->h,chem->cp,T[r+j],chem);

      /* get specific heat of the mixture */
      flame->cpmix[r+j] = calcMixtureCp(nSpecies,chem->cp,&Y[row+col]);
      flame->R[r+j]     = flame->RGAS/flame->Mw[r+j];
      flame->gamma[r+j] = flame->cpmix[r+j]/(flame->cpmix[r+j]-flame->R[r+j]);
      
      /* compute the mixture viscosity */
      if (driver->VISC_MODE)
        flame->mu[r+j] = calcMixViscosity(nSpecies,T[r+j],&Y[row+col],chem);

      /* compute and store concentrations of all species at local point */
      for (i=0; i<nSpecies; i++)
        chem->c[i] = flame->rho[r+j]*Y[row+col+i]/chem->W[i];
      
      /* get rates from mechanism */
      chemComputeProductionRates(chem->cdot,chem->w,chem->K,chem->c,chem->M,
                                 T[r+j],P,chem);
      
      flame->srcH[r+j] = flame->sumCp[r+j] = flame->sumCp2[r+j] = 0.0;
      /* loop over species */
      for (i=0; i<nSpecies; i++) {
        /* calculate the gradients in normalised space */
        Dx_Y  = grad_x(jst,  k, j, i, flame->nVars, Y, flame);
        Dxy_Y = grad_x(jst,  0, j, i, flame->nVars, flame->dY.y, flame);
        Dxx_Y = grad_xx(jst, k, j, i, flame->nVars, Y, flame);
        Dyy_Y = grad_yy(kst, k, j, i, flame->nVars, Y, flame);
        
        /* ----- compute diffusion term -----*/
        diff =   0.5*flame->chiA[r+j]*Dxx_Y + flame->chiB[r+j]*Dxy_Y
               + 0.5*flame->chiC[r+j]*Dyy_Y;
        
        /* ----- compute convection term ----- */
        conv =   (flame->v1[r+j] + vH1)*Dx_Y
               + (flame->v2[r+j] + vH2)*flame->dY.y[col+i];
          
        /* ----- compute chemical source term ----- */
        src  = chem->cdot[i]*chem->W[i]/flame->rho[r+j];
          
        /* ----- compute summation terms required for energy equation ----- */
        if (driver->ENERGY)
          flame->srcH[r+j] += chem->h[i]*Y[row+col+i];
        else { /* temperature */
          flame->sumCp[r+j]  += cp[i]*(flame->chiA[r+j]*Dx_Y + 
                                       flame->chiB[r+j]*flame->dY.y[col+i]);
          flame->sumCp2[r+j] += cp[i]*(flame->chiB[r+j]*Dx_Y + 
                                       flame->chiC[r+j]*flame->dY.y[col+i]);
          flame->srcH[r+j]   += chem->h[i]*src;  /* for chemical source term */
        }
          
        /* ----- update RHS of the species equations ----- */
        RHS[rowL+col+i] = diff - driver->CONV_MODE*conv + src;
        
        /* store chemical source term for heat release computation */
        flame->chemSrcY_tmp[(r+j)*flame->nVars+i] = src;
      } /* end loop over species */
    } /* end loop over columns */
  } /* end loop over rows */

#ifdef ACCESS_MPI
  /* if parallel, make sure that that the mixture cp is updated in the
     ghost cells */
  for (k=flame->k0m; k<=flame->k0; k++) {
    row = k*flame->n1*flame->nVars;
    r   = k*flame->n1;
    for (j=flame->j0; j<=flame->jN; j++) {
      col = j*flame->nVars;
      chemComputeThermoData(chem->h,chem->cp,T[r+j],chem);
      flame->cpmix[r+j] = calcMixtureCp(nSpecies,chem->cp,&Y[row+col]);
    }
  }

  for (k=flame->kN+1; k<=flame->kNp; k++) {
    row = k*flame->n1*flame->nVars;
    r   = k*flame->n1;
    for (j=flame->j0; j<=flame->jN; j++) {
      col = j*flame->nVars;
      chemComputeThermoData(chem->h,chem->cp,T[r+j],chem);
      flame->cpmix[r+j] = calcMixtureCp(nSpecies,chem->cp,&Y[row+col]);
    }
  }
#endif /* ACCESS_MPI */
}
/* -------------------------------------------------------------------------- */

void flameletT_RHS_2D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                      flowStruct *flow, double *Y, double *RHS)
{
  int    i,j,k,jst,kst,r,row,col,rowL;
  double Dx_Y,Dxx_Y,Dyy_Y,Dxy_Y; /* spatial derivatives */
  double Dx_cpmix,Dy_cpmix;       
  double vT1,vT2,vH1,vH2;        /* additional velocities in temperature eq. */
  double src,diff,conv;          /* source, diffusion, and convection terms */
  
  /* ----- solve the temperature equation ----- */
  for (k=flame->k0; k<=flame->kN; k++) {
    row  = k*(flame->n1*flame->nVars); /* start of current row */
    rowL = (k-flame->k0)*(flame->n1*flame->nVars);
    r    = k*flame->n1;
    kst  = stencil_y(k, flame);        /* start of stencil */

    /* precompute the derivative in the y direction to enable cross-derivative
       computation in column loop */
    for (j=flame->j0; j<=flame->jN; j++) {
      col = j*flame->nVars + flame->sT;
      flame->dY.y[col] = grad_y(kst,k,j,flame->sT,flame->nVars,Y,flame);
    }

    for (j=flame->j0; j<=flame->jN; j++) {
      col  = j*flame->nVars + flame->sT;  /* column index */
      jst  = stencil_x(j, flame);         /* start of stencil */
      
      /* compute gradients in x1 direction */
      Dx_Y     = grad_x(jst, k,j,flame->sT,flame->nVars,Y,flame);
      Dxx_Y    = grad_xx(jst,k,j,flame->sT,flame->nVars,Y,flame);
      Dxy_Y    = grad_x(jst, 0,j,flame->sT,flame->nVars,flame->dY.y,flame);
      Dyy_Y    = grad_yy(kst,k,j,flame->sT,flame->nVars,Y,flame);
      Dx_cpmix = grad_x(jst, k,j,0,1,flame->cpmix,flame);
      Dy_cpmix = grad_y(kst, k,j,0,1,flame->cpmix,flame);
      
      /* ----- compute diffusion term ----- */
      diff =   0.5*flame->chiA[r+j]*Dxx_Y + flame->chiB[r+j]*Dxy_Y
             + 0.5*flame->chiC[r+j]*Dyy_Y;

      /* ----- compute convection term ----- */
      /* add additional terms implicitly to convective velocity if current
         direction is enthalpy */
      vH1 = vH2 = 0.0;
      if (driver->H_DIM != 0) {
        if (driver->DIM == driver->H_DIM && flow->D_phi1 != 0.0)
          vH1 = (flame->srcP[r+j] - flow->avgSdot[j])/flow->D_phi1;
        else if (driver->DIM != driver->H_DIM && flow->D_phi2 != 0.0)
          vH2 = (flame->srcP[r+j] - flow->avgSdot[k])/flow->D_phi2;
      }

      /* compute additional velocities for temperature equation (only if 
         not comparing to DNS solution) */
      vT1 = vT2 = 0.0;
      if (!driver->DNS) {
        vT1 = - 0.5/flame->cpmix[r+j]*(flame->chiA[r+j]*Dx_cpmix + 
                                       flame->chiB[r+j]*Dy_cpmix +
                                       flame->sumCp[r+j]);
        vT2 = - 0.5/flame->cpmix[r+j]*(flame->chiB[r+j]*Dx_cpmix + 
                                       flame->chiC[r+j]*Dy_cpmix +
                                       flame->sumCp2[r+j]);
      }

      /* compute total contribution */
      conv = (flame->v1[r+j] + vT1 + vH1)*Dx_Y +
             (flame->v2[r+j] + vT2 + vH2)*flame->dY.y[col];
      
      /* ----- compute chemical source term ----- */
      src  = 1.0/flame->cpmix[r+j]*(flame->srcP[r+j] - flame->srcH[r+j]);

      /* ----- copy sources over if needed ----- */
      if (driver->output_tsrc) {
        flame->diff_temp[r+j] = diff;
        flame->conv_temp[r+j] = conv;
        flame->chem_temp[r+j] = src;
      }
      
      /* ----- update RHS of temperature equation ----- */
      RHS[rowL+col] = diff - driver->CONV_MODE*conv + src;
    } /* end loop over columns */
  } /* end loop over rows */
}
/* -------------------------------------------------------------------------- */
