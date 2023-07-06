#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "flameletInit.h"
#include "flamelet1D.h"
#include "flamelet2D.h"
#include "solve1D.h"
#include "chemMech.h"
#include "transport.h"
#include "accessChem.h"
#include "metric.h"
#include "solveUtils.h"
#include "utils.h"
#include "fileIO.h"
#include "display.h"
#include "sootHMOM.h"
#include "engine.h"
#include "engineUtils.h"
#include "utils.h"
#include "accessMDUC.h"
#include "rk4solver.h"

void solveFlamelet1D(double dt,double P,double phimin,double phimax,double *chi,
                     int inert,void *mem)
{
  int myid = 0;
  int ierr;
  

  /* dereference internal flamelet solver structs */
  mducStruct   *mduc   = (mducStruct*) mem;
  ioStruct     *io     = mduc->io;
  driveStruct  *driver = mduc->driver;
  chemStruct   *chem   = mduc->chem;
  timeStruct   *time   = mduc->time;
  flowStruct   *flow   = mduc->flow;
  flameStruct  *flame  = mduc->flame;
  sdialStruct  *sdial  = mduc->sdial;
  rkStruct     *rk     = mduc->rk;
  engineStruct *engine = mduc->engine;
#ifdef ACCESS_MPI
  mpiStruct *mpi = mduc->mpi; /* dereference mpi info */
  myid = mpi->myid;
  realtype *Y;
#endif /* ACCESS_MPI */

  int  i, j, k, ind, jMaxT;
  char outputFileName[256];
  char rateFile[256];
  char errorFileName[256];
  
  /* set CVODE variables */
  int      flag, flag_global;        /* flags for CVODE return */
  realtype tEnd_,tEnd, dtPre = dt; /* final time to integrate to and previous time step */
  long int nst;         /* current time step number */

  /* check whether Sundials needs to be initialized */
  if (sdial->mem == NULL) {
    if (driver->SDIAL == 0)
      cvodeInitFlamelet(time->t, flame->Y, sdial);
    else if (driver->SDIAL == 1)
      idaInitFlamelet(time->t, flame->Y, sdial);
    else if (driver->SDIAL == -1)
      rk4initFlamelet(time->t, flame->Y, mem);
  }

  flowBound *phiRaw, *phiHat, *Dt_phiRaw, *Dt_phiHat, *DtHat_phi;
  phiRaw    = flow->phiRaw[1];
  phiHat    = flow->phiHat[1];
  Dt_phiRaw = flow->Dt_phiRaw[1];
  Dt_phiHat = flow->Dt_phiHat[1];
  DtHat_phi = flow->DtHat_phi[1];

  checkMonitor(io);

  /* update time step to take */
  time->dt = dt;
  tEnd = time->t_n = time->t + dt;
  
  /* check to make sure the flow field has been initialised for the flamelet */
  if (flow->P_m1 == 0.0)
    mducErr("initial pressure has not been set","solveFlamelet1D");

  /* TODO: move the filter somewhere else to clean up */
  /* store the previous bounds */
  for (j=flow->fhat.nW-1; j>0; j--) {
    phiRaw[j].min = phiRaw[j-1].min;
    phiRaw[j].max = phiRaw[j-1].max;
    phiRaw[j].del = phiRaw[j-1].del;
    phiHat[j].min = phiHat[j-1].min;
    phiHat[j].max = phiHat[j-1].max;
    phiHat[j].del = phiHat[j-1].del;

    /* store previous derivatives */
    Dt_phiRaw[j].min = Dt_phiRaw[j-1].min;
    Dt_phiRaw[j].max = Dt_phiRaw[j-1].max;
    Dt_phiRaw[j].del = Dt_phiRaw[j-1].del;

    Dt_phiHat[j].min = Dt_phiHat[j-1].min;
    Dt_phiHat[j].max = Dt_phiHat[j-1].max;
    Dt_phiHat[j].del = Dt_phiHat[j-1].del;

    DtHat_phi[j].min = DtHat_phi[j-1].min;
    DtHat_phi[j].max = DtHat_phi[j-1].max;
    DtHat_phi[j].del = DtHat_phi[j-1].del;
  }
  /* update the current bounds */
  phiRaw[0].min = phimin;
  phiRaw[0].max = phimax;
  phiRaw[0].del = phimax - phimin;

  /* filter the bounds */
  phiHat[0].min = phiHat[0].max = 0.0;
  for (j=0; j<flow->fhat.nW; j++) {
    phiHat[0].min += phiRaw[j].min*flow->fhat.w[j];
    phiHat[0].max += phiRaw[j].max*flow->fhat.w[j];
  }

  /* ----- compute the derivatives ----- */
  /* dPhi/dt */
  Dt_phiRaw[0].min = (phiRaw[0].min - phiRaw[1].min)/dt;
  Dt_phiRaw[0].max = (phiRaw[0].max - phiRaw[1].max)/dt;
  /* d<Phi>/dt*/
  Dt_phiHat[0].min = (phiHat[0].min - phiHat[1].min)/dt;
  Dt_phiHat[0].max = (phiHat[0].max - phiHat[1].max)/dt;
  /* compute the filtered derivative <dPhi/dt> */
  DtHat_phi[0].min = DtHat_phi[0].max = 0.0;
  for (j=0; j<flow->fhat.nW; j++) {
    DtHat_phi[0].min += Dt_phiRaw[j].min*flow->fhat.w[j];
    DtHat_phi[0].max += Dt_phiRaw[j].max*flow->fhat.w[j];
  }
  
  /* assign flow variables to struct */
  flow->P       = P;
  if (driver->PHI_HAT == 0) { /* no filtering */
    flow->phi1min = phimin;
    flow->phi1max = phimax;
    flow->phi     = flow->phiRaw;
    flow->Dt_phi  = flow->Dt_phiRaw;
  }
  else if (driver->PHI_HAT == 1) { /* filtering */
    flow->phi1min = phiHat[0].min;
    flow->phi1max = phiHat[0].max;
    flow->phi     = flow->phiHat;
    flow->Dt_phi  = flow->Dt_phiHat;
  }
  for (j=0; j<flame->n1; j++) flow->chi1[j] = chi[j];

  /* compute pressure rate */
  flow->Dt_P = (flow->P - flow->P_m1)/time->dt;

  /* get normalisation for variables */
  flow->D_phi1     = flow->phi1max - flow->phi1min;
  flow->D_phi1_m1  = flow->phi1max_m1 - flow->phi1min_m1;
  /* Get difference from corrected value */
  if (driver->H_DIM != 0 && driver->DNS == 1) {
    fprintf(stderr, "  using DNS\n");
    flow->Dt_phi1min = (flow->phi1min_corr - flow->phi1min_m1)/time->dt;
    flow->Dt_phi1max = (flow->phi1max_corr - flow->phi1max_m1)/time->dt;
  }
  else {
    flow->Dt_phi1min = (flow->phi1min - flow->phi1min_m1)/time->dt;
    flow->Dt_phi1max = (flow->phi1max - flow->phi1max_m1)/time->dt;
  }

  if (io->debug) {
    fprintf(stderr, "");
    fprintf(stderr, "%10s | %31s | %31s\n", "", "Phi", "<Phi>");
    fprintf(stderr, "%10s | %15s | %15s | %15s | %15s\n", "", "min", "max",
            "min", "max");
    for (j=0; j<flow->fhat.nW; j++) {
      fprintf(stderr, "%5s - %2i | %15.6e | %15.6e | %15.6e | %15.6e\n", "n", j,
              phiRaw[j].min, phiRaw[j].max, phiHat[j].min, phiHat[j].max);
    }
    
    fprintf(stderr, "%10s | %31s | %31s | %31s\n", "", "dPhi/dt", "d<Phi>/dt",
            "<dPhi/dt>");
    fprintf(stderr, "%10s | %15s | %15s | %15s | %15s | %15s | %15s\n", "",
            "min", "max", "min", "max", "min", "max");
    for (j=0; j<flow->fhat.nW; j++) {
      fprintf(stderr, "%5s - %2i | %15.6e | %15.6e | %15.6e | %15.6e | "
              "%15.6e | %15.6e\n","n",j,Dt_phiRaw[j].min,Dt_phiRaw[j].max,
              Dt_phiHat[j].min,Dt_phiHat[j].max, DtHat_phi[j].min,
              DtHat_phi[j].max);
    }
    
    fprintf(stderr, "  phi min:      %15.6e\n", flow->phi1min);
    fprintf(stderr, "  phi min(n-1): %15.6e\n", flow->phi1min_m1);
    fprintf(stderr, "  phi max:      %15.6e\n", flow->phi1max);
    fprintf(stderr, "  phi max(n-1): %15.6e\n", flow->phi1max_m1);
    fprintf(stderr, "  dPhiMin/dt:   %15.6e\n", flow->Dt_phi1min);
    fprintf(stderr, "  dPhiMax/dt:   %15.6e\n", flow->Dt_phi1max);
    
    fprintf(stderr, "Boundary states:\n");
    fprintf(stderr, "  rho0: %15.6e\n", flame->rho[0]);
    fprintf(stderr, "  T0:   %15.6e\n", flame->Y[flame->sT]);
    fprintf(stderr, "  rho1: %15.6e\n", flame->rho[flame->n1-1]);
    fprintf(stderr, "  T1:   %15.6e\n", flame->Y[(flame->n1-1)*flame->nVars+flame->sT]);
  }

  /* check to make sure that bounds are not changing faster than dP/dt */
  if (driver->H_DIM != 0) {
    if (driver->ENGINE) {
      flow->Dt_phi1min = 1.0/flame->rho[0]*flow->Dt_P;
      flow->Dt_phi1max = 1.0/flame->rho[flame->n1-1]*flow->Dt_P;
    }
    else if (flow->Dt_phi1min < 1.0/flame->rho[0]*flow->Dt_P) {
      fprintf(stderr, "  Clipping dPhi_min/dt: %15.6e %15.6e\n",
              flow->Dt_phi1min, 1.0/flame->rho[0]*flow->Dt_P);
      fprintf(stderr, "flow->Dt_phi1min: %15.6e\n", flow->Dt_phi1min);
      fprintf(stderr, "flame->rho[0]:    %15.6e\n", flame->rho[0]);
      fprintf(stderr, "flow->Dt_P:       %15.6e\n", flow->Dt_P);
      flow->Dt_phi1min = 1.0/flame->rho[0]*flow->Dt_P;
    }
    if (flow->Dt_phi1max > 1.0/flame->rho[flame->n1-1]*flow->Dt_P) {
      fprintf(stderr, "  Clipping dPhi_max/dt: %15.6e %15.6e\n",
              flow->Dt_phi1max, 1.0/flame->rho[flame->n1-1]*flow->Dt_P);
      flow->Dt_phi1max = 1.0/flame->rho[flame->n1-1]*flow->Dt_P;
    }
  }

  /* write the monitor header if first time step*/
#ifdef ACCESS_MPI
  if (time->nt == time->nt0 && myid == mpi->host) {
#else
  if (time->nt == time->nt0) {
#endif /* ACCESS_MPI */
    fprintf(io->monitor, "# nst<1> time[s]<2> dt[s]<3> Tmax[K]<4> dt_max<5>\n");
    fflush(io->monitor);
  }

  /* the pre-step does the following:
     - store solution of Y(N) to Y(N-1)
     - computes transformation variables, v1 and chiA
     (these are invariant over the time step)
     - compute initial RHS */
  flameletPrestep1D(driver, flame, chem, flow, io, time);
  if (inert) { /* ----- update based on new boundary temperature ----- */
    int    k_pre;
    int    iX1    = (flame->n1-1)*flame->nVars; /* start of fuel boundary */
    double *hox   = (double*) calloc(chem->nSpecies, sizeof(double));
    double *hfuel = (double*) calloc(chem->nSpecies, sizeof(double));

    /* calculate the thermo data at the boundaries.  Assumes that boundary
       temperatures have been updated with MDUCsetTempBounds previously */
    chemComputeThermoData(hox,chem->cp,flame->Y[flame->sT], chem);
    chemComputeThermoData(hfuel,chem->cp,flame->Y[iX1+flame->sT],chem);

    if (!flame->H)
      flame->H = (double*) calloc(flame->n1, sizeof(double));

    flame->H[0] = flame->H[flame->n1-1] = 0.0; 

    for (i=0; i<chem->nSpecies; i++) {
      flame->H[0]           += hox[i]*flame->Y[i];
      flame->H[flame->n1-1] += hfuel[i]*flame->Y[iX1+i];
    }
    free(hox); free(hfuel);
    
    /* update the interior profile */
    for (j=1; j < flame->n1-1; j++) {
      k = j*flame->nVars;
      k_pre = (j-1)*flame->nVars;
      
      /* linearly interpolate the enthalpy */
      flame->H[j] = linInterp(flame->x1[j],flame->x1[0],flame->H[0],
                              flame->x1[flame->n1-1],flame->H[flame->n1-1]);
      
      /* iterate for temperature based on the enthalpy and species */
      flame->Y[k+flame->sT] = calcTfromH(chem->nSpecies,chem->nSpecSolve,
                                         flame->H[j],&flame->Y[k],
                                         flame->Y[k_pre+flame->sT],chem);
    } /* end loop over columns */
  
    /* reinitialise CVODE */
    MDUCresetODE(mem);
  }
  else { /* ----- call 1D flamelet solver ----- */
    /* set such that it does not integrate past the end of time step */
    if (driver->SDIAL >= 0)
       flag = SDialSetStopTime(sdial->mem, tEnd);
    else if (driver->SDIAL == -1)
       flag = rk4setStopTime(mem,tEnd);
	 

    io->errorFileWritten = 0;
    /* ----- call SUNDIALS to integrate to next time ----- */
    do {
      time->dt_max = tEnd - time->t;
      flameletPreiteration1D(driver,flame,chem,flow,time,mpi,time->t);
      if (driver->SDIAL == 0) {
        CVodeSetMaxStep(sdial->mem,time->dt_max);
        flag = CVode(sdial->mem, tEnd, sdial->Y, &time->t, sdial->itask);
      }
      else if (driver->SDIAL == 1)
        flag = IDASolve(sdial->mem, tEnd, &time->t, sdial->Y, sdial->Yp,
                        sdial->itask);
      else if (driver->SDIAL == -1) {
        rk4setMaxStep(mem,time->dt_max);
        flag = rk4solve(mem);
      }
      if (flag < 0) {
#ifdef ACCESS_MPI
  if (myid == mpi->host) 
#endif
        fprintf(stderr, "SUNDIALS execution failed: flag = %i\n", flag);
 
        if (io->dumpErrorFile && !io->errorFileWritten ) {
          flameletUpdateSoot1D(chem,flame,driver,io,flow,inert);
#ifdef ACCESS_MPI
          flameletAllgatherv(chem,flame,driver,io,mpi);
#endif /*ACCESS_MPI*/
          sprintf(errorFileName,"%s/%s",io->dataDir,"debug_.rif"); 
          createDataFileName(errorFileName, io, time, outputFileName);
          writeFlamelet(outputFileName, io, driver, time, flame, flow);
        }
#ifdef ACCESS_MPI
  if (myid == mpi->host) 
#endif
        fprintf(stderr, "Resetting CVODE...\n");

        /* restore old solution */
        for (j=0; j < flame->n1; j++) 
          for (i=0; i < flame->nVars; i++)   
	     flame->Y[j*flame->nVars+i] = fmax(0.0,flame->Y_m1[j*flame->nVars+i]);
 
        if (driver->output_w) {
           for (i = 0; i < chem->nRxnTotal; i++)
             chem->w_pre[i] = 0.0;
        } 
        /* reset time */
	time->t = time->t_n - time->dt;
        
        //flameletPreiteration1D(driver,flame,chem,flow,time,mpi,time->t);
        updateFlow1D(time->t, driver, time, flame, flow);
        /* reinitialise CVODE */
        MDUCresetODE(mem);

        if (io->dumpErrorFile && !io->errorFileWritten ) {
          flameletUpdateSoot1D(chem,flame,driver,io,flow,inert);
#ifdef ACCESS_MPI
          flameletAllgatherv(chem,flame,driver,io,mpi);
#endif /*ACCESS_MPI*/
          sprintf(errorFileName,"%s/%s",io->dataDir,"debug2_.rif"); 
          createDataFileName(errorFileName, io, time, outputFileName);
          writeFlamelet(outputFileName, io, driver, time, flame, flow);
	  io->errorFileWritten = 1;
        }
	
      }

      /* get solver data about last time step */
      if (driver->SDIAL >= 0) {
         SDialGetNumSteps(sdial->mem, &nst);
         SDialGetLastStep(sdial->mem, &dtPre);
         sdial->totalIt += nst; 
      }
      else if (driver->SDIAL == -1) {
         nst   = rk4getStep(mem);
         dtPre = rk4getStepSize(mem);
 
         if (nst%io->printFreq == 0)
#ifdef ACCESS_MPI
      	    if (myid == mpi->host)
#endif /* ACCESS_MPI */
		if (rk->iCrit < chem->nSpecies)
 	           fprintf(stdout,"Time step limiting: %s on grid point %i\n", 
                     chem->species[rk->iCrit], rk->jCrit);
	        else
 	           fprintf(stdout,"Time step limiting: Var %i on grid point %i\n", 
                     rk->iCrit, rk->jCrit);
      }
      
      /* update sources */
      if (io->outputHR || driver->output_tsrc) {
        for (j=0; j<flame->n1; j++) {
          if (io->outputHR)
            for (i = 0; i < chem->nSpecSolve; i++)
              flame->chemSrcY[j*flame->nVars+i] +=
                flame->chemSrcY_tmp[j*flame->nVars+i]*dtPre;
          
          if (driver->output_tsrc) {
            flame->chem_src[j] += flame->chem_temp[j]*dtPre;
            flame->conv_src[j] += flame->conv_temp[j]*dtPre;
            flame->diff_src[j] += flame->diff_temp[j]*dtPre;
          }
        }
      }
      if (driver->output_w && flag >= 0) {
         for (i = 0; i < chem->nRxnTotal; i++)
           chem->w_pre[i] += chem->w_tmp[i]*dtPre;
      }
    
      if (engine->spark) {
        engine->qPrem += engine->qdotPrem*dtPre;
      }
      
      /* print monitor */
#ifdef ACCESS_MPI
      if (nst%io->printFreq == 0 && myid == mpi->host)
#else
      time->dt_max = time->dt_max_;
      if (nst%io->printFreq == 0)
#endif /*ACCESS_MPI */
#ifndef RK
        printMonitor(nst, time->t, dtPre, io, sdial->Ydata, flame,time->dt_max);
#ifdef ACCESS_MPI
      
      // copy current solution to flamelet data struct and communicate
      Y   = NV_DATA_P(sdial->Y);
      bCommFlamelet1D(Y, sdial);
#endif /* ACCESS_MPI */
#else
        printMonitor(nst, time->t, dtPre, io, flame->Y, flame,time->dt_max);
#endif /* RK */
      
    } while (flag != sdial->iret);
    
  } /* end solve */

//  if (!inert) {
//  if (driver->output_w) {
//      realtype tmp;
//      flag = CVodeGetQuad(sdial->mem, &tmp, sdial->Q);
#ifdef ACCESS_MPI
//      MPI_Allreduce(sdial->Qdata,chem->w_int, chem->nRxnTotal, MPI_DOUBLE_PRECISION, MPI_SUM, mpi->comm);
#endif
//  } /* output_w */
//  } /* inert */
  flameletUpdateSoot1D(chem,flame,driver,io,flow,inert);
#ifdef ACCESS_MPI
  flameletAllgatherv(chem,flame,driver,io,mpi);
#endif /*ACCESS_MPI*/

  for (j=0;j<flame->n1;j++) {
    flame->sumYi[j] = 0.0;
    for (i=0;i<chem->nSpecies;i++)
      flame->sumYi[j] += flame->Y[j*flame->nVars+i];
  }

  /* update normalisation for next time step */
  flow->P_m1       = flow->P;
  flow->P_pre      = flow->P;
  flow->phi1min_m1 = flow->phi1min;
  flow->phi1max_m1 = flow->phi1max;
  copyArray(flame->n1, flow->chi1, flow->chi1_m1);

  /* update flow time */
  time->t_m1   = time->t;
  time->dt_pre = dt;
  time->dt_cvode_pre = dtPre;
  time->nt++;

  if (driver->ENGINE) {
    if (driver->PRES_MODE == 3)      /* implicit pressure */
      flow->P = flame->Y[flame->sP];
    else if (driver->PRES_MODE == 4) /*  explicit pressure */
      flow->P = calcEnginePressure(time->t, engine, flow, flame);
    
    /* compute the premixed */
    if (engine->spark) {
      /* engine->qdotPrem = (engine->qPrem - engine->qPrem_m1)/time->dt_pre; */
      /* engine->qPrem_m1 = engine->qPrem; */
    }
    
    writeEngineMonitor(time->t_m1, io, chem, engine, flame);
  }
  
  if (driver->ENERGY) {
    /* copy most current enthalpy from solution array for next time step */
    extractVarArray(flame->nVars,flame->n1,1,flame->sT,flame->Y,flame->H);
    
    /* make sure the temperature is up-to-date */
    for (j=0; j<flame->n1; j++) {
      k = j*flame->nVars;
      flame->T[j] = calcTfromH(chem->nSpecies, chem->nSpecSolve,
                               flame->Y[k+flame->sT], &flame->Y[k], 1000, chem);
    }
  }
  else { /* compute the total enthalpy if output to file */
    if (io->outputHtot)
      calcTotalEnthalpy(flame->n1,flame->n2,flame->nVars,flame->sT,flame->Y,
                        flame->H,chem);
  }
 
  /* update production rates */
  if (driver->output_w) {
     for (i = 0; i < chem->nRxnTotal; i++)
       chem->w_int_pre[i] += chem->w_pre[i];
  }
#ifdef ACCESS_MPI
  if (myid == mpi->host) {
#else
  {
#endif /* ACCESS_MPI */
    /* finish averaging of heat release rate over time step */
    if (io->outputHR)
      MDUCcalcHRrate(mem);
    
    /* write solution file at this step (change file name here if desired) */
    if (flame->nDim > 0) {
      if (io->freqType == 0 && time->nt%io->writeFreq == 0) {
        createDataFileName(io->outputData, io, time, outputFileName);
        writeFlamelet(outputFileName, io, driver, time, flame, flow);
        if (driver->output_w) {
        sprintf(rateFile,"%s/%s",io->dataDir,"rif_rates_.rates"); 
        createDataFileName(rateFile, io, time, outputFileName);
        writeRates(outputFileName, chem);
      }
      }
      else if (io->freqType == 1 && (fmod(time->t,io->tFreq) < 0.5*time->dt ||
                                     time->t > io->tWrite+io->tFreq) ) {
        createDataFileName(io->outputData,io,time,outputFileName);
        writeFlamelet(outputFileName,io,driver,time,flame,flow);
        if (driver->output_w) {
        sprintf(rateFile,"%s/%s",io->dataDir,"rif_rates_.rates"); 
        createDataFileName(rateFile, io, time, outputFileName);
        writeRates(outputFileName, chem);
        }
        io->tWrite = time->t;
      }
    }
  }
}
/* -------------------------------------------------------------------------- */

int cvodeFlamelet1D(realtype t, N_Vector u, N_Vector udot, void *data)
{
  int ierr;
  realtype *Y, *RHS;
  /* assign pointers to beginning of Sundials vector memory */
#ifdef ACCESS_MPI
  Y   = NV_DATA_P(u);
  RHS = NV_DATA_P(udot);

  /* communicate data from neighbouring processors */
  bCommFlamelet1D(Y, data);

  sdialStruct *sdial  = (sdialStruct*) data;
  flameStruct *flame  = sdial->args[3];

  /* use the data now copied/communicated in flame->Y for calculating RHS */
  Y = flame->Y;
#else /* SERIAL */
  Y   = NV_DATA_S(u);
  RHS = NV_DATA_S(udot);
#endif /* ACCESS_MPI */

  /* compute the current flamelet RHS */
  ierr = flameletRHS_1D(t, Y, RHS, data);

  return 0;
}

int cvodeFlameletQ1D(realtype t, N_Vector u, N_Vector udot, void *data)
{
  realtype *RHS,*Y;
  sdialStruct *sdial  = (sdialStruct*) data;
  /* assign pointers to beginning of Sundials vector memory */
#ifdef ACCESS_MPI
  RHS = NV_DATA_P(udot);
  Y = NV_DATA_P(u);
#else /* SERIAL */
  RHS = NV_DATA_S(udot);
  Y = NV_DATA_S(u);
#endif /* ACCESS_MPI */

  driveStruct  *driver = (driveStruct*)  sdial->args[0];
  chemStruct   *chem   = (chemStruct*)   sdial->args[1];

  /* compute the current flamelet RHS */
  flameletQ_RHS_1D(sdial, driver, chem, RHS);

  return 0;
}
/* -------------------------------------------------------------------------- */

#ifdef ACCESS_MPI
#ifdef SUNDIALS_FROM_2_7
int cvodeFlamelet1D_local(long int Nlocal, realtype t, N_Vector u, N_Vector udot,
                          void *data)
#else
#ifdef SUNDIALS_FROM_3_0
int cvodeFlamelet1D_local(sunindextype Nlocal, realtype t, N_Vector u, N_Vector udot,
                          void *data)
#else /* SUNDIALS LOWER THAN 2.7 */
int cvodeFlamelet1D_local(int Nlocal, realtype t, N_Vector u, N_Vector udot,
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
  copyArray(mpi->jCount[mpi->myid], Y, &flame->Y[mpi->jStart[mpi->myid]]);
  Y = flame->Y;

  /* compute the current flamelet RHS */
  flameletRHS_1D(t,Y, RHS, data);
  return 0;
}
/* -------------------------------------------------------------------------- */

int bCommFlamelet1D(realtype *Y, void *data)
{
  sdialStruct *sdial  = (sdialStruct*) data;
  flameStruct *flame  = sdial->args[3];
  mpiStruct   *mpi    = sdial->args[7];

  int i, myid=mpi->myid;
  MPI_Status  status;
  MPI_Request request[4];
  MPI_Request rrequest[4];

  /* copy current solution to flamelet memory for RHS computation */
  copyArray(mpi->jCount[myid], Y, &flame->Y[mpi->jStart[myid]]);
#ifndef BLOCK
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
#else
  for (i=0;i<mpi->nSRcomm;i++) {
    if (mpi->myid == mpi->srSender[i]) /* send */
       MPI_Send(&flame->Y[mpi->srIndex[i]],mpi->srCount[i],MPI_DOUBLE_PRECISION,
                mpi->srRecver[i],mpi->srComm[i],mpi->comm);
    else /* receive */
       MPI_Recv(&flame->Y[mpi->srIndex[i]],mpi->srCount[i],MPI_DOUBLE_PRECISION,
                mpi->srSender[i],mpi->srComm[i],mpi->comm,&status);
  }
#endif /* NONBLOCK */
  return 0;
}
/* -------------------------------------------------------------------------- */
#endif /* ACCESS_MPI */

int idaFlamelet1D(realtype t, N_Vector y, N_Vector yp, N_Vector r, void *data)
{
  int      i,j,k,ierr;
  double   dt;
  realtype *Y, *Yprime, *res, *RHS;

  /* dereference the necessary structs */
  sdialStruct  *sdial  = (sdialStruct*)  data;
  driveStruct  *driver = (driveStruct*)  sdial->args[0];
  flowStruct   *flow   = (flowStruct*)   sdial->args[2];
  flameStruct  *flame  = (flameStruct*)  sdial->args[3];
  timeStruct   *time   = (timeStruct*)   sdial->args[4];

  /* assign pointers to beginning of Sundials vector memory */
  Y      = NV_DATA_S(y);  /* current value of dependent variable vector y(t) */
  Yprime = NV_DATA_S(yp); /* current value of y'(t) */
  res    = NV_DATA_S(r);  /* output residual vector */
 
  /* assign memory for the right-hand-side */
  RHS = (double*) calloc(sdial->nEq, sizeof(double));

  /* get current location from start of time step */
  dt = t - time->t_m1;

  /* compute the current flamelet RHS */
  ierr = flameletRHS_1D(t, Y, RHS, data);

  /* update enthalpy equation if necessary */
  if (driver->ENERGY) {
    for (j=0; j<flame->n1; j++) {
      k = j*flame->nVars;
      Yprime[k+flame->sT] = flame->H[j] + 1.0/flame->rho[j]*flow->Dt_P*dt;
    }
  }

  /* compute the residual */
  for (i=0; i<sdial->nEq; i++)
    res[i] = Yprime[i] - RHS[i];

  free(RHS);
  return ierr;
}
/* -------------------------------------------------------------------------- */

int flameletRHS_1D(double t, double *Y, double *RHS, void *args)
{
  int    i,j,k,iCAD;
  double *T; /* local temperature */
  double P,m,dmdt,vol_u,dVdt,rho,qdot;
  double CAD,CAD_rel;
  int    nanFound;
  nanFound = 0;

  /* dereference the necessary structs */
  sdialStruct  *sdial  = (sdialStruct*)  args;
  driveStruct  *driver = (driveStruct*)  sdial->args[0];
  chemStruct   *chem   = (chemStruct*)   sdial->args[1];
  flowStruct   *flow   = (flowStruct*)   sdial->args[2];
  flameStruct  *flame  = (flameStruct*)  sdial->args[3];
  timeStruct   *time   = (timeStruct*)   sdial->args[4];
  engineStruct *engine = (engineStruct*) sdial->args[5];

  /* ----- update quantities before computing derivatives ----- */
  /* get temperature based on form of energy equation */
  if (driver->ENERGY) {
    T = flame->T;
    /* update to most current temperature and enthalpy */
    for (j=0; j<flame->n1; j++) {
      k = j*flame->nVars;
      T[j] = calcTfromH(chem->nSpecies, chem->nSpecSolve,
                        Y[k+flame->sT], &Y[k], 1000, chem);
    } /* end loop over grid */
  }
  else { /* temperature is part of solution */
    T = (double*) calloc(flame->n1, sizeof(double));
    extractVarArray(flame->nVars,flame->n1,1,flame->sT,Y,T);
  }

  /* ----- check for unphysical values first ---- */
  for (j=flame->j0; j<=flame->jN; j++) {
        if ((T[j]<100.0) || (T[j]>5000.0)) {
#ifdef DEBUG
           fprintf(stderr,"Clipping unphysical temperature! %g - %i\n", T[j], j);
#endif /* DEBUG */ 
	   T[j] = fmin(fmax(100.0,T[j]),5000.0);
        }
        if (isnan(T[j]) || isinf(T[j])) {
#ifdef DEBUG
           fprintf(stderr,"NaN or inf temperature!\n");
#endif /* DEBUG */
	   nanFound = 1;
           for (i=0;i<flame->nVars;i++)
              RHS[(j-flame->j0)*flame->nVars+i] = 1.0/0.0;
        }

        flame->sumYi[j] = 0.0;
        for (i=0;i<chem->nSpecies;i++)
	    flame->sumYi[j] += Y[j*flame->nVars+i];

        if (flame->sumYi[j] <= 0.0) {
           for (i=0;i<flame->nVars;i++)
              RHS[(j-flame->j0)*flame->nVars+i] = 1.0/0.0;

	   nanFound = 1; 
        }
  }

  if (nanFound == 1) return 0;


  /* compute a coupled equation for the scalar dissipation rate */
#ifndef CHIPREITER
  if (driver->CHI_MODE == 2) {
    flameletChi_RHS_1D(driver,chem,flame,flow,T,Y,RHS);
    extractVarArray(flame->nVars,flame->n1,1,flame->sChi,Y,flame->chiA);
  }
  else /* interpolate flow variables (chi) for current t */
    updateFlow1D(t, driver, time, flame, flow);
#endif /* CHIPREITER */


  if (driver->ENGINE) {
    /* ----- pre-compute some engine quantities ----- */
    CAD     = time_to_cad(t,engine);
    CAD_rel = fmod(CAD, 720);

    if (driver->PRES_MODE == 3)
      P = Y[flame->sP];
    else if (driver->PRES_MODE == 4)
      P = flow->P;

    /* volume and rate */
    vol_u = engine->vol = calcCylVol(time_to_cad(t,engine),engine);
    iCAD  = gridLocate(engine->nCAD,engine->CAD_vec,CAD_rel);
    dVdt  = linInterp(CAD_rel,engine->CAD_vec[iCAD],engine->dVdt[iCAD],
                      engine->CAD_vec[iCAD+1],engine->dVdt[iCAD+1]);

    if (chem->wiebeMech) {
      wiebeStruct *wiebe = (wiebeStruct*) chem->args[0];
      wiebe->t = CAD_rel;
      /* ensure relative crank angle is defined centred around combustion at
         0CAD, i.e. -360 < CAD_rel < 360 */
      if (CAD_rel > 360.0)
        wiebe->t -= 720;
      wiebe->V    = vol_u;
      wiebe->dVdt = dVdt;
    }
  }

  /* ----- compute species derivatives ----- */
  flameletY_RHS_1D(sdial, driver, chem, flame, flow, T, Y, RHS);

  /* ----- compute energy equation ----- */
  if (driver->ENERGY) {
    /* update the energy based on the species current species and enthalpy
       (calculated in flameletY_RHS_1D) */
    for (j=flame->j0; j<=flame->jN; j++)
      RHS[(j-flame->j0)*flame->nVars+flame->sT] = flame->srcH[j] + flame->visc_src[j];
  }
  else {
    flameletT_RHS_1D(driver, chem, flame, flow, Y, RHS);
    free(T);
  }

  if (chem->sootMode == 1) {
    /* apply the computed soot moment source terms to the equations */
    for (j=flame->j0; j<=flame->jN; j++) {
      k = (j-flame->j0)*flame->nVars + flame->sSoot;
      for (i=0; i<chem->nSootMom; i++) {
        RHS[k+i] = 0.0;
        if (chem->sootMomOn[i]) {
	  RHS[k+i] = (flame->sootMomSrc[j*chem->nSootMom+i]+
                      flame->sootMomTrans[j*chem->nSootMom+i])/flame->rho[j];

        }
      }
    }
  }

  if (driver->ENGINE) {
     /* check mass flow */
    calcEngineRHS_0D(t, time, chem, flow, flame, engine, Y, RHS);

    /* reconstruct the mass terms */
    m = dmdt = 0.0;
    for (i=0; i<3; i++) {
      m    += Y[flame->sM + i];
      dmdt += RHS[flame->sM + i];
    }

    /* calculate the heat losses */
    if (driver->Q_LOSS) {
      for (j=0; j<flame->n1; j++) {
        qdot = calcWoschniQdot(t,P,T[j],flame->gamma[j],engine);
        
        /* convert to correct units [W/(m-K-kg)] and add to source */
        qdot *= flame->R[j]*T[j]/(P*engine->vol);
        RHS[j*flame->nVars+flame->sT] += qdot/flame->cpmix[j];
      }
    }
      
    /* ----- compute the pressure rate ----- */
    if (engine->spark) {
      calcEngineSpark(t,driver,time,chem,flow,flame,engine,Y,RHS);
    }

    if (driver->PRES_MODE == 3) { /* implicit pressure treatment*/
      double dMw_dt = 0.0, dTdt_mean = 0.0;

      double Rmean = 0.0, Tmean = 0.0;
      for (j=0; j<flame->n1; j++) {
        Rmean += flame->R[j]*flow->pdf[j];
        Tmean += T[j]*flow->pdf[j];

        for (i=0; i<chem->nSpecies; i++)
          dMw_dt+=RHS[j*flame->nVars+i]/chem->W[j*flame->nVars+i]*flow->pdf[j];
        
        dTdt_mean += RHS[j*flame->nVars+flame->sT]*flow->pdf[j];
      }
      
      flow->meanRho = rho = m/vol_u;
      RHS[flame->sP] =  Rmean*Tmean/vol_u*(dmdt - rho*dVdt)
                      + rho*(Rmean*dTdt_mean + T[j]*flame->RGAS*dMw_dt);
    }
  } /* end engine */
  /* ----- check for unphysical values again ---- */
  int error=0;
  for (j=flame->j0; j<=flame->jN; j++) {
     for (i=0;i<flame->nVars;i++) {
 	if (isnan(RHS[(j-flame->j0)*flame->nVars+i])||
	    isinf(RHS[(j-flame->j0)*flame->nVars+i])) {
#ifdef DEBUG
 	       fprintf(stderr,"NaN or inf in RHS\n");
               fprintf(stderr,"Var:%i gridPoint:%i value:%15.6e\n",i,j,Y[j*flame->nVars+i]);
#endif /* DEBUG */
	       error += 1;
            }
        }
     if (error > 0 ) {
#ifdef DEBUG
        for (i=0;i<flame->nVars;i++) 
            fprintf(stderr,"Var:%i value:%15.6e\n",i,Y[j*flame->nVars+i]);

        fprintf(stderr,"rhoSigma: %15.6e\n",flame->rhoSigma[j]);
        fprintf(stderr,"mu      : %15.6e\n",flame->mu[j]);
        fprintf(stderr,"D       : %15.6e\n",flame->D[j]);
#endif /* DEBUG */

      //  return 1;
     }
  }
  return 0;
}
/* -------------------------------------------------------------------------- */

void flameletY_RHS_1D(sdialStruct *sdial, driveStruct *driver, chemStruct *chem, flameStruct *flame,
                      flowStruct *flow, double *T, double *Y, double *RHS)
{
  int    i, j, k, jst, kM1, kP1, kL, iS;
  int 	 s,sM1,sP1,sM2;
  int    nSpecies;
  double Dx_Y, Dxx_Y;       /* spatial derivatives */
  double Dx_M, Dx_psiChi, Dx_psiD;
  double src, diff, conv;   /* source, diffusion, and convection terms */
  double vH;                /* additional velocity for enthalpy equation */
  double P, dPdt;           /* pressure terms */
  double vMxmyp,vA,vT,vMxy,vDL;  /* convection velocity for soot moments */
  double Dxx_T,Dx_T, Dx_xi, Dx_chi, Dx_D;  /* Addional gradients required for thermophoresis term */
  double Dxx_Mxmyp, Dxx_N0,Dx_rho,Dx_Mxmyp, Dxx_Mxy;
  double S0,V0,V0xS0y,DN0_DZ,N0,M10,M01,dp02,dpL2;
  int kMxmyp; 
  double Df = 1.8; /* Fractal dimension of soot aggregates Blanquart 2009 */
  DN0_DZ   = driver->sootDiffusivity; /* Diffusivity of small soot particles divided by DZ */
  V0       = 2.0*NuclVol;
  S0       = pow(V0,2.0/3.0);
  nSpecies = chem->nSpecSolve; /* solve only non-steady state species */

  /* set default pressure properties */
  P    = flow->P;
  dPdt = flow->Dt_P;

  if (driver->output_w) {
    for (i = 0; i < chem->nRxnTotal; i++)
        chem->w_tmp[i] = 0.0;
  }
#ifndef PREITER
  /* Update convection velocities */
  flameletConvVel1D(driver,chem,flame,flow,Y,T);
#else
  /* loop over mesh points to solve species equations */
  for (j=flame->j0m; j<=flame->jNp; j++) {
    jst = stencil_x(j, flame); /* start of stencil */
    k   = j*flame->nVars;      /* array index for start of solution */

    /* ----- compute local properties ----- */
    flame->Mw[j]     = calcMixtureW(nSpecies, &Y[k], chem->W);
    /* compute mixture density */
    flame->rho[j]    = flame->Mw[j]*P/(flame->RGAS*T[j]);

    flame->R[j]      = flame->RGAS/flame->Mw[j];
    /* get enthalpy and specific heat from chemical mechanism */
    chemComputeThermoData(chem->h, chem->cp, T[j], chem);
    /* get specific heat of the mixture */
    flame->cpmix[j] = calcMixtureCp(nSpecies, chem->cp, &Y[k]);

    if (driver->detailedDiff) {
       /* Update species Lewis numbers here */
       for (i=0;i<chem->nSpecSolve;i++) {
          k = j*chem->nSpecSolve + i;
          flame->LeZ_LeI[k] = 1.0;
          // flame->LeZ_LeI[k] = calcSpecDiffusivity(nSpecies, T[j], &Y[j*flame->nVars],chem);
       }
    }
  }
#endif

  /* loop over mesh points to solve species equations */
  for (j=flame->j0; j<=flame->jN; j++) {
#ifdef COUNTERFLOW
    if (!(j == 0 || j == flame->n1-1)) 
#endif /* COUNTERFLOW */
    {
    kM1 = (j-1)*flame->nVars;
    kP1 = (j+1)*flame->nVars;
    k   = j*flame->nVars;      /* array index for start of solution */
    kL  = (j-flame->j0)*flame->nVars;
    jst = stencil_x(j, flame); /* start of stencil */
    
    /* get enthalpy and specific heat from chemical mechanism */
    chemComputeThermoData(chem->h, chem->cp, T[j], chem);

    /* get specific heat of the mixture */
    //flame->cpmix[j] = calcMixtureCp(nSpecies, chem->cp, &Y[k]);
    flame->gamma[j] = flame->cpmix[j]/(flame->cpmix[j] - flame->R[j]);

#ifndef FM_VS_MDUC
    /* compute the mixture viscosity */
    if (driver->VISC_MODE)
       if (chem->readTrans) 
          flame->mu[j] = calcMixViscosity(nSpecies, T[j], &Y[k], chem);
#endif /* FM_VS_MDUC */

    /* compute and store concentrations of all species at local point */

    for (i = 0; i < nSpecies; i++)
      chem->c[i] = flame->rho[j]*Y[k+i]/chem->W[i];
    
    /* get rates from mechanism */
    /* Deactivate chemical reactions for rich mixtures that are not present in the domain*/
    int reactive = 1;
    if (flow->chi1[j] <= 0.0 && flow->chi1_m1[j] <= 0.0) 
	  reactive = 0;

    if (reactive > 0) {
       chemComputeProductionRates(chem->cdot,chem->w,chem->K,chem->c,chem->M,
                               T[j],P,chem);
    }
    else {
      /* zero out reaction rates */
         for (i = 0; i < nSpecies; i++)
           chem->cdot[i] = 0.0;
         for (i = 0; i < chem->nRxn; i++)
           chem->w[i] = 0.0;
    }

    /* add any sources due to soot model */
    if (chem->soot > 0) { /* FIXME: default to HMOM for now */
      /* copy the species necessary for soot reactions */
      for (iS=0; iS<chem->nSootSp; iS++)
        chem->sootRxnConc[iS] = chem->c[chem->iSootSp[iS]];

      /* convert from M/rho to M  for source term evaluation*/
      k = j*flame->nVars + flame->sSoot;
      for (i=0; i<chem->nSootMom; i++) 
        flame->sootMom[j*chem->nSootMom+i] = Y[k+i]*flame->rho[j];

      /* get the soot source terms */
      if (chem->sootMode <= 3) {
         if (reactive>0) {
            sootHMOM(T[j],flame->mu[j],flame->rho[j],flame->Mw[j],
                     chem->sootRxnConc,&flame->sootMom[j*chem->nSootMom],
                     &flame->sootMomSrc[j*chem->nSootMom],chem->sootRxnRates,
                     &flame->rhoDotSoot[j],&flame->sootMomSrcPart[j*chem->nSootMom*chem->nSootSrc],
                     chem->sootSrcOn,driver->SABSTOL);
                   
            /* ----- add the rates to the source terms ----- */
            /* recompute rho dot */
            flame->rhoDotSoot[j] = 0.0;

            for (i=0; i<nPAH; i++) {
             chem->cdot[chem->sPAH[i]] -= 2.0*chem->sootRxnRates[i]; 
             flame->rhoDotSoot[j] -= 2.0* chem->sootRxnRates[i] * chem->W[chem->sPAH[i]];
            } 

            chem->cdot[chem->sOH]   -= chem->sootRxnRates[nPAH] + chem->sootRxnRates[nPAH+5];
            chem->cdot[chem->sH]    += chem->sootRxnRates[nPAH+2] - chem->sootRxnRates[nPAH+1];
            chem->cdot[chem->sC2H2] -= chem->sootRxnRates[nPAH+3];
            chem->cdot[chem->sO2]   -= chem->sootRxnRates[nPAH+4];
            chem->cdot[chem->sH2O]  += chem->sootRxnRates[nPAH];
            chem->cdot[chem->sH2]   += chem->sootRxnRates[nPAH+1];
            chem->cdot[chem->sCO]   += 2.0*chem->sootRxnRates[nPAH+4]+chem->sootRxnRates[nPAH+5];

            flame->rhoDotSoot[j]  -= (chem->sootRxnRates[nPAH] + chem->sootRxnRates[nPAH+5])* chem->W[chem->sOH];
            flame->rhoDotSoot[j]  += (chem->sootRxnRates[nPAH+2] - chem->sootRxnRates[nPAH+1])* chem->W[chem->sH];
            flame->rhoDotSoot[j]  -= (chem->sootRxnRates[nPAH+3])* chem->W[chem->sC2H2];
            flame->rhoDotSoot[j]  -= (chem->sootRxnRates[nPAH+4])* chem->W[chem->sO2];
            flame->rhoDotSoot[j]  += (chem->sootRxnRates[nPAH])* chem->W[chem->sH2O];
            flame->rhoDotSoot[j]  += (chem->sootRxnRates[nPAH+1])* chem->W[chem->sH2];
            flame->rhoDotSoot[j]  += (2.0*chem->sootRxnRates[nPAH+4]+chem->sootRxnRates[nPAH+5])* chem->W[chem->sCO];

            if (driver->output_w) { 
               for (i=0; i<chem->nSootRxn;i++) 
	           chem->w_tmp[i+chem->sSoot] += flame->gPDF[j]/flame->rho[j]*chem->sootRxnRates[i];
            }
         }
         else {
           for (i=0; i<chem->nSootMom; i++) 
           	flame->sootMomSrc[j*chem->nSootMom+i] = 0.0;
         }
      } /* sootMode <= 3 */
    } /* soot */

    /* loop over species */
    for (i=0; i<nSpecies; i++) {
      /* calculate the gradients in normalised space */
      Dxx_Y = grad_xx(jst, 0, j, i, flame->nVars, Y, flame);
      if (j == 0)
        Dx_Y = (Y[kP1+i]-Y[k+i])/(flame->x1[j+1]-flame->x1[j]);
      else if (j == flame->n1-1)
        Dx_Y = (Y[k+i]-Y[kM1+i])/(flame->x1[j]-flame->x1[j-1]);
      else {
        //if (flame->v1[j] + vH >= 0.0)
        if (flame->v1[j] <= 0.0)
          Dx_Y = (Y[k+i]-Y[kM1+i])/(flame->x1[j]-flame->x1[j-1]);
        else
          Dx_Y = (Y[kP1+i]-Y[k+i])/(flame->x1[j+1]-flame->x1[j]);
      }

      /* ----- compute diffusion term -----*/
      diff = 0.5*flame->chiA[j]*Dxx_Y;
      
      /* ----- compute convection term ----- */
      conv = (flame->v1[j])*Dx_Y; 
      // TODO: convection due to mass source
      // TODO: fuel source term from rhoDotSpray
      flame->convCoeff[j] = flame->v1[j];
      
      /* ----- compute chemical source term ----- */
      src  = chem->cdot[i]*chem->W[i]/flame->rho[j];
      
      /* store chemical source term for heat release computation */
      flame->chemSrcY_tmp[k+i] = src;

      /* ----- soot source term ----- */
#ifdef RHODOT
      src -= (flame->rhoDotSpray[j] + flame->rhoDotSoot[j]) * Y[k+i] / flame->rho[j];
#else
      src -= (flame->rhoDotSoot[j]) * Y[k+i] / flame->rho[j];
#endif /* RHODOT */

      /* DEBUG: force convective velocity to be zero at the boundary */
      if ((j == 0) && ((flame->v1[j]) > 0.0))
        conv = 0.0;
      else if  ((j == flame->n1-1) && ((flame->v1[j]) < 0.0))
        conv = 0.0;
      
      /* ----- update RHS of the species equations ----- */
#ifndef SOOTONLY
      RHS[kL+i] = driver->DIFF_MODE[0]*diff - driver->CONV_MODE*conv + src;
#else
      RHS[kL+i] = 0.0;
#endif /* SOOTONLY */
    } /* end loop over non-steady state species */

    /* remove soot reactions from heat release */
    for (i=0; i<nPAH; i++) {
     chem->cdot[chem->sPAH[i]] += 2.0*chem->sootRxnRates[i]; 
    } 

    chem->cdot[chem->sOH]   += chem->sootRxnRates[nPAH] + chem->sootRxnRates[nPAH+5];
    chem->cdot[chem->sH]    -= chem->sootRxnRates[nPAH+2] - chem->sootRxnRates[nPAH+1];
    chem->cdot[chem->sC2H2] += chem->sootRxnRates[nPAH+3];
    chem->cdot[chem->sO2]   += chem->sootRxnRates[nPAH+4];
    chem->cdot[chem->sH2O]  -= chem->sootRxnRates[nPAH];
    chem->cdot[chem->sH2]   -= chem->sootRxnRates[nPAH+1];
    chem->cdot[chem->sCO]   -= 2.0*chem->sootRxnRates[nPAH+4]+chem->sootRxnRates[nPAH+5];

    flame->srcH[j] =  0.0;

    /* loop over species */
    for (i=0; i<nSpecies; i++) {
      /* ----- compute chemical source term ----- */
      src  = chem->cdot[i]*chem->W[i]/flame->rho[j];
      
      /* store chemical source term for heat release computation */
      flame->chemSrcY_tmp[k+i] = src;

      /* ----- compute summation terms required for energy equation ----- */
      if (driver->ENERGY) /* solve enthalpy */
        flame->srcH[j] += chem->h[i]*Y[k+i];
      else { /* solve temperature */
        flame->srcH[j]  += chem->h[i]*src;   /* for chemical source term */
      }
    }

    /* ----- copy production rates if needed ----- */
    if (driver->output_w) {
      for (i = 0; i < chem->nRxn; i++) {
          chem->w_tmp[i] += flame->gPDF[j]/flame->rho[j]*chem->w[i];
      }
    }
  }
#ifdef COUNTERFLOW
  else {
    for (i=0;i<flame->nVars;i++) {
        kL  = (j-flame->j0)*flame->nVars;
        RHS[kL+i] = 0.0;
    }
    if (chem->sootMode < 3) 
       for (i=0; i<chem->nSootMom; i++) 
 	flame->sootMomSrc[j*chem->nSootMom+i] = 0.0;
  }
#endif
  } /* end mesh loop for species equations */

  if (chem->sootMode == 1) {
    /* if implicit SOOT equation, make sure moments are up-to-date
       NOTE: the solved equations are not used directly */
    for (j=flame->j0m; j<=flame->jNp; j++) {
      k = j*flame->nVars + flame->sSoot;
      for (i=0; i<chem->nSootMom; i++) {
        flame->sootMom[j*chem->nSootMom+i] = Y[k+i];
      }
      kMxmyp = j* (chem->nSootMom-1);
#ifndef NFREEPATH
      flame->sootMomxmyp[kMxmyp]   = FracMomDensity(4.0/Df-2.0,-6.0/Df+2.0,&Y[k],flame->rho[j])/flame->rho[j];
      flame->sootMomxmyp[kMxmyp+1] = FracMomDensity(4.0/Df-1.0,-6.0/Df+2.0,&Y[k],flame->rho[j])/flame->rho[j];
      flame->sootMomxmyp[kMxmyp+2] = FracMomDensity(4.0/Df-2.0,-6.0/Df+3.0,&Y[k],flame->rho[j])/flame->rho[j];
#else
      flame->sootMomxmyp[kMxmyp]   = FracMomDensity(2.0/Df-1.0,-3.0/Df+1.0,&Y[k],flame->rho[j])/flame->rho[j];
      flame->sootMomxmyp[kMxmyp+1] = FracMomDensity(2.0/Df    ,-3.0/Df+1.0,&Y[k],flame->rho[j])/flame->rho[j];
      flame->sootMomxmyp[kMxmyp+2] = FracMomDensity(2.0/Df-1.0,-3.0/Df+2.0,&Y[k],flame->rho[j])/flame->rho[j];
#endif /*NFREEPATH */
    }
  }

  if (chem->soot>0 && chem->sootMode <= 2) {
    /* loop over mesh points to solve species equations */
    for (j=flame->j0; j<=flame->jN; j++) {
      if (j != 0 && j !=flame->n1-1) {
      jst = stencil_x(j, flame); /* start of stencil */

      /* Compute required gradients */
      Dx_chi = grad_x(jst, 0, j, 0, 1 , flame->chiA, flame);
      Dx_D   = grad_x(jst, 0, j, 0, 1 , flame->D, flame);
      //Dx_chi = (flame->chiA[j]-flame->chiA[j-1])/(flame->x1[j]-flame->x1[j-1]);
      //Dx_D   = (flame->D[j]-flame->D[j-1])/(flame->x1[j]-flame->x1[j-1]);

      if (driver->thermophoresis) {
        Dx_T  = grad_x(jst, 0, j, flame->sT, flame->nVars, Y, flame);
        Dxx_T = grad_xx(jst, 0, j, flame->sT, flame->nVars, Y, flame);
        Dx_xi = grad_x(jst, 0, j, 0, 1 , flame->xi, flame);
      }

      vMxy = flame->vMxy[j];
#ifndef NFREEPATH
      vMxmyp = flame->vMxmyp[j]*DN0_DZ*pow(S0,6.0/Df-2.0)*pow(V0,2.0-4.0/Df);
#else
      vMxmyp = flame->vMxmyp[j]*DN0_DZ*pow(S0,3.0/Df-1.0)*pow(V0,1.0-2.0/Df);
#endif /*NFREEPATH*/
      /* ---- compute soot diffusion in mixture fraciton space ---- */
      for (i=0; i<chem->nSootMom; i++) 
      {
         Dx_Mxmyp = 0.0;
         if (i == chem->nSootMom-1) {
           vMxy += flame->vMxmyp[j]*DN0_DZ;
           Dx_Mxmyp = 0.0;
         }
         else {
           if (driver->sootDiff) {
              s   = j    *(chem->nSootMom-1) + i;
              sP1 = (j+1)*(chem->nSootMom-1) + i;
              sM1 = (j-1)*(chem->nSootMom-1) + i;
              /* use a simple upwinding scheme */
              if (vMxmyp >= 0.0) {
                  Dx_Mxmyp  = (flame->sootMomxmyp[s] - flame->sootMomxmyp[sM1])
                  	/(flame->x1[j]-flame->x1[j-1]);
              }
              else {
                  Dx_Mxmyp  = (flame->sootMomxmyp[sP1] - flame->sootMomxmyp[s])
                  	/(flame->x1[j+1] - flame->x1[j]);
              }
           }
         }

         /* use a simple upwinding scheme */
         s   =  j   *chem->nSootMom + i;
         sP1 = (j+1)*chem->nSootMom + i;
         sM1 = (j-1)*chem->nSootMom + i;

         if (vMxy >= 0.0) {
             Dx_M  = (flame->sootMom[s] - flame->sootMom[sM1])
             	/(flame->x1[j]-flame->x1[j-1]);
         }
         else {
             Dx_M  = (flame->sootMom[sP1] - flame->sootMom[s])
             	/(flame->x1[j+1] - flame->x1[j]);
         }
 
         /* convective term */ 
         flame->sootMomTrans[s]  = -vMxy *Dx_M;
         if (i < chem->nSootMom-1) /* Mxy */
            flame->sootMomTrans[s]  -= vMxmyp*Dx_Mxmyp;

         flame->sootMomTransConv[s] = flame->sootMomTrans[s];
#ifdef RHODOT
         TODO     
#endif /* RHODOT */
         /* source terms due to thermophoresis */
         if (driver->thermophoresis) {
             flame->sootMomTrans[s] += flame->sootMom[s] * (
          		       flame->xi[j]*Dx_T*(Dx_chi + flame->chiA[j]/flame->D[j]*Dx_D) +
          		       2.0*flame->chiA[j]*Dx_xi *Dx_T +
          		       2.0*flame->chiA[j]*flame->xi[j]*Dxx_T);
         }
         if (driver->sootDiff) {
             if (i==chem->nSootMom-1) { /* N0 */
                /* soot diffustion contribution from large particles */
                Dxx_N0 = grad_xx(jst, 0, j, i, chem->nSootMom, flame->sootMom, flame);
                flame->sootMomTrans[s] += 0.5*flame->rho[j]*flame->chiA[j]*DN0_DZ*Dxx_N0;
                flame->sootMomTransDiff[s] = 0.5*flame->rho[j]*flame->chiA[j]*DN0_DZ*Dxx_N0;
             }
             else { /*Mxmyp*/
                /* soot diffustion contribution from large particles */
                int nSootMomxmyp = chem->nSootMom - 1;
                Dxx_Mxmyp = grad_xx(jst, 0, j, i, nSootMomxmyp, flame->sootMomxmyp, flame);
#ifndef NFREEPATH
                flame->sootMomTrans[s] += 
		     0.5*flame->rho[j]*flame->chiA[j]*DN0_DZ*pow(S0,6.0/Df-2.0)*pow(V0,2.0-4.0/Df)
                     *Dxx_Mxmyp;
                flame->sootMomTransDiff[s] = 
		     0.5*flame->rho[j]*flame->chiA[j]*DN0_DZ*pow(S0,6.0/Df-2.0)*pow(V0,2.0-4.0/Df)
                     *Dxx_Mxmyp;
#else
                flame->sootMomTrans[s] += 
		     0.5*flame->rho[j]*flame->chiA[j]*DN0_DZ*pow(S0,3.0/Df-1.0)*pow(V0,1.0-2.0/Df)
                     *Dxx_Mxmyp;
                flame->sootMomTransDiff[s] = 
		     0.5*flame->rho[j]*flame->chiA[j]*DN0_DZ*pow(S0,3.0/Df-1.0)*pow(V0,1.0-2.0/Df)
                     *Dxx_Mxmyp;
#endif
             }
         }
#ifdef DEBUG
         if (isnan(flame->sootMomTrans[s]) || isinf(flame->sootMomTrans[s])) {
            fprintf(stderr,"NaN on Soot moment source after transport! %i\n", i);
            fprintf(stderr,"D,rhoSigma,xi,vMxmyp,vMxy,sootMom: %g %g %g %g %g %g\n", 
          	flame->D[j], flame->rhoSigma[j], flame->xi[j], vMxmyp, vMxy ,flame->sootMom[s]);
         }
#endif /* DEBUG */
        } /* Soot moments */
      } 
      else { /* j == 1 || j == flame->n1 - 1 */
        for (i=0; i<chem->nSootMom; i++) {
           s = j*chem->nSootMom + i;
           flame->sootMomTrans[s]     = 0.0;
           flame->sootMomTransDiff[s] = 0.0;
           flame->sootMomTransConv[s] = 0.0;
        }
      }
    } /* end mesh loop for species equations */
  } /* soot mode <= 2 */
}
void flameletConvVel1D(driveStruct *driver, chemStruct *chem, flameStruct *flame, 
			      flowStruct *flow,double *Y, double *T)
{
  int    i, j, k,jst,kM1,kP1,kMxmyp;
  double vA,vT,vH,vDL,vMxmyp,vMxy;
  double Dx_T, Dx_xi, Dx_chi, Dx_D, Dx_rho,Dx_Y;
  double Dx_cpmix;
  double DN0_DZ,V0,S0;
  double Df = 1.8;
  DN0_DZ   = driver->sootDiffusivity; /* Diffusivity of small soot particles divided by DZ */
  V0       = 2.0*NuclVol;
  S0       = pow(V0,2.0/3.0);
  /* loop over mesh points to solve and ghost cells and update properties */
  for (j=flame->j0m; j<=flame->jNp; j++) {
    k   = j*flame->nVars;      /* array index for start of solution */

    /* ----- compute local properties ----- */
    flame->Mw[j]     = calcMixtureW(chem->nSpecies, &Y[k], chem->W);
    /* compute mixture density */
    flame->rho[j]    = flame->Mw[j]*flow->P/(flame->RGAS*T[j]);

    flame->R[j]      = flame->RGAS/flame->Mw[j];
    /* get enthalpy and specific heat from chemical mechanism */
    chemComputeThermoData(chem->h, chem->cp, T[j], chem);
    /* get specific heat of the mixture */
    flame->cpmix[j] = calcMixtureCp(chem->nSpecies, chem->cp, &Y[k]);

    if (driver->detailedDiff) {
       /* Update species Lewis numbers here */
       for (i=0;i<chem->nSpecSolve;i++) {
          k = j*chem->nSpecSolve + i;
          flame->LeZ_LeI[k] = 1.0;
       }
    }

    /* Evaluate species terms */
    flame->sumCp[j] = 0.0;   /* for chemical source term */
    kM1 = (j-1)*flame->nVars;
    k   =  j   *flame->nVars;
    kP1 = (j+1)*flame->nVars;

    for (i=0; i<chem->nSpecies; i++) {
      /* calculate the gradients in normalised space */
      if (j == 0)
        Dx_Y = (Y[kP1+i]-Y[k+i])/(flame->x1[j+1]-flame->x1[j]);
      else if (j == flame->n1-1)
        Dx_Y = (Y[k+i]-Y[kM1+i])/(flame->x1[j]-flame->x1[j-1]);
      else {
        if (flame->v1[j] <= 0.0)
          Dx_Y = (Y[k+i]-Y[kM1+i])/(flame->x1[j]-flame->x1[j-1]);
        else
          Dx_Y = (Y[kP1+i]-Y[k+i])/(flame->x1[j+1]-flame->x1[j]);
      }
      
      if (!driver->ENERGY) /* temperature equation */
        flame->sumCp[j] += chem->cp[i]*Dx_Y; /* for diffusion term */

    } /* end loop over non-steady state species */
  }

  /* loop over mesh points to solve and ghost cells and update properties */
  if (chem->soot) {
    for (j=flame->j0m; j<=flame->jNp; j++) {
      k   = j*flame->nVars;      /* array index for start of solution */
      /* compute heat conductivity according to SMOOKE */
      // lambda = 2.58e-5 * pow(T/298.0,0.7) * cp 
      /* compute rhoD from Le = 1 assumption --> rhoD = lambda/cp */
#ifdef FM_VS_MDUC
      flame->D[j]  = 2.58e-5 * pow(T[j]/298.0,0.7)/flame->rho[j];
      flame->mu[j] = 0.75*flame->rho[j]*flame->D[j];
#else  
      if (chem->readTrans) {
         flame->mu[j] = calcMixViscosity(chem->nSpecies, T[j], &Y[k], chem);
      }
      else {
          fprintf(stderr, "Transport data required!!!\n");
      }
      /* Invoking Pr = 0.75 /TODO: Detailed computation of lambda/cp */
      flame->D[j]  = flame->mu[j] / (flame->rho[j]*0.75);
#endif /* FM_VS_MDUC */
      if (driver->thermophoresis) 
         flame->xi[j] = 0.55*0.25*flame->mu[j]/flame->D[j]/T[j];

    }
  }

  /* Update convection velocity for energy equation */
  flame->CFL = 1.0E-60;
  for (j=flame->j0; j<=flame->jN; j++) {
    jst = stencil_x(j, flame);        /* start of stencil */
    
    /* ----- compute convection term ----- */
    Dx_cpmix = grad_x(jst, 0, j, 0, 1, flame->cpmix, flame);
    vT = 0.0;
    if (!driver->DNS) /* ignore cpmix terms if comparing with DNS */
      vT = - 0.5*flame->chiA[j]/flame->cpmix[j]*(Dx_cpmix + flame->sumCp[j]);

    flame->vT[j] = flame->v1[j] + vT;
 
    /* ----- compute CFL number */
    if (j != 0 && j != flame->n1-1) 
      flame->CFL = fmax(flame->CFL, flame->vT[j]*(flame->x1[j]-flame->x1[j-1]));
  }

  /* Update convection velocity for soot moments */
  if (chem->soot) {
    for (j=flame->j0; j<=flame->jN; j++) {
      if (j == 0 || j == flame->n1-1) {
          flame->vMxmyp[j] = 0.0;
          flame->vMxy[j]  = 0.0;
      }
      else { 
        jst = stencil_x(j, flame); /* start of stencil */

        /* Compute required gradients */
        Dx_chi = grad_x(jst, 0, j, 0, 1 , flame->chiA, flame);
        Dx_rho = grad_x(jst, 0, j, 0, 1 , flame->rho, flame);
        Dx_D   = grad_x(jst, 0, j, 0, 1 , flame->D, flame);
        //Dx_chi = (flame->chiA[j]-flame->chiA[j-1])/(flame->x1[j]-flame->x1[j-1]);
        //Dx_rho = (flame->rho[j]-flame->rho[j-1])/(flame->x1[j]-flame->x1[j-1]);
        //Dx_D   = (flame->D[j]-flame->D[j-1])/(flame->x1[j]-flame->x1[j-1]);
        vA = - flame->rhoDotSoot[j]*flame->x1[j] 
#ifdef RHODOT
         + flame->rhoDotSpray[j] * (1.0 - flame->x1[j]);
#endif /* RHODOT */
           + 0.25 * (flame->rho[j]*Dx_chi 
          	   + 2.0 * flame->chiA[j] * Dx_rho 
          	   + flame->chiA[j] / flame->D[j]*flame->rho[j]*Dx_D);

        if (driver->thermophoresis) {
          Dx_T  = grad_x(jst, 0, j, flame->sT, flame->nVars, Y, flame);
          Dx_xi = grad_x(jst, 0, j, 0, 1 , flame->xi, flame);
          vT    = - 2.0*flame->chiA[j]*flame->xi[j]*Dx_T;
        }
        else vT = 0.0;

        /* convective velocity in mixture fraction space caused by 
           vA - Advection
           vT - Thermophoresis */
        flame->vMxy[j] = vA + vT;
        flame->vMxmyp[j] = - 0.25 * (flame->rho[j]*Dx_chi
                   + 2.0 * flame->chiA[j] * Dx_rho
                   + flame->chiA[j] / flame->D[j]*flame->rho[j]*Dx_D);; // multiplication with DN0_DZ and V0/S0 will be done later

        /* ----- compute CFL number */
        flame->CFL = fmax(flame->CFL, flame->vMxy[j]/(flame->x1[j]-flame->x1[j-1]));
        flame->CFL = fmax(flame->CFL, flame->vMxmyp[j]/(flame->x1[j]-flame->x1[j-1]));
      }
    } /* end mesh loop for species equations */
  }
  /* Diffusive CFL number */
  for (j=flame->j0; j<=flame->jN; j++) {
    jst = stencil_x(j, flame);        /* start of stencil */
    
    if (j != 0 && j != flame->n1-1) {
      flame->CFL = fmax(flame->CFL, 0.5*flame->chiA[j]/
      		pow(flame->x1[j]-flame->x1[j-1],2));

      if (chem->soot) 
#ifndef NFREEPATH
         flame->CFL = fmax(flame->CFL, 0.5*flame->chiA[j]*DN0_DZ*pow(S0,6.0/Df-2.0)*pow(V0,2.0-4.0/Df)
      		*pow(flame->x1[j]-flame->x1[j-1],2));
#else
         flame->CFL = fmax(flame->CFL, 0.5*flame->chiA[j]*DN0_DZ*pow(S0,3.0/Df-1.0)*pow(V0,1.0-2.0/Df)
      		*pow(flame->x1[j]-flame->x1[j-1],2));
#endif /* NFREEPATH */
    }
  }
}
/* -------------------------------------------------------------------------- */
void flameletQ_RHS_1D(sdialStruct *sdial, driveStruct *driver, chemStruct *chem, double *RHS)
{
  int i;
  if (driver->output_w) {
    for (i = 0; i < chem->nRxnTotal; i++) {
        RHS[i] = chem->w_tmp[i];
    }
  }
}

void flameletT_RHS_1D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                      flowStruct *flow, double *Y, double *RHS)
{
  int    i, j, k, jst, kM1, kP1, kL;
  double Dx_T, Dxx_T, Dx_cpmix; /* spatial derivatives */
  double src, diff, conv;       /* source, diffusion, and convection terms */
  double vT;
  double *cpmix = flame->cpmix;
  const double BOLTZMANN=5.67051E-8; // [W/(m^2 K^4)]

  /* ----- solve the temperature equation ----- */
  for (j=flame->j0; j<=flame->jN; j++) {
#ifdef COUNTERFLOW
    if (!(j == 0 || j == flame->n1-1)) 
#endif /* COUNTERFLOW */
    {
    kM1 = (j-1)*flame->nVars + flame->sT;
    kP1 = (j+1)*flame->nVars + flame->sT;
    k   = j*flame->nVars + flame->sT; /* array index temperature */
    kL  = (j-flame->j0)*flame->nVars + flame->sT;
    jst = stencil_x(j, flame);        /* start of stencil */
    
    /* calculate the gradients in normalised space */
    Dxx_T    = grad_xx(jst, 0, j, flame->sT, flame->nVars, Y, flame);
    
    /* ----- compute diffusion term ----- */
    diff = 0.5*flame->chiA[j]*Dxx_T;
    
//#ifdef UPWIND
    if (j == 0) {
      Dx_T     = (Y[kP1]     - Y[k]    )/(flame->x1[j+1]-flame->x1[j]);
    }
    else if (j == flame->n1-1) {
      Dx_T     = (Y[k]    - Y[kM1]    )/(flame->x1[j]-flame->x1[j-1]);
    }
    /* DEBUG: use an upwinding scheme */
    else {
      if ((flame->vT[j]) >= 0.0) {
        Dx_T     = (Y[k]     - Y[kM1]    )/(flame->x1[j]-flame->x1[j-1]);
      }
      else {
        Dx_T     = (Y[kP1]     - Y[k]    )/(flame->x1[j+1] - flame->x1[j]);
      }
    }
//#else
//    Dx_T    = grad_x(jst, 0, j, flame->sT, flame->nVars, Y, flame);
//#endif /*UPWIND */
    
    conv = Dx_T*flame->vT[j];

    /* compute pressure source term */
    flame->srcP[j] = 1.0/flame->rho[j]*flow->Dt_P;
    
    /* ----- compute chemical source term ----- */
    src  = 1.0/flame->cpmix[j]*(flame->srcP[j] - flame->srcH[j] + flame->visc_src[j]);
    
    /* DEBUG: force convective velocity to be zero at the boundary */
    if ((j == 0) && ((flame->vT[j]) > 0.0))
      conv = 0.0;
    else if  ((j == flame->n1-1) && ((flame->vT[j]) < 0.0))
      conv = 0.0;

    /* ----- copy sources over if needed ----- */
    if (driver->output_tsrc) {
      flame->diff_temp[j] = diff;
      flame->conv_temp[j] = conv;
      flame->chem_temp[j] = src;
    }
    // TODO: Put this into a seperate subroutine!
    if (driver->radiation) {
       /* Gas phase radiation (only CO2 / H2O / CH4 / CO considered) */
       // Grosshandler, W. L., RADCAL: A Narrow-Band Model for Radiation Calculations in a Combustion Environment, NIST technical note 1402, 1993.
       /* 0: H2O 1:CO2 2:CH4 3:CO */
       /* compute partial pressures of radiating species */ 
       int is=0;
       for (i = 0; i < chem->nRad; i++) {
         is = chem->iRad[i];
         chem->pRad[i]  = Y[j*flame->nVars+is]/chem->W[is]*flame->Mw[j]*flow->P;
       }
       
       /* compute emission coefficients */
       double T = fmax(Y[k],1.0);    
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

       double fv=0.0;
       /* Soot radiation */
       if (chem->soot>0) {
          if (chem->sootMode<4) {
             /* compute soot volume fraction in [1/m^3] */
             /* NOTE: for soot mode == 1 we are solving M/rho */
             if (chem->sootMode == 1) 
                   fv = fmax(0.0,flame->Y[j*flame->nVars+flame->sSoot+1]*flame->rho[j]*SootMolarMass/SootDensity);
             else
                   fv = fmax(0.0,flame->sootMom[j*chem->nSootMom+1]*SootMolarMass/SootDensity);
          } /* soot mode < 4 */
       } /* soot */
       double alphas = -3.75e5 + 1735.0 * T;
       if (alphas < 0.0) alphas = 0.0;

       rad += alphas * fv;
       flame->radiation[j] = rad * 2.0*BOLTZMANN * (T2*T2-pow(chem->TuRad,4));
       src -= flame->radiation[j]/(flame->rho[j]*flame->cpmix[j]);
    }
    /* ----- update RHS of temperature equation ----- */
#ifndef SOOTONLY
    RHS[kL] = driver->DIFF_MODE[0]*diff - driver->CONV_MODE*conv 
              + src - flame->rhoDotSoot[j]*Y[k]/flame->rho[j];
#else
    RHS[kL] = 0.0;
#endif /* SOOTONLY */
  }
#ifdef COUNTERFLOW
  else {
    kL  = (j-flame->j0)*flame->nVars + flame->sT;
    RHS[kL] = 0.0;
  }
#endif
  } /* end mesh loop for temperature equation */
}
/* -------------------------------------------------------------------------- */

void flameletChi_RHS_1D(driveStruct *driver,chemStruct *chem,flameStruct *flame,
                        flowStruct *flow, double *T, double *Y, double *RHS)
{
  int    i, j, k, jst;
  double Dx_Y, Dxx_Y, Dx_Sdot;     /* spatial derivatives */
  double wmix, nSpecies;
  double src, diff, conv, strain;  /* individual terms */

  nSpecies = chem->nSpecSolve; /* solve only non-steady state species */

  /* enforce boundary conditions */
  Y[flame->sChi] = 0.0;
  Y[(flame->n1-1)*flame->nVars + flame->sChi] = 0.0;

  /* ----- solve the scalar dissipation rate equation ----- */
  for (j=1; j<flame->n1-1; j++) {
    k   = j*flame->nVars + flame->sChi; /* array index for chi */
    jst = stencil_x(j, flame);          /* start of stencil */
    
    /* ----- compute local properties ----- */
    wmix = calcMixtureW(nSpecies, &Y[j*flame->nVars], chem->W);
    flame->rho[j] = wmix*flow->P/(flame->RGAS*T[j]);
    
    /* calculate the gradients in normalised space */
    Dx_Y     = grad_x(jst, 0, j, flame->sChi, flame->nVars, Y, flame);
    Dxx_Y    = grad_xx(jst, 0, j, flame->sChi, flame->nVars, Y, flame);
    Dx_Sdot  = grad_x(jst, 0, j, 0, 1, flow->avgSdot_M, flame);
    
    /* ----- compute diffusion term ----- */
    diff = 0.5*Y[k]*Dxx_Y;
    
    /* ----- compute convection term ----- */
    conv = 0.25*Dx_Y*Dx_Y;
    
    /* ----- strain rate term ----- */
    strain = -2*flow->strainCondZ[j]*Y[k];
    
    /* ----- compute evaporative source term ----- */
    src = 2*Y[k]/flame->rho[j]*Dx_Sdot
      - flow->avgSdot_M[j]/flame->rho[j]*Dx_Y;
    
    /* ----- update RHS of temperature equation ----- */
    RHS[k] = driver->DIFF_MODE[0]*diff - driver->CONV_MODE*conv - strain +src;
  } /* end mesh loop for temperature equation */
}
/* -------------------------------------------------------------------------- */

int getTimeIndex(ioStruct *io, double tIn)
{
  int iLow = 0, iUp = io->nFiles-1;
  int step = 0;

  if (io->chiData[iLow].tChi > tIn)
   return -1;
  else if (io->chiData[iUp].tChi < tIn) {
   return -2;
   fprintf(stderr, "Warning: solution time is greater than time of the last file!");   
  }

  while(iUp-iLow > 1 && step < 1000) {
    int iMean = ceil((iLow+iUp)*0.5);
    if(io->chiData[iMean].tChi < tIn)
      iLow = iMean; 
    else
      iUp = iMean; 

    step++;
  }
  
  if (step >= 1000){
    fprintf(stderr, "ERROR: Bisection failed");
    exit(2);
  }

  return iLow;
}
/* -------------------------------------------------------------------------- */

void flameletPrestep1D(driveStruct *driver,flameStruct *flame, chemStruct *chem, flowStruct *flow, ioStruct *io, timeStruct *time)
{      
  int i, j, k;
  double S0,V0;
  V0       = 2.0*NuclVol;
  S0       = pow(V0,2.0/3.0);
  
  flowBound *phiStar = &flow->phiStar[1];
  flowBound *Dt_phi  = flow->Dt_phi[1];

  /* set the current delta(phi) to be the one at n-1 */
  phiStar->del = flow->D_phi1;

  if (driver->CHI_MODE == 3) {
    
    /* check if grid is identical */
    int i_m1 = getTimeIndex(io, time->t_m1);
    int i_n = getTimeIndex(io, time->t_n);

    for (j=0; j < flame->n1; j++) {
      // interpolation for previous time
      if (i_m1 >= 0) {
        if (io->chiData[i_m1].chiVal != NULL)
        flow->chi1_m1[j] = linInterp(time->t_m1, io->chiData[i_m1].tChi, io->chiData[i_m1].chiVal[j],
		 io->chiData[i_m1+1].tChi, io->chiData[i_m1+1].chiVal[j]); 
        if (io->chiData[i_m1].pdfVal != NULL)
           flame->gPDF[j] = linInterp(time->t_n, io->chiData[i_n].tChi, io->chiData[i_n].pdfVal[j],
		 io->chiData[i_n+1].tChi, io->chiData[i_n+1].pdfVal[j]); 
        if (io->chiData[i_m1].hSrcVal != NULL)
           flame->visc_src[j] = linInterp(time->t_n, io->chiData[i_n].tChi, io->chiData[i_n].hSrcVal[j],
		 io->chiData[i_n+1].tChi, io->chiData[i_n+1].hSrcVal[j]); 
        if (io->chiData[i_m1].rhoDotVal != NULL)
           flame->rhoDotSpray[j] = linInterp(time->t_n, io->chiData[i_n].tChi, 
		io->chiData[i_n].rhoDotVal[j],io->chiData[i_n+1].tChi, io->chiData[i_n+1].rhoDotVal[j]); 
//        if (io->chiData[i_m1].TDotVal != NULL)
//           flame->TDotSpray[j] = linInterp(time->t_n, io->chiData[i_n].tChi, 
//		io->chiData[i_n].TDotVal[j],io->chiData[i_n+1].tChi, io->chiData[i_n+1].TDotVal[j]); 

      } else if (i_m1 == -1) {
        if (io->chiData[0].chiVal != NULL)
        flow->chi1_m1[j] = io->chiData[0].chiVal[j]; 
        if (io->chiData[0].pdfVal != NULL)
           flame->gPDF[j] = io->chiData[0].pdfVal[j]; 
        if (io->chiData[0].hSrcVal != NULL)
           flame->visc_src[j] = io->chiData[0].hSrcVal[j]; 
        if (io->chiData[0].rhoDotVal != NULL)
           flame->rhoDotSpray[j] = io->chiData[0].rhoDotVal[j]; 
//        if (io->chiData[0].TDotVal != NULL)
//           flame->TDotSpray[j] = io->chiData[0].TDotVal[j]; 
      } else if (i_m1 == -2) {
        if (io->chiData[io->nFiles-1].chiVal != NULL)
        flow->chi1_m1[j] = io->chiData[io->nFiles-1].chiVal[j]; 
        if (io->chiData[io->nFiles-1].pdfVal != NULL)
           flame->gPDF[j] = io->chiData[io->nFiles-1].pdfVal[j]; 
        if (io->chiData[io->nFiles-1].hSrcVal != NULL)
           flame->visc_src[j] = io->chiData[io->nFiles-1].hSrcVal[j]; 
        if (io->chiData[io->nFiles-1].rhoDotVal != NULL)
           flame->rhoDotSpray[j] = io->chiData[io->nFiles-1].rhoDotVal[j]; 
//        if (io->chiData[io->nFiles-1].TDotVal != NULL)
//           flame->TDotSpray[j] = io->chiData[io->nFiles-1].TDotVal[j]; 
      }

      // interpolation for end time of step
      if (i_n >= 0) {
        flow->chi1[j] = linInterp(time->t_n, io->chiData[i_n].tChi, io->chiData[i_n].chiVal[j],
		 io->chiData[i_n+1].tChi, io->chiData[i_n+1].chiVal[j]); 
      } else if (i_n == -1) {
        flow->chi1[j] = io->chiData[0].chiVal[j]; 
      } else if (i_n == -2) {
        flow->chi1[j] = io->chiData[io->nFiles-1].chiVal[j]; 
      }
      //flow->chi1[j] = flow->chi1[j];
      //flow->chi1_m1[j] = flow->chi1_m1[j];
    }
  }

  /* loop over mesh points */
  for (j=0; j<flame->n1; j++) {
    k = j*(flame->nVars);
    
    /* calculate the convective velocity for case with mixture fraction
       NOTE: pressure source term must be added for enthalpy coordinate.
       Is added in calcRHS */
    if (flow->D_phi1_m1 > 0.0) {
      flame->v1[j] = -1.0/flow->D_phi1_m1*( flow->Dt_phi1min + 
                     flame->x1[j]*(flow->Dt_phi1max-flow->Dt_phi1min) );

      /* compute the transformed scalar dissipation rate */
      if (driver->CHI_MODE != 2) {
        if (j == 0 || j == flame->n1-1) /* enforce zero at boundaries */
          flame->chiA[j] = 0.0;
        else
          flame->chiA[j]=1.0/(flow->D_phi1_m1*flow->D_phi1_m1)*flow->chi1_m1[j];
      }
    }
    else {
      flame->chiA[j] = 0.0;
      flame->v1[j] = 0.0;
    }
    
    /* ----- store the previous solution ----- */
    /* loop over species */
    for (i=0; i<flame->nVars; i++) {
      flame->Y_m2[k+i] = flame->Y_m1[k+i];
      flame->Y_m1[k+i] = flame->Y[k+i];
      
      flame->chemSrcY[k+i] = flame->chemSrcY_tmp[k+i] = 0.0;
    } /* end species loop */
         
    /* initialize heat release and source terms to zero if specified */
    flame->HRrate[j] = 0.0;
    if (driver->output_tsrc) 
      flame->chem_src[j] = flame->conv_src[j] = flame->diff_src[j] = 0.0;

  } /* end mesh loop */

  if (driver->output_w) {
     for (i = 0; i < chem->nRxnTotal; i++)
       chem->w_pre[i] = 0.0;
  } 

}
void flameletPreiteration1D(driveStruct *driver,flameStruct *flame, chemStruct *chem, flowStruct *flow,
			    timeStruct *time, mpiStruct *mpi, double t)
{      
  int i,j,k;
  double *T; /* local temperature */

  time->dt_max = 1.0E60;
//#ifdef CHIPREITER
  /* compute a coupled equation for the scalar dissipation rate */
  if (driver->CHI_MODE != 2) 
    updateFlow1D(t, driver, time, flame, flow);
//#endif /* CHIPREITER */
#ifdef PREITER
  /* get temperature based on form of energy equation */
  if (driver->ENERGY) {
    T = flame->T;
    /* update to most current temperature and enthalpy */
    for (j=0; j<flame->n1; j++) {
      k = j*flame->nVars;
      T[j] = calcTfromH(chem->nSpecies, chem->nSpecSolve,
                        flame->Y[k+flame->sT], &flame->Y[k], 1000.0, chem);
    } /* end loop over grid */
  }
  else { /* temperature is part of solution */
    T = (double*) calloc(flame->n1, sizeof(double));
    extractVarArray(flame->nVars,flame->n1,1,flame->sT,flame->Y,T);
  }
  /* Update convective velocities */
  flameletConvVel1D(driver,chem,flame,flow,flame->Y,T);
  free(T);

  /* Update maximum time-step for CFL<0.9 */
  time->dt_max_ = 0.9/flame->CFL;
#ifdef ACCESS_MPI
  MPI_Allreduce(&time->dt_max_,&time->dt_max, 1, MPI_DOUBLE_PRECISION, MPI_MIN, mpi->comm);
#endif /* ACCESS_MPI */
#endif /* PREITER */
}
/* -------------------------------------------------------------------------- */

void updateEnergy1D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                    flowStruct *flow, double dt, double *Y, double *Yp)
{
  int j, k;

  /* update the current enthalpy to H + 1/rho*dPdt - q */
  for (j=0; j<flame->n1; j++) {
    k = j*flame->nVars;
    
    /* make sure the temperature is up-to-date */
    flame->T[j] = calcTfromH(chem->nSpecies, chem->nSpecSolve,
                             Y[k+flame->sT], &Y[k], 1000, chem);
    
    /* update the enthalpy */
    Yp[k+flame->sT] = flame->H[j] + 1.0/flame->rho[j]*flow->Dt_P*dt;
  }
}
/* -------------------------------------------------------------------------- */

void updateFlow1D(double t, driveStruct *driver, timeStruct *time,
                  flameStruct *flame, flowStruct *flow)
{
  int i, j, k;
  double tmp = 1.0;
  
  flowBound *phiStar = &flow->phiStar[1];
  flowBound *phi     = flow->phi[1];
  flowBound *Dt_phi  = flow->Dt_phi[1];

  if (flow->D_phi1 != 0) 
    tmp = linInterp(t, time->t_m1, flow->D_phi1_m1, time->t_n, flow->D_phi1);

  /* compute the current delta phi */
  phiStar->del = flow->D_phi1;

  if (phiStar->del != 0) tmp = phiStar->del;

  /* loop over mesh points */
  for (j=0; j < flame->n1; j++) {
    k = j*(flame->nVars);
    
    /* should update P here as well */
    
    flame->v1[j] = flame->chiA[j] = 0.0;
    if (flow->D_phi1 != 0) 
      /* compute convective velocity due to grid normalization */
      flame->v1[j] = -1.0/tmp*( flow->Dt_phi1min +
                             flame->x1[j]*(flow->Dt_phi1max-flow->Dt_phi1min) );

    /* add pressure source to convective velocity if enthalpy is direction */
    /* compute pressure source term */
    flame->srcP[j] = 1.0/flame->rho[j]*flow->Dt_P;
    if (flow->phiStar[1].del > 0.0)
      flame->v1[j] += driver->H_DIM*flame->srcP[j]/flow->phiStar[1].del;

    /* compute the transformed scalar dissipation rate */
    if (j == 0 || j == flame->n1-1) /* enforce zero at boundaries */
      flame->chiA[j] = 0.0;
    else
#ifndef CONSTANTCHI
      flame->chiA[j] = 1.0/(tmp*tmp)*linInterp(t, time->t_m1,
                           flow->chi1_m1[j], time->t_n, flow->chi1[j]);
#else
      flame->chiA[j] = flow->chi1[j];
#endif
  } /* end mesh loop */
}
/* -------------------------------------------------------------------------- */

void flamelet1D_monitorHeader(ioStruct *io, sdialStruct *sdial)
{
  int i;
  fprintf(io->monitor, "# it<1> subit<2> time[s]<3> dt[s]<4> Tmax[K]<5> "
          "Z|Tmax[-]<6> ");
  i = 7;
  fprintf(io->monitor, "RHS_evals<%i> ", i); i++;
  fprintf(io->monitor, "Jac_evals<%i> ", i); i++;
  fprintf(io->monitor, "Jac_evals<%i> ", i); i++;
  fprintf(io->monitor, "\n");
  fflush(io->monitor);
}
/* -------------------------------------------------------------------------- */

void flamelet1D_monitor(int nst, double t, double dtPre, sdialStruct *sdial,
                        ioStruct *io, flowStruct *flow, flameStruct *flame)
{
  int i, j, ind, jMaxT;
  double Tmax, TmaxZ;

  if (sdial->totalIt == 0)
    flamelet1D_monitorHeader(io,sdial);

  /* get the maximum temperature in the domain at locataion */
  if (flame->T) {
    Tmax  = getMaxValue(flame->n1,flame->n2,0,1,flame->T);
    jMaxT = getIndexMax(flame->n1,1,flame->T);
  }
  else {
    Tmax  = getMaxValue(flame->n1,flame->n2,flame->sT,flame->nVars,flame->Y);
    jMaxT = getIndexMax(flame->n1,flame->sT,flame->Y);
  }
  TmaxZ = flame->x1[jMaxT];

  // print the progress to the screen along with any monitor species
  fprintf(io->monitor, "%6i %6li %15.6e %15.6e %15.6e %15.6e\n", nst,
          sdial->totalIt, t, dtPre, Tmax, TmaxZ);

  if (io->nMonitorY > 0) {
    /* print table header */
    fprintf(io->monitor, "%5s %8s %8s", "Z", "T", "chiA");
    if (flame->nDim == 2)
      fprintf(io->monitor, "%8s %8s","chiB","chiC");
    
    for (i=0; i< io->nMonitorY; ++i)
      fprintf(io->monitor, " %13s", io->monitorSpecies[i]);
    fprintf(io->monitor, "\n");
    
    /* print table */
    for (j=0; j < io->nMonitorZ; ++j) {
      ind = floor((flame->n1-1)*j/(io->nMonitorZ-1)+0.49);
      if (flame->T)
        fprintf(io->monitor, "%5.3f %8.1f %8.3g", flame->x1[ind],
                flame->T[ind], flame->chiA[ind]);
      else
        fprintf(io->monitor, "%5.3f %8.1f %8.3g", flame->x1[ind],
                sdial->Ydata[ind*flame->nVars+flame->sT], flame->chiA[ind]);
      
      if (flame->nDim == 2)
        fprintf(io->monitor, "%8.3g %8.3g", flame->chiB[ind],
                flame->chiC[ind]);
      
      for (i=0; i < io->nMonitorY; ++i)
        fprintf(io->monitor, " %13.3e",
                sdial->Ydata[ind*flame->nVars + io->monitorIndex[i]]);
      
      fprintf(io->monitor, "\n");
    }
    fprintf(io->monitor, "\n");
  } /* end table */
  fflush(io->monitor);
}


void flameletUpdateSoot1D(chemStruct *chem,flameStruct *flame,driveStruct *driver,ioStruct *io,
		          flowStruct *flow,int inert)
{
  int i,j,k,iS;
  double sgCoeff, oxCoeff,o2Coeff, dimProdRate;
  if (io->outputMu || chem->soot > 0) {
    for (j=flame->j0; j<=flame->jN; j++) {
      k = j*flame->nVars;
        if (chem->readTrans)
           flame->mu[j] = calcMixViscosity(chem->nSpecies, flame->Y[k+flame->sT],
                                      &flame->Y[k], chem);
    }

    if (chem->sootMode == 1) { /* copy solution into arrays for output */
      for (j=flame->j0; j<=flame->jN; j++)
        for (i=0; i<flame->nSootMom; i++) {
          flame->sootMom[j*flame->nSootMom+i] = 
            flame->Y[j*flame->nVars+flame->sSoot+i];
        }
    }
    else if (chem->sootMode == 2) {
      /* update the soot moments based on the source terms */
      /* I think this should be done in sundials loop */
      //for (j=flame->j0; j<=flame->jN; j++)
      //  for (i=0; i<flame->nSootMom; i++)
      //    flame->sootMom[j*flame->nSootMom+i] += 
      //      flame->sootMomSrc[j*flame->nSootMom+i]*dt;
    }
    else if (chem->sootMode >= 3) {
        /* compute hmom turbulence-soot interaction coeffcients */ 
        for (j=flame->j0; j<=flame->jN; j++) {
          k = j*flame->nVars;
	  /* compute rho (not up to date in case of inert flamelet */
          if (inert == 0) {
             flame->rho[j] = calcMixtureW(chem->nSpecSolve,&flame->Y[j*flame->nVars], chem->W)*
				flow->P/(flame->RGAS*flame->Y[k+flame->sT]);
          } /*inert*/
          /* compute and store concentrations of all species at local point */
          for (iS=0; iS<chem->nSootSp; iS++)
            chem->sootRxnConc[iS] = flame->rho[j]*flame->Y[k+chem->iSootSp[iS]]/chem->W[chem->iSootSp[iS]];

          sootHMOM_GasPhase_Coeff(flame->Y[k+flame->sT], flame->mu[j],
          	                flame->rho[j], flame->Mw[j], chem->sootRxnConc, 
          			&sgCoeff, &oxCoeff,&o2Coeff, &dimProdRate);

	  flame->sgCoeff[j] 	    = sgCoeff;
	  flame->oxCoeff[j]  	    = oxCoeff;
	  flame->o2Coeff[j] 	    = o2Coeff;
	  flame->dimProdRate[j]     = dimProdRate;
          flame->sqrtT[j] 	    = sqrt(flame->Y[k+flame->sT]);
          flame->T_mu[j]  	    = flame->Y[k+flame->sT]/flame->mu[j];
          flame->muSqrtMw_rhoSqrtT[j] = flame->mu[j]*sqrt(flame->Mw[j]/
          			(1000.0*flame->Y[k+flame->sT]*flame->rho[j]*flame->rho[j]));
        } /* grid loop */
      } /* soot mode 3 */
    } /* soot */
}

void flameletAllgatherv(chemStruct *chem,flameStruct *flame,driveStruct *driver,ioStruct *io,mpiStruct *mpi)
{
     int i,j,k;
     double *sendbuf;
     int myid = mpi->myid; 
     /*TODO: Ideally we would communicate everything in a single buffer!*/
     sendbuf = (double*) calloc(mpi->jCount[myid], sizeof(double));
     for (i=0; i<mpi->jCount[myid]; i++) sendbuf[i] = flame->Y[mpi->jStart[myid]+i];

     MPI_Allgatherv(sendbuf,mpi->jCount[myid],
                 MPI_DOUBLE_PRECISION,flame->Y,mpi->jCount,mpi->jStart,
                 MPI_DOUBLE_PRECISION,mpi->comm);

     /*TODO: is flame->Yp already upto date? */
     if (driver->SDIAL == 1) {
       for (i=0; i<mpi->jCount[myid]; i++) sendbuf[i] = flame->Yp[mpi->jStart[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->jCount[myid],
                   MPI_DOUBLE_PRECISION,flame->Yp,mpi->jCount,mpi->jStart,
                   MPI_DOUBLE_PRECISION,mpi->comm);
     }

     /* retrieve source terms for heat release computation */
     if (io->outputHR) {
       for (i=0; i<mpi->jCount[myid]; i++) sendbuf[i] = flame->chemSrcY[mpi->jStart[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->jCount[myid],
                   MPI_DOUBLE_PRECISION,flame->chemSrcY,mpi->jCount,mpi->jStart,
                   MPI_DOUBLE_PRECISION,mpi->comm);
     }
     free(sendbuf);

     if (driver->output_w) {
       sendbuf = (double*) calloc(chem->nRxnTotal, sizeof(double));
       for (i=0; i<chem->nRxnTotal; i++) sendbuf[i] = chem->w_pre[i];
       MPI_Allreduce(sendbuf,chem->w_pre,chem->nRxnTotal,MPI_DOUBLE_PRECISION,MPI_SUM,mpi->comm);
     }

     /* retrieve properties */
     sendbuf = (double*) calloc(mpi->cCount[myid], sizeof(double));

     for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->cpmix[mpi->cStart[myid]+i];
     MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                 MPI_DOUBLE_PRECISION,flame->cpmix,mpi->cCount,mpi->cStart,
                 MPI_DOUBLE_PRECISION,mpi->comm);

     for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->rho[mpi->cStart[myid]+i];
     MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                 MPI_DOUBLE_PRECISION,flame->rho,mpi->cCount,mpi->cStart,
                 MPI_DOUBLE_PRECISION,mpi->comm);

     for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->Mw[mpi->cStart[myid]+i];
     MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                 MPI_DOUBLE_PRECISION,flame->Mw,mpi->cCount,mpi->cStart,
                 MPI_DOUBLE_PRECISION,mpi->comm);

     for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->gamma[mpi->cStart[myid]+i];
     MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                 MPI_DOUBLE_PRECISION,flame->gamma,mpi->cCount,mpi->cStart,
                 MPI_DOUBLE_PRECISION,mpi->comm);

     if (driver->ENERGY) {
       for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->H[mpi->cStart[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                   MPI_DOUBLE_PRECISION,flame->H,mpi->cCount,mpi->cStart,
                   MPI_DOUBLE_PRECISION,mpi->comm);
     }
     if (chem->soot > 0 || io->outputMu) {
       for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->mu[mpi->cStart[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                   MPI_DOUBLE_PRECISION,flame->mu,mpi->cCount,mpi->cStart,
                   MPI_DOUBLE_PRECISION,mpi->comm);
     }
     free(sendbuf);
     /* soot */

     if (chem->soot > 0 && chem->sootMode <= 3) {
       sendbuf = (double*) calloc(mpi->sootN[myid], sizeof(double));

       for (i=0; i<mpi->sootN[myid]; i++) sendbuf[i] = flame->sootMom[mpi->soot0[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->sootN[myid],
                   MPI_DOUBLE_PRECISION,flame->sootMom,mpi->sootN,mpi->soot0,
                   MPI_DOUBLE_PRECISION,mpi->comm);

       for (i=0; i<mpi->sootN[myid]; i++) sendbuf[i] = flame->sootMomSrc[mpi->soot0[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->sootN[myid],
                   MPI_DOUBLE_PRECISION,flame->sootMomSrc,mpi->sootN,mpi->soot0,
                   MPI_DOUBLE_PRECISION,mpi->comm);

       for (i=0; i<mpi->sootN[myid]; i++) sendbuf[i] = flame->sootMomTrans[mpi->soot0[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->sootN[myid],
                   MPI_DOUBLE_PRECISION,flame->sootMomTrans,mpi->sootN,mpi->soot0,
                   MPI_DOUBLE_PRECISION,mpi->comm);

       for (i=0; i<mpi->sootN[myid]; i++) sendbuf[i] = flame->sootMomTransDiff[mpi->soot0[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->sootN[myid],
                   MPI_DOUBLE_PRECISION,flame->sootMomTransDiff,mpi->sootN,mpi->soot0,
                   MPI_DOUBLE_PRECISION,mpi->comm);

       for (i=0; i<mpi->sootN[myid]; i++) sendbuf[i] = flame->sootMomTransConv[mpi->soot0[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->sootN[myid],
                   MPI_DOUBLE_PRECISION,flame->sootMomTransConv,mpi->sootN,mpi->soot0,
                   MPI_DOUBLE_PRECISION,mpi->comm);

       free(sendbuf);
       if (io->outputSootMxmyp) {
          sendbuf = (double*) calloc(mpi->sootNxmyp[myid], sizeof(double));
          for (i=0; i<mpi->sootNxmyp[myid]; i++) sendbuf[i] = flame->sootMomxmyp[mpi->soot0xmyp[myid]+i];
          MPI_Allgatherv(sendbuf,mpi->sootNxmyp[myid],
                      MPI_DOUBLE_PRECISION,flame->sootMomxmyp,mpi->sootNxmyp,mpi->soot0xmyp,
                      MPI_DOUBLE_PRECISION,mpi->comm);

          free(sendbuf);
       }



       sendbuf = (double*) calloc(mpi->sootSrcPN[myid], sizeof(double));

       for (i=0; i<mpi->sootSrcPN[myid]; i++) sendbuf[i] = flame->sootMomSrcPart[mpi->sootSrcP0[myid]+i];
       MPI_Allgatherv(sendbuf, mpi->sootSrcPN[myid],
	           MPI_DOUBLE_PRECISION,flame->sootMomSrcPart,mpi->sootSrcPN,mpi->sootSrcP0,
		   MPI_DOUBLE_PRECISION,mpi->comm);
       free(sendbuf);

     }
     if (chem->soot > 0 && chem->sootMode >= 3) {
       sendbuf = (double*) calloc(mpi->cCount[myid], sizeof(double));
       for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->sqrtT[mpi->cStart[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                   MPI_DOUBLE_PRECISION,flame->sqrtT,mpi->cCount,mpi->cStart,
                   MPI_DOUBLE_PRECISION,mpi->comm);
       for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->T_mu[mpi->cStart[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                   MPI_DOUBLE_PRECISION,flame->T_mu,mpi->cCount,mpi->cStart,
                   MPI_DOUBLE_PRECISION,mpi->comm);
       for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->muSqrtMw_rhoSqrtT[mpi->cStart[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                   MPI_DOUBLE_PRECISION,flame->muSqrtMw_rhoSqrtT,mpi->cCount,mpi->cStart,
                   MPI_DOUBLE_PRECISION,mpi->comm);
       for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->sgCoeff[mpi->cStart[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                   MPI_DOUBLE_PRECISION,flame->sgCoeff,mpi->cCount,mpi->cStart,
                   MPI_DOUBLE_PRECISION,mpi->comm);
       for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->oxCoeff[mpi->cStart[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                   MPI_DOUBLE_PRECISION,flame->oxCoeff,mpi->cCount,mpi->cStart,
                   MPI_DOUBLE_PRECISION,mpi->comm);
       for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->o2Coeff[mpi->cStart[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                   MPI_DOUBLE_PRECISION,flame->o2Coeff,mpi->cCount,mpi->cStart,
                   MPI_DOUBLE_PRECISION,mpi->comm);
       for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->dimProdRate[mpi->cStart[myid]+i];
       MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                   MPI_DOUBLE_PRECISION,flame->dimProdRate,mpi->cCount,mpi->cStart,
                   MPI_DOUBLE_PRECISION,mpi->comm);
       free(sendbuf);
     }
     if (driver->radiation) {
        sendbuf = (double*) calloc(mpi->cCount[myid], sizeof(double));
        for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->radiation[mpi->cStart[myid]+i];
        MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                    MPI_DOUBLE_PRECISION,flame->radiation,mpi->cCount,mpi->cStart,
                    MPI_DOUBLE_PRECISION,mpi->comm);
        free(sendbuf);
     }
     if (io->outputRhoSigma) {
        sendbuf = (double*) calloc(mpi->cCount[myid], sizeof(double));
        for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->rhoSigma[mpi->cStart[myid]+i];
        MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                    MPI_DOUBLE_PRECISION,flame->rhoSigma,mpi->cCount,mpi->cStart,
                    MPI_DOUBLE_PRECISION,mpi->comm);
        free(sendbuf);
     }
     if (io->outputConvVel) {
        sendbuf = (double*) calloc(mpi->cCount[myid], sizeof(double));
        if (chem->soot>0) {
           for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->vMxmyp[mpi->cStart[myid]+i];
           MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                       MPI_DOUBLE_PRECISION,flame->vMxmyp,mpi->cCount,mpi->cStart,
                       MPI_DOUBLE_PRECISION,mpi->comm);

           for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->vMxy[mpi->cStart[myid]+i];
           MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                       MPI_DOUBLE_PRECISION,flame->vMxy,mpi->cCount,mpi->cStart,
                       MPI_DOUBLE_PRECISION,mpi->comm);
        }
        for (i=0; i<mpi->cCount[myid]; i++) sendbuf[i] = flame->vT[mpi->cStart[myid]+i];
        MPI_Allgatherv(sendbuf,mpi->cCount[myid],
                    MPI_DOUBLE_PRECISION,flame->vT,mpi->cCount,mpi->cStart,
                    MPI_DOUBLE_PRECISION,mpi->comm);
        free(sendbuf);
     }     
}
/* -------------------------------------------------------------------------- */
