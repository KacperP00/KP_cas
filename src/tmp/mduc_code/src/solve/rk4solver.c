#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rk4solver.h"

int rk4init(void *mem) {
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
#ifdef ACCESS_MPI
  mpiStruct    *mpi    = mduc->mpi;
#endif /* ACCESS_MPI */

  /* allocate solution and rhs vector */
  rk->Y   = (double*) calloc(sdial->nEqL, sizeof(double));
  rk->Ytmp= (double*) calloc(sdial->nEqL, sizeof(double));
  rk->k1  = (double*) calloc(sdial->nEqL, sizeof(double));
  rk->k2  = (double*) calloc(sdial->nEqL, sizeof(double));
  rk->k3  = (double*) calloc(sdial->nEqL, sizeof(double));
  rk->k4  = (double*) calloc(sdial->nEqL, sizeof(double));
}

int rk4solve(void *mem) {
  int i,j,kl,kg, ierr, i_,j_;
  double dt, dt_,dY_max,RHS_abs;
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
#ifdef ACCESS_MPI
  mpiStruct    *mpi    = mduc->mpi;
#endif /* ACCESS_MPI */

  /* copy initial values to solution vector */
  for (j=flame->j0;j<=flame->jN;j++) {
     for (i=0;i<flame->nVars;i++) {
        kg    = j * flame->nVars + i;
        kl = (j-flame->j0) * flame->nVars + i;
        rk->Y[kl] =     flame->Y[kg];
     }
  }
  /* communicate */
  bCommFlameletRK(rk->Y,sdial); 

  /* advance  */
  /* compute the current flamelet RHS */
  ierr = flameletRHS_1D(time->t, flame->Y, rk->k1, sdial);

  /* get time step size */
  dt_ = fmin(rk->dt_max, rk->t_end-time->t);
 
  for (j=flame->j0;j<=flame->jN;j++) {
     for (i=0;i<flame->nVars;i++) {
        kl = (j-flame->j0) * flame->nVars + i;
        //dY_max = (fabs(rk->Y[kl])*driver->RELTOL+driver->ABSTOL); 
        dY_max = (fabs(rk->Y[kl])*driver->RELTOL+driver->ABSTOL); 
        RHS_abs = fmax(fabs(rk->k1[kl]),driver->ABSTOL);
        if (dt_ > dY_max/RHS_abs) {
           dt_ = dY_max/RHS_abs;
	   j_  = j;
	   i_  = i;
        }
     }
  }
#ifdef ACCESS_MPI
  MPI_Allreduce(&dt_,&dt, 1, MPI_DOUBLE_PRECISION, MPI_MIN, mpi->comm);
  /* get time step limiting variable */
  j = -1;
  i = -1;
  if (dt == dt_) {
     i = i_;
     j = j_;
  }
  MPI_Allreduce(&i,&rk->iCrit, 1, MPI_INT, MPI_MAX, mpi->comm);
  MPI_Allreduce(&j,&rk->jCrit, 1, MPI_INT, MPI_MAX, mpi->comm);
#else
  dt = dt_;
#endif /*ACCESS_MPI*/
  
     
  //dt = fmax(1.0E-11,dt);
  /* compute k2 */
  for (j=flame->j0;j<=flame->jN;j++) {
     for (i=0;i<flame->nVars;i++) {
        kl = (j-flame->j0) * flame->nVars + i;
//        rk->Ytmp[kl] = fmax(rk->Y[kl] + rk->k1[kl] * dt * 0.5,0.0);
        rk->Ytmp[kl] = rk->Y[kl] + rk->k1[kl] * dt * 0.5;
     }
  }
  bCommFlameletRK(rk->Ytmp,sdial); 

  ierr = flameletRHS_1D(time->t+dt*0.5, flame->Y, rk->k2, sdial);

  /* compute k3 */
  for (j=flame->j0;j<=flame->jN;j++) {
     for (i=0;i<flame->nVars;i++) {
        kl = (j-flame->j0) * flame->nVars + i;
//        rk->Ytmp[kl] = fmax(rk->Y[kl] + rk->k2[kl] * dt * 0.5,0.0);
        rk->Ytmp[kl] = rk->Y[kl] + rk->k2[kl] * dt * 0.5;
     }
  }
  bCommFlameletRK(rk->Ytmp,sdial); 

  ierr = flameletRHS_1D(time->t+dt*0.5, flame->Y, rk->k3, sdial);

  /* compute k4 */
  for (j=flame->j0;j<=flame->jN;j++) {
     for (i=0;i<flame->nVars;i++) {
        kl = (j-flame->j0) * flame->nVars + i;
//        rk->Ytmp[kl] = fmax(rk->Y[kl] + rk->k3[kl] * dt,0.0);
        rk->Ytmp[kl] = rk->Y[kl] + rk->k3[kl] * dt;
     }
  }
  bCommFlameletRK(rk->Ytmp,sdial); 

  ierr = flameletRHS_1D(time->t+dt, flame->Y, rk->k4, sdial);

  /* copy initial values to solution vector */
  for (j=flame->j0;j<=flame->jN;j++) {
     for (i=0;i<flame->nVars;i++) {
        kl = (j-flame->j0) * flame->nVars + i;
//        rk->Y[kl] = fmax( rk->Y[kl] + (1.0/6.0) * 
//            (rk->k1[kl] + 2.0*rk->k2[kl] + 2.0*rk->k3[kl] + rk->k4[kl]) * dt,0.0);
        rk->Y[kl] =  rk->Y[kl] + (1.0/6.0) * 
            (rk->k1[kl] + 2.0*rk->k2[kl] + 2.0*rk->k3[kl] + rk->k4[kl]) * dt;
     }
  }

  /* communicate */
  bCommFlameletRK(rk->Y,sdial); 
 
  time->t += dt;
  rk->step += 1;
  rk->dt = dt;

  if (time->t >= rk->t_end) return 0;
  else 		             return 1;
}

// TODO: replace sdial struct by rk struct
int bCommFlameletRK(double *Y, void *data)
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

int rk4getStep(void *mem)
{
  mducStruct   *mduc   = (mducStruct*) mem;
  rkStruct     *rk     = mduc->rk;
  return rk->step;
}
double rk4getStepSize(void *mem)
{
  mducStruct   *mduc   = (mducStruct*) mem;
  rkStruct     *rk     = mduc->rk;
  return rk->dt;
}

int rk4setStopTime(void *mem, double tEnd)
{
  mducStruct   *mduc   = (mducStruct*) mem;
  rkStruct     *rk     = mduc->rk;
  rk->t_end = tEnd;
}

int rk4setMaxStep(void *mem, double dt_max)
{
  mducStruct   *mduc   = (mducStruct*) mem;
  rkStruct     *rk     = mduc->rk;
  rk->dt_max = dt_max;
}
