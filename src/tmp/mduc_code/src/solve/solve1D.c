#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "solve1D.h"
#include "flamelet1D.h"
#include "metric.h"
#include "solveUtils.h"
#include "utils.h"
#include "engine.h"
#include "engineUtils.h"
#include "fileIO.h"
#include "display.h"
#include "accessChem.h"
#include "chemMech.h"
#include "transport.h"
#include "sootHMOM.h"

int cvodeRHS_0D(realtype t, N_Vector u, N_Vector udot, void *data)
{
  int i,iCAD;
  realtype *Y, *RHS;
  double   P, T, m, dmdt, vol_u, dVdt, rho, qdot;
  double   CAD, CAD_rel;
  sdialStruct *cvode = (sdialStruct*) data;

  // dereference memory assigned for structs
  driveStruct  *driver = (driveStruct*)  cvode->args[0];
  chemStruct   *chem   = (chemStruct*)   cvode->args[1];
  flowStruct   *flow   = (flowStruct*)   cvode->args[2];
  flameStruct  *flame  = (flameStruct*)  cvode->args[3];
  timeStruct   *time   = (timeStruct*)   cvode->args[4];
  engineStruct *engine = (engineStruct*) cvode->args[5];

  // assign pointers to vector memory for solution vector and right-hand-side   
  Y   = NV_DATA_S(u);
  RHS = NV_DATA_S(udot);

  /* get the temperature */
  T = Y[flame->sT];

  if (driver->ENGINE) {
    /* ----- pre-compute some engine quantities ----- */
    CAD     = time_to_cad(t,engine);
    CAD_rel = fmod(CAD, 720);

    if (driver->PRES_MODE == 3)
      P = Y[flame->sP];
    else if (driver->PRES_MODE == 4)
      P = flow->P;

    /* volume and rate */
    vol_u   = engine->vol = calcCylVol(time_to_cad(t,engine),engine);
    iCAD    = gridLocate(engine->nCAD,engine->CAD_vec,CAD_rel);
    dVdt    = linInterp(CAD_rel, engine->CAD_vec[iCAD], engine->dVdt[iCAD],
                        engine->CAD_vec[iCAD+1], engine->dVdt[iCAD+1]);

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

  /* compute the r.h.s. due to reaction */
  calcY_RHS_0D(driver,chem,flame,flow,cvode->row,cvode->col,T,Y,RHS);
  calcT_RHS_0D(driver,chem,flame,flow,cvode->row,cvode->col,Y,RHS);

  /* fprintf(stderr, "mass2: %g\n", Y[flame->sM+2]); */
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
      flame->srcQ[0] = calcWoschniQdot(t,P,T,flame->gamma[0],engine); /* W */

      /* convert to correct units [W/(m-K-kg)] and add to source */
      qdot = flame->srcQ[0]*flame->R[0]*T/(P*engine->vol);
      RHS[flame->sT] += qdot/flame->cpmix[0]; /* K/s */
    }
      
    /* ----- compute the pressure rate ----- */
    if (engine->spark) {
        calcEngineSpark(t, driver, time, chem, flow, flame, engine, Y, RHS);
/*         vol_u -= engine->Vprem; */
/*         dVdt  -= engine->dVbdt; */
    }

    if (driver->PRES_MODE == 3) {
      double dMw_dt = 0.0;
      for (i=0; i<chem->nSpecies; i++)
        dMw_dt += RHS[i]/chem->W[i];
      
      flow->meanRho = rho = m/vol_u;
      RHS[flame->sP] =  flame->R[0]*T/vol_u*(dmdt - rho*dVdt)
                      + rho*(flame->R[0]*RHS[flame->sT] + T*flame->RGAS*dMw_dt);
    }
  }

  return 0;
}
/* -------------------------------------------------------------------------- */

int idaRHS_0D(realtype t, N_Vector y, N_Vector yp, N_Vector r, void *data)
{
  int i, iCAD;
  realtype *Y, *Yprime, *res, *RHS;
  double    dt, H, T, rho, dVdt;
  sdialStruct *ida = (sdialStruct*) data;

  /* dereference pointers to structs */
  driveStruct  *driver = (driveStruct*)  ida->args[0];
  chemStruct   *chem   = (chemStruct*)   ida->args[1];
  flowStruct   *flow   = (flowStruct*)   ida->args[2];
  flameStruct  *flame  = (flameStruct*)  ida->args[3];
  timeStruct   *time   = (timeStruct*)   ida->args[4];
  engineStruct *engine = (engineStruct*) ida->args[5];

  /* assign pointers to beginning of CVODE vector memory */
  Y      = NV_DATA_S(y);  /* current value of dependent variable vector y(t) */
  Yprime = NV_DATA_S(yp); /* current value of y'(t) */
  res    = NV_DATA_S(r);  /* output residual vector */

  /* assign memory for the right-hand-side */
  RHS = (double*) calloc(ida->nEq, sizeof(double));

  /* get current location from start of time step */
  dt = t - time->t_m1;
  
  /* ----- update quantities before computing derivatives ----- */
  if (driver->ENERGY) {
    /* update to most current temperature and enthalpy */
    updateEnergy1D(driver, chem, flame, flow, dt, Y, Yprime);
    T = flame->T[0];
  }
  else { /* temperature is part of solution */
    T = Y[flame->sT];
  }

  /* compute the right-hand-side of species */
  calcY_RHS_0D(driver, chem, flame, flow, ida->row, ida->col, T, Y, RHS);

  /* ----- compute energy equation ----- */
  if (driver->ENERGY) {
    /* update RHS energy to enthalpy that was computed in calcY_RHS */
    RHS[flame->sT] = flame->srcH[0];
  }
  else
    calcT_RHS_0D(driver, chem, flame, flow, ida->row, ida->col, Y, RHS);

  if (driver->ENGINE) {
    fprintf(stderr, "Warning: IDA not tested for engine calculations\n");
    /* get the current engine properties */
    double CAD     = time_to_cad(t,engine);
    double CAD_rel = fmod(CAD, 720);
    double P       = Y[flame->sP];
    engine->vol = calcCylVol(time_to_cad(t,engine),engine);

    /* check mass flow */
    calcEngineRHS_0D(t, time, chem, flow, flame, engine, Y, RHS);

    double m = 0.0, dmdt = 0.0;
    for (i=0, m=0.0; i<3; i++) {
      m    += Y[flame->sM + i];
      dmdt += RHS[flame->sM + i];
    }

    /* calculate the heat losses */
    double qdot;
    qdot = calcWoschniQdot(t,P,T,flame->gamma[0],engine);
    
    /* convert to correct units [W/(m-K-kg)] and add to source */
    qdot *= flame->R[0]*T/(Y[flame->sP]*engine->vol);
    RHS[flame->sT] += qdot/flame->cpmix[0];

    /* ----- compute the pressure rate ----- */
    iCAD = gridLocate(engine->nCAD,engine->CAD_vec,CAD_rel);
    dVdt = linInterp(CAD_rel, engine->CAD_vec[iCAD], engine->dVdt[iCAD],
                     engine->CAD_vec[iCAD+1], engine->dVdt[iCAD+1]);

    double dMw_dt = 0.0;
    for (i=0; i<chem->nSpecies; i++)
      dMw_dt += RHS[i]/chem->W[i];

    flow->meanRho = rho = m/engine->vol;
    RHS[flame->sP] =  flame->R[0]*T/engine->vol*(dmdt - rho*dVdt)
                    + rho*(flame->R[0]*RHS[flame->sT] + T*flame->RGAS*dMw_dt);
  }

  /* compute the residual */
  for (i=0; i<ida->nEq; i++)
    res[i] = Yprime[i] - RHS[i];

  free(RHS);
  return 0;
}
/* -------------------------------------------------------------------------- */

void calcY_RHS_0D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                  flowStruct *flow, int iRow, int iCol, double T, double *Y,
                  double *RHS)
{
  int    i, j, r;
  int    nSpecies;
  double P,dPdt;
  double src;         /* source term */

  if (flame->sP >= 0) {
    P    = Y[flame->sP];
    dPdt = RHS[flame->sP];
  }
  else {
    P    = flow->P;
    dPdt = flow->Dt_P;
  }

  nSpecies = chem->nSpecSolve; /* non-steady state species */
  r  = iRow*flame->n1;         /* define row location indices */
  j  = iCol;

  /* compute local density and pressure source term */
  flame->Mw[r+j]  = calcMixtureW(nSpecies, Y, chem->W);
  flame->rho[r+j] = flame->Mw[r+j]*P/(flame->RGAS*T);
  flame->srcP[j]  = 1.0/flame->rho[r+j]*dPdt;
    
  /* get enthalpy and specific heat from chemical mechanism */
  chemComputeThermoData(chem->h, chem->cp, T, chem);
  /* get specific heat of the mixture */
  flame->cpmix[j] = calcMixtureCp(nSpecies, chem->cp, Y);
  flame->R[j]     = flame->RGAS/flame->Mw[r+j];
  flame->gamma[j] = flame->cpmix[j]/(flame->cpmix[j] - flame->R[j]);
      
  /* compute and store concentrations of all species at local point */
  for (i=0; i<nSpecies; i++)
    chem->c[i] = flame->rho[r+j]*Y[i]/chem->W[i];
      
  /* get rates from mechanism */
  chemComputeProductionRates(chem->cdot, chem->w, chem->K, chem->c, chem->M,
                             T, P, chem);
      
  flame->srcH[j] = 0.0;
  for (i=0; i<nSpecies; i++) {
    /* ----- compute chemical source term ----- */
    src = chem->cdot[i]*chem->W[i]/flame->rho[r+j];
    
    if (driver->ENERGY) /* solve enthalpy */
      flame->srcH[j] += chem->h[i]*Y[i];
    else { /* solve temperature */
      flame->srcH[j] += chem->h[i]*src;   /* for chemical source term */
    }
    
    /* ----- update RHS of the species equations ----- */
    RHS[i] = src;
    
    /* compute the local heat release rate */
    flame->chemSrcY_tmp[(r+j)*flame->nVars+i] = src;
  } /* end loop over species */
}
/* -------------------------------------------------------------------------- */

void calcT_RHS_0D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                  flowStruct *flow, int iRow, int iCol, double *Y, double *RHS)
{
  int j = iCol;
  double src;

  /* ----- solve the temperature equation ----- */
  /* compute chemical source term */
  src  = 1.0/flame->cpmix[j]*(flame->srcP[j] - flame->srcH[j] + flame->visc_src[j]);
  
  if (driver->output_tsrc) {
    flame->diff_temp[j] = 0.0;
    flame->conv_temp[j] = 0.0;
    flame->chem_temp[j] = src;
  }
  
  /* ----- update RHS of temperature equation ----- */
  RHS[flame->sT] = src;  
}
/* -------------------------------------------------------------------------- */
