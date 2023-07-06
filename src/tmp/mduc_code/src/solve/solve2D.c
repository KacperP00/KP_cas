#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "solve2D.h"
#include "accessChem.h"
#include "chemMech.h"
#include "metric.h"
#include "solveUtils.h"
#include "utils.h"
#include "fileIO.h"
#include "display.h"

int cvodeRHS_2D(realtype t, N_Vector u, N_Vector udot, void *data)
{
  int j;
  double *temp;
  realtype *Y, *RHS;
  sdialStruct *cvode = (sdialStruct*) data;

  // dereference memory assigned for structs
  driveStruct *driver = (driveStruct*) cvode->args[0];
  chemStruct  *chem   = (chemStruct*)  cvode->args[1];
  flowStruct  *flow   = (flowStruct*)  cvode->args[2];
  flameStruct *flame  = (flameStruct*) cvode->args[3];
  timeStruct  *time   = (timeStruct*)  cvode->args[4];

  // assign pointers to vector memory for solution vector and right-hand-side   
  Y   = NV_DATA_S(u);
  RHS = NV_DATA_S(udot);

  // update the flow variables for current time
  updateFlow2D(t, cvode->row, driver, time, flame, flow);

  // set the current temperature from the solution array
  temp = (double*) calloc(flame->n1, sizeof(double));
  extractVarArray(flame->nVars, flame->n1, 1, flame->sT, Y, temp);

  // compute the r.h.s. for the species
  calcY_RHS_imp(driver, chem, flame, flow, cvode->row, temp, Y, RHS);

  // compute the r.h.s for the temperature
  calcT_RHS_imp(driver, chem, flame, flow, cvode->row, Y, RHS);
  
  free(temp);
  return 0;
}
/* -------------------------------------------------------------------------- */

int idaRHS_2D(realtype t, N_Vector y, N_Vector yp, N_Vector r, void *data)
{
  int i, k, j;
  realtype *Y, *Yprime, *res, *RHS, *temp, dt;
  sdialStruct *ida = (sdialStruct*) data;

  // dereference pointers to structs
  driveStruct *driver = (driveStruct*) ida->args[0];
  chemStruct  *chem   = (chemStruct*)  ida->args[1];
  flowStruct  *flow   = (flowStruct*)  ida->args[2];
  flameStruct *flame  = (flameStruct*) ida->args[3];
  timeStruct  *time   = (timeStruct*)  ida->args[4];

  // assign pointers to beginning of CVODE vector memory
  Y      = NV_DATA_S(y);  /* current value of dependent variable vector y(t) */
  Yprime = NV_DATA_S(yp); /* current value of y'(t) */
  res    = NV_DATA_S(r);  /* output residual vector */

  /* assign memory for the right-hand-side */
  RHS = (double*) calloc(ida->nEq, sizeof(double));

  // get current location from start of time step
  dt = t - time->t_m1;

  /* ----- interpolate flow variables (chi) for current time ----- */
  updateFlow2D(t, ida->row, driver, time, flame, flow);

  /* ----- update quantities before computing derivatives  ----- */
  if (driver->ENERGY)
    {
      // update to most current temperature and enthalpy
      updateEnergy2D(driver, chem, flame, flow, dt, ida->row, Y, Yprime);
      temp = &flame->T[ida->row*flame->n1];
    }
  else
    { // temperature is part of solution
      temp = (double*) calloc(flame->n1*flame->n2, sizeof(double));
      extractVarArray(flame->nVars, flame->n1, 1, flame->sT, Y, temp);
    }

  /* ----- compute species derivatives ----- */
  calcY_RHS_imp(driver, chem, flame, flow, ida->row, temp, Y, RHS);

  /* ----- compute energy equation ----- */
  if (driver->ENERGY)
    {
      // update RHS energy to enthalpy that was computed in calcY_RHS
      for (j=0; j<flame->n1; j++)
        RHS[(ida->row*flame->n1+j)*flame->nVars+flame->sT] = 
          flame->srcH[ida->row*flame->n1+j];
    }
  else
    calcT_RHS_imp(driver, chem, flame, flow, ida->row, Y, RHS);

  /* compute the residual */
  for (i=0; i<ida->nEq; i++)
    res[i] = Yprime[i] - RHS[i];

  free(RHS);
  if (!driver->ENERGY) free(temp);
  return 0;
}
/* -------------------------------------------------------------------------- */

void updateEnergy2D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                    flowStruct *flow, double dt, int k, double *Y, double *Yp)
{
  int j, r, row, col;

  /* update the current enthalpy to H + 1/rho*dPdt - q */
  
  r   = k*flame->n1;
  row = k*flame->n1*flame->nVars;

  for (j=0; j<flame->n1; j++) {
    col = j*flame->nVars;

    /* make sure the temperature is up-to-date */
    flame->T[r+j] = calcTfromH(chem->nSpecies, chem->nSpecSolve,
                               Y[row+col+flame->sT], &Y[row+col], 1000, chem);
    
    /* update the enthalpy */
    Yp[row+col+flame->sT] = flame->H[r+j] + 1.0/flame->rho[r+j]*flow->Dt_P*dt;
  } /* end loop over columns */
}
/* -------------------------------------------------------------------------- */

void calcExplicit_2D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                     flowStruct *flow, int iRow, double *Cp, double *Y)
{
  int i, j, jst, kst, col;
  
  // stencil bound for this row
  kst = stencil_y(iRow, flame);

  // loop over mesh points to solve the explicit terms
  for (j=0; j<flame->n1; j++)
    {
      col = j*flame->nVars;

      // loop over species and temperature
      for (i=0; i<flame->nVars; i++)
        {
          flame->dY.y[col+i]  = grad_y(kst,iRow,j,i,flame->nVars,Y,flame);
          flame->dY.yy[col+i] = grad_yy(kst,iRow,j,i,flame->nVars,Y,flame);
        }
      
      flame->dCp.y[j] = grad_y(kst, iRow, j, 0, 1, Cp, flame);
    } /* end loop over x1 */

  // compute the cross-derivative by taking the gradient dY.y
  for (j=0; j<flame->n1; j++)
    {
      col = j*flame->nVars;
      jst = stencil_x(j, flame);

      // loop over species and temperature
      for (i=0; i<flame->nVars; i++)
        flame->dY.xy[col+i] = grad_x(jst,0,j,i,flame->nVars,flame->dY.y,flame);
    } /* end loop over x1 */
}
/* -------------------------------------------------------------------------- */

/* compute the 1D solution along a line x1
   derivatives in x1 direction are implicit
   x2 direction is explicit with derivatives stored at time n
 */
void calcY_RHS_imp(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                   flowStruct *flow, int iRow, double *T, double *Y, double *RHS)
{
  int    i, j, jst, r, col;
  double wmix;                // mixture weight
  double Dx_Y, Dxx_Y;         // spatial derivatives
  double vT1, vT2, vH1, vH2;  // additional velocities in temperature eq.
  double src, diff, conv;     // source, diffusion, and convection terms

  // set some local pointers for readability
  int nSpecies = chem->nSpecSolve;
  double *cp   = chem->cp;
  
  r  = iRow*flame->n1;  // define row location indices
  // loop over mesh points to solve the species equations
  for (j = 0; j < flame->n1; j++)
    {
      col = j*flame->nVars;      // column index
      jst = stencil_x(j, flame); // start of stencil
      
      // compute local density and pressure source term
      wmix = calcMixtureW(nSpecies, &Y[col], chem->W);
      flame->rho[r+j] = wmix*flow->P/(flame->RGAS*T[j]);
      flame->srcP[j]  = 1.0/flame->rho[r+j]*flow->Dt_P;

      /* add additional terms implicitly to convective velocity if current
         direction is enthalpy */
      vH1 = vH2 = 0.0;
      if (driver->H_DIM != 0)
        {
          if (driver->DIM == driver->H_DIM && flow->D_phi1 != 0.0)
            vH1 = (flame->srcP[j] - flow->avgSdot[j])/flow->D_phi1;
          else if (driver->DIM != driver->H_DIM && flow->D_phi2 != 0.0)
            vH2 = (flame->srcP[j] - flow->avgSdot[iRow])/flow->D_phi2;
        }
      
      // get enthalpy and specific heat from chemical mechanism
      chemComputeThermoData(chem->h, chem->cp, T[j], chem);
      // get specific heat of the mixture
      flame->cpmix[j] = calcMixtureCp(nSpecies, chem->cp, &Y[col]);
      
      // compute and store concentrations of all species at local point
      for (i = 0; i < nSpecies; i++)
        chem->c[i] = flame->rho[r+j]*Y[col+i]/chem->W[i];
      
      // get rates from mechanism
      chemComputeProductionRates(chem->cdot, chem->w, chem->K, chem->c,
                                 chem->M, T[j], flow->P, chem);
      
      flame->srcH[j] = flame->sumCp[j] = flame->sumCp2[j] = 0.0;
      // loop over species
      for (i = 0; i < nSpecies; i++)
        {
          // calculate the gradients in normalised space
          Dx_Y  = grad_x(jst, 0, j, i, flame->nVars, Y, flame);
          Dxx_Y = grad_xx(jst, 0, j, i, flame->nVars, Y, flame);
          
          /* ----- compute diffusion term -----*/
          diff =   0.5*flame->chiA[r+j]*Dxx_Y
                 + flame->chiB[r+j]*flame->dY.xy[col+i]
                 + 0.5*flame->chiC[r+j]*flame->dY.yy[col+i];

          /* ----- compute convection term ----- */
          conv =   (flame->v1[r+j] + vH1)*Dx_Y
                 + (flame->v2[r+j] + vH2)*flame->dY.y[col+i];
          
          /* ----- compute chemical source term ----- */
          src  = chem->cdot[i]*chem->W[i]/flame->rho[r+j];
          
          /* ----- compute summation terms required for energy equation ----- */
          if (driver->ENERGY)
            flame->srcH[j] += chem->h[i]*Y[col+i];
          else // temperature
            {
              flame->sumCp[j]  += cp[i]*(flame->chiA[r+j]*Dx_Y + 
                                         flame->chiB[r+j]*flame->dY.y[col+i]);
              flame->sumCp2[j] += cp[i]*(flame->chiB[r+j]*Dx_Y + 
                                         flame->chiC[r+j]*flame->dY.y[col+i]);
              flame->srcH[j]   += chem->h[i]*src;  // for chemical source term
            }
          
          /* ----- update RHS of the species equations ----- */
          RHS[col+i] = diff - driver->CONV_MODE*conv + src;

          // store chemical source term for heat release computation
          flame->chemSrcY_tmp[(r+j)*flame->nVars+i] = src;
        } // end loop over species
    } // end loop over mesh
}
/* -------------------------------------------------------------------------- */

void calcT_RHS_imp(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                   flowStruct *flow, int iRow, double *Y, double *RHS)
{
  int    i, j, jst, r, col;
  double Dx_Y, Dxx_Y, Dx_cpmix; // spatial derivatives
  double vT1, vT2, vH1, vH2;    // additional velocities in temperature eq.
  double src, diff, conv;       // source, diffusion, and convection terms
  
  r  = iRow*flame->n1;  // current row location index
  /* ----- solve the temperature equation ----- */
  for (j = 0; j < flame->n1; j++)
    {
      col = j*flame->nVars + flame->sT;  // column index
      jst = stencil_x(j, flame);         // start of stencil
      
      // compute gradients in x1 direction
      Dx_Y     = grad_x(jst, 0, j, flame->sT, flame->nVars, Y, flame);
      Dxx_Y    = grad_xx(jst, 0, j, flame->sT, flame->nVars, Y, flame);
      Dx_cpmix = grad_x(jst, 0, j, 0, 1, flame->cpmix, flame);
      
      /* ----- compute diffusion term -----*/
      diff =   0.5*flame->chiA[r+j]*Dxx_Y
             + flame->chiB[r+j]*flame->dY.xy[col]
             + 0.5*flame->chiC[r+j]*flame->dY.yy[col];

      /* ----- compute convection term ----- */
      /* add additional terms implicitly to convective velocity if current
         direction is enthalpy */
      vH1 = vH2 = 0.0;
      if (driver->H_DIM != 0)
        {
          if (driver->DIM == driver->H_DIM && flow->D_phi1 != 0.0)
            vH1 = (flame->srcP[j] - flow->avgSdot[j])/flow->D_phi1;
          else if (driver->DIM != driver->H_DIM && flow->D_phi2 != 0.0)
            vH2 = (flame->srcP[j] - flow->avgSdot[iRow])/flow->D_phi2;
        }
      /* compute additional velocities for temperature equation (only if 
         not comparing to DNS solution) */
      vT1 = vT2 = 0.0;
      if (!driver->DNS)
        {
          vT1 = - 0.5/flame->cpmix[j]*(flame->chiA[r+j]*Dx_cpmix + 
                                       flame->chiB[r+j]*flame->dCp.y[j] +
                                       flame->sumCp[j]);
          vT2 = - 0.5/flame->cpmix[j]*(flame->chiB[r+j]*Dx_cpmix + 
                                       flame->chiC[r+j]*flame->dCp.y[j] +
                                       flame->sumCp2[j]);
        }

      // compute total contribution
      conv = (flame->v1[r+j] + vT1 + vH1)*Dx_Y +
             (flame->v2[r+j] + vT2 + vH2)*flame->dY.y[col];
      
      /* ----- compute chemical source term ----- */
      src  = 1.0/flame->cpmix[j]*(flame->srcP[j] - flame->srcH[j]);

      /* ----- copy sources over if needed ----- */
      if (driver->output_tsrc)
        {
          flame->diff_temp[r+j] = diff;
          flame->conv_temp[r+j] = conv;
          flame->chem_temp[r+j] = src;
        }
      
      /* ----- update RHS of temperature equation ----- */
      RHS[col] = diff - driver->CONV_MODE*conv + src;
    } // end loop over x1
}
/* -------------------------------------------------------------------------- */

void updateYpre(flameStruct *flame)
{
  int i, j, k;
  int row, col;
  
  // store the solution at current time
  for (k=0; k < flame->n2; k++)
    {
      row = k*(flame->n1*flame->nVars);
      for (j=0; j < flame->n1; j++)
        {
          col = j*flame->nVars;
          for (i=0; i < flame->nVars; i++)
            {
              flame->Y_m2[row+col+i]  = flame->Y_m1[row+col+i];
              flame->Y_m1[row+col+i]  = flame->Y[row+col+i];
              flame->Yp_m1[row+col+i] = flame->Yp[row+col+i];
            } // end loop over species
        } // end loop over columns
    } // end loop over rows
}
/* -------------------------------------------------------------------------- */

double calcMaxTimeStep(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                       flowStruct *flow, double CFL, double *cpmix, double *Y)
{
  int    i, j, k, jst, kst, r, row, col;
  double dx, dy;
  double Dx_Y, Dy_Y, Dx_cpmix, Dy_cpmix;     // spatial derivatives
  double diffA, diffB, diffC;                // diffusion coefficients
  double v1, v2, vT1, vT2, vH1, vH2;         // convective velocities
  double tmp, wmix, temp;

  int nSpecies = chem->nSpecSolve;  
  double dtMax = 1.0; // maximum allowable time step
  tmp = 0.0;

  for (k=0; k<flame->n2; k++)
    {
      row = k*flame->n1*flame->nVars; // row index
      r   = k*flame->n1;              // row index single variable
      kst = stencil_y(k, flame);      // start of stencil

      // compute the grid spacing
      if (k == 0)
        dy = 0.5*(flame->x2[k+1]-flame->x2[k]);
      else if (k == flame->n2-1)
        dy = 0.5*(flame->x2[k]-flame->x2[k-1]);
      else
        dy = 0.5*(flame->x2[k+1]-flame->x2[k-1]);
      
      // loop over mesh points to check the species equations
      for (j=0; j<flame->n1; j++)
        {
          col = j*flame->nVars;      // column index
          jst = stencil_x(j, flame); // start of stencil

          // compute the grid spacing
          if (j == 0)
            dx = 0.5*(flame->x1[j+1]-flame->x1[j]);
          else if (j == flame->n1-1)
            dx = 0.5*(flame->x1[j]-flame->x1[j-1]);
          else
            dx = 0.5*(flame->x1[j+1]-flame->x1[j-1]);
          
          // compute local density and pressure source term
          temp = Y[row+col+flame->sT];
          wmix = calcMixtureW(nSpecies, &Y[row+col], chem->W);
          flame->rho[r+j] = wmix*flow->P/(flame->RGAS*temp);
          flame->srcP[j]  = 1.0/flame->rho[r+j]*flow->Dt_P;

          // compute additional velocity terms for enthalpy direction
          vH1 = vH2 = 0.0;
          if (driver->H_DIM != 0)
            {
              if (driver->H_DIM == 1 && flow->D_phi1 != 0.0)
                vH1 = (flame->srcP[j] - flow->avgSdot[j])/flow->D_phi1;
              else if (driver->H_DIM == 2 && flow->D_phi2 != 0.0)
                vH2 = (flame->srcP[j] - flow->avgSdot[k])/flow->D_phi2;
            }

          /* ----- compute diffusion coefficients -----*/
          diffA = 0.5*flame->chiA[r+j];
          diffB = flame->chiB[r+j];
          diffC = 0.5*flame->chiC[r+j];

          tmp = fmax(diffA/(dx*dx), tmp);
          tmp = fmax(diffB/(dx*dy), tmp);
          tmp = fmax(diffC/(dy*dy), tmp);

          // loop over species
          flame->sumCp[j] = flame->sumCp2[j] = 0.0;
          for (i = 0; i < nSpecies; i++)
            {
              // compute summation terms required for temperature equation
              if (!driver->DNS)
                {
                  // get species gradients
                  Dx_Y = grad_x(jst, k, j, i, flame->nVars, Y, flame);
                  Dy_Y = grad_y(kst, k, j, i, flame->nVars, Y, flame);

                  // compute summation
                  flame->sumCp[j]  += chem->cp[i]*(flame->chiA[r+j]*Dx_Y +
                                                   flame->chiB[r+j]*Dy_Y);
                  flame->sumCp2[j] += chem->cp[i]*(flame->chiB[r+j]*Dx_Y + 
                                                   flame->chiC[r+j]*Dy_Y);
                }
              /* ----- compute convection velocity ----- */
              v1 = flame->v1[r+j] + vH1;
              v2 = flame->v2[r+j] + vH2;
              
              /* ----- check to see if this restricts the time step ----- */
              tmp = fmax(driver->CONV_MODE*fabs(v1)/dx, tmp);
              tmp = fmax(driver->CONV_MODE*fabs(v2)/dy, tmp);
            } /* end loop over species */
        } /* end loop over columns */
      
      /* ----- check the temperature equation ----- */	
      for (j = 0; j < flame->n1; j++)
        {
          jst = stencil_x(j, flame);  // start of stencil
          
          // compute the grid spacing
          if (j == 0)
            dx = 0.5*(flame->x1[j+1]-flame->x1[j]);
          else if (j == flame->n1-1)
            dx = 0.5*(flame->x1[j]-flame->x1[j-1]);
          else
            dx = 0.5*(flame->x1[j+1]-flame->x1[j-1]);

          /* ----- compute convective velocities ----- */
          // compute additional velocity terms for enthalpy direction
          vH1 = vH2 = 0.0;
          if (driver->H_DIM != 0)
            {
              if (driver->H_DIM == 1 && flow->D_phi1 != 0.0)
                vH1 = (flame->srcP[j] - flow->avgSdot[j])/flow->D_phi1;
              else if (driver->H_DIM == 2 && flow->D_phi2 != 0.0)
                vH2 = (flame->srcP[j] - flow->avgSdot[k])/flow->D_phi2;
            }
          // additional temperature contribution
          vT1 = vT2 = 0.0;
          if (!driver->DNS)
            {
              // get gradients of specific heat
              Dx_cpmix = grad_x(jst, k, j, 0, 1, cpmix, flame);
              Dy_cpmix = grad_y(kst, k, j, 0, 1, cpmix, flame);

              vT1 = - 0.5/cpmix[r+j]*(flame->chiA[r+j]*Dx_cpmix +
                                      flame->chiB[r+j]*Dy_cpmix + 
                                      flame->sumCp[j]);
              vT2 = - 0.5/cpmix[r+j]*(flame->chiB[r+j]*Dx_cpmix +
                                      flame->chiC[r+j]*Dy_cpmix +
                                      flame->sumCp2[j]);
            }

          v1 = flame->v1[r+j] + vT1 + vH1;
          v2 = flame->v2[r+j] + vT2 + vH2;

          /* ----- check to see whether this restricts the time step ----- */
          tmp = fmax(driver->CONV_MODE*fabs(v1)/dx, tmp);
          tmp = fmax(driver->CONV_MODE*fabs(v2)/dy, tmp);
        }  /* end loop over columns */
    }  /* end loop over rows */
  
  dtMax = fmin(dtMax, CFL/tmp);

  return dtMax;
}
/* -------------------------------------------------------------------------- */

// compute just for a 2D line for now
void solvePrestep2D(double **chi1, double **chi12, double **chi2,
                    driveStruct *driver, flameStruct *flame, flowStruct *flow)
{
  int    i, j, k, r, row, col;
  double Dphi1, Dphi2, Dphi12;
  double x1_sq, x2_sq;

  double *x1 = flame->x1;
  double *x2 = flame->x2;
  double *chi1_m1  = flow->chi1_m1;
  double *chi2_m1  = flow->chi2_m1;
  double *chi12_m1 = flow->chi12_m1;

  // compute and store inverse change in boundaries
  Dphi1 = Dphi2 = Dphi12 = 1.0;  
  if (flow->D_phi1 != 0.0)
    Dphi1 = flow->D_phi1;

  if (flow->D_phi2 != 0.0)
    Dphi2 = flow->D_phi2;

  if (flow->D_phi1 != 0.0 && flow->D_phi1 != 0.0)
    Dphi12 = flow->D_phi1*flow->D_phi2;
  
  // loop over grid
  for (k=0; k < flame->n2; k++)
    {
      row   = k*flame->n1*flame->nVars;
      r     = k*flame->n1;
      x2_sq = x2[k]*x2[k];

      for (j=0; j < flame->n1; j++)
        {
          col   = j*flame->nVars;
          x1_sq = x1[j]*x1[j];

          /* set the diffusion mode */
          chi1_m1[r+j]  *= driver->DIFF_MODE[0];
          chi12_m1[r+j] *= driver->DIFF_MODE[1];
          chi2_m1[r+j]  *= driver->DIFF_MODE[2];

          // clip cross-dissipation to ensure that it is correctly bounded
          if (chi1_m1[r+j] + 2*chi12_m1[r+j] + chi2_m1[r+j] < 0.0)
            chi12_m1[r+j] = -0.5*(chi1_m1[r+j] + chi2_m1[r+j]);

          // define transformed scalar dissipation rates and convective velocity
          if (driver->H_DIM != 0) // if enthalpy is a direction
            {
              if (j == 0 || j == flame->n1-1)
                flame->chiA[r+j] = 0.0; // enforce zero at left/right bounds
              else
                flame->chiA[r+j] = chi1_m1[r+j]/(Dphi1*Dphi1);

	      // Cross dissipation
              if (j == 0 || j == flame->n1-1)
                flame->chiB[r+j] = 0.0; // enforce zero at left/right bounds
              else if (k == 0 || k == flame->n2-1)
                flame->chiB[r+j] = 0.0; // enforce zero at top/bottom bounds
              else
                flame->chiB[r+j] = chi12_m1[r+j]/(Dphi1*Dphi2);

              if (k == 0 || k == flame->n2-1)
                flame->chiC[r+j] = 0.0; // enforce zero at top/bottom bounds
              else
                flame->chiC[r+j] = chi2_m1[r+j]/(Dphi2*Dphi2);

              /* NOTE: here the normalisation due to boundaries taken care of,
		 additional terms for enthalpy direction added later */
              flame->v1[r+j] = - 1.0/Dphi1*(flow->Dt_phi1min +
                      flame->x1[j]*(flow->Dt_phi1max - flow->Dt_phi1min) );
              flame->v2[r+j] = - 1.0/Dphi2*(flow->Dt_phi2min +
                      flame->x2[k]*(flow->Dt_phi2max - flow->Dt_phi2min) );
            }
          else if (driver->TRANS_Z2 == 1) // standard transform
            {
              if (j == 0 || j == flame->n1-1) // left/right boundaries
                {
                  flame->chiA[r+j] = flame->chiB[r+j] = 0.0;
                  flame->v1[r+j]   = flame->v2[r+j]   = 0.0;
                  if (k == 0 || k == flame->n2 - 1 || j == flame->n1-1)
                    flame->chiC[r+j] = 0.0;
                  else
                    flame->chiC[r+j] = chi2_m1[r+j];
                }
              else if (k == 0 || k == flame->n2-1) // top/bottom boundaries
                {
                  flame->chiB[r+j] = flame->chiC[r+j] = 0.0;
                  flame->v1[r+j]   = flame->v2[r+j]   = 0.0;
                  if (k == 0)
                    flame->chiA[r+j] = chi1_m1[r+j];
                  else
                    flame->chiA[r+j] = 0.0;
                }
              else // interior points
                {
                  flame->chiA[r+j] = chi1_m1[r+j];
                  flame->chiB[r+j] = (chi1_m1[r+j]*x2[k] + chi12_m1[r+j])/
                                     (1.0-x1[j]);
                  flame->chiC[r+j] = (chi1_m1[r+j]*x2_sq +
                                      2*chi12_m1[r+j]*x2[k] + chi2_m1[r+j])/
                                      ((1.0 - x1[j])*(1.0 - x1[j]));

                  flame->v1[r+j]   = 0.0;
                  flame->v2[r+j]   = - flame->chiB[r+j]/(1.0-x1[j]);
                }
            } // end standard transform
          else if (driver->TRANS_Z2 == 2) // Hasse transform
            {
              if (j == 0 || j == flame->n1-1) // left/right boundaries
                {
                  flame->chiA[r+j] = flame->chiB[r+j] = 0.0;
                  flame->chiC[r+j] = 0.0;
                  flame->v1[r+j]   = flame->v2[r+j]   = 0.0;
                }
              else if (k == 0 || k == flame->n2-1) // top/bottom boundaries
                {
                  flame->chiB[r+j] = flame->chiC[r+j] = 0.0;
                  flame->v1[r+j]   = flame->v2[r+j] = 0.0;
                  if (k == 0)
                    flame->chiA[r+j] = chi1_m1[r+j];
                  else if (k == flame->n2-1)
                    flame->chiA[r+j] = chi2_m1[r+j];
                }
              else // interior points
                {
                  flame->chiA[r+j] = chi1_m1[r+j] + chi2_m1[r+j];
                  flame->chiB[r+j] = ( - chi1_m1[r+j]*x2[k]
                                       + chi2[k][j]*(1.0-x2[k]) )/x1[j];
                  flame->chiC[r+j] = (chi1_m1[r+j]*x2_sq + chi2_m1[r+j]*
                                      (1.0-x2[k])*(1.0-x2[k]))/x1_sq;
                  flame->v1[r+j]   = 0.0;
                  flame->v2[r+j]   = flame->chiB[r+j]/x1[j];
                }
            } // end Hasse transform
          
          // store the structured chi in the flow struct for use in output
          flow->chi1[r+j]  = chi1[k][j];
          flow->chi12[r+j] = chi12[k][j];
          flow->chi2[r+j]  = chi2[k][j];

          // zero the variables for heat release computation
          flame->HRrate[r+j] = 0.0;
          for (i=0; i<flame->nVars; i++)
            flame->chemSrcY[row+col+i] = flame->chemSrcY_tmp[row+col+i] = 0.0;

          /* initialize sources to zero if needed*/
          if (driver->output_tsrc)
            {
              flame->chem_src[r+j] = 0.0;
              flame->conv_src[r+j] = 0.0;
              flame->diff_src[r+j] = 0.0;
            }
        } // end loop over x1
    } // end loop over x2
}
/* -------------------------------------------------------------------------- */

// compute just for a 2D line for now
void updateFlow2D(double t, int iRow, driveStruct *driver, timeStruct *time,
                  flameStruct *flame, flowStruct *flow)
{
  int    j, r;
  double Dphi1, Dphi2, Dphi12;
  double x1_sq, x2_sq, *chi1, *chi2, *chi12;

  // declare some local pointers for convenience
  double *x1  = flame->x1;
  double *x2  = flame->x2;
  double t_m1 = time->t_m1;
  double t_n  = time->t_n;

  r     = iRow*flame->n1;
  x2_sq = x2[iRow]*x2[iRow];

  // allocate local memory
  chi1 = (double*) calloc(flame->n1, sizeof(double));
  chi2 = (double*) calloc(flame->n1, sizeof(double));
  chi12 = (double*) calloc(flame->n1, sizeof(double));

  // compute and store inverse change in boundaries
  Dphi1 = Dphi2 = Dphi12 = 1.0;  
  if (flow->D_phi1 != 0.0)
    Dphi1 = flow->D_phi1;

  if (flow->D_phi2 != 0.0)
    Dphi2 = flow->D_phi2;

  if (flow->D_phi1 != 0.0 && flow->D_phi1 != 0.0)
    Dphi12 = flow->D_phi1*flow->D_phi2;

  // linearly interpolate scalar dissipation rate for this row at this time
  for (j=0; j<flame->n1; j++)
    {
      // interpolate the scalar dissipation to the current time
      chi1[j]  = linInterp(t, t_m1, flow->chi1_m1[r+j],  t_n, flow->chi1[r+j]);
      chi12[j] = linInterp(t, t_m1, flow->chi12_m1[r+j], t_n, flow->chi12[r+j]);
      chi2[j]  = linInterp(t, t_m1, flow->chi2_m1[r+j],  t_n, flow->chi2[r+j]);

      // apply the diffusion mode filter
      chi1[j]  *= driver->DIFF_MODE[0];
      chi12[j] *= driver->DIFF_MODE[1];
      chi2[j]  *= driver->DIFF_MODE[2];

      // clip cross-dissipation to ensure that it is correctly bounded
      if (chi1[j] + 2*chi12[j] + chi2[j] < 0.0)
        chi12[j] = -0.5*(chi1[j] + chi2[j]);
    } // end loop over columns
  
  // loop over grid
  for (j=0; j<flame->n1; j++)
    {
      x1_sq = x1[j]*x1[j];

      // if enthalpy is a direction
      if (driver->H_DIM != 0)
        {
          if (j == 0 || j == flame->n1-1)
            flame->chiA[r+j] = 0.0; // enforce zero at left/right bounds
          else
            flame->chiA[r+j] = chi1[j]/(Dphi1*Dphi1);

	  // Cross dissipation
          if (j == 0 || j == flame->n1-1)
            flame->chiB[r+j] = 0.0; // enforce zero at left/right bounds
          else if (iRow == 0 || iRow == flame->n2-1)
            flame->chiB[r+j] = 0.0; // enforce zero at top/bottom bounds
          else
            flame->chiB[r+j] = chi12[j]/(Dphi1*Dphi2);

          if (iRow == 0 || iRow == flame->n2-1)
            flame->chiC[r+j] = 0.0; // enforce zero at top/bottom bounds
          else
            flame->chiC[r+j] = chi2[j]/(Dphi2*Dphi2);
        }
      else if (driver->TRANS_Z2 == 1) // standard transform
        {
          if (j == 0 || j == flame->n1-1) // left/right boundaries
            {
              flame->chiA[r+j] = flame->chiB[r+j] = 0.0;
              flame->v1[r+j]   = flame->v2[r+j]   = 0.0;
              if (iRow == 0 || iRow == flame->n2 - 1 || j == flame->n1-1)
                flame->chiC[r+j] = 0.0;
              else
                flame->chiC[r+j] = chi2[j];
            }
          else if (iRow == 0 || iRow == flame->n2-1) // top/bottom boundaries
            {
              flame->chiB[r+j] = flame->chiC[r+j] = 0.0;
              flame->v1[r+j]   = flame->v2[r+j]   = 0.0;
              if (iRow == 0)
                flame->chiA[r+j] = chi1[j];
              else
                flame->chiA[r+j] = 0.0;
            }
          else // interior points
            {
              if (driver->DIM == 1)
                {
                  flame->chiA[r+j] = chi1[j];
                  flame->chiB[r+j] = (chi1[j]*x2[iRow]+chi12[j])/(1.0-x1[j]);
                  flame->chiC[r+j] = (chi1[j]*x2_sq + 2*chi12[j]*x2[iRow] +
                                      chi2[j])/((1.0 - x1[j])*(1.0 - x1[j]));
                  flame->v1[r+j]   = 0.0;
                  flame->v2[r+j]   = - flame->chiB[r+j]/(1.0-x1[j]);
                }
              else if (driver->DIM == 2)
                {
                  flame->chiA[r+j] = (chi2[j]*x1_sq + 2*chi12[j]*x1[j] +
                                      chi1[j])/
                                      ((1.0 - x2[iRow])*(1.0 - x2[iRow]));
                  flame->chiB[r+j] = (chi2[j]*x1[j]+chi12[j])/(1.0-x2[iRow]);
                  flame->chiC[r+j] = chi2[j];
                  flame->v1[r+j]   = - flame->chiB[r+j]/(1.0-x2[iRow]);
                  flame->v2[r+j]   = 0.0;
                }
            }
        } // end standard transform
      else if (driver->TRANS_Z2 == 2) // Hasse transform
        {
          if (j == 0 || j == flame->n1-1) // left/right boundaries
            {
              flame->chiA[r+j] = flame->chiB[r+j] = flame->chiC[r+j]=0.0;
              flame->v1[r+j]   = flame->v2[r+j]   = 0.0;
            }
          else if (iRow == 0 || iRow == flame->n2-1) // top/bottom boundaries
            {
              flame->chiB[r+j] = flame->chiC[r+j] = 0.0;
              flame->v1[r+j]   = flame->v2[r+j]   = 0.0;
              if (iRow == 0)
                flame->chiA[r+j] = chi1[j];
              else if (iRow == flame->n2-1)
                flame->chiA[r+j] = chi2[j];
            }
          else // interior points
            {
              if (driver->DIM == 1)
                {
                  flame->chiA[r+j] = chi1[j] + chi2[j];
                  flame->chiB[r+j] = ( - chi1[j]*x2[iRow]
                                       + chi2[j]*(1.0-x2[iRow]) )/x1[j];
                  flame->chiC[r+j] = (chi1[j]*x2_sq + chi2[j]*
                                      (1.0-x2[iRow])*(1.0-x2[iRow]))/x1_sq;
                  flame->v1[r+j]   = 0.0;
                  flame->v2[r+j]   = flame->chiB[r+j]/x1[j];
                }
              else if (driver->DIM == 2)
                {
                  flame->chiA[r+j] = (chi2[j]*x1_sq + chi1[j]*
                                      (1.0-x1[j])*(1.0-x1[j]))/x2_sq;
                  flame->chiB[r+j] = ( - chi2[j]*x1[j]
                                       + chi1[j]*(1.0-x1[j]))/x2[iRow];
                  flame->chiC[r+j] = chi1[j] + chi2[j];
                  flame->v1[r+j]   = flame->chiB[r+j]/x2[iRow];
                  flame->v2[r+j]   = 0.0;
                }
            }
        } // end Hasse transform
    } // end loop over x1

  free(chi1);
  free(chi2);
  free(chi12);
}
/* -------------------------------------------------------------------------- */

void calcCpmix2D(chemStruct *chem, flameStruct *flame, double *Y, double *cpmix)
{
  int i, j, k;
  int r, row, col;
  double temp;
  
  // store the solution at current time
  for (k=0; k < flame->n2; k++)
    {
      row = k*(flame->n1*flame->nVars);
      r   = k*flame->n1;
      for (j=0; j < flame->n1; j++)
        {
          col = j*flame->nVars;
          
          temp = Y[row + col + flame->sT]; // current temperature
          // get enthalpy and specific heat from chemical mechanism
          chemComputeThermoData(chem->h, chem->cp, temp, chem);
          // get specific heat of the mixture
          cpmix[r+j] = calcMixtureCp(chem->nSpecSolve, chem->cp, &Y[row+col]);
        } // end loop over columns
    } // end loop over rows
}
/* -------------------------------------------------------------------------- */

void setFlameCol(flameStruct *flame, flameStruct *flameCol)
{
  // variables that need to be transposed
  flameCol->n1      = flame->n2;
  flameCol->n2      = flame->n1;
  flameCol->x1      = flame->x2;
  flameCol->x2      = flame->x1;
  flameCol->dx1     = flame->dx2;
  flameCol->dx2     = flame->dx1;
  flameCol->chiA    = flame->chiC;
  flameCol->chiB    = flame->chiB;
  flameCol->chiC    = flame->chiA;
  flameCol->v1      = flame->v2;
  flameCol->v2      = flame->v1;
  flameCol->grad.x  = flame->grad.y;
  flameCol->grad.xx = flame->grad.yy;
  flameCol->grad.y  = flame->grad.x;
  flameCol->grad.yy = flame->grad.xx;

  // pointers that simply need to point to the same memory as rows
  flameCol->stw       = flame->stw;
  flameCol->order     = flame->order;
  flameCol->Y         = flame->Y;
  flameCol->Yp        = flame->Yp;
  flameCol->Y_m1      = flame->Y_m1;
  flameCol->Y_m2      = flame->Y_m2;
  flameCol->Yp_m1     = flame->Yp_m1;
  flameCol->nVars     = flame->nVars;
  flameCol->sT        = flame->sT;
  flameCol->sP        = flame->sP;
  flameCol->RGAS      = flame->RGAS;
  flameCol->cpmix     = flame->cpmix;
  flameCol->Mw        = flame->Mw;
  flameCol->R         = flame->R;
  flameCol->gamma     = flame->gamma;
  flameCol->rho       = flame->rho;
  flameCol->srcP      = flame->srcP;
  flameCol->srcH      = flame->srcH;
  flameCol->srcQ      = flame->srcQ;
  flameCol->sumCp     = flame->sumCp;
  flameCol->sumCp2    = flame->sumCp2;
  flameCol->HRrate    = flame->HRrate;
  flameCol->chemSrcY  = flame->chemSrcY;
  flameCol->chemSrcY_tmp  = flame->chemSrcY_tmp;
  flameCol->chem_src  = flame->chem_src;
  flameCol->conv_src  = flame->conv_src;
  flameCol->diff_src  = flame->diff_src;
  flameCol->chem_temp = flame->chem_temp;
  flameCol->conv_temp = flame->conv_temp;
  flameCol->diff_temp = flame->diff_temp;
  flameCol->dY.y      = flame->dY.y;
  flameCol->dY.yy     = flame->dY.yy;
  flameCol->dY.xy     = flame->dY.xy;
  flameCol->dCp.y     = flame->dCp.y;
}
/* -------------------------------------------------------------------------- */

void setFlowCol(flowStruct *flow, flowStruct *flowCol)
{
  // variables that need to be transposed
  flowCol->D_phi1     = flow->D_phi2;
  flowCol->D_phi2     = flow->D_phi1;
  flowCol->phi1min    = flow->phi2min;
  flowCol->phi1max    = flow->phi2max;
  flowCol->phi2min    = flow->phi1min;
  flowCol->phi2max    = flow->phi1max;
  flowCol->Dt_phi1min = flow->Dt_phi2min;
  flowCol->Dt_phi1max = flow->Dt_phi2max;
  flowCol->Dt_phi2min = flow->Dt_phi1min;
  flowCol->Dt_phi2max = flow->Dt_phi1max;
  flowCol->chi1       = flow->chi2;
  flowCol->chi2       = flow->chi1;
  flowCol->chi12      = flow->chi12;
  flowCol->chi1_m1    = flow->chi2_m1;
  flowCol->chi2_m1    = flow->chi1_m1;
  flowCol->chi12_m1   = flow->chi12_m1;

  // pointers that simply need to point to the same memory as rows
  flowCol->P          = flow->P;
  flowCol->P_m1       = flow->P_m1;
  flowCol->Dt_P       = flow->Dt_P;
  flowCol->avgSdot    = flow->avgSdot;
}
/* -------------------------------------------------------------------------- */
