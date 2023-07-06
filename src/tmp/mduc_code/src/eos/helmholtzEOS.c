#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "helmholtzEOS.h"
#include "solveUtils.h"

double HH_calcP_rhoT(double rho, double T, eosFluid *fluid)
{
  eosState state;

  HH_setStateRhoT(rho,T,fluid,&state);

  return state.P;
}
/* -------------------------------------------------------------------------- */

double HH_calcT_rhoP(double rho, double P, eosFluid *fluid)
{

}
/* -------------------------------------------------------------------------- */

double HH_calcRho_PT(double P,   double T, eosFluid *fluid)
{

}
/* -------------------------------------------------------------------------- */

void HH_setStateRhoT(double rho, double T, eosFluid *fluid, eosState *state)
{
  double tau, delta;
  double numTerm, denTerm; /* variables for temporary storage */
  double RT = fluid->R*T;
  
  double alpha;   /* dimensionless Helmholtz Energy */
  hhGrad alpha0;  /* ideal gas contribution */
  hhGrad alphaR;  /* real fluid contribution */
  
  /* ----- compute normalized temperature and density ----- */
  delta = rho/fluid->rhoc;
  tau   = fluid->Tc/T;
  
  /* ----- compute the derivative ideal and real contributions ----- */
  HH_calcAlphaIdeal(tau, delta, &alpha0, fluid);
  HH_calcAlphaReal(tau, delta, &alphaR, fluid);
  
  /* ----- total dimensionless Helmholtz energy ----- */
  alpha = alpha0.O + alphaR.O;
  
  /* ----- compute state properties ----- */
  /* some storage variables */
  numTerm = 1 +   delta*alphaR.d - delta*tau*alphaR.dt;
  denTerm = 1 + 2*delta*alphaR.d + delta*delta*alphaR.dd;
  
  /* compute non-dimensional first */
  state->Z   = 1 + delta*alphaR.d;
  state->u   = tau*(alpha0.t + alphaR.t);
  state->h   = state->u + state->Z;
  state->s   = state->u - alpha;
  state->g   = alpha + state->Z;
  state->cv  = -(tau*tau)*(alpha0.tt + alphaR.tt);
  state->cp  = state->cv + numTerm*numTerm/denTerm;
  state->w   = denTerm + numTerm*numTerm/state->cv;
  state->phi = exp(state->Z - 1 - log(state->Z) + alphaR.O);
  
  /* compute gradients */
  state->dP_drho = fluid->R*T*denTerm;
  state->dP_dT   = fluid->R*rho*numTerm;
  
  /* apply R and T to put in dimensional form */
  state->a   = alpha*RT;
  state->u   = state->u*RT;
  state->h   = state->h*RT;
  state->s   = state->s*fluid->R;
  state->g   = state->g*RT;
  state->cv  = state->cv*fluid->R;
  state->cp  = state->cp*fluid->R;
  state->w   = sqrt(state->w*RT);
  state->P   = rho*RT*state->Z;
  state->T   = T;
  state->rho = rho;
}
/* -------------------------------------------------------------------------- */

void HH_setStateTP(double T, double P, eosFluid *fluid, eosState *state)
{
}
/* -------------------------------------------------------------------------- */

void HH_setStateRhoP(double rho, double P, eosFluid *fluid, eosState *state)
{
}
/* -------------------------------------------------------------------------- */

void HH_calcAlphaIdeal(double tau,double delta,hhGrad *alpha0,eosFluid *fluid)
{
  int k;
  double test;
  double expTerm, t_plus;

  /* ----- compute the ideal contribution to the Helmholtz energy ----- */
  alpha0->O =    fluid->a[0] + fluid->a[1]*tau + log(delta)
              + (fluid->c[0]-1)*log(tau);
  if (fluid->c[1] > 0 && fluid->c[2] > 0)
    alpha0->O -= fluid->c[1]*pow(fluid->Tc,fluid->c[2])/
                 (fluid->c[2]*(fluid->c[2]+1))*pow(tau,-fluid->c[2]);
  alpha0->t  =   fluid->a[1] + (fluid->c[0]-1)/tau;
  alpha0->tt = -(fluid->c[0]-1)/(tau*tau);

  /* loop over terms with Einstein coefficients */
  for (k=0; k<4; k++) {
    if (fluid->u[k] > 0) {
      expTerm     = exp(-fluid->u[k]*tau/fluid->Tc);
      t_plus      = fluid->v[k]*fluid->u[k]/fluid->Tc*expTerm/(1-expTerm);
      alpha0->O  += fluid->v[k]*log(1-expTerm);
      alpha0->t  += t_plus;
      alpha0->tt -= t_plus*fluid->u[k]/fluid->Tc*(1 + expTerm/(1-expTerm));
    }
  }
}
/* -------------------------------------------------------------------------- */

void HH_calcAlphaReal(double tau,double delta,hhGrad *alphaR,eosFluid *fluid) 
{
  int i;
  double O_plus, t_plus, d_plus;
  alphaR->O   = 0;
  alphaR->t   = 0;
  alphaR->tt  = 0;
  alphaR->d   = 0;
  alphaR->dd  = 0;
  alphaR->ddd = 0;
  alphaR->dt  = 0;
  
  /* ----- compute real gas contribution ----- */
  /* loop over terms with no exponential */
  for (i=0; i<fluid->nPol; i++) {
    O_plus      =  fluid->n[i]*pow(delta,fluid->dExp[i])*pow(tau,fluid->tExp[i]);
    t_plus      =  O_plus*fluid->tExp[i]/tau;
    d_plus      =  O_plus*fluid->dExp[i]/delta;
    alphaR->O   += O_plus;
    alphaR->t   += t_plus;
    alphaR->tt  += t_plus*(fluid->tExp[i]-1)/tau;
    alphaR->d   += d_plus;
    alphaR->dd  += d_plus*(fluid->dExp[i]-1)/delta;
    /*alphaR->ddd += alphaR->dd*(fluid->dExp[i]-2)/delta; */
    alphaR->dt  += t_plus*fluid->dExp[i]/delta;
  }
  
  /* loop over terms with exponential */
  for (i=fluid->nPol; i<fluid->nT; i++) {
    O_plus      = fluid->n[i]*pow(delta,fluid->dExp[i])*pow(tau,fluid->tExp[i])*
                  exp(-pow(delta,fluid->powExp[i]));
    t_plus      = O_plus*fluid->tExp[i]/tau;
    d_plus      = O_plus/delta*(fluid->dExp[i] - 
                                fluid->powExp[i]*pow(delta,fluid->powExp[i]));
    alphaR->O  += O_plus;
    alphaR->t  += t_plus;
    alphaR->tt += t_plus*(fluid->tExp[i]-1)/tau;
    alphaR->d  += d_plus;
    alphaR->dd += d_plus/delta*(fluid->dExp[i] - 
                                fluid->powExp[i]*pow(delta,fluid->powExp[i])-1);
    alphaR->dd -= O_plus*pow(fluid->powExp[i]/delta,2)*
                         pow(delta,fluid->powExp[i]);
    alphaR->dt += d_plus*fluid->tExp[i]/tau;
  }
}
/* -------------------------------------------------------------------------- */

void HH_calcAlphaD(double tau,double delta,hhGrad *alphaR,eosFluid *fluid)
{
  int i;
  double O_plus, d_plus;
  alphaR->O   = 0;
  alphaR->t   = 0;
  alphaR->tt  = 0;
  alphaR->d   = 0;
  alphaR->dd  = 0;
  alphaR->ddd = 0;
  alphaR->dt  = 0;
  
  /* ----- compute real gas contribution ----- */
  /* loop over terms with no exponential */
  for (i=0; i<fluid->nPol; i++) {
    O_plus      =  fluid->n[i]*pow(delta,fluid->dExp[i])*pow(tau,fluid->tExp[i]);
    d_plus      =  O_plus*fluid->dExp[i]/delta;
    alphaR->O   += O_plus;
    alphaR->d   += d_plus;
    alphaR->dd  += d_plus*(fluid->dExp[i]-1)/delta;
  }
  
  /* loop over terms with exponential */
  for (i=fluid->nPol; i<fluid->nT; i++) {
    O_plus      = fluid->n[i]*pow(delta,fluid->dExp[i])*pow(tau,fluid->tExp[i])*
                  exp(-pow(delta,fluid->powExp[i]));
    d_plus      = O_plus/delta*(fluid->dExp[i]-
                                fluid->powExp[i]*pow(delta,fluid->powExp[i]));
    alphaR->O  += O_plus;
    alphaR->d  += d_plus;
    alphaR->dd += d_plus/delta*(fluid->dExp[i] -
                                fluid->powExp[i]*pow(delta,fluid->powExp[i])-1);
    alphaR->dd -= O_plus*pow(fluid->powExp[i]/delta,2)*
                  pow(delta,fluid->powExp[i]);
  }
}
/* -------------------------------------------------------------------------- */


void HH_calcSatT(double T,eosFluid *fluid,eosState *stateL,eosState *stateV)
{
  int its, maxIts = 100;
  int check;
  double TOL = 1e-4;
  double delta[2];
  double temp;
  eosSat sat;
  
  /* set delta at for liquid (0) and vapour (1) */
  delta[0] = fluid->rhoMax/fluid->rhoc;
  delta[1] = 0.1/fluid->rhoc;
  
  /* initialise arguments for function */
  sat.T     = T;
  sat.fluid = fluid;
  
  if (T > fluid->Tc) {
    fprintf(stderr,"T = %g > Tc = %g\n", T, fluid->Tc);
    fprintf(stderr,"Fluid is supercitical! - no saturation solution\n");
    exit(1);
  }
	
  /* ensure that the initial guess is not too close to the critical point */
  if (T > 0.95*fluid->Tc) {
    temp = 0.92*T;
    sat.T = temp;
    newt(delta, 2, &check, &sat, HH_resSatT, HH_jacSatT);
  }
  
  sat.T = temp;
  newt(delta, 2, &check, &sat, HH_resSatT, HH_jacSatT);
  
  /* compute the liquid and vapour states */
  HH_setStateRhoT(delta[0]*fluid->rhoc, T, fluid, stateL);
  HH_setStateRhoT(delta[1]*fluid->rhoc, T, fluid, stateV);
  its = 1;
  while (fabs(stateL->P-stateV->P) > fluid->satTOL) {
    TOL = TOL/10;
    newt(delta, 2, &check, &sat, HH_resSatT, HH_jacSatT);
    HH_setStateRhoT(delta[0]*fluid->rhoc, T, fluid, stateL);
    HH_setStateRhoT(delta[1]*fluid->rhoc, T, fluid, stateV);
    its++;
    if (its > maxIts) {
      fprintf(stderr,"Saturation pressure did not converge within satTOL\n");
      exit(1);
    }
  }
}
/* ------------------------------------------------------------------------- */

void HH_calcSatP(double P,eosFluid *fluid,eosState *stateL,eosState *stateV)
{
  double Tsat;
  eosSat sat;
  /* eosSat *sat = (eosSat*) calloc(1,sizeof(eosSat)); */
  
  if (P >= fluid->Pc) {
    fprintf(stderr,"P = %e > Pc = %e\n",P,fluid->Pc);
    fprintf(stderr,"Fluid is supercritical - no saturation solution\n");
    exit(1);
  }

  sat.fluid = fluid;
  sat.liq   = stateL;
  sat.vap   = stateV;
  sat.P     = P;
  
  /* !!!!! There should be a call to an Antoine solver here for initial guess */
  Tsat = antoineT(P,fluid->antoine);
  Tsat = brent(HH_brentSatP, Tsat-10, Tsat+10, 1e-2, &sat);
  HH_calcSatT(Tsat, fluid, stateL, stateV);
  /* free(sat); */
}
/* -------------------------------------------------------------------------- */

void HH_calcSatRho(double rho,eosFluid *fluid,eosState *state)
{
  int **ARGS;
  int phase;
  double Tsat;
  /* eosSat *sat = (eosSat*) calloc(1,sizeof(eosSat)); */
  eosSat   sat;
  eosState stateL, stateV;
  
  sat.fluid = fluid;
  sat.liq   = &stateL;
  sat.vap   = &stateV;
  /* TODO: phase is not defined yet!*/
  /*sat.phase = phase;*/
  sat.rho   = rho;
  
  if(rho > fluid->rhoc) { /* liquid */
    phase = 1;
    Tsat = brent(HH_brentSatRho, fluid->Ttp+1, fluid->Tc-1, 1e-2, &sat);
    HH_calcSatT(Tsat, fluid, &stateL, &stateV);
    eosCopyState(state, &stateL);
  }
  else {
    phase = 0;
    Tsat = brent(HH_brentSatRho, fluid->Ttp+1, fluid->Tc-1, 1e-2, &sat);
    HH_calcSatT(Tsat, fluid, &stateL, &stateV);
    eosCopyState(state, &stateV);
  }
}
/* -------------------------------------------------------------------------- */

double HH_brentSatP(double T, void *args)
{
  /* dereference parameters */
  eosSat *sat = (eosSat*) args;
  
  HH_calcSatT(T,sat->fluid,sat->liq,sat->vap);
  return 0.5*(sat->liq->P + sat->vap->P) - sat->P;
}
/* -------------------------------------------------------------------------- */

double HH_brentSatRho(double T, void *args)
{
  /* dereference parameters */
  eosSat *sat = (eosSat*) args;
  
  HH_calcSatT(T,sat->fluid,sat->liq,sat->vap);
  if (sat->phase == 1)
    return sat->liq->rho - sat->rho;
  else
    return sat->vap->rho - sat->rho;
}
/* -------------------------------------------------------------------------- */

void HH_resSatT(int n,double *delta,double *res, void *args)
{
  double tau;

  /* dereference parameters */
  eosSat *sat = (eosSat*) args;
  
  /* [0] = liquid, [1] = vapour */
  hhGrad alpha[2];
  tau = sat->fluid->Tc/(sat->T);
  
  /* calculate derivatives of alpha */
  HH_calcAlphaD(tau, delta[0], &alpha[0], sat->fluid);
  HH_calcAlphaD(tau, delta[1], &alpha[1], sat->fluid);
  
  /* compute the residual of the system */
  res[0] = delta[0]*(1+delta[0]*alpha[0].d) - delta[1]*(1+delta[1]*alpha[1].d);
  res[1] = alpha[0].O - alpha[1].O + log(delta[0]/delta[1]) + 
           delta[0]*alpha[0].d - delta[1]*alpha[1].d;
}
/* ------------------------------------------------------------------------- */

void HH_jacSatT(int n,double *delta,double **jac,void *args)
{
  double tau;

  /* dereference parameters */
  eosSat *sat = (eosSat*) args;
  
  /* [0] = liquid, [1] = vapour */
  hhGrad alpha[2];
  tau = sat->fluid->Tc/(sat->T);
  
  /* calculate derivatives of alpha */
  HH_calcAlphaD(tau, delta[0], &alpha[0], sat->fluid);
  HH_calcAlphaD(tau, delta[1], &alpha[1], sat->fluid);
  
  /* compute Jacobian */
  jac[0][0] =  1 + 2*delta[0]*alpha[0].d + pow(delta[0],2)*alpha[0].dd;
  jac[0][1] = -1 - 2*delta[1]*alpha[1].d - pow(delta[1],2)*alpha[1].dd;
  jac[1][0] =  1/delta[0] + delta[0]*alpha[0].dd + 2*alpha[0].d;
  jac[1][1] = -1/delta[1] - delta[1]*alpha[1].dd - 2*alpha[1].d;
}
/* ------------------------------------------------------------------------- */
