#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "pengRobinson.h"
#include "solveUtils.h"

eosPR * PR_init(eosFluid *f)
{
  eosPR * pr = (eosPR*) calloc(1, sizeof(eosPR));

  pr->nFl = 1;
  pr->f1  = f;

  /* compute the acentric factor */
  f->omega = eosCalcAcentricFac(-1.0, f);
  
  /* compute the constant for the correction */
  pr->kappa = 0.37464 + 1.54226*f->omega - 0.26992*f->omega*f->omega;

  return pr;
}
/* -------------------------------------------------------------------------- */

double PR_calcP_rhoT(double rho, double T, eosFluid *f)
{
  double a,b;
  double nu = f->Mw/rho;
  PR_calcParams(T,&a,&b,f);

  return f->Ru*1000.0*T/(nu-b) - a/(nu*(nu+b) + b*(nu-b));
}
/* -------------------------------------------------------------------------- */

double PR_calcT_rhoP(double rho, double P, eosFluid *fluid)
{

}
/* -------------------------------------------------------------------------- */

double PR_calcRho_PT(double P, double T, eosFluid *fluid)
{

}
/* -------------------------------------------------------------------------- */

/* double PR_brentRes_T_rhoP(double T, void *args) */
/* { */
/*   double a,b; */
/*   double nu = f->Mw/rho; */
/*   PR_calcParams(T,&a,&b,f); */

  
/* } */
/* -------------------------------------------------------------------------- */

void PR_calcParams(double T, double *a, double *b, eosFluid *fluid)
{
  eosPR *pr = (eosPR*)(fluid->eos);
  double ac,Tr,tmp,alpha;
  
  /* calculate the critical parameters */
  tmp  = fluid->Ru*fluid->Tc/fluid->Pc;
  ac   = 0.45724*fluid->Ru*fluid->Tc*tmp;
  b[0] = 0.07780*tmp;
  
  /* calculate the reduced temperatures for each phase */
  Tr = T/fluid->Tc;
  
  /* calculate the real correction to a */
  tmp    = 1.0 + pr->kappa*(1.0 - pow(Tr,0.5));
  alpha  = tmp*tmp;
  a[0]   = ac*alpha;
}
/* -------------------------------------------------------------------------- */

double PR_calcSaturation(double T, eosPR *pr)
{

}
/* -------------------------------------------------------------------------- */

double PR_fugacity(double T, eosPR *pr)
{

}
/* -------------------------------------------------------------------------- */

