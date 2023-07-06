#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "eosFluid.h"
#include "pengRobinson.h"
#include "helmholtzEOS.h"

eosFluid * fluidInit(char *name,char *eos,parseStruct *p)
{
  /* get the fluid properties from the parsed file */
  eosFluid *fluid = paramGetFluid(name,p);

  /* convert densities from (mol/dm^3) to (kg/m^3) */
  fluid->satTOL = 10;
  fluid->rhoc   = fluid->rhoMc*fluid->Mw;
  fluid->Ru     = 8.314472;
  fluid->R      = fluid->Ru/fluid->Mw*1e3;

  strcpy(fluid->eosName,eos);
  if (!strcmp(eos,"ideal")) {
    eosCalcP_rhoT = idealCalcP_rhoT;
    eosCalcT_rhoP = idealCalcT_rhoP;
    eosCalcRho_PT = idealCalcRho_PT;
  }
  else if (!strncmp(eos,"peng-rob",8)) {
    eosCalcP_rhoT = PR_calcP_rhoT;
    eosCalcT_rhoP = PR_calcT_rhoP;
    eosCalcRho_PT = PR_calcRho_PT;
    eosPR *pr = PR_init(fluid);
    fluid->eos = pr;
  }
  else if (!strcmp(eos,"helmholtz")) {
    eosCalcP_rhoT = HH_calcP_rhoT;
    eosCalcT_rhoP = HH_calcT_rhoP;
    eosCalcRho_PT = HH_calcRho_PT;
  }
  else {
    fprintf(stderr,"Warning: no equation of state set, using ideal gas\n");
    strcpy(fluid->eosName,"ideal");
    eosCalcP_rhoT = idealCalcP_rhoT;
    eosCalcT_rhoP = idealCalcT_rhoP;
    eosCalcRho_PT = idealCalcRho_PT;
  }

  return fluid;
}
/* -------------------------------------------------------------------------- */

double idealCalcP_rhoT(double rho, double T, eosFluid *fluid)
{
  return rho*fluid->R*T;
}
/* -------------------------------------------------------------------------- */

double idealCalcT_rhoP(double rho, double P, eosFluid *fluid)
{
  return P/(fluid->R*rho);
}
/* -------------------------------------------------------------------------- */

double idealCalcRho_PT(double P, double T, eosFluid *fluid)
{
  return P/(fluid->R*T);
}
/* -------------------------------------------------------------------------- */

double antoineP(double T, double *ant)
{
  return 1e5*pow(10,ant[0] - ant[1]/(T+ant[2]));
}
/* -------------------------------------------------------------------------- */

double antoineT(double P, double *ant)
{
  return ant[1]/(ant[0] - log10(P/1e5)) - ant[2];
}
/* -------------------------------------------------------------------------- */

void eosPrintFluid(eosFluid *f)
{
  fprintf(stderr," ===== Fluid Properties ===== \n");
  fprintf(stderr,"  Fluid:     %s\n",f->name);
  fprintf(stderr,"  Type:      %s\n",f->eosName);
  if (f->polar == 'p')
    fprintf(stderr,"  Polarity:  %s\n","polar");
  else
    fprintf(stderr,"  Polarity:  %s\n","non-polar");
  fprintf(stderr,"  Mw:        %10.3f [g/mol]\n",f->Mw);
  fprintf(stderr,"  Pc:        %10.3f [bar]\n",f->Pc*1e-5);
  fprintf(stderr,"  Tc:        %10.3f [K]\n",f->Tc);
  fprintf(stderr,"  rhoc:      %10.3f [mol/L]\n",f->rhoMc);
  fprintf(stderr,"             %10.3f [kg/m^3]\n",f->rhoc);
  if (f->antoine) {
    fprintf(stderr,"\n");
    fprintf(stderr,"  Antoine parameters:\n");
    fprintf(stderr,"   A:   %10.3f [-]\n",f->antoine[0]);
    fprintf(stderr,"   B:   %10.3f [-]\n",f->antoine[1]);
    fprintf(stderr,"   C:   %10.3f [-]\n",f->antoine[2]);
  }
  fprintf(stderr," ============================ \n");
}
/* -------------------------------------------------------------------------- */

eosFluid * paramGetFluid(char *name, parseStruct *p)
{
  char fName[1024];
  int iP,nP;
  char **raw,tag[64],buffer[256],*pS;

  sprintf(fName,"fluid-%s\0",name);
  
  iP = paramGetIndex(fName,p);
  if (iP < 0) {
    fprintf(stderr,"could not find fluid in input file");
    return NULL;
  }
  nP  = p->n[iP];
  raw = p->vals[iP];

  /* allocate the memory  */
  eosFluid *f = (eosFluid*) calloc(1, sizeof(eosFluid));

  /* set defaults */
  f->antoine = NULL;
  strcpy(f->name, name);

  f->Mw    = sparamGetReal(nP,"molecular-weight",-1.0,raw);
  f->Pc    = sparamGetReal(nP,"critical-pressure",-1.0,raw);
  f->Tc    = sparamGetReal(nP,"critical-temperature",-1.0,raw);
  f->rhoMc = sparamGetReal(nP,"critical-density",-1.0,raw);

  f->Tbp   = sparamGetReal(nP,"boiling-temperature",-1.0,raw);

  /* get Antoine parameters */
  int nA = sparamGetRealArray(nP,"Antoine",&f->antoine,raw);
  if (nA < 3) {
    fprintf(stderr,"Warning: not enough Antoine parameters (%i)\n",nA);
  }

  return f;
}
/* -------------------------------------------------------------------------- */

double eosCalcAcentricFac(double pSat, eosFluid *f)
{
  if (pSat < 0.0) {
    /* calculate? */
    if (f->antoine) {
      pSat = antoineP(0.7*f->Tc,f->antoine);
    }
    else {
      fprintf(stderr,"Warning: no Antoine parameters for calculating "
              "acentric factor\n");
      return 0.0;
    }
  }
  
  double Pr = pSat/f->Pc;
  return -log10(Pr) - 1.0;
}
/* -------------------------------------------------------------------------- */

void eosCopyState(eosState *state1, eosState *state2)
{
  state1->a       = state2->a;
  state1->u       = state2->u;
  state1->h       = state2->h;
  state1->s       = state2->s;
  state1->g       = state2->g;
  state1->cv      = state2->cv;
  state1->cp      = state2->cp;
  state1->w       = state2->w;
  state1->P       = state2->P;
  state1->T       = state2->T;
  state1->rho     = state2->rho;
  state1->dP_drho = state2->dP_drho;
  state1->dP_dT   = state2->dP_dT;
  state1->phi     = state2->phi;
}
/* -------------------------------------------------------------------------- */
