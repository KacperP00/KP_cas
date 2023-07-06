#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "wiebe.h"
#include "chemStruct.h"

void ComputeProductionRatesWiebe(double *cdot, double *w, double *k, double *c,
                                 double *M, double temp, double pressure,
                                 chemStruct *chem)
{
  int    i;
  double src;
  double V, sumC;
  wiebeStruct *wb = (wiebeStruct*) chem->args[0]; /* Wiebe Parameters */

  /* compute the volume from knowing the initial number of moles */
  sumC = 0.0;
  for (i=0; i<chem->nSpecies; i++)
    sumC += c[i]; /* kmol/m^3 */
  sumC *= 1e3;    /* mol/m^3 */
  V    = wb->mol0/sumC;  /* m^3 */

  /* compute the Wiebe source term */
  src  = chemCalcWiebe_dt(wb->t, wb[0]); /* 1/CAD */
  src *= wb->tConv;                      /* 1/s */
  src *= wb->molF;                       /* kmol/s */
  src /= wb->V;                          /* kmol/(m^3-s) */

  /* compute the consumption rate of fuel */
  for (i=0; i<chem->nFuel; i++) {
    if (c[chem->sFuel[i]] > 0.0)
      cdot[chem->sFuel[i]] = -chem->fuelX[i]*src;
  }

  /* compute the consumption rate of O2 */
  cdot[chem->sO2]  = -(chem->fuelC+0.25*chem->fuelH-0.5*chem->fuelO)*src;

  /* compute the prodcution of products */
  cdot[chem->sCO2] = chem->fuelC*src;
  cdot[chem->sH2O] = 0.5*chem->fuelH*src;
  cdot[chem->sN2]  = 0.0;
}
/* ------------------------------------------------------------------------- */

int initWiebeSpecies(int *nA, char ***atoms, int *nSpMech, char ***species, 
                     chemStruct *chem)
{
  int i;

  /* assign the atoms */
  nA[0] = 4;
  atoms[0] = (char**) calloc(nA[0], sizeof(char*));
  for (i=0; i<nA[0]; i++)
    atoms[0][i] = (char*) calloc(3, sizeof(char));
  strcpy(atoms[0][0], "N");
  strcpy(atoms[0][1], "C");
  strcpy(atoms[0][2], "H");
  strcpy(atoms[0][3], "O");

  /* assign the species */
  if (chem->nFuel == 0) {
    fprintf(stderr, "Error: wiebe function requires fuel species to be set "
            "in input file\n  Use variable: fuel-species\n");
  }
    
  nSpMech[0] = chem->nFuel + 4;
  species[0] = (char**) calloc(nSpMech[0], sizeof(char*));
  for (i=0; i<nSpMech[0]; i++)
    species[0][i] = (char*) calloc(32, sizeof(char));
  for (i=0; i<chem->nFuel; i++)
    strcpy(species[0][i], chem->fuel[i]);
  strcpy(species[0][chem->nFuel],   "O2");
  strcpy(species[0][chem->nFuel+1], "N2");
  strcpy(species[0][chem->nFuel+2], "CO2");
  strcpy(species[0][chem->nFuel+3], "H2O");

  return nSpMech[0];
}
/* ------------------------------------------------------------------------- */

double chemCalcWiebe(double t, wiebeStruct w)
{
  if (t > w.t0)
    return 1.0-exp(-w.a*pow((t-w.t0)/w.dt, w.m+1));
  else
    return 0.0;
}
/* ------------------------------------------------------------------------- */

double chemCalcWiebe_dt(double t, wiebeStruct w)
{
  if (t > w.t0)
    return (1.0-chemCalcWiebe(t,w))*w.a*(w.m+1)/w.dt*pow((t-w.t0)/w.dt,w.m);
  else 
    return 0.0;
}
/* ------------------------------------------------------------------------- */
