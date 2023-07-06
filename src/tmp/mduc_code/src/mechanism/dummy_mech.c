#include <stdlib.h>
#include <stdio.h>

/* mechanism sizes */
int  GetNSpecies()   { return 0; }
int  GetNSpecs()     { return 0; }
int  GetNReactions() { return 0; }

/* dummy functions */
void GetSpeciesNames(char **names)
{
  fprintf(stderr, "Warning: compiled GetSpeciesName not present\n");
}
/* -------------------------------------------------------------------------- */

void GetMolarMass(double *W)
{
  fprintf(stderr, "Warning: compiled GetMolarMass not present\n");
}
/* -------------------------------------------------------------------------- */

void ComputeThermoData(double *h, double *cp, double T)
{
  fprintf(stderr, "Warning: compiled ComputeThermoData not present\n");
}
/* -------------------------------------------------------------------------- */

void ComputeProductionRates(double *cdot, double *w, double *k, double *c, 
                            double *M, double temp, double pressure)
{
  fprintf(stderr, "Warning: compiled ComputeProductionRates not pressent\n");
}
/* -------------------------------------------------------------------------- */
