#include <stdlib.h>
#include "sootHMOM.h"

void soothmom_(double *Temp, double *Viscosity, double *Density,
               double *MolarMass, double *C, double *moments,
               double *src, double *rxnrates, double *rhodot, double *ABSTOL)
{
  sootHMOM(*Temp,*Viscosity,*Density,*MolarMass,C,moments,src,rxnrates,rhodot,
           NULL,NULL,*ABSTOL);
}
void soothmom__(double *Temp, double *Viscosity, double *Density,
                double *MolarMass, double *C, double *moments,
                double *src, double *rxnrates, double *rhodot, double *ABSTOL)
{ soothmom_(Temp,Viscosity,Density,MolarMass,C,moments,src,rxnrates,rhodot,ABSTOL); }

double continuitysource_(double *rxnrates)
{ return ContinuitySource(rxnrates); }
double continuitysource__(double *rxnrates)
{ return continuitysource_(rxnrates); }

double fracmom_(double *VolOrd, double *SurfOrd, double *moments)
{ return FracMom(*VolOrd,*SurfOrd,moments); }
double fracmom__(double *VolOrd, double *SurfOrd, double *moments)
{ return fracmom_(VolOrd,SurfOrd,moments); }

double fracmomlarge_(double *VolOrd, double *SurfOrd, double *moments)
{ return FracMomLarge(*VolOrd,*SurfOrd,moments); }
double fracmomlarge__(double *VolOrd, double *SurfOrd, double *moments)
{ return fracmomlarge_(VolOrd,SurfOrd,moments); }



