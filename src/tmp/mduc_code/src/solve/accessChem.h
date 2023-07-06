#ifndef ACCESSCHEM_H
#define ACCESSCHEM_H

/* header file to provide access functions to chemical mechanism */

/* ----- return the total number of species in the mechanism ----- */
int GetNSpecies();

/* ----- return the number of non-steady-state species ----- */
int GetNSpecs();

/* ----- return the number of elementary reactions ----- */
int GetNReactions();

/* ----- return array of species properties ----- */
void   GetSpeciesNames(char **names);
void   GetMolarMass(double *W);
void   ComputeThermoData(double *h, double *cp, double T);
void   ComputeProductionRates(double *cdot, double *w, double *k, double *c, 
                              double *M, double temp, double pressure);
double GetLindRateCoeff(double temp, double pressure, double k0, double kInf, 
                        double fc, double conc);

#endif // ACCESSCHEM_H
