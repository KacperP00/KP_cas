#ifndef MDUC_TRANSPORT_H
#define MDUC_TRANSPORT_H

#include "chemStruct.h"

/* Calculate the viscosity from estimate basen on Chapman-Enskog Equation
     mu = 5/16*sqrt(pi*Mw*k_b*T)/(pi*sigma^2*Omega^(2,2))

   Variables:
     k_b           - Boltzmann's constant (1.3806488e-23)            [J/K]
     sigma         - Lennard-Jones collision diameter                [m]
     Omega(2,2)    - collision integral
     Mw            - species molecular weight                        [g/mol]
     T             - temperature                                     [K]

   Dependencies:
     eps_k         - ratio of L-J well depth/Boltzmann const.        [K]
*/
double calcViscosity(double Mw, double eps_k, double sigma, double T);

/* Computes the constant part of the viscosity approximation */
double calcViscCoeff(double Mw, double sigma);

/* ----- compute mixture averaged properties ----- */
double calcMixViscosity(int nSpecies, double T, double *Y, chemStruct *chem);

/* calculate the collision integral Omega(2,2)* using the quadratic
   interpolation based on Stockmayer potentials.
     Monchick, L. and Mason, E.A., Transport Properties of Polar Gases,
     Journal of Chemical Physics 35:1676 (1961)

   Inputs:
     Tr - reduced temperature, k_b*T/epsilon [-]
*/
double calcOmegaMu(double Tr);

/* calculate the collision integral Omega(1,1)* using the quadratic
   interpolation based on Stockmayer potentials.
     Monchick, L. and Mason, E.A., Transport Properties of Polar Gases,
     Journal of Chemical Physics 35:1676 (1961)

   Inputs:
     Tr - reduced temperature, k_b*T/epsilon [-]
*/
double calcOmegaD(double Tr);

#endif /* MDUC_TRANSPORT_H */
