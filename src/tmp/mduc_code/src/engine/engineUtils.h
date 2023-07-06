#ifndef MDUC_ENGINE_TOOLS_H
#define MDUC_ENGINE_TOOLS_H

#include "engineStruct.h"

/* conversion between crank angle and absolute time */
double cad_to_time(double CAD, engineStruct *engine);
double time_to_cad(double t, engineStruct *engine);

/* functions for computing the cylinder geometry */
double calcCylLength(double CAD, engineStruct *engine);
double calcPistonVelocity(double CAD, engineStruct *engine);
double calcCylVol(double cad, engineStruct *engine);

/* calculate the valve area for a given lift */
double calcValveArea(double lift, valveStruct *v);

/* calculate mass flow rate through valve assuming isentropic steady flow */
double calcValveSteadyMassFlow(double gamma, double rhoi, double Pu, double Pd,
                               double Af);

double calcEngine_dPdt(double t1, double t2, double P1, double gamma,
                       engineStruct *eng);

/* heat transfer models */
double calcWoschniQdot(double t, double Pc, double Tc, double gamma, 
                       engineStruct *engine);

/* functions for estimating EGR composition give equivalence ratio inputs */
void calcEGR(double phi, double psi, double EGR, double M_fuel, int nC, int nH,
             int nO, double *nb, double *nu);
double calcEGR_K(double Tb);

double calcFlameSpeedHeywood(char *f,double phi,double xb,double P,double T);

#endif /* MDUC_ENGINE_TOOLS_H */
