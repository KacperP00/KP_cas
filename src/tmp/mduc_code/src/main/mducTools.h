#ifndef MDUC_TOOLS_H
#define MDUC_TOOLS_H

#include "mducMain.h"
#include "mducStructs.h"

double calcLiquidFuelTemp(double T_L, double h_v, double cp_L, double T_ref,
                          int nSpecies, int nFuel, int *fuelIndex, double *Y,
                          chemStruct *chem);

void interpolateSolution(optStruct *opts);
void convert1D_2D(optStruct *opts);

/* ----- functions for converting data to ensight format ----- */
/* create an ensight formatted time series from a monitor file */
void genEnsightTime(optStruct *opts);
/* convert flamelet files to ensight case format */
void mduc2ensight(optStruct *opts);
/* generate a ensight constant file from a monitor file */
void writeEnsightConst(optStruct *opts, int n, int freq);
int  readEnsightCaseTime(char *fileName, double **time);
void convertEnsightParticle(optStruct *opts);

// compute mean mixture fractions given a mass-averaged pdf from flow file
void getZfromPDF(optStruct *opts);

void cleanMonitorFile(optStruct *opts);
void combineMonitorFiles(optStruct *opts);
int  readMonitorFile(char *fileName, int iT, double **time, int iV,double **var);
void genEnsightTime(optStruct *opts);
void checkGrid(optStruct *opts);

void chiProfile(optStruct *opts);
void checkChi2D(optStruct *opts, void *mem);
void postProcessChi(optStruct *opts);

void heatReleaseBudget(optStruct *opts);

void transformFlow(optStruct *opts);
void testMetric(optStruct *opts);

void checkThermo(optStruct *opts, chemStruct *chem, flameStruct *flame);

void genFireSpeciesList(optStruct *opts);
void integrateHeatRelease(optStruct *opts);

void testPolyA();

void flameletFileAddWmol(optStruct *opts);

#endif /* MDUC_TOOLS_H */
