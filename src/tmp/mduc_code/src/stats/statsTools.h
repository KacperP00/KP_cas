#ifndef MDUC_STATS_TOOLS_H
#define MDUC_STATS_TOOLS_H

#include "mducStats.h"

void runTest1D(optStruct *opts);
void errorTest1D(optStruct *opts);
void errorTest2D(optStruct *opts);

void regime(optStruct *opts);
void regime2D(optStruct *opts);

double getFileTime(FILE *pFile);

void readOptionsFile(optStruct *opts);

int readMomentFile(optStruct *opts, int **iter, double **time, double ***mean,
                   double ***var);

int interpPDFfromFile(char *fileName, int *n, double **x, double **pdf);

void writeStatsHeader(char *fileName, int dim, int *n, double **x, double t,
                      double *mean, double *var, int nPDF, char *transform);


void printStatsError(FILE *pFile, FILE *err, char *method, int dim, double *mean,
                     double *var, double *meanPDF, double *varPDF, double *fDiv);

#endif /* MDUC_STATS_TOOLS_H */
