#ifndef MDUC_STATS_H
#define MDUC_STATS_H

/* standard c libraries */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <string.h>

#include "mducStructs.h"

/* struct to hold option read from command line */
struct OPTSTRUCT {
  // flags
  int nOpts;
  char *args;
  char file[100], optFile[100], momentFile[100], errorFile[100];
  char **method;
  char transformName[10];
  char ident[100];
  char prefix[100];
  char **names;

  int dim, nMethod;
  int runTest, errorTest;
  int inputMean, inputBeta, inputFile, interpFile, readOptFile;
  int transform, scaled;
  int verbose, regime;
  int dataSet;
  int meanLin; // flag for computing mean with linear boundaries
  int extract; // extract a pdf from a data file
  int range[2];

  // grid parameters
  int    n1, n2, nPiv1, nPiv2;
  double *gPiv1, *gPerc1, *gPiv2, *gPerc2;

  // variables to read in
  double mean1, mean2;
  double var1, var2, var12;
  double max1, max2;

  double beta1, beta2, beta3;
  
};
typedef struct OPTSTRUCT optStruct;

void readList(int n, int iStart, char ***var, char **argv);
void readOptions(int argc, char **argv, optStruct *opts);

/* ----- miscellaneous functions ----- */

#endif /* MDUC_STATS_H */
