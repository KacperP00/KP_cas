#ifndef UDF_SPRAY_PROFILE_H
#define UDF_SPRAY_PROFILE_H

/* include standard libraries */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* include Fluent required libraries */
#include "udf.h"
#include "prf.h"

typedef struct UDF_DPM_STRUCT {
  real penetration;
  real D32;
};
typedef struct UDF_DPM_STRUCT dpmStruct;

/* function to read the profile data file */
int DPM_read_profile(char *fileName, double **t, double **var);

/* utility functions for getting info about the spray */
real DPM_integrate_profile(int n, real *t, real *data);
void DPM_print_profile(char *name, int n, real *t, real *data);

/* functions for post-processing data from the spray */
real DPM_penetration_length(Injection *I, real *max_pen, real *pen_05);
real DPM_get_SauterD(Injection *I);

/* simple functions */
int UDFsprayParseFluentString(char *input, char ***var);
void UDFsprayPrefixString(char *pre, char *str);
double MAX_D(double x, double y);
double MIN_D(double x, double y);
double UDF_lnGamma(double xx);

#endif /* UDF_SPRAY_PROFILE_H */
