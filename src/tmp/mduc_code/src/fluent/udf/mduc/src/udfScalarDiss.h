#ifndef MDUC_UDF_SCALAR_DISS_H
#define MDUC_UDF_SCALAR_DISS_H

#include "mducUDF.h"

/* ----- initialisation functions ----- */
void UDFinitScalarDiss(udfStruct *udf);
int  UDFsetChiMode(char *mode);
int  UDFsetChiForm(char *form);
int  UDFsetChiScaling(char *scaling);

#if !RP_HOST
/* ----- functions for calculating scalar dissipation rate ----- */
/* compute 1D chi based on normalised distribution */
void UDFcalcChi1D(Domain *d, int n1, double *x1, double *chi);

/* do simple binning retrieve scalar dissipation */
void UDFbinChi1D(Domain *d, int n1, double *x1, double *chi);

/* compute 1D chi profile using a reference value computed from domain
   and a functional form */
void UDFcalcChiPitsch(Domain *d, int n1, double *x1, double Z1ref, double *chi);

/* compute 2D chi based on normalised distribution */
void UDFcalcChi2D(Domain *d, int *nMesh, double **xMesh, double **chi1,
                  double **chi12, double **chi2);
void UDFcalcChi2DAvg(Domain *d, int *nMesh, double **xMesh, double **chi1,
                     double **chi12, double **chi2);

/* do weighted binning in 2D space based on the mean values */
void UDFbinChi2D(Domain *d, int n1, int n2, double *x1, double *x2,
                 double **chi1, double **chi2);
#endif /* !RP_HOST */

/* functional forms for chi computations */
double (*UDFcalcChiForm)(double, double, double);

#endif /* MDUC_UDF_SCALAR_DISS_H */
