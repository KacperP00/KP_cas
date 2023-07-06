#ifndef MDUC_UDF_TOOLS_H
#define MDUC_UDF_TOOLS_H

#include "mducUDF.h"
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define SQR(c)    ((c) * (c))

/* initialisation functions */
void initFluentOptions(udfStruct *udf);
int  RP_is_Defined(char *var_name);
void readFluentString(int n, char *input, char ***species);
int  parseFluentString(char *input, char ***var);
void UDFstripStringBraces(char *name);
int  parseFluentReal(char *input, double **var);
void readFluentReal(int n, char *input, double **val);


int    UDFgetThreads(int n, char **names, Thread **threads);

#if !RP_HOST
/* ----- domain operations ----- */
/* general functions for extracting flow properties */
void   initThreadSolve(Domain *d);
void   checkThreadSolve(Domain *d);
void   UDFupdateSolverThreadNames(Domain *d);
double calcAvgPressure(Domain *d);
double calcAvgTemp(Domain *d);
void   UDFupdateEnthalpy(Domain *d, void *mem);
void   clipVariance(Domain *d);
void   clipEnthalpyVariance(Domain *d);
void   getPhiBounds(Domain *d);
void   getEnthalpyBounds(Domain *d);

/* functions for changing domain properties */
void adjustMass(double mSet, Domain *d);
void adjustTemp(double Tset, Domain *d);

void UDFwritePDF_2D(Domain *d, int local);

#endif /* !RP_HOST */


/* I/O functions */
void write1DArray(char *fileName,char *header,int n1,double *x1,double *var);
void write2DArray(char *fileName, char *header, int n1, int n2, double *x1, 
                  double *x2, double **var);

/* math functions */
int    calcFactorial(int x);

/* root solve functions */
int rootBrac(double (*func)(double, void *), double *x1, double *x2, void *args);
int rootBrak(double (*func)(double, void *), double x1, double x2, int n, 
             double *xb1, double *xb2, int *nb, void *args);

double rootBisect(double (*func)(double, void *), double x1, double x2,
                  double tol, void *args);

double rootBrent(double (*funcEval)(double, void *), double x1, double x2,
                 double tol, void *args);

/* sort an array of integers and return the index */
void indexIntArray(int n, int *arr, int *index);

/* miscellaneous */
int UDFgridLocate(int n1, double *x1, double x);
double MAX_D(double x, double y);
double MIN_D(double x, double y);
double ABS_D(double x);
int MAX_I(int x, int y);
int MIN_I(int x, int y);
double SGN_D(double x);

#endif /* MDUC_UDF_TOOLS_H */
