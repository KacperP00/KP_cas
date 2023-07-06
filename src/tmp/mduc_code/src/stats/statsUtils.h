#ifndef MDUC_STATS_UTILS_H
#define MDUC_STATS_UTILS_H

#include <time.h>

/* define struct for passing information required during function evaluation */
struct STATSTRUCT {
  int    nMom;             // number of moments: 1D (3), 2D (6)
  int    n1, n2;           // size of grid in each direction
  int    transform;        // transformation to use
  double *x1, *x2;         // grid for each direction
  double **q;              // a priori pdf (for SML distribution)
  double mean[3], var[4];  // input means and variances
  double moments[6];       // ordered list of moments

  /* pointers to residual and Jacobian functions if desired */  
  void   (*fcn)(int, double *, double *, void *);
  void   (*jac)(int, double *, double **, void *);

#ifdef MDUC_JBETA_ADAPT
  /* parameters for adaptive pdf */
  double *a;
  double x,y,x0,y0,P_00,P_0y,dA;
  double xi, eta;
  int    bound;
  int    iX, nxG, *nyG;
  double *xG, *yG;
#endif /* MDUC_JBETA_ADAPT */
};
typedef struct STATSTRUCT statStruct;

/* include other files that define specific pdfs */
#include "betaPDF.h"
#include "gaussPDF.h"
#include "smlPDF.h"

/* ----- pointer to function for computing pdfs ----- */
int (*calcPDF_1D)(int, double *, double, double, double *, int *);
int (*calcPDF_2D)(int *, double **, double, double, double, double, double,
                  double **, int **, int);

/* ----- functions for clipping the moments and treating special cases ----- */
void clipMoments1D(double *mean, double *var);
void clipMoments2D(double *mean1, double *var1, double *mean2, double *var2,
                   double *var, double *covar);

void getPDFrowBounds(int n, double *x, double *bounds, double tol);

// functions for treating special cases of a pdf (e.g. delta functions)
int  checkPDF_1D(int n, double *x, double mean, double var, double *pdf,
                 int *pdfBound);

int  checkPDF_2D(int *n, double **x, double mean1, double var1, double mean2, 
                 double var2, double var, double **pdf, int **pdfBound,
                 int transform, int method);

/* ----- functions for computing the means ----- */
double calcLinMeanY(int index, int offset, double alpha, double beta, int n, 
                    double *x, double *Y, double *pdf, int *bounds, int pdfFlag);

void calcLinAllMeanY(int offset, double alpha, double beta, int n, 
                    double *x, double *Y,double *Ymean, double *pdf, int *bounds, int pdfFlag);

double calcLinMeanX(int index, int offset, double alpha, double beta, int n, 
                    double *x, double *Y, double *pdf, int *bounds, int pdfFlag,
                    int nSpecies, double *W);

double calcLinMeanY_2D(int index, int offset, int *n, double **x, double *Y,
                       double beta1, double beta2, double beta3, double **pdf,
                       int **pdfBound, int pdfFlag, int transform);

// for computing variances
double calcVar(int *n, int index, int offset, double mean1, double mean2, 
               double *Y1, double *Y2, double **pdf);

// compute the marginal pdf in direction <dir>
void calcMarginalPDF(int dir, int *n, double **x,double *marginal, double **pdf);

/* distance measures for comparison of two distributions */
double calcKullbackLeibler(int n1, int n2, double **P, double **Q);
double calcJensonShannon(int n1, int n2, double **P, double **Q);
double calcHellinger(int n1, int n2, double **P, double **Q);
double calcChiSqDiv(int n1, int n2, double **P, double **Q);
double calcTVD(int n1, int n2, double **P, double **Q);

/* ----- functions for binning pdf of a variable from data ----- */
int binPDF_1D(int nx, int ny, int nz, double ***var, int n1, double *x1,
              double *pdf);

int binPDF_2D(int nx, int ny, int nz, double ***var1, double ***var2,
              int *n, double **x, double **pdf);

/* ----- miscallaneous functions ----- */
// computes the pdf unscaled from the grid
void unscalePDF(int str, int transform, int n1, int n2, double *x1, double *x2,
                double **pdf);

// returns the value of a normal distribution given a mean and variance
double normalDistribution(double x, double mean, double var);

// calculate ln(Gamma) for a value xx: from Numerical Recipes in C, section 6-1
double gamma(double xx);
double gammaProd(int n, int d, double *num, double *den);
double lnGamma(double xx);

#endif // MDUC_STATS_UTILS_H
