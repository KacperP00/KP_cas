#ifndef MDUC_GAUSS_PDF_H
#define MDUC_GAUSS_PDF_H

#include "statsUtils.h"

/* computes a 1D Gaussian pdf given a mean and variance.
     P = exp(a0 + a1*x + a2*x^2);
 */
int calcGaussPDF(int n, double *x, double mean, double var, double *pdf,
                 int *pdfBound);
int calcLogNormPDF(int n, double *x, double mean, double var, double *pdf,
                   int *pdfBound);

/* computes a 2D statistically-most-likely pdf given first five moments by
   maximising the entropy (shown above).  The resulting form is

     P = exp(a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y);
*/
int calcGaussPDF2D(int *n, double **x, double mean1, double var1, double mean2, 
                   double var2, double var, double **pdf, int **pdfBound,
                   int transform);

int calcLogNormPDF2D(int *n,double **x, double mean1, double var1, double mean2, 
                     double var2, double var, double **pdf, int **pdfBound,
                     int transform);

/* definition of distribution */
double gaussDist1D(double x, double a0, double a1, double a2, double a3);
double gaussDist2D(double x1, double x2, double a0, double a1, double a2,
                   double a3, double a4, double a5, double a6);


#endif // MDUC_GAUSS_PDF_H
