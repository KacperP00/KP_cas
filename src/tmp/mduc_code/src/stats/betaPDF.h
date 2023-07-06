#ifndef MDUC_BETA_PDF_H
#define MDUC_BETA_PDF_H

#include "statsUtils.h"

#ifdef MDUC_JBETA_ADAPT
#include "sdialStruct.h"
#endif /* MDUC_JBETA_ADAPT */

/* compute a beta distribution given a mean and a variance
   inputs
     mean, variance
 
   f(x; a,b) = G(a+b)/(G(a)*G(b))*x^(a-1)*(1-x)^(b-1)
 
   The value of f(z) is returned, with (a,b) being defined using the mean
   and variance as:
 
    a = mean*[mean*(1-mean)/var - 1]
    b = (1-mean)*[mean*(1-mean)/var - 1]
 
   NOTE: The gamma function is computed as the log of the gamma function,
         thus G(a+b)/[G(a)+G(b)] == exp(ln[G(a+B)] - ln[G(a)] - ln[G(b)] )
 */
int calcBetaPDF(int n, double *x, double mean, double var, double *pdf,
                int *pdfBound);

int calcBetaPDFall(int n, double *x, double mean, double var, double *pdf,
                int *pdfBound);

void betaDistAll(double *x, double *x1work, double*x2work, double alpha, double beta, int n, double *result);
/* function to compute a Dirichlet distribution for mixing of three scalars
   inputs:
     xmean, ymean
     xvar, yvar, xyvar
 
   f(x,y; b1,b2,b3) = G(b1+b2+b3)/(G(b1)*G(b2))*x^(b1-1)*y^(b2-1)*(1-x-y)^(b3-1)
 
   The value of f(x,y) is returned, with (b1,b2,b3) being defined using the mean
   and variance as:
 
     b1 = xmean*[(1-S)/Q - 1]
     b2 = ymean*[(1-S)/Q - 1]
     b2 = (1-xmean-ymean)*[(1-S)/Q - 1]

   where:
    S = xmean^2 + ymean^2 + (1-xmean-ymean)^2
    Q = xvar + yvar + xyvar
 
   NOTE: The gamma function is computed as the log of the gamma function,
   thus G(b1+b2+b3)/[G(b1)+G(b2)+G(b3)] = exp(ln[G(b1+b2+b3)] 
                                            - ln[G(b1)]-ln[G(b2)]-ln[G(b3)] )
 */
int calcDirichletPDF(int *n,double **x, double mean1, double var1, double mean2,
                     double var2, double var, double **pdf, int **pdfBound,
                     int transform);

/* function to compute a joint beta PDF given the first five moments:
   inputs:
     xmean, ymean
     xvar, yvar, xyvar   

   f(x,y; b) = b0*x^(b1-1)*y^(b2-2)*(1-x)^(b3-1)*(1-y)^(b4-1)*(1-x-y)^(b5-1)

   This form allows preservation of all moments
 */
int calcJointBetaPDF(int *n,double **x, double mean1, double var1, double mean2,
                     double var2, double var, double **pdf, int **pdfBound,
                     int transform);

/* ----- functions required for beta distribution ----- */
double betaDist(double x, double alpha, double beta);

/* 6 parameter joint-beta defined such that:
     f(x,y) = b0*x^(b1-1)*y^(b2-1)*(1-x)^(b3-1)*(1-y)^(b4-1)*(1-x-y)^(b5-1)

   NOTE: properties of this distribution have not been well defined
 */
double jointBetaDist(double x1, double x2, double *b);

/* ----- functions required for Dirichlet distribution ----- */
// computes the distribution at (x1, x2)
double dirichletDist(double x1, double x2, double b1, double b2, double b3);

/* computes the distribution in transformed space: (xi, eta)
     xi = x1,   eta = x2/(1-x1) 
 */
double dirichletDistDoran(double x1,double x2,double b1,double b2,double b3);

/* computes the distribution in transformed space: (xi, eta)
     xi = x1+x2,   eta = x2/(x1+x2) 
 */
double dirichletDistHasse(double x1,double x2,double b1,double b2,double b3);

/* function definitions and analytic Jacobians used by Newton/Broyden method */
void jointBetaRes(int nF, double *a, double *f, void *args);
void jointBetaFunc(int nF, double *a, double *f, void *args);
void jointBetaJac(int nF, double *a, double **r, void *args);

void jointBetaLevMar(double *a, double *y, double **dyda,
                     int na, void *args);

void jointBetaFuncLM(double *p, double *hx, int m, int n, void *args);
void jointBetaJacLM(double *p, double *j, int m, int n, void *args);

/* for computing the joint beta distribution using Appell series */
int calcAppellBetaPDF(int *n, double **x, double mean1,double var1,double mean2,
                      double var2, double var, double **pdf, int **pdfBound,
                      int transform);
void appellF3BetaFunc(int nF, double *a, double *f, void *args);
void appellF2BetaFunc(int nF, double *a, double *f, void *args);

#ifdef MDUC_PDF_KINSOL
#include "sdialStruct.h"

/* functions for KINSOL Newton solver */
void initKinJbeta(sdialStruct *kin, statStruct *sml, double *a);
void destroyKinJbeta(sdialStruct *kin);

int  kinJbetaRes(N_Vector u, N_Vector fval, void *f_data);
#if SUNDIALS_2_3
int kinJbetaJac(long int N, DenseMat J, N_Vector u, N_Vector fu,
                void *jac_data, N_Vector tmp1, N_Vector tmp2);
#else
int kinJbetaJac(int N, N_Vector u, N_Vector fu, DlsMat J, void *jac_data,
                N_Vector tmp1, N_Vector tmp2);
#endif /* SUNDIAL_2_3 */
#elif MDUC_PDF_LM
#include "levmar.h"
#endif /* MDUC_PDF_KINSOL */

#ifdef MDUC_JBETA_ADAPT
void jointBetaFuncODE(int nF, double *a, double *f, void *args);
void jointBetaJacODE(int nF, double *a, double **r, void *args);
void jointBetaODE_RHS_X(double x, double *y, double *dydx, void *args);
void jointBetaODE_RHS_Xgrid(double x, double *y, double *dydx, void *args);
void jointBetaODE_RHS_Y(double x, double *y, double *dydx, void *args);
void jointBetaJacODE_RHS_X(double x, double *y, double *dydx, void *args);
void jointBetaJacODE_RHS_Y(double x, double *y, double *dydx, void *args);

int jointBetaCVODE_dX(realtype t, N_Vector u, N_Vector udot, void *data);
int jointBetaCVODE_dY(realtype t, N_Vector u, N_Vector udot, void *data);
int jointBetaJacCVODE_dX(realtype t, N_Vector u, N_Vector udot, void *data);
int jointBetaJacCVODE_dY(realtype t, N_Vector u, N_Vector udot, void *data);
void jointBetaODE_RHS(double x, double *y, double *dydx, void *args);
#endif /* MDUC_JBETA_ADAPT */


#endif /* MDUC_BETA_PDF_H */
