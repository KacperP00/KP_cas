#ifndef MDUC_SML_PDF_H
#define MDUC_SML_PDF_H

#include "statsUtils.h"

/* computes a 1D statistically-most-likely pdf given a mean and variance.
   The SML distribution is such that it maximises the entropy

     H = -int{ ln(P/q)*q*dx )

   resulting in an analytical form
   
     P = exp(a0 + a1*x + a2*x^2);

 */
int calcSML_PDF(int n, double *x, double mean, double var, double *pdf,
                int *pdfBound);

/* computes a 2D statistically-most-likely pdf given first five moments by
   maximising the entropy (shown above).  The resulting form is

     P = exp(a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y);
*/
int calcSML_PDF2D(int *n, double **x, double mean1, double var1, double mean2, 
                  double var2, double var, double **pdf, int **pdfBound,
                  int transform);

/* definition of distribution */
double smlDist1D(double x, double *a);
double smlDist2D(double x1, double x2, double *a);

/* function definitions and analytic Jacobians for Newton/Broyden method */
void smlRes1D(int nF, double *a, double *f, void *args);
void smlRes2D(int nF, double *a, double *f, void *args);

void smlFunc1D(int nF, double *a, double *f, void *args);
void smlFunc2D(int nF, double *a, double *f, void *args);

void smlJac1D(int nF, double *a, double **r, void *args);
void smlJac2D(int nF, double *a, double **r, void *args);

/* wrapper functions for external Levenberg-Marquardt library */
void smlFunc1D_LM(double *p, double *hx, int m, int n, void *args);
void smlJac1D_LM(double *p, double *j, int m, int n, void *args);
void smlFunc2D_LM(double *p, double *hx, int m, int n, void *args);
void smlJac2D_LM(double *p, double *j, int m, int n, void *args);

#ifdef MDUC_PDF_KINSOL
#include "sdialStruct.h"

/* functions for KINSOL Newton solver */
void initKinSML(sdialStruct *kin, statStruct *sml, double *a);
void destroyKinSML(sdialStruct *kin);

int  kinSML_Res(N_Vector u, N_Vector fval, void *f_data);
#if SUNDIALS_2_3
int kinSML_Jac(long int N, DenseMat J, N_Vector u, N_Vector fu,
               void *jac_data, N_Vector tmp1, N_Vector tmp2);
#else
int kinSML_Jac(int N, N_Vector u, N_Vector fu, DlsMat J, void *jac_data,
               N_Vector tmp1, N_Vector tmp2);
#endif /* SUNDIAL_2_3 */
#elif MDUC_PDF_LM
#include "levmar.h"
#endif /* MDUC_PDF_KINSOL */

#endif /* MDUC_SML_PDF_H */
