#ifndef MDUC_SOLVE_UTILS_H
#define MDUC_SOLVE_UTILS_H

struct LAPINTSTRUCT {
  int    n1, n2, *mask;
  double *x1, *x2;
};
/* -------------------------------------------------------------------------- */
typedef struct LAPINTSTRUCT LapIntStruct;

/*
 Functions to calculate the first derivative
 */
double cGrad(double y_im1, double y_i, double y_ip1, double dx1, double dx2);
double fGrad(double y_i, double y_ip1, double dx);
double bGrad(double y_i, double y_im1, double dx);

/*
 Functions to calculate the second derivative
*/
double cGrad2(double y_im1, double y_i, double y_ip1, double dx1, double dx2);
double fGrad2(double y_i, double y_ip1, double y_ip2, double dx1, double dx2);
double bGrad2(double y_i, double y_im1, double y_im2, double dx1, double dx2);

/*
 Functions for solving linear equations
 */
void solveTriDiag(int n, double *a, double *b, double *c, double *u, double *r);
void triDiagNR(int n, double *a, double *b, double *c, double *u, double *r);

/* matrix inversion */
void invertMatrix(int n, double **a, double **a_inv);

/* for LU decomposition */
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double *b);

/* for QR decomposition */
void qrdcmp(double **a, int n, double *c, double *d, int *sing);
void qrsolve(double **a, int n, double *c, double *d, double *b);
void rsolve(double **a, int n, double *d, double *b);
void qrupdate(double **r, double **qt, int n, double *u, double *v);
void rotate(double **r, double **qt, int n, int i, double a, double b);

/* for SVD decomposition */
void   svdcmp(double **a, int m, int n, double *w, double **v);
void   svbksb(double **u, double *w, double **v, int m, int n, double *b,
              double *x);
double pythag(double a, double b);

/* globally convergent root finding routines */
double brent(double (*funcEval)(double, void *), double x1, double x2,
             double tol, void *args);
// Newton method with linesearch
void newt(double *x, int n, int *check, void *args, 
          void (*vecfunc)(int, double *, double *, void *),
          void (*jacFunc)(int, double *, double **, void *));

void freeNewt(int n, int *a, double *b, double *c, double *d, double *e,
              double **aa);

// Broyden's method (a secant method)
void broyden(double *x, int n, int *check, void *args,
             void (*vecFunc)(int, double *, double *, void *),
             void (*jacFunc)(int, double *, double **, void *));

void freeBroyden(int n, double *a, double *b, double *c, double *d, double *e,
                 double *f, double *g, double *h, double *i, double *j,
                 double **aa, double **bb);

// line search
void lnsrch(int n, double *xold, double fold, double *g, double *p,
            double *x, double *f, double stpmax, int *check, double *fvec,
            void *args, void (*vecFunc)(int, double *, double *, void *), 
            double (*func)(int, double *, double *, void *,
                           void(*)(int, double*, double *, void *)));

// function to numerically compute jacobian if no analytic solution available
void fdjac(int n, double *x, double *fvec, double **df, void *args,
           void (*vecFunc)(int, double *, double *, void *));

// function error to minimise
double dotF(int n, double *x, double *fvec, void *args,
            void(*vecFunc)(int, double *, double *, void *));

/* linear least squares fitting routines */
// solve using standard normal equations
void lfit(double *x, double *y, double *sig, int ndat, double *a, int *ia,
          int ma, double **covar, double *chisq,
          void (*funcs)(double, double *, int));
// solve using singular value decomposition
void svdfit(double *x, double *y, double *sig, int ndata, double *a, int ma,
            double **u, double **v, double *w, double *chisq,
            void (*funcs)(double, double *, int));
void fitPoly(double x, double *p, int np); // basis function for polynomial

// Levenberg-Marquardt Method
void mrqmin(double *x, double *y, double *sig, int ndat, double *a, int *ia,
            int ma, double **covar, double **alpha, double *chisq, void *args,
            void (*funcs)(double *, double *, double **, int, void *));
void mrqcof(double *x, double *y, double *sig, int ndat, double *a, int *ia, 
            int ma, double **alpha, double *beta, double *chisq, void *args,
            void (*funcs)(double *, double *, double **, int, void *));
void covsrt(double **covar, int ma, int *ia, int mfit);

void gaussj(double **a, int n, double **b, int m);

/* sparse linear matrix solve using biconjugate gradient method */
void Linbcg(int n, double *b, double *x, int itol, double tol, int itmax,
            int *iter, double *err, void *args, 
            void (*atimes)(int, double* , double *, int, void *), 
            void (*asolve)(int, double* , double *, int, void *));
double snrm(int n, double *sx, int itol);

// functions for Laplace interpolation
void LaplaceInterp_asolve(int n, double *b, double *x, int itrnsp, void *args);
void LaplaceInterp_atimes_1D(int n, double *u, double *r, int itrnsp,void *args);
void LaplaceInterp_atimes_2D(int n, double *u, double *r, int itrnsp,void *args);

#endif // MDUC_SOLVE_UTILS_H
