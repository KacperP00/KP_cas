#ifndef MDUC_HYPERGEOMETRIC_H
#define MDUC_HYPERGEOMETRIC_H

/* ----- Hypergeometric series ----- */
/* Wrapper function to compute generalised hypergeometric series
     p: number of numerator coefficients
     q: number of denominator coefficients
     a: nuamerator coefficients
     b: denominator coefficients */
double hyper(int p, int q, double *a, double *b, double x);

// specific function calls to treat special cases
double hyper_2F1(double a, double b, double c, double x);
double hyper_3F2(double a, double b, double c, double e, double f, double x);

/* ----- Appell Hypergeometric Series ----- */
double appell_F1(double a, double b, double bp, double g, double x, double y,
                 int nmax, double eps);
double appell_F2(double a, double b, double bp, double g, double gp,
                 double x, double y, int nmax, double eps);
double appell_F3(double a, double ap, double b, double bp, double g,
                 double x, double y, int nmax, double eps);
double appell_F4(double a, double b, double g, double gp, double x, double y,
                 int nmax, double eps);

/* ----- summation tools ----- */
// direct summation of the hypergeometric series
double hypersum(int p, int q, double *a, double *b, double x, int nmax,
                double eps);

/* Uses epsilon algorithm with Shanks transformation to accelerate the 
   hypergeometric series.  This should only be used when the series converges
   linearly, which in the hypergeometric case is when 0 < x < 1 */
double hyperShanks(int p, int q, double *a, double *b, double x, int nmax,
                   double eps);

/* Uses Levin algorithm to accelerate the series.  This is used when x == 1
   since the convergence there is logarithmic and Shanks will not be suited. */
double hyperLevin(int p, int q, double *a, double *b, double x, int nmax,
                  double eps);

#endif /*  MDUC_HYPERGEOMETRIC_H */
