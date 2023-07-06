#ifndef MDUC_SERIES_H
#define MDUC_SERIES_H

/* Computes the next term in a series for the epsilon algorithm using a 
   Shanks transformation */
double ShanksNext(double sum, int *n, int *ncv, double eps, double *e,
                  double *lastval, double *lasteps, int *cnvgd);

/* Uses Levin transformation for series acceleration (non-linear) */
double LevinNext(double sum, double omega, double beta, int *n, int *ncv,
                 double eps, double *num, double *den, double *lastval,
                 double *lasteps, int *cnvgd);

#endif /* MDUC_SERIES_H */
