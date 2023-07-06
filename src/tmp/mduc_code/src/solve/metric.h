#ifndef MDUC_METRIC_H
#define MDUC_METRIC_H

#include "flameStruct.h"

/* initialises the memory and coefficients for the metric operators */
void initMetric(flameStruct *flame);

/* stencil functions: 
     returns the start index of the stencil given the current point
*/
int stencil_x(int j, flameStruct *flame);
int stencil_y(int k, flameStruct *flame);

/* gradient operators:
      st - left/bottom most point in the stencel
      k  - row index
      j  - column index
      i  - species/variable index
      nV - total number of variables in array
      Y  - quantity
 */
double grad_x(int st,int k,int j,int i,int nV,double *Y,flameStruct *flame);
double grad_xx(int st,int k,int j,int i,int nV,double *Y,flameStruct *flame);
double grad_y(int st,int k,int j,int i,int nV,double *Y,flameStruct *flame);
double grad_yy(int st,int k,int j,int i,int nV,double *Y,flameStruct *flame);

/* finite difference coefficients for arbitrary high order */
void fd_coeff(int deriv, int order, double *x1, double x, double *coeff);

#endif /* MDUC_METRIC_H */
