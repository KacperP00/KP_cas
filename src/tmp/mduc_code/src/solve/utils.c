#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "utils.h"
#include "display.h"
#include "accessChem.h" /* should remove dependence on this */

void getZfromGrid(int trans, double x1, double x2, double *Z1, double *Z2)
{
  if (trans == 0) {
    *Z1 = x1;
    *Z2 = x2;
  }
  else if (trans == 1) {
    *Z1 = x1;
    *Z2 = x2*(1.0-x1);
  }
  else if (trans == 2) {
    *Z1 = x1*(1.0-x2);
    *Z2 = x2*x1;
  }
}
/* -------------------------------------------------------------------------- */

void getGridfromZ(int trans, double Z1, double Z2, double *x1, double *x2)
{
  if (trans == 0) {
    *x1 = Z1;
    *x2 = Z2;
  }
  else if (trans == 1) {
    *x1 = Z1;
    *x2 = Z2/(1.0-Z1);
  }
  else if (trans == 2) {
    *x1 = Z1 + Z2;
    *x2 = Z2/(Z1+Z2);
  }
}
/* -------------------------------------------------------------------------- */

void transformZ(double Z1, double Z2, double *x1, double *x2)
{
  *x1 = Z1;
  *x2 = Z2/(1-Z1);
}
/* -------------------------------------------------------------------------- */

void retrieveZ(double x1, double x2, double *Z1, double *Z2)
{
  *Z1 = x1;
  *Z2 = x2*(1-x1);
}
/* -------------------------------------------------------------------------- */

double clip(double x, double min, double max)
{
  return fmax(fmin(x, max), min);
}
/* -------------------------------------------------------------------------- */

double percError(double val, double base)
{
  return fabs(fabs(base-val)/base)*100;
}
/* -------------------------------------------------------------------------- */

int getSpeciesIndex(int nSpecies, char **species, char *name)
{
  int i;
  for (i=0; i < nSpecies; i++)
    if (!strcmp(species[i], name))
      return i;
  
/*   printf("\nWARNING: Species not found in mechanism: %s\n", name); */
  return -1;
}
/* -------------------------------------------------------------------------- */

void convertMassToMoleFrac(int nSpecies, double *W, double *Y, double *X)
{
  int i;
  double moles;

  /* compute the total moles of each species and total moles
     NOTE: assumes mass is 1 g since absolute value is irrelevant */
  for (i=0, moles=0.0; i<nSpecies; i++)
    {
      X[i]   = Y[i]/W[i];
      moles += X[i];
    }

  for (i=0; i<nSpecies; i++) X[i] /= moles;
}
/* -------------------------------------------------------------------------- */

void convertMoleToMassFrac(int nSpecies, double *W, double *X, double *Y)
{
  int i;
  double mass;

  /* compute the total mass of each species and total mass
     NOTE: assumes 1 mole since absolute value is irrelevant */
  for (i=0, mass=0.0; i<nSpecies; i++)
    {
      Y[i]  = X[i]*W[i];
      mass += Y[i];
    }

  /* normalize */
  for (i=0; i<nSpecies; i++) Y[i] /= mass;
}
/* -------------------------------------------------------------------------- */

double calcMixtureW(int nSpecies, double *Y, double *W)
{
  int i;
  double Yi, wmix=0.0, ysum=0.0;

  for (i = 0; i < nSpecies; i++) {
    Yi = fmin(fmax(Y[i],0.0),1.0);
    wmix += Yi/W[i];
    ysum += Yi;
  }
  return ysum/wmix;
}
/* -------------------------------------------------------------------------- */

/* function to calculate the mixture specific heat given:
 number of species in mechanism, mass fraction and specific heat of each
 species */
double calcMixtureCp(int nSpecies, double *cp, double *Y)
{
  int i;
  double cpmix = 0.0;
  
  for (i = 0; i < nSpecies; i++)
    cpmix += Y[i]*cp[i];
  
  return cpmix;
}
/* -------------------------------------------------------------------------- */

int emptyString(char *a)
{
  int i;
  if (a == NULL) return 1;
  for (i=0; i<strlen(a); i++)
    if (!isspace(a[i]))
      return 0;
  return 1;
}
/* -------------------------------------------------------------------------- */

int countChars(char *str, char a)
{
  int i,count = 0;
  if (str == NULL) return 0;
  for (i=0; i<strlen(str); i++)
    if (str[i] == a)
      count++;
  return count;
}
/* -------------------------------------------------------------------------- */

void stripString(char *str)
{
  int pos, len;
  char *p;
  char *loc = (char*) calloc(strlen(str)+1,sizeof(char));

  strcpy(loc, str);
  pos = strlen(loc) - 1;
  if (pos < 0)
    loc[0] = '\0'; /* empty string */
  else {
    while (loc[pos] == ' ' || loc[pos] == '\t' || loc[pos] == '\n' &&
           pos > 0) --pos;
    loc[pos+1] = '\0';
  }
  
  p = loc;
  while (p[0] == ' ') p++;

  /* assign original string to modified local string */
  strcpy(str,p);
  free(loc);
}
/* -------------------------------------------------------------------------- */

void prefixString(char *pre, char *str)
{
  char buffer[200];
  strcpy(buffer, str);
  strcpy(str, pre);
  strcat(str, buffer);
}
/* -------------------------------------------------------------------------- */

int getIndexMax(int n, int offset, double *vec)
{
  int i, k, index = 0;
  double max = vec[0];
  for (i=1; i<n; ++i) {
    k = i*offset;
    if (vec[k] > max) {
      max = vec[k];
      index = i;
    }
  }
  return index;
}
/* -------------------------------------------------------------------------- */

double getMinValue(int n1, int n2, int index, int offset, double *vec)
{
  int j, k, row, col;
  double min = 1.0e20;
  for (k=0; k<n2; k++) {
    row = k*n1*offset;
    for (j=0; j<n1; j++) {
      col = j*offset;
      if (vec[row+col+index] < min)
        min = vec[row+col+index];
    }
  }
  return min;
}
/* -------------------------------------------------------------------------- */

double getMaxValue(int n1, int n2, int index, int offset, double *vec)
{
  int j, k, row, col;
  double max = -1.0e20;
  for (k=0; k<n2; k++) {
    row = k*n1*offset;
    for (j=0; j<n1; j++) {
      col = j*offset;
      if (vec[row+col+index] > max)
        max = vec[row+col+index];
    }
  }
  return max;
}
/* -------------------------------------------------------------------------- */

double getArrayValueStr(char *val, int n1, int n2, int n3, double ***vec)
{
  int i, j, k;
  double min = 1.0e20;
  double max = -1.0e20;
  for (k=0; k<n3; k++) {
    for (j=0; j<n2; j++) {
      for (i=0; i<n1; i++) {
        if (vec[k][j][i] < min)
          min = vec[k][j][i];         
        if (vec[k][j][i] > max)
          max = vec[k][j][i];
      }
    }
  }
  
  if (!strcmp(val, "min"))
    return min;
  else if (!strcmp(val, "max"))
    return max;
  else {
    fprintf(stderr, "Error (getArrayValueStr): invalid type, %s\n", val);
    exit(2);
  }
  return 0.0;
}
/* -------------------------------------------------------------------------- */

double linInterp(double x, double x0, double y0, double xN, double yN)
{
  return y0 + (yN-y0)/(xN-x0)*(x-x0);
}
/* -------------------------------------------------------------------------- */

double biLinInterp(double x, double y, double x0, double x1, double y0,
                   double y1, double Y00, double Y01, double Y10, double Y11)
{
  double dx, dy, f, g, q;

  f = (Y01-Y00)/(x1-x0);
  g = (Y10-Y00)/(y1-y0);
  q = (Y11-Y10-Y01+Y00)/(x1-x0)/(y1-y0);

  return Y00 + f*(x-x0) + g*(y-y0) + q*(x-x0)*(y-y0);
}
/* -------------------------------------------------------------------------- */

double lookupProfile(double xv, int n, double *x, double *y)
{
  int i;
  i = gridLocate(n, x, xv);
  return linInterp(xv, x[i], y[i], x[i+1], y[i+1]);
}
/* -------------------------------------------------------------------------- */

double calcMean(int x0, int xN, int index, int offset, double *pdf, double *Y)
{
  int j;
  double sum = 0.0;
  
  for (j=x0; j <= xN; ++j)
    sum += pdf[j]*Y[j*offset + index];
  
  return sum;
}
/* -------------------------------------------------------------------------- */

double calcMeanX(int x0, int xN, int index, int offset, double *pdf, double *Y,
                 int nSpecies, double *W)
{
  int j, k;
  double sum = 0.0;
  double WMix = 0.0;
  
  for (j = x0; j <= xN; ++j) {
    k = j*offset;
    WMix = calcMixtureW(nSpecies, &Y[k], W);
    sum += pdf[j]*Y[k+index]*WMix/W[index];
  }
  
  return sum;
}
/* -------------------------------------------------------------------------- */

double calcGridGrowth(int n, double a0, double L)
{
  int count = 0;
  double f, fp, r = 1.1;

  do {
    f = (1.0 - pow(r,n))/(1.0-r) - L/a0;
    fp = -(n)*pow(r,n-1)/(1.0-r) + (1.0-pow(r,n))/((1.0-r)*(1.0-r));
    r -= f/fp;
    count++;
    if (++count > 1000) mducErr("grid not converged\n","calcGridGrowth");
  } while (fabs(f/fp) > 1.0e-4*a0);

  return r;
}
/* -------------------------------------------------------------------------- */

void copyArray(int n, double *x, double *y)
{
  int i;
  for (i=0; i<n; i++)
    y[i] = x[i];
}
/* -------------------------------------------------------------------------- */

void convertArrayToFloat(int n, double *x, float *y)
{
  int i;
  for (i=0; i<n; i++)
    y[i] = (float)x[i];
}
/* -------------------------------------------------------------------------- */

void extractVarArray(int n1, int n2, int n3, int index, double *Y, double *var)
{
  int j, k;
  for (k=0; k<n3; k++)
    for (j=0; j<n2; j++)
      var[k*n2 + j] = Y[(k*n2 + j)*n1 + index];
}
/* -------------------------------------------------------------------------- */

void setVarArray(int n1, int n2, int n3, int index, double *Y, double *var)
{
  int j, k;
  for (k=0; k<n3; k++)
    for (j=0; j<n2; j++)
      Y[(k*n2 + j)*n1 + index] = var[k*n2 + j];
}
/* -------------------------------------------------------------------------- */

void convertArray1D_MD(int n1, int n2, int n3, double *YU, double ***YS)
{
  int i, j, k;
  for (k=0; k<n3; k++)
    for (j=0; j<n2; j++)
      for (i=0; i<n1; i++)
        YS[k][j][i] = YU[(k*n2 + j)*n1 + i];
}
/* -------------------------------------------------------------------------- */

void convertArrayMD_1D(int n1, int n2, int n3, double ***YS, double *YU)
{
  int i, j, k;
  for (k=0; k<n3; k++)
    for (j=0; j<n2; j++)
      for (i=0; i<n1; i++)
        YU[(k*n2 + j)*n1 + i] = YS[k][j][i];
}
/* -------------------------------------------------------------------------- */

void transposeRowCol(int n1, int n2, int nVars, double *Y)
{
  int i, j, k;
  int row, col, rowTemp, colTemp;
  
  /* allocate temporary storage */
  double *Ytemp = (double*) calloc(n1*n2*nVars, sizeof(double));
  
  for (k=0; k < n2; k++) {
    row     = k*n1*nVars;
    colTemp = k*nVars;
    
    for (j=0; j < n1; j++) {
      col     = j*nVars;
      rowTemp = j*n2*nVars;
      
      for (i=0; i < nVars; i++)
        Ytemp[rowTemp+colTemp+i] = Y[row+col+i];
    } /* end loop over columns */
  } /* end loop over rows */
  
  /* copy the temporary solution into the solution vector */
  for (i=0; i < n1*n2*nVars; i++)
    Y[i] = Ytemp[i];
  
  free(Ytemp);
}
/* -------------------------------------------------------------------------- */

void interpSolution1D(int n1, double *x, int nI, double *xI, int nVars,
                      double *Y, double *YI)
{
  int i, j, col, jR, jL;

  /* loop over all points on the interpolated mesh */
  for (j=0; j<nI; j++) {
    col = j*nVars;
    
    /* find the coarse mesh points either side of current fine mesh */
    jR = 0;
    do { jR++; } while ( x[jR] < xI[j] && jR < n1 );
    jL = jR - 1; /* assume for within 1 grid point */
    
    /* linearly interpolate all variables */
    for (i=0; i<nVars; i++)
      YI[col+i] = linInterp(xI[j],x[jL],Y[jL*nVars+i],x[jR],Y[jR*nVars+i]);
  }
}
/* -------------------------------------------------------------------------- */

void interpSolution2D(int n1, int n2, double *x1, double *x2, int n1I, int n2I,
                      double *x1I, double *x2I, int nVars, double *Y, double *YI)
{
  int i, j, k, row, col, jR, jL, kA, kB;
  int rowA, rowB, colL, colR;

  /* set the interpolated mesh to zero */
  for (k=0; k<n1I*n2I*nVars; k++) YI[k] = 0.0;  

  /* loop over all points on the interpolated mesh */
  for (k=0; k<n2I; k++) {
    row = k*n1I*nVars;
    
    /* find the coarse mesh points above and below current row of fine mesh */
    kB = findAdjacentGridPoint(n2, x2, x2I[k]);
    kA = kB + 1;
    
    rowA = kA*n1*nVars;
    rowB = kB*n1*nVars;
    
    for (j=0; j<n1I; j++) {
      col = j*nVars;
      
      /* find the coarse mesh points either side of current fine mesh */
      jL = findAdjacentGridPoint(n1, x1, x1I[j]);
      jR = jL + 1;
      
      colL = jL*nVars;
      colR = jR*nVars;
      
      /* bi-linearly interpolate all variables */
      for (i=0; i<nVars; i++)
        if (jL >=0 && kB >= 0)
          YI[row+col+i] = biLinInterp(x1I[j], x2I[k], x1[jL], x1[jR], x2[kB],
                                      x2[kA], Y[rowB+colL+i], Y[rowB+colR+i],
                                      Y[rowA+colL+i], Y[rowA+colR+i]);
    } /* end loop over columns */
  } /* end loop over rows */
}
/* -------------------------------------------------------------------------- */

void interpDNS_1D(int n, double *x, int norm, double xmax, double *var)
{
  int j, jL = 0, jR = 0, jMax = n-1;

  for (j=1; j<n-1; j++) {
    /* keep track of the previous non-zero point */
    if (var[j] > 0.0)
      jL = j;
    else {
      /* find the next non-zero point.  If no non-zero points found,
         will result in the maximum boundary */
      if (j > jR) {
        jR = j;
        while (var[jR] <= 0.0 && jR < jMax)
          jR++;
      }
      
      /* interpolate */
      if (jR < jMax)
        var[j] = linInterp(x[j], x[jL], var[jL], x[jR], var[jR]);
      else if (norm)
        var[j] = linInterp(x[j],x[jL],var[jL],xmax,0.0);
      else
        var[j] = 0.0;
    }
  } /* end loop over columns */
}
/* -------------------------------------------------------------------------- */

int findAdjacentGridPoint(int n1, double *x1, double x)
{
  int jR, jL;

  /* bound the value of x correctly */
  x = fmin(fmax(x, x1[0]),x1[n1-1]);

  /* find the coarse mesh points either side of current fine mesh */
  jR = 0;
  do { jR++; } while ( x1[jR] < x && jR < n1-1);
  jL = jR - 1; /* assume for within 1 grid point */

  return jL;
}
/* -------------------------------------------------------------------------- */

int gridLocate(int n1, double *x1, double x)
{
  int jL = 0, jR = n1-1, jM;
  int ascnd = (x1[n1-1] >= x1[0]);  /* whether grid is ascending or not */

  /* use bisection to find nearest grid point */
  while (jR - jL > 1) {
    /* compute a midpoint */
    jM = (jL+jR)/2;
    if (x >= x1[jM] == ascnd)
      jL = jM;
    else
      jR = jM;
  }
  
  return IMAX(0, IMIN(n1-2, jL));
}
/* -------------------------------------------------------------------------- */

int factorial(int x)
{
  int i, fact;
  
  fact = x;
  for (i=x-1; i>1; --i) fact *= i;
  return fact;
}
/* -------------------------------------------------------------------------- */

int IMIN(int a, int b)
{
  if (a < b) return a;
  return b;
}
/* -------------------------------------------------------------------------- */

int IMAX(int a, int b)
{
  if (a > b) return a;
  return b;
}
/* -------------------------------------------------------------------------- */
