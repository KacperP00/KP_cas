#ifndef MDUC_UDF_WIEBE_H
#define MDUC_UDF_WIEBE_H

#include "mducUDF.h"

typedef struct MDUC_WIEBE_STRUCT {
  /* I/O */
  FILE *oFile;        /* pointer to timing file */
  char outFile[256];  /* file to write Wiebe source terms */

  /* function parameters */
  double a;   /* leading coefficient of Wiebe function */
  double m;   /* exponent of Wiebe function */

  /* duration parameters */
  double alpha; /* duration constant */
  double beta;  /* duration linear coefficient */
  double SOC;   /* start of combustion */
  double dt;    /* combustion duration */

  /* parameters for the overall hydrocarbon reaction (C_a H_b) */
  double nC, nH;
  double src;   /* common source term */

  double mf_n;

  /* coefficients for the species source term equations */
  double aFuel, aO2, aH2O, aCO2;

  /* keep track of integrated src terms */
  double sFuel, sO2, sH2O, sCO2, sHTOT;
  double bFuel, bO2, bH2O, bCO2, bHTOT;
};
typedef struct MDUC_WIEBE_STRUCT wiebeStruct;

/* called during initialization */
void UDFinitWiebe();
void UDFinitWiebeFile(wiebeStruct* w);

/* function to be called during ADJUST */
void WiebeCalcSrc(Domain *d);
void writeWiebe(Domain *d);

/* functions for computing the Wiebe function */
double UDFcalcWiebe(double a, double m, double t0, double dt, double t);
double UDFcalcWiebe_dt(double a, double m, double t0, double dt, double t);

#endif /* MDUC_UDF_WIEBE_H */
