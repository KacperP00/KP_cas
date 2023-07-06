#ifndef MDUC_UDF_PREMIXED_H
#define MDUC_UDF_PREMIXED_H

#include "udf.h"

#define C_PROG(c,t)    C_UDSI(c,t,udsOffset+udf.iProg)
#define C_PROGVAR(c,t) C_UDSI(c,t,udsOffset+udf.iProgV)
#define C_PROG_G(c,t)  C_UDSI_G(c,t,udsOffset+udf.iProg)

/* initialization functions */
void UDFinitPremixed();
void UDFpremixedUpdateCFD(Domain *d);
void igniteFlameKernel();
void flameKernelDev();
void UDFpremixedUpdateProgSrc(Domain *d);

/* Laminar flame speed models */
double calcKeckLaminarFlameSpeed(double phi, double Tu, double Pu);

/* Turbulent flame speed models */
double (*calcTurbFlameSpeed)(double, cell_t, Thread*);
double calcNoTurbFlameSpeed(double sL, cell_t c, Thread *t);
double calcZimontTurbFlameSpeed(double sL, cell_t c, Thread *t);
double calcPetersTurbFlameSpeed(double sL, cell_t c, Thread *t);
double calcPetersTurbFlameSpeedIgn(double sL, cell_t c, Thread *t);

#endif /* MDUC_UDF_PREMIXED_H */
