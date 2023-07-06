#ifndef MDUC_UDF_RIF_H
#define MDUC_UDF_RIF_H

#include "mducUDF.h"

/* check status of injection to determine which solver to use */
int  checkInjection(Domain *d);

/* called during main initialization */
void UDFinitRIF();

/* main decision tree function to solve the RIF */
void UDFsolveRIF(Domain *d, int injection, double PAvg);

/* ===== 1D functions ===== */
/* main function for solving 1D representative flamelet */
void UDF_rif_1D(Domain *d, int inert, double PAvg, void *mduc);

/* function to update the CFD domain from a 1D flamelet solution */
void updateCFD_1D(Domain *d, int nYFav, int *indexYFav, int n1, double *x1,
                  double **Y, double *Temp, double *pdfGlobal);

/* ===== 2D functions ===== */
/* main function for solving 2D representative flamelet */
void UDF_rif_2D(Domain *d, int inert, double PAvg, void *mduc);

/* function to update the CFD domain from a 2D flamelet solution */
void updateCFD_2D(Domain *d, int nYFav, int *indexYFav, int *n, double **x,
                  double ***Y, double **Temp, double **pdfGlobal);

void adjust_T_cond_dTdZ(Domain *d, int n1, double *x1);
void computeModH_src(Domain *d);

/* miscellaneous functions */
void checkRIF_coord(char *name, int *iMean, int *iVar);

void UDFupdateProgSRC(Domain *d);
void UDFupdate(Domain *d);
#endif /* MDUC_UDF_RIF_H */
