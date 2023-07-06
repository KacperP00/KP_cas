#ifndef MDUC_ACCESS_H
#define MDUC_ACCESS_H

#ifdef ACCESS_MPI
#include "mpiStruct.h"
#endif /* ACCESS_MPI */

/* ----- functions required for initialisation and solution ----- */

void *MDUCinit(char *inputFileName);
void *MDUCinitChem(char *inputFileName);
void MDUCdestroy(void *mem);

void MDUCinitFlow(double t, double P, double phi1min, double phi1max,
                  double phi2min, double phi2max, void *mem);
void MDUCinit2Dfrom_BC_and_1D(int nGrid2,int nFuel2,char **fuelName,double *Y,
                              double Tfuel2,void *mem);
void MDUCinit2Dfrom_1D_and_1D(void *mem2, void *mem1);

int MDUCisParallel(void);
#ifdef ACCESS_MPI
void MDUCtempPrintProc(int printproc);
void MDUCsetPrintProc(MPI_Comm GLOBAL_COMM);
void MDUCsetParallel(MPI_Comm GLOBAL_COMM, void *mem);
void MDUCsetParallelmRIF(MPI_Comm GLOBAL_COMM, void *mem, int nLocal, int hostLocal);
void initParallel1D(void *mem);
void initParallel2D(void *mem);
void initParallel2D_ADI(void *mem);
#else
void MDUCsetPrintProc(void *dum);
#endif /* ACCESS_MPI */

void MDUCupdateTime(double t, double dt, int timeStep, void *mem);
void MDUCupdateTimestep(int timeStep, void *mem);

/* function to set conditional mean enthalpy source term */
void MDUCsetAvgSdot(int m, int n, double *Sdot, void *mem);
/* function to set average H for initialization */
void MDUCsetAvgH(int n, double *AvgH, void *mem);
/* Set grid bounds */
void MDUCsetGridBound(double x1min, double x1max, double x2min,
                      double x2max, void *mem);
/* set the scalar dissipation rates */
void MDUCsetScalarDiss(int n1, int n2, double **chi1, double **chi12,
                       double **chi2, void *mem);
void MDUCsetPDF(int n1, int n2, double **pdf, void *mem);
void MDUCsetPressure(int n, double P, void *mem);
void MDUCsetTempBounds(double T1, double T2, double T3, void *mem);

/* update the solution based on new boundary conditions without chemistry */
void MDUCupdateInert1D(double dt, double P, double phimin, double phimax,
                       double Tox, double Tfuel, double *chi, double **Y,
                       double *Temp, void *mem);
void MDUCupdateInert2D(double dt, double P, double phi1min, double phi1max,
                       double phi2min, double phi2max,double Tox,double Tfuel1,
                       double Tfuel2, double **chi1, double **chi12,
                       double **chi2, double ***Y, double **Temp, void *mem);

/* update the temperature based on the dP/dt term (if not solved) */
void MDUCupdateTempDPDT(double P, double dPdt, double **T, void *mem);

void MDUCflamelet1D(double dt,double P,double phimin,double phimax,double *chi,
                    int inert,double **Y,double *Temp,void *mem);
void MDUCflamelet2D(double dt,double P,double phi1min,double phi1max,
                    double phi2min,double phi2max,double **chi1,double **chi12,
                    double **chi2,int inert,double ***Y,double **Temp,
                    void *mem);
void MDUCsolve1D(double dt, double P, double phimin, double phimax,double *chi,
                 double **Y, double *Temp, void *mem);
void MDUCsolve2D(double dt, double P, double phi1min, double phi1max,
                 double phi2min, double phi2max, double **chi1, double **chi12,
                 double **chi2, double ***Y, double **Temp, void *mem);

double MDUCcalcResidual(void *mem);

/* ----- functions for retrieving info about the combustion solver ----- */
/* return the grid properties, either of solution space or physical */
int    MDUCgetGridSize(int dir, void *mem);
int    MDUCgetHdim(void *mem);
int    MDUCgetNdim(void *mem);
int    MDUCgetGrid(double **grid, void *mem);
int  * MDUCgetGrid2D(double ***grid, void *mem);
void   MDUCgetZGrid(double ***Z1grid, double ***Z2grid, void *mem);
void   MDUCgetGridBound(double *x1min, double *x1max, double *x2min,
                        double *x2max, void *mem);

/* retrieve info about the chemical mechanism */
int    MDUCgetNspecies(void *mem);
int    MDUCgetNspecs(void *mem);
int    MDUCgetSpeciesIndex(char *species, void *mem);
void   MDUCcomputeThermoData(double *h, double *cp, double T, void *mem);
double MDUCgetMolWeight(char *species, void *mem);
void   MDUCgetMolWeightAll(double *W, void *mem);
void   MDUCgetSpeciesNames(char ***names, int buf_len, void *mem);

/* chemtable functions */
void  *MDUCloadChemTable(char *chemTabFile, int nDim, char *type);
double MDUClookupChemTable(char *name, double *x, double *vars, void *tab);
double MDUClookupChemTableIndex(int *index, double *x, double *vars, void *tab);
double MDUCchemTableVarMin(char *name, void *tab);
double MDUCchemTableVarMax(char *name, void *tab);
double MDUCchemTableProgMax(double *x, void *tab);
double MDUCchemTableProgMin(double *x, void *tab);
int    MDUCgetChemTabNdim(void *tab);
int    MDUCgetChemTabNvar(void *tab);
int    MDUCgetChemTabNspecies(void *tab);
double MDUCgetChemTabPref(void *tab);
double MDUCgetChemTabTref(void *tab);
void   MDUCgetChemTabVarNames(char **names, void *tab);
void   MDUCgetChemTabSpecies(char **species, void *tab);
double MDUCchemTableProgScaleP(double P, double *x, void *tab);

/* retrieve solution */
double *MDUCgetZrefV(void *mem);
double *MDUCgetChiRefV(void *mem);
double MDUCgetZref(int dir, void *mem);
double MDUCgetChiZ1max(void *mem);
void   MDUCgetZmoments(int dir, double *mean, double *var, void *mem);
double MDUCgetChiRef(int dir, void *mem);
double MDUCgetStrainRate(int dir, void *mem);
void   MDUCgetPointY(double x1, double x2, double *Y, void *mem);
double MDUCgetPointT(double x1, double x2, void *mem);
double MDUCgetPointH(double x1, double x2, void *mem);
void   MDUCgetSolution(double ***Y, double **Temp, void *mem);
void   MDUCgetSootCoeff(double **sqrtT, double **T_mu, double **muSqrtMw_rhoSqrtT,
		      double **sgCoeff, double **oxCoeff, double **o2Coeff,
		      double **dimProdRate, void *mem);
void   MDUCgetSootPAH(double **yPAH, double **srcPAHpos, double **srcPAHneg,
		      double **rhoDot, void *mem);
void   MDUCgetVar(char *name, double **var, void *mem);
void   MDUCcalcHRrate(void *mem);
double MDUCgetHRrate(double volume, void *mem);
double MDUCgetHeatLoss(void *mem);

/* retrieve info about the flow parameters */
double MDUCgetPressure(int n, void *mem);
double MDUCgetHbounds(int n, void *mem);
double MDUCgetRhobounds(int n, void *mem);
double MDUCgetTime(int n, void *mem);
double MDUCgetDPDT(void *mem);

void MDUCstoreCvodePre(void *mem);
void MDUCsetSolution(double t, double P, double ***Y, double **Temp,void *mem);
void MDUCresetODE(void *mem);

/* updates the oxidiser boundary temperature based on a given total enthalpy */
double MDUCupdateOxTemp(double H, double guess, void *mem);

/* calculate the total enthalpy at point j (mesh point) using
   flamelet solution */
double MDUCcalcHtot(int iRow, int iCol, double T, void *mem);

/* iterates for temperature from enthalpy and species
   ***** NOTE: right now only for ox boundary, should be for any point */
double MDUCcalcTfromH(int nSolve, double H, double *Y, double guess, void *mem);

/* ----- functions to display, read or write data ----- */
void MDUCprintStream(void *mem); /* lists fuel and oxidiser b.c. */
void MDUClistSpecies(void *mem); /* lists the species in mechanism */

/* set indices for species tracked in CFD field */
void MDUCsetSpeciesIndex(int nY, char **species, int **index, void *mem);

void MDUCsetDataDir(char *directory, void *mem);
void MDUCsetFilePrefix(char *prefix, void *mem);
/* enables setting of the monitor directly from calling source code */
void MDUCsetMonitor(int nY, int nZ, char **monitorSpecies, void *mem);

/* write a standard flamelet file */
void MDUCwriteFlamelet(char *fileName, void *mem);
/* write a reaction pathway file */
void MDUCwriteRates(char *fileName, void *mem);

/* add a variable to a flamelet file with current name */
void MDUCwriteVariableToFlamelet(char *varName, double **var, void *mem);

/* temporary functions for enabling output of PDF */
double MDUCreadFlowData(double **chi1, double **chi12, double **chi2,
                        double **pdf, double *sdot, double timeval, void *mem);
void MDUCwriteFlowData(double **chi1, double **chi2, double **pdf1,
                       double *sdot, void *mem);
void MDUCwriteFlowDataMem(char *fileName, void *mem);
void MDUCtransformVar(int s, int trans, int n1, int n2, double *x1, double *x2,
                      double **var);
void MDUCreadVariable(char *fileName,char *varName,int n1,int n2,double **var);

void MDUCtransformVar(int s, int trans, int n1, int n2, double *x1, double *x2,
                      double **var);

/* ----- functions for evaluating statistics ----- */
int MDUCcalcPDF(double mean1, double var1, double mean2, double var2,double var,
                double **pdf, int **pdfBound, void *mem);

int MDUCcalcPDF_1D(int method, int n, double *x, double mean, double var,
                   double *pdf, int *pdfBound);

int MDUCcalcPDF_2D(int method, int *n, double **x, double mean1, double var1,
                   double mean2, double var2, double var, double **pdf,
                   int **pdfBound, int transform);

int MDUCcalcMeanY(double mean1, double var1, double mean2, double var2,
                  double var, double ***Y, double *Ymean, double **pdf,
                  int **pdfBound, void *mem);

int MDUCcalcMeanY_1D(int method, int n, double *x, int nSpecies, double mean,
                     double var, double **Y, double *Ymean, double *pdf,
                     int *pdfBound);

void MDUCsetSootMom(int n1, int nSootMom, double **sootMom, void *mem);
void MDUCsetViscSrc(int n1, double *viscSrc, void *mem);
void MDUCsetGlobalPDF(int n1, double *gPDF, void *mem);

int MDUCcalcMeanSootSrc(double mean1, double var1, double mean2, double var2,
                        double var, double *meanSootSrc, void *mem);

void MDUCgetSootSolution(double **Ssrc, void *mem);
void MDUCgetRadiation1D(double *S, void *mem);

void MDUCsootHMOM(int nMom,double T_SI, double Visc_SI, double Rho_SI, double W_SI,
                  double *c_spec_SI, double *mom_SI, double *src_SI, double *rxnrates_SI,
                  double *rhodot_SI, double *srcP_SI, int *srcActive);

int MDUCcalcMeanY_2D(int method, int *n, double **x, int nSpecies,
                     double mean1, double var1, double mean2, double var2,
                     double var, double ***Y, double *Ymean, double **pdf,
                     int **pdfBound, int transform);

int MDUCcalcMomentsPDF_2D(int *n, double **x, double **X, double *meanPDF,
                          double *varPDF, double **pdf, int **pdfBound,
                          int pdfFlag, int transform);

double MDUCcalcMeanYi(int n, double *x, int iY, double mean, double var,
                      double **Y, double *pdf, int *pdfBound, int pdfFlag);

double MDUCcalcMeanVar(char *varName, double mean1, double var1, double mean2,
                       double var2, double var0, double **pdf, int **pdfBound,
                       int pdfFlag, void *mem);

double MDUCcalcMeanDensity(double P, double **TIn, double **pdfIn,
                           int **pdfBound, int pdfFlag, void *mem);

double MDUCcalcLocalMeanVar(char *varName, double *mean, double *var,
                            double **pdf, int **pdfBound,int pdfFlag,void *mem);

double MDUCcalcMeanGamma(double **pdfIn, int **pdfBound, int pdfFlag,
                         double *Rmean, void *mem);

void MDUCupdateY_pre(void *mem);
void MDUCupdateHRrate(double **HR, double **H, void *mem);

/* compute Favre averaged values of all species based on supplied PDF and write
   to file */
void MDUCwriteScalars(double curTime, double *mean,double *var, double **pdf,
                      int **pdfBound, int pdfFlag, void *mem);

void MDUCcalcAvgSdot(int nCFD1, int nCFD2, int nCFD3, double ***Zmean,
                     double ***Zvar, double ***Hmean, double ***Hvar,
                     double ***rhoMean, void *mem);

/* ----- functions for scalar dissipation rate ----- */
/* for interpolating raw data from DNS */
void MDUCinterpChiDNS(double **chi1, double **chi12, double **chi2, void *mem);
void MDUCfitChi1D(int n1, double xmin, double xmax, double *x1, double *chi1);
void MDUClaplaceChi1D(int n1,double xmin,double xmax,double *x1,double *chi1);
void MDUCinterpChiLap(double **chi1, double **chi12, double **chi2, void *mem);

/* compute the functional form */
double MDUCcalcChiLogForm(double Z, double Zref, double Zmax);
double MDUCcalcChiErfcForm(double Z, double Zref, double Zmax);
double MDUCcalcChiQuadForm(double Z, double Zref, double Zmax);
int  MDUCgetNsootMom(void *mem);
void MDUCsetRhoDot(int n1,double *rhoDot, void *mem);
void MDUCgetSootSrc(double **Ssrc, void *mem);

/* ----- functions for operating on the flow field ----- */
#ifdef ACCESS_CFD
/* check whether an injection has occurred by looking
   for mean mixture fraction */
int MDUCcheckInjection(int nCFD1, int nCFD2, int nCFD3, double threshold,
                       double ***Z1mean, double ***Z2mean);

void MDUCupdateCFD_1D(int nCFD1, int nCFD2, int nCFD3, int nYFav,
                      int *indexYFav, double ***C_Y, double ***C_T,
                      double ***C_H, double ***C_PHIMEAN, double ***C_PHIVAR,
                      int dim, int n1, double *x, double **Y, double *T);
/* Errorfunction */
double MDUCerf(double x);
#endif /* ACCESS_CFD */

#endif /* MDUC_ACCESS_H */
