#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "accessMDUC.h"
#include "accessChem.h"
#include "utils.h"

#ifdef ACCESS_MPI
#include "mpiStruct.h"
#endif /* ACCESS_MPI */

/* Local functions for converting between c and fortran arrays
   for .f to .c, specify n1, n2, n3 in the reverse order as allocated in fortran
     e.g. allocate( var(n1,n2,n3) ) -> *ArrayFtoC(n3,n2,n1,...)
   Special functions exist for converting a species/scalar array, as the order
   is not a direct translation of the fortran order.
   If the array is 2D, set n3 to 1
 */
void MDUCintArrayFtoC(int n1, int n2, int n3, int *fval, int ***cval);
void MDUCintArrayCtoF(int n1, int n2, int n3, int ***cval, int *fval);
void MDUCrealArrayFtoC(int n1, int n2, int n3, double *fval, double ***cval);
void MDUCrealArrayCtoF(int n1, int n2, int n3, double ***cval, double *fval);
void MDUCspecArrayFtoC(int n1, int n2, int n3, double *fval, double ***cval);
void MDUCspecArrayCtoF(int n1, int n2, int n3, double ***cval, double *fval);

/* ========================================================================== */
/*                      initialisation/memory functions                       */
/* ========================================================================== */

void mductempprintproc_(int printproc) { MDUCtempPrintProc(printproc); }
void mductempprintproc__(int printproc) { MDUCtempPrintProc(printproc); }

void mducinit_(char *fileName, int **mem) { *mem = MDUCinit(fileName); }
void mducinit__(char *fileName, int **mem) { mducinit_(fileName, mem); }

void mducinitflow_(double *t, double *P, double *phi1min, double *phi1max,
                   double *phi2min, double *phi2max, int **mem)
{ MDUCinitFlow(*t, *P, *phi1min, *phi1max, *phi2min, *phi2max, *mem); }
void mducinitflow__(double *t, double *P, double *phi1min, double *phi1max,
                    double *phi2min, double *phi2max, int **mem)
{ mducinitflow_(t,P,phi1min,phi1max,phi2min,phi2max,mem); }

void mducsetgridbound_(double *x1min, double *x1max, double *x2min,
                      double *x2max, int **mem)
{ MDUCsetGridBound(*x1min, *x1max, *x2min, *x2max, *mem);}
void mducsetgridbound__(double *x1min, double *x1max, double *x2min,
                      double *x2max, int **mem)
{ mducsetgridbound_(x1min, x1max, x2min, x2max, mem);}

void mducinit2dfrom_bc_and_1d_(int *nGrid2,int *nFuel2,int *len,char *fuelName,
                               double *Y,double *Tfuel2,int **mem)
{
  int i, j;
  char **nameL = (char**) calloc(*nFuel2, sizeof(char*));
  for (i=0; i<*nFuel2; i++)
    {
      nameL[i] = (char*) calloc(*len, sizeof(char));
      for (j=0; j<*len; j++)
        nameL[i][j] = fuelName[i+j*(*len)];
    }

  MDUCinit2Dfrom_BC_and_1D(*nGrid2, *nFuel2, nameL, Y, *Tfuel2, *mem);
  for (i=0; i<*nFuel2; i++) free(nameL[i]);
  free(nameL);
}
void mducinit2dfrom_bc_and_1d__(int *nGrid2,int *nFuel2,int *len,char *fuelName,
                                double *Y,double *Tfuel2,int **mem)
{ mducinit2dfrom_bc_and_1d_(nGrid2,nFuel2,len,fuelName,Y,Tfuel2,mem); }

void mducinit2dfrom_1d_and_1d_(int **mem2, int **mem1)
{ MDUCinit2Dfrom_1D_and_1D(*mem2, *mem1); }
void mducinit2dfrom_1d_and_1d__(int **mem2, int **mem1)
{ mducinit2dfrom_1d_and_1d_(mem2, mem1); }

int mducisparallel_() {return MDUCisParallel();}
int mducisparallel__() {return mducisparallel_();}

#ifdef ACCESS_MPI
void mducsetparallel_(MPI_Fint *f_handle, int **mem) {
  MPI_Comm GLOBAL_COMM = MPI_Comm_f2c(*f_handle);
  MDUCsetParallel(GLOBAL_COMM, *mem);
}
void mducsetparallel__(MPI_Fint *f_handle, int **mem)
{ mducsetparallel_(f_handle, mem); }

void mducsetparallelmrif_(MPI_Fint *f_handle, int **mem, int *nLocal, int *hostLocal) {
  MPI_Comm GLOBAL_COMM = MPI_Comm_f2c(*f_handle);
  MDUCsetParallelmRIF(GLOBAL_COMM, *mem, *nLocal, *hostLocal);
}
void mducsetparallelmrif__(MPI_Fint *f_handle, int **mem, int *nLocal, int *hostLocal)
{ mducsetparallelmrif_(f_handle, mem, nLocal, hostLocal); }
#else
void mducsetparallel_(int *dummy, int **mem) {};
void mducsetparallel__(int *dummy, int **mem) {};
void mducsetparallelmrif_(int *dummy, int **mem, int *dummy2,int *dummy3) {};
void mducsetparallelmrif__(int *dummy, int **mem, int *dummy2, int *dummy3) {};
#endif /* ACCESS_MPI */

void mducsetdatadir_(char *directory,int **mem)
{ MDUCsetDataDir(directory, *mem);}
void mducsetdatadir__(char *directory,int **mem)
{ mducsetdatadir_(directory,mem);}

void mducdestroy_(int **mem) { MDUCdestroy(*mem); }
void mducdestroy__(int **mem) { mducdestroy_(mem); }

/* ========================================================================== */
/*                             data retrieval                                 */
/* ========================================================================== */
int mducgetgridsize_(int *dir,int **mem)  { return MDUCgetGridSize(*dir,*mem);}
int mducgetgridsize__(int *dir,int **mem) { return mducgetgridsize_(dir,mem); }

int mducgethdim_(int **mem)  { return MDUCgetHdim(*mem);}
int mducgethdim__(int **mem) { return mducgethdim_(mem); }

int mducgetndim_(int **mem)  { return MDUCgetNdim(*mem); }
int mducgetndim__(int **mem) { return mducgetndim_(mem); }

int mducgetspeciesindex_(char species[25],int **mem)
{ 
  int index = MDUCgetSpeciesIndex(species,*mem);
  if (index >= 0) index++; // convert to Fortran index
  return index;
}
int mducgetspeciesindex__(char species[25], int **mem)
{ return mducgetspeciesindex_(species,mem); }

void mducgetspeciesnames_(char *species, int *len, int **mem){
  // note that c end character needs to be removed in fortran routine afterwards
  char **names;
  int i;
  int imax = MDUCgetNspecies(*mem);

  MDUCgetSpeciesNames(&names,*len, *mem);

  for (i=0;i<imax;i++){
    strncpy(&(species)[i*(*len)],names[i],(*len));
  }
  free(names);
}
void mducgetspeciesnames__(char *species, int *len, int **mem) 
{mducgetspeciesnames_(species, len, mem);}

void mducgetgrid_(double *phi, int **mem) {
  int i, n;
  double *x;
  n = MDUCgetGrid(&x, *mem);
  for (i=0; i<n; i++) phi[i] = x[i];
  free(x);
}
void mducgetgrid__(double *phi, int **mem) { mducgetgrid_(phi, mem); }

void mducgetgrid2d_(double *phi1, double *phi2, int **mem) {
  int i, *n;
  double **x;
  n = MDUCgetGrid2D(&x, *mem);
  for (i=0; i<n[0]; i++) phi1[i] = x[0][i];
  for (i=0; i<n[1]; i++) phi2[i] = x[1][i];
  free(x[0]); free(x[1]); free(x); free(n);
}
void mducgetgrid2d__(double *phi1, double *phi2, int **mem)
{ mducgetgrid2d_(phi1, phi2, mem); }

double mducgetpressure_(int *n,int **mem)  { return MDUCgetPressure(*n,*mem); }
double mducgetpressure__(int *n,int **mem) { return mducgetpressure_(n,mem); }

double mducgettime_(int *n,int **mem)  { return MDUCgetTime(*n,*mem); }
double mducgettime__(int *n,int **mem) { return mducgettime_(n,mem); }

double mducgetzref_(int *dir,int **mem)  { return MDUCgetZref(*dir,*mem); }
double mducgetzref__(int *dir,int **mem) { return mducgetzref_(dir,mem); }

double mducgetchiref_(int *dir,int **mem)  { return MDUCgetChiRef(*dir,*mem); }
double mducgetchiref__(int *dir,int **mem) { return mducgetchiref_(dir,mem); }

double mducgetpointh_(double *x1, double *x2, int **mem)
{ return MDUCgetPointH(*x1,*x2,*mem); }
double mducgetpointh__(double *x1, double *x2, int **mem)
{ return mducgetpointh_(x1,x2,mem); }

double mducgetpointt_(double *x1, double *x2, int **mem)
{ return MDUCgetPointT(*x1,*x2,*mem); }
double mducgetpointt__(double *x1, double *x2, int **mem)
{ return mducgetpointt_(x1,x2,mem); }

void mducgetpointy_(double *x1, double *x2, double *Y, int **mem)
{ MDUCgetPointY(*x1,*x2,Y,*mem); }
void mducgetpointy__(double *x1, double *x2, double *Y, int **mem)
{ mducgetpointy_(x1,x2,Y,mem); }

void mducgetvar_(char *name, double *var, int **mem)
{
  int i, j, k;
  int n1       = MDUCgetGridSize(1, *mem);
  int n2       = MDUCgetGridSize(2, *mem);
  double  **varL = (double**)  calloc(n2, sizeof(double*));
  for (k=0; k<n2; k++)
    {
     varL[k] = (double*)  calloc(n1, sizeof(double*));
    }
  MDUCgetVar(name, varL, *mem);
  // convert arrays to fortran and free local memory
  convertArrayMD_1D(n1, n2, 1, &varL, var);
  for (k=0; k<n2; k++)
    free(varL[k]);
  free(varL);
}
void mducgetvar__(char *name, double *var, int **mem)
{ mducgetvar_(name, var, mem); }

void mducgetsootcoeff1d_(double **mem, int *use_pah, double *sqrtT, double *T_mu, double *muSqrtMw_rhoSqrtT,
			 double *sgCoeff, double *oxCoeff, double *o2Coeff, double *dimProdRate,
			 double *yPAH, double *srcPAHpos, double *srcPAHneg, double *rhoDot)
{
  int i, j;
  int n1 = MDUCgetGridSize(1, *mem);

  MDUCgetSootCoeff(&sqrtT, &T_mu, &muSqrtMw_rhoSqrtT, &sgCoeff, &oxCoeff, &o2Coeff,&dimProdRate,*mem);

  if (*use_pah == 1) 
     MDUCgetSootPAH(&yPAH, &srcPAHpos, &srcPAHneg,&rhoDot, *mem);

}
void mducgetsootcoeff1d__(double **mem, int *use_pah, double *sqrtT, double *T_mu, double *muSqrtMw_rhoSqrtT,
			 double *sgCoeff, double *oxCoeff, double *o2Coeff, double *dimProdRate,
			 double *yPah, double *srcPAHpos, double *srcPAHneg, double *rhoDot)
{ mducgetsootcoeff1d_(mem,use_pah,sqrtT,T_mu,muSqrtMw_rhoSqrtT,sgCoeff,oxCoeff,o2Coeff,dimProdRate,yPah,srcPAHpos,srcPAHneg,rhoDot); }

void mducgetsootsolution1d_(int *nMom, double *fS, int **mem)
{
  int i, j;
  int n1       = MDUCgetGridSize(1, *mem);
  int nSoot    = MDUCgetNsootMom(*mem);

  if (nSoot != *nMom)
     fprintf(stderr,"Number of soot moments not consisent\n");

  // allocate local memory for species
  double **SootLocal = (double**) calloc(n1, sizeof(double*));
  for (j=0; j<n1; j++) SootLocal[j] = (double*) calloc(nSoot, sizeof(double));

  MDUCgetSootSolution(SootLocal, *mem);

  // convert arrays to fortran and free local memory
  MDUCspecArrayCtoF(nSoot, n1, 1, &SootLocal, fS);
  for (j=0; j<n1; j++) free(SootLocal[j]); free(SootLocal);
}

void mducgetsootsolution1d__(int *nMom, double *fS, int **mem)
{ mducgetsootsolution1d_(nMom,fS,mem); }
void mducgetradiation1d_(double *fS, int **mem)
{
  //int n1       = MDUCgetGridSize(1, *mem);
  MDUCgetRadiation1D(fS, *mem);
}

void mducgetradiation1d__(double *fS, int **mem)
{ mducgetradiation1d_(fS,mem); }
void mducgetsolution1d_(double *Y, double *Temp, int **mem)
{
  int i, j;
  int n1       = MDUCgetGridSize(1, *mem);
  int nSpecies = MDUCgetNspecs(*mem);

  // allocate local memory for species
  double **YLocal = (double**) calloc(n1, sizeof(double*));
  for (j=0; j<n1; j++) YLocal[j] = (double*) calloc(nSpecies, sizeof(double));

  MDUCgetSolution(&YLocal, &Temp, *mem);

  // convert arrays to fortran and free local memory
  MDUCspecArrayCtoF(nSpecies, n1, 1, &YLocal, Y);
  for (j=0; j<n1; j++) free(YLocal[j]); free(YLocal);
}
void mducgetsolution1d__(double *Y, double *Temp, int **mem)
{ mducgetsolution1d_(Y,Temp,mem); }

void mducgetsolution2d_(double *Y, double *Temp, int **mem)
{
  int i, j, k;
  int n1       = MDUCgetGridSize(1, *mem);
  int n2       = MDUCgetGridSize(2, *mem);
  int nSpecies = MDUCgetNspecs(*mem);

  // allocate local memory for species
  double ***YLocal = (double***) calloc(n2, sizeof(double**));
  double  **TLocal = (double**)  calloc(n2, sizeof(double*));
  for (k=0; k<n2; k++)
    {
      YLocal[k] = (double**) calloc(n1, sizeof(double*));
      TLocal[k] = (double*)  calloc(n1, sizeof(double*));
      for (j=0; j<n1; j++)
        YLocal[k][j] = (double*) calloc(nSpecies, sizeof(double));
    }

  // retrieve solution from flamelet memory
  MDUCgetSolution(YLocal, TLocal, *mem);

  // convert arrays to fortran and free local memory
  MDUCspecArrayCtoF(nSpecies, n1, n2, YLocal, Y);
  convertArrayMD_1D(n1, n2, 1, &TLocal, Temp);
  for (k=0; k<n2; k++)
    {
      for (j=0; j<n1; j++) free(YLocal[k][j]);
      free(YLocal[k]); free(TLocal[k]);
    }
  free(YLocal); free(TLocal);
}
void mducgetsolution2d__(double *Y, double *Temp, int **mem)
{ mducgetsolution2d_(Y,Temp,mem); }

/* ========================================================================== */
/*                              data modification                             */
/* ========================================================================== */
void mducupdatetime_(double *t, double *dt, int *iter, int **mem)
{ MDUCupdateTime(*t,*dt,*iter,*mem); }
void mducupdatetime__(double *t, double *dt, int *iter, int **mem)
{ mducupdatetime_(t,dt,iter,mem); }

void mducupdatetimestep_(int *iter, int **mem)
{ MDUCupdateTimestep(*iter,*mem); }
void mducupdatetimestep__(int *iter, int **mem)
{ mducupdatetimestep_(iter,mem); }

void mducsetsolution1d_(double *t, double *P, double *Y, double *T, int **mem)
{
  int i, j;
  int n1       = MDUCgetGridSize(1, *mem);
  int nSpecies = MDUCgetNspecs(*mem);

  // allocate local memory for species and convert to C order
  double **YLocal = (double**) calloc(n1, sizeof(double*));
  for (j=0; j<n1; j++) YLocal[j] = (double*) calloc(nSpecies, sizeof(double));
  MDUCspecArrayFtoC(nSpecies, n1, 1, Y, &YLocal);

  MDUCsetSolution(*t, *P, &YLocal, &T, *mem);
  for (j=0; j<n1; j++) free(YLocal[j]); free(YLocal);
}
void mducsetsolution1d__(double *t, double *P, double *Y, double *T, int **mem)
{ mducsetsolution1d_(t,P,Y,T,mem); }

void mducsettempbounds_(double *T1,double *T2,double *T3,int **mem)
{ MDUCsetTempBounds(*T1,*T2,*T3,*mem); }
void mducsettempbounds__(double *T1,double *T2,double *T3,int **mem)
{ mducsettempbounds_(T1,T2,T3,mem); }

void mducresetode_(int **mem)  { MDUCresetODE(*mem); };
void mducresetode__(int **mem) { mducresetode_(mem); };

void mducsetsolution2d_(double *t, double *P, double *Y, double *T, int **mem)
{
  int i, j, k;
  int n1       = MDUCgetGridSize(1, *mem);
  int n2       = MDUCgetGridSize(2, *mem);
  int nSpecies = MDUCgetNspecs(*mem);

  // allocate local memory for species and convert to C order
  double ***YLocal = (double***) calloc(n2, sizeof(double**));
  double  **TLocal = (double**)  calloc(n2, sizeof(double*));
  for (k=0; k<n2; k++)
    {
      YLocal[k] = (double**) calloc(n1, sizeof(double*));
      TLocal[k] = (double*)  calloc(n1, sizeof(double*));
      for (j=0; j<n1; j++)
        YLocal[k][j] = (double*) calloc(nSpecies, sizeof(double));
    }
  MDUCspecArrayFtoC(nSpecies, n1, n2, Y, YLocal);
  convertArray1D_MD(n1, n2, 1, T, &TLocal);

  // set solution
  MDUCsetSolution(*t, *P, YLocal, TLocal, *mem);
  for (k=0; k<n2; k++)
    {
      for (j=0; j<n1; j++) free(YLocal[k][j]);
      free(YLocal[k]); free(TLocal[k]);
    }
  free(YLocal); free(TLocal);
}
void mducsetsolution2d__(double *t, double *P, double *Y, double *T, int **mem)
{ mducsetsolution2d_(t,P,Y,T,mem); }

void mducsetscalardiss1d_(int *n1, double *chi1, int **mem)
{ MDUCsetScalarDiss(*n1, 1, &chi1, NULL, NULL, *mem); }
void mducsetscalardiss1d__(int *n1, double *chi1, int **mem)
{ mducsetscalardiss1d_(n1,chi1,mem); }

void mducsetscalardiss2d_(int *n1, int *n2, double *chi1, double *chi12,
                          double *chi2, int **mem)
{ 
  int j, k;
  double **chi1L  = (double**) calloc(*n2, sizeof(double*));
  double **chi12L = (double**) calloc(*n2, sizeof(double*));
  double **chi2L  = (double**) calloc(*n2, sizeof(double*));
  for (k=0; k<*n2; k++)
    {
      chi1L[k]  = (double*) calloc(*n1, sizeof(double));
      chi12L[k] = (double*) calloc(*n1, sizeof(double));
      chi2L[k]  = (double*) calloc(*n1, sizeof(double));
    }

  /* scalar dissipation rates should be allocated in fortran as
       chiXX = allocate(n1,n2) */
  convertArray1D_MD(*n1, *n2, 1, chi1, &chi1L);
  convertArray1D_MD(*n1, *n2, 1, chi12, &chi12L);
  convertArray1D_MD(*n1, *n2, 1, chi2, &chi2L);

  MDUCsetScalarDiss(*n1, *n2, chi1L, chi12L, chi2L, *mem);
  for (k=0; k<*n2; k++) { free(chi1L[k]); free(chi12L[k]); free(chi2L[k]); }
  free(chi1L); free(chi12L); free(chi2L);
}
void mducsetscalardiss2d__(int *n1, int *n2, double *chi1, double *chi12,
                           double *chi2, int **mem)
{ mducsetscalardiss2d_(n1,n2,chi1,chi12,chi2,mem); }

void mducsetpressure_(int *n, double *P, int **mem)
{ MDUCsetPressure(*n, *P, *mem); }
void mducsetpressure__(int *n, double *P, int **mem)
{ mducsetpressure_(n,P,mem); }

void mducsetavgsdot_(int *m, int *n, double *Sdot, int **mem)
{ MDUCsetAvgSdot(*m, *n, Sdot, *mem); }
void mducsetavgsdot__(int *m, int *n, double *Sdot, int **mem)
{ mducsetavgsdot_(m,n,Sdot,mem); }

void mducsetavgh_(int *n, double *avgH, int **mem)
{ MDUCsetAvgH(*n, avgH, *mem); }
void mducsetavgh__(int *n, double *avgH, int **mem)
{ mducsetavgh_(n, avgH,mem); }

/* ========================================================================== */
/*                          chemistry functions                               */
/* ========================================================================== */
int  mducgetnspecies_(int **mem)  { return MDUCgetNspecies(*mem); }
int  mducgetnspecies__(int **mem) { return mducgetnspecies_(mem); }

int  mducgetnspecs_(int **mem)  { return MDUCgetNspecs(*mem); }
int  mducgetnspecs__(int **mem) { return mducgetnspecs_(mem); }

void mduccomputethermodata_(double *h, double *cp, double *T, int **mem)
{ MDUCcomputeThermoData(h, cp, *T, *mem); }
void mduccomputethermodata__(double *h, double *cp, double *T, int **mem)
{ mduccomputethermodata_(h, cp, T, mem); }

double mducgetmolweight_(char *species, int **mem)
{ return MDUCgetMolWeight(species, *mem); }
double mducgetmolweight__(char *species, int **mem)
{ return MDUCgetMolWeight(species, *mem); }

void mducgetmolweightall_(double *W,int **mem)  { MDUCgetMolWeightAll(W,*mem); }
void mducgetmolweightall__(double *W,int **mem) { mducgetmolweightall_(W,mem); }

void mducsetspeciesindex_(int nY, char **species, int *index, int **mem)
{ 
  int i;
  MDUCsetSpeciesIndex(nY, species, &index, *mem);
  for (i=0; i<nY; i++) index[i]++; // set to fortran index
}
void mducsetspeciesindex__(int nY, char **species, int *index, int **mem)
{ mducsetspeciesindex_(nY, species, index, mem); }

double mduccalctfromh_(int *nSpecies, double *H, double *Y, double *guess,
                       int **mem)
{ return MDUCcalcTfromH(*nSpecies,*H,Y,*guess,*mem); }
double mduccalctfromh__(int *nSpecies, double *H, double *Y, double *guess,
                        int **mem)
{ return mduccalctfromh_(nSpecies,H,Y,guess,mem); }

void mducsetsootmom_(int *n1, int *nSootMom, double *sootMom, int **mem)
{
  int i, j;
  
  /* allocate local memory for moments and convert to C order */
  double **mLocal = (double**) calloc(*n1, sizeof(double*));
  for (j=0; j<*n1; j++) mLocal[j] = (double*) calloc(*nSootMom, sizeof(double));
  MDUCspecArrayFtoC(*nSootMom, *n1, 1, sootMom, &mLocal);
  
  MDUCsetSootMom(*n1, *nSootMom, mLocal, *mem);
  for (j=0; j<*n1; j++) free(mLocal[j]); free(mLocal); 
}
void mducsetsootmom__(int *n1, int *nSootMom, double *sootMom, int **mem)
{ mducsetsootmom_(n1,nSootMom,sootMom,mem); }

void mducsetrhodot_(int *n1, double *rhoDot, int **mem)
{
  MDUCsetRhoDot(*n1, rhoDot, *mem);
}
void mducsetrhodot__(int *n1, double *rhoDot, int **mem)
{ mducsetrhodot_(n1,rhoDot,mem); }

void mducsetviscsrc_(int *n1, double *visc_src, int **mem)
{
  MDUCsetViscSrc(*n1, visc_src, *mem);
}
void mducsetviscsrc__(int *n1, double *visc_src, int **mem)
{ mducsetviscsrc_(n1,visc_src,mem); }
void mducsetglobalpdf_(int *n1, double *gpdf, int **mem)
{
  MDUCsetGlobalPDF(*n1, gpdf, *mem);
}
void mducsetglobalpdf__(int *n1, double *gpdf, int **mem)
{ mducsetglobalpdf_(n1,gpdf,mem); }

void mducmeansootsrc_(double *mean1, double *var1, double *mean2, double *var2,
                      double *var, double *meanSootSrc, int **mem)
{ MDUCcalcMeanSootSrc(*mean1,*var1,*mean2,*var2,*var,meanSootSrc,*mem); }
void mducmeansootsrc__(double *mean1, double *var1, double *mean2, double *var2,
                       double *var, double *meanSootSrc, int **mem)
{ mducmeansootsrc_(mean1,var1,mean2,var2,var,meanSootSrc,mem); }

void mducgetsootsrc_(int *nSootMom,double *Src, int **mem)
{
  int i, j;
  int n1       = MDUCgetGridSize(1, *mem);

  // allocate local memory for species
  double **SrcLocal = (double**) calloc(n1, sizeof(double*));
  for (j=0; j<n1; j++) SrcLocal[j] = (double*) calloc(*nSootMom, sizeof(double));

  MDUCgetSootSrc(SrcLocal, *mem);

  // convert arrays to fortran and free local memory
  MDUCspecArrayCtoF(*nSootMom, n1, 1, &SrcLocal, Src);
  for (j=0; j<n1; j++) free(SrcLocal[j]); free(SrcLocal);
}
void mducgetsootsolution__(int *nSootMom,double *Src, int **mem)
{ mducgetsootsrc_(nSootMom,Src,mem); }

void mducsoothmom_(int *nMom,double *T_SI, double *Visc_SI, double *Rho_SI, double *W_SI,
                  double *c_spec_SI, double *mom_SI, double *src_SI, double *rxnrates_SI,
                  double *rhodot_SI, double *srcP_SI, int *srcActive)
{

  MDUCsootHMOM(*nMom,*T_SI, *Visc_SI, *Rho_SI, *W_SI, c_spec_SI, mom_SI, src_SI, rxnrates_SI,
               rhodot_SI, srcP_SI, srcActive);
}

void mducsoothmom__(int *nMom,double *T_SI, double *Visc_SI, double *Rho_SI, double *W_SI,
                  double *c_spec_SI, double *mom_SI, double *src_SI, double *rxnrates_SI,
                  double *rhodot_SI, double *srcP_SI, int *srcActive)
{
  mducsoothmom_(nMom, T_SI, Visc_SI, Rho_SI, W_SI, c_spec_SI, mom_SI, src_SI, rxnrates_SI,
               rhodot_SI, srcP_SI, srcActive);
}
/* ========================================================================== */
/*                         main solver functions                              */
/* ========================================================================== */

void mducupdateinert1d_(double *dt, double *P, double *phimin, double *phimax,
                        double *Tox, double *Tfuel, double *chi, double *Y,
                        double *Temp, int **mem)
{
  int i, j;
  int n1       = MDUCgetGridSize(1, *mem);
  int nSpecies = MDUCgetNspecs(*mem);

  // allocate local memory for species
  double **YLocal = (double**) calloc(n1, sizeof(double*));
  for (j=0; j<n1; j++) YLocal[j] = (double*) calloc(nSpecies, sizeof(double));

  // call C solver
  MDUCupdateInert1D(*dt,*P,*phimin,*phimax,*Tox,*Tfuel,chi,YLocal,Temp,*mem);

  // convert arrays to fortran and free local memory
  MDUCspecArrayCtoF(nSpecies, n1, 1, &YLocal, Y);
  for (j=0; j<n1; j++) free(YLocal[j]); free(YLocal);
}
void mducupdateinert1d__(double *dt, double *P, double *phimin, double *phimax,
                         double *Tox, double *Tfuel, double *chi, double *Y,
                         double *Temp, int **mem)
{ mducupdateinert1d_(dt,P,phimin,phimax,Tox,Tfuel,chi,Y,Temp,mem); }

void mducflamelet1d_(double *dt,double *P,double *phimin,double *phimax,
                     double *chi,int *inert,double *Y,double *Temp,int **mem)
{
  int i, j;
  int n1       = MDUCgetGridSize(1, *mem);
  int nSpecies = MDUCgetNspecs(*mem);

  // allocate local memory for species
  double **YLocal = (double**) calloc(n1, sizeof(double*));
  for (j=0; j<n1; j++) YLocal[j] = (double*) calloc(nSpecies, sizeof(double));

  // call C solver
  MDUCflamelet1D(*dt, *P, *phimin, *phimax, chi, *inert, YLocal, Temp, *mem);

  // convert arrays to fortran and free local memory
  MDUCspecArrayCtoF(nSpecies, n1, 1, &YLocal, Y);
  for (j=0; j<n1; j++) free(YLocal[j]); free(YLocal);
}
void mducflamelet1d__(double *dt,double *P,double *phimin,double *phimax,
                     double *chi,int *inert,double *Y,double *Temp,int **mem)
{ mducflamelet1d_(dt, P, phimin, phimax, chi, inert, Y, Temp, mem); }

void mducsolve1d_(double *dt, double *P, double *phimin, double *phimax,
                  double *chi, double *Y, double *Temp, int **mem)
{ /* keep for compatibility */
  int inert = 0;
  mducflamelet1d_(dt, P, phimin, phimax, chi, &inert, Y, Temp, mem);
}
void mducsolve1d__(double *dt, double *P, double *phimin, double *phimax,
                   double *chi, double *Y, double *Temp, int **mem)
{ mducsolve1d_(dt, P, phimin, phimax, chi, Y, Temp, mem); }


void mducupdateinert2d_(double *dt, double *P, double *phi1min, double *phi1max,
                        double *phi2min, double *phi2max, double *Tox,
                        double *Tfuel1, double *Tfuel2, double *chi1,
                        double *chi12, double *chi2, double *Y, double *Temp,
                        int **mem)
{
  int    i, j, k;
  int    n1       = MDUCgetGridSize(1, *mem);
  int    n2       = MDUCgetGridSize(2, *mem);
  int    nSpecies = MDUCgetNspecs(*mem);

  // allocate local memory for C and copy values from 
  double ***YL    = (double***) calloc(n2, sizeof(double**));
  double **TempL  = (double**)  calloc(n2, sizeof(double*));
  double **chi1L  = (double**)  calloc(n2, sizeof(double*));
  double **chi12L = (double**)  calloc(n2, sizeof(double*));
  double **chi2L  = (double**)  calloc(n2, sizeof(double*));
  for (k=0; k<n2; k++)
    {
      YL[k]     = (double**) calloc(n1, sizeof(double*));
      TempL[k]  = (double*)  calloc(n1, sizeof(double));
      chi1L[k]  = (double*)  calloc(n1, sizeof(double));
      chi12L[k] = (double*)  calloc(n1, sizeof(double));
      chi2L[k]  = (double*)  calloc(n1, sizeof(double));
      for (j=0; j<n1; j++)
        YL[k][j] = (double*) calloc(nSpecies, sizeof(double));
    }

  /* scalar dissipation rates should be allocated in fortran as
       chiXX = allocate(n1,n2) */
  convertArray1D_MD(n1, n2, 1, chi1, &chi1L);
  convertArray1D_MD(n1, n2, 1, chi12, &chi12L);
  convertArray1D_MD(n1, n2, 1, chi2, &chi2L);

  // call 2D solver
  MDUCupdateInert2D(*dt,*P,*phi1min,*phi1max,*phi2min,*phi2max,*Tox,*Tfuel1,
                    *Tfuel2,chi1L,chi12L,chi2L,YL,TempL,*mem);

  FILE *pFile = fopen("c_species.dat", "w");

  for (i=0; i<nSpecies; i++) fprintf(pFile, "%15.6e\n", YL[4][9][i]);

  fclose(pFile);


  // convert arrays to fortran and free local memory
  MDUCspecArrayCtoF(nSpecies, n1, n2, YL, Y);
  convertArrayMD_1D(n1, n2, 1, &TempL, Temp);
  for (k=0; k<n2; k++)
    {
      for (j=0; j<n1; j++) free(YL[k][j]);
      free(YL[k]);free(TempL[k]);free(chi1L[k]);free(chi12L[k]);free(chi2L[k]);
    }
  free(YL); free(TempL); free(chi1L); free(chi12L); free(chi2L);
}
void mducupdateinert2d__(double *dt, double *P, double *phi1min,double *phi1max,
                         double *phi2min, double *phi2max, double *Tox,
                         double *Tfuel1, double *Tfuel2, double *chi1,
                         double *chi12, double *chi2, double *Y, double *Temp,
                         int **mem)
{
  mducupdateinert2d_(dt,P,phi1min,phi1max,phi2min,phi2max,Tox,Tfuel1,Tfuel2,
                     chi1,chi12,chi2,Y,Temp,mem);
}

void mducflamelet2d_(double *dt,double *P,double *phi1min,double *phi1max,
                     double *phi2min,double *phi2max,double *chi1,double *chi12,
                     double *chi2,int *inert,double *Y,double *Temp,int **mem)
{
  int    i, j, k;
  int    n1       = MDUCgetGridSize(1, *mem);
  int    n2       = MDUCgetGridSize(2, *mem);
  int    nSpecies = MDUCgetNspecs(*mem);

  // allocate local memory for C and copy values from 
  double ***YL    = (double***) calloc(n2, sizeof(double**));
  double **TempL  = (double**)  calloc(n2, sizeof(double*));
  double **chi1L  = (double**)  calloc(n2, sizeof(double*));
  double **chi12L = (double**)  calloc(n2, sizeof(double*));
  double **chi2L  = (double**)  calloc(n2, sizeof(double*));
  for (k=0; k<n2; k++)
    {
      YL[k]     = (double**) calloc(n1, sizeof(double*));
      TempL[k]  = (double*)  calloc(n1, sizeof(double));
      chi1L[k]  = (double*)  calloc(n1, sizeof(double));
      chi12L[k] = (double*)  calloc(n1, sizeof(double));
      chi2L[k]  = (double*)  calloc(n1, sizeof(double));
      for (j=0; j<n1; j++)
        YL[k][j] = (double*) calloc(nSpecies, sizeof(double));
    }

  /* scalar dissipation rates should be allocated in fortran as
       chiXX = allocate(n1,n2) */
  convertArray1D_MD(n1, n2, 1, chi1, &chi1L);
  convertArray1D_MD(n1, n2, 1, chi12, &chi12L);
  convertArray1D_MD(n1, n2, 1, chi2, &chi2L);

  // call 2D solver
  MDUCflamelet2D(*dt,*P,*phi1min,*phi1max,*phi2min,*phi2max,chi1L,chi12L,chi2L,
                 *inert,YL, TempL, *mem);

  // convert arrays to fortran and free local memory
  MDUCspecArrayCtoF(nSpecies, n1, n2, YL, Y);
  convertArrayMD_1D(n1, n2, 1, &TempL, Temp);
  for (k=0; k<n2; k++)
    {
      for (j=0; j<n1; j++) free(YL[k][j]);
      free(YL[k]);free(TempL[k]);free(chi1L[k]);free(chi12L[k]);free(chi2L[k]);
    }
  free(YL); free(TempL); free(chi1L); free(chi12L); free(chi2L);
}
void mducflamelet2d__(double *dt,double *P,double *phi1min,double *phi1max,
                      double *phi2min,double *phi2max,double *chi1,
                      double *chi12,double *chi2,int *inert,double *Y,
                      double *Temp,int **mem)
{  
  mducflamelet2d_(dt,P,phi1min,phi1max,phi2min,phi2max,chi1,chi12,chi2,inert,
                  Y,Temp,mem);
}

void mducsolve2d_(double *dt, double *P, double *phi1min, double *phi1max,
                  double *phi2min, double *phi2max, double *chi1, double *chi12,
                  double *chi2, double *Y, double *Temp, int **mem)
{
  mducflamelet2d_(dt,P,phi1min,phi1max,phi2min,phi2max,chi1,chi12,chi2,0,
                  Y,Temp,mem);
}
void mducsolve2d__(double *dt, double *P, double *phi1min, double *phi1max,
                   double *phi2min, double *phi2max, double *chi1,double *chi12,
                   double *chi2, double *Y, double *Temp, int **mem)
{
  mducsolve2d_(dt,P,phi1min,phi1max,phi2min,phi2max,chi1,chi12,chi2,Y,Temp,mem);
}

/* ========================================================================== */
/*                         statistics functions                               */
/* ========================================================================== */
int mduccalcpdf_(double *mean1, double *var1, double *mean2, double *var2,
                 double *var, double *pdf, int *pdfBound, int **mem)
{
  int j, k, pdfFlag;
  int n1 = MDUCgetGridSize(1, *mem);
  int n2 = MDUCgetGridSize(2, *mem);
  // allocate local 2D arrays for use with C call
  int    **pdfBoundL = (int**)    calloc(n2, sizeof(int*));
  double **pdfLocal  = (double**) calloc(n2, sizeof(double*));
  for (k=0; k<n2; k++)
    {
      pdfBoundL[k] = (int*)    calloc(2,    sizeof(int));
      pdfLocal[k]  = (double*) calloc(n2, sizeof(double));
    }

  pdfFlag = MDUCcalcPDF(*mean1,*var1,*mean2,*var2,*var,pdfLocal,pdfBoundL,*mem);

  // convert bounds to fortran index
  for (k=0; k<n2; k++) { pdfBoundL[0]++; pdfBoundL[1]++; }
  // convert 2D arrays to 2D fortran order
  convertArrayMD_1D(n1, n2, 1, &pdfLocal, pdf);
  /* TODO: write utility function for converting integer array */
  /* convertArrayMD_1D(2, n2, 1, &pdfBoundL, pdfBound);   */
  // free local memory
  for (k=0; k<n2; k++) { free(pdfLocal[k]); free(pdfBoundL[k]); }
  free(pdfLocal); free(pdfBoundL);
  return pdfFlag;
}
int mduccalcpdf__(double *mean1, double *var1, double *mean2, double *var2,
                  double *var, double *pdf, int *pdfBound, int **mem)
{ return mduccalcpdf_(mean1,var1,mean2,var2,var,pdf,pdfBound,mem); }

int mduccalcpdf_1d_(int *method, int *n, double *x, double *mean, double *var,
                    double *pdf, int *pdfBound)
{
  int pdfFlag;
  pdfFlag = MDUCcalcPDF_1D(*method, *n, x, *mean, *var, pdf, pdfBound);
  pdfBound[0]++; pdfBound[1]++; // convert indexing to fortran
  return pdfFlag;
}
int mduccalcpdf_1d__(int *method, int *n, double *x, double *mean, double *var,
                     double *pdf, int *pdfBound)
{ return mduccalcpdf_1d_(method,n,x,mean,var,pdf,pdfBound); }

int mduccalcmeany_1d_(int *method, int *n, double *x, int *nSpecies,
                      double *mean, double *var, double *Y, double *Ymean,
                      double *pdf, int *pdfBound)
{
  int j, pdfFlag;

  // allocate local memory for species and set to those input
  double **YLocal = (double**) calloc(*n, sizeof(double*));
  for (j=0; j<*n; j++) YLocal[j] = (double*) calloc(*nSpecies, sizeof(double));
  MDUCspecArrayFtoC(*nSpecies, *n, 1, Y, &YLocal);

  // compute the mean
  pdfFlag = MDUCcalcMeanY_1D(*method,*n,x,*nSpecies,*mean,*var,YLocal,Ymean,
                             pdf,pdfBound);
  // convert indexing to fortran
  pdfBound[0]++; pdfBound[1]++;

  // convert arrays to fortran and free local memory
  for (j=0; j<*n; j++) free(YLocal[j]); free(YLocal);
  return pdfFlag;
}
int mduccalcmeany_1d__(int *method, int *n, double *x, int *nSpecies,
                       double *mean, double *var, double *Y, double *Ymean,
                       double *pdf, int *pdfBound)
{ return mduccalcmeany_1d_(method,n,x,nSpecies,mean,var,Y,Ymean,pdf,pdfBound); }

int mduccalcmeany_2d_(int *method, int *n1, double *x1, int *n2, double *x2,
                      int *nSpecies, double *mean1, double *var1,
                      double *mean2, double *var2, double *var,
                      double *Y, double *Ymean, double *pdf, int *transform)
{
  int j, k, pdfFlag;

  // allocate local variables for grid
/*   int *n     = (int*)     calloc(2, sizeof(int)); */
  int n[2];
  double **x = (double**) calloc(2, sizeof(double*));
  n[0] = *n1; n[1] = *n2;
  x[0] = x1;  x[1] = x2;

  // allocate local memory for species and set to those input
  double ***YLocal = (double***) calloc(*n2, sizeof(double**));
  for (k=0; k<*n2; k++) {
    YLocal[k] = (double**) calloc(*n1, sizeof(double*));
    for (j=0; j<*n1; j++) {
      YLocal[k][j] = (double*) calloc(*nSpecies, sizeof(double));
    }
  }
  MDUCspecArrayFtoC(*nSpecies, *n1, *n2, Y, YLocal);

  // local pdf memory
  int **pdfBound = (int**) calloc(*n2, sizeof(int*));
  double **pdfL = (double**) calloc(*n2, sizeof(double*));
  for (k=0; k<*n2; k++) {
    pdfL[k]     = (double*) calloc(*n1, sizeof(double));
    pdfBound[k] = (int*)    calloc(2, sizeof(int));
  }

  // compute the mean species
  pdfFlag = MDUCcalcMeanY_2D(*method,n,x,*nSpecies,*mean1,*var1,*mean2,*var2,
                             *var,YLocal,Ymean,pdfL,pdfBound,*transform);
  // convert indexing to fortran
  for (k=0; k<*n2; k++) { pdfBound[k][0]++; pdfBound[k][1]++; }

  // convert arrays to fortran and free local memory
  convertArrayMD_1D(*n1, *n2, 1, &pdfL, pdf);

  for (k=0; k<*n2; k++) {
    for (j=0; j<*n1; j++) free(YLocal[k][j]);
    free(YLocal[k]); free(pdfL[k]); free(pdfBound[k]);
  }
  free(YLocal); free(pdfL); free(pdfBound); free(x);
  return pdfFlag;
}
int mduccalcmeany_2d__(int *method, int *n1, double *x1, int *n2, double *x2,
                      int *nSpecies, double *mean1, double *var1,
                      double *mean2, double *var2, double *var,
                      double *Y, double *Ymean, double *pdf, int *transform)
{
  return mduccalcmeany_2d_(method,n1,x1,n2,x2,nSpecies,mean1,var1,mean2,var2,
                           var,Y,Ymean,pdf,transform);
}

int mduccalcpdf_2d_(int *method, int *n1, double *x1, int *n2, double *x2,
                    double *mean1, double *var1, double *mean2, double *var2,
                    double *var, double *pdf, int *transform)
{
  int j, k, pdfFlag;

  // assign grid parameters to form needed for C call
  int     *n = (int*) calloc(2, sizeof(int));
  double **x = (double**) calloc(2, sizeof(double*));
  n[0] = *n1; n[1] = *n2;
  x[0] = x1;  x[1] = x2;

  // allocate local 2D arrays for use with C call
  int    **pdfBoundL = (int**)    calloc(n[1], sizeof(int*));
  double **pdfLocal  = (double**) calloc(n[1], sizeof(double*));
  for (k=0; k<n[1]; k++)
    {
      pdfBoundL[k] = (int*)    calloc(2,    sizeof(int));
      pdfLocal[k]  = (double*) calloc(n[0], sizeof(double));
    }

  // C call
  pdfFlag = MDUCcalcPDF_2D(*method, n, x, *mean1, *var1, *mean2, *var2, *var,
			   pdfLocal, pdfBoundL, *transform);
  // convert to fortran index
  for (k=0; k<n[1]; k++) { pdfBoundL[k][0]++; pdfBoundL[k][1]++; }

  // convert 2D arrays to 2D fortran order
  convertArrayMD_1D(*n1, *n2, 1, &pdfLocal, pdf);

  // free local memory
  for (k=0; k<n[1]; k++) { free(pdfLocal[k]); free(pdfBoundL[k]); }
  free(pdfLocal); free(pdfBoundL);
  free(n); free(x);
  return pdfFlag;
}
int mduccalcpdf_2d__(int *method, int *n1, double *x1, int *n2, double *x2,
                     double *mean1, double *var1, double *mean2, double *var2,
                     double *var, double *pdf, int *transform)
{
  return mduccalcpdf_2d_(method,n1,x1,n2,x2,mean1,var1,mean2,var2,var,pdf,
                         transform);
}

void mducinterpchidns_(int *p1, double *chi1, int *p2, double *chi12,
		       int *p3, double *chi2, int **mem)
{
  int i, j, k;
  int n1       = MDUCgetGridSize(1, *mem);
  int n2       = MDUCgetGridSize(2, *mem);
  // allocate local memory and fill with values
  double **chi1L = NULL, **chi12L = NULL, **chi2L = NULL;
  if (*p1)
    {
      chi1L = (double**) calloc(n2, sizeof(double*));
      for (k=0; k<n2; k++)
        chi1L[k] = (double*) calloc(n1, sizeof(double));

      convertArray1D_MD(n1, n2, 1, chi1, &chi1L);
    }
  if (*p2)
    {
      chi12L = (double**) calloc(n2, sizeof(double*));
      for (k=0; k<n2; k++)
        chi12L[k] = (double*) calloc(n1, sizeof(double));

      convertArray1D_MD(n1, n2, 1, chi12, &chi12L);
    }
  if (*p3)
    {
      chi2L = (double**) calloc(n2, sizeof(double*));
      for (k=0; k<n2; k++)
        chi2L[k] = (double*) calloc(n1, sizeof(double));

      convertArray1D_MD(n1, n2, 1, chi2, &chi2L);
    }

  MDUCinterpChiDNS(chi1L, chi12L, chi2L, *mem);
  
  // convert arrays to fortran and free local memory
  convertArrayMD_1D(n1, n2, 1, &chi1L, chi1);
  for (k=0; k<n2; k++)
    {
      if (chi1L) free(chi1L[k]);
      if (chi2L) free(chi2L[k]);
    }
  if (chi1L) free(chi1L);
  if (chi2L) free(chi2L);
 
}
void mducinterpchidns__(int *p1, double *chi1, int *p2, double *chi12,
		       int *p3, double *chi2, int **mem)
{
  mducinterpchidns_(p1, chi1, p2, chi12, p3, chi2, mem);
}

/* ========================================================================== */
/*                              I/O functions                                 */
/* ========================================================================== */
void mducwriteflamelet_(char *file,int **mem)  { MDUCwriteFlamelet(file,*mem); }
void mducwriteflamelet__(char *file,int **mem) { mducwriteflamelet_(file,mem); }

void mducwriterates_(char *file,int **mem)  { MDUCwriteRates(file,*mem); }
void mducwriterates__(char *file,int **mem) { mducwriterates_(file,mem); }

void mducwriteflowdatamem_(char *file,int **mem)  { MDUCwriteFlowDataMem(file,*mem); }
void mducwriteflowdatamem__(char *file,int **mem) { mducwriteflowdatamem_(file,mem); }

void mducwriteflowdata_(int *p1, double *chi1, int *p2, double *chi2, int *p3,
                        double *pdf, int *p4, double *sdot, int **mem)
{
  int j, k;
  int n1 = MDUCgetGridSize(1, *mem);
  int n2 = MDUCgetGridSize(2, *mem);

  // allocate local memory and fill with values
  double **chi1L = NULL, **chi2L = NULL, **pdfL = NULL, *sdotL = NULL;
  if (*p1)
    {
      chi1L = (double**) calloc(n2, sizeof(double*));
      for (k=0; k<n2; k++)
        chi1L[k] = (double*) calloc(n1, sizeof(double));

      convertArray1D_MD(n1, n2, 1, chi1, &chi1L);
    }
  /* TODO: add chi12 */
  if (*p2)
    {
      chi2L = (double**) calloc(n2, sizeof(double*));
      for (k=0; k<n2; k++)
        chi2L[k] = (double*) calloc(n1, sizeof(double));

      convertArray1D_MD(n1, n2, 1, chi2, &chi2L);
    }
  if (*p3)
    {
      pdfL = (double**) calloc(n2, sizeof(double*));
      for (k=0; k<n2; k++)
        pdfL[k] = (double*) calloc(n1, sizeof(double));

      convertArray1D_MD(n1, n2, 1, pdf, &pdfL);
    }
  if (*p4) sdotL = sdot;

  MDUCwriteFlowData(chi1L, chi2L, pdfL, sdotL, *mem);

  for (k=0; k<n2; k++)
    {
      if (chi1L) free(chi1L[k]);
      if (chi2L) free(chi2L[k]);
      if (pdfL)  free(pdfL[k]);
    }
  if (chi1L) free(chi1L);
  if (chi2L) free(chi2L);
  if (pdfL)  free(pdfL);
}

void mducwritescalars_(double *ctime, double *pdf, int **mem) 
{
  int j, k;
  int n1 = MDUCgetGridSize(1, *mem);
  int n2 = MDUCgetGridSize(2, *mem);

  if (n2 == 1)
    { /* 1D solution */
      int *pdfBound = (int*) calloc(2, sizeof(int));
      pdfBound[0] = 0; pdfBound[1] = n1-1;
      
      MDUCwriteScalars(*ctime, NULL, NULL, &pdf, &pdfBound, -1, *mem);
      free(pdfBound);
    }
  else /* 2D solution */
    {
      int **pdfBound = (int**)    calloc(n2, sizeof(int*));
      double **pdfL  = (double**) calloc(n2, sizeof(double*));
      for (k=0; k<n2; k++)
        {
          pdfBound[k] = (int*)    calloc(2,  sizeof(int));
          pdfL[k]     = (double*) calloc(n1, sizeof(double));
        }

      convertArray1D_MD(n1, n2, 1, pdf, &pdfL);

      MDUCwriteScalars(*ctime, NULL, NULL, pdfL, pdfBound, -1, *mem);

      for (k=0; k<n2; k++) { free(pdfBound[k]); free(pdfL[k]); }
      free(pdfBound); free(pdfL);
    } /* end 2D */
}
void mducwritescalars__(double *ctime, double *pdf, int **mem)
    { mducwritescalars_(ctime,pdf,mem); }

/* ========================================================================== */
/*                           Miscellaneous functions                          */
/* ========================================================================== */

double mduccalcchilogform_(double *Z, double *Zref, double *Zmax)
{ return MDUCcalcChiLogForm(*Z, *Zref, *Zmax); }
double mduccalcchilogform__(double *Z, double *Zref, double *Zmax)
{ return mduccalcchilogform_(Z,Zref,Zmax); }

double mduccalcchierfcform_(double *Z, double *Zref, double *Zmax)
{ return MDUCcalcChiErfcForm(*Z, *Zref, *Zmax); }
double mduccalcchierfcform__(double *Z, double *Zref, double *Zmax)
{ return mduccalcchierfcform_(Z,Zref,Zmax); }

/* ========================================================================== */

void MDUCintArrayFtoC(int n1, int n2, int n3, int *fval, int ***cval)
{
  int i, j, k;
  for (i=0; i<n1; i++)
    for (j=0; j<n2; j++)
      for (k=0; k<n3; k++)
        cval[i][j][k] = fval[i+n1*(j+n2*k)];
}
/* -------------------------------------------------------------------------- */

void MDUCintArrayCtoF(int n1, int n2, int n3, int ***cval, int *fval)
{
  int i, j, k;
  for (i=0; i<n1; i++)
    for (j=0; j<n2; j++)
      for (k=0; k<n3; k++)
        fval[i+n1*(j+n2*k)] = cval[i][j][k];
}
/* -------------------------------------------------------------------------- */

void MDUCrealArrayFtoC(int n1, int n2, int n3, double *fval, double ***cval)
{
  int i, j, k;
  for (i=0; i<n1; i++)
    for (j=0; j<n2; j++)
      for (k=0; k<n3; k++)
        cval[i][j][k] = fval[i+n1*(j+n2*k)];
}
/* -------------------------------------------------------------------------- */

void MDUCrealArrayCtoF(int n1, int n2, int n3, double ***cval, double *fval)
{
  int i, j, k;
  for (i=0; i<n1; i++)
    for (j=0; j<n2; j++)
      for (k=0; k<n3; k++)
        fval[i+n1*(j+n2*k)] = cval[i][j][k];
}
/* -------------------------------------------------------------------------- */

void MDUCspecArrayFtoC(int nSp, int nZ1, int nZ2, double *fval, double ***cval)
{
  int Z1, Z2, sp;
  for (Z2=0; Z2<nZ2; Z2++)
    for (Z1=0; Z1<nZ1; Z1++)
      for (sp=0; sp<nSp; sp++)
        cval[Z2][Z1][sp] = fval[Z1+nZ1*(Z2+nZ2*sp)];
}
/* -------------------------------------------------------------------------- */

void MDUCspecArrayCtoF(int nSp, int nZ1, int nZ2, double ***cval, double *fval)
{
  int Z1, Z2, sp;
  for (Z2=0; Z2<nZ2; Z2++)
    for (Z1=0; Z1<nZ1; Z1++)
      for (sp=0; sp<nSp; sp++)
        fval[Z1+nZ1*(Z2+nZ2*sp)] = cval[Z2][Z1][sp];
}
/* -------------------------------------------------------------------------- */
