#ifndef MDUC_POSTPROCESS_H
#define MDUC_POSTPROCESS_H

#include "mducStructs.h"

#ifdef ACCESS_MPI
#include "mpiStruct.h"
#endif // ACCESS_MPI

/* ----- functions for post-processing data file ----- */
void postProcess1D(char *postFile, ioStruct *io, driveStruct *driver,
                   void *workMDUC);
void postProcess2D(char *postFile, ioStruct *io, driveStruct *driver,
                   void *workMDUC);
void postProcessMulti1D(char *postFile, ioStruct *io, driveStruct *driver, 
                        void *mem, void *mem1, void *mem2);

int inputPostProcess(int **nt, double **tIn, double **P, double ***chiR,
                     double ***ZR, double ***Zmean, double ***Zvar,
                     double ***xmin, double ***xmax, void *mem);

/* ----- I/O functions ----- */
int readFluentData(char *dataFile, int nYFav, char **species, double **P,
                   double **RHO, double ***YFav, double **HFav, double **K,
                   double **D, double **Z1Fav, double **Z1Var, double **Z2Fav,
                   double **Z2Var, double **Z12Var, double **V);

int readFileList(char *fileForm, int iStart, int iStop, int **nt, double **tIn,
                 double **P, double ***chiR, double ***ZR, double ***Zmean,
                 double ***Zvar, double ***xmin, double ***xmax, void *mem);

int readScalarFile(char *file, int **nt, double **tIn, double **P, double ***T,
                   double ***chiR, double ***ZR, double ***Zmean,
                   double ***Zvar, double ***xmin, double ***xmax, void *mem);

// write a file comparing values with DNS
void writeErrorFileHeader(char *fileName, int dim);
void writeErrorFile(char *fileName, double t, double P, double P_DNS,
                    double dns1min, double dns1max, double dns2min,
                    double dns2max, void *mem);

#endif /* MDUC_POSTPROCESS_H */
