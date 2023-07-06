#ifndef MDUC_MAIN_H
#define MDUC_MAIN_H

#include "mducStructs.h"

struct OPTSTRUCT {
  /* ----- for getting command line options and storing ----- */
  int  nOpts;            // total number of options
  int  nFiles;           // number of files to process
  char *args;            // command line arguments
  char mducFile[50];     // mduc input file to read from
  char **dataFile;       // data files
  char postProcFile[50]; // file with post-process info
  char monitorFile[100]; // monitor file
  char outFile[100];     // general file for output
  char ensDir[100];      // ensight directory
  char ensCase[100];     // existing ensight case file
  int  range[2];         // range for files
  int  nChiMethod;
  char **chiMethod;
  int fileMin, fileMax;

  /* ----- flags for determining what to do ----- */
  int postProcess;     // postProcess a data file
  int fuelTemp;        // compute the fuel temperature for liquid fuels
  int interp;          // interpolate the solution read from the file
  int init2D;
  int convert1D, multi1D;
  int writeEnsight, ensConst, ensTime, ensConv, freq;
  int checkGrid, checkChi, postChi, chiProfile, checkThermo;
  int integratePDF, tFlow, hrBudg, clean, combine;
  int iT, iV;
  int metric, metricOrder;
  int polyA;
  int listSpecies, fireSpecies, induction, intQ;
  int checkVTK;
  int addWmol;
  int iFile;
  int testParse;

  /* ----- options read in ----- */
  int     nPivot, nGrid;
  double *pivot, *gridPerc;

  /* ----- variables read in from command line ----- */
  int dir, n1, n2;
  double *gen;      // generic vector to read in variables
};
typedef struct OPTSTRUCT optStruct;

/* functions to retrieve internal structs from solver */
ioStruct     * MDUCgetIOStruct(void *mem);
driveStruct  * MDUCgetDriveStruct(void *mem);
timeStruct   * MDUCgetTimeStruct(void *mem);
chemStruct   * MDUCgetChemStruct(void *mem);
flameStruct  * MDUCgetFlameStruct(void *mem);
flowStruct   * MDUCgetFlowStruct(void *mem);
engineStruct * MDUCgetEngineStruct(void *mem);
void getMDUCStructs(ioStruct **io, driveStruct **driver, chemStruct **chem,
                    timeStruct **time, flowStruct **flow, flameStruct **flame,
                    void *mem);

#ifdef ACCESS_MPI
void getMPIStruct(mpiStruct **mpi, void *mem);
#endif /* ACCESS_MPI */

// read options from the command line
void readOptions(int argc, char **argv, optStruct *opts);
void readList(int n, int iStart, char ***var, char **argv);

#endif /* MDUCMAIN_H */
