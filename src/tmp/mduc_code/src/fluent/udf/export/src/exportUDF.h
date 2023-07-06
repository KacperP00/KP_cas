#ifndef MDUC_UDF_EXPORT_H
#define MDUC_UDF_EXPORT_H

#include "udf.h"
#include "prf.h"
#include "mpi.h"
#include "version.h"
#include "surf.h"
#include "cxsurf.h"
#include "cxiface.h"

/* struct for storing UDF library options */
typedef struct MDUC_EXPORT_STRUCT {
  /* control flags */
  int single_zone;     /* use single cell part? */
  int single_dir;      /* use single directory? */
  int dump_velocity;   /* export velocity? */
  int dump_species;    /* export species mass fractions? */
  int dump_sdm;        /* export single droplet parameters? */
  int dump_particle;   /* export particles? */
  int dump_uds;        /* export user-defined-scalars? */
  int dump_udm;        /* export user-defined-memory? */
  int binary;          /* write binary format? */
  int poly_mesh;       /* force a polyhedra mesh*/
  int byteOrder;       /* 0: little endian, 1: big endian */
  int restart;         /* is this a restart? */
  int init;            /* is this initial time-step? */
  int export;          /* should data be exported? */
  int format;          /* 0: ensight, 1: vtk */
  int inLine;          /* for VTK inline or appended */
  char vtkFormat[32];  

  /* info about grouping zones together */
  int nZone, *nZoneThreads;
  char **zones, ***zThreadNames;
  Thread ***zThreads;
  int nSurf, *nSurfThreads, nAllSurfThreads, *sThreadIDs;
  char **surfs, ***sThreadNames;
  Thread ***sThreads;
  int nCxSurf, *nCxSurfThreads, nAllCxThreads, *cxThreadIDs;
  char **cxSurfs, ***cxThreadNames;
  int  nCxRef, *cxRefIds;
  char **cxRefNames;

  /* memory for flamelet solver */
  int     nSc, nY, nCells, nWalls, nParts, freq, nInj, nNoz;
  Thread  **cell_th, **wall_th;
  Injection ***streams, **all_streams;
  int    *part_id;
  char  **part_name, ***part_type;
  int    *part_nv, **part_verts;
  char  **inject_id;

  char dir[1000], prefix[1000];

  int fileNo; /* timestep index for ensight files */

  /* for scalar identification */
  int   *sc_indx, *sc_type, *Y_indx;
  char **sc_names, **sc_abbrv, **Y_names;
};
typedef struct MDUC_EXPORT_STRUCT expStruct;

/* defined functions */
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;

/* functions to call from outside */
void UDFinitExport(void);

/* functions for writing geometry in parts */
int initThreads(Thread ***tF, int **ids, char ***name, char ****type, int **nv,
                int ***verts);
int getFluidThreadNames(char ***names);
int getWallThreadNames(char ***names);
void checkCurrentZones(int force);
int countFluidThreads();
int getNumberCellTypes(Thread *t, int **verts, char ***types);
void destroyExportStruct(expStruct *e);

/* functions for writing binary strings */
#if PARALLEL && MDUC_MPI_WRITE
void writeBinStr(MPI_File *pFile, int len, char *str);
void writeBinInt(MPI_File *pFile, int len, int *var);
void writeBinReal(MPI_File *pFile, int len, float *var);
#else /* SERIAL OR ROOT WRITE */
void (*writeDataStr)(FILE*, int, char*);
void (*writeDataInt)(FILE*, int, int*, int);
void (*writeDataReal)(FILE*, int, float*, int);
int  (*writeBinaryData)(void *buffer, int size, int count, FILE *pFile);
int  BINARY_FILE_WRITE(void *buffer, int size, int count, FILE *pFile);
void writeBinStr(FILE *pFile, int len, char *str);
void writeBinInt(FILE *pFile, int len, int *var, int perLine);
void writeBinReal(FILE *pFile, int len, float *var, int perLine);
void writeAsciiStr(FILE *pFile, int len, char *str);
void writeAsciiInt(FILE *pFile, int len, int *var, int perLine);
void writeAsciiReal(FILE *pFile, int len, float *var, int perLine);
#endif /* PARALLEL */

int parseFluentString_EX(char *input, char ***var);
void getInjection(char *name, int nNozzle, Injection ***I);
void indexIntArray(int n, int *arr, int *indx);
int gridLocateIntIndx(int n1, int *x1, int *indx, int x);
void EXstripString(char *str);
void EXstripStringBraces(char *name);
int EXparseThreadNames(int nRaw, char **raw, int *nT, char ***names);
int EX_getCxId(char *name);
int EX_RP_is_Defined(char *var_name);
int EX_checkEndian();
int IMIN(int x, int y);
int IMAX(int x, int y);
double ABS_D(double x);
#endif /* MDUC_UDF_EXPORT_H */
