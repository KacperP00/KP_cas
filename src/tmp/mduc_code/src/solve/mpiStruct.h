#ifndef MPI_STRUCT_H
#define MPI_STRUCT_H

#include "mpi.h"

struct MPISTRUCT {
  MPI_Comm world;   /* global communicator */
  MPI_Comm comm;    /* local communicator */
  MPI_Group worldG; /* world group */
  MPI_Group localG; /* reduced set of processors */

  int nWorld; 
  int nProc;  /* number of processors available */
  int myid;   /* ident of current process */
  int iWorld;
  int host;   /* ident of host node */
  int rxn;    /* last node used for comp of rxn */

  /* info about neighbouring processors */
  int nBcomm, *bComm, *bSend, *bRcv, *bCount;

  /* info about sendrecv comm */
  int nSRcomm, *srComm, *srSender, *srRecver, *srCount,*srIndex;

  /* number of rows/columns to be solved by process */
  int *nRowLocal, *nColLocal;

  /* start and ending rows/columns */
  int *rowStart, *rowEnd, *colStart, *colEnd;

  /* ----- indices/length for arrays with  with one variable ----- */
  /* for rows/columns assigned to processor (i.e., for gathering) */
  int *rStart, *rCount; /* row start index and length */
  int *cStart, *cCount; /* column start index and length */

  /* when plus/minus a row/column is needed for explicit calculations */
  int *rPlus, *rCountPlus;
  int *cPlus, *cCountPlus;

  /* ----- indices/length for solution arrays (multiple variables) ----- */
  /* for row/column solutions assigned to processor (i.e. for gathering) */
  int *kStart, *kCount; /* row start index and length */
  int *jStart, *jCount; /* column start index and length */

  /* when plus/minus a row/column is needed for explicit calculations */
  int *kPlus, *kCountPlus;
  int *jPlus, *jCountPlus;

  /* for soot moments 
     TODO: this should maybe be eliminated */
  int *soot0, *sootN;
  int *soot0xmyp, *sootNxmyp;
  int *sootSrcP0, *sootSrcPN;
};
typedef struct MPISTRUCT mpiStruct;

#endif /* MPI_STRUCT_H */
