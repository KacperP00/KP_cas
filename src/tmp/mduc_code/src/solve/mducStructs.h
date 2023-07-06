#ifndef MDUC_STRUCT_H
#define MDUC_STRUCT_H

#include "ioStruct.h"
#include "driveStruct.h"
#include "chemStruct.h"
#include "timeStruct.h"
#include "flowStruct.h"
#include "flameStruct.h"
#include "engineStruct.h"
#include "sdialStruct.h"
#include "rk4struct.h"

#ifdef ACCESS_MPI
#include "mpiStruct.h"
#endif // ACCESS_MPI

struct MDUCSTRUCT {
  ioStruct     *io;      /* input/output storage */
  driveStruct  *driver;  /* program options */
  chemStruct   *chem;    /* chemical mechanism storage */
  timeStruct   *time;    /* time information */
  flowStruct   *flow;    /* physical flow variables */
  flameStruct  *flame;   /* flamelet solution */
  sdialStruct  *sdial;   /* SUNDIALS solver storage */
  engineStruct *engine;  /* engine parameters */
  rkStruct    *rk;
#ifdef ACCESS_MPI
  mpiStruct    *mpi;     /* parallel processing information */
#endif /* ACCESS_MPI */
};
typedef struct MDUCSTRUCT mducStruct;

#endif /* MDUC_STRUCT_H */
