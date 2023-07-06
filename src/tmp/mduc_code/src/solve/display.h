#ifndef MDUC_DISPLAY_H
#define MDUC_DISPLAY_H

#include "mducStructs.h"

int MDUC_PRINT_PROC;

typedef struct monStruct_ {
  char  id[256];
  FILE  *p;
  int    n;
  int   *type;  /* 0: integer, -1 exponential, > 0, real precision */
  char **names;
  char **units;
  void **vals;
} monStruct;

/* ----- error/output functions ----- */
void mducHeader(char *head);
void mducInfo(char *message, char *func);
void mducErr(char *message, char *func);
void mducWarn(char *message, char *func);

/* ----- for writing monitor files ----- */
monStruct * monitorInit(char *id);
void monitorRegister(char *name,char *unit,void *val,int type,monStruct *mon);
void monitorWriteHead(monStruct *mon);
void monitorWriteLine(monStruct *mon);
void destroyMonitor(monStruct *mon);

void printInput(ioStruct *io, driveStruct *driver, chemStruct *chem,
                timeStruct *time, flameStruct *flame);

void printArray(int n, int index, int offset, double *y);

void printSpeciesList(chemStruct *chem);

void printStreamProperties(ioStruct *io);

void printMonitor(long int nst, double t, double dtPre, ioStruct *io, double *Y,
                  flameStruct *flame, double dt_max);

#endif /* MDUC_DISPLAY_H */
