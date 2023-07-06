#ifndef MDUC_EOS_H
#define MDUC_EOS_H

#include "eosFluid.h"

typedef struct optStruct_ {
  /* ----- for getting command line options and storing ----- */
  int  nOpts;            /* total number of options */
  int  nFiles;           /* number of files to process */
  char *args;            /* command line arguments */
  char inFile[50];       /* input file to read from */
  char outFile[50];      /* output file to write to */

  /* ----- problem type ----- */
  int PT, rhoT, rhoP;

  /* ----- store initial properties ----- */
  

} optStruct;

void readOptions(int argc, char **argv, optStruct *opts);

#endif /* MDUC_EOS_H */
