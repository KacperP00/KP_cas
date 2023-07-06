#ifndef MDUC_CHEM_BIN_H
#define MDUC_CHEM_BIN_H

#include "chemStruct.h"

struct OPTSTRUCT {
  /* ----- options ----- */
  int  nOpts;            /* total number of options */
  char *args;            /* command line arguments */

  /* ----- data/files ----- */
  int  nFiles;           /* size of data file list */
  char **dataFile;       /* list of data files */
  char inputFile[128];   /* chemistry input file */

  /* ----- mechanism info ----- */
  char mechFile[128];    /* mechanism file */
  char mechBase[128];    /* base name of mechanism */
  char mechSuffix[16];   /* mechanism suffix */
  char thermFile[128];   /* thermodynamic data file */
  char transFile[128];   /* species transport data file */
  char format[32];       /* mechanism format */
  char output[128];      /* output file base name */

  int  info;             /* write mech info */
  int  convert;          /* convert mechanism format */
  char convFormat[32];   /* format to convert to */
  int  chemTable;        /* chemtable functionality */
  char chemTabFile[256], chemTabType[32]; 
  int  chemTabDim;
  int  chemTabInfo;
  int  chemTabLookup;
  int  chemTabLookupExp;
  int  cea;
  double *chemTabX;      /* location to lookup chemtable values */
};
typedef struct OPTSTRUCT optStruct;

void readOptions(int argc, char **argv, optStruct *opts);

/* convert between mechanism formats */
void convertMechanism(optStruct *opts);

/* chemtable generation */
int chemTable_fitP(optStruct *opts);

#endif /* MDUC_CHEM_BIN_H */
