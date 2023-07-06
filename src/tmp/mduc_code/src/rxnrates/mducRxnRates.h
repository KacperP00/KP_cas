#ifndef MDUC_POSTREAC_H
#define MDUC_POSTREAC_H

#include "mducStructs.h"

/* struct to hold option read from command line */
struct OPTSTRUCT {
  /* ----- for getting command line options and storing ----- */
  int  nOpts;          // total number of options
  int  nFiles;         // number of files to process
  int  nPAH;	       // number of PAH molecules parsed from cmd line
  char *args;          // command line arguments
  char suffix[10];     // file suffix to use
  char **postFileName; // list of filenames to process
  char **PAHlist; // list of PAH molecules parsed from cmd line
  char outFile[100];   // output file

};
typedef struct OPTSTRUCT optStruct;

// function write the flamelet solution read in to a tabulated format
void writePostFile(char *fileName, optStruct *opts, chemStruct *chem);


// read options from the command line
void readOptions(int argc, char **argv, optStruct *opts);
void readList(int n, int iStart, char ***var, char **argv);

// generate a file name with different suffix
void createOutputFileName(char *prefix, char *fileName, char *ext,
                          char *outputFile);


#endif /* MDUC_POSTREAC_H */
