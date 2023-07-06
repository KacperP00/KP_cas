#ifndef MDUC_PARSER_H
#define MDUC_PARSER_H

typedef struct parseStruct_ {
  int  nIn;
  char *comment;
  char *delims;
  int  *n;
  char **tags;
  char ***vals;
} parseStruct;

parseStruct * parseFile(char *fileName);
void   destroyParser(parseStruct *p);
int    parsePrint(parseStruct *p);
int    paramGetIndex(char *name, parseStruct *p);
int    paramGetInt(char *name, int def, parseStruct *p);
double paramGetReal(char *name, double def, parseStruct *p);
void   paramGetString(char *val, char *name, char *def, parseStruct *p);
int    paramGetStringArray(char ***val, char *name, parseStruct *p);
int    paramGetIntArray(int **out, char *name, parseStruct *p);
int    paramGetRealArray(double **out, char *name, parseStruct *p);
int    sparamGetIndex(int n, char *name, char **vals);
int    sparamGetInt(int n, char *name, int def, char **vals);
double sparamGetReal(int n, char *name, double def, char **vals);
int    sparamGetIntArray(int n, char *name, int **out, char **vals);
int    sparamGetRealArray(int n, char *name, double **out, char **vals);
int    sparamGetStringArray(int n, char *name, char ***out, char **vals);

/* 
 functions to read in input stream from a file
 function should be self-explanatory
*/
/* int    sread(char *str, char *fmt, int len, int n, void *out); */
int    sreadInt(char *str, int len, int n, int *out);
int    sreadReal(char *str, int len, int n, double *out);
int    sreadString(char *str, int len, char *out);
int    sreadStringArray(char *str, int len, int n, char **out);
int    readInteger(FILE *pFile);
double readReal(FILE *pFile);
void   readString(FILE *pFile, char *var);
void   readRealArray(FILE *pFile, int n, double **var);
void   readIntArray(FILE *pFile, int n, int var[]);
void   readOptArray(FILE *pFile, int n, double var[]);

/*
 Gets the number of boundary conditions as defined in an input or flamelet file
 determines number of species define as boundary conditions
 NOTE: does not count temperature
*/
int  getNumBC(FILE *pFile, int *fTemp, int *fH);
int  getNumFuelBC(FILE *pFile, int *fTemp, int *fH, int *fHvap, int *fCpL,
                  int *fTref, int *fSat);

/*
 Reads in the boundary conditions defined in input or flamelet data file.
 reads in temperature plus n species, which are determined by getNumBC.
 dynamically allocates proper memory
 stores species names, mass fractions, and temperature
*/
char readBC(FILE *pFile, int *nComp, char ***ident, double *T, double **Y,
            double *H);
char readFuelBC(FILE *pFile, int *nComp, char ***ident, double *T, double **Y,
                double *H, double *Hvap, double *cpL, double *Tref, char *data);

/*
 Reads in and stores a list of species names to be output in the data files
 dynamically allocates proper memory
*/
int  readNames(FILE *pFile, char ***ident);
int  readRealArraySize(FILE *pFile, double **var);

/*
 set of functions used for old format of Dave's input file
 NOTE: most likely no longer necessary
 */
int  getNumComponents(FILE *pFile);
//void readComposition(FILE *pFile, int *nComp, char ***ident, double **y_species);

#endif /* MDUC_PARSER_H */
