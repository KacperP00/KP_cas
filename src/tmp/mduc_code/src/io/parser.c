#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
/* #include <cstdarg.h> */
#include "parser.h"
#include "utils.h"
#include "display.h"
#include "fileIO.h"

parseStruct * parseFile(char *fileName)
{
  int i,j,n,iC;
  int multiLine, isFlamelet, isQuote;
  if (MDUC_PRINT_PROC)
    fprintf(stderr, "Parsing input file: %s\n", fileName);

  /* strings to store the tag for each line and comments (if any) */
  char tag[64], line[256], tmpL[256], *pch, *pq;

  /* initalize the parser struct */
  parseStruct *p = (parseStruct*) calloc(1, sizeof(parseStruct));
  p->nIn  = 0;
  p->comment = (char*) calloc(8, sizeof(char));
  p->delims  = (char*) calloc(8, sizeof(char));
  strcpy(p->comment, "#!");
  strcpy(p->delims,  " :=+\n");
  p->n    = (int*)    calloc(1, sizeof(int));
  p->tags = (char**)  calloc(1, sizeof(char**));
  p->vals = (char***) calloc(1, sizeof(char***));

  FILE *pF = fopen(fileName, "r"); /* open input file */

  if (!pF) {
    fprintf(stderr, "Warning: cannot open input file %s\n", fileName);
    return p;
  }
  
  /* ----- loop until end of file is reached  ----- */
  while (!feof(pF) && fgets(line,256,pF)) { 
    /* first check to see if end of input is desired */
    if (!strcmp(line, "ENDINPUT"))
      break;
    /* check to see if it is a comment */
    stripString(line);
    if (strspn(line, p->comment) > 0 || !strcmp(line, "\n"))
      continue;
    
    /* check whether line is empty */
    if (emptyString(line)) continue;
    if (!strcmp(line,"body")) continue;

    /* check whether input is across multiple lines */
    multiLine  = countChars(line, '{');
    multiLine -= countChars(line, '}');

    /* reallocate parser arrays */
    p->nIn++;
    i = p->nIn - 1;
    p->n    = (int*)    realloc(p->n,    p->nIn*sizeof(int));
    p->tags = (char**)  realloc(p->tags, p->nIn*sizeof(char**));
    p->vals = (char***) realloc(p->vals, p->nIn*sizeof(char***));
    p->vals[i] = (char**) calloc(1, sizeof(char**));

    /* split the line into parts */
    /* get the tag */
    pch = strtok(line, p->delims);
    p->tags[i] = (char*) calloc(strlen(pch)+1, sizeof(char));
    strcpy(p->tags[i], pch);

    /* backward compatability for flamelets */
    isFlamelet = 0;
    if (!strcmp(p->tags[i], "FLAMELET"))
      isFlamelet = 1;

    /* start getting the values */
    p->n[i] = 0;
    pch = strtok(NULL, p->delims);
    if (pch == NULL && !isFlamelet) {
      fprintf(stderr, "Warning: %s has no inputs\n", p->tags[i]);
      continue;
    }

    do {
      if (strpbrk(pch, p->comment)) break;
      isQuote = countChars(pch, '"');
      p->n[i]++;
      j = p->n[i] - 1;
      n = strlen(pch)+1;
      p->vals[i] = (char**) realloc(p->vals[i], p->n[i]*sizeof(char**));
      p->vals[i][j] = (char*) calloc(n, sizeof(char));
      strcpy(p->vals[i][j], pch);
      if (isQuote > 0) {
        /* locate the next quote */
        pch = strtok(NULL, "\"");
        n += strlen(pch)+4;
        strcpy(tmpL, p->vals[i][j]);
        p->vals[i][j] = (char*) realloc(p->vals[i][j],n*sizeof(char));
        sprintf(p->vals[i][j], "%s %s\0", &tmpL[1], pch);
      }
      pch = strtok(NULL, p->delims);
    } while (pch != NULL);
    
    while ((multiLine > 0 || isFlamelet > 0) &&
           !feof(pF) && fgets(line,256,pF)) {
      if (emptyString(line)) continue;
      multiLine += countChars(line, '{');
      multiLine -= countChars(line, '}');

      pch = strtok(line, p->delims);
      while (pch != NULL) {
        if (strpbrk(pch, p->comment)) break;
        if (!strcmp(pch, "ENDFLAMELET"))
          isFlamelet = 0;
        isQuote = countChars(pch, '"');

        p->n[i]++;
        j = p->n[i] - 1;
        n = strlen(pch)+1;
        p->vals[i] = (char**) realloc(p->vals[i], p->n[i]*sizeof(char**));
        p->vals[i][j] = (char*) calloc(n, sizeof(char));
        strcpy(p->vals[i][j], pch);
        if (isQuote > 0) {
          /* locate the next quote */
          pch = strtok(NULL, "\"");
          n += strlen(pch)+4;
          strcpy(tmpL, p->vals[i][j]);
          p->vals[i][j] = (char*) realloc(p->vals[i][j],n*sizeof(char));
          sprintf(p->vals[i][j], "%s %s\0", &tmpL[1], pch);
        }
        pch = strtok(NULL, p->delims);
      };
    }

    /* read in until end delimiter is found */
  };

  /* do some cleaning up */
  for (i=0; i<p->nIn; i++) {
    /* remove leading braces */
    if (!strcmp(p->vals[i][0], "{")) { /* stand-alone */
      for (j=1; j<p->n[i]; j++)
        strcpy(p->vals[i][j-1], p->vals[i][j]);
      p->n[i]--;
      free(p->vals[i][p->n[i]]);
      p->vals[i] = (char**) realloc(p->vals[i], p->n[i]*sizeof(char**));
    }
    else if (p->vals[i][0][0] == '{') { /* part of first string */
      n = strlen(p->vals[i][0]);
      for (j=1; j<n; j++)
        p->vals[i][0][j-1] = p->vals[i][0][j];
      p->vals[i][0][n-1] = '\0';
    }

    /* remove and trailing braces */
    if (!strcmp(p->vals[i][p->n[i]-1], "}") ||
        !strcmp(p->vals[i][p->n[i]-1], "ENDFLAMELET")) { /* stand-alone */
      p->n[i]--;
      free(p->vals[i][p->n[i]]);
      p->vals[i] = (char**) realloc(p->vals[i], p->n[i]*sizeof(char**));
    }
    else if (strpbrk(p->vals[i][p->n[i]-1], "}")) {
      n = strlen(p->vals[i][p->n[i]-1]);
      p->vals[i][p->n[i]-1][n-1] = '\0';
    }


    /* TODO: can't handle braces in middle of word */
  }

  fclose(pF);
  return p;
}
/* -------------------------------------------------------------------------- */

int paramGetIndex(char *name, parseStruct *p)
{
  int i;
  for (i=0; i<p->nIn; i++)
    if (!strcmp(name, p->tags[i]))
      return i;
  return -1;
}
/* -------------------------------------------------------------------------- */

int paramGetInt(char *name, int def, parseStruct *p)
{
  int i;
  for (i=0; i<p->nIn; i++)
    if (!strcmp(name, p->tags[i]))
      return atoi(p->vals[i][0]);
  return def;
}
/* -------------------------------------------------------------------------- */

double paramGetReal(char *name, double def, parseStruct *p)
{
  int i;
  for (i=0; i<p->nIn; i++)
    if (!strcmp(name, p->tags[i]))
      return atof(p->vals[i][0]);
  return def;
}
/* -------------------------------------------------------------------------- */

void paramGetString(char *val, char *name, char *def, parseStruct *p)
{
  int i,len;
  char buffer[256], *pch = NULL;
  
  len = strlen(name);
  strcpy(buffer,name);
  pch = strchr(buffer,'*');
  if (pch != NULL)
    len = strlen(buffer) - strlen(pch);

  for (i=0; i<p->nIn; i++) {
    if (!strncmp(name, p->tags[i],len)) {
      strcpy(val, p->vals[i][0]);
      return;
    }
  }
  strcpy(val, def);
}
/* -------------------------------------------------------------------------- */

int paramGetStringArray(char ***val, char *name, parseStruct *p)
{
  int i,j,len,nn;
  char buffer[256], *pch = NULL;

  /* check whether to recognize wildcard after specified length */
  len = strlen(name);
  strcpy(buffer,name);
  pch = strchr(buffer,'*');
  if (pch != NULL)
    len = strlen(buffer) - strlen(pch);

  for (i=0; i<p->nIn; i++) {
    if (!strncmp(name, p->tags[i], len)) {
      if (val[0] == NULL)
        val[0] = (char**) calloc(p->n[i], sizeof(char*));
      for (j=0; j<p->n[i]; j++) {
        nn = strlen(p->vals[i][j])+1;
        val[0][j] = (char*) calloc(nn, sizeof(char));
        strcpy(val[0][j], p->vals[i][j]);
      }
      return p->n[i];
    }
  }
  return 0;
}
/* -------------------------------------------------------------------------- */

int paramGetIntArray(int **out, char *name, parseStruct *p)
{
  int i,nA;
  char **arr = NULL;
  nA = paramGetStringArray(&arr,name,p);
  if (nA == 0) {
    out[0] = NULL;
    return 0;
  }
  out[0] = (int*) calloc(nA, sizeof(int));
  for (i=0; i<nA; i++) {
    out[0][i] = atoi(arr[i]);
    free(arr[i]);
  }
  free(arr);
  return nA;
}
/* -------------------------------------------------------------------------- */

int paramGetRealArray(double **out, char *name, parseStruct *p)
{
  int i,nA;
  char **arr = NULL;
  nA = paramGetStringArray(&arr,name,p);
  if (nA == 0) {
    out[0] = NULL;
    return 0;
  }
  out[0] = (double*) calloc(nA, sizeof(double));
  for (i=0; i<nA; i++) {
    out[0][i] = atof(arr[i]);
    free(arr[i]);
  }
  free(arr);
  return nA;
}
/* -------------------------------------------------------------------------- */

int sparamGetIndex(int n, char *name, char **vals)
{
  int i;
  for (i=0; i<n; i++)
    if (!strcmp(name, vals[i]))
      return i;
  return -1;
}
/* -------------------------------------------------------------------------- */

int sparamGetInt(int n, char *name, int def, char **vals)
{
  int i;
  for (i=0; i<n; i++)
    if (!strcmp(name, vals[i]))
      return atoi(vals[i+1]);
  return def;
}
/* -------------------------------------------------------------------------- */

double sparamGetReal(int n, char *name, double def, char **vals)
{
  int i;
  for (i=0; i<n; i++)
    if (!strcmp(name, vals[i]))
      return atof(vals[i+1]);
  return def;
}
/* -------------------------------------------------------------------------- */

int sparamGetIntArray(int n, char *name, int **out, char **vals)
{
  int i,nA;
  char **arr = NULL;
  nA = sparamGetStringArray(n,name,&arr,vals);
  if (nA == 0) {
    out[0] = NULL;
    return 0;
  }
  out[0] = (int*) calloc(nA, sizeof(int));
  for (i=0; i<nA; i++) {
    out[0][i] = atoi(arr[i]);
    free(arr[i]);
  }
  free(arr);
  return nA;
}
/* -------------------------------------------------------------------------- */

int sparamGetRealArray(int n, char *name, double **out, char **vals)
{
  int i,nA;
  char **arr = NULL;
  nA = sparamGetStringArray(n,name,&arr,vals);
  if (nA == 0) {
    out[0] = NULL;
    return 0;
  }
  out[0] = (double*) calloc(nA, sizeof(double));
  for (i=0; i<nA; i++) {
    out[0][i] = atof(arr[i]);
    free(arr[i]);
  }
  free(arr);
  return nA;
}
/* -------------------------------------------------------------------------- */

int sparamGetStringArray(int n, char *name, char ***out, char **vals)
{
  int i,j,nA = 0,iS;
  char aTmp[32];
  iS = sparamGetIndex(n,name,vals);
  if (iS < 0) {
    out[0] = NULL;
    return 0;
  }
  out[0] = (char**) calloc(1, sizeof(char*));
  i = iS+1;
  /* assume single value if brace not found */
  if (strchr(vals[i],'{') == NULL) {
    out[0][0] = (char*) calloc(64, sizeof(char));
    strcpy(out[0][0],vals[i]);
    return 1;
  }
  /* if single value in braces without spaces */
  if (strchr(vals[i],'{') && strchr(vals[i],'}')) {
    out[0][0] = (char*) calloc(64, sizeof(char));
    strncpy(out[0][0],&vals[i][1],strlen(vals[i])-2);
    return 1;
  }
  int done = 0;
  while (strcmp(vals[i],"}")) {
    if (!strcmp(vals[i],"{")) {
      i++;
      continue;
    }
    nA++;
    out[0] = (char**) realloc(out[0], nA*sizeof(char*));
    out[0][nA-1] = (char*) calloc(64, sizeof(char));

    if (vals[i][0] == '{')
      strcpy(out[0][nA-1],&vals[i][1]);
    else if (vals[i][strlen(vals[i])-1] == '}') {
      strncpy(out[0][nA-1],vals[i],strlen(vals[i])-1);
      done = 1;
    }
    else
      strcpy(out[0][nA-1],vals[i]);

    if (done) break;
    i++;
  }
  return nA;
}
/* -------------------------------------------------------------------------- */

int parsePrint(parseStruct *p)
{
  int i,j;
  fprintf(stderr, "Parser values:\n");
  fprintf(stderr, "  Number of inputs: %i\n", p->nIn);
  for (i=0; i<p->nIn; i++) {
    fprintf(stderr, " %3i: %s -> (%2i) ", i, p->tags[i], p->n[i]);
    for (j=0; j<p->n[i]; j++) {
      fprintf(stderr, "%s . ", p->vals[i][j]);
    }
    fprintf(stderr, " |\n");
  }
  return 0;
}
/* -------------------------------------------------------------------------- */

int sreadInt(char *str, int len, int n, int *out)
{
  int i;
  char buff[128];
  for (i=0; i<n; i++) {
    strncpy(buff,&str[i*len],len);
    buff[len] = '\0';
    out[i] = atoi(buff);
  }
  return 0;
}
/* -------------------------------------------------------------------------- */

int sreadReal(char *str, int len, int n, double *out)
{
  int i;
  char buff[128];
  for (i=0; i<n; i++) {
    strncpy(buff,&str[i*len],len);
    buff[len] = '\0';
    out[i] = atof(buff);
  }
  return 0;
}
/* -------------------------------------------------------------------------- */

int sreadString(char *str, int len, char *out)
{
  strncpy(out,str,len);
  out[len] = '\0';
  return 0;
}
/* -------------------------------------------------------------------------- */

int sreadStringArray(char *str, int len, int n, char **out)
{
  int i;
  char buff[128];
  for (i=0; i<n; i++) {
    strncpy(buff,&str[i*len],len);
    buff[len] = '\0';
    strcpy(out[i],buff);
  }
  return 0;
}
/* -------------------------------------------------------------------------- */

int readInteger(FILE *pFile)
{
  int i;
  char tag[100];
  fpos_t pos;
  
  fgetpos(pFile, &pos);  // get the current position of the buffer
  // check to see if an = is present
  fscanf(pFile, "%s", tag);
  if (strcmp(tag, "="))
    fsetpos(pFile, &pos); // reset stream to before value
  
  fscanf(pFile, "%i\n", &i);
  
  return i;
}
/* -------------------------------------------------------------------------- */

double readReal(FILE *pFile)
{
  long double r;
  char tag[100];
  fpos_t pos;
  
  fgetpos(pFile, &pos);  // get the current position of the buffer
  // check to see if an = is present
  fscanf(pFile, "%s", tag);
  if (strcmp(tag, "="))
    fsetpos(pFile, &pos); // reset stream to before value
  
  fscanf(pFile, "%Lf\n", &r);

  /* ----- check for units ----- */
  fgetpos(pFile, &pos);  // get the current position of the buffer
  fscanf(pFile, "%s", tag);
  if (tag[0] != '[')
    fsetpos(pFile, &pos); // reset stream to before value

  return (double)(r);
}
/* -------------------------------------------------------------------------- */

void readString(FILE *pFile, char *var)
{
  int  len;
  char s_tmp[200], line[200], quote[200], *tmp;
  
  fscanf(pFile, "%s\n", s_tmp);
  
  // account for an equals sign if present
  if (!strcmp(s_tmp, "="))
    fscanf(pFile, "%s\n", s_tmp);
  
  // if the string is in quotes, read entire string and strip quotes
  if (s_tmp[0] == '\"') {
    if (s_tmp[1] == '\"') {
      strcpy(s_tmp, ""); /* empty string */
    }
    else {
      fgets(line, 200, pFile);    // read in the entire line
      len = strlen(line);
      tmp = strtok(s_tmp, "\"");  // get rid of leading quote
      strcpy(quote, tmp);
      strcat(quote, " ");
      strncat(quote, line, len-2); // add line without trailing quote
      strcpy(s_tmp, quote);
    }
  }
  
  strcpy(var,s_tmp);
}
/* -------------------------------------------------------------------------- */

void readRealArray(FILE *pFile, int n, double **var)
{
  int i;
  char tag[100];
  fpos_t pos;
  
  // allocate the appropriate memory
  if (var[0] == NULL)
    var[0] = (double*) calloc(n, sizeof(double));
  
  fgetpos(pFile, &pos);
  // check to see if units are present
  fscanf(pFile, "%s", tag);
  if (tag[0] != '[')
    if ( (isdigit(tag[0])) | (isdigit(tag[1])) ) 
      fsetpos(pFile, &pos);
  
  // read in the array
  for (i=0; i < n; i++)
    var[0][i] = readReal(pFile);
}
/* -------------------------------------------------------------------------- */

int readRealArraySize(FILE *pFile, double **var)
{
  int    i, n;
  char   name[25];
  fpos_t pos;
  
  // get the number of components
  n = getNumComponents(pFile);
  
  // allocate memory for names array
  var[0] = (double*) calloc(n, sizeof(double));

  // loop back over data to read in name and mass fraction
  for (i=0; i<n; i++) var[0][i] = readReal(pFile);
  
  fscanf(pFile, "%s\n", name); // read in }
  return n;
}
/* -------------------------------------------------------------------------- */

void readOptArray(FILE *pFile, int n, double var[])
{
  int i, nc;
  char tag[100];
  fpos_t pos;
  
  nc = getNumComponents(pFile);

  fgetpos(pFile, &pos);
  // check to see if units are present
  fscanf(pFile, "%s", tag);
  if (tag[0] != '{')
    if ( (isdigit(tag[0])) | (isdigit(tag[1])) ) 
      fsetpos(pFile, &pos);
  
  if (nc == 1) // only one variable, set all elements to single value
    {
      var[0] = readReal(pFile);
      for (i=1; i<n; i++)
        var[i] = var[0];
    }
  else if (nc <= n) // read in the array
    {
      for (i=0; i < n; i++)
        var[i] = readReal(pFile);
    }
  else
    {
      fprintf(stderr, "Error: options array exceeds bounds.\n");
      exit(2);
    }

  fscanf(pFile, "%s\n", tag); // read in }
}
/* -------------------------------------------------------------------------- */

void readIntArray(FILE *pFile, int n, int var[])
{
  int i;
  char tag[100];
  fpos_t pos;
  
  fgetpos(pFile, &pos);
  // check to see if units are present
  fscanf(pFile, "%s", tag);
  if (tag[0] != '{')
    if ( (isdigit(tag[0])) | (isdigit(tag[1])) ) 
      fsetpos(pFile, &pos);
  
  // read in the array
  for (i=0; i < n; i++)
    var[i] = readInteger(pFile);

  fscanf(pFile, "%s\n", tag); // read in }
}
/* -------------------------------------------------------------------------- */

int getNumBC(FILE *pFile, int *fTemp, int *fH)
{
  int    nBC;
  double r_tmp;
  char   name[100];
  fpos_t pos;
  
  fscanf(pFile, "\n");   // go to next line
  fgetpos(pFile, &pos);  // get the current position of the buffer
  
  // get the number of components
  nBC=0;
  fscanf(pFile, "%s", name);
  do {
    if (!strcmp(name, "end"))
      break;
    if (!strcmp(name, "Temperature"))
      *fTemp = nBC;
    if (!strcmp(name, "Enthalpy"))
      *fH = nBC;
    fscanf(pFile, "%s", name); // read in "="
    r_tmp = readReal(pFile);
    nBC++;
    fscanf(pFile, "%s", name);
    if (name[0] == '[') // read in units if present
      fscanf(pFile, "%s", name);
  } while (strcmp(name,"end"));
  
  fsetpos(pFile, &pos); // reset stream to beginning of composition read
  
  return nBC;
}
/* -------------------------------------------------------------------------- */

char readBC(FILE *pFile, int *nComp, char ***ident, double *T, double **Y,
            double *H)
{
  int i, nBC, iComp, fTemp = -1, fH = -1;
  char name[50];
  
  // read in the begin statement ("begin", "{", etc.)
  fscanf(pFile, "%s\n", name);
  
  // get the number of components (including temperature)
  nBC = getNumBC(pFile, &fTemp, &fH);

  *nComp = nBC;
  
  if (fTemp >= 0) (*nComp)--;
  if (fH >= 0)    (*nComp)--;

  if (fTemp >= 0 && fH >=0)
    fprintf(stderr, "Warning: specifying both temperature and enthalpy BC!\n");
  
  // allocate memory for each (assume temperature has been read)
  Y[0]     = (double*) calloc(*nComp, sizeof(double));
  ident[0] = (char**)  calloc(*nComp, sizeof(char*));
  for (i = 0; i < *nComp; i++)
    ident[0][i] = (char*) calloc(50, sizeof(char));
  
  // loop back over data to read in name and mass fraction
  iComp = 0;
  for (i = 0; i < nBC; i++)
    {
      readString(pFile, ident[0][iComp]);
      if (ident[0][iComp][0] == '[')
        readString(pFile, ident[0][iComp]);
      fscanf(pFile, "%s", name); // read in equals
      if (i == fTemp)
        *T = readReal(pFile);
      else if (i == fH)
        *H = readReal(pFile);
      else
        {
          //stripHeader(ident[0][iComp]); // this should only be if from flamelet
          Y[0][iComp] = readReal(pFile);
          iComp++;
        }
    }
  fscanf(pFile, "%s\n", name); // read in "end"

  if (fTemp >= 0)
    return 'T';
  else if (fH >= 0)
    return 'H';
  else
    return 'U';
}
/* -------------------------------------------------------------------------- */

int getNumFuelBC(FILE *pFile, int *fTemp, int *fH, int *fHvap, int *fCpL,
                 int *fTref, int *fSat)
{
  int    nBC;
  double r_tmp;
  char   name[100];
  fpos_t pos;
  
  fscanf(pFile, "\n");   // go to next line
  fgetpos(pFile, &pos);  // get the current position of the buffer
  
  // get the number of components
  nBC=0;
  fscanf(pFile, "%s", name);
  do {
    if (!strcmp(name, "end"))
      break;
    if (!strcmp(name, "Temperature"))
      *fTemp = nBC;
    if (!strcmp(name, "Enthalpy"))
      *fH = nBC;
    if (!strcmp(name, "Latent-heat"))
      *fHvap = nBC;
    if (!strcmp(name, "Liquid-cp"))
      *fCpL = nBC;
    if (!strcmp(name, "Temp-ref"))
      *fTref = nBC;
    if (!strcmp(name, "sat-data"))
      *fSat = nBC;
    fscanf(pFile, "%s", name); // read in "="
    r_tmp = readReal(pFile);
    nBC++;
    fscanf(pFile, "%s", name);
    if (name[0] == '[') // read in units if present
      fscanf(pFile, "%s", name);
  } while (strcmp(name,"end"));
  
  fsetpos(pFile, &pos); // reset stream to beginning of composition read
  
  return nBC;
}
/* -------------------------------------------------------------------------- */

char readFuelBC(FILE *pFile, int *nComp, char ***ident, double *T, double **Y,
                double *H, double *Hvap, double *cpL, double *Tref, char *data)
{
  int i, nBC, iComp, fTemp = -1, fH = -1, fTref = -1, fHvap = -1, fCpL = -1;
  int fSat = -1;
  char name[50];
  
  // read in the begin statement ("begin", "{", etc.)
  fscanf(pFile, "%s\n", name);
  
  // get the number of components (including temperature)
  nBC = getNumFuelBC(pFile, &fTemp, &fH, &fHvap, &fCpL, &fTref, &fSat);

  *nComp = nBC; // set number of components read

  // decrease the number of components for each non-species component read
  if (fTemp >= 0) (*nComp)--;
  if (fH >= 0)    (*nComp)--;
  if (fTref >= 0) (*nComp)--;
  if (fHvap >= 0) (*nComp)--;
  if (fCpL >= 0)  (*nComp)--;
  if (fSat >= 0)  (*nComp)--;

  if (fTemp >= 0 && fH >=0)
    fprintf(stderr, "Warning: specifying both temperature and enthalpy BC!\n");
  
  // allocate memory for each (assume temperature has been read)
  Y[0]     = (double*) calloc(*nComp, sizeof(double));
  ident[0] = (char**)  calloc(*nComp, sizeof(char*));
  for (i = 0; i < *nComp; i++)
    ident[0][i] = (char*) calloc(50, sizeof(char));
  
  // loop back over data to read in name and mass fraction
  iComp = 0;
  for (i = 0; i < nBC; i++)
    {
      readString(pFile, ident[0][iComp]);
      if (ident[0][iComp][0] == '[')
        readString(pFile, ident[0][iComp]);
      fscanf(pFile, "%s", name); // read in equals
      if (i == fTemp)
        *T = readReal(pFile);
      else if (i == fH)
        *H = readReal(pFile);
      else if (i == fTref)
        *Tref = readReal(pFile);
      else if (i == fHvap)
        *Hvap = readReal(pFile);
      else if (i == fCpL)
        *cpL = readReal(pFile);
      else if (i == fSat)
        readString(pFile, data);
      else
        {
          Y[0][iComp] = readReal(pFile);
          iComp++;
        }
    }
  fscanf(pFile, "%s\n", name); // read in "end"

  if (fCpL >= 0)
    return 'L';
  else if (fSat >= 0)
    return 'D';
  else if (fTemp >= 0)
    return 'T';
  else if (fH >= 0)
    return 'H';
  else
    return 'U';
}
/* -------------------------------------------------------------------------- */

int readNames(FILE *pFile, char ***ident)
{
  int    i, n;
  char   name[25];
  fpos_t pos;
  
  // get the number of components
  n = getNumComponents(pFile);
  
  // allocate memory for names array
  ident[0] = (char**) calloc(n, sizeof(char*));
  for (i=0; i<n; i++)
    ident[0][i] = (char*) calloc(25, sizeof(char));
  
  // loop back over data to read in name and mass fraction
  for (i = 0; i < n; i++)
    readString(pFile, ident[0][i]);
  
  fscanf(pFile, "%s\n", name); // read in }
  return n;
}
/* -------------------------------------------------------------------------- */

int getNumComponents(FILE *pFile)
{
  int    numComp;
  float  r_temp;
  char   name[200];
  fpos_t pos;

  fgetpos(pFile, &pos);  // get the current position of the buffer

  // check for leading brace {
  fscanf(pFile, "%s\n", name);
  if (strcmp(name, "{")) 
    {
      fsetpos(pFile, &pos);
      return 1;
    }
  fgetpos(pFile, &pos);
  
  // get the number of components
  numComp=0;
  do {
    fscanf(pFile, "%s\n", name);
    numComp++;
  } while (strcmp(name,"}"));
  
  fsetpos(pFile, &pos); // reset stream to beginning of composition read
  
  return --numComp;     // do loop always counts one too many
}
/* -------------------------------------------------------------------------- */

/*
void readComposition(FILE *pFile, int *nComp, char ***ident, double **y_species)
{
  int    i;
  char   name[25];
  
  // get the number of components
  *nComp = getNumComponents(pFile);
  
  // allocate memory for each
  y_species[0] = (double*) calloc(*nComp, sizeof(double));
  ident[0]     = (char**)  calloc(*nComp, sizeof(char*));
  for (i = 0; i < *nComp; i++)
    ident[0][i] = (char*) calloc(25, sizeof(char));
  
  // loop back over data to read in name and mass fraction
  for (i = 0; i < *nComp; i++)
    {
      readString(pFile, ident[0][i]);
      y_species[0][i] = readReal(pFile);
    }
  fscanf(pFile, "%s\n", name); // read in }
} */
/* -------------------------------------------------------------------------- */

void destroyParser(parseStruct *p)
{
  int i,j;
  if(p->comment) free(p->comment);
  if(p->delims)  free(p->delims);
  for (i=0; i<p->nIn; i++) {
    if (p->tags[i]) free(p->tags[i]);
    for (j=0; j<p->n[i]; j++)
      if (p->vals[i][j]) free(p->vals[i][j]);
    if (p->vals[i]) free(p->vals[i]);
  }
  if (p->n)    free(p->n);
  if (p->tags) free(p->tags);
  if (p->vals) free(p->vals);
}
/* -------------------------------------------------------------------------- */
