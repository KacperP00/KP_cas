#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "chemTable.h"
#include "chemMech.h"
#include "utils.h"

void *loadChemTable(char *chemTabFile, int nDim, char *type)
{
  int i,j,k,l,p,q,lenData,tType,iSp;
  long double r_tmp;

  fprintf(stderr, "  Loading %s %iD chemtable: %s\n", type, nDim, chemTabFile);
  FILE *chemFile = NULL;
  chemFile = fopen(chemTabFile, "rb");

  if (!chemFile) {
    fprintf(stderr, "Error: could not locate chemtable file %s\n", chemTabFile);
    exit(2);
  }

  /* allocate the chemtable memory */
  chemTab *cT = (chemTab*) calloc(1, sizeof(chemTab));

  /* check the table type */
  if (!strcmp(type, "non-premixed"))
    tType = 0;
  else if (!strcmp(type, "premixed"))
    tType = 1;
  else if (!strcmp(type, "pressure"))
    tType = 9;

  if (nDim > 0)
    cT->nDim = nDim;
  else {
    cT->nDim = 4;   /* ----- TODO add 4th dimension ----- default for now */
    if (tType == 9) cT->nDim = 3; 
  }
  /* ----- read in the table header ----- */
  cT->n = (int*) calloc(cT->nDim, sizeof(int));
  int *n = cT->n;

  /* get length of each co-ordinate */
  fread(cT->n, sizeof(int), cT->nDim, chemFile);

  /* get number of variables */
  fread(&cT->nVar, sizeof(int), 1, chemFile);

  /* ----- allocate the memory ----- */
  /* get the data sizes */
  lenData = 1;
  for (i=0; i<cT->nDim; i++) lenData *= n[i];
  cT->nData = lenData*cT->nVar;

  /* allocate the core data array */
  cT->data = (double*) calloc(cT->nData, sizeof(double));

  /* allocate the independent co-ordinate arrays/masks */
  cT->x = (double**) calloc(cT->nDim, sizeof(double*));
  for (i=0; i<cT->nDim; i++)
    cT->x[i] = (double*) calloc(cT->n[i], sizeof(double));
  cT->mask = (int*) calloc(lenData, sizeof(int));

  /* allocate name arrays */
  cT->names = (char**) calloc(cT->nVar, sizeof(char*));
  for (i=0; i<cT->nVar; i++)
    cT->names[i] = (char*) calloc(64, sizeof(char));

  /* ----- read in data ----- */
  /* read in the mapping variables */
  for (i=0; i<cT->nDim; i++)
    fread(cT->x[i], sizeof(double), cT->n[i], chemFile);

  /* read in the mask */
  fread(cT->mask, sizeof(int), lenData, chemFile);

  /* read in combustion model */
  fread(&cT->model, sizeof(char), 64, chemFile);

  if (tType == 9) {
    fread(&cT->Pref, sizeof(double), 1, chemFile);
    fread(&cT->Tref, sizeof(double), 1, chemFile);
  }

  /* read in variable names */
  cT->nSpecies = 0;
  cT->species  = (char**) calloc(1, sizeof(char*));
  cT->iSp      = (int*)   calloc(1, sizeof(int));
  for (i=0; i<cT->nVar; i++) {
    fread(cT->names[i], sizeof(char), 64, chemFile);
    cT->names[i][63] = '\0';
    stripString(cT->names[i]);
    
    /* check whether the variable is a species */
    if (!strncmp(cT->names[i], "Y_", 2)) {
      iSp = cT->nSpecies;
      cT->nSpecies++;
      cT->species = (char**) realloc(cT->species, cT->nSpecies*sizeof(char*));
      cT->species[iSp] = (char*) calloc(64, sizeof(char));
      cT->iSp     = (int*) realloc(cT->iSp, cT->nSpecies*sizeof(int));
      strcpy(cT->species[iSp], &cT->names[i][2]);
      cT->iSp[iSp] = i;
      if (i > 0 && cT->iSp[iSp-1] != i-1) {
        fprintf(stderr, "Warning: species are non-contiguous in chemtable\n");
        fprintf(stderr, "         Currently not supported\n");
      }
    }
  }

  /* read the mapped variables */
  cT->iRho = cT->iProg = cT->iHR = -1;
  for (q=0; q<cT->nVar; q++) {
    fread(&cT->data[q*lenData], sizeof(double), lenData, chemFile);
    
    if      (!strncmp(cT->names[q], "RHO", 3))  cT->iRho  = q;
    else if (!strncmp(cT->names[q], "PROG", 4)) cT->iProg = q;
    else if (!strncmp(cT->names[q], "HR", 2))   cT->iHR   = q;
  }

  /* ----- post-process the data ----- */
  /* co-ordinates */
  cT->xmin = (double*) calloc(cT->nDim, sizeof(double));
  cT->xmax = (double*) calloc(cT->nDim, sizeof(double));
  for (i=0; i<cT->nDim; i++) {
    cT->xmin[i] =  1.0e20;
    cT->xmax[i] = -1.0e20;
    for (j=0; j<cT->n[i]; j++) {
      cT->xmin[i] = fmin(cT->xmin[i], cT->x[i][j]);
      cT->xmax[i] = fmax(cT->xmax[i], cT->x[i][j]);
    }
  }

  /* variables */
  cT->datMin = (double*) calloc(cT->nVar, sizeof(double));
  cT->datMax = (double*) calloc(cT->nVar, sizeof(double));
  for (i=0; i<cT->nVar; i++) {
    cT->datMin[i] =  1.0e20;
    cT->datMax[i] = -1.0e20;
    for (j=0; j<lenData; j++) {
      cT->datMin[i] = fmin(cT->datMin[i], cT->data[i*lenData+j]);
      cT->datMax[i] = fmax(cT->datMax[i], cT->data[i*lenData+j]);
    }
  }

  if (tType == 1) { /* compute min/max progress variable */
    int nProg = 1;
    for (i=0; i<cT->nDim-1; i++) nProg *= n[i];
    cT->progMin = (double*) calloc(nProg, sizeof(double));
    cT->progMax = (double*) calloc(nProg, sizeof(double));
    for (i=0; i<nProg; i++) {
      cT->progMin[i] =  1.0e20;
      cT->progMax[i] = -1.0e20;
    }

    if (cT->iProg < 0) {
      fprintf(stderr, "Error: progress variable not found in premixed table\n");
      exit(2);
    }
    
    /* TODO: this assumes progress variable is always last co-ordinate */
    if (cT->nDim == 4) {
      for (i=0; i<n[0]; i++) {
        for (j=0; j<n[1]; j++) {
          for (k=0; k<n[2]; k++) {
            for (l=0; l<n[3] ; l++) {
              cT->progMin[cti3(i,j,k,0,n)] =
                fmin(cT->progMin[cti3(i,j,k,0,n)],
                     cT->data[cti4(i,j,k,l,cT->iProg,n)]);
              cT->progMax[cti3(i,j,k,0,n)] =
                fmax(cT->progMax[cti3(i,j,k,0,n)],
                     cT->data[cti4(i,j,k,l,cT->iProg,n)]);
            }
          }
        }
      }
    }
    else {
      for (i=0; i<n[0]; i++) {
        for (j=0; j<n[1]; j++) {
          for (k=0; k<n[2]; k++) {
            cT->progMin[cti2(i,j,0,n)] =
              fmin(cT->progMin[cti2(i,j,0,n)],
                   cT->data[cti3(i,j,k,cT->iProg,n)]);
            cT->progMax[cti2(i,j,0,n)] =
              fmax(cT->progMax[cti2(i,j,0,n)],
                   cT->data[cti3(i,j,k,cT->iProg,n)]);
          }
        }
      }
    }
  } 
  fclose(chemFile);
  return (void*)cT;

}
/* ------------------------------------------------------------------------- */

double lookupChemTableIndex(int *index, double *x, double *vars, void *tab)
{
   /* dereference struct */
   chemTab *cT = (chemTab*) tab;
   if (cT->nDim == 4) {
     return chemTableInterp4DIndex(index, x, vars, tab);
   }
   else if (cT->nDim == 3) {
     return chemTableInterp3DIndex(index, x, vars, tab);
   }

   fprintf(stderr,"Error: interpolation of %iD table not available\n",cT->nDim);
   return 0.0;
}
/* ------------------------------------------------------------------------- */

double lookupChemTable(char *name, double *x, double *vars, void *tab)
{
  /* dereference struct */
  chemTab *cT = (chemTab*) tab;

  /* use specific interpolators for now */
  if (cT->nDim == 2) {
    return chemTableInterp2D(name, x, vars, tab);
  }
  else if (cT->nDim == 3) {
    return chemTableInterp3D(name, x, vars, tab);
  }
  else if (cT->nDim == 4) {
    return chemTableInterp4D(name, x, vars, tab);
  }

  fprintf(stderr,"Error: interpolation of %iD table not available\n",cT->nDim);
  return 0.0;
}
/* ------------------------------------------------------------------------- */

double lookupChemTablePT(char *name, double *x, double *vars, void *tab,
                         void *tabEx)
{
  /* assumes that x is length tab.nDim + 2
       P = x[tab.nDim]
       T = x[tab.nDim] 
     If either P or T are zero, only one expansion is needed */
  
  int i,iV;
  chemTab *cT   = (chemTab*) tab;
  chemTab *cTEx = (chemTab*) tabEx;
  int npar = 11; /* TODO: this should be in the table, not hardcoded */

  double P        = x[cT->nDim];
  double T        = x[cT->nDim+1];
  double Cprog    = x[3]; /* TODO: hardcoded */
  double Pstar    = log(P/cT->Pref);
  double PstarSq  = Pstar*Pstar;
  double Tstar    = T/cT->Tref;
  double TstarSq  = Tstar*Tstar;
  double logTstar = log(Tstar);
  double logTsq   = log(Tstar)*log(Tstar);
  double *xRef = (double*) calloc(cT->nDim, sizeof(double));
  double *xEx  = (double*) calloc(cTEx->nDim, sizeof(double));
  double *cPT  = (double*) calloc(cT->nVar*npar, sizeof(double));
  double CM,Ca,Cb,fP,fT,cMaxRef,cMinRef,V;
  int    all,nV,iV0,offset;

  /* reference coordinates */
  xEx[0] = x[0];
  xEx[1] = x[1];
  if (cTEx->nDim > 2) /* no variance for table at this point */
    xEx[2] = 0.0;

  /* look up the expansion coefficients */
  lookupChemTable(name,xEx,cPT,tab);

  /* copy the current co-ordinates to the reference */
  for (i=0; i<cT->nDim; i++) xRef[i] = x[i];

  all = chemTableCheckLookup(name,&nV,&iV0,&offset,cT);
  if (all < 0)
    fprintf(stderr, "ChemTable Warning: %s not found in table!\n", name);

  /* modify all the variables */
  for (iV=0; iV<nV; iV++) {
    fP = exp(cPT[i*npar+0]*PstarSq + cPT[i*npar+1]*Pstar);
    CM = fmax(cPT[i*npar+3]*(Tstar-1.0)+cPT[i*npar+4]*TstarSq+1.0-cPT[i*npar+4],
              1e-10);
    Ca = fmax(cPT[i*npar+5]*(Tstar-1.0)+1.0, 1e-10);
    Cb = fmax(cPT[i*npar+6]*(Tstar-1.0)+1.0, 1e-10);
    fT = exp(cPT[i*npar+8]*logTsq + cPT[i*npar+7]*logTstar);
    cMaxRef = cPT[i*npar+9];
    cMinRef = cPT[i*npar+10];

    /* transform the progress variable to the reference profile */
    if (Cprog - cMaxRef*CM < 0.0)
      xRef[3] = (Cprog - cMaxRef*CM)/Ca + cMaxRef;
    else
      xRef[3] = (Cprog - cMaxRef*CM)/Cb + cMaxRef;   

    /* look up the variable from the reference profile */
    V = lookupChemTable(cT->names[i],xRef,NULL,tab);

    /* scale the profile */
    V = fP*((V-cMinRef)*fT+cMinRef);

    if (!all) {
      free(xRef); free(xEx); free(cPT);
      return V;
    }
    vars[iV-offset] = V;
  }

  free(xRef); free(xEx); free(cPT);
  return V;
}
/* ------------------------------------------------------------------------- */

double chemTableVarMin(char *name, void *tab)
{
  /* dereference struct */
  chemTab *cT = (chemTab*) tab;
  return cT->datMin[getSpeciesIndex(cT->nVar, cT->names, name)];
}
/* ------------------------------------------------------------------------- */

double chemTableVarMax(char *name, void *tab)
{
  /* dereference struct */
  chemTab *cT = (chemTab*) tab;
  return cT->datMax[getSpeciesIndex(cT->nVar, cT->names, name)];
}
/* ------------------------------------------------------------------------- */

double chemTableProgMax(double *x, void *tab)
{
  int    i,j,l;
  double V,w11,w12,w21,w22,w31,w32;
  /* dereference struct */
  chemTab *cT = (chemTab*) tab;

  /* look up co-ordinate location indices */
  i = gridLocate(cT->n[0], cT->x[0], x[0]); 
  j = gridLocate(cT->n[1], cT->x[1], x[1]); 
  if (cT->nDim == 4) {  /* check dimension */
    /* ----- look up co-ordinate location indices -----  */
    l = gridLocate(cT->n[2], cT->x[2], x[2]);
    
    /* ----- compute the weights ----- */
    w11 = (cT->x[0][i+1] - x[0])/(cT->x[0][i+1] - cT->x[0][i]);
    w21 = (cT->x[1][j+1] - x[1])/(cT->x[1][j+1] - cT->x[1][j]);
    w31 = (cT->x[2][l+1] - x[2])/(cT->x[2][l+1] - cT->x[2][l]);
    w11 = fmin(fmax(w11, 0.0), 1.0);
    w21 = fmin(fmax(w21, 0.0), 1.0);
    w31 = fmin(fmax(w31, 0.0), 1.0);
    w12 = 1.0 - w11;
    w22 = 1.0 - w21;
    w32 = 1.0 - w31;
    
    /* ----- interpolate variable ----- */
    V = w31*( w21*( w11*cT->progMax[cti3(i,  j,  l,  0,cT->n)]
                   +w12*cT->progMax[cti3(i+1,j,  l,  0,cT->n)] )
             +w22*( w11*cT->progMax[cti3(i,  j+1,l,  0,cT->n)]
                   +w12*cT->progMax[cti3(i+1,j+1,l,  0,cT->n)] ) )
       +w32*( w21*( w11*cT->progMax[cti3(i,  j,  l+1,0,cT->n)]
                   +w12*cT->progMax[cti3(i+1,j,  l+1,0,cT->n)] )
             +w22*( w11*cT->progMax[cti3(i,  j+1,l+1,0,cT->n)]
                   +w12*cT->progMax[cti3(i+1,j+1,l+1,0,cT->n)] ) );
  }
  else if (cT->nDim == 3) { 
    /* ----- compute the weights ----- */
    w11 = (cT->x[0][i+1] - x[0])/(cT->x[0][i+1] - cT->x[0][i]);
    w21 = (cT->x[1][j+1] - x[1])/(cT->x[1][j+1] - cT->x[1][j]);
    w11 = fmin(fmax(w11, 0.0), 1.0);
    w21 = fmin(fmax(w21, 0.0), 1.0);
    w12 = 1.0 - w11;
    w22 = 1.0 - w21;
    
    /* ----- interpolate variable ----- */
    V = w21*( w11*cT->progMax[cti2(i,  j,  0,cT->n)]
             +w12*cT->progMax[cti2(i+1,j,  0,cT->n)] )
       +w22*( w11*cT->progMax[cti2(i,  j+1,0,cT->n)]
             +w12*cT->progMax[cti2(i+1,j+1,0,cT->n)] );
  }
  return V;
}
/* ------------------------------------------------------------------------- */

double chemTableProgMin(double *x, void *tab)
{
  int    i,j,l;
  double V,w11,w12,w21,w22,w31,w32;
  /* dereference struct */
  chemTab *cT = (chemTab*) tab;

   /* ----- look up co-ordinate location indices ----- */
  i = gridLocate(cT->n[0], cT->x[0], x[0]); 
  j = gridLocate(cT->n[1], cT->x[1], x[1]); 
  if (cT->nDim == 4) { /* check dimension */
    /* ----- look up co-ordinate location indices ----- */
    l = gridLocate(cT->n[2], cT->x[2], x[2]);
    
    /* ----- compute the weights ----- */
    w11 = (cT->x[0][i+1] - x[0])/(cT->x[0][i+1] - cT->x[0][i]);
    w21 = (cT->x[1][j+1] - x[1])/(cT->x[1][j+1] - cT->x[1][j]);
    w31 = (cT->x[2][l+1] - x[2])/(cT->x[2][l+1] - cT->x[2][l]);
    w11 = fmin(fmax(w11, 0.0), 1.0);
    w21 = fmin(fmax(w21, 0.0), 1.0);
    w31 = fmin(fmax(w31, 0.0), 1.0);
    w12 = 1.0 - w11;
    w22 = 1.0 - w21;
    w32 = 1.0 - w31;
    
    /* ----- interpolate variable ----- */
    V = w31*( w21*( w11*cT->progMin[cti3(i,  j,  l,  0,cT->n)]
                   +w12*cT->progMin[cti3(i+1,j,  l,  0,cT->n)] )
             +w22*( w11*cT->progMin[cti3(i,  j+1,l,  0,cT->n)]
                   +w12*cT->progMin[cti3(i+1,j+1,l,  0,cT->n)] ) )
       +w32*( w21*( w11*cT->progMin[cti3(i,  j,  l+1,0,cT->n)]
                   +w12*cT->progMin[cti3(i+1,j,  l+1,0,cT->n)] )
             +w22*( w11*cT->progMin[cti3(i,  j+1,l+1,0,cT->n)]
                   +w12*cT->progMin[cti3(i+1,j+1,l+1,0,cT->n)] ) );
  }
  else if (cT->nDim == 3) { 
    /* ----- compute the weights ----- */
    w11 = (cT->x[0][i+1] - x[0])/(cT->x[0][i+1] - cT->x[0][i]);
    w21 = (cT->x[1][j+1] - x[1])/(cT->x[1][j+1] - cT->x[1][j]);
    w11 = fmin(fmax(w11, 0.0), 1.0);
    w21 = fmin(fmax(w21, 0.0), 1.0);
    w12 = 1.0 - w11;
    w22 = 1.0 - w21;
    
    /* ----- interpolate variable ----- */
    V = w21*( w11*cT->progMin[cti2(i,  j,  0,cT->n)]
              +w12*cT->progMin[cti2(i+1,j,  0,cT->n)] )
      +w22*( w11*cT->progMin[cti2(i,  j+1,0,cT->n)]
             +w12*cT->progMin[cti2(i+1,j+1,0,cT->n)] );
  }
  return V;
}
/* ------------------------------------------------------------------------- */

double chemTableProgScaleP(double P, double *x, void *tab)
{
  /* dereference struct */
  chemTab *cT = (chemTab*) tab;

  double sc_prog, prog, ratio, tmp, f[3];

  /* store the last co-ordinate  */
  prog = x[2];

  /* look up the factors */
  x[2] = 1.0;
  f[1] = lookupChemTable("PROG",x,NULL,tab);
  x[2] = 2.0;
  f[2] = lookupChemTable("PROG",x,NULL,tab);

  /* scale the pressure */
  ratio = log(P/cT->Pref);
  tmp   = exp(f[1]*ratio*ratio + f[2]*ratio);

  /* restore the original co-ordinate */
  x[2] = prog;
  if (tmp > 0.0)
    x[2] /= tmp;

  return x[2];
}
/* ------------------------------------------------------------------------- */

int chemTableCheckLookup(char *name,int *nV,int *iV0,int *offset,chemTab *cT)
{
  if (!strcmp(name, "all")) { /* lookup all variables */
    *iV0 = 0;
    *nV  = cT->nVar;
    return 1;
  }
  else if (!strcmp(name, "species")) {
    /* TODO: assumes species are contiguous */
    *iV0    = cT->iSp[0];
    *nV     = cT->nSpecies;
    *offset = cT->iSp[0];
    return 2;
  }
  else { /* look up specific variable */
    /* TODO: rename this function? */
    *iV0 = getSpeciesIndex(cT->nVar,cT->names,name);
    *nV  = *iV0+1;
    /* check to see if variable exists */
    if (*iV0 < 0) return -1;
    return 0;
  }
  return -1;
}
/* ------------------------------------------------------------------------- */

double chemTableInterp2D(char *name, double *x, double *vars, chemTab *cT)
{
  int    i,j,iV,iV0,nV,all,offset = 0;
  double V;
  double w11,w12,w21,w22;

  if (!strcmp(name, "all")) { /* lookup all variables */
    all = 1;
    iV0 = 0;
    nV  = cT->nVar;
  }
  else if (!strcmp(name, "species")) {
    /* TODO: assumes species are contiguous */
    all = 2;
    iV0 = cT->iSp[0];
    nV  = cT->nSpecies;
    offset = cT->iSp[0];
  }
  else { /* look up specific variable */
    all = 0;
    /* TODO: rename this function? */
    iV0 = getSpeciesIndex(cT->nVar,cT->names,name);    
    nV  = iV0+1;
    if (iV0 < 0) {
      fprintf(stderr, "ChemTable Warning: %s not found in table!\n", name);
      return 0.0;
    }
  }
  
  /* ----- look up co-ordinate location indices ----- */
  i = gridLocate(cT->n[0], cT->x[0], x[0]); 
  j = gridLocate(cT->n[1], cT->x[1], x[1]); 
  
  /* ----- compute the weights ----- */
  w11 = (cT->x[0][i+1] - x[0])/(cT->x[0][i+1] - cT->x[0][i]);
  w21 = (cT->x[1][j+1] - x[1])/(cT->x[1][j+1] - cT->x[1][j]);
  w11 = fmin(fmax(w11, 0.0), 1.0);
  w21 = fmin(fmax(w21, 0.0), 1.0);
  w12 = 1.0 - w11;
  w22 = 1.0 - w21;
  
  for (iV=iV0; iV<nV; iV++) {
    /* ----- interpolate variable ----- */
    V = w21*( w11*cT->data[cti2(i,  j,  iV,cT->n)]
             +w12*cT->data[cti2(i+1,j,  iV,cT->n)] )
       +w22*( w11*cT->data[cti2(i,  j+1,iV,cT->n)]
             +w12*cT->data[cti2(i+1,j+1,iV,cT->n)] );

    if (!all) return V;
    vars[iV-offset] = V;
  } /* end loop over variables */
  return V;
}
/* ------------------------------------------------------------------------- */

double chemTableInterp3D(char *name, double *x, double *vars, chemTab *cT)
{
  int    i,j,k,iV,iV0,nV,all,offset = 0;
  double V;
  double w11,w12,w21,w22,w31,w32;

  if (!strcmp(name, "all")) { /* lookup all variables */
    all = 1;
    iV0 = 0;
    nV  = cT->nVar;
  }
  else if (!strcmp(name, "species")) {
    /* TODO: assumes species are contiguous */
    all = 2;
    iV0 = cT->iSp[0];
    nV  = cT->nSpecies;
    offset = cT->iSp[0];
  }
  else { /* look up specific variable */
    all = 0;
    /* TODO: rename this function? */
    iV0 = getSpeciesIndex(cT->nVar,cT->names,name);    
    nV  = iV0+1;
    if (iV0 < 0) {
      fprintf(stderr, "ChemTable Warning: %s not found in table!\n", name);
      return 0.0;
    }
  }
  
  /* ----- look up co-ordinate location indices ----- */
  i = gridLocate(cT->n[0], cT->x[0], x[0]); 
  j = gridLocate(cT->n[1], cT->x[1], x[1]); 
  k = gridLocate(cT->n[2], cT->x[2], x[2]); 
  
  /* ----- compute the weights ----- */
  w11 = (cT->x[0][i+1] - x[0])/(cT->x[0][i+1] - cT->x[0][i]);
  w21 = (cT->x[1][j+1] - x[1])/(cT->x[1][j+1] - cT->x[1][j]);
  w31 = (cT->x[2][k+1] - x[2])/(cT->x[2][k+1] - cT->x[2][k]);
  w11 = fmin(fmax(w11, 0.0), 1.0);
  w21 = fmin(fmax(w21, 0.0), 1.0);
  w31 = fmin(fmax(w31, 0.0), 1.0);
  w12 = 1.0 - w11;
  w22 = 1.0 - w21;
  w32 = 1.0 - w31;

  /* ----- interpolate variable ----- */
  for (iV=iV0; iV<nV; iV++) {
    V = w31*( w21*( w11*cT->data[cti3(i,  j,  k,  iV,cT->n)]
                   +w12*cT->data[cti3(i+1,j,  k,  iV,cT->n)] )
             +w22*( w11*cT->data[cti3(i,  j+1,k,  iV,cT->n)]
                   +w12*cT->data[cti3(i+1,j+1,k,  iV,cT->n)] ) )
       +w32*( w21*( w11*cT->data[cti3(i,  j,  k+1,iV,cT->n)]
                   +w12*cT->data[cti3(i+1,j,  k+1,iV,cT->n)] )
             +w22*( w11*cT->data[cti3(i,  j+1,k+1,iV,cT->n)]
                   +w12*cT->data[cti3(i+1,j+1,k+1,iV,cT->n)] ) );

    if (!all) return V;
    vars[iV-offset] = V;
  } /* end loop over variables */

  return V;
}
/* ------------------------------------------------------------------------- */

double chemTableInterp4D(char *name, double *x, double *vars, chemTab *cT)
{
  int    i,j,k,l,ig,jg,kg,lg,iV,iV0,nV,all,offset = 0;
  double V;
  double w11,w12,w21,w22,w31,w32,w41,w42,w,wtest;
  if (!strcmp(name, "all")) { /* lookup all variables */
    all = 1;
    iV0 = 0;
    nV  = cT->nVar;
  }
  else if (!strcmp(name, "species")) {
    /* TODO: assumes species are contiguous */
    all = 2;
    iV0 = cT->iSp[0];
    nV  = cT->nSpecies;
    offset = cT->iSp[0];
  }
  else { /* look up specific variable */
    all = 0;
    /* TODO: rename this function? */
    iV0 = getSpeciesIndex(cT->nVar,cT->names,name);    
    nV  = iV0+1;
    if (iV0 < 0) {
      fprintf(stderr, "ChemTable Warning: %s not found in table!\n", name);
      return 0.0;
    }
  }
  
  /* ----- look up co-ordinate location indices ----- */
  i = gridLocate(cT->n[0], cT->x[0], x[0]); 
  j = gridLocate(cT->n[1], cT->x[1], x[1]); 
  k = gridLocate(cT->n[2], cT->x[2], x[2]); 
  l = gridLocate(cT->n[3], cT->x[3], x[3]);

  /* ----- compute the weights ----- */
  w11 = fmin(fmax((cT->x[0][i+1] - x[0])/(cT->x[0][i+1] - cT->x[0][i]),0),1);
  w21 = fmin(fmax((cT->x[1][j+1] - x[1])/(cT->x[1][j+1] - cT->x[1][j]),0),1);
  w31 = fmin(fmax((cT->x[2][k+1] - x[2])/(cT->x[2][k+1] - cT->x[2][k]),0),1);
  w41 = fmin(fmax((cT->x[3][l+1] - x[3])/(cT->x[3][l+1] - cT->x[3][l]),0),1);
  w12 = 1.0 - w11;
  w22 = 1.0 - w21;
  w32 = 1.0 - w31;
  w42 = 1.0 - w41; 

  /* ----- interpolate variable ----- */
  for (iV=iV0; iV<nV; iV++) {
    V = w41 * (w31*( w21*( w11*cT->data[cti4(i,  j,  k,  l,  iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j,  k,  l,  iV,cT->n)] )
                    +w22*( w11*cT->data[cti4(i,  j+1,k,  l,  iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j+1,k,  l,  iV,cT->n)] ) )
              +w32*( w21*( w11*cT->data[cti4(i,  j,  k+1,l,  iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j,  k+1,l,  iV,cT->n)] )
                    +w22*( w11*cT->data[cti4(i,  j+1,k+1,l,  iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j+1,k+1,l,  iV,cT->n)] ) ))
       + w42* (w31*( w21*( w11*cT->data[cti4(i,  j,  k,  l+1,iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j,  k,  l+1,iV,cT->n)] )
                    +w22*( w11*cT->data[cti4(i,  j+1,k,  l+1,iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j+1,k,  l+1,iV,cT->n)] ) )
              +w32*( w21*( w11*cT->data[cti4(i,  j,  k+1,l+1,iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j,  k+1,l+1,iV,cT->n)] )
                    +w22*( w11*cT->data[cti4(i,  j+1,k+1,l+1,iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j+1,k+1,l+1,iV,cT->n)] ) ));
    if (!all) return V;
    vars[iV-offset] = V;
  } /* end loop over variables */

  return V;
}
/* -------------------------------------------------------------------------- */

double chemTableInterp4DIndex(int *index, double *x, double *vars, chemTab *cT)
{
  int    i,j,k,l,ig,jg,kg,lg,iV,iV0,nV,all,offset = 0;
  double V;
  double w11,w12,w21,w22,w31,w32,w41,w42,w,wtest;
    all = 0;
    iV0 = *index;    
    nV  = iV0+1;
  
  /* ----- look up co-ordinate location indices ----- */
  i = gridLocate(cT->n[0], cT->x[0], x[0]); 
  j = gridLocate(cT->n[1], cT->x[1], x[1]); 
  k = gridLocate(cT->n[2], cT->x[2], x[2]); 
  l = gridLocate(cT->n[3], cT->x[3], x[3]);

  /* ----- compute the weights ----- */
  w11 = (cT->x[0][i+1] - x[0])/(cT->x[0][i+1] - cT->x[0][i]);
  w21 = (cT->x[1][j+1] - x[1])/(cT->x[1][j+1] - cT->x[1][j]);
  w31 = (cT->x[2][k+1] - x[2])/(cT->x[2][k+1] - cT->x[2][k]);
  w41 = (cT->x[3][l+1] - x[3])/(cT->x[3][l+1] - cT->x[3][l]);
  w11 = fmin(fmax(w11, 0.0), 1.0);
  w21 = fmin(fmax(w21, 0.0), 1.0);
  w31 = fmin(fmax(w31, 0.0), 1.0);
  w41 = fmin(fmax(w41, 0.0), 1.0);
  w12 = 1.0 - w11;
  w22 = 1.0 - w21;
  w32 = 1.0 - w31;
  w42 = 1.0 - w41;
  
  /* ----- interpolate variable ----- */
  for (iV=iV0; iV<nV; iV++) {
    V = w41 * (w31*( w21*( w11*cT->data[cti4(i,  j,  k,  l,  iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j,  k,  l,  iV,cT->n)] )
                    +w22*( w11*cT->data[cti4(i,  j+1,k,  l,  iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j+1,k,  l,  iV,cT->n)] ) )
              +w32*( w21*( w11*cT->data[cti4(i,  j,  k+1,l,  iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j,  k+1,l,  iV,cT->n)] )
                    +w22*( w11*cT->data[cti4(i,  j+1,k+1,l,  iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j+1,k+1,l,  iV,cT->n)] ) ))
       + w42* (w31*( w21*( w11*cT->data[cti4(i,  j,  k,  l+1,iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j,  k,  l+1,iV,cT->n)] )
                    +w22*( w11*cT->data[cti4(i,  j+1,k,  l+1,iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j+1,k,  l+1,iV,cT->n)] ) )
              +w32*( w21*( w11*cT->data[cti4(i,  j,  k+1,l+1,iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j,  k+1,l+1,iV,cT->n)] )
                    +w22*( w11*cT->data[cti4(i,  j+1,k+1,l+1,iV,cT->n)]
                          +w12*cT->data[cti4(i+1,j+1,k+1,l+1,iV,cT->n)] ) ));
    if (!all) return V;
    vars[iV-offset] = V;
  } /* end loop over variables */

  return V;
}
/* -------------------------------------------------------------------------- */

double chemTableInterp3DIndex(int *index, double *x, double *vars, chemTab *cT)
{
  int    i,j,k,ig,jg,kg,lg,iV,iV0,nV,all,offset = 0;
  double V;
  double w11,w12,w21,w22,w31,w32,w,wtest;
    all = 0;
    iV0 = *index;    
    nV  = iV0+1;
  
  /* ----- look up co-ordinate location indices ----- */
  i = gridLocate(cT->n[0], cT->x[0], x[0]); 
  j = gridLocate(cT->n[1], cT->x[1], x[1]); 
  k = gridLocate(cT->n[2], cT->x[2], x[2]); 

  /* ----- compute the weights ----- */
  w11 = (cT->x[0][i+1] - x[0])/(cT->x[0][i+1] - cT->x[0][i]);
  w21 = (cT->x[1][j+1] - x[1])/(cT->x[1][j+1] - cT->x[1][j]);
  w31 = (cT->x[2][k+1] - x[2])/(cT->x[2][k+1] - cT->x[2][k]);
  w11 = fmin(fmax(w11, 0.0), 1.0);
  w21 = fmin(fmax(w21, 0.0), 1.0);
  w31 = fmin(fmax(w31, 0.0), 1.0);
  w12 = 1.0 - w11;
  w22 = 1.0 - w21;
  w32 = 1.0 - w31;
  
  /* ----- interpolate variable ----- */
  for (iV=iV0; iV<nV; iV++) {
    V =  w31*( w21*( w11*cT->data[cti3(i,  j,  k,   iV,cT->n)]
                    +w12*cT->data[cti3(i+1,j,  k,   iV,cT->n)] )
              +w22*( w11*cT->data[cti3(i,  j+1,k,   iV,cT->n)]
                    +w12*cT->data[cti3(i+1,j+1,k,   iV,cT->n)] ) )
        +w32*( w21*( w11*cT->data[cti3(i,  j,  k+1, iV,cT->n)]
                    +w12*cT->data[cti3(i+1,j,  k+1, iV,cT->n)] )
              +w22*( w11*cT->data[cti3(i,  j+1,k+1, iV,cT->n)]
                    +w12*cT->data[cti3(i+1,j+1,k+1, iV,cT->n)] ) );
    if (!all) return V;
    vars[iV-offset] = V;
  } /* end loop over variables */

  return V;
}
/* ------------------------------------------------------------------------- */

int  getChemTabNdim(void *tab)
{
  chemTab *cT = (chemTab*) tab;
  return cT->nDim;
}
/* ------------------------------------------------------------------------- */

int  getChemTabNvar(void *tab)
{
  chemTab *cT = (chemTab*) tab;
  return cT->nVar;
}
/* ------------------------------------------------------------------------- */

int  getChemTabNspecies(void *tab)
{
  chemTab *cT = (chemTab*) tab;
  return cT->nSpecies;
}
/* ------------------------------------------------------------------------- */

double getChemTabPref(void *tab)
{
  chemTab *cT = (chemTab*) tab;
  return cT->Pref;
}
/* ------------------------------------------------------------------------- */

double getChemTabTref(void *tab)
{
  chemTab *cT = (chemTab*) tab;
  return cT->Tref;
}
/* ------------------------------------------------------------------------- */

void getChemTabVarNames(char **names, void *tab)
{
  chemTab *cT = (chemTab*) tab;
  int i;
  for (i=0; i<cT->nVar; i++) strcpy(names[i], cT->names[i]);
}
/* ------------------------------------------------------------------------- */

void getChemTabSpecies(char **species, void *tab)
{
  chemTab *cT = (chemTab*) tab;
  int i;
  for (i=0; i<cT->nSpecies; i++) strcpy(species[i], cT->species[i]);
}
/* ------------------------------------------------------------------------- */

void printChemTableInfo(void *tab)
{
  int i;
  /* dereference struct */
  chemTab *cT = (chemTab*) tab;

  for (i=0; i<31; i++) fprintf(stderr, "-");
  fprintf(stderr, " Chemtable Summary ");
  for (i=0; i<30; i++) fprintf(stderr, "-");
  fprintf(stderr, "\n");

  fprintf(stderr, "  Model:        %s\n", cT->model);
  fprintf(stderr, "  Format:       %s\n", cT->format);
  fprintf(stderr, "  Dimensions:   %i\n", cT->nDim);
  fprintf(stderr, "  Co-ordinates:   \n");
  fprintf(stderr, "%7s%10s | %6s | %10s | %10s\n","","name","size","min","max");
  for (i=0; i<cT->nDim; i++)
    fprintf(stderr, "%4i | %10s | %6i | %10.2f | %10.2f\n", i, "", cT->n[i],
            cT->xmin[i], cT->xmax[i]);
  fprintf(stderr, "\n");
  fprintf(stderr, "  Variables (%i) min/max:\n", cT->nVar);
  for (i=0; i<cT->nVar; i++)
    fprintf(stderr, "%3s%-25s %15.6e %15.6e\n", "",cT->names[i],cT->datMin[i],
            cT->datMax[i]);
  fprintf(stderr, "\n  Species (%i):\n", cT->nSpecies);
  for (i=0; i<cT->nSpecies; i++) {
    fprintf(stderr, "%3s%-24s ", "", cT->species[i]);
    if ((i+1)%3 == 0)
      fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  for (i=0; i<80; i++) fprintf(stderr, "-");
  fprintf(stderr, "\n");
}
/* ------------------------------------------------------------------------- */

int indxC(int nx, int *i, int *n)
{
  int k, loc = i[0]*n[1];
  for (k=1; k<nx-2; k++) {
    loc += i[k];
    loc *= n[k+1];
  }
  loc += i[nx-1];
  return loc;
}
/* ------------------------------------------------------------------------- */

int indxF(int nx, int *i, int *n)
{
  int k, loc = i[nx-1]*n[nx-2];
  for (k=nx-2; k>=1; k--) {
    loc += i[k];
    loc *= n[k-1];
  }
  loc += i[0];
  return loc;
}
/* ------------------------------------------------------------------------- */

int cti2(int i, int j, int k, int *n)
{
  return i + n[0]*(j + n[1]*k);
}  
/* ------------------------------------------------------------------------- */

int cti3(int i, int j, int k, int p, int *n)
{
  return i + n[0]*(j + n[1]*(k + n[2]*p));
}  
/* ------------------------------------------------------------------------- */

int cti4(int i, int j, int k,int l, int p, int *n)
{
  return i + n[0]*(j + n[1]*(k + n[2]*(l + n[3]*p)));
}  
/* ------------------------------------------------------------------------- */

int cTab_indx(int i, int j, int k, int p, int q, int *n)
{
  return i + n[0]*(j + n[1]*(k + n[2]*(p + n[3]*q)));
}  
/* ------------------------------------------------------------------------- */

int cTab_indxf(int i, int j, int k, int p, int q, int *n)
{
  return i-1 + n[0]*(j-1 + n[1]*(k-1 + n[2]*(p-1 + n[3]*(q-1))));
}
/* ------------------------------------------------------------------------- */
