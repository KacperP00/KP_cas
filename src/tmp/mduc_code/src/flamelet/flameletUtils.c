#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "flameletUtils.h"
#include "chemUtils.h"
#include "parser.h"
#include "display.h"
#include "utils.h"

flameStruct *paramGetFlamelet(char *name, driveStruct *driver, chemStruct *chem,
                              flowStruct *flow, parseStruct *p)
{
  int i,j,k,nP,iP,iB,iSp,nC;
  char **raw,tag[64],buffer[256],*pS;

  /* error checking */
  if (driver == NULL)
    mducErr("driver must be initialized before flamelet","flameletInit");
  if (chem == NULL)
    mducErr("chemistry must be initialized before flamelet","flameletInit");
  if (flow == NULL)
    mducErr("flow must be initialized before flamelet","flameletInit");
  
  /* find the reaction in the input */
  iP  = paramGetIndex("FLAMELET",p);
  if (iP < 0) {
    mducErr("could not find flamelet in input file","flameletInit");
    return NULL;
  }
  nP  = p->n[iP];
  raw = p->vals[iP];

  /* allocate the memory to return species and coefficients */
  flameStruct *flame = (flameStruct*) calloc(1, sizeof(flameStruct));

  /* set global defaults */
  flame->RGAS    = 8314.34; /* universal gas constant J/kg-K */
  flame->intSrcQ = 0.0;

  /* initialize arrays to NULL */
  flame->x1 = flame->x2 = flame->dx1 = flame->dx2 = NULL;
  flame->v1 = flame->v2 = flame->chiA = flame->chiB = flame->chiC = NULL;
  flame->Y = flame->Yp = flame->Yp_m1 = flame->Y_m1 = flame->Y_m2 = NULL;
  flame->Y_pre = NULL;
  flame->convCoeff = flame->chemSrcY = flame->chemSrcY_tmp = NULL;
  flame->rho = flame->cpmix = flame->sumCp = flame->sumCp2 = flame->T = NULL;
  flame->Mw = flame->R = flame->mu = flame->D = flame->gamma = NULL;
  flame->srcH = flame->srcP = flame->srcQ = flame->HRrate = flame->H = NULL;
  flame->chem_src = flame->diff_src = flame->conv_src = flame->visc_src = NULL;
  flame->chem_temp = flame->diff_temp = flame->conv_temp = flame->gPDF = NULL;
  flame->dY.x = flame->dY.y = flame->dY.xx = flame->dY.yy = flame->dY.xy = NULL;
  flame->dCp.x  = flame->dCp.y   = NULL;
  flame->dCp.xx = flame->dCp.yy  = flame->dCp.xy = NULL;
  flame->grad.x = flame->grad.xx = flame->grad.y = flame->grad.yy = NULL;
  flame->sootMom = flame->sootMomSrc = flame->sootMomSrcPart = NULL;

  /* get the basic numbers */
  flame->nDim = sparamGetInt(nP,"dimensions",1,raw);
  flame->n1   = sparamGetInt(nP,"gridPoints1",101,raw);
  flame->n2   = sparamGetInt(nP,"gridPoints2",1,raw);
  flame->nPiv1 = sparamGetRealArray(nP,"gridPivot1",&flame->gPiv1,raw);
  nC           = sparamGetRealArray(nP,"gridPercent1",&flame->gPc1,raw);
  if (flame->nPiv1 != nC)
    fprintf(stderr, "Warning: pivots/percentages given do not match (%i/%i)\n",
            flame->nPiv1,nC);
  flame->nPiv2 = sparamGetRealArray(nP,"gridPivot2",&flame->gPiv2,raw);
  nC           = sparamGetRealArray(nP,"gridPercent2",&flame->gPc2,raw);
  if (flame->nPiv2 != nC)
    fprintf(stderr, "Warning: pivots/percentages given do not match (%i/%i)\n",
            flame->nPiv2,nC);

  /* read in any global quantities (e.g. if 0D) */
  flame->phi   = sparamGetReal(nP,"phi",1.0,raw);
  flame->Tmean = sparamGetReal(nP,"Temperature",1000,raw);
  
  /* get the boundary conditions */
  flame->stream[0] = paramGetBC("OxidizerSide",nP,raw);
  flame->stream[1] = paramGetBC("FuelSide1",nP,raw);
  flame->stream[2] = paramGetBC("FuelSide2",nP,raw);

  /* assign appropriately */
  /* TODO: enable this in input file? */
  flame->ox    = flame->stream[0];
  flame->fuel1 = flame->stream[1];
  flame->fuel2 = flame->stream[2];

  return flame;
}
/* -------------------------------------------------------------------------- */

bndStruct * paramGetBC(char *name,int n,char **inp)
{
  int iB,j,j0,jN,iSp;
  char **bcs;
  /* determine max value to read */
  iB = sparamGetIndex(n,name,inp);
  if (iB < 0) {
    fprintf(stderr, "Warning: can't find boundary condition %s\n",name);
    return NULL;
  }
  j0 = iB;
  while (strcmp(inp[j0],"{") && strcmp(inp[j0],"begin")) j0++;
  jN = j0+1;
  while (strcmp(inp[j0],"}") && strcmp(inp[jN],"end")) jN++;
  n = jN-j0;
  bcs = &inp[j0];

  /* allocate and set defaults */
  bndStruct *bnd = (bndStruct*) calloc(1, sizeof(bndStruct));
  strcpy(bnd->name, name);
  strcpy(bnd->sat,"none");
  bnd->type  = 'U'; 
  bnd->T     = -1;
  bnd->H     = -1e20;
  bnd->Hvap  = -1e20;
  bnd->cpLiq = -1;
  bnd->Tref  = -1;
  bnd->nSp   = 0;
  bnd->iSp   = (int*)    calloc(1, sizeof(int));
  bnd->sp    = (char**)  calloc(1, sizeof(char*));
  bnd->Y     = (double*) calloc(1, sizeof(double));
  bnd->X     = (double*) calloc(1, sizeof(double));

  /* process boundary condition */
  for (j=j0+1; j<jN-1; j++) {
    if (inp[j][0] == '[') /* unit */
      continue;
    if (!strncmp(inp[j], "Temperature",10)) {
      bnd->T = atof(inp[j+1]);
      bnd->type = 'T';
    }
    else if (!strcmp(inp[j],"Enthalpy")) {
      bnd->H = atof(inp[j+1]);
      bnd->type = 'H';
    }
    else if (!strcmp(inp[j],"Latent-heat")) {
      bnd->Hvap = atof(inp[j+1]);
      bnd->type = 'L';
    }
    else if (!strcmp(inp[j],"Liquid-cp"))
      bnd->cpLiq = atof(inp[j+1]);
    else if (!strcmp(inp[j],"Temp-ref"))
      bnd->Tref = atof(inp[j+1]);
    else if (!strcmp(inp[j],"sat-data")) {
      strcpy(bnd->sat,inp[j+1]);
      bnd->type = 'D';
    }
    else { /* species */
      iSp = bnd->nSp;
      bnd->nSp++;
      bnd->iSp     = (int*)   realloc(bnd->iSp, bnd->nSp*sizeof(int));
      bnd->sp      = (char**) realloc(bnd->sp,  bnd->nSp*sizeof(char*));
      bnd->sp[iSp] = (char*)  calloc(32, sizeof(char));
      bnd->Y = (double*) realloc(bnd->Y, bnd->nSp*sizeof(double));
      bnd->X = (double*) realloc(bnd->X, bnd->nSp*sizeof(double));
      if (!strncmp(inp[j],"X->",3)) {
        strcpy(bnd->sp[iSp], &inp[j][2]);
        bnd->X[iSp] = atof(inp[j+1]);
      }
      else if (!strncmp(inp[j],"Y->",3)) {
        strcpy(bnd->sp[iSp], &inp[j][2]);
        bnd->Y[iSp] = atof(inp[j+1]);
      }
      else {
        strcpy(bnd->sp[iSp], inp[j]);
        bnd->Y[iSp] = atof(inp[j+1]);
      }
    }

    j++;
  }
  return bnd;
}
/* -------------------------------------------------------------------------- */

void writeFlameletBC(FILE *pF,bndStruct *bnd)
{
  int i;
  if (bnd == NULL) return;
  fprintf(pF,"%s {\n",bnd->name);
  if (bnd->type == 'T')
    fprintf(pF,"%20s = %15.2f [K]\n","Temperature",bnd->T);
  else if (bnd->type == 'H')
    fprintf(pF,"%20s = %15.6e [J/kg]\n","Enthalpy",bnd->H);
  else if (bnd->type == 'L') {
    fprintf(pF,"%20s = %15.2f [K]\n","Temperature",bnd->T);
    fprintf(pF,"%20s = %15.6e [J/kg]\n","Latent-heat",bnd->Hvap);
    fprintf(pF,"%20s = %15.2f [J/kg-K]\n","Liquid-cp",bnd->cpLiq);
    fprintf(pF,"%20s = %15.2f [K]\n","Temp-ref",bnd->Tref);
  }
  else if (bnd->type == 'D')
    fprintf(pF,"%20s %s\n","sat-data",bnd->sat);

  for (i=0; i<bnd->nSp; i++)
    fprintf(pF, "%20s = %15.6f\n",bnd->sp[i],bnd->Y[i]);

  fprintf(pF, "}\n");
}
/* -------------------------------------------------------------------------- */
