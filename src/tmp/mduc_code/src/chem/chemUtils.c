#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "chemUtils.h"
#include "chemMech.h"
#include "parser.h"
#include "utils.h"
#include "display.h"

rxnStruct * paramGetReaction(char *name, parseStruct *p)
{
  int i,j,k,iP,iSp,nSp = 0;
  int    rSide = 1;
  double lastNu = -1.0, *nu;
  char tag[64], nuC[16], buffer[256], *pS, splits[16],**raw,**sp;

  char alphabet[52] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

  /* allocate the memory to return species and coefficients */
  rxnStruct *rxn = (rxnStruct*) calloc(1, sizeof(rxnStruct));
  setRxnDefault(rxn);
  
  /* find the reaction in the input */
  iP  = paramGetIndex(name,p);
  if (iP < 0) return NULL;
  raw = p->vals[iP];

  /* loop over the inputs */
  rxn->nS  = rxn->nR = rxn->nP = 0;
  rxn->spR = (char**)  calloc(1, sizeof(char*));
  rxn->spP = (char**)  calloc(1, sizeof(char*));
  rxn->nuR = (double*) calloc(1, sizeof(double));
  rxn->nuP = (double*) calloc(1, sizeof(double));
  for (j=0; j<p->n[iP]; j++) {
    strcpy(tag, raw[j]);

    /* check whether reaction is reversible */
    /* NOTE: this requires a space, i.e. R <=> P */
    if (!strcmp(tag, "<")) {
      rxn->reverse = 1;
      rxn->dir = 'f';
      continue;
    }
    if (!strcmp(tag, ">")) {
      rSide = 0;
      continue;
    }

    /* get the coefficient for the species */
    i = 0;
    while (isdigit(tag[i]) || tag[i] == '.' && i < strlen(tag)) {
      i++;
    };
    
    if (i == strlen(tag)) { /* this was only a digit */
      lastNu = atof(tag);
      continue;
    }

    /* count the species and copy necessary values */
    if (rSide) {
      rxn->nR++;
      nSp = rxn->nR;
      iSp = nSp - 1;
      sp = rxn->spR;
      nu = rxn->nuR;
    }
    else {
      rxn->nP++;
      nSp = rxn->nP;
      iSp = nSp - 1;
      sp = rxn->spP;
      nu = rxn->nuP;
    }
        
    sp      = (char**)  realloc(sp, nSp*sizeof(char*));
    sp[iSp] = (char*)    calloc(32, sizeof(char));
    nu      = (double*) realloc(nu, nSp*sizeof(double));

    if (lastNu > 0.0) {
      strcpy(sp[iSp], tag);
      nu[iSp] = lastNu;
      lastNu = -1.0;
    }
    else {
      /* NOTE: i is set above and is first nondigit character */
      if (i == 0 || (tag[i] == '-')) { /* no leading digit found */
        if      (tag[0] == '1') tag[0] = 'P';
        else if (tag[0] == '2') tag[0] = 'S';
        else if (tag[0] == '3') tag[0] = 'T';
        for (k=0; k<strlen(tag); k++) tag[k] = toupper(tag[k]);
        strcpy(sp[iSp], tag);
        nu[iSp] = 1;
      }
      else { /* leading digit found */
        strcpy(sp[iSp], &tag[i]);
        strncpy(nuC, tag, i);
        nuC[i] = '\0';
        nu[iSp] = atof(nuC);
      }
    }
  }

  rxn->nR = removeDuplicateSpecies(rxn->nR,&rxn->nuR,&rxn->spR);
  rxn->nP = removeDuplicateSpecies(rxn->nP,&rxn->nuP,&rxn->spP);

  /* process the overall species list */
  rxn->nS = rxn->nR + rxn->nP;
  rxn->sp  = (char**)  calloc(rxn->nS, sizeof(char*));
  rxn->nu  = (double*) calloc(rxn->nS, sizeof(double));
  for (i=0; i<rxn->nR; i++) {
    rxn->sp[i] = (char*) calloc(32, sizeof(char));
    strcpy(rxn->sp[i],rxn->spR[i]);
    rxn->nu[i] = -rxn->nuR[i];
  }
  for (i=0,j=rxn->nR; i<rxn->nP; i++,j++) {
    rxn->sp[j] = (char*) calloc(32, sizeof(char));
    strcpy(rxn->sp[j],rxn->spP[i]);
    rxn->nu[j] = rxn->nuP[i];
  }
  rxn->nS = removeDuplicateSpecies(rxn->nS,&rxn->nu,&rxn->sp);

  return rxn;
}
/* -------------------------------------------------------------------------- */

void initChemFuel(chemStruct *chem)
{
  int i;
  int convertYtoX = -1;
  double sum;

  chem->sFuel = (int*) calloc(chem->nFuel, sizeof(int));
  /* check whether mole or mass fractions specified */
  if (chem->fuelY == NULL) {
    chem->fuelY = (double*) calloc(chem->nFuel, sizeof(double));
    convertYtoX += 1;
    /* check if sum of mole fractions are equal to 1 */
    sum 	= 0.0;
    for (i=0; i<chem->nFuel; i++)
	sum += chem->fuelX[i];
    if (abs(sum-1)>1e-3)
	fprintf(stderr, "Warning: sum of fuel mole fractions: %6.5e \n", sum);
  }
  if (chem->fuelX == NULL) {
    chem->fuelX = (double*) calloc(chem->nFuel, sizeof(double));
    convertYtoX += 2;
    /* check if sum of mass fractions are equal to 1 */
    sum 	= 0.0;
    for (i=0; i<chem->nFuel; i++)
	sum += chem->fuelY[i];
    if (abs(sum-1)>1e-3)
	fprintf(stderr, "Warning: sum of fuel mass fractions: %6.5e \n", sum);
	
  }
  
  /* get the species index from the mechanism */
  for (i=0; i<chem->nFuel; i++) {
    chem->sFuel[i] = getSpeciesIndex(chem->nSpecies, chem->species,
                                     chem->fuel[i]);
  }

  if (convertYtoX > 1) {
    fprintf(stderr, "Warning: fuel composition not defined\n");
    return;
  }

  /* convert mass/mole fractions */
  double *W = calloc(chem->nFuel, sizeof(double));
  for (i=0; i<chem->nFuel; i++)
    W[i] = chem->W[chem->sFuel[i]];

  if (convertYtoX) 
    convertMassToMoleFrac(chem->nFuel, W, chem->fuelY, chem->fuelX);
  else
    convertMoleToMassFrac(chem->nFuel, W, chem->fuelX, chem->fuelY);

  chem->fuelMw = 0.0;
  for (i=0; i<chem->nFuel; i++)
    chem->fuelMw += chem->fuelY[i]*W[i];

  if (MDUC_PRINT_PROC) {
    fprintf(stdout,"--------- Stoichiometric properties ---------\n");
    if (!chem->compiled) {
      /* get the stoichiometric properties */
      chem->stoichFA = calcStoichMixRatio(chem);
      fprintf(stdout,"\n  stoich. F/A: %8.3f\n", chem->stoichFA);
      fprintf(stdout,"  stoich. A/F: %8.3f\n", 1.0/chem->stoichFA);
      fprintf(stdout,"  H/C ratio:   %8.3f\n", chem->fuelH/chem->fuelC);
      fprintf(stdout,"  LHV:         %8.3f [MJ/kg]\n\n", chem->fuelLHV*1e-6);
    }

    fprintf(stdout,"---------------------------------------------\n");
    fprintf(stdout," %-20s |    Y     |    X     |\n","species");
    fprintf(stdout,"---------------------------------------------\n");
    /* fprintf(stderr, "# species mass-fraction mole fraction\n"); */
    for (i=0; i<chem->nFuel; i++) {
      /* if (chem->fuelY[i] > 1.0e-6) */
      fprintf(stdout," %-20s | %8.6f | %8.6f |\n",chem->fuel[i],chem->fuelY[i],
              chem->fuelX[i]);
    }
    fprintf(stdout,"---------------------------------------------\n");

  }

  free(W);
}
/* -------------------------------------------------------------------------- */

double calcStoichMixRatio(chemStruct *chem)
{
  int i, j;
  double N_O = 3.77418; /* ratio of (air-O2)/O2  */
  double aC, aH, aO, nTot;
  double m, mFuel, mAir, mO2;
  molStruct *sp;

  /* zero out */
  for (i=0; i<chem->nSpecies; i++) chem->stoichX[i] = 0.0;

  aC = aH = aO = 0;
  mFuel = 0.0;
  for (i=0; i<chem->nFuel; i++) {
    sp = &chem->specs[chem->sFuel[i]];
    for (j=0; j<sp->nAtoms; j++) {
      if (!strcmp(sp->atom[j]->symb, "C"))
        aC += chem->fuelX[i]*sp->nA[j];
      else if (!strcmp(sp->atom[j]->symb, "H"))
        aH += chem->fuelX[i]*sp->nA[j];
      else if (!strcmp(sp->atom[j]->symb, "O"))
        aO += chem->fuelX[i]*sp->nA[j];
      else
        fprintf(stderr, "Stoichiometry warning: fuel atom %s not recognized\n",
                sp->atom[j]->symb);
    } /* end loop over species atoms */

    chem->stoichX[chem->sFuel[i]] = chem->fuelX[i];
    mFuel += chem->fuelX[i]*chem->W[chem->sFuel[i]];
  } /* end loop over fuel species */
  
  /* store the number of hydrocarbon elements in each fuel */
  chem->fuelC = aC;
  chem->fuelH = aH;
  chem->fuelO = aO;

  /* calculate stoichiometric moles */
  chem->stoichX[chem->sCO2] = aC;
  chem->stoichX[chem->sH2O] = 0.5*aH;
  chem->stoichX[chem->sO2]  = 0.5*(2*aC + aH - aO);
  chem->stoichX[chem->sN2]  = chem->stoichX[chem->sO2]*N_O;

  mO2   = chem->stoichX[chem->sO2]*chem->W[chem->sO2];
  mAir  = mO2 + chem->stoichX[chem->sN2]*chem->W[chem->sN2];

  if (chem->fuelLHV < -1e20) { /* allow forcing from input file */
    /* calculate the lower heating value of the fuel in J/(kg fuel) */
    chemComputeThermoData(chem->h, chem->cp, 298.15, chem);
    /* products - reactants */
    chem->fuelLHV =
      chem->stoichX[chem->sCO2]*chem->h[chem->sCO2]*chem->W[chem->sCO2]
      + chem->stoichX[chem->sH2O]*chem->h[chem->sH2O]*chem->W[chem->sH2O]
      - chem->stoichX[chem->sO2]*chem->h[chem->sO2]*chem->W[chem->sO2]
      - chem->stoichX[chem->sN2]*chem->h[chem->sN2]*chem->W[chem->sN2];
    chem->fuelLHV *= 1e-3; /* J */
    for (i=0; i<chem->nFuel; i++) {
      chem->fuelLHV -=
        chem->fuelX[i]*chem->h[chem->sFuel[i]]*chem->W[chem->sFuel[i]]*1e-3;
    }
    chem->fuelLHV *= -1.0;  /* convention */
    chem->fuelLHV /= mFuel; /* J/g */
    chem->fuelLHV *= 1000;  /* J/kg */
  }

  return mFuel/mAir;
}
/* -------------------------------------------------------------------------- */

int getSpeciesIndices(int nCheck, char **spCheck, int *iSp, int nM, char **spM)
{
  if (nCheck < 1) return 0;

  int i, allFound=0;
  for (i=0; i<nCheck; i++) {
    iSp[i] = getSpeciesIndex(nM, spM, spCheck[i]);
    if (iSp >= 0)
      allFound++;
  }
  return allFound - nCheck;
}
/* -------------------------------------------------------------------------- */

void copyRxn(rxnStruct *rc, rxnStruct *r)
{
  int i;

  /* clean up memory first and reallocate in this function */
  destroyRxn(rc);

  /* simple info */
  rc->no      = r->no;
  rc->dir     = r->dir;
  rc->reverse = r->reverse;
  rc->orderF[0] = r->orderF[0];
  rc->orderF[1] = r->orderF[1];
  rc->orderR[0] = r->orderR[0];
  rc->orderR[1] = r->orderR[1];

  /* copy and allocate memory */
  rc->nS = r->nS;
  rc->nR = r->nR;
  rc->nP = r->nP;
  rc->nA = r->nA;
  rc->iS  = (int*)    calloc(rc->nS, sizeof(int));
  rc->iR  = (int*)    calloc(rc->nR, sizeof(int));
  rc->iP  = (int*)    calloc(rc->nP, sizeof(int));
  rc->nu  = (double*) calloc(rc->nS, sizeof(double));
  rc->nuR = (double*) calloc(rc->nR, sizeof(double));
  rc->nuP = (double*) calloc(rc->nP, sizeof(double));
  rc->sp  = (char**)  calloc(rc->nS, sizeof(char*));
  rc->spR = (char**)  calloc(rc->nR, sizeof(char*));
  rc->spP = (char**)  calloc(rc->nP, sizeof(char*));
  for (i=0; i<rc->nS; i++) {
    rc->sp[i] = (char*) calloc(32, sizeof(char));
    rc->nu[i] = r->nu[i];
    rc->iS[i] = r->iS[i];
    strcpy(rc->sp[i], r->sp[i]);
  }
  for (i=0; i<rc->nR; i++) {
    rc->spR[i] = (char*) calloc(32, sizeof(char));
    rc->nuR[i] = r->nuR[i];
    rc->iR[i]  = r->iR[i];
    strcpy(rc->spR[i], r->spR[i]);
  }
  for (i=0; i<rc->nP; i++) {
    rc->spP[i] = (char*) calloc(32, sizeof(char));
    rc->nuP[i] = r->nuP[i];
    rc->iP[i]  = r->iP[i];
    strcpy(rc->spP[i], r->spP[i]);
  }

/*   rc->atoms = (char**) calloc(rc->nA, sizeof(char*)); */
/*   for (i=0; i<rc->nA; i++) { */
/*     rc->atoms[i] = (char*) calloc(3, sizeof(char)); */
/*     strcpy(rc->atoms[i], r->atoms[i]); */
/*   } */

  rc->A      = r->A;
  rc->b      = r->b;
  rc->E      = r->E;
  rc->Ar     = r->Ar;
  rc->br     = r->br;
  rc->Er     = r->Er;
  rc->A_low  = r->A_low;
  rc->b_low  = r->b_low;
  rc->E_low  = r->E_low;
  rc->Ar_low = r->Ar_low;
  rc->br_low = r->br_low;
  rc->Er_low = r->Er_low;

  rc->thirdBody = r->thirdBody;
  rc->pressure  = r->pressure;
  rc->lindemann = r->lindemann;
  rc->troe      = r->troe;
  rc->sri       = r->sri;
  for (i=0; i<3; i++) rc->troeA[i] = r->troeA[i];
  for (i=0; i<3; i++) rc->troeT[i] = r->troeT[i];
  for (i=0; i<5; i++) rc->sriA[i]  = r->sriA[i];

  rc->nEnhance = r->nEnhance;
  rc->spEnh    = (char**)  calloc(rc->nEnhance, sizeof(char*));
  rc->factEnh  = (double*) calloc(rc->nEnhance, sizeof(double));
  for (i=0; i<rc->nEnhance; i++) {
    rc->spEnh[i] = (char*) calloc(32, sizeof(char));
    strcpy(rc->spEnh[i], r->spEnh[i]);
    rc->factEnh[i] = r->factEnh[i];
  }
  rc->nSpecies = r->nSpecies;
  if (rc->pressure > 0) {
    rc->concF = (double*) calloc(rc->nSpecies, sizeof(double));
    for (i=0; i<rc->nSpecies; i++)
      rc->concF[i] = r->concF[i];
  }
}
/* -------------------------------------------------------------------------- */

void destroySpecies(molStruct *sp)
{
  if (sp->nA)    free(sp->nA);
  if (sp->atom)  free(sp->atom);
  if (sp->iRxn)  free(sp->iRxn);
  if (sp->nuRxn) free(sp->nuRxn);
}
/* -------------------------------------------------------------------------- */

void setRxnDefault(rxnStruct *r)
{
  r->thirdBody = 0;
  r->pressure = r->lindemann = r->troe = r->sri = 0;
  r->b = r->b_low = r->br = r->br_low = 0.0;
  r->reverse  = 0;
  r->nEnhance = 0;
  r->nSpecies = 0;
  r->nS  = r->nR  = r->nP  = 0;
  r->iS  = r->iR  = r->iP  = NULL;
  r->nu  = r->nuR = r->nuP = r->factEnh = r->concF = NULL;
  r->sp  = r->spR = r->spP = r->atoms   = r->spEnh = NULL;
}
/* -------------------------------------------------------------------------- */

void destroyRxn(rxnStruct *r)
{
  int i;
  if (r->nu)  free(r->nu);
  if (r->nuR) free(r->nuR);
  if (r->nuP) free(r->nuP);
  if (r->iS)  free(r->iS);
  if (r->iR)  free(r->iR);
  if (r->iP)  free(r->iP);
  if (r->sp) {
    for (i=0; i<r->nS; i++) if (r->sp[i]) free(r->sp[i]);
    free(r->sp);
  }
  if (r->spR) {
    for (i=0; i<r->nR; i++) if (r->spR[i]) free(r->spR[i]);
    free(r->spR);
  }
  if (r->spP) {
    for (i=0; i<r->nP; i++) if (r->spP[i]) free(r->spP[i]);
    free(r->spP);
  }
  if (r->atoms) {
    for (i=0; i<r->nA; i++) free(r->atoms[i]);
    free(r->atoms);
  }

  if (r->spEnh) {
    for (i=0; i<r->nEnhance; i++) if (r->spEnh[i]) free(r->spEnh[i]);
    free(r->spEnh);
  }
  if (r->factEnh) free(r->factEnh);
  if (r->concF)   free(r->concF);
}
/* -------------------------------------------------------------------------- */

void destroyChem(chemStruct *chem)
{
  int i;
  for (i=0; i<chem->nSpecies; ++i) if (chem->species[i]) free(chem->species[i]);
  if (chem->species) free(chem->species);

  if (chem->pt) free(chem->pt);
  
  if (chem->W)    free(chem->W);
  if (chem->Wi)    free(chem->Wi);
  if (chem->cp)   free(chem->cp);
  if (chem->h)    free(chem->h);
  if (chem->cdot) free(chem->cdot);
  if (chem->c)    free(chem->c);
  if (chem->K)    free(chem->K);
  if (chem->w)    free(chem->w);
  if (chem->w_tmp) free(chem->w_tmp);
  if (chem->w_pre) free(chem->w_pre);
  if (chem->w_int) free(chem->w_int);
  if (chem->w_int_pre) free(chem->w_int_pre);
  if (chem->M)    free(chem->M);
  if (chem->kF)   free(chem->kF);
  if (chem->kR)   free(chem->kR);
  if (chem->wF)   free(chem->wF);
  if (chem->wR)   free(chem->wR);

  if (chem->iRad) free(chem->iRad);
  if (chem->pRad) free(chem->pRad);
  if (chem->alpha) free(chem->alpha);

  if (chem->transport) free(chem->transport);
  if (chem->eps_k)     free(chem->eps_k);
  if (chem->sigma)     free(chem->sigma);
  if (chem->coeffMu)   free(chem->coeffMu);
  if (chem->coeffD)    free(chem->coeffD);
   
  if (chem->readTrans) {
     for (i=0; i<chem->nSpecies; ++i) {
        if (chem->phiMixMu_a[i])   free(chem->phiMixMu_a[i]);
        if (chem->WjOverWi[i])        free(chem->WjOverWi[i]);  
        if (chem->Wij_D[i])           free(chem->Wij_D[i]);     
  }
  if (chem->phiMixMu_a)      free(chem->phiMixMu_a);
  if (chem->WjOverWi)        free(chem->WjOverWi);  
  if (chem->Wij_D)           free(chem->Wij_D);     
  if (chem->mu05)	     free(chem->mu05);	     
  if (chem->mu05i)           free(chem->mu05i);     
  }
//  for (i=0; i<chem->nSootSp;++i)
//      if (chem->sootRxnSpec[i]) free(chem->sootRxnSpec[i]);
//     
//  if (chem->sootRxnSpec)     free(chem->sootRxnSpec);       
  if (chem->sootRxnConc)     free(chem->sootRxnConc);        
  if (chem->sootRxnRates)    free(chem->sootRxnRates);       

  for (i=0; i<chem->nSootSrc;++i) {
     if (chem->sootSrcNames[i])    free(chem->sootSrcNames[i]);
//     if (chem->sootSrcSelect[i])   free(chem->sootSrcSelect[i]);
  }
  if (chem->sootSrcNames)    free(chem->sootSrcNames);
//  if (chem->sootSrcSelect)   free(chem->sootSrcSelect);     
  if (chem->sootSrcOn)       free(chem->sootSrcOn);          

//  for (i=0; i<chem->nSootMom;++i) {
//     if (chem->sootMomNames[i])    free(chem->sootMomNames[i]);
//     if (chem->sootMomAbbrv[i])    free(chem->sootMomAbbrv[i]);
//     if (chem->sootMomSelect[i])   free(chem->sootMomSelect[i]);
//  }
//  if (chem->sootMomNames)    free(chem->sootMomNames);
//  if (chem->sootMomAbbrv)    free(chem->sootMomAbbrv);
//  if (chem->sootMomSelect)   free(chem->sootMomSelect);
  if (chem->sootMomOn)       free(chem->sootMomOn);          

  for (i=0; i<chem->nSootPAH;++i) if (chem->sootPAHSpec[i]) free(chem->sootPAHSpec[i]);
  if (chem->sootPAHSpec)     free(chem->sootPAHSpec);       

  for (i=0; i<chem->nFuel; ++i) if (chem->fuel[i]) free(chem->fuel[i]);
  if (chem->fuel)    free(chem->fuel);
  if (chem->sFuel)   free(chem->sFuel);
  if (chem->fuelY)   free(chem->fuelY);
  if (chem->fuelX)   free(chem->fuelX);
  if (chem->stoichX) free(chem->stoichX);

  if (chem->specs) {
    for (i=0; i<chem->nSpecies; ++i) destroySpecies(&chem->specs[i]);
    free(chem->specs);
  }
  if (chem->rxns) {
    for (i=0; i<chem->nRxn; ++i) destroyRxn(&chem->rxns[i]);
    free(chem->rxns);
  }
}
/* -------------------------------------------------------------------------- */
