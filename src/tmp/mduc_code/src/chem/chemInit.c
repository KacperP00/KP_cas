#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "chemInit.h"
#include "chemMech.h"
#include "transport.h"
#include "soot.h"
#include "chemUtils.h"
#include "flameletUtils.h"
#include "accessChem.h"
#include "utils.h"
#include "fileIO.h"
#include "chemKinIO.h"
#include "flameMasterIO.h"
#include "utils.h"
#include "fileIO.h"
#include "display.h"


int setChemDefaults(chemStruct *chem)
{
  strcpy(chem->transFile, "");
  chem->wiebeMech = 0;
  chem->nFuel = 0;
  chem->fuel  = NULL;
  chem->sFuel = NULL;
  chem->fuelY = chem->fuelX = NULL;
  chem->fuelC = chem->fuelH = chem->fuelO = -1;
  chem->fuelLHV = -1e30;
  
  chem->transport =  NULL;
  chem->atomGeom  = NULL;
  chem->eps_k = chem->sigma = NULL;
  chem->muDipole = chem->alphaPolar = chem->Zrot = NULL;
  chem->coeffMu = chem->coeffD = NULL;
  chem->phiMixMu_a = chem->WjOverWi = NULL;
  chem->mu05 = chem->mu05i = NULL;

  /* ----- set constants ----- */
  chem->logPoverR = log(101325/8314.14); /* P_atm/Rgas for chem equil. const. */

  /* ----- mechanism options ----- */
  strcpy(chem->reverseRxn, "fit");
  chem->fitReverse = 1;
  chem->revRxnDBG  = 0;
  strcpy(chem->unitsA,"cm-mol-s-K");
  strcpy(chem->unitsE,"CAL/MOL");

  /* ----- soot model defaults ----- */
  strcpy(chem->sootModel, "");
  chem->soot       = 0;
  chem->sootMode   = 0;
  chem->sootOrder  = 1;
  chem->nSootSrcOn = -1;
  chem->nSootMomOn = -1;
  chem->nSootPAH = -1;
  chem->nSootRxn = 0;
  chem->iSootSp = chem->sootSrcOn = chem->sootMomOn = NULL;
  chem->sootRxnRates = chem->sootRxnConc = NULL;
  chem->sootSrcNames = chem->sootMomNames = chem->sootMomAbbrv = NULL;
  chem->sootSrcSelect = chem->sootMomSelect = chem->sootRxnSpec = NULL;
  chem->sootPAHSpec = NULL;

  /* ----- total number if reactions ----- */
  chem->nRxnTotal = -1;

  return 1;
}
/* -------------------------------------------------------------------------- */

int chemInitPeriodicTable(chemStruct *chem)
{
  /* ----- initialize the atoms ----- */
  const double atom_M[120] = {
    1.00794,     4.002602,  6.941,      9.012182, 10.811,         12.0107,
    14.00674,    15.9994,   18.9984032, 20.1797,  22.989769,      24.305,
    26.9815386,  28.0855,   30.973762,  32.065,   35.453,         39.948,
    39.0983,     40.078,    44.955912,  47.867,   50.9415,        51.9961,
    54.938045,   55.845,    58.933195,  58.6934,  63.546,         65.38,
    69.723,      72.63,     74.9216,    78.96,    79.904,         83.798,
    85.4678,     87.62,     88.90585,   91.224,   92.90638,       95.96,
    98.0,        101.07,    102.9055,   106.42,   107.8682,       112.411,
    114.818,     118.71,    121.76,     127.60,   126.90447,      131.293,
    132.9054519, 137.327,   138.90547,  140.116,  140.90765,      144.242,
    145.0,       150.36,    151.964,    157.25,   158.92535,      162.50,
    164.93032,   167.259,   168.93421,  173.054,  174.9668,       178.49,
    180.94788,   183.84,    186.207,    190.23,   192.217,        195.084,
    196.966569,  200.59,    204.3833,   207.2,    208.9804,       209.0,
    210.0,       222.0,     223.0,      226.0,    227.0,          232.03806,
    231.03588,   238.02891, 237.0,      244.0,    243.0,          247.0,
    247.0,       251.0,     252.0,      257.0,    258.0,          258.0,
    262.0,       267.0,     268.0,      271.0,    272.0,          270.0,
    276.0,       281.0,     280.0,      285.0,    284.0,          289.0,
    288.0,       293.0,     294.0,      294.0,    5.48579883e-07, 2.01410178};
  
  /* symbols of periodic table according to LLNL (http://periodic.lanl.gov/ */
  const char atom_symb[120][3] = {
    "H",                                                                                 "He",
    "Li","Be",                                                  "B", "C", "N", "O", "F", "Ne",
    "Na","Mg",                                                  "Al","Si","P", "S", "Cl","Ar",
    "K", "Ca","Sc","Ti","V", "Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y", "Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I", "Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
                   "Hf","Ta","W", "Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
    "Fr","Ra","Ac","Th","Pa","U", "Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
                   "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Ut","Uq","Up","Uh","Us","Uo",
    "E" ,"D" }; /* E = electron and D = deuterium*/

  const char atom_name[120][25] = {
    "hydrogen","helium","lithium","beryllium","boron","carbon","nitrogen",
    "oxygen","fluorine","neon","sodium","magnesium","aluminum","silicon",
    "phosphorous","sulfur","chlorine","argon","potassium","calcium","scandium",
    "titanium","vanadium","chromium","manganese","iron","cobalt","nickel",
    "copper","zinc","gallium","germanium","arsenic","selenium","bromine",
    "krypton","rubidium","strontium","yttrium","zirconium","niobium",
    "molybdenum","technetium","ruthenium","rhodium","palladium","silver",
    "cadmium","indium","tin","antimony","tellurium","iodine","xenon","caesium",
    "barium","lanthanum","cerium","praseodymium","neodymium","promethium",
    "samarium","europium","gadolinium","terbium","dysprosium","holmium",
    "erbium","thulium","ytterbium","lutetium","hafnium","tantalum","tungsten",
    "rhenium","osmium","iridium","platinum","gold","mercury","thallium","lead",
    "bismuth","polonium","astatine","radon","francium","radium","actinium",
    "thorium","protactinium","uranium","neptunium","plutonium","americium",
    "curium","berkelium","californium","einsteinium","fermium","mendelevium",
    "nobelium","lawrencium","rutherfordium","dubnium","seaborgium","bohrium",
    "hassium","meitnerium","darmstadtium","roentgenium","copernicium",
    "ununtrium","ununquadium","ununpentium","ununhexium","ununseptium",
    "ununoctium","electron","deuterium"
  };

  chem->nPT = 120;
  chem->pt  = (atomStruct*) calloc(chem->nPT, sizeof(atomStruct));
  int i;
  for (i=0; i<chem->nPT; i++) {
    chem->pt[i].no = i+1;
    chem->pt[i].m  = atom_M[i];
    strcpy(chem->pt[i].symb, atom_symb[i]);
  }
  return 120;
}
/* -------------------------------------------------------------------------- */

chemStruct *chemInit(char *mech, parseStruct *p)
{
  int i,j,k,iSp,iR,nRxn,nMatch,rMatch,matchLen;
  FILE *pF;
  char tmpFile[256], buffer[256];
  molStruct *spec;

  /* allocate chemistry struct */
  chemStruct *chem = (chemStruct*) calloc(1, sizeof(chemStruct));
  setChemDefaults(chem);
  
  /* init basic species info */
  chemInitPeriodicTable(chem);

  /* ----- read in the chemical mechanism ----- */
  chem->compiled = 0;

  /* check the mechanism type */
  if (strlen(mech) == 0)
    strcpy(mech, "compiled");
  paramGetString(chem->mechName,"mechanism",mech,p);
  if (!strcmp(chem->mechName, "compiled")) {
    strcpy(chem->mechType, "compiled");
    chem->compiled = 1;
  }
  else if (!strcmp(chem->mechName, "wiebe")) {
    strcpy(chem->mechType, "wiebe");
    chem->wiebeMech = 1;
    fprintf(stderr, "wiebeMech: %i\n", chem->wiebeMech);
  }
  else
    strcpy(chem->mechType, "unknown");

  paramGetString(chem->reverseRxn,"reverse-reaction-rates","exact",p);
  chem->revRxnDBG = paramGetInt("reverse-reaction-debug",1,p);
  if (!strcmp(chem->reverseRxn, "fit"))
    chem->fitReverse = 1;
  else if (!strncmp(chem->reverseRxn, "exact", 3))
    chem->fitReverse = 0;

  /* try to recognize the chemical mechanism format */
  int  nMechTypes = 3;
  char mechSuffix[3][32] = {"mech", "flm", "flmech"};
  char mechTypes[3][32]  = {"chemkin", "flamemaster", "flamemaster"};
  i = 0;
  do {
    sprintf(tmpFile, "%s.%s", chem->mechName, mechSuffix[i]);
    if (checkFileExists(tmpFile)) {
      strcpy(chem->mechType, mechTypes[i]);
    }
    i++;
  } while (!strcmp(chem->mechType, "unknown") && i < nMechTypes);

  /* error checking */
  if (!strcmp(chem->mechType, "unknown")) {
    fprintf(stderr, "Error: mechanism format not recognized\n");
    fprintf(stderr, "       %s\n", chem->mechName);
    exit(2);
  }

  int nSpMech = 0, nAMech;
  char **mechSp, **mechAtoms;
  rxnStruct *r, *rB, *rF;

  if (chem->compiled) {
    /* point the functions appropriately */
    chemComputeThermoData      = ComputeThermoDataComp;
    chemComputeProductionRates = ComputeProductionRatesComp;

    /* ----- get the species information from the mechanism ----- */
    chem->nSpecies   = GetNSpecies();   // total number of species
    chem->nSpecSolve = GetNSpecs();     // number of non-steady state species
    chem->nRxn       = GetNReactions(); // number of elementary reactions
    chem->nRxnTotal  = chem->nRxn;     
    
    /* ----- allocate memory for array of names ----- */
    chem->species = (char**) calloc(chem->nSpecies, sizeof(char*));
    for (i=0; i<chem->nSpecies; i++)
      chem->species[i] = (char*) calloc(32, sizeof(char));
    
    /* get species names from the mechanism file */
    GetSpeciesNames(chem->species);
    
    /* strip all white-space from species names */
    for (i=0; i<chem->nSpecies; i++)
      stripString(chem->species[i]);
  }
  else { /* mechanism was interpreted from file directly */
    /* point the functions appropriately */
    chemComputeThermoData      = ComputeThermoDataInterp;
    chemComputeProductionRates = ComputeProductionRatesInterp;

    /* working array for reactions */
    rxnStruct *rxns;
    /* ----- read in the mechanism ----- */
    if (!strcmp(chem->mechType, "chemkin")) {
      nRxn = readChemKinMech(tmpFile, &nAMech, &mechAtoms, &nSpMech, &mechSp,
                             &rxns,chem);
    }
    else if (!strcmp(chem->mechType, "flamemaster")) {
      nRxn = readFlameMasterMech(tmpFile,&nAMech,&mechAtoms,&nSpMech,&mechSp,
                                 &rxns,chem);
    }
    else if (!strcmp(chem->mechType, "wiebe")) {
      nSpMech = initWiebeSpecies(&nAMech, &mechAtoms, &nSpMech, &mechSp, chem);
      chem->nRxn = 0; /* to be updated later to global reaction */
      chemComputeProductionRates = ComputeProductionRatesWiebe;

      /* allocate arguments to be added to chemStruct for wiebe function */
      chem->args = (void**) calloc(1, sizeof(void*));
    }

    if (nAMech > 0) {
      /* assign the elements from the table */
      chem->nElem = nAMech;
      chem->elems = (atomStruct**) calloc(nAMech, sizeof(atomStruct*));
      for (i=0; i<nAMech; i++) {
        j = getElementIndex(mechAtoms[i], chem);
        chem->elems[i] = &chem->pt[j];
      }
    }

    /* ----- process the reactions ----- */
    for (iR=0; iR<nRxn; iR++) {
      r = &rxns[iR];

      if (r == NULL) {
        fprintf(stderr, "Error: reaction pointer is NULL\n");
        exit(2);
      }

      /* check whether this reaction is only a reverse of another */
      if (r->reverse == 2) {
        if (MDUC_PRINT_PROC) fprintf(stderr, "chem: checking reverse reaction\n");
        /* first check the reaction number */
        for (i=0; i<nRxn; i++) {
          rF = &rxns[i];
          if (i != iR && r->no == rF->no) {
            if (r->nR != rF->nP && r->nP != rF->nR) {
              if (MDUC_PRINT_PROC) fprintf(stderr,"Warning: no. of reacts/prods doesn't match\n");
            }

            /* check whether the reactants/products are the same */
            rMatch = 0;
            for (j=0, nMatch=0; j<r->nR; j++) {
              for (k=0; k<rF->nP; k++) {
                if (!strcmp(r->spR[j], rF->spP[k]) && r->nuR[j] == rF->nuP[k])
                  nMatch++;
              }
            }
            if (nMatch != r->nR) {
              if (MDUC_PRINT_PROC) fprintf(stderr, "Warning: backward reaction doesn't match\n");
            }
            else
              rMatch++;

            for (j=0, nMatch=0; j<r->nP; j++) {
              for (k=0; k<rF->nR; k++) {
                if (!strcmp(r->spP[j], rF->spR[k]) && r->nuP[j] == rF->nuR[k])
                  nMatch++;
              }
            }
            if (nMatch != r->nP) {
              if (MDUC_PRINT_PROC) fprintf(stderr, "Warning: backward reaction doesn't match\n");
            }
            else
              rMatch++;

            /* reaction matches */
            if (rMatch == 2) {
              if (MDUC_PRINT_PROC) fprintf(stderr, "Using reaction %i (%i) as reverse of "
                      " reaction %i (%i)\n",r->no, iR, rF->no, i);
              r->reverse = -1;
              rF->reverse = 3;
              rF->Ar = r->A;
              rF->br = r->b;
              rF->Er = r->E;
              if (r->pressure > 0) {
                if (MDUC_PRINT_PROC) fprintf(stderr, "Warning: reverse pressure dependent "
                        "reactions not yet implemented\n");
              }
            }
          }
        } /* end loop over reactions */
      }

      /* check whether any species are repeated */
      r->nR = removeDuplicateSpecies(r->nR, &r->nuR, &r->spR);
      r->nP = removeDuplicateSpecies(r->nP, &r->nuP, &r->spP);

      /* check the fall off reactions */
      if (r->pressure > 0) {
        for (i=0, j=0; i<r->nR; i++)
          if (r->spR[i][0] == 'M')
            j = 1;
        if (j == 0) {
          r->nR = addSpecies(1.0, "M", r->nR, &r->nuR, &r->spR);
        }
        for (i=0, j=0; i<r->nP; i++)
          if (r->spP[i][0] == 'M')
            j = 1;
        if (j == 0) {
          r->nP = addSpecies(1.0, "M", r->nP, &r->nuP, &r->spP);
        }
      }

      /* get the species indices */
      r->iR = (int*) calloc(r->nR, sizeof(int));
      for (i=0; i<r->nR; i++)
        r->iR[i] = getSpeciesIndex(nSpMech, mechSp, r->spR[i]);
      r->iP = (int*) calloc(r->nP, sizeof(int));
      for (i=0; i<r->nP; i++)
        r->iP[i] = getSpeciesIndex(nSpMech, mechSp, r->spP[i]);

      /* build a total species list */
      r->nS = r->nR; /* at least as long as reactants so far */
      if (r->pressure > 0) r->nS--;
      r->iS = (int*)   calloc(r->nS, sizeof(int));
      r->sp = (char**) calloc(r->nS, sizeof(char**));
      iSp = 0;
      for (i=0; i<r->nR; i++) {
        if (strcmp(r->spR[i], "M")) {
          r->sp[iSp] = (char*) calloc(32, sizeof(char));
          strcpy(r->sp[iSp], r->spR[i]);
          r->iS[iSp] = r->iR[i];
          iSp++;
        }
      }
      iSp = r->nS;
      for (i=0; i<r->nP; i++) {
        rMatch = 0; /* reactant match */
        for (j=0; j<r->nR; j++)
          if (!strcmp(r->spP[i], r->spR[j]))
            rMatch = 1;
        if (!rMatch) {
          r->nS++;
          r->iS      = (int*)   realloc(r->iS, r->nS*sizeof(int));
          r->sp      = (char**) realloc(r->sp, r->nS*sizeof(char*));
          r->sp[iSp] = (char*)   calloc(32,          sizeof(char));
          strcpy(r->sp[iSp], r->spP[i]);
          r->iS[iSp] = r->iP[i];
          iSp++;
        }
      }

      /* build the coefficients for each species */
      r->nu = (double*) calloc(r->nS, sizeof(double));
      for (i=0; i<r->nS; i++) {
        r->nu[i] = 0.0;
        /* check reactants */
        for (j=0; j<r->nR; j++) {
          if (!strcmp(r->sp[i], r->spR[j])) {
            r->nu[i] -= r->nuR[j];
          }
        }

        /* check products */
        for (j=0; j<r->nP; j++) {
          if (!strcmp(r->sp[i], r->spP[j])) {
            r->nu[i] += r->nuP[j];
          }
        }
      }

      /* compute the reaction order */
      r->orderF[0] = r->orderF[1] = r->orderR[0] = r->orderR[1] = 0;
      for (j=0; j<r->nR; j++)
        r->orderF[1] += r->nuR[j];
      for (j=0; j<r->nP; j++)
        r->orderR[1] += r->nuP[j];
      if (r->pressure > 1) {
        /* if in fall-off region, high pressure limit is zero-order in M */
        r->orderF[0] = r->orderF[1];
        r->orderF[1]--;
        r->orderR[0] = r->orderR[1];
        r->orderR[1]--;
      }

      /* assign the enhancement factors */
      if (r->pressure > 0) {
        r->nSpecies = nSpMech;
        r->concF = (double*) calloc(nSpMech, sizeof(double));
        for (i=0; i<nSpMech; i++) r->concF[i] = 0.0;

        /* use total concentration if no enhancement factor */
        if (r->nEnhance == 0) {
          for (i=0; i<nSpMech; i++) r->concF[i] = 1.0;
        }

        /* check whether to include all in concentration */
        for (i=0; i<r->nEnhance; i++) {
          if (!strcmp(r->spEnh[i], "OTHER")) {
            for (i=0; i<nSpMech; i++) r->concF[i] = 1.0;
          }
        }

        /* find the rest of the species*/
        for (i=0; i<r->nEnhance; i++) {
          j = getSpeciesIndex(nSpMech, mechSp, r->spEnh[i]);
          if (j >= 0)
            r->concF[j] = r->factEnh[i];
          else if (strcmp(r->spEnh[i], "OTHER"))
            if (MDUC_PRINT_PROC) fprintf(stderr,"Warning: enhancement factor species %s not found\n",
                    r->spEnh[i]);
        }
      }
    } /* end loop over reactions */

    /* ----- allocate and copy into chemStruct memory ----- */
    chem->nRxn = nRxn;
    /* remove any reverse reactions that were processed */
    for (iR=0; iR<nRxn; iR++) {
      if (rxns[iR].reverse == -1)
        chem->nRxn--;
    } /* end loop over reactions */
    chem->rxns = (rxnStruct*) calloc(chem->nRxn, sizeof(rxnStruct));
    for (iR=0,i=0; iR<nRxn; iR++) {
      if (rxns[iR].reverse > -1) {
        copyRxn(&chem->rxns[i], &rxns[iR]);
        i++;
      }
      destroyRxn(&rxns[iR]);
    }
    free(rxns);
    /* compute the total number of reaction pathways */
    chem->nRxnPath = chem->nRxn;
    for (iR=0; iR<chem->nRxn; iR++) {
      if (chem->rxns[iR].reverse > 0)
        chem->nRxnPath++;
    }
    if (MDUC_PRINT_PROC) fprintf(stderr, "  Total number of reactions: %i\n", chem->nRxn);
    if (MDUC_PRINT_PROC) fprintf(stderr, "  Total number of pathways:  %i\n", chem->nRxnPath);

    /* ----- convert units ----- */
    double Aconv,Econv;
    
    for (iR=0; iR<chem->nRxn; iR++) {
      r = &chem->rxns[iR];
      /* units:
           n is order of reaction
           b is tempeature exponent 
           A has [1/(s*conc^n*K^b)]: cgs -> [cm^(3(n-1)/(s*mol^(n-1)*K^b]
                                     SI  -> [ m^(3(n-1)/(s*mol^(n-1)*K^b]
           E has [energy/mole]:      SI  -> [J/mol]
                                     Br  -> [cal/mol]
           right now, internal units are SI except use kmol
      */
      /* TODO: right now assumes input
               E -> kJ/mol  
               A -> cm, s, mol, K
           SI:
               E -> J/mol
               A -> m, s, mol, K
           needs:
               E -> J/kmol
               A -> m, s, kmol, K
      */
      /* 1e-6 m^3/cm^3 * 1e3 mol/kmol */
      Aconv = Econv = 1.0;
      strcpy(buffer,chem->unitsE);
      for (i=0; i<strlen(buffer); i++) buffer[i] = tolower(buffer[i]);
      /* check numerator */
      if (strstr(buffer,"kj"))
        Econv *= 1.0e3;
      else if (strstr(buffer,"kcal"))
        Econv *= 4.184e3;
      else if (strstr(buffer,"cal"))
        Econv *= 4.184;
      /* denominator */
      if (strstr(buffer,"kmol"))
        Econv *= 1.0;
      else
        Econv *= 1.0e3;

      strcpy(buffer,chem->unitsA);
      for (i=0; i<strlen(buffer); i++) buffer[i] = tolower(buffer[i]);
      /* check concentration units: assumes if not cm then m, and if not kmol
         then mol */
      if (strstr(buffer,"cm"))
        Aconv *= 1.0e-6;
      if (strstr(buffer,"kmol"))
        Aconv *= 1.0;
      else
        Aconv *= 1.0e3;
      r->A      *= pow(Aconv, r->orderF[1]-1);
      r->Ar     *= pow(Aconv, r->orderR[1]-1);
      r->A_low  *= pow(Aconv, r->orderF[0]-1);
      r->Ar_low *= pow(Aconv, r->orderR[0]-1);
      r->E      *= Econv;
      r->Er     *= Econv;
      r->E_low  *= Econv;
      r->Er_low *= Econv;
    } /* end loop over reactions */

    /* allocate the species */
    chem->nSpecies   = nSpMech;
    chem->nSpecSolve = nSpMech; /* TODO: enable steady state species */
    chem->species = (char**)     calloc(chem->nSpecies, sizeof(char*));
    chem->specs   = (molStruct*) calloc(chem->nSpecies, sizeof(molStruct));
    for (i=0; i<chem->nSpecies; i++) {
      spec = &chem->specs[i];
      chem->species[i] = (char*) calloc(32, sizeof(char));
      strcpy(chem->species[i], mechSp[i]);
      /* assign the structs as well, ensuring same order */
      strcpy(chem->specs[i].name, chem->species[i]);
      
      /* determine which reactions contain this species */
      spec->nRxn  = 0;
      spec->iRxn  = (int*)    calloc(1, sizeof(int));
      spec->nuRxn = (double*) calloc(1, sizeof(double));
      for (iR=0; iR<chem->nRxn; iR++) {
        r = &chem->rxns[iR];
        /* check reactants */
        for (k=0; k<r->nR; k++) {
          if (!strcmp(spec->name, r->spR[k])) {
            spec->nRxn++;
            spec->iRxn  = (int*) realloc(spec->iRxn, spec->nRxn*sizeof(int));
            spec->nuRxn = (double*) realloc(spec->nuRxn,
                                            spec->nRxn*sizeof(double));
            spec->iRxn[spec->nRxn-1]  = iR;
            spec->nuRxn[spec->nRxn-1] = -fabs(r->nuR[k]);
          }
        } /* end loop over reactants */

        /* check products */
        for (k=0; k<r->nP; k++) {
          if (!strcmp(spec->name, r->spP[k])) {
            spec->nRxn++;
            spec->iRxn  = (int*) realloc(spec->iRxn, spec->nRxn*sizeof(int));
            spec->nuRxn = (double*) realloc(spec->nuRxn,
                                            spec->nRxn*sizeof(double));
            spec->iRxn[spec->nRxn-1]  = iR;
            spec->nuRxn[spec->nRxn-1] = fabs(r->nuP[k]);
          }
        } /* end loop over products */
      } /* end loop over reactions */
    } /* end loop over species */
    
    if (chem->nRxn == 0 && !strcmp(chem->mechType, "wiebe")) {
      
    }

    /* summarize the species */
    int line, perLine = 8;
    for (i=0; i<chem->nSpecies; i++) {
      spec = &chem->specs[i];
      fprintf(stderr, "%-20s in %3i reactions\n", chem->species[i], spec->nRxn);
      for (j=0, line=1; j<spec->nRxn; j++,line++) {
        fprintf(stderr, " %4i:%+.1f", spec->iRxn[j]+1, spec->nuRxn[j]);
        if (line > 0 && line%perLine == 0)
          if (MDUC_PRINT_PROC) fprintf(stderr, "\n");
      }
      fprintf(stderr, "\n");
    } /* end loop over species */
  }

  /* ----- store location of major species in the mechanism ----- */
  chem->sO2  = getSpeciesIndex(chem->nSpecies, chem->species, "O2");
  chem->sN2  = getSpeciesIndex(chem->nSpecies, chem->species, "N2");
  chem->sCO2 = getSpeciesIndex(chem->nSpecies, chem->species, "CO2");
  chem->sH2O = getSpeciesIndex(chem->nSpecies, chem->species, "H2O");
  chem->sCO  = getSpeciesIndex(chem->nSpecies, chem->species, "CO");
  chem->sH2  = getSpeciesIndex(chem->nSpecies, chem->species, "H2");
  chem->radiation = paramGetInt("radiation",0,p);
  if (chem->radiation) {
     chem->sCH4  = getSpeciesIndex(chem->nSpecies, chem->species, "CH4");
     chem->TuRad = paramGetReal("ambient-temperature",298.0,p);
     chem->nRad = 4;
     chem->iRad = (int*) calloc(chem->nRad,   sizeof(int));
     chem->iRad[0] = chem->sH2O;
     chem->iRad[1] = chem->sCO2;
     chem->iRad[2] = chem->sCH4;
     chem->iRad[3] = chem->sCO;
     chem->pRad    = (double*) calloc(chem->nRad,   sizeof(double));
     chem->alpha   = (double*) calloc(chem->nRad,   sizeof(double));
  }
  
  /* read fuel composition from input file */
  paramGetString(chem->igniteFuel,"ignite-fuel","",p);
  chem->iFuel  = getSpeciesIndex(chem->nSpecies, chem->species, chem->igniteFuel);
  /* TODO: get those information from input file */
  chem->fuelC = paramGetReal("ignite-fuel-C",0.0,p);
  chem->fuelH = paramGetReal("ignite-fuel-H",2.0,p);
  chem->fuelO = paramGetReal("ignite-fuel-O",0.0,p);

  /* radicals */
  chem->sOH  = getSpeciesIndex(chem->nSpecies, chem->species, "OH");
  chem->sH   = getSpeciesIndex(chem->nSpecies, chem->species, "H");
  chem->sO   = getSpeciesIndex(chem->nSpecies, chem->species, "O");
  
  /* ----- allocate memory for thermo properties ----- */
  chem->W    = (double*) calloc(chem->nSpecies,   sizeof(double));
  chem->Wi   = (double*) calloc(chem->nSpecies,   sizeof(double));
  chem->cp   = (double*) calloc(chem->nSpecies,   sizeof(double));
  chem->h    = (double*) calloc(chem->nSpecies,   sizeof(double));
  chem->cdot = (double*) calloc(chem->nSpecSolve, sizeof(double));
  chem->c    = (double*) calloc(chem->nSpecies,   sizeof(double));
  chem->K    = (double*) calloc(chem->nRxn,       sizeof(double));
  chem->w    = (double*) calloc(chem->nRxn,       sizeof(double));
  chem->M    = (double*) calloc(chem->nRxn,       sizeof(double));
  chem->kF   = (double*) calloc(chem->nRxn,       sizeof(double));
  chem->kR   = (double*) calloc(chem->nRxn,       sizeof(double));
  chem->wF   = (double*) calloc(chem->nRxn,       sizeof(double));
  chem->wR   = (double*) calloc(chem->nRxn,       sizeof(double));

  /* ----- allocate memory for transport properties ----- */
  chem->transport = (int*)    calloc(chem->nSpecies, sizeof(double));
  chem->eps_k     = (double*) calloc(chem->nSpecies, sizeof(double));
  chem->sigma     = (double*) calloc(chem->nSpecies, sizeof(double));
  chem->coeffMu   = (double*) calloc(chem->nSpecies, sizeof(double));
  chem->coeffD    = (double*) calloc(chem->nSpecies, sizeof(double));

  /* ----- allocate additional memory ----- */
  chem->stoichX = (double*) calloc(chem->nSpecies, sizeof(double));

  if (chem->compiled) {
    /* store the molecular weights of each species from the mechanism */
    GetMolarMass(chem->W);
    for (i=0;i<chem->nSpecies; i++) chem->Wi[i] = 1.0 / chem->W[i]; 

     
  }
  else {
    /* ----- read in the thermo data ----- */
    int nTh, *nAtoms = NULL, **nAt = NULL;
    char **thSp = NULL, **thInfo = NULL, ***thAtomID = NULL, *thPhase = NULL;
    double Trange[3], **spT = NULL, **hCoeff = NULL;
    
    sprintf(tmpFile, "%s.therm", chem->mechName);
    nTh = readChemKinThermo(tmpFile, Trange, &thSp, &thInfo, &thPhase, &spT,
                            &nAtoms, &nAt, &thAtomID, &hCoeff);
    
    if (nTh != chem->nSpecies) {
      if (MDUC_PRINT_PROC) fprintf(stderr, "Warning: %i species in thermo file / %i in mechanism\n",
              nTh, chem->nSpecies);
    }
    
    for (i=0; i<3; i++) chem->Trange[i] = Trange[i];
    
    /* strip the whitespace from the species names read in */
    for (i=0; i<nTh; i++) {
      stripString(thSp[i]);
      if      (thSp[i][0] == '1') thSp[i][0] = 'P';
      else if (thSp[i][0] == '2') thSp[i][0] = 'S';
      else if (thSp[i][0] == '3') thSp[i][0] = 'T';
      for (j=0; j<nAtoms[i]; j++)
        stripString(thAtomID[i][j]);
    }
    
    /* initialize the thermo data for checking later */
    for (i=0; i<chem->nSpecies; i++)
      for (j=0; j<7; j++)
        chem->specs[i].therm[0][j] = chem->specs[i].therm[1][j] = -1.0E30;
    
    /* assign the thermo data appropriately */
    /* char **matches = (char**) calloc(1, sizeof(char*)); */
    for (i=0; i<chem->nSpecies; i++) {
      j = getSpeciesIndex(nTh, thSp, chem->species[i]);
      if (j < 0) {
        /* test whether it exists partial */
        for (k=0, nMatch=0, matchLen=0; k<nTh; k++) {
          if (!strncmp(chem->species[i], thSp[k], strlen(thSp[k]))) {
            nMatch++;
            if (strlen(thSp[k]) > matchLen) {
              matchLen = strlen(thSp[nMatch-1]);
              j = k;
            }
          }
        }
        
        if (nMatch > 1) {
          if (MDUC_PRINT_PROC) fprintf(stderr, "Warning: %i potential thermo matches for %18s --> ",
                  nMatch, chem->species[i]);
          if (MDUC_PRINT_PROC) fprintf(stderr, "taking longest species: %s\n", thSp[j]);
        }
      
        if (j < 0) { /* still not found */
          fprintf(stderr, "Error: thermo data for %s not found (%i)\n",
                  chem->species[i], i);
          exit(2);
        }
        
      }
      
      /* ----- copy the info from the species ----- */
      strcpy(chem->specs[i].info, thInfo[j]);
      /* elements */
      chem->specs[i].nAtoms = nAtoms[j];
      chem->specs[i].nA   = (int*)         calloc(nAtoms[j],sizeof(int));
      chem->specs[i].atom = (atomStruct**) calloc(nAtoms[j],sizeof(atomStruct));
      chem->W[i] = 0.0;
      for (k=0; k<nAtoms[j]; k++) {
        chem->specs[i].nA[k] = nAt[j][k];
        iSp = getElementIndex(thAtomID[j][k], chem);
        chem->specs[i].atom[k] = &chem->pt[iSp];
        
        /* compute the molecular weights */
        chem->W[i] += nAt[j][k]*chem->pt[iSp].m;
      }
      chem->Wi[i] = 1.0 / chem->W[i]; 
      
      /* phase */
      chem->specs[i].phase = thPhase[j];
      
      /* temperature ranges */
      if (spT[j][0] == 300.0)
        chem->specs[i].Trange[0] = 298; /* low temp */
      else
        chem->specs[i].Trange[0] = spT[j][0]; /* low temp */
      if (spT[j][2] > 0.0)
        chem->specs[i].Trange[1] = spT[j][2]; /* mid temp */
      else
        chem->specs[i].Trange[1] = Trange[1];
      chem->specs[i].Trange[2] = spT[j][1]; /* high temp */
      
      /* store the raw coefficients for temperature ranges */
      /* high temperature */
      for (k=0; k<7; k++)  chem->specs[i].therm[1][k]   = hCoeff[j][k];
      /* low temperature */
      for (k=7; k<14; k++) chem->specs[i].therm[0][k-7] = hCoeff[j][k];

      /* pre-compute coefficients for different property polynomials */
      /* specific heat coefficients */
      for (k=0; k<5; k++) {
        chem->specs[i].aCp[0][k] = chem->specs[i].therm[0][k];
        chem->specs[i].aCp[1][k] = chem->specs[i].therm[1][k];
      }
      /* enthalpy coefficients */
      for (k=0; k<2; k++) {
        chem->specs[i].aH[k][0] = chem->specs[i].therm[k][5];
        chem->specs[i].aH[k][1] = chem->specs[i].therm[k][0];
        chem->specs[i].aH[k][2] = chem->specs[i].therm[k][1]/2.0;
        chem->specs[i].aH[k][3] = chem->specs[i].therm[k][2]/3.0;
        chem->specs[i].aH[k][4] = chem->specs[i].therm[k][3]/4.0;
        chem->specs[i].aH[k][5] = chem->specs[i].therm[k][4]/5.0;
      }
      /* entropy coefficients */
      for (k=0; k<2; k++) {
        chem->specs[i].aS[k][0] = chem->specs[i].therm[k][0]; /* for ln(T) */
        chem->specs[i].aS[k][1] = chem->specs[i].therm[k][6];
        chem->specs[i].aS[k][2] = chem->specs[i].therm[k][1];
        chem->specs[i].aS[k][3] = chem->specs[i].therm[k][2]/2.0;
        chem->specs[i].aS[k][4] = chem->specs[i].therm[k][3]/3.0;
        chem->specs[i].aS[k][5] = chem->specs[i].therm[k][4]/4.0;
      }
    } /* finished species */

    /* compute fits for backward rates */
    if (chem->fitReverse) {
      for (iR=0; iR<chem->nRxn; iR++) {
        r = &chem->rxns[iR];
        if (r->reverse == 1) {
          chemFitReverseRxnRate(r, chem);
          r->reverse = 3;
        }
      } /* end loop over reactions (reverse rate calc) */
    }
    
    /* free temporary memory */
    for (i=0; i<nTh; i++) {
      for (j=0; j<nAtoms[i]; j++) free(thAtomID[i][j]);
      free(thSp[i]); free(thInfo[i]); free(thAtomID[i]); free(nAt[i]);
      free(spT[i]); free(hCoeff[i]);
    }
    free(thSp); free(thInfo); free(thAtomID); free(nAt);
    free(spT); free(hCoeff);
    free(thPhase); free(nAtoms);

    /* TODO: add step here for processing atoms in reactions? */
  }

  /* ----- read in all the transport properties ----- */
  chem->readTrans = paramGetInt("read-transport-data",0,p);
  paramGetString(chem->transFile,"transport-data-file","",p);
  if (strlen(chem->transFile) > 0)
    chem->readTrans = 1;

  if (chem->readTrans) {
    int    nTrans, iMatch;
    char   **transSp = NULL;
    double *eps_k, *sigma;
    eps_k = sigma = NULL;
    if (strlen(chem->transFile) == 0)
      sprintf(chem->transFile, "%s.trans", chem->mechName);

    nTrans = readChemKinTrans(chem->transFile, &transSp, &eps_k, &sigma);

    /* strip the whitespace from the species names read in */
    for (i=0; i<nTrans; i++) {
      stripString(transSp[i]);
      if      (transSp[i][0] == '1') transSp[i][0] = 'P';
      else if (transSp[i][0] == '2') transSp[i][0] = 'S';
      else if (transSp[i][0] == '3') transSp[i][0] = 'T';
    }
    
    if (nTrans != chem->nSpecies) {
      if (MDUC_PRINT_PROC) fprintf(stderr, "Warning: number of species in transport file differs\n");
      if (MDUC_PRINT_PROC) fprintf(stderr, "   mechanism species: %i\n", chem->nSpecies);
      if (MDUC_PRINT_PROC) fprintf(stderr, "   transport species: %i\n", nTrans);
    }
    
    for (i=0; i<chem->nSpecies; i++) {
      iSp = getSpeciesIndex(nTrans, transSp, chem->species[i]);
      if (iSp < 0) {
        if (MDUC_PRINT_PROC) fprintf(stderr,"Warning: species not found in transport file: %32s ",
                chem->species[i]);
        matchLen = strlen(chem->species[i]);

        /* find the longest matched species */
        iMatch = -1;
        for (j=0; j<nTrans; j++) {
          if (!strncmp(chem->species[i], transSp[j], strlen(transSp[j]))) {
            if (iMatch >= 0) {
              if (strlen(transSp[j]) > strlen(transSp[iMatch]))
                iMatch = j;
            }
            else
              iMatch = j;
          }
        }

        if (iMatch < 0) {
          fprintf(stderr,"Error: no transport data found for %s\n",
                  chem->species[i]);
          exit(2);
        }
        else {
          if (MDUC_PRINT_PROC) fprintf(stderr,"-> using transport data found for %s\n",
                  transSp[iMatch]);
          chem->transport[i] = 1;
          chem->eps_k[i] = eps_k[iMatch];
          chem->sigma[i] = sigma[iMatch];

          /* compute the constant portion of the viscosity relation */
          chem->coeffMu[i] = calcViscCoeff(chem->W[i], chem->sigma[i]);
        }
      }
      else {
        chem->transport[i] = 1;
        chem->eps_k[i] = eps_k[iSp];
        chem->sigma[i] = sigma[iSp];
      
        /* compute the constant portion of the viscosity relation */
        chem->coeffMu[i] = calcViscCoeff(chem->W[i], chem->sigma[i]);
      }
    }

    for (i=0; i<nTrans; i++) { free(transSp[i]); };
    free(transSp); free(eps_k); free(sigma);

    /* compute the properties for computing mixture viscosity */
    chem->mu05   = (double*) calloc(chem->nSpecies, sizeof(double));
    chem->mu05i   = (double*) calloc(chem->nSpecies, sizeof(double));

    chem->WjOverWi   = (double**) calloc(chem->nSpecies, sizeof(double*));
    chem->Wij_D      = (double**) calloc(chem->nSpecies, sizeof(double*));
    chem->phiMixMu_a = (double**) calloc(chem->nSpecies, sizeof(double*));
    for (i=0; i<chem->nSpecies; i++) {
      chem->WjOverWi[i]   = (double*) calloc(chem->nSpecies, sizeof(double));
      chem->Wij_D[i]      = (double*) calloc(chem->nSpecies, sizeof(double));
      chem->phiMixMu_a[i] = (double*) calloc(chem->nSpecies, sizeof(double));
      for (j=0; j<chem->nSpecies; j++) {
        /* ratio of molecular weights Wj/Wi */
        chem->WjOverWi[i][j]   = pow(chem->W[j]/chem->W[i], 0.25);
        /* reduced molecular mass for species pair (for D_ij evaluation) */
        chem->Wij_D[i][j]      = chem->W[i]*chem->W[j]/(chem->W[i]+chem->W[j]);
        /* leading coefficient of mixture viscosity relation
              1/sqrt(8)*(1 + Wi/Wj)^(-1/2) */
        chem->phiMixMu_a[i][j] = 0.35355339*pow(chem->W[i]/chem->W[j]+1.0,-0.5);
      }
    }
    /* /\* DEBUG */
    /* writeChemKinTrans("debug_CK.trans",chem); */
    /* GUBED *\/ */
  }

  /* check if there is a global reaction specified */
  chem->stoich = paramGetReaction("global-reaction",p);
  if (chem->stoich) {
    /* index the species */
    chem->stoich->iS = (int*) calloc(chem->stoich->nS, sizeof(int));
    chem->stoich->iR = (int*) calloc(chem->stoich->nR, sizeof(int));
    chem->stoich->iP = (int*) calloc(chem->stoich->nP, sizeof(int));
    for (i=0; i<chem->stoich->nS; i++)
      chem->stoich->iS[i] = getSpeciesIndex(chem->nSpecies,chem->species,
                                            chem->stoich->sp[i]);
    for (i=0; i<chem->stoich->nR; i++)
      chem->stoich->iR[i] = getSpeciesIndex(chem->nSpecies,chem->species,
                                            chem->stoich->spR[i]);
    for (i=0; i<chem->stoich->nP; i++)
      chem->stoich->iP[i] = getSpeciesIndex(chem->nSpecies,chem->species,
                                            chem->stoich->spP[i]);
  }

  /* try to process the fuel information if present
     otherwise, will try to initialize again from stream in initIO */
  chem->nFuel = paramGetStringArray(&chem->fuel,"fuel-species",p);
  if (chem->nFuel > 0) {
    paramGetRealArray(&chem->fuelX,"fuel-mole-frac*",p);
    paramGetRealArray(&chem->fuelY,"fuel-mass-frac*",p);
    chem->fuelLHV = paramGetReal("fuel-lower-heat*",-1e30,p);
    initChemFuel(chem);
  }
  else { /* assume that fuel side is the fuel */
    int iP = paramGetIndex("FLAMELET",p);
    bndStruct *fuelBC = paramGetBC("FuelSide1",p->n[iP],p->vals[iP]);
    if (fuelBC) {
      chem->nFuel = fuelBC->nSp;
      chem->fuel  = (char**) calloc(chem->nFuel, sizeof(char*));
      /* check that the species mass fractions are defined */
      int Ydef = 0;
      for (i=0; i<fuelBC->nSp; i++) {
        if (fuelBC->Y[i] > 0.0) Ydef++;
        chem->fuel[i] = (char*) calloc(32, sizeof(char));
        strcpy(chem->fuel[i],fuelBC->sp[i]);
      }
      
      if (Ydef < fuelBC->nSp) {
        Ydef = 0;
        for (i=0; i<fuelBC->nSp; i++)
          if (fuelBC->X[i] > 0.0) Ydef++;
        if (Ydef < fuelBC->nSp) {
          if (MDUC_PRINT_PROC) fprintf(stderr, "Warning: fuel species not identified\n");
        }
        else {
          chem->fuelX = (double*) calloc(chem->nFuel, sizeof(double));
          for (i=0; i<fuelBC->nSp; i++)
            chem->fuelX[i] = fuelBC->X[i];
        }
      }
      else {
        chem->fuelY = (double*) calloc(chem->nFuel, sizeof(double));
        for (i=0; i<fuelBC->nSp; i++)
          chem->fuelY[i] = fuelBC->Y[i];
      }
      initChemFuel(chem);
    }
    else
      if (MDUC_PRINT_PROC) fprintf(stderr,"no fuel to initialize with\n");
  }

  /* check to see if the soot model is on */
  paramGetString(chem->sootModel,"soot-model","",p);
  /* TODO: move this initialization into initSoot, and pass the parser */

  if (!strcmp(chem->sootModel, "HMOM"))
    chem->soot = 1;

  if (chem->soot > 0) {
	chem->sootMode = paramGetInt("soot-mode",0,p);

    if (chem->sootMode == 2) {
       fprintf(stderr, "Error: Soot mode 2 currently not supported\n");
       exit(2);
    }

    chem->sootOrder = paramGetInt("soot-moment-order",1,p);
    chem->sootTolerance = paramGetReal("soot-tolerance",1.0E-15,p);
    chem->nSootSrcOn =
      paramGetStringArray(&chem->sootSrcSelect,"soot-active-source*",p);
    chem->nSootMomOn = 
      paramGetStringArray(&chem->sootMomSelect,"soot-active-mom*",p);
    chem->nSootPAH = 
      paramGetStringArray(&chem->sootPAHSpec, "soot-PAH-species*",p);
    initSoot(chem);
  }

  /* allocate arrays for reaction pathway analysis */
  /* TODO: Maybe one should only allocate if reaction rates are actually recoreded 
     in order to save memory */
  chem->w_tmp	= (double*) calloc(chem->nRxnTotal,  sizeof(double));
  chem->w_pre	= (double*) calloc(chem->nRxnTotal,  sizeof(double));
  chem->w_int	= (double*) calloc(chem->nRxnTotal,  sizeof(double));
  chem->w_int_pre	= (double*) calloc(chem->nRxnTotal, sizeof(double));

  /* ----- write summary ----- */
  if (MDUC_PRINT_PROC) {
    fprintf(stderr,"------------- Chemistry Summary -------------\n");
    fprintf(stderr," Mechanism:   %s (%s)\n",chem->mechName,chem->mechType);
    fprintf(stderr,"  Species:    %i\n",chem->nSpecies);
    fprintf(stderr,"  Reactions:  %i\n",chem->nRxn);
    fprintf(stderr,"  Rev. Rxns:  %i\n",chem->nRxnPath-chem->nRxn);
    fprintf(stderr,"  Pathways:   %i\n",chem->nRxnPath);
    fprintf(stderr,"  Rev. Rates: %s\n",chem->reverseRxn);

    if (chem->nFuel > 0) {
      fprintf(stderr," ---------------------------------------------\n");
      fprintf(stdout,"  %-20s |    Y     |    X     |\n","Fuel species");
      fprintf(stderr," ---------------------------------------------\n");
      for (i=0; i<chem->nFuel; i++) {
        fprintf(stdout,"  %-20s | %8.6f | %8.6f |\n",chem->fuel[i],
                chem->fuelY[i],chem->fuelX[i]);
      }
    }
    if (chem->soot > 0) {
      fprintf(stderr," Soot Model:  %s\n",chem->sootModel);
      fprintf(stderr," Soot Order:  %i\n",chem->sootOrder);
      fprintf(stderr," Soot Mode:   %i\n",chem->sootMode);
    }
    fprintf(stderr,"----------------------------------------------\n");
  }

  return chem;
}
/* -------------------------------------------------------------------------- */
