#ifndef MDUC_CHEM_STRUCT_H
#define MDUC_CHEM_STRUCT_H

struct ATOMSTRUCT {
  int    no;          /* atomic number */
  double m;           /* molar mass */
  char   symb[3];     /* atomic symbol */
};
typedef struct ATOMSTRUCT atomStruct;

struct MOLSTRUCT {
  char   name[32];    /* molecule name */
  int    nAtoms;      /* number of unique atoms */
  double W;           /* molecular weight */
  int    *nA;         /* number of each atom */
  atomStruct **atom;  /* atoms */
  int    cond;        /* > 0 if condensed */
  char   phase;       /* (g)as, (l)iquid, or (s)olid */

  /* reaction mechanism info for this species */
  int    nRxn;        /* number of reactions where present */
  int    *iRxn;       /* indices of reactions where present */
  double *nuRxn;      /* coefficient for reactions where present */

  /* thermodynamic data for this species: [0] low, [1] high */
  int    nCpCoeff[9]; /* number of coefficients in Cp fit */
  double therm[9][10]; /* raw coefficients */
  double aCp[9][7];   /* coefficients for specific heat */
  double aH[9][8];    /* coefficients for enthalpy */
  double aS[9][8];    /* coefficients for entropy */
  double Texp[9][8];  /* temperature exponents */
  int    nTrange;     /* number of temperature intervals */
  double Trange[9];   /* [0] low, [1] mid [2] high */
  double condT;       /* temperature for enthalpy of condensed species */
  char   info[8];     /* info from thermo file */
  double dHf;         /* enthalpy of formation @ 298.15K [J/mol] */
  double H0[4];       /* H0(298.15) - H0(0) [J/mol] */

  /* transport data */
  int    geom;        /* atomic structure: 0 monatomic, 1 linear, 2 non-lin */
  double eps_k;       /* Lennard-Jones well depth/Boltzman cst. [K] */
  double sigma;       /* Lennard-Jones collision diameter [m] */
  double mu;          /* dipole moment [debye] */
  double alpha;       /* polarizability [Angstrom^3] */
  double Zrot;        /* rotational relaxation collision number @ 298K */
  double coeffMu;     /* const. coeff. for viscosity rel. [?] */
  double coeffD;      /* const. coeff. for diffusivity rel. [?] */
};
typedef struct MOLSTRUCT molStruct;

struct RXNSTRUCT {
  int  no;            /* reaction number */
  char dir;           /* direction: (f)orward or (r)everse */
  int reverse;        /* reversible? */
  int orderF[2];      /* forward reaction order (0: low, 1: inf) */
  int orderR[2];      /* reverse reaction order (0: low, 1: inf) */

  /* ----- element and species information ----- */
  int nS;             /* number of species */
  int nR;             /* number of reactants */
  int nP;             /* number of products */
  int nA;             /* number of atoms in reaction */
  int   *iS;          /* indices of species */
  int   *iR;          /* indices of reactant species */
  int   *iP;          /* indices of product species */
  char **sp;          /* list of species */
  char **spR;         /* list of reactant species */
  char **spP;         /* list of product species */
  char **atoms;       /* atoms in reaction */
  double *nu;         /* mole balance of species */
  double *nuR;        /* mole number of reactant species */
  double *nuP;        /* mole number of product species */

  /* ----- rate constant parameters ----- */
  double A,A_low;     /* pre-exponential factor */
  double b,b_low;     /* temperature exponent */
  double E,E_low;     /* activation energy */
  /* reverse rates */
  double Ar,Ar_low;
  double br,br_low;
  double Er,Er_low;

  /* ----- third body reactions ----- */
  int thirdBody;      /* reaction has third body */

  /* ----- pressure dependence / auxiliary flags ----- */
  int pressure;       /* 0: independent, 1: low-pressure, 2: fall-off */
  int lindemann;      /* use Lindemann formula for pressure dependence */
  int troe;           /* use Troe's formula for pressure dependence */
  int sri;            /* use SRI formula for pressure dependence */
  double troeA[3];    /* Troe coefficients for F_c formula */
  double troeT[3];    /* Troe parameters T***, T*, T** */
  double sriA[5];     /* SRI parameters a, b, c, d, e */

  /* ----- enhancement factors ----- */
  int    nEnhance;    /* number of species with enhancement factors */
  char   **spEnh;     /* species for enhancement factors */
  double *factEnh;    /* enhancement factors */
  int    nSpecies;    /* number of total species to consider */
  double *concF;      /* concentration */
};
typedef struct RXNSTRUCT rxnStruct;

typedef struct sootStruct_ {
  int     nMom;             /* number of soot moments */
  int     nSrc;             /* number of soot source contributions */
  int     nRxn;             /* number of soot reactions */
  int     nSp;              /* number of species involved in soot rxns */
  int    *iSp;              /* indices of species involved in soot rxns */
  char   **rxnSc;           /* species involved in soot reactions */
  double *rxnRates;         /* reaction rates for soot [mol/(m^3-s)] */
  double *rxnConc;          /* pointer to relevant soot rxn species */
  int    nSrcOn;            /* number of active soot sources */
  int    *srcOn;            /* on/off control for specific soot sources */
  char   **srcNames;
  char   **srcSelect;       
} sootStruct;

typedef struct chemStruct_ {
  /* ----- drivers ----- */
  int  compiled;
  int  readThermo;
  int  readTrans;
  int  readMech;
  int  radiation;
  int  soot;          /* flag for soot model: 0 none, 1 HMOM */
  int  sootMode;      /* 0: soot moments from external
                         1: soot moments part of solution implicit
                         2: soot moments updated explicitly */
  int  sootOrder;     /* order of moments to consider in soot model */
  double sootTolerance;     /* order of moments to consider in soot model */
  char sootModel[32]; /* name of soot model */
  char reverseRxn[8]; /* reverse reaction rate treatment */ 
  int  fitReverse;    /* pre-compute fits for reverse reaction rates */
  int  revRxnDBG;

  /* ----- mechanism info ----- */
  char mechName[128]; /* name of chemical mechanism used */
  char mechType[32];  /* format of the mechanism */
  char transFile[256];
  int  wiebeMech;     /* use simple Wiebe function for chemistry model */

  /* ----- data from the chemical mechanism ----- */
  int  nSpecies;      /* total number of species */
  int  nSpecSolve;    /* number of steady state species */
  int  nRxn;          /* number of reactions */
  int  nRxnPath;      /* number of reaction pathways (forward + reverse) */
  int  nElem;         /* number of elements */
  char **species;     /* array of string containing all species names */
  molStruct  *specs;  /* TEST: struct for each species */
  rxnStruct  *rxns;   /* TEST: struct for each reaction */
  atomStruct **elems; /* TEST: struct for each element */
  char unitsA[64];    /* units for leading factor in Arrhenius rate */
  char unitsE[32];    /* units for activation energy */
  double Trange[3];   /* global thermo temperature ranges */

  /* ----- store location of major species in the mechanism ----- */
  int sO2, sN2, sCO2, sH2O, sCO, sH2, sC2H2, sCH4; //, sA4;
  int sOH, sH, sO;
  int *sPAH;
  int iFuel,nFuel, *sFuel;
  char **fuel;
  char igniteFuel[256];
  double fuelC, fuelH, fuelO, fuelMw, fuelLHV, *fuelY, *fuelX;
  rxnStruct *stoich;
  double *stoichX, stoichFA;
  int nRad; /* number of species considered for radiation */
  int *iRad; /* number of species considered for radiation */
  double *pRad; /* partial pressure of radiating species */ 
  double *alpha; /* absorbtion coefficients */ 
  double TuRad; /* ambient temperature */
  /* ----- arrays needed by chemical mechanism ----- */
  double *W;                    /* molecular weight of all species */
  double *Wi;           /* molecular weight */
  double *cp, *h;               /* specific heat and enthalpy */
  double *cdot, *c, *K, *w, *M; /* production rate variables for mechanism */
  double *kF, *kR;              /* forward/reverse rate coefficients */
  double *wF, *wR;              /* forward/reverse reaction rates */
  double *w_tmp,*w_pre, *w_int,*w_int_pre;        /* store reaction rates for flux analysis */

  /* ----- species transport properties ----- */
  int    *transport; 
  int    *atomGeom;             /* atomic structure: 0 mono, 1 linear, 2 nonl */
  double *eps_k;                /* Lennard-Jones well depth/Boltzman cst. [K] */
  double *sigma;                /* Lennard-Jones collision diameter [m] */
  double *muDipole;             /* dipole moment [debye] */
  double *alphaPolar;           /* polarizability in [Angstrom^3] */
  double *Zrot;                 /* rotational relaxation collision number */
  double *coeffMu;              /* const. coeff. for viscosity rel. [?] */
  double *coeffD;               /* const. coeff. for diffusivity rel. [?] */
  double **phiMixMu_a;          /* leading coefficient for mix. viscosity */
  double **WjOverWi;            /* ratios of mol. weight for mix. viscosity  */
  double **Wij_D;               /* reduced mol. weight for diffusion */
  double *mu05;			/* square root of viscosity */
  double *mu05i;		/* inverse square root of viscosity */

  /* ----- soot properties ----- */
  int     nSootMom;             /* number of soot moments */
  int     nSootSrc;             /* number of soot source contributions */
  int     nSootRxn;             /* number of soot reactions */
  int     nSootSp;              /* number of species involved in soot rxns */
  int     nSootPAH;                 /* number of PAH molecules  */
  int    *iSootSp;              /* indices of species involved in soot rxns */
  char   **sootRxnSpec;         /* species involved in soot reactions */
  double *sootRxnRates;         /* reaction rates for soot [mol/(m^3-s)] */
  double *sootRxnConc;          /* pointer to relevant soot rxn species */
  int    nSootSrcOn;
  int    nSootMomOn;
  char   **sootSrcNames;
  char   **sootMomNames;
  char   **sootMomAbbrv;
  char   **sootSrcSelect;       
  char   **sootMomSelect;
  char   **sootPAHSpec;          /* PAH species involved in soot formation */
  int    *sootSrcOn;            /* on/off control for specific soot sources */
  int    *sootMomOn;            /* on/off control for specific soot moments */

  /* ----- for linking other arguments in the mechanism ----- */
  void **args;

  /* ----- periodic table ----- */
  int        nPT;
  atomStruct *pt;

  int     nRxnTotal;             /* number of gas phase + soot reactions */
  int     sSoot;		 /* offset for soot reactions */

  /* ----- some constants ----- */
  double logPoverR;
} chemStruct;

#endif /* MDUC_CHEM_STRUCT_H */
