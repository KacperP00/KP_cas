#ifndef MDUC_FLAME_STRUCT_H
#define MDUC_FLAME_STRUCT_H

struct GRADSTRUCT {
  /* ----- derivatives with respect to variable ----- */
  double *x, *y;   // first derivatives
  double *xx, *yy; // second derivatives
  double *xy;      // cross derivative
};
/* -------------------------------------------------------------------------- */
typedef struct GRADSTRUCT gradStruct;

struct METRICSTRUCT {
  /* ----- metric coefficients for gradients ----- */
  double **x, **y;
  double **xx, **yy;
};
/* -------------------------------------------------------------------------- */
typedef struct METRICSTRUCT metricStruct;

typedef struct bndStruct_ {
  /* ----- boundary condition struct ----- */
  int    type;
  char   name[32];

  /* species information */
  int     nSp;
  int    *iSp;
  char   **sp;
  double *Y;
  double *X;

  /* energy */
  double T;         /* temperature [K] */
  double H;         /* enthalpy [J/kg] */

  /* variables for liquid fuels */
  double Hvap;      /* enthalpy of vaporization [J/kg] */
  double Tref;      /* reference temperature [K] */
  double cpLiq;     /* liquid specific heat [J/kg-K] */
  double P;         /* ambient pressure [Pa] */
  double Mw;        /* molecular weight [g/mol] */
  char   sat[128];  /* saturation data file */
} bndStruct;
/* -------------------------------------------------------------------------- */

struct FLAMESTRUCT {
  /* ----- driver variables ----- */
  int    type;        /* 0: non-premixed, 1: premixed */
  int    space;       /* 0: physical, 1: mixture fraction */
  int    axisymm;     /* 0: planar, 1: axisymmetric*/

  /* ----- variables for grid properties ----- */
  int    nDim;        /* Number of Dimensions in Flamelet */
  int    nFlamelets;  /* number of members of array of flamelet structures */
  int    n1, n2;      /* number of grid points in 1, 2 */
  int    j0,jN,k0,kN; /* local starting points on grid */
  int    j0m,jNp;
  int    k0m,kNp;
  double *x1, *x2;    /* location in x1 or x2 space at index j */
  double *dx1, *dx2;  /* mesh spacing in x1 or x2: dx(j) = x(j+1)-x(j) */
  int    iX10,iX01,iX11;

  /* stretched grid parameters */
  int    nPiv1,nPiv2;    /* number of pivot points */
  double *gPiv1,*gPiv2;  /* pivot variables */
  double *gPc1,*gPc2;    /* percentage for each range */

  /* ----- boundary conditions ----- */
  bndStruct *stream[3];
  bndStruct *ox,*fuel1,*fuel2;

  /* ----- metric variables ----- */
  int          order, stw; /* spatial gradient order and stencil width */
  metricStruct grad;       /* gradient metric operators */
  
  /* ---- variables due to normalisation ----- */
  double *v1, *v2;            /* Convective velocity in x1, x2 direction  */
  double *vT;		      /* Convective velocity for temperature */
  double *vMxmyp, *vMxy; 	      /* Convective velocity for soot moments */
  double CFL;
  double *chiA, *chiB, *chiC; /* transformed scalar dissipation rates */
  /* ---- required for soot transport ---- */
  double *rhoSigma;
  double *xi;
  double *rhoDotSoot;
  double *rhoDotSpray;
  /* ----- arrays to hold solution variables (species & temperature) ----- */
  int    nVars;    /* number of variables at each mesh point */
  double *Y;       /* current solution, t(N), for all x1 and x2 */
  double *Yp;      /* current value of the derivative (Y prime) */
  double *Y_m1;    /* solution at t(N-1) for all x1 and x2 */     
  double *Y_m2;    /* solution at t(N-2) for all x1 and x2 */
  double *Yp_m1;   /* Yprime at t(N-1) */
  double *Y_pre;   /* actual solution at t(N-1) for all x1 and x2 */
  double *sumYi;
  gradStruct dY;   /* to store derivatives for explicit evaluation */
  /* ----- reference properties ----- */
  double Tmax, flameZ;
  
  /* ----- values from species mechanism ----- */
  /* ??? could just use directly from chem_struc ??? */
  int    sT, sP, sO2, sN2, sChi, sSoot, sRxn, sRxnL;
  int    sM, sMa, sMe, sMf, sMb, sTb, sDF; 
  int    sCont, sMom, sZ;
  
  /* ----- working memory for thermo properties needed in solver ----- */
  double *H;                            /* total enthalpy of the flamelet */
  double *T;                            /* temperature array */
  double *rho, *cpmix, *sumCp, *sumCp2; /* local mixture properties */
  double *R;
  double *srcH, *srcP;                  /* pressure and enthalpy source terms */
  double *srcQ;                         /* loss source term for enthalpy */
  double intSrcQ;
  double *HRrate;                       /* heat release rate*/
  gradStruct dCp;                       /* derivative of mix. specific heat */
  double RGAS;                          /* Universal Gas Constant */

  /* ----- local mixture average properties ----- */
  double *mu;                           /* viscosity */
  double *D;                            /* diffusivity */
  double *lambda;                       /* thermal conductivity */
  double *gamma;                        /* ratio of specific heats */
  double *Mw;                           /* molecular weight */
  double *LeZ_LeI;			/* Lewis Z by Lewis of species i */
  double Tmean;

  /* ----- premixed variables ----- */
  double S_L;                           /* laminar burning velocity */
  double thick;                         /* flame thickness */
  double phi;                           /* equivalence ratio */

  /* ----- viscous heat production ----- */
  double *visc_src;                      /* conditional viscous heat production */ 
  double *gPDF;                      /* global pdf */

  /* ----- soot moments ----- */
  int    nSootMom;
  double *sootMom;
  double *sootMomxmyp;
  double *sootMomSrc, *sootMomSrcPart;
  double *sootMomTrans;
  double *sootMomTransConv;
  double *sootMomTransDiff;
  /* ----- radiation ----- */
  double *radiation;
  /* ----- thermodynamic soot coefficients (turbulence-soot interaction) ----- */
  double *sqrtT;
  double *T_mu;
  double *muSqrtMw_rhoSqrtT;
  double *sgCoeff;
  double *oxCoeff;
  double *o2Coeff;
  double *dimProdRate;

  /* ----- PAH variables ----- */
  double *yPAH;
  double *srcPAHpos;
  double *srcPAHneg;

  /* debugging terms */
  double *convCoeff;

  /* Source terms */
  double *chem_src, *conv_src, *diff_src;
  double *chem_temp, *conv_temp, *diff_temp;
  double *chemSrcY, *chemSrcY_tmp;  /* chem src for all x1 and x2 */
};
/* -------------------------------------------------------------------------- */
typedef struct FLAMESTRUCT flameStruct;

#endif /* MDUC_FLAME_STRUCT_H */
