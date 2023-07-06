#include "accessMDUC.h"
#include "udfScalarDiss.h"
#include "udfTools.h"

extern udfStruct udf;

void UDFinitScalarDiss(udfStruct *udf)
{
  /* set defaults if not defined in scheme file */
  if (udf->chiForm1D == NULL) {
    udf->chiForm1D = (char*) calloc(50, sizeof(char));
    strcpy(udf->chiForm1D, "log");
  }
  if (udf->chiForm2D == NULL) {
    udf->chiForm2D = (char*) calloc(50, sizeof(char));
    strcpy(udf->chiForm2D, "log");
  }
  if (udf->chiMode1D == NULL) {
    udf->chiMode1D = (char*) calloc(50, sizeof(char));
    strcpy(udf->chiMode1D, "distribution");
  }
  if (udf->chiMode2D == NULL) {
    udf->chiMode2D = (char*) calloc(50, sizeof(char));
    strcpy(udf->chiMode2D, "distribution");
  }
  if (udf->chiNorm1D == NULL) {
    udf->chiNorm1D = (char*) calloc(50, sizeof(char));
    strcpy(udf->chiNorm1D, "global-zero");
  }
  if (udf->chiNorm2D == NULL) {
    udf->chiNorm2D = (char*) calloc(50, sizeof(char));
    strcpy(udf->chiNorm2D, "global-zero");
  }

  /* convert string input to appropriate integer values */
  udf->chiForm[1] = UDFsetChiForm(udf->chiForm1D);
  udf->chiForm[2] = UDFsetChiForm(udf->chiForm2D);

  udf->chiMode[1] = UDFsetChiMode(udf->chiMode1D);
  udf->chiMode[2] = UDFsetChiMode(udf->chiMode2D);

  udf->chiNorm[1] = UDFsetChiScaling(udf->chiNorm1D);
  udf->chiNorm[2] = UDFsetChiScaling(udf->chiNorm2D);
}
/* ------------------------------------------------------------------------- */

int UDFsetChiMode(char *mode)
{
  if (!strncmp(mode, "distribution", 4))
    return 1;
  else if (!strcmp(mode, "bin"))
    return 2;
  else if (!strcmp(mode, "global-stoich"))
    return 3;

  return -1;
}
/* ------------------------------------------------------------------------- */

int UDFsetChiForm(char *form)
{
  if (!strncmp(form, "log", 3))
    return 1;
  else if (!strcmp(form, "erfc") || !strcmp(form, "error-function"))
    return 2;
  else if (!strncmp(form, "quad", 4))
    return 3;

  return -1;
}
/* ------------------------------------------------------------------------- */

int UDFsetChiScaling(char *scaling)
{
  if (!strcmp(scaling, "none"))
    return 0;
  else if (!strcmp(scaling, "global"))
    return 1;
  else if (!strcmp(scaling, "local"))
    return 2;
  else if (!strcmp(scaling, "global-zero"))
    return 3;
  else if (!strcmp(scaling, "local-zero"))
    return 4;

  return -1;
}
/* ------------------------------------------------------------------------- */

/* ##### start functions that only work for serial/compute node ##### */
#if !RP_HOST
void UDFcalcChi1D(Domain *d, int n1, double *x1, double *chi)
{
  Message0("Computing 1D %s scalar dissipation with %s scaling\n",
           udf.chiForm1D, udf.chiNorm1D);
  Thread *t;
  cell_t c;

  /* local variables for flow solver */
  int    j, iT;
  double taui;              /* turbulent eddy turnover time */
  double c_chi = 2.0;       /* const. for modelling mean scalar dissipation */
  double chiMean;           /* conditional mean scalar dissipation */
  double TINY  = 1.0e-6;    /* for use in potential divide-by-zero */
  double f1, tmp, den, dx;
  double chiLocal, *chiDen; /* storage for denominator f(x1) */
  double pdf1_st, *pdf;     /* storage for 1D and 2D PDFs */
  int    pdfBound[2];       /* non-zero bounds of pdf */
  double mass, totMass, *rwork;
  double Zmax, Zmin, Zplus, Zminus, DZ, xi, xiRef, sigma;
  int    jZmin, jZmax;
  double phiMean, phiVar;

  /* allocate all local memory */
  pdf    = (double*) calloc(n1, sizeof(double));
  chiDen = (double*) calloc(n1, sizeof(double));
  rwork  = (double*) calloc(n1, sizeof(double));

  /* set functional form */
  if (udf.chiForm[1] == 1)
    UDFcalcChiForm = &MDUCcalcChiLogForm;
  else if (udf.chiForm[1] == 2)
    UDFcalcChiForm = &MDUCcalcChiErfcForm;
  else if (udf.chiForm[1] == 3)
    UDFcalcChiForm = &MDUCcalcChiQuadForm;

  int Z1_d1 = 0;
  double beta;

  if (udf.iPhi[1] == udf.iH) {
    /* if enthalpy, use global scaling by default for now */
    Zmin = 0.0;
    Zmax = 1.0;
  }
  else if (udf.chiNorm[1] == 1 || udf.chiNorm[1] == 3) {
    /* get global min/max mixture fraction in the domain */      
    Zmin =  1.0e20;
    Zmax = -1.0e20;
    for (iT=0; iT<udf.nThreadSolve; iT++) {
      t = udf.solve_threads[iT]; /* current thread */
      if (t == NULL)
        continue;
      
      begin_c_loop_int(c,t) {
        /* normalize the co-ordinate */
        phiMean = (C_PHI1(c,t)  - udf.phiMin[1])/udf.deltaPhi[1];
        phiVar  = C_PHIV1(c,t)/udf.deltaPhi2[1];
        if ( phiMean > 0.0 ) {
          /* see if beta would have non-zero value at Z=1 */
          beta=phiMean*(1.0-phiMean)/(phiVar+1e-6)-1;
          beta *= (1.0-phiMean);
          if (beta < 1.0) Z1_d1 = 1;
          
          /* mixture fraction r.m.s. */
          sigma  = pow(MAX_D(phiVar,0.0), 0.5);
          Zplus  = phiMean + 2*sigma;
          Zminus = phiMean - 2*sigma;
          
          if (Zminus < Zmin) Zmin = MAX_D(Zminus, 0.0);
          if (Zplus  > Zmax) Zmax = MIN_D(Zplus,  1.0);
        }
      } end_c_loop_int(c,t);
    } /* end loop over cell threads */

#if RP_NODE
    /* get the global min/max of all nodes */
    Zmin = PRF_GRLOW1(Zmin);
    Zmax = PRF_GRHIGH1(Zmax);
#endif /* RP_NODE */
    
    Zmin = MAX_D(Zmin, 0.0);
    Zmax = MIN_D(Zmax, 1.0);
    if (udf.chiNorm[1] == 3) /* if set, use zero as minimum everywhere */
      Zmin = 0.0;
    
    /* find the mesh location to the left of Zmin */
    jZmin = UDFgridLocate(n1, x1, Zmin);
    
    /* find the mesh location to the right of Zmax */
    jZmax = UDFgridLocate(n1, x1, Zmax) + 1;
    
    /* set the min/max to the values corresponding to the grid */
    Zmin = x1[jZmin];
    Zmax = x1[jZmax];
    DZ   = Zmax - Zmin;
  } /* end get global min/max */

  /* ensure that the scalar dissipation is set to zero */
  for (j=0; j<n1; j++)
    chi[j] = chiDen[j] = 0.0;

  Message0("  max saturation mixture: %8.6f\n", Zmax);
  Message0("  Z1_d1: %i\n", Z1_d1);

  /* loop over all cells in the domain */
  totMass = 0.0;
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL)
      continue;
    
    begin_c_loop_int(c,t) {
      /* normalize the co-ordinate */
      phiMean = (C_PHI1(c,t)  - udf.phiMin[1])/udf.deltaPhi[1];
      phiVar  = C_PHIV1(c,t)/(udf.deltaPhi2[1]);
      C_CHI1_REF(c,t) = 0.0;
      
      if ( phiMean > udf.phiZero[1] ) {
        mass = C_R(c,t)*C_VOLUME(c,t);
        taui = C_D(c,t)/(C_K(c,t)+TINY); /* turb. eddy turnover time */
        
        /* compute the conditional mean scalar dissipation rate */
        chiMean = c_chi*taui*phiVar;
        
        /* compute the 1D PDF */
        /* NOTE: this function returns the PDF already scaled with dx
           i.e., which integrating equivalent to pdf(x)*dx */
        MDUCcalcPDF_1D(udf.pdf1D, n1, x1, phiMean, phiVar, pdf, pdfBound);

        if (udf.chiNorm[1] == 2 || udf.chiNorm[1] == 4) {
          /* set the bounds to consider for the local pdf
             this is taken for now as 2*sigma in each direction, as the
             cutoff for zero given by the PDF calculation is arbitrarily
             based on the cutoff threshold */
          sigma = pow(MAX_D(phiVar,0.0), 0.5);
          
          /* set local cell min/max mixture fraction */
          Zmin = 0.0;
          if (udf.chiNorm[1] == 2)
            Zmin = MAX_D(phiMean-2*sigma, 0.0);
          Zmax = MIN_D(phiMean+2*sigma, 1.0);
          
          /* find the mesh location to the left of Zmin */
          jZmin = UDFgridLocate(n1, x1, Zmin);
          
          /* find the mesh location to the right of Zmax */
          jZmax = UDFgridLocate(n1, x1, Zmax) + 1;
          
          /* set the min/max to the values corresponding to the grid */
          Zmin = x1[jZmin];
          Zmax = x1[jZmax];
          DZ   = Zmax - Zmin;
        }
        
        /* set bounds based on min/max locations */
        pdfBound[0] = jZmin;
        pdfBound[1] = jZmax;
        
        /* use the mean as the reference for the log functional form */
        xiRef = (phiMean-Zmin)/DZ;
        
        /* compute f1(x) based on assumed functional form */
        den = 0.0;
        for (j=pdfBound[0]; j<=pdfBound[1]; j++) {
          /* compute the functional form for F1(xi) as in chiMean = chi_st*F(xi)
             uses log distribution of Pitsch */
          xi = (x1[j] - Zmin)/DZ;
          /* f1 = UDFcalcChiForm(xi, xiRef, 1.0); */
          f1 = UDFcalcChiForm(x1[j], phiMean, Zmax);
          /* compute contribution to integral */
          den += f1*pdf[j];
        } /* end loop over columns */

        /* compute the conditional mean chi at stoichiometric */
        if (den > TINY)
          C_CHI1_REF(c,t) = chiMean/den;
        
        for (j=pdfBound[0]; j <= pdfBound[1]; j++) {
          xi = (x1[j]-Zmin)/(Zmax-Zmin);
          /* chi in normalised coordinates */
          /* chiLocal = C_CHI1_REF(c,t)*UDFcalcChiForm(xi, xiRef, 1.0); */
          chiLocal = C_CHI1_REF(c,t)*UDFcalcChiForm(x1[j],phiMean,Zmax);
          
          /* add contributions to global chi */
          chi[j]    += mass*chiLocal;
          chiDen[j] += mass;
          
          /* temporary error checking to see if chi[j] is defined */
          if (isnan(chi[j])) {
            Message("Zmean = %15.6e, Zvar = %15.6e\n", phiMean, phiVar);
            Message("chi1 = %15.6e\n", C_CHI1_REF(c,t));
            Message("j = %i, Z = %8.6f, xi = %15.6e\n", j, x1[j], xi);
            Message("den = %15.6e, chiLocal = %15.6e\n", den, chiLocal);
            
            Message("\nZ pdf\n");
            for (j=0; j<n1; j++)
              Message("%8.6f %15.6e\n", x1[j], pdf[j]);
            
            Error("Scalar dissipation has become undefined!\n");
          }
        }
      } /* end operations on cells with mixture fraction */
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */

#if RP_NODE
  /* if parallel, get the contribution to chi from each node */
  PRF_GRSUM(chi, n1, rwork);
  PRF_GRSUM(chiDen, n1, rwork);
#endif /* RP_NODE */

  /* compute chi(Z) over the whole grid */
  for (j=0; j<n1; j++) {
    if (chiDen[j] > 0.0)
      chi[j] = chi[j]/chiDen[j];
  } /* end loop over flamelet grid */

  free(pdf); free(chiDen); free(rwork);
}
/* ------------------------------------------------------------------------- */

void UDFbinChi1D(Domain *d, int n1, double *x1, double *chi)
{
  Message0("Binning 1D scalar dissipation rate...\n");
  Thread *t;
  cell_t c;

  /* local variables for flow solver */
  int    i, j, k, iT;
  double TINY = 1.0e-9; /* for use in potential divide-by-zero */
  double tmp;
  double chiMean;
  double *chiDen;        /* temporary storage for integral of denominator */
  double *rwork;
  double taui;           /* turbulent eddy turnover time */
  double c_chi = 2.0;    /* const. for modelling mean scalar dissipation */

  double dx, mass;
  double w1, w2;         /* weights for distributing chi */

  double phiMean;        /* mean of scalar */
  double phiVar;         /* variance of scalar */

  /* allocate all local memory */
  chiDen = (double*) calloc(n1, sizeof(double));
  rwork  = (double*) calloc(n1, sizeof(double));

  for (j=0; j<n1; j++)
    chi[j] = chiDen[j] = 0.0;

  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL)
      continue;
    
    begin_c_loop_int(c,t) {
      /* normalize the co-ordinate */
      phiMean = (C_PHI1(c,t) - udf.phiMin[1])/udf.deltaPhi[1];
      phiVar  = C_PHIV1(c,t)/udf.deltaPhi2[1];
      
      if ( phiMean > TINY ) {
        mass = C_R(c,t)*C_VOLUME(c,t);
        taui = C_D(c,t)/(C_K(c,t)+TINY); /* turb. eddy turnover time */
        
        /* compute the unconditional mean scalar dissipation rate */
        C_CHI1_REF(c,t) = chiMean = c_chi*taui*phiVar;
        
        /* search for the bin that it fits in */
        j = UDFgridLocate(n1, x1, phiMean);
        dx = x1[j+1] - x1[j];
        
        /* take care of the boundary points */
        if (j == 0) {
          w1 = 0.0;
          w2 = 1.0;
        }
        else if (j == n1-2) {
          w1 = 1.0;
          w2 = 0.0;
        }
        else {
          w1 = (x1[j+1] - phiMean)/dx;
          w2 = (phiMean - x1[j])/dx;
        }
        
        chi[j]    += mass*w1*chiMean;
        chiDen[j] += mass*w1;
        
        chi[j+1]    += mass*w2*chiMean;
        chiDen[j+1] += mass*w2;
      }
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */

#if RP_NODE
  /* if parallel, get the contribution to chi from each node */
  PRF_GRSUM(chi, n1, rwork);
  PRF_GRSUM(chiDen, n1, rwork);
#endif /* RP_NODE */
  
  /* compute chi(Z) over the whole grid */
  for (j=0; j<n1; j++) {
    if (chiDen[j] > 0.0)
      chi[j] = chi[j]/chiDen[j];
  } /* end loop over flamelet grid */

  free(chiDen);
  free(rwork);
  Message0("1D binning is done\n");
}
/* ------------------------------------------------------------------------- */

void UDFcalcChiPitsch(Domain *d, int n1, double *x1, double Z1ref, double *chi)
{
  Message0("Computing scalar dissipation rate according to Pitsch...\n");
  Thread *t;
  cell_t c;

  /* local variables for flow solver */
  int    i, j, k, iT;
  int    jZmax, jRef;
  double TINY  = 1.0e-9; /* for use in potential divide-by-zero */
  double tmp, wRef;
  double chiDen;         /* temporary storage for integral of denominator */
  double pdf1_st, *pdf;  /* storage for PDFs */
  int    pdfBound[2];    /* stores bounds of non-zero PDF values */
  double taui;           /* turbulent eddy turnover time */
  double c_chi = 2.0;    /* const. for modelling mean scalar dissipation */
  double chiMean;        /* unconditional mean scalar dissipation */
  double chiStAvg;       /* domain-averaged conditional mean chi */
  double sigma;          /* standard deviation */
  double Zmax = 1.0;     /* maximum mixture fraction in the domain */

  /* allocate all local memory */
  pdf    = (double*) calloc(n1, sizeof(double));

  if (udf.chiNorm[1] == 1) {
    /* ----- get the max in the domain ----- */
    /* this is defined as the mean + 2*sigma */
    Zmax = 0.0;
    for (iT=0; iT<udf.nThreadSolve; iT++) {
      t = udf.solve_threads[iT]; /* current thread */
      if (t == NULL) continue;

      begin_c_loop_int(c,t) {
        if ( (C_Z1MEAN(c,t) > 0.0) && (C_Z1VAR(c,t) > TINY) ) {
          sigma = pow(MAX_D(C_Z1VAR(c,t),0), 0.5);
          
          if (C_Z1MEAN(c,t) + 2*sigma > Zmax)
            Zmax = C_Z1MEAN(c,t) + 2*sigma;
        }
      } end_c_loop_int(c,t);
    } /* end loop over cell threads */

      /* find the mesh location to the right of Zmax */
    j = 0;
    while (x1[j] < Zmax) j++;
    if (j >= n1-1)
      jZmax = n1-1;
    else
      jZmax = j;
  } /* end set global Zmax */

  /* find the value to the left of stoichiometric for interoplation later */
  jRef = 0;
  do { jRef++; } while ( x1[jRef] < Z1ref && jRef < n1 );
  jRef--; /* assume for within 1 grid point */
  
  wRef = (Z1ref-x1[jRef])/(x1[jRef+1]-x1[jRef]);
  
  /* ----- compute mean scalar dissipation for each cell ----- */
  Message("Computing mean scalar dissipation for cells...\n");
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;

    begin_c_loop_int(c,t) {
      C_CHI1_REF(c,t) = 0.0;
      
      if ( C_Z1MEAN(c,t) > udf.phiZero[1] ) {
        taui = C_D(c,t)/C_K(c,t); /* turb. eddy turnover time */

        /* compute the unconditional mean scalar dissipation rate */
        chiMean = c_chi*taui*C_Z1VAR(c,t);
        
        /* compute the PDF based on beta function*/
        /* NOTE: this function returns the PDF already scaled with dx
           i.e., which integrating equivalent to pdf(x)*dx */
        MDUCcalcPDF_1D(udf.pdf1D,n1,x1,C_Z1MEAN(c,t),C_Z1VAR(c,t),pdf,pdfBound);

        if (udf.chiNorm[1] == 2) {
          /* set the bounds to consider for the local pdf
             this is taken for now as 2*sigma in each direction, as the
             cutoff for zero given by the PDF calculation is arbitrarily
             based on the cutoff threshold */
          sigma = pow(MAX_D(C_Z1VAR(c,t),0.0), 0.5);
          
          /* set local cell min/max mixture fraction */
          Zmax = MIN_D(C_Z1MEAN(c,t)+2*sigma, 1.0);
          
          /* find the mesh location to the right of Zmax */
          j = 0;
          while (x1[j] < Zmax) j++;
          if (j >= n1-1)
            jZmax = n1-1;
          else
            jZmax = j;
          
          /* set the min/max to the values corresponding to the grid */
          Zmax = x1[jZmax];
        } /* end get local cell Zmax*/
        
        /* integrate the denominator using the logarithmic functional
           form of Pitsch */
        chiDen = 0.0;
        for (j=pdfBound[0]; j<=pdfBound[1]; j++)
          chiDen += MDUCcalcChiLogForm(x1[j], Z1ref, Zmax)*pdf[j];
        
        /* compute the conditional mean chi at stoichiometric */
        if (chiDen > TINY)
          C_CHI1_REF(c,t) = chiMean/chiDen;
        /* else */
        /* Message("Warning: chiDen = %15.6e in cell %i\n",chiDen, c); */
      }
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */

  /* ----- compute the reference scalar dissipation rate ----- */
  /* domain average the conditional mean scalar dissipation rate weighted
     by the stoichiometric surface area */
  Message("Computing reference scalar dissipation rate...\n");
  chiStAvg  = 0.0;
  chiDen    = 0.0;
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;
    
    begin_c_loop_int(c,t) {
      /* compute the PDF(Z1_st, 0) */
      MDUCcalcPDF_1D(udf.pdf1D,n1,x1,C_Z1MEAN(c,t),C_Z1VAR(c,t),pdf,pdfBound);

      /* interpolate to find the value at Z1Ref */
      pdf1_st = pdf[jRef] + (pdf[jRef+1]-pdf[jRef])*wRef;
      
      /* add to integral for chi1 */
      if (C_CHI1_REF(c,t) > 0.0) {
        /* sum the contributions to numerator and denominator */
        tmp = pow(C_CHI1_REF(c,t),0.5)*C_R(c,t)*C_VOLUME(c,t)*pdf1_st;
        chiStAvg += C_CHI1_REF(c,t)*tmp;
        chiDen   += tmp;
        
        /* Error checking on average chi */
        if (isnan(chiStAvg)) {
          Message("thread %i, cell %i\n",t,c);
          Message("C_Z1MEAN = %15.6e, C_Z1VAR = %15.6e\n",
                  C_Z1MEAN(c,t),C_Z1VAR(c,t));
          Message("chi1_st = %15.6e, pdf1_st = %15.6e, tmp = %15.6e\n",
                  C_CHI1_REF(c,t), pdf1_st, tmp);
          Error("<chiStoich> is undefined\n");
        } /* end error checking */
      } /* end chi computation for non-zero cell */
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */

#if RP_NODE
  /* if parallel, get the contribution to chi from each node */
  PRF_GRSUM1(chiStAvg);
  PRF_GRSUM1(chiDen);
#endif /* RP_NODE */
  
  chiStAvg = chiStAvg/chiDen; /* complete the averaging */
  
  Message("<chi1_stoich> = %15.6e\n", chiStAvg);

  /* ---- compute the scalar dissipation rate profiles ----- */
  /* for now, compute 1D profiles for each direction based on reference value
     at stoichiometric */
  for (j=0; j<n1; j++)
    chi[j] = chiStAvg*MDUCcalcChiLogForm(x1[j], Z1ref, Zmax);

  free(pdf);
}
/* ------------------------------------------------------------------------- */

void UDFbinChi2D(Domain *d, int n1, int n2, double *x1, double *x2,
                 double **chi1, double **chi2)
{
  Message0("Binning 2D scalar dissipation rate...\n");
  Thread *t;
  cell_t c;

  int    i, j, k, iT;
  double w00, w01, w10, w11;             /* corner weighting factors */
  double A00, A01, A10, A11, Atot, Amod; /* areas for weighting  */
  double TINY = 1.0e-12;
  double mass;
  double tau, c_chi = 2.0;
  double x, y, dx, dy;
  double phiMean[3], phiVar[3];

  double chi1Mean, chi2Mean;
  double **chi1Den, **chi2Den;  /* storage for denominators */
  double *rwork;

  /* allocate local memory */
  chi1Den = (double**) calloc(n2, sizeof(double*));
  chi2Den = (double**) calloc(n2, sizeof(double*));
  rwork   = (double*)  calloc(n1, sizeof(double));

  for (k=0; k < n2; k++) {
    chi1Den[k] = (double*) calloc(n1, sizeof(double));
    chi2Den[k] = (double*) calloc(n1, sizeof(double));
  }

  /* make sure everything is set to zero*/
  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      chi1[k][j] = chi2[k][j] = chi1Den[k][j] = chi2Den[k][j] = 0.0;

   for (iT=0; iT<udf.nThreadSolve; iT++) {
     t = udf.solve_threads[iT]; /* current thread */
     if (t == NULL)
       continue;
     
     begin_c_loop_int(c,t) {
       /* normalize the co-ordinate */
       phiMean[1] = (C_PHI1(c,t)  - udf.phiMin[1])/udf.deltaPhi[1];
       phiVar[1]  = C_PHIV1(c,t)/udf.deltaPhi2[1];
       phiMean[2] = (C_PHI2(c,t)  - udf.phiMin[2])/udf.deltaPhi[2];
       phiVar[2]  = C_PHIV2(c,t)/udf.deltaPhi2[2];
       
       /* only consider a PDF over the cells that have a non-zero mixture */
       if ( (C_PHI1(c,t) > udf.phiZero[1]) || (C_PHI2(c,t) > udf.phiZero[2]) ) {
         mass = C_R(c,t)*C_VOLUME(c,t);     /* cell mass */
         tau = C_D(c,t)/(C_K(c,t)+TINY);    /* turb. eddy turnover time */
         chi1Mean = c_chi*tau*C_PHIV1(c,t);
         chi2Mean = c_chi*tau*C_PHIV2(c,t);
         
         /* select transformation of grid */
         if (udf.transform == 0) { /* no transformation*/
           x = phiMean[1];
           y = phiMean[2];
         }
         else if (udf.transform == 1) {
           x = C_PHI1(c,t);
           y = C_PHI2(c,t)/(1-C_PHI1(c,t));
         }
         else if (udf.transform == 2) { /* Hasse's transformation */
           x = C_PHI1(c,t)+C_PHI2(c,t);
           y = C_PHI2(c,t)/(C_PHI1(c,t)+C_PHI2(c,t));
         }
         
         /* search for the bin that it fits in */
         j = UDFgridLocate(n1, x1, phiMean[1]);
         k = UDFgridLocate(n2, x2, phiMean[2]);
         
         /* add weighted distribution to each of the four points
            surrounding the value */
         /* compute the relative areas */
         Atot = (x1[j+1] - x1[j])*(x2[k+1] - x2[k]); /* dx*dy */
         A00  = (x1[j+1] - x)*(x2[k+1] - y);
         A01  = (x1[j+1] - x)*(y - x2[k]);
         A10  = (x - x1[j])*(x2[k+1] - y);
         A11  = (x - x1[j])*(y - x2[k]);
         
         /* compute weights based on areas */
         w00 = A00/Atot;
         w01 = A01/Atot;
         w10 = A10/Atot;
         w11 = A11/Atot;
         
         /* account for the corner points where know chi = 0 */
         if ( (j == 0) && (k == 0) ) {           /* bottom left corner */
           Amod = Atot - A00;
           w01 += A01/Amod*w00;
           w10 += A10/Amod*w00;
           w11 += A11/Amod*w00;
           w00  = 0.0;
         }
         else if ( (j == 0) && (k == n2-2) ) {   /* top left corner */
           Amod = Atot - A01;
           w00 += A00/Amod*w01;
           w10 += A10/Amod*w01;
           w11 += A11/Amod*w01;
           w01  = 0.0;
         }
         else if ( (j == n1-2) && (k == 0) ) {   /* bottom right corner */
           Amod = Atot - A10;
           w00  = A00/Amod*w10;
           w01  = A01/Amod*w10;
           w11  = A11/Amod*w10;
           w10  = 0.0;
         }
         else if ( (j == n1-2) && (k == n2-2) ) { /* top right corner */
           Amod = Atot - A11;
           w00  = A00/Amod*w11;
           w01  = A01/Amod*w11;
           w10  = A10/Amod*w11;
           w11  = 0.0;
         }
         
         /* ----- add contributions to each corner ----- */
         if (C_PHI1(c,t) > udf.phiZero[1]) {
           if (w00 > TINY) { /* bottom left corner */
             chi1[k][j]    += mass*w00*chi1Mean;
             chi1Den[k][j] += mass*w00;
           }
           if (w01 > TINY) { /* top left corner */
             chi1[k+1][j]    += mass*w01*chi1Mean;
             chi1Den[k+1][j] += mass*w01;
           }
           if (w10 > TINY) { /* bottom right corner */
             chi1[k][j+1]    += mass*w10*chi1Mean;
             chi1Den[k][j+1] += mass*w10;
           }
           if (w11 > TINY) { /* top right corner */
             chi1[k+1][j+1]    += mass*w11*chi1Mean;
             chi1Den[k+1][j+1] += mass*w11;
           }
         }

         /* ----- add contributions to each corner ----- */
         if (C_PHI2(c,t) > udf.phiZero[2]) {
           if (w00 > TINY) { /* bottom left corner */
             chi2[k][j]    += mass*w00*chi2Mean;
             chi2Den[k][j] += mass*w00;
           }
           if (w01 > TINY) { /* top left corner */
             chi2[k+1][j]    += mass*w01*chi2Mean;
             chi2Den[k+1][j] += mass*w01;
           }
           if (w10 > TINY) { /* bottom right corner */
             chi2[k][j+1]    += mass*w10*chi2Mean;
             chi2Den[k][j+1] += mass*w10;
           }
           if (w11 > TINY) { /* top right corner */
             chi2[k+1][j+1]    += mass*w11*chi2Mean;
             chi2Den[k+1][j+1] += mass*w11;
           }
         }
         
         /* error checking */
         if ( isnan(chi1[k+1][j]) ) {
           Message("Z1mean = %15.6e, Z1var = %15.6e\n", C_Z1MEAN(c,t),
                   C_Z1VAR(c,t));
           Message("Z2mean = %15.6e, Z2var = %15.6e\n", C_Z2MEAN(c,t),
                   C_Z2VAR(c,t));
           Message("chi1[%i][%i]    = %15.6e\n", k, j, chi1[k][j]);
           Message("chi2[%i][%i]    = %15.6e\n", k, j, chi2[k][j]);
           Message("chi1Den[%i][%i] = %15.6e\n", k, j, chi1Den[k][j]);
           Message("chi2Den[%i][%i] = %15.6e\n", k, j, chi2Den[k][j]);
           
           Error("scalar dissipation undefined!\n");
         }
         
       } /* end chi computation */
     } end_c_loop_int(c,t);
   } /* end loop over cell threads */
   
#if RP_NODE
  /* if parallel, get the contribution to chi from each node */
  for (k=0; k<n2; k++) {
    PRF_GRSUM(chi1[k], n1, rwork);
    PRF_GRSUM(chi2[k], n1, rwork);
    PRF_GRSUM(chi1Den[k], n1, rwork);
    PRF_GRSUM(chi2Den[k], n1, rwork);
  }
#endif /* RP_NODE */

  for (k=0; k<n2; k++) {
    for (j=0; j<n1; j++) {
      /* finish normalisation for chi1 */
      if (chi1Den[k][j] > 0.0)
        chi1[k][j] = chi1[k][j]/chi1Den[k][j];
      else
        chi1[k][j] = 0.0;
      
      /* finish normalisation for chi2 */
      if (chi2Den[k][j] > 0.0)
        chi2[k][j] = chi2[k][j]/chi2Den[k][j];
      else
        chi2[k][j] = 0.0;
    } /* end loop over columns */
  } /* end loop over rows */
  
  for (k=0; k < n2; k++) {
    free(chi1Den[k]);
    free(chi2Den[k]);
  }
  free(chi1Den);
  free(chi2Den);
  free(rwork);
}
/* ------------------------------------------------------------------------- */

void UDFcalcChi2D(Domain *d, int *nMesh, double **xMesh, double **chi1,
                  double **chi12, double **chi2)
{
  Message0("\n----- Computing 2D scalar dissipation rate ");
  Thread *t;
  cell_t c;

  int    i, j, k, iT;
  int    nZ, nZ1, nZ2, counter, nCount = 500, lineCount = 200 ;
  int    nNotCvg, nZ1L, nZ2L;
  int    largeChi[3];
  double largeChiAvg[3];
  double TINY = 1.0e-6;
  double mass, tau, c_chi = 2.0;

  double chiMean[3], chiLocal[3], den[3];
  double **chi0, **chi0Den, **chi1Den, **chi2Den;
  double sigma1, sigma2;

  double **Z1, **Z2, *Z1min, *Z1max, *Z2min, *Z2max, ZmaxG[3];
  double Zplus, Zminus;
  double **f0, **f1, **f2;
  double tmp, *rwork;
  char   pdf2Mode[50];

  /* storage for the pdf */
  int    pdfFlag, pdfBound1D[2], **pdfBound;
  double **pdf, *pdfX1, *pdfX2;

  /* dereference the mesh for readable notation */
  int     n1   = nMesh[0];
  int     n2   = nMesh[1];
  double *x1   = xMesh[0];
  double *x2   = xMesh[1];
  double dx1_0 = 0.5*(x1[1]-x1[0]);
  double dx2_0 = 0.5*(x2[1]-x2[0]);

  /* set chi to zero everywhere */
  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      chi1[k][j] = chi2[k][j] = chi2[k][j] = 0.0;

  /* if want zero scalar dissipation, set flag in input as negative */
  if (udf.chiScale < 0.0) {
    Message0("Using zero scalar dissipation rate\n");
    return;
  }

  /* set functional form */
  if (udf.chiForm[2] == 1)
    UDFcalcChiForm = &MDUCcalcChiLogForm;
  else if (udf.chiForm[2] == 2)
    UDFcalcChiForm = &MDUCcalcChiErfcForm;
  else if (udf.chiForm[2] == 3)
    UDFcalcChiForm = &MDUCcalcChiQuadForm;

  /* ----- allocate local memory ----- */
  Z1       = (double**) calloc(n2, sizeof(double*));
  Z2       = (double**) calloc(n2, sizeof(double*));
  Z1min    = (double*)  calloc(n2, sizeof(double));
  Z1max    = (double*)  calloc(n2, sizeof(double));
  Z2min    = (double*)  calloc(n1, sizeof(double));
  Z2max    = (double*)  calloc(n1, sizeof(double));
  f0       = (double**) calloc(n2, sizeof(double*));
  f1       = (double**) calloc(n2, sizeof(double*));
  f2       = (double**) calloc(n2, sizeof(double*));
  pdf      = (double**) calloc(n2, sizeof(double*));
  pdfX1    = (double*)  calloc(n1, sizeof(double));
  pdfX2    = (double*)  calloc(n2, sizeof(double));
  pdfBound = (int**)    calloc(n2, sizeof(int*));
  chi0     = (double**) calloc(n2, sizeof(double*));
  chi0Den  = (double**) calloc(n2, sizeof(double*));
  chi1Den  = (double**) calloc(n2, sizeof(double*));
  chi2Den  = (double**) calloc(n2, sizeof(double*));
  rwork    = (double*)  calloc(MAX_I(n1,n2), sizeof(double));

  for (k=0; k<n2; k++) {
    Z1[k]       = (double*) calloc(n1, sizeof(double));
    Z2[k]       = (double*) calloc(n1, sizeof(double));
    f0[k]       = (double*) calloc(n1, sizeof(double));
    f1[k]       = (double*) calloc(n1, sizeof(double));
    f2[k]       = (double*) calloc(n1, sizeof(double));
    pdf[k]      = (double*) calloc(n1, sizeof(double));
    pdfBound[k] = (int*)    calloc(2, sizeof(int));
    chi0[k]     = (double*) calloc(n1, sizeof(double));
    chi0Den[k]  = (double*) calloc(n1, sizeof(double));
    chi1Den[k]  = (double*) calloc(n1, sizeof(double));
    chi2Den[k]  = (double*) calloc(n1, sizeof(double));
  } /* end loop over rows */

  /* set chi, etc. to zero everywhere */
  for (k=0; k<n2; k++) {
    for (j=0; j<n1; j++) {
      pdf[k][j]  = 0.0;
      chi0[k][j] = 0.0;
      chi0Den[k][j] = chi1Den[k][j] = chi2Den[k][j] = 0.0;
    } /* end loop over columns */
  } /* end loop over rows */

  /* fill the non-transformed value grid using the appropriate transformation */
  for (k=0; k<n2; k++) {
    for (j=0; j<n1; j++) {
      /* apply transformation */
      if (udf.transform == 0) {
        Z1[k][j] = x1[j];
        Z2[k][j] = x2[k];
      }
      else if (udf.transform == 1) {
        Z1[k][j] = x1[j];
        Z2[k][j] = x2[k]*(1.0-x1[j]);
      }
      else if (udf.transform == 2) {
        Z1[k][j] = x1[j]*(1.0 - x2[k]);
        Z2[k][j] = x2[k]*x1[j];
      }
    } /* end loop over columns */
  } /* end loop over rows */
  
  int    Z1_d1 = 0, Z2_d1 = 0;
  double beta;

  /* compute the global min/max of each mixture fraction in the domain */
  for (k=0; k<n2; k++) { Z1min[k] = 1.0; Z1max[k] = 0.0; }
  for (j=0; j<n1; j++) { Z2min[j] = 1.0; Z2max[j] = 0.0; }
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL)
      continue;
    
    begin_c_loop_int(c,t) {
      if (C_PHI1(c,t) > 0.0) {
        /* check whether beta distribution would be u-shaped/increasing */
        beta  = C_PHI1(c,t)*(1.0-C_PHI1(c,t))/(C_PHIV1(c,t)+1e-6) - 1;
        beta *= (1.0-C_PHI1(c,t));
        if (beta < 1.0) Z1_d1 = 1;
        
        /* get r.m.s. value */
        sigma1 = pow(MAX_D(C_PHIV1(c,t),0), 0.5);
        
        /* find where to condition */
        k = UDFgridLocate(n2, x2, C_PHI2(c,t));
        if (C_PHI2(c,t) > 0.5*(x2[k]+x2[k+1]))
          k = k+1;
        
        /* assume min/max is +/- 2*sigma */
        Zplus  = C_PHI1(c,t) + 2*sigma1;
        Zminus = C_PHI1(c,t) - 2*sigma1;
        if (Zminus < Z1min[k]) Z1min[k] = MAX_D(Zminus, 0.0);
        if (Zplus  > Z1max[k]) Z1max[k] = MIN_D(Zplus,  1.0);
      }
      
      if (C_PHI2(c,t) > 0.0) {
        /* check whether beta distribution would be u-shaped/increasing */
        beta  = C_PHI2(c,t)*(1.0-C_PHI2(c,t))/(C_PHIV2(c,t)+1e-6) - 1;
        beta *= (1.0-C_PHI2(c,t));
        if (beta < 1.0) Z2_d1 = 1;
        
        /* get r.m.s. value */
        sigma2 = pow(MAX_D(C_PHIV2(c,t),0), 0.5);
        
        /* find where to condition */
        j = UDFgridLocate(n1, x1, C_PHI1(c,t));
        if (C_PHI1(c,t) > 0.5*(x1[j]+x1[j+1]))
          j = j+1;
        
        /* assume min/max is +/- 2*sigma */
        Zplus  = C_PHI2(c,t) + 2*sigma2;
        Zminus = C_PHI2(c,t) - 2*sigma2;
        if (Zminus < Z2min[j]) Z2min[j] = MAX_D(Zminus, 0.0);
        if (Zplus  > Z2max[j]) Z2max[j] = MIN_D(Zplus,  1.0);
      }
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */

#if RP_NODE
  /* accumulate min/max values over all nodes*/
  PRF_GRLOW(Z1min, n2, rwork);
  PRF_GRHIGH(Z1max, n2, rwork);
  PRF_GRLOW(Z2min, n1, rwork);
  PRF_GRHIGH(Z2max, n1, rwork);
#endif /* RP_NODE */

  /* get the global max values for the domain */
  for (j=0, ZmaxG[2]=0.0; j<n1; j++) {
    if (Z2max[j] > ZmaxG[2])
      ZmaxG[2] = Z2max[j];
  }
  for (k=0, ZmaxG[1]=0.0; k<n2; k++) {
    if (Z1max[k] > ZmaxG[1])
      ZmaxG[1] = Z1max[k];
  }
  ZmaxG[0] = 1.0-ZmaxG[1]-ZmaxG[2];

  Message0(" Max Mixture fractions:\n");
  Message0("   Z0: %8.6f\n", ZmaxG[0]);
  Message0("   Z1: %8.6f\n", ZmaxG[1]);
  Message0("   Z2: %8.6f\n", ZmaxG[2]);
  Message0("   Z1_d1 --> %i\n", Z1_d1);
  Message0("   Z2_d1 --> %i\n", Z2_d1);

  Message0(" Starting computation of scalar dissipation rate...\n");
  counter = nZ = nZ1 = nZ2 = 0;
  nNotCvg = nZ1L = nZ2L = 0;
  for (i=0; i<3; i++) largeChi[i]    = 0;
  for (i=0; i<3; i++) largeChiAvg[i] = 0.0;
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL)
      continue;
    
    begin_c_loop_int(c,t) {
      /* ensure all storage memory is reset to zero */
      C_CHI0_REF(c,t) = C_CHI1_REF(c,t) = C_CHI2_REF(c,t) = 0.0;
      for (k=0; k<n2; k++)
        for (j=0; j<n1; j++)
          f0[k][j] = f1[k][j] = f2[k][j] = 0.0;
      
      if ( (C_PHI1(c,t) > udf.phiZero[1]) ) {
        mass = C_R(c,t)*C_VOLUME(c,t);   /* cell mass */
        tau  = C_D(c,t)/C_K(c,t); /* turb. eddy turnover time */
        
        /* check if the cell contains both mixture fractions */
        if (C_PHI2(c,t) > udf.phiZero[2]) {
          /* compute mean scalar dissipation rates */
          chiMean[0] = c_chi*tau*C_PHIV0(c,t);
          chiMean[1] = c_chi*tau*C_PHIV1(c,t);
          chiMean[2] = c_chi*tau*C_PHIV2(c,t);
          
          /* set the default joint pdf treatment */
          strcpy(pdf2Mode, udf.pdf2D);
          
          /* if ( C_PHI1(c,t) < dx1_0 || C_PHI2(c,t) < dx2_0 ) */
          strcpy(pdf2Mode, "dirichlet");
          
          /* compute 2D probability distribution */
          pdfFlag = MDUCcalcPDF_2D(pdf2Mode, nMesh, xMesh,
                                   C_PHI1(c,t), C_PHIV1(c,t),
                                   C_PHI2(c,t), C_PHIV2(c,t),
                                   C_PHIV0(c,t), pdf, pdfBound,
                                   udf.transform);
          
          if (pdfFlag == -10) {
            nNotCvg++;
          }
          else if (pdfFlag == 102) {
            nZ1L++;
          }
          else if (pdfFlag == 202) {
            nZ2L++;
          }
          
          den[0] = den[1] = den[2] = 0.0;
          /* compute the integral of the functional form to get 
             mean scalar dissipation */
          for (k=0; k<n2; k++) {
            for (j=0; j<n1; j++) {
              /* only consider region bounded by max mixture */
              if (Z2[k][j] <= ZmaxG[2]*(1.0-Z1[k][j]/ZmaxG[1])) {
                f0[k][j] = UDFcalcChiForm(Z1[k][j]/ZmaxG[1] + Z2[k][j]/ZmaxG[2],
                                          C_PHI1(c,t)/ZmaxG[1]+
                                          C_PHI2(c,t)/ZmaxG[2], 1.0);
                f1[k][j] = UDFcalcChiForm(Z1[k][j],C_PHI1(c,t),ZmaxG[1]);
                f2[k][j] = UDFcalcChiForm(Z2[k][j],C_PHI2(c,t),ZmaxG[2]);

                den[0] += f0[k][j]*pdf[k][j];
                den[1] += f1[k][j]*pdf[k][j];
                den[2] += f2[k][j]*pdf[k][j];
              }
            } /* end loop over columns */
          } /* end loop over rows */ 
          
          /* compute the reference scalar dissipation rate */
          if (den[0] > TINY) C_CHI0_REF(c,t) = chiMean[0]/den[0];
          if (den[1] > TINY) C_CHI1_REF(c,t) = chiMean[1]/den[1];
          if (den[2] > TINY) C_CHI2_REF(c,t) = chiMean[2]/den[2];
          
          if (C_CHI0_REF(c,t) > 100.0) {
            largeChi[0]++;
            largeChiAvg[0] += C_CHI0_REF(c,t);
          }
          if (C_CHI1_REF(c,t) > 100.0) {
            largeChi[1]++;
            largeChiAvg[1] += C_CHI1_REF(c,t);
          }
          if (C_CHI2_REF(c,t) > 100.0) {
            largeChi[2]++;
            largeChiAvg[2] += C_CHI2_REF(c,t);
          }
          
          C_CHI0_REF(c,t) = MIN_D(C_CHI0_REF(c,t), 100.0);
          C_CHI1_REF(c,t) = MIN_D(C_CHI1_REF(c,t), 100.0);
          C_CHI2_REF(c,t) = MIN_D(C_CHI2_REF(c,t), 100.0);
          
          /* scale the distribution and add contribution to mass avg. */
          for (k=0; k<n2; k++) {
            for (j=0; j<n1; j++) {
              chiLocal[0] = C_CHI0_REF(c,t)*f0[k][j];
              chiLocal[1] = C_CHI1_REF(c,t)*f1[k][j];
              chiLocal[2] = C_CHI2_REF(c,t)*f2[k][j];
              
              /* if (chiLocal[2] > 1000) */
              /*                             fprintf(stderr, "chi huge, Z2 mean = %8.6f\n", C_PHI2(c,t)); */
              
              chi0[k][j]    += mass*chiLocal[0];
              chi0Den[k][j] += mass;
              chi1[k][j]    += mass*chiLocal[1];
              chi1Den[k][j] += mass;
              chi2[k][j]    += mass*chiLocal[2];
              chi2Den[k][j] += mass;
            } /* end loop over columns */
          } /* end loop over rows */ 
          nZ++;
        }
        else {
          /* compute mean scalar dissipation rates */
          chiMean[1] = c_chi*tau*C_PHIV1(c,t);
          
          /* this cell only has Z1 */
          MDUCcalcPDF_1D(udf.pdf1D, n1, x1, C_PHI1(c,t),C_PHIV1(c,t),
                         pdfX1, pdfBound1D);

          /* compute f1(x) based on assumed functional form */
          den[1] = 0.0;
          for (j=0; j<n1; j++) {
            /* compute the functional form */
            f1[0][j] = UDFcalcChiForm(x1[j],C_PHI1(c,t),ZmaxG[1]);
            
            /* compute contribution to integral */
            den[1] += f1[0][j]*pdfX1[j];
          } /* end loop over columns */
          
          /* compute the reference scalar dissipation rate */
          if (den[1] > TINY) C_CHI1_REF(c,t) = chiMean[1]/den[1];
          
          /* scale the distribution and add contribution to mass avg. */
          for (j=0; j<n1; j++) {
            chiLocal[1] = C_CHI1_REF(c,t)*f1[0][j];
            
            /* add contributions to global chi */
            chi0[0][j]    += mass*chiLocal[1];
            chi0Den[0][j] += mass;
            chi1[0][j]    += mass*chiLocal[1];
            chi1Den[0][j] += mass;
          } /* end loop over columns */
          nZ1++;
        }
        
      } /* end Z1 */
      else if (C_PHI2(c,t) > udf.phiZero[2]) { /* cell only has Z2 */
        mass = C_R(c,t)*C_VOLUME(c,t);     /* cell mass */
        tau  = C_D(c,t)/C_K(c,t);   /* turb. eddy turnover time */
        chiMean[2] = c_chi*tau*C_PHIV2(c,t);
        
        /* this cell only has Z2, compute 1D distribution */
        MDUCcalcPDF_1D(udf.pdf1D, n2, x2, C_PHI2(c,t), C_PHIV2(c,t),
                       pdfX2, pdfBound1D);
        
        /* compute f2(y) based on assumed functional form */
        den[2] = 0.0;
        for (k=0; k<n2; k++) {
          /* compute functional form */
          f2[k][0] = UDFcalcChiForm(x2[k],C_PHI2(c,t),ZmaxG[2]);
          
          /* compute contribution to integral */
          den[2] += f2[k][0]*pdfX2[k];
        } /* end loop over rows */
        
        /* compute the reference scalar dissipation rate */
        if (den[2] > TINY)  C_CHI2_REF(c,t) = chiMean[2]/den[2];
        
        /* scale the distribution and add contribution to mass avg. */
        for (k=0; k<n2; k++) {
          chiLocal[2] = C_CHI2_REF(c,t)*f2[k][0];
          
          /* add contributions to global chi */
          chi0[k][0]    += mass*chiLocal[2];
          chi0Den[k][0] += mass;
          chi2[k][0]    += mass*chiLocal[2];
          chi2Den[k][0] += mass;
        } /* end loop over rows */
        nZ2++;
      }
      
      counter++;
    } end_c_loop_int(c,t); /* end loop over all cells */
  } /* end loop over cell threads */

#if RP_NODE
  /* sum up all the counters */
  nZ      = PRF_GISUM1(nZ);
  nZ1     = PRF_GISUM1(nZ1);
  nZ2     = PRF_GISUM1(nZ2);
  nZ1L    = PRF_GISUM1(nZ1L);
  nZ2L    = PRF_GISUM1(nZ2L);
  nNotCvg = PRF_GISUM1(nNotCvg);
  counter = PRF_GISUM1(counter);
  for (i=0; i<3; i++) largeChi[i] = PRF_GISUM1(largeChi[i]);
  for (i=0; i<3; i++) largeChiAvg[i] = PRF_GRSUM1(largeChiAvg[i]);

  /* sum the contribution of all dissipation rates */
  for (k=0; k<n2; k++) {
    PRF_GRSUM(chi0[k], n1, rwork);
    PRF_GRSUM(chi0Den[k], n1, rwork);
    PRF_GRSUM(chi1[k], n1, rwork);
    PRF_GRSUM(chi1Den[k], n1, rwork);
    PRF_GRSUM(chi2[k], n1, rwork);
    PRF_GRSUM(chi2Den[k], n1, rwork);
  }
#endif /* RP_NODE */
 
  /* normalise the distribution */
  for (k=0; k<n2; k++) {
    for (j=0; j<n1; j++) {
      if (chi0Den[k][j] > TINY)
        chi0[k][j] = chi0[k][j]/chi0Den[k][j];
      else
        chi0[k][j] = 0.0;
      
      if (chi1Den[k][j] > TINY)
        chi1[k][j] = chi1[k][j]/chi1Den[k][j];
      else
        chi1[k][j] = 0.0;
      
      if (chi2Den[k][j] > TINY)
        chi2[k][j] = chi2[k][j]/chi2Den[k][j];
      else
        chi2[k][j] = 0.0;
      
      chi12[k][j] = 0.5*(chi0[k][j] - chi1[k][j] - chi2[k][j]);
      
      chi1[k][j]  *= udf.chiScale;
      chi12[k][j] *= udf.chiScale;
      chi2[k][j]  *= udf.chiScale;
    } /* end loop over columns */
  } /* end loop over rows */
  
  /* enforce zero boundary conditions as appropriate */
  for (k=0; k<n2; k++) {
    if (udf.transform == 1) {
      chi1[k][0] = chi1[k][n1-1] = 0.0;
      chi2[k][n1-1] = 0.0;
      chi12[k][0] = chi12[k][n1-1] = 0.0;
    }
    else if (udf.transform == 2) {
      chi1[k][0] = chi2[k][0] = 0.0;
    }
  }

  tmp = 100.0/counter;
  Message0("\n Scalar Dissipation summary:\n");
  Message0(" ----- cells updated -----  number ----- percent -----\n");
  Message0("\t total     \t %8i \t %6.2f\n", nZ1+nZ2+nZ,
           (double)(nZ1+nZ2+nZ)*tmp);
  Message0("\t Z1        \t %8i \t %6.2f\n", nZ1, (double)(nZ1)*tmp);
  Message0("\t Z2        \t %8i \t %6.2f\n", nZ2, (double)(nZ2)*tmp);
  Message0("\t Z         \t %8i \t %6.2f\n", nZ,  (double)(nZ)*tmp);
  tmp = 100.0/nZ;
  for (i=0; i<3; i++)
    Message0("\t chi%i big  \t %8i \t %6.2f --> %15.6e\n",i,largeChi[i],
             (double)(largeChi[i])*tmp,largeChiAvg[i]/largeChi[i]);
  tmp = 100.0/nZ; /* want percentages relative to cells with Z */
  Message0("\t Not cvgd. \t %8i \t %6.2f\n", nNotCvg,(double)(nNotCvg)*tmp);
  Message0("\t nZ1L      \t %8i \t %6.2f\n", nZ1L,(double)(nZ1L)*tmp);
  Message0("\t nZ2L      \t %8i \t %6.2f\n", nZ2L,(double)(nZ2L)*tmp);
  Message0(" -----------------------------------------------------\n");

  /* write the distribution */
  FILE *pFile;
#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    {
      pFile = fopen("chiCurrent.dat", "w");
      
      for (k=0; k<n2; k++) {
        for (j=0; j<n1; j++) {
          fprintf(pFile, "%8.6f %8.6f %15.6e %15.6e %15.6e\n",x1[j],x2[k],
                  chi1[k][j], chi12[k][j], chi2[k][j]);
        }
        fprintf(pFile, "\n");
      }
      
      fclose(pFile);
    }


  for (k=0; k<n2; k++) {
    free(Z1[k]); free(Z2[k]);
    free(f0[k]); free(f1[k]); free(f2[k]);
    free(pdf[k]); free(pdfBound[k]);
    free(chi0[k]); free(chi0Den[k]); free(chi1Den[k]); free(chi2Den[k]);
  }
  free(Z1); free(Z2); free(Z1min); free(Z1max); free(Z2min); free(Z2max);
  free(f0); free(f1); free(f2);
  free(pdf); free(pdfBound); free(pdfX1); free(pdfX2);
  free(chi0); free(chi0Den); free(chi1Den); free(chi2Den); free(rwork);
}
/* ------------------------------------------------------------------------- */
/* #endif /\* !RP_HOST *\/ */
/* ##### end functions that only work for serial/compute node ##### */

void UDFcalcChi2DAvg(Domain *d, int *nMesh, double **xMesh, double **chi1,
                     double **chi12, double **chi2)
{
  Message0("\n----- Computing 2D scalar dissipation rate ");
  Thread *t;
  cell_t c;

  int    i, j, k, iT;
  int    nZ, nZ1, nZ2, counter, nCount = 500, lineCount = 200 ;
  double TINY = 1.0e-6;
  double mass, tau, c_chi = 2.0;

  double chiMean[3], chiLocal[3], den[3];
  double **chi0, **chi0Den, **chi1Den, **chi2Den;
  double sigma1, sigma2;

  double **Z1, **Z2, *Z1min, *Z1max, *Z2min, *Z2max, Zplus, Zminus;
  double **f0, **f1, **f2;
  double tmp, *rwork;

  /* storage for the pdf */
  int    pdfBound1D[2], **pdfBound;
  double **pdf, *pdfX1, *pdfX2;

  /* dereference th mesh for ease of notation */
  int     n1 = nMesh[0];
  int     n2 = nMesh[1];
  double *x1 = xMesh[0];
  double *x2 = xMesh[1];
  double dx1_0 = x1[1] - x1[0];
  double dx2_0 = x2[1] - x2[0];

  double TINY_1 = 1.0e-5; /* 0.5*dx1_0; */
  double TINY_2 = 1.0e-5; /* 0.5*dx2_0; */

  char pdf2Mode[50];

  /* set functional form */
  if (udf.chiForm[2] == 1)
    UDFcalcChiForm = &MDUCcalcChiLogForm;
  else if (udf.chiForm[2] == 2)
    UDFcalcChiForm = &MDUCcalcChiErfcForm;
  else if (udf.chiForm[2] == 3)
    UDFcalcChiForm = &MDUCcalcChiQuadForm;

  /* allocate local memory */
  Z1       = (double**) calloc(n2, sizeof(double*));
  Z2       = (double**) calloc(n2, sizeof(double*));
  Z1min    = (double*)  calloc(n2, sizeof(double));
  Z1max    = (double*)  calloc(n2, sizeof(double));
  Z2min    = (double*)  calloc(n1, sizeof(double));
  Z2max    = (double*)  calloc(n1, sizeof(double));
  f0       = (double**) calloc(n2, sizeof(double*));
  f1       = (double**) calloc(n2, sizeof(double*));
  f2       = (double**) calloc(n2, sizeof(double*));
  pdf      = (double**) calloc(n2, sizeof(double*));
  pdfX1    = (double*)  calloc(n1, sizeof(double));
  pdfX2    = (double*)  calloc(n2, sizeof(double));
  pdfBound = (int**)    calloc(n2, sizeof(int*));
  chi0     = (double**) calloc(n2, sizeof(double*));
  chi0Den  = (double**) calloc(n2, sizeof(double*));
  chi1Den  = (double**) calloc(n2, sizeof(double*));
  chi2Den  = (double**) calloc(n2, sizeof(double*));
  rwork    = (double*)  calloc(MAX_I(n1,n2), sizeof(double));

  for (k=0; k<n2; k++)
    {
      Z1[k]       = (double*) calloc(n1, sizeof(double));
      Z2[k]       = (double*) calloc(n1, sizeof(double));
      f0[k]       = (double*) calloc(n1, sizeof(double));
      f1[k]       = (double*) calloc(n1, sizeof(double));
      f2[k]       = (double*) calloc(n1, sizeof(double));
      pdf[k]      = (double*) calloc(n1, sizeof(double));
      pdfBound[k] = (int*)    calloc(2, sizeof(int));
      chi0[k]     = (double*) calloc(n1, sizeof(double));
      chi0Den[k]  = (double*) calloc(n1, sizeof(double));
      chi1Den[k]  = (double*) calloc(n1, sizeof(double));
      chi2Den[k]  = (double*) calloc(n1, sizeof(double));
    }

  /* set chi to zero everywhere */
  for (k=0; k<n2; k++)
    {
      for (j=0; j<n1; j++)
        {
          pdf[k][j]  = 0.0;
          chi1[k][j] = chi2[k][j] = chi2[k][j] = chi0[k][j] = 0.0;
          chi0Den[k][j] = chi1Den[k][j] = chi2Den[k][j] = 0.0;
        }
    }

  /* fill the non-transformed value grid using the appropriate transformation */
  for (k=0; k<n2; k++)
    {
      Z1min[k] = 1.0;
      Z1max[k] = 0.0;

      for (j=0; j<n1; j++)
        {
          Z2min[j] = 1.0;
          Z2max[j] = 0.0;

          /* apply transformation */
          if (udf.transform == 0)
            {
              Z1[k][j] = x1[j];
              Z2[k][j] = x2[k];
            }
          else if (udf.transform == 1)
            {
              Z1[k][j] = x1[j];
              Z2[k][j] = x2[k]*(1.0-x1[j]);
            }
          else if (udf.transform == 2)
            {
              Z1[k][j] = x1[j]*(1.0 - x2[k]);
              Z2[k][j] = x2[k]*x1[j];
            }
        } /* end loop over columns */
    } /* end loop over rows */

  int oneZ1 = 0;
  int oneZ2 = 0;
  double beta;

  /* find the global max in the domain */
  for (iT=0; iT<udf.nThreadSolve; iT++)
    {
      t = udf.solve_threads[iT]; /* current thread */
      if (t == NULL)
        continue;

      begin_c_loop_int(c,t)
        {
          /* use this if computing global Z1max in the domain */
          if (C_PHI1(c,t) > 0.0)
            {

              beta  = C_PHI1(c,t)*(1.0-C_PHI1(c,t))/(C_PHIV1(c,t)+1e-6) - 1;
              beta *= (1.0-C_PHI1(c,t));
              if (beta < 1.0) oneZ1 = 1;

              sigma1 = pow(MAX_D(C_PHIV1(c,t),0), 0.5);

              /* find where to condition */
              k = UDFgridLocate(n2, x2, C_PHI2(c,t));
              if (C_PHI2(c,t) > 0.5*(x2[k]+x2[k+1]))
                k = k+1;

              Zplus  = C_PHI1(c,t) + 2*sigma1;
              Zminus = C_PHI1(c,t) - 2*sigma1;
              if (Zminus < Z1min[k]) Z1min[k] = MAX_D(Zminus, 0.0);
              if (Zplus  > Z1max[k]) Z1max[k] = MIN_D(Zplus,  1.0);
            }

          /* use this if computing global Z2max in the domain */
          if (C_PHI2(c,t) > 0.0)
            {

              beta  = C_PHI2(c,t)*(1.0-C_PHI2(c,t))/(C_PHIV2(c,t)+1e-6) - 1;
              beta *= (1.0-C_PHI2(c,t));
              if (beta < 1.0) oneZ2 = 1;

              sigma2 = pow(MAX_D(C_PHIV2(c,t),0), 0.5);

              /* find where to condition */
              j = UDFgridLocate(n1, x1, C_PHI1(c,t));
              if (C_PHI1(c,t) > 0.5*(x1[j]+x1[j+1]))
                j = j+1;

              Zplus  = C_PHI2(c,t) + 2*sigma2;
              Zminus = C_PHI2(c,t) - 2*sigma2;
              if (Zminus < Z2min[j]) Z2min[j] = MAX_D(Zminus, 0.0);
              if (Zplus  > Z2max[j]) Z2max[j] = MIN_D(Zplus,  1.0);
            }
        }
      end_c_loop_int(c,t) /* end loop over all cells */
    } /* end loop over cell threads */

#if RP_NODE
  /* accumulate the min*/
  PRF_GRLOW(Z1min, n2, rwork);
  PRF_GRHIGH(Z1max, n2, rwork);
  PRF_GRLOW(Z2min, n1, rwork);
  PRF_GRHIGH(Z2max, n1, rwork);
#endif /* RP_NODE */

  double Z2maxAll = 0.0;
  for (j=0; j<n1; j++)
    {
      if (Z2max[j] > Z2maxAll)
        Z2maxAll = Z2max[j];
    }
  double Z1maxAll = 0.0;
  for (k=0; k<n2; k++)
    {
      if (Z1max[k] > Z1maxAll)
        Z1maxAll = Z1max[k];
    }

  double Z0maxAll = 1-Z1maxAll-Z2maxAll;

  Message0(" Max Mixture fractions:\n");
  Message0("   Z0: %8.6f\n", Z0maxAll);
  Message0("   Z1: %8.6f\n", Z1maxAll);
  Message0("   Z2: %8.6f\n", Z2maxAll);
  Message0("   Z1_1 --> %i\n", oneZ1);
  Message0("   Z2_1 --> %i\n", oneZ2);
  
  char buffer[100];
  Message0("starting loop computation...\n");
  counter = nZ = nZ1 = nZ2 = 0;
  for (iT=0; iT<udf.nThreadSolve; iT++)
    {
      t = udf.solve_threads[iT]; /* current thread */
      if (t == NULL)
        continue;

      begin_c_loop_int(c,t)
        {
          /* ensure all storage memory is reset to zero */
          C_CHI0_REF(c,t) = C_CHI1_REF(c,t) = C_CHI2_REF(c,t) = 0.0;
          for (k=0; k<n2; k++)
            for (j=0; j<n1; j++)
              f0[k][j] = f1[k][j] = f2[k][j] = 0.0;

          if ( (C_PHI1(c,t) > TINY_1) )
            {
              mass = C_R(c,t)*C_VOLUME(c,t);   /* cell mass */
              tau  = C_D(c,t)/(C_K(c,t)+TINY); /* turb. eddy turnover time */

              /* check if the cell contains both mixture fractions */
              if (C_PHI2(c,t) > TINY_2)
                {
                  /* compute mean scalar dissipation rates */
                  chiMean[0] = c_chi*tau*C_PHIV0(c,t);
                  chiMean[1] = c_chi*tau*C_PHIV1(c,t);
                  chiMean[2] = c_chi*tau*C_PHIV2(c,t);

                  /* set the default joint pdf treatment */
                  strcpy(pdf2Mode, udf.pdf2D);

                  if ( C_PHI1(c,t) < 2*dx1_0 || C_PHI2(c,t) < 2*dx2_0 )
                    strcpy(pdf2Mode, "dirichlet");

                  /* compute 2D probability distribution */
                  MDUCcalcPDF_2D(pdf2Mode, nMesh, xMesh, C_PHI1(c,t),
                                 C_PHIV1(c,t), C_PHI2(c,t), C_PHIV2(c,t),
                                 C_PHIV0(c,t), pdf, pdfBound, udf.transform);

                  den[0] = den[1] = den[2] = 0.0;
                  for (k=0; k<n2; k++)
                    {
                      for (j=0; j<n1; j++)
                        {
                          if (Z2[k][j] <= Z2maxAll*(1.0-Z1[k][j]/Z1maxAll))
                            {
                              f0[k][j] =
                                UDFcalcChiForm(Z1[k][j]/Z1maxAll+
                                               Z2[k][j]/Z2maxAll,
                                               C_PHI1(c,t)/Z1maxAll+
                                               C_PHI2(c,t)/Z2maxAll,
                                               1);
                              f1[k][j] =
                                UDFcalcChiForm(Z1[k][j],C_PHI1(c,t),Z1maxAll);
                              f2[k][j] =
                                UDFcalcChiForm(Z2[k][j],C_PHI2(c,t),Z2maxAll);

                              den[0] += f0[k][j]*pdf[k][j];
                              den[1] += f1[k][j]*pdf[k][j];
                              den[2] += f2[k][j]*pdf[k][j];
                            }
                        } /* end loop over columns */
                    } /* end loop over rows */ 

                  if (den[0] > TINY) C_CHI0_REF(c,t) = chiMean[0]/den[0];
                  if (den[1] > TINY) C_CHI1_REF(c,t) = chiMean[1]/den[1];
                  if (den[2] > TINY) C_CHI2_REF(c,t) = chiMean[2]/den[2];

                  for (k=0; k<n2; k++)
                    {
                      for (j=0; j<n1; j++)
                        {
                          chiLocal[0] = C_CHI0_REF(c,t)*f0[k][j];
                          chiLocal[1] = C_CHI1_REF(c,t)*f1[k][j];
                          chiLocal[2] = C_CHI2_REF(c,t)*f2[k][j];

                          chi0[k][j]    += mass*chiLocal[0];
                          chi0Den[k][j] += mass;
                          chi1[k][j]    += mass*chiLocal[1];
                          chi1Den[k][j] += mass;
                          chi2[k][j]    += mass*chiLocal[2];
                          chi2Den[k][j] += mass;
                        } /* end loop over columns */
                    } /* end loop over rows */ 
                  nZ++;
                }
              else
                {
                  /* compute mean scalar dissipation rates */
                  chiMean[1] = c_chi*tau*C_PHIV1(c,t);

                  /* this cell only has Z1 */
                  MDUCcalcPDF_1D(udf.pdf1D, n1, x1, C_PHI1(c,t),C_PHIV1(c,t),
                                 pdfX1, pdfBound1D);

                  /* compute f1(x) based on assumed functional form */
                  den[1] = 0.0;
                  for (j=0; j<n1; j++)
                    {
                      /* compute the functional form for F1(xi1, xi2)
                         and F2(xi1,xi2) chiMean_i = chi_st_i*F_i(xi1, xi2)
                         uses inverse error function */
                      f1[0][j] = UDFcalcChiForm(x1[j],C_PHI1(c,t),Z1maxAll);

                      /* compute contribution to integral */
                      den[1] += f1[0][j]*pdfX1[j];
                    } /* end loop over columns */

                  /* compute the conditional mean chi at stoichiometric */
                  if (den[1] > TINY) C_CHI1_REF(c,t) = chiMean[1]/den[1];

                  for (j=0; j<n1; j++)
                    {
                      /* chi in normalised coordinates */
                      chiLocal[1] = C_CHI1_REF(c,t)*f1[0][j];

                      /* add contributions to global chi */
                      chi0[0][j]    += mass*chiLocal[1];
                      chi0Den[0][j] += mass;
                      chi1[0][j]    += mass*chiLocal[1];
                      chi1Den[0][j] += mass;
                    } /* end loop over columns */
                  nZ1++;
                }

            } /* end Z1 */
          else if (C_PHI2(c,t) > TINY_2)
            {
              /* cell only has Z2 */
              mass = C_R(c,t)*C_VOLUME(c,t);     /* cell mass */
              tau  = C_D(c,t)/(C_K(c,t)+TINY);   /* turb. eddy turnover time */
              chiMean[2] = c_chi*tau*C_PHIV2(c,t);

              /* this cell only has Z2, compute 1D distribution */
              MDUCcalcPDF_1D(udf.pdf1D, n2, x2, C_PHI2(c,t), C_PHIV2(c,t),
                             pdfX2, pdfBound1D);

              /* compute f2(y) based on assumed functional form */
              den[2] = 0.0;
              for (k=0; k<n2; k++)
                {
                  /* compute functional form */
                  f2[k][0] = UDFcalcChiForm(x2[k],C_PHI2(c,t),Z2maxAll);

                  /* compute contribution to integral */
                  den[2] += f2[k][0]*pdfX2[k];
                } /* end loop over rows */

              /* compute the conditional mean chi at stoichiometric */
              if (den[2] > TINY)  C_CHI2_REF(c,t) = chiMean[2]/den[2];

              for (k=0; k<n2; k++)
                {
                  /* chi in normalised coordinates */
                  chiLocal[2] = C_CHI2_REF(c,t)*f2[k][0];

                  /* add contributions to global chi */
                  chi0[k][0]    += mass*chiLocal[2];
                  chi0Den[k][0] += mass;
                  chi2[k][0]    += mass*chiLocal[2];
                  chi2Den[k][0] += mass;
                } /* end loop over rows */
              nZ2++;
            }

          counter++;
        }
      end_c_loop_int(c,t) /* end loop over all cells */
    } /* end loop over cell threads */
  Message0("Finished loop over cells\n");

#if RP_NODE
  nZ      = PRF_GISUM1(nZ);
  nZ1     = PRF_GISUM1(nZ1);
  nZ2     = PRF_GISUM1(nZ2);
  counter = PRF_GISUM1(counter);

  for (k=0; k<n2; k++)
    {
      PRF_GRSUM(chi0[k], n1, rwork);
      PRF_GRSUM(chi0Den[k], n1, rwork);
      PRF_GRSUM(chi1[k], n1, rwork);
      PRF_GRSUM(chi1Den[k], n1, rwork);
      PRF_GRSUM(chi2[k], n1, rwork);
      PRF_GRSUM(chi2Den[k], n1, rwork);
    }
#endif /* RP_NODE */
 
  /* normalise the distribution */
  for (k=0; k<n2; k++)
    {
      for (j=0; j<n1; j++)
        {
          if (chi0Den[k][j] > TINY)
            chi0[k][j] = chi0[k][j]/chi0Den[k][j];
          else
            chi0[k][j] = 0.0;

          if (chi1Den[k][j] > TINY)
            chi1[k][j] = chi1[k][j]/chi1Den[k][j];
          else
            chi1[k][j] = 0.0;

          if (chi2Den[k][j] > TINY)
            chi2[k][j] = chi2[k][j]/chi2Den[k][j];
          else
            chi2[k][j] = 0.0;

          chi12[k][j] = 0.5*(chi0[k][j] - chi1[k][j] - chi2[k][j]);
        } /* end loop over columns */
    } /* end loop over rows */

  /* enforce zero boundary conditions as appropriate */
  for (k=0; k<n2; k++)
    {
      if (udf.transform == 1)
        {
          chi1[k][0] = chi1[k][n1-1] = 0.0;
          chi2[k][n1-1] = 0.0;
          chi12[k][0] = chi12[k][n1-1] = 0.0;
        }
      else if (udf.transform == 2)
        {
          chi1[k][0] = chi2[k][0] = 0.0;
        }
    }

  tmp = 100.0/counter;
  Message0("\n Scalar Dissipation summary:\n");
  Message0(" ----- cells updated -----  number ----- percent -----\n");
  Message0("\t total     \t %8i \t %6.2f\n", nZ1+nZ2+nZ,
           (double)(nZ1+nZ2+nZ)*tmp);
  Message0("\t Z1        \t %8i \t %6.2f\n", nZ1, (double)(nZ1)*tmp);
  Message0("\t Z2        \t %8i \t %6.2f\n", nZ2, (double)(nZ2)*tmp);
  Message0("\t Z         \t %8i \t %6.2f\n", nZ,  (double)(nZ)*tmp);
  Message0(" -----------------------------------------------------\n");

 for (k=0; k<n2; k++)
    {
      free(Z1[k]); free(Z2[k]);
      free(f0[k]); free(f1[k]); free(f2[k]);
      free(pdf[k]); free(pdfBound[k]);
      free(chi0[k]); free(chi0Den[k]); free(chi1Den[k]); free(chi2Den[k]);
    }
  free(Z1); free(Z2); free(Z1min); free(Z1max); free(Z2min); free(Z2max);
  free(f0); free(f1); free(f2);
  free(pdf); free(pdfBound); free(pdfX1); free(pdfX2);
  free(chi0); free(chi0Den); free(chi1Den); free(chi2Den); free(rwork);
}
/* ------------------------------------------------------------------------- */
#endif /* !RP_HOST */
/* ##### end functions that only work for serial/compute node ##### */
