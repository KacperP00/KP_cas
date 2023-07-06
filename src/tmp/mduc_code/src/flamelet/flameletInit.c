#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "flameletInit.h"
#include "flameletUtils.h"
#include "init.h"
#include "flamelet1D.h"
#include "flamelet2D.h"
#include "solve1D.h"
#include "utils.h"
#include "metric.h"
#include "transport.h"
#include "soot.h"
#include "chemMech.h"
#include "chemUtils.h"
#include "solveUtils.h"
#include "display.h"
#include "fileIO.h"
#include "sootHMOM.h"
#include "rk4solver.h"

flameStruct *flameletInit(ioStruct *io, driveStruct *driver, chemStruct *chem,
                          flowStruct *flow)
{
  int i, j, k, k_pre, r, row, col, ii, maxDim, nGrid, nEq, iX1;
  double oppSideY, *initH, H_ref, *energy;

  flameStruct *flame = NULL;

  /* get the flamelet from the I/O */
  flame = paramGetFlamelet("1",driver,chem,flow,io->parser);

  /* set global defaults */
  flame->nFlamelets = 1; /* solve at least one flamelet */

  /* read in additional global flamelet parameters */
  flame->order       = paramGetInt("spatial-order",2,io->parser);
  driver->PHYS_COORD = paramGetInt("physical-coordinate",0,io->parser);

  /* ----- initialise parameters from chemical mechanism ----- */
  flame->nVars = chem->nSpecSolve+1;
  flame->sT    = chem->nSpecSolve;

  if (driver->CHI_MODE == 2) { /* solve equation for chi */
    flame->sChi   = flame->nVars;
    flame->nVars += 1;
  }
  if (chem->sootMode == 1) { /* solve equations for soot moments */
    flame->sSoot  = flame->nVars;
    flame->nVars += chem->nSootMom;
    //fprintf(stderr, "adding soot moments to solution: %i\n", chem->nSootMom);
  }

  /* assign some indices for specific species */
  flame->sO2 = getSpeciesIndex(chem->nSpecies, chem->species, "O2");
  flame->sN2 = getSpeciesIndex(chem->nSpecies, chem->species, "N2");

  maxDim = IMAX(flame->n1, flame->n2);
  nGrid  = flame->n1*flame->n2;
  nEq    = nGrid*flame->nVars;
  flame->j0 = 0;
  flame->jN = flame->n1-1;
  flame->k0 = 0;
  flame->kN = flame->n2-1;

  /* solve additional equations for engine */
  if (driver->ENGINE) {
    if (driver->PRES_MODE == 3) { /* implicit pressure */
      flame->sP  = nEq;
      nEq       += 1;
    }
    flame->sM  = nEq;
    flame->sMa = nEq; nEq += 1;   /* air mass */
    flame->sMf = nEq; nEq += 1;   /* fuel mass */
    flame->sMe = nEq; nEq += 1;   /* egr mass */
    if (driver->ENGINE > 1) { /* SI */
      flame->sMb = nEq; nEq += 1; /* burned mass */
      flame->sTb = nEq; nEq += 1; /* burned zone temperature */
      flame->sDF = nEq; nEq += 1; /* flame diameter */
    }
  }

  /* DEBUG: add reaction rates */
  if (driver->output_w) {
     flame->sRxn = nEq;
  }
  /* ----- allocate solution ----- */
  flame->Y     = (double*) calloc(nEq, sizeof(double));
  flame->Yp    = (double*) calloc(nEq, sizeof(double));
  flame->Y_m1  = (double*) calloc(nEq, sizeof(double));
  flame->Y_m2  = (double*) calloc(nEq, sizeof(double));
  flame->Y_pre = (double*) calloc(nEq, sizeof(double));
  flame->sumYi = (double*) calloc(nEq, sizeof(double));
  flame->Yp_m1 = (double*) calloc(nEq, sizeof(double));
  flame->chemSrcY     = (double*) calloc(nEq, sizeof(double));
  flame->chemSrcY_tmp = (double*) calloc(nEq, sizeof(double));

  if (flame->sP >= 0)
    flame->Y[flame->sP] = flow->P;

  /* ----- allocate local properties ----- */
  flame->chiA     = (double*) calloc(nGrid, sizeof(double));
  flame->rhoDotSoot   = (double*) calloc(nGrid, sizeof(double));
  flame->rhoDotSpray      = (double*) calloc(nGrid, sizeof(double));
  flame->xi          = (double*) calloc(nGrid, sizeof(double));
  flame->D           = (double*) calloc(nGrid, sizeof(double));
  flame->rhoSigma    = (double*) calloc(nGrid, sizeof(double));
  flame->v1       = (double*) calloc(nGrid, sizeof(double));
  flame->vT       = (double*) calloc(nGrid, sizeof(double));
  if (chem->soot) {
     flame->vMxmyp          = (double*) calloc(nGrid, sizeof(double));
     flame->vMxy          = (double*) calloc(nGrid, sizeof(double));
  }
  if (flame->rho == NULL)
  flame->rho      = (double*) calloc(nGrid, sizeof(double));
  flame->T        = (double*) calloc(nGrid, sizeof(double));
  flame->gamma    = (double*) calloc(nGrid, sizeof(double));
  flame->cpmix    = (double*) calloc(nGrid, sizeof(double));
  flame->Mw       = (double*) calloc(nGrid, sizeof(double));
  flame->R        = (double*) calloc(nGrid, sizeof(double));
  flame->mu       = (double*) calloc(nGrid, sizeof(double));
  if (chem->soot) 
  flame->D        = (double*) calloc(nGrid, sizeof(double));

  if (driver->detailedDiff)
     flame->LeZ_LeI = (double*) calloc(nGrid*chem->nSpecSolve, sizeof(double));

  flame->srcP     = (double*) calloc(nGrid, sizeof(double));
  flame->srcH     = (double*) calloc(nGrid, sizeof(double));
  flame->srcQ     = (double*) calloc(nGrid, sizeof(double));
  flame->sumCp    = (double*) calloc(nGrid, sizeof(double));
  flame->HRrate   = (double*) calloc(nGrid, sizeof(double));
  flame->visc_src = (double*) calloc(nGrid, sizeof(double));
  flame->radiation = (double*) calloc(nGrid, sizeof(double));
  flame->gPDF = (double*) calloc(nGrid, sizeof(double));

  flame->convCoeff = (double*) calloc(nGrid, sizeof(double));
  if (driver->output_tsrc) {
    flame->conv_src  = (double*) calloc(nGrid, sizeof(double));
    flame->chem_src  = (double*) calloc(nGrid, sizeof(double));
    flame->diff_src  = (double*) calloc(nGrid, sizeof(double));
    flame->conv_temp = (double*) calloc(nGrid, sizeof(double));
    flame->chem_temp = (double*) calloc(nGrid, sizeof(double));
    flame->diff_temp = (double*) calloc(nGrid, sizeof(double));
  }

  if (chem->soot) {
    driver->sootMode = chem->sootMode;
    if (chem->sootMode < 3) {
       flame->nSootMom = chem->nSootMom;
       flame->sootMom    = (double*) calloc(flame->nSootMom*nGrid, sizeof(double));
       flame->sootMomxmyp = (double*) calloc((flame->nSootMom-1)*nGrid, sizeof(double));
       flame->sootMomSrc = (double*) calloc(flame->nSootMom*nGrid, sizeof(double));
       flame->sootMomSrcPart =
         (double*) calloc(chem->nSootSrc*flame->nSootMom*nGrid, sizeof(double));
       flame->sootMomTrans = (double*) calloc(flame->nSootMom*nGrid, sizeof(double));
       flame->sootMomTransConv = (double*) calloc(flame->nSootMom*nGrid, sizeof(double));
       flame->sootMomTransDiff = (double*) calloc(flame->nSootMom*nGrid, sizeof(double));
    }
    if (chem->sootMode >= 3) {
       flame->sqrtT    		= (double*) calloc(nGrid, sizeof(double));
       flame->T_mu    		= (double*) calloc(nGrid, sizeof(double));
       flame->muSqrtMw_rhoSqrtT = (double*) calloc(nGrid, sizeof(double));
       flame->sgCoeff    	= (double*) calloc(nGrid, sizeof(double));
       flame->oxCoeff  		= (double*) calloc(nGrid, sizeof(double));
       flame->o2Coeff  		= (double*) calloc(nGrid, sizeof(double));
       flame->dimProdRate  	= (double*) calloc(nGrid, sizeof(double));
  }
  }

  /* if flamelet file was read, place data in appropriate format */
  if (io->readData) {
    /* compute dx according to grid read in */
    if (flame->n1 > 1) {
      flame->dx1 = (double*) calloc(flame->n1-1, sizeof(double));
      for (j=0; j<flame->n1-1; j++)
        flame->dx1[j] = flame->x1[j+1] - flame->x1[j];
    }
    
    // set the whole solution vector to zero
    for (i=0; i<flame->nVars*flame->n1*flame->n2; ++i)
      flame->Y[i] = 0.0;
    
    // set the whole solution vector to zero
    for (i=0; i<flame->nVars*flame->n1*flame->n2; ++i) {
      flame->Y[i] = 0.0;
      flame->chemSrcY[i] = 0.0;
      flame->chemSrcY_tmp[i] = 0.0;
    }
    
    if (driver->ENERGY) {
      if (!flame->H) { /* file only had temperature */
        /* compute enthalpy */
        flame->H = (double*) calloc(nGrid, sizeof(double));
        
        double *cp = (double*) calloc(chem->nSpecies, sizeof(double));
        double *h  = (double*) calloc(chem->nSpecies, sizeof(double));
        
        for (k=0; k<flame->n2; k++) {
          r   = k*flame->n1;
          for (j=0; j<flame->n1; j++) {
            chemComputeThermoData(h,cp,io->initT[r+j],chem);
            flame->H[r+j] = 0.0;
            for (i=0; i<chem->nSpecSolve; i++)
              flame->H[r+j] += h[i]*io->initY[i][r+j];
          } /* end loop over columns */
        } /* end loop over rows */
        
        free(cp); free(h);
      }
      for (i=0; i<nGrid; i++) flame->T[i] = io->initT[i];
      energy = flame->H;
    }
    else
      energy = io->initT;
      if (!flame->H)  /* file only had temperature */
        /* compute enthalpy */
        flame->H = (double*) calloc(nGrid, sizeof(double));

    /* initialise the solution vector based on the species
       read from the input file */
    for (k=0; k<flame->n2; ++k) {
      r   = k*flame->n1;
      row = k*flame->n1*flame->nVars;
      for (j=0; j<flame->n1; ++j) {
        col = j*flame->nVars;
        double sum=0.0;
        for (i=0; i<io->nInputs; ++i) {
#ifdef SUNDIALS_FROM_3_2
	     // TODO: Prohibit negative values if contraints are enforced
	     // Negative initial values mostly likely caused by importing 
	     // text file rather that binary
             flame->Y[row+col+io->inputIndex[i]] = fmax(io->initY[i][r+j],0.0);
#else
             flame->Y[row+col+io->inputIndex[i]] = io->initY[i][r+j];
#endif 
          sum += flame->Y[row+col+io->inputIndex[i]];
        }
        if (sum<0.95) 
	  fprintf(stderr,"Sum of mass fraction below 0.95: %15.6e %i\n",sum,j);
        
        flame->Y[row+col+flame->sT] = energy[r+j];//io->initT[r+j];
      } /* end loop over columns */
    } /* end loop over rows */

    if (chem->soot) {
      if (chem->sootMode <= 3) {
         for (k=0; k<flame->n2; k++) {
           r   = k*flame->n1;
           row = k*flame->n1*flame->nSootMom;
           for (j=0; j<flame->n1; j++) {
             col = j*flame->nSootMom;
             flame->sootMom[row+col+0] = io->initSoot[0][r+j];
             flame->sootMom[row+col+1] = io->initSoot[1][r+j];
             flame->sootMom[row+col+2] = io->initSoot[2][r+j];
             if (flame->nSootMom > 4) {
               flame->sootMom[row+col+3] = io->initSoot[3][r+j];
               flame->sootMom[row+col+4] = io->initSoot[4][r+j];
               flame->sootMom[row+col+5] = io->initSoot[5][r+j];
               flame->sootMom[row+col+6] = io->initSoot[6][r+j];
             }
             else
               flame->sootMom[row+col+3] = io->initSoot[6][r+j];

             /* Make sure that the minimum values are not smaller than the
 	      * the solver tolerance */ 
             reinitSootMoments(flame->nSootMom,&flame->sootMom[row+col],1.0E-20);
             //reinitSootMoments(flame->nSootMom,&flame->sootMom[row+col],driver->SABSTOL);

             if (chem->sootMode == 1) {
               for (i=0; i<flame->nSootMom; i++) {
                 flame->Y[(k*flame->n1+j)*flame->nVars+flame->sSoot+i] = 
                   flame->sootMom[row+col+i];
               } /* end loop over soot moments */
             } /* implicit soot mode */

           } /* end loop over columns */
         } /* end loop over rows */
      } /* soot mode <= 3 */
      for (k=0; k<flame->n2; k++) {
        r   = k*flame->n1;
        row = k*flame->n1*flame->nSootMom;
        for (j=0; j<flame->n1; j++) {
          /* make sure viscosity is up-to-date */
          flame->mu[r+j] = calcMixViscosity(chem->nSpecies,
                                            flame->Y[j*flame->nVars+flame->sT],
                                            &flame->Y[j*flame->nVars], chem);
        } /* end loop over columns */
      } /* end loop over rows */
    } /* end soot */
    
    flame->nSootMom = chem->nSootMom;

    if (flame->rho == NULL) 
      flame->rho = (double*) calloc(nGrid, sizeof(double));
    
    // allocate memory specific to 2D flamelet
    if (flame->nDim == 2) {
      // compute dx according to grid read in
      flame->dx2 = (double*) calloc(flame->n2-1, sizeof(double));
      for (j=0; j<flame->n2-1; j++)
        flame->dx2[j] = flame->x2[j+1] - flame->x2[j];
      
      flame->chiB    = (double*) calloc(nGrid, sizeof(double));
      flame->chiC    = (double*) calloc(nGrid, sizeof(double));
      flame->v2      = (double*) calloc(nGrid, sizeof(double));
      flame->sumCp2  = (double*) calloc(nGrid, sizeof(double));
      flame->dY.y    = (double*) calloc(flame->nVars*maxDim,sizeof(double));
      flame->dY.yy   = (double*) calloc(flame->nVars*maxDim,sizeof(double));
      flame->dY.xy   = (double*) calloc(flame->nVars*maxDim,sizeof(double));
      flame->dCp.y   = (double*) calloc(maxDim, sizeof(double));
    }
    
    /* initialise the metric quantities */
    if (flame->nDim > 0)
      initMetric(flame);

    /* compute the current properties */
    double wmix;
    for (k=0; k<flame->n2; k++) {
      row = k*flame->n1*flame->nVars;
      r   = k*flame->n1;
      for (j=0; j<flame->n1; j++) {
        col  = j*flame->nVars;
        wmix = calcMixtureW(chem->nSpecSolve, &flame->Y[row+col], chem->W);
        flame->Mw[r+j] = wmix;
        flame->rho[r+j] = wmix*flow->P/(flame->RGAS*flame->Y[row+col+flame->sT]);
        if (chem->readTrans)
          flame->mu[r+j]  = calcMixViscosity(chem->nSpecies,
                                             flame->Y[row+col+flame->sT],
                                             &flame->Y[row+col], chem);
      } /* end loop over columns */
    } /* end loop over rows */
    
    if (chem->soot) {
//       if (chem->sootMode == 1) {
//         for (k=0; k<flame->n2; k++) {
//           r   = k*flame->n1;
//           row = k*flame->n1*flame->nSootMom;
//           for (j=0; j<flame->n1; j++) {
//               for (i=0; i<flame->nSootMom; i++) {
//                 flame->Y[(k*flame->n1+j)*flame->nVars+flame->sSoot+i] /= flame->rho[r+j]; 
//               } /* end loop over soot moments */
//             } /* end grid loop 1 */
//	 } /* end grid loop 2 */
//       }
       if (chem->sootMode >= 3) {
          for (k=0; k<flame->n2; k++) {
            row = k*flame->n1*flame->nVars;
            r   = k*flame->n1;
            for (j=0; j<flame->n1; j++) {
              col  = j*flame->nVars;
              /*TODO: Compute concentration only for PAH */
              /* compute and store concentrations of all species at local point */
              for (i = 0; i < chem->nSpecies; i++)
                chem->c[i] = flame->rho[r+j]*flame->Y[row+col+i]/chem->W[i];
              /* copy the species necessary for soot reactions */
              int iS;
              for (iS=0; iS<chem->nSootSp; iS++)
             chem->sootRxnConc[iS] = chem->c[chem->iSootSp[iS]];

	      double sgCoeff, oxCoeff,o2Coeff, dimProdRate;
              sootHMOM_GasPhase_Coeff(flame->Y[row+col+flame->sT], flame->mu[r+j],
              	                flame->rho[r+j], flame->Mw[r+j], chem->sootRxnConc, 
              			&sgCoeff, &oxCoeff,&o2Coeff, &dimProdRate);

	      flame->sgCoeff[r+j] 	    = sgCoeff;
	      flame->oxCoeff[r+j]  	    = oxCoeff;
	      flame->o2Coeff[r+j] 	    = o2Coeff;
	      flame->dimProdRate[r+j] 	    = dimProdRate;
              flame->sqrtT[r+j] 	    = sqrt(flame->Y[row+col+flame->sT]);
              flame->T_mu[r+j]  	    = flame->Y[row+col+flame->sT]/flame->mu[r+j];
              flame->muSqrtMw_rhoSqrtT[r+j] = flame->mu[r+j]*sqrt(flame->Mw[r+j]/
              			(1000.0*flame->Y[row+col+flame->sT]*flame->rho[r+j]*flame->rho[r+j]));
            } 
          }
       }
    }
    if (driver->output_w)
      for (i=0; i<nGrid; i++) flame->gPDF[i] = io->gPDF[i];

    return flame;
  } // end init from data file
 
  if (driver->output_w) {
     for (i=0; i<nGrid; i++) flame->gPDF[i] = 1.0;
  }
    

  if (driver->output_w)
      for (i=0; i<nGrid; i++) flame->gPDF[i] = 1.0;

  /* initialise the grid in the x1 direction */
  initGrid(flame->n1,0,1,flame->nPiv1,flame->gPiv1,flame->gPc1,
           &flame->x1,&flame->dx1);

  iX1 = flame->iX10 = (flame->n1-1)*flame->nVars;

  if (io->outputHtot && flame->H == NULL)
    flame->H = (double*) calloc(nGrid, sizeof(double));
  
  /* get the species indexing for the boundaries */
  bndStruct *ox    = flame->ox;
  bndStruct *fuel1 = flame->fuel1;
  bndStruct *fuel2 = flame->fuel2;
  if (ox != NULL)
    getSpeciesIndices(ox->nSp,ox->sp,ox->iSp,chem->nSpecies,chem->species);
  if (fuel1 != NULL)
    getSpeciesIndices(fuel1->nSp,fuel1->sp,fuel1->iSp,chem->nSpecies,
                      chem->species);
  if (fuel2 != NULL)
    getSpeciesIndices(fuel2->nSp,fuel2->sp,fuel2->iSp,chem->nSpecies,
                      chem->species);

  fprintf(stderr, "Oxidizer boundary:\n");
  for (i=0; i<ox->nSp; i++) {
    fprintf(stderr,"%25s %5i %15.6e\n",ox->sp[i],ox->iSp[i],ox->Y[i]);
  }

  fprintf(stderr, "Fuel boundary:\n");
  for (i=0; i<fuel1->nSp; i++) {
    fprintf(stderr,"%25s %5i %15.6e\n",fuel1->sp[i],fuel1->iSp[i],fuel1->Y[i]);
  }

  /* if 0D, interpolate the species from the mass fractions based on the mean */
  if (flame->nDim == 0 && flow->phiMean[1] >= 0.0) {

    double *Y_Z0 = (double*) calloc(chem->nSpecies, sizeof(double));
    double *Y_Z1 = (double*) calloc(chem->nSpecies, sizeof(double));
    
    for (i=0; i<chem->nSpecies; i++)
      Y_Z0[i] = Y_Z1[i]  = 0.0; 
    
    /* fill in the mass fractions */
    for (i=0; i<ox->nSp; i++)    Y_Z0[ox->iSp[i]]    = ox->Y[i];
    for (i=0; i<fuel1->nSp; i++) Y_Z1[fuel1->iSp[i]] = fuel1->Y[i];
    
    if (flow->phiMean[2] >= 0) {
      /* first interpolate the streams */
      double *Y_N  = (double*) calloc(chem->nSpecies, sizeof(double));
      double *Y_Z2 = (double*) calloc(chem->nSpecies, sizeof(double));
      double x1_N, T_N;

      for (i=0; i<chem->nSpecies; i++)
        Y_Z2[i] = Y_N[i] = 0.0; 
      
      for (i=0; i<fuel2->nSp; i++)
        Y_Z2[fuel2->iSp[i]] = fuel2->Y[i];

      /* intersection of mixing line with x1 axis */
      if (flow->phiMean[2] < 1.0)
        x1_N = flow->phiMean[1]/(1.0-flow->phiMean[2]);
      else 
        x1_N = 0.0;

      /* first mix the oxidizer side with the fuel 1 side */
      for (i=0; i<chem->nSpecies; i++)
        Y_N[i] = Y_Z0[i]*(1-x1_N) + Y_Z1[i]*x1_N;

      T_N = ox->T*(1-x1_N) + fuel1->T*x1_N;

      /* interpolate along mixing line to fuel 2 */
      for (i=0; i<chem->nSpecies; i++)
        flame->Y[i] = (1.0-flow->phiMean[2])*Y_N[i] + flow->phiMean[2]*Y_Z2[i];

      flame->Y[flame->sT] = (1.0-flow->phiMean[2])*T_N 
        + flow->phiMean[2]*fuel2->T;

      free(Y_N); free(Y_Z2);
    }
    else { /* 1D interpolation */
      for (i=0; i<chem->nSpecies; i++)
        flame->Y[i] = Y_Z0[i]*(1-flow->phiMean[1]) + Y_Z1[i]*flow->phiMean[1];

      /* linearly mix the temperature */
      flame->Y[flame->sT] = ox->T*(1-flow->phiMean[1])
        + fuel1->T*flow->phiMean[1];
    }

    fprintf(stderr, "initial species info:\n");
    for (i=0; i<chem->nSpecies; i++)
      if (flame->Y[i] > 1.0e-7)
        fprintf(stderr, " %25s: %15.6e\n", chem->species[i], flame->Y[i]);
    
    free(Y_Z0); free(Y_Z1);
    return flame;
  } /* end 0D interpolation */

  /* allocate some local memory */
  initH = (double*) calloc(flame->n1, sizeof(double));
  
  /* initialise the field from the boundary conditions given or from a given
     flamelet solution */
  /* ----- set boundary conditions for oxidizer (left, x1=0) ----- */
  /* check to see if 1D solution should be interpolated */
  double *Y0  = (double*) calloc(chem->nSpecies, sizeof(double));
  double *Y1  = (double*) calloc(chem->nSpecies, sizeof(double));
  double *Yox = (double*) calloc(chem->nSpecies, sizeof(double));
  if (flame->nDim == 1 && flow->phiMean[2] >= 0.0 && fuel2->nSp > 0) {
    for (i=0; i<ox->nSp; i++)    Y0[ox->iSp[i]]    = ox->Y[i];
    for (i=0; i<fuel2->nSp; i++) Y1[fuel2->iSp[i]] = fuel2->Y[i];
    fprintf(stderr, " Initial oxidizer species (interpolated):\n");
    for (i=0; i<chem->nSpecies; i++) {
      Yox[i] = Y0[i]*(1-flow->phiMean[2]) + Y1[i]*flow->phiMean[2];
      if (Yox[i] > 1.0e-7)
        fprintf(stderr, " %25s: %15.6e\n", chem->species[i], Yox[i]);
    }
  }
  else { /* use boundary conditions directly */
    for (i=0; i<ox->nSp; i++) Yox[ox->iSp[i]] = ox->Y[i];
  }

  /* set the species */
  for (i=0; i<chem->nSpecies; i++)
    flame->Y[i] = Yox[i];

  initH[0] = 0.0;
  if (ox->type == 'T') { /* temperature specified */
    flame->Y[0 + flame->sT] = ox->T;

    /*  compute enthalpy at boundary */
    chemComputeThermoData(chem->h, chem->cp, ox->T, chem);
    for (i=0; i<chem->nSpecies; i++)
      initH[0] += flame->Y[i]*chem->h[i];
  }
  else if (ox->type == 'H') { /* enthalpy specified */
    initH[0] = ox->H;
    
    /* compute temperature */
    ox->T = calcTfromH(chem->nSpecies, chem->nSpecSolve, initH[0],
                       &flame->Y[0], 300, chem);
    flame->Y[0 + flame->sT] = ox->T;
  }
  else if (ox->type == 'U' || ox->type == 'L')
    printf("Error: oxidiser boundary conidition specified incorrectly\n");
  
  if (driver->ENERGY) /* set the enthalpy if using diff-alg. eqs. */
    flame->H[0] = flame->Y[0 + flame->sT] = initH[0];

  /* ----- set boundary conditions for fuel 1 (right, Z1=1) ----- */
  for (i=0; i<fuel1->nSp; i++)
    flame->Y[iX1+fuel1->iSp[i]] = fuel1->Y[i];

  initH[flame->n1-1] = 0.0;
  if (fuel1->type == 'T') { /* temperature specified */
    flame->Y[iX1+flame->sT] = fuel1->T;
    
    /* compute enthalpy */
    chemComputeThermoData(chem->h, chem->cp, fuel1->T, chem);
    for (i=0; i<fuel1->nSp; i++)
      initH[flame->n1-1] += fuel1->Y[i]*chem->h[fuel1->iSp[i]];
  }
  else if (fuel1->type == 'H') { /* enthalpy specified */
    initH[flame->n1-1] = fuel1->H;
    
    /* compute temperature */
    fuel1->T = calcTfromH(chem->nSpecies,chem->nSpecSolve,initH[flame->n1-1],
                          &flame->Y[iX1],300,chem);
    flame->Y[iX1+flame->sT] = fuel1->T;
  }
  else if (fuel1->type == 'L') { /* liquid fuel specified */
    // compute the enthalpy of the vapour at the reference temperature
    chemComputeThermoData(chem->h, chem->cp, fuel1->Tref, chem);
    H_ref = 0.0;
    for (i=0; i<fuel1->nSp; i++)
      H_ref += fuel1->Y[i]*chem->h[fuel1->iSp[i]];
    
    /* compute the enthalpy of the vapour at fuel boundary */
    initH[flame->n1-1] =   ox->cpLiq*(fuel1->T - fuel1->Tref)
                         - ox->Hvap + H_ref;
    
    /* compute temperature of vapour at boundary */
    fuel1->T = calcTfromH(chem->nSpecies,chem->nSpecSolve,initH[flame->n1-1],
                          &flame->Y[iX1], 300, chem);
    
    flame->Y[iX1+flame->sT] = fuel1->T;
    
    if (MDUC_PRINT_PROC) printf("Fuel 1 specified as liquid, vapour temperature set to: %8.4f K\n",
           fuel1->T);
  }
  else if (fuel1->type == 'D') { // liquid fuel using droplet
    double hfg, cpG, cpL, Psat, Q, B, Xfs, Yfs, Ts, Hs, Textrap;
    
    /* compute the molecular weight of the fuel and oxidizer */
    for (i=0, fuel1->Mw=0.0; i<fuel1->nSp; i++)
      fuel1->Mw += fuel1->Y[i]*chem->W[fuel1->iSp[i]];
    for (i=0, ox->Mw=0.0; i<ox->nSp; i++)
      ox->Mw    += ox->Y[i]*chem->W[ox->iSp[i]];
    fuel1->P = flow->P;
    
    fprintf(stderr, "Mol. weight fuel:     %8.6f [g/mol]\n", fuel1->Mw);
    fprintf(stderr, "Mol. weight oxidizer: %8.6f [g/mol]\n", ox->Mw);
    fprintf(stderr, "Oxidizer pressure:    %15.6e [Pa]\n", fuel1->P);
    
    /* compute the theoretical surface temperature */
    Ts = brent(MDUCcalcSpaldingFunc, 183, 538, 1.0e-4, io);
    NIST_getThermo(fuel1->sat, &hfg, &cpG, &cpL, &Psat, Ts);
    Q = hfg + cpL*(Ts - fuel1->T);
    B = cpG*(ox->T - Ts)/Q;
    Xfs = Psat/fuel1->P;
    Yfs = (B + 0.0)/(1.0 + B);
    
    fprintf(stderr, "Droplet surface properties\n");
    fprintf(stderr, "    temperature:         %8.3f\n", Ts);
    fprintf(stderr, "    fuel mass fraction:  %8.3f\n", Yfs);
    fprintf(stderr, "    fuel mole fraction:  %8.3f\n", Xfs);
    fprintf(stderr, "Droplet Spalding number: %8.3f\n", B);
    
    double *YoxS   = (double*) calloc(ox->nSp, sizeof(double));
    double *YfuelS = (double*) calloc(fuel1->nSp, sizeof(double));
    for (i=0; i<ox->nSp; i++)
      YoxS[i] = linInterp(Yfs,flame->x1[0],ox->Y[i],flame->x1[flame->n1-1],0);
    for (i=0; i<fuel1->nSp; i++)
      YfuelS[i] = linInterp(Yfs,flame->x1[0],0.0,flame->x1[flame->n1-1],
                            fuel1->Y[i]);
    
    /* compute the enthalpy */
    chemComputeThermoData(chem->h, chem->cp, Ts, chem);
    for (i=0, Hs=0.0; i<ox->nSp; i++)
      Hs += YoxS[i]*chem->h[ox->iSp[i]];
    for (i=0; i<fuel1->nSp; i++)
      Hs += YfuelS[i]*chem->h[fuel1->iSp[i]];
    
    /* extrapolate the enthalpy to the boundary */
    double Hslope = (Hs-initH[0])/Yfs;
    initH[flame->n1-1] = initH[0] + Hslope;
    
    Textrap = (Ts-ox->T)/Yfs + ox->T;
    fprintf(stderr, "Extrapolated temperature: %8.3f K\n", Textrap);
    
    fprintf(stderr, "Fuel enthalpy extrapolated to: %15.6e\n",
            initH[flame->n1-1]);
    
    fuel1->T = calcTfromH(chem->nSpecies,chem->nSpecSolve,initH[flame->n1-1],
                          &flame->Y[iX1],300,chem);
    flame->Y[iX1+flame->sT] = fuel1->T;
    
    fprintf(stderr, "Fuel 1 specified as liquid droplet:\n");
    fprintf(stderr, "   setting vapour temperature to %8.4f K\n", fuel1->T);
  }
  else
    printf("Error: fuel 1 boundary condition specified incorrectly\n");

  if (driver->ENERGY) /* set the enthalpy if using diff-alg. eqs. */
    flame->H[flame->n1-1] = flame->Y[iX1+flame->sT] = initH[flame->n1-1];

  if (driver->H_DIM == 1) { /* set flow bounds to enthalpy, just in case */
    flow->phi1min = initH[0];
    flow->phi1max = initH[flame->n1-1];
  }

  /* Ignite: TODO: make this more general (multiple fuel species)*/
  double W,nuO2,nuH2O,nuCO2;
  int iFuel,iO2,iH2O,iCO2;
  double *X0  = (double*) calloc(chem->nSpecies, sizeof(double));

  if (driver->ignite == 1) {
     /* TODO: get those information from input file */
     nuH2O = 2;
     nuCO2 = 2;
     nuO2  = 3;
 
     iFuel = getSpeciesIndex(chem->nSpecies, chem->species, "C2H4");
     iO2   = getSpeciesIndex(chem->nSpecies, chem->species, "O2");
     iH2O  = getSpeciesIndex(chem->nSpecies, chem->species, "H2O");
     iCO2  = getSpeciesIndex(chem->nSpecies, chem->species, "CO2");

  }

  /* ----- linearly interpolate enthalpy and mass fractions ----- */
  for (j=1; j<flame->n1-1; j++) {
    k     = j*flame->nVars;
    k_pre = (j-1)*flame->nVars;


    /* linearly interpolate the enthalpy */
    initH[j] = linInterp(flame->x1[j], flame->x1[0], initH[0],
                         flame->x1[flame->n1-1], initH[flame->n1-1]);
      
    /* interpolate the species mass fractions */
    for (i=0; i<chem->nSpecies; i++)
      flame->Y[k+i] = linInterp(flame->x1[j], flame->x1[0], flame->Y[i],
                                flame->x1[flame->n1-1], flame->Y[iX1+i]);

    W = calcMixtureW(chem->nSpecSolve,&flame->Y[k], chem->W);

    /* Ignite: TODO: make this more general (multiple fuel species)*/
    if (driver->ignite == 1) {

       /* convert to mole fraction */
       for (i=0;i<chem->nSpecies;i++)
	  X0[i] = flame->Y[k+i] * W / chem->W[i];

       if (X0[iO2] >= nuO2 * X0[iFuel]) {// lean 
	  X0[iO2]  -= X0[iFuel]*nuO2;
	  X0[iH2O] += X0[iFuel]*nuH2O;
	  X0[iCO2] += X0[iFuel]*nuCO2;
          X0[iFuel] = 0.0;
       }
       else {// rich
          X0[iFuel] -= X0[iO2] / nuO2;
	  X0[iH2O]  += X0[iO2] / nuO2*nuH2O;
	  X0[iCO2]  += X0[iO2] / nuO2*nuCO2;
	  X0[iO2]   = 0.0;
       }
       /* renormalize mole fractions and recalculate W*/	
       double sum=0.0;
       for (i=0;i<chem->nSpecies;i++)
	  sum+=X0[i];

       W = 0.0;
       for (i=0;i<chem->nSpecies;i++) {
   	  X0[i]/=sum;
	  W += X0[i] * chem->W[i];
       }
       
       /* convert to mass fraction */       
       for (i=0;i<chem->nSpecies;i++)
	  flame->Y[k+i] = X0[i] * chem->W[i]/W;

    }

    /* get the temperature based on the enthalpy */
    if (driver->ENERGY)
      flame->H[j] = flame->Y[k+flame->sT] = initH[j];
    else
      flame->Y[k+flame->sT] = calcTfromH(chem->nSpecies, chem->nSpecSolve,
                                         initH[j], &flame->Y[k],
                                         flame->Y[k_pre+flame->sT], chem);
  } /* end loop over x1 */

  /* ----- set up the 2D field ----- */
  if (flame->nDim == 2) {
    // initialise the grid in the second dimension
    initGrid(flame->n2,0,1,flame->nPiv2,flame->gPiv2,flame->gPc2,
             &flame->x2,&flame->dx2);
    initX2dimension(driver, io, chem, flame, flow, NULL);
  }

  /* initialise the metric quantities */
  initMetric(flame);

  /* initialize soot quantities if necessary */
  if (chem->soot) {
    if (chem->sootMode <= 3) {
       for (k=0; k<flame->n2; k++) {
         for (j=0; j<flame->n1; j++) {
           initSootMoments(chem->nSootMom,
                           &flame->sootMom[(k*flame->n1+j)*flame->nSootMom],
		           1.0E-20);
        //   initSootMoments(chem->nSootMom,
        //                   &flame->sootMom[(k*flame->n1+j)*flame->nSootMom],
	//	           driver->SABSTOL);

           if (chem->sootMode == 1) {
             /* copy initial soot moments into solver variables */
             for (i=0; i<chem->nSootMom; i++)
               flame->Y[(k*flame->n1+j)*flame->nVars+flame->sSoot+i] = 
                flame->sootMom[(k*flame->n1+j)*flame->nSootMom+i];
           }
         } /* end loop over columns */
       } /* end loop over rows */
    } /* soot mode */
  } /* soot */

  /* compute the current properties */
  double wmix;
       for (k=0; k<flame->n2; k++) {
         row = k*flame->n1*flame->nVars;
         r   = k*flame->n1;
         for (j=0; j<flame->n1; j++) {
           col  = j*flame->nVars;
           wmix = calcMixtureW(chem->nSpecSolve, &flame->Y[row+col], chem->W);
      flame->Mw[r+j] = wmix;
      flame->rho[r+j] = wmix*flow->P/(flame->RGAS*flame->Y[row+col+flame->sT]);
      if (chem->readTrans)
           flame->mu[r+j]  = calcMixViscosity(chem->nSpecies,
                                           flame->Y[row+col+flame->sT],
                                                &flame->Y[row+col], chem);
         } /* end loop over columns */
       } /* end loop over rows */

  if (chem->soot) {
     if (chem->sootMode >= 3) {
       for (k=0; k<flame->n2; k++) {
          row = k*flame->n1*flame->nVars;
          r   = k*flame->n1;
         for (j=0; j<flame->n1; j++) {
            col  = j*flame->nVars;
	    double sgCoeff, oxCoeff,o2Coeff, dimProdRate;
            sootHMOM_GasPhase_Coeff(flame->Y[row+col+flame->sT], flame->mu[r+j],
            	                flame->rho[r+j], flame->Mw[r+j], chem->sootRxnConc, 
            			&sgCoeff, &oxCoeff,&o2Coeff, &dimProdRate);

	    flame->sgCoeff[r+j] 	    = sgCoeff;
	    flame->oxCoeff[r+j]  	    = oxCoeff;
	    flame->o2Coeff[r+j] 	    = o2Coeff;
	    flame->dimProdRate[r+j] 	    = dimProdRate;
            flame->sqrtT[r+j] 	    = sqrt(flame->Y[row+col+flame->sT]);
            flame->T_mu[r+j]  	    = flame->Y[row+col+flame->sT]/flame->mu[r+j];
            flame->muSqrtMw_rhoSqrtT[r+j] = flame->mu[r+j]*sqrt(flame->Mw[r+j]/
            			(1000.0*flame->Y[row+col+flame->sT]*flame->rho[r+j]*flame->rho[r+j]));

	  } 
        }
     }
  }
  if (driver->CHI_MODE == 3) {
    for (j=0; j < flame->n1; j++) 
        flame->chiA[j] = io->chiData[0].chiVal[j]; 
  }

  free(Y0); free(Y1); free(Yox); free(initH);free(X0);

  return flame;
}
/* -------------------------------------------------------------------------- */

/* this function sets parameters that are common to all Sundials solvers */
int sdialFlameletInit(double *Y0, driveStruct *driver, flameStruct *flame,
                      chemStruct *chem, sdialStruct *sdial)
{
  int i;
  sdial->mem = NULL;         /* set the memory to NULL */
  /* set the number of system eqs */
  /* TODO: for steady solver, subtract the boundary points off */
  sdial->nEq = flame->n1*flame->n2*flame->nVars;
  sdial->row = sdial->col = 0;
  sdial->totalIt = 0;

  /* DEBUG: add engine components */
  int engineEqs = 0;
  if (driver->ENGINE) {
    engineEqs += 3;
    if (driver->ENGINE > 1) /* SI */
      engineEqs += 3;
    if (driver->PRES_MODE == 3) /* implicit pressure */
      engineEqs += 1;
  }
  sdial->nEq += engineEqs;

  /* DEBUG: add reaction rates */
  if (driver->output_w) {
    sdial->nQEq = chem->nRxnTotal;
  }
  

#ifdef ACCESS_MPI

  mpiStruct *mpi = (mpiStruct*) sdial->args[7];

  if (driver->output_w) {
    sdial->nQEq *= mpi->nProc;
  }

  sdial->nEqL = mpi->nColLocal[mpi->myid]*flame->nVars;
  if (flame->nDim == 2)
    sdial->nEqL *= mpi->nRowLocal[mpi->myid];

  if (driver->output_w) {
      sdial->nQEqL = chem->nRxnTotal;
  }
  
  /* add any engine equations to the root processor? */
  if (driver->ENGINE && mpi->myid == 0) {
    /* TODO: need to add the engine equations somewhere */
    /* sdial->nEqL += engineEqs; */
  }

  sdial->Y     = N_VNew_Parallel(mpi->comm,sdial->nEqL,sdial->nEq);
  sdial->Ydata = NV_DATA_P(sdial->Y);

  if (driver->output_w) {
     sdial->Q     = N_VNew_Parallel(mpi->comm,sdial->nQEqL,sdial->nQEq);
     sdial->Qdata = NV_DATA_P(sdial->Q);
  }
  /* copy the initial solution into the solver memory */
  if (flame->nDim == 1) {
    copyArray(mpi->jCount[mpi->myid],&Y0[mpi->jStart[mpi->myid]],sdial->Ydata);
#ifdef FORCEBUG
    for (i=0; i<mpi->jStart[mpi->myid]; i++)
      Y0[i] = 0.0;
    for (i=mpi->jStart[mpi->myid]+mpi->jCount[mpi->myid]; i<sdial->nEq; i++)
      Y0[i] = 0.0;
#endif /* FORCEBUG */

    // TODO: Initialize Qdata 
    if (driver->output_w) {
         for (i=0; i<sdial->nQEqL; i++)
           sdial->Qdata[i] = chem->w_int[i]; 
    }
  }
  else if (flame->nDim == 2) {
    copyArray(mpi->kCount[mpi->myid],&Y0[mpi->kStart[mpi->myid]],sdial->Ydata);
    for (i=0; i<mpi->kStart[mpi->myid]; i++)
      Y0[i] = 0.0;
    for (i=mpi->kStart[mpi->myid]+mpi->kCount[mpi->myid]; i<sdial->nEq; i++)
      Y0[i] = 0.0;
  }
#else /* serial */
  /* set the solution vector and initialise from solution Y0
     also assign pointer to first memory location in Y
     NOTE: use compiler directive here later for parallel */
  if (flame->nDim == 1)
    sdial->Y = N_VMake_Serial(sdial->nEq,Y0);
  else if (flame->nDim == 0 && Y0)
    sdial->Y = N_VMake_Serial(sdial->nEq,Y0);
  else if (flame->nDim == 2 && sdial->ADI == 0)
    sdial->Y = N_VMake_Serial(sdial->nEq,Y0);
  else
    sdial->Y = N_VNew_Serial(sdial->nEq);

  // TODO: Initialize Qdata 
  if (driver->output_w) {
     sdial->Q = N_VNew_Serial(sdial->nQEq);
       for (i=0; i<sdial->nQEq; i++)
         sdial->Qdata[i] = chem->w_int[i]; 
  }
  /* assign pointer to access data */
  sdial->Ydata = NV_DATA_S(sdial->Y);
#endif /* ACCESS_MPI */
  return 0;
}
/* -------------------------------------------------------------------------- */

/* initialize a CVODE instance to solve the flamelet equations */
int cvodeInitFlamelet(double t0, double *Y0, sdialStruct *cvode)
{
  int flag;  /* for checking return values from cvode functions */
  int bandw; /* upper/lower bandwith */
  int i,j;

  /* dereference pointers to structs */
  driveStruct  *driver = (driveStruct*)  cvode->args[0];
  chemStruct   *chem   = (chemStruct*)   cvode->args[1];
  flameStruct  *flame  = (flameStruct*)  cvode->args[3];

  /* initialize the problem size */
  sdialFlameletInit(Y0,driver,flame,chem,cvode);

  /* create the memory object for CVODE:
     The system is stiff due to chemistry, thus choose solver appropriately
       CV_BDF:    linear multistep method (Backward Differentiation Formulas) 
       CV_NEWTON: Newton iteration (modified newton iteration for banded) */
#ifdef SUNDIALS_FROM_4_0
  cvode->mem = CVodeCreate(CV_BDF);
#else
  cvode->mem = CVodeCreate(CV_BDF,CV_NEWTON);
#endif /* SUNDIALS_FROM_4_0 */

  /* set right-hand-side function based on flamelet dimension */
  if (flame->nDim == 0)
    cvodeRHS   = cvodeRHS_0D;
  else if (flame->nDim == 1)
    cvodeRHS   = cvodeFlamelet1D;
  else if (flame->nDim == 2)
    cvodeRHS   = cvodeFlamelet2D;

#ifdef ACCESS_MPI
  if (flame->nDim == 0)
    cvodeRHS   = cvodeRHS_0D;
  else if (flame->nDim == 1)
    cvodeRHS_L = cvodeFlamelet1D_local;
  else if (flame->nDim == 2)
    cvodeRHS_L = cvodeFlamelet2D_local;
#endif /* ACCESS_MPI */

  /* set the general functions */
  SDialGetNumSteps    = CVodeGetNumSteps;
  SDialGetLastStep    = CVodeGetLastStep;
  SDialGetNumRhsEvals = CVodeGetNumRhsEvals;
  SDialSetStopTime    = CVodeSetStopTime;
  SDialSetInitStep    = CVodeSetInitStep;
  SDialGetCurrentTime = CVodeGetCurrentTime;
  SDialGetCurrentStep = CVodeGetCurrentStep;
  SDialFree           = CVodeFree;

  /* allocate internal memory object:
     initialise time:tInit
     RHS function:   cvodeRHS
     CV_SS:          scalar relative and absolute tolerances
     reltol:         relative tolerance
     abstol:         absolute tolerance
     
     NOTE: Sundials v2.4.0 is now default, support for v2.3.0 available with
           compiler flag SUNDIALS_2_3 */
#ifdef SUNDIALS_2_3
  flag = CVodeMalloc(cvode->mem,cvodeRHS,t0,cvode->Y,CV_SS,
                     driver->RELTOL,&driver->ABSTOL);

  /* set pointer to parameters (structs) needed by the RHS function */
  flag = CVodeSetFdata(cvode->mem,cvode);

  /* set the solver to return after every internal timestep */
  cvode->itask = CV_ONE_STEP_TSTOP;
  cvode->iret  = CV_TSTOP_RETURN;
#else /* Sundials v2.4.0 (default) */
  flag = CVodeInit(cvode->mem,cvodeRHS,t0,cvode->Y);
#ifdef ACCESS_MPI
  mpiStruct *mpi = (mpiStruct*) cvode->args[7];
  cvode->tmp = N_VNew_Parallel(mpi->comm,cvode->nEqL,cvode->nEq);

  /* set abstol as default */
  //for (j=0; j<mpi->jCount[mpi->myid]; j++) 
  for (j=0; j<cvode->nEqL; j++) {
	NV_Ith_P(cvode->tmp,j) = driver->ABSTOL;
  }

//  /* PAH species were found to require small tolerances if 
//     soot model is used due to nucleation source term */
//  for (j=0; j<mpi->cCount[mpi->myid]; j++)
//    for (i=0; i<chem->nSootPAH; i++)
//        NV_Ith_P(cvode->tmp,j*flame->nVars + chem->sPAH[i]) = 1.0E-30; 

  /* set abstol for T if requested */
  if (driver->TABSTOL>=0.0) {
    for (j=0; j<mpi->cCount[mpi->myid]; j++)
      NV_Ith_P(cvode->tmp,j*flame->nVars + flame->sT) = driver->TABSTOL;
  }
  if (driver->SABSTOL>=0.0) {
    for (j=0; j<mpi->cCount[mpi->myid]; j++)
      for (i=0; i<flame->nSootMom; i++) 
          NV_Ith_P(cvode->tmp,j*flame->nVars + flame->sSoot+i) = driver->SABSTOL;
  }
  else { 
    driver->SABSTOL = driver->ABSTOL;
  }

#else
  cvode->tmp = N_VNew_Serial(cvode->nEq);

  /* set abstol as default */
  for (j=0;j<flame->nVars*flame->n1;j++)
     NV_Ith_S(cvode->tmp,j) = driver->ABSTOL;

  //for (j=0; j<mpi->cCount[mpi->myid]; j++)
  //  for (i=0; i<chem->nSootPAH; i++)
  //      NV_Ith_S(cvode->tmp,j*flame->nVars + chem->sPAH[i]) = 1.0E-20; 

  /* set abstol for T if requested */
  if (driver->TABSTOL>=0.0) {
    for (j=0; j<flame->n1; j++)
       NV_Ith_S(cvode->tmp,j*flame->nVars + flame->sT) = driver->TABSTOL;
  }
  if (driver->SABSTOL>=0.0) {
    for (j=0; j<flame->n1; j++)
      for (i=0; i<flame->nSootMom; i++) 
          NV_Ith_S(cvode->tmp,j*flame->nVars + flame->sSoot+i) = driver->SABSTOL;
  }
  else { 
    driver->SABSTOL = driver->ABSTOL;
  }
#endif /* ACCESS_MPI */

  /* pass tolerance arrays to the solver */
  flag = CVodeSVtolerances(cvode->mem,driver->RELTOL,cvode->tmp);
  /*flag = CVodeSStolerances(cvode->mem,driver->RELTOL,driver->ABSTOL);*/

#ifdef SUNDIALS_FROM_3_2
  /* set constraints */
  ///* 1.0 -> constraint: y>=0 */
  N_VConst(1.0, cvode->tmp); 
  /* 0.0 -> no constraint */
  // N_VConst(0.0, cvode->tmp); 
  /* 2.0 -> constraint: y>0 */
#ifdef ACCESS_MPI


  for (j=0; j<mpi->cCount[mpi->myid]; j++)
    for (i=0; i<flame->nSootMom; i++) 
        NV_Ith_P(cvode->tmp,j*flame->nVars + flame->sSoot+i) = 2.0; 

  for (j=0; j<mpi->cCount[mpi->myid]; j++)
    NV_Ith_P(cvode->tmp,j*flame->nVars + flame->sT) = 2.0;
#else /* SERIAL */
  for (j=0; j<flame->n1; j++)
    for (i=0; i<flame->nSootMom; i++) 
        NV_Ith_S(cvode->tmp,j*flame->nVars + flame->sSoot+i) = 2.0; 

  for (j=0; j<mpi->cCount[mpi->myid]; j++)
    NV_Ith_S(cvode->tmp,j*flame->nVars + flame->sT) = 2.0;
#endif /* ACCESS_MPI */

  /* Call CVodeSetConstraints to initialize constraints */
  flag = CVodeSetConstraints(cvode->mem, cvode->tmp);
#endif /* SUNDIALS_FROM_3_2 */

  /* set pointer to parameters (structs) needed by the RHS function */
  flag = CVodeSetUserData(cvode->mem,cvode);

  /* set the solver to return after every internal timestep */
  cvode->itask = CV_ONE_STEP;
  cvode->iret  = CV_TSTOP_RETURN;
#endif /* SUNDIALS_2_3 */

  if (flag != CV_SUCCESS) {
    if (flag == CV_MEM_NULL)
      fprintf(stderr, "CVODE ERROR: CVodeMalloc requires CVodeCreate");
    printf("exiting: unable to allocate with CVodeMalloc\n");
    exit(1);
  }

  /* defined banded matrix, set number of equations, upper/lower bandwidths */
  bandw = flame->nVars;
  if (flame->nDim == 2 && !cvode->ADI)
    bandw *= flame->n1;
  if (flame->nDim == 0) {
#ifdef SUNDIALS_FROM_3_0
    cvode->A  = SUNDenseMatrix(cvode->nEq,cvode->nEq);
    cvode->LS = SUNDenseLinearSolver(cvode->Y,cvode->A);
    flag      = CVDlsSetLinearSolver(cvode->mem,cvode->LS,cvode->A);
#else /* SUNDIALS LOWER THAN 3.0 */
    flag = CVDense(cvode->mem, cvode->nEq);
#endif /* SUNDIALS_FROM_3_0 */
    cvode->row = cvode->col = 0;
  }
  else {
    if (cvode->iMethod == 0) { /* dense */
      mducWarn("Warning: using dense solver"
               "(this is inefficient in most cases)", "cvodeInitFlamelet");
#ifdef SUNDIALS_FROM_3_0
      cvode->A = SUNDenseMatrix(cvode->nEq,cvode->nEq);
      cvode->LS = SUNDenseLinearSolver(cvode->Y,cvode->A);
      flag      = CVDlsSetLinearSolver(cvode->mem,cvode->LS,cvode->A);
#else /* SUNDIALS LOWER THAN 3.0 */
      flag = CVDense(cvode->mem,cvode->nEq);
#endif /* SUNDIALS_FROM_3_0 */
      SDialGetNumJacEvals = CVDlsGetNumJacEvals;
    }
    else if (cvode->iMethod == 1) { /* banded */
#ifdef SUNDIALS_FROM_3_0
#ifdef SUNDIALS_FROM_4_0
      cvode->A  = SUNBandMatrix(cvode->nEq,flame->nVars,flame->nVars);
#else
      cvode->A  = SUNBandMatrix(cvode->nEq,flame->nVars,flame->nVars,2*flame->nVars);
#endif
      cvode->LS = SUNBandLinearSolver(cvode->Y,cvode->A);
      flag      = CVDlsSetLinearSolver(cvode->mem, cvode->LS, cvode->A);
#else /* SUNDIALS LOWER THAN 3.0 */
      flag = CVBand(cvode->mem,cvode->nEq,flame->nVars,flame->nVars);
#endif /* SUNDIALS_FROM_3_0 */
      SDialGetNumJacEvals = CVDiagGetNumRhsEvals;
    }
    else if (cvode->iMethod == 2) { /* GMRES */
#ifdef SUNDIALS_FROM_3_0
      cvode->LS = SUNSPGMR(cvode->Y,cvode->iPreCnd,0);
      flag = CVSpilsSetLinearSolver(cvode->mem, cvode->LS);
#else /* SUNDIALS LOWER THAN 3.0 */
      flag = CVSpgmr(cvode->mem,cvode->iPreCnd,0);
#endif /* SUNDIALS_FROM_3_0 */
      SDialGetNumJacEvals = CVSpilsGetNumJtimesEvals;
    }
    else if (cvode->iMethod == 3) { /* BICGSTAB */
#ifdef SUNDIALS_FROM_3_0
      cvode->LS = SUNSPBCGS(cvode->Y,cvode->iPreCnd,0);
      flag = CVSpilsSetLinearSolver(cvode->mem,cvode->LS);
#else /* SUNDIALS LOWER THAN 3.0 */
      flag = CVSpbcg(cvode->mem,cvode->iPreCnd,0);
#endif /* SUNDIALS_FROM_3_0 */
      SDialGetNumJacEvals = CVSpilsGetNumJtimesEvals;
    }
    else if (cvode->iMethod == 4) { /* TFQMR */
#ifdef SUNDIALS_FROM_3_0
      cvode->LS = SUNSPTFQMR(cvode->Y,cvode->iPreCnd,0);
      flag =  CVSpilsSetLinearSolver(cvode->mem, cvode->LS);
#else /* SUNDIALS LOWER THAN 3.0 */
      flag = CVSptfqmr(cvode->mem,cvode->iPreCnd,0);
#endif /* SUNDIALS_FROM_3_0 */
      SDialGetNumJacEvals = CVSpilsGetNumJtimesEvals;
    }
    
    if (cvode->iMethod > 1) { /* initialize precondition for Krylov methods */
#ifdef ACCESS_MPI
#ifdef SUNDIALS_FROM_3_0
      cvode->sunIndTypeTmp = flame->nVars;
      flag = CVBBDPrecInit(cvode->mem,cvode->nEqL,cvode->sunIndTypeTmp,cvode->sunIndTypeTmp,
			   cvode->sunIndTypeTmp,cvode->sunIndTypeTmp,0.0,cvodeRHS_L,NULL);
#else /* SUNDIALS LOWER THAN 3.0 */
      flag = CVBBDPrecInit(cvode->mem,cvode->nEqL,flame->nVars,flame->nVars,flame->nVars,flame->nVars,
                           0.0,cvodeRHS_L,NULL);
#endif /*SUNDIALS_FROM_3_0 */
#else /* SERIAL */
      flag = CVBandPrecInit(cvode->mem,cvode->nEq,flame->nVars,flame->nVars);
#endif /* ACCESS_MPI */
    }
//    if (driver->output_w) {
//      flag = CVodeQuadInit(cvode->mem,cvodeFlameletQ1D, cvode->Q);
//      if (driver->QABSTOL > 0 && driver->QRELTOL > 0) {	
//         flag = CVodeQuadSStolerances(cvode->mem, driver->QABSTOL, driver->QRELTOL);
#ifdef SUNDIALS_FROM_3_0
//         flag = CVodeSetQuadErrCon(cvode->mem, 0);
#else
//         flag = CVodeSetQuadErrCon(cvode->mem, TRUE);
#endif /* SUNDIALS_FROM_3_0 */
//      }
//    }
  }
  return 0;
}
/* -------------------------------------------------------------------------- */

int idaInitFlamelet(double t0, double *Y0, sdialStruct *ida)
{
  int j;
  int flag;  /* for checking return values from cvode functions */

  /* dereference pointers to structs */
  driveStruct  *driver = (driveStruct*)  ida->args[0];
  chemStruct   *chem   = (chemStruct*)   ida->args[1];
  flowStruct   *flow   = (flowStruct*)   ida->args[2];
  flameStruct  *flame  = (flameStruct*)  ida->args[3];

  /* initialize the problem size */
  sdialFlameletInit(Y0,driver,flame,chem,ida);

  /* allocate additional memory for y'(t) */
  ida->Yp     = N_VNew_Serial(ida->nEq);
  ida->YpData = NV_DATA_S(ida->Yp);

  /* create the memory object for IDA */
  ida->mem = IDACreate();

  /* set right-hand-side function based on flamelet dimension */
  if (flame->nDim == 0)
    idaRHS = idaRHS_0D;
  else if (flame->nDim == 1)
    idaRHS = idaFlamelet1D;
  else if (flame->nDim == 2)
    idaRHS = idaFlamelet2D;

  /* compute the initial derivative */
  double *T;
  if (flame->nDim == 1) {
    /* ensure temperature is up-to-date */
    if (driver->ENERGY) {
      updateEnergy1D(driver,chem,flame,flow,0,ida->Ydata,ida->YpData);
      T = flame->T;
    }
    else {
      T = (double*) calloc(flame->n1, sizeof(double));
      extractVarArray(flame->nVars,flame->n1,flame->n2,flame->sT,ida->Ydata,T);
    }

    /* get species derivatives */
    flameletY_RHS_1D(ida,driver,chem,flame,flow,T,ida->Ydata,ida->YpData);

    /* energy equation */
    if (driver->ENERGY) {
      for (j=0; j<flame->n1; j++)
        ida->YpData[j*flame->nVars + flame->sT] = flame->H[j];
    }
    else
      flameletT_RHS_1D(driver,chem,flame,flow,ida->Ydata,ida->YpData);
  }

  /* set the general functions */
  SDialGetNumSteps = IDAGetNumSteps;
  SDialGetLastStep = IDAGetLastStep;
  SDialSetStopTime = IDASetStopTime;
  SDialSetInitStep = IDASetInitStep;
  SDialFree        = IDAFree;

  /* allocate internal memory object:
     initial time: tInit

     NOTE: Sundials v2.4.0 is now default, support for v2.3.0 available with
           compiler flag SUNDIALS_2_3 */
#ifdef SUNDIALS_2_3
  flag = IDAMalloc(ida->mem,idaRHS,t0,ida->Y,ida->Yp,IDA_SS,
                   driver->RELTOL,&driver->ABSTOL);

  // set pointer to parameters (structs) needed by the RHS function
  flag = IDASetRdata(ida->mem,ida);

  ida->itask = IDA_ONE_STEP_TSTOP;
  ida->iret  = IDA_TSTOP_RETURN;
#else /* Sundials v2.4.0 (default) */
  flag = IDAInit(ida->mem,idaRHS,t0,ida->Y,ida->Yp);
  flag = IDASStolerances(ida->mem,driver->RELTOL,driver->ABSTOL);

  // set pointer to parameters (structs) needed by the RHS function
  flag = IDASetUserData(ida->mem,ida);

  ida->itask = IDA_ONE_STEP;
  ida->iret  = IDA_TSTOP_RETURN;
#endif /* SUNDIALS_2_3 */

  if (flag != IDA_SUCCESS) {
    if (flag == IDA_MEM_NULL)
      fprintf(stderr, "IDA ERROR: IDAMalloc requires IDACreate");
    printf("exiting: unable to allocate with IDAMalloc\n");
    exit(1);
  }
  
  /* defined banded matrix, set number of equations, upper/lower bandwidths */
  if (flame->nDim == 0) {
#ifdef SUNDIALS_FROM_3_0
    ida->A  = SUNDenseMatrix(ida->nEq,ida->nEq);
    ida->LS = SUNDenseLinearSolver(ida->Y,ida->A);
    flag    = IDADlsSetLinearSolver(ida->mem,ida->LS,ida->A);
#else
    flag = IDADense(ida->mem,ida->nEq);
#endif /* SUNDIALS_3_0 */
  }
  else {
#ifdef SUNDIALS_FROM_3_0
    ida->sunIndTypeTmp = flame->nVars;
#ifdef SUNDIALS_FROM_4_0
    ida->A  = SUNBandMatrix(ida->nEq, ida->sunIndTypeTmp,ida->sunIndTypeTmp);
#else
    ida->A  = SUNBandMatrix(ida->nEq, ida->sunIndTypeTmp,ida->sunIndTypeTmp,
			    ida->sunIndTypeTmp+ida->sunIndTypeTmp);
#endif /* SUNDIALS_FROM_4_0 */
    ida->LS = SUNBandLinearSolver(ida->Y,ida->A);
    flag    = IDADlsSetLinearSolver(ida->mem,ida->LS,ida->A);
#else
    flag = IDABand(ida->mem,ida->nEq,flame->nVars,flame->nVars);
#endif /* SUNDIALS_FROM_3_0 */
  }
  /* TODO: add iterative solvers */
  return 0;
}
int rk4initFlamelet(double t0, double *Y0, void *mem)
{
  /* dereference pointers to structs */
  mducStruct   *mduc   = (mducStruct*) mem;
  driveStruct  *driver = mduc->driver;
  chemStruct   *chem   = mduc->chem;
  flowStruct   *flow   = mduc->flow;
  flameStruct  *flame  = mduc->flame;
  sdialStruct  *sdial  = mduc->sdial;
  rkStruct     *rk     = mduc->rk;

  /* initialize the problem size */
  sdialFlameletInit(Y0,driver,flame,chem,sdial); 
  rk4init(mem);
  return 0;
}
/* -------------------------------------------------------------------------- */
