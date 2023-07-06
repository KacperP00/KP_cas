#include <stdlib.h>
#include <stdio.h>
#include "destroy.h"
#include "parser.h"

void destroyFlame(flameStruct *flame)
{
  int i;

  if (flame->x1)           free(flame->x1);
  if (flame->dx1)          free(flame->dx1);
  if (flame->Y)            free(flame->Y);
  if (flame->Yp)           free(flame->Yp);
  if (flame->Y_m1)         free(flame->Y_m1);
  if (flame->Y_m2)         free(flame->Y_m2);
  if (flame->Yp_m1)        free(flame->Yp_m1);
  if (flame->Y_pre)        free(flame->Y_pre);
  if (flame->sumYi)        free(flame->sumYi);
  if (flame->v1)           free(flame->v1);
  if (flame->xi)           free(flame->xi);
  if (flame->D)            free(flame->D);
  if (flame->rhoSigma)     free(flame->rhoSigma);
  if (flame->chiA)         free(flame->chiA);
  if (flame->rhoDotSoot)   free(flame->rhoDotSoot);
  if (flame->rhoDotSpray)  free(flame->rhoDotSpray);
  if (flame->chiB)         free(flame->chiB);
  if (flame->chiC)         free(flame->chiC);
  if (flame->H)            free(flame->H);
  if (flame->HRrate)       free(flame->HRrate);
  if (flame->rho)          free(flame->rho);
  if (flame->gamma)        free(flame->gamma);
  if (flame->cpmix)        free(flame->cpmix);
  if (flame->Mw)           free(flame->Mw);
  if (flame->R)            free(flame->R);
  if (flame->mu)           free(flame->mu);
  if (flame->srcQ)         free(flame->srcQ);
  if (flame->sumCp)        free(flame->sumCp);
  if (flame->srcP)         free(flame->srcP);
  if (flame->srcH)         free(flame->srcH);
  if (flame->chemSrcY)     free(flame->chemSrcY);
  if (flame->chemSrcY_tmp) free(flame->chemSrcY_tmp);
  if (flame->visc_src)     free(flame->visc_src);
  if (flame->radiation)    free(flame->radiation);
  if (flame->gPDF)         free(flame->gPDF);
  if (flame->chem_src)     free(flame->chem_src);
  if (flame->conv_src)     free(flame->conv_src);
  if (flame->diff_src)     free(flame->diff_src);
  if (flame->chem_temp)    free(flame->chem_temp);
  if (flame->conv_temp)    free(flame->conv_temp);
  if (flame->diff_temp)    free(flame->diff_temp);
  if (flame->convCoeff)    free(flame->convCoeff);
  if (flame->sootMom)      free(flame->sootMom);
  if (flame->sootMomSrc)   free(flame->sootMomSrc);
  if (flame->sootMomTrans) free(flame->sootMomTrans);
  if (flame->sootMomTransDiff)  free(flame->sootMomTransDiff);
  if (flame->sootMomTransConv)  free(flame->sootMomTransConv);

  for (i=0; i<flame->n1; i++)
    {
      if (flame->grad.x) free(flame->grad.x[i]);
      if (flame->grad.xx) free(flame->grad.xx[i]);
    }
  if (flame->grad.x)  free(flame->grad.x);
  if (flame->grad.xx) free(flame->grad.xx);
  
  if (flame->nDim == 2)
    {
      if (flame->x2)     free(flame->x2);
      if (flame->dx2)    free(flame->dx2);
      if (flame->sumCp2) free(flame->sumCp2);
      if (flame->v2)     free(flame->v2);
      destroyGrad(&flame->dY);
      destroyGrad(&flame->dCp);

      for (i=0; i<flame->n2; i++)
        {
          if (flame->grad.y) free(flame->grad.y[i]);
          if (flame->grad.yy) free(flame->grad.yy[i]);
        }
      if (flame->grad.y)  free(flame->grad.y);
      if (flame->grad.yy) free(flame->grad.yy);
    }
}
/* -------------------------------------------------------------------------- */

void destroyGrad(gradStruct *grad)
{
  if (grad->x)  free(grad->x);
  if (grad->y)  free(grad->y);
  if (grad->xx) free(grad->xx);
  if (grad->yy) free(grad->yy);
  if (grad->xy) free(grad->xy);
}
/* -------------------------------------------------------------------------- */

void destroyFlow(flowStruct *flow)
{
  int i;
  if (flow->chi1)        free(flow->chi1);
  if (flow->chi2)        free(flow->chi2);
  if (flow->chi12)       free(flow->chi12);
  if (flow->chi1_m1)     free(flow->chi1_m1);
  if (flow->chi2_m1)     free(flow->chi2_m1);
  if (flow->chi12_m1)    free(flow->chi12_m1);
  if (flow->pdf)         free(flow->pdf);
  if (flow->avgSdot)     free(flow->avgSdot);
  if (flow->avgSdot_m1)  free(flow->avgSdot_m1);
  if (flow->avgSdot_M)   free(flow->avgSdot_M);
  if (flow->strainCondZ) free(flow->strainCondZ);
  if (flow->avgH)        free(flow->avgH);

  if (flow->phiStar)    free(flow->phiStar);
  if (flow->phiRaw) {
    for (i=0; i<3; i++) if (flow->phiRaw[i]) free(flow->phiRaw[i]);
    free(flow->phiRaw);
  }
  if (flow->phiHat) {
    for (i=0; i<3; i++) if (flow->phiHat[i]) free(flow->phiHat[i]);
    free(flow->phiHat);
  }
  if (flow->Dt_phiRaw) {
    for (i=0; i<3; i++) if (flow->Dt_phiRaw[i]) free(flow->Dt_phiRaw[i]);
    free(flow->Dt_phiRaw);
  }
  if (flow->Dt_phiHat) {
    for (i=0; i<3; i++) if (flow->Dt_phiHat[i]) free(flow->Dt_phiHat[i]);
    free(flow->Dt_phiHat);
  }
  if (flow->DtHat_phi) {
    for (i=0; i<3; i++) if (flow->DtHat_phi[i]) free(flow->DtHat_phi[i]);
    free(flow->DtHat_phi);
  }
  if (flow->fhat.w) free(flow->fhat.w);
}
/* -------------------------------------------------------------------------- */

void destroyIO(ioStruct *io)
{
  int i;
  
  if (io->nDim) free(io->nDim);
  if (io->n1)   free(io->n1);
  if (io->n2)   free(io->n2);
  
  for (i=0; i<io->nInputs; ++i) {
    if (io->inputSpecies[i]) free(io->inputSpecies[i]);
    if (io->initY[i]) free(io->initY[i]);
  }
  if (io->inputSpecies) free(io->inputSpecies);
  if (io->inputIndex) free(io->inputIndex);
  if (io->initY) free(io->initY);
  if (io->initT) free(io->initT);
  if (io->gPDF) free(io->gPDF);
  
  for (i=0; i<io->nOutputs; ++i) 
    if (io->outputSpecies[i]) free(io->outputSpecies[i]);
  if (io->outputSpecies) free(io->outputSpecies);
  if (io->outputIndex)   free(io->outputIndex);
  if (io->outputVar) {
    for (i=0; i<io->nOutputVar; i++) free(io->outputVar[i]);
    free(io->outputVar);
  }
  if (io->chemFuel) {
    for (i=0; i<io->nChemFuel; i++) free(io->chemFuel[i]);
    free(io->chemFuel);
  }
  if (io->chemFuelY) free(io->chemFuelY);
  
  for (i=0; i<io->nMonitorY; ++i)
    if (io->monitorSpecies[i]) free(io->monitorSpecies[i]);
  if (io->monitorSpecies) free(io->monitorSpecies);
  if (io->monitorIndex)   free(io->monitorIndex);
  
  for (i=0; i<io->nOx; ++i)    if (io->oxName[i])    free(io->oxName[i]);
  for (i=0; i<io->nFuel1; ++i) if (io->fuelName1[i]) free(io->fuelName1[i]);
  if (io->oxName)     free(io->oxName);
  if (io->oxIndex)    free(io->oxIndex);
  if (io->oxY)        free(io->oxY);
  if (io->fuelName1)  free(io->fuelName1);
  if (io->fuelIndex1) free(io->fuelIndex1);
  if (io->fuelY1)     free(io->fuelY1);
  
  if (io->nFuel2 > 0) {
    for (i=0; i<io->nFuel2; ++i) if (io->fuelName2[i]) free(io->fuelName2[i]);
    if (io->fuelName2)  free(io->fuelName2);
    if (io->fuelIndex2) free(io->fuelIndex2);
    if (io->fuelY2)     free(io->fuelY2);
  }
  
  if (io->nMisc > 0) {
    for (i=0; i<io->nMisc; ++i) {
      if (io->miscIn[i])    free(io->miscIn[i]);
      if (io->miscNames[i]) free(io->miscNames[i]);
    }
    if (io->miscIn)    free(io->miscIn);
    if (io->miscNames) free(io->miscNames);
  }
  
  if (io->gPivot1) free(io->gPivot1);
  if (io->gPivot2) free(io->gPivot2);
  if (io->gPerc1)  free(io->gPerc1);
  if (io->gPerc2)  free(io->gPerc2);

  if (io->parser) destroyParser(io->parser);
  free(io->parser);

  if (io->monitor) fclose(io->monitor);
}
/* -------------------------------------------------------------------------- */

void destroySDial(int solver, sdialStruct *sdial)
{
  // free any solution vectors
  if (sdial->Y != NULL) {
#ifdef ACCESS_MPI
    N_VDestroy_Parallel(sdial->Y);
#else
    N_VDestroy_Serial(sdial->Y);
#endif
    sdial->Y     = NULL;
    sdial->Ydata = NULL;
  }
  if (sdial->Q != NULL) {
#ifdef ACCESS_MPI
    N_VDestroy_Parallel(sdial->Q);
#else
    N_VDestroy_Serial(sdial->Q);
#endif
    sdial->Q     = NULL;
    sdial->Qdata = NULL;
  }
  if (sdial->Yp != NULL) {
#ifdef ACCESS_MPI
    N_VDestroy_Parallel(sdial->Yp);
#else
    N_VDestroy_Serial(sdial->Yp);
#endif
    sdial->Yp     = NULL;
    sdial->YpData = NULL;
  }
  if (sdial->Yscale != NULL) {
#ifdef ACCESS_MPI
    N_VDestroy_Parallel(sdial->Yscale);
#else
    N_VDestroy_Serial(sdial->Yscale);
#endif
    sdial->Yscale     = NULL;
    sdial->YscaleData = NULL;
  }
  if (sdial->Fscale != NULL) {
#ifdef ACCESS_MPI
    N_VDestroy_Parallel(sdial->Fscale);
#else
    N_VDestroy_Serial(sdial->Fscale);
#endif
    sdial->Fscale     = NULL;
    sdial->FscaleData = NULL;
  }
  if (sdial->tmp != NULL) {
#ifdef ACCESS_MPI
    N_VDestroy_Parallel(sdial->tmp);
#else
    N_VDestroy_Serial(sdial->tmp);
#endif
    sdial->tmp = NULL;
  }

  // free the integrator memory
  if (sdial->mem != NULL)
    SDialFree(&sdial->mem);
#ifdef SUNDIALS_FROM_3_0
  if (sdial->A != NULL) 
    SUNMatDestroy(sdial->A);
  if (sdial->LS != NULL)
    SUNLinSolFree(sdial->LS);
#endif /* SUNDIALS_FROM_3_0 */
}
/* -------------------------------------------------------------------------- */
