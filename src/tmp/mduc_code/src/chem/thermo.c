#include <stdlib.h>
#include <math.h>
#include "thermo.h"

void CEAcalcSpeciesThermo(int n,molStruct *sp,double *h,double *cp,double *S,double T)
{
  int    i,iT;
  double *a,*b,*s,R = 8314.34,R_W,logT,Tinv;

  /* constant across species */
  Tinv = 1.0/T;
  logT = log(T);

  for (i=0; i<n; i++) { /* loop over species */
    /* determine the correct temperature range */
    iT = -1;
    if (T >= sp[i]->Trange[0]) iT++;
    if (T >= sp[i]->Trange[1]) iT++;
    if (T >= sp[i]->Trange[2]) iT++;
    if (sp->nTrange == 3)
      if (T > sp[i]->Trange[3]) iT++;
    if (iT < 0 || iT > 3) {
      fprintf(stderr,"Error: temperature (%f) is out of range for species %s\n",
              T, chem->specs[i].name);
      exit(2);
    }
    
    /* assign the appropriate coefficients */
    a = sp[i]->aH[iT];
    b = sp[i]->aCp[iT];
    s = sp[i]->aS[iT];
    
    R_W  = R/sp->W;

    /* ----- compute the enthalpy ----- */
    h[i]  = R_W*(Tinv*a[0]+a[1]*logT+
                 a[2]+T*(a[3]+T*(a[4]+T*(a[5]+T*(a[6]+T*a[7])))));
    
    /* ----- compute the specific heat ----- */
    cp[i] = R_W*(Tinv*(Tinv*b[0]+b[1])+b[2]+T*(b[3]+T*(b[4]+T*(b[5]+T*b[6]))));
    
    /* ----- compute the entropy ----- */
    S[i]  = R_W*(Tinv*(Tinv*s[0]+s[1])+s[2]*log(T)+
                 s[3]+T*(s[4]+T*(s[5]+T*(s[6]+T*s[7]))));
  } /* end loop over species */
}
/* -------------------------------------------------------------------------- */
