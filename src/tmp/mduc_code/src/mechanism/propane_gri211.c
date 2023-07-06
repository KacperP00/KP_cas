#include <math.h>
#include <string.h>
#include "propane_gri211.h"

static double GetLindRateCoeff( double temp, double pressure
			, double k0, double kInf, double fc, double conc );

double MAX_C(double X1, double X2);

void ComputeProductionRates( double *cdot, double *w, double *k
			, double *c, double *M, double temp, double pressure )
{
/*
	This function computes rates of production cdot in [kmole/(m^3s)].
	The parameters w ( reaction rate ), k ( rate coefficient )
	and M ( third body concentrations ) are just work space for this function.
	c contains the concentrations of non steady state species in [kmole/m^3]
	and is workspace for the steady state concentrations, which are computed
	in this function.
	temp is the temperature in [K] and pressure is the pressure in [Pa].
	Called functions are 'GetLindRateCoeff', 'ComputeSteadyStates',
	'CatchZero' and the functions that evaluate the broadening factors
	of the Troe formulation of pressure dependant rate coefficients 'Fc*'
*/

	int	nSpec = 56;
	int	nSpecIn = 56;
	double	kTroe0, kTroeInf, fcTroe;
	double	lgt = log( temp );
	double	rt = 8314.34 * temp;

	M[mM1] = c[sN2] + 0.83 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2.4 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.75 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 3.6 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 15.4 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM2] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM3] = c[sN2] + 0.5 * c[sAR] + c[sO] + 6 * c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 3.5 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM4] =  + c[sO] + c[sH] + c[sOH] + c[sH2]
		 + c[sHO2] + c[sH2O2] + c[sCH] + 0.75 * c[sCO]
		 + c[sCH2] + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3]
		 + c[sCH2O] + c[sCH4] + 1.5 * c[sCO2] + c[sCH2OH]
		 + c[sCH3O] + c[sCH3OH] + c[sC2H] + c[sC2H2]
		 + c[sHCCO] + c[sC2H3] + c[sCH2CO] + c[sC2H4]
		 + c[sC2H5] + 1.5 * c[sC2H6] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM5] = c[sN2] + 0.63 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + c[sHO2] + c[sH2O2]
		 + c[sCH] + c[sCO] + c[sCH2] + c[sHCO]
		 + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O] + 2 * c[sCH4]
		 + c[sCH2OH] + c[sCH3O] + c[sCH3OH] + c[sC2H]
		 + c[sC2H2] + c[sHCCO] + c[sC2H3] + c[sCH2CO]
		 + c[sC2H4] + c[sC2H5] + 3 * c[sC2H6] + c[sC]
		 + c[sHCCOH] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN]
		 + c[sHCNO] + c[sNH3] + c[sC3H3] + c[sC3H4]
		 + c[sC3H5] + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7]
		 + c[sC3H8];

	M[mM6] = c[sN2] + 0.38 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 0.73 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 3.65 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM7] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM8] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM9] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM10] = c[sN2] + c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM11] = c[sN2] + c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM12] = c[sN2] + c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM13] = c[sN2] + c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM14] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM15] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM16] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM17] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM18] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM19] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM20] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM21] = c[sN2] + c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM22] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM23] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM24] = c[sN2] + c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM25] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM26] = c[sN2] + c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + c[sC] + c[sHCCOH] + c[sN]
		 + c[sNO] + c[sN2O] + c[sNO2] + c[sNH]
		 + c[sHNO] + c[sNH2] + c[sNNH] + c[sCN]
		 + c[sNCO] + c[sHCN] + c[sHOCN] + c[sHNCO]
		 + c[sH2CN] + c[sHCNN] + c[sHCNO] + c[sNH3]
		 + c[sC3H3] + c[sC3H4] + c[sC3H5] + c[sC3H6]
		 + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM27] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM28] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM29] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM30] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM31] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM32] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM33] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM34] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM35] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];

	M[mM36] = c[sN2] + 0.7 * c[sAR] + c[sO] + c[sO2]
		 + c[sH] + c[sOH] + 2 * c[sH2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + 1.5 * c[sCO] + c[sCH2]
		 + c[sHCO] + c[sCH2GSGXCH2] + c[sCH3] + c[sCH2O]
		 + 2 * c[sCH4] + 2 * c[sCO2] + c[sCH2OH] + c[sCH3O]
		 + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO]
		 + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5]
		 + 3 * c[sC2H6] + 6 * c[sH2O] + c[sC] + c[sHCCOH]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO]
		 + c[sNH3] + c[sC3H3] + c[sC3H4] + c[sC3H5]
		 + c[sC3H6] + c[sNXC3H7] + c[sIXC3H7] + c[sC3H8];



	k[r1f] = 1.2000000000E+11 * exp( -1 * lgt );
	k[r1b] = 8.4127616577E+15 * exp( -1.12163 * lgt - 497822466 / rt );
	k[r2f] = 5.0000000000E+11 * exp( -1 * lgt );
	k[r2b] = 6.9505130758E+13 * exp( -0.684828 * lgt - 425203810.5 / rt );
	k[r3f] = 5.0000000000E+01 * exp( 2.67 * lgt - 26317000 / rt );
	k[r3b] = 2.9737012188E+01 * exp( 2.63494 * lgt - 18429322.73 / rt );
	k[r4f] = 2.0000000000E+10;
	k[r4b] = 1.9300151207E+09 * exp( 0.32701 * lgt - 220535995.2 / rt );
	k[r5f] = 9.6300000000E+03 * exp( 2 * lgt - 16736000 / rt );
	k[r5b] = 1.9422606141E+01 * exp( 2.5719 * lgt - 75687466.21 / rt );
	k[r6f] = 5.7000000000E+10;
	k[r6b] = 4.0039790897E+12 * exp( -0.0551873 * lgt - 739371737.2 / rt );
	k[r7f] = 8.0000000000E+10;
	k[r7b] = 7.8223107159E+12 * exp( -0.345666 * lgt - 383372794.3 / rt );
	k[r8f] = 1.5000000000E+10;
	k[r8b] = 4.7190290127E+10 * exp( 0.174995 * lgt - 787993698.2 / rt );
	k[r9f] = 1.5000000000E+10;
	k[r9b] = 1.0693498414E+12 * exp( -0.407952 * lgt - 420675782.2 / rt );
	k[r10f] = 8.4300000000E+10;
	k[r10b] = 1.1196932142E+13 * exp( -0.275761 * lgt - 288886347 / rt );
	k[r11f] = 1.0200000000E+06 * exp( 1.5 * lgt - 35982000 / rt );
	k[r11b] = 3.4821276777E+02 * exp( 1.9765 * lgt - 19030609.55 / rt );
	k[r12f] = 6.0200000000E+08 * exp( -12552000 / rt );
	k[r12b] = 4.6182917899E+17 * exp( -0.963346 * lgt - 547255890 / rt );
	k[r13f] = 3.0000000000E+10;
	k[r13b] = 7.8737463365E+08 * exp( 0.547884 * lgt - 359430238.7 / rt );
	k[r14f] = 3.0000000000E+10;
	k[r14b] = 4.3452967068E+15 * exp( -0.730634 * lgt - 468930318.2 / rt );
	k[r15f] = 3.9000000000E+10 * exp( -14811000 / rt );
	k[r15b] = 3.8627520933E+07 * exp( 0.41439 * lgt - 70409403.98 / rt );
	k[r16f] = 1.0000000000E+10;
	k[r16b] = 1.3918415136E+08 * exp( 0.696539 * lgt - 302440224 / rt );
	k[r17f] = 1.0000000000E+10;
	k[r17b] = 1.5736169844E+07 * exp( 0.637464 * lgt - 330967678.1 / rt );
	k[r18f] = 3.8800000000E+02 * exp( 2.5 * lgt - 12970000 / rt );
	k[r18b] = 1.7260888544E+00 * exp( 2.67595 * lgt - 32461855.27 / rt );
	k[r19f] = 1.3000000000E+02 * exp( 2.5 * lgt - 20920000 / rt );
	k[r19b] = 5.1152342447E+00 * exp( 2.73503 * lgt - 11884401.26 / rt );
	k[r20f] = 5.0000000000E+10;
	k[r20b] = 2.3360600807E+07 * exp( 1.01264 * lgt - 325621371.1 / rt );
	k[r21f] = 1.0200000000E+04 * exp( 2 * lgt - 7950000 / rt );
	k[r21b] = 3.1458063947E+04 * exp( 1.81798 * lgt - 89764345.93 / rt );
	k[r22f] = 4.6000000000E+16 * exp( -1.41 * lgt - 121127000 / rt );
	k[r22b] = 9.3604172115E+12 * exp( -0.848524 * lgt + 10152055.94 / rt );
	k[r23f] = 1.0200000000E+04 * exp( 2 * lgt - 7950000 / rt );
	k[r23b] = 2.6543809692E-02 * exp( 3.31672 * lgt - 198861019.4 / rt );
	k[r24f] = 3.0000000000E+10;
	k[r24b] = 1.2732201502E+14 * exp( -0.600835 * lgt - 379074877 / rt );
	k[r25f] = 1.9200000000E+04 * exp( 1.83 * lgt - 920000 / rt );
	k[r25b] = 2.2200607391E-01 * exp( 2.70863 * lgt - 108375950.5 / rt );
	k[r26f] = 1.3200000000E+11;
	k[r26b] = 1.1465114667E+08 * exp( 0.900726 * lgt - 325723724.7 / rt );
	k[r27f] = 8.9800000000E+04 * exp( 1.92 * lgt - 23807000 / rt );
	k[r27b] = 1.0047450020E+01 * exp( 2.44412 * lgt - 27335892.79 / rt );
	k[r28f] = 1.0000000000E+11;
	k[r28b] = 1.5577208903E+03 * exp( 1.38578 * lgt - 426695895.8 / rt );
	k[r29f] = 1.0000000000E+10 * exp( -33472000 / rt );
	k[r29b] = 4.5588427293E+06 * exp( 0.562107 * lgt - 14614041.66 / rt );
	k[r30f] = 1.7500000000E+09 * exp( -5648000 / rt );
	k[r30b] = 3.7150267058E+06 * exp( 0.782321 * lgt - 205386794.6 / rt );
	k[r31f] = 2.5000000000E+09 * exp( -199995000 / rt );
	k[r31b] = 2.7356941842E+13 * exp( -0.84172 * lgt - 236876424.1 / rt );
	k[r32f] = 1.0000000000E+11 * exp( -167360000 / rt );
	k[r32b] = 1.0263642436E+09 * exp( 0.0873801 * lgt - 2422408.767 / rt );
	k[r33f] = 2.8000000000E+12 * exp( -0.86 * lgt );
	k[r33b] = 4.0334267652E+15 * exp( -0.871839 * lgt - 204667815.3 / rt );
	k[r34f] = 3.0000000000E+14 * exp( -1.72 * lgt );
	k[r34b] = 4.3215286770E+17 * exp( -1.73184 * lgt - 204667815.3 / rt );
	k[r35f] = 9.3800000000E+12 * exp( -0.76 * lgt );
	k[r35b] = 1.3511979663E+16 * exp( -0.771839 * lgt - 204667815.3 / rt );
	k[r36f] = 3.7500000000E+14 * exp( -1.72 * lgt );
	k[r36b] = 5.4019108462E+17 * exp( -1.73184 * lgt - 204667815.3 / rt );
	k[r37f] = 7.0000000000E+11 * exp( -0.8 * lgt );
	k[r37b] = 1.0083566913E+15 * exp( -0.811839 * lgt - 204667815.3 / rt );
	k[r38f] = 8.3000000000E+10 * exp( -60304000 / rt );
	k[r38b] = 1.6457642104E+08 * exp( 0.436798 * lgt + 12314655.44 / rt );
	k[r39f] = 1.0000000000E+12 * exp( -1 * lgt );
	k[r39b] = 2.3373273118E+14 * exp( -0.649765 * lgt - 433091487.8 / rt );
	k[r40f] = 9.0000000000E+10 * exp( -0.6 * lgt );
	k[r40b] = 2.1035945806E+13 * exp( -0.249765 * lgt - 433091487.8 / rt );
	k[r41f] = 6.0000000000E+13 * exp( -1.25 * lgt );
	k[r41b] = 1.4023963871E+16 * exp( -0.899765 * lgt - 433091487.8 / rt );
	k[r42f] = 5.5000000000E+14 * exp( -2 * lgt );
	k[r42b] = 1.2855300215E+17 * exp( -1.64977 * lgt - 433091487.8 / rt );
	k[r43f] = 2.2000000000E+16 * exp( -2 * lgt );
	k[r43b] = 5.7024497877E+19 * exp( -1.7546 * lgt - 497934272.2 / rt );
	k[r44f] = 3.9700000000E+09 * exp( -2807000 / rt );
	k[r44b] = 1.4164546439E+07 * exp( 0.694036 * lgt - 223454801.5 / rt );
	k[r45f] = 2.8000000000E+10 * exp( -4469000 / rt );
	k[r45b] = 4.5431954492E+09 * exp( 0.362073 * lgt - 232892672.5 / rt );
	k[r46f] = 1.3400000000E+11 * exp( -2657000 / rt );
	k[r46b] = 2.5640414139E+07 * exp( 0.763808 * lgt - 150574339.8 / rt );
	k[r47f] = 1.2100000000E+04 * exp( 2 * lgt - 21757000 / rt );
	k[r47b] = 4.1033565943E+01 * exp( 2.60696 * lgt - 88596143.48 / rt );
	k[r48f] = 1.0000000000E+10 * exp( -15062000 / rt );
	k[r48b] = 7.1960369081E+04 * exp( 1.26593 * lgt - 294661267.7 / rt );
	k[r49f] = 1.1000000000E+11;
	k[r49b] = 5.6099788271E+10 * exp( 0.23442 * lgt - 97396261.59 / rt );
	kTroe0 = 3.2000000000E+21 * exp( -3.14 * lgt - 5146000 / rt );
	kTroeInf = 2.5000000000E+13 * exp( -0.8 * lgt );
	fcTroe = 0.32 * exp( -temp / 78 ) + 0.68 * exp( -temp / 1995 ) + 1 * exp( -5590 / temp );
	k[r50f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM7] );
	kTroe0 = 3.3062692160E+26 * exp( -3.30912 * lgt - 469237853.8 / rt );
	kTroeInf = 2.5830228250E+18 * exp( -0.969123 * lgt - 464091853.8 / rt );
	k[r50b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM7] );
	k[r51f] = 3.0000000000E+10;
	k[r51b] = 1.3435867056E+09 * exp( 0.230182 * lgt - 48621961.01 / rt );
	kTroe0 = 2.4770000000E+27 * exp( -4.76 * lgt - 10209000 / rt );
	kTroeInf = 1.2700000000E+13 * exp( -0.63 * lgt - 1602000 / rt );
	fcTroe = 0.217 * exp( -temp / 74 ) + 0.783 * exp( -temp / 2941 ) + 1 * exp( -6964 / temp );
	k[r52f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM8] );
	kTroe0 = 1.0086217929E+33 * exp( -4.92133 * lgt - 452364201 / rt );
	kTroeInf = 5.1713753613E+18 * exp( -0.791325 * lgt - 443757201 / rt );
	k[r52b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM8] );
	k[r53f] = 6.6000000000E+05 * exp( 1.62 * lgt - 45355000 / rt );
	k[r53b] = 3.7884462370E+02 * exp( 2.13156 * lgt - 36291286.82 / rt );
	kTroe0 = 1.3500000000E+18 * exp( -2.57 * lgt - 5962000 / rt );
	kTroeInf = 1.0900000000E+09 * exp( 0.48 * lgt + 1088000 / rt );
	fcTroe = 0.2176 * exp( -temp / 271 ) + 0.7824 * exp( -temp / 2755 ) + 1 * exp( -6570 / temp );
	k[r54f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM9] );
	kTroe0 = 1.8947346586E+23 * exp( -2.66922 * lgt - 375567406.5 / rt );
	kTroeInf = 1.5298227985E+14 * exp( 0.380781 * lgt - 368517406.5 / rt );
	k[r54b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM9] );
	k[r55f] = 7.3400000000E+10;
	k[r55b] = 3.2391338749E+09 * exp( 0.582947 * lgt - 367317915.9 / rt );
	kTroe0 = 1.2700000000E+26 * exp( -4.82 * lgt - 27322000 / rt );
	kTroeInf = 5.4000000000E+08 * exp( 0.454 * lgt - 15062000 / rt );
	fcTroe = 0.2813 * exp( -temp / 103 ) + 0.7187 * exp( -temp / 1291 ) + 1 * exp( -4160 / temp );
	k[r56f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM10] );
	kTroe0 = 1.2684133243E+30 * exp( -5.20137 * lgt - 150085586.5 / rt );
	kTroeInf = 5.3932535051E+12 * exp( 0.0726333 * lgt - 137825586.5 / rt );
	k[r56b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM10] );
	kTroe0 = 2.2000000000E+24 * exp( -4.8 * lgt - 23263000 / rt );
	kTroeInf = 5.4000000000E+08 * exp( 0.454 * lgt - 10878000 / rt );
	fcTroe = 0.242 * exp( -temp / 94 ) + 0.758 * exp( -temp / 1555 ) + 1 * exp( -4200 / temp );
	k[r57f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM11] );
	kTroe0 = 1.9434371792E+29 * exp( -5.12229 * lgt - 117499132.5 / rt );
	kTroeInf = 4.7702548944E+13 * exp( 0.131708 * lgt - 105114132.5 / rt );
	k[r57b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM11] );
	k[r58f] = 2.3000000000E+07 * exp( 1.05 * lgt - 13703000 / rt );
	k[r58b] = 3.8302995462E+04 * exp( 1.49945 * lgt - 77189081.25 / rt );
	kTroe0 = 3.0000000000E+25 * exp( -4.8 * lgt - 13807000 / rt );
	kTroeInf = 1.8000000000E+10;
	fcTroe = 0.2321 * exp( -temp / 338 ) + 0.7679 * exp( -temp / 1812 ) + 1 * exp( -5081 / temp );
	k[r59f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM12] );
	kTroe0 = 9.3742534746E+29 * exp( -4.66078 * lgt - 419518955.3 / rt );
	kTroeInf = 5.6245520847E+14 * exp( 0.139219 * lgt - 405711955.3 / rt );
	k[r59b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM12] );
	k[r60f] = 2.0000000000E+10;
	k[r60b] = 4.6805022133E+08 * exp( 0.731602 * lgt - 310327901.3 / rt );
	k[r61f] = 1.2000000000E+10;
	k[r61b] = 1.2574755812E+06 * exp( 0.972299 * lgt - 13553877.06 / rt );
	k[r62f] = 6.0000000000E+09;
	k[r62b] = 2.1633990968E+04 * exp( 1.44911 * lgt - 10093307.5 / rt );
	kTroe0 = 8.6000000000E+22 * exp( -4 * lgt - 12657000 / rt );
	kTroeInf = 5.0000000000E+10;
	fcTroe = 0.1098 * exp( -temp / 144 ) + 0.8902 * exp( -temp / 2838 ) + 1 * exp( -45570 / temp );
	k[r63f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM13] );
	kTroe0 = 3.0382474183E+26 * exp( -3.91986 * lgt - 446896409.3 / rt );
	kTroeInf = 1.7664229176E+14 * exp( 0.0801444 * lgt - 434239409.3 / rt );
	k[r63b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM13] );
	k[r64f] = 3.4000000000E+03 * exp( 1.6 * lgt );
	k[r64b] = 3.8440423674E+02 * exp( 1.54093 * lgt - 28527454.01 / rt );
	k[r65f] = 2.0000000000E+10;
	k[r65b] = 5.2917790613E+07 * exp( 0.672527 * lgt - 338855355.3 / rt );
	k[r66f] = 3.2000000000E+10;
	k[r66b] = 3.7912073804E+05 * exp( 0.913225 * lgt - 42081331.07 / rt );
	k[r67f] = 1.6000000000E+10;
	k[r67b] = 6.5225080672E+03 * exp( 1.39003 * lgt - 38620761.51 / rt );
	k[r68f] = 1.7000000000E+04 * exp( 2.1 * lgt - 20376000 / rt );
	k[r68b] = 1.2716073149E+02 * exp( 2.31102 * lgt - 47755532.54 / rt );
	k[r69f] = 4.2000000000E+03 * exp( 2.1 * lgt - 20376000 / rt );
	k[r69b] = 2.7787158476E+02 * exp( 2.37009 * lgt - 19228078.53 / rt );
	kTroe0 = 3.7500000000E+27 * exp( -4.8 * lgt - 7950000 / rt );
	kTroeInf = 1.0000000000E+14 * exp( -1 * lgt );
	fcTroe = 0.3536 * exp( -temp / 132 ) + 0.6464 * exp( -temp / 1315 ) + 1 * exp( -5566 / temp );
	k[r70f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM14] );
	kTroe0 = 2.5617736443E+33 * exp( -5.0463 * lgt - 564432866.5 / rt );
	kTroeInf = 6.8313963848E+19 * exp( -1.2463 * lgt - 556482866.5 / rt );
	k[r70b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM14] );
	kTroe0 = 3.8000000000E+34 * exp( -7.27 * lgt - 30208000 / rt );
	kTroeInf = 5.6000000000E+09 * exp( -10042000 / rt );
	fcTroe = 0.2493 * exp( -temp / 98.5 ) + 0.7507 * exp( -temp / 1302 ) + 1 * exp( -4167 / temp );
	k[r71f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM15] );
	kTroe0 = 8.4202583646E+36 * exp( -7.09812 * lgt - 177009237.8 / rt );
	kTroeInf = 1.2408801800E+12 * exp( 0.171883 * lgt - 156843237.8 / rt );
	k[r71b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM15] );
	kTroe0 = 1.4000000000E+24 * exp( -3.86 * lgt - 13891000 / rt );
	kTroeInf = 6.0800000000E+09 * exp( 0.27 * lgt - 1172000 / rt );
	fcTroe = 0.218 * exp( -temp / 207.5 ) + 0.782 * exp( -temp / 2663 ) + 1 * exp( -6095 / temp );
	k[r72f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM16] );
	kTroe0 = 7.7814512533E+29 * exp( -3.99563 * lgt - 480410256.8 / rt );
	kTroeInf = 3.3793731157E+15 * exp( 0.134367 * lgt - 467691256.8 / rt );
	k[r72b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM16] );
	k[r73f] = 3.0000000000E+10;
	k[r73b] = 3.1644553172E+10 * exp( 0.178352 * lgt - 286290250 / rt );
	kTroe0 = 1.2000000000E+36 * exp( -7.62 * lgt - 29162000 / rt );
	kTroeInf = 1.0800000000E+09 * exp( 0.454 * lgt - 7615000 / rt );
	fcTroe = 0.0247 * exp( -temp / 210 ) + 0.9753 * exp( -temp / 984 ) + 1 * exp( -4374 / temp );
	k[r74f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM17] );
	kTroe0 = 2.2421002822E+39 * exp( -7.74131 * lgt - 180499632.3 / rt );
	kTroeInf = 2.0178902539E+12 * exp( 0.332687 * lgt - 158952632.3 / rt );
	k[r74b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM17] );
	k[r75f] = 1.3250000000E+03 * exp( 2.53 * lgt - 51212000 / rt );
	k[r75b] = 5.5718940109E-01 * exp( 3.01587 * lgt - 17784230.95 / rt );
	kTroe0 = 1.9900000000E+35 * exp( -7.08 * lgt - 27970000 / rt );
	kTroeInf = 5.2100000000E+14 * exp( -0.99 * lgt - 6611000 / rt );
	fcTroe = 0.1578 * exp( -temp / 125 ) + 0.8422 * exp( -temp / 2219 ) + 1 * exp( -6882 / temp );
	k[r76f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM18] );
	kTroe0 = 2.4724095869E+41 * exp( -7.28894 * lgt - 449644917.7 / rt );
	kTroeInf = 6.4729919335E+20 * exp( -1.19894 * lgt - 428285917.7 / rt );
	k[r76b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM18] );
	k[r77f] = 2.0000000000E+09;
	k[r77b] = 2.5019333849E+08 * exp( 0.471548 * lgt - 281753855.5 / rt );
	k[r78f] = 1.1500000000E+05 * exp( 1.9 * lgt - 31506000 / rt );
	k[r78b] = 2.1634657872E+01 * exp( 2.45918 * lgt - 42922570.06 / rt );
	k[r79f] = 1.0000000000E+11;
	k[r79b] = 1.1573048940E+05 * exp( 1.56102 * lgt - 71793685.47 / rt );
	k[r80f] = 5.0000000000E+10 * exp( -33472000 / rt );
	k[r80b] = 3.8326334707E+07 * exp( 0.59717 * lgt - 22501718.93 / rt );
	k[r81f] = 1.1300000000E+10 * exp( -14343000 / rt );
	k[r81b] = 3.2307674042E+03 * exp( 1.57654 * lgt - 143469758.4 / rt );
	k[r82f] = 1.0000000000E+10;
	k[r82b] = 3.7197289740E+08 * exp( 0.411006 * lgt - 125789549.2 / rt );
	kTroe0 = 5.0700000000E+21 * exp( -3.42 * lgt - 352920000 / rt );
	kTroeInf = 4.3000000000E+04 * exp( 1.5 * lgt - 333046000 / rt );
	fcTroe = 0.068 * exp( -temp / 197 ) + 0.932 * exp( -temp / 1540 ) + 1 * exp( -10300 / temp );
	k[r83f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM19] );
	kTroe0 = 1.6124629783E+28 * exp( -4.10217 * lgt - 355207490.6 / rt );
	kTroeInf = 1.3675721512E+11 * exp( 0.817835 * lgt - 335333490.6 / rt );
	k[r83b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM19] );
	k[r84f] = 2.1600000000E+05 * exp( 1.51 * lgt - 14351000 / rt );
	k[r84b] = 2.3953722075E+06 * exp( 1.40516 * lgt - 79193784.45 / rt );
	kTroe0 = 2.3000000000E+12 * exp( -0.9 * lgt + 7113000 / rt );
	kTroeInf = 7.4000000000E+10 * exp( -0.37 * lgt );
	fcTroe = 0.2654 * exp( -temp / 94 ) + 0.7346 * exp( -temp / 1756 ) + 1 * exp( -5182 / temp );
	k[r85f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM20] );
	kTroe0 = 8.2846323983E+20 * exp( -1.92053 * lgt - 211222004.5 / rt );
	kTroeInf = 2.6654904238E+19 * exp( -1.39053 * lgt - 218335004.5 / rt );
	k[r85b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM20] );
	k[r86f] = 3.5700000000E+01 * exp( 2.4 * lgt + 8828000 / rt );
	k[r86b] = 6.6567177782E+02 * exp( 2.33023 * lgt - 63902461.72 / rt );
	k[r87f] = 2.9000000000E+10 * exp( 2092000 / rt );
	k[r87b] = 5.2181990618E+10 * exp( 0.257238 * lgt - 291174456.9 / rt );
	k[r88f] = 1.7500000000E+09 * exp( -1339000 / rt );
	k[r88b] = 6.5812925088E+07 * exp( 0.502125 * lgt - 133020927.9 / rt );
	k[r89f] = 5.8000000000E+11 * exp( -39999000 / rt );
	k[r89b] = 2.1812283744E+10 * exp( 0.502125 * lgt - 171680927.9 / rt );
	k[r90f] = 5.0000000000E+10;
	k[r90b] = 1.1579534097E+13 * exp( -0.254545 * lgt - 649863152.8 / rt );
	k[r91f] = 3.0000000000E+10;
	k[r91b] = 8.0293064785E+13 * exp( -0.603071 * lgt - 379941498.5 / rt );
	k[r92f] = 2.0000000000E+10;
	k[r92b] = 1.9744350048E+15 * exp( -0.760056 * lgt - 327774390.3 / rt );
	k[r93f] = 1.1300000000E+04 * exp( 2 * lgt - 12552000 / rt );
	k[r93b] = 7.6976589861E+03 * exp( 2.18763 * lgt - 88713757.49 / rt );
	k[r94f] = 3.0000000000E+10;
	k[r94b] = 2.1593228251E+15 * exp( -0.822343 * lgt - 365077378.3 / rt );
	kTroe0 = 2.7000000000E+32 * exp( -6.3 * lgt - 12970000 / rt );
	kTroeInf = 6.3000000000E+10;
	fcTroe = 0.7895 * exp( -temp / 83.5 ) + 0.2105 * exp( -temp / 5398 ) + 1 * exp( -8370 / temp );
	k[r95f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM21] );
	kTroe0 = 8.0512050525E+40 * exp( -7.13308 * lgt - 405128078.2 / rt );
	kTroeInf = 1.8786145123E+19 * exp( -0.833081 * lgt - 392158078.2 / rt );
	k[r95b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM21] );
	k[r96f] = 5.6000000000E+04 * exp( 1.6 * lgt - 22677000 / rt );
	k[r96b] = 1.4048779004E+03 * exp( 2.01452 * lgt - 56519418.41 / rt );
	k[r97f] = 2.5010000000E+10;
	k[r97b] = 8.6055923823E+08 * exp( 0.476809 * lgt + 3460569.564 / rt );
	k[r98f] = 1.0000000000E+05 * exp( 1.6 * lgt - 13054000 / rt );
	k[r98b] = 6.3655575377E+02 * exp( 2.00672 * lgt - 68833071.27 / rt );
	k[r99f] = 4.7600000000E+04 * exp( 1.228 * lgt - 293000 / rt );
	k[r99b] = 2.6269086456E+11 * exp( -0.0505181 * lgt - 109793079.5 / rt );
	k[r100f] = 5.0000000000E+10;
	k[r100b] = 2.4469331101E+10 * exp( 0.478112 * lgt - 432160700.4 / rt );
	k[r101f] = 3.4300000000E+06 * exp( 1.18 * lgt + 1870000 / rt );
	k[r101b] = 6.3345866067E+04 * exp( 1.52462 * lgt - 126458865.7 / rt );
	k[r102f] = 5.0000000000E+09;
	k[r102b] = 1.2976325138E+09 * exp( 0.626766 * lgt - 375170685.8 / rt );
	k[r103f] = 5.0000000000E+09;
	k[r103b] = 1.4671042236E+08 * exp( 0.567692 * lgt - 403698139.8 / rt );
	k[r104f] = 1.4400000000E+03 * exp( 2 * lgt + 3515000 / rt );
	k[r104b] = 1.1944991455E+02 * exp( 2.10618 * lgt - 88707316.99 / rt );
	k[r105f] = 6.3000000000E+03 * exp( 2 * lgt - 6276000 / rt );
	k[r105b] = 4.6222629959E+03 * exp( 2.16526 * lgt - 69970862.98 / rt );
	k[r106f] = 2.0000000000E+10;
	k[r106b] = 3.0312687340E+14 * exp( -0.743493 * lgt - 213093401.9 / rt );
	k[r107f] = 2.1800000000E-07 * exp( 4.5 * lgt + 4184000 / rt );
	k[r107b] = 1.4748019480E-03 * exp( 3.75588 * lgt - 96488304.27 / rt );
	k[r108f] = 5.0400000000E+02 * exp( 2.3 * lgt - 56484000 / rt );
	k[r108b] = 9.1663502510E+07 * exp( 1.14487 * lgt - 31366755.02 / rt );
	k[r109f] = 3.3700000000E+04 * exp( 2 * lgt - 58576000 / rt );
	k[r109b] = 1.2786721446E+02 * exp( 2.4917 * lgt - 27405.78234 / rt );
	k[r110f] = 4.8300000000E-07 * exp( 4 * lgt + 8368000 / rt );
	k[r110b] = 9.3422522348E-10 * exp( 4.83242 * lgt - 221431062.7 / rt );
	k[r111f] = 5.0000000000E+09;
	k[r111b] = 5.8488027152E+10 * exp( 0.0735165 * lgt - 351133034.5 / rt );
	k[r112f] = 3.6000000000E+03 * exp( 2 * lgt - 10460000 / rt );
	k[r112b] = 1.6788377430E+01 * exp( 2.38103 * lgt - 41875015.4 / rt );
	k[r113f] = 3.5400000000E+03 * exp( 2.12 * lgt - 3640000 / rt );
	k[r113b] = 7.3854116736E+00 * exp( 2.57434 * lgt - 79899354.5 / rt );
	k[r114f] = 7.5000000000E+09 * exp( -8368000 / rt );
	k[r114b] = 6.3754053454E+07 * exp( 0.492334 * lgt - 62240503.38 / rt );
	k[r115f] = 1.3000000000E+08 * exp( 6820000 / rt );
	k[r115b] = 6.2200353343E+09 * exp( -0.244887 * lgt - 154764529 / rt );
	k[r116f] = 4.2000000000E+11 * exp( -50208000 / rt );
	k[r116b] = 2.0095498773E+13 * exp( -0.244887 * lgt - 211792529 / rt );
	k[r117f] = 2.0000000000E+10;
	k[r117b] = 3.7780097921E+11 * exp( 0.00375206 * lgt - 475691730.1 / rt );
	k[r118f] = 1.0000000000E+09;
	k[r118b] = 2.8267421607E+11 * exp( -0.149487 * lgt - 237487385.7 / rt );
	k[r119f] = 2.0000000000E+10;
	k[r119b] = 3.2301498489E+11 * exp( -0.149417 * lgt - 105836008.7 / rt );
	k[r120f] = 1.5000000000E+11 * exp( -98742000 / rt );
	k[r120b] = 1.5839793423E+14 * exp( -0.51471 * lgt - 356159419.3 / rt );
	k[r121f] = 1.0000000000E+09 * exp( -33472000 / rt );
	k[r121b] = 4.9107860468E+08 * exp( -0.157507 * lgt - 30118937.77 / rt );
	k[r122f] = 5.8000000000E+10 * exp( -2410000 / rt );
	k[r122b] = 2.6634134985E+10 * exp( 0.182253 * lgt - 579654497.4 / rt );
	k[r123f] = 5.0000000000E+10;
	k[r123b] = 9.0545188562E+14 * exp( -1.00978 * lgt - 327673077.5 / rt );
	k[r124f] = 5.0000000000E+10;
	k[r124b] = 5.9866880368E+15 * exp( -1.08697 * lgt - 420064090.2 / rt );
	k[r125f] = 3.3000000000E+10;
	k[r125b] = 1.7512987663E+11 * exp( -0.166273 * lgt - 307322843.1 / rt );
	k[r126f] = 1.1070000000E+05 * exp( 1.79 * lgt - 6987000 / rt );
	k[r126b] = 1.8021322225E+06 * exp( 1.49753 * lgt + 4331973.037 / rt );
	k[r127f] = 1.7130000000E+10 * exp( 3159000 / rt );
	k[r127b] = 2.4825041623E+15 * exp( -0.947689 * lgt - 248453632.8 / rt );
	k[r128f] = 4.0000000000E+10;
	k[r128b] = 1.0797267305E+18 * exp( -1.3719 * lgt - 548460717.8 / rt );
	k[r129f] = 3.0000000000E+10;
	k[r129b] = 1.7367155877E+15 * exp( -1.0309 * lgt - 231170101.7 / rt );
	k[r130f] = 6.0000000000E+10;
	k[r130b] = 4.7412070602E+15 * exp( -1.0052 * lgt - 255534157.6 / rt );
	kTroe0 = 2.6900000000E+22 * exp( -3.74 * lgt - 8100000 / rt );
	kTroeInf = 5.0000000000E+10;
	fcTroe = 0.4243 * exp( -temp / 237 ) + 0.5757 * exp( -temp / 1652 ) + 1 * exp( -5069 / temp );
	k[r131f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM22] );
	kTroe0 = 1.2130524545E+32 * exp( -5.18097 * lgt - 320775841.3 / rt );
	kTroeInf = 2.2547443392E+20 * exp( -1.44097 * lgt - 312675841.3 / rt );
	k[r131b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM22] );
	k[r132f] = 3.4000000000E+09 * exp( -2887000 / rt );
	k[r132b] = 1.6489127663E+06 * exp( 0.675447 * lgt - 273328419 / rt );
	k[r133f] = 9.4600000000E+10 * exp( 2155000 / rt );
	k[r133b] = 1.7498844055E+17 * exp( -1.35597 * lgt - 319203631.8 / rt );
	k[r134f] = 5.0000000000E+10;
	k[r134b] = 1.1388201171E+12 * exp( 0.12683 * lgt - 657557391.2 / rt );
	k[r135f] = 1.3200000000E+10 * exp( -6276000 / rt );
	k[r135b] = 2.5592253472E+09 * exp( 0.0911322 * lgt - 317030138.8 / rt );
	k[r136f] = 5.0000000000E+02 * exp( 2 * lgt - 30250000 / rt );
	k[r136b] = 2.2102362917E+05 * exp( 1.48064 * lgt - 61250366.04 / rt );
	k[r137f] = 3.2000000000E+10;
	k[r137b] = 5.3059702311E+16 * exp( -1.07943 * lgt - 559779690.8 / rt );
	k[r138f] = 4.0000000000E+10;
	k[r138b] = 3.3825283043E+17 * exp( -1.2243 * lgt - 275916843.8 / rt );
	k[r139f] = 2.4600000000E+03 * exp( 2 * lgt - 34602000 / rt );
	k[r139b] = 6.2419602881E+02 * exp( 1.9922 * lgt - 56538652.86 / rt );
	kTroe0 = 2.6900000000E+27 * exp( -5.11 * lgt - 29685000 / rt );
	kTroeInf = 8.1000000000E+08 * exp( 0.5 * lgt - 18870000 / rt );
	fcTroe = 0.4093 * exp( -temp / 275 ) + 0.5907 * exp( -temp / 1226 ) + 1 * exp( -5185 / temp );
	k[r140f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM23] );
	kTroe0 = 9.7210519966E+38 * exp( -6.85567 * lgt - 364650095.4 / rt );
	kTroeInf = 2.9271569209E+20 * exp( -1.24567 * lgt - 353835095.4 / rt );
	k[r140b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM23] );
	k[r141f] = 3.0000000000E+10;
	k[r141b] = 3.9791438217E+10 * exp( 0.240947 * lgt - 382586114.2 / rt );
	k[r142f] = 1.5000000000E+10 * exp( -2510000 / rt );
	k[r142b] = 1.0936408744E+10 * exp( -0.0622863 * lgt - 39812987.98 / rt );
	k[r143f] = 9.0000000000E+09 * exp( -2510000 / rt );
	k[r143b] = 6.5618452463E+09 * exp( -0.0622863 * lgt - 39812987.98 / rt );
	k[r144f] = 2.8000000000E+10;
	k[r144b] = 7.4729040691E+05 * exp( 0.261558 * lgt - 282283554.9 / rt );
	k[r145f] = 1.2000000000E+10;
	k[r145b] = 8.3014013421E+08 * exp( 0.506957 * lgt - 780217827.2 / rt );
	k[r146f] = 7.0000000000E+10;
	k[r146b] = 2.2560577672E+13 * exp( -0.581645 * lgt - 68303354.01 / rt );
	kTroe0 = 2.7000000000E+32 * exp( -6.3 * lgt - 12970000 / rt );
	kTroeInf = 2.0000000000E+10;
	fcTroe = 0.8493 * exp( -temp / 134 ) + 0.1507 * exp( -temp / 2383 ) + 1 * exp( -7265 / temp );
	k[r147f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM24] );
	kTroe0 = 2.3398812007E+42 * exp( -7.60989 * lgt - 408588647.8 / rt );
	kTroeInf = 1.7332453339E+20 * exp( -1.30989 * lgt - 395618647.8 / rt );
	k[r147b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM24] );
	k[r148f] = 3.0000000000E+10;
	k[r148b] = 2.1872817488E+10 * exp( -0.0622863 * lgt - 37302987.98 / rt );
	k[r149f] = 1.2000000000E+10 * exp( 2385000 / rt );
	k[r149b] = 7.3985424247E+16 * exp( -1.28658 * lgt - 310834831.8 / rt );
	k[r150f] = 1.6000000000E+10 * exp( 2385000 / rt );
	k[r150b] = 2.9599839164E+09 * exp( -0.0700846 * lgt - 56854640.83 / rt );
	k[r151f] = 9.0000000000E+09;
	k[r151b] = 6.5618452463E+09 * exp( -0.0622863 * lgt - 37302987.98 / rt );
	k[r152f] = 7.0000000000E+09;
	k[r152b] = 5.1036574138E+09 * exp( -0.0622863 * lgt - 37302987.98 / rt );
	k[r153f] = 1.4000000000E+10;
	k[r153b] = 1.8259393135E+08 * exp( 0.456176 * lgt - 255577298.7 / rt );
	k[r154f] = 4.0000000000E+10 * exp( 2301000 / rt );
	k[r154b] = 2.4252938103E+09 * exp( -0.0224647 * lgt - 77418924.07 / rt );
	k[r155f] = 2.6750000000E+10 * exp( -120499000 / rt );
	k[r155b] = 4.4769860885E+12 * exp( -0.476427 * lgt - 5799013.492 / rt );
	k[r156f] = 3.6000000000E+07 * exp( -37405000 / rt );
	k[r156b] = 9.4812040571E+06 * exp( 0.161037 * lgt - 253672691.5 / rt );
	k[r157f] = 2.4500000000E+01 * exp( 2.47 * lgt - 21673000 / rt );
	k[r157b] = 1.4474473694E+02 * exp( 2.5654 * lgt - 97575856.67 / rt );
	kTroe0 = 1.7700000000E+44 * exp( -9.67 * lgt - 26024000 / rt );
	kTroeInf = 2.1200000000E+13 * exp( -0.97 * lgt - 2594000 / rt );
	fcTroe = 0.4675 * exp( -temp / 151 ) + 0.5325 * exp( -temp / 1038 ) + 1 * exp( -4970 / temp );
	k[r158f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM25] );
	kTroe0 = 3.3628508749E+55 * exp( -11.0554 * lgt - 410861540 / rt );
	kTroeInf = 4.0278213868E+24 * exp( -2.35543 * lgt - 387431540 / rt );
	k[r158b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM25] );
	k[r159f] = 4.9900000000E+09 * exp( 0.1 * lgt - 44350000 / rt );
	k[r159b] = 7.6307556299E+14 * exp( -1.07649 * lgt - 7512622.297 / rt );
	k[r160f] = 2.6480000000E+10;
	k[r160b] = 2.0357930432E+12 * exp( 0.0713867 * lgt - 376381629.1 / rt );
	k[r161f] = 3.3200000000E+00 * exp( 2.81 * lgt - 24518000 / rt );
	k[r161b] = 9.6322067040E+00 * exp( 2.74789 * lgt - 97067794.44 / rt );
	k[r162f] = 3.0000000000E+04 * exp( 1.5 * lgt - 41589000 / rt );
	k[r162b] = 3.9093824409E+05 * exp( 1.19946 * lgt - 78032245.72 / rt );
	k[r163f] = 1.0000000000E+04 * exp( 1.5 * lgt - 41589000 / rt );
	k[r163b] = 1.1525974507E+06 * exp( 1.25853 * lgt - 49504791.71 / rt );
	k[r164f] = 2.2700000000E+02 * exp( 2 * lgt - 38493000 / rt );
	k[r164b] = 1.6630129560E+02 * exp( 1.97431 * lgt - 14128944.14 / rt );
	k[r165f] = 6.1400000000E+03 * exp( 1.74 * lgt - 43723000 / rt );
	k[r165b] = 2.0123493121E+03 * exp( 1.78762 * lgt - 64203283.24 / rt );
	k[r166f] = 2.2440000000E+15 * exp( -1 * lgt - 71128000 / rt );
	k[r166b] = 4.2367823753E+11 * exp( -0.767288 * lgt - 5354428.13 / rt );
	k[r167f] = 1.8700000000E+14 * exp( -1 * lgt - 71128000 / rt );
	k[r167b] = 3.5306519794E+10 * exp( -0.767288 * lgt - 5354428.13 / rt );
	k[r168f] = 7.6000000000E+09 * exp( -1674000 / rt );
	k[r168b] = 2.0670122051E+09 * exp( 0.220873 * lgt - 140568243.4 / rt );
	k[r169f] = 1.8000000000E+10 * exp( -3766000 / rt );
	k[r169b] = 2.5961607220E+09 * exp( 0.369528 * lgt - 85670228.82 / rt );
	k[r170f] = 4.2800000000E-16 * exp( 7.6 * lgt + 14770000 / rt );
	k[r170b] = 6.9793035518E-18 * exp( 7.91045 * lgt - 95661682.83 / rt );
	k[r171f] = 5.0000000000E+10 * exp( -6276000 / rt );
	k[r171b] = 1.2397391326E+08 * exp( 0.846372 * lgt - 639220214.1 / rt );
	k[r172f] = 4.0700000000E+02 * exp( 2.4 * lgt - 837000 / rt );
	k[r172b] = 1.1895545457E+06 * exp( 1.80346 * lgt - 124228378.7 / rt );
	k[r173f] = 3.9800000000E+09 * exp( 1004000 / rt );
	k[r173b] = 4.8461114526E+07 * exp( 0.588862 * lgt - 364035088.3 / rt );
	kTroe0 = 7.0000000000E+47 * exp( -9.31 * lgt - 417814000 / rt );
	kTroeInf = 8.0000000000E+12 * exp( 0.44 * lgt - 371414000 / rt );
	fcTroe = 0.2655 * exp( -temp / 180 ) + 0.7345 * exp( -temp / 1035 ) + 1 * exp( -5417 / temp );
	k[r174f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM27] );
	kTroe0 = 1.3284438040E+42 * exp( -8.99601 * lgt - 237584993.2 / rt );
	kTroeInf = 1.5182214903E+07 * exp( 0.753985 * lgt - 191184993.2 / rt );
	k[r174b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM27] );
	k[r175f] = 8.4000000000E+08 * exp( -16213000 / rt );
	k[r175b] = 6.4762207428E+08 * exp( 0.109475 * lgt - 69543182.98 / rt );
	k[r176f] = 1.6000000000E+09 * exp( -3573000 / rt );
	k[r176b] = 4.9419591153E-02 * exp( 1.82258 * lgt - 357650240.4 / rt );
	k[r177f] = 1.0000000000E+10;
	k[r177b] = 5.0507727073E+01 * exp( 1.5678 * lgt - 344881549.9 / rt );
	k[r178f] = 3.5000000000E+10 * exp( -1381000 / rt );
	k[r178b] = 5.5666593174E+10 * exp( 0.125293 * lgt - 315579803.4 / rt );
	k[r179f] = 2.6500000000E+09 * exp( -26778000 / rt );
	k[r179b] = 2.6196114916E+08 * exp( 0.0908078 * lgt - 158468102.7 / rt );
	k[r180f] = 7.3330000000E+10 * exp( -4686000 / rt );
	k[r180b] = 3.6558062974E+12 * exp( -0.34599 * lgt - 208994758.1 / rt );
	k[r181f] = 1.4000000000E+09 * exp( -45229000 / rt );
	k[r181b] = 1.9933960680E+05 * exp( 1.00025 * lgt - 372592819.5 / rt );
	k[r182f] = 2.9000000000E+10 * exp( -96860000 / rt );
	k[r182b] = 2.5664214343E+05 * exp( 0.965769 * lgt - 241715118.8 / rt );
	k[r183f] = 4.4000000000E+11 * exp( -78994000 / rt );
	k[r183b] = 1.2422464360E+05 * exp( 1.43705 * lgt - 333739164.1 / rt );
	k[r184f] = 2.0000000000E+09 * exp( -88115000 / rt );
	k[r184b] = 2.9509703195E+06 * exp( 0.673244 * lgt - 194942824.3 / rt );
	kTroe0 = 6.2000000000E+11 * exp( -234722000 / rt );
	kTroeInf = 1.3000000000E+11 * exp( -249450000 / rt );
	fcTroe = 1 * exp( -0 / temp );
	k[r185f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM28] );
	kTroe0 = 1.2592150686E+03 * exp( 1.12188 * lgt - 64263353.53 / rt );
	kTroeInf = 2.6402896599E+02 * exp( 1.12188 * lgt - 78991353.53 / rt );
	k[r185b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM28] );
	k[r186f] = 2.1100000000E+09 * exp( 2008000 / rt );
	k[r186b] = 7.2973454994E+10 * exp( -0.253856 * lgt - 28607484.22 / rt );
	k[r187f] = 1.0600000000E+14 * exp( -1.41 * lgt );
	k[r187b] = 2.6632684356E+21 * exp( -2.11249 * lgt - 307901955 / rt );
	k[r188f] = 3.9000000000E+09 * exp( 1004000 / rt );
	k[r188b] = 1.0882103930E+07 * exp( 0.580867 * lgt - 188916511 / rt );
	k[r189f] = 1.3200000000E+11 * exp( -1506000 / rt );
	k[r189b] = 7.3031751256E+05 * exp( 1.01766 * lgt - 118807855.6 / rt );
	k[r190f] = 5.0000000000E+10;
	k[r190b] = 2.7076746651E+12 * exp( -0.269105 * lgt - 298056501.2 / rt );
	k[r191f] = 3.2000000000E+10 * exp( -1381000 / rt );
	k[r191b] = 5.8445036912E+10 * exp( 0.111948 * lgt - 103016420.4 / rt );
	k[r192f] = 2.0000000000E+10;
	k[r192b] = 1.9456650491E+14 * exp( -0.735716 * lgt - 76392793.5 / rt );
	k[r193f] = 2.0000000000E+06 * exp( 1.2 * lgt );
	k[r193b] = 4.0508569759E+07 * exp( 1.20711 * lgt - 166478204.9 / rt );
	k[r194f] = 4.6100000000E+02 * exp( 2 * lgt - 27196000 / rt );
	k[r194b] = 8.8925953097E+03 * exp( 1.70108 * lgt - 30970138.07 / rt );
	k[r195f] = 1.2800000000E+03 * exp( 1.5 * lgt - 418000 / rt );
	k[r195b] = 1.3744405767E+02 * exp( 1.66769 * lgt - 225855845.8 / rt );
	k[r196f] = 1.5000000000E+10;
	k[r196b] = 1.2919459202E+12 * exp( -0.143813 * lgt - 612255304.6 / rt );
	k[r197f] = 2.0000000000E+10 * exp( -57948000 / rt );
	k[r197b] = 1.7544816179E+13 * exp( -0.630881 * lgt - 69498009.06 / rt );
	k[r198f] = 2.1600000000E+10 * exp( -0.23 * lgt );
	k[r198b] = 3.7316880802E+10 * exp( -0.0278224 * lgt - 407946546.5 / rt );
	k[r199f] = 4.1600000000E+11 * exp( -0.45 * lgt );
	k[r199b] = 2.5455980630E+18 * exp( -1.68487 * lgt - 153201382.5 / rt );
	k[r200f] = 7.0000000000E+09;
	k[r200b] = 3.0310363517E+08 * exp( 0.250072 * lgt - 42944039.28 / rt );
	k[r201f] = 4.6000000000E+10;
	k[r201b] = 1.9377110617E+13 * exp( -0.485644 * lgt - 119336832.8 / rt );
	k[r202f] = 4.0000000000E+10 * exp( -15272000 / rt );
	k[r202b] = 2.9122306596E+09 * exp( 0.285135 * lgt - 66103716.55 / rt );
	k[r203f] = 9.0000000000E+04 * exp( 1.5 * lgt + 1925000 / rt );
	k[r203b] = 7.2665378999E+04 * exp( 1.6803 * lgt - 113749501 / rt );
	k[r204f] = 3.3000000000E+08;
	k[r204b] = 9.6125078621E+04 * exp( 0.283779 * lgt - 31168478.84 / rt );
	k[r205f] = 1.3000000000E+11 * exp( -0.11 * lgt - 20836000 / rt );
	k[r205b] = 3.7867455214E+07 * exp( 0.173779 * lgt - 52004478.84 / rt );
	k[r206f] = 5.0000000000E+09;
	k[r206b] = 2.0980165850E+09 * exp( 0.271941 * lgt - 235836294.1 / rt );
	k[r207f] = 2.5000000000E+10;
	k[r207b] = 1.0123009331E+09 * exp( 0.598951 * lgt - 456372289.4 / rt );
	k[r208f] = 7.0000000000E+10;
	k[r208b] = 1.6406505346E+09 * exp( 0.396774 * lgt - 48425742.82 / rt );
	k[r209f] = 5.0000000000E+10;
	k[r209b] = 3.4041783577E+09 * exp( 0.634014 * lgt - 464259966.6 / rt );
	k[r210f] = 2.0000000000E+10;
	k[r210b] = 1.5100507829E+10 * exp( 0.529179 * lgt - 529102751.1 / rt );
	k[r211f] = 2.5000000000E+10;
	k[r211b] = 2.9652759674E+12 * exp( 0.122454 * lgt - 473323679.8 / rt );
	k[r212f] = 8.9500000000E+13 * exp( -1.32 * lgt - 3096000 / rt );
	k[r212b] = 2.2350204463E+18 * exp( -1.47144 * lgt - 206636102.8 / rt );
	k[r213f] = 2.5000000000E+10;
	k[r213b] = 1.3916448087E+08 * exp( 0.466611 * lgt - 221663707.7 / rt );
	k[r214f] = 4.5000000000E+08 * exp( 0.72 * lgt - 2761000 / rt );
	k[r214b] = 4.2118566583E+06 * exp( 1.22167 * lgt - 232312385 / rt );
	k[r215f] = 1.3000000000E+04 * exp( 1.9 * lgt + 3975000 / rt );
	k[r215b] = 1.3493470881E+03 * exp( 2.29684 * lgt - 290419169.5 / rt );
	k[r216f] = 1.0000000000E+10 * exp( -54392000 / rt );
	k[r216b] = 5.7684306980E+08 * exp( 0.1396 * lgt - 55519712.52 / rt );
	k[r217f] = 7.7000000000E+10;
	k[r217b] = 2.7037695910E+11 * exp( 0.0464021 * lgt - 325702214.2 / rt );
	k[r218f] = 4.0000000000E+10;
	k[r218b] = 1.0396721037E+17 * exp( -1.42788 * lgt - 133363208.4 / rt );
	k[r219f] = 8.0000000000E+09 * exp( -31213000 / rt );
	k[r219b] = 2.4158642845E+12 * exp( -0.735333 * lgt - 55744749.09 / rt );
	k[r220f] = 6.1400000000E+09 * exp( 1841000 / rt );
	k[r220b] = 3.1644212506E+13 * exp( -0.99108 * lgt - 58903552.94 / rt );
	k[r221f] = 2.1000000000E+10 * exp( -19707000 / rt );
	k[r221b] = 7.0326838802E+13 * exp( -0.840168 * lgt - 109081533.5 / rt );
	k[r222f] = 2.3500000000E+10;
	k[r222b] = 1.5827484909E+06 * exp( 1.12829 * lgt - 396647763.9 / rt );
	k[r223f] = 5.4000000000E+10;
	k[r223b] = 6.7160099146E+04 * exp( 1.3974 * lgt - 98591262.68 / rt );
	k[r224f] = 2.5000000000E+09;
	k[r224b] = 1.2112594941E+03 * exp( 0.813118 * lgt + 28556046.6 / rt );
	k[r225f] = 2.0000000000E+10;
	k[r225b] = 2.1424003974E+06 * exp( 1.25358 * lgt - 710846567.3 / rt );
	k[r226f] = 2.0000000000E+09 * exp( -83680000 / rt );
	k[r226b] = 1.4740139035E+09 * exp( 0.28657 * lgt - 517209188 / rt );
	k[r227f] = 8.8000000000E+13 * exp( -0.5 * lgt - 200832000 / rt );
	k[r227b] = 8.5522101884E+05 * exp( 0.659108 * lgt + 32032804.7 / rt );
	k[r228f] = 2.8500000000E+14 * exp( -1.52 * lgt - 3096000 / rt );
	k[r228b] = 2.1689968935E+15 * exp( -1.35748 * lgt - 254888645.2 / rt );
	k[r229f] = 5.7000000000E+15 * exp( -2 * lgt - 3347000 / rt );
	k[r229b] = 6.7589934526E+16 * exp( -1.67895 * lgt - 619384888.7 / rt );
	k[r230f] = 1.0400000000E+26 * exp( -3.3 * lgt - 529694000 / rt );
	k[r230b] = 1.3286543200E+20 * exp( -2.81007 * lgt - 7227978.666 / rt );
	k[r231f] = 1.1070000000E+01 * exp( 2.64 * lgt - 20836000 / rt );
	k[r231b] = 5.1098664034E+03 * exp( 2.01723 * lgt - 56936997.57 / rt );
	k[r232f] = 2.7670000000E+00 * exp( 2.64 * lgt - 20836000 / rt );
	k[r232b] = 1.5885052436E-03 * exp( 3.41462 * lgt - 155528260.2 / rt );
	k[r233f] = 2.1340000000E+06 * exp( 1.58 * lgt - 111294000 / rt );
	k[r233b] = 3.7898318393E+02 * exp( 2.38511 * lgt - 14031789.19 / rt );
	k[r234f] = 1.1000000000E+03 * exp( 2.03 * lgt - 55940000 / rt );
	k[r234b] = 4.6430276446E+06 * exp( 1.33782 * lgt - 26668814.91 / rt );
	k[r235f] = 4.4000000000E+00 * exp( 2.26 * lgt - 26778000 / rt );
	k[r235b] = 2.5580724323E+05 * exp( 1.22653 * lgt - 105221801.3 / rt );
	k[r236f] = 1.6000000000E-01 * exp( 2.56 * lgt - 37656000 / rt );
	k[r236b] = 2.1213209026E-03 * exp( 3.08455 * lgt - 129404221 / rt );
	k[r237f] = 1.4000000000E+20 * exp( -3.4 * lgt - 7950000 / rt );
	k[r237b] = 6.9763425470E+23 * exp( -3.52088 * lgt - 114266547.5 / rt );
	k[r238f] = 6.0000000000E+10 * exp( -1674000 / rt );
	k[r238b] = 3.2249699964E+07 * exp( 0.864641 * lgt - 324705795 / rt );
	k[r239f] = 6.3000000000E+10 * exp( -192548000 / rt );
	k[r239b] = 5.2402811303E+10 * exp( -0.0802491 * lgt + 1798622.854 / rt );
	k[r240f] = 2.8570000000E+05 * exp( 1.1 * lgt - 85354000 / rt );
	k[r240b] = 4.0587716300E+08 * exp( 0.414003 * lgt - 77778172.28 / rt );
	kTroe0 = 1.3000000000E+19 * exp( -3.16 * lgt - 3096000 / rt );
	kTroeInf = 3.1000000000E+09 * exp( 0.15 * lgt );
	fcTroe = 0.333 * exp( -temp / 235 ) + 0.667 * exp( -temp / 2117 ) + 1 * exp( -4536 / temp );
	k[r241f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM35] );
	kTroe0 = 8.2668740778E+28 * exp( -4.81452 * lgt - 142330827 / rt );
	kTroeInf = 1.9713315109E+19 * exp( -1.50452 * lgt - 139234827 / rt );
	k[r241b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM35] );
	k[r242f] = 1.0000000000E+10 * exp( -309616000 / rt );
	k[r242b] = 4.7780169452E+11 * exp( -0.505477 * lgt - 211723724.9 / rt );
	k[r243f] = 1.0000000000E+08 * exp( -271960000 / rt );
	k[r243b] = 3.4836230865E+09 * exp( -0.567763 * lgt - 211370712.9 / rt );
	k[r244f] = 1.9000000000E+10;
	k[r244b] = 2.5135888245E+10 * exp( 0.0450435 * lgt - 119852180.5 / rt );
	k[r245f] = 2.9000000000E+10;
	k[r245b] = 1.3471550636E+11 * exp( 0.0914455 * lgt - 445554394.7 / rt );
	k[r246f] = 5.0000000000E+10;
	k[r246b] = 1.1297464329E+14 * exp( -0.560704 * lgt - 306622975.7 / rt );
	k[r247f] = 2.0000000000E+10;
	k[r247b] = 2.0859452002E+16 * exp( -1.18348 * lgt - 342723973.2 / rt );
	k[r248f] = 3.0000000000E+10;
	k[r248b] = 1.6105586461E+12 * exp( -0.257081 * lgt - 175632740.4 / rt );
	k[r249f] = 3.1000000000E+14 * exp( -1.38 * lgt - 5314000 / rt );
	k[r249b] = 1.4877197367E+21 * exp( -2.71677 * lgt - 393812072.8 / rt );
	k[r250f] = 2.9000000000E+11 * exp( -0.69 * lgt - 3180000 / rt );
	k[r250b] = 2.3938518404E+13 * exp( -0.993299 * lgt - 313234271.4 / rt );
	k[r251f] = 3.8000000000E+10 * exp( -0.36 * lgt - 2427000 / rt );
	k[r251b] = 6.9498400218E+17 * exp( -1.9287 * lgt - 102420234.8 / rt );
	k[r252f] = 3.1000000000E+14 * exp( -1.38 * lgt - 5314000 / rt );
	k[r252b] = 1.0846874091E+21 * exp( -2.77906 * lgt - 431115060.8 / rt );
	k[r253f] = 2.9000000000E+11 * exp( -0.69 * lgt - 3180000 / rt );
	k[r253b] = 1.7453428133E+13 * exp( -1.05559 * lgt - 350537259.4 / rt );
	k[r254f] = 3.8000000000E+10 * exp( -0.36 * lgt - 2427000 / rt );
	k[r254b] = 5.0670860788E+17 * exp( -1.99098 * lgt - 139723222.8 / rt );
	k[r255f] = 9.6000000000E+10 * exp( -120499000 / rt );
	k[r255b] = 1.9880214074E+11 * exp( 0.111224 * lgt - 464395689.8 / rt );
	k[r256f] = 1.0000000000E+09 * exp( -91002000 / rt );
	k[r256b] = 3.9811706690E+09 * exp( -0.255053 * lgt - 43280965.12 / rt );
	k[r257f] = 2.2000000000E+10;
	k[r257b] = 2.4301979813E+02 * exp( 1.59933 * lgt - 600136910.1 / rt );
	k[r258f] = 2.0000000000E+09;
	k[r258b] = 3.0964548596E+06 * exp( 0.937705 * lgt - 482701913.9 / rt );
	k[r259f] = 1.2000000000E+10;
	k[r259b] = 1.0014507077E+01 * exp( 1.48825 * lgt - 168088016 / rt );
	k[r260f] = 1.2000000000E+10;
	k[r260b] = 5.0505660659E+03 * exp( 1.05145 * lgt - 240706671.5 / rt );
	k[r261f] = 1.0000000000E+11;
	k[r261b] = 5.9835730142E+04 * exp( 1.71229 * lgt - 282537687.7 / rt );
	k[r262f] = 9.8000000000E+04 * exp( 1.41 * lgt - 35564000 / rt );
	k[r262b] = 5.3405185205E+03 * exp( 1.93958 * lgt - 201312538.4 / rt );
	k[r263f] = 1.5000000000E+05 * exp( 1.57 * lgt - 184096000 / rt );
	k[r263b] = 1.4409492112E+01 * exp( 2.64238 * lgt - 316737252.4 / rt );
	k[r264f] = 2.2000000000E+03 * exp( 2.11 * lgt - 47698000 / rt );
	k[r264b] = 1.7467246323E+01 * exp( 2.5207 * lgt - 5355196.22 / rt );
	k[r265f] = 2.2500000000E+04 * exp( 1.7 * lgt - 15899000 / rt );
	k[r265b] = 5.1310795266E-03 * exp( 3.25802 * lgt - 29203419.61 / rt );
	k[r266f] = 1.0500000000E+02 * exp( 2.5 * lgt - 55647000 / rt );
	k[r266b] = 1.4017279608E+00 * exp( 2.94576 * lgt - 21191873.49 / rt );
	k[r267f] = 4.6500000000E+09 * exp( -28660000 / rt );
	k[r267b] = 6.8840916793E+08 * exp( 0.340929 * lgt - 59047657.94 / rt );
	k[r268f] = 1.5500000000E+09 * exp( -28660000 / rt );
	k[r268b] = 1.9507245204E+09 * exp( 0.279506 * lgt - 151464499.1 / rt );
	k[r269f] = 1.1800000000E+13 * exp( -354468000 / rt );
	k[r269b] = 8.3821316893E+02 * exp( 1.49292 * lgt + 14487351.62 / rt );
	k[r270f] = 2.1000000000E+12 * exp( -0.69 * lgt - 11924000 / rt );
	k[r270b] = 5.5104556881E+11 * exp( -0.458075 * lgt - 300428838 / rt );
	k[r271f] = 2.7000000000E+08 * exp( 0.18 * lgt - 8870000 / rt );
	k[r271b] = 1.2186298809E+03 * exp( 1.4454 * lgt - 218931036.6 / rt );
	k[r272f] = 1.7000000000E+11 * exp( -0.75 * lgt - 12092000 / rt );
	k[r272b] = 1.0172867051E+04 * exp( 1.03995 * lgt - 313901257.6 / rt );
	k[r273f] = 2.0000000000E+04 * exp( 2 * lgt - 8368000 / rt );
	k[r273b] = 2.7547460710E+05 * exp( 1.6587 * lgt - 116082986.4 / rt );
	k[r274f] = 2.3500000000E+10;
	k[r274b] = 3.6265221753E+11 * exp( -0.0699648 * lgt - 209089908.3 / rt );
	k[r275f] = 6.1000000000E+11 * exp( -0.31 * lgt - 1213000 / rt );
	k[r275b] = 1.2107155556E+14 * exp( -0.911043 * lgt - 157800723.2 / rt );
	k[r276f] = 3.7000000000E+09 * exp( 0.15 * lgt + 377000 / rt );
	k[r276b] = 3.4445602388E+10 * exp( 0.0200691 * lgt - 482985663.5 / rt );
	k[r277f] = 5.4000000000E+02 * exp( 2.4 * lgt - 41484000 / rt );
	k[r277b] = 2.8686505207E+00 * exp( 2.76177 * lgt - 18443175.52 / rt );
	k[r278f] = 5.0000000000E+04 * exp( 1.6 * lgt - 3996000 / rt );
	k[r278b] = 2.9455957348E+03 * exp( 1.85693 * lgt - 45797959.97 / rt );
	k[r279f] = 9.4000000000E+03 * exp( 1.94 * lgt - 27029000 / rt );
	k[r279b] = 2.9698811025E+01 * exp( 2.2667 * lgt + 3899501.745 / rt );
	k[r280] = 6.0000000000E+09;
	k[r281] = 3.8000000000E+10;
	k[r282f] = 5.0000000000E+14 * exp( -370000000 / rt );
	k[r282b] = 1.7000000000E+10 * exp( -19880000 / rt );
	k[r283] = 1.0000000000E+09;
	k[r284] = 1.0000000000E+09;
	k[r285] = 1.0000000000E+09;
	k[r286] = 1.0000000000E+09;
	k[r287f] = 3.9800000000E+13 * exp( -293100000 / rt );
	k[r287b] = 1.2670000000E+10 * exp( -32480000 / rt );
	k[r288] = 1.0000000000E+10;
	k[r289f] = 3.1500000000E+15 * exp( -359300000 / rt );
	k[r289b] = 2.5110000000E+09 * exp( 34690000 / rt );
	k[r290] = 5.0000000000E+09 * exp( -6300000 / rt );
	k[r291] = 9.6000000000E+13 * exp( -129800000 / rt );
	k[r292f] = 1.2500000000E+14 * exp( -154900000 / rt );
	k[r292b] = 4.6090000000E+11 * exp( -21490000 / rt );
	k[r293] = 6.3000000000E+13 * exp( -154500000 / rt );
	k[r294] = 1.0000000000E+09 * exp( -20900000 / rt );
	k[r295] = 1.3000000000E+11 * exp( -40600000 / rt );
	k[r296] = 1.0000000000E+11 * exp( -34900000 / rt );
	k[r297] = 3.0000000000E+10 * exp( -24100000 / rt );
	k[r298] = 2.6000000000E+10 * exp( -18700000 / rt );
	k[r299] = 3.7000000000E+09 * exp( -6900000 / rt );
	k[r300] = 2.8000000000E+09 * exp( -3600000 / rt );


	w[r1f] = k[r1f] * c[sO] * c[sO] * M[mM1];
	w[r1b] = k[r1b] * c[sO2] * M[mM1];
	w[r2f] = k[r2f] * c[sO] * c[sH] * M[mM2];
	w[r2b] = k[r2b] * c[sOH] * M[mM2];
	w[r3f] = k[r3f] * c[sO] * c[sH2];
	w[r3b] = k[r3b] * c[sOH] * c[sH];
	w[r4f] = k[r4f] * c[sO] * c[sHO2];
	w[r4b] = k[r4b] * c[sO2] * c[sOH];
	w[r5f] = k[r5f] * c[sO] * c[sH2O2];
	w[r5b] = k[r5b] * c[sHO2] * c[sOH];
	w[r6f] = k[r6f] * c[sO] * c[sCH];
	w[r6b] = k[r6b] * c[sCO] * c[sH];
	w[r7f] = k[r7f] * c[sO] * c[sCH2];
	w[r7b] = k[r7b] * c[sHCO] * c[sH];
	w[r8f] = k[r8f] * c[sO] * c[sCH2GSGXCH2];
	w[r8b] = k[r8b] * c[sCO] * c[sH2];
	w[r9f] = k[r9f] * c[sO] * c[sCH2GSGXCH2];
	w[r9b] = k[r9b] * c[sHCO] * c[sH];
	w[r10f] = k[r10f] * c[sO] * c[sCH3];
	w[r10b] = k[r10b] * c[sCH2O] * c[sH];
	w[r11f] = k[r11f] * c[sO] * c[sCH4];
	w[r11b] = k[r11b] * c[sCH3] * c[sOH];
	w[r12f] = k[r12f] * c[sO] * c[sCO] * M[mM3];
	w[r12b] = k[r12b] * c[sCO2] * M[mM3];
	w[r13f] = k[r13f] * c[sO] * c[sHCO];
	w[r13b] = k[r13b] * c[sCO] * c[sOH];
	w[r14f] = k[r14f] * c[sO] * c[sHCO];
	w[r14b] = k[r14b] * c[sCO2] * c[sH];
	w[r15f] = k[r15f] * c[sO] * c[sCH2O];
	w[r15b] = k[r15b] * c[sHCO] * c[sOH];
	w[r16f] = k[r16f] * c[sO] * c[sCH2OH];
	w[r16b] = k[r16b] * c[sCH2O] * c[sOH];
	w[r17f] = k[r17f] * c[sO] * c[sCH3O];
	w[r17b] = k[r17b] * c[sCH2O] * c[sOH];
	w[r18f] = k[r18f] * c[sO] * c[sCH3OH];
	w[r18b] = k[r18b] * c[sCH2OH] * c[sOH];
	w[r19f] = k[r19f] * c[sO] * c[sCH3OH];
	w[r19b] = k[r19b] * c[sCH3O] * c[sOH];
	w[r20f] = k[r20f] * c[sO] * c[sC2H];
	w[r20b] = k[r20b] * c[sCO] * c[sCH];
	w[r21f] = k[r21f] * c[sO] * c[sC2H2];
	w[r21b] = k[r21b] * c[sHCCO] * c[sH];
	w[r22f] = k[r22f] * c[sO] * c[sC2H2];
	w[r22b] = k[r22b] * c[sC2H] * c[sOH];
	w[r23f] = k[r23f] * c[sO] * c[sC2H2];
	w[r23b] = k[r23b] * c[sCH2] * c[sCO];
	w[r24f] = k[r24f] * c[sO] * c[sC2H3];
	w[r24b] = k[r24b] * c[sCH2CO] * c[sH];
	w[r25f] = k[r25f] * c[sO] * c[sC2H4];
	w[r25b] = k[r25b] * c[sHCO] * c[sCH3];
	w[r26f] = k[r26f] * c[sO] * c[sC2H5];
	w[r26b] = k[r26b] * c[sCH2O] * c[sCH3];
	w[r27f] = k[r27f] * c[sO] * c[sC2H6];
	w[r27b] = k[r27b] * c[sC2H5] * c[sOH];
	w[r28f] = k[r28f] * c[sO] * c[sHCCO];
	w[r28b] = k[r28b] * c[sCO] * c[sCO] * c[sH];
	w[r29f] = k[r29f] * c[sO] * c[sCH2CO];
	w[r29b] = k[r29b] * c[sHCCO] * c[sOH];
	w[r30f] = k[r30f] * c[sO] * c[sCH2CO];
	w[r30b] = k[r30b] * c[sCO2] * c[sCH2];
	w[r31f] = k[r31f] * c[sO2] * c[sCO];
	w[r31b] = k[r31b] * c[sCO2] * c[sO];
	w[r32f] = k[r32f] * c[sO2] * c[sCH2O];
	w[r32b] = k[r32b] * c[sHCO] * c[sHO2];
	w[r33f] = k[r33f] * c[sH] * c[sO2] * M[mM4];
	w[r33b] = k[r33b] * c[sHO2] * M[mM4];
	w[r34f] = k[r34f] * c[sH] * c[sO2] * c[sO2];
	w[r34b] = k[r34b] * c[sO2] * c[sHO2];
	w[r35f] = k[r35f] * c[sH] * c[sO2] * c[sH2O];
	w[r35b] = k[r35b] * c[sH2O] * c[sHO2];
	w[r36f] = k[r36f] * c[sH] * c[sO2] * c[sN2];
	w[r36b] = k[r36b] * c[sN2] * c[sHO2];
	w[r37f] = k[r37f] * c[sH] * c[sO2] * c[sAR];
	w[r37b] = k[r37b] * c[sAR] * c[sHO2];
	w[r38f] = k[r38f] * c[sH] * c[sO2];
	w[r38b] = k[r38b] * c[sOH] * c[sO];
	w[r39f] = k[r39f] * c[sH] * c[sH] * M[mM5];
	w[r39b] = k[r39b] * c[sH2] * M[mM5];
	w[r40f] = k[r40f] * c[sH] * c[sH] * c[sH2];
	w[r40b] = k[r40b] * c[sH2] * c[sH2];
	w[r41f] = k[r41f] * c[sH] * c[sH] * c[sH2O];
	w[r41b] = k[r41b] * c[sH2O] * c[sH2];
	w[r42f] = k[r42f] * c[sH] * c[sH] * c[sCO2];
	w[r42b] = k[r42b] * c[sCO2] * c[sH2];
	w[r43f] = k[r43f] * c[sH] * c[sOH] * M[mM6];
	w[r43b] = k[r43b] * c[sH2O] * M[mM6];
	w[r44f] = k[r44f] * c[sH] * c[sHO2];
	w[r44b] = k[r44b] * c[sH2O] * c[sO];
	w[r45f] = k[r45f] * c[sH] * c[sHO2];
	w[r45b] = k[r45b] * c[sH2] * c[sO2];
	w[r46f] = k[r46f] * c[sH] * c[sHO2];
	w[r46b] = k[r46b] * c[sOH] * c[sOH];
	w[r47f] = k[r47f] * c[sH] * c[sH2O2];
	w[r47b] = k[r47b] * c[sH2] * c[sHO2];
	w[r48f] = k[r48f] * c[sH] * c[sH2O2];
	w[r48b] = k[r48b] * c[sH2O] * c[sOH];
	w[r49f] = k[r49f] * c[sH] * c[sCH];
	w[r49b] = k[r49b] * c[sH2] * c[sC];
	w[r50f] = k[r50f] * c[sH] * c[sCH2];
	w[r50b] = k[r50b] * c[sCH3];
	w[r51f] = k[r51f] * c[sH] * c[sCH2GSGXCH2];
	w[r51b] = k[r51b] * c[sH2] * c[sCH];
	w[r52f] = k[r52f] * c[sH] * c[sCH3];
	w[r52b] = k[r52b] * c[sCH4];
	w[r53f] = k[r53f] * c[sH] * c[sCH4];
	w[r53b] = k[r53b] * c[sH2] * c[sCH3];
	w[r54f] = k[r54f] * c[sH] * c[sHCO];
	w[r54b] = k[r54b] * c[sCH2O];
	w[r55f] = k[r55f] * c[sH] * c[sHCO];
	w[r55b] = k[r55b] * c[sCO] * c[sH2];
	w[r56f] = k[r56f] * c[sH] * c[sCH2O];
	w[r56b] = k[r56b] * c[sCH2OH];
	w[r57f] = k[r57f] * c[sH] * c[sCH2O];
	w[r57b] = k[r57b] * c[sCH3O];
	w[r58f] = k[r58f] * c[sH] * c[sCH2O];
	w[r58b] = k[r58b] * c[sH2] * c[sHCO];
	w[r59f] = k[r59f] * c[sH] * c[sCH2OH];
	w[r59b] = k[r59b] * c[sCH3OH];
	w[r60f] = k[r60f] * c[sH] * c[sCH2OH];
	w[r60b] = k[r60b] * c[sCH2O] * c[sH2];
	w[r61f] = k[r61f] * c[sH] * c[sCH2OH];
	w[r61b] = k[r61b] * c[sCH3] * c[sOH];
	w[r62f] = k[r62f] * c[sH] * c[sCH2OH];
	w[r62b] = k[r62b] * c[sH2O] * c[sCH2GSGXCH2];
	w[r63f] = k[r63f] * c[sH] * c[sCH3O];
	w[r63b] = k[r63b] * c[sCH3OH];
	w[r64f] = k[r64f] * c[sH] * c[sCH3O];
	w[r64b] = k[r64b] * c[sCH2OH] * c[sH];
	w[r65f] = k[r65f] * c[sH] * c[sCH3O];
	w[r65b] = k[r65b] * c[sCH2O] * c[sH2];
	w[r66f] = k[r66f] * c[sH] * c[sCH3O];
	w[r66b] = k[r66b] * c[sCH3] * c[sOH];
	w[r67f] = k[r67f] * c[sH] * c[sCH3O];
	w[r67b] = k[r67b] * c[sH2O] * c[sCH2GSGXCH2];
	w[r68f] = k[r68f] * c[sH] * c[sCH3OH];
	w[r68b] = k[r68b] * c[sH2] * c[sCH2OH];
	w[r69f] = k[r69f] * c[sH] * c[sCH3OH];
	w[r69b] = k[r69b] * c[sH2] * c[sCH3O];
	w[r70f] = k[r70f] * c[sH] * c[sC2H];
	w[r70b] = k[r70b] * c[sC2H2];
	w[r71f] = k[r71f] * c[sH] * c[sC2H2];
	w[r71b] = k[r71b] * c[sC2H3];
	w[r72f] = k[r72f] * c[sH] * c[sC2H3];
	w[r72b] = k[r72b] * c[sC2H4];
	w[r73f] = k[r73f] * c[sH] * c[sC2H3];
	w[r73b] = k[r73b] * c[sC2H2] * c[sH2];
	w[r74f] = k[r74f] * c[sH] * c[sC2H4];
	w[r74b] = k[r74b] * c[sC2H5];
	w[r75f] = k[r75f] * c[sH] * c[sC2H4];
	w[r75b] = k[r75b] * c[sH2] * c[sC2H3];
	w[r76f] = k[r76f] * c[sH] * c[sC2H5];
	w[r76b] = k[r76b] * c[sC2H6];
	w[r77f] = k[r77f] * c[sH] * c[sC2H5];
	w[r77b] = k[r77b] * c[sC2H4] * c[sH2];
	w[r78f] = k[r78f] * c[sH] * c[sC2H6];
	w[r78b] = k[r78b] * c[sH2] * c[sC2H5];
	w[r79f] = k[r79f] * c[sH] * c[sHCCO];
	w[r79b] = k[r79b] * c[sCO] * c[sCH2GSGXCH2];
	w[r80f] = k[r80f] * c[sH] * c[sCH2CO];
	w[r80b] = k[r80b] * c[sH2] * c[sHCCO];
	w[r81f] = k[r81f] * c[sH] * c[sCH2CO];
	w[r81b] = k[r81b] * c[sCO] * c[sCH3];
	w[r82f] = k[r82f] * c[sH] * c[sHCCOH];
	w[r82b] = k[r82b] * c[sCH2CO] * c[sH];
	w[r83f] = k[r83f] * c[sH2] * c[sCO];
	w[r83b] = k[r83b] * c[sCH2O];
	w[r84f] = k[r84f] * c[sOH] * c[sH2];
	w[r84b] = k[r84b] * c[sH2O] * c[sH];
	w[r85f] = k[r85f] * c[sOH] * c[sOH];
	w[r85b] = k[r85b] * c[sH2O2];
	w[r86f] = k[r86f] * c[sOH] * c[sOH];
	w[r86b] = k[r86b] * c[sH2O] * c[sO];
	w[r87f] = k[r87f] * c[sOH] * c[sHO2];
	w[r87b] = k[r87b] * c[sH2O] * c[sO2];
	w[r88f] = k[r88f] * c[sOH] * c[sH2O2];
	w[r88b] = k[r88b] * c[sH2O] * c[sHO2];
	w[r89f] = k[r89f] * c[sOH] * c[sH2O2];
	w[r89b] = k[r89b] * c[sH2O] * c[sHO2];
	w[r90f] = k[r90f] * c[sOH] * c[sC];
	w[r90b] = k[r90b] * c[sCO] * c[sH];
	w[r91f] = k[r91f] * c[sOH] * c[sCH];
	w[r91b] = k[r91b] * c[sHCO] * c[sH];
	w[r92f] = k[r92f] * c[sOH] * c[sCH2];
	w[r92b] = k[r92b] * c[sCH2O] * c[sH];
	w[r93f] = k[r93f] * c[sOH] * c[sCH2];
	w[r93b] = k[r93b] * c[sH2O] * c[sCH];
	w[r94f] = k[r94f] * c[sOH] * c[sCH2GSGXCH2];
	w[r94b] = k[r94b] * c[sCH2O] * c[sH];
	w[r95f] = k[r95f] * c[sOH] * c[sCH3];
	w[r95b] = k[r95b] * c[sCH3OH];
	w[r96f] = k[r96f] * c[sOH] * c[sCH3];
	w[r96b] = k[r96b] * c[sH2O] * c[sCH2];
	w[r97f] = k[r97f] * c[sOH] * c[sCH3];
	w[r97b] = k[r97b] * c[sH2O] * c[sCH2GSGXCH2];
	w[r98f] = k[r98f] * c[sOH] * c[sCH4];
	w[r98b] = k[r98b] * c[sH2O] * c[sCH3];
	w[r99f] = k[r99f] * c[sOH] * c[sCO];
	w[r99b] = k[r99b] * c[sCO2] * c[sH];
	w[r100f] = k[r100f] * c[sOH] * c[sHCO];
	w[r100b] = k[r100b] * c[sCO] * c[sH2O];
	w[r101f] = k[r101f] * c[sOH] * c[sCH2O];
	w[r101b] = k[r101b] * c[sH2O] * c[sHCO];
	w[r102f] = k[r102f] * c[sOH] * c[sCH2OH];
	w[r102b] = k[r102b] * c[sCH2O] * c[sH2O];
	w[r103f] = k[r103f] * c[sOH] * c[sCH3O];
	w[r103b] = k[r103b] * c[sCH2O] * c[sH2O];
	w[r104f] = k[r104f] * c[sOH] * c[sCH3OH];
	w[r104b] = k[r104b] * c[sH2O] * c[sCH2OH];
	w[r105f] = k[r105f] * c[sOH] * c[sCH3OH];
	w[r105b] = k[r105b] * c[sH2O] * c[sCH3O];
	w[r106f] = k[r106f] * c[sOH] * c[sC2H];
	w[r106b] = k[r106b] * c[sHCCO] * c[sH];
	w[r107f] = k[r107f] * c[sOH] * c[sC2H2];
	w[r107b] = k[r107b] * c[sCH2CO] * c[sH];
	w[r108f] = k[r108f] * c[sOH] * c[sC2H2];
	w[r108b] = k[r108b] * c[sHCCOH] * c[sH];
	w[r109f] = k[r109f] * c[sOH] * c[sC2H2];
	w[r109b] = k[r109b] * c[sH2O] * c[sC2H];
	w[r110f] = k[r110f] * c[sOH] * c[sC2H2];
	w[r110b] = k[r110b] * c[sCO] * c[sCH3];
	w[r111f] = k[r111f] * c[sOH] * c[sC2H3];
	w[r111b] = k[r111b] * c[sC2H2] * c[sH2O];
	w[r112f] = k[r112f] * c[sOH] * c[sC2H4];
	w[r112b] = k[r112b] * c[sH2O] * c[sC2H3];
	w[r113f] = k[r113f] * c[sOH] * c[sC2H6];
	w[r113b] = k[r113b] * c[sH2O] * c[sC2H5];
	w[r114f] = k[r114f] * c[sOH] * c[sCH2CO];
	w[r114b] = k[r114b] * c[sH2O] * c[sHCCO];
	w[r115f] = k[r115f] * c[sHO2] * c[sHO2];
	w[r115b] = k[r115b] * c[sH2O2] * c[sO2];
	w[r116f] = k[r116f] * c[sHO2] * c[sHO2];
	w[r116b] = k[r116b] * c[sH2O2] * c[sO2];
	w[r117f] = k[r117f] * c[sHO2] * c[sCH2];
	w[r117b] = k[r117b] * c[sCH2O] * c[sOH];
	w[r118f] = k[r118f] * c[sHO2] * c[sCH3];
	w[r118b] = k[r118b] * c[sCH4] * c[sO2];
	w[r119f] = k[r119f] * c[sHO2] * c[sCH3];
	w[r119b] = k[r119b] * c[sCH3O] * c[sOH];
	w[r120f] = k[r120f] * c[sHO2] * c[sCO];
	w[r120b] = k[r120b] * c[sCO2] * c[sOH];
	w[r121f] = k[r121f] * c[sHO2] * c[sCH2O];
	w[r121b] = k[r121b] * c[sH2O2] * c[sHCO];
	w[r122f] = k[r122f] * c[sC] * c[sO2];
	w[r122b] = k[r122b] * c[sCO] * c[sO];
	w[r123f] = k[r123f] * c[sC] * c[sCH2];
	w[r123b] = k[r123b] * c[sC2H] * c[sH];
	w[r124f] = k[r124f] * c[sC] * c[sCH3];
	w[r124b] = k[r124b] * c[sC2H2] * c[sH];
	w[r125f] = k[r125f] * c[sCH] * c[sO2];
	w[r125b] = k[r125b] * c[sHCO] * c[sO];
	w[r126f] = k[r126f] * c[sCH] * c[sH2];
	w[r126b] = k[r126b] * c[sCH2] * c[sH];
	w[r127f] = k[r127f] * c[sCH] * c[sH2O];
	w[r127b] = k[r127b] * c[sCH2O] * c[sH];
	w[r128f] = k[r128f] * c[sCH] * c[sCH2];
	w[r128b] = k[r128b] * c[sC2H2] * c[sH];
	w[r129f] = k[r129f] * c[sCH] * c[sCH3];
	w[r129b] = k[r129b] * c[sC2H3] * c[sH];
	w[r130f] = k[r130f] * c[sCH] * c[sCH4];
	w[r130b] = k[r130b] * c[sC2H4] * c[sH];
	w[r131f] = k[r131f] * c[sCH] * c[sCO];
	w[r131b] = k[r131b] * c[sHCCO];
	w[r132f] = k[r132f] * c[sCH] * c[sCO2];
	w[r132b] = k[r132b] * c[sCO] * c[sHCO];
	w[r133f] = k[r133f] * c[sCH] * c[sCH2O];
	w[r133b] = k[r133b] * c[sCH2CO] * c[sH];
	w[r134f] = k[r134f] * c[sCH] * c[sHCCO];
	w[r134b] = k[r134b] * c[sC2H2] * c[sCO];
	w[r135f] = k[r135f] * c[sCH2] * c[sO2];
	w[r135b] = k[r135b] * c[sHCO] * c[sOH];
	w[r136f] = k[r136f] * c[sCH2] * c[sH2];
	w[r136b] = k[r136b] * c[sCH3] * c[sH];
	w[r137f] = k[r137f] * c[sCH2] * c[sCH2];
	w[r137b] = k[r137b] * c[sC2H2] * c[sH2];
	w[r138f] = k[r138f] * c[sCH2] * c[sCH3];
	w[r138b] = k[r138b] * c[sC2H4] * c[sH];
	w[r139f] = k[r139f] * c[sCH2] * c[sCH4];
	w[r139b] = k[r139b] * c[sCH3] * c[sCH3];
	w[r140f] = k[r140f] * c[sCH2] * c[sCO];
	w[r140b] = k[r140b] * c[sCH2CO];
	w[r141f] = k[r141f] * c[sCH2] * c[sHCCO];
	w[r141b] = k[r141b] * c[sCO] * c[sC2H3];
	w[r142f] = k[r142f] * c[sCH2GSGXCH2] * c[sN2];
	w[r142b] = k[r142b] * c[sN2] * c[sCH2];
	w[r143f] = k[r143f] * c[sCH2GSGXCH2] * c[sAR];
	w[r143b] = k[r143b] * c[sAR] * c[sCH2];
	w[r144f] = k[r144f] * c[sCH2GSGXCH2] * c[sO2];
	w[r144b] = k[r144b] * c[sCO] * c[sOH] * c[sH];
	w[r145f] = k[r145f] * c[sCH2GSGXCH2] * c[sO2];
	w[r145b] = k[r145b] * c[sH2O] * c[sCO];
	w[r146f] = k[r146f] * c[sCH2GSGXCH2] * c[sH2];
	w[r146b] = k[r146b] * c[sH] * c[sCH3];
	w[r147f] = k[r147f] * c[sCH2GSGXCH2] * c[sH2O];
	w[r147b] = k[r147b] * c[sCH3OH];
	w[r148f] = k[r148f] * c[sCH2GSGXCH2] * c[sH2O];
	w[r148b] = k[r148b] * c[sH2O] * c[sCH2];
	w[r149f] = k[r149f] * c[sCH2GSGXCH2] * c[sCH3];
	w[r149b] = k[r149b] * c[sC2H4] * c[sH];
	w[r150f] = k[r150f] * c[sCH2GSGXCH2] * c[sCH4];
	w[r150b] = k[r150b] * c[sCH3] * c[sCH3];
	w[r151f] = k[r151f] * c[sCH2GSGXCH2] * c[sCO];
	w[r151b] = k[r151b] * c[sCO] * c[sCH2];
	w[r152f] = k[r152f] * c[sCH2GSGXCH2] * c[sCO2];
	w[r152b] = k[r152b] * c[sCO2] * c[sCH2];
	w[r153f] = k[r153f] * c[sCH2GSGXCH2] * c[sCO2];
	w[r153b] = k[r153b] * c[sCH2O] * c[sCO];
	w[r154f] = k[r154f] * c[sCH2GSGXCH2] * c[sC2H6];
	w[r154b] = k[r154b] * c[sC2H5] * c[sCH3];
	w[r155f] = k[r155f] * c[sCH3] * c[sO2];
	w[r155b] = k[r155b] * c[sCH3O] * c[sO];
	w[r156f] = k[r156f] * c[sCH3] * c[sO2];
	w[r156b] = k[r156b] * c[sCH2O] * c[sOH];
	w[r157f] = k[r157f] * c[sCH3] * c[sH2O2];
	w[r157b] = k[r157b] * c[sCH4] * c[sHO2];
	w[r158f] = k[r158f] * c[sCH3] * c[sCH3];
	w[r158b] = k[r158b] * c[sC2H6];
	w[r159f] = k[r159f] * c[sCH3] * c[sCH3];
	w[r159b] = k[r159b] * c[sC2H5] * c[sH];
	w[r160f] = k[r160f] * c[sCH3] * c[sHCO];
	w[r160b] = k[r160b] * c[sCO] * c[sCH4];
	w[r161f] = k[r161f] * c[sCH3] * c[sCH2O];
	w[r161b] = k[r161b] * c[sCH4] * c[sHCO];
	w[r162f] = k[r162f] * c[sCH3] * c[sCH3OH];
	w[r162b] = k[r162b] * c[sCH4] * c[sCH2OH];
	w[r163f] = k[r163f] * c[sCH3] * c[sCH3OH];
	w[r163b] = k[r163b] * c[sCH4] * c[sCH3O];
	w[r164f] = k[r164f] * c[sCH3] * c[sC2H4];
	w[r164b] = k[r164b] * c[sCH4] * c[sC2H3];
	w[r165f] = k[r165f] * c[sCH3] * c[sC2H6];
	w[r165b] = k[r165b] * c[sCH4] * c[sC2H5];
	w[r166f] = k[r166f] * c[sHCO] * c[sH2O];
	w[r166b] = k[r166b] * c[sH2O] * c[sCO] * c[sH];
	w[r167f] = k[r167f] * c[sHCO] * M[mM26];
	w[r167b] = k[r167b] * c[sCO] * c[sH] * M[mM26];
	w[r168f] = k[r168f] * c[sHCO] * c[sO2];
	w[r168b] = k[r168b] * c[sCO] * c[sHO2];
	w[r169f] = k[r169f] * c[sCH2OH] * c[sO2];
	w[r169b] = k[r169b] * c[sCH2O] * c[sHO2];
	w[r170f] = k[r170f] * c[sCH3O] * c[sO2];
	w[r170b] = k[r170b] * c[sCH2O] * c[sHO2];
	w[r171f] = k[r171f] * c[sC2H] * c[sO2];
	w[r171b] = k[r171b] * c[sCO] * c[sHCO];
	w[r172f] = k[r172f] * c[sC2H] * c[sH2];
	w[r172b] = k[r172b] * c[sC2H2] * c[sH];
	w[r173f] = k[r173f] * c[sC2H3] * c[sO2];
	w[r173b] = k[r173b] * c[sCH2O] * c[sHCO];
	w[r174f] = k[r174f] * c[sC2H4];
	w[r174b] = k[r174b] * c[sC2H2] * c[sH2];
	w[r175f] = k[r175f] * c[sC2H5] * c[sO2];
	w[r175b] = k[r175b] * c[sC2H4] * c[sHO2];
	w[r176f] = k[r176f] * c[sHCCO] * c[sO2];
	w[r176b] = k[r176b] * c[sCO] * c[sCO] * c[sOH];
	w[r177f] = k[r177f] * c[sHCCO] * c[sHCCO];
	w[r177b] = k[r177b] * c[sC2H2] * c[sCO] * c[sCO];
	w[r178f] = k[r178f] * c[sN] * c[sNO];
	w[r178b] = k[r178b] * c[sO] * c[sN2];
	w[r179f] = k[r179f] * c[sN] * c[sO2];
	w[r179b] = k[r179b] * c[sO] * c[sNO];
	w[r180f] = k[r180f] * c[sN] * c[sOH];
	w[r180b] = k[r180b] * c[sH] * c[sNO];
	w[r181f] = k[r181f] * c[sN2O] * c[sO];
	w[r181b] = k[r181b] * c[sO2] * c[sN2];
	w[r182f] = k[r182f] * c[sN2O] * c[sO];
	w[r182b] = k[r182b] * c[sNO] * c[sNO];
	w[r183f] = k[r183f] * c[sN2O] * c[sH];
	w[r183b] = k[r183b] * c[sOH] * c[sN2];
	w[r184f] = k[r184f] * c[sN2O] * c[sOH];
	w[r184b] = k[r184b] * c[sHO2] * c[sN2];
	w[r185f] = k[r185f] * c[sN2O];
	w[r185b] = k[r185b] * c[sO] * c[sN2];
	w[r186f] = k[r186f] * c[sHO2] * c[sNO];
	w[r186b] = k[r186b] * c[sOH] * c[sNO2];
	w[r187f] = k[r187f] * c[sNO] * c[sO] * M[mM29];
	w[r187b] = k[r187b] * c[sNO2] * M[mM29];
	w[r188f] = k[r188f] * c[sNO2] * c[sO];
	w[r188b] = k[r188b] * c[sO2] * c[sNO];
	w[r189f] = k[r189f] * c[sNO2] * c[sH];
	w[r189b] = k[r189b] * c[sOH] * c[sNO];
	w[r190f] = k[r190f] * c[sNH] * c[sO];
	w[r190b] = k[r190b] * c[sH] * c[sNO];
	w[r191f] = k[r191f] * c[sNH] * c[sH];
	w[r191b] = k[r191b] * c[sH2] * c[sN];
	w[r192f] = k[r192f] * c[sNH] * c[sOH];
	w[r192b] = k[r192b] * c[sH] * c[sHNO];
	w[r193f] = k[r193f] * c[sNH] * c[sOH];
	w[r193b] = k[r193b] * c[sH2O] * c[sN];
	w[r194f] = k[r194f] * c[sNH] * c[sO2];
	w[r194b] = k[r194b] * c[sO] * c[sHNO];
	w[r195f] = k[r195f] * c[sNH] * c[sO2];
	w[r195b] = k[r195b] * c[sOH] * c[sNO];
	w[r196f] = k[r196f] * c[sNH] * c[sN];
	w[r196b] = k[r196b] * c[sH] * c[sN2];
	w[r197f] = k[r197f] * c[sNH] * c[sH2O];
	w[r197b] = k[r197b] * c[sH2] * c[sHNO];
	w[r198f] = k[r198f] * c[sNH] * c[sNO];
	w[r198b] = k[r198b] * c[sOH] * c[sN2];
	w[r199f] = k[r199f] * c[sNH] * c[sNO];
	w[r199b] = k[r199b] * c[sH] * c[sN2O];
	w[r200f] = k[r200f] * c[sNH2] * c[sO];
	w[r200b] = k[r200b] * c[sNH] * c[sOH];
	w[r201f] = k[r201f] * c[sNH2] * c[sO];
	w[r201b] = k[r201b] * c[sHNO] * c[sH];
	w[r202f] = k[r202f] * c[sNH2] * c[sH];
	w[r202b] = k[r202b] * c[sH2] * c[sNH];
	w[r203f] = k[r203f] * c[sNH2] * c[sOH];
	w[r203b] = k[r203b] * c[sH2O] * c[sNH];
	w[r204f] = k[r204f] * c[sNNH];
	w[r204b] = k[r204b] * c[sH] * c[sN2];
	w[r205f] = k[r205f] * c[sNNH] * M[mM30];
	w[r205b] = k[r205b] * c[sH] * c[sN2] * M[mM30];
	w[r206f] = k[r206f] * c[sNNH] * c[sO2];
	w[r206b] = k[r206b] * c[sN2] * c[sHO2];
	w[r207f] = k[r207f] * c[sNNH] * c[sO];
	w[r207b] = k[r207b] * c[sN2] * c[sOH];
	w[r208f] = k[r208f] * c[sNNH] * c[sO];
	w[r208b] = k[r208b] * c[sNO] * c[sNH];
	w[r209f] = k[r209f] * c[sNNH] * c[sH];
	w[r209b] = k[r209b] * c[sN2] * c[sH2];
	w[r210f] = k[r210f] * c[sNNH] * c[sOH];
	w[r210b] = k[r210b] * c[sN2] * c[sH2O];
	w[r211f] = k[r211f] * c[sNNH] * c[sCH3];
	w[r211b] = k[r211b] * c[sN2] * c[sCH4];
	w[r212f] = k[r212f] * c[sH] * c[sNO] * M[mM31];
	w[r212b] = k[r212b] * c[sHNO] * M[mM31];
	w[r213f] = k[r213f] * c[sHNO] * c[sO];
	w[r213b] = k[r213b] * c[sOH] * c[sNO];
	w[r214f] = k[r214f] * c[sHNO] * c[sH];
	w[r214b] = k[r214b] * c[sNO] * c[sH2];
	w[r215f] = k[r215f] * c[sHNO] * c[sOH];
	w[r215b] = k[r215b] * c[sH2O] * c[sNO];
	w[r216f] = k[r216f] * c[sHNO] * c[sO2];
	w[r216b] = k[r216b] * c[sNO] * c[sHO2];
	w[r217f] = k[r217f] * c[sCN] * c[sO];
	w[r217b] = k[r217b] * c[sN] * c[sCO];
	w[r218f] = k[r218f] * c[sCN] * c[sOH];
	w[r218b] = k[r218b] * c[sH] * c[sNCO];
	w[r219f] = k[r219f] * c[sCN] * c[sH2O];
	w[r219b] = k[r219b] * c[sOH] * c[sHCN];
	w[r220f] = k[r220f] * c[sCN] * c[sO2];
	w[r220b] = k[r220b] * c[sO] * c[sNCO];
	w[r221f] = k[r221f] * c[sCN] * c[sH2];
	w[r221b] = k[r221b] * c[sH] * c[sHCN];
	w[r222f] = k[r222f] * c[sNCO] * c[sO];
	w[r222b] = k[r222b] * c[sCO] * c[sNO];
	w[r223f] = k[r223f] * c[sNCO] * c[sH];
	w[r223b] = k[r223b] * c[sCO] * c[sNH];
	w[r224f] = k[r224f] * c[sNCO] * c[sOH];
	w[r224b] = k[r224b] * c[sCO] * c[sH] * c[sNO];
	w[r225f] = k[r225f] * c[sNCO] * c[sN];
	w[r225b] = k[r225b] * c[sCO] * c[sN2];
	w[r226f] = k[r226f] * c[sNCO] * c[sO2];
	w[r226b] = k[r226b] * c[sCO2] * c[sNO];
	w[r227f] = k[r227f] * c[sNCO] * M[mM32];
	w[r227b] = k[r227b] * c[sCO] * c[sN] * M[mM32];
	w[r228f] = k[r228f] * c[sNCO] * c[sNO];
	w[r228b] = k[r228b] * c[sCO] * c[sN2O];
	w[r229f] = k[r229f] * c[sNCO] * c[sNO];
	w[r229b] = k[r229b] * c[sCO2] * c[sN2];
	w[r230f] = k[r230f] * c[sHCN] * M[mM33];
	w[r230b] = k[r230b] * c[sCN] * c[sH] * M[mM33];
	w[r231f] = k[r231f] * c[sHCN] * c[sO];
	w[r231b] = k[r231b] * c[sH] * c[sNCO];
	w[r232f] = k[r232f] * c[sHCN] * c[sO];
	w[r232b] = k[r232b] * c[sCO] * c[sNH];
	w[r233f] = k[r233f] * c[sHCN] * c[sO];
	w[r233b] = k[r233b] * c[sOH] * c[sCN];
	w[r234f] = k[r234f] * c[sHCN] * c[sOH];
	w[r234b] = k[r234b] * c[sH] * c[sHOCN];
	w[r235f] = k[r235f] * c[sHCN] * c[sOH];
	w[r235b] = k[r235b] * c[sH] * c[sHNCO];
	w[r236f] = k[r236f] * c[sHCN] * c[sOH];
	w[r236b] = k[r236b] * c[sCO] * c[sNH2];
	w[r237f] = k[r237f] * c[sH] * c[sHCN] * M[mM34];
	w[r237b] = k[r237b] * c[sH2CN] * M[mM34];
	w[r238f] = k[r238f] * c[sH2CN] * c[sN];
	w[r238b] = k[r238b] * c[sCH2] * c[sN2];
	w[r239f] = k[r239f] * c[sC] * c[sN2];
	w[r239b] = k[r239b] * c[sN] * c[sCN];
	w[r240f] = k[r240f] * c[sCH] * c[sN2];
	w[r240b] = k[r240b] * c[sN] * c[sHCN];
	w[r241f] = k[r241f] * c[sCH] * c[sN2];
	w[r241b] = k[r241b] * c[sHCNN];
	w[r242f] = k[r242f] * c[sCH2] * c[sN2];
	w[r242b] = k[r242b] * c[sNH] * c[sHCN];
	w[r243f] = k[r243f] * c[sCH2GSGXCH2] * c[sN2];
	w[r243b] = k[r243b] * c[sHCN] * c[sNH];
	w[r244f] = k[r244f] * c[sC] * c[sNO];
	w[r244b] = k[r244b] * c[sO] * c[sCN];
	w[r245f] = k[r245f] * c[sC] * c[sNO];
	w[r245b] = k[r245b] * c[sN] * c[sCO];
	w[r246f] = k[r246f] * c[sCH] * c[sNO];
	w[r246b] = k[r246b] * c[sO] * c[sHCN];
	w[r247f] = k[r247f] * c[sCH] * c[sNO];
	w[r247b] = k[r247b] * c[sNCO] * c[sH];
	w[r248f] = k[r248f] * c[sCH] * c[sNO];
	w[r248b] = k[r248b] * c[sHCO] * c[sN];
	w[r249f] = k[r249f] * c[sCH2] * c[sNO];
	w[r249b] = k[r249b] * c[sHNCO] * c[sH];
	w[r250f] = k[r250f] * c[sCH2] * c[sNO];
	w[r250b] = k[r250b] * c[sHCN] * c[sOH];
	w[r251f] = k[r251f] * c[sCH2] * c[sNO];
	w[r251b] = k[r251b] * c[sHCNO] * c[sH];
	w[r252f] = k[r252f] * c[sCH2GSGXCH2] * c[sNO];
	w[r252b] = k[r252b] * c[sHNCO] * c[sH];
	w[r253f] = k[r253f] * c[sCH2GSGXCH2] * c[sNO];
	w[r253b] = k[r253b] * c[sHCN] * c[sOH];
	w[r254f] = k[r254f] * c[sCH2GSGXCH2] * c[sNO];
	w[r254b] = k[r254b] * c[sHCNO] * c[sH];
	w[r255f] = k[r255f] * c[sCH3] * c[sNO];
	w[r255b] = k[r255b] * c[sH2O] * c[sHCN];
	w[r256f] = k[r256f] * c[sCH3] * c[sNO];
	w[r256b] = k[r256b] * c[sOH] * c[sH2CN];
	w[r257f] = k[r257f] * c[sHCNN] * c[sO];
	w[r257b] = k[r257b] * c[sN2] * c[sH] * c[sCO];
	w[r258f] = k[r258f] * c[sHCNN] * c[sO];
	w[r258b] = k[r258b] * c[sNO] * c[sHCN];
	w[r259f] = k[r259f] * c[sHCNN] * c[sO2];
	w[r259b] = k[r259b] * c[sN2] * c[sHCO] * c[sO];
	w[r260f] = k[r260f] * c[sHCNN] * c[sOH];
	w[r260b] = k[r260b] * c[sN2] * c[sHCO] * c[sH];
	w[r261f] = k[r261f] * c[sHCNN] * c[sH];
	w[r261b] = k[r261b] * c[sN2] * c[sCH2];
	w[r262f] = k[r262f] * c[sHNCO] * c[sO];
	w[r262b] = k[r262b] * c[sCO2] * c[sNH];
	w[r263f] = k[r263f] * c[sHNCO] * c[sO];
	w[r263b] = k[r263b] * c[sCO] * c[sHNO];
	w[r264f] = k[r264f] * c[sHNCO] * c[sO];
	w[r264b] = k[r264b] * c[sOH] * c[sNCO];
	w[r265f] = k[r265f] * c[sHNCO] * c[sH];
	w[r265b] = k[r265b] * c[sCO] * c[sNH2];
	w[r266f] = k[r266f] * c[sHNCO] * c[sH];
	w[r266b] = k[r266b] * c[sNCO] * c[sH2];
	w[r267f] = k[r267f] * c[sHNCO] * c[sOH];
	w[r267b] = k[r267b] * c[sH2O] * c[sNCO];
	w[r268f] = k[r268f] * c[sHNCO] * c[sOH];
	w[r268b] = k[r268b] * c[sCO2] * c[sNH2];
	w[r269f] = k[r269f] * c[sHNCO] * M[mM36];
	w[r269b] = k[r269b] * c[sCO] * c[sNH] * M[mM36];
	w[r270f] = k[r270f] * c[sHCNO] * c[sH];
	w[r270b] = k[r270b] * c[sHNCO] * c[sH];
	w[r271f] = k[r271f] * c[sHCNO] * c[sH];
	w[r271b] = k[r271b] * c[sHCN] * c[sOH];
	w[r272f] = k[r272f] * c[sHCNO] * c[sH];
	w[r272b] = k[r272b] * c[sCO] * c[sNH2];
	w[r273f] = k[r273f] * c[sHOCN] * c[sH];
	w[r273b] = k[r273b] * c[sHNCO] * c[sH];
	w[r274f] = k[r274f] * c[sHCCO] * c[sNO];
	w[r274b] = k[r274b] * c[sCO] * c[sHCNO];
	w[r275f] = k[r275f] * c[sCH3] * c[sN];
	w[r275b] = k[r275b] * c[sH] * c[sH2CN];
	w[r276f] = k[r276f] * c[sCH3] * c[sN];
	w[r276b] = k[r276b] * c[sH2] * c[sHCN];
	w[r277f] = k[r277f] * c[sNH3] * c[sH];
	w[r277b] = k[r277b] * c[sH2] * c[sNH2];
	w[r278f] = k[r278f] * c[sNH3] * c[sOH];
	w[r278b] = k[r278b] * c[sH2O] * c[sNH2];
	w[r279f] = k[r279f] * c[sNH3] * c[sO];
	w[r279b] = k[r279b] * c[sOH] * c[sNH2];
	w[r280] = k[r280] * c[sC3H3] * c[sO2];
	w[r281] = k[r281] * c[sC3H3] * c[sO];
	w[r282f] = k[r282f] * c[sC3H4];
	w[r282b] = k[r282b] * c[sC3H3] * c[sH];
	w[r283] = k[r283] * c[sC3H4] * c[sO];
	w[r284] = k[r284] * c[sC3H4] * c[sO];
	w[r285] = k[r285] * c[sC3H4] * c[sOH];
	w[r286] = k[r286] * c[sC3H4] * c[sOH];
	w[r287f] = k[r287f] * c[sC3H5];
	w[r287b] = k[r287b] * c[sC3H4] * c[sH];
	w[r288] = k[r288] * c[sC3H5] * c[sH];
	w[r289f] = k[r289f] * c[sC3H6];
	w[r289b] = k[r289b] * c[sC2H3] * c[sCH3];
	w[r290] = k[r290] * c[sC3H6] * c[sH];
	w[r291] = k[r291] * c[sNXC3H7];
	w[r292f] = k[r292f] * c[sNXC3H7];
	w[r292b] = k[r292b] * c[sC3H6] * c[sH];
	w[r293] = k[r293] * c[sIXC3H7];
	w[r294] = k[r294] * c[sIXC3H7] * c[sO2];
	w[r295] = k[r295] * c[sC3H8] * c[sH];
	w[r296] = k[r296] * c[sC3H8] * c[sH];
	w[r297] = k[r297] * c[sC3H8] * c[sO];
	w[r298] = k[r298] * c[sC3H8] * c[sO];
	w[r299] = k[r299] * c[sC3H8] * c[sOH];
	w[r300] = k[r300] * c[sC3H8] * c[sOH];


	cdot[sN2] = - w[r36f] + w[r36f] - w[r36b] + w[r36b]
		- w[r142f] + w[r142f] - w[r142b] + w[r142b]
		+ w[r178f] - w[r178b] + w[r181f] - w[r181b]
		+ w[r183f] - w[r183b] + w[r184f] - w[r184b]
		+ w[r185f] - w[r185b] + w[r196f] - w[r196b]
		+ w[r198f] - w[r198b] + w[r204f] - w[r204b]
		+ w[r205f] - w[r205b] + w[r206f] - w[r206b]
		+ w[r207f] - w[r207b] + w[r209f] - w[r209b]
		+ w[r210f] - w[r210b] + w[r211f] - w[r211b]
		+ w[r225f] - w[r225b] + w[r229f] - w[r229b]
		+ w[r238f] - w[r238b] - w[r239f] + w[r239b]
		- w[r240f] + w[r240b] - w[r241f] + w[r241b]
		- w[r242f] + w[r242b] - w[r243f] + w[r243b]
		+ w[r257f] - w[r257b] + w[r259f] - w[r259b]
		+ w[r260f] - w[r260b] + w[r261f] - w[r261b];

	cdot[sAR] = - w[r37f] + w[r37f] - w[r37b] + w[r37b]
		- w[r143f] + w[r143f] - w[r143b] + w[r143b];

	cdot[sO] = - 2 * w[r1f] + 2 * w[r1b] - w[r2f] + w[r2b]
		- w[r3f] + w[r3b] - w[r4f] + w[r4b]
		- w[r5f] + w[r5b] - w[r6f] + w[r6b]
		- w[r7f] + w[r7b] - w[r8f] + w[r8b]
		- w[r9f] + w[r9b] - w[r10f] + w[r10b]
		- w[r11f] + w[r11b] - w[r12f] + w[r12b]
		- w[r13f] + w[r13b] - w[r14f] + w[r14b]
		- w[r15f] + w[r15b] - w[r16f] + w[r16b]
		- w[r17f] + w[r17b] - w[r18f] + w[r18b]
		- w[r19f] + w[r19b] - w[r20f] + w[r20b]
		- w[r21f] + w[r21b] - w[r22f] + w[r22b]
		- w[r23f] + w[r23b] - w[r24f] + w[r24b]
		- w[r25f] + w[r25b] - w[r26f] + w[r26b]
		- w[r27f] + w[r27b] - w[r28f] + w[r28b]
		- w[r29f] + w[r29b] - w[r30f] + w[r30b]
		+ w[r31f] - w[r31b] + w[r38f] - w[r38b]
		+ w[r44f] - w[r44b] + w[r86f] - w[r86b]
		+ w[r122f] - w[r122b] + w[r125f] - w[r125b]
		+ w[r155f] - w[r155b] + w[r178f] - w[r178b]
		+ w[r179f] - w[r179b] - w[r181f] + w[r181b]
		- w[r182f] + w[r182b] + w[r185f] - w[r185b]
		- w[r187f] + w[r187b] - w[r188f] + w[r188b]
		- w[r190f] + w[r190b] + w[r194f] - w[r194b]
		- w[r200f] + w[r200b] - w[r201f] + w[r201b]
		- w[r207f] + w[r207b] - w[r208f] + w[r208b]
		- w[r213f] + w[r213b] - w[r217f] + w[r217b]
		+ w[r220f] - w[r220b] - w[r222f] + w[r222b]
		- w[r231f] + w[r231b] - w[r232f] + w[r232b]
		- w[r233f] + w[r233b] + w[r244f] - w[r244b]
		+ w[r246f] - w[r246b] - w[r257f] + w[r257b]
		- w[r258f] + w[r258b] + w[r259f] - w[r259b]
		- w[r262f] + w[r262b] - w[r263f] + w[r263b]
		- w[r264f] + w[r264b] - w[r279f] + w[r279b]
		- w[r281] - w[r283] - w[r284] - w[r297]
		- w[r298];

	cdot[sO2] = w[r1f] - w[r1b] + w[r4f] - w[r4b]
		- w[r31f] + w[r31b] - w[r32f] + w[r32b]
		- w[r33f] + w[r33b] - 2 * w[r34f] + w[r34f]
		- w[r34b] + 2 * w[r34b] - w[r35f] + w[r35b]
		- w[r36f] + w[r36b] - w[r37f] + w[r37b]
		- w[r38f] + w[r38b] + w[r45f] - w[r45b]
		+ w[r87f] - w[r87b] + w[r115f] - w[r115b]
		+ w[r116f] - w[r116b] + w[r118f] - w[r118b]
		- w[r122f] + w[r122b] - w[r125f] + w[r125b]
		- w[r135f] + w[r135b] - w[r144f] + w[r144b]
		- w[r145f] + w[r145b] - w[r155f] + w[r155b]
		- w[r156f] + w[r156b] - w[r168f] + w[r168b]
		- w[r169f] + w[r169b] - w[r170f] + w[r170b]
		- w[r171f] + w[r171b] - w[r173f] + w[r173b]
		- w[r175f] + w[r175b] - w[r176f] + w[r176b]
		- w[r179f] + w[r179b] + w[r181f] - w[r181b]
		+ w[r188f] - w[r188b] - w[r194f] + w[r194b]
		- w[r195f] + w[r195b] - w[r206f] + w[r206b]
		- w[r216f] + w[r216b] - w[r220f] + w[r220b]
		- w[r226f] + w[r226b] - w[r259f] + w[r259b]
		- w[r280] - w[r294];

	cdot[sH] = - w[r2f] + w[r2b] + w[r3f] - w[r3b]
		+ w[r6f] - w[r6b] + w[r7f] - w[r7b]
		+ w[r9f] - w[r9b] + w[r10f] - w[r10b]
		+ w[r14f] - w[r14b] + w[r21f] - w[r21b]
		+ w[r24f] - w[r24b] + w[r28f] - w[r28b]
		- w[r33f] + w[r33b] - w[r34f] + w[r34b]
		- w[r35f] + w[r35b] - w[r36f] + w[r36b]
		- w[r37f] + w[r37b] - w[r38f] + w[r38b]
		- 2 * w[r39f] + 2 * w[r39b] - 2 * w[r40f] + 2 * w[r40b]
		- 2 * w[r41f] + 2 * w[r41b] - 2 * w[r42f] + 2 * w[r42b]
		- w[r43f] + w[r43b] - w[r44f] + w[r44b]
		- w[r45f] + w[r45b] - w[r46f] + w[r46b]
		- w[r47f] + w[r47b] - w[r48f] + w[r48b]
		- w[r49f] + w[r49b] - w[r50f] + w[r50b]
		- w[r51f] + w[r51b] - w[r52f] + w[r52b]
		- w[r53f] + w[r53b] - w[r54f] + w[r54b]
		- w[r55f] + w[r55b] - w[r56f] + w[r56b]
		- w[r57f] + w[r57b] - w[r58f] + w[r58b]
		- w[r59f] + w[r59b] - w[r60f] + w[r60b]
		- w[r61f] + w[r61b] - w[r62f] + w[r62b]
		- w[r63f] + w[r63b] - w[r64f] + w[r64f]
		- w[r64b] + w[r64b] - w[r65f] + w[r65b]
		- w[r66f] + w[r66b] - w[r67f] + w[r67b]
		- w[r68f] + w[r68b] - w[r69f] + w[r69b]
		- w[r70f] + w[r70b] - w[r71f] + w[r71b]
		- w[r72f] + w[r72b] - w[r73f] + w[r73b]
		- w[r74f] + w[r74b] - w[r75f] + w[r75b]
		- w[r76f] + w[r76b] - w[r77f] + w[r77b]
		- w[r78f] + w[r78b] - w[r79f] + w[r79b]
		- w[r80f] + w[r80b] - w[r81f] + w[r81b]
		- w[r82f] + w[r82f] - w[r82b] + w[r82b]
		+ w[r84f] - w[r84b] + w[r90f] - w[r90b]
		+ w[r91f] - w[r91b] + w[r92f] - w[r92b]
		+ w[r94f] - w[r94b] + w[r99f] - w[r99b]
		+ w[r106f] - w[r106b] + w[r107f] - w[r107b]
		+ w[r108f] - w[r108b] + w[r123f] - w[r123b]
		+ w[r124f] - w[r124b] + w[r126f] - w[r126b]
		+ w[r127f] - w[r127b] + w[r128f] - w[r128b]
		+ w[r129f] - w[r129b] + w[r130f] - w[r130b]
		+ w[r133f] - w[r133b] + w[r136f] - w[r136b]
		+ w[r138f] - w[r138b] + w[r144f] - w[r144b]
		+ w[r146f] - w[r146b] + w[r149f] - w[r149b]
		+ w[r159f] - w[r159b] + w[r166f] - w[r166b]
		+ w[r167f] - w[r167b] + w[r172f] - w[r172b]
		+ w[r180f] - w[r180b] - w[r183f] + w[r183b]
		- w[r189f] + w[r189b] + w[r190f] - w[r190b]
		- w[r191f] + w[r191b] + w[r192f] - w[r192b]
		+ w[r196f] - w[r196b] + w[r199f] - w[r199b]
		+ w[r201f] - w[r201b] - w[r202f] + w[r202b]
		+ w[r204f] - w[r204b] + w[r205f] - w[r205b]
		- w[r209f] + w[r209b] - w[r212f] + w[r212b]
		- w[r214f] + w[r214b] + w[r218f] - w[r218b]
		+ w[r221f] - w[r221b] - w[r223f] + w[r223b]
		+ w[r224f] - w[r224b] + w[r230f] - w[r230b]
		+ w[r231f] - w[r231b] + w[r234f] - w[r234b]
		+ w[r235f] - w[r235b] - w[r237f] + w[r237b]
		+ w[r247f] - w[r247b] + w[r249f] - w[r249b]
		+ w[r251f] - w[r251b] + w[r252f] - w[r252b]
		+ w[r254f] - w[r254b] + w[r257f] - w[r257b]
		+ w[r260f] - w[r260b] - w[r261f] + w[r261b]
		- w[r265f] + w[r265b] - w[r266f] + w[r266b]
		- w[r270f] + w[r270f] - w[r270b] + w[r270b]
		- w[r271f] + w[r271b] - w[r272f] + w[r272b]
		- w[r273f] + w[r273f] - w[r273b] + w[r273b]
		+ w[r275f] - w[r275b] - w[r277f] + w[r277b]
		+ w[r282f] - w[r282b] + w[r287f] - w[r287b]
		- w[r288] - w[r290] + w[r292f] - w[r292b]
		- w[r295] - w[r296];

	cdot[sOH] = w[r2f] - w[r2b] + w[r3f] - w[r3b]
		+ w[r4f] - w[r4b] + w[r5f] - w[r5b]
		+ w[r11f] - w[r11b] + w[r13f] - w[r13b]
		+ w[r15f] - w[r15b] + w[r16f] - w[r16b]
		+ w[r17f] - w[r17b] + w[r18f] - w[r18b]
		+ w[r19f] - w[r19b] + w[r22f] - w[r22b]
		+ w[r27f] - w[r27b] + w[r29f] - w[r29b]
		+ w[r38f] - w[r38b] - w[r43f] + w[r43b]
		+ 2 * w[r46f] - 2 * w[r46b] + w[r48f] - w[r48b]
		+ w[r61f] - w[r61b] + w[r66f] - w[r66b]
		- w[r84f] + w[r84b] - 2 * w[r85f] + 2 * w[r85b]
		- 2 * w[r86f] + 2 * w[r86b] - w[r87f] + w[r87b]
		- w[r88f] + w[r88b] - w[r89f] + w[r89b]
		- w[r90f] + w[r90b] - w[r91f] + w[r91b]
		- w[r92f] + w[r92b] - w[r93f] + w[r93b]
		- w[r94f] + w[r94b] - w[r95f] + w[r95b]
		- w[r96f] + w[r96b] - w[r97f] + w[r97b]
		- w[r98f] + w[r98b] - w[r99f] + w[r99b]
		- w[r100f] + w[r100b] - w[r101f] + w[r101b]
		- w[r102f] + w[r102b] - w[r103f] + w[r103b]
		- w[r104f] + w[r104b] - w[r105f] + w[r105b]
		- w[r106f] + w[r106b] - w[r107f] + w[r107b]
		- w[r108f] + w[r108b] - w[r109f] + w[r109b]
		- w[r110f] + w[r110b] - w[r111f] + w[r111b]
		- w[r112f] + w[r112b] - w[r113f] + w[r113b]
		- w[r114f] + w[r114b] + w[r117f] - w[r117b]
		+ w[r119f] - w[r119b] + w[r120f] - w[r120b]
		+ w[r135f] - w[r135b] + w[r144f] - w[r144b]
		+ w[r156f] - w[r156b] + w[r176f] - w[r176b]
		- w[r180f] + w[r180b] + w[r183f] - w[r183b]
		- w[r184f] + w[r184b] + w[r186f] - w[r186b]
		+ w[r189f] - w[r189b] - w[r192f] + w[r192b]
		- w[r193f] + w[r193b] + w[r195f] - w[r195b]
		+ w[r198f] - w[r198b] + w[r200f] - w[r200b]
		- w[r203f] + w[r203b] + w[r207f] - w[r207b]
		- w[r210f] + w[r210b] + w[r213f] - w[r213b]
		- w[r215f] + w[r215b] - w[r218f] + w[r218b]
		+ w[r219f] - w[r219b] - w[r224f] + w[r224b]
		+ w[r233f] - w[r233b] - w[r234f] + w[r234b]
		- w[r235f] + w[r235b] - w[r236f] + w[r236b]
		+ w[r250f] - w[r250b] + w[r253f] - w[r253b]
		+ w[r256f] - w[r256b] - w[r260f] + w[r260b]
		+ w[r264f] - w[r264b] - w[r267f] + w[r267b]
		- w[r268f] + w[r268b] + w[r271f] - w[r271b]
		- w[r278f] + w[r278b] + w[r279f] - w[r279b]
		- w[r285] - w[r286] + w[r297] + w[r298]
		- w[r299] - w[r300];

	cdot[sH2] = - w[r3f] + w[r3b] + w[r8f] - w[r8b]
		+ w[r39f] - w[r39b] - w[r40f] + 2 * w[r40f]
		- 2 * w[r40b] + w[r40b] + w[r41f] - w[r41b]
		+ w[r42f] - w[r42b] + w[r45f] - w[r45b]
		+ w[r47f] - w[r47b] + w[r49f] - w[r49b]
		+ w[r51f] - w[r51b] + w[r53f] - w[r53b]
		+ w[r55f] - w[r55b] + w[r58f] - w[r58b]
		+ w[r60f] - w[r60b] + w[r65f] - w[r65b]
		+ w[r68f] - w[r68b] + w[r69f] - w[r69b]
		+ w[r73f] - w[r73b] + w[r75f] - w[r75b]
		+ w[r77f] - w[r77b] + w[r78f] - w[r78b]
		+ w[r80f] - w[r80b] - w[r83f] + w[r83b]
		- w[r84f] + w[r84b] - w[r126f] + w[r126b]
		- w[r136f] + w[r136b] + w[r137f] - w[r137b]
		- w[r146f] + w[r146b] - w[r172f] + w[r172b]
		+ w[r174f] - w[r174b] + w[r191f] - w[r191b]
		+ w[r197f] - w[r197b] + w[r202f] - w[r202b]
		+ w[r209f] - w[r209b] + w[r214f] - w[r214b]
		- w[r221f] + w[r221b] + w[r266f] - w[r266b]
		+ w[r276f] - w[r276b] + w[r277f] - w[r277b]
		+ w[r288] + w[r290] + w[r295] + w[r296];

	cdot[sHO2] = - w[r4f] + w[r4b] + w[r5f] - w[r5b]
		+ w[r32f] - w[r32b] + w[r33f] - w[r33b]
		+ w[r34f] - w[r34b] + w[r35f] - w[r35b]
		+ w[r36f] - w[r36b] + w[r37f] - w[r37b]
		- w[r44f] + w[r44b] - w[r45f] + w[r45b]
		- w[r46f] + w[r46b] + w[r47f] - w[r47b]
		- w[r87f] + w[r87b] + w[r88f] - w[r88b]
		+ w[r89f] - w[r89b] - 2 * w[r115f] + 2 * w[r115b]
		- 2 * w[r116f] + 2 * w[r116b] - w[r117f] + w[r117b]
		- w[r118f] + w[r118b] - w[r119f] + w[r119b]
		- w[r120f] + w[r120b] - w[r121f] + w[r121b]
		+ w[r157f] - w[r157b] + w[r168f] - w[r168b]
		+ w[r169f] - w[r169b] + w[r170f] - w[r170b]
		+ w[r175f] - w[r175b] + w[r184f] - w[r184b]
		- w[r186f] + w[r186b] + w[r206f] - w[r206b]
		+ w[r216f] - w[r216b] + w[r294];

	cdot[sH2O2] = - w[r5f] + w[r5b] - w[r47f] + w[r47b]
		- w[r48f] + w[r48b] + w[r85f] - w[r85b]
		- w[r88f] + w[r88b] - w[r89f] + w[r89b]
		+ w[r115f] - w[r115b] + w[r116f] - w[r116b]
		+ w[r121f] - w[r121b] - w[r157f] + w[r157b];

	cdot[sCH] = - w[r6f] + w[r6b] + w[r20f] - w[r20b]
		- w[r49f] + w[r49b] + w[r51f] - w[r51b]
		- w[r91f] + w[r91b] + w[r93f] - w[r93b]
		- w[r125f] + w[r125b] - w[r126f] + w[r126b]
		- w[r127f] + w[r127b] - w[r128f] + w[r128b]
		- w[r129f] + w[r129b] - w[r130f] + w[r130b]
		- w[r131f] + w[r131b] - w[r132f] + w[r132b]
		- w[r133f] + w[r133b] - w[r134f] + w[r134b]
		- w[r240f] + w[r240b] - w[r241f] + w[r241b]
		- w[r246f] + w[r246b] - w[r247f] + w[r247b]
		- w[r248f] + w[r248b];

	cdot[sCO] = w[r6f] - w[r6b] + w[r8f] - w[r8b]
		- w[r12f] + w[r12b] + w[r13f] - w[r13b]
		+ w[r20f] - w[r20b] + w[r23f] - w[r23b]
		+ 2 * w[r28f] - 2 * w[r28b] - w[r31f] + w[r31b]
		+ w[r55f] - w[r55b] + w[r79f] - w[r79b]
		+ w[r81f] - w[r81b] - w[r83f] + w[r83b]
		+ w[r90f] - w[r90b] - w[r99f] + w[r99b]
		+ w[r100f] - w[r100b] + w[r110f] - w[r110b]
		- w[r120f] + w[r120b] + w[r122f] - w[r122b]
		- w[r131f] + w[r131b] + w[r132f] - w[r132b]
		+ w[r134f] - w[r134b] - w[r140f] + w[r140b]
		+ w[r141f] - w[r141b] + w[r144f] - w[r144b]
		+ w[r145f] - w[r145b] - w[r151f] + w[r151f]
		- w[r151b] + w[r151b] + w[r153f] - w[r153b]
		+ w[r160f] - w[r160b] + w[r166f] - w[r166b]
		+ w[r167f] - w[r167b] + w[r168f] - w[r168b]
		+ w[r171f] - w[r171b] + 2 * w[r176f] - 2 * w[r176b]
		+ 2 * w[r177f] - 2 * w[r177b] + w[r217f] - w[r217b]
		+ w[r222f] - w[r222b] + w[r223f] - w[r223b]
		+ w[r224f] - w[r224b] + w[r225f] - w[r225b]
		+ w[r227f] - w[r227b] + w[r228f] - w[r228b]
		+ w[r232f] - w[r232b] + w[r236f] - w[r236b]
		+ w[r245f] - w[r245b] + w[r257f] - w[r257b]
		+ w[r263f] - w[r263b] + w[r265f] - w[r265b]
		+ w[r269f] - w[r269b] + w[r272f] - w[r272b]
		+ w[r274f] - w[r274b] + w[r281];

	cdot[sCH2] = - w[r7f] + w[r7b] + w[r23f] - w[r23b]
		+ w[r30f] - w[r30b] - w[r50f] + w[r50b]
		- w[r92f] + w[r92b] - w[r93f] + w[r93b]
		+ w[r96f] - w[r96b] - w[r117f] + w[r117b]
		- w[r123f] + w[r123b] + w[r126f] - w[r126b]
		- w[r128f] + w[r128b] - w[r135f] + w[r135b]
		- w[r136f] + w[r136b] - 2 * w[r137f] + 2 * w[r137b]
		- w[r138f] + w[r138b] - w[r139f] + w[r139b]
		- w[r140f] + w[r140b] - w[r141f] + w[r141b]
		+ w[r142f] - w[r142b] + w[r143f] - w[r143b]
		+ w[r148f] - w[r148b] + w[r151f] - w[r151b]
		+ w[r152f] - w[r152b] + w[r238f] - w[r238b]
		- w[r242f] + w[r242b] - w[r249f] + w[r249b]
		- w[r250f] + w[r250b] - w[r251f] + w[r251b]
		+ w[r261f] - w[r261b];

	cdot[sHCO] = w[r7f] - w[r7b] + w[r9f] - w[r9b]
		- w[r13f] + w[r13b] - w[r14f] + w[r14b]
		+ w[r15f] - w[r15b] + w[r25f] - w[r25b]
		+ w[r32f] - w[r32b] - w[r54f] + w[r54b]
		- w[r55f] + w[r55b] + w[r58f] - w[r58b]
		+ w[r91f] - w[r91b] - w[r100f] + w[r100b]
		+ w[r101f] - w[r101b] + w[r121f] - w[r121b]
		+ w[r125f] - w[r125b] + w[r132f] - w[r132b]
		+ w[r135f] - w[r135b] - w[r160f] + w[r160b]
		+ w[r161f] - w[r161b] - w[r166f] + w[r166b]
		- w[r167f] + w[r167b] - w[r168f] + w[r168b]
		+ w[r171f] - w[r171b] + w[r173f] - w[r173b]
		+ w[r248f] - w[r248b] + w[r259f] - w[r259b]
		+ w[r260f] - w[r260b] + w[r284] + w[r286];

	cdot[sCH2GSGXCH2] = - w[r8f] + w[r8b] - w[r9f] + w[r9b]
		- w[r51f] + w[r51b] + w[r62f] - w[r62b]
		+ w[r67f] - w[r67b] + w[r79f] - w[r79b]
		- w[r94f] + w[r94b] + w[r97f] - w[r97b]
		- w[r142f] + w[r142b] - w[r143f] + w[r143b]
		- w[r144f] + w[r144b] - w[r145f] + w[r145b]
		- w[r146f] + w[r146b] - w[r147f] + w[r147b]
		- w[r148f] + w[r148b] - w[r149f] + w[r149b]
		- w[r150f] + w[r150b] - w[r151f] + w[r151b]
		- w[r152f] + w[r152b] - w[r153f] + w[r153b]
		- w[r154f] + w[r154b] - w[r243f] + w[r243b]
		- w[r252f] + w[r252b] - w[r253f] + w[r253b]
		- w[r254f] + w[r254b];

	cdot[sCH3] = - w[r10f] + w[r10b] + w[r11f] - w[r11b]
		+ w[r25f] - w[r25b] + w[r26f] - w[r26b]
		+ w[r50f] - w[r50b] - w[r52f] + w[r52b]
		+ w[r53f] - w[r53b] + w[r61f] - w[r61b]
		+ w[r66f] - w[r66b] + w[r81f] - w[r81b]
		- w[r95f] + w[r95b] - w[r96f] + w[r96b]
		- w[r97f] + w[r97b] + w[r98f] - w[r98b]
		+ w[r110f] - w[r110b] - w[r118f] + w[r118b]
		- w[r119f] + w[r119b] - w[r124f] + w[r124b]
		- w[r129f] + w[r129b] + w[r136f] - w[r136b]
		- w[r138f] + w[r138b] + 2 * w[r139f] - 2 * w[r139b]
		+ w[r146f] - w[r146b] - w[r149f] + w[r149b]
		+ 2 * w[r150f] - 2 * w[r150b] + w[r154f] - w[r154b]
		- w[r155f] + w[r155b] - w[r156f] + w[r156b]
		- w[r157f] + w[r157b] - 2 * w[r158f] + 2 * w[r158b]
		- 2 * w[r159f] + 2 * w[r159b] - w[r160f] + w[r160b]
		- w[r161f] + w[r161b] - w[r162f] + w[r162b]
		- w[r163f] + w[r163b] - w[r164f] + w[r164b]
		- w[r165f] + w[r165b] - w[r211f] + w[r211b]
		- w[r255f] + w[r255b] - w[r256f] + w[r256b]
		- w[r275f] + w[r275b] - w[r276f] + w[r276b]
		+ w[r289f] - w[r289b] + w[r291] + w[r293];

	cdot[sCH2O] = w[r10f] - w[r10b] - w[r15f] + w[r15b]
		+ w[r16f] - w[r16b] + w[r17f] - w[r17b]
		+ w[r26f] - w[r26b] - w[r32f] + w[r32b]
		+ w[r54f] - w[r54b] - w[r56f] + w[r56b]
		- w[r57f] + w[r57b] - w[r58f] + w[r58b]
		+ w[r60f] - w[r60b] + w[r65f] - w[r65b]
		+ w[r83f] - w[r83b] + w[r92f] - w[r92b]
		+ w[r94f] - w[r94b] - w[r101f] + w[r101b]
		+ w[r102f] - w[r102b] + w[r103f] - w[r103b]
		+ w[r117f] - w[r117b] - w[r121f] + w[r121b]
		+ w[r127f] - w[r127b] - w[r133f] + w[r133b]
		+ w[r153f] - w[r153b] + w[r156f] - w[r156b]
		- w[r161f] + w[r161b] + w[r169f] - w[r169b]
		+ w[r170f] - w[r170b] + w[r173f] - w[r173b]
		+ w[r280] + w[r283] + w[r285];

	cdot[sCH4] = - w[r11f] + w[r11b] + w[r52f] - w[r52b]
		- w[r53f] + w[r53b] - w[r98f] + w[r98b]
		+ w[r118f] - w[r118b] - w[r130f] + w[r130b]
		- w[r139f] + w[r139b] - w[r150f] + w[r150b]
		+ w[r157f] - w[r157b] + w[r160f] - w[r160b]
		+ w[r161f] - w[r161b] + w[r162f] - w[r162b]
		+ w[r163f] - w[r163b] + w[r164f] - w[r164b]
		+ w[r165f] - w[r165b] + w[r211f] - w[r211b];

	cdot[sCO2] = w[r12f] - w[r12b] + w[r14f] - w[r14b]
		+ w[r30f] - w[r30b] + w[r31f] - w[r31b]
		- w[r42f] + w[r42f] - w[r42b] + w[r42b]
		+ w[r99f] - w[r99b] + w[r120f] - w[r120b]
		- w[r132f] + w[r132b] - w[r152f] + w[r152f]
		- w[r152b] + w[r152b] - w[r153f] + w[r153b]
		+ w[r226f] - w[r226b] + w[r229f] - w[r229b]
		+ w[r262f] - w[r262b] + w[r268f] - w[r268b];

	cdot[sCH2OH] = - w[r16f] + w[r16b] + w[r18f] - w[r18b]
		+ w[r56f] - w[r56b] - w[r59f] + w[r59b]
		- w[r60f] + w[r60b] - w[r61f] + w[r61b]
		- w[r62f] + w[r62b] + w[r64f] - w[r64b]
		+ w[r68f] - w[r68b] - w[r102f] + w[r102b]
		+ w[r104f] - w[r104b] + w[r162f] - w[r162b]
		- w[r169f] + w[r169b];

	cdot[sCH3O] = - w[r17f] + w[r17b] + w[r19f] - w[r19b]
		+ w[r57f] - w[r57b] - w[r63f] + w[r63b]
		- w[r64f] + w[r64b] - w[r65f] + w[r65b]
		- w[r66f] + w[r66b] - w[r67f] + w[r67b]
		+ w[r69f] - w[r69b] - w[r103f] + w[r103b]
		+ w[r105f] - w[r105b] + w[r119f] - w[r119b]
		+ w[r155f] - w[r155b] + w[r163f] - w[r163b]
		- w[r170f] + w[r170b];

	cdot[sCH3OH] = - w[r18f] + w[r18b] - w[r19f] + w[r19b]
		+ w[r59f] - w[r59b] + w[r63f] - w[r63b]
		- w[r68f] + w[r68b] - w[r69f] + w[r69b]
		+ w[r95f] - w[r95b] - w[r104f] + w[r104b]
		- w[r105f] + w[r105b] + w[r147f] - w[r147b]
		- w[r162f] + w[r162b] - w[r163f] + w[r163b];

	cdot[sC2H] = - w[r20f] + w[r20b] + w[r22f] - w[r22b]
		- w[r70f] + w[r70b] - w[r106f] + w[r106b]
		+ w[r109f] - w[r109b] + w[r123f] - w[r123b]
		- w[r171f] + w[r171b] - w[r172f] + w[r172b];

	cdot[sC2H2] = - w[r21f] + w[r21b] - w[r22f] + w[r22b]
		- w[r23f] + w[r23b] + w[r70f] - w[r70b]
		- w[r71f] + w[r71b] + w[r73f] - w[r73b]
		- w[r107f] + w[r107b] - w[r108f] + w[r108b]
		- w[r109f] + w[r109b] - w[r110f] + w[r110b]
		+ w[r111f] - w[r111b] + w[r124f] - w[r124b]
		+ w[r128f] - w[r128b] + w[r134f] - w[r134b]
		+ w[r137f] - w[r137b] + w[r172f] - w[r172b]
		+ w[r174f] - w[r174b] + w[r177f] - w[r177b]
		+ w[r283];

	cdot[sHCCO] = w[r21f] - w[r21b] - w[r28f] + w[r28b]
		+ w[r29f] - w[r29b] - w[r79f] + w[r79b]
		+ w[r80f] - w[r80b] + w[r106f] - w[r106b]
		+ w[r114f] - w[r114b] + w[r131f] - w[r131b]
		- w[r134f] + w[r134b] - w[r141f] + w[r141b]
		- w[r176f] + w[r176b] - 2 * w[r177f] + 2 * w[r177b]
		- w[r274f] + w[r274b] + w[r280];

	cdot[sC2H3] = - w[r24f] + w[r24b] + w[r71f] - w[r71b]
		- w[r72f] + w[r72b] - w[r73f] + w[r73b]
		+ w[r75f] - w[r75b] - w[r111f] + w[r111b]
		+ w[r112f] - w[r112b] + w[r129f] - w[r129b]
		+ w[r141f] - w[r141b] + w[r164f] - w[r164b]
		- w[r173f] + w[r173b] + w[r281] + w[r284]
		+ w[r285] + w[r289f] - w[r289b];

	cdot[sCH2CO] = w[r24f] - w[r24b] - w[r29f] + w[r29b]
		- w[r30f] + w[r30b] - w[r80f] + w[r80b]
		- w[r81f] + w[r81b] + w[r82f] - w[r82b]
		+ w[r107f] - w[r107b] - w[r114f] + w[r114b]
		+ w[r133f] - w[r133b] + w[r140f] - w[r140b];

	cdot[sC2H4] = - w[r25f] + w[r25b] + w[r72f] - w[r72b]
		- w[r74f] + w[r74b] - w[r75f] + w[r75b]
		+ w[r77f] - w[r77b] - w[r112f] + w[r112b]
		+ w[r130f] - w[r130b] + w[r138f] - w[r138b]
		+ w[r149f] - w[r149b] - w[r164f] + w[r164b]
		- w[r174f] + w[r174b] + w[r175f] - w[r175b]
		+ w[r286] + w[r291] + w[r293];

	cdot[sC2H5] = - w[r26f] + w[r26b] + w[r27f] - w[r27b]
		+ w[r74f] - w[r74b] - w[r76f] + w[r76b]
		- w[r77f] + w[r77b] + w[r78f] - w[r78b]
		+ w[r113f] - w[r113b] + w[r154f] - w[r154b]
		+ w[r159f] - w[r159b] + w[r165f] - w[r165b]
		- w[r175f] + w[r175b];

	cdot[sC2H6] = - w[r27f] + w[r27b] + w[r76f] - w[r76b]
		- w[r78f] + w[r78b] - w[r113f] + w[r113b]
		- w[r154f] + w[r154b] + w[r158f] - w[r158b]
		- w[r165f] + w[r165b];

	cdot[sH2O] = - w[r35f] + w[r35f] - w[r35b] + w[r35b]
		- w[r41f] + w[r41f] - w[r41b] + w[r41b]
		+ w[r43f] - w[r43b] + w[r44f] - w[r44b]
		+ w[r48f] - w[r48b] + w[r62f] - w[r62b]
		+ w[r67f] - w[r67b] + w[r84f] - w[r84b]
		+ w[r86f] - w[r86b] + w[r87f] - w[r87b]
		+ w[r88f] - w[r88b] + w[r89f] - w[r89b]
		+ w[r93f] - w[r93b] + w[r96f] - w[r96b]
		+ w[r97f] - w[r97b] + w[r98f] - w[r98b]
		+ w[r100f] - w[r100b] + w[r101f] - w[r101b]
		+ w[r102f] - w[r102b] + w[r103f] - w[r103b]
		+ w[r104f] - w[r104b] + w[r105f] - w[r105b]
		+ w[r109f] - w[r109b] + w[r111f] - w[r111b]
		+ w[r112f] - w[r112b] + w[r113f] - w[r113b]
		+ w[r114f] - w[r114b] - w[r127f] + w[r127b]
		+ w[r145f] - w[r145b] - w[r147f] + w[r147b]
		- w[r148f] + w[r148f] - w[r148b] + w[r148b]
		- w[r166f] + w[r166f] - w[r166b] + w[r166b]
		+ w[r193f] - w[r193b] - w[r197f] + w[r197b]
		+ w[r203f] - w[r203b] + w[r210f] - w[r210b]
		+ w[r215f] - w[r215b] - w[r219f] + w[r219b]
		+ w[r255f] - w[r255b] + w[r267f] - w[r267b]
		+ w[r278f] - w[r278b] + w[r299] + w[r300];

	cdot[sC] = w[r49f] - w[r49b] - w[r90f] + w[r90b]
		- w[r122f] + w[r122b] - w[r123f] + w[r123b]
		- w[r124f] + w[r124b] - w[r239f] + w[r239b]
		- w[r244f] + w[r244b] - w[r245f] + w[r245b];

	cdot[sHCCOH] = - w[r82f] + w[r82b] + w[r108f] - w[r108b];

	cdot[sN] = - w[r178f] + w[r178b] - w[r179f] + w[r179b]
		- w[r180f] + w[r180b] + w[r191f] - w[r191b]
		+ w[r193f] - w[r193b] - w[r196f] + w[r196b]
		+ w[r217f] - w[r217b] - w[r225f] + w[r225b]
		+ w[r227f] - w[r227b] - w[r238f] + w[r238b]
		+ w[r239f] - w[r239b] + w[r240f] - w[r240b]
		+ w[r245f] - w[r245b] + w[r248f] - w[r248b]
		- w[r275f] + w[r275b] - w[r276f] + w[r276b];

	cdot[sNO] = - w[r178f] + w[r178b] + w[r179f] - w[r179b]
		+ w[r180f] - w[r180b] + 2 * w[r182f] - 2 * w[r182b]
		- w[r186f] + w[r186b] - w[r187f] + w[r187b]
		+ w[r188f] - w[r188b] + w[r189f] - w[r189b]
		+ w[r190f] - w[r190b] + w[r195f] - w[r195b]
		- w[r198f] + w[r198b] - w[r199f] + w[r199b]
		+ w[r208f] - w[r208b] - w[r212f] + w[r212b]
		+ w[r213f] - w[r213b] + w[r214f] - w[r214b]
		+ w[r215f] - w[r215b] + w[r216f] - w[r216b]
		+ w[r222f] - w[r222b] + w[r224f] - w[r224b]
		+ w[r226f] - w[r226b] - w[r228f] + w[r228b]
		- w[r229f] + w[r229b] - w[r244f] + w[r244b]
		- w[r245f] + w[r245b] - w[r246f] + w[r246b]
		- w[r247f] + w[r247b] - w[r248f] + w[r248b]
		- w[r249f] + w[r249b] - w[r250f] + w[r250b]
		- w[r251f] + w[r251b] - w[r252f] + w[r252b]
		- w[r253f] + w[r253b] - w[r254f] + w[r254b]
		- w[r255f] + w[r255b] - w[r256f] + w[r256b]
		+ w[r258f] - w[r258b] - w[r274f] + w[r274b];

	cdot[sN2O] = - w[r181f] + w[r181b] - w[r182f] + w[r182b]
		- w[r183f] + w[r183b] - w[r184f] + w[r184b]
		- w[r185f] + w[r185b] + w[r199f] - w[r199b]
		+ w[r228f] - w[r228b];

	cdot[sNO2] = w[r186f] - w[r186b] + w[r187f] - w[r187b]
		- w[r188f] + w[r188b] - w[r189f] + w[r189b];

	cdot[sNH] = - w[r190f] + w[r190b] - w[r191f] + w[r191b]
		- w[r192f] + w[r192b] - w[r193f] + w[r193b]
		- w[r194f] + w[r194b] - w[r195f] + w[r195b]
		- w[r196f] + w[r196b] - w[r197f] + w[r197b]
		- w[r198f] + w[r198b] - w[r199f] + w[r199b]
		+ w[r200f] - w[r200b] + w[r202f] - w[r202b]
		+ w[r203f] - w[r203b] + w[r208f] - w[r208b]
		+ w[r223f] - w[r223b] + w[r232f] - w[r232b]
		+ w[r242f] - w[r242b] + w[r243f] - w[r243b]
		+ w[r262f] - w[r262b] + w[r269f] - w[r269b];

	cdot[sHNO] = w[r192f] - w[r192b] + w[r194f] - w[r194b]
		+ w[r197f] - w[r197b] + w[r201f] - w[r201b]
		+ w[r212f] - w[r212b] - w[r213f] + w[r213b]
		- w[r214f] + w[r214b] - w[r215f] + w[r215b]
		- w[r216f] + w[r216b] + w[r263f] - w[r263b];

	cdot[sNH2] = - w[r200f] + w[r200b] - w[r201f] + w[r201b]
		- w[r202f] + w[r202b] - w[r203f] + w[r203b]
		+ w[r236f] - w[r236b] + w[r265f] - w[r265b]
		+ w[r268f] - w[r268b] + w[r272f] - w[r272b]
		+ w[r277f] - w[r277b] + w[r278f] - w[r278b]
		+ w[r279f] - w[r279b];

	cdot[sNNH] = - w[r204f] + w[r204b] - w[r205f] + w[r205b]
		- w[r206f] + w[r206b] - w[r207f] + w[r207b]
		- w[r208f] + w[r208b] - w[r209f] + w[r209b]
		- w[r210f] + w[r210b] - w[r211f] + w[r211b];

	cdot[sCN] = - w[r217f] + w[r217b] - w[r218f] + w[r218b]
		- w[r219f] + w[r219b] - w[r220f] + w[r220b]
		- w[r221f] + w[r221b] + w[r230f] - w[r230b]
		+ w[r233f] - w[r233b] + w[r239f] - w[r239b]
		+ w[r244f] - w[r244b];

	cdot[sNCO] = w[r218f] - w[r218b] + w[r220f] - w[r220b]
		- w[r222f] + w[r222b] - w[r223f] + w[r223b]
		- w[r224f] + w[r224b] - w[r225f] + w[r225b]
		- w[r226f] + w[r226b] - w[r227f] + w[r227b]
		- w[r228f] + w[r228b] - w[r229f] + w[r229b]
		+ w[r231f] - w[r231b] + w[r247f] - w[r247b]
		+ w[r264f] - w[r264b] + w[r266f] - w[r266b]
		+ w[r267f] - w[r267b];

	cdot[sHCN] = w[r219f] - w[r219b] + w[r221f] - w[r221b]
		- w[r230f] + w[r230b] - w[r231f] + w[r231b]
		- w[r232f] + w[r232b] - w[r233f] + w[r233b]
		- w[r234f] + w[r234b] - w[r235f] + w[r235b]
		- w[r236f] + w[r236b] - w[r237f] + w[r237b]
		+ w[r240f] - w[r240b] + w[r242f] - w[r242b]
		+ w[r243f] - w[r243b] + w[r246f] - w[r246b]
		+ w[r250f] - w[r250b] + w[r253f] - w[r253b]
		+ w[r255f] - w[r255b] + w[r258f] - w[r258b]
		+ w[r271f] - w[r271b] + w[r276f] - w[r276b];

	cdot[sHOCN] = w[r234f] - w[r234b] - w[r273f] + w[r273b];

	cdot[sHNCO] = w[r235f] - w[r235b] + w[r249f] - w[r249b]
		+ w[r252f] - w[r252b] - w[r262f] + w[r262b]
		- w[r263f] + w[r263b] - w[r264f] + w[r264b]
		- w[r265f] + w[r265b] - w[r266f] + w[r266b]
		- w[r267f] + w[r267b] - w[r268f] + w[r268b]
		- w[r269f] + w[r269b] + w[r270f] - w[r270b]
		+ w[r273f] - w[r273b];

	cdot[sH2CN] = w[r237f] - w[r237b] - w[r238f] + w[r238b]
		+ w[r256f] - w[r256b] + w[r275f] - w[r275b];

	cdot[sHCNN] = w[r241f] - w[r241b] - w[r257f] + w[r257b]
		- w[r258f] + w[r258b] - w[r259f] + w[r259b]
		- w[r260f] + w[r260b] - w[r261f] + w[r261b];

	cdot[sHCNO] = w[r251f] - w[r251b] + w[r254f] - w[r254b]
		- w[r270f] + w[r270b] - w[r271f] + w[r271b]
		- w[r272f] + w[r272b] + w[r274f] - w[r274b];

	cdot[sNH3] = - w[r277f] + w[r277b] - w[r278f] + w[r278b]
		- w[r279f] + w[r279b];

	cdot[sC3H3] = - w[r280] - w[r281] + w[r282f] - w[r282b];

	cdot[sC3H4] = - w[r282f] + w[r282b] - w[r283] - w[r284]
		- w[r285] - w[r286] + w[r287f] - w[r287b]
		+ w[r288];

	cdot[sC3H5] = - w[r287f] + w[r287b] - w[r288] + w[r290];

	cdot[sC3H6] = - w[r289f] + w[r289b] - w[r290] + w[r292f]
		- w[r292b] + w[r294];

	cdot[sNXC3H7] = - w[r291] - w[r292f] + w[r292b] + w[r295]
		+ w[r297] + w[r299];

	cdot[sIXC3H7] = - w[r293] - w[r294] + w[r296] + w[r298]
		+ w[r300];

	cdot[sC3H8] = - w[r295] - w[r296] - w[r297] - w[r298]
		- w[r299] - w[r300];

}

double GetLindRateCoeff( double temp, double pressure
				, double k0, double kInf
				, double fc, double conc )
{
	const double	R = 8314.34;   /* [J / kmole K] */
	double			Ntmp;
	double			kl;
	double			f;
	double			cCoeff, dCoeff, log10kNull;

	int				iTroe = 1;

	if ( conc <= 0.0 ) {
		conc = pressure / ( R * temp );
	}
	Ntmp = 0.75 - 1.27 * log10( fc );
	if ( iTroe ) {
		cCoeff = - 0.4 - 0.67 * log10( fc );
		dCoeff = 0.14;
		k0 *= conc / MAX_C(kInf, 1.0e-60);
		//k0 *= conc / kInf;
		log10kNull = log10( k0 );
		f = ( log10kNull + cCoeff ) / ( Ntmp - dCoeff * ( log10kNull + cCoeff ) );
		f = pow( fc, 1.0 / ( f * f + 1.0 ) );
		kInf *= f * k0 / ( 1.0 + k0 );
	}
	else {
		k0 = k0 * conc / kInf;
		kl = k0 / ( 1.0 + k0 );
		f = log10( k0 ) / Ntmp;
		f = pow( fc, 1.0 / ( f * f + 1.0 ) );
		kInf = kInf * f * kl;

	}
	return kInf;
	
}

/*
double CatchZero( double a )
{
	return ( a == 0.0 ) ? 1.0e-20 : a;
}

*/
void GetMolarMass( double *W )
{
	W[sN2] =  2.80200000e+01;
	W[sAR] =  3.99480000e+01;
	W[sO] =  1.60000000e+01;
	W[sO2] =  3.20000000e+01;
	W[sH] =  1.00800000e+00;
	W[sOH] =  1.70080000e+01;
	W[sH2] =  2.01600000e+00;
	W[sHO2] =  3.30080000e+01;
	W[sH2O2] =  3.40160000e+01;
	W[sCH] =  1.30180000e+01;
	W[sCO] =  2.80100000e+01;
	W[sCH2] =  1.40260000e+01;
	W[sHCO] =  2.90180000e+01;
	W[sCH2GSGXCH2] =  1.40260000e+01;
	W[sCH3] =  1.50340000e+01;
	W[sCH2O] =  3.00260000e+01;
	W[sCH4] =  1.60420000e+01;
	W[sCO2] =  4.40100000e+01;
	W[sCH2OH] =  3.10340000e+01;
	W[sCH3O] =  3.10340000e+01;
	W[sCH3OH] =  3.20420000e+01;
	W[sC2H] =  2.50280000e+01;
	W[sC2H2] =  2.60360000e+01;
	W[sHCCO] =  4.10280000e+01;
	W[sC2H3] =  2.70440000e+01;
	W[sCH2CO] =  4.20360000e+01;
	W[sC2H4] =  2.80520000e+01;
	W[sC2H5] =  2.90600000e+01;
	W[sC2H6] =  3.00680000e+01;
	W[sH2O] =  1.80160000e+01;
	W[sC] =  1.20100000e+01;
	W[sHCCOH] =  4.20360000e+01;
	W[sN] =  1.40100000e+01;
	W[sNO] =  3.00100000e+01;
	W[sN2O] =  4.40200000e+01;
	W[sNO2] =  4.60100000e+01;
	W[sNH] =  1.50180000e+01;
	W[sHNO] =  3.10180000e+01;
	W[sNH2] =  1.60260000e+01;
	W[sNNH] =  2.90280000e+01;
	W[sCN] =  2.60200000e+01;
	W[sNCO] =  4.20200000e+01;
	W[sHCN] =  2.70280000e+01;
	W[sHOCN] =  4.30280000e+01;
	W[sHNCO] =  4.30280000e+01;
	W[sH2CN] =  2.80360000e+01;
	W[sHCNN] =  4.10380000e+01;
	W[sHCNO] =  4.30280000e+01;
	W[sNH3] =  1.70340000e+01;
	W[sC3H3] =  3.90540000e+01;
	W[sC3H4] =  4.00620000e+01;
	W[sC3H5] =  4.10700000e+01;
	W[sC3H6] =  4.20780000e+01;
	W[sNXC3H7] =  4.30860000e+01;
	W[sIXC3H7] =  4.30860000e+01;
	W[sC3H8] =  4.40940000e+01;
}

void GetSpeciesNames( char **names )
{
	strcpy( names[sN2], "N2                  " );
	strcpy( names[sAR], "AR                  " );
	strcpy( names[sO], "O                   " );
	strcpy( names[sO2], "O2                  " );
	strcpy( names[sH], "H                   " );
	strcpy( names[sOH], "OH                  " );
	strcpy( names[sH2], "H2                  " );
	strcpy( names[sHO2], "HO2                 " );
	strcpy( names[sH2O2], "H2O2                " );
	strcpy( names[sCH], "CH                  " );
	strcpy( names[sCO], "CO                  " );
	strcpy( names[sCH2], "CH2                 " );
	strcpy( names[sHCO], "HCO                 " );
	strcpy( names[sCH2GSGXCH2], "CH2GSG-CH2          " );
	strcpy( names[sCH3], "CH3                 " );
	strcpy( names[sCH2O], "CH2O                " );
	strcpy( names[sCH4], "CH4                 " );
	strcpy( names[sCO2], "CO2                 " );
	strcpy( names[sCH2OH], "CH2OH               " );
	strcpy( names[sCH3O], "CH3O                " );
	strcpy( names[sCH3OH], "CH3OH               " );
	strcpy( names[sC2H], "C2H                 " );
	strcpy( names[sC2H2], "C2H2                " );
	strcpy( names[sHCCO], "HCCO                " );
	strcpy( names[sC2H3], "C2H3                " );
	strcpy( names[sCH2CO], "CH2CO               " );
	strcpy( names[sC2H4], "C2H4                " );
	strcpy( names[sC2H5], "C2H5                " );
	strcpy( names[sC2H6], "C2H6                " );
	strcpy( names[sH2O], "H2O                 " );
	strcpy( names[sC], "C                   " );
	strcpy( names[sHCCOH], "HCCOH               " );
	strcpy( names[sN], "N                   " );
	strcpy( names[sNO], "NO                  " );
	strcpy( names[sN2O], "N2O                 " );
	strcpy( names[sNO2], "NO2                 " );
	strcpy( names[sNH], "NH                  " );
	strcpy( names[sHNO], "HNO                 " );
	strcpy( names[sNH2], "NH2                 " );
	strcpy( names[sNNH], "NNH                 " );
	strcpy( names[sCN], "CN                  " );
	strcpy( names[sNCO], "NCO                 " );
	strcpy( names[sHCN], "HCN                 " );
	strcpy( names[sHOCN], "HOCN                " );
	strcpy( names[sHNCO], "HNCO                " );
	strcpy( names[sH2CN], "H2CN                " );
	strcpy( names[sHCNN], "HCNN                " );
	strcpy( names[sHCNO], "HCNO                " );
	strcpy( names[sNH3], "NH3                 " );
	strcpy( names[sC3H3], "C3H3                " );
	strcpy( names[sC3H4], "C3H4                " );
	strcpy( names[sC3H5], "C3H5                " );
	strcpy( names[sC3H6], "C3H6                " );
	strcpy( names[sNXC3H7], "N-C3H7              " );
	strcpy( names[sIXC3H7], "I-C3H7              " );
	strcpy( names[sC3H8], "C3H8                " );
//	strcpy( names[sT], "Temp                " );
}



void ComputeThermoData( double *h, double *cp, double T )
{
/*
	This function computes enthalpy 'h' and heat capacity 'cp' as
	function of temperature 'T' for all non steady state species
	in units [J/kg] and [J/kg K], respectively.
	The parameter h and cp should provide workspace of length 56 */


	int i;
	if ( T > 1000.0 ) {
		h[sN2] =  2.96728765e+02 * (
			T * (  2.92664000e+00 + T * (  7.43988400e-04
			+ T * ( -1.89492000e-07 + T * (  2.52425950e-11
			+ T * -1.35067020e-15 ) ) ) ) -9.22797700e+02 );
		cp[sN2] =  2.96728765e+02 * (
			 2.92664000e+00 + T * (  1.48797680e-03 
			+ T * ( -5.68476000e-07 + T * (  1.00970380e-10
			+ T * -6.75335100e-15 ) ) ) );
		h[sAR] =  2.08129068e+02 * (
			T * (  2.50000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) -7.45375000e+02 );
		cp[sAR] =  2.08129068e+02 * (
			 2.50000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sO] =  5.19646250e+02 * (
			T * (  2.56942078e+00 + T * ( -4.29870569e-05
			+ T * (  1.39828196e-08 + T * ( -2.50444497e-12
			+ T *  2.45667382e-16 ) ) ) ) +  2.92175791e+04 );
		cp[sO] =  5.19646250e+02 * (
			 2.56942078e+00 + T * ( -8.59741137e-05 
			+ T * (  4.19484589e-08 + T * ( -1.00177799e-11
			+ T *  1.22833691e-15 ) ) ) );
		h[sO2] =  2.59823125e+02 * (
			T * (  3.28253784e+00 + T * (  7.41543770e-04
			+ T * ( -2.52655556e-07 + T * (  5.23676387e-11
			+ T * -4.33435588e-15 ) ) ) ) -1.08845772e+03 );
		cp[sO2] =  2.59823125e+02 * (
			 3.28253784e+00 + T * (  1.48308754e-03 
			+ T * ( -7.57966669e-07 + T * (  2.09470555e-10
			+ T * -2.16717794e-14 ) ) ) );
		h[sH] =  8.24835317e+03 * (
			T * (  2.50000001e+00 + T * ( -1.15421486e-11
			+ T * (  5.38539827e-15 + T * ( -1.18378809e-18
			+ T *  9.96394714e-23 ) ) ) ) +  2.54736599e+04 );
		cp[sH] =  8.24835317e+03 * (
			 2.50000001e+00 + T * ( -2.30842973e-11 
			+ T * (  1.61561948e-14 + T * ( -4.73515235e-18
			+ T *  4.98197357e-22 ) ) ) );
		h[sOH] =  4.88848777e+02 * (
			T * (  3.09288767e+00 + T * (  2.74214858e-04
			+ T * (  4.21684093e-08 + T * ( -2.19865389e-11
			+ T *  2.34824752e-15 ) ) ) ) +  3.85865700e+03 );
		cp[sOH] =  4.88848777e+02 * (
			 3.09288767e+00 + T * (  5.48429716e-04 
			+ T * (  1.26505228e-07 + T * ( -8.79461556e-11
			+ T *  1.17412376e-14 ) ) ) );
		h[sH2] =  4.12417659e+03 * (
			T * (  3.33727920e+00 + T * ( -2.47012365e-05
			+ T * (  1.66485593e-07 + T * ( -4.48915985e-11
			+ T *  4.00510752e-15 ) ) ) ) -9.50158922e+02 );
		cp[sH2] =  4.12417659e+03 * (
			 3.33727920e+00 + T * ( -4.94024731e-05 
			+ T * (  4.99456778e-07 + T * ( -1.79566394e-10
			+ T *  2.00255376e-14 ) ) ) );
		h[sHO2] =  2.51888633e+02 * (
			T * (  4.01721090e+00 + T * (  1.11991006e-03
			+ T * ( -2.11219383e-07 + T * (  2.85615925e-11
			+ T * -2.15817070e-15 ) ) ) ) +  1.11856713e+02 );
		cp[sHO2] =  2.51888633e+02 * (
			 4.01721090e+00 + T * (  2.23982013e-03 
			+ T * ( -6.33658150e-07 + T * (  1.14246370e-10
			+ T * -1.07908535e-14 ) ) ) );
		h[sH2O2] =  2.44424389e+02 * (
			T * (  4.16500285e+00 + T * (  2.45415847e-03
			+ T * ( -6.33797417e-07 + T * (  9.27964965e-11
			+ T * -5.75816610e-15 ) ) ) ) -1.78617877e+04 );
		cp[sH2O2] =  2.44424389e+02 * (
			 4.16500285e+00 + T * (  4.90831694e-03 
			+ T * ( -1.90139225e-06 + T * (  3.71185986e-10
			+ T * -2.87908305e-14 ) ) ) );
		h[sCH] =  6.38680289e+02 * (
			T * (  2.87846473e+00 + T * (  4.85456840e-04
			+ T * (  4.81485517e-08 + T * ( -3.26719623e-11
			+ T *  3.52158766e-15 ) ) ) ) +  7.10124364e+04 );
		cp[sCH] =  6.38680289e+02 * (
			 2.87846473e+00 + T * (  9.70913681e-04 
			+ T * (  1.44445655e-07 + T * ( -1.30687849e-10
			+ T *  1.76079383e-14 ) ) ) );
		h[sCO] =  2.96834702e+02 * (
			T * (  2.71518561e+00 + T * (  1.03126372e-03
			+ T * ( -3.32941924e-07 + T * (  5.75132520e-11
			+ T * -4.07295432e-15 ) ) ) ) -1.41518724e+04 );
		cp[sCO] =  2.96834702e+02 * (
			 2.71518561e+00 + T * (  2.06252743e-03 
			+ T * ( -9.98825771e-07 + T * (  2.30053008e-10
			+ T * -2.03647716e-14 ) ) ) );
		h[sCH2] =  5.92780550e+02 * (
			T * (  2.87410113e+00 + T * (  1.82819646e-03
			+ T * ( -4.69648657e-07 + T * (  6.50448872e-11
			+ T * -3.75455134e-15 ) ) ) ) +  4.62636040e+04 );
		cp[sCH2] =  5.92780550e+02 * (
			 2.87410113e+00 + T * (  3.65639292e-03 
			+ T * ( -1.40894597e-06 + T * (  2.60179549e-10
			+ T * -1.87727567e-14 ) ) ) );
		h[sHCO] =  2.86523537e+02 * (
			T * (  2.77217438e+00 + T * (  2.47847763e-03
			+ T * ( -8.28152043e-07 + T * (  1.47290445e-10
			+ T * -1.06701742e-14 ) ) ) ) +  4.01191815e+03 );
		cp[sHCO] =  2.86523537e+02 * (
			 2.77217438e+00 + T * (  4.95695526e-03 
			+ T * ( -2.48445613e-06 + T * (  5.89161778e-10
			+ T * -5.33508711e-14 ) ) ) );
		h[sCH2GSGXCH2] =  5.92780550e+02 * (
			T * (  2.29203842e+00 + T * (  2.32794318e-03
			+ T * ( -6.70639823e-07 + T * (  1.04476500e-10
			+ T * -6.79432730e-15 ) ) ) ) +  5.09259997e+04 );
		cp[sCH2GSGXCH2] =  5.92780550e+02 * (
			 2.29203842e+00 + T * (  4.65588637e-03 
			+ T * ( -2.01191947e-06 + T * (  4.17906000e-10
			+ T * -3.39716365e-14 ) ) ) );
		h[sCH3] =  5.53035786e+02 * (
			T * (  2.28571772e+00 + T * (  3.61995018e-03
			+ T * ( -9.95714493e-07 + T * (  1.48921161e-10
			+ T * -9.34308788e-15 ) ) ) ) +  1.67755843e+04 );
		cp[sCH3] =  5.53035786e+02 * (
			 2.28571772e+00 + T * (  7.23990037e-03 
			+ T * ( -2.98714348e-06 + T * (  5.95684644e-10
			+ T * -4.67154394e-14 ) ) ) );
		h[sCH2O] =  2.76904683e+02 * (
			T * (  1.76069008e+00 + T * (  4.60000041e-03
			+ T * ( -1.47419604e-06 + T * (  2.51603030e-10
			+ T * -1.76771128e-14 ) ) ) ) -1.39958323e+04 );
		cp[sCH2O] =  2.76904683e+02 * (
			 1.76069008e+00 + T * (  9.20000082e-03 
			+ T * ( -4.42258813e-06 + T * (  1.00641212e-09
			+ T * -8.83855640e-14 ) ) ) );
		h[sCH4] =  5.18285750e+02 * (
			T * (  7.48514950e-02 + T * (  6.69547335e-03
			+ T * ( -1.91095270e-06 + T * (  3.05731338e-10
			+ T * -2.03630460e-14 ) ) ) ) -9.46834459e+03 );
		cp[sCH4] =  5.18285750e+02 * (
			 7.48514950e-02 + T * (  1.33909467e-02 
			+ T * ( -5.73285809e-06 + T * (  1.22292535e-09
			+ T * -1.01815230e-13 ) ) ) );
		h[sCO2] =  1.88919337e+02 * (
			T * (  3.85746029e+00 + T * (  2.20718513e-03
			+ T * ( -7.38271347e-07 + T * (  1.30872547e-10
			+ T * -9.44168328e-15 ) ) ) ) -4.87591660e+04 );
		cp[sCO2] =  1.88919337e+02 * (
			 3.85746029e+00 + T * (  4.41437026e-03 
			+ T * ( -2.21481404e-06 + T * (  5.23490188e-10
			+ T * -4.72084164e-14 ) ) ) );
		h[sCH2OH] =  2.67910679e+02 * (
			T * (  3.69266569e+00 + T * (  4.32288399e-03
			+ T * ( -1.25033707e-06 + T * (  1.96808659e-10
			+ T * -1.29710840e-14 ) ) ) ) -3.24250627e+03 );
		cp[sCH2OH] =  2.67910679e+02 * (
			 3.69266569e+00 + T * (  8.64576797e-03 
			+ T * ( -3.75101120e-06 + T * (  7.87234636e-10
			+ T * -6.48554201e-14 ) ) ) );
		h[sCH3O] =  2.67910679e+02 * (
			T * (  3.77079900e+00 + T * (  3.93574850e-03
			+ T * ( -8.85461333e-07 + T * (  9.86107750e-11
			+ T * -4.22523200e-15 ) ) ) ) +  1.27832520e+02 );
		cp[sCH3O] =  2.67910679e+02 * (
			 3.77079900e+00 + T * (  7.87149700e-03 
			+ T * ( -2.65638400e-06 + T * (  3.94443100e-10
			+ T * -2.11261600e-14 ) ) ) );
		h[sCH3OH] =  2.59482554e+02 * (
			T * (  1.78970791e+00 + T * (  7.04691460e-03
			+ T * ( -2.12166945e-06 + T * (  3.45427713e-10
			+ T * -2.34120440e-14 ) ) ) ) -2.53748747e+04 );
		cp[sCH3OH] =  2.59482554e+02 * (
			 1.78970791e+00 + T * (  1.40938292e-02 
			+ T * ( -6.36500835e-06 + T * (  1.38171085e-09
			+ T * -1.17060220e-13 ) ) ) );
		h[sC2H] =  3.32201534e+02 * (
			T * (  3.16780652e+00 + T * (  2.37610951e-03
			+ T * ( -6.12623590e-07 + T * (  7.60475630e-11
			+ T * -3.54465540e-15 ) ) ) ) +  6.71210650e+04 );
		cp[sC2H] =  3.32201534e+02 * (
			 3.16780652e+00 + T * (  4.75221902e-03 
			+ T * ( -1.83787077e-06 + T * (  3.04190252e-10
			+ T * -1.77232770e-14 ) ) ) );
		h[sC2H2] =  3.19340144e+02 * (
			T * (  4.14756964e+00 + T * (  2.98083332e-03
			+ T * ( -7.90982840e-07 + T * (  1.16853043e-10
			+ T * -7.22470426e-15 ) ) ) ) +  2.59359992e+04 );
		cp[sC2H2] =  3.19340144e+02 * (
			 4.14756964e+00 + T * (  5.96166664e-03 
			+ T * ( -2.37294852e-06 + T * (  4.67412171e-10
			+ T * -3.61235213e-14 ) ) ) );
		h[sHCCO] =  2.02650385e+02 * (
			T * (  5.62820580e+00 + T * (  2.04267005e-03
			+ T * ( -5.31151567e-07 + T * (  7.15651300e-11
			+ T * -3.88156640e-15 ) ) ) ) +  1.93272150e+04 );
		cp[sHCCO] =  2.02650385e+02 * (
			 5.62820580e+00 + T * (  4.08534010e-03 
			+ T * ( -1.59345470e-06 + T * (  2.86260520e-10
			+ T * -1.94078320e-14 ) ) ) );
		h[sC2H3] =  3.07437509e+02 * (
			T * (  3.01672400e+00 + T * (  5.16511460e-03
			+ T * ( -1.56027450e-06 + T * (  2.54408220e-10
			+ T * -1.72521408e-14 ) ) ) ) +  3.46128739e+04 );
		cp[sC2H3] =  3.07437509e+02 * (
			 3.01672400e+00 + T * (  1.03302292e-02 
			+ T * ( -4.68082349e-06 + T * (  1.01763288e-09
			+ T * -8.62607041e-14 ) ) ) );
		h[sCH2CO] =  1.97790941e+02 * (
			T * (  4.51129732e+00 + T * (  4.50179872e-03
			+ T * ( -1.38979878e-06 + T * (  2.30836470e-10
			+ T * -1.58967640e-14 ) ) ) ) -7.55105311e+03 );
		cp[sCH2CO] =  1.97790941e+02 * (
			 4.51129732e+00 + T * (  9.00359745e-03 
			+ T * ( -4.16939635e-06 + T * (  9.23345882e-10
			+ T * -7.94838201e-14 ) ) ) );
		h[sC2H4] =  2.96390275e+02 * (
			T * (  2.03611116e+00 + T * (  7.32270755e-03
			+ T * ( -2.23692638e-06 + T * (  3.68057308e-10
			+ T * -2.51412122e-14 ) ) ) ) +  4.93988614e+03 );
		cp[sC2H4] =  2.96390275e+02 * (
			 2.03611116e+00 + T * (  1.46454151e-02 
			+ T * ( -6.71077915e-06 + T * (  1.47222923e-09
			+ T * -1.25706061e-13 ) ) ) );
		h[sC2H5] =  2.86109429e+02 * (
			T * (  1.95465642e+00 + T * (  8.69863610e-03
			+ T * ( -2.66068889e-06 + T * (  4.38044223e-10
			+ T * -2.99283152e-14 ) ) ) ) +  1.28575200e+04 );
		cp[sC2H5] =  2.86109429e+02 * (
			 1.95465642e+00 + T * (  1.73972722e-02 
			+ T * ( -7.98206668e-06 + T * (  1.75217689e-09
			+ T * -1.49641576e-13 ) ) ) );
		h[sC2H6] =  2.76517893e+02 * (
			T * (  1.07188150e+00 + T * (  1.08426339e-02
			+ T * ( -3.34186890e-06 + T * (  5.53530003e-10
			+ T * -3.80005780e-14 ) ) ) ) -1.14263932e+04 );
		cp[sC2H6] =  2.76517893e+02 * (
			 1.07188150e+00 + T * (  2.16852677e-02 
			+ T * ( -1.00256067e-05 + T * (  2.21412001e-09
			+ T * -1.90002890e-13 ) ) ) );
		h[sH2O] =  4.61497558e+02 * (
			T * (  3.03399249e+00 + T * (  1.08845902e-03
			+ T * ( -5.46908393e-08 + T * ( -2.42604967e-11
			+ T *  3.36401984e-15 ) ) ) ) -3.00042971e+04 );
		cp[sH2O] =  4.61497558e+02 * (
			 3.03399249e+00 + T * (  2.17691804e-03 
			+ T * ( -1.64072518e-07 + T * ( -9.70419870e-11
			+ T *  1.68200992e-14 ) ) ) );
		h[sC] =  6.92284763e+02 * (
			T * (  2.49266888e+00 + T * (  2.39944642e-05
			+ T * ( -2.41445007e-08 + T * (  9.35727573e-12
			+ T * -9.74555786e-16 ) ) ) ) +  8.54512953e+04 );
		cp[sC] =  6.92284763e+02 * (
			 2.49266888e+00 + T * (  4.79889284e-05 
			+ T * ( -7.24335020e-08 + T * (  3.74291029e-11
			+ T * -4.87277893e-15 ) ) ) );
		h[sHCCOH] =  1.97790941e+02 * (
			T * (  5.92382910e+00 + T * (  3.39618000e-03
			+ T * ( -8.55285467e-07 + T * (  1.12469603e-10
			+ T * -5.98802020e-15 ) ) ) ) +  7.26462600e+03 );
		cp[sHCCOH] =  1.97790941e+02 * (
			 5.92382910e+00 + T * (  6.79236000e-03 
			+ T * ( -2.56585640e-06 + T * (  4.49878410e-10
			+ T * -2.99401010e-14 ) ) ) );
		h[sN] =  5.93457530e+02 * (
			T * (  2.41594290e+00 + T * (  8.74453250e-05
			+ T * ( -3.96745633e-08 + T * (  7.55656125e-12
			+ T * -4.07219640e-16 ) ) ) ) +  5.61337730e+04 );
		cp[sN] =  5.93457530e+02 * (
			 2.41594290e+00 + T * (  1.74890650e-04 
			+ T * ( -1.19023690e-07 + T * (  3.02262450e-11
			+ T * -2.03609820e-15 ) ) ) );
		h[sNO] =  2.77052316e+02 * (
			T * (  3.26060560e+00 + T * (  5.95552150e-04
			+ T * ( -1.43056827e-07 + T * (  1.73644173e-11
			+ T * -8.06721980e-16 ) ) ) ) +  9.92097460e+03 );
		cp[sNO] =  2.77052316e+02 * (
			 3.26060560e+00 + T * (  1.19110430e-03 
			+ T * ( -4.29170480e-07 + T * (  6.94576690e-11
			+ T * -4.03360990e-15 ) ) ) );
		h[sN2O] =  1.88876420e+02 * (
			T * (  4.82307290e+00 + T * (  1.31351255e-03
			+ T * ( -3.19502913e-07 + T * (  4.00017800e-11
			+ T * -1.95504606e-15 ) ) ) ) +  8.07340480e+03 );
		cp[sN2O] =  1.88876420e+02 * (
			 4.82307290e+00 + T * (  2.62702510e-03 
			+ T * ( -9.58508740e-07 + T * (  1.60007120e-10
			+ T * -9.77523030e-15 ) ) ) );
		h[sNO2] =  1.80707238e+02 * (
			T * (  4.88475420e+00 + T * (  1.08619780e-03
			+ T * ( -2.76023020e-07 + T * (  3.93687750e-11
			+ T * -2.10217900e-15 ) ) ) ) +  2.31649830e+03 );
		cp[sNO2] =  1.80707238e+02 * (
			 4.88475420e+00 + T * (  2.17239560e-03 
			+ T * ( -8.28069060e-07 + T * (  1.57475100e-10
			+ T * -1.05108950e-14 ) ) ) );
		h[sNH] =  5.53624983e+02 * (
			T * (  2.78369280e+00 + T * (  6.64921500e-04
			+ T * ( -1.41593490e-07 + T * (  1.95871253e-11
			+ T * -1.10088940e-15 ) ) ) ) +  4.21208480e+04 );
		cp[sNH] =  5.53624983e+02 * (
			 2.78369280e+00 + T * (  1.32984300e-03 
			+ T * ( -4.24780470e-07 + T * (  7.83485010e-11
			+ T * -5.50444700e-15 ) ) ) );
		h[sHNO] =  2.68048875e+02 * (
			T * (  2.97925090e+00 + T * (  1.74720295e-03
			+ T * ( -2.61832593e-07 + T * (  1.43698985e-11
			+ T * -3.86718320e-17 ) ) ) ) +  1.17505820e+04 );
		cp[sHNO] =  2.68048875e+02 * (
			 2.97925090e+00 + T * (  3.49440590e-03 
			+ T * ( -7.85497780e-07 + T * (  5.74795940e-11
			+ T * -1.93359160e-16 ) ) ) );
		h[sNH2] =  5.18803195e+02 * (
			T * (  2.83474210e+00 + T * (  1.60365410e-03
			+ T * ( -3.11302680e-07 + T * (  3.42573825e-11
			+ T * -1.58412288e-15 ) ) ) ) +  2.21719570e+04 );
		cp[sNH2] =  5.18803195e+02 * (
			 2.83474210e+00 + T * (  3.20730820e-03 
			+ T * ( -9.33908040e-07 + T * (  1.37029530e-10
			+ T * -7.92061440e-15 ) ) ) );
		h[sNNH] =  2.86424831e+02 * (
			T * (  3.76675440e+00 + T * (  1.44575410e-03
			+ T * ( -3.47220667e-07 + T * (  4.21064850e-11
			+ T * -2.01837920e-15 ) ) ) ) +  2.86506970e+04 );
		cp[sNNH] =  2.86424831e+02 * (
			 3.76675440e+00 + T * (  2.89150820e-03 
			+ T * ( -1.04166200e-06 + T * (  1.68425940e-10
			+ T * -1.00918960e-14 ) ) ) );
		h[sCN] =  3.19536510e+02 * (
			T * (  3.74598050e+00 + T * (  2.17253875e-05
			+ T * (  9.90199467e-08 + T * ( -1.71629515e-11
			+ T *  8.82683460e-16 ) ) ) ) +  5.15361880e+04 );
		cp[sCN] =  3.19536510e+02 * (
			 3.74598050e+00 + T * (  4.34507750e-05 
			+ T * (  2.97059840e-07 + T * ( -6.86518060e-11
			+ T *  4.41341730e-15 ) ) ) );
		h[sNCO] =  1.97866254e+02 * (
			T * (  5.15218450e+00 + T * (  1.15258805e-03
			+ T * ( -2.93443843e-07 + T * (  3.69727450e-11
			+ T * -1.81955992e-15 ) ) ) ) +  1.40041230e+04 );
		cp[sNCO] =  1.97866254e+02 * (
			 5.15218450e+00 + T * (  2.30517610e-03 
			+ T * ( -8.80331530e-07 + T * (  1.47890980e-10
			+ T * -9.09779960e-15 ) ) ) );
		h[sHCN] =  3.07619506e+02 * (
			T * (  3.80223920e+00 + T * (  1.57321140e-03
			+ T * ( -3.54406167e-07 + T * (  4.15493925e-11
			+ T * -1.95995140e-15 ) ) ) ) +  1.49105120e+04 );
		cp[sHCN] =  3.07619506e+02 * (
			 3.80223920e+00 + T * (  3.14642280e-03 
			+ T * ( -1.06321850e-06 + T * (  1.66197570e-10
			+ T * -9.79975700e-15 ) ) ) );
		h[sHOCN] =  1.93230919e+02 * (
			T * (  3.76643163e+00 + T * (  3.47495871e-03
			+ T * ( -1.09814582e-06 + T * (  1.41267588e-10
			+ T *  0.00000000e+00 ) ) ) ) -2.82204704e+03 );
		cp[sHOCN] =  1.93230919e+02 * (
			 3.76643163e+00 + T * (  6.94991743e-03 
			+ T * ( -3.29443746e-06 + T * (  5.65070350e-10
			+ T *  0.00000000e+00 ) ) ) );
		h[sHNCO] =  1.93230919e+02 * (
			T * (  3.42336176e+00 + T * (  4.21988380e-03
			+ T * ( -1.43127321e-06 + T * (  1.96235356e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.55662337e+04 );
		cp[sHNCO] =  1.93230919e+02 * (
			 3.42336176e+00 + T * (  8.43976760e-03 
			+ T * ( -4.29381962e-06 + T * (  7.84941423e-10
			+ T *  0.00000000e+00 ) ) ) );
		h[sH2CN] =  2.96559424e+02 * (
			T * (  5.20970300e+00 + T * (  1.48464555e-03
			+ T * ( -9.51863033e-08 + T * ( -4.08887500e-11
			+ T *  6.08651780e-15 ) ) ) ) +  2.76771090e+04 );
		cp[sH2CN] =  2.96559424e+02 * (
			 5.20970300e+00 + T * (  2.96929110e-03 
			+ T * ( -2.85558910e-07 + T * ( -1.63555000e-10
			+ T *  3.04325890e-14 ) ) ) );
		h[sHCNN] =  2.02601004e+02 * (
			T * (  5.89463620e+00 + T * (  1.99479795e-03
			+ T * ( -5.32746000e-07 + T * (  7.31234875e-11
			+ T * -4.01893720e-15 ) ) ) ) +  5.34529410e+04 );
		cp[sHCNN] =  2.02601004e+02 * (
			 5.89463620e+00 + T * (  3.98959590e-03 
			+ T * ( -1.59823800e-06 + T * (  2.92493950e-10
			+ T * -2.00946860e-14 ) ) ) );
		h[sHCNO] =  1.93230919e+02 * (
			T * (  4.41622349e+00 + T * (  3.43873524e-03
			+ T * ( -1.09420184e-06 + T * (  1.41027230e-10
			+ T *  0.00000000e+00 ) ) ) ) +  1.88787247e+04 );
		cp[sHCNO] =  1.93230919e+02 * (
			 4.41622349e+00 + T * (  6.87747049e-03 
			+ T * ( -3.28260552e-06 + T * (  5.64108921e-10
			+ T *  0.00000000e+00 ) ) ) );
		h[sNH3] =  4.88102618e+02 * (
			T * (  2.63445210e+00 + T * (  2.83312800e-03
			+ T * ( -5.75955867e-07 + T * (  5.96679025e-11
			+ T * -2.51575720e-15 ) ) ) ) -6.54469580e+03 );
		cp[sNH3] =  4.88102618e+02 * (
			 2.63445210e+00 + T * (  5.66625600e-03 
			+ T * ( -1.72786760e-06 + T * (  2.38671610e-10
			+ T * -1.25787860e-14 ) ) ) );
		h[sC3H3] =  2.12893430e+02 * (
			T * (  8.83104700e+00 + T * (  2.17859700e-03
			+ T * ( -1.36968867e-07 + T * ( -5.92180750e-11
			+ T *  8.75304000e-15 ) ) ) ) +  3.76322300e+04 );
		cp[sC3H3] =  2.12893430e+02 * (
			 8.83104700e+00 + T * (  4.35719400e-03 
			+ T * ( -4.10906600e-07 + T * ( -2.36872300e-10
			+ T *  4.37652000e-14 ) ) ) );
		h[sC3H4] =  2.07536818e+02 * (
			T * (  5.72914400e+00 + T * (  6.18402250e-03
			+ T * ( -1.60187533e-06 + T * (  2.15034100e-10
			+ T * -1.16256040e-14 ) ) ) ) +  2.01298400e+04 );
		cp[sC3H4] =  2.07536818e+02 * (
			 5.72914400e+00 + T * (  1.23680450e-02 
			+ T * ( -4.80562600e-06 + T * (  8.60136400e-10
			+ T * -5.81280200e-14 ) ) ) );
		h[sC3H5] =  2.02443146e+02 * (
			T * (  9.65153900e+00 + T * (  4.03779800e-03
			+ T * ( -2.65514133e-07 + T * ( -1.16267400e-10
			+ T *  1.72065620e-14 ) ) ) ) +  1.53009550e+04 );
		cp[sC3H5] =  2.02443146e+02 * (
			 9.65153900e+00 + T * (  8.07559600e-03 
			+ T * ( -7.96542400e-07 + T * ( -4.65069600e-10
			+ T *  8.60328100e-14 ) ) ) );
		h[sC3H6] =  1.97593517e+02 * (
			T * (  6.73225700e+00 + T * (  7.45416800e-03
			+ T * ( -1.64996633e-06 + T * (  1.80300550e-10
			+ T * -7.53240800e-15 ) ) ) ) -9.23570300e+02 );
		cp[sC3H6] =  1.97593517e+02 * (
			 6.73225700e+00 + T * (  1.49083360e-02 
			+ T * ( -4.94989900e-06 + T * (  7.21202200e-10
			+ T * -3.76620400e-14 ) ) ) );
		h[sNXC3H7] =  1.92970803e+02 * (
			T * (  7.97829000e+00 + T * (  7.88056700e-03
			+ T * ( -1.72441433e-06 + T * (  1.86097300e-10
			+ T * -7.64995600e-15 ) ) ) ) +  8.26499100e+03 );
		cp[sNXC3H7] =  1.92970803e+02 * (
			 7.97829000e+00 + T * (  1.57611340e-02 
			+ T * ( -5.17324300e-06 + T * (  7.44389200e-10
			+ T * -3.82497800e-14 ) ) ) );
		h[sIXC3H7] =  1.92970803e+02 * (
			T * (  8.06336900e+00 + T * (  7.87243800e-03
			+ T * ( -1.72746367e-06 + T * (  1.86931125e-10
			+ T * -7.70884400e-15 ) ) ) ) +  7.35861400e+03 );
		cp[sIXC3H7] =  1.92970803e+02 * (
			 8.06336900e+00 + T * (  1.57448760e-02 
			+ T * ( -5.18239100e-06 + T * (  7.47724500e-10
			+ T * -3.85442200e-14 ) ) ) );
		h[sC3H8] =  1.88559441e+02 * (
			T * (  7.52521700e+00 + T * (  9.44517000e-03
			+ T * ( -2.09464133e-06 + T * (  2.29484325e-10
			+ T * -9.62482000e-15 ) ) ) ) -1.64645480e+04 );
		cp[sC3H8] =  1.88559441e+02 * (
			 7.52521700e+00 + T * (  1.88903400e-02 
			+ T * ( -6.28392400e-06 + T * (  9.17937300e-10
			+ T * -4.81241000e-14 ) ) ) );
	}
	else if (T >= 299.999999 ) {
		h[sN2] =  2.96728765e+02 * (
			T * (  3.29867700e+00 + T * (  7.04120200e-04
			+ T * ( -1.32107400e-06 + T * (  1.41037875e-09
			+ T * -4.88970800e-13 ) ) ) ) -1.02089990e+03 );
		cp[sN2] =  2.96728765e+02 * (
			 3.29867700e+00 + T * (  1.40824040e-03 
			+ T * ( -3.96322200e-06 + T * (  5.64151500e-09
			+ T * -2.44485400e-12 ) ) ) );
		h[sAR] =  2.08129068e+02 * (
			T * (  2.50000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) -7.45375000e+02 );
		cp[sAR] =  2.08129068e+02 * (
			 2.50000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sO] =  5.19646250e+02 * (
			T * (  3.16826710e+00 + T * ( -1.63965942e-03
			+ T * (  2.21435465e-06 + T * ( -1.53201656e-09
			+ T *  4.22531942e-13 ) ) ) ) +  2.91222592e+04 );
		cp[sO] =  5.19646250e+02 * (
			 3.16826710e+00 + T * ( -3.27931884e-03 
			+ T * (  6.64306396e-06 + T * ( -6.12806624e-09
			+ T *  2.11265971e-12 ) ) ) );
		h[sO2] =  2.59823125e+02 * (
			T * (  3.78245636e+00 + T * ( -1.49836708e-03
			+ T * (  3.28243400e-06 + T * ( -2.42032377e-09
			+ T *  6.48745674e-13 ) ) ) ) -1.06394356e+03 );
		cp[sO2] =  2.59823125e+02 * (
			 3.78245636e+00 + T * ( -2.99673416e-03 
			+ T * (  9.84730201e-06 + T * ( -9.68129509e-09
			+ T *  3.24372837e-12 ) ) ) );
		h[sH] =  8.24835317e+03 * (
			T * (  2.50000000e+00 + T * (  3.52666409e-13
			+ T * ( -6.65306547e-16 + T * (  5.75204080e-19
			+ T * -1.85546466e-22 ) ) ) ) +  2.54736599e+04 );
		cp[sH] =  8.24835317e+03 * (
			 2.50000000e+00 + T * (  7.05332819e-13 
			+ T * ( -1.99591964e-15 + T * (  2.30081632e-18
			+ T * -9.27732332e-22 ) ) ) );
		h[sOH] =  4.88848777e+02 * (
			T * (  3.99201543e+00 + T * ( -1.20065876e-03
			+ T * (  1.53931280e-06 + T * ( -9.70283332e-10
			+ T *  2.72822940e-13 ) ) ) ) +  3.61508056e+03 );
		cp[sOH] =  4.88848777e+02 * (
			 3.99201543e+00 + T * ( -2.40131752e-03 
			+ T * (  4.61793841e-06 + T * ( -3.88113333e-09
			+ T *  1.36411470e-12 ) ) ) );
		h[sH2] =  4.12417659e+03 * (
			T * (  2.34433112e+00 + T * (  3.99026037e-03
			+ T * ( -6.49271700e-06 + T * (  5.03930235e-09
			+ T * -1.47522352e-12 ) ) ) ) -9.17935173e+02 );
		cp[sH2] =  4.12417659e+03 * (
			 2.34433112e+00 + T * (  7.98052075e-03 
			+ T * ( -1.94781510e-05 + T * (  2.01572094e-08
			+ T * -7.37611761e-12 ) ) ) );
		h[sHO2] =  2.51888633e+02 * (
			T * (  4.30179801e+00 + T * ( -2.37456025e-03
			+ T * (  7.05276303e-06 + T * ( -6.06909735e-09
			+ T *  1.85845025e-12 ) ) ) ) +  2.94808040e+02 );
		cp[sHO2] =  2.51888633e+02 * (
			 4.30179801e+00 + T * ( -4.74912051e-03 
			+ T * (  2.11582891e-05 + T * ( -2.42763894e-08
			+ T *  9.29225124e-12 ) ) ) );
		h[sH2O2] =  2.44424389e+02 * (
			T * (  4.27611269e+00 + T * ( -2.71411208e-04
			+ T * (  5.57785670e-06 + T * ( -5.39427032e-09
			+ T *  1.72490873e-12 ) ) ) ) -1.77025821e+04 );
		cp[sH2O2] =  2.44424389e+02 * (
			 4.27611269e+00 + T * ( -5.42822417e-04 
			+ T * (  1.67335701e-05 + T * ( -2.15770813e-08
			+ T *  8.62454363e-12 ) ) ) );
		h[sCH] =  6.38680289e+02 * (
			T * (  3.48981665e+00 + T * (  1.61917771e-04
			+ T * ( -5.62996883e-07 + T * (  7.90543317e-10
			+ T * -2.81218134e-13 ) ) ) ) +  7.07972934e+04 );
		cp[sCH] =  6.38680289e+02 * (
			 3.48981665e+00 + T * (  3.23835541e-04 
			+ T * ( -1.68899065e-06 + T * (  3.16217327e-09
			+ T * -1.40609067e-12 ) ) ) );
		h[sCO] =  2.96834702e+02 * (
			T * (  3.57953347e+00 + T * ( -3.05176840e-04
			+ T * (  3.38938110e-07 + T * (  2.26751471e-10
			+ T * -1.80884900e-13 ) ) ) ) -1.43440860e+04 );
		cp[sCO] =  2.96834702e+02 * (
			 3.57953347e+00 + T * ( -6.10353680e-04 
			+ T * (  1.01681433e-06 + T * (  9.07005884e-10
			+ T * -9.04424499e-13 ) ) ) );
		h[sCH2] =  5.92780550e+02 * (
			T * (  3.76267867e+00 + T * (  4.84436072e-04
			+ T * (  9.31632803e-07 + T * ( -9.62727883e-10
			+ T *  3.37483438e-13 ) ) ) ) +  4.60040401e+04 );
		cp[sCH2] =  5.92780550e+02 * (
			 3.76267867e+00 + T * (  9.68872143e-04 
			+ T * (  2.79489841e-06 + T * ( -3.85091153e-09
			+ T *  1.68741719e-12 ) ) ) );
		h[sHCO] =  2.86523537e+02 * (
			T * (  4.22118584e+00 + T * ( -1.62196266e-03
			+ T * (  4.59331487e-06 + T * ( -3.32860233e-09
			+ T *  8.67537730e-13 ) ) ) ) +  3.83956496e+03 );
		cp[sHCO] =  2.86523537e+02 * (
			 4.22118584e+00 + T * ( -3.24392532e-03 
			+ T * (  1.37799446e-05 + T * ( -1.33144093e-08
			+ T *  4.33768865e-12 ) ) ) );
		h[sCH2GSGXCH2] =  5.92780550e+02 * (
			T * (  4.19860411e+00 + T * ( -1.18330710e-03
			+ T * (  2.74432073e-06 + T * ( -1.67203995e-09
			+ T *  3.88629474e-13 ) ) ) ) +  5.04968163e+04 );
		cp[sCH2GSGXCH2] =  5.92780550e+02 * (
			 4.19860411e+00 + T * ( -2.36661419e-03 
			+ T * (  8.23296220e-06 + T * ( -6.68815981e-09
			+ T *  1.94314737e-12 ) ) ) );
		h[sCH3] =  5.53035786e+02 * (
			T * (  3.67359040e+00 + T * (  1.00547588e-03
			+ T * (  1.91007285e-06 + T * ( -1.71779356e-09
			+ T *  5.08771468e-13 ) ) ) ) +  1.64449988e+04 );
		cp[sCH3] =  5.53035786e+02 * (
			 3.67359040e+00 + T * (  2.01095175e-03 
			+ T * (  5.73021856e-06 + T * ( -6.87117425e-09
			+ T *  2.54385734e-12 ) ) ) );
		h[sCH2O] =  2.76904683e+02 * (
			T * (  4.79372315e+00 + T * ( -4.95416684e-03
			+ T * (  1.24406669e-05 + T * ( -9.48213152e-09
			+ T *  2.63545304e-12 ) ) ) ) -1.43089567e+04 );
		cp[sCH2O] =  2.76904683e+02 * (
			 4.79372315e+00 + T * ( -9.90833369e-03 
			+ T * (  3.73220008e-05 + T * ( -3.79285261e-08
			+ T *  1.31772652e-11 ) ) ) );
		h[sCH4] =  5.18285750e+02 * (
			T * (  5.14987613e+00 + T * ( -6.83548940e-03
			+ T * (  1.63933533e-05 + T * ( -1.21185757e-08
			+ T *  3.33387912e-12 ) ) ) ) -1.02466476e+04 );
		cp[sCH4] =  5.18285750e+02 * (
			 5.14987613e+00 + T * ( -1.36709788e-02 
			+ T * (  4.91800599e-05 + T * ( -4.84743026e-08
			+ T *  1.66693956e-11 ) ) ) );
		h[sCO2] =  1.88919337e+02 * (
			T * (  2.35677352e+00 + T * (  4.49229839e-03
			+ T * ( -2.37452090e-06 + T * (  6.14797555e-10
			+ T * -2.87399096e-14 ) ) ) ) -4.83719697e+04 );
		cp[sCO2] =  1.88919337e+02 * (
			 2.35677352e+00 + T * (  8.98459677e-03 
			+ T * ( -7.12356269e-06 + T * (  2.45919022e-09
			+ T * -1.43699548e-13 ) ) ) );
		h[sCH2OH] =  2.67910679e+02 * (
			T * (  3.86388918e+00 + T * (  2.79836152e-03
			+ T * (  1.97757264e-06 + T * ( -2.61330030e-09
			+ T *  8.73934556e-13 ) ) ) ) -3.19391367e+03 );
		cp[sCH2OH] =  2.67910679e+02 * (
			 3.86388918e+00 + T * (  5.59672304e-03 
			+ T * (  5.93271791e-06 + T * ( -1.04532012e-08
			+ T *  4.36967278e-12 ) ) ) );
		h[sCH3O] =  2.67910679e+02 * (
			T * (  2.10620400e+00 + T * (  3.60829750e-03
			+ T * (  1.77949067e-06 + T * ( -1.84440900e-09
			+ T *  4.15122000e-13 ) ) ) ) +  9.78601100e+02 );
		cp[sCH3O] =  2.67910679e+02 * (
			 2.10620400e+00 + T * (  7.21659500e-03 
			+ T * (  5.33847200e-06 + T * ( -7.37763600e-09
			+ T *  2.07561000e-12 ) ) ) );
		h[sCH3OH] =  2.59482554e+02 * (
			T * (  5.71539582e+00 + T * ( -7.61545645e-03
			+ T * (  2.17480385e-05 + T * ( -1.77701722e-08
			+ T *  5.22705396e-12 ) ) ) ) -2.56427656e+04 );
		cp[sCH3OH] =  2.59482554e+02 * (
			 5.71539582e+00 + T * ( -1.52309129e-02 
			+ T * (  6.52441155e-05 + T * ( -7.10806889e-08
			+ T *  2.61352698e-11 ) ) ) );
		h[sC2H] =  3.32201534e+02 * (
			T * (  2.88965733e+00 + T * (  6.70498055e-03
			+ T * ( -9.49231670e-06 + T * (  7.36977613e-09
			+ T * -2.18663022e-12 ) ) ) ) +  6.68393932e+04 );
		cp[sC2H] =  3.32201534e+02 * (
			 2.88965733e+00 + T * (  1.34099611e-02 
			+ T * ( -2.84769501e-05 + T * (  2.94791045e-08
			+ T * -1.09331511e-11 ) ) ) );
		h[sC2H2] =  3.19340144e+02 * (
			T * (  8.08681094e-01 + T * (  1.16807815e-02
			+ T * ( -1.18390605e-05 + T * (  7.00381092e-09
			+ T * -1.70014595e-12 ) ) ) ) +  2.64289807e+04 );
		cp[sC2H2] =  3.19340144e+02 * (
			 8.08681094e-01 + T * (  2.33615629e-02 
			+ T * ( -3.55171815e-05 + T * (  2.80152437e-08
			+ T * -8.50072974e-12 ) ) ) );
		h[sHCCO] =  2.02650385e+02 * (
			T * (  2.25172140e+00 + T * (  8.82751050e-03
			+ T * ( -7.90970033e-06 + T * (  4.31893975e-09
			+ T * -1.01329622e-12 ) ) ) ) +  2.00594490e+04 );
		cp[sHCCO] =  2.02650385e+02 * (
			 2.25172140e+00 + T * (  1.76550210e-02 
			+ T * ( -2.37291010e-05 + T * (  1.72757590e-08
			+ T * -5.06648110e-12 ) ) ) );
		h[sC2H3] =  3.07437509e+02 * (
			T * (  3.21246645e+00 + T * (  7.57395810e-04
			+ T * (  8.64031373e-06 + T * ( -8.94144617e-09
			+ T *  2.94301746e-12 ) ) ) ) +  3.48598468e+04 );
		cp[sC2H3] =  3.07437509e+02 * (
			 3.21246645e+00 + T * (  1.51479162e-03 
			+ T * (  2.59209412e-05 + T * ( -3.57657847e-08
			+ T *  1.47150873e-11 ) ) ) );
		h[sCH2CO] =  1.97790941e+02 * (
			T * (  2.13583630e+00 + T * (  9.05943605e-03
			+ T * ( -5.79824913e-06 + T * (  2.33599392e-09
			+ T * -4.02915230e-13 ) ) ) ) -7.04291804e+03 );
		cp[sCH2CO] =  1.97790941e+02 * (
			 2.13583630e+00 + T * (  1.81188721e-02 
			+ T * ( -1.73947474e-05 + T * (  9.34397568e-09
			+ T * -2.01457615e-12 ) ) ) );
		h[sC2H4] =  2.96390275e+02 * (
			T * (  3.95920148e+00 + T * ( -3.78526124e-03
			+ T * (  1.90330097e-05 + T * ( -1.72897188e-08
			+ T *  5.39768746e-12 ) ) ) ) +  5.08977593e+03 );
		cp[sC2H4] =  2.96390275e+02 * (
			 3.95920148e+00 + T * ( -7.57052247e-03 
			+ T * (  5.70990292e-05 + T * ( -6.91588753e-08
			+ T *  2.69884373e-11 ) ) ) );
		h[sC2H5] =  2.86109429e+02 * (
			T * (  4.30646568e+00 + T * ( -2.09329446e-03
			+ T * (  1.65714269e-05 + T * ( -1.49781651e-08
			+ T *  4.61018008e-12 ) ) ) ) +  1.28416265e+04 );
		cp[sC2H5] =  2.86109429e+02 * (
			 4.30646568e+00 + T * ( -4.18658892e-03 
			+ T * (  4.97142807e-05 + T * ( -5.99126606e-08
			+ T *  2.30509004e-11 ) ) ) );
		h[sC2H6] =  2.76517893e+02 * (
			T * (  4.29142492e+00 + T * ( -2.75077135e-03
			+ T * (  1.99812763e-05 + T * ( -1.77116571e-08
			+ T *  5.37371542e-12 ) ) ) ) -1.15222055e+04 );
		cp[sC2H6] =  2.76517893e+02 * (
			 4.29142492e+00 + T * ( -5.50154270e-03 
			+ T * (  5.99438288e-05 + T * ( -7.08466285e-08
			+ T *  2.68685771e-11 ) ) ) );
		h[sH2O] =  4.61497558e+02 * (
			T * (  4.19864056e+00 + T * ( -1.01821705e-03
			+ T * (  2.17346737e-06 + T * ( -1.37199266e-09
			+ T *  3.54395634e-13 ) ) ) ) -3.02937267e+04 );
		cp[sH2O] =  4.61497558e+02 * (
			 4.19864056e+00 + T * ( -2.03643410e-03 
			+ T * (  6.52040211e-06 + T * ( -5.48797062e-09
			+ T *  1.77197817e-12 ) ) ) );
		h[sC] =  6.92284763e+02 * (
			T * (  2.55423955e+00 + T * ( -1.60768862e-04
			+ T * (  2.44597415e-07 + T * ( -1.83058722e-10
			+ T *  5.33042892e-14 ) ) ) ) +  8.54438832e+04 );
		cp[sC] =  6.92284763e+02 * (
			 2.55423955e+00 + T * ( -3.21537724e-04 
			+ T * (  7.33792245e-07 + T * ( -7.32234889e-10
			+ T *  2.66521446e-13 ) ) ) );
		h[sHCCOH] =  1.97790941e+02 * (
			T * (  1.24237330e+00 + T * (  1.55361005e-02
			+ T * ( -1.69556213e-05 + T * (  1.07842828e-08
			+ T * -2.80291880e-12 ) ) ) ) +  8.03161430e+03 );
		cp[sHCCOH] =  1.97790941e+02 * (
			 1.24237330e+00 + T * (  3.10722010e-02 
			+ T * ( -5.08668640e-05 + T * (  4.31371310e-08
			+ T * -1.40145940e-11 ) ) ) );
		h[sN] =  5.93457530e+02 * (
			T * (  2.50000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  5.61046370e+04 );
		cp[sN] =  5.93457530e+02 * (
			 2.50000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sNO] =  2.77052316e+02 * (
			T * (  4.21847630e+00 + T * ( -2.31948800e-03
			+ T * (  3.68034067e-06 + T * ( -2.33403385e-09
			+ T *  5.60715400e-13 ) ) ) ) +  9.84462300e+03 );
		cp[sNO] =  2.77052316e+02 * (
			 4.21847630e+00 + T * ( -4.63897600e-03 
			+ T * (  1.10410220e-05 + T * ( -9.33613540e-09
			+ T *  2.80357700e-12 ) ) ) );
		h[sN2O] =  1.88876420e+02 * (
			T * (  2.25715020e+00 + T * (  5.65236400e-03
			+ T * ( -4.55710633e-06 + T * (  2.42049515e-09
			+ T * -5.86143640e-13 ) ) ) ) +  8.74177440e+03 );
		cp[sN2O] =  1.88876420e+02 * (
			 2.25715020e+00 + T * (  1.13047280e-02 
			+ T * ( -1.36713190e-05 + T * (  9.68198060e-09
			+ T * -2.93071820e-12 ) ) ) );
		h[sNO2] =  1.80707238e+02 * (
			T * (  3.94403120e+00 + T * ( -7.92714500e-04
			+ T * (  5.55260400e-06 + T * ( -5.11885650e-09
			+ T *  1.56701128e-12 ) ) ) ) +  2.89661790e+03 );
		cp[sNO2] =  1.80707238e+02 * (
			 3.94403120e+00 + T * ( -1.58542900e-03 
			+ T * (  1.66578120e-05 + T * ( -2.04754260e-08
			+ T *  7.83505640e-12 ) ) ) );
		h[sNH] =  5.53624983e+02 * (
			T * (  3.49290850e+00 + T * (  1.55895990e-04
			+ T * ( -4.96349467e-07 + T * (  6.20411050e-10
			+ T * -2.07139340e-13 ) ) ) ) +  4.18806290e+04 );
		cp[sNH] =  5.53624983e+02 * (
			 3.49290850e+00 + T * (  3.11791980e-04 
			+ T * ( -1.48904840e-06 + T * (  2.48164420e-09
			+ T * -1.03569670e-12 ) ) ) );
		h[sHNO] =  2.68048875e+02 * (
			T * (  4.53349160e+00 + T * ( -2.83480855e-03
			+ T * (  6.15773567e-06 + T * ( -4.28427350e-09
			+ T *  1.10909146e-12 ) ) ) ) +  1.15482970e+04 );
		cp[sHNO] =  2.68048875e+02 * (
			 4.53349160e+00 + T * ( -5.66961710e-03 
			+ T * (  1.84732070e-05 + T * ( -1.71370940e-08
			+ T *  5.54545730e-12 ) ) ) );
		h[sNH2] =  5.18803195e+02 * (
			T * (  4.20400290e+00 + T * ( -1.05306925e-03
			+ T * (  2.36894493e-06 + T * ( -1.40287992e-09
			+ T *  3.28814340e-13 ) ) ) ) +  2.18859100e+04 );
		cp[sNH2] =  5.18803195e+02 * (
			 4.20400290e+00 + T * ( -2.10613850e-03 
			+ T * (  7.10683480e-06 + T * ( -5.61151970e-09
			+ T *  1.64407170e-12 ) ) ) );
		h[sNNH] =  2.86424831e+02 * (
			T * (  4.34469270e+00 + T * ( -2.42485360e-03
			+ T * (  6.68648633e-06 + T * ( -5.43161600e-09
			+ T *  1.58939078e-12 ) ) ) ) +  2.87919730e+04 );
		cp[sNNH] =  2.86424831e+02 * (
			 4.34469270e+00 + T * ( -4.84970720e-03 
			+ T * (  2.00594590e-05 + T * ( -2.17264640e-08
			+ T *  7.94695390e-12 ) ) ) );
		h[sCN] =  3.19536510e+02 * (
			T * (  3.61293510e+00 + T * ( -4.77756635e-04
			+ T * (  7.14765900e-07 + T * ( -7.87908075e-11
			+ T * -9.28607120e-14 ) ) ) ) +  5.17083400e+04 );
		cp[sCN] =  3.19536510e+02 * (
			 3.61293510e+00 + T * ( -9.55513270e-04 
			+ T * (  2.14429770e-06 + T * ( -3.15163230e-10
			+ T * -4.64303560e-13 ) ) ) );
		h[sNCO] =  1.97866254e+02 * (
			T * (  2.82693080e+00 + T * (  4.40258440e-03
			+ T * ( -2.79553780e-06 + T * (  1.20042410e-09
			+ T * -2.66271900e-13 ) ) ) ) +  1.46824770e+04 );
		cp[sNCO] =  1.97866254e+02 * (
			 2.82693080e+00 + T * (  8.80516880e-03 
			+ T * ( -8.38661340e-06 + T * (  4.80169640e-09
			+ T * -1.33135950e-12 ) ) ) );
		h[sHCN] =  3.07619506e+02 * (
			T * (  2.25898860e+00 + T * (  5.02558500e-03
			+ T * ( -4.45058767e-06 + T * (  2.52308725e-09
			+ T * -6.01780560e-13 ) ) ) ) +  1.52158530e+04 );
		cp[sHCN] =  3.07619506e+02 * (
			 2.25898860e+00 + T * (  1.00511700e-02 
			+ T * ( -1.33517630e-05 + T * (  1.00923490e-08
			+ T * -3.00890280e-12 ) ) ) );
		h[sHOCN] =  1.93230919e+02 * (
			T * (  3.78604952e+00 + T * (  3.44333961e-03
			+ T * ( -1.07162621e-06 + T * (  1.29298942e-10
			+ T *  2.38721576e-15 ) ) ) ) -2.82698400e+03 );
		cp[sHOCN] =  1.93230919e+02 * (
			 3.78604952e+00 + T * (  6.88667922e-03 
			+ T * ( -3.21487864e-06 + T * (  5.17195767e-10
			+ T *  1.19360788e-14 ) ) ) );
		h[sHNCO] =  1.93230919e+02 * (
			T * (  3.63096317e+00 + T * (  3.65141179e-03
			+ T * ( -7.60166677e-07 + T * ( -1.65317825e-10
			+ T *  7.24471504e-14 ) ) ) ) -1.55873636e+04 );
		cp[sHNCO] =  1.93230919e+02 * (
			 3.63096317e+00 + T * (  7.30282357e-03 
			+ T * ( -2.28050003e-06 + T * ( -6.61271298e-10
			+ T *  3.62235752e-13 ) ) ) );
		h[sH2CN] =  2.96559424e+02 * (
			T * (  2.85166100e+00 + T * (  2.84761655e-03
			+ T * (  3.57046667e-07 + T * ( -4.05653000e-10
			+ T * -4.70221620e-14 ) ) ) ) +  2.86378200e+04 );
		cp[sH2CN] =  2.96559424e+02 * (
			 2.85166100e+00 + T * (  5.69523310e-03 
			+ T * (  1.07114000e-06 + T * ( -1.62261200e-09
			+ T * -2.35110810e-13 ) ) ) );
		h[sHCNN] =  2.02601004e+02 * (
			T * (  2.52431940e+00 + T * (  7.98030950e-03
			+ T * ( -6.27211800e-06 + T * (  3.03138500e-09
			+ T * -6.47147560e-13 ) ) ) ) +  5.42619840e+04 );
		cp[sHCNN] =  2.02601004e+02 * (
			 2.52431940e+00 + T * (  1.59606190e-02 
			+ T * ( -1.88163540e-05 + T * (  1.21255400e-08
			+ T * -3.23573780e-12 ) ) ) );
		h[sHCNO] =  1.93230919e+02 * (
			T * (  2.64727989e+00 + T * (  6.37526710e-03
			+ T * ( -3.49314120e-06 + T * (  1.10358209e-09
			+ T * -1.51504293e-13 ) ) ) ) +  1.92990252e+04 );
		cp[sHCNO] =  1.93230919e+02 * (
			 2.64727989e+00 + T * (  1.27505342e-02 
			+ T * ( -1.04794236e-05 + T * (  4.41432836e-09
			+ T * -7.57521466e-13 ) ) ) );
		h[sNH3] =  4.88102618e+02 * (
			T * (  4.28602740e+00 + T * ( -2.33026150e-03
			+ T * (  7.23950433e-06 + T * ( -5.70222175e-09
			+ T *  1.65276092e-12 ) ) ) ) -6.74172850e+03 );
		cp[sNH3] =  4.88102618e+02 * (
			 4.28602740e+00 + T * ( -4.66052300e-03 
			+ T * (  2.17185130e-05 + T * ( -2.28088870e-08
			+ T *  8.26380460e-12 ) ) ) );
		h[sC3H3] =  2.12893430e+02 * (
			T * (  4.75419900e+00 + T * (  5.54013850e-03
			+ T * (  9.31107667e-08 + T * ( -1.36980300e-09
			+ T *  3.89925800e-13 ) ) ) ) +  3.90468300e+04 );
		cp[sC3H3] =  2.12893430e+02 * (
			 4.75419900e+00 + T * (  1.10802770e-02 
			+ T * (  2.79332300e-07 + T * ( -5.47921200e-09
			+ T *  1.94962900e-12 ) ) ) );
		h[sC3H4] =  2.07536818e+02 * (
			T * ( -2.13196800e-01 + T * (  1.67935650e-02
			+ T * ( -1.26829000e-05 + T * (  6.86459500e-09
			+ T * -1.73800880e-12 ) ) ) ) +  2.16204800e+04 );
		cp[sC3H4] =  2.07536818e+02 * (
			-2.13196800e-01 + T * (  3.35871300e-02 
			+ T * ( -3.80487000e-05 + T * (  2.74583800e-08
			+ T * -8.69004400e-12 ) ) ) );
		h[sC3H5] =  2.02443146e+02 * (
			T * (  2.27648600e+00 + T * (  9.92782000e-03
			+ T * (  3.74614033e-07 + T * ( -2.53643925e-09
			+ T *  6.88268400e-13 ) ) ) ) +  1.78949600e+04 );
		cp[sC3H5] =  2.02443146e+02 * (
			 2.27648600e+00 + T * (  1.98556400e-02 
			+ T * (  1.12384210e-06 + T * ( -1.01457570e-08
			+ T *  3.44134200e-12 ) ) ) );
		h[sC3H6] =  1.97593517e+02 * (
			T * (  1.49330710e+00 + T * (  1.04625850e-02
			+ T * (  1.49559800e-06 + T * ( -4.17228025e-09
			+ T *  1.43162920e-12 ) ) ) ) +  1.07482640e+03 );
		cp[sC3H6] =  1.97593517e+02 * (
			 1.49330710e+00 + T * (  2.09251700e-02 
			+ T * (  4.48679400e-06 + T * ( -1.66891210e-08
			+ T *  7.15814600e-12 ) ) ) );
		h[sNXC3H7] =  1.92970803e+02 * (
			T * (  1.92253600e+00 + T * (  1.23946350e-02
			+ T * (  6.03416333e-07 + T * ( -4.45816250e-09
			+ T *  1.71659920e-12 ) ) ) ) +  1.03988710e+04 );
		cp[sNXC3H7] =  1.92970803e+02 * (
			 1.92253600e+00 + T * (  2.47892700e-02 
			+ T * (  1.81024900e-06 + T * ( -1.78326500e-08
			+ T *  8.58299600e-12 ) ) ) );
		h[sIXC3H7] =  1.92970803e+02 * (
			T * (  1.71329900e+00 + T * (  1.27130800e-02
			+ T * (  5.26936100e-07 + T * ( -4.55321500e-09
			+ T *  1.76554200e-12 ) ) ) ) +  9.58055400e+03 );
		cp[sIXC3H7] =  1.92970803e+02 * (
			 1.71329900e+00 + T * (  2.54261600e-02 
			+ T * (  1.58080830e-06 + T * ( -1.82128600e-08
			+ T *  8.82771000e-12 ) ) ) );
		h[sC3H8] =  1.88559441e+02 * (
			T * (  8.96920800e-01 + T * (  1.33449300e-02
			+ T * (  1.81047500e-06 + T * ( -5.31500000e-09
			+ T *  1.84866600e-12 ) ) ) ) -1.39549180e+04 );
		cp[sC3H8] =  1.88559441e+02 * (
			 8.96920800e-01 + T * (  2.66898600e-02 
			+ T * (  5.43142500e-06 + T * ( -2.12600000e-08
			+ T *  9.24333000e-12 ) ) ) );
	}
	else {
		ComputeThermoData( h, cp, 300.0 );
		for (i = 0; i < sEnd; i++) {
			h[i] = (T-300.)*cp[i] + h[i];
		}
	}
}
int GetNSpecies( void )
{
	return sEnd;
}

int GetNSpecs( void )
{
	return sEnd;
}

int GetNReactions( void )
{
	return rEnd;
}

double MAX_C(double X1, double X2)
{
	return ( (X1 > X2) ? X1 : X2 );
}

