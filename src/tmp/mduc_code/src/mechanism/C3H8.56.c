

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "C3H8.56.h"

/* #include "functions.h" */

static double GetLindRateCoeff( double temp, double pressure
			, double k0, double kInf, double fc, double conc );

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
	double  RGAS    = 8314.340;
	double	lgt;
	double	rt;
        double k9f0, k9fInf, fc9, k9b0, k9bInf;
        double k16f0, k16fInf, fc16, k16b0, k16bInf;
        double k63f0, k63fInf, fc63, k63b0, k63bInf;
        double k86f0, k86fInf, fc86, k86b0, k86bInf;
        double k98f0, k98fInf, fc98, k98b0, k98bInf;
        double k113f0, k113fInf, fc113, k113b0, k113bInf;
        double k144f0, k144fInf, fc144, k144b0, k144bInf;
        double k145f0, k145fInf, fc145, k145b0, k145bInf;
        double k156f0, k156fInf, fc156, k156b0, k156bInf;
        double k204f0, k204fInf, fc204, k204b0, k204bInf;
        double k233f0, k233fInf, fc233, k233b0, k233bInf;
        double k400f0, k400fInf, fc400, k400b0, k400bInf;
	
	rt = RGAS * temp;
	lgt = log( temp );
	  
	M[mM1] = c[sN2] + c[sAR] + c[sH] + c[sO2]
		 + c[sO] + c[sOH] + 2.5 * c[sH2] + 12 * c[sH2O]
		 + c[sHO2] + c[sH2O2] + 1.9 * c[sCO] + 3.8 * c[sCO2]
		 + c[sHCO] + c[sCH3] + c[sCH4] + c[sCH2O]
		 + c[sHOCHO] + c[sCH3O] + c[sCH3OH] + c[sCH3O2]
		 + c[sCH3O2H] + c[sC2H6] + c[sC2H5] + c[sC2H4]
		 + c[sCH3CHO] + c[sCH3CO] + c[sCH3CO3] + c[sCH3CO3H]
		 + c[sCH2CHO] + c[sCH3CO2] + c[sC2H3] + c[sC2H3CHO]
		 + c[sC2H3CO] + c[sC3H8] + c[sNXC3H7] + c[sIXC3H7]
		 + c[sC3H5XAXC3H5] + c[sC3H6] + c[sIXC3H7O2] + c[sIXC3H7O2H]
		 + c[sC3H5XTXC3H5] + c[sC3H5O] + c[sCH3COCH2] + c[sIXC3H7O]
		 + c[sNXC3H7O2] + c[sC3H6OOH1X2XC3H7O2] + c[sC3H6OOH1X3XC3H7O2] + c[sC3H6OOH2X1XC3H7O2]
		 + c[sC3H6OOH1X2O2XC3H7O4] + c[sC3H6OOH1X3O2XC3H7O4] + c[sC3H6OOH2X1O2XC3H7O4] + c[sC3KET12XC3H6O3]
		 + c[sC3KET13XC3H6O3] + c[sC3KET21XC3H6O3] + c[sC3H52X1J3OOHXC3H7O4] + c[sAXC3H5OOH];

	M[mM2] = c[sN2] + c[sAR] + c[sH] + c[sO2]
		 + c[sO] + c[sOH] + 2.5 * c[sH2] + 12 * c[sH2O]
		 + c[sHO2] + c[sH2O2] + 1.9 * c[sCO] + 3.8 * c[sCO2]
		 + c[sHCO] + c[sCH3] + 2 * c[sCH4] + c[sCH2O]
		 + c[sHOCHO] + c[sCH3O] + c[sCH3OH] + c[sCH3O2]
		 + c[sCH3O2H] + 3 * c[sC2H6] + c[sC2H5] + c[sC2H4]
		 + c[sCH3CHO] + c[sCH3CO] + c[sCH3CO3] + c[sCH3CO3H]
		 + c[sCH2CHO] + c[sCH3CO2] + c[sC2H3] + c[sC2H3CHO]
		 + c[sC2H3CO] + c[sC3H8] + c[sNXC3H7] + c[sIXC3H7]
		 + c[sC3H5XAXC3H5] + c[sC3H6] + c[sIXC3H7O2] + c[sIXC3H7O2H]
		 + c[sC3H5XTXC3H5] + c[sC3H5O] + c[sCH3COCH2] + c[sIXC3H7O]
		 + c[sNXC3H7O2] + c[sC3H6OOH1X2XC3H7O2] + c[sC3H6OOH1X3XC3H7O2] + c[sC3H6OOH2X1XC3H7O2]
		 + c[sC3H6OOH1X2O2XC3H7O4] + c[sC3H6OOH1X3O2XC3H7O4] + c[sC3H6OOH2X1O2XC3H7O4] + c[sC3KET12XC3H6O3]
		 + c[sC3KET13XC3H6O3] + c[sC3KET21XC3H6O3] + c[sC3H52X1J3OOHXC3H7O4] + c[sAXC3H5OOH];

	M[mM3] = c[sN2] + c[sAR] + c[sH] + c[sO2]
		 + c[sO] + c[sOH] + 2.5 * c[sH2] + 12 * c[sH2O]
		 + c[sHO2] + c[sH2O2] + 1.5 * c[sCO] + 2 * c[sCO2]
		 + c[sHCO] + c[sCH3] + 2 * c[sCH4] + c[sCH2O]
		 + c[sHOCHO] + c[sCH3O] + c[sCH3OH] + c[sCH3O2]
		 + c[sCH3O2H] + 3 * c[sC2H6] + c[sC2H5] + c[sC2H4]
		 + c[sCH3CHO] + c[sCH3CO] + c[sCH3CO3] + c[sCH3CO3H]
		 + c[sCH2CHO] + c[sCH3CO2] + c[sC2H3] + c[sC2H3CHO]
		 + c[sC2H3CO] + c[sC3H8] + c[sNXC3H7] + c[sIXC3H7]
		 + c[sC3H5XAXC3H5] + c[sC3H6] + c[sIXC3H7O2] + c[sIXC3H7O2H]
		 + c[sC3H5XTXC3H5] + c[sC3H5O] + c[sCH3COCH2] + c[sIXC3H7O]
		 + c[sNXC3H7O2] + c[sC3H6OOH1X2XC3H7O2] + c[sC3H6OOH1X3XC3H7O2] + c[sC3H6OOH2X1XC3H7O2]
		 + c[sC3H6OOH1X2O2XC3H7O4] + c[sC3H6OOH1X3O2XC3H7O4] + c[sC3H6OOH2X1O2XC3H7O4] + c[sC3KET12XC3H6O3]
		 + c[sC3KET13XC3H6O3] + c[sC3KET21XC3H6O3] + c[sC3H52X1J3OOHXC3H7O4] + c[sAXC3H5OOH];

	M[mM4] = c[sN2] + c[sAR] + c[sH] + c[sO2]
		 + c[sO] + c[sOH] + 0.73 * c[sH2] + 12 * c[sH2O]
		 + c[sHO2] + c[sH2O2] + c[sCO] + c[sCO2]
		 + c[sHCO] + c[sCH3] + 2 * c[sCH4] + c[sCH2O]
		 + c[sHOCHO] + c[sCH3O] + c[sCH3OH] + c[sCH3O2]
		 + c[sCH3O2H] + 3 * c[sC2H6] + c[sC2H5] + c[sC2H4]
		 + c[sCH3CHO] + c[sCH3CO] + c[sCH3CO3] + c[sCH3CO3H]
		 + c[sCH2CHO] + c[sCH3CO2] + c[sC2H3] + c[sC2H3CHO]
		 + c[sC2H3CO] + c[sC3H8] + c[sNXC3H7] + c[sIXC3H7]
		 + c[sC3H5XAXC3H5] + c[sC3H6] + c[sIXC3H7O2] + c[sIXC3H7O2H]
		 + c[sC3H5XTXC3H5] + c[sC3H5O] + c[sCH3COCH2] + c[sIXC3H7O]
		 + c[sNXC3H7O2] + c[sC3H6OOH1X2XC3H7O2] + c[sC3H6OOH1X3XC3H7O2] + c[sC3H6OOH2X1XC3H7O2]
		 + c[sC3H6OOH1X2O2XC3H7O4] + c[sC3H6OOH1X3O2XC3H7O4] + c[sC3H6OOH2X1O2XC3H7O4] + c[sC3KET12XC3H6O3]
		 + c[sC3KET13XC3H6O3] + c[sC3KET21XC3H6O3] + c[sC3H52X1J3OOHXC3H7O4] + c[sAXC3H5OOH];

	M[mM5] = c[sN2] + c[sAR] + c[sH] + c[sO2]
		 + c[sO] + c[sOH] + 1.3 * c[sH2] + 14 * c[sH2O]
		 + c[sHO2] + c[sH2O2] + 1.9 * c[sCO] + 3.8 * c[sCO2]
		 + c[sHCO] + c[sCH3] + 2 * c[sCH4] + c[sCH2O]
		 + c[sHOCHO] + c[sCH3O] + c[sCH3OH] + c[sCH3O2]
		 + c[sCH3O2H] + 3 * c[sC2H6] + c[sC2H5] + c[sC2H4]
		 + c[sCH3CHO] + c[sCH3CO] + c[sCH3CO3] + c[sCH3CO3H]
		 + c[sCH2CHO] + c[sCH3CO2] + c[sC2H3] + c[sC2H3CHO]
		 + c[sC2H3CO] + c[sC3H8] + c[sNXC3H7] + c[sIXC3H7]
		 + c[sC3H5XAXC3H5] + c[sC3H6] + c[sIXC3H7O2] + c[sIXC3H7O2H]
		 + c[sC3H5XTXC3H5] + c[sC3H5O] + c[sCH3COCH2] + c[sIXC3H7O]
		 + c[sNXC3H7O2] + c[sC3H6OOH1X2XC3H7O2] + c[sC3H6OOH1X3XC3H7O2] + c[sC3H6OOH2X1XC3H7O2]
		 + c[sC3H6OOH1X2O2XC3H7O4] + c[sC3H6OOH1X3O2XC3H7O4] + c[sC3H6OOH2X1O2XC3H7O4] + c[sC3KET12XC3H6O3]
		 + c[sC3KET13XC3H6O3] + c[sC3KET21XC3H6O3] + c[sC3H52X1J3OOHXC3H7O4] + c[sAXC3H5OOH];

	M[mM6] = c[sN2] + c[sAR] + c[sH] + c[sO2]
		 + c[sO] + c[sOH] + 2.5 * c[sH2] + 12 * c[sH2O]
		 + c[sHO2] + c[sH2O2] + 1.9 * c[sCO] + 3.8 * c[sCO2]
		 + c[sHCO] + c[sCH3] + 2 * c[sCH4] + c[sCH2O]
		 + c[sHOCHO] + c[sCH3O] + c[sCH3OH] + c[sCH3O2]
		 + c[sCH3O2H] + 3 * c[sC2H6] + c[sC2H5] + c[sC2H4]
		 + c[sCH3CHO] + c[sCH3CO] + c[sCH3CO3] + c[sCH3CO3H]
		 + c[sCH2CHO] + c[sCH3CO2] + c[sC2H3] + c[sC2H3CHO]
		 + c[sC2H3CO] + c[sC3H8] + c[sNXC3H7] + c[sIXC3H7]
		 + c[sC3H5XAXC3H5] + c[sC3H6] + c[sIXC3H7O2] + c[sIXC3H7O2H]
		 + c[sC3H5XTXC3H5] + c[sC3H5O] + c[sCH3COCH2] + c[sIXC3H7O]
		 + c[sNXC3H7O2] + c[sC3H6OOH1X2XC3H7O2] + c[sC3H6OOH1X3XC3H7O2] + c[sC3H6OOH2X1XC3H7O2]
		 + c[sC3H6OOH1X2O2XC3H7O4] + c[sC3H6OOH1X3O2XC3H7O4] + c[sC3H6OOH2X1O2XC3H7O4] + c[sC3KET12XC3H6O3]
		 + c[sC3KET13XC3H6O3] + c[sC3KET21XC3H6O3] + c[sC3H52X1J3OOHXC3H7O4] + c[sAXC3H5OOH];

	M[mM8] = c[sN2] + c[sAR] + c[sH] + 12 * c[sO2]
		 + c[sO] + c[sOH] + 2 * c[sH2] + 6 * c[sH2O]
		 + c[sHO2] + c[sH2O2] + 1.5 * c[sCO] + 2 * c[sCO2]
		 + c[sHCO] + c[sCH3] + 2 * c[sCH4] + c[sCH2O]
		 + c[sHOCHO] + c[sCH3O] + c[sCH3OH] + c[sCH3O2]
		 + c[sCH3O2H] + 3 * c[sC2H6] + c[sC2H5] + c[sC2H4]
		 + c[sCH3CHO] + c[sCH3CO] + c[sCH3CO3] + c[sCH3CO3H]
		 + c[sCH2CHO] + c[sCH3CO2] + c[sC2H3] + c[sC2H3CHO]
		 + c[sC2H3CO] + c[sC3H8] + c[sNXC3H7] + c[sIXC3H7]
		 + c[sC3H5XAXC3H5] + c[sC3H6] + c[sIXC3H7O2] + c[sIXC3H7O2H]
		 + c[sC3H5XTXC3H5] + c[sC3H5O] + c[sCH3COCH2] + c[sIXC3H7O]
		 + c[sNXC3H7O2] + c[sC3H6OOH1X2XC3H7O2] + c[sC3H6OOH1X3XC3H7O2] + c[sC3H6OOH2X1XC3H7O2]
		 + c[sC3H6OOH1X2O2XC3H7O4] + c[sC3H6OOH1X3O2XC3H7O4] + c[sC3H6OOH2X1O2XC3H7O4] + c[sC3KET12XC3H6O3]
		 + c[sC3KET13XC3H6O3] + c[sC3KET21XC3H6O3] + c[sC3H52X1J3OOHXC3H7O4] + c[sAXC3H5OOH];

	M[mM0] = c[sN2] + c[sAR] + c[sH] + c[sO2]
		 + c[sO] + c[sOH] + c[sH2] + c[sH2O]
		 + c[sHO2] + c[sH2O2] + c[sCO] + c[sCO2]
		 + c[sHCO] + c[sCH3] + c[sCH4] + c[sCH2O]
		 + c[sHOCHO] + c[sCH3O] + c[sCH3OH] + c[sCH3O2]
		 + c[sCH3O2H] + c[sC2H6] + c[sC2H5] + c[sC2H4]
		 + c[sCH3CHO] + c[sCH3CO] + c[sCH3CO3] + c[sCH3CO3H]
		 + c[sCH2CHO] + c[sCH3CO2] + c[sC2H3] + c[sC2H3CHO]
		 + c[sC2H3CO] + c[sC3H8] + c[sNXC3H7] + c[sIXC3H7]
		 + c[sC3H5XAXC3H5] + c[sC3H6] + c[sIXC3H7O2] + c[sIXC3H7O2H]
		 + c[sC3H5XTXC3H5] + c[sC3H5O] + c[sCH3COCH2] + c[sIXC3H7O]
		 + c[sNXC3H7O2] + c[sC3H6OOH1X2XC3H7O2] + c[sC3H6OOH1X3XC3H7O2] + c[sC3H6OOH2X1XC3H7O2]
		 + c[sC3H6OOH1X2O2XC3H7O4] + c[sC3H6OOH1X3O2XC3H7O4] + c[sC3H6OOH2X1O2XC3H7O4] + c[sC3KET12XC3H6O3]
		 + c[sC3KET13XC3H6O3] + c[sC3KET21XC3H6O3] + c[sC3H52X1J3OOHXC3H7O4] + c[sAXC3H5OOH];

	M[mM9] = c[sN2] + c[sAR] + c[sH] + c[sO2]
		 + c[sO] + c[sOH] + 2 * c[sH2] + 6 * c[sH2O]
		 + c[sHO2] + c[sH2O2] + 1.5 * c[sCO] + 2 * c[sCO2]
		 + c[sHCO] + c[sCH3] + 3 * c[sCH4] + c[sCH2O]
		 + c[sHOCHO] + c[sCH3O] + c[sCH3OH] + c[sCH3O2]
		 + c[sCH3O2H] + 3 * c[sC2H6] + c[sC2H5] + c[sC2H4]
		 + c[sCH3CHO] + c[sCH3CO] + c[sCH3CO3] + c[sCH3CO3H]
		 + c[sCH2CHO] + c[sCH3CO2] + c[sC2H3] + c[sC2H3CHO]
		 + c[sC2H3CO] + c[sC3H8] + c[sNXC3H7] + c[sIXC3H7]
		 + c[sC3H5XAXC3H5] + c[sC3H6] + c[sIXC3H7O2] + c[sIXC3H7O2H]
		 + c[sC3H5XTXC3H5] + c[sC3H5O] + c[sCH3COCH2] + c[sIXC3H7O]
		 + c[sNXC3H7O2] + c[sC3H6OOH1X2XC3H7O2] + c[sC3H6OOH1X3XC3H7O2] + c[sC3H6OOH2X1XC3H7O2]
		 + c[sC3H6OOH1X2O2XC3H7O4] + c[sC3H6OOH1X3O2XC3H7O4] + c[sC3H6OOH2X1O2XC3H7O4] + c[sC3KET12XC3H6O3]
		 + c[sC3KET13XC3H6O3] + c[sC3KET21XC3H6O3] + c[sC3H52X1J3OOHXC3H7O4] + c[sAXC3H5OOH];

	M[mM10] = c[sN2] + c[sAR] + c[sH] + c[sO2]
		 + c[sO] + c[sOH] + 2.5 * c[sH2] + 6.3 * c[sH2O]
		 + c[sHO2] + c[sH2O2] + 1.9 * c[sCO] + 3.8 * c[sCO2]
		 + c[sHCO] + c[sCH3] + 2 * c[sCH4] + c[sCH2O]
		 + c[sHOCHO] + c[sCH3O] + c[sCH3OH] + c[sCH3O2]
		 + c[sCH3O2H] + 3 * c[sC2H6] + c[sC2H5] + c[sC2H4]
		 + c[sCH3CHO] + c[sCH3CO] + c[sCH3CO3] + c[sCH3CO3H]
		 + c[sCH2CHO] + c[sCH3CO2] + c[sC2H3] + c[sC2H3CHO]
		 + c[sC2H3CO] + c[sC3H8] + c[sNXC3H7] + c[sIXC3H7]
		 + c[sC3H5XAXC3H5] + c[sC3H6] + c[sIXC3H7O2] + c[sIXC3H7O2H]
		 + c[sC3H5XTXC3H5] + c[sC3H5O] + c[sCH3COCH2] + c[sIXC3H7O]
		 + c[sNXC3H7O2] + c[sC3H6OOH1X2XC3H7O2] + c[sC3H6OOH1X3XC3H7O2] + c[sC3H6OOH2X1XC3H7O2]
		 + c[sC3H6OOH1X2O2XC3H7O4] + c[sC3H6OOH1X3O2XC3H7O4] + c[sC3H6OOH2X1O2XC3H7O4] + c[sC3KET12XC3H6O3]
		 + c[sC3KET13XC3H6O3] + c[sC3KET21XC3H6O3] + c[sC3H52X1J3OOHXC3H7O4] + c[sAXC3H5OOH];

	M[mM11] = c[sN2] + c[sAR] + c[sH] + c[sO2]
		 + c[sO] + c[sOH] + 2 * c[sH2] + 6 * c[sH2O]
		 + c[sHO2] + c[sH2O2] + 1.5 * c[sCO] + 2 * c[sCO2]
		 + c[sHCO] + c[sCH3] + 2 * c[sCH4] + c[sCH2O]
		 + c[sHOCHO] + c[sCH3O] + c[sCH3OH] + c[sCH3O2]
		 + c[sCH3O2H] + 3 * c[sC2H6] + c[sC2H5] + c[sC2H4]
		 + c[sCH3CHO] + c[sCH3CO] + c[sCH3CO3] + c[sCH3CO3H]
		 + c[sCH2CHO] + c[sCH3CO2] + c[sC2H3] + c[sC2H3CHO]
		 + c[sC2H3CO] + c[sC3H8] + c[sNXC3H7] + c[sIXC3H7]
		 + c[sC3H5XAXC3H5] + c[sC3H6] + c[sIXC3H7O2] + c[sIXC3H7O2H]
		 + c[sC3H5XTXC3H5] + c[sC3H5O] + c[sCH3COCH2] + c[sIXC3H7O]
		 + c[sNXC3H7O2] + c[sC3H6OOH1X2XC3H7O2] + c[sC3H6OOH1X3XC3H7O2] + c[sC3H6OOH2X1XC3H7O2]
		 + c[sC3H6OOH1X2O2XC3H7O4] + c[sC3H6OOH1X3O2XC3H7O4] + c[sC3H6OOH2X1O2XC3H7O4] + c[sC3KET12XC3H6O3]
		 + c[sC3KET13XC3H6O3] + c[sC3KET21XC3H6O3] + c[sC3H52X1J3OOHXC3H7O4] + c[sAXC3H5OOH];

	M[mM12] = c[sN2] + c[sAR] + c[sH] + c[sO2]
		 + c[sO] + c[sOH] + 2 * c[sH2] + 6 * c[sH2O]
		 + c[sHO2] + c[sH2O2] + 1.5 * c[sCO] + 2 * c[sCO2]
		 + c[sHCO] + c[sCH3] + 2 * c[sCH4] + c[sCH2O]
		 + c[sHOCHO] + c[sCH3O] + c[sCH3OH] + c[sCH3O2]
		 + c[sCH3O2H] + 3 * c[sC2H6] + c[sC2H5] + c[sC2H4]
		 + c[sCH3CHO] + c[sCH3CO] + c[sCH3CO3] + c[sCH3CO3H]
		 + c[sCH2CHO] + c[sCH3CO2] + c[sC2H3] + c[sC2H3CHO]
		 + c[sC2H3CO] + c[sC3H8] + c[sNXC3H7] + c[sIXC3H7]
		 + c[sC3H5XAXC3H5] + c[sC3H6] + c[sIXC3H7O2] + c[sIXC3H7O2H]
		 + c[sC3H5XTXC3H5] + c[sC3H5O] + c[sCH3COCH2] + c[sIXC3H7O]
		 + c[sNXC3H7O2] + c[sC3H6OOH1X2XC3H7O2] + c[sC3H6OOH1X3XC3H7O2] + c[sC3H6OOH2X1XC3H7O2]
		 + c[sC3H6OOH1X2O2XC3H7O4] + c[sC3H6OOH1X3O2XC3H7O4] + c[sC3H6OOH2X1O2XC3H7O4] + c[sC3KET12XC3H6O3]
		 + c[sC3KET13XC3H6O3] + c[sC3KET21XC3H6O3] + c[sC3H52X1J3OOHXC3H7O4] + c[sAXC3H5OOH];

	M[mM13] = c[sN2] + c[sAR] + c[sH] + c[sO2]
		 + c[sO] + c[sOH] + 2 * c[sH2] + 6 * c[sH2O]
		 + c[sHO2] + c[sH2O2] + 1.5 * c[sCO] + 2 * c[sCO2]
		 + c[sHCO] + c[sCH3] + 2 * c[sCH4] + c[sCH2O]
		 + c[sHOCHO] + c[sCH3O] + c[sCH3OH] + c[sCH3O2]
		 + c[sCH3O2H] + 3 * c[sC2H6] + c[sC2H5] + c[sC2H4]
		 + c[sCH3CHO] + c[sCH3CO] + c[sCH3CO3] + c[sCH3CO3H]
		 + c[sCH2CHO] + c[sCH3CO2] + c[sC2H3] + c[sC2H3CHO]
		 + c[sC2H3CO] + c[sC3H8] + c[sNXC3H7] + c[sIXC3H7]
		 + c[sC3H5XAXC3H5] + c[sC3H6] + c[sIXC3H7O2] + c[sIXC3H7O2H]
		 + c[sC3H5XTXC3H5] + c[sC3H5O] + c[sCH3COCH2] + c[sIXC3H7O]
		 + c[sNXC3H7O2] + c[sC3H6OOH1X2XC3H7O2] + c[sC3H6OOH1X3XC3H7O2] + c[sC3H6OOH2X1XC3H7O2]
		 + c[sC3H6OOH1X2O2XC3H7O4] + c[sC3H6OOH1X3O2XC3H7O4] + c[sC3H6OOH2X1O2XC3H7O4] + c[sC3KET12XC3H6O3]
		 + c[sC3KET13XC3H6O3] + c[sC3KET21XC3H6O3] + c[sC3H52X1J3OOHXC3H7O4] + c[sAXC3H5OOH];

	M[mM21] = c[sN2] + c[sAR] + c[sH] + c[sO2]
		 + c[sO] + c[sOH] + 2 * c[sH2] + 6 * c[sH2O]
		 + c[sHO2] + c[sH2O2] + 1.5 * c[sCO] + 2 * c[sCO2]
		 + c[sHCO] + c[sCH3] + 2 * c[sCH4] + c[sCH2O]
		 + c[sHOCHO] + c[sCH3O] + c[sCH3OH] + c[sCH3O2]
		 + c[sCH3O2H] + 3 * c[sC2H6] + c[sC2H5] + c[sC2H4]
		 + c[sCH3CHO] + c[sCH3CO] + c[sCH3CO3] + c[sCH3CO3H]
		 + c[sCH2CHO] + c[sCH3CO2] + c[sC2H3] + c[sC2H3CHO]
		 + c[sC2H3CO] + c[sC3H8] + c[sNXC3H7] + c[sIXC3H7]
		 + c[sC3H5XAXC3H5] + c[sC3H6] + c[sIXC3H7O2] + c[sIXC3H7O2H]
		 + c[sC3H5XTXC3H5] + c[sC3H5O] + c[sCH3COCH2] + c[sIXC3H7O]
		 + c[sNXC3H7O2] + c[sC3H6OOH1X2XC3H7O2] + c[sC3H6OOH1X3XC3H7O2] + c[sC3H6OOH2X1XC3H7O2]
		 + c[sC3H6OOH1X2O2XC3H7O4] + c[sC3H6OOH1X3O2XC3H7O4] + c[sC3H6OOH2X1O2XC3H7O4] + c[sC3KET12XC3H6O3]
		 + c[sC3KET13XC3H6O3] + c[sC3KET21XC3H6O3] + c[sC3H52X1J3OOHXC3H7O4] + c[sAXC3H5OOH];


/* 	printf("M[0] = %e, %e, %e, %e, %e, %e\n", M[0], M[1], M[2], M[3], M[4], M[5]); */

	k[r1f] = 1.9150000000E+11 * exp( -68785000 / rt );
	k[r1b] = 5.4810000000E+08 * exp( 0.39 * lgt + 1226000 / rt );
	k[r2f] = 5.0800000000E+01 * exp( 2.67 * lgt - 26326000 / rt );
	k[r2b] = 2.6670000000E+01 * exp( 2.65 * lgt - 20418000 / rt );
	k[r3f] = 2.1600000000E+05 * exp( 1.51 * lgt - 14351000 / rt );
	k[r3b] = 2.2980000000E+06 * exp( 1.4 * lgt - 76651000 / rt );
	k[r4f] = 2.9700000000E+03 * exp( 2.02 * lgt - 56066000 / rt );
	k[r4b] = 1.4650000000E+02 * exp( 2.11 * lgt + 12150000 / rt );
	k[r5f] = 4.5770000000E+16 * exp( -1.4 * lgt - 436810000 / rt );
	k[r5b] = 1.1460000000E+14 * exp( -1.68 * lgt - 3431000 / rt );
	k[r6f] = 4.5150000000E+14 * exp( -0.64 * lgt - 497478000 / rt );
	k[r6b] = 6.1650000000E+09 * exp( -0.5 * lgt );
	k[r7f] = 9.8800000000E+14 * exp( -0.74 * lgt - 427186000 / rt );
	k[r7b] = 4.7140000000E+12 * exp( -1 * lgt );
	k[r8f] = 1.9120000000E+20 * exp( -1.83 * lgt - 495804000 / rt );
	k[r8b] = 4.5000000000E+16 * exp( -2 * lgt );
	 	k9f0 = 3.4820000000E+10 * exp( -0.411 * lgt + 4665000 / rt );
	 	k9fInf = 1.4750000000E+09 * exp( 0.6 * lgt );
	 	fc9 = 0.5 * exp( -temp / 1e-30 ) + 0.5 * exp( -temp / 1e+30 ) + 1 * exp( -0.1 / temp );
	k[r9f] = GetLindRateCoeff( temp, pressure, k9f0, k9fInf
				, fc9, M[mM5] );
	 	k9b0 = 4.9387649140E+13 * exp( -0.423907 * lgt - 199998348.7 / rt );
	 	k9bInf = 2.0920959931E+12 * exp( 0.587093 * lgt - 204663348.7 / rt );
	k[r9b] = GetLindRateCoeff( temp, pressure, k9b0, k9bInf
				, fc9, M[mM5] );
	k[r10] = 1.6600000000E+10 * exp( -3443000 / rt );
	k[r11f] = 7.0790000000E+10 * exp( -1234000 / rt );
	k[r11b] = 2.0270000000E+07 * exp( 0.72 * lgt - 154139000 / rt );
	k[r12] = 3.2500000000E+10;
	k[r13f] = 2.8900000000E+10 * exp( 2079000 / rt );
	k[r13b] = 5.8610000000E+10 * exp( 0.24 * lgt - 289031000 / rt );
	k[r14f] = 4.6340000000E+13 * exp( -0.35 * lgt - 212003000 / rt );
	k[r14b] = 4.2000000000E+11 * exp( -50124000 / rt );
	k[r15f] = 1.4340000000E+10 * exp( -0.35 * lgt - 155059000 / rt );
	k[r15b] = 1.3000000000E+08 * exp( 6816000 / rt );
	 	k16f0 = 1.2020000000E+14 * exp( -190372000 / rt );
	 	k16fInf = 2.9510000000E+14 * exp( -202631000 / rt );
	 	fc16 = 0.5 * exp( -temp / 1e-30 ) + 0.5 * exp( -temp / 1e+30 ) + 1 * exp( -0.1 / temp );
	k[r16f] = GetLindRateCoeff( temp, pressure, k16f0, k16fInf
				, fc16, M[mM6] );
	 	k16b0 = 1.0903484865E+05 * exp( 1.16428 * lgt + 24235773.55 / rt );
	 	k16bInf = 2.6768871743E+05 * exp( 1.16428 * lgt + 11976773.55 / rt );
	k[r16b] = GetLindRateCoeff( temp, pressure, k16b0, k16bInf
				, fc16, M[mM6] );
	k[r17] = 2.4100000000E+10 * exp( -16610000 / rt );
	k[r18f] = 6.0250000000E+10 * exp( -33263000 / rt );
	k[r18b] = 1.0410000000E+08 * exp( 0.7 * lgt - 100207000 / rt );
	k[r19] = 9.5500000000E+03 * exp( 2 * lgt - 16610000 / rt );
	k[r20] = 1.0000000000E+09;
	k[r21f] = 5.8000000000E+11 * exp( -39986000 / rt );
	k[r21b] = 1.0660000000E+10 * exp( 0.59 * lgt - 169243000 / rt );
	k[r23f] = 1.0680000000E-18 * exp( 7.13 * lgt - 55731000 / rt );
	k[r23b] = 3.5750000000E-16 * exp( 6.69 * lgt - 87153000 / rt );
	k[r24f] = 1.4000000000E+02 * exp( 1.95 * lgt + 5636000 / rt );
	k[r24b] = 1.6370000000E+07 * exp( 1.13 * lgt - 95814000 / rt );
	k[r25] = 3.0100000000E+10 * exp( -96232000 / rt );
	k[r26] = 1.8600000000E+14 * exp( -1 * lgt - 71128000 / rt );
	k[r27] = 7.5800000000E+09 * exp( -1715000 / rt );
	k[r28] = 7.3400000000E+10;
	k[r31] = 1.0200000000E+11;
	k[r32] = 1.2100000000E+11;
	k[r33f] = 1.1680000000E+02 * exp( 2.01 * lgt + 19736000 / rt );
	k[r33b] = 6.0600000000E+04 * exp( 1.72 * lgt - 155812000 / rt );
	k[r34] = 3.0000000000E+10;
	k[r41] = 2.8910000000E+24 * exp( -2.93 * lgt - 381832000 / rt );
	k[r43] = 3.4300000000E+06 * exp( 1.18 * lgt + 1870000 / rt );
	k[r44f] = 5.7400000000E+04 * exp( 1.9 * lgt - 11464000 / rt );
	k[r44b] = 2.1090000000E+01 * exp( 2.54 * lgt - 64726000 / rt );
	k[r45] = 6.2600000000E+06 * exp( 1.15 * lgt - 9456000 / rt );
	k[r46] = 3.6360000000E-09 * exp( 5.42 * lgt - 4176000 / rt );
	k[r47f] = 5.8200000000E-06 * exp( 4.53 * lgt - 27434000 / rt );
	k[r47b] = 1.2380000000E-06 * exp( 4.47 * lgt - 13757000 / rt );
	k[r50f] = 2.3000000000E+13 * exp( -209200000 / rt );
	k[r50b] = 3.3460000000E+04 * exp( 1.18 * lgt - 186899000 / rt );
	k[r51f] = 1.5000000000E+16 * exp( -238488000 / rt );
	k[r51b] = 2.3990000000E+11 * exp( 0.46 * lgt - 255308000 / rt );
	k[r55] = 1.8500000000E+04 * exp( 1.51 * lgt + 4025000 / rt );
	k[r60] = 1.0000000000E+09 * exp( -49873000 / rt );
	 	k63f0 = 1.8670000000E+22 * exp( -3 * lgt - 101700000 / rt );
	 	k63fInf = 6.8000000000E+13 * exp( -109495000 / rt );
	 	fc63 = 0.1 * exp( -temp / 2500 ) + 0.9 * exp( -temp / 1300 ) + 1 * exp( -0.1 / temp );
	k[r63f] = GetLindRateCoeff( temp, pressure, k63f0, k63fInf
				, fc63, M[mM0] );
	 	k63b0 = 2.1529707283E+17 * exp( -2.67828 * lgt - 14765645.8 / rt );
	 	k63bInf = 7.8415645164E+08 * exp( 0.321724 * lgt - 22560645.8 / rt );
	k[r63b] = GetLindRateCoeff( temp, pressure, k63b0, k63bInf
				, fc63, M[mM0] );
	k[r64] = 4.3800000000E-22 * exp( 9.5 * lgt + 23016000 / rt );
	k[r65] = 6.6200000000E+08 * exp( -9598000 / rt );
	k[r67] = 1.2000000000E+10;
	k[r68] = 2.0000000000E+10;
	k[r69f] = 3.0100000000E+08;
	k[r69b] = 1.0580000000E+09 * exp( -0.03 * lgt - 286646000 / rt );
	 	k86f0 = 2.9500000000E+41 * exp( -7.35 * lgt - 399405000 / rt );
	 	k86fInf = 1.9000000000E+16 * exp( -383798000 / rt );
	 	fc86 = 0.586 * exp( -temp / 279 ) + 0.414 * exp( -temp / 5459 ) + 1 * exp( -0.1 / temp );
	k[r86f] = GetLindRateCoeff( temp, pressure, k86f0, k86fInf
				, fc86, M[mM0] );
	 	k86b0 = 1.2770094813E+33 * exp( -6.55495 * lgt - 10542309.78 / rt );
	 	k86bInf = 8.2248068290E+07 * exp( 0.795053 * lgt + 5064690.219 / rt );
	k[r86b] = GetLindRateCoeff( temp, pressure, k86b0, k86bInf
				, fc86, M[mM0] );
	k[r88f] = 3.6000000000E+09 * exp( -25501000 / rt );
	k[r88b] = 3.4360000000E+07 * exp( 0.61 * lgt - 24623000 / rt );
	k[r91f] = 5.1300000000E+02 * exp( 2.13 * lgt - 10251000 / rt );
	k[r91b] = 5.2100000000E+01 * exp( 2.64 * lgt - 71672000 / rt );
	k[r97f] = 3.8350000000E+10 * exp( 0.05 * lgt - 354468000 / rt );
	k[r97b] = 6.0300000000E+10;
	 	k98f0 = 3.3100000000E+24 * exp( -4 * lgt - 8820000 / rt );
	 	k98fInf = 2.1400000000E+12 * exp( -0.4 * lgt );
	 	fc98 = 1 * exp( -temp / 1e-15 ) + 1 * exp( -40 / temp );
	k[r98f] = GetLindRateCoeff( temp, pressure, k98f0, k98fInf
				, fc98, M[mM9] );
	 	k98b0 = 2.4823314815E+29 * exp( -3.92933 * lgt - 448865305.6 / rt );
	 	k98bInf = 1.6048910485E+17 * exp( -0.329328 * lgt - 440045305.6 / rt );
	k[r98b] = GetLindRateCoeff( temp, pressure, k98b0, k98bInf
				, fc98, M[mM9] );
	k[r99f] = 3.8600000000E+03 * exp( 2.11 * lgt - 32426000 / rt );
	k[r99b] = 6.7990000000E-01 * exp( 2.75 * lgt - 22301000 / rt );
	k[r100f] = 4.1600000000E+01 * exp( 2.64 * lgt - 8786000 / rt );
	k[r100b] = 7.7960000000E-02 * exp( 3.17 * lgt - 60961000 / rt );
	k[r101f] = 6.9200000000E+05 * exp( 1.56 * lgt - 35501000 / rt );
	k[r101b] = 6.3980000000E+01 * exp( 2.18 * lgt - 19468000 / rt );
	k[r102f] = 9.0400000000E+09 * exp( -103428000 / rt );
	k[r102b] = 9.2170000000E+08 * exp( -0.06 * lgt - 26355000 / rt );
	k[r108f] = 1.1500000000E+10;
	k[r108b] = 7.8270000000E+10 * exp( -0.13 * lgt - 99286000 / rt );
	k[r109f] = 3.0000000000E+09;
	k[r109b] = 3.2470000000E+12 * exp( -0.29 * lgt - 238948000 / rt );
	k[r110] = 8.0000000000E+10;
	k[r111] = 6.0820000000E+04 * exp( 1.54 * lgt - 110416000 / rt );
	k[r112f] = 6.8650000000E-02 * exp( 2.86 * lgt - 42551000 / rt );
	k[r112b] = 1.4880000000E-02 * exp( 3.05 * lgt - 266604000 / rt );
	 	k113f0 = 3.8160000000E+25 * exp( -4.89 * lgt - 14359000 / rt );
	 	k113fInf = 1.0060000000E+05 * exp( 1.63 * lgt );
	 	fc113 = 0.955 * exp( -temp / 880.1 ) + 0.045 * exp( -temp / 2.5e+09 ) + 1 * exp( -1.786e+09 / temp );
	k[r113f] = GetLindRateCoeff( temp, pressure, k113f0, k113fInf
				, fc113, M[mM0] );
	 	k113b0 = 8.2867599918E+33 * exp( -5.93794 * lgt - 154216737.5 / rt );
	 	k113bInf = 2.1846123039E+13 * exp( 0.58206 * lgt - 139857737.5 / rt );
	k[r113b] = GetLindRateCoeff( temp, pressure, k113b0, k113bInf
				, fc113, M[mM0] );
	k[r114f] = 1.9900000000E+09 * exp( -48827000 / rt );
	k[r114b] = 8.5040000000E+09 * exp( -0.5 * lgt - 29326000 / rt );
	k[r115f] = 1.8100000000E+08 * exp( -77320000 / rt );
	k[r115b] = 3.7080000000E+08 * exp( -0.5 * lgt + 5552000 / rt );
	k[r117f] = 7.5000000000E+09 * exp( 4184000 / rt );
	k[r117b] = 3.1830000000E+13 * exp( -0.93 * lgt - 118449000 / rt );
	k[r118] = 1.7500000000E+08 * exp( 7159000 / rt );
	k[r119] = 3.1100000000E+11 * exp( -1.61 * lgt + 4397000 / rt );
	k[r120] = 1.4000000000E+13 * exp( -1.61 * lgt - 7782000 / rt );
	k[r121f] = 9.6000000000E+10;
	k[r121b] = 1.7140000000E+10 * exp( -0.08 * lgt - 176272000 / rt );
	k[r122] = 3.6000000000E+10;
	k[r123] = 6.0000000000E+10;
	k[r124] = 6.3100000000E+14 * exp( -176983000 / rt );
	 	k144f0 = 3.4000000000E+35 * exp( -7.03 * lgt - 11560000 / rt );
	 	k144fInf = 6.7700000000E+13 * exp( -1.18 * lgt - 2736000 / rt );
	 	fc144 = 0.381 * exp( -temp / 73.2 ) + 0.619 * exp( -temp / 1180 ) + 1 * exp( -9999 / temp );
	k[r144f] = GetLindRateCoeff( temp, pressure, k144f0, k144fInf
				, fc144, M[mM10] );
	 	k144b0 = 1.4511995024E+47 * exp( -8.54245 * lgt - 396160526.7 / rt );
	 	k144bInf = 2.8895943033E+25 * exp( -2.69245 * lgt - 387336526.7 / rt );
	k[r144b] = GetLindRateCoeff( temp, pressure, k144b0, k144bInf
				, fc144, M[mM10] );
	 	k145f0 = 1.9900000000E+35 * exp( -7.08 * lgt - 27970000 / rt );
	 	k145fInf = 5.2700000000E+14 * exp( -0.99 * lgt - 6611000 / rt );
	 	fc145 = 0.158 * exp( -temp / 125 ) + 0.842 * exp( -temp / 2219 ) + 1 * exp( -6882 / temp );
	k[r145f] = GetLindRateCoeff( temp, pressure, k145f0, k145fInf
				, fc145, M[mM11] );
	 	k145b0 = 4.3268672779E+42 * exp( -7.68048 * lgt - 450132234.3 / rt );
	 	k145bInf = 1.1458588218E+22 * exp( -1.59048 * lgt - 428773234.3 / rt );
	k[r145b] = GetLindRateCoeff( temp, pressure, k145b0, k145bInf
				, fc145, M[mM11] );
	k[r146f] = 5.5400000000E-01 * exp( 3.5 * lgt - 21619000 / rt );
	k[r146b] = 6.2350000000E-07 * exp( 4.7 * lgt - 28941000 / rt );
	k[r147] = 3.5500000000E+03 * exp( 2.4 * lgt - 24393000 / rt );
	k[r148] = 5.8000000000E+04 * exp( 1.73 * lgt - 4853000 / rt );
	k[r150] = 1.5100000000E-10 * exp( 6 * lgt - 25301000 / rt );
	k[r151] = 1.9040000000E+09 * exp( -66484000 / rt );
	k[r152] = 1.3200000000E+10 * exp( -85646000 / rt );
	 	k156f0 = 4.7150000000E+12 * exp( -3159000 / rt );
	 	k156fInf = 3.9750000000E+06 * exp( 1.28 * lgt - 5406000 / rt );
	 	fc156 = 0.24 * exp( -temp / 40 ) + 0.76 * exp( -temp / 1025 ) + 1 * exp( -0.1 / temp );
	k[r156f] = GetLindRateCoeff( temp, pressure, k156f0, k156fInf
				, fc156, M[mM12] );
	 	k156b0 = 1.2019475679E+15 * exp( 0.122938 * lgt - 156070232.3 / rt );
	 	k156bInf = 1.0133068044E+09 * exp( 1.40294 * lgt - 158317232.3 / rt );
	k[r156b] = GetLindRateCoeff( temp, pressure, k156b0, k156bInf
				, fc156, M[mM12] );
	k[r157f] = 3.0000000000E+17 * exp( -2.86 * lgt - 28288000 / rt );
	k[r157b] = 4.1150000000E+17 * exp( -2.77 * lgt - 82927000 / rt );
	k[r158f] = 1.6000000000E+11 * exp( -1.17 * lgt - 43472000 / rt );
	k[r158b] = 2.8950000000E+10 * exp( -0.87 * lgt - 287441000 / rt );
	k[r161f] = 1.9500000000E+10 * exp( -0.5 * lgt );
	k[r161b] = 2.8950000000E+13 * exp( -0.7 * lgt - 293591000 / rt );
	k[r162f] = 9.8540000000E+08 * exp( 0.49 * lgt + 18577000 / rt );
	k[r162b] = 6.8400000000E+09 * exp( 0.1 * lgt - 44350000 / rt );
	k[r163] = 2.0000000000E+09;
	k[r164] = 1.1000000000E+11;
	k[r193f] = 2.2870000000E+15 * exp( 0.15 * lgt - 337021000 / rt );
	k[r193b] = 1.7500000000E+10;
	k[r194f] = 1.1100000000E+10 * exp( -13012000 / rt );
	k[r194b] = 7.9880000000E+06 * exp( 0.63 * lgt - 73471000 / rt );
	k[r196] = 2.0000000000E+03 * exp( 1.8 * lgt - 5439000 / rt );
	k[r197] = 3.0100000000E+10 * exp( -163804000 / rt );
	k[r198f] = 1.7600000000E+00 * exp( 2.79 * lgt - 20711000 / rt );
	k[r198b] = 7.1910000000E+00 * exp( 2.78 * lgt - 91295000 / rt );
	k[r199] = 3.0100000000E+09 * exp( -49915000 / rt );
	k[r200f] = 3.0100000000E+09 * exp( -49915000 / rt );
	k[r200b] = 2.5190000000E+10 * exp( -0.51 * lgt - 37618000 / rt );
	k[r201] = 3.0100000000E+09 * exp( -49915000 / rt );
	k[r202f] = 3.0000000000E+12 * exp( -1.08 * lgt );
	k[r202b] = 5.0650000000E+15 * exp( -1.87 * lgt - 116022000 / rt );
	k[r203f] = 1.7200000000E+02 * exp( 2.4 * lgt - 3410000 / rt );
	k[r203b] = 1.2160000000E+02 * exp( 2.52 * lgt - 104349000 / rt );
	 	k204f0 = 1.2000000000E+12 * exp( -52375000 / rt );
	 	k204fInf = 3.0000000000E+12 * exp( -69956000 / rt );
	 	fc204 = 1 * exp( -0 / temp );
	k[r204f] = GetLindRateCoeff( temp, pressure, k204f0, k204fInf
				, fc204, M[mM0] );
	 	k204b0 = 7.2799877061E+02 * exp( 1.37732 * lgt + 10376912.91 / rt );
	 	k204bInf = 1.8199969265E+03 * exp( 1.37732 * lgt - 7204087.087 / rt );
	k[r204b] = GetLindRateCoeff( temp, pressure, k204b0, k204bInf
				, fc204, M[mM0] );
	k[r208f] = 4.6250000000E+19 * exp( -1.9 * lgt - 165519000 / rt );
	k[r208b] = 1.2000000000E+08 * exp( 4602000 / rt );
	k[r209] = 1.7500000000E+07 * exp( 13703000 / rt );
	k[r210f] = 2.4100000000E+09 * exp( -41572000 / rt );
	k[r210b] = 3.8450000000E+09 * exp( 0.05 * lgt - 53179000 / rt );
	k[r212f] = 1.9900000000E+09 * exp( -48827000 / rt );
	k[r212b] = 6.7520000000E+08 * exp( -46735000 / rt );
	k[r213] = 1.7000000000E+10 * exp( -85605000 / rt );
	k[r214] = 5.0100000000E+14 * exp( -167988000 / rt );
	k[r215f] = 4.4000000000E+12 * exp( -43932000 / rt );
	k[r215b] = 9.6480000000E+06 * exp( 0.2 * lgt - 88031000 / rt );
	k[r217] = 2.0000000000E+10 * exp( -17573000 / rt );
	 	k233f0 = 9.8000000000E+23 * exp( -3.86 * lgt - 13891000 / rt );
	 	k233fInf = 6.1000000000E+09 * exp( 0.27 * lgt - 1172000 / rt );
	 	fc233 = 0.218 * exp( -temp / 208 ) + 0.782 * exp( -temp / 2663 ) + 1 * exp( -6095 / temp );
	k[r233f] = GetLindRateCoeff( temp, pressure, k233f0, k233fInf
				, fc233, M[mM13] );
	 	k233b0 = 1.4917429652E+30 * exp( -4.21453 * lgt - 466842993.9 / rt );
	 	k233bInf = 9.2853388652E+15 * exp( -0.0845323 * lgt - 454123993.9 / rt );
	k[r233b] = GetLindRateCoeff( temp, pressure, k233b0, k233bInf
				, fc233, M[mM13] );
	k[r235f] = 8.4200000000E-06 * exp( 4.62 * lgt - 10807000 / rt );
	k[r235b] = 2.6330000000E-06 * exp( 4.43 * lgt - 5410000 / rt );
	k[r236] = 1.0200000000E+04 * exp( 1.88 * lgt - 749000 / rt );
	k[r237] = 3.3900000000E+03 * exp( 1.88 * lgt - 749000 / rt );
	k[r238] = 2.0200000000E+10 * exp( -24836000 / rt );
	k[r239f] = 6.6200000000E-03 * exp( 3.7 * lgt - 39748000 / rt );
	k[r239b] = 1.4400000000E-03 * exp( 4.02 * lgt - 22895000 / rt );
	k[r240f] = 4.0000000000E+10 * exp( -243509000 / rt );
	k[r240b] = 6.5630000000E+10 * exp( -0.54 * lgt - 9293000 / rt );
	k[r241] = 1.2000000000E+08 * exp( -28242000 / rt );
	k[r242] = 2.2300000000E+09 * exp( -71923000 / rt );
	k[r244f] = 1.1300000000E+10 * exp( -127319000 / rt );
	k[r244b] = 3.2640000000E+12 * exp( -0.83 * lgt - 66567000 / rt );
	k[r252] = 4.5700000000E+13 * exp( -1.39 * lgt - 4238000 / rt );
	k[r253f] = 3.0300000000E+08 * exp( -0.29 * lgt - 42000 / rt );
	k[r253b] = 8.4830000000E+07 * exp( 0.01 * lgt - 10456000 / rt );
	k[r314f] = 4.4570000000E+17 * exp( -0.43 * lgt - 394886000 / rt );
	k[r314b] = 1.8100000000E+10;
	k[r315f] = 1.3400000000E+10 * exp( -13807000 / rt );
	k[r315b] = 3.3090000000E+07 * exp( 0.61 * lgt - 94893000 / rt );
	k[r316] = 5.9400000000E+09 * exp( -7816000 / rt );
	k[r317] = 9.2400000000E+03 * exp( 1.5 * lgt + 4025000 / rt );
	k[r318] = 1.0050000000E+10 * exp( -170289000 / rt );
	k[r319f] = 3.0100000000E+09 * exp( -49915000 / rt );
	k[r319b] = 4.3030000000E+09 * exp( -0.08 * lgt - 64015000 / rt );
	k[r320] = 2.6080000000E+03 * exp( 1.78 * lgt - 24732000 / rt );
	k[r322] = 1.0000000000E+09 * exp( -13807000 / rt );
	k[r323f] = 3.0100000000E+09 * exp( -49915000 / rt );
	k[r323b] = 8.6450000000E+10 * exp( -0.53 * lgt - 58241000 / rt );
	k[r324] = 3.0430000000E+14 * exp( -0.46 * lgt - 127654000 / rt );
	 	k400f0 = 2.7100000000E+68 * exp( -16.82 * lgt - 54664000 / rt );
	 	k400fInf = 1.8800000000E+10;
	 	fc400 = 0.8473 * exp( -temp / 291 ) + 0.1527 * exp( -temp / 2742 ) + 1 * exp( -7748 / temp );
	k[r400f] = GetLindRateCoeff( temp, pressure, k400f0, k400fInf
				, fc400, M[mM21] );
	 	k400b0 = 3.4568521433E+82 * exp( -18.9698 * lgt - 430310132.3 / rt );
	 	k400bInf = 2.3981114500E+24 * exp( -2.14976 * lgt - 375646132.3 / rt );
	k[r400b] = GetLindRateCoeff( temp, pressure, k400b0, k400bInf
				, fc400, M[mM21] );
	k[r401f] = 3.7190000000E+17 * exp( -0.36 * lgt - 423421000 / rt );
	k[r401b] = 1.0000000000E+11;
	k[r402f] = 2.3460000000E+18 * exp( -0.67 * lgt - 412835000 / rt );
	k[r402b] = 1.0000000000E+11;
	k[r403f] = 2.0000000000E+10 * exp( -207694000 / rt );
	k[r403b] = 1.7860000000E+06 * exp( 0.6 * lgt + 665000 / rt );
	k[r404f] = 6.0000000000E+10 * exp( -218781000 / rt );
	k[r404b] = 3.3790000000E+07 * exp( 0.28 * lgt + 247000 / rt );
	k[r405f] = 1.3000000000E+03 * exp( 2.4 * lgt - 18707000 / rt );
	k[r405b] = 2.2120000000E-02 * exp( 3.35 * lgt - 39166000 / rt );
	k[r406f] = 1.3300000000E+03 * exp( 2.54 * lgt - 28284000 / rt );
	k[r406b] = 1.4280000000E-01 * exp( 3.17 * lgt - 38074000 / rt );
	k[r407f] = 5.4900000000E+02 * exp( 2.5 * lgt - 13138000 / rt );
	k[r407b] = 4.9040000000E-03 * exp( 3.43 * lgt - 27690000 / rt );
	k[r408f] = 3.7100000000E+03 * exp( 2.4 * lgt - 23033000 / rt );
	k[r408b] = 2.0910000000E-01 * exp( 3.01 * lgt - 26916000 / rt );
	k[r409] = 1.0540000000E+07 * exp( 0.97 * lgt - 6636000 / rt );
	k[r410] = 4.6700000000E+04 * exp( 1.61 * lgt + 146000 / rt );
	k[r411f] = 5.6000000000E+09 * exp( -74057000 / rt );
	k[r411b] = 4.1600000000E+08 * exp( -31070000 / rt );
	k[r412f] = 1.6800000000E+10 * exp( -85479000 / rt );
	k[r412b] = 2.3300000000E+09 * exp( -41112000 / rt );
	k[r413f] = 6.4000000000E+01 * exp( 2.17 * lgt - 31464000 / rt );
	k[r413b] = 6.1840000000E+00 * exp( 2.48 * lgt - 62049000 / rt );
	k[r414] = 9.0400000000E-04 * exp( 3.65 * lgt - 29932000 / rt );
	k[r415f] = 3.0000000000E+07 * exp( -53974000 / rt );
	k[r415b] = 3.0000000000E+07 * exp( -53974000 / rt );
	k[r416] = 1.0000000000E+08 * exp( -43514000 / rt );
	k[r419] = 1.0000000000E+08 * exp( -43514000 / rt );
	k[r420f] = 7.9400000000E+08 * exp( -85772000 / rt );
	k[r420b] = 5.3720000000E+13 * exp( -1.33 * lgt - 56066000 / rt );
	k[r421f] = 7.9400000000E+08 * exp( -67781000 / rt );
	k[r421b] = 5.3720000000E+13 * exp( -1.33 * lgt - 38053000 / rt );
	k[r422f] = 3.0000000000E+08 * exp( -29288000 / rt );
	k[r422b] = 1.2200000000E+07 * exp( -38417000 / rt );
	k[r423f] = 3.0000000000E+08 * exp( -29288000 / rt );
	k[r423b] = 1.2200000000E+07 * exp( -38417000 / rt );
	k[r424] = 1.7000000000E+10 * exp( -85605000 / rt );
	k[r425] = 2.0000000000E+09 * exp( -71128000 / rt );
	k[r430] = 1.7000000000E+10 * exp( -85605000 / rt );
	k[r431] = 2.0000000000E+09 * exp( -71128000 / rt );
	k[r432f] = 2.0000000000E+09 * exp( -71128000 / rt );
	k[r432b] = 3.1430000000E+07 * exp( 0.3 * lgt - 36233000 / rt );
	k[r433f] = 1.7000000000E+10 * exp( -85605000 / rt );
	k[r433b] = 1.6850000000E+09 * exp( -0.01 * lgt - 40041000 / rt );
	k[r436f] = 1.7300000000E+50 * exp( -10.98 * lgt - 227066000 / rt );
	k[r436b] = 2.0900000000E+43 * exp( -9.9 * lgt - 69831000 / rt );
	k[r437] = 2.0000000000E+10;
	k[r439] = 2.4100000000E+10;
	k[r441f] = 4.8180000000E+10;
	k[r441b] = 1.2790000000E+08 * exp( 0.8 * lgt - 361832000 / rt );
	k[r442f] = 9.9700000000E+40 * exp( -8.6 * lgt - 173343000 / rt );
	k[r442b] = 1.6300000000E+36 * exp( -8.36 * lgt - 84433000 / rt );
	k[r443f] = 8.7800000000E+39 * exp( -8.1 * lgt - 194891000 / rt );
	k[r443b] = 2.0700000000E+34 * exp( -7.39 * lgt - 50334000 / rt );
	k[r448f] = 2.7300000000E+62 * exp( -13.28 * lgt - 515469000 / rt );
	k[r448b] = 4.7120000000E+56 * exp( -13.19 * lgt - 123595000 / rt );
	k[r449f] = 2.0100000000E+61 * exp( -13.26 * lgt - 495804000 / rt );
	k[r449b] = 1.0620000000E+56 * exp( -12.89 * lgt - 125604000 / rt );
	k[r452] = 1.5800000000E+04 * exp( 1.76 * lgt + 5088000 / rt );
	k[r455] = 5.2400000000E+08 * exp( 0.7 * lgt - 24619000 / rt );
	k[r457f] = 6.0300000000E+07 * exp( 0.7 * lgt - 31932000 / rt );
	k[r457b] = 9.9290000000E+01 * exp( 2 * lgt + 59000 / rt );
	k[r458f] = 3.1200000000E+03 * exp( 2 * lgt + 1247000 / rt );
	k[r458b] = 7.0040000000E+01 * exp( 2.54 * lgt - 124181000 / rt );
	k[r460f] = 1.1100000000E+03 * exp( 2 * lgt - 6071000 / rt );
	k[r460b] = 3.7050000000E-02 * exp( 3.22 * lgt - 42300000 / rt );
	k[r461] = 1.5000000000E+08 * exp( -59371000 / rt );
	k[r463f] = 3.0000000000E+06 * exp( -41547000 / rt );
	k[r463b] = 9.2340000000E+00 * exp( 1.27 * lgt + 58827000 / rt );
	k[r464f] = 1.7300000000E+02 * exp( 2.5 * lgt - 10427000 / rt );
	k[r464b] = 3.6500000000E-01 * exp( 3.15 * lgt - 73555000 / rt );
	k[r466f] = 4.0500000000E+02 * exp( 2.5 * lgt - 40978000 / rt );
	k[r466b] = 1.2700000000E-03 * exp( 3.82 * lgt - 14895000 / rt );
	k[r467f] = 4.8300000000E+31 * exp( -5.81 * lgt - 77404000 / rt );
	k[r467b] = 1.0640000000E+29 * exp( -5.26 * lgt - 124181000 / rt );
	k[r468f] = 4.0000000000E+09 * exp( -166942000 / rt );
	k[r468b] = 4.4280000000E+07 * exp( 0.3 * lgt - 1243000 / rt );
	k[r470] = 1.4000000000E+09 * exp( -253969000 / rt );
	k[r471] = 2.2100000000E-03 * exp( 3.5 * lgt - 23744000 / rt );
	k[r473] = 8.4000000000E-04 * exp( 3.5 * lgt - 48785000 / rt );
	k[r475] = 3.2400000000E+08 * exp( -62342000 / rt );
	k[r476] = 3.2400000000E+08 * exp( -62342000 / rt );
	k[r480] = 3.2400000000E+08 * exp( -62342000 / rt );
	k[r486] = 7.0000000000E+09 * exp( 4184000 / rt );
	k[r487] = 7.0000000000E+09 * exp( 4184000 / rt );
	k[r491f] = 4.0000000000E+08;
	k[r491b] = 6.9370000000E+13 * exp( -1.33 * lgt - 220915000 / rt );
	k[r495f] = 7.1400000000E+12 * exp( -1.21 * lgt - 88073000 / rt );
	k[r495b] = 4.9440000000E+13 * exp( -1.4 * lgt - 370786000 / rt );
	k[r496f] = 2.4700000000E+10 * exp( -0.44 * lgt - 96316000 / rt );
	k[r496b] = 3.9620000000E+10 * exp( -0.61 * lgt - 313465000 / rt );
	k[r507f] = 3.8100000000E+14 * exp( -1.36 * lgt - 23347000 / rt );
	k[r507b] = 2.0000000000E+08 * exp( -73220000 / rt );
	k[r508] = 3.7100000000E+22 * exp( -3.96 * lgt - 29468000 / rt );
	k[r558] = 7.0000000000E+09 * exp( 4184000 / rt );
	k[r560] = 7.0000000000E+09 * exp( 4184000 / rt );
	k[r561f] = 3.0630000000E+20 * exp( -1.64 * lgt - 149411000 / rt );
	k[r561b] = 4.5200000000E+09;
	k[r562f] = 3.9780000000E+22 * exp( -2.19 * lgt - 158574000 / rt );
	k[r562b] = 7.5400000000E+09;
	k[r565] = 5.6000000000E+09 * exp( -56902000 / rt );
	k[r566] = 2.8000000000E+09 * exp( -56902000 / rt );
	k[r568] = 1.7500000000E+07 * exp( 13703000 / rt );
	k[r570] = 1.1300000000E+10 * exp( -127319000 / rt );
	k[r574f] = 2.8000000000E+09 * exp( -56902000 / rt );
	k[r574b] = 1.0000000000E+09 * exp( -41840000 / rt );
	k[r576] = 1.1200000000E+10 * exp( -103094000 / rt );
	k[r578] = 1.4000000000E+13 * exp( -1.61 * lgt - 7782000 / rt );
	k[r580f] = 3.0100000000E+10 * exp( -108910000 / rt );
	k[r580b] = 4.8000000000E+10 * exp( -33263000 / rt );
	k[r581f] = 1.7000000000E+10 * exp( -85605000 / rt );
	k[r581b] = 5.0000000000E+08 * exp( -27196000 / rt );
	k[r585] = 1.4000000000E+13 * exp( -1.61 * lgt - 7782000 / rt );
	k[r589] = 1.4000000000E+13 * exp( -1.61 * lgt - 7782000 / rt );
	k[r592] = 7.0000000000E+09 * exp( 4184000 / rt );
	k[r594] = 7.0000000000E+09 * exp( 4184000 / rt );
	k[r596] = 7.0000000000E+09 * exp( 4184000 / rt );
	k[r603] = 9.4500000000E+15 * exp( -178238000 / rt );
	k[r606f] = 2.1670000000E+17 * exp( -1.25 * lgt - 69371000 / rt );
	k[r606b] = 8.5000000000E+07 * exp( -38493000 / rt );
	k[r609f] = 2.0000000000E+12 * exp( -112340000 / rt );
	k[r609b] = 3.7220000000E+08 * exp( 0.58 * lgt - 49036000 / rt );
	k[r610f] = 3.7500000000E+11 * exp( -102090000 / rt );
	k[r610b] = 9.0520000000E+11 * exp( -0.51 * lgt - 37388000 / rt );
	k[r611f] = 6.0000000000E+12 * exp( -123010000 / rt );
	k[r611b] = 3.7410000000E+10 * exp( 0.12 * lgt - 49413000 / rt );
	k[r616f] = 1.0150000000E+14 * exp( -0.67 * lgt - 53304000 / rt );
	k[r616b] = 1.2500000000E+08 * exp( -35982000 / rt );
	k[r617f] = 4.1950000000E+16 * exp( -1.37 * lgt - 62802000 / rt );
	k[r617b] = 1.2500000000E+08 * exp( -35982000 / rt );
	k[r618] = 3.0350000000E+15 * exp( -0.79 * lgt - 114642000 / rt );
	k[r620] = 1.3100000000E+33 * exp( -7.01 * lgt - 201334000 / rt );
	k[r622f] = 2.3910000000E+25 * exp( -2.94 * lgt - 167778000 / rt );
	k[r622b] = 5.0000000000E+09;
	k[r623f] = 2.8530000000E+20 * exp( -1.63 * lgt - 149327000 / rt );
	k[r623b] = 4.5200000000E+09;
	k[r624f] = 5.2270000000E+22 * exp( -2.24 * lgt - 158239000 / rt );
	k[r624b] = 4.5200000000E+09;
	k[r625] = 2.0000000000E+12 * exp( -110458000 / rt );
	k[r626] = 2.5000000000E+11 * exp( -89538000 / rt );
	k[r627] = 1.0000000000E+12 * exp( -99788000 / rt );
	k[r631f] = 3.0680000000E+09 * exp( 0.26 * lgt - 70500000 / rt );
	k[r631b] = 2.0000000000E+12 * exp( -112340000 / rt );
	k[r632f] = 1.1500000000E+14 * exp( -0.63 * lgt - 72174000 / rt );
	k[r632b] = 1.5640000000E+08 * exp( 0.12 * lgt - 42677000 / rt );
	k[r633] = 9.4500000000E+15 * exp( -179912000 / rt );
	k[r634] = 1.0000000000E+16 * exp( -179912000 / rt );
	k[r635] = 1.0000000000E+16 * exp( -179912000 / rt );
	k[r636f] = 3.8800000000E+19 * exp( -1.46 * lgt - 189828000 / rt );
	k[r636b] = 2.0000000000E+10;
	k[r637] = 1.0000000000E+14 * exp( -121754000 / rt );
	k[r638f] = 2.0280000000E+12 * exp( 0.09 * lgt - 98575000 / rt );
	k[r638b] = 1.5000000000E+08 * exp( -44350000 / rt );
	k[r639f] = 1.0000000000E+09 * exp( -25104000 / rt );
	k[r639b] = 1.2880000000E+08 * exp( -133888000 / rt );


	w[r1f] = k[r1f] * c[sH] * c[sO2];
	w[r1b] = k[r1b] * c[sO] * c[sOH];
	w[r2f] = k[r2f] * c[sO] * c[sH2];
	w[r2b] = k[r2b] * c[sH] * c[sOH];
	w[r3f] = k[r3f] * c[sOH] * c[sH2];
	w[r3b] = k[r3b] * c[sH] * c[sH2O];
	w[r4f] = k[r4f] * c[sO] * c[sH2O];
	w[r4b] = k[r4b] * c[sOH] * c[sOH];
	w[r5f] = k[r5f] * c[sH2] * M[mM1];
	w[r5b] = k[r5b] * c[sH] * c[sH] * M[mM1];
	w[r6f] = k[r6f] * c[sO2] * M[mM2];
	w[r6b] = k[r6b] * c[sO] * c[sO] * M[mM2];
	w[r7f] = k[r7f] * c[sOH] * M[mM3];
	w[r7b] = k[r7b] * c[sO] * c[sH] * M[mM3];
	w[r8f] = k[r8f] * c[sH2O] * M[mM4];
	w[r8b] = k[r8b] * c[sH] * c[sOH] * M[mM4];
	w[r9f] = k[r9f] * c[sH] * c[sO2];
	w[r9b] = k[r9b] * c[sHO2];
	w[r10] = k[r10] * c[sHO2] * c[sH];
	w[r11f] = k[r11f] * c[sHO2] * c[sH];
	w[r11b] = k[r11b] * c[sOH] * c[sOH];
	w[r12] = k[r12] * c[sHO2] * c[sO];
	w[r13f] = k[r13f] * c[sHO2] * c[sOH];
	w[r13b] = k[r13b] * c[sH2O] * c[sO2];
	w[r14f] = k[r14f] * c[sH2O2] * c[sO2];
	w[r14b] = k[r14b] * c[sHO2] * c[sHO2];
	w[r15f] = k[r15f] * c[sH2O2] * c[sO2];
	w[r15b] = k[r15b] * c[sHO2] * c[sHO2];
	w[r16f] = k[r16f] * c[sH2O2];
	w[r16b] = k[r16b] * c[sOH] * c[sOH];
	w[r17] = k[r17] * c[sH2O2] * c[sH];
	w[r18f] = k[r18f] * c[sH2O2] * c[sH];
	w[r18b] = k[r18b] * c[sH2] * c[sHO2];
	w[r19] = k[r19] * c[sH2O2] * c[sO];
	w[r20] = k[r20] * c[sH2O2] * c[sOH];
	w[r21f] = k[r21f] * c[sH2O2] * c[sOH];
	w[r21b] = k[r21b] * c[sH2O] * c[sHO2];
	w[r23f] = k[r23f] * c[sCO] * c[sO2];
	w[r23b] = k[r23b] * c[sCO2] * c[sO];
	w[r24f] = k[r24f] * c[sCO] * c[sOH];
	w[r24b] = k[r24b] * c[sCO2] * c[sH];
	w[r25] = k[r25] * c[sCO] * c[sHO2];
	w[r26] = k[r26] * c[sHCO] * M[mM8];
	w[r27] = k[r27] * c[sHCO] * c[sO2];
	w[r28] = k[r28] * c[sHCO] * c[sH];
	w[r31] = k[r31] * c[sHCO] * c[sOH];
	w[r32] = k[r32] * c[sHCO] * c[sCH3];
	w[r33f] = k[r33f] * c[sHCO] * c[sHO2];
	w[r33b] = k[r33b] * c[sCH2O] * c[sO2];
	w[r34] = k[r34] * c[sHCO] * c[sHO2];
	w[r41] = k[r41] * c[sCH2O] * M[mM0];
	w[r43] = k[r43] * c[sCH2O] * c[sOH];
	w[r44f] = k[r44f] * c[sCH2O] * c[sH];
	w[r44b] = k[r44b] * c[sHCO] * c[sH2];
	w[r45] = k[r45] * c[sCH2O] * c[sO];
	w[r46] = k[r46] * c[sCH2O] * c[sCH3];
	w[r47f] = k[r47f] * c[sCH2O] * c[sHO2];
	w[r47b] = k[r47b] * c[sHCO] * c[sH2O2];
	w[r50f] = k[r50f] * c[sHOCHO];
	w[r50b] = k[r50b] * c[sCO] * c[sH2O];
	w[r51f] = k[r51f] * c[sHOCHO];
	w[r51b] = k[r51b] * c[sCO2] * c[sH2];
	w[r55] = k[r55] * c[sHOCHO] * c[sOH];
	w[r60] = k[r60] * c[sHOCHO] * c[sHO2];
	w[r63f] = k[r63f] * c[sCH3O];
	w[r63b] = k[r63b] * c[sH] * c[sCH2O];
	w[r64] = k[r64] * c[sCH3O] * c[sO2];
	w[r65] = k[r65] * c[sCH2O] * c[sCH3O];
	w[r67] = k[r67] * c[sCH3O] * c[sCH3];
	w[r68] = k[r68] * c[sCH3O] * c[sH];
	w[r69f] = k[r69f] * c[sCH3O] * c[sHO2];
	w[r69b] = k[r69b] * c[sCH2O] * c[sH2O2];
	w[r86f] = k[r86f] * c[sCH3OH];
	w[r86b] = k[r86b] * c[sOH] * c[sCH3];
	w[r88f] = k[r88f] * c[sCH3OH] * c[sH];
	w[r88b] = k[r88b] * c[sCH3O] * c[sH2];
	w[r91f] = k[r91f] * c[sCH3OH] * c[sOH];
	w[r91b] = k[r91b] * c[sCH3O] * c[sH2O];
	w[r97f] = k[r97f] * c[sCH3OH] * c[sCH2O];
	w[r97b] = k[r97b] * c[sCH3O] * c[sCH3O];
	w[r98f] = k[r98f] * c[sCH3] * c[sH];
	w[r98b] = k[r98b] * c[sCH4];
	w[r99f] = k[r99f] * c[sCH4] * c[sH];
	w[r99b] = k[r99b] * c[sCH3] * c[sH2];
	w[r100f] = k[r100f] * c[sCH4] * c[sOH];
	w[r100b] = k[r100b] * c[sCH3] * c[sH2O];
	w[r101f] = k[r101f] * c[sCH4] * c[sO];
	w[r101b] = k[r101b] * c[sCH3] * c[sOH];
	w[r102f] = k[r102f] * c[sCH4] * c[sHO2];
	w[r102b] = k[r102b] * c[sCH3] * c[sH2O2];
	w[r108f] = k[r108f] * c[sCH3] * c[sHO2];
	w[r108b] = k[r108b] * c[sCH3O] * c[sOH];
	w[r109f] = k[r109f] * c[sCH3] * c[sHO2];
	w[r109b] = k[r109b] * c[sCH4] * c[sO2];
	w[r110] = k[r110] * c[sCH3] * c[sO];
	w[r111] = k[r111] * c[sCH3] * c[sO2];
	w[r112f] = k[r112f] * c[sCH3] * c[sO2];
	w[r112b] = k[r112b] * c[sCH2O] * c[sOH];
	w[r113f] = k[r113f] * c[sCH3] * c[sO2];
	w[r113b] = k[r113b] * c[sCH3O2];
	w[r114f] = k[r114f] * c[sCH3O2] * c[sCH2O];
	w[r114b] = k[r114b] * c[sCH3O2H] * c[sHCO];
	w[r115f] = k[r115f] * c[sCH4] * c[sCH3O2];
	w[r115b] = k[r115b] * c[sCH3] * c[sCH3O2H];
	w[r117f] = k[r117f] * c[sCH3O2] * c[sCH3];
	w[r117b] = k[r117b] * c[sCH3O] * c[sCH3O];
	w[r118] = k[r118] * c[sCH3O2] * c[sHO2];
	w[r119] = k[r119] * c[sCH3O2] * c[sCH3O2];
	w[r120] = k[r120] * c[sCH3O2] * c[sCH3O2];
	w[r121f] = k[r121f] * c[sCH3O2] * c[sH];
	w[r121b] = k[r121b] * c[sCH3O] * c[sOH];
	w[r122] = k[r122] * c[sCH3O2] * c[sO];
	w[r123] = k[r123] * c[sCH3O2] * c[sOH];
	w[r124] = k[r124] * c[sCH3O2H];
	w[r144f] = k[r144f] * c[sCH3] * c[sCH3];
	w[r144b] = k[r144b] * c[sC2H6];
	w[r145f] = k[r145f] * c[sC2H5] * c[sH];
	w[r145b] = k[r145b] * c[sC2H6];
	w[r146f] = k[r146f] * c[sC2H6] * c[sH];
	w[r146b] = k[r146b] * c[sC2H5] * c[sH2];
	w[r147] = k[r147] * c[sC2H6] * c[sO];
	w[r148] = k[r148] * c[sC2H6] * c[sOH];
	w[r150] = k[r150] * c[sC2H6] * c[sCH3];
	w[r151] = k[r151] * c[sC2H6] * c[sHO2];
	w[r152] = k[r152] * c[sC2H6] * c[sCH3O2];
	w[r156f] = k[r156f] * c[sC2H4] * c[sH];
	w[r156b] = k[r156b] * c[sC2H5];
	w[r157f] = k[r157f] * c[sC2H5] * c[sO2];
	w[r157b] = k[r157b] * c[sC2H4] * c[sHO2];
	w[r158f] = k[r158f] * c[sC2H5] * c[sO2];
	w[r158b] = k[r158b] * c[sCH3CHO] * c[sOH];
	w[r161f] = k[r161f] * c[sCH3] * c[sC2H5];
	w[r161b] = k[r161b] * c[sCH4] * c[sC2H4];
	w[r162f] = k[r162f] * c[sC2H5] * c[sH];
	w[r162b] = k[r162b] * c[sCH3] * c[sCH3];
	w[r163] = k[r163] * c[sC2H5] * c[sH];
	w[r164] = k[r164] * c[sC2H5] * c[sO];
	w[r193f] = k[r193f] * c[sCH3CHO];
	w[r193b] = k[r193b] * c[sCH3] * c[sHCO];
	w[r194f] = k[r194f] * c[sCH3CHO] * c[sH];
	w[r194b] = k[r194b] * c[sCH3CO] * c[sH2];
	w[r196] = k[r196] * c[sCH3CHO] * c[sOH];
	w[r197] = k[r197] * c[sCH3CHO] * c[sO2];
	w[r198f] = k[r198f] * c[sCH3CHO] * c[sCH3];
	w[r198b] = k[r198b] * c[sCH3CO] * c[sCH4];
	w[r199] = k[r199] * c[sCH3CHO] * c[sHO2];
	w[r200f] = k[r200f] * c[sCH3O2] * c[sCH3CHO];
	w[r200b] = k[r200b] * c[sCH3O2H] * c[sCH3CO];
	w[r201] = k[r201] * c[sCH3CHO] * c[sCH3CO3];
	w[r202f] = k[r202f] * c[sCH3CHO] * c[sOH];
	w[r202b] = k[r202b] * c[sCH3] * c[sHOCHO];
	w[r203f] = k[r203f] * c[sCH3CHO] * c[sOH];
	w[r203b] = k[r203b] * c[sCH2CHO] * c[sH2O];
	w[r204f] = k[r204f] * c[sCH3CO];
	w[r204b] = k[r204b] * c[sCO] * c[sCH3];
	w[r208f] = k[r208f] * c[sCH3CO3];
	w[r208b] = k[r208b] * c[sCH3CO] * c[sO2];
	w[r209] = k[r209] * c[sCH3CO3] * c[sHO2];
	w[r210f] = k[r210f] * c[sH2O2] * c[sCH3CO3];
	w[r210b] = k[r210b] * c[sHO2] * c[sCH3CO3H];
	w[r212f] = k[r212f] * c[sCH2O] * c[sCH3CO3];
	w[r212b] = k[r212b] * c[sHCO] * c[sCH3CO3H];
	w[r213] = k[r213] * c[sC2H6] * c[sCH3CO3];
	w[r214] = k[r214] * c[sCH3CO3H];
	w[r215f] = k[r215f] * c[sCH3CO2] * M[mM0];
	w[r215b] = k[r215b] * c[sCH3] * c[sCO2] * M[mM0];
	w[r217] = k[r217] * c[sCH2CHO] * c[sO2];
	w[r233f] = k[r233f] * c[sC2H3] * c[sH];
	w[r233b] = k[r233b] * c[sC2H4];
	w[r235f] = k[r235f] * c[sC2H4] * c[sH];
	w[r235b] = k[r235b] * c[sC2H3] * c[sH2];
	w[r236] = k[r236] * c[sC2H4] * c[sO];
	w[r237] = k[r237] * c[sC2H4] * c[sO];
	w[r238] = k[r238] * c[sC2H4] * c[sOH];
	w[r239f] = k[r239f] * c[sC2H4] * c[sCH3];
	w[r239b] = k[r239b] * c[sC2H3] * c[sCH4];
	w[r240f] = k[r240f] * c[sC2H4] * c[sO2];
	w[r240b] = k[r240b] * c[sC2H3] * c[sHO2];
	w[r241] = k[r241] * c[sC2H4] * c[sCH3O];
	w[r242] = k[r242] * c[sC2H4] * c[sCH3O2];
	w[r244f] = k[r244f] * c[sC2H4] * c[sCH3CO3];
	w[r244b] = k[r244b] * c[sC2H3] * c[sCH3CO3H];
	w[r252] = k[r252] * c[sC2H3] * c[sO2];
	w[r253f] = k[r253f] * c[sC2H3] * c[sO2];
	w[r253b] = k[r253b] * c[sCH2CHO] * c[sO];
	w[r314f] = k[r314f] * c[sC2H3CHO];
	w[r314b] = k[r314b] * c[sC2H3] * c[sHCO];
	w[r315f] = k[r315f] * c[sC2H3CHO] * c[sH];
	w[r315b] = k[r315b] * c[sC2H3CO] * c[sH2];
	w[r316] = k[r316] * c[sC2H3CHO] * c[sO];
	w[r317] = k[r317] * c[sC2H3CHO] * c[sOH];
	w[r318] = k[r318] * c[sC2H3CHO] * c[sO2];
	w[r319f] = k[r319f] * c[sC2H3CHO] * c[sHO2];
	w[r319b] = k[r319b] * c[sC2H3CO] * c[sH2O2];
	w[r320] = k[r320] * c[sC2H3CHO] * c[sCH3];
	w[r322] = k[r322] * c[sC2H3CHO] * c[sCH3O];
	w[r323f] = k[r323f] * c[sC2H3CHO] * c[sCH3O2];
	w[r323b] = k[r323b] * c[sC2H3CO] * c[sCH3O2H];
	w[r324] = k[r324] * c[sC2H3CO];
	w[r400f] = k[r400f] * c[sCH3] * c[sC2H5];
	w[r400b] = k[r400b] * c[sC3H8];
	w[r401f] = k[r401f] * c[sC3H8];
	w[r401b] = k[r401b] * c[sNXC3H7] * c[sH];
	w[r402f] = k[r402f] * c[sC3H8];
	w[r402b] = k[r402b] * c[sIXC3H7] * c[sH];
	w[r403f] = k[r403f] * c[sC3H8] * c[sO2];
	w[r403b] = k[r403b] * c[sIXC3H7] * c[sHO2];
	w[r404f] = k[r404f] * c[sC3H8] * c[sO2];
	w[r404b] = k[r404b] * c[sNXC3H7] * c[sHO2];
	w[r405f] = k[r405f] * c[sH] * c[sC3H8];
	w[r405b] = k[r405b] * c[sH2] * c[sIXC3H7];
	w[r406f] = k[r406f] * c[sH] * c[sC3H8];
	w[r406b] = k[r406b] * c[sH2] * c[sNXC3H7];
	w[r407f] = k[r407f] * c[sC3H8] * c[sO];
	w[r407b] = k[r407b] * c[sIXC3H7] * c[sOH];
	w[r408f] = k[r408f] * c[sC3H8] * c[sO];
	w[r408b] = k[r408b] * c[sNXC3H7] * c[sOH];
	w[r409] = k[r409] * c[sC3H8] * c[sOH];
	w[r410] = k[r410] * c[sC3H8] * c[sOH];
	w[r411f] = k[r411f] * c[sC3H8] * c[sHO2];
	w[r411b] = k[r411b] * c[sIXC3H7] * c[sH2O2];
	w[r412f] = k[r412f] * c[sC3H8] * c[sHO2];
	w[r412b] = k[r412b] * c[sNXC3H7] * c[sH2O2];
	w[r413f] = k[r413f] * c[sCH3] * c[sC3H8];
	w[r413b] = k[r413b] * c[sCH4] * c[sIXC3H7];
	w[r414] = k[r414] * c[sCH3] * c[sC3H8];
	w[r415f] = k[r415f] * c[sIXC3H7] * c[sC3H8];
	w[r415b] = k[r415b] * c[sNXC3H7] * c[sC3H8];
	w[r416] = k[r416] * c[sC2H3] * c[sC3H8];
	w[r419] = k[r419] * c[sC2H5] * c[sC3H8];
	w[r420f] = k[r420f] * c[sC3H8] * c[sC3H5XAXC3H5];
	w[r420b] = k[r420b] * c[sNXC3H7] * c[sC3H6];
	w[r421f] = k[r421f] * c[sC3H8] * c[sC3H5XAXC3H5];
	w[r421b] = k[r421b] * c[sIXC3H7] * c[sC3H6];
	w[r422f] = k[r422f] * c[sC3H8] * c[sCH3O];
	w[r422b] = k[r422b] * c[sNXC3H7] * c[sCH3OH];
	w[r423f] = k[r423f] * c[sC3H8] * c[sCH3O];
	w[r423b] = k[r423b] * c[sIXC3H7] * c[sCH3OH];
	w[r424] = k[r424] * c[sCH3O2] * c[sC3H8];
	w[r425] = k[r425] * c[sCH3O2] * c[sC3H8];
	w[r430] = k[r430] * c[sIXC3H7O2] * c[sC3H8];
	w[r431] = k[r431] * c[sIXC3H7O2] * c[sC3H8];
	w[r432f] = k[r432f] * c[sC3H8] * c[sCH3CO3];
	w[r432b] = k[r432b] * c[sIXC3H7] * c[sCH3CO3H];
	w[r433f] = k[r433f] * c[sC3H8] * c[sCH3CO3];
	w[r433b] = k[r433b] * c[sNXC3H7] * c[sCH3CO3H];
	w[r436f] = k[r436f] * c[sIXC3H7];
	w[r436b] = k[r436b] * c[sH] * c[sC3H6];
	w[r437] = k[r437] * c[sIXC3H7] * c[sH];
	w[r439] = k[r439] * c[sIXC3H7] * c[sOH];
	w[r441f] = k[r441f] * c[sIXC3H7] * c[sO];
	w[r441b] = k[r441b] * c[sCH3CHO] * c[sCH3];
	w[r442f] = k[r442f] * c[sNXC3H7];
	w[r442b] = k[r442b] * c[sCH3] * c[sC2H4];
	w[r443f] = k[r443f] * c[sNXC3H7];
	w[r443b] = k[r443b] * c[sH] * c[sC3H6];
	w[r448f] = k[r448f] * c[sC3H6];
	w[r448b] = k[r448b] * c[sC2H3] * c[sCH3];
	w[r449f] = k[r449f] * c[sC3H6];
	w[r449b] = k[r449b] * c[sC3H5XAXC3H5] * c[sH];
	w[r452] = k[r452] * c[sC3H6] * c[sO];
	w[r455] = k[r455] * c[sC3H6] * c[sO];
	w[r457f] = k[r457f] * c[sC3H6] * c[sO];
	w[r457b] = k[r457b] * c[sC3H5XTXC3H5] * c[sOH];
	w[r458f] = k[r458f] * c[sC3H6] * c[sOH];
	w[r458b] = k[r458b] * c[sC3H5XAXC3H5] * c[sH2O];
	w[r460f] = k[r460f] * c[sC3H6] * c[sOH];
	w[r460b] = k[r460b] * c[sC3H5XTXC3H5] * c[sH2O];
	w[r461] = k[r461] * c[sC3H6] * c[sHO2];
	w[r463f] = k[r463f] * c[sC3H6] * c[sHO2];
	w[r463b] = k[r463b] * c[sC3H5XTXC3H5] * c[sH2O2];
	w[r464f] = k[r464f] * c[sC3H6] * c[sH];
	w[r464b] = k[r464b] * c[sC3H5XAXC3H5] * c[sH2];
	w[r466f] = k[r466f] * c[sC3H6] * c[sH];
	w[r466b] = k[r466b] * c[sC3H5XTXC3H5] * c[sH2];
	w[r467f] = k[r467f] * c[sC3H6] * c[sH];
	w[r467b] = k[r467b] * c[sC2H4] * c[sCH3];
	w[r468f] = k[r468f] * c[sC3H6] * c[sO2];
	w[r468b] = k[r468b] * c[sC3H5XAXC3H5] * c[sHO2];
	w[r470] = k[r470] * c[sC3H6] * c[sO2];
	w[r471] = k[r471] * c[sC3H6] * c[sCH3];
	w[r473] = k[r473] * c[sC3H6] * c[sCH3];
	w[r475] = k[r475] * c[sC3H6] * c[sCH3CO3];
	w[r476] = k[r476] * c[sC3H6] * c[sCH3O2];
	w[r480] = k[r480] * c[sC3H6] * c[sIXC3H7O2];
	w[r486] = k[r486] * c[sC3H5XAXC3H5] * c[sHO2];
	w[r487] = k[r487] * c[sC3H5XAXC3H5] * c[sCH3O2];
	w[r491f] = k[r491f] * c[sC3H5XAXC3H5] * c[sC2H5];
	w[r491b] = k[r491b] * c[sC2H4] * c[sC3H6];
	w[r495f] = k[r495f] * c[sC3H5XAXC3H5] * c[sO2];
	w[r495b] = k[r495b] * c[sCH2CHO] * c[sCH2O];
	w[r496f] = k[r496f] * c[sC3H5XAXC3H5] * c[sO2];
	w[r496b] = k[r496b] * c[sC2H3CHO] * c[sOH];
	w[r507f] = k[r507f] * c[sC3H5XTXC3H5] * c[sO2];
	w[r507b] = k[r507b] * c[sCH3COCH2] * c[sO];
	w[r508] = k[r508] * c[sC3H5XTXC3H5] * c[sO2];
	w[r558] = k[r558] * c[sIXC3H7] * c[sHO2];
	w[r560] = k[r560] * c[sCH3O2] * c[sIXC3H7];
	w[r561f] = k[r561f] * c[sNXC3H7O2];
	w[r561b] = k[r561b] * c[sNXC3H7] * c[sO2];
	w[r562f] = k[r562f] * c[sIXC3H7O2];
	w[r562b] = k[r562b] * c[sIXC3H7] * c[sO2];
	w[r565] = k[r565] * c[sIXC3H7O2] * c[sCH2O];
	w[r566] = k[r566] * c[sIXC3H7O2] * c[sCH3CHO];
	w[r568] = k[r568] * c[sIXC3H7O2] * c[sHO2];
	w[r570] = k[r570] * c[sC2H4] * c[sIXC3H7O2];
	w[r574f] = k[r574f] * c[sC2H3CHO] * c[sIXC3H7O2];
	w[r574b] = k[r574b] * c[sC2H3CO] * c[sIXC3H7O2H];
	w[r576] = k[r576] * c[sCH4] * c[sIXC3H7O2];
	w[r578] = k[r578] * c[sIXC3H7O2] * c[sCH3O2];
	w[r580f] = k[r580f] * c[sH2] * c[sIXC3H7O2];
	w[r580b] = k[r580b] * c[sH] * c[sIXC3H7O2H];
	w[r581f] = k[r581f] * c[sIXC3H7O2] * c[sC2H6];
	w[r581b] = k[r581b] * c[sIXC3H7O2H] * c[sC2H5];
	w[r585] = k[r585] * c[sIXC3H7O2] * c[sCH3CO3];
	w[r589] = k[r589] * c[sIXC3H7O2] * c[sIXC3H7O2];
	w[r592] = k[r592] * c[sIXC3H7O2] * c[sCH3];
	w[r594] = k[r594] * c[sIXC3H7O2] * c[sIXC3H7];
	w[r596] = k[r596] * c[sIXC3H7O2] * c[sC3H5XAXC3H5];
	w[r603] = k[r603] * c[sIXC3H7O2H];
	w[r606f] = k[r606f] * c[sIXC3H7O];
	w[r606b] = k[r606b] * c[sCH3] * c[sCH3CHO];
	w[r609f] = k[r609f] * c[sNXC3H7O2];
	w[r609b] = k[r609b] * c[sC3H6OOH1X2XC3H7O2];
	w[r610f] = k[r610f] * c[sNXC3H7O2];
	w[r610b] = k[r610b] * c[sC3H6OOH1X3XC3H7O2];
	w[r611f] = k[r611f] * c[sIXC3H7O2];
	w[r611b] = k[r611b] * c[sC3H6OOH2X1XC3H7O2];
	w[r616f] = k[r616f] * c[sC3H6OOH1X2XC3H7O2];
	w[r616b] = k[r616b] * c[sC3H6] * c[sHO2];
	w[r617f] = k[r617f] * c[sC3H6OOH2X1XC3H7O2];
	w[r617b] = k[r617b] * c[sC3H6] * c[sHO2];
	w[r618] = k[r618] * c[sC3H6OOH1X3XC3H7O2];
	w[r620] = k[r620] * c[sC3H6OOH1X2XC3H7O2];
	w[r622f] = k[r622f] * c[sC3H6OOH1X2O2XC3H7O4];
	w[r622b] = k[r622b] * c[sC3H6OOH1X2XC3H7O2] * c[sO2];
	w[r623f] = k[r623f] * c[sC3H6OOH1X3O2XC3H7O4];
	w[r623b] = k[r623b] * c[sC3H6OOH1X3XC3H7O2] * c[sO2];
	w[r624f] = k[r624f] * c[sC3H6OOH2X1O2XC3H7O4];
	w[r624b] = k[r624b] * c[sC3H6OOH2X1XC3H7O2] * c[sO2];
	w[r625] = k[r625] * c[sC3H6OOH1X2O2XC3H7O4];
	w[r626] = k[r626] * c[sC3H6OOH1X3O2XC3H7O4];
	w[r627] = k[r627] * c[sC3H6OOH2X1O2XC3H7O4];
	w[r631f] = k[r631f] * c[sC3H52X1J3OOHXC3H7O4];
	w[r631b] = k[r631b] * c[sC3H6OOH1X3O2XC3H7O4];
	w[r632f] = k[r632f] * c[sC3H52X1J3OOHXC3H7O4];
	w[r632b] = k[r632b] * c[sAXC3H5OOH] * c[sHO2];
	w[r633] = k[r633] * c[sC3KET12XC3H6O3];
	w[r634] = k[r634] * c[sC3KET13XC3H6O3];
	w[r635] = k[r635] * c[sC3KET21XC3H6O3];
	w[r636f] = k[r636f] * c[sAXC3H5OOH];
	w[r636b] = k[r636b] * c[sC3H5O] * c[sOH];
	w[r637] = k[r637] * c[sC3H5O];
	w[r638f] = k[r638f] * c[sC3H5O];
	w[r638b] = k[r638b] * c[sC2H3] * c[sCH2O];
	w[r639f] = k[r639f] * c[sC3H5O] * c[sO2];
	w[r639b] = k[r639b] * c[sC2H3CHO] * c[sHO2];


	cdot[sN2] = 0.0;

	cdot[sAR] = 0.0;

	cdot[sH] = - w[r1f] + w[r1b] + w[r2f] - w[r2b]
		+ w[r3f] - w[r3b] + 2 * w[r5f] - 2 * w[r5b]
		+ w[r7f] - w[r7b] + w[r8f] - w[r8b]
		- w[r9f] + w[r9b] - w[r10] - w[r11f]
		+ w[r11b] - w[r17] - w[r18f] + w[r18b]
		+ w[r24f] - w[r24b] + w[r26] - w[r28]
		+ w[r34] + w[r41] - w[r44f] + w[r44b]
		+ w[r63f] - w[r63b] - w[r68] - w[r88f]
		+ w[r88b] - w[r98f] + w[r98b] - w[r99f]
		+ w[r99b] + w[r110] - w[r121f] + w[r121b]
		- w[r145f] + w[r145b] - w[r146f] + w[r146b]
		- w[r156f] + w[r156b] - w[r162f] + w[r162b]
		- w[r163] + w[r164] - w[r194f] + w[r194b]
		- w[r233f] + w[r233b] - w[r235f] + w[r235b]
		+ w[r237] - w[r315f] + w[r315b] + w[r401f]
		- w[r401b] + w[r402f] - w[r402b] - w[r405f]
		+ w[r405b] - w[r406f] + w[r406b] + w[r436f]
		- w[r436b] - w[r437] + w[r443f] - w[r443b]
		+ w[r449f] - w[r449b] - w[r464f] + w[r464b]
		- w[r466f] + w[r466b] - w[r467f] + w[r467b]
		+ w[r580f] - w[r580b] + w[r637];

	cdot[sO2] = - w[r1f] + w[r1b] - w[r6f] + w[r6b]
		- w[r9f] + w[r9b] + w[r10] + w[r12]
		+ w[r13f] - w[r13b] - w[r14f] + w[r14b]
		- w[r15f] + w[r15b] - w[r23f] + w[r23b]
		- w[r27] + w[r33f] - w[r33b] - w[r64]
		+ w[r109f] - w[r109b] - w[r111] - w[r112f]
		+ w[r112b] - w[r113f] + w[r113b] + w[r118]
		+ w[r119] + w[r120] + w[r122] + w[r123]
		- w[r157f] + w[r157b] - w[r158f] + w[r158b]
		- w[r197] + w[r208f] - w[r208b] + w[r209]
		- w[r217] - w[r240f] + w[r240b] - w[r252]
		- w[r253f] + w[r253b] - w[r318] - w[r403f]
		+ w[r403b] - w[r404f] + w[r404b] - w[r468f]
		+ w[r468b] - w[r470] - w[r495f] + w[r495b]
		- w[r496f] + w[r496b] - w[r507f] + w[r507b]
		- w[r508] + w[r561f] - w[r561b] + w[r562f]
		- w[r562b] + w[r568] + w[r578] + w[r585]
		+ w[r589] + w[r622f] - w[r622b] + w[r623f]
		- w[r623b] + w[r624f] - w[r624b] - w[r639f]
		+ w[r639b];

	cdot[sO] = w[r1f] - w[r1b] - w[r2f] + w[r2b]
		- w[r4f] + w[r4b] + 2 * w[r6f] - 2 * w[r6b]
		+ w[r7f] - w[r7b] - w[r12] - w[r19]
		+ w[r23f] - w[r23b] - w[r45] - w[r101f]
		+ w[r101b] - w[r110] + w[r111] - w[r122]
		- w[r147] - w[r164] - w[r236] - w[r237]
		+ w[r253f] - w[r253b] - w[r316] - w[r407f]
		+ w[r407b] - w[r408f] + w[r408b] - w[r441f]
		+ w[r441b] - w[r452] - w[r455] - w[r457f]
		+ w[r457b] + w[r507f] - w[r507b];

	cdot[sOH] = w[r1f] - w[r1b] + w[r2f] - w[r2b]
		- w[r3f] + w[r3b] + 2 * w[r4f] - 2 * w[r4b]
		- w[r7f] + w[r7b] + w[r8f] - w[r8b]
		+ 2 * w[r11f] - 2 * w[r11b] + w[r12] - w[r13f]
		+ w[r13b] + 2 * w[r16f] - 2 * w[r16b] + w[r17]
		+ w[r19] - w[r20] - w[r21f] + w[r21b]
		- w[r24f] + w[r24b] + w[r25] - w[r31]
		+ w[r34] - w[r43] + w[r45] - w[r55]
		+ w[r55] + w[r60] + w[r86f] - w[r86b]
		- w[r91f] + w[r91b] - w[r100f] + w[r100b]
		+ w[r101f] - w[r101b] + w[r108f] - w[r108b]
		+ w[r112f] - w[r112b] + w[r121f] - w[r121b]
		- w[r123] + w[r124] + w[r147] - w[r148]
		+ w[r158f] - w[r158b] - w[r196] - w[r202f]
		+ w[r202b] - w[r203f] + w[r203b] + w[r214]
		+ w[r217] - w[r238] + w[r316] - w[r317]
		+ w[r407f] - w[r407b] + w[r408f] - w[r408b]
		- w[r409] - w[r410] - w[r439] + w[r455]
		+ w[r457f] - w[r457b] - w[r458f] + w[r458b]
		- w[r460f] + w[r460b] + w[r486] + w[r496f]
		- w[r496b] + w[r558] + w[r603] + w[r618]
		+ w[r620] + w[r625] + w[r626] + w[r627]
		+ w[r633] + w[r634] + w[r635] + w[r636f]
		- w[r636b];

	cdot[sH2] = - w[r2f] + w[r2b] - w[r3f] + w[r3b]
		- w[r5f] + w[r5b] + w[r10] + w[r18f]
		- w[r18b] + w[r28] + w[r44f] - w[r44b]
		+ w[r51f] - w[r51b] + w[r68] + w[r88f]
		- w[r88b] + w[r99f] - w[r99b] + w[r146f]
		- w[r146b] + w[r163] + w[r194f] - w[r194b]
		+ w[r235f] - w[r235b] + w[r315f] - w[r315b]
		+ w[r405f] - w[r405b] + w[r406f] - w[r406b]
		+ w[r464f] - w[r464b] + w[r466f] - w[r466b]
		- w[r580f] + w[r580b];

	cdot[sH2O] = w[r3f] - w[r3b] - w[r4f] + w[r4b]
		- w[r8f] + w[r8b] + w[r13f] - w[r13b]
		+ w[r17] + w[r20] + w[r21f] - w[r21b]
		+ w[r31] + w[r43] + w[r50f] - w[r50b]
		+ w[r55] + w[r91f] - w[r91b] + w[r100f]
		- w[r100b] + w[r148] + w[r196] + w[r203f]
		- w[r203b] + w[r238] + w[r317] + w[r409]
		+ w[r410] + w[r439] + w[r458f] - w[r458b]
		+ w[r460f] - w[r460b];

	cdot[sHO2] = w[r9f] - w[r9b] - w[r10] - w[r11f]
		+ w[r11b] - w[r12] - w[r13f] + w[r13b]
		+ 2 * w[r14f] - 2 * w[r14b] + 2 * w[r15f] - 2 * w[r15b]
		+ w[r18f] - w[r18b] + w[r19] + w[r20]
		+ w[r21f] - w[r21b] - w[r25] + w[r27]
		- w[r33f] + w[r33b] - w[r34] - w[r47f]
		+ w[r47b] - w[r60] + w[r64] - w[r69f]
		+ w[r69b] - w[r102f] + w[r102b] - w[r108f]
		+ w[r108b] - w[r109f] + w[r109b] - w[r118]
		- w[r151] + w[r157f] - w[r157b] + w[r197]
		- w[r199] - w[r209] + w[r210f] - w[r210b]
		+ w[r240f] - w[r240b] + w[r318] - w[r319f]
		+ w[r319b] + w[r403f] - w[r403b] + w[r404f]
		- w[r404b] - w[r411f] + w[r411b] - w[r412f]
		+ w[r412b] - w[r461] - w[r463f] + w[r463b]
		+ w[r468f] - w[r468b] + w[r470] - w[r486]
		- w[r558] - w[r568] + w[r616f] - w[r616b]
		+ w[r617f] - w[r617b] + w[r632f] - w[r632b]
		+ w[r639f] - w[r639b];

	cdot[sH2O2] = - w[r14f] + w[r14b] - w[r15f] + w[r15b]
		- w[r16f] + w[r16b] - w[r17] - w[r18f]
		+ w[r18b] - w[r19] - w[r20] - w[r21f]
		+ w[r21b] + w[r47f] - w[r47b] + w[r60]
		+ w[r69f] - w[r69b] + w[r102f] - w[r102b]
		+ w[r151] + w[r199] - w[r210f] + w[r210b]
		+ w[r319f] - w[r319b] + w[r411f] - w[r411b]
		+ w[r412f] - w[r412b] + w[r461] + w[r463f]
		- w[r463b];

	cdot[sCO] = - w[r23f] + w[r23b] - w[r24f] + w[r24b]
		- w[r25] + w[r26] + w[r27] + w[r28]
		+ w[r31] + w[r32] + w[r50f] - w[r50b]
		+ w[r55] + w[r60] + w[r204f] - w[r204b]
		+ w[r217] + w[r324];

	cdot[sCO2] = w[r23f] - w[r23b] + w[r24f] - w[r24b]
		+ w[r25] + w[r34] + w[r51f] - w[r51b]
		+ w[r215f] - w[r215b];

	cdot[sHCO] = - w[r26] - w[r27] - w[r28] - w[r31]
		- w[r32] - w[r33f] + w[r33b] - w[r34]
		+ w[r41] + w[r43] + w[r44f] - w[r44b]
		+ w[r45] + w[r46] + w[r47f] - w[r47b]
		+ w[r65] + w[r114f] - w[r114b] + w[r193f]
		- w[r193b] + w[r212f] - w[r212b] + w[r236]
		+ w[r252] + w[r314f] - w[r314b] + w[r452]
		+ w[r565] + w[r633];

	cdot[sCH3] = - w[r32] - w[r46] - w[r67] + w[r86f]
		- w[r86b] - w[r98f] + w[r98b] + w[r99f]
		- w[r99b] + w[r100f] - w[r100b] + w[r101f]
		- w[r101b] + w[r102f] - w[r102b] - w[r108f]
		+ w[r108b] - w[r109f] + w[r109b] - w[r110]
		- w[r111] - w[r112f] + w[r112b] - w[r113f]
		+ w[r113b] + w[r115f] - w[r115b] - w[r117f]
		+ w[r117b] - 2 * w[r144f] + 2 * w[r144b] - w[r150]
		- w[r161f] + w[r161b] + 2 * w[r162f] - 2 * w[r162b]
		+ w[r193f] - w[r193b] - w[r198f] + w[r198b]
		+ w[r202f] - w[r202b] + w[r204f] - w[r204b]
		+ w[r215f] - w[r215b] + w[r236] - w[r239f]
		+ w[r239b] - w[r320] - w[r400f] + w[r400b]
		- w[r413f] + w[r413b] - w[r414] + w[r437]
		+ w[r441f] - w[r441b] + w[r442f] - w[r442b]
		+ w[r448f] - w[r448b] + w[r467f] - w[r467b]
		- w[r471] - w[r473] + w[r576] - w[r592]
		+ w[r606f] - w[r606b];

	cdot[sCH4] = w[r32] + w[r46] + w[r67] + w[r98f]
		- w[r98b] - w[r99f] + w[r99b] - w[r100f]
		+ w[r100b] - w[r101f] + w[r101b] - w[r102f]
		+ w[r102b] + w[r109f] - w[r109b] - w[r115f]
		+ w[r115b] + w[r150] + w[r161f] - w[r161b]
		+ w[r198f] - w[r198b] + w[r239f] - w[r239b]
		+ w[r320] + w[r413f] - w[r413b] + w[r414]
		+ w[r471] + w[r473] - w[r576];

	cdot[sCH2O] = w[r33f] - w[r33b] - w[r41] - w[r43]
		- w[r44f] + w[r44b] - w[r45] - w[r46]
		- w[r47f] + w[r47b] + w[r63f] - w[r63b]
		+ w[r64] - w[r65] + w[r67] + w[r68]
		+ w[r69f] - w[r69b] - w[r97f] + w[r97b]
		+ w[r110] + w[r112f] - w[r112b] - w[r114f]
		+ w[r114b] + w[r119] - w[r212f] + w[r212b]
		+ w[r217] + w[r252] + w[r495f] - w[r495b]
		+ w[r508] - w[r565] + w[r618] + w[r620]
		+ w[r634] + w[r635] + w[r638f] - w[r638b];

	cdot[sHOCHO] = - w[r50f] + w[r50b] - w[r51f] + w[r51b]
		- w[r55] - w[r60] + w[r202f] - w[r202b];

	cdot[sCH3O] = - w[r63f] + w[r63b] - w[r64] - w[r65]
		- w[r67] - w[r68] - w[r69f] + w[r69b]
		+ w[r88f] - w[r88b] + w[r91f] - w[r91b]
		+ 2 * w[r97f] - 2 * w[r97b] + w[r108f] - w[r108b]
		+ w[r111] + 2 * w[r117f] - 2 * w[r117b] + 2 * w[r120]
		+ w[r121f] - w[r121b] + w[r122] + w[r124]
		- w[r241] - w[r322] - w[r422f] + w[r422b]
		- w[r423f] + w[r423b] + w[r487] + w[r560]
		+ w[r578] + w[r592];

	cdot[sCH3OH] = w[r65] - w[r86f] + w[r86b] - w[r88f]
		+ w[r88b] - w[r91f] + w[r91b] - w[r97f]
		+ w[r97b] + w[r119] + w[r123] + w[r241]
		+ w[r322] + w[r422f] - w[r422b] + w[r423f]
		- w[r423b];

	cdot[sCH3O2] = w[r113f] - w[r113b] - w[r114f] + w[r114b]
		- w[r115f] + w[r115b] - w[r117f] + w[r117b]
		- w[r118] - 2 * w[r119] - 2 * w[r120] - w[r121f]
		+ w[r121b] - w[r122] - w[r123] - w[r152]
		- w[r200f] + w[r200b] - w[r242] - w[r323f]
		+ w[r323b] - w[r424] - w[r425] - w[r476]
		- w[r487] - w[r560] - w[r578];

	cdot[sCH3O2H] = w[r114f] - w[r114b] + w[r115f] - w[r115b]
		+ w[r118] - w[r124] + w[r152] + w[r200f]
		- w[r200b] + w[r242] + w[r323f] - w[r323b]
		+ w[r424] + w[r425] + w[r476];

	cdot[sC2H6] = w[r144f] - w[r144b] + w[r145f] - w[r145b]
		- w[r146f] + w[r146b] - w[r147] - w[r148]
		- w[r150] - w[r151] - w[r152] - w[r213]
		+ w[r419] - w[r581f] + w[r581b];

	cdot[sC2H5] = - w[r145f] + w[r145b] + w[r146f] - w[r146b]
		+ w[r147] + w[r148] + w[r150] + w[r151]
		+ w[r152] + w[r156f] - w[r156b] - w[r157f]
		+ w[r157b] - w[r158f] + w[r158b] - w[r161f]
		+ w[r161b] - w[r162f] + w[r162b] - w[r163]
		- w[r164] + w[r213] - w[r400f] + w[r400b]
		- w[r419] + w[r437] + w[r452] - w[r491f]
		+ w[r491b] + w[r581f] - w[r581b];

	cdot[sC2H4] = - w[r156f] + w[r156b] + w[r157f] - w[r157b]
		+ w[r161f] - w[r161b] + w[r163] + w[r233f]
		- w[r233b] - w[r235f] + w[r235b] - w[r236]
		- w[r237] - w[r238] - w[r239f] + w[r239b]
		- w[r240f] + w[r240b] - w[r241] - w[r242]
		- w[r244f] + w[r244b] + w[r416] + w[r442f]
		- w[r442b] + w[r467f] - w[r467b] + w[r491f]
		- w[r491b] - w[r570] + w[r618] + w[r620];

	cdot[sCH3CHO] = w[r158f] - w[r158b] + w[r164] - w[r193f]
		+ w[r193b] - w[r194f] + w[r194b] - w[r196]
		- w[r197] - w[r198f] + w[r198b] - w[r199]
		- w[r200f] + w[r200b] - w[r201] - w[r202f]
		+ w[r202b] - w[r203f] + w[r203b] + w[r441f]
		- w[r441b] - w[r566] + w[r606f] - w[r606b]
		+ w[r633];

	cdot[sCH3CO] = w[r194f] - w[r194b] + w[r196] + w[r197]
		+ w[r198f] - w[r198b] + w[r199] + w[r200f]
		- w[r200b] + w[r201] - w[r204f] + w[r204b]
		+ w[r208f] - w[r208b] + w[r508] + w[r566]
		+ w[r635];

	cdot[sCH3CO3] = - w[r201] - w[r208f] + w[r208b] - w[r209]
		- w[r210f] + w[r210b] - w[r212f] + w[r212b]
		- w[r213] - w[r244f] + w[r244b] - w[r432f]
		+ w[r432b] - w[r433f] + w[r433b] - w[r475]
		- w[r585];

	cdot[sCH3CO3H] = w[r201] + w[r209] + w[r210f] - w[r210b]
		+ w[r212f] - w[r212b] + w[r213] - w[r214]
		+ w[r244f] - w[r244b] + w[r432f] - w[r432b]
		+ w[r433f] - w[r433b] + w[r475];

	cdot[sCH2CHO] = w[r203f] - w[r203b] - w[r217] + w[r237]
		+ w[r253f] - w[r253b] + w[r495f] - w[r495b]
		+ w[r634];

	cdot[sCH3CO2] = w[r214] - w[r215f] + w[r215b] + w[r585];

	cdot[sC2H3] = - w[r233f] + w[r233b] + w[r235f] - w[r235b]
		+ w[r238] + w[r239f] - w[r239b] + w[r240f]
		- w[r240b] + w[r241] + w[r242] + w[r244f]
		- w[r244b] - w[r252] - w[r253f] + w[r253b]
		+ w[r314f] - w[r314b] + w[r324] - w[r416]
		+ w[r448f] - w[r448b] + w[r570] + w[r638f]
		- w[r638b];

	cdot[sC2H3CHO] = - w[r314f] + w[r314b] - w[r315f] + w[r315b]
		- w[r316] - w[r317] - w[r318] - w[r319f]
		+ w[r319b] - w[r320] - w[r322] - w[r323f]
		+ w[r323b] + w[r496f] - w[r496b] - w[r574f]
		+ w[r574b] + w[r637] + w[r639f] - w[r639b];

	cdot[sC2H3CO] = w[r315f] - w[r315b] + w[r316] + w[r317]
		+ w[r318] + w[r319f] - w[r319b] + w[r320]
		+ w[r322] + w[r323f] - w[r323b] - w[r324]
		+ w[r574f] - w[r574b];

	cdot[sC3H8] = w[r400f] - w[r400b] - w[r401f] + w[r401b]
		- w[r402f] + w[r402b] - w[r403f] + w[r403b]
		- w[r404f] + w[r404b] - w[r405f] + w[r405b]
		- w[r406f] + w[r406b] - w[r407f] + w[r407b]
		- w[r408f] + w[r408b] - w[r409] - w[r410]
		- w[r411f] + w[r411b] - w[r412f] + w[r412b]
		- w[r413f] + w[r413b] - w[r414] - w[r415f]
		+ w[r415f] - w[r415b] + w[r415b] - w[r416]
		- w[r419] - w[r420f] + w[r420b] - w[r421f]
		+ w[r421b] - w[r422f] + w[r422b] - w[r423f]
		+ w[r423b] - w[r424] - w[r425] - w[r430]
		- w[r431] - w[r432f] + w[r432b] - w[r433f]
		+ w[r433b];

	cdot[sNXC3H7] = w[r401f] - w[r401b] + w[r404f] - w[r404b]
		+ w[r406f] - w[r406b] + w[r408f] - w[r408b]
		+ w[r409] + w[r412f] - w[r412b] + w[r414]
		+ w[r415f] - w[r415b] + w[r419] + w[r420f]
		- w[r420b] + w[r422f] - w[r422b] + w[r424]
		+ w[r430] + w[r433f] - w[r433b] - w[r442f]
		+ w[r442b] - w[r443f] + w[r443b] + w[r561f]
		- w[r561b];

	cdot[sIXC3H7] = w[r402f] - w[r402b] + w[r403f] - w[r403b]
		+ w[r405f] - w[r405b] + w[r407f] - w[r407b]
		+ w[r410] + w[r411f] - w[r411b] + w[r413f]
		- w[r413b] - w[r415f] + w[r415b] + w[r416]
		+ w[r421f] - w[r421b] + w[r423f] - w[r423b]
		+ w[r425] + w[r431] + w[r432f] - w[r432b]
		- w[r436f] + w[r436b] - w[r437] - w[r439]
		- w[r441f] + w[r441b] - w[r558] - w[r560]
		+ w[r562f] - w[r562b] - w[r594];

	cdot[sC3H5XAXC3H5] = - w[r420f] + w[r420b] - w[r421f] + w[r421b]
		+ w[r449f] - w[r449b] + w[r455] + w[r458f]
		- w[r458b] + w[r461] + w[r464f] - w[r464b]
		+ w[r468f] - w[r468b] + w[r471] + w[r475]
		+ w[r476] + w[r480] - w[r486] - w[r487]
		- w[r491f] + w[r491b] - w[r495f] + w[r495b]
		- w[r496f] + w[r496b] - w[r596];

	cdot[sC3H6] = w[r420f] - w[r420b] + w[r421f] - w[r421b]
		+ w[r436f] - w[r436b] + w[r439] + w[r443f]
		- w[r443b] - w[r448f] + w[r448b] - w[r449f]
		+ w[r449b] - w[r452] - w[r455] - w[r457f]
		+ w[r457b] - w[r458f] + w[r458b] - w[r460f]
		+ w[r460b] - w[r461] - w[r463f] + w[r463b]
		- w[r464f] + w[r464b] - w[r466f] + w[r466b]
		- w[r467f] + w[r467b] - w[r468f] + w[r468b]
		- w[r470] - w[r471] - w[r473] - w[r475]
		- w[r476] - w[r480] + w[r491f] - w[r491b]
		+ w[r616f] - w[r616b] + w[r617f] - w[r617b];

	cdot[sIXC3H7O2] = - w[r430] - w[r431] - w[r480] - w[r562f]
		+ w[r562b] - w[r565] - w[r566] - w[r568]
		- w[r570] - w[r574f] + w[r574b] - w[r576]
		- w[r578] - w[r580f] + w[r580b] - w[r581f]
		+ w[r581b] - w[r585] - 2 * w[r589] - w[r592]
		- w[r594] - w[r596] - w[r611f] + w[r611b];

	cdot[sIXC3H7O2H] = w[r430] + w[r431] + w[r480] + w[r565]
		+ w[r566] + w[r568] + w[r570] + w[r574f]
		- w[r574b] + w[r576] + w[r580f] - w[r580b]
		+ w[r581f] - w[r581b] - w[r603];

	cdot[sC3H5XTXC3H5] = w[r457f] - w[r457b] + w[r460f] - w[r460b]
		+ w[r463f] - w[r463b] + w[r466f] - w[r466b]
		+ w[r470] + w[r473] - w[r507f] + w[r507b]
		- w[r508];

	cdot[sC3H5O] = w[r486] + w[r487] + w[r596] + w[r636f]
		- w[r636b] - w[r637] - w[r638f] + w[r638b]
		- w[r639f] + w[r639b];

	cdot[sCH3COCH2] = w[r507f] - w[r507b];

	cdot[sIXC3H7O] = w[r558] + w[r560] + w[r578] + w[r585]
		+ 2 * w[r589] + w[r592] + 2 * w[r594] + w[r596]
		+ w[r603] - w[r606f] + w[r606b];

	cdot[sNXC3H7O2] = - w[r561f] + w[r561b] - w[r609f] + w[r609b]
		- w[r610f] + w[r610b];

	cdot[sC3H6OOH1X2XC3H7O2] = w[r609f] - w[r609b] - w[r616f] + w[r616b]
		- w[r620] + w[r622f] - w[r622b];

	cdot[sC3H6OOH1X3XC3H7O2] = w[r610f] - w[r610b] - w[r618] + w[r623f]
		- w[r623b];

	cdot[sC3H6OOH2X1XC3H7O2] = w[r611f] - w[r611b] - w[r617f] + w[r617b]
		+ w[r624f] - w[r624b];

	cdot[sC3H6OOH1X2O2XC3H7O4] = - w[r622f] + w[r622b] - w[r625];

	cdot[sC3H6OOH1X3O2XC3H7O4] = - w[r623f] + w[r623b] - w[r626] + w[r631f]
		- w[r631b];

	cdot[sC3H6OOH2X1O2XC3H7O4] = - w[r624f] + w[r624b] - w[r627];

	cdot[sC3KET12XC3H6O3] = w[r625] - w[r633];

	cdot[sC3KET13XC3H6O3] = w[r626] - w[r634];

	cdot[sC3KET21XC3H6O3] = w[r627] - w[r635];

	cdot[sC3H52X1J3OOHXC3H7O4] = - w[r631f] + w[r631b] - w[r632f] + w[r632b];

	cdot[sAXC3H5OOH] = w[r632f] - w[r632b] - w[r636f] + w[r636b];

}

double GetLindRateCoeff( double temp, double pressure
				, double k0, double kInf
				, double fc, double conc )
{
	const double	R = 8314.34;   /* [J / kmole K] */
	double			Ntmp;
	double			kl;
	double			f;
	double	cCoeff, dCoeff, log10kNull;
	/* double			conc; */

/* 	Flag			iTroe = TRUE; */

	if ( conc <= 0.0 ) {
		conc = pressure / ( R * temp );
	}
	Ntmp = 0.75 - 1.27 * log10( fc );
	if ( 1 ) {	
		cCoeff = - 0.4 - 0.67 * log10( fc );
		dCoeff = 0.14;
		//k0 *= conc / MAX_D(kInf, 1.0e-60);
		k0 *= conc / kInf;
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
	W[sH] =  1.00800000e+00;
	W[sO2] =  3.20000000e+01;
	W[sO] =  1.60000000e+01;
	W[sOH] =  1.70080000e+01;
	W[sH2] =  2.01600000e+00;
	W[sH2O] =  1.80160000e+01;
	W[sHO2] =  3.30080000e+01;
	W[sH2O2] =  3.40160000e+01;
	W[sCO] =  2.80100000e+01;
	W[sCO2] =  4.40100000e+01;
	W[sHCO] =  2.90180000e+01;
	W[sCH3] =  1.50340000e+01;
	W[sCH4] =  1.60420000e+01;
	W[sCH2O] =  3.00260000e+01;
	W[sHOCHO] =  4.60260000e+01;
	W[sCH3O] =  3.10340000e+01;
	W[sCH3OH] =  3.20420000e+01;
	W[sCH3O2] =  4.70340000e+01;
	W[sCH3O2H] =  4.80420000e+01;
	W[sC2H6] =  3.00680000e+01;
	W[sC2H5] =  2.90600000e+01;
	W[sC2H4] =  2.80520000e+01;
	W[sCH3CHO] =  4.40520000e+01;
	W[sCH3CO] =  4.30440000e+01;
	W[sCH3CO3] =  7.50440000e+01;
	W[sCH3CO3H] =  7.60520000e+01;
	W[sCH2CHO] =  4.30440000e+01;
	W[sCH3CO2] =  5.90440000e+01;
	W[sC2H3] =  2.70440000e+01;
	W[sC2H3CHO] =  5.60620000e+01;
	W[sC2H3CO] =  5.50540000e+01;
	W[sC3H8] =  4.40940000e+01;
	W[sNXC3H7] =  4.30860000e+01;
	W[sIXC3H7] =  4.30860000e+01;
	W[sC3H5XAXC3H5] =  4.10700000e+01;
	W[sC3H6] =  4.20780000e+01;
	W[sIXC3H7O2] =  7.50860000e+01;
	W[sIXC3H7O2H] =  7.60940000e+01;
	W[sC3H5XTXC3H5] =  4.10700000e+01;
	W[sC3H5O] =  5.70700000e+01;
	W[sCH3COCH2] =  5.70700000e+01;
	W[sIXC3H7O] =  5.90860000e+01;
	W[sNXC3H7O2] =  7.50860000e+01;
	W[sC3H6OOH1X2XC3H7O2] =  7.50860000e+01;
	W[sC3H6OOH1X3XC3H7O2] =  7.50860000e+01;
	W[sC3H6OOH2X1XC3H7O2] =  7.50860000e+01;
	W[sC3H6OOH1X2O2XC3H7O4] =  1.07086000e+02;
	W[sC3H6OOH1X3O2XC3H7O4] =  1.07086000e+02;
	W[sC3H6OOH2X1O2XC3H7O4] =  1.07086000e+02;
	W[sC3KET12XC3H6O3] =  9.00780000e+01;
	W[sC3KET13XC3H6O3] =  9.00780000e+01;
	W[sC3KET21XC3H6O3] =  9.00780000e+01;
	W[sC3H52X1J3OOHXC3H7O4] =  1.07086000e+02;
	W[sAXC3H5OOH] =  7.40780000e+01;
}

void GetSpeciesNames( char **names )
{
	strcpy( names[sN2], "N2                  " );
	strcpy( names[sAR], "AR                  " );
	strcpy( names[sH], "H                   " );
	strcpy( names[sO2], "O2                  " );
	strcpy( names[sO], "O                   " );
	strcpy( names[sOH], "OH                  " );
	strcpy( names[sH2], "H2                  " );
	strcpy( names[sH2O], "H2O                 " );
	strcpy( names[sHO2], "HO2                 " );
	strcpy( names[sH2O2], "H2O2                " );
	strcpy( names[sCO], "CO                  " );
	strcpy( names[sCO2], "CO2                 " );
	strcpy( names[sHCO], "HCO                 " );
	strcpy( names[sCH3], "CH3                 " );
	strcpy( names[sCH4], "CH4                 " );
	strcpy( names[sCH2O], "CH2O                " );
	strcpy( names[sHOCHO], "HOCHO               " );
	strcpy( names[sCH3O], "CH3O                " );
	strcpy( names[sCH3OH], "CH3OH               " );
	strcpy( names[sCH3O2], "CH3O2               " );
	strcpy( names[sCH3O2H], "CH3O2H              " );
	strcpy( names[sC2H6], "C2H6                " );
	strcpy( names[sC2H5], "C2H5                " );
	strcpy( names[sC2H4], "C2H4                " );
	strcpy( names[sCH3CHO], "CH3CHO              " );
	strcpy( names[sCH3CO], "CH3CO               " );
	strcpy( names[sCH3CO3], "CH3CO3              " );
	strcpy( names[sCH3CO3H], "CH3CO3H             " );
	strcpy( names[sCH2CHO], "CH2CHO              " );
	strcpy( names[sCH3CO2], "CH3CO2              " );
	strcpy( names[sC2H3], "C2H3                " );
	strcpy( names[sC2H3CHO], "C2H3CHO             " );
	strcpy( names[sC2H3CO], "C2H3CO              " );
	strcpy( names[sC3H8], "C3H8                " );
	strcpy( names[sNXC3H7], "N-C3H7              " );
	strcpy( names[sIXC3H7], "I-C3H7              " );
	strcpy( names[sC3H5XAXC3H5], "C3H5-A-C3H5         " );
	strcpy( names[sC3H6], "C3H6                " );
	strcpy( names[sIXC3H7O2], "I-C3H7O2            " );
	strcpy( names[sIXC3H7O2H], "I-C3H7O2H           " );
	strcpy( names[sC3H5XTXC3H5], "C3H5-T-C3H5         " );
	strcpy( names[sC3H5O], "C3H5O               " );
	strcpy( names[sCH3COCH2], "CH3COCH2            " );
	strcpy( names[sIXC3H7O], "I-C3H7O             " );
	strcpy( names[sNXC3H7O2], "N-C3H7O2            " );
	strcpy( names[sC3H6OOH1X2XC3H7O2], "C3H6OOH1-2-C3H7O2   " );
	strcpy( names[sC3H6OOH1X3XC3H7O2], "C3H6OOH1-3-C3H7O2   " );
	strcpy( names[sC3H6OOH2X1XC3H7O2], "C3H6OOH2-1-C3H7O2   " );
	strcpy( names[sC3H6OOH1X2O2XC3H7O4], "C3H6OOH1-2O2-C3H7O4 " );
	strcpy( names[sC3H6OOH1X3O2XC3H7O4], "C3H6OOH1-3O2-C3H7O4 " );
	strcpy( names[sC3H6OOH2X1O2XC3H7O4], "C3H6OOH2-1O2-C3H7O4 " );
	strcpy( names[sC3KET12XC3H6O3], "C3KET12-C3H6O3      " );
	strcpy( names[sC3KET13XC3H6O3], "C3KET13-C3H6O3      " );
	strcpy( names[sC3KET21XC3H6O3], "C3KET21-C3H6O3      " );
	strcpy( names[sC3H52X1J3OOHXC3H7O4], "C3H52-1J3OOH-C3H7O4 " );
	strcpy( names[sAXC3H5OOH], "A-C3H5OOH           " );
//	strcpy( names[sT], "Temp                " );
//	strcpy( names[sR], "Density             " );
//	strcpy( names[sU], "U-Velocity          " );
//	strcpy( names[sV], "V-Velocity          " );
}



void ComputeThermoData( double *h, double *cp, double T )
{
/*
	This function computes enthalpy 'h' and heat capacity 'cp' as
	function of temperature 'T' for all non steady state species
	in units [J/kg] and [J/kg K], respectively.
	The parameter h and cp should provide workspace of length 56 */


	if ( T > 1000.0 ) {
		h[sN2] =  2.96728765e+02 * (
			T * (  2.92664000e+00 + T * (  7.43988500e-04
			+ T * ( -1.89492033e-07 + T * (  2.52426000e-11
			+ T * -1.35067020e-15 ) ) ) ) -9.22797700e+02 );
		cp[sN2] =  2.96728765e+02 * (
			 2.92664000e+00 + T * (  1.48797700e-03 
			+ T * ( -5.68476100e-07 + T * (  1.00970400e-10
			+ T * -6.75335100e-15 ) ) ) );
		h[sAR] =  2.08129068e+02 * (
			T * (  2.50000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) -7.45375000e+02 );
		cp[sAR] =  2.08129068e+02 * (
			 2.50000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sH] =  8.24835317e+03 * (
			T * (  2.50000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  2.54716300e+04 );
		cp[sH] =  8.24835317e+03 * (
			 2.50000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sO2] =  2.59823125e+02 * (
			T * (  3.69757800e+00 + T * (  3.06759850e-04
			+ T * ( -4.19614000e-08 + T * (  4.43820250e-12
			+ T * -2.27287000e-16 ) ) ) ) -1.23393000e+03 );
		cp[sO2] =  2.59823125e+02 * (
			 3.69757800e+00 + T * (  6.13519700e-04 
			+ T * ( -1.25884200e-07 + T * (  1.77528100e-11
			+ T * -1.13643500e-15 ) ) ) );
		h[sO] =  5.19646250e+02 * (
			T * (  2.54206000e+00 + T * ( -1.37753100e-05
			+ T * ( -1.03426767e-09 + T * (  1.13776675e-12
			+ T * -8.73610400e-17 ) ) ) ) +  2.92308000e+04 );
		cp[sO] =  5.19646250e+02 * (
			 2.54206000e+00 + T * ( -2.75506200e-05 
			+ T * ( -3.10280300e-09 + T * (  4.55106700e-12
			+ T * -4.36805200e-16 ) ) ) );
		h[sOH] =  4.88848777e+02 * (
			T * (  3.72342593e+00 + T * ( -3.21499498e-04
			+ T * (  2.67673786e-07 + T * ( -4.73277323e-11
			+ T *  0.00000000e+00 ) ) ) ) +  3.37726070e+03 );
		cp[sOH] =  4.88848777e+02 * (
			 3.72342593e+00 + T * ( -6.42998996e-04 
			+ T * (  8.03021358e-07 + T * ( -1.89310929e-10
			+ T *  0.00000000e+00 ) ) ) );
		h[sH2] =  4.12417659e+03 * (
			T * (  2.99142300e+00 + T * (  3.50032200e-04
			+ T * ( -1.87794300e-08 + T * ( -2.30789450e-12
			+ T *  3.16550400e-16 ) ) ) ) -8.35034000e+02 );
		cp[sH2] =  4.12417659e+03 * (
			 2.99142300e+00 + T * (  7.00064400e-04 
			+ T * ( -5.63382900e-08 + T * ( -9.23157800e-12
			+ T *  1.58275200e-15 ) ) ) );
		h[sH2O] =  4.61497558e+02 * (
			T * (  2.67214600e+00 + T * (  1.52814650e-03
			+ T * ( -2.91008667e-07 + T * (  3.00249000e-11
			+ T * -1.27832360e-15 ) ) ) ) -2.98992100e+04 );
		cp[sH2O] =  4.61497558e+02 * (
			 2.67214600e+00 + T * (  3.05629300e-03 
			+ T * ( -8.73026000e-07 + T * (  1.20099600e-10
			+ T * -6.39161800e-15 ) ) ) );
		h[sHO2] =  2.51888633e+02 * (
			T * (  4.01721090e+00 + T * (  1.11991006e-03
			+ T * ( -2.11219383e-07 + T * (  2.85615925e-11
			+ T * -2.15817070e-15 ) ) ) ) +  1.11856713e+02 );
		cp[sHO2] =  2.51888633e+02 * (
			 4.01721090e+00 + T * (  2.23982013e-03 
			+ T * ( -6.33658150e-07 + T * (  1.14246370e-10
			+ T * -1.07908535e-14 ) ) ) );
		h[sH2O2] =  2.44424389e+02 * (
			T * (  4.57316700e+00 + T * (  2.16806800e-03
			+ T * ( -4.91563000e-07 + T * (  5.87226000e-11
			+ T * -2.86330800e-15 ) ) ) ) -1.80069600e+04 );
		cp[sH2O2] =  2.44424389e+02 * (
			 4.57316700e+00 + T * (  4.33613600e-03 
			+ T * ( -1.47468900e-06 + T * (  2.34890400e-10
			+ T * -1.43165400e-14 ) ) ) );
		h[sCO] =  2.96834702e+02 * (
			T * (  3.02507800e+00 + T * (  7.21344500e-04
			+ T * ( -1.87694267e-07 + T * (  2.54645250e-11
			+ T * -1.38219040e-15 ) ) ) ) -1.42683500e+04 );
		cp[sCO] =  2.96834702e+02 * (
			 3.02507800e+00 + T * (  1.44268900e-03 
			+ T * ( -5.63082800e-07 + T * (  1.01858100e-10
			+ T * -6.91095200e-15 ) ) ) );
		h[sCO2] =  1.88919337e+02 * (
			T * (  4.45362300e+00 + T * (  1.57008450e-03
			+ T * ( -4.26137000e-07 + T * (  5.98499250e-11
			+ T * -3.33806600e-15 ) ) ) ) -4.89669600e+04 );
		cp[sCO2] =  1.88919337e+02 * (
			 4.45362300e+00 + T * (  3.14016900e-03 
			+ T * ( -1.27841100e-06 + T * (  2.39399700e-10
			+ T * -1.66903300e-14 ) ) ) );
		h[sHCO] =  2.86523537e+02 * (
			T * (  3.55727100e+00 + T * (  1.67278650e-03
			+ T * ( -4.45002000e-07 + T * (  6.17643250e-11
			+ T * -3.42770200e-15 ) ) ) ) +  3.91632400e+03 );
		cp[sHCO] =  2.86523537e+02 * (
			 3.55727100e+00 + T * (  3.34557300e-03 
			+ T * ( -1.33500600e-06 + T * (  2.47057300e-10
			+ T * -1.71385100e-14 ) ) ) );
		h[sCH3] =  5.53035786e+02 * (
			T * (  2.84405200e+00 + T * (  3.06898700e-03
			+ T * ( -7.43448333e-07 + T * (  9.46290250e-11
			+ T * -4.90431800e-15 ) ) ) ) +  1.64378100e+04 );
		cp[sCH3] =  5.53035786e+02 * (
			 2.84405200e+00 + T * (  6.13797400e-03 
			+ T * ( -2.23034500e-06 + T * (  3.78516100e-10
			+ T * -2.45215900e-14 ) ) ) );
		h[sCH4] =  5.18285750e+02 * (
			T * (  1.68347900e+00 + T * (  5.11862000e-03
			+ T * ( -1.29170967e-06 + T * (  1.69639625e-10
			+ T * -9.00684600e-15 ) ) ) ) -1.00807900e+04 );
		cp[sCH4] =  5.18285750e+02 * (
			 1.68347900e+00 + T * (  1.02372400e-02 
			+ T * ( -3.87512900e-06 + T * (  6.78558500e-10
			+ T * -4.50342300e-14 ) ) ) );
		h[sCH2O] =  2.76904683e+02 * (
			T * (  2.99560600e+00 + T * (  3.34066050e-03
			+ T * ( -8.76318333e-07 + T * (  1.18428825e-10
			+ T * -6.42503400e-15 ) ) ) ) -1.53203700e+04 );
		cp[sCH2O] =  2.76904683e+02 * (
			 2.99560600e+00 + T * (  6.68132100e-03 
			+ T * ( -2.62895500e-06 + T * (  4.73715300e-10
			+ T * -3.21251700e-14 ) ) ) );
		h[sHOCHO] =  1.80644418e+02 * (
			T * (  2.71806922e+00 + T * (  6.15775995e-03
			+ T * ( -2.00992155e-06 + T * (  2.65342645e-10
			+ T *  0.00000000e+00 ) ) ) ) -4.67812848e+04 );
		cp[sHOCHO] =  1.80644418e+02 * (
			 2.71806922e+00 + T * (  1.23155199e-02 
			+ T * ( -6.02976465e-06 + T * (  1.06137058e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH3O] =  2.67910679e+02 * (
			T * (  3.77080000e+00 + T * (  3.93574850e-03
			+ T * ( -8.85461333e-07 + T * (  9.86107750e-11
			+ T * -4.22523200e-15 ) ) ) ) +  1.27832500e+02 );
		cp[sCH3O] =  2.67910679e+02 * (
			 3.77080000e+00 + T * (  7.87149700e-03 
			+ T * ( -2.65638400e-06 + T * (  3.94443100e-10
			+ T * -2.11261600e-14 ) ) ) );
		h[sCH3OH] =  2.59482554e+02 * (
			T * (  4.02906100e+00 + T * (  4.68829650e-03
			+ T * ( -1.01675133e-06 + T * (  1.08969825e-10
			+ T * -4.44944600e-15 ) ) ) ) -2.61579100e+04 );
		cp[sCH3OH] =  2.59482554e+02 * (
			 4.02906100e+00 + T * (  9.37659300e-03 
			+ T * ( -3.05025400e-06 + T * (  4.35879300e-10
			+ T * -2.22472300e-14 ) ) ) );
		h[sCH3O2] =  1.76772973e+02 * (
			T * (  3.72224006e+00 + T * (  5.77152465e-03
			+ T * ( -1.51686254e-06 + T * (  1.67577132e-10
			+ T *  0.00000000e+00 ) ) ) ) -5.34692750e+02 );
		cp[sCH3O2] =  1.76772973e+02 * (
			 3.72224006e+00 + T * (  1.15430493e-02 
			+ T * ( -4.55058763e-06 + T * (  6.70308528e-10
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH3O2H] =  1.73063986e+02 * (
			T * (  4.46640989e+00 + T * (  7.55106630e-03
			+ T * ( -2.26701845e-06 + T * (  2.83131345e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.80236075e+04 );
		cp[sCH3O2H] =  1.73063986e+02 * (
			 4.46640989e+00 + T * (  1.51021326e-02 
			+ T * ( -6.80105535e-06 + T * (  1.13252538e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC2H6] =  2.76517893e+02 * (
			T * (  2.78247515e-01 + T * (  1.15844838e-02
			+ T * ( -3.41064787e-06 + T * (  4.19248350e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.13124605e+04 );
		cp[sC2H6] =  2.76517893e+02 * (
			 2.78247515e-01 + T * (  2.31689676e-02 
			+ T * ( -1.02319436e-05 + T * (  1.67699340e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC2H5] =  2.86109429e+02 * (
			T * (  7.19048000e+00 + T * (  3.24203850e-03
			+ T * ( -2.14268833e-07 + T * ( -5.86969750e-11
			+ T *  7.76175400e-15 ) ) ) ) +  1.06745500e+04 );
		cp[sC2H5] =  2.86109429e+02 * (
			 7.19048000e+00 + T * (  6.48407700e-03 
			+ T * ( -6.42806500e-07 + T * ( -2.34787900e-10
			+ T *  3.88087700e-14 ) ) ) );
		h[sC2H4] =  2.96390275e+02 * (
			T * (  3.52841900e+00 + T * (  5.74259000e-03
			+ T * ( -1.47279500e-06 + T * (  1.96115025e-10
			+ T * -1.05336960e-14 ) ) ) ) +  4.42828900e+03 );
		cp[sC2H4] =  2.96390275e+02 * (
			 3.52841900e+00 + T * (  1.14851800e-02 
			+ T * ( -4.41838500e-06 + T * (  7.84460100e-10
			+ T * -5.26684800e-14 ) ) ) );
		h[sCH3CHO] =  1.88739217e+02 * (
			T * (  5.86865000e+00 + T * (  5.39712000e-03
			+ T * ( -1.21517667e-06 + T * (  1.35322800e-10
			+ T * -5.79368800e-15 ) ) ) ) -2.26456900e+04 );
		cp[sCH3CHO] =  1.88739217e+02 * (
			 5.86865000e+00 + T * (  1.07942400e-02 
			+ T * ( -3.64553000e-06 + T * (  5.41291200e-10
			+ T * -2.89684400e-14 ) ) ) );
		h[sCH3CO] =  1.93159093e+02 * (
			T * (  5.61227900e+00 + T * (  4.22494300e-03
			+ T * ( -9.51382333e-07 + T * (  1.05959400e-10
			+ T * -4.53680800e-15 ) ) ) ) -5.18786300e+03 );
		cp[sCH3CO] =  1.93159093e+02 * (
			 5.61227900e+00 + T * (  8.44988600e-03 
			+ T * ( -2.85414700e-06 + T * (  4.23837600e-10
			+ T * -2.26840400e-14 ) ) ) );
		h[sCH3CO3] =  1.10792868e+02 * (
			T * (  6.99421069e+00 + T * (  7.88423790e-03
			+ T * ( -2.36106506e-06 + T * (  2.92470817e-10
			+ T *  0.00000000e+00 ) ) ) ) -2.42270741e+04 );
		cp[sCH3CO3] =  1.10792868e+02 * (
			 6.99421069e+00 + T * (  1.57684758e-02 
			+ T * ( -7.08319519e-06 + T * (  1.16988327e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH3CO3H] =  1.09324410e+02 * (
			T * (  6.38601528e+00 + T * (  1.02502835e-02
			+ T * ( -3.24616256e-06 + T * (  4.21010380e-10
			+ T *  0.00000000e+00 ) ) ) ) -4.34818847e+04 );
		cp[sCH3CO3H] =  1.09324410e+02 * (
			 6.38601528e+00 + T * (  2.05005671e-02 
			+ T * ( -9.73848768e-06 + T * (  1.68404152e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH2CHO] =  1.93159093e+02 * (
			T * (  5.97567000e+00 + T * (  4.06529550e-03
			+ T * ( -9.14541333e-07 + T * (  1.01757600e-10
			+ T * -4.35203400e-15 ) ) ) ) +  4.90321800e+02 );
		cp[sCH2CHO] =  1.93159093e+02 * (
			 5.97567000e+00 + T * (  8.13059100e-03 
			+ T * ( -2.74362400e-06 + T * (  4.07030400e-10
			+ T * -2.17601700e-14 ) ) ) );
		h[sCH3CO2] =  1.40816002e+02 * (
			T * (  4.04855337e+00 + T * (  8.20920375e-03
			+ T * ( -2.51916794e-06 + T * (  3.21024970e-10
			+ T *  0.00000000e+00 ) ) ) ) -2.78926423e+04 );
		cp[sCH3CO2] =  1.40816002e+02 * (
			 4.04855337e+00 + T * (  1.64184075e-02 
			+ T * ( -7.55750381e-06 + T * (  1.28409988e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC2H3] =  3.07437509e+02 * (
			T * (  5.93346800e+00 + T * (  2.00887300e-03
			+ T * ( -1.32224667e-07 + T * ( -3.60316750e-11
			+ T *  4.75728800e-15 ) ) ) ) +  3.18543500e+04 );
		cp[sC2H3] =  3.07437509e+02 * (
			 5.93346800e+00 + T * (  4.01774600e-03 
			+ T * ( -3.96674000e-07 + T * ( -1.44126700e-10
			+ T *  2.37864400e-14 ) ) ) );
		h[sC2H3CHO] =  1.48306161e+02 * (
			T * (  5.56154592e+00 + T * (  8.96479185e-03
			+ T * ( -2.67821586e-06 + T * (  3.30738438e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.29035886e+04 );
		cp[sC2H3CHO] =  1.48306161e+02 * (
			 5.56154592e+00 + T * (  1.79295837e-02 
			+ T * ( -8.03464758e-06 + T * (  1.32295375e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC2H3CO] =  1.51021542e+02 * (
			T * (  7.49438492e+00 + T * (  5.37732150e-03
			+ T * ( -1.33085471e-06 + T * (  1.37113989e-10
			+ T *  0.00000000e+00 ) ) ) ) +  2.82679210e+03 );
		cp[sC2H3CO] =  1.51021542e+02 * (
			 7.49438492e+00 + T * (  1.07546430e-02 
			+ T * ( -3.99256412e-06 + T * (  5.48455956e-10
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3H8] =  1.88559441e+02 * (
			T * (  1.23646898e+00 + T * (  1.58034151e-02
			+ T * ( -4.69289980e-06 + T * (  5.81750253e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.44967440e+04 );
		cp[sC3H8] =  1.88559441e+02 * (
			 1.23646898e+00 + T * (  3.16068302e-02 
			+ T * ( -1.40786994e-05 + T * (  2.32700101e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sNXC3H7] =  1.92970803e+02 * (
			T * (  1.77024723e+00 + T * (  1.38805670e-02
			+ T * ( -4.12295397e-06 + T * (  5.10936707e-10
			+ T *  0.00000000e+00 ) ) ) ) +  1.01768072e+04 );
		cp[sNXC3H7] =  1.92970803e+02 * (
			 1.77024723e+00 + T * (  2.77611340e-02 
			+ T * ( -1.23688619e-05 + T * (  2.04374683e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIXC3H7] =  1.92970803e+02 * (
			T * ( -3.55917828e-01 + T * (  1.56030339e-02
			+ T * ( -4.80156740e-06 + T * (  6.13303987e-10
			+ T *  0.00000000e+00 ) ) ) ) +  9.50205709e+03 );
		cp[sIXC3H7] =  1.92970803e+02 * (
			-3.55917828e-01 + T * (  3.12060678e-02 
			+ T * ( -1.44047022e-05 + T * (  2.45321595e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3H5XAXC3H5] =  2.02443146e+02 * (
			T * (  3.66851913e+00 + T * (  9.81016580e-03
			+ T * ( -2.84310523e-06 + T * (  3.45283085e-10
			+ T *  0.00000000e+00 ) ) ) ) +  1.83798091e+04 );
		cp[sC3H5XAXC3H5] =  2.02443146e+02 * (
			 3.66851913e+00 + T * (  1.96203316e-02 
			+ T * ( -8.52931569e-06 + T * (  1.38113234e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3H6] =  1.97593517e+02 * (
			T * (  1.68134957e+00 + T * (  1.24845359e-02
			+ T * ( -3.70920137e-06 + T * (  4.61033815e-10
			+ T *  0.00000000e+00 ) ) ) ) +  7.40726095e+02 );
		cp[sC3H6] =  1.97593517e+02 * (
			 1.68134957e+00 + T * (  2.49690718e-02 
			+ T * ( -1.11276041e-05 + T * (  1.84413526e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIXC3H7O2] =  1.10730895e+02 * (
			T * (  6.12023736e+00 + T * (  1.42615975e-02
			+ T * ( -4.10011727e-06 + T * (  4.89504662e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.13282549e+04 );
		cp[sIXC3H7O2] =  1.10730895e+02 * (
			 6.12023736e+00 + T * (  2.85231950e-02 
			+ T * ( -1.23003518e-05 + T * (  1.95801865e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIXC3H7O2H] =  1.09264068e+02 * (
			T * (  6.94298974e+00 + T * (  1.59467816e-02
			+ T * ( -4.80474350e-06 + T * (  5.97383883e-10
			+ T *  0.00000000e+00 ) ) ) ) -2.86314978e+04 );
		cp[sIXC3H7O2H] =  1.09264068e+02 * (
			 6.94298974e+00 + T * (  3.18935633e-02 
			+ T * ( -1.44142305e-05 + T * (  2.38953553e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3H5XTXC3H5] =  2.02443146e+02 * (
			T * (  1.94476358e+00 + T * (  1.07086382e-02
			+ T * ( -3.18351890e-06 + T * (  3.95877935e-10
			+ T *  0.00000000e+00 ) ) ) ) +  2.95295307e+04 );
		cp[sC3H5XTXC3H5] =  2.02443146e+02 * (
			 1.94476358e+00 + T * (  2.14172765e-02 
			+ T * ( -9.55055670e-06 + T * (  1.58351174e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3H5O] =  1.45686701e+02 * (
			T * (  3.39074577e+00 + T * (  1.20650810e-02
			+ T * ( -3.78836313e-06 + T * (  4.94752345e-10
			+ T *  0.00000000e+00 ) ) ) ) +  9.00757452e+03 );
		cp[sC3H5O] =  1.45686701e+02 * (
			 3.39074577e+00 + T * (  2.41301620e-02 
			+ T * ( -1.13650894e-05 + T * (  1.97900938e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH3COCH2] =  1.45686701e+02 * (
			T * (  4.27637519e+00 + T * (  1.11174643e-02
			+ T * ( -3.36394730e-06 + T * (  4.22492645e-10
			+ T *  0.00000000e+00 ) ) ) ) -5.97887937e+03 );
		cp[sCH3COCH2] =  1.45686701e+02 * (
			 4.27637519e+00 + T * (  2.22349285e-02 
			+ T * ( -1.00918419e-05 + T * (  1.68997058e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIXC3H7O] =  1.40715906e+02 * (
			T * (  4.63651596e+00 + T * (  1.37315629e-02
			+ T * ( -4.00138490e-06 + T * (  4.86885295e-10
			+ T *  0.00000000e+00 ) ) ) ) -9.57914716e+03 );
		cp[sIXC3H7O] =  1.40715906e+02 * (
			 4.63651596e+00 + T * (  2.74631259e-02 
			+ T * ( -1.20041547e-05 + T * (  1.94754118e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sNXC3H7O2] =  1.10730895e+02 * (
			T * (  4.71102090e+00 + T * (  1.53071024e-02
			+ T * ( -4.48109563e-06 + T * (  5.43703638e-10
			+ T *  0.00000000e+00 ) ) ) ) -8.53286429e+03 );
		cp[sNXC3H7O2] =  1.10730895e+02 * (
			 4.71102090e+00 + T * (  3.06142048e-02 
			+ T * ( -1.34432869e-05 + T * (  2.17481455e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3H6OOH1X2XC3H7O2] =  1.10730895e+02 * (
			T * (  6.11209376e+00 + T * (  1.35997173e-02
			+ T * ( -3.93796473e-06 + T * (  4.82472052e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.58345092e+03 );
		cp[sC3H6OOH1X2XC3H7O2] =  1.10730895e+02 * (
			 6.11209376e+00 + T * (  2.71994345e-02 
			+ T * ( -1.18138942e-05 + T * (  1.92988821e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3H6OOH1X3XC3H7O2] =  1.10730895e+02 * (
			T * (  6.08114425e+00 + T * (  1.50669910e-02
			+ T * ( -4.62095230e-06 + T * (  5.82264653e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.16933310e+03 );
		cp[sC3H6OOH1X3XC3H7O2] =  1.10730895e+02 * (
			 6.08114425e+00 + T * (  3.01339820e-02 
			+ T * ( -1.38628569e-05 + T * (  2.32905861e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3H6OOH2X1XC3H7O2] =  1.10730895e+02 * (
			T * (  9.51872798e+00 + T * (  1.12169631e-02
			+ T * ( -3.08417739e-06 + T * (  3.61057348e-10
			+ T *  0.00000000e+00 ) ) ) ) -3.66915364e+03 );
		cp[sC3H6OOH2X1XC3H7O2] =  1.10730895e+02 * (
			 9.51872798e+00 + T * (  2.24339261e-02 
			+ T * ( -9.25253217e-06 + T * (  1.44422939e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3H6OOH1X2O2XC3H7O4] =  7.76417085e+01 * (
			T * (  1.05243923e+01 + T * (  1.53970017e-02
			+ T * ( -4.60295637e-06 + T * (  5.64612330e-10
			+ T *  0.00000000e+00 ) ) ) ) -2.27195952e+04 );
		cp[sC3H6OOH1X2O2XC3H7O4] =  7.76417085e+01 * (
			 1.05243923e+01 + T * (  3.07940033e-02 
			+ T * ( -1.38088691e-05 + T * (  2.25844932e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3H6OOH1X3O2XC3H7O4] =  7.76417085e+01 * (
			T * (  9.22535194e+00 + T * (  1.63403797e-02
			+ T * ( -4.94110023e-06 + T * (  6.12092480e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.99593555e+04 );
		cp[sC3H6OOH1X3O2XC3H7O4] =  7.76417085e+01 * (
			 9.22535194e+00 + T * (  3.26807594e-02 
			+ T * ( -1.48233007e-05 + T * (  2.44836992e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3H6OOH2X1O2XC3H7O4] =  7.76417085e+01 * (
			T * (  1.05243923e+01 + T * (  1.53970017e-02
			+ T * ( -4.60295637e-06 + T * (  5.64612330e-10
			+ T *  0.00000000e+00 ) ) ) ) -2.27195952e+04 );
		cp[sC3H6OOH2X1O2XC3H7O4] =  7.76417085e+01 * (
			 1.05243923e+01 + T * (  3.07940033e-02 
			+ T * ( -1.38088691e-05 + T * (  2.25844932e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3KET12XC3H6O3] =  9.23015609e+01 * (
			T * (  8.56015723e+00 + T * (  1.40662572e-02
			+ T * ( -4.40872123e-06 + T * (  5.63532810e-10
			+ T *  0.00000000e+00 ) ) ) ) -3.81432115e+04 );
		cp[sC3KET12XC3H6O3] =  9.23015609e+01 * (
			 8.56015723e+00 + T * (  2.81325143e-02 
			+ T * ( -1.32261637e-05 + T * (  2.25413124e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3KET13XC3H6O3] =  9.23015609e+01 * (
			T * (  7.36879167e+00 + T * (  1.44949036e-02
			+ T * ( -4.46309407e-06 + T * (  5.61978412e-10
			+ T *  0.00000000e+00 ) ) ) ) -3.57832882e+04 );
		cp[sC3KET13XC3H6O3] =  9.23015609e+01 * (
			 7.36879167e+00 + T * (  2.89898073e-02 
			+ T * ( -1.33892822e-05 + T * (  2.24791365e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3KET21XC3H6O3] =  9.23015609e+01 * (
			T * (  5.66290215e+00 + T * (  1.60890021e-02
			+ T * ( -5.12006937e-06 + T * (  6.62461047e-10
			+ T *  0.00000000e+00 ) ) ) ) -3.89423019e+04 );
		cp[sC3KET21XC3H6O3] =  9.23015609e+01 * (
			 5.66290215e+00 + T * (  3.21780042e-02 
			+ T * ( -1.53602081e-05 + T * (  2.64984419e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3H52X1J3OOHXC3H7O4] =  7.76417085e+01 * (
			T * (  9.67696819e+00 + T * (  1.68174876e-02
			+ T * ( -5.36589597e-06 + T * (  6.94307162e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.36030398e+04 );
		cp[sC3H52X1J3OOHXC3H7O4] =  7.76417085e+01 * (
			 9.67696819e+00 + T * (  3.36349752e-02 
			+ T * ( -1.60976879e-05 + T * (  2.77722865e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sAXC3H5OOH] =  1.12237641e+02 * (
			T * (  7.07396238e+00 + T * (  1.25148535e-02
			+ T * ( -3.74323233e-06 + T * (  4.65006210e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.05510508e+04 );
		cp[sAXC3H5OOH] =  1.12237641e+02 * (
			 7.07396238e+00 + T * (  2.50297070e-02 
			+ T * ( -1.12296970e-05 + T * (  1.86002484e-09
			+ T *  0.00000000e+00 ) ) ) );
	}
	else {
		h[sN2] =  2.96728765e+02 * (
			T * (  3.29867700e+00 + T * (  7.04120000e-04
			+ T * ( -1.32107400e-06 + T * (  1.41037875e-09
			+ T * -4.88971000e-13 ) ) ) ) -1.02090000e+03 );
		cp[sN2] =  2.96728765e+02 * (
			 3.29867700e+00 + T * (  1.40824000e-03 
			+ T * ( -3.96322200e-06 + T * (  5.64151500e-09
			+ T * -2.44485500e-12 ) ) ) );
		h[sAR] =  2.08129068e+02 * (
			T * (  2.50000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) -7.45375000e+02 );
		cp[sAR] =  2.08129068e+02 * (
			 2.50000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sH] =  8.24835317e+03 * (
			T * (  2.50000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  2.54716300e+04 );
		cp[sH] =  8.24835317e+03 * (
			 2.50000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sO2] =  2.59823125e+02 * (
			T * (  3.21293600e+00 + T * (  5.63743000e-04
			+ T * ( -1.91871667e-07 + T * (  3.28469250e-10
			+ T * -1.75371080e-13 ) ) ) ) -1.00524900e+03 );
		cp[sO2] =  2.59823125e+02 * (
			 3.21293600e+00 + T * (  1.12748600e-03 
			+ T * ( -5.75615000e-07 + T * (  1.31387700e-09
			+ T * -8.76855400e-13 ) ) ) );
		h[sO] =  5.19646250e+02 * (
			T * (  2.94642900e+00 + T * ( -8.19083000e-04
			+ T * (  8.07010667e-07 + T * ( -4.00710750e-10
			+ T *  7.78139200e-14 ) ) ) ) +  2.91476400e+04 );
		cp[sO] =  5.19646250e+02 * (
			 2.94642900e+00 + T * ( -1.63816600e-03 
			+ T * (  2.42103200e-06 + T * ( -1.60284300e-09
			+ T *  3.89069600e-13 ) ) ) );
		h[sOH] =  4.88848777e+02 * (
			T * (  3.41896226e+00 + T * (  1.59627901e-04
			+ T * ( -1.02764239e-07 + T * (  9.11018735e-11
			+ T * -2.00390958e-14 ) ) ) ) +  3.45264448e+03 );
		cp[sOH] =  4.88848777e+02 * (
			 3.41896226e+00 + T * (  3.19255801e-04 
			+ T * ( -3.08292717e-07 + T * (  3.64407494e-10
			+ T * -1.00195479e-13 ) ) ) );
		h[sH2] =  4.12417659e+03 * (
			T * (  3.29812400e+00 + T * (  4.12472100e-04
			+ T * ( -2.71433833e-07 + T * ( -2.36885850e-11
			+ T *  8.26974400e-14 ) ) ) ) -1.01252100e+03 );
		cp[sH2] =  4.12417659e+03 * (
			 3.29812400e+00 + T * (  8.24944200e-04 
			+ T * ( -8.14301500e-07 + T * ( -9.47543400e-11
			+ T *  4.13487200e-13 ) ) ) );
		h[sH2O] =  4.61497558e+02 * (
			T * (  3.38684200e+00 + T * (  1.73749100e-03
			+ T * ( -2.11823200e-06 + T * (  1.74214525e-09
			+ T * -5.01317600e-13 ) ) ) ) -3.02081100e+04 );
		cp[sH2O] =  4.61497558e+02 * (
			 3.38684200e+00 + T * (  3.47498200e-03 
			+ T * ( -6.35469600e-06 + T * (  6.96858100e-09
			+ T * -2.50658800e-12 ) ) ) );
		h[sHO2] =  2.51888633e+02 * (
			T * (  4.30179801e+00 + T * ( -2.37456025e-03
			+ T * (  7.05276303e-06 + T * ( -6.06909735e-09
			+ T *  1.85845025e-12 ) ) ) ) +  2.94808040e+02 );
		cp[sHO2] =  2.51888633e+02 * (
			 4.30179801e+00 + T * ( -4.74912051e-03 
			+ T * (  2.11582891e-05 + T * ( -2.42763894e-08
			+ T *  9.29225124e-12 ) ) ) );
		h[sH2O2] =  2.44424389e+02 * (
			T * (  3.38875400e+00 + T * (  3.28461300e-03
			+ T * ( -4.95004333e-08 + T * ( -1.15645150e-09
			+ T *  4.94303000e-13 ) ) ) ) -1.76631500e+04 );
		cp[sH2O2] =  2.44424389e+02 * (
			 3.38875400e+00 + T * (  6.56922600e-03 
			+ T * ( -1.48501300e-07 + T * ( -4.62580600e-09
			+ T *  2.47151500e-12 ) ) ) );
		h[sCO] =  2.96834702e+02 * (
			T * (  3.26245200e+00 + T * (  7.55970500e-04
			+ T * ( -1.29391833e-06 + T * (  1.39548600e-09
			+ T * -4.94990200e-13 ) ) ) ) -1.43105400e+04 );
		cp[sCO] =  2.96834702e+02 * (
			 3.26245200e+00 + T * (  1.51194100e-03 
			+ T * ( -3.88175500e-06 + T * (  5.58194400e-09
			+ T * -2.47495100e-12 ) ) ) );
		h[sCO2] =  1.88919337e+02 * (
			T * (  2.27572500e+00 + T * (  4.96103600e-03
			+ T * ( -3.46970333e-06 + T * (  1.71667175e-09
			+ T * -4.23456000e-13 ) ) ) ) -4.83731400e+04 );
		cp[sCO2] =  1.88919337e+02 * (
			 2.27572500e+00 + T * (  9.92207200e-03 
			+ T * ( -1.04091100e-05 + T * (  6.86668700e-09
			+ T * -2.11728000e-12 ) ) ) );
		h[sHCO] =  2.86523537e+02 * (
			T * (  2.89833000e+00 + T * (  3.09957350e-03
			+ T * ( -3.20769467e-06 + T * (  2.72456250e-09
			+ T * -9.14977000e-13 ) ) ) ) +  4.15992200e+03 );
		cp[sHCO] =  2.86523537e+02 * (
			 2.89833000e+00 + T * (  6.19914700e-03 
			+ T * ( -9.62308400e-06 + T * (  1.08982500e-08
			+ T * -4.57488500e-12 ) ) ) );
		h[sCH3] =  5.53035786e+02 * (
			T * (  2.43044300e+00 + T * (  5.56205000e-03
			+ T * ( -5.60073333e-06 + T * (  4.05457250e-09
			+ T * -1.17299060e-12 ) ) ) ) +  1.64237800e+04 );
		cp[sCH3] =  5.53035786e+02 * (
			 2.43044300e+00 + T * (  1.11241000e-02 
			+ T * ( -1.68022000e-05 + T * (  1.62182900e-08
			+ T * -5.86495300e-12 ) ) ) );
		h[sCH4] =  5.18285750e+02 * (
			T * (  7.78741500e-01 + T * (  8.73834000e-03
			+ T * ( -9.27803000e-06 + T * (  7.62427000e-09
			+ T * -2.44786200e-12 ) ) ) ) -9.82522900e+03 );
		cp[sCH4] =  5.18285750e+02 * (
			 7.78741500e-01 + T * (  1.74766800e-02 
			+ T * ( -2.78340900e-05 + T * (  3.04970800e-08
			+ T * -1.22393100e-11 ) ) ) );
		h[sCH2O] =  2.76904683e+02 * (
			T * (  1.65273100e+00 + T * (  6.31572000e-03
			+ T * ( -6.29389333e-06 + T * (  5.12507750e-09
			+ T * -1.68264740e-12 ) ) ) ) -1.48654000e+04 );
		cp[sCH2O] =  2.76904683e+02 * (
			 1.65273100e+00 + T * (  1.26314400e-02 
			+ T * ( -1.88816800e-05 + T * (  2.05003100e-08
			+ T * -8.41323700e-12 ) ) ) );
		h[sHOCHO] =  1.80644418e+02 * (
			T * (  1.43548185e+00 + T * (  8.16815080e-03
			+ T * ( -3.54191403e-06 + T * (  8.30332443e-10
			+ T * -8.04352206e-14 ) ) ) ) -4.64616504e+04 );
		cp[sHOCHO] =  1.80644418e+02 * (
			 1.43548185e+00 + T * (  1.63363016e-02 
			+ T * ( -1.06257421e-05 + T * (  3.32132977e-09
			+ T * -4.02176103e-13 ) ) ) );
		h[sCH3O] =  2.67910679e+02 * (
			T * (  2.10620400e+00 + T * (  3.60829750e-03
			+ T * (  1.77949067e-06 + T * ( -1.84440900e-09
			+ T *  4.15122200e-13 ) ) ) ) +  9.78601100e+02 );
		cp[sCH3O] =  2.67910679e+02 * (
			 2.10620400e+00 + T * (  7.21659500e-03 
			+ T * (  5.33847200e-06 + T * ( -7.37763600e-09
			+ T *  2.07561100e-12 ) ) ) );
		h[sCH3OH] =  2.59482554e+02 * (
			T * (  2.66011500e+00 + T * (  3.67075400e-03
			+ T * (  2.39001700e-06 + T * ( -2.19829850e-09
			+ T *  4.78114000e-13 ) ) ) ) -2.53534800e+04 );
		cp[sCH3OH] =  2.59482554e+02 * (
			 2.66011500e+00 + T * (  7.34150800e-03 
			+ T * (  7.17005100e-06 + T * ( -8.79319400e-09
			+ T *  2.39057000e-12 ) ) ) );
		h[sCH3O2] =  1.76772973e+02 * (
			T * (  4.26146906e+00 + T * (  5.04367995e-03
			+ T * ( -1.07168728e-06 + T * (  5.23523168e-11
			+ T *  8.36678206e-15 ) ) ) ) -6.84394259e+02 );
		cp[sCH3O2] =  1.76772973e+02 * (
			 4.26146906e+00 + T * (  1.00873599e-02 
			+ T * ( -3.21506184e-06 + T * (  2.09409267e-10
			+ T *  4.18339103e-14 ) ) ) );
		h[sCH3O2H] =  1.73063986e+02 * (
			T * (  3.23442817e+00 + T * (  9.50648835e-03
			+ T * ( -3.77954290e-06 + T * (  8.50766632e-10
			+ T * -8.23660444e-14 ) ) ) ) -1.77197926e+04 );
		cp[sCH3O2H] =  1.73063986e+02 * (
			 3.23442817e+00 + T * (  1.90129767e-02 
			+ T * ( -1.13386287e-05 + T * (  3.40306653e-09
			+ T * -4.11830222e-13 ) ) ) );
		h[sC2H6] =  2.76517893e+02 * (
			T * ( -2.52854344e-02 + T * (  1.20382377e-02
			+ T * ( -3.72978240e-06 + T * (  5.20852252e-10
			+ T * -1.05973723e-14 ) ) ) ) -1.12345534e+04 );
		cp[sC2H6] =  2.76517893e+02 * (
			-2.52854344e-02 + T * (  2.40764754e-02 
			+ T * ( -1.11893472e-05 + T * (  2.08340901e-09
			+ T * -5.29868616e-14 ) ) ) );
		h[sC2H5] =  2.86109429e+02 * (
			T * (  2.69070200e+00 + T * (  4.35956650e-03
			+ T * (  1.47327967e-06 + T * (  2.33467575e-10
			+ T * -7.85554600e-13 ) ) ) ) +  1.28704000e+04 );
		cp[sC2H5] =  2.86109429e+02 * (
			 2.69070200e+00 + T * (  8.71913300e-03 
			+ T * (  4.41983900e-06 + T * (  9.33870300e-10
			+ T * -3.92777300e-12 ) ) ) );
		h[sC2H4] =  2.96390275e+02 * (
			T * ( -8.61488000e-01 + T * (  1.39808150e-02
			+ T * ( -1.12955900e-05 + T * (  6.96288000e-09
			+ T * -1.94757580e-12 ) ) ) ) +  5.57304600e+03 );
		cp[sC2H4] =  2.96390275e+02 * (
			-8.61488000e-01 + T * (  2.79616300e-02 
			+ T * ( -3.38867700e-05 + T * (  2.78515200e-08
			+ T * -9.73787900e-12 ) ) ) );
		h[sCH3CHO] =  1.88739217e+02 * (
			T * (  2.50569500e+00 + T * (  6.68495500e-03
			+ T * (  1.55731767e-06 + T * ( -2.82035000e-09
			+ T *  8.52713200e-13 ) ) ) ) -2.12458900e+04 );
		cp[sCH3CHO] =  1.88739217e+02 * (
			 2.50569500e+00 + T * (  1.33699100e-02 
			+ T * (  4.67195300e-06 + T * ( -1.12814000e-08
			+ T *  4.26356600e-12 ) ) ) );
		h[sCH3CO] =  1.93159093e+02 * (
			T * (  3.12527800e+00 + T * (  4.88911000e-03
			+ T * (  1.50714933e-06 + T * ( -2.25236550e-09
			+ T *  6.38743600e-13 ) ) ) ) -4.10850800e+03 );
		cp[sCH3CO] =  1.93159093e+02 * (
			 3.12527800e+00 + T * (  9.77822000e-03 
			+ T * (  4.52144800e-06 + T * ( -9.00946200e-09
			+ T *  3.19371800e-12 ) ) ) );
		h[sCH3CO3] =  1.10792868e+02 * (
			T * (  3.60373432e+00 + T * (  1.35040171e-02
			+ T * ( -6.94311460e-06 + T * (  2.12635276e-09
			+ T * -2.87692220e-13 ) ) ) ) -2.34205171e+04 );
		cp[sCH3CO3] =  1.10792868e+02 * (
			 3.60373432e+00 + T * (  2.70080341e-02 
			+ T * ( -2.08293438e-05 + T * (  8.50541104e-09
			+ T * -1.43846110e-12 ) ) ) );
		h[sCH3CO3H] =  1.09324410e+02 * (
			T * (  2.24135876e+00 + T * (  1.68981757e-02
			+ T * ( -8.46291607e-06 + T * (  2.41895897e-09
			+ T * -2.98532314e-13 ) ) ) ) -4.24677831e+04 );
		cp[sCH3CO3H] =  1.09324410e+02 * (
			 2.24135876e+00 + T * (  3.37963514e-02 
			+ T * ( -2.53887482e-05 + T * (  9.67583587e-09
			+ T * -1.49266157e-12 ) ) ) );
		h[sCH2CHO] =  1.93159093e+02 * (
			T * (  3.40906200e+00 + T * (  5.36928500e-03
			+ T * (  6.30497333e-07 + T * ( -1.78964575e-09
			+ T *  5.73477000e-13 ) ) ) ) +  1.52147700e+03 );
		cp[sCH2CHO] =  1.93159093e+02 * (
			 3.40906200e+00 + T * (  1.07385700e-02 
			+ T * (  1.89149200e-06 + T * ( -7.15858300e-09
			+ T *  2.86738500e-12 ) ) ) );
		h[sCH3CO2] =  1.40816002e+02 * (
			T * (  1.37440768e+00 + T * (  1.24557802e-02
			+ T * ( -5.81029647e-06 + T * (  1.56199877e-09
			+ T * -1.81903367e-13 ) ) ) ) -2.72330150e+04 );
		cp[sCH3CO2] =  1.40816002e+02 * (
			 1.37440768e+00 + T * (  2.49115604e-02 
			+ T * ( -1.74308894e-05 + T * (  6.24799508e-09
			+ T * -9.09516835e-13 ) ) ) );
		h[sC2H3] =  3.07437509e+02 * (
			T * (  2.45927600e+00 + T * (  3.68573800e-03
			+ T * (  7.03291000e-07 + T * ( -3.30410500e-10
			+ T * -2.36956800e-13 ) ) ) ) +  3.33522500e+04 );
		cp[sC2H3] =  3.07437509e+02 * (
			 2.45927600e+00 + T * (  7.37147600e-03 
			+ T * (  2.10987300e-06 + T * ( -1.32164200e-09
			+ T * -1.18478400e-12 ) ) ) );
		h[sC2H3CHO] =  1.48306161e+02 * (
			T * (  2.92355162e-01 + T * (  1.77160709e-02
			+ T * ( -9.83121080e-06 + T * (  3.20250310e-09
			+ T * -4.52288216e-13 ) ) ) ) -1.16521584e+04 );
		cp[sC2H3CHO] =  1.48306161e+02 * (
			 2.92355162e-01 + T * (  3.54321417e-02 
			+ T * ( -2.94936324e-05 + T * (  1.28100124e-08
			+ T * -2.26144108e-12 ) ) ) );
		h[sC2H3CO] =  1.51021542e+02 * (
			T * (  1.36242013e+00 + T * (  1.57636986e-02
			+ T * ( -1.00072978e-05 + T * (  3.70417780e-09
			+ T * -5.75943060e-13 ) ) ) ) +  4.25770215e+03 );
		cp[sC2H3CO] =  1.51021542e+02 * (
			 1.36242013e+00 + T * (  3.15273972e-02 
			+ T * ( -3.00218935e-05 + T * (  1.48167112e-08
			+ T * -2.87971530e-12 ) ) ) );
		h[sC3H8] =  1.88559441e+02 * (
			T * ( -2.93546760e-01 + T * (  1.81935127e-02
			+ T * ( -6.50165410e-06 + T * (  1.23975106e-09
			+ T * -9.11839576e-14 ) ) ) ) -1.41148884e+04 );
		cp[sC3H8] =  1.88559441e+02 * (
			-2.93546760e-01 + T * (  3.63870253e-02 
			+ T * ( -1.95049623e-05 + T * (  4.95900425e-09
			+ T * -4.55919788e-13 ) ) ) );
		h[sNXC3H7] =  1.92970803e+02 * (
			T * (  4.89680103e-01 + T * (  1.58872692e-02
			+ T * ( -5.64735837e-06 + T * (  1.06805725e-09
			+ T * -7.76212142e-14 ) ) ) ) +  1.04955773e+04 );
		cp[sNXC3H7] =  1.92970803e+02 * (
			 4.89680103e-01 + T * (  3.17745383e-02 
			+ T * ( -1.69420751e-05 + T * (  4.27222899e-09
			+ T * -3.88106071e-13 ) ) ) );
		h[sIXC3H7] =  1.92970803e+02 * (
			T * (  1.56896579e+00 + T * (  1.20745612e-02
			+ T * ( -1.60491576e-06 + T * ( -8.13848015e-10
			+ T *  2.50142984e-13 ) ) ) ) +  9.08600357e+03 );
		cp[sIXC3H7] =  1.92970803e+02 * (
			 1.56896579e+00 + T * (  2.41491224e-02 
			+ T * ( -4.81474729e-06 + T * ( -3.25539206e-09
			+ T *  1.25071492e-12 ) ) ) );
		h[sC3H5XAXC3H5] =  2.02443146e+02 * (
			T * ( -5.29131958e-01 + T * (  1.67279550e-02
			+ T * ( -8.44670090e-06 + T * (  2.57164385e-09
			+ T * -3.46516680e-13 ) ) ) ) +  1.93834226e+04 );
		cp[sC3H5XAXC3H5] =  2.02443146e+02 * (
			-5.29131958e-01 + T * (  3.34559100e-02 
			+ T * ( -2.53401027e-05 + T * (  1.02865754e-08
			+ T * -1.73258340e-12 ) ) ) );
		h[sC3H6] =  1.97593517e+02 * (
			T * (  3.94615444e-01 + T * (  1.44553831e-02
			+ T * ( -5.16289360e-06 + T * (  9.72035522e-10
			+ T * -6.75780704e-14 ) ) ) ) +  1.06688164e+03 );
		cp[sC3H6] =  1.97593517e+02 * (
			 3.94615444e-01 + T * (  2.89107662e-02 
			+ T * ( -1.54886808e-05 + T * (  3.88814209e-09
			+ T * -3.37890352e-13 ) ) ) );
		h[sIXC3H7O2] =  1.10730895e+02 * (
			T * (  1.49941639e+00 + T * (  2.21540602e-02
			+ T * ( -1.07471485e-05 + T * (  3.24217840e-09
			+ T * -4.46741138e-13 ) ) ) ) -1.02587980e+04 );
		cp[sIXC3H7O2] =  1.10730895e+02 * (
			 1.49941639e+00 + T * (  4.43081205e-02 
			+ T * ( -3.22414456e-05 + T * (  1.29687136e-08
			+ T * -2.23370569e-12 ) ) ) );
		h[sIXC3H7O2H] =  1.09264068e+02 * (
			T * (  5.19265570e-01 + T * (  2.66055614e-02
			+ T * ( -1.35052297e-05 + T * (  4.08366782e-09
			+ T * -5.47502466e-13 ) ) ) ) -2.71048486e+04 );
		cp[sIXC3H7O2H] =  1.09264068e+02 * (
			 5.19265570e-01 + T * (  5.32111228e-02 
			+ T * ( -4.05156892e-05 + T * (  1.63346713e-08
			+ T * -2.73751233e-12 ) ) ) );
		h[sC3H5XTXC3H5] =  2.02443146e+02 * (
			T * (  2.17916644e+00 + T * (  1.01913311e-02
			+ T * ( -2.63804611e-06 + T * (  1.19226547e-10
			+ T *  5.40797072e-14 ) ) ) ) +  2.94895338e+04 );
		cp[sC3H5XTXC3H5] =  2.02443146e+02 * (
			 2.17916644e+00 + T * (  2.03826623e-02 
			+ T * ( -7.91413834e-06 + T * (  4.76906187e-10
			+ T *  2.70398536e-13 ) ) ) );
		h[sC3H5O] =  1.45686701e+02 * (
			T * (  1.19822582e+00 + T * (  1.52789918e-02
			+ T * ( -6.02100920e-06 + T * (  1.21537508e-09
			+ T * -8.39709124e-14 ) ) ) ) +  9.58217784e+03 );
		cp[sC3H5O] =  1.45686701e+02 * (
			 1.19822582e+00 + T * (  3.05579837e-02 
			+ T * ( -1.80630276e-05 + T * (  4.86150033e-09
			+ T * -4.19854562e-13 ) ) ) );
		h[sCH3COCH2] =  1.45686701e+02 * (
			T * (  1.80339187e+00 + T * (  1.50703542e-02
			+ T * ( -6.45018507e-06 + T * (  1.59549759e-09
			+ T * -1.73220636e-13 ) ) ) ) -5.37233261e+03 );
		cp[sCH3COCH2] =  1.45686701e+02 * (
			 1.80339187e+00 + T * (  3.01407085e-02 
			+ T * ( -1.93505552e-05 + T * (  6.38199034e-09
			+ T * -8.66103180e-13 ) ) ) );
		h[sIXC3H7O] =  1.40715906e+02 * (
			T * ( -1.00369576e-01 + T * (  2.15559464e-02
			+ T * ( -1.03543939e-05 + T * (  3.01670910e-09
			+ T * -3.94429934e-13 ) ) ) ) -8.44902998e+03 );
		cp[sIXC3H7O] =  1.40715906e+02 * (
			-1.00369576e-01 + T * (  4.31118928e-02 
			+ T * ( -3.10631817e-05 + T * (  1.20668364e-08
			+ T * -1.97214967e-12 ) ) ) );
		h[sNXC3H7O2] =  1.10730895e+02 * (
			T * (  2.10731492e+00 + T * (  1.98082493e-02
			+ T * ( -8.31638663e-06 + T * (  2.14862575e-09
			+ T * -2.62480660e-13 ) ) ) ) -7.93745567e+03 );
		cp[sNXC3H7O2] =  1.10730895e+02 * (
			 2.10731492e+00 + T * (  3.96164986e-02 
			+ T * ( -2.49491599e-05 + T * (  8.59450300e-09
			+ T * -1.31240330e-12 ) ) ) );
		h[sC3H6OOH1X2XC3H7O2] =  1.10730895e+02 * (
			T * (  2.87774562e+00 + T * (  1.87083500e-02
			+ T * ( -7.86860210e-06 + T * (  1.94982715e-09
			+ T * -2.12085124e-13 ) ) ) ) -7.82368119e+02 );
		cp[sC3H6OOH1X2XC3H7O2] =  1.10730895e+02 * (
			 2.87774562e+00 + T * (  3.74166999e-02 
			+ T * ( -2.36058063e-05 + T * (  7.79930860e-09
			+ T * -1.06042562e-12 ) ) ) );
		h[sC3H6OOH1X3XC3H7O2] =  1.10730895e+02 * (
			T * (  1.88331465e+00 + T * (  2.20078026e-02
			+ T * ( -1.02619487e-05 + T * (  2.83037155e-09
			+ T * -3.50646368e-13 ) ) ) ) -1.68779168e+02 );
		cp[sC3H6OOH1X3XC3H7O2] =  1.10730895e+02 * (
			 1.88331465e+00 + T * (  4.40156051e-02 
			+ T * ( -3.07858462e-05 + T * (  1.13214862e-08
			+ T * -1.75323184e-12 ) ) ) );
		h[sC3H6OOH2X1XC3H7O2] =  1.10730895e+02 * (
			T * (  2.09193950e+00 + T * (  2.34610197e-02
			+ T * ( -1.30093610e-05 + T * (  4.30953633e-09
			+ T * -6.15937958e-13 ) ) ) ) -1.89377918e+03 );
		cp[sC3H6OOH2X1XC3H7O2] =  1.10730895e+02 * (
			 2.09193950e+00 + T * (  4.69220394e-02 
			+ T * ( -3.90280831e-05 + T * (  1.72381453e-08
			+ T * -3.07968979e-12 ) ) ) );
		h[sC3H6OOH1X2O2XC3H7O4] =  7.76417085e+01 * (
			T * (  2.65196584e+00 + T * (  2.87319075e-02
			+ T * ( -1.57396956e-05 + T * (  5.13978893e-09
			+ T * -7.37574774e-13 ) ) ) ) -2.08829371e+04 );
		cp[sC3H6OOH1X2O2XC3H7O4] =  7.76417085e+01 * (
			 2.65196584e+00 + T * (  5.74638149e-02 
			+ T * ( -4.72190867e-05 + T * (  2.05591557e-08
			+ T * -3.68787387e-12 ) ) ) );
		h[sC3H6OOH1X3O2XC3H7O4] =  7.76417085e+01 * (
			T * (  3.14864588e+00 + T * (  2.66771285e-02
			+ T * ( -1.36110204e-05 + T * (  4.18855550e-09
			+ T * -5.78576652e-13 ) ) ) ) -1.85473645e+04 );
		cp[sC3H6OOH1X3O2XC3H7O4] =  7.76417085e+01 * (
			 3.14864588e+00 + T * (  5.33542571e-02 
			+ T * ( -4.08330611e-05 + T * (  1.67542220e-08
			+ T * -2.89288326e-12 ) ) ) );
		h[sC3H6OOH2X1O2XC3H7O4] =  7.76417085e+01 * (
			T * (  2.65196584e+00 + T * (  2.87319075e-02
			+ T * ( -1.57396956e-05 + T * (  5.13978893e-09
			+ T * -7.37574774e-13 ) ) ) ) -2.08829371e+04 );
		cp[sC3H6OOH2X1O2XC3H7O4] =  7.76417085e+01 * (
			 2.65196584e+00 + T * (  5.74638149e-02 
			+ T * ( -4.72190867e-05 + T * (  2.05591557e-08
			+ T * -3.68787387e-12 ) ) ) );
		h[sC3KET12XC3H6O3] =  9.23015609e+01 * (
			T * (  1.10507238e+00 + T * (  2.63698353e-02
			+ T * ( -1.43935258e-05 + T * (  4.54045730e-09
			+ T * -6.21071136e-13 ) ) ) ) -3.63627536e+04 );
		cp[sC3KET12XC3H6O3] =  9.23015609e+01 * (
			 1.10507238e+00 + T * (  5.27396706e-02 
			+ T * ( -4.31805774e-05 + T * (  1.81618292e-08
			+ T * -3.10535568e-12 ) ) ) );
		h[sC3KET13XC3H6O3] =  9.23015609e+01 * (
			T * (  3.55241022e+00 + T * (  2.09360135e-02
			+ T * ( -9.81834567e-06 + T * (  2.74957250e-09
			+ T * -3.50091954e-13 ) ) ) ) -3.48902671e+04 );
		cp[sC3KET13XC3H6O3] =  9.23015609e+01 * (
			 3.55241022e+00 + T * (  4.18720270e-02 
			+ T * ( -2.94550370e-05 + T * (  1.09982900e-08
			+ T * -1.75045977e-12 ) ) ) );
		h[sC3KET21XC3H6O3] =  9.23015609e+01 * (
			T * (  4.55686367e+00 + T * (  1.78538418e-02
			+ T * ( -6.49040180e-06 + T * (  1.17673858e-09
			+ T * -7.39507614e-14 ) ) ) ) -3.86710975e+04 );
		cp[sC3KET21XC3H6O3] =  9.23015609e+01 * (
			 4.55686367e+00 + T * (  3.57076837e-02 
			+ T * ( -1.94712054e-05 + T * (  4.70695431e-09
			+ T * -3.69753807e-13 ) ) ) );
		h[sC3H52X1J3OOHXC3H7O4] =  7.76417085e+01 * (
			T * (  4.12253742e+00 + T * (  2.59776805e-02
			+ T * ( -1.27911242e-05 + T * (  3.64629093e-09
			+ T * -4.59641072e-13 ) ) ) ) -1.22759164e+04 );
		cp[sC3H52X1J3OOHXC3H7O4] =  7.76417085e+01 * (
			 4.12253742e+00 + T * (  5.19553611e-02 
			+ T * ( -3.83733727e-05 + T * (  1.45851637e-08
			+ T * -2.29820536e-12 ) ) ) );
		h[sAXC3H5OOH] =  1.12237641e+02 * (
			T * (  2.43934647e+00 + T * (  2.01035319e-02
			+ T * ( -9.84408930e-06 + T * (  2.86789000e-09
			+ T * -3.70341022e-13 ) ) ) ) -9.43679906e+03 );
		cp[sAXC3H5OOH] =  1.12237641e+02 * (
			 2.43934647e+00 + T * (  4.02070638e-02 
			+ T * ( -2.95322679e-05 + T * (  1.14715600e-08
			+ T * -1.85170511e-12 ) ) ) );
	}
}


int    Get_Num_Species_in_Chem_Mech_Ngi()
{
  return 56;
}
int    Get_Num_Species_in_Chem_Mech_Ni()
{
  return 59;
}
int    Get_Num_Elementry_RXN_in_Chem_Mech_Ner()
{
  return 357;
}

