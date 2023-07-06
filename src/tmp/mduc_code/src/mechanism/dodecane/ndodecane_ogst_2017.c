#include "ndodecane_ogst_2017.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
static double GetLindRateCoeff( double temp, double pressure
			, double k0, double kInf, double fc, double conc );

static double GetPlogRateCoeff( double temp, double pressure, double lgt, double rt_inv
			, double *PlogP, double *PlogA, double *PlogB, double *PlogE, int np );

void ComputeProductionRates( double *cdot, double *w, double *k
			, double *c, double *M, double temp, double pressure );
double MAX_C(double X1, double X2);
double MIN_C(double X1, double X2);
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

	int	nSpec = 155;
	int	nSpecIn = 155;
        int     i;
	double	kTroe0, kTroeInf, fcTroe;
	double	RGAS = 8314.34;
	double	lgt = log( temp );
	double	rt_inv = 1.0/(RGAS * temp);
	double PlogA[14],PlogB[14],PlogE[14],PlogP[14];
	int	np;

        if (temp<300.0) {

	   for (i=0; i<nSpec;i++) cdot[i] = 0.0;

	   return;
        }

	M[mM15] = c[sN2] + c[sAR] + c[sO] + 2 * c[sH2]
		 + c[sH] + c[sOH] + 14 * c[sH2O] + 0.78 * c[sO2]
		 + c[sHO2] + c[sH2O2] + c[sCH] + c[sTXCH2]
		 + c[sCH2O] + c[sHCO] + c[sCH3] + 3.8 * c[sCO2]
		 + 1.9 * c[sCO] + c[sSXCH2] + c[sC2H4] + c[sCH3O]
		 + c[sCH3O2] + c[sC2H5] + c[sCH2CO] + c[sO2CHO]
		 + c[sHOCHO] + c[sOCHO] + c[sC2H2] + c[sC2H3]
		 + c[sHCCO] + c[sCH2CHO] + c[sC2H5O] + c[sAXC3H5]
		 + c[sC3H6] + c[sCH3COCH2] + c[sNXC3H7] + c[sPXC4H8]
		 + c[sPXC4H9] + c[sC5H11] + c[sC6H12X1XC6H12] + c[sC6H13]
		 + c[sC7H14] + c[sC7H15] + c[sC8H17] + c[sC9H19]
		 + c[sNXC12H26] + c[sC12H25] + c[sC12H24] + c[sC12H25O2]
		 + c[sC12OOHXT6XC12H25O2] + c[sC12OOHO2XT6XC12H25O4] + c[sC12KETXT6XC12H24O3] + c[sNXC8H17CHO]
		 + c[sC7H15CO] + c[sC3H7CHO] + c[sNXC7H15CHO] + c[sC7H15O2X1XC7H15O2]
		 + c[sC7H14OOH1X3XC7H15O2] + c[sC5H4CH2] + c[sA1XC6H6] + c[sNXC4H5]
		 + c[sIXC4H5] + c[sA1XXC6H5] + c[sNXC4H3] + c[sC3H3]
		 + c[sA1C2H2XC8H7] + c[sA1C2H3XC8H8] + c[sA1C2HXC8H6] + c[sA1C2HYXC8H5]
		 + c[sA1C2H3YXC8H7] + c[sA2XXC10H7] + c[sA2XC10H8] + c[sC4H4]
		 + c[sA2YXC10H7] + c[sA2C2H2AXC12H9] + c[sA2C2H2BXC12H9] + c[sA2C2HAXC12H8]
		 + c[sA2C2HBXC12H8] + c[sA2C2HAYXC12H7] + c[sA2C2HBYXC12H7] + c[sA2R5XC12H8]
		 + c[sA2R5XXC12H7] + c[sA2R5C2H2XC14H9] + c[sA2R5C2HXC14H8] + c[sA2R5C2HYXC14H7]
		 + c[sP2XC12H10] + c[sP2XXC12H9] + c[sA3XXC14H9] + c[sA3XC14H10]
		 + c[sA3YXC14H9] + c[sA3R5XXC16H9] + c[sA3R5XC16H10] + c[sA4XC16H10]
		 + c[sA4XXC16H9] + c[sA4R5XC18H10] + c[sFLTNXC16H10] + c[sC5H6]
		 + c[sC5H5] + c[sSXC3H5] + c[sCH4] + c[sC4H6]
		 + c[sTXC5H5O] + c[sC5H4O] + c[sSXC5H5O] + c[sC9H8]
		 + c[sC9H7] + c[sA1CH2XC7H7] + c[sC9H6O] + c[sOXC6H4]
		 + c[sA1CH3XC7H8] + c[sA1OHXC6H6O] + c[sHOA1CH3XC7H8O] + c[sOA1CH3XC7H7O]
		 + c[sA1CH2OXC7H7O] + c[sA1CH2OHXC7H8O] + c[sA1CHOXC7H6O] + c[sA1OXC6H5O]
		 + c[sA1CH3YXC7H7] + c[sPXC3H4] + c[sA1C2H4XC8H9] + c[sA1C2H5XC8H10]
		 + c[sC8H9O2] + c[sC8H8OOH] + c[sOC8H7OOH] + c[sH2C2]
		 + c[sA1CH3CH3XC8H10] + c[sA1CH3CH2XC8H9] + c[sA1CH3CHOXC8H8O] + c[sA2CH3XC11H10]
		 + c[sA1CHOCH2XC8H7O] + c[sA1CHOCHOXC8H6O2] + c[sA2OHXC10H8O] + c[sA2CH2XC11H9]
		 + c[sA2CH2OXC11H9O] + c[sA2CHOXC11H8O] + c[sA2OXC10H7O] + c[sC4H2]
		 + c[sOC6H4O] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sC]
		 + c[sHCNN] + c[sHCNO] + c[sNH3];

	M[mM16] = 1.5 * c[sN2] + c[sAR] + c[sO] + 3.7 * c[sH2]
		 + c[sH] + c[sOH] + 7.5 * c[sH2O] + 1.2 * c[sO2]
		 + c[sHO2] + 7.7 * c[sH2O2] + c[sCH] + c[sTXCH2]
		 + c[sCH2O] + c[sHCO] + c[sCH3] + 1.6 * c[sCO2]
		 + 2.8 * c[sCO] + c[sSXCH2] + c[sC2H4] + c[sCH3O]
		 + c[sCH3O2] + c[sC2H5] + c[sCH2CO] + c[sO2CHO]
		 + c[sHOCHO] + c[sOCHO] + c[sC2H2] + c[sC2H3]
		 + c[sHCCO] + c[sCH2CHO] + c[sC2H5O] + c[sAXC3H5]
		 + c[sC3H6] + c[sCH3COCH2] + c[sNXC3H7] + c[sPXC4H8]
		 + c[sPXC4H9] + c[sC5H11] + c[sC6H12X1XC6H12] + c[sC6H13]
		 + c[sC7H14] + c[sC7H15] + c[sC8H17] + c[sC9H19]
		 + c[sNXC12H26] + c[sC12H25] + c[sC12H24] + c[sC12H25O2]
		 + c[sC12OOHXT6XC12H25O2] + c[sC12OOHO2XT6XC12H25O4] + c[sC12KETXT6XC12H24O3] + c[sNXC8H17CHO]
		 + c[sC7H15CO] + c[sC3H7CHO] + c[sNXC7H15CHO] + c[sC7H15O2X1XC7H15O2]
		 + c[sC7H14OOH1X3XC7H15O2] + c[sC5H4CH2] + c[sA1XC6H6] + c[sNXC4H5]
		 + c[sIXC4H5] + c[sA1XXC6H5] + c[sNXC4H3] + c[sC3H3]
		 + c[sA1C2H2XC8H7] + c[sA1C2H3XC8H8] + c[sA1C2HXC8H6] + c[sA1C2HYXC8H5]
		 + c[sA1C2H3YXC8H7] + c[sA2XXC10H7] + c[sA2XC10H8] + c[sC4H4]
		 + c[sA2YXC10H7] + c[sA2C2H2AXC12H9] + c[sA2C2H2BXC12H9] + c[sA2C2HAXC12H8]
		 + c[sA2C2HBXC12H8] + c[sA2C2HAYXC12H7] + c[sA2C2HBYXC12H7] + c[sA2R5XC12H8]
		 + c[sA2R5XXC12H7] + c[sA2R5C2H2XC14H9] + c[sA2R5C2HXC14H8] + c[sA2R5C2HYXC14H7]
		 + c[sP2XC12H10] + c[sP2XXC12H9] + c[sA3XXC14H9] + c[sA3XC14H10]
		 + c[sA3YXC14H9] + c[sA3R5XXC16H9] + c[sA3R5XC16H10] + c[sA4XC16H10]
		 + c[sA4XXC16H9] + c[sA4R5XC18H10] + c[sFLTNXC16H10] + c[sC5H6]
		 + c[sC5H5] + c[sSXC3H5] + c[sCH4] + c[sC4H6]
		 + c[sTXC5H5O] + c[sC5H4O] + c[sSXC5H5O] + c[sC9H8]
		 + c[sC9H7] + c[sA1CH2XC7H7] + c[sC9H6O] + c[sOXC6H4]
		 + c[sA1CH3XC7H8] + c[sA1OHXC6H6O] + c[sHOA1CH3XC7H8O] + c[sOA1CH3XC7H7O]
		 + c[sA1CH2OXC7H7O] + c[sA1CH2OHXC7H8O] + c[sA1CHOXC7H6O] + c[sA1OXC6H5O]
		 + c[sA1CH3YXC7H7] + c[sPXC3H4] + c[sA1C2H4XC8H9] + c[sA1C2H5XC8H10]
		 + c[sC8H9O2] + c[sC8H8OOH] + c[sOC8H7OOH] + c[sH2C2]
		 + c[sA1CH3CH3XC8H10] + c[sA1CH3CH2XC8H9] + c[sA1CH3CHOXC8H8O] + c[sA2CH3XC11H10]
		 + c[sA1CHOCH2XC8H7O] + c[sA1CHOCHOXC8H6O2] + c[sA2OHXC10H8O] + c[sA2CH2XC11H9]
		 + c[sA2CH2OXC11H9O] + c[sA2CHOXC11H8O] + c[sA2OXC10H7O] + c[sC4H2]
		 + c[sOC6H4O] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sC]
		 + c[sHCNN] + c[sHCNO] + c[sNH3];

	M[mM3] = c[sN2] + c[sAR] + c[sO] + 2 * c[sH2]
		 + c[sH] + c[sOH] + 12 * c[sH2O] + c[sO2]
		 + c[sHO2] + c[sH2O2] + c[sCH] + c[sTXCH2]
		 + c[sCH2O] + c[sHCO] + c[sCH3] + 3.6 * c[sCO2]
		 + 1.75 * c[sCO] + c[sSXCH2] + c[sC2H4] + c[sCH3O]
		 + c[sCH3O2] + c[sC2H5] + c[sCH2CO] + c[sO2CHO]
		 + c[sHOCHO] + c[sOCHO] + c[sC2H2] + c[sC2H3]
		 + c[sHCCO] + c[sCH2CHO] + c[sC2H5O] + c[sAXC3H5]
		 + c[sC3H6] + c[sCH3COCH2] + c[sNXC3H7] + c[sPXC4H8]
		 + c[sPXC4H9] + c[sC5H11] + c[sC6H12X1XC6H12] + c[sC6H13]
		 + c[sC7H14] + c[sC7H15] + c[sC8H17] + c[sC9H19]
		 + c[sNXC12H26] + c[sC12H25] + c[sC12H24] + c[sC12H25O2]
		 + c[sC12OOHXT6XC12H25O2] + c[sC12OOHO2XT6XC12H25O4] + c[sC12KETXT6XC12H24O3] + c[sNXC8H17CHO]
		 + c[sC7H15CO] + c[sC3H7CHO] + c[sNXC7H15CHO] + c[sC7H15O2X1XC7H15O2]
		 + c[sC7H14OOH1X3XC7H15O2] + c[sC5H4CH2] + c[sA1XC6H6] + c[sNXC4H5]
		 + c[sIXC4H5] + c[sA1XXC6H5] + c[sNXC4H3] + c[sC3H3]
		 + c[sA1C2H2XC8H7] + c[sA1C2H3XC8H8] + c[sA1C2HXC8H6] + c[sA1C2HYXC8H5]
		 + c[sA1C2H3YXC8H7] + c[sA2XXC10H7] + c[sA2XC10H8] + c[sC4H4]
		 + c[sA2YXC10H7] + c[sA2C2H2AXC12H9] + c[sA2C2H2BXC12H9] + c[sA2C2HAXC12H8]
		 + c[sA2C2HBXC12H8] + c[sA2C2HAYXC12H7] + c[sA2C2HBYXC12H7] + c[sA2R5XC12H8]
		 + c[sA2R5XXC12H7] + c[sA2R5C2H2XC14H9] + c[sA2R5C2HXC14H8] + c[sA2R5C2HYXC14H7]
		 + c[sP2XC12H10] + c[sP2XXC12H9] + c[sA3XXC14H9] + c[sA3XC14H10]
		 + c[sA3YXC14H9] + c[sA3R5XXC16H9] + c[sA3R5XC16H10] + c[sA4XC16H10]
		 + c[sA4XXC16H9] + c[sA4R5XC18H10] + c[sFLTNXC16H10] + c[sC5H6]
		 + c[sC5H5] + c[sSXC3H5] + c[sCH4] + c[sC4H6]
		 + c[sTXC5H5O] + c[sC5H4O] + c[sSXC5H5O] + c[sC9H8]
		 + c[sC9H7] + c[sA1CH2XC7H7] + c[sC9H6O] + c[sOXC6H4]
		 + c[sA1CH3XC7H8] + c[sA1OHXC6H6O] + c[sHOA1CH3XC7H8O] + c[sOA1CH3XC7H7O]
		 + c[sA1CH2OXC7H7O] + c[sA1CH2OHXC7H8O] + c[sA1CHOXC7H6O] + c[sA1OXC6H5O]
		 + c[sA1CH3YXC7H7] + c[sPXC3H4] + c[sA1C2H4XC8H9] + c[sA1C2H5XC8H10]
		 + c[sC8H9O2] + c[sC8H8OOH] + c[sOC8H7OOH] + c[sH2C2]
		 + c[sA1CH3CH3XC8H10] + c[sA1CH3CH2XC8H9] + c[sA1CH3CHOXC8H8O] + c[sA2CH3XC11H10]
		 + c[sA1CHOCH2XC8H7O] + c[sA1CHOCHOXC8H6O2] + c[sA2OHXC10H8O] + c[sA2CH2XC11H9]
		 + c[sA2CH2OXC11H9O] + c[sA2CHOXC11H8O] + c[sA2OXC10H7O] + c[sC4H2]
		 + c[sOC6H4O] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sC]
		 + c[sHCNN] + c[sHCNO] + c[sNH3];

	M[mM0] = c[sN2] + c[sAR] + c[sO] + c[sH2]
		 + c[sH] + c[sOH] + c[sH2O] + c[sO2]
		 + c[sHO2] + c[sH2O2] + c[sCH] + c[sTXCH2]
		 + c[sCH2O] + c[sHCO] + c[sCH3] + c[sCO2]
		 + c[sCO] + c[sSXCH2] + c[sC2H4] + c[sCH3O]
		 + c[sCH3O2] + c[sC2H5] + c[sCH2CO] + c[sO2CHO]
		 + c[sHOCHO] + c[sOCHO] + c[sC2H2] + c[sC2H3]
		 + c[sHCCO] + c[sCH2CHO] + c[sC2H5O] + c[sAXC3H5]
		 + c[sC3H6] + c[sCH3COCH2] + c[sNXC3H7] + c[sPXC4H8]
		 + c[sPXC4H9] + c[sC5H11] + c[sC6H12X1XC6H12] + c[sC6H13]
		 + c[sC7H14] + c[sC7H15] + c[sC8H17] + c[sC9H19]
		 + c[sNXC12H26] + c[sC12H25] + c[sC12H24] + c[sC12H25O2]
		 + c[sC12OOHXT6XC12H25O2] + c[sC12OOHO2XT6XC12H25O4] + c[sC12KETXT6XC12H24O3] + c[sNXC8H17CHO]
		 + c[sC7H15CO] + c[sC3H7CHO] + c[sNXC7H15CHO] + c[sC7H15O2X1XC7H15O2]
		 + c[sC7H14OOH1X3XC7H15O2] + c[sC5H4CH2] + c[sA1XC6H6] + c[sNXC4H5]
		 + c[sIXC4H5] + c[sA1XXC6H5] + c[sNXC4H3] + c[sC3H3]
		 + c[sA1C2H2XC8H7] + c[sA1C2H3XC8H8] + c[sA1C2HXC8H6] + c[sA1C2HYXC8H5]
		 + c[sA1C2H3YXC8H7] + c[sA2XXC10H7] + c[sA2XC10H8] + c[sC4H4]
		 + c[sA2YXC10H7] + c[sA2C2H2AXC12H9] + c[sA2C2H2BXC12H9] + c[sA2C2HAXC12H8]
		 + c[sA2C2HBXC12H8] + c[sA2C2HAYXC12H7] + c[sA2C2HBYXC12H7] + c[sA2R5XC12H8]
		 + c[sA2R5XXC12H7] + c[sA2R5C2H2XC14H9] + c[sA2R5C2HXC14H8] + c[sA2R5C2HYXC14H7]
		 + c[sP2XC12H10] + c[sP2XXC12H9] + c[sA3XXC14H9] + c[sA3XC14H10]
		 + c[sA3YXC14H9] + c[sA3R5XXC16H9] + c[sA3R5XC16H10] + c[sA4XC16H10]
		 + c[sA4XXC16H9] + c[sA4R5XC18H10] + c[sFLTNXC16H10] + c[sC5H6]
		 + c[sC5H5] + c[sSXC3H5] + c[sCH4] + c[sC4H6]
		 + c[sTXC5H5O] + c[sC5H4O] + c[sSXC5H5O] + c[sC9H8]
		 + c[sC9H7] + c[sA1CH2XC7H7] + c[sC9H6O] + c[sOXC6H4]
		 + c[sA1CH3XC7H8] + c[sA1OHXC6H6O] + c[sHOA1CH3XC7H8O] + c[sOA1CH3XC7H7O]
		 + c[sA1CH2OXC7H7O] + c[sA1CH2OHXC7H8O] + c[sA1CHOXC7H6O] + c[sA1OXC6H5O]
		 + c[sA1CH3YXC7H7] + c[sPXC3H4] + c[sA1C2H4XC8H9] + c[sA1C2H5XC8H10]
		 + c[sC8H9O2] + c[sC8H8OOH] + c[sOC8H7OOH] + c[sH2C2]
		 + c[sA1CH3CH3XC8H10] + c[sA1CH3CH2XC8H9] + c[sA1CH3CHOXC8H8O] + c[sA2CH3XC11H10]
		 + c[sA1CHOCH2XC8H7O] + c[sA1CHOCHOXC8H6O2] + c[sA2OHXC10H8O] + c[sA2CH2XC11H9]
		 + c[sA2CH2OXC11H9O] + c[sA2CHOXC11H8O] + c[sA2OXC10H7O] + c[sC4H2]
		 + c[sOC6H4O] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sC]
		 + c[sHCNN] + c[sHCNO] + c[sNH3];

	M[mM20] = c[sN2] + c[sAR] + c[sO] + 2 * c[sH2]
		 + c[sH] + c[sOH] + 12 * c[sH2O] + c[sO2]
		 + c[sHO2] + c[sH2O2] + c[sCH] + c[sTXCH2]
		 + c[sCH2O] + c[sHCO] + c[sCH3] + 3.6 * c[sCO2]
		 + 1.75 * c[sCO] + c[sSXCH2] + c[sC2H4] + c[sCH3O]
		 + c[sCH3O2] + c[sC2H5] + c[sCH2CO] + c[sO2CHO]
		 + c[sHOCHO] + c[sOCHO] + c[sC2H2] + c[sC2H3]
		 + c[sHCCO] + c[sCH2CHO] + c[sC2H5O] + c[sAXC3H5]
		 + c[sC3H6] + c[sCH3COCH2] + c[sNXC3H7] + c[sPXC4H8]
		 + c[sPXC4H9] + c[sC5H11] + c[sC6H12X1XC6H12] + c[sC6H13]
		 + c[sC7H14] + c[sC7H15] + c[sC8H17] + c[sC9H19]
		 + c[sNXC12H26] + c[sC12H25] + c[sC12H24] + c[sC12H25O2]
		 + c[sC12OOHXT6XC12H25O2] + c[sC12OOHO2XT6XC12H25O4] + c[sC12KETXT6XC12H24O3] + c[sNXC8H17CHO]
		 + c[sC7H15CO] + c[sC3H7CHO] + c[sNXC7H15CHO] + c[sC7H15O2X1XC7H15O2]
		 + c[sC7H14OOH1X3XC7H15O2] + c[sC5H4CH2] + c[sA1XC6H6] + c[sNXC4H5]
		 + c[sIXC4H5] + c[sA1XXC6H5] + c[sNXC4H3] + c[sC3H3]
		 + c[sA1C2H2XC8H7] + c[sA1C2H3XC8H8] + c[sA1C2HXC8H6] + c[sA1C2HYXC8H5]
		 + c[sA1C2H3YXC8H7] + c[sA2XXC10H7] + c[sA2XC10H8] + c[sC4H4]
		 + c[sA2YXC10H7] + c[sA2C2H2AXC12H9] + c[sA2C2H2BXC12H9] + c[sA2C2HAXC12H8]
		 + c[sA2C2HBXC12H8] + c[sA2C2HAYXC12H7] + c[sA2C2HBYXC12H7] + c[sA2R5XC12H8]
		 + c[sA2R5XXC12H7] + c[sA2R5C2H2XC14H9] + c[sA2R5C2HXC14H8] + c[sA2R5C2HYXC14H7]
		 + c[sP2XC12H10] + c[sP2XXC12H9] + c[sA3XXC14H9] + c[sA3XC14H10]
		 + c[sA3YXC14H9] + c[sA3R5XXC16H9] + c[sA3R5XC16H10] + c[sA4XC16H10]
		 + c[sA4XXC16H9] + c[sA4R5XC18H10] + c[sFLTNXC16H10] + c[sC5H6]
		 + c[sC5H5] + c[sSXC3H5] + c[sCH4] + c[sC4H6]
		 + c[sTXC5H5O] + c[sC5H4O] + c[sSXC5H5O] + c[sC9H8]
		 + c[sC9H7] + c[sA1CH2XC7H7] + c[sC9H6O] + c[sOXC6H4]
		 + c[sA1CH3XC7H8] + c[sA1OHXC6H6O] + c[sHOA1CH3XC7H8O] + c[sOA1CH3XC7H7O]
		 + c[sA1CH2OXC7H7O] + c[sA1CH2OHXC7H8O] + c[sA1CHOXC7H6O] + c[sA1OXC6H5O]
		 + c[sA1CH3YXC7H7] + c[sPXC3H4] + c[sA1C2H4XC8H9] + c[sA1C2H5XC8H10]
		 + c[sC8H9O2] + c[sC8H8OOH] + c[sOC8H7OOH] + c[sH2C2]
		 + c[sA1CH3CH3XC8H10] + c[sA1CH3CH2XC8H9] + c[sA1CH3CHOXC8H8O] + c[sA2CH3XC11H10]
		 + c[sA1CHOCH2XC8H7O] + c[sA1CHOCHOXC8H6O2] + c[sA2OHXC10H8O] + c[sA2CH2XC11H9]
		 + c[sA2CH2OXC11H9O] + c[sA2CHOXC11H8O] + c[sA2OXC10H7O] + c[sC4H2]
		 + c[sOC6H4O] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sC]
		 + c[sHCNN] + c[sHCNO] + c[sNH3];

	M[mM21] = c[sN2] + c[sAR] + c[sO] + 2 * c[sH2]
		 + c[sH] + c[sOH] + c[sO2] + c[sHO2]
		 + c[sH2O2] + c[sCH] + c[sTXCH2] + c[sCH2O]
		 + c[sHCO] + c[sCH3] + 3.6 * c[sCO2] + 1.75 * c[sCO]
		 + c[sSXCH2] + c[sC2H4] + c[sCH3O] + c[sCH3O2]
		 + c[sC2H5] + c[sCH2CO] + c[sO2CHO] + c[sHOCHO]
		 + c[sOCHO] + c[sC2H2] + c[sC2H3] + c[sHCCO]
		 + c[sCH2CHO] + c[sC2H5O] + c[sAXC3H5] + c[sC3H6]
		 + c[sCH3COCH2] + c[sNXC3H7] + c[sPXC4H8] + c[sPXC4H9]
		 + c[sC5H11] + c[sC6H12X1XC6H12] + c[sC6H13] + c[sC7H14]
		 + c[sC7H15] + c[sC8H17] + c[sC9H19] + c[sNXC12H26]
		 + c[sC12H25] + c[sC12H24] + c[sC12H25O2] + c[sC12OOHXT6XC12H25O2]
		 + c[sC12OOHO2XT6XC12H25O4] + c[sC12KETXT6XC12H24O3] + c[sNXC8H17CHO] + c[sC7H15CO]
		 + c[sC3H7CHO] + c[sNXC7H15CHO] + c[sC7H15O2X1XC7H15O2] + c[sC7H14OOH1X3XC7H15O2]
		 + c[sC5H4CH2] + c[sA1XC6H6] + c[sNXC4H5] + c[sIXC4H5]
		 + c[sA1XXC6H5] + c[sNXC4H3] + c[sC3H3] + c[sA1C2H2XC8H7]
		 + c[sA1C2H3XC8H8] + c[sA1C2HXC8H6] + c[sA1C2HYXC8H5] + c[sA1C2H3YXC8H7]
		 + c[sA2XXC10H7] + c[sA2XC10H8] + c[sC4H4] + c[sA2YXC10H7]
		 + c[sA2C2H2AXC12H9] + c[sA2C2H2BXC12H9] + c[sA2C2HAXC12H8] + c[sA2C2HBXC12H8]
		 + c[sA2C2HAYXC12H7] + c[sA2C2HBYXC12H7] + c[sA2R5XC12H8] + c[sA2R5XXC12H7]
		 + c[sA2R5C2H2XC14H9] + c[sA2R5C2HXC14H8] + c[sA2R5C2HYXC14H7] + c[sP2XC12H10]
		 + c[sP2XXC12H9] + c[sA3XXC14H9] + c[sA3XC14H10] + c[sA3YXC14H9]
		 + c[sA3R5XXC16H9] + c[sA3R5XC16H10] + c[sA4XC16H10] + c[sA4XXC16H9]
		 + c[sA4R5XC18H10] + c[sFLTNXC16H10] + c[sC5H6] + c[sC5H5]
		 + c[sSXC3H5] + c[sCH4] + c[sC4H6] + c[sTXC5H5O]
		 + c[sC5H4O] + c[sSXC5H5O] + c[sC9H8] + c[sC9H7]
		 + c[sA1CH2XC7H7] + c[sC9H6O] + c[sOXC6H4] + c[sA1CH3XC7H8]
		 + c[sA1OHXC6H6O] + c[sHOA1CH3XC7H8O] + c[sOA1CH3XC7H7O] + c[sA1CH2OXC7H7O]
		 + c[sA1CH2OHXC7H8O] + c[sA1CHOXC7H6O] + c[sA1OXC6H5O] + c[sA1CH3YXC7H7]
		 + c[sPXC3H4] + c[sA1C2H4XC8H9] + c[sA1C2H5XC8H10] + c[sC8H9O2]
		 + c[sC8H8OOH] + c[sOC8H7OOH] + c[sH2C2] + c[sA1CH3CH3XC8H10]
		 + c[sA1CH3CH2XC8H9] + c[sA1CH3CHOXC8H8O] + c[sA2CH3XC11H10] + c[sA1CHOCH2XC8H7O]
		 + c[sA1CHOCHOXC8H6O2] + c[sA2OHXC10H8O] + c[sA2CH2XC11H9] + c[sA2CH2OXC11H9O]
		 + c[sA2CHOXC11H8O] + c[sA2OXC10H7O] + c[sC4H2] + c[sOC6H4O]
		 + c[sN] + c[sNO] + c[sN2O] + c[sNO2]
		 + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH]
		 + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN]
		 + c[sHNCO] + c[sH2CN] + c[sC] + c[sHCNN]
		 + c[sHCNO] + c[sNH3];

	M[mM32] = c[sN2] + c[sAR] + c[sO] + 2 * c[sH2]
		 + c[sH] + c[sOH] + 6 * c[sH2O] + c[sO2]
		 + c[sHO2] + c[sH2O2] + c[sCH] + c[sTXCH2]
		 + c[sCH2O] + c[sHCO] + c[sCH3] + 2 * c[sCO2]
		 + 1.5 * c[sCO] + c[sSXCH2] + c[sC2H4] + c[sCH3O]
		 + c[sCH3O2] + c[sC2H5] + c[sCH2CO] + c[sO2CHO]
		 + c[sHOCHO] + c[sOCHO] + c[sC2H2] + c[sC2H3]
		 + c[sHCCO] + c[sCH2CHO] + c[sC2H5O] + c[sAXC3H5]
		 + c[sC3H6] + c[sCH3COCH2] + c[sNXC3H7] + c[sPXC4H8]
		 + c[sPXC4H9] + c[sC5H11] + c[sC6H12X1XC6H12] + c[sC6H13]
		 + c[sC7H14] + c[sC7H15] + c[sC8H17] + c[sC9H19]
		 + c[sNXC12H26] + c[sC12H25] + c[sC12H24] + c[sC12H25O2]
		 + c[sC12OOHXT6XC12H25O2] + c[sC12OOHO2XT6XC12H25O4] + c[sC12KETXT6XC12H24O3] + c[sNXC8H17CHO]
		 + c[sC7H15CO] + c[sC3H7CHO] + c[sNXC7H15CHO] + c[sC7H15O2X1XC7H15O2]
		 + c[sC7H14OOH1X3XC7H15O2] + c[sC5H4CH2] + c[sA1XC6H6] + c[sNXC4H5]
		 + c[sIXC4H5] + c[sA1XXC6H5] + c[sNXC4H3] + c[sC3H3]
		 + c[sA1C2H2XC8H7] + c[sA1C2H3XC8H8] + c[sA1C2HXC8H6] + c[sA1C2HYXC8H5]
		 + c[sA1C2H3YXC8H7] + c[sA2XXC10H7] + c[sA2XC10H8] + c[sC4H4]
		 + c[sA2YXC10H7] + c[sA2C2H2AXC12H9] + c[sA2C2H2BXC12H9] + c[sA2C2HAXC12H8]
		 + c[sA2C2HBXC12H8] + c[sA2C2HAYXC12H7] + c[sA2C2HBYXC12H7] + c[sA2R5XC12H8]
		 + c[sA2R5XXC12H7] + c[sA2R5C2H2XC14H9] + c[sA2R5C2HXC14H8] + c[sA2R5C2HYXC14H7]
		 + c[sP2XC12H10] + c[sP2XXC12H9] + c[sA3XXC14H9] + c[sA3XC14H10]
		 + c[sA3YXC14H9] + c[sA3R5XXC16H9] + c[sA3R5XC16H10] + c[sA4XC16H10]
		 + c[sA4XXC16H9] + c[sA4R5XC18H10] + c[sFLTNXC16H10] + c[sC5H6]
		 + c[sC5H5] + c[sSXC3H5] + 2 * c[sCH4] + c[sC4H6]
		 + c[sTXC5H5O] + c[sC5H4O] + c[sSXC5H5O] + c[sC9H8]
		 + c[sC9H7] + c[sA1CH2XC7H7] + c[sC9H6O] + c[sOXC6H4]
		 + c[sA1CH3XC7H8] + c[sA1OHXC6H6O] + c[sHOA1CH3XC7H8O] + c[sOA1CH3XC7H7O]
		 + c[sA1CH2OXC7H7O] + c[sA1CH2OHXC7H8O] + c[sA1CHOXC7H6O] + c[sA1OXC6H5O]
		 + c[sA1CH3YXC7H7] + c[sPXC3H4] + c[sA1C2H4XC8H9] + c[sA1C2H5XC8H10]
		 + c[sC8H9O2] + c[sC8H8OOH] + c[sOC8H7OOH] + c[sH2C2]
		 + c[sA1CH3CH3XC8H10] + c[sA1CH3CH2XC8H9] + c[sA1CH3CHOXC8H8O] + c[sA2CH3XC11H10]
		 + c[sA1CHOCH2XC8H7O] + c[sA1CHOCHOXC8H6O2] + c[sA2OHXC10H8O] + c[sA2CH2XC11H9]
		 + c[sA2CH2OXC11H9O] + c[sA2CHOXC11H8O] + c[sA2OXC10H7O] + c[sC4H2]
		 + c[sOC6H4O] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sC]
		 + c[sHCNN] + c[sHCNO] + c[sNH3];



	k[r3f] = -80207000 *rt_inv;
	k[r3b] = -0.0350345 * lgt - 74440909.81 * rt_inv ;
	k[r4f] = 1.51 * lgt - 14351000 * rt_inv ;
	k[r4b] = 1.41638 * lgt - 77059304.43 * rt_inv ;
	k[r5f] = 2.42 * lgt + 8075000 * rt_inv ;
	k[r5b] = 2.36141 * lgt - 60399394.61 * rt_inv ;
	k[r1f] = -63957000 *rt_inv;
	k[r1b] = 0.434958 * lgt + 6604131.16 * rt_inv ;
	k[r16] = 2.09 * lgt + 6071000 * rt_inv ;
	k[r21] = 6817000 *rt_inv;
	k[r20] = -50133000 *rt_inv;
	k[r19] = 2079000 *rt_inv;
	k[r18] = 1 * lgt + 3029000 * rt_inv ;
	k[r17] = -1234000 *rt_inv;
	k[r26] = -1331000 *rt_inv;
	k[r23] = -16610000 *rt_inv;
	k[r39] = -30418000 *rt_inv;
	k[rG06f] = -13010000 *rt_inv;
	k[rG06b] = -0.291033 * lgt - 1692273.875 * rt_inv ;
	k[rG08] = 3160000 *rt_inv;
	k[rG19] = -6280000 *rt_inv;
	k[rG20] = -6280000 *rt_inv;
	k[rG17] = 2 * lgt - 12550000 * rt_inv ;
	k[rG21] = -6280000 *rt_inv;
	k[rG38b] = -0.0622863 * lgt - 37302987.98 * rt_inv ;
	k[rG34b] = -0.573354 * lgt - 68461920.07 * rt_inv ;
	k[rG28f] = -2510000 *rt_inv;
	k[rG28b] = -0.0622863 * lgt - 39812987.98 * rt_inv ;
	k[rG55] = 1.6 * lgt - 22680000 * rt_inv ;
	k[rG57f] = -1.34 * lgt - 5930000 * rt_inv ;
	k[rG57b] = -0.86027 * lgt - 176384.3555 * rt_inv ;
	k[rG59] = -57910000 *rt_inv;
	k[rG74f] = 0.1 * lgt - 44350000 * rt_inv ;
	k[rG74b] = -1.0902 * lgt - 7214861.361 * rt_inv ;
	k[r29] = -0.664 * lgt - 1388000 * rt_inv ;
	k[r31] = 2.18 * lgt - 75072000 * rt_inv ;
	k[r28f] = 2.053 * lgt + 1488000 * rt_inv ;
	k[r28b] = 0.773803 * lgt - 105981755.2 * rt_inv ;
	k[r38] = -1715000 *rt_inv;
	k[r36f] = 0.659 * lgt - 62232000 * rt_inv ;
	k[r36b] = 0.891456 * lgt + 3517584.76 * rt_inv ;
	k[rDD68] = 4602000 *rt_inv;
	k[rG46] = -14810000 *rt_inv;
	k[rG47] = 1.18 * lgt + 1870000 * rt_inv ;
	k[rG45] = 1.9 * lgt - 11470000 * rt_inv ;
	k[rG49] = 2 * lgt - 50210000 * rt_inv ;
	k[rLPDD62] =  -1.1 * lgt;
	k[rLPG89] = -0.168 * lgt - 61525000 * rt_inv ;
	k[rDD44] = -0.353 * lgt - 73555000 * rt_inv ;
	k[rDD47] = 1.51 * lgt + 4025000 * rt_inv ;
	k[rDD50] = 2.1 * lgt - 20368000 * rt_inv ;
	k[rDD48] = -49873000 *rt_inv;
	k[rDD53] = -1.542 * lgt - 463169000 * rt_inv ;
	k[rDD52] = 2.06 * lgt - 3833000 * rt_inv ;
	k[rLPG63] = -0.019 * lgt + 6396000 * rt_inv ;
	k[rLPDD59] = 18.783 * lgt + 37518000 * rt_inv ;
	k[rLPG61] = -0.019 * lgt + 4846000 * rt_inv ;
	k[rLPG64] = -0.019 * lgt - 49004000 * rt_inv ;
	k[rDD66] = -48785000 *rt_inv;
	k[rDD67] = -1.126 * lgt - 171544000 * rt_inv ;
	k[rG116] = 2 * lgt - 7950000 * rt_inv ;
	k[rG117] = 2 * lgt - 7950000 * rt_inv ;
	k[rLPG120f] = 3.703 * lgt + 7342000 * rt_inv ;
	k[rLPG120b] = -4.52 * lgt - 266751000 * rt_inv ;
	k[rG122] = 0.73 * lgt - 10790000 * rt_inv ;
	k[rG135] = -1.39 * lgt - 4250000 * rt_inv ;
	k[rLPG134] = 0.29 * lgt - 50000 * rt_inv ;
	k[rG133] = 1.61 * lgt + 1610000 * rt_inv ;
	k[rLPG131] = 0.21 * lgt + 1790000 * rt_inv ;
	k[rG160] = 4.2 * lgt + 3600000 * rt_inv ;
	k[rG158] = 2.47 * lgt - 3890000 * rt_inv ;
	k[rG159] = 1.36 * lgt - 3710000 * rt_inv ;
	k[rLPG157] = 0.88 * lgt - 4770000 * rt_inv ;
	k[rG156] = 2.75 * lgt - 48740000 * rt_inv ;
	k[rG161] = -7.8 * lgt - 29540000 * rt_inv ;
	k[rLPDD15] = -0.019 * lgt + 4010000 * rt_inv ;
	k[rG171] = 1.02 * lgt + 8510000 * rt_inv ;
	k[rG111] = -3570000 *rt_inv;
	k[rLPG124] = -0.048 * lgt - 15621000 * rt_inv ;
	k[rLPR028] = -1.428 * lgt - 13051000 * rt_inv ;
	k[rLPG123f] = -1.428 * lgt - 46521000 * rt_inv ;
	k[rLPG123b] = 0.596 * lgt - 22499000 * rt_inv ;
	k[rLPG127f] = -1.428 * lgt - 21421000 * rt_inv ;
	k[rLPG127b] = 0.502 * lgt - 60107000 * rt_inv ;
	k[rLPR025] = -1.446 * lgt - 13299000 * rt_inv ;
	k[rLPR030] = -0.509 * lgt - 119000 * rt_inv ;
	k[rLPG139] = -0.009 * lgt - 119000 * rt_inv ;
	k[rLPG140] = -0.009 * lgt - 119000 * rt_inv ;
	k[rLPG137] = -6.879 * lgt - 197579000 * rt_inv ;
	k[rLPG138] = 0.021 * lgt + 1531000 * rt_inv ;
	k[rLPG136] = -6.579 * lgt - 207051000 * rt_inv ;
	k[rLPG144] = -0.009 * lgt - 128000 * rt_inv ;
	k[rG168] = -2.02 * lgt - 86820000 * rt_inv ;
	k[rLPR015] = -3.805 * lgt - 106695000 * rt_inv ;
	k[rLPR116] = -0.619 * lgt - 5490000 * rt_inv ;
	k[rLPR115] = -0.087 * lgt - 726000 * rt_inv ;
	k[rLPR133] = 4.164 * lgt - 41693000 * rt_inv ;
	k[rLPLPR118] = -0.508 * lgt - 248000 * rt_inv ;
	k[rLPR016] = -12.379 * lgt - 104020000 * rt_inv ;
	k[rLPLPR113] = -0.589 * lgt - 3840000 * rt_inv ;
	k[rLPR111] = -3.319 * lgt - 110010000 * rt_inv ;
	k[rLPR007] = -6.861 * lgt - 224777000 * rt_inv ;
	k[rLPLPR117] = -1.559 * lgt - 16050000 * rt_inv ;
	k[rR141] = -2.39 * lgt - 46780000 * rt_inv ;
	k[rR143] = 1.6 * lgt + 4070000 * rt_inv ;
	k[rLPR153] = 2.926 * lgt - 27350000 * rt_inv ;
	k[rLPR142] = 1.6 * lgt - 1370000 * rt_inv ;
	k[rLPR151] = 2.776 * lgt + 6471000 * rt_inv ;
	k[rG184] = -1.63 * lgt - 14300000 * rt_inv ;
	k[rLPHP55] = -0.493 * lgt - 969000 * rt_inv ;
	k[rLPHP43] = -0.434 * lgt - 297422000 * rt_inv ;
	k[rLPHP54] = -0.493 * lgt - 4529000 * rt_inv ;
	k[rLPHP52] = -0.493 * lgt - 4529000 * rt_inv ;
	k[rLPHP53] = 1.847 * lgt + 3421000 * rt_inv ;
	k[rLPDD514] = 0.015 * lgt - 5525000 * rt_inv ;
	k[rLPHP56] = -0.493 * lgt - 969000 * rt_inv ;
	k[rLPHP39] = -1.4 * lgt - 162800000 * rt_inv ;
	k[rHP41] = -0.64 * lgt - 112570000 * rt_inv ;
	k[rHP29] = -2.61 * lgt - 134000000 * rt_inv ;
	k[rHP32] = 8.84 * lgt - 29730000 * rt_inv ;
	k[rDD199] = 16736000 *rt_inv;
	k[rDD212] = 0.086 * lgt - 6452000 * rt_inv ;
	k[rLPDD205] = -297064000 *rt_inv;
	k[rDD214] = -3.686 * lgt - 142868000 * rt_inv ;
	k[rDD208] = 0.279 * lgt - 128848000 * rt_inv ;
	k[rDD43] = 16736000 *rt_inv;
	k[rHP26] = -1.65 * lgt + 7630000 * rt_inv ;
	k[rLPHP22] = -297064000 *rt_inv;
	k[rHP15] = 0.45 * lgt - 8430000 * rt_inv ;
	k[rLPHP12] = -1.65 * lgt - 132560000 * rt_inv ;
	k[rHP11] = -1.75 * lgt - 133780000 * rt_inv ;
	k[rHP10] = 0.02 * lgt - 116250000 * rt_inv ;
	k[rLPDD286] = -1.413 * lgt - 110770000 * rt_inv ;
	k[rDD282] = -0.833 * lgt - 118041000 * rt_inv ;
	k[rDD288] = 0.336 * lgt - 100547000 * rt_inv ;
	k[rDD287] = -2.405 * lgt - 160016000 * rt_inv ;
	k[rDD323] = -1.902 * lgt - 118943000 * rt_inv ;
	k[rDD324] = -2.324 * lgt - 159951000 * rt_inv ;
	k[rDD319] = 0.693 * lgt - 101293000 * rt_inv ;
	k[rDD210] = -2.246 * lgt - 370781000 * rt_inv ;
	k[rDD421] = -2.223 * lgt - 371839000 * rt_inv ;
	k[rDD436] = -2.147 * lgt - 370881000 * rt_inv ;
	k[rDD425] = -2.229 * lgt - 371576000 * rt_inv ;
	k[rDD422] = -207694000 *rt_inv;
	k[rDD434] = 1.544 * lgt - 556000 * rt_inv ;
	k[rDD429] = 3.557 * lgt - 56660000 * rt_inv ;
	k[rDD424] = 2.52 * lgt - 12080000 * rt_inv ;
	k[rDD432] = 2.521 * lgt - 18533000 * rt_inv ;
	k[rDD427] = -74015000 *rt_inv;
	k[rLPDD431] = 3.73 * lgt - 61347000 * rt_inv ;
	k[rDD405f] = 0.309 * lgt - 152536000 * rt_inv ;
	k[rDD405b] = 2.187 * lgt + 7611000 * rt_inv ;
	k[rDD418] = 0.51 * lgt - 113311000 * rt_inv ;
	k[rDD415] = -1.265 * lgt - 131866000 * rt_inv ;
	k[rDD416] = -0.397 * lgt - 121439000 * rt_inv ;
	k[rLPDD406] = -1.223 * lgt - 128372000 * rt_inv ;
	k[rDD408f] = -0.324 * lgt - 2235000 * rt_inv ;
	k[rDD408b] = -2.989 * lgt - 165486000 * rt_inv ;
	k[rDD472] = -7.955 * lgt - 166778000 * rt_inv ;
	k[rLPDD465f] = 1.188 * lgt - 77083000 * rt_inv ;
	k[rLPDD465b] = 0.937827 * lgt - 23205657.49 * rt_inv ;
	k[rLPDD473] = 4.946 * lgt - 24970000 * rt_inv ;
	k[rLPDD486f] = -1.558 * lgt - 14198000 * rt_inv ;
	k[rLPDD486b] = -4.19719 * lgt - 177495502.4 * rt_inv ;
	k[rLPDD502] = -0.852 * lgt - 81442000 * rt_inv ;
	k[rLPLPDD495] = -1.292 * lgt - 171906000 * rt_inv ;
	k[rLPDD499] = -2.155 * lgt - 178393000 * rt_inv ;
	k[rLPDD490] = -1.679 * lgt - 172679000 * rt_inv ;
	k[rLPR140] = -14.011 * lgt - 294551000 * rt_inv ;
	k[rLPLPDD107] = 9.495 * lgt + 41108000 * rt_inv ;
	k[rDD332N] = 0.76 * lgt + 1423000 * rt_inv ;
	k[rDD294] = 0.76 * lgt + 1423000 * rt_inv ;
	k[rDD290] = -119179000 *rt_inv;
	k[rLPDD293] = -3.723 * lgt - 152644000 * rt_inv ;
	k[rHP60] = 0.76 * lgt + 1420000 * rt_inv ;
	k[rDD239] = -1.716 * lgt - 14670000 * rt_inv ;
	k[rDD263] = -2.149 * lgt - 18957000 * rt_inv ;
	k[rDD260f] = -85563000 *rt_inv;
	k[rDD260b] = -0.23424 * lgt - 33307702.91 * rt_inv ;
	k[rDD259] = 1.23 * lgt - 127068000 * rt_inv ;
	k[rP000f] = 0.5 * lgt - 8370000 * rt_inv ;
	k[rP000b] = 0.794062 * lgt - 150584803.3 * rt_inv ;
	k[rP001f] = -0.89 * lgt - 38250000 * rt_inv ;
	k[rP001b] = -2.30141 * lgt - 194882567.4 * rt_inv ;
	k[rP002f] = -3.45 * lgt - 85090000 * rt_inv ;
	k[rP002b] = -4.51606 * lgt - 195705788.7 * rt_inv ;
	k[rP003f] = -8.9 * lgt - 405860000 * rt_inv ;
	k[rP003b] = -8.60594 * lgt - 548074803.3 * rt_inv ;
	k[rP004f] = -14.65 * lgt - 596540000 * rt_inv ;
	k[rP004b] = -14.1444 * lgt - 262723022 * rt_inv ;
	k[rP005f] = -17.77 * lgt - 130960000 * rt_inv ;
	k[rP005b] = -19.121 * lgt - 575688821.6 * rt_inv ;
	k[rP006f] = -1 * lgt - 37240000 * rt_inv ;
	k[rP006b] = -2.11735 * lgt - 336087370.7 * rt_inv ;
	k[rP007f] = -3.3 * lgt - 104430000 * rt_inv ;
	k[rP007b] = -4.07199 * lgt - 357260592 * rt_inv ;
	k[rP008] = -7.74 * lgt - 99800000 * rt_inv ;
	k[rP009f] = -10.1 * lgt - 70960000 * rt_inv ;
	k[rP009b] = -11.8517 * lgt - 535008868.8 * rt_inv ;
	k[rP010f] = -9.57 * lgt - 71190000 * rt_inv ;
	k[rP010b] = -11.0277 * lgt - 677453672.1 * rt_inv ;
	k[rP011f] = -7 * lgt - 131820000 * rt_inv ;
	k[rP011b] = -8.24617 * lgt - 262051890.7 * rt_inv ;
	k[rK012f] = 2.05 * lgt - 13230000 * rt_inv ;
	k[rK012b] = 0.27021 * lgt - 191386104.5 * rt_inv ;
	k[rP013b] = -2.06437 * lgt - 497693591.2 * rt_inv ;
	k[rP014f] = 1.47 * lgt - 23150000 * rt_inv ;
	k[rP014b] = -0.382874 * lgt - 44811809.88 * rt_inv ;
	k[rP015f] = 4.47 * lgt - 18710000 * rt_inv ;
	k[rP015b] = 4.39414 * lgt - 28239401.86 * rt_inv ;
	k[rP016f] = -12.48 * lgt - 619590000 * rt_inv ;
	k[rP016b] = -12.2493 * lgt - 132929473 * rt_inv ;
	k[rK017f] = 1.89 * lgt - 73550000 * rt_inv ;
	k[rK017b] = 2.46953 * lgt - 19962836.49 * rt_inv ;
	k[rP018f] = 2.68 * lgt - 3070000 * rt_inv ;
	k[rP018b] = 3.16591 * lgt - 12191140.91 * rt_inv ;
	k[rK019f] = 0.54 * lgt - 115340000 * rt_inv ;
	k[rK019b] = 0.647444 * lgt - 97958978.78 * rt_inv ;
	k[rK020f] = 1.02 * lgt - 161810000 * rt_inv ;
	k[rK020b] = 1.32078 * lgt - 30401375.45 * rt_inv ;
	k[rK021f] = 0.49 * lgt - 44480000 * rt_inv ;
	k[rK021b] = 1.13958 * lgt - 346144739 * rt_inv ;
	k[rP022b] = 0.555961 * lgt - 364373043.4 * rt_inv ;
	k[rP023f] = -12.48 * lgt - 619590000 * rt_inv ;
	k[rP023b] = -12.2599 * lgt - 135887131.6 * rt_inv ;
	k[rK024f] = 2.36 * lgt - 70780000 * rt_inv ;
	k[rK024b] = 2.92894 * lgt - 20150495.14 * rt_inv ;
	k[rP025f] = 3.33 * lgt - 6090000 * rt_inv ;
	k[rP025b] = 3.80532 * lgt - 18168799.56 * rt_inv ;
	k[rP026f] = 0.34 * lgt - 465490000 * rt_inv ;
	k[rP026b] = 0.452698 * lgt + 831847.1725 * rt_inv ;
	k[rP027f] = 2.75 * lgt - 48740000 * rt_inv ;
	k[rP027b] = 3.2115 * lgt - 15491516.35 * rt_inv ;
	k[rP028f] = 4.2 * lgt + 3600000 * rt_inv ;
	k[rP028b] = 4.56787 * lgt - 25859820.78 * rt_inv ;
	k[rK100f] = 2.5 * lgt - 5370000 * rt_inv ;
	k[rK100b] = 0.979368 * lgt - 386582470.9 * rt_inv ;
	k[rK102f] = 2.55 * lgt - 13310000 * rt_inv ;
	k[rK102b] = 1.2202 * lgt - 292022356 * rt_inv ;
	k[rP104f] = 1.62 * lgt - 18570000 * rt_inv ;
	k[rP104b] = 0.397642 * lgt - 279901334.8 * rt_inv ;
	k[rP105f] = -1.44 * lgt - 65930000 * rt_inv ;
	k[rP105b] = -3.13503 * lgt - 311885598.9 * rt_inv ;
	k[rP106f] = -4.24 * lgt - 99850000 * rt_inv ;
	k[rP106b] = -6.03013 * lgt - 365963746.5 * rt_inv ;
	k[rP107f] = 2.61 * lgt - 6000000 * rt_inv ;
	k[rP107b] = 0.803868 * lgt - 279700375 * rt_inv ;
	k[rP108f] = -12.48 * lgt - 619550000 * rt_inv ;
	k[rP108b] = -12.2468 * lgt - 121361984.6 * rt_inv ;
	k[rK109f] = 1.87 * lgt - 71530000 * rt_inv ;
	k[rK109b] = 2.45204 * lgt - 6415348.089 * rt_inv ;
	k[rK110f] = 3.02 * lgt - 18300000 * rt_inv ;
	k[rK110b] = 3.50842 * lgt - 15893652.52 * rt_inv ;
	k[rP111f] = -12.48 * lgt - 619550000 * rt_inv ;
	k[rP111b] = -12.2625 * lgt - 122653241.3 * rt_inv ;
	k[rK112f] = 1.87 * lgt - 71530000 * rt_inv ;
	k[rK112b] = 2.43625 * lgt - 7706604.8 * rt_inv ;
	k[rK113f] = 3.02 * lgt - 18300000 * rt_inv ;
	k[rK113b] = 3.49263 * lgt - 17184909.23 * rt_inv ;
	k[rK114f] = 2.05 * lgt - 8080000 * rt_inv ;
	k[rK114b] = 0.192174 * lgt - 181300070.2 * rt_inv ;
	k[rK115f] = 2.05 * lgt - 13230000 * rt_inv ;
	k[rK115b] = 0.287919 * lgt - 186404982.5 * rt_inv ;
	k[rP116b] = -2.02971 * lgt - 26435709.78 * rt_inv ;
	k[rP117b] = -1.93396 * lgt - 26390622.05 * rt_inv ;
	k[rP118f] = -1.44 * lgt - 65930000 * rt_inv ;
	k[rP118b] = -2.88767 * lgt - 27251057.87 * rt_inv ;
	k[rP119f] = -1.44 * lgt - 65930000 * rt_inv ;
	k[rP119b] = -2.80771 * lgt - 28497226.85 * rt_inv ;
	k[rP120f] = -4.24 * lgt - 99850000 * rt_inv ;
	k[rP120b] = -5.78528 * lgt - 92856693.82 * rt_inv ;
	k[rP121f] = -4.24 * lgt - 99850000 * rt_inv ;
	k[rP121b] = -5.68953 * lgt - 92811606.1 * rt_inv ;
	k[rK122f] = 1.02 * lgt - 161810000 * rt_inv ;
	k[rK122b] = 1.38384 * lgt - 60866185.59 * rt_inv ;
	k[rK123f] = 0.49 * lgt - 44480000 * rt_inv ;
	k[rK123b] = 1.20264 * lgt - 376609549.1 * rt_inv ;
	k[rP124b] = 0.619019 * lgt - 394837853.5 * rt_inv ;
	k[rK125f] = 1.02 * lgt - 161810000 * rt_inv ;
	k[rK125b] = 1.29787 * lgt - 58106247.65 * rt_inv ;
	k[rK126f] = 0.49 * lgt - 44480000 * rt_inv ;
	k[rK126b] = 1.11667 * lgt - 373849611.2 * rt_inv ;
	k[rP127b] = 0.533047 * lgt - 392077915.6 * rt_inv ;
	k[rP128f] = -12.48 * lgt - 619550000 * rt_inv ;
	k[rP128b] = -12.2367 * lgt - 111810655.9 * rt_inv ;
	k[rK129f] = 1.88 * lgt - 70380000 * rt_inv ;
	k[rK129b] = 2.47205 * lgt + 4285980.552 * rt_inv ;
	k[rP130f] = 3.33 * lgt - 6090000 * rt_inv ;
	k[rP130b] = 3.82843 * lgt + 5867676.125 * rt_inv ;
	k[rP131f] = -12.48 * lgt - 619550000 * rt_inv ;
	k[rP131b] = -12.225 * lgt - 110360350.9 * rt_inv ;
	k[rK132f] = 1.88 * lgt - 70380000 * rt_inv ;
	k[rK132b] = 2.48384 * lgt + 5736285.547 * rt_inv ;
	k[rP133f] = 3.33 * lgt - 6090000 * rt_inv ;
	k[rP133b] = 3.84022 * lgt + 7317981.12 * rt_inv ;
	k[rK200f] = 0.23 * lgt - 71240000 * rt_inv ;
	k[rK200b] = 0.811179 * lgt - 85017550.86 * rt_inv ;
	k[rP201f] = -55900000 *rt_inv;
	k[rP201b] = 0.217336 * lgt - 170621365.3 * rt_inv ;
	k[rP202f] = -12.48 * lgt - 619550000 * rt_inv ;
	k[rP202b] = -12.2259 * lgt - 125068499.8 * rt_inv ;
	k[rK203f] = 1.87 * lgt - 71530000 * rt_inv ;
	k[rK203b] = 2.47285 * lgt - 10121863.29 * rt_inv ;
	k[rK204f] = 3.02 * lgt - 18300000 * rt_inv ;
	k[rK204b] = 3.52923 * lgt - 19600167.72 * rt_inv ;
	k[rK205f] = 2.05 * lgt - 13230000 * rt_inv ;
	k[rK205b] = -0.322738 * lgt - 188836244.6 * rt_inv ;
	k[rP206b] = -2.54462 * lgt - 28821884.2 * rt_inv ;
	k[rP207f] = -1.44 * lgt - 65930000 * rt_inv ;
	k[rP207b] = -3.38177 * lgt - 33343747.49 * rt_inv ;
	k[rP208f] = -4.24 * lgt - 99850000 * rt_inv ;
	k[rP208b] = -6.30019 * lgt - 95242868.24 * rt_inv ;
	k[rP209f] = -12.48 * lgt - 619550000 * rt_inv ;
	k[rP209b] = -12.2371 * lgt - 118358316.7 * rt_inv ;
	k[rK210f] = 1.88 * lgt - 70380000 * rt_inv ;
	k[rK210b] = 2.47167 * lgt - 2261680.225 * rt_inv ;
	k[rP211f] = 3.33 * lgt - 6090000 * rt_inv ;
	k[rP211b] = 3.82804 * lgt - 679984.652 * rt_inv ;
	k[rK212f] = 1.02 * lgt - 161810000 * rt_inv ;
	k[rK212b] = 1.88558 * lgt - 52723839.08 * rt_inv ;
	k[rK213f] = 0.49 * lgt - 44480000 * rt_inv ;
	k[rK213b] = 1.70438 * lgt - 368467202.6 * rt_inv ;
	k[rP214b] = 1.12075 * lgt - 386695507 * rt_inv ;
	k[rP301f] = -9070000 *rt_inv;
	k[rP301b] = -1.57654 * lgt - 32778677.98 * rt_inv ;
	k[rP302f] = -470000 *rt_inv;
	k[rP302b] = -1.78804 * lgt - 500210459.3 * rt_inv ;
	k[rP304f] = -12.48 * lgt - 619550000 * rt_inv ;
	k[rP304b] = -12.4557 * lgt - 126953664.6 * rt_inv ;
	k[rP305f] = 1.8 * lgt - 68420000 * rt_inv ;
	k[rP305b] = 2.17308 * lgt - 8897028.137 * rt_inv ;
	k[rP306f] = 3.33 * lgt - 6090000 * rt_inv ;
	k[rP306b] = 3.60945 * lgt - 9275332.564 * rt_inv ;
	k[rK401f] = 2.5 * lgt - 5370000 * rt_inv ;
	k[rK401b] = 0.97143 * lgt - 409641149.6 * rt_inv ;
	k[rK403f] = 2.5 * lgt - 5370000 * rt_inv ;
	k[rK403b] = 0.96566 * lgt - 412605223.6 * rt_inv ;
	k[rP405f] = 1.62 * lgt - 18570000 * rt_inv ;
	k[rP405b] = 0.473091 * lgt - 323520364.2 * rt_inv ;
	k[rP406f] = 1.62 * lgt - 18570000 * rt_inv ;
	k[rP406b] = 0.393137 * lgt - 322274195.2 * rt_inv ;
	k[rP407f] = 2.05 * lgt - 13230000 * rt_inv ;
	k[rP407b] = 0.444306 * lgt - 286506547.6 * rt_inv ;
	k[rP408b] = -1.92589 * lgt - 766849162.3 * rt_inv ;
	k[rP409b] = -1.93166 * lgt - 769813236.2 * rt_inv ;
	k[rP410f] = -1.44 * lgt - 65930000 * rt_inv ;
	k[rP410b] = -3.12264 * lgt - 325039818.2 * rt_inv ;
	k[rP411f] = -1.44 * lgt - 65930000 * rt_inv ;
	k[rP411b] = -3.11662 * lgt - 326553587.1 * rt_inv ;
	k[rP412f] = -4.24 * lgt - 99850000 * rt_inv ;
	k[rP412b] = -6.03025 * lgt - 400196782.8 * rt_inv ;
	k[rP413f] = -4.24 * lgt - 99850000 * rt_inv ;
	k[rP413b] = -6.03603 * lgt - 403160856.7 * rt_inv ;
	k[rP414f] = 2.61 * lgt - 6000000 * rt_inv ;
	k[rP414b] = 0.801281 * lgt - 318383370.2 * rt_inv ;
	k[rP415f] = 2.61 * lgt - 6000000 * rt_inv ;
	k[rP415b] = 0.817072 * lgt - 317092113.5 * rt_inv ;
	k[rP416f] = -9070000 *rt_inv;
	k[rP416b] = -1.89045 * lgt - 242743191.6 * rt_inv ;
	k[rP417f] = -9070000 *rt_inv;
	k[rP417b] = -1.90968 * lgt - 253371937.3 * rt_inv ;
	k[rP418f] = -470000 *rt_inv;
	k[rP418b] = -2.12118 * lgt - 720803718.6 * rt_inv ;
	k[rP419f] = -12.48 * lgt - 619550000 * rt_inv ;
	k[rP419b] = -12.2546 * lgt - 110187626.9 * rt_inv ;
	k[rK420f] = 1.8 * lgt - 68420000 * rt_inv ;
	k[rK420b] = 2.37424 * lgt + 7869009.54 * rt_inv ;
	k[rK421f] = 3.33 * lgt - 6090000 * rt_inv ;
	k[rK421b] = 3.81061 * lgt + 7490705.113 * rt_inv ;
	k[rP422f] = -12.48 * lgt - 619550000 * rt_inv ;
	k[rP422b] = -12.2546 * lgt - 110187626.9 * rt_inv ;
	k[rK423f] = 1.87 * lgt - 71530000 * rt_inv ;
	k[rK423b] = 2.44424 * lgt + 4759009.54 * rt_inv ;
	k[rK424f] = 3.02 * lgt - 18300000 * rt_inv ;
	k[rK424b] = 3.50061 * lgt - 4719294.887 * rt_inv ;
	k[rP425f] = 1.08 * lgt - 294550000 * rt_inv ;
	k[rP425b] = 2.8367 * lgt - 18258059.5 * rt_inv ;
	k[rP501f] = 2.5 * lgt - 5370000 * rt_inv ;
	k[rP501b] = 0.969451 * lgt - 413991985.3 * rt_inv ;
	k[rP502f] = 1.62 * lgt - 18570000 * rt_inv ;
	k[rP502b] = 0.969942 * lgt - 321579572.7 * rt_inv ;
	k[rP503f] = 2.05 * lgt - 8080000 * rt_inv ;
	k[rP503b] = 0.783907 * lgt - 210403876.8 * rt_inv ;
	k[rP504b] = -1.08918 * lgt - 488612879.9 * rt_inv ;
	k[rP505f] = 2.61 * lgt - 6000000 * rt_inv ;
	k[rP505b] = 0.78322 * lgt - 318828753.1 * rt_inv ;
	k[rP506f] = -12.48 * lgt - 619550000 * rt_inv ;
	k[rP506b] = -12.252 * lgt - 114884568.4 * rt_inv ;
	k[rK507f] = 1.87 * lgt - 71530000 * rt_inv ;
	k[rK507b] = 2.44675 * lgt + 62068.03303 * rt_inv ;
	k[rK508f] = 3.02 * lgt - 18300000 * rt_inv ;
	k[rK508b] = 3.50313 * lgt - 9416236.394 * rt_inv ;
	k[rK600f] = 2.05 * lgt - 8080000 * rt_inv ;
	k[rK600b] = 0.882089 * lgt - 289602697.5 * rt_inv ;
	k[rP601f] = -12.48 * lgt - 619550000 * rt_inv ;
	k[rP601b] = -12.2341 * lgt - 85298616.63 * rt_inv ;
	k[rK602f] = 1.87 * lgt - 71530000 * rt_inv ;
	k[rK602b] = 2.46473 * lgt + 29648019.84 * rt_inv ;
	k[rK603f] = 3.02 * lgt - 18300000 * rt_inv ;
	k[rK603b] = 3.52111 * lgt + 20169715.41 * rt_inv ;
	k[rK700f] = 2.05 * lgt - 8080000 * rt_inv ;
	k[rK700b] = 0.76659 * lgt - 230176187.5 * rt_inv ;
	k[rK701f] = 2.05 * lgt - 8080000 * rt_inv ;
	k[rK701b] = 0.882751 * lgt - 279789056.4 * rt_inv ;
	k[rP800] = -9070000 *rt_inv;
	k[rP801] = -9070000 *rt_inv;
	k[rP802] = -470000 *rt_inv;
	k[rCP01f] = -15.16 * lgt - 486900000 * rt_inv ;
	k[rCP01b] = -14.464 * lgt - 137434558.3 * rt_inv ;
	k[rCP02f] = -9450000 *rt_inv;
	k[rCP02b] = 1.04483 * lgt - 93057921.84 * rt_inv ;
	k[rCP03f] = -51650000 *rt_inv;
	k[rCP03b] = 1.13526 * lgt - 2981345.18 * rt_inv ;
	k[rCP04] = -51650000 *rt_inv;
	k[rCP05f] = 2.71 * lgt - 4630000 * rt_inv ;
	k[rCP05b] = 3.7198 * lgt - 82471831.66 * rt_inv ;
	k[rCP06f] =  2 * lgt;
	k[rCP06b] = 2.95121 * lgt - 146316226.3 * rt_inv ;
	k[rCP07f] = -155440000 *rt_inv;
	k[rCP07b] = 0.683127 * lgt - 10637907.05 * rt_inv ;
	k[rCP08f] = 2.6 * lgt - 53970000 * rt_inv ;
	k[rCP08b] = 2.89385 * lgt - 71243579.16 * rt_inv ;
	k[rCP09f] =  4 * lgt;
	k[rCP09b] = 4.52606 * lgt - 92521076.66 * rt_inv ;
	k[rCP10b] = 0.560401 * lgt - 117036937.8 * rt_inv ;
	k[rCP11b] = -0.2405 * lgt - 123947539.5 * rt_inv ;
	k[rCP12f] = 0.88 * lgt - 4770000 * rt_inv ;
	k[rCP12b] = -0.124313 * lgt - 111518622.2 * rt_inv ;
	k[rCP13] = 1.36 * lgt - 3710000 * rt_inv ;
	k[rCP14] = -7.8 * lgt - 29540000 * rt_inv ;
	k[rCP15] = -7.8 * lgt - 29540000 * rt_inv ;
	k[rCP16f] = -12.5 * lgt - 176000000 * rt_inv ;
	k[rCP16b] = -13.5715 * lgt - 487565243.7 * rt_inv ;
	k[rCP17] = -4.03 * lgt - 147300000 * rt_inv ;
	k[rCP18] = -4.85 * lgt - 103650000 * rt_inv ;
	k[rCP19b] = -1.03661 * lgt - 240967609.8 * rt_inv ;
	k[rCP20f] = 1.3 * lgt - 73920000 * rt_inv ;
	k[rCP20b] = 0.698346 * lgt - 244326478.6 * rt_inv ;
	k[rCP21b] = -0.692564 * lgt - 18654658.39 * rt_inv ;
	k[rCP23f] = -125520000 *rt_inv;
	k[rCP23b] = 0.103818 * lgt - 72608471.51 * rt_inv ;
	k[rCP24] = -150620000 *rt_inv;
	k[rCP28] = 7.6 * lgt + 14770000 * rt_inv ;
	k[rCP30] = 1.02 * lgt + 8510000 * rt_inv ;
	k[rCP31] = -8 * lgt - 454700000 * rt_inv ;
	k[rCP32f] = 1.46 * lgt - 5670000 * rt_inv ;
	k[rCP32b] = 0.796263 * lgt - 220916454.1 * rt_inv ;
	k[rCP33f] = -8370000 *rt_inv;
	k[rCP33b] = 1.12696 * lgt - 406300873.9 * rt_inv ;
	k[rCP34] = 3.07 * lgt - 23970000 * rt_inv ;
	k[rI00b] = -2.03247 * lgt - 521300387.2 * rt_inv ;
	k[rI01f] = -15.16 * lgt - 486900000 * rt_inv ;
	k[rI01b] = -14.9678 * lgt - 147532832.2 * rt_inv ;
	k[rI02f] = -9450000 *rt_inv;
	k[rI02b] = 0.540969 * lgt - 103156195.8 * rt_inv ;
	k[rI03f] = 2.48 * lgt - 46280000 * rt_inv ;
	k[rI03b] = 1.36023 * lgt - 103605736.2 * rt_inv ;
	k[rI05f] = 2.71 * lgt - 4630000 * rt_inv ;
	k[rI05b] = 3.21593 * lgt - 92570105.57 * rt_inv ;
	k[rI06f] =  2 * lgt;
	k[rI06b] = 2.44735 * lgt - 156414500.2 * rt_inv ;
	k[rI07f] = -155440000 *rt_inv;
	k[rI07b] = 0.179262 * lgt - 20736180.96 * rt_inv ;
	k[rI08f] = 2.6 * lgt - 53970000 * rt_inv ;
	k[rI08b] = 2.38999 * lgt - 81341853.07 * rt_inv ;
	k[rI09f] =  4 * lgt;
	k[rI09b] = 4.0222 * lgt - 102619350.6 * rt_inv ;
	k[rI12] = 0.88 * lgt - 4770000 * rt_inv ;
	k[rI15] = -7.8 * lgt - 29540000 * rt_inv ;
	k[rI17] = -4.03 * lgt - 147300000 * rt_inv ;
	k[rI18] = -4.85 * lgt - 103650000 * rt_inv ;
	k[rI19b] = -0.525451 * lgt - 246670384.8 * rt_inv ;
	k[rI20f] = 1.3 * lgt - 73920000 * rt_inv ;
	k[rI20b] = 1.20951 * lgt - 250029253.6 * rt_inv ;
	k[rI23] = -7.74 * lgt - 99800000 * rt_inv ;
	k[rI31] = -8 * lgt - 454700000 * rt_inv ;
	k[rI32] = 1.46 * lgt - 5670000 * rt_inv ;
	k[rT01f] = 2.17 * lgt - 17420000 * rt_inv ;
	k[rT01b] = 3.57723 * lgt - 50936358.17 * rt_inv ;
	k[rT02f] = 0.68 * lgt - 373240000 * rt_inv ;
	k[rT02b] = 0.254459 * lgt + 1763440.913 * rt_inv ;
	k[rT03f] = -1.73 * lgt - 436010000 * rt_inv ;
	k[rT03b] = -0.111273 * lgt + 6505423.169 * rt_inv ;
	k[rT04f] = -14.15 * lgt - 285890000 * rt_inv ;
	k[rT04b] = -12.1057 * lgt - 218378017.7 * rt_inv ;
	k[rT05f] = -337550000 *rt_inv;
	k[rT05b] = 2.2454 * lgt - 46091406.02 * rt_inv ;
	k[rT06f] = 2.5 * lgt - 192650000 * rt_inv ;
	k[rT06b] = 2.06155 * lgt - 22309907.82 * rt_inv ;
	k[rT07f] = 3.98 * lgt - 14160000 * rt_inv ;
	k[rT07b] = 3.90326 * lgt - 72229922.61 * rt_inv ;
	k[rT08f] = 2.39 * lgt + 2520000 * rt_inv ;
	k[rT08b] = 2.21964 * lgt - 118258227 * rt_inv ;
	k[rT09f] = 2.88 * lgt - 13480000 * rt_inv ;
	k[rT09b] = 2.76966 * lgt - 48452858.8 * rt_inv ;
	k[rT10f] = 3.37 * lgt - 19750000 * rt_inv ;
	k[rT10b] = 1.86503 * lgt - 21643114.7 * rt_inv ;
	k[rT11f] = 4.09 * lgt - 10650000 * rt_inv ;
	k[rT11b] = 3.97822 * lgt - 62953832.43 * rt_inv ;
	k[rT12f] = 0.3 * lgt - 18420000 * rt_inv ;
	k[rT12b] = -0.891202 * lgt - 447620388 * rt_inv ;
	k[rT13f] = 1.8 * lgt - 16630000 * rt_inv ;
	k[rT13b] = 0.787767 * lgt - 79387362.49 * rt_inv ;
	k[rT14f] = -93120000 *rt_inv;
	k[rT14b] = -0.595515 * lgt - 160103077.4 * rt_inv ;
	k[rT15f] = 2.5 * lgt - 61440000 * rt_inv ;
	k[rT15b] = 1.67228 * lgt - 53175579.93 * rt_inv ;
	k[rT16f] = -50000000 *rt_inv;
	k[rT16b] = -0.637038 * lgt - 151028340.4 * rt_inv ;
	k[rT17b] = -0.6588 * lgt - 335402917 * rt_inv ;
	k[rT18b] = -0.492206 * lgt - 341422792.2 * rt_inv ;
	k[rT19f] = 1.03 * lgt + 9410000 * rt_inv ;
	k[rT19b] = 0.819065 * lgt - 50768437.11 * rt_inv ;
	k[rT20] = -7.74 * lgt - 99800000 * rt_inv ;
	k[rT21f] = -1.55 * lgt - 47370000 * rt_inv ;
	k[rT21b] = -1.51645 * lgt - 253641155 * rt_inv ;
	k[rT22f] = -8.86 * lgt - 69380000 * rt_inv ;
	k[rT22b] = -8.19414 * lgt - 329708198.2 * rt_inv ;
	k[rT24f] = 2.1 * lgt - 20380000 * rt_inv ;
	k[rT24b] = 1.96844 * lgt - 20126215 * rt_inv ;
	k[rT25f] = 2.5 * lgt - 20920000 * rt_inv ;
	k[rT25b] = 2.33341 * lgt - 14900124.81 * rt_inv ;
	k[rT26f] = 2 * lgt - 6280000 * rt_inv ;
	k[rT26b] = 1.77482 * lgt - 68734519.42 * rt_inv ;
	k[rT27f] = 1.5 * lgt - 41590000 * rt_inv ;
	k[rT27b] = 0.849667 * lgt - 50249369.82 * rt_inv ;
	k[rT28f] = 2.19 * lgt - 160000 * rt_inv ;
	k[rT28b] = 2.07772 * lgt - 106603252.2 * rt_inv ;
	k[rT29f] = -5.08 * lgt - 93090000 * rt_inv ;
	k[rT29b] = -4.82261 * lgt - 34519369.16 * rt_inv ;
	k[rT30f] = -6.21 * lgt - 154180000 * rt_inv ;
	k[rT30b] = -3.78987 * lgt - 40040180.59 * rt_inv ;
	k[rT31f] = -6.1 * lgt - 120540000 * rt_inv ;
	k[rT31b] = -3.79215 * lgt - 112843432.7 * rt_inv ;
	k[rT32b] = 0.606189 * lgt - 374502732.7 * rt_inv ;
	k[rT33b] = 0.571154 * lgt - 368736642.5 * rt_inv ;
	k[rT34b] = 0.512565 * lgt - 437211037.1 * rt_inv ;
	k[rT35f] = 7.6 * lgt + 14770000 * rt_inv ;
	k[rT35b] = 7.84448 * lgt - 131322717.9 * rt_inv ;
	k[rT36] = -342000000 *rt_inv;
	k[rT37] = 1.16 * lgt - 10060000 * rt_inv ;
	k[rT38] = -7570000 *rt_inv;
	k[rT39] = 1.35 * lgt + 6580000 * rt_inv ;
	k[rT40] = 2.5 * lgt - 157130000 * rt_inv ;
	k[rT41] = 2.5 * lgt - 42690000 * rt_inv ;
	k[rT42] = 6.21 * lgt - 6820000 * rt_inv ;
	k[rT43f] = -15.92 * lgt - 522120000 * rt_inv ;
	k[rT43b] = -15.741 * lgt - 155676974.5 * rt_inv ;
	k[rT44f] = -51870000 *rt_inv;
	k[rT44b] = 0.527769 * lgt - 118500338 * rt_inv ;
	k[rT45f] = -30760000 *rt_inv;
	k[rT45b] = 0.492735 * lgt - 91624247.79 * rt_inv ;
	k[rT46f] = 2 * lgt + 5490000 * rt_inv ;
	k[rT46b] = 2.43415 * lgt - 123848642.4 * rt_inv ;
	k[rT47] = -152400000 *rt_inv;
	k[rT50f] = -12.48 * lgt - 619590000 * rt_inv ;
	k[rT50b] = -12.2592 * lgt - 143886133 * rt_inv ;
	k[rT51f] = 1.8 * lgt - 68420000 * rt_inv ;
	k[rT51b] = 2.36964 * lgt - 25789496.49 * rt_inv ;
	k[rT52f] = -61500000 *rt_inv;
	k[rT52b] = 0.53461 * lgt - 13103406.31 * rt_inv ;
	k[rT53f] = 2.69 * lgt - 2590000 * rt_inv ;
	k[rT53b] = 3.16602 * lgt - 22667800.92 * rt_inv ;
	k[rT54f] = 4.46 * lgt - 57060000 * rt_inv ;
	k[rT54b] = 4.51087 * lgt - 23342651.31 * rt_inv ;
	k[rT55b] = -1.23308 * lgt - 538461229.5 * rt_inv ;
	k[rT56b] = -1.54684 * lgt - 111153956.2 * rt_inv ;
	k[rT57b] = -0.785212 * lgt - 263236749.6 * rt_inv ;
	k[rT58f] = -2.27 * lgt - 30080000 * rt_inv ;
	k[rT58b] = -3.38189 * lgt - 70672825.02 * rt_inv ;
	k[rT59] = -0.44 * lgt + 6900000 * rt_inv ;
	k[rT60] = -0.44 * lgt + 6900000 * rt_inv ;
	k[rE01b] = -0.723732 * lgt - 430033992.4 * rt_inv ;
	k[rE02b] = -1.4418 * lgt - 329900840.3 * rt_inv ;
	k[rE03b] = -2.29587 * lgt - 434547961.2 * rt_inv ;
	k[rE04f] = 2.17 * lgt - 17420000 * rt_inv ;
	k[rE04b] = 4.25437 * lgt - 58903820.16 * rt_inv ;
	k[rE05f] = 2.88 * lgt - 13480000 * rt_inv ;
	k[rE05b] = 3.4468 * lgt - 56420320.79 * rt_inv ;
	k[rE06f] = 4.71 * lgt - 25990000 * rt_inv ;
	k[rE06b] = 5.78253 * lgt - 29029371.1 * rt_inv ;
	k[rE07f] = 4.09 * lgt - 10650000 * rt_inv ;
	k[rE07b] = 5.1275 * lgt - 7923280.916 * rt_inv ;
	k[rE08f] = 2.01 * lgt - 1530000 * rt_inv ;
	k[rE08b] = 2.98891 * lgt - 67277675.53 * rt_inv ;
	k[rE09f] = 2.6 * lgt - 58200000 * rt_inv ;
	k[rE09b] = 2.92155 * lgt + 5094971.583 * rt_inv ;
	k[rE10f] = 3.65 * lgt - 29930000 * rt_inv ;
	k[rE10b] = 4.20376 * lgt - 41882525.92 * rt_inv ;
	k[rE11f] = 0.78 * lgt - 161940000 * rt_inv ;
	k[rE11b] = 2.47345 * lgt - 6105276.279 * rt_inv ;
	k[rE12f] = 2.08 * lgt - 134330000 * rt_inv ;
	k[rE12b] = 1.84471 * lgt - 9686488.018 * rt_inv ;
	k[rE13b] = 0.113509 * lgt - 308429851.5 * rt_inv ;
	k[rE14b] = 0.0198847 * lgt - 371138156 * rt_inv ;
	k[rE15f] = 3220000 *rt_inv;
	k[rE15b] = -0.405265 * lgt - 314123006.4 * rt_inv ;
	k[rE16f] = -1.63 * lgt - 14300000 * rt_inv ;
	k[rE16b] = -1.8782 * lgt - 94319836.75 * rt_inv ;
	k[rE17f] = 0.03 * lgt + 1650000 * rt_inv ;
	k[rE17b] = 0.465129 * lgt - 387258232 * rt_inv ;
	k[rE18f] = 0.03 * lgt + 1650000 * rt_inv ;
	k[rE18b] = 0.346652 * lgt - 375315438.3 * rt_inv ;
	k[rE19] = 4180000 *rt_inv;
	k[rE30] = -8.23 * lgt - 21620000 * rt_inv ;
	k[rE31] = 1.84 * lgt + 2420000 * rt_inv ;
	k[rE32] = -10.15 * lgt - 170770000 * rt_inv ;
	k[rE33] = -4.48 * lgt - 136430000 * rt_inv ;
	k[rE34] = -1.12 * lgt - 113030000 * rt_inv ;
	k[rE35] = 3.57 * lgt - 67350000 * rt_inv ;
	k[rE37] = -180100000 *rt_inv;
	k[rST01] = -326900000 *rt_inv;
	k[rST02f] = 1.66 * lgt - 2750000 * rt_inv ;
	k[rST02b] = 2.42964 * lgt - 146713214.8 * rt_inv ;
	k[rST03f] = 2 * lgt - 38490000 * rt_inv ;
	k[rST03b] = 1.94272 * lgt - 14154671.17 * rt_inv ;
	k[rST04f] = -7.8 * lgt - 29540000 * rt_inv ;
	k[rST04b] = -7.56182 * lgt - 103841676.9 * rt_inv ;
	k[rST05f] = -7.8 * lgt - 29540000 * rt_inv ;
	k[rST05b] = -7.44334 * lgt - 115784470.6 * rt_inv ;
	k[rST06f] = 2.88 * lgt - 13480000 * rt_inv ;
	k[rST06b] = 3.2153 * lgt + 6725309.257 * rt_inv ;
	k[rST10f] = -61500000 *rt_inv;
	k[rST10b] = 0.533908 * lgt - 5104404.95 * rt_inv ;
	k[rST11f] = 0.21 * lgt + 1790000 * rt_inv ;
	k[rST11b] = 1.0994 * lgt - 542745477.2 * rt_inv ;
	k[rST12f] = 1.61 * lgt + 1610000 * rt_inv ;
	k[rST12b] = 1.89788 * lgt - 71644724.18 * rt_inv ;
	k[rST13] = 0.29 * lgt - 50000 * rt_inv ;
	k[rST14f] = -1.39 * lgt - 4250000 * rt_inv ;
	k[rST14b] = -1.01328 * lgt - 393498943.6 * rt_inv ;
	k[rST00b] = 0.771407 * lgt - 83837426.03 * rt_inv ;
	k[rXY00f] = -0.6 * lgt - 396590000 * rt_inv ;
	k[rXY00b] = -1.02266 * lgt - 21983285.59 * rt_inv ;
	k[rXY01f] = -3.58 * lgt - 460930000 * rt_inv ;
	k[rXY01b] = -1.95607 * lgt - 18728819.02 * rt_inv ;
	k[rXY02f] = 3.98 * lgt - 14160000 * rt_inv ;
	k[rXY02b] = 3.90614 * lgt - 72626649.11 * rt_inv ;
	k[rXY03f] = 4.09 * lgt - 10650000 * rt_inv ;
	k[rXY03b] = 3.98111 * lgt - 63350558.93 * rt_inv ;
	k[rXY04f] = 2.39 * lgt + 2520000 * rt_inv ;
	k[rXY04b] = 2.22252 * lgt - 118654953.5 * rt_inv ;
	k[rXY05f] = 2.5 * lgt - 192650000 * rt_inv ;
	k[rXY05b] = 2.06443 * lgt - 22706634.32 * rt_inv ;
	k[rXY06f] = 2.5 * lgt - 61440000 * rt_inv ;
	k[rXY06b] = 1.67516 * lgt - 53572306.43 * rt_inv ;
	k[rXY07f] = -93120000 *rt_inv;
	k[rXY07b] = -0.592634 * lgt - 160499803.9 * rt_inv ;
	k[rXY09f] = 2.17 * lgt - 17420000 * rt_inv ;
	k[rXY09b] = 3.57309 * lgt - 50922686.05 * rt_inv ;
	k[rXY10f] = 2.88 * lgt - 13480000 * rt_inv ;
	k[rXY10b] = 2.77812 * lgt - 48875800.75 * rt_inv ;
	k[rXY11] = 1.55 * lgt - 12930000 * rt_inv ;
	k[rXY12] = -337550000 *rt_inv;
	k[rXY13f] = -14.15 * lgt - 285890000 * rt_inv ;
	k[rXY13b] = -12.1034 * lgt - 218295533.4 * rt_inv ;
	k[rXY14f] = -1.34 * lgt - 6660000 * rt_inv ;
	k[rXY14b] = -1.74104 * lgt - 283859870.9 * rt_inv ;
	k[rXY15f] = -2.47 * lgt - 67750000 * rt_inv ;
	k[rXY15b] = -0.706342 * lgt - 288930613.3 * rt_inv ;
	k[rXY16f] = -2.36 * lgt - 34110000 * rt_inv ;
	k[rXY16b] = -0.717967 * lgt - 361405951.1 * rt_inv ;
	k[rXY18f] = 2.42 * lgt - 31130000 * rt_inv ;
	k[rXY18b] = 2.45392 * lgt - 237768739.7 * rt_inv ;
	k[rXY19f] = 1.87 * lgt - 20930000 * rt_inv ;
	k[rXY19b] = 2.52177 * lgt - 282703438.3 * rt_inv ;
	k[rXY201] = -0.31 * lgt + 2750000 * rt_inv ;
	k[rXY202] = -1.44 * lgt - 58340000 * rt_inv ;
	k[rXY203] = -1.33 * lgt - 24700000 * rt_inv ;
	k[rXY22] = -7.74 * lgt - 99800000 * rt_inv ;
	k[rXY23f] = -0.6 * lgt - 396590000 * rt_inv ;
	k[rXY23b] = -1.02766 * lgt - 19452677.75 * rt_inv ;
	k[rXY24] = -3.58 * lgt - 460930000 * rt_inv ;
	k[rXY25] = -342000000 *rt_inv;
	k[rXY26f] = 3.98 * lgt - 14160000 * rt_inv ;
	k[rXY26b] = 3.90114 * lgt - 70096041.28 * rt_inv ;
	k[rXY27f] = 4.09 * lgt - 10650000 * rt_inv ;
	k[rXY27b] = 3.9761 * lgt - 60819951.09 * rt_inv ;
	k[rXY28f] = 2.39 * lgt + 2520000 * rt_inv ;
	k[rXY28b] = 2.21751 * lgt - 116124345.7 * rt_inv ;
	k[rXY29f] = 2.5 * lgt - 192650000 * rt_inv ;
	k[rXY29b] = 2.05943 * lgt - 20176026.48 * rt_inv ;
	k[rXY30f] = 2.5 * lgt - 61440000 * rt_inv ;
	k[rXY30b] = 1.67015 * lgt - 51041698.59 * rt_inv ;
	k[rXY31f] = -93120000 *rt_inv;
	k[rXY31b] = -0.597639 * lgt - 157969196.1 * rt_inv ;
	k[rXY33] = 1.16 * lgt - 10060000 * rt_inv ;
	k[rXY34] = -7570000 *rt_inv;
	k[rXY35] = 1.35 * lgt + 6580000 * rt_inv ;
	k[rXY36] = 2.5 * lgt - 157130000 * rt_inv ;
	k[rXY37] = 2.5 * lgt - 42690000 * rt_inv ;
	k[rXY38] = 6.21 * lgt - 6820000 * rt_inv ;
	k[rXY39f] = 2.17 * lgt - 17420000 * rt_inv ;
	k[rXY39b] = 4.21307 * lgt - 67516080.26 * rt_inv ;
	k[rXY40f] = 2.17 * lgt - 17420000 * rt_inv ;
	k[rXY40b] = 3.56984 * lgt - 50158374.84 * rt_inv ;
	k[rXY41f] = 2.88 * lgt - 13480000 * rt_inv ;
	k[rXY41b] = 3.41811 * lgt - 65469194.96 * rt_inv ;
	k[rXY42] = 2.88 * lgt - 13480000 * rt_inv ;
	k[rXY43f] = -1.34 * lgt - 6660000 * rt_inv ;
	k[rXY43b] = -1.77587 * lgt - 280931537.4 * rt_inv ;
	k[rXY44] = -2.47 * lgt - 67750000 * rt_inv ;
	k[rXY45f] = -2.36 * lgt - 34110000 * rt_inv ;
	k[rXY45b] = -0.716213 * lgt - 363172247.8 * rt_inv ;
	k[rXY47f] = 2.42 * lgt - 31980000 * rt_inv ;
	k[rXY47b] = 2.41908 * lgt - 235690406.2 * rt_inv ;
	k[rXY48] = 1.03 * lgt + 9410000 * rt_inv ;
	k[rXY50] = -342000000 *rt_inv;
	k[rXY51] = 1.16 * lgt - 10060000 * rt_inv ;
	k[rXY52] = -7570000 *rt_inv;
	k[rXY53] = 1.35 * lgt + 6580000 * rt_inv ;
	k[rXY54] = 2.5 * lgt - 157130000 * rt_inv ;
	k[rXY55] = 2.5 * lgt - 42690000 * rt_inv ;
	k[rXY56] = 6.21 * lgt - 6820000 * rt_inv ;
	k[rXY57f] = 2.17 * lgt - 17420000 * rt_inv ;
	k[rXY57b] = 4.24966 * lgt - 72210710.39 * rt_inv ;
	k[rXY58] = 2.88 * lgt - 13480000 * rt_inv ;
	k[rN01f] = 2.17 * lgt - 17420000 * rt_inv ;
	k[rN01b] = 4.14699 * lgt - 50392663.12 * rt_inv ;
	k[rN02f] = 2.88 * lgt - 13480000 * rt_inv ;
	k[rN02b] = 3.74052 * lgt - 49438848.61 * rt_inv ;
	k[rN04f] = -0.6 * lgt - 396590000 * rt_inv ;
	k[rN04b] = -0.435322 * lgt - 9868756.914 * rt_inv ;
	k[rN05f] = -5.02 * lgt - 478030000 * rt_inv ;
	k[rN05b] = -2.80977 * lgt - 12814647.68 * rt_inv ;
	k[rN06f] = -14.15 * lgt - 285890000 * rt_inv ;
	k[rN06b] = -12.1044 * lgt - 207395890.8 * rt_inv ;
	k[rN07f] = -337550000 *rt_inv;
	k[rN07b] = 1.84259 * lgt - 103227017.7 * rt_inv ;
	k[rN08f] = 3.98 * lgt - 14160000 * rt_inv ;
	k[rN08b] = 4.49348 * lgt - 60512120.44 * rt_inv ;
	k[rN09f] = 4.09 * lgt - 10650000 * rt_inv ;
	k[rN09b] = 4.56844 * lgt - 51236030.25 * rt_inv ;
	k[rN10f] = 2.39 * lgt + 2520000 * rt_inv ;
	k[rN10b] = 2.80985 * lgt - 106540424.9 * rt_inv ;
	k[rN11f] = 2.5 * lgt - 192650000 * rt_inv ;
	k[rN11b] = 2.65177 * lgt - 10592105.65 * rt_inv ;
	k[rN12f] = -93120000 *rt_inv;
	k[rN12b] = -0.00529554 * lgt - 148385275.3 * rt_inv ;
	k[rN13f] = 2.5 * lgt - 61440000 * rt_inv ;
	k[rN13b] = 2.2625 * lgt - 41457777.75 * rt_inv ;
	k[rN14] = -18960000 *rt_inv;
	k[rN15] = -18960000 *rt_inv;
	k[rN16b] = -0.655882 * lgt - 318374553.9 * rt_inv ;
	k[rN18f] = 1.03 * lgt + 9410000 * rt_inv ;
	k[rN18b] = 0.821983 * lgt - 33740074.02 * rt_inv ;
	k[rN19] = -7.74 * lgt - 99800000 * rt_inv ;
	k[rN20f] = -1.55 * lgt - 47370000 * rt_inv ;
	k[rN20b] = -1.93741 * lgt - 261565117.8 * rt_inv ;
	k[rN21f] = 0.37 * lgt - 70750000 * rt_inv ;
	k[rN21b] = 1.33906 * lgt - 347371636.2 * rt_inv ;
	k[rN22f] = -5.08 * lgt - 93090000 * rt_inv ;
	k[rN22b] = -5.24648 * lgt - 59471695 * rt_inv ;
	k[rN23f] = -6.21 * lgt - 154180000 * rt_inv ;
	k[rN23b] = -3.7915 * lgt - 46086416.7 * rt_inv ;
	k[rN24b] = 0.182318 * lgt - 399455058.5 * rt_inv ;
	k[rN25b] = 0.147283 * lgt - 393688968.3 * rt_inv ;
	k[rN26b] = 0.0886938 * lgt - 462163363 * rt_inv ;
	k[rN27f] = 7.6 * lgt + 14770000 * rt_inv ;
	k[rN27b] = 7.42061 * lgt - 156275043.7 * rt_inv ;
	k[rN28] = -342000000 *rt_inv;
	k[rN29] = 1.16 * lgt - 10060000 * rt_inv ;
	k[rN30] = -7570000 *rt_inv;
	k[rN31] = 1.35 * lgt + 6580000 * rt_inv ;
	k[rN32] = 2.5 * lgt - 157130000 * rt_inv ;
	k[rN33] = 2.5 * lgt - 42690000 * rt_inv ;
	k[rN34] = 6.21 * lgt - 6820000 * rt_inv ;
	k[rOX00f] = -12.48 * lgt - 619590000 * rt_inv ;
	k[rOX00b] = -12.2685 * lgt - 143558218.7 * rt_inv ;
	k[rOX02f] = -19.31 * lgt - 284180000 * rt_inv ;
	k[rOX02b] = -20.3484 * lgt - 537844017.3 * rt_inv ;
	k[rOX03f] = 1.8 * lgt - 68420000 * rt_inv ;
	k[rOX03b] = 2.3603 * lgt - 25461582.19 * rt_inv ;
	k[rOX04f] = 2.68 * lgt - 3070000 * rt_inv ;
	k[rOX04b] = 3.14667 * lgt - 22819886.61 * rt_inv ;
	k[rOX05f] = 3.25 * lgt - 23390000 * rt_inv ;
	k[rOX05b] = 1.73243 * lgt - 24846500.63 * rt_inv ;
	k[rOX06f] = -269000000 *rt_inv;
	k[rOX06b] = 0.19859 * lgt + 2368432.605 * rt_inv ;
	k[rOX99f] = 1.55 * lgt - 12930000 * rt_inv ;
	k[rOX99b] = 0.54483 * lgt - 73831723.77 * rt_inv ;
	k[rOX07f] = 1.8 * lgt - 16630000 * rt_inv ;
	k[rOX07b] = 0.79483 * lgt - 77531723.77 * rt_inv ;
	k[rOX08f] = 0.3 * lgt - 18420000 * rt_inv ;
	k[rOX08b] = -0.903806 * lgt - 447183774 * rt_inv ;
	k[rOX09f] = -0.77 * lgt - 63980000 * rt_inv ;
	k[rOX09b] = 0.100471 * lgt - 368646648 * rt_inv ;
	k[rOX10f] = -61500000 *rt_inv;
	k[rOX10b] = 0.525263 * lgt - 12775492 * rt_inv ;
	k[rOX11f] = -2.27 * lgt - 30080000 * rt_inv ;
	k[rOX11b] = -3.36548 * lgt - 69145100.61 * rt_inv ;
	k[rOX12f] = -37580000 *rt_inv;
	k[rOX12b] = -1.69733 * lgt - 249952673.3 * rt_inv ;
	k[rOX13b] = -1.21667 * lgt - 536933505.1 * rt_inv ;
	k[rOX14b] = -1.53043 * lgt - 109626231.8 * rt_inv ;
	k[rOX15b] = -0.768803 * lgt - 261709025.2 * rt_inv ;
	k[rOX16f] = 4.57 * lgt - 22000000 * rt_inv ;
	k[rOX16b] = 4.52848 * lgt - 56045263 * rt_inv ;
	k[rOX17f] = -0.61 * lgt - 310110000 * rt_inv ;
	k[rOX17b] = 1.46428 * lgt - 186012874.1 * rt_inv ;
	k[rOX18f] = -15.92 * lgt - 522120000 * rt_inv ;
	k[rOX18b] = -15.7214 * lgt - 154257949.8 * rt_inv ;
	k[rOX19f] = 3.4 * lgt - 30260000 * rt_inv ;
	k[rOX19b] = 3.94744 * lgt - 95471313.33 * rt_inv ;
	k[rOX20f] = 3.4 * lgt + 4780000 * rt_inv ;
	k[rOX20b] = 3.85381 * lgt - 123139617.8 * rt_inv ;
	k[rOX21f] = 4.7 * lgt - 20200000 * rt_inv ;
	k[rOX21b] = 4.72866 * lgt - 94324468.15 * rt_inv ;
	k[rOX22f] = -173220000 *rt_inv;
	k[rOX22b] = 0.185728 * lgt - 10021298.54 * rt_inv ;
	k[rOX23f] = -152400000 *rt_inv;
	k[rOX23b] = 2.57168 * lgt - 46699482.58 * rt_inv ;
	k[rOX24f] = 0.03 * lgt + 1650000 * rt_inv ;
	k[rOX24b] = -0.571852 * lgt - 171657572.7 * rt_inv ;
	k[rOX25f] = 1.3 * lgt - 73920000 * rt_inv ;
	k[rOX25b] = 1.13311 * lgt - 176666441.5 * rt_inv ;
	k[rOX26f] = -246860000 *rt_inv;
	k[rOX26b] = 2.13692 * lgt - 208819519.7 * rt_inv ;
	k[rOX27f] = 1.45 * lgt - 16180000 * rt_inv ;
	k[rOX27b] = 2.92318 * lgt - 193385973.8 * rt_inv ;
	k[rOX28] = -20920000 *rt_inv;
	k[rOX30f] = -18950000 *rt_inv;
	k[rOX30b] = -0.681513 * lgt - 84971054.68 * rt_inv ;
	k[rOX31f] = -18950000 *rt_inv;
	k[rOX31b] = -0.802706 * lgt - 449243458.8 * rt_inv ;
	k[rOX32f] = 3.25 * lgt - 23390000 * rt_inv ;
	k[rOX32b] = 2.13353 * lgt - 26376185.49 * rt_inv ;
	k[rOX33f] = -2.27 * lgt - 30080000 * rt_inv ;
	k[rOX33b] = -3.06357 * lgt - 96420665.62 * rt_inv ;
	k[rOX34f] = -2.27 * lgt - 30080000 * rt_inv ;
	k[rOX34b] = -3.04777 * lgt - 95129408.91 * rt_inv ;
	k[rOX35] = -37580000 *rt_inv;
	k[rOX36] = -37580000 *rt_inv;
	k[rOX37b] = -0.914757 * lgt - 564209070.1 * rt_inv ;
	k[rOX38b] = -0.898966 * lgt - 562917813.4 * rt_inv ;
	k[rOX39b] = -1.22852 * lgt - 136901796.8 * rt_inv ;
	k[rOX40b] = -1.21273 * lgt - 135610540.1 * rt_inv ;
	k[rOX41f] = -0.61 * lgt - 310110000 * rt_inv ;
	k[rOX41b] = 1.18469 * lgt - 220346419.8 * rt_inv ;
	k[rOX42f] = -15.92 * lgt - 522120000 * rt_inv ;
	k[rOX42b] = -15.7988 * lgt - 157847595.9 * rt_inv ;
	k[rOX43f] = 3.4 * lgt - 30260000 * rt_inv ;
	k[rOX43b] = 3.86999 * lgt - 99060959.38 * rt_inv ;
	k[rOX44f] = 3.4 * lgt + 4780000 * rt_inv ;
	k[rOX44b] = 3.77637 * lgt - 126729263.8 * rt_inv ;
	k[rOX45f] = 4.7 * lgt - 20200000 * rt_inv ;
	k[rOX45b] = 4.65122 * lgt - 97914114.2 * rt_inv ;
	k[rOX46f] = -152400000 *rt_inv;
	k[rOX46b] = 1.86567 * lgt - 87541656.19 * rt_inv ;
	k[rOX48] = 1.3 * lgt - 73920000 * rt_inv ;
	k[rOX50] = -2.27 * lgt - 30080000 * rt_inv ;
	k[rOX51] = -2.27 * lgt - 30080000 * rt_inv ;
	k[rOX52] = -2.27 * lgt - 30080000 * rt_inv ;
	k[rOX53] = -2.27 * lgt - 30080000 * rt_inv ;
	k[rOX54] = -2.27 * lgt - 30080000 * rt_inv ;
	k[rOX60] = 3.25 * lgt - 23390000 * rt_inv ;
	k[rOX61] = 3.25 * lgt - 23390000 * rt_inv ;
	k[rOX62] = 3.25 * lgt - 23390000 * rt_inv ;
	k[rOX63] = 3.25 * lgt - 23390000 * rt_inv ;
	k[rOX64] = 3.25 * lgt - 23390000 * rt_inv ;
	k[rOX65] = 3.25 * lgt - 23390000 * rt_inv ;
	k[rOX66] = 3.25 * lgt - 23390000 * rt_inv ;
	k[rNO178f] = -1485460 *rt_inv;
	k[rNO178b] = 0.125609 * lgt - 315655094.7 * rt_inv ;
	k[rNO179f] = 1 * lgt - 27198600 * rt_inv ;
	k[rNO179b] = 1.09006 * lgt - 158871938.8 * rt_inv ;
	k[rNO180f] = -1611000 *rt_inv;
	k[rNO180b] = -0.344901 * lgt - 203845469.9 * rt_inv ;
	k[rNO181f] = -45233400 *rt_inv;
	k[rNO181b] = 1.001 * lgt - 372613989.4 * rt_inv ;
	k[rNO182f] = -96868900 *rt_inv;
	k[rNO182b] = 0.965452 * lgt - 241753193.4 * rt_inv ;
	k[rNO183f] = -79001500 *rt_inv;
	k[rNO183b] = 1.43596 * lgt - 335820958.2 * rt_inv ;
	k[rNO184f] = -88123500 *rt_inv;
	k[rNO184b] = 0.674332 * lgt - 192860164.8 * rt_inv ;
	k[rNO186f] = 2008510 *rt_inv;
	k[rNO186b] = -0.254945 * lgt - 30698139.73 * rt_inv ;
	k[rNO187f] =  -1.41 * lgt;
	k[rNO187b] = -2.11281 * lgt - 307931129.6 * rt_inv ;
	k[rNO188f] = 1004260 *rt_inv;
	k[rNO188b] = 0.581618 * lgt - 188933014.9 * rt_inv ;
	k[rNO189f] = -1506390 *rt_inv;
	k[rNO189b] = 1.01658 * lgt - 120882533.7 * rt_inv ;
	k[rNO190b] = -0.269423 * lgt - 298102553.3 * rt_inv ;
	k[rNO191f] = -1380850 *rt_inv;
	k[rNO191b] = 0.110513 * lgt - 103015023.5 * rt_inv ;
	k[rNO192b] = -0.734627 * lgt - 74318505.35 * rt_inv ;
	k[rNO193f] =  1.2 * lgt;
	k[rNO193b] = 1.21689 * lgt - 164342477.9 * rt_inv ;
	k[rNO194f] = 2 * lgt - 27198600 * rt_inv ;
	k[rNO194b] = 1.70033 * lgt - 30955974.19 * rt_inv ;
	k[rNO195f] = 1.5 * lgt - 418440 * rt_inv ;
	k[rNO195b] = 1.66553 * lgt - 227959862.1 * rt_inv ;
	k[rNO196b] = -0.143813 * lgt - 612272188 * rt_inv ;
	k[rNO197f] = -57954000 *rt_inv;
	k[rNO197b] = -0.641003 * lgt - 69564200.92 * rt_inv ;
	k[rNO198f] =  -0.23 * lgt;
	k[rNO198b] = -0.0289121 * lgt - 410037718 * rt_inv ;
	k[rNO199f] =  -0.45 * lgt;
	k[rNO199b] = -1.68487 * lgt - 153218259.8 * rt_inv ;
	k[rNO200b] = 0.248665 * lgt - 45064379.46 * rt_inv ;
	k[rNO201b] = -0.485962 * lgt - 119382884.8 * rt_inv ;
	k[rNO202f] = -15273100 *rt_inv;
	k[rNO202b] = 0.2837 * lgt - 66103569.64 * rt_inv ;
	k[rNO203f] = 1.5 * lgt + 1924830 * rt_inv ;
	k[rNO203b] = 1.69008 * lgt - 111613944.1 * rt_inv ;
	k[rNO204b] = 0.283778 * lgt - 31185362.18 * rt_inv ;
	k[rNO205f] = -0.11 * lgt - 20838300 * rt_inv ;
	k[rNO205b] = 0.173778 * lgt - 52023662.18 * rt_inv ;
	k[rNO206b] = 0.270871 * lgt - 235848710.9 * rt_inv ;
	k[rNO207b] = 0.597544 * lgt - 458492635.5 * rt_inv ;
	k[rNO208b] = 0.396456 * lgt - 48454917.49 * rt_inv ;
	k[rNO209b] = 0.632578 * lgt - 464258725.7 * rt_inv ;
	k[rNO210b] = 0.538954 * lgt - 526967030.1 * rt_inv ;
	k[rNO211b] = 0.113805 * lgt - 473171880.5 * rt_inv ;
	k[rNO212f] = -1.32 * lgt - 3096460 * rt_inv ;
	k[rNO212b] = -1.47144 * lgt - 206619685.4 * rt_inv ;
	k[rNO213b] = 0.465204 * lgt - 223784047.9 * rt_inv ;
	k[rNO214f] = 0.72 * lgt - 2761710 * rt_inv ;
	k[rNO214b] = 1.22024 * lgt - 232311848.1 * rt_inv ;
	k[rNO215f] = 1.9 * lgt + 3975180 * rt_inv ;
	k[rNO215b] = 2.30661 * lgt - 288283262.5 * rt_inv ;
	k[rNO216f] = -54397200 *rt_inv;
	k[rNO216b] = 0.138531 * lgt - 55537323.31 * rt_inv ;
	k[rNO217b] = 0.0458286 * lgt - 325738498.6 * rt_inv ;
	k[rNO218b] = -1.42679 * lgt - 131288920.2 * rt_inv ;
	k[rNO219f] = -31215700 *rt_inv;
	k[rNO219b] = -0.745109 * lgt - 62067118.2 * rt_inv ;
	k[rNO220f] = 1841140 *rt_inv;
	k[rNO220b] = -0.991831 * lgt - 58886649.06 * rt_inv ;
	k[rNO221f] = 2.45 * lgt - 9373060 * rt_inv ;
	k[rNO221b] = 1.61127 * lgt - 102932782.6 * rt_inv ;
	k[rNO222b] = 1.12772 * lgt - 396684048.3 * rt_inv ;
	k[rNO223b] = 1.39714 * lgt - 98581495.08 * rt_inv ;
	k[rNO224b] = 0.813951 * lgt + 30623225 * rt_inv ;
	k[rNO225b] = 1.25333 * lgt - 710853683.1 * rt_inv ;
	k[rNO226f] = -83688100 *rt_inv;
	k[rNO226b] = 0.283477 * lgt - 517280772.3 * rt_inv ;
	k[rNO227f] = -226167000 *rt_inv;
	k[rNO227b] = 1.15885 * lgt + 6690694.945 * rt_inv ;
	k[rNO228f] = -1.52 * lgt - 3096460 * rt_inv ;
	k[rNO228b] = -1.35774 * lgt - 254896214.9 * rt_inv ;
	k[rNO229f] = -2 * lgt - 3347520 * rt_inv ;
	k[rNO229b] = -1.68097 * lgt - 619436488.3 * rt_inv ;
	k[rNO230f] = -3.3 * lgt - 529746000 * rt_inv ;
	k[rNO230b] = -2.81007 * lgt - 3112913.846 * rt_inv ;
	k[rNO231f] = 2.64 * lgt - 20838300 * rt_inv ;
	k[rNO231b] = 2.01691 * lgt - 52801407.41 * rt_inv ;
	k[rNO232f] = 2.64 * lgt - 20838300 * rt_inv ;
	k[rNO232b] = 3.41405 * lgt - 151382902.5 * rt_inv ;
	k[rNO233f] = 1.58 * lgt - 111305000 * rt_inv ;
	k[rNO233b] = 2.3837 * lgt - 11979187.19 * rt_inv ;
	k[rNO234f] = 2.03 * lgt - 55945500 * rt_inv ;
	k[rNO234b] = 1.3399 * lgt - 20408332.66 * rt_inv ;
	k[rNO235] = 2.26 * lgt - 26780200 * rt_inv ;
	k[rNO236f] = 2.56 * lgt - 37659600 * rt_inv ;
	k[rNO236b] = 3.08538 * lgt - 123139823 * rt_inv ;
	k[rNO238f] = -1673760 *rt_inv;
	k[rNO238b] = 0.86464 * lgt - 324705561 * rt_inv ;
	k[rNO239f] = -192566000 *rt_inv;
	k[rNO239b] = -0.0802483 * lgt + 1780628.844 * rt_inv ;
	k[rNO240f] = 0.88 * lgt - 84232000 * rt_inv ;
	k[rNO240b] = 0.194004 * lgt - 80840108.46 * rt_inv ;
	k[rNO242f] = -309646000 *rt_inv;
	k[rNO242b] = -0.505476 * lgt - 215937661.1 * rt_inv ;
	k[rNO243f] = -271986000 *rt_inv;
	k[rNO243b] = -0.567762 * lgt - 215580649.1 * rt_inv ;
	k[rNO244b] = 0.045361 * lgt - 119823005.9 * rt_inv ;
	k[rNO245b] = 0.0911896 * lgt - 445561504.5 * rt_inv ;
	k[rNO246b] = -0.560387 * lgt - 310777743.2 * rt_inv ;
	k[rNO247b] = -1.18348 * lgt - 342740850.6 * rt_inv ;
	k[rNO248b] = -0.257081 * lgt - 175632740.4 * rt_inv ;
	k[rNO249f] = -1.38 * lgt - 5314190 * rt_inv ;
	k[rNO249b] = -2.71396 * lgt - 393811209.4 * rt_inv ;
	k[rNO250f] = -0.69 * lgt - 3180150 * rt_inv ;
	k[rNO250b] = -0.994388 * lgt - 319509529.1 * rt_inv ;
	k[rNO251f] = -0.36 * lgt - 2426950 * rt_inv ;
	k[rNO251b] = -1.92712 * lgt - 102425371.3 * rt_inv ;
	k[rNO252f] = -1.38 * lgt - 5314190 * rt_inv ;
	k[rNO252b] = -2.77625 * lgt - 431114197.4 * rt_inv ;
	k[rNO253f] = -0.69 * lgt - 3180150 * rt_inv ;
	k[rNO253b] = -1.05667 * lgt - 356812517.1 * rt_inv ;
	k[rNO254f] = -0.36 * lgt - 2426950 * rt_inv ;
	k[rNO254b] = -1.98941 * lgt - 139728359.3 * rt_inv ;
	k[rNO255f] = -120511000 *rt_inv;
	k[rNO255b] = 0.113056 * lgt - 468389751.4 * rt_inv ;
	k[rNO256f] = -91010800 *rt_inv;
	k[rNO256b] = -0.262998 * lgt - 45223611.48 * rt_inv ;
	k[rNO257b] = 1.59876 * lgt - 600190077.9 * rt_inv ;
	k[rNO258b] = 0.937388 * lgt - 486915030.7 * rt_inv ;
	k[rNO259b] = 1.48749 * lgt - 168071258.2 * rt_inv ;
	k[rNO260b] = 1.05254 * lgt - 238632389.3 * rt_inv ;
	k[rNO261b] = 1.71229 * lgt - 282520816.4 * rt_inv ;
	k[rNO262f] = 1.41 * lgt - 35567400 * rt_inv ;
	k[rNO262b] = 1.93443 * lgt - 201414117.4 * rt_inv ;
	k[rNO263f] = 1.57 * lgt - 184114000 * rt_inv ;
	k[rNO263b] = 2.639 * lgt - 316809467.6 * rt_inv ;
	k[rNO264f] = 2.11 * lgt - 47702200 * rt_inv ;
	k[rNO264b] = 2.51648 * lgt - 7497667.142 * rt_inv ;
	k[rNO265f] = 1.7 * lgt - 15900700 * rt_inv ;
	k[rNO265b] = 3.25496 * lgt - 29213282.77 * rt_inv ;
	k[rNO266f] = 2.5 * lgt - 55652600 * rt_inv ;
	k[rNO266b] = 2.94152 * lgt - 21214157.33 * rt_inv ;
	k[rNO267f] = 1.5 * lgt - 15063900 * rt_inv ;
	k[rNO267b] = 1.8479 * lgt - 43333761.76 * rt_inv ;
	k[rNO268f] = 1.5 * lgt - 15063900 * rt_inv ;
	k[rNO268b] = 1.77576 * lgt - 135846237.9 * rt_inv ;
	k[rNO269f] = -354503000 *rt_inv;
	k[rNO269b] = 1.48986 * lgt + 14427311.12 * rt_inv ;
	k[rNO270f] = -0.69 * lgt - 11925600 * rt_inv ;
	k[rNO270b] = -0.456838 * lgt - 300424198.1 * rt_inv ;
	k[rNO271f] = 0.18 * lgt - 8870940 * rt_inv ;
	k[rNO271b] = 1.44274 * lgt - 225201897.8 * rt_inv ;
	k[rNO272f] = -0.75 * lgt - 12092900 * rt_inv ;
	k[rNO272b] = 1.03812 * lgt - 313904080.8 * rt_inv ;
	k[rNO273] = 2 * lgt - 8368810 * rt_inv ;
	k[rNO274b] = -0.0686479 * lgt - 209085327.2 * rt_inv ;
	k[rNO275f] = -0.31 * lgt - 1213480 * rt_inv ;
	k[rNO275b] = -0.917899 * lgt - 157660761.4 * rt_inv ;
	k[rNO276f] = 0.15 * lgt + 376596 * rt_inv ;
	k[rNO276b] = 0.0117785 * lgt - 487028321 * rt_inv ;
	k[rNO277f] = 2.4 * lgt - 41488400 * rt_inv ;
	k[rNO277b] = 2.76033 * lgt - 18446328.61 * rt_inv ;
	k[rNO278f] = 1.6 * lgt - 3996110 * rt_inv ;
	k[rNO278b] = 1.86671 * lgt - 43662343.04 * rt_inv ;
	k[rNO279f] = 1.94 * lgt - 27031200 * rt_inv ;
	k[rNO279b] = 2.2653 * lgt + 1776961.574 * rt_inv ;
	k[rNO280f] = -60046200 *rt_inv;
	k[rNO280b] = 0.54457 * lgt - 26894950.2 * rt_inv ;
	k[rNO281f] = -0.752 * lgt - 1443620 * rt_inv ;
	k[rNO281b] = -1.16221 * lgt - 252108683.9 * rt_inv ;
	k[rNO282f] = 2950000 *rt_inv;
	k[rNO282b] = -0.100357 * lgt - 475695653.8 * rt_inv ;
	k[rNO283f] = -47283800 *rt_inv;
	k[rNO283b] = 0.934296 * lgt - 142048514.8 * rt_inv ;
	k[r15f] = -1.72 * lgt - 2196000 * rt_inv ;
	w[r15f] =  0.44 * lgt;
	k[r15b] = -1.73291 * lgt - 206859348.7 * rt_inv ;
	w[r15b] = 0.427093 * lgt - 204663348.7 * rt_inv ;
	k[r22] = -2.3 * lgt - 203970000 * rt_inv ;
	w[r22] = 0.9 * lgt - 203966000 * rt_inv ;
	k[rG14] = -2.76 * lgt - 6690000 * rt_inv ;
	w[rG14] = 0.0;
	k[rLPG60f] = -4.89 * lgt - 14360000 * rt_inv ;
	w[rLPG60f] =  1.63 * lgt;
	k[rLPG60b] = -6.659 * lgt - 158341000 * rt_inv ;
	w[rLPG60b] = -0.139 * lgt - 143981000 * rt_inv ;
	k[r27] = -2.79 * lgt - 17540000 * rt_inv ;
	w[r27] = -9980000 *rt_inv;
	k[rLPG24] = -5.11 * lgt - 29690000 * rt_inv ;
	w[rLPG24] = 0.5 * lgt - 18870000 * rt_inv ;
	k[rLPG44] = -13.872 * lgt - 189916000 * rt_inv ;
	w[rLPG44] = -8.622 * lgt - 177536000 * rt_inv ;
	k[rG115f] = -4.66 * lgt - 15820000 * rt_inv ;
	w[rG115f] = 1.27 * lgt - 11330000 * rt_inv ;
	k[rG115b] = -4.48812 * lgt - 162604360.4 * rt_inv ;
	w[rG115b] = 1.44188 * lgt - 158114360.4 * rt_inv ;
	k[rG155f] = -6.64 * lgt - 24140000 * rt_inv ;
	w[rG155f] = 1.46 * lgt - 5670000 * rt_inv ;
	k[rG155b] = -6.76131 * lgt - 175460755 * rt_inv ;
	w[rG155b] = 1.33869 * lgt - 156990755 * rt_inv ;
	k[rLPDD76f] = 20.949 * lgt + 187371000 * rt_inv ;
	w[rLPDD76f] = 20.949 * lgt + 169799000 * rt_inv ;
	k[rLPDD76b] = 1.465 * lgt + 3246000 * rt_inv ;
	w[rLPDD76b] = 1.465 * lgt - 14326000 * rt_inv ;
	k[rG183] = -6.796 * lgt - 166004000 * rt_inv ;
	w[rG183] = -0.506 * lgt - 153584000 * rt_inv ;
	k[rG163] = -16.267 * lgt - 179210000 * rt_inv ;
	w[rG163] = -0.067 * lgt - 127040000 * rt_inv ;
	k[rT23f] = -7.44 * lgt - 58910000 * rt_inv ;
	w[rT23f] = 0.52 * lgt - 210000 * rt_inv ;
	k[rT23b] = -6.95964 * lgt - 492237148.5 * rt_inv ;
	w[rT23b] = 1.00036 * lgt - 433537148.5 * rt_inv ;
	k[rOX01f] = -18.87 * lgt - 376980000 * rt_inv ;
	w[rOX01f] = 0.62 * lgt - 323430000 * rt_inv ;
	k[rOX01b] = -18.4977 * lgt - 50925964.37 * rt_inv ;
	w[rOX01b] = 0.992264 * lgt + 2624035.634 * rt_inv ;
	k[rNO185f] = -237005000 *rt_inv;
	w[rNO185f] = -234410000 *rt_inv;
	k[rNO185b] = 1.1222 * lgt - 66517184.86 * rt_inv ;
	w[rNO185b] = 1.1222 * lgt - 63922184.86 * rt_inv ;
	k[rNO237f] = -3.4 * lgt - 7950370 * rt_inv ;
	w[rNO237f] = 0.0;
	k[rNO237b] = -3.52088 * lgt - 110066098 * rt_inv ;
	w[rNO237b] = -0.120877 * lgt - 102115728 * rt_inv ;
	k[rNO241f] = -3.16 * lgt - 3096460 * rt_inv ;
	w[rNO241f] =  0.15 * lgt;
	k[rNO241b] = -4.81452 * lgt - 142331281 * rt_inv ;
	w[rNO241b] = -1.50452 * lgt - 139234821 * rt_inv ;
	for (i=0;i<=rNO283b;i++) k[i] = exp(k[i]);
	for (i=r15f;i<=rNO241b;i++) {
		k[i] = exp(k[i]);
		w[i] = exp(w[i]);
	}
	k[r3f] *= 8.6310000000E+11;
	k[r3b] *= 5.0623635033E+11;
	k[r4f] *= 2.0660000000E+05;
	k[r4b] *= 2.1217331924E+06;
	k[r5f] *= 3.3200000000E+01;
	k[r5b] *= 5.8130803668E+02;
	k[r1f] *= 1.0740000000E+11;
	k[r1b] *= 2.1355916977E+08;
	k[r16] *= 2.8120000000E+03;
	k[r21] *= 1.4170000000E+08;
	k[r20] *= 9.4060000000E+11;
	k[r19] *= 8.2200000000E+08;
	k[r18] *= 1.7200000000E+07;
	k[r17] *= 6.2450000000E+10;
	k[r26] *= 3.7230000000E+08;
	k[r23] *= 2.4100000000E+10;
	k[r39] *= 1.9790000000E+10;
	k[rG06f] *= 1.0800000000E+11;
	k[rG06b] *= 1.7379827487E+12;
	k[rG08] *= 6.2810000000E+09;
	k[rG19] *= 9.9490000000E+09;
	k[rG20] *= 2.4000000000E+09;
	k[rG17] *= 1.1300000000E+04;
	k[rG21] *= 5.0000000000E+09;
	k[rG38b] *= 2.1872817488E+10;
	k[rG34b] *= 2.1163554877E+13;
	k[rG28f] *= 1.5000000000E+10;
	k[rG28b] *= 1.0936408744E+10;
	k[rG55] *= 5.6000000000E+04;
	k[rG57f] *= 7.6150000000E+14;
	k[rG57b] *= 2.5866625919E+13;
	k[rG59] *= 7.2770000000E+08;
	k[rG74f] *= 5.5100000000E+09;
	k[rG74b] *= 9.4829147064E+14;
	k[r29] *= 5.7150000000E+09;
	k[r31] *= 1.5700000000E+02;
	k[r28f] *= 6.9880000000E+01;
	k[r28b] *= 3.9069941448E+08;
	k[r38] *= 3.4280000000E+10;
	k[r36f] *= 4.3760000000E+08;
	k[r36b] *= 8.4837182791E+04;
	k[rDD68] *= 1.3830000000E+08;
	k[rG46] *= 3.9000000000E+10;
	k[rG47] *= 5.8080000000E+07;
	k[rG45] *= 5.6760000000E+04;
	k[rG49] *= 1.4000000000E+03;
	k[rLPDD62] *= 1.7910000000E+13;
	k[rLPG89] *= 5.9680000000E+12;
	k[rDD44] *= 5.3180000000E+11;
	k[rDD47] *= 1.8500000000E+04;
	k[rDD50] *= 3.6680000000E+03;
	k[rDD48] *= 1.0000000000E+09;
	k[rDD53] *= 3.4710000000E+22;
	k[rDD52] *= 2.6200000000E+03;
	k[rLPG63] *= 1.1520000000E+08;
	k[rLPDD59] *= 4.2510000000E-55;
	k[rLPG61] *= 1.5280000000E+10;
	k[rLPG64] *= 3.6340000000E+09;
	k[rDD66] *= 1.9900000000E+09;
	k[rDD67] *= 1.4540000000E+16;
	k[rG116] *= 8.1000000000E+03;
	k[rG117] *= 9.1850000000E+03;
	k[rLPG120f] *= 3.7780000000E-04;
	k[rLPG120b] *= 2.0570000000E+31;
	k[rG122] *= 1.2800000000E+06;
	k[rG135] *= 2.7550000000E+14;
	k[rLPG134] *= 2.5380000000E+08;
	k[rG133] *= 1.1790000000E+03;
	k[rLPG131] *= 1.0300000000E+10;
	k[rG160] *= 4.1840000000E-05;
	k[rG158] *= 7.1500000000E+01;
	k[rG159] *= 3.8900000000E+05;
	k[rLPG157] *= 7.3880000000E+06;
	k[rG156] *= 3.4820000000E+02;
	k[rG161] *= 3.5970000000E+33;
	k[rLPDD15] *= 4.3030000000E+09;
	k[rG171] *= 6.3900000000E+03;
	k[rG111] *= 6.0920000000E+07;
	k[rLPG124] *= 1.4570000000E+11;
	k[rLPR028] *= 8.7400000000E+12;
	k[rLPG123f] *= 4.8550000000E+15;
	k[rLPG123b] *= 3.8770000000E+07;
	k[rLPG127f] *= 4.1370000000E+15;
	k[rLPG127b] *= 3.3920000000E+08;
	k[rLPR025] *= 3.7050000000E+14;
	k[rLPR030] *= 4.6260000000E+11;
	k[rLPG139] *= 1.7080000000E+07;
	k[rLPG140] *= 1.5990000000E+07;
	k[rLPG137] *= 7.0040000000E+34;
	k[rLPG138] *= 3.0250000000E+10;
	k[rLPG136] *= 1.4210000000E+34;
	k[rLPG144] *= 3.2560000000E+10;
	k[rG168] *= 1.3200000000E+20;
	k[rLPR015] *= 1.1770000000E+26;
	k[rLPR116] *= 2.0330000000E+11;
	k[rLPR115] *= 1.1950000000E+11;
	k[rLPR133] *= 8.2140000000E-06;
	k[rLPLPR118] *= 2.8880000000E+11;
	k[rLPR016] *= 2.2850000000E+54;
	k[rLPLPR113] *= 4.5190000000E+12;
	k[rLPR111] *= 1.3780000000E+20;
	k[rLPR007] *= 1.1950000000E+34;
	k[rLPLPR117] *= 1.7170000000E+15;
	k[rR141] *= 4.3550000000E+18;
	k[rR143] *= 3.1760000000E+04;
	k[rLPR153] *= 8.8680000000E+01;
	k[rLPR142] *= 1.4940000000E+05;
	k[rLPR151] *= 4.4630000000E+00;
	k[rG184] *= 3.0820000000E+13;
	k[rLPHP55] *= 3.1670000000E+10;
	k[rLPHP43] *= 9.7530000000E+16;
	k[rLPHP54] *= 4.3310000000E+11;
	k[rLPHP52] *= 4.1170000000E+11;
	k[rLPHP53] *= 1.4860000000E+04;
	k[rLPDD514] *= 1.3160000000E+10;
	k[rLPHP56] *= 3.2380000000E+10;
	k[rLPHP39] *= 8.6100000000E+17;
	k[rHP41] *= 8.3900000000E+14;
	k[rHP29] *= 6.8490000000E+21;
	k[rHP32] *= 4.1220000000E-19;
	k[rDD199] *= 2.3670000000E+08;
	k[rDD212] *= 9.3780000000E+09;
	k[rLPDD205] *= 7.1310000000E+15;
	k[rDD214] *= 5.5420000000E+25;
	k[rDD208] *= 7.1930000000E+11;
	k[rDD43] *= 1.0000000000E+08;
	k[rHP26] *= 1.9030000000E+13;
	k[rLPHP22] *= 1.0030000000E+16;
	k[rHP15] *= 6.6600000000E+08;
	k[rLPHP12] *= 6.7480000000E+18;
	k[rHP11] *= 2.9740000000E+19;
	k[rHP10] *= 1.8870000000E+12;
	k[rLPDD286] *= 4.6000000000E+16;
	k[rDD282] *= 7.3220000000E+14;
	k[rDD288] *= 9.4310000000E+10;
	k[rDD287] *= 1.1680000000E+22;
	k[rDD323] *= 5.8390000000E+18;
	k[rDD324] *= 5.9620000000E+21;
	k[rDD319] *= 5.9340000000E+09;
	k[rDD210] *= 1.5380000000E+24;
	k[rDD421] *= 1.8190000000E+24;
	k[rDD436] *= 8.5070000000E+23;
	k[rDD425] *= 1.2680000000E+24;
	k[rDD422] *= 1.9440000000E+11;
	k[rDD434] *= 2.0710000000E+05;
	k[rDD429] *= 7.1390000000E-02;
	k[rDD424] *= 3.9120000000E+03;
	k[rDD432] *= 5.7120000000E+04;
	k[rDD427] *= 1.3330000000E+11;
	k[rLPDD431] *= 2.0700000000E-02;
	k[rDD405f] *= 2.3720000000E+12;
	k[rDD405b] *= 5.6190000000E+02;
	k[rDD418] *= 7.2920000000E+10;
	k[rDD415] *= 1.4520000000E+17;
	k[rDD416] *= 3.8460000000E+14;
	k[rLPDD406] *= 3.2050000000E+16;
	k[rDD408f] *= 1.7070000000E+11;
	k[rDD408b] *= 4.2560000000E+25;
	k[rDD472] *= 8.7100000000E+38;
	k[rLPDD465f] *= 2.4200000000E+07;
	k[rLPDD465b] *= 7.2484982708E+06;
	k[rLPDD473] *= 1.5350000000E-06;
	k[rLPDD486f] *= 3.1810000000E+14;
	k[rLPDD486b] *= 7.1597693926E+28;
	k[rLPDD502] *= 1.1130000000E+14;
	k[rLPLPDD495] *= 3.8030000000E+19;
	k[rLPDD499] *= 4.3110000000E+23;
	k[rLPDD490] *= 9.4840000000E+20;
	k[rLPR140] *= 1.1290000000E+63;
	k[rLPLPDD107] *= 6.2000000000E-25;
	k[rDD332N] *= 1.0990000000E+08;
	k[rDD294] *= 1.5280000000E+07;
	k[rDD290] *= 8.2480000000E+12;
	k[rLPDD293] *= 1.4280000000E+26;
	k[rHP60] *= 7.7880000000E+06;
	k[rDD239] *= 3.4920000000E+15;
	k[rDD263] *= 1.0440000000E+17;
	k[rDD260f] *= 2.5000000000E+10;
	k[rDD260b] *= 6.7155278570E+09;
	k[rDD259] *= 7.2030000000E+09;
	k[rP000f] *= 1.5000000000E+09;
	k[rP000b] *= 6.1199043605E+09;
	k[rP001f] *= 4.6200000000E+12;
	k[rP001b] *= 2.4680465396E+22;
	k[rP002f] *= 6.8000000000E+21;
	k[rP002b] *= 6.5879507878E+29;
	k[rP003f] *= 1.4500000000E+45;
	k[rP003b] *= 5.9159075488E+45;
	k[rP004f] *= 2.2400000000E+68;
	k[rP004b] *= 5.7144070433E+62;
	k[rP005f] *= 9.6000000000E+67;
	k[rP005b] *= 3.7352568967E+81;
	k[rP006f] *= 1.3800000000E+13;
	k[rP006b] *= 3.0077645189E+23;
	k[rP007f] *= 1.6700000000E+20;
	k[rP007b] *= 6.6010235311E+28;
	k[rP008] *= 2.1600000000E+36;
	k[rP009f] *= 8.2500000000E+43;
	k[rP009b] *= 2.9151378138E+59;
	k[rP010f] *= 1.0700000000E+42;
	k[rP010b] *= 1.5425608195E+58;
	k[rP011f] *= 5.7700000000E+34;
	k[rP011b] *= 5.2012066553E+44;
	k[rK012f] *= 3.2900000000E+03;
	k[rK012b] *= 2.8151364037E+14;
	k[rP013b] *= 3.8029710467E+23;
	k[rP014f] *= 1.8700000000E+04;
	k[rP014b] *= 7.4111081854E+11;
	k[rP015f] *= 9.4500000000E-06;
	k[rP015b] *= 2.7559045191E-05;
	k[rP016f] *= 4.3000000000E+60;
	k[rP016b] *= 6.5130890750E+54;
	k[rK017f] *= 1.2900000000E+05;
	k[rK017b] *= 4.5583789859E+01;
	k[rP018f] *= 7.8000000000E+00;
	k[rP018b] *= 2.8305822277E-02;
	k[rK019f] *= 5.9000000000E+10;
	k[rK019b] *= 1.3748876245E+10;
	k[rK020f] *= 7.1800000000E+10;
	k[rK020b] *= 4.1942662262E+06;
	k[rK021f] *= 1.6500000000E+08;
	k[rK021b] *= 2.2486281887E+06;
	k[rP022b] *= 3.4989212542E+08;
	k[rP023f] *= 4.3000000000E+60;
	k[rP023b] *= 6.1874015827E+54;
	k[rK024f] *= 5.2300000000E+03;
	k[rK024b] *= 1.7556732603E+00;
	k[rP025f] *= 1.3400000000E-01;
	k[rP025b] *= 4.6196306864E-04;
	k[rP026f] *= 3.0100000000E+14;
	k[rP026b] *= 1.8586223399E+09;
	k[rP027f] *= 6.3500000000E+01;
	k[rP027b] *= 9.1474604387E-02;
	k[rP028f] *= 6.5500000000E-05;
	k[rP028b] *= 9.6901068957E-07;
	k[rK100f] *= 1.3400000000E+01;
	k[rK100b] *= 3.3117446183E+14;
	k[rK102f] *= 3.0200000000E+00;
	k[rK102b] *= 3.6667703577E+10;
	k[rP104f] *= 3.8000000000E+04;
	k[rP104b] *= 1.0751660854E+14;
	k[rP105f] *= 3.6000000000E+14;
	k[rP105b] *= 7.9754609459E+25;
	k[rP106f] *= 3.6200000000E+25;
	k[rP106b] *= 9.6548479060E+36;
	k[rP107f] *= 1.2600000000E+01;
	k[rP107b] *= 7.7106058842E+12;
	k[rP108f] *= 8.6000000000E+60;
	k[rP108b] *= 6.6467392724E+54;
	k[rK109f] *= 2.6500000000E+05;
	k[rK109b] *= 4.7781328800E+01;
	k[rK110f] *= 9.6300000000E-01;
	k[rK110b] *= 1.7831960352E-03;
	k[rP111f] *= 8.6000000000E+60;
	k[rP111b] *= 7.5406406999E+54;
	k[rK112f] *= 2.6500000000E+05;
	k[rK112b] *= 5.4207306452E+01;
	k[rK113f] *= 9.6300000000E-01;
	k[rK113b] *= 2.0230130968E-03;
	k[rK114f] *= 1.2800000000E+03;
	k[rK114b] *= 2.2606983071E+14;
	k[rK115f] *= 3.2900000000E+03;
	k[rK115b] *= 2.8326023455E+14;
	k[rP116b] *= 4.8470418645E+18;
	k[rP117b] *= 2.3628374402E+18;
	k[rP118f] *= 7.2000000000E+14;
	k[rP118b] *= 1.0487461179E+20;
	k[rP119f] *= 7.2000000000E+14;
	k[rP119b] *= 5.7999873154E+19;
	k[rP120f] *= 3.6200000000E+25;
	k[rP120b] *= 1.2440515728E+31;
	k[rP121f] *= 3.6200000000E+25;
	k[rP121b] *= 6.0645063853E+30;
	k[rK122f] *= 7.1800000000E+10;
	k[rK122b] *= 2.6267266215E+06;
	k[rK123f] *= 1.6500000000E+08;
	k[rK123b] *= 1.4082395363E+06;
	k[rP124b] *= 2.1912556594E+08;
	k[rK125f] *= 7.1800000000E+10;
	k[rK125b] *= 4.3354949417E+06;
	k[rK126f] *= 1.6500000000E+08;
	k[rK126b] *= 2.3243436664E+06;
	k[rP127b] *= 3.6167364161E+08;
	k[rP128f] *= 2.1500000000E+60;
	k[rP128b] *= 5.9464917059E+54;
	k[rK129f] *= 1.3200000000E+05;
	k[rK129b] *= 8.5172313916E+01;
	k[rP130f] *= 6.7200000000E-02;
	k[rP130b] *= 4.4530158534E-04;
	k[rP131f] *= 2.1500000000E+60;
	k[rP131b] *= 5.1910717013E+54;
	k[rK132f] *= 1.3200000000E+05;
	k[rK132b] *= 7.4352342586E+01;
	k[rP133f] *= 6.7200000000E-02;
	k[rP133b] *= 3.8873214199E-04;
	k[rK200f] *= 2.8800000000E+11;
	k[rK200b] *= 2.1483543223E+08;
	k[rP201f] *= 3.5200000000E+09;
	k[rP201b] *= 7.1773752755E+10;
	k[rP202f] *= 8.6000000000E+60;
	k[rP202b] *= 1.2560982049E+55;
	k[rK203f] *= 2.6500000000E+05;
	k[rK203b] *= 9.0296969496E+01;
	k[rK204f] *= 9.6300000000E-01;
	k[rK204b] *= 3.3698769382E-03;
	k[rK205f] *= 3.2900000000E+03;
	k[rK205b] *= 2.4444249624E+16;
	k[rP206b] *= 2.0390362347E+20;
	k[rP207f] *= 3.6000000000E+14;
	k[rP207b] *= 4.1687273816E+21;
	k[rP208f] *= 3.6200000000E+25;
	k[rP208b] *= 5.2334316593E+32;
	k[rP209f] *= 2.1500000000E+60;
	k[rP209b] *= 5.7147903002E+54;
	k[rK210f] *= 1.3200000000E+05;
	k[rK210b] *= 8.1853626892E+01;
	k[rP211f] *= 6.7200000000E-02;
	k[rP211b] *= 4.2795068192E-04;
	k[rK212f] *= 7.1800000000E+10;
	k[rK212b] *= 6.1509446396E+04;
	k[rK213f] *= 1.6500000000E+08;
	k[rK213b] *= 3.2976417705E+04;
	k[rP214b] *= 5.1312124153E+06;
	k[rP301f] *= 9.5500000000E+08;
	k[rP301b] *= 1.0546908314E+16;
	k[rP302f] *= 1.3900000000E+10;
	k[rP302b] *= 2.4550866312E+23;
	k[rP304f] *= 8.6000000000E+60;
	k[rP304b] *= 8.1347711535E+55;
	k[rP305f] *= 4.0100000000E+05;
	k[rP305b] *= 8.8489841788E+02;
	k[rP306f] *= 2.6900000000E-01;
	k[rP306b] *= 6.0962361138E-03;
	k[rK401f] *= 1.3400000000E+01;
	k[rK401b] *= 3.2105425876E+14;
	k[rK403f] *= 1.3400000000E+01;
	k[rK403b] *= 4.0290414867E+14;
	k[rP405f] *= 3.8000000000E+04;
	k[rP405b] *= 4.8091959456E+13;
	k[rP406f] *= 3.8000000000E+04;
	k[rP406b] *= 8.6959251872E+13;
	k[rP407f] *= 3.2900000000E+03;
	k[rP407b] *= 1.9744900653E+13;
	k[rP408b] *= 3.4325764465E+26;
	k[rP409b] *= 4.3076808769E+26;
	k[rP410f] *= 1.8000000000E+14;
	k[rP410b] *= 2.8481563769E+25;
	k[rP411f] *= 1.8000000000E+14;
	k[rP411b] *= 3.1202065530E+25;
	k[rP412f] *= 3.6200000000E+25;
	k[rP412b] *= 3.7764146576E+36;
	k[rP413f] *= 3.6200000000E+25;
	k[rP413b] *= 4.7391775412E+36;
	k[rP414f] *= 1.2600000000E+01;
	k[rP414b] *= 7.1189266062E+12;
	k[rP415f] *= 1.2600000000E+01;
	k[rP415b] *= 6.2750170621E+12;
	k[rP416f] *= 3.1800000000E+08;
	k[rP416b] *= 6.3789616436E+19;
	k[rP417f] *= 2.3900000000E+08;
	k[rP417b] *= 1.9791193614E+19;
	k[rP418f] *= 1.3900000000E+10;
	k[rP418b] *= 1.8408529969E+27;
	k[rP419f] *= 4.3000000000E+60;
	k[rP419b] *= 8.2369419565E+54;
	k[rK420f] *= 2.0000000000E+05;
	k[rK420b] *= 8.9377808763E+01;
	k[rK421f] *= 1.3400000000E-01;
	k[rK421b] *= 6.1498561740E-04;
	k[rP422f] *= 4.3000000000E+60;
	k[rP422b] *= 8.2369419565E+54;
	k[rK423f] *= 2.6500000000E+05;
	k[rK423b] *= 1.1842559661E+02;
	k[rK424f] *= 9.6300000000E-01;
	k[rK424b] *= 4.4196354445E-03;
	k[rP425f] *= 1.3000000000E+11;
	k[rP425b] *= 5.8424565666E-02;
	k[rP501f] *= 1.3400000000E+01;
	k[rP501b] *= 3.5811034648E+14;
	k[rP502f] *= 3.8000000000E+04;
	k[rP502b] *= 6.9008878969E+11;
	k[rP503f] *= 1.2800000000E+03;
	k[rP503b] *= 7.7618330600E+10;
	k[rP504b] *= 3.8824030263E+17;
	k[rP505f] *= 1.2600000000E+01;
	k[rP505b] *= 4.2972952007E+12;
	k[rP506f] *= 2.1500000000E+60;
	k[rP506b] *= 7.2045750124E+54;
	k[rK507f] *= 2.6500000000E+05;
	k[rK507b] *= 2.0716574153E+02;
	k[rK508f] *= 9.6300000000E-01;
	k[rK508b] *= 7.7314117922E-03;
	k[rK600f] *= 1.2800000000E+03;
	k[rK600b] *= 4.3852288122E+11;
	k[rP601f] *= 1.7200000000E+60;
	k[rP601b] *= 1.2604271375E+54;
	k[rK602f] *= 5.3000000000E+05;
	k[rK602b] *= 9.0608162916E+01;
	k[rK603f] *= 1.9300000000E+00;
	k[rK603b] *= 3.3885134819E-03;
	k[rK700f] *= 1.2800000000E+03;
	k[rK700b] *= 8.9132609449E+10;
	k[rK701f] *= 1.2800000000E+03;
	k[rK701b] *= 1.1012450591E+11;
	k[rP800] *= 6.3700000000E+08;
	k[rP801] *= 9.5500000000E+08;
	k[rP802] *= 1.3900000000E+10;
	k[rCP01f] *= 1.7300000000E+68;
	k[rCP01b] *= 7.6332526441E+61;
	k[rCP02f] *= 2.8000000000E+10;
	k[rCP02b] *= 2.8821974020E+06;
	k[rCP03f] *= 3.3000000000E+11;
	k[rCP03b] *= 1.1622581454E+04;
	k[rCP04] *= 3.3000000000E+11;
	k[rCP05f] *= 4.7700000000E+01;
	k[rCP05b] *= 2.8798925207E-03;
	k[rCP06f] *= 3.0800000000E+03;
	k[rCP06b] *= 3.2559435064E+00;
	k[rCP07f] *= 1.0000000000E+11;
	k[rCP07b] *= 6.2582530698E+07;
	k[rCP08f] *= 1.1000000000E+01;
	k[rCP08b] *= 1.0420300995E+00;
	k[rCP09f] *= 1.8000000000E-04;
	k[rCP09b] *= 3.4048615419E-05;
	k[rCP10b] *= 1.4518195750E+09;
	k[rCP11b] *= 1.5162412910E+12;
	k[rCP12f] *= 7.6600000000E+06;
	k[rCP12b] *= 2.2014248883E+10;
	k[rCP13] *= 3.8900000000E+05;
	k[rCP14] *= 3.7500000000E+33;
	k[rCP15] *= 3.7500000000E+33;
	k[rCP16f] *= 6.8700000000E+52;
	k[rCP16b] *= 5.5097990876E+64;
	k[rCP17] *= 6.3900000000E+26;
	k[rCP18] *= 4.9100000000E+28;
	k[rCP19b] *= 1.2179628270E+15;
	k[rCP20f] *= 4.3400000000E+04;
	k[rCP20b] *= 1.5015495409E+06;
	k[rCP21b] *= 1.5338425964E+13;
	k[rCP23f] *= 2.0000000000E+13;
	k[rCP23b] *= 9.8598712442E+08;
	k[rCP24] *= 1.0000000000E+12;
	k[rCP28] *= 1.4300000000E-16;
	k[rCP30] *= 1.2800000000E+04;
	k[rCP31] *= 3.3700000000E+44;
	k[rCP32f] *= 2.7400000000E+06;
	k[rCP32b] *= 1.0257136572E+12;
	k[rCP33f] *= 1.0000000000E+10;
	k[rCP33b] *= 2.2459405120E+04;
	k[rCP34] *= 7.8600000000E-04;
	k[rI00b] *= 5.2485026205E+25;
	k[rI01f] *= 1.7300000000E+68;
	k[rI01b] *= 1.8449825236E+63;
	k[rI02f] *= 2.8000000000E+10;
	k[rI02b] *= 6.9663668743E+07;
	k[rI03f] *= 3.1600000000E+01;
	k[rI03b] *= 4.5880603280E+08;
	k[rI05f] *= 4.7700000000E+01;
	k[rI05b] *= 6.9607958996E-02;
	k[rI06f] *= 3.0800000000E+03;
	k[rI06b] *= 7.8697236255E+01;
	k[rI07f] *= 1.0000000000E+11;
	k[rI07b] *= 1.5126405584E+09;
	k[rI08f] *= 1.1000000000E+01;
	k[rI08b] *= 2.5186213693E+01;
	k[rI09f] *= 1.8000000000E-04;
	k[rI09b] *= 8.2296634648E-04;
	k[rI12] *= 3.8300000000E+06;
	k[rI15] *= 1.8800000000E+33;
	k[rI17] *= 2.5600000000E+26;
	k[rI18] *= 1.9600000000E+28;
	k[rI19b] *= 2.1780081605E+13;
	k[rI20f] *= 1.7400000000E+04;
	k[rI20b] *= 2.6913158063E+04;
	k[rI23] *= 4.3200000000E+36;
	k[rI31] *= 3.3700000000E+44;
	k[rI32] *= 1.3700000000E+06;
	k[rT01f] *= 2.3100000000E+03;
	k[rT01b] *= 4.4820102636E-03;
	k[rT02f] *= 1.5600000000E+13;
	k[rT02b] *= 1.2986502050E+10;
	k[rT03f] *= 4.3500000000E+22;
	k[rT03b] *= 5.2773986505E+10;
	k[rT04f] *= 5.8300000000E+64;
	k[rT04b] *= 8.4963346413E+55;
	k[rT05f] *= 8.2000000000E+14;
	k[rT05b] *= 4.3614434871E+01;
	k[rT06f] *= 2.1800000000E+04;
	k[rT06b] *= 2.5740303959E+04;
	k[rT07f] *= 6.4700000000E-03;
	k[rT07b] *= 1.2565333739E-03;
	k[rT08f] *= 1.7700000000E+02;
	k[rT08b] *= 3.5302337866E+02;
	k[rT09f] *= 7.8300000000E-01;
	k[rT09b] *= 6.8369388406E-01;
	k[rT10f] *= 3.1400000000E-02;
	k[rT10b] *= 1.7210570482E+04;
	k[rT11f] *= 1.1800000000E-03;
	k[rT11b] *= 1.3441383262E-04;
	k[rT12f] *= 1.6900000000E+09;
	k[rT12b] *= 1.2674962700E+17;
	k[rT13f] *= 1.6600000000E+04;
	k[rT13b] *= 1.9399287458E+07;
	k[rT14f] *= 4.2200000000E+11;
	k[rT14b] *= 1.5060641068E+14;
	k[rT15f] *= 9.3300000000E+01;
	k[rT15b] *= 1.6675285524E+04;
	k[rT16f] *= 7.9400000000E+10;
	k[rT16b] *= 1.0571054566E+14;
	k[rT17b] *= 8.0100757473E+17;
	k[rT18b] *= 4.8247335590E+17;
	k[rT19f] *= 1.1900000000E+06;
	k[rT19b] *= 5.8610185036E+06;
	k[rT20] *= 4.3200000000E+36;
	k[rT21f] *= 3.7600000000E+12;
	k[rT21b] *= 3.1822436227E+12;
	k[rT22f] *= 6.2600000000E+34;
	k[rT22b] *= 2.6747811416E+32;
	k[rT24f] *= 4.2000000000E+03;
	k[rT24b] *= 1.0428348047E+03;
	k[rT25f] *= 1.3000000000E+02;
	k[rT25b] *= 1.8932230751E+01;
	k[rT26f] *= 6.3000000000E+03;
	k[rT26b] *= 1.6064500624E+04;
	k[rT27f] *= 1.0000000000E+04;
	k[rT27b] *= 4.5627585431E+06;
	k[rT28f] *= 8.5500000000E+01;
	k[rT28b] *= 9.8745294140E+02;
	k[rT29f] *= 5.2600000000E+28;
	k[rT29b] *= 6.3725845392E+24;
	k[rT30f] *= 7.2100000000E+33;
	k[rT30b] *= 4.1991617492E+20;
	k[rT31f] *= 2.3700000000E+32;
	k[rT31b] *= 1.5941382341E+20;
	k[rT32b] *= 3.7590981019E+08;
	k[rT33b] *= 1.1057324691E+08;
	k[rT34b] *= 9.6657755389E+08;
	k[rT35f] *= 2.8500000000E-16;
	k[rT35b] *= 4.8973857002E-17;
	k[rT36] *= 2.1000000000E+16;
	k[rT37] *= 4.0900000000E+06;
	k[rT38] *= 5.8400000000E+09;
	k[rT39] *= 2.8900000000E+05;
	k[rT40] *= 1.2000000000E+02;
	k[rT41] *= 4.0900000000E+01;
	k[rT42] *= 3.4900000000E-11;
	k[rT43f] *= 1.0100000000E+71;
	k[rT43b] *= 1.5737620878E+66;
	k[rT44f] *= 1.1500000000E+11;
	k[rT44b] *= 4.1803986880E+08;
	k[rT45f] *= 2.8100000000E+10;
	k[rT45b] *= 5.9912631081E+07;
	k[rT46f] *= 2.9500000000E+03;
	k[rT46b] *= 1.1012916128E+02;
	k[rT47] *= 2.9000000000E+10;
	k[rT50f] *= 8.3900000000E+60;
	k[rT50b] *= 5.9493864399E+54;
	k[rT51f] *= 3.9100000000E+05;
	k[rT51b] *= 6.4682847912E+01;
	k[rT52f] *= 2.6000000000E+10;
	k[rT52b] *= 2.5227716017E+06;
	k[rT53f] *= 1.3600000000E+01;
	k[rT53b] *= 2.3105306687E-02;
	k[rT54f] *= 1.7900000000E-05;
	k[rT54b] *= 5.4416006531E-06;
	k[rT55b] *= 1.6480389513E+20;
	k[rT56b] *= 3.6132196175E+17;
	k[rT57b] *= 6.9312742321E+13;
	k[rT58f] *= 8.5700000000E+17;
	k[rT58b] *= 2.0524277087E+22;
	k[rT59] *= 2.5500000000E+10;
	k[rT60] *= 2.5500000000E+10;
	k[rE01b] *= 9.7891110020E+17;
	k[rE02b] *= 3.0096212144E+21;
	k[rE03b] *= 1.1999372607E+25;
	k[rE04f] *= 2.3100000000E+03;
	k[rE04b] *= 6.1576309766E-06;
	k[rE05f] *= 7.8300000000E-01;
	k[rE05b] *= 9.3929607283E-04;
	k[rE06f] *= 4.8300000000E-05;
	k[rE06b] *= 4.1554038409E-10;
	k[rE07f] *= 1.9600000000E-03;
	k[rE07b] *= 9.8904120951E-09;
	k[rE08f] *= 4.4700000000E+03;
	k[rE08b] *= 3.9494269248E-01;
	k[rE09f] *= 8.0300000000E+00;
	k[rE09b] *= 6.3577513555E-02;
	k[rE10f] *= 7.5300000000E-04;
	k[rE10b] *= 1.1904806541E-05;
	k[rE11f] *= 1.7200000000E+11;
	k[rE11b] *= 4.2169346325E+00;
	k[rE12f] *= 3.7900000000E+06;
	k[rE12b] *= 1.0739437552E+04;
	k[rE13b] *= 1.1039793641E+09;
	k[rE14b] *= 1.6361457188E+11;
	k[rE15f] *= 3.3100000000E+09;
	k[rE15b] *= 4.0209968985E+12;
	k[rE16f] *= 3.7000000000E+13;
	k[rE16b] *= 1.4870776226E+14;
	k[rE17f] *= 3.1700000000E+10;
	k[rE17b] *= 8.0200953199E+08;
	k[rE18f] *= 3.1700000000E+10;
	k[rE18b] *= 2.4313367079E+09;
	k[rE19] *= 7.0000000000E+09;
	k[rE30] *= 4.8200000000E+33;
	k[rE31] *= 3.8800000000E+01;
	k[rE32] *= 1.2100000000E+45;
	k[rE33] *= 2.4500000000E+25;
	k[rE34] *= 1.2300000000E+13;
	k[rE35] *= 5.4500000000E-02;
	k[rE37] *= 6.3000000000E+14;
	k[rST01] *= 2.4000000000E+14;
	k[rST02f] *= 5.1000000000E+04;
	k[rST02b] *= 3.2882754328E+00;
	k[rST03f] *= 1.1400000000E+02;
	k[rST03b] *= 3.0178171755E+02;
	k[rST04f] *= 1.8800000000E+33;
	k[rST04b] *= 3.7188357721E+32;
	k[rST05f] *= 1.8800000000E+33;
	k[rST05b] *= 1.2267086362E+32;
	k[rST06f] *= 7.8300000000E-01;
	k[rST06b] *= 8.8911775253E-03;
	k[rST10f] *= 1.3300000000E+10;
	k[rST10b] *= 2.6187006539E+06;
	k[rST11f] *= 1.0300000000E+10;
	k[rST11b] *= 2.0850645484E+07;
	k[rST12f] *= 6.7000000000E+02;
	k[rST12b] *= 5.5513158673E+01;
	k[rST13] *= 3.0300000000E+08;
	k[rST14f] *= 4.5800000000E+13;
	k[rST14b] *= 2.9579446639E+12;
	k[rST00b] *= 7.2340731681E+07;
	k[rXY00f] *= 2.5000000000E+18;
	k[rXY00b] *= 1.0584848633E+15;
	k[rXY01f] *= 4.3200000000E+29;
	k[rXY01b] *= 2.4885817320E+17;
	k[rXY02f] *= 1.2900000000E-02;
	k[rXY02b] *= 1.2741968022E-03;
	k[rXY03f] *= 2.3600000000E-03;
	k[rXY03b] *= 1.3672597006E-04;
	k[rXY04f] *= 3.5400000000E+02;
	k[rXY04b] *= 3.5909595731E+02;
	k[rXY05f] *= 4.3600000000E+04;
	k[rXY05b] *= 2.6183079225E+04;
	k[rXY06f] *= 1.8700000000E+02;
	k[rXY06b] *= 1.6998488070E+04;
	k[rXY07f] *= 8.4400000000E+11;
	k[rXY07b] *= 1.5319708690E+14;
	k[rXY09f] *= 4.6200000000E+03;
	k[rXY09b] *= 3.7531845431E-03;
	k[rXY10f] *= 1.5700000000E+00;
	k[rXY10b] *= 6.9907410294E-01;
	k[rXY11] *= 1.8200000000E+05;
	k[rXY12] *= 8.2000000000E+14;
	k[rXY13f] *= 5.8300000000E+64;
	k[rXY13b] *= 7.9321711003E+55;
	k[rXY14f] *= 4.3700000000E+15;
	k[rXY14b] *= 2.3757277888E+18;
	k[rXY15f] *= 5.9900000000E+20;
	k[rXY15b] *= 1.1442380922E+14;
	k[rXY16f] *= 1.9700000000E+19;
	k[rXY16b] *= 3.8323543459E+13;
	k[rXY18f] *= 1.3800000000E-01;
	k[rXY18b] *= 1.4917919852E-01;
	k[rXY19f] *= 6.5700000000E+00;
	k[rXY19b] *= 3.0056748401E-02;
	k[rXY201] *= 2.2800000000E+10;
	k[rXY202] *= 3.1300000000E+15;
	k[rXY203] *= 1.0300000000E+14;
	k[rXY22] *= 4.3200000000E+36;
	k[rXY23f] *= 1.2500000000E+18;
	k[rXY23b] *= 1.0731519576E+15;
	k[rXY24] *= 2.1600000000E+29;
	k[rXY25] *= 2.1000000000E+16;
	k[rXY26f] *= 6.4700000000E-03;
	k[rXY26b] *= 1.2958586951E-03;
	k[rXY27f] *= 1.1800000000E-03;
	k[rXY27b] *= 1.3862053915E-04;
	k[rXY28f] *= 1.7700000000E+02;
	k[rXY28b] *= 3.6407183790E+02;
	k[rXY29f] *= 2.1800000000E+04;
	k[rXY29b] *= 2.6545889980E+04;
	k[rXY30f] *= 9.3300000000E+01;
	k[rXY30b] *= 1.7197166576E+04;
	k[rXY31f] *= 4.2200000000E+11;
	k[rXY31b] *= 1.5531989112E+14;
	k[rXY33] *= 4.0900000000E+06;
	k[rXY34] *= 5.8400000000E+09;
	k[rXY35] *= 2.8900000000E+05;
	k[rXY36] *= 1.2000000000E+02;
	k[rXY37] *= 4.0900000000E+01;
	k[rXY38] *= 3.4900000000E-11;
	k[rXY39f] *= 2.3100000000E+03;
	k[rXY39b] *= 8.2660146550E-06;
	k[rXY40f] *= 2.3100000000E+03;
	k[rXY40b] *= 2.8887454010E-03;
	k[rXY41f] *= 7.8300000000E-01;
	k[rXY41b] *= 1.5357183416E-03;
	k[rXY42] *= 7.8300000000E-01;
	k[rXY43f] *= 4.3700000000E+15;
	k[rXY43b] *= 4.9836580785E+18;
	k[rXY44] *= 5.9900000000E+20;
	k[rXY45f] *= 1.9700000000E+19;
	k[rXY45b] *= 2.9093665114E+13;
	k[rXY47f] *= 1.3800000000E-01;
	k[rXY47b] *= 3.1293910075E-01;
	k[rXY48] *= 1.1900000000E+06;
	k[rXY50] *= 4.2000000000E+16;
	k[rXY51] *= 8.1800000000E+06;
	k[rXY52] *= 1.1700000000E+10;
	k[rXY53] *= 5.7800000000E+05;
	k[rXY54] *= 2.4000000000E+02;
	k[rXY55] *= 8.1800000000E+01;
	k[rXY56] *= 6.9800000000E-11;
	k[rXY57f] *= 4.6200000000E+03;
	k[rXY57b] *= 5.9828399303E-06;
	k[rXY58] *= 1.5700000000E+00;
	k[rN01f] *= 2.3100000000E+03;
	k[rN01b] *= 2.5952295683E-05;
	k[rN02f] *= 7.8300000000E-01;
	k[rN02b] *= 1.0967302857E-03;
	k[rN04f] *= 1.2500000000E+18;
	k[rN04b] *= 7.7463001274E+12;
	k[rN05f] *= 3.2000000000E+34;
	k[rN05b] *= 2.7785868184E+20;
	k[rN06f] *= 5.8300000000E+64;
	k[rN06b] *= 8.1687995828E+55;
	k[rN07f] *= 8.2000000000E+14;
	k[rN07b] *= 6.5241355800E+02;
	k[rN08f] *= 6.4700000000E-03;
	k[rN08b] *= 9.3538573948E-06;
	k[rN09f] *= 1.1800000000E-03;
	k[rN09b] *= 1.0006004205E-06;
	k[rN10f] *= 1.7700000000E+02;
	k[rN10b] *= 2.6279686713E+00;
	k[rN11f] *= 2.1800000000E+04;
	k[rN11b] *= 1.9161539004E+02;
	k[rN12f] *= 4.2200000000E+11;
	k[rN12b] *= 1.1211408447E+12;
	k[rN13f] *= 9.3300000000E+01;
	k[rN13b] *= 1.2413378431E+02;
	k[rN14] *= 1.1000000000E+10;
	k[rN15] *= 1.4700000000E+10;
	k[rN16b] *= 1.0052196131E+18;
	k[rN18f] *= 1.1900000000E+06;
	k[rN18b] *= 7.3552497358E+06;
	k[rN19] *= 4.3200000000E+36;
	k[rN20f] *= 3.7600000000E+12;
	k[rN20b] *= 1.7159707068E+14;
	k[rN21f] *= 3.0800000000E+06;
	k[rN21b] *= 1.9552148774E+03;
	k[rN22f] *= 5.2600000000E+28;
	k[rN22b] *= 2.7382161014E+26;
	k[rN23f] *= 7.2100000000E+33;
	k[rN23b] *= 3.2171025483E+20;
	k[rN24b] *= 1.6152352136E+10;
	k[rN25b] *= 4.7511875788E+09;
	k[rN26b] *= 4.1532571362E+10;
	k[rN27f] *= 2.8500000000E-16;
	k[rN27b] *= 2.1043424841E-15;
	k[rN28] *= 2.1000000000E+16;
	k[rN29] *= 4.0900000000E+06;
	k[rN30] *= 5.8400000000E+09;
	k[rN31] *= 2.8900000000E+05;
	k[rN32] *= 1.2000000000E+02;
	k[rN33] *= 4.0900000000E+01;
	k[rN34] *= 3.4900000000E-11;
	k[rOX00f] *= 1.2900000000E+61;
	k[rOX00b] *= 8.0660235973E+54;
	k[rOX02f] *= 5.0000000000E+75;
	k[rOX02b] *= 2.7802587455E+86;
	k[rOX03f] *= 6.0200000000E+05;
	k[rOX03b] *= 8.7814943551E+01;
	k[rOX04f] *= 2.3400000000E+01;
	k[rOX04b] *= 3.5054860727E-02;
	k[rOX05f] *= 1.3200000000E-01;
	k[rOX05b] *= 5.9403631765E+04;
	k[rOX06f] *= 4.3400000000E+11;
	k[rOX06b] *= 3.8490104908E+08;
	k[rOX99f] *= 2.7400000000E+05;
	k[rOX99b] *= 2.4619859479E+08;
	k[rOX07f] *= 1.9900000000E+04;
	k[rOX07b] *= 1.7880846848E+07;
	k[rOX08f] *= 2.0300000000E+09;
	k[rOX08b] *= 1.2500560065E+17;
	k[rOX09f] *= 5.8000000000E+13;
	k[rOX09b] *= 7.5724507824E+08;
	k[rOX10f] *= 4.0000000000E+10;
	k[rOX10b] *= 3.4223477473E+06;
	k[rOX11f] *= 8.5700000000E+17;
	k[rOX11b] *= 1.7896397866E+22;
	k[rOX12f] *= 3.0000000000E+10;
	k[rOX12b] *= 5.2188245827E+17;
	k[rOX13b] *= 1.4370279961E+20;
	k[rOX14b] *= 3.1505916423E+17;
	k[rOX15b] *= 6.0438105009E+13;
	k[rOX16f] *= 3.8900000000E-06;
	k[rOX16b] *= 1.4511636996E-05;
	k[rOX17f] *= 6.5900000000E+15;
	k[rOX17b] *= 1.3972060163E+03;
	k[rOX18f] *= 1.0100000000E+71;
	k[rOX18b] *= 1.4737470456E+66;
	k[rOX19f] *= 6.8300000000E-02;
	k[rOX19b] *= 2.3250079098E-04;
	k[rOX20f] *= 1.7300000000E-02;
	k[rOX20b] *= 6.0479791809E-04;
	k[rOX21f] *= 3.7000000000E-07;
	k[rOX21b] *= 2.3145507889E-06;
	k[rOX22f] *= 7.3200000000E+10;
	k[rOX22b] *= 1.5149637533E+09;
	k[rOX23f] *= 2.9000000000E+10;
	k[rOX23b] *= 1.8592380137E-04;
	k[rOX24f] *= 6.3400000000E+10;
	k[rOX24b] *= 5.2814831459E+13;
	k[rOX25f] *= 6.5100000000E+04;
	k[rOX25b] *= 1.0783545407E+05;
	k[rOX26f] *= 7.4000000000E+11;
	k[rOX26b] *= 9.9091879383E-02;
	k[rOX27f] *= 4.3000000000E+06;
	k[rOX27b] *= 2.1555116600E-01;
	k[rOX28] *= 3.0000000000E+10;
	k[rOX30f] *= 1.8300000000E+10;
	k[rOX30b] *= 3.1407237352E+12;
	k[rOX31f] *= 1.8300000000E+10;
	k[rOX31b] *= 3.1219053405E+17;
	k[rOX32f] *= 2.2000000000E-01;
	k[rOX32b] *= 2.7428177703E+04;
	k[rOX33f] *= 8.5700000000E+17;
	k[rOX33b] *= 2.7654693660E+21;
	k[rOX34f] *= 8.5700000000E+17;
	k[rOX34b] *= 2.4376382025E+21;
	k[rOX35] *= 3.0000000000E+10;
	k[rOX36] *= 3.0000000000E+10;
	k[rOX37b] *= 2.2205903841E+19;
	k[rOX38b] *= 1.9573516232E+19;
	k[rOX39b] *= 4.8685018832E+16;
	k[rOX40b] *= 4.2913677966E+16;
	k[rOX41f] *= 8.6200000000E+15;
	k[rOX41b] *= 5.2489324223E+03;
	k[rOX42f] *= 1.0100000000E+71;
	k[rOX42b] *= 1.0160881342E+66;
	k[rOX43f] *= 6.8300000000E-02;
	k[rOX43b] *= 1.6029975810E-04;
	k[rOX44f] *= 1.7300000000E-02;
	k[rOX44b] *= 4.1698335546E-04;
	k[rOX45f] *= 3.7000000000E-07;
	k[rOX45b] *= 1.5957878251E-06;
	k[rOX46f] *= 2.9000000000E+10;
	k[rOX46b] *= 1.8719655890E-02;
	k[rOX48] *= 6.5100000000E+04;
	k[rOX50] *= 8.5700000000E+17;
	k[rOX51] *= 8.5700000000E+17;
	k[rOX52] *= 8.5700000000E+17;
	k[rOX53] *= 8.5700000000E+17;
	k[rOX54] *= 8.5700000000E+17;
	k[rOX60] *= 1.1000000000E-01;
	k[rOX61] *= 1.1000000000E-01;
	k[rOX62] *= 2.2000000000E-01;
	k[rOX63] *= 1.7600000000E-01;
	k[rOX64] *= 2.2000000000E-01;
	k[rOX65] *= 2.2000000000E-01;
	k[rOX66] *= 2.2000000000E-01;
	k[rNO178f] *= 2.7000000000E+10;
	k[rNO178b] *= 4.3388802917E+10;
	k[rNO179f] *= 9.0000000000E+06;
	k[rNO179b] *= 8.9706941899E+05;
	k[rNO180f] *= 3.3600000000E+10;
	k[rNO180b] *= 1.6842590024E+12;
	k[rNO181f] *= 1.4000000000E+09;
	k[rNO181b] *= 1.9769877959E+05;
	k[rNO182f] *= 2.9000000000E+10;
	k[rNO182b] *= 2.5400577397E+05;
	k[rNO183f] *= 3.8700000000E+11;
	k[rNO183b] *= 1.0866779635E+05;
	k[rNO184f] *= 2.0000000000E+09;
	k[rNO184b] *= 2.9275287484E+06;
	k[rNO186f] *= 2.1100000000E+09;
	k[rNO186b] *= 7.3558269264E+10;
	k[rNO187f] *= 1.0600000000E+14;
	k[rNO187b] *= 2.6359098753E+21;
	k[rNO188f] *= 3.9000000000E+09;
	k[rNO188b] *= 1.0792457321E+07;
	k[rNO189f] *= 1.3200000000E+11;
	k[rNO189b] *= 7.2634609588E+05;
	k[rNO190b] *= 2.1728852178E+12;
	k[rNO191f] *= 3.2000000000E+10;
	k[rNO191b] *= 5.9124157137E+10;
	k[rNO192b] *= 1.9563032926E+14;
	k[rNO193f] *= 2.0000000000E+06;
	k[rNO193b] *= 3.7949445387E+07;
	k[rNO194f] *= 4.6100000000E+02;
	k[rNO194b] *= 8.9664608808E+03;
	k[rNO195f] *= 1.2800000000E+03;
	k[rNO195b] *= 1.3826132234E+02;
	k[rNO196b] *= 1.3094290066E+12;
	k[rNO197f] *= 2.0000000000E+10;
	k[rNO197b] *= 1.9049155742E+13;
	k[rNO198f] *= 2.1600000000E+10;
	k[rNO198b] *= 3.7616193806E+10;
	k[rNO199f] *= 3.6500000000E+11;
	k[rNO199b] *= 2.2637270869E+18;
	k[rNO200b] *= 1.2959748343E+08;
	k[rNO201b] *= 1.6479578930E+13;
	k[rNO202f] *= 4.0000000000E+10;
	k[rNO202b] *= 2.9460702266E+09;
	k[rNO203f] *= 9.0000000000E+04;
	k[rNO203b] *= 6.8074751794E+04;
	k[rNO204b] *= 9.7425878465E+04;
	k[rNO205f] *= 1.3000000000E+11;
	k[rNO205b] *= 3.8379891517E+07;
	k[rNO206b] *= 2.0937266561E+09;
	k[rNO207b] *= 1.0099381235E+09;
	k[rNO208b] *= 1.6237968686E+09;
	k[rNO209b] *= 3.4437574222E+09;
	k[rNO210b] *= 1.4146630067E+10;
	k[rNO211b] *= 3.1641992010E+12;
	k[rNO212f] *= 4.4800000000E+13;
	k[rNO212b] *= 1.1038289431E+18;
	k[rNO213b] *= 1.3883872366E+08;
	k[rNO214f] *= 9.0000000000E+08;
	k[rNO214b] *= 8.5215952651E+06;
	k[rNO215f] *= 1.3000000000E+04;
	k[rNO215b] *= 1.2641022364E+03;
	k[rNO216f] *= 1.0000000000E+10;
	k[rNO216b] *= 5.7565969606E+08;
	k[rNO217b] *= 2.7111042230E+11;
	k[rNO218b] *= 1.0453566818E+17;
	k[rNO219f] *= 8.0000000000E+09;
	k[rNO219b] *= 2.5787783168E+12;
	k[rNO220f] *= 6.1400000000E+09;
	k[rNO220b] *= 3.1907062411E+13;
	k[rNO221f] *= 2.9500000000E+02;
	k[rNO221b] *= 9.7657700109E+05;
	k[rNO222b] *= 1.5870420808E+06;
	k[rNO223b] *= 6.7133230431E+04;
	k[rNO224b] *= 1.2338609554E+03;
	k[rNO225b] *= 2.1705234378E+06;
	k[rNO226f] *= 2.0000000000E+09;
	k[rNO226b] *= 1.5016004241E+09;
	k[rNO227f] *= 3.1000000000E+11;
	k[rNO227b] *= 3.0522373799E+03;
	k[rNO228f] *= 1.9000000000E+14;
	k[rNO228b] *= 1.4649695236E+15;
	k[rNO229f] *= 3.8000000000E+15;
	k[rNO229b] *= 4.5997936024E+16;
	k[rNO230f] *= 1.0400000000E+26;
	k[rNO230b] *= 1.3466251020E+20;
	k[rNO231f] *= 2.0300000000E+01;
	k[rNO231b] *= 9.3995763477E+03;
	k[rNO232f] *= 5.0700000000E+00;
	k[rNO232b] *= 2.9185288262E-03;
	k[rNO233f] *= 3.9100000000E+06;
	k[rNO233b] *= 6.9276269181E+02;
	k[rNO234f] *= 1.1000000000E+03;
	k[rNO234b] *= 4.6317021066E+06;
	k[rNO235] *= 4.4000000000E+00;
	k[rNO236f] *= 1.6000000000E-01;
	k[rNO236b] *= 2.1320662570E-03;
	k[rNO238f] *= 6.0000000000E+10;
	k[rNO238b] *= 3.2249916840E+07;
	k[rNO239f] *= 6.3000000000E+10;
	k[rNO239b] *= 5.2402458901E+10;
	k[rNO240f] *= 3.1200000000E+06;
	k[rNO240b] *= 4.4323704327E+09;
	k[rNO242f] *= 1.0000000000E+10;
	k[rNO242b] *= 4.7779848137E+11;
	k[rNO243f] *= 1.0000000000E+08;
	k[rNO243b] *= 3.4835996596E+09;
	k[rNO244b] *= 2.5396777936E+10;
	k[rNO245b] *= 1.3648298329E+11;
	k[rNO246b] *= 9.3600725705E+13;
	k[rNO247b] *= 1.7124685946E+16;
	k[rNO248b] *= 1.3206580898E+12;
	k[rNO249f] *= 3.1000000000E+14;
	k[rNO249b] *= 1.4756130541E+21;
	k[rNO250f] *= 2.9000000000E+11;
	k[rNO250b] *= 2.4130363332E+13;
	k[rNO251f] *= 3.8000000000E+10;
	k[rNO251b] *= 6.9568023699E+17;
	k[rNO252f] *= 3.1000000000E+14;
	k[rNO252b] *= 1.0758605005E+21;
	k[rNO253f] *= 2.9000000000E+11;
	k[rNO253b] *= 1.7593301102E+13;
	k[rNO254f] *= 3.8000000000E+10;
	k[rNO254b] *= 5.0721622845E+17;
	k[rNO255f] *= 9.6000000000E+10;
	k[rNO255b] *= 1.9782922307E+11;
	k[rNO256f] *= 1.0000000000E+09;
	k[rNO256b] *= 4.2288439657E+09;
	k[rNO257b] *= 2.4697660271E+02;
	k[rNO258b] *= 3.0646463697E+06;
	k[rNO259b] *= 1.0097759616E+01;
	k[rNO260b] *= 5.0782150163E+03;
	k[rNO261b] *= 5.9037616232E+04;
	k[rNO262f] *= 9.8000000000E+04;
	k[rNO262b] *= 5.5134956480E+03;
	k[rNO263f] *= 1.5000000000E+05;
	k[rNO263b] *= 1.4764153289E+01;
	k[rNO264f] *= 2.2000000000E+03;
	k[rNO264b] *= 1.7806968531E+01;
	k[rNO265f] *= 2.2500000000E+04;
	k[rNO265b] *= 5.2410499750E-03;
	k[rNO266f] *= 1.0500000000E+02;
	k[rNO266b] *= 1.4489867013E+00;
	k[rNO267f] *= 3.3000000000E+04;
	k[rNO267b] *= 4.6768074920E+03;
	k[rNO268f] *= 3.3000000000E+03;
	k[rNO268b] *= 4.2977345395E+03;
	k[rNO269f] *= 1.1800000000E+13;
	k[rNO269b] *= 8.6775822032E+02;
	k[rNO270f] *= 2.1000000000E+12;
	k[rNO270b] *= 5.4601430394E+11;
	k[rNO271f] *= 2.7000000000E+08;
	k[rNO271b] *= 1.2271666890E+03;
	k[rNO272f] *= 1.7000000000E+11;
	k[rNO272b] *= 1.0296021231E+04;
	k[rNO273] *= 2.0000000000E+04;
	k[rNO274b] *= 1.3897160079E+11;
	k[rNO275f] *= 6.1000000000E+11;
	k[rNO275b] *= 1.2930672014E+14;
	k[rNO276f] *= 3.7000000000E+09;
	k[rNO276b] *= 3.7216033418E+10;
	k[rNO277f] *= 5.4000000000E+02;
	k[rNO277b] *= 2.9019836947E+00;
	k[rNO278f] *= 5.0000000000E+04;
	k[rNO278b] *= 2.7595080532E+03;
	k[rNO279f] *= 9.4000000000E+03;
	k[rNO279b] *= 2.9629291835E+01;
	k[rNO280f] *= 1.0000000000E+10;
	k[rNO280b] *= 1.7495095849E+07;
	k[rNO281f] *= 6.1600000000E+12;
	k[rNO281b] *= 8.8583920205E+13;
	k[rNO282f] *= 3.2500000000E+09;
	k[rNO282b] *= 7.7093522462E+11;
	k[rNO283f] *= 3.0000000000E+09;
	k[rNO283b] *= 2.6896823347E+04;
	k[rG09] = 6.7100000000E+10;
	k[rG15] = 8.0000000000E+10;
	k[rG16] = 2.0000000000E+10;
	k[rG36] = 7.7590000000E+09;
	k[rG35] = 4.5590000000E+10;
	k[rG38f] = 3.0000000000E+10;
	k[rG34f] = 7.0000000000E+10;
	k[rG52] = 6.3410000000E+10;
	k[rG53] = 2.1410000000E+10;
	k[rG72] = 1.0000000000E+11;
	k[rLPG65] = 1.4500000000E+09;
	k[r32] = 1.4700000000E+11;
	k[r35] = 3.0200000000E+10;
	k[rDD54] = 2.0600000000E+11;
	k[rG130] = 1.6140000000E+10;
	k[rR037] = 4.9170000000E+11;
	k[rG109] = 1.0000000000E+11;
	k[rR022] = 2.7500000000E+10;
	k[rG110] = 1.0000000000E+11;
	k[rR024] = 1.0000000000E+10;
	k[rP013f] = 6.0000000000E+09;
	k[rP022f] = 2.5000000000E+09;
	k[rP116f] = 6.0000000000E+09;
	k[rP117f] = 6.0000000000E+09;
	k[rP124f] = 2.5000000000E+09;
	k[rP127f] = 2.5000000000E+09;
	k[rP206f] = 6.0000000000E+09;
	k[rP214f] = 2.5000000000E+09;
	k[rP408f] = 6.0000000000E+09;
	k[rP409f] = 6.0000000000E+09;
	k[rP504f] = 6.0000000000E+09;
	k[rCP10f] = 6.0000000000E+09;
	k[rCP11f] = 6.0000000000E+09;
	k[rCP19f] = 7.0000000000E+10;
	k[rCP21f] = 3.1000000000E+10;
	k[rCP22] = 1.0200000000E+10;
	k[rCP25] = 6.6700000000E+09;
	k[rCP26] = 3.3300000000E+09;
	k[rCP27] = 1.6700000000E+09;
	k[rCP29] = 1.3300000000E+09;
	k[rI00f] = 1.0000000000E+10;
	k[rI19f] = 2.8000000000E+10;
	k[rI21] = 1.2400000000E+10;
	k[rI22] = 4.0800000000E+09;
	k[rT17f] = 2.2800000000E+11;
	k[rT18f] = 2.0000000000E+10;
	k[rT32f] = 1.3300000000E+10;
	k[rT33f] = 6.6700000000E+09;
	k[rT34f] = 3.3300000000E+09;
	k[rT55f] = 1.0000000000E+11;
	k[rT56f] = 3.0000000000E+10;
	k[rT57f] = 3.0000000000E+10;
	k[rE01f] = 3.6100000000E+10;
	k[rE02f] = 2.0000000000E+10;
	k[rE03f] = 2.0000000000E+10;
	k[rE13f] = 1.6700000000E+09;
	k[rE14f] = 2.4100000000E+10;
	k[rE36] = 1.0000000000E+08;
	k[rST00f] = 1.0000000000E+10;
	k[rXY17] = 2.0000000000E+10;
	k[rXY46] = 2.0000000000E+10;
	k[rN16f] = 2.2800000000E+11;
	k[rN17] = 2.0000000000E+10;
	k[rN24f] = 1.3300000000E+10;
	k[rN25f] = 6.6700000000E+09;
	k[rN26f] = 3.3300000000E+09;
	k[rOX13f] = 1.0000000000E+11;
	k[rOX14f] = 3.0000000000E+10;
	k[rOX15f] = 3.0000000000E+10;
	k[rOX37f] = 1.0000000000E+11;
	k[rOX38f] = 1.0000000000E+11;
	k[rOX39f] = 3.0000000000E+10;
	k[rOX40f] = 3.0000000000E+10;
	k[rOX47] = 1.6800000000E+11;
	k[rNO190f] = 4.0000000000E+10;
	k[rNO192f] = 2.0000000000E+10;
	k[rNO196f] = 1.5000000000E+10;
	k[rNO200f] = 3.0000000000E+09;
	k[rNO201f] = 3.9000000000E+10;
	k[rNO204f] = 3.3000000000E+08;
	k[rNO206f] = 5.0000000000E+09;
	k[rNO207f] = 2.5000000000E+10;
	k[rNO208f] = 7.0000000000E+10;
	k[rNO209f] = 5.0000000000E+10;
	k[rNO210f] = 2.0000000000E+10;
	k[rNO211f] = 2.5000000000E+10;
	k[rNO213f] = 2.5000000000E+10;
	k[rNO217f] = 7.7000000000E+10;
	k[rNO218f] = 4.0000000000E+10;
	k[rNO222f] = 2.3500000000E+10;
	k[rNO223f] = 5.4000000000E+10;
	k[rNO224f] = 2.5000000000E+09;
	k[rNO225f] = 2.0000000000E+10;
	k[rNO244f] = 1.9000000000E+10;
	k[rNO245f] = 2.9000000000E+10;
	k[rNO246f] = 4.1000000000E+10;
	k[rNO247f] = 1.6200000000E+10;
	k[rNO248f] = 2.4600000000E+10;
	k[rNO257f] = 2.2000000000E+10;
	k[rNO258f] = 2.0000000000E+09;
	k[rNO259f] = 1.2000000000E+10;
	k[rNO260f] = 1.2000000000E+10;
	k[rNO261f] = 1.0000000000E+11;
	k[rNO274f] = 9.0000000000E+09;
	kTroe0 = 6.4390000000E+14	 * k[r15f];
	kTroeInf = 4.7040000000E+09	 * w[r15f];
	fcTroe = 0.5 * exp( -temp * 1e+30 ) + 0.5 * exp( -temp * 1e-30 );

	k[r15f] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM15] );
	kTroe0 = 9.1328854916E+17	 * k[r15b];
	kTroeInf = 6.6720132556E+12	 * w[r15b];

	k[r15b] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM15] );
	kTroe0 = 2.1760000000E+22	 * k[r22];
	kTroeInf = 1.7470000000E+13	 * w[r22];
	fcTroe = 0.57 * exp( -temp * 1e+30 ) + 0.43 * exp( -temp * 1e-30 );

	k[r22] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM16] );
	kTroe0 = 1.0400000000E+20	 * k[rG14];
	kTroeInf = 6.0000000000E+11	 * w[rG14];
	fcTroe = 0.438 * exp( -temp * 0.010989 ) + 0.562 * exp( -temp * 0.00017135 ) + 1 * exp( -8552 * RGAS* rt_inv );

	k[rG14] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM3] );
	kTroe0 = 3.8200000000E+25	 * k[rLPG60f];
	kTroeInf = 1.0100000000E+05	 * w[rLPG60f];
	fcTroe = 0.955 * exp( -temp * 0.00113623 ) + 0.045 * exp( -temp * 4e-10 ) + 1 * exp( -1.786e+09 * RGAS* rt_inv );

	k[rLPG60f] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM0] );
	kTroe0 = 1.6260000000E+36	 * k[rLPG60b];
	kTroeInf = 4.2980000000E+15	 * w[rLPG60b];

	k[rLPG60b] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM0] );
	kTroe0 = 1.5500000000E+18	 * k[r27];
	kTroeInf = 1.8000000000E+07	 * w[r27];
	fcTroe = 1 * exp( -0 * RGAS* rt_inv );

	k[r27] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM20] );
	kTroe0 = 2.6900000000E+27	 * k[rLPG24];
	kTroeInf = 8.1000000000E+08	 * w[rLPG24];
	fcTroe = 0.4093 * exp( -temp * 0.00363636 ) + 0.5907 * exp( -temp * 0.000815661 ) + 1 * exp( -5185 * RGAS* rt_inv );

	k[rLPG24] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM3] );
	kTroe0 = 6.3030000000E+58	 * k[rLPG44];
	kTroeInf = 1.5470000000E+43	 * w[rLPG44];
	fcTroe = 0.242 * exp( -temp * 0.0106383 ) + 0.758 * exp( -temp * 0.000643087 ) + 1 * exp( -4200 * RGAS* rt_inv );

	k[rLPG44] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM3] );
	kTroe0 = 4.3040000000E+25	 * k[rG115f];
	kTroeInf = 1.1610000000E+07	 * w[rG115f];
	fcTroe = 0.2122 * exp( -temp * -9.79432e-05 );

	k[rG115f] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM3] );
	kTroe0 = 9.4097780916E+27	 * k[rG115b];
	kTroeInf = 2.5382788951E+09	 * w[rG115b];

	k[rG115b] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM3] );
	kTroe0 = 2.3320000000E+33	 * k[rG155f];
	kTroeInf = 1.5740000000E+06	 * w[rG155f];
	fcTroe = 1.569 * exp( -temp * 0.00334448 ) + -0.569 * exp( -temp * -0.000109325 ) + 1 * exp( -152.4 * RGAS* rt_inv );

	k[rG155f] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM3] );
	kTroe0 = 4.2990018456E+36	 * k[rG155b];
	kTroeInf = 2.9016418975E+09	 * w[rG155b];

	k[rG155b] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM3] );
	kTroe0 = 3.1690000000E-68	 * k[rLPDD76f];
	kTroeInf = 7.9230000000E-68	 * w[rLPDD76f];
	fcTroe = 1 * exp( -0 * RGAS* rt_inv );

	k[rLPDD76f] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM0] );
	kTroe0 = 5.1040000000E+02	 * k[rLPDD76b];
	kTroeInf = 1.2760000000E+03	 * w[rLPDD76b];

	k[rLPDD76b] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM0] );
	kTroe0 = 3.9900000000E+36	 * k[rG183];
	kTroeInf = 1.9510000000E+15	 * w[rG183];
	fcTroe = 1 * exp( -temp * 0.000763359 ) + 1 * exp( -48100 * RGAS* rt_inv );

	k[rG183] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM3] );
	kTroe0 = 1.6290000000E+67	 * k[rG163];
	kTroeInf = 1.3850000000E+13	 * w[rG163];
	fcTroe = 0.8106 * exp( -temp * 0.00361011 ) + 0.1894 * exp( -temp * 0.000114312 ) + 1 * exp( -7891 * RGAS* rt_inv );

	k[rG163] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM3] );
	kTroe0 = 4.6600000000E+35	 * k[rT23f];
	kTroeInf = 2.4300000000E+09	 * w[rT23f];
	fcTroe = 0.3 * exp( -temp * 0.01 ) + 0.7 * exp( -temp * 1.11111e-05 ) + 1 * exp( -10000 * RGAS* rt_inv );

	k[rT23f] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM3] );
	kTroe0 = 4.3784646437E+38	 * k[rT23b];
	kTroeInf = 2.2831907906E+12	 * w[rT23b];

	k[rT23b] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM3] );
	kTroe0 = 1.0000000000E+81	 * k[rOX01f];
	kTroeInf = 4.3000000000E+12	 * w[rOX01f];
	fcTroe = 0.098 * exp( -temp * 0.00143678 ) + 0.902 * exp( -temp * 0.0027933 ) + 1 * exp( -3856 * RGAS* rt_inv );

	k[rOX01f] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM3] );
	kTroe0 = 4.2026355221E+75	 * k[rOX01b];
	kTroeInf = 1.8071332745E+07	 * w[rOX01b];

	k[rOX01b] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM3] );
	kTroe0 = 6.3700000000E+11	 * k[rNO185f];
	kTroeInf = 7.9100000000E+10	 * w[rNO185f];
	fcTroe = 1 * exp( -0 * RGAS* rt_inv );

	k[rNO185f] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM32] );
	kTroe0 = 1.3071786932E+03	 * k[rNO185b];
	kTroeInf = 1.6231999157E+02	 * w[rNO185b];

	k[rNO185b] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM32] );
	kTroe0 = 1.4000000000E+20	 * k[rNO237f];
	kTroeInf = 3.3000000000E+10	 * w[rNO237f];
	fcTroe = 1 * exp( -0 * RGAS* rt_inv );

	k[rNO237f] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM32] );
	kTroe0 = 6.8832428931E+23	 * k[rNO237b];
	kTroeInf = 1.6224786819E+14	 * w[rNO237b];

	k[rNO237b] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM32] );
	kTroe0 = 1.3000000000E+19	 * k[rNO241f];
	kTroeInf = 3.1000000000E+09	 * w[rNO241f];
	fcTroe = 0.333 * exp( -temp * 0.00425532 ) + 0.667 * exp( -temp * 0.000472367 ) + 1 * exp( -4536 * RGAS* rt_inv );

	k[rNO241f] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM32] );
	kTroe0 = 8.2668184843E+28	 * k[rNO241b];
	kTroeInf = 1.9713182539E+19	 * w[rNO241b];

	k[rNO241b] = GetLindRateCoeff( temp, pressure,kTroe0 , kTroeInf
				, fcTroe, M[mM32] );


	w[r3f] = k[r3f] * c[sO] * c[sH2];
	w[r3b] = k[r3b] * c[sOH] * c[sH];
	w[r4f] = k[r4f] * c[sH2] * c[sOH];
	w[r4b] = k[r4b] * c[sH] * c[sH2O];
	w[r5f] = k[r5f] * c[sOH] * c[sOH];
	w[r5b] = k[r5b] * c[sH2O] * c[sO];
	w[r1f] = k[r1f] * c[sH] * c[sO2];
	w[r1b] = k[r1b] * c[sOH] * c[sO];
	w[r16] = k[r16] * c[sHO2] * c[sH];
	w[r21] = k[r21] * c[sHO2] * c[sHO2];
	w[r20] = k[r20] * c[sHO2] * c[sHO2];
	w[r19] = k[r19] * c[sHO2] * c[sOH];
	w[r18] = k[r18] * c[sHO2] * c[sO];
	w[r17] = k[r17] * c[sHO2] * c[sH];
	w[r26] = k[r26] * c[sH2O2] * c[sOH];
	w[r23] = k[r23] * c[sH2O2] * c[sH];
	w[r39] = k[r39] * c[sH2O2] * c[sOH];
	w[rG06f] = k[rG06f] * c[sCH] * c[sH2];
	w[rG06b] = k[rG06b] * c[sH] * c[sTXCH2];
	w[rG08] = k[rG08] * c[sCH] * c[sH2O];
	w[rG19] = k[rG19] * c[sTXCH2] * c[sO2];
	w[rG20] = k[rG20] * c[sTXCH2] * c[sO2];
	w[rG17] = k[rG17] * c[sTXCH2] * c[sOH];
	w[rG21] = k[rG21] * c[sTXCH2] * c[sO2];
	w[rG38b] = k[rG38b] * c[sH2O] * c[sTXCH2];
	w[rG34b] = k[rG34b] * c[sH] * c[sCH3];
	w[rG28f] = k[rG28f] * c[sSXCH2] * c[sN2];
	w[rG28b] = k[rG28b] * c[sN2] * c[sTXCH2];
	w[rG55] = k[rG55] * c[sCH3] * c[sOH];
	w[rG57f] = k[rG57f] * c[sCH3] * c[sOH];
	w[rG57b] = k[rG57b] * c[sH2O] * c[sSXCH2];
	w[rG59] = k[rG59] * c[sCH3] * c[sO2];
	w[rG74f] = k[rG74f] * c[sCH3] * c[sCH3];
	w[rG74b] = k[rG74b] * c[sH] * c[sC2H5];
	w[r29] = k[r29] * c[sCO] * c[sOH];
	w[r31] = k[r31] * c[sCO] * c[sHO2];
	w[r28f] = k[r28f] * c[sCO] * c[sOH];
	w[r28b] = k[r28b] * c[sH] * c[sCO2];
	w[r38] = k[r38] * c[sHCO] * c[sO2];
	w[r36f] = k[r36f] * c[sHCO] * M[mM21];
	w[r36b] = k[r36b] * c[sH] * c[sCO] * M[mM21];
	w[rDD68] = k[rDD68] * c[sHCO] * c[sO2];
	w[rG46] = k[rG46] * c[sCH2O] * c[sO];
	w[rG47] = k[rG47] * c[sCH2O] * c[sOH];
	w[rG45] = k[rG45] * c[sCH2O] * c[sH];
	w[rG49] = k[rG49] * c[sCH2O] * c[sHO2];
	w[rLPDD62] = k[rLPDD62] * c[sCH2O] * c[sOH];
	w[rLPG89] = k[rLPG89] * c[sCH3O] * c[sO2];
	w[rDD44] = k[rDD44] * c[sOCHO] * M[mM0];
	w[rDD47] = k[rDD47] * c[sHOCHO] * c[sOH];
	w[rDD50] = k[rDD50] * c[sHOCHO] * c[sH];
	w[rDD48] = k[rDD48] * c[sHOCHO] * c[sHO2];
	w[rDD53] = k[rDD53] * c[sHOCHO];
	w[rDD52] = k[rDD52] * c[sHOCHO] * c[sOH];
	w[rLPG63] = k[rLPG63] * c[sCH3O2] * c[sHO2];
	w[rLPDD59] = k[rLPDD59] * c[sCH3O2];
	w[rLPG61] = k[rLPG61] * c[sCH3O2] * c[sCH3];
	w[rLPG64] = k[rLPG64] * c[sCH3O2] * c[sCH2O];
	w[rDD66] = k[rDD66] * c[sCH2O] * c[sO2CHO];
	w[rDD67] = k[rDD67] * c[sO2CHO];
	w[rG116] = k[rG116] * c[sC2H2] * c[sO];
	w[rG117] = k[rG117] * c[sC2H2] * c[sO];
	w[rLPG120f] = k[rLPG120f] * c[sC2H2] * c[sOH];
	w[rLPG120b] = k[rLPG120b] * c[sH] * c[sCH2CO];
	w[rG122] = k[rG122] * c[sC2H2] * c[sOH];
	w[rG135] = k[rG135] * c[sC2H3] * c[sO2];
	w[rLPG134] = k[rLPG134] * c[sC2H3] * c[sO2];
	w[rG133] = k[rG133] * c[sC2H3] * c[sO2];
	w[rLPG131] = k[rLPG131] * c[sC2H3] * c[sO];
	w[rG160] = k[rG160] * c[sC2H4] * c[sOH];
	w[rG158] = k[rG158] * c[sC2H4] * c[sO];
	w[rG159] = k[rG159] * c[sC2H4] * c[sO];
	w[rLPG157] = k[rLPG157] * c[sC2H4] * c[sO];
	w[rG156] = k[rG156] * c[sC2H4] * c[sH];
	w[rG161] = k[rG161] * c[sC2H4] * c[sOH];
	w[rLPDD15] = k[rLPDD15] * c[sCH3O2] * c[sC2H5];
	w[rG171] = k[rG171] * c[sC2H5] * c[sO2];
	w[rG111] = k[rG111] * c[sHCCO] * c[sO2];
	w[rLPG124] = k[rLPG124] * c[sCH2CO] * c[sH];
	w[rLPR028] = k[rLPR028] * c[sCH2CO] * c[sCH3];
	w[rLPG123f] = k[rLPG123f] * c[sCH2CO] * c[sH];
	w[rLPG123b] = k[rLPG123b] * c[sH2] * c[sHCCO];
	w[rLPG127f] = k[rLPG127f] * c[sCH2CO] * c[sOH];
	w[rLPG127b] = k[rLPG127b] * c[sH2O] * c[sHCCO];
	w[rLPR025] = k[rLPR025] * c[sCH2CO] * c[sOH];
	w[rLPR030] = k[rLPR030] * c[sCH2CHO] * c[sCH3];
	w[rLPG139] = k[rLPG139] * c[sCH2CHO] * c[sO2];
	w[rLPG140] = k[rLPG140] * c[sCH2CHO] * c[sO2];
	w[rLPG137] = k[rLPG137] * c[sCH2CHO];
	w[rLPG138] = k[rLPG138] * c[sCH2CHO] * c[sO];
	w[rLPG136] = k[rLPG136] * c[sCH2CHO];
	w[rLPG144] = k[rLPG144] * c[sCH2CHO] * c[sOH];
	w[rG168] = k[rG168] * c[sC2H5O];
	w[rLPR015] = k[rLPR015] * c[sH] * c[sAXC3H5];
	w[rLPR116] = k[rLPR116] * c[sAXC3H5] * c[sHO2];
	w[rLPR115] = k[rLPR115] * c[sAXC3H5] * c[sHCO];
	w[rLPR133] = k[rLPR133] * c[sAXC3H5] * c[sO2];
	w[rLPLPR118] = k[rLPLPR118] * c[sAXC3H5] * c[sCH3O2];
	w[rLPR016] = k[rLPR016] * c[sAXC3H5] * c[sH];
	w[rLPLPR113] = k[rLPLPR113] * c[sAXC3H5] * c[sO];
	w[rLPR111] = k[rLPR111] * c[sAXC3H5] * c[sO2];
	w[rLPR007] = k[rLPR007] * c[sAXC3H5];
	w[rLPLPR117] = k[rLPLPR117] * c[sAXC3H5] * c[sHO2];
	w[rR141] = k[rR141] * c[sC3H6] * c[sH];
	w[rR143] = k[rR143] * c[sC3H6] * c[sO];
	w[rLPR153] = k[rLPR153] * c[sC3H6] * c[sH];
	w[rLPR142] = k[rLPR142] * c[sC3H6] * c[sO];
	w[rLPR151] = k[rLPR151] * c[sC3H6] * c[sOH];
	w[rG184] = k[rG184] * c[sNXC3H7] * c[sO2];
	w[rLPHP55] = k[rLPHP55] * c[sPXC4H8] * c[sOH];
	w[rLPHP43] = k[rLPHP43] * c[sPXC4H8];
	w[rLPHP54] = k[rLPHP54] * c[sPXC4H8] * c[sO];
	w[rLPHP52] = k[rLPHP52] * c[sPXC4H8] * c[sO];
	w[rLPHP53] = k[rLPHP53] * c[sPXC4H8] * c[sO];
	w[rLPDD514] = k[rLPDD514] * c[sPXC4H8] * c[sH];
	w[rLPHP56] = k[rLPHP56] * c[sPXC4H8] * c[sOH];
	w[rLPHP39] = k[rLPHP39] * c[sPXC4H9];
	w[rHP41] = k[rHP41] * c[sPXC4H9];
	w[rHP29] = k[rHP29] * c[sC5H11];
	w[rHP32] = k[rHP32] * c[sC5H11];
	w[rDD199] = k[rDD199] * c[sC6H12X1XC6H12] * c[sOH];
	w[rDD212] = k[rDD212] * c[sC6H12X1XC6H12] * c[sH];
	w[rLPDD205] = k[rLPDD205] * c[sC6H12X1XC6H12];
	w[rDD214] = k[rDD214] * c[sC6H13];
	w[rDD208] = k[rDD208] * c[sC6H13];
	w[rDD43] = k[rDD43] * c[sC7H14] * c[sOH];
	w[rHP26] = k[rHP26] * c[sC7H14] * c[sOH];
	w[rLPHP22] = k[rLPHP22] * c[sC7H14];
	w[rHP15] = k[rHP15] * c[sC7H14] * c[sH];
	w[rLPHP12] = k[rLPHP12] * c[sC7H15];
	w[rHP11] = k[rHP11] * c[sC7H15];
	w[rHP10] = k[rHP10] * c[sC7H15];
	w[rLPDD286] = k[rLPDD286] * c[sC8H17];
	w[rDD282] = k[rDD282] * c[sC8H17];
	w[rDD288] = k[rDD288] * c[sC8H17];
	w[rDD287] = k[rDD287] * c[sC8H17];
	w[rDD323] = k[rDD323] * c[sC9H19];
	w[rDD324] = k[rDD324] * c[sC9H19];
	w[rDD319] = k[rDD319] * c[sC9H19];
	w[rDD210] = k[rDD210] * c[sNXC12H26];
	w[rDD421] = k[rDD421] * c[sNXC12H26];
	w[rDD436] = k[rDD436] * c[sNXC12H26];
	w[rDD425] = k[rDD425] * c[sNXC12H26];
	w[rDD422] = k[rDD422] * c[sNXC12H26] * c[sO2];
	w[rDD434] = k[rDD434] * c[sNXC12H26] * c[sOH];
	w[rDD429] = k[rDD429] * c[sNXC12H26] * c[sHO2];
	w[rDD424] = k[rDD424] * c[sNXC12H26] * c[sO];
	w[rDD432] = k[rDD432] * c[sNXC12H26] * c[sH];
	w[rDD427] = k[rDD427] * c[sNXC12H26] * c[sO2CHO];
	w[rLPDD431] = k[rLPDD431] * c[sNXC12H26] * c[sCH3O2];
	w[rDD405f] = k[rDD405f] * c[sC12H25];
	w[rDD405b] = k[rDD405b] * c[sH] * c[sC12H24];
	w[rDD418] = k[rDD418] * c[sC12H25];
	w[rDD415] = k[rDD415] * c[sC12H25];
	w[rDD416] = k[rDD416] * c[sC12H25];
	w[rLPDD406] = k[rLPDD406] * c[sC12H25];
	w[rDD408f] = k[rDD408f] * c[sC12H25] * c[sO2];
	w[rDD408b] = k[rDD408b] * c[sC12H25O2];
	w[rDD472] = k[rDD472] * c[sC12H25O2];
	w[rLPDD465f] = k[rLPDD465f] * c[sC12H25O2];
	w[rLPDD465b] = k[rLPDD465b] * c[sC12OOHXT6XC12H25O2];
	w[rLPDD473] = k[rLPDD473] * c[sC12OOHXT6XC12H25O2];
	w[rLPDD486f] = k[rLPDD486f] * c[sC12OOHXT6XC12H25O2] * c[sO2];
	w[rLPDD486b] = k[rLPDD486b] * c[sC12OOHO2XT6XC12H25O4];
	w[rLPDD502] = k[rLPDD502] * c[sC12OOHO2XT6XC12H25O4];
	w[rLPLPDD495] = k[rLPLPDD495] * c[sC12KETXT6XC12H24O3];
	w[rLPDD499] = k[rLPDD499] * c[sC12KETXT6XC12H24O3];
	w[rLPDD490] = k[rLPDD490] * c[sC12KETXT6XC12H24O3];
	w[rLPR140] = k[rLPR140] * c[sCH3COCH2];
	w[rLPLPDD107] = k[rLPLPDD107] * c[sCH3COCH2];
	w[rDD332N] = k[rDD332N] * c[sNXC8H17CHO] * c[sOH];
	w[rDD294] = k[rDD294] * c[sNXC7H15CHO] * c[sOH];
	w[rDD290] = k[rDD290] * c[sC7H15CO];
	w[rLPDD293] = k[rLPDD293] * c[sC7H15CO];
	w[rHP60] = k[rHP60] * c[sC3H7CHO] * c[sOH];
	w[rDD239] = k[rDD239] * c[sC6H13] * c[sO2];
	w[rDD263] = k[rDD263] * c[sC7H15] * c[sO2];
	w[rDD260f] = k[rDD260f] * c[sC7H15O2X1XC7H15O2];
	w[rDD260b] = k[rDD260b] * c[sC7H14OOH1X3XC7H15O2];
	w[rDD259] = k[rDD259] * c[sC7H14OOH1X3XC7H15O2];
	w[rP000f] = k[rP000f] * c[sC5H4CH2] * c[sH];
	w[rP000b] = k[rP000b] * c[sH] * c[sA1XC6H6];
	w[rP001f] = k[rP001f] * c[sNXC4H5] * c[sC2H2];
	w[rP001b] = k[rP001b] * c[sH] * c[sC5H4CH2];
	w[rP002f] = k[rP002f] * c[sIXC4H5] * c[sC2H2];
	w[rP002b] = k[rP002b] * c[sH] * c[sC5H4CH2];
	w[rP003f] = k[rP003f] * c[sC5H4CH2];
	w[rP003b] = k[rP003b] * c[sA1XC6H6];
	w[rP004f] = k[rP004f] * c[sC5H4CH2];
	w[rP004b] = k[rP004b] * c[sH] * c[sA1XXC6H5];
	w[rP005f] = k[rP005f] * c[sNXC4H3] * c[sC2H2];
	w[rP005b] = k[rP005b] * c[sA1XXC6H5];
	w[rP006f] = k[rP006f] * c[sNXC4H5] * c[sC2H2];
	w[rP006b] = k[rP006b] * c[sH] * c[sA1XC6H6];
	w[rP007f] = k[rP007f] * c[sIXC4H5] * c[sC2H2];
	w[rP007b] = k[rP007b] * c[sH] * c[sA1XC6H6];
	w[rP008] = k[rP008] * c[sAXC3H5] * c[sC3H3];
	w[rP009f] = k[rP009f] * c[sC3H3] * c[sC3H3];
	w[rP009b] = k[rP009b] * c[sC5H4CH2];
	w[rP010f] = k[rP010f] * c[sC3H3] * c[sC3H3];
	w[rP010b] = k[rP010b] * c[sA1XC6H6];
	w[rP011f] = k[rP011f] * c[sC3H3] * c[sC3H3];
	w[rP011b] = k[rP011b] * c[sH] * c[sA1XXC6H5];
	w[rK012f] = k[rK012f] * c[sA1XXC6H5] * c[sC2H2];
	w[rK012b] = k[rK012b] * c[sA1C2H2XC8H7];
	w[rP013b] = k[rP013b] * c[sA1C2H3XC8H8];
	w[rP014f] = k[rP014f] * c[sA1XC6H6] * c[sC2H3];
	w[rP014b] = k[rP014b] * c[sH] * c[sA1C2H3XC8H8];
	w[rP015f] = k[rP015f] * c[sA1XXC6H5] * c[sC2H4];
	w[rP015b] = k[rP015b] * c[sC2H3] * c[sA1XC6H6];
	w[rP016f] = k[rP016f] * c[sA1C2HXC8H6];
	w[rP016b] = k[rP016b] * c[sH] * c[sA1C2HYXC8H5];
	w[rK017f] = k[rK017f] * c[sA1C2HXC8H6] * c[sH];
	w[rK017b] = k[rK017b] * c[sH2] * c[sA1C2HYXC8H5];
	w[rP018f] = k[rP018f] * c[sA1C2HXC8H6] * c[sOH];
	w[rP018b] = k[rP018b] * c[sH2O] * c[sA1C2HYXC8H5];
	w[rK019f] = k[rK019f] * c[sA1C2H2XC8H7];
	w[rK019b] = k[rK019b] * c[sA1C2H3YXC8H7];
	w[rK020f] = k[rK020f] * c[sA1C2H2XC8H7];
	w[rK020b] = k[rK020b] * c[sH] * c[sA1C2HXC8H6];
	w[rK021f] = k[rK021f] * c[sA1C2H2XC8H7] * c[sH];
	w[rK021b] = k[rK021b] * c[sH2] * c[sA1C2HXC8H6];
	w[rP022b] = k[rP022b] * c[sH2O] * c[sA1C2HXC8H6];
	w[rP023f] = k[rP023f] * c[sA1C2H3XC8H8];
	w[rP023b] = k[rP023b] * c[sH] * c[sA1C2H3YXC8H7];
	w[rK024f] = k[rK024f] * c[sA1C2H3XC8H8] * c[sH];
	w[rK024b] = k[rK024b] * c[sH2] * c[sA1C2H3YXC8H7];
	w[rP025f] = k[rP025f] * c[sA1C2H3XC8H8] * c[sOH];
	w[rP025b] = k[rP025b] * c[sH2O] * c[sA1C2H3YXC8H7];
	w[rP026f] = k[rP026f] * c[sA1C2H3XC8H8];
	w[rP026b] = k[rP026b] * c[sH] * c[sA1C2H2XC8H7];
	w[rP027f] = k[rP027f] * c[sA1C2H3XC8H8] * c[sH];
	w[rP027b] = k[rP027b] * c[sH2] * c[sA1C2H2XC8H7];
	w[rP028f] = k[rP028f] * c[sA1C2H3XC8H8] * c[sOH];
	w[rP028b] = k[rP028b] * c[sH2O] * c[sA1C2H2XC8H7];
	w[rK100f] = k[rK100f] * c[sA1C2HYXC8H5] * c[sC2H2];
	w[rK100b] = k[rK100b] * c[sA2XXC10H7];
	w[rK102f] = k[rK102f] * c[sA1C2H3YXC8H7] * c[sC2H2];
	w[rK102b] = k[rK102b] * c[sH] * c[sA2XC10H8];
	w[rP104f] = k[rP104f] * c[sA1C2H2XC8H7] * c[sC2H2];
	w[rP104b] = k[rP104b] * c[sH] * c[sA2XC10H8];
	w[rP105f] = k[rP105f] * c[sA1C2HXC8H6] * c[sC2H3];
	w[rP105b] = k[rP105b] * c[sH] * c[sA2XC10H8];
	w[rP106f] = k[rP106f] * c[sA1C2HYXC8H5] * c[sC2H4];
	w[rP106b] = k[rP106b] * c[sH] * c[sA2XC10H8];
	w[rP107f] = k[rP107f] * c[sA1XXC6H5] * c[sC4H4];
	w[rP107b] = k[rP107b] * c[sH] * c[sA2XC10H8];
	w[rP108f] = k[rP108f] * c[sA2XC10H8];
	w[rP108b] = k[rP108b] * c[sH] * c[sA2XXC10H7];
	w[rK109f] = k[rK109f] * c[sA2XC10H8] * c[sH];
	w[rK109b] = k[rK109b] * c[sH2] * c[sA2XXC10H7];
	w[rK110f] = k[rK110f] * c[sA2XC10H8] * c[sOH];
	w[rK110b] = k[rK110b] * c[sH2O] * c[sA2XXC10H7];
	w[rP111f] = k[rP111f] * c[sA2XC10H8];
	w[rP111b] = k[rP111b] * c[sH] * c[sA2YXC10H7];
	w[rK112f] = k[rK112f] * c[sA2XC10H8] * c[sH];
	w[rK112b] = k[rK112b] * c[sH2] * c[sA2YXC10H7];
	w[rK113f] = k[rK113f] * c[sA2XC10H8] * c[sOH];
	w[rK113b] = k[rK113b] * c[sH2O] * c[sA2YXC10H7];
	w[rK114f] = k[rK114f] * c[sA2XXC10H7] * c[sC2H2];
	w[rK114b] = k[rK114b] * c[sA2C2H2AXC12H9];
	w[rK115f] = k[rK115f] * c[sA2YXC10H7] * c[sC2H2];
	w[rK115b] = k[rK115b] * c[sA2C2H2BXC12H9];
	w[rP116b] = k[rP116b] * c[sH] * c[sA2C2H2AXC12H9];
	w[rP117b] = k[rP117b] * c[sH] * c[sA2C2H2BXC12H9];
	w[rP118f] = k[rP118f] * c[sA2XC10H8] * c[sC2H3];
	w[rP118b] = k[rP118b] * c[sH2] * c[sA2C2H2AXC12H9];
	w[rP119f] = k[rP119f] * c[sA2XC10H8] * c[sC2H3];
	w[rP119b] = k[rP119b] * c[sH2] * c[sA2C2H2BXC12H9];
	w[rP120f] = k[rP120f] * c[sA2XXC10H7] * c[sC2H4];
	w[rP120b] = k[rP120b] * c[sH2] * c[sA2C2H2AXC12H9];
	w[rP121f] = k[rP121f] * c[sA2YXC10H7] * c[sC2H4];
	w[rP121b] = k[rP121b] * c[sH2] * c[sA2C2H2BXC12H9];
	w[rK122f] = k[rK122f] * c[sA2C2H2AXC12H9];
	w[rK122b] = k[rK122b] * c[sH] * c[sA2C2HAXC12H8];
	w[rK123f] = k[rK123f] * c[sA2C2H2AXC12H9] * c[sH];
	w[rK123b] = k[rK123b] * c[sH2] * c[sA2C2HAXC12H8];
	w[rP124b] = k[rP124b] * c[sH2O] * c[sA2C2HAXC12H8];
	w[rK125f] = k[rK125f] * c[sA2C2H2BXC12H9];
	w[rK125b] = k[rK125b] * c[sH] * c[sA2C2HBXC12H8];
	w[rK126f] = k[rK126f] * c[sA2C2H2BXC12H9] * c[sH];
	w[rK126b] = k[rK126b] * c[sH2] * c[sA2C2HBXC12H8];
	w[rP127b] = k[rP127b] * c[sH2O] * c[sA2C2HBXC12H8];
	w[rP128f] = k[rP128f] * c[sA2C2HAXC12H8];
	w[rP128b] = k[rP128b] * c[sH] * c[sA2C2HAYXC12H7];
	w[rK129f] = k[rK129f] * c[sA2C2HAXC12H8] * c[sH];
	w[rK129b] = k[rK129b] * c[sH2] * c[sA2C2HAYXC12H7];
	w[rP130f] = k[rP130f] * c[sA2C2HAXC12H8] * c[sOH];
	w[rP130b] = k[rP130b] * c[sH2O] * c[sA2C2HAYXC12H7];
	w[rP131f] = k[rP131f] * c[sA2C2HBXC12H8];
	w[rP131b] = k[rP131b] * c[sH] * c[sA2C2HBYXC12H7];
	w[rK132f] = k[rK132f] * c[sA2C2HBXC12H8] * c[sH];
	w[rK132b] = k[rK132b] * c[sH2] * c[sA2C2HBYXC12H7];
	w[rP133f] = k[rP133f] * c[sA2C2HBXC12H8] * c[sOH];
	w[rP133b] = k[rP133b] * c[sH2O] * c[sA2C2HBYXC12H7];
	w[rK200f] = k[rK200f] * c[sA2C2H2AXC12H9];
	w[rK200b] = k[rK200b] * c[sH] * c[sA2R5XC12H8];
	w[rP201f] = k[rP201f] * c[sA2C2HAXC12H8] * c[sH];
	w[rP201b] = k[rP201b] * c[sH] * c[sA2R5XC12H8];
	w[rP202f] = k[rP202f] * c[sA2R5XC12H8];
	w[rP202b] = k[rP202b] * c[sH] * c[sA2R5XXC12H7];
	w[rK203f] = k[rK203f] * c[sA2R5XC12H8] * c[sH];
	w[rK203b] = k[rK203b] * c[sH2] * c[sA2R5XXC12H7];
	w[rK204f] = k[rK204f] * c[sA2R5XC12H8] * c[sOH];
	w[rK204b] = k[rK204b] * c[sH2O] * c[sA2R5XXC12H7];
	w[rK205f] = k[rK205f] * c[sA2R5XXC12H7] * c[sC2H2];
	w[rK205b] = k[rK205b] * c[sA2R5C2H2XC14H9];
	w[rP206b] = k[rP206b] * c[sH] * c[sA2R5C2H2XC14H9];
	w[rP207f] = k[rP207f] * c[sA2R5XC12H8] * c[sC2H3];
	w[rP207b] = k[rP207b] * c[sH2] * c[sA2R5C2H2XC14H9];
	w[rP208f] = k[rP208f] * c[sA2R5XXC12H7] * c[sC2H4];
	w[rP208b] = k[rP208b] * c[sH2] * c[sA2R5C2H2XC14H9];
	w[rP209f] = k[rP209f] * c[sA2R5C2HXC14H8];
	w[rP209b] = k[rP209b] * c[sH] * c[sA2R5C2HYXC14H7];
	w[rK210f] = k[rK210f] * c[sA2R5C2HXC14H8] * c[sH];
	w[rK210b] = k[rK210b] * c[sH2] * c[sA2R5C2HYXC14H7];
	w[rP211f] = k[rP211f] * c[sA2R5C2HXC14H8] * c[sOH];
	w[rP211b] = k[rP211b] * c[sH2O] * c[sA2R5C2HYXC14H7];
	w[rK212f] = k[rK212f] * c[sA2R5C2H2XC14H9];
	w[rK212b] = k[rK212b] * c[sH] * c[sA2R5C2HXC14H8];
	w[rK213f] = k[rK213f] * c[sA2R5C2H2XC14H9] * c[sH];
	w[rK213b] = k[rK213b] * c[sH2] * c[sA2R5C2HXC14H8];
	w[rP214b] = k[rP214b] * c[sH2O] * c[sA2R5C2HXC14H8];
	w[rP301f] = k[rP301f] * c[sA1XC6H6] * c[sA1XXC6H5];
	w[rP301b] = k[rP301b] * c[sH] * c[sP2XC12H10];
	w[rP302f] = k[rP302f] * c[sA1XXC6H5] * c[sA1XXC6H5];
	w[rP302b] = k[rP302b] * c[sP2XC12H10];
	w[rP304f] = k[rP304f] * c[sP2XC12H10];
	w[rP304b] = k[rP304b] * c[sH] * c[sP2XXC12H9];
	w[rP305f] = k[rP305f] * c[sP2XC12H10] * c[sH];
	w[rP305b] = k[rP305b] * c[sH2] * c[sP2XXC12H9];
	w[rP306f] = k[rP306f] * c[sP2XC12H10] * c[sOH];
	w[rP306b] = k[rP306b] * c[sH2O] * c[sP2XXC12H9];
	w[rK401f] = k[rK401f] * c[sA2C2HAYXC12H7] * c[sC2H2];
	w[rK401b] = k[rK401b] * c[sA3XXC14H9];
	w[rK403f] = k[rK403f] * c[sA2C2HBYXC12H7] * c[sC2H2];
	w[rK403b] = k[rK403b] * c[sA3XXC14H9];
	w[rP405f] = k[rP405f] * c[sA2C2H2AXC12H9] * c[sC2H2];
	w[rP405b] = k[rP405b] * c[sH] * c[sA3XC14H10];
	w[rP406f] = k[rP406f] * c[sA2C2H2BXC12H9] * c[sC2H2];
	w[rP406b] = k[rP406b] * c[sH] * c[sA3XC14H10];
	w[rP407f] = k[rP407f] * c[sP2XXC12H9] * c[sC2H2];
	w[rP407b] = k[rP407b] * c[sH] * c[sA3XC14H10];
	w[rP408b] = k[rP408b] * c[sA3XC14H10];
	w[rP409b] = k[rP409b] * c[sA3XC14H10];
	w[rP410f] = k[rP410f] * c[sA2C2HAXC12H8] * c[sC2H3];
	w[rP410b] = k[rP410b] * c[sH] * c[sA3XC14H10];
	w[rP411f] = k[rP411f] * c[sA2C2HBXC12H8] * c[sC2H3];
	w[rP411b] = k[rP411b] * c[sH] * c[sA3XC14H10];
	w[rP412f] = k[rP412f] * c[sA2C2HAYXC12H7] * c[sC2H4];
	w[rP412b] = k[rP412b] * c[sH] * c[sA3XC14H10];
	w[rP413f] = k[rP413f] * c[sA2C2HBYXC12H7] * c[sC2H4];
	w[rP413b] = k[rP413b] * c[sH] * c[sA3XC14H10];
	w[rP414f] = k[rP414f] * c[sA2XXC10H7] * c[sC4H4];
	w[rP414b] = k[rP414b] * c[sH] * c[sA3XC14H10];
	w[rP415f] = k[rP415f] * c[sA2YXC10H7] * c[sC4H4];
	w[rP415b] = k[rP415b] * c[sH] * c[sA3XC14H10];
	w[rP416f] = k[rP416f] * c[sA1C2HXC8H6] * c[sA1XXC6H5];
	w[rP416b] = k[rP416b] * c[sH] * c[sA3XC14H10];
	w[rP417f] = k[rP417f] * c[sA1C2HYXC8H5] * c[sA1XC6H6];
	w[rP417b] = k[rP417b] * c[sH] * c[sA3XC14H10];
	w[rP418f] = k[rP418f] * c[sA1C2HYXC8H5] * c[sA1XXC6H5];
	w[rP418b] = k[rP418b] * c[sA3XC14H10];
	w[rP419f] = k[rP419f] * c[sA3XC14H10];
	w[rP419b] = k[rP419b] * c[sH] * c[sA3XXC14H9];
	w[rK420f] = k[rK420f] * c[sA3XC14H10] * c[sH];
	w[rK420b] = k[rK420b] * c[sH2] * c[sA3XXC14H9];
	w[rK421f] = k[rK421f] * c[sA3XC14H10] * c[sOH];
	w[rK421b] = k[rK421b] * c[sH2O] * c[sA3XXC14H9];
	w[rP422f] = k[rP422f] * c[sA3XC14H10];
	w[rP422b] = k[rP422b] * c[sH] * c[sA3YXC14H9];
	w[rK423f] = k[rK423f] * c[sA3XC14H10] * c[sH];
	w[rK423b] = k[rK423b] * c[sH2] * c[sA3YXC14H9];
	w[rK424f] = k[rK424f] * c[sA3XC14H10] * c[sOH];
	w[rK424b] = k[rK424b] * c[sH2O] * c[sA3YXC14H9];
	w[rP425f] = k[rP425f] * c[sA3XXC14H9];
	w[rP425b] = k[rP425b] * c[sC2H2] * c[sA2R5XXC12H7];
	w[rP501f] = k[rP501f] * c[sA2R5C2HYXC14H7] * c[sC2H2];
	w[rP501b] = k[rP501b] * c[sA3R5XXC16H9];
	w[rP502f] = k[rP502f] * c[sA2R5C2H2XC14H9] * c[sC2H2];
	w[rP502b] = k[rP502b] * c[sH] * c[sA3R5XC16H10];
	w[rP503f] = k[rP503f] * c[sA3YXC14H9] * c[sC2H2];
	w[rP503b] = k[rP503b] * c[sH] * c[sA3R5XC16H10];
	w[rP504b] = k[rP504b] * c[sH2] * c[sA3R5XC16H10];
	w[rP505f] = k[rP505f] * c[sA2R5XXC12H7] * c[sC4H4];
	w[rP505b] = k[rP505b] * c[sH] * c[sA3R5XC16H10];
	w[rP506f] = k[rP506f] * c[sA3R5XC16H10];
	w[rP506b] = k[rP506b] * c[sH] * c[sA3R5XXC16H9];
	w[rK507f] = k[rK507f] * c[sA3R5XC16H10] * c[sH];
	w[rK507b] = k[rK507b] * c[sH2] * c[sA3R5XXC16H9];
	w[rK508f] = k[rK508f] * c[sA3R5XC16H10] * c[sOH];
	w[rK508b] = k[rK508b] * c[sH2O] * c[sA3R5XXC16H9];
	w[rK600f] = k[rK600f] * c[sA3XXC14H9] * c[sC2H2];
	w[rK600b] = k[rK600b] * c[sH] * c[sA4XC16H10];
	w[rP601f] = k[rP601f] * c[sA4XC16H10];
	w[rP601b] = k[rP601b] * c[sH] * c[sA4XXC16H9];
	w[rK602f] = k[rK602f] * c[sA4XC16H10] * c[sH];
	w[rK602b] = k[rK602b] * c[sH2] * c[sA4XXC16H9];
	w[rK603f] = k[rK603f] * c[sA4XC16H10] * c[sOH];
	w[rK603b] = k[rK603b] * c[sH2O] * c[sA4XXC16H9];
	w[rK700f] = k[rK700f] * c[sA4XXC16H9] * c[sC2H2];
	w[rK700b] = k[rK700b] * c[sH] * c[sA4R5XC18H10];
	w[rK701f] = k[rK701f] * c[sA3R5XXC16H9] * c[sC2H2];
	w[rK701b] = k[rK701b] * c[sH] * c[sA4R5XC18H10];
	w[rP800] = k[rP800] * c[sA2XC10H8] * c[sA1XXC6H5];
	w[rP801] = k[rP801] * c[sA2XXC10H7] * c[sA1XC6H6];
	w[rP802] = k[rP802] * c[sA2XXC10H7] * c[sA1XXC6H5];
	w[rCP01f] = k[rCP01f] * c[sC5H6];
	w[rCP01b] = k[rCP01b] * c[sH] * c[sC5H5];
	w[rCP02f] = k[rCP02f] * c[sC5H6] * c[sH];
	w[rCP02b] = k[rCP02b] * c[sH2] * c[sC5H5];
	w[rCP03f] = k[rCP03f] * c[sC5H6] * c[sH];
	w[rCP03b] = k[rCP03b] * c[sC2H2] * c[sAXC3H5];
	w[rCP04] = k[rCP04] * c[sC5H6] * c[sH];
	w[rCP05f] = k[rCP05f] * c[sC5H6] * c[sO];
	w[rCP05b] = k[rCP05b] * c[sOH] * c[sC5H5];
	w[rCP06f] = k[rCP06f] * c[sC5H6] * c[sOH];
	w[rCP06b] = k[rCP06b] * c[sH2O] * c[sC5H5];
	w[rCP07f] = k[rCP07f] * c[sC5H6] * c[sO2];
	w[rCP07b] = k[rCP07b] * c[sHO2] * c[sC5H5];
	w[rCP08f] = k[rCP08f] * c[sC5H6] * c[sHO2];
	w[rCP08b] = k[rCP08b] * c[sH2O2] * c[sC5H5];
	w[rCP09f] = k[rCP09f] * c[sC5H6] * c[sCH3];
	w[rCP09b] = k[rCP09b] * c[sCH4] * c[sC5H5];
	w[rCP10b] = k[rCP10b] * c[sC2H4] * c[sC5H5];
	w[rCP11b] = k[rCP11b] * c[sC4H6] * c[sC5H5];
	w[rCP12f] = k[rCP12f] * c[sC5H6] * c[sO];
	w[rCP12b] = k[rCP12b] * c[sH] * c[sTXC5H5O];
	w[rCP13] = k[rCP13] * c[sC5H6] * c[sO];
	w[rCP14] = k[rCP14] * c[sC5H6] * c[sOH];
	w[rCP15] = k[rCP15] * c[sC5H6] * c[sOH];
	w[rCP16f] = k[rCP16f] * c[sC3H3] * c[sC2H2];
	w[rCP16b] = k[rCP16b] * c[sC5H5];
	w[rCP17] = k[rCP17] * c[sC5H5] * c[sC5H5];
	w[rCP18] = k[rCP18] * c[sC5H5] * c[sCH3];
	w[rCP19b] = k[rCP19b] * c[sH] * c[sC5H4O];
	w[rCP20f] = k[rCP20f] * c[sC5H5] * c[sO2];
	w[rCP20b] = k[rCP20b] * c[sOH] * c[sC5H4O];
	w[rCP21b] = k[rCP21b] * c[sOH] * c[sSXC5H5O];
	w[rCP23f] = k[rCP23f] * c[sSXC5H5O];
	w[rCP23b] = k[rCP23b] * c[sH] * c[sC5H4O];
	w[rCP24] = k[rCP24] * c[sTXC5H5O];
	w[rCP28] = k[rCP28] * c[sSXC5H5O] * c[sO2];
	w[rCP30] = k[rCP30] * c[sTXC5H5O] * c[sO2];
	w[rCP31] = k[rCP31] * c[sC5H4O];
	w[rCP32f] = k[rCP32f] * c[sC5H4O] * c[sH];
	w[rCP32b] = k[rCP32b] * c[sTXC5H5O];
	w[rCP33f] = k[rCP33f] * c[sC5H4O] * c[sO];
	w[rCP33b] = k[rCP33b] * c[sCO2] * c[sC4H4];
	w[rCP34] = k[rCP34] * c[sC5H5] * c[sC5H6];
	w[rI00b] = k[rI00b] * c[sC9H8];
	w[rI01f] = k[rI01f] * c[sC9H8];
	w[rI01b] = k[rI01b] * c[sH] * c[sC9H7];
	w[rI02f] = k[rI02f] * c[sC9H8] * c[sH];
	w[rI02b] = k[rI02b] * c[sH2] * c[sC9H7];
	w[rI03f] = k[rI03f] * c[sA1CH2XC7H7] * c[sC2H2];
	w[rI03b] = k[rI03b] * c[sH] * c[sC9H8];
	w[rI05f] = k[rI05f] * c[sC9H8] * c[sO];
	w[rI05b] = k[rI05b] * c[sOH] * c[sC9H7];
	w[rI06f] = k[rI06f] * c[sC9H8] * c[sOH];
	w[rI06b] = k[rI06b] * c[sH2O] * c[sC9H7];
	w[rI07f] = k[rI07f] * c[sC9H8] * c[sO2];
	w[rI07b] = k[rI07b] * c[sHO2] * c[sC9H7];
	w[rI08f] = k[rI08f] * c[sC9H8] * c[sHO2];
	w[rI08b] = k[rI08b] * c[sH2O2] * c[sC9H7];
	w[rI09f] = k[rI09f] * c[sC9H8] * c[sCH3];
	w[rI09b] = k[rI09b] * c[sCH4] * c[sC9H7];
	w[rI12] = k[rI12] * c[sC9H8] * c[sO];
	w[rI15] = k[rI15] * c[sC9H8] * c[sOH];
	w[rI17] = k[rI17] * c[sC9H7] * c[sC5H5];
	w[rI18] = k[rI18] * c[sC9H7] * c[sCH3];
	w[rI19b] = k[rI19b] * c[sH] * c[sC9H6O];
	w[rI20f] = k[rI20f] * c[sC9H7] * c[sO2];
	w[rI20b] = k[rI20b] * c[sOH] * c[sC9H6O];
	w[rI23] = k[rI23] * c[sC9H7] * c[sC3H3];
	w[rI31] = k[rI31] * c[sC9H6O];
	w[rI32] = k[rI32] * c[sC9H6O] * c[sH];
	w[rT01f] = k[rT01f] * c[sA1CH3XC7H8] * c[sH];
	w[rT01b] = k[rT01b] * c[sCH3] * c[sA1XC6H6];
	w[rT02f] = k[rT02f] * c[sA1CH3XC7H8];
	w[rT02b] = k[rT02b] * c[sH] * c[sA1CH2XC7H7];
	w[rT03f] = k[rT03f] * c[sA1CH3XC7H8];
	w[rT03b] = k[rT03b] * c[sCH3] * c[sA1XXC6H5];
	w[rT04f] = k[rT04f] * c[sA1CH2XC7H7] * c[sH];
	w[rT04b] = k[rT04b] * c[sCH3] * c[sA1XXC6H5];
	w[rT05f] = k[rT05f] * c[sA1CH2XC7H7];
	w[rT05b] = k[rT05b] * c[sC2H2] * c[sC5H5];
	w[rT06f] = k[rT06f] * c[sA1CH3XC7H8] * c[sO2];
	w[rT06b] = k[rT06b] * c[sHO2] * c[sA1CH2XC7H7];
	w[rT07f] = k[rT07f] * c[sA1CH3XC7H8] * c[sH];
	w[rT07b] = k[rT07b] * c[sH2] * c[sA1CH2XC7H7];
	w[rT08f] = k[rT08f] * c[sA1CH3XC7H8] * c[sOH];
	w[rT08b] = k[rT08b] * c[sH2O] * c[sA1CH2XC7H7];
	w[rT09f] = k[rT09f] * c[sA1CH3XC7H8] * c[sOH];
	w[rT09b] = k[rT09b] * c[sCH3] * c[sA1OHXC6H6O];
	w[rT10f] = k[rT10f] * c[sA1CH3XC7H8] * c[sOH];
	w[rT10b] = k[rT10b] * c[sH] * c[sHOA1CH3XC7H8O];
	w[rT11f] = k[rT11f] * c[sA1CH3XC7H8] * c[sO];
	w[rT11b] = k[rT11b] * c[sOH] * c[sA1CH2XC7H7];
	w[rT12f] = k[rT12f] * c[sA1CH3XC7H8] * c[sO];
	w[rT12b] = k[rT12b] * c[sHOA1CH3XC7H8O];
	w[rT13f] = k[rT13f] * c[sA1CH3XC7H8] * c[sO];
	w[rT13b] = k[rT13b] * c[sH] * c[sOA1CH3XC7H7O];
	w[rT14f] = k[rT14f] * c[sA1CH3XC7H8] * c[sCH3];
	w[rT14b] = k[rT14b] * c[sCH4] * c[sA1CH2XC7H7];
	w[rT15f] = k[rT15f] * c[sA1CH3XC7H8] * c[sHO2];
	w[rT15b] = k[rT15b] * c[sH2O2] * c[sA1CH2XC7H7];
	w[rT16f] = k[rT16f] * c[sA1CH3XC7H8] * c[sA1XXC6H5];
	w[rT16b] = k[rT16b] * c[sA1XC6H6] * c[sA1CH2XC7H7];
	w[rT17b] = k[rT17b] * c[sA1CH2OXC7H7O];
	w[rT18b] = k[rT18b] * c[sA1CH2OHXC7H8O];
	w[rT19f] = k[rT19f] * c[sA1CH2XC7H7] * c[sHO2];
	w[rT19b] = k[rT19b] * c[sOH] * c[sA1CH2OXC7H7O];
	w[rT20] = k[rT20] * c[sA1CH2XC7H7] * c[sC3H3];
	w[rT21f] = k[rT21f] * c[sA1CH2XC7H7] * c[sO2];
	w[rT21b] = k[rT21b] * c[sOH] * c[sA1CHOXC7H6O];
	w[rT22f] = k[rT22f] * c[sA1CH2XC7H7] * c[sO2];
	w[rT22b] = k[rT22b] * c[sCH2O] * c[sA1OXC6H5O];
	w[rT24f] = k[rT24f] * c[sA1CH2OHXC7H8O] * c[sH];
	w[rT24b] = k[rT24b] * c[sH2] * c[sA1CH2OXC7H7O];
	w[rT25f] = k[rT25f] * c[sA1CH2OHXC7H8O] * c[sO];
	w[rT25b] = k[rT25b] * c[sOH] * c[sA1CH2OXC7H7O];
	w[rT26f] = k[rT26f] * c[sA1CH2OHXC7H8O] * c[sOH];
	w[rT26b] = k[rT26b] * c[sH2O] * c[sA1CH2OXC7H7O];
	w[rT27f] = k[rT27f] * c[sA1CH2OHXC7H8O] * c[sCH3];
	w[rT27b] = k[rT27b] * c[sCH4] * c[sA1CH2OXC7H7O];
	w[rT28f] = k[rT28f] * c[sA1XXC6H5] * c[sCH2O];
	w[rT28b] = k[rT28b] * c[sHCO] * c[sA1XC6H6];
	w[rT29f] = k[rT29f] * c[sA1CH2OXC7H7O];
	w[rT29b] = k[rT29b] * c[sH] * c[sA1CHOXC7H6O];
	w[rT30f] = k[rT30f] * c[sA1CH2OXC7H7O];
	w[rT30b] = k[rT30b] * c[sCH2O] * c[sA1XXC6H5];
	w[rT31f] = k[rT31f] * c[sA1CH2OXC7H7O];
	w[rT31b] = k[rT31b] * c[sHCO] * c[sA1XC6H6];
	w[rT32b] = k[rT32b] * c[sH2] * c[sA1CHOXC7H6O];
	w[rT33b] = k[rT33b] * c[sOH] * c[sA1CHOXC7H6O];
	w[rT34b] = k[rT34b] * c[sH2O] * c[sA1CHOXC7H6O];
	w[rT35f] = k[rT35f] * c[sA1CH2OXC7H7O] * c[sO2];
	w[rT35b] = k[rT35b] * c[sHO2] * c[sA1CHOXC7H6O];
	w[rT36] = k[rT36] * c[sA1CHOXC7H6O];
	w[rT37] = k[rT37] * c[sA1CHOXC7H6O] * c[sH];
	w[rT38] = k[rT38] * c[sA1CHOXC7H6O] * c[sO];
	w[rT39] = k[rT39] * c[sA1CHOXC7H6O] * c[sOH];
	w[rT40] = k[rT40] * c[sA1CHOXC7H6O] * c[sO2];
	w[rT41] = k[rT41] * c[sA1CHOXC7H6O] * c[sHO2];
	w[rT42] = k[rT42] * c[sA1CHOXC7H6O] * c[sCH3];
	w[rT43f] = k[rT43f] * c[sHOA1CH3XC7H8O];
	w[rT43b] = k[rT43b] * c[sH] * c[sOA1CH3XC7H7O];
	w[rT44f] = k[rT44f] * c[sHOA1CH3XC7H8O] * c[sH];
	w[rT44b] = k[rT44b] * c[sH2] * c[sOA1CH3XC7H7O];
	w[rT45f] = k[rT45f] * c[sHOA1CH3XC7H8O] * c[sO];
	w[rT45b] = k[rT45b] * c[sOH] * c[sOA1CH3XC7H7O];
	w[rT46f] = k[rT46f] * c[sHOA1CH3XC7H8O] * c[sOH];
	w[rT46b] = k[rT46b] * c[sH2O] * c[sOA1CH3XC7H7O];
	w[rT47] = k[rT47] * c[sOA1CH3XC7H7O];
	w[rT50f] = k[rT50f] * c[sA1CH3XC7H8];
	w[rT50b] = k[rT50b] * c[sH] * c[sA1CH3YXC7H7];
	w[rT51f] = k[rT51f] * c[sA1CH3XC7H8] * c[sH];
	w[rT51b] = k[rT51b] * c[sH2] * c[sA1CH3YXC7H7];
	w[rT52f] = k[rT52f] * c[sA1CH3XC7H8] * c[sO];
	w[rT52b] = k[rT52b] * c[sOH] * c[sA1CH3YXC7H7];
	w[rT53f] = k[rT53f] * c[sA1CH3XC7H8] * c[sOH];
	w[rT53b] = k[rT53b] * c[sH2O] * c[sA1CH3YXC7H7];
	w[rT54f] = k[rT54f] * c[sA1CH3XC7H8] * c[sCH3];
	w[rT54b] = k[rT54b] * c[sCH4] * c[sA1CH3YXC7H7];
	w[rT55b] = k[rT55b] * c[sOA1CH3XC7H7O];
	w[rT56b] = k[rT56b] * c[sH] * c[sOA1CH3XC7H7O];
	w[rT57b] = k[rT57b] * c[sOH] * c[sOA1CH3XC7H7O];
	w[rT58f] = k[rT58f] * c[sA1CH3YXC7H7] * c[sO2];
	w[rT58b] = k[rT58b] * c[sO] * c[sOA1CH3XC7H7O];
	w[rT59] = k[rT59] * c[sA1CH3YXC7H7] * c[sO2];
	w[rT60] = k[rT60] * c[sA1CH3YXC7H7] * c[sO2];
	w[rE01b] = k[rE01b] * c[sA1C2H5XC8H10];
	w[rE02b] = k[rE02b] * c[sA1C2H5XC8H10];
	w[rE03b] = k[rE03b] * c[sA1C2H5XC8H10];
	w[rE04f] = k[rE04f] * c[sA1C2H5XC8H10] * c[sH];
	w[rE04b] = k[rE04b] * c[sC2H5] * c[sA1XC6H6];
	w[rE05f] = k[rE05f] * c[sA1C2H5XC8H10] * c[sOH];
	w[rE05b] = k[rE05b] * c[sC2H5] * c[sA1OHXC6H6O];
	w[rE06f] = k[rE06f] * c[sA1C2H5XC8H10] * c[sH];
	w[rE06b] = k[rE06b] * c[sH2] * c[sA1C2H4XC8H9];
	w[rE07f] = k[rE07f] * c[sA1C2H5XC8H10] * c[sO];
	w[rE07b] = k[rE07b] * c[sOH] * c[sA1C2H4XC8H9];
	w[rE08f] = k[rE08f] * c[sA1C2H5XC8H10] * c[sOH];
	w[rE08b] = k[rE08b] * c[sH2O] * c[sA1C2H4XC8H9];
	w[rE09f] = k[rE09f] * c[sA1C2H5XC8H10] * c[sHO2];
	w[rE09b] = k[rE09b] * c[sH2O2] * c[sA1C2H4XC8H9];
	w[rE10f] = k[rE10f] * c[sA1C2H5XC8H10] * c[sCH3];
	w[rE10b] = k[rE10b] * c[sCH4] * c[sA1C2H4XC8H9];
	w[rE11f] = k[rE11f] * c[sA1C2H4XC8H9];
	w[rE11b] = k[rE11b] * c[sC2H4] * c[sA1XXC6H5];
	w[rE12f] = k[rE12f] * c[sA1C2H4XC8H9];
	w[rE12b] = k[rE12b] * c[sH] * c[sA1C2H3XC8H8];
	w[rE13b] = k[rE13b] * c[sH2] * c[sA1C2H3XC8H8];
	w[rE14b] = k[rE14b] * c[sH2O] * c[sA1C2H3XC8H8];
	w[rE15f] = k[rE15f] * c[sA1C2H4XC8H9] * c[sCH3];
	w[rE15b] = k[rE15b] * c[sCH4] * c[sA1C2H3XC8H8];
	w[rE16f] = k[rE16f] * c[sA1C2H4XC8H9] * c[sO2];
	w[rE16b] = k[rE16b] * c[sHO2] * c[sA1C2H3XC8H8];
	w[rE17f] = k[rE17f] * c[sA1C2H4XC8H9] * c[sO];
	w[rE17b] = k[rE17b] * c[sCH2O] * c[sA1CH2XC7H7];
	w[rE18f] = k[rE18f] * c[sA1C2H4XC8H9] * c[sO];
	w[rE18b] = k[rE18b] * c[sCH3] * c[sA1CHOXC7H6O];
	w[rE19] = k[rE19] * c[sA1C2H4XC8H9] * c[sHO2];
	w[rE30] = k[rE30] * c[sA1C2H4XC8H9] * c[sO2];
	w[rE31] = k[rE31] * c[sA1C2H4XC8H9] * c[sO2];
	w[rE32] = k[rE32] * c[sC8H9O2];
	w[rE33] = k[rE33] * c[sC8H9O2];
	w[rE34] = k[rE34] * c[sC8H9O2];
	w[rE35] = k[rE35] * c[sC8H9O2];
	w[rE37] = k[rE37] * c[sOC8H7OOH];
	w[rST01] = k[rST01] * c[sA1C2H3XC8H8];
	w[rST02f] = k[rST02f] * c[sA1C2H3XC8H8] * c[sO];
	w[rST02b] = k[rST02b] * c[sHCO] * c[sA1CH2XC7H7];
	w[rST03f] = k[rST03f] * c[sA1C2H3XC8H8] * c[sCH3];
	w[rST03b] = k[rST03b] * c[sCH4] * c[sA1C2H2XC8H7];
	w[rST04f] = k[rST04f] * c[sA1C2H3XC8H8] * c[sOH];
	w[rST04b] = k[rST04b] * c[sCH3] * c[sA1CHOXC7H6O];
	w[rST05f] = k[rST05f] * c[sA1C2H3XC8H8] * c[sOH];
	w[rST05b] = k[rST05b] * c[sCH2O] * c[sA1CH2XC7H7];
	w[rST06f] = k[rST06f] * c[sA1C2H3XC8H8] * c[sOH];
	w[rST06b] = k[rST06b] * c[sA1OHXC6H6O] * c[sC2H3];
	w[rST10f] = k[rST10f] * c[sA1C2H3XC8H8] * c[sO];
	w[rST10b] = k[rST10b] * c[sOH] * c[sA1C2H3YXC8H7];
	w[rST11f] = k[rST11f] * c[sA1C2H2XC8H7] * c[sO];
	w[rST11b] = k[rST11b] * c[sCO] * c[sA1CH2XC7H7];
	w[rST12f] = k[rST12f] * c[sA1C2H2XC8H7] * c[sO2];
	w[rST12b] = k[rST12b] * c[sHO2] * c[sA1C2HXC8H6];
	w[rST13] = k[rST13] * c[sA1C2H2XC8H7] * c[sO2];
	w[rST14f] = k[rST14f] * c[sA1C2H2XC8H7] * c[sO2];
	w[rST14b] = k[rST14b] * c[sHCO] * c[sA1CHOXC7H6O];
	w[rST00b] = k[rST00b] * c[sHCCO] * c[sA1XC6H6];
	w[rXY00f] = k[rXY00f] * c[sA1CH3CH3XC8H10];
	w[rXY00b] = k[rXY00b] * c[sH] * c[sA1CH3CH2XC8H9];
	w[rXY01f] = k[rXY01f] * c[sA1CH3CH3XC8H10];
	w[rXY01b] = k[rXY01b] * c[sCH3] * c[sA1CH3YXC7H7];
	w[rXY02f] = k[rXY02f] * c[sA1CH3CH3XC8H10] * c[sH];
	w[rXY02b] = k[rXY02b] * c[sH2] * c[sA1CH3CH2XC8H9];
	w[rXY03f] = k[rXY03f] * c[sA1CH3CH3XC8H10] * c[sO];
	w[rXY03b] = k[rXY03b] * c[sOH] * c[sA1CH3CH2XC8H9];
	w[rXY04f] = k[rXY04f] * c[sA1CH3CH3XC8H10] * c[sOH];
	w[rXY04b] = k[rXY04b] * c[sH2O] * c[sA1CH3CH2XC8H9];
	w[rXY05f] = k[rXY05f] * c[sA1CH3CH3XC8H10] * c[sO2];
	w[rXY05b] = k[rXY05b] * c[sHO2] * c[sA1CH3CH2XC8H9];
	w[rXY06f] = k[rXY06f] * c[sA1CH3CH3XC8H10] * c[sHO2];
	w[rXY06b] = k[rXY06b] * c[sH2O2] * c[sA1CH3CH2XC8H9];
	w[rXY07f] = k[rXY07f] * c[sA1CH3CH3XC8H10] * c[sCH3];
	w[rXY07b] = k[rXY07b] * c[sCH4] * c[sA1CH3CH2XC8H9];
	w[rXY09f] = k[rXY09f] * c[sA1CH3CH3XC8H10] * c[sH];
	w[rXY09b] = k[rXY09b] * c[sCH3] * c[sA1CH3XC7H8];
	w[rXY10f] = k[rXY10f] * c[sA1CH3CH3XC8H10] * c[sOH];
	w[rXY10b] = k[rXY10b] * c[sCH3] * c[sHOA1CH3XC7H8O];
	w[rXY11] = k[rXY11] * c[sA1CH3CH3XC8H10] * c[sO];
	w[rXY12] = k[rXY12] * c[sA1CH3CH2XC8H9];
	w[rXY13f] = k[rXY13f] * c[sA1CH3CH2XC8H9] * c[sH];
	w[rXY13b] = k[rXY13b] * c[sCH3] * c[sA1CH3YXC7H7];
	w[rXY14f] = k[rXY14f] * c[sA1CH3CH2XC8H9] * c[sO];
	w[rXY14b] = k[rXY14b] * c[sH] * c[sA1CH3CHOXC8H8O];
	w[rXY15f] = k[rXY15f] * c[sA1CH3CH2XC8H9] * c[sO];
	w[rXY15b] = k[rXY15b] * c[sCH2O] * c[sA1CH3YXC7H7];
	w[rXY16f] = k[rXY16f] * c[sA1CH3CH2XC8H9] * c[sO];
	w[rXY16b] = k[rXY16b] * c[sHCO] * c[sA1CH3XC7H8];
	w[rXY18f] = k[rXY18f] * c[sA1CH3CH2XC8H9] * c[sO2];
	w[rXY18b] = k[rXY18b] * c[sOH] * c[sA1CH3CHOXC8H8O];
	w[rXY19f] = k[rXY19f] * c[sA1CH3CH2XC8H9] * c[sO2];
	w[rXY19b] = k[rXY19b] * c[sCH2O] * c[sOA1CH3XC7H7O];
	w[rXY201] = k[rXY201] * c[sA1CH3CH2XC8H9] * c[sHO2];
	w[rXY202] = k[rXY202] * c[sA1CH3CH2XC8H9] * c[sHO2];
	w[rXY203] = k[rXY203] * c[sA1CH3CH2XC8H9] * c[sHO2];
	w[rXY22] = k[rXY22] * c[sA1CH3CH2XC8H9] * c[sC3H3];
	w[rXY23f] = k[rXY23f] * c[sA1CH3CHOXC8H8O];
	w[rXY23b] = k[rXY23b] * c[sH] * c[sA1CHOCH2XC8H7O];
	w[rXY24] = k[rXY24] * c[sA1CH3CHOXC8H8O];
	w[rXY25] = k[rXY25] * c[sA1CH3CHOXC8H8O];
	w[rXY26f] = k[rXY26f] * c[sA1CH3CHOXC8H8O] * c[sH];
	w[rXY26b] = k[rXY26b] * c[sH2] * c[sA1CHOCH2XC8H7O];
	w[rXY27f] = k[rXY27f] * c[sA1CH3CHOXC8H8O] * c[sO];
	w[rXY27b] = k[rXY27b] * c[sOH] * c[sA1CHOCH2XC8H7O];
	w[rXY28f] = k[rXY28f] * c[sA1CH3CHOXC8H8O] * c[sOH];
	w[rXY28b] = k[rXY28b] * c[sH2O] * c[sA1CHOCH2XC8H7O];
	w[rXY29f] = k[rXY29f] * c[sA1CH3CHOXC8H8O] * c[sO2];
	w[rXY29b] = k[rXY29b] * c[sHO2] * c[sA1CHOCH2XC8H7O];
	w[rXY30f] = k[rXY30f] * c[sA1CH3CHOXC8H8O] * c[sHO2];
	w[rXY30b] = k[rXY30b] * c[sH2O2] * c[sA1CHOCH2XC8H7O];
	w[rXY31f] = k[rXY31f] * c[sA1CH3CHOXC8H8O] * c[sCH3];
	w[rXY31b] = k[rXY31b] * c[sCH4] * c[sA1CHOCH2XC8H7O];
	w[rXY33] = k[rXY33] * c[sA1CH3CHOXC8H8O] * c[sH];
	w[rXY34] = k[rXY34] * c[sA1CH3CHOXC8H8O] * c[sO];
	w[rXY35] = k[rXY35] * c[sA1CH3CHOXC8H8O] * c[sOH];
	w[rXY36] = k[rXY36] * c[sA1CH3CHOXC8H8O] * c[sO2];
	w[rXY37] = k[rXY37] * c[sA1CH3CHOXC8H8O] * c[sHO2];
	w[rXY38] = k[rXY38] * c[sA1CH3CHOXC8H8O] * c[sCH3];
	w[rXY39f] = k[rXY39f] * c[sA1CH3CHOXC8H8O] * c[sH];
	w[rXY39b] = k[rXY39b] * c[sHCO] * c[sA1CH3XC7H8];
	w[rXY40f] = k[rXY40f] * c[sA1CH3CHOXC8H8O] * c[sH];
	w[rXY40b] = k[rXY40b] * c[sCH3] * c[sA1CHOXC7H6O];
	w[rXY41f] = k[rXY41f] * c[sA1CH3CHOXC8H8O] * c[sOH];
	w[rXY41b] = k[rXY41b] * c[sHCO] * c[sHOA1CH3XC7H8O];
	w[rXY42] = k[rXY42] * c[sA1CH3CHOXC8H8O] * c[sOH];
	w[rXY43f] = k[rXY43f] * c[sA1CHOCH2XC8H7O] * c[sO];
	w[rXY43b] = k[rXY43b] * c[sH] * c[sA1CHOCHOXC8H6O2];
	w[rXY44] = k[rXY44] * c[sA1CHOCH2XC8H7O] * c[sO];
	w[rXY45f] = k[rXY45f] * c[sA1CHOCH2XC8H7O] * c[sO];
	w[rXY45b] = k[rXY45b] * c[sHCO] * c[sA1CHOXC7H6O];
	w[rXY47f] = k[rXY47f] * c[sA1CHOCH2XC8H7O] * c[sO2];
	w[rXY47b] = k[rXY47b] * c[sOH] * c[sA1CHOCHOXC8H6O2];
	w[rXY48] = k[rXY48] * c[sA1CHOCH2XC8H7O] * c[sHO2];
	w[rXY50] = k[rXY50] * c[sA1CHOCHOXC8H6O2];
	w[rXY51] = k[rXY51] * c[sA1CHOCHOXC8H6O2] * c[sH];
	w[rXY52] = k[rXY52] * c[sA1CHOCHOXC8H6O2] * c[sO];
	w[rXY53] = k[rXY53] * c[sA1CHOCHOXC8H6O2] * c[sOH];
	w[rXY54] = k[rXY54] * c[sA1CHOCHOXC8H6O2] * c[sO2];
	w[rXY55] = k[rXY55] * c[sA1CHOCHOXC8H6O2] * c[sHO2];
	w[rXY56] = k[rXY56] * c[sA1CHOCHOXC8H6O2] * c[sCH3];
	w[rXY57f] = k[rXY57f] * c[sA1CHOCHOXC8H6O2] * c[sH];
	w[rXY57b] = k[rXY57b] * c[sHCO] * c[sA1CHOXC7H6O];
	w[rXY58] = k[rXY58] * c[sA1CHOCHOXC8H6O2] * c[sOH];
	w[rN01f] = k[rN01f] * c[sA2CH3XC11H10] * c[sH];
	w[rN01b] = k[rN01b] * c[sCH3] * c[sA2XC10H8];
	w[rN02f] = k[rN02f] * c[sA2CH3XC11H10] * c[sOH];
	w[rN02b] = k[rN02b] * c[sCH3] * c[sA2OHXC10H8O];
	w[rN04f] = k[rN04f] * c[sA2CH3XC11H10];
	w[rN04b] = k[rN04b] * c[sH] * c[sA2CH2XC11H9];
	w[rN05f] = k[rN05f] * c[sA2CH3XC11H10];
	w[rN05b] = k[rN05b] * c[sCH3] * c[sA2XXC10H7];
	w[rN06f] = k[rN06f] * c[sA2CH2XC11H9] * c[sH];
	w[rN06b] = k[rN06b] * c[sCH3] * c[sA2XXC10H7];
	w[rN07f] = k[rN07f] * c[sA2CH2XC11H9];
	w[rN07b] = k[rN07b] * c[sC2H2] * c[sC9H7];
	w[rN08f] = k[rN08f] * c[sA2CH3XC11H10] * c[sH];
	w[rN08b] = k[rN08b] * c[sH2] * c[sA2CH2XC11H9];
	w[rN09f] = k[rN09f] * c[sA2CH3XC11H10] * c[sO];
	w[rN09b] = k[rN09b] * c[sOH] * c[sA2CH2XC11H9];
	w[rN10f] = k[rN10f] * c[sA2CH3XC11H10] * c[sOH];
	w[rN10b] = k[rN10b] * c[sH2O] * c[sA2CH2XC11H9];
	w[rN11f] = k[rN11f] * c[sA2CH3XC11H10] * c[sO2];
	w[rN11b] = k[rN11b] * c[sHO2] * c[sA2CH2XC11H9];
	w[rN12f] = k[rN12f] * c[sA2CH3XC11H10] * c[sCH3];
	w[rN12b] = k[rN12b] * c[sCH4] * c[sA2CH2XC11H9];
	w[rN13f] = k[rN13f] * c[sA2CH3XC11H10] * c[sHO2];
	w[rN13b] = k[rN13b] * c[sH2O2] * c[sA2CH2XC11H9];
	w[rN14] = k[rN14] * c[sA2CH3XC11H10] * c[sO];
	w[rN15] = k[rN15] * c[sA2CH3XC11H10] * c[sO];
	w[rN16b] = k[rN16b] * c[sA2CH2OXC11H9O];
	w[rN18f] = k[rN18f] * c[sA2CH2XC11H9] * c[sHO2];
	w[rN18b] = k[rN18b] * c[sOH] * c[sA2CH2OXC11H9O];
	w[rN19] = k[rN19] * c[sA2CH2XC11H9] * c[sC3H3];
	w[rN20f] = k[rN20f] * c[sA2CH2XC11H9] * c[sO2];
	w[rN20b] = k[rN20b] * c[sOH] * c[sA2CHOXC11H8O];
	w[rN21f] = k[rN21f] * c[sA2CH2XC11H9] * c[sO2];
	w[rN21b] = k[rN21b] * c[sCH2O] * c[sA2OXC10H7O];
	w[rN22f] = k[rN22f] * c[sA2CH2OXC11H9O];
	w[rN22b] = k[rN22b] * c[sH] * c[sA2CHOXC11H8O];
	w[rN23f] = k[rN23f] * c[sA2CH2OXC11H9O];
	w[rN23b] = k[rN23b] * c[sCH2O] * c[sA2XXC10H7];
	w[rN24b] = k[rN24b] * c[sH2] * c[sA2CHOXC11H8O];
	w[rN25b] = k[rN25b] * c[sOH] * c[sA2CHOXC11H8O];
	w[rN26b] = k[rN26b] * c[sH2O] * c[sA2CHOXC11H8O];
	w[rN27f] = k[rN27f] * c[sA2CH2OXC11H9O] * c[sO2];
	w[rN27b] = k[rN27b] * c[sHO2] * c[sA2CHOXC11H8O];
	w[rN28] = k[rN28] * c[sA2CHOXC11H8O];
	w[rN29] = k[rN29] * c[sA2CHOXC11H8O] * c[sH];
	w[rN30] = k[rN30] * c[sA2CHOXC11H8O] * c[sO];
	w[rN31] = k[rN31] * c[sA2CHOXC11H8O] * c[sOH];
	w[rN32] = k[rN32] * c[sA2CHOXC11H8O] * c[sO2];
	w[rN33] = k[rN33] * c[sA2CHOXC11H8O] * c[sHO2];
	w[rN34] = k[rN34] * c[sA2CHOXC11H8O] * c[sCH3];
	w[rOX00f] = k[rOX00f] * c[sA1XC6H6];
	w[rOX00b] = k[rOX00b] * c[sH] * c[sA1XXC6H5];
	w[rOX02f] = k[rOX02f] * c[sC4H2] * c[sC2H2];
	w[rOX02b] = k[rOX02b] * c[sOXC6H4];
	w[rOX03f] = k[rOX03f] * c[sA1XC6H6] * c[sH];
	w[rOX03b] = k[rOX03b] * c[sH2] * c[sA1XXC6H5];
	w[rOX04f] = k[rOX04f] * c[sA1XC6H6] * c[sOH];
	w[rOX04b] = k[rOX04b] * c[sH2O] * c[sA1XXC6H5];
	w[rOX05f] = k[rOX05f] * c[sA1XC6H6] * c[sOH];
	w[rOX05b] = k[rOX05b] * c[sH] * c[sA1OHXC6H6O];
	w[rOX06f] = k[rOX06f] * c[sA1XC6H6] * c[sO2];
	w[rOX06b] = k[rOX06b] * c[sHO2] * c[sA1XXC6H5];
	w[rOX99f] = k[rOX99f] * c[sA1XC6H6] * c[sO];
	w[rOX99b] = k[rOX99b] * c[sH] * c[sA1OXC6H5O];
	w[rOX07f] = k[rOX07f] * c[sA1XC6H6] * c[sO];
	w[rOX07b] = k[rOX07b] * c[sH] * c[sA1OXC6H5O];
	w[rOX08f] = k[rOX08f] * c[sA1XC6H6] * c[sO];
	w[rOX08b] = k[rOX08b] * c[sA1OHXC6H6O];
	w[rOX09f] = k[rOX09f] * c[sA1XC6H6] * c[sO];
	w[rOX09b] = k[rOX09b] * c[sCO] * c[sC5H6];
	w[rOX10f] = k[rOX10f] * c[sA1XC6H6] * c[sO];
	w[rOX10b] = k[rOX10b] * c[sOH] * c[sA1XXC6H5];
	w[rOX11f] = k[rOX11f] * c[sA1XXC6H5] * c[sO2];
	w[rOX11b] = k[rOX11b] * c[sO] * c[sA1OXC6H5O];
	w[rOX12f] = k[rOX12f] * c[sA1XXC6H5] * c[sO2];
	w[rOX12b] = k[rOX12b] * c[sH] * c[sOC6H4O];
	w[rOX13b] = k[rOX13b] * c[sA1OXC6H5O];
	w[rOX14b] = k[rOX14b] * c[sH] * c[sA1OXC6H5O];
	w[rOX15b] = k[rOX15b] * c[sOH] * c[sA1OXC6H5O];
	w[rOX16f] = k[rOX16f] * c[sA1XXC6H5] * c[sCH4];
	w[rOX16b] = k[rOX16b] * c[sCH3] * c[sA1XC6H6];
	w[rOX17f] = k[rOX17f] * c[sA1OHXC6H6O];
	w[rOX17b] = k[rOX17b] * c[sCO] * c[sC5H6];
	w[rOX18f] = k[rOX18f] * c[sA1OHXC6H6O];
	w[rOX18b] = k[rOX18b] * c[sH] * c[sA1OXC6H5O];
	w[rOX19f] = k[rOX19f] * c[sA1OHXC6H6O] * c[sH];
	w[rOX19b] = k[rOX19b] * c[sH2] * c[sA1OXC6H5O];
	w[rOX20f] = k[rOX20f] * c[sA1OHXC6H6O] * c[sOH];
	w[rOX20b] = k[rOX20b] * c[sH2O] * c[sA1OXC6H5O];
	w[rOX21f] = k[rOX21f] * c[sA1OHXC6H6O] * c[sCH3];
	w[rOX21b] = k[rOX21b] * c[sCH4] * c[sA1OXC6H5O];
	w[rOX22f] = k[rOX22f] * c[sA1OHXC6H6O] * c[sO2];
	w[rOX22b] = k[rOX22b] * c[sHO2] * c[sA1OXC6H5O];
	w[rOX23f] = k[rOX23f] * c[sA1OXC6H5O];
	w[rOX23b] = k[rOX23b] * c[sC5H5] * c[sCO];
	w[rOX24f] = k[rOX24f] * c[sA1OXC6H5O] * c[sO];
	w[rOX24b] = k[rOX24b] * c[sH] * c[sOC6H4O];
	w[rOX25f] = k[rOX25f] * c[sA1OXC6H5O] * c[sO2];
	w[rOX25b] = k[rOX25b] * c[sOH] * c[sOC6H4O];
	w[rOX26f] = k[rOX26f] * c[sOC6H4O];
	w[rOX26b] = k[rOX26b] * c[sCO] * c[sC5H4O];
	w[rOX27f] = k[rOX27f] * c[sOC6H4O] * c[sH];
	w[rOX27b] = k[rOX27b] * c[sCO] * c[sTXC5H5O];
	w[rOX28] = k[rOX28] * c[sOC6H4O] * c[sO];
	w[rOX30f] = k[rOX30f] * c[sA2XC10H8] * c[sO];
	w[rOX30b] = k[rOX30b] * c[sH] * c[sA2OXC10H7O];
	w[rOX31f] = k[rOX31f] * c[sA2XC10H8] * c[sO];
	w[rOX31b] = k[rOX31b] * c[sA2OHXC10H8O];
	w[rOX32f] = k[rOX32f] * c[sA2XC10H8] * c[sOH];
	w[rOX32b] = k[rOX32b] * c[sH] * c[sA2OHXC10H8O];
	w[rOX33f] = k[rOX33f] * c[sA2XXC10H7] * c[sO2];
	w[rOX33b] = k[rOX33b] * c[sO] * c[sA2OXC10H7O];
	w[rOX34f] = k[rOX34f] * c[sA2YXC10H7] * c[sO2];
	w[rOX34b] = k[rOX34b] * c[sO] * c[sA2OXC10H7O];
	w[rOX35] = k[rOX35] * c[sA2XXC10H7] * c[sO2];
	w[rOX36] = k[rOX36] * c[sA2YXC10H7] * c[sO2];
	w[rOX37b] = k[rOX37b] * c[sA2OXC10H7O];
	w[rOX38b] = k[rOX38b] * c[sA2OXC10H7O];
	w[rOX39b] = k[rOX39b] * c[sH] * c[sA2OXC10H7O];
	w[rOX40b] = k[rOX40b] * c[sH] * c[sA2OXC10H7O];
	w[rOX41f] = k[rOX41f] * c[sA2OHXC10H8O];
	w[rOX41b] = k[rOX41b] * c[sCO] * c[sC9H8];
	w[rOX42f] = k[rOX42f] * c[sA2OHXC10H8O];
	w[rOX42b] = k[rOX42b] * c[sH] * c[sA2OXC10H7O];
	w[rOX43f] = k[rOX43f] * c[sA2OHXC10H8O] * c[sH];
	w[rOX43b] = k[rOX43b] * c[sH2] * c[sA2OXC10H7O];
	w[rOX44f] = k[rOX44f] * c[sA2OHXC10H8O] * c[sOH];
	w[rOX44b] = k[rOX44b] * c[sH2O] * c[sA2OXC10H7O];
	w[rOX45f] = k[rOX45f] * c[sA2OHXC10H8O] * c[sCH3];
	w[rOX45b] = k[rOX45b] * c[sCH4] * c[sA2OXC10H7O];
	w[rOX46f] = k[rOX46f] * c[sA2OXC10H7O];
	w[rOX46b] = k[rOX46b] * c[sCO] * c[sC9H7];
	w[rOX48] = k[rOX48] * c[sA2OXC10H7O] * c[sO2];
	w[rOX50] = k[rOX50] * c[sA3XXC14H9] * c[sO2];
	w[rOX51] = k[rOX51] * c[sA3YXC14H9] * c[sO2];
	w[rOX52] = k[rOX52] * c[sA4XXC16H9] * c[sO2];
	w[rOX53] = k[rOX53] * c[sA2R5XXC12H7] * c[sO2];
	w[rOX54] = k[rOX54] * c[sA3R5XXC16H9] * c[sO2];
	w[rOX60] = k[rOX60] * c[sA3XC14H10] * c[sOH];
	w[rOX61] = k[rOX61] * c[sA3XC14H10] * c[sOH];
	w[rOX62] = k[rOX62] * c[sA4XC16H10] * c[sOH];
	w[rOX63] = k[rOX63] * c[sA2R5XC12H8] * c[sOH];
	w[rOX64] = k[rOX64] * c[sA3R5XC16H10] * c[sOH];
	w[rOX65] = k[rOX65] * c[sA4R5XC18H10] * c[sOH];
	w[rOX66] = k[rOX66] * c[sFLTNXC16H10] * c[sOH];
	w[rNO178f] = k[rNO178f] * c[sN] * c[sNO];
	w[rNO178b] = k[rNO178b] * c[sO] * c[sN2];
	w[rNO179f] = k[rNO179f] * c[sN] * c[sO2];
	w[rNO179b] = k[rNO179b] * c[sO] * c[sNO];
	w[rNO180f] = k[rNO180f] * c[sN] * c[sOH];
	w[rNO180b] = k[rNO180b] * c[sH] * c[sNO];
	w[rNO181f] = k[rNO181f] * c[sN2O] * c[sO];
	w[rNO181b] = k[rNO181b] * c[sO2] * c[sN2];
	w[rNO182f] = k[rNO182f] * c[sN2O] * c[sO];
	w[rNO182b] = k[rNO182b] * c[sNO] * c[sNO];
	w[rNO183f] = k[rNO183f] * c[sN2O] * c[sH];
	w[rNO183b] = k[rNO183b] * c[sOH] * c[sN2];
	w[rNO184f] = k[rNO184f] * c[sN2O] * c[sOH];
	w[rNO184b] = k[rNO184b] * c[sHO2] * c[sN2];
	w[rNO186f] = k[rNO186f] * c[sHO2] * c[sNO];
	w[rNO186b] = k[rNO186b] * c[sOH] * c[sNO2];
	w[rNO187f] = k[rNO187f] * c[sNO] * c[sO] * M[mM32];
	w[rNO187b] = k[rNO187b] * c[sNO2] * M[mM32];
	w[rNO188f] = k[rNO188f] * c[sNO2] * c[sO];
	w[rNO188b] = k[rNO188b] * c[sO2] * c[sNO];
	w[rNO189f] = k[rNO189f] * c[sNO2] * c[sH];
	w[rNO189b] = k[rNO189b] * c[sOH] * c[sNO];
	w[rNO190b] = k[rNO190b] * c[sH] * c[sNO];
	w[rNO191f] = k[rNO191f] * c[sNH] * c[sH];
	w[rNO191b] = k[rNO191b] * c[sH2] * c[sN];
	w[rNO192b] = k[rNO192b] * c[sH] * c[sHNO];
	w[rNO193f] = k[rNO193f] * c[sNH] * c[sOH];
	w[rNO193b] = k[rNO193b] * c[sH2O] * c[sN];
	w[rNO194f] = k[rNO194f] * c[sNH] * c[sO2];
	w[rNO194b] = k[rNO194b] * c[sO] * c[sHNO];
	w[rNO195f] = k[rNO195f] * c[sNH] * c[sO2];
	w[rNO195b] = k[rNO195b] * c[sOH] * c[sNO];
	w[rNO196b] = k[rNO196b] * c[sH] * c[sN2];
	w[rNO197f] = k[rNO197f] * c[sNH] * c[sH2O];
	w[rNO197b] = k[rNO197b] * c[sH2] * c[sHNO];
	w[rNO198f] = k[rNO198f] * c[sNH] * c[sNO];
	w[rNO198b] = k[rNO198b] * c[sOH] * c[sN2];
	w[rNO199f] = k[rNO199f] * c[sNH] * c[sNO];
	w[rNO199b] = k[rNO199b] * c[sH] * c[sN2O];
	w[rNO200b] = k[rNO200b] * c[sNH] * c[sOH];
	w[rNO201b] = k[rNO201b] * c[sHNO] * c[sH];
	w[rNO202f] = k[rNO202f] * c[sNH2] * c[sH];
	w[rNO202b] = k[rNO202b] * c[sH2] * c[sNH];
	w[rNO203f] = k[rNO203f] * c[sNH2] * c[sOH];
	w[rNO203b] = k[rNO203b] * c[sH2O] * c[sNH];
	w[rNO204b] = k[rNO204b] * c[sH] * c[sN2];
	w[rNO205f] = k[rNO205f] * c[sNNH] * M[mM32];
	w[rNO205b] = k[rNO205b] * c[sH] * c[sN2] * M[mM32];
	w[rNO206b] = k[rNO206b] * c[sN2] * c[sHO2];
	w[rNO207b] = k[rNO207b] * c[sN2] * c[sOH];
	w[rNO208b] = k[rNO208b] * c[sNO] * c[sNH];
	w[rNO209b] = k[rNO209b] * c[sN2] * c[sH2];
	w[rNO210b] = k[rNO210b] * c[sN2] * c[sH2O];
	w[rNO211b] = k[rNO211b] * c[sN2] * c[sCH4];
	w[rNO212f] = k[rNO212f] * c[sH] * c[sNO] * M[mM32];
	w[rNO212b] = k[rNO212b] * c[sHNO] * M[mM32];
	w[rNO213b] = k[rNO213b] * c[sOH] * c[sNO];
	w[rNO214f] = k[rNO214f] * c[sHNO] * c[sH];
	w[rNO214b] = k[rNO214b] * c[sNO] * c[sH2];
	w[rNO215f] = k[rNO215f] * c[sHNO] * c[sOH];
	w[rNO215b] = k[rNO215b] * c[sH2O] * c[sNO];
	w[rNO216f] = k[rNO216f] * c[sHNO] * c[sO2];
	w[rNO216b] = k[rNO216b] * c[sNO] * c[sHO2];
	w[rNO217b] = k[rNO217b] * c[sN] * c[sCO];
	w[rNO218b] = k[rNO218b] * c[sH] * c[sNCO];
	w[rNO219f] = k[rNO219f] * c[sCN] * c[sH2O];
	w[rNO219b] = k[rNO219b] * c[sOH] * c[sHCN];
	w[rNO220f] = k[rNO220f] * c[sCN] * c[sO2];
	w[rNO220b] = k[rNO220b] * c[sO] * c[sNCO];
	w[rNO221f] = k[rNO221f] * c[sCN] * c[sH2];
	w[rNO221b] = k[rNO221b] * c[sH] * c[sHCN];
	w[rNO222b] = k[rNO222b] * c[sCO] * c[sNO];
	w[rNO223b] = k[rNO223b] * c[sCO] * c[sNH];
	w[rNO224b] = k[rNO224b] * c[sCO] * c[sH] * c[sNO];
	w[rNO225b] = k[rNO225b] * c[sCO] * c[sN2];
	w[rNO226f] = k[rNO226f] * c[sNCO] * c[sO2];
	w[rNO226b] = k[rNO226b] * c[sCO2] * c[sNO];
	w[rNO227f] = k[rNO227f] * c[sNCO] * M[mM32];
	w[rNO227b] = k[rNO227b] * c[sCO] * c[sN] * M[mM32];
	w[rNO228f] = k[rNO228f] * c[sNCO] * c[sNO];
	w[rNO228b] = k[rNO228b] * c[sCO] * c[sN2O];
	w[rNO229f] = k[rNO229f] * c[sNCO] * c[sNO];
	w[rNO229b] = k[rNO229b] * c[sCO2] * c[sN2];
	w[rNO230f] = k[rNO230f] * c[sHCN] * M[mM32];
	w[rNO230b] = k[rNO230b] * c[sCN] * c[sH] * M[mM32];
	w[rNO231f] = k[rNO231f] * c[sHCN] * c[sO];
	w[rNO231b] = k[rNO231b] * c[sH] * c[sNCO];
	w[rNO232f] = k[rNO232f] * c[sHCN] * c[sO];
	w[rNO232b] = k[rNO232b] * c[sCO] * c[sNH];
	w[rNO233f] = k[rNO233f] * c[sHCN] * c[sO];
	w[rNO233b] = k[rNO233b] * c[sOH] * c[sCN];
	w[rNO234f] = k[rNO234f] * c[sHCN] * c[sOH];
	w[rNO234b] = k[rNO234b] * c[sH] * c[sHOCN];
	w[rNO235] = k[rNO235] * c[sHCN] * c[sOH];
	w[rNO236f] = k[rNO236f] * c[sHCN] * c[sOH];
	w[rNO236b] = k[rNO236b] * c[sCO] * c[sNH2];
	w[rNO238f] = k[rNO238f] * c[sH2CN] * c[sN];
	w[rNO238b] = k[rNO238b] * c[sTXCH2] * c[sN2];
	w[rNO239f] = k[rNO239f] * c[sC] * c[sN2];
	w[rNO239b] = k[rNO239b] * c[sN] * c[sCN];
	w[rNO240f] = k[rNO240f] * c[sCH] * c[sN2];
	w[rNO240b] = k[rNO240b] * c[sN] * c[sHCN];
	w[rNO242f] = k[rNO242f] * c[sTXCH2] * c[sN2];
	w[rNO242b] = k[rNO242b] * c[sNH] * c[sHCN];
	w[rNO243f] = k[rNO243f] * c[sSXCH2] * c[sN2];
	w[rNO243b] = k[rNO243b] * c[sHCN] * c[sNH];
	w[rNO244b] = k[rNO244b] * c[sO] * c[sCN];
	w[rNO245b] = k[rNO245b] * c[sN] * c[sCO];
	w[rNO246b] = k[rNO246b] * c[sO] * c[sHCN];
	w[rNO247b] = k[rNO247b] * c[sNCO] * c[sH];
	w[rNO248b] = k[rNO248b] * c[sHCO] * c[sN];
	w[rNO249f] = k[rNO249f] * c[sTXCH2] * c[sNO];
	w[rNO249b] = k[rNO249b] * c[sHNCO] * c[sH];
	w[rNO250f] = k[rNO250f] * c[sTXCH2] * c[sNO];
	w[rNO250b] = k[rNO250b] * c[sHCN] * c[sOH];
	w[rNO251f] = k[rNO251f] * c[sTXCH2] * c[sNO];
	w[rNO251b] = k[rNO251b] * c[sHCNO] * c[sH];
	w[rNO252f] = k[rNO252f] * c[sSXCH2] * c[sNO];
	w[rNO252b] = k[rNO252b] * c[sHNCO] * c[sH];
	w[rNO253f] = k[rNO253f] * c[sSXCH2] * c[sNO];
	w[rNO253b] = k[rNO253b] * c[sHCN] * c[sOH];
	w[rNO254f] = k[rNO254f] * c[sSXCH2] * c[sNO];
	w[rNO254b] = k[rNO254b] * c[sHCNO] * c[sH];
	w[rNO255f] = k[rNO255f] * c[sCH3] * c[sNO];
	w[rNO255b] = k[rNO255b] * c[sH2O] * c[sHCN];
	w[rNO256f] = k[rNO256f] * c[sCH3] * c[sNO];
	w[rNO256b] = k[rNO256b] * c[sOH] * c[sH2CN];
	w[rNO257b] = k[rNO257b] * c[sN2] * c[sH] * c[sCO];
	w[rNO258b] = k[rNO258b] * c[sNO] * c[sHCN];
	w[rNO259b] = k[rNO259b] * c[sN2] * c[sHCO] * c[sO];
	w[rNO260b] = k[rNO260b] * c[sN2] * c[sHCO] * c[sH];
	w[rNO261b] = k[rNO261b] * c[sN2] * c[sTXCH2];
	w[rNO262f] = k[rNO262f] * c[sHNCO] * c[sO];
	w[rNO262b] = k[rNO262b] * c[sCO2] * c[sNH];
	w[rNO263f] = k[rNO263f] * c[sHNCO] * c[sO];
	w[rNO263b] = k[rNO263b] * c[sCO] * c[sHNO];
	w[rNO264f] = k[rNO264f] * c[sHNCO] * c[sO];
	w[rNO264b] = k[rNO264b] * c[sOH] * c[sNCO];
	w[rNO265f] = k[rNO265f] * c[sHNCO] * c[sH];
	w[rNO265b] = k[rNO265b] * c[sCO] * c[sNH2];
	w[rNO266f] = k[rNO266f] * c[sHNCO] * c[sH];
	w[rNO266b] = k[rNO266b] * c[sNCO] * c[sH2];
	w[rNO267f] = k[rNO267f] * c[sHNCO] * c[sOH];
	w[rNO267b] = k[rNO267b] * c[sH2O] * c[sNCO];
	w[rNO268f] = k[rNO268f] * c[sHNCO] * c[sOH];
	w[rNO268b] = k[rNO268b] * c[sCO2] * c[sNH2];
	w[rNO269f] = k[rNO269f] * c[sHNCO] * M[mM32];
	w[rNO269b] = k[rNO269b] * c[sCO] * c[sNH] * M[mM32];
	w[rNO270f] = k[rNO270f] * c[sHCNO] * c[sH];
	w[rNO270b] = k[rNO270b] * c[sHNCO] * c[sH];
	w[rNO271f] = k[rNO271f] * c[sHCNO] * c[sH];
	w[rNO271b] = k[rNO271b] * c[sHCN] * c[sOH];
	w[rNO272f] = k[rNO272f] * c[sHCNO] * c[sH];
	w[rNO272b] = k[rNO272b] * c[sCO] * c[sNH2];
	w[rNO273] = k[rNO273] * c[sHOCN] * c[sH];
	w[rNO274b] = k[rNO274b] * c[sCO] * c[sHCNO];
	w[rNO275f] = k[rNO275f] * c[sCH3] * c[sN];
	w[rNO275b] = k[rNO275b] * c[sH] * c[sH2CN];
	w[rNO276f] = k[rNO276f] * c[sCH3] * c[sN];
	w[rNO276b] = k[rNO276b] * c[sH2] * c[sHCN];
	w[rNO277f] = k[rNO277f] * c[sNH3] * c[sH];
	w[rNO277b] = k[rNO277b] * c[sH2] * c[sNH2];
	w[rNO278f] = k[rNO278f] * c[sNH3] * c[sOH];
	w[rNO278b] = k[rNO278b] * c[sH2O] * c[sNH2];
	w[rNO279f] = k[rNO279f] * c[sNH3] * c[sO];
	w[rNO279b] = k[rNO279b] * c[sOH] * c[sNH2];
	w[rNO280f] = k[rNO280f] * c[sNH] * c[sCO2];
	w[rNO280b] = k[rNO280b] * c[sCO] * c[sHNO];
	w[rNO281f] = k[rNO281f] * c[sCN] * c[sNO2];
	w[rNO281b] = k[rNO281b] * c[sNO] * c[sNCO];
	w[rNO282f] = k[rNO282f] * c[sNCO] * c[sNO2];
	w[rNO282b] = k[rNO282b] * c[sCO2] * c[sN2O];
	w[rNO283f] = k[rNO283f] * c[sN] * c[sCO2];
	w[rNO283b] = k[rNO283b] * c[sCO] * c[sNO];
	w[rG09] = k[rG09] * c[sCH] * c[sO2];
	w[rG15] = k[rG15] * c[sTXCH2] * c[sO];
	w[rG16] = k[rG16] * c[sTXCH2] * c[sOH];
	w[rG36] = k[rG36] * c[sSXCH2] * c[sO2];
	w[rG35] = k[rG35] * c[sSXCH2] * c[sO2];
	w[rG38f] = k[rG38f] * c[sSXCH2] * c[sH2O];
	w[rG34f] = k[rG34f] * c[sSXCH2] * c[sH2];
	w[rG52] = k[rG52] * c[sCH3] * c[sO];
	w[rG53] = k[rG53] * c[sCH3] * c[sO];
	w[rG72] = k[rG72] * c[sCH3] * c[sTXCH2];
	w[rLPG65] = k[rLPG65] * c[sCH3] * c[sHO2];
	w[r32] = k[r32] * c[sHCO] * c[sH];
	w[r35] = k[r35] * c[sHCO] * c[sOH];
	w[rDD54] = k[rDD54] * c[sHCO] * c[sOH];
	w[rG130] = k[rG130] * c[sC2H3] * c[sH];
	w[rR037] = k[rR037] * c[sC2H5] * c[sHO2];
	w[rG109] = k[rG109] * c[sHCCO] * c[sH];
	w[rR022] = k[rR022] * c[sHCCO] * c[sCH3];
	w[rG110] = k[rG110] * c[sHCCO] * c[sO];
	w[rR024] = k[rR024] * c[sHCCO] * c[sOH];
	w[rP013f] = k[rP013f] * c[sA1XXC6H5] * c[sC2H3];
	w[rP022f] = k[rP022f] * c[sA1C2H2XC8H7] * c[sOH];
	w[rP116f] = k[rP116f] * c[sA2XXC10H7] * c[sC2H3];
	w[rP117f] = k[rP117f] * c[sA2YXC10H7] * c[sC2H3];
	w[rP124f] = k[rP124f] * c[sA2C2H2AXC12H9] * c[sOH];
	w[rP127f] = k[rP127f] * c[sA2C2H2BXC12H9] * c[sOH];
	w[rP206f] = k[rP206f] * c[sA2R5XXC12H7] * c[sC2H3];
	w[rP214f] = k[rP214f] * c[sA2R5C2H2XC14H9] * c[sOH];
	w[rP408f] = k[rP408f] * c[sA2C2HAYXC12H7] * c[sC2H3];
	w[rP409f] = k[rP409f] * c[sA2C2HBYXC12H7] * c[sC2H3];
	w[rP504f] = k[rP504f] * c[sA3YXC14H9] * c[sC2H3];
	w[rCP10f] = k[rCP10f] * c[sC5H6] * c[sC2H3];
	w[rCP11f] = k[rCP11f] * c[sC5H6] * c[sNXC4H5];
	w[rCP19f] = k[rCP19f] * c[sC5H5] * c[sO];
	w[rCP21f] = k[rCP21f] * c[sC5H5] * c[sHO2];
	w[rCP22] = k[rCP22] * c[sC5H5] * c[sOH];
	w[rCP25] = k[rCP25] * c[sSXC5H5O] * c[sH];
	w[rCP26] = k[rCP26] * c[sSXC5H5O] * c[sO];
	w[rCP27] = k[rCP27] * c[sSXC5H5O] * c[sOH];
	w[rCP29] = k[rCP29] * c[sTXC5H5O] * c[sH];
	w[rI00f] = k[rI00f] * c[sA1XXC6H5] * c[sC3H3];
	w[rI19f] = k[rI19f] * c[sC9H7] * c[sO];
	w[rI21] = k[rI21] * c[sC9H7] * c[sHO2];
	w[rI22] = k[rI22] * c[sC9H7] * c[sOH];
	w[rT17f] = k[rT17f] * c[sA1CH2XC7H7] * c[sO];
	w[rT18f] = k[rT18f] * c[sA1CH2XC7H7] * c[sOH];
	w[rT32f] = k[rT32f] * c[sA1CH2OXC7H7O] * c[sH];
	w[rT33f] = k[rT33f] * c[sA1CH2OXC7H7O] * c[sO];
	w[rT34f] = k[rT34f] * c[sA1CH2OXC7H7O] * c[sOH];
	w[rT55f] = k[rT55f] * c[sA1CH3YXC7H7] * c[sO];
	w[rT56f] = k[rT56f] * c[sA1CH3YXC7H7] * c[sOH];
	w[rT57f] = k[rT57f] * c[sA1CH3YXC7H7] * c[sHO2];
	w[rE01f] = k[rE01f] * c[sA1C2H4XC8H9] * c[sH];
	w[rE02f] = k[rE02f] * c[sA1CH2XC7H7] * c[sCH3];
	w[rE03f] = k[rE03f] * c[sA1XXC6H5] * c[sC2H5];
	w[rE13f] = k[rE13f] * c[sA1C2H4XC8H9] * c[sH];
	w[rE14f] = k[rE14f] * c[sA1C2H4XC8H9] * c[sOH];
	w[rE36] = k[rE36] * c[sC8H8OOH] * c[sO2];
	w[rST00f] = k[rST00f] * c[sA1C2HXC8H6] * c[sOH];
	w[rXY17] = k[rXY17] * c[sA1CH3CH2XC8H9] * c[sOH];
	w[rXY46] = k[rXY46] * c[sA1CHOCH2XC8H7O] * c[sOH];
	w[rN16f] = k[rN16f] * c[sA2CH2XC11H9] * c[sO];
	w[rN17] = k[rN17] * c[sA2CH2XC11H9] * c[sOH];
	w[rN24f] = k[rN24f] * c[sA2CH2OXC11H9O] * c[sH];
	w[rN25f] = k[rN25f] * c[sA2CH2OXC11H9O] * c[sO];
	w[rN26f] = k[rN26f] * c[sA2CH2OXC11H9O] * c[sOH];
	w[rOX13f] = k[rOX13f] * c[sA1XXC6H5] * c[sO];
	w[rOX14f] = k[rOX14f] * c[sA1XXC6H5] * c[sOH];
	w[rOX15f] = k[rOX15f] * c[sA1XXC6H5] * c[sHO2];
	w[rOX37f] = k[rOX37f] * c[sA2XXC10H7] * c[sO];
	w[rOX38f] = k[rOX38f] * c[sA2YXC10H7] * c[sO];
	w[rOX39f] = k[rOX39f] * c[sA2XXC10H7] * c[sOH];
	w[rOX40f] = k[rOX40f] * c[sA2YXC10H7] * c[sOH];
	w[rOX47] = k[rOX47] * c[sA2OXC10H7O] * c[sO];
	w[rNO190f] = k[rNO190f] * c[sNH] * c[sO];
	w[rNO192f] = k[rNO192f] * c[sNH] * c[sOH];
	w[rNO196f] = k[rNO196f] * c[sNH] * c[sN];
	w[rNO200f] = k[rNO200f] * c[sNH2] * c[sO];
	w[rNO201f] = k[rNO201f] * c[sNH2] * c[sO];
	w[rNO204f] = k[rNO204f] * c[sNNH];
	w[rNO206f] = k[rNO206f] * c[sNNH] * c[sO2];
	w[rNO207f] = k[rNO207f] * c[sNNH] * c[sO];
	w[rNO208f] = k[rNO208f] * c[sNNH] * c[sO];
	w[rNO209f] = k[rNO209f] * c[sNNH] * c[sH];
	w[rNO210f] = k[rNO210f] * c[sNNH] * c[sOH];
	w[rNO211f] = k[rNO211f] * c[sNNH] * c[sCH3];
	w[rNO213f] = k[rNO213f] * c[sHNO] * c[sO];
	w[rNO217f] = k[rNO217f] * c[sCN] * c[sO];
	w[rNO218f] = k[rNO218f] * c[sCN] * c[sOH];
	w[rNO222f] = k[rNO222f] * c[sNCO] * c[sO];
	w[rNO223f] = k[rNO223f] * c[sNCO] * c[sH];
	w[rNO224f] = k[rNO224f] * c[sNCO] * c[sOH];
	w[rNO225f] = k[rNO225f] * c[sNCO] * c[sN];
	w[rNO244f] = k[rNO244f] * c[sC] * c[sNO];
	w[rNO245f] = k[rNO245f] * c[sC] * c[sNO];
	w[rNO246f] = k[rNO246f] * c[sCH] * c[sNO];
	w[rNO247f] = k[rNO247f] * c[sCH] * c[sNO];
	w[rNO248f] = k[rNO248f] * c[sCH] * c[sNO];
	w[rNO257f] = k[rNO257f] * c[sHCNN] * c[sO];
	w[rNO258f] = k[rNO258f] * c[sHCNN] * c[sO];
	w[rNO259f] = k[rNO259f] * c[sHCNN] * c[sO2];
	w[rNO260f] = k[rNO260f] * c[sHCNN] * c[sOH];
	w[rNO261f] = k[rNO261f] * c[sHCNN] * c[sH];
	w[rNO274f] = k[rNO274f] * c[sHCCO] * c[sNO];
	w[r15f] = k[r15f] * c[sH] * c[sO2];
	w[r15b] = k[r15b] * c[sHO2];
	w[r22] = k[r22] * c[sH2O2];
	w[rG14] = k[rG14] * c[sTXCH2] * c[sH];
	w[rLPG60f] = k[rLPG60f] * c[sCH3] * c[sO2];
	w[rLPG60b] = k[rLPG60b] * c[sCH3O2];
	w[r27] = k[r27] * c[sCO] * c[sO];
	w[rLPG24] = k[rLPG24] * c[sTXCH2] * c[sCO];
	w[rLPG44] = k[rLPG44] * c[sCH3O];
	w[rG115f] = k[rG115f] * c[sC2H2] * c[sH];
	w[rG115b] = k[rG115b] * c[sC2H3];
	w[rG155f] = k[rG155f] * c[sC2H4] * c[sH];
	w[rG155b] = k[rG155b] * c[sC2H5];
	w[rLPDD76f] = k[rLPDD76f] * c[sCH2CHO];
	w[rLPDD76b] = k[rLPDD76b] * c[sCO] * c[sCH3];
	w[rG183] = k[rG183] * c[sNXC3H7];
	w[rG163] = k[rG163] * c[sNXC3H7];
	w[rT23f] = k[rT23f] * c[sA1CH2OXC7H7O] * c[sH];
	w[rT23b] = k[rT23b] * c[sA1CH2OHXC7H8O];
	w[rOX01f] = k[rOX01f] * c[sA1XXC6H5];
	w[rOX01b] = k[rOX01b] * c[sH] * c[sOXC6H4];
	w[rNO185f] = k[rNO185f] * c[sN2O];
	w[rNO185b] = k[rNO185b] * c[sO] * c[sN2];
	w[rNO237f] = k[rNO237f] * c[sH] * c[sHCN];
	w[rNO237b] = k[rNO237b] * c[sH2CN];
	w[rNO241f] = k[rNO241f] * c[sCH] * c[sN2];
	w[rNO241b] = k[rNO241b] * c[sHCNN];


	cdot[sN2] = - w[rG28f] + w[rG28f] - w[rG28b] + w[rG28b]
		+ w[rNO178f] - w[rNO178b] + w[rNO181f] - w[rNO181b]
		+ w[rNO183f] - w[rNO183b] + w[rNO184f] - w[rNO184b]
		- w[rNO196b] + w[rNO198f] - w[rNO198b] - w[rNO204b]
		+ w[rNO205f] - w[rNO205b] - w[rNO206b] - w[rNO207b]
		- w[rNO209b] - w[rNO210b] - w[rNO211b] - w[rNO225b]
		+ w[rNO229f] - w[rNO229b] + w[rNO238f] - w[rNO238b]
		- w[rNO239f] + w[rNO239b] - w[rNO240f] + w[rNO240b]
		- w[rNO242f] + w[rNO242b] - w[rNO243f] + w[rNO243b]
		- w[rNO257b] - w[rNO259b] - w[rNO260b] - w[rNO261b]
		+ w[rNO196f] + w[rNO204f] + w[rNO206f] + w[rNO207f]
		+ w[rNO209f] + w[rNO210f] + w[rNO211f] + w[rNO225f]
		+ w[rNO257f] + w[rNO259f] + w[rNO260f] + w[rNO261f]
		+ w[rNO185f] - w[rNO185b] - w[rNO241f] + w[rNO241b];

	cdot[sAR] = 0.0;

	cdot[sO] = - w[r3f] + w[r3b] + w[r5f] - w[r5b]
		+ w[r1f] - w[r1b] - w[r18] + w[rG20]
		- w[rG46] - w[rG116] - w[rG117] + w[rLPG134]
		- w[rLPG131] - w[rG158] - w[rG159] - w[rLPG157]
		- w[rLPG138] + w[rLPR133] - w[rLPLPR113] - w[rR143]
		- w[rLPR142] - w[rLPHP54] - w[rLPHP52] - w[rLPHP53]
		- w[rDD424] - w[rCP05f] + w[rCP05b] - w[rCP12f]
		+ w[rCP12b] - w[rCP13] + w[rCP19b] - w[rCP33f]
		+ w[rCP33b] - w[rI05f] + w[rI05b] - w[rI12]
		+ w[rI19b] - w[rT11f] + w[rT11b] - w[rT12f]
		+ w[rT12b] - w[rT13f] + w[rT13b] + w[rT17b]
		- w[rT25f] + w[rT25b] + w[rT33b] - w[rT38]
		- w[rT45f] + w[rT45b] - w[rT52f] + w[rT52b]
		+ w[rT55b] + w[rT58f] - w[rT58b] - w[rE07f]
		+ w[rE07b] - w[rE17f] + w[rE17b] - w[rE18f]
		+ w[rE18b] - w[rST02f] + w[rST02b] - w[rST10f]
		+ w[rST10b] - w[rST11f] + w[rST11b] + w[rST13]
		- w[rXY03f] + w[rXY03b] - w[rXY11] - w[rXY14f]
		+ w[rXY14b] - w[rXY15f] + w[rXY15b] - w[rXY16f]
		+ w[rXY16b] - w[rXY27f] + w[rXY27b] - w[rXY34]
		- w[rXY43f] + w[rXY43b] - w[rXY44] - w[rXY45f]
		+ w[rXY45b] - w[rXY52] + w[rXY52] - w[rN09f]
		+ w[rN09b] - w[rN14] - w[rN15] + w[rN16b]
		+ w[rN25b] - w[rN30] - w[rOX99f] + w[rOX99b]
		- w[rOX07f] + w[rOX07b] - w[rOX08f] + w[rOX08b]
		- w[rOX09f] + w[rOX09b] - w[rOX10f] + w[rOX10b]
		+ w[rOX11f] - w[rOX11b] + w[rOX13b] - w[rOX24f]
		+ w[rOX24b] - w[rOX28] - w[rOX30f] + w[rOX30b]
		- w[rOX31f] + w[rOX31b] + w[rOX33f] - w[rOX33b]
		+ w[rOX34f] - w[rOX34b] + w[rOX37b] + w[rOX38b]
		+ w[rNO178f] - w[rNO178b] + w[rNO179f] - w[rNO179b]
		- w[rNO181f] + w[rNO181b] - w[rNO182f] + w[rNO182b]
		- w[rNO187f] + w[rNO187b] - w[rNO188f] + w[rNO188b]
		+ w[rNO190b] + w[rNO194f] - w[rNO194b] + w[rNO200b]
		+ w[rNO201b] + w[rNO207b] + w[rNO208b] + w[rNO213b]
		+ w[rNO217b] + w[rNO220f] - w[rNO220b] + w[rNO222b]
		- w[rNO231f] + w[rNO231b] - w[rNO232f] + w[rNO232b]
		- w[rNO233f] + w[rNO233b] - w[rNO244b] - w[rNO246b]
		+ w[rNO257b] + w[rNO258b] - w[rNO259b] - w[rNO262f]
		+ w[rNO262b] - w[rNO263f] + w[rNO263b] - w[rNO264f]
		+ w[rNO264b] - w[rNO279f] + w[rNO279b] + w[rG09]
		- w[rG15] - w[rG52] - w[rG53] - w[rG110]
		- w[rCP19f] - w[rCP26] - w[rI19f] - w[rT17f]
		- w[rT33f] - w[rT55f] - w[rN16f] - w[rN25f]
		- w[rOX13f] - w[rOX37f] - w[rOX38f] - w[rOX47]
		- w[rNO190f] - w[rNO200f] - w[rNO201f] - w[rNO207f]
		- w[rNO208f] - w[rNO213f] - w[rNO217f] - w[rNO222f]
		+ w[rNO244f] + w[rNO246f] - w[rNO257f] - w[rNO258f]
		+ w[rNO259f] - w[r27] + w[rNO185f] - w[rNO185b];

	cdot[sH2] = - w[r3f] + w[r3b] - w[r4f] + w[r4b]
		+ w[r16] - w[rG06f] + w[rG06b] + w[rG34b]
		+ w[rG45] + w[rDD50] + w[rG156] + w[rLPG123f]
		- w[rLPG123b] + w[rLPR153] + w[rDD432] + w[rK017f]
		- w[rK017b] + w[rK021f] - w[rK021b] + w[rK024f]
		- w[rK024b] + w[rP027f] - w[rP027b] + w[rK109f]
		- w[rK109b] + w[rK112f] - w[rK112b] + w[rP118f]
		- w[rP118b] + w[rP119f] - w[rP119b] + w[rP120f]
		- w[rP120b] + w[rP121f] - w[rP121b] + w[rK123f]
		- w[rK123b] + w[rK126f] - w[rK126b] + w[rK129f]
		- w[rK129b] + w[rK132f] - w[rK132b] + w[rK203f]
		- w[rK203b] + w[rP207f] - w[rP207b] + w[rP208f]
		- w[rP208b] + w[rK210f] - w[rK210b] + w[rK213f]
		- w[rK213b] + w[rP305f] - w[rP305b] + w[rK420f]
		- w[rK420b] + w[rK423f] - w[rK423b] - w[rP504b]
		+ w[rK507f] - w[rK507b] + w[rK602f] - w[rK602b]
		+ w[rP800] + w[rP801] + w[rP802] + w[rCP02f]
		- w[rCP02b] + w[rI02f] - w[rI02b] + w[rT07f]
		- w[rT07b] + w[rT24f] - w[rT24b] - w[rT32b]
		+ w[rT37] + w[rT44f] - w[rT44b] + w[rT51f]
		- w[rT51b] + w[rE06f] - w[rE06b] - w[rE13b]
		+ w[rXY02f] - w[rXY02b] + w[rXY26f] - w[rXY26b]
		+ w[rXY33] + w[rN08f] - w[rN08b] - w[rN24b]
		+ w[rN29] + w[rOX03f] - w[rOX03b] + w[rOX19f]
		- w[rOX19b] + w[rOX43f] - w[rOX43b] + w[rNO191f]
		- w[rNO191b] + w[rNO197f] - w[rNO197b] + w[rNO202f]
		- w[rNO202b] - w[rNO209b] + w[rNO214f] - w[rNO214b]
		- w[rNO221f] + w[rNO221b] + w[rNO266f] - w[rNO266b]
		+ w[rNO276f] - w[rNO276b] + w[rNO277f] - w[rNO277b]
		- w[rG34f] + w[rG53] + w[r32] + w[rG130]
		+ w[rP504f] + w[rCP25] + w[rCP29] + w[rT32f]
		+ w[rE13f] + w[rXY17] + w[rXY46] + w[rN24f]
		+ w[rNO209f];

	cdot[sH] = w[r3f] - w[r3b] + w[r4f] - w[r4b]
		- w[r1f] + w[r1b] - w[r16] - w[r17]
		- w[r23] + w[rG06f] - w[rG06b] + w[rG08]
		+ 2 * w[rG19] + w[rG21] - w[rG34b] + w[rG74f]
		- w[rG74b] + w[r29] + w[r28f] - w[r28b]
		+ w[r36f] - w[r36b] - w[rG45] + w[rDD44]
		- w[rDD50] + w[rDD50] + w[rDD52] + w[rLPDD59]
		+ w[rG116] + w[rLPG120f] - w[rLPG120b] + w[rLPG157]
		- w[rG156] - w[rLPG124] - w[rLPG123f] + w[rLPG123b]
		+ w[rLPG136] - w[rLPR015] - w[rLPR016] - w[rR141]
		- w[rLPR153] - w[rLPDD514] + w[rLPHP39] - w[rDD212]
		- w[rHP15] - w[rDD432] + w[rDD405f] - w[rDD405b]
		- w[rP000f] + w[rP000f] - w[rP000b] + w[rP000b]
		+ w[rP001f] - w[rP001b] + w[rP002f] - w[rP002b]
		+ w[rP004f] - w[rP004b] + w[rP006f] - w[rP006b]
		+ w[rP007f] - w[rP007b] + 2 * w[rP008] + w[rP011f]
		- w[rP011b] + w[rP014f] - w[rP014b] + w[rP016f]
		- w[rP016b] - w[rK017f] + w[rK017b] + w[rK020f]
		- w[rK020b] - w[rK021f] + w[rK021b] + w[rP023f]
		- w[rP023b] - w[rK024f] + w[rK024b] + w[rP026f]
		- w[rP026b] - w[rP027f] + w[rP027b] + w[rK102f]
		- w[rK102b] + w[rP104f] - w[rP104b] + w[rP105f]
		- w[rP105b] + w[rP106f] - w[rP106b] + w[rP107f]
		- w[rP107b] + w[rP108f] - w[rP108b] - w[rK109f]
		+ w[rK109b] + w[rP111f] - w[rP111b] - w[rK112f]
		+ w[rK112b] - w[rP116b] - w[rP117b] + w[rK122f]
		- w[rK122b] - w[rK123f] + w[rK123b] + w[rK125f]
		- w[rK125b] - w[rK126f] + w[rK126b] + w[rP128f]
		- w[rP128b] - w[rK129f] + w[rK129b] + w[rP131f]
		- w[rP131b] - w[rK132f] + w[rK132b] + w[rK200f]
		- w[rK200b] - w[rP201f] + w[rP201f] - w[rP201b]
		+ w[rP201b] + w[rP202f] - w[rP202b] - w[rK203f]
		+ w[rK203b] - w[rP206b] + w[rP209f] - w[rP209b]
		- w[rK210f] + w[rK210b] + w[rK212f] - w[rK212b]
		- w[rK213f] + w[rK213b] + w[rP301f] - w[rP301b]
		+ w[rP304f] - w[rP304b] - w[rP305f] + w[rP305b]
		+ w[rP405f] - w[rP405b] + w[rP406f] - w[rP406b]
		+ w[rP407f] - w[rP407b] + w[rP410f] - w[rP410b]
		+ w[rP411f] - w[rP411b] + w[rP412f] - w[rP412b]
		+ w[rP413f] - w[rP413b] + w[rP414f] - w[rP414b]
		+ w[rP415f] - w[rP415b] + w[rP416f] - w[rP416b]
		+ w[rP417f] - w[rP417b] + w[rP419f] - w[rP419b]
		- w[rK420f] + w[rK420b] + w[rP422f] - w[rP422b]
		- w[rK423f] + w[rK423b] + w[rP502f] - w[rP502b]
		+ w[rP503f] - w[rP503b] + w[rP505f] - w[rP505b]
		+ w[rP506f] - w[rP506b] - w[rK507f] + w[rK507b]
		+ w[rK600f] - w[rK600b] + w[rP601f] - w[rP601b]
		- w[rK602f] + w[rK602b] + w[rK700f] - w[rK700b]
		+ w[rK701f] - w[rK701b] + w[rP800] + w[rP801]
		+ w[rCP01f] - w[rCP01b] - w[rCP02f] + w[rCP02b]
		- w[rCP03f] + w[rCP03b] - w[rCP04] + w[rCP12f]
		- w[rCP12b] + 2 * w[rCP17] + 2 * w[rCP18] - w[rCP19b]
		+ w[rCP23f] - w[rCP23b] - w[rCP32f] + w[rCP32b]
		+ w[rI01f] - w[rI01b] - w[rI02f] + w[rI02b]
		+ w[rI03f] - w[rI03b] + 2 * w[rI12] + 2 * w[rI17]
		+ 2 * w[rI18] - w[rI19b] + 2 * w[rI23] - w[rI32]
		- w[rT01f] + w[rT01b] + w[rT02f] - w[rT02b]
		- w[rT04f] + w[rT04b] - w[rT07f] + w[rT07b]
		+ w[rT10f] - w[rT10b] + w[rT13f] - w[rT13b]
		+ 2 * w[rT20] - w[rT24f] + w[rT24b] + w[rT29f]
		- w[rT29b] + w[rT32b] + w[rT36] - w[rT37]
		+ w[rT43f] - w[rT43b] - w[rT44f] + w[rT44b]
		+ w[rT47] + w[rT50f] - w[rT50b] - w[rT51f]
		+ w[rT51b] - w[rT56b] + w[rT59] + w[rE01b]
		- w[rE04f] + w[rE04b] - w[rE06f] + w[rE06b]
		+ w[rE12f] - w[rE12b] + w[rE13b] + w[rXY00f]
		- w[rXY00b] - w[rXY02f] + w[rXY02b] - w[rXY09f]
		+ w[rXY09b] + 2 * w[rXY11] + w[rXY12] - w[rXY13f]
		+ w[rXY13b] + w[rXY14f] - w[rXY14b] + w[rXY201]
		+ 2 * w[rXY22] + w[rXY23f] - w[rXY23b] + w[rXY25]
		- w[rXY26f] + w[rXY26b] - w[rXY33] - w[rXY39f]
		+ w[rXY39b] - w[rXY40f] + w[rXY40b] + w[rXY42]
		+ w[rXY43f] - w[rXY43b] + w[rXY48] + w[rXY50]
		- w[rXY51] + w[rXY51] - w[rXY57f] + w[rXY57b]
		+ w[rXY58] - w[rN01f] + w[rN01b] + w[rN04f]
		- w[rN04b] - w[rN06f] + w[rN06b] - w[rN08f]
		+ w[rN08b] + 2 * w[rN14] + 2 * w[rN19] + w[rN22f]
		- w[rN22b] + w[rN24b] + w[rN28] - w[rN29]
		+ w[rOX00f] - w[rOX00b] - w[rOX03f] + w[rOX03b]
		+ w[rOX05f] - w[rOX05b] + w[rOX99f] - w[rOX99b]
		+ w[rOX07f] - w[rOX07b] + w[rOX12f] - w[rOX12b]
		- w[rOX14b] + w[rOX18f] - w[rOX18b] - w[rOX19f]
		+ w[rOX19b] + w[rOX24f] - w[rOX24b] - w[rOX27f]
		+ w[rOX27b] + w[rOX30f] - w[rOX30b] + w[rOX32f]
		- w[rOX32b] + w[rOX35] + w[rOX36] - w[rOX39b]
		- w[rOX40b] + w[rOX42f] - w[rOX42b] - w[rOX43f]
		+ w[rOX43b] + w[rNO180f] - w[rNO180b] - w[rNO183f]
		+ w[rNO183b] - w[rNO189f] + w[rNO189b] - w[rNO190b]
		- w[rNO191f] + w[rNO191b] - w[rNO192b] - w[rNO196b]
		+ w[rNO199f] - w[rNO199b] - w[rNO201b] - w[rNO202f]
		+ w[rNO202b] - w[rNO204b] + w[rNO205f] - w[rNO205b]
		+ w[rNO209b] - w[rNO212f] + w[rNO212b] - w[rNO214f]
		+ w[rNO214b] - w[rNO218b] + w[rNO221f] - w[rNO221b]
		+ w[rNO223b] - w[rNO224b] + w[rNO230f] - w[rNO230b]
		+ w[rNO231f] - w[rNO231b] + w[rNO234f] - w[rNO234b]
		+ w[rNO235] - w[rNO247b] + w[rNO249f] - w[rNO249b]
		+ w[rNO251f] - w[rNO251b] + w[rNO252f] - w[rNO252b]
		+ w[rNO254f] - w[rNO254b] - w[rNO257b] - w[rNO260b]
		+ w[rNO261b] - w[rNO265f] + w[rNO265b] - w[rNO266f]
		+ w[rNO266b] - w[rNO270f] + w[rNO270f] - w[rNO270b]
		+ w[rNO270b] - w[rNO271f] + w[rNO271b] - w[rNO272f]
		+ w[rNO272b] - w[rNO273] + w[rNO273] + w[rNO275f]
		- w[rNO275b] - w[rNO277f] + w[rNO277b] + w[rG15]
		+ w[rG16] + w[rG35] + w[rG34f] + w[rG52]
		+ w[rG53] + w[rG72] - w[r32] - w[rG130]
		- w[rG109] + w[rG110] + w[rP116f] + w[rP117f]
		+ w[rP206f] + w[rCP19f] + w[rCP22] - w[rCP25]
		- w[rCP29] + w[rI19f] + 2 * w[rI22] - w[rT32f]
		+ w[rT56f] - w[rE01f] - w[rE13f] + w[rN17]
		- w[rN24f] + w[rOX14f] + w[rOX39f] + w[rOX40f]
		+ w[rOX47] + w[rNO190f] + w[rNO192f] + w[rNO196f]
		+ w[rNO201f] + w[rNO204f] - w[rNO209f] + w[rNO218f]
		- w[rNO223f] + w[rNO224f] + w[rNO247f] + w[rNO257f]
		+ w[rNO260f] - w[rNO261f] - w[r15f] + w[r15b]
		- w[rG14] + w[rLPG44] - w[rG115f] + w[rG115b]
		- w[rG155f] + w[rG155b] + w[rG183] - w[rT23f]
		+ w[rT23b] + w[rOX01f] - w[rOX01b] - w[rNO237f]
		+ w[rNO237b];

	cdot[sOH] = w[r3f] - w[r3b] - w[r4f] + w[r4b]
		- 2 * w[r5f] + 2 * w[r5b] + w[r1f] - w[r1b]
		- w[r19] + w[r18] + 2 * w[r17] - w[r26]
		+ w[r23] - w[r39] - w[rG17] + w[rG21]
		- w[rG55] - w[rG57f] + w[rG57b] + w[rG59]
		- w[r29] + w[r31] - w[r28f] + w[r28b]
		+ w[rG46] - w[rG47] - w[rLPDD62] - w[rDD47]
		+ w[rDD47] + w[rDD48] + w[rDD53] - w[rDD52]
		+ w[rLPG63] + w[rLPG64] + w[rDD66] - w[rLPG120f]
		+ w[rLPG120b] - w[rG122] - w[rG160] - w[rG161]
		+ w[rG111] - w[rLPG127f] + w[rLPG127b] - w[rLPR025]
		+ w[rLPG139] + w[rLPG140] - w[rLPG144] + w[rLPR111]
		+ w[rLPLPR117] - w[rLPR151] - w[rLPHP55] - w[rLPHP56]
		- w[rDD199] - w[rDD43] - w[rHP26] - w[rDD434]
		+ w[rDD424] + w[rDD427] + w[rLPDD431] + w[rLPDD502]
		+ w[rLPLPDD495] + w[rLPDD499] + w[rLPDD490] - w[rDD332N]
		- w[rDD294] - w[rHP60] + w[rDD259] - w[rP018f]
		+ w[rP018b] + w[rP022b] - w[rP025f] + w[rP025b]
		- w[rP028f] + w[rP028b] - w[rK110f] + w[rK110b]
		- w[rK113f] + w[rK113b] + w[rP124b] + w[rP127b]
		- w[rP130f] + w[rP130b] - w[rP133f] + w[rP133b]
		- w[rK204f] + w[rK204b] - w[rP211f] + w[rP211b]
		+ w[rP214b] - w[rP306f] + w[rP306b] - w[rK421f]
		+ w[rK421b] - w[rK424f] + w[rK424b] - w[rK508f]
		+ w[rK508b] - w[rK603f] + w[rK603b] + w[rCP05f]
		- w[rCP05b] - w[rCP06f] + w[rCP06b] - w[rCP14]
		- w[rCP15] + w[rCP20f] - w[rCP20b] - w[rCP21b]
		+ w[rI05f] - w[rI05b] - w[rI06f] + w[rI06b]
		- w[rI15] + w[rI20f] - w[rI20b] - w[rT08f]
		+ w[rT08b] - w[rT09f] + w[rT09b] - w[rT10f]
		+ w[rT10b] + w[rT11f] - w[rT11b] + w[rT18b]
		+ w[rT19f] - w[rT19b] + w[rT21f] - w[rT21b]
		+ w[rT25f] - w[rT25b] - w[rT26f] + w[rT26b]
		- w[rT33b] + w[rT34b] + w[rT38] - w[rT39]
		+ w[rT45f] - w[rT45b] - w[rT46f] + w[rT46b]
		+ w[rT52f] - w[rT52b] - w[rT53f] + w[rT53b]
		+ w[rT56b] - w[rT57b] - w[rE05f] + w[rE05b]
		+ w[rE07f] - w[rE07b] - w[rE08f] + w[rE08b]
		+ w[rE14b] + w[rE19] + w[rE31] + w[rE34]
		+ w[rE37] - w[rST04f] + w[rST04b] - w[rST05f]
		+ w[rST05b] - w[rST06f] + w[rST06b] + w[rST10f]
		- w[rST10b] + w[rST00b] + w[rXY03f] - w[rXY03b]
		- w[rXY04f] + w[rXY04b] - w[rXY10f] + w[rXY10b]
		+ w[rXY18f] - w[rXY18b] + w[rXY201] + w[rXY202]
		+ w[rXY203] + w[rXY27f] - w[rXY27b] - w[rXY28f]
		+ w[rXY28b] + w[rXY34] - w[rXY35] - w[rXY41f]
		+ w[rXY41b] - w[rXY42] + w[rXY47f] - w[rXY47b]
		+ w[rXY48] - w[rXY53] + w[rXY53] - w[rXY58]
		- w[rN02f] + w[rN02b] + w[rN09f] - w[rN09b]
		- w[rN10f] + w[rN10b] + w[rN18f] - w[rN18b]
		+ w[rN20f] - w[rN20b] - w[rN25b] + w[rN26b]
		+ w[rN30] - w[rN31] - w[rOX04f] + w[rOX04b]
		- w[rOX05f] + w[rOX05b] + w[rOX10f] - w[rOX10b]
		+ w[rOX14b] - w[rOX15b] - w[rOX20f] + w[rOX20b]
		+ w[rOX25f] - w[rOX25b] - w[rOX32f] + w[rOX32b]
		+ w[rOX39b] + w[rOX40b] - w[rOX44f] + w[rOX44b]
		+ w[rOX48] - w[rOX60] - w[rOX61] - w[rOX62]
		- w[rOX63] - w[rOX64] - w[rOX65] - w[rOX66]
		- w[rNO180f] + w[rNO180b] + w[rNO183f] - w[rNO183b]
		- w[rNO184f] + w[rNO184b] + w[rNO186f] - w[rNO186b]
		+ w[rNO189f] - w[rNO189b] + w[rNO192b] - w[rNO193f]
		+ w[rNO193b] + w[rNO195f] - w[rNO195b] + w[rNO198f]
		- w[rNO198b] - w[rNO200b] - w[rNO203f] + w[rNO203b]
		- w[rNO207b] + w[rNO210b] - w[rNO213b] - w[rNO215f]
		+ w[rNO215b] + w[rNO218b] + w[rNO219f] - w[rNO219b]
		+ w[rNO224b] + w[rNO233f] - w[rNO233b] - w[rNO234f]
		+ w[rNO234b] - w[rNO235] - w[rNO236f] + w[rNO236b]
		+ w[rNO250f] - w[rNO250b] + w[rNO253f] - w[rNO253b]
		+ w[rNO256f] - w[rNO256b] + w[rNO260b] + w[rNO264f]
		- w[rNO264b] - w[rNO267f] + w[rNO267b] - w[rNO268f]
		+ w[rNO268b] + w[rNO271f] - w[rNO271b] - w[rNO278f]
		+ w[rNO278b] + w[rNO279f] - w[rNO279b] - w[rG16]
		+ w[rG35] + w[rLPG65] - w[r35] - w[rDD54]
		+ w[rR037] - w[rR024] - w[rP022f] - w[rP124f]
		- w[rP127f] - w[rP214f] + w[rCP21f] - w[rCP22]
		+ w[rCP26] - w[rCP27] - w[rI22] - w[rT18f]
		+ w[rT33f] - w[rT34f] - w[rT56f] + w[rT57f]
		- w[rE14f] + w[rE36] - w[rST00f] - w[rXY17]
		- w[rXY46] - w[rN17] + w[rN25f] - w[rN26f]
		- w[rOX14f] + w[rOX15f] - w[rOX39f] - w[rOX40f]
		- w[rNO192f] + w[rNO200f] + w[rNO207f] - w[rNO210f]
		+ w[rNO213f] - w[rNO218f] - w[rNO224f] - w[rNO260f]
		+ 2 * w[r22];

	cdot[sH2O] = w[r4f] - w[r4b] + w[r5f] - w[r5b]
		+ w[r19] + w[r26] + w[r23] + w[r39]
		- w[rG08] + w[rG17] - w[rG38b] + w[rG38b]
		+ w[rG55] + w[rG57f] - w[rG57b] + w[rG47]
		+ w[rDD47] + w[rDD52] + w[rG160] + w[rLPG127f]
		- w[rLPG127b] + w[rLPR151] + w[rDD434] + w[rDD332N]
		+ w[rDD294] + w[rHP60] + w[rP018f] - w[rP018b]
		- w[rP022b] + w[rP025f] - w[rP025b] + w[rP028f]
		- w[rP028b] + w[rK110f] - w[rK110b] + w[rK113f]
		- w[rK113b] - w[rP124b] - w[rP127b] + w[rP130f]
		- w[rP130b] + w[rP133f] - w[rP133b] + w[rK204f]
		- w[rK204b] + w[rP211f] - w[rP211b] - w[rP214b]
		+ w[rP306f] - w[rP306b] + w[rK421f] - w[rK421b]
		+ w[rK424f] - w[rK424b] + w[rK508f] - w[rK508b]
		+ w[rK603f] - w[rK603b] + w[rCP06f] - w[rCP06b]
		+ w[rI06f] - w[rI06b] + w[rT08f] - w[rT08b]
		+ w[rT26f] - w[rT26b] - w[rT34b] + w[rT39]
		+ w[rT46f] - w[rT46b] + w[rT53f] - w[rT53b]
		+ w[rE08f] - w[rE08b] - w[rE14b] + w[rXY04f]
		- w[rXY04b] + w[rXY28f] - w[rXY28b] + w[rXY35]
		+ w[rN10f] - w[rN10b] - w[rN26b] + w[rN31]
		+ w[rOX04f] - w[rOX04b] + w[rOX20f] - w[rOX20b]
		+ w[rOX44f] - w[rOX44b] + w[rNO193f] - w[rNO193b]
		- w[rNO197f] + w[rNO197b] + w[rNO203f] - w[rNO203b]
		- w[rNO210b] + w[rNO215f] - w[rNO215b] - w[rNO219f]
		+ w[rNO219b] + w[rNO255f] - w[rNO255b] + w[rNO267f]
		- w[rNO267b] + w[rNO278f] - w[rNO278b] + w[rG36]
		- w[rG38f] + w[rG38f] + w[r35] + w[rP022f]
		+ w[rP124f] + w[rP127f] + w[rP214f] + w[rCP27]
		+ w[rI21] + w[rT34f] + w[rE14f] + w[rN26f]
		+ w[rNO210f];

	cdot[sO2] = - w[r1f] + w[r1b] + w[r16] + w[r21]
		+ w[r20] + w[r19] + w[r18] - w[rG19]
		- w[rG20] - w[rG21] - w[rG59] - w[r38]
		- w[rDD68] - w[rLPG89] + w[rLPG63] + w[rDD67]
		- w[rG135] - w[rLPG134] - w[rG133] - w[rG171]
		- w[rG111] - w[rLPG139] - w[rLPG140] + w[rLPR116]
		- w[rLPR133] - w[rLPR111] - w[rG184] - w[rDD422]
		- w[rDD408f] + w[rDD408b] - w[rLPDD486f] + w[rLPDD486b]
		- w[rDD239] - w[rDD263] - w[rCP07f] + w[rCP07b]
		- w[rCP20f] + w[rCP20b] - w[rCP28] - w[rCP30]
		- w[rI07f] + w[rI07b] - w[rI20f] + w[rI20b]
		- w[rT06f] + w[rT06b] - w[rT21f] + w[rT21b]
		- w[rT22f] + w[rT22b] - w[rT35f] + w[rT35b]
		- w[rT40] - w[rT58f] + w[rT58b] - w[rT59]
		- w[rT60] - w[rE16f] + w[rE16b] - w[rE30]
		- w[rE31] + w[rE32] - w[rST12f] + w[rST12b]
		- w[rST13] - w[rST14f] + w[rST14b] - w[rXY05f]
		+ w[rXY05b] - w[rXY18f] + w[rXY18b] - w[rXY19f]
		+ w[rXY19b] - w[rXY29f] + w[rXY29b] - w[rXY36]
		- w[rXY47f] + w[rXY47b] - w[rXY54] + w[rXY54]
		- w[rN11f] + w[rN11b] - w[rN20f] + w[rN20b]
		- w[rN21f] + w[rN21b] - w[rN27f] + w[rN27b]
		- w[rN32] - w[rOX06f] + w[rOX06b] - w[rOX11f]
		+ w[rOX11b] - w[rOX12f] + w[rOX12b] - w[rOX22f]
		+ w[rOX22b] - w[rOX25f] + w[rOX25b] - w[rOX33f]
		+ w[rOX33b] - w[rOX34f] + w[rOX34b] - w[rOX35]
		- w[rOX36] - w[rOX48] - w[rOX50] - w[rOX51]
		- w[rOX52] - w[rOX53] - w[rOX54] - w[rNO179f]
		+ w[rNO179b] + w[rNO181f] - w[rNO181b] + w[rNO188f]
		- w[rNO188b] - w[rNO194f] + w[rNO194b] - w[rNO195f]
		+ w[rNO195b] + w[rNO206b] - w[rNO216f] + w[rNO216b]
		- w[rNO220f] + w[rNO220b] - w[rNO226f] + w[rNO226b]
		+ w[rNO259b] - w[rG09] - w[rG36] - w[rG35]
		- w[rE36] - w[rNO206f] - w[rNO259f] - w[r15f]
		+ w[r15b] - w[rLPG60f] + w[rLPG60b];

	cdot[sHO2] = - w[r16] - 2 * w[r21] - 2 * w[r20] - w[r19]
		- w[r18] - w[r17] + w[r26] + w[r39]
		- w[r31] + w[r38] - w[rG49] + w[rLPG89]
		- w[rDD48] - w[rLPG63] + w[rG133] + w[rG171]
		- w[rLPR116] - w[rLPLPR117] + w[rG184] + w[rDD422]
		- w[rDD429] + w[rDD472] + w[rLPDD473] + w[rDD239]
		+ w[rCP07f] - w[rCP07b] - w[rCP08f] + w[rCP08b]
		+ w[rCP21b] + w[rCP28] + w[rCP30] + w[rI07f]
		- w[rI07b] - w[rI08f] + w[rI08b] + w[rT06f]
		- w[rT06b] - w[rT15f] + w[rT15b] - w[rT19f]
		+ w[rT19b] + w[rT35f] - w[rT35b] + w[rT40]
		- w[rT41] + w[rT57b] - w[rE09f] + w[rE09b]
		+ w[rE16f] - w[rE16b] - w[rE19] + w[rE33]
		+ w[rST12f] - w[rST12b] + w[rXY05f] - w[rXY05b]
		- w[rXY06f] + w[rXY06b] - w[rXY201] - w[rXY202]
		- w[rXY203] + w[rXY29f] - w[rXY29b] - w[rXY30f]
		+ w[rXY30b] + w[rXY36] - w[rXY37] - w[rXY48]
		- w[rXY55] + w[rXY55] + w[rN11f] - w[rN11b]
		- w[rN13f] + w[rN13b] - w[rN18f] + w[rN18b]
		+ w[rN27f] - w[rN27b] + w[rN32] - w[rN33]
		+ w[rOX06f] - w[rOX06b] + w[rOX15b] + w[rOX22f]
		- w[rOX22b] + w[rNO184f] - w[rNO184b] - w[rNO186f]
		+ w[rNO186b] - w[rNO206b] + w[rNO216f] - w[rNO216b]
		- w[rLPG65] - w[rR037] - w[rCP21f] - w[rI21]
		- w[rT57f] - w[rOX15f] + w[rNO206f] + w[r15f]
		- w[r15b];

	cdot[sH2O2] = w[r21] + w[r20] - w[r26] - w[r23]
		- w[r39] + w[rG49] + w[rDD48] + w[rDD429]
		+ w[rCP08f] - w[rCP08b] + w[rI08f] - w[rI08b]
		+ w[rT15f] - w[rT15b] + w[rT41] + w[rE09f]
		- w[rE09b] + w[rXY06f] - w[rXY06b] + w[rXY30f]
		- w[rXY30b] + w[rXY37] + w[rN13f] - w[rN13b]
		+ w[rN33] - w[r22];

	cdot[sCH] = - w[rG06f] + w[rG06b] - w[rG08] + w[rG17]
		- w[rNO240f] + w[rNO240b] + w[rNO246b] + w[rNO247b]
		+ w[rNO248b] - w[rG09] - w[rNO246f] - w[rNO247f]
		- w[rNO248f] - w[rNO241f] + w[rNO241b];

	cdot[sTXCH2] = w[rG06f] - w[rG06b] - w[rG19] - w[rG20]
		- w[rG17] - w[rG21] - w[rG38b] + w[rG28f]
		- w[rG28b] + w[rG55] + w[rG117] + w[rG158]
		+ w[rLPHP54] + w[rNO238f] - w[rNO238b] - w[rNO242f]
		+ w[rNO242b] - w[rNO249f] + w[rNO249b] - w[rNO250f]
		+ w[rNO250b] - w[rNO251f] + w[rNO251b] - w[rNO261b]
		- w[rG15] - w[rG16] + w[rG38f] - w[rG72]
		+ w[rNO261f] - w[rG14] - w[rLPG24];

	cdot[sCH2O] = w[rG08] + w[rG20] + w[rG59] - w[rG46]
		- w[rG47] - w[rG45] - w[rG49] - w[rLPDD62]
		+ w[rLPG89] - w[rLPG64] - w[rDD66] + w[rG135]
		+ w[rG158] + w[rLPG139] + w[rLPG138] + w[rG168]
		+ w[rLPR111] + w[rLPHP55] + w[rLPHP53] + w[rDD199]
		+ w[rDD43] + w[rLPR140] + w[rDD259] + w[rT22f]
		- w[rT22b] - w[rT28f] + w[rT28b] + w[rT30f]
		- w[rT30b] + w[rE17f] - w[rE17b] + w[rE19]
		+ w[rE37] + w[rST05f] - w[rST05b] + w[rXY15f]
		- w[rXY15b] + w[rXY19f] - w[rXY19b] + w[rXY202]
		+ w[rXY44] + w[rN21f] - w[rN21b] + w[rN23f]
		- w[rN23b] + w[rG16] + w[rG52] + w[rLPG44];

	cdot[sHCO] = - w[r38] - w[r36f] + w[r36b] - w[rDD68]
		+ w[rG46] + w[rG47] + w[rG45] + w[rG49]
		+ w[rDD53] + w[rLPG64] + w[rDD66] + w[rDD67]
		+ w[rG135] + w[rG159] + w[rLPR030] + 2 * w[rLPG140]
		+ w[rLPG138] + w[rLPG144] - w[rLPR115] + w[rR143]
		+ w[rLPHP54] + w[rLPHP56] + w[rHP26] + w[rCP14]
		+ w[rCP15] + w[rI15] + w[rT28f] - w[rT28b]
		+ w[rT31f] - w[rT31b] + w[rE31] + w[rE34]
		+ w[rST02f] - w[rST02b] + w[rST14f] - w[rST14b]
		+ w[rXY16f] - w[rXY16b] + w[rXY203] + w[rXY39f]
		- w[rXY39b] + w[rXY41f] - w[rXY41b] + w[rXY45f]
		- w[rXY45b] + w[rXY57f] - w[rXY57b] + w[rXY58]
		- w[rNO248b] - w[rNO259b] - w[rNO260b] + w[rG09]
		+ w[rG15] - w[r32] - w[r35] - w[rDD54]
		+ 2 * w[rR024] + w[rNO248f] + w[rNO259f] + w[rNO260f];

	cdot[sCH3] = - w[rG34b] - w[rG55] - w[rG57f] + w[rG57b]
		- w[rG59] - 2 * w[rG74f] + 2 * w[rG74b] - w[rLPG61]
		+ w[rG122] + w[rG159] + w[rLPG124] - w[rLPR028]
		- w[rLPR030] + w[rLPG137] + w[rG168] + w[rLPR015]
		+ w[rLPR133] + w[rLPR007] + w[rR141] + w[rLPR142]
		+ w[rLPHP43] + w[rLPHP52] + w[rLPHP56] + w[rDD282]
		+ w[rLPLPDD107] - w[rCP09f] + w[rCP09b] - w[rCP18]
		+ w[rCP34] - w[rI09f] + w[rI09b] - w[rI18]
		+ w[rT01f] - w[rT01b] + w[rT03f] - w[rT03b]
		+ w[rT04f] - w[rT04b] + w[rT09f] - w[rT09b]
		- w[rT14f] + w[rT14b] - w[rT27f] + w[rT27b]
		- w[rT42] - w[rT54f] + w[rT54b] + w[rE02b]
		- w[rE10f] + w[rE10b] - w[rE15f] + w[rE15b]
		+ w[rE18f] - w[rE18b] - w[rST03f] + w[rST03b]
		+ w[rST04f] - w[rST04b] + w[rXY01f] - w[rXY01b]
		- w[rXY07f] + w[rXY07b] + w[rXY09f] - w[rXY09b]
		+ w[rXY10f] - w[rXY10b] + w[rXY13f] - w[rXY13b]
		+ w[rXY24] - w[rXY31f] + w[rXY31b] - w[rXY38]
		+ w[rXY40f] - w[rXY40b] + w[rXY42] - w[rXY56]
		+ w[rXY56] + w[rN01f] - w[rN01b] + w[rN02f]
		- w[rN02b] + w[rN05f] - w[rN05b] + w[rN06f]
		- w[rN06b] - w[rN12f] + w[rN12b] + w[rN15]
		- w[rN34] + w[rOX16f] - w[rOX16b] - w[rOX21f]
		+ w[rOX21b] - w[rOX45f] + w[rOX45b] + w[rOX60]
		+ w[rOX61] + w[rNO211b] - w[rNO255f] + w[rNO255b]
		- w[rNO256f] + w[rNO256b] - w[rNO275f] + w[rNO275b]
		- w[rNO276f] + w[rNO276b] + w[rG34f] - w[rG52]
		- w[rG53] - w[rG72] - w[rLPG65] - w[rR022]
		- w[rE02f] - w[rNO211f] + w[rG14] - w[rLPG60f]
		+ w[rLPG60b] + w[rLPDD76f] - w[rLPDD76b] + w[rG163];

	cdot[sCO2] = w[rG19] + w[r29] + w[r31] + w[r28f]
		- w[r28b] + w[rDD44] + w[rDD50] + w[rDD52]
		+ w[rCP33f] - w[rCP33b] + w[rT59] + w[rNO226f]
		- w[rNO226b] + w[rNO229f] - w[rNO229b] + w[rNO262f]
		- w[rNO262b] + w[rNO268f] - w[rNO268b] - w[rNO280f]
		+ w[rNO280b] + w[rNO282f] - w[rNO282b] - w[rNO283f]
		+ w[rNO283b] + w[r27];

	cdot[sCO] = w[rG21] - w[r29] - w[r31] - w[r28f]
		+ w[r28b] + w[r38] + w[r36f] - w[r36b]
		+ w[rDD47] + w[rDD48] + w[rG117] + w[rG122]
		+ 2 * w[rG111] + w[rLPG124] + w[rLPR028] + w[rLPR025]
		+ w[rLPG139] + w[rLPG137] + w[rLPR115] + w[rLPHP52]
		+ w[rDD332N] + w[rDD290] + w[rHP60] + w[rCP13]
		+ w[rCP24] + w[rCP31] + w[rI31] + w[rI32]
		+ w[rT36] + w[rT37] + w[rT38] + w[rT39]
		+ w[rT40] + w[rT41] + w[rT42] + w[rT47]
		+ 2 * w[rT60] + w[rE37] + w[rST11f] - w[rST11b]
		+ w[rST13] + w[rXY11] + w[rXY24] + w[rXY25]
		+ w[rXY33] + w[rXY34] + w[rXY35] + w[rXY36]
		+ w[rXY37] + w[rXY38] + w[rXY42] + w[rXY44]
		+ 2 * w[rXY50] + w[rXY51] + w[rXY52] + w[rXY53]
		+ w[rXY54] + w[rXY55] + w[rXY56] + w[rXY58]
		+ w[rN14] + w[rN15] + w[rN28] + w[rN29]
		+ w[rN30] + w[rN31] + w[rN32] + w[rN33]
		+ w[rN34] + w[rOX09f] - w[rOX09b] + w[rOX17f]
		- w[rOX17b] + w[rOX23f] - w[rOX23b] + w[rOX26f]
		- w[rOX26b] + w[rOX27f] - w[rOX27b] + 2 * w[rOX28]
		+ w[rOX35] + w[rOX36] + w[rOX41f] - w[rOX41b]
		+ w[rOX46f] - w[rOX46b] + w[rOX48] + 2 * w[rOX50]
		+ 2 * w[rOX51] + 2 * w[rOX52] + 2 * w[rOX53] + 2 * w[rOX54]
		+ w[rOX60] + w[rOX61] - w[rNO217b] - w[rNO222b]
		- w[rNO223b] - w[rNO224b] - w[rNO225b] + w[rNO227f]
		- w[rNO227b] + w[rNO228f] - w[rNO228b] + w[rNO232f]
		- w[rNO232b] + w[rNO236f] - w[rNO236b] - w[rNO245b]
		- w[rNO257b] + w[rNO263f] - w[rNO263b] + w[rNO265f]
		- w[rNO265b] + w[rNO269f] - w[rNO269b] + w[rNO272f]
		- w[rNO272b] - w[rNO274b] + w[rNO280f] - w[rNO280b]
		+ w[rNO283f] - w[rNO283b] + w[rG36] + w[rG35]
		+ w[rG53] + w[r32] + w[r35] + w[rG109]
		+ w[rR022] + 2 * w[rG110] + w[rOX47] + w[rNO217f]
		+ w[rNO222f] + w[rNO223f] + w[rNO224f] + w[rNO225f]
		+ w[rNO245f] + w[rNO257f] + w[rNO274f] - w[r27]
		- w[rLPG24] + w[rLPDD76f] - w[rLPDD76b];

	cdot[sSXCH2] = w[rG38b] + w[rG34b] - w[rG28f] + w[rG28b]
		+ w[rG57f] - w[rG57b] - w[rNO243f] + w[rNO243b]
		- w[rNO252f] + w[rNO252b] - w[rNO253f] + w[rNO253b]
		- w[rNO254f] + w[rNO254b] - w[rG36] - w[rG35]
		- w[rG38f] - w[rG34f] + w[rG109];

	cdot[sC2H4] = - w[rG160] - w[rG158] - w[rG159] - w[rLPG157]
		- w[rG156] - w[rG161] + w[rG171] + w[rR141]
		+ w[rHP41] + w[rHP29] + w[rDD214] + w[rHP10]
		+ w[rDD287] + w[rDD324] - w[rP015f] + w[rP015b]
		- w[rP106f] + w[rP106b] - w[rP120f] + w[rP120b]
		- w[rP121f] + w[rP121b] - w[rP208f] + w[rP208b]
		- w[rP412f] + w[rP412b] - w[rP413f] + w[rP413b]
		- w[rCP10b] + w[rCP13] + w[rCP15] + w[rI15]
		+ w[rE11f] - w[rE11b] + w[rG72] + w[rR022]
		+ w[rCP10f] - w[rG155f] + w[rG155b] + w[rG163];

	cdot[sCH3O] = - w[rLPG89] + w[rLPG63] + 2 * w[rLPG61] + w[rLPG64]
		+ w[rLPDD15] + w[rLPR025] + w[rLPG144] + w[rLPLPR118]
		+ w[rLPDD431] + w[rLPG65] - w[rLPG44];

	cdot[sCH3O2] = w[rLPDD62] - w[rLPG63] - w[rLPDD59] - w[rLPG61]
		- w[rLPG64] - w[rLPDD15] - w[rLPLPR118] - w[rLPDD431]
		+ w[rLPG60f] - w[rLPG60b];

	cdot[sC2H5] = w[rG74f] - w[rG74b] - w[rLPDD15] - w[rG171]
		+ w[rLPR028] + w[rLPR030] + w[rR143] + w[rLPHP54]
		+ w[rLPHP52] + w[rLPHP56] + w[rHP41] + w[rHP32]
		+ w[rHP26] + w[rDD288] + w[rDD319] + w[rLPDD490]
		+ w[rE03b] + w[rE04f] - w[rE04b] + w[rE05f]
		- w[rE05b] - w[rR037] - w[rE03f] + w[rG155f]
		- w[rG155b];

	cdot[sCH2CO] = w[rLPG120f] - w[rLPG120b] - w[rLPG124] - w[rLPR028]
		- w[rLPG123f] + w[rLPG123b] - w[rLPG127f] + w[rLPG127b]
		- w[rLPR025] + w[rLPG136] + w[rLPR133] + w[rLPDD490]
		+ w[rLPLPDD107] + w[rLPDD293] + w[rOX28] + w[rLPG24];

	cdot[sO2CHO] = w[rDD68] - w[rDD66] - w[rDD67] - w[rDD427];

	cdot[sHOCHO] = - w[rDD47] - w[rDD50] - w[rDD48] - w[rDD53]
		- w[rDD52] + w[rLPDD59] + w[rDD54];

	cdot[sOCHO] = - w[rDD44] + w[rDD66] + w[rDD427];

	cdot[sC2H2] = - w[rG116] - w[rG117] - w[rLPG120f] + w[rLPG120b]
		- w[rG122] + w[rG133] + w[rLPR111] + w[rLPR007]
		- w[rP001f] + w[rP001b] - w[rP002f] + w[rP002b]
		- w[rP005f] + w[rP005b] - w[rP006f] + w[rP006b]
		- w[rP007f] + w[rP007b] - w[rK012f] + w[rK012b]
		- w[rK100f] + w[rK100b] - w[rK102f] + w[rK102b]
		- w[rP104f] + w[rP104b] - w[rK114f] + w[rK114b]
		- w[rK115f] + w[rK115b] - w[rK205f] + w[rK205b]
		- w[rK401f] + w[rK401b] - w[rK403f] + w[rK403b]
		- w[rP405f] + w[rP405b] - w[rP406f] + w[rP406b]
		- w[rP407f] + w[rP407b] + w[rP425f] - w[rP425b]
		- w[rP501f] + w[rP501b] - w[rP502f] + w[rP502b]
		- w[rP503f] + w[rP503b] - w[rK600f] + w[rK600b]
		- w[rK700f] + w[rK700b] - w[rK701f] + w[rK701b]
		+ w[rCP03f] - w[rCP03b] + w[rCP04] + w[rCP13]
		+ w[rCP15] - w[rCP16f] + w[rCP16b] + 2 * w[rCP31]
		- w[rI03f] + w[rI03b] + w[rI31] + w[rT05f]
		- w[rT05b] + w[rXY12] + w[rN07f] - w[rN07b]
		- w[rOX02f] + w[rOX02b] + w[rOX28] + w[rG130]
		- w[rG115f] + w[rG115b];

	cdot[sC2H3] = - w[rG135] - w[rLPG134] - w[rG133] - w[rLPG131]
		+ w[rG160] + w[rG156] + w[rLPR015] + w[rLPR140]
		+ w[rP013b] - w[rP014f] + w[rP014b] + w[rP015f]
		- w[rP015b] - w[rP105f] + w[rP105b] + w[rP116b]
		+ w[rP117b] - w[rP118f] + w[rP118b] - w[rP119f]
		+ w[rP119b] + w[rP206b] - w[rP207f] + w[rP207b]
		+ w[rP408b] + w[rP409b] - w[rP410f] + w[rP410b]
		- w[rP411f] + w[rP411b] + w[rP504b] + w[rCP10b]
		+ w[rT60] + w[rST06f] - w[rST06b] - w[rG130]
		- w[rP013f] - w[rP116f] - w[rP117f] - w[rP206f]
		- w[rP408f] - w[rP409f] - w[rP504f] - w[rCP10f]
		+ w[rG115f] - w[rG115b];

	cdot[sHCCO] = w[rG116] - w[rG111] + w[rLPG123f] - w[rLPG123b]
		+ w[rLPG127f] - w[rLPG127b] - w[rST00b] + w[rOX62]
		+ w[rOX63] + w[rOX64] + w[rOX65] + w[rOX66]
		+ w[rNO274b] - w[rG109] - w[rR022] - w[rG110]
		- w[rR024] + w[rST00f] - w[rNO274f];

	cdot[sCH2CHO] = w[rLPG134] + w[rLPG131] + w[rLPG157] - w[rLPR030]
		- w[rLPG139] - w[rLPG140] - w[rLPG137] - w[rLPG138]
		- w[rLPG136] - w[rLPG144] + w[rLPR142] - w[rLPDD76f]
		+ w[rLPDD76b];

	cdot[sC2H5O] = w[rG161] + w[rLPDD15] - w[rG168] + w[rR037];

	cdot[sAXC3H5] = - w[rLPR015] - w[rLPR116] - w[rLPR115] - w[rLPR133]
		- w[rLPLPR118] - w[rLPR016] - w[rLPLPR113] - w[rLPR111]
		- w[rLPR007] - w[rLPLPR117] + w[rLPR153] + w[rLPR151]
		+ w[rLPHP43] + w[rLPDD205] + w[rLPHP22] - w[rP008]
		+ w[rCP03f] - w[rCP03b];

	cdot[sC3H6] = w[rLPR116] + w[rLPR115] + w[rLPR016] - w[rR141]
		- w[rR143] - w[rLPR153] - w[rLPR142] - w[rLPR151]
		+ w[rG184] + w[rLPHP53] + w[rHP32] + w[rDD208]
		+ w[rHP11] + w[rDD415] + w[rG183];

	cdot[sCH3COCH2] = w[rLPLPR118] + w[rLPLPR113] + w[rLPLPR117] + w[rLPLPDD495]
		- w[rLPR140] - w[rLPLPDD107];

	cdot[sNXC3H7] = - w[rG184] + w[rLPHP55] + w[rHP29] + w[rLPDD205]
		+ w[rDD208] + w[rLPHP12] + w[rDD323] + w[rDD421]
		+ w[rHP60] - w[rG183] - w[rG163];

	cdot[sPXC4H8] = - w[rLPHP55] - w[rLPHP43] - w[rLPHP54] - w[rLPHP52]
		- w[rLPHP53] - w[rLPDD514] - w[rLPHP56] + w[rLPHP39]
		+ w[rLPHP12] + w[rLPDD286] + w[rLPDD406];

	cdot[sPXC4H9] = w[rLPDD514] - w[rLPHP39] - w[rHP41] + w[rDD214]
		+ w[rHP26] + w[rLPHP22] + w[rHP11] + w[rLPDD286]
		+ w[rDD425];

	cdot[sC5H11] = - w[rHP29] - w[rHP32] + w[rDD199] + w[rHP10]
		+ w[rDD436] + w[rDD416];

	cdot[sC6H12X1XC6H12] = - w[rDD199] - w[rDD212] - w[rLPDD205] + w[rDD288]
		+ w[rDD323] + w[rDD418] + w[rDD239] + w[rDD259];

	cdot[sC6H13] = w[rDD212] - w[rDD214] - w[rDD208] + w[rDD43]
		+ w[rDD287] + 2 * w[rDD210] + w[rDD418] + w[rLPDD293]
		- w[rDD239];

	cdot[sC7H14] = - w[rDD43] - w[rHP26] - w[rLPHP22] - w[rHP15]
		+ w[rDD282] + w[rDD319] + w[rDD416];

	cdot[sC7H15] = w[rHP15] - w[rLPHP12] - w[rHP11] - w[rHP10]
		+ w[rDD324] + w[rDD436] + w[rDD290] - w[rDD263];

	cdot[sC8H17] = - w[rLPDD286] - w[rDD282] - w[rDD288] - w[rDD287]
		+ w[rDD425] + w[rLPDD406] + w[rDD332N];

	cdot[sC9H19] = - w[rDD323] - w[rDD324] - w[rDD319] + w[rDD421]
		+ w[rDD415];

	cdot[sNXC12H26] = - w[rDD210] - w[rDD421] - w[rDD436] - w[rDD425]
		- w[rDD422] - w[rDD434] - w[rDD429] - w[rDD424]
		- w[rDD432] - w[rDD427] - w[rLPDD431];

	cdot[sC12H25] = w[rDD422] + w[rDD434] + w[rDD429] + w[rDD424]
		+ w[rDD432] + w[rDD427] + w[rLPDD431] - w[rDD405f]
		+ w[rDD405b] - w[rDD418] - w[rDD415] - w[rDD416]
		- w[rLPDD406] - w[rDD408f] + w[rDD408b];

	cdot[sC12H24] = w[rDD405f] - w[rDD405b] + w[rDD472] + w[rLPDD473];

	cdot[sC12H25O2] = w[rDD408f] - w[rDD408b] - w[rDD472] - w[rLPDD465f]
		+ w[rLPDD465b];

	cdot[sC12OOHXT6XC12H25O2] = w[rLPDD465f] - w[rLPDD465b] - w[rLPDD473] - w[rLPDD486f]
		+ w[rLPDD486b];

	cdot[sC12OOHO2XT6XC12H25O4] = w[rLPDD486f] - w[rLPDD486b] - w[rLPDD502];

	cdot[sC12KETXT6XC12H24O3] = w[rLPDD502] - w[rLPLPDD495] - w[rLPDD499] - w[rLPDD490];

	cdot[sNXC8H17CHO] = w[rLPLPDD495] - w[rDD332N];

	cdot[sC7H15CO] = w[rLPDD499] + w[rDD294] - w[rDD290] - w[rLPDD293];

	cdot[sC3H7CHO] = w[rLPDD499] - w[rHP60];

	cdot[sNXC7H15CHO] = w[rLPDD490] - w[rDD294];

	cdot[sC7H15O2X1XC7H15O2] = w[rDD263] - w[rDD260f] + w[rDD260b];

	cdot[sC7H14OOH1X3XC7H15O2] = w[rDD260f] - w[rDD260b] - w[rDD259];

	cdot[sC5H4CH2] = - w[rP000f] + w[rP000b] + w[rP001f] - w[rP001b]
		+ w[rP002f] - w[rP002b] - w[rP003f] + w[rP003b]
		- w[rP004f] + w[rP004b] + w[rP008] + w[rP009f]
		- w[rP009b] + w[rCP18] + w[rT47] + w[rT59];

	cdot[sA1XC6H6] = w[rP000f] - w[rP000b] + w[rP003f] - w[rP003b]
		+ w[rP006f] - w[rP006b] + w[rP007f] - w[rP007b]
		+ w[rP010f] - w[rP010b] - w[rP014f] + w[rP014b]
		+ w[rP015f] - w[rP015b] - w[rP301f] + w[rP301b]
		- w[rP417f] + w[rP417b] - w[rP801] + w[rT01f]
		- w[rT01b] + w[rT16f] - w[rT16b] + w[rT28f]
		- w[rT28b] + w[rT31f] - w[rT31b] + w[rE04f]
		- w[rE04b] + w[rST01] - w[rST00b] + w[rXY12]
		- w[rOX00f] + w[rOX00b] - w[rOX03f] + w[rOX03b]
		- w[rOX04f] + w[rOX04b] - w[rOX05f] + w[rOX05b]
		- w[rOX06f] + w[rOX06b] - w[rOX99f] + w[rOX99b]
		- w[rOX07f] + w[rOX07b] - w[rOX08f] + w[rOX08b]
		- w[rOX09f] + w[rOX09b] - w[rOX10f] + w[rOX10b]
		+ w[rOX16f] - w[rOX16b] + w[rST00f];

	cdot[sNXC4H5] = - w[rP001f] + w[rP001b] - w[rP006f] + w[rP006b]
		+ w[rCP11b] + w[rCP24] - w[rCP11f];

	cdot[sIXC4H5] = - w[rP002f] + w[rP002b] - w[rP007f] + w[rP007b];

	cdot[sA1XXC6H5] = w[rP004f] - w[rP004b] + w[rP005f] - w[rP005b]
		+ w[rP011f] - w[rP011b] - w[rK012f] + w[rK012b]
		+ w[rP013b] - w[rP015f] + w[rP015b] - w[rP107f]
		+ w[rP107b] - w[rP301f] + w[rP301b] - 2 * w[rP302f]
		+ 2 * w[rP302b] - w[rP416f] + w[rP416b] - w[rP418f]
		+ w[rP418b] - w[rP800] - w[rP802] + w[rI00b]
		+ w[rT03f] - w[rT03b] + w[rT04f] - w[rT04b]
		- w[rT16f] + w[rT16b] - w[rT28f] + w[rT28b]
		+ w[rT30f] - w[rT30b] + w[rT36] + w[rT37]
		+ w[rT38] + w[rT39] + w[rT40] + w[rT41]
		+ w[rT42] + w[rE03b] + w[rE11f] - w[rE11b]
		+ w[rE37] + w[rXY24] + w[rXY44] + w[rXY50]
		+ w[rOX00f] - w[rOX00b] + w[rOX03f] - w[rOX03b]
		+ w[rOX04f] - w[rOX04b] + w[rOX06f] - w[rOX06b]
		+ w[rOX10f] - w[rOX10b] - w[rOX11f] + w[rOX11b]
		- w[rOX12f] + w[rOX12b] + w[rOX13b] + w[rOX14b]
		+ w[rOX15b] - w[rOX16f] + w[rOX16b] - w[rP013f]
		- w[rI00f] - w[rE03f] - w[rOX13f] - w[rOX14f]
		- w[rOX15f] - w[rOX01f] + w[rOX01b];

	cdot[sNXC4H3] = - w[rP005f] + w[rP005b];

	cdot[sC3H3] = - w[rP008] - 2 * w[rP009f] + 2 * w[rP009b] - 2 * w[rP010f]
		+ 2 * w[rP010b] - 2 * w[rP011f] + 2 * w[rP011b] - w[rCP16f]
		+ w[rCP16b] + w[rI00b] - w[rI23] - w[rT20]
		- w[rXY22] - w[rN19] - w[rI00f];

	cdot[sA1C2H2XC8H7] = w[rK012f] - w[rK012b] - w[rK019f] + w[rK019b]
		- w[rK020f] + w[rK020b] - w[rK021f] + w[rK021b]
		+ w[rP022b] + w[rP026f] - w[rP026b] + w[rP027f]
		- w[rP027b] + w[rP028f] - w[rP028b] - w[rP104f]
		+ w[rP104b] + w[rST03f] - w[rST03b] - w[rST11f]
		+ w[rST11b] - w[rST12f] + w[rST12b] - w[rST13]
		- w[rST14f] + w[rST14b] - w[rP022f];

	cdot[sA1C2H3XC8H8] = - w[rP013b] + w[rP014f] - w[rP014b] - w[rP023f]
		+ w[rP023b] - w[rK024f] + w[rK024b] - w[rP025f]
		+ w[rP025b] - w[rP026f] + w[rP026b] - w[rP027f]
		+ w[rP027b] - w[rP028f] + w[rP028b] + w[rE12f]
		- w[rE12b] - w[rE13b] - w[rE14b] + w[rE15f]
		- w[rE15b] + w[rE16f] - w[rE16b] + w[rE33]
		- w[rST01] - w[rST02f] + w[rST02b] - w[rST03f]
		+ w[rST03b] - w[rST04f] + w[rST04b] - w[rST05f]
		+ w[rST05b] - w[rST06f] + w[rST06b] - w[rST10f]
		+ w[rST10b] + w[rP013f] + w[rE13f] + w[rE14f];

	cdot[sA1C2HXC8H6] = - w[rP016f] + w[rP016b] - w[rK017f] + w[rK017b]
		- w[rP018f] + w[rP018b] + w[rK020f] - w[rK020b]
		+ w[rK021f] - w[rK021b] - w[rP022b] - w[rP105f]
		+ w[rP105b] - w[rP416f] + w[rP416b] + w[rST12f]
		- w[rST12b] + w[rST00b] + w[rP022f] - w[rST00f];

	cdot[sA1C2HYXC8H5] = w[rP016f] - w[rP016b] + w[rK017f] - w[rK017b]
		+ w[rP018f] - w[rP018b] - w[rK100f] + w[rK100b]
		- w[rP106f] + w[rP106b] - w[rP417f] + w[rP417b]
		- w[rP418f] + w[rP418b];

	cdot[sA1C2H3YXC8H7] = w[rK019f] - w[rK019b] + w[rP023f] - w[rP023b]
		+ w[rK024f] - w[rK024b] + w[rP025f] - w[rP025b]
		- w[rK102f] + w[rK102b] + w[rI32] + w[rST10f]
		- w[rST10b];

	cdot[sA2XXC10H7] = w[rK100f] - w[rK100b] + w[rP108f] - w[rP108b]
		+ w[rK109f] - w[rK109b] + w[rK110f] - w[rK110b]
		- w[rK114f] + w[rK114b] + w[rP116b] - w[rP120f]
		+ w[rP120b] - w[rP414f] + w[rP414b] - w[rP801]
		- w[rP802] + w[rN05f] - w[rN05b] + w[rN06f]
		- w[rN06b] + w[rN23f] - w[rN23b] + w[rN28]
		+ w[rN29] + w[rN30] + w[rN31] + w[rN32]
		+ w[rN33] + w[rN34] - w[rOX33f] + w[rOX33b]
		- w[rOX35] + w[rOX37b] + w[rOX39b] + w[rOX53]
		- w[rP116f] - w[rOX37f] - w[rOX39f];

	cdot[sA2XC10H8] = w[rK102f] - w[rK102b] + w[rP104f] - w[rP104b]
		+ w[rP105f] - w[rP105b] + w[rP106f] - w[rP106b]
		+ w[rP107f] - w[rP107b] - w[rP108f] + w[rP108b]
		- w[rK109f] + w[rK109b] - w[rK110f] + w[rK110b]
		- w[rP111f] + w[rP111b] - w[rK112f] + w[rK112b]
		- w[rK113f] + w[rK113b] - w[rP118f] + w[rP118b]
		- w[rP119f] + w[rP119b] - w[rP800] + w[rCP17]
		+ w[rI18] + w[rT20] + w[rN01f] - w[rN01b]
		+ w[rN14] - w[rOX30f] + w[rOX30b] - w[rOX31f]
		+ w[rOX31b] - w[rOX32f] + w[rOX32b] + w[rOX63];

	cdot[sC4H4] = - w[rP107f] + w[rP107b] - w[rP414f] + w[rP414b]
		- w[rP415f] + w[rP415b] - w[rP505f] + w[rP505b]
		+ w[rCP33f] - w[rCP33b];

	cdot[sA2YXC10H7] = w[rP111f] - w[rP111b] + w[rK112f] - w[rK112b]
		+ w[rK113f] - w[rK113b] - w[rK115f] + w[rK115b]
		+ w[rP117b] - w[rP121f] + w[rP121b] - w[rP415f]
		+ w[rP415b] - w[rOX34f] + w[rOX34b] - w[rOX36]
		+ w[rOX38b] + w[rOX40b] - w[rP117f] - w[rOX38f]
		- w[rOX40f];

	cdot[sA2C2H2AXC12H9] = w[rK114f] - w[rK114b] - w[rP116b] + w[rP118f]
		- w[rP118b] + w[rP120f] - w[rP120b] - w[rK122f]
		+ w[rK122b] - w[rK123f] + w[rK123b] + w[rP124b]
		- w[rK200f] + w[rK200b] - w[rP405f] + w[rP405b]
		+ w[rOX51] + w[rP116f] - w[rP124f];

	cdot[sA2C2H2BXC12H9] = w[rK115f] - w[rK115b] - w[rP117b] + w[rP119f]
		- w[rP119b] + w[rP121f] - w[rP121b] - w[rK125f]
		+ w[rK125b] - w[rK126f] + w[rK126b] + w[rP127b]
		- w[rP406f] + w[rP406b] + w[rOX50] + w[rP117f]
		- w[rP127f];

	cdot[sA2C2HAXC12H8] = w[rK122f] - w[rK122b] + w[rK123f] - w[rK123b]
		- w[rP124b] - w[rP128f] + w[rP128b] - w[rK129f]
		+ w[rK129b] - w[rP130f] + w[rP130b] - w[rP201f]
		+ w[rP201b] - w[rP410f] + w[rP410b] + w[rOX60]
		+ w[rP124f];

	cdot[sA2C2HBXC12H8] = w[rK125f] - w[rK125b] + w[rK126f] - w[rK126b]
		- w[rP127b] - w[rP131f] + w[rP131b] - w[rK132f]
		+ w[rK132b] - w[rP133f] + w[rP133b] - w[rP411f]
		+ w[rP411b] + w[rOX61] + w[rP127f];

	cdot[sA2C2HAYXC12H7] = w[rP128f] - w[rP128b] + w[rK129f] - w[rK129b]
		+ w[rP130f] - w[rP130b] - w[rK401f] + w[rK401b]
		+ w[rP408b] - w[rP412f] + w[rP412b] - w[rP408f];

	cdot[sA2C2HBYXC12H7] = w[rP131f] - w[rP131b] + w[rK132f] - w[rK132b]
		+ w[rP133f] - w[rP133b] - w[rK403f] + w[rK403b]
		+ w[rP409b] - w[rP413f] + w[rP413b] - w[rP409f];

	cdot[sA2R5XC12H8] = w[rK200f] - w[rK200b] + w[rP201f] - w[rP201b]
		- w[rP202f] + w[rP202b] - w[rK203f] + w[rK203b]
		- w[rK204f] + w[rK204b] - w[rP207f] + w[rP207b]
		+ w[rI23] - w[rOX63];

	cdot[sA2R5XXC12H7] = w[rP202f] - w[rP202b] + w[rK203f] - w[rK203b]
		+ w[rK204f] - w[rK204b] - w[rK205f] + w[rK205b]
		+ w[rP206b] - w[rP208f] + w[rP208b] + w[rP425f]
		- w[rP425b] - w[rP505f] + w[rP505b] - w[rOX53]
		- w[rP206f];

	cdot[sA2R5C2H2XC14H9] = w[rK205f] - w[rK205b] - w[rP206b] + w[rP207f]
		- w[rP207b] + w[rP208f] - w[rP208b] - w[rK212f]
		+ w[rK212b] - w[rK213f] + w[rK213b] + w[rP214b]
		- w[rP502f] + w[rP502b] + w[rP206f] - w[rP214f];

	cdot[sA2R5C2HXC14H8] = - w[rP209f] + w[rP209b] - w[rK210f] + w[rK210b]
		- w[rP211f] + w[rP211b] + w[rK212f] - w[rK212b]
		+ w[rK213f] - w[rK213b] - w[rP214b] + w[rP214f];

	cdot[sA2R5C2HYXC14H7] = w[rP209f] - w[rP209b] + w[rK210f] - w[rK210b]
		+ w[rP211f] - w[rP211b] - w[rP501f] + w[rP501b];

	cdot[sP2XC12H10] = w[rP301f] - w[rP301b] + w[rP302f] - w[rP302b]
		- w[rP304f] + w[rP304b] - w[rP305f] + w[rP305b]
		- w[rP306f] + w[rP306b];

	cdot[sP2XXC12H9] = w[rP304f] - w[rP304b] + w[rP305f] - w[rP305b]
		+ w[rP306f] - w[rP306b] - w[rP407f] + w[rP407b];

	cdot[sA3XXC14H9] = w[rK401f] - w[rK401b] + w[rK403f] - w[rK403b]
		+ w[rP419f] - w[rP419b] + w[rK420f] - w[rK420b]
		+ w[rK421f] - w[rK421b] - w[rP425f] + w[rP425b]
		- w[rK600f] + w[rK600b] - w[rOX50] + w[rOX52]
		+ w[rOX54];

	cdot[sA3XC14H10] = w[rP405f] - w[rP405b] + w[rP406f] - w[rP406b]
		+ w[rP407f] - w[rP407b] - w[rP408b] - w[rP409b]
		+ w[rP410f] - w[rP410b] + w[rP411f] - w[rP411b]
		+ w[rP412f] - w[rP412b] + w[rP413f] - w[rP413b]
		+ w[rP414f] - w[rP414b] + w[rP415f] - w[rP415b]
		+ w[rP416f] - w[rP416b] + w[rP417f] - w[rP417b]
		+ w[rP418f] - w[rP418b] - w[rP419f] + w[rP419b]
		- w[rK420f] + w[rK420b] - w[rK421f] + w[rK421b]
		- w[rP422f] + w[rP422b] - w[rK423f] + w[rK423b]
		- w[rK424f] + w[rK424b] + w[rI17] + w[rN19]
		- w[rOX60] - w[rOX61] + w[rOX62] + w[rOX64]
		+ w[rOX66] + w[rP408f] + w[rP409f];

	cdot[sA3YXC14H9] = w[rP422f] - w[rP422b] + w[rK423f] - w[rK423b]
		+ w[rK424f] - w[rK424b] - w[rP503f] + w[rP503b]
		+ w[rP504b] - w[rOX51] - w[rP504f];

	cdot[sA3R5XXC16H9] = w[rP501f] - w[rP501b] + w[rP506f] - w[rP506b]
		+ w[rK507f] - w[rK507b] + w[rK508f] - w[rK508b]
		- w[rK701f] + w[rK701b] - w[rOX54];

	cdot[sA3R5XC16H10] = w[rP502f] - w[rP502b] + w[rP503f] - w[rP503b]
		- w[rP504b] + w[rP505f] - w[rP505b] - w[rP506f]
		+ w[rP506b] - w[rK507f] + w[rK507b] - w[rK508f]
		+ w[rK508b] - w[rOX64] + w[rP504f];

	cdot[sA4XC16H10] = w[rK600f] - w[rK600b] - w[rP601f] + w[rP601b]
		- w[rK602f] + w[rK602b] - w[rK603f] + w[rK603b]
		- w[rOX62] + w[rOX65];

	cdot[sA4XXC16H9] = w[rP601f] - w[rP601b] + w[rK602f] - w[rK602b]
		+ w[rK603f] - w[rK603b] - w[rK700f] + w[rK700b]
		- w[rOX52];

	cdot[sA4R5XC18H10] = w[rK700f] - w[rK700b] + w[rK701f] - w[rK701b]
		- w[rOX65];

	cdot[sFLTNXC16H10] = w[rP800] + w[rP801] + w[rP802] - w[rOX66];

	cdot[sC5H6] = - w[rCP01f] + w[rCP01b] - w[rCP02f] + w[rCP02b]
		- w[rCP03f] + w[rCP03b] - w[rCP04] - w[rCP05f]
		+ w[rCP05b] - w[rCP06f] + w[rCP06b] - w[rCP07f]
		+ w[rCP07b] - w[rCP08f] + w[rCP08b] - w[rCP09f]
		+ w[rCP09b] + w[rCP10b] + w[rCP11b] - w[rCP12f]
		+ w[rCP12b] - w[rCP13] - w[rCP14] - w[rCP15]
		- w[rCP34] + w[rOX09f] - w[rOX09b] + w[rOX17f]
		- w[rOX17b] - w[rCP10f] - w[rCP11f];

	cdot[sC5H5] = w[rCP01f] - w[rCP01b] + w[rCP02f] - w[rCP02b]
		+ w[rCP05f] - w[rCP05b] + w[rCP06f] - w[rCP06b]
		+ w[rCP07f] - w[rCP07b] + w[rCP08f] - w[rCP08b]
		+ w[rCP09f] - w[rCP09b] - w[rCP10b] - w[rCP11b]
		+ w[rCP16f] - w[rCP16b] - 2 * w[rCP17] - w[rCP18]
		+ w[rCP19b] - w[rCP20f] + w[rCP20b] + w[rCP21b]
		- w[rCP34] - w[rI17] + w[rT05f] - w[rT05b]
		+ w[rOX23f] - w[rOX23b] + w[rCP10f] + w[rCP11f]
		- w[rCP19f] - w[rCP21f] - w[rCP22];

	cdot[sSXC3H5] = w[rCP04];

	cdot[sCH4] = w[rCP09f] - w[rCP09b] + w[rI09f] - w[rI09b]
		+ w[rT14f] - w[rT14b] + w[rT27f] - w[rT27b]
		+ w[rT42] + w[rT54f] - w[rT54b] + w[rE10f]
		- w[rE10b] + w[rE15f] - w[rE15b] + w[rST03f]
		- w[rST03b] + w[rXY07f] - w[rXY07b] + w[rXY31f]
		- w[rXY31b] + w[rXY38] + w[rN12f] - w[rN12b]
		+ w[rN34] - w[rOX16f] + w[rOX16b] + w[rOX21f]
		- w[rOX21b] + w[rOX45f] - w[rOX45b] - w[rNO211b]
		+ w[rNO211f];

	cdot[sC4H6] = - w[rCP11b] + w[rCP14] + w[rCP11f];

	cdot[sTXC5H5O] = w[rCP12f] - w[rCP12b] - w[rCP24] - w[rCP30]
		+ w[rCP32f] - w[rCP32b] + w[rOX27f] - w[rOX27b]
		- w[rCP29];

	cdot[sC5H4O] = - w[rCP19b] + w[rCP20f] - w[rCP20b] + w[rCP23f]
		- w[rCP23b] + w[rCP28] + w[rCP30] - w[rCP31]
		- w[rCP32f] + w[rCP32b] - w[rCP33f] + w[rCP33b]
		+ w[rOX26f] - w[rOX26b] + w[rCP19f] + w[rCP25]
		+ w[rCP26] + w[rCP27] + w[rCP29];

	cdot[sSXC5H5O] = - w[rCP21b] - w[rCP23f] + w[rCP23b] - w[rCP28]
		+ w[rCP21f] + w[rCP22] - w[rCP25] - w[rCP26]
		- w[rCP27];

	cdot[sC9H8] = w[rCP34] - w[rI00b] - w[rI01f] + w[rI01b]
		- w[rI02f] + w[rI02b] + w[rI03f] - w[rI03b]
		- w[rI05f] + w[rI05b] - w[rI06f] + w[rI06b]
		- w[rI07f] + w[rI07b] - w[rI08f] + w[rI08b]
		- w[rI09f] + w[rI09b] - w[rI12] - w[rI15]
		+ w[rOX41f] - w[rOX41b] + w[rI00f];

	cdot[sC9H7] = w[rI01f] - w[rI01b] + w[rI02f] - w[rI02b]
		+ w[rI05f] - w[rI05b] + w[rI06f] - w[rI06b]
		+ w[rI07f] - w[rI07b] + w[rI08f] - w[rI08b]
		+ w[rI09f] - w[rI09b] - w[rI17] - w[rI18]
		+ w[rI19b] - w[rI20f] + w[rI20b] - w[rI23]
		+ w[rN07f] - w[rN07b] + w[rN15] + w[rOX46f]
		- w[rOX46b] - w[rI19f] - w[rI21] - w[rI22];

	cdot[sA1CH2XC7H7] = - w[rI03f] + w[rI03b] + w[rT02f] - w[rT02b]
		- w[rT04f] + w[rT04b] - w[rT05f] + w[rT05b]
		+ w[rT06f] - w[rT06b] + w[rT07f] - w[rT07b]
		+ w[rT08f] - w[rT08b] + w[rT11f] - w[rT11b]
		+ w[rT14f] - w[rT14b] + w[rT15f] - w[rT15b]
		+ w[rT16f] - w[rT16b] + w[rT17b] + w[rT18b]
		- w[rT19f] + w[rT19b] - w[rT20] - w[rT21f]
		+ w[rT21b] - w[rT22f] + w[rT22b] + w[rE02b]
		+ w[rE17f] - w[rE17b] + w[rE19] + w[rE31]
		+ w[rE34] + w[rST02f] - w[rST02b] + w[rST05f]
		- w[rST05b] + w[rST11f] - w[rST11b] + w[rST13]
		- w[rT17f] - w[rT18f] - w[rE02f];

	cdot[sC9H6O] = w[rI12] - w[rI19b] + w[rI20f] - w[rI20b]
		- w[rI31] - w[rI32] + w[rOX35] + w[rOX36]
		+ w[rOX48] + w[rI19f] + w[rI21] + w[rI22]
		+ w[rOX47];

	cdot[sOXC6H4] = w[rI15] + w[rI31] + w[rOX02f] - w[rOX02b]
		+ w[rOX01f] - w[rOX01b];

	cdot[sA1CH3XC7H8] = - w[rT01f] + w[rT01b] - w[rT02f] + w[rT02b]
		- w[rT03f] + w[rT03b] - w[rT06f] + w[rT06b]
		- w[rT07f] + w[rT07b] - w[rT08f] + w[rT08b]
		- w[rT09f] + w[rT09b] - w[rT10f] + w[rT10b]
		- w[rT11f] + w[rT11b] - w[rT12f] + w[rT12b]
		- w[rT13f] + w[rT13b] - w[rT14f] + w[rT14b]
		- w[rT15f] + w[rT15b] - w[rT16f] + w[rT16b]
		- w[rT50f] + w[rT50b] - w[rT51f] + w[rT51b]
		- w[rT52f] + w[rT52b] - w[rT53f] + w[rT53b]
		- w[rT54f] + w[rT54b] + w[rXY09f] - w[rXY09b]
		+ w[rXY11] + w[rXY16f] - w[rXY16b] + w[rXY203]
		+ w[rXY39f] - w[rXY39b];

	cdot[sA1OHXC6H6O] = w[rT09f] - w[rT09b] + w[rE05f] - w[rE05b]
		+ w[rST06f] - w[rST06b] + w[rOX05f] - w[rOX05b]
		+ w[rOX08f] - w[rOX08b] - w[rOX17f] + w[rOX17b]
		- w[rOX18f] + w[rOX18b] - w[rOX19f] + w[rOX19b]
		- w[rOX20f] + w[rOX20b] - w[rOX21f] + w[rOX21b]
		- w[rOX22f] + w[rOX22b];

	cdot[sHOA1CH3XC7H8O] = w[rT10f] - w[rT10b] + w[rT12f] - w[rT12b]
		- w[rT43f] + w[rT43b] - w[rT44f] + w[rT44b]
		- w[rT45f] + w[rT45b] - w[rT46f] + w[rT46b]
		+ w[rXY10f] - w[rXY10b] + w[rXY41f] - w[rXY41b];

	cdot[sOA1CH3XC7H7O] = w[rT13f] - w[rT13b] + w[rT43f] - w[rT43b]
		+ w[rT44f] - w[rT44b] + w[rT45f] - w[rT45b]
		+ w[rT46f] - w[rT46b] - w[rT47] - w[rT55b]
		- w[rT56b] - w[rT57b] + w[rT58f] - w[rT58b]
		+ w[rXY19f] - w[rXY19b] + w[rT55f] + w[rT56f]
		+ w[rT57f];

	cdot[sA1CH2OXC7H7O] = - w[rT17b] + w[rT19f] - w[rT19b] + w[rT24f]
		- w[rT24b] + w[rT25f] - w[rT25b] + w[rT26f]
		- w[rT26b] + w[rT27f] - w[rT27b] - w[rT29f]
		+ w[rT29b] - w[rT30f] + w[rT30b] - w[rT31f]
		+ w[rT31b] + w[rT32b] + w[rT33b] + w[rT34b]
		- w[rT35f] + w[rT35b] + w[rT17f] - w[rT32f]
		- w[rT33f] - w[rT34f] - w[rT23f] + w[rT23b];

	cdot[sA1CH2OHXC7H8O] = - w[rT18b] - w[rT24f] + w[rT24b] - w[rT25f]
		+ w[rT25b] - w[rT26f] + w[rT26b] - w[rT27f]
		+ w[rT27b] + w[rT18f] + w[rT23f] - w[rT23b];

	cdot[sA1CHOXC7H6O] = w[rT21f] - w[rT21b] + w[rT29f] - w[rT29b]
		- w[rT32b] - w[rT33b] - w[rT34b] + w[rT35f]
		- w[rT35b] - w[rT36] - w[rT37] - w[rT38]
		- w[rT39] - w[rT40] - w[rT41] - w[rT42]
		+ w[rE18f] - w[rE18b] + w[rST04f] - w[rST04b]
		+ w[rST14f] - w[rST14b] + w[rXY40f] - w[rXY40b]
		+ w[rXY45f] - w[rXY45b] + w[rXY51] + w[rXY52]
		+ w[rXY53] + w[rXY54] + w[rXY55] + w[rXY56]
		+ w[rXY57f] - w[rXY57b] + w[rT32f] + w[rT33f]
		+ w[rT34f];

	cdot[sA1OXC6H5O] = w[rT22f] - w[rT22b] + w[rXY42] + w[rXY58]
		+ w[rOX99f] - w[rOX99b] + w[rOX07f] - w[rOX07b]
		+ w[rOX11f] - w[rOX11b] - w[rOX13b] - w[rOX14b]
		- w[rOX15b] + w[rOX18f] - w[rOX18b] + w[rOX19f]
		- w[rOX19b] + w[rOX20f] - w[rOX20b] + w[rOX21f]
		- w[rOX21b] + w[rOX22f] - w[rOX22b] - w[rOX23f]
		+ w[rOX23b] - w[rOX24f] + w[rOX24b] - w[rOX25f]
		+ w[rOX25b] + w[rOX13f] + w[rOX14f] + w[rOX15f];

	cdot[sA1CH3YXC7H7] = w[rT50f] - w[rT50b] + w[rT51f] - w[rT51b]
		+ w[rT52f] - w[rT52b] + w[rT53f] - w[rT53b]
		+ w[rT54f] - w[rT54b] + w[rT55b] + w[rT56b]
		+ w[rT57b] - w[rT58f] + w[rT58b] - w[rT59]
		- w[rT60] + w[rXY01f] - w[rXY01b] + w[rXY13f]
		- w[rXY13b] + w[rXY15f] - w[rXY15b] + w[rXY202]
		+ w[rXY25] + w[rXY33] + w[rXY34] + w[rXY35]
		+ w[rXY36] + w[rXY37] + w[rXY38] - w[rT55f]
		- w[rT56f] - w[rT57f];

	cdot[sPXC3H4] = w[rT60];

	cdot[sA1C2H4XC8H9] = w[rE01b] + w[rE06f] - w[rE06b] + w[rE07f]
		- w[rE07b] + w[rE08f] - w[rE08b] + w[rE09f]
		- w[rE09b] + w[rE10f] - w[rE10b] - w[rE11f]
		+ w[rE11b] - w[rE12f] + w[rE12b] + w[rE13b]
		+ w[rE14b] - w[rE15f] + w[rE15b] - w[rE16f]
		+ w[rE16b] - w[rE17f] + w[rE17b] - w[rE18f]
		+ w[rE18b] - w[rE19] - w[rE30] - w[rE31]
		+ w[rE32] - w[rE01f] - w[rE13f] - w[rE14f];

	cdot[sA1C2H5XC8H10] = - w[rE01b] - w[rE02b] - w[rE03b] - w[rE04f]
		+ w[rE04b] - w[rE05f] + w[rE05b] - w[rE06f]
		+ w[rE06b] - w[rE07f] + w[rE07b] - w[rE08f]
		+ w[rE08b] - w[rE09f] + w[rE09b] - w[rE10f]
		+ w[rE10b] + w[rE01f] + w[rE02f] + w[rE03f];

	cdot[sC8H9O2] = w[rE30] - w[rE32] - w[rE33] - w[rE34]
		- w[rE35];

	cdot[sC8H8OOH] = w[rE35] - w[rE36];

	cdot[sOC8H7OOH] = - w[rE37] + w[rE36];

	cdot[sH2C2] = w[rST01];

	cdot[sA1CH3CH3XC8H10] = - w[rXY00f] + w[rXY00b] - w[rXY01f] + w[rXY01b]
		- w[rXY02f] + w[rXY02b] - w[rXY03f] + w[rXY03b]
		- w[rXY04f] + w[rXY04b] - w[rXY05f] + w[rXY05b]
		- w[rXY06f] + w[rXY06b] - w[rXY07f] + w[rXY07b]
		- w[rXY09f] + w[rXY09b] - w[rXY10f] + w[rXY10b]
		- w[rXY11];

	cdot[sA1CH3CH2XC8H9] = w[rXY00f] - w[rXY00b] + w[rXY02f] - w[rXY02b]
		+ w[rXY03f] - w[rXY03b] + w[rXY04f] - w[rXY04b]
		+ w[rXY05f] - w[rXY05b] + w[rXY06f] - w[rXY06b]
		+ w[rXY07f] - w[rXY07b] - w[rXY12] - w[rXY13f]
		+ w[rXY13b] - w[rXY14f] + w[rXY14b] - w[rXY15f]
		+ w[rXY15b] - w[rXY16f] + w[rXY16b] - w[rXY18f]
		+ w[rXY18b] - w[rXY19f] + w[rXY19b] - w[rXY201]
		- w[rXY202] - w[rXY203] - w[rXY22] - w[rXY17];

	cdot[sA1CH3CHOXC8H8O] = w[rXY14f] - w[rXY14b] + w[rXY18f] - w[rXY18b]
		+ w[rXY201] - w[rXY23f] + w[rXY23b] - w[rXY24]
		- w[rXY25] - w[rXY26f] + w[rXY26b] - w[rXY27f]
		+ w[rXY27b] - w[rXY28f] + w[rXY28b] - w[rXY29f]
		+ w[rXY29b] - w[rXY30f] + w[rXY30b] - w[rXY31f]
		+ w[rXY31b] - w[rXY33] - w[rXY34] - w[rXY35]
		- w[rXY36] - w[rXY37] - w[rXY38] - w[rXY39f]
		+ w[rXY39b] - w[rXY40f] + w[rXY40b] - w[rXY41f]
		+ w[rXY41b] - w[rXY42] + w[rXY17];

	cdot[sA2CH3XC11H10] = w[rXY22] - w[rN01f] + w[rN01b] - w[rN02f]
		+ w[rN02b] - w[rN04f] + w[rN04b] - w[rN05f]
		+ w[rN05b] - w[rN08f] + w[rN08b] - w[rN09f]
		+ w[rN09b] - w[rN10f] + w[rN10b] - w[rN11f]
		+ w[rN11b] - w[rN12f] + w[rN12b] - w[rN13f]
		+ w[rN13b] - w[rN14] - w[rN15];

	cdot[sA1CHOCH2XC8H7O] = w[rXY23f] - w[rXY23b] + w[rXY26f] - w[rXY26b]
		+ w[rXY27f] - w[rXY27b] + w[rXY28f] - w[rXY28b]
		+ w[rXY29f] - w[rXY29b] + w[rXY30f] - w[rXY30b]
		+ w[rXY31f] - w[rXY31b] - w[rXY43f] + w[rXY43b]
		- w[rXY44] - w[rXY45f] + w[rXY45b] - w[rXY47f]
		+ w[rXY47b] - w[rXY48] - w[rXY46];

	cdot[sA1CHOCHOXC8H6O2] = w[rXY43f] - w[rXY43b] + w[rXY47f] - w[rXY47b]
		+ w[rXY48] - w[rXY50] - w[rXY51] - w[rXY52]
		- w[rXY53] - w[rXY54] - w[rXY55] - w[rXY56]
		- w[rXY57f] + w[rXY57b] - w[rXY58] + w[rXY46];

	cdot[sA2OHXC10H8O] = w[rN02f] - w[rN02b] + w[rOX31f] - w[rOX31b]
		+ w[rOX32f] - w[rOX32b] - w[rOX41f] + w[rOX41b]
		- w[rOX42f] + w[rOX42b] - w[rOX43f] + w[rOX43b]
		- w[rOX44f] + w[rOX44b] - w[rOX45f] + w[rOX45b];

	cdot[sA2CH2XC11H9] = w[rN04f] - w[rN04b] - w[rN06f] + w[rN06b]
		- w[rN07f] + w[rN07b] + w[rN08f] - w[rN08b]
		+ w[rN09f] - w[rN09b] + w[rN10f] - w[rN10b]
		+ w[rN11f] - w[rN11b] + w[rN12f] - w[rN12b]
		+ w[rN13f] - w[rN13b] + w[rN16b] - w[rN18f]
		+ w[rN18b] - w[rN19] - w[rN20f] + w[rN20b]
		- w[rN21f] + w[rN21b] - w[rN16f] - w[rN17];

	cdot[sA2CH2OXC11H9O] = - w[rN16b] + w[rN18f] - w[rN18b] - w[rN22f]
		+ w[rN22b] - w[rN23f] + w[rN23b] + w[rN24b]
		+ w[rN25b] + w[rN26b] - w[rN27f] + w[rN27b]
		+ w[rN16f] + w[rN17] - w[rN24f] - w[rN25f]
		- w[rN26f];

	cdot[sA2CHOXC11H8O] = w[rN20f] - w[rN20b] + w[rN22f] - w[rN22b]
		- w[rN24b] - w[rN25b] - w[rN26b] + w[rN27f]
		- w[rN27b] - w[rN28] - w[rN29] - w[rN30]
		- w[rN31] - w[rN32] - w[rN33] - w[rN34]
		+ w[rN24f] + w[rN25f] + w[rN26f];

	cdot[sA2OXC10H7O] = w[rN21f] - w[rN21b] + w[rOX30f] - w[rOX30b]
		+ w[rOX33f] - w[rOX33b] + w[rOX34f] - w[rOX34b]
		- w[rOX37b] - w[rOX38b] - w[rOX39b] - w[rOX40b]
		+ w[rOX42f] - w[rOX42b] + w[rOX43f] - w[rOX43b]
		+ w[rOX44f] - w[rOX44b] + w[rOX45f] - w[rOX45b]
		- w[rOX46f] + w[rOX46b] - w[rOX48] + w[rOX37f]
		+ w[rOX38f] + w[rOX39f] + w[rOX40f] - w[rOX47];

	cdot[sC4H2] = - w[rOX02f] + w[rOX02b];

	cdot[sOC6H4O] = w[rOX12f] - w[rOX12b] + w[rOX24f] - w[rOX24b]
		+ w[rOX25f] - w[rOX25b] - w[rOX26f] + w[rOX26b]
		- w[rOX27f] + w[rOX27b] - w[rOX28];

	cdot[sN] = - w[rNO178f] + w[rNO178b] - w[rNO179f] + w[rNO179b]
		- w[rNO180f] + w[rNO180b] + w[rNO191f] - w[rNO191b]
		+ w[rNO193f] - w[rNO193b] + w[rNO196b] - w[rNO217b]
		+ w[rNO225b] + w[rNO227f] - w[rNO227b] - w[rNO238f]
		+ w[rNO238b] + w[rNO239f] - w[rNO239b] + w[rNO240f]
		- w[rNO240b] - w[rNO245b] - w[rNO248b] - w[rNO275f]
		+ w[rNO275b] - w[rNO276f] + w[rNO276b] - w[rNO283f]
		+ w[rNO283b] - w[rNO196f] + w[rNO217f] - w[rNO225f]
		+ w[rNO245f] + w[rNO248f];

	cdot[sNO] = - w[rNO178f] + w[rNO178b] + w[rNO179f] - w[rNO179b]
		+ w[rNO180f] - w[rNO180b] + 2 * w[rNO182f] - 2 * w[rNO182b]
		- w[rNO186f] + w[rNO186b] - w[rNO187f] + w[rNO187b]
		+ w[rNO188f] - w[rNO188b] + w[rNO189f] - w[rNO189b]
		- w[rNO190b] + w[rNO195f] - w[rNO195b] - w[rNO198f]
		+ w[rNO198b] - w[rNO199f] + w[rNO199b] - w[rNO208b]
		- w[rNO212f] + w[rNO212b] - w[rNO213b] + w[rNO214f]
		- w[rNO214b] + w[rNO215f] - w[rNO215b] + w[rNO216f]
		- w[rNO216b] - w[rNO222b] - w[rNO224b] + w[rNO226f]
		- w[rNO226b] - w[rNO228f] + w[rNO228b] - w[rNO229f]
		+ w[rNO229b] + w[rNO244b] + w[rNO245b] + w[rNO246b]
		+ w[rNO247b] + w[rNO248b] - w[rNO249f] + w[rNO249b]
		- w[rNO250f] + w[rNO250b] - w[rNO251f] + w[rNO251b]
		- w[rNO252f] + w[rNO252b] - w[rNO253f] + w[rNO253b]
		- w[rNO254f] + w[rNO254b] - w[rNO255f] + w[rNO255b]
		- w[rNO256f] + w[rNO256b] - w[rNO258b] + w[rNO274b]
		+ w[rNO281f] - w[rNO281b] + w[rNO283f] - w[rNO283b]
		+ w[rNO190f] + w[rNO208f] + w[rNO213f] + w[rNO222f]
		+ w[rNO224f] - w[rNO244f] - w[rNO245f] - w[rNO246f]
		- w[rNO247f] - w[rNO248f] + w[rNO258f] - w[rNO274f];

	cdot[sN2O] = - w[rNO181f] + w[rNO181b] - w[rNO182f] + w[rNO182b]
		- w[rNO183f] + w[rNO183b] - w[rNO184f] + w[rNO184b]
		+ w[rNO199f] - w[rNO199b] + w[rNO228f] - w[rNO228b]
		+ w[rNO282f] - w[rNO282b] - w[rNO185f] + w[rNO185b];

	cdot[sNO2] = w[rNO186f] - w[rNO186b] + w[rNO187f] - w[rNO187b]
		- w[rNO188f] + w[rNO188b] - w[rNO189f] + w[rNO189b]
		- w[rNO281f] + w[rNO281b] - w[rNO282f] + w[rNO282b];

	cdot[sNH] = w[rNO190b] - w[rNO191f] + w[rNO191b] + w[rNO192b]
		- w[rNO193f] + w[rNO193b] - w[rNO194f] + w[rNO194b]
		- w[rNO195f] + w[rNO195b] + w[rNO196b] - w[rNO197f]
		+ w[rNO197b] - w[rNO198f] + w[rNO198b] - w[rNO199f]
		+ w[rNO199b] - w[rNO200b] + w[rNO202f] - w[rNO202b]
		+ w[rNO203f] - w[rNO203b] - w[rNO208b] - w[rNO223b]
		+ w[rNO232f] - w[rNO232b] + w[rNO242f] - w[rNO242b]
		+ w[rNO243f] - w[rNO243b] + w[rNO262f] - w[rNO262b]
		+ w[rNO269f] - w[rNO269b] - w[rNO280f] + w[rNO280b]
		- w[rNO190f] - w[rNO192f] - w[rNO196f] + w[rNO200f]
		+ w[rNO208f] + w[rNO223f];

	cdot[sHNO] = - w[rNO192b] + w[rNO194f] - w[rNO194b] + w[rNO197f]
		- w[rNO197b] - w[rNO201b] + w[rNO212f] - w[rNO212b]
		+ w[rNO213b] - w[rNO214f] + w[rNO214b] - w[rNO215f]
		+ w[rNO215b] - w[rNO216f] + w[rNO216b] + w[rNO263f]
		- w[rNO263b] + w[rNO280f] - w[rNO280b] + w[rNO192f]
		+ w[rNO201f] - w[rNO213f];

	cdot[sNH2] = w[rNO200b] + w[rNO201b] - w[rNO202f] + w[rNO202b]
		- w[rNO203f] + w[rNO203b] + w[rNO236f] - w[rNO236b]
		+ w[rNO265f] - w[rNO265b] + w[rNO268f] - w[rNO268b]
		+ w[rNO272f] - w[rNO272b] + w[rNO277f] - w[rNO277b]
		+ w[rNO278f] - w[rNO278b] + w[rNO279f] - w[rNO279b]
		- w[rNO200f] - w[rNO201f];

	cdot[sNNH] = w[rNO204b] - w[rNO205f] + w[rNO205b] + w[rNO206b]
		+ w[rNO207b] + w[rNO208b] + w[rNO209b] + w[rNO210b]
		+ w[rNO211b] - w[rNO204f] - w[rNO206f] - w[rNO207f]
		- w[rNO208f] - w[rNO209f] - w[rNO210f] - w[rNO211f];

	cdot[sCN] = w[rNO217b] + w[rNO218b] - w[rNO219f] + w[rNO219b]
		- w[rNO220f] + w[rNO220b] - w[rNO221f] + w[rNO221b]
		+ w[rNO230f] - w[rNO230b] + w[rNO233f] - w[rNO233b]
		+ w[rNO239f] - w[rNO239b] - w[rNO244b] - w[rNO281f]
		+ w[rNO281b] - w[rNO217f] - w[rNO218f] + w[rNO244f];

	cdot[sNCO] = - w[rNO218b] + w[rNO220f] - w[rNO220b] + w[rNO222b]
		+ w[rNO223b] + w[rNO224b] + w[rNO225b] - w[rNO226f]
		+ w[rNO226b] - w[rNO227f] + w[rNO227b] - w[rNO228f]
		+ w[rNO228b] - w[rNO229f] + w[rNO229b] + w[rNO231f]
		- w[rNO231b] - w[rNO247b] + w[rNO264f] - w[rNO264b]
		+ w[rNO266f] - w[rNO266b] + w[rNO267f] - w[rNO267b]
		+ w[rNO281f] - w[rNO281b] - w[rNO282f] + w[rNO282b]
		+ w[rNO218f] - w[rNO222f] - w[rNO223f] - w[rNO224f]
		- w[rNO225f] + w[rNO247f];

	cdot[sHCN] = w[rNO219f] - w[rNO219b] + w[rNO221f] - w[rNO221b]
		- w[rNO230f] + w[rNO230b] - w[rNO231f] + w[rNO231b]
		- w[rNO232f] + w[rNO232b] - w[rNO233f] + w[rNO233b]
		- w[rNO234f] + w[rNO234b] - w[rNO235] - w[rNO236f]
		+ w[rNO236b] + w[rNO240f] - w[rNO240b] + w[rNO242f]
		- w[rNO242b] + w[rNO243f] - w[rNO243b] - w[rNO246b]
		+ w[rNO250f] - w[rNO250b] + w[rNO253f] - w[rNO253b]
		+ w[rNO255f] - w[rNO255b] - w[rNO258b] + w[rNO271f]
		- w[rNO271b] + w[rNO276f] - w[rNO276b] + w[rNO246f]
		+ w[rNO258f] - w[rNO237f] + w[rNO237b];

	cdot[sHOCN] = w[rNO234f] - w[rNO234b] - w[rNO273];

	cdot[sHNCO] = w[rNO235] + w[rNO249f] - w[rNO249b] + w[rNO252f]
		- w[rNO252b] - w[rNO262f] + w[rNO262b] - w[rNO263f]
		+ w[rNO263b] - w[rNO264f] + w[rNO264b] - w[rNO265f]
		+ w[rNO265b] - w[rNO266f] + w[rNO266b] - w[rNO267f]
		+ w[rNO267b] - w[rNO268f] + w[rNO268b] - w[rNO269f]
		+ w[rNO269b] + w[rNO270f] - w[rNO270b] + w[rNO273];

	cdot[sH2CN] = - w[rNO238f] + w[rNO238b] + w[rNO256f] - w[rNO256b]
		+ w[rNO275f] - w[rNO275b] + w[rNO237f] - w[rNO237b];

	cdot[sC] = - w[rNO239f] + w[rNO239b] + w[rNO244b] + w[rNO245b]
		- w[rNO244f] - w[rNO245f];

	cdot[sHCNN] = w[rNO257b] + w[rNO258b] + w[rNO259b] + w[rNO260b]
		+ w[rNO261b] - w[rNO257f] - w[rNO258f] - w[rNO259f]
		- w[rNO260f] - w[rNO261f] + w[rNO241f] - w[rNO241b];

	cdot[sHCNO] = w[rNO251f] - w[rNO251b] + w[rNO254f] - w[rNO254b]
		- w[rNO270f] + w[rNO270b] - w[rNO271f] + w[rNO271b]
		- w[rNO272f] + w[rNO272b] - w[rNO274b] + w[rNO274f];

	cdot[sNH3] = - w[rNO277f] + w[rNO277b] - w[rNO278f] + w[rNO278b]
		- w[rNO279f] + w[rNO279b];

}

double GetLindRateCoeff( double temp, double pressure
				, double k0, double kInf
				, double fc, double conc )
{
	const double	R = 8314.34;   /* [J / kmole K] */
	double		Ntmp;
	double		kl;
	double		f;
	double		cCoeff, dCoeff, log10kNull;
	int		iTroe = 1;

	if ( conc <= 0.0 ) {
		/*conc = pressure / ( R * temp );*/
		return 0.0;
	}

	Ntmp = 0.75 - 1.27 * log10(MAX_C(fc,1.0e-60));
	cCoeff = - 0.4 - 0.67 * log10(MAX_C(fc,1.0e-60));
	dCoeff = 0.14;

	k0 *= conc / MAX_C(kInf, 1.0e-60);
	log10kNull = log10(MAX_C(k0,1.0e-60));

	f = ( log10kNull + cCoeff ) / ( Ntmp - dCoeff * ( log10kNull + cCoeff ) );
	f = pow( MAX_C(fc,1.0e-60), 1.0 / ( f * f + 1.0 ) );

	kInf *= f * k0 / ( 1.0 + k0 );

	if (isnan(kInf)) {
		fprintf(stdout,"kInf is not a number\n");
		fprintf(stdout,"temp: %15.6e\n",temp);
		fprintf(stdout,"pressure: %15.6e\n",pressure);
		fprintf(stdout,"k0: %15.6e\n",k0);
		fprintf(stdout,"kInf: %15.6e\n",kInf);
		fprintf(stdout,"fc: %15.6e\n",fc);
		fprintf(stdout,"conc: %15.6e\n",conc);
		fprintf(stdout,"f: %15.6e\n",f);
		fprintf(stdout,"Ntmp: %15.6e\n",Ntmp);
		fprintf(stdout,"cCoeff: %15.6e\n",cCoeff);
		fprintf(stdout,"log10kNull: %15.6e\n",log10kNull);
	}
	if (isinf(kInf)) {
		fprintf(stdout,"kInf is not a infinity\n");
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
	W[sH2] =  2.01600000e+00;
	W[sH] =  1.00800000e+00;
	W[sOH] =  1.70080000e+01;
	W[sH2O] =  1.80160000e+01;
	W[sO2] =  3.20000000e+01;
	W[sHO2] =  3.30080000e+01;
	W[sH2O2] =  3.40160000e+01;
	W[sCH] =  1.30180000e+01;
	W[sTXCH2] =  1.40260000e+01;
	W[sCH2O] =  3.00260000e+01;
	W[sHCO] =  2.90180000e+01;
	W[sCH3] =  1.50340000e+01;
	W[sCO2] =  4.40100000e+01;
	W[sCO] =  2.80100000e+01;
	W[sSXCH2] =  1.40260000e+01;
	W[sC2H4] =  2.80520000e+01;
	W[sCH3O] =  3.10340000e+01;
	W[sCH3O2] =  4.70340000e+01;
	W[sC2H5] =  2.90600000e+01;
	W[sCH2CO] =  4.20360000e+01;
	W[sO2CHO] =  6.10180000e+01;
	W[sHOCHO] =  4.60260000e+01;
	W[sOCHO] =  4.50180000e+01;
	W[sC2H2] =  2.60360000e+01;
	W[sC2H3] =  2.70440000e+01;
	W[sHCCO] =  4.10280000e+01;
	W[sCH2CHO] =  4.30440000e+01;
	W[sC2H5O] =  4.50600000e+01;
	W[sAXC3H5] =  4.10700000e+01;
	W[sC3H6] =  4.20780000e+01;
	W[sCH3COCH2] =  5.70700000e+01;
	W[sNXC3H7] =  4.30860000e+01;
	W[sPXC4H8] =  5.61040000e+01;
	W[sPXC4H9] =  5.71120000e+01;
	W[sC5H11] =  7.11380000e+01;
	W[sC6H12X1XC6H12] =  8.41560000e+01;
	W[sC6H13] =  8.51640000e+01;
	W[sC7H14] =  9.81820000e+01;
	W[sC7H15] =  9.91900000e+01;
	W[sC8H17] =  1.13216000e+02;
	W[sC9H19] =  1.27242000e+02;
	W[sNXC12H26] =  1.70328000e+02;
	W[sC12H25] =  1.69320000e+02;
	W[sC12H24] =  1.68312000e+02;
	W[sC12H25O2] =  2.01320000e+02;
	W[sC12OOHXT6XC12H25O2] =  2.01320000e+02;
	W[sC12OOHO2XT6XC12H25O4] =  2.33320000e+02;
	W[sC12KETXT6XC12H24O3] =  2.16312000e+02;
	W[sNXC8H17CHO] =  1.42234000e+02;
	W[sC7H15CO] =  1.27200000e+02;
	W[sC3H7CHO] =  7.21040000e+01;
	W[sNXC7H15CHO] =  1.28208000e+02;
	W[sC7H15O2X1XC7H15O2] =  1.31190000e+02;
	W[sC7H14OOH1X3XC7H15O2] =  1.31190000e+02;
	W[sC5H4CH2] =  7.81080000e+01;
	W[sA1XC6H6] =  7.81080000e+01;
	W[sNXC4H5] =  5.30800000e+01;
	W[sIXC4H5] =  5.30800000e+01;
	W[sA1XXC6H5] =  7.71000000e+01;
	W[sNXC4H3] =  5.10640000e+01;
	W[sC3H3] =  3.90540000e+01;
	W[sA1C2H2XC8H7] =  1.03136000e+02;
	W[sA1C2H3XC8H8] =  1.04144000e+02;
	W[sA1C2HXC8H6] =  1.02128000e+02;
	W[sA1C2HYXC8H5] =  1.01120000e+02;
	W[sA1C2H3YXC8H7] =  1.03136000e+02;
	W[sA2XXC10H7] =  1.27156000e+02;
	W[sA2XC10H8] =  1.28164000e+02;
	W[sC4H4] =  5.20720000e+01;
	W[sA2YXC10H7] =  1.27156000e+02;
	W[sA2C2H2AXC12H9] =  1.53192000e+02;
	W[sA2C2H2BXC12H9] =  1.53192000e+02;
	W[sA2C2HAXC12H8] =  1.52184000e+02;
	W[sA2C2HBXC12H8] =  1.52184000e+02;
	W[sA2C2HAYXC12H7] =  1.51176000e+02;
	W[sA2C2HBYXC12H7] =  1.51176000e+02;
	W[sA2R5XC12H8] =  1.52184000e+02;
	W[sA2R5XXC12H7] =  1.51176000e+02;
	W[sA2R5C2H2XC14H9] =  1.77212000e+02;
	W[sA2R5C2HXC14H8] =  1.76204000e+02;
	W[sA2R5C2HYXC14H7] =  1.75196000e+02;
	W[sP2XC12H10] =  1.54200000e+02;
	W[sP2XXC12H9] =  1.53192000e+02;
	W[sA3XXC14H9] =  1.77212000e+02;
	W[sA3XC14H10] =  1.78220000e+02;
	W[sA3YXC14H9] =  1.77212000e+02;
	W[sA3R5XXC16H9] =  2.01232000e+02;
	W[sA3R5XC16H10] =  2.02240000e+02;
	W[sA4XC16H10] =  2.02240000e+02;
	W[sA4XXC16H9] =  2.01232000e+02;
	W[sA4R5XC18H10] =  2.26260000e+02;
	W[sFLTNXC16H10] =  2.02240000e+02;
	W[sC5H6] =  6.60980000e+01;
	W[sC5H5] =  6.50900000e+01;
	W[sSXC3H5] =  4.10700000e+01;
	W[sCH4] =  1.60420000e+01;
	W[sC4H6] =  5.40880000e+01;
	W[sTXC5H5O] =  8.10900000e+01;
	W[sC5H4O] =  8.00820000e+01;
	W[sSXC5H5O] =  8.10900000e+01;
	W[sC9H8] =  1.16154000e+02;
	W[sC9H7] =  1.15146000e+02;
	W[sA1CH2XC7H7] =  9.11260000e+01;
	W[sC9H6O] =  1.30138000e+02;
	W[sOXC6H4] =  7.60920000e+01;
	W[sA1CH3XC7H8] =  9.21340000e+01;
	W[sA1OHXC6H6O] =  9.41080000e+01;
	W[sHOA1CH3XC7H8O] =  1.08134000e+02;
	W[sOA1CH3XC7H7O] =  1.07126000e+02;
	W[sA1CH2OXC7H7O] =  1.07126000e+02;
	W[sA1CH2OHXC7H8O] =  1.08134000e+02;
	W[sA1CHOXC7H6O] =  1.06118000e+02;
	W[sA1OXC6H5O] =  9.31000000e+01;
	W[sA1CH3YXC7H7] =  9.11260000e+01;
	W[sPXC3H4] =  4.00620000e+01;
	W[sA1C2H4XC8H9] =  1.05152000e+02;
	W[sA1C2H5XC8H10] =  1.06160000e+02;
	W[sC8H9O2] =  1.37152000e+02;
	W[sC8H8OOH] =  1.37152000e+02;
	W[sOC8H7OOH] =  1.52144000e+02;
	W[sH2C2] =  2.60360000e+01;
	W[sA1CH3CH3XC8H10] =  1.06160000e+02;
	W[sA1CH3CH2XC8H9] =  1.05152000e+02;
	W[sA1CH3CHOXC8H8O] =  1.20144000e+02;
	W[sA2CH3XC11H10] =  1.42190000e+02;
	W[sA1CHOCH2XC8H7O] =  1.19136000e+02;
	W[sA1CHOCHOXC8H6O2] =  1.34128000e+02;
	W[sA2OHXC10H8O] =  1.44164000e+02;
	W[sA2CH2XC11H9] =  1.41182000e+02;
	W[sA2CH2OXC11H9O] =  1.57182000e+02;
	W[sA2CHOXC11H8O] =  1.56174000e+02;
	W[sA2OXC10H7O] =  1.43156000e+02;
	W[sC4H2] =  5.00560000e+01;
	W[sOC6H4O] =  1.08092000e+02;
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
	W[sC] =  1.20100000e+01;
	W[sHCNN] =  4.10380000e+01;
	W[sHCNO] =  4.30280000e+01;
	W[sNH3] =  1.70340000e+01;
}

void GetSpeciesNames( char **names )
{
	strcpy( names[sN2], "N2                  " );
	strcpy( names[sAR], "AR                  " );
	strcpy( names[sO], "O                   " );
	strcpy( names[sH2], "H2                  " );
	strcpy( names[sH], "H                   " );
	strcpy( names[sOH], "OH                  " );
	strcpy( names[sH2O], "H2O                 " );
	strcpy( names[sO2], "O2                  " );
	strcpy( names[sHO2], "HO2                 " );
	strcpy( names[sH2O2], "H2O2                " );
	strcpy( names[sCH], "CH                  " );
	strcpy( names[sTXCH2], "T-CH2               " );
	strcpy( names[sCH2O], "CH2O                " );
	strcpy( names[sHCO], "HCO                 " );
	strcpy( names[sCH3], "CH3                 " );
	strcpy( names[sCO2], "CO2                 " );
	strcpy( names[sCO], "CO                  " );
	strcpy( names[sSXCH2], "S-CH2               " );
	strcpy( names[sC2H4], "C2H4                " );
	strcpy( names[sCH3O], "CH3O                " );
	strcpy( names[sCH3O2], "CH3O2               " );
	strcpy( names[sC2H5], "C2H5                " );
	strcpy( names[sCH2CO], "CH2CO               " );
	strcpy( names[sO2CHO], "O2CHO               " );
	strcpy( names[sHOCHO], "HOCHO               " );
	strcpy( names[sOCHO], "OCHO                " );
	strcpy( names[sC2H2], "C2H2                " );
	strcpy( names[sC2H3], "C2H3                " );
	strcpy( names[sHCCO], "HCCO                " );
	strcpy( names[sCH2CHO], "CH2CHO              " );
	strcpy( names[sC2H5O], "C2H5O               " );
	strcpy( names[sAXC3H5], "A-C3H5              " );
	strcpy( names[sC3H6], "C3H6                " );
	strcpy( names[sCH3COCH2], "CH3COCH2            " );
	strcpy( names[sNXC3H7], "N-C3H7              " );
	strcpy( names[sPXC4H8], "P-C4H8              " );
	strcpy( names[sPXC4H9], "P-C4H9              " );
	strcpy( names[sC5H11], "C5H11               " );
	strcpy( names[sC6H12X1XC6H12], "C6H12-1-C6H12       " );
	strcpy( names[sC6H13], "C6H13               " );
	strcpy( names[sC7H14], "C7H14               " );
	strcpy( names[sC7H15], "C7H15               " );
	strcpy( names[sC8H17], "C8H17               " );
	strcpy( names[sC9H19], "C9H19               " );
	strcpy( names[sNXC12H26], "N-C12H26            " );
	strcpy( names[sC12H25], "C12H25              " );
	strcpy( names[sC12H24], "C12H24              " );
	strcpy( names[sC12H25O2], "C12H25O2            " );
	strcpy( names[sC12OOHXT6XC12H25O2], "C12OOH-T6-C12H25O2  " );
	strcpy( names[sC12OOHO2XT6XC12H25O4], "C12OOHO2-T6-C12H25O4" );
	strcpy( names[sC12KETXT6XC12H24O3], "C12KET-T6-C12H24O3  " );
	strcpy( names[sNXC8H17CHO], "N-C8H17CHO          " );
	strcpy( names[sC7H15CO], "C7H15CO             " );
	strcpy( names[sC3H7CHO], "C3H7CHO             " );
	strcpy( names[sNXC7H15CHO], "N-C7H15CHO          " );
	strcpy( names[sC7H15O2X1XC7H15O2], "C7H15O2-1-C7H15O2   " );
	strcpy( names[sC7H14OOH1X3XC7H15O2], "C7H14OOH1-3-C7H15O2 " );
	strcpy( names[sC5H4CH2], "C5H4CH2             " );
	strcpy( names[sA1XC6H6], "A1-C6H6             " );
	strcpy( names[sNXC4H5], "N-C4H5              " );
	strcpy( names[sIXC4H5], "I-C4H5              " );
	strcpy( names[sA1XXC6H5], "A1--C6H5            " );
	strcpy( names[sNXC4H3], "N-C4H3              " );
	strcpy( names[sC3H3], "C3H3                " );
	strcpy( names[sA1C2H2XC8H7], "A1C2H2-C8H7         " );
	strcpy( names[sA1C2H3XC8H8], "A1C2H3-C8H8         " );
	strcpy( names[sA1C2HXC8H6], "A1C2H-C8H6          " );
	strcpy( names[sA1C2HYXC8H5], "A1C2H*-C8H5         " );
	strcpy( names[sA1C2H3YXC8H7], "A1C2H3*-C8H7        " );
	strcpy( names[sA2XXC10H7], "A2--C10H7           " );
	strcpy( names[sA2XC10H8], "A2-C10H8            " );
	strcpy( names[sC4H4], "C4H4                " );
	strcpy( names[sA2YXC10H7], "A2*-C10H7           " );
	strcpy( names[sA2C2H2AXC12H9], "A2C2H2A-C12H9       " );
	strcpy( names[sA2C2H2BXC12H9], "A2C2H2B-C12H9       " );
	strcpy( names[sA2C2HAXC12H8], "A2C2HA-C12H8        " );
	strcpy( names[sA2C2HBXC12H8], "A2C2HB-C12H8        " );
	strcpy( names[sA2C2HAYXC12H7], "A2C2HA*-C12H7       " );
	strcpy( names[sA2C2HBYXC12H7], "A2C2HB*-C12H7       " );
	strcpy( names[sA2R5XC12H8], "A2R5-C12H8          " );
	strcpy( names[sA2R5XXC12H7], "A2R5--C12H7         " );
	strcpy( names[sA2R5C2H2XC14H9], "A2R5C2H2-C14H9      " );
	strcpy( names[sA2R5C2HXC14H8], "A2R5C2H-C14H8       " );
	strcpy( names[sA2R5C2HYXC14H7], "A2R5C2H*-C14H7      " );
	strcpy( names[sP2XC12H10], "P2-C12H10           " );
	strcpy( names[sP2XXC12H9], "P2--C12H9           " );
	strcpy( names[sA3XXC14H9], "A3--C14H9           " );
	strcpy( names[sA3XC14H10], "A3-C14H10           " );
	strcpy( names[sA3YXC14H9], "A3*-C14H9           " );
	strcpy( names[sA3R5XXC16H9], "A3R5--C16H9         " );
	strcpy( names[sA3R5XC16H10], "A3R5-C16H10         " );
	strcpy( names[sA4XC16H10], "A4-C16H10           " );
	strcpy( names[sA4XXC16H9], "A4--C16H9           " );
	strcpy( names[sA4R5XC18H10], "A4R5-C18H10         " );
	strcpy( names[sFLTNXC16H10], "FLTN-C16H10         " );
	strcpy( names[sC5H6], "C5H6                " );
	strcpy( names[sC5H5], "C5H5                " );
	strcpy( names[sSXC3H5], "S-C3H5              " );
	strcpy( names[sCH4], "CH4                 " );
	strcpy( names[sC4H6], "C4H6                " );
	strcpy( names[sTXC5H5O], "T-C5H5O             " );
	strcpy( names[sC5H4O], "C5H4O               " );
	strcpy( names[sSXC5H5O], "S-C5H5O             " );
	strcpy( names[sC9H8], "C9H8                " );
	strcpy( names[sC9H7], "C9H7                " );
	strcpy( names[sA1CH2XC7H7], "A1CH2-C7H7          " );
	strcpy( names[sC9H6O], "C9H6O               " );
	strcpy( names[sOXC6H4], "O-C6H4              " );
	strcpy( names[sA1CH3XC7H8], "A1CH3-C7H8          " );
	strcpy( names[sA1OHXC6H6O], "A1OH-C6H6O          " );
	strcpy( names[sHOA1CH3XC7H8O], "HOA1CH3-C7H8O       " );
	strcpy( names[sOA1CH3XC7H7O], "OA1CH3-C7H7O        " );
	strcpy( names[sA1CH2OXC7H7O], "A1CH2O-C7H7O        " );
	strcpy( names[sA1CH2OHXC7H8O], "A1CH2OH-C7H8O       " );
	strcpy( names[sA1CHOXC7H6O], "A1CHO-C7H6O         " );
	strcpy( names[sA1OXC6H5O], "A1O-C6H5O           " );
	strcpy( names[sA1CH3YXC7H7], "A1CH3*-C7H7         " );
	strcpy( names[sPXC3H4], "P-C3H4              " );
	strcpy( names[sA1C2H4XC8H9], "A1C2H4-C8H9         " );
	strcpy( names[sA1C2H5XC8H10], "A1C2H5-C8H10        " );
	strcpy( names[sC8H9O2], "C8H9O2              " );
	strcpy( names[sC8H8OOH], "C8H8OOH             " );
	strcpy( names[sOC8H7OOH], "OC8H7OOH            " );
	strcpy( names[sH2C2], "H2C2                " );
	strcpy( names[sA1CH3CH3XC8H10], "A1CH3CH3-C8H10      " );
	strcpy( names[sA1CH3CH2XC8H9], "A1CH3CH2-C8H9       " );
	strcpy( names[sA1CH3CHOXC8H8O], "A1CH3CHO-C8H8O      " );
	strcpy( names[sA2CH3XC11H10], "A2CH3-C11H10        " );
	strcpy( names[sA1CHOCH2XC8H7O], "A1CHOCH2-C8H7O      " );
	strcpy( names[sA1CHOCHOXC8H6O2], "A1CHOCHO-C8H6O2     " );
	strcpy( names[sA2OHXC10H8O], "A2OH-C10H8O         " );
	strcpy( names[sA2CH2XC11H9], "A2CH2-C11H9         " );
	strcpy( names[sA2CH2OXC11H9O], "A2CH2O-C11H9O       " );
	strcpy( names[sA2CHOXC11H8O], "A2CHO-C11H8O        " );
	strcpy( names[sA2OXC10H7O], "A2O-C10H7O          " );
	strcpy( names[sC4H2], "C4H2                " );
	strcpy( names[sOC6H4O], "OC6H4O              " );
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
	strcpy( names[sC], "C                   " );
	strcpy( names[sHCNN], "HCNN                " );
	strcpy( names[sHCNO], "HCNO                " );
	strcpy( names[sNH3], "NH3                 " );
}



void ComputeThermoData( double *h, double *cp, double T )
{
/*
	This function computes enthalpy 'h' and heat capacity 'cp' as
	function of temperature 'T' for all non steady state species
	in units [J/kg] and [J/kg K], respectively.
	The parameter h and cp should provide workspace of length 155 */


	int i;
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
		h[sO] =  5.19646250e+02 * (
			T * (  2.54206000e+00 + T * ( -1.37753100e-05
			+ T * ( -1.03426767e-09 + T * (  1.13776675e-12
			+ T * -8.73610400e-17 ) ) ) ) +  2.92308000e+04 );
		cp[sO] =  5.19646250e+02 * (
			 2.54206000e+00 + T * ( -2.75506200e-05 
			+ T * ( -3.10280300e-09 + T * (  4.55106700e-12
			+ T * -4.36805200e-16 ) ) ) );
		h[sH2] =  4.12417659e+03 * (
			T * (  2.99142300e+00 + T * (  3.50032200e-04
			+ T * ( -1.87794300e-08 + T * ( -2.30789450e-12
			+ T *  3.16550400e-16 ) ) ) ) -8.35034000e+02 );
		cp[sH2] =  4.12417659e+03 * (
			 2.99142300e+00 + T * (  7.00064400e-04 
			+ T * ( -5.63382900e-08 + T * ( -9.23157800e-12
			+ T *  1.58275200e-15 ) ) ) );
		h[sH] =  8.24835317e+03 * (
			T * (  2.50000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  2.54716300e+04 );
		cp[sH] =  8.24835317e+03 * (
			 2.50000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sOH] =  4.88848777e+02 * (
			T * (  2.86472886e+00 + T * (  5.28252240e-04
			+ T * ( -8.63609193e-08 + T * (  7.63046685e-12
			+ T * -2.66391752e-16 ) ) ) ) +  3.68362875e+03 );
		cp[sOH] =  4.88848777e+02 * (
			 2.86472886e+00 + T * (  1.05650448e-03 
			+ T * ( -2.59082758e-07 + T * (  3.05218674e-11
			+ T * -1.33195876e-15 ) ) ) );
		h[sH2O] =  4.61497558e+02 * (
			T * (  2.67214600e+00 + T * (  1.52814650e-03
			+ T * ( -2.91008667e-07 + T * (  3.00249000e-11
			+ T * -1.27832360e-15 ) ) ) ) -2.98992100e+04 );
		cp[sH2O] =  4.61497558e+02 * (
			 2.67214600e+00 + T * (  3.05629300e-03 
			+ T * ( -8.73026000e-07 + T * (  1.20099600e-10
			+ T * -6.39161800e-15 ) ) ) );
		h[sO2] =  2.59823125e+02 * (
			T * (  3.69757800e+00 + T * (  3.06759850e-04
			+ T * ( -4.19614000e-08 + T * (  4.43820250e-12
			+ T * -2.27287000e-16 ) ) ) ) -1.23393000e+03 );
		cp[sO2] =  2.59823125e+02 * (
			 3.69757800e+00 + T * (  6.13519700e-04 
			+ T * ( -1.25884200e-07 + T * (  1.77528100e-11
			+ T * -1.13643500e-15 ) ) ) );
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
		h[sCH] =  6.38680289e+02 * (
			T * (  2.87846473e+00 + T * (  4.85456840e-04
			+ T * (  4.81485517e-08 + T * ( -3.26719623e-11
			+ T *  3.52158766e-15 ) ) ) ) +  7.10124364e+04 );
		cp[sCH] =  6.38680289e+02 * (
			 2.87846473e+00 + T * (  9.70913681e-04 
			+ T * (  1.44445655e-07 + T * ( -1.30687849e-10
			+ T *  1.76079383e-14 ) ) ) );
		h[sTXCH2] =  5.92780550e+02 * (
			T * (  2.87410113e+00 + T * (  1.82819646e-03
			+ T * ( -4.69648657e-07 + T * (  6.50448872e-11
			+ T * -3.75455134e-15 ) ) ) ) +  4.62636040e+04 );
		cp[sTXCH2] =  5.92780550e+02 * (
			 2.87410113e+00 + T * (  3.65639292e-03 
			+ T * ( -1.40894597e-06 + T * (  2.60179549e-10
			+ T * -1.87727567e-14 ) ) ) );
		h[sCH2O] =  2.76904683e+02 * (
			T * (  1.76069008e+00 + T * (  4.60000041e-03
			+ T * ( -1.47419604e-06 + T * (  2.51603030e-10
			+ T * -1.76771128e-14 ) ) ) ) -1.39958323e+04 );
		cp[sCH2O] =  2.76904683e+02 * (
			 1.76069008e+00 + T * (  9.20000082e-03 
			+ T * ( -4.42258813e-06 + T * (  1.00641212e-09
			+ T * -8.83855640e-14 ) ) ) );
		h[sHCO] =  2.86523537e+02 * (
			T * (  2.77217438e+00 + T * (  2.47847763e-03
			+ T * ( -8.28152043e-07 + T * (  1.47290445e-10
			+ T * -1.06701742e-14 ) ) ) ) +  4.01191815e+03 );
		cp[sHCO] =  2.86523537e+02 * (
			 2.77217438e+00 + T * (  4.95695526e-03 
			+ T * ( -2.48445613e-06 + T * (  5.89161778e-10
			+ T * -5.33508711e-14 ) ) ) );
		h[sCH3] =  5.53035786e+02 * (
			T * (  2.97812060e+00 + T * (  2.89892600e-03
			+ T * ( -6.58526667e-07 + T * (  7.68244750e-11
			+ T * -3.58348320e-15 ) ) ) ) +  1.65095130e+04 );
		cp[sCH3] =  5.53035786e+02 * (
			 2.97812060e+00 + T * (  5.79785200e-03 
			+ T * ( -1.97558000e-06 + T * (  3.07297900e-10
			+ T * -1.79174160e-14 ) ) ) );
		h[sCO2] =  1.88919337e+02 * (
			T * (  4.45362300e+00 + T * (  1.57008450e-03
			+ T * ( -4.26137000e-07 + T * (  5.98499250e-11
			+ T * -3.33806600e-15 ) ) ) ) -4.89669600e+04 );
		cp[sCO2] =  1.88919337e+02 * (
			 4.45362300e+00 + T * (  3.14016900e-03 
			+ T * ( -1.27841100e-06 + T * (  2.39399700e-10
			+ T * -1.66903300e-14 ) ) ) );
		h[sCO] =  2.96834702e+02 * (
			T * (  3.02507800e+00 + T * (  7.21344500e-04
			+ T * ( -1.87694267e-07 + T * (  2.54645250e-11
			+ T * -1.38219040e-15 ) ) ) ) -1.42683500e+04 );
		cp[sCO] =  2.96834702e+02 * (
			 3.02507800e+00 + T * (  1.44268900e-03 
			+ T * ( -5.63082800e-07 + T * (  1.01858100e-10
			+ T * -6.91095200e-15 ) ) ) );
		h[sSXCH2] =  5.92780550e+02 * (
			T * (  2.29203842e+00 + T * (  2.32794318e-03
			+ T * ( -6.70639823e-07 + T * (  1.04476500e-10
			+ T * -6.79432730e-15 ) ) ) ) +  5.09259997e+04 );
		cp[sSXCH2] =  5.92780550e+02 * (
			 2.29203842e+00 + T * (  4.65588637e-03 
			+ T * ( -2.01191947e-06 + T * (  4.17906000e-10
			+ T * -3.39716365e-14 ) ) ) );
		h[sC2H4] =  2.96390275e+02 * (
			T * (  2.03611116e+00 + T * (  7.32270755e-03
			+ T * ( -2.23692638e-06 + T * (  3.68057308e-10
			+ T * -2.51412122e-14 ) ) ) ) +  4.93988614e+03 );
		cp[sC2H4] =  2.96390275e+02 * (
			 2.03611116e+00 + T * (  1.46454151e-02 
			+ T * ( -6.71077915e-06 + T * (  1.47222923e-09
			+ T * -1.25706061e-13 ) ) ) );
		h[sCH3O] =  2.67910679e+02 * (
			T * (  3.77079900e+00 + T * (  3.93574850e-03
			+ T * ( -8.85461333e-07 + T * (  9.86107750e-11
			+ T * -4.22523200e-15 ) ) ) ) +  1.27832520e+02 );
		cp[sCH3O] =  2.67910679e+02 * (
			 3.77079900e+00 + T * (  7.87149700e-03 
			+ T * ( -2.65638400e-06 + T * (  3.94443100e-10
			+ T * -2.11261600e-14 ) ) ) );
		h[sCH3O2] =  1.76772973e+02 * (
			T * (  5.92505819e+00 + T * (  4.50097271e-03
			+ T * ( -1.08084770e-06 + T * (  1.31090679e-10
			+ T * -6.28526006e-15 ) ) ) ) -1.53258958e+03 );
		cp[sCH3O2] =  1.76772973e+02 * (
			 5.92505819e+00 + T * (  9.00194542e-03 
			+ T * ( -3.24254309e-06 + T * (  5.24362718e-10
			+ T * -3.14263003e-14 ) ) ) );
		h[sC2H5] =  2.86109429e+02 * (
			T * (  1.95465642e+00 + T * (  8.69863610e-03
			+ T * ( -2.66068889e-06 + T * (  4.38044223e-10
			+ T * -2.99283152e-14 ) ) ) ) +  1.28575200e+04 );
		cp[sC2H5] =  2.86109429e+02 * (
			 1.95465642e+00 + T * (  1.73972722e-02 
			+ T * ( -7.98206668e-06 + T * (  1.75217689e-09
			+ T * -1.49641576e-13 ) ) ) );
		h[sCH2CO] =  1.97790941e+02 * (
			T * (  4.51129732e+00 + T * (  4.50179872e-03
			+ T * ( -1.38979878e-06 + T * (  2.30836470e-10
			+ T * -1.58967640e-14 ) ) ) ) -7.55105311e+03 );
		cp[sCH2CO] =  1.97790941e+02 * (
			 4.51129732e+00 + T * (  9.00359745e-03 
			+ T * ( -4.16939635e-06 + T * (  9.23345882e-10
			+ T * -7.94838201e-14 ) ) ) );
		h[sO2CHO] =  1.36260448e+02 * (
			T * (  4.05631532e+00 + T * (  5.14596660e-03
			+ T * ( -1.63298733e-06 + T * (  2.10704868e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.73830665e+04 );
		cp[sO2CHO] =  1.36260448e+02 * (
			 4.05631532e+00 + T * (  1.02919332e-02 
			+ T * ( -4.89896198e-06 + T * (  8.42819470e-10
			+ T *  0.00000000e+00 ) ) ) );
		h[sHOCHO] =  1.80644418e+02 * (
			T * (  1.24573687e+00 + T * (  8.36210310e-03
			+ T * ( -3.07392626e-06 + T * (  4.41095550e-10
			+ T *  0.00000000e+00 ) ) ) ) -4.65097525e+04 );
		cp[sHOCHO] =  1.80644418e+02 * (
			 1.24573687e+00 + T * (  1.67242062e-02 
			+ T * ( -9.22177878e-06 + T * (  1.76438220e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sOCHO] =  1.84689235e+02 * (
			T * (  1.69314163e+00 + T * (  6.21628800e-03
			+ T * ( -2.29557633e-06 + T * (  3.30393593e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.63316485e+04 );
		cp[sOCHO] =  1.84689235e+02 * (
			 1.69314163e+00 + T * (  1.24325760e-02 
			+ T * ( -6.88672900e-06 + T * (  1.32157437e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC2H2] =  3.19340144e+02 * (
			T * (  4.14756964e+00 + T * (  2.98083332e-03
			+ T * ( -7.90982840e-07 + T * (  1.16853043e-10
			+ T * -7.22470426e-15 ) ) ) ) +  2.59359992e+04 );
		cp[sC2H2] =  3.19340144e+02 * (
			 4.14756964e+00 + T * (  5.96166664e-03 
			+ T * ( -2.37294852e-06 + T * (  4.67412171e-10
			+ T * -3.61235213e-14 ) ) ) );
		h[sC2H3] =  3.07437509e+02 * (
			T * (  3.01672400e+00 + T * (  5.16511460e-03
			+ T * ( -1.56027450e-06 + T * (  2.54408220e-10
			+ T * -1.72521408e-14 ) ) ) ) +  3.46128739e+04 );
		cp[sC2H3] =  3.07437509e+02 * (
			 3.01672400e+00 + T * (  1.03302292e-02 
			+ T * ( -4.68082349e-06 + T * (  1.01763288e-09
			+ T * -8.62607041e-14 ) ) ) );
		h[sHCCO] =  2.02650385e+02 * (
			T * (  5.62820580e+00 + T * (  2.04267005e-03
			+ T * ( -5.31151567e-07 + T * (  7.15651300e-11
			+ T * -3.88156640e-15 ) ) ) ) +  1.93272150e+04 );
		cp[sHCCO] =  2.02650385e+02 * (
			 5.62820580e+00 + T * (  4.08534010e-03 
			+ T * ( -1.59345470e-06 + T * (  2.86260520e-10
			+ T * -1.94078320e-14 ) ) ) );
		h[sCH2CHO] =  1.93159093e+02 * (
			T * (  2.42606357e+00 + T * (  8.62000105e-03
			+ T * ( -3.25710706e-06 + T * (  6.66389180e-10
			+ T * -5.64240156e-14 ) ) ) ) +  8.33106990e+02 );
		cp[sCH2CHO] =  1.93159093e+02 * (
			 2.42606357e+00 + T * (  1.72400021e-02 
			+ T * ( -9.77132119e-06 + T * (  2.66555672e-09
			+ T * -2.82120078e-13 ) ) ) );
		h[sC2H5O] =  1.84517088e+02 * (
			T * (  2.46262349e+00 + T * (  1.04751979e-02
			+ T * ( -3.13097250e-06 + T * (  3.91101567e-10
			+ T *  0.00000000e+00 ) ) ) ) -3.83932658e+03 );
		cp[sC2H5O] =  1.84517088e+02 * (
			 2.46262349e+00 + T * (  2.09503959e-02 
			+ T * ( -9.39291750e-06 + T * (  1.56440627e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sAXC3H5] =  2.02443146e+02 * (
			T * (  2.28794927e+00 + T * (  1.18200788e-02
			+ T * ( -4.26304833e-06 + T * (  8.42096350e-10
			+ T * -6.94898898e-14 ) ) ) ) +  1.83033514e+04 );
		cp[sAXC3H5] =  2.02443146e+02 * (
			 2.28794927e+00 + T * (  2.36401575e-02 
			+ T * ( -1.27891450e-05 + T * (  3.36838540e-09
			+ T * -3.47449449e-13 ) ) ) );
		h[sC3H6] =  1.97593517e+02 * (
			T * (  4.71697982e-01 + T * (  1.44756535e-02
			+ T * ( -5.22006063e-06 + T * (  1.02860800e-09
			+ T * -8.46150282e-14 ) ) ) ) +  1.12603387e+03 );
		cp[sC3H6] =  1.97593517e+02 * (
			 4.71697982e-01 + T * (  2.89513070e-02 
			+ T * ( -1.56601819e-05 + T * (  4.11443199e-09
			+ T * -4.23075141e-13 ) ) ) );
		h[sCH3COCH2] =  1.45686701e+02 * (
			T * (  4.12743010e+00 + T * (  1.16865282e-02
			+ T * ( -3.66800960e-06 + T * (  4.73988545e-10
			+ T *  0.00000000e+00 ) ) ) ) -7.31939257e+03 );
		cp[sCH3COCH2] =  1.45686701e+02 * (
			 4.12743010e+00 + T * (  2.33730564e-02 
			+ T * ( -1.10040288e-05 + T * (  1.89595418e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sNXC3H7] =  1.92970803e+02 * (
			T * (  7.70404050e+00 + T * (  8.02077000e-03
			+ T * ( -1.76053223e-06 + T * (  1.90636007e-10
			+ T * -7.87069240e-15 ) ) ) ) +  8.29795310e+03 );
		cp[sNXC3H7] =  1.92970803e+02 * (
			 7.70404050e+00 + T * (  1.60415400e-02 
			+ T * ( -5.28159670e-06 + T * (  7.62544030e-10
			+ T * -3.93534620e-14 ) ) ) );
		h[sPXC4H8] =  1.48195138e+02 * (
			T * (  3.04470367e+00 + T * (  1.63725883e-02
			+ T * ( -4.84544123e-06 + T * (  5.99360043e-10
			+ T *  0.00000000e+00 ) ) ) ) -2.52177534e+03 );
		cp[sPXC4H8] =  1.48195138e+02 * (
			 3.04470367e+00 + T * (  3.27451765e-02 
			+ T * ( -1.45363237e-05 + T * (  2.39744017e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sPXC4H9] =  1.45579563e+02 * (
			T * (  3.81812330e+00 + T * (  1.70395244e-02
			+ T * ( -4.97118110e-06 + T * (  6.08013118e-10
			+ T *  0.00000000e+00 ) ) ) ) +  6.65901304e+03 );
		cp[sPXC4H9] =  1.45579563e+02 * (
			 3.81812330e+00 + T * (  3.40790489e-02 
			+ T * ( -1.49135433e-05 + T * (  2.43205247e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC5H11] =  1.16876212e+02 * (
			T * (  4.88920629e+00 + T * (  2.11417269e-02
			+ T * ( -6.19477000e-06 + T * (  7.60311908e-10
			+ T *  0.00000000e+00 ) ) ) ) +  3.43475468e+03 );
		cp[sC5H11] =  1.16876212e+02 * (
			 4.88920629e+00 + T * (  4.22834537e-02 
			+ T * ( -1.85843100e-05 + T * (  3.04124763e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC6H12X1XC6H12] =  9.87967584e+01 * (
			T * (  4.97075871e+00 + T * (  2.48318661e-02
			+ T * ( -7.43158390e-06 + T * (  9.29005450e-10
			+ T *  0.00000000e+00 ) ) ) ) -8.89572273e+03 );
		cp[sC6H12X1XC6H12] =  9.87967584e+01 * (
			 4.97075871e+00 + T * (  4.96637321e-02 
			+ T * ( -2.22947517e-05 + T * (  3.71602180e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC6H13] =  9.76274012e+01 * (
			T * (  4.66412018e+00 + T * (  2.65476469e-02
			+ T * ( -7.98363293e-06 + T * (  1.00158224e-09
			+ T *  0.00000000e+00 ) ) ) ) +  6.21004203e+02 );
		cp[sC6H13] =  9.76274012e+01 * (
			 4.66412018e+00 + T * (  5.30952938e-02 
			+ T * ( -2.39508988e-05 + T * (  4.00632895e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC7H14] =  8.46829358e+01 * (
			T * (  5.45858240e+00 + T * (  2.93078813e-02
			+ T * ( -8.76964720e-06 + T * (  1.09504683e-09
			+ T *  0.00000000e+00 ) ) ) ) -1.33463299e+04 );
		cp[sC7H14] =  8.46829358e+01 * (
			 5.45858240e+00 + T * (  5.86157625e-02 
			+ T * ( -2.63089416e-05 + T * (  4.38018732e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC7H15] =  8.38223611e+01 * (
			T * (  3.74721159e+00 + T * (  3.24672581e-02
			+ T * ( -1.00447008e-05 + T * (  1.29354536e-09
			+ T *  0.00000000e+00 ) ) ) ) -3.37018357e+03 );
		cp[sC7H15] =  8.38223611e+01 * (
			 3.74721159e+00 + T * (  6.49345162e-02 
			+ T * ( -3.01341025e-05 + T * (  5.17418142e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC8H17] =  7.34378533e+01 * (
			T * (  6.64595645e+00 + T * (  3.49584911e-02
			+ T * ( -1.05512069e-05 + T * (  1.32847881e-09
			+ T *  0.00000000e+00 ) ) ) ) -5.77046319e+03 );
		cp[sC8H17] =  7.34378533e+01 * (
			 6.64595645e+00 + T * (  6.99169822e-02 
			+ T * ( -3.16536208e-05 + T * (  5.31391524e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC9H19] =  6.53427327e+01 * (
			T * (  7.61279482e+00 + T * (  3.91887732e-02
			+ T * ( -1.18454948e-05 + T * (  1.49353250e-09
			+ T *  0.00000000e+00 ) ) ) ) -8.95900903e+03 );
		cp[sC9H19] =  6.53427327e+01 * (
			 7.61279482e+00 + T * (  7.83775464e-02 
			+ T * ( -3.55364844e-05 + T * (  5.97413001e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sNXC12H26] =  4.88137006e+01 * (
			T * (  9.97283422e+00 + T * (  5.38076480e-02
			+ T * ( -1.63000417e-05 + T * (  2.05974186e-09
			+ T *  0.00000000e+00 ) ) ) ) -4.31954165e+04 );
		cp[sNXC12H26] =  4.88137006e+01 * (
			 9.97283422e+00 + T * (  1.07615296e-01 
			+ T * ( -4.89001251e-05 + T * (  8.23896743e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC12H25] =  4.91042996e+01 * (
			T * (  1.05038381e+01 + T * (  5.18907195e-02
			+ T * ( -1.57336507e-05 + T * (  1.98955892e-09
			+ T *  0.00000000e+00 ) ) ) ) -1.85218187e+04 );
		cp[sC12H25] =  4.91042996e+01 * (
			 1.05038381e+01 + T * (  1.03781439e-01 
			+ T * ( -4.72009520e-05 + T * (  7.95823570e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC12H24] =  4.93983792e+01 * (
			T * (  9.89916240e+00 + T * (  5.08066350e-02
			+ T * ( -1.54033264e-05 + T * (  1.94782579e-09
			+ T *  0.00000000e+00 ) ) ) ) -2.92166047e+04 );
		cp[sC12H24] =  4.93983792e+01 * (
			 9.89916240e+00 + T * (  1.01613270e-01 
			+ T * ( -4.62099791e-05 + T * (  7.79130314e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC12H25O2] =  4.12991258e+01 * (
			T * (  1.23894782e+01 + T * (  5.40244640e-02
			+ T * ( -1.62884320e-05 + T * (  2.03850321e-09
			+ T *  0.00000000e+00 ) ) ) ) -3.76523809e+04 );
		cp[sC12H25O2] =  4.12991258e+01 * (
			 1.23894782e+01 + T * (  1.08048928e-01 
			+ T * ( -4.88652961e-05 + T * (  8.15401282e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC12OOHXT6XC12H25O2] =  4.12991258e+01 * (
			T * (  1.27466803e+01 + T * (  5.45462995e-02
			+ T * ( -1.68126854e-05 + T * (  2.15175494e-09
			+ T *  0.00000000e+00 ) ) ) ) -3.13603674e+04 );
		cp[sC12OOHXT6XC12H25O2] =  4.12991258e+01 * (
			 1.27466803e+01 + T * (  1.09092599e-01 
			+ T * ( -5.04380562e-05 + T * (  8.60701975e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC12OOHO2XT6XC12H25O4] =  3.56349220e+01 * (
			T * (  1.51871349e+01 + T * (  5.79965270e-02
			+ T * ( -1.82215220e-05 + T * (  2.35498310e-09
			+ T *  0.00000000e+00 ) ) ) ) -5.13250050e+04 );
		cp[sC12OOHO2XT6XC12H25O4] =  3.56349220e+01 * (
			 1.51871349e+01 + T * (  1.15993054e-01 
			+ T * ( -5.46645659e-05 + T * (  9.41993241e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC12KETXT6XC12H24O3] =  3.84367950e+01 * (
			T * (  1.68711079e+01 + T * (  5.14824840e-02
			+ T * ( -1.55478239e-05 + T * (  1.95013682e-09
			+ T *  0.00000000e+00 ) ) ) ) -6.72051106e+04 );
		cp[sC12KETXT6XC12H24O3] =  3.84367950e+01 * (
			 1.68711079e+01 + T * (  1.02964968e-01 
			+ T * ( -4.66434716e-05 + T * (  7.80054730e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sNXC8H17CHO] =  5.84553623e+01 * (
			T * (  8.73373357e+00 + T * (  3.87105596e-02
			+ T * ( -1.16666906e-05 + T * (  1.46272009e-09
			+ T *  0.00000000e+00 ) ) ) ) -4.35042303e+04 );
		cp[sNXC8H17CHO] =  5.84553623e+01 * (
			 8.73373357e+00 + T * (  7.74211192e-02 
			+ T * ( -3.50000718e-05 + T * (  5.85088036e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC7H15CO] =  6.53643082e+01 * (
			T * (  8.20398940e+00 + T * (  3.27046298e-02
			+ T * ( -9.89618267e-06 + T * (  1.24500341e-09
			+ T *  0.00000000e+00 ) ) ) ) -2.17195608e+04 );
		cp[sC7H15CO] =  6.53643082e+01 * (
			 8.20398940e+00 + T * (  6.54092597e-02 
			+ T * ( -2.96885480e-05 + T * (  4.98001365e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3H7CHO] =  1.15310385e+02 * (
			T * (  3.99143562e+00 + T * (  1.76801297e-02
			+ T * ( -5.37657257e-06 + T * (  6.80257197e-10
			+ T *  0.00000000e+00 ) ) ) ) -2.76351912e+04 );
		cp[sC3H7CHO] =  1.15310385e+02 * (
			 3.99143562e+00 + T * (  3.53602593e-02 
			+ T * ( -1.61297177e-05 + T * (  2.72102879e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sNXC7H15CHO] =  6.48503994e+01 * (
			T * (  7.80908620e+00 + T * (  3.44731614e-02
			+ T * ( -1.03975574e-05 + T * (  1.30488637e-09
			+ T *  0.00000000e+00 ) ) ) ) -4.03376810e+04 );
		cp[sNXC7H15CHO] =  6.48503994e+01 * (
			 7.80908620e+00 + T * (  6.89463227e-02 
			+ T * ( -3.11926721e-05 + T * (  5.21954548e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC7H15O2X1XC7H15O2] =  6.33763244e+01 * (
			T * (  8.40579521e+00 + T * (  3.22445627e-02
			+ T * ( -9.55938713e-06 + T * (  1.17625910e-09
			+ T *  0.00000000e+00 ) ) ) ) -2.09993977e+04 );
		cp[sC7H15O2X1XC7H15O2] =  6.33763244e+01 * (
			 8.40579521e+00 + T * (  6.44891254e-02 
			+ T * ( -2.86781614e-05 + T * (  4.70503640e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC7H14OOH1X3XC7H15O2] =  6.33763244e+01 * (
			T * (  8.32815125e+00 + T * (  3.33022445e-02
			+ T * ( -1.03508678e-05 + T * (  1.33053336e-09
			+ T *  0.00000000e+00 ) ) ) ) -1.47896752e+04 );
		cp[sC7H14OOH1X3XC7H15O2] =  6.33763244e+01 * (
			 8.32815125e+00 + T * (  6.66044890e-02 
			+ T * ( -3.10526033e-05 + T * (  5.32213346e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC5H4CH2] =  1.06446715e+02 * (
			T * (  2.78194214e+00 + T * (  2.03161008e-02
			+ T * ( -7.83431090e-06 + T * (  1.62764004e-09
			+ T * -1.39361817e-13 ) ) ) ) +  2.43155607e+04 );
		cp[sC5H4CH2] =  1.06446715e+02 * (
			 2.78194214e+00 + T * (  4.06322016e-02 
			+ T * ( -2.35029327e-05 + T * (  6.51056017e-09
			+ T * -6.96809087e-13 ) ) ) );
		h[sA1XC6H6] =  1.06446715e+02 * (
			T * ( -2.06240612e-01 + T * (  2.32061220e-02
			+ T * ( -9.25511787e-06 + T * (  1.97227634e-09
			+ T * -1.72073052e-13 ) ) ) ) +  8.09883905e+03 );
		cp[sA1XC6H6] =  1.06446715e+02 * (
			-2.06240612e-01 + T * (  4.64122440e-02 
			+ T * ( -2.77653536e-05 + T * (  7.88910537e-09
			+ T * -8.60365259e-13 ) ) ) );
		h[sNXC4H5] =  1.56637905e+02 * (
			T * (  4.87674639e+00 + T * (  1.13767150e-02
			+ T * ( -3.92382327e-06 + T * (  7.38128637e-10
			+ T * -5.82913132e-14 ) ) ) ) +  4.11081097e+04 );
		cp[sNXC4H5] =  1.56637905e+02 * (
			 4.87674639e+00 + T * (  2.27534299e-02 
			+ T * ( -1.17714698e-05 + T * (  2.95251455e-09
			+ T * -2.91456566e-13 ) ) ) );
		h[sIXC4H5] =  1.56637905e+02 * (
			T * (  4.34643669e+00 + T * (  1.22880720e-02
			+ T * ( -4.36512283e-06 + T * (  8.47120312e-10
			+ T * -6.87039266e-14 ) ) ) ) +  3.58709780e+04 );
		cp[sIXC4H5] =  1.56637905e+02 * (
			 4.34643669e+00 + T * (  2.45761440e-02 
			+ T * ( -1.30953685e-05 + T * (  3.38848125e-09
			+ T * -3.43519633e-13 ) ) ) );
		h[sA1XXC6H5] =  1.07838392e+02 * (
			T * (  1.38016336e+00 + T * (  2.02016004e-02
			+ T * ( -8.07502950e-06 + T * (  1.72180830e-09
			+ T * -1.50192160e-13 ) ) ) ) +  3.86973520e+04 );
		cp[sA1XXC6H5] =  1.07838392e+02 * (
			 1.38016336e+00 + T * (  4.04032009e-02 
			+ T * ( -2.42250885e-05 + T * (  6.88723321e-09
			+ T * -7.50960802e-13 ) ) ) );
		h[sNXC4H3] =  1.62821949e+02 * (
			T * (  7.25330164e+00 + T * (  5.97904230e-03
			+ T * ( -1.75571892e-06 + T * (  2.74954688e-10
			+ T * -1.76803350e-14 ) ) ) ) +  6.28977574e+04 );
		cp[sNXC4H3] =  1.62821949e+02 * (
			 7.25330164e+00 + T * (  1.19580846e-02 
			+ T * ( -5.26715675e-06 + T * (  1.09981875e-09
			+ T * -8.84016751e-14 ) ) ) );
		h[sC3H3] =  2.12893430e+02 * (
			T * (  6.14915291e+00 + T * (  4.67031583e-03
			+ T * ( -1.25018451e-06 + T * (  1.72539079e-10
			+ T * -9.21649988e-15 ) ) ) ) +  3.83854848e+04 );
		cp[sC3H3] =  2.12893430e+02 * (
			 6.14915291e+00 + T * (  9.34063166e-03 
			+ T * ( -3.75055354e-06 + T * (  6.90156316e-10
			+ T * -4.60824994e-14 ) ) ) );
		h[sA1C2H2XC8H7] =  8.06153041e+01 * (
			T * (  5.98044803e+00 + T * (  2.34715873e-02
			+ T * ( -8.91261580e-06 + T * (  1.82444488e-09
			+ T * -1.54218606e-13 ) ) ) ) +  4.32864172e+04 );
		cp[sA1C2H2XC8H7] =  8.06153041e+01 * (
			 5.98044803e+00 + T * (  4.69431747e-02 
			+ T * ( -2.67378474e-05 + T * (  7.29777950e-09
			+ T * -7.71093028e-13 ) ) ) );
		h[sA1C2H3XC8H8] =  7.98350361e+01 * (
			T * (  5.40554217e-01 + T * (  3.08651181e-02
			+ T * ( -1.24649102e-05 + T * (  2.67616467e-09
			+ T * -2.34609968e-13 ) ) ) ) +  1.50413170e+04 );
		cp[sA1C2H3XC8H8] =  7.98350361e+01 * (
			 5.40554217e-01 + T * (  6.17302362e-02 
			+ T * ( -3.73947305e-05 + T * (  1.07046587e-08
			+ T * -1.17304984e-12 ) ) ) );
		h[sA1C2HXC8H6] =  8.14109745e+01 * (
			T * (  5.81520488e+00 + T * (  2.20436466e-02
			+ T * ( -8.40179527e-06 + T * (  1.72568807e-09
			+ T * -1.46275782e-13 ) ) ) ) +  3.30271906e+04 );
		cp[sA1C2HXC8H6] =  8.14109745e+01 * (
			 5.81520488e+00 + T * (  4.40872933e-02 
			+ T * ( -2.52053858e-05 + T * (  6.90275228e-09
			+ T * -7.31378908e-13 ) ) ) );
		h[sA1C2HYXC8H5] =  8.22225079e+01 * (
			T * (  7.23812069e+00 + T * (  1.91906054e-02
			+ T * ( -7.29502437e-06 + T * (  1.49290312e-09
			+ T * -1.26070293e-13 ) ) ) ) +  6.49528135e+04 );
		cp[sA1C2HYXC8H5] =  8.22225079e+01 * (
			 7.23812069e+00 + T * (  3.83812109e-02 
			+ T * ( -2.18850731e-05 + T * (  5.97161247e-09
			+ T * -6.30351467e-13 ) ) ) );
		h[sA1C2H3YXC8H7] =  8.06153041e+01 * (
			T * (  3.90114779e+00 + T * (  2.57947010e-02
			+ T * ( -1.01693507e-05 + T * (  2.13977724e-09
			+ T * -1.84609351e-13 ) ) ) ) +  4.59935428e+04 );
		cp[sA1C2H3YXC8H7] =  8.06153041e+01 * (
			 3.90114779e+00 + T * (  5.15894020e-02 
			+ T * ( -3.05080522e-05 + T * (  8.55910896e-09
			+ T * -9.23046757e-13 ) ) ) );
		h[sA2XXC10H7] =  6.53869263e+01 * (
			T * (  3.22892303e+00 + T * (  3.15632243e-02
			+ T * ( -1.26860794e-05 + T * (  2.71135172e-09
			+ T * -2.36685024e-13 ) ) ) ) +  4.78400840e+04 );
		cp[sA2XXC10H7] =  6.53869263e+01 * (
			 3.22892303e+00 + T * (  6.31264486e-02 
			+ T * ( -3.80582381e-05 + T * (  1.08454069e-08
			+ T * -1.18342512e-12 ) ) ) );
		h[sA2XC10H8] =  6.48726632e+01 * (
			T * (  1.76826275e+00 + T * (  3.44571753e-02
			+ T * ( -1.38107392e-05 + T * (  2.94785772e-09
			+ T * -2.57194122e-13 ) ) ) ) +  1.45412795e+04 );
		cp[sA2XC10H8] =  6.48726632e+01 * (
			 1.76826275e+00 + T * (  6.89143506e-02 
			+ T * ( -4.14322176e-05 + T * (  1.17914309e-08
			+ T * -1.28597061e-12 ) ) ) );
		h[sC4H4] =  1.59670072e+02 * (
			T * (  4.97237210e+00 + T * (  9.65699520e-03
			+ T * ( -3.27065503e-06 + T * (  6.07512635e-10
			+ T * -4.74199476e-14 ) ) ) ) +  3.30561454e+04 );
		cp[sC4H4] =  1.59670072e+02 * (
			 4.97237210e+00 + T * (  1.93139904e-02 
			+ T * ( -9.81196508e-06 + T * (  2.43005054e-09
			+ T * -2.37099738e-13 ) ) ) );
		h[sA2YXC10H7] =  6.53869263e+01 * (
			T * (  3.29950506e+00 + T * (  3.15066683e-02
			+ T * ( -1.26586694e-05 + T * (  2.70451890e-09
			+ T * -2.36015394e-13 ) ) ) ) +  4.76658373e+04 );
		cp[sA2YXC10H7] =  6.53869263e+01 * (
			 3.29950506e+00 + T * (  6.30133365e-02 
			+ T * ( -3.79760083e-05 + T * (  1.08180756e-08
			+ T * -1.18007697e-12 ) ) ) );
		h[sA2C2H2AXC12H9] =  5.42739830e+01 * (
			T * (  8.53385239e+00 + T * (  3.43771399e-02
			+ T * ( -1.33916700e-05 + T * (  2.78702865e-09
			+ T * -2.38124486e-13 ) ) ) ) +  5.27583345e+04 );
		cp[sA2C2H2AXC12H9] =  5.42739830e+01 * (
			 8.53385239e+00 + T * (  6.87542797e-02 
			+ T * ( -4.01750101e-05 + T * (  1.11481146e-08
			+ T * -1.19062243e-12 ) ) ) );
		h[sA2C2H2BXC12H9] =  5.42739830e+01 * (
			T * (  7.59341132e+00 + T * (  3.51388460e-02
			+ T * ( -1.36912519e-05 + T * (  2.85023082e-09
			+ T * -2.43721726e-13 ) ) ) ) +  5.29378620e+04 );
		cp[sA2C2H2BXC12H9] =  5.42739830e+01 * (
			 7.59341132e+00 + T * (  7.02776920e-02 
			+ T * ( -4.10737558e-05 + T * (  1.14009233e-08
			+ T * -1.21860863e-12 ) ) ) );
		h[sA2C2HAXC12H8] =  5.46334700e+01 * (
			T * (  7.55690939e+00 + T * (  3.35536637e-02
			+ T * ( -1.30965256e-05 + T * (  2.73699948e-09
			+ T * -2.34921860e-13 ) ) ) ) +  3.91372224e+04 );
		cp[sA2C2HAXC12H8] =  5.46334700e+01 * (
			 7.55690939e+00 + T * (  6.71073273e-02 
			+ T * ( -3.92895768e-05 + T * (  1.09479979e-08
			+ T * -1.17460930e-12 ) ) ) );
		h[sA2C2HBXC12H8] =  5.46334700e+01 * (
			T * (  7.63899557e+00 + T * (  3.34667928e-02
			+ T * ( -1.30500037e-05 + T * (  2.72499707e-09
			+ T * -2.33731998e-13 ) ) ) ) +  3.92947046e+04 );
		cp[sA2C2HBXC12H8] =  5.46334700e+01 * (
			 7.63899557e+00 + T * (  6.69335855e-02 
			+ T * ( -3.91500110e-05 + T * (  1.08999883e-08
			+ T * -1.16865999e-12 ) ) ) );
		h[sA2C2HAYXC12H7] =  5.49977510e+01 * (
			T * (  8.88459555e+00 + T * (  3.07883585e-02
			+ T * ( -1.20329594e-05 + T * (  2.51481775e-09
			+ T * -2.15735058e-13 ) ) ) ) +  7.36259656e+04 );
		cp[sA2C2HAYXC12H7] =  5.49977510e+01 * (
			 8.88459555e+00 + T * (  6.15767170e-02 
			+ T * ( -3.60988783e-05 + T * (  1.00592710e-08
			+ T * -1.07867529e-12 ) ) ) );
		h[sA2C2HBYXC12H7] =  5.49977510e+01 * (
			T * (  8.88789581e+00 + T * (  3.07752080e-02
			+ T * ( -1.20249504e-05 + T * (  2.51275178e-09
			+ T * -2.15535646e-13 ) ) ) ) +  7.39803629e+04 );
		cp[sA2C2HBYXC12H7] =  5.49977510e+01 * (
			 8.88789581e+00 + T * (  6.15504161e-02 
			+ T * ( -3.60748511e-05 + T * (  1.00510071e-08
			+ T * -1.07767823e-12 ) ) ) );
		h[sA2R5XC12H8] =  5.46334700e+01 * (
			T * (  3.65432884e+00 + T * (  3.76323618e-02
			+ T * ( -1.51621650e-05 + T * (  3.24488352e-09
			+ T * -2.83461654e-13 ) ) ) ) +  2.65223472e+04 );
		cp[sA2R5XC12H8] =  5.46334700e+01 * (
			 3.65432884e+00 + T * (  7.52647236e-02 
			+ T * ( -4.54864951e-05 + T * (  1.29795341e-08
			+ T * -1.41730827e-12 ) ) ) );
		h[sA2R5XXC12H7] =  5.49977510e+01 * (
			T * (  4.90108932e+00 + T * (  3.49465809e-02
			+ T * ( -1.41408620e-05 + T * (  3.03365740e-09
			+ T * -2.65367420e-13 ) ) ) ) +  5.94391140e+04 );
		cp[sA2R5XXC12H7] =  5.49977510e+01 * (
			 4.90108932e+00 + T * (  6.98931618e-02 
			+ T * ( -4.24225860e-05 + T * (  1.21346296e-08
			+ T * -1.32683710e-12 ) ) ) );
		h[sA2R5C2H2XC14H9] =  4.69174774e+01 * (
			T * (  7.80126948e+00 + T * (  4.07307673e-02
			+ T * ( -1.62550125e-05 + T * (  3.45333605e-09
			+ T * -3.00000482e-13 ) ) ) ) +  6.50870635e+04 );
		cp[sA2R5C2H2XC14H9] =  4.69174774e+01 * (
			 7.80126948e+00 + T * (  8.14615345e-02 
			+ T * ( -4.87650376e-05 + T * (  1.38133442e-08
			+ T * -1.50000241e-12 ) ) ) );
		h[sA2R5C2HXC14H8] =  4.71858755e+01 * (
			T * (  9.29417050e+00 + T * (  3.68848112e-02
			+ T * ( -1.45277531e-05 + T * (  3.05383987e-09
			+ T * -2.63098636e-13 ) ) ) ) +  5.14110057e+04 );
		cp[sA2R5C2HXC14H8] =  4.71858755e+01 * (
			 9.29417050e+00 + T * (  7.37696223e-02 
			+ T * ( -4.35832594e-05 + T * (  1.22153595e-08
			+ T * -1.31549318e-12 ) ) ) );
		h[sA2R5C2HYXC14H7] =  4.74573620e+01 * (
			T * (  1.06117453e+01 + T * (  3.41337752e-02
			+ T * ( -1.34734130e-05 + T * (  2.83431000e-09
			+ T * -2.44192732e-13 ) ) ) ) +  8.51154522e+04 );
		cp[sA2R5C2HYXC14H7] =  4.74573620e+01 * (
			 1.06117453e+01 + T * (  6.82675505e-02 
			+ T * ( -4.04202391e-05 + T * (  1.13372400e-08
			+ T * -1.22096366e-12 ) ) ) );
		h[sP2XC12H10] =  5.39191958e+01 * (
			T * (  5.73686527e+00 + T * (  3.77329723e-02
			+ T * ( -1.46228290e-05 + T * (  3.04040273e-09
			+ T * -2.60024582e-13 ) ) ) ) +  1.66022411e+04 );
		cp[sP2XC12H10] =  5.39191958e+01 * (
			 5.73686527e+00 + T * (  7.54659445e-02 
			+ T * ( -4.38684869e-05 + T * (  1.21616109e-08
			+ T * -1.30012291e-12 ) ) ) );
		h[sP2XXC12H9] =  5.42739830e+01 * (
			T * (  3.97670430e+00 + T * (  3.89716035e-02
			+ T * ( -1.58160723e-05 + T * (  3.40940205e-09
			+ T * -2.99715664e-13 ) ) ) ) +  5.02964162e+04 );
		cp[sP2XXC12H9] =  5.42739830e+01 * (
			 3.97670430e+00 + T * (  7.79432070e-02 
			+ T * ( -4.74482168e-05 + T * (  1.36376082e-08
			+ T * -1.49857832e-12 ) ) ) );
		h[sA3XXC14H9] =  4.69174774e+01 * (
			T * (  4.71264594e+00 + T * (  4.33615001e-02
			+ T * ( -1.75337211e-05 + T * (  3.76148872e-09
			+ T * -3.29127952e-13 ) ) ) ) +  5.37939644e+04 );
		cp[sA3XXC14H9] =  4.69174774e+01 * (
			 4.71264594e+00 + T * (  8.67230002e-02 
			+ T * ( -5.26011632e-05 + T * (  1.50459549e-08
			+ T * -1.64563976e-12 ) ) ) );
		h[sA3XC14H10] =  4.66521154e+01 * (
			T * (  3.38725839e+00 + T * (  4.60942802e-02
			+ T * ( -1.85762331e-05 + T * (  3.97809798e-09
			+ T * -3.47768516e-13 ) ) ) ) +  1.91061794e+04 );
		cp[sA3XC14H10] =  4.66521154e+01 * (
			 3.38725839e+00 + T * (  9.21885604e-02 
			+ T * ( -5.57286994e-05 + T * (  1.59123919e-08
			+ T * -1.73884258e-12 ) ) ) );
		h[sA3YXC14H9] =  4.69174774e+01 * (
			T * (  4.71264594e+00 + T * (  4.33615001e-02
			+ T * ( -1.75337211e-05 + T * (  3.76148872e-09
			+ T * -3.29127952e-13 ) ) ) ) +  5.37939644e+04 );
		cp[sA3YXC14H9] =  4.69174774e+01 * (
			 4.71264594e+00 + T * (  8.67230002e-02 
			+ T * ( -5.26011632e-05 + T * (  1.50459549e-08
			+ T * -1.64563976e-12 ) ) ) );
		h[sA3R5XXC16H9] =  4.13171861e+01 * (
			T * (  6.47777347e+00 + T * (  4.66668858e-02
			+ T * ( -1.89542299e-05 + T * (  4.07616817e-09
			+ T * -3.57133286e-13 ) ) ) ) +  6.47489852e+04 );
		cp[sA3R5XXC16H9] =  4.13171861e+01 * (
			 6.47777347e+00 + T * (  9.33337716e-02 
			+ T * ( -5.68626896e-05 + T * (  1.63046727e-08
			+ T * -1.78566643e-12 ) ) ) );
		h[sA3R5XC16H10] =  4.11112540e+01 * (
			T * (  5.07024731e+00 + T * (  4.94958652e-02
			+ T * ( -2.00466753e-05 + T * (  4.30511463e-09
			+ T * -3.76952302e-13 ) ) ) ) +  3.06528296e+04 );
		cp[sA3R5XC16H10] =  4.11112540e+01 * (
			 5.07024731e+00 + T * (  9.89917305e-02 
			+ T * ( -6.01400259e-05 + T * (  1.72204585e-08
			+ T * -1.88476151e-12 ) ) ) );
		h[sA4XC16H10] =  4.11112540e+01 * (
			T * (  4.54060055e+00 + T * (  4.99057604e-02
			+ T * ( -2.02125434e-05 + T * (  4.33937548e-09
			+ T * -3.79804636e-13 ) ) ) ) +  2.12755890e+04 );
		cp[sA4XC16H10] =  4.11112540e+01 * (
			 4.54060055e+00 + T * (  9.98115207e-02 
			+ T * ( -6.06376301e-05 + T * (  1.73575019e-08
			+ T * -1.89902318e-12 ) ) ) );
		h[sA4XXC16H9] =  4.13171861e+01 * (
			T * (  5.85098138e+00 + T * (  4.71615553e-02
			+ T * ( -1.91632717e-05 + T * (  4.12147133e-09
			+ T * -3.61083018e-13 ) ) ) ) +  5.89572568e+04 );
		cp[sA4XXC16H9] =  4.13171861e+01 * (
			 5.85098138e+00 + T * (  9.43231105e-02 
			+ T * ( -5.74898152e-05 + T * (  1.64858853e-08
			+ T * -1.80541509e-12 ) ) ) );
		h[sA4R5XC18H10] =  3.67468399e+01 * (
			T * (  6.20190827e+00 + T * (  5.33281375e-02
			+ T * ( -2.16943100e-05 + T * (  4.66940325e-09
			+ T * -4.09294700e-13 ) ) ) ) +  3.34439422e+04 );
		cp[sA4R5XC18H10] =  3.67468399e+01 * (
			 6.20190827e+00 + T * (  1.06656275e-01 
			+ T * ( -6.50829301e-05 + T * (  1.86776130e-08
			+ T * -2.04647350e-12 ) ) ) );
		h[sFLTNXC16H10] =  4.11112540e+01 * (
			T * (  4.54792547e+00 + T * (  4.99994870e-02
			+ T * ( -2.02900657e-05 + T * (  4.36322555e-09
			+ T * -3.82399156e-13 ) ) ) ) +  2.54780117e+04 );
		cp[sFLTNXC16H10] =  4.11112540e+01 * (
			 4.54792547e+00 + T * (  9.99989740e-02 
			+ T * ( -6.08701972e-05 + T * (  1.74529022e-08
			+ T * -1.91199578e-12 ) ) ) );
		h[sC5H6] =  1.25788072e+02 * (
			T * (  2.30537462e-01 + T * (  2.04785913e-02
			+ T * ( -8.05296527e-06 + T * (  1.69940870e-09
			+ T * -1.47274884e-13 ) ) ) ) +  1.43779465e+04 );
		cp[sC5H6] =  1.25788072e+02 * (
			 2.30537462e-01 + T * (  4.09571826e-02 
			+ T * ( -2.41588958e-05 + T * (  6.79763480e-09
			+ T * -7.36374421e-13 ) ) ) );
		h[sC5H5] =  1.27736058e+02 * (
			T * (  4.21464919e+00 + T * (  1.35917364e-02
			+ T * ( -4.43910697e-06 + T * (  7.72450297e-10
			+ T * -5.55759746e-14 ) ) ) ) +  2.88952416e+04 );
		cp[sC5H5] =  1.27736058e+02 * (
			 4.21464919e+00 + T * (  2.71834728e-02 
			+ T * ( -1.33173209e-05 + T * (  3.08980119e-09
			+ T * -2.77879873e-13 ) ) ) );
		h[sSXC3H5] =  2.02443146e+02 * (
			T * (  2.02509360e+00 + T * (  1.17756625e-02
			+ T * ( -4.27515187e-06 + T * (  8.48948055e-10
			+ T * -7.03589448e-14 ) ) ) ) +  3.11812042e+04 );
		cp[sSXC3H5] =  2.02443146e+02 * (
			 2.02509360e+00 + T * (  2.35513249e-02 
			+ T * ( -1.28254556e-05 + T * (  3.39579222e-09
			+ T * -3.51794724e-13 ) ) ) );
		h[sCH4] =  5.18285750e+02 * (
			T * (  1.65326226e+00 + T * (  5.01315495e-03
			+ T * ( -1.10553746e-06 + T * (  1.34120785e-10
			+ T * -6.29393516e-15 ) ) ) ) -1.00095936e+04 );
		cp[sCH4] =  5.18285750e+02 * (
			 1.65326226e+00 + T * (  1.00263099e-02 
			+ T * ( -3.31661238e-06 + T * (  5.36483138e-10
			+ T * -3.14696758e-14 ) ) ) );
		h[sC4H6] =  1.53718755e+02 * (
			T * ( -8.99531092e+00 + T * (  3.00857534e-02
			+ T * ( -1.40019253e-05 + T * (  3.33325140e-09
			+ T * -3.14847380e-13 ) ) ) ) +  1.49296107e+04 );
		cp[sC4H6] =  1.53718755e+02 * (
			-8.99531092e+00 + T * (  6.01715069e-02 
			+ T * ( -4.20057758e-05 + T * (  1.33330056e-08
			+ T * -1.57423690e-12 ) ) ) );
		h[sTXC5H5O] =  1.02532248e+02 * (
			T * (  1.26065350e+01 + T * (  8.37353350e-03
			+ T * ( -2.03658623e-06 + T * (  2.49186440e-10
			+ T * -1.20223668e-14 ) ) ) ) +  1.41146570e+03 );
		cp[sTXC5H5O] =  1.02532248e+02 * (
			 1.26065350e+01 + T * (  1.67470670e-02 
			+ T * ( -6.10975870e-06 + T * (  9.96745760e-10
			+ T * -6.01118340e-14 ) ) ) );
		h[sC5H4O] =  1.03822832e+02 * (
			T * (  4.25344911e+00 + T * (  1.56819909e-02
			+ T * ( -6.09546950e-06 + T * (  1.27102091e-09
			+ T * -1.08969098e-13 ) ) ) ) +  3.87579835e+03 );
		cp[sC5H4O] =  1.03822832e+02 * (
			 4.25344911e+00 + T * (  3.13639818e-02 
			+ T * ( -1.82864085e-05 + T * (  5.08408365e-09
			+ T * -5.44845492e-13 ) ) ) );
		h[sSXC5H5O] =  1.02532248e+02 * (
			T * (  8.54053120e+00 + T * (  1.14947550e-02
			+ T * ( -3.18125210e-06 + T * (  4.26540300e-10
			+ T * -1.94918720e-14 ) ) ) ) +  2.22636990e+04 );
		cp[sSXC5H5O] =  1.02532248e+02 * (
			 8.54053120e+00 + T * (  2.29895100e-02 
			+ T * ( -9.54375630e-06 + T * (  1.70616120e-09
			+ T * -9.74593600e-14 ) ) ) );
		h[sC9H8] =  7.15803158e+01 * (
			T * (  1.15459802e+00 + T * (  3.27112098e-02
			+ T * ( -1.30835036e-05 + T * (  2.78922102e-09
			+ T * -2.43185350e-13 ) ) ) ) +  1.68166108e+04 );
		cp[sC9H8] =  7.15803158e+01 * (
			 1.15459802e+00 + T * (  6.54224196e-02 
			+ T * ( -3.92505107e-05 + T * (  1.11568841e-08
			+ T * -1.21592675e-12 ) ) ) );
		h[sC9H7] =  7.22069373e+01 * (
			T * (  3.65597547e+00 + T * (  2.87404231e-02
			+ T * ( -1.14290200e-05 + T * (  2.42569698e-09
			+ T * -2.10772824e-13 ) ) ) ) +  3.06843457e+04 );
		cp[sC9H7] =  7.22069373e+01 * (
			 3.65597547e+00 + T * (  5.74808463e-02 
			+ T * ( -3.42870600e-05 + T * (  9.70278793e-09
			+ T * -1.05386412e-12 ) ) ) );
		h[sA1CH2XC7H7] =  9.12400413e+01 * (
			T * (  3.30049696e+00 + T * (  2.40027670e-02
			+ T * ( -9.28143407e-06 + T * (  1.93092839e-09
			+ T * -1.65430827e-13 ) ) ) ) +  2.17498572e+04 );
		cp[sA1CH2XC7H7] =  9.12400413e+01 * (
			 3.30049696e+00 + T * (  4.80055340e-02 
			+ T * ( -2.78443022e-05 + T * (  7.72371356e-09
			+ T * -8.27154136e-13 ) ) ) );
		h[sC9H6O] =  6.38886413e+01 * (
			T * (  4.65659248e+00 + T * (  2.85027911e-02
			+ T * ( -1.14391400e-05 + T * (  2.44044361e-09
			+ T * -2.12668074e-13 ) ) ) ) +  4.57857140e+03 );
		cp[sC9H6O] =  6.38886413e+01 * (
			 4.65659248e+00 + T * (  5.70055822e-02 
			+ T * ( -3.43174199e-05 + T * (  9.76177442e-09
			+ T * -1.06334037e-12 ) ) ) );
		h[sOXC6H4] =  1.09266940e+02 * (
			T * (  2.98618725e+00 + T * (  1.68818922e-02
			+ T * ( -6.67461303e-06 + T * (  1.40963421e-09
			+ T * -1.22000829e-13 ) ) ) ) +  5.12231321e+04 );
		cp[sOXC6H4] =  1.09266940e+02 * (
			 2.98618725e+00 + T * (  3.37637843e-02 
			+ T * ( -2.00238391e-05 + T * (  5.63853682e-09
			+ T * -6.10004145e-13 ) ) ) );
		h[sA1CH3XC7H8] =  9.02418217e+01 * (
			T * ( -1.01117220e+00 + T * (  2.92650956e-02
			+ T * ( -1.15865023e-05 + T * (  2.45545248e-09
			+ T * -2.13361740e-13 ) ) ) ) +  3.99363395e+03 );
		cp[sA1CH3XC7H8] =  9.02418217e+01 * (
			-1.01117220e+00 + T * (  5.85301912e-02 
			+ T * ( -3.47595069e-05 + T * (  9.82180993e-09
			+ T * -1.06680870e-12 ) ) ) );
		h[sA1OHXC6H6O] =  8.83489183e+01 * (
			T * (  9.33151850e-01 + T * (  2.53298752e-02
			+ T * ( -1.05872946e-05 + T * (  2.34389723e-09
			+ T * -2.10588162e-13 ) ) ) ) -1.37575260e+04 );
		cp[sA1OHXC6H6O] =  8.83489183e+01 * (
			 9.33151850e-01 + T * (  5.06597504e-02 
			+ T * ( -3.17618838e-05 + T * (  9.37558892e-09
			+ T * -1.05294081e-12 ) ) ) );
		h[sHOA1CH3XC7H8O] =  7.68892300e+01 * (
			T * ( -6.71722270e-02 + T * (  3.15933480e-02
			+ T * ( -1.30255983e-05 + T * (  2.85416352e-09
			+ T * -2.54530110e-13 ) ) ) ) -1.78554350e+04 );
		cp[sHOA1CH3XC7H8O] =  7.68892300e+01 * (
			-6.71722270e-02 + T * (  6.31866960e-02 
			+ T * ( -3.90767950e-05 + T * (  1.14166541e-08
			+ T * -1.27265055e-12 ) ) ) );
		h[sOA1CH3XC7H7O] =  7.76127177e+01 * (
			T * (  4.14521668e+00 + T * (  2.51392063e-02
			+ T * ( -9.69283660e-06 + T * (  2.00710898e-09
			+ T * -1.71158850e-13 ) ) ) ) -1.22524065e+03 );
		cp[sOA1CH3XC7H7O] =  7.76127177e+01 * (
			 4.14521668e+00 + T * (  5.02784126e-02 
			+ T * ( -2.90785098e-05 + T * (  8.02843592e-09
			+ T * -8.55794251e-13 ) ) ) );
		h[sA1CH2OXC7H7O] =  7.76127177e+01 * (
			T * (  2.07930551e+00 + T * (  2.76359457e-02
			+ T * ( -1.10434375e-05 + T * (  2.35169868e-09
			+ T * -2.04829704e-13 ) ) ) ) +  1.17599254e+04 );
		cp[sA1CH2OXC7H7O] =  7.76127177e+01 * (
			 2.07930551e+00 + T * (  5.52718914e-02 
			+ T * ( -3.31303125e-05 + T * (  9.40679473e-09
			+ T * -1.02414852e-12 ) ) ) );
		h[sA1CH2OHXC7H8O] =  7.68892300e+01 * (
			T * (  1.51623145e+01 + T * (  1.34685184e-02
			+ T * ( -3.26763210e-06 + T * (  3.98205190e-10
			+ T * -1.91421567e-14 ) ) ) ) -1.88226234e+04 );
		cp[sA1CH2OHXC7H8O] =  7.68892300e+01 * (
			 1.51623145e+01 + T * (  2.69370369e-02 
			+ T * ( -9.80289631e-06 + T * (  1.59282076e-09
			+ T * -9.57107837e-14 ) ) ) );
		h[sA1CHOXC7H6O] =  7.83499501e+01 * (
			T * (  1.87355756e+00 + T * (  2.63115775e-02
			+ T * ( -1.05881654e-05 + T * (  2.26600767e-09
			+ T * -1.98061225e-13 ) ) ) ) -7.23603865e+03 );
		cp[sA1CHOXC7H6O] =  7.83499501e+01 * (
			 1.87355756e+00 + T * (  5.26231551e-02 
			+ T * ( -3.17644962e-05 + T * (  9.06403069e-09
			+ T * -9.90306123e-13 ) ) ) );
		h[sA1OXC6H5O] =  8.93054780e+01 * (
			T * (  3.39256520e+00 + T * (  2.08689845e-02
			+ T * ( -8.32789277e-06 + T * (  1.77206751e-09
			+ T * -1.54268601e-13 ) ) ) ) +  3.60336039e+03 );
		cp[sA1OXC6H5O] =  8.93054780e+01 * (
			 3.39256520e+00 + T * (  4.17379690e-02 
			+ T * ( -2.49836783e-05 + T * (  7.08827005e-09
			+ T * -7.71343006e-13 ) ) ) );
		h[sA1CH3YXC7H7] =  9.12400413e+01 * (
			T * (  5.19780822e-01 + T * (  2.63109877e-02
			+ T * ( -1.04327811e-05 + T * (  2.21204344e-09
			+ T * -1.92250440e-13 ) ) ) ) +  3.45681510e+04 );
		cp[sA1CH3YXC7H7] =  9.12400413e+01 * (
			 5.19780822e-01 + T * (  5.26219754e-02 
			+ T * ( -3.12983433e-05 + T * (  8.84817377e-09
			+ T * -9.61252202e-13 ) ) ) );
		h[sPXC3H4] =  2.07536818e+02 * (
			T * (  2.81460543e+00 + T * (  9.27622480e-03
			+ T * ( -3.18342256e-06 + T * (  5.99878425e-10
			+ T * -4.74970514e-14 ) ) ) ) +  2.07010771e+04 );
		cp[sPXC3H4] =  2.07536818e+02 * (
			 2.81460543e+00 + T * (  1.85524496e-02 
			+ T * ( -9.55026768e-06 + T * (  2.39951370e-09
			+ T * -2.37485257e-13 ) ) ) );
		h[sA1C2H4XC8H9] =  7.90697276e+01 * (
			T * (  1.61326962e+01 + T * (  1.41452137e-02
			+ T * ( -3.39339587e-06 + T * (  4.10441592e-10
			+ T * -1.96275066e-14 ) ) ) ) +  2.08791061e+04 );
		cp[sA1C2H4XC8H9] =  7.90697276e+01 * (
			 1.61326962e+01 + T * (  2.82904273e-02 
			+ T * ( -1.01801876e-05 + T * (  1.64176637e-09
			+ T * -9.81375329e-14 ) ) ) );
		h[sA1C2H5XC8H10] =  7.83189525e+01 * (
			T * (  1.56901336e+01 + T * (  1.61831537e-02
			+ T * ( -3.89548593e-06 + T * (  4.72473905e-10
			+ T * -2.26403582e-14 ) ) ) ) -4.38669907e+03 );
		cp[sA1C2H5XC8H10] =  7.83189525e+01 * (
			 1.56901336e+01 + T * (  3.23663075e-02 
			+ T * ( -1.16864578e-05 + T * (  1.88989562e-09
			+ T * -1.13201791e-13 ) ) ) );
		h[sC8H9O2] =  6.06213544e+01 * (
			T * (  1.74946727e+00 + T * (  3.66934387e-02
			+ T * ( -1.47880348e-05 + T * (  3.17289810e-09
			+ T * -2.78067684e-13 ) ) ) ) +  1.71120874e+04 );
		cp[sC8H9O2] =  6.06213544e+01 * (
			 1.74946727e+00 + T * (  7.33868774e-02 
			+ T * ( -4.43641045e-05 + T * (  1.26915924e-08
			+ T * -1.39033842e-12 ) ) ) );
		h[sC8H8OOH] =  6.06213544e+01 * (
			T * (  1.74946727e+00 + T * (  3.66934387e-02
			+ T * ( -1.47880348e-05 + T * (  3.17289810e-09
			+ T * -2.78067684e-13 ) ) ) ) +  1.71120874e+04 );
		cp[sC8H8OOH] =  6.06213544e+01 * (
			 1.74946727e+00 + T * (  7.33868774e-02 
			+ T * ( -4.43641045e-05 + T * (  1.26915924e-08
			+ T * -1.39033842e-12 ) ) ) );
		h[sOC8H7OOH] =  5.46478336e+01 * (
			T * (  5.90515216e+00 + T * (  3.35364725e-02
			+ T * ( -1.34427129e-05 + T * (  2.87414260e-09
			+ T * -2.51281962e-13 ) ) ) ) -2.17516302e+04 );
		cp[sOC8H7OOH] =  5.46478336e+01 * (
			 5.90515216e+00 + T * (  6.70729450e-02 
			+ T * ( -4.03281386e-05 + T * (  1.14965704e-08
			+ T * -1.25640981e-12 ) ) ) );
		h[sH2C2] =  3.19340144e+02 * (
			T * (  4.27803400e+00 + T * (  2.37814020e-03
			+ T * ( -5.43366967e-07 + T * (  6.36570150e-11
			+ T * -2.97727580e-15 ) ) ) ) +  4.83166880e+04 );
		cp[sH2C2] =  3.19340144e+02 * (
			 4.27803400e+00 + T * (  4.75628040e-03 
			+ T * ( -1.63010090e-06 + T * (  2.54628060e-10
			+ T * -1.48863790e-14 ) ) ) );
		h[sA1CH3CH3XC8H10] =  7.83189525e+01 * (
			T * ( -1.95577967e+00 + T * (  3.54776361e-02
			+ T * ( -1.39991811e-05 + T * (  2.95930647e-09
			+ T * -2.56679432e-13 ) ) ) ) -6.91883225e+01 );
		cp[sA1CH3CH3XC8H10] =  7.83189525e+01 * (
			-1.95577967e+00 + T * (  7.09552723e-02 
			+ T * ( -4.19975432e-05 + T * (  1.18372259e-08
			+ T * -1.28339716e-12 ) ) ) );
		h[sA1CH3CH2XC8H9] =  7.90697276e+01 * (
			T * (  2.36833258e+00 + T * (  3.01967514e-02
			+ T * ( -1.16835702e-05 + T * (  2.43202571e-09
			+ T * -2.08472266e-13 ) ) ) ) +  1.76344521e+04 );
		cp[sA1CH3CH2XC8H9] =  7.90697276e+01 * (
			 2.36833258e+00 + T * (  6.03935028e-02 
			+ T * ( -3.50507105e-05 + T * (  9.72810282e-09
			+ T * -1.04236133e-12 ) ) ) );
		h[sA1CH3CHOXC8H8O] =  6.92031229e+01 * (
			T * (  9.22059379e-01 + T * (  3.25311564e-02
			+ T * ( -1.30067510e-05 + T * (  2.77167152e-09
			+ T * -2.41572622e-13 ) ) ) ) -1.13895901e+04 );
		cp[sA1CH3CHOXC8H8O] =  6.92031229e+01 * (
			 9.22059379e-01 + T * (  6.50623128e-02 
			+ T * ( -3.90202529e-05 + T * (  1.10866861e-08
			+ T * -1.20786311e-12 ) ) ) );
		h[sA2CH3XC11H10] =  5.84734510e+01 * (
			T * (  1.43553166e+00 + T * (  4.06226654e-02
			+ T * ( -1.61966967e-05 + T * (  3.44489735e-09
			+ T * -2.99849226e-13 ) ) ) ) +  1.04031044e+04 );
		cp[sA2CH3XC11H10] =  5.84734510e+01 * (
			 1.43553166e+00 + T * (  8.12453307e-02 
			+ T * ( -4.85900900e-05 + T * (  1.37795894e-08
			+ T * -1.49924613e-12 ) ) ) );
		h[sA1CHOCH2XC8H7O] =  6.97886449e+01 * (
			T * (  5.25955102e+00 + T * (  2.72439788e-02
			+ T * ( -1.06895200e-05 + T * (  2.24419442e-09
			+ T * -1.93360743e-13 ) ) ) ) +  6.61511721e+03 );
		cp[sA1CHOCH2XC8H7O] =  6.97886449e+01 * (
			 5.25955102e+00 + T * (  5.44879576e-02 
			+ T * ( -3.20685601e-05 + T * (  8.97677768e-09
			+ T * -9.66803716e-13 ) ) ) );
		h[sA1CHOCHOXC8H6O2] =  6.19881009e+01 * (
			T * (  3.92023646e+00 + T * (  2.95004798e-02
			+ T * ( -1.19741621e-05 + T * (  2.57441813e-09
			+ T * -2.25560976e-13 ) ) ) ) -2.20802686e+04 );
		cp[sA1CHOCHOXC8H6O2] =  6.19881009e+01 * (
			 3.92023646e+00 + T * (  5.90009596e-02 
			+ T * ( -3.59224864e-05 + T * (  1.02976725e-08
			+ T * -1.12780488e-12 ) ) ) );
		h[sA2OHXC10H8O] =  5.76727893e+01 * (
			T * (  2.08930252e+01 + T * (  1.55280033e-02
			+ T * ( -3.81358540e-06 + T * (  4.69682165e-10
			+ T * -2.27647762e-14 ) ) ) ) -1.35886443e+04 );
		cp[sA2OHXC10H8O] =  5.76727893e+01 * (
			 2.08930252e+01 + T * (  3.10560066e-02 
			+ T * ( -1.14407562e-05 + T * (  1.87872866e-09
			+ T * -1.13823881e-13 ) ) ) );
		h[sA2CH2XC11H9] =  5.88909351e+01 * (
			T * (  4.97463689e+00 + T * (  3.55734585e-02
			+ T * ( -1.40082277e-05 + T * (  2.95063847e-09
			+ T * -2.54923676e-13 ) ) ) ) +  2.96267836e+04 );
		cp[sA2CH2XC11H9] =  5.88909351e+01 * (
			 4.97463689e+00 + T * (  7.11469171e-02 
			+ T * ( -4.20246831e-05 + T * (  1.18025539e-08
			+ T * -1.27461838e-12 ) ) ) );
		h[sA2CH2OXC11H9O] =  5.28962604e+01 * (
			T * (  3.88859655e+00 + T * (  3.90441740e-02
			+ T * ( -1.56825396e-05 + T * (  3.34906190e-09
			+ T * -2.92138322e-13 ) ) ) ) +  2.16424127e+04 );
		cp[sA2CH2OXC11H9O] =  5.28962604e+01 * (
			 3.88859655e+00 + T * (  7.80883481e-02 
			+ T * ( -4.70476187e-05 + T * (  1.33962476e-08
			+ T * -1.46069161e-12 ) ) ) );
		h[sA2CHOXC11H8O] =  5.32376708e+01 * (
			T * ( -5.55792190e-01 + T * (  4.36874222e-02
			+ T * ( -1.87272784e-05 + T * (  4.18083215e-09
			+ T * -3.76185702e-13 ) ) ) ) +  1.03455228e+03 );
		cp[sA2CHOXC11H8O] =  5.32376708e+01 * (
			-5.55792190e-01 + T * (  8.73748443e-02 
			+ T * ( -5.61818352e-05 + T * (  1.67233286e-08
			+ T * -1.88092851e-12 ) ) ) );
		h[sA2OXC10H7O] =  5.80788790e+01 * (
			T * (  2.10591364e+01 + T * (  1.41281535e-02
			+ T * ( -3.44428953e-06 + T * (  4.22167585e-10
			+ T * -2.03949534e-14 ) ) ) ) +  4.09143507e+03 );
		cp[sA2OXC10H7O] =  5.80788790e+01 * (
			 2.10591364e+01 + T * (  2.82563070e-02 
			+ T * ( -1.03328686e-05 + T * (  1.68867034e-09
			+ T * -1.01974767e-13 ) ) ) );
		h[sC4H2] =  1.66100767e+02 * (
			T * (  9.75839793e+00 + T * (  1.89436611e-03
			+ T * (  1.02047338e-07 + T * ( -1.58413756e-10
			+ T *  2.25860064e-14 ) ) ) ) +  5.22698696e+04 );
		cp[sC4H2] =  1.66100767e+02 * (
			 9.75839793e+00 + T * (  3.78873223e-03 
			+ T * (  3.06142015e-07 + T * ( -6.33655024e-10
			+ T *  1.12930032e-13 ) ) ) );
		h[sOC6H4O] =  7.69191059e+01 * (
			T * (  5.70290193e+00 + T * (  1.92522500e-02
			+ T * ( -7.62482393e-06 + T * (  1.60822969e-09
			+ T * -1.38851568e-13 ) ) ) ) -1.40768967e+04 );
		cp[sOC6H4O] =  7.69191059e+01 * (
			 5.70290193e+00 + T * (  3.85045001e-02 
			+ T * ( -2.28744718e-05 + T * (  6.43291874e-09
			+ T * -6.94257842e-13 ) ) ) );
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
			+ T * -1.95995140e-15 ) ) ) ) +  1.44072920e+04 );
		cp[sHCN] =  3.07619506e+02 * (
			 3.80223920e+00 + T * (  3.14642280e-03 
			+ T * ( -1.06321850e-06 + T * (  1.66197570e-10
			+ T * -9.79975700e-15 ) ) ) );
		h[sHOCN] =  1.93230919e+02 * (
			T * (  5.89784885e+00 + T * (  1.58394696e-03
			+ T * ( -3.72670213e-07 + T * (  4.43107860e-11
			+ T * -2.08678354e-15 ) ) ) ) -3.70653331e+03 );
		cp[sHOCN] =  1.93230919e+02 * (
			 5.89784885e+00 + T * (  3.16789393e-03 
			+ T * ( -1.11801064e-06 + T * (  1.77243144e-10
			+ T * -1.04339177e-14 ) ) ) );
		h[sHNCO] =  1.93230919e+02 * (
			T * (  6.22395134e+00 + T * (  1.58932002e-03
			+ T * ( -3.64595850e-07 + T * (  4.26837908e-11
			+ T * -1.99004391e-15 ) ) ) ) -1.66599344e+04 );
		cp[sHNCO] =  1.93230919e+02 * (
			 6.22395134e+00 + T * (  3.17864004e-03 
			+ T * ( -1.09378755e-06 + T * (  1.70735163e-10
			+ T * -9.95021955e-15 ) ) ) );
		h[sH2CN] =  2.96559424e+02 * (
			T * (  5.20970300e+00 + T * (  1.48464555e-03
			+ T * ( -9.51863033e-08 + T * ( -4.08887500e-11
			+ T *  6.08651780e-15 ) ) ) ) +  2.76771090e+04 );
		cp[sH2CN] =  2.96559424e+02 * (
			 5.20970300e+00 + T * (  2.96929110e-03 
			+ T * ( -2.85558910e-07 + T * ( -1.63555000e-10
			+ T *  3.04325890e-14 ) ) ) );
		h[sC] =  6.92284763e+02 * (
			T * (  2.49266888e+00 + T * (  2.39944642e-05
			+ T * ( -2.41445007e-08 + T * (  9.35727573e-12
			+ T * -9.74555786e-16 ) ) ) ) +  8.54512953e+04 );
		cp[sC] =  6.92284763e+02 * (
			 2.49266888e+00 + T * (  4.79889284e-05 
			+ T * ( -7.24335020e-08 + T * (  3.74291029e-11
			+ T * -4.87277893e-15 ) ) ) );
		h[sHCNN] =  2.02601004e+02 * (
			T * (  5.89463620e+00 + T * (  1.99479795e-03
			+ T * ( -5.32746000e-07 + T * (  7.31234875e-11
			+ T * -4.01893720e-15 ) ) ) ) +  5.34529410e+04 );
		cp[sHCNN] =  2.02601004e+02 * (
			 5.89463620e+00 + T * (  3.98959590e-03 
			+ T * ( -1.59823800e-06 + T * (  2.92493950e-10
			+ T * -2.00946860e-14 ) ) ) );
		h[sHCNO] =  1.93230919e+02 * (
			T * (  6.59860456e+00 + T * (  1.51389313e-03
			+ T * ( -3.59014487e-07 + T * (  4.29166320e-11
			+ T * -2.02878782e-15 ) ) ) ) +  1.79661339e+04 );
		cp[sHCNO] =  1.93230919e+02 * (
			 6.59860456e+00 + T * (  3.02778626e-03 
			+ T * ( -1.07704346e-06 + T * (  1.71666528e-10
			+ T * -1.01439391e-14 ) ) ) );
		h[sNH3] =  4.88102618e+02 * (
			T * (  2.63445210e+00 + T * (  2.83312800e-03
			+ T * ( -5.75955867e-07 + T * (  5.96679025e-11
			+ T * -2.51575720e-15 ) ) ) ) -6.54469580e+03 );
		cp[sNH3] =  4.88102618e+02 * (
			 2.63445210e+00 + T * (  5.66625600e-03 
			+ T * ( -1.72786760e-06 + T * (  2.38671610e-10
			+ T * -1.25787860e-14 ) ) ) );
	}
	else if (T >= 299.999999 ) {
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
		h[sO] =  5.19646250e+02 * (
			T * (  2.94642900e+00 + T * ( -8.19083000e-04
			+ T * (  8.07010667e-07 + T * ( -4.00710750e-10
			+ T *  7.78139200e-14 ) ) ) ) +  2.91476400e+04 );
		cp[sO] =  5.19646250e+02 * (
			 2.94642900e+00 + T * ( -1.63816600e-03 
			+ T * (  2.42103200e-06 + T * ( -1.60284300e-09
			+ T *  3.89069600e-13 ) ) ) );
		h[sH2] =  4.12417659e+03 * (
			T * (  3.29812400e+00 + T * (  4.12472100e-04
			+ T * ( -2.71433833e-07 + T * ( -2.36885850e-11
			+ T *  8.26974400e-14 ) ) ) ) -1.01252100e+03 );
		cp[sH2] =  4.12417659e+03 * (
			 3.29812400e+00 + T * (  8.24944200e-04 
			+ T * ( -8.14301500e-07 + T * ( -9.47543400e-11
			+ T *  4.13487200e-13 ) ) ) );
		h[sH] =  8.24835317e+03 * (
			T * (  2.50000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  2.54716300e+04 );
		cp[sH] =  8.24835317e+03 * (
			 2.50000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sOH] =  4.88848777e+02 * (
			T * (  4.12530561e+00 + T * ( -1.61272470e-03
			+ T * (  2.17588230e-06 + T * ( -1.44963411e-09
			+ T *  4.12474758e-13 ) ) ) ) +  3.34630913e+03 );
		cp[sOH] =  4.88848777e+02 * (
			 4.12530561e+00 + T * ( -3.22544939e-03 
			+ T * (  6.52764691e-06 + T * ( -5.79853643e-09
			+ T *  2.06237379e-12 ) ) ) );
		h[sH2O] =  4.61497558e+02 * (
			T * (  3.38684200e+00 + T * (  1.73749100e-03
			+ T * ( -2.11823200e-06 + T * (  1.74214525e-09
			+ T * -5.01317600e-13 ) ) ) ) -3.02081100e+04 );
		cp[sH2O] =  4.61497558e+02 * (
			 3.38684200e+00 + T * (  3.47498200e-03 
			+ T * ( -6.35469600e-06 + T * (  6.96858100e-09
			+ T * -2.50658800e-12 ) ) ) );
		h[sO2] =  2.59823125e+02 * (
			T * (  3.21293600e+00 + T * (  5.63743000e-04
			+ T * ( -1.91871667e-07 + T * (  3.28469250e-10
			+ T * -1.75371080e-13 ) ) ) ) -1.00524900e+03 );
		cp[sO2] =  2.59823125e+02 * (
			 3.21293600e+00 + T * (  1.12748600e-03 
			+ T * ( -5.75615000e-07 + T * (  1.31387700e-09
			+ T * -8.76855400e-13 ) ) ) );
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
		h[sCH] =  6.38680289e+02 * (
			T * (  3.48981665e+00 + T * (  1.61917771e-04
			+ T * ( -5.62996883e-07 + T * (  7.90543317e-10
			+ T * -2.81218134e-13 ) ) ) ) +  7.07972934e+04 );
		cp[sCH] =  6.38680289e+02 * (
			 3.48981665e+00 + T * (  3.23835541e-04 
			+ T * ( -1.68899065e-06 + T * (  3.16217327e-09
			+ T * -1.40609067e-12 ) ) ) );
		h[sTXCH2] =  5.92780550e+02 * (
			T * (  3.76267867e+00 + T * (  4.84436072e-04
			+ T * (  9.31632803e-07 + T * ( -9.62727883e-10
			+ T *  3.37483438e-13 ) ) ) ) +  4.60040401e+04 );
		cp[sTXCH2] =  5.92780550e+02 * (
			 3.76267867e+00 + T * (  9.68872143e-04 
			+ T * (  2.79489841e-06 + T * ( -3.85091153e-09
			+ T *  1.68741719e-12 ) ) ) );
		h[sCH2O] =  2.76904683e+02 * (
			T * (  4.79372315e+00 + T * ( -4.95416684e-03
			+ T * (  1.24406669e-05 + T * ( -9.48213152e-09
			+ T *  2.63545304e-12 ) ) ) ) -1.43089567e+04 );
		cp[sCH2O] =  2.76904683e+02 * (
			 4.79372315e+00 + T * ( -9.90833369e-03 
			+ T * (  3.73220008e-05 + T * ( -3.79285261e-08
			+ T *  1.31772652e-11 ) ) ) );
		h[sHCO] =  2.86523537e+02 * (
			T * (  4.22118584e+00 + T * ( -1.62196266e-03
			+ T * (  4.59331487e-06 + T * ( -3.32860233e-09
			+ T *  8.67537730e-13 ) ) ) ) +  3.83956496e+03 );
		cp[sHCO] =  2.86523537e+02 * (
			 4.22118584e+00 + T * ( -3.24392532e-03 
			+ T * (  1.37799446e-05 + T * ( -1.33144093e-08
			+ T *  4.33768865e-12 ) ) ) );
		h[sCH3] =  5.53035786e+02 * (
			T * (  3.65717970e+00 + T * (  1.06329895e-03
			+ T * (  1.81946277e-06 + T * ( -1.65452507e-09
			+ T *  4.93141480e-13 ) ) ) ) +  1.64227160e+04 );
		cp[sCH3] =  5.53035786e+02 * (
			 3.65717970e+00 + T * (  2.12659790e-03 
			+ T * (  5.45838830e-06 + T * ( -6.61810030e-09
			+ T *  2.46570740e-12 ) ) ) );
		h[sCO2] =  1.88919337e+02 * (
			T * (  2.27572500e+00 + T * (  4.96103600e-03
			+ T * ( -3.46970333e-06 + T * (  1.71667175e-09
			+ T * -4.23456000e-13 ) ) ) ) -4.83731400e+04 );
		cp[sCO2] =  1.88919337e+02 * (
			 2.27572500e+00 + T * (  9.92207200e-03 
			+ T * ( -1.04091100e-05 + T * (  6.86668700e-09
			+ T * -2.11728000e-12 ) ) ) );
		h[sCO] =  2.96834702e+02 * (
			T * (  3.26245200e+00 + T * (  7.55970500e-04
			+ T * ( -1.29391833e-06 + T * (  1.39548600e-09
			+ T * -4.94990200e-13 ) ) ) ) -1.43105400e+04 );
		cp[sCO] =  2.96834702e+02 * (
			 3.26245200e+00 + T * (  1.51194100e-03 
			+ T * ( -3.88175500e-06 + T * (  5.58194400e-09
			+ T * -2.47495100e-12 ) ) ) );
		h[sSXCH2] =  5.92780550e+02 * (
			T * (  4.19860411e+00 + T * ( -1.18330710e-03
			+ T * (  2.74432073e-06 + T * ( -1.67203995e-09
			+ T *  3.88629474e-13 ) ) ) ) +  5.04968163e+04 );
		cp[sSXCH2] =  5.92780550e+02 * (
			 4.19860411e+00 + T * ( -2.36661419e-03 
			+ T * (  8.23296220e-06 + T * ( -6.68815981e-09
			+ T *  1.94314737e-12 ) ) ) );
		h[sC2H4] =  2.96390275e+02 * (
			T * (  3.95920148e+00 + T * ( -3.78526124e-03
			+ T * (  1.90330097e-05 + T * ( -1.72897188e-08
			+ T *  5.39768746e-12 ) ) ) ) +  5.08977593e+03 );
		cp[sC2H4] =  2.96390275e+02 * (
			 3.95920148e+00 + T * ( -7.57052247e-03 
			+ T * (  5.70990292e-05 + T * ( -6.91588753e-08
			+ T *  2.69884373e-11 ) ) ) );
		h[sCH3O] =  2.67910679e+02 * (
			T * (  2.10620400e+00 + T * (  3.60829750e-03
			+ T * (  1.77949067e-06 + T * ( -1.84440900e-09
			+ T *  4.15122000e-13 ) ) ) ) +  9.78601100e+02 );
		cp[sCH3O] =  2.67910679e+02 * (
			 2.10620400e+00 + T * (  7.21659500e-03 
			+ T * (  5.33847200e-06 + T * ( -7.37763600e-09
			+ T *  2.07561000e-12 ) ) ) );
		h[sCH3O2] =  1.76772973e+02 * (
			T * (  4.76597792e+00 + T * ( -1.75538574e-03
			+ T * (  1.51464717e-05 + T * ( -1.41690932e-08
			+ T *  4.43182964e-12 ) ) ) ) -4.82401289e+02 );
		cp[sCH3O2] =  1.76772973e+02 * (
			 4.76597792e+00 + T * ( -3.51077148e-03 
			+ T * (  4.54394152e-05 + T * ( -5.66763729e-08
			+ T *  2.21591482e-11 ) ) ) );
		h[sC2H5] =  2.86109429e+02 * (
			T * (  4.30646568e+00 + T * ( -2.09329446e-03
			+ T * (  1.65714269e-05 + T * ( -1.49781651e-08
			+ T *  4.61018008e-12 ) ) ) ) +  1.28416265e+04 );
		cp[sC2H5] =  2.86109429e+02 * (
			 4.30646568e+00 + T * ( -4.18658892e-03 
			+ T * (  4.97142807e-05 + T * ( -5.99126606e-08
			+ T *  2.30509004e-11 ) ) ) );
		h[sCH2CO] =  1.97790941e+02 * (
			T * (  2.13583630e+00 + T * (  9.05943605e-03
			+ T * ( -5.79824913e-06 + T * (  2.33599392e-09
			+ T * -4.02915230e-13 ) ) ) ) -7.04291804e+03 );
		cp[sCH2CO] =  1.97790941e+02 * (
			 2.13583630e+00 + T * (  1.81188721e-02 
			+ T * ( -1.73947474e-05 + T * (  9.34397568e-09
			+ T * -2.01457615e-12 ) ) ) );
		h[sO2CHO] =  1.36260448e+02 * (
			T * (  3.96059309e+00 + T * (  5.30011395e-03
			+ T * ( -1.75237784e-06 + T * (  2.54291815e-10
			+ T * -5.74975204e-15 ) ) ) ) -1.73599383e+04 );
		cp[sO2CHO] =  1.36260448e+02 * (
			 3.96059309e+00 + T * (  1.06002279e-02 
			+ T * ( -5.25713351e-06 + T * (  1.01716726e-09
			+ T * -2.87487602e-14 ) ) ) );
		h[sHOCHO] =  1.80644418e+02 * (
			T * (  1.28069021e+00 + T * (  7.64438790e-03
			+ T * ( -1.88050159e-06 + T * ( -3.07421998e-10
			+ T *  1.62854647e-13 ) ) ) ) -4.64347524e+04 );
		cp[sHOCHO] =  1.80644418e+02 * (
			 1.28069021e+00 + T * (  1.52887758e-02 
			+ T * ( -5.64150476e-06 + T * ( -1.22968799e-09
			+ T *  8.14273233e-13 ) ) ) );
		h[sOCHO] =  1.84689235e+02 * (
			T * (  1.42991854e+00 + T * (  6.14016420e-03
			+ T * ( -1.76153699e-06 + T * ( -9.32566377e-11
			+ T *  1.01590714e-13 ) ) ) ) -1.62042814e+04 );
		cp[sOCHO] =  1.84689235e+02 * (
			 1.42991854e+00 + T * (  1.22803284e-02 
			+ T * ( -5.28461096e-06 + T * ( -3.73026551e-10
			+ T *  5.07953571e-13 ) ) ) );
		h[sC2H2] =  3.19340144e+02 * (
			T * (  8.08681094e-01 + T * (  1.16807815e-02
			+ T * ( -1.18390605e-05 + T * (  7.00381092e-09
			+ T * -1.70014595e-12 ) ) ) ) +  2.64289807e+04 );
		cp[sC2H2] =  3.19340144e+02 * (
			 8.08681094e-01 + T * (  2.33615629e-02 
			+ T * ( -3.55171815e-05 + T * (  2.80152437e-08
			+ T * -8.50072974e-12 ) ) ) );
		h[sC2H3] =  3.07437509e+02 * (
			T * (  3.21246645e+00 + T * (  7.57395810e-04
			+ T * (  8.64031373e-06 + T * ( -8.94144617e-09
			+ T *  2.94301746e-12 ) ) ) ) +  3.48598468e+04 );
		cp[sC2H3] =  3.07437509e+02 * (
			 3.21246645e+00 + T * (  1.51479162e-03 
			+ T * (  2.59209412e-05 + T * ( -3.57657847e-08
			+ T *  1.47150873e-11 ) ) ) );
		h[sHCCO] =  2.02650385e+02 * (
			T * (  2.25172140e+00 + T * (  8.82751050e-03
			+ T * ( -7.90970033e-06 + T * (  4.31893975e-09
			+ T * -1.01329622e-12 ) ) ) ) +  2.00594490e+04 );
		cp[sHCCO] =  2.02650385e+02 * (
			 2.25172140e+00 + T * (  1.76550210e-02 
			+ T * ( -2.37291010e-05 + T * (  1.72757590e-08
			+ T * -5.06648110e-12 ) ) ) );
		h[sCH2CHO] =  1.93159093e+02 * (
			T * (  1.09685733e+00 + T * (  1.10114398e-02
			+ T * ( -4.81944813e-06 + T * (  7.51948945e-10
			+ T *  1.21798575e-13 ) ) ) ) +  1.06943322e+03 );
		cp[sCH2CHO] =  1.93159093e+02 * (
			 1.09685733e+00 + T * (  2.20228796e-02 
			+ T * ( -1.44583444e-05 + T * (  3.00779578e-09
			+ T *  6.08992877e-13 ) ) ) );
		h[sC2H5O] =  1.84517088e+02 * (
			T * (  4.94420708e-01 + T * (  1.35887217e-02
			+ T * ( -5.53030033e-06 + T * (  1.28801050e-09
			+ T * -1.29699383e-13 ) ) ) ) -3.35252925e+03 );
		cp[sC2H5O] =  1.84517088e+02 * (
			 4.94420708e-01 + T * (  2.71774434e-02 
			+ T * ( -1.65909010e-05 + T * (  5.15204200e-09
			+ T * -6.48496915e-13 ) ) ) );
		h[sAXC3H5] =  2.02443146e+02 * (
			T * ( -1.03516444e+00 + T * (  1.87521683e-02
			+ T * ( -1.08793747e-05 + T * (  3.69156533e-09
			+ T * -4.87482308e-13 ) ) ) ) +  1.88792254e+04 );
		cp[sAXC3H5] =  2.02443146e+02 * (
			-1.03516444e+00 + T * (  3.75043366e-02 
			+ T * ( -3.26381242e-05 + T * (  1.47662613e-08
			+ T * -2.43741154e-12 ) ) ) );
		h[sC3H6] =  1.97593517e+02 * (
			T * ( -2.29261670e-03 + T * (  1.55130533e-02
			+ T * ( -5.57171827e-06 + T * (  4.73985425e-10
			+ T *  2.49915830e-13 ) ) ) ) +  1.13437406e+03 );
		cp[sC3H6] =  1.97593517e+02 * (
			-2.29261670e-03 + T * (  3.10261065e-02 
			+ T * ( -1.67151548e-05 + T * (  1.89594170e-09
			+ T *  1.24957915e-12 ) ) ) );
		h[sCH3COCH2] =  1.45686701e+02 * (
			T * (  1.22337251e+00 + T * (  1.62273371e-02
			+ T * ( -7.11808393e-06 + T * (  1.74194434e-09
			+ T * -1.79832060e-13 ) ) ) ) -6.59419324e+03 );
		cp[sCH3COCH2] =  1.45686701e+02 * (
			 1.22337251e+00 + T * (  3.24546742e-02 
			+ T * ( -2.13542518e-05 + T * (  6.96777735e-09
			+ T * -8.99160299e-13 ) ) ) );
		h[sNXC3H7] =  1.92970803e+02 * (
			T * (  1.04754730e+00 + T * (  1.30038970e-02
			+ T * (  7.85408400e-07 + T * ( -4.89807925e-09
			+ T *  1.87360232e-12 ) ) ) ) +  1.06326370e+04 );
		cp[sNXC3H7] =  1.92970803e+02 * (
			 1.04754730e+00 + T * (  2.60077940e-02 
			+ T * (  2.35622520e-06 + T * ( -1.95923170e-08
			+ T *  9.36801160e-12 ) ) ) );
		h[sPXC4H8] =  1.48195138e+02 * (
			T * ( -8.31372089e-01 + T * (  2.26290489e-02
			+ T * ( -9.78861863e-06 + T * (  2.50551090e-09
			+ T * -2.86383360e-13 ) ) ) ) -1.57875035e+03 );
		cp[sPXC4H8] =  1.48195138e+02 * (
			-8.31372089e-01 + T * (  4.52580978e-02 
			+ T * ( -2.93658559e-05 + T * (  1.00220436e-08
			+ T * -1.43191680e-12 ) ) ) );
		h[sPXC4H9] =  1.45579563e+02 * (
			T * ( -4.37779725e-01 + T * (  2.39486182e-02
			+ T * ( -1.04674386e-05 + T * (  2.74466180e-09
			+ T * -3.24021328e-13 ) ) ) ) +  7.68945248e+03 );
		cp[sPXC4H9] =  1.45579563e+02 * (
			-4.37779725e-01 + T * (  4.78972364e-02 
			+ T * ( -3.14023159e-05 + T * (  1.09786472e-08
			+ T * -1.62010664e-12 ) ) ) );
		h[sC5H11] =  1.16876212e+02 * (
			T * ( -9.05255912e-01 + T * (  3.05316426e-02
			+ T * ( -1.36497275e-05 + T * (  3.65233675e-09
			+ T * -4.37719230e-13 ) ) ) ) +  4.83995303e+03 );
		cp[sC5H11] =  1.16876212e+02 * (
			-9.05255912e-01 + T * (  6.10632852e-02 
			+ T * ( -4.09491825e-05 + T * (  1.46093470e-08
			+ T * -2.18859615e-12 ) ) ) );
		h[sC6H12X1XC6H12] =  9.87967584e+01 * (
			T * ( -1.35275205e+00 + T * (  3.49327713e-02
			+ T * ( -1.53136007e-05 + T * (  3.92418358e-09
			+ T * -4.42592350e-13 ) ) ) ) -7.34368617e+03 );
		cp[sC6H12X1XC6H12] =  9.87967584e+01 * (
			-1.35275205e+00 + T * (  6.98655426e-02 
			+ T * ( -4.59408022e-05 + T * (  1.56967343e-08
			+ T * -2.21296175e-12 ) ) ) );
		h[sC6H13] =  9.76274012e+01 * (
			T * ( -2.04871465e-01 + T * (  3.41900636e-02
			+ T * ( -1.38149304e-05 + T * (  3.15389505e-09
			+ T * -3.06240116e-13 ) ) ) ) +  1.83280393e+03 );
		cp[sC6H13] =  9.76274012e+01 * (
			-2.04871465e-01 + T * (  6.83801272e-02 
			+ T * ( -4.14447912e-05 + T * (  1.26155802e-08
			+ T * -1.53120058e-12 ) ) ) );
		h[sC7H14] =  8.46829358e+01 * (
			T * ( -2.03026994e+00 + T * (  4.13162189e-02
			+ T * ( -1.81838157e-05 + T * (  4.69264555e-09
			+ T * -5.35142440e-13 ) ) ) ) -1.15141029e+04 );
		cp[sC7H14] =  8.46829358e+01 * (
			-2.03026994e+00 + T * (  8.26324377e-02 
			+ T * ( -5.45514471e-05 + T * (  1.87705822e-08
			+ T * -2.67571220e-12 ) ) ) );
		h[sC7H15] =  8.38223611e+01 * (
			T * ( -3.79155767e-02 + T * (  3.78363285e-02
			+ T * ( -1.35824545e-05 + T * (  2.33169736e-09
			+ T * -9.84721490e-14 ) ) ) ) -2.35605303e+03 );
		cp[sC7H15] =  8.38223611e+01 * (
			-3.79155767e-02 + T * (  7.56726570e-02 
			+ T * ( -4.07473634e-05 + T * (  9.32678943e-09
			+ T * -4.92360745e-13 ) ) ) );
		h[sC8H17] =  7.34378533e+01 * (
			T * ( -7.72759438e-01 + T * (  4.66274853e-02
			+ T * ( -1.94815748e-05 + T * (  4.63925535e-09
			+ T * -4.74254966e-13 ) ) ) ) -3.92689511e+03 );
		cp[sC8H17] =  7.34378533e+01 * (
			-7.72759438e-01 + T * (  9.32549705e-02 
			+ T * ( -5.84447245e-05 + T * (  1.85570214e-08
			+ T * -2.37127483e-12 ) ) ) );
		h[sC9H19] =  6.53427327e+01 * (
			T * ( -1.04387292e+00 + T * (  5.28086415e-02
			+ T * ( -2.22733324e-05 + T * (  5.36215415e-09
			+ T * -5.54808550e-13 ) ) ) ) -6.80818512e+03 );
		cp[sC9H19] =  6.53427327e+01 * (
			-1.04387292e+00 + T * (  1.05617283e-01 
			+ T * ( -6.68199971e-05 + T * (  2.14486166e-08
			+ T * -2.77404275e-12 ) ) ) );
		h[sNXC12H26] =  4.88137006e+01 * (
			T * ( -2.62181594e+00 + T * (  7.36188555e-02
			+ T * ( -3.14656757e-05 + T * (  7.68603170e-09
			+ T * -8.07204460e-13 ) ) ) ) -4.00654253e+04 );
		cp[sNXC12H26] =  4.88137006e+01 * (
			-2.62181594e+00 + T * (  1.47237711e-01 
			+ T * ( -9.43970271e-05 + T * (  3.07441268e-08
			+ T * -4.03602230e-12 ) ) ) );
		h[sC12H25] =  4.91042996e+01 * (
			T * ( -1.85028741e+00 + T * (  7.13353540e-02
			+ T * ( -3.06305518e-05 + T * (  7.52208480e-09
			+ T * -7.94908600e-13 ) ) ) ) -1.54530435e+04 );
		cp[sC12H25] =  4.91042996e+01 * (
			-1.85028741e+00 + T * (  1.42670708e-01 
			+ T * ( -9.18916555e-05 + T * (  3.00883392e-08
			+ T * -3.97454300e-12 ) ) ) );
		h[sC12H24] =  4.93983792e+01 * (
			T * ( -2.43874931e+00 + T * (  7.02564230e-02
			+ T * ( -3.03344326e-05 + T * (  7.50771467e-09
			+ T * -8.01580282e-13 ) ) ) ) -2.61556835e+04 );
		cp[sC12H24] =  4.93983792e+01 * (
			-2.43874931e+00 + T * (  1.40512846e-01 
			+ T * ( -9.10032977e-05 + T * (  3.00308587e-08
			+ T * -4.00790141e-12 ) ) ) );
		h[sC12H25O2] =  4.12991258e+01 * (
			T * (  1.58304260e+00 + T * (  7.13840480e-02
			+ T * ( -2.99213936e-05 + T * (  7.25500520e-09
			+ T * -7.75971112e-13 ) ) ) ) -3.50130986e+04 );
		cp[sC12H25O2] =  4.12991258e+01 * (
			 1.58304260e+00 + T * (  1.42768096e-01 
			+ T * ( -8.97641809e-05 + T * (  2.90200208e-08
			+ T * -3.87985556e-12 ) ) ) );
		h[sC12OOHXT6XC12H25O2] =  4.12991258e+01 * (
			T * (  1.42078475e+00 + T * (  7.20197365e-02
			+ T * ( -2.98554576e-05 + T * (  6.83165328e-09
			+ T * -6.42451046e-13 ) ) ) ) -2.85025839e+04 );
		cp[sC12OOHXT6XC12H25O2] =  4.12991258e+01 * (
			 1.42078475e+00 + T * (  1.44039473e-01 
			+ T * ( -8.95663728e-05 + T * (  2.73266131e-08
			+ T * -3.21225523e-12 ) ) ) );
		h[sC12OOHO2XT6XC12H25O4] =  3.56349220e+01 * (
			T * (  3.14841764e+00 + T * (  7.66022995e-02
			+ T * ( -3.21410091e-05 + T * (  7.36427198e-09
			+ T * -6.90304458e-13 ) ) ) ) -4.82915577e+04 );
		cp[sC12OOHO2XT6XC12H25O4] =  3.56349220e+01 * (
			 3.14841764e+00 + T * (  1.53204599e-01 
			+ T * ( -9.64230273e-05 + T * (  2.94570879e-08
			+ T * -3.45152229e-12 ) ) ) );
		h[sC12KETXT6XC12H24O3] =  3.84367950e+01 * (
			T * (  2.09319905e+00 + T * (  7.53591275e-02
			+ T * ( -3.44404850e-05 + T * (  9.25207215e-09
			+ T * -1.10102720e-12 ) ) ) ) -6.36120922e+04 );
		cp[sC12KETXT6XC12H24O3] =  3.84367950e+01 * (
			 2.09319905e+00 + T * (  1.50718255e-01 
			+ T * ( -1.03321455e-04 + T * (  3.70082886e-08
			+ T * -5.50513600e-12 ) ) ) );
		h[sNXC8H17CHO] =  5.84553623e+01 * (
			T * (  7.56488136e-01 + T * (  5.14860075e-02
			+ T * ( -2.16636960e-05 + T * (  5.27230180e-09
			+ T * -5.64192190e-13 ) ) ) ) -4.15508169e+04 );
		cp[sNXC8H17CHO] =  5.84553623e+01 * (
			 7.56488136e-01 + T * (  1.02972015e-01 
			+ T * ( -6.49910881e-05 + T * (  2.10892072e-08
			+ T * -2.82096095e-12 ) ) ) );
		h[sC7H15CO] =  6.53643082e+01 * (
			T * (  1.73697916e+00 + T * (  4.30131311e-02
			+ T * ( -1.79164136e-05 + T * (  4.27960817e-09
			+ T * -4.45543698e-13 ) ) ) ) -2.01298819e+04 );
		cp[sC7H15CO] =  6.53643082e+01 * (
			 1.73697916e+00 + T * (  8.60262622e-02 
			+ T * ( -5.37492409e-05 + T * (  1.71184327e-08
			+ T * -2.22771849e-12 ) ) ) );
		h[sC3H7CHO] =  1.15310385e+02 * (
			T * (  1.87415959e+00 + T * (  2.09620158e-02
			+ T * ( -7.83829263e-06 + T * (  1.56728418e-09
			+ T * -1.21888782e-13 ) ) ) ) -2.71032194e+04 );
		cp[sC3H7CHO] =  1.15310385e+02 * (
			 1.87415959e+00 + T * (  4.19240315e-02 
			+ T * ( -2.35148779e-05 + T * (  6.26913673e-09
			+ T * -6.09443908e-13 ) ) ) );
		h[sNXC7H15CHO] =  6.48503994e+01 * (
			T * (  9.66850312e-01 + T * (  4.54031898e-02
			+ T * ( -1.89236823e-05 + T * (  4.54126703e-09
			+ T * -4.76993758e-13 ) ) ) ) -3.86587355e+04 );
		cp[sNXC7H15CHO] =  6.48503994e+01 * (
			 9.66850312e-01 + T * (  9.08063795e-02 
			+ T * ( -5.67710469e-05 + T * (  1.81650681e-08
			+ T * -2.38496879e-12 ) ) ) );
		h[sC7H15O2X1XC7H15O2] =  6.33763244e+01 * (
			T * (  1.20888552e+00 + T * (  4.42030853e-02
			+ T * ( -1.93287087e-05 + T * (  5.08936723e-09
			+ T * -6.12920672e-13 ) ) ) ) -1.92918765e+04 );
		cp[sC7H15O2X1XC7H15O2] =  6.33763244e+01 * (
			 1.20888552e+00 + T * (  8.84061706e-02 
			+ T * ( -5.79861260e-05 + T * (  2.03574689e-08
			+ T * -3.06460336e-12 ) ) ) );
		h[sC7H14OOH1X3XC7H15O2] =  6.33763244e+01 * (
			T * (  1.11146206e+00 + T * (  4.45797696e-02
			+ T * ( -1.89117461e-05 + T * (  4.47225235e-09
			+ T * -4.44520412e-13 ) ) ) ) -1.29868314e+04 );
		cp[sC7H14OOH1X3XC7H15O2] =  6.33763244e+01 * (
			 1.11146206e+00 + T * (  8.91595392e-02 
			+ T * ( -5.67352382e-05 + T * (  1.78890094e-08
			+ T * -2.22260206e-12 ) ) ) );
		h[sC5H4CH2] =  1.06446715e+02 * (
			T * ( -5.34007612e+00 + T * (  3.58641913e-02
			+ T * ( -2.15274819e-05 + T * (  6.96727892e-09
			+ T * -7.90002910e-13 ) ) ) ) +  2.58936616e+04 );
		cp[sC5H4CH2] =  1.06446715e+02 * (
			-5.34007612e+00 + T * (  7.17283827e-02 
			+ T * ( -6.45824457e-05 + T * (  2.78691157e-08
			+ T * -3.95001455e-12 ) ) ) );
		h[sA1XC6H6] =  1.06446715e+02 * (
			T * ( -5.51558393e+00 + T * (  3.22726613e-02
			+ T * ( -1.47134309e-05 + T * (  1.86928040e-09
			+ T *  6.20564508e-13 ) ) ) ) +  9.11031457e+03 );
		cp[sA1XC6H6] =  1.06446715e+02 * (
			-5.51558393e+00 + T * (  6.45453225e-02 
			+ T * ( -4.41402928e-05 + T * (  7.47712161e-09
			+ T *  3.10282254e-12 ) ) ) );
		h[sNXC4H5] =  1.56637905e+02 * (
			T * ( -1.16849950e+00 + T * (  2.39503037e-02
			+ T * ( -1.70792334e-05 + T * (  7.65610660e-09
			+ T * -1.51981393e-12 ) ) ) ) +  4.22787216e+04 );
		cp[sNXC4H5] =  1.56637905e+02 * (
			-1.16849950e+00 + T * (  4.79006074e-02 
			+ T * ( -5.12377002e-05 + T * (  3.06244264e-08
			+ T * -7.59906965e-12 ) ) ) );
		h[sIXC4H5] =  1.56637905e+02 * (
			T * ( -3.31905498e-01 + T * (  2.20081938e-02
			+ T * ( -1.42563415e-05 + T * (  5.78210790e-09
			+ T * -1.03434304e-12 ) ) ) ) +  3.67510686e+04 );
		cp[sIXC4H5] =  1.56637905e+02 * (
			-3.31905498e-01 + T * (  4.40163876e-02 
			+ T * ( -4.27690246e-05 + T * (  2.31284316e-08
			+ T * -5.17171519e-12 ) ) ) );
		h[sA1XXC6H5] =  1.07838392e+02 * (
			T * ( -4.87654845e+00 + T * (  3.13402891e-02
			+ T * ( -1.62467429e-05 + T * (  3.52805717e-09
			+ T *  1.03703662e-13 ) ) ) ) +  3.99269438e+04 );
		cp[sA1XXC6H5] =  1.07838392e+02 * (
			-4.87654845e+00 + T * (  6.26805782e-02 
			+ T * ( -4.87402286e-05 + T * (  1.41122287e-08
			+ T *  5.18518312e-13 ) ) ) );
		h[sNXC4H3] =  1.62821949e+02 * (
			T * ( -3.55175031e-02 + T * (  2.15254251e-02
			+ T * ( -1.91909716e-05 + T * (  1.03970785e-08
			+ T * -2.41501714e-12 ) ) ) ) +  6.43506593e+04 );
		cp[sNXC4H3] =  1.62821949e+02 * (
			-3.55175031e-02 + T * (  4.30508503e-02 
			+ T * ( -5.75729147e-05 + T * (  4.15883142e-08
			+ T * -1.20750857e-11 ) ) ) );
		h[sC3H3] =  2.12893430e+02 * (
			T * (  1.40299238e+00 + T * (  1.50886664e-02
			+ T * ( -1.32816458e-05 + T * (  7.33836572e-09
			+ T * -1.74110916e-12 ) ) ) ) +  3.93108220e+04 );
		cp[sC3H3] =  2.12893430e+02 * (
			 1.40299238e+00 + T * (  3.01773327e-02 
			+ T * ( -3.98449373e-05 + T * (  2.93534629e-08
			+ T * -8.70554579e-12 ) ) ) );
		h[sA1C2H2XC8H7] =  8.06153041e+01 * (
			T * ( -6.30997035e+00 + T * (  4.75453915e-02
			+ T * ( -3.18566445e-05 + T * (  1.24202002e-08
			+ T * -2.03584362e-12 ) ) ) ) +  4.57329298e+04 );
		cp[sA1C2H2XC8H7] =  8.06153041e+01 * (
			-6.30997035e+00 + T * (  9.50907829e-02 
			+ T * ( -9.55699336e-05 + T * (  4.96808010e-08
			+ T * -1.01792181e-11 ) ) ) );
		h[sA1C2H3XC8H8] =  7.98350361e+01 * (
			T * ( -5.38499941e+00 + T * (  4.10182578e-02
			+ T * ( -1.78153959e-05 + T * (  1.39773752e-09
			+ T *  1.12227810e-12 ) ) ) ) +  1.60857559e+04 );
		cp[sA1C2H3XC8H8] =  7.98350361e+01 * (
			-5.38499941e+00 + T * (  8.20365155e-02 
			+ T * ( -5.34461878e-05 + T * (  5.59095007e-09
			+ T *  5.61139050e-12 ) ) ) );
		h[sA1C2HXC8H6] =  8.14109745e+01 * (
			T * ( -5.21036925e+00 + T * (  4.32775972e-02
			+ T * ( -2.81669161e-05 + T * (  1.05480176e-08
			+ T * -1.63353233e-12 ) ) ) ) +  3.52488620e+04 );
		cp[sA1C2HXC8H6] =  8.14109745e+01 * (
			-5.21036925e+00 + T * (  8.65551944e-02 
			+ T * ( -8.45007483e-05 + T * (  4.21920706e-08
			+ T * -8.16766167e-12 ) ) ) );
		h[sA1C2HYXC8H5] =  8.22225079e+01 * (
			T * ( -4.42757639e+00 + T * (  4.18334322e-02
			+ T * ( -2.90035454e-05 + T * (  1.17571415e-08
			+ T * -2.03633970e-12 ) ) ) ) +  6.73302359e+04 );
		cp[sA1C2HYXC8H5] =  8.22225079e+01 * (
			-4.42757639e+00 + T * (  8.36668645e-02 
			+ T * ( -8.70106362e-05 + T * (  4.70285661e-08
			+ T * -1.01816985e-11 ) ) ) );
		h[sA1C2H3YXC8H7] =  8.06153041e+01 * (
			T * ( -5.36214520e+00 + T * (  4.33516648e-02
			+ T * ( -2.51432653e-05 + T * (  7.52849635e-09
			+ T * -6.81362836e-13 ) ) ) ) +  4.77818209e+04 );
		cp[sA1C2H3YXC8H7] =  8.06153041e+01 * (
			-5.36214520e+00 + T * (  8.67033297e-02 
			+ T * ( -7.54297960e-05 + T * (  3.01139854e-08
			+ T * -3.40681418e-12 ) ) ) );
		h[sA2XXC10H7] =  6.53869263e+01 * (
			T * ( -8.02718034e+00 + T * (  5.14622590e-02
			+ T * ( -2.78090670e-05 + T * (  6.80338458e-09
			+ T * -1.44911911e-13 ) ) ) ) +  5.01363344e+04 );
		cp[sA2XXC10H7] =  6.53869263e+01 * (
			-8.02718034e+00 + T * (  1.02924518e-01 
			+ T * ( -8.34272010e-05 + T * (  2.72135383e-08
			+ T * -7.24559554e-13 ) ) ) );
		h[sA2XC10H8] =  6.48726632e+01 * (
			T * ( -8.72434585e+00 + T * (  5.26880040e-02
			+ T * ( -2.67236897e-05 + T * (  5.46364935e-09
			+ T *  2.84133212e-13 ) ) ) ) +  1.66588912e+04 );
		cp[sA2XC10H8] =  6.48726632e+01 * (
			-8.72434585e+00 + T * (  1.05376008e-01 
			+ T * ( -8.01710690e-05 + T * (  2.18545974e-08
			+ T *  1.42066606e-12 ) ) ) );
		h[sC4H4] =  1.59670072e+02 * (
			T * ( -2.31343354e-01 + T * (  2.05907248e-02
			+ T * ( -1.49208019e-05 + T * (  6.88585392e-09
			+ T * -1.41275363e-12 ) ) ) ) +  3.40632704e+04 );
		cp[sC4H4] =  1.59670072e+02 * (
			-2.31343354e-01 + T * (  4.11814497e-02 
			+ T * ( -4.47624056e-05 + T * (  2.75434157e-08
			+ T * -7.06376813e-12 ) ) ) );
		h[sA2YXC10H7] =  6.53869263e+01 * (
			T * ( -8.00768796e+00 + T * (  5.15206445e-02
			+ T * ( -2.79396999e-05 + T * (  6.91229315e-09
			+ T * -1.77768442e-13 ) ) ) ) +  4.99740633e+04 );
		cp[sA2YXC10H7] =  6.53869263e+01 * (
			-8.00768796e+00 + T * (  1.03041289e-01 
			+ T * ( -8.38190998e-05 + T * (  2.76491726e-08
			+ T * -8.88842208e-13 ) ) ) );
		h[sA2C2H2AXC12H9] =  5.42739830e+01 * (
			T * ( -9.26784872e+00 + T * (  6.75215650e-02
			+ T * ( -4.32639473e-05 + T * (  1.56805083e-08
			+ T * -2.32697154e-12 ) ) ) ) +  5.64832554e+04 );
		cp[sA2C2H2AXC12H9] =  5.42739830e+01 * (
			-9.26784872e+00 + T * (  1.35043130e-01 
			+ T * ( -1.29791842e-04 + T * (  6.27220331e-08
			+ T * -1.16348577e-11 ) ) ) );
		h[sA2C2H2BXC12H9] =  5.42739830e+01 * (
			T * ( -9.38518818e+00 + T * (  6.74128980e-02
			+ T * ( -4.29874243e-05 + T * (  1.53620107e-08
			+ T * -2.18934312e-12 ) ) ) ) +  5.63724233e+04 );
		cp[sA2C2H2BXC12H9] =  5.42739830e+01 * (
			-9.38518818e+00 + T * (  1.34825796e-01 
			+ T * ( -1.28962273e-04 + T * (  6.14480428e-08
			+ T * -1.09467156e-11 ) ) ) );
		h[sA2C2HAXC12H8] =  5.46334700e+01 * (
			T * ( -8.23047877e+00 + T * (  6.30265880e-02
			+ T * ( -3.91663747e-05 + T * (  1.34197745e-08
			+ T * -1.77092462e-12 ) ) ) ) +  4.23747629e+04 );
		cp[sA2C2HAXC12H8] =  5.46334700e+01 * (
			-8.23047877e+00 + T * (  1.26053176e-01 
			+ T * ( -1.17499124e-04 + T * (  5.36790980e-08
			+ T * -8.85462308e-12 ) ) ) );
		h[sA2C2HBXC12H8] =  5.46334700e+01 * (
			T * ( -8.22579974e+00 + T * (  6.31237755e-02
			+ T * ( -3.93802473e-05 + T * (  1.35996841e-08
			+ T * -1.82516934e-12 ) ) ) ) +  4.25495114e+04 );
		cp[sA2C2HBXC12H8] =  5.46334700e+01 * (
			-8.22579974e+00 + T * (  1.26247551e-01 
			+ T * ( -1.18140742e-04 + T * (  5.43987363e-08
			+ T * -9.12584671e-12 ) ) ) );
		h[sA2C2HAYXC12H7] =  5.49977510e+01 * (
			T * ( -7.36585075e+00 + T * (  6.12789955e-02
			+ T * ( -3.95790920e-05 + T * (  1.43485070e-08
			+ T * -2.10117262e-12 ) ) ) ) +  7.69836555e+04 );
		cp[sA2C2HAYXC12H7] =  5.49977510e+01 * (
			-7.36585075e+00 + T * (  1.22557991e-01 
			+ T * ( -1.18737276e-04 + T * (  5.73940282e-08
			+ T * -1.05058631e-11 ) ) ) );
		h[sA2C2HBYXC12H7] =  5.49977510e+01 * (
			T * ( -7.35510706e+00 + T * (  6.12626680e-02
			+ T * ( -3.95804270e-05 + T * (  1.43577795e-08
			+ T * -2.10465530e-12 ) ) ) ) +  7.73354743e+04 );
		cp[sA2C2HBYXC12H7] =  5.49977510e+01 * (
			-7.35510706e+00 + T * (  1.22525336e-01 
			+ T * ( -1.18741281e-04 + T * (  5.74311178e-08
			+ T * -1.05232765e-11 ) ) ) );
		h[sA2R5XC12H8] =  5.46334700e+01 * (
			T * ( -1.05497902e+01 + T * (  6.27683950e-02
			+ T * ( -3.45486817e-05 + T * (  8.82472825e-09
			+ T * -3.29016768e-13 ) ) ) ) +  2.94426605e+04 );
		cp[sA2R5XC12H8] =  5.46334700e+01 * (
			-1.05497902e+01 + T * (  1.25536790e-01 
			+ T * ( -1.03646045e-04 + T * (  3.52989130e-08
			+ T * -1.64508384e-12 ) ) ) );
		h[sA2R5XXC12H7] =  5.49977510e+01 * (
			T * ( -9.79699017e+00 + T * (  6.11386065e-02
			+ T * ( -3.49775030e-05 + T * (  9.69866022e-09
			+ T * -6.32579322e-13 ) ) ) ) +  6.24840181e+04 );
		cp[sA2R5XXC12H7] =  5.49977510e+01 * (
			-9.79699017e+00 + T * (  1.22277213e-01 
			+ T * ( -1.04932509e-04 + T * (  3.87946409e-08
			+ T * -3.16289661e-12 ) ) ) );
		h[sA2R5C2H2XC14H9] =  4.69174774e+01 * (
			T * ( -9.79888742e+00 + T * (  7.27002575e-02
			+ T * ( -4.27866497e-05 + T * (  1.28967134e-08
			+ T * -1.20348472e-12 ) ) ) ) +  6.87094743e+04 );
		cp[sA2R5C2H2XC14H9] =  4.69174774e+01 * (
			-9.79888742e+00 + T * (  1.45400515e-01 
			+ T * ( -1.28359949e-04 + T * (  5.15868534e-08
			+ T * -6.01742362e-12 ) ) ) );
		h[sA2R5C2HXC14H8] =  4.71858755e+01 * (
			T * ( -9.95199604e+00 + T * (  7.27067810e-02
			+ T * ( -4.64002687e-05 + T * (  1.63758248e-08
			+ T * -2.27673192e-12 ) ) ) ) +  5.53993662e+04 );
		cp[sA2R5C2HXC14H8] =  4.71858755e+01 * (
			-9.95199604e+00 + T * (  1.45413562e-01 
			+ T * ( -1.39200806e-04 + T * (  6.55032994e-08
			+ T * -1.13836596e-11 ) ) ) );
		h[sA2R5C2HYXC14H7] =  4.74573620e+01 * (
			T * ( -9.09090029e+00 + T * (  7.09590820e-02
			+ T * ( -4.67904557e-05 + T * (  1.72745595e-08
			+ T * -2.59576030e-12 ) ) ) ) +  8.92211518e+04 );
		cp[sA2R5C2HYXC14H7] =  4.74573620e+01 * (
			-9.09090029e+00 + T * (  1.41918164e-01 
			+ T * ( -1.40371367e-04 + T * (  6.90982381e-08
			+ T * -1.29788015e-11 ) ) ) );
		h[sP2XC12H10] =  5.39191958e+01 * (
			T * ( -1.19438051e+01 + T * (  7.10815795e-02
			+ T * ( -4.44991497e-05 + T * (  1.55126429e-08
			+ T * -2.11533328e-12 ) ) ) ) +  2.01936932e+04 );
		cp[sP2XC12H10] =  5.39191958e+01 * (
			-1.19438051e+01 + T * (  1.42163159e-01 
			+ T * ( -1.33497449e-04 + T * (  6.20505718e-08
			+ T * -1.05766664e-11 ) ) ) );
		h[sP2XXC12H9] =  5.42739830e+01 * (
			T * ( -9.50091731e+00 + T * (  6.26052515e-02
			+ T * ( -3.27906074e-05 + T * (  6.89835438e-09
			+ T *  3.35908690e-13 ) ) ) ) +  5.29903482e+04 );
		cp[sP2XXC12H9] =  5.42739830e+01 * (
			-9.50091731e+00 + T * (  1.25210503e-01 
			+ T * ( -9.83718223e-05 + T * (  2.75934175e-08
			+ T *  1.67954345e-12 ) ) ) );
		h[sA3XXC14H9] =  4.69174774e+01 * (
			T * ( -1.08881743e+01 + T * (  7.05935385e-02
			+ T * ( -3.78438093e-05 + T * (  8.97938392e-09
			+ T * -9.00424064e-14 ) ) ) ) +  5.70158539e+04 );
		cp[sA3XXC14H9] =  4.69174774e+01 * (
			-1.08881743e+01 + T * (  1.41187077e-01 
			+ T * ( -1.13531428e-04 + T * (  3.59175357e-08
			+ T * -4.50212032e-13 ) ) ) );
		h[sA3XC14H10] =  4.66521154e+01 * (
			T * ( -1.15461369e+01 + T * (  7.18790815e-02
			+ T * ( -3.69563747e-05 + T * (  7.80450353e-09
			+ T *  2.91950464e-13 ) ) ) ) +  2.21687904e+04 );
		cp[sA3XC14H10] =  4.66521154e+01 * (
			-1.15461369e+01 + T * (  1.43758163e-01 
			+ T * ( -1.10869124e-04 + T * (  3.12180141e-08
			+ T *  1.45975232e-12 ) ) ) );
		h[sA3YXC14H9] =  4.69174774e+01 * (
			T * ( -1.08881743e+01 + T * (  7.05935385e-02
			+ T * ( -3.78438093e-05 + T * (  8.97938392e-09
			+ T * -9.00424064e-14 ) ) ) ) +  5.70158539e+04 );
		cp[sA3YXC14H9] =  4.69174774e+01 * (
			-1.08881743e+01 + T * (  1.41187077e-01 
			+ T * ( -1.13531428e-04 + T * (  3.59175357e-08
			+ T * -4.50212032e-13 ) ) ) );
		h[sA3R5XXC16H9] =  4.13171861e+01 * (
			T * ( -1.25419370e+01 + T * (  8.01855080e-02
			+ T * ( -4.50096677e-05 + T * (  1.19011023e-08
			+ T * -5.87324606e-13 ) ) ) ) +  6.87107686e+04 );
		cp[sA3R5XXC16H9] =  4.13171861e+01 * (
			-1.25419370e+01 + T * (  1.60371016e-01 
			+ T * ( -1.35029003e-04 + T * (  4.76044093e-08
			+ T * -2.93662303e-12 ) ) ) );
		h[sA3R5XC16H10] =  4.11112540e+01 * (
			T * ( -1.32241574e+01 + T * (  8.14313235e-02
			+ T * ( -4.39890187e-05 + T * (  1.06098466e-08
			+ T * -1.70634345e-13 ) ) ) ) +  3.44430693e+04 );
		cp[sA3R5XC16H10] =  4.11112540e+01 * (
			-1.32241574e+01 + T * (  1.62862647e-01 
			+ T * ( -1.31967056e-04 + T * (  4.24393865e-08
			+ T * -8.53171724e-13 ) ) ) );
		h[sA4XC16H10] =  4.11112540e+01 * (
			T * ( -1.31524443e+01 + T * (  8.04394215e-02
			+ T * ( -4.25732390e-05 + T * (  9.77297245e-09
			+ T *  1.48798225e-14 ) ) ) ) +  2.49673872e+04 );
		cp[sA4XC16H10] =  4.11112540e+01 * (
			-1.31524443e+01 + T * (  1.60878843e-01 
			+ T * ( -1.27719717e-04 + T * (  3.90918898e-08
			+ T *  7.43991125e-14 ) ) ) );
		h[sA4XXC16H9] =  4.13171861e+01 * (
			T * ( -1.23671835e+01 + T * (  7.88287575e-02
			+ T * ( -4.30812110e-05 + T * (  1.07157418e-08
			+ T * -3.08983758e-13 ) ) ) ) +  6.27797890e+04 );
		cp[sA4XXC16H9] =  4.13171861e+01 * (
			-1.23671835e+01 + T * (  1.57657515e-01 
			+ T * ( -1.29243633e-04 + T * (  4.28629673e-08
			+ T * -1.54491879e-12 ) ) ) );
		h[sA4R5XC18H10] =  3.67468399e+01 * (
			T * ( -1.47695663e+01 + T * (  8.98290110e-02
			+ T * ( -4.93966333e-05 + T * (  1.24430526e-08
			+ T * -4.12874670e-13 ) ) ) ) +  3.78467972e+04 );
		cp[sA4R5XC18H10] =  3.67468399e+01 * (
			-1.47695663e+01 + T * (  1.79658022e-01 
			+ T * ( -1.48189900e-04 + T * (  4.97722102e-08
			+ T * -2.06437335e-12 ) ) ) );
		h[sFLTNXC16H10] =  4.11112540e+01 * (
			T * ( -1.29396091e+01 + T * (  8.01598610e-02
			+ T * ( -4.21507513e-05 + T * (  9.40122770e-09
			+ T *  1.36967841e-13 ) ) ) ) +  2.91084888e+04 );
		cp[sFLTNXC16H10] =  4.11112540e+01 * (
			-1.29396091e+01 + T * (  1.60319722e-01 
			+ T * ( -1.26452254e-04 + T * (  3.76049108e-08
			+ T *  6.84839205e-13 ) ) ) );
		h[sC5H6] =  1.25788072e+02 * (
			T * ( -5.13691194e+00 + T * (  3.03476727e-02
			+ T * ( -1.53517612e-05 + T * (  3.21143002e-09
			+ T *  1.48242970e-13 ) ) ) ) +  1.53675713e+04 );
		cp[sC5H6] =  1.25788072e+02 * (
			-5.13691194e+00 + T * (  6.06953453e-02 
			+ T * ( -4.60552837e-05 + T * (  1.28457201e-08
			+ T *  7.41214852e-13 ) ) ) );
		h[sC5H5] =  1.27736058e+02 * (
			T * ( -7.37844042e+00 + T * (  4.86195909e-02
			+ T * ( -5.65263793e-05 + T * (  3.79546668e-08
			+ T * -1.02415096e-11 ) ) ) ) +  3.05514662e+04 );
		cp[sC5H5] =  1.27736058e+02 * (
			-7.37844042e+00 + T * (  9.72391818e-02 
			+ T * ( -1.69579138e-04 + T * (  1.51818667e-07
			+ T * -5.12075479e-11 ) ) ) );
		h[sSXC3H5] =  2.02443146e+02 * (
			T * (  3.13106581e-01 + T * (  1.59384831e-02
			+ T * ( -8.44733377e-06 + T * (  2.57497682e-09
			+ T * -2.70603708e-13 ) ) ) ) +  3.13767683e+04 );
		cp[sSXC3H5] =  2.02443146e+02 * (
			 3.13106581e-01 + T * (  3.18769663e-02 
			+ T * ( -2.53420013e-05 + T * (  1.02999073e-08
			+ T * -1.35301854e-12 ) ) ) );
		h[sCH4] =  5.18285750e+02 * (
			T * (  5.14911468e+00 + T * ( -6.83110045e-03
			+ T * (  1.63817974e-05 + T * ( -1.21061692e-08
			+ T *  3.33206882e-12 ) ) ) ) -1.02465983e+04 );
		cp[sCH4] =  5.18285750e+02 * (
			 5.14911468e+00 + T * ( -1.36622009e-02 
			+ T * (  4.91453921e-05 + T * ( -4.84246767e-08
			+ T *  1.66603441e-11 ) ) ) );
		h[sC4H6] =  1.53718755e+02 * (
			T * (  4.01336263e+00 + T * (  2.22313425e-03
			+ T * (  2.60227673e-05 + T * ( -2.79185323e-08
			+ T *  9.21507692e-12 ) ) ) ) +  1.14807231e+04 );
		cp[sC4H6] =  1.53718755e+02 * (
			 4.01336263e+00 + T * (  4.44626850e-03 
			+ T * (  7.80683019e-05 + T * ( -1.11674129e-07
			+ T *  4.60753846e-11 ) ) ) );
		h[sTXC5H5O] =  1.02532248e+02 * (
			T * (  2.30436010e-01 + T * (  1.61612860e-02
			+ T * (  9.63363600e-06 + T * ( -1.76701532e-08
			+ T *  6.68143480e-12 ) ) ) ) +  5.55547240e+03 );
		cp[sTXC5H5O] =  1.02532248e+02 * (
			 2.30436010e-01 + T * (  3.23225720e-02 
			+ T * (  2.89009080e-05 + T * ( -7.06806130e-08
			+ T *  3.34071740e-11 ) ) ) );
		h[sC5H4O] =  1.03822832e+02 * (
			T * ( -3.64380971e+00 + T * (  3.07164598e-02
			+ T * ( -1.97383079e-05 + T * (  7.08083390e-09
			+ T * -1.00545227e-12 ) ) ) ) +  5.46809680e+03 );
		cp[sC5H4O] =  1.03822832e+02 * (
			-3.64380971e+00 + T * (  6.14329196e-02 
			+ T * ( -5.92149236e-05 + T * (  2.83233356e-08
			+ T * -5.02726134e-12 ) ) ) );
		h[sSXC5H5O] =  1.02532248e+02 * (
			T * ( -3.07776000e+00 + T * (  2.62908395e-02
			+ T * ( -9.61883767e-06 + T * ( -8.47136975e-10
			+ T *  1.26722798e-12 ) ) ) ) +  2.55104550e+04 );
		cp[sSXC5H5O] =  1.02532248e+02 * (
			-3.07776000e+00 + T * (  5.25816790e-02 
			+ T * ( -2.88565130e-05 + T * ( -3.38854790e-09
			+ T *  6.33613990e-12 ) ) ) );
		h[sC9H8] =  7.15803158e+01 * (
			T * ( -8.12447817e+00 + T * (  4.88828533e-02
			+ T * ( -2.43478658e-05 + T * (  4.70737525e-09
			+ T *  3.68066426e-13 ) ) ) ) +  1.86589996e+04 );
		cp[sC9H8] =  7.15803158e+01 * (
			-8.12447817e+00 + T * (  9.77657067e-02 
			+ T * ( -7.30435974e-05 + T * (  1.88295010e-08
			+ T *  1.84033213e-12 ) ) ) );
		h[sC9H7] =  7.22069373e+01 * (
			T * ( -8.73685384e+00 + T * (  5.17108180e-02
			+ T * ( -3.07807798e-05 + T * (  9.39057395e-09
			+ T * -8.81210540e-13 ) ) ) ) +  3.31641009e+04 );
		cp[sC9H7] =  7.22069373e+01 * (
			-8.73685384e+00 + T * (  1.03421636e-01 
			+ T * ( -9.23423393e-05 + T * (  3.75622958e-08
			+ T * -4.40605270e-12 ) ) ) );
		h[sA1CH2XC7H7] =  9.12400413e+01 * (
			T * ( -6.07053038e+00 + T * (  4.17600754e-02
			+ T * ( -2.47233361e-05 + T * (  7.82884618e-09
			+ T * -8.47341736e-13 ) ) ) ) +  2.35894712e+04 );
		cp[sA1CH2XC7H7] =  9.12400413e+01 * (
			-6.07053038e+00 + T * (  8.35201507e-02 
			+ T * ( -7.41700083e-05 + T * (  3.13153847e-08
			+ T * -4.23670868e-12 ) ) ) );
		h[sC9H6O] =  6.38886413e+01 * (
			T * ( -6.53928778e+00 + T * (  4.84661643e-02
			+ T * ( -2.72566219e-05 + T * (  7.41748685e-09
			+ T * -4.49986784e-13 ) ) ) ) +  6.88883578e+03 );
		cp[sC9H6O] =  6.38886413e+01 * (
			-6.53928778e+00 + T * (  9.69323286e-02 
			+ T * ( -8.17698656e-05 + T * (  2.96699474e-08
			+ T * -2.24993392e-12 ) ) ) );
		h[sOXC6H4] =  1.09266940e+02 * (
			T * ( -3.46229657e+00 + T * (  2.87008288e-02
			+ T * ( -1.64328123e-05 + T * (  4.76701207e-09
			+ T * -3.90861418e-13 ) ) ) ) +  5.25223614e+04 );
		cp[sOXC6H4] =  1.09266940e+02 * (
			-3.46229657e+00 + T * (  5.74016575e-02 
			+ T * ( -4.92984369e-05 + T * (  1.90680483e-08
			+ T * -1.95430709e-12 ) ) ) );
		h[sA1CH3XC7H8] =  9.02418217e+01 * (
			T * ( -4.54072038e+00 + T * (  3.42713573e-02
			+ T * ( -1.19037675e-05 + T * ( -1.04849410e-09
			+ T *  1.48355959e-12 ) ) ) ) +  4.64121087e+03 );
		cp[sA1CH3XC7H8] =  9.02418217e+01 * (
			-4.54072038e+00 + T * (  6.85427145e-02 
			+ T * ( -3.57113024e-05 + T * ( -4.19397642e-09
			+ T *  7.41779795e-12 ) ) ) );
		h[sA1OHXC6H6O] =  8.83489183e+01 * (
			T * ( -3.56571190e+00 + T * (  3.30067717e-02
			+ T * ( -1.30985939e-05 + T * ( -9.05634825e-10
			+ T *  1.72483122e-12 ) ) ) ) -1.31101467e+04 );
		cp[sA1OHXC6H6O] =  8.83489183e+01 * (
			-3.56571190e+00 + T * (  6.60135435e-02 
			+ T * ( -3.92957818e-05 + T * ( -3.62253930e-09
			+ T *  8.62415610e-12 ) ) ) );
		h[sHOA1CH3XC7H8O] =  7.68892300e+01 * (
			T * ( -2.49882920e+00 + T * (  3.45763628e-02
			+ T * ( -9.60285530e-06 + T * ( -4.32281625e-09
			+ T *  2.72653336e-12 ) ) ) ) -1.76336196e+04 );
		cp[sHOA1CH3XC7H8O] =  7.68892300e+01 * (
			-2.49882920e+00 + T * (  6.91527256e-02 
			+ T * ( -2.88085659e-05 + T * ( -1.72912650e-08
			+ T *  1.36326668e-11 ) ) ) );
		h[sOA1CH3XC7H7O] =  7.76127177e+01 * (
			T * ( -3.88641950e+00 + T * (  3.92129808e-02
			+ T * ( -2.06935634e-05 + T * (  5.42917367e-09
			+ T * -3.31557060e-13 ) ) ) ) +  4.71681386e+02 );
		cp[sOA1CH3XC7H7O] =  7.76127177e+01 * (
			-3.88641950e+00 + T * (  7.84259616e-02 
			+ T * ( -6.20806903e-05 + T * (  2.17166947e-08
			+ T * -1.65778530e-12 ) ) ) );
		h[sA1CH2OXC7H7O] =  7.76127177e+01 * (
			T * ( -4.75332952e+00 + T * (  3.92916108e-02
			+ T * ( -1.85706249e-05 + T * (  2.99161712e-09
			+ T *  5.03808950e-13 ) ) ) ) +  1.31155256e+04 );
		cp[sA1CH2OXC7H7O] =  7.76127177e+01 * (
			-4.75332952e+00 + T * (  7.85832217e-02 
			+ T * ( -5.57118748e-05 + T * (  1.19664685e-08
			+ T *  2.51904475e-12 ) ) ) );
		h[sA1CH2OHXC7H8O] =  7.68892300e+01 * (
			T * (  2.85739935e+00 + T * (  1.19385310e-02
			+ T * (  2.80169339e-05 + T * ( -3.34963965e-08
			+ T *  1.13987775e-11 ) ) ) ) -1.37956049e+04 );
		cp[sA1CH2OHXC7H8O] =  7.68892300e+01 * (
			 2.85739935e+00 + T * (  2.38770620e-02 
			+ T * (  8.40508017e-05 + T * ( -1.33985586e-07
			+ T *  5.69938876e-11 ) ) ) );
		h[sA1CHOXC7H6O] =  7.83499501e+01 * (
			T * ( -3.47171048e+00 + T * (  3.46445945e-02
			+ T * ( -1.44201170e-05 + T * (  8.59677740e-10
			+ T *  9.62020522e-13 ) ) ) ) -6.14558774e+03 );
		cp[sA1CHOXC7H6O] =  7.83499501e+01 * (
			-3.47171048e+00 + T * (  6.92891889e-02 
			+ T * ( -4.32603509e-05 + T * (  3.43871096e-09
			+ T *  4.81010261e-12 ) ) ) );
		h[sA1OXC6H5O] =  8.93054780e+01 * (
			T * ( -4.51502441e+00 + T * (  3.51975764e-02
			+ T * ( -1.97518737e-05 + T * (  5.32274058e-09
			+ T * -2.90337352e-13 ) ) ) ) +  5.19173466e+03 );
		cp[sA1OXC6H5O] =  8.93054780e+01 * (
			-4.51502441e+00 + T * (  7.03951529e-02 
			+ T * ( -5.92556211e-05 + T * (  2.12909623e-08
			+ T * -1.45168676e-12 ) ) ) );
		h[sA1CH3YXC7H7] =  9.12400413e+01 * (
			T * ( -3.91657299e+00 + T * (  3.32935047e-02
			+ T * ( -1.33241679e-05 + T * (  5.11701993e-10
			+ T *  9.97118780e-13 ) ) ) ) +  3.54243469e+04 );
		cp[sA1CH3YXC7H7] =  9.12400413e+01 * (
			-3.91657299e+00 + T * (  6.65870094e-02 
			+ T * ( -3.99725038e-05 + T * (  2.04680797e-09
			+ T *  4.98559390e-12 ) ) ) );
		h[sPXC3H4] =  2.07536818e+02 * (
			T * (  1.46175323e+00 + T * (  1.23013301e-02
			+ T * ( -6.34064650e-06 + T * (  2.15090855e-09
			+ T * -3.33458480e-13 ) ) ) ) +  2.09209793e+04 );
		cp[sPXC3H4] =  2.07536818e+02 * (
			 1.46175323e+00 + T * (  2.46026602e-02 
			+ T * ( -1.90219395e-05 + T * (  8.60363422e-09
			+ T * -1.66729240e-12 ) ) ) );
		h[sA1C2H4XC8H9] =  7.90697276e+01 * (
			T * (  7.33299107e-01 + T * (  2.29526579e-02
			+ T * (  1.26085744e-05 + T * ( -2.28091853e-08
			+ T *  8.51179356e-12 ) ) ) ) +  2.61572945e+04 );
		cp[sA1C2H4XC8H9] =  7.90697276e+01 * (
			 7.33299107e-01 + T * (  4.59053158e-02 
			+ T * (  3.78257231e-05 + T * ( -9.12367411e-08
			+ T *  4.25589678e-11 ) ) ) );
		h[sA1C2H5XC8H10] =  7.83189525e+01 * (
			T * (  1.24076722e+00 + T * (  1.79566415e-02
			+ T * (  2.51407491e-05 + T * ( -3.29760752e-08
			+ T *  1.14949361e-11 ) ) ) ) +  1.18391719e+03 );
		cp[sA1C2H5XC8H10] =  7.83189525e+01 * (
			 1.24076722e+00 + T * (  3.59132829e-02 
			+ T * (  7.54222474e-05 + T * ( -1.31904301e-07
			+ T *  5.74746803e-11 ) ) ) );
		h[sC8H9O2] =  6.06213544e+01 * (
			T * ( -3.87929417e+00 + T * (  4.49066033e-02
			+ T * ( -1.64561836e-05 + T * ( -1.01586420e-09
			+ T *  1.91431787e-12 ) ) ) ) +  1.81922098e+04 );
		cp[sC8H9O2] =  6.06213544e+01 * (
			-3.87929417e+00 + T * (  8.98132066e-02 
			+ T * ( -4.93685508e-05 + T * ( -4.06345681e-09
			+ T *  9.57158935e-12 ) ) ) );
		h[sC8H8OOH] =  6.06213544e+01 * (
			T * ( -3.87929417e+00 + T * (  4.49066033e-02
			+ T * ( -1.64561836e-05 + T * ( -1.01586420e-09
			+ T *  1.91431787e-12 ) ) ) ) +  1.81922098e+04 );
		cp[sC8H8OOH] =  6.06213544e+01 * (
			-3.87929417e+00 + T * (  8.98132066e-02 
			+ T * ( -4.93685508e-05 + T * ( -4.06345681e-09
			+ T *  9.57158935e-12 ) ) ) );
		h[sOC8H7OOH] =  5.46478336e+01 * (
			T * ( -1.93576547e+00 + T * (  4.62450358e-02
			+ T * ( -2.03243198e-05 + T * (  1.97462569e-09
			+ T *  1.08205393e-12 ) ) ) ) -2.01714879e+04 );
		cp[sOC8H7OOH] =  5.46478336e+01 * (
			-1.93576547e+00 + T * (  9.24900716e-02 
			+ T * ( -6.09729593e-05 + T * (  7.89850275e-09
			+ T *  5.41026965e-12 ) ) ) );
		h[sH2C2] =  3.19340144e+02 * (
			T * (  3.28154830e+00 + T * (  3.48823955e-03
			+ T * ( -7.95174800e-07 + T * ( -3.02610800e-10
			+ T *  1.96379090e-13 ) ) ) ) +  4.86217940e+04 );
		cp[sH2C2] =  3.19340144e+02 * (
			 3.28154830e+00 + T * (  6.97647910e-03 
			+ T * ( -2.38552440e-06 + T * ( -1.21044320e-09
			+ T *  9.81895450e-13 ) ) ) );
		h[sA1CH3CH3XC8H10] =  7.83189525e+01 * (
			T * ( -3.46066830e+00 + T * (  3.58894658e-02
			+ T * ( -8.52036773e-06 + T * ( -4.37176937e-09
			+ T *  2.45713912e-12 ) ) ) ) +  1.62314629e+02 );
		cp[sA1CH3CH3XC8H10] =  7.83189525e+01 * (
			-3.46066830e+00 + T * (  7.17789316e-02 
			+ T * ( -2.55611032e-05 + T * ( -1.74870775e-08
			+ T *  1.22856956e-11 ) ) ) );
		h[sA1CH3CH2XC8H9] =  7.90697276e+01 * (
			T * ( -5.06171538e+00 + T * (  4.35354132e-02
			+ T * ( -2.15249935e-05 + T * (  4.61621605e-09
			+ T *  9.95743576e-14 ) ) ) ) +  1.90750247e+04 );
		cp[sA1CH3CH2XC8H9] =  7.90697276e+01 * (
			-5.06171538e+00 + T * (  8.70708264e-02 
			+ T * ( -6.45749806e-05 + T * (  1.84648642e-08
			+ T *  4.97871788e-13 ) ) ) );
		h[sA1CH3CHOXC8H8O] =  6.92031229e+01 * (
			T * ( -2.28640538e+00 + T * (  3.60390418e-02
			+ T * ( -1.07859591e-05 + T * ( -2.61287950e-09
			+ T *  1.97213191e-12 ) ) ) ) -1.07389561e+04 );
		cp[sA1CH3CHOXC8H8O] =  6.92031229e+01 * (
			-2.28640538e+00 + T * (  7.20780836e-02 
			+ T * ( -3.23578774e-05 + T * ( -1.04515180e-08
			+ T *  9.86065956e-12 ) ) ) );
		h[sA2CH3XC11H10] =  5.84734510e+01 * (
			T * ( -7.20788596e+00 + T * (  5.46223345e-02
			+ T * ( -2.37103682e-05 + T * (  2.35931837e-09
			+ T *  1.20563270e-12 ) ) ) ) +  1.21406213e+04 );
		cp[sA2CH3XC11H10] =  5.84734510e+01 * (
			-7.20788596e+00 + T * (  1.09244669e-01 
			+ T * ( -7.11311046e-05 + T * (  9.43727348e-09
			+ T *  6.02816350e-12 ) ) ) );
		h[sA1CHOCH2XC8H7O] =  6.97886449e+01 * (
			T * ( -3.91654934e+00 + T * (  4.37771316e-02
			+ T * ( -2.39035516e-05 + T * (  6.44161300e-09
			+ T * -4.00917698e-13 ) ) ) ) +  8.48223476e+03 );
		cp[sA1CHOCH2XC8H7O] =  6.97886449e+01 * (
			-3.91654934e+00 + T * (  8.75542631e-02 
			+ T * ( -7.17106548e-05 + T * (  2.57664520e-08
			+ T * -2.00458849e-12 ) ) ) );
		h[sA1CHOCHOXC8H6O2] =  6.19881009e+01 * (
			T * ( -1.05921174e+00 + T * (  3.62775300e-02
			+ T * ( -1.32481793e-05 + T * ( -7.00683277e-10
			+ T *  1.44399999e-12 ) ) ) ) -2.09983131e+04 );
		cp[sA1CHOCHOXC8H6O2] =  6.19881009e+01 * (
			-1.05921174e+00 + T * (  7.25550600e-02 
			+ T * ( -3.97445378e-05 + T * ( -2.80273311e-09
			+ T *  7.21999995e-12 ) ) ) );
		h[sA2OHXC10H8O] =  5.76727893e+01 * (
			T * ( -2.08768263e+00 + T * (  3.84049753e-02
			+ T * ( -5.11976743e-06 + T * ( -1.01164408e-08
			+ T *  4.67519558e-12 ) ) ) ) -6.29056385e+03 );
		cp[sA2OHXC10H8O] =  5.76727893e+01 * (
			-2.08768263e+00 + T * (  7.68099506e-02 
			+ T * ( -1.53593023e-05 + T * ( -4.04657632e-08
			+ T *  2.33759779e-11 ) ) ) );
		h[sA2CH2XC11H9] =  5.88909351e+01 * (
			T * ( -9.33279584e+00 + T * (  6.19204430e-02
			+ T * ( -3.60778607e-05 + T * (  1.08420439e-08
			+ T * -1.00357549e-12 ) ) ) ) +  3.25141112e+04 );
		cp[sA2CH2XC11H9] =  5.88909351e+01 * (
			-9.33279584e+00 + T * (  1.23840886e-01 
			+ T * ( -1.08233582e-04 + T * (  4.33681755e-08
			+ T * -5.01787746e-12 ) ) ) );
		h[sA2CH2OXC11H9O] =  5.28962604e+01 * (
			T * ( -7.78907165e+00 + T * (  5.91291930e-02
			+ T * ( -2.97536620e-05 + T * (  6.00588557e-09
			+ T *  3.26602332e-13 ) ) ) ) +  2.40306202e+04 );
		cp[sA2CH2OXC11H9O] =  5.28962604e+01 * (
			-7.78907165e+00 + T * (  1.18258386e-01 
			+ T * ( -8.92609860e-05 + T * (  2.40235423e-08
			+ T *  1.63301166e-12 ) ) ) );
		h[sA2CHOXC11H8O] =  5.32376708e+01 * (
			T * ( -4.87929110e+00 + T * (  4.91546815e-02
			+ T * ( -1.51615113e-05 + T * ( -5.04520860e-09
			+ T *  3.54298268e-12 ) ) ) ) +  1.63189711e+03 );
		cp[sA2CHOXC11H8O] =  5.32376708e+01 * (
			-4.87929110e+00 + T * (  9.83093630e-02 
			+ T * ( -4.54845339e-05 + T * ( -2.01808344e-08
			+ T *  1.77149134e-11 ) ) ) );
		h[sA2OXC10H7O] =  5.80788790e+01 * (
			T * ( -1.15176448e+00 + T * (  3.05677256e-02
			+ T * (  1.06717028e-05 + T * ( -2.48571322e-08
			+ T *  9.59980086e-12 ) ) ) ) +  1.14058756e+04 );
		cp[sA2OXC10H7O] =  5.80788790e+01 * (
			-1.15176448e+00 + T * (  6.11354512e-02 
			+ T * (  3.20151083e-05 + T * ( -9.94285290e-08
			+ T *  4.79990043e-11 ) ) ) );
		h[sC4H2] =  1.66100767e+02 * (
			T * (  1.73325212e-01 + T * (  2.26974515e-02
			+ T * ( -2.43374610e-05 + T * (  1.48812934e-08
			+ T * -3.74969432e-12 ) ) ) ) +  5.42239385e+04 );
		cp[sC4H2] =  1.66100767e+02 * (
			 1.73325212e-01 + T * (  4.53949030e-02 
			+ T * ( -7.30123830e-05 + T * (  5.95251736e-08
			+ T * -1.87484716e-11 ) ) ) );
		h[sOC6H4O] =  7.69191059e+01 * (
			T * ( -2.04371804e+00 + T * (  3.30482234e-02
			+ T * ( -1.89325813e-05 + T * (  5.62252577e-09
			+ T * -5.34699342e-13 ) ) ) ) -1.24369410e+04 );
		cp[sOC6H4O] =  7.69191059e+01 * (
			-2.04371804e+00 + T * (  6.60964467e-02 
			+ T * ( -5.67977439e-05 + T * (  2.24901031e-08
			+ T * -2.67349671e-12 ) ) ) );
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
			+ T * -6.01780560e-13 ) ) ) ) +  1.47126330e+04 );
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
		h[sC] =  6.92284763e+02 * (
			T * (  2.55423955e+00 + T * ( -1.60768862e-04
			+ T * (  2.44597415e-07 + T * ( -1.83058722e-10
			+ T *  5.33042892e-14 ) ) ) ) +  8.54438832e+04 );
		cp[sC] =  6.92284763e+02 * (
			 2.55423955e+00 + T * ( -3.21537724e-04 
			+ T * (  7.33792245e-07 + T * ( -7.32234889e-10
			+ T *  2.66521446e-13 ) ) ) );
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
	return 155;
}

int GetNReactions( void )
{
	return rEnd;
}

double MAX_C(double X1, double X2)
{
  return ( (X1 > X2) ? X1 : X2 );
}
double MIN_C(double X1, double X2)
{
  return ( (X2 > X1) ? X1 : X2 );
}
double GetPlogRateCoeff( double temp, double pressure, double lgt, double rt_inv, double *PlogP, double *PlogA, double *PlogB, double *PlogE, int np)
{
	double kR,kR_l,kR_r;
	int i;
	if (pressure <= PlogP[0]) 
		kR = PlogA[0] * exp(PlogB[0]*lgt-PlogE[0]*rt_inv);
	else if (pressure >= PlogP[np-1]) 
		kR = PlogA[np-1] * exp(PlogB[np-1]*lgt-PlogE[np-1]*rt_inv);
	else {
		/* interpolate */
		for (i=0;i<np;i++) {
			if (pressure <= PlogP[i]) 
				break;
		}

		kR_l = PlogA[i-1] * exp(PlogB[i-1]*lgt-PlogE[i-1]*rt_inv); 
		kR_r = PlogA[i] * exp(PlogB[i]*lgt-PlogE[i]*rt_inv); 
		kR = exp(log(MAX_C(kR_l,1.0e-60)) + 
			(log(MAX_C(kR_r,1.0e-60)) - log(MAX_C(kR_l,1.0e-60)))
			/(log(PlogP[i])-log(PlogP[i-1]))*(log(pressure)-log(PlogP[i-1]))); 
	}
  return MIN_C(kR,DBL_MAX);
}
