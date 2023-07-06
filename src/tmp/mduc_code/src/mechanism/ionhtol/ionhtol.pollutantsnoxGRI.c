#include "ionhtol.pollutantsnoxGRI.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
static double GetLindRateCoeff( double temp, double pressure
			, double k0, double kInf, double fc, double conc );

void ComputeProductionRates( double *cdot, double *w, double *k
			, double *c, double *M, double temp, double pressure );
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

	int	nSpec = 158;
	int	nSpecIn = 158;
	double	kTroe0, kTroeInf, fcTroe;
	double	RGAS = 8314.34;
	double	lgt = log( temp );
	double	rt = RGAS * temp;

	M[mM0] = c[sN2] + c[sAR] + c[sC6H5CH3] + c[sHO2]
		 + c[sOH] + c[sH2O] + c[sO2] + c[sC2H4]
		 + c[sH] + c[sC2H3] + c[sH2] + c[sCH2O]
		 + c[sHCO] + c[sC2H6] + c[sC2H5] + c[sCH3O2]
		 + c[sCH3] + c[sCH3O] + c[sCH3O2H] + c[sCH4]
		 + c[sCO] + c[sCO2] + c[sO] + c[sH2O2]
		 + c[sC2H2] + c[sCH2] + c[sCH2OH] + c[sHCCO]
		 + c[sCH2CO] + c[sCH3CO] + c[sAXC3H5] + c[sC3H6]
		 + c[sC3H5O] + c[sCH3OH] + c[sCH3CHO] + c[sNXC3H7]
		 + c[sC5H9] + c[sPXC4H9] + c[s1XC5H10] + c[s1XC5H11]
		 + c[sCH3COCH2O2] + c[sCH3COCH2] + c[sCH3COCH2O2H] + c[sCH3COCH2O]
		 + c[sCH2CHO] + c[sNXC7H16] + c[s1XC7H15] + c[s2XC7H15]
		 + c[s3XC7H15] + c[s4XC7H15] + c[s1XC7H14] + c[s2XC7H15O2]
		 + c[s2X3XC7H14O2H] + c[s2X4XC7H14O2H] + c[s2X4XC7H14OOHO2] + c[sNC7KET24XC7H14O3]
		 + c[sCH] + c[sC3H5XTXC3H5] + c[sIXC3H7] + c[sIXC4H8]
		 + c[sIXC4H7] + c[sTXC4H9] + c[sIXC4H9] + c[sCH3COCH3]
		 + c[sC2H5CO] + c[sNEOXC5H11] + c[sIXC4H7O] + c[sIXC4H8O]
		 + c[sIC4H8OOHXIO2XC4H9O4] + c[sIC4H8OOHXIXC4H9O2] + c[sIC4H8OOHXTO2XC4H9O4] + c[sIC4H8OOHXTXC4H9O2]
		 + c[sIC4KETIIXC4H8O3] + c[sTC4H8OOHXIXC4H9O2] + c[sTXC4H9O] + c[sIXC4H9O2]
		 + c[sTXC4H9O2] + c[sTXC3H6CHO] + c[sTXC3H6OCHO] + c[sIXC3H5CHO]
		 + c[sIXC3H7CHO] + c[sIXC4H7OOH] + c[sIXC3H5CO] + c[sYXC7H15]
		 + c[sYXC7H14] + c[sXXC7H13] + c[sIXC8H18] + c[sAXC8H17]
		 + c[sBXC8H17] + c[sCXC8H17] + c[sDXC8H17] + c[sAXC8H17O2]
		 + c[sAXC8H17O2H] + c[sCXC8H17O2] + c[sCXC8H17O2H] + c[sAXC8H17O]
		 + c[sCXC8H17O] + c[sAC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXBXC8H17O2] + c[sAC8H16OOHXCXC8H17O2]
		 + c[sCC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXAO2XC8H17O4] + c[sAC8H16OOHXBO2XC8H17O4] + c[sAC8H16OOHXCO2XC8H17O4]
		 + c[sCC8H16OOHXAO2XC8H17O4] + c[sIC8KETABXC8H16O3] + c[sC6H5CH2] + c[sC6H5]
		 + c[sC6H6] + c[sC6H5OH] + c[sC6H5O] + c[sC6H5CHO]
		 + c[sC6H5CO] + c[sC4H5XIXC4H5] + c[sC4H6] + c[sC6H5CH2O]
		 + c[sC6H5C2H5] + c[sC6H5CH2OH] + c[sC5H5] + c[sC6H4O2]
		 + c[sC5H6] + c[sC5H4O] + c[sC5H5O] + c[sC4H5XNXC4H5]
		 + c[sCH2CHCHCO] + c[sCHCHCHO] + c[sHCCCO] + c[sC2H5OH]
		 + c[sC2H4OH] + c[sCH3CHOH] + c[sCH3CH2O] + c[sCH3HCO]
		 + c[sCH2HCO] + c[sHOC2H4O2] + c[sC2H] + c[sH2CC]
		 + c[sC3H3] + c[sPXC3H4] + c[sAXC3H4] + c[s1XCH2]
		 + c[sC] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN]
		 + c[sHCNO] + c[sNH3];

	M[mM22] = c[sN2] + c[sAR] + c[sC6H5CH3] + c[sHO2]
		 + c[sOH] + 12 * c[sH2O] + c[sO2] + c[sC2H4]
		 + c[sH] + c[sC2H3] + 2.5 * c[sH2] + c[sCH2O]
		 + c[sHCO] + c[sC2H6] + c[sC2H5] + c[sCH3O2]
		 + c[sCH3] + c[sCH3O] + c[sCH3O2H] + c[sCH4]
		 + 1.9 * c[sCO] + 3.8 * c[sCO2] + c[sO] + c[sH2O2]
		 + c[sC2H2] + c[sCH2] + c[sCH2OH] + c[sHCCO]
		 + c[sCH2CO] + c[sCH3CO] + c[sAXC3H5] + c[sC3H6]
		 + c[sC3H5O] + c[sCH3OH] + c[sCH3CHO] + c[sNXC3H7]
		 + c[sC5H9] + c[sPXC4H9] + c[s1XC5H10] + c[s1XC5H11]
		 + c[sCH3COCH2O2] + c[sCH3COCH2] + c[sCH3COCH2O2H] + c[sCH3COCH2O]
		 + c[sCH2CHO] + c[sNXC7H16] + c[s1XC7H15] + c[s2XC7H15]
		 + c[s3XC7H15] + c[s4XC7H15] + c[s1XC7H14] + c[s2XC7H15O2]
		 + c[s2X3XC7H14O2H] + c[s2X4XC7H14O2H] + c[s2X4XC7H14OOHO2] + c[sNC7KET24XC7H14O3]
		 + c[sCH] + c[sC3H5XTXC3H5] + c[sIXC3H7] + c[sIXC4H8]
		 + c[sIXC4H7] + c[sTXC4H9] + c[sIXC4H9] + c[sCH3COCH3]
		 + c[sC2H5CO] + c[sNEOXC5H11] + c[sIXC4H7O] + c[sIXC4H8O]
		 + c[sIC4H8OOHXIO2XC4H9O4] + c[sIC4H8OOHXIXC4H9O2] + c[sIC4H8OOHXTO2XC4H9O4] + c[sIC4H8OOHXTXC4H9O2]
		 + c[sIC4KETIIXC4H8O3] + c[sTC4H8OOHXIXC4H9O2] + c[sTXC4H9O] + c[sIXC4H9O2]
		 + c[sTXC4H9O2] + c[sTXC3H6CHO] + c[sTXC3H6OCHO] + c[sIXC3H5CHO]
		 + c[sIXC3H7CHO] + c[sIXC4H7OOH] + c[sIXC3H5CO] + c[sYXC7H15]
		 + c[sYXC7H14] + c[sXXC7H13] + c[sIXC8H18] + c[sAXC8H17]
		 + c[sBXC8H17] + c[sCXC8H17] + c[sDXC8H17] + c[sAXC8H17O2]
		 + c[sAXC8H17O2H] + c[sCXC8H17O2] + c[sCXC8H17O2H] + c[sAXC8H17O]
		 + c[sCXC8H17O] + c[sAC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXBXC8H17O2] + c[sAC8H16OOHXCXC8H17O2]
		 + c[sCC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXAO2XC8H17O4] + c[sAC8H16OOHXBO2XC8H17O4] + c[sAC8H16OOHXCO2XC8H17O4]
		 + c[sCC8H16OOHXAO2XC8H17O4] + c[sIC8KETABXC8H16O3] + c[sC6H5CH2] + c[sC6H5]
		 + c[sC6H6] + c[sC6H5OH] + c[sC6H5O] + c[sC6H5CHO]
		 + c[sC6H5CO] + c[sC4H5XIXC4H5] + c[sC4H6] + c[sC6H5CH2O]
		 + c[sC6H5C2H5] + c[sC6H5CH2OH] + c[sC5H5] + c[sC6H4O2]
		 + c[sC5H6] + c[sC5H4O] + c[sC5H5O] + c[sC4H5XNXC4H5]
		 + c[sCH2CHCHCO] + c[sCHCHCHO] + c[sHCCCO] + c[sC2H5OH]
		 + c[sC2H4OH] + c[sCH3CHOH] + c[sCH3CH2O] + c[sCH3HCO]
		 + c[sCH2HCO] + c[sHOC2H4O2] + c[sC2H] + c[sH2CC]
		 + c[sC3H3] + c[sPXC3H4] + c[sAXC3H4] + c[s1XCH2]
		 + c[sC] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN]
		 + c[sHCNO] + c[sNH3];

	M[mM1] = c[sN2] + c[sAR] + c[sC6H5CH3] + c[sHO2]
		 + c[sOH] + 5 * c[sH2O] + c[sO2] + c[sC2H4]
		 + c[sH] + c[sC2H3] + 2 * c[sH2] + c[sCH2O]
		 + c[sHCO] + c[sC2H6] + c[sC2H5] + c[sCH3O2]
		 + c[sCH3] + c[sCH3O] + c[sCH3O2H] + c[sCH4]
		 + 2 * c[sCO] + 3 * c[sCO2] + c[sO] + c[sH2O2]
		 + c[sC2H2] + c[sCH2] + c[sCH2OH] + c[sHCCO]
		 + c[sCH2CO] + c[sCH3CO] + c[sAXC3H5] + c[sC3H6]
		 + c[sC3H5O] + c[sCH3OH] + c[sCH3CHO] + c[sNXC3H7]
		 + c[sC5H9] + c[sPXC4H9] + c[s1XC5H10] + c[s1XC5H11]
		 + c[sCH3COCH2O2] + c[sCH3COCH2] + c[sCH3COCH2O2H] + c[sCH3COCH2O]
		 + c[sCH2CHO] + c[sNXC7H16] + c[s1XC7H15] + c[s2XC7H15]
		 + c[s3XC7H15] + c[s4XC7H15] + c[s1XC7H14] + c[s2XC7H15O2]
		 + c[s2X3XC7H14O2H] + c[s2X4XC7H14O2H] + c[s2X4XC7H14OOHO2] + c[sNC7KET24XC7H14O3]
		 + c[sCH] + c[sC3H5XTXC3H5] + c[sIXC3H7] + c[sIXC4H8]
		 + c[sIXC4H7] + c[sTXC4H9] + c[sIXC4H9] + c[sCH3COCH3]
		 + c[sC2H5CO] + c[sNEOXC5H11] + c[sIXC4H7O] + c[sIXC4H8O]
		 + c[sIC4H8OOHXIO2XC4H9O4] + c[sIC4H8OOHXIXC4H9O2] + c[sIC4H8OOHXTO2XC4H9O4] + c[sIC4H8OOHXTXC4H9O2]
		 + c[sIC4KETIIXC4H8O3] + c[sTC4H8OOHXIXC4H9O2] + c[sTXC4H9O] + c[sIXC4H9O2]
		 + c[sTXC4H9O2] + c[sTXC3H6CHO] + c[sTXC3H6OCHO] + c[sIXC3H5CHO]
		 + c[sIXC3H7CHO] + c[sIXC4H7OOH] + c[sIXC3H5CO] + c[sYXC7H15]
		 + c[sYXC7H14] + c[sXXC7H13] + c[sIXC8H18] + c[sAXC8H17]
		 + c[sBXC8H17] + c[sCXC8H17] + c[sDXC8H17] + c[sAXC8H17O2]
		 + c[sAXC8H17O2H] + c[sCXC8H17O2] + c[sCXC8H17O2H] + c[sAXC8H17O]
		 + c[sCXC8H17O] + c[sAC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXBXC8H17O2] + c[sAC8H16OOHXCXC8H17O2]
		 + c[sCC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXAO2XC8H17O4] + c[sAC8H16OOHXBO2XC8H17O4] + c[sAC8H16OOHXCO2XC8H17O4]
		 + c[sCC8H16OOHXAO2XC8H17O4] + c[sIC8KETABXC8H16O3] + c[sC6H5CH2] + c[sC6H5]
		 + c[sC6H6] + c[sC6H5OH] + c[sC6H5O] + c[sC6H5CHO]
		 + c[sC6H5CO] + c[sC4H5XIXC4H5] + c[sC4H6] + c[sC6H5CH2O]
		 + c[sC6H5C2H5] + c[sC6H5CH2OH] + c[sC5H5] + c[sC6H4O2]
		 + c[sC5H6] + c[sC5H4O] + c[sC5H5O] + c[sC4H5XNXC4H5]
		 + c[sCH2CHCHCO] + c[sCHCHCHO] + c[sHCCCO] + c[sC2H5OH]
		 + c[sC2H4OH] + c[sCH3CHOH] + c[sCH3CH2O] + c[sCH3HCO]
		 + c[sCH2HCO] + c[sHOC2H4O2] + c[sC2H] + c[sH2CC]
		 + c[sC3H3] + c[sPXC3H4] + c[sAXC3H4] + c[s1XCH2]
		 + c[sC] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN]
		 + c[sHCNO] + c[sNH3];

	M[mM3] = c[sN2] + c[sAR] + c[sC6H5CH3] + c[sHO2]
		 + c[sOH] + 16 * c[sH2O] + c[sO2] + c[sC2H4]
		 + c[sH] + c[sC2H3] + 2 * c[sH2] + c[sCH2O]
		 + c[sHCO] + c[sC2H6] + c[sC2H5] + c[sCH3O2]
		 + c[sCH3] + c[sCH3O] + c[sCH3O2H] + c[sCH4]
		 + 2 * c[sCO] + 3 * c[sCO2] + c[sO] + c[sH2O2]
		 + c[sC2H2] + c[sCH2] + c[sCH2OH] + c[sHCCO]
		 + c[sCH2CO] + c[sCH3CO] + c[sAXC3H5] + c[sC3H6]
		 + c[sC3H5O] + c[sCH3OH] + c[sCH3CHO] + c[sNXC3H7]
		 + c[sC5H9] + c[sPXC4H9] + c[s1XC5H10] + c[s1XC5H11]
		 + c[sCH3COCH2O2] + c[sCH3COCH2] + c[sCH3COCH2O2H] + c[sCH3COCH2O]
		 + c[sCH2CHO] + c[sNXC7H16] + c[s1XC7H15] + c[s2XC7H15]
		 + c[s3XC7H15] + c[s4XC7H15] + c[s1XC7H14] + c[s2XC7H15O2]
		 + c[s2X3XC7H14O2H] + c[s2X4XC7H14O2H] + c[s2X4XC7H14OOHO2] + c[sNC7KET24XC7H14O3]
		 + c[sCH] + c[sC3H5XTXC3H5] + c[sIXC3H7] + c[sIXC4H8]
		 + c[sIXC4H7] + c[sTXC4H9] + c[sIXC4H9] + c[sCH3COCH3]
		 + c[sC2H5CO] + c[sNEOXC5H11] + c[sIXC4H7O] + c[sIXC4H8O]
		 + c[sIC4H8OOHXIO2XC4H9O4] + c[sIC4H8OOHXIXC4H9O2] + c[sIC4H8OOHXTO2XC4H9O4] + c[sIC4H8OOHXTXC4H9O2]
		 + c[sIC4KETIIXC4H8O3] + c[sTC4H8OOHXIXC4H9O2] + c[sTXC4H9O] + c[sIXC4H9O2]
		 + c[sTXC4H9O2] + c[sTXC3H6CHO] + c[sTXC3H6OCHO] + c[sIXC3H5CHO]
		 + c[sIXC3H7CHO] + c[sIXC4H7OOH] + c[sIXC3H5CO] + c[sYXC7H15]
		 + c[sYXC7H14] + c[sXXC7H13] + c[sIXC8H18] + c[sAXC8H17]
		 + c[sBXC8H17] + c[sCXC8H17] + c[sDXC8H17] + c[sAXC8H17O2]
		 + c[sAXC8H17O2H] + c[sCXC8H17O2] + c[sCXC8H17O2H] + c[sAXC8H17O]
		 + c[sCXC8H17O] + c[sAC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXBXC8H17O2] + c[sAC8H16OOHXCXC8H17O2]
		 + c[sCC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXAO2XC8H17O4] + c[sAC8H16OOHXBO2XC8H17O4] + c[sAC8H16OOHXCO2XC8H17O4]
		 + c[sCC8H16OOHXAO2XC8H17O4] + c[sIC8KETABXC8H16O3] + c[sC6H5CH2] + c[sC6H5]
		 + c[sC6H6] + c[sC6H5OH] + c[sC6H5O] + c[sC6H5CHO]
		 + c[sC6H5CO] + c[sC4H5XIXC4H5] + c[sC4H6] + c[sC6H5CH2O]
		 + c[sC6H5C2H5] + c[sC6H5CH2OH] + c[sC5H5] + c[sC6H4O2]
		 + c[sC5H6] + c[sC5H4O] + c[sC5H5O] + c[sC4H5XNXC4H5]
		 + c[sCH2CHCHCO] + c[sCHCHCHO] + c[sHCCCO] + c[sC2H5OH]
		 + c[sC2H4OH] + c[sCH3CHOH] + c[sCH3CH2O] + c[sCH3HCO]
		 + c[sCH2HCO] + c[sHOC2H4O2] + c[sC2H] + c[sH2CC]
		 + c[sC3H3] + c[sPXC3H4] + c[sAXC3H4] + c[s1XCH2]
		 + c[sC] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN]
		 + c[sHCNO] + c[sNH3];

	M[mM41] = c[sN2] + c[sAR] + c[sC6H5CH3] + c[sHO2]
		 + c[sOH] + 6 * c[sH2O] + c[sO2] + c[sC2H4]
		 + c[sH] + c[sC2H3] + 2 * c[sH2] + c[sCH2O]
		 + c[sHCO] + 3 * c[sC2H6] + c[sC2H5] + c[sCH3O2]
		 + c[sCH3] + c[sCH3O] + c[sCH3O2H] + 2 * c[sCH4]
		 + 1.5 * c[sCO] + 2 * c[sCO2] + c[sO] + c[sH2O2]
		 + c[sC2H2] + c[sCH2] + c[sCH2OH] + c[sHCCO]
		 + c[sCH2CO] + c[sCH3CO] + c[sAXC3H5] + c[sC3H6]
		 + c[sC3H5O] + c[sCH3OH] + c[sCH3CHO] + c[sNXC3H7]
		 + c[sC5H9] + c[sPXC4H9] + c[s1XC5H10] + c[s1XC5H11]
		 + c[sCH3COCH2O2] + c[sCH3COCH2] + c[sCH3COCH2O2H] + c[sCH3COCH2O]
		 + c[sCH2CHO] + c[sNXC7H16] + c[s1XC7H15] + c[s2XC7H15]
		 + c[s3XC7H15] + c[s4XC7H15] + c[s1XC7H14] + c[s2XC7H15O2]
		 + c[s2X3XC7H14O2H] + c[s2X4XC7H14O2H] + c[s2X4XC7H14OOHO2] + c[sNC7KET24XC7H14O3]
		 + c[sCH] + c[sC3H5XTXC3H5] + c[sIXC3H7] + c[sIXC4H8]
		 + c[sIXC4H7] + c[sTXC4H9] + c[sIXC4H9] + c[sCH3COCH3]
		 + c[sC2H5CO] + c[sNEOXC5H11] + c[sIXC4H7O] + c[sIXC4H8O]
		 + c[sIC4H8OOHXIO2XC4H9O4] + c[sIC4H8OOHXIXC4H9O2] + c[sIC4H8OOHXTO2XC4H9O4] + c[sIC4H8OOHXTXC4H9O2]
		 + c[sIC4KETIIXC4H8O3] + c[sTC4H8OOHXIXC4H9O2] + c[sTXC4H9O] + c[sIXC4H9O2]
		 + c[sTXC4H9O2] + c[sTXC3H6CHO] + c[sTXC3H6OCHO] + c[sIXC3H5CHO]
		 + c[sIXC3H7CHO] + c[sIXC4H7OOH] + c[sIXC3H5CO] + c[sYXC7H15]
		 + c[sYXC7H14] + c[sXXC7H13] + c[sIXC8H18] + c[sAXC8H17]
		 + c[sBXC8H17] + c[sCXC8H17] + c[sDXC8H17] + c[sAXC8H17O2]
		 + c[sAXC8H17O2H] + c[sCXC8H17O2] + c[sCXC8H17O2H] + c[sAXC8H17O]
		 + c[sCXC8H17O] + c[sAC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXBXC8H17O2] + c[sAC8H16OOHXCXC8H17O2]
		 + c[sCC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXAO2XC8H17O4] + c[sAC8H16OOHXBO2XC8H17O4] + c[sAC8H16OOHXCO2XC8H17O4]
		 + c[sCC8H16OOHXAO2XC8H17O4] + c[sIC8KETABXC8H16O3] + c[sC6H5CH2] + c[sC6H5]
		 + c[sC6H6] + c[sC6H5OH] + c[sC6H5O] + c[sC6H5CHO]
		 + c[sC6H5CO] + c[sC4H5XIXC4H5] + c[sC4H6] + c[sC6H5CH2O]
		 + c[sC6H5C2H5] + c[sC6H5CH2OH] + c[sC5H5] + c[sC6H4O2]
		 + c[sC5H6] + c[sC5H4O] + c[sC5H5O] + c[sC4H5XNXC4H5]
		 + c[sCH2CHCHCO] + c[sCHCHCHO] + c[sHCCCO] + c[sC2H5OH]
		 + c[sC2H4OH] + c[sCH3CHOH] + c[sCH3CH2O] + c[sCH3HCO]
		 + c[sCH2HCO] + c[sHOC2H4O2] + c[sC2H] + c[sH2CC]
		 + c[sC3H3] + c[sPXC3H4] + c[sAXC3H4] + c[s1XCH2]
		 + c[sC] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN]
		 + c[sHCNO] + c[sNH3];

	M[mM14] = c[sN2] + c[sAR] + c[sC6H5CH3] + c[sHO2]
		 + c[sOH] + 6 * c[sH2O] + c[sO2] + c[sC2H4]
		 + c[sH] + c[sC2H3] + 2 * c[sH2] + c[sCH2O]
		 + c[sHCO] + c[sC2H6] + c[sC2H5] + c[sCH3O2]
		 + c[sCH3] + c[sCH3O] + c[sCH3O2H] + 2 * c[sCH4]
		 + 1.5 * c[sCO] + 2 * c[sCO2] + c[sO] + c[sH2O2]
		 + c[sC2H2] + c[sCH2] + c[sCH2OH] + c[sHCCO]
		 + c[sCH2CO] + c[sCH3CO] + c[sAXC3H5] + c[sC3H6]
		 + c[sC3H5O] + c[sCH3OH] + c[sCH3CHO] + c[sNXC3H7]
		 + c[sC5H9] + c[sPXC4H9] + c[s1XC5H10] + c[s1XC5H11]
		 + c[sCH3COCH2O2] + c[sCH3COCH2] + c[sCH3COCH2O2H] + c[sCH3COCH2O]
		 + c[sCH2CHO] + c[sNXC7H16] + c[s1XC7H15] + c[s2XC7H15]
		 + c[s3XC7H15] + c[s4XC7H15] + c[s1XC7H14] + c[s2XC7H15O2]
		 + c[s2X3XC7H14O2H] + c[s2X4XC7H14O2H] + c[s2X4XC7H14OOHO2] + c[sNC7KET24XC7H14O3]
		 + c[sCH] + c[sC3H5XTXC3H5] + c[sIXC3H7] + c[sIXC4H8]
		 + c[sIXC4H7] + c[sTXC4H9] + c[sIXC4H9] + c[sCH3COCH3]
		 + c[sC2H5CO] + c[sNEOXC5H11] + c[sIXC4H7O] + c[sIXC4H8O]
		 + c[sIC4H8OOHXIO2XC4H9O4] + c[sIC4H8OOHXIXC4H9O2] + c[sIC4H8OOHXTO2XC4H9O4] + c[sIC4H8OOHXTXC4H9O2]
		 + c[sIC4KETIIXC4H8O3] + c[sTC4H8OOHXIXC4H9O2] + c[sTXC4H9O] + c[sIXC4H9O2]
		 + c[sTXC4H9O2] + c[sTXC3H6CHO] + c[sTXC3H6OCHO] + c[sIXC3H5CHO]
		 + c[sIXC3H7CHO] + c[sIXC4H7OOH] + c[sIXC3H5CO] + c[sYXC7H15]
		 + c[sYXC7H14] + c[sXXC7H13] + c[sIXC8H18] + c[sAXC8H17]
		 + c[sBXC8H17] + c[sCXC8H17] + c[sDXC8H17] + c[sAXC8H17O2]
		 + c[sAXC8H17O2H] + c[sCXC8H17O2] + c[sCXC8H17O2H] + c[sAXC8H17O]
		 + c[sCXC8H17O] + c[sAC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXBXC8H17O2] + c[sAC8H16OOHXCXC8H17O2]
		 + c[sCC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXAO2XC8H17O4] + c[sAC8H16OOHXBO2XC8H17O4] + c[sAC8H16OOHXCO2XC8H17O4]
		 + c[sCC8H16OOHXAO2XC8H17O4] + c[sIC8KETABXC8H16O3] + c[sC6H5CH2] + c[sC6H5]
		 + c[sC6H6] + c[sC6H5OH] + c[sC6H5O] + c[sC6H5CHO]
		 + c[sC6H5CO] + c[sC4H5XIXC4H5] + c[sC4H6] + c[sC6H5CH2O]
		 + c[sC6H5C2H5] + c[sC6H5CH2OH] + c[sC5H5] + c[sC6H4O2]
		 + c[sC5H6] + c[sC5H4O] + c[sC5H5O] + c[sC4H5XNXC4H5]
		 + c[sCH2CHCHCO] + c[sCHCHCHO] + c[sHCCCO] + c[sC2H5OH]
		 + c[sC2H4OH] + c[sCH3CHOH] + c[sCH3CH2O] + c[sCH3HCO]
		 + c[sCH2HCO] + c[sHOC2H4O2] + c[sC2H] + c[sH2CC]
		 + c[sC3H3] + c[sPXC3H4] + c[sAXC3H4] + c[s1XCH2]
		 + c[sC] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN]
		 + c[sHCNO] + c[sNH3];

	M[mM15] = c[sN2] + c[sAR] + c[sC6H5CH3] + c[sHO2]
		 + c[sOH] + 5 * c[sH2O] + c[sO2] + c[sC2H4]
		 + c[sH] + c[sC2H3] + 2 * c[sH2] + c[sCH2O]
		 + c[sHCO] + c[sC2H6] + c[sC2H5] + c[sCH3O2]
		 + c[sCH3] + c[sCH3O] + c[sCH3O2H] + c[sCH4]
		 + 2 * c[sCO] + 3 * c[sCO2] + c[sO] + c[sH2O2]
		 + c[sC2H2] + c[sCH2] + c[sCH2OH] + c[sHCCO]
		 + c[sCH2CO] + c[sCH3CO] + c[sAXC3H5] + c[sC3H6]
		 + c[sC3H5O] + c[sCH3OH] + c[sCH3CHO] + c[sNXC3H7]
		 + c[sC5H9] + c[sPXC4H9] + c[s1XC5H10] + c[s1XC5H11]
		 + c[sCH3COCH2O2] + c[sCH3COCH2] + c[sCH3COCH2O2H] + c[sCH3COCH2O]
		 + c[sCH2CHO] + c[sNXC7H16] + c[s1XC7H15] + c[s2XC7H15]
		 + c[s3XC7H15] + c[s4XC7H15] + c[s1XC7H14] + c[s2XC7H15O2]
		 + c[s2X3XC7H14O2H] + c[s2X4XC7H14O2H] + c[s2X4XC7H14OOHO2] + c[sNC7KET24XC7H14O3]
		 + c[sCH] + c[sC3H5XTXC3H5] + c[sIXC3H7] + c[sIXC4H8]
		 + c[sIXC4H7] + c[sTXC4H9] + c[sIXC4H9] + c[sCH3COCH3]
		 + c[sC2H5CO] + c[sNEOXC5H11] + c[sIXC4H7O] + c[sIXC4H8O]
		 + c[sIC4H8OOHXIO2XC4H9O4] + c[sIC4H8OOHXIXC4H9O2] + c[sIC4H8OOHXTO2XC4H9O4] + c[sIC4H8OOHXTXC4H9O2]
		 + c[sIC4KETIIXC4H8O3] + c[sTC4H8OOHXIXC4H9O2] + c[sTXC4H9O] + c[sIXC4H9O2]
		 + c[sTXC4H9O2] + c[sTXC3H6CHO] + c[sTXC3H6OCHO] + c[sIXC3H5CHO]
		 + c[sIXC3H7CHO] + c[sIXC4H7OOH] + c[sIXC3H5CO] + c[sYXC7H15]
		 + c[sYXC7H14] + c[sXXC7H13] + c[sIXC8H18] + c[sAXC8H17]
		 + c[sBXC8H17] + c[sCXC8H17] + c[sDXC8H17] + c[sAXC8H17O2]
		 + c[sAXC8H17O2H] + c[sCXC8H17O2] + c[sCXC8H17O2H] + c[sAXC8H17O]
		 + c[sCXC8H17O] + c[sAC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXBXC8H17O2] + c[sAC8H16OOHXCXC8H17O2]
		 + c[sCC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXAO2XC8H17O4] + c[sAC8H16OOHXBO2XC8H17O4] + c[sAC8H16OOHXCO2XC8H17O4]
		 + c[sCC8H16OOHXAO2XC8H17O4] + c[sIC8KETABXC8H16O3] + c[sC6H5CH2] + c[sC6H5]
		 + c[sC6H6] + c[sC6H5OH] + c[sC6H5O] + c[sC6H5CHO]
		 + c[sC6H5CO] + c[sC4H5XIXC4H5] + c[sC4H6] + c[sC6H5CH2O]
		 + c[sC6H5C2H5] + c[sC6H5CH2OH] + c[sC5H5] + c[sC6H4O2]
		 + c[sC5H6] + c[sC5H4O] + c[sC5H5O] + c[sC4H5XNXC4H5]
		 + c[sCH2CHCHCO] + c[sCHCHCHO] + c[sHCCCO] + c[sC2H5OH]
		 + c[sC2H4OH] + c[sCH3CHOH] + c[sCH3CH2O] + c[sCH3HCO]
		 + c[sCH2HCO] + c[sHOC2H4O2] + c[sC2H] + c[sH2CC]
		 + c[sC3H3] + c[sPXC3H4] + c[sAXC3H4] + c[s1XCH2]
		 + c[sC] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN]
		 + c[sHCNO] + c[sNH3];

	M[mM16] = c[sN2] + c[sAR] + c[sC6H5CH3] + c[sHO2]
		 + c[sOH] + 5 * c[sH2O] + c[sO2] + c[sC2H4]
		 + c[sH] + c[sC2H3] + 2 * c[sH2] + c[sCH2O]
		 + c[sHCO] + c[sC2H6] + c[sC2H5] + c[sCH3O2]
		 + c[sCH3] + c[sCH3O] + c[sCH3O2H] + c[sCH4]
		 + 2 * c[sCO] + 3 * c[sCO2] + c[sO] + c[sH2O2]
		 + c[sC2H2] + c[sCH2] + c[sCH2OH] + c[sHCCO]
		 + c[sCH2CO] + c[sCH3CO] + c[sAXC3H5] + c[sC3H6]
		 + c[sC3H5O] + c[sCH3OH] + c[sCH3CHO] + c[sNXC3H7]
		 + c[sC5H9] + c[sPXC4H9] + c[s1XC5H10] + c[s1XC5H11]
		 + c[sCH3COCH2O2] + c[sCH3COCH2] + c[sCH3COCH2O2H] + c[sCH3COCH2O]
		 + c[sCH2CHO] + c[sNXC7H16] + c[s1XC7H15] + c[s2XC7H15]
		 + c[s3XC7H15] + c[s4XC7H15] + c[s1XC7H14] + c[s2XC7H15O2]
		 + c[s2X3XC7H14O2H] + c[s2X4XC7H14O2H] + c[s2X4XC7H14OOHO2] + c[sNC7KET24XC7H14O3]
		 + c[sCH] + c[sC3H5XTXC3H5] + c[sIXC3H7] + c[sIXC4H8]
		 + c[sIXC4H7] + c[sTXC4H9] + c[sIXC4H9] + c[sCH3COCH3]
		 + c[sC2H5CO] + c[sNEOXC5H11] + c[sIXC4H7O] + c[sIXC4H8O]
		 + c[sIC4H8OOHXIO2XC4H9O4] + c[sIC4H8OOHXIXC4H9O2] + c[sIC4H8OOHXTO2XC4H9O4] + c[sIC4H8OOHXTXC4H9O2]
		 + c[sIC4KETIIXC4H8O3] + c[sTC4H8OOHXIXC4H9O2] + c[sTXC4H9O] + c[sIXC4H9O2]
		 + c[sTXC4H9O2] + c[sTXC3H6CHO] + c[sTXC3H6OCHO] + c[sIXC3H5CHO]
		 + c[sIXC3H7CHO] + c[sIXC4H7OOH] + c[sIXC3H5CO] + c[sYXC7H15]
		 + c[sYXC7H14] + c[sXXC7H13] + c[sIXC8H18] + c[sAXC8H17]
		 + c[sBXC8H17] + c[sCXC8H17] + c[sDXC8H17] + c[sAXC8H17O2]
		 + c[sAXC8H17O2H] + c[sCXC8H17O2] + c[sCXC8H17O2H] + c[sAXC8H17O]
		 + c[sCXC8H17O] + c[sAC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXBXC8H17O2] + c[sAC8H16OOHXCXC8H17O2]
		 + c[sCC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXAO2XC8H17O4] + c[sAC8H16OOHXBO2XC8H17O4] + c[sAC8H16OOHXCO2XC8H17O4]
		 + c[sCC8H16OOHXAO2XC8H17O4] + c[sIC8KETABXC8H16O3] + c[sC6H5CH2] + c[sC6H5]
		 + c[sC6H6] + c[sC6H5OH] + c[sC6H5O] + c[sC6H5CHO]
		 + c[sC6H5CO] + c[sC4H5XIXC4H5] + c[sC4H6] + c[sC6H5CH2O]
		 + c[sC6H5C2H5] + c[sC6H5CH2OH] + c[sC5H5] + c[sC6H4O2]
		 + c[sC5H6] + c[sC5H4O] + c[sC5H5O] + c[sC4H5XNXC4H5]
		 + c[sCH2CHCHCO] + c[sCHCHCHO] + c[sHCCCO] + c[sC2H5OH]
		 + c[sC2H4OH] + c[sCH3CHOH] + c[sCH3CH2O] + c[sCH3HCO]
		 + c[sCH2HCO] + c[sHOC2H4O2] + c[sC2H] + c[sH2CC]
		 + c[sC3H3] + c[sPXC3H4] + c[sAXC3H4] + c[s1XCH2]
		 + c[sC] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN]
		 + c[sHCNO] + c[sNH3];

	M[mM17] = c[sN2] + c[sAR] + c[sC6H5CH3] + c[sHO2]
		 + c[sOH] + 5 * c[sH2O] + c[sO2] + c[sC2H4]
		 + c[sH] + c[sC2H3] + c[sH2] + c[sCH2O]
		 + c[sHCO] + c[sC2H6] + c[sC2H5] + c[sCH3O2]
		 + c[sCH3] + c[sCH3O] + c[sCH3O2H] + c[sCH4]
		 + c[sCO] + c[sCO2] + c[sO] + c[sH2O2]
		 + c[sC2H2] + c[sCH2] + c[sCH2OH] + c[sHCCO]
		 + c[sCH2CO] + c[sCH3CO] + c[sAXC3H5] + c[sC3H6]
		 + c[sC3H5O] + c[sCH3OH] + c[sCH3CHO] + c[sNXC3H7]
		 + c[sC5H9] + c[sPXC4H9] + c[s1XC5H10] + c[s1XC5H11]
		 + c[sCH3COCH2O2] + c[sCH3COCH2] + c[sCH3COCH2O2H] + c[sCH3COCH2O]
		 + c[sCH2CHO] + c[sNXC7H16] + c[s1XC7H15] + c[s2XC7H15]
		 + c[s3XC7H15] + c[s4XC7H15] + c[s1XC7H14] + c[s2XC7H15O2]
		 + c[s2X3XC7H14O2H] + c[s2X4XC7H14O2H] + c[s2X4XC7H14OOHO2] + c[sNC7KET24XC7H14O3]
		 + c[sCH] + c[sC3H5XTXC3H5] + c[sIXC3H7] + c[sIXC4H8]
		 + c[sIXC4H7] + c[sTXC4H9] + c[sIXC4H9] + c[sCH3COCH3]
		 + c[sC2H5CO] + c[sNEOXC5H11] + c[sIXC4H7O] + c[sIXC4H8O]
		 + c[sIC4H8OOHXIO2XC4H9O4] + c[sIC4H8OOHXIXC4H9O2] + c[sIC4H8OOHXTO2XC4H9O4] + c[sIC4H8OOHXTXC4H9O2]
		 + c[sIC4KETIIXC4H8O3] + c[sTC4H8OOHXIXC4H9O2] + c[sTXC4H9O] + c[sIXC4H9O2]
		 + c[sTXC4H9O2] + c[sTXC3H6CHO] + c[sTXC3H6OCHO] + c[sIXC3H5CHO]
		 + c[sIXC3H7CHO] + c[sIXC4H7OOH] + c[sIXC3H5CO] + c[sYXC7H15]
		 + c[sYXC7H14] + c[sXXC7H13] + c[sIXC8H18] + c[sAXC8H17]
		 + c[sBXC8H17] + c[sCXC8H17] + c[sDXC8H17] + c[sAXC8H17O2]
		 + c[sAXC8H17O2H] + c[sCXC8H17O2] + c[sCXC8H17O2H] + c[sAXC8H17O]
		 + c[sCXC8H17O] + c[sAC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXBXC8H17O2] + c[sAC8H16OOHXCXC8H17O2]
		 + c[sCC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXAO2XC8H17O4] + c[sAC8H16OOHXBO2XC8H17O4] + c[sAC8H16OOHXCO2XC8H17O4]
		 + c[sCC8H16OOHXAO2XC8H17O4] + c[sIC8KETABXC8H16O3] + c[sC6H5CH2] + c[sC6H5]
		 + c[sC6H6] + c[sC6H5OH] + c[sC6H5O] + c[sC6H5CHO]
		 + c[sC6H5CO] + c[sC4H5XIXC4H5] + c[sC4H6] + c[sC6H5CH2O]
		 + c[sC6H5C2H5] + c[sC6H5CH2OH] + c[sC5H5] + c[sC6H4O2]
		 + c[sC5H6] + c[sC5H4O] + c[sC5H5O] + c[sC4H5XNXC4H5]
		 + c[sCH2CHCHCO] + c[sCHCHCHO] + c[sHCCCO] + c[sC2H5OH]
		 + c[sC2H4OH] + c[sCH3CHOH] + c[sCH3CH2O] + c[sCH3HCO]
		 + c[sCH2HCO] + c[sHOC2H4O2] + c[sC2H] + c[sH2CC]
		 + c[sC3H3] + c[sPXC3H4] + c[sAXC3H4] + c[s1XCH2]
		 + c[sC] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN]
		 + c[sHCNO] + c[sNH3];

	M[mM29] = c[sN2] + c[sAR] + c[sC6H5CH3] + c[sHO2]
		 + c[sOH] + 6 * c[sH2O] + c[sO2] + 2.5 * c[sC2H4]
		 + c[sH] + c[sC2H3] + 2 * c[sH2] + c[sCH2O]
		 + c[sHCO] + 3 * c[sC2H6] + c[sC2H5] + c[sCH3O2]
		 + c[sCH3] + c[sCH3O] + c[sCH3O2H] + 2 * c[sCH4]
		 + 1.5 * c[sCO] + 2 * c[sCO2] + c[sO] + c[sH2O2]
		 + 2.5 * c[sC2H2] + c[sCH2] + c[sCH2OH] + c[sHCCO]
		 + c[sCH2CO] + c[sCH3CO] + c[sAXC3H5] + c[sC3H6]
		 + c[sC3H5O] + c[sCH3OH] + c[sCH3CHO] + c[sNXC3H7]
		 + c[sC5H9] + c[sPXC4H9] + c[s1XC5H10] + c[s1XC5H11]
		 + c[sCH3COCH2O2] + c[sCH3COCH2] + c[sCH3COCH2O2H] + c[sCH3COCH2O]
		 + c[sCH2CHO] + c[sNXC7H16] + c[s1XC7H15] + c[s2XC7H15]
		 + c[s3XC7H15] + c[s4XC7H15] + c[s1XC7H14] + c[s2XC7H15O2]
		 + c[s2X3XC7H14O2H] + c[s2X4XC7H14O2H] + c[s2X4XC7H14OOHO2] + c[sNC7KET24XC7H14O3]
		 + c[sCH] + c[sC3H5XTXC3H5] + c[sIXC3H7] + c[sIXC4H8]
		 + c[sIXC4H7] + c[sTXC4H9] + c[sIXC4H9] + c[sCH3COCH3]
		 + c[sC2H5CO] + c[sNEOXC5H11] + c[sIXC4H7O] + c[sIXC4H8O]
		 + c[sIC4H8OOHXIO2XC4H9O4] + c[sIC4H8OOHXIXC4H9O2] + c[sIC4H8OOHXTO2XC4H9O4] + c[sIC4H8OOHXTXC4H9O2]
		 + c[sIC4KETIIXC4H8O3] + c[sTC4H8OOHXIXC4H9O2] + c[sTXC4H9O] + c[sIXC4H9O2]
		 + c[sTXC4H9O2] + c[sTXC3H6CHO] + c[sTXC3H6OCHO] + c[sIXC3H5CHO]
		 + c[sIXC3H7CHO] + c[sIXC4H7OOH] + c[sIXC3H5CO] + c[sYXC7H15]
		 + c[sYXC7H14] + c[sXXC7H13] + c[sIXC8H18] + c[sAXC8H17]
		 + c[sBXC8H17] + c[sCXC8H17] + c[sDXC8H17] + c[sAXC8H17O2]
		 + c[sAXC8H17O2H] + c[sCXC8H17O2] + c[sCXC8H17O2H] + c[sAXC8H17O]
		 + c[sCXC8H17O] + c[sAC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXBXC8H17O2] + c[sAC8H16OOHXCXC8H17O2]
		 + c[sCC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXAO2XC8H17O4] + c[sAC8H16OOHXBO2XC8H17O4] + c[sAC8H16OOHXCO2XC8H17O4]
		 + c[sCC8H16OOHXAO2XC8H17O4] + c[sIC8KETABXC8H16O3] + c[sC6H5CH2] + c[sC6H5]
		 + c[sC6H6] + c[sC6H5OH] + c[sC6H5O] + c[sC6H5CHO]
		 + c[sC6H5CO] + c[sC4H5XIXC4H5] + c[sC4H6] + c[sC6H5CH2O]
		 + c[sC6H5C2H5] + c[sC6H5CH2OH] + c[sC5H5] + c[sC6H4O2]
		 + c[sC5H6] + c[sC5H4O] + c[sC5H5O] + c[sC4H5XNXC4H5]
		 + c[sCH2CHCHCO] + c[sCHCHCHO] + c[sHCCCO] + c[sC2H5OH]
		 + c[sC2H4OH] + c[sCH3CHOH] + c[sCH3CH2O] + c[sCH3HCO]
		 + c[sCH2HCO] + c[sHOC2H4O2] + c[sC2H] + c[sH2CC]
		 + c[sC3H3] + c[sPXC3H4] + c[sAXC3H4] + c[s1XCH2]
		 + c[sC] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN]
		 + c[sHCNO] + c[sNH3];

	M[mM35] = c[sN2] + 0.7 * c[sAR] + c[sC6H5CH3] + c[sHO2]
		 + c[sOH] + 6 * c[sH2O] + c[sO2] + 3 * c[sC2H4]
		 + c[sH] + c[sC2H3] + 2 * c[sH2] + c[sCH2O]
		 + c[sHCO] + 3 * c[sC2H6] + c[sC2H5] + c[sCH3O2]
		 + c[sCH3] + c[sCH3O] + c[sCH3O2H] + 2 * c[sCH4]
		 + 1.5 * c[sCO] + 2 * c[sCO2] + c[sO] + c[sH2O2]
		 + 3 * c[sC2H2] + c[sCH2] + c[sCH2OH] + c[sHCCO]
		 + c[sCH2CO] + c[sCH3CO] + c[sAXC3H5] + c[sC3H6]
		 + c[sC3H5O] + c[sCH3OH] + c[sCH3CHO] + c[sNXC3H7]
		 + c[sC5H9] + c[sPXC4H9] + c[s1XC5H10] + c[s1XC5H11]
		 + c[sCH3COCH2O2] + c[sCH3COCH2] + c[sCH3COCH2O2H] + c[sCH3COCH2O]
		 + c[sCH2CHO] + c[sNXC7H16] + c[s1XC7H15] + c[s2XC7H15]
		 + c[s3XC7H15] + c[s4XC7H15] + c[s1XC7H14] + c[s2XC7H15O2]
		 + c[s2X3XC7H14O2H] + c[s2X4XC7H14O2H] + c[s2X4XC7H14OOHO2] + c[sNC7KET24XC7H14O3]
		 + c[sCH] + c[sC3H5XTXC3H5] + c[sIXC3H7] + c[sIXC4H8]
		 + c[sIXC4H7] + c[sTXC4H9] + c[sIXC4H9] + c[sCH3COCH3]
		 + c[sC2H5CO] + c[sNEOXC5H11] + c[sIXC4H7O] + c[sIXC4H8O]
		 + c[sIC4H8OOHXIO2XC4H9O4] + c[sIC4H8OOHXIXC4H9O2] + c[sIC4H8OOHXTO2XC4H9O4] + c[sIC4H8OOHXTXC4H9O2]
		 + c[sIC4KETIIXC4H8O3] + c[sTC4H8OOHXIXC4H9O2] + c[sTXC4H9O] + c[sIXC4H9O2]
		 + c[sTXC4H9O2] + c[sTXC3H6CHO] + c[sTXC3H6OCHO] + c[sIXC3H5CHO]
		 + c[sIXC3H7CHO] + c[sIXC4H7OOH] + c[sIXC3H5CO] + c[sYXC7H15]
		 + c[sYXC7H14] + c[sXXC7H13] + c[sIXC8H18] + c[sAXC8H17]
		 + c[sBXC8H17] + c[sCXC8H17] + c[sDXC8H17] + c[sAXC8H17O2]
		 + c[sAXC8H17O2H] + c[sCXC8H17O2] + c[sCXC8H17O2H] + c[sAXC8H17O]
		 + c[sCXC8H17O] + c[sAC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXBXC8H17O2] + c[sAC8H16OOHXCXC8H17O2]
		 + c[sCC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXAO2XC8H17O4] + c[sAC8H16OOHXBO2XC8H17O4] + c[sAC8H16OOHXCO2XC8H17O4]
		 + c[sCC8H16OOHXAO2XC8H17O4] + c[sIC8KETABXC8H16O3] + c[sC6H5CH2] + c[sC6H5]
		 + c[sC6H6] + c[sC6H5OH] + c[sC6H5O] + c[sC6H5CHO]
		 + c[sC6H5CO] + c[sC4H5XIXC4H5] + c[sC4H6] + c[sC6H5CH2O]
		 + c[sC6H5C2H5] + c[sC6H5CH2OH] + c[sC5H5] + c[sC6H4O2]
		 + c[sC5H6] + c[sC5H4O] + c[sC5H5O] + c[sC4H5XNXC4H5]
		 + c[sCH2CHCHCO] + c[sCHCHCHO] + c[sHCCCO] + c[sC2H5OH]
		 + c[sC2H4OH] + c[sCH3CHOH] + c[sCH3CH2O] + c[sCH3HCO]
		 + c[sCH2HCO] + c[sHOC2H4O2] + c[sC2H] + c[sH2CC]
		 + c[sC3H3] + c[sPXC3H4] + c[sAXC3H4] + c[s1XCH2]
		 + c[sC] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN]
		 + c[sHCNO] + c[sNH3];

	M[mM36] = c[sN2] + 0.7 * c[sAR] + c[sC6H5CH3] + c[sHO2]
		 + c[sOH] + 6 * c[sH2O] + c[sO2] + 3 * c[sC2H4]
		 + c[sH] + c[sC2H3] + 2 * c[sH2] + c[sCH2O]
		 + c[sHCO] + 3 * c[sC2H6] + c[sC2H5] + c[sCH3O2]
		 + c[sCH3] + c[sCH3O] + c[sCH3O2H] + 2 * c[sCH4]
		 + 1.5 * c[sCO] + 2 * c[sCO2] + c[sO] + c[sH2O2]
		 + 3 * c[sC2H2] + c[sCH2] + c[sCH2OH] + c[sHCCO]
		 + c[sCH2CO] + c[sCH3CO] + c[sAXC3H5] + c[sC3H6]
		 + c[sC3H5O] + c[sCH3OH] + c[sCH3CHO] + c[sNXC3H7]
		 + c[sC5H9] + c[sPXC4H9] + c[s1XC5H10] + c[s1XC5H11]
		 + c[sCH3COCH2O2] + c[sCH3COCH2] + c[sCH3COCH2O2H] + c[sCH3COCH2O]
		 + c[sCH2CHO] + c[sNXC7H16] + c[s1XC7H15] + c[s2XC7H15]
		 + c[s3XC7H15] + c[s4XC7H15] + c[s1XC7H14] + c[s2XC7H15O2]
		 + c[s2X3XC7H14O2H] + c[s2X4XC7H14O2H] + c[s2X4XC7H14OOHO2] + c[sNC7KET24XC7H14O3]
		 + c[sCH] + c[sC3H5XTXC3H5] + c[sIXC3H7] + c[sIXC4H8]
		 + c[sIXC4H7] + c[sTXC4H9] + c[sIXC4H9] + c[sCH3COCH3]
		 + c[sC2H5CO] + c[sNEOXC5H11] + c[sIXC4H7O] + c[sIXC4H8O]
		 + c[sIC4H8OOHXIO2XC4H9O4] + c[sIC4H8OOHXIXC4H9O2] + c[sIC4H8OOHXTO2XC4H9O4] + c[sIC4H8OOHXTXC4H9O2]
		 + c[sIC4KETIIXC4H8O3] + c[sTC4H8OOHXIXC4H9O2] + c[sTXC4H9O] + c[sIXC4H9O2]
		 + c[sTXC4H9O2] + c[sTXC3H6CHO] + c[sTXC3H6OCHO] + c[sIXC3H5CHO]
		 + c[sIXC3H7CHO] + c[sIXC4H7OOH] + c[sIXC3H5CO] + c[sYXC7H15]
		 + c[sYXC7H14] + c[sXXC7H13] + c[sIXC8H18] + c[sAXC8H17]
		 + c[sBXC8H17] + c[sCXC8H17] + c[sDXC8H17] + c[sAXC8H17O2]
		 + c[sAXC8H17O2H] + c[sCXC8H17O2] + c[sCXC8H17O2H] + c[sAXC8H17O]
		 + c[sCXC8H17O] + c[sAC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXBXC8H17O2] + c[sAC8H16OOHXCXC8H17O2]
		 + c[sCC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXAO2XC8H17O4] + c[sAC8H16OOHXBO2XC8H17O4] + c[sAC8H16OOHXCO2XC8H17O4]
		 + c[sCC8H16OOHXAO2XC8H17O4] + c[sIC8KETABXC8H16O3] + c[sC6H5CH2] + c[sC6H5]
		 + c[sC6H6] + c[sC6H5OH] + c[sC6H5O] + c[sC6H5CHO]
		 + c[sC6H5CO] + c[sC4H5XIXC4H5] + c[sC4H6] + c[sC6H5CH2O]
		 + c[sC6H5C2H5] + c[sC6H5CH2OH] + c[sC5H5] + c[sC6H4O2]
		 + c[sC5H6] + c[sC5H4O] + c[sC5H5O] + c[sC4H5XNXC4H5]
		 + c[sCH2CHCHCO] + c[sCHCHCHO] + c[sHCCCO] + c[sC2H5OH]
		 + c[sC2H4OH] + c[sCH3CHOH] + c[sCH3CH2O] + c[sCH3HCO]
		 + c[sCH2HCO] + c[sHOC2H4O2] + c[sC2H] + c[sH2CC]
		 + c[sC3H3] + c[sPXC3H4] + c[sAXC3H4] + c[s1XCH2]
		 + c[sC] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN]
		 + c[sHCNO] + c[sNH3];

	M[mM40] = c[sN2] + c[sAR] + c[sC6H5CH3] + c[sHO2]
		 + c[sOH] + 6 * c[sH2O] + c[sO2] + c[sC2H4]
		 + c[sH] + c[sC2H3] + 2 * c[sH2] + c[sCH2O]
		 + c[sHCO] + 3 * c[sC2H6] + c[sC2H5] + c[sCH3O2]
		 + c[sCH3] + c[sCH3O] + c[sCH3O2H] + 2 * c[sCH4]
		 + 1.5 * c[sCO] + 2 * c[sCO2] + c[sO] + c[sH2O2]
		 + c[sC2H2] + c[sCH2] + c[sCH2OH] + c[sHCCO]
		 + c[sCH2CO] + c[sCH3CO] + c[sAXC3H5] + c[sC3H6]
		 + c[sC3H5O] + c[sCH3OH] + c[sCH3CHO] + c[sNXC3H7]
		 + c[sC5H9] + c[sPXC4H9] + c[s1XC5H10] + c[s1XC5H11]
		 + c[sCH3COCH2O2] + c[sCH3COCH2] + c[sCH3COCH2O2H] + c[sCH3COCH2O]
		 + c[sCH2CHO] + c[sNXC7H16] + c[s1XC7H15] + c[s2XC7H15]
		 + c[s3XC7H15] + c[s4XC7H15] + c[s1XC7H14] + c[s2XC7H15O2]
		 + c[s2X3XC7H14O2H] + c[s2X4XC7H14O2H] + c[s2X4XC7H14OOHO2] + c[sNC7KET24XC7H14O3]
		 + c[sCH] + c[sC3H5XTXC3H5] + c[sIXC3H7] + c[sIXC4H8]
		 + c[sIXC4H7] + c[sTXC4H9] + c[sIXC4H9] + c[sCH3COCH3]
		 + c[sC2H5CO] + c[sNEOXC5H11] + c[sIXC4H7O] + c[sIXC4H8O]
		 + c[sIC4H8OOHXIO2XC4H9O4] + c[sIC4H8OOHXIXC4H9O2] + c[sIC4H8OOHXTO2XC4H9O4] + c[sIC4H8OOHXTXC4H9O2]
		 + c[sIC4KETIIXC4H8O3] + c[sTC4H8OOHXIXC4H9O2] + c[sTXC4H9O] + c[sIXC4H9O2]
		 + c[sTXC4H9O2] + c[sTXC3H6CHO] + c[sTXC3H6OCHO] + c[sIXC3H5CHO]
		 + c[sIXC3H7CHO] + c[sIXC4H7OOH] + c[sIXC3H5CO] + c[sYXC7H15]
		 + c[sYXC7H14] + c[sXXC7H13] + c[sIXC8H18] + c[sAXC8H17]
		 + c[sBXC8H17] + c[sCXC8H17] + c[sDXC8H17] + c[sAXC8H17O2]
		 + c[sAXC8H17O2H] + c[sCXC8H17O2] + c[sCXC8H17O2H] + c[sAXC8H17O]
		 + c[sCXC8H17O] + c[sAC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXBXC8H17O2] + c[sAC8H16OOHXCXC8H17O2]
		 + c[sCC8H16OOHXAXC8H17O2] + c[sAC8H16OOHXAO2XC8H17O4] + c[sAC8H16OOHXBO2XC8H17O4] + c[sAC8H16OOHXCO2XC8H17O4]
		 + c[sCC8H16OOHXAO2XC8H17O4] + c[sIC8KETABXC8H16O3] + c[sC6H5CH2] + c[sC6H5]
		 + c[sC6H6] + c[sC6H5OH] + c[sC6H5O] + c[sC6H5CHO]
		 + c[sC6H5CO] + c[sC4H5XIXC4H5] + c[sC4H6] + c[sC6H5CH2O]
		 + c[sC6H5C2H5] + c[sC6H5CH2OH] + c[sC5H5] + c[sC6H4O2]
		 + c[sC5H6] + c[sC5H4O] + c[sC5H5O] + c[sC4H5XNXC4H5]
		 + c[sCH2CHCHCO] + c[sCHCHCHO] + c[sHCCCO] + c[sC2H5OH]
		 + c[sC2H4OH] + c[sCH3CHOH] + c[sCH3CH2O] + c[sCH3HCO]
		 + c[sCH2HCO] + c[sHOC2H4O2] + c[sC2H] + c[sH2CC]
		 + c[sC3H3] + c[sPXC3H4] + c[sAXC3H4] + c[s1XCH2]
		 + c[sC] + c[sN] + c[sNO] + c[sN2O]
		 + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2]
		 + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN]
		 + c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN]
		 + c[sHCNO] + c[sNH3];



	k[r49] = 1.4500000000E+13 * exp( -1 * lgt );
	k[r61f] = 1.5000000000E+04 * exp( 2 * lgt - 25104000 / rt );
	k[r61b] = 1.7500000000E+03 * exp( 2 * lgt - 21380000 / rt );
	k[r62f] = 2.0200000000E+10 * exp( -24916000 / rt );
	k[r62b] = 1.0150000000E+10 * exp( -84600000 / rt );
	k[r112] = 4.0000000000E+09 * exp( 1046000 / rt );
	kTroe0 = 6.9180000000E+39 * exp( -6.431 * lgt - 448525000 / rt );
	kTroeInf = 8.8510000000E+20 * exp( -1.22 * lgt - 427605000 / rt );
	fcTroe = 1 * exp( -0 / temp );
	k[r113f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM0] );
	kTroe0 = 3.2017114613E+32 * exp( -5.83156 * lgt - 26336853.18 / rt );
	kTroeInf = 4.0963209228E+13 * exp( -0.62056 * lgt - 5416853.181 / rt );
	k[r113b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM0] );
	k[r533] = 1.9000000000E+09 * exp( 5021000 / rt );
	k[r535] = 4.6000000000E+07 * exp( 10878000 / rt );
	kTroe0 = 3.5000000000E+10 * exp( -0.41 * lgt + 4669000 / rt );
	kTroeInf = 1.4750000000E+09 * exp( 0.6 * lgt );
	fcTroe = 0.5 * exp( -temp / 1e-30 ) + 0.5 * exp( -temp / 1e+30 ) + 1 * exp( -0.1 / temp );
	k[ri26f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM22] );
	kTroe0 = 4.9170882993E+13 * exp( -0.423969 * lgt - 196710750.7 / rt );
	kTroeInf = 2.0722014976E+12 * exp( 0.586031 * lgt - 201379750.7 / rt );
	k[ri26b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM22] );
	k[ri5f] = 1.5100000000E-10 * exp( 6 * lgt - 25301000 / rt );
	k[ri5b] = 9.6490000000E-13 * exp( 6.56 * lgt - 42760000 / rt );
	k[ri7f] = 1.4000000000E+02 * exp( 1.95 * lgt + 5636000 / rt );
	k[ri7b] = 1.5680000000E+04 * exp( 1.95 * lgt - 87822000 / rt );
	k[ri8f] = 1.9700000000E+11 * exp( -69203000 / rt );
	k[ri8b] = 1.5550000000E+10 * exp( -1778000 / rt );
	k[ri10f] = 2.9700000000E+03 * exp( 2.02 * lgt - 56066000 / rt );
	k[ri10b] = 3.0130000000E+02 * exp( 2.02 * lgt + 16108000 / rt );
	k[ri13P1f] = 1.0000000000E+09;
	k[ri13P1b] = 1.6850000000E+08 * exp( 0.33 * lgt - 131629000 / rt );
	k[ri17f] = 5.5400000000E-01 * exp( 3.5 * lgt - 21619000 / rt );
	k[ri17b] = 1.3550000000E-04 * exp( 4.06 * lgt - 37058000 / rt );
	k[ri22f] = 1.1000000000E+10;
	k[ri22b] = 4.7800000000E+11 * exp( -0.35 * lgt - 102717000 / rt );
	k[ri23f] = 3.0100000000E+10 * exp( -96232000 / rt );
	k[ri23b] = 6.4350000000E+12 * exp( -0.33 * lgt - 354008000 / rt );
	kTroe0 = 1.3500000000E+18 * exp( -2.788 * lgt - 17535000 / rt );
	kTroeInf = 1.8000000000E+07 * exp( -9975000 / rt );
	fcTroe = 1 * exp( -0 / temp );
	k[ri27f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM22] );
	kTroe0 = 1.0328042105E+27 * exp( -3.75343 * lgt - 552312028.5 / rt );
	kTroeInf = 1.3770722807E+16 * exp( -0.965431 * lgt - 544752028.5 / rt );
	k[ri27b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM22] );
	k[ri29f] = 7.3400000000E+10;
	k[ri29b] = 4.8130000000E+11 * exp( -376560000 / rt );
	k[ri31f] = 6.2830000000E+26 * exp( -3.57 * lgt - 389949000 / rt );
	k[ri31b] = 2.6600000000E+18 * exp( -2.57 * lgt - 1787000 / rt );
	k[ri33f] = 9.3340000000E+05 * exp( 1.5 * lgt - 12452000 / rt );
	k[ri33b] = 7.4530000000E+04 * exp( 1.5 * lgt - 73848000 / rt );
	k[ri36f] = 8.0000000000E+10;
	k[ri36b] = 1.0550000000E+12 * exp( -291332000 / rt );
	k[ri37f] = 1.9950000000E+15 * exp( -1.57 * lgt - 122215000 / rt );
	k[ri37b] = 3.5850000000E+15 * exp( -1.59 * lgt + 6824000 / rt );
	kTroe0 = 2.3440000000E+22 * exp( -2.7 * lgt - 128030000 / rt );
	kTroeInf = 5.4500000000E+13 * exp( -56484000 / rt );
	fcTroe = 1 * exp( -0 / temp );
	k[ri40] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM0] );
	k[ri42f] = 3.2500000000E+10;
	k[ri42b] = 7.8570000000E+11 * exp( -0.33 * lgt - 231752000 / rt );
	k[ri43f] = 2.9740000000E+07 * exp( 0.33 * lgt + 16154000 / rt );
	k[ri43b] = 2.0500000000E+10 * exp( -162967000 / rt );
	k[ri44] = 5.5000000000E+07 * exp( -10142000 / rt );
	k[ri46] = 7.5800000000E+09 * exp( -1715000 / rt );
	k[ri47] = 7.0800000000E+10 * exp( -1255000 / rt );
	k[ri48] = 1.6600000000E+10 * exp( -3431000 / rt );
	k[ri52] = 2.4100000000E+10 * exp( -16610000 / rt );
	k[ri54] = 5.8200000000E-06 * exp( 4.53 * lgt - 27434000 / rt );
	kTroe0 = 2.2540000000E+34 * exp( -7.269 * lgt - 27518000 / rt );
	kTroeInf = 3.1100000000E+08 * exp( 0.58 * lgt - 10832000 / rt );
	fcTroe = 1 * exp( -temp / 675 ) + 1 * exp( -1e+15 / temp );
	k[ri60f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM1] );
	kTroe0 = 1.1636981944E+36 * exp( -6.82121 * lgt - 186153480.3 / rt );
	kTroeInf = 1.6056350420E+10 * exp( 1.02779 * lgt - 169467480.3 / rt );
	k[ri60b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM1] );
	k[ri69] = 6.1200000000E+03 * exp( 2 * lgt - 7950000 / rt );
	k[ri72] = 1.2900000000E+17 * exp( -3.3 * lgt - 1188000 / rt );
	k[ri83] = 3.8100000000E+03 * exp( 2 * lgt - 6866000 / rt );
	kTroe0 = 6.0100000000E+30 * exp( -5.39 * lgt - 151461000 / rt );
	kTroeInf = 2.8000000000E+14 * exp( -0.73 * lgt - 137319000 / rt );
	fcTroe = 0.04 * exp( -temp / 67.6 ) + 0.96 * exp( -temp / 1855 ) + 1 * exp( -7543 / temp );
	k[ri84f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM0] );
	kTroe0 = 2.9231660646E+26 * exp( -4.86324 * lgt - 32941650.37 / rt );
	kTroeInf = 1.3618743728E+10 * exp( -0.203239 * lgt - 18799650.37 / rt );
	k[ri84b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM0] );
	k[ri87] = 1.4300000000E+04 * exp( 2 * lgt - 7950000 / rt );
	k[ri94] = 1.0000000000E+10 * exp( -8368000 / rt );
	k[ri111f] = 3.6100000000E+10;
	k[ri111b] = 5.4460000000E+13 * exp( -1.03 * lgt - 71044000 / rt );
	kTroe0 = 1.2000000000E+12 * exp( -52375000 / rt );
	kTroeInf = 3.0000000000E+12 * exp( -69956000 / rt );
	fcTroe = 1 * exp( -0 / temp );
	k[ri136] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM0] );
	k[ri146] = 2.3970000000E+48 * exp( -9.9 * lgt - 343423000 / rt );
	k[ri159] = 3.1200000000E+03 * exp( 2 * lgt + 1247000 / rt );
	k[ri177] = 1.7300000000E+02 * exp( 2.5 * lgt - 10427000 / rt );
	k[ri521f] = 2.3430000000E+24 * exp( -3.42 * lgt - 127486000 / rt );
	k[ri521b] = 1.0440000000E+20 * exp( -3.3 * lgt );
	k[ri522] = 6.3100000000E+14 * exp( -176983000 / rt );
	k[ri568] = 2.0280000000E+12 * exp( 0.09 * lgt - 98575000 / rt );
	k[ri570] = 7.0000000000E+09 * exp( 4184000 / rt );
	k[ri1455] = 2.4000000000E+08 * exp( 3573000 / rt );
	k[ri1457] = 3.7300000000E+09 * exp( 4238000 / rt );
	kTroe0 = 5.7480000000E+35 * exp( -7.4 * lgt - 2619000 / rt );
	kTroeInf = 5.6490000000E+10 * exp( 0.1 * lgt );
	fcTroe = 0.975 * exp( -temp / 1e-15 ) + 0.025 * exp( -temp / 8000 ) + 1 * exp( -3000 / temp );
	k[r16f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM3] );
	kTroe0 = 1.3090853649E+44 * exp( -8.19715 * lgt - 393901606.2 / rt );
	kTroeInf = 1.2865384876E+19 * exp( -0.697155 * lgt - 391282606.2 / rt );
	k[r16b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM3] );
	k[r63] = 1.5100000000E+04 * exp( 1.91 * lgt - 15631000 / rt );
	k[r65] = 4.0000000000E+09 * exp( -117152000 / rt );
	k[r85] = 1.0000000000E+09;
	k[r91] = 2.8000000000E+10;
	k[r104f] = 1.5330000000E+09 * exp( -332921000 / rt );
	k[r104b] = 3.0000000000E+10;
	k[r135] = 6.6200000000E-03 * exp( 3.7 * lgt - 39748000 / rt );
	k[r137] = 7.0000000000E+15 * exp( -342126000 / rt );
	k[r142] = 1.7000000000E+09 * exp( -44769000 / rt );
	k[r169] = 9.4680000000E+13 * exp( -0.55 * lgt - 121336000 / rt );
	k[r170f] = 4.2680000000E+15 * exp( -0.64 * lgt - 155728000 / rt );
	k[r170b] = 1.6000000000E+10 * exp( -12134000 / rt );
	k[r171] = 2.0000000000E+09 * exp( -20920000 / rt );
	k[r174] = 1.0000000000E+09;
	k[r218] = 2.5000000000E+13 * exp( -125520000 / rt );
	k[r255] = 3.8620000000E+17 * exp( -1.41 * lgt - 126273000 / rt );
	k[r330] = 1.0000000000E+16 * exp( -298738000 / rt );
	k[r331] = 2.8000000000E+10 * exp( -16736000 / rt );
	k[r335] = 6.8000000000E+10 * exp( -12803000 / rt );
	k[r464] = 4.1070000000E+17 * exp( -1.44 * lgt - 127152000 / rt );
	k[r465f] = 5.5720000000E+15 * exp( -0.66 * lgt - 158490000 / rt );
	k[r465b] = 1.6000000000E+10 * exp( -12134000 / rt );
	k[r536] = 1.0000000000E+08;
	k[r538] = 7.6000000000E+07 * exp( 2092000 / rt );
	k[r1014f] = 3.1620000000E+16 * exp( -1 * lgt - 138072000 / rt );
	k[r1014b] = 1.0000000000E+07;
	k[r1017] = 1.0000000000E+09;
	k[r1018] = 3.9810000000E+15 * exp( -167360000 / rt );
	k[r1019] = 1.0000000000E+15 * exp( -92048000 / rt );
	k[r1190f] = 7.4900000000E+16 * exp( -0.47 * lgt - 286813000 / rt );
	k[r1190b] = 1.0000000000E+11;
	k[r1330] = 1.4000000000E+08 * exp( -36819000 / rt );
	k[r1335] = 1.0000000000E+09;
	k[r1400f] = 1.6000000000E+13 * exp( -146440000 / rt );
	k[r1400b] = 1.0000000000E+11 * exp( -66944000 / rt );
	k[r1401f] = 1.8100000000E+08 * exp( -7699000 / rt );
	k[r1401b] = 7.1510000000E+07 * exp( -14163000 / rt );
	k[r1452] = 1.7300000000E+09 * exp( 4226000 / rt );
	k[r2078] = 7.8100000000E+22 * exp( -1.57 * lgt - 368443000 / rt );
	k[r2079] = 8.3040000000E+22 * exp( -1.6 * lgt - 369322000 / rt );
	k[r2080] = 5.6000000000E+04 * exp( 2 * lgt - 32217000 / rt );
	k[r2081] = 1.8200000000E+04 * exp( 2 * lgt - 20920000 / rt );
	k[r2082] = 1.8200000000E+04 * exp( 2 * lgt - 20920000 / rt );
	k[r2083] = 9.1000000000E+03 * exp( 2 * lgt - 20920000 / rt );
	k[r2088] = 1.0500000000E+07 * exp( 0.97 * lgt - 6653000 / rt );
	k[r2089] = 9.4000000000E+04 * exp( 1.61 * lgt - 167000 / rt );
	k[r2090] = 9.4000000000E+04 * exp( 1.61 * lgt - 167000 / rt );
	k[r2091] = 4.7000000000E+04 * exp( 1.61 * lgt - 167000 / rt );
	k[r2092] = 1.0040000000E+10 * exp( -81170000 / rt );
	k[r2093] = 6.8800000000E+09 * exp( -71128000 / rt );
	k[r2094] = 6.8800000000E+09 * exp( -71128000 / rt );
	k[r2095] = 3.4400000000E+09 * exp( -71128000 / rt );
	k[r2100f] = 2.5000000000E+10 * exp( -205016000 / rt );
	k[r2100b] = 2.5000000000E+09;
	k[r2101f] = 4.0000000000E+10 * exp( -199158000 / rt );
	k[r2101b] = 4.0700000000E+09;
	k[r2102f] = 4.0000000000E+10 * exp( -199158000 / rt );
	k[r2102b] = 4.0700000000E+09;
	k[r2103f] = 2.0000000000E+10 * exp( -199158000 / rt );
	k[r2103b] = 4.0000000000E+09;
	k[r2116] = 1.2100000000E+10 * exp( -85479000 / rt );
	k[r2117] = 8.0640000000E+09 * exp( -74057000 / rt );
	k[r2118] = 8.0640000000E+09 * exp( -74057000 / rt );
	k[r2119] = 4.0320000000E+09 * exp( -74057000 / rt );
	k[r2142] = 2.7870000000E+18 * exp( -1.42 * lgt - 127361000 / rt );
	k[r2144] = 3.8850000000E+16 * exp( -0.89 * lgt - 124390000 / rt );
	k[r2147] = 2.7950000000E+11 * exp( 0.19 * lgt - 96274000 / rt );
	k[r2151] = 4.6120000000E+16 * exp( -0.89 * lgt - 126315000 / rt );
	k[r2159f] = 2.0000000000E+11 * exp( -46442000 / rt );
	k[r2159b] = 3.0000000000E+11 * exp( -58994000 / rt );
	k[r2160f] = 2.0000000000E+11 * exp( -75730000 / rt );
	k[r2160b] = 6.0000000000E+11 * exp( -88282000 / rt );
	k[r2161f] = 2.0000000000E+11 * exp( -75730000 / rt );
	k[r2161b] = 2.0000000000E+11 * exp( -75730000 / rt );
	k[r2162f] = 2.0000000000E+11 * exp( -75730000 / rt );
	k[r2162b] = 3.0000000000E+11 * exp( -88282000 / rt );
	k[r2185] = 1.0000000000E+16 * exp( -297064000 / rt );
	k[r2188f] = 7.9940000000E+18 * exp( -1.64 * lgt - 141252000 / rt );
	k[r2188b] = 2.0000000000E+09;
	k[r2189f] = 1.0680000000E+22 * exp( -2.34 * lgt - 149578000 / rt );
	k[r2189b] = 2.0000000000E+09;
	k[r2190f] = 9.8310000000E+21 * exp( -2.34 * lgt - 149578000 / rt );
	k[r2190b] = 2.0000000000E+09;
	k[r2191f] = 2.0430000000E+21 * exp( -2.34 * lgt - 149578000 / rt );
	k[r2191b] = 2.0000000000E+09;
	k[r2212P1f] = 8.9000000000E+12 * exp( -124265000 / rt );
	k[r2212P1b] = 2.7040000000E+13 * exp( -0.54 * lgt - 70835000 / rt );
	k[r2212P2f] = 5.9500000000E+12 * exp( -116734000 / rt );
	k[r2212P2b] = 3.3590000000E+11 * exp( -0.07 * lgt - 73555000 / rt );
	k[r2214f] = 1.0960000000E+12 * exp( -92676000 / rt );
	k[r2214b] = 1.4000000000E+10 * exp( -0.07 * lgt - 49497000 / rt );
	k[r2258f] = 4.6420000000E+19 * exp( -2.46 * lgt - 75145000 / rt );
	k[r2258b] = 7.0000000000E+07 * exp( -32635000 / rt );
	k[r2281] = 5.0000000000E+13 * exp( -106692000 / rt );
	k[r2290f] = 3.9510000000E+22 * exp( -2.49 * lgt - 149871000 / rt );
	k[r2290b] = 2.0000000000E+09;
	k[r2291f] = 3.9510000000E+22 * exp( -2.49 * lgt - 149871000 / rt );
	k[r2291b] = 2.0000000000E+09;
	k[r2309] = 1.2380000000E+11 * exp( -80124000 / rt );
	k[r2327] = 1.0000000000E+16 * exp( -179912000 / rt );
	kTroe0 = 2.4770000000E+27 * exp( -4.76 * lgt - 10209900 / rt );
	kTroeInf = 1.2700000000E+13 * exp( -0.63 * lgt - 1602630 / rt );
	fcTroe = 0.217 * exp( -temp / 74 ) + 0.783 * exp( -temp / 2941 ) + 1 * exp( -6964 / temp );
	k[rnog52f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM41] );
	kTroe0 = 1.8460168359E+32 * exp( -4.68828 * lgt - 450281118.1 / rt );
	kTroeInf = 9.4648420732E+17 * exp( -0.558284 * lgt - 441673848.1 / rt );
	k[rnog52b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM41] );
	k[ri2f] = 1.7270000000E+01 * exp( 3 * lgt - 34409000 / rt );
	k[ri2b] = 6.6100000000E-01 * exp( 3 * lgt - 32401000 / rt );
	k[ri3f] = 1.9300000000E+02 * exp( 2.4 * lgt - 8812000 / rt );
	k[ri3b] = 4.8200000000E-01 * exp( 2.9 * lgt - 62174000 / rt );
	k[ri4f] = 2.1300000000E+03 * exp( 2.21 * lgt - 27112000 / rt );
	k[ri4b] = 3.5570000000E+01 * exp( 2.21 * lgt - 16401000 / rt );
	k[ri6f] = 1.0200000000E+11;
	k[ri6b] = 2.8960000000E+12 * exp( -440157000 / rt );
	k[ri9f] = 5.0800000000E+01 * exp( 2.67 * lgt - 26326000 / rt );
	k[ri9b] = 2.2310000000E+01 * exp( 2.67 * lgt - 17560000 / rt );
	k[ri11f] = 2.1600000000E+05 * exp( 1.51 * lgt - 14351000 / rt );
	k[ri11b] = 9.3520000000E+05 * exp( 1.51 * lgt - 77739000 / rt );
	k[ri12f] = 1.8600000000E+14 * exp( -1 * lgt - 71128000 / rt );
	k[ri12b] = 6.4670000000E+07 * exp( 1849000 / rt );
	k[ri14f] = 1.0200000000E+04 * exp( 1.88 * lgt - 749000 / rt );
	k[ri14b] = 2.8510000000E+05 * exp( 1.05 * lgt - 132926000 / rt );
	kTroe0 = 1.1120000000E+28 * exp( -5 * lgt - 18610000 / rt );
	kTroeInf = 1.0810000000E+09 * exp( 0.45 * lgt - 7623000 / rt );
	fcTroe = 1 * exp( -temp / 95 ) + 1 * exp( -200 / temp );
	k[ri15f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM1] );
	kTroe0 = 2.8169976375E+30 * exp( -4.87602 * lgt - 171547144.8 / rt );
	kTroeInf = 2.7384662285E+11 * exp( 0.573982 * lgt - 160560144.8 / rt );
	k[ri15b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM1] );
	k[ri19f] = 1.2200000000E+27 * exp( -5.76 * lgt - 42258000 / rt );
	k[ri19b] = 1.2590000000E+27 * exp( -5.63 * lgt - 93345000 / rt );
	k[ri20f] = 5.1250000000E+03 * exp( 2.06 * lgt - 3577000 / rt );
	k[ri20b] = 1.0100000000E+04 * exp( 2.06 * lgt - 96148000 / rt );
	k[ri21f] = 1.1300000000E+11 * exp( -32844000 / rt );
	k[ri21b] = 2.0800000000E+10 * exp( -53220000 / rt );
	kTroe0 = 1.1350000000E+30 * exp( -5.246 * lgt - 7134000 / rt );
	kTroeInf = 9.2140000000E+13 * exp( -1.17 * lgt - 2660000 / rt );
	fcTroe = 0.595 * exp( -temp / 1120 ) + 0.405 * exp( -temp / 69.6 ) + 1 * exp( -1e+15 / temp );
	k[ri24f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM1] );
	kTroe0 = 4.8444453978E+41 * exp( -6.75845 * lgt - 391734526.7 / rt );
	kTroeInf = 3.9327506516E+25 * exp( -2.68245 * lgt - 387260526.7 / rt );
	k[ri24b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM1] );
	k[ri25f] = 1.8370000000E+24 * exp( -3 * lgt - 512958000 / rt );
	k[ri25b] = 2.2500000000E+16 * exp( -2 * lgt );
	k[ri28f] = 1.0680000000E-18 * exp( 7.13 * lgt - 55731000 / rt );
	k[ri28b] = 9.4440000000E-18 * exp( 7.13 * lgt - 81755000 / rt );
	k[ri30f] = 3.0200000000E+10;
	k[ri30b] = 8.6970000000E+10 * exp( -367774000 / rt );
	k[ri32f] = 3.4300000000E+06 * exp( 1.18 * lgt + 1870000 / rt );
	k[ri32b] = 1.1860000000E+06 * exp( 1.18 * lgt - 122926000 / rt );
	k[ri34f] = 4.1600000000E+08 * exp( 0.57 * lgt - 11556000 / rt );
	k[ri34b] = 1.4590000000E+07 * exp( 0.57 * lgt - 64183000 / rt );
	k[ri35f] = 2.2500000000E+10 * exp( -17991000 / rt );
	k[ri35b] = 6.7560000000E+11 * exp( -318110000 / rt );
	k[ri38] = 3.6360000000E-09 * exp( 5.42 * lgt - 4176000 / rt );
	k[ri45] = 1.5000000000E+09;
	k[ri50P1f] = 5.9420000000E+14 * exp( -0.66 * lgt - 222380000 / rt );
	k[ri50P1b] = 4.2000000000E+11 * exp( -50124000 / rt );
	kTroe0 = 3.0410000000E+24 * exp( -4.63 * lgt - 8573000 / rt );
	kTroeInf = 1.2360000000E+11 * exp( -0.37 * lgt );
	fcTroe = 0.53 * exp( -temp / 100 ) + 0.47 * exp( -temp / 2000 ) + 1 * exp( -1e+15 / temp );
	k[ri51f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM22] );
	kTroe0 = 3.2583834458E+33 * exp( -5.79849 * lgt - 228020605.5 / rt );
	kTroeInf = 1.3243544686E+20 * exp( -1.53849 * lgt - 219447605.5 / rt );
	k[ri51b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM22] );
	k[ri53f] = 3.4200000000E+08 * exp( -80709000 / rt );
	k[ri53b] = 3.3650000000E+08 * exp( -0.33 * lgt - 10468000 / rt );
	k[ri66] = 7.2800000000E+16 * exp( -2.54 * lgt - 7569000 / rt );
	k[ri74f] = 1.0000000000E+15 * exp( -1.56 * lgt );
	k[ri74b] = 7.0260000000E+14 * exp( -1.56 * lgt - 12510000 / rt );
	k[ri75f] = 1.1300000000E+04 * exp( 2 * lgt - 12552000 / rt );
	k[ri75b] = 3.4370000000E+04 * exp( 2 * lgt - 88492000 / rt );
	k[ri76] = 3.2900000000E+18 * exp( -3.3 * lgt - 12000000 / rt );
	k[ri77f] = 3.3000000000E+10;
	k[ri77b] = 4.4020000000E+10 * exp( -301206000 / rt );
	k[ri88] = 2.1900000000E-07 * exp( 4.5 * lgt + 4184000 / rt );
	k[ri89] = 1.1000000000E+10 * exp( -14226000 / rt );
	k[ri90] = 1.7500000000E+09 * exp( -5648000 / rt );
	k[ri91] = 3.2900000000E+18 * exp( -3.3 * lgt - 12000000 / rt );
	kTroe0 = 3.6000000000E+12 * exp( -247986000 / rt );
	kTroeInf = 3.0000000000E+14 * exp( -296980000 / rt );
	fcTroe = 1 * exp( -0 / temp );
	k[ri92] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM0] );
	k[ri95] = 2.0000000000E+11 * exp( -33472000 / rt );
	k[ri96] = 1.0000000000E+10;
	k[ri98] = 8.0000000000E+10;
	k[ri101] = 1.0100000000E+18 * exp( -3.3 * lgt - 6309000 / rt );
	k[ri147] = 2.7300000000E+62 * exp( -13.28 * lgt - 515469000 / rt );
	k[ri151f] = 2.0100000000E+61 * exp( -13.26 * lgt - 495804000 / rt );
	k[ri151b] = 4.8870000000E+53 * exp( -12.25 * lgt - 117487000 / rt );
	k[ri153] = 1.5800000000E+04 * exp( 1.76 * lgt + 5088000 / rt );
	k[ri173] = 1.1100000000E+03 * exp( 2 * lgt - 6071000 / rt );
	k[ri174] = 5.2400000000E+08 * exp( 0.7 * lgt - 24619000 / rt );
	k[ri176] = 6.0300000000E+07 * exp( 0.7 * lgt - 31932000 / rt );
	k[ri179] = 4.8300000000E+30 * exp( -5.81 * lgt - 77404000 / rt );
	k[ri180f] = 8.5690000000E+18 * exp( -1.57 * lgt - 168783000 / rt );
	k[ri180b] = 1.3000000000E+10 * exp( -6527000 / rt );
	k[ri181] = 2.0000000000E+10;
	k[ri182] = 4.5000000000E+08 * exp( -21004000 / rt );
	k[ri206] = 1.9200000000E+66 * exp( -14.22 * lgt - 535970000 / rt );
	k[ri207f] = 3.0700000000E+55 * exp( -11.49 * lgt - 478231000 / rt );
	k[ri207b] = 3.3000000000E+49 * exp( -11.1 * lgt - 102341000 / rt );
	k[ri210] = 3.3300000000E+04 * exp( 1.76 * lgt - 318000 / rt );
	k[ri211] = 3.4000000000E+02 * exp( 2.5 * lgt - 10427000 / rt );
	k[ri213] = 5.2000000000E+03 * exp( 2 * lgt + 1247000 / rt );
	k[ri214] = 1.5800000000E+04 * exp( 1.76 * lgt + 5088000 / rt );
	k[ri215] = 1.9280000000E+01 * exp( 2.6 * lgt - 58199000 / rt );
	k[ri217] = 1.9280000000E+01 * exp( 2.6 * lgt - 58199000 / rt );
	k[ri259f] = 4.6500000000E+46 * exp( -9.83 * lgt - 230455000 / rt );
	k[ri259b] = 5.8890000000E+41 * exp( -9.42 * lgt - 71044000 / rt );
	k[ri262] = 1.6400000000E+37 * exp( -7.4 * lgt - 161795000 / rt );
	k[ri263] = 4.9800000000E+32 * exp( -6.23 * lgt - 167653000 / rt );
	k[ri308] = 1.2190000000E+23 * exp( -1.99 * lgt - 351247000 / rt );
	k[ri309] = 1.0540000000E+07 * exp( 0.97 * lgt - 6636000 / rt );
	k[ri310] = 5.6280000000E+04 * exp( 2 * lgt - 32217000 / rt );
	k[ri314] = 1.0000000000E+14 * exp( -129704000 / rt );
	k[ri317] = 1.8340000000E+15 * exp( -0.73 * lgt - 54015000 / rt );
	k[ri340f] = 4.8200000000E+10 * exp( -33263000 / rt );
	k[ri340b] = 1.8750000000E+09 * exp( 0.33 * lgt - 101504000 / rt );
	k[ri343f] = 9.0000000000E+10 * exp( -63178000 / rt );
	k[ri343b] = 1.8180000000E+10 * exp( -43514000 / rt );
	k[ri344f] = 3.0000000000E+03 * exp( 2 * lgt - 10460000 / rt );
	k[ri344b] = 2.6230000000E+03 * exp( 2 * lgt - 54225000 / rt );
	k[ri347] = 2.0000000000E+10;
	k[ri493] = 3.0560000000E+17 * exp( -1.2 * lgt - 135101000 / rt );
	k[ri50P2f] = 1.8390000000E+11 * exp( -0.66 * lgt - 165477000 / rt );
	k[ri50P2b] = 1.3000000000E+08 * exp( 6816000 / rt );
	k[ri525] = 1.9900000000E+09 * exp( -48827000 / rt );
	k[ri13P2f] = 5.8000000000E+11 * exp( -39999000 / rt );
	k[ri13P2b] = 9.7710000000E+10 * exp( 0.33 * lgt - 171628000 / rt );
	k[ri540] = 1.4000000000E+13 * exp( -1.61 * lgt - 7782000 / rt );
	k[ri710] = 7.0000000000E+09 * exp( 4184000 / rt );
	k[ri711f] = 1.2900000000E+09 * exp( -55815000 / rt );
	k[ri711b] = 1.0000000000E+09 * exp( -31380000 / rt );
	k[ri776f] = 7.0660000000E+18 * exp( -1.34 * lgt - 150247000 / rt );
	k[ri776b] = 2.2600000000E+09;
	k[ri778f] = 6.2770000000E+21 * exp( -1.85 * lgt - 153302000 / rt );
	k[ri778b] = 1.4100000000E+10;
	k[ri779] = 2.5000000000E+10 * exp( -89538000 / rt );
	k[ri780] = 2.0000000000E+11 * exp( -110458000 / rt );
	k[ri790] = 2.2230000000E+20 * exp( -2.26 * lgt - 94851000 / rt );
	k[ri791] = 7.6120000000E+18 * exp( -1.93 * lgt - 92550000 / rt );
	k[ri793] = 4.0000000000E+11 * exp( -92048000 / rt );
	k[ri794] = 4.0000000000E+11 * exp( -92048000 / rt );
	k[ri797] = 8.4510000000E+15 * exp( -0.68 * lgt - 122047000 / rt );
	k[ri832] = 4.8650000000E+14 * exp( -0.25 * lgt - 63764000 / rt );
	k[ri901f] = 7.5000000000E+10 * exp( -102090000 / rt );
	k[ri901b] = 2.7510000000E+11 * exp( -0.57 * lgt - 39279000 / rt );
	k[ri902f] = 9.0000000000E+11 * exp( -123010000 / rt );
	k[ri902b] = 2.2010000000E+11 * exp( -0.38 * lgt - 59496000 / rt );
	k[ri903f] = 5.0000000000E+10 * exp( -100834000 / rt );
	k[ri903b] = 3.9640000000E+10 * exp( -0.2 * lgt - 57656000 / rt );
	k[ri1193] = 1.0130000000E+18 * exp( -1.45 * lgt - 129913000 / rt );
	k[ri1347] = 3.8100000000E+14 * exp( -1.36 * lgt - 23347000 / rt );
	k[ri1393f] = 1.4010000000E+23 * exp( -2.14 * lgt - 155226000 / rt );
	k[ri1393b] = 1.4100000000E+10;
	k[ri1394] = 7.0000000000E+09 * exp( 4184000 / rt );
	k[ri1396f] = 1.2080000000E+19 * exp( -1.42 * lgt - 150540000 / rt );
	k[ri1396b] = 2.2600000000E+09;
	k[ri1412] = 1.5000000000E+16 * exp( -175728000 / rt );
	k[ri1458] = 7.4700000000E+08 * exp( -59622000 / rt );
	k[ri1460] = 3.7100000000E+22 * exp( -3.96 * lgt - 29468000 / rt );
	k[ri1461] = 9.6400000000E+09;
	k[ri1462] = 3.9800000000E+13 * exp( -40585000 / rt );
	k[ri1463] = 2.8790000000E+16 * exp( -0.63 * lgt - 172716000 / rt );
	k[ri1466] = 2.3040000000E+18 * exp( -0.91 * lgt - 384928000 / rt );
	k[ri1468] = 6.4000000000E+15 * exp( -190581000 / rt );
	k[ri1508] = 3.6750000000E+09 * exp( -5481000 / rt );
	k[ri1510] = 1.0700000000E+10 * exp( 1456000 / rt );
	k[ri1511] = 1.0000000000E+09 * exp( -49873000 / rt );
	k[ri1516] = 4.7800000000E+12 * exp( -114934000 / rt );
	k[ri1517] = 3.0000000000E+07 * exp( -6899000 / rt );
	k[ri1536f] = 6.0000000000E+09 * exp( -166942000 / rt );
	k[ri1536b] = 2.2090000000E+09 * exp( -0.28 * lgt - 126000 / rt );
	k[ri1546] = 5.0000000000E+13 * exp( -121754000 / rt );
	k[ri2720] = 9.9720000000E+14 * exp( -0.52 * lgt - 128700000 / rt );
	k[ri2722f] = 6.7230000000E+15 * exp( -0.59 * lgt - 160373000 / rt );
	k[ri2722b] = 2.6000000000E+10 * exp( -10460000 / rt );
	k[ri2727f] = 3.0000000000E-32 * exp( -12552000 / rt );
	k[ri2727b] = 2.0000000000E-32 * exp( -73220000 / rt );
	k[ri2731f] = 5.3000000000E+03 * exp( 2.2 * lgt );
	k[ri2731b] = 1.8340000000E+04 * exp( 2.14 * lgt - 93094000 / rt );
	k[ri2735f] = 3.0700000000E-01 * exp( 3.2 * lgt + 14644000 / rt );
	k[ri2735b] = 4.6000000000E+00 * exp( 3.14 * lgt - 141838000 / rt );
	k[ri2743] = 8.6820000000E+13 * exp( -0.96 * lgt - 160415000 / rt );
	k[ri3134] = 1.8300000000E+55 * exp( -11.43 * lgt - 512540000 / rt );
	k[ri3135] = 1.8200000000E+55 * exp( -11.55 * lgt - 501243000 / rt );
	k[ri3136] = 3.5600000000E+54 * exp( -11.38 * lgt - 484507000 / rt );
	k[ri3137] = 3.7600000000E+54 * exp( -11.36 * lgt - 509193000 / rt );
	k[ri3138] = 7.5100000000E+48 * exp( -9.38 * lgt - 403045000 / rt );
	k[ri3140] = 1.9400000000E+57 * exp( -11.84 * lgt - 414132000 / rt );
	k[ri3141] = 5.7500000000E+49 * exp( -9.66 * lgt - 410199000 / rt );
	k[ri3142] = 7.3410000000E+02 * exp( 2.77 * lgt - 34087000 / rt );
	k[ri3143] = 5.7360000000E+02 * exp( 2.49 * lgt - 17255000 / rt );
	k[ri3144] = 6.0200000000E+02 * exp( 2.4 * lgt - 10807000 / rt );
	k[ri3145] = 1.8800000000E+02 * exp( 2.75 * lgt - 26276000 / rt );
	k[ri3150] = 2.6300000000E+04 * exp( 1.8 * lgt - 1164000 / rt );
	k[ri3151] = 9.0000000000E+02 * exp( 2 * lgt + 4740000 / rt );
	k[ri3152] = 2.4000000000E+03 * exp( 1.9 * lgt + 6071000 / rt );
	k[ri3153] = 1.7800000000E+04 * exp( 1.8 * lgt - 5987000 / rt );
	k[ri3158] = 2.5200000000E+10 * exp( -85521000 / rt );
	k[ri3159] = 5.6000000000E+09 * exp( -74015000 / rt );
	k[ri3160] = 2.8000000000E+09 * exp( -66986000 / rt );
	k[ri3161] = 1.6800000000E+10 * exp( -85521000 / rt );
	k[ri3166] = 6.3000000000E+10 * exp( -212380000 / rt );
	k[ri3167] = 1.4000000000E+10 * exp( -201711000 / rt );
	k[ri3168] = 7.0000000000E+09 * exp( -192715000 / rt );
	k[ri3169] = 4.2000000000E+10 * exp( -212380000 / rt );
	k[ri3199] = 2.5200000000E+10 * exp( -85521000 / rt );
	k[ri3200] = 5.6000000000E+09 * exp( -74015000 / rt );
	k[ri3201] = 2.8000000000E+09 * exp( -66986000 / rt );
	k[ri3202] = 1.6800000000E+10 * exp( -85521000 / rt );
	k[ri3203] = 1.8140000000E+10 * exp( -85479000 / rt );
	k[ri3205] = 1.8140000000E+10 * exp( -85479000 / rt );
	k[ri3207] = 4.0320000000E+09 * exp( -74057000 / rt );
	k[ri3209] = 4.0320000000E+09 * exp( -74057000 / rt );
	k[ri3211] = 2.0000000000E+09 * exp( -66944000 / rt );
	k[ri3213] = 2.0000000000E+09 * exp( -66944000 / rt );
	k[ri3354] = 1.2100000000E+10 * exp( -85479000 / rt );
	k[ri3356] = 1.2100000000E+10 * exp( -85479000 / rt );
	k[ri3178] = 2.7130000000E+17 * exp( -1.27 * lgt - 124265000 / rt );
	k[ri3179] = 1.9220000000E+22 * exp( -2.66 * lgt - 120750000 / rt );
	k[ri3180] = 1.1530000000E+17 * exp( -1.09 * lgt - 124851000 / rt );
	k[ri3184] = 3.7840000000E+22 * exp( -2.66 * lgt - 136231000 / rt );
	k[ri3189] = 2.6190000000E+19 * exp( -1.75 * lgt - 134348000 / rt );
	k[ri3364f] = 1.3900000000E+11 * exp( -64434000 / rt );
	k[ri3364b] = 4.1630000000E+11 * exp( -67781000 / rt );
	k[ri3365f] = 3.7080000000E+11 * exp( -85354000 / rt );
	k[ri3365b] = 1.8590000000E+10 * exp( 0.58 * lgt - 110416000 / rt );
	k[ri3215f] = 1.8310000000E+18 * exp( -1.43 * lgt - 150582000 / rt );
	k[ri3215b] = 3.7500000000E+08;
	k[ri3216f] = 5.2050000000E+19 * exp( -1.91 * lgt - 161461000 / rt );
	k[ri3216b] = 7.5000000000E+08;
	k[ri3217f] = 8.0650000000E+21 * exp( -2.15 * lgt - 148155000 / rt );
	k[ri3217b] = 1.5000000000E+09;
	k[ri3218f] = 3.0530000000E+17 * exp( -1.43 * lgt - 147235000 / rt );
	k[ri3218b] = 6.2500000000E+08;
	k[ri3257] = 1.2600000000E+16 * exp( -177820000 / rt );
	k[ri3259] = 7.5000000000E+15 * exp( -173803000 / rt );
	k[ri3261] = 1.9700000000E+21 * exp( -2.26 * lgt - 62007000 / rt );
	k[ri3263] = 9.6240000000E+19 * exp( -1.73 * lgt - 53681000 / rt );
	k[ri3265f] = 1.8750000000E+10 * exp( -102090000 / rt );
	k[ri3265b] = 4.6260000000E+10 * exp( -0.51 * lgt - 42426000 / rt );
	k[ri3266f] = 6.2500000000E+09 * exp( -87236000 / rt );
	k[ri3266b] = 5.1450000000E+08 * exp( -0.07 * lgt - 34309000 / rt );
	k[ri3267f] = 3.9080000000E+08 * exp( -71337000 / rt );
	k[ri3267b] = 1.6370000000E+07 * exp( 0.07 * lgt - 32552000 / rt );
	k[ri3272P1f] = 3.5150000000E+09 * exp( -93512000 / rt );
	k[ri3272P1b] = 2.9110000000E+09 * exp( -0.51 * lgt - 33890000 / rt );
	k[ri3272P2] = 5.0000000000E+10 * exp( -112340000 / rt );
	k[ri3272P3] = 1.5000000000E+11 * exp( -123010000 / rt );
	k[ri3297f] = 3.4690000000E+19 * exp( -1.81 * lgt - 148825000 / rt );
	k[ri3297b] = 3.7500000000E+08;
	k[ri3298f] = 1.8820000000E+22 * exp( -2.37 * lgt - 160080000 / rt );
	k[ri3298b] = 7.5000000000E+08;
	k[ri3299f] = 5.6580000000E+23 * exp( -2.59 * lgt - 150833000 / rt );
	k[ri3299b] = 1.0000000000E+09;
	k[ri3306f] = 8.4230000000E+19 * exp( -1.88 * lgt - 149034000 / rt );
	k[ri3306b] = 3.7500000000E+08;
	k[ri3311] = 6.2500000000E+09 * exp( -89538000 / rt );
	k[ri3312] = 6.2500000000E+09 * exp( -89538000 / rt );
	k[ri3313] = 2.8130000000E+08 * exp( -80960000 / rt );
	k[ri3367] = 1.0500000000E+16 * exp( -174054000 / rt );
	k[ri3377] = 9.0870000000E+17 * exp( -1.26 * lgt - 119579000 / rt );
	k[ri3378] = 1.2520000000E+17 * exp( -1.08 * lgt - 118031000 / rt );
	k[ri3498] = 1.0290000000E+24 * exp( -3.08 * lgt - 112299000 / rt );
	k[ri3506] = 1.0370000000E+26 * exp( -4.08 * lgt - 108617000 / rt );
	k[r1021102f] = 3.1000000000E+15 * exp( -372999000 / rt );
	k[r1021102b] = 1.6779496382E+12 * exp( -0.48179 * lgt + 5544325.333 / rt );
	k[r1021103] = 1.3000000000E+15 * exp( -372999000 / rt );
	k[r1021104f] = 3.0000000000E+11 * exp( -180000000 / rt );
	k[r1021104b] = 2.2812792256E+11 * exp( -0.49576 * lgt - 2836425.38 / rt );
	k[r1021105] = 1.0200000000E+01 * exp( 2.5 * lgt - 51628000 / rt );
	k[r1021106] = 6.1200000000E+10 * exp( -34455000 / rt );
	k[r1021107] = 6.1200000000E+09 * exp( -21539000 / rt );
	k[r1021108] = 1.2600000000E+10 * exp( -10795000 / rt );
	k[r1021109] = 6.3000000000E+08;
	k[r1021110] = 3.1600000000E+08 * exp( -39748000 / rt );
	k[r1021111] = 2.7000000000E-04 * exp( 4.02 * lgt + 6653000 / rt );
	k[r1021112f] = 1.0500000000E+08 * exp( -39748000 / rt );
	k[r1021112b] = 5.6599349003E+06 * exp( 0.693323 * lgt - 57166252.31 / rt );
	k[r1021113] = 2.7700000000E+00 * exp( 2.81 * lgt - 24154000 / rt );
	k[r1021114f] = 6.0000000000E+09;
	k[r1021114b] = 6.6868969779E+10 * exp( 0.131025 * lgt - 2500528.728 / rt );
	k[r1021115] = 1.1200000000E+08 * exp( -0.33 * lgt - 76567000 / rt );
	k[r1021116] = 5.0000000000E+09 * exp( -1.07 * lgt - 45355000 / rt );
	k[r1021117] = 2.5000000000E+11;
	k[r1021118] = 8.0000000000E+10;
	k[r1021119] = 1.1500000000E+10;
	k[r1021121] = 2.8600000000E+13 * exp( -1.41 * lgt - 20066000 / rt );
	k[r1021124f] = 2.0000000000E+15 * exp( -304177000 / rt );
	k[r1021124b] = 3.8462273852E+04 * exp( 1.32211 * lgt + 27937312.59 / rt );
	k[r1021125f] = 6.0000000000E+10;
	k[r1021125b] = 7.7036740770E+17 * exp( -0.402997 * lgt - 353330463.1 / rt );
	k[r1021126f] = 2.0000000000E+11 * exp( -173218000 / rt );
	k[r1021126b] = 1.3541974923E+10 * exp( 0.0041436 * lgt + 59734778.98 / rt );
	k[r1021127] = 8.4300000000E+09 * exp( -10807000 / rt );
	k[r1021128] = 8.0000000000E+10 * exp( -34455000 / rt );
	k[r1021133] = 1.0200000000E+10 * exp( -162967000 / rt );
	k[r1021134] = 1.7100000000E+06 * exp( 1.18 * lgt + 1870000 / rt );
	k[r1021135] = 5.0000000000E+10 * exp( -20619000 / rt );
	k[r1021136] = 1.2000000000E+10 * exp( -21539000 / rt );
	k[r1021137] = 9.0400000000E+09 * exp( -12887000 / rt );
	k[r1021138] = 2.7700000000E+00 * exp( 2.81 * lgt - 24154000 / rt );
	k[r1021140f] = 8.5100000000E+08 * exp( -12527000 / rt );
	k[r1021140b] = 5.6295794657E+19 * exp( -1.91033 * lgt - 125636063.8 / rt );
	kTroe0 = 6.6000000000E+69 * exp( -16.3 * lgt - 29288000 / rt );
	kTroeInf = 1.0000000000E+11;
	fcTroe = 1 * exp( -temp / 584.9 ) + 1 * exp( -6113 / temp );
	k[r1021142] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM14] );
	k[r1021143f] = 6.3000000000E+10 * exp( -251040000 / rt );
	k[r1021143b] = 2.3005872575E+09 * exp( 1.02521 * lgt - 264314015.8 / rt );
	k[r1021144] = 2.1100000000E+10 * exp( -19121000 / rt );
	k[r1021146] = 1.1000000000E+10 * exp( -18954000 / rt );
	k[r1021147] = 1.1000000000E+10 * exp( -18954000 / rt );
	k[r1021150] = 8.5500000000E+01 * exp( 2.19 * lgt - 159000 / rt );
	k[r1021153] = 3.0000000000E+10;
	k[r1021154f] = 2.6000000000E+10 * exp( -25606000 / rt );
	k[r1021154b] = 9.6191587410E+14 * exp( -1.13132 * lgt - 59819021.27 / rt );
	k[r1021155] = 3.0000000000E+10 * exp( -37572000 / rt );
	k[r1021156] = 3.7600000000E+54 * exp( -12.06 * lgt - 310871000 / rt );
	k[r1021157f] = 1.0000000000E+10 * exp( -158992000 / rt );
	k[r1021157b] = 4.0990133035E+08 * exp( 0.197564 * lgt + 753322.3065 / rt );
	k[r1021158] = 1.0000000000E+10;
	k[r1021159] = 8.5000000000E+10;
	k[r1021160] = 1.0000000000E+10;
	k[r1021161] = 8.0000000000E+09;
	kTroe0 = 1.0000000000E+88 * exp( -21.84 * lgt - 58074000 / rt );
	kTroeInf = 2.5000000000E+11;
	fcTroe = 0.957 * exp( -temp / 304.2 ) + 0.043 * exp( -temp / 60000 ) + 1 * exp( -5896 / temp );
	k[r1021162] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM0] );
	k[r1021163] = 7.4000000000E+11 * exp( -246856000 / rt );
	k[r1021164] = 3.0000000000E+10 * exp( -20920000 / rt );
	k[r1021165] = 4.3000000000E+06 * exp( 1.45 * lgt - 16318000 / rt );
	k[r1021166] = 2.9500000000E+03 * exp( 2 * lgt + 5489000 / rt );
	k[r1021167f] = 1.3000000000E+10 * exp( -44350000 / rt );
	k[r1021167b] = 2.1216529261E+17 * exp( -0.742015 * lgt - 324513267.1 / rt );
	k[r1021168] = 1.1500000000E+11 * exp( -51882000 / rt );
	k[r1021169] = 2.8100000000E+10 * exp( -30761000 / rt );
	k[r1021170] = 6.0000000000E+09;
	k[r1021171f] = 6.0000000000E+09;
	k[r1021171b] = 7.2403279291E+06 * exp( 1.53149 * lgt - 67040441.46 / rt );
	k[r1021172] = 4.9100000000E+09 * exp( -18410000 / rt );
	k[r1021173f] = 6.0000000000E+09;
	k[r1021173b] = 3.6045144362E+09 * exp( 0.824348 * lgt - 19918781.04 / rt );
	k[r1021174f] = 1.5000000000E+11;
	k[r1021174b] = 6.9084532247E+12 * exp( 0.837446 * lgt - 344098705.7 / rt );
	k[r1021175f] = 4.0000000000E+10 * exp( -155436000 / rt );
	k[r1021175b] = 1.2201420424E+12 * exp( -0.851415 * lgt - 12717044.97 / rt );
	k[r1021177] = 3.0800000000E+03 * exp( 2 * lgt );
	k[r1021178] = 7.2000000000E+10 * exp( -14644000 / rt );
	k[r1021179] = 4.8000000000E+01 * exp( 2.71 * lgt - 4602000 / rt );
	k[r1021180f] = 1.2000000000E-04 * exp( 4 * lgt );
	k[r1021180b] = 3.9412816925E+00 * exp( 2.80907 * lgt - 108879200.6 / rt );
	k[r1021182] = 3.1600000000E+08 * exp( -33472000 / rt );
	k[r1021185f] = 1.5000000000E-04 * exp( 4 * lgt );
	k[r1021185b] = 6.0170617535E-03 * exp( 3.64434 * lgt - 34444619.59 / rt );
	k[r1021186] = 3.2000000000E+10 * exp( -0.17 * lgt - 1841000 / rt );
	k[r1021187] = 5.8000000000E+10 * exp( -0.02 * lgt - 84000 / rt );
	k[r1021191] = 6.3000000000E+26 * exp( -4.69 * lgt - 48744000 / rt );
	k[r1021192] = 8.0000000000E+14 * exp( -2.48 * lgt - 45898000 / rt );
	k[r1021193] = 6.6000000000E+10 * exp( -0.02 * lgt - 11464000 / rt );
	k[r1021194] = 3.0000000000E+05 * exp( 1.45 * lgt + 3598000 / rt );
	k[r1021195] = 3.0000000000E+09;
	k[r1021205] = 2.9000000000E+32 * exp( -6.5 * lgt - 88784000 / rt );
	k[r1021206] = 1.1000000000E+79 * exp( -19.62 * lgt - 277190000 / rt );
	k[r1021208] = 6.2000000000E+41 * exp( -7.87 * lgt - 412961000 / rt );
	k[r1021209] = 2.1000000000E+58 * exp( -13.27 * lgt - 170749000 / rt );
	k[r1021210] = 3.9800000000E+11 * exp( 0.7 * lgt - 176816000 / rt );
	k[r1021211] = 6.0000000000E+09;
	k[r1021212f] = 2.0000000000E+13 * exp( -223844000 / rt );
	k[r1021212b] = 9.9567712169E+15 * exp( -0.707144 * lgt - 176722339.6 / rt );
	k[r1021218] = 1.0000000000E+09;
	k[r1021244] = 1.4000000000E+06 * exp( 1 * lgt );
	k[r1021245f] = 1.0000000000E+11;
	k[r1021245b] = 8.1083168944E+05 * exp( 1.76859 * lgt - 119229998.5 / rt );
	k[r1021246] = 3.0000000000E+09;
	k[r1021247] = 1.0000000000E+14 * exp( -138072000 / rt );
	kTroe0 = 2.8800000000E+82 * exp( -18.9 * lgt - 459863000 / rt );
	kTroeInf = 5.9400000000E+23 * exp( -1.68 * lgt - 381426000 / rt );
	fcTroe = 0.5 * exp( -temp / 200 ) + 0.5 * exp( -temp / 890 ) + 1 * exp( -4600 / temp );
	k[reth139f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM15] );
	kTroe0 = 3.5650422488E+71 * exp( -17.6231 * lgt - 89996936.09 / rt );
	kTroeInf = 7.3528996382E+12 * exp( -0.403132 * lgt - 11559936.09 / rt );
	k[reth139b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM15] );
	kTroe0 = 3.2520000000E+82 * exp( -18.81 * lgt - 480867000 / rt );
	kTroeInf = 1.2500000000E+23 * exp( -1.54 * lgt - 401685000 / rt );
	fcTroe = 0.5 * exp( -temp / 300 ) + 0.5 * exp( -temp / 900 ) + 1 * exp( -5000 / temp );
	k[reth140f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM16] );
	kTroe0 = 4.0032767005E+71 * exp( -17.3297 * lgt - 84762598.33 / rt );
	kTroeInf = 1.5387748695E+12 * exp( -0.0597353 * lgt - 5580598.325 / rt );
	k[reth140b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM16] );
	kTroe0 = 2.5700000000E+80 * exp( -18.85 * lgt - 361715000 / rt );
	kTroeInf = 2.7900000000E+13 * exp( 0.09 * lgt - 276713000 / rt );
	fcTroe = 0.3 * exp( -temp / 350 ) + 0.7 * exp( -temp / 800 ) + 1 * exp( -3800 / temp );
	k[reth141] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM17] );
	k[reth143f] = 1.7400000000E+08 * exp( 0.27 * lgt - 2510000 / rt );
	k[reth143b] = 5.2909480971E+06 * exp( 0.541485 * lgt - 80644708.51 / rt );
	k[reth144] = 4.6400000000E+08 * exp( 0.15 * lgt );
	k[reth145] = 7.4600000000E+08 * exp( 0.3 * lgt - 6837000 / rt );
	k[reth146f] = 1.2300000000E+04 * exp( 1.8 * lgt - 21330000 / rt );
	k[reth146b] = 3.6988486123E+01 * exp( 2.16849 * lgt - 34386627.52 / rt );
	k[reth147] = 2.5800000000E+04 * exp( 1.65 * lgt - 11828000 / rt );
	k[reth148f] = 1.5000000000E+04 * exp( 1.6 * lgt - 12711000 / rt );
	k[reth148b] = 2.0986284902E+03 * exp( 1.91153 * lgt - 10554907.16 / rt );
	k[reth149f] = 9.4100000000E+04 * exp( 1.7 * lgt - 22840000 / rt );
	k[reth149b] = 1.7069672358E+02 * exp( 2.03475 * lgt - 27708935.83 / rt );
	k[reth150] = 1.8800000000E+04 * exp( 1.85 * lgt - 7632000 / rt );
	k[reth151f] = 1.5800000000E+04 * exp( 2 * lgt - 18610000 / rt );
	k[reth151b] = 1.3334462327E+03 * exp( 2.27779 * lgt - 8266215.47 / rt );
	k[reth152f] = 2.1900000000E-01 * exp( 3.18 * lgt - 40258000 / rt );
	k[reth152b] = 2.1303811690E-01 * exp( 3.26932 * lgt - 60260657.14 / rt );
	k[reth153] = 7.2800000000E-01 * exp( 2.99 * lgt - 33254000 / rt );
	k[reth154f] = 1.4500000000E-01 * exp( 2.99 * lgt - 32003000 / rt );
	k[reth154b] = 6.5624206032E+00 * exp( 3.02235 * lgt - 36792936.78 / rt );
	k[reth155] = 8.2000000000E+00 * exp( 2.55 * lgt - 44978000 / rt );
	k[reth156f] = 1.2300000000E+01 * exp( 2.55 * lgt - 65898000 / rt );
	k[reth156b] = 3.4366718767E+01 * exp( 2.16857 * lgt - 15903882.86 / rt );
	k[reth157] = 2.5000000000E+09 * exp( -100416000 / rt );
	k[reth158] = 1.1600000000E+32 * exp( -5.89 * lgt - 105746000 / rt );
	k[reth159] = 1.3500000000E+35 * exp( -6.96 * lgt - 99579000 / rt );
	k[reth160f] = 4.0000000000E+07 * exp( -4602000 / rt );
	k[reth160b] = 6.1340666442E+04 * exp( 0.435067 * lgt - 134918351.3 / rt );
	k[reth161f] = 4.6800000000E-01 * exp( 3.16 * lgt - 22510000 / rt );
	k[reth161b] = 5.2224793874E-02 * exp( 3.39705 * lgt - 171526411.3 / rt );
	k[reth163] = 3.0000000000E+10;
	k[reth164] = 1.0000000000E+10;
	k[reth165] = 4.8200000000E+11 * exp( -20991000 / rt );
	k[reth167] = 1.0000000000E+10;
	k[reth169] = 3.0000000000E+10;
	k[reth171] = 4.0000000000E+10;
	k[reth172f] = 5.0000000000E+09;
	k[reth172b] = 1.9211183859E+08 * exp( 0.360418 * lgt - 406771572 / rt );
	k[reth173] = 1.0000000000E+11 * exp( -104600000 / rt );
	k[reth174] = 9.2400000000E+03 * exp( 1.5 * lgt + 4025000 / rt );
	k[reth175f] = 1.7200000000E+02 * exp( 2.4 * lgt - 3410000 / rt );
	k[reth175b] = 1.8009678112E+00 * exp( 2.84673 * lgt - 91814102.05 / rt );
	k[reth177] = 1.7700000000E+15 * exp( -1.9 * lgt - 12447000 / rt );
	k[reth179] = 4.6600000000E+10 * exp( -0.35 * lgt - 12502000 / rt );
	k[reth180f] = 1.8500000000E+09 * exp( 0.4 * lgt - 22422000 / rt );
	k[reth180b] = 1.9156952595E+06 * exp( 0.943736 * lgt - 45748021.07 / rt );
	k[reth181] = 3.9000000000E-10 * exp( 5.8 * lgt - 9205000 / rt );
	k[reth182] = 2.4500000000E-02 * exp( 3.15 * lgt - 23962000 / rt );
	k[reth183] = 2.4000000000E+16 * exp( -2.2 * lgt - 58702000 / rt );
	k[reth184] = 2.3200000000E+08 * exp( 0.4 * lgt - 62191000 / rt );
	k[reth185f] = 1.0000000000E+11 * exp( -176565000 / rt );
	k[reth185b] = 2.9665373901E+08 * exp( 0.332827 * lgt - 15981058.31 / rt );
	k[reth200f] = 1.0000000000E+09 * exp( 4602000 / rt );
	k[reth200b] = 5.1274620564E+15 * exp( -0.470019 * lgt - 122845559.2 / rt );
	k[reth201f] = 6.0000000000E+10 * exp( -102508000 / rt );
	k[reth201b] = 1.4635785336E-07 * exp( 2.42188 * lgt - 128881545.3 / rt );
	k[reth234] = 5.0000000000E+10;
	k[reth235] = 2.0000000000E+10;
	k[reth236f] = 1.0000000000E+11;
	k[reth236b] = 1.8964600615E+07 * exp( 1.02735 * lgt - 344375900.3 / rt );
	k[reth237] = 3.0000000000E+10;
	k[reth238] = 3.0000000000E+07;
	k[reth239] = 4.9000000000E+11 * exp( -0.5 * lgt );
	k[reth240] = 7.0000000000E+09;
	k[reth241f] = 3.0000000000E+09;
	k[reth241b] = 4.7509941111E+11 * exp( -0.178878 * lgt - 208419416.7 / rt );
	k[reth242f] = 1.1700000000E+43 * exp( -9.83 * lgt - 183075000 / rt );
	k[reth242b] = 3.3346351640E+33 * exp( -8.29873 * lgt - 168158562.1 / rt );
	k[reth243] = 1.8100000000E+43 * exp( -9.61 * lgt - 191912000 / rt );
	k[reth260f] = 4.0900000000E+02 * exp( 2.39 * lgt - 3616000 / rt );
	k[reth260b] = 4.4325383135E+06 * exp( 1.48044 * lgt - 127721493.7 / rt );
	k[reth261] = 5.0000000000E+10;
	k[reth262f] = 2.0000000000E+10;
	k[reth262b] = 5.0269250498E+14 * exp( -1.09908 * lgt - 212040585.7 / rt );
	k[reth263f] = 9.0400000000E+09 * exp( 1912000 / rt );
	k[reth263b] = 1.1780430570E+04 * exp( 0.809132 * lgt - 564376413.8 / rt );
	k[reth267] = 2.9500000000E+10 * exp( -4657000 / rt );
	k[reth271] = 1.0000000000E+10;
	kTroe0 = 2.4500000000E+12 * exp( -0.64 * lgt - 207945000 / rt );
	kTroeInf = 8.0000000000E+14 * exp( -0.52 * lgt - 212338000 / rt );
	fcTroe = 1 * exp( -0 / temp );
	k[rhai193f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM29] );
	kTroe0 = 3.7555644894E+09 * exp( -0.0527248 * lgt - 18054967.59 / rt );
	kTroeInf = 1.2263067720E+12 * exp( 0.0672752 * lgt - 22447967.59 / rt );
	k[rhai193b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM29] );
	k[rhai194f] = 1.0000000000E+10;
	k[rhai194b] = 8.3437759661E+10 * exp( -0.0125926 * lgt - 421399988.7 / rt );
	k[rhai195f] = 1.0000000000E+10;
	k[rhai195b] = 8.3622473155E+07 * exp( 0.317696 * lgt - 327394710.9 / rt );
	k[rhai348f] = 4.2200000000E+10 * exp( -259826000 / rt );
	k[rhai348b] = 3.9192829801E+07 * exp( 0.339519 * lgt - 8227844.371 / rt );
	k[rhai349f] = 4.9000000000E+09 * exp( 0.42 * lgt - 317147000 / rt );
	k[rhai349b] = 1.3229478695E+05 * exp( 0.670417 * lgt - 396974381.2 / rt );
	kTroe0 = 7.0000000000E+47 * exp( -9.31 * lgt - 417814000 / rt );
	kTroeInf = 8.0000000000E+12 * exp( 0.44 * lgt - 371414000 / rt );
	fcTroe = 0.2655 * exp( -temp / 180 ) + 0.7345 * exp( -temp / 1035 ) + 1 * exp( -5417 / temp );
	k[rhai350f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM35] );
	kTroe0 = 3.1654308814E+39 * exp( -8.46614 * lgt - 49435769.44 / rt );
	kTroeInf = 3.6176352930E+04 * exp( 1.28386 * lgt - 3035769.436 / rt );
	k[rhai350b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM35] );
	k[rhai352f] = 2.5000000000E+09;
	k[rhai352b] = 8.1261795561E+13 * exp( 0.182642 * lgt - 368480745.9 / rt );
	k[rhai353f] = 2.5000000000E+09;
	k[rhai353b] = 1.8193180871E+14 * exp( 0.100232 * lgt - 370800081 / rt );
	k[rhai354f] = 3.0000000000E+08;
	k[rhai354b] = 6.9410903482E+09 * exp( 0.196611 * lgt - 167100995.2 / rt );
	k[rhai355f] = 3.0000000000E+08;
	k[rhai355b] = 1.5539960848E+10 * exp( 0.114201 * lgt - 169420330.3 / rt );
	k[rhai356f] = 4.0000000000E+11 * exp( -175310000 / rt );
	k[rhai356b] = 5.9567053916E-03 * exp( 2.21623 * lgt - 280322939.2 / rt );
	k[rhai357f] = 3.1200000000E+55 * exp( -13.93 * lgt - 381162000 / rt );
	k[rhai357b] = 6.9851673752E+55 * exp( -14.0124 * lgt - 383481335.1 / rt );
	k[rnog8f] = 1.5000000000E+10;
	k[rnog8b] = 4.8515577882E+10 * exp( 0.172986 * lgt - 788045613 / rt );
	k[rnog9f] = 1.5000000000E+10;
	k[rnog9b] = 1.0436059884E+12 * exp( -0.40932 * lgt - 419143898 / rt );
	k[rnog51f] = 3.0000000000E+10;
	k[rnog51b] = 1.7669826979E+09 * exp( 0.195649 * lgt - 51987116.22 / rt );
	k[rnog62f] = 6.0000000000E+09;
	k[rnog62b] = 9.4919169392E+03 * exp( 1.60586 * lgt - 7026119.562 / rt );
	k[rnog67f] = 1.6000000000E+10;
	k[rnog67b] = 6.0389242625E+03 * exp( 1.39977 * lgt - 38585202.52 / rt );
	k[rnog79f] = 1.0000000000E+11;
	k[rnog79b] = 1.8286756321E+05 * exp( 1.64778 * lgt - 72258311.24 / rt );
	k[rnog94f] = 3.0000000000E+10;
	k[rnog94b] = 2.1487124305E+15 * exp( -0.825211 * lgt - 372657800.3 / rt );
	k[rnog97f] = 2.5010000000E+10;
	k[rnog97b] = 7.8591092956E+08 * exp( 0.489444 * lgt + 4323162.781 / rt );
	k[rnog142f] = 1.5000000000E+10 * exp( -2510640 / rt );
	k[rnog142b] = 6.7453611699E+09 * exp( -0.00559423 * lgt - 45106260.69 / rt );
	k[rnog144f] = 2.8000000000E+10;
	k[rnog144b] = 7.9385589251E+05 * exp( 0.260437 * lgt - 281989516.8 / rt );
	k[rnog145f] = 1.2000000000E+10;
	k[rnog145b] = 7.9258453041E+08 * exp( 0.51432 * lgt - 780192786.3 / rt );
	k[rnog146f] = 7.0000000000E+10;
	k[rnog146b] = 2.2524813868E+13 * exp( -0.586449 * lgt - 69401243.77 / rt );
	kTroe0 = 2.7000000000E+32 * exp( -6.3 * lgt - 12971700 / rt );
	kTroeInf = 2.0000000000E+10;
	fcTroe = 0.8493 * exp( -temp / 134 ) + 0.1507 * exp( -temp / 2383 ) + 1 * exp( -7265 / temp );
	k[rnog147f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM40] );
	kTroe0 = 1.9568400826E+42 * exp( -7.5866 * lgt - 408577469 / rt );
	kTroeInf = 1.4495111723E+20 * exp( -1.2866 * lgt - 395605769 / rt );
	k[rnog147b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM40] );
	k[rnog148f] = 3.0000000000E+10;
	k[rnog148b] = 1.3490722340E+10 * exp( -0.00559423 * lgt - 42595620.69 / rt );
	k[rnog149f] = 1.2000000000E+10 * exp( 2385110 / rt );
	k[rnog149b] = 6.9369801954E+16 * exp( -1.27256 * lgt - 309616557.3 / rt );
	k[rnog150f] = 1.6000000000E+10 * exp( 2385110 / rt );
	k[rnog150b] = 1.5915907230E+10 * exp( -0.307277 * lgt - 60070104.15 / rt );
	k[rnog151f] = 9.0000000000E+09;
	k[rnog151b] = 4.0472167020E+09 * exp( -0.00559423 * lgt - 42595620.69 / rt );
	k[rnog152f] = 7.0000000000E+09;
	k[rnog152b] = 3.1478352126E+09 * exp( -0.00559423 * lgt - 42595620.69 / rt );
	k[rnog153f] = 1.4000000000E+10;
	k[rnog153b] = 1.8215154484E+08 * exp( 0.457368 * lgt - 262818268.6 / rt );
	k[rnog154f] = 4.0000000000E+10 * exp( 2301420 / rt );
	k[rnog154b] = 1.3724062641E+08 * exp( 0.363879 * lgt - 78036865.42 / rt );
	k[rnog49f] = 1.1000000000E+11;
	k[rnog49b] = 4.3110693870E+10 * exp( 0.268172 * lgt - 94080437.51 / rt );
	k[rnog90f] = 5.0000000000E+10;
	k[rnog90b] = 1.1614017272E+13 * exp( -0.257094 * lgt - 650165751 / rt );
	k[rnog122f] = 5.8000000000E+10 * exp( -2410220 / rt );
	k[rnog122b] = 2.7207713449E+10 * exp( 0.181246 * lgt - 579645063.3 / rt );
	k[rnog123f] = 5.0000000000E+10;
	k[rnog123b] = 5.6190050992E+14 * exp( -0.800194 * lgt - 323271233.4 / rt );
	k[rnog124f] = 5.0000000000E+10;
	k[rnog124b] = 8.5101927116E+15 * exp( -1.1289 * lgt - 420571103.9 / rt );
	k[rnog178f] = 3.5000000000E+10 * exp( -1380850 / rt );
	k[rnog178b] = 5.6244744521E+10 * exp( 0.125609 * lgt - 315550484.7 / rt );
	k[rnog179f] = 2.6500000000E+09 * exp( -26780200 / rt );
	k[rnog179b] = 2.6413710670E+08 * exp( 0.0900563 * lgt - 158453538.8 / rt );
	k[rnog180f] = 7.3330000000E+10 * exp( -4686530 / rt );
	k[rnog180b] = 3.6192088004E+12 * exp( -0.348283 * lgt - 209290776.5 / rt );
	k[rnog181f] = 1.4000000000E+09 * exp( -45233400 / rt );
	k[rnog181b] = 1.9769877959E+05 * exp( 1.001 * lgt - 372613989.4 / rt );
	k[rnog182f] = 2.9000000000E+10 * exp( -96868900 / rt );
	k[rnog182b] = 2.5400577397E+05 * exp( 0.965452 * lgt - 241753193.4 / rt );
	k[rnog183f] = 4.4000000000E+11 * exp( -79001500 / rt );
	k[rnog183b] = 1.2548164902E+05 * exp( 1.43934 * lgt - 333451181.6 / rt );
	k[rnog184f] = 2.0000000000E+09 * exp( -88123500 / rt );
	k[rnog184b] = 2.8550511833E+06 * exp( 0.669888 * lgt - 191946343.3 / rt );
	kTroe0 = 6.2000000000E+11 * exp( -234745000 / rt );
	kTroeInf = 1.3000000000E+11 * exp( -249474000 / rt );
	fcTroe = 1 * exp( -0 / temp );
	k[rnog185f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM40] );
	kTroe0 = 1.2722932335E+03 * exp( 1.1222 * lgt - 64257184.86 / rt );
	kTroeInf = 2.6677116187E+02 * exp( 1.1222 * lgt - 78986184.86 / rt );
	k[rnog185b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM40] );
	k[rnog186f] = 2.1100000000E+09 * exp( 2008510 / rt );
	k[rnog186b] = 7.5425599797E+10 * exp( -0.250502 * lgt - 31611961.19 / rt );
	k[rnog187f] = 1.0600000000E+14 * exp( -1.41 * lgt );
	k[rnog187b] = 2.6359098753E+21 * exp( -2.11281 * lgt - 307931129.6 / rt );
	k[rnog188f] = 3.9000000000E+09 * exp( 1004260 / rt );
	k[rnog188b] = 1.0792457320E+07 * exp( 0.581618 * lgt - 188933014.9 / rt );
	k[rnog189f] = 1.3200000000E+11 * exp( -1506390 / rt );
	k[rnog189b] = 7.3770249488E+05 * exp( 1.01996 * lgt - 118512757.2 / rt );
	k[rnog190f] = 5.0000000000E+10;
	k[rnog190b] = 2.7331843723E+12 * exp( -0.270467 * lgt - 298076640.8 / rt );
	k[rnog191f] = 3.2000000000E+10 * exp( -1380850 / rt );
	k[rnog191b] = 5.8754729629E+10 * exp( 0.111557 * lgt - 103040936 / rt );
	k[rnog192f] = 2.0000000000E+10;
	k[rnog192b] = 1.9261874113E+14 * exp( -0.738009 * lgt - 76688281.91 / rt );
	k[rnog193f] = 2.0000000000E+06 * exp( 1.2 * lgt );
	k[rnog193b] = 3.7131770463E+07 * exp( 1.21455 * lgt - 166738167 / rt );
	k[rnog194f] = 4.6100000000E+02 * exp( 2 * lgt - 27198600 / rt );
	k[rnog194b] = 8.9664608808E+03 * exp( 1.70033 * lgt - 30955974.19 / rt );
	k[rnog195f] = 1.2800000000E+03 * exp( 1.5 * lgt - 418440 / rt );
	k[rnog195b] = 1.4130596080E+02 * exp( 1.66787 * lgt - 225564173.1 / rt );
	k[rnog196f] = 1.5000000000E+10;
	k[rnog196b] = 1.3176622007E+12 * exp( -0.144858 * lgt - 612246275.5 / rt );
	k[rnog197f] = 2.0000000000E+10 * exp( -57954000 / rt );
	k[rnog197b] = 1.9049155742E+13 * exp( -0.641003 * lgt - 69564200.92 / rt );
	k[rnog198f] = 2.1600000000E+10 * exp( -0.23 * lgt );
	k[rnog198b] = 3.8444536166E+10 * exp( -0.0265746 * lgt - 407642029 / rt );
	k[rnog199f] = 4.1600000000E+11 * exp( -0.45 * lgt );
	k[rnog199b] = 2.5962509240E+18 * exp( -1.68592 * lgt - 153192347.4 / rt );
	k[rnog200f] = 7.0000000000E+09;
	k[rnog200b] = 3.0905311820E+08 * exp( 0.251003 * lgt - 42668690.42 / rt );
	k[rnog201f] = 4.6000000000E+10;
	k[rnog201b] = 1.9559667416E+13 * exp( -0.487006 * lgt - 119356972.3 / rt );
	k[rnog202f] = 4.0000000000E+10 * exp( -15273100 / rt );
	k[rnog202b] = 2.9276621945E+09 * exp( 0.284744 * lgt - 66129482.11 / rt );
	k[rnog203f] = 9.0000000000E+04 * exp( 1.5 * lgt + 1924830 / rt );
	k[rnog203b] = 6.6607984180E+04 * exp( 1.68774 * lgt - 114009633.1 / rt );
	k[rnog204f] = 3.3000000000E+08;
	k[rnog204b] = 9.8038455525E+04 * exp( 0.282734 * lgt - 31159449.71 / rt );
	k[rnog205f] = 1.3000000000E+11 * exp( -0.11 * lgt - 20838300 / rt );
	k[rnog205b] = 3.8621209752E+07 * exp( 0.172734 * lgt - 51997749.71 / rt );
	k[rnog206f] = 5.0000000000E+09;
	k[rnog206b] = 2.0868560283E+09 * exp( 0.268765 * lgt - 232539200.4 / rt );
	k[rnog207f] = 2.5000000000E+10;
	k[rnog207b] = 1.0321778678E+09 * exp( 0.599881 * lgt - 456096946.5 / rt );
	k[rnog208f] = 7.0000000000E+10;
	k[rnog208b] = 1.6237968686E+09 * exp( 0.396456 * lgt - 48454917.49 / rt );
	k[rnog209f] = 5.0000000000E+10;
	k[rnog209b] = 3.4222396741E+09 * exp( 0.633623 * lgt - 464284638.2 / rt );
	k[rnog210f] = 2.0000000000E+10;
	k[rnog210b] = 1.3841820746E+10 * exp( 0.536617 * lgt - 529362719.2 / rt );
	k[rnog211f] = 2.5000000000E+10;
	k[rnog211b] = 5.5351855085E+11 * exp( 0.35445 * lgt - 471230667.8 / rt );
	k[rnog212f] = 8.9500000000E+13 * exp( -1.32 * lgt - 3096460 / rt );
	k[rnog212b] = 2.1914151927E+18 * exp( -1.47039 * lgt - 206645597.9 / rt );
	k[rnog213f] = 2.5000000000E+10;
	k[rnog213b] = 1.4189607700E+08 * exp( 0.467542 * lgt - 221388358.9 / rt );
	k[rnog214f] = 4.5000000000E+08 * exp( 0.72 * lgt - 2761710 / rt );
	k[rnog214b] = 4.2341747438E+06 * exp( 1.22128 * lgt - 232337760.6 / rt );
	k[rnog215f] = 1.3000000000E+04 * exp( 1.9 * lgt + 3975180 / rt );
	k[rnog215b] = 1.2368653509E+03 * exp( 2.30428 * lgt - 290678951.6 / rt );
	k[rnog216f] = 1.0000000000E+10 * exp( -54397200 / rt );
	k[rnog216b] = 5.7377065124E+08 * exp( 0.136425 * lgt - 52227812.81 / rt );
	k[rnog217f] = 7.7000000000E+10;
	k[rnog217b] = 2.7111042230E+11 * exp( 0.0458286 * lgt - 325738498.6 / rt );
	k[rnog218f] = 4.0000000000E+10;
	k[rnog218b] = 1.0292641680E+17 * exp( -1.43017 * lgt - 133658696.8 / rt );
	k[rnog219f] = 8.0000000000E+09 * exp( -31215700 / rt );
	k[rnog219b] = 2.6355653306E+12 * exp( -0.742772 * lgt - 55487486.99 / rt );
	k[rnog220f] = 6.1400000000E+09 * exp( 1841140 / rt );
	k[rnog220b] = 3.1907062411E+13 * exp( -0.991831 * lgt - 58886649.06 / rt );
	k[rnog221f] = 2.1000000000E+10 * exp( -19708500 / rt );
	k[rnog221b] = 6.9956150179E+13 * exp( -0.839778 * lgt - 109058368 / rt );
	k[rnog222f] = 2.3500000000E+10;
	k[rnog222b] = 1.5870420808E+06 * exp( 1.12772 * lgt - 396684048.3 / rt );
	k[rnog223f] = 5.4000000000E+10;
	k[rnog223b] = 6.6713759554E+04 * exp( 1.39818 * lgt - 98607407.55 / rt );
	k[rnog224f] = 2.5000000000E+09;
	k[rnog224b] = 1.2148665539E+03 * exp( 0.810569 * lgt + 28253448.44 / rt );
	k[rnog225f] = 2.0000000000E+10;
	k[rnog225b] = 2.1705234378E+06 * exp( 1.25333 * lgt - 710853683.1 / rt );
	k[rnog226f] = 2.0000000000E+09 * exp( -83688100 / rt );
	k[rnog226b] = 1.5016004241E+09 * exp( 0.283477 * lgt - 517280772.3 / rt );
	k[rnog227f] = 8.8000000000E+13 * exp( -0.5 * lgt - 200851000 / rt );
	k[rnog227b] = 8.6644157880E+05 * exp( 0.658852 * lgt + 32006694.95 / rt );
	k[rnog228f] = 2.8500000000E+14 * exp( -1.52 * lgt - 3096460 / rt );
	k[rnog228b] = 2.1974542854E+15 * exp( -1.35774 * lgt - 254896214.9 / rt );
	k[rnog229f] = 5.7000000000E+15 * exp( -2 * lgt - 3347520 / rt );
	k[rnog229b] = 6.8996904036E+16 * exp( -1.68097 * lgt - 619436488.3 / rt );
	k[rnog230f] = 1.0400000000E+26 * exp( -3.3 * lgt - 529746000 / rt );
	k[rnog230b] = 1.3550921711E+20 * exp( -2.81111 * lgt - 7270943.549 / rt );
	k[rnog231f] = 1.1070000000E+01 * exp( 2.64 * lgt - 20838300 / rt );
	k[rnog231b] = 5.1580077839E+03 * exp( 2.01587 * lgt - 56959437.11 / rt );
	k[rnog232f] = 2.7670000000E+00 * exp( 2.64 * lgt - 20838300 / rt );
	k[rnog232b] = 1.5928144501E-03 * exp( 3.41405 * lgt - 155566844.7 / rt );
	k[rnog233f] = 2.1340000000E+06 * exp( 1.58 * lgt - 111305000 / rt );
	k[rnog233b] = 3.8642207202E+02 * exp( 2.38604 * lgt - 13767440.33 / rt );
	k[rnog234f] = 1.1000000000E+03 * exp( 2.03 * lgt - 55945500 / rt );
	k[rnog234b] = 4.5604003860E+06 * exp( 1.33652 * lgt - 26962051.4 / rt );
	k[rnog235] = 4.4000000000E+00 * exp( 2.26 * lgt - 26780200 / rt );
	k[rnog236f] = 1.6000000000E-01 * exp( 2.56 * lgt - 37659600 / rt );
	k[rnog236b] = 2.0861278488E-03 * exp( 3.08305 * lgt - 129719454.2 / rt );
	k[rnog237f] = 1.4000000000E+20 * exp( -3.4 * lgt - 7950370 / rt );
	k[rnog237b] = 6.8402340894E+23 * exp( -3.51983 * lgt - 114275952.6 / rt );
	k[rnog238f] = 6.0000000000E+10 * exp( -1673760 / rt );
	k[rnog238b] = 1.9891112510E+07 * exp( 0.921332 * lgt - 329998193.7 / rt );
	k[rnog239f] = 6.3000000000E+10 * exp( -192566000 / rt );
	k[rnog239b] = 5.2402458901E+10 * exp( -0.0802483 * lgt + 1780628.844 / rt );
	k[rnog240f] = 1.4570000000E+05 * exp( 1.1 * lgt - 85361800 / rt );
	k[rnog240b] = 1.5822299564E+08 * exp( 0.448146 * lgt - 74445476.64 / rt );
	kTroe0 = 1.3000000000E+19 * exp( -3.16 * lgt - 3096460 / rt );
	kTroeInf = 3.1000000000E+09 * exp( 0.15 * lgt );
	fcTroe = 0.333 * exp( -temp / 235 ) + 0.667 * exp( -temp / 2117 ) + 1 * exp( -4536 / temp );
	k[rnog241f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM40] );
	kTroe0 = 6.3192712531E+28 * exp( -4.78038 * lgt - 138990791.4 / rt );
	kTroeInf = 1.5069031450E+19 * exp( -1.47038 * lgt - 135894331.4 / rt );
	k[rnog241b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM40] );
	k[rnog242f] = 1.0000000000E+10 * exp( -309646000 / rt );
	k[rnog242b] = 7.7466563433E+11 * exp( -0.562168 * lgt - 206461086.2 / rt );
	k[rnog243f] = 1.0000000000E+08 * exp( -271986000 / rt );
	k[rnog243b] = 3.4835996596E+09 * exp( -0.567762 * lgt - 211396706.9 / rt );
	k[rnog244f] = 1.9000000000E+10;
	k[rnog244b] = 2.5396777936E+10 * exp( 0.045361 * lgt - 119823005.9 / rt );
	k[rnog245f] = 2.9000000000E+10;
	k[rnog245b] = 1.3648298329E+11 * exp( 0.0911896 * lgt - 445561504.5 / rt );
	k[rnog246f] = 5.0000000000E+10;
	k[rnog246b] = 8.7255730632E+13 * exp( -0.526245 * lgt - 303253311.4 / rt );
	k[rnog247f] = 2.0000000000E+10;
	k[rnog247b] = 1.6262537951E+16 * exp( -1.15038 * lgt - 339374448.5 / rt );
	k[rnog248f] = 3.0000000000E+10;
	k[rnog248b] = 1.1902808466E+12 * exp( -0.222945 * lgt - 170740226.9 / rt );
	k[rnog249f] = 3.1000000000E+14 * exp( -1.38 * lgt - 5314190 / rt );
	k[rnog249b] = 2.4074880728E+21 * exp( -2.7717 * lgt - 388492664.2 / rt );
	k[rnog250f] = 2.9000000000E+11 * exp( -0.69 * lgt - 3180150 / rt );
	k[rnog250b] = 3.9984637446E+13 * exp( -1.04874 * lgt - 307637265.2 / rt );
	k[rnog251f] = 3.8000000000E+10 * exp( -0.36 * lgt - 2426950 / rt );
	k[rnog251b] = 1.1350142698E+18 * exp( -1.98486 * lgt - 97106826.14 / rt );
	k[rnog252f] = 3.1000000000E+14 * exp( -1.38 * lgt - 5314190 / rt );
	k[rnog252b] = 1.0826251042E+21 * exp( -2.77729 * lgt - 431088284.9 / rt );
	k[rnog253f] = 2.9000000000E+11 * exp( -0.69 * lgt - 3180150 / rt );
	k[rnog253b] = 1.7980721388E+13 * exp( -1.05434 * lgt - 350232885.9 / rt );
	k[rnog254f] = 3.8000000000E+10 * exp( -0.36 * lgt - 2426950 / rt );
	k[rnog254b] = 5.1040541218E+17 * exp( -1.99046 * lgt - 139702446.8 / rt );
	k[rnog255f] = 9.6000000000E+10 * exp( -120511000 / rt );
	k[rnog255b] = 1.8704236432E+11 * exp( 0.125107 * lgt - 463240573.1 / rt );
	k[rnog256f] = 1.0000000000E+09 * exp( -91010800 / rt );
	k[rnog256b] = 4.0863068469E+09 * exp( -0.248609 * lgt - 41862686.28 / rt );
	k[rnog257f] = 2.2000000000E+10;
	k[rnog257b] = 2.4852949814E+02 * exp( 1.59771 * lgt - 600164165.4 / rt );
	k[rnog258f] = 2.0000000000E+09;
	k[rnog258b] = 3.0646463697E+06 * exp( 0.937388 * lgt - 482731088.5 / rt );
	k[rnog259f] = 1.2000000000E+10;
	k[rnog259b] = 9.7626862770E+00 * exp( 1.48749 * lgt - 166519234.3 / rt );
	k[rnog260f] = 1.2000000000E+10;
	k[rnog260b] = 4.8341235444E+03 * exp( 1.04915 * lgt - 239450142 / rt );
	k[rnog261f] = 1.0000000000E+11;
	k[rnog261b] = 3.6185714969E+04 * exp( 1.77002 * lgt - 287839361.6 / rt );
	k[rnog262f] = 9.8000000000E+04 * exp( 1.41 * lgt - 35567400 / rt );
	k[rnog262b] = 5.5134956480E+03 * exp( 1.93443 * lgt - 201414117.4 / rt );
	k[rnog263f] = 1.5000000000E+05 * exp( 1.57 * lgt - 184114000 / rt );
	k[rnog263b] = 1.4764153289E+01 * exp( 2.639 * lgt - 316809467.6 / rt );
	k[rnog264f] = 2.2000000000E+03 * exp( 2.11 * lgt - 47702200 / rt );
	k[rnog264b] = 1.8199093965E+01 * exp( 2.51882 * lgt - 5101978.111 / rt );
	k[rnog265f] = 2.2500000000E+04 * exp( 1.7 * lgt - 15900700 / rt );
	k[rnog265b] = 5.2083021419E-03 * exp( 3.256 * lgt - 29239195.24 / rt );
	k[rnog266f] = 1.0500000000E+02 * exp( 2.5 * lgt - 55652600 / rt );
	k[rnog266b] = 1.4399329478E+00 * exp( 2.94256 * lgt - 21240069.8 / rt );
	k[rnog267f] = 4.6500000000E+09 * exp( -28663200 / rt );
	k[rnog267b] = 6.4480549612E+08 * exp( 0.345558 * lgt - 59328750.79 / rt );
	k[rnog268f] = 1.5500000000E+09 * exp( -28663200 / rt );
	k[rnog268b] = 1.9751385652E+09 * exp( 0.273424 * lgt - 151841226.9 / rt );
	k[rnog269f] = 1.1800000000E+13 * exp( -354503000 / rt );
	k[rnog269b] = 8.6775822032E+02 * exp( 1.48986 * lgt + 14427311.12 / rt );
	k[rnog270f] = 2.1000000000E+12 * exp( -0.69 * lgt - 11925600 / rt );
	k[rnog270b] = 5.4601430394E+11 * exp( -0.456838 * lgt - 300424198.1 / rt );
	k[rnog271f] = 2.7000000000E+08 * exp( 0.18 * lgt - 8870940 / rt );
	k[rnog271b] = 1.2463534027E+03 * exp( 1.44612 * lgt - 218648179.1 / rt );
	k[rnog272f] = 1.7000000000E+11 * exp( -0.75 * lgt - 12092900 / rt );
	k[rnog272b] = 1.0231688247E+04 * exp( 1.03916 * lgt - 313929993.3 / rt );
	k[rnog273] = 2.0000000000E+04 * exp( 2 * lgt - 8368810 / rt );
	k[rnog274f] = 2.3500000000E+10;
	k[rnog274b] = 5.7721314694E+11 * exp( 0.017323 * lgt - 209533808.1 / rt );
	k[rnog275f] = 6.1000000000E+11 * exp( -0.31 * lgt - 1213480 / rt );
	k[rnog275b] = 1.2302482746E+14 * exp( -0.906892 * lgt - 156669612.8 / rt );
	k[rnog276f] = 3.7000000000E+09 * exp( 0.15 * lgt + 376596 / rt );
	k[rnog276b] = 3.5186787743E+10 * exp( 0.0238295 * lgt - 481879142.6 / rt );
	k[rnog277f] = 5.4000000000E+02 * exp( 2.4 * lgt - 41488400 / rt );
	k[rnog277b] = 2.8838511300E+00 * exp( 2.76138 * lgt - 18472241.08 / rt );
	k[rnog278f] = 5.0000000000E+04 * exp( 1.6 * lgt - 3996110 / rt );
	k[rnog278b] = 2.7000505166E+03 * exp( 1.86437 * lgt - 46058032.07 / rt );
	k[rnog279f] = 9.4000000000E+03 * exp( 1.94 * lgt - 27031200 / rt );
	k[rnog279b] = 3.0281755443E+01 * exp( 2.26763 * lgt + 4172650.606 / rt );
	k[rnog280f] = 4.0950565564E+12 * exp( -0.2 * lgt );
	k[rnog280b] = 1.9054048295E+12 * exp( -0.0903729 * lgt - 74098001.86 / rt );


	w[r49] = k[r49] * c[sHO2] * c[sOH];
	w[r61f] = k[r61f] * c[sC2H4] * c[sH];
	w[r61b] = k[r61b] * c[sC2H3] * c[sH2];
	w[r62f] = k[r62f] * c[sC2H4] * c[sOH];
	w[r62b] = k[r62b] * c[sC2H3] * c[sH2O];
	w[r112] = k[r112] * c[sC2H3] * c[sO2];
	w[r113f] = k[r113f] * c[sC2H6];
	w[r113b] = k[r113b] * c[sH] * c[sC2H5];
	w[r533] = k[r533] * c[sCH3O2] * c[sCH3];
	w[r535] = k[r535] * c[sCH3O2] * c[sHO2];
	w[ri26f] = k[ri26f] * c[sH] * c[sO2];
	w[ri26b] = k[ri26b] * c[sHO2];
	w[ri5f] = k[ri5f] * c[sC2H6] * c[sCH3];
	w[ri5b] = k[ri5b] * c[sC2H5] * c[sCH4];
	w[ri7f] = k[ri7f] * c[sCO] * c[sOH];
	w[ri7b] = k[ri7b] * c[sCO2] * c[sH];
	w[ri8f] = k[ri8f] * c[sH] * c[sO2];
	w[ri8b] = k[ri8b] * c[sO] * c[sOH];
	w[ri10f] = k[ri10f] * c[sO] * c[sH2O];
	w[ri10b] = k[ri10b] * c[sOH] * c[sOH];
	w[ri13P1f] = k[ri13P1f] * c[sH2O2] * c[sOH];
	w[ri13P1b] = k[ri13P1b] * c[sH2O] * c[sHO2];
	w[ri17f] = k[ri17f] * c[sC2H6] * c[sH];
	w[ri17b] = k[ri17b] * c[sC2H5] * c[sH2];
	w[ri22f] = k[ri22f] * c[sCH3] * c[sHO2];
	w[ri22b] = k[ri22b] * c[sCH3O] * c[sOH];
	w[ri23f] = k[ri23f] * c[sCO] * c[sHO2];
	w[ri23b] = k[ri23b] * c[sCO2] * c[sOH];
	w[ri27f] = k[ri27f] * c[sCO] * c[sO];
	w[ri27b] = k[ri27b] * c[sCO2];
	w[ri29f] = k[ri29f] * c[sHCO] * c[sH];
	w[ri29b] = k[ri29b] * c[sCO] * c[sH2];
	w[ri31f] = k[ri31f] * c[sCH2O] * M[mM0];
	w[ri31b] = k[ri31b] * c[sHCO] * c[sH] * M[mM0];
	w[ri33f] = k[ri33f] * c[sCH2O] * c[sH];
	w[ri33b] = k[ri33b] * c[sHCO] * c[sH2];
	w[ri36f] = k[ri36f] * c[sCH3] * c[sO];
	w[ri36b] = k[ri36b] * c[sCH2O] * c[sH];
	w[ri37f] = k[ri37f] * c[sCH3] * c[sO2];
	w[ri37b] = k[ri37b] * c[sCH3O] * c[sO];
	w[ri40] = k[ri40] * c[sCH3O];
	w[ri42f] = k[ri42f] * c[sHO2] * c[sO];
	w[ri42b] = k[ri42b] * c[sOH] * c[sO2];
	w[ri43f] = k[ri43f] * c[sHCO] * c[sHO2];
	w[ri43b] = k[ri43b] * c[sCH2O] * c[sO2];
	w[ri44] = k[ri44] * c[sCH3O] * c[sO2];
	w[ri46] = k[ri46] * c[sHCO] * c[sO2];
	w[ri47] = k[ri47] * c[sHO2] * c[sH];
	w[ri48] = k[ri48] * c[sHO2] * c[sH];
	w[ri52] = k[ri52] * c[sH2O2] * c[sH];
	w[ri54] = k[ri54] * c[sCH2O] * c[sHO2];
	w[ri60f] = k[ri60f] * c[sC2H2] * c[sH];
	w[ri60b] = k[ri60b] * c[sC2H3];
	w[ri69] = k[ri69] * c[sC2H2] * c[sO];
	w[ri72] = k[ri72] * c[sCH2] * c[sO2];
	w[ri83] = k[ri83] * c[sCH2OH] * c[sO2];
	w[ri84f] = k[ri84f] * c[sCH2OH];
	w[ri84b] = k[ri84b] * c[sH] * c[sCH2O];
	w[ri87] = k[ri87] * c[sC2H2] * c[sO];
	w[ri94] = k[ri94] * c[sCH2CO] * c[sOH];
	w[ri111f] = k[ri111f] * c[sC2H5] * c[sH];
	w[ri111b] = k[ri111b] * c[sCH3] * c[sCH3];
	w[ri136] = k[ri136] * c[sCH3CO];
	w[ri146] = k[ri146] * c[sAXC3H5];
	w[ri159] = k[ri159] * c[sC3H6] * c[sOH];
	w[ri177] = k[ri177] * c[sC3H6] * c[sH];
	w[ri521f] = k[ri521f] * c[sCH3O2] * M[mM0];
	w[ri521b] = k[ri521b] * c[sCH3] * c[sO2] * M[mM0];
	w[ri522] = k[ri522] * c[sCH3O2H];
	w[ri568] = k[ri568] * c[sC3H5O];
	w[ri570] = k[ri570] * c[sAXC3H5] * c[sHO2];
	w[ri1455] = k[ri1455] * c[sHCCO] * c[sO2];
	w[ri1457] = k[ri1457] * c[sCH2CO] * c[sOH];
	w[r16f] = k[r16f] * c[sCH3] * c[sOH];
	w[r16b] = k[r16b] * c[sCH3OH];
	w[r63] = k[r63] * c[sC2H4] * c[sO];
	w[r65] = k[r65] * c[sC2H2] * c[sO2];
	w[r85] = k[r85] * c[sC2H3] * c[sO2];
	w[r91] = k[r91] * c[sCH2CO] * c[sOH];
	w[r104f] = k[r104f] * c[sCH3OH] * c[sCH2O];
	w[r104b] = k[r104b] * c[sCH3O] * c[sCH3O];
	w[r135] = k[r135] * c[sC2H4] * c[sCH3];
	w[r137] = k[r137] * c[sCH3CHO];
	w[r142] = k[r142] * c[sCH3CHO] * c[sHO2];
	w[r169] = k[r169] * c[sNXC3H7];
	w[r170f] = k[r170f] * c[sNXC3H7];
	w[r170b] = k[r170b] * c[sH] * c[sC3H6];
	w[r171] = k[r171] * c[sNXC3H7] * c[sO2];
	w[r174] = k[r174] * c[sC3H6] * c[sOH];
	w[r218] = k[r218] * c[sC5H9];
	w[r255] = k[r255] * c[sPXC4H9];
	w[r330] = k[r330] * c[s1XC5H10];
	w[r331] = k[r331] * c[s1XC5H10] * c[sH];
	w[r335] = k[r335] * c[s1XC5H10] * c[sOH];
	w[r464] = k[r464] * c[s1XC5H11];
	w[r465f] = k[r465f] * c[s1XC5H11];
	w[r465b] = k[r465b] * c[sH] * c[s1XC5H10];
	w[r536] = k[r536] * c[sCH3O2] * c[sHO2];
	w[r538] = k[r538] * c[sCH3O2] * c[sCH3O2];
	w[r1014f] = k[r1014f] * c[sCH3COCH2O2];
	w[r1014b] = k[r1014b] * c[sCH3COCH2] * c[sO2];
	w[r1017] = k[r1017] * c[sHO2] * c[sCH3COCH2O2];
	w[r1018] = k[r1018] * c[sCH3COCH2O2H];
	w[r1019] = k[r1019] * c[sCH3COCH2O];
	w[r1190f] = k[r1190f] * c[sHO2] * M[mM0];
	w[r1190b] = k[r1190b] * c[sO] * c[sOH] * M[mM0];
	w[r1330] = k[r1330] * c[sC3H6] * c[sCH3];
	w[r1335] = k[r1335] * c[sAXC3H5] * c[sHO2];
	w[r1400f] = k[r1400f] * c[sCH2CHO];
	w[r1400b] = k[r1400b] * c[sCH2CO] * c[sH];
	w[r1401f] = k[r1401f] * c[sCH2CHO] * c[sO2];
	w[r1401b] = k[r1401b] * c[sCH2CO] * c[sHO2];
	w[r1452] = k[r1452] * c[sCH2CO] * c[sOH];
	w[r2078] = k[r2078] * c[sNXC7H16];
	w[r2079] = k[r2079] * c[sNXC7H16];
	w[r2080] = k[r2080] * c[sNXC7H16] * c[sH];
	w[r2081] = k[r2081] * c[sNXC7H16] * c[sH];
	w[r2082] = k[r2082] * c[sNXC7H16] * c[sH];
	w[r2083] = k[r2083] * c[sNXC7H16] * c[sH];
	w[r2088] = k[r2088] * c[sNXC7H16] * c[sOH];
	w[r2089] = k[r2089] * c[sNXC7H16] * c[sOH];
	w[r2090] = k[r2090] * c[sNXC7H16] * c[sOH];
	w[r2091] = k[r2091] * c[sNXC7H16] * c[sOH];
	w[r2092] = k[r2092] * c[sNXC7H16] * c[sHO2];
	w[r2093] = k[r2093] * c[sNXC7H16] * c[sHO2];
	w[r2094] = k[r2094] * c[sNXC7H16] * c[sHO2];
	w[r2095] = k[r2095] * c[sNXC7H16] * c[sHO2];
	w[r2100f] = k[r2100f] * c[sNXC7H16] * c[sO2];
	w[r2100b] = k[r2100b] * c[s1XC7H15] * c[sHO2];
	w[r2101f] = k[r2101f] * c[sNXC7H16] * c[sO2];
	w[r2101b] = k[r2101b] * c[s2XC7H15] * c[sHO2];
	w[r2102f] = k[r2102f] * c[sNXC7H16] * c[sO2];
	w[r2102b] = k[r2102b] * c[s3XC7H15] * c[sHO2];
	w[r2103f] = k[r2103f] * c[sNXC7H16] * c[sO2];
	w[r2103b] = k[r2103b] * c[s4XC7H15] * c[sHO2];
	w[r2116] = k[r2116] * c[sNXC7H16] * c[sCH3O2];
	w[r2117] = k[r2117] * c[sNXC7H16] * c[sCH3O2];
	w[r2118] = k[r2118] * c[sNXC7H16] * c[sCH3O2];
	w[r2119] = k[r2119] * c[sNXC7H16] * c[sCH3O2];
	w[r2142] = k[r2142] * c[s1XC7H15];
	w[r2144] = k[r2144] * c[s2XC7H15];
	w[r2147] = k[r2147] * c[s3XC7H15];
	w[r2151] = k[r2151] * c[s4XC7H15];
	w[r2159f] = k[r2159f] * c[s1XC7H15];
	w[r2159b] = k[r2159b] * c[s3XC7H15];
	w[r2160f] = k[r2160f] * c[s1XC7H15];
	w[r2160b] = k[r2160b] * c[s4XC7H15];
	w[r2161f] = k[r2161f] * c[s2XC7H15];
	w[r2161b] = k[r2161b] * c[s3XC7H15];
	w[r2162f] = k[r2162f] * c[s1XC7H15];
	w[r2162b] = k[r2162b] * c[s2XC7H15];
	w[r2185] = k[r2185] * c[s1XC7H14];
	w[r2188f] = k[r2188f] * c[s2XC7H15O2];
	w[r2188b] = k[r2188b] * c[s1XC7H15] * c[sO2];
	w[r2189f] = k[r2189f] * c[s2XC7H15O2];
	w[r2189b] = k[r2189b] * c[s2XC7H15] * c[sO2];
	w[r2190f] = k[r2190f] * c[s2XC7H15O2];
	w[r2190b] = k[r2190b] * c[s3XC7H15] * c[sO2];
	w[r2191f] = k[r2191f] * c[s2XC7H15O2];
	w[r2191b] = k[r2191b] * c[s4XC7H15] * c[sO2];
	w[r2212P1f] = k[r2212P1f] * c[s2XC7H15O2];
	w[r2212P1b] = k[r2212P1b] * c[s2X3XC7H14O2H];
	w[r2212P2f] = k[r2212P2f] * c[s2XC7H15O2];
	w[r2212P2b] = k[r2212P2b] * c[s2X3XC7H14O2H];
	w[r2214f] = k[r2214f] * c[s2XC7H15O2];
	w[r2214b] = k[r2214b] * c[s2X4XC7H14O2H];
	w[r2258f] = k[r2258f] * c[s2X3XC7H14O2H];
	w[r2258b] = k[r2258b] * c[s1XC7H14] * c[sHO2];
	w[r2281] = k[r2281] * c[s2X4XC7H14O2H];
	w[r2290f] = k[r2290f] * c[s2X4XC7H14OOHO2];
	w[r2290b] = k[r2290b] * c[s2X3XC7H14O2H] * c[sO2];
	w[r2291f] = k[r2291f] * c[s2X4XC7H14OOHO2];
	w[r2291b] = k[r2291b] * c[s2X4XC7H14O2H] * c[sO2];
	w[r2309] = k[r2309] * c[s2X4XC7H14OOHO2];
	w[r2327] = k[r2327] * c[sNC7KET24XC7H14O3];
	w[rnog52f] = k[rnog52f] * c[sH] * c[sCH3];
	w[rnog52b] = k[rnog52b] * c[sCH4];
	w[ri2f] = k[ri2f] * c[sCH4] * c[sH];
	w[ri2b] = k[ri2b] * c[sCH3] * c[sH2];
	w[ri3f] = k[ri3f] * c[sCH4] * c[sOH];
	w[ri3b] = k[ri3b] * c[sCH3] * c[sH2O];
	w[ri4f] = k[ri4f] * c[sCH4] * c[sO];
	w[ri4b] = k[ri4b] * c[sCH3] * c[sOH];
	w[ri6f] = k[ri6f] * c[sHCO] * c[sOH];
	w[ri6b] = k[ri6b] * c[sCO] * c[sH2O];
	w[ri9f] = k[ri9f] * c[sO] * c[sH2];
	w[ri9b] = k[ri9b] * c[sH] * c[sOH];
	w[ri11f] = k[ri11f] * c[sOH] * c[sH2];
	w[ri11b] = k[ri11b] * c[sH] * c[sH2O];
	w[ri12f] = k[ri12f] * c[sHCO] * M[mM22];
	w[ri12b] = k[ri12b] * c[sH] * c[sCO] * M[mM22];
	w[ri14f] = k[ri14f] * c[sC2H4] * c[sO];
	w[ri14b] = k[ri14b] * c[sCH3] * c[sHCO];
	w[ri15f] = k[ri15f] * c[sH] * c[sC2H4];
	w[ri15b] = k[ri15b] * c[sC2H5];
	w[ri19f] = k[ri19f] * c[sC2H5] * c[sO2];
	w[ri19b] = k[ri19b] * c[sC2H4] * c[sHO2];
	w[ri20f] = k[ri20f] * c[sC2H6] * c[sOH];
	w[ri20b] = k[ri20b] * c[sC2H5] * c[sH2O];
	w[ri21f] = k[ri21f] * c[sC2H6] * c[sO];
	w[ri21b] = k[ri21b] * c[sC2H5] * c[sOH];
	w[ri24f] = k[ri24f] * c[sCH3] * c[sCH3];
	w[ri24b] = k[ri24b] * c[sC2H6];
	w[ri25f] = k[ri25f] * c[sH2O] * M[mM22];
	w[ri25b] = k[ri25b] * c[sH] * c[sOH] * M[mM22];
	w[ri28f] = k[ri28f] * c[sCO] * c[sO2];
	w[ri28b] = k[ri28b] * c[sCO2] * c[sO];
	w[ri30f] = k[ri30f] * c[sHCO] * c[sO];
	w[ri30b] = k[ri30b] * c[sCO] * c[sOH];
	w[ri32f] = k[ri32f] * c[sCH2O] * c[sOH];
	w[ri32b] = k[ri32b] * c[sHCO] * c[sH2O];
	w[ri34f] = k[ri34f] * c[sCH2O] * c[sO];
	w[ri34b] = k[ri34b] * c[sHCO] * c[sOH];
	w[ri35f] = k[ri35f] * c[sCH3] * c[sOH];
	w[ri35b] = k[ri35b] * c[sCH2O] * c[sH2];
	w[ri38] = k[ri38] * c[sCH2O] * c[sCH3];
	w[ri45] = k[ri45] * c[sCH3] * c[sHO2];
	w[ri50P1f] = k[ri50P1f] * c[sH2O2] * c[sO2];
	w[ri50P1b] = k[ri50P1b] * c[sHO2] * c[sHO2];
	w[ri51f] = k[ri51f] * c[sOH] * c[sOH];
	w[ri51b] = k[ri51b] * c[sH2O2];
	w[ri53f] = k[ri53f] * c[sCH4] * c[sHO2];
	w[ri53b] = k[ri53b] * c[sCH3] * c[sH2O2];
	w[ri66] = k[ri66] * c[sCH2] * c[sO2];
	w[ri74f] = k[ri74f] * c[sCH2] * c[sH];
	w[ri74b] = k[ri74b] * c[sCH] * c[sH2];
	w[ri75f] = k[ri75f] * c[sCH2] * c[sOH];
	w[ri75b] = k[ri75b] * c[sCH] * c[sH2O];
	w[ri76] = k[ri76] * c[sCH2] * c[sO2];
	w[ri77f] = k[ri77f] * c[sCH] * c[sO2];
	w[ri77b] = k[ri77b] * c[sHCO] * c[sO];
	w[ri88] = k[ri88] * c[sC2H2] * c[sOH];
	w[ri89] = k[ri89] * c[sCH2CO] * c[sH];
	w[ri90] = k[ri90] * c[sCH2CO] * c[sO];
	w[ri91] = k[ri91] * c[sCH2] * c[sO2];
	w[ri92] = k[ri92] * c[sCH2CO];
	w[ri95] = k[ri95] * c[sCH2CO] * c[sH];
	w[ri96] = k[ri96] * c[sHCCO] * c[sOH];
	w[ri98] = k[ri98] * c[sHCCO] * c[sO];
	w[ri101] = k[ri101] * c[sCH2] * c[sO2];
	w[ri147] = k[ri147] * c[sC3H6];
	w[ri151f] = k[ri151f] * c[sC3H6];
	w[ri151b] = k[ri151b] * c[sAXC3H5] * c[sH];
	w[ri153] = k[ri153] * c[sC3H6] * c[sO];
	w[ri173] = k[ri173] * c[sC3H6] * c[sOH];
	w[ri174] = k[ri174] * c[sC3H6] * c[sO];
	w[ri176] = k[ri176] * c[sC3H6] * c[sO];
	w[ri179] = k[ri179] * c[sC3H6] * c[sH];
	w[ri180f] = k[ri180f] * c[sIXC3H7];
	w[ri180b] = k[ri180b] * c[sH] * c[sC3H6];
	w[ri181] = k[ri181] * c[sIXC3H7] * c[sH];
	w[ri182] = k[ri182] * c[sIXC3H7] * c[sO2];
	w[ri206] = k[ri206] * c[sIXC4H8];
	w[ri207f] = k[ri207f] * c[sIXC4H8];
	w[ri207b] = k[ri207b] * c[sIXC4H7] * c[sH];
	w[ri210] = k[ri210] * c[sIXC4H8] * c[sO];
	w[ri211] = k[ri211] * c[sIXC4H8] * c[sH];
	w[ri213] = k[ri213] * c[sIXC4H8] * c[sOH];
	w[ri214] = k[ri214] * c[sIXC4H8] * c[sO];
	w[ri215] = k[ri215] * c[sIXC4H8] * c[sCH3O2];
	w[ri217] = k[ri217] * c[sIXC4H8] * c[sHO2];
	w[ri259f] = k[ri259f] * c[sTXC4H9];
	w[ri259b] = k[ri259b] * c[sH] * c[sIXC4H8];
	w[ri262] = k[ri262] * c[sIXC4H9];
	w[ri263] = k[ri263] * c[sIXC4H9];
	w[ri308] = k[ri308] * c[sCH3COCH3];
	w[ri309] = k[ri309] * c[sCH3COCH3] * c[sOH];
	w[ri310] = k[ri310] * c[sCH3COCH3] * c[sH];
	w[ri314] = k[ri314] * c[sCH3COCH2];
	w[ri317] = k[ri317] * c[sC2H5CO];
	w[ri340f] = k[ri340f] * c[sH2O2] * c[sH];
	w[ri340b] = k[ri340b] * c[sH2] * c[sHO2];
	w[ri343f] = k[ri343f] * c[sCH3] * c[sH];
	w[ri343b] = k[ri343b] * c[sCH2] * c[sH2];
	w[ri344f] = k[ri344f] * c[sCH3] * c[sOH];
	w[ri344b] = k[ri344b] * c[sCH2] * c[sH2O];
	w[ri347] = k[ri347] * c[sCH3CO] * c[sH];
	w[ri493] = k[ri493] * c[sNEOXC5H11];
	w[ri50P2f] = k[ri50P2f] * c[sH2O2] * c[sO2];
	w[ri50P2b] = k[ri50P2b] * c[sHO2] * c[sHO2];
	w[ri525] = k[ri525] * c[sCH3O2] * c[sCH2O];
	w[ri13P2f] = k[ri13P2f] * c[sH2O2] * c[sOH];
	w[ri13P2b] = k[ri13P2b] * c[sH2O] * c[sHO2];
	w[ri540] = k[ri540] * c[sCH3O2] * c[sCH3O2];
	w[ri710] = k[ri710] * c[sIXC4H7] * c[sHO2];
	w[ri711f] = k[ri711f] * c[sIXC4H8] * c[sHO2];
	w[ri711b] = k[ri711b] * c[sIXC4H8O] * c[sOH];
	w[ri776f] = k[ri776f] * c[sIC4H8OOHXIO2XC4H9O4];
	w[ri776b] = k[ri776b] * c[sIC4H8OOHXIXC4H9O2] * c[sO2];
	w[ri778f] = k[ri778f] * c[sIC4H8OOHXTO2XC4H9O4];
	w[ri778b] = k[ri778b] * c[sIC4H8OOHXTXC4H9O2] * c[sO2];
	w[ri779] = k[ri779] * c[sIC4H8OOHXIO2XC4H9O4];
	w[ri780] = k[ri780] * c[sIC4H8OOHXTO2XC4H9O4];
	w[ri790] = k[ri790] * c[sTC4H8OOHXIXC4H9O2];
	w[ri791] = k[ri791] * c[sIC4H8OOHXTXC4H9O2];
	w[ri793] = k[ri793] * c[sIC4H8OOHXTXC4H9O2];
	w[ri794] = k[ri794] * c[sTC4H8OOHXIXC4H9O2];
	w[ri797] = k[ri797] * c[sIC4H8OOHXIXC4H9O2];
	w[ri832] = k[ri832] * c[sTXC4H9O];
	w[ri901f] = k[ri901f] * c[sIXC4H9O2];
	w[ri901b] = k[ri901b] * c[sIC4H8OOHXIXC4H9O2];
	w[ri902f] = k[ri902f] * c[sTXC4H9O2];
	w[ri902b] = k[ri902b] * c[sTC4H8OOHXIXC4H9O2];
	w[ri903f] = k[ri903f] * c[sIXC4H9O2];
	w[ri903b] = k[ri903b] * c[sIC4H8OOHXTXC4H9O2];
	w[ri1193] = k[ri1193] * c[sIXC4H7O];
	w[ri1347] = k[ri1347] * c[sC3H5XTXC3H5] * c[sO2];
	w[ri1393f] = k[ri1393f] * c[sTXC4H9O2];
	w[ri1393b] = k[ri1393b] * c[sTXC4H9] * c[sO2];
	w[ri1394] = k[ri1394] * c[sTXC4H9] * c[sHO2];
	w[ri1396f] = k[ri1396f] * c[sIXC4H9O2];
	w[ri1396b] = k[ri1396b] * c[sIXC4H9] * c[sO2];
	w[ri1412] = k[ri1412] * c[sIC4KETIIXC4H8O3];
	w[ri1458] = k[ri1458] * c[sCH3] * c[sO2];
	w[ri1460] = k[ri1460] * c[sC3H5XTXC3H5] * c[sO2];
	w[ri1461] = k[ri1461] * c[sTXC3H6CHO] * c[sHO2];
	w[ri1462] = k[ri1462] * c[sTXC3H6OCHO];
	w[ri1463] = k[ri1463] * c[sTXC3H6CHO];
	w[ri1466] = k[ri1466] * c[sIXC3H7CHO];
	w[ri1468] = k[ri1468] * c[sIXC4H7OOH];
	w[ri1508] = k[ri1508] * c[sTXC3H6CHO] * c[sHO2];
	w[ri1510] = k[ri1510] * c[sIXC3H5CHO] * c[sOH];
	w[ri1511] = k[ri1511] * c[sIXC3H5CHO] * c[sHO2];
	w[ri1516] = k[ri1516] * c[sIXC3H5CO];
	w[ri1517] = k[ri1517] * c[sIXC4H7O] * c[sO2];
	w[ri1536f] = k[ri1536f] * c[sIXC4H8] * c[sO2];
	w[ri1536b] = k[ri1536b] * c[sIXC4H7] * c[sHO2];
	w[ri1546] = k[ri1546] * c[sIXC4H7O];
	w[ri2720] = k[ri2720] * c[sYXC7H15];
	w[ri2722f] = k[ri2722f] * c[sYXC7H15];
	w[ri2722b] = k[ri2722b] * c[sYXC7H14] * c[sH];
	w[ri2727f] = k[ri2727f] * c[sYXC7H15] * c[sO2];
	w[ri2727b] = k[ri2727b] * c[sYXC7H14] * c[sHO2];
	w[ri2731f] = k[ri2731f] * c[sYXC7H14] * c[sH];
	w[ri2731b] = k[ri2731b] * c[sXXC7H13] * c[sH2];
	w[ri2735f] = k[ri2735f] * c[sYXC7H14] * c[sOH];
	w[ri2735b] = k[ri2735b] * c[sXXC7H13] * c[sH2O];
	w[ri2743] = k[ri2743] * c[sXXC7H13];
	w[ri3134] = k[ri3134] * c[sIXC8H18];
	w[ri3135] = k[ri3135] * c[sIXC8H18];
	w[ri3136] = k[ri3136] * c[sIXC8H18];
	w[ri3137] = k[ri3137] * c[sIXC8H18];
	w[ri3138] = k[ri3138] * c[sIXC8H18];
	w[ri3140] = k[ri3140] * c[sIXC8H18];
	w[ri3141] = k[ri3141] * c[sIXC8H18];
	w[ri3142] = k[ri3142] * c[sIXC8H18] * c[sH];
	w[ri3143] = k[ri3143] * c[sIXC8H18] * c[sH];
	w[ri3144] = k[ri3144] * c[sIXC8H18] * c[sH];
	w[ri3145] = k[ri3145] * c[sIXC8H18] * c[sH];
	w[ri3150] = k[ri3150] * c[sIXC8H18] * c[sOH];
	w[ri3151] = k[ri3151] * c[sIXC8H18] * c[sOH];
	w[ri3152] = k[ri3152] * c[sIXC8H18] * c[sOH];
	w[ri3153] = k[ri3153] * c[sIXC8H18] * c[sOH];
	w[ri3158] = k[ri3158] * c[sIXC8H18] * c[sHO2];
	w[ri3159] = k[ri3159] * c[sIXC8H18] * c[sHO2];
	w[ri3160] = k[ri3160] * c[sIXC8H18] * c[sHO2];
	w[ri3161] = k[ri3161] * c[sIXC8H18] * c[sHO2];
	w[ri3166] = k[ri3166] * c[sIXC8H18] * c[sO2];
	w[ri3167] = k[ri3167] * c[sIXC8H18] * c[sO2];
	w[ri3168] = k[ri3168] * c[sIXC8H18] * c[sO2];
	w[ri3169] = k[ri3169] * c[sIXC8H18] * c[sO2];
	w[ri3199] = k[ri3199] * c[sIXC8H18] * c[sCH3O2];
	w[ri3200] = k[ri3200] * c[sIXC8H18] * c[sCH3O2];
	w[ri3201] = k[ri3201] * c[sIXC8H18] * c[sCH3O2];
	w[ri3202] = k[ri3202] * c[sIXC8H18] * c[sCH3O2];
	w[ri3203] = k[ri3203] * c[sIXC8H18] * c[sAXC8H17O2];
	w[ri3205] = k[ri3205] * c[sIXC8H18] * c[sCXC8H17O2];
	w[ri3207] = k[ri3207] * c[sIXC8H18] * c[sAXC8H17O2];
	w[ri3209] = k[ri3209] * c[sIXC8H18] * c[sCXC8H17O2];
	w[ri3211] = k[ri3211] * c[sIXC8H18] * c[sAXC8H17O2];
	w[ri3213] = k[ri3213] * c[sIXC8H18] * c[sCXC8H17O2];
	w[ri3354] = k[ri3354] * c[sIXC8H18] * c[sAXC8H17O2];
	w[ri3356] = k[ri3356] * c[sIXC8H18] * c[sCXC8H17O2];
	w[ri3178] = k[ri3178] * c[sAXC8H17];
	w[ri3179] = k[ri3179] * c[sAXC8H17];
	w[ri3180] = k[ri3180] * c[sBXC8H17];
	w[ri3184] = k[ri3184] * c[sCXC8H17];
	w[ri3189] = k[ri3189] * c[sDXC8H17];
	w[ri3364f] = k[ri3364f] * c[sAXC8H17];
	w[ri3364b] = k[ri3364b] * c[sDXC8H17];
	w[ri3365f] = k[ri3365f] * c[sAXC8H17];
	w[ri3365b] = k[ri3365b] * c[sCXC8H17];
	w[ri3215f] = k[ri3215f] * c[sAXC8H17O2];
	w[ri3215b] = k[ri3215b] * c[sAXC8H17] * c[sO2];
	w[ri3216f] = k[ri3216f] * c[sAXC8H17O2];
	w[ri3216b] = k[ri3216b] * c[sBXC8H17] * c[sO2];
	w[ri3217f] = k[ri3217f] * c[sCXC8H17O2];
	w[ri3217b] = k[ri3217b] * c[sCXC8H17] * c[sO2];
	w[ri3218f] = k[ri3218f] * c[sAXC8H17O2];
	w[ri3218b] = k[ri3218b] * c[sDXC8H17] * c[sO2];
	w[ri3257] = k[ri3257] * c[sAXC8H17O2H];
	w[ri3259] = k[ri3259] * c[sCXC8H17O2H];
	w[ri3261] = k[ri3261] * c[sAXC8H17O];
	w[ri3263] = k[ri3263] * c[sCXC8H17O];
	w[ri3265f] = k[ri3265f] * c[sAXC8H17O2];
	w[ri3265b] = k[ri3265b] * c[sAC8H16OOHXAXC8H17O2];
	w[ri3266f] = k[ri3266f] * c[sAXC8H17O2];
	w[ri3266b] = k[ri3266b] * c[sAC8H16OOHXBXC8H17O2];
	w[ri3267f] = k[ri3267f] * c[sAXC8H17O2];
	w[ri3267b] = k[ri3267b] * c[sAC8H16OOHXCXC8H17O2];
	w[ri3272P1f] = k[ri3272P1f] * c[sCXC8H17O2];
	w[ri3272P1b] = k[ri3272P1b] * c[sCC8H16OOHXAXC8H17O2];
	w[ri3272P2] = k[ri3272P2] * c[sCXC8H17O2];
	w[ri3272P3] = k[ri3272P3] * c[sCXC8H17O2];
	w[ri3297f] = k[ri3297f] * c[sAC8H16OOHXAO2XC8H17O4];
	w[ri3297b] = k[ri3297b] * c[sAC8H16OOHXAXC8H17O2] * c[sO2];
	w[ri3298f] = k[ri3298f] * c[sAC8H16OOHXBO2XC8H17O4];
	w[ri3298b] = k[ri3298b] * c[sAC8H16OOHXBXC8H17O2] * c[sO2];
	w[ri3299f] = k[ri3299f] * c[sAC8H16OOHXCO2XC8H17O4];
	w[ri3299b] = k[ri3299b] * c[sAC8H16OOHXCXC8H17O2] * c[sO2];
	w[ri3306f] = k[ri3306f] * c[sCC8H16OOHXAO2XC8H17O4];
	w[ri3306b] = k[ri3306b] * c[sCC8H16OOHXAXC8H17O2] * c[sO2];
	w[ri3311] = k[ri3311] * c[sAC8H16OOHXAO2XC8H17O4];
	w[ri3312] = k[ri3312] * c[sAC8H16OOHXBO2XC8H17O4];
	w[ri3313] = k[ri3313] * c[sAC8H16OOHXCO2XC8H17O4];
	w[ri3367] = k[ri3367] * c[sIC8KETABXC8H16O3];
	w[ri3377] = k[ri3377] * c[sAC8H16OOHXAXC8H17O2];
	w[ri3378] = k[ri3378] * c[sAC8H16OOHXBXC8H17O2];
	w[ri3498] = k[ri3498] * c[sAC8H16OOHXAXC8H17O2];
	w[ri3506] = k[ri3506] * c[sAC8H16OOHXCXC8H17O2];
	w[r1021102f] = k[r1021102f] * c[sC6H5CH3];
	w[r1021102b] = k[r1021102b] * c[sH] * c[sC6H5CH2];
	w[r1021103] = k[r1021103] * c[sC6H5CH3];
	w[r1021104f] = k[r1021104f] * c[sC6H5CH3] * c[sO2];
	w[r1021104b] = k[r1021104b] * c[sHO2] * c[sC6H5CH2];
	w[r1021105] = k[r1021105] * c[sC6H5CH3] * c[sHO2];
	w[r1021106] = k[r1021106] * c[sC6H5CH3] * c[sH];
	w[r1021107] = k[r1021107] * c[sC6H5CH3] * c[sH];
	w[r1021108] = k[r1021108] * c[sC6H5CH3] * c[sOH];
	w[r1021109] = k[r1021109] * c[sC6H5CH3] * c[sO];
	w[r1021110] = k[r1021110] * c[sC6H5CH3] * c[sCH3];
	w[r1021111] = k[r1021111] * c[sC6H5CH3] * c[sC6H5];
	w[r1021112f] = k[r1021112f] * c[sC6H5OH] * c[sC6H5CH2];
	w[r1021112b] = k[r1021112b] * c[sC6H5CH3] * c[sC6H5O];
	w[r1021113] = k[r1021113] * c[sC6H5CH2] * c[sC6H5CHO];
	w[r1021114f] = k[r1021114f] * c[sC4H5XIXC4H5] * c[sC6H5CH3];
	w[r1021114b] = k[r1021114b] * c[sC6H5CH2] * c[sC4H6];
	w[r1021115] = k[r1021115] * c[sC6H5CH2] * c[sO2];
	w[r1021116] = k[r1021116] * c[sC6H5CH2] * c[sO2];
	w[r1021117] = k[r1021117] * c[sC6H5CH2] * c[sO];
	w[r1021118] = k[r1021118] * c[sC6H5CH2] * c[sO];
	w[r1021119] = k[r1021119] * c[sC6H5CH2] * c[sHO2];
	w[r1021121] = k[r1021121] * c[sC6H5CH2O];
	w[r1021124f] = k[r1021124f] * c[sC6H5C2H5];
	w[r1021124b] = k[r1021124b] * c[sCH3] * c[sC6H5CH2];
	w[r1021125f] = k[r1021125f] * c[sC6H5CH2] * c[sOH];
	w[r1021125b] = k[r1021125b] * c[sC6H5CH2OH];
	w[r1021126f] = k[r1021126f] * c[sC6H5CH2OH] * c[sO2];
	w[r1021126b] = k[r1021126b] * c[sHO2] * c[sC6H5CH2O];
	w[r1021127] = k[r1021127] * c[sC6H5CH2OH] * c[sOH];
	w[r1021128] = k[r1021128] * c[sC6H5CH2OH] * c[sH];
	w[r1021133] = k[r1021133] * c[sC6H5CHO] * c[sO2];
	w[r1021134] = k[r1021134] * c[sC6H5CHO] * c[sOH];
	w[r1021135] = k[r1021135] * c[sC6H5CHO] * c[sH];
	w[r1021136] = k[r1021136] * c[sC6H5CHO] * c[sH];
	w[r1021137] = k[r1021137] * c[sC6H5CHO] * c[sO];
	w[r1021138] = k[r1021138] * c[sCH3] * c[sC6H5CHO];
	w[r1021140f] = k[r1021140f] * c[sC6H5] * c[sCO];
	w[r1021140b] = k[r1021140b] * c[sC6H5CO];
	w[r1021142] = k[r1021142] * c[sC6H5] * c[sH];
	w[r1021143f] = k[r1021143f] * c[sC6H6] * c[sO2];
	w[r1021143b] = k[r1021143b] * c[sHO2] * c[sC6H5];
	w[r1021144] = k[r1021144] * c[sC6H6] * c[sOH];
	w[r1021146] = k[r1021146] * c[sC6H6] * c[sO];
	w[r1021147] = k[r1021147] * c[sC6H6] * c[sO];
	w[r1021150] = k[r1021150] * c[sC6H5] * c[sCH2O];
	w[r1021153] = k[r1021153] * c[sC6H5] * c[sHO2];
	w[r1021154f] = k[r1021154f] * c[sC6H5] * c[sO2];
	w[r1021154b] = k[r1021154b] * c[sO] * c[sC6H5O];
	w[r1021155] = k[r1021155] * c[sC6H5] * c[sO2];
	w[r1021156] = k[r1021156] * c[sC6H5O];
	w[r1021157f] = k[r1021157f] * c[sC6H5OH] * c[sO2];
	w[r1021157b] = k[r1021157b] * c[sHO2] * c[sC6H5O];
	w[r1021158] = k[r1021158] * c[sC6H5O] * c[sH];
	w[r1021159] = k[r1021159] * c[sC6H5O] * c[sO];
	w[r1021160] = k[r1021160] * c[sC6H5O] * c[sO];
	w[r1021161] = k[r1021161] * c[sC6H5O] * c[sHO2];
	w[r1021162] = k[r1021162] * c[sC6H5O] * c[sH];
	w[r1021163] = k[r1021163] * c[sC6H4O2];
	w[r1021164] = k[r1021164] * c[sC6H4O2] * c[sO];
	w[r1021165] = k[r1021165] * c[sC6H4O2] * c[sH];
	w[r1021166] = k[r1021166] * c[sC6H5OH] * c[sOH];
	w[r1021167f] = k[r1021167f] * c[sC6H6] * c[sOH];
	w[r1021167b] = k[r1021167b] * c[sH] * c[sC6H5OH];
	w[r1021168] = k[r1021168] * c[sC6H5OH] * c[sH];
	w[r1021169] = k[r1021169] * c[sC6H5OH] * c[sO];
	w[r1021170] = k[r1021170] * c[sC2H3] * c[sC6H5OH];
	w[r1021171f] = k[r1021171f] * c[sC4H5XNXC4H5] * c[sC6H5OH];
	w[r1021171b] = k[r1021171b] * c[sC6H5O] * c[sC4H6];
	w[r1021172] = k[r1021172] * c[sC6H5] * c[sC6H5OH];
	w[r1021173f] = k[r1021173f] * c[sC4H5XIXC4H5] * c[sC6H5OH];
	w[r1021173b] = k[r1021173b] * c[sC6H5O] * c[sC4H6];
	w[r1021174f] = k[r1021174f] * c[sC5H5] * c[sH];
	w[r1021174b] = k[r1021174b] * c[sC5H6];
	w[r1021175f] = k[r1021175f] * c[sC5H6] * c[sO2];
	w[r1021175b] = k[r1021175b] * c[sHO2] * c[sC5H5];
	w[r1021177] = k[r1021177] * c[sC5H6] * c[sOH];
	w[r1021178] = k[r1021178] * c[sC5H6] * c[sH];
	w[r1021179] = k[r1021179] * c[sC5H6] * c[sO];
	w[r1021180f] = k[r1021180f] * c[sC5H6] * c[sC2H3];
	w[r1021180b] = k[r1021180b] * c[sC2H4] * c[sC5H5];
	w[r1021182] = k[r1021182] * c[sC5H6] * c[sC6H5O];
	w[r1021185f] = k[r1021185f] * c[sC5H6] * c[sC6H5CH2];
	w[r1021185b] = k[r1021185b] * c[sC6H5CH3] * c[sC5H5];
	w[r1021186] = k[r1021186] * c[sC5H5] * c[sO];
	w[r1021187] = k[r1021187] * c[sC5H5] * c[sO];
	w[r1021191] = k[r1021191] * c[sC5H5] * c[sHO2];
	w[r1021192] = k[r1021192] * c[sC5H5] * c[sO2];
	w[r1021193] = k[r1021193] * c[sCH2CHCHCO] * c[sH];
	w[r1021194] = k[r1021194] * c[sCH2CHCHCO] * c[sO];
	w[r1021195] = k[r1021195] * c[sCH2CHCHCO] * c[sOH];
	w[r1021205] = k[r1021205] * c[sC5H5O];
	w[r1021206] = k[r1021206] * c[sC5H5O];
	w[r1021208] = k[r1021208] * c[sC5H4O];
	w[r1021209] = k[r1021209] * c[sC5H4O] * c[sH];
	w[r1021210] = k[r1021210] * c[sC4H5XNXC4H5];
	w[r1021211] = k[r1021211] * c[sC5H6] * c[sC4H5XIXC4H5];
	w[r1021212f] = k[r1021212f] * c[sC4H5XIXC4H5];
	w[r1021212b] = k[r1021212b] * c[sC4H5XNXC4H5];
	w[r1021218] = k[r1021218] * c[sC4H5XNXC4H5] * c[sO2];
	w[r1021244] = k[r1021244] * c[sHCCCO] * c[sO2];
	w[r1021245f] = k[r1021245f] * c[sHCCCO] * c[sH];
	w[r1021245b] = k[r1021245b] * c[sCO] * c[sC2H2];
	w[r1021246] = k[r1021246] * c[sCHCHCHO] * c[sO2];
	w[r1021247] = k[r1021247] * c[sCHCHCHO];
	w[reth139f] = k[reth139f] * c[sC2H5OH];
	w[reth139b] = k[reth139b] * c[sCH2OH] * c[sCH3];
	w[reth140f] = k[reth140f] * c[sC2H5OH];
	w[reth140b] = k[reth140b] * c[sOH] * c[sC2H5];
	w[reth141] = k[reth141] * c[sC2H5OH];
	w[reth143f] = k[reth143f] * c[sC2H5OH] * c[sOH];
	w[reth143b] = k[reth143b] * c[sH2O] * c[sC2H4OH];
	w[reth144] = k[reth144] * c[sC2H5OH] * c[sOH];
	w[reth145] = k[reth145] * c[sC2H5OH] * c[sOH];
	w[reth146f] = k[reth146f] * c[sC2H5OH] * c[sH];
	w[reth146b] = k[reth146b] * c[sH2] * c[sC2H4OH];
	w[reth147] = k[reth147] * c[sC2H5OH] * c[sH];
	w[reth148f] = k[reth148f] * c[sC2H5OH] * c[sH];
	w[reth148b] = k[reth148b] * c[sH2] * c[sCH3CH2O];
	w[reth149f] = k[reth149f] * c[sC2H5OH] * c[sO];
	w[reth149b] = k[reth149b] * c[sOH] * c[sC2H4OH];
	w[reth150] = k[reth150] * c[sC2H5OH] * c[sO];
	w[reth151f] = k[reth151f] * c[sC2H5OH] * c[sO];
	w[reth151b] = k[reth151b] * c[sOH] * c[sCH3CH2O];
	w[reth152f] = k[reth152f] * c[sC2H5OH] * c[sCH3];
	w[reth152b] = k[reth152b] * c[sCH4] * c[sC2H4OH];
	w[reth153] = k[reth153] * c[sC2H5OH] * c[sCH3];
	w[reth154f] = k[reth154f] * c[sC2H5OH] * c[sCH3];
	w[reth154b] = k[reth154b] * c[sCH4] * c[sCH3CH2O];
	w[reth155] = k[reth155] * c[sC2H5OH] * c[sHO2];
	w[reth156f] = k[reth156f] * c[sC2H5OH] * c[sHO2];
	w[reth156b] = k[reth156b] * c[sH2O2] * c[sC2H4OH];
	w[reth157] = k[reth157] * c[sC2H5OH] * c[sHO2];
	w[reth158] = k[reth158] * c[sCH3CH2O] * M[mM0];
	w[reth159] = k[reth159] * c[sCH3CH2O] * M[mM0];
	w[reth160f] = k[reth160f] * c[sCH3CH2O] * c[sO2];
	w[reth160b] = k[reth160b] * c[sHO2] * c[sCH3HCO];
	w[reth161f] = k[reth161f] * c[sCH3CH2O] * c[sCO];
	w[reth161b] = k[reth161b] * c[sCO2] * c[sC2H5];
	w[reth163] = k[reth163] * c[sCH3CH2O] * c[sH];
	w[reth164] = k[reth164] * c[sCH3CH2O] * c[sOH];
	w[reth165] = k[reth165] * c[sCH3CHOH] * c[sO2];
	w[reth167] = k[reth167] * c[sCH3CHOH] * c[sCH3];
	w[reth169] = k[reth169] * c[sCH3CHOH] * c[sH];
	w[reth171] = k[reth171] * c[sCH3CHOH] * c[sHO2];
	w[reth172f] = k[reth172f] * c[sCH3CHOH] * c[sOH];
	w[reth172b] = k[reth172b] * c[sH2O] * c[sCH3HCO];
	w[reth173] = k[reth173] * c[sCH3CHOH] * M[mM0];
	w[reth174] = k[reth174] * c[sCH3HCO] * c[sOH];
	w[reth175f] = k[reth175f] * c[sCH3HCO] * c[sOH];
	w[reth175b] = k[reth175b] * c[sH2O] * c[sCH2HCO];
	w[reth177] = k[reth177] * c[sCH3HCO] * c[sO];
	w[reth179] = k[reth179] * c[sCH3HCO] * c[sH];
	w[reth180f] = k[reth180f] * c[sCH3HCO] * c[sH];
	w[reth180b] = k[reth180b] * c[sH2] * c[sCH2HCO];
	w[reth181] = k[reth181] * c[sCH3HCO] * c[sCH3];
	w[reth182] = k[reth182] * c[sCH3HCO] * c[sCH3];
	w[reth183] = k[reth183] * c[sCH3HCO] * c[sHO2];
	w[reth184] = k[reth184] * c[sCH3HCO] * c[sHO2];
	w[reth185f] = k[reth185f] * c[sCH3HCO] * c[sO2];
	w[reth185b] = k[reth185b] * c[sHO2] * c[sCH3CO];
	w[reth200f] = k[reth200f] * c[sC2H4OH] * c[sO2];
	w[reth200b] = k[reth200b] * c[sHOC2H4O2];
	w[reth201f] = k[reth201f] * c[sHOC2H4O2];
	w[reth201b] = k[reth201b] * c[sOH] * c[sCH2O] * c[sCH2O];
	w[reth234] = k[reth234] * c[sCH2HCO] * c[sH];
	w[reth235] = k[reth235] * c[sCH2HCO] * c[sH];
	w[reth236f] = k[reth236f] * c[sCH2HCO] * c[sO];
	w[reth236b] = k[reth236b] * c[sHCO] * c[sCH2O];
	w[reth237] = k[reth237] * c[sCH2HCO] * c[sOH];
	w[reth238] = k[reth238] * c[sCH2HCO] * c[sO2];
	w[reth239] = k[reth239] * c[sCH2HCO] * c[sCH3];
	w[reth240] = k[reth240] * c[sCH2HCO] * c[sHO2];
	w[reth241f] = k[reth241f] * c[sCH2HCO] * c[sHO2];
	w[reth241b] = k[reth241b] * c[sO2] * c[sCH3HCO];
	w[reth242f] = k[reth242f] * c[sCH2HCO];
	w[reth242b] = k[reth242b] * c[sCO] * c[sCH3];
	w[reth243] = k[reth243] * c[sCH2HCO];
	w[reth260f] = k[reth260f] * c[sC2H] * c[sH2];
	w[reth260b] = k[reth260b] * c[sH] * c[sC2H2];
	w[reth261] = k[reth261] * c[sC2H] * c[sO];
	w[reth262f] = k[reth262f] * c[sC2H] * c[sOH];
	w[reth262b] = k[reth262b] * c[sH] * c[sHCCO];
	w[reth263f] = k[reth263f] * c[sC2H] * c[sO2];
	w[reth263b] = k[reth263b] * c[sH] * c[sCO] * c[sCO];
	w[reth267] = k[reth267] * c[sHCCO] * c[sO];
	w[reth271] = k[reth271] * c[sHCCO] * c[sHCCO];
	w[rhai193f] = k[rhai193f] * c[sC2H2];
	w[rhai193b] = k[rhai193b] * c[sH2CC];
	w[rhai194f] = k[rhai194f] * c[sH2CC] * c[sO2];
	w[rhai194b] = k[rhai194b] * c[sCO2] * c[sCH2];
	w[rhai195f] = k[rhai195f] * c[sH2CC] * c[sO2];
	w[rhai195b] = k[rhai195b] * c[sHCO] * c[sHCO];
	w[rhai348f] = k[rhai348f] * c[sC2H4] * c[sO2];
	w[rhai348b] = k[rhai348b] * c[sHO2] * c[sC2H3];
	w[rhai349f] = k[rhai349f] * c[sC2H4] * c[sO2];
	w[rhai349b] = k[rhai349b] * c[sH] * c[sCO2] * c[sCH3];
	w[rhai350f] = k[rhai350f] * c[sC2H4];
	w[rhai350b] = k[rhai350b] * c[sH2] * c[sH2CC];
	w[rhai352f] = k[rhai352f] * c[sC3H3] * c[sH];
	w[rhai352b] = k[rhai352b] * c[sPXC3H4];
	w[rhai353f] = k[rhai353f] * c[sC3H3] * c[sH];
	w[rhai353b] = k[rhai353b] * c[sAXC3H4];
	w[rhai354f] = k[rhai354f] * c[sC3H3] * c[sHO2];
	w[rhai354b] = k[rhai354b] * c[sO2] * c[sPXC3H4];
	w[rhai355f] = k[rhai355f] * c[sC3H3] * c[sHO2];
	w[rhai355b] = k[rhai355b] * c[sO2] * c[sAXC3H4];
	w[rhai356f] = k[rhai356f] * c[sPXC3H4] * c[sO2];
	w[rhai356b] = k[rhai356b] * c[sCO] * c[sHCO] * c[sCH3];
	w[rhai357f] = k[rhai357f] * c[sPXC3H4] * M[mM36];
	w[rhai357b] = k[rhai357b] * c[sAXC3H4] * M[mM36];
	w[rnog8f] = k[rnog8f] * c[sO] * c[s1XCH2];
	w[rnog8b] = k[rnog8b] * c[sCO] * c[sH2];
	w[rnog9f] = k[rnog9f] * c[sO] * c[s1XCH2];
	w[rnog9b] = k[rnog9b] * c[sHCO] * c[sH];
	w[rnog51f] = k[rnog51f] * c[sH] * c[s1XCH2];
	w[rnog51b] = k[rnog51b] * c[sH2] * c[sCH];
	w[rnog62f] = k[rnog62f] * c[sH] * c[sCH2OH];
	w[rnog62b] = k[rnog62b] * c[sH2O] * c[s1XCH2];
	w[rnog67f] = k[rnog67f] * c[sH] * c[sCH3O];
	w[rnog67b] = k[rnog67b] * c[sH2O] * c[s1XCH2];
	w[rnog79f] = k[rnog79f] * c[sH] * c[sHCCO];
	w[rnog79b] = k[rnog79b] * c[sCO] * c[s1XCH2];
	w[rnog94f] = k[rnog94f] * c[sOH] * c[s1XCH2];
	w[rnog94b] = k[rnog94b] * c[sCH2O] * c[sH];
	w[rnog97f] = k[rnog97f] * c[sOH] * c[sCH3];
	w[rnog97b] = k[rnog97b] * c[sH2O] * c[s1XCH2];
	w[rnog142f] = k[rnog142f] * c[s1XCH2] * c[sN2];
	w[rnog142b] = k[rnog142b] * c[sN2] * c[sCH2];
	w[rnog144f] = k[rnog144f] * c[s1XCH2] * c[sO2];
	w[rnog144b] = k[rnog144b] * c[sCO] * c[sOH] * c[sH];
	w[rnog145f] = k[rnog145f] * c[s1XCH2] * c[sO2];
	w[rnog145b] = k[rnog145b] * c[sH2O] * c[sCO];
	w[rnog146f] = k[rnog146f] * c[s1XCH2] * c[sH2];
	w[rnog146b] = k[rnog146b] * c[sH] * c[sCH3];
	w[rnog147f] = k[rnog147f] * c[s1XCH2] * c[sH2O];
	w[rnog147b] = k[rnog147b] * c[sCH3OH];
	w[rnog148f] = k[rnog148f] * c[s1XCH2] * c[sH2O];
	w[rnog148b] = k[rnog148b] * c[sH2O] * c[sCH2];
	w[rnog149f] = k[rnog149f] * c[s1XCH2] * c[sCH3];
	w[rnog149b] = k[rnog149b] * c[sC2H4] * c[sH];
	w[rnog150f] = k[rnog150f] * c[s1XCH2] * c[sCH4];
	w[rnog150b] = k[rnog150b] * c[sCH3] * c[sCH3];
	w[rnog151f] = k[rnog151f] * c[s1XCH2] * c[sCO];
	w[rnog151b] = k[rnog151b] * c[sCO] * c[sCH2];
	w[rnog152f] = k[rnog152f] * c[s1XCH2] * c[sCO2];
	w[rnog152b] = k[rnog152b] * c[sCO2] * c[sCH2];
	w[rnog153f] = k[rnog153f] * c[s1XCH2] * c[sCO2];
	w[rnog153b] = k[rnog153b] * c[sCH2O] * c[sCO];
	w[rnog154f] = k[rnog154f] * c[s1XCH2] * c[sC2H6];
	w[rnog154b] = k[rnog154b] * c[sC2H5] * c[sCH3];
	w[rnog49f] = k[rnog49f] * c[sH] * c[sCH];
	w[rnog49b] = k[rnog49b] * c[sH2] * c[sC];
	w[rnog90f] = k[rnog90f] * c[sOH] * c[sC];
	w[rnog90b] = k[rnog90b] * c[sCO] * c[sH];
	w[rnog122f] = k[rnog122f] * c[sC] * c[sO2];
	w[rnog122b] = k[rnog122b] * c[sCO] * c[sO];
	w[rnog123f] = k[rnog123f] * c[sC] * c[sCH2];
	w[rnog123b] = k[rnog123b] * c[sC2H] * c[sH];
	w[rnog124f] = k[rnog124f] * c[sC] * c[sCH3];
	w[rnog124b] = k[rnog124b] * c[sC2H2] * c[sH];
	w[rnog178f] = k[rnog178f] * c[sN] * c[sNO];
	w[rnog178b] = k[rnog178b] * c[sO] * c[sN2];
	w[rnog179f] = k[rnog179f] * c[sN] * c[sO2];
	w[rnog179b] = k[rnog179b] * c[sO] * c[sNO];
	w[rnog180f] = k[rnog180f] * c[sN] * c[sOH];
	w[rnog180b] = k[rnog180b] * c[sH] * c[sNO];
	w[rnog181f] = k[rnog181f] * c[sN2O] * c[sO];
	w[rnog181b] = k[rnog181b] * c[sO2] * c[sN2];
	w[rnog182f] = k[rnog182f] * c[sN2O] * c[sO];
	w[rnog182b] = k[rnog182b] * c[sNO] * c[sNO];
	w[rnog183f] = k[rnog183f] * c[sN2O] * c[sH];
	w[rnog183b] = k[rnog183b] * c[sOH] * c[sN2];
	w[rnog184f] = k[rnog184f] * c[sN2O] * c[sOH];
	w[rnog184b] = k[rnog184b] * c[sHO2] * c[sN2];
	w[rnog185f] = k[rnog185f] * c[sN2O];
	w[rnog185b] = k[rnog185b] * c[sO] * c[sN2];
	w[rnog186f] = k[rnog186f] * c[sHO2] * c[sNO];
	w[rnog186b] = k[rnog186b] * c[sOH] * c[sNO2];
	w[rnog187f] = k[rnog187f] * c[sNO] * c[sO] * M[mM40];
	w[rnog187b] = k[rnog187b] * c[sNO2] * M[mM40];
	w[rnog188f] = k[rnog188f] * c[sNO2] * c[sO];
	w[rnog188b] = k[rnog188b] * c[sO2] * c[sNO];
	w[rnog189f] = k[rnog189f] * c[sNO2] * c[sH];
	w[rnog189b] = k[rnog189b] * c[sOH] * c[sNO];
	w[rnog190f] = k[rnog190f] * c[sNH] * c[sO];
	w[rnog190b] = k[rnog190b] * c[sH] * c[sNO];
	w[rnog191f] = k[rnog191f] * c[sNH] * c[sH];
	w[rnog191b] = k[rnog191b] * c[sH2] * c[sN];
	w[rnog192f] = k[rnog192f] * c[sNH] * c[sOH];
	w[rnog192b] = k[rnog192b] * c[sH] * c[sHNO];
	w[rnog193f] = k[rnog193f] * c[sNH] * c[sOH];
	w[rnog193b] = k[rnog193b] * c[sH2O] * c[sN];
	w[rnog194f] = k[rnog194f] * c[sNH] * c[sO2];
	w[rnog194b] = k[rnog194b] * c[sO] * c[sHNO];
	w[rnog195f] = k[rnog195f] * c[sNH] * c[sO2];
	w[rnog195b] = k[rnog195b] * c[sOH] * c[sNO];
	w[rnog196f] = k[rnog196f] * c[sNH] * c[sN];
	w[rnog196b] = k[rnog196b] * c[sH] * c[sN2];
	w[rnog197f] = k[rnog197f] * c[sNH] * c[sH2O];
	w[rnog197b] = k[rnog197b] * c[sH2] * c[sHNO];
	w[rnog198f] = k[rnog198f] * c[sNH] * c[sNO];
	w[rnog198b] = k[rnog198b] * c[sOH] * c[sN2];
	w[rnog199f] = k[rnog199f] * c[sNH] * c[sNO];
	w[rnog199b] = k[rnog199b] * c[sH] * c[sN2O];
	w[rnog200f] = k[rnog200f] * c[sNH2] * c[sO];
	w[rnog200b] = k[rnog200b] * c[sNH] * c[sOH];
	w[rnog201f] = k[rnog201f] * c[sNH2] * c[sO];
	w[rnog201b] = k[rnog201b] * c[sHNO] * c[sH];
	w[rnog202f] = k[rnog202f] * c[sNH2] * c[sH];
	w[rnog202b] = k[rnog202b] * c[sH2] * c[sNH];
	w[rnog203f] = k[rnog203f] * c[sNH2] * c[sOH];
	w[rnog203b] = k[rnog203b] * c[sH2O] * c[sNH];
	w[rnog204f] = k[rnog204f] * c[sNNH];
	w[rnog204b] = k[rnog204b] * c[sH] * c[sN2];
	w[rnog205f] = k[rnog205f] * c[sNNH] * M[mM40];
	w[rnog205b] = k[rnog205b] * c[sH] * c[sN2] * M[mM40];
	w[rnog206f] = k[rnog206f] * c[sNNH] * c[sO2];
	w[rnog206b] = k[rnog206b] * c[sN2] * c[sHO2];
	w[rnog207f] = k[rnog207f] * c[sNNH] * c[sO];
	w[rnog207b] = k[rnog207b] * c[sN2] * c[sOH];
	w[rnog208f] = k[rnog208f] * c[sNNH] * c[sO];
	w[rnog208b] = k[rnog208b] * c[sNO] * c[sNH];
	w[rnog209f] = k[rnog209f] * c[sNNH] * c[sH];
	w[rnog209b] = k[rnog209b] * c[sN2] * c[sH2];
	w[rnog210f] = k[rnog210f] * c[sNNH] * c[sOH];
	w[rnog210b] = k[rnog210b] * c[sN2] * c[sH2O];
	w[rnog211f] = k[rnog211f] * c[sNNH] * c[sCH3];
	w[rnog211b] = k[rnog211b] * c[sN2] * c[sCH4];
	w[rnog212f] = k[rnog212f] * c[sH] * c[sNO] * M[mM40];
	w[rnog212b] = k[rnog212b] * c[sHNO] * M[mM40];
	w[rnog213f] = k[rnog213f] * c[sHNO] * c[sO];
	w[rnog213b] = k[rnog213b] * c[sOH] * c[sNO];
	w[rnog214f] = k[rnog214f] * c[sHNO] * c[sH];
	w[rnog214b] = k[rnog214b] * c[sNO] * c[sH2];
	w[rnog215f] = k[rnog215f] * c[sHNO] * c[sOH];
	w[rnog215b] = k[rnog215b] * c[sH2O] * c[sNO];
	w[rnog216f] = k[rnog216f] * c[sHNO] * c[sO2];
	w[rnog216b] = k[rnog216b] * c[sNO] * c[sHO2];
	w[rnog217f] = k[rnog217f] * c[sCN] * c[sO];
	w[rnog217b] = k[rnog217b] * c[sN] * c[sCO];
	w[rnog218f] = k[rnog218f] * c[sCN] * c[sOH];
	w[rnog218b] = k[rnog218b] * c[sH] * c[sNCO];
	w[rnog219f] = k[rnog219f] * c[sCN] * c[sH2O];
	w[rnog219b] = k[rnog219b] * c[sOH] * c[sHCN];
	w[rnog220f] = k[rnog220f] * c[sCN] * c[sO2];
	w[rnog220b] = k[rnog220b] * c[sO] * c[sNCO];
	w[rnog221f] = k[rnog221f] * c[sCN] * c[sH2];
	w[rnog221b] = k[rnog221b] * c[sH] * c[sHCN];
	w[rnog222f] = k[rnog222f] * c[sNCO] * c[sO];
	w[rnog222b] = k[rnog222b] * c[sCO] * c[sNO];
	w[rnog223f] = k[rnog223f] * c[sNCO] * c[sH];
	w[rnog223b] = k[rnog223b] * c[sCO] * c[sNH];
	w[rnog224f] = k[rnog224f] * c[sNCO] * c[sOH];
	w[rnog224b] = k[rnog224b] * c[sCO] * c[sH] * c[sNO];
	w[rnog225f] = k[rnog225f] * c[sNCO] * c[sN];
	w[rnog225b] = k[rnog225b] * c[sCO] * c[sN2];
	w[rnog226f] = k[rnog226f] * c[sNCO] * c[sO2];
	w[rnog226b] = k[rnog226b] * c[sCO2] * c[sNO];
	w[rnog227f] = k[rnog227f] * c[sNCO] * M[mM40];
	w[rnog227b] = k[rnog227b] * c[sCO] * c[sN] * M[mM40];
	w[rnog228f] = k[rnog228f] * c[sNCO] * c[sNO];
	w[rnog228b] = k[rnog228b] * c[sCO] * c[sN2O];
	w[rnog229f] = k[rnog229f] * c[sNCO] * c[sNO];
	w[rnog229b] = k[rnog229b] * c[sCO2] * c[sN2];
	w[rnog230f] = k[rnog230f] * c[sHCN] * M[mM40];
	w[rnog230b] = k[rnog230b] * c[sCN] * c[sH] * M[mM40];
	w[rnog231f] = k[rnog231f] * c[sHCN] * c[sO];
	w[rnog231b] = k[rnog231b] * c[sH] * c[sNCO];
	w[rnog232f] = k[rnog232f] * c[sHCN] * c[sO];
	w[rnog232b] = k[rnog232b] * c[sCO] * c[sNH];
	w[rnog233f] = k[rnog233f] * c[sHCN] * c[sO];
	w[rnog233b] = k[rnog233b] * c[sOH] * c[sCN];
	w[rnog234f] = k[rnog234f] * c[sHCN] * c[sOH];
	w[rnog234b] = k[rnog234b] * c[sH] * c[sHOCN];
	w[rnog235] = k[rnog235] * c[sHCN] * c[sOH];
	w[rnog236f] = k[rnog236f] * c[sHCN] * c[sOH];
	w[rnog236b] = k[rnog236b] * c[sCO] * c[sNH2];
	w[rnog237f] = k[rnog237f] * c[sH] * c[sHCN] * M[mM40];
	w[rnog237b] = k[rnog237b] * c[sH2CN] * M[mM40];
	w[rnog238f] = k[rnog238f] * c[sH2CN] * c[sN];
	w[rnog238b] = k[rnog238b] * c[sCH2] * c[sN2];
	w[rnog239f] = k[rnog239f] * c[sC] * c[sN2];
	w[rnog239b] = k[rnog239b] * c[sN] * c[sCN];
	w[rnog240f] = k[rnog240f] * c[sCH] * c[sN2];
	w[rnog240b] = k[rnog240b] * c[sN] * c[sHCN];
	w[rnog241f] = k[rnog241f] * c[sCH] * c[sN2];
	w[rnog241b] = k[rnog241b] * c[sHCNN];
	w[rnog242f] = k[rnog242f] * c[sCH2] * c[sN2];
	w[rnog242b] = k[rnog242b] * c[sNH] * c[sHCN];
	w[rnog243f] = k[rnog243f] * c[s1XCH2] * c[sN2];
	w[rnog243b] = k[rnog243b] * c[sHCN] * c[sNH];
	w[rnog244f] = k[rnog244f] * c[sC] * c[sNO];
	w[rnog244b] = k[rnog244b] * c[sO] * c[sCN];
	w[rnog245f] = k[rnog245f] * c[sC] * c[sNO];
	w[rnog245b] = k[rnog245b] * c[sN] * c[sCO];
	w[rnog246f] = k[rnog246f] * c[sCH] * c[sNO];
	w[rnog246b] = k[rnog246b] * c[sO] * c[sHCN];
	w[rnog247f] = k[rnog247f] * c[sCH] * c[sNO];
	w[rnog247b] = k[rnog247b] * c[sNCO] * c[sH];
	w[rnog248f] = k[rnog248f] * c[sCH] * c[sNO];
	w[rnog248b] = k[rnog248b] * c[sHCO] * c[sN];
	w[rnog249f] = k[rnog249f] * c[sCH2] * c[sNO];
	w[rnog249b] = k[rnog249b] * c[sHNCO] * c[sH];
	w[rnog250f] = k[rnog250f] * c[sCH2] * c[sNO];
	w[rnog250b] = k[rnog250b] * c[sHCN] * c[sOH];
	w[rnog251f] = k[rnog251f] * c[sCH2] * c[sNO];
	w[rnog251b] = k[rnog251b] * c[sHCNO] * c[sH];
	w[rnog252f] = k[rnog252f] * c[s1XCH2] * c[sNO];
	w[rnog252b] = k[rnog252b] * c[sHNCO] * c[sH];
	w[rnog253f] = k[rnog253f] * c[s1XCH2] * c[sNO];
	w[rnog253b] = k[rnog253b] * c[sHCN] * c[sOH];
	w[rnog254f] = k[rnog254f] * c[s1XCH2] * c[sNO];
	w[rnog254b] = k[rnog254b] * c[sHCNO] * c[sH];
	w[rnog255f] = k[rnog255f] * c[sCH3] * c[sNO];
	w[rnog255b] = k[rnog255b] * c[sH2O] * c[sHCN];
	w[rnog256f] = k[rnog256f] * c[sCH3] * c[sNO];
	w[rnog256b] = k[rnog256b] * c[sOH] * c[sH2CN];
	w[rnog257f] = k[rnog257f] * c[sHCNN] * c[sO];
	w[rnog257b] = k[rnog257b] * c[sN2] * c[sH] * c[sCO];
	w[rnog258f] = k[rnog258f] * c[sHCNN] * c[sO];
	w[rnog258b] = k[rnog258b] * c[sNO] * c[sHCN];
	w[rnog259f] = k[rnog259f] * c[sHCNN] * c[sO2];
	w[rnog259b] = k[rnog259b] * c[sN2] * c[sHCO] * c[sO];
	w[rnog260f] = k[rnog260f] * c[sHCNN] * c[sOH];
	w[rnog260b] = k[rnog260b] * c[sN2] * c[sHCO] * c[sH];
	w[rnog261f] = k[rnog261f] * c[sHCNN] * c[sH];
	w[rnog261b] = k[rnog261b] * c[sN2] * c[sCH2];
	w[rnog262f] = k[rnog262f] * c[sHNCO] * c[sO];
	w[rnog262b] = k[rnog262b] * c[sCO2] * c[sNH];
	w[rnog263f] = k[rnog263f] * c[sHNCO] * c[sO];
	w[rnog263b] = k[rnog263b] * c[sCO] * c[sHNO];
	w[rnog264f] = k[rnog264f] * c[sHNCO] * c[sO];
	w[rnog264b] = k[rnog264b] * c[sOH] * c[sNCO];
	w[rnog265f] = k[rnog265f] * c[sHNCO] * c[sH];
	w[rnog265b] = k[rnog265b] * c[sCO] * c[sNH2];
	w[rnog266f] = k[rnog266f] * c[sHNCO] * c[sH];
	w[rnog266b] = k[rnog266b] * c[sNCO] * c[sH2];
	w[rnog267f] = k[rnog267f] * c[sHNCO] * c[sOH];
	w[rnog267b] = k[rnog267b] * c[sH2O] * c[sNCO];
	w[rnog268f] = k[rnog268f] * c[sHNCO] * c[sOH];
	w[rnog268b] = k[rnog268b] * c[sCO2] * c[sNH2];
	w[rnog269f] = k[rnog269f] * c[sHNCO] * M[mM40];
	w[rnog269b] = k[rnog269b] * c[sCO] * c[sNH] * M[mM40];
	w[rnog270f] = k[rnog270f] * c[sHCNO] * c[sH];
	w[rnog270b] = k[rnog270b] * c[sHNCO] * c[sH];
	w[rnog271f] = k[rnog271f] * c[sHCNO] * c[sH];
	w[rnog271b] = k[rnog271b] * c[sHCN] * c[sOH];
	w[rnog272f] = k[rnog272f] * c[sHCNO] * c[sH];
	w[rnog272b] = k[rnog272b] * c[sCO] * c[sNH2];
	w[rnog273] = k[rnog273] * c[sHOCN] * c[sH];
	w[rnog274f] = k[rnog274f] * c[sHCCO] * c[sNO];
	w[rnog274b] = k[rnog274b] * c[sCO] * c[sHCNO];
	w[rnog275f] = k[rnog275f] * c[sCH3] * c[sN];
	w[rnog275b] = k[rnog275b] * c[sH] * c[sH2CN];
	w[rnog276f] = k[rnog276f] * c[sCH3] * c[sN];
	w[rnog276b] = k[rnog276b] * c[sH2] * c[sHCN];
	w[rnog277f] = k[rnog277f] * c[sNH3] * c[sH];
	w[rnog277b] = k[rnog277b] * c[sH2] * c[sNH2];
	w[rnog278f] = k[rnog278f] * c[sNH3] * c[sOH];
	w[rnog278b] = k[rnog278b] * c[sH2O] * c[sNH2];
	w[rnog279f] = k[rnog279f] * c[sNH3] * c[sO];
	w[rnog279b] = k[rnog279b] * c[sOH] * c[sNH2];
	w[rnog280f] = k[rnog280f] * c[sCH3] * c[sNO2];
	w[rnog280b] = k[rnog280b] * c[sNO] * c[sCH3O];


	cdot[sN2] = - w[rnog142f] + w[rnog142f] - w[rnog142b] + w[rnog142b]
		+ w[rnog178f] - w[rnog178b] + w[rnog181f] - w[rnog181b]
		+ w[rnog183f] - w[rnog183b] + w[rnog184f] - w[rnog184b]
		+ w[rnog185f] - w[rnog185b] + w[rnog196f] - w[rnog196b]
		+ w[rnog198f] - w[rnog198b] + w[rnog204f] - w[rnog204b]
		+ w[rnog205f] - w[rnog205b] + w[rnog206f] - w[rnog206b]
		+ w[rnog207f] - w[rnog207b] + w[rnog209f] - w[rnog209b]
		+ w[rnog210f] - w[rnog210b] + w[rnog211f] - w[rnog211b]
		+ w[rnog225f] - w[rnog225b] + w[rnog229f] - w[rnog229b]
		+ w[rnog238f] - w[rnog238b] - w[rnog239f] + w[rnog239b]
		- w[rnog240f] + w[rnog240b] - w[rnog241f] + w[rnog241b]
		- w[rnog242f] + w[rnog242b] - w[rnog243f] + w[rnog243b]
		+ w[rnog257f] - w[rnog257b] + w[rnog259f] - w[rnog259b]
		+ w[rnog260f] - w[rnog260b] + w[rnog261f] - w[rnog261b];

	cdot[sAR] = 0.0;

	cdot[sC6H5CH3] = - w[r1021102f] + w[r1021102b] - w[r1021103] - w[r1021104f]
		+ w[r1021104b] - w[r1021105] - w[r1021106] - w[r1021107]
		- w[r1021108] - w[r1021109] - w[r1021110] - w[r1021111]
		+ w[r1021112f] - w[r1021112b] + w[r1021113] - w[r1021114f]
		+ w[r1021114b] + w[r1021185f] - w[r1021185b];

	cdot[sHO2] = - w[r49] - w[r535] + w[ri26f] - w[ri26b]
		+ w[ri13P1f] - w[ri13P1b] - w[ri22f] + w[ri22b]
		- w[ri23f] + w[ri23b] - w[ri42f] + w[ri42b]
		- w[ri43f] + w[ri43b] + w[ri44] + w[ri46]
		- w[ri47] - w[ri48] - w[ri54] + w[ri83]
		- w[ri570] + w[r85] - w[r142] + w[r171]
		- w[r536] - w[r1017] - w[r1190f] + w[r1190b]
		- w[r1335] + w[r1401f] - w[r1401b] - w[r2092]
		- w[r2093] - w[r2094] - w[r2095] + w[r2100f]
		- w[r2100b] + w[r2101f] - w[r2101b] + w[r2102f]
		- w[r2102b] + w[r2103f] - w[r2103b] + w[r2258f]
		- w[r2258b] + w[ri19f] - w[ri19b] - w[ri45]
		+ 2 * w[ri50P1f] - 2 * w[ri50P1b] - w[ri53f] + w[ri53b]
		+ w[ri182] - w[ri217] + w[ri340f] - w[ri340b]
		+ 2 * w[ri50P2f] - 2 * w[ri50P2b] + w[ri13P2f] - w[ri13P2b]
		- w[ri710] - w[ri711f] + w[ri711b] + w[ri790]
		+ w[ri791] - w[ri1394] - w[ri1461] - w[ri1508]
		- w[ri1511] + w[ri1517] + w[ri1536f] - w[ri1536b]
		+ w[ri2727f] - w[ri2727b] - w[ri3158] - w[ri3159]
		- w[ri3160] - w[ri3161] + w[ri3166] + w[ri3167]
		+ w[ri3168] + w[ri3169] + w[r1021104f] - w[r1021104b]
		- w[r1021105] - w[r1021119] + w[r1021126f] - w[r1021126b]
		+ w[r1021133] + w[r1021143f] - w[r1021143b] - w[r1021153]
		+ w[r1021157f] - w[r1021157b] - w[r1021161] + w[r1021175f]
		- w[r1021175b] - w[r1021191] + w[r1021246] - w[reth155]
		- w[reth156f] + w[reth156b] - w[reth157] + w[reth160f]
		- w[reth160b] + w[reth165] - w[reth171] - w[reth183]
		- w[reth184] + w[reth185f] - w[reth185b] - w[reth240]
		- w[reth241f] + w[reth241b] + w[rhai348f] - w[rhai348b]
		- w[rhai354f] + w[rhai354b] - w[rhai355f] + w[rhai355b]
		+ w[rnog184f] - w[rnog184b] - w[rnog186f] + w[rnog186b]
		+ w[rnog206f] - w[rnog206b] + w[rnog216f] - w[rnog216b];

	cdot[sOH] = - w[r49] - w[r62f] + w[r62b] - w[ri7f]
		+ w[ri7b] + w[ri8f] - w[ri8b] + 2 * w[ri10f]
		- 2 * w[ri10b] - w[ri13P1f] + w[ri13P1b] + w[ri22f]
		- w[ri22b] + w[ri23f] - w[ri23b] + w[ri42f]
		- w[ri42b] + 2 * w[ri47] + w[ri52] + w[ri72]
		- w[ri94] - w[ri159] + w[ri522] + w[ri570]
		- w[ri1457] - w[r16f] + w[r16b] + w[r63]
		- w[r91] - w[r174] - w[r335] + w[r536]
		+ w[r1018] + w[r1190f] - w[r1190b] + w[r1335]
		- w[r1452] - w[r2088] - w[r2089] - w[r2090]
		- w[r2091] + w[r2281] + w[r2309] + w[r2327]
		- w[ri3f] + w[ri3b] + w[ri4f] - w[ri4b]
		- w[ri6f] + w[ri6b] + w[ri9f] - w[ri9b]
		- w[ri11f] + w[ri11b] - w[ri20f] + w[ri20b]
		+ w[ri21f] - w[ri21b] + w[ri25f] - w[ri25b]
		+ w[ri30f] - w[ri30b] - w[ri32f] + w[ri32b]
		+ w[ri34f] - w[ri34b] - w[ri35f] + w[ri35b]
		- 2 * w[ri51f] + 2 * w[ri51b] - w[ri75f] + w[ri75b]
		- w[ri88] - w[ri96] - w[ri173] + w[ri174]
		+ w[ri176] - w[ri213] - w[ri309] - w[ri344f]
		+ w[ri344b] - w[ri13P2f] + w[ri13P2b] + w[ri710]
		+ w[ri711f] - w[ri711b] + w[ri779] + w[ri780]
		+ w[ri793] + w[ri794] + w[ri797] + w[ri1394]
		+ w[ri1412] + w[ri1458] + w[ri1461] + w[ri1468]
		- w[ri1510] - w[ri2735f] + w[ri2735b] - w[ri3150]
		- w[ri3151] - w[ri3152] - w[ri3153] + w[ri3257]
		+ w[ri3259] + w[ri3311] + w[ri3312] + w[ri3313]
		+ w[ri3367] + w[ri3377] + w[ri3378] - w[r1021108]
		+ w[r1021109] + w[r1021115] + w[r1021119] - w[r1021125f]
		+ w[r1021125b] - w[r1021127] - w[r1021134] + w[r1021137]
		- w[r1021144] + w[r1021153] + w[r1021161] - w[r1021166]
		- w[r1021167f] + w[r1021167b] + w[r1021169] - w[r1021177]
		+ w[r1021179] + w[r1021191] - w[r1021195] + w[reth140f]
		- w[reth140b] - w[reth143f] + w[reth143b] - w[reth144]
		- w[reth145] + w[reth149f] - w[reth149b] + w[reth150]
		+ w[reth151f] - w[reth151b] - w[reth164] + 2 * w[reth171]
		- w[reth172f] + w[reth172b] - w[reth174] - w[reth175f]
		+ w[reth175b] + w[reth177] + w[reth201f] - w[reth201b]
		- w[reth237] + w[reth238] + w[reth240] - w[reth262f]
		+ w[reth262b] - w[rnog94f] + w[rnog94b] - w[rnog97f]
		+ w[rnog97b] + w[rnog144f] - w[rnog144b] - w[rnog90f]
		+ w[rnog90b] - w[rnog180f] + w[rnog180b] + w[rnog183f]
		- w[rnog183b] - w[rnog184f] + w[rnog184b] + w[rnog186f]
		- w[rnog186b] + w[rnog189f] - w[rnog189b] - w[rnog192f]
		+ w[rnog192b] - w[rnog193f] + w[rnog193b] + w[rnog195f]
		- w[rnog195b] + w[rnog198f] - w[rnog198b] + w[rnog200f]
		- w[rnog200b] - w[rnog203f] + w[rnog203b] + w[rnog207f]
		- w[rnog207b] - w[rnog210f] + w[rnog210b] + w[rnog213f]
		- w[rnog213b] - w[rnog215f] + w[rnog215b] - w[rnog218f]
		+ w[rnog218b] + w[rnog219f] - w[rnog219b] - w[rnog224f]
		+ w[rnog224b] + w[rnog233f] - w[rnog233b] - w[rnog234f]
		+ w[rnog234b] - w[rnog235] - w[rnog236f] + w[rnog236b]
		+ w[rnog250f] - w[rnog250b] + w[rnog253f] - w[rnog253b]
		+ w[rnog256f] - w[rnog256b] - w[rnog260f] + w[rnog260b]
		+ w[rnog264f] - w[rnog264b] - w[rnog267f] + w[rnog267b]
		- w[rnog268f] + w[rnog268b] + w[rnog271f] - w[rnog271b]
		- w[rnog278f] + w[rnog278b] + w[rnog279f] - w[rnog279b];

	cdot[sH2O] = w[r49] + w[r62f] - w[r62b] - w[ri10f]
		+ w[ri10b] + w[ri13P1f] - w[ri13P1b] + w[ri52]
		+ w[ri94] + w[ri159] + w[r335] + w[r2088]
		+ w[r2089] + w[r2090] + w[r2091] + w[ri3f]
		- w[ri3b] + w[ri6f] - w[ri6b] + w[ri11f]
		- w[ri11b] + w[ri20f] - w[ri20b] - w[ri25f]
		+ w[ri25b] + w[ri32f] - w[ri32b] + w[ri66]
		+ w[ri75f] - w[ri75b] + w[ri173] + w[ri213]
		+ w[ri309] + w[ri344f] - w[ri344b] + w[ri13P2f]
		- w[ri13P2b] + w[ri1510] + w[ri2735f] - w[ri2735b]
		+ w[ri3150] + w[ri3151] + w[ri3152] + w[ri3153]
		+ w[r1021108] + w[r1021127] + w[r1021134] + w[r1021144]
		+ w[r1021166] + w[r1021177] + w[reth141] + w[reth143f]
		- w[reth143b] + w[reth144] + w[reth145] + w[reth163]
		+ w[reth164] + w[reth167] + w[reth169] + w[reth172f]
		- w[reth172b] + w[reth174] + w[reth175f] - w[reth175b]
		+ w[reth237] + w[rnog62f] - w[rnog62b] + w[rnog67f]
		- w[rnog67b] + w[rnog97f] - w[rnog97b] + w[rnog145f]
		- w[rnog145b] - w[rnog147f] + w[rnog147b] - w[rnog148f]
		+ w[rnog148f] - w[rnog148b] + w[rnog148b] + w[rnog193f]
		- w[rnog193b] - w[rnog197f] + w[rnog197b] + w[rnog203f]
		- w[rnog203b] + w[rnog210f] - w[rnog210b] + w[rnog215f]
		- w[rnog215b] - w[rnog219f] + w[rnog219b] + w[rnog255f]
		- w[rnog255b] + w[rnog267f] - w[rnog267b] + w[rnog278f]
		- w[rnog278b];

	cdot[sO2] = w[r49] - w[r112] + w[r535] - w[ri26f]
		+ w[ri26b] - w[ri8f] + w[ri8b] - w[ri37f]
		+ w[ri37b] + w[ri42f] - w[ri42b] + w[ri43f]
		- w[ri43b] - w[ri44] - w[ri46] + w[ri48]
		- w[ri72] - w[ri83] + w[ri521f] - w[ri521b]
		- w[ri1455] - w[r65] - w[r85] - w[r171]
		+ w[r536] + w[r538] + w[r1014f] - w[r1014b]
		+ w[r1017] - w[r1401f] + w[r1401b] - w[r2100f]
		+ w[r2100b] - w[r2101f] + w[r2101b] - w[r2102f]
		+ w[r2102b] - w[r2103f] + w[r2103b] + w[r2188f]
		- w[r2188b] + w[r2189f] - w[r2189b] + w[r2190f]
		- w[r2190b] + w[r2191f] - w[r2191b] + w[r2290f]
		- w[r2290b] + w[r2291f] - w[r2291b] - w[ri19f]
		+ w[ri19b] - w[ri28f] + w[ri28b] + w[ri45]
		- w[ri50P1f] + w[ri50P1b] - w[ri66] - w[ri76]
		- w[ri77f] + w[ri77b] - w[ri91] - w[ri101]
		- w[ri182] - w[ri50P2f] + w[ri50P2b] + w[ri540]
		+ w[ri776f] - w[ri776b] + w[ri778f] - w[ri778b]
		- w[ri1347] + w[ri1393f] - w[ri1393b] + w[ri1396f]
		- w[ri1396b] - w[ri1458] - w[ri1460] + w[ri1508]
		- w[ri1517] - w[ri1536f] + w[ri1536b] - w[ri2727f]
		+ w[ri2727b] - w[ri3166] - w[ri3167] - w[ri3168]
		- w[ri3169] + w[ri3215f] - w[ri3215b] + w[ri3216f]
		- w[ri3216b] + w[ri3217f] - w[ri3217b] + w[ri3218f]
		- w[ri3218b] + w[ri3297f] - w[ri3297b] + w[ri3298f]
		- w[ri3298b] + w[ri3299f] - w[ri3299b] + w[ri3306f]
		- w[ri3306b] - w[r1021104f] + w[r1021104b] - w[r1021115]
		- w[r1021116] - w[r1021126f] + w[r1021126b] - w[r1021133]
		- w[r1021143f] + w[r1021143b] - w[r1021154f] + w[r1021154b]
		- w[r1021155] - w[r1021157f] + w[r1021157b] - w[r1021175f]
		+ w[r1021175b] - w[r1021192] - w[r1021218] - w[r1021244]
		- w[r1021246] - w[reth160f] + w[reth160b] - w[reth165]
		- w[reth185f] + w[reth185b] - w[reth200f] + w[reth200b]
		- w[reth238] + w[reth241f] - w[reth241b] - w[reth263f]
		+ w[reth263b] - w[rhai194f] + w[rhai194b] - w[rhai195f]
		+ w[rhai195b] - w[rhai348f] + w[rhai348b] - w[rhai349f]
		+ w[rhai349b] + w[rhai354f] - w[rhai354b] + w[rhai355f]
		- w[rhai355b] - w[rhai356f] + w[rhai356b] - w[rnog144f]
		+ w[rnog144b] - w[rnog145f] + w[rnog145b] - w[rnog122f]
		+ w[rnog122b] - w[rnog179f] + w[rnog179b] + w[rnog181f]
		- w[rnog181b] + w[rnog188f] - w[rnog188b] - w[rnog194f]
		+ w[rnog194b] - w[rnog195f] + w[rnog195b] - w[rnog206f]
		+ w[rnog206b] - w[rnog216f] + w[rnog216b] - w[rnog220f]
		+ w[rnog220b] - w[rnog226f] + w[rnog226b] - w[rnog259f]
		+ w[rnog259b];

	cdot[sC2H4] = - w[r61f] + w[r61b] - w[r62f] + w[r62b]
		- w[r63] - w[r135] + w[r169] + w[r255]
		+ w[r464] + w[r2142] + w[r2327] - w[ri14f]
		+ w[ri14b] - w[ri15f] + w[ri15b] + w[ri19f]
		- w[ri19b] + w[ri179] + w[r1021170] + w[r1021180f]
		- w[r1021180b] + w[reth141] + w[reth163] + w[reth169]
		- w[rhai348f] + w[rhai348b] - w[rhai349f] + w[rhai349b]
		- w[rhai350f] + w[rhai350b] + w[rnog149f] - w[rnog149b];

	cdot[sH] = - w[r61f] + w[r61b] + w[r113f] - w[r113b]
		- w[ri26f] + w[ri26b] + w[ri7f] - w[ri7b]
		- w[ri8f] + w[ri8b] - w[ri17f] + w[ri17b]
		- w[ri29f] + w[ri29b] + w[ri31f] - w[ri31b]
		- w[ri33f] + w[ri33b] + w[ri36f] - w[ri36b]
		+ w[ri40] - w[ri47] - w[ri48] - w[ri52]
		- w[ri60f] + w[ri60b] + w[ri84f] - w[ri84b]
		+ w[ri87] - w[ri111f] + w[ri111b] - w[ri177]
		+ w[r170f] - w[r170b] - w[r331] + w[r465f]
		- w[r465b] + w[r1400f] - w[r1400b] - w[r2080]
		- w[r2081] - w[r2082] - w[r2083] - w[rnog52f]
		+ w[rnog52b] - w[ri2f] + w[ri2b] + w[ri9f]
		- w[ri9b] + w[ri11f] - w[ri11b] + w[ri12f]
		- w[ri12b] - w[ri15f] + w[ri15b] + w[ri25f]
		- w[ri25b] - w[ri74f] + w[ri74b] + 2 * w[ri76]
		+ w[ri88] - w[ri89] - w[ri95] + w[ri98]
		+ w[ri151f] - w[ri151b] - w[ri179] + w[ri180f]
		- w[ri180b] - w[ri181] + w[ri207f] - w[ri207b]
		- w[ri211] + w[ri259f] - w[ri259b] + w[ri263]
		- w[ri310] - w[ri340f] + w[ri340b] - w[ri343f]
		+ w[ri343b] - w[ri347] + w[ri1463] + w[ri1466]
		+ w[ri1546] + w[ri2722f] - w[ri2722b] - w[ri2731f]
		+ w[ri2731b] + w[ri3134] + w[ri3135] + w[ri3136]
		+ w[ri3137] - w[ri3142] - w[ri3143] - w[ri3144]
		- w[ri3145] + w[r1021102f] - w[r1021102b] - w[r1021106]
		- w[r1021107] + w[r1021117] + w[r1021121] - w[r1021128]
		- w[r1021135] - w[r1021136] - w[r1021142] + w[r1021155]
		- w[r1021158] + w[r1021159] - w[r1021162] - w[r1021165]
		+ w[r1021167f] - w[r1021167b] - w[r1021168] - w[r1021174f]
		+ w[r1021174b] - w[r1021178] + w[r1021187] - w[r1021193]
		+ w[r1021205] - w[r1021209] - w[r1021245f] + w[r1021245b]
		- w[reth146f] + w[reth146b] - w[reth147] - w[reth148f]
		+ w[reth148b] + w[reth158] - w[reth163] - w[reth169]
		+ w[reth173] - w[reth179] - w[reth180f] + w[reth180b]
		- w[reth234] - w[reth235] + w[reth239] + w[reth243]
		+ w[reth260f] - w[reth260b] + w[reth262f] - w[reth262b]
		+ w[reth263f] - w[reth263b] + w[rhai349f] - w[rhai349b]
		- w[rhai352f] + w[rhai352b] - w[rhai353f] + w[rhai353b]
		+ w[rnog9f] - w[rnog9b] - w[rnog51f] + w[rnog51b]
		- w[rnog62f] + w[rnog62b] - w[rnog67f] + w[rnog67b]
		- w[rnog79f] + w[rnog79b] + w[rnog94f] - w[rnog94b]
		+ w[rnog144f] - w[rnog144b] + w[rnog146f] - w[rnog146b]
		+ w[rnog149f] - w[rnog149b] - w[rnog49f] + w[rnog49b]
		+ w[rnog90f] - w[rnog90b] + w[rnog123f] - w[rnog123b]
		+ w[rnog124f] - w[rnog124b] + w[rnog180f] - w[rnog180b]
		- w[rnog183f] + w[rnog183b] - w[rnog189f] + w[rnog189b]
		+ w[rnog190f] - w[rnog190b] - w[rnog191f] + w[rnog191b]
		+ w[rnog192f] - w[rnog192b] + w[rnog196f] - w[rnog196b]
		+ w[rnog199f] - w[rnog199b] + w[rnog201f] - w[rnog201b]
		- w[rnog202f] + w[rnog202b] + w[rnog204f] - w[rnog204b]
		+ w[rnog205f] - w[rnog205b] - w[rnog209f] + w[rnog209b]
		- w[rnog212f] + w[rnog212b] - w[rnog214f] + w[rnog214b]
		+ w[rnog218f] - w[rnog218b] + w[rnog221f] - w[rnog221b]
		- w[rnog223f] + w[rnog223b] + w[rnog224f] - w[rnog224b]
		+ w[rnog230f] - w[rnog230b] + w[rnog231f] - w[rnog231b]
		+ w[rnog234f] - w[rnog234b] + w[rnog235] - w[rnog237f]
		+ w[rnog237b] + w[rnog247f] - w[rnog247b] + w[rnog249f]
		- w[rnog249b] + w[rnog251f] - w[rnog251b] + w[rnog252f]
		- w[rnog252b] + w[rnog254f] - w[rnog254b] + w[rnog257f]
		- w[rnog257b] + w[rnog260f] - w[rnog260b] - w[rnog261f]
		+ w[rnog261b] - w[rnog265f] + w[rnog265b] - w[rnog266f]
		+ w[rnog266b] - w[rnog270f] + w[rnog270f] - w[rnog270b]
		+ w[rnog270b] - w[rnog271f] + w[rnog271b] - w[rnog272f]
		+ w[rnog272b] - w[rnog273] + w[rnog273] + w[rnog275f]
		- w[rnog275b] - w[rnog277f] + w[rnog277b];

	cdot[sC2H3] = w[r61f] - w[r61b] + w[r62f] - w[r62b]
		- w[r112] + w[ri60f] - w[ri60b] + w[ri568]
		+ w[r63] - w[r85] + w[r135] + w[r218]
		+ w[r1335] + w[ri147] - w[r1021170] - w[r1021180f]
		+ w[r1021180b] + w[r1021210] + w[rhai348f] - w[rhai348b];

	cdot[sH2] = w[r61f] - w[r61b] + w[ri17f] - w[ri17b]
		+ w[ri29f] - w[ri29b] + w[ri33f] - w[ri33b]
		+ w[ri48] + w[ri177] + w[r331] + w[r2080]
		+ w[r2081] + w[r2082] + w[r2083] + w[ri2f]
		- w[ri2b] - w[ri9f] + w[ri9b] - w[ri11f]
		+ w[ri11b] + w[ri35f] - w[ri35b] + w[ri74f]
		- w[ri74b] + w[ri95] + w[ri101] + w[ri181]
		+ w[ri211] + w[ri310] + w[ri340f] - w[ri340b]
		+ w[ri343f] - w[ri343b] + w[ri347] + w[ri2731f]
		- w[ri2731b] + w[ri3142] + w[ri3143] + w[ri3144]
		+ w[ri3145] + w[r1021106] + w[r1021128] + w[r1021135]
		+ w[r1021168] + w[r1021178] + w[reth146f] - w[reth146b]
		+ w[reth147] + w[reth148f] - w[reth148b] + w[reth179]
		+ w[reth180f] - w[reth180b] + w[reth235] - w[reth260f]
		+ w[reth260b] + w[rhai350f] - w[rhai350b] + w[rnog8f]
		- w[rnog8b] + w[rnog51f] - w[rnog51b] - w[rnog146f]
		+ w[rnog146b] + w[rnog49f] - w[rnog49b] + w[rnog191f]
		- w[rnog191b] + w[rnog197f] - w[rnog197b] + w[rnog202f]
		- w[rnog202b] + w[rnog209f] - w[rnog209b] + w[rnog214f]
		- w[rnog214b] - w[rnog221f] + w[rnog221b] + w[rnog266f]
		- w[rnog266b] + w[rnog276f] - w[rnog276b] + w[rnog277f]
		- w[rnog277b];

	cdot[sCH2O] = w[r112] - w[ri31f] + w[ri31b] - w[ri33f]
		+ w[ri33b] + w[ri36f] - w[ri36b] + w[ri40]
		+ w[ri43f] - w[ri43b] + w[ri44] - w[ri54]
		+ w[ri83] + w[ri84f] - w[ri84b] + w[ri568]
		+ w[r91] - w[r104f] + w[r104b] + w[r174]
		+ w[r538] + w[r1019] + w[r1335] - w[ri32f]
		+ w[ri32b] - w[ri34f] + w[ri34b] + w[ri35f]
		- w[ri35b] - w[ri38] + w[ri91] - w[ri525]
		+ w[ri797] + w[ri1193] + w[ri1412] + w[ri1458]
		+ w[ri1460] + w[ri3261] + w[ri3377] + w[ri3378]
		+ w[r1021116] + w[r1021118] - w[r1021150] + w[r1021218]
		+ w[reth159] + 2 * w[reth201f] - 2 * w[reth201b] + w[reth236f]
		- w[reth236b] + w[reth238] + w[reth240] + w[rnog94f]
		- w[rnog94b] + w[rnog153f] - w[rnog153b];

	cdot[sHCO] = w[r112] - w[ri29f] + w[ri29b] + w[ri31f]
		- w[ri31b] + w[ri33f] - w[ri33b] - w[ri43f]
		+ w[ri43b] - w[ri46] + w[ri54] + w[ri72]
		+ w[ri1455] + 2 * w[r65] + w[r91] + w[r137]
		- w[ri6f] + w[ri6b] - w[ri12f] + w[ri12b]
		+ w[ri14f] - w[ri14b] - w[ri30f] + w[ri30b]
		+ w[ri32f] - w[ri32b] + w[ri34f] - w[ri34b]
		+ w[ri38] + w[ri77f] - w[ri77b] + 2 * w[ri96]
		+ w[ri153] + w[ri214] + w[ri525] + w[ri1462]
		+ w[r1021136] + w[r1021146] + w[r1021150] + w[r1021192]
		+ w[r1021244] + w[r1021247] + w[reth234] + w[reth236f]
		- w[reth236b] + w[reth240] + 2 * w[rhai195f] - 2 * w[rhai195b]
		+ w[rhai356f] - w[rhai356b] + w[rnog9f] - w[rnog9b]
		+ w[rnog248f] - w[rnog248b] + w[rnog259f] - w[rnog259b]
		+ w[rnog260f] - w[rnog260b];

	cdot[sC2H6] = - w[r113f] + w[r113b] - w[ri5f] + w[ri5b]
		- w[ri17f] + w[ri17b] - w[ri20f] + w[ri20b]
		- w[ri21f] + w[ri21b] + w[ri24f] - w[ri24b]
		- w[rnog154f] + w[rnog154b];

	cdot[sC2H5] = w[r113f] - w[r113b] + w[ri5f] - w[ri5b]
		+ w[ri17f] - w[ri17b] - w[ri111f] + w[ri111b]
		+ w[r174] + w[r255] + w[r330] + w[r2078]
		+ w[r2151] + w[ri15f] - w[ri15b] - w[ri19f]
		+ w[ri19b] + w[ri20f] - w[ri20b] + w[ri21f]
		- w[ri21b] + w[ri153] + w[ri317] + w[reth140f]
		- w[reth140b] + w[reth161f] - w[reth161b] + w[reth239]
		+ w[rnog154f] - w[rnog154b];

	cdot[sCH3O2] = - w[r533] - w[r535] - w[ri521f] + w[ri521b]
		- w[r536] - 2 * w[r538] - w[r2116] - w[r2117]
		- w[r2118] - w[r2119] - w[ri215] - w[ri525]
		- 2 * w[ri540] - w[ri3199] - w[ri3200] - w[ri3201]
		- w[ri3202];

	cdot[sCH3] = - w[r533] - w[ri5f] + w[ri5b] - w[ri22f]
		+ w[ri22b] - w[ri36f] + w[ri36b] - w[ri37f]
		+ w[ri37b] + 2 * w[ri111f] - 2 * w[ri111b] + w[ri136]
		+ w[ri146] + w[ri521f] - w[ri521b] - w[r16f]
		+ w[r16b] - w[r135] + w[r137] + w[r169]
		- w[r1330] + w[r1452] - w[rnog52f] + w[rnog52b]
		+ w[ri2f] - w[ri2b] + w[ri3f] - w[ri3b]
		+ w[ri4f] - w[ri4b] + w[ri14f] - w[ri14b]
		- 2 * w[ri24f] + 2 * w[ri24b] - w[ri35f] + w[ri35b]
		- w[ri38] - w[ri45] + w[ri53f] - w[ri53b]
		+ w[ri89] + w[ri147] + w[ri179] + w[ri206]
		+ 2 * w[ri210] + w[ri262] + w[ri308] + w[ri314]
		- w[ri343f] + w[ri343b] - w[ri344f] + w[ri344b]
		+ w[ri493] + w[ri832] - w[ri1458] + w[ri3138]
		+ w[ri3178] + w[ri3180] + w[r1021103] + w[r1021107]
		- w[r1021110] + w[r1021124f] - w[r1021124b] - w[r1021138]
		+ w[reth139f] - w[reth139b] - w[reth152f] + w[reth152b]
		- w[reth153] - w[reth154f] + w[reth154b] + w[reth159]
		- w[reth167] - w[reth181] - w[reth182] + w[reth234]
		- w[reth239] + w[reth242f] - w[reth242b] + w[rhai349f]
		- w[rhai349b] + w[rhai356f] - w[rhai356b] - w[rnog97f]
		+ w[rnog97b] + w[rnog146f] - w[rnog146b] - w[rnog149f]
		+ w[rnog149b] + 2 * w[rnog150f] - 2 * w[rnog150b] + w[rnog154f]
		- w[rnog154b] - w[rnog124f] + w[rnog124b] - w[rnog211f]
		+ w[rnog211b] - w[rnog255f] + w[rnog255b] - w[rnog256f]
		+ w[rnog256b] - w[rnog275f] + w[rnog275b] - w[rnog276f]
		+ w[rnog276b] - w[rnog280f] + w[rnog280b];

	cdot[sCH3O] = 2 * w[r533] + w[ri22f] - w[ri22b] + w[ri37f]
		- w[ri37b] - w[ri40] - w[ri44] + w[ri522]
		+ 2 * w[r104f] - 2 * w[r104b] + w[r536] + 2 * w[ri540]
		- w[rnog67f] + w[rnog67b] + w[rnog280f] - w[rnog280b];

	cdot[sCH3O2H] = w[r535] - w[ri522] + w[r2116] + w[r2117]
		+ w[r2118] + w[r2119] + w[ri215] + w[ri525]
		+ w[ri3199] + w[ri3200] + w[ri3201] + w[ri3202];

	cdot[sCH4] = w[ri5f] - w[ri5b] + w[r135] + w[r1330]
		+ w[rnog52f] - w[rnog52b] - w[ri2f] + w[ri2b]
		- w[ri3f] + w[ri3b] - w[ri4f] + w[ri4b]
		+ w[ri38] + w[ri45] - w[ri53f] + w[ri53b]
		+ w[r1021110] + w[r1021138] + w[reth152f] - w[reth152b]
		+ w[reth153] + w[reth154f] - w[reth154b] + w[reth181]
		+ w[reth182] - w[rnog150f] + w[rnog150b] + w[rnog211f]
		- w[rnog211b];

	cdot[sCO] = - w[ri7f] + w[ri7b] - w[ri23f] + w[ri23b]
		- w[ri27f] + w[ri27b] + w[ri29f] - w[ri29b]
		+ w[ri46] + w[ri69] + w[ri136] + w[ri1457]
		+ w[ri6f] - w[ri6b] + w[ri12f] - w[ri12b]
		- w[ri28f] + w[ri28b] + w[ri30f] - w[ri30b]
		+ w[ri66] + w[ri89] + w[ri92] + 2 * w[ri98]
		+ w[ri317] + w[ri1516] - w[r1021140f] + w[r1021140b]
		+ w[r1021156] + w[r1021158] + w[r1021163] + 2 * w[r1021164]
		+ w[r1021165] + w[r1021186] + w[r1021193] + w[r1021206]
		+ w[r1021208] + w[r1021209] + 2 * w[r1021244] + w[r1021245f]
		- w[r1021245b] + w[r1021246] - w[reth161f] + w[reth161b]
		+ w[reth238] + w[reth239] + w[reth242f] - w[reth242b]
		+ w[reth261] + 2 * w[reth263f] - 2 * w[reth263b] + 2 * w[reth271]
		+ w[rhai356f] - w[rhai356b] + w[rnog8f] - w[rnog8b]
		+ w[rnog79f] - w[rnog79b] + w[rnog144f] - w[rnog144b]
		+ w[rnog145f] - w[rnog145b] - w[rnog151f] + w[rnog151f]
		- w[rnog151b] + w[rnog151b] + w[rnog153f] - w[rnog153b]
		+ w[rnog90f] - w[rnog90b] + w[rnog122f] - w[rnog122b]
		+ w[rnog217f] - w[rnog217b] + w[rnog222f] - w[rnog222b]
		+ w[rnog223f] - w[rnog223b] + w[rnog224f] - w[rnog224b]
		+ w[rnog225f] - w[rnog225b] + w[rnog227f] - w[rnog227b]
		+ w[rnog228f] - w[rnog228b] + w[rnog232f] - w[rnog232b]
		+ w[rnog236f] - w[rnog236b] + w[rnog245f] - w[rnog245b]
		+ w[rnog257f] - w[rnog257b] + w[rnog263f] - w[rnog263b]
		+ w[rnog265f] - w[rnog265b] + w[rnog269f] - w[rnog269b]
		+ w[rnog272f] - w[rnog272b] + w[rnog274f] - w[rnog274b];

	cdot[sCO2] = w[ri7f] - w[ri7b] + w[ri23f] - w[ri23b]
		+ w[ri27f] - w[ri27b] + w[ri1455] + w[r1452]
		+ w[ri28f] - w[ri28b] + w[ri76] + w[ri90]
		+ w[ri101] + w[r1021160] + w[r1021161] + w[r1021195]
		+ w[reth161f] - w[reth161b] + w[reth267] + w[rhai194f]
		- w[rhai194b] + w[rhai349f] - w[rhai349b] - w[rnog152f]
		+ w[rnog152f] - w[rnog152b] + w[rnog152b] - w[rnog153f]
		+ w[rnog153b] + w[rnog226f] - w[rnog226b] + w[rnog229f]
		- w[rnog229b] + w[rnog262f] - w[rnog262b] + w[rnog268f]
		- w[rnog268b];

	cdot[sO] = w[ri8f] - w[ri8b] - w[ri10f] + w[ri10b]
		- w[ri27f] + w[ri27b] - w[ri36f] + w[ri36b]
		+ w[ri37f] - w[ri37b] - w[ri42f] + w[ri42b]
		- w[ri69] - w[ri87] - w[r63] + w[r1190f]
		- w[r1190b] - w[ri4f] + w[ri4b] - w[ri9f]
		+ w[ri9b] - w[ri14f] + w[ri14b] - w[ri21f]
		+ w[ri21b] + w[ri28f] - w[ri28b] - w[ri30f]
		+ w[ri30b] - w[ri34f] + w[ri34b] + w[ri77f]
		- w[ri77b] - w[ri90] + w[ri91] - w[ri98]
		- w[ri153] - w[ri174] - w[ri176] - w[ri210]
		- w[ri214] + w[ri1347] - w[r1021109] - w[r1021117]
		- w[r1021118] - w[r1021137] - w[r1021146] - w[r1021147]
		+ w[r1021154f] - w[r1021154b] - w[r1021159] - w[r1021160]
		- w[r1021164] - w[r1021169] - w[r1021179] - w[r1021186]
		- w[r1021187] - w[r1021194] - w[reth149f] + w[reth149b]
		- w[reth150] - w[reth151f] + w[reth151b] - w[reth177]
		- w[reth236f] + w[reth236b] - w[reth261] - w[reth267]
		- w[rnog8f] + w[rnog8b] - w[rnog9f] + w[rnog9b]
		+ w[rnog122f] - w[rnog122b] + w[rnog178f] - w[rnog178b]
		+ w[rnog179f] - w[rnog179b] - w[rnog181f] + w[rnog181b]
		- w[rnog182f] + w[rnog182b] + w[rnog185f] - w[rnog185b]
		- w[rnog187f] + w[rnog187b] - w[rnog188f] + w[rnog188b]
		- w[rnog190f] + w[rnog190b] + w[rnog194f] - w[rnog194b]
		- w[rnog200f] + w[rnog200b] - w[rnog201f] + w[rnog201b]
		- w[rnog207f] + w[rnog207b] - w[rnog208f] + w[rnog208b]
		- w[rnog213f] + w[rnog213b] - w[rnog217f] + w[rnog217b]
		+ w[rnog220f] - w[rnog220b] - w[rnog222f] + w[rnog222b]
		- w[rnog231f] + w[rnog231b] - w[rnog232f] + w[rnog232b]
		- w[rnog233f] + w[rnog233b] + w[rnog244f] - w[rnog244b]
		+ w[rnog246f] - w[rnog246b] - w[rnog257f] + w[rnog257b]
		- w[rnog258f] + w[rnog258b] + w[rnog259f] - w[rnog259b]
		- w[rnog262f] + w[rnog262b] - w[rnog263f] + w[rnog263b]
		- w[rnog264f] + w[rnog264b] - w[rnog279f] + w[rnog279b];

	cdot[sH2O2] = - w[ri13P1f] + w[ri13P1b] - w[ri52] + w[ri54]
		+ w[r142] + w[r2092] + w[r2093] + w[r2094]
		+ w[r2095] - w[ri50P1f] + w[ri50P1b] + w[ri51f]
		- w[ri51b] + w[ri53f] - w[ri53b] + w[ri217]
		- w[ri340f] + w[ri340b] - w[ri50P2f] + w[ri50P2b]
		- w[ri13P2f] + w[ri13P2b] + w[ri1511] + w[ri3158]
		+ w[ri3159] + w[ri3160] + w[ri3161] + w[r1021105]
		+ w[reth155] + w[reth156f] - w[reth156b] + w[reth157]
		+ w[reth183] + w[reth184];

	cdot[sC2H2] = - w[ri60f] + w[ri60b] - w[ri69] - w[ri87]
		+ w[ri146] - w[r65] + w[r85] - w[ri88]
		+ w[r1021164] + 2 * w[r1021208] + w[r1021210] + w[r1021245f]
		- w[r1021245b] + w[r1021246] + w[r1021247] + w[reth260f]
		- w[reth260b] + w[reth271] - w[rhai193f] + w[rhai193b]
		+ w[rnog124f] - w[rnog124b];

	cdot[sCH2] = w[ri69] - w[ri72] - w[ri66] - w[ri74f]
		+ w[ri74b] - w[ri75f] + w[ri75b] - w[ri76]
		+ w[ri90] - w[ri91] + w[ri92] - w[ri101]
		+ w[ri343f] - w[ri343b] + w[ri344f] - w[ri344b]
		+ w[rhai194f] - w[rhai194b] + w[rnog142f] - w[rnog142b]
		+ w[rnog148f] - w[rnog148b] + w[rnog151f] - w[rnog151b]
		+ w[rnog152f] - w[rnog152b] - w[rnog123f] + w[rnog123b]
		+ w[rnog238f] - w[rnog238b] - w[rnog242f] + w[rnog242b]
		- w[rnog249f] + w[rnog249b] - w[rnog250f] + w[rnog250b]
		- w[rnog251f] + w[rnog251b] + w[rnog261f] - w[rnog261b];

	cdot[sCH2OH] = - w[ri83] - w[ri84f] + w[ri84b] + w[ri1457]
		+ w[reth139f] - w[reth139b] - w[rnog62f] + w[rnog62b];

	cdot[sHCCO] = w[ri87] + w[ri94] - w[ri1455] + w[ri95]
		- w[ri96] - w[ri98] + w[r1021194] + w[reth262f]
		- w[reth262b] - w[reth267] - 2 * w[reth271] - w[rnog79f]
		+ w[rnog79b] - w[rnog274f] + w[rnog274b];

	cdot[sCH2CO] = - w[ri94] - w[ri1457] - w[r91] + w[r1400f]
		- w[r1400b] + w[r1401f] - w[r1401b] - w[r1452]
		+ w[ri88] - w[ri89] - w[ri90] - w[ri92]
		- w[ri95] + w[ri210] + w[ri314] + w[ri347]
		+ w[r1021164] + w[reth235] + w[reth237] + w[reth243];

	cdot[sCH3CO] = - w[ri136] + w[r142] + w[r1019] + w[ri308]
		- w[ri347] + w[ri1460] + w[reth174] + w[reth177]
		+ w[reth179] + w[reth181] + w[reth183] + w[reth185f]
		- w[reth185b];

	cdot[sAXC3H5] = - w[ri146] + w[ri159] + w[ri177] - w[ri570]
		+ w[r330] + w[r1330] - w[r1335] + w[r2185]
		+ w[ri151f] - w[ri151b] + w[ri174] + w[r1021193]
		+ w[r1021195];

	cdot[sC3H6] = - w[ri159] - w[ri177] + w[r170f] - w[r170b]
		+ w[r171] - w[r174] + w[r218] - w[r1330]
		+ w[r2144] + w[r2147] - w[ri147] - w[ri151f]
		+ w[ri151b] - w[ri153] - w[ri173] - w[ri174]
		- w[ri176] - w[ri179] + w[ri180f] - w[ri180b]
		+ w[ri181] + w[ri182] + w[ri262] + w[ri797]
		+ w[ri2743] + w[ri3189] + w[reth167];

	cdot[sC3H5O] = - w[ri568] + w[ri570];

	cdot[sCH3OH] = w[r16f] - w[r16b] - w[r104f] + w[r104b]
		+ w[r538] + w[rnog147f] - w[rnog147b];

	cdot[sCH3CHO] = - w[r137] - w[r142] + w[r2281] + w[r2327];

	cdot[sNXC3H7] = - w[r169] - w[r170f] + w[r170b] - w[r171]
		+ w[r464] + w[r2079];

	cdot[sC5H9] = - w[r218] + w[r331] + w[r335];

	cdot[sPXC4H9] = - w[r255] + w[r2079] + w[r2144] + w[r2147]
		+ w[r2185];

	cdot[s1XC5H10] = - w[r330] - w[r331] - w[r335] + w[r465f]
		- w[r465b] + w[r2151] + w[r2281];

	cdot[s1XC5H11] = - w[r464] - w[r465f] + w[r465b] + w[r2078]
		+ w[r2142];

	cdot[sCH3COCH2O2] = - w[r1014f] + w[r1014b] - w[r1017];

	cdot[sCH3COCH2] = w[r1014f] - w[r1014b] + w[r2327] + w[ri309]
		+ w[ri310] - w[ri314] + w[ri1347];

	cdot[sCH3COCH2O2H] = w[r1017] - w[r1018];

	cdot[sCH3COCH2O] = w[r1018] - w[r1019];

	cdot[sCH2CHO] = - w[r1400f] + w[r1400b] - w[r1401f] + w[r1401b]
		+ w[r1021194];

	cdot[sNXC7H16] = - w[r2078] - w[r2079] - w[r2080] - w[r2081]
		- w[r2082] - w[r2083] - w[r2088] - w[r2089]
		- w[r2090] - w[r2091] - w[r2092] - w[r2093]
		- w[r2094] - w[r2095] - w[r2100f] + w[r2100b]
		- w[r2101f] + w[r2101b] - w[r2102f] + w[r2102b]
		- w[r2103f] + w[r2103b] - w[r2116] - w[r2117]
		- w[r2118] - w[r2119];

	cdot[s1XC7H15] = w[r2080] + w[r2088] + w[r2092] + w[r2100f]
		- w[r2100b] + w[r2116] - w[r2142] - w[r2159f]
		+ w[r2159b] - w[r2160f] + w[r2160b] - w[r2162f]
		+ w[r2162b] + w[r2188f] - w[r2188b];

	cdot[s2XC7H15] = w[r2081] + w[r2089] + w[r2093] + w[r2101f]
		- w[r2101b] + w[r2117] - w[r2144] - w[r2161f]
		+ w[r2161b] + w[r2162f] - w[r2162b] + w[r2189f]
		- w[r2189b];

	cdot[s3XC7H15] = w[r2082] + w[r2090] + w[r2094] + w[r2102f]
		- w[r2102b] + w[r2118] - w[r2147] + w[r2159f]
		- w[r2159b] + w[r2161f] - w[r2161b] + w[r2190f]
		- w[r2190b];

	cdot[s4XC7H15] = w[r2083] + w[r2091] + w[r2095] + w[r2103f]
		- w[r2103b] + w[r2119] - w[r2151] + w[r2160f]
		- w[r2160b] + w[r2191f] - w[r2191b];

	cdot[s1XC7H14] = - w[r2185] + w[r2258f] - w[r2258b];

	cdot[s2XC7H15O2] = - w[r2188f] + w[r2188b] - w[r2189f] + w[r2189b]
		- w[r2190f] + w[r2190b] - w[r2191f] + w[r2191b]
		- w[r2212P1f] + w[r2212P1b] - w[r2212P2f] + w[r2212P2b]
		- w[r2214f] + w[r2214b];

	cdot[s2X3XC7H14O2H] = w[r2212P1f] - w[r2212P1b] + w[r2212P2f] - w[r2212P2b]
		- w[r2258f] + w[r2258b] + w[r2290f] - w[r2290b];

	cdot[s2X4XC7H14O2H] = w[r2214f] - w[r2214b] - w[r2281] + w[r2291f]
		- w[r2291b];

	cdot[s2X4XC7H14OOHO2] = - w[r2290f] + w[r2290b] - w[r2291f] + w[r2291b]
		- w[r2309];

	cdot[sNC7KET24XC7H14O3] = w[r2309] - w[r2327];

	cdot[sCH] = w[ri74f] - w[ri74b] + w[ri75f] - w[ri75b]
		- w[ri77f] + w[ri77b] + w[reth261] + w[reth267]
		+ w[rnog51f] - w[rnog51b] - w[rnog49f] + w[rnog49b]
		- w[rnog240f] + w[rnog240b] - w[rnog241f] + w[rnog241b]
		- w[rnog246f] + w[rnog246b] - w[rnog247f] + w[rnog247b]
		- w[rnog248f] + w[rnog248b];

	cdot[sC3H5XTXC3H5] = w[ri173] + w[ri176] + w[ri206] + w[ri1193]
		- w[ri1347] - w[ri1460] + w[ri1516];

	cdot[sIXC3H7] = - w[ri180f] + w[ri180b] - w[ri181] - w[ri182]
		+ w[ri214] + w[ri2720] + w[ri3141];

	cdot[sIXC4H8] = - w[ri206] - w[ri207f] + w[ri207b] - w[ri210]
		- w[ri211] - w[ri213] - w[ri214] - w[ri215]
		- w[ri217] + w[ri259f] - w[ri259b] + w[ri263]
		+ w[ri493] - w[ri711f] + w[ri711b] + w[ri790]
		+ w[ri791] - w[ri1536f] + w[ri1536b] + w[ri2720]
		+ w[ri3179] + w[ri3184] + w[ri3506];

	cdot[sIXC4H7] = w[ri207f] - w[ri207b] + w[ri211] + w[ri213]
		+ w[ri215] + w[ri217] - w[ri710] + w[ri1536f]
		- w[ri1536b] + w[ri2743];

	cdot[sTXC4H9] = - w[ri259f] + w[ri259b] + w[ri1393f] - w[ri1393b]
		- w[ri1394] + w[ri3140] + w[ri3184] + w[ri3498];

	cdot[sIXC4H9] = - w[ri262] - w[ri263] + w[ri1396f] - w[ri1396b]
		+ w[ri3140] + w[ri3179];

	cdot[sCH3COCH3] = - w[ri308] - w[ri309] - w[ri310] + w[ri832]
		+ w[ri1462] + w[ri3263];

	cdot[sC2H5CO] = - w[ri317] + w[ri1412];

	cdot[sNEOXC5H11] = - w[ri493] + w[ri3141] + w[ri3189] + w[ri3263];

	cdot[sIXC4H7O] = w[ri710] - w[ri1193] + w[ri1468] - w[ri1517]
		- w[ri1546];

	cdot[sIXC4H8O] = w[ri711f] - w[ri711b] + w[ri793] + w[ri794];

	cdot[sIC4H8OOHXIO2XC4H9O4] = - w[ri776f] + w[ri776b] - w[ri779];

	cdot[sIC4H8OOHXIXC4H9O2] = w[ri776f] - w[ri776b] - w[ri797] + w[ri901f]
		- w[ri901b];

	cdot[sIC4H8OOHXTO2XC4H9O4] = - w[ri778f] + w[ri778b] - w[ri780];

	cdot[sIC4H8OOHXTXC4H9O2] = w[ri778f] - w[ri778b] - w[ri791] - w[ri793]
		+ w[ri903f] - w[ri903b] + w[ri3506];

	cdot[sIC4KETIIXC4H8O3] = w[ri779] + w[ri780] - w[ri1412];

	cdot[sTC4H8OOHXIXC4H9O2] = - w[ri790] - w[ri794] + w[ri902f] - w[ri902b];

	cdot[sTXC4H9O] = - w[ri832] + w[ri1394];

	cdot[sIXC4H9O2] = - w[ri901f] + w[ri901b] - w[ri903f] + w[ri903b]
		- w[ri1396f] + w[ri1396b];

	cdot[sTXC4H9O2] = - w[ri902f] + w[ri902b] - w[ri1393f] + w[ri1393b];

	cdot[sTXC3H6CHO] = - w[ri1461] - w[ri1463] + w[ri1466] - w[ri1508]
		+ w[ri3367];

	cdot[sTXC3H6OCHO] = w[ri1461] - w[ri1462];

	cdot[sIXC3H5CHO] = w[ri1463] - w[ri1510] - w[ri1511] + w[ri1517]
		+ w[ri1546];

	cdot[sIXC3H7CHO] = - w[ri1466] + w[ri1508] + w[ri3367];

	cdot[sIXC4H7OOH] = - w[ri1468] + w[ri3498];

	cdot[sIXC3H5CO] = w[ri1510] + w[ri1511] - w[ri1516];

	cdot[sYXC7H15] = - w[ri2720] - w[ri2722f] + w[ri2722b] - w[ri2727f]
		+ w[ri2727b] + w[ri3138] + w[ri3261];

	cdot[sYXC7H14] = w[ri2722f] - w[ri2722b] + w[ri2727f] - w[ri2727b]
		- w[ri2731f] + w[ri2731b] - w[ri2735f] + w[ri2735b]
		+ w[ri3178] + w[ri3180] + w[ri3377] + w[ri3378];

	cdot[sXXC7H13] = w[ri2731f] - w[ri2731b] + w[ri2735f] - w[ri2735b]
		- w[ri2743];

	cdot[sIXC8H18] = - w[ri3134] - w[ri3135] - w[ri3136] - w[ri3137]
		- w[ri3138] - w[ri3140] - w[ri3141] - w[ri3142]
		- w[ri3143] - w[ri3144] - w[ri3145] - w[ri3150]
		- w[ri3151] - w[ri3152] - w[ri3153] - w[ri3158]
		- w[ri3159] - w[ri3160] - w[ri3161] - w[ri3166]
		- w[ri3167] - w[ri3168] - w[ri3169] - w[ri3199]
		- w[ri3200] - w[ri3201] - w[ri3202] - w[ri3203]
		- w[ri3205] - w[ri3207] - w[ri3209] - w[ri3211]
		- w[ri3213] - w[ri3354] - w[ri3356];

	cdot[sAXC8H17] = w[ri3134] + w[ri3142] + w[ri3150] + w[ri3158]
		+ w[ri3166] + w[ri3199] + w[ri3203] + w[ri3205]
		- w[ri3178] - w[ri3179] - w[ri3364f] + w[ri3364b]
		- w[ri3365f] + w[ri3365b] + w[ri3215f] - w[ri3215b];

	cdot[sBXC8H17] = w[ri3135] + w[ri3143] + w[ri3151] + w[ri3159]
		+ w[ri3167] + w[ri3200] + w[ri3207] + w[ri3209]
		- w[ri3180] + w[ri3216f] - w[ri3216b];

	cdot[sCXC8H17] = w[ri3136] + w[ri3144] + w[ri3152] + w[ri3160]
		+ w[ri3168] + w[ri3201] + w[ri3211] + w[ri3213]
		- w[ri3184] + w[ri3365f] - w[ri3365b] + w[ri3217f]
		- w[ri3217b];

	cdot[sDXC8H17] = w[ri3137] + w[ri3145] + w[ri3153] + w[ri3161]
		+ w[ri3169] + w[ri3202] + w[ri3354] + w[ri3356]
		- w[ri3189] + w[ri3364f] - w[ri3364b] + w[ri3218f]
		- w[ri3218b];

	cdot[sAXC8H17O2] = - w[ri3203] - w[ri3207] - w[ri3211] - w[ri3354]
		- w[ri3215f] + w[ri3215b] - w[ri3216f] + w[ri3216b]
		- w[ri3218f] + w[ri3218b] - w[ri3265f] + w[ri3265b]
		- w[ri3266f] + w[ri3266b] - w[ri3267f] + w[ri3267b];

	cdot[sAXC8H17O2H] = w[ri3203] + w[ri3207] + w[ri3211] + w[ri3354]
		- w[ri3257];

	cdot[sCXC8H17O2] = - w[ri3205] - w[ri3209] - w[ri3213] - w[ri3356]
		- w[ri3217f] + w[ri3217b] - w[ri3272P1f] + w[ri3272P1b]
		- w[ri3272P2] - w[ri3272P3];

	cdot[sCXC8H17O2H] = w[ri3205] + w[ri3209] + w[ri3213] + w[ri3356]
		- w[ri3259];

	cdot[sAXC8H17O] = w[ri3257] - w[ri3261];

	cdot[sCXC8H17O] = w[ri3259] - w[ri3263];

	cdot[sAC8H16OOHXAXC8H17O2] = w[ri3265f] - w[ri3265b] + w[ri3297f] - w[ri3297b]
		- w[ri3377] - w[ri3498];

	cdot[sAC8H16OOHXBXC8H17O2] = w[ri3266f] - w[ri3266b] + w[ri3298f] - w[ri3298b]
		- w[ri3378];

	cdot[sAC8H16OOHXCXC8H17O2] = w[ri3267f] - w[ri3267b] + w[ri3299f] - w[ri3299b]
		- w[ri3506];

	cdot[sCC8H16OOHXAXC8H17O2] = w[ri3272P1f] - w[ri3272P1b] + w[ri3272P2] + w[ri3272P3]
		+ w[ri3306f] - w[ri3306b];

	cdot[sAC8H16OOHXAO2XC8H17O4] = - w[ri3297f] + w[ri3297b] - w[ri3311];

	cdot[sAC8H16OOHXBO2XC8H17O4] = - w[ri3298f] + w[ri3298b] - w[ri3312];

	cdot[sAC8H16OOHXCO2XC8H17O4] = - w[ri3299f] + w[ri3299b] - w[ri3313];

	cdot[sCC8H16OOHXAO2XC8H17O4] = - w[ri3306f] + w[ri3306b];

	cdot[sIC8KETABXC8H16O3] = w[ri3311] + w[ri3312] + w[ri3313] - w[ri3367];

	cdot[sC6H5CH2] = w[r1021102f] - w[r1021102b] + w[r1021104f] - w[r1021104b]
		+ w[r1021105] + w[r1021106] + w[r1021108] + w[r1021109]
		+ w[r1021110] + w[r1021111] - w[r1021112f] + w[r1021112b]
		- w[r1021113] + w[r1021114f] - w[r1021114b] - w[r1021115]
		- w[r1021116] - w[r1021117] - w[r1021118] - w[r1021119]
		+ w[r1021124f] - w[r1021124b] - w[r1021125f] + w[r1021125b]
		- w[r1021185f] + w[r1021185b];

	cdot[sC6H5] = w[r1021103] - w[r1021111] + w[r1021118] - w[r1021140f]
		+ w[r1021140b] - w[r1021142] + w[r1021143f] - w[r1021143b]
		+ w[r1021144] - w[r1021150] - w[r1021153] - w[r1021154f]
		+ w[r1021154b] - w[r1021155] - w[r1021172];

	cdot[sC6H6] = w[r1021107] + w[r1021111] + w[r1021136] + w[r1021142]
		- w[r1021143f] + w[r1021143b] - w[r1021144] - w[r1021146]
		- w[r1021147] + w[r1021150] - w[r1021167f] + w[r1021167b]
		+ w[r1021172];

	cdot[sC6H5OH] = - w[r1021112f] + w[r1021112b] + w[r1021147] - w[r1021157f]
		+ w[r1021157b] + w[r1021162] - w[r1021166] + w[r1021167f]
		- w[r1021167b] - w[r1021168] - w[r1021169] - w[r1021170]
		- w[r1021171f] + w[r1021171b] - w[r1021172] - w[r1021173f]
		+ w[r1021173b] + w[r1021182];

	cdot[sC6H5O] = w[r1021112f] - w[r1021112b] + w[r1021116] + w[r1021153]
		+ w[r1021154f] - w[r1021154b] - w[r1021156] + w[r1021157f]
		- w[r1021157b] - w[r1021158] - w[r1021159] - w[r1021160]
		- w[r1021161] - w[r1021162] + w[r1021166] + w[r1021168]
		+ w[r1021169] + w[r1021170] + w[r1021171f] - w[r1021171b]
		+ w[r1021172] + w[r1021173f] - w[r1021173b] - w[r1021182];

	cdot[sC6H5CHO] = - w[r1021113] + w[r1021115] + w[r1021117] + w[r1021121]
		- w[r1021133] - w[r1021134] - w[r1021135] - w[r1021136]
		- w[r1021137] - w[r1021138];

	cdot[sC6H5CO] = w[r1021113] + w[r1021133] + w[r1021134] + w[r1021135]
		+ w[r1021137] + w[r1021138] + w[r1021140f] - w[r1021140b];

	cdot[sC4H5XIXC4H5] = - w[r1021114f] + w[r1021114b] - w[r1021173f] + w[r1021173b]
		- w[r1021211] - w[r1021212f] + w[r1021212b];

	cdot[sC4H6] = w[r1021114f] - w[r1021114b] + w[r1021171f] - w[r1021171b]
		+ w[r1021173f] - w[r1021173b] + w[r1021211];

	cdot[sC6H5CH2O] = w[r1021119] - w[r1021121] + w[r1021126f] - w[r1021126b]
		+ w[r1021127] + w[r1021128];

	cdot[sC6H5C2H5] = - w[r1021124f] + w[r1021124b];

	cdot[sC6H5CH2OH] = w[r1021125f] - w[r1021125b] - w[r1021126f] + w[r1021126b]
		- w[r1021127] - w[r1021128];

	cdot[sC5H5] = w[r1021146] + w[r1021156] + w[r1021160] + w[r1021161]
		- w[r1021174f] + w[r1021174b] + w[r1021175f] - w[r1021175b]
		+ w[r1021177] + w[r1021178] + w[r1021179] + w[r1021180f]
		- w[r1021180b] + w[r1021182] + w[r1021185f] - w[r1021185b]
		- w[r1021186] - w[r1021187] - w[r1021191] - w[r1021192]
		+ w[r1021211];

	cdot[sC6H4O2] = w[r1021155] + w[r1021159] - w[r1021163] - w[r1021164]
		- w[r1021165];

	cdot[sC5H6] = w[r1021158] + w[r1021174f] - w[r1021174b] - w[r1021175f]
		+ w[r1021175b] - w[r1021177] - w[r1021178] - w[r1021179]
		- w[r1021180f] + w[r1021180b] - w[r1021182] - w[r1021185f]
		+ w[r1021185b] - w[r1021211];

	cdot[sC5H4O] = w[r1021163] + w[r1021187] + w[r1021205] - w[r1021208]
		- w[r1021209];

	cdot[sC5H5O] = w[r1021165] + w[r1021191] - w[r1021205] - w[r1021206];

	cdot[sC4H5XNXC4H5] = - w[r1021171f] + w[r1021171b] + w[r1021186] + w[r1021206]
		+ w[r1021209] - w[r1021210] + w[r1021212f] - w[r1021212b]
		- w[r1021218];

	cdot[sCH2CHCHCO] = w[r1021192] - w[r1021193] - w[r1021194] - w[r1021195];

	cdot[sCHCHCHO] = w[r1021218] - w[r1021246] - w[r1021247];

	cdot[sHCCCO] = - w[r1021244] - w[r1021245f] + w[r1021245b];

	cdot[sC2H5OH] = - w[reth139f] + w[reth139b] - w[reth140f] + w[reth140b]
		- w[reth141] - w[reth143f] + w[reth143b] - w[reth144]
		- w[reth145] - w[reth146f] + w[reth146b] - w[reth147]
		- w[reth148f] + w[reth148b] - w[reth149f] + w[reth149b]
		- w[reth150] - w[reth151f] + w[reth151b] - w[reth152f]
		+ w[reth152b] - w[reth153] - w[reth154f] + w[reth154b]
		- w[reth155] - w[reth156f] + w[reth156b] - w[reth157];

	cdot[sC2H4OH] = w[reth143f] - w[reth143b] + w[reth146f] - w[reth146b]
		+ w[reth149f] - w[reth149b] + w[reth152f] - w[reth152b]
		+ w[reth156f] - w[reth156b] - w[reth200f] + w[reth200b];

	cdot[sCH3CHOH] = w[reth144] + w[reth147] + w[reth150] + w[reth153]
		+ w[reth155] - w[reth165] - w[reth167] - w[reth169]
		- w[reth171] - w[reth172f] + w[reth172b] - w[reth173];

	cdot[sCH3CH2O] = w[reth145] + w[reth148f] - w[reth148b] + w[reth151f]
		- w[reth151b] + w[reth154f] - w[reth154b] + w[reth157]
		- w[reth158] - w[reth159] - w[reth160f] + w[reth160b]
		- w[reth161f] + w[reth161b] - w[reth163] - w[reth164];

	cdot[sCH3HCO] = w[reth158] + w[reth160f] - w[reth160b] + w[reth164]
		+ w[reth165] + w[reth171] + w[reth172f] - w[reth172b]
		+ w[reth173] - w[reth174] - w[reth175f] + w[reth175b]
		- w[reth177] - w[reth179] - w[reth180f] + w[reth180b]
		- w[reth181] - w[reth182] - w[reth183] - w[reth184]
		- w[reth185f] + w[reth185b] + w[reth241f] - w[reth241b];

	cdot[sCH2HCO] = w[reth175f] - w[reth175b] + w[reth180f] - w[reth180b]
		+ w[reth182] + w[reth184] - w[reth234] - w[reth235]
		- w[reth236f] + w[reth236b] - w[reth237] - w[reth238]
		- w[reth239] - w[reth240] - w[reth241f] + w[reth241b]
		- w[reth242f] + w[reth242b] - w[reth243];

	cdot[sHOC2H4O2] = w[reth200f] - w[reth200b] - w[reth201f] + w[reth201b];

	cdot[sC2H] = - w[reth260f] + w[reth260b] - w[reth261] - w[reth262f]
		+ w[reth262b] - w[reth263f] + w[reth263b] + w[rnog123f]
		- w[rnog123b];

	cdot[sH2CC] = w[rhai193f] - w[rhai193b] - w[rhai194f] + w[rhai194b]
		- w[rhai195f] + w[rhai195b] + w[rhai350f] - w[rhai350b];

	cdot[sC3H3] = - w[rhai352f] + w[rhai352b] - w[rhai353f] + w[rhai353b]
		- w[rhai354f] + w[rhai354b] - w[rhai355f] + w[rhai355b];

	cdot[sPXC3H4] = w[rhai352f] - w[rhai352b] + w[rhai354f] - w[rhai354b]
		- w[rhai356f] + w[rhai356b] - w[rhai357f] + w[rhai357b];

	cdot[sAXC3H4] = w[rhai353f] - w[rhai353b] + w[rhai355f] - w[rhai355b]
		+ w[rhai357f] - w[rhai357b];

	cdot[s1XCH2] = - w[rnog8f] + w[rnog8b] - w[rnog9f] + w[rnog9b]
		- w[rnog51f] + w[rnog51b] + w[rnog62f] - w[rnog62b]
		+ w[rnog67f] - w[rnog67b] + w[rnog79f] - w[rnog79b]
		- w[rnog94f] + w[rnog94b] + w[rnog97f] - w[rnog97b]
		- w[rnog142f] + w[rnog142b] - w[rnog144f] + w[rnog144b]
		- w[rnog145f] + w[rnog145b] - w[rnog146f] + w[rnog146b]
		- w[rnog147f] + w[rnog147b] - w[rnog148f] + w[rnog148b]
		- w[rnog149f] + w[rnog149b] - w[rnog150f] + w[rnog150b]
		- w[rnog151f] + w[rnog151b] - w[rnog152f] + w[rnog152b]
		- w[rnog153f] + w[rnog153b] - w[rnog154f] + w[rnog154b]
		- w[rnog243f] + w[rnog243b] - w[rnog252f] + w[rnog252b]
		- w[rnog253f] + w[rnog253b] - w[rnog254f] + w[rnog254b];

	cdot[sC] = w[rnog49f] - w[rnog49b] - w[rnog90f] + w[rnog90b]
		- w[rnog122f] + w[rnog122b] - w[rnog123f] + w[rnog123b]
		- w[rnog124f] + w[rnog124b] - w[rnog239f] + w[rnog239b]
		- w[rnog244f] + w[rnog244b] - w[rnog245f] + w[rnog245b];

	cdot[sN] = - w[rnog178f] + w[rnog178b] - w[rnog179f] + w[rnog179b]
		- w[rnog180f] + w[rnog180b] + w[rnog191f] - w[rnog191b]
		+ w[rnog193f] - w[rnog193b] - w[rnog196f] + w[rnog196b]
		+ w[rnog217f] - w[rnog217b] - w[rnog225f] + w[rnog225b]
		+ w[rnog227f] - w[rnog227b] - w[rnog238f] + w[rnog238b]
		+ w[rnog239f] - w[rnog239b] + w[rnog240f] - w[rnog240b]
		+ w[rnog245f] - w[rnog245b] + w[rnog248f] - w[rnog248b]
		- w[rnog275f] + w[rnog275b] - w[rnog276f] + w[rnog276b];

	cdot[sNO] = - w[rnog178f] + w[rnog178b] + w[rnog179f] - w[rnog179b]
		+ w[rnog180f] - w[rnog180b] + 2 * w[rnog182f] - 2 * w[rnog182b]
		- w[rnog186f] + w[rnog186b] - w[rnog187f] + w[rnog187b]
		+ w[rnog188f] - w[rnog188b] + w[rnog189f] - w[rnog189b]
		+ w[rnog190f] - w[rnog190b] + w[rnog195f] - w[rnog195b]
		- w[rnog198f] + w[rnog198b] - w[rnog199f] + w[rnog199b]
		+ w[rnog208f] - w[rnog208b] - w[rnog212f] + w[rnog212b]
		+ w[rnog213f] - w[rnog213b] + w[rnog214f] - w[rnog214b]
		+ w[rnog215f] - w[rnog215b] + w[rnog216f] - w[rnog216b]
		+ w[rnog222f] - w[rnog222b] + w[rnog224f] - w[rnog224b]
		+ w[rnog226f] - w[rnog226b] - w[rnog228f] + w[rnog228b]
		- w[rnog229f] + w[rnog229b] - w[rnog244f] + w[rnog244b]
		- w[rnog245f] + w[rnog245b] - w[rnog246f] + w[rnog246b]
		- w[rnog247f] + w[rnog247b] - w[rnog248f] + w[rnog248b]
		- w[rnog249f] + w[rnog249b] - w[rnog250f] + w[rnog250b]
		- w[rnog251f] + w[rnog251b] - w[rnog252f] + w[rnog252b]
		- w[rnog253f] + w[rnog253b] - w[rnog254f] + w[rnog254b]
		- w[rnog255f] + w[rnog255b] - w[rnog256f] + w[rnog256b]
		+ w[rnog258f] - w[rnog258b] - w[rnog274f] + w[rnog274b]
		+ w[rnog280f] - w[rnog280b];

	cdot[sN2O] = - w[rnog181f] + w[rnog181b] - w[rnog182f] + w[rnog182b]
		- w[rnog183f] + w[rnog183b] - w[rnog184f] + w[rnog184b]
		- w[rnog185f] + w[rnog185b] + w[rnog199f] - w[rnog199b]
		+ w[rnog228f] - w[rnog228b];

	cdot[sNO2] = w[rnog186f] - w[rnog186b] + w[rnog187f] - w[rnog187b]
		- w[rnog188f] + w[rnog188b] - w[rnog189f] + w[rnog189b]
		- w[rnog280f] + w[rnog280b];

	cdot[sNH] = - w[rnog190f] + w[rnog190b] - w[rnog191f] + w[rnog191b]
		- w[rnog192f] + w[rnog192b] - w[rnog193f] + w[rnog193b]
		- w[rnog194f] + w[rnog194b] - w[rnog195f] + w[rnog195b]
		- w[rnog196f] + w[rnog196b] - w[rnog197f] + w[rnog197b]
		- w[rnog198f] + w[rnog198b] - w[rnog199f] + w[rnog199b]
		+ w[rnog200f] - w[rnog200b] + w[rnog202f] - w[rnog202b]
		+ w[rnog203f] - w[rnog203b] + w[rnog208f] - w[rnog208b]
		+ w[rnog223f] - w[rnog223b] + w[rnog232f] - w[rnog232b]
		+ w[rnog242f] - w[rnog242b] + w[rnog243f] - w[rnog243b]
		+ w[rnog262f] - w[rnog262b] + w[rnog269f] - w[rnog269b];

	cdot[sHNO] = w[rnog192f] - w[rnog192b] + w[rnog194f] - w[rnog194b]
		+ w[rnog197f] - w[rnog197b] + w[rnog201f] - w[rnog201b]
		+ w[rnog212f] - w[rnog212b] - w[rnog213f] + w[rnog213b]
		- w[rnog214f] + w[rnog214b] - w[rnog215f] + w[rnog215b]
		- w[rnog216f] + w[rnog216b] + w[rnog263f] - w[rnog263b];

	cdot[sNH2] = - w[rnog200f] + w[rnog200b] - w[rnog201f] + w[rnog201b]
		- w[rnog202f] + w[rnog202b] - w[rnog203f] + w[rnog203b]
		+ w[rnog236f] - w[rnog236b] + w[rnog265f] - w[rnog265b]
		+ w[rnog268f] - w[rnog268b] + w[rnog272f] - w[rnog272b]
		+ w[rnog277f] - w[rnog277b] + w[rnog278f] - w[rnog278b]
		+ w[rnog279f] - w[rnog279b];

	cdot[sNNH] = - w[rnog204f] + w[rnog204b] - w[rnog205f] + w[rnog205b]
		- w[rnog206f] + w[rnog206b] - w[rnog207f] + w[rnog207b]
		- w[rnog208f] + w[rnog208b] - w[rnog209f] + w[rnog209b]
		- w[rnog210f] + w[rnog210b] - w[rnog211f] + w[rnog211b];

	cdot[sCN] = - w[rnog217f] + w[rnog217b] - w[rnog218f] + w[rnog218b]
		- w[rnog219f] + w[rnog219b] - w[rnog220f] + w[rnog220b]
		- w[rnog221f] + w[rnog221b] + w[rnog230f] - w[rnog230b]
		+ w[rnog233f] - w[rnog233b] + w[rnog239f] - w[rnog239b]
		+ w[rnog244f] - w[rnog244b];

	cdot[sNCO] = w[rnog218f] - w[rnog218b] + w[rnog220f] - w[rnog220b]
		- w[rnog222f] + w[rnog222b] - w[rnog223f] + w[rnog223b]
		- w[rnog224f] + w[rnog224b] - w[rnog225f] + w[rnog225b]
		- w[rnog226f] + w[rnog226b] - w[rnog227f] + w[rnog227b]
		- w[rnog228f] + w[rnog228b] - w[rnog229f] + w[rnog229b]
		+ w[rnog231f] - w[rnog231b] + w[rnog247f] - w[rnog247b]
		+ w[rnog264f] - w[rnog264b] + w[rnog266f] - w[rnog266b]
		+ w[rnog267f] - w[rnog267b];

	cdot[sHCN] = w[rnog219f] - w[rnog219b] + w[rnog221f] - w[rnog221b]
		- w[rnog230f] + w[rnog230b] - w[rnog231f] + w[rnog231b]
		- w[rnog232f] + w[rnog232b] - w[rnog233f] + w[rnog233b]
		- w[rnog234f] + w[rnog234b] - w[rnog235] - w[rnog236f]
		+ w[rnog236b] - w[rnog237f] + w[rnog237b] + w[rnog240f]
		- w[rnog240b] + w[rnog242f] - w[rnog242b] + w[rnog243f]
		- w[rnog243b] + w[rnog246f] - w[rnog246b] + w[rnog250f]
		- w[rnog250b] + w[rnog253f] - w[rnog253b] + w[rnog255f]
		- w[rnog255b] + w[rnog258f] - w[rnog258b] + w[rnog271f]
		- w[rnog271b] + w[rnog276f] - w[rnog276b];

	cdot[sHOCN] = w[rnog234f] - w[rnog234b] - w[rnog273];

	cdot[sHNCO] = w[rnog235] + w[rnog249f] - w[rnog249b] + w[rnog252f]
		- w[rnog252b] - w[rnog262f] + w[rnog262b] - w[rnog263f]
		+ w[rnog263b] - w[rnog264f] + w[rnog264b] - w[rnog265f]
		+ w[rnog265b] - w[rnog266f] + w[rnog266b] - w[rnog267f]
		+ w[rnog267b] - w[rnog268f] + w[rnog268b] - w[rnog269f]
		+ w[rnog269b] + w[rnog270f] - w[rnog270b] + w[rnog273];

	cdot[sH2CN] = w[rnog237f] - w[rnog237b] - w[rnog238f] + w[rnog238b]
		+ w[rnog256f] - w[rnog256b] + w[rnog275f] - w[rnog275b];

	cdot[sHCNN] = w[rnog241f] - w[rnog241b] - w[rnog257f] + w[rnog257b]
		- w[rnog258f] + w[rnog258b] - w[rnog259f] + w[rnog259b]
		- w[rnog260f] + w[rnog260b] - w[rnog261f] + w[rnog261b];

	cdot[sHCNO] = w[rnog251f] - w[rnog251b] + w[rnog254f] - w[rnog254b]
		- w[rnog270f] + w[rnog270b] - w[rnog271f] + w[rnog271b]
		- w[rnog272f] + w[rnog272b] + w[rnog274f] - w[rnog274b];

	cdot[sNH3] = - w[rnog277f] + w[rnog277b] - w[rnog278f] + w[rnog278b]
		- w[rnog279f] + w[rnog279b];

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
	W[sC6H5CH3] =  9.21340000e+01;
	W[sHO2] =  3.30080000e+01;
	W[sOH] =  1.70080000e+01;
	W[sH2O] =  1.80160000e+01;
	W[sO2] =  3.20000000e+01;
	W[sC2H4] =  2.80520000e+01;
	W[sH] =  1.00800000e+00;
	W[sC2H3] =  2.70440000e+01;
	W[sH2] =  2.01600000e+00;
	W[sCH2O] =  3.00260000e+01;
	W[sHCO] =  2.90180000e+01;
	W[sC2H6] =  3.00680000e+01;
	W[sC2H5] =  2.90600000e+01;
	W[sCH3O2] =  4.70340000e+01;
	W[sCH3] =  1.50340000e+01;
	W[sCH3O] =  3.10340000e+01;
	W[sCH3O2H] =  4.80420000e+01;
	W[sCH4] =  1.60420000e+01;
	W[sCO] =  2.80100000e+01;
	W[sCO2] =  4.40100000e+01;
	W[sO] =  1.60000000e+01;
	W[sH2O2] =  3.40160000e+01;
	W[sC2H2] =  2.60360000e+01;
	W[sCH2] =  1.40260000e+01;
	W[sCH2OH] =  3.10340000e+01;
	W[sHCCO] =  4.10280000e+01;
	W[sCH2CO] =  4.20360000e+01;
	W[sCH3CO] =  4.30440000e+01;
	W[sAXC3H5] =  4.10700000e+01;
	W[sC3H6] =  4.20780000e+01;
	W[sC3H5O] =  5.70700000e+01;
	W[sCH3OH] =  3.20420000e+01;
	W[sCH3CHO] =  4.40520000e+01;
	W[sNXC3H7] =  4.30860000e+01;
	W[sC5H9] =  6.91220000e+01;
	W[sPXC4H9] =  5.71120000e+01;
	W[s1XC5H10] =  7.01300000e+01;
	W[s1XC5H11] =  7.11380000e+01;
	W[sCH3COCH2O2] =  8.90700000e+01;
	W[sCH3COCH2] =  5.70700000e+01;
	W[sCH3COCH2O2H] =  9.00780000e+01;
	W[sCH3COCH2O] =  7.30700000e+01;
	W[sCH2CHO] =  4.30440000e+01;
	W[sNXC7H16] =  1.00198000e+02;
	W[s1XC7H15] =  9.91900000e+01;
	W[s2XC7H15] =  9.91900000e+01;
	W[s3XC7H15] =  9.91900000e+01;
	W[s4XC7H15] =  9.91900000e+01;
	W[s1XC7H14] =  9.81820000e+01;
	W[s2XC7H15O2] =  1.31190000e+02;
	W[s2X3XC7H14O2H] =  1.31190000e+02;
	W[s2X4XC7H14O2H] =  1.31190000e+02;
	W[s2X4XC7H14OOHO2] =  1.63190000e+02;
	W[sNC7KET24XC7H14O3] =  1.46182000e+02;
	W[sCH] =  1.30180000e+01;
	W[sC3H5XTXC3H5] =  4.10700000e+01;
	W[sIXC3H7] =  4.30860000e+01;
	W[sIXC4H8] =  5.61040000e+01;
	W[sIXC4H7] =  5.50960000e+01;
	W[sTXC4H9] =  5.71120000e+01;
	W[sIXC4H9] =  5.71120000e+01;
	W[sCH3COCH3] =  5.80780000e+01;
	W[sC2H5CO] =  5.70700000e+01;
	W[sNEOXC5H11] =  7.11380000e+01;
	W[sIXC4H7O] =  7.10960000e+01;
	W[sIXC4H8O] =  7.21040000e+01;
	W[sIC4H8OOHXIO2XC4H9O4] =  1.21112000e+02;
	W[sIC4H8OOHXIXC4H9O2] =  8.91120000e+01;
	W[sIC4H8OOHXTO2XC4H9O4] =  1.21112000e+02;
	W[sIC4H8OOHXTXC4H9O2] =  8.91120000e+01;
	W[sIC4KETIIXC4H8O3] =  1.04104000e+02;
	W[sTC4H8OOHXIXC4H9O2] =  8.91120000e+01;
	W[sTXC4H9O] =  7.31120000e+01;
	W[sIXC4H9O2] =  8.91120000e+01;
	W[sTXC4H9O2] =  8.91120000e+01;
	W[sTXC3H6CHO] =  7.10960000e+01;
	W[sTXC3H6OCHO] =  8.70960000e+01;
	W[sIXC3H5CHO] =  7.00880000e+01;
	W[sIXC3H7CHO] =  7.21040000e+01;
	W[sIXC4H7OOH] =  8.81040000e+01;
	W[sIXC3H5CO] =  6.90800000e+01;
	W[sYXC7H15] =  9.91900000e+01;
	W[sYXC7H14] =  9.81820000e+01;
	W[sXXC7H13] =  9.71740000e+01;
	W[sIXC8H18] =  1.14224000e+02;
	W[sAXC8H17] =  1.13216000e+02;
	W[sBXC8H17] =  1.13216000e+02;
	W[sCXC8H17] =  1.13216000e+02;
	W[sDXC8H17] =  1.13216000e+02;
	W[sAXC8H17O2] =  1.45216000e+02;
	W[sAXC8H17O2H] =  1.46224000e+02;
	W[sCXC8H17O2] =  1.45216000e+02;
	W[sCXC8H17O2H] =  1.46224000e+02;
	W[sAXC8H17O] =  1.29216000e+02;
	W[sCXC8H17O] =  1.29216000e+02;
	W[sAC8H16OOHXAXC8H17O2] =  1.45216000e+02;
	W[sAC8H16OOHXBXC8H17O2] =  1.45216000e+02;
	W[sAC8H16OOHXCXC8H17O2] =  1.45216000e+02;
	W[sCC8H16OOHXAXC8H17O2] =  1.45216000e+02;
	W[sAC8H16OOHXAO2XC8H17O4] =  1.77216000e+02;
	W[sAC8H16OOHXBO2XC8H17O4] =  1.77216000e+02;
	W[sAC8H16OOHXCO2XC8H17O4] =  1.77216000e+02;
	W[sCC8H16OOHXAO2XC8H17O4] =  1.77216000e+02;
	W[sIC8KETABXC8H16O3] =  1.60208000e+02;
	W[sC6H5CH2] =  9.11260000e+01;
	W[sC6H5] =  7.71000000e+01;
	W[sC6H6] =  7.81080000e+01;
	W[sC6H5OH] =  9.41080000e+01;
	W[sC6H5O] =  9.31000000e+01;
	W[sC6H5CHO] =  1.06118000e+02;
	W[sC6H5CO] =  1.05110000e+02;
	W[sC4H5XIXC4H5] =  5.30800000e+01;
	W[sC4H6] =  5.40880000e+01;
	W[sC6H5CH2O] =  1.07126000e+02;
	W[sC6H5C2H5] =  1.06160000e+02;
	W[sC6H5CH2OH] =  1.08134000e+02;
	W[sC5H5] =  6.50900000e+01;
	W[sC6H4O2] =  1.08092000e+02;
	W[sC5H6] =  6.60980000e+01;
	W[sC5H4O] =  8.00820000e+01;
	W[sC5H5O] =  8.10900000e+01;
	W[sC4H5XNXC4H5] =  5.30800000e+01;
	W[sCH2CHCHCO] =  6.80720000e+01;
	W[sCHCHCHO] =  5.50540000e+01;
	W[sHCCCO] =  5.30380000e+01;
	W[sC2H5OH] =  4.60680000e+01;
	W[sC2H4OH] =  4.50600000e+01;
	W[sCH3CHOH] =  4.50600000e+01;
	W[sCH3CH2O] =  4.50600000e+01;
	W[sCH3HCO] =  4.40520000e+01;
	W[sCH2HCO] =  4.30440000e+01;
	W[sHOC2H4O2] =  7.70600000e+01;
	W[sC2H] =  2.50280000e+01;
	W[sH2CC] =  2.60360000e+01;
	W[sC3H3] =  3.90540000e+01;
	W[sPXC3H4] =  4.00620000e+01;
	W[sAXC3H4] =  4.00620000e+01;
	W[s1XCH2] =  1.40260000e+01;
	W[sC] =  1.20100000e+01;
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
}

void GetSpeciesNames( char **names )
{
	strcpy( names[sN2], "N2                  " );
	strcpy( names[sAR], "AR                  " );
	strcpy( names[sC6H5CH3], "C6H5CH3             " );
	strcpy( names[sHO2], "HO2                 " );
	strcpy( names[sOH], "OH                  " );
	strcpy( names[sH2O], "H2O                 " );
	strcpy( names[sO2], "O2                  " );
	strcpy( names[sC2H4], "C2H4                " );
	strcpy( names[sH], "H                   " );
	strcpy( names[sC2H3], "C2H3                " );
	strcpy( names[sH2], "H2                  " );
	strcpy( names[sCH2O], "CH2O                " );
	strcpy( names[sHCO], "HCO                 " );
	strcpy( names[sC2H6], "C2H6                " );
	strcpy( names[sC2H5], "C2H5                " );
	strcpy( names[sCH3O2], "CH3O2               " );
	strcpy( names[sCH3], "CH3                 " );
	strcpy( names[sCH3O], "CH3O                " );
	strcpy( names[sCH3O2H], "CH3O2H              " );
	strcpy( names[sCH4], "CH4                 " );
	strcpy( names[sCO], "CO                  " );
	strcpy( names[sCO2], "CO2                 " );
	strcpy( names[sO], "O                   " );
	strcpy( names[sH2O2], "H2O2                " );
	strcpy( names[sC2H2], "C2H2                " );
	strcpy( names[sCH2], "CH2                 " );
	strcpy( names[sCH2OH], "CH2OH               " );
	strcpy( names[sHCCO], "HCCO                " );
	strcpy( names[sCH2CO], "CH2CO               " );
	strcpy( names[sCH3CO], "CH3CO               " );
	strcpy( names[sAXC3H5], "A-C3H5              " );
	strcpy( names[sC3H6], "C3H6                " );
	strcpy( names[sC3H5O], "C3H5O               " );
	strcpy( names[sCH3OH], "CH3OH               " );
	strcpy( names[sCH3CHO], "CH3CHO              " );
	strcpy( names[sNXC3H7], "N-C3H7              " );
	strcpy( names[sC5H9], "C5H9                " );
	strcpy( names[sPXC4H9], "P-C4H9              " );
	strcpy( names[s1XC5H10], "1-C5H10             " );
	strcpy( names[s1XC5H11], "1-C5H11             " );
	strcpy( names[sCH3COCH2O2], "CH3COCH2O2          " );
	strcpy( names[sCH3COCH2], "CH3COCH2            " );
	strcpy( names[sCH3COCH2O2H], "CH3COCH2O2H         " );
	strcpy( names[sCH3COCH2O], "CH3COCH2O           " );
	strcpy( names[sCH2CHO], "CH2CHO              " );
	strcpy( names[sNXC7H16], "N-C7H16             " );
	strcpy( names[s1XC7H15], "1-C7H15             " );
	strcpy( names[s2XC7H15], "2-C7H15             " );
	strcpy( names[s3XC7H15], "3-C7H15             " );
	strcpy( names[s4XC7H15], "4-C7H15             " );
	strcpy( names[s1XC7H14], "1-C7H14             " );
	strcpy( names[s2XC7H15O2], "2-C7H15O2           " );
	strcpy( names[s2X3XC7H14O2H], "2-3-C7H14O2H        " );
	strcpy( names[s2X4XC7H14O2H], "2-4-C7H14O2H        " );
	strcpy( names[s2X4XC7H14OOHO2], "2-4-C7H14OOHO2      " );
	strcpy( names[sNC7KET24XC7H14O3], "NC7KET24-C7H14O3    " );
	strcpy( names[sCH], "CH                  " );
	strcpy( names[sC3H5XTXC3H5], "C3H5-T-C3H5         " );
	strcpy( names[sIXC3H7], "I-C3H7              " );
	strcpy( names[sIXC4H8], "I-C4H8              " );
	strcpy( names[sIXC4H7], "I-C4H7              " );
	strcpy( names[sTXC4H9], "T-C4H9              " );
	strcpy( names[sIXC4H9], "I-C4H9              " );
	strcpy( names[sCH3COCH3], "CH3COCH3            " );
	strcpy( names[sC2H5CO], "C2H5CO              " );
	strcpy( names[sNEOXC5H11], "NEO-C5H11           " );
	strcpy( names[sIXC4H7O], "I-C4H7O             " );
	strcpy( names[sIXC4H8O], "I-C4H8O             " );
	strcpy( names[sIC4H8OOHXIO2XC4H9O4], "IC4H8OOH-IO2-C4H9O4 " );
	strcpy( names[sIC4H8OOHXIXC4H9O2], "IC4H8OOH-I-C4H9O2   " );
	strcpy( names[sIC4H8OOHXTO2XC4H9O4], "IC4H8OOH-TO2-C4H9O4 " );
	strcpy( names[sIC4H8OOHXTXC4H9O2], "IC4H8OOH-T-C4H9O2   " );
	strcpy( names[sIC4KETIIXC4H8O3], "IC4KETII-C4H8O3     " );
	strcpy( names[sTC4H8OOHXIXC4H9O2], "TC4H8OOH-I-C4H9O2   " );
	strcpy( names[sTXC4H9O], "T-C4H9O             " );
	strcpy( names[sIXC4H9O2], "I-C4H9O2            " );
	strcpy( names[sTXC4H9O2], "T-C4H9O2            " );
	strcpy( names[sTXC3H6CHO], "T-C3H6CHO           " );
	strcpy( names[sTXC3H6OCHO], "T-C3H6OCHO          " );
	strcpy( names[sIXC3H5CHO], "I-C3H5CHO           " );
	strcpy( names[sIXC3H7CHO], "I-C3H7CHO           " );
	strcpy( names[sIXC4H7OOH], "I-C4H7OOH           " );
	strcpy( names[sIXC3H5CO], "I-C3H5CO            " );
	strcpy( names[sYXC7H15], "Y-C7H15             " );
	strcpy( names[sYXC7H14], "Y-C7H14             " );
	strcpy( names[sXXC7H13], "X-C7H13             " );
	strcpy( names[sIXC8H18], "I-C8H18             " );
	strcpy( names[sAXC8H17], "A-C8H17             " );
	strcpy( names[sBXC8H17], "B-C8H17             " );
	strcpy( names[sCXC8H17], "C-C8H17             " );
	strcpy( names[sDXC8H17], "D-C8H17             " );
	strcpy( names[sAXC8H17O2], "A-C8H17O2           " );
	strcpy( names[sAXC8H17O2H], "A-C8H17O2H          " );
	strcpy( names[sCXC8H17O2], "C-C8H17O2           " );
	strcpy( names[sCXC8H17O2H], "C-C8H17O2H          " );
	strcpy( names[sAXC8H17O], "A-C8H17O            " );
	strcpy( names[sCXC8H17O], "C-C8H17O            " );
	strcpy( names[sAC8H16OOHXAXC8H17O2], "AC8H16OOH-A-C8H17O2 " );
	strcpy( names[sAC8H16OOHXBXC8H17O2], "AC8H16OOH-B-C8H17O2 " );
	strcpy( names[sAC8H16OOHXCXC8H17O2], "AC8H16OOH-C-C8H17O2 " );
	strcpy( names[sCC8H16OOHXAXC8H17O2], "CC8H16OOH-A-C8H17O2 " );
	strcpy( names[sAC8H16OOHXAO2XC8H17O4], "AC8H16OOH-AO2-C8H17O4" );
	strcpy( names[sAC8H16OOHXBO2XC8H17O4], "AC8H16OOH-BO2-C8H17O4" );
	strcpy( names[sAC8H16OOHXCO2XC8H17O4], "AC8H16OOH-CO2-C8H17O4" );
	strcpy( names[sCC8H16OOHXAO2XC8H17O4], "CC8H16OOH-AO2-C8H17O4" );
	strcpy( names[sIC8KETABXC8H16O3], "IC8KETAB-C8H16O3    " );
	strcpy( names[sC6H5CH2], "C6H5CH2             " );
	strcpy( names[sC6H5], "C6H5                " );
	strcpy( names[sC6H6], "C6H6                " );
	strcpy( names[sC6H5OH], "C6H5OH              " );
	strcpy( names[sC6H5O], "C6H5O               " );
	strcpy( names[sC6H5CHO], "C6H5CHO             " );
	strcpy( names[sC6H5CO], "C6H5CO              " );
	strcpy( names[sC4H5XIXC4H5], "C4H5-I-C4H5         " );
	strcpy( names[sC4H6], "C4H6                " );
	strcpy( names[sC6H5CH2O], "C6H5CH2O            " );
	strcpy( names[sC6H5C2H5], "C6H5C2H5            " );
	strcpy( names[sC6H5CH2OH], "C6H5CH2OH           " );
	strcpy( names[sC5H5], "C5H5                " );
	strcpy( names[sC6H4O2], "C6H4O2              " );
	strcpy( names[sC5H6], "C5H6                " );
	strcpy( names[sC5H4O], "C5H4O               " );
	strcpy( names[sC5H5O], "C5H5O               " );
	strcpy( names[sC4H5XNXC4H5], "C4H5-N-C4H5         " );
	strcpy( names[sCH2CHCHCO], "CH2CHCHCO           " );
	strcpy( names[sCHCHCHO], "CHCHCHO             " );
	strcpy( names[sHCCCO], "HCCCO               " );
	strcpy( names[sC2H5OH], "C2H5OH              " );
	strcpy( names[sC2H4OH], "C2H4OH              " );
	strcpy( names[sCH3CHOH], "CH3CHOH             " );
	strcpy( names[sCH3CH2O], "CH3CH2O             " );
	strcpy( names[sCH3HCO], "CH3HCO              " );
	strcpy( names[sCH2HCO], "CH2HCO              " );
	strcpy( names[sHOC2H4O2], "HOC2H4O2            " );
	strcpy( names[sC2H], "C2H                 " );
	strcpy( names[sH2CC], "H2CC                " );
	strcpy( names[sC3H3], "C3H3                " );
	strcpy( names[sPXC3H4], "P-C3H4              " );
	strcpy( names[sAXC3H4], "A-C3H4              " );
	strcpy( names[s1XCH2], "1-CH2               " );
	strcpy( names[sC], "C                   " );
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
}



void ComputeThermoData( double *h, double *cp, double T )
{
/*
	This function computes enthalpy 'h' and heat capacity 'cp' as
	function of temperature 'T' for all non steady state species
	in units [J/kg] and [J/kg K], respectively.
	The parameter h and cp should provide workspace of length 158 */


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
		h[sC6H5CH3] =  9.02418217e+01 * (
			T * (  1.29400340e+01 + T * (  1.33456435e-02
			+ T * ( -3.22795017e-06 + T * (  3.93465725e-10
			+ T * -1.89327202e-14 ) ) ) ) -6.97649080e+02 );
		cp[sC6H5CH3] =  9.02418217e+01 * (
			 1.29400340e+01 + T * (  2.66912870e-02 
			+ T * ( -9.68385050e-06 + T * (  1.57386290e-09
			+ T * -9.46636010e-14 ) ) ) );
		h[sHO2] =  2.51888633e+02 * (
			T * (  3.00201579e+00 + T * (  2.19276212e-03
			+ T * ( -6.58647397e-07 + T * (  8.30225242e-11
			+ T *  0.00000000e+00 ) ) ) ) +  8.44944001e+02 );
		cp[sHO2] =  2.51888633e+02 * (
			 3.00201579e+00 + T * (  4.38552424e-03 
			+ T * ( -1.97594219e-06 + T * (  3.32090097e-10
			+ T *  0.00000000e+00 ) ) ) );
		h[sOH] =  4.88848777e+02 * (
			T * (  3.69039275e+00 + T * ( -2.90401159e-04
			+ T * (  2.54434419e-07 + T * ( -4.50655148e-11
			+ T *  0.00000000e+00 ) ) ) ) +  3.67872158e+03 );
		cp[sOH] =  4.88848777e+02 * (
			 3.69039275e+00 + T * ( -5.80802318e-04 
			+ T * (  7.63303256e-07 + T * ( -1.80262059e-10
			+ T *  0.00000000e+00 ) ) ) );
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
		h[sC2H4] =  2.96390275e+02 * (
			T * (  3.52841900e+00 + T * (  5.74259000e-03
			+ T * ( -1.47279500e-06 + T * (  1.96115025e-10
			+ T * -1.05336960e-14 ) ) ) ) +  4.42828900e+03 );
		cp[sC2H4] =  2.96390275e+02 * (
			 3.52841900e+00 + T * (  1.14851800e-02 
			+ T * ( -4.41838500e-06 + T * (  7.84460100e-10
			+ T * -5.26684800e-14 ) ) ) );
		h[sH] =  8.24835317e+03 * (
			T * (  2.50104422e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  2.54747466e+04 );
		cp[sH] =  8.24835317e+03 * (
			 2.50104422e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sC2H3] =  3.07437509e+02 * (
			T * (  5.93346800e+00 + T * (  2.00887300e-03
			+ T * ( -1.32224667e-07 + T * ( -3.60316750e-11
			+ T *  4.75728800e-15 ) ) ) ) +  3.18543500e+04 );
		cp[sC2H3] =  3.07437509e+02 * (
			 5.93346800e+00 + T * (  4.01774600e-03 
			+ T * ( -3.96674000e-07 + T * ( -1.44126700e-10
			+ T *  2.37864400e-14 ) ) ) );
		h[sH2] =  4.12417659e+03 * (
			T * (  2.99142300e+00 + T * (  3.50032200e-04
			+ T * ( -1.87794300e-08 + T * ( -2.30789450e-12
			+ T *  3.16550400e-16 ) ) ) ) -8.35034000e+02 );
		cp[sH2] =  4.12417659e+03 * (
			 2.99142300e+00 + T * (  7.00064400e-04 
			+ T * ( -5.63382900e-08 + T * ( -9.23157800e-12
			+ T *  1.58275200e-15 ) ) ) );
		h[sCH2O] =  2.76904683e+02 * (
			T * (  2.99560600e+00 + T * (  3.34066050e-03
			+ T * ( -8.76318333e-07 + T * (  1.18428825e-10
			+ T * -6.42503400e-15 ) ) ) ) -1.53203700e+04 );
		cp[sCH2O] =  2.76904683e+02 * (
			 2.99560600e+00 + T * (  6.68132100e-03 
			+ T * ( -2.62895500e-06 + T * (  4.73715300e-10
			+ T * -3.21251700e-14 ) ) ) );
		h[sHCO] =  2.86523537e+02 * (
			T * (  3.55727100e+00 + T * (  1.67278650e-03
			+ T * ( -4.45002000e-07 + T * (  6.17643250e-11
			+ T * -3.42770200e-15 ) ) ) ) +  3.91632400e+03 );
		cp[sHCO] =  2.86523537e+02 * (
			 3.55727100e+00 + T * (  3.34557300e-03 
			+ T * ( -1.33500600e-06 + T * (  2.47057300e-10
			+ T * -1.71385100e-14 ) ) ) );
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
		h[sCH3O2] =  1.76772973e+02 * (
			T * (  3.72224006e+00 + T * (  5.77152465e-03
			+ T * ( -1.51686254e-06 + T * (  1.67577132e-10
			+ T *  0.00000000e+00 ) ) ) ) -5.34692750e+02 );
		cp[sCH3O2] =  1.76772973e+02 * (
			 3.72224006e+00 + T * (  1.15430493e-02 
			+ T * ( -4.55058763e-06 + T * (  6.70308528e-10
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH3] =  5.53035786e+02 * (
			T * (  2.84405200e+00 + T * (  3.06898700e-03
			+ T * ( -7.43448333e-07 + T * (  9.46290250e-11
			+ T * -4.90431800e-15 ) ) ) ) +  1.64378100e+04 );
		cp[sCH3] =  5.53035786e+02 * (
			 2.84405200e+00 + T * (  6.13797400e-03 
			+ T * ( -2.23034500e-06 + T * (  3.78516100e-10
			+ T * -2.45215900e-14 ) ) ) );
		h[sCH3O] =  2.67910679e+02 * (
			T * (  3.77080000e+00 + T * (  3.93574850e-03
			+ T * ( -8.85461333e-07 + T * (  9.86107750e-11
			+ T * -4.22523200e-15 ) ) ) ) +  1.27832500e+02 );
		cp[sCH3O] =  2.67910679e+02 * (
			 3.77080000e+00 + T * (  7.87149700e-03 
			+ T * ( -2.65638400e-06 + T * (  3.94443100e-10
			+ T * -2.11261600e-14 ) ) ) );
		h[sCH3O2H] =  1.73063986e+02 * (
			T * (  4.46640989e+00 + T * (  7.55106630e-03
			+ T * ( -2.26701845e-06 + T * (  2.83131345e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.80236075e+04 );
		cp[sCH3O2H] =  1.73063986e+02 * (
			 4.46640989e+00 + T * (  1.51021326e-02 
			+ T * ( -6.80105535e-06 + T * (  1.13252538e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH4] =  5.18285750e+02 * (
			T * (  1.68347900e+00 + T * (  5.11862000e-03
			+ T * ( -1.29170967e-06 + T * (  1.69639625e-10
			+ T * -9.00684600e-15 ) ) ) ) -1.00807900e+04 );
		cp[sCH4] =  5.18285750e+02 * (
			 1.68347900e+00 + T * (  1.02372400e-02 
			+ T * ( -3.87512900e-06 + T * (  6.78558500e-10
			+ T * -4.50342300e-14 ) ) ) );
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
		h[sO] =  5.19646250e+02 * (
			T * (  2.54206000e+00 + T * ( -1.37753100e-05
			+ T * ( -1.03426767e-09 + T * (  1.13776675e-12
			+ T * -8.73610400e-17 ) ) ) ) +  2.92308000e+04 );
		cp[sO] =  5.19646250e+02 * (
			 2.54206000e+00 + T * ( -2.75506200e-05 
			+ T * ( -3.10280300e-09 + T * (  4.55106700e-12
			+ T * -4.36805200e-16 ) ) ) );
		h[sH2O2] =  2.44424389e+02 * (
			T * (  4.57316700e+00 + T * (  2.16806800e-03
			+ T * ( -4.91563000e-07 + T * (  5.87226000e-11
			+ T * -2.86330800e-15 ) ) ) ) -1.80069600e+04 );
		cp[sH2O2] =  2.44424389e+02 * (
			 4.57316700e+00 + T * (  4.33613600e-03 
			+ T * ( -1.47468900e-06 + T * (  2.34890400e-10
			+ T * -1.43165400e-14 ) ) ) );
		h[sC2H2] =  3.19340144e+02 * (
			T * (  4.43677000e+00 + T * (  2.68801950e-03
			+ T * ( -6.37605667e-07 + T * (  8.21594750e-11
			+ T * -4.31342000e-15 ) ) ) ) +  2.56676600e+04 );
		cp[sC2H2] =  3.19340144e+02 * (
			 4.43677000e+00 + T * (  5.37603900e-03 
			+ T * ( -1.91281700e-06 + T * (  3.28637900e-10
			+ T * -2.15671000e-14 ) ) ) );
		h[sCH2] =  5.92780550e+02 * (
			T * (  3.63640800e+00 + T * (  9.66528500e-04
			+ T * ( -5.62338667e-08 + T * ( -2.52474750e-11
			+ T *  3.61651200e-15 ) ) ) ) +  4.53413400e+04 );
		cp[sCH2] =  5.92780550e+02 * (
			 3.63640800e+00 + T * (  1.93305700e-03 
			+ T * ( -1.68701600e-07 + T * ( -1.00989900e-10
			+ T *  1.80825600e-14 ) ) ) );
		h[sCH2OH] =  2.67910679e+02 * (
			T * (  6.32752000e+00 + T * (  1.80413550e-03
			+ T * ( -1.06718233e-07 + T * ( -4.84687500e-11
			+ T *  7.01941000e-15 ) ) ) ) -4.47450900e+03 );
		cp[sCH2OH] =  2.67910679e+02 * (
			 6.32752000e+00 + T * (  3.60827100e-03 
			+ T * ( -3.20154700e-07 + T * ( -1.93875000e-10
			+ T *  3.50970500e-14 ) ) ) );
		h[sHCCO] =  2.02650385e+02 * (
			T * (  6.75807300e+00 + T * (  1.00020000e-03
			+ T * ( -6.75869000e-08 + T * ( -2.60283000e-11
			+ T *  3.93033000e-15 ) ) ) ) +  1.90151300e+04 );
		cp[sHCCO] =  2.02650385e+02 * (
			 6.75807300e+00 + T * (  2.00040000e-03 
			+ T * ( -2.02760700e-07 + T * ( -1.04113200e-10
			+ T *  1.96516500e-14 ) ) ) );
		h[sCH2CO] =  1.97790941e+02 * (
			T * (  6.03881700e+00 + T * (  2.90242000e-03
			+ T * ( -6.40318000e-07 + T * (  6.98621250e-11
			+ T * -2.91773600e-15 ) ) ) ) -8.58340200e+03 );
		cp[sCH2CO] =  1.97790941e+02 * (
			 6.03881700e+00 + T * (  5.80484000e-03 
			+ T * ( -1.92095400e-06 + T * (  2.79448500e-10
			+ T * -1.45886800e-14 ) ) ) );
		h[sCH3CO] =  1.93159093e+02 * (
			T * (  5.61227900e+00 + T * (  4.22494300e-03
			+ T * ( -9.51382333e-07 + T * (  1.05959400e-10
			+ T * -4.53680800e-15 ) ) ) ) -5.18786300e+03 );
		cp[sCH3CO] =  1.93159093e+02 * (
			 5.61227900e+00 + T * (  8.44988600e-03 
			+ T * ( -2.85414700e-06 + T * (  4.23837600e-10
			+ T * -2.26840400e-14 ) ) ) );
		h[sAXC3H5] =  2.02443146e+02 * (
			T * (  9.65153900e+00 + T * (  4.03779800e-03
			+ T * ( -2.65514133e-07 + T * ( -1.16267400e-10
			+ T *  1.72065600e-14 ) ) ) ) +  1.53009600e+04 );
		cp[sAXC3H5] =  2.02443146e+02 * (
			 9.65153900e+00 + T * (  8.07559600e-03 
			+ T * ( -7.96542400e-07 + T * ( -4.65069600e-10
			+ T *  8.60328000e-14 ) ) ) );
		h[sC3H6] =  1.97593517e+02 * (
			T * (  1.68134957e+00 + T * (  1.24845359e-02
			+ T * ( -3.70920137e-06 + T * (  4.61033815e-10
			+ T *  0.00000000e+00 ) ) ) ) +  7.40726095e+02 );
		cp[sC3H6] =  1.97593517e+02 * (
			 1.68134957e+00 + T * (  2.49690718e-02 
			+ T * ( -1.11276041e-05 + T * (  1.84413526e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC3H5O] =  1.45686701e+02 * (
			T * (  3.39074577e+00 + T * (  1.20650810e-02
			+ T * ( -3.78836313e-06 + T * (  4.94752345e-10
			+ T *  0.00000000e+00 ) ) ) ) +  9.00757452e+03 );
		cp[sC3H5O] =  1.45686701e+02 * (
			 3.39074577e+00 + T * (  2.41301620e-02 
			+ T * ( -1.13650894e-05 + T * (  1.97900938e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH3OH] =  2.59482554e+02 * (
			T * (  4.02906100e+00 + T * (  4.68829650e-03
			+ T * ( -1.01675133e-06 + T * (  1.08969825e-10
			+ T * -4.44944600e-15 ) ) ) ) -2.61579100e+04 );
		cp[sCH3OH] =  2.59482554e+02 * (
			 4.02906100e+00 + T * (  9.37659300e-03 
			+ T * ( -3.05025400e-06 + T * (  4.35879300e-10
			+ T * -2.22472300e-14 ) ) ) );
		h[sCH3CHO] =  1.88739217e+02 * (
			T * (  5.86865000e+00 + T * (  5.39712000e-03
			+ T * ( -1.21517667e-06 + T * (  1.35322800e-10
			+ T * -5.79368800e-15 ) ) ) ) -2.26456900e+04 );
		cp[sCH3CHO] =  1.88739217e+02 * (
			 5.86865000e+00 + T * (  1.07942400e-02 
			+ T * ( -3.64553000e-06 + T * (  5.41291200e-10
			+ T * -2.89684400e-14 ) ) ) );
		h[sNXC3H7] =  1.92970803e+02 * (
			T * (  7.97829100e+00 + T * (  7.88056500e-03
			+ T * ( -1.72441433e-06 + T * (  1.86097300e-10
			+ T * -7.64995600e-15 ) ) ) ) +  7.57940200e+03 );
		cp[sNXC3H7] =  1.92970803e+02 * (
			 7.97829100e+00 + T * (  1.57611300e-02 
			+ T * ( -5.17324300e-06 + T * (  7.44389200e-10
			+ T * -3.82497800e-14 ) ) ) );
		h[sC5H9] =  1.20285003e+02 * (
			T * (  3.78447384e+00 + T * (  1.96499757e-02
			+ T * ( -5.93368333e-06 + T * (  7.46485712e-10
			+ T *  0.00000000e+00 ) ) ) ) +  1.12889042e+04 );
		cp[sC5H9] =  1.20285003e+02 * (
			 3.78447384e+00 + T * (  3.92999515e-02 
			+ T * ( -1.78010500e-05 + T * (  2.98594285e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sPXC4H9] =  1.45579563e+02 * (
			T * (  3.81812330e+00 + T * (  1.70395244e-02
			+ T * ( -4.97118110e-06 + T * (  6.08013118e-10
			+ T *  0.00000000e+00 ) ) ) ) +  6.65901304e+03 );
		cp[sPXC4H9] =  1.45579563e+02 * (
			 3.81812330e+00 + T * (  3.40790489e-02 
			+ T * ( -1.49135433e-05 + T * (  2.43205247e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[s1XC5H10] =  1.18556110e+02 * (
			T * (  3.98580522e+00 + T * (  2.06214993e-02
			+ T * ( -6.14634990e-06 + T * (  7.65388102e-10
			+ T *  0.00000000e+00 ) ) ) ) -5.70112071e+03 );
		cp[s1XC5H10] =  1.18556110e+02 * (
			 3.98580522e+00 + T * (  4.12429986e-02 
			+ T * ( -1.84390497e-05 + T * (  3.06155241e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[s1XC5H11] =  1.16876212e+02 * (
			T * (  4.88920629e+00 + T * (  2.11417269e-02
			+ T * ( -6.19477000e-06 + T * (  7.60311908e-10
			+ T *  0.00000000e+00 ) ) ) ) +  3.43475468e+03 );
		cp[s1XC5H11] =  1.16876212e+02 * (
			 4.88920629e+00 + T * (  4.22834537e-02 
			+ T * ( -1.85843100e-05 + T * (  3.04124763e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH3COCH2O2] =  9.33461323e+01 * (
			T * (  6.39329449e+00 + T * (  1.25899650e-02
			+ T * ( -3.65240633e-06 + T * (  4.39950632e-10
			+ T *  0.00000000e+00 ) ) ) ) -2.07465812e+04 );
		cp[sCH3COCH2O2] =  9.33461323e+01 * (
			 6.39329449e+00 + T * (  2.51799300e-02 
			+ T * ( -1.09572190e-05 + T * (  1.75980253e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH3COCH2] =  1.45686701e+02 * (
			T * (  4.27637519e+00 + T * (  1.11174643e-02
			+ T * ( -3.36394730e-06 + T * (  4.22492645e-10
			+ T *  0.00000000e+00 ) ) ) ) -5.97887937e+03 );
		cp[sCH3COCH2] =  1.45686701e+02 * (
			 4.27637519e+00 + T * (  2.22349285e-02 
			+ T * ( -1.00918419e-05 + T * (  1.68997058e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH3COCH2O2H] =  9.23015609e+01 * (
			T * (  7.24028399e+00 + T * (  1.42563754e-02
			+ T * ( -4.35283280e-06 + T * (  5.47627642e-10
			+ T *  0.00000000e+00 ) ) ) ) -3.94177289e+04 );
		cp[sCH3COCH2O2H] =  9.23015609e+01 * (
			 7.24028399e+00 + T * (  2.85127508e-02 
			+ T * ( -1.30584984e-05 + T * (  2.19051057e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH3COCH2O] =  1.13785959e+02 * (
			T * (  7.47517424e+00 + T * (  1.41867268e-02
			+ T * ( -5.12626773e-06 + T * (  7.26424122e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.34371904e+04 );
		cp[sCH3COCH2O] =  1.13785959e+02 * (
			 7.47517424e+00 + T * (  2.83734536e-02 
			+ T * ( -1.53788032e-05 + T * (  2.90569649e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH2CHO] =  1.93159093e+02 * (
			T * (  5.97567000e+00 + T * (  4.06529550e-03
			+ T * ( -9.14541333e-07 + T * (  1.01757600e-10
			+ T * -4.35203400e-15 ) ) ) ) +  4.90321800e+02 );
		cp[sCH2CHO] =  1.93159093e+02 * (
			 5.97567000e+00 + T * (  8.13059100e-03 
			+ T * ( -2.74362400e-06 + T * (  4.07030400e-10
			+ T * -2.17601700e-14 ) ) ) );
		h[sNXC7H16] =  8.29791014e+01 * (
			T * (  5.14079241e+00 + T * (  3.26539335e-02
			+ T * ( -9.82758747e-06 + T * (  1.23431681e-09
			+ T *  0.00000000e+00 ) ) ) ) -2.72533890e+04 );
		cp[sNXC7H16] =  8.29791014e+01 * (
			 5.14079241e+00 + T * (  6.53078671e-02 
			+ T * ( -2.94827624e-05 + T * (  4.93726726e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[s1XC7H15] =  8.38223611e+01 * (
			T * (  5.67638200e+00 + T * (  3.07314223e-02
			+ T * ( -9.25845533e-06 + T * (  1.16367723e-09
			+ T *  0.00000000e+00 ) ) ) ) -2.58110033e+03 );
		cp[s1XC7H15] =  8.38223611e+01 * (
			 5.67638200e+00 + T * (  6.14628445e-02 
			+ T * ( -2.77753660e-05 + T * (  4.65470891e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[s2XC7H15] =  8.38223611e+01 * (
			T * (  3.74721159e+00 + T * (  3.24672581e-02
			+ T * ( -1.00447008e-05 + T * (  1.29354536e-09
			+ T *  0.00000000e+00 ) ) ) ) -3.37018357e+03 );
		cp[s2XC7H15] =  8.38223611e+01 * (
			 3.74721159e+00 + T * (  6.49345162e-02 
			+ T * ( -3.01341025e-05 + T * (  5.17418142e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[s3XC7H15] =  8.38223611e+01 * (
			T * (  3.74721159e+00 + T * (  3.24672581e-02
			+ T * ( -1.00447008e-05 + T * (  1.29354536e-09
			+ T *  0.00000000e+00 ) ) ) ) -3.37018357e+03 );
		cp[s3XC7H15] =  8.38223611e+01 * (
			 3.74721159e+00 + T * (  6.49345162e-02 
			+ T * ( -3.01341025e-05 + T * (  5.17418142e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[s4XC7H15] =  8.38223611e+01 * (
			T * (  3.74721159e+00 + T * (  3.24672581e-02
			+ T * ( -1.00447008e-05 + T * (  1.29354536e-09
			+ T *  0.00000000e+00 ) ) ) ) -3.37018357e+03 );
		cp[s4XC7H15] =  8.38223611e+01 * (
			 3.74721159e+00 + T * (  6.49345162e-02 
			+ T * ( -3.01341025e-05 + T * (  5.17418142e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[s1XC7H14] =  8.46829358e+01 * (
			T * (  5.97612676e+00 + T * (  2.90245262e-02
			+ T * ( -8.70960413e-06 + T * (  1.09150687e-09
			+ T *  0.00000000e+00 ) ) ) ) -1.20977632e+04 );
		cp[s1XC7H14] =  8.46829358e+01 * (
			 5.97612676e+00 + T * (  5.80490524e-02 
			+ T * ( -2.61288124e-05 + T * (  4.36602746e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[s2XC7H15O2] =  6.33763244e+01 * (
			T * (  9.57252117e+00 + T * (  3.10811618e-02
			+ T * ( -9.14153807e-06 + T * (  1.12196996e-09
			+ T *  0.00000000e+00 ) ) ) ) -2.39185444e+04 );
		cp[s2XC7H15O2] =  6.33763244e+01 * (
			 9.57252117e+00 + T * (  6.21623236e-02 
			+ T * ( -2.74246142e-05 + T * (  4.48787986e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[s2X3XC7H14O2H] =  6.33763244e+01 * (
			T * (  9.35686979e+00 + T * (  3.23298594e-02
			+ T * ( -1.00497718e-05 + T * (  1.29834269e-09
			+ T *  0.00000000e+00 ) ) ) ) -1.76700972e+04 );
		cp[s2X3XC7H14O2H] =  6.33763244e+01 * (
			 9.35686979e+00 + T * (  6.46597189e-02 
			+ T * ( -3.01493154e-05 + T * (  5.19337077e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[s2X4XC7H14O2H] =  6.33763244e+01 * (
			T * (  9.35686979e+00 + T * (  3.23298594e-02
			+ T * ( -1.00497718e-05 + T * (  1.29834269e-09
			+ T *  0.00000000e+00 ) ) ) ) -1.76700972e+04 );
		cp[s2X4XC7H14O2H] =  6.33763244e+01 * (
			 9.35686979e+00 + T * (  6.46597189e-02 
			+ T * ( -3.01493154e-05 + T * (  5.19337077e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[s2X4XC7H14OOHO2] =  5.09488326e+01 * (
			T * (  1.26922886e+01 + T * (  3.46523221e-02
			+ T * ( -1.08959172e-05 + T * (  1.41029433e-09
			+ T *  0.00000000e+00 ) ) ) ) -3.76787486e+04 );
		cp[s2X4XC7H14OOHO2] =  5.09488326e+01 * (
			 1.26922886e+01 + T * (  6.93046441e-02 
			+ T * ( -3.26877516e-05 + T * (  5.64117730e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sNC7KET24XC7H14O3] =  5.68766332e+01 * (
			T * (  1.14177280e+01 + T * (  3.09063655e-02
			+ T * ( -9.36374100e-06 + T * (  1.17520390e-09
			+ T *  0.00000000e+00 ) ) ) ) -5.49015428e+04 );
		cp[sNC7KET24XC7H14O3] =  5.68766332e+01 * (
			 1.14177280e+01 + T * (  6.18127309e-02 
			+ T * ( -2.80912230e-05 + T * (  4.70081561e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH] =  6.38680289e+02 * (
			T * (  2.19622300e+00 + T * (  1.17019050e-03
			+ T * ( -2.35273367e-07 + T * (  2.25189550e-11
			+ T * -7.71008000e-16 ) ) ) ) +  7.08672300e+04 );
		cp[sCH] =  6.38680289e+02 * (
			 2.19622300e+00 + T * (  2.34038100e-03 
			+ T * ( -7.05820100e-07 + T * (  9.00758200e-11
			+ T * -3.85504000e-15 ) ) ) );
		h[sC3H5XTXC3H5] =  2.02443146e+02 * (
			T * (  1.94476358e+00 + T * (  1.07086382e-02
			+ T * ( -3.18351890e-06 + T * (  3.95877935e-10
			+ T *  0.00000000e+00 ) ) ) ) +  2.95295307e+04 );
		cp[sC3H5XTXC3H5] =  2.02443146e+02 * (
			 1.94476358e+00 + T * (  2.14172765e-02 
			+ T * ( -9.55055670e-06 + T * (  1.58351174e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIXC3H7] =  1.92970803e+02 * (
			T * (  8.06336900e+00 + T * (  7.87244000e-03
			+ T * ( -1.72746400e-06 + T * (  1.86931125e-10
			+ T * -7.70884400e-15 ) ) ) ) +  5.31387100e+03 );
		cp[sIXC3H7] =  1.92970803e+02 * (
			 8.06336900e+00 + T * (  1.57448800e-02 
			+ T * ( -5.18239200e-06 + T * (  7.47724500e-10
			+ T * -3.85442200e-14 ) ) ) );
		h[sIXC4H8] =  1.48195138e+02 * (
			T * (  2.86958571e+00 + T * (  1.64824603e-02
			+ T * ( -4.88104033e-06 + T * (  6.04074228e-10
			+ T *  0.00000000e+00 ) ) ) ) -4.22675236e+03 );
		cp[sIXC4H8] =  1.48195138e+02 * (
			 2.86958571e+00 + T * (  3.29649207e-02 
			+ T * ( -1.46431210e-05 + T * (  2.41629691e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIXC4H7] =  1.50906418e+02 * (
			T * (  4.86718299e+00 + T * (  1.37705580e-02
			+ T * ( -3.99386733e-06 + T * (  4.84570495e-10
			+ T *  0.00000000e+00 ) ) ) ) +  1.34120718e+04 );
		cp[sIXC4H7] =  1.50906418e+02 * (
			 4.86718299e+00 + T * (  2.75411161e-02 
			+ T * ( -1.19816020e-05 + T * (  1.93828198e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sTXC4H9] =  1.45579563e+02 * (
			T * ( -1.58631603e+00 + T * (  2.11173753e-02
			+ T * ( -6.44271080e-06 + T * (  8.14133765e-10
			+ T *  0.00000000e+00 ) ) ) ) +  4.56608047e+03 );
		cp[sTXC4H9] =  1.45579563e+02 * (
			-1.58631603e+00 + T * (  4.22347506e-02 
			+ T * ( -1.93281324e-05 + T * (  3.25653506e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIXC4H9] =  1.45579563e+02 * (
			T * (  3.77508718e+00 + T * (  1.71621240e-02
			+ T * ( -5.00002717e-06 + T * (  6.08502173e-10
			+ T *  0.00000000e+00 ) ) ) ) +  5.51336849e+03 );
		cp[sIXC4H9] =  1.45579563e+02 * (
			 3.77508718e+00 + T * (  3.43242480e-02 
			+ T * ( -1.50000815e-05 + T * (  2.43400869e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH3COCH3] =  1.43158167e+02 * (
			T * (  1.07057510e+00 + T * (  1.52820507e-02
			+ T * ( -4.79495923e-06 + T * (  6.18846660e-10
			+ T *  0.00000000e+00 ) ) ) ) -2.75121809e+04 );
		cp[sCH3COCH3] =  1.43158167e+02 * (
			 1.07057510e+00 + T * (  3.05641014e-02 
			+ T * ( -1.43848777e-05 + T * (  2.47538664e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC2H5CO] =  1.45686701e+02 * (
			T * (  3.47466548e+00 + T * (  1.16103793e-02
			+ T * ( -3.55000970e-06 + T * (  4.50665342e-10
			+ T *  0.00000000e+00 ) ) ) ) -5.84231374e+03 );
		cp[sC2H5CO] =  1.45686701e+02 * (
			 3.47466548e+00 + T * (  2.32207586e-02 
			+ T * ( -1.06500291e-05 + T * (  1.80266137e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sNEOXC5H11] =  1.16876212e+02 * (
			T * (  5.77200574e+00 + T * (  2.10271610e-02
			+ T * ( -6.28626133e-06 + T * (  7.85677565e-10
			+ T *  0.00000000e+00 ) ) ) ) +  5.14995970e+02 );
		cp[sNEOXC5H11] =  1.16876212e+02 * (
			 5.77200574e+00 + T * (  4.20543220e-02 
			+ T * ( -1.88587840e-05 + T * (  3.14271026e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIXC4H7O] =  1.16945257e+02 * (
			T * (  4.14496432e+00 + T * (  1.10818501e-02
			+ T * ( -3.07324644e-06 + T * (  3.84453565e-10
			+ T *  0.00000000e+00 ) ) ) ) +  2.20157645e+03 );
		cp[sIXC4H7O] =  1.16945257e+02 * (
			 4.14496432e+00 + T * (  2.21637003e-02 
			+ T * ( -9.21973931e-06 + T * (  1.53781426e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIXC4H8O] =  1.15310385e+02 * (
			T * (  6.04731897e+00 + T * (  1.64117506e-02
			+ T * ( -4.84914563e-06 + T * (  5.97657432e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.94905574e+04 );
		cp[sIXC4H8O] =  1.15310385e+02 * (
			 6.04731897e+00 + T * (  3.28235012e-02 
			+ T * ( -1.45474369e-05 + T * (  2.39062973e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIC4H8OOHXIO2XC4H9O4] =  6.86500099e+01 * (
			T * (  9.96114693e+00 + T * (  2.04295144e-02
			+ T * ( -6.20924330e-06 + T * (  7.82442710e-10
			+ T *  0.00000000e+00 ) ) ) ) -2.47969947e+04 );
		cp[sIC4H8OOHXIO2XC4H9O4] =  6.86500099e+01 * (
			 9.96114693e+00 + T * (  4.08590288e-02 
			+ T * ( -1.86277299e-05 + T * (  3.12977084e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIC4H8OOHXIXC4H9O2] =  9.33021366e+01 * (
			T * (  6.91173411e+00 + T * (  1.94716722e-02
			+ T * ( -6.03886310e-06 + T * (  7.73428673e-10
			+ T *  0.00000000e+00 ) ) ) ) -5.89539486e+03 );
		cp[sIC4H8OOHXIXC4H9O2] =  9.33021366e+01 * (
			 6.91173411e+00 + T * (  3.89433443e-02 
			+ T * ( -1.81165893e-05 + T * (  3.09371469e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIC4H8OOHXTO2XC4H9O4] =  6.86500099e+01 * (
			T * (  1.24802464e+01 + T * (  1.80460041e-02
			+ T * ( -5.23893583e-06 + T * (  6.36004885e-10
			+ T *  0.00000000e+00 ) ) ) ) -2.82744672e+04 );
		cp[sIC4H8OOHXTO2XC4H9O4] =  6.86500099e+01 * (
			 1.24802464e+01 + T * (  3.60920081e-02 
			+ T * ( -1.57168075e-05 + T * (  2.54401954e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIC4H8OOHXTXC4H9O2] =  9.33021366e+01 * (
			T * (  2.90219569e+00 + T * (  2.28751103e-02
			+ T * ( -7.36397300e-06 + T * (  9.69660102e-10
			+ T *  0.00000000e+00 ) ) ) ) -6.99529787e+03 );
		cp[sIC4H8OOHXTXC4H9O2] =  9.33021366e+01 * (
			 2.90219569e+00 + T * (  4.57502205e-02 
			+ T * ( -2.20919190e-05 + T * (  3.87864041e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIC4KETIIXC4H8O3] =  7.98657112e+01 * (
			T * (  7.48905109e+00 + T * (  2.00241579e-02
			+ T * ( -6.43434440e-06 + T * (  8.45643625e-10
			+ T *  0.00000000e+00 ) ) ) ) -4.02203586e+04 );
		cp[sIC4KETIIXC4H8O3] =  7.98657112e+01 * (
			 7.48905109e+00 + T * (  4.00483158e-02 
			+ T * ( -1.93030332e-05 + T * (  3.38257450e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sTC4H8OOHXIXC4H9O2] =  9.33021366e+01 * (
			T * (  9.75631688e+00 + T * (  1.63648944e-02
			+ T * ( -4.69674377e-06 + T * (  5.64285470e-10
			+ T *  0.00000000e+00 ) ) ) ) -9.37988130e+03 );
		cp[sTC4H8OOHXIXC4H9O2] =  9.33021366e+01 * (
			 9.75631688e+00 + T * (  3.27297888e-02 
			+ T * ( -1.40902313e-05 + T * (  2.25714188e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sTXC4H9O] =  1.13720593e+02 * (
			T * (  6.29676884e+00 + T * (  1.74474283e-02
			+ T * ( -5.10250397e-06 + T * (  6.23539407e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.53396915e+04 );
		cp[sTXC4H9O] =  1.13720593e+02 * (
			 6.29676884e+00 + T * (  3.48948567e-02 
			+ T * ( -1.53075119e-05 + T * (  2.49415763e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIXC4H9O2] =  9.33021366e+01 * (
			T * (  5.89612230e+00 + T * (  1.92046343e-02
			+ T * ( -5.65367667e-06 + T * (  6.94097068e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.31450848e+04 );
		cp[sIXC4H9O2] =  9.33021366e+01 * (
			 5.89612230e+00 + T * (  3.84092686e-02 
			+ T * ( -1.69610300e-05 + T * (  2.77638827e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sTXC4H9O2] =  9.33021366e+01 * (
			T * (  8.51208693e+00 + T * (  1.67335993e-02
			+ T * ( -4.64891537e-06 + T * (  5.42574382e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.66528502e+04 );
		cp[sTXC4H9O2] =  9.33021366e+01 * (
			 8.51208693e+00 + T * (  3.34671986e-02 
			+ T * ( -1.39467461e-05 + T * (  2.17029753e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sTXC3H6CHO] =  1.16945257e+02 * (
			T * (  3.91099601e+00 + T * (  1.68059067e-02
			+ T * ( -5.25798680e-06 + T * (  6.65921487e-10
			+ T *  0.00000000e+00 ) ) ) ) -9.38390232e+03 );
		cp[sTXC3H6CHO] =  1.16945257e+02 * (
			 3.91099601e+00 + T * (  3.36118134e-02 
			+ T * ( -1.57739604e-05 + T * (  2.66368595e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sTXC3H6OCHO] =  9.54617893e+01 * (
			T * (  6.53671659e+00 + T * (  1.73906247e-02
			+ T * ( -5.60683163e-06 + T * (  7.42276282e-10
			+ T *  0.00000000e+00 ) ) ) ) -2.33929526e+04 );
		cp[sTXC3H6OCHO] =  9.54617893e+01 * (
			 6.53671659e+00 + T * (  3.47812495e-02 
			+ T * ( -1.68204949e-05 + T * (  2.96910513e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIXC3H5CHO] =  1.18627154e+02 * (
			T * (  7.19597854e+00 + T * (  1.24978145e-02
			+ T * ( -3.68171107e-06 + T * (  4.52301075e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.72891601e+04 );
		cp[sIXC3H5CHO] =  1.18627154e+02 * (
			 7.19597854e+00 + T * (  2.49956291e-02 
			+ T * ( -1.10451332e-05 + T * (  1.80920430e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIXC3H7CHO] =  1.15310385e+02 * (
			T * (  3.80412072e+00 + T * (  1.80965247e-02
			+ T * ( -5.56004940e-06 + T * (  7.08189018e-10
			+ T *  0.00000000e+00 ) ) ) ) -2.86218337e+04 );
		cp[sIXC3H7CHO] =  1.15310385e+02 * (
			 3.80412072e+00 + T * (  3.61930495e-02 
			+ T * ( -1.66801482e-05 + T * (  2.83275607e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIXC4H7OOH] =  9.43696086e+01 * (
			T * (  8.26065337e+00 + T * (  1.65318550e-02
			+ T * ( -4.92782483e-06 + T * (  6.10472037e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.63753273e+04 );
		cp[sIXC4H7OOH] =  9.43696086e+01 * (
			 8.26065337e+00 + T * (  3.30637100e-02 
			+ T * ( -1.47834745e-05 + T * (  2.44188815e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIXC3H5CO] =  1.20358135e+02 * (
			T * (  7.57359694e+00 + T * (  1.07440472e-02
			+ T * ( -3.18437445e-06 + T * (  3.92260472e-10
			+ T *  0.00000000e+00 ) ) ) ) +  1.33531033e+03 );
		cp[sIXC3H5CO] =  1.20358135e+02 * (
			 7.57359694e+00 + T * (  2.14880944e-02 
			+ T * ( -9.55312336e-06 + T * (  1.56904189e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sYXC7H15] =  8.38223611e+01 * (
			T * (  1.32921118e+00 + T * (  3.37901245e-02
			+ T * ( -1.02602831e-05 + T * (  1.29187030e-09
			+ T *  0.00000000e+00 ) ) ) ) -6.14108024e+03 );
		cp[sYXC7H15] =  8.38223611e+01 * (
			 1.32921118e+00 + T * (  6.75802490e-02 
			+ T * ( -3.07808493e-05 + T * (  5.16748120e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sYXC7H14] =  8.46829358e+01 * (
			T * (  5.51949916e+00 + T * (  2.91994107e-02
			+ T * ( -8.66953550e-06 + T * (  1.07295509e-09
			+ T *  0.00000000e+00 ) ) ) ) -1.61420030e+04 );
		cp[sYXC7H14] =  8.46829358e+01 * (
			 5.51949916e+00 + T * (  5.83988214e-02 
			+ T * ( -2.60086065e-05 + T * (  4.29182037e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sXXC7H13] =  8.55613642e+01 * (
			T * (  6.53246968e+00 + T * (  2.64904270e-02
			+ T * ( -7.59724337e-06 + T * (  9.08533427e-10
			+ T *  0.00000000e+00 ) ) ) ) -5.02592988e+02 );
		cp[sXXC7H13] =  8.55613642e+01 * (
			 6.53246968e+00 + T * (  5.29808540e-02 
			+ T * ( -2.27917301e-05 + T * (  3.63413371e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIXC8H18] =  7.27897815e+01 * (
			T * (  8.11399812e+00 + T * (  3.59469073e-02
			+ T * ( -1.08572515e-05 + T * (  1.36609200e-09
			+ T *  0.00000000e+00 ) ) ) ) -3.29583815e+04 );
		cp[sIXC8H18] =  7.27897815e+01 * (
			 8.11399812e+00 + T * (  7.18938145e-02 
			+ T * ( -3.25717545e-05 + T * (  5.46436799e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sAXC8H17] =  7.34378533e+01 * (
			T * (  8.65029385e+00 + T * (  3.40202160e-02
			+ T * ( -1.02866931e-05 + T * (  1.29533016e-09
			+ T *  0.00000000e+00 ) ) ) ) -8.28527571e+03 );
		cp[sAXC8H17] =  7.34378533e+01 * (
			 8.65029385e+00 + T * (  6.80404319e-02 
			+ T * ( -3.08600794e-05 + T * (  5.18132064e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sBXC8H17] =  7.34378533e+01 * (
			T * (  6.98293343e+00 + T * (  3.54557074e-02
			+ T * ( -1.09409074e-05 + T * (  1.40465082e-09
			+ T *  0.00000000e+00 ) ) ) ) -9.15032631e+03 );
		cp[sBXC8H17] =  7.34378533e+01 * (
			 6.98293343e+00 + T * (  7.09114148e-02 
			+ T * ( -3.28227222e-05 + T * (  5.61860330e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCXC8H17] =  7.34378533e+01 * (
			T * (  5.37958758e+00 + T * (  3.63998835e-02
			+ T * ( -1.11696659e-05 + T * (  1.42412168e-09
			+ T *  0.00000000e+00 ) ) ) ) -1.03210414e+04 );
		cp[sCXC8H17] =  7.34378533e+01 * (
			 5.37958758e+00 + T * (  7.27997671e-02 
			+ T * ( -3.35089976e-05 + T * (  5.69648672e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sDXC8H17] =  7.34378533e+01 * (
			T * (  8.65029385e+00 + T * (  3.40202160e-02
			+ T * ( -1.02866931e-05 + T * (  1.29533016e-09
			+ T *  0.00000000e+00 ) ) ) ) -8.68789272e+03 );
		cp[sDXC8H17] =  7.34378533e+01 * (
			 8.65029385e+00 + T * (  6.80404319e-02 
			+ T * ( -3.08600794e-05 + T * (  5.18132064e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sAXC8H17O2] =  5.72549857e+01 * (
			T * (  1.15102096e+01 + T * (  3.52642283e-02
			+ T * ( -1.05720902e-05 + T * (  1.32059416e-09
			+ T *  0.00000000e+00 ) ) ) ) -2.71417088e+04 );
		cp[sAXC8H17O2] =  5.72549857e+01 * (
			 1.15102096e+01 + T * (  7.05284566e-02 
			+ T * ( -3.17162705e-05 + T * (  5.28237664e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sAXC8H17O2H] =  5.68602965e+01 * (
			T * (  1.23498105e+01 + T * (  3.69404708e-02
			+ T * ( -1.12779752e-05 + T * (  1.42926277e-09
			+ T *  0.00000000e+00 ) ) ) ) -4.46533044e+04 );
		cp[sAXC8H17O2H] =  5.68602965e+01 * (
			 1.23498105e+01 + T * (  7.38809415e-02 
			+ T * ( -3.38339255e-05 + T * (  5.71705110e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCXC8H17O2] =  5.72549857e+01 * (
			T * (  1.40288813e+01 + T * (  3.28815408e-02
			+ T * ( -9.60228217e-06 + T * (  1.17423962e-09
			+ T *  0.00000000e+00 ) ) ) ) -3.06207677e+04 );
		cp[sCXC8H17O2] =  5.72549857e+01 * (
			 1.40288813e+01 + T * (  6.57630815e-02 
			+ T * ( -2.88068465e-05 + T * (  4.69695847e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCXC8H17O2H] =  5.68602965e+01 * (
			T * (  1.48376064e+01 + T * (  3.45871570e-02
			+ T * ( -1.03191748e-05 + T * (  1.28445539e-09
			+ T *  0.00000000e+00 ) ) ) ) -4.81195410e+04 );
		cp[sCXC8H17O2H] =  5.68602965e+01 * (
			 1.48376064e+01 + T * (  6.91743141e-02 
			+ T * ( -3.09575245e-05 + T * (  5.13782157e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sAXC8H17O] =  6.43445084e+01 * (
			T * (  1.04791401e+01 + T * (  3.45965196e-02
			+ T * ( -1.04618950e-05 + T * (  1.31870540e-09
			+ T *  0.00000000e+00 ) ) ) ) -2.55317132e+04 );
		cp[sAXC8H17O] =  6.43445084e+01 * (
			 1.04791401e+01 + T * (  6.91930393e-02 
			+ T * ( -3.13856849e-05 + T * (  5.27482162e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCXC8H17O] =  6.43445084e+01 * (
			T * (  1.20991457e+01 + T * (  3.33418731e-02
			+ T * ( -9.95639703e-06 + T * (  1.24032280e-09
			+ T *  0.00000000e+00 ) ) ) ) -2.94126400e+04 );
		cp[sCXC8H17O] =  6.43445084e+01 * (
			 1.20991457e+01 + T * (  6.66837462e-02 
			+ T * ( -2.98691911e-05 + T * (  4.96129118e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sAC8H16OOHXAXC8H17O2] =  5.72549857e+01 * (
			T * (  1.28765889e+01 + T * (  3.50242149e-02
			+ T * ( -1.07096994e-05 + T * (  1.35854679e-09
			+ T *  0.00000000e+00 ) ) ) ) -2.03801280e+04 );
		cp[sAC8H16OOHXAXC8H17O2] =  5.72549857e+01 * (
			 1.28765889e+01 + T * (  7.00484298e-02 
			+ T * ( -3.21290982e-05 + T * (  5.43418715e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sAC8H16OOHXBXC8H17O2] =  5.72549857e+01 * (
			T * (  1.13274233e+01 + T * (  3.65604868e-02
			+ T * ( -1.15334440e-05 + T * (  1.50815428e-09
			+ T *  0.00000000e+00 ) ) ) ) -2.09293541e+04 );
		cp[sAC8H16OOHXBXC8H17O2] =  5.72549857e+01 * (
			 1.13274233e+01 + T * (  7.31209737e-02 
			+ T * ( -3.46003320e-05 + T * (  6.03261713e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sAC8H16OOHXCXC8H17O2] =  5.72549857e+01 * (
			T * (  9.68145324e+00 + T * (  3.73374700e-02
			+ T * ( -1.15677184e-05 + T * (  1.48377604e-09
			+ T *  0.00000000e+00 ) ) ) ) -2.19394580e+04 );
		cp[sAC8H16OOHXCXC8H17O2] =  5.72549857e+01 * (
			 9.68145324e+00 + T * (  7.46749400e-02 
			+ T * ( -3.47031551e-05 + T * (  5.93510415e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCC8H16OOHXAXC8H17O2] =  5.72549857e+01 * (
			T * (  1.53715418e+01 + T * (  3.26595881e-02
			+ T * ( -9.74415223e-06 + T * (  1.21244411e-09
			+ T *  0.00000000e+00 ) ) ) ) -2.38510335e+04 );
		cp[sCC8H16OOHXAXC8H17O2] =  5.72549857e+01 * (
			 1.53715418e+01 + T * (  6.53191761e-02 
			+ T * ( -2.92324567e-05 + T * (  4.84977645e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sAC8H16OOHXAO2XC8H17O4] =  4.69164184e+01 * (
			T * (  1.37592150e+01 + T * (  3.90705930e-02
			+ T * ( -1.23019228e-05 + T * (  1.59436672e-09
			+ T *  0.00000000e+00 ) ) ) ) -3.83968741e+04 );
		cp[sAC8H16OOHXAO2XC8H17O4] =  4.69164184e+01 * (
			 1.37592150e+01 + T * (  7.81411860e-02 
			+ T * ( -3.69057684e-05 + T * (  6.37746686e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sAC8H16OOHXBO2XC8H17O4] =  4.69164184e+01 * (
			T * (  1.50085341e+01 + T * (  3.80285230e-02
			+ T * ( -1.19008975e-05 + T * (  1.53564711e-09
			+ T *  0.00000000e+00 ) ) ) ) -4.09454874e+04 );
		cp[sAC8H16OOHXBO2XC8H17O4] =  4.69164184e+01 * (
			 1.50085341e+01 + T * (  7.60570460e-02 
			+ T * ( -3.57026925e-05 + T * (  6.14258843e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sAC8H16OOHXCO2XC8H17O4] =  4.69164184e+01 * (
			T * (  1.61942108e+01 + T * (  3.68961862e-02
			+ T * ( -1.14364193e-05 + T * (  1.46516074e-09
			+ T *  0.00000000e+00 ) ) ) ) -4.18648526e+04 );
		cp[sAC8H16OOHXCO2XC8H17O4] =  4.69164184e+01 * (
			 1.61942108e+01 + T * (  7.37923724e-02 
			+ T * ( -3.43092578e-05 + T * (  5.86064297e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCC8H16OOHXAO2XC8H17O4] =  4.69164184e+01 * (
			T * (  1.61942108e+01 + T * (  3.68961862e-02
			+ T * ( -1.14364193e-05 + T * (  1.46516074e-09
			+ T *  0.00000000e+00 ) ) ) ) -4.18648526e+04 );
		cp[sCC8H16OOHXAO2XC8H17O4] =  4.69164184e+01 * (
			 1.61942108e+01 + T * (  7.37923724e-02 
			+ T * ( -3.43092578e-05 + T * (  5.86064297e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sIC8KETABXC8H16O3] =  5.18971587e+01 * (
			T * (  1.46142912e+01 + T * (  3.44417520e-02
			+ T * ( -1.06282252e-05 + T * (  1.35703653e-09
			+ T *  0.00000000e+00 ) ) ) ) -5.76012840e+04 );
		cp[sIC8KETABXC8H16O3] =  5.18971587e+01 * (
			 1.46142912e+01 + T * (  6.88835040e-02 
			+ T * ( -3.18846757e-05 + T * (  5.42814612e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC6H5CH2] =  9.12400413e+01 * (
			T * (  1.40439800e+01 + T * (  1.17469365e-02
			+ T * ( -2.84584557e-06 + T * (  3.47271025e-10
			+ T * -1.67228840e-14 ) ) ) ) +  1.85642030e+04 );
		cp[sC6H5CH2] =  9.12400413e+01 * (
			 1.40439800e+01 + T * (  2.34938730e-02 
			+ T * ( -8.53753670e-06 + T * (  1.38908410e-09
			+ T * -8.36144200e-14 ) ) ) );
		h[sC6H5] =  1.07838392e+02 * (
			T * (  1.14314180e+01 + T * (  8.50952250e-03
			+ T * ( -1.94624137e-06 + T * (  2.20236717e-10
			+ T * -9.61008340e-15 ) ) ) ) +  3.39423480e+04 );
		cp[sC6H5] =  1.07838392e+02 * (
			 1.14314180e+01 + T * (  1.70190450e-02 
			+ T * ( -5.83872410e-06 + T * (  8.80946870e-10
			+ T * -4.80504170e-14 ) ) ) );
		h[sC6H6] =  1.06446715e+02 * (
			T * (  9.77760890e+00 + T * (  1.29849920e-02
			+ T * ( -4.20744533e-06 + T * (  7.46770125e-10
			+ T * -5.58391360e-14 ) ) ) ) +  3.83029950e+04 );
		cp[sC6H6] =  1.06446715e+02 * (
			 9.77760890e+00 + T * (  2.59699840e-02 
			+ T * ( -1.26223360e-05 + T * (  2.98708050e-09
			+ T * -2.79195680e-13 ) ) ) );
		h[sC6H5OH] =  8.83489183e+01 * (
			T * (  1.49120730e+01 + T * (  9.18906750e-03
			+ T * ( -2.06610427e-06 + T * (  2.29958053e-10
			+ T * -9.84191300e-15 ) ) ) ) -1.83751990e+04 );
		cp[sC6H5OH] =  8.83489183e+01 * (
			 1.49120730e+01 + T * (  1.83781350e-02 
			+ T * ( -6.19831280e-06 + T * (  9.19832210e-10
			+ T * -4.92095650e-14 ) ) ) );
		h[sC6H5O] =  8.93054780e+01 * (
			T * (  1.38339840e+01 + T * (  8.80920150e-03
			+ T * ( -2.02320857e-06 + T * (  2.29970432e-10
			+ T * -1.00898362e-14 ) ) ) ) -6.92125490e+02 );
		cp[sC6H5O] =  8.93054780e+01 * (
			 1.38339840e+01 + T * (  1.76184030e-02 
			+ T * ( -6.06962570e-06 + T * (  9.19881730e-10
			+ T * -5.04491810e-14 ) ) ) );
		h[sC6H5CHO] =  7.83499501e+01 * (
			T * (  3.82504540e+00 + T * (  2.19042381e-02
			+ T * ( -7.11241910e-06 + T * (  9.42002285e-10
			+ T *  0.00000000e+00 ) ) ) ) -7.68300125e+03 );
		cp[sC6H5CHO] =  7.83499501e+01 * (
			 3.82504540e+00 + T * (  4.38084761e-02 
			+ T * ( -2.13372573e-05 + T * (  3.76800914e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC6H5CO] =  7.91013224e+01 * (
			T * (  5.30322109e+00 + T * (  1.93744705e-02
			+ T * ( -6.35730250e-06 + T * (  8.49035522e-10
			+ T *  0.00000000e+00 ) ) ) ) +  9.55521205e+03 );
		cp[sC6H5CO] =  7.91013224e+01 * (
			 5.30322109e+00 + T * (  3.87489410e-02 
			+ T * ( -1.90719075e-05 + T * (  3.39614209e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC4H5XIXC4H5] =  1.56637905e+02 * (
			T * (  1.19977600e+01 + T * (  3.99529000e-03
			+ T * ( -2.69939100e-07 + T * ( -1.14218325e-10
			+ T *  1.72738200e-14 ) ) ) ) +  3.22849300e+04 );
		cp[sC4H5XIXC4H5] =  1.56637905e+02 * (
			 1.19977600e+01 + T * (  7.99058000e-03 
			+ T * ( -8.09817300e-07 + T * ( -4.56873300e-10
			+ T *  8.63691000e-14 ) ) ) );
		h[sC4H6] =  1.53718755e+02 * (
			T * (  8.04658300e+00 + T * (  8.24262500e-03
			+ T * ( -1.84074233e-06 + T * (  2.03089825e-10
			+ T * -8.59015600e-15 ) ) ) ) +  1.37013000e+04 );
		cp[sC4H6] =  1.53718755e+02 * (
			 8.04658300e+00 + T * (  1.64852500e-02 
			+ T * ( -5.52222700e-06 + T * (  8.12359300e-10
			+ T * -4.29507800e-14 ) ) ) );
		h[sC6H5CH2O] =  7.76127177e+01 * (
			T * (  4.73095622e+00 + T * (  2.23394748e-02
			+ T * ( -7.04831573e-06 + T * (  9.14245510e-10
			+ T *  0.00000000e+00 ) ) ) ) +  1.02920692e+04 );
		cp[sC6H5CH2O] =  7.76127177e+01 * (
			 4.73095622e+00 + T * (  4.46789497e-02 
			+ T * ( -2.11449472e-05 + T * (  3.65698204e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC6H5C2H5] =  7.83189525e+01 * (
			T * (  5.63102551e+00 + T * (  2.61394630e-02
			+ T * ( -8.08143927e-06 + T * (  1.03338710e-09
			+ T *  0.00000000e+00 ) ) ) ) -8.82227909e+02 );
		cp[sC6H5C2H5] =  7.83189525e+01 * (
			 5.63102551e+00 + T * (  5.22789261e-02 
			+ T * ( -2.42443178e-05 + T * (  4.13354840e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC6H5CH2OH] =  7.68892300e+01 * (
			T * (  5.03378845e+00 + T * (  2.31785659e-02
			+ T * ( -7.21369517e-06 + T * (  9.27658675e-10
			+ T *  0.00000000e+00 ) ) ) ) -1.59791417e+04 );
		cp[sC6H5CH2OH] =  7.68892300e+01 * (
			 5.03378845e+00 + T * (  4.63571318e-02 
			+ T * ( -2.16410855e-05 + T * (  3.71063470e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC5H5] =  1.27736058e+02 * (
			T * (  7.47439380e+00 + T * (  8.00636650e-03
			+ T * ( -2.16077160e-09 + T * ( -8.95492575e-10
			+ T *  1.84730142e-13 ) ) ) ) +  2.81916630e+04 );
		cp[sC5H5] =  1.27736058e+02 * (
			 7.47439380e+00 + T * (  1.60127330e-02 
			+ T * ( -6.48231480e-09 + T * ( -3.58197030e-09
			+ T *  9.23650710e-13 ) ) ) );
		h[sC6H4O2] =  7.69191059e+01 * (
			T * (  1.43886174e+01 + T * (  9.08121050e-03
			+ T * ( -2.23311559e-06 + T * (  2.75244700e-10
			+ T * -1.33474453e-14 ) ) ) ) -2.12444054e+04 );
		cp[sC6H4O2] =  7.69191059e+01 * (
			 1.43886174e+01 + T * (  1.81624210e-02 
			+ T * ( -6.69934678e-06 + T * (  1.10097880e-09
			+ T * -6.67372266e-14 ) ) ) );
		h[sC5H6] =  1.25788072e+02 * (
			T * (  9.97578480e+00 + T * (  9.45277150e-03
			+ T * ( -2.28038203e-06 + T * (  2.77483500e-10
			+ T * -1.33360472e-14 ) ) ) ) +  1.10816930e+04 );
		cp[sC5H6] =  1.25788072e+02 * (
			 9.97578480e+00 + T * (  1.89055430e-02 
			+ T * ( -6.84114610e-06 + T * (  1.10993400e-09
			+ T * -6.66802360e-14 ) ) ) );
		h[sC5H4O] =  1.03822832e+02 * (
			T * (  6.68707870e+00 + T * (  1.27579043e-02
			+ T * ( -3.99828073e-06 + T * (  5.11626438e-10
			+ T *  0.00000000e+00 ) ) ) ) +  8.80393103e+01 );
		cp[sC5H4O] =  1.03822832e+02 * (
			 6.68707870e+00 + T * (  2.55158085e-02 
			+ T * ( -1.19948422e-05 + T * (  2.04650575e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC5H5O] =  1.02532248e+02 * (
			T * (  4.82915206e+00 + T * (  1.60799429e-02
			+ T * ( -5.18694580e-06 + T * (  6.81266442e-10
			+ T *  0.00000000e+00 ) ) ) ) +  1.86243938e+04 );
		cp[sC5H5O] =  1.02532248e+02 * (
			 4.82915206e+00 + T * (  3.21598858e-02 
			+ T * ( -1.55608374e-05 + T * (  2.72506577e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC4H5XNXC4H5] =  1.56637905e+02 * (
			T * (  1.28659700e+01 + T * (  3.97168450e-03
			+ T * ( -2.87548867e-07 + T * ( -1.16390875e-10
			+ T *  1.79022620e-14 ) ) ) ) +  3.78355200e+04 );
		cp[sC4H5XNXC4H5] =  1.56637905e+02 * (
			 1.28659700e+01 + T * (  7.94336900e-03 
			+ T * ( -8.62646600e-07 + T * ( -4.65563500e-10
			+ T *  8.95113100e-14 ) ) ) );
		h[sCH2CHCHCO] =  1.22140381e+02 * (
			T * (  1.06577517e+01 + T * (  6.68661305e-03
			+ T * ( -1.68578436e-06 + T * (  1.74796677e-10
			+ T *  0.00000000e+00 ) ) ) ) -3.38617306e+03 );
		cp[sCH2CHCHCO] =  1.22140381e+02 * (
			 1.06577517e+01 + T * (  1.33732261e-02 
			+ T * ( -5.05735308e-06 + T * (  6.99186707e-10
			+ T *  0.00000000e+00 ) ) ) );
		h[sCHCHCHO] =  1.51021542e+02 * (
			T * ( -7.38411152e+01 + T * (  1.04360651e-01
			+ T * ( -5.19962990e-05 + T * (  8.99662125e-09
			+ T *  0.00000000e+00 ) ) ) ) +  2.26497218e+04 );
		cp[sCHCHCHO] =  1.51021542e+02 * (
			-7.38411152e+01 + T * (  2.08721303e-01 
			+ T * ( -1.55988897e-04 + T * (  3.59864850e-08
			+ T *  0.00000000e+00 ) ) ) );
		h[sHCCCO] =  1.56761944e+02 * (
			T * (  8.92529842e+00 + T * (  2.02539283e-03
			+ T * ( -6.56048020e-07 + T * (  1.08013303e-10
			+ T *  0.00000000e+00 ) ) ) ) -6.55268915e+02 );
		cp[sHCCCO] =  1.56761944e+02 * (
			 8.92529842e+00 + T * (  4.05078565e-03 
			+ T * ( -1.96814406e-06 + T * (  4.32053214e-10
			+ T *  0.00000000e+00 ) ) ) );
		h[sC2H5OH] =  1.80479726e+02 * (
			T * (  2.40390887e+00 + T * (  1.14894971e-02
			+ T * ( -3.34607570e-06 + T * (  4.10315382e-10
			+ T *  0.00000000e+00 ) ) ) ) -3.00556140e+04 );
		cp[sC2H5OH] =  1.80479726e+02 * (
			 2.40390887e+00 + T * (  2.29789942e-02 
			+ T * ( -1.00382271e-05 + T * (  1.64126153e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sC2H4OH] =  1.84517088e+02 * (
			T * (  4.53468150e+00 + T * (  8.60506330e-03
			+ T * ( -2.54582098e-06 + T * (  3.20857480e-10
			+ T *  0.00000000e+00 ) ) ) ) -6.39836604e+03 );
		cp[sC2H4OH] =  1.84517088e+02 * (
			 4.53468150e+00 + T * (  1.72101266e-02 
			+ T * ( -7.63746295e-06 + T * (  1.28342992e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH3CHOH] =  1.84517088e+02 * (
			T * (  1.30193879e+00 + T * (  1.14403872e-02
			+ T * ( -3.65066547e-06 + T * (  4.83084037e-10
			+ T *  0.00000000e+00 ) ) ) ) -5.99986135e+03 );
		cp[sCH3CHOH] =  1.84517088e+02 * (
			 1.30193879e+00 + T * (  2.28807743e-02 
			+ T * ( -1.09519964e-05 + T * (  1.93233615e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH3CH2O] =  1.84517088e+02 * (
			T * (  2.81457039e+00 + T * (  1.02691869e-02
			+ T * ( -3.16582511e-06 + T * (  4.10898748e-10
			+ T *  0.00000000e+00 ) ) ) ) -3.95403690e+03 );
		cp[sCH3CH2O] =  1.84517088e+02 * (
			 2.81457039e+00 + T * (  2.05383738e-02 
			+ T * ( -9.49747532e-06 + T * (  1.64359499e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sCH3HCO] =  1.88739217e+02 * (
			T * (  5.86865000e+00 + T * (  5.39712000e-03
			+ T * ( -1.21517667e-06 + T * (  1.35322800e-10
			+ T * -5.79368800e-15 ) ) ) ) -2.26456900e+04 );
		cp[sCH3HCO] =  1.88739217e+02 * (
			 5.86865000e+00 + T * (  1.07942400e-02 
			+ T * ( -3.64553000e-06 + T * (  5.41291200e-10
			+ T * -2.89684400e-14 ) ) ) );
		h[sCH2HCO] =  1.93159093e+02 * (
			T * (  5.97567000e+00 + T * (  4.06529550e-03
			+ T * ( -9.14541333e-07 + T * (  1.01757600e-10
			+ T * -4.35203400e-15 ) ) ) ) +  4.90321800e+02 );
		cp[sCH2HCO] =  1.93159093e+02 * (
			 5.97567000e+00 + T * (  8.13059100e-03 
			+ T * ( -2.74362400e-06 + T * (  4.07030400e-10
			+ T * -2.17601700e-14 ) ) ) );
		h[sHOC2H4O2] =  1.07894368e+02 * (
			T * (  7.28358795e+00 + T * (  8.24244810e-03
			+ T * ( -1.80079794e-06 + T * (  1.56175649e-10
			+ T *  0.00000000e+00 ) ) ) ) -2.25435240e+04 );
		cp[sHOC2H4O2] =  1.07894368e+02 * (
			 7.28358795e+00 + T * (  1.64848962e-02 
			+ T * ( -5.40239381e-06 + T * (  6.24702597e-10
			+ T *  0.00000000e+00 ) ) ) );
		h[sC2H] =  3.32201534e+02 * (
			T * (  3.98636700e+00 + T * (  1.57156150e-03
			+ T * ( -4.22414333e-07 + T * (  7.31090750e-11
			+ T * -5.43264000e-15 ) ) ) ) +  6.65588400e+04 );
		cp[sC2H] =  3.32201534e+02 * (
			 3.98636700e+00 + T * (  3.14312300e-03 
			+ T * ( -1.26724300e-06 + T * (  2.92436300e-10
			+ T * -2.71632000e-14 ) ) ) );
		h[sH2CC] =  3.19340144e+02 * (
			T * (  4.27803400e+00 + T * (  2.37814020e-03
			+ T * ( -5.43366967e-07 + T * (  6.36570150e-11
			+ T * -2.97727580e-15 ) ) ) ) +  4.83166880e+04 );
		cp[sH2CC] =  3.19340144e+02 * (
			 4.27803400e+00 + T * (  4.75628040e-03 
			+ T * ( -1.63010090e-06 + T * (  2.54628060e-10
			+ T * -1.48863790e-14 ) ) ) );
		h[sC3H3] =  2.12893430e+02 * (
			T * (  8.83104700e+00 + T * (  2.17859700e-03
			+ T * ( -1.36968867e-07 + T * ( -5.92180750e-11
			+ T *  8.75304000e-15 ) ) ) ) +  3.76322300e+04 );
		cp[sC3H3] =  2.12893430e+02 * (
			 8.83104700e+00 + T * (  4.35719400e-03 
			+ T * ( -4.10906600e-07 + T * ( -2.36872300e-10
			+ T *  4.37652000e-14 ) ) ) );
		h[sPXC3H4] =  2.07536818e+02 * (
			T * (  1.76237742e+00 + T * (  1.02784310e-02
			+ T * ( -3.37305433e-06 + T * (  4.59840425e-10
			+ T *  0.00000000e+00 ) ) ) ) +  2.17000925e+04 );
		cp[sPXC3H4] =  2.07536818e+02 * (
			 1.76237742e+00 + T * (  2.05568619e-02 
			+ T * ( -1.01191630e-05 + T * (  1.83936170e-09
			+ T *  0.00000000e+00 ) ) ) );
		h[sAXC3H4] =  2.07536818e+02 * (
			T * (  5.72914400e+00 + T * (  6.18402250e-03
			+ T * ( -1.60187533e-06 + T * (  2.15034100e-10
			+ T * -1.16256040e-14 ) ) ) ) +  2.01298400e+04 );
		cp[sAXC3H4] =  2.07536818e+02 * (
			 5.72914400e+00 + T * (  1.23680450e-02 
			+ T * ( -4.80562600e-06 + T * (  8.60136400e-10
			+ T * -5.81280200e-14 ) ) ) );
		h[s1XCH2] =  5.92780550e+02 * (
			T * (  2.29203842e+00 + T * (  2.32794318e-03
			+ T * ( -6.70639823e-07 + T * (  1.04476500e-10
			+ T * -6.79432730e-15 ) ) ) ) +  5.09259997e+04 );
		cp[s1XCH2] =  5.92780550e+02 * (
			 2.29203842e+00 + T * (  4.65588637e-03 
			+ T * ( -2.01191947e-06 + T * (  4.17906000e-10
			+ T * -3.39716365e-14 ) ) ) );
		h[sC] =  6.92284763e+02 * (
			T * (  2.49266888e+00 + T * (  2.39944642e-05
			+ T * ( -2.41445007e-08 + T * (  9.35727573e-12
			+ T * -9.74555786e-16 ) ) ) ) +  8.54512953e+04 );
		cp[sC] =  6.92284763e+02 * (
			 2.49266888e+00 + T * (  4.79889284e-05 
			+ T * ( -7.24335020e-08 + T * (  3.74291029e-11
			+ T * -4.87277893e-15 ) ) ) );
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
		h[sC6H5CH3] =  9.02418217e+01 * (
			T * (  1.61526630e+00 + T * (  1.05497190e-02
			+ T * (  2.84553393e-05 + T * ( -3.31526650e-08
			+ T *  1.11913208e-11 ) ) ) ) +  4.07563000e+03 );
		cp[sC6H5CH3] =  9.02418217e+01 * (
			 1.61526630e+00 + T * (  2.10994380e-02 
			+ T * (  8.53660180e-05 + T * ( -1.32610660e-07
			+ T *  5.59566040e-11 ) ) ) );
		h[sHO2] =  2.51888633e+02 * (
			T * (  3.18310656e+00 + T * (  1.83383975e-03
			+ T * ( -3.10795041e-07 + T * ( -8.14632298e-11
			+ T *  3.02279824e-14 ) ) ) ) +  8.09181013e+02 );
		cp[sHO2] =  2.51888633e+02 * (
			 3.18310656e+00 + T * (  3.66767950e-03 
			+ T * ( -9.32385122e-07 + T * ( -3.25852919e-10
			+ T *  1.51139912e-13 ) ) ) );
		h[sOH] =  4.88848777e+02 * (
			T * (  3.43586219e+00 + T * (  1.01117902e-04
			+ T * ( -3.78488040e-08 + T * (  6.06112872e-11
			+ T * -1.48730206e-14 ) ) ) ) +  3.74321252e+03 );
		cp[sOH] =  4.88848777e+02 * (
			 3.43586219e+00 + T * (  2.02235804e-04 
			+ T * ( -1.13546412e-07 + T * (  2.42445149e-10
			+ T * -7.43651031e-14 ) ) ) );
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
		h[sC2H4] =  2.96390275e+02 * (
			T * ( -8.61488000e-01 + T * (  1.39808150e-02
			+ T * ( -1.12955900e-05 + T * (  6.96288000e-09
			+ T * -1.94757580e-12 ) ) ) ) +  5.57304600e+03 );
		cp[sC2H4] =  2.96390275e+02 * (
			-8.61488000e-01 + T * (  2.79616300e-02 
			+ T * ( -3.38867700e-05 + T * (  2.78515200e-08
			+ T * -9.73787900e-12 ) ) ) );
		h[sH] =  8.24835317e+03 * (
			T * (  2.50104422e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  2.54747466e+04 );
		cp[sH] =  8.24835317e+03 * (
			 2.50104422e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sC2H3] =  3.07437509e+02 * (
			T * (  2.45927600e+00 + T * (  3.68573800e-03
			+ T * (  7.03291000e-07 + T * ( -3.30410500e-10
			+ T * -2.36956800e-13 ) ) ) ) +  3.33522500e+04 );
		cp[sC2H3] =  3.07437509e+02 * (
			 2.45927600e+00 + T * (  7.37147600e-03 
			+ T * (  2.10987300e-06 + T * ( -1.32164200e-09
			+ T * -1.18478400e-12 ) ) ) );
		h[sH2] =  4.12417659e+03 * (
			T * (  3.29812400e+00 + T * (  4.12472100e-04
			+ T * ( -2.71433833e-07 + T * ( -2.36885850e-11
			+ T *  8.26974400e-14 ) ) ) ) -1.01252100e+03 );
		cp[sH2] =  4.12417659e+03 * (
			 3.29812400e+00 + T * (  8.24944200e-04 
			+ T * ( -8.14301500e-07 + T * ( -9.47543400e-11
			+ T *  4.13487200e-13 ) ) ) );
		h[sCH2O] =  2.76904683e+02 * (
			T * (  1.65273100e+00 + T * (  6.31572000e-03
			+ T * ( -6.29389333e-06 + T * (  5.12507750e-09
			+ T * -1.68264740e-12 ) ) ) ) -1.48654000e+04 );
		cp[sCH2O] =  2.76904683e+02 * (
			 1.65273100e+00 + T * (  1.26314400e-02 
			+ T * ( -1.88816800e-05 + T * (  2.05003100e-08
			+ T * -8.41323700e-12 ) ) ) );
		h[sHCO] =  2.86523537e+02 * (
			T * (  2.89833000e+00 + T * (  3.09957350e-03
			+ T * ( -3.20769467e-06 + T * (  2.72456250e-09
			+ T * -9.14977000e-13 ) ) ) ) +  4.15992200e+03 );
		cp[sHCO] =  2.86523537e+02 * (
			 2.89833000e+00 + T * (  6.19914700e-03 
			+ T * ( -9.62308400e-06 + T * (  1.08982500e-08
			+ T * -4.57488500e-12 ) ) ) );
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
		h[sCH3O2] =  1.76772973e+02 * (
			T * (  4.26146906e+00 + T * (  5.04367995e-03
			+ T * ( -1.07168728e-06 + T * (  5.23523168e-11
			+ T *  8.36678206e-15 ) ) ) ) -6.84394259e+02 );
		cp[sCH3O2] =  1.76772973e+02 * (
			 4.26146906e+00 + T * (  1.00873599e-02 
			+ T * ( -3.21506184e-06 + T * (  2.09409267e-10
			+ T *  4.18339103e-14 ) ) ) );
		h[sCH3] =  5.53035786e+02 * (
			T * (  2.43044300e+00 + T * (  5.56205000e-03
			+ T * ( -5.60073333e-06 + T * (  4.05457250e-09
			+ T * -1.17299060e-12 ) ) ) ) +  1.64237800e+04 );
		cp[sCH3] =  5.53035786e+02 * (
			 2.43044300e+00 + T * (  1.11241000e-02 
			+ T * ( -1.68022000e-05 + T * (  1.62182900e-08
			+ T * -5.86495300e-12 ) ) ) );
		h[sCH3O] =  2.67910679e+02 * (
			T * (  2.10620400e+00 + T * (  3.60829750e-03
			+ T * (  1.77949067e-06 + T * ( -1.84440900e-09
			+ T *  4.15122200e-13 ) ) ) ) +  9.78601100e+02 );
		cp[sCH3O] =  2.67910679e+02 * (
			 2.10620400e+00 + T * (  7.21659500e-03 
			+ T * (  5.33847200e-06 + T * ( -7.37763600e-09
			+ T *  2.07561100e-12 ) ) ) );
		h[sCH3O2H] =  1.73063986e+02 * (
			T * (  3.23442817e+00 + T * (  9.50648835e-03
			+ T * ( -3.77954290e-06 + T * (  8.50766632e-10
			+ T * -8.23660444e-14 ) ) ) ) -1.77197926e+04 );
		cp[sCH3O2H] =  1.73063986e+02 * (
			 3.23442817e+00 + T * (  1.90129767e-02 
			+ T * ( -1.13386287e-05 + T * (  3.40306653e-09
			+ T * -4.11830222e-13 ) ) ) );
		h[sCH4] =  5.18285750e+02 * (
			T * (  7.78741500e-01 + T * (  8.73834000e-03
			+ T * ( -9.27803000e-06 + T * (  7.62427000e-09
			+ T * -2.44786200e-12 ) ) ) ) -9.82522900e+03 );
		cp[sCH4] =  5.18285750e+02 * (
			 7.78741500e-01 + T * (  1.74766800e-02 
			+ T * ( -2.78340900e-05 + T * (  3.04970800e-08
			+ T * -1.22393100e-11 ) ) ) );
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
		h[sO] =  5.19646250e+02 * (
			T * (  2.94642900e+00 + T * ( -8.19083000e-04
			+ T * (  8.07010667e-07 + T * ( -4.00710750e-10
			+ T *  7.78139200e-14 ) ) ) ) +  2.91476400e+04 );
		cp[sO] =  5.19646250e+02 * (
			 2.94642900e+00 + T * ( -1.63816600e-03 
			+ T * (  2.42103200e-06 + T * ( -1.60284300e-09
			+ T *  3.89069600e-13 ) ) ) );
		h[sH2O2] =  2.44424389e+02 * (
			T * (  3.38875400e+00 + T * (  3.28461300e-03
			+ T * ( -4.95004333e-08 + T * ( -1.15645150e-09
			+ T *  4.94303000e-13 ) ) ) ) -1.76631500e+04 );
		cp[sH2O2] =  2.44424389e+02 * (
			 3.38875400e+00 + T * (  6.56922600e-03 
			+ T * ( -1.48501300e-07 + T * ( -4.62580600e-09
			+ T *  2.47151500e-12 ) ) ) );
		h[sC2H2] =  3.19340144e+02 * (
			T * (  2.01356200e+00 + T * (  7.59522500e-03
			+ T * ( -5.38773000e-06 + T * (  2.26974800e-09
			+ T * -3.82549200e-13 ) ) ) ) +  2.61244400e+04 );
		cp[sC2H2] =  3.19340144e+02 * (
			 2.01356200e+00 + T * (  1.51904500e-02 
			+ T * ( -1.61631900e-05 + T * (  9.07899200e-09
			+ T * -1.91274600e-12 ) ) ) );
		h[sCH2] =  5.92780550e+02 * (
			T * (  3.76223700e+00 + T * (  5.79909500e-04
			+ T * (  8.29861667e-08 + T * (  2.20020900e-10
			+ T * -1.46648700e-13 ) ) ) ) +  4.53679100e+04 );
		cp[sCH2] =  5.92780550e+02 * (
			 3.76223700e+00 + T * (  1.15981900e-03 
			+ T * (  2.48958500e-07 + T * (  8.80083600e-10
			+ T * -7.33243500e-13 ) ) ) );
		h[sCH2OH] =  2.67910679e+02 * (
			T * (  2.86262800e+00 + T * (  5.00763500e-03
			+ T * ( -1.76181200e-07 + T * ( -1.28463500e-09
			+ T *  4.49208200e-13 ) ) ) ) -3.34967900e+03 );
		cp[sCH2OH] =  2.67910679e+02 * (
			 2.86262800e+00 + T * (  1.00152700e-02 
			+ T * ( -5.28543600e-07 + T * ( -5.13854000e-09
			+ T *  2.24604100e-12 ) ) ) );
		h[sHCCO] =  2.02650385e+02 * (
			T * (  5.04796500e+00 + T * (  2.22673900e-03
			+ T * (  7.56094333e-08 + T * ( -3.70523750e-10
			+ T *  4.50148400e-14 ) ) ) ) +  1.96589200e+04 );
		cp[sHCCO] =  2.02650385e+02 * (
			 5.04796500e+00 + T * (  4.45347800e-03 
			+ T * (  2.26828300e-07 + T * ( -1.48209500e-09
			+ T *  2.25074200e-13 ) ) ) );
		h[sCH2CO] =  1.97790941e+02 * (
			T * (  2.97497100e+00 + T * (  6.05935500e-03
			+ T * ( -7.81682000e-07 + T * ( -1.61667125e-09
			+ T *  7.81129800e-13 ) ) ) ) -7.63263700e+03 );
		cp[sCH2CO] =  1.97790941e+02 * (
			 2.97497100e+00 + T * (  1.21187100e-02 
			+ T * ( -2.34504600e-06 + T * ( -6.46668500e-09
			+ T *  3.90564900e-12 ) ) ) );
		h[sCH3CO] =  1.93159093e+02 * (
			T * (  3.12527800e+00 + T * (  4.88911000e-03
			+ T * (  1.50714933e-06 + T * ( -2.25236550e-09
			+ T *  6.38743600e-13 ) ) ) ) -4.10850800e+03 );
		cp[sCH3CO] =  1.93159093e+02 * (
			 3.12527800e+00 + T * (  9.77822000e-03 
			+ T * (  4.52144800e-06 + T * ( -9.00946200e-09
			+ T *  3.19371800e-12 ) ) ) );
		h[sAXC3H5] =  2.02443146e+02 * (
			T * (  2.27648600e+00 + T * (  9.92782000e-03
			+ T * (  3.74614000e-07 + T * ( -2.53644000e-09
			+ T *  6.88268400e-13 ) ) ) ) +  1.78949700e+04 );
		cp[sAXC3H5] =  2.02443146e+02 * (
			 2.27648600e+00 + T * (  1.98556400e-02 
			+ T * (  1.12384200e-06 + T * ( -1.01457600e-08
			+ T *  3.44134200e-12 ) ) ) );
		h[sC3H6] =  1.97593517e+02 * (
			T * (  3.94615444e-01 + T * (  1.44553831e-02
			+ T * ( -5.16289360e-06 + T * (  9.72035522e-10
			+ T * -6.75780704e-14 ) ) ) ) +  1.06688164e+03 );
		cp[sC3H6] =  1.97593517e+02 * (
			 3.94615444e-01 + T * (  2.89107662e-02 
			+ T * ( -1.54886808e-05 + T * (  3.88814209e-09
			+ T * -3.37890352e-13 ) ) ) );
		h[sC3H5O] =  1.45686701e+02 * (
			T * (  1.19822582e+00 + T * (  1.52789918e-02
			+ T * ( -6.02100920e-06 + T * (  1.21537508e-09
			+ T * -8.39709124e-14 ) ) ) ) +  9.58217784e+03 );
		cp[sC3H5O] =  1.45686701e+02 * (
			 1.19822582e+00 + T * (  3.05579837e-02 
			+ T * ( -1.80630276e-05 + T * (  4.86150033e-09
			+ T * -4.19854562e-13 ) ) ) );
		h[sCH3OH] =  2.59482554e+02 * (
			T * (  2.66011500e+00 + T * (  3.67075400e-03
			+ T * (  2.39001700e-06 + T * ( -2.19829850e-09
			+ T *  4.78114000e-13 ) ) ) ) -2.53534800e+04 );
		cp[sCH3OH] =  2.59482554e+02 * (
			 2.66011500e+00 + T * (  7.34150800e-03 
			+ T * (  7.17005100e-06 + T * ( -8.79319400e-09
			+ T *  2.39057000e-12 ) ) ) );
		h[sCH3CHO] =  1.88739217e+02 * (
			T * (  2.50569500e+00 + T * (  6.68495500e-03
			+ T * (  1.55731767e-06 + T * ( -2.82035000e-09
			+ T *  8.52713200e-13 ) ) ) ) -2.12458900e+04 );
		cp[sCH3CHO] =  1.88739217e+02 * (
			 2.50569500e+00 + T * (  1.33699100e-02 
			+ T * (  4.67195300e-06 + T * ( -1.12814000e-08
			+ T *  4.26356600e-12 ) ) ) );
		h[sNXC3H7] =  1.92970803e+02 * (
			T * (  1.92253700e+00 + T * (  1.23946350e-02
			+ T * (  6.03416333e-07 + T * ( -4.45816500e-09
			+ T *  1.71659920e-12 ) ) ) ) +  9.71328100e+03 );
		cp[sNXC3H7] =  1.92970803e+02 * (
			 1.92253700e+00 + T * (  2.47892700e-02 
			+ T * (  1.81024900e-06 + T * ( -1.78326600e-08
			+ T *  8.58299600e-12 ) ) ) );
		h[sC5H9] =  1.20285003e+02 * (
			T * ( -1.38013950e+00 + T * (  2.78804243e-02
			+ T * ( -1.23381309e-05 + T * (  3.17209752e-09
			+ T * -3.57077670e-13 ) ) ) ) +  1.25589824e+04 );
		cp[sC5H9] =  1.20285003e+02 * (
			-1.38013950e+00 + T * (  5.57608487e-02 
			+ T * ( -3.70143928e-05 + T * (  1.26883901e-08
			+ T * -1.78538835e-12 ) ) ) );
		h[sPXC4H9] =  1.45579563e+02 * (
			T * ( -4.37779725e-01 + T * (  2.39486182e-02
			+ T * ( -1.04674386e-05 + T * (  2.74466180e-09
			+ T * -3.24021328e-13 ) ) ) ) +  7.68945248e+03 );
		cp[sPXC4H9] =  1.45579563e+02 * (
			-4.37779725e-01 + T * (  4.78972364e-02 
			+ T * ( -3.14023159e-05 + T * (  1.09786472e-08
			+ T * -1.62010664e-12 ) ) ) );
		h[s1XC5H10] =  1.18556110e+02 * (
			T * ( -1.06223481e+00 + T * (  2.87109147e-02
			+ T * ( -1.24828963e-05 + T * (  3.18412472e-09
			+ T * -3.59219578e-13 ) ) ) ) -4.46546666e+03 );
		cp[s1XC5H10] =  1.18556110e+02 * (
			-1.06223481e+00 + T * (  5.74218294e-02 
			+ T * ( -3.74486890e-05 + T * (  1.27364989e-08
			+ T * -1.79609789e-12 ) ) ) );
		h[s1XC5H11] =  1.16876212e+02 * (
			T * ( -9.05255912e-01 + T * (  3.05316426e-02
			+ T * ( -1.36497275e-05 + T * (  3.65233675e-09
			+ T * -4.37719230e-13 ) ) ) ) +  4.83995303e+03 );
		cp[s1XC5H11] =  1.16876212e+02 * (
			-9.05255912e-01 + T * (  6.10632852e-02 
			+ T * ( -4.09491825e-05 + T * (  1.46093470e-08
			+ T * -2.18859615e-12 ) ) ) );
		h[sCH3COCH2O2] =  9.33461323e+01 * (
			T * (  5.95535468e+00 + T * (  1.35127603e-02
			+ T * ( -4.57950103e-06 + T * (  8.84339628e-10
			+ T * -8.07845114e-14 ) ) ) ) -2.06679464e+04 );
		cp[sCH3COCH2O2] =  9.33461323e+01 * (
			 5.95535468e+00 + T * (  2.70255205e-02 
			+ T * ( -1.37385031e-05 + T * (  3.53735851e-09
			+ T * -4.03922557e-13 ) ) ) );
		h[sCH3COCH2] =  1.45686701e+02 * (
			T * (  1.80339187e+00 + T * (  1.50703542e-02
			+ T * ( -6.45018507e-06 + T * (  1.59549759e-09
			+ T * -1.73220636e-13 ) ) ) ) -5.37233261e+03 );
		cp[sCH3COCH2] =  1.45686701e+02 * (
			 1.80339187e+00 + T * (  3.01407085e-02 
			+ T * ( -1.93505552e-05 + T * (  6.38199034e-09
			+ T * -8.66103180e-13 ) ) ) );
		h[sCH3COCH2O2H] =  9.23015609e+01 * (
			T * (  4.94789761e+00 + T * (  1.80237216e-02
			+ T * ( -7.39066443e-06 + T * (  1.74574218e-09
			+ T * -1.84253852e-13 ) ) ) ) -3.88687178e+04 );
		cp[sCH3COCH2O2H] =  9.23015609e+01 * (
			 4.94789761e+00 + T * (  3.60474432e-02 
			+ T * ( -2.21719933e-05 + T * (  6.98296874e-09
			+ T * -9.21269260e-13 ) ) ) );
		h[sCH3COCH2O] =  1.13785959e+02 * (
			T * (  4.18402754e+00 + T * (  1.86951530e-02
			+ T * ( -7.93248273e-06 + T * (  1.46425396e-09
			+ T * -5.16760300e-14 ) ) ) ) -1.25344087e+04 );
		cp[sCH3COCH2O] =  1.13785959e+02 * (
			 4.18402754e+00 + T * (  3.73903061e-02 
			+ T * ( -2.37974482e-05 + T * (  5.85701585e-09
			+ T * -2.58380150e-13 ) ) ) );
		h[sCH2CHO] =  1.93159093e+02 * (
			T * (  3.40906200e+00 + T * (  5.36928500e-03
			+ T * (  6.30497333e-07 + T * ( -1.78964575e-09
			+ T *  5.73477000e-13 ) ) ) ) +  1.52147700e+03 );
		cp[sCH2CHO] =  1.93159093e+02 * (
			 3.40906200e+00 + T * (  1.07385700e-02 
			+ T * (  1.89149200e-06 + T * ( -7.15858300e-09
			+ T *  2.86738500e-12 ) ) ) );
		h[sNXC7H16] =  8.29791014e+01 * (
			T * ( -1.26836187e+00 + T * (  4.27177910e-02
			+ T * ( -1.75115595e-05 + T * (  4.07364302e-09
			+ T * -4.04789850e-13 ) ) ) ) -2.56586565e+04 );
		cp[sNXC7H16] =  8.29791014e+01 * (
			-1.26836187e+00 + T * (  8.54355820e-02 
			+ T * ( -5.25346786e-05 + T * (  1.62945721e-08
			+ T * -2.02394925e-12 ) ) ) );
		h[s1XC7H15] =  8.38223611e+01 * (
			T * ( -4.99570406e-01 + T * (  4.04413234e-02
			+ T * ( -1.66844251e-05 + T * (  3.91373270e-09
			+ T * -3.93232454e-13 ) ) ) ) -1.04590223e+03 );
		cp[s1XC7H15] =  8.38223611e+01 * (
			-4.99570406e-01 + T * (  8.08826467e-02 
			+ T * ( -5.00532754e-05 + T * (  1.56549308e-08
			+ T * -1.96616227e-12 ) ) ) );
		h[s2XC7H15] =  8.38223611e+01 * (
			T * ( -3.79155767e-02 + T * (  3.78363285e-02
			+ T * ( -1.35824545e-05 + T * (  2.33169736e-09
			+ T * -9.84721490e-14 ) ) ) ) -2.35605303e+03 );
		cp[s2XC7H15] =  8.38223611e+01 * (
			-3.79155767e-02 + T * (  7.56726570e-02 
			+ T * ( -4.07473634e-05 + T * (  9.32678943e-09
			+ T * -4.92360745e-13 ) ) ) );
		h[s3XC7H15] =  8.38223611e+01 * (
			T * ( -3.79155767e-02 + T * (  3.78363285e-02
			+ T * ( -1.35824545e-05 + T * (  2.33169736e-09
			+ T * -9.84721490e-14 ) ) ) ) -2.35605303e+03 );
		cp[s3XC7H15] =  8.38223611e+01 * (
			-3.79155767e-02 + T * (  7.56726570e-02 
			+ T * ( -4.07473634e-05 + T * (  9.32678943e-09
			+ T * -4.92360745e-13 ) ) ) );
		h[s4XC7H15] =  8.38223611e+01 * (
			T * ( -3.79155767e-02 + T * (  3.78363285e-02
			+ T * ( -1.35824545e-05 + T * (  2.33169736e-09
			+ T * -9.84721490e-14 ) ) ) ) -2.35605303e+03 );
		cp[s4XC7H15] =  8.38223611e+01 * (
			-3.79155767e-02 + T * (  7.56726570e-02 
			+ T * ( -4.07473634e-05 + T * (  9.32678943e-09
			+ T * -4.92360745e-13 ) ) ) );
		h[s1XC7H14] =  8.46829358e+01 * (
			T * ( -1.67720549e+00 + T * (  4.12305800e-02
			+ T * ( -1.82168036e-05 + T * (  4.69655757e-09
			+ T * -5.31475966e-13 ) ) ) ) -1.02168601e+04 );
		cp[s1XC7H14] =  8.46829358e+01 * (
			-1.67720549e+00 + T * (  8.24611601e-02 
			+ T * ( -5.46504108e-05 + T * (  1.87862303e-08
			+ T * -2.65737983e-12 ) ) ) );
		h[s2XC7H15O2] =  6.33763244e+01 * (
			T * (  1.51378168e+00 + T * (  4.42786373e-02
			+ T * ( -1.97485716e-05 + T * (  5.29504655e-09
			+ T * -6.41483444e-13 ) ) ) ) -2.19818400e+04 );
		cp[s2XC7H15O2] =  6.33763244e+01 * (
			 1.51378168e+00 + T * (  8.85572745e-02 
			+ T * ( -5.92457147e-05 + T * (  2.11801862e-08
			+ T * -3.20741722e-12 ) ) ) );
		h[s2X3XC7H14O2H] =  6.33763244e+01 * (
			T * (  1.62083964e+00 + T * (  4.41818108e-02
			+ T * ( -1.88175745e-05 + T * (  4.40814997e-09
			+ T * -4.20738684e-13 ) ) ) ) -1.57072842e+04 );
		cp[s2X3XC7H14O2H] =  6.33763244e+01 * (
			 1.62083964e+00 + T * (  8.83636215e-02 
			+ T * ( -5.64527235e-05 + T * (  1.76325999e-08
			+ T * -2.10369342e-12 ) ) ) );
		h[s2X4XC7H14O2H] =  6.33763244e+01 * (
			T * (  1.62083964e+00 + T * (  4.41818108e-02
			+ T * ( -1.88175745e-05 + T * (  4.40814997e-09
			+ T * -4.20738684e-13 ) ) ) ) -1.57072842e+04 );
		cp[s2X4XC7H14O2H] =  6.33763244e+01 * (
			 1.62083964e+00 + T * (  8.83636215e-02 
			+ T * ( -5.64527235e-05 + T * (  1.76325999e-08
			+ T * -2.10369342e-12 ) ) ) );
		h[s2X4XC7H14OOHO2] =  5.09488326e+01 * (
			T * (  3.09603698e+00 + T * (  4.96591120e-02
			+ T * ( -2.23015309e-05 + T * (  5.60478175e-09
			+ T * -5.95687324e-13 ) ) ) ) -3.52824732e+04 );
		cp[s2X4XC7H14OOHO2] =  5.09488326e+01 * (
			 3.09603698e+00 + T * (  9.93182239e-02 
			+ T * ( -6.69045928e-05 + T * (  2.24191270e-08
			+ T * -2.97843662e-12 ) ) ) );
		h[sNC7KET24XC7H14O3] =  5.68766332e+01 * (
			T * (  2.69842062e+00 + T * (  4.50729276e-02
			+ T * ( -2.06471129e-05 + T * (  5.56972598e-09
			+ T * -6.68357896e-13 ) ) ) ) -5.27915898e+04 );
		cp[sNC7KET24XC7H14O3] =  5.68766332e+01 * (
			 2.69842062e+00 + T * (  9.01458551e-02 
			+ T * ( -6.19413386e-05 + T * (  2.22789039e-08
			+ T * -3.34178948e-12 ) ) ) );
		h[sCH] =  6.38680289e+02 * (
			T * (  3.20020200e+00 + T * (  1.03643800e-03
			+ T * ( -1.71147700e-06 + T * (  1.43347250e-09
			+ T * -3.91106600e-13 ) ) ) ) +  7.04525900e+04 );
		cp[sCH] =  6.38680289e+02 * (
			 3.20020200e+00 + T * (  2.07287600e-03 
			+ T * ( -5.13443100e-06 + T * (  5.73389000e-09
			+ T * -1.95553300e-12 ) ) ) );
		h[sC3H5XTXC3H5] =  2.02443146e+02 * (
			T * (  2.17916644e+00 + T * (  1.01913311e-02
			+ T * ( -2.63804611e-06 + T * (  1.19226547e-10
			+ T *  5.40797072e-14 ) ) ) ) +  2.94895338e+04 );
		cp[sC3H5XTXC3H5] =  2.02443146e+02 * (
			 2.17916644e+00 + T * (  2.03826623e-02 
			+ T * ( -7.91413834e-06 + T * (  4.76906187e-10
			+ T *  2.70398536e-13 ) ) ) );
		h[sIXC3H7] =  1.92970803e+02 * (
			T * (  1.71330000e+00 + T * (  1.27130800e-02
			+ T * (  5.26936000e-07 + T * ( -4.55321500e-09
			+ T *  1.76554200e-12 ) ) ) ) +  7.53580900e+03 );
		cp[sIXC3H7] =  1.92970803e+02 * (
			 1.71330000e+00 + T * (  2.54261600e-02 
			+ T * (  1.58080800e-06 + T * ( -1.82128600e-08
			+ T *  8.82771000e-12 ) ) ) );
		h[sIXC4H8] =  1.48195138e+02 * (
			T * (  9.38433173e-01 + T * (  1.95273643e-02
			+ T * ( -7.21457160e-06 + T * (  1.46816769e-09
			+ T * -1.22887096e-13 ) ) ) ) -3.74817891e+03 );
		cp[sIXC4H8] =  1.48195138e+02 * (
			 9.38433173e-01 + T * (  3.90547287e-02 
			+ T * ( -2.16437148e-05 + T * (  5.87267077e-09
			+ T * -6.14435479e-13 ) ) ) );
		h[sIXC4H7] =  1.50906418e+02 * (
			T * ( -7.20881697e-04 + T * (  2.18247865e-02
			+ T * ( -1.05461959e-05 + T * (  3.09962457e-09
			+ T * -4.08756720e-13 ) ) ) ) +  1.45717785e+04 );
		cp[sIXC4H7] =  1.50906418e+02 * (
			-7.20881697e-04 + T * (  4.36495730e-02 
			+ T * ( -3.16385877e-05 + T * (  1.23984983e-08
			+ T * -2.04378360e-12 ) ) ) );
		h[sTXC4H9] =  1.45579563e+02 * (
			T * ( -2.73729203e+00 + T * (  2.27695172e-02
			+ T * ( -7.54636463e-06 + T * (  1.14237763e-09
			+ T * -3.10643590e-14 ) ) ) ) +  4.87138887e+03 );
		cp[sTXC4H9] =  1.45579563e+02 * (
			-2.73729203e+00 + T * (  4.55390345e-02 
			+ T * ( -2.26390939e-05 + T * (  4.56951052e-09
			+ T * -1.55321795e-13 ) ) ) );
		h[sIXC4H9] =  1.45579563e+02 * (
			T * ( -1.00837470e+00 + T * (  2.50543502e-02
			+ T * ( -1.13980403e-05 + T * (  3.15064458e-09
			+ T * -3.95104106e-13 ) ) ) ) +  6.65557909e+03 );
		cp[sIXC4H9] =  1.45579563e+02 * (
			-1.00837470e+00 + T * (  5.01087003e-02 
			+ T * ( -3.41941210e-05 + T * (  1.26025783e-08
			+ T * -1.97552053e-12 ) ) ) );
		h[sCH3COCH3] =  1.43158167e+02 * (
			T * (  1.30767163e+00 + T * (  1.46010871e-02
			+ T * ( -3.96818723e-06 + T * (  1.63037522e-10
			+ T *  9.35502406e-14 ) ) ) ) -2.75328269e+04 );
		cp[sCH3COCH3] =  1.43158167e+02 * (
			 1.30767163e+00 + T * (  2.92021742e-02 
			+ T * ( -1.19045617e-05 + T * (  6.52150087e-10
			+ T *  4.67751203e-13 ) ) ) );
		h[sC2H5CO] =  1.45686701e+02 * (
			T * (  2.93313946e+00 + T * (  1.23713956e-02
			+ T * ( -4.04073343e-06 + T * (  5.86752620e-10
			+ T * -1.05368668e-14 ) ) ) ) -5.69663063e+03 );
		cp[sC2H5CO] =  1.45686701e+02 * (
			 2.93313946e+00 + T * (  2.47427911e-02 
			+ T * ( -1.21222003e-05 + T * (  2.34701048e-09
			+ T * -5.26843338e-14 ) ) ) );
		h[sNEOXC5H11] =  1.16876212e+02 * (
			T * ( -1.58140132e+00 + T * (  3.28587533e-02
			+ T * ( -1.56040105e-05 + T * (  4.36831982e-09
			+ T * -5.37419850e-13 ) ) ) ) +  2.30933742e+03 );
		cp[sNEOXC5H11] =  1.16876212e+02 * (
			-1.58140132e+00 + T * (  6.57175067e-02 
			+ T * ( -4.68120314e-05 + T * (  1.74732793e-08
			+ T * -2.68709925e-12 ) ) ) );
		h[sIXC4H7O] =  1.16945257e+02 * (
			T * (  1.32646878e+00 + T * (  1.57581752e-02
			+ T * ( -6.88479067e-06 + T * (  1.90577798e-09
			+ T * -2.36963878e-13 ) ) ) ) +  2.87093070e+03 );
		cp[sIXC4H7O] =  1.16945257e+02 * (
			 1.32646878e+00 + T * (  3.15163503e-02 
			+ T * ( -2.06543720e-05 + T * (  7.62311192e-09
			+ T * -1.18481939e-12 ) ) ) );
		h[sIXC4H8O] =  1.15310385e+02 * (
			T * ( -2.97373741e+00 + T * (  3.12309319e-02
			+ T * ( -1.68116070e-05 + T * (  5.33363400e-09
			+ T * -7.34765648e-13 ) ) ) ) -1.73274318e+04 );
		cp[sIXC4H8O] =  1.15310385e+02 * (
			-2.97373741e+00 + T * (  6.24618637e-02 
			+ T * ( -5.04348211e-05 + T * (  2.13345360e-08
			+ T * -3.67382824e-12 ) ) ) );
		h[sIC4H8OOHXIO2XC4H9O4] =  6.86500099e+01 * (
			T * (  5.11703412e+00 + T * (  2.81609880e-02
			+ T * ( -1.22353513e-05 + T * (  3.06850708e-09
			+ T * -3.36953588e-13 ) ) ) ) -2.36073583e+04 );
		cp[sIC4H8OOHXIO2XC4H9O4] =  6.86500099e+01 * (
			 5.11703412e+00 + T * (  5.63219760e-02 
			+ T * ( -3.67060538e-05 + T * (  1.22740283e-08
			+ T * -1.68476794e-12 ) ) ) );
		h[sIC4H8OOHXIXC4H9O2] =  9.33021366e+01 * (
			T * (  2.14188719e+00 + T * (  2.69638383e-02
			+ T * ( -1.17639704e-05 + T * (  2.89268730e-09
			+ T * -3.03239606e-13 ) ) ) ) -4.70862582e+03 );
		cp[sIC4H8OOHXIXC4H9O2] =  9.33021366e+01 * (
			 2.14188719e+00 + T * (  5.39276766e-02 
			+ T * ( -3.52919112e-05 + T * (  1.15707492e-08
			+ T * -1.51619803e-12 ) ) ) );
		h[sIC4H8OOHXTO2XC4H9O4] =  6.86500099e+01 * (
			T * (  4.68534433e+00 + T * (  3.09050930e-02
			+ T * ( -1.56669145e-05 + T * (  4.78424157e-09
			+ T * -6.46457298e-13 ) ) ) ) -2.64124548e+04 );
		cp[sIC4H8OOHXTO2XC4H9O4] =  6.86500099e+01 * (
			 4.68534433e+00 + T * (  6.18101860e-02 
			+ T * ( -4.70007436e-05 + T * (  1.91369663e-08
			+ T * -3.23228649e-12 ) ) ) );
		h[sIC4H8OOHXTXC4H9O2] =  9.33021366e+01 * (
			T * (  4.65203922e+00 + T * (  1.99516401e-02
			+ T * ( -4.96749033e-06 + T * (  1.02978369e-11
			+ T *  1.49019577e-13 ) ) ) ) -7.40781120e+03 );
		cp[sIC4H8OOHXTXC4H9O2] =  9.33021366e+01 * (
			 4.65203922e+00 + T * (  3.99032801e-02 
			+ T * ( -1.49024710e-05 + T * (  4.11913477e-11
			+ T *  7.45097887e-13 ) ) ) );
		h[sIC4KETIIXC4H8O3] =  7.98657112e+01 * (
			T * (  2.25615224e+00 + T * (  2.80658977e-02
			+ T * ( -1.24083426e-05 + T * (  2.97693002e-09
			+ T * -2.90746368e-13 ) ) ) ) -3.88957414e+04 );
		cp[sIC4KETIIXC4H8O3] =  7.98657112e+01 * (
			 2.25615224e+00 + T * (  5.61317954e-02 
			+ T * ( -3.72250277e-05 + T * (  1.19077201e-08
			+ T * -1.45373184e-12 ) ) ) );
		h[sTC4H8OOHXIXC4H9O2] =  9.33021366e+01 * (
			T * (  2.39936083e+00 + T * (  2.85793727e-02
			+ T * ( -1.46740864e-05 + T * (  4.56553797e-09
			+ T * -6.28996546e-13 ) ) ) ) -7.63231691e+03 );
		cp[sTC4H8OOHXIXC4H9O2] =  9.33021366e+01 * (
			 2.39936083e+00 + T * (  5.71587455e-02 
			+ T * ( -4.40222592e-05 + T * (  1.82621519e-08
			+ T * -3.14498273e-12 ) ) ) );
		h[sTXC4H9O] =  1.13720593e+02 * (
			T * ( -7.70464068e-01 + T * (  2.90963330e-02
			+ T * ( -1.45399307e-05 + T * (  4.37355548e-09
			+ T * -5.83672080e-13 ) ) ) ) -1.36502805e+04 );
		cp[sTXC4H9O] =  1.13720593e+02 * (
			-7.70464068e-01 + T * (  5.81926660e-02 
			+ T * ( -4.36197921e-05 + T * (  1.74942219e-08
			+ T * -2.91836040e-12 ) ) ) );
		h[sIXC4H9O2] =  9.33021366e+01 * (
			T * (  2.16441645e+00 + T * (  2.53129504e-02
			+ T * ( -1.05578196e-05 + T * (  2.61954310e-09
			+ T * -2.94856350e-13 ) ) ) ) -1.22481418e+04 );
		cp[sIXC4H9O2] =  9.33021366e+01 * (
			 2.16441645e+00 + T * (  5.06259009e-02 
			+ T * ( -3.16734589e-05 + T * (  1.04781724e-08
			+ T * -1.47428175e-12 ) ) ) );
		h[sTXC4H9O2] =  9.33021366e+01 * (
			T * (  1.67446795e+00 + T * (  2.82304004e-02
			+ T * ( -1.41712979e-05 + T * (  4.41856535e-09
			+ T * -6.18559950e-13 ) ) ) ) -1.50470809e+04 );
		cp[sTXC4H9O2] =  9.33021366e+01 * (
			 1.67446795e+00 + T * (  5.64608009e-02 
			+ T * ( -4.25138936e-05 + T * (  1.76742614e-08
			+ T * -3.09279975e-12 ) ) ) );
		h[sTXC3H6CHO] =  1.16945257e+02 * (
			T * (  1.87052762e+00 + T * (  2.07434839e-02
			+ T * ( -8.89385670e-06 + T * (  2.25382902e-09
			+ T * -2.55741266e-13 ) ) ) ) -8.97730744e+03 );
		cp[sTXC3H6CHO] =  1.16945257e+02 * (
			 1.87052762e+00 + T * (  4.14869677e-02 
			+ T * ( -2.66815701e-05 + T * (  9.01531610e-09
			+ T * -1.27870633e-12 ) ) ) );
		h[sTXC3H6OCHO] =  9.54617893e+01 * (
			T * (  3.70830259e-01 + T * (  2.69237831e-02
			+ T * ( -1.27492522e-05 + T * (  3.32205592e-09
			+ T * -3.58457460e-13 ) ) ) ) -2.18391262e+04 );
		cp[sTXC3H6OCHO] =  9.54617893e+01 * (
			 3.70830259e-01 + T * (  5.38475661e-02 
			+ T * ( -3.82477565e-05 + T * (  1.32882237e-08
			+ T * -1.79228730e-12 ) ) ) );
		h[sIXC3H5CHO] =  1.18627154e+02 * (
			T * (  6.27183793e-01 + T * (  2.33390127e-02
			+ T * ( -1.24810210e-05 + T * (  3.95826355e-09
			+ T * -5.47904310e-13 ) ) ) ) -1.57203117e+04 );
		cp[sIXC3H5CHO] =  1.18627154e+02 * (
			 6.27183793e-01 + T * (  4.66780254e-02 
			+ T * ( -3.74430631e-05 + T * (  1.58330542e-08
			+ T * -2.73952155e-12 ) ) ) );
		h[sIXC3H7CHO] =  1.15310385e+02 * (
			T * ( -2.73021382e-01 + T * (  2.44848154e-02
			+ T * ( -1.04256683e-05 + T * (  2.50132363e-09
			+ T * -2.55024148e-13 ) ) ) ) -2.76054737e+04 );
		cp[sIXC3H7CHO] =  1.15310385e+02 * (
			-2.73021382e-01 + T * (  4.89696307e-02 
			+ T * ( -3.12770049e-05 + T * (  1.00052945e-08
			+ T * -1.27512074e-12 ) ) ) );
		h[sIXC4H7OOH] =  9.43696086e+01 * (
			T * (  2.99117402e+00 + T * (  2.51674639e-02
			+ T * ( -1.18760020e-05 + T * (  3.34880385e-09
			+ T * -4.22106818e-13 ) ) ) ) -1.51095046e+04 );
		cp[sIXC4H7OOH] =  9.43696086e+01 * (
			 2.99117402e+00 + T * (  5.03349278e-02 
			+ T * ( -3.56280061e-05 + T * (  1.33952154e-08
			+ T * -2.11053409e-12 ) ) ) );
		h[sIXC3H5CO] =  1.20358135e+02 * (
			T * (  1.40475057e+00 + T * (  2.09346201e-02
			+ T * ( -1.14633650e-05 + T * (  3.69373333e-09
			+ T * -5.16243830e-13 ) ) ) ) +  2.80734533e+03 );
		cp[sIXC3H5CO] =  1.20358135e+02 * (
			 1.40475057e+00 + T * (  4.18692403e-02 
			+ T * ( -3.43900951e-05 + T * (  1.47749333e-08
			+ T * -2.58121915e-12 ) ) ) );
		h[sYXC7H15] =  8.38223611e+01 * (
			T * ( -3.88343746e+00 + T * (  4.19967208e-02
			+ T * ( -1.65431512e-05 + T * (  3.61892672e-09
			+ T * -3.32033050e-13 ) ) ) ) -4.84718312e+03 );
		cp[sYXC7H15] =  8.38223611e+01 * (
			-3.88343746e+00 + T * (  8.39934415e-02 
			+ T * ( -4.96294536e-05 + T * (  1.44757069e-08
			+ T * -1.66016525e-12 ) ) ) );
		h[sYXC7H14] =  8.46829358e+01 * (
			T * ( -8.42232644e-01 + T * (  3.94899149e-02
			+ T * ( -1.68191492e-05 + T * (  4.22335120e-09
			+ T * -4.74403974e-13 ) ) ) ) -1.45971538e+04 );
		cp[sYXC7H14] =  8.46829358e+01 * (
			-8.42232644e-01 + T * (  7.89798297e-02 
			+ T * ( -5.04574475e-05 + T * (  1.68934048e-08
			+ T * -2.37201987e-12 ) ) ) );
		h[sXXC7H13] =  8.55613642e+01 * (
			T * ( -1.02563925e+00 + T * (  3.93160690e-02
			+ T * ( -1.83228719e-05 + T * (  5.31662222e-09
			+ T * -7.09728908e-13 ) ) ) ) +  1.25714263e+03 );
		cp[sXXC7H13] =  8.55613642e+01 * (
			-1.02563925e+00 + T * (  7.86321380e-02 
			+ T * ( -5.49686158e-05 + T * (  2.12664889e-08
			+ T * -3.54864454e-12 ) ) ) );
		h[sIXC8H18] =  7.27897815e+01 * (
			T * ( -4.20868893e+00 + T * (  5.57202905e-02
			+ T * ( -2.63782194e-05 + T * (  7.31015605e-09
			+ T * -8.87486382e-13 ) ) ) ) -2.99446875e+04 );
		cp[sIXC8H18] =  7.27897815e+01 * (
			-4.20868893e+00 + T * (  1.11440581e-01 
			+ T * ( -7.91346582e-05 + T * (  2.92406242e-08
			+ T * -4.43743191e-12 ) ) ) );
		h[sAXC8H17] =  7.34378533e+01 * (
			T * ( -3.41944741e+00 + T * (  5.34015945e-02
			+ T * ( -2.55137188e-05 + T * (  7.13354205e-09
			+ T * -8.72957298e-13 ) ) ) ) -5.33514196e+03 );
		cp[sAXC8H17] =  7.34378533e+01 * (
			-3.41944741e+00 + T * (  1.06803189e-01 
			+ T * ( -7.65411563e-05 + T * (  2.85341682e-08
			+ T * -4.36478649e-12 ) ) ) );
		h[sBXC8H17] =  7.34378533e+01 * (
			T * ( -3.09104262e+00 + T * (  5.11594480e-02
			+ T * ( -2.28286291e-05 + T * (  5.75459850e-09
			+ T * -6.14026160e-13 ) ) ) ) -6.62829069e+03 );
		cp[sBXC8H17] =  7.34378533e+01 * (
			-3.09104262e+00 + T * (  1.02318896e-01 
			+ T * ( -6.84858873e-05 + T * (  2.30183940e-08
			+ T * -3.07013080e-12 ) ) ) );
		h[sCXC8H17] =  7.34378533e+01 * (
			T * ( -9.73159697e-02 + T * (  4.46326862e-02
			+ T * ( -1.70957938e-05 + T * (  3.44101320e-09
			+ T * -2.55576792e-13 ) ) ) ) -8.91212728e+03 );
		cp[sCXC8H17] =  7.34378533e+01 * (
			-9.73159697e-02 + T * (  8.92653724e-02 
			+ T * ( -5.12873814e-05 + T * (  1.37640528e-08
			+ T * -1.27788396e-12 ) ) ) );
		h[sDXC8H17] =  7.34378533e+01 * (
			T * ( -3.41944741e+00 + T * (  5.34015945e-02
			+ T * ( -2.55137188e-05 + T * (  7.13354205e-09
			+ T * -8.72957298e-13 ) ) ) ) -5.73775897e+03 );
		cp[sDXC8H17] =  7.34378533e+01 * (
			-3.41944741e+00 + T * (  1.06803189e-01 
			+ T * ( -7.65411563e-05 + T * (  2.85341682e-08
			+ T * -4.36478649e-12 ) ) ) );
		h[sAXC8H17O2] =  5.72549857e+01 * (
			T * ( -4.17981776e-01 + T * (  5.45165365e-02
			+ T * ( -2.57893830e-05 + T * (  7.19626250e-09
			+ T * -8.85443976e-13 ) ) ) ) -2.42387571e+04 );
		cp[sAXC8H17O2] =  5.72549857e+01 * (
			-4.17981776e-01 + T * (  1.09033073e-01 
			+ T * ( -7.73681490e-05 + T * (  2.87850500e-08
			+ T * -4.42721988e-12 ) ) ) );
		h[sAXC8H17O2H] =  5.68602965e+01 * (
			T * ( -1.41990836e+00 + T * (  5.90138665e-02
			+ T * ( -2.85849639e-05 + T * (  8.04937108e-09
			+ T * -9.87307934e-13 ) ) ) ) -4.12827930e+04 );
		cp[sAXC8H17O2H] =  5.68602965e+01 * (
			-1.41990836e+00 + T * (  1.18027733e-01 
			+ T * ( -8.57548917e-05 + T * (  3.21974843e-08
			+ T * -4.93653967e-12 ) ) ) );
		h[sCXC8H17O2] =  5.72549857e+01 * (
			T * ( -8.81398118e-01 + T * (  5.73144075e-02
			+ T * ( -2.92656130e-05 + T * (  8.93021525e-09
			+ T * -1.19787280e-12 ) ) ) ) -2.70381270e+04 );
		cp[sCXC8H17O2] =  5.72549857e+01 * (
			-8.81398118e-01 + T * (  1.14628815e-01 
			+ T * ( -8.77968391e-05 + T * (  3.57208610e-08
			+ T * -5.98936398e-12 ) ) ) );
		h[sCXC8H17O2H] =  5.68602965e+01 * (
			T * ( -1.80339987e+00 + T * (  6.16548540e-02
			+ T * ( -3.19177854e-05 + T * (  9.72183205e-09
			+ T * -1.28961253e-12 ) ) ) ) -4.40953853e+04 );
		cp[sCXC8H17O2H] =  5.68602965e+01 * (
			-1.80339987e+00 + T * (  1.23309708e-01 
			+ T * ( -9.57533561e-05 + T * (  3.88873282e-08
			+ T * -6.44806264e-12 ) ) ) );
		h[sAXC8H17O] =  6.43445084e+01 * (
			T * ( -3.80491312e+00 + T * (  5.75923045e-02
			+ T * ( -2.85871395e-05 + T * (  8.29694817e-09
			+ T * -1.04895079e-12 ) ) ) ) -2.20474923e+04 );
		cp[sAXC8H17O] =  6.43445084e+01 * (
			-3.80491312e+00 + T * (  1.15184609e-01 
			+ T * ( -8.57614185e-05 + T * (  3.31877927e-08
			+ T * -5.24475393e-12 ) ) ) );
		h[sCXC8H17O] =  6.43445084e+01 * (
			T * ( -3.83187637e+00 + T * (  5.92907770e-02
			+ T * ( -3.06963249e-05 + T * (  9.35788083e-09
			+ T * -1.24344684e-12 ) ) ) ) -2.55647052e+04 );
		cp[sCXC8H17O] =  6.43445084e+01 * (
			-3.83187637e+00 + T * (  1.18581554e-01 
			+ T * ( -9.20889748e-05 + T * (  3.74315233e-08
			+ T * -6.21723419e-12 ) ) ) );
		h[sAC8H16OOHXAXC8H17O2] =  5.72549857e+01 * (
			T * ( -6.38717337e-01 + T * (  5.67115800e-02
			+ T * ( -2.77353200e-05 + T * (  7.88123630e-09
			+ T * -9.74664034e-13 ) ) ) ) -1.70745917e+04 );
		cp[sAC8H16OOHXAXC8H17O2] =  5.72549857e+01 * (
			-6.38717337e-01 + T * (  1.13423160e-01 
			+ T * ( -8.32059600e-05 + T * (  3.15249452e-08
			+ T * -4.87332017e-12 ) ) ) );
		h[sAC8H16OOHXBXC8H17O2] =  5.72549857e+01 * (
			T * ( -6.93031347e-01 + T * (  5.51507470e-02
			+ T * ( -2.54658905e-05 + T * (  6.54188850e-09
			+ T * -6.99532598e-13 ) ) ) ) -1.79009147e+04 );
		cp[sAC8H16OOHXBXC8H17O2] =  5.72549857e+01 * (
			-6.93031347e-01 + T * (  1.10301494e-01 
			+ T * ( -7.63976715e-05 + T * (  2.61675540e-08
			+ T * -3.49766299e-12 ) ) ) );
		h[sAC8H16OOHXCXC8H17O2] =  5.72549857e+01 * (
			T * (  2.58309817e+00 + T * (  4.81702185e-02
			+ T * ( -1.95346384e-05 + T * (  4.28338858e-09
			+ T * -3.72966342e-13 ) ) ) ) -2.01335775e+04 );
		cp[sAC8H16OOHXCXC8H17O2] =  5.72549857e+01 * (
			 2.58309817e+00 + T * (  9.63404369e-02 
			+ T * ( -5.86039153e-05 + T * (  1.71335543e-08
			+ T * -1.86483171e-12 ) ) ) );
		h[sCC8H16OOHXAXC8H17O2] =  5.72549857e+01 * (
			T * ( -1.09111518e+00 + T * (  5.94785170e-02
			+ T * ( -3.11836806e-05 + T * (  9.60600030e-09
			+ T * -1.28616809e-12 ) ) ) ) -1.98751651e+04 );
		cp[sCC8H16OOHXAXC8H17O2] =  5.72549857e+01 * (
			-1.09111518e+00 + T * (  1.18957034e-01 
			+ T * ( -9.35510419e-05 + T * (  3.84240012e-08
			+ T * -6.43084043e-12 ) ) ) );
		h[sAC8H16OOHXAO2XC8H17O4] =  4.69164184e+01 * (
			T * (  2.15186702e+00 + T * (  5.72483635e-02
			+ T * ( -2.61442177e-05 + T * (  6.69850128e-09
			+ T * -7.27569306e-13 ) ) ) ) -3.55015670e+04 );
		cp[sAC8H16OOHXAO2XC8H17O4] =  4.69164184e+01 * (
			 2.15186702e+00 + T * (  1.14496727e-01 
			+ T * ( -7.84326532e-05 + T * (  2.67940051e-08
			+ T * -3.63784653e-12 ) ) ) );
		h[sAC8H16OOHXBO2XC8H17O4] =  4.69164184e+01 * (
			T * (  1.32679918e+00 + T * (  5.97158645e-02
			+ T * ( -2.86727328e-05 + T * (  7.84412273e-09
			+ T * -9.22268938e-13 ) ) ) ) -3.75654653e+04 );
		cp[sAC8H16OOHXBO2XC8H17O4] =  4.69164184e+01 * (
			 1.32679918e+00 + T * (  1.19431729e-01 
			+ T * ( -8.60181983e-05 + T * (  3.13764909e-08
			+ T * -4.61134469e-12 ) ) ) );
		h[sAC8H16OOHXCO2XC8H17O4] =  4.69164184e+01 * (
			T * (  1.69476791e+00 + T * (  6.00660575e-02
			+ T * ( -2.95349242e-05 + T * (  8.35755413e-09
			+ T * -1.02287168e-12 ) ) ) ) -3.83062977e+04 );
		cp[sAC8H16OOHXCO2XC8H17O4] =  4.69164184e+01 * (
			 1.69476791e+00 + T * (  1.20132115e-01 
			+ T * ( -8.86047726e-05 + T * (  3.34302165e-08
			+ T * -5.11435840e-12 ) ) ) );
		h[sCC8H16OOHXAO2XC8H17O4] =  4.69164184e+01 * (
			T * (  1.69476791e+00 + T * (  6.00660575e-02
			+ T * ( -2.95349242e-05 + T * (  8.35755413e-09
			+ T * -1.02287168e-12 ) ) ) ) -3.83062977e+04 );
		cp[sCC8H16OOHXAO2XC8H17O4] =  4.69164184e+01 * (
			 1.69476791e+00 + T * (  1.20132115e-01 
			+ T * ( -8.86047726e-05 + T * (  3.34302165e-08
			+ T * -5.11435840e-12 ) ) ) );
		h[sIC8KETABXC8H16O3] =  5.18971587e+01 * (
			T * ( -1.15294347e+00 + T * (  5.97987060e-02
			+ T * ( -3.05910575e-05 + T * (  9.03308852e-09
			+ T * -1.15247693e-12 ) ) ) ) -5.37517462e+04 );
		cp[sIC8KETABXC8H16O3] =  5.18971587e+01 * (
			-1.15294347e+00 + T * (  1.19597412e-01 
			+ T * ( -9.17731724e-05 + T * (  3.61323541e-08
			+ T * -5.76238463e-12 ) ) ) );
		h[sC6H5CH2] =  9.12400413e+01 * (
			T * (  4.81115400e-01 + T * (  1.92564160e-02
			+ T * (  1.09538307e-05 + T * ( -1.92431803e-08
			+ T *  7.08461360e-12 ) ) ) ) +  2.33070270e+04 );
		cp[sC6H5CH2] =  9.12400413e+01 * (
			 4.81115400e-01 + T * (  3.85128320e-02 
			+ T * (  3.28614920e-05 + T * ( -7.69727210e-08
			+ T *  3.54230680e-11 ) ) ) );
		h[sC6H5] =  1.07838392e+02 * (
			T * ( -2.34050750e+00 + T * (  2.13801525e-02
			+ T * ( -8.50605533e-07 + T * ( -7.66717900e-09
			+ T *  3.24910380e-12 ) ) ) ) +  3.83767340e+04 );
		cp[sC6H5] =  1.07838392e+02 * (
			-2.34050750e+00 + T * (  4.27603050e-02 
			+ T * ( -2.55181660e-06 + T * ( -3.06687160e-08
			+ T *  1.62455190e-11 ) ) ) );
		h[sC6H6] =  1.06446715e+02 * (
			T * ( -3.52109420e+00 + T * (  4.03026715e-02
			+ T * ( -3.36291400e-05 + T * (  1.71465552e-08
			+ T * -3.79215020e-12 ) ) ) ) +  4.10422350e+04 );
		cp[sC6H6] =  1.06446715e+02 * (
			-3.52109420e+00 + T * (  8.06053430e-02 
			+ T * ( -1.00887420e-04 + T * (  6.85862210e-08
			+ T * -1.89607510e-11 ) ) ) );
		h[sC6H5OH] =  8.83489183e+01 * (
			T * ( -1.69565390e+00 + T * (  2.61356495e-02
			+ T * ( -2.40080167e-06 + T * ( -8.96490075e-09
			+ T *  4.08981460e-12 ) ) ) ) -1.32841210e+04 );
		cp[sC6H5OH] =  8.83489183e+01 * (
			-1.69565390e+00 + T * (  5.22712990e-02 
			+ T * ( -7.20240500e-06 + T * ( -3.58596030e-08
			+ T *  2.04490730e-11 ) ) ) );
		h[sC6H5O] =  8.93054780e+01 * (
			T * ( -1.82194330e+00 + T * (  2.40612550e-02
			+ T * ( -1.55974340e-06 + T * ( -8.50464850e-09
			+ T *  3.72992740e-12 ) ) ) ) +  4.24291800e+03 );
		cp[sC6H5O] =  8.93054780e+01 * (
			-1.82194330e+00 + T * (  4.81225100e-02 
			+ T * ( -4.67923020e-06 + T * ( -3.40185940e-08
			+ T *  1.86496370e-11 ) ) ) );
		h[sC6H5CHO] =  7.83499501e+01 * (
			T * ( -2.37082285e+00 + T * (  3.14421564e-02
			+ T * ( -1.42153585e-05 + T * (  3.48540207e-09
			+ T * -3.48949898e-13 ) ) ) ) -6.11656186e+03 );
		cp[sC6H5CHO] =  7.83499501e+01 * (
			-2.37082285e+00 + T * (  6.28843128e-02 
			+ T * ( -4.26460754e-05 + T * (  1.39416083e-08
			+ T * -1.74474949e-12 ) ) ) );
		h[sC6H5CO] =  7.91013224e+01 * (
			T * ( -1.56729686e+00 + T * (  2.99908756e-02
			+ T * ( -1.43051195e-05 + T * (  3.71667135e-09
			+ T * -3.97876906e-13 ) ) ) ) +  1.12873830e+04 );
		cp[sC6H5CO] =  7.91013224e+01 * (
			-1.56729686e+00 + T * (  5.99817511e-02 
			+ T * ( -4.29153585e-05 + T * (  1.48666854e-08
			+ T * -1.98938453e-12 ) ) ) );
		h[sC4H5XIXC4H5] =  1.56637905e+02 * (
			T * (  3.87944300e+00 + T * (  9.98832000e-03
			+ T * (  6.24259000e-07 + T * ( -2.32673825e-09
			+ T *  4.77223200e-13 ) ) ) ) +  3.52685900e+04 );
		cp[sC4H5XIXC4H5] =  1.56637905e+02 * (
			 3.87944300e+00 + T * (  1.99766400e-02 
			+ T * (  1.87277700e-06 + T * ( -9.30695300e-09
			+ T *  2.38611600e-12 ) ) ) );
		h[sC4H6] =  1.53718755e+02 * (
			T * (  3.19710800e+00 + T * (  1.01279600e-02
			+ T * (  2.17006400e-06 + T * ( -4.14610500e-09
			+ T *  1.28005640e-12 ) ) ) ) +  1.57152000e+04 );
		cp[sC4H6] =  1.53718755e+02 * (
			 3.19710800e+00 + T * (  2.02559200e-02 
			+ T * (  6.51019200e-06 + T * ( -1.65844200e-08
			+ T *  6.40028200e-12 ) ) ) );
		h[sC6H5CH2O] =  7.76127177e+01 * (
			T * ( -4.77736690e+00 + T * (  3.75524654e-02
			+ T * ( -1.89510944e-05 + T * (  5.45725072e-09
			+ T * -6.76268596e-13 ) ) ) ) +  1.26234438e+04 );
		cp[sC6H5CH2O] =  7.76127177e+01 * (
			-4.77736690e+00 + T * (  7.51049308e-02 
			+ T * ( -5.68532831e-05 + T * (  2.18290029e-08
			+ T * -3.38134298e-12 ) ) ) );
		h[sC6H5C2H5] =  7.83189525e+01 * (
			T * ( -5.93600077e+00 + T * (  4.47309208e-02
			+ T * ( -2.27076837e-05 + T * (  6.65259335e-09
			+ T * -8.42796204e-13 ) ) ) ) +  1.94317497e+03 );
		cp[sC6H5C2H5] =  7.83189525e+01 * (
			-5.93600077e+00 + T * (  8.94618416e-02 
			+ T * ( -6.81230510e-05 + T * (  2.66103734e-08
			+ T * -4.21398102e-12 ) ) ) );
		h[sC6H5CH2OH] =  7.68892300e+01 * (
			T * ( -3.72475901e+00 + T * (  3.71187522e-02
			+ T * ( -1.80493929e-05 + T * (  5.02947335e-09
			+ T * -6.04398106e-13 ) ) ) ) -1.38224993e+04 );
		cp[sC6H5CH2OH] =  7.68892300e+01 * (
			-3.72475901e+00 + T * (  7.42375044e-02 
			+ T * ( -5.41481788e-05 + T * (  2.01178934e-08
			+ T * -3.02199053e-12 ) ) ) );
		h[sC5H5] =  1.27736058e+02 * (
			T * (  9.83498220e-01 + T * (  1.68257380e-02
			+ T * ( -3.68472700e-08 + T * ( -9.18584850e-09
			+ T *  4.62823680e-12 ) ) ) ) +  2.97340990e+04 );
		cp[sC5H5] =  1.27736058e+02 * (
			 9.83498220e-01 + T * (  3.36514760e-02 
			+ T * ( -1.10541810e-07 + T * ( -3.67433940e-08
			+ T *  2.31411840e-11 ) ) ) );
		h[sC6H4O2] =  7.69191059e+01 * (
			T * (  3.79867882e+00 + T * (  1.25838285e-02
			+ T * (  1.26615639e-05 + T * ( -1.76694379e-08
			+ T *  6.12253146e-12 ) ) ) ) -1.72429606e+04 );
		cp[sC6H4O2] =  7.69191059e+01 * (
			 3.79867882e+00 + T * (  2.51676569e-02 
			+ T * (  3.79846917e-05 + T * ( -7.06777516e-08
			+ T *  3.06126573e-11 ) ) ) );
		h[sC5H6] =  1.25788072e+02 * (
			T * (  8.61089570e-01 + T * (  7.40201550e-03
			+ T * (  2.40362983e-05 + T * ( -2.83451375e-08
			+ T *  9.73799440e-12 ) ) ) ) +  1.48017550e+04 );
		cp[sC5H6] =  1.25788072e+02 * (
			 8.61089570e-01 + T * (  1.48040310e-02 
			+ T * (  7.21088950e-05 + T * ( -1.13380550e-07
			+ T *  4.86899720e-11 ) ) ) );
		h[sC5H4O] =  1.03822832e+02 * (
			T * ( -4.15438982e+00 + T * (  3.05987543e-02
			+ T * ( -1.84327408e-05 + T * (  6.24350998e-09
			+ T * -8.92877106e-13 ) ) ) ) +  2.68411141e+03 );
		cp[sC5H4O] =  1.03822832e+02 * (
			-4.15438982e+00 + T * (  6.11975086e-02 
			+ T * ( -5.52982224e-05 + T * (  2.49740399e-08
			+ T * -4.46438553e-12 ) ) ) );
		h[sC5H5O] =  1.02532248e+02 * (
			T * ( -2.83112840e+00 + T * (  2.83638644e-02
			+ T * ( -1.48252434e-05 + T * (  4.37311117e-09
			+ T * -5.52009694e-13 ) ) ) ) +  2.04992154e+04 );
		cp[sC5H5O] =  1.02532248e+02 * (
			-2.83112840e+00 + T * (  5.67277287e-02 
			+ T * ( -4.44757303e-05 + T * (  1.74924447e-08
			+ T * -2.76004847e-12 ) ) ) );
		h[sC4H5XNXC4H5] =  1.56637905e+02 * (
			T * (  2.99524000e+00 + T * (  1.14422800e-02
			+ T * (  6.58490333e-07 + T * ( -2.87061250e-09
			+ T *  6.39564800e-13 ) ) ) ) +  4.14221800e+04 );
		cp[sC4H5XNXC4H5] =  1.56637905e+02 * (
			 2.99524000e+00 + T * (  2.28845600e-02 
			+ T * (  1.97547100e-06 + T * ( -1.14824500e-08
			+ T *  3.19782400e-12 ) ) ) );
		h[sCH2CHCHCO] =  1.22140381e+02 * (
			T * ( -1.74329010e-01 + T * (  2.49941346e-02
			+ T * ( -1.69451194e-05 + T * (  6.43459755e-09
			+ T * -1.00883216e-12 ) ) ) ) -8.53247629e+02 );
		cp[sCH2CHCHCO] =  1.22140381e+02 * (
			-1.74329010e-01 + T * (  4.99882693e-02 
			+ T * ( -5.08353583e-05 + T * (  2.57383902e-08
			+ T * -5.04416081e-12 ) ) ) );
		h[sCHCHCHO] =  1.51021542e+02 * (
			T * (  2.03583477e+00 + T * (  1.27869923e-02
			+ T * ( -6.05711080e-06 + T * (  1.57137786e-09
			+ T * -1.69244471e-13 ) ) ) ) +  1.73077219e+00 );
		cp[sCHCHCHO] =  1.51021542e+02 * (
			 2.03583477e+00 + T * (  2.55739847e-02 
			+ T * ( -1.81713324e-05 + T * (  6.28551145e-09
			+ T * -8.46222353e-13 ) ) ) );
		h[sHCCCO] =  1.56761944e+02 * (
			T * (  6.40938042e+00 + T * (  5.78575580e-03
			+ T * ( -3.42133077e-06 + T * (  1.13131095e-09
			+ T * -1.60430065e-13 ) ) ) ) +  2.70126564e+00 );
		cp[sHCCCO] =  1.56761944e+02 * (
			 6.40938042e+00 + T * (  1.15715116e-02 
			+ T * ( -1.02639923e-05 + T * (  4.52524382e-09
			+ T * -8.02150325e-13 ) ) ) );
		h[sC2H5OH] =  1.80479726e+02 * (
			T * (  4.23149045e-01 + T * (  1.46429084e-02
			+ T * ( -5.79483663e-06 + T * (  1.33409973e-09
			+ T * -1.34983450e-13 ) ) ) ) -2.95683054e+04 );
		cp[sC2H5OH] =  1.80479726e+02 * (
			 4.23149045e-01 + T * (  2.92858167e-02 
			+ T * ( -1.73845099e-05 + T * (  5.33639892e-09
			+ T * -6.74917252e-13 ) ) ) );
		h[sC2H4OH] =  1.84517088e+02 * (
			T * (  7.76157606e-01 + T * (  1.45501257e-02
			+ T * ( -7.13268563e-06 + T * (  2.04208474e-09
			+ T * -2.51183200e-13 ) ) ) ) -5.46808396e+03 );
		cp[sC2H4OH] =  1.84517088e+02 * (
			 7.76157606e-01 + T * (  2.91002514e-02 
			+ T * ( -2.13980569e-05 + T * (  8.16833894e-09
			+ T * -1.25591600e-12 ) ) ) );
		h[sCH3CHOH] =  1.84517088e+02 * (
			T * (  2.10026806e+00 + T * (  9.63133395e-03
			+ T * ( -1.72479145e-06 + T * ( -4.90476455e-10
			+ T *  1.87279412e-13 ) ) ) ) -6.12873036e+03 );
		cp[sCH3CHOH] =  1.84517088e+02 * (
			 2.10026806e+00 + T * (  1.92626679e-02 
			+ T * ( -5.17437436e-06 + T * ( -1.96190582e-09
			+ T *  9.36397060e-13 ) ) ) );
		h[sCH3CH2O] =  1.84517088e+02 * (
			T * ( -2.71296378e-01 + T * (  1.49419906e-02
			+ T * ( -6.56968493e-06 + T * (  1.59334973e-09
			+ T * -1.55593011e-13 ) ) ) ) -3.16397196e+03 );
		cp[sCH3CH2O] =  1.84517088e+02 * (
			-2.71296378e-01 + T * (  2.98839812e-02 
			+ T * ( -1.97090548e-05 + T * (  6.37339893e-09
			+ T * -7.77965054e-13 ) ) ) );
		h[sCH3HCO] =  1.88739217e+02 * (
			T * (  2.50569500e+00 + T * (  6.68495500e-03
			+ T * (  1.55731767e-06 + T * ( -2.82035000e-09
			+ T *  8.52713200e-13 ) ) ) ) -2.12458900e+04 );
		cp[sCH3HCO] =  1.88739217e+02 * (
			 2.50569500e+00 + T * (  1.33699100e-02 
			+ T * (  4.67195300e-06 + T * ( -1.12814000e-08
			+ T *  4.26356600e-12 ) ) ) );
		h[sCH2HCO] =  1.93159093e+02 * (
			T * (  3.40906200e+00 + T * (  5.36928500e-03
			+ T * (  6.30497333e-07 + T * ( -1.78964575e-09
			+ T *  5.73477000e-13 ) ) ) ) +  1.52147700e+03 );
		cp[sCH2HCO] =  1.93159093e+02 * (
			 3.40906200e+00 + T * (  1.07385700e-02 
			+ T * (  1.89149200e-06 + T * ( -7.15858300e-09
			+ T *  2.86738500e-12 ) ) ) );
		h[sHOC2H4O2] =  1.07894368e+02 * (
			T * (  3.58811609e+00 + T * (  1.49707080e-02
			+ T * ( -7.83855163e-06 + T * (  2.81588813e-09
			+ T * -4.57327372e-13 ) ) ) ) -2.17409435e+04 );
		cp[sHOC2H4O2] =  1.07894368e+02 * (
			 3.58811609e+00 + T * (  2.99414161e-02 
			+ T * ( -2.35156549e-05 + T * (  1.12635525e-08
			+ T * -2.28663686e-12 ) ) ) );
		h[sC2H] =  3.32201534e+02 * (
			T * (  2.73770400e+00 + T * (  4.02422300e-03
			+ T * ( -3.08143667e-06 + T * (  1.63131475e-09
			+ T * -3.87916000e-13 ) ) ) ) +  6.68381300e+04 );
		cp[sC2H] =  3.32201534e+02 * (
			 2.73770400e+00 + T * (  8.04844600e-03 
			+ T * ( -9.24431000e-06 + T * (  6.52525900e-09
			+ T * -1.93958000e-12 ) ) ) );
		h[sH2CC] =  3.19340144e+02 * (
			T * (  3.28154830e+00 + T * (  3.48823955e-03
			+ T * ( -7.95174800e-07 + T * ( -3.02610800e-10
			+ T *  1.96379090e-13 ) ) ) ) +  4.86217940e+04 );
		cp[sH2CC] =  3.19340144e+02 * (
			 3.28154830e+00 + T * (  6.97647910e-03 
			+ T * ( -2.38552440e-06 + T * ( -1.21044320e-09
			+ T *  9.81895450e-13 ) ) ) );
		h[sC3H3] =  2.12893430e+02 * (
			T * (  4.75419900e+00 + T * (  5.54013850e-03
			+ T * (  9.31107667e-08 + T * ( -1.36980300e-09
			+ T *  3.89925800e-13 ) ) ) ) +  3.90468300e+04 );
		cp[sC3H3] =  2.12893430e+02 * (
			 4.75419900e+00 + T * (  1.10802770e-02 
			+ T * (  2.79332300e-07 + T * ( -5.47921200e-09
			+ T *  1.94962900e-12 ) ) ) );
		h[sPXC3H4] =  2.07536818e+02 * (
			T * (  3.02973000e+00 + T * (  7.49480500e-03
			+ T * ( -4.66166667e-07 + T * ( -9.92404750e-10
			+ T *  2.77643400e-13 ) ) ) ) +  2.14840800e+04 );
		cp[sPXC3H4] =  2.07536818e+02 * (
			 3.02973000e+00 + T * (  1.49896100e-02 
			+ T * ( -1.39850000e-06 + T * ( -3.96961900e-09
			+ T *  1.38821700e-12 ) ) ) );
		h[sAXC3H4] =  2.07536818e+02 * (
			T * ( -2.13196800e-01 + T * (  1.67935650e-02
			+ T * ( -1.26829000e-05 + T * (  6.86459500e-09
			+ T * -1.73800880e-12 ) ) ) ) +  2.16204800e+04 );
		cp[sAXC3H4] =  2.07536818e+02 * (
			-2.13196800e-01 + T * (  3.35871300e-02 
			+ T * ( -3.80487000e-05 + T * (  2.74583800e-08
			+ T * -8.69004400e-12 ) ) ) );
		h[s1XCH2] =  5.92780550e+02 * (
			T * (  4.19860411e+00 + T * ( -1.18330710e-03
			+ T * (  2.74432073e-06 + T * ( -1.67203995e-09
			+ T *  3.88629474e-13 ) ) ) ) +  5.04968163e+04 );
		cp[s1XCH2] =  5.92780550e+02 * (
			 4.19860411e+00 + T * ( -2.36661419e-03 
			+ T * (  8.23296220e-06 + T * ( -6.68815981e-09
			+ T *  1.94314737e-12 ) ) ) );
		h[sC] =  6.92284763e+02 * (
			T * (  2.55423955e+00 + T * ( -1.60768862e-04
			+ T * (  2.44597415e-07 + T * ( -1.83058722e-10
			+ T *  5.33042892e-14 ) ) ) ) +  8.54438832e+04 );
		cp[sC] =  6.92284763e+02 * (
			 2.55423955e+00 + T * ( -3.21537724e-04 
			+ T * (  7.33792245e-07 + T * ( -7.32234889e-10
			+ T *  2.66521446e-13 ) ) ) );
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
	return 158;
}

int GetNReactions( void )
{
	return rEnd;
}

double MAX_C(double X1, double X2)
{
  return ( (X1 > X2) ? X1 : X2 );
}
