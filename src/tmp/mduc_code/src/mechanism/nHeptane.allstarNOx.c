#include "nHeptane.allstarNOx.h"

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

	int	nSpec = 62;
	int	nSpecIn = 62;
	double	kTroe0, kTroeInf, fcTroe;
	double	RGAS = 8314.34;
	double	lgt = log( temp );
	double	rt = RGAS * temp;

	M[mM1] = 0.4 * c[sN2] + 0.4 * c[sO2] + c[sH] + c[sOH]
		 + c[sO] + c[sH2] + 6.5 * c[sH2O] + c[sHO2]
		 + c[sH2O2] + 0.75 * c[sCO] + c[sCO2] + c[sCH]
		 + c[sHCO] + c[sCH2OH] + c[s3XCH2] + c[sC2H2]
		 + c[sCH3] + c[sC2H4] + c[s1XCH2] + c[sCH2O]
		 + 3 * c[sCH4] + c[sC2H6] + c[sHCCO] + c[sC2H3]
		 + c[sC2H5] + c[sC3H4] + c[sC3H5] + c[sC3H6]
		 + c[sNXC3H7] + c[s1XC4H8] + c[sPXC4H9] + c[s1XC5H11]
		 + c[sC6H11] + c[s1XC6H12] + c[s1XC7H15] + c[s2XC7H15]
		 + 3 * c[sNXC7H16] + c[s1C7H15O2XC7H15O2] + c[s1HEOOHX2XC7H15O2] + c[s2OO7OOH1XO2C7H15O2]
		 + c[s1X2HEOOHXO2C7H14OOH] + c[sOC7OOHXOC7H13OOH] + c[s1X2HEPOOXOC7H13O];

	M[mM41] = c[sN2] + c[sO2] + c[sH] + c[sOH]
		 + c[sO] + 2 * c[sH2] + 6 * c[sH2O] + c[sHO2]
		 + c[sH2O2] + 1.5 * c[sCO] + 2 * c[sCO2] + c[sCH]
		 + c[sHCO] + c[sCH2OH] + c[s3XCH2] + c[sC2H2]
		 + c[sCH3] + c[sC2H4] + c[s1XCH2] + c[sCH2O]
		 + 2 * c[sCH4] + 3 * c[sC2H6] + c[sHCCO] + c[sC2H3]
		 + c[sC2H5] + c[sC3H4] + c[sC3H5] + c[sC3H6]
		 + c[sNXC3H7] + c[s1XC4H8] + c[sPXC4H9] + c[s1XC5H11]
		 + c[sC6H11] + c[s1XC6H12] + c[s1XC7H15] + c[s2XC7H15]
		 + c[sNXC7H16] + c[s1C7H15O2XC7H15O2] + c[s1HEOOHX2XC7H15O2] + c[s2OO7OOH1XO2C7H15O2]
		 + c[s1X2HEOOHXO2C7H14OOH] + c[sOC7OOHXOC7H13OOH] + c[s1X2HEPOOXOC7H13O] + c[sCH2]
		 + c[sCH3OH] + c[sC] + c[sC2H] + c[sN]
		 + c[sNO] + c[sN2O] + c[sNO2] + c[sNH]
		 + c[sHNO] + c[sNH2] + c[sNNH] + c[sCN]
		 + c[sNCO] + c[sHCN] + c[sHOCN] + c[sH2CN]
		 + c[sHCNN] + c[sNH3];

	M[mM40] = c[sN2] + c[sO2] + c[sH] + c[sOH]
		 + c[sO] + 2 * c[sH2] + 6 * c[sH2O] + c[sHO2]
		 + c[sH2O2] + 1.5 * c[sCO] + 2 * c[sCO2] + c[sCH]
		 + c[sHCO] + c[sCH2OH] + c[s3XCH2] + c[sC2H2]
		 + c[sCH3] + c[sC2H4] + c[s1XCH2] + c[sCH2O]
		 + 2 * c[sCH4] + 3 * c[sC2H6] + c[sHCCO] + c[sC2H3]
		 + c[sC2H5] + c[sC3H4] + c[sC3H5] + c[sC3H6]
		 + c[sNXC3H7] + c[s1XC4H8] + c[sPXC4H9] + c[s1XC5H11]
		 + c[sC6H11] + c[s1XC6H12] + c[s1XC7H15] + c[s2XC7H15]
		 + c[sNXC7H16] + c[s1C7H15O2XC7H15O2] + c[s1HEOOHX2XC7H15O2] + c[s2OO7OOH1XO2C7H15O2]
		 + c[s1X2HEOOHXO2C7H14OOH] + c[sOC7OOHXOC7H13OOH] + c[s1X2HEPOOXOC7H13O] + c[sCH2]
		 + c[sCH3OH] + c[sC] + c[sC2H] + c[sN]
		 + c[sNO] + c[sN2O] + c[sNO2] + c[sNH]
		 + c[sHNO] + c[sNH2] + c[sNNH] + c[sCN]
		 + c[sNCO] + c[sHCN] + c[sHOCN] + c[sH2CN]
		 + c[sHCNN] + c[sNH3];



	k[r1f] = 2.0000000000E+11 * exp( -70300000 / rt );
	k[r1b] = 4.0674859656E+08 * exp( 0.433091 * lgt + 1978606.19 / rt );
	k[r2f] = 5.0600000000E+01 * exp( 2.67 * lgt - 26300000 / rt );
	k[r2b] = 2.9161028451E+01 * exp( 2.63868 * lgt - 18754035.03 / rt );
	k[r3f] = 1.0000000000E+05 * exp( 1.6 * lgt - 13800000 / rt );
	k[r3b] = 1.0374701530E+06 * exp( 1.50368 * lgt - 78257627.95 / rt );
	k[r4f] = 1.5000000000E+06 * exp( 1.14 * lgt - 420000 / rt );
	k[r4b] = 2.7003157570E+07 * exp( 1.075 * lgt - 72423592.92 / rt );
	k[r5f] = 1.8000000000E+12 * exp( -1 * lgt );
	k[r5b] = 4.2896464801E+14 * exp( -0.65384 * lgt - 433084285.8 / rt );
	k[r6f] = 2.9000000000E+11 * exp( -1 * lgt );
	k[r6b] = 1.9584083292E+16 * exp( -1.11825 * lgt - 497816927 / rt );
	k[r7f] = 2.2000000000E+16 * exp( -2 * lgt );
	k[r7b] = 5.4393535658E+19 * exp( -1.75016 * lgt - 497541913.7 / rt );
	k[r8f] = 2.3000000000E+12 * exp( -0.8 * lgt );
	k[r8b] = 3.5261428498E+15 * exp( -0.799081 * lgt - 196197039.3 / rt );
	k[r9f] = 1.5000000000E+11 * exp( -4200000 / rt );
	k[r9b] = 2.7328549235E+07 * exp( 0.747015 * lgt - 161262675.3 / rt );
	k[r10f] = 2.5000000000E+10 * exp( -2900000 / rt );
	k[r10b] = 3.8861265555E+09 * exp( 0.345241 * lgt - 239787246.4 / rt );
	k[r11f] = 3.0000000000E+10 * exp( -7200000 / rt );
	k[r11b] = 9.8394282821E+07 * exp( 0.682016 * lgt - 236266268.2 / rt );
	k[r12f] = 1.8000000000E+10 * exp( 1700000 / rt );
	k[r12b] = 1.6125075469E+09 * exp( 0.313924 * lgt - 227641281.4 / rt );
	k[r13f] = 6.0000000000E+10;
	k[r13b] = 9.6761767494E+10 * exp( 0.248925 * lgt - 301344874.4 / rt );
	k[r14f] = 2.5000000000E+08 * exp( 5200000 / rt );
	k[r14b] = 3.2944353960E+10 * exp( -0.419617 * lgt - 173817058 / rt );
	k[r15f] = 3.2500000000E+16 * exp( -2 * lgt );
	k[r15b] = 3.6038864236E+25 * exp( -3.16571 * lgt - 218151422.1 / rt );
	k[r19f] = 5.4000000000E+09 * exp( -4200000 / rt );
	k[r19b] = 6.6085368415E+07 * exp( 0.668541 * lgt - 126527816.4 / rt );
	k[r20f] = 6.0000000000E+03 * exp( 1.5 * lgt + 3100000 / rt );
	k[r20b] = 3.3391130375E+10 * exp( 0.219853 * lgt - 106117412.1 / rt );
	k[r21f] = 1.5000000000E+11 * exp( -98700000 / rt );
	k[r21b] = 1.5208852508E+14 * exp( -0.533132 * lgt - 364980087.4 / rt );
	k[r22f] = 7.1000000000E+07 * exp( 19000000 / rt );
	k[r22b] = 5.4267465235E+16 * exp( -0.965305 * lgt - 515755732.9 / rt );
	k[r25] = 6.0000000000E+10;
	k[r26] = 3.4000000000E+09 * exp( -2900000 / rt );
	k[r27f] = 5.7000000000E+09 * exp( 3200000 / rt );
	k[r27b] = 2.9168083631E+18 * exp( -1.12914 * lgt - 370052367.4 / rt );
	k[r28f] = 1.5660000000E+11 * exp( -65930000 / rt );
	k[r28b] = 3.2470280470E+07 * exp( 0.227488 * lgt - 1742296.835 / rt );
	k[r33f] = 3.0000000000E+09;
	k[r33b] = 9.5364690411E+08 * exp( 0.228407 * lgt - 132009336.2 / rt );
	k[r35f] = 6.0000000000E+09 * exp( 7500000 / rt );
	k[r35b] = 8.0785783149E+08 * exp( 0.197553 * lgt - 1876479.045 / rt );
	k[r37f] = 1.2000000000E+10 * exp( -3400000 / rt );
	k[r37b] = 7.6647993862E+16 * exp( -1.243 * lgt - 554259617.7 / rt );
	k[r39f] = 4.2000000000E+10;
	k[r39b] = 5.3662222907E+17 * exp( -1.26594 * lgt - 269432043 / rt );
	k[r40f] = 1.3000000000E+10 * exp( -6200000 / rt );
	k[r40b] = 8.1440410754E+05 * exp( 0.262989 * lgt - 246276743.9 / rt );
	k[r41f] = 1.2000000000E+10 * exp( -6200000 / rt );
	k[r41b] = 9.9702653624E+14 * exp( -0.670998 * lgt - 788578441.7 / rt );
	k[r42f] = 1.2000000000E+10;
	k[r42b] = 5.1878402018E+09 * exp( -0.0101671 * lgt - 37408739.09 / rt );
	k[r43f] = 3.1000000000E+10;
	k[r43b] = 8.3958300935E+05 * exp( 0.252822 * lgt - 277485482.9 / rt );
	k[r44f] = 7.2000000000E+10;
	k[r44b] = 2.1669400497E+13 * exp( -0.587347 * lgt - 64229316.31 / rt );
	k[r45f] = 5.0000000000E+13 * exp( -320000000 / rt );
	k[r45b] = 3.3588207512E+08 * exp( 0.104758 * lgt + 58435338.47 / rt );
	k[r46f] = 2.3000000000E+07 * exp( 1.05 * lgt - 13700000 / rt );
	k[r46b] = 3.6820837009E+04 * exp( 1.50092 * lgt - 68348947.29 / rt );
	k[r47f] = 4.1500000000E+08 * exp( 0.57 * lgt - 11600000 / rt );
	k[r47b] = 3.8288313488E+05 * exp( 0.9896 * lgt - 58702982.32 / rt );
	k[r48f] = 3.4000000000E+06 * exp( 1.2 * lgt + 1900000 / rt );
	k[r48b] = 5.6470333036E+04 * exp( 1.5546 * lgt - 117206575.2 / rt );
	k[r49f] = 3.0000000000E+09 * exp( -54700000 / rt );
	k[r49b] = 4.0714689870E+09 * exp( -0.31394 * lgt - 51478758.86 / rt );
	k[r53f] = 8.4300000000E+10;
	k[r53b] = 1.1237221766E+13 * exp( -0.271329 * lgt - 295073528.7 / rt );
	kTroe0 = 6.2570000000E+17 * exp( -1.8 * lgt );
	kTroeInf = 2.1080000000E+11;
	fcTroe = 0.577 * exp( -temp / 2370 );
	k[ra34f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, -1.0 );
	kTroe0 = 4.6913273812E+22 * exp( -1.7293 * lgt - 440045149.2 / rt );
	kTroeInf = 1.5805207159E+16 * exp( 0.0707032 * lgt - 440045149.2 / rt );
	k[ra34b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, -1.0 );
	k[r55] = 2.2600000000E+11 * exp( -64800000 / rt );
	k[r57] = 3.3000000000E+08 * exp( -37400000 / rt );
	k[r58f] = 1.8000000000E+10;
	k[r58b] = 1.4129673015E+10 * exp( -0.0498658 * lgt - 144244286.7 / rt );
	k[r59f] = 3.6000000000E+09;
	k[r59b] = 1.7605972553E+11 * exp( 0.0697843 * lgt - 243848109.8 / rt );
	k[r60] = 1.0000000000E+13 * exp( -134000000 / rt );
	kTroe0 = 1.2720000000E+35 * exp( -7 * lgt - 11560000 / rt );
	kTroeInf = 1.8130000000E+10;
	fcTroe = 0.38 * exp( -temp / 73 ) + 0.62 * exp( -temp / 1180 );
	k[r36f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, -1.0 );
	kTroe0 = 1.7651615539E+46 * exp( -8.33318 * lgt - 393858761.4 / rt );
	kTroeInf = 2.5159102966E+21 * exp( -1.33318 * lgt - 382298761.4 / rt );
	k[r36b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, -1.0 );
	k[r62f] = 5.0000000000E+10 * exp( -105000000 / rt );
	k[r62b] = 1.1263303224E+07 * exp( 0.21071 * lgt + 12646403.51 / rt );
	k[r63f] = 1.8000000000E+10;
	k[r63b] = 9.6631178059E+08 * exp( 0.55687 * lgt - 315437882.3 / rt );
	k[r64f] = 4.0000000000E+07 * exp( -8900000 / rt );
	k[r64b] = 1.3814266479E+07 * exp( 0.211628 * lgt - 87450635.84 / rt );
	k[r69f] = 5.0000000000E+10 * exp( -105000000 / rt );
	k[r69b] = 1.1263303224E+07 * exp( 0.21071 * lgt + 12646403.51 / rt );
	k[r70f] = 3.0000000000E+10;
	k[r70b] = 1.6105196343E+09 * exp( 0.55687 * lgt - 315437882.3 / rt );
	k[r71f] = 1.0000000000E+10 * exp( -30000000 / rt );
	k[r71b] = 3.4535666198E+09 * exp( 0.211628 * lgt - 108550635.8 / rt );
	k[r84f] = 1.3000000000E+01 * exp( 3 * lgt - 33600000 / rt );
	k[r84b] = 4.1320233178E-02 * exp( 3.27546 * lgt - 26639136.6 / rt );
	k[r86f] = 1.6000000000E+04 * exp( 1.83 * lgt - 11600000 / rt );
	k[r86b] = 5.2761241400E+02 * exp( 2.00914 * lgt - 69096764.55 / rt );
	k[r105f] = 1.5000000000E+11;
	k[r105b] = 1.2413080129E+05 * exp( 1.64114 * lgt - 114828016.3 / rt );
	k[r106] = 9.6000000000E+10;
	k[r109f] = 2.0000000000E+05 * exp( 1.5 * lgt - 126000000 / rt );
	k[r109b] = 5.6489502505E+02 * exp( 1.71102 * lgt - 133473395.6 / rt );
	k[r111f] = 1.7200000000E+01 * exp( 2.8 * lgt - 2100000 / rt );
	k[r111b] = 1.9767804328E-05 * exp( 4.21906 * lgt - 196680018.1 / rt );
	k[r112f] = 1.7200000000E+01 * exp( 2.8 * lgt - 2100000 / rt );
	k[r112b] = 2.3887468852E+01 * exp( 2.57792 * lgt - 81852001.8 / rt );
	kTroe0 = 1.1870000000E+39 * exp( -7.5 * lgt - 190400000 / rt );
	kTroeInf = 2.0000000000E+14 * exp( -166290000 / rt );
	fcTroe = 0.35 * exp( -0 / temp );
	k[ra51f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, -1.0 );
	kTroe0 = 2.2850342328E+37 * exp( -7.94676 * lgt - 31790514.87 / rt );
	kTroeInf = 3.8500998026E+12 * exp( -0.446758 * lgt - 7680514.875 / rt );
	k[ra51b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, -1.0 );
	k[ra125] = 5.4200000000E+09;
	k[r129f] = 2.5000000000E+14 * exp( -319800000 / rt );
	k[r129b] = 7.5332457193E+08 * exp( 0.253959 * lgt - 141322711.7 / rt );
	k[r131f] = 1.7000000000E+12 * exp( -62900000 / rt );
	k[r131b] = 2.6610256106E+08 * exp( 0.700717 * lgt - 43032196.86 / rt );
	k[r134f] = 6.5000000000E+10 * exp( -24900000 / rt );
	k[r134b] = 1.0555750123E+08 * exp( 0.6044 * lgt - 69489824.81 / rt );
	kTroe0 = 1.0000000000E+13 * exp( -126000000 / rt );
	kTroeInf = 1.3000000000E+13 * exp( -167000000 / rt );
	fcTroe = 1 * exp( -temp / 422.8 ) + 0.411 * exp( -73.4 / temp );
	k[ra58f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, -1.0 );
	kTroe0 = 3.9227577100E+10 * exp( -0.122937 * lgt + 26911158.33 / rt );
	kTroeInf = 5.0995850230E+10 * exp( -0.122937 * lgt - 14088841.67 / rt );
	k[ra58b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, -1.0 );
	k[r146f] = 3.0000000000E+10;
	k[r146b] = 1.5280987618E+06 * exp( 0.911984 * lgt - 37561661.66 / rt );
	k[r149f] = 1.1000000000E+07 * exp( 6300000 / rt );
	k[r149b] = 6.6154019371E+07 * exp( -0.122018 * lgt - 36985881.02 / rt );
	k[r164] = 1.4000000000E+06 * exp( 1.5 * lgt - 31100000 / rt );
	k[r166] = 7.2000000000E+03 * exp( 2 * lgt - 3600000 / rt );
	k[r170] = 1.5000000000E-10 * exp( 6 * lgt - 25400000 / rt );
	k[r189] = 1.0000000000E+09;
	k[r190] = 1.0000000000E+09;
	k[r191f] = 3.9800000000E+13 * exp( -293100000 / rt );
	k[r191b] = 4.7045354409E+06 * exp( 0.920325 * lgt - 19705269.47 / rt );
	k[r192f] = 5.0000000000E+09;
	k[r192b] = 1.4084864527E+05 * exp( 1.26648 * lgt - 159689555.2 / rt );
	k[r193f] = 6.0000000000E+08 * exp( -41900000 / rt );
	k[r193b] = 1.0873190303E+05 * exp( 0.921243 * lgt + 35297691.18 / rt );
	k[r195f] = 3.1500000000E+15 * exp( -359000000 / rt );
	k[r195b] = 9.4467695801E+02 * exp( 1.73191 * lgt + 58477088.52 / rt );
	k[r196f] = 5.0000000000E+09 * exp( -6300000 / rt );
	k[r196b] = 4.1940112662E+09 * exp( -0.292037 * lgt - 110274697.9 / rt );
	k[r200f] = 7.9000000000E+09;
	k[r200b] = 8.9245699864E+05 * exp( 0.914116 * lgt - 57921366.61 / rt );
	k[r202f] = 4.0000000000E+09;
	k[r202b] = 3.4809292081E+10 * exp( -0.388353 * lgt - 168432325.9 / rt );
	k[r203] = 8.9600000000E+09 * exp( -35600000 / rt );
	k[r205f] = 9.6000000000E+13 * exp( -129800000 / rt );
	k[r205b] = 6.1932489287E+02 * exp( 1.76946 * lgt - 21962779.28 / rt );
	k[r206f] = 1.2500000000E+14 * exp( -154900000 / rt );
	k[r206b] = 1.7661829307E+09 * exp( 0.392108 * lgt - 11587778.91 / rt );
	k[r207f] = 1.0000000000E+09 * exp( -20900000 / rt );
	k[r207b] = 2.1661959347E+07 * exp( 0.393027 * lgt - 73784818.25 / rt );
	k[r255f] = 8.0000000000E+16 * exp( -307400000 / rt );
	k[r255b] = 7.8514484060E+07 * exp( 1.00703 * lgt - 23524837.58 / rt );
	k[r265f] = 6.5000000000E+09;
	k[r265b] = 1.2449931145E+08 * exp( 0.666947 * lgt - 58081854.44 / rt );
	k[r292] = 2.5000000000E+13 * exp( -120600000 / rt );
	k[r349] = 3.2000000000E+13 * exp( -118900000 / rt );
	k[r353] = 2.5000000000E+13 * exp( -125600000 / rt );
	k[r355] = 2.5000000000E+16 * exp( -297700000 / rt );
	k[r357] = 5.0000000000E+09;
	k[r385] = 5.0000000000E+09;
	k[rn192] = 2.5000000000E+13 * exp( -120600000 / rt );
	k[rn193] = 1.6000000000E+13 * exp( -118500000 / rt );
	k[rn194] = 4.0000000000E+13 * exp( -138200000 / rt );
	k[r477] = 2.0000000000E+11 * exp( -75800000 / rt );
	k[r478] = 3.0000000000E+11 * exp( -88400000 / rt );
	k[rn196] = 3.1600000000E+16 * exp( -339200000 / rt );
	k[rn197] = 7.3000000000E+04 * exp( 2 * lgt - 32200000 / rt );
	k[rn198] = 3.5000000000E+04 * exp( 2 * lgt - 20900000 / rt );
	k[rn201] = 1.0560000000E+07 * exp( 1.1 * lgt - 7600000 / rt );
	k[rn202] = 5.2000000000E+06 * exp( 1.3 * lgt - 2900000 / rt );
	k[rn209] = 1.7900000000E+10 * exp( -81200000 / rt );
	k[rn210] = 1.3400000000E+10 * exp( -71200000 / rt );
	k[rn217] = 5.5000000000E+10 * exp( -205200000 / rt );
	k[rn218] = 8.0000000000E+10 * exp( -199300000 / rt );
	k[ri66f] = 2.0000000000E+09;
	k[ri66b] = 1.7500000000E+15 * exp( -117300000 / rt );
	k[ri67f] = 2.0000000000E+09;
	k[ri67b] = 1.7500000000E+15 * exp( -117300000 / rt );
	k[ri68] = 6.0000000000E+11 * exp( -85600000 / rt );
	k[ri69] = 5.0000000000E+08;
	k[ri70] = 2.0000000000E+11 * exp( -71200000 / rt );
	k[ri71] = 1.0000000000E+09 * exp( -31400000 / rt );
	k[ri72] = 8.4000000000E+14 * exp( -180100000 / rt );
	k[ri73] = 2.0000000000E+13 * exp( -62800000 / rt );
	kTroe0 = 2.4770000000E+27 * exp( -4.76 * lgt - 10209900 / rt );
	kTroeInf = 1.2700000000E+13 * exp( -0.63 * lgt - 1602630 / rt );
	fcTroe = 0.217 * exp( -temp / 74 ) + 0.783 * exp( -temp / 2941 ) + 1 * exp( -6964 / temp );
	k[rnog52f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM41] );
	kTroe0 = 1.8571868185E+32 * exp( -4.6893 * lgt - 450255049.2 / rt );
	kTroeInf = 9.5221124725E+17 * exp( -0.559297 * lgt - 441647779.2 / rt );
	k[rnog52b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM41] );
	k[rnog8f] = 1.5000000000E+10;
	k[rnog8b] = 4.7604294830E+10 * exp( 0.16589 * lgt - 782848374.9 / rt );
	k[rnog9f] = 1.5000000000E+10;
	k[rnog9b] = 9.6339160312E+11 * exp( -0.407758 * lgt - 413951792.3 / rt );
	k[rnog51f] = 3.0000000000E+10;
	k[rnog51b] = 1.7462655565E+09 * exp( 0.187386 * lgt - 46785218.14 / rt );
	k[rnog62f] = 6.0000000000E+09;
	k[rnog62b] = 4.8004082511E+04 * exp( 1.28791 * lgt - 13046700.23 / rt );
	k[rnog67f] = 1.6000000000E+10;
	k[rnog67b] = 1.2801088670E+05 * exp( 1.28791 * lgt - 13046700.23 / rt );
	k[rnog79f] = 1.0000000000E+11;
	k[rnog79b] = 1.9141807991E+05 * exp( 1.6513 * lgt - 77419277.2 / rt );
	k[rnog94f] = 3.0000000000E+10;
	k[rnog94b] = 2.0884049407E+15 * exp( -0.827358 * lgt - 366848810 / rt );
	k[rnog97f] = 2.5010000000E+10;
	k[rnog97b] = 8.6213425898E+08 * exp( 0.491031 * lgt - 228311.6391 / rt );
	k[rnog142f] = 1.5000000000E+10 * exp( -2510640 / rt );
	k[rnog142b] = 6.4851138161E+09 * exp( -0.010174 * lgt - 39919365.22 / rt );
	k[rnog145f] = 1.2000000000E+10;
	k[rnog145b] = 8.0354055760E+08 * exp( 0.502665 * lgt - 775027396.7 / rt );
	kTroe0 = 2.7000000000E+32 * exp( -6.3 * lgt - 12971700 / rt );
	kTroeInf = 2.0000000000E+10;
	fcTroe = 0.8493 * exp( -temp / 134 ) + 0.1507 * exp( -temp / 2383 ) + 1 * exp( -7265 / temp );
	k[rnog147f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM40] );
	kTroe0 = 1.8146352403E+42 * exp( -7.58679 * lgt - 403377906.8 / rt );
	kTroeInf = 1.3441742521E+20 * exp( -1.28679 * lgt - 390406206.8 / rt );
	k[rnog147b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM40] );
	k[rnog148f] = 3.0000000000E+10;
	k[rnog148b] = 1.2970227632E+10 * exp( -0.010174 * lgt - 37408725.22 / rt );
	k[rnog149f] = 1.2000000000E+10 * exp( 2385110 / rt );
	k[rnog149b] = 6.6283580313E+16 * exp( -1.27611 * lgt - 304455672.1 / rt );
	k[rnog150f] = 1.6000000000E+10 * exp( 2385110 / rt );
	k[rnog150b] = 1.5305721048E+10 * exp( -0.31189 * lgt - 54883342.91 / rt );
	k[rnog151f] = 9.0000000000E+09;
	k[rnog151b] = 3.8910682896E+09 * exp( -0.010174 * lgt - 37408725.22 / rt );
	k[rnog152f] = 7.0000000000E+09;
	k[rnog152b] = 3.0263864475E+09 * exp( -0.010174 * lgt - 37408725.22 / rt );
	k[rnog153f] = 1.4000000000E+10;
	k[rnog153b] = 1.7512236837E+08 * exp( 0.452789 * lgt - 257631397.9 / rt );
	k[rnog154f] = 4.0000000000E+10 * exp( 2301420 / rt );
	k[rnog154b] = 4.0587833585E+08 * exp( 0.180011 * lgt - 75151759.06 / rt );
	k[rnog49f] = 1.1000000000E+11;
	k[rnog49b] = 4.4315662097E+10 * exp( 0.264487 * lgt - 94065509.37 / rt );
	k[rnog90f] = 5.0000000000E+10;
	k[rnog90b] = 1.1741375337E+13 * exp( -0.254665 * lgt - 649543612.4 / rt );
	k[rnog122f] = 5.8000000000E+10 * exp( -2410220 / rt );
	k[rnog122b] = 2.7699570052E+10 * exp( 0.178426 * lgt - 579675226.2 / rt );
	k[rnog123f] = 5.0000000000E+10;
	k[rnog123b] = 5.5838949387E+14 * exp( -0.79915 * lgt - 323297170.7 / rt );
	k[rnog124f] = 5.0000000000E+10;
	k[rnog124b] = 8.4573099810E+15 * exp( -1.12786 * lgt - 420597052.1 / rt );
	k[rnog178f] = 3.5000000000E+10 * exp( -1380850 / rt );
	k[rnog178b] = 5.6259447347E+10 * exp( 0.125962 * lgt - 315575625.1 / rt );
	k[rnog179f] = 2.6500000000E+09 * exp( -26780200 / rt );
	k[rnog179b] = 2.6890262107E+08 * exp( 0.0872414 * lgt - 158483692.6 / rt );
	k[rnog180f] = 7.3330000000E+10 * exp( -4686530 / rt );
	k[rnog180b] = 3.6587672165E+12 * exp( -0.34585 * lgt - 208668628.8 / rt );
	k[rnog181f] = 1.4000000000E+09 * exp( -45233400 / rt );
	k[rnog181b] = 1.9406281486E+05 * exp( 1.0043 * lgt - 372587652.2 / rt );
	k[rnog182f] = 2.9000000000E+10 * exp( -96868900 / rt );
	k[rnog182b] = 2.5376634987E+05 * exp( 0.96558 * lgt - 241731869.7 / rt );
	k[rnog183f] = 4.4000000000E+11 * exp( -79001500 / rt );
	k[rnog183b] = 1.2404036478E+05 * exp( 1.43739 * lgt - 334077146 / rt );
	k[rnog184f] = 2.0000000000E+09 * exp( -88123500 / rt );
	k[rnog184b] = 3.0946749225E+06 * exp( 0.690378 * lgt - 186136470.7 / rt );
	kTroe0 = 6.2000000000E+11 * exp( -234745000 / rt );
	kTroeInf = 1.3000000000E+11 * exp( -249474000 / rt );
	fcTroe = 1 * exp( -0 / temp );
	k[rnog185f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM40] );
	kTroe0 = 1.2726258212E+03 * exp( 1.12255 * lgt - 64282325.2 / rt );
	kTroeInf = 2.6684089798E+02 * exp( 1.12255 * lgt - 79011325.2 / rt );
	k[rnog185b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM40] );
	k[rnog186f] = 2.1100000000E+09 * exp( 2008510 / rt );
	k[rnog186b] = 6.9537903393E+10 * exp( -0.27051 * lgt - 37425650.36 / rt );
	k[rnog187f] = 1.0600000000E+14 * exp( -1.41 * lgt );
	k[rnog187b] = 2.6334252847E+21 * exp( -2.11268 * lgt - 307909805.9 / rt );
	k[rnog188f] = 3.9000000000E+09 * exp( 1004260 / rt );
	k[rnog188b] = 1.0601192504E+07 * exp( 0.584433 * lgt - 188902861.1 / rt );
	k[rnog189f] = 1.3200000000E+11 * exp( -1506390 / rt );
	k[rnog189b] = 7.2972649080E+05 * exp( 1.01752 * lgt - 119134904.9 / rt );
	k[rnog190f] = 5.0000000000E+10;
	k[rnog190b] = 2.7135459805E+12 * exp( -0.269294 * lgt - 298081254.4 / rt );
	k[rnog191f] = 3.2000000000E+10 * exp( -1380850 / rt );
	k[rnog191b] = 6.0396496297E+10 * exp( 0.107873 * lgt - 103025970.6 / rt );
	k[rnog192f] = 2.0000000000E+10;
	k[rnog192b] = 1.9472408866E+14 * exp( -0.735576 * lgt - 76066134.17 / rt );
	k[rnog193f] = 2.0000000000E+06 * exp( 1.2 * lgt );
	k[rnog193b] = 3.9162226410E+07 * exp( 1.21156 * lgt - 166102748.5 / rt );
	k[rnog194f] = 4.6100000000E+02 * exp( 2 * lgt - 27198600 / rt );
	k[rnog194b] = 9.1282321619E+03 * exp( 1.69752 * lgt - 30986127.98 / rt );
	k[rnog195f] = 1.2800000000E+03 * exp( 1.5 * lgt - 418440 / rt );
	k[rnog195b] = 1.4127757047E+02 * exp( 1.6638 * lgt - 226221088.2 / rt );
	k[rnog196f] = 1.5000000000E+10;
	k[rnog196b] = 1.3085365476E+12 * exp( -0.143332 * lgt - 612276029.5 / rt );
	k[rnog197f] = 2.0000000000E+10 * exp( -57954000 / rt );
	k[rnog197b] = 1.8769126812E+13 * exp( -0.639259 * lgt - 69562506.22 / rt );
	k[rnog198f] = 2.1600000000E+10 * exp( -0.23 * lgt );
	k[rnog198b] = 3.7765501403E+10 * exp( -0.027482 * lgt - 408293930.7 / rt );
	k[rnog199f] = 4.1600000000E+11 * exp( -0.45 * lgt );
	k[rnog199b] = 2.5800283392E+18 * exp( -1.68487 * lgt - 153218284.7 / rt );
	k[rnog200f] = 7.0000000000E+09;
	k[rnog200b] = 3.0351506078E+08 * exp( 0.249742 * lgt - 43295451.77 / rt );
	k[rnog201f] = 4.6000000000E+10;
	k[rnog201b] = 1.9419127898E+13 * exp( -0.485833 * lgt - 119361585.9 / rt );
	k[rnog202f] = 4.0000000000E+10 * exp( -15273100 / rt );
	k[rnog202b] = 3.0094690250E+09 * exp( 0.28106 * lgt - 66114516.74 / rt );
	k[rnog203f] = 9.0000000000E+04 * exp( 1.5 * lgt + 1924830 / rt );
	k[rnog203b] = 7.0250271522E+04 * exp( 1.68474 * lgt - 113374214.7 / rt );
	k[rnog204f] = 3.3000000000E+08;
	k[rnog204b] = 9.7359476541E+04 * exp( 0.28426 * lgt - 31189203.66 / rt );
	k[rnog205f] = 1.3000000000E+11 * exp( -0.11 * lgt - 20838300 / rt );
	k[rnog205b] = 3.8353733183E+07 * exp( 0.17426 * lgt - 52027503.66 / rt );
	k[rnog206f] = 5.0000000000E+09;
	k[rnog206b] = 2.2615508700E+09 * exp( 0.285179 * lgt - 227386243 / rt );
	k[rnog207f] = 2.5000000000E+10;
	k[rnog207b] = 1.0129910682E+09 * exp( 0.599102 * lgt - 456727524.4 / rt );
	k[rnog208f] = 7.0000000000E+10;
	k[rnog208b] = 1.6222662888E+09 * exp( 0.396584 * lgt - 48433593.76 / rt );
	k[rnog209f] = 5.0000000000E+10;
	k[rnog209b] = 3.5154691570E+09 * exp( 0.63042 * lgt - 464273489.4 / rt );
	k[rnog210f] = 2.0000000000E+10;
	k[rnog210b] = 1.4588777297E+10 * exp( 0.534103 * lgt - 528731117.4 / rt );
	k[rnog211f] = 2.5000000000E+10;
	k[rnog211b] = 5.5301114644E+11 * exp( 0.354963 * lgt - 471234352.8 / rt );
	k[rnog212f] = 8.9500000000E+13 * exp( -1.32 * lgt - 3096460 / rt );
	k[rnog212b] = 2.2051942734E+18 * exp( -1.47144 * lgt - 206619660.6 / rt );
	k[rnog213f] = 2.5000000000E+10;
	k[rnog213b] = 1.3935337940E+08 * exp( 0.466281 * lgt - 222015120.2 / rt );
	k[rnog214f] = 4.5000000000E+08 * exp( 0.72 * lgt - 2761710 / rt );
	k[rnog214b] = 4.3524890822E+06 * exp( 1.2176 * lgt - 232322795.2 / rt );
	k[rnog215f] = 1.3000000000E+04 * exp( 1.9 * lgt + 3975180 / rt );
	k[rnog215b] = 1.3045001708E+03 * exp( 2.30128 * lgt - 290043533.1 / rt );
	k[rnog216f] = 1.0000000000E+10 * exp( -54397200 / rt );
	k[rnog216b] = 6.2222613072E+08 * exp( 0.152358 * lgt - 47071038.78 / rt );
	k[rnog217f] = 7.7000000000E+10;
	k[rnog217b] = 2.7086445768E+11 * exp( 0.0459522 * lgt - 325717184 / rt );
	k[rnog218f] = 4.0000000000E+10;
	k[rnog218b] = 1.0405141573E+17 * exp( -1.42774 * lgt - 133036549 / rt );
	k[rnog219f] = 8.0000000000E+09 * exp( -31215700 / rt );
	k[rnog219b] = 2.4989183677E+12 * exp( -0.739777 * lgt - 56122905.4 / rt );
	k[rnog220f] = 6.1400000000E+09 * exp( 1841140 / rt );
	k[rnog220b] = 3.2482723916E+13 * exp( -0.994646 * lgt - 58916802.85 / rt );
	k[rnog221f] = 2.1000000000E+10 * exp( -19708500 / rt );
	k[rnog221b] = 6.8054522061E+13 * exp( -0.836093 * lgt - 109073333.4 / rt );
	k[rnog222f] = 2.3500000000E+10;
	k[rnog222b] = 1.5856022387E+06 * exp( 1.12784 * lgt - 396662733.7 / rt );
	k[rnog223f] = 5.4000000000E+10;
	k[rnog223b] = 6.7135614311E+04 * exp( 1.39713 * lgt - 98581479.32 / rt );
	k[rnog224f] = 2.5000000000E+09;
	k[rnog224b] = 1.2281886500E+03 * exp( 0.812998 * lgt + 28875587.07 / rt );
	k[rnog225f] = 2.0000000000E+10;
	k[rnog225b] = 2.1691211103E+06 * exp( 1.2538 * lgt - 710857508.8 / rt );
	k[rnog226f] = 2.0000000000E+09 * exp( -83688100 / rt );
	k[rnog226b] = 1.5273280298E+09 * exp( 0.280784 * lgt - 517289639.6 / rt );
	k[rnog227f] = 8.8000000000E+13 * exp( -0.5 * lgt - 200851000 / rt );
	k[rnog227b] = 8.6647223146E+05 * exp( 0.658847 * lgt + 32006685.84 / rt );
	k[rnog228f] = 2.8500000000E+14 * exp( -1.52 * lgt - 3096460 / rt );
	k[rnog228b] = 2.1975320261E+15 * exp( -1.35774 * lgt - 254896224 / rt );
	k[rnog229f] = 5.7000000000E+15 * exp( -2 * lgt - 3347520 / rt );
	k[rnog229b] = 6.8953360682E+16 * exp( -1.68049 * lgt - 619440342.1 / rt );
	k[rnog230f] = 1.0400000000E+26 * exp( -3.3 * lgt - 529746000 / rt );
	k[rnog230b] = 1.3466249242E+20 * exp( -2.81007 * lgt - 7296880.883 / rt );
	k[rnog231f] = 1.1070000000E+01 * exp( 2.64 * lgt - 20838300 / rt );
	k[rnog231b] = 5.1209466259E+03 * exp( 2.01704 * lgt - 56964050.72 / rt );
	k[rnog232f] = 2.7670000000E+00 * exp( 2.64 * lgt - 20838300 / rt );
	k[rnog232b] = 1.5913693710E-03 * exp( 3.41417 * lgt - 155545530 / rt );
	k[rnog233f] = 2.1340000000E+06 * exp( 1.58 * lgt - 111305000 / rt );
	k[rnog233b] = 3.7949760661E+02 * exp( 2.38478 * lgt - 14394201.67 / rt );
	k[rnog234f] = 1.1000000000E+03 * exp( 2.03 * lgt - 55945500 / rt );
	k[rnog234b] = 4.6102461468E+06 * exp( 1.33895 * lgt - 26339903.66 / rt );
	k[rnog235] = 4.4000000000E+00 * exp( 2.26 * lgt - 26780200 / rt );
	k[rnog236f] = 1.6000000000E-01 * exp( 2.56 * lgt - 37659600 / rt );
	k[rnog236b] = 2.1222650026E-03 * exp( 3.08443 * lgt - 129071378.3 / rt );
	k[rnog237f] = 1.4000000000E+20 * exp( -3.4 * lgt - 7950370 / rt );
	k[rnog237b] = 6.8832438020E+23 * exp( -3.52088 * lgt - 114250015.3 / rt );
	k[rnog238f] = 6.0000000000E+10 * exp( -1673760 / rt );
	k[rnog238b] = 1.9877558079E+07 * exp( 0.921814 * lgt - 330002010.3 / rt );
	k[rnog239f] = 6.3000000000E+10 * exp( -192566000 / rt );
	k[rnog239b] = 5.2438191940E+10 * exp( -0.0807297 * lgt + 1784445.46 / rt );
	k[rnog240f] = 1.4570000000E+05 * exp( 1.1 * lgt - 85361800 / rt );
	k[rnog240b] = 1.5833209632E+08 * exp( 0.447664 * lgt - 74441697.26 / rt );
	kTroe0 = 1.3000000000E+19 * exp( -3.16 * lgt - 3096460 / rt );
	kTroeInf = 3.1000000000E+09 * exp( 0.15 * lgt );
	fcTroe = 0.333 * exp( -temp / 235 ) + 0.667 * exp( -temp / 2117 ) + 1 * exp( -4536 / temp );
	k[rnog241f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM40] );
	kTroe0 = 6.3236286273E+28 * exp( -4.78086 * lgt - 138987012 / rt );
	kTroeInf = 1.5079422111E+19 * exp( -1.47086 * lgt - 135890552 / rt );
	k[rnog241b] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, M[mM40] );
	k[rnog242f] = 1.0000000000E+10 * exp( -309646000 / rt );
	k[rnog242b] = 7.7519387589E+11 * exp( -0.562649 * lgt - 206457269.6 / rt );
	k[rnog243f] = 1.0000000000E+08 * exp( -271986000 / rt );
	k[rnog243b] = 3.3514803431E+09 * exp( -0.572823 * lgt - 206205994.8 / rt );
	k[rnog244f] = 1.9000000000E+10;
	k[rnog244b] = 2.5420739351E+10 * exp( 0.0452326 * lgt - 119844329.6 / rt );
	k[rnog245f] = 2.9000000000E+10;
	k[rnog245b] = 1.3648781173E+11 * exp( 0.0911848 * lgt - 445561513.6 / rt );
	k[rnog246f] = 5.0000000000E+10;
	k[rnog246b] = 8.7338721798E+13 * exp( -0.526374 * lgt - 303274672.3 / rt );
	k[rnog247f] = 2.0000000000E+10;
	k[rnog247b] = 1.6161045445E+16 * exp( -1.14934 * lgt - 339400423 / rt );
	k[rnog248f] = 3.0000000000E+10;
	k[rnog248b] = 1.1511671958E+12 * exp( -0.217977 * lgt - 170730440.4 / rt );
	k[rnog249f] = 3.1000000000E+14 * exp( -1.38 * lgt - 5314190 / rt );
	k[rnog249b] = 1.7609416243E+20 * exp( -2.43118 * lgt - 280813793.9 / rt );
	k[rnog250f] = 2.9000000000E+11 * exp( -0.69 * lgt - 3180150 / rt );
	k[rnog250b] = 3.9305184111E+13 * exp( -1.05013 * lgt - 308285350.3 / rt );
	k[rnog251f] = 3.8000000000E+10 * exp( -0.36 * lgt - 2426950 / rt );
	k[rnog251b] = 2.1585736040E+16 * exp( -1.41118 * lgt - 277926553.9 / rt );
	k[rnog252f] = 3.1000000000E+14 * exp( -1.38 * lgt - 5314190 / rt );
	k[rnog252b] = 7.6132712379E+19 * exp( -2.44135 * lgt - 318222519.1 / rt );
	k[rnog253f] = 2.9000000000E+11 * exp( -0.69 * lgt - 3180150 / rt );
	k[rnog253b] = 1.6993239501E+13 * exp( -1.06031 * lgt - 345694075.5 / rt );
	k[rnog254f] = 3.8000000000E+10 * exp( -0.36 * lgt - 2426950 / rt );
	k[rnog254b] = 9.3323970013E+15 * exp( -1.42135 * lgt - 315335279.1 / rt );
	k[rnog255f] = 9.6000000000E+10 * exp( -120511000 / rt );
	k[rnog255b] = 1.9391465190E+11 * exp( 0.120725 * lgt - 463253237.1 / rt );
	k[rnog256f] = 1.0000000000E+09 * exp( -91010800 / rt );
	k[rnog256b] = 4.0167971473E+09 * exp( -0.249996 * lgt - 42510768.71 / rt );
	k[rnog257f] = 2.2000000000E+10;
	k[rnog257b] = 2.4658435595E+02 * exp( 1.59936 * lgt - 600172604.8 / rt );
	k[rnog258f] = 2.0000000000E+09;
	k[rnog258b] = 3.0617576551E+06 * exp( 0.937516 * lgt - 482709764.8 / rt );
	k[rnog259f] = 1.2000000000E+10;
	k[rnog259b] = 9.6056012366E+00 * exp( 1.49012 * lgt - 166543380.9 / rt );
	k[rnog260f] = 1.2000000000E+10;
	k[rnog260b] = 4.7231146305E+03 * exp( 1.05703 * lgt - 238821987.1 / rt );
	k[rnog261f] = 1.0000000000E+11;
	k[rnog261b] = 3.6388428729E+04 * exp( 1.76946 * lgt - 287817240.8 / rt );
	k[rnog262f] = 9.8000000000E+04 * exp( 1.41 * lgt - 35567400 / rt );
	k[rnog262b] = 7.4840443777E+04 * exp( 1.59507 * lgt - 309097638.5 / rt );
	k[rnog263f] = 1.5000000000E+05 * exp( 1.57 * lgt - 184114000 / rt );
	k[rnog263b] = 2.0040629117E+02 * exp( 2.29964 * lgt - 424492960.6 / rt );
	k[rnog264f] = 2.2000000000E+03 * exp( 2.11 * lgt - 47702200 / rt );
	k[rnog264b] = 2.4282527034E+02 * exp( 2.17808 * lgt - 113433547.1 / rt );
	k[rnog265f] = 2.2500000000E+04 * exp( 1.7 * lgt - 15900700 / rt );
	k[rnog265b] = 7.1208316678E-02 * exp( 2.91548 * lgt - 136918074.6 / rt );
	k[rnog266f] = 1.0500000000E+02 * exp( 2.5 * lgt - 55652600 / rt );
	k[rnog266b] = 2.0109819819E+01 * exp( 2.5994 * lgt - 128929912 / rt );
	k[rnog267f] = 4.6500000000E+09 * exp( -28663200 / rt );
	k[rnog267b] = 9.2394781882E+09 * exp( 0.00308589 * lgt - 166398140 / rt );
	k[rnog268f] = 1.5500000000E+09 * exp( -28663200 / rt );
	k[rnog268b] = 2.7299819172E+10 * exp( -0.0646704 * lgt - 258897986.7 / rt );
	k[rnog269f] = 1.1800000000E+13 * exp( -354503000 / rt );
	k[rnog269b] = 1.1789926501E+04 * exp( 1.15038 * lgt - 93277505.59 / rt );
	k[rnog270f] = 2.1000000000E+12 * exp( -0.69 * lgt - 11925600 / rt );
	k[rnog270b] = 2.1000000000E+12 * exp( -0.69 * lgt - 11925600 / rt );
	k[rnog271f] = 2.7000000000E+08 * exp( 0.18 * lgt - 8870940 / rt );
	k[rnog271b] = 6.4421722950E+04 * exp( 0.871047 * lgt - 38476536.34 / rt );
	k[rnog272f] = 1.7000000000E+11 * exp( -0.75 * lgt - 12092900 / rt );
	k[rnog272b] = 5.3801839268E+05 * exp( 0.465477 * lgt - 133110274.6 / rt );
	k[rnog273] = 2.0000000000E+04 * exp( 2 * lgt - 8368810 / rt );
	k[rnog274f] = 2.3500000000E+10;
	k[rnog274b] = 1.1047408842E+10 * exp( 0.58995 * lgt - 390327606.3 / rt );
	k[rnog275f] = 6.1000000000E+11 * exp( -0.31 * lgt - 1213480 / rt );
	k[rnog275b] = 1.2225393001E+14 * exp( -0.905846 * lgt - 156695547.5 / rt );
	k[rnog276f] = 3.7000000000E+09 * exp( 0.15 * lgt + 376596 / rt );
	k[rnog276b] = 3.5943353873E+10 * exp( 0.0211918 * lgt - 481890111.9 / rt );
	k[rnog277f] = 5.4000000000E+02 * exp( 2.4 * lgt - 41488400 / rt );
	k[rnog277b] = 2.9644337604E+00 * exp( 2.75769 * lgt - 18457275.71 / rt );
	k[rnog278f] = 5.0000000000E+04 * exp( 1.6 * lgt - 3996110 / rt );
	k[rnog278b] = 2.8476958769E+03 * exp( 1.86137 * lgt - 45422613.66 / rt );
	k[rnog279f] = 9.4000000000E+03 * exp( 1.94 * lgt - 27031200 / rt );
	k[rnog279b] = 2.9739123479E+01 * exp( 2.26637 * lgt + 3545889.261 / rt );
	k[rnog280f] = 4.0950565564E+12 * exp( -0.2 * lgt );
	k[rnog280b] = 9.7539466006E+10 * exp( 0.0206437 * lgt - 104810126.3 / rt );


	w[r1f] = k[r1f] * c[sO2] * c[sH];
	w[r1b] = k[r1b] * c[sO] * c[sOH];
	w[r2f] = k[r2f] * c[sH2] * c[sO];
	w[r2b] = k[r2b] * c[sH] * c[sOH];
	w[r3f] = k[r3f] * c[sH2] * c[sOH];
	w[r3b] = k[r3b] * c[sH] * c[sH2O];
	w[r4f] = k[r4f] * c[sOH] * c[sOH];
	w[r4b] = k[r4b] * c[sO] * c[sH2O];
	w[r5f] = k[r5f] * c[sH] * c[sH] * M[mM1];
	w[r5b] = k[r5b] * c[sH2] * M[mM1];
	w[r6f] = k[r6f] * c[sO] * c[sO] * M[mM1];
	w[r6b] = k[r6b] * c[sO2] * M[mM1];
	w[r7f] = k[r7f] * c[sH] * c[sOH] * M[mM1];
	w[r7b] = k[r7b] * c[sH2O] * M[mM1];
	w[r8f] = k[r8f] * c[sH] * c[sO2] * M[mM1];
	w[r8b] = k[r8b] * c[sHO2] * M[mM1];
	w[r9f] = k[r9f] * c[sHO2] * c[sH];
	w[r9b] = k[r9b] * c[sOH] * c[sOH];
	w[r10f] = k[r10f] * c[sHO2] * c[sH];
	w[r10b] = k[r10b] * c[sO2] * c[sH2];
	w[r11f] = k[r11f] * c[sHO2] * c[sH];
	w[r11b] = k[r11b] * c[sO] * c[sH2O];
	w[r12f] = k[r12f] * c[sHO2] * c[sO];
	w[r12b] = k[r12b] * c[sO2] * c[sOH];
	w[r13f] = k[r13f] * c[sHO2] * c[sOH];
	w[r13b] = k[r13b] * c[sO2] * c[sH2O];
	w[r14f] = k[r14f] * c[sHO2] * c[sHO2];
	w[r14b] = k[r14b] * c[sO2] * c[sH2O2];
	w[r15f] = k[r15f] * c[sOH] * c[sOH] * M[mM1];
	w[r15b] = k[r15b] * c[sH2O2] * M[mM1];
	w[r19f] = k[r19f] * c[sH2O2] * c[sOH];
	w[r19b] = k[r19b] * c[sHO2] * c[sH2O];
	w[r20f] = k[r20f] * c[sCO] * c[sOH];
	w[r20b] = k[r20b] * c[sH] * c[sCO2];
	w[r21f] = k[r21f] * c[sCO] * c[sHO2];
	w[r21b] = k[r21b] * c[sOH] * c[sCO2];
	w[r22f] = k[r22f] * c[sCO] * c[sO] * M[mM1];
	w[r22b] = k[r22b] * c[sCO2] * M[mM1];
	w[r25] = k[r25] * c[sCH] * c[sO2];
	w[r26] = k[r26] * c[sCH] * c[sCO2];
	w[r27f] = k[r27f] * c[sCH] * c[sH2O];
	w[r27b] = k[r27b] * c[sCH2OH];
	w[r28f] = k[r28f] * c[sHCO] * M[mM1];
	w[r28b] = k[r28b] * c[sH] * c[sCO] * M[mM1];
	w[r33f] = k[r33f] * c[sHCO] * c[sO2];
	w[r33b] = k[r33b] * c[sHO2] * c[sCO];
	w[r35f] = k[r35f] * c[s3XCH2] * c[sH];
	w[r35b] = k[r35b] * c[sH2] * c[sCH];
	w[r37f] = k[r37f] * c[s3XCH2] * c[s3XCH2];
	w[r37b] = k[r37b] * c[sH2] * c[sC2H2];
	w[r39f] = k[r39f] * c[s3XCH2] * c[sCH3];
	w[r39b] = k[r39b] * c[sH] * c[sC2H4];
	w[r40f] = k[r40f] * c[s3XCH2] * c[sO2];
	w[r40b] = k[r40b] * c[sH] * c[sOH] * c[sCO];
	w[r41f] = k[r41f] * c[s3XCH2] * c[sO2];
	w[r41b] = k[r41b] * c[sH2] * c[sCO2];
	w[r42f] = k[r42f] * c[s1XCH2] * M[mM1];
	w[r42b] = k[r42b] * c[s3XCH2] * M[mM1];
	w[r43f] = k[r43f] * c[s1XCH2] * c[sO2];
	w[r43b] = k[r43b] * c[sH] * c[sOH] * c[sCO];
	w[r44f] = k[r44f] * c[s1XCH2] * c[sH2];
	w[r44b] = k[r44b] * c[sH] * c[sCH3];
	w[r45f] = k[r45f] * c[sCH2O] * M[mM1];
	w[r45b] = k[r45b] * c[sH] * c[sHCO] * M[mM1];
	w[r46f] = k[r46f] * c[sCH2O] * c[sH];
	w[r46b] = k[r46b] * c[sH2] * c[sHCO];
	w[r47f] = k[r47f] * c[sCH2O] * c[sO];
	w[r47b] = k[r47b] * c[sOH] * c[sHCO];
	w[r48f] = k[r48f] * c[sCH2O] * c[sOH];
	w[r48b] = k[r48b] * c[sH2O] * c[sHCO];
	w[r49f] = k[r49f] * c[sCH2O] * c[sHO2];
	w[r49b] = k[r49b] * c[sH2O2] * c[sHCO];
	w[r53f] = k[r53f] * c[sCH3] * c[sO];
	w[r53b] = k[r53b] * c[sH] * c[sCH2O];
	w[ra34f] = k[ra34f] * c[sCH3] * c[sH];
	w[ra34b] = k[ra34b] * c[sCH4];
	w[r55] = k[r55] * c[sCH3] * c[sOH];
	w[r57] = k[r57] * c[sCH3] * c[sO2];
	w[r58f] = k[r58f] * c[sCH3] * c[sHO2];
	w[r58b] = k[r58b] * c[sOH] * c[sCH2OH];
	w[r59f] = k[r59f] * c[sCH3] * c[sHO2];
	w[r59b] = k[r59b] * c[sO2] * c[sCH4];
	w[r60] = k[r60] * c[sCH3] * c[sCH3];
	w[r36f] = k[r36f] * c[sCH3] * c[sCH3];
	w[r36b] = k[r36b] * c[sC2H6];
	w[r62f] = k[r62f] * c[sCH2OH] * M[mM1];
	w[r62b] = k[r62b] * c[sH] * c[sCH2O] * M[mM1];
	w[r63f] = k[r63f] * c[sCH2OH] * c[sH];
	w[r63b] = k[r63b] * c[sH2] * c[sCH2O];
	w[r64f] = k[r64f] * c[sCH2OH] * c[sO2];
	w[r64b] = k[r64b] * c[sHO2] * c[sCH2O];
	w[r69f] = k[r69f] * c[sCH2OH] * M[mM1];
	w[r69b] = k[r69b] * c[sH] * c[sCH2O] * M[mM1];
	w[r70f] = k[r70f] * c[sCH2OH] * c[sH];
	w[r70b] = k[r70b] * c[sH2] * c[sCH2O];
	w[r71f] = k[r71f] * c[sCH2OH] * c[sO2];
	w[r71b] = k[r71b] * c[sHO2] * c[sCH2O];
	w[r84f] = k[r84f] * c[sCH4] * c[sH];
	w[r84b] = k[r84b] * c[sCH3] * c[sH2];
	w[r86f] = k[r86f] * c[sCH4] * c[sOH];
	w[r86b] = k[r86b] * c[sCH3] * c[sH2O];
	w[r105f] = k[r105f] * c[sHCCO] * c[sH];
	w[r105b] = k[r105b] * c[sCO] * c[s3XCH2];
	w[r106] = k[r106] * c[sHCCO] * c[sO];
	w[r109f] = k[r109f] * c[sC2H2] * c[sO2];
	w[r109b] = k[r109b] * c[sOH] * c[sHCCO];
	w[r111f] = k[r111f] * c[sC2H2] * c[sO];
	w[r111b] = k[r111b] * c[sCO] * c[s3XCH2];
	w[r112f] = k[r112f] * c[sC2H2] * c[sO];
	w[r112b] = k[r112b] * c[sH] * c[sHCCO];
	w[ra51f] = k[ra51f] * c[sC2H3];
	w[ra51b] = k[ra51b] * c[sH] * c[sC2H2];
	w[ra125] = k[ra125] * c[sC2H3] * c[sO2];
	w[r129f] = k[r129f] * c[sC2H4] * M[mM1];
	w[r129b] = k[r129b] * c[sH2] * c[sC2H2] * M[mM1];
	w[r131f] = k[r131f] * c[sC2H4] * c[sH];
	w[r131b] = k[r131b] * c[sH2] * c[sC2H3];
	w[r134f] = k[r134f] * c[sC2H4] * c[sOH];
	w[r134b] = k[r134b] * c[sH2O] * c[sC2H3];
	w[ra58f] = k[ra58f] * c[sC2H5];
	w[ra58b] = k[ra58b] * c[sH] * c[sC2H4];
	w[r146f] = k[r146f] * c[sC2H5] * c[sH];
	w[r146b] = k[r146b] * c[sCH3] * c[sCH3];
	w[r149f] = k[r149f] * c[sC2H5] * c[sO2];
	w[r149b] = k[r149b] * c[sHO2] * c[sC2H4];
	w[r164] = k[r164] * c[sC2H6] * c[sH];
	w[r166] = k[r166] * c[sC2H6] * c[sOH];
	w[r170] = k[r170] * c[sC2H6] * c[sCH3];
	w[r189] = k[r189] * c[sC3H4] * c[sOH];
	w[r190] = k[r190] * c[sC3H4] * c[sOH];
	w[r191f] = k[r191f] * c[sC3H5];
	w[r191b] = k[r191b] * c[sH] * c[sC3H4];
	w[r192f] = k[r192f] * c[sC3H5] * c[sH];
	w[r192b] = k[r192b] * c[sH2] * c[sC3H4];
	w[r193f] = k[r193f] * c[sC3H5] * c[sO2];
	w[r193b] = k[r193b] * c[sHO2] * c[sC3H4];
	w[r195f] = k[r195f] * c[sC3H6];
	w[r195b] = k[r195b] * c[sCH3] * c[sC2H3];
	w[r196f] = k[r196f] * c[sC3H6] * c[sH];
	w[r196b] = k[r196b] * c[sH2] * c[sC3H5];
	w[r200f] = k[r200f] * c[sC3H6] * c[sOH];
	w[r200b] = k[r200b] * c[sCH2O] * c[sC2H5];
	w[r202f] = k[r202f] * c[sC3H6] * c[sOH];
	w[r202b] = k[r202b] * c[sH2O] * c[sC3H5];
	w[r203] = k[r203] * c[sC3H6] * c[sCH3];
	w[r205f] = k[r205f] * c[sNXC3H7];
	w[r205b] = k[r205b] * c[sC2H4] * c[sCH3];
	w[r206f] = k[r206f] * c[sNXC3H7];
	w[r206b] = k[r206b] * c[sC3H6] * c[sH];
	w[r207f] = k[r207f] * c[sNXC3H7] * c[sO2];
	w[r207b] = k[r207b] * c[sHO2] * c[sC3H6];
	w[r255f] = k[r255f] * c[s1XC4H8];
	w[r255b] = k[r255b] * c[sCH3] * c[sC3H5];
	w[r265f] = k[r265f] * c[s1XC4H8] * c[sOH];
	w[r265b] = k[r265b] * c[sCH2O] * c[sNXC3H7];
	w[r292] = k[r292] * c[sPXC4H9];
	w[r349] = k[r349] * c[s1XC5H11];
	w[r353] = k[r353] * c[sC6H11];
	w[r355] = k[r355] * c[s1XC6H12];
	w[r357] = k[r357] * c[s1XC6H12] * c[sH];
	w[r385] = k[r385] * c[s1XC6H12] * c[sOH];
	w[rn192] = k[rn192] * c[s1XC7H15];
	w[rn193] = k[rn193] * c[s2XC7H15];
	w[rn194] = k[rn194] * c[s2XC7H15];
	w[r477] = k[r477] * c[s1XC7H15];
	w[r478] = k[r478] * c[s2XC7H15];
	w[rn196] = k[rn196] * c[sNXC7H16];
	w[rn197] = k[rn197] * c[sNXC7H16] * c[sH];
	w[rn198] = k[rn198] * c[sNXC7H16] * c[sH];
	w[rn201] = k[rn201] * c[sNXC7H16] * c[sOH];
	w[rn202] = k[rn202] * c[sNXC7H16] * c[sOH];
	w[rn209] = k[rn209] * c[sNXC7H16] * c[sHO2];
	w[rn210] = k[rn210] * c[sNXC7H16] * c[sHO2];
	w[rn217] = k[rn217] * c[sNXC7H16] * c[sO2];
	w[rn218] = k[rn218] * c[sNXC7H16] * c[sO2];
	w[ri66f] = k[ri66f] * c[s1XC7H15] * c[sO2];
	w[ri66b] = k[ri66b] * c[s1C7H15O2XC7H15O2];
	w[ri67f] = k[ri67f] * c[s2XC7H15] * c[sO2];
	w[ri67b] = k[ri67b] * c[s1C7H15O2XC7H15O2];
	w[ri68] = k[ri68] * c[s1C7H15O2XC7H15O2];
	w[ri69] = k[ri69] * c[s1HEOOHX2XC7H15O2] * c[sO2];
	w[ri70] = k[ri70] * c[s2OO7OOH1XO2C7H15O2];
	w[ri71] = k[ri71] * c[s1X2HEOOHXO2C7H14OOH];
	w[ri72] = k[ri72] * c[sOC7OOHXOC7H13OOH];
	w[ri73] = k[ri73] * c[s1X2HEPOOXOC7H13O];
	w[rnog52f] = k[rnog52f] * c[sH] * c[sCH3];
	w[rnog52b] = k[rnog52b] * c[sCH4];
	w[rnog8f] = k[rnog8f] * c[sO] * c[s1XCH2];
	w[rnog8b] = k[rnog8b] * c[sCO] * c[sH2];
	w[rnog9f] = k[rnog9f] * c[sO] * c[s1XCH2];
	w[rnog9b] = k[rnog9b] * c[sHCO] * c[sH];
	w[rnog51f] = k[rnog51f] * c[sH] * c[s1XCH2];
	w[rnog51b] = k[rnog51b] * c[sH2] * c[sCH];
	w[rnog62f] = k[rnog62f] * c[sH] * c[sCH2OH];
	w[rnog62b] = k[rnog62b] * c[sH2O] * c[s1XCH2];
	w[rnog67f] = k[rnog67f] * c[sH] * c[sCH2OH];
	w[rnog67b] = k[rnog67b] * c[sH2O] * c[s1XCH2];
	w[rnog79f] = k[rnog79f] * c[sH] * c[sHCCO];
	w[rnog79b] = k[rnog79b] * c[sCO] * c[s1XCH2];
	w[rnog94f] = k[rnog94f] * c[sOH] * c[s1XCH2];
	w[rnog94b] = k[rnog94b] * c[sCH2O] * c[sH];
	w[rnog97f] = k[rnog97f] * c[sOH] * c[sCH3];
	w[rnog97b] = k[rnog97b] * c[sH2O] * c[s1XCH2];
	w[rnog142f] = k[rnog142f] * c[s1XCH2] * c[sN2];
	w[rnog142b] = k[rnog142b] * c[sN2] * c[sCH2];
	w[rnog145f] = k[rnog145f] * c[s1XCH2] * c[sO2];
	w[rnog145b] = k[rnog145b] * c[sH2O] * c[sCO];
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
	w[rnog249b] = k[rnog249b] * c[sHOCN] * c[sH];
	w[rnog250f] = k[rnog250f] * c[sCH2] * c[sNO];
	w[rnog250b] = k[rnog250b] * c[sHCN] * c[sOH];
	w[rnog251f] = k[rnog251f] * c[sCH2] * c[sNO];
	w[rnog251b] = k[rnog251b] * c[sHOCN] * c[sH];
	w[rnog252f] = k[rnog252f] * c[s1XCH2] * c[sNO];
	w[rnog252b] = k[rnog252b] * c[sHOCN] * c[sH];
	w[rnog253f] = k[rnog253f] * c[s1XCH2] * c[sNO];
	w[rnog253b] = k[rnog253b] * c[sHCN] * c[sOH];
	w[rnog254f] = k[rnog254f] * c[s1XCH2] * c[sNO];
	w[rnog254b] = k[rnog254b] * c[sHOCN] * c[sH];
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
	w[rnog262f] = k[rnog262f] * c[sHOCN] * c[sO];
	w[rnog262b] = k[rnog262b] * c[sCO2] * c[sNH];
	w[rnog263f] = k[rnog263f] * c[sHOCN] * c[sO];
	w[rnog263b] = k[rnog263b] * c[sCO] * c[sHNO];
	w[rnog264f] = k[rnog264f] * c[sHOCN] * c[sO];
	w[rnog264b] = k[rnog264b] * c[sOH] * c[sNCO];
	w[rnog265f] = k[rnog265f] * c[sHOCN] * c[sH];
	w[rnog265b] = k[rnog265b] * c[sCO] * c[sNH2];
	w[rnog266f] = k[rnog266f] * c[sHOCN] * c[sH];
	w[rnog266b] = k[rnog266b] * c[sNCO] * c[sH2];
	w[rnog267f] = k[rnog267f] * c[sHOCN] * c[sOH];
	w[rnog267b] = k[rnog267b] * c[sH2O] * c[sNCO];
	w[rnog268f] = k[rnog268f] * c[sHOCN] * c[sOH];
	w[rnog268b] = k[rnog268b] * c[sCO2] * c[sNH2];
	w[rnog269f] = k[rnog269f] * c[sHOCN] * M[mM40];
	w[rnog269b] = k[rnog269b] * c[sCO] * c[sNH] * M[mM40];
	w[rnog270f] = k[rnog270f] * c[sHOCN] * c[sH];
	w[rnog270b] = k[rnog270b] * c[sHOCN] * c[sH];
	w[rnog271f] = k[rnog271f] * c[sHOCN] * c[sH];
	w[rnog271b] = k[rnog271b] * c[sHCN] * c[sOH];
	w[rnog272f] = k[rnog272f] * c[sHOCN] * c[sH];
	w[rnog272b] = k[rnog272b] * c[sCO] * c[sNH2];
	w[rnog273] = k[rnog273] * c[sHOCN] * c[sH];
	w[rnog274f] = k[rnog274f] * c[sHCCO] * c[sNO];
	w[rnog274b] = k[rnog274b] * c[sCO] * c[sHOCN];
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
	w[rnog280b] = k[rnog280b] * c[sNO] * c[sCH2OH];


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

	cdot[sO2] = - w[r1f] + w[r1b] + w[r6f] - w[r6b]
		- w[r8f] + w[r8b] + w[r10f] - w[r10b]
		+ w[r12f] - w[r12b] + w[r13f] - w[r13b]
		+ w[r14f] - w[r14b] - w[r25] - w[r33f]
		+ w[r33b] - w[r40f] + w[r40b] - w[r41f]
		+ w[r41b] - w[r43f] + w[r43b] - w[r57]
		+ w[r59f] - w[r59b] - w[r64f] + w[r64b]
		- w[r71f] + w[r71b] - w[r109f] + w[r109b]
		- w[ra125] - w[r149f] + w[r149b] - w[r193f]
		+ w[r193b] - w[r207f] + w[r207b] - w[rn217]
		- w[rn218] - w[ri66f] + w[ri66b] - w[ri67f]
		+ w[ri67b] - w[ri69] - w[rnog145f] + w[rnog145b]
		- w[rnog122f] + w[rnog122b] - w[rnog179f] + w[rnog179b]
		+ w[rnog181f] - w[rnog181b] + w[rnog188f] - w[rnog188b]
		- w[rnog194f] + w[rnog194b] - w[rnog195f] + w[rnog195b]
		- w[rnog206f] + w[rnog206b] - w[rnog216f] + w[rnog216b]
		- w[rnog220f] + w[rnog220b] - w[rnog226f] + w[rnog226b]
		- w[rnog259f] + w[rnog259b];

	cdot[sH] = - w[r1f] + w[r1b] + w[r2f] - w[r2b]
		+ w[r3f] - w[r3b] - 2 * w[r5f] + 2 * w[r5b]
		- w[r7f] + w[r7b] - w[r8f] + w[r8b]
		- w[r9f] + w[r9b] - w[r10f] + w[r10b]
		- w[r11f] + w[r11b] + w[r20f] - w[r20b]
		+ w[r28f] - w[r28b] - w[r35f] + w[r35b]
		+ w[r39f] - w[r39b] + w[r40f] - w[r40b]
		+ w[r43f] - w[r43b] + w[r44f] - w[r44b]
		+ w[r45f] - w[r45b] - w[r46f] + w[r46b]
		+ w[r53f] - w[r53b] - w[ra34f] + w[ra34b]
		+ w[r55] + w[r62f] - w[r62b] - w[r63f]
		+ w[r63b] + w[r69f] - w[r69b] - w[r70f]
		+ w[r70b] - w[r84f] + w[r84b] - w[r105f]
		+ w[r105b] + w[r106] + w[r112f] - w[r112b]
		+ w[ra51f] - w[ra51b] - w[r131f] + w[r131b]
		+ w[ra58f] - w[ra58b] - w[r146f] + w[r146b]
		- w[r164] + w[r191f] - w[r191b] - w[r192f]
		+ w[r192b] - w[r196f] + w[r196b] + w[r206f]
		- w[r206b] - w[r357] - w[rn197] - w[rn198]
		- w[rnog52f] + w[rnog52b] + w[rnog9f] - w[rnog9b]
		- w[rnog51f] + w[rnog51b] - w[rnog62f] + w[rnog62b]
		- w[rnog67f] + w[rnog67b] - w[rnog79f] + w[rnog79b]
		+ w[rnog94f] - w[rnog94b] + w[rnog149f] - w[rnog149b]
		- w[rnog49f] + w[rnog49b] + w[rnog90f] - w[rnog90b]
		+ w[rnog123f] - w[rnog123b] + w[rnog124f] - w[rnog124b]
		+ w[rnog180f] - w[rnog180b] - w[rnog183f] + w[rnog183b]
		- w[rnog189f] + w[rnog189b] + w[rnog190f] - w[rnog190b]
		- w[rnog191f] + w[rnog191b] + w[rnog192f] - w[rnog192b]
		+ w[rnog196f] - w[rnog196b] + w[rnog199f] - w[rnog199b]
		+ w[rnog201f] - w[rnog201b] - w[rnog202f] + w[rnog202b]
		+ w[rnog204f] - w[rnog204b] + w[rnog205f] - w[rnog205b]
		- w[rnog209f] + w[rnog209b] - w[rnog212f] + w[rnog212b]
		- w[rnog214f] + w[rnog214b] + w[rnog218f] - w[rnog218b]
		+ w[rnog221f] - w[rnog221b] - w[rnog223f] + w[rnog223b]
		+ w[rnog224f] - w[rnog224b] + w[rnog230f] - w[rnog230b]
		+ w[rnog231f] - w[rnog231b] + w[rnog234f] - w[rnog234b]
		+ w[rnog235] - w[rnog237f] + w[rnog237b] + w[rnog247f]
		- w[rnog247b] + w[rnog249f] - w[rnog249b] + w[rnog251f]
		- w[rnog251b] + w[rnog252f] - w[rnog252b] + w[rnog254f]
		- w[rnog254b] + w[rnog257f] - w[rnog257b] + w[rnog260f]
		- w[rnog260b] - w[rnog261f] + w[rnog261b] - w[rnog265f]
		+ w[rnog265b] - w[rnog266f] + w[rnog266b] - w[rnog270f]
		+ w[rnog270f] - w[rnog270b] + w[rnog270b] - w[rnog271f]
		+ w[rnog271b] - w[rnog272f] + w[rnog272b] - w[rnog273]
		+ w[rnog273] + w[rnog275f] - w[rnog275b] - w[rnog277f]
		+ w[rnog277b];

	cdot[sOH] = w[r1f] - w[r1b] + w[r2f] - w[r2b]
		- w[r3f] + w[r3b] - 2 * w[r4f] + 2 * w[r4b]
		- w[r7f] + w[r7b] + 2 * w[r9f] - 2 * w[r9b]
		+ w[r12f] - w[r12b] - w[r13f] + w[r13b]
		- 2 * w[r15f] + 2 * w[r15b] - w[r19f] + w[r19b]
		- w[r20f] + w[r20b] + w[r21f] - w[r21b]
		+ w[r40f] - w[r40b] + w[r43f] - w[r43b]
		+ w[r47f] - w[r47b] - w[r48f] + w[r48b]
		- w[r55] + w[r57] + w[r58f] - w[r58b]
		- w[r86f] + w[r86b] + w[r109f] - w[r109b]
		- w[r134f] + w[r134b] - w[r166] - w[r189]
		- w[r190] - w[r200f] + w[r200b] - w[r202f]
		+ w[r202b] - w[r265f] + w[r265b] - w[r385]
		- w[rn201] - w[rn202] + w[ri71] + w[ri72]
		- w[rnog94f] + w[rnog94b] - w[rnog97f] + w[rnog97b]
		- w[rnog90f] + w[rnog90b] - w[rnog180f] + w[rnog180b]
		+ w[rnog183f] - w[rnog183b] - w[rnog184f] + w[rnog184b]
		+ w[rnog186f] - w[rnog186b] + w[rnog189f] - w[rnog189b]
		- w[rnog192f] + w[rnog192b] - w[rnog193f] + w[rnog193b]
		+ w[rnog195f] - w[rnog195b] + w[rnog198f] - w[rnog198b]
		+ w[rnog200f] - w[rnog200b] - w[rnog203f] + w[rnog203b]
		+ w[rnog207f] - w[rnog207b] - w[rnog210f] + w[rnog210b]
		+ w[rnog213f] - w[rnog213b] - w[rnog215f] + w[rnog215b]
		- w[rnog218f] + w[rnog218b] + w[rnog219f] - w[rnog219b]
		- w[rnog224f] + w[rnog224b] + w[rnog233f] - w[rnog233b]
		- w[rnog234f] + w[rnog234b] - w[rnog235] - w[rnog236f]
		+ w[rnog236b] + w[rnog250f] - w[rnog250b] + w[rnog253f]
		- w[rnog253b] + w[rnog256f] - w[rnog256b] - w[rnog260f]
		+ w[rnog260b] + w[rnog264f] - w[rnog264b] - w[rnog267f]
		+ w[rnog267b] - w[rnog268f] + w[rnog268b] + w[rnog271f]
		- w[rnog271b] - w[rnog278f] + w[rnog278b] + w[rnog279f]
		- w[rnog279b];

	cdot[sO] = w[r1f] - w[r1b] - w[r2f] + w[r2b]
		+ w[r4f] - w[r4b] - 2 * w[r6f] + 2 * w[r6b]
		+ w[r11f] - w[r11b] - w[r12f] + w[r12b]
		- w[r22f] + w[r22b] + w[r25] - w[r47f]
		+ w[r47b] - w[r53f] + w[r53b] - w[r106]
		- w[r111f] + w[r111b] - w[r112f] + w[r112b]
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

	cdot[sH2] = - w[r2f] + w[r2b] - w[r3f] + w[r3b]
		+ w[r5f] - w[r5b] + w[r10f] - w[r10b]
		+ w[r35f] - w[r35b] + w[r37f] - w[r37b]
		+ w[r41f] - w[r41b] - w[r44f] + w[r44b]
		+ w[r46f] - w[r46b] + w[r60] + w[r63f]
		- w[r63b] + w[r70f] - w[r70b] + w[r84f]
		- w[r84b] + w[r129f] - w[r129b] + w[r131f]
		- w[r131b] + w[r164] + w[r192f] - w[r192b]
		+ w[r196f] - w[r196b] + w[r357] + w[rn197]
		+ w[rn198] + w[rnog8f] - w[rnog8b] + w[rnog51f]
		- w[rnog51b] + w[rnog49f] - w[rnog49b] + w[rnog191f]
		- w[rnog191b] + w[rnog197f] - w[rnog197b] + w[rnog202f]
		- w[rnog202b] + w[rnog209f] - w[rnog209b] + w[rnog214f]
		- w[rnog214b] - w[rnog221f] + w[rnog221b] + w[rnog266f]
		- w[rnog266b] + w[rnog276f] - w[rnog276b] + w[rnog277f]
		- w[rnog277b];

	cdot[sH2O] = w[r3f] - w[r3b] + w[r4f] - w[r4b]
		+ w[r7f] - w[r7b] + w[r11f] - w[r11b]
		+ w[r13f] - w[r13b] + w[r19f] - w[r19b]
		- w[r27f] + w[r27b] + w[r48f] - w[r48b]
		+ w[r86f] - w[r86b] + w[r134f] - w[r134b]
		+ w[r166] + w[r202f] - w[r202b] + w[r385]
		+ w[rn201] + w[rn202] + w[rnog62f] - w[rnog62b]
		+ w[rnog67f] - w[rnog67b] + w[rnog97f] - w[rnog97b]
		+ w[rnog145f] - w[rnog145b] - w[rnog147f] + w[rnog147b]
		- w[rnog148f] + w[rnog148f] - w[rnog148b] + w[rnog148b]
		+ w[rnog193f] - w[rnog193b] - w[rnog197f] + w[rnog197b]
		+ w[rnog203f] - w[rnog203b] + w[rnog210f] - w[rnog210b]
		+ w[rnog215f] - w[rnog215b] - w[rnog219f] + w[rnog219b]
		+ w[rnog255f] - w[rnog255b] + w[rnog267f] - w[rnog267b]
		+ w[rnog278f] - w[rnog278b];

	cdot[sHO2] = w[r8f] - w[r8b] - w[r9f] + w[r9b]
		- w[r10f] + w[r10b] - w[r11f] + w[r11b]
		- w[r12f] + w[r12b] - w[r13f] + w[r13b]
		- 2 * w[r14f] + 2 * w[r14b] + w[r19f] - w[r19b]
		- w[r21f] + w[r21b] + w[r33f] - w[r33b]
		- w[r49f] + w[r49b] - w[r58f] + w[r58b]
		- w[r59f] + w[r59b] + w[r64f] - w[r64b]
		+ w[r71f] - w[r71b] + w[r149f] - w[r149b]
		+ w[r193f] - w[r193b] + w[r207f] - w[r207b]
		- w[rn209] - w[rn210] + w[rn217] + w[rn218]
		+ w[rnog184f] - w[rnog184b] - w[rnog186f] + w[rnog186b]
		+ w[rnog206f] - w[rnog206b] + w[rnog216f] - w[rnog216b];

	cdot[sH2O2] = w[r14f] - w[r14b] + w[r15f] - w[r15b]
		- w[r19f] + w[r19b] + w[r49f] - w[r49b]
		+ w[rn209] + w[rn210];

	cdot[sCO] = - w[r20f] + w[r20b] - w[r21f] + w[r21b]
		- w[r22f] + w[r22b] + w[r26] + w[r28f]
		- w[r28b] + w[r33f] - w[r33b] + w[r40f]
		- w[r40b] + w[r43f] - w[r43b] + w[r105f]
		- w[r105b] + 2 * w[r106] + w[r111f] - w[r111b]
		+ w[ri73] + w[rnog8f] - w[rnog8b] + w[rnog79f]
		- w[rnog79b] + w[rnog145f] - w[rnog145b] - w[rnog151f]
		+ w[rnog151f] - w[rnog151b] + w[rnog151b] + w[rnog153f]
		- w[rnog153b] + w[rnog90f] - w[rnog90b] + w[rnog122f]
		- w[rnog122b] + w[rnog217f] - w[rnog217b] + w[rnog222f]
		- w[rnog222b] + w[rnog223f] - w[rnog223b] + w[rnog224f]
		- w[rnog224b] + w[rnog225f] - w[rnog225b] + w[rnog227f]
		- w[rnog227b] + w[rnog228f] - w[rnog228b] + w[rnog232f]
		- w[rnog232b] + w[rnog236f] - w[rnog236b] + w[rnog245f]
		- w[rnog245b] + w[rnog257f] - w[rnog257b] + w[rnog263f]
		- w[rnog263b] + w[rnog265f] - w[rnog265b] + w[rnog269f]
		- w[rnog269b] + w[rnog272f] - w[rnog272b] + w[rnog274f]
		- w[rnog274b];

	cdot[sCO2] = w[r20f] - w[r20b] + w[r21f] - w[r21b]
		+ w[r22f] - w[r22b] - w[r26] + w[r41f]
		- w[r41b] - w[rnog152f] + w[rnog152f] - w[rnog152b]
		+ w[rnog152b] - w[rnog153f] + w[rnog153b] + w[rnog226f]
		- w[rnog226b] + w[rnog229f] - w[rnog229b] + w[rnog262f]
		- w[rnog262b] + w[rnog268f] - w[rnog268b];

	cdot[sCH] = - w[r25] - w[r26] - w[r27f] + w[r27b]
		+ w[r35f] - w[r35b] + w[rnog51f] - w[rnog51b]
		- w[rnog49f] + w[rnog49b] - w[rnog240f] + w[rnog240b]
		- w[rnog241f] + w[rnog241b] - w[rnog246f] + w[rnog246b]
		- w[rnog247f] + w[rnog247b] - w[rnog248f] + w[rnog248b];

	cdot[sHCO] = w[r25] + w[r26] - w[r28f] + w[r28b]
		- w[r33f] + w[r33b] + w[r45f] - w[r45b]
		+ w[r46f] - w[r46b] + w[r47f] - w[r47b]
		+ w[r48f] - w[r48b] + w[r49f] - w[r49b]
		+ w[ra125] + w[r190] + w[rnog9f] - w[rnog9b]
		+ w[rnog248f] - w[rnog248b] + w[rnog259f] - w[rnog259b]
		+ w[rnog260f] - w[rnog260b];

	cdot[sCH2OH] = w[r27f] - w[r27b] + w[r55] + w[r58f]
		- w[r58b] - w[r62f] + w[r62b] - w[r63f]
		+ w[r63b] - w[r64f] + w[r64b] - w[r69f]
		+ w[r69b] - w[r70f] + w[r70b] - w[r71f]
		+ w[r71b] - w[rnog62f] + w[rnog62b] - w[rnog67f]
		+ w[rnog67b] + w[rnog280f] - w[rnog280b];

	cdot[s3XCH2] = - w[r35f] + w[r35b] - 2 * w[r37f] + 2 * w[r37b]
		- w[r39f] + w[r39b] - w[r40f] + w[r40b]
		- w[r41f] + w[r41b] + w[r42f] - w[r42b]
		+ w[r105f] - w[r105b] + w[r111f] - w[r111b];

	cdot[sC2H2] = w[r37f] - w[r37b] - w[r109f] + w[r109b]
		- w[r111f] + w[r111b] - w[r112f] + w[r112b]
		+ w[ra51f] - w[ra51b] + w[r129f] - w[r129b]
		+ w[rnog124f] - w[rnog124b];

	cdot[sCH3] = - w[r39f] + w[r39b] + w[r44f] - w[r44b]
		- w[r53f] + w[r53b] - w[ra34f] + w[ra34b]
		- w[r55] - w[r57] - w[r58f] + w[r58b]
		- w[r59f] + w[r59b] - 2 * w[r60] - 2 * w[r36f]
		+ 2 * w[r36b] + w[r84f] - w[r84b] + w[r86f]
		- w[r86b] + 2 * w[r146f] - 2 * w[r146b] - w[r170]
		+ w[r195f] - w[r195b] - w[r203] + w[r205f]
		- w[r205b] + w[r255f] - w[r255b] + w[rn194]
		- w[rnog52f] + w[rnog52b] - w[rnog97f] + w[rnog97b]
		- w[rnog149f] + w[rnog149b] + 2 * w[rnog150f] - 2 * w[rnog150b]
		+ w[rnog154f] - w[rnog154b] - w[rnog124f] + w[rnog124b]
		- w[rnog211f] + w[rnog211b] - w[rnog255f] + w[rnog255b]
		- w[rnog256f] + w[rnog256b] - w[rnog275f] + w[rnog275b]
		- w[rnog276f] + w[rnog276b] - w[rnog280f] + w[rnog280b];

	cdot[sC2H4] = w[r39f] - w[r39b] + w[r60] - w[r129f]
		+ w[r129b] - w[r131f] + w[r131b] - w[r134f]
		+ w[r134b] + w[ra58f] - w[ra58b] + w[r149f]
		- w[r149b] + w[r190] + w[r205f] - w[r205b]
		+ w[r292] + w[r349] + w[rn192] + w[rnog149f]
		- w[rnog149b];

	cdot[s1XCH2] = - w[r42f] + w[r42b] - w[r43f] + w[r43b]
		- w[r44f] + w[r44b] - w[rnog8f] + w[rnog8b]
		- w[rnog9f] + w[rnog9b] - w[rnog51f] + w[rnog51b]
		+ w[rnog62f] - w[rnog62b] + w[rnog67f] - w[rnog67b]
		+ w[rnog79f] - w[rnog79b] - w[rnog94f] + w[rnog94b]
		+ w[rnog97f] - w[rnog97b] - w[rnog142f] + w[rnog142b]
		- w[rnog145f] + w[rnog145b] - w[rnog147f] + w[rnog147b]
		- w[rnog148f] + w[rnog148b] - w[rnog149f] + w[rnog149b]
		- w[rnog150f] + w[rnog150b] - w[rnog151f] + w[rnog151b]
		- w[rnog152f] + w[rnog152b] - w[rnog153f] + w[rnog153b]
		- w[rnog154f] + w[rnog154b] - w[rnog243f] + w[rnog243b]
		- w[rnog252f] + w[rnog252b] - w[rnog253f] + w[rnog253b]
		- w[rnog254f] + w[rnog254b];

	cdot[sCH2O] = - w[r45f] + w[r45b] - w[r46f] + w[r46b]
		- w[r47f] + w[r47b] - w[r48f] + w[r48b]
		- w[r49f] + w[r49b] + w[r53f] - w[r53b]
		+ w[r57] + w[r62f] - w[r62b] + w[r63f]
		- w[r63b] + w[r64f] - w[r64b] + w[r69f]
		- w[r69b] + w[r70f] - w[r70b] + w[r71f]
		- w[r71b] + w[ra125] + w[r189] + w[r200f]
		- w[r200b] + w[r265f] - w[r265b] + w[ri73]
		+ w[rnog94f] - w[rnog94b] + w[rnog153f] - w[rnog153b];

	cdot[sCH4] = w[ra34f] - w[ra34b] + w[r59f] - w[r59b]
		- w[r84f] + w[r84b] - w[r86f] + w[r86b]
		+ w[r170] + w[r203] + w[rnog52f] - w[rnog52b]
		- w[rnog150f] + w[rnog150b] + w[rnog211f] - w[rnog211b];

	cdot[sC2H6] = w[r36f] - w[r36b] - w[r164] - w[r166]
		- w[r170] - w[rnog154f] + w[rnog154b];

	cdot[sHCCO] = - w[r105f] + w[r105b] - w[r106] + w[r109f]
		- w[r109b] + w[r112f] - w[r112b] - w[rnog79f]
		+ w[rnog79b] - w[rnog274f] + w[rnog274b];

	cdot[sC2H3] = - w[ra51f] + w[ra51b] - w[ra125] + w[r131f]
		- w[r131b] + w[r134f] - w[r134b] + w[r189]
		+ w[r195f] - w[r195b];

	cdot[sC2H5] = - w[ra58f] + w[ra58b] - w[r146f] + w[r146b]
		- w[r149f] + w[r149b] + w[r164] + w[r166]
		+ w[r170] + w[r200f] - w[r200b] + w[r292]
		+ w[rnog154f] - w[rnog154b];

	cdot[sC3H4] = - w[r189] - w[r190] + w[r191f] - w[r191b]
		+ w[r192f] - w[r192b] + w[r193f] - w[r193b];

	cdot[sC3H5] = - w[r191f] + w[r191b] - w[r192f] + w[r192b]
		- w[r193f] + w[r193b] + w[r196f] - w[r196b]
		+ w[r202f] - w[r202b] + w[r203] + w[r255f]
		- w[r255b] + w[r353] + w[r355];

	cdot[sC3H6] = - w[r195f] + w[r195b] - w[r196f] + w[r196b]
		- w[r200f] + w[r200b] - w[r202f] + w[r202b]
		- w[r203] + w[r206f] - w[r206b] + w[r207f]
		- w[r207b] + w[r353] + w[rn193];

	cdot[sNXC3H7] = - w[r205f] + w[r205b] - w[r206f] + w[r206b]
		- w[r207f] + w[r207b] + w[r265f] - w[r265b]
		+ w[r349] + w[r355] + w[rn196];

	cdot[s1XC4H8] = - w[r255f] + w[r255b] - w[r265f] + w[r265b];

	cdot[sPXC4H9] = - w[r292] + w[rn193] + w[rn196];

	cdot[s1XC5H11] = - w[r349] + w[rn192] + w[ri73];

	cdot[sC6H11] = - w[r353] + w[r357] + w[r385];

	cdot[s1XC6H12] = - w[r355] - w[r357] - w[r385] + w[rn194];

	cdot[s1XC7H15] = - w[rn192] - w[r477] + w[r478] + w[rn197]
		+ w[rn201] + w[rn209] + w[rn217] - w[ri66f]
		+ w[ri66b];

	cdot[s2XC7H15] = - w[rn193] - w[rn194] + w[r477] - w[r478]
		+ w[rn198] + w[rn202] + w[rn210] + w[rn218]
		- w[ri67f] + w[ri67b];

	cdot[sNXC7H16] = - w[rn196] - w[rn197] - w[rn198] - w[rn201]
		- w[rn202] - w[rn209] - w[rn210] - w[rn217]
		- w[rn218];

	cdot[s1C7H15O2XC7H15O2] = w[ri66f] - w[ri66b] + w[ri67f] - w[ri67b]
		- w[ri68];

	cdot[s1HEOOHX2XC7H15O2] = w[ri68] - w[ri69];

	cdot[s2OO7OOH1XO2C7H15O2] = w[ri69] - w[ri70];

	cdot[s1X2HEOOHXO2C7H14OOH] = w[ri70] - w[ri71];

	cdot[sOC7OOHXOC7H13OOH] = w[ri71] - w[ri72];

	cdot[s1X2HEPOOXOC7H13O] = w[ri72] - w[ri73];

	cdot[sCH2] = w[rnog142f] - w[rnog142b] + w[rnog148f] - w[rnog148b]
		+ w[rnog151f] - w[rnog151b] + w[rnog152f] - w[rnog152b]
		- w[rnog123f] + w[rnog123b] + w[rnog238f] - w[rnog238b]
		- w[rnog242f] + w[rnog242b] - w[rnog249f] + w[rnog249b]
		- w[rnog250f] + w[rnog250b] - w[rnog251f] + w[rnog251b]
		+ w[rnog261f] - w[rnog261b];

	cdot[sCH3OH] = w[rnog147f] - w[rnog147b];

	cdot[sC] = w[rnog49f] - w[rnog49b] - w[rnog90f] + w[rnog90b]
		- w[rnog122f] + w[rnog122b] - w[rnog123f] + w[rnog123b]
		- w[rnog124f] + w[rnog124b] - w[rnog239f] + w[rnog239b]
		- w[rnog244f] + w[rnog244b] - w[rnog245f] + w[rnog245b];

	cdot[sC2H] = w[rnog123f] - w[rnog123b];

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

	cdot[sHOCN] = w[rnog234f] - w[rnog234b] + w[rnog235] + w[rnog249f]
		- w[rnog249b] + w[rnog251f] - w[rnog251b] + w[rnog252f]
		- w[rnog252b] + w[rnog254f] - w[rnog254b] - w[rnog262f]
		+ w[rnog262b] - w[rnog263f] + w[rnog263b] - w[rnog264f]
		+ w[rnog264b] - w[rnog265f] + w[rnog265b] - w[rnog266f]
		+ w[rnog266b] - w[rnog267f] + w[rnog267b] - w[rnog268f]
		+ w[rnog268b] - w[rnog269f] + w[rnog269b] - w[rnog270f]
		+ w[rnog270f] - w[rnog270b] + w[rnog270b] - w[rnog271f]
		+ w[rnog271b] - w[rnog272f] + w[rnog272b] - w[rnog273]
		+ w[rnog273] + w[rnog274f] - w[rnog274b];

	cdot[sH2CN] = w[rnog237f] - w[rnog237b] - w[rnog238f] + w[rnog238b]
		+ w[rnog256f] - w[rnog256b] + w[rnog275f] - w[rnog275b];

	cdot[sHCNN] = w[rnog241f] - w[rnog241b] - w[rnog257f] + w[rnog257b]
		- w[rnog258f] + w[rnog258b] - w[rnog259f] + w[rnog259b]
		- w[rnog260f] + w[rnog260b] - w[rnog261f] + w[rnog261b];

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
	W[sO2] =  3.20000000e+01;
	W[sH] =  1.00800000e+00;
	W[sOH] =  1.70080000e+01;
	W[sO] =  1.60000000e+01;
	W[sH2] =  2.01600000e+00;
	W[sH2O] =  1.80160000e+01;
	W[sHO2] =  3.30080000e+01;
	W[sH2O2] =  3.40160000e+01;
	W[sCO] =  2.80100000e+01;
	W[sCO2] =  4.40100000e+01;
	W[sCH] =  1.30180000e+01;
	W[sHCO] =  2.90180000e+01;
	W[sCH2OH] =  3.10340000e+01;
	W[s3XCH2] =  1.40260000e+01;
	W[sC2H2] =  2.60360000e+01;
	W[sCH3] =  1.50340000e+01;
	W[sC2H4] =  2.80520000e+01;
	W[s1XCH2] =  1.40260000e+01;
	W[sCH2O] =  3.00260000e+01;
	W[sCH4] =  1.60420000e+01;
	W[sC2H6] =  3.00680000e+01;
	W[sHCCO] =  4.10280000e+01;
	W[sC2H3] =  2.70440000e+01;
	W[sC2H5] =  2.90600000e+01;
	W[sC3H4] =  4.00620000e+01;
	W[sC3H5] =  4.10700000e+01;
	W[sC3H6] =  4.20780000e+01;
	W[sNXC3H7] =  4.30860000e+01;
	W[s1XC4H8] =  5.61040000e+01;
	W[sPXC4H9] =  5.71120000e+01;
	W[s1XC5H11] =  7.11380000e+01;
	W[sC6H11] =  8.31480000e+01;
	W[s1XC6H12] =  8.41560000e+01;
	W[s1XC7H15] =  9.91900000e+01;
	W[s2XC7H15] =  9.91900000e+01;
	W[sNXC7H16] =  1.00198000e+02;
	W[s1C7H15O2XC7H15O2] =  1.31190000e+02;
	W[s1HEOOHX2XC7H15O2] =  1.31190000e+02;
	W[s2OO7OOH1XO2C7H15O2] =  1.63190000e+02;
	W[s1X2HEOOHXO2C7H14OOH] =  1.63190000e+02;
	W[sOC7OOHXOC7H13OOH] =  1.46182000e+02;
	W[s1X2HEPOOXOC7H13O] =  1.29174000e+02;
	W[sCH2] =  1.40260000e+01;
	W[sCH3OH] =  3.20420000e+01;
	W[sC] =  1.20100000e+01;
	W[sC2H] =  2.50280000e+01;
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
	W[sH2CN] =  2.80360000e+01;
	W[sHCNN] =  4.10380000e+01;
	W[sNH3] =  1.70340000e+01;
}

void GetSpeciesNames( char **names )
{
	strcpy( names[sN2], "N2                  " );
	strcpy( names[sO2], "O2                  " );
	strcpy( names[sH], "H                   " );
	strcpy( names[sOH], "OH                  " );
	strcpy( names[sO], "O                   " );
	strcpy( names[sH2], "H2                  " );
	strcpy( names[sH2O], "H2O                 " );
	strcpy( names[sHO2], "HO2                 " );
	strcpy( names[sH2O2], "H2O2                " );
	strcpy( names[sCO], "CO                  " );
	strcpy( names[sCO2], "CO2                 " );
	strcpy( names[sCH], "CH                  " );
	strcpy( names[sHCO], "HCO                 " );
	strcpy( names[sCH2OH], "CH2OH               " );
	strcpy( names[s3XCH2], "3-CH2               " );
	strcpy( names[sC2H2], "C2H2                " );
	strcpy( names[sCH3], "CH3                 " );
	strcpy( names[sC2H4], "C2H4                " );
	strcpy( names[s1XCH2], "1-CH2               " );
	strcpy( names[sCH2O], "CH2O                " );
	strcpy( names[sCH4], "CH4                 " );
	strcpy( names[sC2H6], "C2H6                " );
	strcpy( names[sHCCO], "HCCO                " );
	strcpy( names[sC2H3], "C2H3                " );
	strcpy( names[sC2H5], "C2H5                " );
	strcpy( names[sC3H4], "C3H4                " );
	strcpy( names[sC3H5], "C3H5                " );
	strcpy( names[sC3H6], "C3H6                " );
	strcpy( names[sNXC3H7], "N-C3H7              " );
	strcpy( names[s1XC4H8], "1-C4H8              " );
	strcpy( names[sPXC4H9], "P-C4H9              " );
	strcpy( names[s1XC5H11], "1-C5H11             " );
	strcpy( names[sC6H11], "C6H11               " );
	strcpy( names[s1XC6H12], "1-C6H12             " );
	strcpy( names[s1XC7H15], "1-C7H15             " );
	strcpy( names[s2XC7H15], "2-C7H15             " );
	strcpy( names[sNXC7H16], "N-C7H16             " );
/* 	strcpy( names[s1C7H15O2XC7H15O2], "1C7H15O2-C7H15O2    " ); */
/* 	strcpy( names[s1HEOOHX2XC7H15O2], "1HEOOH-2-C7H15O2    " ); */
/* 	strcpy( names[s2OO7OOH1XO2C7H15O2], "2OO7OOH1-O2C7H15O2  " ); */
/* 	strcpy( names[s1X2HEOOHXO2C7H14OOH], "1-2HEOOH-O2C7H14OOH " ); */
/* 	strcpy( names[sOC7OOHXOC7H13OOH], "OC7OOH-OC7H13OOH    " ); */
/* 	strcpy( names[s1X2HEPOOXOC7H13O], "1-2HEPOO-OC7H13O    " ); */
	strcpy( names[s1C7H15O2XC7H15O2], "1C7H15O2    " );
	strcpy( names[s1HEOOHX2XC7H15O2], "1HEOOH-2    " );
	strcpy( names[s2OO7OOH1XO2C7H15O2], "2OO7OOH1  " );
	strcpy( names[s1X2HEOOHXO2C7H14OOH], "1-2HEOOH " );
	strcpy( names[sOC7OOHXOC7H13OOH], "OC7OOH    " );
	strcpy( names[s1X2HEPOOXOC7H13O], "1-2HEPOO    " );
	strcpy( names[sCH2], "CH2                 " );
	strcpy( names[sCH3OH], "CH3OH               " );
	strcpy( names[sC], "C                   " );
	strcpy( names[sC2H], "C2H                 " );
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
	strcpy( names[sH2CN], "H2CN                " );
	strcpy( names[sHCNN], "HCNN                " );
	strcpy( names[sNH3], "NH3                 " );
}



void ComputeThermoData( double *h, double *cp, double T )
{
/*
	This function computes enthalpy 'h' and heat capacity 'cp' as
	function of temperature 'T' for all non steady state species
	in units [J/kg] and [J/kg K], respectively.
	The parameter h and cp should provide workspace of length 62 */


	int i;
	if ( T > 1000.0 ) {
		h[sN2] =  2.96728765e+02 * (
			T * (  2.85328990e+00 + T * (  8.01106400e-04
			+ T * ( -2.09789643e-07 + T * (  2.86025550e-11
			+ T * -1.56114930e-15 ) ) ) ) -8.90080930e+02 );
		cp[sN2] =  2.96728765e+02 * (
			 2.85328990e+00 + T * (  1.60221280e-03 
			+ T * ( -6.29368930e-07 + T * (  1.14410220e-10
			+ T * -7.80574650e-15 ) ) ) );
		h[sO2] =  2.59823125e+02 * (
			T * (  3.61221390e+00 + T * (  3.74265830e-04
			+ T * ( -6.60688233e-08 + T * (  8.43725200e-12
			+ T * -4.78147480e-16 ) ) ) ) -1.19781510e+03 );
		cp[sO2] =  2.59823125e+02 * (
			 3.61221390e+00 + T * (  7.48531660e-04 
			+ T * ( -1.98206470e-07 + T * (  3.37490080e-11
			+ T * -2.39073740e-15 ) ) ) );
		h[sH] =  8.24835317e+03 * (
			T * (  2.50000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  2.54716270e+04 );
		cp[sH] =  8.24835317e+03 * (
			 2.50000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sOH] =  4.88848777e+02 * (
			T * (  2.88978140e+00 + T * (  5.00293950e-04
			+ T * ( -7.34960233e-08 + T * (  5.04782200e-12
			+ T * -7.88196620e-17 ) ) ) ) +  3.88570420e+03 );
		cp[sOH] =  4.88848777e+02 * (
			 2.88978140e+00 + T * (  1.00058790e-03 
			+ T * ( -2.20488070e-07 + T * (  2.01912880e-11
			+ T * -3.94098310e-16 ) ) ) );
		h[sO] =  5.19646250e+02 * (
			T * (  2.53429610e+00 + T * ( -6.23908500e-06
			+ T * ( -4.18757467e-09 + T * (  1.72574655e-12
			+ T * -1.27594190e-16 ) ) ) ) +  2.92311080e+04 );
		cp[sO] =  5.19646250e+02 * (
			 2.53429610e+00 + T * ( -1.24781700e-05 
			+ T * ( -1.25627240e-08 + T * (  6.90298620e-12
			+ T * -6.37970950e-16 ) ) ) );
		h[sH2] =  4.12417659e+03 * (
			T * (  3.06670950e+00 + T * (  2.87368775e-04
			+ T * (  4.64610633e-09 + T * ( -6.37087950e-12
			+ T *  5.81971480e-16 ) ) ) ) -8.65474120e+02 );
		cp[sH2] =  4.12417659e+03 * (
			 3.06670950e+00 + T * (  5.74737550e-04 
			+ T * (  1.39383190e-08 + T * ( -2.54835180e-11
			+ T *  2.90985740e-15 ) ) ) );
		h[sH2O] =  4.61497558e+02 * (
			T * (  2.61104720e+00 + T * (  1.57815650e-03
			+ T * ( -3.09951460e-07 + T * (  3.33288450e-11
			+ T * -1.49378702e-15 ) ) ) ) -2.98681670e+04 );
		cp[sH2O] =  4.61497558e+02 * (
			 2.61104720e+00 + T * (  3.15631300e-03 
			+ T * ( -9.29854380e-07 + T * (  1.33315380e-10
			+ T * -7.46893510e-15 ) ) ) );
		h[sHO2] =  2.51888633e+02 * (
			T * (  3.78662800e+00 + T * (  1.39427020e-03
			+ T * ( -3.38956933e-07 + T * (  4.29598650e-11
			+ T * -2.20437040e-15 ) ) ) ) +  1.18885000e+03 );
		cp[sHO2] =  2.51888633e+02 * (
			 3.78662800e+00 + T * (  2.78854040e-03 
			+ T * ( -1.01687080e-06 + T * (  1.71839460e-10
			+ T * -1.10218520e-14 ) ) ) );
		h[sH2O2] =  2.44424389e+02 * (
			T * (  4.57316670e+00 + T * (  2.16806815e-03
			+ T * ( -4.91562933e-07 + T * (  5.87225925e-11
			+ T * -2.86330720e-15 ) ) ) ) -1.80069610e+04 );
		cp[sH2O2] =  2.44424389e+02 * (
			 4.57316670e+00 + T * (  4.33613630e-03 
			+ T * ( -1.47468880e-06 + T * (  2.34890370e-10
			+ T * -1.43165360e-14 ) ) ) );
		h[sCO] =  2.96834702e+02 * (
			T * (  3.02507810e+00 + T * (  7.21344250e-04
			+ T * ( -1.87694260e-07 + T * (  2.54645325e-11
			+ T * -1.38219032e-15 ) ) ) ) -1.42683500e+04 );
		cp[sCO] =  2.96834702e+02 * (
			 3.02507810e+00 + T * (  1.44268850e-03 
			+ T * ( -5.63082780e-07 + T * (  1.01858130e-10
			+ T * -6.91095160e-15 ) ) ) );
		h[sCO2] =  1.88919337e+02 * (
			T * (  4.45362280e+00 + T * (  1.57008435e-03
			+ T * ( -4.26136833e-07 + T * (  5.98499175e-11
			+ T * -3.33806640e-15 ) ) ) ) -4.89669610e+04 );
		cp[sCO2] =  1.88919337e+02 * (
			 4.45362280e+00 + T * (  3.14016870e-03 
			+ T * ( -1.27841050e-06 + T * (  2.39399670e-10
			+ T * -1.66903320e-14 ) ) ) );
		h[sCH] =  6.38680289e+02 * (
			T * (  2.19622330e+00 + T * (  1.17019050e-03
			+ T * ( -2.35273377e-07 + T * (  2.25189555e-11
			+ T * -7.71008020e-16 ) ) ) ) +  7.08672340e+04 );
		cp[sCH] =  6.38680289e+02 * (
			 2.19622330e+00 + T * (  2.34038100e-03 
			+ T * ( -7.05820130e-07 + T * (  9.00758220e-11
			+ T * -3.85504010e-15 ) ) ) );
		h[sHCO] =  2.86523537e+02 * (
			T * (  3.47383480e+00 + T * (  1.71851135e-03
			+ T * ( -4.54422133e-07 + T * (  6.23216125e-11
			+ T * -3.40886620e-15 ) ) ) ) +  3.95940050e+03 );
		cp[sHCO] =  2.86523537e+02 * (
			 3.47383480e+00 + T * (  3.43702270e-03 
			+ T * ( -1.36326640e-06 + T * (  2.49286450e-10
			+ T * -1.70443310e-14 ) ) ) );
		h[sCH2OH] =  2.67910679e+02 * (
			T * (  4.72350410e+00 + T * (  3.05102205e-03
			+ T * ( -6.37736467e-07 + T * (  6.90185675e-11
			+ T * -2.90967340e-15 ) ) ) ) -3.93291650e+03 );
		cp[sCH2OH] =  2.67910679e+02 * (
			 4.72350410e+00 + T * (  6.10204410e-03 
			+ T * ( -1.91320940e-06 + T * (  2.76074270e-10
			+ T * -1.45483670e-14 ) ) ) );
		h[s3XCH2] =  5.92780550e+02 * (
			T * (  3.63640790e+00 + T * (  9.66528300e-04
			+ T * ( -5.62338767e-08 + T * ( -2.52474850e-11
			+ T *  3.61651160e-15 ) ) ) ) +  4.53413400e+04 );
		cp[s3XCH2] =  5.92780550e+02 * (
			 3.63640790e+00 + T * (  1.93305660e-03 
			+ T * ( -1.68701630e-07 + T * ( -1.00989940e-10
			+ T *  1.80825580e-14 ) ) ) );
		h[sC2H2] =  3.19340144e+02 * (
			T * (  4.43677040e+00 + T * (  2.68801955e-03
			+ T * ( -6.37605567e-07 + T * (  8.21594725e-11
			+ T * -4.31341900e-15 ) ) ) ) +  2.56676640e+04 );
		cp[sC2H2] =  3.19340144e+02 * (
			 4.43677040e+00 + T * (  5.37603910e-03 
			+ T * ( -1.91281670e-06 + T * (  3.28637890e-10
			+ T * -2.15670950e-14 ) ) ) );
		h[sCH3] =  5.53035786e+02 * (
			T * (  2.84405160e+00 + T * (  3.06898705e-03
			+ T * ( -7.43448400e-07 + T * (  9.46290200e-11
			+ T * -4.90431800e-15 ) ) ) ) +  1.64378090e+04 );
		cp[sCH3] =  5.53035786e+02 * (
			 2.84405160e+00 + T * (  6.13797410e-03 
			+ T * ( -2.23034520e-06 + T * (  3.78516080e-10
			+ T * -2.45215900e-14 ) ) ) );
		h[sC2H4] =  2.96390275e+02 * (
			T * (  3.52841880e+00 + T * (  5.74259250e-03
			+ T * ( -1.47279510e-06 + T * (  1.96115013e-10
			+ T * -1.05336970e-14 ) ) ) ) +  4.42828860e+03 );
		cp[sC2H4] =  2.96390275e+02 * (
			 3.52841880e+00 + T * (  1.14851850e-02 
			+ T * ( -4.41838530e-06 + T * (  7.84460050e-10
			+ T * -5.26684850e-14 ) ) ) );
		h[s1XCH2] =  5.92780550e+02 * (
			T * (  3.55288860e+00 + T * (  1.03339415e-03
			+ T * ( -6.38038667e-08 + T * ( -2.76168350e-11
			+ T *  4.04269920e-15 ) ) ) ) +  4.98497540e+04 );
		cp[s1XCH2] =  5.92780550e+02 * (
			 3.55288860e+00 + T * (  2.06678830e-03 
			+ T * ( -1.91411600e-07 + T * ( -1.10467340e-10
			+ T *  2.02134960e-14 ) ) ) );
		h[sCH2O] =  2.76904683e+02 * (
			T * (  2.99560620e+00 + T * (  3.34066060e-03
			+ T * ( -8.76318233e-07 + T * (  1.18428823e-10
			+ T * -6.42503500e-15 ) ) ) ) -1.53203690e+04 );
		cp[sCH2O] =  2.76904683e+02 * (
			 2.99560620e+00 + T * (  6.68132120e-03 
			+ T * ( -2.62895470e-06 + T * (  4.73715290e-10
			+ T * -3.21251750e-14 ) ) ) );
		h[sCH4] =  5.18285750e+02 * (
			T * (  1.68347880e+00 + T * (  5.11861800e-03
			+ T * ( -1.29170953e-06 + T * (  1.69639622e-10
			+ T * -9.00684620e-15 ) ) ) ) -1.00807870e+04 );
		cp[sCH4] =  5.18285750e+02 * (
			 1.68347880e+00 + T * (  1.02372360e-02 
			+ T * ( -3.87512860e-06 + T * (  6.78558490e-10
			+ T * -4.50342310e-14 ) ) ) );
		h[sC2H6] =  2.76517893e+02 * (
			T * (  4.82593820e+00 + T * (  6.92021450e-03
			+ T * ( -1.51908627e-06 + T * (  1.68124180e-10
			+ T * -7.19632280e-15 ) ) ) ) -1.27177930e+04 );
		cp[sC2H6] =  2.76517893e+02 * (
			 4.82593820e+00 + T * (  1.38404290e-02 
			+ T * ( -4.55725880e-06 + T * (  6.72496720e-10
			+ T * -3.59816140e-14 ) ) ) );
		h[sHCCO] =  2.02650385e+02 * (
			T * (  6.75807290e+00 + T * (  1.00020015e-03
			+ T * ( -6.75869100e-08 + T * ( -2.60282950e-11
			+ T *  3.93032940e-15 ) ) ) ) +  1.90151330e+04 );
		cp[sHCCO] =  2.02650385e+02 * (
			 6.75807290e+00 + T * (  2.00040030e-03 
			+ T * ( -2.02760730e-07 + T * ( -1.04113180e-10
			+ T *  1.96516470e-14 ) ) ) );
		h[sC2H3] =  3.07437509e+02 * (
			T * (  5.93346790e+00 + T * (  2.00887280e-03
			+ T * ( -1.32224650e-07 + T * ( -3.60316650e-11
			+ T *  4.75728700e-15 ) ) ) ) +  3.18543460e+04 );
		cp[sC2H3] =  3.07437509e+02 * (
			 5.93346790e+00 + T * (  4.01774560e-03 
			+ T * ( -3.96673950e-07 + T * ( -1.44126660e-10
			+ T *  2.37864350e-14 ) ) ) );
		h[sC2H5] =  2.86109429e+02 * (
			T * (  7.19048020e+00 + T * (  3.24203865e-03
			+ T * ( -2.14268823e-07 + T * ( -5.86969850e-11
			+ T *  7.76175460e-15 ) ) ) ) +  1.06745490e+04 );
		cp[sC2H5] =  2.86109429e+02 * (
			 7.19048020e+00 + T * (  6.48407730e-03 
			+ T * ( -6.42806470e-07 + T * ( -2.34787940e-10
			+ T *  3.88087730e-14 ) ) ) );
		h[sC3H4] =  2.07536818e+02 * (
			T * (  6.53446390e+00 + T * (  5.16887200e-03
			+ T * ( -1.13336863e-06 + T * (  1.22728437e-10
			+ T * -5.07274960e-15 ) ) ) ) +  1.94417070e+04 );
		cp[sC3H4] =  2.07536818e+02 * (
			 6.53446390e+00 + T * (  1.03377440e-02 
			+ T * ( -3.40010590e-06 + T * (  4.90913750e-10
			+ T * -2.53637480e-14 ) ) ) );
		h[sC3H5] =  2.02443146e+02 * (
			T * (  7.90919780e+00 + T * (  6.05762750e-03
			+ T * ( -1.37252877e-06 + T * (  1.53916990e-10
			+ T * -6.64714660e-15 ) ) ) ) +  1.23541560e+04 );
		cp[sC3H5] =  2.02443146e+02 * (
			 7.90919780e+00 + T * (  1.21152550e-02 
			+ T * ( -4.11758630e-06 + T * (  6.15667960e-10
			+ T * -3.32357330e-14 ) ) ) );
		h[sC3H6] =  1.97593517e+02 * (
			T * (  6.73225690e+00 + T * (  7.45416800e-03
			+ T * ( -1.64996647e-06 + T * (  1.80300552e-10
			+ T * -7.53240860e-15 ) ) ) ) -9.23570310e+02 );
		cp[sC3H6] =  1.97593517e+02 * (
			 6.73225690e+00 + T * (  1.49083360e-02 
			+ T * ( -4.94989940e-06 + T * (  7.21202210e-10
			+ T * -3.76620430e-14 ) ) ) );
		h[sNXC3H7] =  1.92970803e+02 * (
			T * (  7.97829060e+00 + T * (  7.88056650e-03
			+ T * ( -1.72441437e-06 + T * (  1.86097305e-10
			+ T * -7.64995640e-15 ) ) ) ) +  7.57940230e+03 );
		cp[sNXC3H7] =  1.92970803e+02 * (
			 7.97829060e+00 + T * (  1.57611330e-02 
			+ T * ( -5.17324310e-06 + T * (  7.44389220e-10
			+ T * -3.82497820e-14 ) ) ) );
		h[s1XC4H8] =  1.48195138e+02 * (
			T * (  2.05358410e+00 + T * (  1.71752535e-02
			+ T * ( -5.29439900e-06 + T * (  8.27241550e-10
			+ T * -5.07220900e-14 ) ) ) ) -2.13972310e+03 );
		cp[s1XC4H8] =  1.48195138e+02 * (
			 2.05358410e+00 + T * (  3.43505070e-02 
			+ T * ( -1.58831970e-05 + T * (  3.30896620e-09
			+ T * -2.53610450e-13 ) ) ) );
		h[sPXC4H9] =  1.45579563e+02 * (
			T * (  9.00491050e+00 + T * (  1.20016595e-02
			+ T * ( -2.98506227e-06 + T * (  3.83064850e-10
			+ T * -1.96494080e-14 ) ) ) ) +  3.46301120e+03 );
		cp[sPXC4H9] =  1.45579563e+02 * (
			 9.00491050e+00 + T * (  2.40033190e-02 
			+ T * ( -8.95518680e-06 + T * (  1.53225940e-09
			+ T * -9.82470400e-14 ) ) ) );
		h[s1XC5H11] =  1.16876212e+02 * (
			T * (  1.06551910e+01 + T * (  1.54696330e-02
			+ T * ( -3.93019900e-06 + T * (  5.13945775e-10
			+ T * -2.68096960e-14 ) ) ) ) +  1.91752170e+01 );
		cp[s1XC5H11] =  1.16876212e+02 * (
			 1.06551910e+01 + T * (  3.09392660e-02 
			+ T * ( -1.17905970e-05 + T * (  2.05578310e-09
			+ T * -1.34048480e-13 ) ) ) );
		h[sC6H11] =  9.99944677e+01 * (
			T * (  1.59277710e+01 + T * (  1.18720645e-02
			+ T * ( -2.30289060e-06 + T * (  2.02744433e-10
			+ T * -5.36624520e-15 ) ) ) ) -8.64265620e+03 );
		cp[sC6H11] =  9.99944677e+01 * (
			 1.59277710e+01 + T * (  2.37441290e-02 
			+ T * ( -6.90867180e-06 + T * (  8.10977730e-10
			+ T * -2.68312260e-14 ) ) ) );
		h[s1XC6H12] =  9.87967584e+01 * (
			T * (  1.86636350e+01 + T * (  1.04857255e-02
			+ T * ( -1.03609363e-06 + T * ( -1.71629045e-10
			+ T *  3.20472160e-14 ) ) ) ) -1.35908950e+04 );
		cp[s1XC6H12] =  9.87967584e+01 * (
			 1.86636350e+01 + T * (  2.09714510e-02 
			+ T * ( -3.10828090e-06 + T * ( -6.86516180e-10
			+ T *  1.60236080e-13 ) ) ) );
		h[s1XC7H15] =  8.38223611e+01 * (
			T * (  1.49188600e+01 + T * (  2.16028655e-02
			+ T * ( -5.54310600e-06 + T * (  7.31940400e-10
			+ T * -3.85228520e-14 ) ) ) ) -7.17518750e+03 );
		cp[s1XC7H15] =  8.38223611e+01 * (
			 1.49188600e+01 + T * (  4.32057310e-02 
			+ T * ( -1.66293180e-05 + T * (  2.92776160e-09
			+ T * -1.92614260e-13 ) ) ) );
		h[s2XC7H15] =  8.38223611e+01 * (
			T * (  1.49188600e+01 + T * (  2.16028655e-02
			+ T * ( -5.54310600e-06 + T * (  7.31940400e-10
			+ T * -3.85228520e-14 ) ) ) ) -7.17518750e+03 );
		cp[s2XC7H15] =  8.38223611e+01 * (
			 1.49188600e+01 + T * (  4.32057310e-02 
			+ T * ( -1.66293180e-05 + T * (  2.92776160e-09
			+ T * -1.92614260e-13 ) ) ) );
		h[sNXC7H16] =  8.29791014e+01 * (
			T * (  2.28188930e+01 + T * (  1.62717270e-02
			+ T * ( -3.70668033e-06 + T * (  4.28293575e-10
			+ T * -1.92424202e-14 ) ) ) ) -3.36787380e+04 );
		cp[sNXC7H16] =  8.29791014e+01 * (
			 2.28188930e+01 + T * (  3.25434540e-02 
			+ T * ( -1.11200410e-05 + T * (  1.71317430e-09
			+ T * -9.62121010e-14 ) ) ) );
		h[s1C7H15O2XC7H15O2] =  6.33763244e+01 * (
			T * (  2.85328990e+00 + T * (  8.01106400e-04
			+ T * ( -2.09789643e-07 + T * (  2.86025550e-11
			+ T * -1.56114930e-15 ) ) ) ) -8.90080930e+02 );
		cp[s1C7H15O2XC7H15O2] =  6.33763244e+01 * (
			 2.85328990e+00 + T * (  1.60221280e-03 
			+ T * ( -6.29368930e-07 + T * (  1.14410220e-10
			+ T * -7.80574650e-15 ) ) ) );
		h[s1HEOOHX2XC7H15O2] =  6.33763244e+01 * (
			T * (  2.85328990e+00 + T * (  8.01106400e-04
			+ T * ( -2.09789643e-07 + T * (  2.86025550e-11
			+ T * -1.56114930e-15 ) ) ) ) -8.90080930e+02 );
		cp[s1HEOOHX2XC7H15O2] =  6.33763244e+01 * (
			 2.85328990e+00 + T * (  1.60221280e-03 
			+ T * ( -6.29368930e-07 + T * (  1.14410220e-10
			+ T * -7.80574650e-15 ) ) ) );
		h[s2OO7OOH1XO2C7H15O2] =  5.09488326e+01 * (
			T * (  2.85328990e+00 + T * (  8.01106400e-04
			+ T * ( -2.09789643e-07 + T * (  2.86025550e-11
			+ T * -1.56114930e-15 ) ) ) ) -8.90080930e+02 );
		cp[s2OO7OOH1XO2C7H15O2] =  5.09488326e+01 * (
			 2.85328990e+00 + T * (  1.60221280e-03 
			+ T * ( -6.29368930e-07 + T * (  1.14410220e-10
			+ T * -7.80574650e-15 ) ) ) );
		h[s1X2HEOOHXO2C7H14OOH] =  5.09488326e+01 * (
			T * (  2.85328990e+00 + T * (  8.01106400e-04
			+ T * ( -2.09789643e-07 + T * (  2.86025550e-11
			+ T * -1.56114930e-15 ) ) ) ) -8.90080930e+02 );
		cp[s1X2HEOOHXO2C7H14OOH] =  5.09488326e+01 * (
			 2.85328990e+00 + T * (  1.60221280e-03 
			+ T * ( -6.29368930e-07 + T * (  1.14410220e-10
			+ T * -7.80574650e-15 ) ) ) );
		h[sOC7OOHXOC7H13OOH] =  5.68766332e+01 * (
			T * (  2.85328990e+00 + T * (  8.01106400e-04
			+ T * ( -2.09789643e-07 + T * (  2.86025550e-11
			+ T * -1.56114930e-15 ) ) ) ) -8.90080930e+02 );
		cp[sOC7OOHXOC7H13OOH] =  5.68766332e+01 * (
			 2.85328990e+00 + T * (  1.60221280e-03 
			+ T * ( -6.29368930e-07 + T * (  1.14410220e-10
			+ T * -7.80574650e-15 ) ) ) );
		h[s1X2HEPOOXOC7H13O] =  6.43654296e+01 * (
			T * (  2.85328990e+00 + T * (  8.01106400e-04
			+ T * ( -2.09789643e-07 + T * (  2.86025550e-11
			+ T * -1.56114930e-15 ) ) ) ) -8.90080930e+02 );
		cp[s1X2HEPOOXOC7H13O] =  6.43654296e+01 * (
			 2.85328990e+00 + T * (  1.60221280e-03 
			+ T * ( -6.29368930e-07 + T * (  1.14410220e-10
			+ T * -7.80574650e-15 ) ) ) );
		h[sCH2] =  5.92780550e+02 * (
			T * (  3.63640800e+00 + T * (  9.66528500e-04
			+ T * ( -5.62338667e-08 + T * ( -2.52474750e-11
			+ T *  3.61651200e-15 ) ) ) ) +  4.53413400e+04 );
		cp[sCH2] =  5.92780550e+02 * (
			 3.63640800e+00 + T * (  1.93305700e-03 
			+ T * ( -1.68701600e-07 + T * ( -1.00989900e-10
			+ T *  1.80825600e-14 ) ) ) );
		h[sCH3OH] =  2.59482554e+02 * (
			T * (  4.02906100e+00 + T * (  4.68829650e-03
			+ T * ( -1.01675133e-06 + T * (  1.08969825e-10
			+ T * -4.44944600e-15 ) ) ) ) -2.61579100e+04 );
		cp[sCH3OH] =  2.59482554e+02 * (
			 4.02906100e+00 + T * (  9.37659300e-03 
			+ T * ( -3.05025400e-06 + T * (  4.35879300e-10
			+ T * -2.22472300e-14 ) ) ) );
		h[sC] =  6.92284763e+02 * (
			T * (  2.49266888e+00 + T * (  2.39944642e-05
			+ T * ( -2.41445007e-08 + T * (  9.35727573e-12
			+ T * -9.74555786e-16 ) ) ) ) +  8.54512953e+04 );
		cp[sC] =  6.92284763e+02 * (
			 2.49266888e+00 + T * (  4.79889284e-05 
			+ T * ( -7.24335020e-08 + T * (  3.74291029e-11
			+ T * -4.87277893e-15 ) ) ) );
		h[sC2H] =  3.32201534e+02 * (
			T * (  3.98636700e+00 + T * (  1.57156150e-03
			+ T * ( -4.22414333e-07 + T * (  7.31090750e-11
			+ T * -5.43264000e-15 ) ) ) ) +  6.65588400e+04 );
		cp[sC2H] =  3.32201534e+02 * (
			 3.98636700e+00 + T * (  3.14312300e-03 
			+ T * ( -1.26724300e-06 + T * (  2.92436300e-10
			+ T * -2.71632000e-14 ) ) ) );
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
			T * (  3.70441770e+00 + T * ( -7.10937650e-04
			+ T * (  9.55679733e-07 + T * ( -3.00722125e-10
			+ T * -2.79093540e-15 ) ) ) ) -1.06407950e+03 );
		cp[sN2] =  2.96728765e+02 * (
			 3.70441770e+00 + T * ( -1.42187530e-03 
			+ T * (  2.86703920e-06 + T * ( -1.20288850e-09
			+ T * -1.39546770e-14 ) ) ) );
		h[sO2] =  2.59823125e+02 * (
			T * (  3.78371350e+00 + T * ( -1.51168170e-03
			+ T * (  3.31642503e-06 + T * ( -2.45472752e-09
			+ T *  6.60636500e-13 ) ) ) ) -1.06381070e+03 );
		cp[sO2] =  2.59823125e+02 * (
			 3.78371350e+00 + T * ( -3.02336340e-03 
			+ T * (  9.94927510e-06 + T * ( -9.81891010e-09
			+ T *  3.30318250e-12 ) ) ) );
		h[sH] =  8.24835317e+03 * (
			T * (  2.50000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  2.54716270e+04 );
		cp[sH] =  8.24835317e+03 * (
			 2.50000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sOH] =  4.88848777e+02 * (
			T * (  3.87373000e+00 + T * ( -6.69688600e-04
			+ T * (  5.44945033e-07 + T * ( -1.30334098e-10
			+ T *  8.36539480e-15 ) ) ) ) +  3.58023480e+03 );
		cp[sOH] =  4.88848777e+02 * (
			 3.87373000e+00 + T * ( -1.33937720e-03 
			+ T * (  1.63483510e-06 + T * ( -5.21336390e-10
			+ T *  4.18269740e-14 ) ) ) );
		h[sO] =  5.19646250e+02 * (
			T * (  3.03094010e+00 + T * ( -1.12629265e-03
			+ T * (  1.32748467e-06 + T * ( -8.15123025e-10
			+ T *  2.03040700e-13 ) ) ) ) +  2.91365260e+04 );
		cp[sO] =  5.19646250e+02 * (
			 3.03094010e+00 + T * ( -2.25258530e-03 
			+ T * (  3.98245400e-06 + T * ( -3.26049210e-09
			+ T *  1.01520350e-12 ) ) ) );
		h[sH2] =  4.12417659e+03 * (
			T * (  3.35535140e+00 + T * (  2.50680720e-04
			+ T * ( -7.66896933e-08 + T * ( -1.19763310e-10
			+ T *  9.70451700e-14 ) ) ) ) -1.01916260e+03 );
		cp[sH2] =  4.12417659e+03 * (
			 3.35535140e+00 + T * (  5.01361440e-04 
			+ T * ( -2.30069080e-07 + T * ( -4.79053240e-10
			+ T *  4.85225850e-13 ) ) ) );
		h[sH2O] =  4.61497558e+02 * (
			T * (  4.16772340e+00 + T * ( -9.05748500e-04
			+ T * (  1.98237627e-06 + T * ( -1.21730053e-09
			+ T *  3.05839820e-13 ) ) ) ) -3.02899690e+04 );
		cp[sH2O] =  4.61497558e+02 * (
			 4.16772340e+00 + T * ( -1.81149700e-03 
			+ T * (  5.94712880e-06 + T * ( -4.86920210e-09
			+ T *  1.52919910e-12 ) ) ) );
		h[sHO2] =  2.51888633e+02 * (
			T * (  3.50948500e+00 + T * (  5.74983500e-04
			+ T * (  1.95947530e-06 + T * ( -1.94488798e-09
			+ T *  5.92157660e-13 ) ) ) ) +  1.38033310e+03 );
		cp[sHO2] =  2.51888633e+02 * (
			 3.50948500e+00 + T * (  1.14996700e-03 
			+ T * (  5.87842590e-06 + T * ( -7.77955190e-09
			+ T *  2.96078830e-12 ) ) ) );
		h[sH2O2] =  2.44424389e+02 * (
			T * (  3.38875360e+00 + T * (  3.28461300e-03
			+ T * ( -4.95004200e-08 + T * ( -1.15645138e-09
			+ T *  4.94302940e-13 ) ) ) ) -1.76631470e+04 );
		cp[sH2O2] =  2.44424389e+02 * (
			 3.38875360e+00 + T * (  6.56922600e-03 
			+ T * ( -1.48501260e-07 + T * ( -4.62580550e-09
			+ T *  2.47151470e-12 ) ) ) );
		h[sCO] =  2.96834702e+02 * (
			T * (  3.26245170e+00 + T * (  7.55970450e-04
			+ T * ( -1.29391840e-06 + T * (  1.39548605e-09
			+ T * -4.94990240e-13 ) ) ) ) -1.43105390e+04 );
		cp[sCO] =  2.96834702e+02 * (
			 3.26245170e+00 + T * (  1.51194090e-03 
			+ T * ( -3.88175520e-06 + T * (  5.58194420e-09
			+ T * -2.47495120e-12 ) ) ) );
		h[sCO2] =  1.88919337e+02 * (
			T * (  2.27572470e+00 + T * (  4.96103615e-03
			+ T * ( -3.46970433e-06 + T * (  1.71667170e-09
			+ T * -4.23456020e-13 ) ) ) ) -4.83731410e+04 );
		cp[sCO2] =  1.88919337e+02 * (
			 2.27572470e+00 + T * (  9.92207230e-03 
			+ T * ( -1.04091130e-05 + T * (  6.86668680e-09
			+ T * -2.11728010e-12 ) ) ) );
		h[sCH] =  6.38680289e+02 * (
			T * (  3.20020250e+00 + T * (  1.03643780e-03
			+ T * ( -1.71147713e-06 + T * (  1.43347258e-09
			+ T * -3.91106640e-13 ) ) ) ) +  7.04525940e+04 );
		cp[sCH] =  6.38680289e+02 * (
			 3.20020250e+00 + T * (  2.07287560e-03 
			+ T * ( -5.13443140e-06 + T * (  5.73389030e-09
			+ T * -1.95553320e-12 ) ) ) );
		h[sHCO] =  2.86523537e+02 * (
			T * (  3.88401920e+00 + T * ( -4.14872240e-04
			+ T * (  2.59669363e-06 + T * ( -1.76542405e-09
			+ T *  3.99434600e-13 ) ) ) ) +  4.05638600e+03 );
		cp[sHCO] =  2.86523537e+02 * (
			 3.88401920e+00 + T * ( -8.29744480e-04 
			+ T * (  7.79008090e-06 + T * ( -7.06169620e-09
			+ T *  1.99717300e-12 ) ) ) );
		h[sCH2OH] =  2.67910679e+02 * (
			T * (  3.33684060e+00 + T * (  3.29406185e-03
			+ T * (  9.99310933e-07 + T * ( -1.46799285e-09
			+ T *  4.24591440e-13 ) ) ) ) -3.31682670e+03 );
		cp[sCH2OH] =  2.67910679e+02 * (
			 3.33684060e+00 + T * (  6.58812370e-03 
			+ T * (  2.99793280e-06 + T * ( -5.87197140e-09
			+ T *  2.12295720e-12 ) ) ) );
		h[s3XCH2] =  5.92780550e+02 * (
			T * (  3.76223710e+00 + T * (  5.79909550e-04
			+ T * (  8.29861800e-08 + T * (  2.20020890e-10
			+ T * -1.46648708e-13 ) ) ) ) +  4.53679060e+04 );
		cp[s3XCH2] =  5.92780550e+02 * (
			 3.76223710e+00 + T * (  1.15981910e-03 
			+ T * (  2.48958540e-07 + T * (  8.80083560e-10
			+ T * -7.33243540e-13 ) ) ) );
		h[sC2H2] =  3.19340144e+02 * (
			T * (  2.01356220e+00 + T * (  7.59522300e-03
			+ T * ( -5.38772967e-06 + T * (  2.26974795e-09
			+ T * -3.82549200e-13 ) ) ) ) +  2.61244430e+04 );
		cp[sC2H2] =  3.19340144e+02 * (
			 2.01356220e+00 + T * (  1.51904460e-02 
			+ T * ( -1.61631890e-05 + T * (  9.07899180e-09
			+ T * -1.91274600e-12 ) ) ) );
		h[sCH3] =  5.53035786e+02 * (
			T * (  2.43044280e+00 + T * (  5.56204950e-03
			+ T * ( -5.60073433e-06 + T * (  4.05457175e-09
			+ T * -1.17299052e-12 ) ) ) ) +  1.64237810e+04 );
		cp[sCH3] =  5.53035786e+02 * (
			 2.43044280e+00 + T * (  1.11240990e-02 
			+ T * ( -1.68022030e-05 + T * (  1.62182870e-08
			+ T * -5.86495260e-12 ) ) ) );
		h[sC2H4] =  2.96390275e+02 * (
			T * ( -8.61487980e-01 + T * (  1.39808145e-02
			+ T * ( -1.12955907e-05 + T * (  6.96288050e-09
			+ T * -1.94757578e-12 ) ) ) ) +  5.57304590e+03 );
		cp[sC2H4] =  2.96390275e+02 * (
			-8.61487980e-01 + T * (  2.79616290e-02 
			+ T * ( -3.38867720e-05 + T * (  2.78515220e-08
			+ T * -9.73787890e-12 ) ) ) );
		h[s1XCH2] =  5.92780550e+02 * (
			T * (  3.97126510e+00 + T * ( -8.49544400e-05
			+ T * (  3.41789633e-07 + T * (  6.23137700e-10
			+ T * -3.96253260e-13 ) ) ) ) +  4.98936760e+04 );
		cp[s1XCH2] =  5.92780550e+02 * (
			 3.97126510e+00 + T * ( -1.69908880e-04 
			+ T * (  1.02536890e-06 + T * (  2.49255080e-09
			+ T * -1.98126630e-12 ) ) ) );
		h[sCH2O] =  2.76904683e+02 * (
			T * (  1.65273120e+00 + T * (  6.31571950e-03
			+ T * ( -6.29389500e-06 + T * (  5.12507850e-09
			+ T * -1.68264742e-12 ) ) ) ) -1.48654040e+04 );
		cp[sCH2O] =  2.76904683e+02 * (
			 1.65273120e+00 + T * (  1.26314390e-02 
			+ T * ( -1.88816850e-05 + T * (  2.05003140e-08
			+ T * -8.41323710e-12 ) ) ) );
		h[sCH4] =  5.18285750e+02 * (
			T * (  7.78741480e-01 + T * (  8.73834200e-03
			+ T * ( -9.27803000e-06 + T * (  7.62427000e-09
			+ T * -2.44786140e-12 ) ) ) ) -9.82522850e+03 );
		cp[sCH4] =  5.18285750e+02 * (
			 7.78741480e-01 + T * (  1.74766840e-02 
			+ T * ( -2.78340900e-05 + T * (  3.04970800e-08
			+ T * -1.22393070e-11 ) ) ) );
		h[sC2H6] =  2.76517893e+02 * (
			T * (  1.46253870e+00 + T * (  7.74733350e-03
			+ T * (  1.92683577e-06 + T * ( -3.14457975e-09
			+ T *  9.17253420e-13 ) ) ) ) -1.12391760e+04 );
		cp[sC2H6] =  2.76517893e+02 * (
			 1.46253870e+00 + T * (  1.54946670e-02 
			+ T * (  5.78050730e-06 + T * ( -1.25783190e-08
			+ T *  4.58626710e-12 ) ) ) );
		h[sHCCO] =  2.02650385e+02 * (
			T * (  5.04796510e+00 + T * (  2.22673895e-03
			+ T * (  7.56094267e-08 + T * ( -3.70523650e-10
			+ T *  4.50148300e-14 ) ) ) ) +  1.96589180e+04 );
		cp[sHCCO] =  2.02650385e+02 * (
			 5.04796510e+00 + T * (  4.45347790e-03 
			+ T * (  2.26828280e-07 + T * ( -1.48209460e-09
			+ T *  2.25074150e-13 ) ) ) );
		h[sC2H3] =  3.07437509e+02 * (
			T * (  2.45927640e+00 + T * (  3.68573820e-03
			+ T * (  7.03290967e-07 + T * ( -3.30410525e-10
			+ T * -2.36956760e-13 ) ) ) ) +  3.33522500e+04 );
		cp[sC2H3] =  3.07437509e+02 * (
			 2.45927640e+00 + T * (  7.37147640e-03 
			+ T * (  2.10987290e-06 + T * ( -1.32164210e-09
			+ T * -1.18478380e-12 ) ) ) );
		h[sC2H5] =  2.86109429e+02 * (
			T * (  2.69070170e+00 + T * (  4.35956660e-03
			+ T * (  1.47327957e-06 + T * (  2.33467577e-10
			+ T * -7.85554700e-13 ) ) ) ) +  1.28704040e+04 );
		cp[sC2H5] =  2.86109429e+02 * (
			 2.69070170e+00 + T * (  8.71913320e-03 
			+ T * (  4.41983870e-06 + T * (  9.33870310e-10
			+ T * -3.92777350e-12 ) ) ) );
		h[sC3H4] =  2.07536818e+02 * (
			T * (  2.66811280e+00 + T * (  8.08889400e-03
			+ T * (  3.24122083e-07 + T * ( -2.88316400e-09
			+ T *  1.13041972e-12 ) ) ) ) +  2.08009570e+04 );
		cp[sC3H4] =  2.07536818e+02 * (
			 2.66811280e+00 + T * (  1.61777880e-02 
			+ T * (  9.72366250e-07 + T * ( -1.15326560e-08
			+ T *  5.65209860e-12 ) ) ) );
		h[sC3H5] =  2.02443146e+02 * (
			T * ( -5.41004000e-01 + T * (  1.36420505e-02
			+ T * ( -3.21217763e-07 + T * ( -4.78236550e-09
			+ T *  1.96788350e-12 ) ) ) ) +  1.51303950e+04 );
		cp[sC3H5] =  2.02443146e+02 * (
			-5.41004000e-01 + T * (  2.72841010e-02 
			+ T * ( -9.63653290e-07 + T * ( -1.91294620e-08
			+ T *  9.83941750e-12 ) ) ) );
		h[sC3H6] =  1.97593517e+02 * (
			T * (  1.49330710e+00 + T * (  1.04625875e-02
			+ T * (  1.49559793e-06 + T * ( -4.17228025e-09
			+ T *  1.43162930e-12 ) ) ) ) +  1.07482640e+03 );
		cp[sC3H6] =  1.97593517e+02 * (
			 1.49330710e+00 + T * (  2.09251750e-02 
			+ T * (  4.48679380e-06 + T * ( -1.66891210e-08
			+ T *  7.15814650e-12 ) ) ) );
		h[sNXC3H7] =  1.92970803e+02 * (
			T * (  1.92253690e+00 + T * (  1.23946370e-02
			+ T * (  6.03416400e-07 + T * ( -4.45816450e-09
			+ T *  1.71659926e-12 ) ) ) ) +  9.71328130e+03 );
		cp[sNXC3H7] =  1.92970803e+02 * (
			 1.92253690e+00 + T * (  2.47892740e-02 
			+ T * (  1.81024920e-06 + T * ( -1.78326580e-08
			+ T *  8.58299630e-12 ) ) ) );
		h[s1XC4H8] =  1.48195138e+02 * (
			T * (  1.18113800e+00 + T * (  1.54266900e-02
			+ T * (  1.69550823e-06 + T * ( -6.16372200e-09
			+ T *  2.22203860e-12 ) ) ) ) -1.79040040e+03 );
		cp[s1XC4H8] =  1.48195138e+02 * (
			 1.18113800e+00 + T * (  3.08533800e-02 
			+ T * (  5.08652470e-06 + T * ( -2.46548880e-08
			+ T *  1.11101930e-11 ) ) ) );
		h[sPXC4H9] =  1.45579563e+02 * (
			T * (  2.18641570e+00 + T * (  1.65024830e-02
			+ T * (  1.01004027e-06 + T * ( -6.01018000e-09
			+ T *  2.26127540e-12 ) ) ) ) +  5.89791800e+03 );
		cp[sPXC4H9] =  1.45579563e+02 * (
			 2.18641570e+00 + T * (  3.30049660e-02 
			+ T * (  3.03012080e-06 + T * ( -2.40407200e-08
			+ T *  1.13063770e-11 ) ) ) );
		h[s1XC5H11] =  1.16876212e+02 * (
			T * (  2.19509120e+00 + T * (  2.06457545e-02
			+ T * (  2.19508063e-06 + T * ( -8.43787700e-09
			+ T *  3.08107980e-12 ) ) ) ) +  3.02183180e+03 );
		cp[s1XC5H11] =  1.16876212e+02 * (
			 2.19509120e+00 + T * (  4.12915090e-02 
			+ T * (  6.58524190e-06 + T * ( -3.37515080e-08
			+ T *  1.54053990e-11 ) ) ) );
		h[sC6H11] =  9.99944677e+01 * (
			T * ( -1.39422800e+00 + T * (  2.36034655e-02
			+ T * (  3.98680633e-06 + T * ( -1.04072395e-08
			+ T *  3.48067140e-12 ) ) ) ) -2.21779000e+03 );
		cp[sC6H11] =  9.99944677e+01 * (
			-1.39422800e+00 + T * (  4.72069310e-02 
			+ T * (  1.19604190e-05 + T * ( -4.16289580e-08
			+ T *  1.74033570e-11 ) ) ) );
		h[s1XC6H12] =  9.87967584e+01 * (
			T * (  1.96862030e+00 + T * (  2.38281155e-02
			+ T * (  2.20051243e-06 + T * ( -9.28704325e-09
			+ T *  3.38449260e-12 ) ) ) ) -7.71187890e+03 );
		cp[s1XC6H12] =  9.87967584e+01 * (
			 1.96862030e+00 + T * (  4.76562310e-02 
			+ T * (  6.60153730e-06 + T * ( -3.71481730e-08
			+ T *  1.69224630e-11 ) ) ) );
		h[s1XC7H15] =  8.38223611e+01 * (
			T * (  2.46798900e+00 + T * (  2.94440055e-02
			+ T * (  3.16893967e-06 + T * ( -1.21711593e-08
			+ T *  4.41048580e-12 ) ) ) ) -2.82337940e+03 );
		cp[s1XC7H15] =  8.38223611e+01 * (
			 2.46798900e+00 + T * (  5.88880110e-02 
			+ T * (  9.50681900e-06 + T * ( -4.86846370e-08
			+ T *  2.20524290e-11 ) ) ) );
		h[s2XC7H15] =  8.38223611e+01 * (
			T * (  2.46798900e+00 + T * (  2.94440055e-02
			+ T * (  3.16893967e-06 + T * ( -1.21711593e-08
			+ T *  4.41048580e-12 ) ) ) ) -2.82337940e+03 );
		cp[s2XC7H15] =  8.38223611e+01 * (
			 2.46798900e+00 + T * (  5.88880110e-02 
			+ T * (  9.50681900e-06 + T * ( -4.86846370e-08
			+ T *  2.20524290e-11 ) ) ) );
		h[sNXC7H16] =  8.29791014e+01 * (
			T * (  3.01495460e+00 + T * (  2.72286015e-02
			+ T * (  7.27089367e-06 + T * ( -1.35585278e-08
			+ T *  4.16174600e-12 ) ) ) ) -2.60033790e+04 );
		cp[sNXC7H16] =  8.29791014e+01 * (
			 3.01495460e+00 + T * (  5.44572030e-02 
			+ T * (  2.18126810e-05 + T * ( -5.42341110e-08
			+ T *  2.08087300e-11 ) ) ) );
		h[s1C7H15O2XC7H15O2] =  6.33763244e+01 * (
			T * (  3.70441770e+00 + T * ( -7.10937650e-04
			+ T * (  9.55679733e-07 + T * ( -3.00722125e-10
			+ T * -2.79093540e-15 ) ) ) ) -1.06407950e+03 );
		cp[s1C7H15O2XC7H15O2] =  6.33763244e+01 * (
			 3.70441770e+00 + T * ( -1.42187530e-03 
			+ T * (  2.86703920e-06 + T * ( -1.20288850e-09
			+ T * -1.39546770e-14 ) ) ) );
		h[s1HEOOHX2XC7H15O2] =  6.33763244e+01 * (
			T * (  3.70441770e+00 + T * ( -7.10937650e-04
			+ T * (  9.55679733e-07 + T * ( -3.00722125e-10
			+ T * -2.79093540e-15 ) ) ) ) -1.06407950e+03 );
		cp[s1HEOOHX2XC7H15O2] =  6.33763244e+01 * (
			 3.70441770e+00 + T * ( -1.42187530e-03 
			+ T * (  2.86703920e-06 + T * ( -1.20288850e-09
			+ T * -1.39546770e-14 ) ) ) );
		h[s2OO7OOH1XO2C7H15O2] =  5.09488326e+01 * (
			T * (  3.70441770e+00 + T * ( -7.10937650e-04
			+ T * (  9.55679733e-07 + T * ( -3.00722125e-10
			+ T * -2.79093540e-15 ) ) ) ) -1.06407950e+03 );
		cp[s2OO7OOH1XO2C7H15O2] =  5.09488326e+01 * (
			 3.70441770e+00 + T * ( -1.42187530e-03 
			+ T * (  2.86703920e-06 + T * ( -1.20288850e-09
			+ T * -1.39546770e-14 ) ) ) );
		h[s1X2HEOOHXO2C7H14OOH] =  5.09488326e+01 * (
			T * (  3.70441770e+00 + T * ( -7.10937650e-04
			+ T * (  9.55679733e-07 + T * ( -3.00722125e-10
			+ T * -2.79093540e-15 ) ) ) ) -1.06407950e+03 );
		cp[s1X2HEOOHXO2C7H14OOH] =  5.09488326e+01 * (
			 3.70441770e+00 + T * ( -1.42187530e-03 
			+ T * (  2.86703920e-06 + T * ( -1.20288850e-09
			+ T * -1.39546770e-14 ) ) ) );
		h[sOC7OOHXOC7H13OOH] =  5.68766332e+01 * (
			T * (  3.70441770e+00 + T * ( -7.10937650e-04
			+ T * (  9.55679733e-07 + T * ( -3.00722125e-10
			+ T * -2.79093540e-15 ) ) ) ) -1.06407950e+03 );
		cp[sOC7OOHXOC7H13OOH] =  5.68766332e+01 * (
			 3.70441770e+00 + T * ( -1.42187530e-03 
			+ T * (  2.86703920e-06 + T * ( -1.20288850e-09
			+ T * -1.39546770e-14 ) ) ) );
		h[s1X2HEPOOXOC7H13O] =  6.43654296e+01 * (
			T * (  3.70441770e+00 + T * ( -7.10937650e-04
			+ T * (  9.55679733e-07 + T * ( -3.00722125e-10
			+ T * -2.79093540e-15 ) ) ) ) -1.06407950e+03 );
		cp[s1X2HEPOOXOC7H13O] =  6.43654296e+01 * (
			 3.70441770e+00 + T * ( -1.42187530e-03 
			+ T * (  2.86703920e-06 + T * ( -1.20288850e-09
			+ T * -1.39546770e-14 ) ) ) );
		h[sCH2] =  5.92780550e+02 * (
			T * (  3.76223700e+00 + T * (  5.79909500e-04
			+ T * (  8.29861667e-08 + T * (  2.20020900e-10
			+ T * -1.46648700e-13 ) ) ) ) +  4.53679100e+04 );
		cp[sCH2] =  5.92780550e+02 * (
			 3.76223700e+00 + T * (  1.15981900e-03 
			+ T * (  2.48958500e-07 + T * (  8.80083600e-10
			+ T * -7.33243500e-13 ) ) ) );
		h[sCH3OH] =  2.59482554e+02 * (
			T * (  2.66011500e+00 + T * (  3.67075400e-03
			+ T * (  2.39001700e-06 + T * ( -2.19829850e-09
			+ T *  4.78114000e-13 ) ) ) ) -2.53534800e+04 );
		cp[sCH3OH] =  2.59482554e+02 * (
			 2.66011500e+00 + T * (  7.34150800e-03 
			+ T * (  7.17005100e-06 + T * ( -8.79319400e-09
			+ T *  2.39057000e-12 ) ) ) );
		h[sC] =  6.92284763e+02 * (
			T * (  2.55423955e+00 + T * ( -1.60768862e-04
			+ T * (  2.44597415e-07 + T * ( -1.83058722e-10
			+ T *  5.33042892e-14 ) ) ) ) +  8.54438832e+04 );
		cp[sC] =  6.92284763e+02 * (
			 2.55423955e+00 + T * ( -3.21537724e-04 
			+ T * (  7.33792245e-07 + T * ( -7.32234889e-10
			+ T *  2.66521446e-13 ) ) ) );
		h[sC2H] =  3.32201534e+02 * (
			T * (  2.73770400e+00 + T * (  4.02422300e-03
			+ T * ( -3.08143667e-06 + T * (  1.63131475e-09
			+ T * -3.87916000e-13 ) ) ) ) +  6.68381300e+04 );
		cp[sC2H] =  3.32201534e+02 * (
			 2.73770400e+00 + T * (  8.04844600e-03 
			+ T * ( -9.24431000e-06 + T * (  6.52525900e-09
			+ T * -1.93958000e-12 ) ) ) );
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
	return 62;
}

int GetNReactions( void )
{
	return rEnd;
}

double MAX_C(double X1, double X2)
{
  return ( (X1 > X2) ? X1 : X2 );
}
