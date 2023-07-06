#include "nHeptane.allstarnew_oksred.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
static double GetLindRateCoeff( double temp, double pressure
			, double k0, double kInf, double fc, double conc );

void ComputeProductionRates( double *cdot, double *w, double *k
			, double *c, double *M, double temp, double pressure );
void ComputeSteadyStates( double *k, double *c, double *M );
double MAX_C(double X1, double X2);
double CatchZero( double a );
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

	int	nSpec = 44;
	int	nSpecIn = 25;
	double	kTroe0, kTroeInf, fcTroe;
	double	RGAS = 8314.34;
	double	lgt = log( temp );
	double	rt = RGAS * temp;

	M[mM1] = 0.4 * c[sN2] + 0.4 * c[sO2] + c[sH] + c[sOH]
		 + c[sO] + c[sH2] + 6.5 * c[sH2O] + c[sHO2]
		 + c[sH2O2] + 0.75 * c[sCO] + c[sCO2] + c[s3XCH2]
		 + c[sC2H2] + c[sCH3] + c[sC2H4] + c[sCH2O]
		 + 3 * c[sCH4] + c[sC3H4] + c[sC3H5] + c[sC3H6]
		 + c[s1XC4H8] + c[s1XC6H12] + 3 * c[sNXC7H16] + c[s1C7H15O2XC7H15O2]
		 + c[sOC7OOHXOC7H13OOH];



	k[r1f] = 2.0000000000E+11 * exp( -70300000 / rt );
	k[r1b] = 4.0674859656E+08 * exp( 0.433091 * lgt + 1978606.19 / rt );
	k[r2f] = 5.0600000000E+01 * exp( 2.67 * lgt - 26300000 / rt );
	k[r2b] = 2.9161028451E+01 * exp( 2.63868 * lgt - 18754035.03 / rt );
	k[r3f] = 1.0000000000E+05 * exp( 1.6 * lgt - 13800000 / rt );
	k[r3b] = 1.0374701530E+06 * exp( 1.50368 * lgt - 78257627.95 / rt );
	k[r4f] = 1.5000000000E+06 * exp( 1.14 * lgt - 420000 / rt );
	k[r4b] = 2.7003157570E+07 * exp( 1.075 * lgt - 72423592.92 / rt );
	k[r5f] = 1.8000000000E+12 * exp( -1 * lgt );
	k[r5b] = 4.2896464802E+14 * exp( -0.65384 * lgt - 433084285.8 / rt );
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
	k[r22b] = 5.4267465236E+16 * exp( -0.965305 * lgt - 515755732.9 / rt );
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
	k[r46b] = 3.6820837010E+04 * exp( 1.50092 * lgt - 68348947.29 / rt );
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
	k[r58b] = 2.7601829058E+11 * exp( -0.160882 * lgt - 113532166.5 / rt );
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
	k[r62b] = 5.7658060013E+05 * exp( 0.321726 * lgt - 18065716.66 / rt );
	k[r63f] = 1.8000000000E+10;
	k[r63b] = 4.9466538837E+07 * exp( 0.667886 * lgt - 346150002.4 / rt );
	k[r64f] = 4.0000000000E+07 * exp( -8900000 / rt );
	k[r64b] = 7.0716715146E+05 * exp( 0.322645 * lgt - 118162756 / rt );
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
	k[r111b] = 1.9767804329E-05 * exp( 4.21906 * lgt - 196680018.1 / rt );
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


	ComputeSteadyStates( k, c, M );

	double	cTot = pressure / ( RGAS * temp );
	int iNow;
	for ( iNow = nSpecIn; iNow < nSpec; ++iNow ) {
		if ( c[iNow] > cTot ) {
			c[iNow] = cTot;
		}
	}
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
	w[r58b] = k[r58b] * c[sOH] * c[sCH3O];
	w[r59f] = k[r59f] * c[sCH3] * c[sHO2];
	w[r59b] = k[r59b] * c[sO2] * c[sCH4];
	w[r60] = k[r60] * c[sCH3] * c[sCH3];
	w[r36f] = k[r36f] * c[sCH3] * c[sCH3];
	w[r36b] = k[r36b] * c[sC2H6];
	w[r62f] = k[r62f] * c[sCH3O] * M[mM1];
	w[r62b] = k[r62b] * c[sH] * c[sCH2O] * M[mM1];
	w[r63f] = k[r63f] * c[sCH3O] * c[sH];
	w[r63b] = k[r63b] * c[sH2] * c[sCH2O];
	w[r64f] = k[r64f] * c[sCH3O] * c[sO2];
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


	cdot[sN2] = 0.0;

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
		+ w[ri67b] - w[ri69];

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
		- w[r206b] - w[r357] - w[rn197] - w[rn198];

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
		- w[rn201] - w[rn202] + w[ri71] + w[ri72];

	cdot[sO] = w[r1f] - w[r1b] - w[r2f] + w[r2b]
		+ w[r4f] - w[r4b] - 2 * w[r6f] + 2 * w[r6b]
		+ w[r11f] - w[r11b] - w[r12f] + w[r12b]
		- w[r22f] + w[r22b] + w[r25] - w[r47f]
		+ w[r47b] - w[r53f] + w[r53b] - w[r106]
		- w[r111f] + w[r111b] - w[r112f] + w[r112b];

	cdot[sH2] = - w[r2f] + w[r2b] - w[r3f] + w[r3b]
		+ w[r5f] - w[r5b] + w[r10f] - w[r10b]
		+ w[r35f] - w[r35b] + w[r37f] - w[r37b]
		+ w[r41f] - w[r41b] - w[r44f] + w[r44b]
		+ w[r46f] - w[r46b] + w[r60] + w[r63f]
		- w[r63b] + w[r70f] - w[r70b] + w[r84f]
		- w[r84b] + w[r129f] - w[r129b] + w[r131f]
		- w[r131b] + w[r164] + w[r192f] - w[r192b]
		+ w[r196f] - w[r196b] + w[r357] + w[rn197]
		+ w[rn198];

	cdot[sH2O] = w[r3f] - w[r3b] + w[r4f] - w[r4b]
		+ w[r7f] - w[r7b] + w[r11f] - w[r11b]
		+ w[r13f] - w[r13b] + w[r19f] - w[r19b]
		- w[r27f] + w[r27b] + w[r48f] - w[r48b]
		+ w[r86f] - w[r86b] + w[r134f] - w[r134b]
		+ w[r166] + w[r202f] - w[r202b] + w[r385]
		+ w[rn201] + w[rn202];

	cdot[sHO2] = w[r8f] - w[r8b] - w[r9f] + w[r9b]
		- w[r10f] + w[r10b] - w[r11f] + w[r11b]
		- w[r12f] + w[r12b] - w[r13f] + w[r13b]
		- 2 * w[r14f] + 2 * w[r14b] + w[r19f] - w[r19b]
		- w[r21f] + w[r21b] + w[r33f] - w[r33b]
		- w[r49f] + w[r49b] - w[r58f] + w[r58b]
		- w[r59f] + w[r59b] + w[r64f] - w[r64b]
		+ w[r71f] - w[r71b] + w[r149f] - w[r149b]
		+ w[r193f] - w[r193b] + w[r207f] - w[r207b]
		- w[rn209] - w[rn210] + w[rn217] + w[rn218];

	cdot[sH2O2] = w[r14f] - w[r14b] + w[r15f] - w[r15b]
		- w[r19f] + w[r19b] + w[r49f] - w[r49b]
		+ w[rn209] + w[rn210];

	cdot[sCO] = - w[r20f] + w[r20b] - w[r21f] + w[r21b]
		- w[r22f] + w[r22b] + w[r26] + w[r28f]
		- w[r28b] + w[r33f] - w[r33b] + w[r40f]
		- w[r40b] + w[r43f] - w[r43b] + w[r105f]
		- w[r105b] + 2 * w[r106] + w[r111f] - w[r111b]
		+ w[ri73];

	cdot[sCO2] = w[r20f] - w[r20b] + w[r21f] - w[r21b]
		+ w[r22f] - w[r22b] - w[r26] + w[r41f]
		- w[r41b];

	cdot[s3XCH2] = - w[r35f] + w[r35b] - 2 * w[r37f] + 2 * w[r37b]
		- w[r39f] + w[r39b] - w[r40f] + w[r40b]
		- w[r41f] + w[r41b] + w[r42f] - w[r42b]
		+ w[r105f] - w[r105b] + w[r111f] - w[r111b];

	cdot[sC2H2] = w[r37f] - w[r37b] - w[r109f] + w[r109b]
		- w[r111f] + w[r111b] - w[r112f] + w[r112b]
		+ w[ra51f] - w[ra51b] + w[r129f] - w[r129b];

	cdot[sCH3] = - w[r39f] + w[r39b] + w[r44f] - w[r44b]
		- w[r53f] + w[r53b] - w[ra34f] + w[ra34b]
		- w[r55] - w[r57] - w[r58f] + w[r58b]
		- w[r59f] + w[r59b] - 2 * w[r60] - 2 * w[r36f]
		+ 2 * w[r36b] + w[r84f] - w[r84b] + w[r86f]
		- w[r86b] + 2 * w[r146f] - 2 * w[r146b] - w[r170]
		+ w[r195f] - w[r195b] - w[r203] + w[r205f]
		- w[r205b] + w[r255f] - w[r255b] + w[rn194];

	cdot[sC2H4] = w[r39f] - w[r39b] + w[r60] - w[r129f]
		+ w[r129b] - w[r131f] + w[r131b] - w[r134f]
		+ w[r134b] + w[ra58f] - w[ra58b] + w[r149f]
		- w[r149b] + w[r190] + w[r205f] - w[r205b]
		+ w[r292] + w[r349] + w[rn192];

	cdot[sCH2O] = - w[r45f] + w[r45b] - w[r46f] + w[r46b]
		- w[r47f] + w[r47b] - w[r48f] + w[r48b]
		- w[r49f] + w[r49b] + w[r53f] - w[r53b]
		+ w[r57] + w[r62f] - w[r62b] + w[r63f]
		- w[r63b] + w[r64f] - w[r64b] + w[r69f]
		- w[r69b] + w[r70f] - w[r70b] + w[r71f]
		- w[r71b] + w[ra125] + w[r189] + w[r200f]
		- w[r200b] + w[r265f] - w[r265b] + w[ri73];

	cdot[sCH4] = w[ra34f] - w[ra34b] + w[r59f] - w[r59b]
		- w[r84f] + w[r84b] - w[r86f] + w[r86b]
		+ w[r170] + w[r203];

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

	cdot[s1XC4H8] = - w[r255f] + w[r255b] - w[r265f] + w[r265b];

	cdot[s1XC6H12] = - w[r355] - w[r357] - w[r385] + w[rn194];

	cdot[sNXC7H16] = - w[rn196] - w[rn197] - w[rn198] - w[rn201]
		- w[rn202] - w[rn209] - w[rn210] - w[rn217]
		- w[rn218];

	cdot[s1C7H15O2XC7H15O2] = w[ri66f] - w[ri66b] + w[ri67f] - w[ri67b]
		- w[ri68];

	cdot[sOC7OOHXOC7H13OOH] = w[ri71] - w[ri72];

}

double CatchZero( double a )
{
	return ( a == 0.0 ) ? 1.0e-20 : a;
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
	W[s3XCH2] =  1.40260000e+01;
	W[sC2H2] =  2.60360000e+01;
	W[sCH3] =  1.50340000e+01;
	W[sC2H4] =  2.80520000e+01;
	W[sCH2O] =  3.00260000e+01;
	W[sCH4] =  1.60420000e+01;
	W[sC3H4] =  4.00620000e+01;
	W[sC3H5] =  4.10700000e+01;
	W[sC3H6] =  4.20780000e+01;
	W[s1XC4H8] =  5.61040000e+01;
	W[s1XC6H12] =  8.41560000e+01;
	W[sNXC7H16] =  1.00198000e+02;
	W[s1C7H15O2XC7H15O2] =  1.31190000e+02;
	W[sOC7OOHXOC7H13OOH] =  1.46182000e+02;
	W[sCH] =  1.30180000e+01;
	W[sHCO] =  2.90180000e+01;
	W[sCH2OH] =  3.10340000e+01;
	W[s1XCH2] =  1.40260000e+01;
	W[sCH3O] =  3.10340000e+01;
	W[sC2H6] =  3.00680000e+01;
	W[sHCCO] =  4.10280000e+01;
	W[sC2H3] =  2.70440000e+01;
	W[sC2H5] =  2.90600000e+01;
	W[sNXC3H7] =  4.30860000e+01;
	W[sPXC4H9] =  5.71120000e+01;
	W[s1XC5H11] =  7.11380000e+01;
	W[sC6H11] =  8.31480000e+01;
	W[s1XC7H15] =  9.91900000e+01;
	W[s2XC7H15] =  9.91900000e+01;
	W[s1HEOOHX2XC7H15O2] =  1.31190000e+02;
	W[s2OO7OOH1XO2C7H15O2] =  1.63190000e+02;
	W[s1X2HEOOHXO2C7H14OOH] =  1.63190000e+02;
	W[s1X2HEPOOXOC7H13O] =  1.29174000e+02;
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
	strcpy( names[s3XCH2], "3-CH2               " );
	strcpy( names[sC2H2], "C2H2                " );
	strcpy( names[sCH3], "CH3                 " );
	strcpy( names[sC2H4], "C2H4                " );
	strcpy( names[sCH2O], "CH2O                " );
	strcpy( names[sCH4], "CH4                 " );
	strcpy( names[sC3H4], "C3H4                " );
	strcpy( names[sC3H5], "C3H5                " );
	strcpy( names[sC3H6], "C3H6                " );
	strcpy( names[s1XC4H8], "1-C4H8              " );
	strcpy( names[s1XC6H12], "1-C6H12             " );
	strcpy( names[sNXC7H16], "N-C7H16             " );
	strcpy( names[s1C7H15O2XC7H15O2], "1C7H15O2-C7H15O2    " );
	strcpy( names[sOC7OOHXOC7H13OOH], "OC7OOH-OC7H13OOH    " );
	strcpy( names[sCH], "CH                  " );
	strcpy( names[sHCO], "HCO                 " );
	strcpy( names[sCH2OH], "CH2OH               " );
	strcpy( names[s1XCH2], "1-CH2               " );
	strcpy( names[sCH3O], "CH3O                " );
	strcpy( names[sC2H6], "C2H6                " );
	strcpy( names[sHCCO], "HCCO                " );
	strcpy( names[sC2H3], "C2H3                " );
	strcpy( names[sC2H5], "C2H5                " );
	strcpy( names[sNXC3H7], "N-C3H7              " );
	strcpy( names[sPXC4H9], "P-C4H9              " );
	strcpy( names[s1XC5H11], "1-C5H11             " );
	strcpy( names[sC6H11], "C6H11               " );
	strcpy( names[s1XC7H15], "1-C7H15             " );
	strcpy( names[s2XC7H15], "2-C7H15             " );
	strcpy( names[s1HEOOHX2XC7H15O2], "1HEOOH-2-C7H15O2    " );
	strcpy( names[s2OO7OOH1XO2C7H15O2], "2OO7OOH1-O2C7H15O2  " );
	strcpy( names[s1X2HEOOHXO2C7H14OOH], "1-2HEOOH-O2C7H14OOH " );
	strcpy( names[s1X2HEPOOXOC7H13O], "1-2HEPOO-OC7H13O    " );
}


void ComputeSteadyStates(double *k, double *c, double *M)
{
// c[s1HEOOHX2XC7H15O2] (uncoupled) --------------------

		c[s1HEOOHX2XC7H15O2] = ( 0
			+ k[ri68] * c[s1C7H15O2XC7H15O2]
			) / ( CatchZero(
			+ k[ri69] * c[sO2]
			) );


// c[s1XCH2] (uncoupled) --------------------

		c[s1XCH2] = ( 0
			+ k[r42b] * c[s3XCH2] * M[mM1]
			+ k[r43b] * c[sH] * c[sOH] * c[sCO]
			+ k[r44b] * c[sH] * c[sCH3]
			) / ( CatchZero(
			+ k[r42f] * M[mM1]
			+ k[r43f] * c[sO2]
			+ k[r44f] * c[sH2]
			) );


// c[sCH3O] (uncoupled) --------------------

		c[sCH3O] = ( 0
			+ k[r55] * c[sCH3] * c[sOH]
			+ k[r58f] * c[sCH3] * c[sHO2]
			+ k[r62b] * c[sH] * c[sCH2O] * M[mM1]
			+ k[r63b] * c[sH2] * c[sCH2O]
			+ k[r64b] * c[sHO2] * c[sCH2O]
			) / ( CatchZero(
			+ k[r58b] * c[sOH]
			+ k[r62f] * M[mM1]
			+ k[r63f] * c[sH]
			+ k[r64f] * c[sO2]
			) );


// c[sC2H6] (uncoupled) --------------------

		c[sC2H6] = ( 0
			+ k[r36f] * c[sCH3] * c[sCH3]
			) / ( CatchZero(
			+ k[r36b]
			+ k[r164] * c[sH]
			+ k[r166] * c[sOH]
			+ k[r170] * c[sCH3]
			) );


// c[sHCCO] (uncoupled) --------------------

		c[sHCCO] = ( 0
			+ k[r105b] * c[sCO] * c[s3XCH2]
			+ k[r109f] * c[sC2H2] * c[sO2]
			+ k[r112f] * c[sC2H2] * c[sO]
			) / ( CatchZero(
			+ k[r105f] * c[sH]
			+ k[r106] * c[sO]
			+ k[r109b] * c[sOH]
			+ k[r112b] * c[sH]
			) );


// c[sC2H3] (uncoupled) --------------------

		c[sC2H3] = ( 0
			+ k[ra51b] * c[sH] * c[sC2H2]
			+ k[r131f] * c[sC2H4] * c[sH]
			+ k[r134f] * c[sC2H4] * c[sOH]
			+ k[r189] * c[sC3H4] * c[sOH]
			+ k[r195f] * c[sC3H6]
			) / ( CatchZero(
			+ k[ra51f]
			+ k[ra125] * c[sO2]
			+ k[r131b] * c[sH2]
			+ k[r134b] * c[sH2O]
			+ k[r195b] * c[sCH3]
			) );


// c[s1XC7H15] c[s2XC7H15] (coupled)  --------------------

// Primary denominators-----------------------

	 double X1XC7H15_denom1 = CatchZero( 0
			+ k[rn192]
			+ k[r477]
			+ k[ri66f] * c[sO2]
			);

	 double X2XC7H15_denom1 = CatchZero( 0
			+ k[rn193]
			+ k[rn194]
			+ k[r478]
			+ k[ri67f] * c[sO2]
			);


// Primary constant parts -----------------------

	 double X1XC7H15_ct1 = ( 0
			+ k[rn197] * c[sNXC7H16] * c[sH]
			+ k[rn201] * c[sNXC7H16] * c[sOH]
			+ k[rn209] * c[sNXC7H16] * c[sHO2]
			+ k[rn217] * c[sNXC7H16] * c[sO2]
			+ k[ri66b] * c[s1C7H15O2XC7H15O2]
			);

	 double X2XC7H15_ct1 = ( 0
			+ k[rn198] * c[sNXC7H16] * c[sH]
			+ k[rn202] * c[sNXC7H16] * c[sOH]
			+ k[rn210] * c[sNXC7H16] * c[sHO2]
			+ k[rn218] * c[sNXC7H16] * c[sO2]
			+ k[ri67b] * c[s1C7H15O2XC7H15O2]
			);


// X1XC7H15 ---------------------------------------

	 double X1XC7H15_denom2 = CatchZero( X1XC7H15_denom1
			);

	 double X1XC7H15_ct2 = ( X1XC7H15_ct1
			) / X1XC7H15_denom2;



	 double X1XC7H15_X2XC7H15 = ( 0 + k[r478] 
			) / X1XC7H15_denom2 ;



	 double X1XC7H15_X2XC7H15_coeff = ( 0 + k[r477] 
			);


// X2XC7H15 ---------------------------------------

	 double X2XC7H15_denom2 = CatchZero( X2XC7H15_denom1
			- X1XC7H15_X2XC7H15_coeff * X1XC7H15_X2XC7H15
			);

	 double X2XC7H15_ct2 = ( X2XC7H15_ct1
			+ X1XC7H15_X2XC7H15_coeff * X1XC7H15_ct2
			) / X2XC7H15_denom2;






// Reconstruction ------------------------------------

	c[s2XC7H15] = ( X2XC7H15_ct2
			) ;

	c[s1XC7H15] = ( X1XC7H15_ct2
			+ X1XC7H15_X2XC7H15 * c[s2XC7H15]
			) ;

// c[sPXC4H9] (uncoupled) --------------------

		c[sPXC4H9] = ( 0
			+ k[rn193] * c[s2XC7H15]
			+ k[rn196] * c[sNXC7H16]
			) / ( CatchZero(
			+ k[r292]
			) );


// c[sCH2OH] c[sCH] (coupled)  --------------------

// Primary denominators-----------------------

	 double CH2OH_denom1 = CatchZero( 0
			+ k[r27b]
			+ k[r69f] * M[mM1]
			+ k[r70f] * c[sH]
			+ k[r71f] * c[sO2]
			);

	 double CH_denom1 = CatchZero( 0
			+ k[r25] * c[sO2]
			+ k[r26] * c[sCO2]
			+ k[r27f] * c[sH2O]
			+ k[r35b] * c[sH2]
			);


// Primary constant parts -----------------------

	 double CH2OH_ct1 = ( 0
			+ k[r69b] * c[sH] * c[sCH2O] * M[mM1]
			+ k[r70b] * c[sH2] * c[sCH2O]
			+ k[r71b] * c[sHO2] * c[sCH2O]
			);

	 double CH_ct1 = ( 0
			+ k[r35f] * c[s3XCH2] * c[sH]
			);


// CH2OH ---------------------------------------

	 double CH2OH_denom2 = CatchZero( CH2OH_denom1
			);

	 double CH2OH_ct2 = ( CH2OH_ct1
			) / CH2OH_denom2;



	 double CH2OH_CH = ( 0 + k[r27f]  * c[sH2O]
			) / CH2OH_denom2 ;



	 double CH2OH_CH_coeff = ( 0 + k[r27b] 
			);


// CH ---------------------------------------

	 double CH_denom2 = CatchZero( CH_denom1
			- CH2OH_CH_coeff * CH2OH_CH
			);

	 double CH_ct2 = ( CH_ct1
			+ CH2OH_CH_coeff * CH2OH_ct2
			) / CH_denom2;






// Reconstruction ------------------------------------

	c[sCH] = ( CH_ct2
			) ;

	c[sCH2OH] = ( CH2OH_ct2
			+ CH2OH_CH * c[sCH]
			) ;

// c[sC6H11] (uncoupled) --------------------

		c[sC6H11] = ( 0
			+ k[r357] * c[s1XC6H12] * c[sH]
			+ k[r385] * c[s1XC6H12] * c[sOH]
			) / ( CatchZero(
			+ k[r353]
			) );


// c[sC2H5] (uncoupled) --------------------

		c[sC2H5] = ( 0
			+ k[ra58b] * c[sH] * c[sC2H4]
			+ k[r146b] * c[sCH3] * c[sCH3]
			+ k[r149b] * c[sHO2] * c[sC2H4]
			+ k[r164] * c[sC2H6] * c[sH]
			+ k[r166] * c[sC2H6] * c[sOH]
			+ k[r170] * c[sC2H6] * c[sCH3]
			+ k[r200f] * c[sC3H6] * c[sOH]
			+ k[r292] * c[sPXC4H9]
			) / ( CatchZero(
			+ k[ra58f]
			+ k[r146f] * c[sH]
			+ k[r149f] * c[sO2]
			+ k[r200b] * c[sCH2O]
			) );


// c[s2OO7OOH1XO2C7H15O2] (uncoupled) --------------------

		c[s2OO7OOH1XO2C7H15O2] = ( 0
			+ k[ri69] * c[s1HEOOHX2XC7H15O2] * c[sO2]
			) / ( CatchZero(
			+ k[ri70]
			) );


// c[s1X2HEOOHXO2C7H14OOH] (uncoupled) --------------------

		c[s1X2HEOOHXO2C7H14OOH] = ( 0
			+ k[ri70] * c[s2OO7OOH1XO2C7H15O2]
			) / ( CatchZero(
			+ k[ri71]
			) );


// c[sHCO] (uncoupled) --------------------

		c[sHCO] = ( 0
			+ k[r25] * c[sCH] * c[sO2]
			+ k[r26] * c[sCH] * c[sCO2]
			+ k[r28b] * c[sH] * c[sCO] * M[mM1]
			+ k[r33b] * c[sHO2] * c[sCO]
			+ k[r45f] * c[sCH2O] * M[mM1]
			+ k[r46f] * c[sCH2O] * c[sH]
			+ k[r47f] * c[sCH2O] * c[sO]
			+ k[r48f] * c[sCH2O] * c[sOH]
			+ k[r49f] * c[sCH2O] * c[sHO2]
			+ k[ra125] * c[sC2H3] * c[sO2]
			+ k[r190] * c[sC3H4] * c[sOH]
			) / ( CatchZero(
			+ k[r28f] * M[mM1]
			+ k[r33f] * c[sO2]
			+ k[r45b] * c[sH] * M[mM1]
			+ k[r46b] * c[sH2]
			+ k[r47b] * c[sOH]
			+ k[r48b] * c[sH2O]
			+ k[r49b] * c[sH2O2]
			) );


// c[s1X2HEPOOXOC7H13O] (uncoupled) --------------------

		c[s1X2HEPOOXOC7H13O] = ( 0
			+ k[ri72] * c[sOC7OOHXOC7H13OOH]
			) / ( CatchZero(
			+ k[ri73]
			) );


// c[s1XC5H11] (uncoupled) --------------------

		c[s1XC5H11] = ( 0
			+ k[rn192] * c[s1XC7H15]
			+ k[ri73] * c[s1X2HEPOOXOC7H13O]
			) / ( CatchZero(
			+ k[r349]
			) );


// c[sNXC3H7] (uncoupled) --------------------

		c[sNXC3H7] = ( 0
			+ k[r205b] * c[sC2H4] * c[sCH3]
			+ k[r206b] * c[sC3H6] * c[sH]
			+ k[r207b] * c[sHO2] * c[sC3H6]
			+ k[r265f] * c[s1XC4H8] * c[sOH]
			+ k[r349] * c[s1XC5H11]
			+ k[r355] * c[s1XC6H12]
			+ k[rn196] * c[sNXC7H16]
			) / ( CatchZero(
			+ k[r205f]
			+ k[r206f]
			+ k[r207f] * c[sO2]
			+ k[r265b] * c[sCH2O]
			) );


}

void ComputeThermoData( double *h, double *cp, double T )
{
/*
	This function computes enthalpy 'h' and heat capacity 'cp' as
	function of temperature 'T' for all non steady state species
	in units [J/kg] and [J/kg K], respectively.
	The parameter h and cp should provide workspace of length 25 */


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
		h[s1XC4H8] =  1.48195138e+02 * (
			T * (  2.05358410e+00 + T * (  1.71752535e-02
			+ T * ( -5.29439900e-06 + T * (  8.27241550e-10
			+ T * -5.07220900e-14 ) ) ) ) -2.13972310e+03 );
		cp[s1XC4H8] =  1.48195138e+02 * (
			 2.05358410e+00 + T * (  3.43505070e-02 
			+ T * ( -1.58831970e-05 + T * (  3.30896620e-09
			+ T * -2.53610450e-13 ) ) ) );
		h[s1XC6H12] =  9.87967584e+01 * (
			T * (  1.86636350e+01 + T * (  1.04857255e-02
			+ T * ( -1.03609363e-06 + T * ( -1.71629045e-10
			+ T *  3.20472160e-14 ) ) ) ) -1.35908950e+04 );
		cp[s1XC6H12] =  9.87967584e+01 * (
			 1.86636350e+01 + T * (  2.09714510e-02 
			+ T * ( -3.10828090e-06 + T * ( -6.86516180e-10
			+ T *  1.60236080e-13 ) ) ) );
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
		h[sOC7OOHXOC7H13OOH] =  5.68766332e+01 * (
			T * (  2.85328990e+00 + T * (  8.01106400e-04
			+ T * ( -2.09789643e-07 + T * (  2.86025550e-11
			+ T * -1.56114930e-15 ) ) ) ) -8.90080930e+02 );
		cp[sOC7OOHXOC7H13OOH] =  5.68766332e+01 * (
			 2.85328990e+00 + T * (  1.60221280e-03 
			+ T * ( -6.29368930e-07 + T * (  1.14410220e-10
			+ T * -7.80574650e-15 ) ) ) );
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
		h[s1XC4H8] =  1.48195138e+02 * (
			T * (  1.18113800e+00 + T * (  1.54266900e-02
			+ T * (  1.69550823e-06 + T * ( -6.16372200e-09
			+ T *  2.22203860e-12 ) ) ) ) -1.79040040e+03 );
		cp[s1XC4H8] =  1.48195138e+02 * (
			 1.18113800e+00 + T * (  3.08533800e-02 
			+ T * (  5.08652470e-06 + T * ( -2.46548880e-08
			+ T *  1.11101930e-11 ) ) ) );
		h[s1XC6H12] =  9.87967584e+01 * (
			T * (  1.96862030e+00 + T * (  2.38281155e-02
			+ T * (  2.20051243e-06 + T * ( -9.28704325e-09
			+ T *  3.38449260e-12 ) ) ) ) -7.71187890e+03 );
		cp[s1XC6H12] =  9.87967584e+01 * (
			 1.96862030e+00 + T * (  4.76562310e-02 
			+ T * (  6.60153730e-06 + T * ( -3.71481730e-08
			+ T *  1.69224630e-11 ) ) ) );
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
		h[sOC7OOHXOC7H13OOH] =  5.68766332e+01 * (
			T * (  3.70441770e+00 + T * ( -7.10937650e-04
			+ T * (  9.55679733e-07 + T * ( -3.00722125e-10
			+ T * -2.79093540e-15 ) ) ) ) -1.06407950e+03 );
		cp[sOC7OOHXOC7H13OOH] =  5.68766332e+01 * (
			 3.70441770e+00 + T * ( -1.42187530e-03 
			+ T * (  2.86703920e-06 + T * ( -1.20288850e-09
			+ T * -1.39546770e-14 ) ) ) );
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
	return 25;
}

int GetNReactions( void )
{
	return rEnd;
}

double MAX_C(double X1, double X2)
{
  return ( (X1 > X2) ? X1 : X2 );
}
