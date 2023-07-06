#include "nDaM.allstar4.h"

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

	int	nSpec = 93;
	int	nSpecIn = 93;
	double	kTroe0, kTroeInf, fcTroe;
	double	RGAS = 8314.34;
	double	lgt = log( temp );
	double	rt = RGAS * temp;

	M[mM1] = 0.4 * c[sN2] + 0.4 * c[sO2] + c[sH2] + 6.5 * c[sH2O]
		 + 0.75 * c[sCO] + 3 * c[sCH4] + 3 * c[sNXC10H22] + 3 * c[sA2CH3XC11H10];



	k[r1f] = 2.0000000000E+11 * exp( -70300000 / rt );
	k[r2f] = 5.0600000000E+01 * exp( 2.67 * lgt - 26300000 / rt );
	k[r3f] = 1.0000000000E+05 * exp( 1.6 * lgt - 13800000 / rt );
	k[r4f] = 1.5000000000E+06 * exp( 1.14 * lgt - 420000 / rt );
	k[r5f] = 2.3000000000E+12 * exp( -0.8 * lgt );
	k[r6f] = 1.5000000000E+11 * exp( -4200000 / rt );
	k[r7f] = 2.5000000000E+10 * exp( -2900000 / rt );
	k[r8f] = 6.0000000000E+10;
	k[r9f] = 3.0000000000E+10 * exp( -7200000 / rt );
	k[r10f] = 1.8000000000E+10 * exp( 1700000 / rt );
	k[r11f] = 2.5000000000E+08 * exp( 5200000 / rt );
	k[r12f] = 2.0250000000E+17 * exp( -2 * lgt );
	k[r13f] = 1.7000000000E+09 * exp( -15700000 / rt );
	k[r14f] = 1.0000000000E+10 * exp( -15000000 / rt );
	k[r15f] = 2.8000000000E+10 * exp( -26800000 / rt );
	k[r16f] = 5.4000000000E+09 * exp( -4200000 / rt );
	k[r17f] = 1.8000000000E+12 * exp( -1 * lgt );
	k[r18f] = 2.2000000000E+16 * exp( -2 * lgt );
	k[r19f] = 2.9000000000E+11 * exp( -1 * lgt );
	k[r20f] = 6.0000000000E+03 * exp( 1.5 * lgt + 3100000 / rt );
	k[r21] = 1.5000000000E+11 * exp( -98700000 / rt );
	k[r22] = 7.1000000000E+07 * exp( 19000000 / rt );
	k[r24] = 3.3000000000E+10;
	k[r25] = 3.0000000000E+10;
	k[r26] = 3.4000000000E+09 * exp( -2900000 / rt );
	k[r27f] = 5.7000000000E+09 * exp( 3200000 / rt );
	k[r28f] = 7.1000000000E+11 * exp( -70300000 / rt );
	k[r29] = 9.0000000000E+10;
	k[r32] = 1.0000000000E+11;
	k[r35] = 3.0000000000E+09;
	k[r37f] = 6.0000000000E+09 * exp( 7500000 / rt );
	k[r39] = 1.3000000000E+10 * exp( -6200000 / rt );
	k[r40] = 1.2000000000E+10 * exp( -6200000 / rt );
	k[r42] = 1.2000000000E+10 * exp( -3400000 / rt );
	k[r43] = 1.1000000000E+11 * exp( -3400000 / rt );
	k[r44] = 4.2000000000E+10;
	k[r45f] = 1.2000000000E+10;
	k[r46f] = 7.2000000000E+10;
	k[r48] = 5.0000000000E+13 * exp( -320000000 / rt );
	k[r49] = 2.3000000000E+07 * exp( 1.05 * lgt - 13700000 / rt );
	k[r50] = 4.1500000000E+08 * exp( 0.57 * lgt - 11600000 / rt );
	k[r51] = 3.4000000000E+06 * exp( 1.2 * lgt + 1900000 / rt );
	k[r54f] = 6.0000000000E+10 * exp( -170700000 / rt );
	k[r56f] = 6.0000000000E+10 * exp( -63000000 / rt );
	k[r57] = 8.4300000000E+10;
	k[r58] = 3.3000000000E+08 * exp( -37400000 / rt );
	k[r59] = 2.2600000000E+11 * exp( -64800000 / rt );
	k[r60] = 4.7500000000E+13 * exp( -0.13 * lgt - 88000000 / rt );
	k[r61] = 1.8000000000E+10;
	kTroe0 = 6.2570000000E+17 * exp( -1.8 * lgt );
	kTroeInf = 2.1080000000E+11;
	fcTroe = 0.577 * exp( -temp / 2370 );
	k[r62f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, -1.0 );
	kTroe0 = 1.2720000000E+35 * exp( -7 * lgt - 11560000 / rt );
	kTroeInf = 1.8130000000E+10;
	fcTroe = 0.38 * exp( -temp / 73 ) + 0.62 * exp( -temp / 1180 );
	k[r63] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, -1.0 );
	k[r64] = 1.0000000000E+13 * exp( -134000000 / rt );
	k[r65] = 3.6000000000E+09;
	k[r66] = 5.0000000000E+10 * exp( -105000000 / rt );
	k[r67] = 1.8000000000E+10;
	k[r68f] = 4.0000000000E+07 * exp( -8900000 / rt );
	k[r71f] = 1.3000000000E+01 * exp( 3 * lgt - 33600000 / rt );
	k[r73f] = 1.6000000000E+04 * exp( 1.83 * lgt - 11600000 / rt );
	k[r76] = 3.0000000000E+10 * exp( 1700000 / rt );
	k[rm69] = 5.0000000000E+10 * exp( -105000000 / rt );
	k[rm70] = 3.0000000000E+10;
	k[r78] = 3.0000000000E+09;
	k[r79f] = 1.5000000000E+11;
	k[r80] = 9.6000000000E+10;
	k[r84f] = 6.0220000000E+10 * exp( -116400000 / rt );
	k[r85f] = 1.2000000000E+11 * exp( -54200000 / rt );
	k[r86] = 1.7200000000E+01 * exp( 2.8 * lgt - 2100000 / rt );
	k[r87] = 1.7200000000E+01 * exp( 2.8 * lgt - 2100000 / rt );
	k[r88] = 1.0000000000E+13 * exp( -248000000 / rt );
	k[r89f] = 4.2000000000E+10 * exp( -16100000 / rt );
	k[r90f] = 2.3000000000E+09 * exp( -2900000 / rt );
	k[r91f] = 1.0000000000E+10;
	kTroe0 = 1.1870000000E+39 * exp( -7.5 * lgt - 190400000 / rt );
	kTroeInf = 2.0000000000E+14 * exp( -166290000 / rt );
	fcTroe = 0.35 * exp( -0 / temp );
	k[r92f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, -1.0 );
	k[r125] = 1.2100000000E+08;
	k[ra125] = 5.4200000000E+09;
	k[r113f] = 2.5000000000E+14 * exp( -319800000 / rt );
	k[r115f] = 5.4000000000E+11 * exp( -62900000 / rt );
	k[r116f] = 6.5000000000E+10 * exp( -24900000 / rt );
	k[r117f] = 4.1500000000E+09 * exp( -46500000 / rt );
	kTroe0 = 2.4610000000E+13 * exp( -126000000 / rt );
	kTroeInf = 3.2000000000E+13 * exp( -167000000 / rt );
	fcTroe = 1 * exp( -temp / 422.8 ) + 0.411 * exp( -73.4 / temp );
	k[r119f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, -1.0 );
	k[r121] = 1.4000000000E+09;
	k[r122] = 1.1000000000E+07 * exp( 6300000 / rt );
	k[r126f] = 1.4000000000E+06 * exp( 1.5 * lgt - 31100000 / rt );
	k[r128] = 7.2000000000E+03 * exp( 2 * lgt - 3600000 / rt );
	k[r129f] = 1.5000000000E-10 * exp( 6 * lgt - 25400000 / rt );
	k[r141f] = 1.8000000000E+11;
	k[r142f] = 1.2000000000E+10 * exp( -27600000 / rt );
	k[r143] = 1.0000000000E+08 * exp( -12470000 / rt );
	k[r144] = 3.0000000000E+07 * exp( -15440000 / rt );
	k[r145] = 3.8000000000E+10;
	k[r147f] = 1.0000000000E+14 * exp( -293000000 / rt );
	k[r148] = 1.0000000000E+09 * exp( -6280000 / rt );
	k[r149] = 2.0000000000E+09 * exp( -32200000 / rt );
	k[r153] = 1.0000000000E+09;
	k[r155f] = 2.0000000000E+13 * exp( -259600000 / rt );
	k[r157f] = 2.0000000000E+10 * exp( -10000000 / rt );
	k[r158f] = 1.0000000000E+14 * exp( -293000000 / rt );
	k[r159] = 1.0000000000E+09 * exp( -6280000 / rt );
	k[r164] = 1.0000000000E+09;
	k[r166] = 3.9800000000E+13 * exp( -293100000 / rt );
	k[r167] = 1.0000000000E+10;
	k[r168] = 1.0000000000E+09;
	k[r193] = 6.0000000000E+08 * exp( -41900000 / rt );
	k[r169f] = 3.1500000000E+15 * exp( -359300000 / rt );
	k[r170f] = 5.0000000000E+09 * exp( -6300000 / rt );
	k[r172f] = 8.9600000000E+09 * exp( -35600000 / rt );
	k[r176] = 7.9000000000E+09;
	k[r178f] = 9.6000000000E+13 * exp( -129800000 / rt );
	k[r179f] = 1.2500000000E+14 * exp( -154900000 / rt );
	k[r180f] = 1.0000000000E+09 * exp( -20900000 / rt );
	k[r207f] = 3.5000000000E+10;
	k[r208] = 7.5000000000E+09;
	k[r211f] = 2.0000000000E+10 * exp( -225900000 / rt );
	k[r213f] = 1.0000000000E+13 * exp( -250000000 / rt );
	k[r215f] = 2.0000000000E+10;
	k[r221f] = 2.0000000000E+09 * exp( -21000000 / rt );
	k[r224f] = 1.5000000000E+11 * exp( -42700000 / rt );
	k[r225f] = 1.5000000000E+11 * exp( -42700000 / rt );
	k[r226f] = 7.0000000000E+10 * exp( -12600000 / rt );
	k[r227] = 7.0000000000E+10 * exp( -12600000 / rt );
	kTroe0 = 5.6010000000E+27 * exp( -5.84 * lgt - 19680000 / rt );
	kTroeInf = 6.7000000000E+04 * exp( 1.43 * lgt - 15500000 / rt );
	fcTroe = 1 * exp( -0 / temp );
	k[r234f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, -1.0 );
	k[r235f] = 5.5000000000E+09 * exp( -10000000 / rt );
	k[r236f] = 5.5000000000E+09 * exp( -10000000 / rt );
	k[r237] = 2.0000000000E+10;
	k[r238] = 2.0000000000E+10;
	k[r239f] = 1.0000000000E+11;
	k[r240] = 5.1000000000E+18 * exp( -3.24 * lgt - 23700000 / rt );
	k[r242] = 5.1000000000E+18 * exp( -3.24 * lgt - 23700000 / rt );
	k[r243] = 4.4800000000E+23 * exp( -4.55 * lgt - 22900000 / rt );
	k[r245f] = 1.0000000000E+09;
	k[r248f] = 3.0000000000E+04 * exp( 2 * lgt - 25100000 / rt );
	k[r250] = 2.0000000000E+04 * exp( 2 * lgt - 8400000 / rt );
	k[r241f] = 1.2000000000E+14 * exp( -206400000 / rt );
	k[rm242] = 1.0000000000E+11 * exp( -154900000 / rt );
	k[r244] = 1.0000000000E+08;
	k[r252] = 4.0000000000E+10;
	k[r255f] = 8.0000000000E+16 * exp( -307400000 / rt );
	k[r258] = 5.0000000000E+10 * exp( -16300000 / rt );
	k[r265] = 6.5000000000E+09;
	k[r269] = 1.0000000000E+08 * exp( -71400000 / rt );
	k[r279] = 2.2650000000E+11 * exp( -12800000 / rt );
	k[r292] = 2.5000000000E+13 * exp( -120600000 / rt );
	k[r294] = 1.0000000000E+09 * exp( -8400000 / rt );
	k[r333] = 1.2990000000E+13 * exp( -120000000 / rt );
	k[r334] = 1.2990000000E+13 * exp( -120000000 / rt );
	k[r335] = 1.2990000000E+13 * exp( -120000000 / rt );
	k[r336] = 1.2990000000E+13 * exp( -120000000 / rt );
	k[r347] = 3.1590000000E+16 * exp( -338580000 / rt );
	k[r354] = 3.2500000000E+06 * exp( 1 * lgt - 2934000 / rt );
	k[r364] = 8.4000000000E+09 * exp( -71060000 / rt );
	k[r369] = 3.0000000000E+09 * exp( -39710000 / rt );
	k[r374] = 2.2500000000E+04 * exp( 2 * lgt - 20900000 / rt );
	k[r381] = 3.6990000000E+13 * exp( -120000000 / rt );
	k[r383] = 1.2990000000E+13 * exp( -120000000 / rt );
	k[r384] = 1.2990000000E+13 * exp( -120000000 / rt );
	k[r385] = 1.2990000000E+13 * exp( -120000000 / rt );
	k[r386] = 1.2990000000E+13 * exp( -120000000 / rt );
	k[r387] = 1.2990000000E+13 * exp( -120000000 / rt );
	k[r400] = 1.5790000000E+16 * exp( -338580000 / rt );
	k[r401] = 1.5790000000E+16 * exp( -338580000 / rt );
	k[r409] = 3.9000000000E+06 * exp( 1 * lgt - 2934000 / rt );
	k[r421] = 1.0000000000E+10 * exp( -71060000 / rt );
	k[r427] = 3.6000000000E+09 * exp( -39710000 / rt );
	k[r433] = 2.7000000000E+04 * exp( 2 * lgt - 20900000 / rt );
	k[r444] = 3.6990000000E+13 * exp( -120000000 / rt );
	k[r456] = 3.6990000000E+13 * exp( -120000000 / rt );
	k[rc395] = 2.5000000000E+16 * exp( -297700000 / rt );
	k[rc396] = 1.0000000000E+10;
	k[rc424] = 8.0000000000E+03 * exp( 2 * lgt - 20900000 / rt );
	k[rc425] = 1.3000000000E+06 * exp( 1.3 * lgt - 2900000 / rt );
	k[r468] = 3.6990000000E+13 * exp( -120000000 / rt );
	k[r481] = 3.6990000000E+13 * exp( -120000000 / rt );
	k[r482] = 3.6990000000E+13 * exp( -120000000 / rt );
	k[r483] = 3.6990000000E+13 * exp( -120000000 / rt );
	k[r484] = 3.6990000000E+13 * exp( -120000000 / rt );
	k[r485] = 3.6990000000E+13 * exp( -120000000 / rt );
	k[r496] = 3.1590000000E+16 * exp( -338580000 / rt );
	k[r497] = 2.5000000000E+10 * exp( -205200000 / rt );
	k[r498] = 1.6000000000E+11 * exp( -199300000 / rt );
	k[r499] = 1.6000000000E+11 * exp( -199300000 / rt );
	k[r500] = 1.6000000000E+11 * exp( -199300000 / rt );
	k[r501] = 1.6000000000E+11 * exp( -199300000 / rt );
	k[r502] = 8.6100000000E+06 * exp( 1.05 * lgt - 7600000 / rt );
	k[r503] = 1.0400000000E+07 * exp( 1.25 * lgt - 2900000 / rt );
	k[r504] = 1.0400000000E+07 * exp( 1.25 * lgt - 2900000 / rt );
	k[r505] = 1.0400000000E+07 * exp( 1.25 * lgt - 2900000 / rt );
	k[r506] = 1.0400000000E+07 * exp( 1.25 * lgt - 2900000 / rt );
	k[r512] = 1.1200000000E+10 * exp( -81200000 / rt );
	k[r513] = 2.6900000000E+10 * exp( -71200000 / rt );
	k[r514] = 2.6900000000E+10 * exp( -71200000 / rt );
	k[r515] = 2.6900000000E+10 * exp( -71200000 / rt );
	k[r516] = 2.6900000000E+10 * exp( -71200000 / rt );
	k[r522] = 5.6000000000E+04 * exp( 2 * lgt - 32200000 / rt );
	k[r523] = 7.2000000000E+04 * exp( 2 * lgt - 20900000 / rt );
	k[r524] = 7.2000000000E+04 * exp( 2 * lgt - 20900000 / rt );
	k[r525] = 7.2000000000E+04 * exp( 2 * lgt - 20900000 / rt );
	k[r526] = 7.2000000000E+04 * exp( 2 * lgt - 20900000 / rt );
	k[rl10f] = 2.0000000000E+09;
	k[rl10b] = 5.0000000000E+14 * exp( -117100000 / rt );
	k[rl20f] = 2.0000000000E+09;
	k[rl20b] = 5.0000000000E+14 * exp( -117100000 / rt );
	k[rl30f] = 2.0000000000E+09;
	k[rl30b] = 5.0000000000E+14 * exp( -117100000 / rt );
	k[rl40f] = 2.0000000000E+09;
	k[rl40b] = 5.0000000000E+14 * exp( -117100000 / rt );
	k[rl50f] = 2.0000000000E+09;
	k[rl50b] = 5.0000000000E+14 * exp( -117100000 / rt );
	k[rl60f] = 6.4000000000E+12 * exp( -61100000 / rt );
	k[rl60b] = 8.0000000000E+11 * exp( -27600000 / rt );
	k[rl70] = 3.0000000000E+08;
	k[rl80] = 3.1620000000E+12 * exp( -82840000 / rt );
	k[rl90] = 1.1000000000E+09 * exp( -31400000 / rt );
	k[rl100] = 1.0000000000E+15 * exp( -170000000 / rt );
	k[rl230] = 2.0000000000E+13 * exp( -62800000 / rt );
	kTroe0 = 8.0010000000E+17 * exp( -2.66 * lgt - 14200000 / rt );
	kTroeInf = 2.9000000000E+04 * exp( 1.4 * lgt - 9600000 / rt );
	fcTroe = 1 * exp( -0 / temp );
	k[rs257f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, -1.0 );
	k[rs260f] = 2.0000000000E+10;
	k[rs261f] = 1.0000000000E+08;
	k[rs266f] = 5.5000000000E+09 * exp( -10000000 / rt );
	kTroe0 = 4.4990000000E+20 * exp( -3.28 * lgt - 42700000 / rt );
	kTroeInf = 3.7000000000E+04 * exp( 1.33 * lgt - 8800000 / rt );
	fcTroe = 1 * exp( -0 / temp );
	k[rs268f] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf
				, fcTroe, -1.0 );
	k[rs271] = 5.0000000000E+08;
	k[rs272] = 2.0000000000E+08 * exp( -8370000 / rt );
	k[rs273] = 1.0000000000E+10;
	k[rs274f] = 3.0000000000E+13 * exp( -305000000 / rt );
	k[rs276] = 4.5000000000E+15 * exp( -449500000 / rt );
	k[rsA276] = 4.5000000000E+13 * exp( -303430000 / rt );
	k[rs277f] = 4.5000000000E+15 * exp( -449500000 / rt );
	k[rs278f] = 2.5000000000E+11 * exp( -66940000 / rt );
	k[rs279f] = 1.6300000000E+05 * exp( 1.42 * lgt - 6070000 / rt );
	k[rs280f] = 3.5530000000E-02 * exp( 3.8 * lgt - 3930000 / rt );
	k[rsA280f] = 6.3000000000E+10 * exp( -251040000 / rt );
	k[rs281f] = 2.2300000000E+10 * exp( -20000000 / rt );
	k[rs2A281f] = 5.0000000000E+10 * exp( -4180000 / rt );
	k[rs275] = 8.0900000000E+08 * exp( -31250000 / rt );
	k[rs3A281f] = 2.5300000000E+11;
	k[rs284f] = 1.1500000000E+11 * exp( -51880000 / rt );
	k[rsA284f] = 2.2300000000E+10 * exp( -33170000 / rt );
	k[rs285f] = 6.0200000000E+09;
	k[rsA285f] = 1.2600000000E+10 * exp( -12130000 / rt );
	k[rsA33f] = 2.5100000000E+11 * exp( -183700000 / rt );
	k[rsA39f] = 1.0000000000E+11;
	k[rsA40f] = 1.0000000000E+11;
	k[rsA41f] = 1.0000000000E+10;
	k[rsA43f] = 3.0000000000E+10;
	k[rsA44f] = 2.1900000000E+05 * exp( 1.77 * lgt - 12550000 / rt );
	k[rsA45f] = 1.8100000000E+10 * exp( -12890000 / rt );
	k[rsA46f] = 3.4300000000E+06 * exp( 1.18 * lgt + 1870000 / rt );
	k[rsA49f] = 1.0000000000E+10 * exp( -86660000 / rt );
	k[rsA50f] = 2.5100000000E+11 * exp( -183680000 / rt );
	k[rs288f] = 2.5000000000E+11 * exp( -0.48 * lgt - 6276000 / rt );
	k[rs290] = 1.0000000000E+10;
	k[rs297f] = 4.2000000000E+08 * exp( 0.48 * lgt + 300000 / rt );
	k[rs298f] = 5.0000000000E+10 * exp( -0.48 * lgt - 6276000 / rt );
	k[rs299f] = 2.1000000000E+10 * exp( -19100000 / rt );
	k[rs300f] = 2.0000000000E+10;
	k[rs303f] = 1.4000000000E+29 * exp( -6.28 * lgt - 14800000 / rt );
	k[rAs303f] = 1.1000000000E+59 * exp( -14.56 * lgt - 33100000 / rt );
	k[rs304f] = 1.0000000000E+10;
	k[rs306f] = 4.5000000000E+15 * exp( -449500000 / rt );
	k[rs307f] = 2.5000000000E+11 * exp( -66940000 / rt );
	k[rs308f] = 1.6300000000E+05 * exp( 1.42 * lgt - 6070000 / rt );
	k[rs310f] = 2.0000000000E+10;
	k[rs313f] = 1.0000000000E+10;
	k[rs315f] = 1.0000000000E+11;
	k[rs318f] = 2.0000000000E+10;
	k[rs319f] = 1.0000000000E+11;
	k[rs322f] = 1.0000000000E+11;
	k[rs325f] = 2.0000000000E+10;
	k[rs326f] = 1.0000000000E+10;
	k[rs327f] = 1.0000000000E+10;
	k[rs329f] = 1.0000000000E+11;
	k[rs332f] = 2.0000000000E+10;
	k[rs333f] = 1.0000000000E+10;
	k[rs334] = 1.0000000000E+10;
	k[rs335] = 1.0000000000E+10;
	k[rs339] = 1.0000000000E+10;
	k[rs342f] = 1.0000000000E+10;
	k[rs342b] = 2.5090000000E+14 * exp( -220790000 / rt );
	k[rF1f] = 3.0000000000E+11 * exp( -173000000 / rt );
	k[rAF1f] = 4.0000000000E+08 * exp( -58900000 / rt );
	k[rAF2f] = 1.2000000000E+10 * exp( -5379000 / rt );
	k[rF3f] = 5.1900000000E+05 * exp( 1 * lgt - 3658000 / rt );
	k[rF4f] = 1.2600000000E+11 * exp( -35150000 / rt );
	k[rF7f] = 1.6500000000E+11;
	k[rF9] = 2.5000000000E+11;
	k[rAF12f] = 6.0600000000E+10;
	k[rAF13f] = 2.5100000000E+08 * exp( 0.4 * lgt );
	k[rAF13b] = 1.3800000000E+18 * exp( -0.5 * lgt - 281000000 / rt );


	w[r1f] = k[r1f] * c[sO2] * c[sH];
	w[r2f] = k[r2f] * c[sH2] * c[sO];
	w[r3f] = k[r3f] * c[sH2] * c[sOH];
	w[r4f] = k[r4f] * c[sOH] * c[sOH];
	w[r5f] = k[r5f] * c[sH] * c[sO2] * M[mM1];
	w[r6f] = k[r6f] * c[sHO2] * c[sH];
	w[r7f] = k[r7f] * c[sHO2] * c[sH];
	w[r8f] = k[r8f] * c[sHO2] * c[sOH];
	w[r9f] = k[r9f] * c[sHO2] * c[sH];
	w[r10f] = k[r10f] * c[sHO2] * c[sO];
	w[r11f] = k[r11f] * c[sHO2] * c[sHO2];
	w[r12f] = k[r12f] * c[sOH] * c[sOH] * M[mM1];
	w[r13f] = k[r13f] * c[sH2O2] * c[sH];
	w[r14f] = k[r14f] * c[sH2O2] * c[sH];
	w[r15f] = k[r15f] * c[sH2O2] * c[sO];
	w[r16f] = k[r16f] * c[sH2O2] * c[sOH];
	w[r17f] = k[r17f] * c[sH] * c[sH] * M[mM1];
	w[r18f] = k[r18f] * c[sH] * c[sOH] * M[mM1];
	w[r19f] = k[r19f] * c[sO] * c[sO] * M[mM1];
	w[r20f] = k[r20f] * c[sCO] * c[sOH];
	w[r21] = k[r21] * c[sCO] * c[sHO2];
	w[r22] = k[r22] * c[sCO] * c[sO] * M[mM1];
	w[r24] = k[r24] * c[sCH] * c[sO];
	w[r25] = k[r25] * c[sCH] * c[sO2];
	w[r26] = k[r26] * c[sCH] * c[sCO2];
	w[r27f] = k[r27f] * c[sCH] * c[sH2O];
	w[r28f] = k[r28f] * c[sHCO] * M[mM1];
	w[r29] = k[r29] * c[sHCO] * c[sH];
	w[r32] = k[r32] * c[sHCO] * c[sOH];
	w[r35] = k[r35] * c[sHCO] * c[sO2];
	w[r37f] = k[r37f] * c[s3XCH2] * c[sH];
	w[r39] = k[r39] * c[s3XCH2] * c[sO2];
	w[r40] = k[r40] * c[s3XCH2] * c[sO2];
	w[r42] = k[r42] * c[s3XCH2] * c[s3XCH2];
	w[r43] = k[r43] * c[s3XCH2] * c[s3XCH2];
	w[r44] = k[r44] * c[s3XCH2] * c[sCH3];
	w[r45f] = k[r45f] * c[s1XCH2] * M[mM1];
	w[r46f] = k[r46f] * c[s1XCH2] * c[sH2];
	w[r48] = k[r48] * c[sCH2O] * M[mM1];
	w[r49] = k[r49] * c[sCH2O] * c[sH];
	w[r50] = k[r50] * c[sCH2O] * c[sO];
	w[r51] = k[r51] * c[sCH2O] * c[sOH];
	w[r54f] = k[r54f] * c[sCH2O] * c[sO2];
	w[r56f] = k[r56f] * c[sCH3] * c[sH];
	w[r57] = k[r57] * c[sCH3] * c[sO];
	w[r58] = k[r58] * c[sCH3] * c[sO2];
	w[r59] = k[r59] * c[sCH3] * c[sOH];
	w[r60] = k[r60] * c[sCH3O] * c[sH];
	w[r61] = k[r61] * c[sCH3] * c[sHO2];
	w[r62f] = k[r62f] * c[sCH3] * c[sH];
	w[r63] = k[r63] * c[sCH3] * c[sCH3];
	w[r64] = k[r64] * c[sCH3] * c[sCH3];
	w[r65] = k[r65] * c[sCH3] * c[sHO2];
	w[r66] = k[r66] * c[sCH3O] * M[mM1];
	w[r67] = k[r67] * c[sCH3O] * c[sH];
	w[r68f] = k[r68f] * c[sCH3O] * c[sO2];
	w[r71f] = k[r71f] * c[sCH4] * c[sH];
	w[r73f] = k[r73f] * c[sCH4] * c[sOH];
	w[r76] = k[r76] * c[sCH4] * c[sCH];
	w[rm69] = k[rm69] * c[sCH2OH] * M[mM1];
	w[rm70] = k[rm70] * c[sCH2OH] * c[sH];
	w[r78] = k[r78] * c[sC2H] * c[sO2];
	w[r79f] = k[r79f] * c[sHCCO] * c[sH];
	w[r80] = k[r80] * c[sHCCO] * c[sO];
	w[r84f] = k[r84f] * c[sC2H2] * c[sH];
	w[r85f] = k[r85f] * c[sC2H2] * c[sOH];
	w[r86] = k[r86] * c[sC2H2] * c[sO];
	w[r87] = k[r87] * c[sC2H2] * c[sO];
	w[r88] = k[r88] * c[sCH2CO] * M[mM1];
	w[r89f] = k[r89f] * c[sCH2CO] * c[sH];
	w[r90f] = k[r90f] * c[sCH2CO] * c[sO];
	w[r91f] = k[r91f] * c[sCH2CO] * c[sOH];
	w[r92f] = k[r92f] * c[sC2H3];
	w[r125] = k[r125] * c[sC2H3] * c[sO2];
	w[ra125] = k[ra125] * c[sC2H3] * c[sO2];
	w[r113f] = k[r113f] * c[sC2H4] * M[mM1];
	w[r115f] = k[r115f] * c[sC2H4] * c[sH];
	w[r116f] = k[r116f] * c[sC2H4] * c[sOH];
	w[r117f] = k[r117f] * c[sC2H4] * c[sCH3];
	w[r119f] = k[r119f] * c[sC2H5];
	w[r121] = k[r121] * c[sC2H5] * c[sC2H5];
	w[r122] = k[r122] * c[sC2H5] * c[sO2];
	w[r126f] = k[r126f] * c[sC2H6] * c[sH];
	w[r128] = k[r128] * c[sC2H6] * c[sOH];
	w[r129f] = k[r129f] * c[sC2H6] * c[sCH3];
	w[r141f] = k[r141f] * c[sC2H2] * c[s1XCH2];
	w[r142f] = k[r142f] * c[sC2H2] * c[s3XCH2];
	w[r143] = k[r143] * c[sC2H2] * c[sHCCO];
	w[r144] = k[r144] * c[sC3H3] * c[sO2];
	w[r145] = k[r145] * c[sC3H3] * c[sO];
	w[r147f] = k[r147f] * c[sC3H4] * M[mM1];
	w[r148] = k[r148] * c[sC3H4] * c[sH];
	w[r149] = k[r149] * c[sC3H4] * c[sCH3];
	w[r153] = k[r153] * c[sC3H4] * c[sOH];
	w[r155f] = k[r155f] * c[sC3H4];
	w[r157f] = k[r157f] * c[sPXC3H4] * c[sH];
	w[r158f] = k[r158f] * c[sPXC3H4] * M[mM1];
	w[r159] = k[r159] * c[sPXC3H4] * c[sH];
	w[r164] = k[r164] * c[sPXC3H4] * c[sOH];
	w[r166] = k[r166] * c[sC3H5];
	w[r167] = k[r167] * c[sC3H5] * c[sH];
	w[r168] = k[r168] * c[sC3H5] * c[sCH3];
	w[r193] = k[r193] * c[sC3H5] * c[sO2];
	w[r169f] = k[r169f] * c[sC3H6];
	w[r170f] = k[r170f] * c[sC3H6] * c[sH];
	w[r172f] = k[r172f] * c[sC3H6] * c[sCH3];
	w[r176] = k[r176] * c[sC3H6] * c[sOH];
	w[r178f] = k[r178f] * c[sNXC3H7];
	w[r179f] = k[r179f] * c[sNXC3H7];
	w[r180f] = k[r180f] * c[sNXC3H7] * c[sO2];
	w[r207f] = k[r207f] * c[sC2H2] * c[sC2H];
	w[r208] = k[r208] * c[sC4H2] * c[sOH];
	w[r211f] = k[r211f] * c[sC2H2] * c[sC2H2];
	w[r213f] = k[r213f] * c[sSXC4H3] * M[mM1];
	w[r215f] = k[r215f] * c[sSXC4H3] * c[sH];
	w[r221f] = k[r221f] * c[sC2H2] * c[sC2H3];
	w[r224f] = k[r224f] * c[sC4H4] * c[sH];
	w[r225f] = k[r225f] * c[sC4H4] * c[sH];
	w[r226f] = k[r226f] * c[sC4H4] * c[sOH];
	w[r227] = k[r227] * c[sC4H4] * c[sOH];
	w[r234f] = k[r234f] * c[sC2H2] * c[sC2H3];
	w[r235f] = k[r235f] * c[sC4H4] * c[sH];
	w[r236f] = k[r236f] * c[sC4H4] * c[sH];
	w[r237] = k[r237] * c[sSXC4H5] * c[sH];
	w[r238] = k[r238] * c[sUXC4H5] * c[sH];
	w[r239f] = k[r239f] * c[sUXC4H5] * c[sH];
	w[r240] = k[r240] * c[sUXC4H5] * c[sO2];
	w[r242] = k[r242] * c[sSXC4H5] * c[sO2];
	w[r243] = k[r243] * c[sSXC4H5] * c[sO2];
	w[r245f] = k[r245f] * c[sC3H3] * c[sCH3];
	w[r248f] = k[r248f] * c[sC4H6] * c[sH];
	w[r250] = k[r250] * c[sC4H6] * c[sOH];
	w[r241f] = k[r241f] * c[sC4H7];
	w[rm242] = k[rm242] * c[sC4H7];
	w[r244] = k[r244] * c[sC4H7] * c[sO2];
	w[r252] = k[r252] * c[sC4H7] * c[sC3H5];
	w[r255f] = k[r255f] * c[s1XC4H8];
	w[r258] = k[r258] * c[s1XC4H8] * c[sH];
	w[r265] = k[r265] * c[s1XC4H8] * c[sOH];
	w[r269] = k[r269] * c[s1XC4H8] * c[sHO2];
	w[r279] = k[r279] * c[s1XC4H8] * c[sOH];
	w[r292] = k[r292] * c[sPXC4H9];
	w[r294] = k[r294] * c[sPXC4H9] * c[sO2];
	w[r333] = k[r333] * c[s1XC5H9];
	w[r334] = k[r334] * c[s1XC5H9];
	w[r335] = k[r335] * c[s1XC5H9];
	w[r336] = k[r336] * c[s1XC5H9];
	w[r347] = k[r347] * c[s1XC5H10];
	w[r354] = k[r354] * c[s1XC5H10] * c[sOH];
	w[r364] = k[r364] * c[s1XC5H10] * c[sHO2];
	w[r369] = k[r369] * c[s1XC5H10] * c[sCH3];
	w[r374] = k[r374] * c[s1XC5H10] * c[sH];
	w[r381] = k[r381] * c[s1XC5H11];
	w[r383] = k[r383] * c[sC6H11];
	w[r384] = k[r384] * c[sC6H11];
	w[r385] = k[r385] * c[sC6H11];
	w[r386] = k[r386] * c[sC6H11];
	w[r387] = k[r387] * c[sC6H11];
	w[r400] = k[r400] * c[s1XC6H12];
	w[r401] = k[r401] * c[s1XC6H12];
	w[r409] = k[r409] * c[s1XC6H12] * c[sOH];
	w[r421] = k[r421] * c[s1XC6H12] * c[sHO2];
	w[r427] = k[r427] * c[s1XC6H12] * c[sCH3];
	w[r433] = k[r433] * c[s1XC6H12] * c[sH];
	w[r444] = k[r444] * c[s1XC6H13];
	w[r456] = k[r456] * c[s1XC7H15];
	w[rc395] = k[rc395] * c[s1XC7H14];
	w[rc396] = k[rc396] * c[sPXC4H9] * c[sC3H5];
	w[rc424] = k[rc424] * c[s1XC7H14] * c[sH];
	w[rc425] = k[rc425] * c[s1XC7H14] * c[sOH];
	w[r468] = k[r468] * c[s1XC8H17];
	w[r481] = k[r481] * c[s2XC10H21];
	w[r482] = k[r482] * c[s3XC10H21];
	w[r483] = k[r483] * c[s4XC10H21];
	w[r484] = k[r484] * c[s5XC10H21];
	w[r485] = k[r485] * c[s1XC10H21];
	w[r496] = k[r496] * c[sNXC10H22];
	w[r497] = k[r497] * c[sNXC10H22] * c[sO2];
	w[r498] = k[r498] * c[sNXC10H22] * c[sO2];
	w[r499] = k[r499] * c[sNXC10H22] * c[sO2];
	w[r500] = k[r500] * c[sNXC10H22] * c[sO2];
	w[r501] = k[r501] * c[sNXC10H22] * c[sO2];
	w[r502] = k[r502] * c[sNXC10H22] * c[sOH];
	w[r503] = k[r503] * c[sNXC10H22] * c[sOH];
	w[r504] = k[r504] * c[sNXC10H22] * c[sOH];
	w[r505] = k[r505] * c[sNXC10H22] * c[sOH];
	w[r506] = k[r506] * c[sNXC10H22] * c[sOH];
	w[r512] = k[r512] * c[sNXC10H22] * c[sHO2];
	w[r513] = k[r513] * c[sNXC10H22] * c[sHO2];
	w[r514] = k[r514] * c[sNXC10H22] * c[sHO2];
	w[r515] = k[r515] * c[sNXC10H22] * c[sHO2];
	w[r516] = k[r516] * c[sNXC10H22] * c[sHO2];
	w[r522] = k[r522] * c[sNXC10H22] * c[sH];
	w[r523] = k[r523] * c[sNXC10H22] * c[sH];
	w[r524] = k[r524] * c[sNXC10H22] * c[sH];
	w[r525] = k[r525] * c[sNXC10H22] * c[sH];
	w[r526] = k[r526] * c[sNXC10H22] * c[sH];
	w[rl10f] = k[rl10f] * c[s1XC10H21] * c[sO2];
	w[rl10b] = k[rl10b] * c[s2XRO2XC10H21O2];
	w[rl20f] = k[rl20f] * c[s2XC10H21] * c[sO2];
	w[rl20b] = k[rl20b] * c[s2XRO2XC10H21O2];
	w[rl30f] = k[rl30f] * c[s3XC10H21] * c[sO2];
	w[rl30b] = k[rl30b] * c[s2XRO2XC10H21O2];
	w[rl40f] = k[rl40f] * c[s4XC10H21] * c[sO2];
	w[rl40b] = k[rl40b] * c[s2XRO2XC10H21O2];
	w[rl50f] = k[rl50f] * c[s5XC10H21] * c[sO2];
	w[rl50b] = k[rl50b] * c[s2XRO2XC10H21O2];
	w[rl60f] = k[rl60f] * c[s2XRO2XC10H21O2];
	w[rl60b] = k[rl60b] * c[s3X2XRO2HXC10H20OOH];
	w[rl70] = k[rl70] * c[s3X2XRO2HXC10H20OOH] * c[sO2];
	w[rl80] = k[rl80] * c[s32O2RO2HXO2C10H20OOH];
	w[rl90] = k[rl90] * c[s3X2HOROHXHO2C10H19OOH];
	w[rl100] = k[rl100] * c[s3X2ORO2HXOC10H19OOH];
	w[rl230] = k[rl230] * c[s2X1OROXOC10H19O];
	w[rs257f] = k[rs257f] * c[sUXC4H3] * c[sC2H2];
	w[rs260f] = k[rs260f] * c[sC3H3] * c[sC3H3];
	w[rs261f] = k[rs261f] * c[sC3H4] * c[sC3H3];
	w[rs266f] = k[rs266f] * c[sC6H6] * c[sH];
	w[rs268f] = k[rs268f] * c[sUXC4H5] * c[sC2H2];
	w[rs271] = k[rs271] * c[sC3H3] * c[sC3H3];
	w[rs272] = k[rs272] * c[sC3H4] * c[sC3H3];
	w[rs273] = k[rs273] * c[sUXC6H7];
	w[rs274f] = k[rs274f] * c[sA1XXC6H5];
	w[rs276] = k[rs276] * c[sA1XC6H6];
	w[rsA276] = k[rsA276] * c[sA1XXC6H5];
	w[rs277f] = k[rs277f] * c[sA1XC6H6];
	w[rs278f] = k[rs278f] * c[sA1XC6H6] * c[sH];
	w[rs279f] = k[rs279f] * c[sA1XC6H6] * c[sOH];
	w[rs280f] = k[rs280f] * c[sA1XC6H6] * c[sO];
	w[rsA280f] = k[rsA280f] * c[sA1XC6H6] * c[sO2];
	w[rs281f] = k[rs281f] * c[sA1XC6H6] * c[sO];
	w[rs2A281f] = k[rs2A281f] * c[sA1XXC6H5] * c[sHO2];
	w[rs275] = k[rs275] * c[sA1XXC6H5] * c[sO2];
	w[rs3A281f] = k[rs3A281f] * c[sC6H5O] * c[sH];
	w[rs284f] = k[rs284f] * c[sC6H5OH] * c[sH];
	w[rsA284f] = k[rsA284f] * c[sC6H5OH] * c[sH];
	w[rs285f] = k[rs285f] * c[sC6H5OH] * c[sOH];
	w[rsA285f] = k[rsA285f] * c[sC6H5OH] * c[sO];
	w[rsA33f] = k[rsA33f] * c[sC6H5O];
	w[rsA39f] = k[rsA39f] * c[sC5H5] * c[sH];
	w[rsA40f] = k[rsA40f] * c[sC5H5] * c[sO];
	w[rsA41f] = k[rsA41f] * c[sC5H5] * c[sO];
	w[rsA43f] = k[rsA43f] * c[sC5H5] * c[sHO2];
	w[rsA44f] = k[rsA44f] * c[sC5H6] * c[sH];
	w[rsA45f] = k[rsA45f] * c[sC5H6] * c[sO];
	w[rsA46f] = k[rsA46f] * c[sC5H6] * c[sOH];
	w[rsA49f] = k[rsA49f] * c[sC5H6] * c[sO2];
	w[rsA50f] = k[rsA50f] * c[sC5H5O];
	w[rs288f] = k[rs288f] * c[sA1C2HXC8H6] * c[sH];
	w[rs290] = k[rs290] * c[sA1C2HXC8H6] * c[sOH];
	w[rs297f] = k[rs297f] * c[sA1C2HYXC8H5] * c[sH];
	w[rs298f] = k[rs298f] * c[sA1C2HXC8H6] * c[sH];
	w[rs299f] = k[rs299f] * c[sA1C2HXC8H6] * c[sOH];
	w[rs300f] = k[rs300f] * c[sA1C2HXC8H6] * c[sC2H];
	w[rs303f] = k[rs303f] * c[sA1C2HYXC8H5] * c[sC2H2];
	w[rAs303f] = k[rAs303f] * c[sA1C2HYXC8H5] * c[sC2H2];
	w[rs304f] = k[rs304f] * c[sA1C2HACXC10H7];
	w[rs306f] = k[rs306f] * c[sA2XC10H8];
	w[rs307f] = k[rs307f] * c[sA2XC10H8] * c[sH];
	w[rs308f] = k[rs308f] * c[sA2XC10H8] * c[sOH];
	w[rs310f] = k[rs310f] * c[sA2XC10H8] * c[sCH3];
	w[rs313f] = k[rs313f] * c[sA2XXXC10H7] * c[sC2H2];
	w[rs315f] = k[rs315f] * c[sA2R5XC12H8] * c[sH];
	w[rs318f] = k[rs318f] * c[sA2R5XC12H8] * c[sCH3];
	w[rs319f] = k[rs319f] * c[sA2R5C2HXC14H8] * c[sH];
	w[rs322f] = k[rs322f] * c[sA2R5C2HXC14H8] * c[sH];
	w[rs325f] = k[rs325f] * c[sA2R5C2HXC14H8] * c[sCH3];
	w[rs326f] = k[rs326f] * c[sA2R5C2HYXC14H7] * c[sC2H2];
	w[rs327f] = k[rs327f] * c[sANC2HACXC16H9];
	w[rs329f] = k[rs329f] * c[sA3R5XC16H10] * c[sH];
	w[rs332f] = k[rs332f] * c[sA3R5XC16H10] * c[sCH3];
	w[rs333f] = k[rs333f] * c[sA3R5XXC16H9] * c[sC2H2];
	w[rs334] = k[rs334] * c[sA2R5XC12H8] * c[sOH];
	w[rs335] = k[rs335] * c[sA2R5XXC12H7] * c[sO2];
	w[rs339] = k[rs339] * c[sA3R5XC16H10] * c[sOH];
	w[rs342f] = k[rs342f] * c[sA3R5ACXC18H11];
	w[rs342b] = k[rs342b] * c[sA4XC18H10] * c[sH];
	w[rF1f] = k[rF1f] * c[sA2CH3XC11H10] * c[sO2];
	w[rAF1f] = k[rAF1f] * c[sA2CH3XC11H10] * c[sHO2];
	w[rAF2f] = k[rAF2f] * c[sA2CH3XC11H10] * c[sO];
	w[rF3f] = k[rF3f] * c[sA2CH3XC11H10] * c[sOH];
	w[rF4f] = k[rF4f] * c[sA2CH3XC11H10] * c[sH];
	w[rF7f] = k[rF7f] * c[sA2CH2XC11H9] * c[sO];
	w[rF9] = k[rF9] * c[sA2CH2XC11H9] * c[sHO2];
	w[rAF12f] = k[rAF12f] * c[sA2CH2XC11H9] * c[sOH];
	w[rAF13f] = k[rAF13f] * c[sA2CH2XC11H9] * c[sA2CH2XC11H9];
	w[rAF13b] = k[rAF13b] * c[sBIBENZYLXC22H18];


	cdot[sN2] = 0.0;

	cdot[sO2] = - w[r1f] - w[r5f] + w[r7f] + w[r8f]
		+ w[r10f] + w[r11f] + w[r19f] - w[r25]
		- w[r35] - w[r39] - w[r40] - w[r54f]
		- w[r58] + w[r65] - w[r68f] - w[r78]
		- w[r125] - w[ra125] - w[r122] - w[r144]
		- w[r193] - w[r180f] - w[r240] - w[r242]
		- w[r243] - w[r244] - w[r294] - w[r497]
		- w[r498] - w[r499] - w[r500] - w[r501]
		- w[rl10f] + w[rl10b] - w[rl20f] + w[rl20b]
		- w[rl30f] + w[rl30b] - w[rl40f] + w[rl40b]
		- w[rl50f] + w[rl50b] - w[rl70] - w[rsA280f]
		- w[rs275] - w[rsA49f] - w[rs335] - w[rF1f];

	cdot[sH] = - w[r1f] + w[r2f] + w[r3f] - w[r5f]
		- w[r6f] - w[r7f] - w[r9f] - w[r13f]
		- w[r14f] - 2 * w[r17f] - w[r18f] + w[r20f]
		+ w[r24] + w[r28f] - w[r29] - w[r37f]
		+ w[r39] + 2 * w[r43] + w[r44] + w[r46f]
		+ w[r48] - w[r49] - w[r56f] + w[r57]
		+ w[r59] - w[r60] - w[r62f] + w[r66]
		- w[r67] - w[r71f] + w[r76] + w[rm69]
		- w[rm70] - w[r79f] + w[r80] - w[r84f]
		+ w[r87] - w[r89f] + w[r92f] - w[r115f]
		+ w[r119f] - w[r126f] + w[r141f] + w[r142f]
		+ w[r147f] - w[r148] - w[r157f] + w[r158f]
		- w[r159] + w[r166] - w[r167] - w[r170f]
		+ w[r179f] + w[r207f] + w[r211f] + w[r213f]
		- w[r215f] + w[r221f] - w[r224f] - w[r225f]
		- w[r235f] - w[r236f] - w[r237] - w[r238]
		- w[r239f] + w[r239f] - w[r248f] + w[r241f]
		- w[r258] - w[r374] - w[r433] - w[rc424]
		- w[r522] - w[r523] - w[r524] - w[r525]
		- w[r526] + w[rs261f] - w[rs266f] + w[rs272]
		+ w[rs273] + w[rs277f] - w[rs278f] + w[rs281f]
		- w[rs3A281f] - w[rs284f] - w[rsA284f] - w[rsA39f]
		- w[rsA44f] - w[rs288f] - w[rs297f] - w[rs298f]
		+ w[rs306f] - w[rs307f] + w[rs313f] - w[rs315f]
		- w[rs319f] - w[rs322f] - w[rs329f] + w[rs342f]
		- w[rs342b] - w[rF4f] + w[rF7f];

	cdot[sOH] = w[r1f] + w[r2f] - w[r3f] - 2 * w[r4f]
		+ 2 * w[r6f] - w[r8f] + w[r10f] - 2 * w[r12f]
		+ w[r14f] + w[r15f] - w[r16f] - w[r18f]
		- w[r20f] + w[r21] - w[r32] + w[r39]
		+ w[r50] - w[r51] + w[r58] - w[r59]
		+ w[r60] + w[r61] - w[r73f] - w[r85f]
		- w[r91f] - w[r116f] - w[r128] - w[r153]
		- w[r164] - w[r176] - w[r208] - w[r226f]
		- w[r227] - w[r250] - w[r265] - w[r279]
		- w[r354] - w[r409] - w[rc425] - w[r502]
		- w[r503] - w[r504] - w[r505] - w[r506]
		+ w[rl90] + w[rl100] - w[rs279f] + w[rs280f]
		+ w[rs2A281f] + w[rsA284f] - w[rs285f] + w[rsA285f]
		+ w[rsA43f] + w[rsA45f] - w[rsA46f] + w[rsA49f]
		- w[rs290] - w[rs299f] - w[rs308f] - w[rs334]
		- w[rs339] + w[rAF2f] - w[rF3f] - w[rAF12f];

	cdot[sO] = w[r1f] - w[r2f] + w[r4f] + w[r9f]
		- w[r10f] - w[r15f] - 2 * w[r19f] - w[r22]
		- w[r24] + w[r25] - w[r50] - w[r57]
		+ w[r78] - w[r80] - w[r86] - w[r87]
		- w[r90f] - w[r145] - w[rs280f] - w[rs281f]
		+ w[rs275] - w[rsA285f] - w[rsA40f] - w[rsA41f]
		- w[rsA45f] - w[rAF2f] - w[rF7f];

	cdot[sH2] = - w[r2f] - w[r3f] + w[r7f] + w[r13f]
		+ w[r17f] + w[r29] + w[r37f] + w[r40]
		+ w[r42] - w[r46f] + w[r49] + w[r56f]
		+ w[r64] + w[r67] + w[r71f] + w[rm70]
		+ w[r84f] + w[r113f] + w[r115f] + w[r126f]
		+ w[r148] + w[r159] + w[r167] + w[r170f]
		+ w[r215f] + w[r224f] + w[r225f] + w[r237]
		+ w[r238] + w[r248f] + w[r258] + w[r374]
		+ w[r433] + w[rc424] + w[r522] + w[r523]
		+ w[r524] + w[r525] + w[r526] + w[rs278f]
		+ w[rs284f] + w[rsA44f] + w[rs298f] + w[rs307f]
		+ w[rs315f] + w[rs322f] + w[rs329f] + w[rF4f]
		+ w[rAF12f];

	cdot[sH2O] = w[r3f] + w[r4f] + w[r8f] + w[r9f]
		+ w[r14f] + w[r16f] + w[r18f] - w[r27f]
		+ w[r32] + w[r51] + w[r73f] + w[r85f]
		+ w[r116f] + w[r128] + w[r226f] + w[r227]
		+ w[r250] + w[r279] + w[r354] + w[r409]
		+ w[rc425] + w[r502] + w[r503] + w[r504]
		+ w[r505] + w[r506] + w[rs279f] + w[rs285f]
		+ w[rsA46f] + w[rs299f] + w[rs308f] + w[rF3f]
		+ w[rF9];

	cdot[sHO2] = w[r5f] - w[r6f] - w[r7f] - w[r8f]
		- w[r9f] - w[r10f] - 2 * w[r11f] + w[r13f]
		+ w[r15f] + w[r16f] - w[r21] + w[r35]
		+ w[r54f] - w[r61] - w[r65] + w[r68f]
		+ w[r125] + w[r122] + w[r193] + w[r180f]
		+ w[r240] + w[r242] + w[r244] - w[r269]
		+ w[r294] - w[r364] - w[r421] + w[r497]
		+ w[r498] + w[r499] + w[r500] + w[r501]
		- w[r512] - w[r513] - w[r514] - w[r515]
		- w[r516] + w[rsA280f] - w[rs2A281f] - w[rsA43f]
		+ w[rF1f] - w[rAF1f] - w[rF9];

	cdot[sH2O2] = w[r11f] + w[r12f] - w[r13f] - w[r14f]
		- w[r15f] - w[r16f] + w[r269] + w[r364]
		+ w[r421] + w[r512] + w[r513] + w[r514]
		+ w[r515] + w[r516] + w[rAF1f];

	cdot[sCO] = - w[r20f] - w[r21] - w[r22] + w[r24]
		+ w[r26] + w[r28f] + w[r29] + w[r32]
		+ w[r35] + w[r39] + w[r79f] + 2 * w[r80]
		+ w[r86] + w[r88] + w[r89f] + w[r143]
		+ w[rsA33f] + w[rsA40f] + w[rsA50f] + 2 * w[rs335];

	cdot[sCO2] = w[r20f] + w[r21] + w[r22] - w[r26]
		+ w[r40];

	cdot[sCH] = - w[r24] - w[r25] - w[r26] - w[r27f]
		+ w[r37f] - w[r76];

	cdot[sHCO] = w[r25] + w[r26] - w[r28f] - w[r29]
		- w[r32] - w[r35] + w[r48] + w[r49]
		+ w[r50] + w[r51] + w[r54f] + 2 * w[r90f]
		+ w[r91f] + w[ra125] + w[r144] + 2 * w[r243]
		+ w[rs339];

	cdot[sCH2OH] = w[r27f] - w[rm69] - w[rm70];

	cdot[s3XCH2] = - w[r37f] - w[r39] - w[r40] - 2 * w[r42]
		- 2 * w[r43] - w[r44] + w[r45f] + w[r56f]
		+ w[r79f] + w[r86] + w[r88] - w[r142f]
		+ w[rs339];

	cdot[sC2H2] = w[r42] + w[r43] - w[r84f] - w[r85f]
		- w[r86] - w[r87] + w[r92f] + w[r125]
		+ w[r113f] - w[r141f] - w[r142f] - w[r143]
		+ w[r157f] - w[r207f] + w[r208] - 2 * w[r211f]
		- w[r221f] - w[r234f] - w[rs257f] - w[rs268f]
		+ w[rs276] + w[rsA276] + w[rs288f] + w[rs300f]
		- w[rs303f] - w[rAs303f] - w[rs313f] + w[rs319f]
		- w[rs326f] - w[rs333f];

	cdot[sCH3] = - w[r44] + w[r46f] - w[r56f] - w[r57]
		- w[r58] - w[r59] + w[r60] - w[r61]
		- w[r62f] - 2 * w[r63] - 2 * w[r64] - w[r65]
		+ w[r71f] + w[r73f] + w[r89f] - w[r117f]
		- w[r129f] - w[r149] + w[r157f] - w[r168]
		+ w[r169f] - w[r172f] + w[r178f] - w[r245f]
		+ w[r255f] + w[r335] - w[r369] - w[r427]
		- w[rs310f] - w[rs318f] - w[rs325f] - w[rs332f];

	cdot[sC2H4] = w[r44] + w[r64] + w[r76] - w[r113f]
		- w[r115f] - w[r116f] - w[r117f] + w[r119f]
		+ w[r121] + w[r122] + w[r178f] + w[rm242]
		+ w[r292] + w[r336] + w[r381] + w[r387]
		+ w[r444] + w[r456] + w[r468] + w[r485];

	cdot[s1XCH2] = - w[r45f] - w[r46f] - w[r141f];

	cdot[sCH2O] = - w[r48] - w[r49] - w[r50] - w[r51]
		- w[r54f] + w[r57] + w[r58] + w[r66]
		+ w[r67] + w[r68f] + w[rm69] + w[rm70]
		+ w[r91f] + w[ra125] + w[r145] + w[r153]
		+ w[r164] + w[r176] + w[r265] + w[rl230];

	cdot[sCH3O] = w[r59] - w[r60] + w[r61] - w[r66]
		- w[r67] - w[r68f];

	cdot[sCH4] = w[r62f] + w[r65] - w[r71f] - w[r73f]
		- w[r76] + w[r117f] + w[r129f] + w[r149]
		+ w[r168] + w[r172f] + w[r369] + w[r427]
		+ w[rs310f] + w[rs318f] + w[rs325f] + w[rs332f];

	cdot[sC2H6] = w[r63] + w[r121] - w[r126f] - w[r128]
		- w[r129f];

	cdot[sC2H] = - w[r78] + w[r84f] + w[r85f] + w[r145]
		- w[r207f] - w[rs300f];

	cdot[sHCCO] = w[r78] - w[r79f] - w[r80] + w[r87]
		- w[r143] + w[r208] + w[rs290] + w[rs334];

	cdot[sCH2CO] = - w[r88] - w[r89f] - w[r90f] - w[r91f]
		+ w[r144] + w[rl230];

	cdot[sC2H3] = - w[r92f] - w[r125] - w[ra125] + w[r115f]
		+ w[r116f] + w[r117f] + w[r153] + w[r164]
		+ w[r169f] - w[r221f] - w[r234f] + w[r243]
		+ w[rm242] + w[r334] + w[r385];

	cdot[sC2H5] = - w[r119f] - 2 * w[r121] - w[r122] + w[r126f]
		+ w[r128] + w[r129f] + w[r176] + w[r292]
		+ w[r333] + w[r347] + w[r386] + w[r401];

	cdot[sC3H3] = w[r141f] + w[r142f] + w[r143] - w[r144]
		- w[r145] + w[r147f] + w[r148] + w[r149]
		+ w[r158f] + w[r159] - w[r245f] - 2 * w[rs260f]
		- w[rs261f] - 2 * w[rs271] - w[rs272];

	cdot[sC3H4] = - w[r147f] - w[r148] - w[r149] - w[r153]
		- w[r155f] + w[r166] + w[r167] + w[r168]
		+ w[r193] + w[r333] + w[r383] - w[rs261f]
		- w[rs272];

	cdot[sPXC3H4] = w[r155f] - w[r157f] - w[r158f] - w[r159]
		- w[r164];

	cdot[sC3H5] = - w[r166] - w[r167] - w[r168] - w[r193]
		+ w[r170f] + w[r172f] - w[r252] + w[r255f]
		+ w[r336] + w[r347] + w[r384] + w[r400]
		+ w[rc395] - w[rc396] + w[rc424] + w[rc425];

	cdot[sC3H6] = - w[r169f] - w[r170f] - w[r172f] - w[r176]
		+ w[r179f] + w[r180f] + w[r252] + w[r334]
		+ w[r384] + w[r481];

	cdot[sNXC3H7] = - w[r178f] - w[r179f] - w[r180f] + w[r265]
		+ w[r381] + w[r383] + w[r400];

	cdot[sC4H2] = w[r207f] - w[r208] + w[r213f] + w[r215f];

	cdot[sUXC4H3] = w[r211f] + w[r224f] + w[r226f] - w[rs257f]
		+ w[rsA276];

	cdot[sSXC4H3] = - w[r213f] - w[r215f] + w[r225f] + w[r227];

	cdot[sC4H4] = w[r221f] - w[r224f] - w[r225f] - w[r226f]
		- w[r227] - w[r235f] - w[r236f] + w[r237]
		+ w[r238] + w[r240] + w[r242] + w[rs276];

	cdot[sUXC4H5] = w[r234f] + w[r236f] - w[r238] - w[r239f]
		- w[r240] - w[rs268f] + w[rsA40f] + w[rsA50f];

	cdot[sSXC4H5] = w[r235f] - w[r237] + w[r239f] - w[r242]
		- w[r243] + w[r248f] + w[r250];

	cdot[sC4H6] = w[r245f] - w[r248f] - w[r250] + w[r241f]
		+ w[r244] + w[r252] + w[r335] + w[r386];

	cdot[sC4H7] = - w[r241f] - w[rm242] - w[r244] - w[r252]
		+ w[r258] + w[r269] + w[r279] + w[r387]
		+ w[r401];

	cdot[s1XC4H8] = - w[r255f] - w[r258] - w[r265] - w[r269]
		- w[r279] + w[r294] + w[r385] + w[rc424]
		+ w[rc425] + w[r482];

	cdot[sPXC4H9] = - w[r292] - w[r294] + w[r444] + w[rc395]
		- w[rc396] + w[r484];

	cdot[s1XC5H9] = - w[r333] - w[r334] - w[r335] - w[r336]
		+ w[r354] + w[r364] + w[r369] + w[r374];

	cdot[s1XC5H10] = - w[r347] - w[r354] - w[r364] - w[r369]
		- w[r374] + w[r483];

	cdot[s1XC5H11] = - w[r381] + w[r456] + w[r483] + 2 * w[r496];

	cdot[sC6H11] = - w[r383] - w[r384] - w[r385] - w[r386]
		- w[r387] + w[r409] + w[r421] + w[r427]
		+ w[r433];

	cdot[s1XC6H12] = - w[r400] - w[r401] - w[r409] - w[r421]
		- w[r427] - w[r433] + w[r484];

	cdot[s1XC6H13] = - w[r444] + w[r468] + w[r482];

	cdot[s1XC7H15] = - w[r456] + w[r481] + w[rl230];

	cdot[s1XC7H14] = - w[rc395] + w[rc396] - w[rc424] - w[rc425];

	cdot[s1XC8H17] = - w[r468] + w[r485];

	cdot[s2XC10H21] = - w[r481] + w[r501] + w[r506] + w[r516]
		+ w[r526] - w[rl20f] + w[rl20b];

	cdot[s3XC10H21] = - w[r482] + w[r500] + w[r505] + w[r515]
		+ w[r525] - w[rl30f] + w[rl30b];

	cdot[s4XC10H21] = - w[r483] + w[r499] + w[r504] + w[r514]
		+ w[r524] - w[rl40f] + w[rl40b];

	cdot[s5XC10H21] = - w[r484] + w[r498] + w[r503] + w[r513]
		+ w[r523] - w[rl50f] + w[rl50b];

	cdot[s1XC10H21] = - w[r485] + w[r497] + w[r502] + w[r512]
		+ w[r522] - w[rl10f] + w[rl10b];

	cdot[sNXC10H22] = - w[r496] - w[r497] - w[r498] - w[r499]
		- w[r500] - w[r501] - w[r502] - w[r503]
		- w[r504] - w[r505] - w[r506] - w[r512]
		- w[r513] - w[r514] - w[r515] - w[r516]
		- w[r522] - w[r523] - w[r524] - w[r525]
		- w[r526];

	cdot[s2XRO2XC10H21O2] = w[rl10f] - w[rl10b] + w[rl20f] - w[rl20b]
		+ w[rl30f] - w[rl30b] + w[rl40f] - w[rl40b]
		+ w[rl50f] - w[rl50b] - w[rl60f] + w[rl60b];

	cdot[s3X2XRO2HXC10H20OOH] = w[rl60f] - w[rl60b] - w[rl70];

	cdot[s32O2RO2HXO2C10H20OOH] = w[rl70] - w[rl80];

	cdot[s3X2HOROHXHO2C10H19OOH] = w[rl80] - w[rl90];

	cdot[s3X2ORO2HXOC10H19OOH] = w[rl90] - w[rl100];

	cdot[s2X1OROXOC10H19O] = w[rl100] - w[rl230];

	cdot[sUXC6H5] = w[rs257f] + w[rs274f];

	cdot[sC6H6] = w[rs260f] + w[rs261f] - w[rs266f];

	cdot[sUXC6H7] = w[rs266f] + w[rs268f] - w[rs273];

	cdot[sA1XC6H6] = w[rs271] + w[rs272] + w[rs273] - w[rs276]
		- w[rs277f] - w[rs278f] - w[rs279f] - w[rs280f]
		- w[rsA280f] - w[rs281f] + w[rsA284f] + w[rs290];

	cdot[sA1XXC6H5] = - w[rs274f] - w[rsA276] + w[rs277f] + w[rs278f]
		+ w[rs279f] + w[rs280f] + w[rsA280f] - w[rs2A281f]
		- w[rs275] + w[rs288f];

	cdot[sC6H5O] = w[rs281f] + w[rs2A281f] + w[rs275] - w[rs3A281f]
		+ w[rs284f] + w[rs285f] + w[rsA285f] - w[rsA33f];

	cdot[sC6H5OH] = w[rs3A281f] - w[rs284f] - w[rsA284f] - w[rs285f]
		- w[rsA285f];

	cdot[sC5H5] = w[rsA33f] - w[rsA39f] - w[rsA40f] - w[rsA41f]
		- w[rsA43f] + w[rsA44f] + w[rsA45f] + w[rsA46f];

	cdot[sC5H6] = w[rsA39f] - w[rsA44f] - w[rsA45f] - w[rsA46f]
		- w[rsA49f];

	cdot[sC5H5O] = w[rsA41f] + w[rsA43f] + w[rsA49f] - w[rsA50f];

	cdot[sA1C2HXC8H6] = - w[rs288f] - w[rs290] + w[rs297f] - w[rs298f]
		- w[rs299f] - w[rs300f];

	cdot[sA1C2HYXC8H5] = - w[rs297f] + w[rs298f] + w[rs299f] + w[rs300f]
		- w[rs303f] - w[rAs303f];

	cdot[sA1C2HACXC10H7] = w[rs303f] - w[rs304f];

	cdot[sA2XXXC10H7] = w[rAs303f] + w[rs304f] + w[rs306f] + w[rs307f]
		+ w[rs308f] + w[rs310f] - w[rs313f] + w[rs335];

	cdot[sA2XC10H8] = - w[rs306f] - w[rs307f] - w[rs308f] - w[rs310f]
		+ w[rs334];

	cdot[sA2R5XC12H8] = w[rs313f] - w[rs315f] - w[rs318f] - w[rs334];

	cdot[sA2R5XXC12H7] = w[rs315f] + w[rs318f] + w[rs319f] - w[rs335];

	cdot[sA2R5C2HXC14H8] = - w[rs319f] - w[rs322f] - w[rs325f] + w[rs339];

	cdot[sA2R5C2HYXC14H7] = w[rs322f] + w[rs325f] - w[rs326f];

	cdot[sANC2HACXC16H9] = w[rs326f] - w[rs327f];

	cdot[sA3R5XXC16H9] = w[rs327f] + w[rs329f] + w[rs332f] - w[rs333f];

	cdot[sA3R5XC16H10] = - w[rs329f] - w[rs332f] - w[rs339];

	cdot[sA3R5ACXC18H11] = w[rs333f] - w[rs342f] + w[rs342b];

	cdot[sA4XC18H10] = w[rs342f] - w[rs342b];

	cdot[sA2CH3XC11H10] = - w[rF1f] - w[rAF1f] - w[rAF2f] - w[rF3f]
		- w[rF4f];

	cdot[sA2CH2XC11H9] = w[rF1f] + w[rAF1f] + w[rAF2f] + w[rF3f]
		+ w[rF4f] - w[rF7f] - w[rF9] - w[rAF12f]
		- 2 * w[rAF13f] + 2 * w[rAF13b];

	cdot[sA2CHOXC11H8O] = w[rF7f] + w[rF9] + w[rAF12f];

	cdot[sBIBENZYLXC22H18] = w[rAF13f] - w[rAF13b];

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
	W[sN2] =  2.80140000e+01;
	W[sO2] =  3.20000000e+01;
	W[sH] =  1.00800000e+00;
	W[sOH] =  1.70080000e+01;
	W[sO] =  1.60000000e+01;
	W[sH2] =  2.01600000e+00;
	W[sH2O] =  1.80160000e+01;
	W[sHO2] =  3.30080000e+01;
	W[sH2O2] =  3.40160000e+01;
	W[sCO] =  2.80110000e+01;
	W[sCO2] =  4.40110000e+01;
	W[sCH] =  1.30190000e+01;
	W[sHCO] =  2.90190000e+01;
	W[sCH2OH] =  3.10270000e+01;
	W[s3XCH2] =  1.40270000e+01;
	W[sC2H2] =  2.60460000e+01;
	W[sCH3] =  1.50350000e+01;
	W[sC2H4] =  2.80620000e+01;
	W[s1XCH2] =  1.40270000e+01;
	W[sCH2O] =  3.00270000e+01;
	W[sCH3O] =  3.10340000e+01;
	W[sCH4] =  1.60430000e+01;
	W[sC2H6] =  3.00780000e+01;
	W[sC2H] =  2.50380000e+01;
	W[sHCCO] =  4.10300000e+01;
	W[sCH2CO] =  4.20380000e+01;
	W[sC2H3] =  2.70540000e+01;
	W[sC2H5] =  2.90700000e+01;
	W[sC3H3] =  3.90650000e+01;
	W[sC3H4] =  4.00730000e+01;
	W[sPXC3H4] =  4.00730000e+01;
	W[sC3H5] =  4.10810000e+01;
	W[sC3H6] =  4.20500000e+01;
	W[sNXC3H7] =  4.30600000e+01;
	W[sC4H2] =  5.00680000e+01;
	W[sUXC4H3] =  5.10640000e+01;
	W[sSXC4H3] =  5.10760000e+01;
	W[sC4H4] =  5.20840000e+01;
	W[sUXC4H5] =  5.30800000e+01;
	W[sSXC4H5] =  5.30920000e+01;
	W[sC4H6] =  5.40980000e+01;
	W[sC4H7] =  5.51040000e+01;
	W[s1XC4H8] =  5.60600000e+01;
	W[sPXC4H9] =  5.71120000e+01;
	W[s1XC5H9] =  6.91220000e+01;
	W[s1XC5H10] =  7.01300000e+01;
	W[s1XC5H11] =  7.11380000e+01;
	W[sC6H11] =  8.31480000e+01;
	W[s1XC6H12] =  8.41560000e+01;
	W[s1XC6H13] =  8.51640000e+01;
	W[s1XC7H15] =  9.91900000e+01;
	W[s1XC7H14] =  9.81820000e+01;
	W[s1XC8H17] =  1.13216000e+02;
	W[s2XC10H21] =  1.41268000e+02;
	W[s3XC10H21] =  1.41268000e+02;
	W[s4XC10H21] =  1.41268000e+02;
	W[s5XC10H21] =  1.41268000e+02;
	W[s1XC10H21] =  1.41268000e+02;
	W[sNXC10H22] =  1.42276000e+02;
	W[s2XRO2XC10H21O2] =  1.73268000e+02;
	W[s3X2XRO2HXC10H20OOH] =  1.73268000e+02;
	W[s32O2RO2HXO2C10H20OOH] =  2.05268000e+02;
	W[s3X2HOROHXHO2C10H19OOH] =  2.05268000e+02;
	W[s3X2ORO2HXOC10H19OOH] =  1.88260000e+02;
	W[s2X1OROXOC10H19O] =  1.71252000e+02;
	W[sUXC6H5] =  7.71000000e+01;
	W[sC6H6] =  7.81080000e+01;
	W[sUXC6H7] =  7.91160000e+01;
	W[sA1XC6H6] =  7.81080000e+01;
	W[sA1XXC6H5] =  7.71000000e+01;
	W[sC6H5O] =  9.31000000e+01;
	W[sC6H5OH] =  9.41080000e+01;
	W[sC5H5] =  6.50900000e+01;
	W[sC5H6] =  6.60980000e+01;
	W[sC5H5O] =  8.10900000e+01;
	W[sA1C2HXC8H6] =  1.02128000e+02;
	W[sA1C2HYXC8H5] =  1.01120000e+02;
	W[sA1C2HACXC10H7] =  1.27156000e+02;
	W[sA2XXXC10H7] =  1.27156000e+02;
	W[sA2XC10H8] =  1.28164000e+02;
	W[sA2R5XC12H8] =  1.52184000e+02;
	W[sA2R5XXC12H7] =  1.51176000e+02;
	W[sA2R5C2HXC14H8] =  1.76204000e+02;
	W[sA2R5C2HYXC14H7] =  1.75196000e+02;
	W[sANC2HACXC16H9] =  2.01232000e+02;
	W[sA3R5XXC16H9] =  2.01232000e+02;
	W[sA3R5XC16H10] =  2.02240000e+02;
	W[sA3R5ACXC18H11] =  2.27268000e+02;
	W[sA4XC18H10] =  2.26260000e+02;
	W[sA2CH3XC11H10] =  1.42190000e+02;
	W[sA2CH2XC11H9] =  1.41182000e+02;
	W[sA2CHOXC11H8O] =  1.56174000e+02;
	W[sBIBENZYLXC22H18] =  2.82364000e+02;
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
	strcpy( names[sCH3O], "CH3O                " );
	strcpy( names[sCH4], "CH4                 " );
	strcpy( names[sC2H6], "C2H6                " );
	strcpy( names[sC2H], "C2H                 " );
	strcpy( names[sHCCO], "HCCO                " );
	strcpy( names[sCH2CO], "CH2CO               " );
	strcpy( names[sC2H3], "C2H3                " );
	strcpy( names[sC2H5], "C2H5                " );
	strcpy( names[sC3H3], "C3H3                " );
	strcpy( names[sC3H4], "C3H4                " );
	strcpy( names[sPXC3H4], "P-C3H4              " );
	strcpy( names[sC3H5], "C3H5                " );
	strcpy( names[sC3H6], "C3H6                " );
	strcpy( names[sNXC3H7], "N-C3H7              " );
	strcpy( names[sC4H2], "C4H2                " );
	strcpy( names[sUXC4H3], "U-C4H3              " );
	strcpy( names[sSXC4H3], "S-C4H3              " );
	strcpy( names[sC4H4], "C4H4                " );
	strcpy( names[sUXC4H5], "U-C4H5              " );
	strcpy( names[sSXC4H5], "S-C4H5              " );
	strcpy( names[sC4H6], "C4H6                " );
	strcpy( names[sC4H7], "C4H7                " );
	strcpy( names[s1XC4H8], "1-C4H8              " );
	strcpy( names[sPXC4H9], "P-C4H9              " );
	strcpy( names[s1XC5H9], "1-C5H9              " );
	strcpy( names[s1XC5H10], "1-C5H10             " );
	strcpy( names[s1XC5H11], "1-C5H11             " );
	strcpy( names[sC6H11], "C6H11               " );
	strcpy( names[s1XC6H12], "1-C6H12             " );
	strcpy( names[s1XC6H13], "1-C6H13             " );
	strcpy( names[s1XC7H15], "1-C7H15             " );
	strcpy( names[s1XC7H14], "1-C7H14             " );
	strcpy( names[s1XC8H17], "1-C8H17             " );
	strcpy( names[s2XC10H21], "2-C10H21            " );
	strcpy( names[s3XC10H21], "3-C10H21            " );
	strcpy( names[s4XC10H21], "4-C10H21            " );
	strcpy( names[s5XC10H21], "5-C10H21            " );
	strcpy( names[s1XC10H21], "1-C10H21            " );
	strcpy( names[sNXC10H22], "N-C10H22            " );
	strcpy( names[s2XRO2XC10H21O2], "2-RO2-C10H21O2      " );
	strcpy( names[s3X2XRO2HXC10H20OOH], "3-2-RO2H-C10H20OOH  " );
	strcpy( names[s32O2RO2HXO2C10H20OOH], "32O2RO2H-O2C10H20OOH" );
	strcpy( names[s3X2HOROHXHO2C10H19OOH], "3-2HOROH-HO2C10H19OOH" );
	strcpy( names[s3X2ORO2HXOC10H19OOH], "3-2ORO2H-OC10H19OOH " );
	strcpy( names[s2X1OROXOC10H19O], "2-1ORO-OC10H19O     " );
	strcpy( names[sUXC6H5], "U-C6H5              " );
	strcpy( names[sC6H6], "C6H6                " );
	strcpy( names[sUXC6H7], "U-C6H7              " );
	strcpy( names[sA1XC6H6], "A1-C6H6             " );
	strcpy( names[sA1XXC6H5], "A1--C6H5            " );
	strcpy( names[sC6H5O], "C6H5O               " );
	strcpy( names[sC6H5OH], "C6H5OH              " );
	strcpy( names[sC5H5], "C5H5                " );
	strcpy( names[sC5H6], "C5H6                " );
	strcpy( names[sC5H5O], "C5H5O               " );
	strcpy( names[sA1C2HXC8H6], "A1C2H-C8H6          " );
	strcpy( names[sA1C2HYXC8H5], "A1C2H*-C8H5         " );
	strcpy( names[sA1C2HACXC10H7], "A1C2HAC-C10H7       " );
	strcpy( names[sA2XXXC10H7], "A2-X-C10H7          " );
	strcpy( names[sA2XC10H8], "A2-C10H8            " );
	strcpy( names[sA2R5XC12H8], "A2R5-C12H8          " );
	strcpy( names[sA2R5XXC12H7], "A2R5--C12H7         " );
	strcpy( names[sA2R5C2HXC14H8], "A2R5C2H-C14H8       " );
	strcpy( names[sA2R5C2HYXC14H7], "A2R5C2H*-C14H7      " );
	strcpy( names[sANC2HACXC16H9], "ANC2HAC-C16H9       " );
	strcpy( names[sA3R5XXC16H9], "A3R5--C16H9         " );
	strcpy( names[sA3R5XC16H10], "A3R5-C16H10         " );
	strcpy( names[sA3R5ACXC18H11], "A3R5AC-C18H11       " );
	strcpy( names[sA4XC18H10], "A4-C18H10           " );
	strcpy( names[sA2CH3XC11H10], "A2CH3-C11H10        " );
	strcpy( names[sA2CH2XC11H9], "A2CH2-C11H9         " );
	strcpy( names[sA2CHOXC11H8O], "A2CHO-C11H8O        " );
	strcpy( names[sBIBENZYLXC22H18], "BIBENZYL-C22H18     " );
}



void ComputeThermoData( double *h, double *cp, double T )
{
/*
	This function computes enthalpy 'h' and heat capacity 'cp' as
	function of temperature 'T' for all non steady state species
	in units [J/kg] and [J/kg K], respectively.
	The parameter h and cp should provide workspace of length 93 */


	int i;
	if ( T > 1000.0 ) {
		h[sN2] =  2.96792318e+02 * (
			T * (  2.92664000e+00 + T * (  7.43988400e-04
			+ T * ( -1.89492000e-07 + T * (  2.52425950e-11
			+ T * -1.35067020e-15 ) ) ) ) -9.22797700e+02 );
		cp[sN2] =  2.96792318e+02 * (
			 2.92664000e+00 + T * (  1.48797680e-03 
			+ T * ( -5.68476000e-07 + T * (  1.00970380e-10
			+ T * -6.75335100e-15 ) ) ) );
		h[sO2] =  2.59823125e+02 * (
			T * (  3.69757800e+00 + T * (  3.06759850e-04
			+ T * ( -4.19614000e-08 + T * (  4.43820250e-12
			+ T * -2.27287080e-16 ) ) ) ) -1.23393010e+03 );
		cp[sO2] =  2.59823125e+02 * (
			 3.69757800e+00 + T * (  6.13519700e-04 
			+ T * ( -1.25884200e-07 + T * (  1.77528100e-11
			+ T * -1.13643540e-15 ) ) ) );
		h[sH] =  8.24835317e+03 * (
			T * (  2.50000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  2.54716200e+04 );
		cp[sH] =  8.24835317e+03 * (
			 2.50000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sOH] =  4.88848777e+02 * (
			T * (  2.88273000e+00 + T * (  5.06987150e-04
			+ T * ( -7.58959000e-08 + T * (  5.43670750e-12
			+ T * -1.02526100e-16 ) ) ) ) +  3.88688800e+03 );
		cp[sOH] =  4.88848777e+02 * (
			 2.88273000e+00 + T * (  1.01397430e-03 
			+ T * ( -2.27687700e-07 + T * (  2.17468300e-11
			+ T * -5.12630500e-16 ) ) ) );
		h[sO] =  5.19646250e+02 * (
			T * (  2.54205900e+00 + T * ( -1.37753050e-05
			+ T * ( -1.03426767e-09 + T * (  1.13776675e-12
			+ T * -8.73610200e-17 ) ) ) ) +  2.92308000e+04 );
		cp[sO] =  5.19646250e+02 * (
			 2.54205900e+00 + T * ( -2.75506100e-05 
			+ T * ( -3.10280300e-09 + T * (  4.55106700e-12
			+ T * -4.36805100e-16 ) ) ) );
		h[sH2] =  4.12417659e+03 * (
			T * (  2.99142300e+00 + T * (  3.50032200e-04
			+ T * ( -1.87794267e-08 + T * ( -2.30789450e-12
			+ T *  3.16550380e-16 ) ) ) ) -8.35034000e+02 );
		cp[sH2] =  4.12417659e+03 * (
			 2.99142300e+00 + T * (  7.00064400e-04 
			+ T * ( -5.63382800e-08 + T * ( -9.23157800e-12
			+ T *  1.58275190e-15 ) ) ) );
		h[sH2O] =  4.61497558e+02 * (
			T * (  2.67214500e+00 + T * (  1.52814650e-03
			+ T * ( -2.91008667e-07 + T * (  3.00249100e-11
			+ T * -1.27832360e-15 ) ) ) ) -2.98992100e+04 );
		cp[sH2O] =  4.61497558e+02 * (
			 2.67214500e+00 + T * (  3.05629300e-03 
			+ T * ( -8.73026000e-07 + T * (  1.20099640e-10
			+ T * -6.39161800e-15 ) ) ) );
		h[sHO2] =  2.51888633e+02 * (
			T * (  4.07219100e+00 + T * (  1.06564800e-03
			+ T * ( -1.76938167e-07 + T * (  1.52806725e-11
			+ T * -5.68232800e-16 ) ) ) ) -1.57972700e+02 );
		cp[sHO2] =  2.51888633e+02 * (
			 4.07219100e+00 + T * (  2.13129600e-03 
			+ T * ( -5.30814500e-07 + T * (  6.11226900e-11
			+ T * -2.84116400e-15 ) ) ) );
		h[sH2O2] =  2.44424389e+02 * (
			T * (  4.57316700e+00 + T * (  2.16806800e-03
			+ T * ( -4.91562933e-07 + T * (  5.87225750e-11
			+ T * -2.86330720e-15 ) ) ) ) -1.80069600e+04 );
		cp[sH2O2] =  2.44424389e+02 * (
			 4.57316700e+00 + T * (  4.33613600e-03 
			+ T * ( -1.47468880e-06 + T * (  2.34890300e-10
			+ T * -1.43165360e-14 ) ) ) );
		h[sCO] =  2.96824105e+02 * (
			T * (  3.02507800e+00 + T * (  7.21344250e-04
			+ T * ( -1.87694233e-07 + T * (  2.54645325e-11
			+ T * -1.38219020e-15 ) ) ) ) -1.42683500e+04 );
		cp[sCO] =  2.96824105e+02 * (
			 3.02507800e+00 + T * (  1.44268850e-03 
			+ T * ( -5.63082700e-07 + T * (  1.01858130e-10
			+ T * -6.91095100e-15 ) ) ) );
		h[sCO2] =  1.88915044e+02 * (
			T * (  4.45362300e+00 + T * (  1.57008400e-03
			+ T * ( -4.26136833e-07 + T * (  5.98499000e-11
			+ T * -3.33806660e-15 ) ) ) ) -4.89669600e+04 );
		cp[sCO2] =  1.88915044e+02 * (
			 4.45362300e+00 + T * (  3.14016800e-03 
			+ T * ( -1.27841050e-06 + T * (  2.39399600e-10
			+ T * -1.66903330e-14 ) ) ) );
		h[sCH] =  6.38631231e+02 * (
			T * (  2.19622300e+00 + T * (  1.17019050e-03
			+ T * ( -2.35273367e-07 + T * (  2.25189550e-11
			+ T * -7.71008000e-16 ) ) ) ) +  7.08672300e+04 );
		cp[sCH] =  6.38631231e+02 * (
			 2.19622300e+00 + T * (  2.34038100e-03 
			+ T * ( -7.05820100e-07 + T * (  9.00758200e-11
			+ T * -3.85504000e-15 ) ) ) );
		h[sHCO] =  2.86513663e+02 * (
			T * (  3.55727100e+00 + T * (  1.67278600e-03
			+ T * ( -4.45002000e-07 + T * (  6.17643000e-11
			+ T * -3.42770000e-15 ) ) ) ) +  3.91632400e+03 );
		cp[sHCO] =  2.86513663e+02 * (
			 3.55727100e+00 + T * (  3.34557200e-03 
			+ T * ( -1.33500600e-06 + T * (  2.47057200e-10
			+ T * -1.71385000e-14 ) ) ) );
		h[sCH2OH] =  2.67971122e+02 * (
			T * (  6.32752000e+00 + T * (  1.80413500e-03
			+ T * ( -1.06718233e-07 + T * ( -4.84687500e-11
			+ T *  7.01940800e-15 ) ) ) ) -4.47450900e+03 );
		cp[sCH2OH] =  2.67971122e+02 * (
			 6.32752000e+00 + T * (  3.60827000e-03 
			+ T * ( -3.20154700e-07 + T * ( -1.93875000e-10
			+ T *  3.50970400e-14 ) ) ) );
		h[s3XCH2] =  5.92738290e+02 * (
			T * (  3.63640700e+00 + T * (  9.66528000e-04
			+ T * ( -5.62338667e-08 + T * ( -2.52474850e-11
			+ T *  3.61651000e-15 ) ) ) ) +  4.53413400e+04 );
		cp[s3XCH2] =  5.92738290e+02 * (
			 3.63640700e+00 + T * (  1.93305600e-03 
			+ T * ( -1.68701600e-07 + T * ( -1.00989940e-10
			+ T *  1.80825500e-14 ) ) ) );
		h[sC2H2] =  3.19217538e+02 * (
			T * (  4.43677000e+00 + T * (  2.68801950e-03
			+ T * ( -6.37605333e-07 + T * (  8.21594750e-11
			+ T * -4.31341800e-15 ) ) ) ) +  2.56676600e+04 );
		cp[sC2H2] =  3.19217538e+02 * (
			 4.43677000e+00 + T * (  5.37603900e-03 
			+ T * ( -1.91281600e-06 + T * (  3.28637900e-10
			+ T * -2.15670900e-14 ) ) ) );
		h[sCH3] =  5.52999002e+02 * (
			T * (  2.84405100e+00 + T * (  3.06898700e-03
			+ T * ( -7.43448333e-07 + T * (  9.46290250e-11
			+ T * -4.90431800e-15 ) ) ) ) +  1.64378090e+04 );
		cp[sCH3] =  5.52999002e+02 * (
			 2.84405100e+00 + T * (  6.13797400e-03 
			+ T * ( -2.23034500e-06 + T * (  3.78516100e-10
			+ T * -2.45215900e-14 ) ) ) );
		h[sC2H4] =  2.96284655e+02 * (
			T * (  3.52841800e+00 + T * (  5.74259250e-03
			+ T * ( -1.47279500e-06 + T * (  1.96115000e-10
			+ T * -1.05336960e-14 ) ) ) ) +  4.42828800e+03 );
		cp[sC2H4] =  2.96284655e+02 * (
			 3.52841800e+00 + T * (  1.14851850e-02 
			+ T * ( -4.41838500e-06 + T * (  7.84460000e-10
			+ T * -5.26684800e-14 ) ) ) );
		h[s1XCH2] =  5.92738290e+02 * (
			T * (  3.55288800e+00 + T * (  1.03339400e-03
			+ T * ( -6.38038667e-08 + T * ( -2.76168325e-11
			+ T *  4.04269800e-15 ) ) ) ) +  4.98497500e+04 );
		cp[s1XCH2] =  5.92738290e+02 * (
			 3.55288800e+00 + T * (  2.06678800e-03 
			+ T * ( -1.91411600e-07 + T * ( -1.10467330e-10
			+ T *  2.02134900e-14 ) ) ) );
		h[sCH2O] =  2.76895461e+02 * (
			T * (  2.99560600e+00 + T * (  3.34066050e-03
			+ T * ( -8.76318000e-07 + T * (  1.18428825e-10
			+ T * -6.42503400e-15 ) ) ) ) -1.53203690e+04 );
		cp[sCH2O] =  2.76895461e+02 * (
			 2.99560600e+00 + T * (  6.68132100e-03 
			+ T * ( -2.62895400e-06 + T * (  4.73715300e-10
			+ T * -3.21251700e-14 ) ) ) );
		h[sCH3O] =  2.67910679e+02 * (
			T * (  3.77079900e+00 + T * (  3.93574850e-03
			+ T * ( -8.85461333e-07 + T * (  9.86107750e-11
			+ T * -4.22523200e-15 ) ) ) ) +  1.27832520e+02 );
		cp[sCH3O] =  2.67910679e+02 * (
			 3.77079900e+00 + T * (  7.87149700e-03 
			+ T * ( -2.65638400e-06 + T * (  3.94443100e-10
			+ T * -2.11261600e-14 ) ) ) );
		h[sCH4] =  5.18253444e+02 * (
			T * (  1.68347800e+00 + T * (  5.11861800e-03
			+ T * ( -1.29170933e-06 + T * (  1.69639625e-10
			+ T * -9.00684600e-15 ) ) ) ) -1.00807870e+04 );
		cp[sCH4] =  5.18253444e+02 * (
			 1.68347800e+00 + T * (  1.02372360e-02 
			+ T * ( -3.87512800e-06 + T * (  6.78558500e-10
			+ T * -4.50342300e-14 ) ) ) );
		h[sC2H6] =  2.76425959e+02 * (
			T * (  4.82593800e+00 + T * (  6.92021450e-03
			+ T * ( -1.51908600e-06 + T * (  1.68124175e-10
			+ T * -7.19632200e-15 ) ) ) ) -1.27177930e+04 );
		cp[sC2H6] =  2.76425959e+02 * (
			 4.82593800e+00 + T * (  1.38404290e-02 
			+ T * ( -4.55725800e-06 + T * (  6.72496700e-10
			+ T * -3.59816100e-14 ) ) ) );
		h[sC2H] =  3.32068855e+02 * (
			T * (  4.42768800e+00 + T * (  1.10813400e-03
			+ T * ( -2.01631733e-07 + T * (  2.47062925e-11
			+ T * -1.47023580e-15 ) ) ) ) +  6.59041500e+04 );
		cp[sC2H] =  3.32068855e+02 * (
			 4.42768800e+00 + T * (  2.21626800e-03 
			+ T * ( -6.04895200e-07 + T * (  9.88251700e-11
			+ T * -7.35117900e-15 ) ) ) );
		h[sHCCO] =  2.02640507e+02 * (
			T * (  6.75807300e+00 + T * (  1.00020000e-03
			+ T * ( -6.75869000e-08 + T * ( -2.60282950e-11
			+ T *  3.93032800e-15 ) ) ) ) +  1.90151300e+04 );
		cp[sHCCO] =  2.02640507e+02 * (
			 6.75807300e+00 + T * (  2.00040000e-03 
			+ T * ( -2.02760700e-07 + T * ( -1.04113180e-10
			+ T *  1.96516400e-14 ) ) ) );
		h[sCH2CO] =  1.97781531e+02 * (
			T * (  6.03881700e+00 + T * (  2.90242000e-03
			+ T * ( -6.40317667e-07 + T * (  6.98621000e-11
			+ T * -2.91773520e-15 ) ) ) ) -8.58340200e+03 );
		cp[sCH2CO] =  1.97781531e+02 * (
			 6.03881700e+00 + T * (  5.80484000e-03 
			+ T * ( -1.92095300e-06 + T * (  2.79448400e-10
			+ T * -1.45886760e-14 ) ) ) );
		h[sC2H3] =  3.07323871e+02 * (
			T * (  5.93346800e+00 + T * (  2.00887250e-03
			+ T * ( -1.32224633e-07 + T * ( -3.60316650e-11
			+ T *  4.75728600e-15 ) ) ) ) +  3.18543400e+04 );
		cp[sC2H3] =  3.07323871e+02 * (
			 5.93346800e+00 + T * (  4.01774500e-03 
			+ T * ( -3.96673900e-07 + T * ( -1.44126660e-10
			+ T *  2.37864300e-14 ) ) ) );
		h[sC2H5] =  2.86011008e+02 * (
			T * (  7.19048000e+00 + T * (  3.24203850e-03
			+ T * ( -2.14268800e-07 + T * ( -5.86969750e-11
			+ T *  7.76175400e-15 ) ) ) ) +  1.06745490e+04 );
		cp[sC2H5] =  2.86011008e+02 * (
			 7.19048000e+00 + T * (  6.48407700e-03 
			+ T * ( -6.42806400e-07 + T * ( -2.34787900e-10
			+ T *  3.88087700e-14 ) ) ) );
		h[sC3H3] =  2.12833483e+02 * (
			T * (  8.83104700e+00 + T * (  2.17859700e-03
			+ T * ( -1.36968867e-07 + T * ( -5.92180750e-11
			+ T *  8.75304000e-15 ) ) ) ) +  3.76322300e+04 );
		cp[sC3H3] =  2.12833483e+02 * (
			 8.83104700e+00 + T * (  4.35719400e-03 
			+ T * ( -4.10906600e-07 + T * ( -2.36872300e-10
			+ T *  4.37652000e-14 ) ) ) );
		h[sC3H4] =  2.07479849e+02 * (
			T * (  5.72914400e+00 + T * (  6.18402250e-03
			+ T * ( -1.60187533e-06 + T * (  2.15034100e-10
			+ T * -1.16256040e-14 ) ) ) ) +  2.01298400e+04 );
		cp[sC3H4] =  2.07479849e+02 * (
			 5.72914400e+00 + T * (  1.23680450e-02 
			+ T * ( -4.80562600e-06 + T * (  8.60136400e-10
			+ T * -5.81280200e-14 ) ) ) );
		h[sPXC3H4] =  2.07479849e+02 * (
			T * (  5.51103400e+00 + T * (  6.23478100e-03
			+ T * ( -1.60472133e-06 + T * (  2.14344250e-10
			+ T * -1.15431220e-14 ) ) ) ) +  1.99204000e+04 );
		cp[sPXC3H4] =  2.07479849e+02 * (
			 5.51103400e+00 + T * (  1.24695620e-02 
			+ T * ( -4.81416400e-06 + T * (  8.57377000e-10
			+ T * -5.77156100e-14 ) ) ) );
		h[sC3H5] =  2.02388939e+02 * (
			T * (  9.65153900e+00 + T * (  4.03779800e-03
			+ T * ( -2.65514133e-07 + T * ( -1.16267400e-10
			+ T *  1.72065620e-14 ) ) ) ) +  1.53009550e+04 );
		cp[sC3H5] =  2.02388939e+02 * (
			 9.65153900e+00 + T * (  8.07559600e-03 
			+ T * ( -7.96542400e-07 + T * ( -4.65069600e-10
			+ T *  8.60328100e-14 ) ) ) );
		h[sC3H6] =  1.97725089e+02 * (
			T * (  6.73225700e+00 + T * (  7.45416800e-03
			+ T * ( -1.64996633e-06 + T * (  1.80300550e-10
			+ T * -7.53240800e-15 ) ) ) ) -9.23570300e+02 );
		cp[sC3H6] =  1.97725089e+02 * (
			 6.73225700e+00 + T * (  1.49083360e-02 
			+ T * ( -4.94989900e-06 + T * (  7.21202200e-10
			+ T * -3.76620400e-14 ) ) ) );
		h[sNXC3H7] =  1.93087320e+02 * (
			T * (  7.97829000e+00 + T * (  7.88056700e-03
			+ T * ( -1.72441433e-06 + T * (  1.86097300e-10
			+ T * -7.64995600e-15 ) ) ) ) +  8.26499100e+03 );
		cp[sNXC3H7] =  1.93087320e+02 * (
			 7.97829000e+00 + T * (  1.57611340e-02 
			+ T * ( -5.17324300e-06 + T * (  7.44389200e-10
			+ T * -3.82497800e-14 ) ) ) );
		h[sC4H2] =  1.66060957e+02 * (
			T * (  9.03140740e+00 + T * (  3.02362630e-03
			+ T * ( -6.49596267e-07 + T * (  6.88715750e-11
			+ T * -2.77121600e-15 ) ) ) ) +  4.98634000e+04 );
		cp[sC4H2] =  1.66060957e+02 * (
			 9.03140740e+00 + T * (  6.04725260e-03 
			+ T * ( -1.94878880e-06 + T * (  2.75486300e-10
			+ T * -1.38560800e-14 ) ) ) );
		h[sUXC4H3] =  1.62821949e+02 * (
			T * (  1.07527380e+01 + T * (  2.69057650e-03
			+ T * ( -1.84987900e-07 + T * ( -7.63066500e-11
			+ T *  1.15234800e-14 ) ) ) ) +  5.88687600e+04 );
		cp[sUXC4H3] =  1.62821949e+02 * (
			 1.07527380e+01 + T * (  5.38115300e-03 
			+ T * ( -5.54963700e-07 + T * ( -3.05226600e-10
			+ T *  5.76174000e-14 ) ) ) );
		h[sSXC4H3] =  1.62783695e+02 * (
			T * (  1.13140950e+01 + T * (  2.50720700e-03
			+ T * ( -1.78348133e-07 + T * ( -7.06327250e-11
			+ T *  1.08065580e-14 ) ) ) ) +  5.18121100e+04 );
		cp[sSXC4H3] =  1.62783695e+02 * (
			 1.13140950e+01 + T * (  5.01441400e-03 
			+ T * ( -5.35044400e-07 + T * ( -2.82530900e-10
			+ T *  5.40327900e-14 ) ) ) );
		h[sC4H4] =  1.59633285e+02 * (
			T * (  1.06977730e+01 + T * (  3.49100700e-03
			+ T * ( -2.18924900e-07 + T * ( -9.71129250e-11
			+ T *  1.44018920e-14 ) ) ) ) +  3.03480300e+04 );
		cp[sC4H4] =  1.59633285e+02 * (
			 1.06977730e+01 + T * (  6.98201400e-03 
			+ T * ( -6.56774700e-07 + T * ( -3.88451700e-10
			+ T *  7.20094600e-14 ) ) ) );
		h[sUXC4H5] =  1.56637905e+02 * (
			T * (  1.28659710e+01 + T * (  3.97168450e-03
			+ T * ( -2.87548867e-07 + T * ( -1.16390875e-10
			+ T *  1.79022620e-14 ) ) ) ) +  3.68546400e+04 );
		cp[sUXC4H5] =  1.56637905e+02 * (
			 1.28659710e+01 + T * (  7.94336900e-03 
			+ T * ( -8.62646600e-07 + T * ( -4.65563500e-10
			+ T *  8.95113100e-14 ) ) ) );
		h[sSXC4H5] =  1.56602501e+02 * (
			T * (  1.19977620e+01 + T * (  3.99529000e-03
			+ T * ( -2.69939067e-07 + T * ( -1.14218325e-10
			+ T *  1.72738220e-14 ) ) ) ) +  3.22849300e+04 );
		cp[sSXC4H5] =  1.56602501e+02 * (
			 1.19977620e+01 + T * (  7.99058000e-03 
			+ T * ( -8.09817200e-07 + T * ( -4.56873300e-10
			+ T *  8.63691100e-14 ) ) ) );
		h[sC4H6] =  1.53690340e+02 * (
			T * (  1.25443660e+01 + T * (  4.79826250e-03
			+ T * ( -3.06233733e-07 + T * ( -1.35741000e-10
			+ T *  2.01072720e-14 ) ) ) ) +  8.57873500e+03 );
		cp[sC4H6] =  1.53690340e+02 * (
			 1.25443660e+01 + T * (  9.59652500e-03 
			+ T * ( -9.18701200e-07 + T * ( -5.42964000e-10
			+ T *  1.00536360e-13 ) ) ) );
		h[sC4H7] =  1.50884509e+02 * (
			T * (  5.44066200e+00 + T * (  1.37230150e-02
			+ T * ( -4.55165333e-06 + T * (  8.35088000e-10
			+ T * -6.58010000e-14 ) ) ) ) +  2.09074600e+04 );
		cp[sC4H7] =  1.50884509e+02 * (
			 5.44066200e+00 + T * (  2.74460300e-02 
			+ T * ( -1.36549600e-05 + T * (  3.34035200e-09
			+ T * -3.29005000e-13 ) ) ) );
		h[s1XC4H8] =  1.48311452e+02 * (
			T * (  2.05358410e+00 + T * (  1.71752535e-02
			+ T * ( -5.29439867e-06 + T * (  8.27241550e-10
			+ T * -5.07220900e-14 ) ) ) ) -2.13972310e+03 );
		cp[s1XC4H8] =  1.48311452e+02 * (
			 2.05358410e+00 + T * (  3.43505070e-02 
			+ T * ( -1.58831960e-05 + T * (  3.30896620e-09
			+ T * -2.53610450e-13 ) ) ) );
		h[sPXC4H9] =  1.45579563e+02 * (
			T * (  9.00491050e+00 + T * (  1.20016595e-02
			+ T * ( -2.98506227e-06 + T * (  3.83064850e-10
			+ T * -1.96494080e-14 ) ) ) ) +  3.46301120e+03 );
		cp[sPXC4H9] =  1.45579563e+02 * (
			 9.00491050e+00 + T * (  2.40033190e-02 
			+ T * ( -8.95518680e-06 + T * (  1.53225940e-09
			+ T * -9.82470400e-14 ) ) ) );
		h[s1XC5H9] =  1.20285003e+02 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s1XC5H9] =  1.20285003e+02 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s1XC5H10] =  1.18556110e+02 * (
			T * (  1.25414540e+01 + T * (  1.19705350e-02
			+ T * ( -2.58151763e-06 + T * (  2.78949250e-10
			+ T * -1.15936232e-14 ) ) ) ) -8.71753520e+03 );
		cp[s1XC5H10] =  1.18556110e+02 * (
			 1.25414540e+01 + T * (  2.39410700e-02 
			+ T * ( -7.74455290e-06 + T * (  1.11579700e-09
			+ T * -5.79681160e-14 ) ) ) );
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
		h[s1XC6H13] =  9.76274012e+01 * (
			T * (  1.27597700e+01 + T * (  1.85671395e-02
			+ T * ( -4.75203667e-06 + T * (  6.26313900e-10
			+ T * -3.29194140e-14 ) ) ) ) -3.56937500e+03 );
		cp[s1XC6H13] =  9.76274012e+01 * (
			 1.27597700e+01 + T * (  3.71342790e-02 
			+ T * ( -1.42561100e-05 + T * (  2.50525560e-09
			+ T * -1.64597070e-13 ) ) ) );
		h[s1XC7H15] =  8.38223611e+01 * (
			T * (  1.49188600e+01 + T * (  2.16028655e-02
			+ T * ( -5.54310600e-06 + T * (  7.31940400e-10
			+ T * -3.85228520e-14 ) ) ) ) -7.17518750e+03 );
		cp[s1XC7H15] =  8.38223611e+01 * (
			 1.49188600e+01 + T * (  4.32057310e-02 
			+ T * ( -1.66293180e-05 + T * (  2.92776160e-09
			+ T * -1.92614260e-13 ) ) ) );
		h[s1XC7H14] =  8.46829358e+01 * (
			T * (  2.20515900e+01 + T * (  1.21791940e-02
			+ T * ( -1.18553937e-06 + T * ( -2.05597705e-10
			+ T *  3.78767540e-14 ) ) ) ) -1.76170040e+04 );
		cp[s1XC7H14] =  8.46829358e+01 * (
			 2.20515900e+01 + T * (  2.43583880e-02 
			+ T * ( -3.55661810e-06 + T * ( -8.22390820e-10
			+ T *  1.89383770e-13 ) ) ) );
		h[s1XC8H17] =  7.34378533e+01 * (
			T * (  1.70271150e+01 + T * (  2.46903135e-02
			+ T * ( -6.35649000e-06 + T * (  8.41807725e-10
			+ T * -4.44200780e-14 ) ) ) ) -1.07635270e+04 );
		cp[s1XC8H17] =  7.34378533e+01 * (
			 1.70271150e+01 + T * (  4.93806270e-02 
			+ T * ( -1.90694700e-05 + T * (  3.36723090e-09
			+ T * -2.22100390e-13 ) ) ) );
		h[s2XC10H21] =  5.88550840e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s2XC10H21] =  5.88550840e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s3XC10H21] =  5.88550840e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s3XC10H21] =  5.88550840e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s4XC10H21] =  5.88550840e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s4XC10H21] =  5.88550840e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s5XC10H21] =  5.88550840e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s5XC10H21] =  5.88550840e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s1XC10H21] =  5.88550840e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s1XC10H21] =  5.88550840e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sNXC10H22] =  5.84381062e+01 * (
			T * (  2.28188930e+01 + T * (  1.62717270e-02
			+ T * ( -3.70668033e-06 + T * (  4.28293575e-10
			+ T * -1.92424202e-14 ) ) ) ) -3.36787380e+04 );
		cp[sNXC10H22] =  5.84381062e+01 * (
			 2.28188930e+01 + T * (  3.25434540e-02 
			+ T * ( -1.11200410e-05 + T * (  1.71317430e-09
			+ T * -9.62121010e-14 ) ) ) );
		h[s2XRO2XC10H21O2] =  4.79854330e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s2XRO2XC10H21O2] =  4.79854330e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s3X2XRO2HXC10H20OOH] =  4.79854330e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s3X2XRO2HXC10H20OOH] =  4.79854330e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s32O2RO2HXO2C10H20OOH] =  4.05048035e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s32O2RO2HXO2C10H20OOH] =  4.05048035e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s3X2HOROHXHO2C10H19OOH] =  4.05048035e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s3X2HOROHXHO2C10H19OOH] =  4.05048035e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s3X2ORO2HXOC10H19OOH] =  4.41641347e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s3X2ORO2HXOC10H19OOH] =  4.41641347e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s2X1OROXOC10H19O] =  4.85503235e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s2X1OROXOC10H19O] =  4.85503235e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sUXC6H5] =  1.07838392e+02 * (
			T * (  1.72154000e+01 + T * (  4.31053400e-03
			+ T * ( -2.74044667e-07 + T * ( -1.18804100e-10
			+ T *  1.76881720e-14 ) ) ) ) +  6.40440300e+04 );
		cp[sUXC6H5] =  1.07838392e+02 * (
			 1.72154000e+01 + T * (  8.62106800e-03 
			+ T * ( -8.22134000e-07 + T * ( -4.75216400e-10
			+ T *  8.84408600e-14 ) ) ) );
		h[sC6H6] =  1.06446715e+02 * (
			T * (  9.77760890e+00 + T * (  1.29849920e-02
			+ T * ( -4.20744533e-06 + T * (  7.46770125e-10
			+ T * -5.58391360e-14 ) ) ) ) +  3.83029950e+04 );
		cp[sC6H6] =  1.06446715e+02 * (
			 9.77760890e+00 + T * (  2.59699840e-02 
			+ T * ( -1.26223360e-05 + T * (  2.98708050e-09
			+ T * -2.79195680e-13 ) ) ) );
		h[sUXC6H7] =  1.05090500e+02 * (
			T * (  2.18819190e+01 + T * (  3.02735450e-03
			+ T * (  8.79012400e-07 + T * ( -4.48246875e-10
			+ T *  4.85432300e-14 ) ) ) ) +  4.10291760e+04 );
		cp[sUXC6H7] =  1.05090500e+02 * (
			 2.18819190e+01 + T * (  6.05470900e-03 
			+ T * (  2.63703720e-06 + T * ( -1.79298750e-09
			+ T *  2.42716150e-13 ) ) ) );
		h[sA1XC6H6] =  1.06446715e+02 * (
			T * (  4.87752300e+00 + T * (  1.71098080e-02
			+ T * ( -5.89377267e-06 + T * (  1.09893858e-09
			+ T * -8.54067420e-14 ) ) ) ) +  6.42381650e+03 );
		cp[sA1XC6H6] =  1.06446715e+02 * (
			 4.87752300e+00 + T * (  3.42196160e-02 
			+ T * ( -1.76813180e-05 + T * (  4.39575430e-09
			+ T * -4.27033710e-13 ) ) ) );
		h[sA1XXC6H5] =  1.07838392e+02 * (
			T * (  5.21730470e+00 + T * (  1.53063115e-02
			+ T * ( -5.31589900e-06 + T * (  9.97514125e-10
			+ T * -7.79098880e-14 ) ) ) ) +  3.58955480e+04 );
		cp[sA1XXC6H5] =  1.07838392e+02 * (
			 5.21730470e+00 + T * (  3.06126230e-02 
			+ T * ( -1.59476970e-05 + T * (  3.99005650e-09
			+ T * -3.89549440e-13 ) ) ) );
		h[sC6H5O] =  8.93054780e+01 * (
			T * (  1.38339840e+01 + T * (  8.80920150e-03
			+ T * ( -2.02320857e-06 + T * (  2.29970432e-10
			+ T * -1.00898362e-14 ) ) ) ) -6.92125490e+02 );
		cp[sC6H5O] =  8.93054780e+01 * (
			 1.38339840e+01 + T * (  1.76184030e-02 
			+ T * ( -6.06962570e-06 + T * (  9.19881730e-10
			+ T * -5.04491810e-14 ) ) ) );
		h[sC6H5OH] =  8.83489183e+01 * (
			T * (  1.49435680e+01 + T * (  9.15892250e-03
			+ T * ( -2.05354930e-06 + T * (  2.27621310e-10
			+ T * -9.68216120e-15 ) ) ) ) -1.83868280e+04 );
		cp[sC6H5OH] =  8.83489183e+01 * (
			 1.49435680e+01 + T * (  1.83178450e-02 
			+ T * ( -6.16064790e-06 + T * (  9.10485240e-10
			+ T * -4.84108060e-14 ) ) ) );
		h[sC5H5] =  1.27736058e+02 * (
			T * (  9.68981460e+00 + T * (  9.19131000e-03
			+ T * ( -2.08829473e-06 + T * (  2.34834430e-10
			+ T * -1.01754162e-14 ) ) ) ) +  1.10212420e+04 );
		cp[sC5H5] =  1.27736058e+02 * (
			 9.68981460e+00 + T * (  1.83826200e-02 
			+ T * ( -6.26488420e-06 + T * (  9.39337720e-10
			+ T * -5.08770810e-14 ) ) ) );
		h[sC5H6] =  1.25788072e+02 * (
			T * (  9.68981460e+00 + T * (  9.19131000e-03
			+ T * ( -2.08829473e-06 + T * (  2.34834430e-10
			+ T * -1.01754162e-14 ) ) ) ) +  1.10212420e+04 );
		cp[sC5H6] =  1.25788072e+02 * (
			 9.68981460e+00 + T * (  1.83826200e-02 
			+ T * ( -6.26488420e-06 + T * (  9.39337720e-10
			+ T * -5.08770810e-14 ) ) ) );
		h[sC5H5O] =  1.02532248e+02 * (
			T * (  1.38339840e+01 + T * (  8.80920150e-03
			+ T * ( -2.02320857e-06 + T * (  2.29970432e-10
			+ T * -1.00898362e-14 ) ) ) ) -6.92125490e+02 );
		cp[sC5H5O] =  1.02532248e+02 * (
			 1.38339840e+01 + T * (  1.76184030e-02 
			+ T * ( -6.06962570e-06 + T * (  9.19881730e-10
			+ T * -5.04491810e-14 ) ) ) );
		h[sA1C2HXC8H6] =  8.14109745e+01 * (
			T * (  9.09223910e+00 + T * (  1.80699605e-02
			+ T * ( -6.20154033e-06 + T * (  1.15291007e-09
			+ T * -8.93769520e-14 ) ) ) ) +  3.28858950e+04 );
		cp[sA1C2HXC8H6] =  8.14109745e+01 * (
			 9.09223910e+00 + T * (  3.61399210e-02 
			+ T * ( -1.86046210e-05 + T * (  4.61164030e-09
			+ T * -4.46884760e-13 ) ) ) );
		h[sA1C2HYXC8H5] =  8.22225079e+01 * (
			T * (  9.43202080e+00 + T * (  1.62664645e-02
			+ T * ( -5.62366633e-06 + T * (  1.05148562e-09
			+ T * -8.18800960e-14 ) ) ) ) +  6.23576270e+04 );
		cp[sA1C2HYXC8H5] =  8.22225079e+01 * (
			 9.43202080e+00 + T * (  3.25329290e-02 
			+ T * ( -1.68709990e-05 + T * (  4.20594250e-09
			+ T * -4.09400480e-13 ) ) ) );
		h[sA1C2HACXC10H7] =  6.53869263e+01 * (
			T * (  1.21998720e+01 + T * (  2.12108275e-02
			+ T * ( -7.27111867e-06 + T * (  1.34979190e-09
			+ T * -1.04474232e-13 ) ) ) ) +  6.75357020e+04 );
		cp[sA1C2HACXC10H7] =  6.53869263e+01 * (
			 1.21998720e+01 + T * (  4.24216550e-02 
			+ T * ( -2.18133560e-05 + T * (  5.39916760e-09
			+ T * -5.22371160e-13 ) ) ) );
		h[sA2XXXC10H7] =  6.53869263e+01 * (
			T * (  9.14717570e+00 + T * (  2.40291180e-02
			+ T * ( -8.43824033e-06 + T * (  1.59573013e-09
			+ T * -1.25301978e-13 ) ) ) ) +  4.17725310e+04 );
		cp[sA2XXXC10H7] =  6.53869263e+01 * (
			 9.14717570e+00 + T * (  4.80582360e-02 
			+ T * ( -2.53147210e-05 + T * (  6.38292050e-09
			+ T * -6.26509890e-13 ) ) ) );
		h[sA2XC10H8] =  6.48726632e+01 * (
			T * (  8.80739400e+00 + T * (  2.58326140e-02
			+ T * ( -9.01611433e-06 + T * (  1.69715458e-09
			+ T * -1.32798832e-13 ) ) ) ) +  1.23007990e+04 );
		cp[sA2XC10H8] =  6.48726632e+01 * (
			 8.80739400e+00 + T * (  5.16652280e-02 
			+ T * ( -2.70483430e-05 + T * (  6.78861830e-09
			+ T * -6.63994160e-13 ) ) ) );
		h[sA2R5XC12H8] =  5.46334700e+01 * (
			T * (  1.16241070e+01 + T * (  2.79811160e-02
			+ T * ( -9.80209567e-06 + T * (  1.84608017e-09
			+ T * -1.44315888e-13 ) ) ) ) +  2.40156570e+04 );
		cp[sA2R5XC12H8] =  5.46334700e+01 * (
			 1.16241070e+01 + T * (  5.59622320e-02 
			+ T * ( -2.94062870e-05 + T * (  7.38432070e-09
			+ T * -7.21579440e-13 ) ) ) );
		h[sA2R5XXC12H7] =  5.49977510e+01 * (
			T * (  1.19638890e+01 + T * (  2.61776200e-02
			+ T * ( -9.22422167e-06 + T * (  1.74465572e-09
			+ T * -1.36819034e-13 ) ) ) ) +  5.34873890e+04 );
		cp[sA2R5XXC12H7] =  5.49977510e+01 * (
			 1.19638890e+01 + T * (  5.23552400e-02 
			+ T * ( -2.76726650e-05 + T * (  6.97862290e-09
			+ T * -6.84095170e-13 ) ) ) );
		h[sA2R5C2HXC14H8] =  4.71858755e+01 * (
			T * (  1.58388230e+01 + T * (  2.89412685e-02
			+ T * ( -1.01098633e-05 + T * (  1.90005167e-09
			+ T * -1.48286096e-13 ) ) ) ) +  5.04777350e+04 );
		cp[sA2R5C2HXC14H8] =  4.71858755e+01 * (
			 1.58388230e+01 + T * (  5.78825370e-02 
			+ T * ( -3.03295900e-05 + T * (  7.60020670e-09
			+ T * -7.41430480e-13 ) ) ) );
		h[sA2R5C2HYXC14H7] =  4.74573620e+01 * (
			T * (  1.61786050e+01 + T * (  2.71377725e-02
			+ T * ( -9.53198933e-06 + T * (  1.79862722e-09
			+ T * -1.40789242e-13 ) ) ) ) +  7.99494670e+04 );
		cp[sA2R5C2HYXC14H7] =  4.74573620e+01 * (
			 1.61786050e+01 + T * (  5.42755450e-02 
			+ T * ( -2.85959680e-05 + T * (  7.19450890e-09
			+ T * -7.03946210e-13 ) ) ) );
		h[sANC2HACXC16H9] =  4.13171861e+01 * (
			T * (  1.89464560e+01 + T * (  3.20821360e-02
			+ T * ( -1.11794417e-05 + T * (  2.09693350e-09
			+ T * -1.63383376e-13 ) ) ) ) +  8.51275420e+04 );
		cp[sANC2HACXC16H9] =  4.13171861e+01 * (
			 1.89464560e+01 + T * (  6.41642720e-02 
			+ T * ( -3.35383250e-05 + T * (  8.38773400e-09
			+ T * -8.16916880e-13 ) ) ) );
		h[sA3R5XXC16H9] =  4.13171861e+01 * (
			T * (  1.58937590e+01 + T * (  3.49004260e-02
			+ T * ( -1.23465633e-05 + T * (  2.34287170e-09
			+ T * -1.84211124e-13 ) ) ) ) +  5.82573190e+04 );
		cp[sA3R5XXC16H9] =  4.13171861e+01 * (
			 1.58937590e+01 + T * (  6.98008520e-02 
			+ T * ( -3.70396900e-05 + T * (  9.37148680e-09
			+ T * -9.21055620e-13 ) ) ) );
		h[sA3R5XC16H10] =  4.11112540e+01 * (
			T * (  1.55539780e+01 + T * (  3.67039220e-02
			+ T * ( -1.29244370e-05 + T * (  2.44429615e-09
			+ T * -1.91707978e-13 ) ) ) ) +  2.87855880e+04 );
		cp[sA3R5XC16H10] =  4.11112540e+01 * (
			 1.55539780e+01 + T * (  7.34078440e-02 
			+ T * ( -3.87733110e-05 + T * (  9.77718460e-09
			+ T * -9.58539890e-13 ) ) ) );
		h[sA3R5ACXC18H11] =  3.65838569e+01 * (
			T * (  1.86616100e+01 + T * (  3.98447890e-02
			+ T * ( -1.39940153e-05 + T * (  2.64117800e-09
			+ T * -2.06805260e-13 ) ) ) ) +  6.34353960e+04 );
		cp[sA3R5ACXC18H11] =  3.65838569e+01 * (
			 1.86616100e+01 + T * (  7.96895780e-02 
			+ T * ( -4.19820460e-05 + T * (  1.05647120e-08
			+ T * -1.03402630e-12 ) ) ) );
		h[sA4XC18H10] =  3.67468399e+01 * (
			T * (  1.86348930e+01 + T * (  3.88005750e-02
			+ T * ( -1.37573170e-05 + T * (  2.61631575e-09
			+ T * -2.06118560e-13 ) ) ) ) +  2.99008240e+04 );
		cp[sA4XC18H10] =  3.67468399e+01 * (
			 1.86348930e+01 + T * (  7.76011500e-02 
			+ T * ( -4.12719510e-05 + T * (  1.04652630e-08
			+ T * -1.03059280e-12 ) ) ) );
		h[sA2CH3XC11H10] =  5.84734510e+01 * (
			T * (  2.55377270e+01 + T * (  7.88844950e-03
			+ T * (  2.11806143e-06 + T * (  7.07124350e-11
			+ T * -4.75411680e-13 ) ) ) ) +  2.80964450e+03 );
		cp[sA2CH3XC11H10] =  5.84734510e+01 * (
			 2.55377270e+01 + T * (  1.57768990e-02 
			+ T * (  6.35418430e-06 + T * (  2.82849740e-10
			+ T * -2.37705840e-12 ) ) ) );
		h[sA2CH2XC11H9] =  5.88909351e+01 * (
			T * (  2.53690800e+01 + T * (  7.56497750e-03
			+ T * (  1.96068580e-06 + T * (  4.39942775e-11
			+ T * -4.44976440e-13 ) ) ) ) +  1.99007250e+04 );
		cp[sA2CH2XC11H9] =  5.88909351e+01 * (
			 2.53690800e+01 + T * (  1.51299550e-02 
			+ T * (  5.88205740e-06 + T * (  1.75977110e-10
			+ T * -2.22488220e-12 ) ) ) );
		h[sA2CHOXC11H8O] =  5.32376708e+01 * (
			T * (  1.62272000e+01 + T * (  6.20089500e-03
			+ T * (  2.07737123e-06 + T * (  3.03401025e-10
			+ T * -2.90436360e-13 ) ) ) ) -7.73170090e+01 );
		cp[sA2CHOXC11H8O] =  5.32376708e+01 * (
			 1.62272000e+01 + T * (  1.24017900e-02 
			+ T * (  6.23211370e-06 + T * (  1.21360410e-09
			+ T * -1.45218180e-12 ) ) ) );
		h[sBIBENZYLXC22H18] =  2.94454676e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[sBIBENZYLXC22H18] =  2.94454676e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
	}
	else if (T >= 299.999999 ) {
		h[sN2] =  2.96792318e+02 * (
			T * (  3.29867700e+00 + T * (  7.04120200e-04
			+ T * ( -1.32107400e-06 + T * (  1.41037875e-09
			+ T * -4.88970800e-13 ) ) ) ) -1.02089990e+03 );
		cp[sN2] =  2.96792318e+02 * (
			 3.29867700e+00 + T * (  1.40824040e-03 
			+ T * ( -3.96322200e-06 + T * (  5.64151500e-09
			+ T * -2.44485400e-12 ) ) ) );
		h[sO2] =  2.59823125e+02 * (
			T * (  3.21293600e+00 + T * (  5.63743200e-04
			+ T * ( -1.91871667e-07 + T * (  3.28469325e-10
			+ T * -1.75371080e-13 ) ) ) ) -1.00524900e+03 );
		cp[sO2] =  2.59823125e+02 * (
			 3.21293600e+00 + T * (  1.12748640e-03 
			+ T * ( -5.75615000e-07 + T * (  1.31387730e-09
			+ T * -8.76855400e-13 ) ) ) );
		h[sH] =  8.24835317e+03 * (
			T * (  2.50000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  2.54716200e+04 );
		cp[sH] =  8.24835317e+03 * (
			 2.50000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sOH] =  4.88848777e+02 * (
			T * (  3.63726600e+00 + T * (  9.25455000e-05
			+ T * ( -5.58721533e-07 + T * (  5.96800500e-10
			+ T * -1.68628840e-13 ) ) ) ) +  3.60678100e+03 );
		cp[sOH] =  4.88848777e+02 * (
			 3.63726600e+00 + T * (  1.85091000e-04 
			+ T * ( -1.67616460e-06 + T * (  2.38720200e-09
			+ T * -8.43144200e-13 ) ) ) );
		h[sO] =  5.19646250e+02 * (
			T * (  2.94642800e+00 + T * ( -8.19083250e-04
			+ T * (  8.07010333e-07 + T * ( -4.00710775e-10
			+ T *  7.78139200e-14 ) ) ) ) +  2.91476400e+04 );
		cp[sO] =  5.19646250e+02 * (
			 2.94642800e+00 + T * ( -1.63816650e-03 
			+ T * (  2.42103100e-06 + T * ( -1.60284310e-09
			+ T *  3.89069600e-13 ) ) ) );
		h[sH2] =  4.12417659e+03 * (
			T * (  3.29812400e+00 + T * (  4.12472050e-04
			+ T * ( -2.71433833e-07 + T * ( -2.36885850e-11
			+ T *  8.26974400e-14 ) ) ) ) -1.01252090e+03 );
		cp[sH2] =  4.12417659e+03 * (
			 3.29812400e+00 + T * (  8.24944100e-04 
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
			T * (  2.97996300e+00 + T * (  2.49834850e-03
			+ T * ( -1.26366567e-06 + T * (  5.88548000e-10
			+ T * -1.61780480e-13 ) ) ) ) +  1.76227300e+02 );
		cp[sHO2] =  2.51888633e+02 * (
			 2.97996300e+00 + T * (  4.99669700e-03 
			+ T * ( -3.79099700e-06 + T * (  2.35419200e-09
			+ T * -8.08902400e-13 ) ) ) );
		h[sH2O2] =  2.44424389e+02 * (
			T * (  3.38875300e+00 + T * (  3.28461300e-03
			+ T * ( -4.95004167e-08 + T * ( -1.15645125e-09
			+ T *  4.94302800e-13 ) ) ) ) -1.76631400e+04 );
		cp[sH2O2] =  2.44424389e+02 * (
			 3.38875300e+00 + T * (  6.56922600e-03 
			+ T * ( -1.48501250e-07 + T * ( -4.62580500e-09
			+ T *  2.47151400e-12 ) ) ) );
		h[sCO] =  2.96824105e+02 * (
			T * (  3.26245100e+00 + T * (  7.55970450e-04
			+ T * ( -1.29391833e-06 + T * (  1.39548600e-09
			+ T * -4.94990200e-13 ) ) ) ) -1.43105390e+04 );
		cp[sCO] =  2.96824105e+02 * (
			 3.26245100e+00 + T * (  1.51194090e-03 
			+ T * ( -3.88175500e-06 + T * (  5.58194400e-09
			+ T * -2.47495100e-12 ) ) ) );
		h[sCO2] =  1.88915044e+02 * (
			T * (  2.27572400e+00 + T * (  4.96103600e-03
			+ T * ( -3.46970433e-06 + T * (  1.71667150e-09
			+ T * -4.23456000e-13 ) ) ) ) -4.83731400e+04 );
		cp[sCO2] =  1.88915044e+02 * (
			 2.27572400e+00 + T * (  9.92207200e-03 
			+ T * ( -1.04091130e-05 + T * (  6.86668600e-09
			+ T * -2.11728000e-12 ) ) ) );
		h[sCH] =  6.38631231e+02 * (
			T * (  3.20020200e+00 + T * (  1.03643750e-03
			+ T * ( -1.71147700e-06 + T * (  1.43347250e-09
			+ T * -3.91106600e-13 ) ) ) ) +  7.04525900e+04 );
		cp[sCH] =  6.38631231e+02 * (
			 3.20020200e+00 + T * (  2.07287500e-03 
			+ T * ( -5.13443100e-06 + T * (  5.73389000e-09
			+ T * -1.95553300e-12 ) ) ) );
		h[sHCO] =  2.86513663e+02 * (
			T * (  2.89832900e+00 + T * (  3.09957300e-03
			+ T * ( -3.20769467e-06 + T * (  2.72456225e-09
			+ T * -9.14977000e-13 ) ) ) ) +  4.15992200e+03 );
		cp[sHCO] =  2.86513663e+02 * (
			 2.89832900e+00 + T * (  6.19914600e-03 
			+ T * ( -9.62308400e-06 + T * (  1.08982490e-08
			+ T * -4.57488500e-12 ) ) ) );
		h[sCH2OH] =  2.67971122e+02 * (
			T * (  2.86262800e+00 + T * (  5.00763650e-03
			+ T * ( -1.76181167e-07 + T * ( -1.28463475e-09
			+ T *  4.49208200e-13 ) ) ) ) -3.34967800e+03 );
		cp[sCH2OH] =  2.67971122e+02 * (
			 2.86262800e+00 + T * (  1.00152730e-02 
			+ T * ( -5.28543500e-07 + T * ( -5.13853900e-09
			+ T *  2.24604100e-12 ) ) ) );
		h[s3XCH2] =  5.92738290e+02 * (
			T * (  3.76223700e+00 + T * (  5.79909550e-04
			+ T * (  8.29861667e-08 + T * (  2.20020900e-10
			+ T * -1.46648700e-13 ) ) ) ) +  4.53679000e+04 );
		cp[s3XCH2] =  5.92738290e+02 * (
			 3.76223700e+00 + T * (  1.15981910e-03 
			+ T * (  2.48958500e-07 + T * (  8.80083600e-10
			+ T * -7.33243500e-13 ) ) ) );
		h[sC2H2] =  3.19217538e+02 * (
			T * (  2.01356200e+00 + T * (  7.59522300e-03
			+ T * ( -5.38772967e-06 + T * (  2.26974800e-09
			+ T * -3.82549200e-13 ) ) ) ) +  2.61244400e+04 );
		cp[sC2H2] =  3.19217538e+02 * (
			 2.01356200e+00 + T * (  1.51904460e-02 
			+ T * ( -1.61631890e-05 + T * (  9.07899200e-09
			+ T * -1.91274600e-12 ) ) ) );
		h[sCH3] =  5.52999002e+02 * (
			T * (  2.43044200e+00 + T * (  5.56204950e-03
			+ T * ( -5.60073333e-06 + T * (  4.05457200e-09
			+ T * -1.17299040e-12 ) ) ) ) +  1.64237810e+04 );
		cp[sCH3] =  5.52999002e+02 * (
			 2.43044200e+00 + T * (  1.11240990e-02 
			+ T * ( -1.68022000e-05 + T * (  1.62182880e-08
			+ T * -5.86495200e-12 ) ) ) );
		h[sC2H4] =  2.96284655e+02 * (
			T * ( -8.61488000e-01 + T * (  1.39808100e-02
			+ T * ( -1.12955900e-05 + T * (  6.96288000e-09
			+ T * -1.94757580e-12 ) ) ) ) +  5.57304600e+03 );
		cp[sC2H4] =  2.96284655e+02 * (
			-8.61488000e-01 + T * (  2.79616200e-02 
			+ T * ( -3.38867700e-05 + T * (  2.78515200e-08
			+ T * -9.73787900e-12 ) ) ) );
		h[s1XCH2] =  5.92738290e+02 * (
			T * (  3.97126500e+00 + T * ( -8.49544000e-05
			+ T * (  3.41789633e-07 + T * (  6.23137500e-10
			+ T * -3.96253200e-13 ) ) ) ) +  4.98936700e+04 );
		cp[s1XCH2] =  5.92738290e+02 * (
			 3.97126500e+00 + T * ( -1.69908800e-04 
			+ T * (  1.02536890e-06 + T * (  2.49255000e-09
			+ T * -1.98126600e-12 ) ) ) );
		h[sCH2O] =  2.76895461e+02 * (
			T * (  1.65273110e+00 + T * (  6.31571950e-03
			+ T * ( -6.29389333e-06 + T * (  5.12507750e-09
			+ T * -1.68264740e-12 ) ) ) ) -1.48654040e+04 );
		cp[sCH2O] =  2.76895461e+02 * (
			 1.65273110e+00 + T * (  1.26314390e-02 
			+ T * ( -1.88816800e-05 + T * (  2.05003100e-08
			+ T * -8.41323700e-12 ) ) ) );
		h[sCH3O] =  2.67910679e+02 * (
			T * (  2.10620400e+00 + T * (  3.60829750e-03
			+ T * (  1.77949067e-06 + T * ( -1.84440900e-09
			+ T *  4.15122000e-13 ) ) ) ) +  9.78601100e+02 );
		cp[sCH3O] =  2.67910679e+02 * (
			 2.10620400e+00 + T * (  7.21659500e-03 
			+ T * (  5.33847200e-06 + T * ( -7.37763600e-09
			+ T *  2.07561000e-12 ) ) ) );
		h[sCH4] =  5.18253444e+02 * (
			T * (  7.78741500e-01 + T * (  8.73834000e-03
			+ T * ( -9.27803000e-06 + T * (  7.62427000e-09
			+ T * -2.44786140e-12 ) ) ) ) -9.82522900e+03 );
		cp[sCH4] =  5.18253444e+02 * (
			 7.78741500e-01 + T * (  1.74766800e-02 
			+ T * ( -2.78340900e-05 + T * (  3.04970800e-08
			+ T * -1.22393070e-11 ) ) ) );
		h[sC2H6] =  2.76425959e+02 * (
			T * (  1.46253880e+00 + T * (  7.74733350e-03
			+ T * (  1.92683567e-06 + T * ( -3.14457975e-09
			+ T *  9.17253400e-13 ) ) ) ) -1.12391760e+04 );
		cp[sC2H6] =  2.76425959e+02 * (
			 1.46253880e+00 + T * (  1.54946670e-02 
			+ T * (  5.78050700e-06 + T * ( -1.25783190e-08
			+ T *  4.58626700e-12 ) ) ) );
		h[sC2H] =  3.32068855e+02 * (
			T * (  3.05066700e+00 + T * (  3.02583700e-03
			+ T * ( -1.65221133e-06 + T * (  7.01039750e-10
			+ T * -1.63866640e-13 ) ) ) ) +  6.63001100e+04 );
		cp[sC2H] =  3.32068855e+02 * (
			 3.05066700e+00 + T * (  6.05167400e-03 
			+ T * ( -4.95663400e-06 + T * (  2.80415900e-09
			+ T * -8.19333200e-13 ) ) ) );
		h[sHCCO] =  2.02640507e+02 * (
			T * (  5.04796500e+00 + T * (  2.22673900e-03
			+ T * (  7.56094000e-08 + T * ( -3.70523625e-10
			+ T *  4.50148200e-14 ) ) ) ) +  1.96589100e+04 );
		cp[sHCCO] =  2.02640507e+02 * (
			 5.04796500e+00 + T * (  4.45347800e-03 
			+ T * (  2.26828200e-07 + T * ( -1.48209450e-09
			+ T *  2.25074100e-13 ) ) ) );
		h[sCH2CO] =  1.97781531e+02 * (
			T * (  2.97497000e+00 + T * (  6.05935600e-03
			+ T * ( -7.81681667e-07 + T * ( -1.61667125e-09
			+ T *  7.81129800e-13 ) ) ) ) -7.63263600e+03 );
		cp[sCH2CO] =  1.97781531e+02 * (
			 2.97497000e+00 + T * (  1.21187120e-02 
			+ T * ( -2.34504500e-06 + T * ( -6.46668500e-09
			+ T *  3.90564900e-12 ) ) ) );
		h[sC2H3] =  3.07323871e+02 * (
			T * (  2.45927600e+00 + T * (  3.68573800e-03
			+ T * (  7.03290667e-07 + T * ( -3.30410525e-10
			+ T * -2.36956760e-13 ) ) ) ) +  3.33522500e+04 );
		cp[sC2H3] =  3.07323871e+02 * (
			 2.45927600e+00 + T * (  7.37147600e-03 
			+ T * (  2.10987200e-06 + T * ( -1.32164210e-09
			+ T * -1.18478380e-12 ) ) ) );
		h[sC2H5] =  2.86011008e+02 * (
			T * (  2.69070100e+00 + T * (  4.35956650e-03
			+ T * (  1.47327933e-06 + T * (  2.33467575e-10
			+ T * -7.85554600e-13 ) ) ) ) +  1.28704040e+04 );
		cp[sC2H5] =  2.86011008e+02 * (
			 2.69070100e+00 + T * (  8.71913300e-03 
			+ T * (  4.41983800e-06 + T * (  9.33870300e-10
			+ T * -3.92777300e-12 ) ) ) );
		h[sC3H3] =  2.12833483e+02 * (
			T * (  4.75419900e+00 + T * (  5.54013850e-03
			+ T * (  9.31107667e-08 + T * ( -1.36980300e-09
			+ T *  3.89925800e-13 ) ) ) ) +  3.90468300e+04 );
		cp[sC3H3] =  2.12833483e+02 * (
			 4.75419900e+00 + T * (  1.10802770e-02 
			+ T * (  2.79332300e-07 + T * ( -5.47921200e-09
			+ T *  1.94962900e-12 ) ) ) );
		h[sC3H4] =  2.07479849e+02 * (
			T * ( -2.13196800e-01 + T * (  1.67935650e-02
			+ T * ( -1.26829000e-05 + T * (  6.86459500e-09
			+ T * -1.73800880e-12 ) ) ) ) +  2.16204800e+04 );
		cp[sC3H4] =  2.07479849e+02 * (
			-2.13196800e-01 + T * (  3.35871300e-02 
			+ T * ( -3.80487000e-05 + T * (  2.74583800e-08
			+ T * -8.69004400e-12 ) ) ) );
		h[sPXC3H4] =  2.07479849e+02 * (
			T * (  6.27144700e-01 + T * (  1.55808950e-02
			+ T * ( -1.24922100e-05 + T * (  7.41029250e-09
			+ T * -1.99747620e-12 ) ) ) ) +  2.11356000e+04 );
		cp[sPXC3H4] =  2.07479849e+02 * (
			 6.27144700e-01 + T * (  3.11617900e-02 
			+ T * ( -3.74766300e-05 + T * (  2.96411700e-08
			+ T * -9.98738100e-12 ) ) ) );
		h[sC3H5] =  2.02388939e+02 * (
			T * (  2.27648600e+00 + T * (  9.92782000e-03
			+ T * (  3.74614033e-07 + T * ( -2.53643925e-09
			+ T *  6.88268400e-13 ) ) ) ) +  1.78949600e+04 );
		cp[sC3H5] =  2.02388939e+02 * (
			 2.27648600e+00 + T * (  1.98556400e-02 
			+ T * (  1.12384210e-06 + T * ( -1.01457570e-08
			+ T *  3.44134200e-12 ) ) ) );
		h[sC3H6] =  1.97725089e+02 * (
			T * (  1.49330710e+00 + T * (  1.04625850e-02
			+ T * (  1.49559800e-06 + T * ( -4.17228025e-09
			+ T *  1.43162920e-12 ) ) ) ) +  1.07482640e+03 );
		cp[sC3H6] =  1.97725089e+02 * (
			 1.49330710e+00 + T * (  2.09251700e-02 
			+ T * (  4.48679400e-06 + T * ( -1.66891210e-08
			+ T *  7.15814600e-12 ) ) ) );
		h[sNXC3H7] =  1.93087320e+02 * (
			T * (  1.92253600e+00 + T * (  1.23946350e-02
			+ T * (  6.03416333e-07 + T * ( -4.45816250e-09
			+ T *  1.71659920e-12 ) ) ) ) +  1.03988710e+04 );
		cp[sNXC3H7] =  1.93087320e+02 * (
			 1.92253600e+00 + T * (  2.47892700e-02 
			+ T * (  1.81024900e-06 + T * ( -1.78326500e-08
			+ T *  8.58299600e-12 ) ) ) );
		h[sC4H2] =  1.66060957e+02 * (
			T * (  4.00519180e+00 + T * (  9.90500100e-03
			+ T * ( -3.28862583e-06 + T * ( -1.65878955e-09
			+ T *  1.21548258e-12 ) ) ) ) +  5.11566930e+04 );
		cp[sC4H2] =  1.66060957e+02 * (
			 4.00519180e+00 + T * (  1.98100020e-02 
			+ T * ( -9.86587750e-06 + T * ( -6.63515820e-09
			+ T *  6.07741290e-12 ) ) ) );
		h[sUXC4H3] =  1.62821949e+02 * (
			T * (  4.15388100e+00 + T * (  8.63143500e-03
			+ T * ( -7.96458000e-08 + T * ( -2.54675000e-09
			+ T *  8.68100800e-13 ) ) ) ) +  6.10352600e+04 );
		cp[sUXC4H3] =  1.62821949e+02 * (
			 4.15388100e+00 + T * (  1.72628700e-02 
			+ T * ( -2.38937400e-07 + T * ( -1.01870000e-08
			+ T *  4.34050400e-12 ) ) ) );
		h[sSXC4H3] =  1.62783695e+02 * (
			T * (  6.54579900e+00 + T * (  6.21238400e-03
			+ T * (  1.86774200e-07 + T * ( -1.40778525e-09
			+ T *  3.33043660e-13 ) ) ) ) +  5.35250200e+04 );
		cp[sSXC4H3] =  1.62783695e+02 * (
			 6.54579900e+00 + T * (  1.24247680e-02 
			+ T * (  5.60322600e-07 + T * ( -5.63114100e-09
			+ T *  1.66521830e-12 ) ) ) );
		h[sC4H4] =  1.59633285e+02 * (
			T * (  3.23389300e+00 + T * (  9.32817000e-03
			+ T * (  4.23440167e-07 + T * ( -2.35252400e-09
			+ T *  5.91222000e-13 ) ) ) ) +  3.30109700e+04 );
		cp[sC4H4] =  1.59633285e+02 * (
			 3.23389300e+00 + T * (  1.86563400e-02 
			+ T * (  1.27032050e-06 + T * ( -9.41009600e-09
			+ T *  2.95611000e-12 ) ) ) );
		h[sUXC4H5] =  1.56637905e+02 * (
			T * (  2.99524000e+00 + T * (  1.14422800e-02
			+ T * (  6.58490333e-07 + T * ( -2.87061350e-09
			+ T *  6.39564600e-13 ) ) ) ) +  4.04413400e+04 );
		cp[sUXC4H5] =  1.56637905e+02 * (
			 2.99524000e+00 + T * (  2.28845600e-02 
			+ T * (  1.97547100e-06 + T * ( -1.14824540e-08
			+ T *  3.19782300e-12 ) ) ) );
		h[sSXC4H5] =  1.56602501e+02 * (
			T * (  3.87944300e+00 + T * (  9.98831500e-03
			+ T * (  6.24259000e-07 + T * ( -2.32673825e-09
			+ T *  4.77223200e-13 ) ) ) ) +  3.52685900e+04 );
		cp[sSXC4H5] =  1.56602501e+02 * (
			 3.87944300e+00 + T * (  1.99766300e-02 
			+ T * (  1.87277700e-06 + T * ( -9.30695300e-09
			+ T *  2.38611600e-12 ) ) ) );
		h[sC4H6] =  1.53690340e+02 * (
			T * (  1.93162400e+00 + T * (  1.23951500e-02
			+ T * (  1.00602367e-06 + T * ( -2.88671400e-09
			+ T *  5.17324600e-13 ) ) ) ) +  1.25360850e+04 );
		cp[sC4H6] =  1.53690340e+02 * (
			 1.93162400e+00 + T * (  2.47903000e-02 
			+ T * (  3.01807100e-06 + T * ( -1.15468560e-08
			+ T *  2.58662300e-12 ) ) ) );
		h[sC4H7] =  1.50884509e+02 * (
			T * (  9.27815100e+00 + T * ( -1.47815250e-02
			+ T * (  5.44288667e-05 + T * ( -5.00144000e-08
			+ T *  1.58692100e-11 ) ) ) ) +  2.15084700e+04 );
		cp[sC4H7] =  1.50884509e+02 * (
			 9.27815100e+00 + T * ( -2.95630500e-02 
			+ T * (  1.63286600e-04 + T * ( -2.00057600e-07
			+ T *  7.93460500e-11 ) ) ) );
		h[s1XC4H8] =  1.48311452e+02 * (
			T * (  1.18113800e+00 + T * (  1.54266900e-02
			+ T * (  1.69550823e-06 + T * ( -6.16372200e-09
			+ T *  2.22203860e-12 ) ) ) ) -1.79040040e+03 );
		cp[s1XC4H8] =  1.48311452e+02 * (
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
		h[s1XC5H9] =  1.20285003e+02 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s1XC5H9] =  1.20285003e+02 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s1XC5H10] =  1.18556110e+02 * (
			T * (  1.76278590e+00 + T * (  1.94513945e-02
			+ T * (  1.94602580e-06 + T * ( -7.63679700e-09
			+ T *  2.76789300e-12 ) ) ) ) -4.81098830e+03 );
		cp[s1XC5H10] =  1.18556110e+02 * (
			 1.76278590e+00 + T * (  3.89027890e-02 
			+ T * (  5.83807740e-06 + T * ( -3.05471880e-08
			+ T *  1.38394650e-11 ) ) ) );
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
		h[s1XC6H13] =  9.76274012e+01 * (
			T * (  2.30623530e+00 + T * (  2.51000225e-02
			+ T * (  2.64625910e-06 + T * ( -1.03151638e-08
			+ T *  3.75887380e-12 ) ) ) ) +  1.02697110e+02 );
		cp[s1XC6H13] =  9.76274012e+01 * (
			 2.30623530e+00 + T * (  5.02000450e-02 
			+ T * (  7.93877730e-06 + T * ( -4.12606550e-08
			+ T *  1.87943690e-11 ) ) ) );
		h[s1XC7H15] =  8.38223611e+01 * (
			T * (  2.46798900e+00 + T * (  2.94440055e-02
			+ T * (  3.16893967e-06 + T * ( -1.21711593e-08
			+ T *  4.41048580e-12 ) ) ) ) -2.82337940e+03 );
		cp[s1XC7H15] =  8.38223611e+01 * (
			 2.46798900e+00 + T * (  5.88880110e-02 
			+ T * (  9.50681900e-06 + T * ( -4.86846370e-08
			+ T *  2.20524290e-11 ) ) ) );
		h[s1XC7H14] =  8.46829358e+01 * (
			T * (  2.29726600e+00 + T * (  2.79542730e-02
			+ T * (  2.71505207e-06 + T * ( -1.11203588e-08
			+ T *  4.07020260e-12 ) ) ) ) -1.06558790e+04 );
		cp[s1XC7H14] =  8.46829358e+01 * (
			 2.29726600e+00 + T * (  5.59085460e-02 
			+ T * (  8.14515620e-06 + T * ( -4.44814350e-08
			+ T *  2.03510130e-11 ) ) ) );
		h[s1XC8H17] =  7.34378533e+01 * (
			T * (  2.61413190e+00 + T * (  3.38039100e-02
			+ T * (  3.73002267e-06 + T * ( -1.41030850e-08
			+ T *  5.09678420e-12 ) ) ) ) -5.74694140e+03 );
		cp[s1XC8H17] =  7.34378533e+01 * (
			 2.61413190e+00 + T * (  6.76078200e-02 
			+ T * (  1.11900680e-05 + T * ( -5.64123400e-08
			+ T *  2.54839210e-11 ) ) ) );
		h[s2XC10H21] =  5.88550840e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s2XC10H21] =  5.88550840e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s3XC10H21] =  5.88550840e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s3XC10H21] =  5.88550840e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s4XC10H21] =  5.88550840e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s4XC10H21] =  5.88550840e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s5XC10H21] =  5.88550840e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s5XC10H21] =  5.88550840e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s1XC10H21] =  5.88550840e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s1XC10H21] =  5.88550840e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sNXC10H22] =  5.84381062e+01 * (
			T * (  3.01495460e+00 + T * (  2.72286015e-02
			+ T * (  7.27089367e-06 + T * ( -1.35585278e-08
			+ T *  4.16174600e-12 ) ) ) ) -2.60033790e+04 );
		cp[sNXC10H22] =  5.84381062e+01 * (
			 3.01495460e+00 + T * (  5.44572030e-02 
			+ T * (  2.18126810e-05 + T * ( -5.42341110e-08
			+ T *  2.08087300e-11 ) ) ) );
		h[s2XRO2XC10H21O2] =  4.79854330e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s2XRO2XC10H21O2] =  4.79854330e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s3X2XRO2HXC10H20OOH] =  4.79854330e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s3X2XRO2HXC10H20OOH] =  4.79854330e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s32O2RO2HXO2C10H20OOH] =  4.05048035e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s32O2RO2HXO2C10H20OOH] =  4.05048035e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s3X2HOROHXHO2C10H19OOH] =  4.05048035e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s3X2HOROHXHO2C10H19OOH] =  4.05048035e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s3X2ORO2HXOC10H19OOH] =  4.41641347e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s3X2ORO2HXOC10H19OOH] =  4.41641347e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[s2X1OROXOC10H19O] =  4.85503235e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[s2X1OROXOC10H19O] =  4.85503235e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
		h[sUXC6H5] =  1.07838392e+02 * (
			T * (  4.85426800e+00 + T * (  1.51582950e-02
			+ T * (  5.80964000e-08 + T * ( -4.52752500e-09
			+ T *  1.47850220e-12 ) ) ) ) +  6.81732300e+04 );
		cp[sUXC6H5] =  1.07838392e+02 * (
			 4.85426800e+00 + T * (  3.03165900e-02 
			+ T * (  1.74289200e-07 + T * ( -1.81101000e-08
			+ T *  7.39251100e-12 ) ) ) );
		h[sC6H6] =  1.06446715e+02 * (
			T * ( -3.52109420e+00 + T * (  4.03026715e-02
			+ T * ( -3.36291400e-05 + T * (  1.71465552e-08
			+ T * -3.79215020e-12 ) ) ) ) +  4.10422350e+04 );
		cp[sC6H6] =  1.06446715e+02 * (
			-3.52109420e+00 + T * (  8.06053430e-02 
			+ T * ( -1.00887420e-04 + T * (  6.85862210e-08
			+ T * -1.89607510e-11 ) ) ) );
		h[sUXC6H7] =  1.05090500e+02 * (
			T * (  5.78835200e+00 + T * (  9.85662000e-03
			+ T * (  1.29601040e-05 + T * ( -1.27202075e-08
			+ T *  3.10446500e-12 ) ) ) ) +  4.74284760e+04 );
		cp[sUXC6H7] =  1.05090500e+02 * (
			 5.78835200e+00 + T * (  1.97132400e-02 
			+ T * (  3.88803120e-05 + T * ( -5.08808300e-08
			+ T *  1.55223250e-11 ) ) ) );
		h[sA1XC6H6] =  1.06446715e+02 * (
			T * ( -5.15447510e+00 + T * (  3.04725240e-02
			+ T * ( -1.21158027e-05 + T * (  1.67725935e-10
			+ T *  1.05215958e-12 ) ) ) ) +  9.10877510e+03 );
		cp[sA1XC6H6] =  1.06446715e+02 * (
			-5.15447510e+00 + T * (  6.09450480e-02 
			+ T * ( -3.63474080e-05 + T * (  6.70903740e-10
			+ T *  5.26079790e-12 ) ) ) );
		h[sA1XXC6H5] =  1.07838392e+02 * (
			T * ( -4.14205870e+00 + T * (  2.77593785e-02
			+ T * ( -1.11881920e-05 + T * (  2.56120050e-10
			+ T *  9.27311400e-13 ) ) ) ) +  3.84103110e+04 );
		cp[sA1XXC6H5] =  1.07838392e+02 * (
			-4.14205870e+00 + T * (  5.55187570e-02 
			+ T * ( -3.35645760e-05 + T * (  1.02448020e-09
			+ T *  4.63655700e-12 ) ) ) );
		h[sC6H5O] =  8.93054780e+01 * (
			T * ( -1.82194330e+00 + T * (  2.40612550e-02
			+ T * ( -1.55974340e-06 + T * ( -8.50464850e-09
			+ T *  3.72992740e-12 ) ) ) ) +  4.24291800e+03 );
		cp[sC6H5O] =  8.93054780e+01 * (
			-1.82194330e+00 + T * (  4.81225100e-02 
			+ T * ( -4.67923020e-06 + T * ( -3.40185940e-08
			+ T *  1.86496370e-11 ) ) ) );
		h[sC6H5OH] =  8.83489183e+01 * (
			T * ( -1.66776560e+00 + T * (  2.60455850e-02
			+ T * ( -2.26316690e-06 + T * ( -9.06427900e-09
			+ T *  4.11724820e-12 ) ) ) ) -1.32875390e+04 );
		cp[sC6H5OH] =  8.83489183e+01 * (
			-1.66776560e+00 + T * (  5.20911700e-02 
			+ T * ( -6.78950070e-06 + T * ( -3.62571160e-08
			+ T *  2.05862410e-11 ) ) ) );
		h[sC5H5] =  1.27736058e+02 * (
			T * ( -3.19673920e+00 + T * (  2.04068050e-02
			+ T * (  2.27216843e-07 + T * ( -7.84364750e-09
			+ T *  3.15444620e-12 ) ) ) ) +  1.52906760e+04 );
		cp[sC5H5] =  1.27736058e+02 * (
			-3.19673920e+00 + T * (  4.08136100e-02 
			+ T * (  6.81650530e-07 + T * ( -3.13745900e-08
			+ T *  1.57722310e-11 ) ) ) );
		h[sC5H6] =  1.25788072e+02 * (
			T * ( -3.19673920e+00 + T * (  2.04068050e-02
			+ T * (  2.27216843e-07 + T * ( -7.84364750e-09
			+ T *  3.15444620e-12 ) ) ) ) +  1.52906760e+04 );
		cp[sC5H6] =  1.25788072e+02 * (
			-3.19673920e+00 + T * (  4.08136100e-02 
			+ T * (  6.81650530e-07 + T * ( -3.13745900e-08
			+ T *  1.57722310e-11 ) ) ) );
		h[sC5H5O] =  1.02532248e+02 * (
			T * ( -1.82194330e+00 + T * (  2.40612550e-02
			+ T * ( -1.55974340e-06 + T * ( -8.50464850e-09
			+ T *  3.72992740e-12 ) ) ) ) +  4.24291800e+03 );
		cp[sC5H5O] =  1.02532248e+02 * (
			-1.82194330e+00 + T * (  4.81225100e-02 
			+ T * ( -4.67923020e-06 + T * ( -3.40185940e-08
			+ T *  1.86496370e-11 ) ) ) );
		h[sA1C2HXC8H6] =  8.14109745e+01 * (
			T * ( -5.57356520e+00 + T * (  4.33654820e-02
			+ T * ( -2.80397020e-05 + T * (  1.04617595e-08
			+ T * -1.62112044e-12 ) ) ) ) +  3.63172340e+04 );
		cp[sA1C2HXC8H6] =  8.14109745e+01 * (
			-5.57356520e+00 + T * (  8.67309640e-02 
			+ T * ( -8.41191060e-05 + T * (  4.18470380e-08
			+ T * -8.10560220e-12 ) ) ) );
		h[sA1C2HYXC8H5] =  8.22225079e+01 * (
			T * ( -4.56114890e+00 + T * (  4.06523370e-02
			+ T * ( -2.71120913e-05 + T * (  1.05501538e-08
			+ T * -1.74596862e-12 ) ) ) ) +  6.56187700e+04 );
		cp[sA1C2HYXC8H5] =  8.22225079e+01 * (
			-4.56114890e+00 + T * (  8.13046740e-02 
			+ T * ( -8.13362740e-05 + T * (  4.22006150e-08
			+ T * -8.72984310e-12 ) ) ) );
		h[sA1C2HACXC10H7] =  6.53869263e+01 * (
			T * ( -2.87461430e+00 + T * (  4.67993190e-02
			+ T * ( -2.90542303e-05 + T * (  1.05309350e-08
			+ T * -1.60301288e-12 ) ) ) ) +  7.11222060e+04 );
		cp[sA1C2HACXC10H7] =  6.53869263e+01 * (
			-2.87461430e+00 + T * (  9.35986380e-02 
			+ T * ( -8.71626910e-05 + T * (  4.21237400e-08
			+ T * -8.01506440e-12 ) ) ) );
		h[sA2XXXC10H7] =  6.53869263e+01 * (
			T * ( -6.02056560e+00 + T * (  4.52001605e-02
			+ T * ( -2.05431377e-05 + T * (  2.84409750e-09
			+ T *  7.00488980e-13 ) ) ) ) +  4.57999680e+04 );
		cp[sA2XXXC10H7] =  6.53869263e+01 * (
			-6.02056560e+00 + T * (  9.04003210e-02 
			+ T * ( -6.16294130e-05 + T * (  1.13763900e-08
			+ T *  3.50244490e-12 ) ) ) );
		h[sA2XC10H8] =  6.48726632e+01 * (
			T * ( -7.03298200e+00 + T * (  4.79133060e-02
			+ T * ( -2.14707483e-05 + T * (  2.75570325e-09
			+ T *  8.25337180e-13 ) ) ) ) +  1.64984320e+04 );
		cp[sA2XC10H8] =  6.48726632e+01 * (
			-7.03298200e+00 + T * (  9.58266120e-02 
			+ T * ( -6.44122450e-05 + T * (  1.10228130e-08
			+ T *  4.12668590e-12 ) ) ) );
		h[sA2R5XC12H8] =  5.46334700e+01 * (
			T * ( -8.55541730e+00 + T * (  5.94867900e-02
			+ T * ( -3.34815967e-05 + T * (  1.01576480e-08
			+ T * -1.15698604e-12 ) ) ) ) +  2.90701080e+04 );
		cp[sA2R5XC12H8] =  5.46334700e+01 * (
			-8.55541730e+00 + T * (  1.18973580e-01 
			+ T * ( -1.00444790e-04 + T * (  4.06305920e-08
			+ T * -5.78493020e-12 ) ) ) );
		h[sA2R5XXC12H7] =  5.49977510e+01 * (
			T * ( -7.54300090e+00 + T * (  5.67736450e-02
			+ T * ( -3.25539857e-05 + T * (  1.02460423e-08
			+ T * -1.28183422e-12 ) ) ) ) +  5.83716440e+04 );
		cp[sA2R5XXC12H7] =  5.49977510e+01 * (
			-7.54300090e+00 + T * (  1.13547290e-01 
			+ T * ( -9.76619570e-05 + T * (  4.09841690e-08
			+ T * -6.40917110e-12 ) ) ) );
		h[sA2R5C2HXC14H8] =  4.71858755e+01 * (
			T * ( -8.97450750e+00 + T * (  7.23797500e-02
			+ T * ( -4.94054967e-05 + T * (  2.04516818e-08
			+ T * -3.83026600e-12 ) ) ) ) +  5.62785670e+04 );
		cp[sA2R5C2HXC14H8] =  4.71858755e+01 * (
			-8.97450750e+00 + T * (  1.44759500e-01 
			+ T * ( -1.48216490e-04 + T * (  8.18067270e-08
			+ T * -1.91513300e-11 ) ) ) );
		h[sA2R5C2HYXC14H7] =  4.74573620e+01 * (
			T * ( -7.96209110e+00 + T * (  6.96666050e-02
			+ T * ( -4.84778833e-05 + T * (  2.05400758e-08
			+ T * -3.95511420e-12 ) ) ) ) +  8.55801030e+04 );
		cp[sA2R5C2HYXC14H7] =  4.74573620e+01 * (
			-7.96209110e+00 + T * (  1.39333210e-01 
			+ T * ( -1.45433650e-04 + T * (  8.21603030e-08
			+ T * -1.97755710e-11 ) ) ) );
		h[sANC2HACXC16H9] =  4.13171861e+01 * (
			T * ( -6.27555650e+00 + T * (  7.58135850e-02
			+ T * ( -5.04200233e-05 + T * (  2.05208570e-08
			+ T * -3.81215840e-12 ) ) ) ) +  9.10835390e+04 );
		cp[sANC2HACXC16H9] =  4.13171861e+01 * (
			-6.27555650e+00 + T * (  1.51627170e-01 
			+ T * ( -1.51260070e-04 + T * (  8.20834280e-08
			+ T * -1.90607920e-11 ) ) ) );
		h[sA3R5XXC16H9] =  4.13171861e+01 * (
			T * ( -9.42150780e+00 + T * (  7.42144300e-02
			+ T * ( -4.19089300e-05 + T * (  1.28340195e-08
			+ T * -1.50865664e-12 ) ) ) ) +  6.46542500e+04 );
		cp[sA3R5XXC16H9] =  4.13171861e+01 * (
			-9.42150780e+00 + T * (  1.48428860e-01 
			+ T * ( -1.25726790e-04 + T * (  5.13360780e-08
			+ T * -7.54328320e-12 ) ) ) );
		h[sA3R5XC16H10] =  4.11112540e+01 * (
			T * ( -1.04339240e+01 + T * (  7.69275750e-02
			+ T * ( -4.28365400e-05 + T * (  1.27456255e-08
			+ T * -1.38380846e-12 ) ) ) ) +  3.53527140e+04 );
		cp[sA3R5XC16H10] =  4.11112540e+01 * (
			-1.04339240e+01 + T * (  1.53855150e-01 
			+ T * ( -1.28509620e-04 + T * (  5.09825020e-08
			+ T * -6.91904230e-12 ) ) ) );
		h[sA3R5ACXC18H11] =  3.65838569e+01 * (
			T * ( -7.73497320e+00 + T * (  8.03614100e-02
			+ T * ( -4.38510700e-05 + T * (  1.28148007e-08
			+ T * -1.36570086e-12 ) ) ) ) +  7.01576860e+04 );
		cp[sA3R5ACXC18H11] =  3.65838569e+01 * (
			-7.73497320e+00 + T * (  1.60722820e-01 
			+ T * ( -1.31553210e-04 + T * (  5.12592030e-08
			+ T * -6.82850430e-12 ) ) ) );
		h[sA4XC18H10] =  3.67468399e+01 * (
			T * ( -1.42954020e+01 + T * (  9.08867350e-02
			+ T * ( -5.35984000e-05 + T * (  1.68837257e-08
			+ T * -1.97088082e-12 ) ) ) ) +  3.80833950e+04 );
		cp[sA4XC18H10] =  3.67468399e+01 * (
			-1.42954020e+01 + T * (  1.81773470e-01 
			+ T * ( -1.60795200e-04 + T * (  6.75349030e-08
			+ T * -9.85440410e-12 ) ) ) );
		h[sA2CH3XC11H10] =  5.84734510e+01 * (
			T * ( -7.50660900e-02 + T * (  3.47692260e-02
			+ T * ( -1.70670320e-06 + T * ( -8.90050350e-09
			+ T *  3.36666780e-12 ) ) ) ) +  1.04955660e+04 );
		cp[sA2CH3XC11H10] =  5.84734510e+01 * (
			-7.50660900e-02 + T * (  6.95384520e-02 
			+ T * ( -5.12010960e-06 + T * ( -3.56020140e-08
			+ T *  1.68333390e-11 ) ) ) );
		h[sA2CH2XC11H9] =  5.88909351e+01 * (
			T * (  3.75132530e-01 + T * (  3.45229840e-02
			+ T * ( -2.22163000e-06 + T * ( -9.12558650e-09
			+ T *  3.61566340e-12 ) ) ) ) +  2.72279280e+04 );
		cp[sA2CH2XC11H9] =  5.88909351e+01 * (
			 3.75132530e-01 + T * (  6.90459680e-02 
			+ T * ( -6.66489000e-06 + T * ( -3.65023460e-08
			+ T *  1.80783170e-11 ) ) ) );
		h[sA2CHOXC11H8O] =  5.32376708e+01 * (
			T * ( -7.29367490e+00 + T * (  3.52177990e-02
			+ T * ( -4.15382433e-06 + T * ( -1.03542845e-08
			+ T *  5.07184220e-12 ) ) ) ) +  5.95325930e+03 );
		cp[sA2CHOXC11H8O] =  5.32376708e+01 * (
			-7.29367490e+00 + T * (  7.04355980e-02 
			+ T * ( -1.24614730e-05 + T * ( -4.14171380e-08
			+ T *  2.53592110e-11 ) ) ) );
		h[sBIBENZYLXC22H18] =  2.94454676e+01 * (
			T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) ) +  0.00000000e+00 );
		cp[sBIBENZYLXC22H18] =  2.94454676e+01 * (
			 0.00000000e+00 + T * (  0.00000000e+00 
			+ T * (  0.00000000e+00 + T * (  0.00000000e+00
			+ T *  0.00000000e+00 ) ) ) );
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
	return 93;
}

int GetNReactions( void )
{
	return rEnd;
}

double MAX_C(double X1, double X2)
{
  return ( (X1 > X2) ? X1 : X2 );
}
