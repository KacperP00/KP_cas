#ifndef MDUC_UDF_MEMORY_H
#define MDUC_UDF_MEMORY_H

/* ----- DEFINE MACROS FOR USER-DEFINED VARIABLES ----- */
/* ----- user-defined memory  ----- */
static int udmOffset = 0; /* UDM_UNRESERVED; */

/* storage for reference scalar dissipation rates */
#define C_CHI0_REF(c,t)     C_UDMI(c,t,udmOffset+udf.iCHIR[0])
#define C_CHI1_REF(c,t)     C_UDMI(c,t,udmOffset+udf.iCHIR[1])
#define C_CHI2_REF(c,t)     C_UDMI(c,t,udmOffset+udf.iCHIR[2])

#define C_HR_RATE(c,t)      C_UDMI(c,t,udmOffset+udf.iHRrate)
#define C_HR_PREM(c,t)      C_UDMI(c,t,udmOffset+udf.iHRprem)

/* additional memory require for single droplet evaporation model (SDM) */
#ifdef MDUC_SPRAY 
/* the rate of fuel volume evaporated in a cell */
#define C_FUEL1_VOL(c,t)    C_UDMI(c,t,udmOffset+udf.iFuelV[1])
#define C_FUEL2_VOL(c,t)    C_UDMI(c,t,udmOffset+udf.iFuelV[2])

/* source terms for variance equation */
#define C_Z0VAR_SRC(c,t)    C_UDMI(c,t,udmOffset+udf.iZV_SRC[0])
#define C_Z1VAR_SRC(c,t)    C_UDMI(c,t,udmOffset+udf.iZV_SRC[1])
#define C_Z2VAR_SRC(c,t)    C_UDMI(c,t,udmOffset+udf.iZV_SRC[2])

/* dissipation terms for variance equation */
#define C_Z0VAR_EPS(c,t)    C_UDMI(c,t,udmOffset+udf.iZV_EPS[0])
#define C_Z1VAR_EPS(c,t)    C_UDMI(c,t,udmOffset+udf.iZV_EPS[1])
#define C_Z2VAR_EPS(c,t)    C_UDMI(c,t,udmOffset+udf.iZV_EPS[2])

/* turbulent diffusivity correction factor */
#define C_FDT_0(c,t)        C_UDMI(c,t,udmOffset+udf.iFDT[0])
#define C_FDT_1(c,t)        C_UDMI(c,t,udmOffset+udf.iFDT[1])
#define C_FDT_2(c,t)        C_UDMI(c,t,udmOffset+udf.iFDT[2])

/* premixed combustion variables */
#define C_REGIME(c,t)       C_UDMI(c,t,udmOffset+udf.iRegime)
#define C_PROGSRC(c,t)      C_UDMI(c,t,udmOffset+udf.iProgSrc)
#define C_RIFSRC(c,t)       C_UDMI(c,t,udmOffset+udf.iRIFSrc)

#define C_PSCALE(c,t)       C_UDMI(c,t,udmOffset+udf.iPscale)
#define C_SL(c,t)           C_UDMI(c,t,udmOffset+udf.iSL)
#define C_ST(c,t)           C_UDMI(c,t,udmOffset+udf.iST)
#define C_CRIF(c,t)         C_UDMI(c,t,udmOffset+udf.iCrif)
#define C_PROGMAX(c,t)      C_UDMI(c,t,udmOffset+udf.iProgMax)
#define C_PROGMIN(c,t)      C_UDMI(c,t,udmOffset+udf.iProgMin)
#define C_DIFFSC(c,t)       C_UDMI(c,t,udmOffset+udf.iDiffSc)

#ifdef DEBUG_SDM
/* constants/exponents for SDM model */
#define C_XI_0(c,t)         C_UDMI(c,t,udmOffset+udf.iXI[0])
#define C_XI_1(c,t)         C_UDMI(c,t,udmOffset+udf.iXI[1])
#define C_XI_2(c,t)         C_UDMI(c,t,udmOffset+udf.iXI[2])
#define C_ALPHA_B0(c,t)     C_UDMI(c,t,udmOffset+udf.iAlphaB[0])
#define C_ALPHA_B1(c,t)     C_UDMI(c,t,udmOffset+udf.iAlphaB[1])
#define C_ALPHA_B2(c,t)     C_UDMI(c,t,udmOffset+udf.iAlphaB[2])

/* Spalding number */
#define C_SPALD_NO(c,t)     C_UDMI(c,t,udmOffset+udf.iSPALD)

#define C_Z0VAR_PRD(c,t)    C_UDMI(c,t,udmOffset+udf.iZV_PRD[0])
#define C_Z1VAR_PRD(c,t)    C_UDMI(c,t,udmOffset+udf.iZV_PRD[1])
#define C_Z2VAR_PRD(c,t)    C_UDMI(c,t,udmOffset+udf.iZV_PRD[2])

#define C_Z0VAR_D2(c,t)     C_UDMI(c,t,udmOffset+udf.iZV_D2[0])
#define C_Z1VAR_D2(c,t)     C_UDMI(c,t,udmOffset+udf.iZV_D2[1])
#define C_Z2VAR_D2(c,t)     C_UDMI(c,t,udmOffset+udf.iZV_D2[2])

#define C_Z0VAR_SRC_R(c,t)  C_UDMI(c,t,udmOffset+udf.iZV_SRC_R[0])
#define C_Z1VAR_SRC_R(c,t)  C_UDMI(c,t,udmOffset+udf.iZV_SRC_R[1])
#define C_Z2VAR_SRC_R(c,t)  C_UDMI(c,t,udmOffset+udf.iZV_SRC_R[2])
#define C_Z0VAR_EPS_R(c,t)  C_UDMI(c,t,udmOffset+udf.iZV_EPS_R[0])
#define C_Z1VAR_EPS_R(c,t)  C_UDMI(c,t,udmOffset+udf.iZV_EPS_R[1])
#define C_Z2VAR_EPS_R(c,t)  C_UDMI(c,t,udmOffset+udf.iZV_EPS_R[2])

#define C_Z0VAR_SRC_L(c,t)  C_UDMI(c,t,udmOffset+udf.iZV_SRC_L[0])
#define C_Z1VAR_SRC_L(c,t)  C_UDMI(c,t,udmOffset+udf.iZV_SRC_L[1])
#define C_Z2VAR_SRC_L(c,t)  C_UDMI(c,t,udmOffset+udf.iZV_SRC_L[2])
#define C_RE_TURB(c,t)      C_UDMI(c,t,udmOffset+udf.iReTURB)
#endif /* DEBUG_SDM */

#endif /* MDUC_SPRAY */

void UDFsetMemoryIndex();
int  UDFassignMem(int index, char *id, int *var, char **name);

#endif /* MDUC_UDF_MEMORY_H */
