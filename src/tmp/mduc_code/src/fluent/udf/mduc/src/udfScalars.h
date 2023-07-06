#ifndef MDUC_UDF_SCALARS_H
#define MDUC_UDF_SCALARS_H

/* ----- define macros for user defined scalars  ----- */
static int udsOffset = 0; /* UDS_UNRESERVED; */
#define C_Z1MEAN(c,t) C_UDSI(c,t,udsOffset+udf.iZ1)
#define C_Z1VAR(c,t)  C_UDSI(c,t,udsOffset+udf.iZ1V)
#define C_Z2MEAN(c,t) C_UDSI(c,t,udsOffset+udf.iZ2)
#define C_Z2VAR(c,t)  C_UDSI(c,t,udsOffset+udf.iZ2V)
#define C_Z0VAR(c,t)  C_UDSI(c,t,udsOffset+udf.iZ0V)
#define C_HTOT(c,t)   C_UDSI(c,t,udsOffset+udf.iH)
#define C_HVAR(c,t)   C_UDSI(c,t,udsOffset+udf.iHV)
#define C_HMOD(c,t)   C_UDSI(c,t,udsOffset+udf.iHMOD)
#define C_EGR(c,t)    C_UDSI(c,t,udsOffset+udf.iEGR)
#define C_EGRVAR(c,t) C_UDSI(c,t,udsOffset+udf.iEGRV)
#define C_PROGTOTAL(c,t)  C_UDSI(c,t,udsOffset+udf.iProgTotal)

/* scalar gradients */
#define C_Z1MEAN_G(c,t) C_UDSI_G(c,t,udsOffset+udf.iZ1)
#define C_Z1VAR_G(c,t)  C_UDSI_G(c,t,udsOffset+udf.iZ1V)
#define C_Z2MEAN_G(c,t) C_UDSI_G(c,t,udsOffset+udf.iZ2)
#define C_Z2VAR_G(c,t)  C_UDSI_G(c,t,udsOffset+udf.iZ2V)
#define C_Z0VAR_G(c,t)  C_UDSI_G(c,t,udsOffset+udf.iZ0V)
#define C_HTOT_G(c,t)   C_UDSI_G(c,t,udsOffset+udf.iH)
#define C_HVAR_G(c,t)   C_UDSI_G(c,t,udsOffset+udf.iHV)
#define C_HMOD_G(c,t)   C_UDSI_G(c,t,udsOffset+udf.iHMOD)
#define C_EGR_G(c,t)    C_UDSI_G(c,t,udsOffset+udf.iEGR)
#define C_EGRVAR_G(c,t) C_UDSI_G(c,t,udsOffset+udf.iEGRV)
#define C_PROGTOTAL_G(c,t)  C_UDSI_G(c,t,udsOffset+udf.iProgTotal)

/* general IDs */
#define C_PHI1(c,t)     C_UDSI(c,t,udsOffset+udf.iPhi[1])
#define C_PHIV1(c,t)    C_UDSI(c,t,udsOffset+udf.iPhiV[1])
#define C_PHI2(c,t)     C_UDSI(c,t,udsOffset+udf.iPhi[2])
#define C_PHIV2(c,t)    C_UDSI(c,t,udsOffset+udf.iPhiV[2])
#define C_PHIV0(c,t)    C_UDSI(c,t,udsOffset+udf.iPhiV[0])

/* function to set default values at beginning of iteration */
void UDFsetScalarDefaultSources(Domain *d);

/* functions for mapping scalars correctly in library */
void UDFsetScalarIndex();

#endif /* MDUC_UDF_SCALARS_H */
