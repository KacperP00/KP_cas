#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include "iFileIO.h"

iFileStruct * iFileRead(char *fileName)
{
  int i, j, k;
  int debug = 0;
  char *f;
  long lDum;
  long endAPB, endDGB, endHeader;
  long int len;

  iFileStruct *iFile = (iFileStruct*) calloc(1, sizeof(iFileStruct));
  iFileOP *op = (iFileOP*) calloc(1, sizeof(iFileOP));

  /* open the iFile for reading */
  FILE *fP = fopen(fileName, "rb");
  if (fP == NULL) {
    fprintf(stderr, "Error: couldn't open iFile %s\n", fileName);
    return NULL;
  }

  char *buff = (char*) calloc(1024, sizeof(char));
  char cTmp;

  // [ a, b, c, bt, br, Pext, datah] = iFileReadMat_v4('M271SCRE_HCCI_baseline_map3.350')

  /* allocate enough memory for hard coded operational parameters */
  iFile->apb.nOp = 28; /* hard coded OPs */
  iFile->apb.op = (iFileOP*) calloc(iFile->apb.nOp, sizeof(iFileOP));

  /* ----- read in the general parameters ----- */
  fread(&iFile->apb.parlng,     sizeof(short),         1, fP);
  fread(&iFile->apb.grpanz,     sizeof(short),         1, fP);
  fread(iFile->apb.filkom,      sizeof(char),         80, fP);
  fread(iFile->apb.mesdat,      sizeof(char),         24, fP);
  fread(iFile->apb.parfil,      sizeof(char),         18, fP);
  fread(iFile->apb.pronam,      sizeof(char),         18, fP);
  fread(iFile->apb.prfstd,      sizeof(char),         10, fP);
  fread(iFile->apb.motnam,      sizeof(char),         18, fP);
  fread(&iFile->apb.mottyp,     sizeof(short),         1, fP);
  fread(&iFile->apb.tktanz,     sizeof(short),         1, fP);
  fread(&iFile->apb.geoein,     sizeof(short),         1, fP);
  fread(&iFile->apb.hublng,     sizeof(double),        1, fP);
  fread(&iFile->apb.pleuel,     sizeof(double),        1, fP);
  fread(&iFile->apb.bohrng,     sizeof(double),        1, fP);
  fread(&iFile->apb.kompre,     sizeof(double),        1, fP);
  fread(&iFile->apb.desaxi,     sizeof(double),        1, fP);
  fread(&iFile->apb.grpflg,     sizeof(int32_t),       1, fP);
  fread(&iFile->apb.datahandle, sizeof(int32_t),       1, fP);
  fread(&iFile->apb.version,    sizeof(short),         1, fP);
  fread(&iFile->apb.desaxi2,    sizeof(double),        1, fP);
  fread(&iFile->apb.pExt,       sizeof(uint32_t),      1, fP);

  /* ensure that strings end with null characters */
  iFile->apb.filkom[79] = '\0';
  iFile->apb.mesdat[23] = '\0';
  iFile->apb.parfil[17] = '\0';
  iFile->apb.pronam[17] = '\0';
  iFile->apb.prfstd[9]  = '\0';
  iFile->apb.motnam[17] = '\0';

  /* determine whether the file needs an endian swap */
  bool swap_end;
  if (iFile->apb.parlng == 1024)
    swap_end = false;
  else if (iFile->apb.parlng == 4)
    swap_end = true;
  else {
    fprintf(stderr, "Error: Endianess of iFile not recognized!\n");
    exit(1);
  }

  if (swap_end)
    fprintf(stderr, "#Info: Swapping endian in iFile\n");

  if (iFile->apb.grpanz == 0) {
    fprintf(stderr, "Error: iFile is empty\n");
    exit(2);
  }

  if (debug) {
    fprintf(stderr, "iFile parameter block:\n");
    fprintf(stderr, "  parlng: %hi\n", iFile->apb.parlng);
    fprintf(stderr, "  grpanz: %hi\n", iFile->apb.grpanz);
    fprintf(stderr, "  filkom: %s\n", iFile->apb.filkom);
    fprintf(stderr, "  mesdat: %s\n", iFile->apb.mesdat);
    fprintf(stderr, "  parfil: %s\n", iFile->apb.parfil);
    fprintf(stderr, "  pronam: %s\n", iFile->apb.pronam);
    fprintf(stderr, "  prfstd: %s\n", iFile->apb.prfstd);
    fprintf(stderr, "  motnam: %s\n", iFile->apb.motnam);
    fprintf(stderr, "  mottyp: %hi\n", iFile->apb.mottyp);
    fprintf(stderr, "  tktanz: %hi\n", iFile->apb.tktanz);
    fprintf(stderr, "  geoein: %hi\n", iFile->apb.geoein);
    fprintf(stderr, "  hublng: %g\n", iFile->apb.hublng);
    fprintf(stderr, "  pleuel: %g\n", iFile->apb.pleuel);
    fprintf(stderr, "  bohrng: %g\n", iFile->apb.bohrng);
    fprintf(stderr, "  kompre: %g\n", iFile->apb.kompre);
    fprintf(stderr, "  desaxi: %g\n", iFile->apb.desaxi);
    fprintf(stderr, "  grpflg: %i\n", iFile->apb.grpflg);
    fprintf(stderr, "  datahandle: %i\n", iFile->apb.datahandle);
    fprintf(stderr, "  version: %i\n", iFile->apb.version);
    fprintf(stderr, "  desaxi2: %f\n", iFile->apb.desaxi2);
    fprintf(stderr, "  pExt:    %i\n", iFile->apb.pExt);
  }

  /* read hard-coded operating parameters */
  for (i=0; i<iFile->apb.nOp; i++) {
    fread(iFile->apb.op[i].name,   sizeof(char),   10, fP);
    fread(iFile->apb.op[i].unit,   sizeof(char),   10, fP);
    fread(&iFile->apb.op[i].value, sizeof(double),  1, fP);

    iFile->apb.op[i].name[9] = '\0';
    iFile->apb.op[i].unit[9] = '\0';
  } /* end loop over hard-coded operating parameters */

  /* check whether the header size read in was correct */
  endAPB = ftell(fP);
  if (endAPB != iFile->apb.parlng) {
    fprintf(stderr, "Warning: Header was read as %li bytes long!\n", endAPB);
    fprintf(stderr, "         file wants %hi bytes\n", iFile->apb.parlng);
  }

  if (swap_end) {
    swap_endian(&iFile->apb.parlng,     sizeof(short));
    swap_endian(&iFile->apb.grpanz,     sizeof(short));
    swap_endian(&iFile->apb.mottyp,     sizeof(short));
    swap_endian(&iFile->apb.tktanz,     sizeof(short));
    swap_endian(&iFile->apb.geoein,     sizeof(short));
    swap_endian(&iFile->apb.hublng,     sizeof(double));
    swap_endian(&iFile->apb.pleuel,     sizeof(double));
    swap_endian(&iFile->apb.bohrng,     sizeof(double));
    swap_endian(&iFile->apb.kompre,     sizeof(double));
    swap_endian(&iFile->apb.desaxi,     sizeof(double));
    swap_endian(&iFile->apb.grpflg,     sizeof(int32_t));
    swap_endian(&iFile->apb.datahandle, sizeof(int32_t));
    swap_endian(&iFile->apb.version,    sizeof(short));
    swap_endian(&iFile->apb.desaxi2,    sizeof(double));
    swap_endian(&iFile->apb.pExt,       sizeof(uint32_t));
  }

  /* read in extended operating parameters */
  int nExOp;
  if (iFile->apb.pExt > 0) {
    if (debug)
      fprintf(stderr, " Reading extended operating parameters\n");
    nExOp = iFileReadExtendedOPs(fP, iFile->apb.pExt, &iFile->apb.op);
    iFile->apb.nOp += nExOp;
  }

  /* output the variables */
  fprintf(stderr, "=====================================================\n");
  fprintf(stderr, "|                    Variables                      |\n");
  fprintf(stderr, "+-----+--------------+-----------------+------------+\n");
  fprintf(stderr, "| Var | %-12s | %-15s | %-10s |\n", "Name","Value","Units");
  fprintf(stderr, "+-----+--------------+-----------------+------------+\n");
  for (i=0; i<iFile->apb.nOp; i++) {
    fprintf(stderr, "| %3i | %-12s | %15g | %-10s |\n", i,
            iFile->apb.op[i].name, iFile->apb.op[i].value,
            iFile->apb.op[i].unit);
  }
  fprintf(stderr, "+-----+--------------+-----------------+------------+\n\n");

  /* ----- read DGB ----- */
  iFileDATGRP *iDGB;
  iFileDAT *pDI;
  iFile->dgb = (iFileDATGRP*) calloc(iFile->apb.grpanz, sizeof(iFileDATGRP));
  for (i=0; i<iFile->apb.grpanz; i++) {
    fread(&iFile->dgb[i].datart,    sizeof(short),   1, fP);
    fread(&iFile->dgb[i].absart,    sizeof(short),   1, fP);
    fread(iFile->dgb[i].absein,     sizeof(char),   10, fP);
    fread(&iFile->dgb[i].dltphi,    sizeof(double),  1, fP);
    fread(&iFile->dgb[i].dltzei,    sizeof(double),  1, fP);
    fread(&iFile->dgb[i].kananz,    sizeof(short),   1, fP);
    fread(&iFile->dgb[i].absanf,    sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].absend,    sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].fortyp,    sizeof(short),   1, fP);
    fread(&iFile->dgb[i].beranz,    sizeof(short),   1, fP);
    fread(&iFile->dgb[i].blkanz,    sizeof(short),   1, fP);
    fread(&iFile->dgb[i].zykanz,    sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].zykanf,    sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].zykend,    sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].zyklng,    sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].zyksta,    sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].zyksum,    sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].datfor,    sizeof(short),   1, fP);
    fread(&iFile->dgb[i].mitanz,    sizeof(short),   1, fP);
    fread(&iFile->dgb[i].mitkor,    sizeof(short),   1, fP);
    fread(&iFile->dgb[i].datofs,    sizeof(short),   1, fP);
    fread(&iFile->dgb[i].diradr,                4,   1, fP);
    fread(&iFile->dgb[i].absadr,                4,   1, fP);
    fread(&iFile->dgb[i].mpladr,                4,   1, fP);
    fread(&iFile->dgb[i].aztadr,                4,   1, fP);
    fread(&iFile->dgb[i].rztadr,                4,   1, fP);
    fread(&iFile->dgb[i].datadr,                4,   1, fP);

    iFile->dgb[i].absein[9] = '\0';

    if (debug) {
      fprintf(stderr, " datart: %hi\n", iFile->dgb[i].datart);
      fprintf(stderr, " absart: %hi\n", iFile->dgb[i].absart);
      fprintf(stderr, " absein: %s\n", iFile->dgb[i].absein);
      fprintf(stderr, " dltphi: %f\n", iFile->dgb[i].dltphi);
      fprintf(stderr, " dltzei: %f\n", iFile->dgb[i].dltzei);
      fprintf(stderr, " kananz: %hi\n", iFile->dgb[i].kananz);
      fprintf(stderr, " absanf: %i\n", iFile->dgb[i].absanf);
      fprintf(stderr, " absend: %i\n", iFile->dgb[i].absend);
      fprintf(stderr, " fortyp: %hi\n", iFile->dgb[i].fortyp);
      fprintf(stderr, " beranz: %hi\n", iFile->dgb[i].beranz);
      fprintf(stderr, " blkanz: %hi\n", iFile->dgb[i].blkanz);
      fprintf(stderr, " zykanz: %i\n", iFile->dgb[i].zykanz);
      fprintf(stderr, " zykanf: %i\n", iFile->dgb[i].zykanf);
      fprintf(stderr, " zykend: %i\n", iFile->dgb[i].zykend);
      fprintf(stderr, " zyklng: %i\n", iFile->dgb[i].zyklng);
      fprintf(stderr, " zyksta: %i\n", iFile->dgb[i].zyksta);
      fprintf(stderr, " zyksum: %i\n", iFile->dgb[i].zyksum);
      fprintf(stderr, " datfor: %hi\n", iFile->dgb[i].datfor);
      fprintf(stderr, " mitanz: %hi\n", iFile->dgb[i].mitanz);
      fprintf(stderr, " mitkor: %hi\n", iFile->dgb[i].mitkor);
      fprintf(stderr, " datofs: %hi\n", iFile->dgb[i].datofs);
      fprintf(stderr, " diradr: 0x%p\n", iFile->dgb[i].diradr);
      fprintf(stderr, " absadr: 0x%p\n", iFile->dgb[i].absadr);
      fprintf(stderr, " mpladr: 0x%p\n", iFile->dgb[i].mpladr);
      fprintf(stderr, " aztadr: 0x%p\n", iFile->dgb[i].aztadr);
      fprintf(stderr, " rztadr: 0x%p\n", iFile->dgb[i].rztadr);
      fprintf(stderr, " datadr: 0x%p\n", iFile->dgb[i].datadr);
    }

    fread(&iFile->dgb[i].refanf,    sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].thekx1,    sizeof(double),  1, fP);
    fread(&iFile->dgb[i].thekx2,    sizeof(double),  1, fP);
    fread(&iFile->dgb[i].polexp,    sizeof(double),  1, fP);
    fread(&iFile->dgb[i].thefac,    sizeof(double),  1, fP);
    fread(iFile->dgb[i].refzyl,     sizeof(char),   10, fP);
    fread(iFile->dgb[i].refsau,     sizeof(char),   10, fP);
    fread(iFile->dgb[i].refaus,     sizeof(char),   10, fP);
    fread(&iFile->dgb[i].aztlen,    sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].memhandle, sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].abslen,    sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].refzylmit, sizeof(short),   1, fP);
    fread(&iFile->dgb[i].refsaumit, sizeof(short),   1, fP);
    fread(&iFile->dgb[i].refausmit, sizeof(short),   1, fP);
    fread(&iFile->dgb[i].azttyp,    sizeof(short),   1, fP);
    fread(&iFile->dgb[i].aztofs,    sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].blkofs,    sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].blkdur,    sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].blknum,    sizeof(int32_t), 1, fP);
    fread(&iFile->dgb[i].synadr,                 4,  1, fP);
    fread(&iFile->dgb[i].synlen,    sizeof(int32_t), 1, fP);
    fread(iFile->dgb[i].reserve,    sizeof(char),   44, fP);
    
    iFile->dgb[i].refzyl[9]   = '\0';
    iFile->dgb[i].refsau[9]   = '\0';
    iFile->dgb[i].refaus[9]   = '\0';
    iFile->dgb[i].reserve[43] = '\0';

    if (debug) {
      fprintf(stderr, " refanf: %i\n", iFile->dgb[i].refanf);
      fprintf(stderr, " thekx1: %f\n",  iFile->dgb[i].thekx1);
      fprintf(stderr, " thekx2: %f\n",  iFile->dgb[i].thekx2);
      fprintf(stderr, " polexp: %f\n",  iFile->dgb[i].polexp);
      fprintf(stderr, " thefac: %f\n",  iFile->dgb[i].thefac);
      fprintf(stderr, " refzyl: %s\n",  iFile->dgb[i].refzyl);
      fprintf(stderr, " refsau: %s\n",  iFile->dgb[i].refsau);
      fprintf(stderr, " refaus: %s\n",  iFile->dgb[i].refaus);
      fprintf(stderr, " aztlen: %i\n", iFile->dgb[i].aztlen);
      fprintf(stderr, " memhandle: %i\n", iFile->dgb[i].memhandle);
      fprintf(stderr, " abslen: %i\n", iFile->dgb[i].abslen);
      fprintf(stderr, " refzylmit: %hi\n", iFile->dgb[i].refzylmit);
      fprintf(stderr, " refsaumit: %hi\n", iFile->dgb[i].refsaumit);
      fprintf(stderr, " refausmit: %hi\n", iFile->dgb[i].refausmit);
      fprintf(stderr, " azttyp: %hi\n", iFile->dgb[i].azttyp);
      fprintf(stderr, " aztofs: %i\n", iFile->dgb[i].aztofs);
      fprintf(stderr, " blkofs: %i\n", iFile->dgb[i].blkofs);
      fprintf(stderr, " blkdur: %i\n", iFile->dgb[i].blkdur);
      fprintf(stderr, " blknum: %i\n", iFile->dgb[i].blknum);
      fprintf(stderr, " synadr: 0x%p\n",  iFile->dgb[i].synadr);
      fprintf(stderr, " synlen: %i\n", iFile->dgb[i].synlen);
      fprintf(stderr, " reserve: %s\n",  iFile->dgb[i].reserve);
    }

    if (swap_end) {
      swap_endian(&iFile->dgb[i].datart,    sizeof(short));
      swap_endian(&iFile->dgb[i].absart,    sizeof(short));
      swap_endian(&iFile->dgb[i].dltphi,    sizeof(double));
      swap_endian(&iFile->dgb[i].dltzei,    sizeof(double));
      swap_endian(&iFile->dgb[i].kananz,    sizeof(short));
      swap_endian(&iFile->dgb[i].absanf,    sizeof(int32_t));
      swap_endian(&iFile->dgb[i].absend,    sizeof(int32_t));
      swap_endian(&iFile->dgb[i].fortyp,    sizeof(short));
      swap_endian(&iFile->dgb[i].beranz,    sizeof(short));
      swap_endian(&iFile->dgb[i].blkanz,    sizeof(short));
      swap_endian(&iFile->dgb[i].zykanz,    sizeof(int32_t));
      swap_endian(&iFile->dgb[i].zykanf,    sizeof(int32_t));
      swap_endian(&iFile->dgb[i].zykend,    sizeof(int32_t));
      swap_endian(&iFile->dgb[i].zyklng,    sizeof(int32_t));
      swap_endian(&iFile->dgb[i].zyksta,    sizeof(int32_t));
      swap_endian(&iFile->dgb[i].zyksum,    sizeof(int32_t));
      swap_endian(&iFile->dgb[i].datfor,    sizeof(short));
      swap_endian(&iFile->dgb[i].mitanz,    sizeof(short));
      swap_endian(&iFile->dgb[i].mitkor,    sizeof(short));
      swap_endian(&iFile->dgb[i].datofs,    sizeof(short));
      swap_endian(&iFile->dgb[i].diradr,    4);
      swap_endian(&iFile->dgb[i].absadr,    4);
      swap_endian(&iFile->dgb[i].mpladr,    4);
      swap_endian(&iFile->dgb[i].aztadr,    4);
      swap_endian(&iFile->dgb[i].rztadr,    4);
      swap_endian(&iFile->dgb[i].datadr,    4);
      swap_endian(&iFile->dgb[i].refanf,    sizeof(int32_t));
      swap_endian(&iFile->dgb[i].thekx1,    sizeof(double));
      swap_endian(&iFile->dgb[i].thekx2,    sizeof(double));
      swap_endian(&iFile->dgb[i].polexp,    sizeof(double));
      swap_endian(&iFile->dgb[i].thefac,    sizeof(double));
      swap_endian(&iFile->dgb[i].aztlen,    sizeof(int32_t));
      swap_endian(&iFile->dgb[i].memhandle, sizeof(int32_t));
      swap_endian(&iFile->dgb[i].abslen,    sizeof(int32_t));
      swap_endian(&iFile->dgb[i].refzylmit, sizeof(short));
      swap_endian(&iFile->dgb[i].refsaumit, sizeof(short));
      swap_endian(&iFile->dgb[i].refausmit, sizeof(short));
      swap_endian(&iFile->dgb[i].azttyp,    sizeof(short));
      swap_endian(&iFile->dgb[i].aztofs,    sizeof(int32_t));
      swap_endian(&iFile->dgb[i].blkofs,    sizeof(int32_t));
      swap_endian(&iFile->dgb[i].blkdur,    sizeof(int32_t));
      swap_endian(&iFile->dgb[i].blknum,    sizeof(int32_t));
      swap_endian(&iFile->dgb[i].synadr,    4);
      swap_endian(&iFile->dgb[i].synlen,    sizeof(int32_t));
    }
  } /* end loop over groups */
  endDGB = ftell(fP);

  unsigned int mpSize, aztlen;
  iFile->nC = 0;
  /* ----- read in the data groups ----- */
  for (i=0; i<iFile->apb.grpanz; i++) {
    /* pointer to current data group */
    iDGB = &iFile->dgb[i];

    /* only read in data that is CAD based for now */
/*     if (iDGB->datart != 1) */
/*       continue; */

    /* go to the address of the data directory */
    fseek(fP, (long)iDGB->diradr, SEEK_SET);

    iDGB->diradr = (iFileDAT*) malloc(iDGB->kananz*sizeof(iFileDAT));

    /* loop over channels */
    for (j=0; j<iDGB->kananz; j++) {
      pDI = &iDGB->diradr[j];

      fread(pDI->signam,   sizeof(char),   10, fP);
      fread(pDI->sigein,   sizeof(char),   10, fP);
      fread(pDI->sigkom,   sizeof(char),    6, fP);
      fread(&pDI->minimum, sizeof(float),   1, fP);
      fread(&pDI->maximum, sizeof(float),   1, fP);
      fread(&pDI->dstov,   sizeof(double),  1, fP);
      fread(&pDI->kstov,   sizeof(double),  1, fP);
      fread(&pDI->abskor,  sizeof(int32_t), 1, fP);
      fread(&pDI->desori,  sizeof(short),   1, fP);
      fread(&pDI->zykofs,  sizeof(int32_t), 1, fP);
      fread(&pDI->kaltyp,  sizeof(short),   1, fP);
      fread(&pDI->nultyp,  sizeof(short),   1, fP);
      fread(&pDI->sigtyp,  sizeof(short),   1, fP);
      fread(&pDI->bitnum,  sizeof(short),   1, fP);
      fread(&pDI->mitofs,  sizeof(int32_t), 1, fP);
      fread(&pDI->nulwrt,  sizeof(double),  1, fP);
      fread(pDI->kalfak,   sizeof(double),  2, fP);
      fread(pDI->intern,   sizeof(int32_t), 6, fP);
      fread(&pDI->refwrt,  sizeof(double),  1, fP);

      pDI->signam[9] = '\0';
      pDI->sigein[9] = '\0';
      pDI->sigkom[5] = '\0';

      if (swap_end) {
        swap_endian(&pDI->minimum, sizeof(float));
        swap_endian(&pDI->maximum, sizeof(float));
        swap_endian(&pDI->dstov,   sizeof(double));
        swap_endian(&pDI->kstov,   sizeof(double));
        swap_endian(&pDI->abskor,  sizeof(int32_t));
        swap_endian(&pDI->desori,  sizeof(short));
        swap_endian(&pDI->zykofs,  sizeof(int32_t));
        swap_endian(&pDI->kaltyp,  sizeof(short));
        swap_endian(&pDI->nultyp,  sizeof(short));
        swap_endian(&pDI->sigtyp,  sizeof(short));
        swap_endian(&pDI->bitnum,  sizeof(short));
        swap_endian(&pDI->mitofs,  sizeof(int32_t));
        swap_endian(&pDI->nulwrt,  sizeof(double));
        for (k=0; k<2; k++)
          swap_endian(pDI->kalfak, sizeof(double));
        for (k=0; k<6; k++)
          swap_endian(pDI->intern, sizeof(int32_t));
        swap_endian(&pDI->refwrt,  sizeof(double));
      }

      iFile->nC++;
    } /* end loop over channels */

    /* output the channel info */
    fprintf(stderr,"====================================================");
    fprintf(stderr,"===========================\n");
    fprintf(stderr,"| %20s Group %3i : %3i Channels (type %2i) %19s |\n", "",
            i, iDGB->kananz, iDGB->datart, "");
    fprintf(stderr,"+-----+------------+------------+------+------------");
    fprintf(stderr,"+---------------+---------+\n");
    fprintf(stderr,"| Chn | Name       | Units      | Type | Slope      ");
    fprintf(stderr,"| Offset        | Samples |\n");
		fprintf(stderr,"+-----+------------+------------+------+------------");
    fprintf(stderr,"+---------------+---------+\n");
    for (j=0; j<iDGB->kananz; j++) {
      pDI = &iDGB->diradr[j];
      /* write the channel info */
      fprintf(stderr, "");
      fprintf(stderr, "| %3u | %-10s | %-10s | %4u | %10.4f | %13.4f | %7i |\n",
              j, pDI->signam, pDI->sigein, pDI->sigtyp, pDI->kalfak[1],
              pDI->kalfak[0], pDI->zykofs);
    }
		fprintf(stderr,"+-----+------------+------------+------+------------");
    fprintf(stderr,"+---------------+---------+\n\n");

    /* ----- load the measurement table (required) ----- */
    fseek(fP, (long)iDGB->mpladr, SEEK_SET);
    mpSize = iDGB->kananz*iDGB->beranz;
    iDGB->mpladr = (iFileMPTAB*) malloc(mpSize*sizeof(iFileMPTAB));
    
    for (j=0; j<mpSize; j++) {
      fread(&iDGB->mpladr[j].beranf, sizeof(int32_t), 1, fP);
      fread(&iDGB->mpladr[j].wrtanz, sizeof(int32_t), 1, fP);
      fread(&iDGB->mpladr[j].absint, sizeof(int32_t), 1, fP);
      fread(&iDGB->mpladr[j].beradr, sizeof(int32_t), 1, fP);
      fread(&iDGB->mpladr[j].adrint, sizeof(int32_t), 1, fP);
    }

    if (swap_end) {
      swap_endian(&iDGB->mpladr[j].beranf, sizeof(int32_t));
      swap_endian(&iDGB->mpladr[j].wrtanz, sizeof(int32_t));
      swap_endian(&iDGB->mpladr[j].absint, sizeof(int32_t));
      swap_endian(&iDGB->mpladr[j].beradr, sizeof(int32_t));
      swap_endian(&iDGB->mpladr[j].adrint, sizeof(int32_t));
    }

    /* ----- load the Absolute Time Table (AZT) ----- */
    if (iDGB->aztadr > 0) {
      fseek(fP, (long)iDGB->aztadr, SEEK_SET);
      aztlen = (unsigned int)(iDGB->aztlen/4);

      iDGB->aztadr = (int*) malloc(aztlen*sizeof(int));
      
      fread(iDGB->aztadr, sizeof(int), aztlen, fP);
    } /* end AZT read */

    /* ----- load the Relative Time Table (RZT) ----- */
    if (iDGB->aztadr > 0) {
      fseek(fP, (long)iDGB->rztadr, SEEK_SET);

      iDGB->rztadr = (short*) calloc(iDGB->zykanz, sizeof(short));
      
      fread(iDGB->rztadr, sizeof(short), iDGB->zykanz, fP);
    } /* end AZT read */
         
  } /* end loop over groups for data read */
  endHeader = ftell(fP);

  fprintf(stderr, "Number of channels: %u\n", iFile->nC);
	fprintf(stderr, "End of APB filepos: %li\n", endAPB);
	fprintf(stderr, "End of DGB filepos: %li\n", endDGB);
	fprintf(stderr, "End of DI/MT/AZT/RZT of CA data: %li\n", endHeader);

  /* ----- Allocate/read the channel data ----- */
  int iCyl, cyc, iM, iC = 0;
  size_t nRead;
  double CA0, dCA;
	unsigned int mtStart;
  iFileMPTAB *pMT;
  iFileChannel *ch;
  FILE *fDbg;
  char outFile[64];
  iFile->cdat = (iFileChannel*) calloc(iFile->nC, sizeof(iFileChannel));

  /* set channel type to unknown */
  for (i=0; i<iFile->nC; i++)
    iFile->cdat[iC].type = 0;

  for (i=0; i<iFile->apb.grpanz; i++) {
    /* pointer to current data group */
    iDGB = &iFile->dgb[i];

    /* loop over channels for this data group */
    for (j=0; j<iDGB->kananz; j++) {
      /* current channel */
      ch = &iFile->cdat[iC];
      pDI = &iDGB->diradr[j];

      /* error checking */
      if (pDI == NULL) {
        fprintf(stderr, "Error: channel %i of group %i is invalid\n", j, i);
        exit(2);
      }

      ch->no        = iC;
      ch->type      = iDGB->datart;
      ch->axisType  = iDGB->datart;
      ch->cal_off   = pDI->kalfak[0];
      ch->cal_slope = pDI->kalfak[1];
      ch->tdc_off   = pDI->abskor*iDGB->dltphi;
      strcpy(ch->name,    pDI->signam);
      strcpy(ch->unit[0], pDI->sigein);
      strcpy(ch->unit[1], iDGB->absein);

      /* Get the cylinder number 
         Assumes that the last single digit number is cylinder number */
      len = strlen(ch->name);
      iCyl = 0;
      if (isdigit(ch->name[len-1]))
        iCyl = atoi(&ch->name[len-1]);

      if (iCyl > 0)
        ch->cyl = iCyl;

      /* determine the data type */
      switch (iDGB->datart) {
      case 1:  /* crank angle */
        break;
      case 2:  /* time */
        break;
      case 3:  /* realtime (RTP) */
        break;
      case 4:  /* results */
        break;
      case 5:  /* asynchronous UTC data */
        break;
      case 10: /* cycle*/
        break;
      default:
        fprintf(stderr, "Warning: channel data not recognized\n");
        fprintf(stderr, "  %i-%i: %s type %i\n",i,j,pDI->signam,iDGB->datart);
      } /* end switch datart */

      /* determine the signal type */
      /* if (pDI->sigtyp == 0) */
/*         ch->type = getIFileSigType(ch->name); */
/*       else  */
/*         ch->type = pDI->sigtyp; */

      switch (pDI->sigtyp) {
      case 1:  /* cylinder pressure */
        break;
      case 2:  /* needle lift signal */
        break;
      case 6:  /* dT ? (time/degree) */
        break;
      case 10:  /* injection pressure */
        break;
      case 11:  /* ignition angle */
        break;
      case 14:  /* knock boss acceleration */
        break;
      case 16:  /* intake manifold pressure */
        break;
      case 17:  /* exhaust manifold pressure */
        break;
      }

      /* configure the offset */

      /* get the measurement tables */
      fprintf(stderr,"\n+----------------------------------------------+\n");
      fprintf(stderr,"|       %8s (%3i) Measurement Tables      |\n",
              ch->name, iC);
      fprintf(stderr,"+-----+---------+---------+-------+------------+\n");
      fprintf(stderr,"| No  | Start   | Samples | Res   | Data Start |\n");
      fprintf(stderr,"+-----+---------+---------+-------+------------+\n");
      for (iM=0; iM<iDGB->beranz; iM++) {
        pMT = &iDGB->mpladr[iM];
      } /* end loop over measurement tables */

      ch->nCyc  = iDGB->zykanz;
      ch->nData = ch->nMT = 0;
      for (iM=0; iM<iDGB->beranz; iM++) {
        pMT = &iDGB->mpladr[iM];
        fprintf(stderr, "| %2u  |", iM);

        if (pMT->wrtanz > 0 && pMT->adrint > 0 && pMT->absint > 0) {
          ch->nData = pMT->wrtanz;

          /* set up the data axis */
          ch->axis    = (double*) calloc(ch->nData, sizeof(double));
          ch->dx      = pMT->absint;
          ch->axis[0] = pMT->beranf;
          if (ch->type == 1) {
            ch->dx      *= iDGB->dltphi;
            ch->axis[0] *= iDGB->dltphi;
            fprintf(stderr," %7.2f | %7d | %5.2f |",ch->axis[0],ch->nData,
                    ch->dx);

          }
          else
            fprintf(stderr," %7d | %7d | %5d |",ch->axis[0],ch->nData,
                    ch->dx);
          fprintf(stderr, " %10d |\n", pMT->beradr);
          mtStart = pMT->beradr;
          
          for (k=1; k<ch->nData; k++)
            ch->axis[k] = ch->axis[k-1] + ch->dx;
          ch->nMT++;
        }
        else
          fprintf(stderr, " %7d | %7d | %5d | %10d | (Invalid)\n", pMT->beranf,
                  pMT->wrtanz, pMT->absint, pMT->beradr);
      }
      fprintf(stderr,"+----------------------------------------------+\n");
      fprintf(stderr, "Number of valid measurement tables: %i\n", ch->nMT);

      /* load the data */
      fprintf(stderr, "loading data %i %i\n", ch->nCyc, ch->nData);
      int nDat = ch->nCyc*ch->nData;
      float *dataBuf = (float*) calloc(ch->nCyc*ch->nData, sizeof(float));
      short *shortBuf = (short*) calloc(ch->nCyc*ch->nData, sizeof(short));
      float *convBuf;
      short *shortP = NULL;
      char *rawDat = NULL;
      int vDat, iP;
      double correctDat;
      ch->data = (double**) calloc(ch->nCyc, sizeof(double*));
      for (k=0; k<ch->nCyc; k++)
        ch->data[k] = (double*) calloc(ch->nData, sizeof(double));

      fprintf(stderr, "short size %i\n", sizeof(short));

      if (pDI->zykofs > ch->nData) { /* channel data is interleaved */
        fprintf(stderr, "data is interleaved %i %i %i %i %i\n", pDI->zykofs, ch->nData, iDGB->datfor, iDGB->fortyp, mtStart);
        for (cyc=0; cyc<ch->nCyc; cyc++) {
          fseek(fP,mtStart+cyc*pDI->zykofs*iDGB->datfor,SEEK_SET);
/*           if (iDGB->datfor == 2) */
/*           nRead = fread((void*)((char*)dataBuf+cyc*ch->nData*2), */
/*                         2,ch->nData,fP); */
            nRead = fread(&shortBuf[cyc*ch->nData],iDGB->datfor,ch->nData,fP);
/*           else if (iDGB->datfor == 4) */
/*             nRead = fread(&dataBuf[cyc*ch->nData],iDGB->datfor,ch->nData,fP); */

          if (nRead != ch->nData)
            fprintf(stderr, "error: data read not correct\n");
        } /* end loop over cycles */


      }
      else { /* channel data is sequential */
 /*        fprintf(stderr, "data is sequential %i %i\n", pDI->zykofs, ch->nData); */
/*         for (cyc=0; cyc<ch->nCyc; cyc++) { */
/*           fseek(fP,mtStart+cyc*ch->nData,SEEK_SET); */
/*           nRead = fread(&ch->data[cyc],iDGB->datfor,ch->nData,fP); */
/*           if (nRead != ch->nData) */
/*             fprintf(stderr, "error: data read not correct\n"); */
/*         } */
      }

      
      /* swap endian if necessary */
      /* if (iDGB->datfor == 2) */
/*         rawDat = (char*)shortBuf; */
/*       else if (iDGB->datfor == 4) */
/*         rawDat = (char*)dataBuf; */
      rawDat = (char*)shortBuf;

        if (swap_end) {
          fprintf(stderr, "swapping endian\n");
          for (iP=0; iP<nDat; iP++) {
            swap_endian(&rawDat[iP*2], iDGB->datfor);
          }
        }
/*       if (swap_end) { */
/*         fprintf(stderr, "swap endian data\n"); */
/*         for (iP=0; iP<nDat; iP++) { */
/*           shortP = (short*)(rawDat+iP*d); */
/*           swap_endian(shortP, iDGB->datfor); */
/*         } */
/*       } */

      int dp;
      if (iDGB->fortyp == 0) {
        convBuf = (float*) calloc(nDat, sizeof(float));
        fprintf(stderr, "converting %i-bit fixed point data\n", 8*iDGB->datfor);
        for (iP=0; iP<nDat; iP++) {
          dp = nDat - iP - 1;
          shortP = (short*)(rawDat+dp*2);
/*           vDat = (int)(rawDat+iP*iDGB->datfor); */
/*           convBuf[iP] = (float)(vDat); */
          vDat = (int)(*shortP);
          convBuf[dp] = (float)(vDat);
        }
      }
      else 
        convBuf = dataBuf;
            
      /* apply the offset correction */
      for (cyc=0; cyc<ch->nCyc; cyc++)
        for (k=0; k<ch->nData; k++)
/*           ch->data[cyc][k] = convBuf[cyc*ch->nData+k]; */
          ch->data[cyc][k] = ch->cal_slope*convBuf[cyc*ch->nData+k]+ch->cal_off;

      if (ch->type == 1) { /* write pressure data */
        sprintf(outFile, "debug_%s.dat", ch->name);
        fDbg = fopen(outFile, "w");

        for (iP=0; iP<nDat; iP++)
          fprintf(fDbg, "%15.6e\n", convBuf[iP]);
        fprintf(fDbg, "\n\n");

        for (k=0; k<ch->nData; k++) {
          fprintf(fDbg, "%15.6e ", ch->axis[k]);
          for (cyc=0; cyc<ch->nCyc; cyc++)
            fprintf(fDbg, "%15.6e ", ch->data[cyc][k]);
          fprintf(fDbg, "\n");
        }
        
        fclose(fDbg);
      }

      iC++;
    } /* end loop over channels */

  } /* end loop over data groups */

  fclose(fP);
  return iFile;
}
/* -------------------------------------------------------------------------- */

int getIFileSigType(char *sigName)
{
  int i;
  char name[64];
  /* convert to all lower case */
  strcpy(name, sigName);
  for (i=0; i<strlen(name); i++) name[i] = tolower(name[i]);
  
  if (!strncmp(name, "cylpr", 5) ||
      !strncmp(name, "pcyl", 4)  || !strncmp(name, "pzyl", 4)  ||
      !strncmp(name, "p_cyl", 5) || !strncmp(name, "p_zyl", 5))
    return 1;
  else if (!strncmp(name, "nlift", 5) || !strncmp(name, "needle", 6))
    return 2;
  else if (!strncmp(name, "plin", 4) || !strncmp(name, "inj", 3))
    return 10;
  else if (!strcmp(name, "spark") || !strcmp(name, "sprk") ||
           !strcmp(name, "spk") || !strncmp(name, "ign", 3) )
    return 11;
  else if (!strncmp(name, "pin", 3) || !strncmp(name, "pman", 4))
    return 16;
  else if (!strncmp(name, "pex", 3))
    return 17;
}
/* -------------------------------------------------------------------------- */

int  iFileReadExtendedOPs(FILE *fP, uint32_t pExt, iFileOP **op)
{
  if (pExt == 0) {
    fprintf(stderr, "no extended operating parameters\n");
    return 0;
  }
  
  int i, j, k, nK, nL, iK, nExOp, iOp;
  int32_t oldpos, pD;
  char *ssKey, **longNames, *paramFile;

  uint32_t nNumberOfObjects, ulKeyLength, ulObjSize, ulNameLength;
  int32_t eType, nFirstObject, nLastObject;
  int32_t nextExtObj, keyFilePos, ulObjFilePos, ulUnitLength;

  /* store file positions */
  j = k = 1;
  oldpos = ftell(fP);
  fseek(fP, pExt, SEEK_SET); 
  pD = ftell(fP);

  /* read in the info about the data objects */
  fread(&nNumberOfObjects, sizeof(uint32_t), 1, fP);
  fread(&nFirstObject,     sizeof(int32_t),  1, fP);
  fread(&nLastObject,      sizeof(int32_t),  1, fP);

  nextExtObj = nFirstObject;

  /* local memory */
  ssKey     = (char*)  calloc(1, sizeof(char));
  longNames = (char**) calloc(1, sizeof(char*));
  paramFile = (char*)  calloc(1, sizeof(char));

  /* loop over extended objects */
  nExOp = 0;
  nL = nK = 0;
  for (i=0; i<nNumberOfObjects; i++) {
    fseek(fP, pExt + nextExtObj, SEEK_SET);

    /* read info about this object */
    fread(&eType,        sizeof(int32_t),  1, fP);
    fread(&ulKeyLength,  sizeof(uint32_t), 1, fP);
    fread(&keyFilePos,   sizeof(int32_t),  1, fP);
    fread(&ulObjSize,    sizeof(uint32_t), 1, fP);
    fread(&ulObjFilePos, sizeof(int32_t),  1, fP);
    fread(&nextExtObj,   sizeof(int32_t),  1, fP);

    /* read in the key */
    fseek(fP, pExt+keyFilePos, SEEK_SET);
    ssKey = (char*) realloc(ssKey, ulKeyLength*sizeof(char));
    fread(ssKey, sizeof(char), ulKeyLength, fP);

    /* move to object */
    fseek(fP, pD + ulObjFilePos, SEEK_SET);
    switch (eType) {
    case 1: { /* RTP type */
      break;
    }
    case 2: { /* eensor sensitivity, TDC value */
      break;
    }
    case 3: { /* SourceSignal, AmplifierType, Filter, DriftCompensation */
              /* OutputOffset, AmplifierRangeUnit, SensorCalibrationDate */
      break;
    }
    case 5: { /* extended signal name */
      fread(&ulNameLength, sizeof(uint32_t), 1, fP);
      nL++;
      longNames       = (char**) realloc(longNames,    nL*sizeof(char*));
      longNames[nL-1] = (char*)   calloc(ulNameLength,    sizeof(char));
      fread(longNames[nL-1], sizeof(char), ulNameLength, fP);
      /* fprintf(stderr, "%5i - Long Name: %s\n", i, longNames[nL-1]); */
    }
    case 7: { /* YRangeMin */
      break;
    }
    case 8: { /* YRangeMax */
      break;
    }
    case 15: { /* parameter file */
      paramFile = (char*) realloc(paramFile, ulObjSize*sizeof(char));
      fread(paramFile, sizeof(char), ulObjSize, fP);
      /* fprintf(stderr, "%5i - Embedded parameter file: %s\n", i, paramFile); */
      break;
    }
    case 16: { /* CY' SYNC' */
      break;
    }
    case 17: { /* SYS' iFileName software firmware */
      break;
    }
    case 19: { /* extended parameter */
      nK++;
      nExOp++;
      iK = nK-1+28;
      op[0] = (iFileOP*) realloc(op[0], (nK+28)*sizeof(iFileOP));

      /* read in the OP name */
      fread(&ulNameLength, sizeof(uint32_t), 1, fP);
      fread(op[0][iK].name, sizeof(char), ulNameLength, fP);

      /* read in the length of the units */
      fread(&ulUnitLength, sizeof(int32_t), 1, fP);
      if (ulUnitLength > 1)
        fread(op[0][iK].unit, sizeof(char), ulUnitLength, fP);
      else
        strcpy(op[0][iK].unit, "");

      /* read in the value */
      fread(&op[0][iK].value, sizeof(double), 1, fP);

      /* fprintf(stderr, "%5i - %s --> %s: %f [%s]\n", i, ssKey, op[0][iK].name, */
/*               op[0][iK].value, op[0][iK].unit); */
      break;
    }
    case 22: { /* SYS' options licence hardware serial no. */
      break;
    }
    default: {
      fprintf(stderr, "Warning: eType %i not supported for %s\n", eType, ssKey);
      break;
    }
    } /* end eType switch */

  } /* end loop over number of objects */

  for (i=0; i<nL; i++) free(longNames[i]);
  free(longNames);
  free(ssKey); free(paramFile);

  /* go back to old file position */
  fseek(fP, oldpos, SEEK_SET);
  fprintf(stderr, "finished reading extended objects\n");
  return nExOp;
}
/* -------------------------------------------------------------------------- */

char *readFile2ptr(char *fileName, int32_t *fileLength)
{
  char *f;
  int32_t lengthRead;

  /* open the file for reading */
  FILE *pFile = fopen(fileName, "rb");
  if ( !(pFile = fopen(fileName,"rb")) ) {
    fprintf(stderr, "\nError on opening file\n");
    fprintf(stderr, "File is already open or does not exist\n");
    exit(2);
  }

  /* determine the fileLength of the file */
  fseek(pFile, 0L, SEEK_END);
  *fileLength = ftell(pFile);
  fseek(pFile, 0L, SEEK_SET);

  if (*fileLength == -1) {
    fprintf(stderr, "\n File is empty \n");
  }

  /* allocate memory for the whole file */
  f = (char*) calloc(*fileLength, sizeof(char));

  /* read in the file */
  lengthRead = fread(f, sizeof(char), *fileLength, pFile);

  if ( *fileLength != lengthRead ) {
    fprintf(stderr, "File-Length = %i\n", *fileLength);
    fprintf(stderr, "\nEndOfFile encountered while reading file!\n");
  }

  fclose(pFile);
  return(f);
}
/* -------------------------------------------------------------------------- */

void swap_endian(void *input, int n)
{
  int i,j;
  char *src = (char*)input;
  char *dst = (char*) calloc(n, sizeof(char));

  for (i=0, j=n-1; i<n; i++, j--)
    dst[i] = src[j];
  for (i=0; i<n; i++) src[i] = dst[i];
  free(dst);
}
/* -------------------------------------------------------------------------- */

double double_swap(double d)
{
  int i, j;
  double a;
  unsigned char *dst = (unsigned char *)&a;
  unsigned char *src = (unsigned char *)&d;

  for (i=0, j=7; i<8; i++, j--)
    dst[i] = src[j];

  return a;
}
/* -------------------------------------------------------------------------- */

long long_swap(long d)
{
  int i, j;
  long a;
  unsigned char *dst = (unsigned char *)&a;
  unsigned char *src = (unsigned char *)&d;

  for (i=0, j=3; i<4; i++, j--)
    dst[i] = src[j];

  return a;
}
/* -------------------------------------------------------------------------- */

short short_swap(short d)
{
  int i, j;
  short a;
  unsigned char *dst = (unsigned char *)&a;
  unsigned char *src = (unsigned char *)&d;

  for (i=0, j=1; i<2; i++, j--)
    dst[i] = src[j];

  return a;
}
/* -------------------------------------------------------------------------- */


//! Byte swap unsigned short
/* uint16_t swap_uint16( uint16_t val )  */
/* { */
/*     return (val << 8) | (val >> 8 ); */
/* } */

/* //! Byte swap short */
/* int16_t swap_int16( int16_t val )  */
/* { */
/*     return (val << 8) | ((val >> 8) & 0xFF); */
/* } */

/* //! Byte swap unsigned int */
/* uint32_t swap_uint32( uint32_t val ) */
/* { */
/*     val = ((val << 8) & 0xFF00FF00 ) | ((val >> 8) & 0xFF00FF );  */
/*     return (val << 16) | (val >> 16); */
/* } */

/* //! Byte swap int */
/* int32_t swap_int32( int32_t val ) */
/* { */
/*     val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF );  */
/*     return (val << 16) | ((val >> 16) & 0xFFFF); */
/* } */

/* void ChangeMemEndianness(uint64_t *mem, size_t size)  */
/* { */
/* uint64_t m1 = 0xFF00FF00FF00FF00ULL, m2 = m1 >> CHAR_BIT; */

/* size = (size + (sizeof (uint64_t) - 1)) / sizeof (uint64_t); */
/* for(; size; size--, mem++) */
/*   *mem = ((*mem & m1) >> CHAR_BIT) | ((*mem & m2) << CHAR_BIT); */
/* } */
