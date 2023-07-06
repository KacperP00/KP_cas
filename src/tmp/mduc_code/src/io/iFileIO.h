#ifndef IFILE_READ_H
#define IFILE_READ_H

#include <stdint.h>

/* info for operating parameter */
typedef struct {
  char       name[10];
  char       unit[10];
  double     value;
} iFileOP;  /* 28 Byte */

/* struct for header information */
typedef struct {
  short      parlng;      /* header size in bytes */
  short      grpanz;      /* number of data groups in file (Gruppenanzahl) */
  char       filkom[80];  /* file comment */
  char       mesdat[24];  /* measurement date (Messdate)*/
  char       parfil[18];  /* measurement file used */
  char       pronam[18];
  char       prfstd[10];
  char       motnam[18];  /* engine name (Motorname) */
  short      mottyp;      /* engine type */
  short      tktanz;      /* strokes per cycle (0: 4 stroke, 1: 2 stroke */
  short      geoein;      /* geometry measurement basis (0: SI, 1: inches )*/
  double     hublng;      /* stroke [mm]         (Hublaenge) */
  double     pleuel;      /* rod length [mm]     (Pleullaenge) */
  double     bohrng;      /* bore [mm]           (Bohrung) */
  double     kompre;      /* compression ratio */
  double     desaxi;      /* axial pin offset [mm] */
  int32_t    grpflg;
  int32_t    datahandle;
  short      version;     /* iFile version */
  double     desaxi2;     /* axial pin offset 2 [mm] */
  uint32_t   pExt;        /* IFile Datenobjekte IndiCom 1.3 */
  short      nOp;         /* number of operating points */
  iFileOP    *op;         /* default operating parameters */
} iFileAPB;

typedef struct {
  char       signam[10];  /* signal name */
  char       sigein[10];  /* signal units */
  char       sigkom[6];
  float      minimum;
  float      maximum;
  double     dstov;
  double     kstov;
  int32_t    abskor;
  short      desori;
  int32_t    zykofs;      /* number of samples */
  short      kaltyp;
  short      nultyp;
  short      sigtyp;      /* signal type */
  short      bitnum;
  int32_t    mitofs;      /* Adressoffset fuer Var-Funktion */
  double     nulwrt;
  double     kalfak[2];   /* calibration curve ([0] offset, [1] slope) */
  int32_t    intern[6];
  double     refwrt;
} iFileDAT;

typedef struct {
  int32_t    beranf;      /* Anzahl der Messbereiche mit konst. Messauflosung */
  int32_t    wrtanz;      /* Anzahl der Messwerte */
  int32_t    absint;      /* Messaufloesung */
  int32_t    beradr;
  int32_t    adrint;      /* address interval */
} iFileMPTAB;

typedef struct {
  int        no;          /* channel number */
  int        cyl;         /* cylinder number */
  char       name[64];
  char       unit[2][10];
  int        signal;      /* signal type */
  int        filter;
  char       filterType[16];

  /* */
  int        type;
  int        axisType;
  int        nCyc;        /* number of cycles */
  int        nMT;
  int        nData;       /* number of data points per cycle */
  double     **data;
  double     *axis;
  double     dx;
  double     cal_off;     /* offset of calibration curve */
  double     cal_slope;   /* slope of calibration curve */
  double     tdc_off;
  
  /* statistics */
  double     min,max;
  double     *mean;
  double     *var;
} iFileChannel;

typedef struct {
  short      datart;      /* data type   1: Crank Angle 
                                         2: Timebase
                                         3: realtime results 
                                         4: results
                                         5: asynchronous UTC data */
  short      absart;      /* abscissa type (Abszissenart) */
  char       absein[10];  /* abscissa unit (Abszisseneinheit) */
  double     dltphi;      /* delta angle [CAD] */
  double     dltzei;      /* dt [ms] (delZeit) */
  short      kananz;      /* number of channels (Kanalanzahl) */
  int32_t    absanf;
  int32_t    absend;
  short      fortyp;      /* data type */
  short      beranz;      /* number of measurements per cycle */
  short      blkanz;
  int32_t    zykanz;      /* number of cycles */
  int32_t    zykanf;      /* first cycle */
  int32_t    zykend;      /* last cycle */
  int32_t    zyklng;      /* cycle length */
  int32_t    zyksta;
  int32_t    zyksum;
  short      datfor;      /* data length */
  short      mitanz;
  short      mitkor;
  short      datofs;
  iFileDAT   *diradr;     /* address of the data directory */
  uint8_t    *absadr;     /* address of the absolute time */
  iFileMPTAB *mpladr;     /* address of the measurement table */
  int32_t    *aztadr;     /* address of AZT (Absolutzeittabelle) */
  int16_t    *rztadr;     /* address of RZT (Relativzeittabelle) */
  uint8_t    *datadr;
  int32_t    refanf;
  double     thekx1;      /* Angle 1 */
  double     thekx2;      /* Angle 2 */
  double     polexp;      /* polytropic exponent */
  double     thefac;
  char       refzyl[10];  /* zero level channel */
  char       refsau[10];  /* intake channel */
  char       refaus[10];  /* exhaust channel */
  int32_t    aztlen;      /* Laenge der AZT */
  int32_t    memhandle;   /* Memory handle fuer DGR auf PC-Version */
  int32_t    abslen;      /* abscissa length (Laenge der Abszissentabelle) */
  short      refzylmit;
  short      refsaumit;
  short      refausmit;
  short      azttyp;
  int32_t    aztofs;
  int32_t    blkofs;
  int32_t    blkdur;
  int32_t    blknum;
  char       *synadr;     /* address of CDM? */
  int32_t    synlen;      /* length of CDM? */
  char       reserve[44];
} iFileDATGRP;

typedef struct {
  int32_t    eType;
  uint32_t   ulKeyLength;
  int32_t    lKeyFilePos;
  uint32_t   ulObjectSize;
  int32_t    lObjectFilePos;
  int32_t    lNextExtensionObjectFilePos;
} iFileExtensionObject;

typedef struct {
  uint32_t   ulNameLength;
  char       *pcName;
  int32_t    ulUnitLenght;
  char       *pcUnit;
  double     dValue;
} iFileExtensionOP;
typedef struct {
  uint32_t   lNumberOfObjects;
  int32_t    lFirstExtensionObjectFilePos;
  int32_t    lLastExtensionObjectFilePos;
} iFileExtension;

struct IFILESTRUCT {
  iFileAPB     apb;   /* Allgemeinen Parameterblock / General paramters */
  iFileDATGRP  *dgb;  /* Datengruppenbeschreibung / Data description */
  iFileDAT     di;    /* Data struct */

  /* channel data */
  int          nC;    /* number of channels */
  iFileChannel *cdat; /* channel data */
  iFileChannel **Pcyl; /* pointer to cylinder pressure */
  iFileChannel *Pint; /* pointer to intake manifold pressure */
  iFileChannel *Pexh; /* pointer to exhaust manifold pressure */
  iFileChannel *L_in; /* pointer to intake lift */
  iFileChannel *L_ex; /* pointer to exhaust lift */
};
/* -------------------------------------------------------------------------- */
typedef struct IFILESTRUCT iFileStruct;

iFileStruct * iFileRead(char *fileName);

char *readFile2ptr(char *fileName, int32_t *fileLength);
int  iFileReadExtendedOPs(FILE *fP, uint32_t pExt, iFileOP **op);

int getIFileSigType(char *sigName);
void   swap_endian(void *input, int n);
short  short_swap(short d);
long   long_swap(long d);
double double_swap(double d);

#endif /* IFILE_READ_H */
