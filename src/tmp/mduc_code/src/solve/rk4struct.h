#ifndef RK_STRUCT_H
#define RK_STRUCT_H

struct RKSTRUCT {
  double *k1;
  double *k2;
  double *k3;
  double *k4;
  
  double *Y;
  double *Ytmp;
  
  double dt;
  double dt_max;
  double t_end;
  int    step;
  int    iCrit;
  int    jCrit;
};
typedef struct RKSTRUCT rkStruct;

#endif /* RK_STRUCT_H */
