#ifndef MDUC_POLY_INTERSECT_H
#define MDUC_POLY_INTERSECT_H

typedef struct{double x; double y;} gPoint;
typedef double real; // could be double.

real polyIntersect(int na, gPoint *a, int nb, gPoint *b);

#endif /* MDUC_POLY_INTERSECT_H */
