#include "commonlib.h"  // for REAL
#include "sparselib.h"  // for sparseMatrix, sparseVector


#ifdef __cplusplus
extern "C" {
#endif

void LUmod( int mode, int n, int krow, int kcol,
            sparseMatrix *L, sparseMatrix *U,
            REAL *y, REAL *z, REAL *w );
void Lprod ( int mode, int n, sparseMatrix *L, REAL *y, REAL *z );
void Usolve ( int mode, int n, sparseMatrix *U, REAL *y );
void LUforw ( int first, int last, int n, int nu, 
              REAL eps, sparseMatrix *L, sparseMatrix *U, REAL *y );
void LUback ( int first, int *last, int n, int nu,
              REAL eps, sparseMatrix *L, sparseMatrix *U, REAL *y, REAL *z );

void elm1 ( int first, int last, sparseVector *x, REAL *y, REAL cs, REAL sn );
void elm2 ( int first, int last, REAL *x, sparseVector *y, REAL cs, REAL sn );
void elm3 ( int first, int last, sparseVector *x, sparseVector *y, REAL cs, REAL sn );
void elmgen ( REAL *x, REAL *y, REAL eps, REAL *cs, REAL *sn );

#ifdef __cplusplus
}
#endif
