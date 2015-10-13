

void LUmod_dense ( int mode, int maxmod, int n, int krow, int kcol,
             double *L, double *U, double *y, double *z, double *w );
void Lprod_dense ( int mode, int maxmod, int n,
             double *L, double *y, double *z );
void LUforw_dense ( int first, int last, int n, int nu, int maxmod,
              double eps, double *L, double *U, double *y );
void LUback_dense ( int first, int *last, int n, int nu,
              int maxmod, double eps,
              double *L, double *U, double *y, double *z );
void Usolve_dense ( int mode, int maxmod, int n,
              double *U, double *y );
void Lprod_dense ( int mode, int maxmod, int n,
             double *L, double *y, double *z );
void elm_dense ( int first, int last, double *x, double *y, double cs, double sn );
void elmgen_dense ( double *x, double *y, double eps, double *cs, double *sn );

