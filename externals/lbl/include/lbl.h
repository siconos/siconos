/*
 * A generic interface to the MA27 and MA57 wrappers.
 *
 * 10 Sept 08: First version.
 *             Michael P. Friedlander
 *             Dominique Orban
 *
 */

#ifndef _LBL_LIB_H
#define _LBL_LIB_H

#include "ma27.h"
#include "ma57.h"

#define LBL_I_PIV_SELECTION    1
#define LBL_I_PIV_NUMERICAL    2
#define LBL_I_SCALING          3
#define LBL_D_PIV_THRESH       4
#define LBL_D_PIV_NONZERO      5

typedef struct LBL_Data {
    /* LBL problem object.*/
    /* ----------------------------------------------------------------- */

    Ma27_Data *ma27;
    /* Pointer to an MA27 problem context. */

    Ma57_Data *ma57;
    /* Pointer to an MA57 problem context. */

    int lblsolver;
    /* Specifies which solver is being used for the current problem context. */

    int *irn, *jcn;
    /* Convenience pointers to sparsity patters that are stored in ma27/ma57. */

} LBL_Data;

/* Interfaces to the LBL routines. */

LBL_Data  * LBL_Initialize( int nz, int n, FILE *outfile, int lblsolver );
int LBL_Analyze( LBL_Data *data, int iflag );
int LBL_Factorize(LBL_Data *data, double A[] );
int LBL_Solve( LBL_Data *data, double x[] );
int LBL_Refine( LBL_Data *data, double x[], double rhs[],
                double A[], double tol, int maxitref, int job );
void LBL_Finalize( LBL_Data *data );
void   LBL_set_int_parm( LBL_Data *data, int parm, int val );
void   LBL_set_real_parm( LBL_Data *data, int parm, double val );
int    LBL_get_int_parm( LBL_Data *data, int parm );
double LBL_get_real_parm( LBL_Data *data, int parm );

void *LBL_Calloc(int nmemb, size_t size);
void *LBL_Realloc(void *ptr, int nmemb, size_t size);
void LBL_Free_Object(void **ptr);
#define LBL_Free(obj) LBL_Free_Object( (void**)(&(obj)) )

/* Errors and warnings */

#define ALLOCATION_ERROR   -1
#define REALLOCATION_ERROR -2
#define DEALLOCATION_ERROR -3

#define SETERRQ(n,s) {                                                \
  fprintf( stderr, "\n" );                                            \
  fprintf( stderr, "  LBL Error::        Code   = %d\n", n  );        \
  fprintf( stderr, "                     Msg   :: %s\n", s );         \
  fprintf( stderr, "  Error occured in   function %s\n", __FUNCT__ ); \
  fprintf( stderr, "                     file     %s\n", __FILE__ );  \
  fprintf( stderr, "                     line     %d\n", __LINE__ );  \
  exit( n );                                                        \
}

#define SETERRQi(n,s,i) {                                             \
  fprintf( stderr, "\n" );                                            \
  fprintf( stderr, "  LBL Error::        Code   = %d\n", n );         \
  fprintf( stderr, "                     Msg   :: %s\n", s );         \
  fprintf( stderr, "                     Value  : %d\n", i );         \
  fprintf( stderr, "  Error occured in   function %s\n", __FUNCT__ ); \
  fprintf( stderr, "                     file     %s\n", __FILE__ );  \
  fprintf( stderr, "                     line     %d\n", __LINE__ );  \
  exit( n );                                                        \
}

#define SETERRQg(n,s,x) {                                             \
  fprintf( stderr, "\n" );                                            \
  fprintf( stderr, "  LBL Error::        Code   = %d\n", n );         \
  fprintf( stderr, "                     Msg   :: %s\n", s );         \
  fprintf( stderr, "                     Value  : %g\n", x );         \
  fprintf( stderr, "  Error occured in   function %s\n", __FUNCT__ ); \
  fprintf( stderr, "                     file     %s\n", __FILE__ );  \
  fprintf( stderr, "                     line     %d\n", __LINE__ );  \
  exit( n );                                                        \
}

#define SETWARNQ(n,s) {                                               \
  fprintf( stderr, "\n" );                                            \
  fprintf( stderr, "  LBL Warning::      Code   = %d\n", n );         \
  fprintf( stderr, "                     Msg   :: %s\n", s );         \
  fprintf( stderr, "  Warning occured in function %s\n", __FUNCT__ ); \
  fprintf( stderr, "                     file     %s\n", __FILE__ );  \
  fprintf( stderr, "                     line     %d\n", __LINE__ );  \
}

#define SETWARNQi(n,s,i) {                                            \
  fprintf( stderr, "\n" );                                            \
  fprintf( stderr, "  LBL Warning::      Code   = %d\n", n );         \
  fprintf( stderr, "                     Msg   :: %s\n", s );         \
  fprintf( stderr, "                     Value  : %d\n", i );         \
  fprintf( stderr, "  Error occured in   function %s\n", __FUNCT__ ); \
  fprintf( stderr, "                     file     %s\n", __FILE__ );  \
  fprintf( stderr, "                     line     %d\n", __LINE__ );  \
}

#define SETWARNQg(n,s,x) {                                            \
  fprintf( stderr, "\n" );                                            \
  fprintf( stderr, "  LBL Warning::      Code   = %d\n", n );         \
  fprintf( stderr, "                     Msg   :: %s\n", s );         \
  fprintf( stderr, "                     Value  : %g\n", x );         \
  fprintf( stderr, "  Error occured in   function %s\n", __FUNCT__ ); \
  fprintf( stderr, "                     file     %s\n", __FILE__ );  \
  fprintf( stderr, "                     line     %d\n", __LINE__ );  \
}

#define SETWARNQs(n,s1,s2) {                                          \
  fprintf( stderr, "\n" );                                            \
  fprintf( stderr, "  LBL Warning::      Code   = %d\n", n );         \
  fprintf( stderr, "                     Msg   :: %s\n", s1 );        \
  fprintf( stderr, "                     Value  : %s\n", s2 );        \
  fprintf( stderr, "  Error occured in   function %s\n", __FUNCT__ ); \
  fprintf( stderr, "                     file     %s\n", __FILE__ );  \
  fprintf( stderr, "                     line     %d\n", __LINE__ );  \
}

#endif /* _LBL_LIB_H */
