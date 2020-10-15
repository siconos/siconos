/*
 * ============================================
 * A wrapper around the Harwell subroutine MA27
 * to dynamically allocate and ajust workspace.
 *
 * D. Orban            Montreal, September 2003
 * ============================================
 */

#ifndef _MA27_H
#define _MA27_H

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define imax(a,b) ((a) > (b) ? (a) : (b))
#define imin(a,b) ((a) < (b) ? (a) : (b))
#define ZERO       0.0

#define PRINT(...) if (ma27->outfile) fprintf(ma27->outfile, __VA_ARGS__)

#ifdef _AIX
#define FUNDERSCORE(a) a
#else
#define FUNDERSCORE(a) a##_
#endif

#define MA27ID        FUNDERSCORE(ma27id)
#define MA27AD        FUNDERSCORE(ma27ad)
#define MA27BD        FUNDERSCORE(ma27bd)
#define MA27CD        FUNDERSCORE(ma27cd)
#define MA27FACTORS   FUNDERSCORE(ma27factors)
#define MA27QDEMASC   FUNDERSCORE(ma27qdemasc)

typedef struct Ma27_Data {
    int     n, nz;               /* Order and #nonzeros */
    int     icntl[31], info[21];
    double  cntl[6];
    int    *irn, *jcn;           /* Sparsity pattern    */
    int    *iw, liw;             /* Integer workspace   */
    int    *ikeep;               /* Pivot sequence      */
    int    *iw1;                 /* Integer workspace   */
    int     nsteps;
    int     iflag;               /* Pivot selection     */
    double  ops;                 /* Operation count     */

    int     la;
    double *factors;             /* Matrix factors      */
    int     maxfrt;
    double *w;                   /* Real workspace      */

    double *residual;            /* = b - Ax            */

    char    fetched;             /* Factors have been fetched
                                  * Used for de-allocation */

    FILE   *outfile;             /* File for log output */
} Ma27_Data;

/* Below I indicate arrays of fixed length by specifying it
 * explicitly, e.g. icntl[30], arrays of variable length by
 * specifying it implicitly, e.g. iw[]. The remaining variables
 * are either variables of the type indicated, or arrays whose
 * size is fixed, but depends on the value of other parameters.
 */

extern void MA27ID( int icntl[30], double cntl[5] );
extern void MA27AD( int *n, int *nz, int *irn, int *jcn, int iw[],
		    int *liw, int *ikeep, int *iw1, int *nsteps,
		    int *iflag, int icntl[30], double cntl[5],
		    int info[20], double *ops );
extern void MA27BD( int *n, int *nz, int *irn, int *jcn,
		    double *a, int *la, int iw[], int *liw,
		    int *ikeep, int *nsteps, int *maxfrt, int *iw1,
		    int icntl[30], double cntl[5], int info[20] );
extern void MA27CD( int *n, double *a, int *la, int iw[], int *liw,
		    double *w, int *maxfrt, double *rhs, int *iw1,
		    int *nsteps, int icntl[30], int info[20] );
//extern void MA27FACTORS( int *n, double a[], int *la, int iw[],
//                         int *liw, int *maxfrt,
//                         int iw2[], int *nblk, int *latop,
//                         int icntl[30], int colrhs[], int *nnzD,
//                         int id[], int jd[], double d[],
//                         int *nnzL, int il[], int jl[], double l[] );
//extern void MA27QDEMASC( int *n, int iw[], int *liwm1, int iw2[],
//                         int *nblk, int *latop, int icntl[30] );


/* Interfaces to the above MA27 subroutines */

Ma27_Data * MA27_Initialize(    int nz, int n, FILE *outfile );
int         MA27_Analyze(       Ma27_Data *data, int iflag );
int         MA27_Factorize(     Ma27_Data *data, double A[] );
int         MA27_Solve(         Ma27_Data *data, double x[] );
int         MA27_Refine(        Ma27_Data *data, double x[], double rhs[],
                                double A[], double tol, int maxitref );
void        MA27_Finalize(      Ma27_Data *data );

#define LIW_MIN    500
#define PIV_MIN   -0.5
#define PIV_MAX    0.5

#endif /* _MA27_H */
