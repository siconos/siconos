// Demo program for the LBL library.
// D. Orban, 2010.

#include <stdio.h>
#include <stdlib.h>
#include "lbl.h"

int main(void) {

    LBL_Data *lbl;      // Main data structure.
    int lblsolver = 1;  // MA27=0, MA57=1.
    int n, nnz;         // Matrix order and number of nonzeros.
    double *val, *rhs;  // Matrix values and right-hand side.
    FILE *logfile, *mat;

    // Local variables.
    int i, error;

    // Open stream for logging.
    logfile = fopen("lbl.log", "w");
    
    // Read matrix and right-hand side from file.
    mat = fopen("mat.txt", "r");
    fscanf(mat, "%d%d", &n, &nnz);

    // ... Initialize LBL data structure.
    lbl = LBL_Initialize(nnz, n, logfile, lblsolver);

    // ... Allocate memory.
    val = (double *)LBL_Calloc(nnz, sizeof(double));
    rhs = (double *)LBL_Calloc(n, sizeof(double));

    // ... Read in values.
    for(i = 0; i < nnz; i++)
        fscanf(mat, "%d%d%lf", &lbl->irn[i], &lbl->jcn[i], &val[i]);
    for(i = 0; i < n; i++)
        fscanf(mat, "%lf", &rhs[i]);
    fclose(mat);

    // Display input.
    printf("Input matrix:\n");
    for(i = 0; i < nnz; i++)
        printf("%d  %d  %lf\n", lbl->irn[i], lbl->jcn[i], val[i]);
    printf("Input right-hand side:\n");
    for(i = 0; i < n; i++)
        printf("%lf ", rhs[i]);
    printf("\n");

    // Optionally, we may set some specific parameters (MA57 only).
    LBL_set_int_parm(lbl, LBL_I_SCALING, 0);        // No scaling.
    LBL_set_int_parm(lbl, LBL_I_PIV_NUMERICAL, 1);  // Do pivoting.
    LBL_set_real_parm(lbl, LBL_D_PIV_THRESH, 0.5);  // Pivot threshold.

    // Analyze.
    error = LBL_Analyze(lbl, 0);  // iflag=0: automatic pivot choice.
    if(error) {
        fprintf(stderr, "Error return from Analyze: %d\n", error);
        goto fail;
    }

    // Factorize.
    error = LBL_Factorize(lbl, val);
    if(error) {
        fprintf(stderr, "Error return from Factorize: %d\n", error);
        goto fail;
    }

    // Solve.
    error = LBL_Solve(lbl, rhs);
    if(error) {
        fprintf(stderr, "Error return from Solve: %d\n", error);
        goto fail;
    }

    // LBL_Solve() overwrites rhs with solution vector.
    printf("Solution:\n");
    for(i = 0; i < n; i++)
        printf("%lf ", rhs[i]);
    printf("\n");

fail:
    // Terminate.
    LBL_Finalize(lbl);

    // Close logging stream.
    fclose(logfile);

    // Free locally allocated memory.
    LBL_Free(val);
    LBL_Free(rhs);

    return error;
}
