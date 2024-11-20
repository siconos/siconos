static char help[] = "Color a matrix, returning set sizes and indices. \n\n";

/*
    Example:
        ./ex16 -f <matrix file> -a_mat_view draw -draw_pause -1
        ./ex16 -f <matrix file> -a_mat_view ascii::ascii_info
*/

#include "graph_tools_petsc.h"
int color_graph_petsc(int n, double *M, int *n_colors, int **set_sizes, int ***set_indices) {
    Mat A;
    PetscScalar val;
    PetscReal dtol = 1.e-16;
    PetscLogDouble time_start, time_end;

    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(NULL, NULL, NULL, help));

    // Create PETSC sparse matrix
    PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
    PetscCall(MatSetType(A, MATAIJ)); // CSR format
    PetscCall(MatSetSizes(A, n, n, PETSC_DECIDE, PETSC_DECIDE));

    // Fill PETSC sparse matrix
    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            val = PetscAbsScalar(M[row + n * col]);
            if (val > dtol) {
                PetscCall(MatSetValue(A, row, col, val, INSERT_VALUES));
            }
        }
    }
    PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

    // Possible to fill the sparse matrix faster ?
    // The line below stores zero entries because M is stored as dense.
    // If M is stored as sparse we could use it
    // PetscCall(MatSetValues(A, n, idx, n, idx, val, INSERT_VALUES));

    /* View matrix to check
    PetscCall(MatView(A, PETSC_VIEWER_STDOUT_WORLD)); */

    /*          *
     * COLORING *
     *          */
    PetscCall(PetscTime(&time_start));
    MatColoring mc;
    ISColoring iscoloring;

    PetscCall(MatColoringCreate(A, &mc));
    PetscCall(MatColoringSetDistance(mc, 1));

    PetscCall(MatColoringSetType(mc, MATCOLORINGJP)); // Coloring algorithm 
    // PetscCall(MatColoringSetType(mc, MATCOLORINGGREEDY));

    PetscCall(MatColoringSetFromOptions(mc));
    PetscCall(MatColoringApply(mc, &iscoloring));
    // PetscCall(MatColoringView(mc, PETSC_VIEWER_STDOUT_WORLD)); // View coloring
    // PetscCall(ISColoringView(iscoloring, PETSC_VIEWER_STDOUT_WORLD)); // View IScoloring

    /* Get index sets for each color */
    PetscInt nn;
    IS *is;
    int *size = NULL; // Array of sizes of each color set
    int **indexes = NULL; // Array of pointers to index sets
    const PetscInt *idxin = NULL;
    PetscCall(ISColoringGetIS(iscoloring, PETSC_USE_POINTER, &nn, &is)); // Get index sets
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "n_colors = %d\n", nn));

    size = (int *)malloc(nn * sizeof(int));
    indexes = (int **)malloc(nn * sizeof(int *));

    for (int i = 0; i < nn; i++) {
        PetscCall(ISGetLocalSize(is[i], &size[i])); // conversion PetscInt to int ?
        indexes[i] = (int *)malloc(size[i] * sizeof(int)); // allocate indexes
        PetscCall(ISGetIndices(is[i], &idxin)); // Get indices for i-th color

        /*
        COPYING INDICES IN NEW ARRAY
        TO USE IT OUTSIDE THIS FUNCTION
        WITHOUT PETSC

        I COULD ONLY USE POINTERS IF I WROTE BOTH COLORING AND COMPUTING IN PETSC???
        */
        for (int j = 0; j < size[i]; j++) {
            indexes[i][j] = idxin[j];
        }

        PetscCall(ISRestoreIndices(is[i], &idxin));
    }

    // Call this because of option PETSC_USE_POINTER (see https://petsc.org/release/manualpages/IS/ISColoringGetIS/)
    PetscCall(ISColoringRestoreIS(iscoloring, PETSC_USE_POINTER, &is));

    /* for (int i = 0; i < nn; i++)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%d : [", i));
        for (int j = 0; j < size[i]; j++)
        {
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, " %d ", indexes[i][j]));
        }
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "]\n"));
    } */

    PetscCall(MatColoringDestroy(&mc));
    PetscCall(ISColoringDestroy(&iscoloring));
    PetscCall(MatDestroy(&A));

    /* for (int i = 0; i < nn; i++)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "size[%d] = %d\n", i, size[i]));
    } */

    *n_colors = nn;
    *set_sizes = size;
    *set_indices = indexes;

    PetscCall(PetscTime(&time_end));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Time to color: %f\n", time_end - time_start));

    PetscCall(PetscFinalize());

    return 0;

}
