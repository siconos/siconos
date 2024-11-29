// static char help[] = "Color a matrix, returning set sizes and indices. \n\n";

/*
    Example:
        ./ex16 -f <matrix file> -a_mat_view draw -draw_pause -1
        ./ex16 -f <matrix file> -a_mat_view ascii::ascii_info
*/

#include "graph_tools_petsc.h"
int color_graph_petsc(int n, NumericsMatrix *M, long int *n_colors, size_t **set_sizes, size_t ***set_indices) {
    Mat A;
    // PetscLogDouble time_start, time_end;

    PetscFunctionBeginUser;

    // PetscCall(PetscTime(&time_start));
    switch (M->storageType) {
        case NM_DENSE: {
            // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "DENSE\n"));
            PetscCall(MatCreateSeqDense(PETSC_COMM_SELF, n, n, M->matrix0, &A));
            // MatCreateSeqDense does not work because coloring requires MATSEQAIJ matrix
            PetscCall(MatConvert(A, MATSEQAIJ, MAT_INPLACE_MATRIX, &A));

            /* PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
            PetscCall(MatSetType(A, MATAIJ)); // CSR format
            PetscCall(MatSetSizes(A, n, n, PETSC_DECIDE, PETSC_DECIDE));

            // Fill PETSC sparse matrix
            double *dense = M->matrix0;
            for (int row = 0; row < n; row++) {
                for (int col = 0; col < n; col++) {
                    val = PetscAbsScalar(dense[row + n * col]);
                    if (val > dtol) {
                        PetscCall(MatSetValue(A, row, col, val, INSERT_VALUES));
                    }
                }
            } */
            PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
            PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
            break;
        }
        case NM_SPARSE: {
            // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "SPARSE\n"));

            /* I copied this from Numericsmatrix.c function NM_row_prod_no_diag1x1 */
            CSparseMatrix* sparse;
            if (M->matrix2->origin == NSM_CSR) {
                sparse = NM_csr(M);
            } else {
                sparse = NM_csc_trans(M);
            }

            int64_t* Mp = sparse->p;
            int64_t* Mi = sparse->i;
            double* Mx = sparse->x;

            PetscCall(MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, n, n, Mp, Mi, Mx, &A));
            PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
            PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
            break;
        }
        case NM_SPARSE_BLOCK: {
            fprintf(stderr, "color_graph_petsc :: matrix storage not supported yet %d", M->storageType);
            exit(EXIT_FAILURE);
        }
        default: {
            fprintf(stderr, "color_graph_petsc :: unknown matrix storage %d", M->storageType);
            exit(EXIT_FAILURE);
        }
    }
    // PetscCall(PetscTime(&time_end));
    // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Time to build matrix: %f\n", time_end - time_start));
    

    // Create PETSC sparse matrix
    /* PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    switch (M->storageType) {
        case NM_DENSE: {
            // PetscCall(MatCreateSeqDense(PETSC_COMM_WORLD, n, n, M->matrix0, &A));
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Creating PETSC matrix...\n"));
            PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
            MatSetType(A, MATAIJ);
            MatSetSizes(A, n, n, PETSC_DECIDE, PETSC_DECIDE);
            double* dense = M->matrix0;
            for (int row = 0; row < n; row++) {
                for (int col = 0; col < n; col++) {
                    val = PetscAbsScalar(dense[row + n * col]);
                    if (val > dtol) {
                        PetscCall(MatSetValue(A, row, col, val, INSERT_VALUES));
                    }
                }
            }
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "...PETSC matrix created.\n"));
        }
        case NM_SPARSE: {
            fprintf(stderr, "color_graph_petsc :: storage type not supported yet %d", M->storageType);
            exit(EXIT_FAILURE);
            // PetscCall(MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, n, n, ))

        }
        case NM_SPARSE_BLOCK: {
            fprintf(stderr, "color_graph_petsc :: storage type not supported yet %d", M->storageType);
            exit(EXIT_FAILURE);
        }
        default: {
            fprintf(stderr, "color_graph_petsc :: unknown matrix storage %d", M->storageType);
            exit(EXIT_FAILURE);
    }
    }
    PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
 */

    /* View matrix to check
    PetscCall(MatView(A, PETSC_VIEWER_STDOUT_WORLD)); */

    /*          *
     * COLORING *
     *          */
    // PetscCall(PetscTime(&time_start));
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
    size_t *size = NULL; // Array of sizes of each color set
    size_t **indexes = NULL; // Array of pointers to index sets
    const PetscInt *idxin = NULL;
    PetscCall(ISColoringGetIS(iscoloring, PETSC_USE_POINTER, &nn, &is)); // Get index sets
    // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "n_colors = %ld\n", nn));

    PetscInt size_petsc;

    size = (size_t *)malloc((size_t)nn * sizeof(size_t));
    indexes = (size_t **)malloc((size_t)nn * sizeof(size_t *));

    for (int i = 0; i < (int)nn; i++) {
        PetscCall(ISGetLocalSize(is[i], &size_petsc));
        size[i] = (size_t)size_petsc;
        // PetscCall(ISGetLocalSize(is[i], &size[i])); // This gave me a conversion warning
        indexes[i] = (size_t *)malloc(size[i] * sizeof(size_t)); // allocate indexes
        PetscCall(ISGetIndices(is[i], &idxin)); // Get indices for i-th color

        /*
        COPYING INDICES IN NEW ARRAY
        TO USE IT OUTSIDE THIS FUNCTION
        WITHOUT PETSC

        I COULD ONLY USE POINTERS IF I WROTE BOTH COLORING AND COMPUTING IN PETSC???
        */
        for (int j = 0; j < (int)size[i]; j++) {
            indexes[i][j] = (size_t)idxin[j];
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

    // PetscCall(PetscTime(&time_end));
    // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Time to color: %f\n", time_end - time_start));

    // PetscCall(PetscFinalize());

    return 0;

}
