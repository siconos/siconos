// static char help[] = "Color a matrix, returning set sizes and indices. \n\n";

/*
    Example:
        ./ex16 -f <matrix file> -a_mat_view draw -draw_pause -1
        ./ex16 -f <matrix file> -a_mat_view ascii::ascii_info
*/

#include "graph_tools.h"
int color_graph(int n, NumericsMatrix *M, size_t *n_colors, size_t **set_sizes, size_t ***set_indices) {
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

    /* TEST IF MATRIX IS SYMMETRIC */
    /* PetscReal eps = 1e-2;
    PetscBool is_symmetric;
    PetscCall(MatIsSymmetric(A, eps, &is_symmetric));

    if (!is_symmetric) {
        PetscCall(MatDestroy(&A));
        // printf("NOT SYMMETRIC\n");
        return 1;
    } */
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

    *n_colors = (size_t)nn;
    *set_sizes = size;
    *set_indices = indexes;

    // PetscCall(PetscTime(&time_end));
    // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Time to color: %f\n", time_end - time_start));

    // PetscCall(PetscFinalize());

    return 0;

}


int color_graph_permut(int n, NumericsMatrix *M, size_t *n_colors, size_t **sum_sizes, size_t *inv_permutation) {
    
    /*                            *
     *   BUILD PETSC MATRIX       *
     *                            */
    Mat A;

    PetscFunctionBeginUser;
    switch (M->storageType) {
        case NM_DENSE: {
            PetscCall(MatCreateSeqDense(PETSC_COMM_SELF, n, n, M->matrix0, &A));
            PetscCall(MatConvert(A, MATSEQAIJ, MAT_INPLACE_MATRIX, &A));
            PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
            PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
            break;
        }
        case NM_SPARSE: {
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

    /* TEST IF MATRIX IS SYMMETRIC */
    /* PetscReal eps = 1e-2;
    PetscBool is_symmetric;
    PetscCall(MatIsSymmetric(A, eps, &is_symmetric));

    if (!is_symmetric) {
        PetscCall(MatDestroy(&A));
        // printf("NOT SYMMETRIC\n");
        return 1;
    } */

    /*          *
     * COLORING *
     *          */
    MatColoring mc;
    ISColoring iscoloring;

    PetscCall(MatColoringCreate(A, &mc));
    PetscCall(MatColoringSetDistance(mc, 1));

    PetscCall(MatColoringSetType(mc, MATCOLORINGJP)); // Coloring algorithm 
    // PetscCall(MatColoringSetType(mc, MATCOLORINGGREEDY));

    PetscCall(MatColoringSetFromOptions(mc));
    PetscCall(MatColoringApply(mc, &iscoloring));

    /* Get index sets for each color */
    PetscInt nn;
    IS *is;
    size_t *sum_size = NULL; // Array of sum of sizes of each color set
    const PetscInt *idxin = NULL;
    PetscCall(ISColoringGetIS(iscoloring, PETSC_USE_POINTER, &nn, &is)); // Get index sets

    PetscInt size_petsc;

    sum_size = (size_t *)malloc((size_t)(nn + 1) * sizeof(size_t)); // sum of sizes before color i
    sum_size[0] = 0;

    int k = 0;

    for (int i = 0; i < (int)nn; i++) {
        PetscCall(ISGetLocalSize(is[i], &size_petsc));
        sum_size[i + 1] = (size_t)size_petsc + sum_size[i];
        PetscCall(ISGetIndices(is[i], &idxin)); // Get indices for i-th color

        /*
        COPYING INDICES IN NEW ARRAY
        TO USE IT OUTSIDE THIS FUNCTION
        WITHOUT PETSC

        I COULD ONLY USE POINTERS IF I WROTE BOTH COLORING AND COMPUTING IN PETSC???
        */
        for (int j = 0; j < (int)size_petsc; j++) {
            inv_permutation[k] = (size_t)idxin[j];
            k++;
        }

        PetscCall(ISRestoreIndices(is[i], &idxin));
    }

    // Call this because of option PETSC_USE_POINTER (see https://petsc.org/release/manualpages/IS/ISColoringGetIS/)
    PetscCall(ISColoringRestoreIS(iscoloring, PETSC_USE_POINTER, &is));

    PetscCall(MatColoringDestroy(&mc));
    PetscCall(ISColoringDestroy(&iscoloring));
    PetscCall(MatDestroy(&A));

    *n_colors = (size_t)nn; // conversion from PetscInt (= long int) to size_t (= unsigned long). No issue if nn >= 0, which is the case
    *sum_sizes = sum_size;

    return 0;
}

/* EQUITABLE GRAPH COLORING */


typedef struct node {
    size_t val;
    struct node *next;
} node_t;

node_t *create_node(size_t val) {
    node_t *new_node = (node_t *)malloc(sizeof(node_t));
    new_node->val = val;
    new_node->next = NULL;

    return new_node;
}

void push_new_node(node_t **head_node, size_t val) {
    node_t *new_node = create_node(val);

    new_node->next = *head_node;
    *head_node = new_node;
}

void push_existing_node(node_t **head_node, node_t *node) {
    node->next = *head_node;
    *head_node = node;
}

void pop(node_t **head_node, int position) {
    if (position == 0) {
        *head_node = (*head_node)->next;
    }
    else {
        int i = 1;
        node_t *current_node = (*head_node)->next;
        node_t *previous_node = *head_node;

        while (current_node->next) {
            if (i == position) {
                previous_node->next = current_node->next;
                break;
            }
            else {
                previous_node = current_node;
                current_node = current_node->next;
                i += 1;
            }
        }
    }
}

void free_list(node_t **head) {
    node_t *current_node = *head;
    node_t *next_node;
    do {
      next_node = current_node->next; // to avoid warning of using a pointer after freeing it
      free(current_node);
      current_node = next_node;

    } while (current_node);
}

void print_list(node_t **head)
{
    node_t *current_node = *head;

    printf("[");
    while (current_node->next)
    {
        printf(" %ld ", current_node->val);
        current_node = current_node->next;
    }
    printf("]\n");
}

typedef struct element {
    size_t value;
    size_t index;
} element_t;

int compare(const void *a, const void *b) {

    element_t *elA = (element_t *)a;
    element_t *elB = (element_t *)b;

    return (int)(elB->value - elA->value);
}

void create_adjacency_lists(int n, NumericsMatrix *M, node_t **adjacency_lists) {
    switch (M->storageType) {
        case NM_DENSE: {
            assert(M->matrix0);
            double *dense = M->matrix0;
            
            for (size_t row = 0; row < (size_t)n; row++) {
                for (size_t col = 0; col < (size_t)n; col++) {
                    if (dense[row + col * (size_t)M->size0] > DBL_EPSILON)
                        push_new_node(&adjacency_lists[row], col);
                }
            }
            break;
        }
        case NM_SPARSE: {
            CSparseMatrix* sparse;
            if (M->matrix2->origin == NSM_CSR) {
                sparse = NM_csr(M);
            } else {
                sparse = NM_csc_trans(M);
            }

            int64_t* Mp = sparse->p;
            int64_t* Mi = sparse->i;

            for (int row = 0; row < n; row++) {
                adjacency_lists[row] = create_node(0);
                for (CS_INT p = Mp[row]; p < Mp[row + 1]; ++p) {
                    // push_new_node(&adjacency_lists[row], (size_t)Mi[p]);
                    if (row != Mi[p]) push_new_node(&adjacency_lists[row], (size_t)Mi[p]);
                }
            }
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
}

int color_graph_permut_equitable(int n, NumericsMatrix *M, size_t *n_colors, size_t **sum_sizes, size_t *inv_permutation) {
    
    int err_not_symmetric;

    err_not_symmetric = color_graph_permut(n, M, n_colors, sum_sizes, inv_permutation);

    /* Matrix not symmetric */
    if (err_not_symmetric == 1) return 1;

    // printf("color_graph_permut_equitable:: err_not_symmetric = %d\n", err_not_symmetric);

    size_t tmp_n_colors = *n_colors;
    size_t *tmp_sum_sizes = *sum_sizes;

    int sorted = 1;

    element_t *sizes_and_indexes = (element_t *)malloc(tmp_n_colors * 2 * sizeof(element_t));
    sizes_and_indexes[0].value = tmp_sum_sizes[1] - tmp_sum_sizes[0];
    sizes_and_indexes[0].index = 0;
    for (size_t i = 1; i < tmp_n_colors; i++)
    {
        sizes_and_indexes[i].value = tmp_sum_sizes[i + 1] - tmp_sum_sizes[i];
        sizes_and_indexes[i].index = i;

        if (sizes_and_indexes[i].value > sizes_and_indexes[i - 1].value)
            sorted = 0;
    }
    for (size_t i = tmp_n_colors; i < 2 * tmp_n_colors; i++)
    {
        sizes_and_indexes[i].value = 0;
        sizes_and_indexes[i].index = i;
    }

    /* for (int i = 0; i < 2 * n_colors; i++)
    {
        printf("%d %d\n", sizes_and_indexes[i].value, sizes_and_indexes[i].index);
    } */

    if (sorted == 0)
        qsort(sizes_and_indexes, tmp_n_colors, sizeof(element_t), compare);

    /* printf("Sorted = %d\n", sorted);
    printf("Number of colors = %d\n", tmp_n_colors);
    printf("Max size = %d, min size %d\n", sizes_and_indexes[0].value, sizes_and_indexes[tmp_n_colors - 1].value);
 */
    /* for (int i = 0; i < 2 * tmp_n_colors; i++)
    {
        printf("%d %d\n", sizes_and_indexes[i].value, sizes_and_indexes[i].index);
    } */

    /* ACTUALLY IT IS ALREADY SORTED BY PETSC FROM BIGGEST TO SMALLEST SET
       SO I DON'T REALLY NEED ALL OF THIS
    */

    /* Array containing lists of vertices for each color */
    node_t **vertices_array = (node_t **)malloc(tmp_n_colors * 2 * sizeof(node_t *));
    for (size_t i = 0; i < tmp_n_colors * 2; i++)
    {
        vertices_array[i] = create_node(0);
        if (i < tmp_n_colors)
        {
            for (size_t j = tmp_sum_sizes[i]; j < tmp_sum_sizes[i + 1]; j++)
            {
                push_new_node(&vertices_array[i], inv_permutation[j]);
            }
        }

        /* printf("Color %d. List = ", i);
        print_list(&vertices_array[i]); */
    }

    /* Better adjacency lists */
    node_t **adjacency_lists = (node_t **)malloc((size_t)n * sizeof(node_t *));
    create_adjacency_lists(n, M, adjacency_lists);

    /* for (int row = 0; row < n; row++)
    {
        adjacency_lists[row] = create_node(0);
        for (CS_INT p = Mp[row]; p < Mp[row + 1]; ++p)
        {
            push_new_node(&adjacency_lists[row], Mi[p]);
            // if (row != Mi[p]) push_new_node(&adjacency_lists[row], Mi[p]);
        }

        printf("adjacency_list[%d] = ", row);
        print_list(&adjacency_lists[row]);
    }
    */

    // Color of each vertex
    size_t *colors = (size_t *)malloc((size_t)n * sizeof(size_t));
    size_t current_color = 0;
    for (size_t v = 0; v < (size_t)n; v++)
    {
        if (v == tmp_sum_sizes[current_color + 1])
            current_color += 1;
        colors[inv_permutation[v]] = current_color;
    }

    // printf("color_graph_permut_equitable:: colors done\n");

    /* printf("Colors = [");
    for (int i = 0; i < n; i++)
        printf(" %d ", colors[i]);
    printf("]\n"); */

    // Recoloring
    size_t current_number_of_colors = tmp_n_colors;
    int position;
    int change_is_possible = 1;
    node_t *current_node;
    node_t *current_node_compare;

    size_t i;

    size_t i_max, i_min;
    element_t tmp;

    // int n_iter = 0;

    // printf("Start recoloring...\n");

    while ((sizes_and_indexes[0].value - sizes_and_indexes[current_number_of_colors - 1].value) > 1)
    {

        /* for (int i = 0; i < tmp_n_colors * 2; i++)
        {
        printf("Color %d. List = ", i);
        print_list(&vertices_array[i]);
        }

        printf("\n(value, index) sizes_and_indexes = [");
        for (int i = 0; i < tmp_n_colors * 2; i++)
            printf(" (%d %d) ", sizes_and_indexes[i].value, sizes_and_indexes[i].index);
        printf("]\n"); */


        // n_iter += 1;
        // if (n_iter > 10) break;


        change_is_possible = 1;
        // See if some vertex in biggest color set can be changed
        i_max = sizes_and_indexes[0].index;
        i_min = sizes_and_indexes[current_number_of_colors - 1].index;
        current_node = vertices_array[i_max];
        position = 0;
        while (current_node->next)
        {
            change_is_possible = 1;
            // vertex = current_node->val;

            // Try to change its color to the smallest color set
            /* is it fatser to compare to the nodes in the target color set ?
               or is it faster to compare to the neighbors ?

            It depends on the number of neighbors and the color set sizes,
            but I think it's better to compare to the neihbors

            */

            // Compare to neighbors
            current_node_compare = adjacency_lists[current_node->val];
            while (current_node_compare->next)
            {
                if (colors[current_node_compare->val] == i_min)
                {
                    change_is_possible = 0;
                }
                current_node_compare = current_node_compare->next;
            }

            // Compare to target color set
            /* current_node_compare = vertices_array[i_min];
            while (current_node_compare->next)
            {
                if (NM_get_value(M, vertex, current_node_compare->val) > 0)
                {
                    change_is_possible = 0;
                }
                current_node_compare = current_node_compare->next;
            } */

            if (change_is_possible == 1)
            {
                pop(&vertices_array[i_max], position);                    // remove vertex from biggest color set
                push_existing_node(&vertices_array[i_min], current_node); // add it to smallest color set
                colors[current_node->val] = i_min;

                /* Check sizes */
                sizes_and_indexes[0].value -= 1;
                sizes_and_indexes[current_number_of_colors - 1].value += 1;

                i = 0;
                while (sizes_and_indexes[i].value < sizes_and_indexes[i + 1].value)
                {
                    tmp = sizes_and_indexes[i];
                    sizes_and_indexes[i] = sizes_and_indexes[i + 1];
                    sizes_and_indexes[i + 1] = tmp;
                    i += 1;
                }

                i = current_number_of_colors - 1;
                while (sizes_and_indexes[i].value > sizes_and_indexes[i - 1].value)
                {
                    tmp = sizes_and_indexes[i];
                    sizes_and_indexes[i] = sizes_and_indexes[i - 1];
                    sizes_and_indexes[i - 1] = tmp;
                    i -= 1;
                }

                break;
            }
            else
            {
                current_node = current_node->next;
                position += 1;
            }
        }

        if (change_is_possible == 0)
        {
            // printf("Adding color...\n");
            /* No change was possible from the biggest color set
               So we create a new color
            */
            current_node = vertices_array[i_max];
            current_number_of_colors += 1;
            i_min = current_number_of_colors - 1;
            colors[vertices_array[i_max]->val] = i_min;
            pop(&vertices_array[i_max], 0); // remove first vertex from biggest color set
            push_existing_node(&vertices_array[i_min], current_node);
            // push_head(&vertices_array[i_min], vertex); // add it to smallest color set

            sizes_and_indexes[0].value -= 1;
            sizes_and_indexes[current_number_of_colors - 1].value += 1;

            i = 0;
            while (sizes_and_indexes[i].value < sizes_and_indexes[i + 1].value)
            {
                tmp = sizes_and_indexes[i];
                sizes_and_indexes[i] = sizes_and_indexes[i + 1];
                sizes_and_indexes[i + 1] = tmp;
                i += 1;
            }

            // printf("...color added.\n");
        }
    }

    // printf("Equitable coloring done.\n");

    /* for (int i = 0; i < 2 * tmp_n_colors; i++)
    {
        printf("%d %d\n", sizes_and_indexes[i].value, sizes_and_indexes[i].index);
    }

    for (int i = 0; i < tmp_n_colors * 2; i++)
    {
        printf("Color %d. List = ", i);
        print_list(&vertices_array[i]);
    } */

    /* for (int color = 0; color < current_number_of_colors; color++)
    {
        node_t *current_node = vertices_array[color];

        printf("Color %d. Vertices = [", color);
        while (current_node->next)
        {
            printf(" %d ", current_node->val);
            current_node = current_node->next;
        }
        printf("]\n");
    } */

    /* Modify inv_permutation and sum_sizes */
    size_t new_n_colors = current_number_of_colors;
    size_t *new_sum_sizes = (size_t *)malloc((new_n_colors + 1) * sizeof(size_t));

    new_sum_sizes[0] = 0;
    size_t color;
    i = 0;
    for (size_t j = 0; j < new_n_colors; j++)
    {
        color = sizes_and_indexes[j].index;
        new_sum_sizes[j + 1] = new_sum_sizes[j] + sizes_and_indexes[j].value;
        current_node = vertices_array[color];
        while (current_node->next)
        {
            inv_permutation[i] = current_node->val;
            i += 1;
            current_node = current_node->next;
        }
    }

    free(sizes_and_indexes);

    for (size_t i = 0; i < tmp_n_colors * 2; i++) free_list(&vertices_array[i]);
    free(vertices_array);

    for (size_t i = 0; i < (size_t)n; i++) free_list(&adjacency_lists[i]);
    free(adjacency_lists);
    
    free(colors);

    free(*sum_sizes);

    /* printf("New number of colors = %d\n", new_n_colors); */

    *n_colors = new_n_colors;
    *sum_sizes = new_sum_sizes;

    return 0;
}

/* Check if a block is full of zeros
Inlined for better performance, is it necessary?
*/
static inline bool block_is_full_of_zeros(double *block, size_t n_row, size_t n_col, size_t col_offset) {
    for (int j = 0; j < n_col; j++) {
        for (int i = 0; i < n_row; i++) {
            /* Perhaps I should compare to a threshold here instead of 0 */
            if (block[i + j * col_offset] != 0) return false;
        }
    }
    return true;
}

int color_graph_2(int nc, NumericsMatrix *M, size_t *n_colors, size_t **set_sizes, size_t ***set_indices) {
    assert(M->size0 == M->size1);
    int n = M->size0;
    int d = n / nc; // dimension of contact space
    
    Mat A;

    PetscFunctionBeginUser;

    // Pointers for sparse storage of matrix to be colored.
    int64_t* Mp = NULL;
    int64_t* Mi = NULL;
    double* Mx = NULL;

    switch (M->storageType) {
        case NM_DENSE: {
            /*
            THIS SHOULD BE TESTED !!!!
            */
            double *M_mat = M->matrix0;

            /* Count non zero blocks */
            size_t nnz = 0;
            for (int contact_i = 0; contact_i < nc; contact_i++) {
                for (int contact_j = 0; contact_j < nc; contact_j++) {
                    /* Check if block is 0 */
                    if (block_is_full_of_zeros(&M_mat[d * contact_i + d * contact_j * n], d, d, n) == false) nnz++;
                }
            }

            Mp = (size_t *)malloc((nc + 1) * sizeof(size_t));
            Mi = (size_t *)malloc(nnz * sizeof(size_t));
            Mx = (double *)malloc(nnz * sizeof(double));

            Mp[0] = 0;
            size_t current_index = 0;

            for (int contact_i = 0; contact_i < nc; contact_i++) {
                Mp[contact_i + 1] = Mp[contact_i];
                for (int contact_j = 0; contact_j < nc; contact_j++) {
                    /* Check if block is 0 */
                    if (block_is_full_of_zeros(&M_mat[d * contact_i + d * contact_j * n], d, d, n) == false) {
                        Mp[contact_i + 1]++;
                        Mi[current_index] = contact_j;
                        Mx[current_index] = 1.;
                        current_index++;
                    }
                }
            }

            /* Is it true that MatColoringCreate only supports MATSEQAIJ format? */
            PetscCall(MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, nc, nc, Mp, Mi, Mx, &A));
            PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
            PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
            break;
        }
        case NM_SPARSE: {
            /*
            THIS SHOULD BE TESTED !!!!
            */
            CSparseMatrix* sparse;
            if (M->matrix2->origin == NSM_CSR) {
                sparse = NM_csr(M);
            } else {
                sparse = NM_csc_trans(M);
            }

            CS_INT* Mp_ = sparse->p;
            CS_INT* Mi_ = sparse->i;
            double* Mx_ = sparse->x;

            bool block_is_zero = true;

            // Can I use sparse->nz ???
            size_t nnz;

            CS_INT *ps = (CS_INT *)malloc(d * sizeof(CS_INT));

            // Count number of nonzero blocks
            for (size_t contact_i = 0; contact_i < nc; contact_i++) {

                // Initialize pointers to go trhough lines
                for (size_t i = 0; i < d; i++) ps[i] = Mp_[d * contact_i + i];

                for (size_t current_contact_col = 0; current_contact_col < nc; current_contact_col++) {
                    
                    block_is_zero = true;

                    for (size_t i = 0; i < d; i++) {

                        while ((Mi_[ps[i]] / d < current_contact_col) || (ps[i] < Mp_[d * (contact_i + 1) + i])) {
                            if (Mx_[ps[i]] != 0) block_is_zero = false;
                            ps[i]++;
                        } 
                    }

                    if (block_is_zero == false) nnz +=1;

                }
            }

            // Initialize sparse storage
            Mp = (size_t *)malloc((nc + 1) * sizeof(size_t));
            Mi = (size_t *)malloc(nnz * sizeof(size_t));
            Mx = (double *)malloc(nnz * sizeof(double));

            Mp[0] = 0;
            size_t current_index = 0;

            for (size_t contact_i = 0; contact_i < nc; contact_i++) {

                // Initialize pointers to go trhough lines
                for (size_t i = 0; i < d; i++) ps[i] = Mp_[d * contact_i + i];

                for (size_t current_contact_col = 0; current_contact_col < nc; current_contact_col++) {

                    block_is_zero = true;
                    
                    for (size_t i = 0; i < d; i++) {

                        while ((Mi_[ps[i]] / d < current_contact_col) || (ps[i] < Mp_[d * (contact_i + 1) + i])) {                           
                            if (Mx_[ps[i]] != 0) block_is_zero = false;
                            ps[i]++;
                        } 
                    }

                    if (block_is_zero == false) {
                        Mp[contact_i + 1]++;
                        Mi[current_index] = current_contact_col;
                        Mx[current_index] = 1.;
                        current_index++;
                    }

                }
            }

            PetscCall(MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, nc, nc, Mp, Mi, Mx, &A));
            PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
            PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

            free(ps);
            break;
        }
        case NM_SPARSE_BLOCK: {
            /* Can I make assumptions ? Like all blocks are squared with size d * d ? */

            SparseBlockStructuredMatrix *sbm = M->matrix1;

            /* TO DO

            Can I improve the current method to only go through the block matrix once and not twice ?
            
            If all blocks are squared and the same size, can I use PETSC block matrix?
            MatCreateSeqBAIJ(PETSC_COMM_WORLD, d, n, n, PETSC_DEFAULT, NULL, &A);
            
            For some reason this does not support blocks of size more than 1 ??? 
            MatCreateSeqBAIJWithArrays(PETSC_COMM_WORLD, d, n, n, sbm->blocksize0, sbm->blocksize1, sbm->block, &A);

            Look at MatCreateSeqSBAIJWithArrays, which create a symmetric matrix (what does it mean?) */

            /* Second method: construct sparse directly from block sparse */
            size_t nnz = 0; // Number of non zero blocks
            
            int *nb_row_blocks = (int *)malloc(sbm->blocknumber0 * sizeof(int));
            int *nb_col_blocks = (int *)malloc(sbm->blocknumber1 * sizeof(int));
            int tmp = 0;

            for (int i = 0; i < sbm->blocknumber0; i++) {
                nb_row_blocks[i] = sbm->blocksize0[i] - tmp;
                tmp = sbm->blocksize0[i];
            }

            tmp = 0;
            for (int i = 0; i < sbm->blocknumber1; i++) {
                nb_col_blocks[i] = sbm->blocksize1[i] - tmp;
                tmp = sbm->blocksize1[i];
            }

            int nb_row_block, nb_col_block;
            int i_start_block, j_start_block;

            // Iterate through lines of blocks
            for (int row_block_number = 0; row_block_number < sbm->blocknumber0; row_block_number++) { 
                nb_row_block = nb_row_blocks[row_block_number];
                // Iterate through blocks of a line of blocks
                for (int block_number = sbm->index1_data[row_block_number]; block_number < sbm->index1_data[row_block_number + 1]; block_number++) {
                    nb_col_block = nb_col_blocks[sbm->index2_data[block_number]];
                    // Inline function in if statement, is it good? 
                    if (block_is_full_of_zeros(sbm->block[block_number], nb_row_block, nb_col_block, nb_col_block) == false) nnz++;  
                }
            }

            // Allocations
            Mp = (size_t *)malloc((nc + 1) * sizeof(size_t));
            Mi = (size_t *)malloc(nnz * sizeof(size_t));
            Mx = (double *)malloc(nnz * sizeof(double));

            Mp[0] = 0;
            size_t current_index = 0;

            for (int row_block_number = 0; row_block_number < sbm->blocknumber0; row_block_number++) { 
                nb_row_block = nb_row_blocks[row_block_number];
                i_start_block = sbm->blocksize0[row_block_number] - nb_row_block;
                Mp[row_block_number + 1] = Mp[row_block_number];
                // Iterate through blocks of a line of blocks
                for (int block_number = sbm->index1_data[row_block_number]; block_number < sbm->index1_data[row_block_number + 1]; block_number++) {
                    nb_col_block = nb_col_blocks[sbm->index2_data[block_number]];
                    j_start_block = sbm->blocksize1[sbm->index2_data[block_number]] - nb_col_block;
                    // Iterate within a block
                    if (block_is_full_of_zeros(sbm->block[block_number], nb_row_block, nb_col_block, nb_col_block) == false) {
                        Mp[row_block_number + 1]++;
                        Mi[current_index] = j_start_block / d;
                        Mx[current_index] = 1.;
                        current_index++;
                    }
                }
            }
            /* Is it true that MatColoringCreate only supports MATSEQAIJ format? */
            PetscCall(MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, nc, nc, Mp, Mi, Mx, &A));
            PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
            PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
            break;
        }
        default: {
            fprintf(stderr, "color_graph_2 :: unknown matrix storage %d", M->storageType);
            exit(EXIT_FAILURE);
        }
    }

    // PetscCall(MatView(A, PETSC_VIEWER_STDOUT_WORLD));

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

    // Free matrix sparse storage
    free(Mp);
    free(Mi);
    free(Mx);

    /* for (int i = 0; i < nn; i++)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "size[%d] = %d\n", i, size[i]));
    } */

    *n_colors = (size_t)nn;
    *set_sizes = size;
    *set_indices = indexes;

    // PetscCall(PetscTime(&time_end));
    // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Time to color: %f\n", time_end - time_start));

    // PetscCall(PetscFinalize());

    return 0;

}