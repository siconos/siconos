#include "graph_tools.h"

void build_graph(int n, double *M, int **all_neighbors, int *degrees, int *s)
{
    *s = n + 1; // Shrinking factor = minimal degree in the graph

    int fill = 0;

    for (int i = 0; i < n; i++)
    {
        fill = 0;
        degrees[i] = 0;
        for (int j = 0; j < n; j++)
        {
            if ((i != j) && ((M[i + n * j] != 0.) || (M[j + n * i] != 0.)))
                degrees[i] += 1;
        };

        all_neighbors[i] = (int *)malloc(degrees[i] * sizeof(int));

        for (int j = 0; j < n; j++)
        {
            if ((i != j) && ((M[i + n * j] != 0.) || (M[j + n * i] != 0.)))
            {
                all_neighbors[i][fill] = j;
                fill += 1;
            }
                
        };

        if ((degrees[i] < *s) && degrees[i] > 0)
            *s = degrees[i];
    }
}

/* void compute_degrees(int n, double *M, int *degrees, int *s)
{
    *s = n + 1; // Shrinking factor = minimal degree in the graph
    for (int i = 0; i < n; i++)
    {
        degrees[i] = 0;
        for (int j = 0; j < n; j++)
        {
            if ((i != j) && ((M[i + n * j] != 0.) || (M[j + n * i] != 0.)))
                degrees[i] += 1;
        };

        if ((degrees[i] < *s) && degrees[i] > 0)
            *s = degrees[i];
    }
}
 */
int choose_color(int *palette, int rank)
{
    int choice = 0;

    while (rank > 0)
    {
        if (palette[choice] == 1)
        {
            rank -= 1;
        }
        choice += 1;
    }

    while (palette[choice] == 0)
    {
        choice += 1;
    }

    return choice;
}

int compare_color(int n, int *neighbors, int degree, int *colors, int color)
{
    for (int j = 0; j < degree; j++)
    {
        if (color == colors[neighbors[j]]) return 0;
    }

    return 1;

    /* for (int j = 0; j < n; j++)
    {
        if ((j != i) && ((M[i + n * j] != 0.) || (M[j + n * i] != 0.)) && (colors[i] == colors[j]))
            return 0;
    }

    return 1; */
}

void remove_color(int n, int *neighbors, int degree, int color, int **palettes, int *palettes_sizes)
{
    for (int j = 0; j < degree; j++)
    {
        if (palettes[neighbors[j]][color] == 1)
        {
            palettes[neighbors[j]][color] = 0;
            palettes_sizes[neighbors[j]] -= 1;
        }
    }

    /* for (int j = 0; j < n; j++)
    {
        if ((j != i) && ((M[i + n * j] != 0.) || (M[j + n * i] != 0.)) && (palettes[j][color] == 1))
        {
            palettes[j][color] = 0;
            palettes_sizes[j] -= 1;
        }
    } */
}

void color_graph(int n, double *M, int *colors, int *n_colors)
{

    srand(time(0));

    double time;

    int no_progress_streak = 0;
    int no_progress = 1;

    int s = 0;
    int *degrees = NULL;
    int **all_neighbors = NULL;
    degrees = (int *)malloc(n * sizeof(int));
    all_neighbors = (int **)malloc(n * sizeof(int *));

    time = omp_get_wtime();
    build_graph(n, M, all_neighbors, degrees, &s);
    time = omp_get_wtime() - time;
    printf("Time to build graph: %fs\n", time);
    printf("Minimal degree of graph: %d\n", s);

    /* if (s > 2) {
        printf("s too big, set it to 2");
        s = 2;
    } */
    int *uncolored = NULL;
    uncolored = (int *)malloc(n * sizeof(int));

    int n_uncolored = n;

    for (int i = 0; i < n; i++)
    {
        uncolored[i] = 1;
    }

    int *palettes_sizes = NULL;
    palettes_sizes = (int *)malloc(n * sizeof(int));

    int *next_colors = NULL;
    next_colors = (int *)malloc(n * sizeof(int));

    int **palettes = NULL;
    palettes = (int **)malloc(n * sizeof(int *));

    for (int i = 0; i < n; i++)
    {
        palettes[i] = (int *)malloc((degrees[i] + 1) * sizeof(int));
        palettes_sizes[i] = 1 + degrees[i] / s;
        // printf("Vertex [%d]: degree = %d, palette_size = %d\n", i, degrees[i], palettes_sizes[i]);
        next_colors[i] = degrees[i] / s + 1;

        for (int j = 0; j < degrees[i] + 1; j++)
        {
            if (j <= degrees[i] / s)
                palettes[i][j] = 1;
            else
                palettes[i][j] = 0;
        }
    }

    // int num_threads = omp_get_max_threads();
    // int size = n / num_threads;

    // #pragma omp parallel default(none) \
    shared(size, n, all_neighbors, degrees, n_uncolored, uncolored, no_progress, next_colors, colors, palettes, palettes_sizes, no_progress_streak)
    // {
        //int thread_id = omp_get_thread_num();

    time = omp_get_wtime();

    int total_no_progress = 0;
    int iter = 0;

    while (n_uncolored > 0)
    {
        // #pragma omp single
        no_progress = 1;

        /* Assign random color to each uncolored vertex */
        // for (int i = thread_id * size; i < (thread_id + 1) * size; i++)
        // #pragma omp single
        for (int i = 0; i < n; i++)
        {
            if (uncolored[i] == 1)
            {
                colors[i] = choose_color(palettes[i], rand() % palettes_sizes[i]);
            }
        }
        // #pragma omp barrier
        
        /* Conflict resolution */
        //for (int i = thread_id * size; i < (thread_id + 1) * size; i++)
        // #pragma omp single
        for (int i = 0; i < n; i++)
        {
            if (uncolored[i] == 1)
            {
                if (compare_color(n, all_neighbors[i], degrees[i], colors, colors[i]) == 1)
                {
                    uncolored[i] = 0;
                    n_uncolored -= 1;
                    remove_color(n, all_neighbors[i], degrees[i], colors[i], palettes, palettes_sizes);
                    no_progress = 0; // progress has been made
                    next_colors[i] = colors[i] + 1; // only useful for n_colors computation
                }
            }
        }
        // #pragma omp barrier

        /* Feed the hungry */
        // for (int i = thread_id * size; i < (thread_id + 1) * size; i++)
        // #pragma omp single
        for (int i = 0; i < n; i++)
        {
            if (uncolored[i] == 1)
            {
                if (palettes_sizes[i] == 0)
                {
                    palettes[i][next_colors[i]] = 1;
                    palettes_sizes[i] += 1;
                    next_colors[i] += 1;
                }
            }
        }
        // #pragma omp barrier

        // #pragma omp single
        if (no_progress == 1)
        {
            no_progress_streak += 1;
            total_no_progress += 1;
            if (no_progress_streak > 5)
            {
                int random_vertex = choose_color(uncolored, rand() % n_uncolored); // select random vertex among uncolored
                palettes[random_vertex][next_colors[random_vertex]] = 1;
                palettes_sizes[random_vertex] += 1;
                next_colors[random_vertex] += 1;
                no_progress_streak = 0;
            }
        }
        // #pragma omp barrier
        iter += 1;
    }
    // }
    time = omp_get_wtime() - time;
    printf("Time in while loop: %fs\n", time);
    printf("Total no progress = %d\n", total_no_progress);
    printf("Total iterations = %d\n", iter);

    *n_colors = 0; // This only works if the starting palette is smaller than the ending one
    for (int i = 0; i < n; i++)
    {
        if (next_colors[i] > *n_colors)
        {
            *n_colors = next_colors[i];
        }
    }

    printf("Number of colors: %d\n", *n_colors);

    /* int check = check_coloring(n, all_neighbors, degrees, colors);
    if (check == 0)
    {
        printf("Bad coloring");
    }
    else if (check == 1)
    {
        printf("Good coloring");
    } */

    /* Free memory */
    free(degrees);
    for (int i = 0; i < n; i++)
    {
        free(all_neighbors[i]);
        free(palettes[i]);
    }
    free(all_neighbors);
    free(palettes);
    free(uncolored);
    free(palettes_sizes);
    free(next_colors);    
}

int check_coloring(int n, int **all_neighbors, int *degrees, int *colors)
{
    for (int i = 0; i < n; i++)
    {
        if (compare_color(n, all_neighbors[i], degrees[i], colors, colors[i]) == 0)
        {
            return 0;
        }
    }

    return 1;
}