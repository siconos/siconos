#ifndef GRAPH_TOOLS_PETSC_H
#define GRAPH_TOOLS_PETSC_H

#include <petscmat.h>
#include <petscsys.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>

int color_graph_petsc(int n, double *M, int *n_colors, int **set_sizes, int ***set_indices);

#endif
