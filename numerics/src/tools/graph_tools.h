#ifndef GRAPH_TOOLS_H
#define GRAPH_TOOLS_H

#include <petscmat.h>
#include <petscsys.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include "NumericsMatrix.h" 

int color_graph(int n, NumericsMatrix *M, long int *n_colors, size_t **set_sizes, size_t ***set_indices);

#endif
