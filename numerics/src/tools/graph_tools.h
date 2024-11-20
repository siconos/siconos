#ifndef GRAPH_TOOLS_H
#define GRAPH_TOOLS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>

void build_graph(int n, double *M, int **all_neighbors, int *degrees, int *s);
// void compute_degrees(int n, double *M, int *degrees, int *s);
int choose_color(int *palette, int rank);
int compare_color(int n, int *neighbors, int degree, int *colors, int color);
void remove_color(int n, int *neighbors, int degree, int color, int **palettes, int *palettes_sizes);
void color_graph(int n, double *M, int *colors, int *n_colors);
int check_coloring(int n, int **all_neighbors, int *degrees, int *colors);

#endif