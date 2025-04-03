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
#include <float.h> // for DBL_EPSILON


int color_graph(int n, NumericsMatrix *M, size_t *n_colors, size_t **set_sizes, size_t ***set_indices);

int color_graph_permut(int n, NumericsMatrix *M, size_t *n_colors, size_t **set_sizes, size_t *inv_permutation);

int color_graph_permut_equitable(int n, NumericsMatrix *M, size_t *n_colors, size_t **sum_sizes, size_t *inv_permutation);

typedef struct node node_t;

node_t *create_node(size_t val);

void push_new_node(node_t **head_node, size_t val);

void push_existing_node(node_t **head_node, node_t *node);

void pop(node_t **head_node, int position);

void print_list(node_t **head);

void free_list(node_t **head);

typedef struct element element_t;

int compare(const void *a, const void *b);

void create_adjacency_lists(int n, NumericsMatrix *M, node_t **adjacency_lists);

#endif
