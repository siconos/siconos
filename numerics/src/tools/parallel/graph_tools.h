/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef GRAPH_TOOLS_H
#define GRAPH_TOOLS_H

#include <float.h>  // for DBL_EPSILON
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "NumericsMatrix.h"
#include "SparseBlockMatrix.h"  // for SparseBlockStructuredMatrix because it's not included in NumericsMatrix.h

/*
Various graph coloring functions
*/

/** Color a graph and returns color sets
 *
 *  \param[in] n number of lines of the matrix M
 *  \param[in] M matrix representing the graph to be colored  (M_ij != 0 iff nodes i and j are
 * connected) \param[out] n_colors number of colors obtained coloring the graph \param[out]
 * set_sizes array of size n_colors containing the size of each color set \param[out]
 * set_indices pointers to arrays containing the elements of each color set \return 0 if
 * succeed
 */
int color_graph(int n, NumericsMatrix *M, size_t *n_colors, size_t **set_sizes,
                size_t ***set_indices);

/** Color a graph and return reordered lines
 *
 *  \param[in] n number of lines of the matrix M
 *  \param[in] M matrix representing the graph to be colored  (M_ij != 0 iff nodes i and j are
 * connected) \param[out] n_colors number of colors obtained coloring the graph \param[out]
 * sum_sizes array of size n_colors + 1 containing the aggregated sizes of each color set, with
 * sum_sizes[0] = 0 \param[out] inv_permutation array of size n such that: color set i = {
 * inv_permutation[j] for sum_sizes[i] <= j < sum_sizes[i + 1] } \return 0 if succeed
 */
int color_graph_permut(int n, NumericsMatrix *M, size_t *n_colors, size_t **sum_sizes,
                       size_t *inv_permutation);

/** EQUITABLY color a graph and return reordered lines. Equitable means the size difference of
 * two color sets is at most 1.
 *
 *  \param[in] n number of lines of the matrix M
 *  \param[in] M matrix representing the graph to be colored  (M_ij != 0 iff nodes i and j are
 * connected) \param[out] n_colors number of colors obtained coloring the graph \param[out]
 * sum_sizes array of size n_colors + 1 containing the aggregated sizes of each color set, with
 * sum_sizes[0] = 0 \param[out] inv_permutation array of size n such that: color set i = {
 * inv_permutation[j] for sum_sizes[i] <= j < sum_sizes[i + 1] } \return 0 if succeed
 */
int color_graph_permut_equitable(int n, NumericsMatrix *M, size_t *n_colors,
                                 size_t **sum_sizes, size_t *inv_permutation);

/** Color a graph defined by matrix blocks and returns color sets
 *
 *  \param[in] nc number of contacts. Must be a divisor of the number of lines of the matrix M
 *  \param[in] M matrix representing the graph to be colored (two nodes are connected iff the
 * block of size `nc`x`nc` at i,j is full of zeros) \param[out] n_colors number of colors
 * obtained coloring the graph \param[out] set_sizes array of size n_colors containing the size
 * of each color set \param[out] set_indices pointers to arrays containing the elements of each
 * color set \return 0 if succeed
 */
int color_graph_block(size_t nc, NumericsMatrix *M, size_t *n_colors, size_t **set_sizes,
                      size_t ***set_indices);

#ifdef __cplusplus
extern "C" {
#endif

int color_graph_block_permut(size_t nc, NumericsMatrix *M, size_t *n_colors,
                             size_t **set_sizes, size_t *inv_permutation);

#ifdef __cplusplus
}
#endif

typedef struct node node_t;

node_t *create_node_(size_t val);

void push_new_node(node_t **head_node, size_t val);

void push_existing_node(node_t **head_node, node_t *node);

void pop(node_t **head_node, int position);

void print_list(node_t **head);

void free_list(node_t **head);

typedef struct element element_t;

int compare(const void *a, const void *b);

void create_adjacency_lists(int n, NumericsMatrix *M, node_t **adjacency_lists);

#endif
