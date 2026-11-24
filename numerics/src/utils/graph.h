// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*!\file op3x3.h
 * \brief linear algebra operations in 3D*/

#ifndef _graph_h_
#define _graph_h_
#include <stddef.h>
/* node in a list of vertices
 */
struct node {
  int vertex;
  struct node* next;
};

// Create a node
struct node* create_node(int v);
struct node* free_node_list(struct node* list);
/* node in a list of connected component
 */
struct connectedcomponent_node {
  struct node* connectedcomponent;
  struct connectedcomponent_node* next;
};

// Create node in connectedcomponent
struct connectedcomponent_node* create_node_connectedcomponent(struct node* v);

size_t len_connectedcomponent(struct node* connectedcomponent);

void print_connectedcomponent(struct node* connectedcomponent);

struct Graph {
  int numVertices;
  int* visited;
  struct node** adjLists;
};

// Create graph
struct Graph* create_graph(int vertices);
struct Graph* free_graph(struct Graph* graph);

// Add edge
void add_edge(struct Graph* graph, int src, int dest);

// Print the graph
void print_graph(struct Graph* graph);

void DFS(struct Graph* graph, int vertex);

void DFS_compute_connectedcomponent(struct Graph* graph, int vertex,
                                    struct node** connected_component);

int compute_number_connectedcomponents(struct Graph* graph);
struct connectedcomponent_node* compute_connectedcomponents(struct Graph* graph);

unsigned int len_connectedcomponents(struct connectedcomponent_node* connectedcomponentList);

struct connectedcomponent_node* free_connectedcomponents(
    struct connectedcomponent_node* connectedcomponentList);

void print_connectedcomponents(struct connectedcomponent_node* connectedcomponentList);

#endif
