/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

/*!\file op3x3.h
 * \brief linear algebra operations in 3D*/

#ifndef _graph_h_
#define _graph_h_




/* node in a list of vertices
 */
struct node {
  int vertex;
  struct node* next;
};

// Create a node
struct node *create_node(int v);
struct node* free_node_list(struct node* list);
/* node in a list of connected component
 */
struct connectedcomponent_node {
  struct node* connectedcomponent;
  struct connectedcomponent_node* next;
};

// Create node in connectedcomponent
struct connectedcomponent_node* create_node_connectedcomponent(struct node* v);

int len_connectedcomponent(struct node* connectedcomponent);

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





void DFS_compute_connectedcomponent(struct Graph* graph, int vertex, struct node**  connected_component);

int compute_number_connectedcomponents(struct Graph* graph);
struct connectedcomponent_node* compute_connectedcomponents(struct Graph* graph);

unsigned int len_connectedcomponents(struct connectedcomponent_node* connectedcomponentList);

struct connectedcomponent_node* free_connectedcomponents(
    struct connectedcomponent_node* connectedcomponentList);


void print_connectedcomponents(struct connectedcomponent_node* connectedcomponentList);

#endif

