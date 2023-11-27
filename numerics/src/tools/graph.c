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

#include "graph.h"

#include <stdio.h>
#include <stdlib.h>


struct node* create_node(int v) {
  struct node* newNode = malloc(sizeof(struct node));
  newNode->vertex = v;
  newNode->next = NULL;
  return newNode;
}



struct connectedcomponent_node*  create_node_connectedcomponent(struct node**  v) {
  struct connectedcomponent_node* newNode = malloc(sizeof(struct connectedcomponent_node));
  newNode->connectedcomponent = v;
  newNode->next = NULL;
  return newNode;
}



int len_connectedcomponent(struct node **connectedcomponent) {
  int len =0;
  struct node * temp = connectedcomponent[0];
  while(temp != NULL)
    {
      temp=temp->next;
      len++;
    }
  return len;
}

void add_node_connectedcomponent(struct node** connectedcomponent, struct node* node) {
  int len = len_connectedcomponent(connectedcomponent);
  connectedcomponent[len] = node;
  if (len>0){
    connectedcomponent[len-1]->next =node;
  }
}
void print_connectedcomponent(struct node** connectedcomponent) {
  struct node* temp = connectedcomponent[0];
  printf("node connectedcomponent : [");
  while (temp != NULL)
    {
      printf("%i, ", temp->vertex);
      temp = temp->next;
    }
  printf("]\n");
  int len = len_connectedcomponent(connectedcomponent);
  for (int i =0; i < len; i ++ )
    {
      printf("(%i", connectedcomponent[i]->vertex);
      if (connectedcomponent[i]->next != NULL)
	printf("--> %i)\t", connectedcomponent[i]->next->vertex);
      else
	printf("--> NULL)\t");
    }
  printf("\n");
}




int len_connectedcomponentList(struct connectedcomponent_node ** connectedcomponentList) {
  int len =0;
  while(connectedcomponentList[len] != NULL)
    {
      len++;
    }
  return len;
}
void add_connectedcomponent_in_connectedcomponentList(struct connectedcomponent_node** connectedcomponentList, struct  connectedcomponent_node* connectedcomponent) {
  int len = len_connectedcomponentList(connectedcomponentList);
  connectedcomponentList[len] = connectedcomponent;
  if (len>0){
    connectedcomponentList[len-1]->next = connectedcomponent;
  }
}



// DFS algo
void DFS(struct Graph* graph, int vertex) {
  struct node* adjList = graph->adjLists[vertex];
  struct node* temp = adjList;

  graph->visited[vertex] = 1;
  //printf("Visited %d \n", vertex);

  while (temp != NULL) {
    int connectedVertex = temp->vertex;

    if (graph->visited[connectedVertex] == 0) {
      DFS(graph, connectedVertex);
    }
    temp = temp->next;
  }
}



void DFS_connected_connectedcomponent(struct Graph* graph, int vertex, struct node**  connected_component) {
  struct node* adjList = graph->adjLists[vertex];
  struct node* temp = adjList;

  graph->visited[vertex] = 1;
  //printf("Visited %d \n", vertex);

  struct node* new_node = create_node(vertex);
  add_node_connectedcomponent(connected_component,new_node);

  while (temp != NULL) {
    int connectedVertex = temp->vertex;
    if (graph->visited[connectedVertex] == 0)
      {
      DFS_connected_connectedcomponent(graph, connectedVertex, connected_component);
    }
    temp = temp->next;
  }
}

struct node** free_adj_list(struct node** adjLists, int size){
  for (int k =0; k< size; k++)
    {
      if (adjLists[k] !=NULL) free(adjLists[k]);
    }
  free(adjLists);
  return NULL;


}

// Create graph
struct Graph*  create_graph(int vertices) {
  struct Graph* graph = malloc(sizeof(struct Graph));
  graph->numVertices = vertices;

  graph->adjLists = malloc(vertices * sizeof(struct node*));

  graph->visited = malloc(vertices * sizeof(int));

  int i;
  for (i = 0; i < vertices; i++) {
    graph->adjLists[i] = NULL;
    graph->visited[i] = 0;
  }
  return graph;
}

struct Graph*  free_graph(struct Graph* graph) {
  graph->adjLists = free_adj_list(graph->adjLists, graph->numVertices);
  free(graph->visited);
  free(graph);
  return NULL;
}




// Add edge
void addEdge(struct Graph* graph, int src, int dest) {
  // Add edge from src to dest
  struct node* newNode = create_node(dest);
  newNode->next = graph->adjLists[src];
  graph->adjLists[src] = newNode;

  // Add edge from dest to src
  newNode = create_node(src);
  newNode->next = graph->adjLists[dest];
  graph->adjLists[dest] = newNode;
}

// Print the graph
void printGraph(struct Graph* graph) {
  int v;
  for (v = 0; v < graph->numVertices; v++) {
    struct node* temp = graph->adjLists[v];
    printf("\n Adjacency list of vertex %d\n ", v);
    while (temp) {
      printf("%d -> ", temp->vertex);
      temp = temp->next;
    }
    printf("\n");
  }
}

int compute_number_connectedcomponents(struct Graph *graph) {
  int n_connectedcomponent =0;
  int n_vertices = graph->numVertices;
  for (int n =0; n < n_vertices; n++)
    {
      if (graph->visited[n] ==0)
	{
	  DFS(graph, n);
	  n_connectedcomponent ++;
	}
    }

  return n_connectedcomponent;
}

struct connectedcomponent_node**  compute_connectedcomponents(struct Graph *graph) {

  int n_vertices=graph->numVertices;
  int n_connectedcomponent =0;
  struct connectedcomponent_node** connectedcomponentList =malloc(n_vertices * sizeof(struct connectedcomponent_node*));
  for (int i = 0; i < n_vertices; i++) {
    connectedcomponentList[i] = NULL;
  }

  for (int n =0; n < n_vertices; n++)
    {
      if (graph->visited[n] ==0)
	{
	  /* printf("new connectedcomponent\n"); */
	  struct node** connectedcomponent=malloc(n_vertices * sizeof(struct node*));
	  for (int i = 0; i < n_vertices; i++) {
	    connectedcomponent[i] = NULL;
	  }

	  DFS_connected_connectedcomponent(graph, n, connectedcomponent);

	  // print connectedcomponent
	  /* print_connectedcomponent(connectedcomponent); */

	  add_connectedcomponent_in_connectedcomponentList(connectedcomponentList, create_node_connectedcomponent(connectedcomponent));

	  n_connectedcomponent ++;
	  /* printf("number of connectedcomponent = %i\n",  n_connectedcomponent);  */
	}
    }


  return connectedcomponentList;
}
struct connectedcomponent_node**  free_connectedcomponents(struct connectedcomponent_node** connectedcomponentList, struct Graph *graph) {


  int n_vertices=graph->numVertices;
  for (int i = 0; i < n_vertices; i++) {

    if  (connectedcomponentList[i] != NULL)
      {
	struct node** connectedcomponent= connectedcomponentList[i]->connectedcomponent;

	for (int k = 0; k < n_vertices; k++) {
	  if (connectedcomponent[k] != NULL) free(connectedcomponent[k]);
	}
	free(connectedcomponent);

	free(connectedcomponentList[i]);
      }
  }
  return NULL;


}

void print_connectedcomponents(struct connectedcomponent_node** connectedcomponentList){

  struct connectedcomponent_node* temp = connectedcomponentList[0];
  int c =0;
  while (temp !=NULL)
    {
      printf("connected component number %i :\n", c);
      print_connectedcomponent(temp->connectedcomponent);
      temp=temp->next;
      c++;
    }


  printf("number of connectedcomponent = %i\n",  c);
}

unsigned int len_connectedcomponents(struct connectedcomponent_node** connectedcomponentList){
  struct connectedcomponent_node* temp = connectedcomponentList[0];
  unsigned int c =0;
  while (temp !=NULL)
    {
      temp=temp->next;
      c++;
    }
  return c;

}
