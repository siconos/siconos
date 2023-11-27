#include <stdio.h>           // for printf, NULL, fclose, fopen, FILE
#include <stdlib.h>          // for free, malloc

#include "graph.h" 

struct Graph* graph_1() {
  int n_vertices =8;
  struct Graph* graph = create_graph(n_vertices);
  addEdge(graph, 0, 1);
  addEdge(graph, 0, 2);
  addEdge(graph, 1, 2);
  addEdge(graph, 2, 3);

  addEdge(graph, 0+4, 1+4);
  addEdge(graph, 0+4, 2+4);
  addEdge(graph, 1+4, 3+4);
  addEdge(graph, 2+4, 3+4);
  printGraph(graph);
  return graph;
}

struct Graph* graph_2() {
  int n_vertices =9;
  struct Graph* graph =  create_graph(n_vertices);
  addEdge(graph, 0, 1);
  addEdge(graph, 0, 2);
  addEdge(graph, 1, 2);
  addEdge(graph, 2, 3);

  addEdge(graph, 0+4, 1+4);
  addEdge(graph, 0+4, 2+4);
  addEdge(graph, 1+4, 3+4);
  addEdge(graph, 2+4, 3+4);

  addEdge(graph, 0, 8);
  //addEdge(graph, 7, 8);

  printGraph(graph);
  return graph;
}






int main(void) {


  struct Graph* graph1 = graph_1();
  int n_connectedcomponent =compute_number_connectedcomponents(graph1);
  printf(" number of connectedcomponent = %i\n",  n_connectedcomponent);
  free_graph(graph1);

  n_connectedcomponent =0;
  struct Graph* graph2 = graph_2();
  struct connectedcomponent_node** connectedcomponentList =  compute_connectedcomponents(graph2);

  struct connectedcomponent_node* temp = connectedcomponentList[0];
  int c =0;
  while (temp !=NULL)
    {
      printf("connected component number %i :\n", c);
      print_node_connectedcomponent(temp->connectedcomponent);
      temp=temp->next;
      c++;
    }
  
  
  printf("number of connectedcomponent = %i\n",  n_connectedcomponent); 

  connectedcomponentList= free_connectedcomponents(connectedcomponentList, graph2);
  
  //free(connectedcomponentList);
  free_graph(graph2);
    
  
  return 0;
}

  
