#include "intersection_union.h"
#include<stdio.h>
#include<stdlib.h>

/* Function prints Intersection of arr1[] and arr2[]
   m is the number of elements in arr1[]
   n is the number of elements in arr2[] */
void printIntersection(int arr1[], int arr2[], int m, int n)
{
  int i = 0, j = 0;
  while (i < m && j < n)
  {
    if (arr1[i] < arr2[j])
      i++;
    else if (arr2[j] < arr1[i])
      j++;
    else /* if arr1[i] == arr2[j] */
    {
      printf(" %d ", arr2[j++]);
      i++;
    }
  }
}

void compute_intersection(int * arr1, int * arr2, int m, int n, int * intersection_set, int * intersection_set_size )
{
  
  int i = 0, j = 0;
  int size= 0;
  while (i < m && j < n)
  {
    if (arr1[i] < arr2[j])
      i++;
    else if (arr2[j] < arr1[i])
      j++;
    else /* if arr1[i] == arr2[j] */
    {
      /* printf(" %d ", arr2[j++]); */
      intersection_set[size] = arr2[j++];
      size++;
      i++;
    }
  }
  *intersection_set_size = size;
}

/* Function prints union of arr1[] and arr2[]
   m is the number of elements in arr1[]
   n is the number of elements in arr2[] */
void printUnion(int arr1[], int arr2[], int m, int n)
{
  int i = 0, j = 0;
  while (i < m && j < n)
  {
    if (arr1[i] < arr2[j])
      printf(" %d ", arr1[i++]);
    else if (arr2[j] < arr1[i])
      printf(" %d ", arr2[j++]);
    else
    {
      printf(" %d ", arr2[j++]);
      i++;
    }
  }
 
  /* Print remaining elements of the larger array */
  while(i < m)
   printf(" %d ", arr1[i++]);
  while(j < n)
   printf(" %d ", arr2[j++]);
}

void compute_union(int arr1[], int arr2[], int m, int n, int * union_set, int * union_set_size)
{
  int i = 0, j = 0;
  int size =0;
  while (i < m && j < n)
  {
    if (arr1[i] < arr2[j])
    {
      //printf(" %d ", arr1[i++]);
      union_set[size]=arr1[i++];
      size++;
    }
    else if (arr2[j] < arr1[i])
    {
      //printf(" %d ", arr2[j++]);
      union_set[size]=arr2[j++];
      size++;
    }
    else
    {
      //printf(" %d ", arr2[j++]);
      union_set[size]=arr2[j++];
      size++;
      i++;
    }
  }
 
  /* Print remaining elements of the larger array */
  while(i < m)
    //printf(" %d ", arr1[i++]);
  {
    union_set[size] = arr1[i++];
    size++;
  }
  while(j < n)
  {
    union_set[size] = arr2[j++];
    size++;
  }
  * union_set_size = size;
  //printf(" %d ", arr2[j++]);
  
}
 
/* /\* Driver program to test above function *\/ */
/* int main() */
/* { */
/*   int arr1[] = {1, 2, 4, 5, 6}; */
/*   int arr2[] = {2, 3, 5, 7}; */
/*   int m = sizeof(arr1)/sizeof(arr1[0]); */
/*   int n = sizeof(arr2)/sizeof(arr2[0]); */
/*   printUnion(arr1, arr2, m, n); */


/*   int max_size=m+n ; */
/*   printf("\n"); */
/*   int * union_set = (int*) malloc(max_size*sizeof(int)); */
/*   int union_set_size; */
/*   compute_union(arr1, arr2, m,  n,  union_set, &union_set_size ); */
/*   printf("union_set_size = %i\n", union_set_size); */
/*   for (int i=0; i < union_set_size; i++) printf("% i ", union_set[i]); */

  
/*   getchar(); */
/*   return 0; */
/* } */


/* /\* Driver program to test above function *\/ */
/* int main() */
/* { */
/*   int arr1[] = {1, 2, 4, 5, 6}; */
/*   int arr2[] = {2, 3, 5, 7}; */
/*   int m = sizeof(arr1)/sizeof(arr1[0]); */
/*   int n = sizeof(arr2)/sizeof(arr2[0]); */
/*   printIntersection(arr1, arr2, m, n); */
/*   int max_size=0; */
/*   if (m<n) */
/*   { */
/*     max_size = m; */
/*   } */
/*   else */
/*     max_size = n; */

/*   printf("\n"); */
/*   int * intersectionset = (int*) malloc(max_size*sizeof(int)); */
/*   int intersection_size; */
/*   intersection(arr1, arr2, m,  n,  intersectionset, &intersection_size ); */
/*   printf("intersection_size = %i\n", intersection_size); */
/*   for (int i=0; i < intersection_size; i++) printf("% i ", intersectionset[i]); */
/*   getchar(); */
/*   return 0; */
/* } */
