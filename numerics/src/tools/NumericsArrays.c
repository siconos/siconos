/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "NumericsArrays.h"

void NA_diffns(int *na, int *a, int *nb, int * b, int *nc, int *c)
{

  int pta, ptb, ptc;
  int aa, i;

  pta = 0;
  ptb = 0;
  ptc = 0;

  if (*nb == 0)
  {

    for (i = 0 ; i < *na ; i++)
      c[i] = a[i];
    *nc  = *na;

  }

  else
  {

    for (i = 0 ; i < *na ; i++)
      c[i] = -1;

    while ((pta < *na) && (ptb < *nb))
    {

      aa  = a[pta];

      if (b[ptb] > aa)
      {

        c[ptc] = aa ;
        ptc    = ptc + 1 ;
        pta = pta + 1;
      }
      else if (b[ptb] == aa)
      {

        pta = pta + 1;

      }
      else
      {

        while ((b[ptb] < aa) && (ptb < *nb))
        {


          ptb = ptb + 1;

          if (ptb >= *nb)
          {

            c[ptc] = aa;
            ptc    = ptc + 1;

            break;

          }
        }

      }



    }



    for (i = pta + 1; i < *na ; i++)
    {


      c[ptc] = a[i];
      ptc = ptc + 1;
    }

    *nc = ptc;

  }

}


size_t NA_rm_duplicate(size_t *arr, size_t len)
{
  size_t prev = 0;
  size_t curr = 1;
  size_t last = len - 1;
  while (curr <= last) {
    for (prev = 0; prev < curr && arr[curr] != arr[prev]; ++prev);
    if (prev == curr) {
      ++curr;
    } else {
      arr[curr] = arr[last];
      --last;
    }
  }
  return curr;
}


// Merge arr1[0..n1-1] and arr2[0..n2-1] into
// arr3[0..n1+n2-1]
void NA_merge_sorted_arrays(size_t * arr1, size_t * arr2, size_t n1,
                            size_t n2, size_t *arr3)
{
    size_t i = 0, j = 0, k = 0;

    // Traverse both array
    while (i<n1 && j <n2)
    {
        // Check if current element of first
        // array is smaller than current element
        // of second array. If yes, store first
        // array element and increment first array
        // index. Otherwise do same with second array
        if (arr1[i] < arr2[j])
            arr3[k++] = arr1[i++];
        else
            arr3[k++] = arr2[j++];
    }

    // Store remaining elements of first array
    while (i < n1)
        arr3[k++] = arr1[i++];

    // Store remaining elements of second array
    while (j < n2)
        arr3[k++] = arr2[j++];
}

static void NA_swap(size_t *xp, size_t *yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void NA_sort_bubble(size_t *arr, size_t len)
{
   size_t i, j;
   for (i = 0; i < len-1; i++)
       for (j = 0; j < len-i-1; j++)
           if (arr[j] > arr[j+1])
              NA_swap(&arr[j], &arr[j+1]);
}


size_t  NA_merge_and_sort_sorted_arrays(size_t * arr1, size_t * arr2, size_t n1,
                                     size_t n2, size_t *arr3)
{
  NA_merge_sorted_arrays(arr1, arr2,  n1, n2, arr3);
  int n3 = NA_rm_duplicate(arr3, n1+n2);
  NA_sort_bubble(arr3, n3);
  
  return n3;
}

void NA_display(size_t * arr1,  size_t n1)
{
  printf("Array display:\t[");
    for(size_t j =0 ; j< n1 ; j++)
      printf("%zu\t",arr1[j]);
  printf("]\n");
}


/* swap two indices */
void uint_swap (unsigned int *a, unsigned int *b)
{
  unsigned int temp = *a;
  *a = *b;
  *b = temp;
}

/* shuffle an unsigned array */
void uint_shuffle (unsigned int *a, unsigned int n) {

  for (unsigned int i = 0; i < n - 1; i++)
  {
    uint_swap  (&a[i], &a[i + rand()%(n - i)]);
  }
}

