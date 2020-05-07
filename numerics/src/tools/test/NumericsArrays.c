/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

/*
  Tests functions for NumericsMatrix structure

 */

#include <stdio.h>           // for printf, size_t
#include "NumericsArrays.h"  // for NA_display, NA_merge_and_sort_sorted_arrays

static int NumericsArrays_merge_test0()
{
  int info =0;
  int n1 =3;
  size_t arr1[]=  {0, 1, 2};
  int n2 =3;
  size_t arr2[]=  {0, 1, 2};
  int n3 =6;
  size_t  arr3[]=  {99, 99,99, 99, 99, 99};

  NA_merge_sorted_arrays(arr1, arr2, n1, n2, arr3);
  NA_display(arr3, n3);

  int n_rm = NA_rm_duplicate(arr3, n3);
  NA_display(arr3, n_rm);

  NA_sort_bubble(arr3, n_rm);
  NA_display(arr3, n_rm);

  int n = NA_merge_and_sort_sorted_arrays(arr1, arr2, n1, n2, arr3);
  NA_display(arr3, n);

  return info;
}



static int NumericsArrays_merge_test1()
{
  int info =0;
  int n1 =3;
  size_t arr1[]=  {1, 1, 1};
  int n2 =3;
  size_t arr2[]=  {1, 1, 1};
  int n3 =6;
  size_t  arr3[]=  {99, 99,99, 99, 99, 99};

  NA_merge_sorted_arrays(arr1, arr2, n1, n2, arr3);
  NA_display(arr3, n3);

  int n_rm = NA_rm_duplicate(arr3, n3);
  NA_display(arr3, n_rm);

  NA_sort_bubble(arr3, n_rm);
  NA_display(arr3, n_rm);

  int n = NA_merge_and_sort_sorted_arrays(arr1, arr2, n1, n2, arr3);
  NA_display(arr3, n);

  return info;
}


int main(void)
{

  printf("========= Starts Numerics tests for NumericsArrays ========= \n");

  int info = NumericsArrays_merge_test0();
  info = NumericsArrays_merge_test1();



  printf("========= End Numerics tests for NumericsArrays ========= \n");
  return info;
}

