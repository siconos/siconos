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

#ifndef NSSTOOLS_H
#define NSSTOOLS_H

/*!\file NSSTools.h
  Header to collect basic tools for integer arrays
*/

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /**
     Input na, a, nb, b
     Output nc, c
     a and b: interger vectors in increasing order
     c : vector of integers of a that are not in b.
  */
  void NA_diffns(int *na, int *a, int *nb, int * b, int *nc, int *c);
  
  /** */
  void NA_sortsn_(int *ddl_i, int *sort, int *n);

  size_t NA_rm_duplicate(size_t *arr, size_t len);

  void NA_sort_bubble(size_t *arr, size_t len);

  
  void NA_merge_sorted_arrays(size_t * arr1, size_t * arr2, size_t n1,
                size_t n2, size_t *arr3);
  
  size_t  NA_merge_and_sort_sorted_arrays(size_t * arr1, size_t * arr2, size_t n1,
                                       size_t n2, size_t *arr3);
  void NA_display(size_t * arr1,  size_t n1);
  
/* swap two indices */
  void uint_swap (unsigned int *a, unsigned int *b);
  /* shuffle an unsigned array */
  void uint_shuffle (unsigned int *a, unsigned int n);
  
  
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

