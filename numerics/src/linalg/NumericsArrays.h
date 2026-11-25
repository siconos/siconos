/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef NUMERICSARRAYS_H
#define NUMERICSARRAYS_H

/*!\file NumericsArrays.h
  Header to collect basic tools for integer arrays
*/

#include <stddef.h>  // for size_t
#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

#if defined(__cplusplus)
extern "C" {
#endif

/**
   Input na, a, nb, b
   Output nc, c
   a and b: interger vectors in increasing order
   c : vector of integers of a that are not in b.
*/
void NA_diffns(int* na, int* a, int* nb, int* b, int* nc, int* c);

/** */
void NA_sortsn_(int* ddl_i, int* sort, int* n);

size_t NA_rm_duplicate(size_t* arr, size_t len);

void NA_sort_bubble(size_t* arr, size_t len);

void NA_merge_sorted_arrays(size_t* arr1, size_t* arr2, size_t n1, size_t n2, size_t* arr3);

size_t NA_merge_and_sort_sorted_arrays(size_t* arr1, size_t* arr2, size_t n1, size_t n2,
                                       size_t* arr3);
void NA_display(size_t* arr1, size_t n1);

size_t NA_max(size_t* arr, size_t len);
size_t NA_min(size_t* arr, size_t len);
/* swap two indices */
void uint_swap(unsigned int* a, unsigned int* b);
/* shuffle an unsigned array */
void uint_shuffle(unsigned int* a, size_t n);

/** @brief swap two indices
 * @param[in,out] a first index
 * @param[in,out] b second index
 */
void size_t_swap(size_t* a, size_t* b);

/** @brief shuffle a size_t array
 * @param[in,out] a the array
 * @param[in] n size of the array
 */
void size_t_shuffle(size_t* a, size_t n);

#if defined(__cplusplus)
}
#endif

#endif
