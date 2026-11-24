/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef ENUM_TOOL_H
#define ENUM_TOOL_H

typedef struct EnumerationStruct {
  unsigned long long int current;
  unsigned long long int counter;
  unsigned long long int nb_cases;
  double progress;
} EnumerationStruct;

EnumerationStruct* enum_init(int M);
int enum_next(int* zw, int size, EnumerationStruct* enum_struct);

/** Compute the total number of cases that should be enumerated
 * \param M the size of the MCLP problem.
 */
unsigned long long int enum_compute_nb_cases(int M);

/* /\** Initialize the enumeration process. */
/*  * \param M the size of the MCLP problem. */
/*  *\/ */
/* void initEnum(int M); */

/* /\** Iterate in the enumeration */
/*  * \param[in,out] the next iterate */
/*  *\/ */
/* int nextEnum(int * W2V, int size); */

#endif
