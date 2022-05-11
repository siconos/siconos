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
#ifndef ENUM_TOOL_H
#define ENUM_TOOL_H


typedef struct EnumerationStruct
{
  unsigned long long int current;
  unsigned long long int counter;
  unsigned long long int nb_cases;
  double progress ;
} EnumerationStruct;


EnumerationStruct * enum_init(int M);
int enum_next(int * zw, int size, EnumerationStruct * enum_struct);



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
