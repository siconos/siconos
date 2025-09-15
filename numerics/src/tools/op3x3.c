/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "op3x3.h"

/** print a matrix
 * \param mat double* a
 */
#include <stdio.h>
void print3x3(double* mat) {
  SET3X3(mat);

  printf("%10.4g ", *mat00);
  printf("%10.4g ", *mat01);
  printf("%10.4g\n", *mat02);

  printf("%10.4g ", *mat10);
  printf("%10.4g ", *mat11);
  printf("%10.4g\n", *mat12);

  printf("%10.4g ", *mat20);
  printf("%10.4g ", *mat21);
  printf("%10.4g\n", *mat22);
}

/** print a vector
 * \param[in] v double*
 */
void print3(double* v) {
  printf("%10.4g\n", *v++);
  printf("%10.4g\n", *v++);
  printf("%10.4g\n", *v);
}
