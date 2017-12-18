/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "SiconosConfig.h"
#include "NumericsMatrix.h"

#ifdef WITH_LPSOLVE
#include "vertex_extraction.h"

int main(void)
{
  double Hdat[] = {1,  0, -1, 0,
                0,  1, 0, -1};

  NumericsMatrix* H = NM_create_from_data(NM_DENSE, 4, 2, Hdat);

  double K[] = {-2, -3, -7, -8};

  polyhedron P = { SICONOS_SET_POLYHEDRON, 4, 0, H, K, NULL, NULL};

  int basis[11] = {0};

  siconos_find_vertex(&P, 2, basis);

  for(unsigned i = 0; i < 11; ++i) printf("%d ", basis[i]);
  printf("\n");
  return 0;
}
#else

#error no LP solver configured

#endif
