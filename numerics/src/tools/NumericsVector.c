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
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "sanitizer.h"
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"

#include "NumericsVector.h"
void NV_display(double * m, int nRow)
{
  int lin;
  printf("vector of size\t%d\t =\n[", nRow);
  if (nRow == 0)
  {
    printf("]\n");
  }
  for (lin = 0; lin < nRow; lin++)
  {
    printf(" %.15e", m[lin]);
    if (lin != nRow - 1)
      printf(", ");
    else
      printf("]\n");
  }

}
void NV_write_in_file_python(double * m,  int nRow, FILE* file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, NV_write_in_file_python  failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(file, "size = %d; \n", nRow);
  fprintf(file, "data= [");
  for (int i = 0; i < nRow; i++)
  {
    fprintf(file, "%32.24e,\t ", m[i]);
  }
  fprintf(file, "]");
}
