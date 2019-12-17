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
#include "NumericsVector.h"
#include <math.h>    // for fabs
#include <stdio.h>   // for fprintf, printf, FILE, stderr
#include <stdlib.h>  // for exit, EXIT_FAILURE
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"   // for DEBUG_PRINTF

void NV_display(double * m, int nRow)
{
  int lin;
  printf("vector of size\t%d\t =\n[", nRow);
  if(nRow == 0)
  {
    printf("]\n");
  }
  for(lin = 0; lin < nRow; lin++)
  {
    printf(" %.15e", m[lin]);
    if(lin != nRow - 1)
      printf(", ");
    else
      printf("]\n");
  }

}
void NV_write_in_file_python(double * m,  int nRow, FILE* file)
{
  if(! m)
  {
    fprintf(stderr, "Numerics, NV_write_in_file_python  failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(file, "size = %d; \n", nRow);
  fprintf(file, "data= [");
  for(int i = 0; i < nRow; i++)
  {
    fprintf(file, "%32.24e,\t ", m[i]);
  }
  fprintf(file, "]");
}

bool NV_equal(double * x, double * y, int n, double tol)
{
  for(int i =0; i< n ; i++)
  {
    if(fabs(x[i] - y[i]) >= tol)
    {
      DEBUG_PRINTF("error %i = %e\n",i, fabs(x[i]) - y[i]);
      return false;
    }
  }
  return true;
}
double NV_max(double * x, int n)
{
  double max = x[0];
  for(int i =1; i< n ; i++)
  {
    max  = fmax(max, x[i]);
  }
  return max;
}
