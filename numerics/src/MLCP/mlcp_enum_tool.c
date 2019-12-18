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

#include "mlcp_enum_tool.h"
#include <stdio.h>
#include "numerics_verbose.h"

static unsigned long long int sCurrentEnum = 0;
static unsigned long long int sCmpEnum = 0;
static unsigned long long int sNbCase = 0;
static double sProgress = 0;
static int sMm = 0;

static void affectW2V(int * W2V);

void initEnum(int M)
{
  int cmp;
  /*  sCurrentEnum = 0;*/

  sCmpEnum = 0;
  sNbCase = 1;
  sMm = M;

  for(cmp = 0; cmp < sMm; cmp++)
    sNbCase = sNbCase << 1;
  sProgress = 0;
}

void affectW2V(int * W2V)
{
  unsigned long  int aux = sCurrentEnum;
  for(int i = 0; i < sMm; i++)
  {
    W2V[i] = aux & 1;
    aux = aux >> 1;
  }
  if(verbose)
  {
    for(int i = 0; i < sMm; i++)
      printf("wv[%d]=%d \t", i, W2V[i]);
    printf("\n");
  }

}

int nextEnum(int * W2V)
{
  if(sCmpEnum == sNbCase)
    return 0;
  if(sCurrentEnum >= sNbCase)
  {
    sCurrentEnum = 0;
  }
  if(verbose)
    printf("try enum :%d\n", (int)sCurrentEnum);
  affectW2V(W2V);
  sCurrentEnum++;
  sCmpEnum++;
  if(verbose && sCmpEnum > (unsigned long int)sProgress * sNbCase)
  {
    sProgress += 0.001;
    printf(" progress %f %d \n", sProgress, (int) sCurrentEnum);
  }

  return 1;
}
