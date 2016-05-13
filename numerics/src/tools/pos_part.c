/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "NSSTools.h"

void pos_part(unsigned n, double* restrict x, double* restrict x_plus)
{

  for (unsigned i = 0 ; i < n ; ++i)
  {
    if (x[i] > 0.0)
    {
      x_plus[i] = x[i];
    }
    else
    {
      x_plus[i] =  0.0;
    }
  }
}
