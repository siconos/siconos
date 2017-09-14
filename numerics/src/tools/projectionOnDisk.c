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
#include <math.h>
#include "projectionOnDisk.h"

void projectionOnDisk(double* r, double  R)
{

  double normTsquare = r[0] * r[0] + r[1] * r[1];
  if (normTsquare <= R * R)
  {
    return ;
  }
  else
  {
    normTsquare = sqrt(normTsquare);
    r[0] = R * r[0] / normTsquare;
    r[1] = R * r[1] / normTsquare;
    return;
  }
}

void projectionOnGeneralDisk(double* r, double  R, int dim)
{

}
