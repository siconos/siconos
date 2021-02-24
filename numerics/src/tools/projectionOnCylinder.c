/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include "projectionOnCylinder.h"

void projectionOnCylinder(double* r, double  R)
{

  double normTsquare = r[1] * r[1] + r[2] * r[2];

  if(r[0] >= 0)
  {
    if(normTsquare <= R * R)
    {
      return ;
    }
    else
    {
      normTsquare = sqrt(normTsquare);
      r[1] = R * r[1] / normTsquare;
      r[2] = R * r[2] / normTsquare;
      return;
    }
  }
  else
  {
    r[0] = 0.0;
    if(0 < normTsquare)
    {

      normTsquare = sqrt(normTsquare);
      r[1] = R * r[1] / normTsquare;
      r[2] = R * r[2] / normTsquare;
      /*    r[1]=0.0; */
      /*    r[2]=0.0; */
      return;
    }
    else
    {
      r[1] = 0.0;
      r[2] = 0.0;
      return;
    }

  }
}
void projectionOnGeneralCylinder(double* r, double  R, int dim)
{

}
