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
#include <math.h>
#include "projectionOnRollingCone.h"

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "debug.h"

unsigned int projectionOnRollingCone(double* r, double  mu, double mur)
{
  double normT = hypot(r[1], r[2]);
  double normMT = hypot(r[3], r[4]);

  if (mu * normT  + mur * normMT <= - r[0])
  {
    r[0] = 0.0;
    r[1] = 0.0;
    r[2] = 0.0;
    r[3] = 0.0;
    r[4] = 0.0;
    return PROJRCONE_DUAL;
  }
  else if ((normT <= mu * r[0]) && (normMT <= mur * r[0]))
  {
    return PROJRCONE_INSIDE;
  }
  else
  {
    double mu2 = mu * mu;
    double mur2 = mur*mur;


    double trial_rn =  (mu * normT + mur * normMT + r[0]) / (mur2+ mu2 + 1.0);
    if ((normT > mu * trial_rn) && (normMT > mur * trial_rn))
    {
      r[0] = trial_rn;
      r[1] = mu * r[0] * r[1] / normT;
      r[2] = mu * r[0] * r[2] / normT;
      r[3] = mur * r[0] * r[3] / normMT;
      r[4] = mur * r[0] * r[4] / normMT;
      return PROJRCONE_BOUNDARY_FRICTION_ROLLING;
    }

    trial_rn = (mu * normT + r[0]) / (mu2 + 1.0);
    if ((normT > mu * trial_rn) && (normMT <= mur * trial_rn))
    {
      r[0] = trial_rn;
      r[1] = mu * r[0] * r[1] / normT;
      r[2] = mu * r[0] * r[2] / normT;
      //r[3] = r[3] ;
      //r[4] = r[4] ;
      return PROJRCONE_BOUNDARY_FRICTION;
    }
    trial_rn =  (mur * normMT + r[0]) / (mur2 + 1.0);
    if ((normT <= mu * trial_rn) && (normMT > mur * trial_rn))
    {
      r[0] = trial_rn;
      //r[1] = r[1] ;
      //r[2] = r[2] ;
      r[3] = mur * r[0] * r[3] / normMT;
      r[4] = mur * r[0] * r[4] / normMT;
      return PROJRCONE_BOUNDARY_ROLLING;
    }
    else
      return 20;
  }
}
unsigned projectionOnDualRollingCone(double* u, double  mu, double mur)
{
  return 0;

}

void display_status_rolling_cone(unsigned int status)
{
  printf("status = %u\n", status);
  if (status == PROJRCONE_INSIDE)
  {
    printf("PROJRCONE_INSIDE reaction was inside the cone\n");
  }
  else if (status == PROJRCONE_DUAL)
  {
    printf("PROJRCONE_DUAL reaction was inside the dual cone\n");
  }
  else if (status == PROJRCONE_BOUNDARY_FRICTION)
  {
    printf("PROJRCONE_BOUNDARY_FRICTION reaction is projected on the friction cone\n");
  }
  else if (status == PROJRCONE_BOUNDARY_ROLLING)
  {
    printf("PROJRCONE_BOUNDARY_ROLLING reaction is projected on the rolling cone\n");
  }
  else if (status == PROJRCONE_BOUNDARY_FRICTION_ROLLING)
  {
    printf("PROJRCONE_BOUNDARY_FRICTION_ROLLING reaction is projected on the both cones\n");
  }
}
