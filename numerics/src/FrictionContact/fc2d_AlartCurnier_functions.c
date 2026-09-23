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

#include "fc2d_AlartCurnier_functions.h"

#include <stdio.h>

// #define DEBUG_MESSAGES
// #include "siconos_debug.h"
#include <math.h>  // for fabs
#define MAX(a, b) ((a) > (b) ? (a) : (b))
/**
 * \brief Compute the 2D Alart--Curnier function and its Jacobians.
 *
 * The implementation follows the Jean--Moreau formulation where the friction
 * radius is \f$\mu r_N\f$.
 */
void computeAlartCurnier2D(double* reaction, double* velocity, double mu, double* rho,
                           double* residual, double* jacobian_u, double* jacobian_r) {
  const double rN = reaction[0];
  const double rT = reaction[1];
  const double uN = velocity[0];
  const double uT = velocity[1];
  const double rhoN = rho[0];
  const double rhoT = rho[1];

  /* DEBUG_PRINTF("rhoN = %e\n", rhoN);   */
  /* DEBUG_PRINTF("rhoT= %e\n", rhoT); */

  
  const double RVN = rN - rhoN * uN;
  const double RVT = rT - rhoT * uT;
  

  /* Normal part. */
  if (RVN >= 0.0) {
    if (residual) residual[0] = rhoN * uN;
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 0] = rhoN;
      jacobian_u[1 + 2 * 0] = 0.0;
      jacobian_r[0 + 2 * 0] = 0.0;
      jacobian_r[1 + 2 * 0] = 0.0;
    }
  } else {
    if (residual) residual[0] = rN;
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 0] = 0.0;
      jacobian_u[1 + 2 * 0] = 0.0;
      jacobian_r[0 + 2 * 0] = 1.0;
      jacobian_r[1 + 2 * 0] = 0.0;
    }
  }
  
  //const double radius = mu * MAX(0,rN); //MJ variant
  const double radius = mu * MAX(0,RVN); //AC variant
  /* Tangential part: projection on [-radius, radius]. */
  if (fabs(RVT) <= radius) {
    if (residual) residual[1] = rhoT * uT;
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 1] = 0.0;
      jacobian_u[1 + 2 * 1] = rhoT;
      jacobian_r[0 + 2 * 1] = 0.0;
      jacobian_r[1 + 2 * 1] = 0.0;
    }
  } else {
    /* Tangential part: projection on radius */    
    if (RVT > radius) {
      if (residual) residual[1] = rT - radius;
      if (jacobian_u && jacobian_r) {
        jacobian_u[0 + 2 * 1] = 0.0;
        jacobian_u[1 + 2 * 1] = 0.0;
        jacobian_r[0 + 2 * 1] = -mu;
        jacobian_r[1 + 2 * 1] = 1.0;
      }
    }
    /* Tangential part: projection on -radius */    
    if (RVT < -radius) {
      if (residual) residual[1] = rT + radius;
      if (jacobian_u && jacobian_r) {
        jacobian_u[0 + 2 * 1] = 0.0;
        jacobian_u[1 + 2 * 1] = 0.0;
        jacobian_r[0 + 2 * 1] = mu;
        jacobian_r[1 + 2 * 1] = 1.0;
      }
    }
  }
}
