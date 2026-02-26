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
#include <assert.h>  // for assert
#include <math.h>    // for sqrt

/* #include "MohrCoulomb2DProblem.h"         // for MohrCoulomb2D... */
/* #include "NumericsFwd.h"                  // for MohrCoulomb2D... */
/* #include "NumericsMatrix.h"               // for NumericsMatrix */
/* #include "mc2d_AlartCurnier_functions.h"  // for mc2d_computeAlartCu... */
/* #include "mc2d_naturalmap_functions.h" */
/* #include "mc2d_onecone_nonsmooth_Newton_solvers.h"  // for mc2d_computeNonsmoo... */
/* #include "numerics_verbose.h"                       // for numerics_printf */
#include "op5x5.h"  // for SET3, eig_3x3
#include "projectionOnRollingCone.h"

// extern computeNonsmoothFunction Function;

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"  // for DEBUG_PRINTF
#include "utils/numerics_errors.h"

void rolling_friction_3D_computeNaturalMap(double R[5], double velocity[5], double mu,
                                           double mur, double * Rho, double F[5],
                                           double A[25], double B[25]) {
  DEBUG_PRINT("rolling_friction_3D_computeNaturalMap\n");
  DEBUG_EXPR_WE(for (int i = 0; i < 5; i++) printf("R[%i]= %12.8e,\t velocity[%i]= %12.8e,\n",
                                                   i, R[i], i, velocity[i]););

  SET5(R);
  SET5(velocity);

  double RV[5]; /* = {0. , 0., 0.}; */
  double rho = *Rho;
  double normVT, normOmegaT;

  normVT = sqrt(*velocity1 * *velocity1 + *velocity2 * *velocity2);
  normOmegaT = sqrt(*velocity3 * *velocity3 + *velocity4 * *velocity4);
  
  RV[0] = *R0 - rho * (*velocity0 + mu * normVT + mur * normOmegaT);
  RV[1] = *R1 - rho * *velocity1;
  RV[2] = *R2 - rho * *velocity2;
  RV[3] = *R3 - rho * *velocity3;
  RV[4] = *R4 - rho * *velocity4;

  cpy5(RV, F);

  DEBUG_PRINTF("mu= %12.8e \n", mu);
  DEBUG_PRINTF("mur= %12.8e \n", mur);
  DEBUG_PRINTF("rho= %12.8e \n", rho);
  unsigned int where = projectionOnRollingCone(F, mu, mur);

  DEBUG_EXPR_WE(printf("projection F\n"); display5(F););

  F[0] = *R0 - F[0];
  F[1] = *R1 - F[1];
  F[2] = *R2 - F[2];
  F[3] = *R3 - F[3];
  F[4] = *R4 - F[4];

  DEBUG_EXPR_WE(display5(F););

  DEBUG_EXPR_WE(if (where == PROJRCONE_DUAL) {printf("We are in the polar cone\n");}
		  else if (where == PROJRCONE_INSIDE) {printf("We are in the cone\n");}
		  else 
                    printf("We are outside the cone and its polar\n"););
  if (A && B) {
    
    SET5X5(A);
    SET5X5(B);
    double s1, s2, s3, s4;
    if (where == PROJRCONE_DUAL) {
      zero5x5(A00);
      eye5x5(B00);
    } else if (where == PROJRCONE_INSIDE) {
      if (normVT <= 0.0) {
        s1 = 1.;
        s2 = 0.;
      } else {
        s1 = *velocity1 / normVT;
        s2 = *velocity2 / normVT;
      }
      if (normOmegaT <= 0.0) {
        s3 = 1.;
        s4 = 0.;
      } else {
        s3 = *velocity3 / normOmegaT;
        s4 = *velocity4 / normOmegaT;
      }
      DEBUG_PRINTF("normVT = %6.4e\t, s1 = %6.4e\t, s2 = %6.4e\n ", normVT, s1, s2);
      DEBUG_PRINTF("normOmegaT = %6.4e\t, s3 = %6.4e\t, s4 = %6.4e\n ", normOmegaT, s3, s4);
      zero5x5(A00);

      *A00 = rho;
      *A01 = rho * mu * s1;
      *A02 = rho * mu * s2;
      *A03 = rho * mur * s3;
      *A04 = rho * mur * s4;

      *A11 = rho;
      *A22 = rho;
      *A33 = rho;
      *A44 = rho;

      zero5x5(B00);
    } else {
      if (normVT <= 0.0) {
        s1 = 1.;
        s2 = 0.;
      } else {
        s1 = *velocity1 / normVT;
        s2 = *velocity2 / normVT;
      }
      if (normOmegaT <= 0.0) {
        s3 = 1.;
        s4 = 0.;
      } else {
        s3 = *velocity3 / normOmegaT;
        s4 = *velocity4 / normOmegaT;
      }
      DEBUG_PRINTF("normVT = %6.4e\t, s1 = %6.4e\t, s2 = %6.4e\n ", normVT, s1, s2);
      DEBUG_PRINTF("normOmegaT = %6.4e\t, s3 = %6.4e\t, s4 = %6.4e\n ", normOmegaT, s3, s4);

      double H[25];
      /* //zero5x5(H); */
      // unsigned int where2 =
      subdifferentialProjectionOnRollingCone(H, RV, mu, mur);
      DEBUG_EXPR_WE(printf("H:"); display5x5(H););

      // A = rho * (I+D) * H

      // B = rho * (I+D) we use B for storage
      zero5x5(B00);

      *B00 = rho;
      *B01 = rho * mu * s1;
      *B02 = rho * mu * s2;
      *B03 = rho * mur * s3;
      *B04 = rho * mur * s4;

      *B11 = rho;
      *B22 = rho;
      *B33 = rho;
      *B44 = rho;

      DEBUG_EXPR_WE(printf("rho * (I+D):"); display5x5(B00););
      mm5x5(H, B00, A00);
      // B = I - H
      eye5x5(B00);
      sub5x5(H, B00);

    }
    DEBUG_EXPR_WE(printf("B"); display5x5(B00););
    DEBUG_EXPR_WE(printf("A:"); display5x5(A00););
  }
}
