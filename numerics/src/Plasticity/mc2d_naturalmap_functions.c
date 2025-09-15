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
#include "mc2d_naturalmap_functions.h"

#include <assert.h>  // for assert
#include <math.h>    // for sqrt

#include "MohrCoulomb2DProblem.h"                   // for MohrCoulomb2D...
#include "NumericsFwd.h"                            // for MohrCoulomb2D...
#include "NumericsMatrix.h"                         // for NumericsMatrix
#include "mc2d_AlartCurnier_functions.h"            // for mc2d_computeAlartCu...
#include "mc2d_onecone_nonsmooth_Newton_solvers.h"  // for mc2d_computeNonsmoo...
#include "numerics_verbose.h"                       // for numerics_printf
#include "op3x3.h"                                  // for SET3, eig_3x3
#include "projectionOnCone.h"

extern computeNonsmoothFunction Function;

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"  // for DEBUG_PRINTF

void mc2d_computeNaturalMap(double R[3], double velocity[3], double eta, double theta,
                            double Rho[3], double F[3], double A[9], double B[9]) {
  DEBUG_PRINT("mc2d_computeNaturalMap starts\n");
  DEBUG_EXPR_WE(for (int i = 0; i < 3; i++) printf("R[%i]= %12.8e,\t velocity[%i]= %12.8e,\n",
                                                   i, R[i], i, velocity[i]););

  SET3(R);
  SET3(velocity);
  SET3(Rho);

  SET3X3(A);
  SET3X3(B);

  double RV[3]; /* = {0. , 0., 0.}; */
  double rho = *Rho0;
  double normVT;

  normVT = sqrt(*velocity1 * *velocity1 + *velocity2 * *velocity2);
  RV[0] = *R0 - rho * (*velocity0 + theta * normVT);
  RV[1] = *R1 - rho * *velocity1;
  RV[2] = *R2 - rho * *velocity2;

  cpy3(RV, F);

  DEBUG_PRINTF("theta= %12.8e \n", theta);
  DEBUG_PRINTF("eta= %12.8e \n", eta);
  DEBUG_PRINTF("rho= %12.8e \n", rho);
  unsigned int where = projectionOnCone(F, eta);

  DEBUG_EXPR_WE(for (int i = 0; i < 3; i++) printf("projection F[%i]= %12.8e \n", i, F[i]));

  F[0] = *R0 - F[0];
  F[1] = *R1 - F[1];
  F[2] = *R2 - F[2];

  DEBUG_EXPR_WE(for (int i = 0; i < 3; i++) printf("F[%i]= %12.8e \n", i, F[i]));

  DEBUG_EXPR_WE(if (where == PROJCONE_DUAL) printf("We are in the polar cone\n");
                if (where == PROJCONE_INSIDE) printf("We are in the cone\n");
                if (where == PROJCONE_BOUNDARY)
                    printf("We are outside the cone and its polar\n"););
  if (A00 && B00) {
    double s1, s2;
    if (where == PROJCONE_DUAL) {
      zero3x3(A00);
      eye3x3(B00);
    } else if (where == PROJCONE_INSIDE) {
      if (normVT <= 0.0) {
        s1 = 1.;
        s2 = 0.;
      } else {
        s1 = *velocity1 / normVT;
        s2 = *velocity2 / normVT;
      }
      DEBUG_PRINTF("normVT = %6.4e\t, s1 = %6.4e\t, s2 = %6.4e\n ", normVT, s1, s2);

      *A00 = rho;
      *A01 = rho * theta * s1;
      *A02 = rho * theta * s2;
      *A10 = 0.;
      *A11 = rho;
      ;
      *A12 = 0.;
      *A20 = 0.;
      *A21 = 0.;
      *A22 = rho;
      zero3x3(B00);

    } else  // where ==  PROJCONE_BOUNDARY
    {
      if (normVT <= 0.0) {
        s1 = 1.;
        s2 = 0.;
      } else {
        s1 = *velocity1 / normVT;
        s2 = *velocity2 / normVT;
      }
      DEBUG_PRINTF("normVT = %6.4e\t, s1 = %6.4e\t, s2 = %6.4e\n ", normVT, s1, s2);
      double H[9];
      /* //zero3x3(H); */
      unsigned int where2 = subdifferentialProjectionOnCone(H, RV, eta);
      DEBUG_EXPR_WE(printf("H:"); display3x3(H););

      // A = rho * (I+D) * H

      // B = rho * (I+D) we use B for storage
      *B00 = rho;
      *B01 = rho * theta * s1;
      *B02 = rho * theta * s2;
      *B10 = 0.;
      *B11 = rho;
      ;
      *B12 = 0.;
      *B20 = 0.;
      *B21 = 0.;
      *B22 = rho;
      DEBUG_EXPR_WE(printf("rho * (I+D):"); display3x3(B00););
      mm3x3(H, B00, A00);
      DEBUG_EXPR_WE(printf("A:"); display3x3(A00););

      // B = I - H
      eye3x3(B00);
      sub3x3(H, B00);
      DEBUG_EXPR_WE(printf("B:"); display3x3(B00););
    }
  }

  // getchar();
}
