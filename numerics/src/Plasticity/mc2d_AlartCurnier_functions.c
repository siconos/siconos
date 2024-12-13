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
#include "mc2d_AlartCurnier_functions.h"  // for mc2d_computeAlartCu...

#include <assert.h>  // for assert
#include <math.h>    // for sqrt

#include "MohrCoulomb2DProblem.h"                   // for MohrCoulomb2D...
#include "NumericsFwd.h"                            // for MohrCoulomb2D...
#include "NumericsMatrix.h"                         // for NumericsMatrix
#include "mc2d_onecone_nonsmooth_Newton_solvers.h"  // for mc2d_computeNonsmoo...
#include "numerics_verbose.h"                       // for numerics_printf
#include "op3x3.h"                                  // for SET3, eig_3x3
#include "siconos_debug.h"                          // for DEBUG_PRINTF

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

extern computeNonsmoothFunction Function;
/* #define VERBOSE_DEBUG */

/* Alart & Curnier version (Radius = mu*max(0,RVN)) */
void mc2d_computeAlartCurnierSTDOld(double R[3], double velocity[3], double eta, double mu,
                                    double rho[3], double F[3], double A[9], double B[9]) {
  double RVN, RVT, RVS;
  double RhoN = rho[0];
  double RhoT = rho[1];
  double Radius, RV, RV1, RV3, GammaTT, GammaTS, GammaST, GammaSS;

  RVN = R[0] - RhoN * velocity[0];
  RVT = R[1] - RhoT * velocity[1];
  RVS = R[2] - RhoT * velocity[2];
  RV = sqrt(RVT * RVT + RVS * RVS);
  // Radius = mu*R[0];

  if (A) {
    A[0 + 3 * 1] = 0.;
    A[0 + 3 * 2] = 0.;
    A[1 + 3 * 0] = 0.;
    A[2 + 3 * 0] = 0.;
  }

  if (B) {
    B[0 + 3 * 1] = 0.;
    B[0 + 3 * 2] = 0.;
    B[1 + 3 * 0] = 0.;
    B[2 + 3 * 0] = 0.;
  }

  if (RVN >= 0.0) {
    DEBUG_PRINT("Normal part in the cone\n");
    Radius = mu * RVN;
    F[0] = RhoN * (velocity[0]);
    if (A && B) {
      A[0 + 3 * 0] = RhoN;
      B[0 + 3 * 0] = 0.0;
    }
  } else {
    DEBUG_PRINT("Normal part out the cone\n");
    Radius = 0.0;
    F[0] = R[0];
    if (A && B) {
      A[0 + 3 * 0] = 0.0;
      B[0 + 3 * 0] = 1.0;
    }
  }

  // Compute the value of the Alart--Curnier Function and its gradient for the tangential part

  DEBUG_PRINTF("Radius=%le\n", Radius);
  DEBUG_PRINTF("RV=%le\n", RV);

  if (RV <= Radius)  // We are in the disk and Radius is positive
  {
    DEBUG_PRINT("We are in the disk\n");
    F[1] = RhoT * (velocity[1]);
    F[2] = RhoT * (velocity[2]);
    if (A && B) {
      A[1 + 3 * 1] = RhoT;
      A[1 + 3 * 2] = 0.0;
      A[2 + 3 * 1] = 0.0;
      A[2 + 3 * 2] = RhoT;
      B[1 + 3 * 0] = 0.0;
      B[1 + 3 * 1] = 0.0;
      B[1 + 3 * 2] = 0.0;
      B[2 + 3 * 0] = 0.0;
      B[2 + 3 * 1] = 0.0;
      B[2 + 3 * 2] = 0.0;
    }
  } else if (RV > Radius)  // We are out the disk and Radius is postive
  {
    if (Radius > 0) {
      DEBUG_PRINT("We are out the disk and Radius is positive\n");
      RV1 = 1.0 / RV;
      F[1] = R[1] - Radius * RVT * RV1;
      F[2] = R[2] - Radius * RVS * RV1;
      if (A && B) {
        RV3 = RV1 * RV1 * RV1;
        GammaTT = RV1 - RVT * RVT * RV3;
        GammaTS = -RVT * RVS * RV3;
        GammaST = GammaTS;
        GammaSS = RV1 - RVS * RVS * RV3;

        A[1 + 3 * 0] = mu * RhoN * RVT * RV1;
        A[2 + 3 * 0] = mu * RhoN * RVS * RV1;

        A[1 + 3 * 1] = GammaTT * RhoT * Radius;

        A[1 + 3 * 2] = GammaTS * RhoT * Radius;
        A[2 + 3 * 1] = GammaST * RhoT * Radius;

        A[2 + 3 * 2] = GammaSS * RhoT * Radius;

        B[1 + 3 * 0] = -mu * RVT * RV1;

        B[1 + 3 * 1] = 1.0 - GammaTT * Radius;
        B[1 + 3 * 2] = -GammaTS * Radius;

        B[2 + 3 * 0] = -mu * RVS * RV1;

        B[2 + 3 * 1] = -GammaST * Radius;
        B[2 + 3 * 2] = 1.0 - GammaSS * Radius;
      }
    } else {
      DEBUG_PRINT("We are out the disk and Radius is zero\n");

      F[1] = R[1];
      F[2] = R[2];
      if (A && B) {
        A[1 + 3 * 1] = 0.0;
        A[1 + 3 * 2] = 0.0;
        A[2 + 3 * 1] = 0.0;
        A[2 + 3 * 2] = 0.0;

        B[1 + 3 * 0] = 0.0;
        B[1 + 3 * 1] = 1.0;
        B[1 + 3 * 2] = 0.0;
        B[2 + 3 * 0] = 0.0;
        B[2 + 3 * 1] = 0.0;
        B[2 + 3 * 2] = 1.0;
      }
    }
  }
  /*   else // We are out the disk and Radius is negative */
  /*     { */
  /* #ifdef VERBOSE_DEBUG */
  /*       printf("We are out the disk and Radius is negative\n"); */
  /* #endif */

  /*       /\*Version original *\/ */
  /*       F[1] = R[1] ; */
  /*       F[2] = R[2] ; */
  /*       if (A && B){ */
  /*    A[1+3*1]=0.0; */
  /*    A[1+3*2]=0.0; */
  /*    A[2+3*1]=0.0; */
  /*    A[2+3*2]=0.0; */

  /*    B[1+3*0]=0.0; */
  /*    B[1+3*1]=1.0; */
  /*    B[1+3*2]=0.0; */
  /*    B[2+3*0]=0.0; */
  /*    B[2+3*1]=0.0; */
  /*    B[2+3*2]=1.0;} */
  /*     } */

  DEBUG_PRINTF("F[0] = %le\t", F[0]);
  DEBUG_PRINTF("F[1] = %le\t", F[1]);
  DEBUG_PRINTF("F[2] = %le\n", F[2]);

  DEBUG_EXPR(if (A) NM_dense_display(A, 3, 3, 3););
  DEBUG_EXPR(if (B) NM_dense_display(B, 3, 3, 3););
  DEBUG_EXPR(if (B) {
    double diago = 0.0;
    for (int l = 0; l < 3; l++) {
      for (int k = 0; k < 3; k++) {
        if (k == l)
          diago = 1.0;
        else
          diago = 0.0;
        printf("I-B[%i+3*%i] = %le\t", l, k, diago - B[l + 3 * k]);
      }
      printf("\n");
    }
  });
}

/* Alart & Curnier version (Radius = mu*max(0,RVN)) */
void mc2d_computeAlartCurnierSTD(double R[3], double velocity[3], double eta, double mu,
                                 double rho[3], double F[3], double A[9], double B[9]) {
  DEBUG_PRINT("mc2d_computeAlartCurnierSTD starts\n");
  DEBUG_EXPR_WE(for (int i = 0; i < 3; i++) printf("R[%i]= %12.8e,\t velocity[%i]= %12.8e,\n",
                                                   i, R[i], i, velocity[i]););

  SET3(R);
  SET3(velocity);
  SET3(rho);
  SET3(F);
  SET3X3(A);
  SET3X3(B);

  double RVN, RVT, RVS;
  double RhoN = *rho0;
  double RhoT = *rho1;
  double Radius, RV, RV1, RV3, GammaTT, GammaTS, GammaST, GammaSS;

  RVN = *R0 - RhoN * *velocity0;
  RVT = *R1 - RhoT * *velocity1;
  RVS = *R2 - RhoT * *velocity2;
  RV = sqrt(RVT * RVT + RVS * RVS);
  // Radius = mu*R[0];

  if (A00) {
    // always true
    *A01 = 0.;
    *A02 = 0.;

    // these should appear under conditions
    *A10 = 0.;
    *A20 = 0.;
  }

  if (B00) {
    // always true
    *B01 = 0.;
    *B02 = 0.;

    // these should appear under conditions
    *B10 = 0.;
    *B20 = 0.;
  }

  if (RVN > 0.0) {
    DEBUG_PRINT("Normal part in the cone\n");

    Radius = mu * RVN;
    *F0 = RhoN * *velocity0;
    if (A00 && B00) {
      *A00 = RhoN;
      *B00 = 0.0;
    }
  } else {
    DEBUG_PRINT("Normal part out the cone\n");
    Radius = 0.0;
    *F0 = *R0;
    if (A00 && B00) {
      *A00 = 0.0;
      *B00 = 1.0;
    }
  }

  // Compute the value of the Alart--Curnier Function and its gradient for the tangential part

  DEBUG_PRINTF("Radius=%le\n", Radius);
  DEBUG_PRINTF("RV=%le\n", RV);

  if (RV <= Radius)  // We are in the disk and Radius is positive
  {
    DEBUG_PRINT("We are in the disk\n");

    *F1 = RhoT * *velocity1;
    *F2 = RhoT * *velocity2;
    if (A00 && B00) {
      *A11 = RhoT;
      *A12 = 0.0;
      *A21 = 0.0;
      *A22 = RhoT;
      *B10 = 0.0;
      *B11 = 0.0;
      *B12 = 0.0;
      *B20 = 0.0;
      *B21 = 0.0;
      *B22 = 0.0;
    }
  } else if (RV > Radius)  // We are out the disk and Radius is positive
  {
    if (Radius > 0) {
      DEBUG_PRINT("We are out the disk and Radius is positive\n");
      RV1 = 1.0 / RV;
      *F1 = *R1 - Radius * RVT * RV1;
      *F2 = *R2 - Radius * RVS * RV1;
      if (A00 && B00) {
        RV3 = RV1 * RV1 * RV1;
        GammaTT = RV1 - RVT * RVT * RV3;
        GammaTS = -RVT * RVS * RV3;
        GammaST = GammaTS;
        GammaSS = RV1 - RVS * RVS * RV3;

        *A10 = mu * RhoN * RVT * RV1;
        *A20 = mu * RhoN * RVS * RV1;

        *A11 = GammaTT * RhoT * Radius;

        *A12 = GammaTS * RhoT * Radius;
        *A21 = GammaST * RhoT * Radius;

        *A22 = GammaSS * RhoT * Radius;

        *B10 = -mu * RVT * RV1;

        *B11 = 1.0 - GammaTT * Radius;
        *B12 = -GammaTS * Radius;

        *B20 = -mu * RVS * RV1;

        *B21 = -GammaST * Radius;
        *B22 = 1.0 - GammaSS * Radius;
      }
    } else {
      DEBUG_PRINT("We are out the disk and Radius is zero\n");

      *F1 = *R1;
      *F2 = *R2;
      if (A00 && B00) {
        *A11 = 0.0;
        *A12 = 0.0;
        *A21 = 0.0;
        *A22 = 0.0;

        *B10 = 0.0;
        *B11 = 1.0;
        *B12 = 0.0;
        *B20 = 0.0;
        *B21 = 0.0;
        *B22 = 1.0;
      }
    }
  }
}

/* Christensen & Pang version (Radius = mu* R[0])*/
void mc2d_computeAlartCurnierJeanMoreau(double R[3], double velocity[3], double mu,
                                        double rho[3], double F[3], double A[9], double B[9]) {
  double RVN, RVT, RVS;
  double RhoN = rho[0];
  double RhoT = rho[1];
  double Radius, RV, RV1, RV3, GammaTT, GammaTS, GammaST, GammaSS;

  RVN = R[0] - RhoN * velocity[0];
  RVT = R[1] - RhoT * velocity[1];
  RVS = R[2] - RhoT * velocity[2];
  RV = sqrt(RVT * RVT + RVS * RVS);
  Radius = mu * R[0];

  // Compute the value of the Alart--Curnier Function and its gradient for the normal part
  DEBUG_PRINTF("[Numerics]  mc2d_computeAlartCurnierJeanMoreau - RVN = %e\n", RVN);
  if (RVN >= 0.0) {
    DEBUG_PRINT("[Numerics]  mc2d_computeAlartCurnierJeanMoreau - Normal part in the cone\n");
    F[0] = RhoN * (velocity[0]);
    if (A && B) {
      A[0 + 3 * 0] = RhoN;
      B[0 + 3 * 0] = 0.0;
    }
  } else {
    DEBUG_PRINT("[Numerics]  mc2d_computeAlartCurnierJeanMoreau - Normal part out the cone\n");
    F[0] = R[0];
    if (A && B) {
      A[0 + 3 * 0] = 0.0;
      B[0 + 3 * 0] = 1.0;
    }
  }

  // Compute the value of the Alart--Curnier Function and its gradient for the tangential part

  DEBUG_PRINTF("[Numerics]  mc2d_computeAlartCurnierJeanMoreau - Radius=%le\n", Radius);
  DEBUG_PRINTF("[Numerics]  mc2d_computeAlartCurnierJeanMoreau - RV=%le\n", RV);
  if (RV < Radius || RV < 1e-20)  // We are in the disk
  {
    DEBUG_PRINT("[Numerics]  mc2d_computeAlartCurnierJeanMoreau - We are in the disk \n");

    F[1] = RhoT * (velocity[1]);
    F[2] = RhoT * (velocity[2]);
    if (A && B) {
      A[1 + 3 * 1] = RhoT;
      A[1 + 3 * 2] = 0.0;
      A[2 + 3 * 1] = 0.0;
      A[2 + 3 * 2] = RhoT;
      B[1 + 3 * 0] = 0.0;
      B[1 + 3 * 1] = 0.0;
      B[1 + 3 * 2] = 0.0;
      B[2 + 3 * 0] = 0.0;
      B[2 + 3 * 1] = 0.0;
      B[2 + 3 * 2] = 0.0;
    }
  } else  // We are out the disk
  {
    DEBUG_PRINT("[Numerics]  mc2d_computeAlartCurnierJeanMoreau - We are out the disk\n");

    /*        RV1 = 1.0/RV; */
    /*        F[1] = R[1] - Radius*RVT*RV1; */
    /*        F[2] = R[2] - Radius*RVS*RV1; */

    /*        RV3 = RV1*RV1*RV1; */
    /*        GammaTT = (RV - RVT*RVT)*RV3; */
    /*        GammaTS =  - RVT*RVS*RV3; */
    /*        GammaST =  GammaTS; */
    /*        GammaSS = (RV - RVS*RVS)*RV3; */

    RV1 = 1.0 / RV;
    F[1] = R[1] - Radius * RVT * RV1;
    F[2] = R[2] - Radius * RVS * RV1;
    if (A && B) {
      RV3 = RV1 * RV1 * RV1;
      GammaTT = RV1 - RVT * RVT * RV3;
      GammaTS = -RVT * RVS * RV3;
      GammaST = GammaTS;
      GammaSS = RV1 - RVS * RVS * RV3;

      // come back to r2146
      //        A[0+3*1]= mu * RhoN * RVT*RV1;
      //        A[0+3*2]= mu * RhoN * RVS*RV1;

      A[1 + 3 * 1] = GammaTT * RhoT * Radius;

      A[1 + 3 * 2] = GammaTS * RhoT * Radius;

      A[2 + 3 * 1] = GammaST * RhoT * Radius;
      A[2 + 3 * 2] = GammaSS * RhoT * Radius;

      B[1 + 3 * 0] = -mu * RVT * RV1;

      B[1 + 3 * 1] = 1.0 - GammaTT * Radius;
      B[1 + 3 * 2] = -GammaTS * Radius;

      B[2 + 3 * 0] = -mu * RVS * RV1;

      B[2 + 3 * 1] = -GammaST * Radius;
      B[2 + 3 * 2] = 1.0 - GammaSS * Radius;
    }
  }

  DEBUG_EXPR(NV_display(F, 3););

  DEBUG_EXPR(if (A) NM_dense_display(A, 3, 3, 3););
  DEBUG_EXPR(if (B) NM_dense_display(B, 3, 3, 3););
}

void mc2d_compute_rho_split_spectral_norm_cond(MohrCoulomb2DProblem* localproblem,
                                               double* rho) {
  double* MLocal = localproblem->M->matrix0;
  assert(MLocal[0 + 0 * 3] > 0);

  DEBUG_EXPR(NM_dense_display(MLocal, 3, 3, 3););
  double sw = MLocal[1 + 1 * 3] + MLocal[2 + 2 * 3];

  double dw = sw * sw - 4.0 * (sw - MLocal[2 + 1 * 3] + MLocal[1 + 2 * 3]);
  DEBUG_PRINTF("dw = %e\n", dw);
  if (dw > 0.0)
    dw = sqrt(dw);
  else
    dw = 0.0;

  rho[0] = 1.0 / MLocal[0 + 0 * 3];
  rho[1] = 2.0 * (sw - dw) / ((sw + dw) * (sw + dw));
  rho[2] = rho[1];

  assert(rho[0] > 0);
  assert(rho[1] > 0);
  assert(rho[2] > 0);

  DEBUG_PRINTF("sw=%le\t  ", sw);
  DEBUG_PRINTF("dw=%le\n ", dw);
  DEBUG_PRINTF("rho[0]=%le\t", rho[0]);
  DEBUG_PRINTF("rho[1]=%le\t", rho[1]);
  DEBUG_PRINTF("rho[2]=%le\n", rho[2]);
}

void mc2d_compute_rho_split_spectral_norm(MohrCoulomb2DProblem* localproblem, double* rho) {
  double* MLocal = localproblem->M->matrix0;
  assert(MLocal[0 + 0 * 3] > 0);

  DEBUG_EXPR(NM_dense_display(MLocal, 3, 3, 3););
  double sw = MLocal[1 + 1 * 3] + MLocal[2 + 2 * 3];

  double dw = sw * sw - 4.0 * (sw - MLocal[2 + 1 * 3] + MLocal[1 + 2 * 3]);
  DEBUG_PRINTF("dw = %e\n", dw);
  if (dw > 0.0)
    dw = sqrt(dw);
  else
    dw = 0.0;

  rho[0] = 1.0 / MLocal[0 + 0 * 3];

  rho[1] = 2.0 / (sw + dw);
  rho[2] = rho[1];

  assert(rho[0] > 0);
  assert(rho[1] > 0);
  assert(rho[2] > 0);

  DEBUG_PRINTF("sw=%le\t  ", sw);
  DEBUG_PRINTF("dw=%le\n ", dw);
  DEBUG_PRINTF("rho[0]=%le\t", rho[0]);
  DEBUG_PRINTF("rho[1]=%le\t", rho[1]);
  DEBUG_PRINTF("rho[2]=%le\n", rho[2]);
}

void mc2d_compute_rho_spectral_norm(MohrCoulomb2DProblem* localproblem, double* rho) {
  double* MLocal = localproblem->M->matrix0;
  double worktmp[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double eig[3] = {0.0, 0.0, 0.0};
  if (eig_3x3(MLocal, worktmp, eig))
    numerics_printf("mc2d_compute_rho_spectral_norm : failed");
  DEBUG_PRINTF("eig[0] = %4.2e, eig[1] = %4.2e, eig[2] = %4.2e", eig[0], eig[1], eig[2]);
  DEBUG_PRINTF("1/eig[0] = %4.2e, 1/eig[1] = %4.2e, 1/eig[2] = %4.2e", 1.0 / eig[0],
               1.0 / eig[1], 1.0 / eig[2]);
  rho[0] = 1.0 / eig[0];
  rho[1] = rho[0];
  rho[2] = rho[0];

  DEBUG_PRINTF("rho[0]=%le\t", rho[0]);
  DEBUG_PRINTF("rho[1]=%le\t", rho[1]);
  DEBUG_PRINTF("rho[2]=%le\n", rho[2]);
}
