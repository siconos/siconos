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
#include "projectionOnRollingCone.h"

#include "op5x5.h"  // for SET3, eig_3x3
                    //
#include <math.h>   // for sqrt
#include <stdio.h>  // for printf

unsigned int projectionOnRollingCone(double* r, double mu, double mur) {
  double normT = sqrt(r[1] * r[1] + r[2] * r[2]);
  double normMT = sqrt(r[3] * r[3] + r[4] * r[4]);
  /* hypot of libm is sure but really slow */
  /* double normT = hypot(r[1], r[2]); */
  /* double normMT = hypot(r[3], r[4]); */

  if (mu * normT + mur * normMT <= -r[0]) {
    r[0] = 0.0;
    r[1] = 0.0;
    r[2] = 0.0;
    r[3] = 0.0;
    r[4] = 0.0;
    return PROJRCONE_DUAL;
  } else if ((normT <= mu * r[0]) && (normMT <= mur * r[0])) {
    return PROJRCONE_INSIDE;
  } else {
    double mu2 = mu * mu;
    double mur2 = mur * mur;

    double trial_rn = (mu * normT + mur * normMT + r[0]) / (mur2 + mu2 + 1.0);
    if ((normT > mu * trial_rn) && (normMT > mur * trial_rn)) {
      r[0] = trial_rn;
      r[1] = mu * r[0] * r[1] / normT;
      r[2] = mu * r[0] * r[2] / normT;
      r[3] = mur * r[0] * r[3] / normMT;
      r[4] = mur * r[0] * r[4] / normMT;
      return PROJRCONE_BOUNDARY_FRICTION_ROLLING;
    }

    trial_rn = (mu * normT + r[0]) / (mu2 + 1.0);
    if ((normT > mu * trial_rn) && (normMT <= mur * trial_rn)) {
      r[0] = trial_rn;
      r[1] = mu * r[0] * r[1] / normT;
      r[2] = mu * r[0] * r[2] / normT;
      // r[3] = r[3] ;
      // r[4] = r[4] ;
      return PROJRCONE_BOUNDARY_FRICTION;
    }
    trial_rn = (mur * normMT + r[0]) / (mur2 + 1.0);
    if ((normT <= mu * trial_rn) && (normMT > mur * trial_rn)) {
      r[0] = trial_rn;
      // r[1] = r[1] ;
      // r[2] = r[2] ;
      r[3] = mur * r[0] * r[3] / normMT;
      r[4] = mur * r[0] * r[4] / normMT;
      return PROJRCONE_BOUNDARY_ROLLING;
    } else
      return 20;
  }
}
unsigned int projectionOn2DRollingCone(double* r, double mu, double mur) {
  double normT = fabs(r[1]);
  double normMT = fabs(r[2]);

  if (mu * normT + mur * normMT <= -r[0]) {
    r[0] = 0.0;
    r[1] = 0.0;
    r[2] = 0.0;
    return PROJRCONE_DUAL;
  } else if ((normT <= mu * r[0]) && (normMT <= mur * r[0])) {
    return PROJRCONE_INSIDE;
  } else {
    double mu2 = mu * mu;
    double mur2 = mur * mur;

    double trial_rn = (mu * normT + mur * normMT + r[0]) / (mur2 + mu2 + 1.0);
    if ((normT > mu * trial_rn) && (normMT > mur * trial_rn)) {
      r[0] = trial_rn;
      r[1] = mu * r[0] * r[1] / normT;
      r[2] = mur * r[0] * r[2] / normMT;
      return PROJRCONE_BOUNDARY_FRICTION_ROLLING;
    }

    trial_rn = (mu * normT + r[0]) / (mu2 + 1.0);
    if ((normT > mu * trial_rn) && (normMT <= mur * trial_rn)) {
      r[0] = trial_rn;
      r[1] = mu * r[0] * r[1] / normT;
      // r[2] = r[2] ;
      return PROJRCONE_BOUNDARY_FRICTION;
    }
    trial_rn = (mur * normMT + r[0]) / (mur2 + 1.0);
    if ((normT <= mu * trial_rn) && (normMT > mur * trial_rn)) {
      r[0] = trial_rn;
      // r[1] = r[1] ;
      r[2] = mur * r[0] * r[2] / normMT;
      return PROJRCONE_BOUNDARY_ROLLING;
    } else
      return 20;
  }
}
unsigned projectionOnDualRollingCone(double* u, double mu, double mur) { return 0; }

void display_status_rolling_cone(unsigned int status) {
  printf("status = %u\n", status);
  if (status == PROJRCONE_INSIDE) {
    printf("PROJRCONE_INSIDE reaction was inside the cone\n");
  } else if (status == PROJRCONE_DUAL) {
    printf("PROJRCONE_DUAL reaction was inside the dual cone\n");
  } else if (status == PROJRCONE_BOUNDARY_FRICTION) {
    printf("PROJRCONE_BOUNDARY_FRICTION reaction is projected on the friction cone\n");
  } else if (status == PROJRCONE_BOUNDARY_ROLLING) {
    printf("PROJRCONE_BOUNDARY_ROLLING reaction is projected on the rolling cone\n");
  } else if (status == PROJRCONE_BOUNDARY_FRICTION_ROLLING) {
    printf("PROJRCONE_BOUNDARY_FRICTION_ROLLING reaction is projected on the both cones\n");
  }
}
unsigned subdifferentialProjectionOnRollingCone(double* H, double* r, double mu, double mur) {
  double normT = sqrt(r[1] * r[1] + r[2] * r[2]);
  double normMT = sqrt(r[3] * r[3] + r[4] * r[4]);
  SET5X5(H);
  zero5x5(H00);

  if (mu * normT + mur * normMT <= -r[0]) {
    // printf("We are in the polar cone\n");
    return PROJRCONE_DUAL;
  } else if ((normT <= mu * r[0]) && (normMT <= mur * r[0])) {
    *H00 = 1.;
    *H11 = 1.;
    *H22 = 1.;
    *H33 = 1.;
    *H44 = 1.;
    // printf("We are in the cone\n");
    return PROJRCONE_INSIDE;
  } else {
    double mu2 = mu * mu;
    double mur2 = mur * mur;

    double trial_rn = (mu * normT + mur * normMT + r[0]) / (mur2 + mu2 + 1.0);
    if ((normT > mu * trial_rn) && (normMT > mur * trial_rn)) {
      double oneoveroneplusmu2 = 1. / (1. + mu2 + mur2);

      double s1 = r[1] / normT;
      double s2 = r[2] / normT;
      *H00 = 1. * oneoveroneplusmu2;

      double muoveroneplusmu2 = mu * oneoveroneplusmu2;
      *H10 = muoveroneplusmu2 * s1;
      *H20 = muoveroneplusmu2 * s2;

      double s3 = r[1] / normMT;
      double s4 = r[2] / normMT;
      double muroveroneplusmu2 = mur * oneoveroneplusmu2;
      *H30 = muroveroneplusmu2 * s3;
      *H40 = muroveroneplusmu2 * s4;

      *H01 = *H10;
      *H02 = *H20;
      *H03 = *H30;
      *H04 = *H40;

      double muoveroneplusmu2overnormT = muoveroneplusmu2 / normT;
      double r0plusmunormT = r[0] + mu * normT;

      *H11 = muoveroneplusmu2overnormT * (r0plusmunormT - r[0] * s1 * s1);
      *H12 = muoveroneplusmu2overnormT * (-r[0] * s1 * s2);
      *H21 = muoveroneplusmu2overnormT * (-r[0] * s1 * s2);
      *H22 = muoveroneplusmu2overnormT * (r0plusmunormT - r[0] * s2 * s2);



      double muroveroneplusmu2overnormMT = muroveroneplusmu2 / normMT;
      double r0plusmunormMT = r[0] + mur * normMT;

      *H33 = muroveroneplusmu2overnormMT * (r0plusmunormMT - r[0] * s3 * s3);
      *H34 = muroveroneplusmu2overnormMT * (-r[0] * s3 * s4);
      *H43 = muroveroneplusmu2overnormMT * (-r[0] * s3 * s4);
      *H44 = muroveroneplusmu2overnormMT * (r0plusmunormMT - r[0] * s4 * s4);



      return PROJRCONE_BOUNDARY_FRICTION_ROLLING;
    }


    trial_rn = (mu * normT + r[0]) / (mu2 + 1.0);
    if ((normT > mu * trial_rn) && (normMT <= mur * trial_rn)) {
      double oneoveroneplusmu2 = 1. / (1. + mu2);

      double s1 = r[1] / normT;
      double s2 = r[2] / normT;
      *H00 = 1. * oneoveroneplusmu2;

      double muoveroneplusmu2 = mu * oneoveroneplusmu2;
      *H10 = muoveroneplusmu2 * s1;
      *H20 = muoveroneplusmu2 * s2;

      *H01 = *H10;
      *H02 = *H20;

      double muoveroneplusmu2overnormT = muoveroneplusmu2 / normT;
      double r0plusmunormT = r[0] + mu * normT;

      *H11 = muoveroneplusmu2overnormT * (r0plusmunormT - r[0] * s1 * s1);
      *H12 = muoveroneplusmu2overnormT * (-r[0] * s1 * s2);
      *H21 = muoveroneplusmu2overnormT * (-r[0] * s1 * s2);
      *H22 = muoveroneplusmu2overnormT * (r0plusmunormT - r[0] * s2 * s2);

      *H33 = 1.;
      *H44 = 1.;

      return PROJRCONE_BOUNDARY_FRICTION;
    }


    trial_rn = (mur * normMT + r[0]) / (mur2 + 1.0);
    if ((normT <= mu * trial_rn) && (normMT > mur * trial_rn)) {
      double oneoveroneplusmur2 = 1. / (1. + mur2);

      *H00 = 1. * oneoveroneplusmur2;

      double s3 = r[3] / normMT;
      double s4 = r[4] / normMT;
      double muroveroneplusmu2 = mur * oneoveroneplusmur2;

      *H30 = muroveroneplusmu2 * s3;
      *H40 = muroveroneplusmu2 * s4;

      *H03 = *H30;
      *H04 = *H40;

      double muroveroneplusmu2overnormMT = muroveroneplusmu2 / normMT;
      double r0plusmunormMT = r[0] + mur * normMT;

      *H33 = muroveroneplusmu2overnormMT * (r0plusmunormMT - r[0] * s3 * s3);
      *H34 = muroveroneplusmu2overnormMT * (-r[0] * s3 * s4);
      *H43 = muroveroneplusmu2overnormMT * (-r[0] * s3 * s4);
      *H44 = muroveroneplusmu2overnormMT * (r0plusmunormMT - r[0] * s4 * s4);

      *H11 = 1.;
      *H22 = 1.;

      return PROJRCONE_BOUNDARY_ROLLING;
    } else
      return 20;
  }
}
