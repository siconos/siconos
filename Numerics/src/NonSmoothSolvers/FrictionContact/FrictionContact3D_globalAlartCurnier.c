/* Siconos-Numerics, Copyright INRIA 2005-2010.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/

#include "LA.h"
#include "op3x3.h"
#include "FrictionContact3D_Solvers.h"
#include "FrictionContactProblem.h"
#include "FrictionContact3D_compute_error.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "Friction_cst.h"


//#define VERBOSE_DEBUG

void printm(unsigned int nl, unsigned int nc, double *m)
{
  for (unsigned int i = 0; i < nl; ++i)
  {
    for (unsigned int j = 0; j < nc; ++j)
    {
      printf("%10.4g ", *(m + j + nc * i));
    }
    printf("\n");
  }
}

double sign(double x)
{
  return copysign(1., x);
}


void frictionContact3D_AlartCurnierFABGenerated(
  double rn,
  double rt1,
  double rt2,
  double un,
  double ut1,
  double ut2,
  double mu,
  double rhon,
  double rhot1,
  double rhot2,
  double *result)
{
  double x0 = rhon * un;
  double x1 = rn - x0;
  double x2 = sign(x1);
  double x3 = rhot1 * ut1;
  double x4 = rt1 - x3;
  double x5 = rn / 2;
  double x6 = fabs(x1);
  double x7 = x6 / 2;
  double x8 = x0 / 2;
  double x9 = x5 + x7;
  double x10 = x9 - x8;
  double x11 = x4 * x4;
  double x12 = rhot2 * ut2;
  double x13 = rt2 - x12;
  double x14 = x13 * x13;
  double x15 = x11 + x14;
  double x16 = mu * x10;
  double x17 = sqrt(x15);
  double x19 = 1. / (x17);
  double x22 = x19 * x19 * x19;
  double x23 = x2 / 2;
  double x24 = -rhon / 2;
  double x25 = x2 * x24;
  double x26 = x24 + x25;
  double x27 = 1.0 / 2.0 + x23;
  double x28 = x3 - rt1;
  double x30 = x12 - rt2;
  double x31 = x16 * x19;
  double x32 = mu * x10 * x22;
  double x33 = mu * x4;
  double x34 = x30 * x32;
  double x35 = x32 * x4;
  double x36 = x22 * x4;
  double x37 = mu * x19;
  double x38 = x28 * x32;
  double x39 = x16 * x22;
  double x40 = x13 * x16 * x36;
  double x41 = x19 * x33;
  double x42 = rhot1 * x39;
  double x43 = x27 * x37;
  double x44 = rhot2 * x39;
  double x45 = x26 * x37;
  double x46 = x13 * x34;
  double x47 = rhot1 * x40;
  double x48 = rhot2 * x40;
  double x49 = x34 * x4;
  double x50 = x13 * x38;
  double x51 = x28 * x35;
  double x52 = x27 * x41;
  double x53 = x26 * x41;
  double x54 = x14 * x44;
  double x55 = x13 * x45;
  double x56 = x13 * x43;
  double x57 = x11 * x42;
  double x58 = -rhot2;
  double x59 = rhon / 2;
  double x60;
  double x61;
  double x62;
  double x63;
  double x64;
  double x65;
  double x66 = rhot1 * x31;
  double x67 = -rhot1;
  double x68;
  double x69;
  double x70 = x13 * x31;
  double x71 = rhon * x23;
  double x72;
  double x73 = x31 * x4;
  double x74;
  double x75 = x31 + x46;
  double x76 = x31 + x51;
  double x77 = x5 + x8;
  double x78 = x59 + x71;
  if (x17 <= x16)
  {
    x60 = 0;
    x61 = 0;
    x62 = 0;
    x63 = -x58;
    x64 = 0;
    x65 = 0;
    x68 = 0;
    x69 = 0;
    x72 = 0;
    x74 = -x67;
    result[1] = rt1 - x4;
    result[16] = 0;
    result[2] = rt2 - x13;
    result[20] = 0;
  };
  if (x16 < x17)
  {
    x60 = -x47;
    x61 = -x48;
    x62 = -x53;
    x63 = -x54 + rhot2 * x31;
    x64 = -x52;
    x65 = -x56;
    x68 = -x49;
    x69 = -x50;
    x72 = -x55;
    x74 = x66 - x57;
    result[1] = rt1 - x73;
    result[16] = 1 - x76;
    result[2] = rt2 - x70;
    result[20] = 1 - x75;
  };
  result[0] = x77 - x7;
  result[3] = x78;
  result[6] = 0;
  result[9] = 0;
  result[12] = 1.0 / 2.0 - x23;
  result[15] = 0;
  result[18] = 0;
  result[4] = x62;
  result[7] = x74;
  result[10] = x61;
  result[13] = x64;
  result[19] = x68;
  result[5] = x72;
  result[8] = x60;
  result[11] = x63;
  result[14] = x65;
  result[17] = x69;
}
void frictionContact3D_AlartCurnierFGenerated(
  double rn,
  double rt1,
  double rt2,
  double un,
  double ut1,
  double ut2,
  double mu,
  double rhon,
  double rhot1,
  double rhot2,
  double *result)
{
  double x0 = rhot1 * ut1;
  double x1 = rt1 - x0;
  double x2 = rn / 2;
  double x3 = rhon * un;
  double x4 = rn - x3;
  double x5 = fabs(x4);
  double x6 = x5 / 2;
  double x7 = x3 / 2;
  double x8 = x2 + x6;
  double x9 = x8 - x7;
  double x10 = x1 * x1;
  double x11 = rhot2 * ut2;
  double x12 = rt2 - x11;
  double x13 = x12 * x12;
  double x14 = x10 + x13;
  double x15 = mu * x9;
  double x16 = sqrt(x14);
  double x18 = 1. / (x16);
  double x20 = x15 * x18;
  double x21 = x12 * x20;
  double x22 = x1 * x20;
  double x23 = x2 + x7;
  if (x16 <= x15)
  {
    result[1] = rt1 - x1;
    result[2] = rt2 - x12;
  };
  if (x15 < x16)
  {
    result[1] = rt1 - x22;
    result[2] = rt2 - x21;
  };
  result[0] = x23 - x6;
}
void frictionContact3D_AlartCurnierABGenerated(
  double rn,
  double rt1,
  double rt2,
  double un,
  double ut1,
  double ut2,
  double mu,
  double rhon,
  double rhot1,
  double rhot2,
  double *result)
{
  double x0 = rhon * un;
  double x1 = rn - x0;
  double x2 = sign(x1);
  double x3 = rhot1 * ut1;
  double x4 = rt1 - x3;
  double x5 = x4 * x4;
  double x6 = rhot2 * ut2;
  double x7 = rt2 - x6;
  double x8 = x7 * x7;
  double x9 = x5 + x8;
  double x10 = rn / 2;
  double x11 = fabs(x1);
  double x12 = x11 / 2;
  double x13 = x10 + x12;
  double x14 = x0 / 2;
  double x15 = x13 - x14;
  double x16 = mu * x15;
  double x17 = sqrt(x9);
  double x19 = 1. / (x17);
  double x22 = x19 * x19 * x19;
  double x23 = x2 / 2;
  double x24 = -rhon / 2;
  double x25 = x2 * x24;
  double x26 = x24 + x25;
  double x27 = 1.0 / 2.0 + x23;
  double x28 = x3 - rt1;
  double x30 = x6 - rt2;
  double x31 = x16 * x19;
  double x32 = mu * x15 * x22;
  double x33 = mu * x4;
  double x34 = x32 * x4;
  double x35 = x32 * x7;
  double x36 = mu * x19;
  double x37 = x16 * x22 * x4 * x7;
  double x38 = x16 * x22;
  double x39 = x19 * x33;
  double x40 = rhot1 * x38;
  double x41 = x36 * x7;
  double x42 = rhot2 * x38;
  double x43 = rhot1 * x37;
  double x44 = rhot2 * x37;
  double x45 = x30 * x34;
  double x46 = x30 * x35;
  double x47 = x28 * x35;
  double x48 = x28 * x34;
  double x49 = x26 * x41;
  double x50 = x27 * x39;
  double x51 = x26 * x39;
  double x52 = x42 * x8;
  double x53 = x40 * x5;
  double x54 = x27 * x41;
  double x55 = -rhot2;
  double x56 = rhon / 2;
  double x57;
  double x58;
  double x59;
  double x60;
  double x61 = rhot2 * x31;
  double x62;
  double x63 = rhot1 * x31;
  double x64;
  double x65 = rhon * x23;
  double x66;
  double x67;
  double x68;
  double x69;
  double x70 = x31 + x48;
  double x71 = x31 + x46;
  double x72 = x56 + x65;
  if (x17 <= x16)
  {
    x57 = 0;
    x58 = 0;
    x59 = 0;
    x60 = 0;
    x62 = 0;
    x64 = rhot1;
    x66 = 0;
    x67 = 0;
    x68 = -x55;
    x69 = 0;
    result[13] = 0;
    result[17] = 0;
  };
  if (x16 < x17)
  {
    x57 = -x51;
    x58 = -x47;
    x59 = -x49;
    x60 = -x43;
    x62 = -x50;
    x64 = x63 - x53;
    x66 = -x54;
    x67 = -x44;
    x68 = x61 - x52;
    x69 = -x45;
    result[13] = 1 - x70;
    result[17] = 1 - x71;
  };
  result[0] = x72;
  result[3] = 0;
  result[6] = 0;
  result[9] = 1.0 / 2.0 - x23;
  result[12] = 0;
  result[15] = 0;
  result[1] = x57;
  result[4] = x64;
  result[7] = x67;
  result[10] = x62;
  result[16] = x69;
  result[2] = x59;
  result[5] = x60;
  result[8] = x68;
  result[11] = x66;
  result[14] = x58;
}

void frictionContact3D_AlartCurnierCKPSFABGenerated(
  double rn,
  double rt1,
  double rt2,
  double un,
  double ut1,
  double ut2,
  double mu,
  double rhon,
  double rhot1,
  double rhot2,
  double *result)
{
  double x0 = rhon * un;
  double x1 = rn - x0;
  double x2 = sign(x1);
  double x3 = rhot1 * ut1;
  double x4 = rt1 - x3;
  double x5 = x4 * x4;
  double x6 = rhot2 * ut2;
  double x7 = rt2 - x6;
  double x8 = x7 * x7;
  double x9 = x5 + x8;
  double x10 = mu * rn;
  double x11 = sqrt(x9);
  double x13 = 1. / (x11);
  double x15 = x13 * x13 * x13;
  double x17 = x3 - rt1;
  double x19 = x10 * x13;
  double x20 = x6 - rt2;
  double x21 = mu * rn * x15;
  double x22 = x20 * x21;
  double x23 = x17 * x21;
  double x24 = x10 * x15 * x4 * x7;
  double x25 = x10 * x15;
  double x26 = mu * x4;
  double x27 = rhot2 * x25;
  double x28 = mu * x7;
  double x29 = rhot1 * x25;
  double x30 = x22 * x4;
  double x31 = x22 * x7;
  double x32 = x23 * x7;
  double x33 = rhot1 * x24;
  double x34 = x23 * x4;
  double x35 = rhot2 * x24;
  double x36 = x27 * x8;
  double x37 = rhon / 2;
  double x38 = x2 / 2;
  double x39 = x29 * x5;
  double x40 = x13 * x28;
  double x41 = x2 * x37;
  double x42 = x13 * x26;
  double x43 = -rhot2;
  double x44;
  double x45 = rhot2 * x19;
  double x46 = x19 * x7;
  double x47;
  double x48;
  double x49;
  double x50 = x19 * x4;
  double x51 = x0 / 2;
  double x52;
  double x53;
  double x54 = fabs(x1) / 2;
  double x55 = rn / 2;
  double x56;
  double x57;
  double x58 = x51 + x55;
  double x59 = x19 + x31;
  double x60 = x37 + x41;
  double x61 = x19 + x34;
  if (x11 <= x10)
  {
    x44 = 0;
    x47 = 0;
    x48 = rhot1;
    x49 = 0;
    x52 = 0;
    x53 = 0;
    x56 = 0;
    x57 = -x43;
    result[1] = rt1 - x4;
    result[16] = 0;
    result[2] = rt2 - x7;
    result[20] = 0;
  };
  if (x10 < x11)
  {
    x44 = -x32;
    x47 = -x33;
    x48 = -x39 + rhot1 * x19;
    x49 = -x35;
    x52 = -x30;
    x53 = -x42;
    x56 = -x40;
    x57 = x45 - x36;
    result[1] = rt1 - x50;
    result[16] = 1 - x61;
    result[2] = rt2 - x46;
    result[20] = 1 - x59;
  };
  result[0] = x58 - x54;
  result[3] = x60;
  result[6] = 0;
  result[9] = 0;
  result[12] = 1.0 / 2.0 - x38;
  result[15] = 0;
  result[18] = 0;
  result[4] = 0;
  result[7] = x48;
  result[10] = x49;
  result[13] = x53;
  result[19] = x52;
  result[5] = 0;
  result[8] = x47;
  result[11] = x57;
  result[14] = x56;
  result[17] = x44;
}
void frictionContact3D_AlartCurnierCKPSFGenerated(
  double rn,
  double rt1,
  double rt2,
  double un,
  double ut1,
  double ut2,
  double mu,
  double rhon,
  double rhot1,
  double rhot2,
  double *result)
{
  double x0 = rhot1 * ut1;
  double x1 = rt1 - x0;
  double x2 = x1 * x1;
  double x3 = rhot2 * ut2;
  double x4 = rt2 - x3;
  double x5 = x4 * x4;
  double x6 = x2 + x5;
  double x7 = mu * rn;
  double x8 = sqrt(x6);
  double x10 = 1. / (x8);
  double x12 = x10 * x7;
  double x13 = rhon * un;
  double x14 = x13 / 2;
  double x15 = x1 * x12;
  double x16 = x12 * x4;
  double x17 = fabs(rn - x13) / 2;
  double x18 = rn / 2;
  double x19 = x14 + x18;
  if (x7 < x8)
  {
    result[1] = rt1 - x15;
    result[2] = rt2 - x16;
  };
  if (x8 <= x7)
  {
    result[1] = rt1 - x1;
    result[2] = rt2 - x4;
  };
  result[0] = x19 - x17;
}
void frictionContact3D_AlartCurnierCKPSABGenerated(
  double rn,
  double rt1,
  double rt2,
  double un,
  double ut1,
  double ut2,
  double mu,
  double rhon,
  double rhot1,
  double rhot2,
  double *result)
{
  double x0 = rhon * un;
  double x1 = rn - x0;
  double x2 = sign(x1);
  double x3 = rhot1 * ut1;
  double x4 = rt1 - x3;
  double x5 = x4 * x4;
  double x6 = rhot2 * ut2;
  double x7 = rt2 - x6;
  double x8 = x7 * x7;
  double x9 = x5 + x8;
  double x10 = mu * rn;
  double x11 = sqrt(x9);
  double x13 = 1. / (x11 * x11 * x11);
  double x16 = pow(x13, (1.0 / 3.0));
  double x17 = x3 - rt1;
  double x19 = x10 * x16;
  double x20 = x6 - rt2;
  double x21 = mu * rn * x13;
  double x22 = x13 * x4;
  double x23 = x10 * x13;
  double x24 = x20 * x21;
  double x25 = x17 * x21;
  double x26 = x10 * x22 * x7;
  double x27 = mu * x4;
  double x28 = rhot2 * x23;
  double x29 = rhot1 * x23;
  double x30 = mu * x7;
  double x31 = x24 * x4;
  double x32 = x24 * x7;
  double x33 = x25 * x4;
  double x34 = rhot2 * x26;
  double x35 = rhot1 * x26;
  double x36 = x25 * x7;
  double x37 = rhon / 2;
  double x38 = x2 / 2;
  double x39 = x29 * x5;
  double x40 = x28 * x8;
  double x41 = x2 * x37;
  double x42 = x16 * x30;
  double x43 = x16 * x27;
  double x44 = -rhot2;
  double x45;
  double x46 = rhot2 * x19;
  double x47;
  double x48;
  double x49 = -rhot1;
  double x50;
  double x51;
  double x52;
  double x53;
  double x54;
  double x55 = x19 + x32;
  double x56 = x37 + x41;
  double x57 = x19 + x33;
  if (x11 <= x10)
  {
    x45 = 0;
    x47 = 0;
    x48 = -x44;
    x50 = 0;
    x51 = 0;
    x52 = 0;
    x53 = -x49;
    x54 = 0;
    result[13] = 0;
    result[17] = 0;
  };
  if (x10 < x11)
  {
    x45 = -x36;
    x47 = -x35;
    x48 = x46 - x40;
    x50 = -x43;
    x51 = -x31;
    x52 = -x34;
    x53 = -x39 + rhot1 * x19;
    x54 = -x42;
    result[13] = 1 - x57;
    result[17] = 1 - x55;
  };
  result[0] = x56;
  result[3] = 0;
  result[6] = 0;
  result[9] = 1.0 / 2.0 - x38;
  result[12] = 0;
  result[15] = 0;
  result[1] = 0;
  result[4] = x53;
  result[7] = x52;
  result[10] = x50;
  result[16] = x51;
  result[2] = 0;
  result[5] = x47;
  result[8] = x48;
  result[11] = x54;
  result[14] = x45;
}

void frictionContact3D_localAlartCurnierFunctionGenerated(
  double *reaction,
  double *velocity,
  double mu,
  double *rho,
  double *f,
  double *A,
  double *B)
{
  double result[21]; //3 + 2 * 9

  assert(reaction);
  assert(velocity);
  assert(rho);

  SET3(reaction);
  SET3(velocity);
  SET3(rho);


  if (f && A && B)
  {

    frictionContact3D_AlartCurnierFABGenerated(
      *reaction0, *reaction1, *reaction2,
      *velocity0, *velocity1, *velocity2,
      mu,
      *rho0, *rho1, *rho2,
      result);
    cpy3(result, f);
    cpy3x3(result + 3, A);
    cpy3x3(result + 12, B);
  }

  else
  {
    if (f)
    {
      frictionContact3D_AlartCurnierFGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3(result, f);
    }

    if (A && B)
    {
      frictionContact3D_AlartCurnierABGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3x3(result, A);
      cpy3x3(result + 9, B);
    }
  }
}

void frictionContact3D_localAlartCurnierCKPSFunctionGenerated(
  double *reaction,
  double *velocity,
  double mu,
  double *rho,
  double *f,
  double *A,
  double *B)
{
  double result[21]; //3 + 2 * 9

  assert(reaction);
  assert(velocity);
  assert(rho);

  SET3(reaction);
  SET3(velocity);
  SET3(rho);


  if (f && A && B)
  {

    frictionContact3D_AlartCurnierCKPSFABGenerated(
      *reaction0, *reaction1, *reaction2,
      *velocity0, *velocity1, *velocity2,
      mu,
      *rho0, *rho1, *rho2,
      result);
    cpy3(result, f);
    cpy3x3(result + 3, A);
    cpy3x3(result + 12, B);
  }

  else
  {
    if (f)
    {
      frictionContact3D_AlartCurnierCKPSFGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3(result, f);
    }

    if (A && B)
    {
      frictionContact3D_AlartCurnierCKPSABGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3x3(result, A);
      cpy3x3(result + 9, B);
    }
  }
}


void frictionContact3D_globalAlartCurnierFunctionGenerated(
  unsigned int problemSize,
  double *reaction,
  double *velocity,
  double *mu,
  double *rho,
  double *result,
  double *A,
  double *B)
{
  assert(reaction);
  assert(velocity);
  assert(rho);
  assert(mu);

  assert(problemSize / 3 > 0);
  assert(problemSize % 3 == 0);

  unsigned int i;
  for (i = 0; i < problemSize; i += 3)
  {

    frictionContact3D_localAlartCurnierFunctionGenerated(reaction,
        velocity,
        *mu,
        rho,
        result, A, B);

    reaction += 3;
    velocity += 3;
    mu++;
    rho += 3;

    if (result)
      result += 3;

    if (A)
      A += 9;

    if (B)
      B += 9;

  }

}

void frictionContact3D_localAlartCurnierFunctionHandMade(
  double *reaction,
  double *velocity,
  double mu,
  double *rho,
  double *ACF,
  double *A,
  double *B)
{
  assert(0);

  /*
  assert (reaction);
  assert (velocity);
  assert (rho);

  SET3(reaction);
  SET3(velocity);
  SET3(rho);

  SET3X3MAYBE(A);
  SET3X3MAYBE(B);
  SET3MAYBE(ACF);

  if (A && B)
  {
  *A01 = 0.;
  *A02 = 0.;
  *B01 = 0.;
  *B02 = 0.;
  }


  double D0,D1,D2,muD0;

  D0 = *reaction0 - *rho0 * *velocity0;
  D1 = *reaction1 - *rho1 * *velocity1;
  D2 = *reaction2 - *rho2 * *velocity2;

  muD0 = mu*D0;

  double hypotD1D2 = hypot(D1,D2);


  if (muD0 <= 0.)
  {

  if (ACF)
  {
    *ACF1 = *reaction1;
    *ACF2 = *reaction2;
  }

  if (A && B)
  {

    *A10 = 0.;
    *A11 = 0.;
    *A12 = 0.;

    *A20 = 0.;
    *A21 = 0.;
    *A22 = 0.;

    *B10 = 0.;
    *B11 = 1.;
    *B12 = 0.;

    *B20 = 0.;
    *B21 = 0.;
    *B22 = 1.;
  }


  };

  if (0<muD0 && hypotD1D2<=muD0)
  {
  if (ACF)
  {
    *ACF1 = *reaction1 - D1;
    *ACF2 = *reaction2 - D2;
  }

  if (A && B)
  {
    *A10 = 0.;
    *A11 = *rho1;
    *A12 = 0.;

    *A20 = 0.;
    *A21 = 0.;
    *A22 = *rho2;

    *B10 = 0.;
    *B11 = 0.;
    *B12 = 0.;

    *B20 = 0.;
    *B21 = 0.;
    *B22 = 0.;
  }

  }

  if (D0<0.)
  {

  if (ACF)
  {
    *ACF0 = *reaction0;
  }

  if (A && B)
  {
    *A00 = 0.;
    *B00 = 1.;
  }

  }

  if (D0>=0.)
  {
  if (ACF)
  {
    *ACF0 = *reaction0 - D0;
  }

  if (A && B)
  {
    *A00 = *rho0;
    *B00 = 0.;
  }

  }

  if (0 < muD0 && muD0 < hypotD1D2)
  {
  double cubehypotD1D2 = hypotD1D2 * hypotD1D2 * hypotD1D2;
  double muD0rho1 = muD0* *rho1;

  if (ACF)
  {
    *ACF1 = *reaction1 - muD0*D1/hypotD1D2;
    *ACF2 = *reaction2 - muD0*D2/hypotD1D2;
  }

  if (A && B)
  {
    *A10 = mu*D1* *rho0/hypotD1D2;
    *A11 = -muD0rho1*D1*D1/cubehypotD1D2 + muD0rho1/hypotD1D2;
    *A12 = - D0*D1*D2*mu**rho2/cubehypotD1D2;

    *A20 = D2*mu* *rho0/hypotD1D2;
    *A21 = - D0*D1*D2*mu* *rho1/cubehypotD1D2;
    *A22 = - muD0* *rho2 *D2*D2/cubehypotD1D2 + muD0* *rho2/hypotD1D2;

    *B10 = - mu*D1/hypotD1D2;
    *B11 = 1 + muD0*D1*D1/cubehypotD1D2 - muD0/hypotD1D2;
    *B12 =  muD0*D1*D2/cubehypotD1D2;

    *B20 = - mu*D2/hypotD1D2;
    *B21 =  muD0*D1*D2/cubehypotD1D2;
    *B22 = 1 + muD0*D2*D2/cubehypotD1D2 - muD0/hypotD1D2;
  }

  }
  */
};


void frictionContact3D_globalAlartCurnierFunctionHandMade(
  unsigned int problemSize,
  double *reaction,
  double *velocity,
  double *mu,
  double *rho,
  double *result,
  double *A,
  double *B)
{
  assert(reaction);
  assert(velocity);
  assert(rho);
  assert(mu);


  assert(problemSize / 3 > 0);
  assert(problemSize % 3 == 0);

  unsigned int i;
  for (i = 0; i < problemSize; i += 3)
  {

    frictionContact3D_localAlartCurnierFunctionHandMade(reaction,
        velocity,
        *mu,
        rho,
        result, A, B);

    scal3(-1., result);


    /*
        computeAlartCurnierSTD(reaction,
                               velocity,
                               *mu,
                               rho,
                               result,A,B);
    */

    /* generated function = hand made */
#ifndef NDEBUG
    double result_g[3];
    double A_g[9];
    double B_g[9];

    frictionContact3D_localAlartCurnierFunctionGenerated(reaction,
        velocity,
        *mu,
        rho,
        result_g, A_g, B_g);
    if (result)
      sub3(result, result_g);

    if (A)
      sub3x3(A, A_g);

    if (B)
      sub3x3(B, B_g);


    assert(hypot3(result_g) < 1e-7);
    assert(hypot9(A_g) < 1e-7);
    assert(hypot9(B_g) < 1e-7);
#endif

    reaction += 3;
    velocity += 3;
    mu++;
    rho += 3;
    result += 3;
    A += 9;
    B += 9;

  }

}

void computeAWpB(
  unsigned int problemSize,
  double *A,
  double *W,
  double *B,
  double *result)
{
  assert(problemSize >= 3);

  double Wij[9], Ai[9], Bi[9], tmp[9];

  for (unsigned int ip3 = 0, ip9 = 0; ip3 < problemSize; ip3 += 3, ip9 += 9)
  {
    assert(ip9 < 3 * problemSize - 8);

    extract3x3(3, ip9, 0, A, Ai);
    extract3x3(3, ip9, 0, B, Bi);

    for (unsigned int jp3 = 0; jp3 < problemSize; jp3 += 3)
    {
      assert(jp3 < problemSize - 2);
      assert(ip3 < problemSize - 2);

      extract3x3(problemSize, ip3, jp3, W, Wij);
      mm3x3(Ai, Wij, tmp);
      if (jp3 == ip3) add3x3(Bi, tmp);
      insert3x3(problemSize, ip3, jp3, result, tmp);

#ifdef VERBOSE_DEBUG_1
      if (jp3 == ip3)
      {
        printf("Ai\n");
        print3x3(Ai);

        printf("Bi\n");
        print3x3(Bi);

        printf("Wij");
        print3x3(Wij);

        printf("result\n");
        print3x3(tmp);
      }
#endif

    }
  }
}

int globalLineSearchGP(
  unsigned int problemSize,
  double *reaction,
  double *velocity,
  double *mu,
  double *rho,
  double *F,
  double *A,
  double *B,
  double *W,
  double *qfree,
  double *AWpB,
  double *direction,
  double *tmp,
  double alpha[1],
  unsigned int maxiter_ls)
{
  double inf = 1e10;
  double alphamin = 0.0;
  double alphamax = inf;

  double m1 = 0.01, m2 = 0.99;

  // Computation of q(t) and q'(t) for t =0

  double q0 = 0.5 * DDOT(problemSize, F, 1, F, 1);

  // useless (already computed)
  computeAWpB(problemSize, A, W, B, AWpB);

  //  tmp <- AWpB * direction
  DGEMV(LA_NOTRANS, problemSize, problemSize, 1., AWpB, problemSize, direction, 1, 0., tmp, 1);

  double dqdt0 = DDOT(problemSize, F, 1, tmp, 1);

  for (int iter = 0; iter < maxiter_ls; ++iter)
  {

    // tmp <- alpha*direction+reaction
    DCOPY(problemSize, reaction, 1, tmp, 1);
    DAXPY(problemSize, alpha[0], direction, 1, tmp, 1);

    // velocity <- W*tmp + qfree
    DCOPY(problemSize, qfree, 1, velocity, 1);
    DGEMV(LA_NOTRANS, problemSize, problemSize, 1.,
          W, problemSize, tmp, 1, 1., velocity, 1);

    frictionContact3D_globalAlartCurnierFunctionGenerated(problemSize, tmp,
        velocity, mu, rho, F,
        NULL, NULL);

    double q  = 0.5 * DDOT(problemSize, F, 1, F, 1);

    assert(q >= 0);

    double slope = (q - q0) / alpha[0];

    int C1 = (slope >= m2 * dqdt0);
    int C2 = (slope <= m1 * dqdt0);

    if (C1 && C2)
    {
      if (verbose > 1)
      {
        printf("global line search success. Number of iteration = %i  alpha = %.10e, q = %.10e\n", iter, alpha[0], q);
      }

      return 0;

    }
    else if (!C1)
    {
      alphamin = alpha[0];
    }
    else
    {
      // not(C2)
      alphamax = alpha[0];
    }

    if (alpha[0] < inf)
    {
      alpha[0] = 0.5 * (alphamin + alphamax);
    }
    else
    {
      alpha[0] = alphamin;
    }

  }
  if (verbose > 1)
  {
    printf("global line search failed. max number of iteration reached  = %i  with alpha = %.10e \n", maxiter_ls, alpha[0]);
  }

  return -1;
}

void frictionContact3D_globalAlartCurnier(
  FrictionContactProblem* problem,
  double *reaction,
  double *velocity,
  int *info,
  SolverOptions *options)
{
  assert(problem);
  assert(reaction);
  assert(velocity);
  assert(info);
  assert(options);

  assert(problem->dimension == 3);

  assert(options->iparam);
  assert(options->dparam);

  assert(problem->q);
  assert(problem->mu);
  assert(problem->M);
  assert(problem->M->matrix0);

  unsigned int problemSize = 3 * problem->numberOfContacts;

  unsigned int iter = 0;
  unsigned int itermax = options->iparam[0];
  unsigned int erritermax = options->iparam[1];

  assert(itermax > 0);

  double tolerance = options->dparam[0];
  assert(tolerance > 0);

  unsigned int problemSize2 = problemSize * problemSize;
  unsigned int _3problemSize = 3 * problemSize;

  void *buffer;

  if (!options->dWork)
    buffer = malloc((14 * problemSize +
                     problemSize2) * sizeof(double) +
                    problemSize * sizeof(int));
  else
    buffer = options->dWork;

  double *F = (double *) buffer; //malloc(problemSize*sizeof(double));
  double *tmp1 = (double *) F + problemSize; //malloc(problemSize*sizeof(double));
  double *tmp2 = (double *) tmp1 + problemSize; //malloc(problemSize*sizeof(double));
  double *A = tmp2 + problemSize; //malloc(3*problemSize*sizeof(double));
  double *B = A + _3problemSize; //malloc(3*problemSize*sizeof(double));
  double *rho = B + _3problemSize; //malloc(problemSize*sizeof(double));
  double *AWpB = rho + problemSize;// malloc(problemSize*problemSize*sizeof(double));
  int *ipiv = (int *)(AWpB + problemSize2);  // malloc(problemSize*sizeof(int));

  double w;

  int LWORK;
  double *WORK = NULL;
  // iparam[2] != 0 => use of DGELS
  if (options->iparam[2])
  {
    int dgelsinfo[1];

    DGELS(problemSize, problemSize,
          1, AWpB, problemSize,
          F, problemSize, &w, -1, dgelsinfo);

    LWORK = (int) w;

    WORK = (double *) malloc(w * sizeof(double));
  }

  for (unsigned int i = 0; i < problemSize; ++i) rho[i] = 1.;

  info[0] = 1;

  // velocity <- M*reaction + qfree
  DCOPY(problemSize, problem->q, 1, velocity, 1);
  DGEMV(LA_NOTRANS, problemSize, problemSize, 1.,
        problem->M->matrix0, problemSize, reaction, 1, 1., velocity, 1);





  while (iter++ < itermax)
  {

    frictionContact3D_globalAlartCurnierFunctionGenerated(problemSize,
        reaction, velocity,
        problem->mu, rho,
        F, A, B);


    // AW + B
    computeAWpB(problemSize, A, problem->M->matrix0, B, AWpB);

    int fail;



    DCOPY(problemSize, F, 1, tmp1, 1);
    DSCAL(problemSize, -1., tmp1, 1);

    if (options->iparam[2])
    {
      assert(WORK);
      DGELS(problemSize, problemSize, 1, AWpB, problemSize,
            tmp1, problemSize, WORK, LWORK, fail);
    }
    else
    {
      DGESV(problemSize, 1, AWpB, problemSize, ipiv,
            tmp1, problemSize, fail);
    }


    assert(fail >= 0);

    if (fail > 0)
      /*if (verbose>0)*/
      printf("GLOBALAC: warning DGESV fail with U(%d,%d) == 0.\n", fail, fail);

    // line search
    double alpha = 1;
    int info_ls = globalLineSearchGP(problemSize, reaction, velocity, problem->mu, rho, F, A, B,
                                     problem->M->matrix0, problem->q, AWpB, tmp1, tmp2, &alpha, 100);

    if (!info_ls)
      DAXPY(problemSize, alpha, tmp1, 1, reaction, 1);
    else
      DAXPY(problemSize, 1, tmp1, 1., reaction, 1);


    // velocity <- M*reaction + qfree
    DCOPY(problemSize, problem->q, 1, velocity, 1);
    DGEMV(LA_NOTRANS, problemSize, problemSize, 1.,
          problem->M->matrix0, problemSize, reaction, 1, 1., velocity, 1);



    options->dparam[1] = INFINITY;

    if (!(iter % erritermax))
    {
      frictionContact3D_globalAlartCurnierFunctionGenerated(problemSize,
          reaction, velocity,
          problem->mu, rho,
          F, NULL, NULL);



      FrictionContact3D_compute_error(problem, reaction, velocity,
                                      tolerance, options, &(options->dparam[1]));


      assert((DNRM2(problemSize, F, 1)
              / (1 + DNRM2(problemSize, problem->q, 1)))
             <= 10 * options->dparam[1]);


    }

    if (verbose > 0)
      printf("GLOBALAC: iteration %d : error=%g\n", iter, options->dparam[1]);

    if (options->dparam[1] < tolerance)
    {
      info[0] = 0;
      break;
    }


  }



  if (verbose > 0)
  {
    if (!info[0])
      printf("GLOBALAC: convergence after %d iterations, error : %g\n",
             iter, options->dparam[1]);
    else
    {
      printf("GLOBALAC: no convergence after %d iterations, error : %g\n",
             iter, options->dparam[1]);
    }
  }

#ifdef DUMP_PROBLEM
  if (info[0])
  {
    static int file_counter = 0;
    char filename[64];
    printf("GLOBALAC: dumping problem\n");
    sprintf(filename, "GLOBALAC_failure%d.dat", file_counter++);
    FILE* file = fopen(filename, "w");
    frictionContact_printInFile(problem, file);
    fclose(file);
  }
#endif

  if (!options->dWork)
  {
    assert(buffer);
    free(buffer);

    if (WORK)
      free(WORK);

  }
  else
  {
    assert(buffer == options->dWork);
  }


}

int frictionContact3D_globalAlartCurnier_setDefaultSolverOptions(
  SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the default solver options for the GLOBALAC Solver\n");
  }

  options->solverId = SICONOS_FRICTION_3D_GLOBALAC;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *) malloc(options->iSize * sizeof(int));
  options->dparam = (double *) malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  for (unsigned int i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 200;
  options->iparam[1] = 1;
  options->dparam[0] = 1e-3;

  options->internalSolvers = NULL;

  return 0;
}
