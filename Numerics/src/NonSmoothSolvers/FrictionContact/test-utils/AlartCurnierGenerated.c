#include <assert.h>
#include <math.h>
#include "op3x3.h"

#define sign(x) copysign(1.,x)

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
  double x5 = rhot1 * ut1;
  double x6 = rt1 - x5;
  double x8 = NAN;
  double x9 = x6 * x6;
  double x10 = rhot2 * ut2;
  double x11 = rt2 - x10;
  double x12 = x11 * x11;
  double x13 = x12 + x9;
  double x14 = NAN;
  double x15 = sqrt(x13);
  double x16 = NAN;
  double x17 = 1. / (x15);
  double x18 = -rhon;
  double x20 = NAN;
  double x21 = NAN;
  double x22 = NAN;
  double x23 = x17 * x17 * x17;
  double x25 = NAN;
  double x26 = x5 - rt1;
  double x27 = NAN;
  double x28 = NAN;
  double x29 = x10 - rt2;
  if (x1 <= 0)
  {
    x8 = 0;
    x14 = mu * x8;
    x16 = x15 <= x14;
    x20 = 0;
    x21 = x14 < x15;
    x25 = 0;
    x28 = x14 * x17;
    result[0] = rn - x8;
    result[3] = -x20;
    result[12] = 1 - x25;
  };
  if (0 < x1)
  {
    x8 = x1;
    x14 = mu * x8;
    x16 = x15 <= x14;
    x20 = x18;
    x21 = x14 < x15;
    x25 = 1;
    x28 = x14 * x17;
    result[0] = rn - x8;
    result[3] = -x20;
    result[12] = 1 - x25;
  };
  if (x15 <= x14)
  {
    x22 = 0;
    x27 = 1;
    x28 = x14 * x17;
    result[0] = rn - x8;
    result[1] = rt1 - x6;
    result[4] = 0;
    result[7] = rhot1;
    result[10] = 0;
    result[13] = 0;
    result[16] = 0;
    result[19] = 0;
    result[2] = rt2 - x11;
    result[5] = 0;
    result[8] = 0;
    result[11] = rhot2;
    result[14] = 0;
    result[17] = 0;
    result[20] = 0;
  };
  if ((x1 <= 0) && (x15 <= x14))
  {
    result[1] = rt1 - x6;
    result[4] = 0;
    result[7] = rhot1;
    result[10] = 0;
    result[13] = 0;
    result[16] = 0;
    result[19] = 0;
    result[2] = rt2 - x11;
    result[5] = 0;
    result[8] = 0;
    result[11] = rhot2;
    result[14] = 0;
    result[17] = 0;
    result[20] = 0;
  };
  if ((x1 <= 0) && (x14 < x15))
  {
    result[1] = rt1 - x28 * x6;
    result[4] = -mu * x17 * x20 * x6;
    result[7] = rhot1 * x28 - rhot1 * x14 * x23 * x9;
    result[10] = -rhot2 * x11 * x14 * x23 * x6;
    result[13] = -mu * x17 * x25 * x6;
    result[16] = 1 - x28 - mu * x23 * x26 * x6 * x8;
    result[19] = -mu * x23 * x29 * x6 * x8;
    result[2] = rt2 - x11 * x28;
    result[5] = -mu * x11 * x17 * x20;
    result[8] = -rhot1 * x11 * x14 * x23 * x6;
    result[11] = rhot2 * x28 - rhot2 * x12 * x14 * x23;
    result[14] = -mu * x11 * x17 * x25;
    result[17] = -mu * x11 * x23 * x26 * x8;
    result[20] = 1 - x28 - mu * x11 * x23 * x29 * x8;
  };
  if ((0 < x1) && (x15 <= x14))
  {
    result[1] = rt1 - x6;
    result[4] = 0;
    result[7] = rhot1;
    result[10] = 0;
    result[13] = 0;
    result[16] = 0;
    result[19] = 0;
    result[2] = rt2 - x11;
    result[5] = 0;
    result[8] = 0;
    result[11] = rhot2;
    result[14] = 0;
    result[17] = 0;
    result[20] = 0;
  };
  if ((0 < x1) && (x14 < x15))
  {
    result[1] = rt1 - x28 * x6;
    result[4] = -mu * x17 * x20 * x6;
    result[7] = rhot1 * x28 - rhot1 * x14 * x23 * x9;
    result[10] = -rhot2 * x11 * x14 * x23 * x6;
    result[13] = -mu * x17 * x25 * x6;
    result[16] = 1 - x28 - mu * x23 * x26 * x6 * x8;
    result[19] = -mu * x23 * x29 * x6 * x8;
    result[2] = rt2 - x11 * x28;
    result[5] = -mu * x11 * x17 * x20;
    result[8] = -rhot1 * x11 * x14 * x23 * x6;
    result[11] = rhot2 * x28 - rhot2 * x12 * x14 * x23;
    result[14] = -mu * x11 * x17 * x25;
    result[17] = -mu * x11 * x23 * x26 * x8;
    result[20] = 1 - x28 - mu * x11 * x23 * x29 * x8;
  };
  result[6] = 0;
  result[9] = 0;
  result[15] = 0;
  result[18] = 0;
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
  double x0 = rhon * un;
  double x1 = rn - x0;
  double x2 = rhot1 * ut1;
  double x3 = rt1 - x2;
  double x8 = NAN;
  double x9 = x3 * x3;
  double x10 = rhot2 * ut2;
  double x11 = rt2 - x10;
  double x12 = x11 * x11;
  double x13 = x12 + x9;
  double x14 = NAN;
  double x15 = sqrt(x13);
  double x16 = NAN;
  double x17 = 1. / (x15);
  double x18 = NAN;
  if (x1 <= 0)
  {
    x8 = 0;
    x14 = mu * x8;
    x16 = x15 <= x14;
    x18 = x14 < x15;
    result[0] = rn - x8;
  };
  if (0 < x1)
  {
    x8 = x1;
    x14 = mu * x8;
    x16 = x15 <= x14;
    x18 = x14 < x15;
    result[0] = rn - x8;
  };
  if ((x1 <= 0) && (x15 <= x14))
  {
    result[1] = rt1 - x3;
    result[2] = rt2 - x11;
  };
  if ((x1 <= 0) && (x14 < x15))
  {
    result[1] = rt1 - x14 * x17 * x3;
    result[2] = rt2 - x11 * x14 * x17;
  };
  if ((0 < x1) && (x15 <= x14))
  {
    result[1] = rt1 - x3;
    result[2] = rt2 - x11;
  };
  if ((0 < x1) && (x14 < x15))
  {
    result[1] = rt1 - x14 * x17 * x3;
    result[2] = rt2 - x11 * x14 * x17;
  };
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
  double x5 = rhot1 * ut1;
  double x6 = rt1 - x5;
  double x7 = x6 * x6;
  double x8 = rhot2 * ut2;
  double x9 = rt2 - x8;
  double x10 = x9 * x9;
  double x11 = x10 + x7;
  double x12 = -rhon;
  double x14 = NAN;
  double x16 = NAN;
  double x17 = NAN;
  double x18 = sqrt(x11);
  double x19 = NAN;
  double x20 = 1. / (x18);
  double x21 = NAN;
  double x22 = NAN;
  double x23 = x20 * x20 * x20;
  double x25 = NAN;
  double x26 = x5 - rt1;
  double x27 = NAN;
  double x28 = NAN;
  double x29 = x8 - rt2;
  if (x1 <= 0)
  {
    x14 = 0;
    x16 = 0;
    x17 = mu * x16;
    x19 = x18 <= x17;
    x21 = x17 < x18;
    x25 = 0;
    x28 = x17 * x20;
    result[0] = -x14;
    result[9] = 1 - x25;
  };
  if (0 < x1)
  {
    x14 = x12;
    x16 = x1;
    x17 = mu * x16;
    x19 = x18 <= x17;
    x21 = x17 < x18;
    x25 = 1;
    x28 = x17 * x20;
    result[0] = -x14;
    result[9] = 1 - x25;
  };
  if (x18 <= x17)
  {
    x22 = 0;
    x27 = 1;
    x28 = x17 * x20;
    result[1] = 0;
    result[4] = rhot1;
    result[7] = 0;
    result[10] = 0;
    result[13] = 0;
    result[16] = 0;
    result[2] = 0;
    result[5] = 0;
    result[8] = rhot2;
    result[11] = 0;
    result[14] = 0;
    result[17] = 0;
  };
  if ((x1 <= 0) && (x18 <= x17))
  {
    result[1] = 0;
    result[4] = rhot1;
    result[7] = 0;
    result[10] = 0;
    result[13] = 0;
    result[16] = 0;
    result[2] = 0;
    result[5] = 0;
    result[8] = rhot2;
    result[11] = 0;
    result[14] = 0;
    result[17] = 0;
  };
  if ((x1 <= 0) && (x17 < x18))
  {
    result[1] = -mu * x14 * x20 * x6;
    result[4] = rhot1 * x28 - rhot1 * x17 * x23 * x7;
    result[7] = -rhot2 * x17 * x23 * x6 * x9;
    result[10] = -mu * x20 * x25 * x6;
    result[13] = 1 - x28 - mu * x16 * x23 * x26 * x6;
    result[16] = -mu * x16 * x23 * x29 * x6;
    result[2] = -mu * x14 * x20 * x9;
    result[5] = -rhot1 * x17 * x23 * x6 * x9;
    result[8] = rhot2 * x28 - rhot2 * x10 * x17 * x23;
    result[11] = -mu * x20 * x25 * x9;
    result[14] = -mu * x16 * x23 * x26 * x9;
    result[17] = 1 - x28 - mu * x16 * x23 * x29 * x9;
  };
  if ((0 < x1) && (x18 <= x17))
  {
    result[1] = 0;
    result[4] = rhot1;
    result[7] = 0;
    result[10] = 0;
    result[13] = 0;
    result[16] = 0;
    result[2] = 0;
    result[5] = 0;
    result[8] = rhot2;
    result[11] = 0;
    result[14] = 0;
    result[17] = 0;
  };
  if ((0 < x1) && (x17 < x18))
  {
    result[1] = -mu * x14 * x20 * x6;
    result[4] = rhot1 * x28 - rhot1 * x17 * x23 * x7;
    result[7] = -rhot2 * x17 * x23 * x6 * x9;
    result[10] = -mu * x20 * x25 * x6;
    result[13] = 1 - x28 - mu * x16 * x23 * x26 * x6;
    result[16] = -mu * x16 * x23 * x29 * x6;
    result[2] = -mu * x14 * x20 * x9;
    result[5] = -rhot1 * x17 * x23 * x6 * x9;
    result[8] = rhot2 * x28 - rhot2 * x10 * x17 * x23;
    result[11] = -mu * x20 * x25 * x9;
    result[14] = -mu * x16 * x23 * x26 * x9;
    result[17] = 1 - x28 - mu * x16 * x23 * x29 * x9;
  };
  result[3] = 0;
  result[6] = 0;
  result[12] = 0;
  result[15] = 0;
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
  double x44 = NAN;
  double x45 = rhot2 * x19;
  double x46 = x19 * x7;
  double x47 = NAN;
  double x48 = NAN;
  double x49 = NAN;
  double x50 = x19 * x4;
  double x51 = x0 / 2;
  double x52 = NAN;
  double x53 = NAN;
  double x54 = fabs(x1) / 2;
  double x55 = rn / 2;
  double x56 = NAN;
  double x57 = NAN;
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
  double x45 = NAN;
  double x46 = rhot2 * x19;
  double x47 = NAN;
  double x48 = NAN;
  double x49 = -rhot1;
  double x50 = NAN;
  double x51 = NAN;
  double x52 = NAN;
  double x53 = NAN;
  double x54 = NAN;
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
