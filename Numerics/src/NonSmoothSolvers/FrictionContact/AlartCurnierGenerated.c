#include <assert.h>
#include <math.h>
#include "op3x3.h"
#include "AlartCurnierGenerated.h"
#define sign(x) copysign(1.,x)

// this file consists of generated code
//#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

void frictionContact3D_AlartCurnierJeanMoreauFABGenerated(
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
  double x0=-rhon*un + rn;
  int    x1=x0 <= 0;
  int    x2=x0 > 0;
  double x3=rhot1*ut1;
  double x4=rt1 - x3;
  double x5=x4*x4;
  double x6=rhot2*ut2;
  double x7=rt2 - x6;
  double x8=x7*x7;
  double x9=x5 + x8;
  double x10=sqrt(x9);
  int    x11=rn <= 0;
  int    x12=rn > 0;
  double x13;
  double x14;
  int x15 = 0;
  double x16=1./(x10);
  double x17;
  int x18 = 0;
  double x19=sqrt(1./(x9*x9*x9));
  double x20;
  double x21=x4*x7;
  double x22;
  double x23;
  double x24=-rt1 + x3;
  double x25;
  double x26=-rt2 + x6;
  double x27;
  if (x11)
  {
    x13=0;
    x23=0;
  };
  if (x12)
  {
    x13=rn;
    x23=mu*x16;
  };
  if (x11 || x12)
  {
    x14=mu*x13;
    x15=x10 <= x14;
    x17=x14*x16;
    x18=x10 > x14;
    x20=mu*rhot1*x13*x19;
    x22=mu*rhot2*x13*x19;
    x25=mu*x13*x19*x4;
    x27=mu*x13*x19*x7;
  };
  if (x1)
  {
    result[0]=rn;
    result[3]=0;
    result[12]=1;
  };
  if (x2)
  {
    result[0]=rn - x0;
    result[3]=rhon;
    result[12]=0;
  };
  if ((x11 || x12) && (x15))
  {
    result[1]=rt1 - x4;
    result[7]=rhot1;
    result[10]=0;
    result[13]=0;
    result[16]=0;
    result[19]=0;
    result[2]=rt2 - x7;
    result[8]=0;
    result[11]=rhot2;
    result[14]=0;
    result[17]=0;
    result[20]=0;
  };
  if ((x11 || x12) && (x18))
  {
    result[1]=rt1 - x17*x4;
    result[7]=rhot1*x17 - x20*x5;
    result[10]=-x21*x22;
    result[13]=-x23*x4;
    result[16]=-x17 - x24*x25 + 1;
    result[19]=-x25*x26;
    result[2]=rt2 - x17*x7;
    result[8]=-x20*x21;
    result[11]=rhot2*x17 - x22*x8;
    result[14]=-x23*x7;
    result[17]=-x24*x27;
    result[20]=-x17 - x26*x27 + 1;
  };
  result[6]=0;
  result[9]=0;
  result[15]=0;
  result[18]=0;
  result[4]=0;
  result[5]=0;
}
void frictionContact3D_AlartCurnierJeanMoreauFGenerated(
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
  double x0=-rhon*un + rn;
  double x1=-rhot1*ut1 + rt1;
  double x2=-rhot2*ut2 + rt2;
  double x3=sqrt(pow(x1, 2) + pow(x2, 2));
  double x4;
  double x5;
  int x6 = 0;
  double x7;
  int x8 = 0;
  if (rn <= 0)
  {
    x4=0;
  };
  if (rn > 0)
  {
    x4=rn;
  };
  if (rn <= 0 || rn > 0)
  {
    x5=mu*x4;
    x6=x3 <= x5;
    x7=mu*x4/x3;
    x8=x3 > x5;
  };
  if (x0 <= 0)
  {
    result[0]=rn;
  };
  if (x0 > 0)
  {
    result[0]=rn - x0;
  };
  if ((rn <= 0 || rn > 0) && (x6))
  {
    result[1]=rt1 - x1;
    result[2]=rt2 - x2;
  };
  if ((rn <= 0 || rn > 0) && (x8))
  {
    result[1]=rt1 - x1*x7;
    result[2]=rt2 - x2*x7;
  };
}
void frictionContact3D_AlartCurnierJeanMoreauABGenerated(
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
  double x0=-rhon*un + rn;
  int    x1=x0 <= 0;
  int    x2=x0 > 0;
  double x3=rhot1*ut1;
  double x4=rt1 - x3;
  double x5=x4*x4;
  double x6=rhot2*ut2;
  double x7=rt2 - x6;
  double x8=x7*x7;
  double x9=x5 + x8;
  double x10=sqrt(x9);
  int    x11=rn <= 0;
  int    x12=rn > 0;
  double x13;
  double x14;
  int x15 = 0;
  double x16=1./(x10);
  double x17;
  double x18=sqrt(1./(x9*x9*x9));
  double x19;
  int x20 = 0;
  double x21=x4*x7;
  double x22;
  double x23;
  double x24=-rt1 + x3;
  double x25;
  double x26=-rt2 + x6;
  double x27;
  if (x11)
  {
    x13=0;
    x23=0;
  };
  if (x12)
  {
    x13=rn;
    x23=mu*x16;
  };
  if (x11 || x12)
  {
    x14=mu*x13;
    x15=x10 <= x14;
    x17=x14*x16;
    x19=mu*rhot1*x13*x18;
    x20=x10 > x14;
    x22=mu*rhot2*x13*x18;
    x25=mu*x13*x18*x4;
    x27=mu*x13*x18*x7;
  };
  if (x1)
  {
    result[0]=0;
    result[9]=1;
  };
  if (x2)
  {
    result[0]=rhon;
    result[9]=0;
  };
  if ((x11 || x12) && (x15))
  {
    result[4]=rhot1;
    result[7]=0;
    result[10]=0;
    result[13]=0;
    result[16]=0;
    result[5]=0;
    result[8]=rhot2;
    result[11]=0;
    result[14]=0;
    result[17]=0;
  };
  if ((x11 || x12) && (x20))
  {
    result[4]=rhot1*x17 - x19*x5;
    result[7]=-x21*x22;
    result[10]=-x23*x4;
    result[13]=-x17 - x24*x25 + 1;
    result[16]=-x25*x26;
    result[5]=-x19*x21;
    result[8]=rhot2*x17 - x22*x8;
    result[11]=-x23*x7;
    result[14]=-x24*x27;
    result[17]=-x17 - x26*x27 + 1;
  };
  result[3]=0;
  result[6]=0;
  result[12]=0;
  result[15]=0;
  result[1]=0;
  result[2]=0;
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
  double x0=-rhon*un + rn;
  int    x1=x0 <= 0;
  int    x2=x0 > 0;
  double x3;
  double x4;
  double x5;
  double x6=rhot1*ut1;
  double x7=rt1 - x6;
  double x8=x7*x7;
  double x9=rhot2*ut2;
  double x10=rt2 - x9;
  double x11=x10*x10;
  double x12=x11 + x8;
  double x13=sqrt(x12);
  double x14;
  int x15 = 0;
  double x16=1./(x13);
  double x17;
  int x18 = 0;
  double x19;
  double x20=sqrt(1./(x12*x12*x12));
  double x21;
  double x22=x10*x7;
  double x23;
  double x24;
  double x25=-rt1 + x6;
  double x26;
  double x27=-rt2 + x9;
  double x28;
  if (x1)
  {
    x3=0;
    x4=0;
    x5=0;
  };
  if (x2)
  {
    x3=x0;
    x4=-rhon;
    x5=1;
  };
  if (x1 || x2)
  {
    x14=mu*x3;
    x15=x13 <= x14;
    x17=x14*x16;
    x18=x13 > x14;
    x19=mu*x16*x4;
    x21=mu*rhot1*x20*x3;
    x23=mu*rhot2*x20*x3;
    x24=mu*x16*x5;
    x26=mu*x20*x3*x7;
    x28=mu*x10*x20*x3;
    result[0]=rn - x3;
    result[3]=-x4;
    result[12]=-x5 + 1;
  };
  if ((x1 || x2) && (x15))
  {
    result[1]=rt1 - x7;
    result[4]=0;
    result[7]=rhot1;
    result[10]=0;
    result[13]=0;
    result[16]=0;
    result[19]=0;
    result[2]=rt2 - x10;
    result[5]=0;
    result[8]=0;
    result[11]=rhot2;
    result[14]=0;
    result[17]=0;
    result[20]=0;
  };
  if ((x1 || x2) && (x18))
  {
    result[1]=rt1 - x17*x7;
    result[4]=-x19*x7;
    result[7]=rhot1*x17 - x21*x8;
    result[10]=-x22*x23;
    result[13]=-x24*x7;
    result[16]=-x17 - x25*x26 + 1;
    result[19]=-x26*x27;
    result[2]=rt2 - x10*x17;
    result[5]=-x10*x19;
    result[8]=-x21*x22;
    result[11]=rhot2*x17 - x11*x23;
    result[14]=-x10*x24;
    result[17]=-x25*x28;
    result[20]=-x17 - x27*x28 + 1;
  };
  result[6]=0;
  result[9]=0;
  result[15]=0;
  result[18]=0;
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
  double x0=-rhon*un + rn;
  double x1;
  double x2=-rhot1*ut1 + rt1;
  double x3=-rhot2*ut2 + rt2;
  double x4=sqrt(pow(x2, 2) + pow(x3, 2));
  double x5;
  int x6 = 0;
  double x7;
  int x8 = 0;
  if (x0 <= 0)
  {
    x1=0;
  };
  if (x0 > 0)
  {
    x1=x0;
  };
  if (x0 <= 0 || x0 > 0)
  {
    x5=mu*x1;
    x6=x4 <= x5;
    x7=mu*x1/x4;
    x8=x4 > x5;
    result[0]=rn - x1;
  };
  if ((x0 <= 0 || x0 > 0) && (x6))
  {
    result[1]=rt1 - x2;
    result[2]=rt2 - x3;
  };
  if ((x0 <= 0 || x0 > 0) && (x8))
  {
    result[1]=rt1 - x2*x7;
    result[2]=rt2 - x3*x7;
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
  double x0=-rhon*un + rn;
  int    x1=x0 <= 0;
  int    x2=x0 > 0;
  double x3;
  double x4;
  double x5=rhot1*ut1;
  double x6=rt1 - x5;
  double x7=x6*x6;
  double x8=rhot2*ut2;
  double x9=rt2 - x8;
  double x10=x9*x9;
  double x11=x10 + x7;
  double x12=sqrt(x11);
  double x13;
  double x14;
  int x15 = 0;
  double x16=1./(x12);
  double x17;
  int x18 = 0;
  double x19;
  double x20=sqrt(1./(x11*x11*x11));
  double x21;
  double x22=x6*x9;
  double x23;
  double x24;
  double x25=-rt1 + x5;
  double x26;
  double x27=-rt2 + x8;
  double x28;
  if (x1)
  {
    x3=0;
    x4=0;
    x13=0;
  };
  if (x2)
  {
    x3=-rhon;
    x4=1;
    x13=x0;
  };
  if (x1 || x2)
  {
    x14=mu*x13;
    x15=x12 <= x14;
    x17=mu*x16*x3;
    x18=x12 > x14;
    x19=x14*x16;
    x21=mu*rhot1*x13*x20;
    x23=mu*rhot2*x13*x20;
    x24=mu*x16*x4;
    x26=mu*x13*x20*x6;
    x28=mu*x13*x20*x9;
    result[0]=-x3;
    result[9]=-x4 + 1;
  };
  if ((x1 || x2) && (x15))
  {
    result[1]=0;
    result[4]=rhot1;
    result[7]=0;
    result[10]=0;
    result[13]=0;
    result[16]=0;
    result[2]=0;
    result[5]=0;
    result[8]=rhot2;
    result[11]=0;
    result[14]=0;
    result[17]=0;
  };
  if ((x1 || x2) && (x18))
  {
    result[1]=-x17*x6;
    result[4]=rhot1*x19 - x21*x7;
    result[7]=-x22*x23;
    result[10]=-x24*x6;
    result[13]=-x19 - x25*x26 + 1;
    result[16]=-x26*x27;
    result[2]=-x17*x9;
    result[5]=-x21*x22;
    result[8]=rhot2*x19 - x10*x23;
    result[11]=-x24*x9;
    result[14]=-x25*x28;
    result[17]=-x19 - x27*x28 + 1;
  };
  result[3]=0;
  result[6]=0;
  result[12]=0;
  result[15]=0;
}


void frictionContact3D_AlartCurnierFunctionGenerated(
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

void frictionContact3D_localAlartCurnierJeanMoreauFunctionGenerated(
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

    frictionContact3D_AlartCurnierJeanMoreauFABGenerated(
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
      frictionContact3D_AlartCurnierJeanMoreauFGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3(result, f);
    }

    if (A && B)
    {
      frictionContact3D_AlartCurnierJeanMoreauABGenerated(
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
