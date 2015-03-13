#include <math.h>
#include <assert.h>
#include <op3x3.h>
#include <stdlib.h>

#define Sign(x) ((x>0) - (x<0))
#define Max fmax
#define Heaviside(x) (.5*Sign(x) + .5)
#define Rand(x) ((double) rand()/ (double) RAND_MAX)

#define random1 .5
#define random2 .5

#ifdef __cplusplus
#include <cmath>
#define CHECK(x)
#define XCHECK(x) assert(isfinite(x))
#else
#define CHECK(x)
#define XCHECK(x) assert(isfinite(x))
#endif

// temporary bug fix for overloaded pow. Sympy generates code with long double
// and it is not clear for all compiler which cast should be applied.
// The real fix is to prevent sympy from adding the type specifier
#ifdef __cplusplus
#define pow(x, y) std::pow(static_cast<double>(x), static_cast<double>(y))
#else
#define pow(x,y) pow((double)x, (double)y)
#endif

#pragma GCC diagnostic ignored "-Wconversion"

// hack, should be prevented in sage/sympy/maple or in code generation
#define sqrt(x) (fabs(x) < 1e-12 ? 0 : sqrt(x))

void frictionContact3D_FischerBurmeisterFABGenerated(
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
{
  double x0=pow(ut1, 2); CHECK(x0);
  double x1=pow(ut2, 2); CHECK(x1);
  double x2=mu*sqrt(x0 + x1) + un; CHECK(x2);
  double x3=mu*rn; CHECK(x3);
  double x4=pow(mu, 2); CHECK(x4);
  double x5=pow(rn, 2)*x4 + pow(rt1, 2) + pow(rt2, 2) + x0*x4 + x1*x4 + pow(x2, 2); CHECK(x5);
  double x6=mu*x2; CHECK(x6);
  double x7=2*sqrt(pow(rt1*x3 + ut1*x6, 2) + pow(rt2*x3 + ut2*x6, 2)); CHECK(x7);
  result[0]=x2 + x3 - 0.5*sqrt(x5 - x7) - 0.5*sqrt(x5 + x7); XCHECK(result[0]);
}
{
  double x0=pow(mu, 2); CHECK(x0);
  double x1=pow(ut1, 2); CHECK(x1);
  double x2=x0*x1; CHECK(x2);
  double x3=pow(ut2, 2); CHECK(x3);
  double x4=x0*x3; CHECK(x4);
  double x5=sqrt(x1 + x3); CHECK(x5);
  double x6=mu*x5 + un; CHECK(x6);
  double x7=pow(x6, 2); CHECK(x7);
  double x8=sqrt(x2 + x4 + x7); CHECK(x8);
  double x9=2.828427124*mu; CHECK(x9);
  double x10=pow(mu, 2); CHECK(x10);
  double x11=7071067810.0*mu; CHECK(x11);
  double x12=sqrt(2500000000.0*x10 - x11 + 7500000000.0); CHECK(x12);
  double x13=sqrt(7500000000.0*x10 + x11 + 2500000000.0); CHECK(x13);
  int    x14=x8 > 0;
  double x15=mu*rn; CHECK(x15);
  double x16=mu*x6; CHECK(x16);
  double x17=rt1*x15 + ut1*x16; CHECK(x17);
  double x18=rt2*x15 + ut2*x16; CHECK(x18);
  double x19=sqrt(pow(x17, 2) + pow(x18, 2)); CHECK(x19);
  int    x20=x19 <= 0;
  int    x21=x19 > 0;
  double x22=pow(rn, 2)*x0 + pow(rt1, 2) + pow(rt2, 2) + x2 + x4 + x7; CHECK(x22);
  double x23=2*x19; CHECK(x23);
  double x24=x22 - x23; CHECK(x24);
  int    x25=x24 <= 0;
  double x26=x22 + x23; CHECK(x26);
  int    x27=x26 <= 0;
  int    x28=x24 > 0;
  double x29=pow(rt1, 2); CHECK(x29);
  double x30=pow(rt2, 2); CHECK(x30);
  double x31=pow(rn, 2)*x10; CHECK(x31);
  double x32=pow(un, 2) + x29 + x30 + x31; CHECK(x32);
  double x33=2.0*sqrt(x29*x31 + x30*x31); CHECK(x33);
  double x34=sqrt(x32 - x33); CHECK(x34);
  double x35=sqrt(x32 + x33); CHECK(x35);
  double x36=0.5*un; CHECK(x36);
  int    x37=x26 > 0;
  double x38=(mu*ut1*x17 + mu*ut2*x18)/x19; CHECK(x38);
  if (x8 <= 0)
  {
    result[3]=1.00000000000000; XCHECK(result[3]);
  };
  if (x14 && (x14 || x20) && (x14 || x21) && (x14 || x25) && (x14 || x27) && (x14 || x28) && (x20 || x21) && (x14 || x20 || x21) && (x14 || x20 || x25) && (x14 || x20 || x27) && (x14 || x20 || x28) && (x14 || x21 || x25) && (x14 || x21 || x27) && (x14 || x21 || x28) && (x14 || x25 || x27) && (x14 || x25 || x28) && (x20 || x21 || x25) && (x20 || x21 || x27) && (x20 || x21 || x28) && (x20 || x25 || x27) && (x20 || x25 || x28))
  {
    result[3]=(-1.6329931624e-5*mu*x12 + 4.0e-10*x12*x13 - 5.773502694e-6*x12 - 1.0e-5*x13)/(sqrt(x10 - x9 + 3.0)*sqrt(3.0*x10 + x9 + 1.0)); XCHECK(result[3]);
  };
  if (x14 && x21 && x28 && x37 && x5 <= 0)
  {
    result[3]=(x34*x35 - x34*x36 - x35*x36)/(x34*x35); XCHECK(result[3]);
  };
  if (x14 && x21 && x28 && x37 && x5 > 0)
  {
    result[3]=1 - 0.5*(x38 + x6)/sqrt(x26) - 0.5*(-x38 + x6)/sqrt(x24); XCHECK(result[3]);
  };
}
{
  double x0=0.707106781*mu; CHECK(x0);
  double x1=pow(mu, 2); CHECK(x1);
  double x2=pow(ut1, 2); CHECK(x2);
  double x3=x1*x2; CHECK(x3);
  double x4=pow(ut2, 2); CHECK(x4);
  double x5=x1*x4; CHECK(x5);
  double x6=sqrt(x2 + x4); CHECK(x6);
  double x7=mu*x6 + un; CHECK(x7);
  double x8=pow(x7, 2); CHECK(x8);
  double x9=sqrt(x3 + x5 + x8); CHECK(x9);
  double x10=pow(mu, 2); CHECK(x10);
  int    x11=x9 > 0;
  double x12=mu*rn; CHECK(x12);
  double x13=mu*ut1; CHECK(x13);
  double x14=rt1*x12 + x13*x7; CHECK(x14);
  double x15=mu*x7; CHECK(x15);
  double x16=rt2*x12 + ut2*x15; CHECK(x16);
  double x17=sqrt(pow(x14, 2) + pow(x16, 2)); CHECK(x17);
  int    x18=x17 <= 0;
  int    x19=x17 > 0;
  double x20=pow(rn, 2)*x1 + pow(rt1, 2) + pow(rt2, 2) + x3 + x5 + x8; CHECK(x20);
  double x21=2*x17; CHECK(x21);
  double x22=x20 - x21; CHECK(x22);
  int    x23=x22 <= 0;
  double x24=x20 + x21; CHECK(x24);
  int    x25=x24 <= 0;
  int    x26=x22 > 0;
  double x27=pow(rt1, 2); CHECK(x27);
  double x28=pow(rn, 2)*x10; CHECK(x28);
  double x29=pow(rt2, 2); CHECK(x29);
  double x30=sqrt(x27*x28 + x28*x29); CHECK(x30);
  double x31=pow(un, 2) + x27 + x28 + x29; CHECK(x31);
  double x32=2.0*x30; CHECK(x32);
  double x33=sqrt(x31 - x32); CHECK(x33);
  double x34=sqrt(x31 + x32); CHECK(x34);
  double x35=0.5*rn*rt1*un*x10; CHECK(x35);
  double x36=0.3535533905*mu*un*x30; CHECK(x36);
  int    x37=x24 > 0;
  double x38=1.0/x6; CHECK(x38);
  double x39=x13*x38; CHECK(x39);
  double x40=ut1*x1; CHECK(x40);
  double x41=x39*x7 + x40; CHECK(x41);
  double x42=(ut2*x16*x38*x40 + (1.0/2.0)*x14*(2*x15 + 2*x3*x38))/x17; CHECK(x42);
  if (x9 <= 0)
  {
    result[6]=x0; XCHECK(result[6]);
  };
  if (x11 && (x11 || x18) && (x11 || x19) && (x11 || x23) && (x11 || x25) && (x11 || x26) && (x18 || x19) && (x11 || x18 || x19) && (x11 || x18 || x23) && (x11 || x18 || x25) && (x11 || x18 || x26) && (x11 || x19 || x23) && (x11 || x19 || x25) && (x11 || x19 || x26) && (x11 || x23 || x25) && (x11 || x23 || x26) && (x18 || x19 || x23) && (x18 || x19 || x25) && (x18 || x19 || x26) && (x18 || x23 || x25) && (x18 || x23 || x26))
  {
    result[6]=(1.414213562e-5*mu*sqrt(7071067810.0*mu + 7500000000.0*x10 + 2500000000.0) - 0.4082482906*mu - 1.154700539*x10)/sqrt(2.828427124*mu + 3.0*x10 + 1.0); XCHECK(result[6]);
  };
  if (x11 && x19 && x26 && x37 && x6 <= 0)
  {
    result[6]=(x0*x30*x33*x34 - x33*x35 - x33*x36 + x34*x35 - x34*x36)/(x30*x33*x34); XCHECK(result[6]);
  };
  if (x11 && x19 && x26 && x37 && x6 > 0)
  {
    result[6]=x39 - 0.5*(x41 + x42)/sqrt(x24) - 0.5*(x41 - x42)/sqrt(x22); XCHECK(result[6]);
  };
}
{
  double x0=0.707106781*mu; CHECK(x0);
  double x1=pow(mu, 2); CHECK(x1);
  double x2=pow(ut1, 2); CHECK(x2);
  double x3=x1*x2; CHECK(x3);
  double x4=pow(ut2, 2); CHECK(x4);
  double x5=x1*x4; CHECK(x5);
  double x6=sqrt(x2 + x4); CHECK(x6);
  double x7=mu*x6 + un; CHECK(x7);
  double x8=pow(x7, 2); CHECK(x8);
  double x9=sqrt(x3 + x5 + x8); CHECK(x9);
  double x10=pow(mu, 2); CHECK(x10);
  int    x11=x9 > 0;
  double x12=mu*rn; CHECK(x12);
  double x13=mu*x7; CHECK(x13);
  double x14=rt1*x12 + ut1*x13; CHECK(x14);
  double x15=mu*ut2; CHECK(x15);
  double x16=rt2*x12 + x15*x7; CHECK(x16);
  double x17=sqrt(pow(x14, 2) + pow(x16, 2)); CHECK(x17);
  int    x18=x17 <= 0;
  int    x19=x17 > 0;
  double x20=pow(rn, 2)*x1 + pow(rt1, 2) + pow(rt2, 2) + x3 + x5 + x8; CHECK(x20);
  double x21=2*x17; CHECK(x21);
  double x22=x20 - x21; CHECK(x22);
  int    x23=x22 <= 0;
  double x24=x20 + x21; CHECK(x24);
  int    x25=x24 <= 0;
  int    x26=x22 > 0;
  double x27=pow(rt1, 2); CHECK(x27);
  double x28=pow(rn, 2)*x10; CHECK(x28);
  double x29=pow(rt2, 2); CHECK(x29);
  double x30=sqrt(x27*x28 + x28*x29); CHECK(x30);
  double x31=pow(un, 2) + x27 + x28 + x29; CHECK(x31);
  double x32=2.0*x30; CHECK(x32);
  double x33=sqrt(x31 - x32); CHECK(x33);
  double x34=sqrt(x31 + x32); CHECK(x34);
  double x35=0.5*rn*rt2*un*x10; CHECK(x35);
  double x36=0.3535533905*mu*un*x30; CHECK(x36);
  int    x37=x24 > 0;
  double x38=1.0/x6; CHECK(x38);
  double x39=x15*x38; CHECK(x39);
  double x40=ut2*x1; CHECK(x40);
  double x41=x39*x7 + x40; CHECK(x41);
  double x42=(ut1*x14*x38*x40 + (1.0/2.0)*x16*(2*x13 + 2*x38*x5))/x17; CHECK(x42);
  if (x9 <= 0)
  {
    result[9]=x0; XCHECK(result[9]);
  };
  if (x11 && (x11 || x18) && (x11 || x19) && (x11 || x23) && (x11 || x25) && (x11 || x26) && (x18 || x19) && (x11 || x18 || x19) && (x11 || x18 || x23) && (x11 || x18 || x25) && (x11 || x18 || x26) && (x11 || x19 || x23) && (x11 || x19 || x25) && (x11 || x19 || x26) && (x11 || x23 || x25) && (x11 || x23 || x26) && (x18 || x19 || x23) && (x18 || x19 || x25) && (x18 || x19 || x26) && (x18 || x23 || x25) && (x18 || x23 || x26))
  {
    result[9]=(1.414213562e-5*mu*sqrt(7071067810.0*mu + 7500000000.0*x10 + 2500000000.0) - 0.4082482906*mu - 1.154700539*x10)/sqrt(2.828427124*mu + 3.0*x10 + 1.0); XCHECK(result[9]);
  };
  if (x11 && x19 && x26 && x37 && x6 <= 0)
  {
    result[9]=(x0*x30*x33*x34 - x33*x35 - x33*x36 + x34*x35 - x34*x36)/(x30*x33*x34); XCHECK(result[9]);
  };
  if (x11 && x19 && x26 && x37 && x6 > 0)
  {
    result[9]=x39 - 0.5*(x41 + x42)/sqrt(x24) - 0.5*(x41 - x42)/sqrt(x22); XCHECK(result[9]);
  };
}
{
  double x0=pow(rt1, 2); CHECK(x0);
  double x1=pow(rt2, 2); CHECK(x1);
  double x2=pow(mu, 2); CHECK(x2);
  double x3=pow(rn, 2)*x2; CHECK(x3);
  double x4=sqrt(x0 + x1 + x3); CHECK(x4);
  double x5=2.828427124*mu; CHECK(x5);
  double x6=pow(mu, 2); CHECK(x6);
  double x7=7071067810.0*mu; CHECK(x7);
  double x8=sqrt(2500000000.0*x6 - x7 + 7500000000.0); CHECK(x8);
  double x9=mu*x8; CHECK(x9);
  double x10=sqrt(7500000000.0*x6 + x7 + 2500000000.0); CHECK(x10);
  double x11=mu*rn; CHECK(x11);
  double x12=pow(ut1, 2); CHECK(x12);
  double x13=pow(ut2, 2); CHECK(x13);
  double x14=mu*sqrt(x12 + x13) + un; CHECK(x14);
  double x15=mu*x14; CHECK(x15);
  double x16=rt1*x11 + ut1*x15; CHECK(x16);
  double x17=rt2*x11 + ut2*x15; CHECK(x17);
  double x18=sqrt(pow(x16, 2) + pow(x17, 2)); CHECK(x18);
  double x19=x0 + x1 + x12*x2 + x13*x2 + pow(x14, 2) + x3; CHECK(x19);
  double x20=2*x18; CHECK(x20);
  double x21=x19 - x20; CHECK(x21);
  double x22=x19 + x20; CHECK(x22);
  double x23=rn*x2; CHECK(x23);
  double x24=(mu*rt1*x16 + mu*rt2*x17)/x18; CHECK(x24);
  if (x4 <= 0)
  {
    result[12]=mu; XCHECK(result[12]);
  };
  if (x18 <= 0 || x21 <= 0 || x22 <= 0)
  {
    result[12]=(1.414213562e-5*mu*x10 - 1.0e-5*x10*x6 + 4.0e-10*x10*x9 - 5.773502694e-6*x6*x8 - 8.164965812e-6*x9)/(sqrt(-x5 + x6 + 3.0)*sqrt(x5 + 3.0*x6 + 1.0)); XCHECK(result[12]);
  };
  if (x18 > 0 && x21 > 0 && x22 > 0 && x4 > 0)
  {
    result[12]=mu - 0.5*(x23 + x24)/sqrt(x22) - 0.5*(x23 - x24)/sqrt(x21); XCHECK(result[12]);
  };
}
{
  double x0=pow(rt1, 2); CHECK(x0);
  double x1=pow(rt2, 2); CHECK(x1);
  double x2=pow(mu, 2); CHECK(x2);
  double x3=pow(rn, 2)*x2; CHECK(x3);
  double x4=sqrt(x0 + x1 + x3); CHECK(x4);
  double x5=2.828427124*mu; CHECK(x5);
  double x6=pow(mu, 2); CHECK(x6);
  double x7=7071067810.0*mu; CHECK(x7);
  double x8=sqrt(2500000000.0*x6 - x7 + 7500000000.0); CHECK(x8);
  double x9=sqrt(7500000000.0*x6 + x7 + 2500000000.0); CHECK(x9);
  double x10=mu*rn; CHECK(x10);
  double x11=pow(ut1, 2); CHECK(x11);
  double x12=pow(ut2, 2); CHECK(x12);
  double x13=mu*sqrt(x11 + x12) + un; CHECK(x13);
  double x14=mu*x13; CHECK(x14);
  double x15=rt1*x10 + ut1*x14; CHECK(x15);
  double x16=sqrt(pow(x15, 2) + pow(rt2*x10 + ut2*x14, 2)); CHECK(x16);
  double x17=x0 + x1 + x11*x2 + x12*x2 + pow(x13, 2) + x3; CHECK(x17);
  double x18=2*x16; CHECK(x18);
  double x19=x17 - x18; CHECK(x19);
  double x20=x17 + x18; CHECK(x20);
  double x21=x10*x15/x16; CHECK(x21);
  if (x4 <= 0)
  {
    result[15]=0.0; XCHECK(result[15]);
  };
  if (x16 <= 0 || x19 <= 0 || x20 <= 0)
  {
    result[15]=(-4.082482906e-6*mu*x8 + 7.07106781e-6*mu*x9 - 5.773502694e-6*x8 - 1.0e-5*x9)/(sqrt(-x5 + x6 + 3.0)*sqrt(x5 + 3.0*x6 + 1.0)); XCHECK(result[15]);
  };
  if (x16 > 0 && x19 > 0 && x20 > 0 && x4 > 0)
  {
    result[15]=-0.5*(rt1 + x21)/sqrt(x20) - 0.5*(rt1 - x21)/sqrt(x19); XCHECK(result[15]);
  };
}
{
  double x0=pow(rt1, 2); CHECK(x0);
  double x1=pow(rt2, 2); CHECK(x1);
  double x2=pow(mu, 2); CHECK(x2);
  double x3=pow(rn, 2)*x2; CHECK(x3);
  double x4=sqrt(x0 + x1 + x3); CHECK(x4);
  double x5=2.828427124*mu; CHECK(x5);
  double x6=pow(mu, 2); CHECK(x6);
  double x7=7071067810.0*mu; CHECK(x7);
  double x8=sqrt(2500000000.0*x6 - x7 + 7500000000.0); CHECK(x8);
  double x9=sqrt(7500000000.0*x6 + x7 + 2500000000.0); CHECK(x9);
  double x10=mu*rn; CHECK(x10);
  double x11=pow(ut1, 2); CHECK(x11);
  double x12=pow(ut2, 2); CHECK(x12);
  double x13=mu*sqrt(x11 + x12) + un; CHECK(x13);
  double x14=mu*x13; CHECK(x14);
  double x15=rt2*x10 + ut2*x14; CHECK(x15);
  double x16=sqrt(pow(x15, 2) + pow(rt1*x10 + ut1*x14, 2)); CHECK(x16);
  double x17=x0 + x1 + x11*x2 + x12*x2 + pow(x13, 2) + x3; CHECK(x17);
  double x18=2*x16; CHECK(x18);
  double x19=x17 - x18; CHECK(x19);
  double x20=x17 + x18; CHECK(x20);
  double x21=x10*x15/x16; CHECK(x21);
  if (x4 <= 0)
  {
    result[18]=0.0; XCHECK(result[18]);
  };
  if (x16 <= 0 || x19 <= 0 || x20 <= 0)
  {
    result[18]=(-4.082482906e-6*mu*x8 + 7.07106781e-6*mu*x9 - 5.773502694e-6*x8 - 1.0e-5*x9)/(sqrt(-x5 + x6 + 3.0)*sqrt(x5 + 3.0*x6 + 1.0)); XCHECK(result[18]);
  };
  if (x16 > 0 && x19 > 0 && x20 > 0 && x4 > 0)
  {
    result[18]=-0.5*(rt2 + x21)/sqrt(x20) - 0.5*(rt2 - x21)/sqrt(x19); XCHECK(result[18]);
  };
}
{
  double x0=mu*ut1; CHECK(x0);
  double x1=pow(mu, 2); CHECK(x1);
  double x2=pow(ut1, 2); CHECK(x2);
  double x3=pow(ut2, 2); CHECK(x3);
  double x4=mu*sqrt(x2 + x3) + un; CHECK(x4);
  double x5=pow(rn, 2)*x1 + pow(rt1, 2) + pow(rt2, 2) + x1*x2 + x1*x3 + pow(x4, 2); CHECK(x5);
  double x6=mu*rn; CHECK(x6);
  double x7=rt1*x6 + x0*x4; CHECK(x7);
  double x8=sqrt(pow(x7, 2) + pow(mu*ut2*x4 + rt2*x6, 2)); CHECK(x8);
  double x9=2*x8; CHECK(x9);
  double x10 = NAN;
  if (x8 > 0)
  {
    x10=0.5*x7/x8; XCHECK(x10);
  };
  if (x8 <= 0)
  {
    x10=0.5*random1/sqrt(pow(random1, 2) + pow(random2, 2)); XCHECK(x10);
  };
  result[1]=rt1 + x0 + x10*sqrt(x5 - x9) - x10*sqrt(x5 + x9); XCHECK(result[1]);
}
{
  double x0=mu*rn; CHECK(x0);
  double x1=pow(ut1, 2); CHECK(x1);
  double x2=pow(ut2, 2); CHECK(x2);
  double x3=sqrt(x1 + x2); CHECK(x3);
  double x4=mu*x3 + un; CHECK(x4);
  double x5=mu*x4; CHECK(x5);
  double x6=rt1*x0 + ut1*x5; CHECK(x6);
  double x7=rt2*x0 + ut2*x5; CHECK(x7);
  double x8=pow(x6, 2) + pow(x7, 2); CHECK(x8);
  double x9=sqrt(x8); CHECK(x9);
  int    x10=x9 <= 0;
  int    x11=x9 > 0;
  double x12=x10 || x11; CHECK(x12);
  double x13=pow(mu, 2); CHECK(x13);
  double x14=x1*x13; CHECK(x14);
  double x15=x13*x2; CHECK(x15);
  double x16=pow(x4, 2); CHECK(x16);
  double x17=sqrt(x14 + x15 + x16); CHECK(x17);
  double x18=7071067810.0*mu; CHECK(x18);
  double x19=pow(mu, 2); CHECK(x19);
  double x20=sqrt(-x18 + 2500000000.0*x19 + 7500000000.0); CHECK(x20);
  double x21=sqrt(x18 + 7500000000.0*x19 + 2500000000.0); CHECK(x21);
  double x22=mu*x20; CHECK(x22);
  double x23=pow(pow(random1, 2) + pow(random2, 2), -1.0/2.0); CHECK(x23);
  double x24=2.828427124*mu; CHECK(x24);
  int    x25=x17 > 0;
  double x26=pow(rn, 2)*x13 + pow(rt1, 2) + pow(rt2, 2) + x14 + x15 + x16; CHECK(x26);
  double x27=2*x9; CHECK(x27);
  double x28=x26 - x27; CHECK(x28);
  int    x29=x28 <= 0;
  double x30=x26 + x27; CHECK(x30);
  int    x31=x30 <= 0;
  int    x32=x28 > 0;
  double x33=pow(rt1, 2); CHECK(x33);
  double x34=pow(rt2, 2); CHECK(x34);
  double x35=pow(rn, 2)*x19; CHECK(x35);
  double x36=pow(un, 2) + x33 + x34 + x35; CHECK(x36);
  double x37=sqrt(x33*x35 + x34*x35); CHECK(x37);
  double x38=2.0*x37; CHECK(x38);
  double x39=sqrt(x36 - x38); CHECK(x39);
  double x40=sqrt(x36 + x38); CHECK(x40);
  double x41=1/(x39*x40); CHECK(x41);
  double x42=0.5*mu*rn*rt1*un; CHECK(x42);
  double x43=0.5*random1*un; CHECK(x43);
  int    x44=x30 > 0;
  double x45=sqrt(x28); CHECK(x45);
  double x46=mu*ut1; CHECK(x46);
  double x47=1.0/x9; CHECK(x47);
  double x48=x46*x6; CHECK(x48);
  double x49=mu*ut2*x7; CHECK(x49);
  double x50 = NAN;
  double x51=sqrt(x30); CHECK(x51);
  double x52 = NAN;
  double x53=x47*(x48 + x49); CHECK(x53);
  if (x11)
  {
    x50=0.5*x46*x47 + 0.5*x6*(-x48 - x49)/pow(x8, 3.0/2.0); XCHECK(x50);
    x52=0.5*x47*x6; XCHECK(x52);
  };
  if (x10)
  {
    x50=0; XCHECK(x50);
    x52=0.5*random1*x23; XCHECK(x52);
  };
  if ((x17 <= 0) && (x12))
  {
    result[4]=0.0; XCHECK(result[4]);
  };
  if ((x12 && x25 && (x10 || x25) && (x11 || x25) && (x25 || x29) && (x25 || x31) && (x25 || x32) && (x10 || x11 || x25) && (x10 || x11 || x29) && (x10 || x11 || x31) && (x10 || x11 || x32) && (x10 || x25 || x29) && (x10 || x25 || x31) && (x10 || x25 || x32) && (x10 || x29 || x31) && (x10 || x29 || x32) && (x11 || x25 || x29) && (x11 || x25 || x31) && (x11 || x25 || x32) && (x25 || x29 || x31) && (x25 || x29 || x32)) && (x11))
  {
    result[4]=(-10206.20726*x20 + 17677.66952*x21 - 28867.51346*x22)/(x20*x21); XCHECK(result[4]);
  };
  if ((x12 && x25 && (x10 || x25) && (x11 || x25) && (x25 || x29) && (x25 || x31) && (x25 || x32) && (x10 || x11 || x25) && (x10 || x11 || x29) && (x10 || x11 || x31) && (x10 || x11 || x32) && (x10 || x25 || x29) && (x10 || x25 || x31) && (x10 || x25 || x32) && (x10 || x29 || x31) && (x10 || x29 || x32) && (x11 || x25 || x29) && (x11 || x25 || x31) && (x11 || x25 || x32) && (x25 || x29 || x31) && (x25 || x29 || x32)) && (x10))
  {
    result[4]=x23*(-5.773502694e-6*random1*x20 + 1.0e-5*random1*x21 - 1.6329931624e-5*random1*x22)/(sqrt(x19 - x24 + 3.0)*sqrt(3.0*x19 + x24 + 1.0)); XCHECK(result[4]);
  };
  if ((x11 && x25 && x32 && x44 && x3 <= 0) && (x11))
  {
    result[4]=x41*(-x39*x42 + x40*x42)/x37; XCHECK(result[4]);
  };
  if ((x11 && x25 && x32 && x44 && x3 <= 0) && (x10))
  {
    result[4]=x23*x41*(-x39*x43 + x40*x43); XCHECK(result[4]);
  };
  if (x11 && x25 && x32 && x44 && x3 > 0)
  {
    result[4]=x45*x50 - x50*x51 - x52*(x4 + x53)/x51 + x52*(x4 - x53)/x45; XCHECK(result[4]);
  };
}
{
  double x0=mu*rn; CHECK(x0);
  double x1=pow(ut1, 2); CHECK(x1);
  double x2=pow(ut2, 2); CHECK(x2);
  double x3=sqrt(x1 + x2); CHECK(x3);
  double x4=mu*x3 + un; CHECK(x4);
  double x5=mu*x4; CHECK(x5);
  double x6=ut1*x5; CHECK(x6);
  double x7=rt1*x0 + x6; CHECK(x7);
  double x8=rt2*x0 + ut2*x5; CHECK(x8);
  double x9=pow(x7, 2) + pow(x8, 2); CHECK(x9);
  double x10=sqrt(x9); CHECK(x10);
  int    x11=x10 <= 0;
  int    x12=x10 > 0;
  double x13=x11 || x12; CHECK(x13);
  double x14=pow(mu, 2); CHECK(x14);
  double x15=x1*x14; CHECK(x15);
  double x16=x14*x2; CHECK(x16);
  double x17=pow(x4, 2); CHECK(x17);
  double x18=sqrt(x15 + x16 + x17); CHECK(x18);
  double x19=2.828427124*mu; CHECK(x19);
  double x20=pow(mu, 2); CHECK(x20);
  double x21=sqrt(x19 + 3.0*x20 + 1.0); CHECK(x21);
  double x22=x21*sqrt(-x19 + x20 + 3.0); CHECK(x22);
  double x23=7071067810.0*mu; CHECK(x23);
  double x24=sqrt(2500000000.0*x20 - x23 + 7500000000.0); CHECK(x24);
  double x25=sqrt(7500000000.0*x20 + x23 + 2500000000.0); CHECK(x25);
  double x26=13258252140000.0*x25; CHECK(x26);
  double x27=mu*x24; CHECK(x27);
  double x28=6250000000000.0*x25; CHECK(x28);
  double x29=x20*x24; CHECK(x29);
  double x30=pow(mu, 3); CHECK(x30);
  double x31=sqrt(pow(random1, 2) + pow(random2, 2)); CHECK(x31);
  double x32=1.0/x31; CHECK(x32);
  double x33=mu*x31; CHECK(x33);
  int    x34=x18 > 0;
  double x35=pow(rn, 2)*x14 + pow(rt1, 2) + pow(rt2, 2) + x15 + x16 + x17; CHECK(x35);
  double x36=2*x10; CHECK(x36);
  double x37=x35 - x36; CHECK(x37);
  int    x38=x37 <= 0;
  double x39=x35 + x36; CHECK(x39);
  int    x40=x39 <= 0;
  int    x41=x37 > 0;
  double x42=pow(rt1, 2); CHECK(x42);
  double x43=pow(rn, 2); CHECK(x43);
  double x44=x20*x43; CHECK(x44);
  double x45=pow(rt2, 2); CHECK(x45);
  double x46=x42*x44 + x44*x45; CHECK(x46);
  double x47=pow(x46, 3.0/2.0); CHECK(x47);
  double x48=pow(un, 2) + x42 + x44 + x45; CHECK(x48);
  double x49=sqrt(x46); CHECK(x49);
  double x50=2.0*x49; CHECK(x50);
  double x51=sqrt(x48 - x50); CHECK(x51);
  double x52=sqrt(x48 + x50); CHECK(x52);
  double x53=4.0*x47*x51*x52; CHECK(x53);
  double x54=x42*x53; CHECK(x54);
  double x55=x45*x53; CHECK(x55);
  double x56=2.0*pow(rt2, 6)*un*x30*x43; CHECK(x56);
  double x57=pow(mu, 4); CHECK(x57);
  double x58=pow(rn, 3); CHECK(x58);
  double x59=1.414213562*pow(rt1, 5)*un*x57*x58; CHECK(x59);
  double x60=pow(mu, 5); CHECK(x60);
  double x61=pow(rn, 4); CHECK(x61);
  double x62=pow(rt2, 4); CHECK(x62);
  double x63=2.0*un*x60*x61*x62; CHECK(x63);
  double x64=1.414213562*rt1*un*x57*x58*x62; CHECK(x64);
  double x65=pow(un, 3); CHECK(x65);
  double x66=2.0*x30*x43*x62*x65; CHECK(x66);
  double x67=4.0*un*x30*x42*x43*x62; CHECK(x67);
  double x68=2.0*pow(rt1, 4)*un*x30*x43*x45; CHECK(x68);
  double x69=2.828427124*pow(rt1, 3)*un*x45*x57*x58; CHECK(x69);
  double x70=2.0*un*x42*x45*x60*x61; CHECK(x70);
  double x71=2.0*x30*x42*x43*x45*x65; CHECK(x71);
  double x72=2.0*mu*un*x42*x47; CHECK(x72);
  double x73=4.0*mu*un*x45*x47; CHECK(x73);
  double x74=0.5*random1*rn*rt1*un*x20; CHECK(x74);
  double x75=0.3535533905*mu*random1*un*x49; CHECK(x75);
  int    x76=x39 > 0;
  double x77=sqrt(x37); CHECK(x77);
  double x78=1.0/x3; CHECK(x78);
  double x79=x15*x78; CHECK(x79);
  double x80=1.0/x10; CHECK(x80);
  double x81=ut1*x14; CHECK(x81);
  double x82=ut2*x78*x8*x81; CHECK(x82);
  double x83=(1.0/2.0)*x7*(2*x5 + 2*x79); CHECK(x83);
  double x84 = NAN;
  double x85=sqrt(x39); CHECK(x85);
  double x86 = NAN;
  double x87=x6*x78 + x81; CHECK(x87);
  double x88=x80*(x82 + x83); CHECK(x88);
  if (x12)
  {
    x84=0.5*x7*(-x82 - x83)/pow(x9, 3.0/2.0) + 0.5*x80*(x5 + x79); XCHECK(x84);
    x86=0.5*x7*x80; XCHECK(x86);
  };
  if (x11)
  {
    x84=0; XCHECK(x84);
    x86=0.5*random1*x32; XCHECK(x86);
  };
  if ((x18 <= 0) && (x13))
  {
    result[7]=mu; XCHECK(result[7]);
  };
  if ((x13 && x34 && (x11 || x34) && (x12 || x34) && (x34 || x38) && (x34 || x40) && (x34 || x41) && (x11 || x12 || x34) && (x11 || x12 || x38) && (x11 || x12 || x40) && (x11 || x12 || x41) && (x11 || x34 || x38) && (x11 || x34 || x40) && (x11 || x34 || x41) && (x11 || x38 || x40) && (x11 || x38 || x41) && (x12 || x34 || x38) && (x12 || x34 || x40) && (x12 || x34 || x41) && (x34 || x38 || x40) && (x34 || x38 || x41)) && (x12))
  {
    result[7]=(mu*x28 - x20*x26 - 61343466120000.0*x24*x30 - 7654655448000.0*x24 + 1000000000.0*x25*x27 + 707106781.0*x25*x29 + x26 - 46909709380000.0*x27 + x28*x30 - 104613624400000.0*x29)/(1.767766952e+18*mu*x22 + 2.5e+18*x22); XCHECK(result[7]);
  };
  if ((x13 && x34 && (x11 || x34) && (x12 || x34) && (x34 || x38) && (x34 || x40) && (x34 || x41) && (x11 || x12 || x34) && (x11 || x12 || x38) && (x11 || x12 || x40) && (x11 || x12 || x41) && (x11 || x34 || x38) && (x11 || x34 || x40) && (x11 || x34 || x41) && (x11 || x38 || x40) && (x11 || x38 || x41) && (x12 || x34 || x38) && (x12 || x34 || x40) && (x12 || x34 || x41) && (x34 || x38 || x40) && (x34 || x38 || x41)) && (x11))
  {
    result[7]=x32*(-0.4082482906*mu*random1 - 1.154700539*random1*x20 + 2.0e-5*x25*x33)/x21; XCHECK(result[7]);
  };
  if ((x12 && x34 && x41 && x76 && x3 <= 0) && (x12))
  {
    result[7]=(mu*x54 + mu*x55 - x51*x56 - x51*x59 - x51*x63 - x51*x64 - x51*x66 - x51*x67 - x51*x68 - x51*x69 - x51*x70 - x51*x71 - x51*x72 - x51*x73 + x52*x56 + x52*x59 + x52*x63 + x52*x64 + x52*x66 + x52*x67 + x52*x68 + x52*x69 + x52*x70 + x52*x71 - x52*x72 - x52*x73)/(x54 + x55); XCHECK(result[7]);
  };
  if ((x12 && x34 && x41 && x76 && x3 <= 0) && (x11))
  {
    result[7]=x32*(x33*x49*x51*x52 - x51*x74 - x51*x75 - x52*x74 + x52*x75)/(x49*x51*x52); XCHECK(result[7]);
  };
  if (x12 && x34 && x41 && x76 && x3 > 0)
  {
    result[7]=mu + x77*x84 - x84*x85 - x86*(x87 + x88)/x85 + x86*(x87 - x88)/x77; XCHECK(result[7]);
  };
}
{
  double x0=mu*rn; CHECK(x0);
  double x1=pow(ut1, 2); CHECK(x1);
  double x2=pow(ut2, 2); CHECK(x2);
  double x3=sqrt(x1 + x2); CHECK(x3);
  double x4=mu*x3 + un; CHECK(x4);
  double x5=mu*x4; CHECK(x5);
  double x6=rt1*x0 + ut1*x5; CHECK(x6);
  double x7=ut2*x5; CHECK(x7);
  double x8=rt2*x0 + x7; CHECK(x8);
  double x9=pow(x6, 2) + pow(x8, 2); CHECK(x9);
  double x10=sqrt(x9); CHECK(x10);
  int    x11=x10 <= 0;
  int    x12=x10 > 0;
  double x13=x11 || x12; CHECK(x13);
  double x14=pow(mu, 2); CHECK(x14);
  double x15=x1*x14; CHECK(x15);
  double x16=x14*x2; CHECK(x16);
  double x17=pow(x4, 2); CHECK(x17);
  double x18=sqrt(x15 + x16 + x17); CHECK(x18);
  double x19=2.828427124*mu; CHECK(x19);
  double x20=pow(mu, 2); CHECK(x20);
  double x21=sqrt(x19 + 3.0*x20 + 1.0); CHECK(x21);
  double x22=x21*sqrt(-x19 + x20 + 3.0); CHECK(x22);
  double x23=7071067810.0*mu; CHECK(x23);
  double x24=sqrt(2500000000.0*x20 - x23 + 7500000000.0); CHECK(x24);
  double x25=sqrt(7500000000.0*x20 + x23 + 2500000000.0); CHECK(x25);
  double x26=5303300858.0*x25; CHECK(x26);
  double x27=2500000000.0*x25; CHECK(x27);
  double x28=pow(mu, 3); CHECK(x28);
  double x29=pow(pow(random1, 2) + pow(random2, 2), -1.0/2.0); CHECK(x29);
  int    x30=x18 > 0;
  double x31=pow(rn, 2)*x14 + pow(rt1, 2) + pow(rt2, 2) + x15 + x16 + x17; CHECK(x31);
  double x32=2*x10; CHECK(x32);
  double x33=x31 - x32; CHECK(x33);
  int    x34=x33 <= 0;
  double x35=x31 + x32; CHECK(x35);
  int    x36=x35 <= 0;
  int    x37=x33 > 0;
  double x38=pow(rt1, 2); CHECK(x38);
  double x39=pow(rn, 2); CHECK(x39);
  double x40=x20*x39; CHECK(x40);
  double x41=pow(rt2, 2); CHECK(x41);
  double x42=x38*x40 + x40*x41; CHECK(x42);
  double x43=pow(x42, 3.0/2.0); CHECK(x43);
  double x44=pow(un, 2) + x38 + x40 + x41; CHECK(x44);
  double x45=sqrt(x42); CHECK(x45);
  double x46=2.0*x45; CHECK(x46);
  double x47=sqrt(x44 - x46); CHECK(x47);
  double x48=sqrt(x44 + x46); CHECK(x48);
  double x49=4.0*x43*x47*x48; CHECK(x49);
  double x50=pow(mu, 4); CHECK(x50);
  double x51=pow(rn, 3); CHECK(x51);
  double x52=pow(rt1, 5); CHECK(x52);
  double x53=1.414213562*un*x50*x51*x52; CHECK(x53);
  double x54=2.0*rt1*pow(rt2, 5)*un*x28*x39; CHECK(x54);
  double x55=2.0*rt2*un*x28*x39*x52; CHECK(x55);
  double x56=1.414213562*rt1*pow(rt2, 4)*un*x50*x51; CHECK(x56);
  double x57=pow(mu, 5); CHECK(x57);
  double x58=pow(rn, 4); CHECK(x58);
  double x59=pow(rt2, 3); CHECK(x59);
  double x60=2.0*rt1*un*x57*x58*x59; CHECK(x60);
  double x61=pow(rt1, 3); CHECK(x61);
  double x62=2.0*rt2*un*x57*x58*x61; CHECK(x62);
  double x63=pow(un, 3); CHECK(x63);
  double x64=2.0*rt1*x28*x39*x59*x63; CHECK(x64);
  double x65=2.0*rt2*x28*x39*x61*x63; CHECK(x65);
  double x66=4.0*un*x28*x39*x59*x61; CHECK(x66);
  double x67=2.828427124*un*x41*x50*x51*x61; CHECK(x67);
  double x68=2.0*mu*rt1*rt2*un*x43; CHECK(x68);
  double x69=0.5*random1*rn*rt2*un*x20; CHECK(x69);
  double x70=0.3535533905*mu*random1*un*x45; CHECK(x70);
  int    x71=x35 > 0;
  double x72=sqrt(x33); CHECK(x72);
  double x73=1.0/x3; CHECK(x73);
  double x74=ut1*ut2*x14*x73; CHECK(x74);
  double x75=1.0/x10; CHECK(x75);
  double x76=x6*x74; CHECK(x76);
  double x77=(1.0/2.0)*x8*(2*x16*x73 + 2*x5); CHECK(x77);
  double x78 = NAN;
  double x79=sqrt(x35); CHECK(x79);
  double x80 = NAN;
  double x81=ut2*x14 + x7*x73; CHECK(x81);
  double x82=x75*(x76 + x77); CHECK(x82);
  if (x12)
  {
    x78=0.5*x6*(-x76 - x77)/pow(x9, 3.0/2.0) + 0.5*x74*x75; XCHECK(x78);
    x80=0.5*x6*x75; XCHECK(x80);
  };
  if (x11)
  {
    x78=0; XCHECK(x78);
    x80=0.5*random1*x29; XCHECK(x80);
  };
  if ((x18 <= 0) && (x13))
  {
    result[10]=0.0; XCHECK(result[10]);
  };
  if ((x13 && x30 && (x11 || x30) && (x12 || x30) && (x30 || x34) && (x30 || x36) && (x30 || x37) && (x11 || x12 || x30) && (x11 || x12 || x34) && (x11 || x12 || x36) && (x11 || x12 || x37) && (x11 || x30 || x34) && (x11 || x30 || x36) && (x11 || x30 || x37) && (x11 || x34 || x36) && (x11 || x34 || x37) && (x12 || x30 || x34) && (x12 || x30 || x36) && (x12 || x30 || x37) && (x30 || x34 || x36) && (x30 || x34 || x37)) && (x12))
  {
    result[10]=(7216878369.0*mu*x24 - mu*x27 + 1020620730.0*x20*x24 + x20*x26 + 1443375674.0*x24*x28 + 3061862179.0*x24 - x26 - x27*x28)/(707106781000000.0*mu*x22 + 1.0e+15*x22); XCHECK(result[10]);
  };
  if ((x13 && x30 && (x11 || x30) && (x12 || x30) && (x30 || x34) && (x30 || x36) && (x30 || x37) && (x11 || x12 || x30) && (x11 || x12 || x34) && (x11 || x12 || x36) && (x11 || x12 || x37) && (x11 || x30 || x34) && (x11 || x30 || x36) && (x11 || x30 || x37) && (x11 || x34 || x36) && (x11 || x34 || x37) && (x12 || x30 || x34) && (x12 || x30 || x36) && (x12 || x30 || x37) && (x30 || x34 || x36) && (x30 || x34 || x37)) && (x11))
  {
    result[10]=x29*(-0.4082482906*mu*random1 - 1.1547005388*random1*x20)/x21; XCHECK(result[10]);
  };
  if ((x12 && x30 && x37 && x71 && x3 <= 0) && (x12))
  {
    result[10]=(-x47*x53 + x47*x54 + x47*x55 - x47*x56 + x47*x60 + x47*x62 + x47*x64 + x47*x65 + x47*x66 - x47*x67 + x47*x68 + x48*x53 - x48*x54 - x48*x55 + x48*x56 - x48*x60 - x48*x62 - x48*x64 - x48*x65 - x48*x66 + x48*x67 + x48*x68)/(x38*x49 + x41*x49); XCHECK(result[10]);
  };
  if ((x12 && x30 && x37 && x71 && x3 <= 0) && (x11))
  {
    result[10]=x29*(-x47*x69 - x47*x70 - x48*x69 + x48*x70)/(x45*x47*x48); XCHECK(result[10]);
  };
  if (x12 && x30 && x37 && x71 && x3 > 0)
  {
    result[10]=x72*x78 - x78*x79 - x80*(x81 + x82)/x79 + x80*(x81 - x82)/x72; XCHECK(result[10]);
  };
}
{
  double x0=mu*rn; CHECK(x0);
  double x1=pow(ut1, 2); CHECK(x1);
  double x2=pow(ut2, 2); CHECK(x2);
  double x3=mu*sqrt(x1 + x2) + un; CHECK(x3);
  double x4=mu*x3; CHECK(x4);
  double x5=rt1*x0 + ut1*x4; CHECK(x5);
  double x6=rt2*x0 + ut2*x4; CHECK(x6);
  double x7=pow(x5, 2) + pow(x6, 2); CHECK(x7);
  double x8=sqrt(x7); CHECK(x8);
  int    x9=x8 <= 0;
  int    x10=x8 > 0;
  double x11=pow(rt1, 2); CHECK(x11);
  double x12=pow(rt2, 2); CHECK(x12);
  double x13=pow(mu, 2); CHECK(x13);
  double x14=pow(rn, 2)*x13; CHECK(x14);
  double x15=sqrt(x11 + x12 + x14); CHECK(x15);
  double x16=pow(pow(random1, 2) + pow(random2, 2), -1.0/2.0); CHECK(x16);
  double x17=2.828427124*mu; CHECK(x17);
  double x18=pow(mu, 2); CHECK(x18);
  double x19=7071067810.0*mu; CHECK(x19);
  double x20=sqrt(2500000000.0*x18 - x19 + 7500000000.0); CHECK(x20);
  double x21=mu*x20; CHECK(x21);
  double x22=sqrt(7500000000.0*x18 + x19 + 2500000000.0); CHECK(x22);
  double x23=mu*x22; CHECK(x23);
  double x24=x18*x20; CHECK(x24);
  double x25=x18*x22; CHECK(x25);
  double x26=x16*(-8.164965812e-6*random1*x21 - 1.414213562e-5*random1*x23 - 5.773502694e-6*random1*x24 + 1.0e-5*random1*x25)/(sqrt(-x17 + x18 + 3.0)*sqrt(x17 + 3.0*x18 + 1.0)); CHECK(x26);
  double x27=x1*x13 + x11 + x12 + x13*x2 + x14 + pow(x3, 2); CHECK(x27);
  double x28=2*x8; CHECK(x28);
  double x29=x27 - x28; CHECK(x29);
  double x30=x27 + x28; CHECK(x30);
  double x31=sqrt(x29); CHECK(x31);
  double x32=mu*rt1; CHECK(x32);
  double x33=1.0/x8; CHECK(x33);
  double x34=x32*x5; CHECK(x34);
  double x35=mu*rt2*x6; CHECK(x35);
  double x36 = NAN;
  double x37=sqrt(x30); CHECK(x37);
  double x38 = NAN;
  double x39=rn*x13; CHECK(x39);
  double x40=x33*(x34 + x35); CHECK(x40);
  if (x10)
  {
    x36=0.5*x32*x33 + 0.5*x5*(-x34 - x35)/pow(x7, 3.0/2.0); XCHECK(x36);
    x38=0.5*x33*x5; XCHECK(x38);
  };
  if (x9)
  {
    x36=0; XCHECK(x36);
    x38=0.5*random1*x16; XCHECK(x38);
    result[13]=x26; XCHECK(result[13]);
  };
  if ((x15 <= 0) && (x10 || x9))
  {
    result[13]=0.0; XCHECK(result[13]);
  };
  if ((x29 <= 0 || x30 <= 0) && (x10))
  {
    result[13]=(-14433.75672*x21 - 24999.99999*x23 - 10206.20726*x24 + 17677.66952*x25)/(x20*x22); XCHECK(result[13]);
  };
  if ((x29 <= 0 || x30 <= 0) && (x9))
  {
    result[13]=x26; XCHECK(result[13]);
  };
  if (x10 && x15 > 0 && x29 > 0 && x30 > 0)
  {
    result[13]=x31*x36 - x36*x37 - x38*(x39 + x40)/x37 + x38*(x39 - x40)/x31; XCHECK(result[13]);
  };
}
{
  double x0=mu*rn; CHECK(x0);
  double x1=pow(ut1, 2); CHECK(x1);
  double x2=pow(ut2, 2); CHECK(x2);
  double x3=mu*sqrt(x1 + x2) + un; CHECK(x3);
  double x4=mu*x3; CHECK(x4);
  double x5=rt1*x0 + ut1*x4; CHECK(x5);
  double x6=pow(x5, 2); CHECK(x6);
  double x7=x6 + pow(rt2*x0 + ut2*x4, 2); CHECK(x7);
  double x8=sqrt(x7); CHECK(x8);
  int    x9=x8 <= 0;
  int    x10=x8 > 0;
  double x11=pow(rt1, 2); CHECK(x11);
  double x12=pow(rt2, 2); CHECK(x12);
  double x13=pow(mu, 2); CHECK(x13);
  double x14=pow(rn, 2)*x13; CHECK(x14);
  double x15=sqrt(x11 + x12 + x14); CHECK(x15);
  double x16=sqrt(pow(random1, 2) + pow(random2, 2)); CHECK(x16);
  double x17=1.0/x16; CHECK(x17);
  double x18=2.828427124*mu; CHECK(x18);
  double x19=pow(mu, 2); CHECK(x19);
  double x20=sqrt(-x18 + x19 + 3.0); CHECK(x20);
  double x21=sqrt(x18 + 3.0*x19 + 1.0); CHECK(x21);
  double x22=7071067810.0*mu; CHECK(x22);
  double x23=sqrt(2500000000.0*x19 - x22 + 7500000000.0); CHECK(x23);
  double x24=sqrt(7500000000.0*x19 + x22 + 2500000000.0); CHECK(x24);
  double x25=mu*x23; CHECK(x25);
  double x26=mu*x24; CHECK(x26);
  double x27=x23*x24; CHECK(x27);
  double x28=x17*(-5.773502694e-6*random1*x23 + 1.0e-5*random1*x24 - 4.082482906e-6*random1*x25 - 7.07106781e-6*random1*x26 + 4.0e-10*x16*x27)/(x20*x21); CHECK(x28);
  double x29=x20*x21; CHECK(x29);
  double x30=x1*x13 + x11 + x12 + x13*x2 + x14 + pow(x3, 2); CHECK(x30);
  double x31=2*x8; CHECK(x31);
  double x32=x30 - x31; CHECK(x32);
  double x33=x30 + x31; CHECK(x33);
  double x34=sqrt(x32); CHECK(x34);
  double x35=1.0/x8; CHECK(x35);
  double x36=x0*x35; CHECK(x36);
  double x37 = NAN;
  double x38=sqrt(x33); CHECK(x38);
  double x39=x36*x5; CHECK(x39);
  double x40 = NAN;
  if (x10)
  {
    x37=-0.5*x0*x6/pow(x7, 3.0/2.0) + 0.5*x36; XCHECK(x37);
    x40=0.5*x35*x5; XCHECK(x40);
  };
  if (x9)
  {
    x37=0; XCHECK(x37);
    x40=0.5*random1*x17; XCHECK(x40);
    result[16]=x28; XCHECK(result[16]);
  };
  if ((x15 <= 0) && (x10 || x9))
  {
    result[16]=1.00000000000000; XCHECK(result[16]);
  };
  if ((x32 <= 0 || x33 <= 0) && (x10))
  {
    result[16]=(-28067069980000.0*x19*x23 - 4419417380000.0*x19*x24 - 17860862700000.0*x23 + 707106781.0*x24*x25 + 30935921680000.0*x24 - 36084391820000.0*x25 - 12500000000000.0*x26 + 1000000000.0*x27)/(1.767766952e+18*mu*x29 + 2.5e+18*x29); XCHECK(result[16]);
  };
  if ((x32 <= 0 || x33 <= 0) && (x9))
  {
    result[16]=x28; XCHECK(result[16]);
  };
  if (x10 && x15 > 0 && x32 > 0 && x33 > 0)
  {
    result[16]=x34*x37 - x37*x38 + 1 - x40*(rt1 + x39)/x38 + x40*(rt1 - x39)/x34; XCHECK(result[16]);
  };
}
{
  double x0=mu*rn; CHECK(x0);
  double x1=pow(ut1, 2); CHECK(x1);
  double x2=pow(ut2, 2); CHECK(x2);
  double x3=mu*sqrt(x1 + x2) + un; CHECK(x3);
  double x4=mu*x3; CHECK(x4);
  double x5=rt1*x0 + ut1*x4; CHECK(x5);
  double x6=rt2*x0 + ut2*x4; CHECK(x6);
  double x7=pow(x5, 2) + pow(x6, 2); CHECK(x7);
  double x8=sqrt(x7); CHECK(x8);
  int    x9=x8 <= 0;
  int    x10=x8 > 0;
  double x11=pow(rt1, 2); CHECK(x11);
  double x12=pow(rt2, 2); CHECK(x12);
  double x13=pow(mu, 2); CHECK(x13);
  double x14=pow(rn, 2)*x13; CHECK(x14);
  double x15=sqrt(x11 + x12 + x14); CHECK(x15);
  double x16=pow(pow(random1, 2) + pow(random2, 2), -1.0/2.0); CHECK(x16);
  double x17=2.828427124*mu; CHECK(x17);
  double x18=pow(mu, 2); CHECK(x18);
  double x19=7071067810.0*mu; CHECK(x19);
  double x20=sqrt(2500000000.0*x18 - x19 + 7500000000.0); CHECK(x20);
  double x21=sqrt(7500000000.0*x18 + x19 + 2500000000.0); CHECK(x21);
  double x22=mu*x20; CHECK(x22);
  double x23=mu*x21; CHECK(x23);
  double x24=x16*(-5.773502694e-6*random1*x20 + 1.0e-5*random1*x21 - 4.082482906e-6*random1*x22 - 7.07106781e-6*random1*x23)/(sqrt(-x17 + x18 + 3.0)*sqrt(x17 + 3.0*x18 + 1.0)); CHECK(x24);
  double x25=x1*x13 + x11 + x12 + x13*x2 + x14 + pow(x3, 2); CHECK(x25);
  double x26=2*x8; CHECK(x26);
  double x27=x25 - x26; CHECK(x27);
  double x28=x25 + x26; CHECK(x28);
  double x29=sqrt(x27); CHECK(x29);
  double x30=mu*rn*x6; CHECK(x30);
  double x31 = NAN;
  double x32=sqrt(x28); CHECK(x32);
  double x33=1.0/x8; CHECK(x33);
  double x34=x30*x33; CHECK(x34);
  double x35 = NAN;
  if (x10)
  {
    x31=-0.5*x30*x5/pow(x7, 3.0/2.0); XCHECK(x31);
    x35=0.5*x33*x5; XCHECK(x35);
  };
  if (x9)
  {
    x31=0; XCHECK(x31);
    x35=0.5*random1*x16; XCHECK(x35);
    result[19]=x24; XCHECK(result[19]);
  };
  if ((x15 <= 0) && (x10 || x9))
  {
    result[19]=0.0; XCHECK(result[19]);
  };
  if ((x27 <= 0 || x28 <= 0) && (x10))
  {
    result[19]=-2000000*(-17860862710.0*x18*x20 + 13258252145.0*x18*x21 + 2551551811.0*x20 - 4419417381.0*x21 - 7216878370.0*x22 - 12499999995.0*x23)/(2000000000000.0*x20*x21 + 1414213562000.0*x21*x22); XCHECK(result[19]);
  };
  if ((x27 <= 0 || x28 <= 0) && (x9))
  {
    result[19]=x24; XCHECK(result[19]);
  };
  if (x10 && x15 > 0 && x27 > 0 && x28 > 0)
  {
    result[19]=x29*x31 - x31*x32 - x35*(rt2 + x34)/x32 + x35*(rt2 - x34)/x29; XCHECK(result[19]);
  };
}
{
  double x0=mu*ut2; CHECK(x0);
  double x1=pow(mu, 2); CHECK(x1);
  double x2=pow(ut1, 2); CHECK(x2);
  double x3=pow(ut2, 2); CHECK(x3);
  double x4=mu*sqrt(x2 + x3) + un; CHECK(x4);
  double x5=pow(rn, 2)*x1 + pow(rt1, 2) + pow(rt2, 2) + x1*x2 + x1*x3 + pow(x4, 2); CHECK(x5);
  double x6=mu*rn; CHECK(x6);
  double x7=rt2*x6 + x0*x4; CHECK(x7);
  double x8=sqrt(pow(x7, 2) + pow(mu*ut1*x4 + rt1*x6, 2)); CHECK(x8);
  double x9=2*x8; CHECK(x9);
  double x10 = NAN;
  if (x8 > 0)
  {
    x10=0.5*x7/x8; XCHECK(x10);
  };
  if (x8 <= 0)
  {
    x10=0.5*random2/sqrt(pow(random1, 2) + pow(random2, 2)); XCHECK(x10);
  };
  result[2]=rt2 + x0 + x10*sqrt(x5 - x9) - x10*sqrt(x5 + x9); XCHECK(result[2]);
}
{
  double x0=mu*rn; CHECK(x0);
  double x1=pow(ut1, 2); CHECK(x1);
  double x2=pow(ut2, 2); CHECK(x2);
  double x3=sqrt(x1 + x2); CHECK(x3);
  double x4=mu*x3 + un; CHECK(x4);
  double x5=mu*x4; CHECK(x5);
  double x6=rt1*x0 + ut1*x5; CHECK(x6);
  double x7=rt2*x0 + ut2*x5; CHECK(x7);
  double x8=pow(x6, 2) + pow(x7, 2); CHECK(x8);
  double x9=sqrt(x8); CHECK(x9);
  int    x10=x9 <= 0;
  int    x11=x9 > 0;
  double x12=x10 || x11; CHECK(x12);
  double x13=pow(mu, 2); CHECK(x13);
  double x14=x1*x13; CHECK(x14);
  double x15=x13*x2; CHECK(x15);
  double x16=pow(x4, 2); CHECK(x16);
  double x17=sqrt(x14 + x15 + x16); CHECK(x17);
  double x18=7071067810.0*mu; CHECK(x18);
  double x19=pow(mu, 2); CHECK(x19);
  double x20=sqrt(-x18 + 2500000000.0*x19 + 7500000000.0); CHECK(x20);
  double x21=sqrt(x18 + 7500000000.0*x19 + 2500000000.0); CHECK(x21);
  double x22=mu*x20; CHECK(x22);
  double x23=pow(pow(random1, 2) + pow(random2, 2), -1.0/2.0); CHECK(x23);
  double x24=2.828427124*mu; CHECK(x24);
  int    x25=x17 > 0;
  double x26=pow(rn, 2)*x13 + pow(rt1, 2) + pow(rt2, 2) + x14 + x15 + x16; CHECK(x26);
  double x27=2*x9; CHECK(x27);
  double x28=x26 - x27; CHECK(x28);
  int    x29=x28 <= 0;
  double x30=x26 + x27; CHECK(x30);
  int    x31=x30 <= 0;
  int    x32=x28 > 0;
  double x33=pow(rt1, 2); CHECK(x33);
  double x34=pow(rt2, 2); CHECK(x34);
  double x35=pow(rn, 2)*x19; CHECK(x35);
  double x36=pow(un, 2) + x33 + x34 + x35; CHECK(x36);
  double x37=sqrt(x33*x35 + x34*x35); CHECK(x37);
  double x38=2.0*x37; CHECK(x38);
  double x39=sqrt(x36 - x38); CHECK(x39);
  double x40=sqrt(x36 + x38); CHECK(x40);
  double x41=1/(x39*x40); CHECK(x41);
  double x42=0.5*mu*rn*rt2*un; CHECK(x42);
  double x43=0.5*random2*un; CHECK(x43);
  int    x44=x30 > 0;
  double x45=sqrt(x28); CHECK(x45);
  double x46=mu*ut2; CHECK(x46);
  double x47=1.0/x9; CHECK(x47);
  double x48=mu*ut1*x6; CHECK(x48);
  double x49=x46*x7; CHECK(x49);
  double x50 = NAN;
  double x51=sqrt(x30); CHECK(x51);
  double x52 = NAN;
  double x53=x47*(x48 + x49); CHECK(x53);
  if (x11)
  {
    x50=0.5*x46*x47 + 0.5*x7*(-x48 - x49)/pow(x8, 3.0/2.0); XCHECK(x50);
    x52=0.5*x47*x7; XCHECK(x52);
  };
  if (x10)
  {
    x50=0; XCHECK(x50);
    x52=0.5*random2*x23; XCHECK(x52);
  };
  if ((x17 <= 0) && (x12))
  {
    result[5]=0.0; XCHECK(result[5]);
  };
  if ((x12 && x25 && (x10 || x25) && (x11 || x25) && (x25 || x29) && (x25 || x31) && (x25 || x32) && (x10 || x11 || x25) && (x10 || x11 || x29) && (x10 || x11 || x31) && (x10 || x11 || x32) && (x10 || x25 || x29) && (x10 || x25 || x31) && (x10 || x25 || x32) && (x10 || x29 || x31) && (x10 || x29 || x32) && (x11 || x25 || x29) && (x11 || x25 || x31) && (x11 || x25 || x32) && (x25 || x29 || x31) && (x25 || x29 || x32)) && (x11))
  {
    result[5]=(-10206.20726*x20 + 17677.66952*x21 - 28867.51346*x22)/(x20*x21); XCHECK(result[5]);
  };
  if ((x12 && x25 && (x10 || x25) && (x11 || x25) && (x25 || x29) && (x25 || x31) && (x25 || x32) && (x10 || x11 || x25) && (x10 || x11 || x29) && (x10 || x11 || x31) && (x10 || x11 || x32) && (x10 || x25 || x29) && (x10 || x25 || x31) && (x10 || x25 || x32) && (x10 || x29 || x31) && (x10 || x29 || x32) && (x11 || x25 || x29) && (x11 || x25 || x31) && (x11 || x25 || x32) && (x25 || x29 || x31) && (x25 || x29 || x32)) && (x10))
  {
    result[5]=x23*(-5.773502694e-6*random2*x20 + 1.0e-5*random2*x21 - 1.6329931624e-5*random2*x22)/(sqrt(x19 - x24 + 3.0)*sqrt(3.0*x19 + x24 + 1.0)); XCHECK(result[5]);
  };
  if ((x11 && x25 && x32 && x44 && x3 <= 0) && (x11))
  {
    result[5]=x41*(-x39*x42 + x40*x42)/x37; XCHECK(result[5]);
  };
  if ((x11 && x25 && x32 && x44 && x3 <= 0) && (x10))
  {
    result[5]=x23*x41*(-x39*x43 + x40*x43); XCHECK(result[5]);
  };
  if (x11 && x25 && x32 && x44 && x3 > 0)
  {
    result[5]=x45*x50 - x50*x51 - x52*(x4 + x53)/x51 + x52*(x4 - x53)/x45; XCHECK(result[5]);
  };
}
{
  double x0=mu*rn; CHECK(x0);
  double x1=pow(ut1, 2); CHECK(x1);
  double x2=pow(ut2, 2); CHECK(x2);
  double x3=sqrt(x1 + x2); CHECK(x3);
  double x4=mu*x3 + un; CHECK(x4);
  double x5=mu*x4; CHECK(x5);
  double x6=ut1*x5; CHECK(x6);
  double x7=rt1*x0 + x6; CHECK(x7);
  double x8=rt2*x0 + ut2*x5; CHECK(x8);
  double x9=pow(x7, 2) + pow(x8, 2); CHECK(x9);
  double x10=sqrt(x9); CHECK(x10);
  int    x11=x10 <= 0;
  int    x12=x10 > 0;
  double x13=x11 || x12; CHECK(x13);
  double x14=pow(mu, 2); CHECK(x14);
  double x15=x1*x14; CHECK(x15);
  double x16=x14*x2; CHECK(x16);
  double x17=pow(x4, 2); CHECK(x17);
  double x18=sqrt(x15 + x16 + x17); CHECK(x18);
  double x19=2.828427124*mu; CHECK(x19);
  double x20=pow(mu, 2); CHECK(x20);
  double x21=sqrt(x19 + 3.0*x20 + 1.0); CHECK(x21);
  double x22=x21*sqrt(-x19 + x20 + 3.0); CHECK(x22);
  double x23=7071067810.0*mu; CHECK(x23);
  double x24=sqrt(2500000000.0*x20 - x23 + 7500000000.0); CHECK(x24);
  double x25=sqrt(7500000000.0*x20 + x23 + 2500000000.0); CHECK(x25);
  double x26=5303300858.0*x25; CHECK(x26);
  double x27=2500000000.0*x25; CHECK(x27);
  double x28=pow(mu, 3); CHECK(x28);
  double x29=pow(pow(random1, 2) + pow(random2, 2), -1.0/2.0); CHECK(x29);
  int    x30=x18 > 0;
  double x31=pow(rn, 2)*x14 + pow(rt1, 2) + pow(rt2, 2) + x15 + x16 + x17; CHECK(x31);
  double x32=2*x10; CHECK(x32);
  double x33=x31 - x32; CHECK(x33);
  int    x34=x33 <= 0;
  double x35=x31 + x32; CHECK(x35);
  int    x36=x35 <= 0;
  int    x37=x33 > 0;
  double x38=pow(rt1, 2); CHECK(x38);
  double x39=pow(rn, 2); CHECK(x39);
  double x40=x20*x39; CHECK(x40);
  double x41=pow(rt2, 2); CHECK(x41);
  double x42=x38*x40 + x40*x41; CHECK(x42);
  double x43=pow(x42, 3.0/2.0); CHECK(x43);
  double x44=pow(un, 2) + x38 + x40 + x41; CHECK(x44);
  double x45=sqrt(x42); CHECK(x45);
  double x46=2.0*x45; CHECK(x46);
  double x47=sqrt(x44 - x46); CHECK(x47);
  double x48=sqrt(x44 + x46); CHECK(x48);
  double x49=4.0*x43*x47*x48; CHECK(x49);
  double x50=pow(mu, 4); CHECK(x50);
  double x51=pow(rn, 3); CHECK(x51);
  double x52=pow(rt2, 5); CHECK(x52);
  double x53=1.414213562*un*x50*x51*x52; CHECK(x53);
  double x54=2.0*rt1*un*x28*x39*x52; CHECK(x54);
  double x55=2.0*pow(rt1, 5)*rt2*un*x28*x39; CHECK(x55);
  double x56=1.414213562*pow(rt1, 4)*rt2*un*x50*x51; CHECK(x56);
  double x57=pow(mu, 5); CHECK(x57);
  double x58=pow(rn, 4); CHECK(x58);
  double x59=pow(rt2, 3); CHECK(x59);
  double x60=2.0*rt1*un*x57*x58*x59; CHECK(x60);
  double x61=pow(rt1, 3); CHECK(x61);
  double x62=2.0*rt2*un*x57*x58*x61; CHECK(x62);
  double x63=pow(un, 3); CHECK(x63);
  double x64=2.0*rt1*x28*x39*x59*x63; CHECK(x64);
  double x65=2.0*rt2*x28*x39*x61*x63; CHECK(x65);
  double x66=4.0*un*x28*x39*x59*x61; CHECK(x66);
  double x67=2.828427124*un*x38*x50*x51*x59; CHECK(x67);
  double x68=2.0*mu*rt1*rt2*un*x43; CHECK(x68);
  double x69=0.5*random2*rn*rt1*un*x20; CHECK(x69);
  double x70=0.3535533905*mu*random2*un*x45; CHECK(x70);
  int    x71=x35 > 0;
  double x72=sqrt(x33); CHECK(x72);
  double x73=1.0/x3; CHECK(x73);
  double x74=ut1*ut2*x14*x73; CHECK(x74);
  double x75=1.0/x10; CHECK(x75);
  double x76=x74*x8; CHECK(x76);
  double x77=(1.0/2.0)*x7*(2*x15*x73 + 2*x5); CHECK(x77);
  double x78 = NAN;
  double x79=sqrt(x35); CHECK(x79);
  double x80 = NAN;
  double x81=ut1*x14 + x6*x73; CHECK(x81);
  double x82=x75*(x76 + x77); CHECK(x82);
  if (x12)
  {
    x78=0.5*x74*x75 + 0.5*x8*(-x76 - x77)/pow(x9, 3.0/2.0); XCHECK(x78);
    x80=0.5*x75*x8; XCHECK(x80);
  };
  if (x11)
  {
    x78=0; XCHECK(x78);
    x80=0.5*random2*x29; XCHECK(x80);
  };
  if ((x18 <= 0) && (x13))
  {
    result[8]=0.0; XCHECK(result[8]);
  };
  if ((x13 && x30 && (x11 || x30) && (x12 || x30) && (x30 || x34) && (x30 || x36) && (x30 || x37) && (x11 || x12 || x30) && (x11 || x12 || x34) && (x11 || x12 || x36) && (x11 || x12 || x37) && (x11 || x30 || x34) && (x11 || x30 || x36) && (x11 || x30 || x37) && (x11 || x34 || x36) && (x11 || x34 || x37) && (x12 || x30 || x34) && (x12 || x30 || x36) && (x12 || x30 || x37) && (x30 || x34 || x36) && (x30 || x34 || x37)) && (x12))
  {
    result[8]=(7216878369.0*mu*x24 - mu*x27 + 1020620730.0*x20*x24 + x20*x26 + 1443375674.0*x24*x28 + 3061862179.0*x24 - x26 - x27*x28)/(707106781000000.0*mu*x22 + 1.0e+15*x22); XCHECK(result[8]);
  };
  if ((x13 && x30 && (x11 || x30) && (x12 || x30) && (x30 || x34) && (x30 || x36) && (x30 || x37) && (x11 || x12 || x30) && (x11 || x12 || x34) && (x11 || x12 || x36) && (x11 || x12 || x37) && (x11 || x30 || x34) && (x11 || x30 || x36) && (x11 || x30 || x37) && (x11 || x34 || x36) && (x11 || x34 || x37) && (x12 || x30 || x34) && (x12 || x30 || x36) && (x12 || x30 || x37) && (x30 || x34 || x36) && (x30 || x34 || x37)) && (x11))
  {
    result[8]=x29*(-0.4082482906*mu*random2 - 1.1547005388*random2*x20)/x21; XCHECK(result[8]);
  };
  if ((x12 && x30 && x37 && x71 && x3 <= 0) && (x12))
  {
    result[8]=(-x47*x53 + x47*x54 + x47*x55 - x47*x56 + x47*x60 + x47*x62 + x47*x64 + x47*x65 + x47*x66 - x47*x67 + x47*x68 + x48*x53 - x48*x54 - x48*x55 + x48*x56 - x48*x60 - x48*x62 - x48*x64 - x48*x65 - x48*x66 + x48*x67 + x48*x68)/(x38*x49 + x41*x49); XCHECK(result[8]);
  };
  if ((x12 && x30 && x37 && x71 && x3 <= 0) && (x11))
  {
    result[8]=x29*(-x47*x69 - x47*x70 - x48*x69 + x48*x70)/(x45*x47*x48); XCHECK(result[8]);
  };
  if (x12 && x30 && x37 && x71 && x3 > 0)
  {
    result[8]=x72*x78 - x78*x79 - x80*(x81 + x82)/x79 + x80*(x81 - x82)/x72; XCHECK(result[8]);
  };
}
{
  double x0=mu*rn; CHECK(x0);
  double x1=pow(ut1, 2); CHECK(x1);
  double x2=pow(ut2, 2); CHECK(x2);
  double x3=sqrt(x1 + x2); CHECK(x3);
  double x4=mu*x3 + un; CHECK(x4);
  double x5=mu*x4; CHECK(x5);
  double x6=rt1*x0 + ut1*x5; CHECK(x6);
  double x7=ut2*x5; CHECK(x7);
  double x8=rt2*x0 + x7; CHECK(x8);
  double x9=pow(x6, 2) + pow(x8, 2); CHECK(x9);
  double x10=sqrt(x9); CHECK(x10);
  int    x11=x10 <= 0;
  int    x12=x10 > 0;
  double x13=x11 || x12; CHECK(x13);
  double x14=pow(mu, 2); CHECK(x14);
  double x15=x1*x14; CHECK(x15);
  double x16=x14*x2; CHECK(x16);
  double x17=pow(x4, 2); CHECK(x17);
  double x18=sqrt(x15 + x16 + x17); CHECK(x18);
  double x19=2.828427124*mu; CHECK(x19);
  double x20=pow(mu, 2); CHECK(x20);
  double x21=sqrt(x19 + 3.0*x20 + 1.0); CHECK(x21);
  double x22=x21*sqrt(-x19 + x20 + 3.0); CHECK(x22);
  double x23=7071067810.0*mu; CHECK(x23);
  double x24=sqrt(2500000000.0*x20 - x23 + 7500000000.0); CHECK(x24);
  double x25=sqrt(7500000000.0*x20 + x23 + 2500000000.0); CHECK(x25);
  double x26=13258252140000.0*x25; CHECK(x26);
  double x27=mu*x24; CHECK(x27);
  double x28=6250000000000.0*x25; CHECK(x28);
  double x29=x20*x24; CHECK(x29);
  double x30=pow(mu, 3); CHECK(x30);
  double x31=sqrt(pow(random1, 2) + pow(random2, 2)); CHECK(x31);
  double x32=1.0/x31; CHECK(x32);
  double x33=mu*x31; CHECK(x33);
  int    x34=x18 > 0;
  double x35=pow(rn, 2)*x14 + pow(rt1, 2) + pow(rt2, 2) + x15 + x16 + x17; CHECK(x35);
  double x36=2*x10; CHECK(x36);
  double x37=x35 - x36; CHECK(x37);
  int    x38=x37 <= 0;
  double x39=x35 + x36; CHECK(x39);
  int    x40=x39 <= 0;
  int    x41=x37 > 0;
  double x42=pow(rt1, 2); CHECK(x42);
  double x43=pow(rn, 2); CHECK(x43);
  double x44=x20*x43; CHECK(x44);
  double x45=pow(rt2, 2); CHECK(x45);
  double x46=x42*x44 + x44*x45; CHECK(x46);
  double x47=pow(x46, 3.0/2.0); CHECK(x47);
  double x48=pow(un, 2) + x42 + x44 + x45; CHECK(x48);
  double x49=sqrt(x46); CHECK(x49);
  double x50=2.0*x49; CHECK(x50);
  double x51=sqrt(x48 - x50); CHECK(x51);
  double x52=sqrt(x48 + x50); CHECK(x52);
  double x53=4.0*x47*x51*x52; CHECK(x53);
  double x54=x42*x53; CHECK(x54);
  double x55=x45*x53; CHECK(x55);
  double x56=2.0*pow(rt1, 6)*un*x30*x43; CHECK(x56);
  double x57=pow(mu, 4); CHECK(x57);
  double x58=pow(rn, 3); CHECK(x58);
  double x59=1.414213562*pow(rt2, 5)*un*x57*x58; CHECK(x59);
  double x60=pow(mu, 5); CHECK(x60);
  double x61=pow(rn, 4); CHECK(x61);
  double x62=pow(rt1, 4); CHECK(x62);
  double x63=2.0*un*x60*x61*x62; CHECK(x63);
  double x64=1.414213562*rt2*un*x57*x58*x62; CHECK(x64);
  double x65=pow(un, 3); CHECK(x65);
  double x66=2.0*x30*x43*x62*x65; CHECK(x66);
  double x67=2.0*pow(rt2, 4)*un*x30*x42*x43; CHECK(x67);
  double x68=4.0*un*x30*x43*x45*x62; CHECK(x68);
  double x69=2.828427124*pow(rt2, 3)*un*x42*x57*x58; CHECK(x69);
  double x70=2.0*un*x42*x45*x60*x61; CHECK(x70);
  double x71=2.0*x30*x42*x43*x45*x65; CHECK(x71);
  double x72=4.0*mu*un*x42*x47; CHECK(x72);
  double x73=2.0*mu*un*x45*x47; CHECK(x73);
  double x74=0.5*random2*rn*rt2*un*x20; CHECK(x74);
  double x75=0.3535533905*mu*random2*un*x49; CHECK(x75);
  int    x76=x39 > 0;
  double x77=sqrt(x37); CHECK(x77);
  double x78=1.0/x3; CHECK(x78);
  double x79=x16*x78; CHECK(x79);
  double x80=1.0/x10; CHECK(x80);
  double x81=ut2*x14; CHECK(x81);
  double x82=ut1*x6*x78*x81; CHECK(x82);
  double x83=(1.0/2.0)*x8*(2*x5 + 2*x79); CHECK(x83);
  double x84 = NAN;
  double x85=sqrt(x39); CHECK(x85);
  double x86 = NAN;
  double x87=x7*x78 + x81; CHECK(x87);
  double x88=x80*(x82 + x83); CHECK(x88);
  if (x12)
  {
    x84=0.5*x8*(-x82 - x83)/pow(x9, 3.0/2.0) + 0.5*x80*(x5 + x79); XCHECK(x84);
    x86=0.5*x8*x80; XCHECK(x86);
  };
  if (x11)
  {
    x84=0; XCHECK(x84);
    x86=0.5*random2*x32; XCHECK(x86);
  };
  if ((x18 <= 0) && (x13))
  {
    result[11]=mu; XCHECK(result[11]);
  };
  if ((x13 && x34 && (x11 || x34) && (x12 || x34) && (x34 || x38) && (x34 || x40) && (x34 || x41) && (x11 || x12 || x34) && (x11 || x12 || x38) && (x11 || x12 || x40) && (x11 || x12 || x41) && (x11 || x34 || x38) && (x11 || x34 || x40) && (x11 || x34 || x41) && (x11 || x38 || x40) && (x11 || x38 || x41) && (x12 || x34 || x38) && (x12 || x34 || x40) && (x12 || x34 || x41) && (x34 || x38 || x40) && (x34 || x38 || x41)) && (x12))
  {
    result[11]=(mu*x28 - x20*x26 - 61343466120000.0*x24*x30 - 7654655448000.0*x24 + 1000000000.0*x25*x27 + 707106781.0*x25*x29 + x26 - 46909709380000.0*x27 + x28*x30 - 104613624400000.0*x29)/(1.767766952e+18*mu*x22 + 2.5e+18*x22); XCHECK(result[11]);
  };
  if ((x13 && x34 && (x11 || x34) && (x12 || x34) && (x34 || x38) && (x34 || x40) && (x34 || x41) && (x11 || x12 || x34) && (x11 || x12 || x38) && (x11 || x12 || x40) && (x11 || x12 || x41) && (x11 || x34 || x38) && (x11 || x34 || x40) && (x11 || x34 || x41) && (x11 || x38 || x40) && (x11 || x38 || x41) && (x12 || x34 || x38) && (x12 || x34 || x40) && (x12 || x34 || x41) && (x34 || x38 || x40) && (x34 || x38 || x41)) && (x11))
  {
    result[11]=x32*(-0.4082482906*mu*random2 - 1.154700539*random2*x20 + 2.0e-5*x25*x33)/x21; XCHECK(result[11]);
  };
  if ((x12 && x34 && x41 && x76 && x3 <= 0) && (x12))
  {
    result[11]=(mu*x54 + mu*x55 - x51*x56 - x51*x59 - x51*x63 - x51*x64 - x51*x66 - x51*x67 - x51*x68 - x51*x69 - x51*x70 - x51*x71 - x51*x72 - x51*x73 + x52*x56 + x52*x59 + x52*x63 + x52*x64 + x52*x66 + x52*x67 + x52*x68 + x52*x69 + x52*x70 + x52*x71 - x52*x72 - x52*x73)/(x54 + x55); XCHECK(result[11]);
  };
  if ((x12 && x34 && x41 && x76 && x3 <= 0) && (x11))
  {
    result[11]=x32*(x33*x49*x51*x52 - x51*x74 - x51*x75 - x52*x74 + x52*x75)/(x49*x51*x52); XCHECK(result[11]);
  };
  if (x12 && x34 && x41 && x76 && x3 > 0)
  {
    result[11]=mu + x77*x84 - x84*x85 - x86*(x87 + x88)/x85 + x86*(x87 - x88)/x77; XCHECK(result[11]);
  };
}
{
  double x0=mu*rn; CHECK(x0);
  double x1=pow(ut1, 2); CHECK(x1);
  double x2=pow(ut2, 2); CHECK(x2);
  double x3=mu*sqrt(x1 + x2) + un; CHECK(x3);
  double x4=mu*x3; CHECK(x4);
  double x5=rt1*x0 + ut1*x4; CHECK(x5);
  double x6=rt2*x0 + ut2*x4; CHECK(x6);
  double x7=pow(x5, 2) + pow(x6, 2); CHECK(x7);
  double x8=sqrt(x7); CHECK(x8);
  int    x9=x8 <= 0;
  int    x10=x8 > 0;
  double x11=pow(rt1, 2); CHECK(x11);
  double x12=pow(rt2, 2); CHECK(x12);
  double x13=pow(mu, 2); CHECK(x13);
  double x14=pow(rn, 2)*x13; CHECK(x14);
  double x15=sqrt(x11 + x12 + x14); CHECK(x15);
  double x16=pow(pow(random1, 2) + pow(random2, 2), -1.0/2.0); CHECK(x16);
  double x17=2.828427124*mu; CHECK(x17);
  double x18=pow(mu, 2); CHECK(x18);
  double x19=7071067810.0*mu; CHECK(x19);
  double x20=sqrt(2500000000.0*x18 - x19 + 7500000000.0); CHECK(x20);
  double x21=mu*x20; CHECK(x21);
  double x22=sqrt(7500000000.0*x18 + x19 + 2500000000.0); CHECK(x22);
  double x23=mu*x22; CHECK(x23);
  double x24=x18*x20; CHECK(x24);
  double x25=x18*x22; CHECK(x25);
  double x26=x16*(-8.164965812e-6*random2*x21 - 1.414213562e-5*random2*x23 - 5.773502694e-6*random2*x24 + 1.0e-5*random2*x25)/(sqrt(-x17 + x18 + 3.0)*sqrt(x17 + 3.0*x18 + 1.0)); CHECK(x26);
  double x27=x1*x13 + x11 + x12 + x13*x2 + x14 + pow(x3, 2); CHECK(x27);
  double x28=2*x8; CHECK(x28);
  double x29=x27 - x28; CHECK(x29);
  double x30=x27 + x28; CHECK(x30);
  double x31=sqrt(x29); CHECK(x31);
  double x32=mu*rt2; CHECK(x32);
  double x33=1.0/x8; CHECK(x33);
  double x34=mu*rt1*x5; CHECK(x34);
  double x35=x32*x6; CHECK(x35);
  double x36 = NAN;
  double x37=sqrt(x30); CHECK(x37);
  double x38 = NAN;
  double x39=rn*x13; CHECK(x39);
  double x40=x33*(x34 + x35); CHECK(x40);
  if (x10)
  {
    x36=0.5*x32*x33 + 0.5*x6*(-x34 - x35)/pow(x7, 3.0/2.0); XCHECK(x36);
    x38=0.5*x33*x6; XCHECK(x38);
  };
  if (x9)
  {
    x36=0; XCHECK(x36);
    x38=0.5*random2*x16; XCHECK(x38);
    result[14]=x26; XCHECK(result[14]);
  };
  if ((x15 <= 0) && (x10 || x9))
  {
    result[14]=0.0; XCHECK(result[14]);
  };
  if ((x29 <= 0 || x30 <= 0) && (x10))
  {
    result[14]=(-14433.75672*x21 - 24999.99999*x23 - 10206.20726*x24 + 17677.66952*x25)/(x20*x22); XCHECK(result[14]);
  };
  if ((x29 <= 0 || x30 <= 0) && (x9))
  {
    result[14]=x26; XCHECK(result[14]);
  };
  if (x10 && x15 > 0 && x29 > 0 && x30 > 0)
  {
    result[14]=x31*x36 - x36*x37 - x38*(x39 + x40)/x37 + x38*(x39 - x40)/x31; XCHECK(result[14]);
  };
}
{
  double x0=mu*rn; CHECK(x0);
  double x1=pow(ut1, 2); CHECK(x1);
  double x2=pow(ut2, 2); CHECK(x2);
  double x3=mu*sqrt(x1 + x2) + un; CHECK(x3);
  double x4=mu*x3; CHECK(x4);
  double x5=rt1*x0 + ut1*x4; CHECK(x5);
  double x6=rt2*x0 + ut2*x4; CHECK(x6);
  double x7=pow(x5, 2) + pow(x6, 2); CHECK(x7);
  double x8=sqrt(x7); CHECK(x8);
  int    x9=x8 <= 0;
  int    x10=x8 > 0;
  double x11=pow(rt1, 2); CHECK(x11);
  double x12=pow(rt2, 2); CHECK(x12);
  double x13=pow(mu, 2); CHECK(x13);
  double x14=pow(rn, 2)*x13; CHECK(x14);
  double x15=sqrt(x11 + x12 + x14); CHECK(x15);
  double x16=pow(pow(random1, 2) + pow(random2, 2), -1.0/2.0); CHECK(x16);
  double x17=2.828427124*mu; CHECK(x17);
  double x18=pow(mu, 2); CHECK(x18);
  double x19=7071067810.0*mu; CHECK(x19);
  double x20=sqrt(2500000000.0*x18 - x19 + 7500000000.0); CHECK(x20);
  double x21=sqrt(7500000000.0*x18 + x19 + 2500000000.0); CHECK(x21);
  double x22=mu*x20; CHECK(x22);
  double x23=mu*x21; CHECK(x23);
  double x24=x16*(-5.773502694e-6*random2*x20 + 1.0e-5*random2*x21 - 4.082482906e-6*random2*x22 - 7.07106781e-6*random2*x23)/(sqrt(-x17 + x18 + 3.0)*sqrt(x17 + 3.0*x18 + 1.0)); CHECK(x24);
  double x25=x1*x13 + x11 + x12 + x13*x2 + x14 + pow(x3, 2); CHECK(x25);
  double x26=2*x8; CHECK(x26);
  double x27=x25 - x26; CHECK(x27);
  double x28=x25 + x26; CHECK(x28);
  double x29=sqrt(x27); CHECK(x29);
  double x30=mu*rn*x5; CHECK(x30);
  double x31 = NAN;
  double x32=sqrt(x28); CHECK(x32);
  double x33=1.0/x8; CHECK(x33);
  double x34=x30*x33; CHECK(x34);
  double x35 = NAN;
  if (x10)
  {
    x31=-0.5*x30*x6/pow(x7, 3.0/2.0); XCHECK(x31);
    x35=0.5*x33*x6; XCHECK(x35);
  };
  if (x9)
  {
    x31=0; XCHECK(x31);
    x35=0.5*random2*x16; XCHECK(x35);
    result[17]=x24; XCHECK(result[17]);
  };
  if ((x15 <= 0) && (x10 || x9))
  {
    result[17]=0.0; XCHECK(result[17]);
  };
  if ((x27 <= 0 || x28 <= 0) && (x10))
  {
    result[17]=-2000000*(-17860862710.0*x18*x20 + 13258252145.0*x18*x21 + 2551551811.0*x20 - 4419417381.0*x21 - 7216878370.0*x22 - 12499999995.0*x23)/(2000000000000.0*x20*x21 + 1414213562000.0*x21*x22); XCHECK(result[17]);
  };
  if ((x27 <= 0 || x28 <= 0) && (x9))
  {
    result[17]=x24; XCHECK(result[17]);
  };
  if (x10 && x15 > 0 && x27 > 0 && x28 > 0)
  {
    result[17]=x29*x31 - x31*x32 - x35*(rt1 + x34)/x32 + x35*(rt1 - x34)/x29; XCHECK(result[17]);
  };
}
{
  double x0=mu*rn; CHECK(x0);
  double x1=pow(ut1, 2); CHECK(x1);
  double x2=pow(ut2, 2); CHECK(x2);
  double x3=mu*sqrt(x1 + x2) + un; CHECK(x3);
  double x4=mu*x3; CHECK(x4);
  double x5=rt2*x0 + ut2*x4; CHECK(x5);
  double x6=pow(x5, 2); CHECK(x6);
  double x7=x6 + pow(rt1*x0 + ut1*x4, 2); CHECK(x7);
  double x8=sqrt(x7); CHECK(x8);
  int    x9=x8 <= 0;
  int    x10=x8 > 0;
  double x11=pow(rt1, 2); CHECK(x11);
  double x12=pow(rt2, 2); CHECK(x12);
  double x13=pow(mu, 2); CHECK(x13);
  double x14=pow(rn, 2)*x13; CHECK(x14);
  double x15=sqrt(x11 + x12 + x14); CHECK(x15);
  double x16=sqrt(pow(random1, 2) + pow(random2, 2)); CHECK(x16);
  double x17=1.0/x16; CHECK(x17);
  double x18=2.828427124*mu; CHECK(x18);
  double x19=pow(mu, 2); CHECK(x19);
  double x20=sqrt(-x18 + x19 + 3.0); CHECK(x20);
  double x21=sqrt(x18 + 3.0*x19 + 1.0); CHECK(x21);
  double x22=7071067810.0*mu; CHECK(x22);
  double x23=sqrt(2500000000.0*x19 - x22 + 7500000000.0); CHECK(x23);
  double x24=sqrt(7500000000.0*x19 + x22 + 2500000000.0); CHECK(x24);
  double x25=mu*x23; CHECK(x25);
  double x26=mu*x24; CHECK(x26);
  double x27=x23*x24; CHECK(x27);
  double x28=x17*(-5.773502694e-6*random2*x23 + 1.0e-5*random2*x24 - 4.082482906e-6*random2*x25 - 7.07106781e-6*random2*x26 + 4.0e-10*x16*x27)/(x20*x21); CHECK(x28);
  double x29=x20*x21; CHECK(x29);
  double x30=x1*x13 + x11 + x12 + x13*x2 + x14 + pow(x3, 2); CHECK(x30);
  double x31=2*x8; CHECK(x31);
  double x32=x30 - x31; CHECK(x32);
  double x33=x30 + x31; CHECK(x33);
  double x34=sqrt(x32); CHECK(x34);
  double x35=1.0/x8; CHECK(x35);
  double x36=x0*x35; CHECK(x36);
  double x37 = NAN;
  double x38=sqrt(x33); CHECK(x38);
  double x39=x36*x5; CHECK(x39);
  double x40 = NAN;
  if (x10)
  {
    x37=-0.5*x0*x6/pow(x7, 3.0/2.0) + 0.5*x36; XCHECK(x37);
    x40=0.5*x35*x5; XCHECK(x40);
  };
  if (x9)
  {
    x37=0; XCHECK(x37);
    x40=0.5*random2*x17; XCHECK(x40);
    result[20]=x28; XCHECK(result[20]);
  };
  if ((x15 <= 0) && (x10 || x9))
  {
    result[20]=1.00000000000000; XCHECK(result[20]);
  };
  if ((x32 <= 0 || x33 <= 0) && (x10))
  {
    result[20]=(-28067069980000.0*x19*x23 - 4419417380000.0*x19*x24 - 17860862700000.0*x23 + 707106781.0*x24*x25 + 30935921680000.0*x24 - 36084391820000.0*x25 - 12500000000000.0*x26 + 1000000000.0*x27)/(1.767766952e+18*mu*x29 + 2.5e+18*x29); XCHECK(result[20]);
  };
  if ((x32 <= 0 || x33 <= 0) && (x9))
  {
    result[20]=x28; XCHECK(result[20]);
  };
  if (x10 && x15 > 0 && x32 > 0 && x33 > 0)
  {
    result[20]=x34*x37 - x37*x38 + 1 - x40*(rt2 + x39)/x38 + x40*(rt2 - x39)/x34; XCHECK(result[20]);
  };
}
}
void frictionContact3D_FischerBurmeisterFGenerated(
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
  double x0=pow(ut1, 2); CHECK(x0);
  double x1=pow(ut2, 2); CHECK(x1);
  double x2=mu*sqrt(x0 + x1) + un; CHECK(x2);
  double x3=mu*rn; CHECK(x3);
  double x4=pow(mu, 2); CHECK(x4);
  double x5=pow(rn, 2)*x4 + pow(rt1, 2) + pow(rt2, 2) + x0*x4 + x1*x4 + pow(x2, 2); CHECK(x5);
  double x6=mu*ut1; CHECK(x6);
  double x7=rt1*x3 + x2*x6; CHECK(x7);
  double x8=mu*ut2; CHECK(x8);
  double x9=rt2*x3 + x2*x8; CHECK(x9);
  double x10=sqrt(pow(x7, 2) + pow(x9, 2)); CHECK(x10);
  double x11=2*x10; CHECK(x11);
  double x12=0.5*sqrt(-x11 + x5); CHECK(x12);
  double x13=0.5*sqrt(x11 + x5); CHECK(x13);
  double x14=1.0/x10; CHECK(x14);
  int    x15=x10 > 0;
  double x16=pow(pow(random1, 2) + pow(random2, 2), -1.0/2.0); CHECK(x16);
  int    x17=x10 <= 0;
  double x18 = NAN;
  double x19 = NAN;
  if (x15)
  {
    x18=x14*x7; XCHECK(x18);
    x19=x14*x9; XCHECK(x19);
  };
  if (x17)
  {
    x18=random1*x16; XCHECK(x18);
    x19=random2*x16; XCHECK(x19);
  };
  result[0]=-x12 - x13 + x2 + x3; XCHECK(result[0]);
  result[1]=rt1 + x12*x18 - x13*x18 + x6; XCHECK(result[1]);
  result[2]=rt2 + x12*x19 - x13*x19 + x8; XCHECK(result[2]);
}
void frictionContact3D_FischerBurmeisterABGenerated(
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
  double x0=pow(mu, 2); CHECK(x0);
  double x1=pow(ut1, 2); CHECK(x1);
  double x2=x0*x1; CHECK(x2);
  double x3=pow(ut2, 2); CHECK(x3);
  double x4=x0*x3; CHECK(x4);
  double x5=sqrt(x1 + x3); CHECK(x5);
  double x6=mu*x5 + un; CHECK(x6);
  double x7=pow(x6, 2); CHECK(x7);
  double x8=sqrt(x2 + x4 + x7); CHECK(x8);
  int    x9=x8 <= 0;
  double x10=2.828427124*mu; CHECK(x10);
  double x11=pow(mu, 2); CHECK(x11);
  double x12=sqrt(-x10 + x11 + 3.0); CHECK(x12);
  double x13=1.0/x12; CHECK(x13);
  double x14=sqrt(x10 + 3.0*x11 + 1.0); CHECK(x14);
  double x15=1.0/x14; CHECK(x15);
  double x16=x13*x15; CHECK(x16);
  double x17=7071067810.0*mu; CHECK(x17);
  double x18=sqrt(2500000000.0*x11 - x17 + 7500000000.0); CHECK(x18);
  double x19=5.773502694e-6*x18; CHECK(x19);
  double x20=sqrt(7500000000.0*x11 + x17 + 2500000000.0); CHECK(x20);
  double x21=1.0e-5*x20; CHECK(x21);
  double x22=-x19 - x21; CHECK(x22);
  double x23=mu*x18; CHECK(x23);
  double x24=1.6329931624e-5*x23; CHECK(x24);
  double x25=x18*x20; CHECK(x25);
  double x26=4.0e-10*x25; CHECK(x26);
  int    x27=x8 > 0;
  double x28=mu*rn; CHECK(x28);
  double x29=mu*x6; CHECK(x29);
  double x30=rt1*x28 + ut1*x29; CHECK(x30);
  double x31=pow(x30, 2); CHECK(x31);
  double x32=rt2*x28 + ut2*x29; CHECK(x32);
  double x33=pow(x32, 2); CHECK(x33);
  double x34=x31 + x33; CHECK(x34);
  double x35=sqrt(x34); CHECK(x35);
  int    x36=x35 <= 0;
  int    x37=x35 > 0;
  double x38=x36 || x37; CHECK(x38);
  double x39=pow(rt1, 2); CHECK(x39);
  double x40=pow(rt2, 2); CHECK(x40);
  double x41=pow(rn, 2)*x0; CHECK(x41);
  double x42=x2 + x39 + x4 + x40 + x41 + x7; CHECK(x42);
  double x43=2*x35; CHECK(x43);
  double x44=x42 - x43; CHECK(x44);
  int    x45=x44 <= 0;
  double x46=x42 + x43; CHECK(x46);
  int    x47=x46 <= 0;
  int    x48=x44 > 0;
  double x49=x36 || x45 || x47; CHECK(x49);
  double x50=x27 && x38 && x49 && (x27 || x36) && (x27 || x37) && (x27 || x45) && (x27 || x47) && (x27 || x48) && (x27 || x36 || x37) && (x27 || x36 || x45) && (x27 || x36 || x47) && (x27 || x36 || x48) && (x27 || x37 || x45) && (x27 || x37 || x47) && (x27 || x37 || x48) && (x27 || x45 || x47) && (x27 || x45 || x48) && (x36 || x37 || x45) && (x36 || x37 || x47) && (x36 || x37 || x48) && (x36 || x45 || x48); CHECK(x50);
  double x51=pow(rt1, 2); CHECK(x51);
  double x52=pow(rt2, 2); CHECK(x52);
  double x53=pow(rn, 2); CHECK(x53);
  double x54=x11*x53; CHECK(x54);
  double x55=pow(un, 2) + x51 + x52 + x54; CHECK(x55);
  double x56=x51*x54 + x52*x54; CHECK(x56);
  double x57=sqrt(x56); CHECK(x57);
  double x58=2.0*x57; CHECK(x58);
  double x59=sqrt(x55 - x58); CHECK(x59);
  double x60=1.0/x59; CHECK(x60);
  double x61=sqrt(x55 + x58); CHECK(x61);
  double x62=1.0/x61; CHECK(x62);
  double x63=0.5*un; CHECK(x63);
  double x64=x59*x63; CHECK(x64);
  double x65=x61*x63; CHECK(x65);
  int    x66=x46 > 0;
  double x67=x27 && x37 && x48 && x66 && x5 <= 0; CHECK(x67);
  double x68=sqrt(x46); CHECK(x68);
  double x69=0.5/x68; CHECK(x69);
  double x70=mu*x30; CHECK(x70);
  double x71=ut1*x70; CHECK(x71);
  double x72=mu*x32; CHECK(x72);
  double x73=ut2*x72; CHECK(x73);
  double x74=1.0/x35; CHECK(x74);
  double x75=x74*(x71 + x73); CHECK(x75);
  double x76=x69*(x6 + x75); CHECK(x76);
  double x77=sqrt(x44); CHECK(x77);
  double x78=0.5/x77; CHECK(x78);
  double x79=x78*(x6 - x75); CHECK(x79);
  double x80=x27 && x37 && x48 && x66 && x5 > 0; CHECK(x80);
  double x81=0.707106781*mu; CHECK(x81);
  double x82=0.4082482906*mu; CHECK(x82);
  double x83=1.154700539*x11; CHECK(x83);
  double x84=mu*x20; CHECK(x84);
  double x85=x15*(-x82 - x83 + 1.414213562e-5*x84); CHECK(x85);
  double x86=1.0/x57; CHECK(x86);
  double x87=x60*x62*x86; CHECK(x87);
  double x88=rn*rt1*x11; CHECK(x88);
  double x89=0.3535533905*mu*un*x57; CHECK(x89);
  double x90=x59*x89; CHECK(x90);
  double x91=x61*x89; CHECK(x91);
  double x92=x57*x59*x61; CHECK(x92);
  double x93=x81*x92 - x90 - x91; CHECK(x93);
  double x94=1.0/x5; CHECK(x94);
  double x95=mu*x94; CHECK(x95);
  double x96=ut1*x95; CHECK(x96);
  double x97=ut1*x0 + x6*x96; CHECK(x97);
  double x98=ut1*ut2*x0*x94; CHECK(x98);
  double x99=x32*x98; CHECK(x99);
  double x100=2*x29; CHECK(x100);
  double x101=x2*x94; CHECK(x101);
  double x102=(1.0/2.0)*x30*(x100 + 2*x101); CHECK(x102);
  double x103=x74*(x102 + x99); CHECK(x103);
  double x104=x69*(x103 + x97); CHECK(x104);
  double x105=x78*(-x103 + x97); CHECK(x105);
  double x106=rn*rt2*x11; CHECK(x106);
  double x107=ut2*x95; CHECK(x107);
  double x108=ut2*x0 + x107*x6; CHECK(x108);
  double x109=x30*x98; CHECK(x109);
  double x110=x4*x94; CHECK(x110);
  double x111=(1.0/2.0)*x32*(x100 + 2*x110); CHECK(x111);
  double x112=x74*(x109 + x111); CHECK(x112);
  double x113=x69*(x108 + x112); CHECK(x113);
  double x114=x78*(x108 - x112); CHECK(x114);
  double x115=sqrt(x39 + x40 + x41); CHECK(x115);
  int    x116=x115 <= 0;
  double x117=8.164965812e-6*x23; CHECK(x117);
  double x118=1.414213562e-5*x84; CHECK(x118);
  double x119=mu*x18*x20; CHECK(x119);
  double x120=rn*x0; CHECK(x120);
  double x121=rt1*x70; CHECK(x121);
  double x122=rt2*x72; CHECK(x122);
  double x123=x74*(x121 + x122); CHECK(x123);
  double x124=x69*(x120 + x123); CHECK(x124);
  double x125=x78*(x120 - x123); CHECK(x125);
  double x126=x37 && x48 && x66 && x115 > 0; CHECK(x126);
  double x127=4.082482906e-6*x23; CHECK(x127);
  double x128=7.07106781e-6*x84; CHECK(x128);
  double x129=x16*(-x127 + x128 + x22); CHECK(x129);
  double x130=x28*x74; CHECK(x130);
  double x131=x130*x30; CHECK(x131);
  double x132=x69*(rt1 + x131); CHECK(x132);
  double x133=x78*(rt1 - x131); CHECK(x133);
  double x134=x130*x32; CHECK(x134);
  double x135=x69*(rt2 + x134); CHECK(x135);
  double x136=x78*(rt2 - x134); CHECK(x136);
  double x137 = NAN;
  double x138=1/(x18*x20); CHECK(x138);
  double x139=10206.20726*x18; CHECK(x139);
  double x140=17677.66952*x20; CHECK(x140);
  double x141=x138*(-x139 + x140 - 28867.51346*x23); CHECK(x141);
  double x142=sqrt(pow(random1, 2) + pow(random2, 2)); CHECK(x142);
  double x143=1.0/x142; CHECK(x143);
  double x144=x13*x143*x15; CHECK(x144);
  double x145=random1*x19; CHECK(x145);
  double x146=random1*x21; CHECK(x146);
  double x147=-x145 + x146; CHECK(x147);
  double x148=mu*rn*rt1; CHECK(x148);
  double x149=random1*x64; CHECK(x149);
  double x150=random1*x65; CHECK(x150);
  double x151=x143*x60*x62; CHECK(x151);
  double x152=0.5*x77; CHECK(x152);
  double x153=mu*x74; CHECK(x153);
  double x154=-x71 - x73; CHECK(x154);
  double x155=pow(x34, -3.0/2.0); CHECK(x155);
  double x156=x155*x30; CHECK(x156);
  double x157 = NAN;
  double x158=0.5*x68; CHECK(x158);
  double x159 = NAN;
  double x160 = NAN;
  double x161=x12*x14; CHECK(x161);
  double x162=mu*x12*x14; CHECK(x162);
  double x163=1.0/(2.5e+18*x161 + 1.767766952e+18*x162); CHECK(x163);
  double x164=13258252140000.0*x20; CHECK(x164);
  double x165=x11*x18; CHECK(x165);
  double x166=pow(mu, 3); CHECK(x166);
  double x167=x166*x18; CHECK(x167);
  double x168=x166*x20; CHECK(x168);
  double x169=x163*(-x11*x164 + 1000000000.0*x119 + x164 + 707106781.0*x165*x20 - 104613624400000.0*x165 - 61343466120000.0*x167 + 6250000000000.0*x168 - 7654655448000.0*x18 - 46909709380000.0*x23 + 6250000000000.0*x84); CHECK(x169);
  double x170=x143*x15; CHECK(x170);
  double x171=-random1*x82; CHECK(x171);
  double x172=2.0e-5*x142*x84; CHECK(x172);
  double x173=pow(x56, 3.0/2.0); CHECK(x173);
  double x174=4.0*x173*x59*x61; CHECK(x174);
  double x175=x174*x51; CHECK(x175);
  double x176=x174*x52; CHECK(x176);
  double x177=1.0/(x175 + x176); CHECK(x177);
  double x178=pow(mu, 5); CHECK(x178);
  double x179=pow(rn, 4); CHECK(x179);
  double x180=2.0*un*x178*x179*x51*x52; CHECK(x180);
  double x181=pow(un, 3); CHECK(x181);
  double x182=2.0*x166*x181*x51*x52*x53; CHECK(x182);
  double x183=mu*x175 + mu*x176 - x180*x59 + x180*x61 - x182*x59 + x182*x61; CHECK(x183);
  double x184=2.0*pow(rt2, 6)*un*x166*x53; CHECK(x184);
  double x185=pow(mu, 4); CHECK(x185);
  double x186=pow(rn, 3); CHECK(x186);
  double x187=pow(rt1, 5); CHECK(x187);
  double x188=1.414213562*un*x185*x186*x187; CHECK(x188);
  double x189=-x188*x59; CHECK(x189);
  double x190=x188*x61; CHECK(x190);
  double x191=pow(rt2, 4); CHECK(x191);
  double x192=2.0*un*x178*x179*x191; CHECK(x192);
  double x193=1.414213562*rt1*un*x185*x186*x191; CHECK(x193);
  double x194=-x193*x59; CHECK(x194);
  double x195=x193*x61; CHECK(x195);
  double x196=2.0*x166*x181*x191*x53; CHECK(x196);
  double x197=un*x166*x191*x51*x53*x59; CHECK(x197);
  double x198=un*x166*x191*x51*x53*x61; CHECK(x198);
  double x199=pow(rt1, 4); CHECK(x199);
  double x200=un*x166*x199*x52*x53*x59; CHECK(x200);
  double x201=un*x166*x199*x52*x53*x61; CHECK(x201);
  double x202=pow(rt1, 3); CHECK(x202);
  double x203=2.828427124*un*x185*x186*x202*x52; CHECK(x203);
  double x204=-x203*x59; CHECK(x204);
  double x205=x203*x61; CHECK(x205);
  double x206=mu*un*x173*x51*x59; CHECK(x206);
  double x207=mu*un*x173*x51*x61; CHECK(x207);
  double x208=mu*un*x173*x52*x59; CHECK(x208);
  double x209=mu*un*x173*x52*x61; CHECK(x209);
  double x210=x143*x60*x62*x86; CHECK(x210);
  double x211=-random1*x90 + random1*x91; CHECK(x211);
  double x212=mu*x142*x92; CHECK(x212);
  double x213=-x102 - x99; CHECK(x213);
  double x214 = NAN;
  double x215=5303300858.0*x20; CHECK(x215);
  double x216=(x11*x215 + 1020620730.0*x165 + 1443375674.0*x167 - 2500000000.0*x168 + 3061862179.0*x18 - x215 + 7216878369.0*x23 - 2500000000.0*x84)/(1.0e+15*x161 + 707106781000000.0*x162); CHECK(x216);
  double x217=1.1547005388*x11; CHECK(x217);
  double x218=pow(rt2, 5); CHECK(x218);
  double x219=2.0*rt1*un*x166*x218*x53; CHECK(x219);
  double x220=2.0*rt2*un*x166*x187*x53; CHECK(x220);
  double x221=pow(rt2, 3); CHECK(x221);
  double x222=2.0*rt1*un*x178*x179*x221; CHECK(x222);
  double x223=2.0*rt2*un*x178*x179*x202; CHECK(x223);
  double x224=2.0*rt1*x166*x181*x221*x53; CHECK(x224);
  double x225=2.0*rt2*x166*x181*x202*x53; CHECK(x225);
  double x226=4.0*un*x166*x202*x221*x53; CHECK(x226);
  double x227=2.0*mu*rt1*rt2*un*x173; CHECK(x227);
  double x228=x219*x59 - x219*x61 + x220*x59 - x220*x61 + x222*x59 - x222*x61 + x223*x59 - x223*x61 + x224*x59 - x224*x61 + x225*x59 - x225*x61 + x226*x59 - x226*x61 + x227*x59 + x227*x61; CHECK(x228);
  double x229=x74*x98; CHECK(x229);
  double x230=-x109 - x111; CHECK(x230);
  double x231 = NAN;
  double x232=x144*(-random1*x117 - random1*x118 - x11*x145 + x11*x146); CHECK(x232);
  double x233=x138*(-x11*x139 + x11*x140 - 14433.75672*x23 - 24999.99999*x84); CHECK(x233);
  double x234=x45 || x47; CHECK(x234);
  double x235=-x121 - x122; CHECK(x235);
  double x236 = NAN;
  double x237 = NAN;
  double x238=x142*x26; CHECK(x238);
  double x239=-random1*x127 - random1*x128 + x147; CHECK(x239);
  double x240=x144*(x238 + x239); CHECK(x240);
  double x241=x11*x20; CHECK(x241);
  double x242=x163*(707106781.0*x119 - 28067069980000.0*x165 - 17860862700000.0*x18 + 30935921680000.0*x20 - 36084391820000.0*x23 - 4419417380000.0*x241 + 1000000000.0*x25 - 12500000000000.0*x84); CHECK(x242);
  double x243=mu*rn*x155; CHECK(x243);
  double x244 = NAN;
  double x245=x144*x239; CHECK(x245);
  double x246=-2000000*(-17860862710.0*x165 + 2551551811.0*x18 - 4419417381.0*x20 - 7216878370.0*x23 + 13258252145.0*x241 - 12499999995.0*x84)/(1414213562000.0*x119 + 2000000000000.0*x25); CHECK(x246);
  double x247 = NAN;
  double x248=-x247*x68 + x247*x77; CHECK(x248);
  double x249=random2*x19; CHECK(x249);
  double x250=random2*x21; CHECK(x250);
  double x251=-x249 + x250; CHECK(x251);
  double x252=mu*rn*rt2; CHECK(x252);
  double x253=random2*x64; CHECK(x253);
  double x254=random2*x65; CHECK(x254);
  double x255=x155*x32; CHECK(x255);
  double x256 = NAN;
  double x257 = NAN;
  double x258=-random2*x82; CHECK(x258);
  double x259=1.414213562*un*x185*x186*x218; CHECK(x259);
  double x260=-x259*x59; CHECK(x260);
  double x261=x259*x61; CHECK(x261);
  double x262=1.414213562*rt2*un*x185*x186*x199; CHECK(x262);
  double x263=-x262*x59; CHECK(x263);
  double x264=x262*x61; CHECK(x264);
  double x265=2.828427124*un*x185*x186*x221*x51; CHECK(x265);
  double x266=-x265*x59; CHECK(x266);
  double x267=x265*x61; CHECK(x267);
  double x268=-random2*x90 + random2*x91; CHECK(x268);
  double x269 = NAN;
  double x270=2.0*pow(rt1, 6)*un*x166*x53; CHECK(x270);
  double x271=2.0*un*x178*x179*x199; CHECK(x271);
  double x272=2.0*x166*x181*x199*x53; CHECK(x272);
  double x273 = NAN;
  double x274=x144*(-random2*x117 - random2*x118 - x11*x249 + x11*x250); CHECK(x274);
  double x275 = NAN;
  double x276=-random2*x127 - random2*x128 + x251; CHECK(x276);
  double x277=x144*x276; CHECK(x277);
  double x278=x144*(x238 + x276); CHECK(x278);
  double x279 = NAN;
  if (x38)
  {
    x137=0.0; XCHECK(x137);
    x160=mu; XCHECK(x160);
    x237=1.00000000000000; XCHECK(x237);
  };
  if (x37)
  {
    x157=ut1*x153 + x154*x156; XCHECK(x157);
    x159=x30*x74; XCHECK(x159);
    x214=x156*x213 + x74*(x101 + x29); XCHECK(x214);
    x231=x156*x230 + x229; XCHECK(x231);
    x236=rt1*x153 + x156*x235; XCHECK(x236);
    x244=x130 - x243*x31; XCHECK(x244);
    x247=-0.5*x243*x30*x32; XCHECK(x247);
    x256=ut2*x153 + x154*x255; XCHECK(x256);
    x257=x32*x74; XCHECK(x257);
    x269=x213*x255 + x229; XCHECK(x269);
    x273=x230*x255 + x74*(x110 + x29); XCHECK(x273);
    x275=rt2*x153 + x235*x255; XCHECK(x275);
    x279=x130 - x243*x33; XCHECK(x279);
  };
  if (x36)
  {
    x157=0; XCHECK(x157);
    x159=random1*x143; XCHECK(x159);
    x214=0; XCHECK(x214);
    x231=0; XCHECK(x231);
    x236=0; XCHECK(x236);
    x244=0; XCHECK(x244);
    x247=0; XCHECK(x247);
    x256=0; XCHECK(x256);
    x257=random2*x143; XCHECK(x257);
    x269=0; XCHECK(x269);
    x273=0; XCHECK(x273);
    x275=0; XCHECK(x275);
    x279=0; XCHECK(x279);
    result[10]=x232; XCHECK(result[10]);
    result[13]=x240; XCHECK(result[13]);
    result[16]=x245; XCHECK(result[16]);
    result[11]=x274; XCHECK(result[11]);
    result[14]=x277; XCHECK(result[14]);
    result[17]=x278; XCHECK(result[17]);
  };
  if (x9)
  {
    result[0]=1.00000000000000; XCHECK(result[0]);
    result[3]=x81; XCHECK(result[3]);
    result[6]=x81; XCHECK(result[6]);
    result[1]=x137; XCHECK(result[1]);
    result[4]=x160; XCHECK(result[4]);
    result[7]=x137; XCHECK(result[7]);
    result[2]=x137; XCHECK(result[2]);
    result[5]=x137; XCHECK(result[5]);
    result[8]=x160; XCHECK(result[8]);
  };
  if (x50)
  {
    result[0]=x16*(x22 - x24 + x26); XCHECK(result[0]);
    result[3]=x85; XCHECK(result[3]);
    result[6]=x85; XCHECK(result[6]);
  };
  if (x67)
  {
    result[0]=x60*x62*(x59*x61 - x64 - x65); XCHECK(result[0]);
    result[3]=x87*(-x64*x88 + x65*x88 + x93); XCHECK(result[3]);
    result[6]=x87*(-x106*x64 + x106*x65 + x93); XCHECK(result[6]);
  };
  if (x80)
  {
    result[0]=-x76 - x79 + 1; XCHECK(result[0]);
    result[3]=-x104 - x105 + x96; XCHECK(result[3]);
    result[6]=x107 - x113 - x114; XCHECK(result[6]);
    result[1]=x152*x157 - x157*x158 - x159*x76 + x159*x79; XCHECK(result[1]);
    result[4]=mu - x104*x159 + x105*x159 + x152*x214 - x158*x214; XCHECK(result[4]);
    result[7]=-x113*x159 + x114*x159 + x152*x231 - x158*x231; XCHECK(result[7]);
    result[2]=x152*x256 - x158*x256 - x257*x76 + x257*x79; XCHECK(result[2]);
    result[5]=-x104*x257 + x105*x257 + x152*x269 - x158*x269; XCHECK(result[5]);
    result[8]=mu - x113*x257 + x114*x257 + x152*x273 - x158*x273; XCHECK(result[8]);
  };
  if (x116)
  {
    result[9]=mu; XCHECK(result[9]);
    result[12]=0.0; XCHECK(result[12]);
    result[15]=0.0; XCHECK(result[15]);
    result[10]=x137; XCHECK(result[10]);
    result[13]=x237; XCHECK(result[13]);
    result[16]=x137; XCHECK(result[16]);
    result[11]=x137; XCHECK(result[11]);
    result[14]=x137; XCHECK(result[14]);
    result[17]=x237; XCHECK(result[17]);
  };
  if (x49)
  {
    result[9]=x16*(-x11*x19 - x11*x21 - x117 + x118 + 4.0e-10*x119); XCHECK(result[9]);
    result[12]=x129; XCHECK(result[12]);
    result[15]=x129; XCHECK(result[15]);
  };
  if (x126)
  {
    result[9]=mu - x124 - x125; XCHECK(result[9]);
    result[12]=-x132 - x133; XCHECK(result[12]);
    result[15]=-x135 - x136; XCHECK(result[15]);
    result[10]=-x124*x159 + x125*x159 + x152*x236 - x158*x236; XCHECK(result[10]);
    result[13]=-x132*x159 + x133*x159 + x152*x244 - x158*x244 + 1; XCHECK(result[13]);
    result[16]=-x135*x159 + x136*x159 + x248; XCHECK(result[16]);
    result[11]=-x124*x257 + x125*x257 + x152*x275 - x158*x275; XCHECK(result[11]);
    result[14]=-x132*x257 + x133*x257 + x248; XCHECK(result[14]);
    result[17]=-x135*x257 + x136*x257 + x152*x279 - x158*x279 + 1; XCHECK(result[17]);
  };
  if ((x50) && (x37))
  {
    result[1]=x141; XCHECK(result[1]);
    result[4]=x169; XCHECK(result[4]);
    result[7]=x216; XCHECK(result[7]);
    result[2]=x141; XCHECK(result[2]);
    result[5]=x216; XCHECK(result[5]);
    result[8]=x169; XCHECK(result[8]);
  };
  if ((x50) && (x36))
  {
    result[1]=x144*(-random1*x24 + x147); XCHECK(result[1]);
    result[4]=x170*(-random1*x83 + x171 + x172); XCHECK(result[4]);
    result[7]=x170*(-random1*x217 + x171); XCHECK(result[7]);
    result[2]=x144*(-random2*x24 + x251); XCHECK(result[2]);
    result[5]=x170*(-random2*x217 + x258); XCHECK(result[5]);
    result[8]=x170*(-random2*x83 + x172 + x258); XCHECK(result[8]);
  };
  if ((x67) && (x37))
  {
    result[1]=x87*(-x148*x64 + x148*x65); XCHECK(result[1]);
    result[4]=x177*(x183 - x184*x59 + x184*x61 + x189 + x190 - x192*x59 + x192*x61 + x194 + x195 - x196*x59 + x196*x61 - 4.0*x197 + 4.0*x198 - 2.0*x200 + 2.0*x201 + x204 + x205 - 2.0*x206 - 2.0*x207 - 4.0*x208 - 4.0*x209); XCHECK(result[4]);
    result[7]=x177*(x189 + x190 + x194 + x195 + x204 + x205 + x228); XCHECK(result[7]);
    result[2]=x87*(-x252*x64 + x252*x65); XCHECK(result[2]);
    result[5]=x177*(x228 + x260 + x261 + x263 + x264 + x266 + x267); XCHECK(result[5]);
    result[8]=x177*(x183 - 2.0*x197 + 2.0*x198 - 4.0*x200 + 4.0*x201 - 4.0*x206 - 4.0*x207 - 2.0*x208 - 2.0*x209 + x260 + x261 + x263 + x264 + x266 + x267 - x270*x59 + x270*x61 - x271*x59 + x271*x61 - x272*x59 + x272*x61); XCHECK(result[8]);
  };
  if ((x67) && (x36))
  {
    result[1]=x151*(-x149 + x150); XCHECK(result[1]);
    result[4]=x210*(-x149*x88 - x150*x88 + x211 + x212); XCHECK(result[4]);
    result[7]=x210*(-x106*x149 - x106*x150 + x211); XCHECK(result[7]);
    result[2]=x151*(-x253 + x254); XCHECK(result[2]);
    result[5]=x210*(-x253*x88 - x254*x88 + x268); XCHECK(result[5]);
    result[8]=x210*(-x106*x253 - x106*x254 + x212 + x268); XCHECK(result[8]);
  };
  if ((x234) && (x37))
  {
    result[10]=x233; XCHECK(result[10]);
    result[13]=x242; XCHECK(result[13]);
    result[16]=x246; XCHECK(result[16]);
    result[11]=x233; XCHECK(result[11]);
    result[14]=x246; XCHECK(result[14]);
    result[17]=x242; XCHECK(result[17]);
  };
  if ((x234) && (x36))
  {
    result[10]=x232; XCHECK(result[10]);
    result[13]=x240; XCHECK(result[13]);
    result[16]=x245; XCHECK(result[16]);
    result[11]=x274; XCHECK(result[11]);
    result[14]=x277; XCHECK(result[14]);
    result[17]=x278; XCHECK(result[17]);
  };
}

void frictionContact3D_FischerBurmeisterFunctionGenerated(
  double *reaction,
  double *velocity,
  double mu,
  double *rho,
  double *f,
  double *A,
  double *B)
{
  double result[21];

  assert(reaction);
  assert(velocity);
  assert(rho);

  SET3(reaction);
  SET3(velocity);
  SET3(rho);


  if (f && A && B)
  {

    frictionContact3D_FischerBurmeisterFABGenerated(
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
      frictionContact3D_FischerBurmeisterFGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3(result, f);
    }

    if (A && B)
    {
      frictionContact3D_FischerBurmeisterABGenerated(
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
