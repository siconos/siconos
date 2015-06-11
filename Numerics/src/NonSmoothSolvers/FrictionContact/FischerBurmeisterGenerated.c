#include <math.h>
#include <assert.h>
#include <op3x3.h>
#include <stdlib.h>

//#define DEBUG_MESSAGES 1
//#include <stdio.h>
#include <debug.h>
#include "FischerBurmeisterGenerated.h"

#define RESULT_CHECK(X)
#define VALUE_CHECK(X)

#define ZERO 1e-13
#define NOT_ZERO(x) fabs(x) > 0
#define IS_NOT_ZERO(x) fabs(x) > 0
#define IS_POSITIVE(x) 1

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
#define sqrt(x) (x < 0 ? 0 : sqrt(x))

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
double x1;
double x2;
double x3;
double x4;
double x8;
double x11;
double x12;
double x13;
double x30;
int x31;


x1=ut1*ut1;
x2=ut2*ut2;
x3=(assert(IS_POSITIVE(x1 + x2)), sqrt(x1 + x2));
x4=mu*x3 + un;
x8=mu*mu;
x11=x1*x8;
x12=x2*x8;
x13=x4*x4;
x30=(assert(IS_POSITIVE(x11 + x12 + x13)), sqrt(x11 + x12 + x13));
x31=x30 <= ZERO;

double x6;
double x7;
double x9;
double x10;
double x14;
double x15;
double x16;
double x17;
double x18;
double x19;
double x20;
double x21;
double x22;
double x23;
double x24;
double x27;
double x41;
int x42;

double x32;
double x33;
double x34;
double x35;
double x36;
double x37;
double x38;
double x39;
double x40;

x6=rt1*rt1;
x7=rt2*rt2;
x9=rn*rn;
x10=x8*x9;
x14=x10 + x11 + x12 + x13 + x6 + x7;
x15=mu*ut1;
x16=mu*rn*rt1 + x15*x4;
x17=x16*x16;
x18=mu*ut2;
x19=mu*rn*rt2 + x18*x4;
x20=x19*x19;
x21=x17 + x20;
x22=(assert(IS_POSITIVE(x21)), sqrt(x21));
x23=2*x22;
x24=x14 - x23;
x27=x14 + x23;
x41=fabs(x24);
x42=x22 <= ZERO || x27 <= ZERO || x41 <= ZERO;

int x54;

double x5;
double x43;
double x44;
double x45;
double x46;
double x47;
double x48;
double x49;
double x50;
double x51;
double x52;
double x53;

x54=x3 <= ZERO;

int x63;
int x64;
int x65;
int x66;
int x67;
int x68;

double x25;
double x28;
double x55;
double x56;
double x57;
double x58;
double x59;
double x60;
double x61;
double x62;

x63=x3 > ZERO;
x64=x30 > ZERO;
x65=x22 > ZERO;
x66=x27 > ZERO;
x67=x41 > ZERO;
x68=x63 && x64 && x65 && x66 && x67;

double x99;
int x100;


x43=x6 + x7;
x99=(assert(IS_POSITIVE(x10 + x43)), sqrt(x10 + x43));
x100=x99 <= ZERO;

int x107;
int x108;

double x103;
double x104;
double x105;
double x106;

x107=x99 > ZERO;
x108=x107 && x65 && x66 && x67;

int x121;

double x26;
double x29;
double x119;
double x120;

x121=x22 > 0;

int x124;

double x122;
double x123;

x124=x22 <= 0;

int x137;

double x133;
double x134;
double x135;
double x136;

x137=x121 && x63 && x64 && x65 && x66 && x67;

int x138;


x138=x124 && x63 && x64 && x65 && x66 && x67;

int x207;

double x205;
double x206;

x207=x107 && x121 && x65 && x66 && x67;

int x208;


x208=x107 && x124 && x65 && x66 && x67;

if (x31)
{

}
else if (x42)
{
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
    x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
    x36=0.5*x32;
    x37=0.5*x33;
    x38=-x36 - x37;
    x39=1.4142135623730951455*x32;
    x40=x32*x33;

}
else if (x54)
{
    x5=mu*rn;
    x44=2.0*(assert(IS_POSITIVE(x43)), sqrt(x43))*x5;
    x45=un*un + x10 + x6 + x7;
    x46=(assert(IS_POSITIVE(-x44 + x45)), sqrt(-x44 + x45));
    x47=1.0/(assert(IS_NOT_ZERO(x46)), x46);
    x48=(assert(IS_POSITIVE(x44 + x45)), sqrt(x44 + x45));
    x49=1.0/(assert(IS_NOT_ZERO(x48)), x48);
    x50=0.5*un;
    x51=x46*x50;
    x52=x48*x50;
    x53=x46*x48;

}
else if (x68)
{
    x25=(assert(IS_POSITIVE(x24)), sqrt(x24));
    x28=(assert(IS_POSITIVE(x27)), sqrt(x27));
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x56=x15*x16;
    x57=x18*x19;
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x59=x58*(x56 + x57);
    x60=x55*(x4 + x59);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x62=x61*(x4 - x59);

}
else if (x100)
{

}
else if (x108)
{
    x25=(assert(IS_POSITIVE(x24)), sqrt(x24));
    x28=(assert(IS_POSITIVE(x27)), sqrt(x27));
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x103=rn*x8;
    x104=x58*(mu*rt1*x16 + mu*rt2*x19);
    x105=x55*(x103 + x104);
    x106=x61*(x103 - x104);

}
else if (x121)
{
    x25=(assert(IS_POSITIVE(x24)), sqrt(x24));
    x26=0.5*x25;
    x28=(assert(IS_POSITIVE(x27)), sqrt(x27));
    x29=0.5*x28;
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x119=mu*ut1 + rt1;
    x120=x16*x58;

}
else if (x124)
{
    x25=(assert(IS_POSITIVE(x24)), sqrt(x24));
    x26=0.5*x25;
    x28=(assert(IS_POSITIVE(x27)), sqrt(x27));
    x29=0.5*x28;
    x119=mu*ut1 + rt1;
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x123=random1*x122;

}
else if (x137)
{
    x25=(assert(IS_POSITIVE(x24)), sqrt(x24));
    x26=0.5*x25;
    x28=(assert(IS_POSITIVE(x27)), sqrt(x27));
    x29=0.5*x28;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x56=x15*x16;
    x57=x18*x19;
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x59=x58*(x56 + x57);
    x60=x55*(x4 + x59);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x62=x61*(x4 - x59);
    x120=x16*x58;
    x133=-x56 - x57;
    x134=pow(x21, -3.0/2.0);
    x135=x134*x16;
    x136=x133*x135 + x15*x58;

}
else if (x138)
{
    x25=(assert(IS_POSITIVE(x24)), sqrt(x24));
    x28=(assert(IS_POSITIVE(x27)), sqrt(x27));
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x56=x15*x16;
    x57=x18*x19;
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x59=x58*(x56 + x57);
    x60=x55*(x4 + x59);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x62=x61*(x4 - x59);
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x123=random1*x122;

}
else if (x207)
{
    x25=(assert(IS_POSITIVE(x24)), sqrt(x24));
    x26=0.5*x25;
    x28=(assert(IS_POSITIVE(x27)), sqrt(x27));
    x29=0.5*x28;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x103=rn*x8;
    x104=x58*(mu*rt1*x16 + mu*rt2*x19);
    x105=x55*(x103 + x104);
    x106=x61*(x103 - x104);
    x120=x16*x58;
    x134=pow(x21, -3.0/2.0);
    x135=x134*x16;
    x205=-mu*rt1*x16 - mu*rt2*x19;
    x206=mu*rt1*x58 + x135*x205;

}
else if (x208)
{
    x25=(assert(IS_POSITIVE(x24)), sqrt(x24));
    x28=(assert(IS_POSITIVE(x27)), sqrt(x27));
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x103=rn*x8;
    x104=x58*(mu*rt1*x16 + mu*rt2*x19);
    x105=x55*(x103 + x104);
    x106=x61*(x103 - x104);
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x123=random1*x122;

}/* Assignment result[0, 0]=-x26 - x29 + x4 + x5 */
x5=mu*rn;
x25=(assert(IS_POSITIVE(x24)), sqrt(x24));
x26=0.5*x25;
x28=(assert(IS_POSITIVE(x27)), sqrt(x27));
x29=0.5*x28;
result[0] = -x26 - x29 + x4 + x5;
/* Assignment result[1, 0]=Piecewise((x119 + x120*x26 - x120*x29, x121), (x119 + x123*x26 - x123*x29, x124)) */
if (x121)
{
    DEBUG_PRINT("Case (x121) is True.\n");
    /* Assignment result[1, 0]=x119 + x120*x26 - x120*x29 */
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x119=mu*ut1 + rt1;
    x120=x16*x58;
    result[1] = x119 + x120*x26 - x120*x29;
}
else if (x124)
{
    DEBUG_PRINT("Case (x124) is True.\n");
    /* Assignment result[1, 0]=x119 + x123*x26 - x123*x29 */
    x119=mu*ut1 + rt1;
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x123=random1*x122;
    result[1] = x119 + x123*x26 - x123*x29;
}
/* Assignment result[2, 0]=Piecewise((x218 + x219*x26 - x219*x29, x121), (x218 + x220*x26 - x220*x29, x124)) */
if (x121)
{
    DEBUG_PRINT("Case (x121) is True.\n");
    /* Assignment result[2, 0]=x218 + x219*x26 - x219*x29 */
    double x218;
    double x219;
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x218=mu*ut2 + rt2;
    x219=x19*x58;
    result[2] = x218 + x219*x26 - x219*x29;
}
else if (x124)
{
    DEBUG_PRINT("Case (x124) is True.\n");
    /* Assignment result[2, 0]=x218 + x220*x26 - x220*x29 */
    double x218;
    double x220;
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x218=mu*ut2 + rt2;
    x220=random2*x122;
    result[2] = x218 + x220*x26 - x220*x29;
}
/* Assignment result[0, 1]=Piecewise((1.00000000000000, x31), (x35*(-mu*x39 + x38 + x40), x42), (x47*x49*(-x51 - x52 + x53), x54), (-x60 - x62 + 1, x68)) */
if (x31)
{
    DEBUG_PRINT("Case (x31) is True.\n");
    /* Assignment result[0, 1]=1.00000000000000 */
    result[3] = 1.0;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[0, 1]=x35*(-mu*x39 + x38 + x40) */
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
    x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
    x36=0.5*x32;
    x37=0.5*x33;
    x38=-x36 - x37;
    x39=1.4142135623730951455*x32;
    x40=x32*x33;
    result[3] = x35*(-mu*x39 + x38 + x40);
}
else if (x54)
{
    DEBUG_PRINT("Case (x54) is True.\n");
    /* Assignment result[0, 1]=x47*x49*(-x51 - x52 + x53) */
    x44=2.0*(assert(IS_POSITIVE(x43)), sqrt(x43))*x5;
    x45=un*un + x10 + x6 + x7;
    x46=(assert(IS_POSITIVE(-x44 + x45)), sqrt(-x44 + x45));
    x47=1.0/(assert(IS_NOT_ZERO(x46)), x46);
    x48=(assert(IS_POSITIVE(x44 + x45)), sqrt(x44 + x45));
    x49=1.0/(assert(IS_NOT_ZERO(x48)), x48);
    x50=0.5*un;
    x51=x46*x50;
    x52=x48*x50;
    x53=x46*x48;
    result[3] = x47*x49*(-x51 - x52 + x53);
}
else if (x68)
{
    DEBUG_PRINT("Case (x68) is True.\n");
    /* Assignment result[0, 1]=-x60 - x62 + 1 */
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x56=x15*x16;
    x57=x18*x19;
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x59=x58*(x56 + x57);
    x60=x55*(x4 + x59);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x62=x61*(x4 - x59);
    result[3] = -x60 - x62 + 1;
}
/* Assignment result[1, 1]=Piecewise((0.0, x31), (x125, x42), (x127*(-x131*x51 + x131*x52 - x132*x51 + x132*x52), x54), (-x120*x60 + x120*x62 + x136*x26 - x136*x29, x137), (-x123*x60 + x123*x62, x138)) */
if (x31)
{
    DEBUG_PRINT("Case (x31) is True.\n");
    /* Assignment result[1, 1]=0.0 */
    result[4] = 0.0;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[1, 1]=x125 */
    double x101;
    double x109;
    double x110;
    double x125;
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
    x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
    x101=mu*x32;
    x109=0.35355339059327378637*x32;
    x110=0.35355339059327378637*x33;
    x125=x35*(-1.0*x101 - x109 + x110);
    result[4] = x125;
}
else if (x54)
{
    DEBUG_PRINT("Case (x54) is True.\n");
    /* Assignment result[1, 1]=x127*(-x131*x51 + x131*x52 - x132*x51 + x132*x52) */
    double x72;
    double x126;
    double x127;
    double x128;
    double x129;
    double x130;
    double x131;
    double x132;
    x44=2.0*(assert(IS_POSITIVE(x43)), sqrt(x43))*x5;
    x45=un*un + x10 + x6 + x7;
    x46=(assert(IS_POSITIVE(-x44 + x45)), sqrt(-x44 + x45));
    x47=1.0/(assert(IS_NOT_ZERO(x46)), x46);
    x48=(assert(IS_POSITIVE(x44 + x45)), sqrt(x44 + x45));
    x49=1.0/(assert(IS_NOT_ZERO(x48)), x48);
    x50=0.5*un;
    x51=x46*x50;
    x52=x48*x50;
    x72=x10*x6 + x10*x7;
    x126=pow(x72, 3.0/2.0);
    x127=(assert(IS_NOT_ZERO(x126)), x47*x49/x126);
    x128=pow(mu, 3);
    x129=pow(rn, 3);
    x130=pow(rt1, 3);
    x131=x128*x129*x130;
    x132=rt1*x128*x129*x7;
    result[4] = x127*(-x131*x51 + x131*x52 - x132*x51 + x132*x52);
}
else if (x137)
{
    DEBUG_PRINT("Case (x137) is True.\n");
    /* Assignment result[1, 1]=-x120*x60 + x120*x62 + x136*x26 - x136*x29 */
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x56=x15*x16;
    x57=x18*x19;
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x59=x58*(x56 + x57);
    x60=x55*(x4 + x59);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x62=x61*(x4 - x59);
    x120=x16*x58;
    x133=-x56 - x57;
    x134=pow(x21, -3.0/2.0);
    x135=x134*x16;
    x136=x133*x135 + x15*x58;
    result[4] = -x120*x60 + x120*x62 + x136*x26 - x136*x29;
}
else if (x138)
{
    DEBUG_PRINT("Case (x138) is True.\n");
    /* Assignment result[1, 1]=-x123*x60 + x123*x62 */
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x56=x15*x16;
    x57=x18*x19;
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x59=x58*(x56 + x57);
    x60=x55*(x4 + x59);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x62=x61*(x4 - x59);
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x123=random1*x122;
    result[4] = -x123*x60 + x123*x62;
}
/* Assignment result[2, 1]=Piecewise((0.0, x31), (x125, x42), (x127*(-x221*x51 + x221*x52 - x222*x51 + x222*x52), x54), (-x219*x60 + x219*x62 + x224*x26 - x224*x29, x137), (-x220*x60 + x220*x62, x138)) */
if (x31)
{
    DEBUG_PRINT("Case (x31) is True.\n");
    /* Assignment result[2, 1]=0.0 */
    result[5] = 0.0;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[2, 1]=x125 */
    double x101;
    double x109;
    double x110;
    double x125;
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
    x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
    x101=mu*x32;
    x109=0.35355339059327378637*x32;
    x110=0.35355339059327378637*x33;
    x125=x35*(-1.0*x101 - x109 + x110);
    result[5] = x125;
}
else if (x54)
{
    DEBUG_PRINT("Case (x54) is True.\n");
    /* Assignment result[2, 1]=x127*(-x221*x51 + x221*x52 - x222*x51 + x222*x52) */
    double x72;
    double x126;
    double x127;
    double x128;
    double x129;
    double x193;
    double x221;
    double x222;
    x44=2.0*(assert(IS_POSITIVE(x43)), sqrt(x43))*x5;
    x45=un*un + x10 + x6 + x7;
    x46=(assert(IS_POSITIVE(-x44 + x45)), sqrt(-x44 + x45));
    x47=1.0/(assert(IS_NOT_ZERO(x46)), x46);
    x48=(assert(IS_POSITIVE(x44 + x45)), sqrt(x44 + x45));
    x49=1.0/(assert(IS_NOT_ZERO(x48)), x48);
    x50=0.5*un;
    x51=x46*x50;
    x52=x48*x50;
    x72=x10*x6 + x10*x7;
    x126=pow(x72, 3.0/2.0);
    x127=(assert(IS_NOT_ZERO(x126)), x47*x49/x126);
    x128=pow(mu, 3);
    x129=pow(rn, 3);
    x193=pow(rt2, 3);
    x221=x128*x129*x193;
    x222=rt2*x128*x129*x6;
    result[5] = x127*(-x221*x51 + x221*x52 - x222*x51 + x222*x52);
}
else if (x137)
{
    DEBUG_PRINT("Case (x137) is True.\n");
    /* Assignment result[2, 1]=-x219*x60 + x219*x62 + x224*x26 - x224*x29 */
    double x219;
    double x223;
    double x224;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x56=x15*x16;
    x57=x18*x19;
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x59=x58*(x56 + x57);
    x60=x55*(x4 + x59);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x62=x61*(x4 - x59);
    x133=-x56 - x57;
    x134=pow(x21, -3.0/2.0);
    x219=x19*x58;
    x223=x134*x19;
    x224=x133*x223 + x18*x58;
    result[5] = -x219*x60 + x219*x62 + x224*x26 - x224*x29;
}
else if (x138)
{
    DEBUG_PRINT("Case (x138) is True.\n");
    /* Assignment result[2, 1]=-x220*x60 + x220*x62 */
    double x220;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x56=x15*x16;
    x57=x18*x19;
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x59=x58*(x56 + x57);
    x60=x55*(x4 + x59);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x62=x61*(x4 - x59);
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x220=random2*x122;
    result[5] = -x220*x60 + x220*x62;
}
/* Assignment result[0, 2]=Piecewise((x69, x31), (x71, x42), (x74*(-x51*x75 + x52*x75 + x77), x54), (x79 - x88 - x89, x68)) */
if (x31)
{
    DEBUG_PRINT("Case (x31) is True.\n");
    /* Assignment result[0, 2]=x69 */
    double x69;
    x69=0.70710678118654757274*mu;
    result[6] = x69;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[0, 2]=x71 */
    double x69;
    double x70;
    double x71;
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
    x69=0.70710678118654757274*mu;
    x70=x33*x69;
    x71=x34*(-x69 + x70 - 2.0*x8);
    result[6] = x71;
}
else if (x54)
{
    DEBUG_PRINT("Case (x54) is True.\n");
    /* Assignment result[0, 2]=x74*(-x51*x75 + x52*x75 + x77) */
    double x72;
    double x73;
    double x74;
    double x75;
    double x76;
    double x77;
    x44=2.0*(assert(IS_POSITIVE(x43)), sqrt(x43))*x5;
    x45=un*un + x10 + x6 + x7;
    x46=(assert(IS_POSITIVE(-x44 + x45)), sqrt(-x44 + x45));
    x47=1.0/(assert(IS_NOT_ZERO(x46)), x46);
    x48=(assert(IS_POSITIVE(x44 + x45)), sqrt(x44 + x45));
    x49=1.0/(assert(IS_NOT_ZERO(x48)), x48);
    x50=0.5*un;
    x51=x46*x50;
    x52=x48*x50;
    x53=x46*x48;
    x72=x10*x6 + x10*x7;
    x73=(assert(IS_POSITIVE(x72)), sqrt(x72));
    x74=(assert(IS_NOT_ZERO(x73)), x47*x49/x73);
    x75=rn*rt1*x8;
    x76=0.35355339059327378637*mu*un*x73;
    x77=0.70710678118654757274*mu*x53*x73 - x46*x76 - x48*x76;
    result[6] = x74*(-x51*x75 + x52*x75 + x77);
}
else if (x68)
{
    DEBUG_PRINT("Case (x68) is True.\n");
    /* Assignment result[0, 2]=x79 - x88 - x89 */
    double x78;
    double x79;
    double x80;
    double x81;
    double x82;
    double x83;
    double x84;
    double x85;
    double x86;
    double x87;
    double x88;
    double x89;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x78=1.0/(assert(IS_NOT_ZERO(x3)), x3);
    x79=x15*x78;
    x80=ut1*x8 + x4*x79;
    x81=ut1*ut2*x78*x8;
    x82=x19*x81;
    x83=mu*x4;
    x84=2*x83;
    x85=x1*x78*x8;
    x86=(1.0/2.0)*x16*(x84 + 2*x85);
    x87=x58*(x82 + x86);
    x88=x55*(x80 + x87);
    x89=x61*(x80 - x87);
    result[6] = x79 - x88 - x89;
}
/* Assignment result[1, 2]=Piecewise((mu, x31), (x151, x42), (x155*(-x149*x178 + x149*x179 + x161 - x162*x46 + x162*x48 + x166 + x167 - x169*x46 + x169*x48 + x171 + x172 - x173*x46 + x173*x48 - x174*x175 + x174*x176 + x181 + x182 - 2.0*x183 - 2.0*x184 - 4.0*x185 - 4.0*x186), x54), (mu - x120*x88 + x120*x89 + x188*x26 - x188*x29, x137), (mu - x123*x88 + x123*x89, x138)) */
if (x31)
{
    DEBUG_PRINT("Case (x31) is True.\n");
    /* Assignment result[1, 2]=mu */
    result[7] = mu;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[1, 2]=x151 */
    double x101;
    double x102;
    double x128;
    double x139;
    double x140;
    double x141;
    double x142;
    double x143;
    double x144;
    double x145;
    double x146;
    double x147;
    double x148;
    double x149;
    double x150;
    double x151;
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x40=x32*x33;
    x101=mu*x32;
    x102=x101*x33;
    x128=pow(mu, 3);
    x139=11.313708498984761164*x102 + 16.0*x40;
    x140=1.0/(assert(IS_NOT_ZERO(x139)), x139);
    x141=4.2426406871192847703*x32;
    x142=4.2426406871192847703*x33;
    x143=mu*x33;
    x144=2.0*x143;
    x145=x32*x8;
    x146=x128*x32;
    x147=x142*x8;
    x148=-x147;
    x149=2.0*x128;
    x150=x149*x33;
    x151=x140*(-26.0*x101 + 16.0*x102 - x141 + x142 + x144 + 11.313708498984761164*x145*x33 - 57.982756057296896302*x145 - 34.0*x146 + x148 + x150);
    result[7] = x151;
}
else if (x54)
{
    DEBUG_PRINT("Case (x54) is True.\n");
    /* Assignment result[1, 2]=x155*(-x149*x178 + x149*x179 + x161 - x162*x46 + x162*x48 + x166 + x167 - x169*x46 + x169*x48 + x171 + x172 - x173*x46 + x173*x48 - x174*x175 + x174*x176 + x181 + x182 - 2.0*x183 - 2.0*x184 - 4.0*x185 - 4.0*x186) */
    double x72;
    double x126;
    double x128;
    double x129;
    double x130;
    double x149;
    double x152;
    double x153;
    double x154;
    double x155;
    double x156;
    double x157;
    double x158;
    double x159;
    double x160;
    double x161;
    double x162;
    double x163;
    double x164;
    double x165;
    double x166;
    double x167;
    double x168;
    double x169;
    double x170;
    double x171;
    double x172;
    double x173;
    double x174;
    double x175;
    double x176;
    double x177;
    double x178;
    double x179;
    double x180;
    double x181;
    double x182;
    double x183;
    double x184;
    double x185;
    double x186;
    x44=2.0*(assert(IS_POSITIVE(x43)), sqrt(x43))*x5;
    x45=un*un + x10 + x6 + x7;
    x46=(assert(IS_POSITIVE(-x44 + x45)), sqrt(-x44 + x45));
    x48=(assert(IS_POSITIVE(x44 + x45)), sqrt(x44 + x45));
    x72=x10*x6 + x10*x7;
    x126=pow(x72, 3.0/2.0);
    x128=pow(mu, 3);
    x129=pow(rn, 3);
    x130=pow(rt1, 3);
    x149=2.0*x128;
    x152=4.0*x126*x46*x48;
    x153=x152*x6;
    x154=x152*x7;
    x155=1.0/(assert(IS_NOT_ZERO((x153 + x154))), (x153 + x154));
    x156=pow(mu, 5);
    x157=pow(rn, 4);
    x158=2.0*un*x156*x157*x6*x7;
    x159=pow(un, 3);
    x160=2.0*x128*x159*x6*x7*x9;
    x161=mu*x153 + mu*x154 - x158*x46 + x158*x48 - x160*x46 + x160*x48;
    x162=2.0*pow(rt2, 6)*un*x128*x9;
    x163=pow(mu, 4);
    x164=pow(rt1, 5);
    x165=1.4142135623730951455*un*x129*x163*x164;
    x166=-x165*x46;
    x167=x165*x48;
    x168=pow(rt2, 4);
    x169=2.0*un*x156*x157*x168;
    x170=1.4142135623730951455*rt1*un*x129*x163*x168;
    x171=-x170*x46;
    x172=x170*x48;
    x173=2.0*x128*x159*x168*x9;
    x174=4.0*x128;
    x175=un*x168*x46*x6*x9;
    x176=un*x168*x48*x6*x9;
    x177=pow(rt1, 4);
    x178=un*x177*x46*x7*x9;
    x179=un*x177*x48*x7*x9;
    x180=2.8284271247461902909*un*x129*x130*x163*x7;
    x181=-x180*x46;
    x182=x180*x48;
    x183=mu*un*x126*x46*x6;
    x184=mu*un*x126*x48*x6;
    x185=mu*un*x126*x46*x7;
    x186=mu*un*x126*x48*x7;
    result[7] = x155*(-x149*x178 + x149*x179 + x161 - x162*x46 + x162*x48 + x166 + x167 - x169*x46 + x169*x48 + x171 + x172 - x173*x46 + x173*x48 - x174*x175 + x174*x176 + x181 + x182 - 2.0*x183 - 2.0*x184 - 4.0*x185 - 4.0*x186);
}
else if (x137)
{
    DEBUG_PRINT("Case (x137) is True.\n");
    /* Assignment result[1, 2]=mu - x120*x88 + x120*x89 + x188*x26 - x188*x29 */
    double x78;
    double x79;
    double x80;
    double x81;
    double x82;
    double x83;
    double x84;
    double x85;
    double x86;
    double x87;
    double x88;
    double x89;
    double x187;
    double x188;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x78=1.0/(assert(IS_NOT_ZERO(x3)), x3);
    x79=x15*x78;
    x80=ut1*x8 + x4*x79;
    x81=ut1*ut2*x78*x8;
    x82=x19*x81;
    x83=mu*x4;
    x84=2*x83;
    x85=x1*x78*x8;
    x86=(1.0/2.0)*x16*(x84 + 2*x85);
    x87=x58*(x82 + x86);
    x88=x55*(x80 + x87);
    x89=x61*(x80 - x87);
    x120=x16*x58;
    x134=pow(x21, -3.0/2.0);
    x135=x134*x16;
    x187=-x82 - x86;
    x188=x135*x187 + x58*(x83 + x85);
    result[7] = mu - x120*x88 + x120*x89 + x188*x26 - x188*x29;
}
else if (x138)
{
    DEBUG_PRINT("Case (x138) is True.\n");
    /* Assignment result[1, 2]=mu - x123*x88 + x123*x89 */
    double x78;
    double x79;
    double x80;
    double x81;
    double x82;
    double x83;
    double x84;
    double x85;
    double x86;
    double x87;
    double x88;
    double x89;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x78=1.0/(assert(IS_NOT_ZERO(x3)), x3);
    x79=x15*x78;
    x80=ut1*x8 + x4*x79;
    x81=ut1*ut2*x78*x8;
    x82=x19*x81;
    x83=mu*x4;
    x84=2*x83;
    x85=x1*x78*x8;
    x86=(1.0/2.0)*x16*(x84 + 2*x85);
    x87=x58*(x82 + x86);
    x88=x55*(x80 + x87);
    x89=x61*(x80 - x87);
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x123=random1*x122;
    result[7] = mu - x123*x88 + x123*x89;
}
/* Assignment result[2, 2]=Piecewise((0.0, x31), (x189, x42), (x155*(x200 + x226 + x227 + x229 + x230 + x232 + x233), x54), (-x219*x88 + x219*x89 + x234*x26 - x234*x29, x137), (-x220*x88 + x220*x89, x138)) */
if (x31)
{
    DEBUG_PRINT("Case (x31) is True.\n");
    /* Assignment result[2, 2]=0.0 */
    result[8] = 0.0;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[2, 2]=x189 */
    double x101;
    double x102;
    double x128;
    double x139;
    double x140;
    double x141;
    double x142;
    double x143;
    double x144;
    double x146;
    double x147;
    double x149;
    double x150;
    double x189;
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x39=1.4142135623730951455*x32;
    x40=x32*x33;
    x101=mu*x32;
    x102=x101*x33;
    x128=pow(mu, 3);
    x139=11.313708498984761164*x102 + 16.0*x40;
    x140=1.0/(assert(IS_NOT_ZERO(x139)), x139);
    x141=4.2426406871192847703*x32;
    x142=4.2426406871192847703*x33;
    x143=mu*x33;
    x144=2.0*x143;
    x146=x128*x32;
    x147=x142*x8;
    x149=2.0*x128;
    x150=x149*x33;
    x189=x140*(10.0*x101 + x141 - x142 - x144 + 2.0*x146 + x147 - x150 + x39*x8);
    result[8] = x189;
}
else if (x54)
{
    DEBUG_PRINT("Case (x54) is True.\n");
    /* Assignment result[2, 2]=x155*(x200 + x226 + x227 + x229 + x230 + x232 + x233) */
    double x72;
    double x126;
    double x128;
    double x129;
    double x130;
    double x152;
    double x153;
    double x154;
    double x155;
    double x156;
    double x157;
    double x159;
    double x163;
    double x164;
    double x177;
    double x190;
    double x191;
    double x192;
    double x193;
    double x194;
    double x195;
    double x196;
    double x197;
    double x198;
    double x199;
    double x200;
    double x225;
    double x226;
    double x227;
    double x228;
    double x229;
    double x230;
    double x231;
    double x232;
    double x233;
    x44=2.0*(assert(IS_POSITIVE(x43)), sqrt(x43))*x5;
    x45=un*un + x10 + x6 + x7;
    x46=(assert(IS_POSITIVE(-x44 + x45)), sqrt(-x44 + x45));
    x48=(assert(IS_POSITIVE(x44 + x45)), sqrt(x44 + x45));
    x72=x10*x6 + x10*x7;
    x126=pow(x72, 3.0/2.0);
    x128=pow(mu, 3);
    x129=pow(rn, 3);
    x130=pow(rt1, 3);
    x152=4.0*x126*x46*x48;
    x153=x152*x6;
    x154=x152*x7;
    x155=1.0/(assert(IS_NOT_ZERO((x153 + x154))), (x153 + x154));
    x156=pow(mu, 5);
    x157=pow(rn, 4);
    x159=pow(un, 3);
    x163=pow(mu, 4);
    x164=pow(rt1, 5);
    x177=pow(rt1, 4);
    x190=pow(rt2, 5);
    x191=2.0*rt1*un*x128*x190*x9;
    x192=2.0*rt2*un*x128*x164*x9;
    x193=pow(rt2, 3);
    x194=2.0*rt1*un*x156*x157*x193;
    x195=2.0*rt2*un*x130*x156*x157;
    x196=2.0*rt1*x128*x159*x193*x9;
    x197=2.0*rt2*x128*x130*x159*x9;
    x198=4.0*un*x128*x130*x193*x9;
    x199=2.0*mu*rt1*rt2*un*x126;
    x200=x191*x46 - x191*x48 + x192*x46 - x192*x48 + x194*x46 - x194*x48 + x195*x46 - x195*x48 + x196*x46 - x196*x48 + x197*x46 - x197*x48 + x198*x46 - x198*x48 + x199*x46 + x199*x48;
    x225=1.4142135623730951455*un*x129*x163*x190;
    x226=-x225*x46;
    x227=x225*x48;
    x228=1.4142135623730951455*rt2*un*x129*x163*x177;
    x229=-x228*x46;
    x230=x228*x48;
    x231=2.8284271247461902909*un*x129*x163*x193*x6;
    x232=-x231*x46;
    x233=x231*x48;
    result[8] = x155*(x200 + x226 + x227 + x229 + x230 + x232 + x233);
}
else if (x137)
{
    DEBUG_PRINT("Case (x137) is True.\n");
    /* Assignment result[2, 2]=-x219*x88 + x219*x89 + x234*x26 - x234*x29 */
    double x78;
    double x79;
    double x80;
    double x81;
    double x82;
    double x83;
    double x84;
    double x85;
    double x86;
    double x87;
    double x88;
    double x89;
    double x187;
    double x201;
    double x219;
    double x223;
    double x234;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x78=1.0/(assert(IS_NOT_ZERO(x3)), x3);
    x79=x15*x78;
    x80=ut1*x8 + x4*x79;
    x81=ut1*ut2*x78*x8;
    x82=x19*x81;
    x83=mu*x4;
    x84=2*x83;
    x85=x1*x78*x8;
    x86=(1.0/2.0)*x16*(x84 + 2*x85);
    x87=x58*(x82 + x86);
    x88=x55*(x80 + x87);
    x89=x61*(x80 - x87);
    x134=pow(x21, -3.0/2.0);
    x187=-x82 - x86;
    x201=x58*x81;
    x219=x19*x58;
    x223=x134*x19;
    x234=x187*x223 + x201;
    result[8] = -x219*x88 + x219*x89 + x234*x26 - x234*x29;
}
else if (x138)
{
    DEBUG_PRINT("Case (x138) is True.\n");
    /* Assignment result[2, 2]=-x220*x88 + x220*x89 */
    double x78;
    double x79;
    double x80;
    double x81;
    double x82;
    double x83;
    double x84;
    double x85;
    double x86;
    double x87;
    double x88;
    double x89;
    double x220;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x78=1.0/(assert(IS_NOT_ZERO(x3)), x3);
    x79=x15*x78;
    x80=ut1*x8 + x4*x79;
    x81=ut1*ut2*x78*x8;
    x82=x19*x81;
    x83=mu*x4;
    x84=2*x83;
    x85=x1*x78*x8;
    x86=(1.0/2.0)*x16*(x84 + 2*x85);
    x87=x58*(x82 + x86);
    x88=x55*(x80 + x87);
    x89=x61*(x80 - x87);
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x220=random2*x122;
    result[8] = -x220*x88 + x220*x89;
}
/* Assignment result[0, 3]=Piecewise((x69, x31), (x71, x42), (x74*(-x51*x90 + x52*x90 + x77), x54), (x91 - x97 - x98, x68)) */
if (x31)
{
    DEBUG_PRINT("Case (x31) is True.\n");
    /* Assignment result[0, 3]=x69 */
    double x69;
    x69=0.70710678118654757274*mu;
    result[9] = x69;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[0, 3]=x71 */
    double x69;
    double x70;
    double x71;
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
    x69=0.70710678118654757274*mu;
    x70=x33*x69;
    x71=x34*(-x69 + x70 - 2.0*x8);
    result[9] = x71;
}
else if (x54)
{
    DEBUG_PRINT("Case (x54) is True.\n");
    /* Assignment result[0, 3]=x74*(-x51*x90 + x52*x90 + x77) */
    double x72;
    double x73;
    double x74;
    double x76;
    double x77;
    double x90;
    x44=2.0*(assert(IS_POSITIVE(x43)), sqrt(x43))*x5;
    x45=un*un + x10 + x6 + x7;
    x46=(assert(IS_POSITIVE(-x44 + x45)), sqrt(-x44 + x45));
    x47=1.0/(assert(IS_NOT_ZERO(x46)), x46);
    x48=(assert(IS_POSITIVE(x44 + x45)), sqrt(x44 + x45));
    x49=1.0/(assert(IS_NOT_ZERO(x48)), x48);
    x50=0.5*un;
    x51=x46*x50;
    x52=x48*x50;
    x53=x46*x48;
    x72=x10*x6 + x10*x7;
    x73=(assert(IS_POSITIVE(x72)), sqrt(x72));
    x74=(assert(IS_NOT_ZERO(x73)), x47*x49/x73);
    x76=0.35355339059327378637*mu*un*x73;
    x77=0.70710678118654757274*mu*x53*x73 - x46*x76 - x48*x76;
    x90=rn*rt2*x8;
    result[9] = x74*(-x51*x90 + x52*x90 + x77);
}
else if (x68)
{
    DEBUG_PRINT("Case (x68) is True.\n");
    /* Assignment result[0, 3]=x91 - x97 - x98 */
    double x78;
    double x81;
    double x83;
    double x84;
    double x91;
    double x92;
    double x93;
    double x94;
    double x95;
    double x96;
    double x97;
    double x98;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x78=1.0/(assert(IS_NOT_ZERO(x3)), x3);
    x81=ut1*ut2*x78*x8;
    x83=mu*x4;
    x84=2*x83;
    x91=x18*x78;
    x92=ut2*x8 + x4*x91;
    x93=x16*x81;
    x94=x2*x78*x8;
    x95=(1.0/2.0)*x19*(x84 + 2*x94);
    x96=x58*(x93 + x95);
    x97=x55*(x92 + x96);
    x98=x61*(x92 - x96);
    result[9] = x91 - x97 - x98;
}
/* Assignment result[1, 3]=Piecewise((0.0, x31), (x189, x42), (x155*(x166 + x167 + x171 + x172 + x181 + x182 + x200), x54), (-x120*x97 + x120*x98 + x203*x26 - x203*x29, x137), (-x123*x97 + x123*x98, x138)) */
if (x31)
{
    DEBUG_PRINT("Case (x31) is True.\n");
    /* Assignment result[1, 3]=0.0 */
    result[10] = 0.0;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[1, 3]=x189 */
    double x101;
    double x102;
    double x128;
    double x139;
    double x140;
    double x141;
    double x142;
    double x143;
    double x144;
    double x146;
    double x147;
    double x149;
    double x150;
    double x189;
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x39=1.4142135623730951455*x32;
    x40=x32*x33;
    x101=mu*x32;
    x102=x101*x33;
    x128=pow(mu, 3);
    x139=11.313708498984761164*x102 + 16.0*x40;
    x140=1.0/(assert(IS_NOT_ZERO(x139)), x139);
    x141=4.2426406871192847703*x32;
    x142=4.2426406871192847703*x33;
    x143=mu*x33;
    x144=2.0*x143;
    x146=x128*x32;
    x147=x142*x8;
    x149=2.0*x128;
    x150=x149*x33;
    x189=x140*(10.0*x101 + x141 - x142 - x144 + 2.0*x146 + x147 - x150 + x39*x8);
    result[10] = x189;
}
else if (x54)
{
    DEBUG_PRINT("Case (x54) is True.\n");
    /* Assignment result[1, 3]=x155*(x166 + x167 + x171 + x172 + x181 + x182 + x200) */
    double x72;
    double x126;
    double x128;
    double x129;
    double x130;
    double x152;
    double x153;
    double x154;
    double x155;
    double x156;
    double x157;
    double x159;
    double x163;
    double x164;
    double x165;
    double x166;
    double x167;
    double x168;
    double x170;
    double x171;
    double x172;
    double x180;
    double x181;
    double x182;
    double x190;
    double x191;
    double x192;
    double x193;
    double x194;
    double x195;
    double x196;
    double x197;
    double x198;
    double x199;
    double x200;
    x44=2.0*(assert(IS_POSITIVE(x43)), sqrt(x43))*x5;
    x45=un*un + x10 + x6 + x7;
    x46=(assert(IS_POSITIVE(-x44 + x45)), sqrt(-x44 + x45));
    x48=(assert(IS_POSITIVE(x44 + x45)), sqrt(x44 + x45));
    x72=x10*x6 + x10*x7;
    x126=pow(x72, 3.0/2.0);
    x128=pow(mu, 3);
    x129=pow(rn, 3);
    x130=pow(rt1, 3);
    x152=4.0*x126*x46*x48;
    x153=x152*x6;
    x154=x152*x7;
    x155=1.0/(assert(IS_NOT_ZERO((x153 + x154))), (x153 + x154));
    x156=pow(mu, 5);
    x157=pow(rn, 4);
    x159=pow(un, 3);
    x163=pow(mu, 4);
    x164=pow(rt1, 5);
    x165=1.4142135623730951455*un*x129*x163*x164;
    x166=-x165*x46;
    x167=x165*x48;
    x168=pow(rt2, 4);
    x170=1.4142135623730951455*rt1*un*x129*x163*x168;
    x171=-x170*x46;
    x172=x170*x48;
    x180=2.8284271247461902909*un*x129*x130*x163*x7;
    x181=-x180*x46;
    x182=x180*x48;
    x190=pow(rt2, 5);
    x191=2.0*rt1*un*x128*x190*x9;
    x192=2.0*rt2*un*x128*x164*x9;
    x193=pow(rt2, 3);
    x194=2.0*rt1*un*x156*x157*x193;
    x195=2.0*rt2*un*x130*x156*x157;
    x196=2.0*rt1*x128*x159*x193*x9;
    x197=2.0*rt2*x128*x130*x159*x9;
    x198=4.0*un*x128*x130*x193*x9;
    x199=2.0*mu*rt1*rt2*un*x126;
    x200=x191*x46 - x191*x48 + x192*x46 - x192*x48 + x194*x46 - x194*x48 + x195*x46 - x195*x48 + x196*x46 - x196*x48 + x197*x46 - x197*x48 + x198*x46 - x198*x48 + x199*x46 + x199*x48;
    result[10] = x155*(x166 + x167 + x171 + x172 + x181 + x182 + x200);
}
else if (x137)
{
    DEBUG_PRINT("Case (x137) is True.\n");
    /* Assignment result[1, 3]=-x120*x97 + x120*x98 + x203*x26 - x203*x29 */
    double x78;
    double x81;
    double x83;
    double x84;
    double x91;
    double x92;
    double x93;
    double x94;
    double x95;
    double x96;
    double x97;
    double x98;
    double x201;
    double x202;
    double x203;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x78=1.0/(assert(IS_NOT_ZERO(x3)), x3);
    x81=ut1*ut2*x78*x8;
    x83=mu*x4;
    x84=2*x83;
    x91=x18*x78;
    x92=ut2*x8 + x4*x91;
    x93=x16*x81;
    x94=x2*x78*x8;
    x95=(1.0/2.0)*x19*(x84 + 2*x94);
    x96=x58*(x93 + x95);
    x97=x55*(x92 + x96);
    x98=x61*(x92 - x96);
    x120=x16*x58;
    x134=pow(x21, -3.0/2.0);
    x135=x134*x16;
    x201=x58*x81;
    x202=-x93 - x95;
    x203=x135*x202 + x201;
    result[10] = -x120*x97 + x120*x98 + x203*x26 - x203*x29;
}
else if (x138)
{
    DEBUG_PRINT("Case (x138) is True.\n");
    /* Assignment result[1, 3]=-x123*x97 + x123*x98 */
    double x78;
    double x81;
    double x83;
    double x84;
    double x91;
    double x92;
    double x93;
    double x94;
    double x95;
    double x96;
    double x97;
    double x98;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x78=1.0/(assert(IS_NOT_ZERO(x3)), x3);
    x81=ut1*ut2*x78*x8;
    x83=mu*x4;
    x84=2*x83;
    x91=x18*x78;
    x92=ut2*x8 + x4*x91;
    x93=x16*x81;
    x94=x2*x78*x8;
    x95=(1.0/2.0)*x19*(x84 + 2*x94);
    x96=x58*(x93 + x95);
    x97=x55*(x92 + x96);
    x98=x61*(x92 - x96);
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x123=random1*x122;
    result[10] = -x123*x97 + x123*x98;
}
/* Assignment result[2, 3]=Piecewise((mu, x31), (x151, x42), (x155*(-x149*x175 + x149*x176 + x161 - x174*x178 + x174*x179 - 4.0*x183 - 4.0*x184 - 2.0*x185 - 2.0*x186 + x226 + x227 + x229 + x230 + x232 + x233 - x235*x46 + x235*x48 - x236*x46 + x236*x48 - x237*x46 + x237*x48), x54), (mu - x219*x97 + x219*x98 + x238*x26 - x238*x29, x137), (mu - x220*x97 + x220*x98, x138)) */
if (x31)
{
    DEBUG_PRINT("Case (x31) is True.\n");
    /* Assignment result[2, 3]=mu */
    result[11] = mu;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[2, 3]=x151 */
    double x101;
    double x102;
    double x128;
    double x139;
    double x140;
    double x141;
    double x142;
    double x143;
    double x144;
    double x145;
    double x146;
    double x147;
    double x148;
    double x149;
    double x150;
    double x151;
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x40=x32*x33;
    x101=mu*x32;
    x102=x101*x33;
    x128=pow(mu, 3);
    x139=11.313708498984761164*x102 + 16.0*x40;
    x140=1.0/(assert(IS_NOT_ZERO(x139)), x139);
    x141=4.2426406871192847703*x32;
    x142=4.2426406871192847703*x33;
    x143=mu*x33;
    x144=2.0*x143;
    x145=x32*x8;
    x146=x128*x32;
    x147=x142*x8;
    x148=-x147;
    x149=2.0*x128;
    x150=x149*x33;
    x151=x140*(-26.0*x101 + 16.0*x102 - x141 + x142 + x144 + 11.313708498984761164*x145*x33 - 57.982756057296896302*x145 - 34.0*x146 + x148 + x150);
    result[11] = x151;
}
else if (x54)
{
    DEBUG_PRINT("Case (x54) is True.\n");
    /* Assignment result[2, 3]=x155*(-x149*x175 + x149*x176 + x161 - x174*x178 + x174*x179 - 4.0*x183 - 4.0*x184 - 2.0*x185 - 2.0*x186 + x226 + x227 + x229 + x230 + x232 + x233 - x235*x46 + x235*x48 - x236*x46 + x236*x48 - x237*x46 + x237*x48) */
    double x72;
    double x126;
    double x128;
    double x129;
    double x149;
    double x152;
    double x153;
    double x154;
    double x155;
    double x156;
    double x157;
    double x158;
    double x159;
    double x160;
    double x161;
    double x163;
    double x168;
    double x174;
    double x175;
    double x176;
    double x177;
    double x178;
    double x179;
    double x183;
    double x184;
    double x185;
    double x186;
    double x190;
    double x193;
    double x225;
    double x226;
    double x227;
    double x228;
    double x229;
    double x230;
    double x231;
    double x232;
    double x233;
    double x235;
    double x236;
    double x237;
    x44=2.0*(assert(IS_POSITIVE(x43)), sqrt(x43))*x5;
    x45=un*un + x10 + x6 + x7;
    x46=(assert(IS_POSITIVE(-x44 + x45)), sqrt(-x44 + x45));
    x48=(assert(IS_POSITIVE(x44 + x45)), sqrt(x44 + x45));
    x72=x10*x6 + x10*x7;
    x126=pow(x72, 3.0/2.0);
    x128=pow(mu, 3);
    x129=pow(rn, 3);
    x149=2.0*x128;
    x152=4.0*x126*x46*x48;
    x153=x152*x6;
    x154=x152*x7;
    x155=1.0/(assert(IS_NOT_ZERO((x153 + x154))), (x153 + x154));
    x156=pow(mu, 5);
    x157=pow(rn, 4);
    x158=2.0*un*x156*x157*x6*x7;
    x159=pow(un, 3);
    x160=2.0*x128*x159*x6*x7*x9;
    x161=mu*x153 + mu*x154 - x158*x46 + x158*x48 - x160*x46 + x160*x48;
    x163=pow(mu, 4);
    x168=pow(rt2, 4);
    x174=4.0*x128;
    x175=un*x168*x46*x6*x9;
    x176=un*x168*x48*x6*x9;
    x177=pow(rt1, 4);
    x178=un*x177*x46*x7*x9;
    x179=un*x177*x48*x7*x9;
    x183=mu*un*x126*x46*x6;
    x184=mu*un*x126*x48*x6;
    x185=mu*un*x126*x46*x7;
    x186=mu*un*x126*x48*x7;
    x190=pow(rt2, 5);
    x193=pow(rt2, 3);
    x225=1.4142135623730951455*un*x129*x163*x190;
    x226=-x225*x46;
    x227=x225*x48;
    x228=1.4142135623730951455*rt2*un*x129*x163*x177;
    x229=-x228*x46;
    x230=x228*x48;
    x231=2.8284271247461902909*un*x129*x163*x193*x6;
    x232=-x231*x46;
    x233=x231*x48;
    x235=2.0*pow(rt1, 6)*un*x128*x9;
    x236=2.0*un*x156*x157*x177;
    x237=2.0*x128*x159*x177*x9;
    result[11] = x155*(-x149*x175 + x149*x176 + x161 - x174*x178 + x174*x179 - 4.0*x183 - 4.0*x184 - 2.0*x185 - 2.0*x186 + x226 + x227 + x229 + x230 + x232 + x233 - x235*x46 + x235*x48 - x236*x46 + x236*x48 - x237*x46 + x237*x48);
}
else if (x137)
{
    DEBUG_PRINT("Case (x137) is True.\n");
    /* Assignment result[2, 3]=mu - x219*x97 + x219*x98 + x238*x26 - x238*x29 */
    double x78;
    double x81;
    double x83;
    double x84;
    double x91;
    double x92;
    double x93;
    double x94;
    double x95;
    double x96;
    double x97;
    double x98;
    double x202;
    double x219;
    double x223;
    double x238;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x78=1.0/(assert(IS_NOT_ZERO(x3)), x3);
    x81=ut1*ut2*x78*x8;
    x83=mu*x4;
    x84=2*x83;
    x91=x18*x78;
    x92=ut2*x8 + x4*x91;
    x93=x16*x81;
    x94=x2*x78*x8;
    x95=(1.0/2.0)*x19*(x84 + 2*x94);
    x96=x58*(x93 + x95);
    x97=x55*(x92 + x96);
    x98=x61*(x92 - x96);
    x134=pow(x21, -3.0/2.0);
    x202=-x93 - x95;
    x219=x19*x58;
    x223=x134*x19;
    x238=x202*x223 + x58*(x83 + x94);
    result[11] = mu - x219*x97 + x219*x98 + x238*x26 - x238*x29;
}
else if (x138)
{
    DEBUG_PRINT("Case (x138) is True.\n");
    /* Assignment result[2, 3]=mu - x220*x97 + x220*x98 */
    double x78;
    double x81;
    double x83;
    double x84;
    double x91;
    double x92;
    double x93;
    double x94;
    double x95;
    double x96;
    double x97;
    double x98;
    double x220;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x78=1.0/(assert(IS_NOT_ZERO(x3)), x3);
    x81=ut1*ut2*x78*x8;
    x83=mu*x4;
    x84=2*x83;
    x91=x18*x78;
    x92=ut2*x8 + x4*x91;
    x93=x16*x81;
    x94=x2*x78*x8;
    x95=(1.0/2.0)*x19*(x84 + 2*x94);
    x96=x58*(x93 + x95);
    x97=x55*(x92 + x96);
    x98=x61*(x92 - x96);
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x220=random2*x122;
    result[11] = mu - x220*x97 + x220*x98;
}
/* Assignment result[0, 4]=Piecewise((mu, x100), (x35*(x102 - x32*x69 - x36*x8 - x37*x8 + x70), x42), (mu - x105 - x106, x108)) */
if (x100)
{
    DEBUG_PRINT("Case (x100) is True.\n");
    /* Assignment result[0, 4]=mu */
    result[12] = mu;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[0, 4]=x35*(x102 - x32*x69 - x36*x8 - x37*x8 + x70) */
    double x69;
    double x70;
    double x101;
    double x102;
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
    x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
    x36=0.5*x32;
    x37=0.5*x33;
    x69=0.70710678118654757274*mu;
    x70=x33*x69;
    x101=mu*x32;
    x102=x101*x33;
    result[12] = x35*(x102 - x32*x69 - x36*x8 - x37*x8 + x70);
}
else if (x108)
{
    DEBUG_PRINT("Case (x108) is True.\n");
    /* Assignment result[0, 4]=mu - x105 - x106 */
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x103=rn*x8;
    x104=x58*(mu*rt1*x16 + mu*rt2*x19);
    x105=x55*(x103 + x104);
    x106=x61*(x103 - x104);
    result[12] = mu - x105 - x106;
}
/* Assignment result[1, 4]=Piecewise((0.0, x100), (x204, x42), (-x105*x120 + x106*x120 + x206*x26 - x206*x29, x207), (-x105*x123 + x106*x123, x208)) */
if (x100)
{
    DEBUG_PRINT("Case (x100) is True.\n");
    /* Assignment result[1, 4]=0.0 */
    result[13] = 0.0;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[1, 4]=x204 */
    double x109;
    double x110;
    double x204;
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
    x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
    x36=0.5*x32;
    x37=0.5*x33;
    x109=0.35355339059327378637*x32;
    x110=0.35355339059327378637*x33;
    x204=x35*(-mu*x36 - mu*x37 - x109*x8 + x110*x8);
    result[13] = x204;
}
else if (x207)
{
    DEBUG_PRINT("Case (x207) is True.\n");
    /* Assignment result[1, 4]=-x105*x120 + x106*x120 + x206*x26 - x206*x29 */
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x103=rn*x8;
    x104=x58*(mu*rt1*x16 + mu*rt2*x19);
    x105=x55*(x103 + x104);
    x106=x61*(x103 - x104);
    x120=x16*x58;
    x134=pow(x21, -3.0/2.0);
    x135=x134*x16;
    x205=-mu*rt1*x16 - mu*rt2*x19;
    x206=mu*rt1*x58 + x135*x205;
    result[13] = -x105*x120 + x106*x120 + x206*x26 - x206*x29;
}
else if (x208)
{
    DEBUG_PRINT("Case (x208) is True.\n");
    /* Assignment result[1, 4]=-x105*x123 + x106*x123 */
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x103=rn*x8;
    x104=x58*(mu*rt1*x16 + mu*rt2*x19);
    x105=x55*(x103 + x104);
    x106=x61*(x103 - x104);
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x123=random1*x122;
    result[13] = -x105*x123 + x106*x123;
}
/* Assignment result[2, 4]=Piecewise((0.0, x100), (x204, x42), (-x105*x219 + x106*x219 + x239*x26 - x239*x29, x207), (-x105*x220 + x106*x220, x208)) */
if (x100)
{
    DEBUG_PRINT("Case (x100) is True.\n");
    /* Assignment result[2, 4]=0.0 */
    result[14] = 0.0;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[2, 4]=x204 */
    double x109;
    double x110;
    double x204;
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
    x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
    x36=0.5*x32;
    x37=0.5*x33;
    x109=0.35355339059327378637*x32;
    x110=0.35355339059327378637*x33;
    x204=x35*(-mu*x36 - mu*x37 - x109*x8 + x110*x8);
    result[14] = x204;
}
else if (x207)
{
    DEBUG_PRINT("Case (x207) is True.\n");
    /* Assignment result[2, 4]=-x105*x219 + x106*x219 + x239*x26 - x239*x29 */
    double x219;
    double x223;
    double x239;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x103=rn*x8;
    x104=x58*(mu*rt1*x16 + mu*rt2*x19);
    x105=x55*(x103 + x104);
    x106=x61*(x103 - x104);
    x134=pow(x21, -3.0/2.0);
    x205=-mu*rt1*x16 - mu*rt2*x19;
    x219=x19*x58;
    x223=x134*x19;
    x239=mu*rt2*x58 + x205*x223;
    result[14] = -x105*x219 + x106*x219 + x239*x26 - x239*x29;
}
else if (x208)
{
    DEBUG_PRINT("Case (x208) is True.\n");
    /* Assignment result[2, 4]=-x105*x220 + x106*x220 */
    double x220;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x103=rn*x8;
    x104=x58*(mu*rt1*x16 + mu*rt2*x19);
    x105=x55*(x103 + x104);
    x106=x61*(x103 - x104);
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x220=random2*x122;
    result[14] = -x105*x220 + x106*x220;
}
/* Assignment result[0, 5]=Piecewise((0.0, x100), (x111, x42), (-x114 - x115, x108)) */
if (x100)
{
    DEBUG_PRINT("Case (x100) is True.\n");
    /* Assignment result[0, 5]=0.0 */
    result[15] = 0.0;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[0, 5]=x111 */
    double x109;
    double x110;
    double x111;
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
    x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
    x36=0.5*x32;
    x37=0.5*x33;
    x38=-x36 - x37;
    x109=0.35355339059327378637*x32;
    x110=0.35355339059327378637*x33;
    x111=x35*(-mu*x109 + mu*x110 + x38);
    result[15] = x111;
}
else if (x108)
{
    DEBUG_PRINT("Case (x108) is True.\n");
    /* Assignment result[0, 5]=-x114 - x115 */
    double x112;
    double x113;
    double x114;
    double x115;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x112=x5*x58;
    x113=x112*x16;
    x114=x55*(rt1 + x113);
    x115=x61*(rt1 - x113);
    result[15] = -x114 - x115;
}
/* Assignment result[1, 5]=Piecewise((1.00000000000000, x100), (x212, x42), (-x114*x120 + x115*x120 + x214*x26 - x214*x29 + 1, x207), (-x114*x123 + x115*x123 + 1, x208)) */
if (x100)
{
    DEBUG_PRINT("Case (x100) is True.\n");
    /* Assignment result[1, 5]=1.00000000000000 */
    result[16] = 1.0;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[1, 5]=x212 */
    double x101;
    double x102;
    double x139;
    double x140;
    double x143;
    double x145;
    double x209;
    double x210;
    double x211;
    double x212;
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x40=x32*x33;
    x101=mu*x32;
    x102=x101*x33;
    x139=11.313708498984761164*x102 + 16.0*x40;
    x140=1.0/(assert(IS_NOT_ZERO(x139)), x139);
    x143=mu*x33;
    x145=x32*x8;
    x209=9.8994949366116653522*x32;
    x210=4.0*x143;
    x211=1.4142135623730951455*x33;
    x212=x140*(-20.0*x101 + x139 - 15.556349186104045046*x145 - x209 - x210 - x211*x8 + 9.8994949366116653522*x33);
    result[16] = x212;
}
else if (x207)
{
    DEBUG_PRINT("Case (x207) is True.\n");
    /* Assignment result[1, 5]=-x114*x120 + x115*x120 + x214*x26 - x214*x29 + 1 */
    double x112;
    double x113;
    double x114;
    double x115;
    double x213;
    double x214;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x112=x5*x58;
    x113=x112*x16;
    x114=x55*(rt1 + x113);
    x115=x61*(rt1 - x113);
    x120=x16*x58;
    x134=pow(x21, -3.0/2.0);
    x213=mu*rn*x134;
    x214=x112 - x17*x213;
    result[16] = -x114*x120 + x115*x120 + x214*x26 - x214*x29 + 1;
}
else if (x208)
{
    DEBUG_PRINT("Case (x208) is True.\n");
    /* Assignment result[1, 5]=-x114*x123 + x115*x123 + 1 */
    double x112;
    double x113;
    double x114;
    double x115;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x112=x5*x58;
    x113=x112*x16;
    x114=x55*(rt1 + x113);
    x115=x61*(rt1 - x113);
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x123=random1*x122;
    result[16] = -x114*x123 + x115*x123 + 1;
}
/* Assignment result[2, 5]=Piecewise((0.0, x100), (x215, x42), (-x114*x219 + x115*x219 + x217, x207), (-x114*x220 + x115*x220, x208)) */
if (x100)
{
    DEBUG_PRINT("Case (x100) is True.\n");
    /* Assignment result[2, 5]=0.0 */
    result[17] = 0.0;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[2, 5]=x215 */
    double x101;
    double x102;
    double x139;
    double x140;
    double x142;
    double x143;
    double x147;
    double x148;
    double x209;
    double x210;
    double x211;
    double x215;
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x39=1.4142135623730951455*x32;
    x40=x32*x33;
    x101=mu*x32;
    x102=x101*x33;
    x139=11.313708498984761164*x102 + 16.0*x40;
    x140=1.0/(assert(IS_NOT_ZERO(x139)), x139);
    x142=4.2426406871192847703*x33;
    x143=mu*x33;
    x147=x142*x8;
    x148=-x147;
    x209=9.8994949366116653522*x32;
    x210=4.0*x143;
    x211=1.4142135623730951455*x33;
    x215=x140*(4.0*x101 + x148 + x209*x8 + x210 + x211 - x39);
    result[17] = x215;
}
else if (x207)
{
    DEBUG_PRINT("Case (x207) is True.\n");
    /* Assignment result[2, 5]=-x114*x219 + x115*x219 + x217 */
    double x112;
    double x113;
    double x114;
    double x115;
    double x216;
    double x217;
    double x219;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x112=x5*x58;
    x113=x112*x16;
    x114=x55*(rt1 + x113);
    x115=x61*(rt1 - x113);
    x134=pow(x21, -3.0/2.0);
    x216=mu*rn*x134*x16*x19;
    x217=-x216*x26 + x216*x29;
    x219=x19*x58;
    result[17] = -x114*x219 + x115*x219 + x217;
}
else if (x208)
{
    DEBUG_PRINT("Case (x208) is True.\n");
    /* Assignment result[2, 5]=-x114*x220 + x115*x220 */
    double x112;
    double x113;
    double x114;
    double x115;
    double x220;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x112=x5*x58;
    x113=x112*x16;
    x114=x55*(rt1 + x113);
    x115=x61*(rt1 - x113);
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x220=random2*x122;
    result[17] = -x114*x220 + x115*x220;
}
/* Assignment result[0, 6]=Piecewise((0.0, x100), (x111, x42), (-x117 - x118, x108)) */
if (x100)
{
    DEBUG_PRINT("Case (x100) is True.\n");
    /* Assignment result[0, 6]=0.0 */
    result[18] = 0.0;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[0, 6]=x111 */
    double x109;
    double x110;
    double x111;
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
    x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
    x36=0.5*x32;
    x37=0.5*x33;
    x38=-x36 - x37;
    x109=0.35355339059327378637*x32;
    x110=0.35355339059327378637*x33;
    x111=x35*(-mu*x109 + mu*x110 + x38);
    result[18] = x111;
}
else if (x108)
{
    DEBUG_PRINT("Case (x108) is True.\n");
    /* Assignment result[0, 6]=-x117 - x118 */
    double x112;
    double x116;
    double x117;
    double x118;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x112=x5*x58;
    x116=x112*x19;
    x117=x55*(rt2 + x116);
    x118=x61*(rt2 - x116);
    result[18] = -x117 - x118;
}
/* Assignment result[1, 6]=Piecewise((0.0, x100), (x215, x42), (-x117*x120 + x118*x120 + x217, x207), (-x117*x123 + x118*x123, x208)) */
if (x100)
{
    DEBUG_PRINT("Case (x100) is True.\n");
    /* Assignment result[1, 6]=0.0 */
    result[19] = 0.0;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[1, 6]=x215 */
    double x101;
    double x102;
    double x139;
    double x140;
    double x142;
    double x143;
    double x147;
    double x148;
    double x209;
    double x210;
    double x211;
    double x215;
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x39=1.4142135623730951455*x32;
    x40=x32*x33;
    x101=mu*x32;
    x102=x101*x33;
    x139=11.313708498984761164*x102 + 16.0*x40;
    x140=1.0/(assert(IS_NOT_ZERO(x139)), x139);
    x142=4.2426406871192847703*x33;
    x143=mu*x33;
    x147=x142*x8;
    x148=-x147;
    x209=9.8994949366116653522*x32;
    x210=4.0*x143;
    x211=1.4142135623730951455*x33;
    x215=x140*(4.0*x101 + x148 + x209*x8 + x210 + x211 - x39);
    result[19] = x215;
}
else if (x207)
{
    DEBUG_PRINT("Case (x207) is True.\n");
    /* Assignment result[1, 6]=-x117*x120 + x118*x120 + x217 */
    double x112;
    double x116;
    double x117;
    double x118;
    double x216;
    double x217;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x112=x5*x58;
    x116=x112*x19;
    x117=x55*(rt2 + x116);
    x118=x61*(rt2 - x116);
    x120=x16*x58;
    x134=pow(x21, -3.0/2.0);
    x216=mu*rn*x134*x16*x19;
    x217=-x216*x26 + x216*x29;
    result[19] = -x117*x120 + x118*x120 + x217;
}
else if (x208)
{
    DEBUG_PRINT("Case (x208) is True.\n");
    /* Assignment result[1, 6]=-x117*x123 + x118*x123 */
    double x112;
    double x116;
    double x117;
    double x118;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x112=x5*x58;
    x116=x112*x19;
    x117=x55*(rt2 + x116);
    x118=x61*(rt2 - x116);
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x123=random1*x122;
    result[19] = -x117*x123 + x118*x123;
}
/* Assignment result[2, 6]=Piecewise((1.00000000000000, x100), (x212, x42), (-x117*x219 + x118*x219 + x240*x26 - x240*x29 + 1, x207), (-x117*x220 + x118*x220 + 1, x208)) */
if (x100)
{
    DEBUG_PRINT("Case (x100) is True.\n");
    /* Assignment result[2, 6]=1.00000000000000 */
    result[20] = 1.0;
}
else if (x42)
{
    DEBUG_PRINT("Case (x42) is True.\n");
    /* Assignment result[2, 6]=x212 */
    double x101;
    double x102;
    double x139;
    double x140;
    double x143;
    double x145;
    double x209;
    double x210;
    double x211;
    double x212;
    x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
    x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
    x40=x32*x33;
    x101=mu*x32;
    x102=x101*x33;
    x139=11.313708498984761164*x102 + 16.0*x40;
    x140=1.0/(assert(IS_NOT_ZERO(x139)), x139);
    x143=mu*x33;
    x145=x32*x8;
    x209=9.8994949366116653522*x32;
    x210=4.0*x143;
    x211=1.4142135623730951455*x33;
    x212=x140*(-20.0*x101 + x139 - 15.556349186104045046*x145 - x209 - x210 - x211*x8 + 9.8994949366116653522*x33);
    result[20] = x212;
}
else if (x207)
{
    DEBUG_PRINT("Case (x207) is True.\n");
    /* Assignment result[2, 6]=-x117*x219 + x118*x219 + x240*x26 - x240*x29 + 1 */
    double x112;
    double x116;
    double x117;
    double x118;
    double x213;
    double x219;
    double x240;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x112=x5*x58;
    x116=x112*x19;
    x117=x55*(rt2 + x116);
    x118=x61*(rt2 - x116);
    x134=pow(x21, -3.0/2.0);
    x213=mu*rn*x134;
    x219=x19*x58;
    x240=x112 - x20*x213;
    result[20] = -x117*x219 + x118*x219 + x240*x26 - x240*x29 + 1;
}
else if (x208)
{
    DEBUG_PRINT("Case (x208) is True.\n");
    /* Assignment result[2, 6]=-x117*x220 + x118*x220 + 1 */
    double x112;
    double x116;
    double x117;
    double x118;
    double x220;
    x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
    x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
    x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
    x112=x5*x58;
    x116=x112*x19;
    x117=x55*(rt2 + x116);
    x118=x61*(rt2 - x116);
    x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x220=random2*x122;
    result[20] = -x117*x220 + x118*x220 + 1;
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
double x1;
double x2;
double x3;
double x6;
double x7;
double x8;
int x15;

double x4;
double x5;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;

x1=ut1*ut1;
x2=ut2*ut2;
x3=mu*(assert(IS_POSITIVE(x1 + x2)), sqrt(x1 + x2)) + un;
x6=mu*rn*rt1 + mu*ut1*x3;
x7=mu*rn*rt2 + mu*ut2*x3;
x8=(assert(IS_POSITIVE(x6*x6 + x7*x7)), sqrt(x6*x6 + x7*x7));
x15=x8 > 0;

int x18;

double x16;
double x17;

x18=x8 <= 0;

if (x15)
{
    x4=mu*mu;
    x5=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x3*x3;
    x9=2*x8;
    x10=0.5*(assert(IS_POSITIVE(x5 - x9)), sqrt(x5 - x9));
    x11=0.5*(assert(IS_POSITIVE(x5 + x9)), sqrt(x5 + x9));
    x12=mu*ut1 + rt1;
    x13=1.0/(assert(IS_NOT_ZERO(x8)), x8);
    x14=x13*x6;

}
else if (x18)
{
    x4=mu*mu;
    x5=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x3*x3;
    x9=2*x8;
    x10=0.5*(assert(IS_POSITIVE(x5 - x9)), sqrt(x5 - x9));
    x11=0.5*(assert(IS_POSITIVE(x5 + x9)), sqrt(x5 + x9));
    x12=mu*ut1 + rt1;
    x16=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x17=random1*x16;

}/* Assignment result[0, 0]=mu*rn - x10 - x11 + x3 */
x4=mu*mu;
x5=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x3*x3;
x9=2*x8;
x10=0.5*(assert(IS_POSITIVE(x5 - x9)), sqrt(x5 - x9));
x11=0.5*(assert(IS_POSITIVE(x5 + x9)), sqrt(x5 + x9));
result[0] = mu*rn - x10 - x11 + x3;
/* Assignment result[1, 0]=Piecewise((x10*x14 - x11*x14 + x12, x15), (x10*x17 - x11*x17 + x12, x18)) */
if (x15)
{
    DEBUG_PRINT("Case (x15) is True.\n");
    /* Assignment result[1, 0]=x10*x14 - x11*x14 + x12 */
    x12=mu*ut1 + rt1;
    x13=1.0/(assert(IS_NOT_ZERO(x8)), x8);
    x14=x13*x6;
    result[1] = x10*x14 - x11*x14 + x12;
}
else if (x18)
{
    DEBUG_PRINT("Case (x18) is True.\n");
    /* Assignment result[1, 0]=x10*x17 - x11*x17 + x12 */
    x12=mu*ut1 + rt1;
    x16=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x17=random1*x16;
    result[1] = x10*x17 - x11*x17 + x12;
}
/* Assignment result[2, 0]=Piecewise((x10*x20 - x11*x20 + x19, x15), (x10*x21 - x11*x21 + x19, x18)) */
if (x15)
{
    DEBUG_PRINT("Case (x15) is True.\n");
    /* Assignment result[2, 0]=x10*x20 - x11*x20 + x19 */
    double x19;
    double x20;
    x13=1.0/(assert(IS_NOT_ZERO(x8)), x8);
    x19=mu*ut2 + rt2;
    x20=x13*x7;
    result[2] = x10*x20 - x11*x20 + x19;
}
else if (x18)
{
    DEBUG_PRINT("Case (x18) is True.\n");
    /* Assignment result[2, 0]=x10*x21 - x11*x21 + x19 */
    double x19;
    double x21;
    x16=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x19=mu*ut2 + rt2;
    x21=random2*x16;
    result[2] = x10*x21 - x11*x21 + x19;
}
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
double x1;
double x2;
double x3;
double x6;
double x7;
double x8;
int x15;

double x4;
double x5;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;

x1=ut1*ut1;
x2=ut2*ut2;
x3=mu*(assert(IS_POSITIVE(x1 + x2)), sqrt(x1 + x2)) + un;
x6=mu*rn*rt1 + mu*ut1*x3;
x7=mu*rn*rt2 + mu*ut2*x3;
x8=(assert(IS_POSITIVE(x6*x6 + x7*x7)), sqrt(x6*x6 + x7*x7));
x15=x8 > 0;

int x18;

double x16;
double x17;

x18=x8 <= 0;

if (x15)
{
    x4=mu*mu;
    x5=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x3*x3;
    x9=2*x8;
    x10=0.5*(assert(IS_POSITIVE(x5 - x9)), sqrt(x5 - x9));
    x11=0.5*(assert(IS_POSITIVE(x5 + x9)), sqrt(x5 + x9));
    x12=mu*ut1 + rt1;
    x13=1.0/(assert(IS_NOT_ZERO(x8)), x8);
    x14=x13*x6;

}
else if (x18)
{
    x4=mu*mu;
    x5=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x3*x3;
    x9=2*x8;
    x10=0.5*(assert(IS_POSITIVE(x5 - x9)), sqrt(x5 - x9));
    x11=0.5*(assert(IS_POSITIVE(x5 + x9)), sqrt(x5 + x9));
    x12=mu*ut1 + rt1;
    x16=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x17=random1*x16;

}/* Assignment result[0, 0]=mu*rn - x10 - x11 + x3 */
x4=mu*mu;
x5=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x3*x3;
x9=2*x8;
x10=0.5*(assert(IS_POSITIVE(x5 - x9)), sqrt(x5 - x9));
x11=0.5*(assert(IS_POSITIVE(x5 + x9)), sqrt(x5 + x9));
result[0] = mu*rn - x10 - x11 + x3;
/* Assignment result[1, 0]=Piecewise((x10*x14 - x11*x14 + x12, x15), (x10*x17 - x11*x17 + x12, x18)) */
if (x15)
{
    DEBUG_PRINT("Case (x15) is True.\n");
    /* Assignment result[1, 0]=x10*x14 - x11*x14 + x12 */
    x12=mu*ut1 + rt1;
    x13=1.0/(assert(IS_NOT_ZERO(x8)), x8);
    x14=x13*x6;
    result[1] = x10*x14 - x11*x14 + x12;
}
else if (x18)
{
    DEBUG_PRINT("Case (x18) is True.\n");
    /* Assignment result[1, 0]=x10*x17 - x11*x17 + x12 */
    x12=mu*ut1 + rt1;
    x16=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x17=random1*x16;
    result[1] = x10*x17 - x11*x17 + x12;
}
/* Assignment result[2, 0]=Piecewise((x10*x20 - x11*x20 + x19, x15), (x10*x21 - x11*x21 + x19, x18)) */
if (x15)
{
    DEBUG_PRINT("Case (x15) is True.\n");
    /* Assignment result[2, 0]=x10*x20 - x11*x20 + x19 */
    double x19;
    double x20;
    x13=1.0/(assert(IS_NOT_ZERO(x8)), x8);
    x19=mu*ut2 + rt2;
    x20=x13*x7;
    result[2] = x10*x20 - x11*x20 + x19;
}
else if (x18)
{
    DEBUG_PRINT("Case (x18) is True.\n");
    /* Assignment result[2, 0]=x10*x21 - x11*x21 + x19 */
    double x19;
    double x21;
    x16=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x19=mu*ut2 + rt2;
    x21=random2*x16;
    result[2] = x10*x21 - x11*x21 + x19;
}
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
