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

// sqrt(DBL_EPSILON)
#define ZERO 1.490116119384766e-08
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
#define sqrt(x) ((x < 0) && ( x > - ZERO) ? 0 : (assert(x>=0),sqrt(x)))

void fc3d_FischerBurmeisterFABGenerated(
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
    x31=x30 <= 0;
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
    x42=x22 <= 0 || x27 <= 0 || x41 <= 0;
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
    x54=x3 <= 0;
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
    x63=x3 > 0;
    x64=x30 > 0;
    x65=x22 > 0;
    x66=x27 > 0;
    x67=x41 > 0;
    x68=x63 && x64 && x65 && x66 && x67;
    double x99;
    int x100;
    x43=x6 + x7;
    x99=(assert(IS_POSITIVE(x10 + x43)), sqrt(x10 + x43));
    x100=x99 <= 0;
    int x107;
    int x108;
    double x103;
    double x104;
    double x105;
    double x106;
    x107=x99 > 0;
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
    /* Assignment result[0, 0]=-x26 - x29 + x4 + x5 */
    x5=mu*rn;
    x25=(assert(IS_POSITIVE(x24)), sqrt(x24));
    x26=0.5*x25;
    x28=(assert(IS_POSITIVE(x27)), sqrt(x27));
    x29=0.5*x28;result[0] = -x26 - x29 + x4 + x5;

    /* Assignment result[1, 0]=Piecewise((x119 + x120*x26 - x120*x29, x121), (x119 + x123*x26 - x123*x29, x124)) */
    if (x121)
    {
        DEBUG_PRINT("Case (x121) is True.\n");
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x119=mu*ut1 + rt1;
        x120=x16*x58;

        /* Assignment result[1, 0]=x119 + x120*x26 - x120*x29 */
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x119=mu*ut1 + rt1;
        x120=x16*x58;result[1] = x119 + x120*x26 - x120*x29;
    }
    else if (x124)
    {
        DEBUG_PRINT("Case (x124) is True.\n");
        x119=mu*ut1 + rt1;
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x123=random1*x122;

        /* Assignment result[1, 0]=x119 + x123*x26 - x123*x29 */
        x119=mu*ut1 + rt1;
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x123=random1*x122;result[1] = x119 + x123*x26 - x123*x29;
    }

    /* Assignment result[2, 0]=Piecewise((x218 + x219*x26 - x219*x29, x121), (x218 + x220*x26 - x220*x29, x124)) */
    double x218;
    double x219;
    double x220;if (x121)
    {
        DEBUG_PRINT("Case (x121) is True.\n");
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x218=mu*ut2 + rt2;
        x219=x19*x58;

        /* Assignment result[2, 0]=x218 + x219*x26 - x219*x29 */
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x218=mu*ut2 + rt2;
        x219=x19*x58;result[2] = x218 + x219*x26 - x219*x29;
    }
    else if (x124)
    {
        DEBUG_PRINT("Case (x124) is True.\n");
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x218=mu*ut2 + rt2;
        x220=random2*x122;

        /* Assignment result[2, 0]=x218 + x220*x26 - x220*x29 */
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x218=mu*ut2 + rt2;
        x220=random2*x122;result[2] = x218 + x220*x26 - x220*x29;
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
        x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
        x36=0.5*x32;
        x37=0.5*x33;
        x38=-x36 - x37;
        x39=1.4142135623730951455*x32;
        x40=x32*x33;

        /* Assignment result[0, 1]=x35*(-mu*x39 + x38 + x40) */
        x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
        x36=0.5*x32;
        x37=0.5*x33;
        x38=-x36 - x37;
        x39=1.4142135623730951455*x32;
        x40=x32*x33;result[3] = x35*(-mu*x39 + x38 + x40);
    }
    else if (x54)
    {
        DEBUG_PRINT("Case (x54) is True.\n");
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
        x53=x46*x48;result[3] = x47*x49*(-x51 - x52 + x53);
    }
    else if (x68)
    {
        DEBUG_PRINT("Case (x68) is True.\n");
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x56=x15*x16;
        x57=x18*x19;
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x59=x58*(x56 + x57);
        x60=x55*(x4 + x59);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x62=x61*(x4 - x59);

        /* Assignment result[0, 1]=-x60 - x62 + 1 */
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x56=x15*x16;
        x57=x18*x19;
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x59=x58*(x56 + x57);
        x60=x55*(x4 + x59);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x62=x61*(x4 - x59);result[3] = -x60 - x62 + 1;
    }

    /* Assignment result[1, 1]=Piecewise((0.0, x31), (x125, x42), (x127*(-x131*x51 + x131*x52 - x132*x51 + x132*x52), x54), (-x120*x60 + x120*x62 + x136*x26 - x136*x29, x137), (-x123*x60 + x123*x62, x138)) */
    double x72;
    double x101;
    double x109;
    double x110;
    double x125;
    double x126;
    double x127;
    double x128;
    double x129;
    double x130;
    double x131;
    double x132;if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");

        /* Assignment result[1, 1]=0.0 */
        result[4] = 0.0;
    }
    else if (x42)
    {
        DEBUG_PRINT("Case (x42) is True.\n");
        x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
        x101=mu*x32;
        x109=0.35355339059327378637*x32;
        x110=0.35355339059327378637*x33;
        x125=x35*(-1.0*x101 - x109 + x110);

        /* Assignment result[1, 1]=x125 */
        x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
        x101=mu*x32;
        x109=0.35355339059327378637*x32;
        x110=0.35355339059327378637*x33;
        x125=x35*(-1.0*x101 - x109 + x110);result[4] = x125;
    }
    else if (x54)
    {
        DEBUG_PRINT("Case (x54) is True.\n");
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

        /* Assignment result[1, 1]=x127*(-x131*x51 + x131*x52 - x132*x51 + x132*x52) */
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
        x132=rt1*x128*x129*x7;result[4] = x127*(-x131*x51 + x131*x52 - x132*x51 + x132*x52);
    }
    else if (x137)
    {
        DEBUG_PRINT("Case (x137) is True.\n");
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
        x136=x133*x135 + x15*x58;result[4] = -x120*x60 + x120*x62 + x136*x26 - x136*x29;
    }
    else if (x138)
    {
        DEBUG_PRINT("Case (x138) is True.\n");
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
        x123=random1*x122;result[4] = -x123*x60 + x123*x62;
    }

    /* Assignment result[2, 1]=Piecewise((0.0, x31), (x125, x42), (x127*(-x221*x51 + x221*x52 - x222*x51 + x222*x52), x54), (-x219*x60 + x219*x62 + x224*x26 - x224*x29, x137), (-x220*x60 + x220*x62, x138)) */
    double x193;
    double x221;
    double x222;
    double x223;
    double x224;if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");

        /* Assignment result[2, 1]=0.0 */
        result[5] = 0.0;
    }
    else if (x42)
    {
        DEBUG_PRINT("Case (x42) is True.\n");
        x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
        x101=mu*x32;
        x109=0.35355339059327378637*x32;
        x110=0.35355339059327378637*x33;
        x125=x35*(-1.0*x101 - x109 + x110);

        /* Assignment result[2, 1]=x125 */
        x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
        x101=mu*x32;
        x109=0.35355339059327378637*x32;
        x110=0.35355339059327378637*x33;
        x125=x35*(-1.0*x101 - x109 + x110);result[5] = x125;
    }
    else if (x54)
    {
        DEBUG_PRINT("Case (x54) is True.\n");
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

        /* Assignment result[2, 1]=x127*(-x221*x51 + x221*x52 - x222*x51 + x222*x52) */
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
        x222=rt2*x128*x129*x6;result[5] = x127*(-x221*x51 + x221*x52 - x222*x51 + x222*x52);
    }
    else if (x137)
    {
        DEBUG_PRINT("Case (x137) is True.\n");
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

        /* Assignment result[2, 1]=-x219*x60 + x219*x62 + x224*x26 - x224*x29 */
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
        x224=x133*x223 + x18*x58;result[5] = -x219*x60 + x219*x62 + x224*x26 - x224*x29;
    }
    else if (x138)
    {
        DEBUG_PRINT("Case (x138) is True.\n");
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

        /* Assignment result[2, 1]=-x220*x60 + x220*x62 */
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x56=x15*x16;
        x57=x18*x19;
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x59=x58*(x56 + x57);
        x60=x55*(x4 + x59);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x62=x61*(x4 - x59);
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x220=random2*x122;result[5] = -x220*x60 + x220*x62;
    }

    /* Assignment result[0, 2]=Piecewise((x69, x31), (x71, x42), (x74*(-x51*x75 + x52*x75 + x77), x54), (x79 - x88 - x89, x68)) */
    double x69;
    double x70;
    double x71;
    double x73;
    double x74;
    double x75;
    double x76;
    double x77;
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
    double x89;if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x69=0.70710678118654757274*mu;

        /* Assignment result[0, 2]=x69 */
        x69=0.70710678118654757274*mu;result[6] = x69;
    }
    else if (x42)
    {
        DEBUG_PRINT("Case (x42) is True.\n");
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x69=0.70710678118654757274*mu;
        x70=x33*x69;
        x71=x34*(-x69 + x70 - 2.0*x8);

        /* Assignment result[0, 2]=x71 */
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x69=0.70710678118654757274*mu;
        x70=x33*x69;
        x71=x34*(-x69 + x70 - 2.0*x8);result[6] = x71;
    }
    else if (x54)
    {
        DEBUG_PRINT("Case (x54) is True.\n");
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

        /* Assignment result[0, 2]=x74*(-x51*x75 + x52*x75 + x77) */
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
        x77=0.70710678118654757274*mu*x53*x73 - x46*x76 - x48*x76;result[6] = x74*(-x51*x75 + x52*x75 + x77);
    }
    else if (x68)
    {
        DEBUG_PRINT("Case (x68) is True.\n");
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

        /* Assignment result[0, 2]=x79 - x88 - x89 */
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
        x89=x61*(x80 - x87);result[6] = x79 - x88 - x89;
    }

    /* Assignment result[1, 2]=Piecewise((mu, x31), (x151, x42), (x155*(-x149*x178 + x149*x179 + x161 - x162*x46 + x162*x48 + x166 + x167 - x169*x46 + x169*x48 + x171 + x172 - x173*x46 + x173*x48 - x174*x175 + x174*x176 + x181 + x182 - 2.0*x183 - 2.0*x184 - 4.0*x185 - 4.0*x186), x54), (mu - x120*x88 + x120*x89 + x188*x26 - x188*x29, x137), (mu - x123*x88 + x123*x89, x138)) */
    double x102;
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
    double x187;
    double x188;if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");

        /* Assignment result[1, 2]=mu */
        result[7] = mu;
    }
    else if (x42)
    {
        DEBUG_PRINT("Case (x42) is True.\n");
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

        /* Assignment result[1, 2]=x151 */
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
        x151=x140*(-26.0*x101 + 16.0*x102 - x141 + x142 + x144 + 11.313708498984761164*x145*x33 - 57.982756057296896302*x145 - 34.0*x146 + x148 + x150);result[7] = x151;
    }
    else if (x54)
    {
        DEBUG_PRINT("Case (x54) is True.\n");
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

        /* Assignment result[1, 2]=x155*(-x149*x178 + x149*x179 + x161 - x162*x46 + x162*x48 + x166 + x167 - x169*x46 + x169*x48 + x171 + x172 - x173*x46 + x173*x48 - x174*x175 + x174*x176 + x181 + x182 - 2.0*x183 - 2.0*x184 - 4.0*x185 - 4.0*x186) */
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
        x186=mu*un*x126*x48*x7;result[7] = x155*(-x149*x178 + x149*x179 + x161 - x162*x46 + x162*x48 + x166 + x167 - x169*x46 + x169*x48 + x171 + x172 - x173*x46 + x173*x48 - x174*x175 + x174*x176 + x181 + x182 - 2.0*x183 - 2.0*x184 - 4.0*x185 - 4.0*x186);
    }
    else if (x137)
    {
        DEBUG_PRINT("Case (x137) is True.\n");
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

        /* Assignment result[1, 2]=mu - x120*x88 + x120*x89 + x188*x26 - x188*x29 */
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
        x188=x135*x187 + x58*(x83 + x85);result[7] = mu - x120*x88 + x120*x89 + x188*x26 - x188*x29;
    }
    else if (x138)
    {
        DEBUG_PRINT("Case (x138) is True.\n");
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

        /* Assignment result[1, 2]=mu - x123*x88 + x123*x89 */
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
        x123=random1*x122;result[7] = mu - x123*x88 + x123*x89;
    }

    /* Assignment result[2, 2]=Piecewise((0.0, x31), (x189, x42), (x155*(x200 + x226 + x227 + x229 + x230 + x232 + x233), x54), (-x219*x88 + x219*x89 + x234*x26 - x234*x29, x137), (-x220*x88 + x220*x89, x138)) */
    double x189;
    double x190;
    double x191;
    double x192;
    double x194;
    double x195;
    double x196;
    double x197;
    double x198;
    double x199;
    double x200;
    double x201;
    double x225;
    double x226;
    double x227;
    double x228;
    double x229;
    double x230;
    double x231;
    double x232;
    double x233;
    double x234;if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");

        /* Assignment result[2, 2]=0.0 */
        result[8] = 0.0;
    }
    else if (x42)
    {
        DEBUG_PRINT("Case (x42) is True.\n");
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

        /* Assignment result[2, 2]=x189 */
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
        x189=x140*(10.0*x101 + x141 - x142 - x144 + 2.0*x146 + x147 - x150 + x39*x8);result[8] = x189;
    }
    else if (x54)
    {
        DEBUG_PRINT("Case (x54) is True.\n");
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

        /* Assignment result[2, 2]=x155*(x200 + x226 + x227 + x229 + x230 + x232 + x233) */
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
        x233=x231*x48;result[8] = x155*(x200 + x226 + x227 + x229 + x230 + x232 + x233);
    }
    else if (x137)
    {
        DEBUG_PRINT("Case (x137) is True.\n");
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

        /* Assignment result[2, 2]=-x219*x88 + x219*x89 + x234*x26 - x234*x29 */
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
        x234=x187*x223 + x201;result[8] = -x219*x88 + x219*x89 + x234*x26 - x234*x29;
    }
    else if (x138)
    {
        DEBUG_PRINT("Case (x138) is True.\n");
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

        /* Assignment result[2, 2]=-x220*x88 + x220*x89 */
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
        x220=random2*x122;result[8] = -x220*x88 + x220*x89;
    }

    /* Assignment result[0, 3]=Piecewise((x69, x31), (x71, x42), (x74*(-x51*x90 + x52*x90 + x77), x54), (x91 - x97 - x98, x68)) */
    double x90;
    double x91;
    double x92;
    double x93;
    double x94;
    double x95;
    double x96;
    double x97;
    double x98;if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x69=0.70710678118654757274*mu;

        /* Assignment result[0, 3]=x69 */
        x69=0.70710678118654757274*mu;result[9] = x69;
    }
    else if (x42)
    {
        DEBUG_PRINT("Case (x42) is True.\n");
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x69=0.70710678118654757274*mu;
        x70=x33*x69;
        x71=x34*(-x69 + x70 - 2.0*x8);

        /* Assignment result[0, 3]=x71 */
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x69=0.70710678118654757274*mu;
        x70=x33*x69;
        x71=x34*(-x69 + x70 - 2.0*x8);result[9] = x71;
    }
    else if (x54)
    {
        DEBUG_PRINT("Case (x54) is True.\n");
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

        /* Assignment result[0, 3]=x74*(-x51*x90 + x52*x90 + x77) */
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
        x90=rn*rt2*x8;result[9] = x74*(-x51*x90 + x52*x90 + x77);
    }
    else if (x68)
    {
        DEBUG_PRINT("Case (x68) is True.\n");
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

        /* Assignment result[0, 3]=x91 - x97 - x98 */
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
        x98=x61*(x92 - x96);result[9] = x91 - x97 - x98;
    }

    /* Assignment result[1, 3]=Piecewise((0.0, x31), (x189, x42), (x155*(x166 + x167 + x171 + x172 + x181 + x182 + x200), x54), (-x120*x97 + x120*x98 + x203*x26 - x203*x29, x137), (-x123*x97 + x123*x98, x138)) */
    double x202;
    double x203;if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");

        /* Assignment result[1, 3]=0.0 */
        result[10] = 0.0;
    }
    else if (x42)
    {
        DEBUG_PRINT("Case (x42) is True.\n");
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

        /* Assignment result[1, 3]=x189 */
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
        x189=x140*(10.0*x101 + x141 - x142 - x144 + 2.0*x146 + x147 - x150 + x39*x8);result[10] = x189;
    }
    else if (x54)
    {
        DEBUG_PRINT("Case (x54) is True.\n");
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

        /* Assignment result[1, 3]=x155*(x166 + x167 + x171 + x172 + x181 + x182 + x200) */
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
        x200=x191*x46 - x191*x48 + x192*x46 - x192*x48 + x194*x46 - x194*x48 + x195*x46 - x195*x48 + x196*x46 - x196*x48 + x197*x46 - x197*x48 + x198*x46 - x198*x48 + x199*x46 + x199*x48;result[10] = x155*(x166 + x167 + x171 + x172 + x181 + x182 + x200);
    }
    else if (x137)
    {
        DEBUG_PRINT("Case (x137) is True.\n");
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

        /* Assignment result[1, 3]=-x120*x97 + x120*x98 + x203*x26 - x203*x29 */
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
        x203=x135*x202 + x201;result[10] = -x120*x97 + x120*x98 + x203*x26 - x203*x29;
    }
    else if (x138)
    {
        DEBUG_PRINT("Case (x138) is True.\n");
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

        /* Assignment result[1, 3]=-x123*x97 + x123*x98 */
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
        x123=random1*x122;result[10] = -x123*x97 + x123*x98;
    }

    /* Assignment result[2, 3]=Piecewise((mu, x31), (x151, x42), (x155*(-x149*x175 + x149*x176 + x161 - x174*x178 + x174*x179 - 4.0*x183 - 4.0*x184 - 2.0*x185 - 2.0*x186 + x226 + x227 + x229 + x230 + x232 + x233 - x235*x46 + x235*x48 - x236*x46 + x236*x48 - x237*x46 + x237*x48), x54), (mu - x219*x97 + x219*x98 + x238*x26 - x238*x29, x137), (mu - x220*x97 + x220*x98, x138)) */
    double x235;
    double x236;
    double x237;
    double x238;if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");

        /* Assignment result[2, 3]=mu */
        result[11] = mu;
    }
    else if (x42)
    {
        DEBUG_PRINT("Case (x42) is True.\n");
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

        /* Assignment result[2, 3]=x151 */
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
        x151=x140*(-26.0*x101 + 16.0*x102 - x141 + x142 + x144 + 11.313708498984761164*x145*x33 - 57.982756057296896302*x145 - 34.0*x146 + x148 + x150);result[11] = x151;
    }
    else if (x54)
    {
        DEBUG_PRINT("Case (x54) is True.\n");
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

        /* Assignment result[2, 3]=x155*(-x149*x175 + x149*x176 + x161 - x174*x178 + x174*x179 - 4.0*x183 - 4.0*x184 - 2.0*x185 - 2.0*x186 + x226 + x227 + x229 + x230 + x232 + x233 - x235*x46 + x235*x48 - x236*x46 + x236*x48 - x237*x46 + x237*x48) */
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
        x237=2.0*x128*x159*x177*x9;result[11] = x155*(-x149*x175 + x149*x176 + x161 - x174*x178 + x174*x179 - 4.0*x183 - 4.0*x184 - 2.0*x185 - 2.0*x186 + x226 + x227 + x229 + x230 + x232 + x233 - x235*x46 + x235*x48 - x236*x46 + x236*x48 - x237*x46 + x237*x48);
    }
    else if (x137)
    {
        DEBUG_PRINT("Case (x137) is True.\n");
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

        /* Assignment result[2, 3]=mu - x219*x97 + x219*x98 + x238*x26 - x238*x29 */
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
        x238=x202*x223 + x58*(x83 + x94);result[11] = mu - x219*x97 + x219*x98 + x238*x26 - x238*x29;
    }
    else if (x138)
    {
        DEBUG_PRINT("Case (x138) is True.\n");
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

        /* Assignment result[2, 3]=mu - x220*x97 + x220*x98 */
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
        x220=random2*x122;result[11] = mu - x220*x97 + x220*x98;
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

        /* Assignment result[0, 4]=x35*(x102 - x32*x69 - x36*x8 - x37*x8 + x70) */
        x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
        x36=0.5*x32;
        x37=0.5*x33;
        x69=0.70710678118654757274*mu;
        x70=x33*x69;
        x101=mu*x32;
        x102=x101*x33;result[12] = x35*(x102 - x32*x69 - x36*x8 - x37*x8 + x70);
    }
    else if (x108)
    {
        DEBUG_PRINT("Case (x108) is True.\n");
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x103=rn*x8;
        x104=x58*(mu*rt1*x16 + mu*rt2*x19);
        x105=x55*(x103 + x104);
        x106=x61*(x103 - x104);

        /* Assignment result[0, 4]=mu - x105 - x106 */
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x103=rn*x8;
        x104=x58*(mu*rt1*x16 + mu*rt2*x19);
        x105=x55*(x103 + x104);
        x106=x61*(x103 - x104);result[12] = mu - x105 - x106;
    }

    /* Assignment result[1, 4]=Piecewise((0.0, x100), (x204, x42), (-x105*x120 + x106*x120 + x206*x26 - x206*x29, x207), (-x105*x123 + x106*x123, x208)) */
    double x204;if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");

        /* Assignment result[1, 4]=0.0 */
        result[13] = 0.0;
    }
    else if (x42)
    {
        DEBUG_PRINT("Case (x42) is True.\n");
        x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
        x36=0.5*x32;
        x37=0.5*x33;
        x109=0.35355339059327378637*x32;
        x110=0.35355339059327378637*x33;
        x204=x35*(-mu*x36 - mu*x37 - x109*x8 + x110*x8);

        /* Assignment result[1, 4]=x204 */
        x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
        x36=0.5*x32;
        x37=0.5*x33;
        x109=0.35355339059327378637*x32;
        x110=0.35355339059327378637*x33;
        x204=x35*(-mu*x36 - mu*x37 - x109*x8 + x110*x8);result[13] = x204;
    }
    else if (x207)
    {
        DEBUG_PRINT("Case (x207) is True.\n");
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
        x206=mu*rt1*x58 + x135*x205;result[13] = -x105*x120 + x106*x120 + x206*x26 - x206*x29;
    }
    else if (x208)
    {
        DEBUG_PRINT("Case (x208) is True.\n");
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x103=rn*x8;
        x104=x58*(mu*rt1*x16 + mu*rt2*x19);
        x105=x55*(x103 + x104);
        x106=x61*(x103 - x104);
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x123=random1*x122;

        /* Assignment result[1, 4]=-x105*x123 + x106*x123 */
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x103=rn*x8;
        x104=x58*(mu*rt1*x16 + mu*rt2*x19);
        x105=x55*(x103 + x104);
        x106=x61*(x103 - x104);
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x123=random1*x122;result[13] = -x105*x123 + x106*x123;
    }

    /* Assignment result[2, 4]=Piecewise((0.0, x100), (x204, x42), (-x105*x219 + x106*x219 + x239*x26 - x239*x29, x207), (-x105*x220 + x106*x220, x208)) */
    double x239;if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");

        /* Assignment result[2, 4]=0.0 */
        result[14] = 0.0;
    }
    else if (x42)
    {
        DEBUG_PRINT("Case (x42) is True.\n");
        x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
        x36=0.5*x32;
        x37=0.5*x33;
        x109=0.35355339059327378637*x32;
        x110=0.35355339059327378637*x33;
        x204=x35*(-mu*x36 - mu*x37 - x109*x8 + x110*x8);

        /* Assignment result[2, 4]=x204 */
        x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
        x36=0.5*x32;
        x37=0.5*x33;
        x109=0.35355339059327378637*x32;
        x110=0.35355339059327378637*x33;
        x204=x35*(-mu*x36 - mu*x37 - x109*x8 + x110*x8);result[14] = x204;
    }
    else if (x207)
    {
        DEBUG_PRINT("Case (x207) is True.\n");
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

        /* Assignment result[2, 4]=-x105*x219 + x106*x219 + x239*x26 - x239*x29 */
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
        x239=mu*rt2*x58 + x205*x223;result[14] = -x105*x219 + x106*x219 + x239*x26 - x239*x29;
    }
    else if (x208)
    {
        DEBUG_PRINT("Case (x208) is True.\n");
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x103=rn*x8;
        x104=x58*(mu*rt1*x16 + mu*rt2*x19);
        x105=x55*(x103 + x104);
        x106=x61*(x103 - x104);
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x220=random2*x122;

        /* Assignment result[2, 4]=-x105*x220 + x106*x220 */
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x103=rn*x8;
        x104=x58*(mu*rt1*x16 + mu*rt2*x19);
        x105=x55*(x103 + x104);
        x106=x61*(x103 - x104);
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x220=random2*x122;result[14] = -x105*x220 + x106*x220;
    }

    /* Assignment result[0, 5]=Piecewise((0.0, x100), (x111, x42), (-x114 - x115, x108)) */
    double x111;
    double x112;
    double x113;
    double x114;
    double x115;if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");

        /* Assignment result[0, 5]=0.0 */
        result[15] = 0.0;
    }
    else if (x42)
    {
        DEBUG_PRINT("Case (x42) is True.\n");
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

        /* Assignment result[0, 5]=x111 */
        x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
        x36=0.5*x32;
        x37=0.5*x33;
        x38=-x36 - x37;
        x109=0.35355339059327378637*x32;
        x110=0.35355339059327378637*x33;
        x111=x35*(-mu*x109 + mu*x110 + x38);result[15] = x111;
    }
    else if (x108)
    {
        DEBUG_PRINT("Case (x108) is True.\n");
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x112=x5*x58;
        x113=x112*x16;
        x114=x55*(rt1 + x113);
        x115=x61*(rt1 - x113);

        /* Assignment result[0, 5]=-x114 - x115 */
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x112=x5*x58;
        x113=x112*x16;
        x114=x55*(rt1 + x113);
        x115=x61*(rt1 - x113);result[15] = -x114 - x115;
    }

    /* Assignment result[1, 5]=Piecewise((1.00000000000000, x100), (x212, x42), (-x114*x120 + x115*x120 + x214*x26 - x214*x29 + 1, x207), (-x114*x123 + x115*x123 + 1, x208)) */
    double x209;
    double x210;
    double x211;
    double x212;
    double x213;
    double x214;if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");

        /* Assignment result[1, 5]=1.00000000000000 */
        result[16] = 1.0;
    }
    else if (x42)
    {
        DEBUG_PRINT("Case (x42) is True.\n");
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

        /* Assignment result[1, 5]=x212 */
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
        x212=x140*(-20.0*x101 + x139 - 15.556349186104045046*x145 - x209 - x210 - x211*x8 + 9.8994949366116653522*x33);result[16] = x212;
    }
    else if (x207)
    {
        DEBUG_PRINT("Case (x207) is True.\n");
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

        /* Assignment result[1, 5]=-x114*x120 + x115*x120 + x214*x26 - x214*x29 + 1 */
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
        x214=x112 - x17*x213;result[16] = -x114*x120 + x115*x120 + x214*x26 - x214*x29 + 1;
    }
    else if (x208)
    {
        DEBUG_PRINT("Case (x208) is True.\n");
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x112=x5*x58;
        x113=x112*x16;
        x114=x55*(rt1 + x113);
        x115=x61*(rt1 - x113);
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x123=random1*x122;

        /* Assignment result[1, 5]=-x114*x123 + x115*x123 + 1 */
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x112=x5*x58;
        x113=x112*x16;
        x114=x55*(rt1 + x113);
        x115=x61*(rt1 - x113);
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x123=random1*x122;result[16] = -x114*x123 + x115*x123 + 1;
    }

    /* Assignment result[2, 5]=Piecewise((0.0, x100), (x215, x42), (-x114*x219 + x115*x219 + x217, x207), (-x114*x220 + x115*x220, x208)) */
    double x215;
    double x216;
    double x217;if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");

        /* Assignment result[2, 5]=0.0 */
        result[17] = 0.0;
    }
    else if (x42)
    {
        DEBUG_PRINT("Case (x42) is True.\n");
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

        /* Assignment result[2, 5]=x215 */
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
        x215=x140*(4.0*x101 + x148 + x209*x8 + x210 + x211 - x39);result[17] = x215;
    }
    else if (x207)
    {
        DEBUG_PRINT("Case (x207) is True.\n");
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

        /* Assignment result[2, 5]=-x114*x219 + x115*x219 + x217 */
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
        x219=x19*x58;result[17] = -x114*x219 + x115*x219 + x217;
    }
    else if (x208)
    {
        DEBUG_PRINT("Case (x208) is True.\n");
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x112=x5*x58;
        x113=x112*x16;
        x114=x55*(rt1 + x113);
        x115=x61*(rt1 - x113);
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x220=random2*x122;

        /* Assignment result[2, 5]=-x114*x220 + x115*x220 */
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x112=x5*x58;
        x113=x112*x16;
        x114=x55*(rt1 + x113);
        x115=x61*(rt1 - x113);
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x220=random2*x122;result[17] = -x114*x220 + x115*x220;
    }

    /* Assignment result[0, 6]=Piecewise((0.0, x100), (x111, x42), (-x117 - x118, x108)) */
    double x116;
    double x117;
    double x118;if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");

        /* Assignment result[0, 6]=0.0 */
        result[18] = 0.0;
    }
    else if (x42)
    {
        DEBUG_PRINT("Case (x42) is True.\n");
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

        /* Assignment result[0, 6]=x111 */
        x32=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x8 + 3.0)), sqrt(-2.8284271247461902909*mu + x8 + 3.0));
        x33=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x8 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x8 + 3.0));
        x34=1.0/(assert(IS_NOT_ZERO(x33)), x33);
        x35=(assert(IS_NOT_ZERO(x32)), x34/x32);
        x36=0.5*x32;
        x37=0.5*x33;
        x38=-x36 - x37;
        x109=0.35355339059327378637*x32;
        x110=0.35355339059327378637*x33;
        x111=x35*(-mu*x109 + mu*x110 + x38);result[18] = x111;
    }
    else if (x108)
    {
        DEBUG_PRINT("Case (x108) is True.\n");
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x112=x5*x58;
        x116=x112*x19;
        x117=x55*(rt2 + x116);
        x118=x61*(rt2 - x116);

        /* Assignment result[0, 6]=-x117 - x118 */
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x112=x5*x58;
        x116=x112*x19;
        x117=x55*(rt2 + x116);
        x118=x61*(rt2 - x116);result[18] = -x117 - x118;
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

        /* Assignment result[1, 6]=x215 */
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
        x215=x140*(4.0*x101 + x148 + x209*x8 + x210 + x211 - x39);result[19] = x215;
    }
    else if (x207)
    {
        DEBUG_PRINT("Case (x207) is True.\n");
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

        /* Assignment result[1, 6]=-x117*x120 + x118*x120 + x217 */
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
        x217=-x216*x26 + x216*x29;result[19] = -x117*x120 + x118*x120 + x217;
    }
    else if (x208)
    {
        DEBUG_PRINT("Case (x208) is True.\n");
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x112=x5*x58;
        x116=x112*x19;
        x117=x55*(rt2 + x116);
        x118=x61*(rt2 - x116);
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x123=random1*x122;

        /* Assignment result[1, 6]=-x117*x123 + x118*x123 */
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x112=x5*x58;
        x116=x112*x19;
        x117=x55*(rt2 + x116);
        x118=x61*(rt2 - x116);
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x123=random1*x122;result[19] = -x117*x123 + x118*x123;
    }

    /* Assignment result[2, 6]=Piecewise((1.00000000000000, x100), (x212, x42), (-x117*x219 + x118*x219 + x240*x26 - x240*x29 + 1, x207), (-x117*x220 + x118*x220 + 1, x208)) */
    double x240;if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");

        /* Assignment result[2, 6]=1.00000000000000 */
        result[20] = 1.0;
    }
    else if (x42)
    {
        DEBUG_PRINT("Case (x42) is True.\n");
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

        /* Assignment result[2, 6]=x212 */
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
        x212=x140*(-20.0*x101 + x139 - 15.556349186104045046*x145 - x209 - x210 - x211*x8 + 9.8994949366116653522*x33);result[20] = x212;
    }
    else if (x207)
    {
        DEBUG_PRINT("Case (x207) is True.\n");
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

        /* Assignment result[2, 6]=-x117*x219 + x118*x219 + x240*x26 - x240*x29 + 1 */
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
        x240=x112 - x20*x213;result[20] = -x117*x219 + x118*x219 + x240*x26 - x240*x29 + 1;
    }
    else if (x208)
    {
        DEBUG_PRINT("Case (x208) is True.\n");
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x112=x5*x58;
        x116=x112*x19;
        x117=x55*(rt2 + x116);
        x118=x61*(rt2 - x116);
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x220=random2*x122;

        /* Assignment result[2, 6]=-x117*x220 + x118*x220 + 1 */
        x55=(assert(IS_NOT_ZERO(x28)), 0.5/x28);
        x58=1.0/(assert(IS_NOT_ZERO(x22)), x22);
        x61=(assert(IS_NOT_ZERO(x25)), 0.5/x25);
        x112=x5*x58;
        x116=x112*x19;
        x117=x55*(rt2 + x116);
        x118=x61*(rt2 - x116);
        x122=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x220=random2*x122;result[20] = -x117*x220 + x118*x220 + 1;
    }
}
void fc3d_FischerBurmeisterFGenerated(
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
    }
    /* Assignment result[0, 0]=mu*rn - x10 - x11 + x3 */
    x4=mu*mu;
    x5=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x3*x3;
    x9=2*x8;
    x10=0.5*(assert(IS_POSITIVE(x5 - x9)), sqrt(x5 - x9));
    x11=0.5*(assert(IS_POSITIVE(x5 + x9)), sqrt(x5 + x9));result[0] = mu*rn - x10 - x11 + x3;

    /* Assignment result[1, 0]=Piecewise((x10*x14 - x11*x14 + x12, x15), (x10*x17 - x11*x17 + x12, x18)) */
    if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");
        x12=mu*ut1 + rt1;
        x13=1.0/(assert(IS_NOT_ZERO(x8)), x8);
        x14=x13*x6;

        /* Assignment result[1, 0]=x10*x14 - x11*x14 + x12 */
        x12=mu*ut1 + rt1;
        x13=1.0/(assert(IS_NOT_ZERO(x8)), x8);
        x14=x13*x6;result[1] = x10*x14 - x11*x14 + x12;
    }
    else if (x18)
    {
        DEBUG_PRINT("Case (x18) is True.\n");
        x12=mu*ut1 + rt1;
        x16=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x17=random1*x16;

        /* Assignment result[1, 0]=x10*x17 - x11*x17 + x12 */
        x12=mu*ut1 + rt1;
        x16=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x17=random1*x16;result[1] = x10*x17 - x11*x17 + x12;
    }

    /* Assignment result[2, 0]=Piecewise((x10*x20 - x11*x20 + x19, x15), (x10*x21 - x11*x21 + x19, x18)) */
    double x19;
    double x20;
    double x21;if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");
        x13=1.0/(assert(IS_NOT_ZERO(x8)), x8);
        x19=mu*ut2 + rt2;
        x20=x13*x7;

        /* Assignment result[2, 0]=x10*x20 - x11*x20 + x19 */
        x13=1.0/(assert(IS_NOT_ZERO(x8)), x8);
        x19=mu*ut2 + rt2;
        x20=x13*x7;result[2] = x10*x20 - x11*x20 + x19;
    }
    else if (x18)
    {
        DEBUG_PRINT("Case (x18) is True.\n");
        x16=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x19=mu*ut2 + rt2;
        x21=random2*x16;

        /* Assignment result[2, 0]=x10*x21 - x11*x21 + x19 */
        x16=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x19=mu*ut2 + rt2;
        x21=random2*x16;result[2] = x10*x21 - x11*x21 + x19;
    }
}
void fc3d_FischerBurmeisterABGenerated(
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
    double x5;
    double x6;
    double x7;
    double x8;
    double x9;
    int x10;
    x1=mu*mu;
    x2=ut1*ut1;
    x3=x1*x2;
    x4=ut2*ut2;
    x5=x1*x4;
    x6=(assert(IS_POSITIVE(x2 + x4)), sqrt(x2 + x4));
    x7=mu*x6 + un;
    x8=x7*x7;
    x9=(assert(IS_POSITIVE(x3 + x5 + x8)), sqrt(x3 + x5 + x8));
    x10=x9 <= 0;
    double x20;
    double x21;
    double x22;
    double x23;
    double x24;
    double x25;
    double x26;
    double x27;
    double x28;
    double x29;
    double x30;
    double x31;
    double x32;
    double x33;
    double x34;
    double x35;
    int x36;
    double x11;
    double x12;
    double x13;
    double x14;
    double x15;
    double x16;
    double x17;
    double x18;
    double x19;
    x20=mu*x7;
    x21=mu*rn*rt1 + ut1*x20;
    x22=x21*x21;
    x23=mu*rn*rt2 + ut2*x20;
    x24=x23*x23;
    x25=x22 + x24;
    x26=(assert(IS_POSITIVE(x25)), sqrt(x25));
    x27=rt1*rt1;
    x28=rt2*rt2;
    x29=rn*rn;
    x30=x1*x29;
    x31=x27 + x28 + x3 + x30 + x5 + x8;
    x32=2*x26;
    x33=x31 + x32;
    x34=x31 - x32;
    x35=fabs(x34);
    x36=x26 <= 0 || x33 <= 0 || x35 <= 0;
    int x49;
    double x37;
    double x38;
    double x39;
    double x40;
    double x41;
    double x42;
    double x43;
    double x44;
    double x45;
    double x46;
    double x47;
    double x48;
    x49=x6 <= 0;
    int x62;
    int x63;
    int x64;
    int x65;
    int x66;
    int x67;
    double x50;
    double x51;
    double x52;
    double x53;
    double x54;
    double x55;
    double x56;
    double x57;
    double x58;
    double x59;
    double x60;
    double x61;
    x62=x6 > 0;
    x63=x9 > 0;
    x64=x26 > 0;
    x65=x33 > 0;
    x66=x35 > 0;
    x67=x62 && x63 && x64 && x65 && x66;
    double x98;
    int x99;
    x38=x27 + x28;
    x98=(assert(IS_POSITIVE(x30 + x38)), sqrt(x30 + x38));
    x99=x98 <= 0;
    int x106;
    int x107;
    double x102;
    double x103;
    double x104;
    double x105;
    x106=x98 > 0;
    x107=x106 && x64 && x65 && x66;
    int x130;
    int x131;
    double x126;
    double x127;
    double x128;
    double x129;
    x130=x26 > 0;
    x131=x130 && x62 && x63 && x64 && x65 && x66;
    int x134;
    int x135;
    double x132;
    double x133;
    x134=x26 <= 0;
    x135=x134 && x62 && x63 && x64 && x65 && x66;
    int x204;
    double x202;
    double x203;
    x204=x106 && x130 && x64 && x65 && x66;
    int x205;
    x205=x106 && x134 && x64 && x65 && x66;
    /* Assignment result[0, 0]=Piecewise((1.00000000000000, x10), (x14*(-mu*x18 + x17 + x19), x36), (x42*x44*(-x46 - x47 + x48), x49), (-x58 - x61 + 1, x67)) */
    if (x10)
    {
        DEBUG_PRINT("Case (x10) is True.\n");

        /* Assignment result[0, 0]=1.00000000000000 */
        result[0] = 1.0;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x15=0.5*x11;
        x16=0.5*x12;
        x17=-x15 - x16;
        x18=1.4142135623730951455*x11;
        x19=x11*x12;

        /* Assignment result[0, 0]=x14*(-mu*x18 + x17 + x19) */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x15=0.5*x11;
        x16=0.5*x12;
        x17=-x15 - x16;
        x18=1.4142135623730951455*x11;
        x19=x11*x12;result[0] = x14*(-mu*x18 + x17 + x19);
    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x42=1.0/(assert(IS_NOT_ZERO(x41)), x41);
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x44=1.0/(assert(IS_NOT_ZERO(x43)), x43);
        x45=0.5*un;
        x46=x41*x45;
        x47=x43*x45;
        x48=x41*x43;

        /* Assignment result[0, 0]=x42*x44*(-x46 - x47 + x48) */
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x42=1.0/(assert(IS_NOT_ZERO(x41)), x41);
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x44=1.0/(assert(IS_NOT_ZERO(x43)), x43);
        x45=0.5*un;
        x46=x41*x45;
        x47=x43*x45;
        x48=x41*x43;result[0] = x42*x44*(-x46 - x47 + x48);
    }
    else if (x67)
    {
        DEBUG_PRINT("Case (x67) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x52=mu*x21;
        x53=ut1*x52;
        x54=mu*x23;
        x55=ut2*x54;
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x57=x56*(x53 + x55);
        x58=x51*(x57 + x7);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x61=x60*(-x57 + x7);

        /* Assignment result[0, 0]=-x58 - x61 + 1 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x52=mu*x21;
        x53=ut1*x52;
        x54=mu*x23;
        x55=ut2*x54;
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x57=x56*(x53 + x55);
        x58=x51*(x57 + x7);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x61=x60*(-x57 + x7);result[0] = -x58 - x61 + 1;
    }

    /* Assignment result[1, 0]=Piecewise((0.0, x10), (x118, x36), (x120*(-x124*x46 + x124*x47 - x125*x46 + x125*x47), x49), (-x128*x50 + x128*x59 - x129*x58 + x129*x61, x131), (-x133*x58 + x133*x61, x135)) */
    double x71;
    double x100;
    double x108;
    double x109;
    double x118;
    double x119;
    double x120;
    double x121;
    double x122;
    double x123;
    double x124;
    double x125;if (x10)
    {
        DEBUG_PRINT("Case (x10) is True.\n");

        /* Assignment result[1, 0]=0.0 */
        result[1] = 0.0;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x100=mu*x11;
        x108=0.35355339059327378637*x11;
        x109=0.35355339059327378637*x12;
        x118=x14*(-1.0*x100 - x108 + x109);

        /* Assignment result[1, 0]=x118 */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x100=mu*x11;
        x108=0.35355339059327378637*x11;
        x109=0.35355339059327378637*x12;
        x118=x14*(-1.0*x100 - x108 + x109);result[1] = x118;
    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x42=1.0/(assert(IS_NOT_ZERO(x41)), x41);
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x44=1.0/(assert(IS_NOT_ZERO(x43)), x43);
        x45=0.5*un;
        x46=x41*x45;
        x47=x43*x45;
        x71=x27*x30 + x28*x30;
        x119=pow(x71, 3.0/2.0);
        x120=(assert(IS_NOT_ZERO(x119)), x42*x44/x119);
        x121=pow(mu, 3);
        x122=pow(rn, 3);
        x123=pow(rt1, 3);
        x124=x121*x122*x123;
        x125=rt1*x121*x122*x28;

        /* Assignment result[1, 0]=x120*(-x124*x46 + x124*x47 - x125*x46 + x125*x47) */
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x42=1.0/(assert(IS_NOT_ZERO(x41)), x41);
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x44=1.0/(assert(IS_NOT_ZERO(x43)), x43);
        x45=0.5*un;
        x46=x41*x45;
        x47=x43*x45;
        x71=x27*x30 + x28*x30;
        x119=pow(x71, 3.0/2.0);
        x120=(assert(IS_NOT_ZERO(x119)), x42*x44/x119);
        x121=pow(mu, 3);
        x122=pow(rn, 3);
        x123=pow(rt1, 3);
        x124=x121*x122*x123;
        x125=rt1*x121*x122*x28;result[1] = x120*(-x124*x46 + x124*x47 - x125*x46 + x125*x47);
    }
    else if (x131)
    {
        DEBUG_PRINT("Case (x131) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x52=mu*x21;
        x53=ut1*x52;
        x54=mu*x23;
        x55=ut2*x54;
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x57=x56*(x53 + x55);
        x58=x51*(x57 + x7);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x61=x60*(-x57 + x7);
        x126=pow(x25, -3.0/2.0);
        x127=-x53 - x55;
        x128=0.5*mu*ut1*x56 + 0.5*x126*x127*x21;
        x129=x21*x56;

        /* Assignment result[1, 0]=-x128*x50 + x128*x59 - x129*x58 + x129*x61 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x52=mu*x21;
        x53=ut1*x52;
        x54=mu*x23;
        x55=ut2*x54;
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x57=x56*(x53 + x55);
        x58=x51*(x57 + x7);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x61=x60*(-x57 + x7);
        x126=pow(x25, -3.0/2.0);
        x127=-x53 - x55;
        x128=0.5*mu*ut1*x56 + 0.5*x126*x127*x21;
        x129=x21*x56;result[1] = -x128*x50 + x128*x59 - x129*x58 + x129*x61;
    }
    else if (x135)
    {
        DEBUG_PRINT("Case (x135) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x52=mu*x21;
        x53=ut1*x52;
        x54=mu*x23;
        x55=ut2*x54;
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x57=x56*(x53 + x55);
        x58=x51*(x57 + x7);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x61=x60*(-x57 + x7);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x133=random1*x132;

        /* Assignment result[1, 0]=-x133*x58 + x133*x61 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x52=mu*x21;
        x53=ut1*x52;
        x54=mu*x23;
        x55=ut2*x54;
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x57=x56*(x53 + x55);
        x58=x51*(x57 + x7);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x61=x60*(-x57 + x7);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x133=random1*x132;result[1] = -x133*x58 + x133*x61;
    }

    /* Assignment result[2, 0]=Piecewise((0.0, x10), (x118, x36), (x120*(-x215*x46 + x215*x47 - x216*x46 + x216*x47), x49), (-x217*x50 + x217*x59 - x218*x58 + x218*x61, x131), (-x219*x58 + x219*x61, x135)) */
    double x190;
    double x215;
    double x216;
    double x217;
    double x218;
    double x219;if (x10)
    {
        DEBUG_PRINT("Case (x10) is True.\n");

        /* Assignment result[2, 0]=0.0 */
        result[2] = 0.0;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x100=mu*x11;
        x108=0.35355339059327378637*x11;
        x109=0.35355339059327378637*x12;
        x118=x14*(-1.0*x100 - x108 + x109);

        /* Assignment result[2, 0]=x118 */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x100=mu*x11;
        x108=0.35355339059327378637*x11;
        x109=0.35355339059327378637*x12;
        x118=x14*(-1.0*x100 - x108 + x109);result[2] = x118;
    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x42=1.0/(assert(IS_NOT_ZERO(x41)), x41);
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x44=1.0/(assert(IS_NOT_ZERO(x43)), x43);
        x45=0.5*un;
        x46=x41*x45;
        x47=x43*x45;
        x71=x27*x30 + x28*x30;
        x119=pow(x71, 3.0/2.0);
        x120=(assert(IS_NOT_ZERO(x119)), x42*x44/x119);
        x121=pow(mu, 3);
        x122=pow(rn, 3);
        x190=pow(rt2, 3);
        x215=x121*x122*x190;
        x216=rt2*x121*x122*x27;

        /* Assignment result[2, 0]=x120*(-x215*x46 + x215*x47 - x216*x46 + x216*x47) */
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x42=1.0/(assert(IS_NOT_ZERO(x41)), x41);
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x44=1.0/(assert(IS_NOT_ZERO(x43)), x43);
        x45=0.5*un;
        x46=x41*x45;
        x47=x43*x45;
        x71=x27*x30 + x28*x30;
        x119=pow(x71, 3.0/2.0);
        x120=(assert(IS_NOT_ZERO(x119)), x42*x44/x119);
        x121=pow(mu, 3);
        x122=pow(rn, 3);
        x190=pow(rt2, 3);
        x215=x121*x122*x190;
        x216=rt2*x121*x122*x27;result[2] = x120*(-x215*x46 + x215*x47 - x216*x46 + x216*x47);
    }
    else if (x131)
    {
        DEBUG_PRINT("Case (x131) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x52=mu*x21;
        x53=ut1*x52;
        x54=mu*x23;
        x55=ut2*x54;
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x57=x56*(x53 + x55);
        x58=x51*(x57 + x7);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x61=x60*(-x57 + x7);
        x126=pow(x25, -3.0/2.0);
        x127=-x53 - x55;
        x217=0.5*mu*ut2*x56 + 0.5*x126*x127*x23;
        x218=x23*x56;

        /* Assignment result[2, 0]=-x217*x50 + x217*x59 - x218*x58 + x218*x61 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x52=mu*x21;
        x53=ut1*x52;
        x54=mu*x23;
        x55=ut2*x54;
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x57=x56*(x53 + x55);
        x58=x51*(x57 + x7);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x61=x60*(-x57 + x7);
        x126=pow(x25, -3.0/2.0);
        x127=-x53 - x55;
        x217=0.5*mu*ut2*x56 + 0.5*x126*x127*x23;
        x218=x23*x56;result[2] = -x217*x50 + x217*x59 - x218*x58 + x218*x61;
    }
    else if (x135)
    {
        DEBUG_PRINT("Case (x135) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x52=mu*x21;
        x53=ut1*x52;
        x54=mu*x23;
        x55=ut2*x54;
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x57=x56*(x53 + x55);
        x58=x51*(x57 + x7);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x61=x60*(-x57 + x7);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x219=random2*x132;

        /* Assignment result[2, 0]=-x219*x58 + x219*x61 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x52=mu*x21;
        x53=ut1*x52;
        x54=mu*x23;
        x55=ut2*x54;
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x57=x56*(x53 + x55);
        x58=x51*(x57 + x7);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x61=x60*(-x57 + x7);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x219=random2*x132;result[2] = -x219*x58 + x219*x61;
    }

    /* Assignment result[0, 1]=Piecewise((x68, x10), (x70, x36), (x73*(-x46*x74 + x47*x74 + x76), x49), (x79 - x87 - x88, x67)) */
    double x68;
    double x69;
    double x70;
    double x72;
    double x73;
    double x74;
    double x75;
    double x76;
    double x77;
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
    double x88;if (x10)
    {
        DEBUG_PRINT("Case (x10) is True.\n");
        x68=0.70710678118654757274*mu;

        /* Assignment result[0, 1]=x68 */
        x68=0.70710678118654757274*mu;result[3] = x68;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x68=0.70710678118654757274*mu;
        x69=x12*x68;
        x70=x13*(-2.0*x1 - x68 + x69);

        /* Assignment result[0, 1]=x70 */
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x68=0.70710678118654757274*mu;
        x69=x12*x68;
        x70=x13*(-2.0*x1 - x68 + x69);result[3] = x70;
    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x42=1.0/(assert(IS_NOT_ZERO(x41)), x41);
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x44=1.0/(assert(IS_NOT_ZERO(x43)), x43);
        x45=0.5*un;
        x46=x41*x45;
        x47=x43*x45;
        x48=x41*x43;
        x71=x27*x30 + x28*x30;
        x72=(assert(IS_POSITIVE(x71)), sqrt(x71));
        x73=(assert(IS_NOT_ZERO(x72)), x42*x44/x72);
        x74=rn*rt1*x1;
        x75=0.35355339059327378637*mu*un*x72;
        x76=0.70710678118654757274*mu*x48*x72 - x41*x75 - x43*x75;

        /* Assignment result[0, 1]=x73*(-x46*x74 + x47*x74 + x76) */
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x42=1.0/(assert(IS_NOT_ZERO(x41)), x41);
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x44=1.0/(assert(IS_NOT_ZERO(x43)), x43);
        x45=0.5*un;
        x46=x41*x45;
        x47=x43*x45;
        x48=x41*x43;
        x71=x27*x30 + x28*x30;
        x72=(assert(IS_POSITIVE(x71)), sqrt(x71));
        x73=(assert(IS_NOT_ZERO(x72)), x42*x44/x72);
        x74=rn*rt1*x1;
        x75=0.35355339059327378637*mu*un*x72;
        x76=0.70710678118654757274*mu*x48*x72 - x41*x75 - x43*x75;result[3] = x73*(-x46*x74 + x47*x74 + x76);
    }
    else if (x67)
    {
        DEBUG_PRINT("Case (x67) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x79=ut1*x78;
        x80=ut1*x1 + x7*x79;
        x81=ut1*ut2*x1*x77;
        x82=x23*x81;
        x83=2*x20;
        x84=x1*x2*x77;
        x85=(1.0/2.0)*x21*(x83 + 2*x84);
        x86=x56*(x82 + x85);
        x87=x51*(x80 + x86);
        x88=x60*(x80 - x86);

        /* Assignment result[0, 1]=x79 - x87 - x88 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x79=ut1*x78;
        x80=ut1*x1 + x7*x79;
        x81=ut1*ut2*x1*x77;
        x82=x23*x81;
        x83=2*x20;
        x84=x1*x2*x77;
        x85=(1.0/2.0)*x21*(x83 + 2*x84);
        x86=x56*(x82 + x85);
        x87=x51*(x80 + x86);
        x88=x60*(x80 - x86);result[3] = x79 - x87 - x88;
    }

    /* Assignment result[1, 1]=Piecewise((mu, x10), (x148, x36), (x152*(-x146*x175 + x146*x176 + x158 - x159*x41 + x159*x43 + x163 + x164 - x166*x41 + x166*x43 + x168 + x169 - x170*x41 + x170*x43 - x171*x172 + x171*x173 + x178 + x179 - 2.0*x180 - 2.0*x181 - 4.0*x182 - 4.0*x183), x49), (mu - x129*x87 + x129*x88 - x185*x50 + x185*x59, x131), (mu - x133*x87 + x133*x88, x135)) */
    double x101;
    double x136;
    double x137;
    double x138;
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
    double x185;if (x10)
    {
        DEBUG_PRINT("Case (x10) is True.\n");

        /* Assignment result[1, 1]=mu */
        result[4] = mu;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x121=pow(mu, 3);
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x138=4.2426406871192847703*x11;
        x139=4.2426406871192847703*x12;
        x140=mu*x12;
        x141=2.0*x140;
        x142=x1*x11;
        x143=x11*x121;
        x144=x1*x139;
        x145=-x144;
        x146=2.0*x121;
        x147=x12*x146;
        x148=x137*(-26.0*x100 + 16.0*x101 + 11.313708498984761164*x12*x142 - x138 + x139 + x141 - 57.982756057296896302*x142 - 34.0*x143 + x145 + x147);

        /* Assignment result[1, 1]=x148 */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x121=pow(mu, 3);
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x138=4.2426406871192847703*x11;
        x139=4.2426406871192847703*x12;
        x140=mu*x12;
        x141=2.0*x140;
        x142=x1*x11;
        x143=x11*x121;
        x144=x1*x139;
        x145=-x144;
        x146=2.0*x121;
        x147=x12*x146;
        x148=x137*(-26.0*x100 + 16.0*x101 + 11.313708498984761164*x12*x142 - x138 + x139 + x141 - 57.982756057296896302*x142 - 34.0*x143 + x145 + x147);result[4] = x148;
    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x71=x27*x30 + x28*x30;
        x119=pow(x71, 3.0/2.0);
        x121=pow(mu, 3);
        x122=pow(rn, 3);
        x123=pow(rt1, 3);
        x146=2.0*x121;
        x149=4.0*x119*x41*x43;
        x150=x149*x27;
        x151=x149*x28;
        x152=1.0/(assert(IS_NOT_ZERO((x150 + x151))), (x150 + x151));
        x153=pow(mu, 5);
        x154=pow(rn, 4);
        x155=2.0*un*x153*x154*x27*x28;
        x156=pow(un, 3);
        x157=2.0*x121*x156*x27*x28*x29;
        x158=mu*x150 + mu*x151 - x155*x41 + x155*x43 - x157*x41 + x157*x43;
        x159=2.0*pow(rt2, 6)*un*x121*x29;
        x160=pow(mu, 4);
        x161=pow(rt1, 5);
        x162=1.4142135623730951455*un*x122*x160*x161;
        x163=-x162*x41;
        x164=x162*x43;
        x165=pow(rt2, 4);
        x166=2.0*un*x153*x154*x165;
        x167=1.4142135623730951455*rt1*un*x122*x160*x165;
        x168=-x167*x41;
        x169=x167*x43;
        x170=2.0*x121*x156*x165*x29;
        x171=4.0*x121;
        x172=un*x165*x27*x29*x41;
        x173=un*x165*x27*x29*x43;
        x174=pow(rt1, 4);
        x175=un*x174*x28*x29*x41;
        x176=un*x174*x28*x29*x43;
        x177=2.8284271247461902909*un*x122*x123*x160*x28;
        x178=-x177*x41;
        x179=x177*x43;
        x180=mu*un*x119*x27*x41;
        x181=mu*un*x119*x27*x43;
        x182=mu*un*x119*x28*x41;
        x183=mu*un*x119*x28*x43;

        /* Assignment result[1, 1]=x152*(-x146*x175 + x146*x176 + x158 - x159*x41 + x159*x43 + x163 + x164 - x166*x41 + x166*x43 + x168 + x169 - x170*x41 + x170*x43 - x171*x172 + x171*x173 + x178 + x179 - 2.0*x180 - 2.0*x181 - 4.0*x182 - 4.0*x183) */
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x71=x27*x30 + x28*x30;
        x119=pow(x71, 3.0/2.0);
        x121=pow(mu, 3);
        x122=pow(rn, 3);
        x123=pow(rt1, 3);
        x146=2.0*x121;
        x149=4.0*x119*x41*x43;
        x150=x149*x27;
        x151=x149*x28;
        x152=1.0/(assert(IS_NOT_ZERO((x150 + x151))), (x150 + x151));
        x153=pow(mu, 5);
        x154=pow(rn, 4);
        x155=2.0*un*x153*x154*x27*x28;
        x156=pow(un, 3);
        x157=2.0*x121*x156*x27*x28*x29;
        x158=mu*x150 + mu*x151 - x155*x41 + x155*x43 - x157*x41 + x157*x43;
        x159=2.0*pow(rt2, 6)*un*x121*x29;
        x160=pow(mu, 4);
        x161=pow(rt1, 5);
        x162=1.4142135623730951455*un*x122*x160*x161;
        x163=-x162*x41;
        x164=x162*x43;
        x165=pow(rt2, 4);
        x166=2.0*un*x153*x154*x165;
        x167=1.4142135623730951455*rt1*un*x122*x160*x165;
        x168=-x167*x41;
        x169=x167*x43;
        x170=2.0*x121*x156*x165*x29;
        x171=4.0*x121;
        x172=un*x165*x27*x29*x41;
        x173=un*x165*x27*x29*x43;
        x174=pow(rt1, 4);
        x175=un*x174*x28*x29*x41;
        x176=un*x174*x28*x29*x43;
        x177=2.8284271247461902909*un*x122*x123*x160*x28;
        x178=-x177*x41;
        x179=x177*x43;
        x180=mu*un*x119*x27*x41;
        x181=mu*un*x119*x27*x43;
        x182=mu*un*x119*x28*x41;
        x183=mu*un*x119*x28*x43;result[4] = x152*(-x146*x175 + x146*x176 + x158 - x159*x41 + x159*x43 + x163 + x164 - x166*x41 + x166*x43 + x168 + x169 - x170*x41 + x170*x43 - x171*x172 + x171*x173 + x178 + x179 - 2.0*x180 - 2.0*x181 - 4.0*x182 - 4.0*x183);
    }
    else if (x131)
    {
        DEBUG_PRINT("Case (x131) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x79=ut1*x78;
        x80=ut1*x1 + x7*x79;
        x81=ut1*ut2*x1*x77;
        x82=x23*x81;
        x83=2*x20;
        x84=x1*x2*x77;
        x85=(1.0/2.0)*x21*(x83 + 2*x84);
        x86=x56*(x82 + x85);
        x87=x51*(x80 + x86);
        x88=x60*(x80 - x86);
        x126=pow(x25, -3.0/2.0);
        x129=x21*x56;
        x184=-x82 - x85;
        x185=0.5*x126*x184*x21 + 0.5*x56*(x20 + x84);

        /* Assignment result[1, 1]=mu - x129*x87 + x129*x88 - x185*x50 + x185*x59 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x79=ut1*x78;
        x80=ut1*x1 + x7*x79;
        x81=ut1*ut2*x1*x77;
        x82=x23*x81;
        x83=2*x20;
        x84=x1*x2*x77;
        x85=(1.0/2.0)*x21*(x83 + 2*x84);
        x86=x56*(x82 + x85);
        x87=x51*(x80 + x86);
        x88=x60*(x80 - x86);
        x126=pow(x25, -3.0/2.0);
        x129=x21*x56;
        x184=-x82 - x85;
        x185=0.5*x126*x184*x21 + 0.5*x56*(x20 + x84);result[4] = mu - x129*x87 + x129*x88 - x185*x50 + x185*x59;
    }
    else if (x135)
    {
        DEBUG_PRINT("Case (x135) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x79=ut1*x78;
        x80=ut1*x1 + x7*x79;
        x81=ut1*ut2*x1*x77;
        x82=x23*x81;
        x83=2*x20;
        x84=x1*x2*x77;
        x85=(1.0/2.0)*x21*(x83 + 2*x84);
        x86=x56*(x82 + x85);
        x87=x51*(x80 + x86);
        x88=x60*(x80 - x86);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x133=random1*x132;

        /* Assignment result[1, 1]=mu - x133*x87 + x133*x88 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x79=ut1*x78;
        x80=ut1*x1 + x7*x79;
        x81=ut1*ut2*x1*x77;
        x82=x23*x81;
        x83=2*x20;
        x84=x1*x2*x77;
        x85=(1.0/2.0)*x21*(x83 + 2*x84);
        x86=x56*(x82 + x85);
        x87=x51*(x80 + x86);
        x88=x60*(x80 - x86);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x133=random1*x132;result[4] = mu - x133*x87 + x133*x88;
    }

    /* Assignment result[2, 1]=Piecewise((0.0, x10), (x186, x36), (x152*(x197 + x221 + x222 + x224 + x225 + x227 + x228), x49), (-x218*x87 + x218*x88 - x229*x50 + x229*x59, x131), (-x219*x87 + x219*x88, x135)) */
    double x186;
    double x187;
    double x188;
    double x189;
    double x191;
    double x192;
    double x193;
    double x194;
    double x195;
    double x196;
    double x197;
    double x198;
    double x220;
    double x221;
    double x222;
    double x223;
    double x224;
    double x225;
    double x226;
    double x227;
    double x228;
    double x229;if (x10)
    {
        DEBUG_PRINT("Case (x10) is True.\n");

        /* Assignment result[2, 1]=0.0 */
        result[5] = 0.0;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x18=1.4142135623730951455*x11;
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x121=pow(mu, 3);
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x138=4.2426406871192847703*x11;
        x139=4.2426406871192847703*x12;
        x140=mu*x12;
        x141=2.0*x140;
        x143=x11*x121;
        x144=x1*x139;
        x146=2.0*x121;
        x147=x12*x146;
        x186=x137*(x1*x18 + 10.0*x100 + x138 - x139 - x141 + 2.0*x143 + x144 - x147);

        /* Assignment result[2, 1]=x186 */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x18=1.4142135623730951455*x11;
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x121=pow(mu, 3);
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x138=4.2426406871192847703*x11;
        x139=4.2426406871192847703*x12;
        x140=mu*x12;
        x141=2.0*x140;
        x143=x11*x121;
        x144=x1*x139;
        x146=2.0*x121;
        x147=x12*x146;
        x186=x137*(x1*x18 + 10.0*x100 + x138 - x139 - x141 + 2.0*x143 + x144 - x147);result[5] = x186;
    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x71=x27*x30 + x28*x30;
        x119=pow(x71, 3.0/2.0);
        x121=pow(mu, 3);
        x122=pow(rn, 3);
        x123=pow(rt1, 3);
        x149=4.0*x119*x41*x43;
        x150=x149*x27;
        x151=x149*x28;
        x152=1.0/(assert(IS_NOT_ZERO((x150 + x151))), (x150 + x151));
        x153=pow(mu, 5);
        x154=pow(rn, 4);
        x156=pow(un, 3);
        x160=pow(mu, 4);
        x161=pow(rt1, 5);
        x174=pow(rt1, 4);
        x187=pow(rt2, 5);
        x188=2.0*rt1*un*x121*x187*x29;
        x189=2.0*rt2*un*x121*x161*x29;
        x190=pow(rt2, 3);
        x191=2.0*rt1*un*x153*x154*x190;
        x192=2.0*rt2*un*x123*x153*x154;
        x193=2.0*rt1*x121*x156*x190*x29;
        x194=2.0*rt2*x121*x123*x156*x29;
        x195=4.0*un*x121*x123*x190*x29;
        x196=2.0*mu*rt1*rt2*un*x119;
        x197=x188*x41 - x188*x43 + x189*x41 - x189*x43 + x191*x41 - x191*x43 + x192*x41 - x192*x43 + x193*x41 - x193*x43 + x194*x41 - x194*x43 + x195*x41 - x195*x43 + x196*x41 + x196*x43;
        x220=1.4142135623730951455*un*x122*x160*x187;
        x221=-x220*x41;
        x222=x220*x43;
        x223=1.4142135623730951455*rt2*un*x122*x160*x174;
        x224=-x223*x41;
        x225=x223*x43;
        x226=2.8284271247461902909*un*x122*x160*x190*x27;
        x227=-x226*x41;
        x228=x226*x43;

        /* Assignment result[2, 1]=x152*(x197 + x221 + x222 + x224 + x225 + x227 + x228) */
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x71=x27*x30 + x28*x30;
        x119=pow(x71, 3.0/2.0);
        x121=pow(mu, 3);
        x122=pow(rn, 3);
        x123=pow(rt1, 3);
        x149=4.0*x119*x41*x43;
        x150=x149*x27;
        x151=x149*x28;
        x152=1.0/(assert(IS_NOT_ZERO((x150 + x151))), (x150 + x151));
        x153=pow(mu, 5);
        x154=pow(rn, 4);
        x156=pow(un, 3);
        x160=pow(mu, 4);
        x161=pow(rt1, 5);
        x174=pow(rt1, 4);
        x187=pow(rt2, 5);
        x188=2.0*rt1*un*x121*x187*x29;
        x189=2.0*rt2*un*x121*x161*x29;
        x190=pow(rt2, 3);
        x191=2.0*rt1*un*x153*x154*x190;
        x192=2.0*rt2*un*x123*x153*x154;
        x193=2.0*rt1*x121*x156*x190*x29;
        x194=2.0*rt2*x121*x123*x156*x29;
        x195=4.0*un*x121*x123*x190*x29;
        x196=2.0*mu*rt1*rt2*un*x119;
        x197=x188*x41 - x188*x43 + x189*x41 - x189*x43 + x191*x41 - x191*x43 + x192*x41 - x192*x43 + x193*x41 - x193*x43 + x194*x41 - x194*x43 + x195*x41 - x195*x43 + x196*x41 + x196*x43;
        x220=1.4142135623730951455*un*x122*x160*x187;
        x221=-x220*x41;
        x222=x220*x43;
        x223=1.4142135623730951455*rt2*un*x122*x160*x174;
        x224=-x223*x41;
        x225=x223*x43;
        x226=2.8284271247461902909*un*x122*x160*x190*x27;
        x227=-x226*x41;
        x228=x226*x43;result[5] = x152*(x197 + x221 + x222 + x224 + x225 + x227 + x228);
    }
    else if (x131)
    {
        DEBUG_PRINT("Case (x131) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x79=ut1*x78;
        x80=ut1*x1 + x7*x79;
        x81=ut1*ut2*x1*x77;
        x82=x23*x81;
        x83=2*x20;
        x84=x1*x2*x77;
        x85=(1.0/2.0)*x21*(x83 + 2*x84);
        x86=x56*(x82 + x85);
        x87=x51*(x80 + x86);
        x88=x60*(x80 - x86);
        x126=pow(x25, -3.0/2.0);
        x184=-x82 - x85;
        x198=0.5*ut1*ut2*x1*x56*x77;
        x218=x23*x56;
        x229=0.5*x126*x184*x23 + x198;

        /* Assignment result[2, 1]=-x218*x87 + x218*x88 - x229*x50 + x229*x59 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x79=ut1*x78;
        x80=ut1*x1 + x7*x79;
        x81=ut1*ut2*x1*x77;
        x82=x23*x81;
        x83=2*x20;
        x84=x1*x2*x77;
        x85=(1.0/2.0)*x21*(x83 + 2*x84);
        x86=x56*(x82 + x85);
        x87=x51*(x80 + x86);
        x88=x60*(x80 - x86);
        x126=pow(x25, -3.0/2.0);
        x184=-x82 - x85;
        x198=0.5*ut1*ut2*x1*x56*x77;
        x218=x23*x56;
        x229=0.5*x126*x184*x23 + x198;result[5] = -x218*x87 + x218*x88 - x229*x50 + x229*x59;
    }
    else if (x135)
    {
        DEBUG_PRINT("Case (x135) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x79=ut1*x78;
        x80=ut1*x1 + x7*x79;
        x81=ut1*ut2*x1*x77;
        x82=x23*x81;
        x83=2*x20;
        x84=x1*x2*x77;
        x85=(1.0/2.0)*x21*(x83 + 2*x84);
        x86=x56*(x82 + x85);
        x87=x51*(x80 + x86);
        x88=x60*(x80 - x86);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x219=random2*x132;

        /* Assignment result[2, 1]=-x219*x87 + x219*x88 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x79=ut1*x78;
        x80=ut1*x1 + x7*x79;
        x81=ut1*ut2*x1*x77;
        x82=x23*x81;
        x83=2*x20;
        x84=x1*x2*x77;
        x85=(1.0/2.0)*x21*(x83 + 2*x84);
        x86=x56*(x82 + x85);
        x87=x51*(x80 + x86);
        x88=x60*(x80 - x86);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x219=random2*x132;result[5] = -x219*x87 + x219*x88;
    }

    /* Assignment result[0, 2]=Piecewise((x68, x10), (x70, x36), (x73*(-x46*x89 + x47*x89 + x76), x49), (x90 - x96 - x97, x67)) */
    double x89;
    double x90;
    double x91;
    double x92;
    double x93;
    double x94;
    double x95;
    double x96;
    double x97;if (x10)
    {
        DEBUG_PRINT("Case (x10) is True.\n");
        x68=0.70710678118654757274*mu;

        /* Assignment result[0, 2]=x68 */
        x68=0.70710678118654757274*mu;result[6] = x68;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x68=0.70710678118654757274*mu;
        x69=x12*x68;
        x70=x13*(-2.0*x1 - x68 + x69);

        /* Assignment result[0, 2]=x70 */
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x68=0.70710678118654757274*mu;
        x69=x12*x68;
        x70=x13*(-2.0*x1 - x68 + x69);result[6] = x70;
    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x42=1.0/(assert(IS_NOT_ZERO(x41)), x41);
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x44=1.0/(assert(IS_NOT_ZERO(x43)), x43);
        x45=0.5*un;
        x46=x41*x45;
        x47=x43*x45;
        x48=x41*x43;
        x71=x27*x30 + x28*x30;
        x72=(assert(IS_POSITIVE(x71)), sqrt(x71));
        x73=(assert(IS_NOT_ZERO(x72)), x42*x44/x72);
        x75=0.35355339059327378637*mu*un*x72;
        x76=0.70710678118654757274*mu*x48*x72 - x41*x75 - x43*x75;
        x89=rn*rt2*x1;

        /* Assignment result[0, 2]=x73*(-x46*x89 + x47*x89 + x76) */
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x42=1.0/(assert(IS_NOT_ZERO(x41)), x41);
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x44=1.0/(assert(IS_NOT_ZERO(x43)), x43);
        x45=0.5*un;
        x46=x41*x45;
        x47=x43*x45;
        x48=x41*x43;
        x71=x27*x30 + x28*x30;
        x72=(assert(IS_POSITIVE(x71)), sqrt(x71));
        x73=(assert(IS_NOT_ZERO(x72)), x42*x44/x72);
        x75=0.35355339059327378637*mu*un*x72;
        x76=0.70710678118654757274*mu*x48*x72 - x41*x75 - x43*x75;
        x89=rn*rt2*x1;result[6] = x73*(-x46*x89 + x47*x89 + x76);
    }
    else if (x67)
    {
        DEBUG_PRINT("Case (x67) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x81=ut1*ut2*x1*x77;
        x83=2*x20;
        x90=ut2*x78;
        x91=ut2*x1 + x7*x90;
        x92=x21*x81;
        x93=x1*x4*x77;
        x94=(1.0/2.0)*x23*(x83 + 2*x93);
        x95=x56*(x92 + x94);
        x96=x51*(x91 + x95);
        x97=x60*(x91 - x95);

        /* Assignment result[0, 2]=x90 - x96 - x97 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x81=ut1*ut2*x1*x77;
        x83=2*x20;
        x90=ut2*x78;
        x91=ut2*x1 + x7*x90;
        x92=x21*x81;
        x93=x1*x4*x77;
        x94=(1.0/2.0)*x23*(x83 + 2*x93);
        x95=x56*(x92 + x94);
        x96=x51*(x91 + x95);
        x97=x60*(x91 - x95);result[6] = x90 - x96 - x97;
    }

    /* Assignment result[1, 2]=Piecewise((0.0, x10), (x186, x36), (x152*(x163 + x164 + x168 + x169 + x178 + x179 + x197), x49), (-x129*x96 + x129*x97 - x200*x50 + x200*x59, x131), (-x133*x96 + x133*x97, x135)) */
    double x199;
    double x200;if (x10)
    {
        DEBUG_PRINT("Case (x10) is True.\n");

        /* Assignment result[1, 2]=0.0 */
        result[7] = 0.0;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x18=1.4142135623730951455*x11;
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x121=pow(mu, 3);
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x138=4.2426406871192847703*x11;
        x139=4.2426406871192847703*x12;
        x140=mu*x12;
        x141=2.0*x140;
        x143=x11*x121;
        x144=x1*x139;
        x146=2.0*x121;
        x147=x12*x146;
        x186=x137*(x1*x18 + 10.0*x100 + x138 - x139 - x141 + 2.0*x143 + x144 - x147);

        /* Assignment result[1, 2]=x186 */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x18=1.4142135623730951455*x11;
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x121=pow(mu, 3);
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x138=4.2426406871192847703*x11;
        x139=4.2426406871192847703*x12;
        x140=mu*x12;
        x141=2.0*x140;
        x143=x11*x121;
        x144=x1*x139;
        x146=2.0*x121;
        x147=x12*x146;
        x186=x137*(x1*x18 + 10.0*x100 + x138 - x139 - x141 + 2.0*x143 + x144 - x147);result[7] = x186;
    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x71=x27*x30 + x28*x30;
        x119=pow(x71, 3.0/2.0);
        x121=pow(mu, 3);
        x122=pow(rn, 3);
        x123=pow(rt1, 3);
        x149=4.0*x119*x41*x43;
        x150=x149*x27;
        x151=x149*x28;
        x152=1.0/(assert(IS_NOT_ZERO((x150 + x151))), (x150 + x151));
        x153=pow(mu, 5);
        x154=pow(rn, 4);
        x156=pow(un, 3);
        x160=pow(mu, 4);
        x161=pow(rt1, 5);
        x162=1.4142135623730951455*un*x122*x160*x161;
        x163=-x162*x41;
        x164=x162*x43;
        x165=pow(rt2, 4);
        x167=1.4142135623730951455*rt1*un*x122*x160*x165;
        x168=-x167*x41;
        x169=x167*x43;
        x177=2.8284271247461902909*un*x122*x123*x160*x28;
        x178=-x177*x41;
        x179=x177*x43;
        x187=pow(rt2, 5);
        x188=2.0*rt1*un*x121*x187*x29;
        x189=2.0*rt2*un*x121*x161*x29;
        x190=pow(rt2, 3);
        x191=2.0*rt1*un*x153*x154*x190;
        x192=2.0*rt2*un*x123*x153*x154;
        x193=2.0*rt1*x121*x156*x190*x29;
        x194=2.0*rt2*x121*x123*x156*x29;
        x195=4.0*un*x121*x123*x190*x29;
        x196=2.0*mu*rt1*rt2*un*x119;
        x197=x188*x41 - x188*x43 + x189*x41 - x189*x43 + x191*x41 - x191*x43 + x192*x41 - x192*x43 + x193*x41 - x193*x43 + x194*x41 - x194*x43 + x195*x41 - x195*x43 + x196*x41 + x196*x43;

        /* Assignment result[1, 2]=x152*(x163 + x164 + x168 + x169 + x178 + x179 + x197) */
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x71=x27*x30 + x28*x30;
        x119=pow(x71, 3.0/2.0);
        x121=pow(mu, 3);
        x122=pow(rn, 3);
        x123=pow(rt1, 3);
        x149=4.0*x119*x41*x43;
        x150=x149*x27;
        x151=x149*x28;
        x152=1.0/(assert(IS_NOT_ZERO((x150 + x151))), (x150 + x151));
        x153=pow(mu, 5);
        x154=pow(rn, 4);
        x156=pow(un, 3);
        x160=pow(mu, 4);
        x161=pow(rt1, 5);
        x162=1.4142135623730951455*un*x122*x160*x161;
        x163=-x162*x41;
        x164=x162*x43;
        x165=pow(rt2, 4);
        x167=1.4142135623730951455*rt1*un*x122*x160*x165;
        x168=-x167*x41;
        x169=x167*x43;
        x177=2.8284271247461902909*un*x122*x123*x160*x28;
        x178=-x177*x41;
        x179=x177*x43;
        x187=pow(rt2, 5);
        x188=2.0*rt1*un*x121*x187*x29;
        x189=2.0*rt2*un*x121*x161*x29;
        x190=pow(rt2, 3);
        x191=2.0*rt1*un*x153*x154*x190;
        x192=2.0*rt2*un*x123*x153*x154;
        x193=2.0*rt1*x121*x156*x190*x29;
        x194=2.0*rt2*x121*x123*x156*x29;
        x195=4.0*un*x121*x123*x190*x29;
        x196=2.0*mu*rt1*rt2*un*x119;
        x197=x188*x41 - x188*x43 + x189*x41 - x189*x43 + x191*x41 - x191*x43 + x192*x41 - x192*x43 + x193*x41 - x193*x43 + x194*x41 - x194*x43 + x195*x41 - x195*x43 + x196*x41 + x196*x43;result[7] = x152*(x163 + x164 + x168 + x169 + x178 + x179 + x197);
    }
    else if (x131)
    {
        DEBUG_PRINT("Case (x131) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x81=ut1*ut2*x1*x77;
        x83=2*x20;
        x90=ut2*x78;
        x91=ut2*x1 + x7*x90;
        x92=x21*x81;
        x93=x1*x4*x77;
        x94=(1.0/2.0)*x23*(x83 + 2*x93);
        x95=x56*(x92 + x94);
        x96=x51*(x91 + x95);
        x97=x60*(x91 - x95);
        x126=pow(x25, -3.0/2.0);
        x129=x21*x56;
        x198=0.5*ut1*ut2*x1*x56*x77;
        x199=-x92 - x94;
        x200=0.5*x126*x199*x21 + x198;

        /* Assignment result[1, 2]=-x129*x96 + x129*x97 - x200*x50 + x200*x59 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x81=ut1*ut2*x1*x77;
        x83=2*x20;
        x90=ut2*x78;
        x91=ut2*x1 + x7*x90;
        x92=x21*x81;
        x93=x1*x4*x77;
        x94=(1.0/2.0)*x23*(x83 + 2*x93);
        x95=x56*(x92 + x94);
        x96=x51*(x91 + x95);
        x97=x60*(x91 - x95);
        x126=pow(x25, -3.0/2.0);
        x129=x21*x56;
        x198=0.5*ut1*ut2*x1*x56*x77;
        x199=-x92 - x94;
        x200=0.5*x126*x199*x21 + x198;result[7] = -x129*x96 + x129*x97 - x200*x50 + x200*x59;
    }
    else if (x135)
    {
        DEBUG_PRINT("Case (x135) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x81=ut1*ut2*x1*x77;
        x83=2*x20;
        x90=ut2*x78;
        x91=ut2*x1 + x7*x90;
        x92=x21*x81;
        x93=x1*x4*x77;
        x94=(1.0/2.0)*x23*(x83 + 2*x93);
        x95=x56*(x92 + x94);
        x96=x51*(x91 + x95);
        x97=x60*(x91 - x95);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x133=random1*x132;

        /* Assignment result[1, 2]=-x133*x96 + x133*x97 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x81=ut1*ut2*x1*x77;
        x83=2*x20;
        x90=ut2*x78;
        x91=ut2*x1 + x7*x90;
        x92=x21*x81;
        x93=x1*x4*x77;
        x94=(1.0/2.0)*x23*(x83 + 2*x93);
        x95=x56*(x92 + x94);
        x96=x51*(x91 + x95);
        x97=x60*(x91 - x95);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x133=random1*x132;result[7] = -x133*x96 + x133*x97;
    }

    /* Assignment result[2, 2]=Piecewise((mu, x10), (x148, x36), (x152*(-x146*x172 + x146*x173 + x158 - x171*x175 + x171*x176 - 4.0*x180 - 4.0*x181 - 2.0*x182 - 2.0*x183 + x221 + x222 + x224 + x225 + x227 + x228 - x230*x41 + x230*x43 - x231*x41 + x231*x43 - x232*x41 + x232*x43), x49), (mu - x218*x96 + x218*x97 - x233*x50 + x233*x59, x131), (mu - x219*x96 + x219*x97, x135)) */
    double x230;
    double x231;
    double x232;
    double x233;if (x10)
    {
        DEBUG_PRINT("Case (x10) is True.\n");

        /* Assignment result[2, 2]=mu */
        result[8] = mu;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x121=pow(mu, 3);
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x138=4.2426406871192847703*x11;
        x139=4.2426406871192847703*x12;
        x140=mu*x12;
        x141=2.0*x140;
        x142=x1*x11;
        x143=x11*x121;
        x144=x1*x139;
        x145=-x144;
        x146=2.0*x121;
        x147=x12*x146;
        x148=x137*(-26.0*x100 + 16.0*x101 + 11.313708498984761164*x12*x142 - x138 + x139 + x141 - 57.982756057296896302*x142 - 34.0*x143 + x145 + x147);

        /* Assignment result[2, 2]=x148 */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x121=pow(mu, 3);
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x138=4.2426406871192847703*x11;
        x139=4.2426406871192847703*x12;
        x140=mu*x12;
        x141=2.0*x140;
        x142=x1*x11;
        x143=x11*x121;
        x144=x1*x139;
        x145=-x144;
        x146=2.0*x121;
        x147=x12*x146;
        x148=x137*(-26.0*x100 + 16.0*x101 + 11.313708498984761164*x12*x142 - x138 + x139 + x141 - 57.982756057296896302*x142 - 34.0*x143 + x145 + x147);result[8] = x148;
    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x71=x27*x30 + x28*x30;
        x119=pow(x71, 3.0/2.0);
        x121=pow(mu, 3);
        x122=pow(rn, 3);
        x146=2.0*x121;
        x149=4.0*x119*x41*x43;
        x150=x149*x27;
        x151=x149*x28;
        x152=1.0/(assert(IS_NOT_ZERO((x150 + x151))), (x150 + x151));
        x153=pow(mu, 5);
        x154=pow(rn, 4);
        x155=2.0*un*x153*x154*x27*x28;
        x156=pow(un, 3);
        x157=2.0*x121*x156*x27*x28*x29;
        x158=mu*x150 + mu*x151 - x155*x41 + x155*x43 - x157*x41 + x157*x43;
        x160=pow(mu, 4);
        x165=pow(rt2, 4);
        x171=4.0*x121;
        x172=un*x165*x27*x29*x41;
        x173=un*x165*x27*x29*x43;
        x174=pow(rt1, 4);
        x175=un*x174*x28*x29*x41;
        x176=un*x174*x28*x29*x43;
        x180=mu*un*x119*x27*x41;
        x181=mu*un*x119*x27*x43;
        x182=mu*un*x119*x28*x41;
        x183=mu*un*x119*x28*x43;
        x187=pow(rt2, 5);
        x190=pow(rt2, 3);
        x220=1.4142135623730951455*un*x122*x160*x187;
        x221=-x220*x41;
        x222=x220*x43;
        x223=1.4142135623730951455*rt2*un*x122*x160*x174;
        x224=-x223*x41;
        x225=x223*x43;
        x226=2.8284271247461902909*un*x122*x160*x190*x27;
        x227=-x226*x41;
        x228=x226*x43;
        x230=2.0*pow(rt1, 6)*un*x121*x29;
        x231=2.0*un*x153*x154*x174;
        x232=2.0*x121*x156*x174*x29;

        /* Assignment result[2, 2]=x152*(-x146*x172 + x146*x173 + x158 - x171*x175 + x171*x176 - 4.0*x180 - 4.0*x181 - 2.0*x182 - 2.0*x183 + x221 + x222 + x224 + x225 + x227 + x228 - x230*x41 + x230*x43 - x231*x41 + x231*x43 - x232*x41 + x232*x43) */
        x37=mu*rn;
        x39=2.0*x37*(assert(IS_POSITIVE(x38)), sqrt(x38));
        x40=un*un + x27 + x28 + x30;
        x41=(assert(IS_POSITIVE(-x39 + x40)), sqrt(-x39 + x40));
        x43=(assert(IS_POSITIVE(x39 + x40)), sqrt(x39 + x40));
        x71=x27*x30 + x28*x30;
        x119=pow(x71, 3.0/2.0);
        x121=pow(mu, 3);
        x122=pow(rn, 3);
        x146=2.0*x121;
        x149=4.0*x119*x41*x43;
        x150=x149*x27;
        x151=x149*x28;
        x152=1.0/(assert(IS_NOT_ZERO((x150 + x151))), (x150 + x151));
        x153=pow(mu, 5);
        x154=pow(rn, 4);
        x155=2.0*un*x153*x154*x27*x28;
        x156=pow(un, 3);
        x157=2.0*x121*x156*x27*x28*x29;
        x158=mu*x150 + mu*x151 - x155*x41 + x155*x43 - x157*x41 + x157*x43;
        x160=pow(mu, 4);
        x165=pow(rt2, 4);
        x171=4.0*x121;
        x172=un*x165*x27*x29*x41;
        x173=un*x165*x27*x29*x43;
        x174=pow(rt1, 4);
        x175=un*x174*x28*x29*x41;
        x176=un*x174*x28*x29*x43;
        x180=mu*un*x119*x27*x41;
        x181=mu*un*x119*x27*x43;
        x182=mu*un*x119*x28*x41;
        x183=mu*un*x119*x28*x43;
        x187=pow(rt2, 5);
        x190=pow(rt2, 3);
        x220=1.4142135623730951455*un*x122*x160*x187;
        x221=-x220*x41;
        x222=x220*x43;
        x223=1.4142135623730951455*rt2*un*x122*x160*x174;
        x224=-x223*x41;
        x225=x223*x43;
        x226=2.8284271247461902909*un*x122*x160*x190*x27;
        x227=-x226*x41;
        x228=x226*x43;
        x230=2.0*pow(rt1, 6)*un*x121*x29;
        x231=2.0*un*x153*x154*x174;
        x232=2.0*x121*x156*x174*x29;result[8] = x152*(-x146*x172 + x146*x173 + x158 - x171*x175 + x171*x176 - 4.0*x180 - 4.0*x181 - 2.0*x182 - 2.0*x183 + x221 + x222 + x224 + x225 + x227 + x228 - x230*x41 + x230*x43 - x231*x41 + x231*x43 - x232*x41 + x232*x43);
    }
    else if (x131)
    {
        DEBUG_PRINT("Case (x131) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x81=ut1*ut2*x1*x77;
        x83=2*x20;
        x90=ut2*x78;
        x91=ut2*x1 + x7*x90;
        x92=x21*x81;
        x93=x1*x4*x77;
        x94=(1.0/2.0)*x23*(x83 + 2*x93);
        x95=x56*(x92 + x94);
        x96=x51*(x91 + x95);
        x97=x60*(x91 - x95);
        x126=pow(x25, -3.0/2.0);
        x199=-x92 - x94;
        x218=x23*x56;
        x233=0.5*x126*x199*x23 + 0.5*x56*(x20 + x93);

        /* Assignment result[2, 2]=mu - x218*x96 + x218*x97 - x233*x50 + x233*x59 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x81=ut1*ut2*x1*x77;
        x83=2*x20;
        x90=ut2*x78;
        x91=ut2*x1 + x7*x90;
        x92=x21*x81;
        x93=x1*x4*x77;
        x94=(1.0/2.0)*x23*(x83 + 2*x93);
        x95=x56*(x92 + x94);
        x96=x51*(x91 + x95);
        x97=x60*(x91 - x95);
        x126=pow(x25, -3.0/2.0);
        x199=-x92 - x94;
        x218=x23*x56;
        x233=0.5*x126*x199*x23 + 0.5*x56*(x20 + x93);result[8] = mu - x218*x96 + x218*x97 - x233*x50 + x233*x59;
    }
    else if (x135)
    {
        DEBUG_PRINT("Case (x135) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x81=ut1*ut2*x1*x77;
        x83=2*x20;
        x90=ut2*x78;
        x91=ut2*x1 + x7*x90;
        x92=x21*x81;
        x93=x1*x4*x77;
        x94=(1.0/2.0)*x23*(x83 + 2*x93);
        x95=x56*(x92 + x94);
        x96=x51*(x91 + x95);
        x97=x60*(x91 - x95);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x219=random2*x132;

        /* Assignment result[2, 2]=mu - x219*x96 + x219*x97 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x77=1.0/(assert(IS_NOT_ZERO(x6)), x6);
        x78=mu*x77;
        x81=ut1*ut2*x1*x77;
        x83=2*x20;
        x90=ut2*x78;
        x91=ut2*x1 + x7*x90;
        x92=x21*x81;
        x93=x1*x4*x77;
        x94=(1.0/2.0)*x23*(x83 + 2*x93);
        x95=x56*(x92 + x94);
        x96=x51*(x91 + x95);
        x97=x60*(x91 - x95);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x219=random2*x132;result[8] = mu - x219*x96 + x219*x97;
    }

    /* Assignment result[0, 3]=Piecewise((mu, x99), (x14*(-x1*x15 - x1*x16 + x101 - x11*x68 + x69), x36), (mu - x104 - x105, x107)) */
    if (x99)
    {
        DEBUG_PRINT("Case (x99) is True.\n");

        /* Assignment result[0, 3]=mu */
        result[9] = mu;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x15=0.5*x11;
        x16=0.5*x12;
        x68=0.70710678118654757274*mu;
        x69=x12*x68;
        x100=mu*x11;
        x101=x100*x12;

        /* Assignment result[0, 3]=x14*(-x1*x15 - x1*x16 + x101 - x11*x68 + x69) */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x15=0.5*x11;
        x16=0.5*x12;
        x68=0.70710678118654757274*mu;
        x69=x12*x68;
        x100=mu*x11;
        x101=x100*x12;result[9] = x14*(-x1*x15 - x1*x16 + x101 - x11*x68 + x69);
    }
    else if (x107)
    {
        DEBUG_PRINT("Case (x107) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x102=rn*x1;
        x103=x56*(mu*rt1*x21 + mu*rt2*x23);
        x104=x51*(x102 + x103);
        x105=x60*(x102 - x103);

        /* Assignment result[0, 3]=mu - x104 - x105 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x102=rn*x1;
        x103=x56*(mu*rt1*x21 + mu*rt2*x23);
        x104=x51*(x102 + x103);
        x105=x60*(x102 - x103);result[9] = mu - x104 - x105;
    }

    /* Assignment result[1, 3]=Piecewise((0.0, x99), (x201, x36), (-x104*x129 + x105*x129 - x203*x50 + x203*x59, x204), (-x104*x133 + x105*x133, x205)) */
    double x201;if (x99)
    {
        DEBUG_PRINT("Case (x99) is True.\n");

        /* Assignment result[1, 3]=0.0 */
        result[10] = 0.0;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x15=0.5*x11;
        x16=0.5*x12;
        x108=0.35355339059327378637*x11;
        x109=0.35355339059327378637*x12;
        x201=x14*(-mu*x15 - mu*x16 - x1*x108 + x1*x109);

        /* Assignment result[1, 3]=x201 */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x15=0.5*x11;
        x16=0.5*x12;
        x108=0.35355339059327378637*x11;
        x109=0.35355339059327378637*x12;
        x201=x14*(-mu*x15 - mu*x16 - x1*x108 + x1*x109);result[10] = x201;
    }
    else if (x204)
    {
        DEBUG_PRINT("Case (x204) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x52=mu*x21;
        x54=mu*x23;
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x102=rn*x1;
        x103=x56*(mu*rt1*x21 + mu*rt2*x23);
        x104=x51*(x102 + x103);
        x105=x60*(x102 - x103);
        x126=pow(x25, -3.0/2.0);
        x129=x21*x56;
        x202=-rt1*x52 - rt2*x54;
        x203=0.5*mu*rt1*x56 + 0.5*x126*x202*x21;

        /* Assignment result[1, 3]=-x104*x129 + x105*x129 - x203*x50 + x203*x59 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x52=mu*x21;
        x54=mu*x23;
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x102=rn*x1;
        x103=x56*(mu*rt1*x21 + mu*rt2*x23);
        x104=x51*(x102 + x103);
        x105=x60*(x102 - x103);
        x126=pow(x25, -3.0/2.0);
        x129=x21*x56;
        x202=-rt1*x52 - rt2*x54;
        x203=0.5*mu*rt1*x56 + 0.5*x126*x202*x21;result[10] = -x104*x129 + x105*x129 - x203*x50 + x203*x59;
    }
    else if (x205)
    {
        DEBUG_PRINT("Case (x205) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x102=rn*x1;
        x103=x56*(mu*rt1*x21 + mu*rt2*x23);
        x104=x51*(x102 + x103);
        x105=x60*(x102 - x103);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x133=random1*x132;

        /* Assignment result[1, 3]=-x104*x133 + x105*x133 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x102=rn*x1;
        x103=x56*(mu*rt1*x21 + mu*rt2*x23);
        x104=x51*(x102 + x103);
        x105=x60*(x102 - x103);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x133=random1*x132;result[10] = -x104*x133 + x105*x133;
    }

    /* Assignment result[2, 3]=Piecewise((0.0, x99), (x201, x36), (-x104*x218 + x105*x218 - x234*x50 + x234*x59, x204), (-x104*x219 + x105*x219, x205)) */
    double x234;if (x99)
    {
        DEBUG_PRINT("Case (x99) is True.\n");

        /* Assignment result[2, 3]=0.0 */
        result[11] = 0.0;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x15=0.5*x11;
        x16=0.5*x12;
        x108=0.35355339059327378637*x11;
        x109=0.35355339059327378637*x12;
        x201=x14*(-mu*x15 - mu*x16 - x1*x108 + x1*x109);

        /* Assignment result[2, 3]=x201 */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x15=0.5*x11;
        x16=0.5*x12;
        x108=0.35355339059327378637*x11;
        x109=0.35355339059327378637*x12;
        x201=x14*(-mu*x15 - mu*x16 - x1*x108 + x1*x109);result[11] = x201;
    }
    else if (x204)
    {
        DEBUG_PRINT("Case (x204) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x52=mu*x21;
        x54=mu*x23;
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x102=rn*x1;
        x103=x56*(mu*rt1*x21 + mu*rt2*x23);
        x104=x51*(x102 + x103);
        x105=x60*(x102 - x103);
        x126=pow(x25, -3.0/2.0);
        x202=-rt1*x52 - rt2*x54;
        x218=x23*x56;
        x234=0.5*mu*rt2*x56 + 0.5*x126*x202*x23;

        /* Assignment result[2, 3]=-x104*x218 + x105*x218 - x234*x50 + x234*x59 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x52=mu*x21;
        x54=mu*x23;
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x102=rn*x1;
        x103=x56*(mu*rt1*x21 + mu*rt2*x23);
        x104=x51*(x102 + x103);
        x105=x60*(x102 - x103);
        x126=pow(x25, -3.0/2.0);
        x202=-rt1*x52 - rt2*x54;
        x218=x23*x56;
        x234=0.5*mu*rt2*x56 + 0.5*x126*x202*x23;result[11] = -x104*x218 + x105*x218 - x234*x50 + x234*x59;
    }
    else if (x205)
    {
        DEBUG_PRINT("Case (x205) is True.\n");
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x102=rn*x1;
        x103=x56*(mu*rt1*x21 + mu*rt2*x23);
        x104=x51*(x102 + x103);
        x105=x60*(x102 - x103);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x219=random2*x132;

        /* Assignment result[2, 3]=-x104*x219 + x105*x219 */
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x102=rn*x1;
        x103=x56*(mu*rt1*x21 + mu*rt2*x23);
        x104=x51*(x102 + x103);
        x105=x60*(x102 - x103);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x219=random2*x132;result[11] = -x104*x219 + x105*x219;
    }

    /* Assignment result[0, 4]=Piecewise((0.0, x99), (x110, x36), (-x113 - x114, x107)) */
    double x110;
    double x111;
    double x112;
    double x113;
    double x114;if (x99)
    {
        DEBUG_PRINT("Case (x99) is True.\n");

        /* Assignment result[0, 4]=0.0 */
        result[12] = 0.0;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x15=0.5*x11;
        x16=0.5*x12;
        x17=-x15 - x16;
        x108=0.35355339059327378637*x11;
        x109=0.35355339059327378637*x12;
        x110=x14*(-mu*x108 + mu*x109 + x17);

        /* Assignment result[0, 4]=x110 */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x15=0.5*x11;
        x16=0.5*x12;
        x17=-x15 - x16;
        x108=0.35355339059327378637*x11;
        x109=0.35355339059327378637*x12;
        x110=x14*(-mu*x108 + mu*x109 + x17);result[12] = x110;
    }
    else if (x107)
    {
        DEBUG_PRINT("Case (x107) is True.\n");
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x112=x111*x21;
        x113=x51*(rt1 + x112);
        x114=x60*(rt1 - x112);

        /* Assignment result[0, 4]=-x113 - x114 */
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x112=x111*x21;
        x113=x51*(rt1 + x112);
        x114=x60*(rt1 - x112);result[12] = -x113 - x114;
    }

    /* Assignment result[1, 4]=Piecewise((1.00000000000000, x99), (x209, x36), (-x113*x129 + x114*x129 - x211*x50 + x211*x59 + 1, x204), (-x113*x133 + x114*x133 + 1, x205)) */
    double x206;
    double x207;
    double x208;
    double x209;
    double x210;
    double x211;if (x99)
    {
        DEBUG_PRINT("Case (x99) is True.\n");

        /* Assignment result[1, 4]=1.00000000000000 */
        result[13] = 1.0;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x140=mu*x12;
        x142=x1*x11;
        x206=9.8994949366116653522*x11;
        x207=4.0*x140;
        x208=1.4142135623730951455*x12;
        x209=x137*(-x1*x208 - 20.0*x100 + 9.8994949366116653522*x12 + x136 - 15.556349186104045046*x142 - x206 - x207);

        /* Assignment result[1, 4]=x209 */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x140=mu*x12;
        x142=x1*x11;
        x206=9.8994949366116653522*x11;
        x207=4.0*x140;
        x208=1.4142135623730951455*x12;
        x209=x137*(-x1*x208 - 20.0*x100 + 9.8994949366116653522*x12 + x136 - 15.556349186104045046*x142 - x206 - x207);result[13] = x209;
    }
    else if (x204)
    {
        DEBUG_PRINT("Case (x204) is True.\n");
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x112=x111*x21;
        x113=x51*(rt1 + x112);
        x114=x60*(rt1 - x112);
        x126=pow(x25, -3.0/2.0);
        x129=x21*x56;
        x210=0.5*mu*rn*x56;
        x211=-0.5*mu*rn*x126*x22 + x210;

        /* Assignment result[1, 4]=-x113*x129 + x114*x129 - x211*x50 + x211*x59 + 1 */
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x112=x111*x21;
        x113=x51*(rt1 + x112);
        x114=x60*(rt1 - x112);
        x126=pow(x25, -3.0/2.0);
        x129=x21*x56;
        x210=0.5*mu*rn*x56;
        x211=-0.5*mu*rn*x126*x22 + x210;result[13] = -x113*x129 + x114*x129 - x211*x50 + x211*x59 + 1;
    }
    else if (x205)
    {
        DEBUG_PRINT("Case (x205) is True.\n");
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x112=x111*x21;
        x113=x51*(rt1 + x112);
        x114=x60*(rt1 - x112);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x133=random1*x132;

        /* Assignment result[1, 4]=-x113*x133 + x114*x133 + 1 */
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x112=x111*x21;
        x113=x51*(rt1 + x112);
        x114=x60*(rt1 - x112);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x133=random1*x132;result[13] = -x113*x133 + x114*x133 + 1;
    }

    /* Assignment result[2, 4]=Piecewise((0.0, x99), (x212, x36), (-x113*x218 + x114*x218 + x214, x204), (-x113*x219 + x114*x219, x205)) */
    double x212;
    double x213;
    double x214;if (x99)
    {
        DEBUG_PRINT("Case (x99) is True.\n");

        /* Assignment result[2, 4]=0.0 */
        result[14] = 0.0;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x18=1.4142135623730951455*x11;
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x139=4.2426406871192847703*x12;
        x140=mu*x12;
        x144=x1*x139;
        x145=-x144;
        x206=9.8994949366116653522*x11;
        x207=4.0*x140;
        x208=1.4142135623730951455*x12;
        x212=x137*(x1*x206 + 4.0*x100 + x145 - x18 + x207 + x208);

        /* Assignment result[2, 4]=x212 */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x18=1.4142135623730951455*x11;
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x139=4.2426406871192847703*x12;
        x140=mu*x12;
        x144=x1*x139;
        x145=-x144;
        x206=9.8994949366116653522*x11;
        x207=4.0*x140;
        x208=1.4142135623730951455*x12;
        x212=x137*(x1*x206 + 4.0*x100 + x145 - x18 + x207 + x208);result[14] = x212;
    }
    else if (x204)
    {
        DEBUG_PRINT("Case (x204) is True.\n");
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x112=x111*x21;
        x113=x51*(rt1 + x112);
        x114=x60*(rt1 - x112);
        x126=pow(x25, -3.0/2.0);
        x213=0.5*mu*rn*x126*x21*x23;
        x214=x213*x50 - x213*x59;
        x218=x23*x56;

        /* Assignment result[2, 4]=-x113*x218 + x114*x218 + x214 */
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x112=x111*x21;
        x113=x51*(rt1 + x112);
        x114=x60*(rt1 - x112);
        x126=pow(x25, -3.0/2.0);
        x213=0.5*mu*rn*x126*x21*x23;
        x214=x213*x50 - x213*x59;
        x218=x23*x56;result[14] = -x113*x218 + x114*x218 + x214;
    }
    else if (x205)
    {
        DEBUG_PRINT("Case (x205) is True.\n");
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x112=x111*x21;
        x113=x51*(rt1 + x112);
        x114=x60*(rt1 - x112);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x219=random2*x132;

        /* Assignment result[2, 4]=-x113*x219 + x114*x219 */
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x112=x111*x21;
        x113=x51*(rt1 + x112);
        x114=x60*(rt1 - x112);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x219=random2*x132;result[14] = -x113*x219 + x114*x219;
    }

    /* Assignment result[0, 5]=Piecewise((0.0, x99), (x110, x36), (-x116 - x117, x107)) */
    double x115;
    double x116;
    double x117;if (x99)
    {
        DEBUG_PRINT("Case (x99) is True.\n");

        /* Assignment result[0, 5]=0.0 */
        result[15] = 0.0;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x15=0.5*x11;
        x16=0.5*x12;
        x17=-x15 - x16;
        x108=0.35355339059327378637*x11;
        x109=0.35355339059327378637*x12;
        x110=x14*(-mu*x108 + mu*x109 + x17);

        /* Assignment result[0, 5]=x110 */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x13=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x14=(assert(IS_NOT_ZERO(x11)), x13/x11);
        x15=0.5*x11;
        x16=0.5*x12;
        x17=-x15 - x16;
        x108=0.35355339059327378637*x11;
        x109=0.35355339059327378637*x12;
        x110=x14*(-mu*x108 + mu*x109 + x17);result[15] = x110;
    }
    else if (x107)
    {
        DEBUG_PRINT("Case (x107) is True.\n");
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x115=x111*x23;
        x116=x51*(rt2 + x115);
        x117=x60*(rt2 - x115);

        /* Assignment result[0, 5]=-x116 - x117 */
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x115=x111*x23;
        x116=x51*(rt2 + x115);
        x117=x60*(rt2 - x115);result[15] = -x116 - x117;
    }

    /* Assignment result[1, 5]=Piecewise((0.0, x99), (x212, x36), (-x116*x129 + x117*x129 + x214, x204), (-x116*x133 + x117*x133, x205)) */
    if (x99)
    {
        DEBUG_PRINT("Case (x99) is True.\n");

        /* Assignment result[1, 5]=0.0 */
        result[16] = 0.0;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x18=1.4142135623730951455*x11;
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x139=4.2426406871192847703*x12;
        x140=mu*x12;
        x144=x1*x139;
        x145=-x144;
        x206=9.8994949366116653522*x11;
        x207=4.0*x140;
        x208=1.4142135623730951455*x12;
        x212=x137*(x1*x206 + 4.0*x100 + x145 - x18 + x207 + x208);

        /* Assignment result[1, 5]=x212 */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x18=1.4142135623730951455*x11;
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x139=4.2426406871192847703*x12;
        x140=mu*x12;
        x144=x1*x139;
        x145=-x144;
        x206=9.8994949366116653522*x11;
        x207=4.0*x140;
        x208=1.4142135623730951455*x12;
        x212=x137*(x1*x206 + 4.0*x100 + x145 - x18 + x207 + x208);result[16] = x212;
    }
    else if (x204)
    {
        DEBUG_PRINT("Case (x204) is True.\n");
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x115=x111*x23;
        x116=x51*(rt2 + x115);
        x117=x60*(rt2 - x115);
        x126=pow(x25, -3.0/2.0);
        x129=x21*x56;
        x213=0.5*mu*rn*x126*x21*x23;
        x214=x213*x50 - x213*x59;

        /* Assignment result[1, 5]=-x116*x129 + x117*x129 + x214 */
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x115=x111*x23;
        x116=x51*(rt2 + x115);
        x117=x60*(rt2 - x115);
        x126=pow(x25, -3.0/2.0);
        x129=x21*x56;
        x213=0.5*mu*rn*x126*x21*x23;
        x214=x213*x50 - x213*x59;result[16] = -x116*x129 + x117*x129 + x214;
    }
    else if (x205)
    {
        DEBUG_PRINT("Case (x205) is True.\n");
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x115=x111*x23;
        x116=x51*(rt2 + x115);
        x117=x60*(rt2 - x115);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x133=random1*x132;

        /* Assignment result[1, 5]=-x116*x133 + x117*x133 */
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x115=x111*x23;
        x116=x51*(rt2 + x115);
        x117=x60*(rt2 - x115);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x133=random1*x132;result[16] = -x116*x133 + x117*x133;
    }

    /* Assignment result[2, 5]=Piecewise((1.00000000000000, x99), (x209, x36), (-x116*x218 + x117*x218 - x235*x50 + x235*x59 + 1, x204), (-x116*x219 + x117*x219 + 1, x205)) */
    double x235;if (x99)
    {
        DEBUG_PRINT("Case (x99) is True.\n");

        /* Assignment result[2, 5]=1.00000000000000 */
        result[17] = 1.0;
    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x140=mu*x12;
        x142=x1*x11;
        x206=9.8994949366116653522*x11;
        x207=4.0*x140;
        x208=1.4142135623730951455*x12;
        x209=x137*(-x1*x208 - 20.0*x100 + 9.8994949366116653522*x12 + x136 - 15.556349186104045046*x142 - x206 - x207);

        /* Assignment result[2, 5]=x209 */
        x11=(assert(IS_POSITIVE(-2.8284271247461902909*mu + x1 + 3.0)), sqrt(-2.8284271247461902909*mu + x1 + 3.0));
        x12=(assert(IS_POSITIVE(8.4852813742385695406*mu + 9.0*x1 + 3.0)), sqrt(8.4852813742385695406*mu + 9.0*x1 + 3.0));
        x19=x11*x12;
        x100=mu*x11;
        x101=x100*x12;
        x136=11.313708498984761164*x101 + 16.0*x19;
        x137=1.0/(assert(IS_NOT_ZERO(x136)), x136);
        x140=mu*x12;
        x142=x1*x11;
        x206=9.8994949366116653522*x11;
        x207=4.0*x140;
        x208=1.4142135623730951455*x12;
        x209=x137*(-x1*x208 - 20.0*x100 + 9.8994949366116653522*x12 + x136 - 15.556349186104045046*x142 - x206 - x207);result[17] = x209;
    }
    else if (x204)
    {
        DEBUG_PRINT("Case (x204) is True.\n");
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x115=x111*x23;
        x116=x51*(rt2 + x115);
        x117=x60*(rt2 - x115);
        x126=pow(x25, -3.0/2.0);
        x210=0.5*mu*rn*x56;
        x218=x23*x56;
        x235=-0.5*mu*rn*x126*x24 + x210;

        /* Assignment result[2, 5]=-x116*x218 + x117*x218 - x235*x50 + x235*x59 + 1 */
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x115=x111*x23;
        x116=x51*(rt2 + x115);
        x117=x60*(rt2 - x115);
        x126=pow(x25, -3.0/2.0);
        x210=0.5*mu*rn*x56;
        x218=x23*x56;
        x235=-0.5*mu*rn*x126*x24 + x210;result[17] = -x116*x218 + x117*x218 - x235*x50 + x235*x59 + 1;
    }
    else if (x205)
    {
        DEBUG_PRINT("Case (x205) is True.\n");
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x115=x111*x23;
        x116=x51*(rt2 + x115);
        x117=x60*(rt2 - x115);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x219=random2*x132;

        /* Assignment result[2, 5]=-x116*x219 + x117*x219 + 1 */
        x37=mu*rn;
        x50=(assert(IS_POSITIVE(x33)), sqrt(x33));
        x51=(assert(IS_NOT_ZERO(x50)), 0.5/x50);
        x56=1.0/(assert(IS_NOT_ZERO(x26)), x26);
        x59=(assert(IS_POSITIVE(x34)), sqrt(x34));
        x60=(assert(IS_NOT_ZERO(x59)), 0.5/x59);
        x111=x37*x56;
        x115=x111*x23;
        x116=x51*(rt2 + x115);
        x117=x60*(rt2 - x115);
        x132=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x219=random2*x132;result[17] = -x116*x219 + x117*x219 + 1;
    }
}

void fc3d_FischerBurmeisterFMeritGenerated(
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

    /* Assignment result=0.5*(mu*rn - x10 - x11 + x3)**2 + 0.5*(mu*ut1 + rt1 + x10*x16 - x11*x16)**2 + 0.5*(mu*ut2 + rt2 + x10*x17 - x11*x17)**2 */
    double x1;
    double x2;
    double x3;
    double x4;
    double x5;
    double x6;
    double x7;
    double x8;
    double x9;
    double x10;
    double x11;
    double x12;
    int x13;
    double x14;
    int x15;
    double x16;
    double x17;x1=ut1*ut1;
    x2=ut2*ut2;
    x3=mu*(assert(IS_POSITIVE(x1 + x2)), sqrt(x1 + x2)) + un;
    x4=mu*mu;
    x5=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x3*x3;
    x6=mu*rn*rt1 + mu*ut1*x3;
    x7=mu*rn*rt2 + mu*ut2*x3;
    x8=(assert(IS_POSITIVE(x6*x6 + x7*x7)), sqrt(x6*x6 + x7*x7));
    x9=2*x8;
    x10=0.5*(assert(IS_POSITIVE(x5 - x9)), sqrt(x5 - x9));
    x11=0.5*(assert(IS_POSITIVE(x5 + x9)), sqrt(x5 + x9));
    x12=1.0/(assert(IS_NOT_ZERO(x8)), x8);
    x13=x8 > 0;
    x14=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
    x15=x8 <= 0;
    x16=((x13) ? (x12*x6): (random1*x14));
    x17=((x13) ? (x12*x7): (random2*x14));
    result[0] = 0.5*(mu*rn - x10 - x11 + x3)*(mu*rn - x10 - x11 + x3) + 0.5*(mu*ut1 + rt1 + x10*x16 - x11*x16)*(mu*ut1 + rt1 + x10*x16 - x11*x16) + 0.5*(mu*ut2 + rt2 + x10*x17 - x11*x17)*(mu*ut2 + rt2 + x10*x17 - x11*x17);
}
void fc3d_FischerBurmeisterGradFMeritGenerated(
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
    double x5;
    double x7;
    double x8;
    double x9;
    double x10;
    double x11;
    double x12;
    int x45;
    double x4;
    double x6;
    double x13;
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
    double x25;
    double x26;
    double x27;
    double x28;
    double x29;
    double x30;
    double x31;
    double x32;
    double x33;
    double x34;
    double x35;
    double x36;
    double x37;
    double x38;
    double x39;
    double x40;
    double x41;
    double x42;
    double x43;
    double x44;
    x1=ut1*ut1;
    x2=ut2*ut2;
    x3=(assert(IS_POSITIVE(x1 + x2)), sqrt(x1 + x2));
    x5=mu*x3 + un;
    x7=mu*rn*rt1 + mu*ut1*x5;
    x8=x7*x7;
    x9=mu*rn*rt2 + mu*ut2*x5;
    x10=x9*x9;
    x11=x10 + x8;
    x12=(assert(IS_POSITIVE(x11)), sqrt(x11));
    x45=x12 > 0;
    int x53;
    double x46;
    double x47;
    double x48;
    double x49;
    double x50;
    double x51;
    double x52;
    x53=x12 <= 0;
    if (x45)
    {
        x4=mu*mu;
        x6=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5;
        x13=2*x12;
        x14=(assert(IS_POSITIVE(-x13 + x6)), sqrt(-x13 + x6));
        x15=(assert(IS_POSITIVE(x13 + x6)), sqrt(x13 + x6));
        x16=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x14 - 0.25*x15;
        x17=(assert(IS_NOT_ZERO(x15)), 1.0/x15);
        x18=rn*x4;
        x19=mu*rt1;
        x20=x19*x7;
        x21=mu*rt2;
        x22=x21*x9;
        x23=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x24=x23*(x20 + x22);
        x25=x17*(x18 + x24);
        x26=1.0/(assert(IS_NOT_ZERO(x14)), x14);
        x27=1.0*x26;
        x28=x18 - x24;
        x29=x16*(2*mu - x25 - x27*x28);
        x30=0.5*rt1;
        x31=0.5*mu*ut1;
        x32=0.25*x14*x23*x7 - 0.25*x15*x23*x7 + x30 + x31;
        x33=pow(x11, -3.0/2.0);
        x34=x33*(-x20 - x22);
        x35=x19*x23 + x34*x7;
        x36=1.0*x15;
        x37=x23*x26*x7;
        x38=x23*x7;
        x39=0.5*rt2;
        x40=0.5*mu*ut2;
        x41=0.25*x14*x23*x9 - 0.25*x15*x23*x9 + x39 + x40;
        x42=x21*x23 + x34*x9;
        x43=x23*x26*x9;
        x44=x23*x9;
    }
    else if (x53)
    {
        x4=mu*mu;
        x6=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5;
        x13=2*x12;
        x14=(assert(IS_POSITIVE(-x13 + x6)), sqrt(-x13 + x6));
        x15=(assert(IS_POSITIVE(x13 + x6)), sqrt(x13 + x6));
        x16=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x14 - 0.25*x15;
        x17=(assert(IS_NOT_ZERO(x15)), 1.0/x15);
        x18=rn*x4;
        x19=mu*rt1;
        x20=x19*x7;
        x21=mu*rt2;
        x22=x21*x9;
        x23=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x24=x23*(x20 + x22);
        x25=x17*(x18 + x24);
        x26=1.0/(assert(IS_NOT_ZERO(x14)), x14);
        x27=1.0*x26;
        x28=x18 - x24;
        x29=x16*(2*mu - x25 - x27*x28);
        x30=0.5*rt1;
        x31=0.5*mu*ut1;
        x39=0.5*rt2;
        x40=0.5*mu*ut2;
        x46=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x47=0.25*random1*x14*x46 - 0.25*random1*x15*x46 + x30 + x31;
        x48=random1*x26*x46;
        x49=random1*x46;
        x50=0.25*random2*x14*x46 - 0.25*random2*x15*x46 + x39 + x40;
        x51=random2*x26*x46;
        x52=random2*x46;
    }
    /* Assignment result[0, 0]=Piecewise((x29 + x32*(x14*x35 - x25*x38 + x28*x37 - x35*x36) + x41*(x14*x42 - x25*x44 + x28*x43 - x36*x42), x45), (x29 + x47*(-x25*x49 + x28*x48) + x50*(-x25*x52 + x28*x51), x53)) */

    if (x45)
    {
        DEBUG_PRINT("Case (x45) is True.\n");
        x4=mu*mu;
        x6=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5;
        x13=2*x12;
        x14=(assert(IS_POSITIVE(-x13 + x6)), sqrt(-x13 + x6));
        x15=(assert(IS_POSITIVE(x13 + x6)), sqrt(x13 + x6));
        x16=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x14 - 0.25*x15;
        x17=(assert(IS_NOT_ZERO(x15)), 1.0/x15);
        x18=rn*x4;
        x19=mu*rt1;
        x20=x19*x7;
        x21=mu*rt2;
        x22=x21*x9;
        x23=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x24=x23*(x20 + x22);
        x25=x17*(x18 + x24);
        x26=1.0/(assert(IS_NOT_ZERO(x14)), x14);
        x27=1.0*x26;
        x28=x18 - x24;
        x29=x16*(2*mu - x25 - x27*x28);
        x30=0.5*rt1;
        x31=0.5*mu*ut1;
        x32=0.25*x14*x23*x7 - 0.25*x15*x23*x7 + x30 + x31;
        x33=pow(x11, -3.0/2.0);
        x34=x33*(-x20 - x22);
        x35=x19*x23 + x34*x7;
        x36=1.0*x15;
        x37=x23*x26*x7;
        x38=x23*x7;
        x39=0.5*rt2;
        x40=0.5*mu*ut2;
        x41=0.25*x14*x23*x9 - 0.25*x15*x23*x9 + x39 + x40;
        x42=x21*x23 + x34*x9;
        x43=x23*x26*x9;
        x44=x23*x9;

        /* Assignment result[0, 0]=x29 + x32*(x14*x35 - x25*x38 + x28*x37 - x35*x36) + x41*(x14*x42 - x25*x44 + x28*x43 - x36*x42) */
        x4=mu*mu;
        x6=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5;
        x13=2*x12;
        x14=(assert(IS_POSITIVE(-x13 + x6)), sqrt(-x13 + x6));
        x15=(assert(IS_POSITIVE(x13 + x6)), sqrt(x13 + x6));
        x16=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x14 - 0.25*x15;
        x17=(assert(IS_NOT_ZERO(x15)), 1.0/x15);
        x18=rn*x4;
        x19=mu*rt1;
        x20=x19*x7;
        x21=mu*rt2;
        x22=x21*x9;
        x23=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x24=x23*(x20 + x22);
        x25=x17*(x18 + x24);
        x26=1.0/(assert(IS_NOT_ZERO(x14)), x14);
        x27=1.0*x26;
        x28=x18 - x24;
        x29=x16*(2*mu - x25 - x27*x28);
        x30=0.5*rt1;
        x31=0.5*mu*ut1;
        x32=0.25*x14*x23*x7 - 0.25*x15*x23*x7 + x30 + x31;
        x33=pow(x11, -3.0/2.0);
        x34=x33*(-x20 - x22);
        x35=x19*x23 + x34*x7;
        x36=1.0*x15;
        x37=x23*x26*x7;
        x38=x23*x7;
        x39=0.5*rt2;
        x40=0.5*mu*ut2;
        x41=0.25*x14*x23*x9 - 0.25*x15*x23*x9 + x39 + x40;
        x42=x21*x23 + x34*x9;
        x43=x23*x26*x9;
        x44=x23*x9;
        result[0] = x29 + x32*(x14*x35 - x25*x38 + x28*x37 - x35*x36) + x41*(x14*x42 - x25*x44 + x28*x43 - x36*x42);
    }
    else if (x53)
    {
        DEBUG_PRINT("Case (x53) is True.\n");
        x4=mu*mu;
        x6=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5;
        x13=2*x12;
        x14=(assert(IS_POSITIVE(-x13 + x6)), sqrt(-x13 + x6));
        x15=(assert(IS_POSITIVE(x13 + x6)), sqrt(x13 + x6));
        x16=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x14 - 0.25*x15;
        x17=(assert(IS_NOT_ZERO(x15)), 1.0/x15);
        x18=rn*x4;
        x19=mu*rt1;
        x20=x19*x7;
        x21=mu*rt2;
        x22=x21*x9;
        x23=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x24=x23*(x20 + x22);
        x25=x17*(x18 + x24);
        x26=1.0/(assert(IS_NOT_ZERO(x14)), x14);
        x27=1.0*x26;
        x28=x18 - x24;
        x29=x16*(2*mu - x25 - x27*x28);
        x30=0.5*rt1;
        x31=0.5*mu*ut1;
        x39=0.5*rt2;
        x40=0.5*mu*ut2;
        x46=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x47=0.25*random1*x14*x46 - 0.25*random1*x15*x46 + x30 + x31;
        x48=random1*x26*x46;
        x49=random1*x46;
        x50=0.25*random2*x14*x46 - 0.25*random2*x15*x46 + x39 + x40;
        x51=random2*x26*x46;
        x52=random2*x46;

        /* Assignment result[0, 0]=x29 + x47*(-x25*x49 + x28*x48) + x50*(-x25*x52 + x28*x51) */
        x4=mu*mu;
        x6=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5;
        x13=2*x12;
        x14=(assert(IS_POSITIVE(-x13 + x6)), sqrt(-x13 + x6));
        x15=(assert(IS_POSITIVE(x13 + x6)), sqrt(x13 + x6));
        x16=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x14 - 0.25*x15;
        x17=(assert(IS_NOT_ZERO(x15)), 1.0/x15);
        x18=rn*x4;
        x19=mu*rt1;
        x20=x19*x7;
        x21=mu*rt2;
        x22=x21*x9;
        x23=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x24=x23*(x20 + x22);
        x25=x17*(x18 + x24);
        x26=1.0/(assert(IS_NOT_ZERO(x14)), x14);
        x27=1.0*x26;
        x28=x18 - x24;
        x29=x16*(2*mu - x25 - x27*x28);
        x30=0.5*rt1;
        x31=0.5*mu*ut1;
        x39=0.5*rt2;
        x40=0.5*mu*ut2;
        x46=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x47=0.25*random1*x14*x46 - 0.25*random1*x15*x46 + x30 + x31;
        x48=random1*x26*x46;
        x49=random1*x46;
        x50=0.25*random2*x14*x46 - 0.25*random2*x15*x46 + x39 + x40;
        x51=random2*x26*x46;
        x52=random2*x46;
        result[0] = x29 + x47*(-x25*x49 + x28*x48) + x50*(-x25*x52 + x28*x51);
    }

    /* Assignment result[0, 1]=Piecewise((x32*(x14*x62 - x36*x62 + x37*x57 - x38*x56 + 2) + x41*(x43*x57 - x44*x56 + x60) + x58, x45), (x47*(x48*x57 - x49*x56 + 2) + x50*(x51*x57 - x52*x56) + x58, x53)) */
    double x54;
    double x55;
    double x56;
    double x57;
    double x58;
    double x59;
    double x60;
    double x61;
    double x62;
    if (x45)
    {
        DEBUG_PRINT("Case (x45) is True.\n");
        x4=mu*mu;
        x6=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5;
        x13=2*x12;
        x14=(assert(IS_POSITIVE(-x13 + x6)), sqrt(-x13 + x6));
        x15=(assert(IS_POSITIVE(x13 + x6)), sqrt(x13 + x6));
        x16=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x14 - 0.25*x15;
        x17=(assert(IS_NOT_ZERO(x15)), 1.0/x15);
        x23=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x26=1.0/(assert(IS_NOT_ZERO(x14)), x14);
        x27=1.0*x26;
        x30=0.5*rt1;
        x31=0.5*mu*ut1;
        x32=0.25*x14*x23*x7 - 0.25*x15*x23*x7 + x30 + x31;
        x33=pow(x11, -3.0/2.0);
        x36=1.0*x15;
        x37=x23*x26*x7;
        x38=x23*x7;
        x39=0.5*rt2;
        x40=0.5*mu*ut2;
        x41=0.25*x14*x23*x9 - 0.25*x15*x23*x9 + x39 + x40;
        x43=x23*x26*x9;
        x44=x23*x9;
        x54=mu*rn*x23;
        x55=x54*x7;
        x56=x17*(rt1 + x55);
        x57=rt1 - x55;
        x58=x16*(-x27*x57 - x56);
        x59=mu*rn*x33*x7*x9;
        x60=-x14*x59 + x36*x59;
        x61=mu*rn*x33;
        x62=x54 - x61*x8;

        /* Assignment result[0, 1]=x32*(x14*x62 - x36*x62 + x37*x57 - x38*x56 + 2) + x41*(x43*x57 - x44*x56 + x60) + x58 */
        x4=mu*mu;
        x6=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5;
        x13=2*x12;
        x14=(assert(IS_POSITIVE(-x13 + x6)), sqrt(-x13 + x6));
        x15=(assert(IS_POSITIVE(x13 + x6)), sqrt(x13 + x6));
        x16=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x14 - 0.25*x15;
        x17=(assert(IS_NOT_ZERO(x15)), 1.0/x15);
        x23=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x26=1.0/(assert(IS_NOT_ZERO(x14)), x14);
        x27=1.0*x26;
        x30=0.5*rt1;
        x31=0.5*mu*ut1;
        x32=0.25*x14*x23*x7 - 0.25*x15*x23*x7 + x30 + x31;
        x33=pow(x11, -3.0/2.0);
        x36=1.0*x15;
        x37=x23*x26*x7;
        x38=x23*x7;
        x39=0.5*rt2;
        x40=0.5*mu*ut2;
        x41=0.25*x14*x23*x9 - 0.25*x15*x23*x9 + x39 + x40;
        x43=x23*x26*x9;
        x44=x23*x9;
        x54=mu*rn*x23;
        x55=x54*x7;
        x56=x17*(rt1 + x55);
        x57=rt1 - x55;
        x58=x16*(-x27*x57 - x56);
        x59=mu*rn*x33*x7*x9;
        x60=-x14*x59 + x36*x59;
        x61=mu*rn*x33;
        x62=x54 - x61*x8;
        result[1] = x32*(x14*x62 - x36*x62 + x37*x57 - x38*x56 + 2) + x41*(x43*x57 - x44*x56 + x60) + x58;
    }
    else if (x53)
    {
        DEBUG_PRINT("Case (x53) is True.\n");
        x4=mu*mu;
        x6=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5;
        x13=2*x12;
        x14=(assert(IS_POSITIVE(-x13 + x6)), sqrt(-x13 + x6));
        x15=(assert(IS_POSITIVE(x13 + x6)), sqrt(x13 + x6));
        x16=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x14 - 0.25*x15;
        x17=(assert(IS_NOT_ZERO(x15)), 1.0/x15);
        x23=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x26=1.0/(assert(IS_NOT_ZERO(x14)), x14);
        x27=1.0*x26;
        x30=0.5*rt1;
        x31=0.5*mu*ut1;
        x39=0.5*rt2;
        x40=0.5*mu*ut2;
        x46=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x47=0.25*random1*x14*x46 - 0.25*random1*x15*x46 + x30 + x31;
        x48=random1*x26*x46;
        x49=random1*x46;
        x50=0.25*random2*x14*x46 - 0.25*random2*x15*x46 + x39 + x40;
        x51=random2*x26*x46;
        x52=random2*x46;
        x54=mu*rn*x23;
        x55=x54*x7;
        x56=x17*(rt1 + x55);
        x57=rt1 - x55;
        x58=x16*(-x27*x57 - x56);

        /* Assignment result[0, 1]=x47*(x48*x57 - x49*x56 + 2) + x50*(x51*x57 - x52*x56) + x58 */
        x4=mu*mu;
        x6=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5;
        x13=2*x12;
        x14=(assert(IS_POSITIVE(-x13 + x6)), sqrt(-x13 + x6));
        x15=(assert(IS_POSITIVE(x13 + x6)), sqrt(x13 + x6));
        x16=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x14 - 0.25*x15;
        x17=(assert(IS_NOT_ZERO(x15)), 1.0/x15);
        x23=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x26=1.0/(assert(IS_NOT_ZERO(x14)), x14);
        x27=1.0*x26;
        x30=0.5*rt1;
        x31=0.5*mu*ut1;
        x39=0.5*rt2;
        x40=0.5*mu*ut2;
        x46=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x47=0.25*random1*x14*x46 - 0.25*random1*x15*x46 + x30 + x31;
        x48=random1*x26*x46;
        x49=random1*x46;
        x50=0.25*random2*x14*x46 - 0.25*random2*x15*x46 + x39 + x40;
        x51=random2*x26*x46;
        x52=random2*x46;
        x54=mu*rn*x23;
        x55=x54*x7;
        x56=x17*(rt1 + x55);
        x57=rt1 - x55;
        x58=x16*(-x27*x57 - x56);
        result[1] = x47*(x48*x57 - x49*x56 + 2) + x50*(x51*x57 - x52*x56) + x58;
    }

    /* Assignment result[0, 2]=Piecewise((x32*(x37*x65 - x38*x64 + x60) + x41*(x14*x67 - x36*x67 + x43*x65 - x44*x64 + 2) + x66, x45), (x47*(x48*x65 - x49*x64) + x50*(x51*x65 - x52*x64 + 2) + x66, x53)) */
    double x63;
    double x64;
    double x65;
    double x66;
    double x67;
    if (x45)
    {
        DEBUG_PRINT("Case (x45) is True.\n");
        x4=mu*mu;
        x6=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5;
        x13=2*x12;
        x14=(assert(IS_POSITIVE(-x13 + x6)), sqrt(-x13 + x6));
        x15=(assert(IS_POSITIVE(x13 + x6)), sqrt(x13 + x6));
        x16=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x14 - 0.25*x15;
        x17=(assert(IS_NOT_ZERO(x15)), 1.0/x15);
        x23=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x26=1.0/(assert(IS_NOT_ZERO(x14)), x14);
        x27=1.0*x26;
        x30=0.5*rt1;
        x31=0.5*mu*ut1;
        x32=0.25*x14*x23*x7 - 0.25*x15*x23*x7 + x30 + x31;
        x33=pow(x11, -3.0/2.0);
        x36=1.0*x15;
        x37=x23*x26*x7;
        x38=x23*x7;
        x39=0.5*rt2;
        x40=0.5*mu*ut2;
        x41=0.25*x14*x23*x9 - 0.25*x15*x23*x9 + x39 + x40;
        x43=x23*x26*x9;
        x44=x23*x9;
        x54=mu*rn*x23;
        x59=mu*rn*x33*x7*x9;
        x60=-x14*x59 + x36*x59;
        x61=mu*rn*x33;
        x63=x54*x9;
        x64=x17*(rt2 + x63);
        x65=rt2 - x63;
        x66=x16*(-x27*x65 - x64);
        x67=-x10*x61 + x54;

        /* Assignment result[0, 2]=x32*(x37*x65 - x38*x64 + x60) + x41*(x14*x67 - x36*x67 + x43*x65 - x44*x64 + 2) + x66 */
        x4=mu*mu;
        x6=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5;
        x13=2*x12;
        x14=(assert(IS_POSITIVE(-x13 + x6)), sqrt(-x13 + x6));
        x15=(assert(IS_POSITIVE(x13 + x6)), sqrt(x13 + x6));
        x16=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x14 - 0.25*x15;
        x17=(assert(IS_NOT_ZERO(x15)), 1.0/x15);
        x23=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x26=1.0/(assert(IS_NOT_ZERO(x14)), x14);
        x27=1.0*x26;
        x30=0.5*rt1;
        x31=0.5*mu*ut1;
        x32=0.25*x14*x23*x7 - 0.25*x15*x23*x7 + x30 + x31;
        x33=pow(x11, -3.0/2.0);
        x36=1.0*x15;
        x37=x23*x26*x7;
        x38=x23*x7;
        x39=0.5*rt2;
        x40=0.5*mu*ut2;
        x41=0.25*x14*x23*x9 - 0.25*x15*x23*x9 + x39 + x40;
        x43=x23*x26*x9;
        x44=x23*x9;
        x54=mu*rn*x23;
        x59=mu*rn*x33*x7*x9;
        x60=-x14*x59 + x36*x59;
        x61=mu*rn*x33;
        x63=x54*x9;
        x64=x17*(rt2 + x63);
        x65=rt2 - x63;
        x66=x16*(-x27*x65 - x64);
        x67=-x10*x61 + x54;
        result[2] = x32*(x37*x65 - x38*x64 + x60) + x41*(x14*x67 - x36*x67 + x43*x65 - x44*x64 + 2) + x66;
    }
    else if (x53)
    {
        DEBUG_PRINT("Case (x53) is True.\n");
        x4=mu*mu;
        x6=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5;
        x13=2*x12;
        x14=(assert(IS_POSITIVE(-x13 + x6)), sqrt(-x13 + x6));
        x15=(assert(IS_POSITIVE(x13 + x6)), sqrt(x13 + x6));
        x16=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x14 - 0.25*x15;
        x17=(assert(IS_NOT_ZERO(x15)), 1.0/x15);
        x23=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x26=1.0/(assert(IS_NOT_ZERO(x14)), x14);
        x27=1.0*x26;
        x30=0.5*rt1;
        x31=0.5*mu*ut1;
        x39=0.5*rt2;
        x40=0.5*mu*ut2;
        x46=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x47=0.25*random1*x14*x46 - 0.25*random1*x15*x46 + x30 + x31;
        x48=random1*x26*x46;
        x49=random1*x46;
        x50=0.25*random2*x14*x46 - 0.25*random2*x15*x46 + x39 + x40;
        x51=random2*x26*x46;
        x52=random2*x46;
        x54=mu*rn*x23;
        x63=x54*x9;
        x64=x17*(rt2 + x63);
        x65=rt2 - x63;
        x66=x16*(-x27*x65 - x64);

        /* Assignment result[0, 2]=x47*(x48*x65 - x49*x64) + x50*(x51*x65 - x52*x64 + 2) + x66 */
        x4=mu*mu;
        x6=rn*rn*x4 + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5;
        x13=2*x12;
        x14=(assert(IS_POSITIVE(-x13 + x6)), sqrt(-x13 + x6));
        x15=(assert(IS_POSITIVE(x13 + x6)), sqrt(x13 + x6));
        x16=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x14 - 0.25*x15;
        x17=(assert(IS_NOT_ZERO(x15)), 1.0/x15);
        x23=1.0/(assert(IS_NOT_ZERO(x12)), x12);
        x26=1.0/(assert(IS_NOT_ZERO(x14)), x14);
        x27=1.0*x26;
        x30=0.5*rt1;
        x31=0.5*mu*ut1;
        x39=0.5*rt2;
        x40=0.5*mu*ut2;
        x46=(assert(IS_POSITIVE(random1*random1 + random2*random2)), pow(random1*random1 + random2*random2, -1.0/2.0));
        x47=0.25*random1*x14*x46 - 0.25*random1*x15*x46 + x30 + x31;
        x48=random1*x26*x46;
        x49=random1*x46;
        x50=0.25*random2*x14*x46 - 0.25*random2*x15*x46 + x39 + x40;
        x51=random2*x26*x46;
        x52=random2*x46;
        x54=mu*rn*x23;
        x63=x54*x9;
        x64=x17*(rt2 + x63);
        x65=rt2 - x63;
        x66=x16*(-x27*x65 - x64);
        result[2] = x47*(x48*x65 - x49*x64) + x50*(x51*x65 - x52*x64 + 2) + x66;
    }
}



void fc3d_FischerBurmeisterFunctionGenerated(
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

    fc3d_FischerBurmeisterFABGenerated(
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
      fc3d_FischerBurmeisterFGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3(result, f);
    }

    if (A && B)
    {
      fc3d_FischerBurmeisterABGenerated(
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

void fc3d_FischerBurmeisterGradMeritFunctionGenerated(
  double *reaction,
  double *velocity,
  double mu,
  double *rho,
  double *gf)
{
  double result[3];

  assert(reaction);
  assert(velocity);
  assert(rho);

  SET3(reaction);
  SET3(velocity);
  SET3(rho);

  fc3d_FischerBurmeisterGradFMeritGenerated(
    *reaction0, *reaction1, *reaction2,
    *velocity0, *velocity1, *velocity2,
    mu,
    *rho0, *rho1, *rho2,
    result);
  cpy3(result, gf);
}
