#include <math.h>
#include <assert.h>
#include <op3x3.h>
#include <stdlib.h>

//#define DEBUG_MESSAGES 1
//#include <stdio.h>
#include <debug.h>
#include "NaturalMapGenerated.h"

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

void fc3d_NaturalMapFABGenerated(
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
    // Not supported in C:
    // Heaviside
    // Heaviside
    // Heaviside
    // Heaviside
    // Heaviside
    // Heaviside
    // Heaviside
    double x3;
    int x31;
    double x1;
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
    /*@ assert (ut1*ut1 + ut2*ut2) >= 0.;*/
    x3=sqrt(ut1*ut1 + ut2*ut2);
    /*@ assert (x3) >= 0.;*/
    x31=x3 <= ZERO;
    int x36;
    double x2;
    double x4;
    double x5;
    double x6;
    double x7;
    double x8;
    double x9;
    double x10;
    double x11;
    double x12;
    double x13;
    double x17;
    double x32;
    double x33;
    double x34;
    double x35;
    x36=x3 > ZERO;
    int x72;
    double x14;
    double x18;
    double x45;
    double x60;
    double x70;
    double x71;
    x5=mu*ut1;
    x6=-rt1 + x5;
    x7=x6*x6;
    x8=mu*ut2;
    x9=-rt2 + x8;
    x10=x9*x9;
    x11=x10 + x7;
    /*@ assert (x11) >= 0.;*/
    x12=sqrt(x11);
    x72=x12 > 0;
    int x73;
    x73=x12 <= 0;
    int x80;
    double x61;
    double x62;
    double x64;
    double x79;
    x80=x36 && x72;
    int x81;
    x81=x36 && x73;
    if (x31)
    {
        x1=mu*rn;
        x20=-1.0*un;
        /*@ assert (-1.0*un) >= 0.;*/
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x26=Heaviside(x20 + x25);
        x27=0.5*x26;
        x28=-1.0*x24;
        x29=Heaviside(x1 + x20 + x28);
        x30=0.5*x29;
    }
    else if (x36)
    {
        x1=mu*rn;
        x2=-un;
        /*@ assert (-un) >= 0.;*/
        x4=-mu*x3 + x1 + x2;
        x13=x12 + x4;
        x17=-x12 + x4;
        x32=Heaviside(x13);
        x33=0.5*x32;
        x34=Heaviside(x17);
        x35=0.5*x34;
    }
    else if (x72)
    {
        x1=mu*rn;
        x2=-un;
        /*@ assert (-un) >= 0.;*/
        x4=-mu*x3 + x1 + x2;
        x13=x12 + x4;
        x14=Max(0, x13);
        x17=-x12 + x4;
        x18=Max(0, x17);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x60=rt1 - x5;
        x70=0.5*x14*x45;
        x71=0.5*x18*x45;
    }
    else if (x80)
    {
        x1=mu*rn;
        x2=-un;
        /*@ assert (-un) >= 0.;*/
        x4=-mu*x3 + x1 + x2;
        x13=x12 + x4;
        x17=-x12 + x4;
        x32=Heaviside(x13);
        x34=Heaviside(x17);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x60=rt1 - x5;
        x61=0.5*x32*x45;
        x62=x60*x61;
        x64=0.5*x34*x45;
        x79=x6*x64;
    }
    /* Assignment result[0, 0]=x1 + x16 - x19 */
    double x15;
    double x16;
    double x19;x1=mu*rn;
    x2=-un;
    /*@ assert (-un) >= 0.;*/
    x4=-mu*x3 + x1 + x2;
    x13=x12 + x4;
    x14=Max(0, x13);
    x15=0.5*x14;
    x16=-x15;
    x17=-x12 + x4;
    x18=Max(0, x17);
    x19=0.5*x18;
    result[0] = x1 + x16 - x19;


    /* Assignment result[1, 0]=Piecewise((rt1 - x6*x71 - x60*x70, x72), (rt1, x73)) */

    if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x60=rt1 - x5;
        x70=0.5*x14*x45;
        x71=0.5*x18*x45;

        /* Assignment result[1, 0]=rt1 - x6*x71 - x60*x70 */
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x60=rt1 - x5;
        x70=0.5*x14*x45;
        x71=0.5*x18*x45;
        result[1] = rt1 - x6*x71 - x60*x70;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[1, 0]=rt1 */

        result[1] = rt1;

    }


    /* Assignment result[2, 0]=Piecewise((rt2 - x66*x70 - x71*x9, x72), (rt2 + x16 + x19, x73)) */
    double x66;
    if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x66=rt2 - x8;
        x70=0.5*x14*x45;
        x71=0.5*x18*x45;

        /* Assignment result[2, 0]=rt2 - x66*x70 - x71*x9 */
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x66=rt2 - x8;
        x70=0.5*x14*x45;
        x71=0.5*x18*x45;
        result[2] = rt2 - x66*x70 - x71*x9;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[2, 0]=rt2 + x16 + x19 */

        result[2] = rt2 + x16 + x19;

    }


    /* Assignment result[0, 1]=Piecewise((x27 + x30, x31), (x33 + x35, x36)) */

    if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x20=-1.0*un;
        /*@ assert (-1.0*un) >= 0.;*/
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x26=Heaviside(x20 + x25);
        x27=0.5*x26;
        x28=-1.0*x24;
        x29=Heaviside(x1 + x20 + x28);
        x30=0.5*x29;

        /* Assignment result[0, 1]=x27 + x30 */
        x20=-1.0*un;
        /*@ assert (-1.0*un) >= 0.;*/
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x26=Heaviside(x20 + x25);
        x27=0.5*x26;
        x28=-1.0*x24;
        x29=Heaviside(x1 + x20 + x28);
        x30=0.5*x29;
        result[3] = x27 + x30;

    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x32=Heaviside(x13);
        x33=0.5*x32;
        x34=Heaviside(x17);
        x35=0.5*x34;

        /* Assignment result[0, 1]=x33 + x35 */
        x32=Heaviside(x13);
        x33=0.5*x32;
        x34=Heaviside(x17);
        x35=0.5*x34;
        result[3] = x33 + x35;

    }


    /* Assignment result[1, 1]=Piecewise((-0.5*x37*x75*(x77 - 1.0*x78), x31), (x62 + x79, x80), (0, x81)) */
    double x37;
    double x74;
    double x75;
    double x76;
    double x77;
    double x78;
    if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x28=-1.0*x24;
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24;
        x74=Heaviside(x24);
        x75=rt1*x74;
        x76=un - 1.0*x1;
        x77=Heaviside(x28 + x76);
        x78=Heaviside(x24 + x76);

        /* Assignment result[1, 1]=-0.5*x37*x75*(x77 - 1.0*x78) */
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x28=-1.0*x24;
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24;
        x74=Heaviside(x24);
        x75=rt1*x74;
        x76=un - 1.0*x1;
        x77=Heaviside(x28 + x76);
        x78=Heaviside(x24 + x76);
        result[4] = -0.5*x37*x75*(x77 - 1.0*x78);

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");
        x32=Heaviside(x13);
        x34=Heaviside(x17);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x60=rt1 - x5;
        x61=0.5*x32*x45;
        x62=x60*x61;
        x64=0.5*x34*x45;
        x79=x6*x64;

        /* Assignment result[1, 1]=x62 + x79 */
        x32=Heaviside(x13);
        x34=Heaviside(x17);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x60=rt1 - x5;
        x61=0.5*x32*x45;
        x62=x60*x61;
        x64=0.5*x34*x45;
        x79=x6*x64;
        result[4] = x62 + x79;

    }
    else if (x81)
    {
        DEBUG_PRINT("Case (x81) is True.\n");

        /* Assignment result[1, 1]=0 */

        result[4] = 0;
        /*@ assert (result[4]) >= 0.;*/
    }


    /* Assignment result[2, 1]=Piecewise((x37*(x130*x77 - x130*x78 + x113*x129 - x115*x129 + x24*x27 - x24*x30), x31), (x131 + x67, x80), (x33 - x35, x81)) */
    double x67;
    double x113;
    double x115;
    double x129;
    double x130;
    double x131;
    if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x20=-1.0*un;
        /*@ assert (-1.0*un) >= 0.;*/
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x26=Heaviside(x20 + x25);
        x27=0.5*x26;
        x28=-1.0*x24;
        x29=Heaviside(x1 + x20 + x28);
        x30=0.5*x29;
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24;
        x74=Heaviside(x24);
        x76=un - 1.0*x1;
        x77=Heaviside(x28 + x76);
        x78=Heaviside(x24 + x76);
        x113=rt2*x78;
        x115=rt2*x77;
        x129=0.5*x74;
        x130=0.5*x24*x74;

        /* Assignment result[2, 1]=x37*(x130*x77 - x130*x78 + x113*x129 - x115*x129 + x24*x27 - x24*x30) */
        x20=-1.0*un;
        /*@ assert (-1.0*un) >= 0.;*/
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x26=Heaviside(x20 + x25);
        x27=0.5*x26;
        x28=-1.0*x24;
        x29=Heaviside(x1 + x20 + x28);
        x30=0.5*x29;
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24;
        x74=Heaviside(x24);
        x76=un - 1.0*x1;
        x77=Heaviside(x28 + x76);
        x78=Heaviside(x24 + x76);
        x113=rt2*x78;
        x115=rt2*x77;
        x129=0.5*x74;
        x130=0.5*x24*x74;
        result[5] = x37*(x130*x77 - x130*x78 + x113*x129 - x115*x129 + x24*x27 - x24*x30);

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");
        x32=Heaviside(x13);
        x34=Heaviside(x17);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x61=0.5*x32*x45;
        x64=0.5*x34*x45;
        x66=rt2 - x8;
        x67=x61*x66;
        x131=x64*x9;

        /* Assignment result[2, 1]=x131 + x67 */
        x32=Heaviside(x13);
        x34=Heaviside(x17);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x61=0.5*x32*x45;
        x64=0.5*x34*x45;
        x66=rt2 - x8;
        x67=x61*x66;
        x131=x64*x9;
        result[5] = x131 + x67;

    }
    else if (x81)
    {
        DEBUG_PRINT("Case (x81) is True.\n");
        x32=Heaviside(x13);
        x33=0.5*x32;
        x34=Heaviside(x17);
        x35=0.5*x34;

        /* Assignment result[2, 1]=x33 - x35 */
        x32=Heaviside(x13);
        x33=0.5*x32;
        x34=Heaviside(x17);
        x35=0.5*x34;
        result[5] = x33 - x35;

    }


    /* Assignment result[0, 2]=Piecewise((x38*(x39*rt1 - x40*rt1 + x42), x31), (x49 - x51, x36)) */
    double x38;
    double x39;
    double x40;
    double x41;
    double x42;
    double x43;
    double x44;
    double x46;
    double x47;
    double x48;
    double x49;
    double x50;
    double x51;
    if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x20=-1.0*un;
        /*@ assert (-1.0*un) >= 0.;*/
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x26=Heaviside(x20 + x25);
        x28=-1.0*x24;
        x29=Heaviside(x1 + x20 + x28);
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24;
        x38=0.25*x37*mu;
        x39=2.0*x26;
        x40=2.0*x29;
        x41=1.4142135623730951455*x24;
        x42=x26*x41 + x29*x41;

        /* Assignment result[0, 2]=x38*(x39*rt1 - x40*rt1 + x42) */
        x20=-1.0*un;
        /*@ assert (-1.0*un) >= 0.;*/
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x26=Heaviside(x20 + x25);
        x28=-1.0*x24;
        x29=Heaviside(x1 + x20 + x28);
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24;
        x38=0.25*x37*mu;
        x39=2.0*x26;
        x40=2.0*x29;
        x41=1.4142135623730951455*x24;
        x42=x26*x41 + x29*x41;
        result[6] = x38*(x39*rt1 - x40*rt1 + x42);

    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x32=Heaviside(x13);
        x33=0.5*x32;
        x34=Heaviside(x17);
        x35=0.5*x34;
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        x44=-x43*x5;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x47=x46*x6;
        x48=x44 + x47;
        x49=-x33*x48;
        x50=x44 - x47;
        x51=x35*x50;

        /* Assignment result[0, 2]=x49 - x51 */
        x32=Heaviside(x13);
        x33=0.5*x32;
        x34=Heaviside(x17);
        x35=0.5*x34;
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        x44=-x43*x5;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x47=x46*x6;
        x48=x44 + x47;
        x49=-x33*x48;
        x50=x44 - x47;
        x51=x35*x50;
        result[6] = x49 - x51;

    }


    /* Assignment result[1, 2]=Piecewise((-x74*x82*(x101*x77 + x101*x78 - x87*x91 + x87*x96 + x100 - x21*x84 - x85*x86 + x85*x89 - x86*x88 + x88*x89 - x93*x94 + x93*x99 - x98), x31), (-x19*(-x103*mu + x46) - x15*(-x105*mu + x104) - x48*x62 - x50*x79, x80), (0, x81)) */
    double x82;
    double x83;
    double x84;
    double x85;
    double x86;
    double x87;
    double x88;
    double x89;
    double x90;
    double x91;
    double x92;
    double x93;
    double x94;
    double x95;
    double x96;
    double x97;
    double x98;
    double x99;
    double x100;
    double x101;
    double x102;
    double x103;
    double x104;
    double x105;
    if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x28=-1.0*x24;
        x74=Heaviside(x24);
        x76=un - 1.0*x1;
        x77=Heaviside(x28 + x76);
        x78=Heaviside(x24 + x76);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=0.25*mu/pow(x23, 5.0/2.0);
        /*@ assert (x23) >= 0.;*/
        x83=pow(x23, 3.0/2.0);
        x84=4.0*x83;
        x85=pow(rt1, 5);
        x86=1.4142135623730951455*x78;
        x87=pow(rt2, 4);
        /*@ assert (x87) >= 0.;*/
        x88=rt1*x87;
        x89=1.4142135623730951455*x77;
        x90=Max(0, x2 + x25);
        x91=2.0*x90;
        x92=pow(rt1, 3);
        x93=x22*x92;
        x94=2.8284271247461902909*x78;
        x95=Max(0, x1 + x2 - x24);
        x96=2.0*x95;
        x97=x21*x22;
        x98=x91*x97;
        x99=2.8284271247461902909*x77;
        x100=x96*x97;
        x101=2.0*x21*x83;

        /* Assignment result[1, 2]=-x74*x82*(x101*x77 + x101*x78 - x87*x91 + x87*x96 + x100 - x21*x84 - x85*x86 + x85*x89 - x86*x88 + x88*x89 - x93*x94 + x93*x99 - x98) */
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x28=-1.0*x24;
        x74=Heaviside(x24);
        x76=un - 1.0*x1;
        x77=Heaviside(x28 + x76);
        x78=Heaviside(x24 + x76);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=0.25*mu/pow(x23, 5.0/2.0);
        /*@ assert (x23) >= 0.;*/
        x83=pow(x23, 3.0/2.0);
        x84=4.0*x83;
        x85=pow(rt1, 5);
        x86=1.4142135623730951455*x78;
        x87=pow(rt2, 4);
        /*@ assert (x87) >= 0.;*/
        x88=rt1*x87;
        x89=1.4142135623730951455*x77;
        x90=Max(0, x2 + x25);
        x91=2.0*x90;
        x92=pow(rt1, 3);
        x93=x22*x92;
        x94=2.8284271247461902909*x78;
        x95=Max(0, x1 + x2 - x24);
        x96=2.0*x95;
        x97=x21*x22;
        x98=x91*x97;
        x99=2.8284271247461902909*x77;
        x100=x96*x97;
        x101=2.0*x21*x83;
        result[7] = -x74*x82*(x101*x77 + x101*x78 - x87*x91 + x87*x96 + x100 - x21*x84 - x85*x86 + x85*x89 - x86*x88 + x88*x89 - x93*x94 + x93*x99 - x98);

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");
        x32=Heaviside(x13);
        x34=Heaviside(x17);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        x44=-x43*x5;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x47=x46*x6;
        x48=x44 + x47;
        x50=x44 - x47;
        x60=rt1 - x5;
        x61=0.5*x32*x45;
        x62=x60*x61;
        x64=0.5*x34*x45;
        x79=x6*x64;
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        x103=x102*x7;
        x104=-x46;
        x105=x102*x6*x60;

        /* Assignment result[1, 2]=-x19*(-x103*mu + x46) - x15*(-x105*mu + x104) - x48*x62 - x50*x79 */
        x32=Heaviside(x13);
        x34=Heaviside(x17);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        x44=-x43*x5;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x47=x46*x6;
        x48=x44 + x47;
        x50=x44 - x47;
        x60=rt1 - x5;
        x61=0.5*x32*x45;
        x62=x60*x61;
        x64=0.5*x34*x45;
        x79=x6*x64;
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        x103=x102*x7;
        x104=-x46;
        x105=x102*x6*x60;
        result[7] = -x19*(-x103*mu + x46) - x15*(-x105*mu + x104) - x48*x62 - x50*x79;

    }
    else if (x81)
    {
        DEBUG_PRINT("Case (x81) is True.\n");

        /* Assignment result[1, 2]=0 */

        result[7] = 0;
        /*@ assert (result[7]) >= 0.;*/
    }


    /* Assignment result[2, 2]=Piecewise((-x82*(-x135*x143 + x136*x143 - x149*x77 - x149*x78 - x106*x75 + x109*x75 - x111*x75 + x114*x75 + x116*x75 - x135*x144 + x136*x144 + x137 + x138*x85 + x138*x88 + x139*x93 - x140*x93 - x141*x93 - x142*x85 - x142*x88 - x145*x85 - x145*x88 - x146*x94 + x146*x99 + x147*x148 - x147*x150 + x151*x86 - x151*x89 - x39*x85 - x39*x88 - x40*x85 - x40*x88), x31), (x152*mu + x121 - x131*x50 - x48*x67, x80), (x49 + x51, x81)) */
    double x106;
    double x107;
    double x108;
    double x109;
    double x111;
    double x112;
    double x114;
    double x116;
    double x118;
    double x120;
    double x121;
    double x128;
    double x132;
    double x133;
    double x134;
    double x135;
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
    if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x20=-1.0*un;
        /*@ assert (-1.0*un) >= 0.;*/
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x26=Heaviside(x20 + x25);
        x28=-1.0*x24;
        x29=Heaviside(x1 + x20 + x28);
        x39=2.0*x26;
        x40=2.0*x29;
        x74=Heaviside(x24);
        x75=rt1*x74;
        x76=un - 1.0*x1;
        x77=Heaviside(x28 + x76);
        x78=Heaviside(x24 + x76);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=0.25*mu/pow(x23, 5.0/2.0);
        /*@ assert (x23) >= 0.;*/
        x83=pow(x23, 3.0/2.0);
        x84=4.0*x83;
        x85=pow(rt1, 5);
        x86=1.4142135623730951455*x78;
        x87=pow(rt2, 4);
        /*@ assert (x87) >= 0.;*/
        x88=rt1*x87;
        x89=1.4142135623730951455*x77;
        x90=Max(0, x2 + x25);
        x91=2.0*x90;
        x92=pow(rt1, 3);
        x93=x22*x92;
        x94=2.8284271247461902909*x78;
        x95=Max(0, x1 + x2 - x24);
        x96=2.0*x95;
        x99=2.8284271247461902909*x77;
        x106=x84*rt2;
        x107=pow(rt1, 4);
        /*@ assert (x107) >= 0.;*/
        x108=pow(rt2, 3);
        x109=x108*x91;
        x111=x108*x96;
        x112=2.0*x83;
        x113=rt2*x78;
        x114=x112*x113;
        x115=rt2*x77;
        x116=x112*x115;
        x132=1.4142135623730951455*x26*x83;
        x133=1.4142135623730951455*x29*x83;
        x134=x21*x83;
        x135=1.4142135623730951455*x74*x78;
        x136=1.4142135623730951455*x74*x77;
        x137=-x132*x21 - x132*x22 + x133*x21 + x133*x22 + x134*x135 - x134*x136;
        x138=4.0*x74;
        x139=8.0*x74;
        x140=4.0*x26;
        x141=4.0*x29;
        x142=2.0*x74*x78;
        x143=pow(rt2, 5);
        x144=rt2*x107;
        x145=2.0*x74*x77;
        x146=x108*x21*x74;
        x147=rt2*x92;
        x148=2.0*x74*x90;
        x149=4.0*x22*x74*x92;
        x150=2.0*x74*x95;
        x151=x22*x74*x83;

        /* Assignment result[2, 2]=-x82*(-x135*x143 + x136*x143 - x149*x77 - x149*x78 - x106*x75 + x109*x75 - x111*x75 + x114*x75 + x116*x75 - x135*x144 + x136*x144 + x137 + x138*x85 + x138*x88 + x139*x93 - x140*x93 - x141*x93 - x142*x85 - x142*x88 - x145*x85 - x145*x88 - x146*x94 + x146*x99 + x147*x148 - x147*x150 + x151*x86 - x151*x89 - x39*x85 - x39*x88 - x40*x85 - x40*x88) */
        x20=-1.0*un;
        /*@ assert (-1.0*un) >= 0.;*/
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x26=Heaviside(x20 + x25);
        x28=-1.0*x24;
        x29=Heaviside(x1 + x20 + x28);
        x39=2.0*x26;
        x40=2.0*x29;
        x74=Heaviside(x24);
        x75=rt1*x74;
        x76=un - 1.0*x1;
        x77=Heaviside(x28 + x76);
        x78=Heaviside(x24 + x76);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=0.25*mu/pow(x23, 5.0/2.0);
        /*@ assert (x23) >= 0.;*/
        x83=pow(x23, 3.0/2.0);
        x84=4.0*x83;
        x85=pow(rt1, 5);
        x86=1.4142135623730951455*x78;
        x87=pow(rt2, 4);
        /*@ assert (x87) >= 0.;*/
        x88=rt1*x87;
        x89=1.4142135623730951455*x77;
        x90=Max(0, x2 + x25);
        x91=2.0*x90;
        x92=pow(rt1, 3);
        x93=x22*x92;
        x94=2.8284271247461902909*x78;
        x95=Max(0, x1 + x2 - x24);
        x96=2.0*x95;
        x99=2.8284271247461902909*x77;
        x106=x84*rt2;
        x107=pow(rt1, 4);
        /*@ assert (x107) >= 0.;*/
        x108=pow(rt2, 3);
        x109=x108*x91;
        x111=x108*x96;
        x112=2.0*x83;
        x113=rt2*x78;
        x114=x112*x113;
        x115=rt2*x77;
        x116=x112*x115;
        x132=1.4142135623730951455*x26*x83;
        x133=1.4142135623730951455*x29*x83;
        x134=x21*x83;
        x135=1.4142135623730951455*x74*x78;
        x136=1.4142135623730951455*x74*x77;
        x137=-x132*x21 - x132*x22 + x133*x21 + x133*x22 + x134*x135 - x134*x136;
        x138=4.0*x74;
        x139=8.0*x74;
        x140=4.0*x26;
        x141=4.0*x29;
        x142=2.0*x74*x78;
        x143=pow(rt2, 5);
        x144=rt2*x107;
        x145=2.0*x74*x77;
        x146=x108*x21*x74;
        x147=rt2*x92;
        x148=2.0*x74*x90;
        x149=4.0*x22*x74*x92;
        x150=2.0*x74*x95;
        x151=x22*x74*x83;
        result[8] = -x82*(-x135*x143 + x136*x143 - x149*x77 - x149*x78 - x106*x75 + x109*x75 - x111*x75 + x114*x75 + x116*x75 - x135*x144 + x136*x144 + x137 + x138*x85 + x138*x88 + x139*x93 - x140*x93 - x141*x93 - x142*x85 - x142*x88 - x145*x85 - x145*x88 - x146*x94 + x146*x99 + x147*x148 - x147*x150 + x151*x86 - x151*x89 - x39*x85 - x39*x88 - x40*x85 - x40*x88);

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");
        x32=Heaviside(x13);
        x34=Heaviside(x17);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        x44=-x43*x5;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x47=x46*x6;
        x48=x44 + x47;
        x50=x44 - x47;
        x61=0.5*x32*x45;
        x64=0.5*x34*x45;
        x66=rt2 - x8;
        x67=x61*x66;
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        x118=0.5*x102*x14;
        x120=x102*x19*x6*x9;
        x121=x120*mu;
        x128=x6*x66;
        x131=x64*x9;
        x152=x118*x128;

        /* Assignment result[2, 2]=x152*mu + x121 - x131*x50 - x48*x67 */
        x32=Heaviside(x13);
        x34=Heaviside(x17);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        x44=-x43*x5;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x47=x46*x6;
        x48=x44 + x47;
        x50=x44 - x47;
        x61=0.5*x32*x45;
        x64=0.5*x34*x45;
        x66=rt2 - x8;
        x67=x61*x66;
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        x118=0.5*x102*x14;
        x120=x102*x19*x6*x9;
        x121=x120*mu;
        x128=x6*x66;
        x131=x64*x9;
        x152=x118*x128;
        result[8] = x152*mu + x121 - x131*x50 - x48*x67;

    }
    else if (x81)
    {
        DEBUG_PRINT("Case (x81) is True.\n");
        x32=Heaviside(x13);
        x33=0.5*x32;
        x34=Heaviside(x17);
        x35=0.5*x34;
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        x44=-x43*x5;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x47=x46*x6;
        x48=x44 + x47;
        x49=-x33*x48;
        x50=x44 - x47;
        x51=x35*x50;

        /* Assignment result[2, 2]=x49 + x51 */
        x32=Heaviside(x13);
        x33=0.5*x32;
        x34=Heaviside(x17);
        x35=0.5*x34;
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        x44=-x43*x5;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x47=x46*x6;
        x48=x44 + x47;
        x49=-x33*x48;
        x50=x44 - x47;
        x51=x35*x50;
        result[8] = x49 + x51;

    }


    /* Assignment result[0, 3]=Piecewise((x38*(x39*rt2 - x40*rt2 + x42), x31), (x55 - x57, x36)) */
    double x52;
    double x53;
    double x54;
    double x55;
    double x56;
    double x57;
    if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x20=-1.0*un;
        /*@ assert (-1.0*un) >= 0.;*/
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x26=Heaviside(x20 + x25);
        x28=-1.0*x24;
        x29=Heaviside(x1 + x20 + x28);
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24;
        x38=0.25*x37*mu;
        x39=2.0*x26;
        x40=2.0*x29;
        x41=1.4142135623730951455*x24;
        x42=x26*x41 + x29*x41;

        /* Assignment result[0, 3]=x38*(x39*rt2 - x40*rt2 + x42) */
        x20=-1.0*un;
        /*@ assert (-1.0*un) >= 0.;*/
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x26=Heaviside(x20 + x25);
        x28=-1.0*x24;
        x29=Heaviside(x1 + x20 + x28);
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24;
        x38=0.25*x37*mu;
        x39=2.0*x26;
        x40=2.0*x29;
        x41=1.4142135623730951455*x24;
        x42=x26*x41 + x29*x41;
        result[9] = x38*(x39*rt2 - x40*rt2 + x42);

    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x32=Heaviside(x13);
        x33=0.5*x32;
        x34=Heaviside(x17);
        x35=0.5*x34;
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x52=-x43*x8;
        x53=x46*x9;
        x54=x52 + x53;
        x55=-x33*x54;
        x56=x52 - x53;
        x57=x35*x56;

        /* Assignment result[0, 3]=x55 - x57 */
        x32=Heaviside(x13);
        x33=0.5*x32;
        x34=Heaviside(x17);
        x35=0.5*x34;
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x52=-x43*x8;
        x53=x46*x9;
        x54=x52 + x53;
        x55=-x33*x54;
        x56=x52 - x53;
        x57=x35*x56;
        result[9] = x55 - x57;

    }


    /* Assignment result[1, 3]=Piecewise((-x74*x82*rt1*(-x107*x86 + x107*x89 - x86*x87 + x87*x89 - x106 + x109 + x110*x91 - x110*x96 - x111 + x114 + x116 - x94*x97 + x97*x99), x31), (x119*mu + x121 - x54*x62 - x56*x79, x80), (0, x81)) */
    double x110;
    double x117;
    double x119;
    if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x28=-1.0*x24;
        x74=Heaviside(x24);
        x76=un - 1.0*x1;
        x77=Heaviside(x28 + x76);
        x78=Heaviside(x24 + x76);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=0.25*mu/pow(x23, 5.0/2.0);
        /*@ assert (x23) >= 0.;*/
        x83=pow(x23, 3.0/2.0);
        x84=4.0*x83;
        x86=1.4142135623730951455*x78;
        x87=pow(rt2, 4);
        /*@ assert (x87) >= 0.;*/
        x89=1.4142135623730951455*x77;
        x90=Max(0, x2 + x25);
        x91=2.0*x90;
        x94=2.8284271247461902909*x78;
        x95=Max(0, x1 + x2 - x24);
        x96=2.0*x95;
        x97=x21*x22;
        x99=2.8284271247461902909*x77;
        x106=x84*rt2;
        x107=pow(rt1, 4);
        /*@ assert (x107) >= 0.;*/
        x108=pow(rt2, 3);
        x109=x108*x91;
        x110=x21*rt2;
        x111=x108*x96;
        x112=2.0*x83;
        x113=rt2*x78;
        x114=x112*x113;
        x115=rt2*x77;
        x116=x112*x115;

        /* Assignment result[1, 3]=-x74*x82*rt1*(-x107*x86 + x107*x89 - x86*x87 + x87*x89 - x106 + x109 + x110*x91 - x110*x96 - x111 + x114 + x116 - x94*x97 + x97*x99) */
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x28=-1.0*x24;
        x74=Heaviside(x24);
        x76=un - 1.0*x1;
        x77=Heaviside(x28 + x76);
        x78=Heaviside(x24 + x76);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=0.25*mu/pow(x23, 5.0/2.0);
        /*@ assert (x23) >= 0.;*/
        x83=pow(x23, 3.0/2.0);
        x84=4.0*x83;
        x86=1.4142135623730951455*x78;
        x87=pow(rt2, 4);
        /*@ assert (x87) >= 0.;*/
        x89=1.4142135623730951455*x77;
        x90=Max(0, x2 + x25);
        x91=2.0*x90;
        x94=2.8284271247461902909*x78;
        x95=Max(0, x1 + x2 - x24);
        x96=2.0*x95;
        x97=x21*x22;
        x99=2.8284271247461902909*x77;
        x106=x84*rt2;
        x107=pow(rt1, 4);
        /*@ assert (x107) >= 0.;*/
        x108=pow(rt2, 3);
        x109=x108*x91;
        x110=x21*rt2;
        x111=x108*x96;
        x112=2.0*x83;
        x113=rt2*x78;
        x114=x112*x113;
        x115=rt2*x77;
        x116=x112*x115;
        result[10] = -x74*x82*rt1*(-x107*x86 + x107*x89 - x86*x87 + x87*x89 - x106 + x109 + x110*x91 - x110*x96 - x111 + x114 + x116 - x94*x97 + x97*x99);

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");
        x32=Heaviside(x13);
        x34=Heaviside(x17);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x52=-x43*x8;
        x53=x46*x9;
        x54=x52 + x53;
        x56=x52 - x53;
        x60=rt1 - x5;
        x61=0.5*x32*x45;
        x62=x60*x61;
        x64=0.5*x34*x45;
        x79=x6*x64;
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        x117=x60*x9;
        x118=0.5*x102*x14;
        x119=x117*x118;
        x120=x102*x19*x6*x9;
        x121=x120*mu;

        /* Assignment result[1, 3]=x119*mu + x121 - x54*x62 - x56*x79 */
        x32=Heaviside(x13);
        x34=Heaviside(x17);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x52=-x43*x8;
        x53=x46*x9;
        x54=x52 + x53;
        x56=x52 - x53;
        x60=rt1 - x5;
        x61=0.5*x32*x45;
        x62=x60*x61;
        x64=0.5*x34*x45;
        x79=x6*x64;
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        x117=x60*x9;
        x118=0.5*x102*x14;
        x119=x117*x118;
        x120=x102*x19*x6*x9;
        x121=x120*mu;
        result[10] = x119*mu + x121 - x54*x62 - x56*x79;

    }
    else if (x81)
    {
        DEBUG_PRINT("Case (x81) is True.\n");

        /* Assignment result[1, 3]=0 */

        result[10] = 0;
        /*@ assert (result[10]) >= 0.;*/
    }


    /* Assignment result[2, 3]=Piecewise((-x82*(-x107*x148 + x107*x150 - 1.17157287525381*x108*x21*x74*x77 - 6.82842712474619*x108*x21*x74*x78 + x138*x143 - x143*x154 - x143*x155 - x143*x39 - x143*x40 + x100*x74 - x22*x74*x84 - x74*x98 + 0.585786437626905*x151*x77 + 3.41421356237309*x151*x78 + x137 + x138*x144 + x139*x153 - x140*x153 - x141*x153 - x144*x154 - x144*x155 - x144*x39 - x144*x40), x31), (-x19*(-x156*mu + x46) - x15*(-x157*mu + x104) - x131*x56 - x54*x67, x80), (x55 + x57, x81)) */
    double x153;
    double x154;
    double x155;
    double x156;
    double x157;
    if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x20=-1.0*un;
        /*@ assert (-1.0*un) >= 0.;*/
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x26=Heaviside(x20 + x25);
        x28=-1.0*x24;
        x29=Heaviside(x1 + x20 + x28);
        x39=2.0*x26;
        x40=2.0*x29;
        x74=Heaviside(x24);
        x76=un - 1.0*x1;
        x77=Heaviside(x28 + x76);
        x78=Heaviside(x24 + x76);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=0.25*mu/pow(x23, 5.0/2.0);
        /*@ assert (x23) >= 0.;*/
        x83=pow(x23, 3.0/2.0);
        x84=4.0*x83;
        x90=Max(0, x2 + x25);
        x91=2.0*x90;
        x95=Max(0, x1 + x2 - x24);
        x96=2.0*x95;
        x97=x21*x22;
        x98=x91*x97;
        x100=x96*x97;
        x107=pow(rt1, 4);
        /*@ assert (x107) >= 0.;*/
        x108=pow(rt2, 3);
        x132=1.4142135623730951455*x26*x83;
        x133=1.4142135623730951455*x29*x83;
        x134=x21*x83;
        x135=1.4142135623730951455*x74*x78;
        x136=1.4142135623730951455*x74*x77;
        x137=-x132*x21 - x132*x22 + x133*x21 + x133*x22 + x134*x135 - x134*x136;
        x138=4.0*x74;
        x139=8.0*x74;
        x140=4.0*x26;
        x141=4.0*x29;
        x143=pow(rt2, 5);
        x144=rt2*x107;
        x148=2.0*x74*x90;
        x150=2.0*x74*x95;
        x151=x22*x74*x83;
        x153=x108*x21;
        x154=3.4142135623730949234*x74*x78;
        x155=0.58578643762690496555*x74*x77;

        /* Assignment result[2, 3]=-x82*(-x107*x148 + x107*x150 - 1.17157287525381*x108*x21*x74*x77 - 6.82842712474619*x108*x21*x74*x78 + x138*x143 - x143*x154 - x143*x155 - x143*x39 - x143*x40 + x100*x74 - x22*x74*x84 - x74*x98 + 0.585786437626905*x151*x77 + 3.41421356237309*x151*x78 + x137 + x138*x144 + x139*x153 - x140*x153 - x141*x153 - x144*x154 - x144*x155 - x144*x39 - x144*x40) */
        x20=-1.0*un;
        /*@ assert (-1.0*un) >= 0.;*/
        x21=rt1*rt1;
        x22=rt2*rt2;
        x23=x21 + x22;
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23);
        x25=x1 + x24;
        x26=Heaviside(x20 + x25);
        x28=-1.0*x24;
        x29=Heaviside(x1 + x20 + x28);
        x39=2.0*x26;
        x40=2.0*x29;
        x74=Heaviside(x24);
        x76=un - 1.0*x1;
        x77=Heaviside(x28 + x76);
        x78=Heaviside(x24 + x76);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=0.25*mu/pow(x23, 5.0/2.0);
        /*@ assert (x23) >= 0.;*/
        x83=pow(x23, 3.0/2.0);
        x84=4.0*x83;
        x90=Max(0, x2 + x25);
        x91=2.0*x90;
        x95=Max(0, x1 + x2 - x24);
        x96=2.0*x95;
        x97=x21*x22;
        x98=x91*x97;
        x100=x96*x97;
        x107=pow(rt1, 4);
        /*@ assert (x107) >= 0.;*/
        x108=pow(rt2, 3);
        x132=1.4142135623730951455*x26*x83;
        x133=1.4142135623730951455*x29*x83;
        x134=x21*x83;
        x135=1.4142135623730951455*x74*x78;
        x136=1.4142135623730951455*x74*x77;
        x137=-x132*x21 - x132*x22 + x133*x21 + x133*x22 + x134*x135 - x134*x136;
        x138=4.0*x74;
        x139=8.0*x74;
        x140=4.0*x26;
        x141=4.0*x29;
        x143=pow(rt2, 5);
        x144=rt2*x107;
        x148=2.0*x74*x90;
        x150=2.0*x74*x95;
        x151=x22*x74*x83;
        x153=x108*x21;
        x154=3.4142135623730949234*x74*x78;
        x155=0.58578643762690496555*x74*x77;
        result[11] = -x82*(-x107*x148 + x107*x150 - 1.1715728752538099311*x108*x21*x74*x77 - 6.8284271247461898469*x108*x21*x74*x78 + x138*x143 - x143*x154 - x143*x155 - x143*x39 - x143*x40 + x100*x74 - x22*x74*x84 - x74*x98 + 0.58578643762690496555*x151*x77 + 3.4142135623730949234*x151*x78 + x137 + x138*x144 + x139*x153 - x140*x153 - x141*x153 - x144*x154 - x144*x155 - x144*x39 - x144*x40);

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");
        x32=Heaviside(x13);
        x34=Heaviside(x17);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x52=-x43*x8;
        x53=x46*x9;
        x54=x52 + x53;
        x56=x52 - x53;
        x61=0.5*x32*x45;
        x64=0.5*x34*x45;
        x66=rt2 - x8;
        x67=x61*x66;
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        x104=-x46;
        x131=x64*x9;
        x156=x10*x102;
        x157=x102*x66*x9;

        /* Assignment result[2, 3]=-x19*(-x156*mu + x46) - x15*(-x157*mu + x104) - x131*x56 - x54*x67 */
        x32=Heaviside(x13);
        x34=Heaviside(x17);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x52=-x43*x8;
        x53=x46*x9;
        x54=x52 + x53;
        x56=x52 - x53;
        x61=0.5*x32*x45;
        x64=0.5*x34*x45;
        x66=rt2 - x8;
        x67=x61*x66;
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        x104=-x46;
        x131=x64*x9;
        x156=x10*x102;
        x157=x102*x66*x9;
        result[11] = -x19*(-x156*mu + x46) - x15*(-x157*mu + x104) - x131*x56 - x54*x67;

    }
    else if (x81)
    {
        DEBUG_PRINT("Case (x81) is True.\n");
        x32=Heaviside(x13);
        x33=0.5*x32;
        x34=Heaviside(x17);
        x35=0.5*x34;
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x52=-x43*x8;
        x53=x46*x9;
        x54=x52 + x53;
        x55=-x33*x54;
        x56=x52 - x53;
        x57=x35*x56;

        /* Assignment result[2, 3]=x55 + x57 */
        x32=Heaviside(x13);
        x33=0.5*x32;
        x34=Heaviside(x17);
        x35=0.5*x34;
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3;
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x52=-x43*x8;
        x53=x46*x9;
        x54=x52 + x53;
        x55=-x33*x54;
        x56=x52 - x53;
        x57=x35*x56;
        result[11] = x55 + x57;

    }


    /* Assignment result[0, 4]=mu + x58 - x59 */
    double x58;
    double x59;x32=Heaviside(x13);
    x33=0.5*x32;
    x34=Heaviside(x17);
    x35=0.5*x34;
    x58=-x33*mu;
    x59=x35*mu;
    result[12] = mu + x58 - x59;


    /* Assignment result[1, 4]=Piecewise((-x122*x60 - x35*x47, x72), (0, x73)) */
    double x122;
    if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x47=x46*x6;
        x60=rt1 - x5;
        x122=0.5*x32*x45*mu;

        /* Assignment result[1, 4]=-x122*x60 - x35*x47 */
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x47=x46*x6;
        x60=rt1 - x5;
        x122=0.5*x32*x45*mu;
        result[13] = -x122*x60 - x35*x47;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[1, 4]=0 */

        result[13] = 0;
        /*@ assert (result[13]) >= 0.;*/
    }


    /* Assignment result[2, 4]=Piecewise((-x122*x66 - x35*x53, x72), (x58 + x59, x73)) */

    if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x53=x46*x9;
        x66=rt2 - x8;
        x122=0.5*x32*x45*mu;

        /* Assignment result[2, 4]=-x122*x66 - x35*x53 */
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12;
        x46=x45*mu;
        x53=x46*x9;
        x66=rt2 - x8;
        x122=0.5*x32*x45*mu;
        result[14] = -x122*x66 - x35*x53;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[2, 4]=x58 + x59 */

        result[14] = x58 + x59;

    }


    /* Assignment result[0, 5]=x63 + x65 */
    double x63;
    double x65;/*@ assert (x12) != 0.;*/
    x45=1.0/x12;
    x60=rt1 - x5;
    x61=0.5*x32*x45;
    x62=x60*x61;
    x63=-x62;
    x64=0.5*x34*x45;
    x65=x60*x64;
    result[15] = x63 + x65;


    /* Assignment result[1, 5]=Piecewise((-x124*x60**2 - x19*(x103 + x126) - x15*(x105 + x45) + 1 + x125*x6*x60, x72), (1, x73)) */
    double x123;
    double x124;
    double x125;
    double x126;
    if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        x103=x102*x7;
        x105=x102*x6*x60;
        /*@ assert (x11) != 0.;*/
        x123=1.0/x11;
        x124=0.5*x123*x32;
        x125=0.5*x123*x34;
        x126=-x45;

        /* Assignment result[1, 5]=-x124*x60**2 - x19*(x103 + x126) - x15*(x105 + x45) + 1 + x125*x6*x60 */
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        x103=x102*x7;
        x105=x102*x6*x60;
        /*@ assert (x11) != 0.;*/
        x123=1.0/x11;
        x124=0.5*x123*x32;
        x125=0.5*x123*x34;
        x126=-x45;
        result[16] = -x124*x60*x60 - x19*(x103 + x126) - x15*(x105 + x45) + 1 + x125*x6*x60;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[1, 5]=1 */

        result[16] = 1;
        /*@ assert (result[16]) >= 0.;*/
        /*@ assert (result[16]) != 0.;*/
    }


    /* Assignment result[2, 5]=Piecewise((x117*x125 + x127 - x152, x72), (x63 - x65, x73)) */
    double x127;
    if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        x66=rt2 - x8;
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        x117=x60*x9;
        x118=0.5*x102*x14;
        x120=x102*x19*x6*x9;
        /*@ assert (x11) != 0.;*/
        x123=1.0/x11;
        x124=0.5*x123*x32;
        x125=0.5*x123*x34;
        x127=-x120 - x124*x60*x66;
        x128=x6*x66;
        x152=x118*x128;

        /* Assignment result[2, 5]=x117*x125 + x127 - x152 */
        x66=rt2 - x8;
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        x117=x60*x9;
        x118=0.5*x102*x14;
        x120=x102*x19*x6*x9;
        /*@ assert (x11) != 0.;*/
        x123=1.0/x11;
        x124=0.5*x123*x32;
        x125=0.5*x123*x34;
        x127=-x120 - x124*x60*x66;
        x128=x6*x66;
        x152=x118*x128;
        result[17] = x117*x125 + x127 - x152;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[2, 5]=x63 - x65 */

        result[17] = x63 - x65;

    }


    /* Assignment result[0, 6]=x68 + x69 */
    double x68;
    double x69;x66=rt2 - x8;
    x67=x61*x66;
    x68=-x67;
    x69=x64*x66;
    result[18] = x68 + x69;


    /* Assignment result[1, 6]=Piecewise((-x119 + x125*x128 + x127, x72), (0, x73)) */

    if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        x117=x60*x9;
        x118=0.5*x102*x14;
        x119=x117*x118;
        x120=x102*x19*x6*x9;
        /*@ assert (x11) != 0.;*/
        x123=1.0/x11;
        x124=0.5*x123*x32;
        x125=0.5*x123*x34;
        x127=-x120 - x124*x60*x66;
        x128=x6*x66;

        /* Assignment result[1, 6]=-x119 + x125*x128 + x127 */
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        x117=x60*x9;
        x118=0.5*x102*x14;
        x119=x117*x118;
        x120=x102*x19*x6*x9;
        /*@ assert (x11) != 0.;*/
        x123=1.0/x11;
        x124=0.5*x123*x32;
        x125=0.5*x123*x34;
        x127=-x120 - x124*x60*x66;
        x128=x6*x66;
        result[19] = -x119 + x125*x128 + x127;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[1, 6]=0 */

        result[19] = 0;
        /*@ assert (result[19]) >= 0.;*/
    }


    /* Assignment result[2, 6]=Piecewise((-x124*x66**2 - x19*(x126 + x156) - x15*(x157 + x45) + 1 + x125*x66*x9, x72), (1 + x68 - x69, x73)) */

    if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        /*@ assert (x11) != 0.;*/
        x123=1.0/x11;
        x124=0.5*x123*x32;
        x125=0.5*x123*x34;
        x126=-x45;
        x156=x10*x102;
        x157=x102*x66*x9;

        /* Assignment result[2, 6]=-x124*x66**2 - x19*(x126 + x156) - x15*(x157 + x45) + 1 + x125*x66*x9 */
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x102=pow(x11, -3.0/2.0);
        /*@ assert (x11) != 0.;*/
        x123=1.0/x11;
        x124=0.5*x123*x32;
        x125=0.5*x123*x34;
        x126=-x45;
        x156=x10*x102;
        x157=x102*x66*x9;
        result[20] = -x124*x66*x66 - x19*(x126 + x156) - x15*(x157 + x45) + 1 + x125*x66*x9;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[2, 6]=1 + x68 - x69 */

        result[20] = 1 + x68 - x69;

    }
}
void fc3d_NaturalMapFGenerated(
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
    double x3;
    double x4;
    double x5;
    double x6;
    double x7;
    int x15;
    double x1;
    double x2;
    double x8;
    double x10;
    double x12;
    double x13;
    double x14;
    x3=mu*ut1;
    x4=-rt1 + x3;
    x5=mu*ut2;
    x6=-rt2 + x5;
    /*@ assert (x6*x6 + x4*x4) >= 0.;*/
    x7=sqrt(x6*x6 + x4*x4);
    x15=x7 > 0;
    int x16;
    x16=x7 <= 0;
    if (x15)
    {
        x1=mu*rn;
        /*@ assert (ut1*ut1 + ut2*ut2) >= 0.;*/
        x2=-mu*sqrt(ut1*ut1 + ut2*ut2) - un + x1;
        x8=Max(0, x2 + x7);
        x10=Max(0, x2 - x7);
        /*@ assert (x7) != 0.;*/
        x12=1.0/x7;
        x13=0.5*x8*x12;
        x14=0.5*x12*x10;
    }
    /* Assignment result[0, 0]=x1 - x11 + x9 */
    double x9;
    double x11;x1=mu*rn;
    /*@ assert (ut1*ut1 + ut2*ut2) >= 0.;*/
    x2=-mu*sqrt(ut1*ut1 + ut2*ut2) - un + x1;
    x8=Max(0, x2 + x7);
    x9=-0.5*x8;
    x10=Max(0, x2 - x7);
    x11=0.5*x10;
    result[0] = x1 - x11 + x9;


    /* Assignment result[1, 0]=Piecewise((rt1 - (rt1 - x3)*x13 - x14*x4, x15), (rt1, x16)) */

    if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");
        /*@ assert (x7) != 0.;*/
        x12=1.0/x7;
        x13=0.5*x8*x12;
        x14=0.5*x12*x10;

        /* Assignment result[1, 0]=rt1 - (rt1 - x3)*x13 - x14*x4 */
        /*@ assert (x7) != 0.;*/
        x12=1.0/x7;
        x13=0.5*x8*x12;
        x14=0.5*x12*x10;
        result[1] = rt1 - (rt1 - x3)*x13 - x14*x4;

    }
    else if (x16)
    {
        DEBUG_PRINT("Case (x16) is True.\n");

        /* Assignment result[1, 0]=rt1 */

        result[1] = rt1;

    }
    /*@ assert (result[1]) >= 0.;*/

    /* Assignment result[2, 0]=Piecewise((rt2 - (rt2 - x5)*x13 - x14*x6, x15), (rt2 + x11 + x9, x16)) */

    if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");
        /*@ assert (x7) != 0.;*/
        x12=1.0/x7;
        x13=0.5*x8*x12;
        x14=0.5*x12*x10;

        /* Assignment result[2, 0]=rt2 - (rt2 - x5)*x13 - x14*x6 */
        /*@ assert (x7) != 0.;*/
        x12=1.0/x7;
        x13=0.5*x8*x12;
        x14=0.5*x12*x10;
        result[2] = rt2 - (rt2 - x5)*x13 - x14*x6;

    }
    else if (x16)
    {
        DEBUG_PRINT("Case (x16) is True.\n");

        /* Assignment result[2, 0]=rt2 + x11 + x9 */

        result[2] = rt2 + x11 + x9;

    }
    /*@ assert (result[2]) >= 0.;*/
}
void fc3d_NaturalMapABGenerated(
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
    // Not supported in C:
    // Heaviside
    // Heaviside
    // Heaviside
    // Heaviside
    // Heaviside
    // Heaviside
    // Heaviside
    double x13;
    int x14;
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
    /*@ assert (ut1*ut1 + ut2*ut2) >= 0.;*/
    x13=sqrt(ut1*ut1 + ut2*ut2);
    /*@ assert (x13) >= 0.;*/
    x14=x13 <= ZERO;
    int x31;
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
    x31=x13 > ZERO;
    int x71;
    int x72;
    double x40;
    double x55;
    double x56;
    double x57;
    double x59;
    double x70;
    x17=mu*ut1;
    x18=-rt1 + x17;
    x19=x18*x18;
    x20=mu*ut2;
    x21=-rt2 + x20;
    x22=x21*x21;
    x23=x19 + x22;
    /*@ assert (x23) >= 0.;*/
    x24=sqrt(x23);
    x71=x24 > 0;
    x72=x31 && x71;
    int x73;
    int x74;
    x73=x24 <= 0;
    x74=x31 && x73;
    double x41;
    double x42;
    double x119;
    if (x14)
    {
        x1=-1.0*un;
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x8=Heaviside(x1 + x7);
        x9=0.5*x8;
        x10=-1.0*x6;
        x11=Heaviside(x2 + x1 + x10);
        x12=0.5*x11;
    }
    else if (x31)
    {
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x27=0.5*x26;
        x28=x16 - x24;
        x29=Heaviside(x28);
        x30=0.5*x29;
    }
    else if (x72)
    {
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x28=x16 - x24;
        x29=Heaviside(x28);
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x55=rt1 - x17;
        x56=0.5*x26*x40;
        x57=x55*x56;
        x59=0.5*x29*x40;
        x70=x18*x59;
    }
    else if (x71)
    {
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x28=x16 - x24;
        x29=Heaviside(x28);
        x30=0.5*x29;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x42=x18*x41;
        x55=rt1 - x17;
        x119=0.5*mu*x26*x40;
    }
    /* Assignment result[0, 0]=Piecewise((x12 + x9, x14), (x27 + x30, x31)) */

    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x1=-1.0*un;
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x8=Heaviside(x1 + x7);
        x9=0.5*x8;
        x10=-1.0*x6;
        x11=Heaviside(x2 + x1 + x10);
        x12=0.5*x11;

        /* Assignment result[0, 0]=x12 + x9 */
        x1=-1.0*un;
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x8=Heaviside(x1 + x7);
        x9=0.5*x8;
        x10=-1.0*x6;
        x11=Heaviside(x2 + x1 + x10);
        x12=0.5*x11;
        result[0] = x12 + x9;

    }
    else if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x27=0.5*x26;
        x28=x16 - x24;
        x29=Heaviside(x28);
        x30=0.5*x29;

        /* Assignment result[0, 0]=x27 + x30 */
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x27=0.5*x26;
        x28=x16 - x24;
        x29=Heaviside(x28);
        x30=0.5*x29;
        result[0] = x27 + x30;

    }
    /*@ assert (result[0]) >= 0.;*/

    /* Assignment result[1, 0]=Piecewise((-0.5*x32*x66*(x68 - 1.0*x69), x14), (x57 + x70, x72), (0, x74)) */
    double x32;
    double x65;
    double x66;
    double x67;
    double x68;
    double x69;
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x10=-1.0*x6;
        /*@ assert (x6) != 0.;*/
        x32=1.0/x6;
        x65=Heaviside(x6);
        x66=rt1*x65;
        x67=un - 1.0*x2;
        x68=Heaviside(x10 + x67);
        x69=Heaviside(x6 + x67);

        /* Assignment result[1, 0]=-0.5*x32*x66*(x68 - 1.0*x69) */
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x10=-1.0*x6;
        /*@ assert (x6) != 0.;*/
        x32=1.0/x6;
        x65=Heaviside(x6);
        x66=rt1*x65;
        x67=un - 1.0*x2;
        x68=Heaviside(x10 + x67);
        x69=Heaviside(x6 + x67);
        result[1] = -0.5*x32*x66*(x68 - 1.0*x69);

    }
    else if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x28=x16 - x24;
        x29=Heaviside(x28);
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x55=rt1 - x17;
        x56=0.5*x26*x40;
        x57=x55*x56;
        x59=0.5*x29*x40;
        x70=x18*x59;

        /* Assignment result[1, 0]=x57 + x70 */
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x28=x16 - x24;
        x29=Heaviside(x28);
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x55=rt1 - x17;
        x56=0.5*x26*x40;
        x57=x55*x56;
        x59=0.5*x29*x40;
        x70=x18*x59;
        result[1] = x57 + x70;

    }
    else if (x74)
    {
        DEBUG_PRINT("Case (x74) is True.\n");

        /* Assignment result[1, 0]=0 */

        result[1] = 0;
        /*@ assert (result[1]) >= 0.;*/
    }
    /*@ assert (result[1]) >= 0.;*/

    /* Assignment result[2, 0]=Piecewise((x32*(x68*x127 - x69*x127 + x126*x110 - x126*x112 - x6*x12 + x6*x9), x14), (x128 + x62, x72), (x27 - x30, x74)) */
    double x61;
    double x62;
    double x110;
    double x112;
    double x126;
    double x127;
    double x128;
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x1=-1.0*un;
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x8=Heaviside(x1 + x7);
        x9=0.5*x8;
        x10=-1.0*x6;
        x11=Heaviside(x2 + x1 + x10);
        x12=0.5*x11;
        /*@ assert (x6) != 0.;*/
        x32=1.0/x6;
        x65=Heaviside(x6);
        x67=un - 1.0*x2;
        x68=Heaviside(x10 + x67);
        x69=Heaviside(x6 + x67);
        x110=rt2*x69;
        x112=rt2*x68;
        x126=0.5*x65;
        x127=0.5*x6*x65;

        /* Assignment result[2, 0]=x32*(x68*x127 - x69*x127 + x126*x110 - x126*x112 - x6*x12 + x6*x9) */
        x1=-1.0*un;
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x8=Heaviside(x1 + x7);
        x9=0.5*x8;
        x10=-1.0*x6;
        x11=Heaviside(x2 + x1 + x10);
        x12=0.5*x11;
        /*@ assert (x6) != 0.;*/
        x32=1.0/x6;
        x65=Heaviside(x6);
        x67=un - 1.0*x2;
        x68=Heaviside(x10 + x67);
        x69=Heaviside(x6 + x67);
        x110=rt2*x69;
        x112=rt2*x68;
        x126=0.5*x65;
        x127=0.5*x6*x65;
        result[2] = x32*(x68*x127 - x69*x127 + x126*x110 - x126*x112 - x6*x12 + x6*x9);

    }
    else if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x28=x16 - x24;
        x29=Heaviside(x28);
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x56=0.5*x26*x40;
        x59=0.5*x29*x40;
        x61=rt2 - x20;
        x62=x61*x56;
        x128=x21*x59;

        /* Assignment result[2, 0]=x128 + x62 */
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x28=x16 - x24;
        x29=Heaviside(x28);
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x56=0.5*x26*x40;
        x59=0.5*x29*x40;
        x61=rt2 - x20;
        x62=x61*x56;
        x128=x21*x59;
        result[2] = x128 + x62;

    }
    else if (x74)
    {
        DEBUG_PRINT("Case (x74) is True.\n");
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x27=0.5*x26;
        x28=x16 - x24;
        x29=Heaviside(x28);
        x30=0.5*x29;

        /* Assignment result[2, 0]=x27 - x30 */
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x27=0.5*x26;
        x28=x16 - x24;
        x29=Heaviside(x28);
        x30=0.5*x29;
        result[2] = x27 - x30;

    }
    /*@ assert (result[2]) >= 0.;*/

    /* Assignment result[0, 1]=Piecewise((x33*(rt1*x34 - rt1*x35 + x37), x14), (x44 - x46, x31)) */
    double x33;
    double x34;
    double x35;
    double x36;
    double x37;
    double x38;
    double x39;
    double x43;
    double x44;
    double x45;
    double x46;
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x1=-1.0*un;
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x8=Heaviside(x1 + x7);
        x10=-1.0*x6;
        x11=Heaviside(x2 + x1 + x10);
        /*@ assert (x6) != 0.;*/
        x32=1.0/x6;
        x33=0.25*x32*mu;
        x34=2.0*x8;
        x35=2.0*x11;
        x36=1.4142135623730951455*x6;
        x37=x36*x11 + x8*x36;

        /* Assignment result[0, 1]=x33*(rt1*x34 - rt1*x35 + x37) */
        x1=-1.0*un;
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x8=Heaviside(x1 + x7);
        x10=-1.0*x6;
        x11=Heaviside(x2 + x1 + x10);
        /*@ assert (x6) != 0.;*/
        x32=1.0/x6;
        x33=0.25*x32*mu;
        x34=2.0*x8;
        x35=2.0*x11;
        x36=1.4142135623730951455*x6;
        x37=x36*x11 + x8*x36;
        result[3] = x33*(rt1*x34 - rt1*x35 + x37);

    }
    else if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x27=0.5*x26;
        x28=x16 - x24;
        x29=Heaviside(x28);
        x30=0.5*x29;
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        x39=-x17*x38;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x42=x18*x41;
        x43=x39 + x42;
        x44=-x27*x43;
        x45=x39 - x42;
        x46=x45*x30;

        /* Assignment result[0, 1]=x44 - x46 */
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x27=0.5*x26;
        x28=x16 - x24;
        x29=Heaviside(x28);
        x30=0.5*x29;
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        x39=-x17*x38;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x42=x18*x41;
        x43=x39 + x42;
        x44=-x27*x43;
        x45=x39 - x42;
        x46=x45*x30;
        result[3] = x44 - x46;

    }
    /*@ assert (result[3]) >= 0.;*/

    /* Assignment result[1, 1]=Piecewise((-x75*x65*(x94*x68 + x94*x69 - x80*x84 + x80*x89 - x3*x77 - x78*x79 + x78*x82 - x81*x79 + x81*x82 - x87*x86 - x91 + x92*x86 + x93), x14), (-x100*(-mu*x102 + x101) - x43*x57 - x45*x70 - x96*(-mu*x98 + x41), x72), (0, x74)) */
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
    double x89;
    double x90;
    double x91;
    double x92;
    double x93;
    double x94;
    double x95;
    double x96;
    double x97;
    double x98;
    double x99;
    double x100;
    double x101;
    double x102;
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x10=-1.0*x6;
        x15=-un;
        x65=Heaviside(x6);
        x67=un - 1.0*x2;
        x68=Heaviside(x10 + x67);
        x69=Heaviside(x6 + x67);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=0.25*mu/pow(x5, 5.0/2.0);
        /*@ assert (x5) >= 0.;*/
        x76=pow(x5, 3.0/2.0);
        x77=4.0*x76;
        x78=pow(rt1, 5);
        x79=1.4142135623730951455*x69;
        x80=pow(rt2, 4);
        /*@ assert (x80) >= 0.;*/
        x81=rt1*x80;
        x82=1.4142135623730951455*x68;
        x83=Max(0, x15 + x7);
        x84=2.0*x83;
        x85=pow(rt1, 3);
        x86=x85*x4;
        x87=2.8284271247461902909*x69;
        x88=Max(0, x2 + x15 - x6);
        x89=2.0*x88;
        x90=x3*x4;
        x91=x90*x84;
        x92=2.8284271247461902909*x68;
        x93=x90*x89;
        x94=2.0*x3*x76;

        /* Assignment result[1, 1]=-x75*x65*(x94*x68 + x94*x69 - x80*x84 + x80*x89 - x3*x77 - x78*x79 + x78*x82 - x81*x79 + x81*x82 - x87*x86 - x91 + x92*x86 + x93) */
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x10=-1.0*x6;
        x15=-un;
        x65=Heaviside(x6);
        x67=un - 1.0*x2;
        x68=Heaviside(x10 + x67);
        x69=Heaviside(x6 + x67);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=0.25*mu/pow(x5, 5.0/2.0);
        /*@ assert (x5) >= 0.;*/
        x76=pow(x5, 3.0/2.0);
        x77=4.0*x76;
        x78=pow(rt1, 5);
        x79=1.4142135623730951455*x69;
        x80=pow(rt2, 4);
        /*@ assert (x80) >= 0.;*/
        x81=rt1*x80;
        x82=1.4142135623730951455*x68;
        x83=Max(0, x15 + x7);
        x84=2.0*x83;
        x85=pow(rt1, 3);
        x86=x85*x4;
        x87=2.8284271247461902909*x69;
        x88=Max(0, x2 + x15 - x6);
        x89=2.0*x88;
        x90=x3*x4;
        x91=x90*x84;
        x92=2.8284271247461902909*x68;
        x93=x90*x89;
        x94=2.0*x3*x76;
        result[4] = -x75*x65*(x94*x68 + x94*x69 - x80*x84 + x80*x89 - x3*x77 - x78*x79 + x78*x82 - x81*x79 + x81*x82 - x87*x86 - x91 + x92*x86 + x93);

    }
    else if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x28=x16 - x24;
        x29=Heaviside(x28);
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        x39=-x17*x38;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x42=x18*x41;
        x43=x39 + x42;
        x45=x39 - x42;
        x55=rt1 - x17;
        x56=0.5*x26*x40;
        x57=x55*x56;
        x59=0.5*x29*x40;
        x70=x18*x59;
        x95=Max(0, x28);
        x96=0.5*x95;
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x98=x19*x97;
        x99=Max(0, x25);
        x100=0.5*x99;
        x101=-x41;
        x102=x18*x55*x97;

        /* Assignment result[1, 1]=-x100*(-mu*x102 + x101) - x43*x57 - x45*x70 - x96*(-mu*x98 + x41) */
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x28=x16 - x24;
        x29=Heaviside(x28);
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        x39=-x17*x38;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x42=x18*x41;
        x43=x39 + x42;
        x45=x39 - x42;
        x55=rt1 - x17;
        x56=0.5*x26*x40;
        x57=x55*x56;
        x59=0.5*x29*x40;
        x70=x18*x59;
        x95=Max(0, x28);
        x96=0.5*x95;
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x98=x19*x97;
        x99=Max(0, x25);
        x100=0.5*x99;
        x101=-x41;
        x102=x18*x55*x97;
        result[4] = -x100*(-mu*x102 + x101) - x43*x57 - x45*x70 - x96*(-mu*x98 + x41);

    }
    else if (x74)
    {
        DEBUG_PRINT("Case (x74) is True.\n");

        /* Assignment result[1, 1]=0 */

        result[4] = 0;
        /*@ assert (result[4]) >= 0.;*/
    }
    /*@ assert (result[4]) >= 0.;*/

    /* Assignment result[2, 1]=Piecewise((-x75*(-x140*x132 + x140*x133 - x146*x68 - x69*x146 - x108*x66 + x134 - x141*x132 + x141*x133 - x143*x87 + x143*x92 + x144*x145 - x144*x147 - x34*x78 - x66*x103 + x66*x106 + x66*x111 + x66*x113 + x78*x135 - x78*x139 - x78*x142 - x78*x35 + x79*x148 + x81*x135 - x81*x139 - x81*x142 - x81*x34 - x81*x35 - x82*x148 + x86*x136 - x86*x137 - x86*x138), x14), (mu*x149 + x118 - x43*x62 - x45*x128, x72), (x44 + x46, x74)) */
    double x103;
    double x104;
    double x105;
    double x106;
    double x108;
    double x109;
    double x111;
    double x113;
    double x115;
    double x117;
    double x118;
    double x125;
    double x129;
    double x130;
    double x131;
    double x132;
    double x133;
    double x134;
    double x135;
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
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x1=-1.0*un;
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x8=Heaviside(x1 + x7);
        x10=-1.0*x6;
        x11=Heaviside(x2 + x1 + x10);
        x15=-un;
        x34=2.0*x8;
        x35=2.0*x11;
        x65=Heaviside(x6);
        x66=rt1*x65;
        x67=un - 1.0*x2;
        x68=Heaviside(x10 + x67);
        x69=Heaviside(x6 + x67);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=0.25*mu/pow(x5, 5.0/2.0);
        /*@ assert (x5) >= 0.;*/
        x76=pow(x5, 3.0/2.0);
        x77=4.0*x76;
        x78=pow(rt1, 5);
        x79=1.4142135623730951455*x69;
        x80=pow(rt2, 4);
        /*@ assert (x80) >= 0.;*/
        x81=rt1*x80;
        x82=1.4142135623730951455*x68;
        x83=Max(0, x15 + x7);
        x84=2.0*x83;
        x85=pow(rt1, 3);
        x86=x85*x4;
        x87=2.8284271247461902909*x69;
        x88=Max(0, x2 + x15 - x6);
        x89=2.0*x88;
        x92=2.8284271247461902909*x68;
        x103=rt2*x77;
        x104=pow(rt1, 4);
        /*@ assert (x104) >= 0.;*/
        x105=pow(rt2, 3);
        x106=x105*x84;
        x108=x105*x89;
        x109=2.0*x76;
        x110=rt2*x69;
        x111=x109*x110;
        x112=rt2*x68;
        x113=x109*x112;
        x129=1.4142135623730951455*x8*x76;
        x130=1.4142135623730951455*x76*x11;
        x131=x3*x76;
        x132=1.4142135623730951455*x65*x69;
        x133=1.4142135623730951455*x65*x68;
        x134=x131*x132 - x131*x133 - x3*x129 + x3*x130 - x4*x129 + x4*x130;
        x135=4.0*x65;
        x136=8.0*x65;
        x137=4.0*x8;
        x138=4.0*x11;
        x139=2.0*x65*x69;
        x140=pow(rt2, 5);
        x141=rt2*x104;
        x142=2.0*x65*x68;
        x143=x3*x105*x65;
        x144=rt2*x85;
        x145=2.0*x65*x83;
        x146=4.0*x85*x4*x65;
        x147=2.0*x65*x88;
        x148=x4*x76*x65;

        /* Assignment result[2, 1]=-x75*(-x140*x132 + x140*x133 - x146*x68 - x69*x146 - x108*x66 + x134 - x141*x132 + x141*x133 - x143*x87 + x143*x92 + x144*x145 - x144*x147 - x34*x78 - x66*x103 + x66*x106 + x66*x111 + x66*x113 + x78*x135 - x78*x139 - x78*x142 - x78*x35 + x79*x148 + x81*x135 - x81*x139 - x81*x142 - x81*x34 - x81*x35 - x82*x148 + x86*x136 - x86*x137 - x86*x138) */
        x1=-1.0*un;
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x8=Heaviside(x1 + x7);
        x10=-1.0*x6;
        x11=Heaviside(x2 + x1 + x10);
        x15=-un;
        x34=2.0*x8;
        x35=2.0*x11;
        x65=Heaviside(x6);
        x66=rt1*x65;
        x67=un - 1.0*x2;
        x68=Heaviside(x10 + x67);
        x69=Heaviside(x6 + x67);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=0.25*mu/pow(x5, 5.0/2.0);
        /*@ assert (x5) >= 0.;*/
        x76=pow(x5, 3.0/2.0);
        x77=4.0*x76;
        x78=pow(rt1, 5);
        x79=1.4142135623730951455*x69;
        x80=pow(rt2, 4);
        /*@ assert (x80) >= 0.;*/
        x81=rt1*x80;
        x82=1.4142135623730951455*x68;
        x83=Max(0, x15 + x7);
        x84=2.0*x83;
        x85=pow(rt1, 3);
        x86=x85*x4;
        x87=2.8284271247461902909*x69;
        x88=Max(0, x2 + x15 - x6);
        x89=2.0*x88;
        x92=2.8284271247461902909*x68;
        x103=rt2*x77;
        x104=pow(rt1, 4);
        /*@ assert (x104) >= 0.;*/
        x105=pow(rt2, 3);
        x106=x105*x84;
        x108=x105*x89;
        x109=2.0*x76;
        x110=rt2*x69;
        x111=x109*x110;
        x112=rt2*x68;
        x113=x109*x112;
        x129=1.4142135623730951455*x8*x76;
        x130=1.4142135623730951455*x76*x11;
        x131=x3*x76;
        x132=1.4142135623730951455*x65*x69;
        x133=1.4142135623730951455*x65*x68;
        x134=x131*x132 - x131*x133 - x3*x129 + x3*x130 - x4*x129 + x4*x130;
        x135=4.0*x65;
        x136=8.0*x65;
        x137=4.0*x8;
        x138=4.0*x11;
        x139=2.0*x65*x69;
        x140=pow(rt2, 5);
        x141=rt2*x104;
        x142=2.0*x65*x68;
        x143=x3*x105*x65;
        x144=rt2*x85;
        x145=2.0*x65*x83;
        x146=4.0*x85*x4*x65;
        x147=2.0*x65*x88;
        x148=x4*x76*x65;
        result[5] = -x75*(-x140*x132 + x140*x133 - x146*x68 - x69*x146 - x108*x66 + x134 - x141*x132 + x141*x133 - x143*x87 + x143*x92 + x144*x145 - x144*x147 - x34*x78 - x66*x103 + x66*x106 + x66*x111 + x66*x113 + x78*x135 - x78*x139 - x78*x142 - x78*x35 + x79*x148 + x81*x135 - x81*x139 - x81*x142 - x81*x34 - x81*x35 - x82*x148 + x86*x136 - x86*x137 - x86*x138);

    }
    else if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x28=x16 - x24;
        x29=Heaviside(x28);
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        x39=-x17*x38;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x42=x18*x41;
        x43=x39 + x42;
        x45=x39 - x42;
        x56=0.5*x26*x40;
        x59=0.5*x29*x40;
        x61=rt2 - x20;
        x62=x61*x56;
        x95=Max(0, x28);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x99=Max(0, x25);
        x115=0.5*x97*x99;
        x117=0.5*x18*x95*x21*x97;
        x118=mu*x117;
        x125=x18*x61;
        x128=x21*x59;
        x149=x125*x115;

        /* Assignment result[2, 1]=mu*x149 + x118 - x43*x62 - x45*x128 */
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x28=x16 - x24;
        x29=Heaviside(x28);
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        x39=-x17*x38;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x42=x18*x41;
        x43=x39 + x42;
        x45=x39 - x42;
        x56=0.5*x26*x40;
        x59=0.5*x29*x40;
        x61=rt2 - x20;
        x62=x61*x56;
        x95=Max(0, x28);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x99=Max(0, x25);
        x115=0.5*x97*x99;
        x117=0.5*x18*x95*x21*x97;
        x118=mu*x117;
        x125=x18*x61;
        x128=x21*x59;
        x149=x125*x115;
        result[5] = mu*x149 + x118 - x43*x62 - x45*x128;

    }
    else if (x74)
    {
        DEBUG_PRINT("Case (x74) is True.\n");
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x27=0.5*x26;
        x28=x16 - x24;
        x29=Heaviside(x28);
        x30=0.5*x29;
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        x39=-x17*x38;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x42=x18*x41;
        x43=x39 + x42;
        x44=-x27*x43;
        x45=x39 - x42;
        x46=x45*x30;

        /* Assignment result[2, 1]=x44 + x46 */
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x27=0.5*x26;
        x28=x16 - x24;
        x29=Heaviside(x28);
        x30=0.5*x29;
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        x39=-x17*x38;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x42=x18*x41;
        x43=x39 + x42;
        x44=-x27*x43;
        x45=x39 - x42;
        x46=x45*x30;
        result[5] = x44 + x46;

    }
    /*@ assert (result[5]) >= 0.;*/

    /* Assignment result[0, 2]=Piecewise((x33*(rt2*x34 - rt2*x35 + x37), x14), (x50 - x52, x31)) */
    double x47;
    double x48;
    double x49;
    double x50;
    double x51;
    double x52;
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x1=-1.0*un;
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x8=Heaviside(x1 + x7);
        x10=-1.0*x6;
        x11=Heaviside(x2 + x1 + x10);
        /*@ assert (x6) != 0.;*/
        x32=1.0/x6;
        x33=0.25*x32*mu;
        x34=2.0*x8;
        x35=2.0*x11;
        x36=1.4142135623730951455*x6;
        x37=x36*x11 + x8*x36;

        /* Assignment result[0, 2]=x33*(rt2*x34 - rt2*x35 + x37) */
        x1=-1.0*un;
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x8=Heaviside(x1 + x7);
        x10=-1.0*x6;
        x11=Heaviside(x2 + x1 + x10);
        /*@ assert (x6) != 0.;*/
        x32=1.0/x6;
        x33=0.25*x32*mu;
        x34=2.0*x8;
        x35=2.0*x11;
        x36=1.4142135623730951455*x6;
        x37=x36*x11 + x8*x36;
        result[6] = x33*(rt2*x34 - rt2*x35 + x37);

    }
    else if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x27=0.5*x26;
        x28=x16 - x24;
        x29=Heaviside(x28);
        x30=0.5*x29;
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x47=-x20*x38;
        x48=x21*x41;
        x49=x47 + x48;
        x50=-x27*x49;
        x51=x47 - x48;
        x52=x51*x30;

        /* Assignment result[0, 2]=x50 - x52 */
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x27=0.5*x26;
        x28=x16 - x24;
        x29=Heaviside(x28);
        x30=0.5*x29;
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x47=-x20*x38;
        x48=x21*x41;
        x49=x47 + x48;
        x50=-x27*x49;
        x51=x47 - x48;
        x52=x51*x30;
        result[6] = x50 - x52;

    }
    /*@ assert (result[6]) >= 0.;*/

    /* Assignment result[1, 2]=Piecewise((-rt1*x75*x65*(-x104*x79 + x104*x82 - x80*x79 + x80*x82 - x103 + x106 + x107*x84 - x107*x89 - x108 + x111 + x113 - x90*x87 + x90*x92), x14), (mu*x116 + x118 - x49*x57 - x51*x70, x72), (0, x74)) */
    double x107;
    double x114;
    double x116;
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x10=-1.0*x6;
        x15=-un;
        x65=Heaviside(x6);
        x67=un - 1.0*x2;
        x68=Heaviside(x10 + x67);
        x69=Heaviside(x6 + x67);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=0.25*mu/pow(x5, 5.0/2.0);
        /*@ assert (x5) >= 0.;*/
        x76=pow(x5, 3.0/2.0);
        x77=4.0*x76;
        x79=1.4142135623730951455*x69;
        x80=pow(rt2, 4);
        /*@ assert (x80) >= 0.;*/
        x82=1.4142135623730951455*x68;
        x83=Max(0, x15 + x7);
        x84=2.0*x83;
        x87=2.8284271247461902909*x69;
        x88=Max(0, x2 + x15 - x6);
        x89=2.0*x88;
        x90=x3*x4;
        x92=2.8284271247461902909*x68;
        x103=rt2*x77;
        x104=pow(rt1, 4);
        /*@ assert (x104) >= 0.;*/
        x105=pow(rt2, 3);
        x106=x105*x84;
        x107=rt2*x3;
        x108=x105*x89;
        x109=2.0*x76;
        x110=rt2*x69;
        x111=x109*x110;
        x112=rt2*x68;
        x113=x109*x112;

        /* Assignment result[1, 2]=-rt1*x75*x65*(-x104*x79 + x104*x82 - x80*x79 + x80*x82 - x103 + x106 + x107*x84 - x107*x89 - x108 + x111 + x113 - x90*x87 + x90*x92) */
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x10=-1.0*x6;
        x15=-un;
        x65=Heaviside(x6);
        x67=un - 1.0*x2;
        x68=Heaviside(x10 + x67);
        x69=Heaviside(x6 + x67);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=0.25*mu/pow(x5, 5.0/2.0);
        /*@ assert (x5) >= 0.;*/
        x76=pow(x5, 3.0/2.0);
        x77=4.0*x76;
        x79=1.4142135623730951455*x69;
        x80=pow(rt2, 4);
        /*@ assert (x80) >= 0.;*/
        x82=1.4142135623730951455*x68;
        x83=Max(0, x15 + x7);
        x84=2.0*x83;
        x87=2.8284271247461902909*x69;
        x88=Max(0, x2 + x15 - x6);
        x89=2.0*x88;
        x90=x3*x4;
        x92=2.8284271247461902909*x68;
        x103=rt2*x77;
        x104=pow(rt1, 4);
        /*@ assert (x104) >= 0.;*/
        x105=pow(rt2, 3);
        x106=x105*x84;
        x107=rt2*x3;
        x108=x105*x89;
        x109=2.0*x76;
        x110=rt2*x69;
        x111=x109*x110;
        x112=rt2*x68;
        x113=x109*x112;
        result[7] = -rt1*x75*x65*(-x104*x79 + x104*x82 - x80*x79 + x80*x82 - x103 + x106 + x107*x84 - x107*x89 - x108 + x111 + x113 - x90*x87 + x90*x92);

    }
    else if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x28=x16 - x24;
        x29=Heaviside(x28);
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x47=-x20*x38;
        x48=x21*x41;
        x49=x47 + x48;
        x51=x47 - x48;
        x55=rt1 - x17;
        x56=0.5*x26*x40;
        x57=x55*x56;
        x59=0.5*x29*x40;
        x70=x18*x59;
        x95=Max(0, x28);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x99=Max(0, x25);
        x114=x55*x21;
        x115=0.5*x97*x99;
        x116=x114*x115;
        x117=0.5*x18*x95*x21*x97;
        x118=mu*x117;

        /* Assignment result[1, 2]=mu*x116 + x118 - x49*x57 - x51*x70 */
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x28=x16 - x24;
        x29=Heaviside(x28);
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x47=-x20*x38;
        x48=x21*x41;
        x49=x47 + x48;
        x51=x47 - x48;
        x55=rt1 - x17;
        x56=0.5*x26*x40;
        x57=x55*x56;
        x59=0.5*x29*x40;
        x70=x18*x59;
        x95=Max(0, x28);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x99=Max(0, x25);
        x114=x55*x21;
        x115=0.5*x97*x99;
        x116=x114*x115;
        x117=0.5*x18*x95*x21*x97;
        x118=mu*x117;
        result[7] = mu*x116 + x118 - x49*x57 - x51*x70;

    }
    else if (x74)
    {
        DEBUG_PRINT("Case (x74) is True.\n");

        /* Assignment result[1, 2]=0 */

        result[7] = 0;
        /*@ assert (result[7]) >= 0.;*/
    }
    /*@ assert (result[7]) >= 0.;*/

    /* Assignment result[2, 2]=Piecewise((-x75*(-x104*x145 + x104*x147 - 1.17157287525381*x3*x105*x65*x68 - 6.82842712474619*x3*x105*x65*x69 + x140*x135 - x140*x151 - x140*x152 - x34*x140 - x140*x35 - x77*x4*x65 - x65*x91 + x65*x93 + 0.585786437626905*x68*x148 + 3.41421356237309*x69*x148 + x134 + x141*x135 - x141*x151 - x141*x152 - x141*x35 + x150*x136 - x150*x137 - x150*x138 - x34*x141), x14), (-x100*(-mu*x154 + x101) - x49*x62 - x51*x128 - x96*(-mu*x153 + x41), x72), (x50 + x52, x74)) */
    double x150;
    double x151;
    double x152;
    double x153;
    double x154;
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x1=-1.0*un;
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x8=Heaviside(x1 + x7);
        x10=-1.0*x6;
        x11=Heaviside(x2 + x1 + x10);
        x15=-un;
        x34=2.0*x8;
        x35=2.0*x11;
        x65=Heaviside(x6);
        x67=un - 1.0*x2;
        x68=Heaviside(x10 + x67);
        x69=Heaviside(x6 + x67);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=0.25*mu/pow(x5, 5.0/2.0);
        /*@ assert (x5) >= 0.;*/
        x76=pow(x5, 3.0/2.0);
        x77=4.0*x76;
        x83=Max(0, x15 + x7);
        x84=2.0*x83;
        x88=Max(0, x2 + x15 - x6);
        x89=2.0*x88;
        x90=x3*x4;
        x91=x90*x84;
        x93=x90*x89;
        x104=pow(rt1, 4);
        /*@ assert (x104) >= 0.;*/
        x105=pow(rt2, 3);
        x129=1.4142135623730951455*x8*x76;
        x130=1.4142135623730951455*x76*x11;
        x131=x3*x76;
        x132=1.4142135623730951455*x65*x69;
        x133=1.4142135623730951455*x65*x68;
        x134=x131*x132 - x131*x133 - x3*x129 + x3*x130 - x4*x129 + x4*x130;
        x135=4.0*x65;
        x136=8.0*x65;
        x137=4.0*x8;
        x138=4.0*x11;
        x140=pow(rt2, 5);
        x141=rt2*x104;
        x145=2.0*x65*x83;
        x147=2.0*x65*x88;
        x148=x4*x76*x65;
        x150=x3*x105;
        x151=3.4142135623730949234*x65*x69;
        x152=0.58578643762690496555*x65*x68;

        /* Assignment result[2, 2]=-x75*(-x104*x145 + x104*x147 - 1.17157287525381*x3*x105*x65*x68 - 6.82842712474619*x3*x105*x65*x69 + x140*x135 - x140*x151 - x140*x152 - x34*x140 - x140*x35 - x77*x4*x65 - x65*x91 + x65*x93 + 0.585786437626905*x68*x148 + 3.41421356237309*x69*x148 + x134 + x141*x135 - x141*x151 - x141*x152 - x141*x35 + x150*x136 - x150*x137 - x150*x138 - x34*x141) */
        x1=-1.0*un;
        x2=mu*rn;
        x3=rt1*rt1;
        x4=rt2*rt2;
        x5=x3 + x4;
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5);
        x7=x2 + x6;
        x8=Heaviside(x1 + x7);
        x10=-1.0*x6;
        x11=Heaviside(x2 + x1 + x10);
        x15=-un;
        x34=2.0*x8;
        x35=2.0*x11;
        x65=Heaviside(x6);
        x67=un - 1.0*x2;
        x68=Heaviside(x10 + x67);
        x69=Heaviside(x6 + x67);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=0.25*mu/pow(x5, 5.0/2.0);
        /*@ assert (x5) >= 0.;*/
        x76=pow(x5, 3.0/2.0);
        x77=4.0*x76;
        x83=Max(0, x15 + x7);
        x84=2.0*x83;
        x88=Max(0, x2 + x15 - x6);
        x89=2.0*x88;
        x90=x3*x4;
        x91=x90*x84;
        x93=x90*x89;
        x104=pow(rt1, 4);
        /*@ assert (x104) >= 0.;*/
        x105=pow(rt2, 3);
        x129=1.4142135623730951455*x8*x76;
        x130=1.4142135623730951455*x76*x11;
        x131=x3*x76;
        x132=1.4142135623730951455*x65*x69;
        x133=1.4142135623730951455*x65*x68;
        x134=x131*x132 - x131*x133 - x3*x129 + x3*x130 - x4*x129 + x4*x130;
        x135=4.0*x65;
        x136=8.0*x65;
        x137=4.0*x8;
        x138=4.0*x11;
        x140=pow(rt2, 5);
        x141=rt2*x104;
        x145=2.0*x65*x83;
        x147=2.0*x65*x88;
        x148=x4*x76*x65;
        x150=x3*x105;
        x151=3.4142135623730949234*x65*x69;
        x152=0.58578643762690496555*x65*x68;
        result[8] = -x75*(-x104*x145 + x104*x147 - 1.1715728752538099311*x3*x105*x65*x68 - 6.8284271247461898469*x3*x105*x65*x69 + x140*x135 - x140*x151 - x140*x152 - x34*x140 - x140*x35 - x77*x4*x65 - x65*x91 + x65*x93 + 0.58578643762690496555*x68*x148 + 3.4142135623730949234*x69*x148 + x134 + x141*x135 - x141*x151 - x141*x152 - x141*x35 + x150*x136 - x150*x137 - x150*x138 - x34*x141);

    }
    else if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x28=x16 - x24;
        x29=Heaviside(x28);
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x47=-x20*x38;
        x48=x21*x41;
        x49=x47 + x48;
        x51=x47 - x48;
        x56=0.5*x26*x40;
        x59=0.5*x29*x40;
        x61=rt2 - x20;
        x62=x61*x56;
        x95=Max(0, x28);
        x96=0.5*x95;
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x99=Max(0, x25);
        x100=0.5*x99;
        x101=-x41;
        x128=x21*x59;
        x153=x22*x97;
        x154=x61*x21*x97;

        /* Assignment result[2, 2]=-x100*(-mu*x154 + x101) - x49*x62 - x51*x128 - x96*(-mu*x153 + x41) */
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x28=x16 - x24;
        x29=Heaviside(x28);
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x47=-x20*x38;
        x48=x21*x41;
        x49=x47 + x48;
        x51=x47 - x48;
        x56=0.5*x26*x40;
        x59=0.5*x29*x40;
        x61=rt2 - x20;
        x62=x61*x56;
        x95=Max(0, x28);
        x96=0.5*x95;
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x99=Max(0, x25);
        x100=0.5*x99;
        x101=-x41;
        x128=x21*x59;
        x153=x22*x97;
        x154=x61*x21*x97;
        result[8] = -x100*(-mu*x154 + x101) - x49*x62 - x51*x128 - x96*(-mu*x153 + x41);

    }
    else if (x74)
    {
        DEBUG_PRINT("Case (x74) is True.\n");
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x27=0.5*x26;
        x28=x16 - x24;
        x29=Heaviside(x28);
        x30=0.5*x29;
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x47=-x20*x38;
        x48=x21*x41;
        x49=x47 + x48;
        x50=-x27*x49;
        x51=x47 - x48;
        x52=x51*x30;

        /* Assignment result[2, 2]=x50 + x52 */
        x2=mu*rn;
        x15=-un;
        x16=-mu*x13 + x2 + x15;
        x25=x16 + x24;
        x26=Heaviside(x25);
        x27=0.5*x26;
        x28=x16 - x24;
        x29=Heaviside(x28);
        x30=0.5*x29;
        /*@ assert (x13) != 0.;*/
        x38=1.0/x13;
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x47=-x20*x38;
        x48=x21*x41;
        x49=x47 + x48;
        x50=-x27*x49;
        x51=x47 - x48;
        x52=x51*x30;
        result[8] = x50 + x52;

    }
    /*@ assert (result[8]) >= 0.;*/

    /* Assignment result[0, 3]=mu + x53 - x54 */
    double x53;
    double x54;x2=mu*rn;
    x15=-un;
    x16=-mu*x13 + x2 + x15;
    x25=x16 + x24;
    x26=Heaviside(x25);
    x27=0.5*x26;
    x28=x16 - x24;
    x29=Heaviside(x28);
    x30=0.5*x29;
    x53=-mu*x27;
    x54=mu*x30;
    result[9] = mu + x53 - x54;


    /* Assignment result[1, 3]=Piecewise((-x42*x30 - x55*x119, x71), (0, x73)) */

    if (x71)
    {
        DEBUG_PRINT("Case (x71) is True.\n");
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x42=x18*x41;
        x55=rt1 - x17;
        x119=0.5*mu*x26*x40;

        /* Assignment result[1, 3]=-x42*x30 - x55*x119 */
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x42=x18*x41;
        x55=rt1 - x17;
        x119=0.5*mu*x26*x40;
        result[10] = -x42*x30 - x55*x119;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[1, 3]=0 */

        result[10] = 0;
        /*@ assert (result[10]) >= 0.;*/
    }
    /*@ assert (result[10]) >= 0.;*/

    /* Assignment result[2, 3]=Piecewise((-x48*x30 - x61*x119, x71), (x53 + x54, x73)) */

    if (x71)
    {
        DEBUG_PRINT("Case (x71) is True.\n");
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x48=x21*x41;
        x61=rt2 - x20;
        x119=0.5*mu*x26*x40;

        /* Assignment result[2, 3]=-x48*x30 - x61*x119 */
        /*@ assert (x24) != 0.;*/
        x40=1.0/x24;
        x41=mu*x40;
        x48=x21*x41;
        x61=rt2 - x20;
        x119=0.5*mu*x26*x40;
        result[11] = -x48*x30 - x61*x119;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[2, 3]=x53 + x54 */

        result[11] = x53 + x54;

    }
    /*@ assert (result[11]) >= 0.;*/

    /* Assignment result[0, 4]=x58 + x60 */
    double x58;
    double x60;/*@ assert (x24) != 0.;*/
    x40=1.0/x24;
    x55=rt1 - x17;
    x56=0.5*x26*x40;
    x57=x55*x56;
    x58=-x57;
    x59=0.5*x29*x40;
    x60=x55*x59;
    result[12] = x58 + x60;


    /* Assignment result[1, 4]=Piecewise((1 - x100*(x102 + x40) + x18*x55*x122 - x55**2*x121 - x96*(x123 + x98), x71), (1, x73)) */
    double x120;
    double x121;
    double x122;
    double x123;
    if (x71)
    {
        DEBUG_PRINT("Case (x71) is True.\n");
        x95=Max(0, x28);
        x96=0.5*x95;
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x98=x19*x97;
        x99=Max(0, x25);
        x100=0.5*x99;
        x102=x18*x55*x97;
        /*@ assert (x23) != 0.;*/
        x120=1.0/x23;
        x121=0.5*x26*x120;
        x122=0.5*x29*x120;
        x123=-x40;

        /* Assignment result[1, 4]=1 - x100*(x102 + x40) + x18*x55*x122 - x55**2*x121 - x96*(x123 + x98) */
        x95=Max(0, x28);
        x96=0.5*x95;
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x98=x19*x97;
        x99=Max(0, x25);
        x100=0.5*x99;
        x102=x18*x55*x97;
        /*@ assert (x23) != 0.;*/
        x120=1.0/x23;
        x121=0.5*x26*x120;
        x122=0.5*x29*x120;
        x123=-x40;
        result[13] = 1 - x100*(x102 + x40) + x18*x55*x122 - x55*x55*x121 - x96*(x123 + x98);

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[1, 4]=1 */

        result[13] = 1;
        /*@ assert (result[13]) >= 0.;*/
        /*@ assert (result[13]) != 0.;*/
    }
    /*@ assert (result[13]) >= 0.;*/

    /* Assignment result[2, 4]=Piecewise((x114*x122 + x124 - x149, x71), (x58 - x60, x73)) */
    double x124;
    if (x71)
    {
        DEBUG_PRINT("Case (x71) is True.\n");
        x61=rt2 - x20;
        x95=Max(0, x28);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x99=Max(0, x25);
        x114=x55*x21;
        x115=0.5*x97*x99;
        x117=0.5*x18*x95*x21*x97;
        /*@ assert (x23) != 0.;*/
        x120=1.0/x23;
        x121=0.5*x26*x120;
        x122=0.5*x29*x120;
        x124=-x117 - x55*x61*x121;
        x125=x18*x61;
        x149=x125*x115;

        /* Assignment result[2, 4]=x114*x122 + x124 - x149 */
        x61=rt2 - x20;
        x95=Max(0, x28);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x99=Max(0, x25);
        x114=x55*x21;
        x115=0.5*x97*x99;
        x117=0.5*x18*x95*x21*x97;
        /*@ assert (x23) != 0.;*/
        x120=1.0/x23;
        x121=0.5*x26*x120;
        x122=0.5*x29*x120;
        x124=-x117 - x55*x61*x121;
        x125=x18*x61;
        x149=x125*x115;
        result[14] = x114*x122 + x124 - x149;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[2, 4]=x58 - x60 */

        result[14] = x58 - x60;

    }
    /*@ assert (result[14]) >= 0.;*/

    /* Assignment result[0, 5]=x63 + x64 */
    double x63;
    double x64;x61=rt2 - x20;
    x62=x61*x56;
    x63=-x62;
    x64=x61*x59;
    result[15] = x63 + x64;


    /* Assignment result[1, 5]=Piecewise((-x116 + x124 + x125*x122, x71), (0, x73)) */

    if (x71)
    {
        DEBUG_PRINT("Case (x71) is True.\n");
        x95=Max(0, x28);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x99=Max(0, x25);
        x114=x55*x21;
        x115=0.5*x97*x99;
        x116=x114*x115;
        x117=0.5*x18*x95*x21*x97;
        /*@ assert (x23) != 0.;*/
        x120=1.0/x23;
        x121=0.5*x26*x120;
        x122=0.5*x29*x120;
        x124=-x117 - x55*x61*x121;
        x125=x18*x61;

        /* Assignment result[1, 5]=-x116 + x124 + x125*x122 */
        x95=Max(0, x28);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x99=Max(0, x25);
        x114=x55*x21;
        x115=0.5*x97*x99;
        x116=x114*x115;
        x117=0.5*x18*x95*x21*x97;
        /*@ assert (x23) != 0.;*/
        x120=1.0/x23;
        x121=0.5*x26*x120;
        x122=0.5*x29*x120;
        x124=-x117 - x55*x61*x121;
        x125=x18*x61;
        result[16] = -x116 + x124 + x125*x122;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[1, 5]=0 */

        result[16] = 0;
        /*@ assert (result[16]) >= 0.;*/
    }
    /*@ assert (result[16]) >= 0.;*/

    /* Assignment result[2, 5]=Piecewise((1 - x100*(x154 + x40) + x61*x21*x122 - x61**2*x121 - x96*(x123 + x153), x71), (1 + x63 - x64, x73)) */

    if (x71)
    {
        DEBUG_PRINT("Case (x71) is True.\n");
        x95=Max(0, x28);
        x96=0.5*x95;
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x99=Max(0, x25);
        x100=0.5*x99;
        /*@ assert (x23) != 0.;*/
        x120=1.0/x23;
        x121=0.5*x26*x120;
        x122=0.5*x29*x120;
        x123=-x40;
        x153=x22*x97;
        x154=x61*x21*x97;

        /* Assignment result[2, 5]=1 - x100*(x154 + x40) + x61*x21*x122 - x61**2*x121 - x96*(x123 + x153) */
        x95=Max(0, x28);
        x96=0.5*x95;
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x97=pow(x23, -3.0/2.0);
        x99=Max(0, x25);
        x100=0.5*x99;
        /*@ assert (x23) != 0.;*/
        x120=1.0/x23;
        x121=0.5*x26*x120;
        x122=0.5*x29*x120;
        x123=-x40;
        x153=x22*x97;
        x154=x61*x21*x97;
        result[17] = 1 - x100*(x154 + x40) + x61*x21*x122 - x61*x61*x121 - x96*(x123 + x153);

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[2, 5]=1 + x63 - x64 */

        result[17] = 1 + x63 - x64;

    }
    /*@ assert (result[17]) >= 0.;*/
}

void fc3d_NaturalMapFunctionGenerated(
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

    fc3d_NaturalMapFABGenerated(
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
      fc3d_NaturalMapFGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3(result, f);
    }

    if (A && B)
    {
      fc3d_NaturalMapABGenerated(
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
