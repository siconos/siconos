#include <math.h>
#include <assert.h>
#include <op3x3.h>
#include <stdlib.h>
#include <stdio.h>
#include "SiconosCompat.h"

//#define DEBUG_MESSAGES 1
//#define DEBUG_WHERE_MESSAGES 1
//#include <stdio.h>
#include <debug.h>
#include "FischerBurmeisterGenerated.h"

#define Sign(x) ((x>0) - (x<0))
#define Max fmax
#define Heaviside(x) (.5*Sign(x) + .5)
#define Rand(x) ((double) rand()/ (double) RAND_MAX)

/* Round-off errors: for values supposed to be >=0,  value = 0 is replaced by value <= ZERO */
/* Note:
 * ZERO=DBL_EPSILON*100 (too small) => test-3D_NSN_FB-NESpheres_30_1 fails (division by 0)
 * ZERO=1e-10 (too big) => test-3D_NSN_FB-OneObject-i100000-499.hdf5 fails (no convergence)
 */
#define ZERO DBL_EPSILON*500

#define ZERO_SQR  ZERO * ZERO

#ifdef FB_DEBUG_POST_CHECK
#define POST_CHECK_POW(x) assert(isfinite((double)x))
#define POST_CHECK_ADD(x) assert(isfinite((double)x))
#define POST_CHECK_MUL(x) assert(isfinite((double)x))
#define POST_CHECK(x) assert(isfinite((double)x))
#else
#define POST_CHECK_POW(x)
#define POST_CHECK_ADD(x)
#define POST_CHECK_MUL(x)
#define POST_CHECK(x)
#endif

#define NOT_ZERO(x) fabs(x) > 0
#define IS_NOT_ZERO(x) fabs(x) > 0
#define IS_POSITIVE(x) 1

#define random1 sqrt(2)/2
#define random2 sqrt(2)/2

#ifdef __cplusplus
#include <cmath>
#define CHECK(x)
#else
#define CHECK(x)
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
#define sqrt(x) ((( (double)x >= 0. && x <= ZERO_SQR) || ((double)x < 0. && (double)x >= -ZERO)) ? 0. : (assert(x>=0.), sqrt(x)))

// ./fb2.py --ccode --ccodefac --ccodeAB --wrapper --merit
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
    double x5;
    double x6;
    double x7;
    double x8;
    double x9;
    double x10;
    double x11;
    double x12;
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
    int x37;
    double x38;
    int x39;
    double x28;
    double x29;
    double x30;
    double x31;
    double x32;
    double x33;
    double x34;
    double x35;
    double x36;
    x1=ut1*ut1; POST_CHECK_POW(x1);
    /*@ assert (x1) >= 0.;*/
    x2=ut2*ut2; POST_CHECK_POW(x2);
    /*@ assert (x2) >= 0.;*/
    /*@ assert (x1 + x2) >= 0.;*/
    x3=sqrt(x1 + x2); POST_CHECK_POW(x3);
    /*@ assert (x3) >= 0.;*/
    x4=x3*mu + un; POST_CHECK_ADD(x4);
    x5=mu*rn; POST_CHECK_MUL(x5);
    /*@ assert (x5) >= 0.;*/
    /*@ assert (x5) != 0.;*/
    x6=rt1*rt1; POST_CHECK_POW(x6);
    /*@ assert (x6) >= 0.;*/
    x7=rt2*rt2; POST_CHECK_POW(x7);
    /*@ assert (x7) >= 0.;*/
    x8=mu*mu; POST_CHECK_POW(x8);
    x9=rn*rn; POST_CHECK_POW(x9);
    /*@ assert (x9) != 0.;*/
    x10=x8*x9; POST_CHECK_MUL(x10);
    x11=x1*x8; POST_CHECK_MUL(x11);
    x12=x2*x8; POST_CHECK_MUL(x12);
    x13=x4*x4 + x6 + x7 + x10 + x11 + x12; POST_CHECK_ADD(x13);
    x14=mu*ut1; POST_CHECK_MUL(x14);
    x15=x5*rt1 + x14*x4; POST_CHECK_ADD(x15);
    x16=x15*x15; POST_CHECK_POW(x16);
    x17=mu*ut2; POST_CHECK_MUL(x17);
    x18=x5*rt2 + x17*x4; POST_CHECK_ADD(x18);
    x19=x18*x18; POST_CHECK_POW(x19);
    x20=x16 + x19; POST_CHECK_ADD(x20);
    /*@ assert (x20) >= 0.;*/
    x21=sqrt(x20); POST_CHECK_POW(x21);
    x22=2*x21; POST_CHECK_MUL(x22);
    x23=x13 - x22; POST_CHECK_ADD(x23);
    x37=x21 <= ZERO; POST_CHECK(x37);
    x38=fabs(x23); POST_CHECK(x38);
    /*@ assert (x38) >= 0.;*/
    x39=x37 || x38 <= ZERO; POST_CHECK(x39);
    int x49;
    double x40;
    double x41;
    double x42;
    double x43;
    double x44;
    double x45;
    double x46;
    double x47;
    double x48;
    x49=x3 <= ZERO; POST_CHECK(x49);
    int x58;
    int x59;
    int x60;
    int x61;
    double x24;
    double x26;
    double x50;
    double x51;
    double x52;
    double x53;
    double x54;
    double x55;
    double x56;
    double x57;
    x58=x3 > 0; POST_CHECK(x58);
    x59=x21 > 0; POST_CHECK(x59);
    x60=x38 > 0; POST_CHECK(x60);
    x61=x58 && x59 && x60; POST_CHECK(x61);
    int x100;
    double x92;
    double x93;
    double x94;
    double x95;
    double x96;
    double x97;
    double x98;
    double x99;
    x100=x59 && x60; POST_CHECK(x100);
    double x25;
    double x27;
    double x111;
    double x112;
    double x113;
    double x114;
    int x125;
    x125=x37 && x58 && x59 && x60; POST_CHECK(x125);
    int x177;
    x177=x37 && x59 && x60; POST_CHECK(x177);
    if (x39)
    {
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x32=0.5*x28; POST_CHECK_MUL(x32);
        x33=0.5*x29; POST_CHECK_MUL(x33);
        x34=-x32 - x33; POST_CHECK_ADD(x34);
        x35=1.4142135623730951454746218587388284504413604736328125*x28; POST_CHECK_MUL(x35);
        x36=x28*x29; POST_CHECK_MUL(x36);
    }
    else if (x49)
    {
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (x42) != 0.;*/
        x43=1.0/x42; POST_CHECK_POW(x43);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        /*@ assert (x44) != 0.;*/
        x45=1.0/x44; POST_CHECK_POW(x45);
        x46=0.5*un; POST_CHECK_MUL(x46);
        /*@ assert (x46) >= 0.;*/
        /*@ assert (x46) != 0.;*/
        x47=x42*x46; POST_CHECK_MUL(x47);
        x48=x42*x44; POST_CHECK_MUL(x48);
    }
    else if (x61)
    {
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); POST_CHECK_POW(x24);
        /*@ assert (x13 + x22) >= 0.;*/
        x26=sqrt(x13 + x22); POST_CHECK_POW(x26);
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        x51=x14*x15; POST_CHECK_MUL(x51);
        x52=x17*x18; POST_CHECK_MUL(x52);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x54=x53*(x51 + x52); POST_CHECK_MUL(x54);
        x55=x50*(x4 + x54); POST_CHECK_MUL(x55);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x57=x56*(x4 - x54); POST_CHECK_MUL(x57);
    }
    else if (x100)
    {
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); POST_CHECK_POW(x24);
        /*@ assert (x13 + x22) >= 0.;*/
        x26=sqrt(x13 + x22); POST_CHECK_POW(x26);
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x92=x8*rn; POST_CHECK_MUL(x92);
        x93=mu*rt1; POST_CHECK_MUL(x93);
        x94=x15*x93; POST_CHECK_MUL(x94);
        x95=mu*rt2; POST_CHECK_MUL(x95);
        x96=x18*x95; POST_CHECK_MUL(x96);
        x97=x53*(x94 + x96); POST_CHECK_MUL(x97);
        x98=x50*(x92 + x97); POST_CHECK_MUL(x98);
        x99=x56*(x92 - x97); POST_CHECK_MUL(x99);
    }
    else if (x59)
    {
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); POST_CHECK_POW(x24);
        x25=0.5*x24; POST_CHECK_MUL(x25);
        /*@ assert (x13 + x22) >= 0.;*/
        x26=sqrt(x13 + x22); POST_CHECK_POW(x26);
        x27=0.5*x26; POST_CHECK_MUL(x27);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x111=rt1 + x14; POST_CHECK_ADD(x111);
        x112=x15*x53; POST_CHECK_MUL(x112);
    }
    else if (x37)
    {
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); POST_CHECK_POW(x24);
        x25=0.5*x24; POST_CHECK_MUL(x25);
        /*@ assert (x13 + x22) >= 0.;*/
        x26=sqrt(x13 + x22); POST_CHECK_POW(x26);
        x27=0.5*x26; POST_CHECK_MUL(x27);
        x111=rt1 + x14; POST_CHECK_ADD(x111);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);
    }
    else if (x125)
    {
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); POST_CHECK_POW(x24);
        /*@ assert (x13 + x22) >= 0.;*/
        x26=sqrt(x13 + x22); POST_CHECK_POW(x26);
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        x51=x14*x15; POST_CHECK_MUL(x51);
        x52=x17*x18; POST_CHECK_MUL(x52);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x54=x53*(x51 + x52); POST_CHECK_MUL(x54);
        x55=x50*(x4 + x54); POST_CHECK_MUL(x55);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x57=x56*(x4 - x54); POST_CHECK_MUL(x57);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);
    }
    else if (x177)
    {
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); POST_CHECK_POW(x24);
        /*@ assert (x13 + x22) >= 0.;*/
        x26=sqrt(x13 + x22); POST_CHECK_POW(x26);
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x92=x8*rn; POST_CHECK_MUL(x92);
        x93=mu*rt1; POST_CHECK_MUL(x93);
        x94=x15*x93; POST_CHECK_MUL(x94);
        x95=mu*rt2; POST_CHECK_MUL(x95);
        x96=x18*x95; POST_CHECK_MUL(x96);
        x97=x53*(x94 + x96); POST_CHECK_MUL(x97);
        x98=x50*(x92 + x97); POST_CHECK_MUL(x98);
        x99=x56*(x92 - x97); POST_CHECK_MUL(x99);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);
    }
    /* Assignment result[0, 0]=-x25 - x27 + x4 + x5 */
    /*@ assert (x23) >= 0.;*/
    x24=sqrt(x23); POST_CHECK_POW(x24);
    x25=0.5*x24; POST_CHECK_MUL(x25);
    /*@ assert (x13 + x22) >= 0.;*/
    x26=sqrt(x13 + x22); POST_CHECK_POW(x26);
    x27=0.5*x26; POST_CHECK_MUL(x27);
    result[0] = -x25 - x27 + x4 + x5;


    /* Assignment result[1, 0]=Piecewise((x111 + x112*x25 - x112*x27, x59), (x111 + x114*x25 - x114*x27, x37)) */

    if (x59)
    {
        DEBUG_PRINT("Case (x59) is True.\n");
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x111=rt1 + x14; POST_CHECK_ADD(x111);
        x112=x15*x53; POST_CHECK_MUL(x112);

        /* Assignment result[1, 0]=x111 + x112*x25 - x112*x27 */
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x111=rt1 + x14; POST_CHECK_ADD(x111);
        x112=x15*x53; POST_CHECK_MUL(x112);
        result[1] = x111 + x112*x25 - x112*x27;

    }
    else if (x37)
    {
        DEBUG_PRINT("Case (x37) is True.\n");
        x111=rt1 + x14; POST_CHECK_ADD(x111);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);

        /* Assignment result[1, 0]=x111 + x114*x25 - x114*x27 */
        x111=rt1 + x14; POST_CHECK_ADD(x111);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);
        result[1] = x111 + x114*x25 - x114*x27;

    }


    /* Assignment result[2, 0]=Piecewise((x187 + x188*x25 - x188*x27, x59), (x187 + x189*x25 - x189*x27, x37)) */
    double x187;
    double x188;
    double x189;
    if (x59)
    {
        DEBUG_PRINT("Case (x59) is True.\n");
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x187=rt2 + x17; POST_CHECK_ADD(x187);
        x188=x18*x53; POST_CHECK_MUL(x188);

        /* Assignment result[2, 0]=x187 + x188*x25 - x188*x27 */
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x187=rt2 + x17; POST_CHECK_ADD(x187);
        x188=x18*x53; POST_CHECK_MUL(x188);
        result[2] = x187 + x188*x25 - x188*x27;

    }
    else if (x37)
    {
        DEBUG_PRINT("Case (x37) is True.\n");
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x187=rt2 + x17; POST_CHECK_ADD(x187);
        x189=random2*x113; POST_CHECK_MUL(x189);

        /* Assignment result[2, 0]=x187 + x189*x25 - x189*x27 */
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x187=rt2 + x17; POST_CHECK_ADD(x187);
        x189=random2*x113; POST_CHECK_MUL(x189);
        result[2] = x187 + x189*x25 - x189*x27;

    }


    /* Assignment result[0, 1]=Piecewise((x31*(-x35*mu + x34 + x36), x39), (x43*x45*(-x44*x46 - x47 + x48), x49), (1 - x55 - x57, x61)) */

    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x32=0.5*x28; POST_CHECK_MUL(x32);
        x33=0.5*x29; POST_CHECK_MUL(x33);
        x34=-x32 - x33; POST_CHECK_ADD(x34);
        x35=1.4142135623730951454746218587388284504413604736328125*x28; POST_CHECK_MUL(x35);
        x36=x28*x29; POST_CHECK_MUL(x36);

        /* Assignment result[0, 1]=x31*(-x35*mu + x34 + x36) */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x32=0.5*x28; POST_CHECK_MUL(x32);
        x33=0.5*x29; POST_CHECK_MUL(x33);
        x34=-x32 - x33; POST_CHECK_ADD(x34);
        x35=1.4142135623730951454746218587388284504413604736328125*x28; POST_CHECK_MUL(x35);
        x36=x28*x29; POST_CHECK_MUL(x36);
        result[3] = x31*(-x35*mu + x34 + x36);

    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (x42) != 0.;*/
        x43=1.0/x42; POST_CHECK_POW(x43);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        /*@ assert (x44) != 0.;*/
        x45=1.0/x44; POST_CHECK_POW(x45);
        x46=0.5*un; POST_CHECK_MUL(x46);
        /*@ assert (x46) >= 0.;*/
        /*@ assert (x46) != 0.;*/
        x47=x42*x46; POST_CHECK_MUL(x47);
        x48=x42*x44; POST_CHECK_MUL(x48);

        /* Assignment result[0, 1]=x43*x45*(-x44*x46 - x47 + x48) */
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (x42) != 0.;*/
        x43=1.0/x42; POST_CHECK_POW(x43);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        /*@ assert (x44) != 0.;*/
        x45=1.0/x44; POST_CHECK_POW(x45);
        x46=0.5*un; POST_CHECK_MUL(x46);
        /*@ assert (x46) >= 0.;*/
        /*@ assert (x46) != 0.;*/
        x47=x42*x46; POST_CHECK_MUL(x47);
        x48=x42*x44; POST_CHECK_MUL(x48);
        result[3] = x43*x45*(-x44*x46 - x47 + x48);

    }
    else if (x61)
    {
        DEBUG_PRINT("Case (x61) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        x51=x14*x15; POST_CHECK_MUL(x51);
        x52=x17*x18; POST_CHECK_MUL(x52);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x54=x53*(x51 + x52); POST_CHECK_MUL(x54);
        x55=x50*(x4 + x54); POST_CHECK_MUL(x55);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x57=x56*(x4 - x54); POST_CHECK_MUL(x57);

        /* Assignment result[0, 1]=1 - x55 - x57 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        x51=x14*x15; POST_CHECK_MUL(x51);
        x52=x17*x18; POST_CHECK_MUL(x52);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x54=x53*(x51 + x52); POST_CHECK_MUL(x54);
        x55=x50*(x4 + x54); POST_CHECK_MUL(x55);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x57=x56*(x4 - x54); POST_CHECK_MUL(x57);
        result[3] = 1 - x55 - x57;

    }


    /* Assignment result[1, 1]=Piecewise((x115, x39), (x117*(0.5*x118*x119*x44*x7*rt1*un - x118*x119*x47*x7*rt1 + 0.5*x118*x119*x120*x44*un - x118*x119*x120*x47), x49), (-x112*x55 + x112*x57 + x124*x25 - x124*x27, x61), (-x114*x55 + x114*x57, x125)) */
    double x65;
    double x90;
    double x101;
    double x102;
    double x115;
    double x116;
    double x117;
    double x118;
    double x119;
    double x120;
    double x121;
    double x122;
    double x123;
    double x124;
    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x101=0.353553390593273786368655464684707112610340118408203125*x28; POST_CHECK_MUL(x101);
        x102=0.353553390593273786368655464684707112610340118408203125*x29; POST_CHECK_MUL(x102);
        x115=x31*(-x101 + x102 - 1.0*x90); POST_CHECK_MUL(x115);

        /* Assignment result[1, 1]=x115 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x101=0.353553390593273786368655464684707112610340118408203125*x28; POST_CHECK_MUL(x101);
        x102=0.353553390593273786368655464684707112610340118408203125*x29; POST_CHECK_MUL(x102);
        x115=x31*(-x101 + x102 - 1.0*x90); POST_CHECK_MUL(x115);
        result[4] = x115;

    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (x42) != 0.;*/
        x43=1.0/x42; POST_CHECK_POW(x43);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        /*@ assert (x44) != 0.;*/
        x45=1.0/x44; POST_CHECK_POW(x45);
        x46=0.5*un; POST_CHECK_MUL(x46);
        /*@ assert (x46) >= 0.;*/
        /*@ assert (x46) != 0.;*/
        x47=x42*x46; POST_CHECK_MUL(x47);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x116=pow(x65, 3.0/2.0); POST_CHECK_POW(x116);
        /*@ assert (x116) != 0.;*/
        x117=x43*x45/x116; POST_CHECK_MUL(x117);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x119=pow(rn, 3); POST_CHECK_POW(x119);
        /*@ assert (x119) >= 0.;*/
        /*@ assert (x119) != 0.;*/
        x120=pow(rt1, 3); POST_CHECK_POW(x120);

        /* Assignment result[1, 1]=x117*(0.5*x118*x119*x44*x7*rt1*un - x118*x119*x47*x7*rt1 + 0.5*x118*x119*x120*x44*un - x118*x119*x120*x47) */
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (x42) != 0.;*/
        x43=1.0/x42; POST_CHECK_POW(x43);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        /*@ assert (x44) != 0.;*/
        x45=1.0/x44; POST_CHECK_POW(x45);
        x46=0.5*un; POST_CHECK_MUL(x46);
        /*@ assert (x46) >= 0.;*/
        /*@ assert (x46) != 0.;*/
        x47=x42*x46; POST_CHECK_MUL(x47);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x116=pow(x65, 3.0/2.0); POST_CHECK_POW(x116);
        /*@ assert (x116) != 0.;*/
        x117=x43*x45/x116; POST_CHECK_MUL(x117);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x119=pow(rn, 3); POST_CHECK_POW(x119);
        /*@ assert (x119) >= 0.;*/
        /*@ assert (x119) != 0.;*/
        x120=pow(rt1, 3); POST_CHECK_POW(x120);
        result[4] = x117*(0.5*x118*x119*x44*x7*rt1*un - x118*x119*x47*x7*rt1 + 0.5*x118*x119*x120*x44*un - x118*x119*x120*x47);

    }
    else if (x61)
    {
        DEBUG_PRINT("Case (x61) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        x51=x14*x15; POST_CHECK_MUL(x51);
        x52=x17*x18; POST_CHECK_MUL(x52);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x54=x53*(x51 + x52); POST_CHECK_MUL(x54);
        x55=x50*(x4 + x54); POST_CHECK_MUL(x55);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x57=x56*(x4 - x54); POST_CHECK_MUL(x57);
        x112=x15*x53; POST_CHECK_MUL(x112);
        x121=-x51 - x52; POST_CHECK_ADD(x121);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x123=x122*x15; POST_CHECK_MUL(x123);
        x124=x121*x123 + x14*x53; POST_CHECK_ADD(x124);

        /* Assignment result[1, 1]=-x112*x55 + x112*x57 + x124*x25 - x124*x27 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        x51=x14*x15; POST_CHECK_MUL(x51);
        x52=x17*x18; POST_CHECK_MUL(x52);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x54=x53*(x51 + x52); POST_CHECK_MUL(x54);
        x55=x50*(x4 + x54); POST_CHECK_MUL(x55);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x57=x56*(x4 - x54); POST_CHECK_MUL(x57);
        x112=x15*x53; POST_CHECK_MUL(x112);
        x121=-x51 - x52; POST_CHECK_ADD(x121);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x123=x122*x15; POST_CHECK_MUL(x123);
        x124=x121*x123 + x14*x53; POST_CHECK_ADD(x124);
        result[4] = -x112*x55 + x112*x57 + x124*x25 - x124*x27;

    }
    else if (x125)
    {
        DEBUG_PRINT("Case (x125) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        x51=x14*x15; POST_CHECK_MUL(x51);
        x52=x17*x18; POST_CHECK_MUL(x52);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x54=x53*(x51 + x52); POST_CHECK_MUL(x54);
        x55=x50*(x4 + x54); POST_CHECK_MUL(x55);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x57=x56*(x4 - x54); POST_CHECK_MUL(x57);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);

        /* Assignment result[1, 1]=-x114*x55 + x114*x57 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        x51=x14*x15; POST_CHECK_MUL(x51);
        x52=x17*x18; POST_CHECK_MUL(x52);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x54=x53*(x51 + x52); POST_CHECK_MUL(x54);
        x55=x50*(x4 + x54); POST_CHECK_MUL(x55);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x57=x56*(x4 - x54); POST_CHECK_MUL(x57);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);
        result[4] = -x114*x55 + x114*x57;

    }


    /* Assignment result[2, 1]=Piecewise((x115, x39), (x117*(0.5*x118*x119*x44*x6*rt2*un - x118*x119*x47*x6*rt2 + 0.5*x118*x119*x169*x44*un - x118*x119*x169*x47), x49), (-x188*x55 + x188*x57 + x191*x25 - x191*x27, x61), (-x189*x55 + x189*x57, x125)) */
    double x169;
    double x190;
    double x191;
    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x101=0.353553390593273786368655464684707112610340118408203125*x28; POST_CHECK_MUL(x101);
        x102=0.353553390593273786368655464684707112610340118408203125*x29; POST_CHECK_MUL(x102);
        x115=x31*(-x101 + x102 - 1.0*x90); POST_CHECK_MUL(x115);

        /* Assignment result[2, 1]=x115 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x101=0.353553390593273786368655464684707112610340118408203125*x28; POST_CHECK_MUL(x101);
        x102=0.353553390593273786368655464684707112610340118408203125*x29; POST_CHECK_MUL(x102);
        x115=x31*(-x101 + x102 - 1.0*x90); POST_CHECK_MUL(x115);
        result[5] = x115;

    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (x42) != 0.;*/
        x43=1.0/x42; POST_CHECK_POW(x43);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        /*@ assert (x44) != 0.;*/
        x45=1.0/x44; POST_CHECK_POW(x45);
        x46=0.5*un; POST_CHECK_MUL(x46);
        /*@ assert (x46) >= 0.;*/
        /*@ assert (x46) != 0.;*/
        x47=x42*x46; POST_CHECK_MUL(x47);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x116=pow(x65, 3.0/2.0); POST_CHECK_POW(x116);
        /*@ assert (x116) != 0.;*/
        x117=x43*x45/x116; POST_CHECK_MUL(x117);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x119=pow(rn, 3); POST_CHECK_POW(x119);
        /*@ assert (x119) >= 0.;*/
        /*@ assert (x119) != 0.;*/
        x169=pow(rt2, 3); POST_CHECK_POW(x169);

        /* Assignment result[2, 1]=x117*(0.5*x118*x119*x44*x6*rt2*un - x118*x119*x47*x6*rt2 + 0.5*x118*x119*x169*x44*un - x118*x119*x169*x47) */
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (x42) != 0.;*/
        x43=1.0/x42; POST_CHECK_POW(x43);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        /*@ assert (x44) != 0.;*/
        x45=1.0/x44; POST_CHECK_POW(x45);
        x46=0.5*un; POST_CHECK_MUL(x46);
        /*@ assert (x46) >= 0.;*/
        /*@ assert (x46) != 0.;*/
        x47=x42*x46; POST_CHECK_MUL(x47);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x116=pow(x65, 3.0/2.0); POST_CHECK_POW(x116);
        /*@ assert (x116) != 0.;*/
        x117=x43*x45/x116; POST_CHECK_MUL(x117);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x119=pow(rn, 3); POST_CHECK_POW(x119);
        /*@ assert (x119) >= 0.;*/
        /*@ assert (x119) != 0.;*/
        x169=pow(rt2, 3); POST_CHECK_POW(x169);
        result[5] = x117*(0.5*x118*x119*x44*x6*rt2*un - x118*x119*x47*x6*rt2 + 0.5*x118*x119*x169*x44*un - x118*x119*x169*x47);

    }
    else if (x61)
    {
        DEBUG_PRINT("Case (x61) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        x51=x14*x15; POST_CHECK_MUL(x51);
        x52=x17*x18; POST_CHECK_MUL(x52);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x54=x53*(x51 + x52); POST_CHECK_MUL(x54);
        x55=x50*(x4 + x54); POST_CHECK_MUL(x55);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x57=x56*(x4 - x54); POST_CHECK_MUL(x57);
        x121=-x51 - x52; POST_CHECK_ADD(x121);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x188=x18*x53; POST_CHECK_MUL(x188);
        x190=x122*x18; POST_CHECK_MUL(x190);
        x191=x121*x190 + x17*x53; POST_CHECK_ADD(x191);

        /* Assignment result[2, 1]=-x188*x55 + x188*x57 + x191*x25 - x191*x27 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        x51=x14*x15; POST_CHECK_MUL(x51);
        x52=x17*x18; POST_CHECK_MUL(x52);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x54=x53*(x51 + x52); POST_CHECK_MUL(x54);
        x55=x50*(x4 + x54); POST_CHECK_MUL(x55);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x57=x56*(x4 - x54); POST_CHECK_MUL(x57);
        x121=-x51 - x52; POST_CHECK_ADD(x121);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x188=x18*x53; POST_CHECK_MUL(x188);
        x190=x122*x18; POST_CHECK_MUL(x190);
        x191=x121*x190 + x17*x53; POST_CHECK_ADD(x191);
        result[5] = -x188*x55 + x188*x57 + x191*x25 - x191*x27;

    }
    else if (x125)
    {
        DEBUG_PRINT("Case (x125) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        x51=x14*x15; POST_CHECK_MUL(x51);
        x52=x17*x18; POST_CHECK_MUL(x52);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x54=x53*(x51 + x52); POST_CHECK_MUL(x54);
        x55=x50*(x4 + x54); POST_CHECK_MUL(x55);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x57=x56*(x4 - x54); POST_CHECK_MUL(x57);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x189=random2*x113; POST_CHECK_MUL(x189);

        /* Assignment result[2, 1]=-x189*x55 + x189*x57 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        x51=x14*x15; POST_CHECK_MUL(x51);
        x52=x17*x18; POST_CHECK_MUL(x52);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        x54=x53*(x51 + x52); POST_CHECK_MUL(x54);
        x55=x50*(x4 + x54); POST_CHECK_MUL(x55);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x57=x56*(x4 - x54); POST_CHECK_MUL(x57);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x189=random2*x113; POST_CHECK_MUL(x189);
        result[5] = -x189*x55 + x189*x57;

    }


    /* Assignment result[0, 2]=Piecewise((x64, x39), (x67*(0.5*x44*x8*rt1*rn*un - x47*x8*rt1*rn + x69), x49), (x71 - x80 - x81, x61)) */
    double x62;
    double x63;
    double x64;
    double x66;
    double x67;
    double x68;
    double x69;
    double x70;
    double x71;
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
    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        x62=0.70710678118654757273731092936941422522068023681640625*mu; POST_CHECK_MUL(x62);
        x63=x29*x62; POST_CHECK_MUL(x63);
        x64=x30*(-x62 + x63 - 2.0*x8); POST_CHECK_MUL(x64);

        /* Assignment result[0, 2]=x64 */
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        x62=0.70710678118654757273731092936941422522068023681640625*mu; POST_CHECK_MUL(x62);
        x63=x29*x62; POST_CHECK_MUL(x63);
        x64=x30*(-x62 + x63 - 2.0*x8); POST_CHECK_MUL(x64);
        result[6] = x64;

    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (x42) != 0.;*/
        x43=1.0/x42; POST_CHECK_POW(x43);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        /*@ assert (x44) != 0.;*/
        x45=1.0/x44; POST_CHECK_POW(x45);
        x46=0.5*un; POST_CHECK_MUL(x46);
        /*@ assert (x46) >= 0.;*/
        /*@ assert (x46) != 0.;*/
        x47=x42*x46; POST_CHECK_MUL(x47);
        x48=x42*x44; POST_CHECK_MUL(x48);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x66=sqrt(x65); POST_CHECK_POW(x66);
        /*@ assert (x66) != 0.;*/
        x67=x43*x45/x66; POST_CHECK_MUL(x67);
        x68=0.353553390593273786368655464684707112610340118408203125*x66*mu*un; POST_CHECK_MUL(x68);
        x69=0.70710678118654757273731092936941422522068023681640625*x48*x66*mu - x42*x68 - x44*x68; POST_CHECK_ADD(x69);

        /* Assignment result[0, 2]=x67*(0.5*x44*x8*rt1*rn*un - x47*x8*rt1*rn + x69) */
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (x42) != 0.;*/
        x43=1.0/x42; POST_CHECK_POW(x43);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        /*@ assert (x44) != 0.;*/
        x45=1.0/x44; POST_CHECK_POW(x45);
        x46=0.5*un; POST_CHECK_MUL(x46);
        /*@ assert (x46) >= 0.;*/
        /*@ assert (x46) != 0.;*/
        x47=x42*x46; POST_CHECK_MUL(x47);
        x48=x42*x44; POST_CHECK_MUL(x48);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x66=sqrt(x65); POST_CHECK_POW(x66);
        /*@ assert (x66) != 0.;*/
        x67=x43*x45/x66; POST_CHECK_MUL(x67);
        x68=0.353553390593273786368655464684707112610340118408203125*x66*mu*un; POST_CHECK_MUL(x68);
        x69=0.70710678118654757273731092936941422522068023681640625*x48*x66*mu - x42*x68 - x44*x68; POST_CHECK_ADD(x69);
        result[6] = x67*(0.5*x44*x8*rt1*rn*un - x47*x8*rt1*rn + x69);

    }
    else if (x61)
    {
        DEBUG_PRINT("Case (x61) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x71=x14*x70; POST_CHECK_MUL(x71);
        x72=x8*ut1 + x4*x71; POST_CHECK_ADD(x72);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x74=x18*x73; POST_CHECK_MUL(x74);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x77=x11*x70; POST_CHECK_MUL(x77);
        x78=(1.0/2.0)*x15*(x76 + 2*x77); POST_CHECK_MUL(x78);
        x79=x53*(x74 + x78); POST_CHECK_MUL(x79);
        x80=x50*(x72 + x79); POST_CHECK_MUL(x80);
        x81=x56*(x72 - x79); POST_CHECK_MUL(x81);

        /* Assignment result[0, 2]=x71 - x80 - x81 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x71=x14*x70; POST_CHECK_MUL(x71);
        x72=x8*ut1 + x4*x71; POST_CHECK_ADD(x72);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x74=x18*x73; POST_CHECK_MUL(x74);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x77=x11*x70; POST_CHECK_MUL(x77);
        x78=(1.0/2.0)*x15*(x76 + 2*x77); POST_CHECK_MUL(x78);
        x79=x53*(x74 + x78); POST_CHECK_MUL(x79);
        x80=x50*(x72 + x79); POST_CHECK_MUL(x80);
        x81=x56*(x72 - x79); POST_CHECK_MUL(x81);
        result[6] = x71 - x80 - x81;

    }


    /* Assignment result[1, 2]=Piecewise((x138, x39), (x142*(-2.0*x118*x147*x42*x9*un + 2.0*x118*x147*x44*x9*un + 4.0*x118*x152*x44*x6*x9*un + 2.0*x118*x157*x44*x7*x9*un - 2.0*x143*x144*x152*x42*un + 2.0*x143*x144*x152*x44*un - 2.0*x118*x145*x152*x42*x9 + 2.0*x118*x145*x152*x44*x9 - x136*x158 + x146 + x150 + x151 + x153 + x154 - x155*x156 + x159 + x160 - 2.0*x161 - 2.0*x162 - 4.0*x163 - 4.0*x164), x49), (mu - x112*x80 + x112*x81 + x166*x25 - x166*x27, x61), (mu - x114*x80 + x114*x81, x125)) */
    double x91;
    double x126;
    double x127;
    double x128;
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
    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x128=4.24264068711928477029005080112256109714508056640625*x28; POST_CHECK_MUL(x128);
        x129=4.24264068711928477029005080112256109714508056640625*x29; POST_CHECK_MUL(x129);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x131=2.0*x130; POST_CHECK_MUL(x131);
        x132=x28*x8; POST_CHECK_MUL(x132);
        x133=x118*x28; POST_CHECK_MUL(x133);
        x134=x129*x8; POST_CHECK_MUL(x134);
        x135=-x134; POST_CHECK_MUL(x135);
        x136=2.0*x118; POST_CHECK_MUL(x136);
        x137=x136*x29; POST_CHECK_MUL(x137);
        x138=x127*(-x128 + x129 + x131 - 57.9827560572968963015227927826344966888427734375*x132 + 11.3137084989847611637969748699106276035308837890625*x132*x29 - 34.0*x133 + x135 + x137 - 26.0*x90 + 16.0*x91); POST_CHECK_MUL(x138);

        /* Assignment result[1, 2]=x138 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x128=4.24264068711928477029005080112256109714508056640625*x28; POST_CHECK_MUL(x128);
        x129=4.24264068711928477029005080112256109714508056640625*x29; POST_CHECK_MUL(x129);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x131=2.0*x130; POST_CHECK_MUL(x131);
        x132=x28*x8; POST_CHECK_MUL(x132);
        x133=x118*x28; POST_CHECK_MUL(x133);
        x134=x129*x8; POST_CHECK_MUL(x134);
        x135=-x134; POST_CHECK_MUL(x135);
        x136=2.0*x118; POST_CHECK_MUL(x136);
        x137=x136*x29; POST_CHECK_MUL(x137);
        x138=x127*(-x128 + x129 + x131 - 57.9827560572968963015227927826344966888427734375*x132 + 11.3137084989847611637969748699106276035308837890625*x132*x29 - 34.0*x133 + x135 + x137 - 26.0*x90 + 16.0*x91); POST_CHECK_MUL(x138);
        result[7] = x138;

    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x116=pow(x65, 3.0/2.0); POST_CHECK_POW(x116);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x119=pow(rn, 3); POST_CHECK_POW(x119);
        /*@ assert (x119) >= 0.;*/
        /*@ assert (x119) != 0.;*/
        x120=pow(rt1, 3); POST_CHECK_POW(x120);
        x136=2.0*x118; POST_CHECK_MUL(x136);
        x139=4.0*x116*x42*x44; POST_CHECK_MUL(x139);
        x140=x139*x6; POST_CHECK_MUL(x140);
        x141=x139*x7; POST_CHECK_MUL(x141);
        /*@ assert (x140 + x141) != 0.;*/
        x142=1.0/(x140 + x141); POST_CHECK_POW(x142);
        x143=pow(mu, 5); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x144=pow(rn, 4); POST_CHECK_POW(x144);
        /*@ assert (x144) >= 0.;*/
        /*@ assert (x144) != 0.;*/
        x145=pow(un, 3); POST_CHECK_POW(x145);
        /*@ assert (x145) >= 0.;*/
        /*@ assert (x145) != 0.;*/
        x146=x140*mu + x141*mu - 2.0*x143*x144*x42*x6*x7*un + 2.0*x143*x144*x44*x6*x7*un - 2.0*x118*x145*x42*x6*x7*x9 + 2.0*x118*x145*x44*x6*x7*x9; POST_CHECK_ADD(x146);
        x147=pow(rt2, 6); POST_CHECK_POW(x147);
        /*@ assert (x147) >= 0.;*/
        x148=pow(mu, 4); POST_CHECK_POW(x148);
        /*@ assert (x148) >= 0.;*/
        /*@ assert (x148) != 0.;*/
        x149=pow(rt1, 5); POST_CHECK_POW(x149);
        x150=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x149*x42*un; POST_CHECK_MUL(x150);
        x151=1.4142135623730951454746218587388284504413604736328125*x119*x148*x149*x44*un; POST_CHECK_MUL(x151);
        x152=pow(rt2, 4); POST_CHECK_POW(x152);
        /*@ assert (x152) >= 0.;*/
        x153=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x152*x42*rt1*un; POST_CHECK_MUL(x153);
        x154=1.4142135623730951454746218587388284504413604736328125*x119*x148*x152*x44*rt1*un; POST_CHECK_MUL(x154);
        x155=4.0*x118; POST_CHECK_MUL(x155);
        x156=x152*x42*x6*x9*un; POST_CHECK_MUL(x156);
        x157=pow(rt1, 4); POST_CHECK_POW(x157);
        /*@ assert (x157) >= 0.;*/
        x158=x157*x42*x7*x9*un; POST_CHECK_MUL(x158);
        x159=-2.828427124746190290949243717477656900882720947265625*x119*x120*x148*x42*x7*un; POST_CHECK_MUL(x159);
        x160=2.828427124746190290949243717477656900882720947265625*x119*x120*x148*x44*x7*un; POST_CHECK_MUL(x160);
        x161=x116*x42*x6*mu*un; POST_CHECK_MUL(x161);
        x162=x116*x44*x6*mu*un; POST_CHECK_MUL(x162);
        x163=x116*x42*x7*mu*un; POST_CHECK_MUL(x163);
        x164=x116*x44*x7*mu*un; POST_CHECK_MUL(x164);

        /* Assignment result[1, 2]=x142*(-2.0*x118*x147*x42*x9*un + 2.0*x118*x147*x44*x9*un + 4.0*x118*x152*x44*x6*x9*un + 2.0*x118*x157*x44*x7*x9*un - 2.0*x143*x144*x152*x42*un + 2.0*x143*x144*x152*x44*un - 2.0*x118*x145*x152*x42*x9 + 2.0*x118*x145*x152*x44*x9 - x136*x158 + x146 + x150 + x151 + x153 + x154 - x155*x156 + x159 + x160 - 2.0*x161 - 2.0*x162 - 4.0*x163 - 4.0*x164) */
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x116=pow(x65, 3.0/2.0); POST_CHECK_POW(x116);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x119=pow(rn, 3); POST_CHECK_POW(x119);
        /*@ assert (x119) >= 0.;*/
        /*@ assert (x119) != 0.;*/
        x120=pow(rt1, 3); POST_CHECK_POW(x120);
        x136=2.0*x118; POST_CHECK_MUL(x136);
        x139=4.0*x116*x42*x44; POST_CHECK_MUL(x139);
        x140=x139*x6; POST_CHECK_MUL(x140);
        x141=x139*x7; POST_CHECK_MUL(x141);
        /*@ assert (x140 + x141) != 0.;*/
        x142=1.0/(x140 + x141); POST_CHECK_POW(x142);
        x143=pow(mu, 5); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x144=pow(rn, 4); POST_CHECK_POW(x144);
        /*@ assert (x144) >= 0.;*/
        /*@ assert (x144) != 0.;*/
        x145=pow(un, 3); POST_CHECK_POW(x145);
        /*@ assert (x145) >= 0.;*/
        /*@ assert (x145) != 0.;*/
        x146=x140*mu + x141*mu - 2.0*x143*x144*x42*x6*x7*un + 2.0*x143*x144*x44*x6*x7*un - 2.0*x118*x145*x42*x6*x7*x9 + 2.0*x118*x145*x44*x6*x7*x9; POST_CHECK_ADD(x146);
        x147=pow(rt2, 6); POST_CHECK_POW(x147);
        /*@ assert (x147) >= 0.;*/
        x148=pow(mu, 4); POST_CHECK_POW(x148);
        /*@ assert (x148) >= 0.;*/
        /*@ assert (x148) != 0.;*/
        x149=pow(rt1, 5); POST_CHECK_POW(x149);
        x150=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x149*x42*un; POST_CHECK_MUL(x150);
        x151=1.4142135623730951454746218587388284504413604736328125*x119*x148*x149*x44*un; POST_CHECK_MUL(x151);
        x152=pow(rt2, 4); POST_CHECK_POW(x152);
        /*@ assert (x152) >= 0.;*/
        x153=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x152*x42*rt1*un; POST_CHECK_MUL(x153);
        x154=1.4142135623730951454746218587388284504413604736328125*x119*x148*x152*x44*rt1*un; POST_CHECK_MUL(x154);
        x155=4.0*x118; POST_CHECK_MUL(x155);
        x156=x152*x42*x6*x9*un; POST_CHECK_MUL(x156);
        x157=pow(rt1, 4); POST_CHECK_POW(x157);
        /*@ assert (x157) >= 0.;*/
        x158=x157*x42*x7*x9*un; POST_CHECK_MUL(x158);
        x159=-2.828427124746190290949243717477656900882720947265625*x119*x120*x148*x42*x7*un; POST_CHECK_MUL(x159);
        x160=2.828427124746190290949243717477656900882720947265625*x119*x120*x148*x44*x7*un; POST_CHECK_MUL(x160);
        x161=x116*x42*x6*mu*un; POST_CHECK_MUL(x161);
        x162=x116*x44*x6*mu*un; POST_CHECK_MUL(x162);
        x163=x116*x42*x7*mu*un; POST_CHECK_MUL(x163);
        x164=x116*x44*x7*mu*un; POST_CHECK_MUL(x164);
        result[7] = x142*(-2.0*x118*x147*x42*x9*un + 2.0*x118*x147*x44*x9*un + 4.0*x118*x152*x44*x6*x9*un + 2.0*x118*x157*x44*x7*x9*un - 2.0*x143*x144*x152*x42*un + 2.0*x143*x144*x152*x44*un - 2.0*x118*x145*x152*x42*x9 + 2.0*x118*x145*x152*x44*x9 - x136*x158 + x146 + x150 + x151 + x153 + x154 - x155*x156 + x159 + x160 - 2.0*x161 - 2.0*x162 - 4.0*x163 - 4.0*x164);

    }
    else if (x61)
    {
        DEBUG_PRINT("Case (x61) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x71=x14*x70; POST_CHECK_MUL(x71);
        x72=x8*ut1 + x4*x71; POST_CHECK_ADD(x72);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x74=x18*x73; POST_CHECK_MUL(x74);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x77=x11*x70; POST_CHECK_MUL(x77);
        x78=(1.0/2.0)*x15*(x76 + 2*x77); POST_CHECK_MUL(x78);
        x79=x53*(x74 + x78); POST_CHECK_MUL(x79);
        x80=x50*(x72 + x79); POST_CHECK_MUL(x80);
        x81=x56*(x72 - x79); POST_CHECK_MUL(x81);
        x112=x15*x53; POST_CHECK_MUL(x112);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x123=x122*x15; POST_CHECK_MUL(x123);
        x165=-x74 - x78; POST_CHECK_ADD(x165);
        x166=x53*(x75 + x77) + x123*x165; POST_CHECK_ADD(x166);

        /* Assignment result[1, 2]=mu - x112*x80 + x112*x81 + x166*x25 - x166*x27 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x71=x14*x70; POST_CHECK_MUL(x71);
        x72=x8*ut1 + x4*x71; POST_CHECK_ADD(x72);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x74=x18*x73; POST_CHECK_MUL(x74);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x77=x11*x70; POST_CHECK_MUL(x77);
        x78=(1.0/2.0)*x15*(x76 + 2*x77); POST_CHECK_MUL(x78);
        x79=x53*(x74 + x78); POST_CHECK_MUL(x79);
        x80=x50*(x72 + x79); POST_CHECK_MUL(x80);
        x81=x56*(x72 - x79); POST_CHECK_MUL(x81);
        x112=x15*x53; POST_CHECK_MUL(x112);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x123=x122*x15; POST_CHECK_MUL(x123);
        x165=-x74 - x78; POST_CHECK_ADD(x165);
        x166=x53*(x75 + x77) + x123*x165; POST_CHECK_ADD(x166);
        result[7] = mu - x112*x80 + x112*x81 + x166*x25 - x166*x27;

    }
    else if (x125)
    {
        DEBUG_PRINT("Case (x125) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x71=x14*x70; POST_CHECK_MUL(x71);
        x72=x8*ut1 + x4*x71; POST_CHECK_ADD(x72);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x74=x18*x73; POST_CHECK_MUL(x74);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x77=x11*x70; POST_CHECK_MUL(x77);
        x78=(1.0/2.0)*x15*(x76 + 2*x77); POST_CHECK_MUL(x78);
        x79=x53*(x74 + x78); POST_CHECK_MUL(x79);
        x80=x50*(x72 + x79); POST_CHECK_MUL(x80);
        x81=x56*(x72 - x79); POST_CHECK_MUL(x81);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);

        /* Assignment result[1, 2]=mu - x114*x80 + x114*x81 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x71=x14*x70; POST_CHECK_MUL(x71);
        x72=x8*ut1 + x4*x71; POST_CHECK_ADD(x72);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x74=x18*x73; POST_CHECK_MUL(x74);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x77=x11*x70; POST_CHECK_MUL(x77);
        x78=(1.0/2.0)*x15*(x76 + 2*x77); POST_CHECK_MUL(x78);
        x79=x53*(x74 + x78); POST_CHECK_MUL(x79);
        x80=x50*(x72 + x79); POST_CHECK_MUL(x80);
        x81=x56*(x72 - x79); POST_CHECK_MUL(x81);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);
        result[7] = mu - x114*x80 + x114*x81;

    }


    /* Assignment result[2, 2]=Piecewise((x167, x39), (x142*(x170 + x192 + x193 + x194 + x195 + x196 + x197), x49), (-x188*x80 + x188*x81 + x198*x25 - x198*x27, x61), (-x189*x80 + x189*x81, x125)) */
    double x167;
    double x168;
    double x170;
    double x171;
    double x192;
    double x193;
    double x194;
    double x195;
    double x196;
    double x197;
    double x198;
    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x35=1.4142135623730951454746218587388284504413604736328125*x28; POST_CHECK_MUL(x35);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x128=4.24264068711928477029005080112256109714508056640625*x28; POST_CHECK_MUL(x128);
        x129=4.24264068711928477029005080112256109714508056640625*x29; POST_CHECK_MUL(x129);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x131=2.0*x130; POST_CHECK_MUL(x131);
        x133=x118*x28; POST_CHECK_MUL(x133);
        x134=x129*x8; POST_CHECK_MUL(x134);
        x136=2.0*x118; POST_CHECK_MUL(x136);
        x137=x136*x29; POST_CHECK_MUL(x137);
        x167=x127*(x128 - x129 - x131 + 2.0*x133 + x134 - x137 + x35*x8 + 10.0*x90); POST_CHECK_MUL(x167);

        /* Assignment result[2, 2]=x167 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x35=1.4142135623730951454746218587388284504413604736328125*x28; POST_CHECK_MUL(x35);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x128=4.24264068711928477029005080112256109714508056640625*x28; POST_CHECK_MUL(x128);
        x129=4.24264068711928477029005080112256109714508056640625*x29; POST_CHECK_MUL(x129);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x131=2.0*x130; POST_CHECK_MUL(x131);
        x133=x118*x28; POST_CHECK_MUL(x133);
        x134=x129*x8; POST_CHECK_MUL(x134);
        x136=2.0*x118; POST_CHECK_MUL(x136);
        x137=x136*x29; POST_CHECK_MUL(x137);
        x167=x127*(x128 - x129 - x131 + 2.0*x133 + x134 - x137 + x35*x8 + 10.0*x90); POST_CHECK_MUL(x167);
        result[8] = x167;

    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x116=pow(x65, 3.0/2.0); POST_CHECK_POW(x116);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x119=pow(rn, 3); POST_CHECK_POW(x119);
        /*@ assert (x119) >= 0.;*/
        /*@ assert (x119) != 0.;*/
        x120=pow(rt1, 3); POST_CHECK_POW(x120);
        x139=4.0*x116*x42*x44; POST_CHECK_MUL(x139);
        x140=x139*x6; POST_CHECK_MUL(x140);
        x141=x139*x7; POST_CHECK_MUL(x141);
        /*@ assert (x140 + x141) != 0.;*/
        x142=1.0/(x140 + x141); POST_CHECK_POW(x142);
        x143=pow(mu, 5); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x144=pow(rn, 4); POST_CHECK_POW(x144);
        /*@ assert (x144) >= 0.;*/
        /*@ assert (x144) != 0.;*/
        x145=pow(un, 3); POST_CHECK_POW(x145);
        /*@ assert (x145) >= 0.;*/
        /*@ assert (x145) != 0.;*/
        x148=pow(mu, 4); POST_CHECK_POW(x148);
        /*@ assert (x148) >= 0.;*/
        /*@ assert (x148) != 0.;*/
        x149=pow(rt1, 5); POST_CHECK_POW(x149);
        x157=pow(rt1, 4); POST_CHECK_POW(x157);
        /*@ assert (x157) >= 0.;*/
        x168=pow(rt2, 5); POST_CHECK_POW(x168);
        x169=pow(rt2, 3); POST_CHECK_POW(x169);
        x170=2.0*x116*x42*rt1*rt2*mu*un + 2.0*x116*x44*rt1*rt2*mu*un + 2.0*x118*x168*x42*x9*rt1*un - 2.0*x118*x168*x44*x9*rt1*un + 2.0*x143*x144*x169*x42*rt1*un - 2.0*x143*x144*x169*x44*rt1*un + 2.0*x118*x145*x169*x42*x9*rt1 - 2.0*x118*x145*x169*x44*x9*rt1 + 2.0*x118*x149*x42*x9*rt2*un - 2.0*x118*x149*x44*x9*rt2*un + 2.0*x120*x143*x144*x42*rt2*un - 2.0*x120*x143*x144*x44*rt2*un + 2.0*x118*x120*x145*x42*x9*rt2 - 2.0*x118*x120*x145*x44*x9*rt2 + 4.0*x118*x120*x169*x42*x9*un - 4.0*x118*x120*x169*x44*x9*un; POST_CHECK_ADD(x170);
        x192=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x168*x42*un; POST_CHECK_MUL(x192);
        x193=1.4142135623730951454746218587388284504413604736328125*x119*x148*x168*x44*un; POST_CHECK_MUL(x193);
        x194=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x157*x42*rt2*un; POST_CHECK_MUL(x194);
        x195=1.4142135623730951454746218587388284504413604736328125*x119*x148*x157*x44*rt2*un; POST_CHECK_MUL(x195);
        x196=-2.828427124746190290949243717477656900882720947265625*x119*x148*x169*x42*x6*un; POST_CHECK_MUL(x196);
        x197=2.828427124746190290949243717477656900882720947265625*x119*x148*x169*x44*x6*un; POST_CHECK_MUL(x197);

        /* Assignment result[2, 2]=x142*(x170 + x192 + x193 + x194 + x195 + x196 + x197) */
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x116=pow(x65, 3.0/2.0); POST_CHECK_POW(x116);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x119=pow(rn, 3); POST_CHECK_POW(x119);
        /*@ assert (x119) >= 0.;*/
        /*@ assert (x119) != 0.;*/
        x120=pow(rt1, 3); POST_CHECK_POW(x120);
        x139=4.0*x116*x42*x44; POST_CHECK_MUL(x139);
        x140=x139*x6; POST_CHECK_MUL(x140);
        x141=x139*x7; POST_CHECK_MUL(x141);
        /*@ assert (x140 + x141) != 0.;*/
        x142=1.0/(x140 + x141); POST_CHECK_POW(x142);
        x143=pow(mu, 5); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x144=pow(rn, 4); POST_CHECK_POW(x144);
        /*@ assert (x144) >= 0.;*/
        /*@ assert (x144) != 0.;*/
        x145=pow(un, 3); POST_CHECK_POW(x145);
        /*@ assert (x145) >= 0.;*/
        /*@ assert (x145) != 0.;*/
        x148=pow(mu, 4); POST_CHECK_POW(x148);
        /*@ assert (x148) >= 0.;*/
        /*@ assert (x148) != 0.;*/
        x149=pow(rt1, 5); POST_CHECK_POW(x149);
        x157=pow(rt1, 4); POST_CHECK_POW(x157);
        /*@ assert (x157) >= 0.;*/
        x168=pow(rt2, 5); POST_CHECK_POW(x168);
        x169=pow(rt2, 3); POST_CHECK_POW(x169);
        x170=2.0*x116*x42*rt1*rt2*mu*un + 2.0*x116*x44*rt1*rt2*mu*un + 2.0*x118*x168*x42*x9*rt1*un - 2.0*x118*x168*x44*x9*rt1*un + 2.0*x143*x144*x169*x42*rt1*un - 2.0*x143*x144*x169*x44*rt1*un + 2.0*x118*x145*x169*x42*x9*rt1 - 2.0*x118*x145*x169*x44*x9*rt1 + 2.0*x118*x149*x42*x9*rt2*un - 2.0*x118*x149*x44*x9*rt2*un + 2.0*x120*x143*x144*x42*rt2*un - 2.0*x120*x143*x144*x44*rt2*un + 2.0*x118*x120*x145*x42*x9*rt2 - 2.0*x118*x120*x145*x44*x9*rt2 + 4.0*x118*x120*x169*x42*x9*un - 4.0*x118*x120*x169*x44*x9*un; POST_CHECK_ADD(x170);
        x192=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x168*x42*un; POST_CHECK_MUL(x192);
        x193=1.4142135623730951454746218587388284504413604736328125*x119*x148*x168*x44*un; POST_CHECK_MUL(x193);
        x194=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x157*x42*rt2*un; POST_CHECK_MUL(x194);
        x195=1.4142135623730951454746218587388284504413604736328125*x119*x148*x157*x44*rt2*un; POST_CHECK_MUL(x195);
        x196=-2.828427124746190290949243717477656900882720947265625*x119*x148*x169*x42*x6*un; POST_CHECK_MUL(x196);
        x197=2.828427124746190290949243717477656900882720947265625*x119*x148*x169*x44*x6*un; POST_CHECK_MUL(x197);
        result[8] = x142*(x170 + x192 + x193 + x194 + x195 + x196 + x197);

    }
    else if (x61)
    {
        DEBUG_PRINT("Case (x61) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x71=x14*x70; POST_CHECK_MUL(x71);
        x72=x8*ut1 + x4*x71; POST_CHECK_ADD(x72);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x74=x18*x73; POST_CHECK_MUL(x74);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x77=x11*x70; POST_CHECK_MUL(x77);
        x78=(1.0/2.0)*x15*(x76 + 2*x77); POST_CHECK_MUL(x78);
        x79=x53*(x74 + x78); POST_CHECK_MUL(x79);
        x80=x50*(x72 + x79); POST_CHECK_MUL(x80);
        x81=x56*(x72 - x79); POST_CHECK_MUL(x81);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x165=-x74 - x78; POST_CHECK_ADD(x165);
        x171=x53*x73; POST_CHECK_MUL(x171);
        x188=x18*x53; POST_CHECK_MUL(x188);
        x190=x122*x18; POST_CHECK_MUL(x190);
        x198=x165*x190 + x171; POST_CHECK_ADD(x198);

        /* Assignment result[2, 2]=-x188*x80 + x188*x81 + x198*x25 - x198*x27 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x71=x14*x70; POST_CHECK_MUL(x71);
        x72=x8*ut1 + x4*x71; POST_CHECK_ADD(x72);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x74=x18*x73; POST_CHECK_MUL(x74);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x77=x11*x70; POST_CHECK_MUL(x77);
        x78=(1.0/2.0)*x15*(x76 + 2*x77); POST_CHECK_MUL(x78);
        x79=x53*(x74 + x78); POST_CHECK_MUL(x79);
        x80=x50*(x72 + x79); POST_CHECK_MUL(x80);
        x81=x56*(x72 - x79); POST_CHECK_MUL(x81);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x165=-x74 - x78; POST_CHECK_ADD(x165);
        x171=x53*x73; POST_CHECK_MUL(x171);
        x188=x18*x53; POST_CHECK_MUL(x188);
        x190=x122*x18; POST_CHECK_MUL(x190);
        x198=x165*x190 + x171; POST_CHECK_ADD(x198);
        result[8] = -x188*x80 + x188*x81 + x198*x25 - x198*x27;

    }
    else if (x125)
    {
        DEBUG_PRINT("Case (x125) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x71=x14*x70; POST_CHECK_MUL(x71);
        x72=x8*ut1 + x4*x71; POST_CHECK_ADD(x72);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x74=x18*x73; POST_CHECK_MUL(x74);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x77=x11*x70; POST_CHECK_MUL(x77);
        x78=(1.0/2.0)*x15*(x76 + 2*x77); POST_CHECK_MUL(x78);
        x79=x53*(x74 + x78); POST_CHECK_MUL(x79);
        x80=x50*(x72 + x79); POST_CHECK_MUL(x80);
        x81=x56*(x72 - x79); POST_CHECK_MUL(x81);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x189=random2*x113; POST_CHECK_MUL(x189);

        /* Assignment result[2, 2]=-x189*x80 + x189*x81 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x71=x14*x70; POST_CHECK_MUL(x71);
        x72=x8*ut1 + x4*x71; POST_CHECK_ADD(x72);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x74=x18*x73; POST_CHECK_MUL(x74);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x77=x11*x70; POST_CHECK_MUL(x77);
        x78=(1.0/2.0)*x15*(x76 + 2*x77); POST_CHECK_MUL(x78);
        x79=x53*(x74 + x78); POST_CHECK_MUL(x79);
        x80=x50*(x72 + x79); POST_CHECK_MUL(x80);
        x81=x56*(x72 - x79); POST_CHECK_MUL(x81);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x189=random2*x113; POST_CHECK_MUL(x189);
        result[8] = -x189*x80 + x189*x81;

    }


    /* Assignment result[0, 3]=Piecewise((x64, x39), (x67*(0.5*x44*x8*rt2*rn*un - x47*x8*rt2*rn + x69), x49), (x82 - x88 - x89, x61)) */
    double x82;
    double x83;
    double x84;
    double x85;
    double x86;
    double x87;
    double x88;
    double x89;
    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        x62=0.70710678118654757273731092936941422522068023681640625*mu; POST_CHECK_MUL(x62);
        x63=x29*x62; POST_CHECK_MUL(x63);
        x64=x30*(-x62 + x63 - 2.0*x8); POST_CHECK_MUL(x64);

        /* Assignment result[0, 3]=x64 */
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        x62=0.70710678118654757273731092936941422522068023681640625*mu; POST_CHECK_MUL(x62);
        x63=x29*x62; POST_CHECK_MUL(x63);
        x64=x30*(-x62 + x63 - 2.0*x8); POST_CHECK_MUL(x64);
        result[9] = x64;

    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (x42) != 0.;*/
        x43=1.0/x42; POST_CHECK_POW(x43);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        /*@ assert (x44) != 0.;*/
        x45=1.0/x44; POST_CHECK_POW(x45);
        x46=0.5*un; POST_CHECK_MUL(x46);
        /*@ assert (x46) >= 0.;*/
        /*@ assert (x46) != 0.;*/
        x47=x42*x46; POST_CHECK_MUL(x47);
        x48=x42*x44; POST_CHECK_MUL(x48);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x66=sqrt(x65); POST_CHECK_POW(x66);
        /*@ assert (x66) != 0.;*/
        x67=x43*x45/x66; POST_CHECK_MUL(x67);
        x68=0.353553390593273786368655464684707112610340118408203125*x66*mu*un; POST_CHECK_MUL(x68);
        x69=0.70710678118654757273731092936941422522068023681640625*x48*x66*mu - x42*x68 - x44*x68; POST_CHECK_ADD(x69);

        /* Assignment result[0, 3]=x67*(0.5*x44*x8*rt2*rn*un - x47*x8*rt2*rn + x69) */
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (x42) != 0.;*/
        x43=1.0/x42; POST_CHECK_POW(x43);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        /*@ assert (x44) != 0.;*/
        x45=1.0/x44; POST_CHECK_POW(x45);
        x46=0.5*un; POST_CHECK_MUL(x46);
        /*@ assert (x46) >= 0.;*/
        /*@ assert (x46) != 0.;*/
        x47=x42*x46; POST_CHECK_MUL(x47);
        x48=x42*x44; POST_CHECK_MUL(x48);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x66=sqrt(x65); POST_CHECK_POW(x66);
        /*@ assert (x66) != 0.;*/
        x67=x43*x45/x66; POST_CHECK_MUL(x67);
        x68=0.353553390593273786368655464684707112610340118408203125*x66*mu*un; POST_CHECK_MUL(x68);
        x69=0.70710678118654757273731092936941422522068023681640625*x48*x66*mu - x42*x68 - x44*x68; POST_CHECK_ADD(x69);
        result[9] = x67*(0.5*x44*x8*rt2*rn*un - x47*x8*rt2*rn + x69);

    }
    else if (x61)
    {
        DEBUG_PRINT("Case (x61) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x82=x17*x70; POST_CHECK_MUL(x82);
        x83=x8*ut2 + x4*x82; POST_CHECK_ADD(x83);
        x84=x15*x73; POST_CHECK_MUL(x84);
        x85=x12*x70; POST_CHECK_MUL(x85);
        x86=(1.0/2.0)*x18*(x76 + 2*x85); POST_CHECK_MUL(x86);
        x87=x53*(x84 + x86); POST_CHECK_MUL(x87);
        x88=x50*(x83 + x87); POST_CHECK_MUL(x88);
        x89=x56*(x83 - x87); POST_CHECK_MUL(x89);

        /* Assignment result[0, 3]=x82 - x88 - x89 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x82=x17*x70; POST_CHECK_MUL(x82);
        x83=x8*ut2 + x4*x82; POST_CHECK_ADD(x83);
        x84=x15*x73; POST_CHECK_MUL(x84);
        x85=x12*x70; POST_CHECK_MUL(x85);
        x86=(1.0/2.0)*x18*(x76 + 2*x85); POST_CHECK_MUL(x86);
        x87=x53*(x84 + x86); POST_CHECK_MUL(x87);
        x88=x50*(x83 + x87); POST_CHECK_MUL(x88);
        x89=x56*(x83 - x87); POST_CHECK_MUL(x89);
        result[9] = x82 - x88 - x89;

    }


    /* Assignment result[1, 3]=Piecewise((x167, x39), (x142*(x150 + x151 + x153 + x154 + x159 + x160 + x170), x49), (-x112*x88 + x112*x89 + x173*x25 - x173*x27, x61), (-x114*x88 + x114*x89, x125)) */
    double x172;
    double x173;
    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x35=1.4142135623730951454746218587388284504413604736328125*x28; POST_CHECK_MUL(x35);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x128=4.24264068711928477029005080112256109714508056640625*x28; POST_CHECK_MUL(x128);
        x129=4.24264068711928477029005080112256109714508056640625*x29; POST_CHECK_MUL(x129);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x131=2.0*x130; POST_CHECK_MUL(x131);
        x133=x118*x28; POST_CHECK_MUL(x133);
        x134=x129*x8; POST_CHECK_MUL(x134);
        x136=2.0*x118; POST_CHECK_MUL(x136);
        x137=x136*x29; POST_CHECK_MUL(x137);
        x167=x127*(x128 - x129 - x131 + 2.0*x133 + x134 - x137 + x35*x8 + 10.0*x90); POST_CHECK_MUL(x167);

        /* Assignment result[1, 3]=x167 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x35=1.4142135623730951454746218587388284504413604736328125*x28; POST_CHECK_MUL(x35);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x128=4.24264068711928477029005080112256109714508056640625*x28; POST_CHECK_MUL(x128);
        x129=4.24264068711928477029005080112256109714508056640625*x29; POST_CHECK_MUL(x129);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x131=2.0*x130; POST_CHECK_MUL(x131);
        x133=x118*x28; POST_CHECK_MUL(x133);
        x134=x129*x8; POST_CHECK_MUL(x134);
        x136=2.0*x118; POST_CHECK_MUL(x136);
        x137=x136*x29; POST_CHECK_MUL(x137);
        x167=x127*(x128 - x129 - x131 + 2.0*x133 + x134 - x137 + x35*x8 + 10.0*x90); POST_CHECK_MUL(x167);
        result[10] = x167;

    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x116=pow(x65, 3.0/2.0); POST_CHECK_POW(x116);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x119=pow(rn, 3); POST_CHECK_POW(x119);
        /*@ assert (x119) >= 0.;*/
        /*@ assert (x119) != 0.;*/
        x120=pow(rt1, 3); POST_CHECK_POW(x120);
        x139=4.0*x116*x42*x44; POST_CHECK_MUL(x139);
        x140=x139*x6; POST_CHECK_MUL(x140);
        x141=x139*x7; POST_CHECK_MUL(x141);
        /*@ assert (x140 + x141) != 0.;*/
        x142=1.0/(x140 + x141); POST_CHECK_POW(x142);
        x143=pow(mu, 5); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x144=pow(rn, 4); POST_CHECK_POW(x144);
        /*@ assert (x144) >= 0.;*/
        /*@ assert (x144) != 0.;*/
        x145=pow(un, 3); POST_CHECK_POW(x145);
        /*@ assert (x145) >= 0.;*/
        /*@ assert (x145) != 0.;*/
        x148=pow(mu, 4); POST_CHECK_POW(x148);
        /*@ assert (x148) >= 0.;*/
        /*@ assert (x148) != 0.;*/
        x149=pow(rt1, 5); POST_CHECK_POW(x149);
        x150=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x149*x42*un; POST_CHECK_MUL(x150);
        x151=1.4142135623730951454746218587388284504413604736328125*x119*x148*x149*x44*un; POST_CHECK_MUL(x151);
        x152=pow(rt2, 4); POST_CHECK_POW(x152);
        /*@ assert (x152) >= 0.;*/
        x153=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x152*x42*rt1*un; POST_CHECK_MUL(x153);
        x154=1.4142135623730951454746218587388284504413604736328125*x119*x148*x152*x44*rt1*un; POST_CHECK_MUL(x154);
        x159=-2.828427124746190290949243717477656900882720947265625*x119*x120*x148*x42*x7*un; POST_CHECK_MUL(x159);
        x160=2.828427124746190290949243717477656900882720947265625*x119*x120*x148*x44*x7*un; POST_CHECK_MUL(x160);
        x168=pow(rt2, 5); POST_CHECK_POW(x168);
        x169=pow(rt2, 3); POST_CHECK_POW(x169);
        x170=2.0*x116*x42*rt1*rt2*mu*un + 2.0*x116*x44*rt1*rt2*mu*un + 2.0*x118*x168*x42*x9*rt1*un - 2.0*x118*x168*x44*x9*rt1*un + 2.0*x143*x144*x169*x42*rt1*un - 2.0*x143*x144*x169*x44*rt1*un + 2.0*x118*x145*x169*x42*x9*rt1 - 2.0*x118*x145*x169*x44*x9*rt1 + 2.0*x118*x149*x42*x9*rt2*un - 2.0*x118*x149*x44*x9*rt2*un + 2.0*x120*x143*x144*x42*rt2*un - 2.0*x120*x143*x144*x44*rt2*un + 2.0*x118*x120*x145*x42*x9*rt2 - 2.0*x118*x120*x145*x44*x9*rt2 + 4.0*x118*x120*x169*x42*x9*un - 4.0*x118*x120*x169*x44*x9*un; POST_CHECK_ADD(x170);

        /* Assignment result[1, 3]=x142*(x150 + x151 + x153 + x154 + x159 + x160 + x170) */
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x116=pow(x65, 3.0/2.0); POST_CHECK_POW(x116);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x119=pow(rn, 3); POST_CHECK_POW(x119);
        /*@ assert (x119) >= 0.;*/
        /*@ assert (x119) != 0.;*/
        x120=pow(rt1, 3); POST_CHECK_POW(x120);
        x139=4.0*x116*x42*x44; POST_CHECK_MUL(x139);
        x140=x139*x6; POST_CHECK_MUL(x140);
        x141=x139*x7; POST_CHECK_MUL(x141);
        /*@ assert (x140 + x141) != 0.;*/
        x142=1.0/(x140 + x141); POST_CHECK_POW(x142);
        x143=pow(mu, 5); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x144=pow(rn, 4); POST_CHECK_POW(x144);
        /*@ assert (x144) >= 0.;*/
        /*@ assert (x144) != 0.;*/
        x145=pow(un, 3); POST_CHECK_POW(x145);
        /*@ assert (x145) >= 0.;*/
        /*@ assert (x145) != 0.;*/
        x148=pow(mu, 4); POST_CHECK_POW(x148);
        /*@ assert (x148) >= 0.;*/
        /*@ assert (x148) != 0.;*/
        x149=pow(rt1, 5); POST_CHECK_POW(x149);
        x150=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x149*x42*un; POST_CHECK_MUL(x150);
        x151=1.4142135623730951454746218587388284504413604736328125*x119*x148*x149*x44*un; POST_CHECK_MUL(x151);
        x152=pow(rt2, 4); POST_CHECK_POW(x152);
        /*@ assert (x152) >= 0.;*/
        x153=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x152*x42*rt1*un; POST_CHECK_MUL(x153);
        x154=1.4142135623730951454746218587388284504413604736328125*x119*x148*x152*x44*rt1*un; POST_CHECK_MUL(x154);
        x159=-2.828427124746190290949243717477656900882720947265625*x119*x120*x148*x42*x7*un; POST_CHECK_MUL(x159);
        x160=2.828427124746190290949243717477656900882720947265625*x119*x120*x148*x44*x7*un; POST_CHECK_MUL(x160);
        x168=pow(rt2, 5); POST_CHECK_POW(x168);
        x169=pow(rt2, 3); POST_CHECK_POW(x169);
        x170=2.0*x116*x42*rt1*rt2*mu*un + 2.0*x116*x44*rt1*rt2*mu*un + 2.0*x118*x168*x42*x9*rt1*un - 2.0*x118*x168*x44*x9*rt1*un + 2.0*x143*x144*x169*x42*rt1*un - 2.0*x143*x144*x169*x44*rt1*un + 2.0*x118*x145*x169*x42*x9*rt1 - 2.0*x118*x145*x169*x44*x9*rt1 + 2.0*x118*x149*x42*x9*rt2*un - 2.0*x118*x149*x44*x9*rt2*un + 2.0*x120*x143*x144*x42*rt2*un - 2.0*x120*x143*x144*x44*rt2*un + 2.0*x118*x120*x145*x42*x9*rt2 - 2.0*x118*x120*x145*x44*x9*rt2 + 4.0*x118*x120*x169*x42*x9*un - 4.0*x118*x120*x169*x44*x9*un; POST_CHECK_ADD(x170);
        result[10] = x142*(x150 + x151 + x153 + x154 + x159 + x160 + x170);

    }
    else if (x61)
    {
        DEBUG_PRINT("Case (x61) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x82=x17*x70; POST_CHECK_MUL(x82);
        x83=x8*ut2 + x4*x82; POST_CHECK_ADD(x83);
        x84=x15*x73; POST_CHECK_MUL(x84);
        x85=x12*x70; POST_CHECK_MUL(x85);
        x86=(1.0/2.0)*x18*(x76 + 2*x85); POST_CHECK_MUL(x86);
        x87=x53*(x84 + x86); POST_CHECK_MUL(x87);
        x88=x50*(x83 + x87); POST_CHECK_MUL(x88);
        x89=x56*(x83 - x87); POST_CHECK_MUL(x89);
        x112=x15*x53; POST_CHECK_MUL(x112);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x123=x122*x15; POST_CHECK_MUL(x123);
        x171=x53*x73; POST_CHECK_MUL(x171);
        x172=-x84 - x86; POST_CHECK_ADD(x172);
        x173=x123*x172 + x171; POST_CHECK_ADD(x173);

        /* Assignment result[1, 3]=-x112*x88 + x112*x89 + x173*x25 - x173*x27 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x82=x17*x70; POST_CHECK_MUL(x82);
        x83=x8*ut2 + x4*x82; POST_CHECK_ADD(x83);
        x84=x15*x73; POST_CHECK_MUL(x84);
        x85=x12*x70; POST_CHECK_MUL(x85);
        x86=(1.0/2.0)*x18*(x76 + 2*x85); POST_CHECK_MUL(x86);
        x87=x53*(x84 + x86); POST_CHECK_MUL(x87);
        x88=x50*(x83 + x87); POST_CHECK_MUL(x88);
        x89=x56*(x83 - x87); POST_CHECK_MUL(x89);
        x112=x15*x53; POST_CHECK_MUL(x112);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x123=x122*x15; POST_CHECK_MUL(x123);
        x171=x53*x73; POST_CHECK_MUL(x171);
        x172=-x84 - x86; POST_CHECK_ADD(x172);
        x173=x123*x172 + x171; POST_CHECK_ADD(x173);
        result[10] = -x112*x88 + x112*x89 + x173*x25 - x173*x27;

    }
    else if (x125)
    {
        DEBUG_PRINT("Case (x125) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x82=x17*x70; POST_CHECK_MUL(x82);
        x83=x8*ut2 + x4*x82; POST_CHECK_ADD(x83);
        x84=x15*x73; POST_CHECK_MUL(x84);
        x85=x12*x70; POST_CHECK_MUL(x85);
        x86=(1.0/2.0)*x18*(x76 + 2*x85); POST_CHECK_MUL(x86);
        x87=x53*(x84 + x86); POST_CHECK_MUL(x87);
        x88=x50*(x83 + x87); POST_CHECK_MUL(x88);
        x89=x56*(x83 - x87); POST_CHECK_MUL(x89);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);

        /* Assignment result[1, 3]=-x114*x88 + x114*x89 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x82=x17*x70; POST_CHECK_MUL(x82);
        x83=x8*ut2 + x4*x82; POST_CHECK_ADD(x83);
        x84=x15*x73; POST_CHECK_MUL(x84);
        x85=x12*x70; POST_CHECK_MUL(x85);
        x86=(1.0/2.0)*x18*(x76 + 2*x85); POST_CHECK_MUL(x86);
        x87=x53*(x84 + x86); POST_CHECK_MUL(x87);
        x88=x50*(x83 + x87); POST_CHECK_MUL(x88);
        x89=x56*(x83 - x87); POST_CHECK_MUL(x89);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);
        result[10] = -x114*x88 + x114*x89;

    }


    /* Assignment result[2, 3]=Piecewise((x138, x39), (x142*(2.0*x118*x152*x44*x6*x9*un + 4.0*x118*x157*x44*x7*x9*un - 2.0*x118*x199*x42*x9*un + 2.0*x118*x199*x44*x9*un - 2.0*x143*x144*x157*x42*un + 2.0*x143*x144*x157*x44*un - 2.0*x118*x145*x157*x42*x9 + 2.0*x118*x145*x157*x44*x9 - x136*x156 + x146 - x155*x158 - 4.0*x161 - 4.0*x162 - 2.0*x163 - 2.0*x164 + x192 + x193 + x194 + x195 + x196 + x197), x49), (mu - x188*x88 + x188*x89 + x200*x25 - x200*x27, x61), (mu - x189*x88 + x189*x89, x125)) */
    double x199;
    double x200;
    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x128=4.24264068711928477029005080112256109714508056640625*x28; POST_CHECK_MUL(x128);
        x129=4.24264068711928477029005080112256109714508056640625*x29; POST_CHECK_MUL(x129);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x131=2.0*x130; POST_CHECK_MUL(x131);
        x132=x28*x8; POST_CHECK_MUL(x132);
        x133=x118*x28; POST_CHECK_MUL(x133);
        x134=x129*x8; POST_CHECK_MUL(x134);
        x135=-x134; POST_CHECK_MUL(x135);
        x136=2.0*x118; POST_CHECK_MUL(x136);
        x137=x136*x29; POST_CHECK_MUL(x137);
        x138=x127*(-x128 + x129 + x131 - 57.9827560572968963015227927826344966888427734375*x132 + 11.3137084989847611637969748699106276035308837890625*x132*x29 - 34.0*x133 + x135 + x137 - 26.0*x90 + 16.0*x91); POST_CHECK_MUL(x138);

        /* Assignment result[2, 3]=x138 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x128=4.24264068711928477029005080112256109714508056640625*x28; POST_CHECK_MUL(x128);
        x129=4.24264068711928477029005080112256109714508056640625*x29; POST_CHECK_MUL(x129);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x131=2.0*x130; POST_CHECK_MUL(x131);
        x132=x28*x8; POST_CHECK_MUL(x132);
        x133=x118*x28; POST_CHECK_MUL(x133);
        x134=x129*x8; POST_CHECK_MUL(x134);
        x135=-x134; POST_CHECK_MUL(x135);
        x136=2.0*x118; POST_CHECK_MUL(x136);
        x137=x136*x29; POST_CHECK_MUL(x137);
        x138=x127*(-x128 + x129 + x131 - 57.9827560572968963015227927826344966888427734375*x132 + 11.3137084989847611637969748699106276035308837890625*x132*x29 - 34.0*x133 + x135 + x137 - 26.0*x90 + 16.0*x91); POST_CHECK_MUL(x138);
        result[11] = x138;

    }
    else if (x49)
    {
        DEBUG_PRINT("Case (x49) is True.\n");
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x116=pow(x65, 3.0/2.0); POST_CHECK_POW(x116);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x119=pow(rn, 3); POST_CHECK_POW(x119);
        /*@ assert (x119) >= 0.;*/
        /*@ assert (x119) != 0.;*/
        x136=2.0*x118; POST_CHECK_MUL(x136);
        x139=4.0*x116*x42*x44; POST_CHECK_MUL(x139);
        x140=x139*x6; POST_CHECK_MUL(x140);
        x141=x139*x7; POST_CHECK_MUL(x141);
        /*@ assert (x140 + x141) != 0.;*/
        x142=1.0/(x140 + x141); POST_CHECK_POW(x142);
        x143=pow(mu, 5); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x144=pow(rn, 4); POST_CHECK_POW(x144);
        /*@ assert (x144) >= 0.;*/
        /*@ assert (x144) != 0.;*/
        x145=pow(un, 3); POST_CHECK_POW(x145);
        /*@ assert (x145) >= 0.;*/
        /*@ assert (x145) != 0.;*/
        x146=x140*mu + x141*mu - 2.0*x143*x144*x42*x6*x7*un + 2.0*x143*x144*x44*x6*x7*un - 2.0*x118*x145*x42*x6*x7*x9 + 2.0*x118*x145*x44*x6*x7*x9; POST_CHECK_ADD(x146);
        x148=pow(mu, 4); POST_CHECK_POW(x148);
        /*@ assert (x148) >= 0.;*/
        /*@ assert (x148) != 0.;*/
        x152=pow(rt2, 4); POST_CHECK_POW(x152);
        /*@ assert (x152) >= 0.;*/
        x155=4.0*x118; POST_CHECK_MUL(x155);
        x156=x152*x42*x6*x9*un; POST_CHECK_MUL(x156);
        x157=pow(rt1, 4); POST_CHECK_POW(x157);
        /*@ assert (x157) >= 0.;*/
        x158=x157*x42*x7*x9*un; POST_CHECK_MUL(x158);
        x161=x116*x42*x6*mu*un; POST_CHECK_MUL(x161);
        x162=x116*x44*x6*mu*un; POST_CHECK_MUL(x162);
        x163=x116*x42*x7*mu*un; POST_CHECK_MUL(x163);
        x164=x116*x44*x7*mu*un; POST_CHECK_MUL(x164);
        x168=pow(rt2, 5); POST_CHECK_POW(x168);
        x169=pow(rt2, 3); POST_CHECK_POW(x169);
        x192=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x168*x42*un; POST_CHECK_MUL(x192);
        x193=1.4142135623730951454746218587388284504413604736328125*x119*x148*x168*x44*un; POST_CHECK_MUL(x193);
        x194=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x157*x42*rt2*un; POST_CHECK_MUL(x194);
        x195=1.4142135623730951454746218587388284504413604736328125*x119*x148*x157*x44*rt2*un; POST_CHECK_MUL(x195);
        x196=-2.828427124746190290949243717477656900882720947265625*x119*x148*x169*x42*x6*un; POST_CHECK_MUL(x196);
        x197=2.828427124746190290949243717477656900882720947265625*x119*x148*x169*x44*x6*un; POST_CHECK_MUL(x197);
        x199=pow(rt1, 6); POST_CHECK_POW(x199);

        /* Assignment result[2, 3]=x142*(2.0*x118*x152*x44*x6*x9*un + 4.0*x118*x157*x44*x7*x9*un - 2.0*x118*x199*x42*x9*un + 2.0*x118*x199*x44*x9*un - 2.0*x143*x144*x157*x42*un + 2.0*x143*x144*x157*x44*un - 2.0*x118*x145*x157*x42*x9 + 2.0*x118*x145*x157*x44*x9 - x136*x156 + x146 - x155*x158 - 4.0*x161 - 4.0*x162 - 2.0*x163 - 2.0*x164 + x192 + x193 + x194 + x195 + x196 + x197) */
        /*@ assert (x6 + x7) >= 0.;*/
        x40=sqrt(x6 + x7); POST_CHECK_POW(x40);
        x41=un*un + x6 + x7 + x10; POST_CHECK_ADD(x41);
        /*@ assert (-2.0*x40*x5 + x41) >= 0.;*/
        x42=sqrt(-2.0*x40*x5 + x41); POST_CHECK_POW(x42);
        /*@ assert (2.0*x40*mu*rn + x41) >= 0.;*/
        x44=sqrt(2.0*x40*mu*rn + x41); POST_CHECK_POW(x44);
        x65=x6*x8*x9 + x7*x8*x9; POST_CHECK_ADD(x65);
        /*@ assert (x65) >= 0.;*/
        x116=pow(x65, 3.0/2.0); POST_CHECK_POW(x116);
        x118=pow(mu, 3); POST_CHECK_POW(x118);
        /*@ assert (x118) >= 0.;*/
        /*@ assert (x118) != 0.;*/
        x119=pow(rn, 3); POST_CHECK_POW(x119);
        /*@ assert (x119) >= 0.;*/
        /*@ assert (x119) != 0.;*/
        x136=2.0*x118; POST_CHECK_MUL(x136);
        x139=4.0*x116*x42*x44; POST_CHECK_MUL(x139);
        x140=x139*x6; POST_CHECK_MUL(x140);
        x141=x139*x7; POST_CHECK_MUL(x141);
        /*@ assert (x140 + x141) != 0.;*/
        x142=1.0/(x140 + x141); POST_CHECK_POW(x142);
        x143=pow(mu, 5); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x144=pow(rn, 4); POST_CHECK_POW(x144);
        /*@ assert (x144) >= 0.;*/
        /*@ assert (x144) != 0.;*/
        x145=pow(un, 3); POST_CHECK_POW(x145);
        /*@ assert (x145) >= 0.;*/
        /*@ assert (x145) != 0.;*/
        x146=x140*mu + x141*mu - 2.0*x143*x144*x42*x6*x7*un + 2.0*x143*x144*x44*x6*x7*un - 2.0*x118*x145*x42*x6*x7*x9 + 2.0*x118*x145*x44*x6*x7*x9; POST_CHECK_ADD(x146);
        x148=pow(mu, 4); POST_CHECK_POW(x148);
        /*@ assert (x148) >= 0.;*/
        /*@ assert (x148) != 0.;*/
        x152=pow(rt2, 4); POST_CHECK_POW(x152);
        /*@ assert (x152) >= 0.;*/
        x155=4.0*x118; POST_CHECK_MUL(x155);
        x156=x152*x42*x6*x9*un; POST_CHECK_MUL(x156);
        x157=pow(rt1, 4); POST_CHECK_POW(x157);
        /*@ assert (x157) >= 0.;*/
        x158=x157*x42*x7*x9*un; POST_CHECK_MUL(x158);
        x161=x116*x42*x6*mu*un; POST_CHECK_MUL(x161);
        x162=x116*x44*x6*mu*un; POST_CHECK_MUL(x162);
        x163=x116*x42*x7*mu*un; POST_CHECK_MUL(x163);
        x164=x116*x44*x7*mu*un; POST_CHECK_MUL(x164);
        x168=pow(rt2, 5); POST_CHECK_POW(x168);
        x169=pow(rt2, 3); POST_CHECK_POW(x169);
        x192=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x168*x42*un; POST_CHECK_MUL(x192);
        x193=1.4142135623730951454746218587388284504413604736328125*x119*x148*x168*x44*un; POST_CHECK_MUL(x193);
        x194=-1.4142135623730951454746218587388284504413604736328125*x119*x148*x157*x42*rt2*un; POST_CHECK_MUL(x194);
        x195=1.4142135623730951454746218587388284504413604736328125*x119*x148*x157*x44*rt2*un; POST_CHECK_MUL(x195);
        x196=-2.828427124746190290949243717477656900882720947265625*x119*x148*x169*x42*x6*un; POST_CHECK_MUL(x196);
        x197=2.828427124746190290949243717477656900882720947265625*x119*x148*x169*x44*x6*un; POST_CHECK_MUL(x197);
        x199=pow(rt1, 6); POST_CHECK_POW(x199);
        result[11] = x142*(2.0*x118*x152*x44*x6*x9*un + 4.0*x118*x157*x44*x7*x9*un - 2.0*x118*x199*x42*x9*un + 2.0*x118*x199*x44*x9*un - 2.0*x143*x144*x157*x42*un + 2.0*x143*x144*x157*x44*un - 2.0*x118*x145*x157*x42*x9 + 2.0*x118*x145*x157*x44*x9 - x136*x156 + x146 - x155*x158 - 4.0*x161 - 4.0*x162 - 2.0*x163 - 2.0*x164 + x192 + x193 + x194 + x195 + x196 + x197);

    }
    else if (x61)
    {
        DEBUG_PRINT("Case (x61) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x82=x17*x70; POST_CHECK_MUL(x82);
        x83=x8*ut2 + x4*x82; POST_CHECK_ADD(x83);
        x84=x15*x73; POST_CHECK_MUL(x84);
        x85=x12*x70; POST_CHECK_MUL(x85);
        x86=(1.0/2.0)*x18*(x76 + 2*x85); POST_CHECK_MUL(x86);
        x87=x53*(x84 + x86); POST_CHECK_MUL(x87);
        x88=x50*(x83 + x87); POST_CHECK_MUL(x88);
        x89=x56*(x83 - x87); POST_CHECK_MUL(x89);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x172=-x84 - x86; POST_CHECK_ADD(x172);
        x188=x18*x53; POST_CHECK_MUL(x188);
        x190=x122*x18; POST_CHECK_MUL(x190);
        x200=x53*(x75 + x85) + x172*x190; POST_CHECK_ADD(x200);

        /* Assignment result[2, 3]=mu - x188*x88 + x188*x89 + x200*x25 - x200*x27 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x82=x17*x70; POST_CHECK_MUL(x82);
        x83=x8*ut2 + x4*x82; POST_CHECK_ADD(x83);
        x84=x15*x73; POST_CHECK_MUL(x84);
        x85=x12*x70; POST_CHECK_MUL(x85);
        x86=(1.0/2.0)*x18*(x76 + 2*x85); POST_CHECK_MUL(x86);
        x87=x53*(x84 + x86); POST_CHECK_MUL(x87);
        x88=x50*(x83 + x87); POST_CHECK_MUL(x88);
        x89=x56*(x83 - x87); POST_CHECK_MUL(x89);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x172=-x84 - x86; POST_CHECK_ADD(x172);
        x188=x18*x53; POST_CHECK_MUL(x188);
        x190=x122*x18; POST_CHECK_MUL(x190);
        x200=x53*(x75 + x85) + x172*x190; POST_CHECK_ADD(x200);
        result[11] = mu - x188*x88 + x188*x89 + x200*x25 - x200*x27;

    }
    else if (x125)
    {
        DEBUG_PRINT("Case (x125) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x82=x17*x70; POST_CHECK_MUL(x82);
        x83=x8*ut2 + x4*x82; POST_CHECK_ADD(x83);
        x84=x15*x73; POST_CHECK_MUL(x84);
        x85=x12*x70; POST_CHECK_MUL(x85);
        x86=(1.0/2.0)*x18*(x76 + 2*x85); POST_CHECK_MUL(x86);
        x87=x53*(x84 + x86); POST_CHECK_MUL(x87);
        x88=x50*(x83 + x87); POST_CHECK_MUL(x88);
        x89=x56*(x83 - x87); POST_CHECK_MUL(x89);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x189=random2*x113; POST_CHECK_MUL(x189);

        /* Assignment result[2, 3]=mu - x189*x88 + x189*x89 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        /*@ assert (x3) != 0.;*/
        x70=1.0/x3; POST_CHECK_POW(x70);
        x73=x70*x8*ut1*ut2; POST_CHECK_MUL(x73);
        x75=x4*mu; POST_CHECK_MUL(x75);
        x76=2*x75; POST_CHECK_MUL(x76);
        x82=x17*x70; POST_CHECK_MUL(x82);
        x83=x8*ut2 + x4*x82; POST_CHECK_ADD(x83);
        x84=x15*x73; POST_CHECK_MUL(x84);
        x85=x12*x70; POST_CHECK_MUL(x85);
        x86=(1.0/2.0)*x18*(x76 + 2*x85); POST_CHECK_MUL(x86);
        x87=x53*(x84 + x86); POST_CHECK_MUL(x87);
        x88=x50*(x83 + x87); POST_CHECK_MUL(x88);
        x89=x56*(x83 - x87); POST_CHECK_MUL(x89);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x189=random2*x113; POST_CHECK_MUL(x189);
        result[11] = mu - x189*x88 + x189*x89;

    }


    /* Assignment result[0, 4]=Piecewise((x31*(-x28*x62 - x32*x8 - x33*x8 + x63 + x91), x39), (mu - x98 - x99, x100)) */

    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x32=0.5*x28; POST_CHECK_MUL(x32);
        x33=0.5*x29; POST_CHECK_MUL(x33);
        x62=0.70710678118654757273731092936941422522068023681640625*mu; POST_CHECK_MUL(x62);
        x63=x29*x62; POST_CHECK_MUL(x63);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);

        /* Assignment result[0, 4]=x31*(-x28*x62 - x32*x8 - x33*x8 + x63 + x91) */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x32=0.5*x28; POST_CHECK_MUL(x32);
        x33=0.5*x29; POST_CHECK_MUL(x33);
        x62=0.70710678118654757273731092936941422522068023681640625*mu; POST_CHECK_MUL(x62);
        x63=x29*x62; POST_CHECK_MUL(x63);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        result[12] = x31*(-x28*x62 - x32*x8 - x33*x8 + x63 + x91);

    }
    else if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x92=x8*rn; POST_CHECK_MUL(x92);
        x93=mu*rt1; POST_CHECK_MUL(x93);
        x94=x15*x93; POST_CHECK_MUL(x94);
        x95=mu*rt2; POST_CHECK_MUL(x95);
        x96=x18*x95; POST_CHECK_MUL(x96);
        x97=x53*(x94 + x96); POST_CHECK_MUL(x97);
        x98=x50*(x92 + x97); POST_CHECK_MUL(x98);
        x99=x56*(x92 - x97); POST_CHECK_MUL(x99);

        /* Assignment result[0, 4]=mu - x98 - x99 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x92=x8*rn; POST_CHECK_MUL(x92);
        x93=mu*rt1; POST_CHECK_MUL(x93);
        x94=x15*x93; POST_CHECK_MUL(x94);
        x95=mu*rt2; POST_CHECK_MUL(x95);
        x96=x18*x95; POST_CHECK_MUL(x96);
        x97=x53*(x94 + x96); POST_CHECK_MUL(x97);
        x98=x50*(x92 + x97); POST_CHECK_MUL(x98);
        x99=x56*(x92 - x97); POST_CHECK_MUL(x99);
        result[12] = mu - x98 - x99;

    }


    /* Assignment result[1, 4]=Piecewise((x174, x39), (-x112*x98 + x112*x99 + x176*x25 - x176*x27, x100), (-x114*x98 + x114*x99, x177)) */
    double x174;
    double x175;
    double x176;
    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x32=0.5*x28; POST_CHECK_MUL(x32);
        x33=0.5*x29; POST_CHECK_MUL(x33);
        x101=0.353553390593273786368655464684707112610340118408203125*x28; POST_CHECK_MUL(x101);
        x102=0.353553390593273786368655464684707112610340118408203125*x29; POST_CHECK_MUL(x102);
        x174=x31*(-x32*mu - x33*mu - x101*x8 + x102*x8); POST_CHECK_MUL(x174);

        /* Assignment result[1, 4]=x174 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x32=0.5*x28; POST_CHECK_MUL(x32);
        x33=0.5*x29; POST_CHECK_MUL(x33);
        x101=0.353553390593273786368655464684707112610340118408203125*x28; POST_CHECK_MUL(x101);
        x102=0.353553390593273786368655464684707112610340118408203125*x29; POST_CHECK_MUL(x102);
        x174=x31*(-x32*mu - x33*mu - x101*x8 + x102*x8); POST_CHECK_MUL(x174);
        result[13] = x174;

    }
    else if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x92=x8*rn; POST_CHECK_MUL(x92);
        x93=mu*rt1; POST_CHECK_MUL(x93);
        x94=x15*x93; POST_CHECK_MUL(x94);
        x95=mu*rt2; POST_CHECK_MUL(x95);
        x96=x18*x95; POST_CHECK_MUL(x96);
        x97=x53*(x94 + x96); POST_CHECK_MUL(x97);
        x98=x50*(x92 + x97); POST_CHECK_MUL(x98);
        x99=x56*(x92 - x97); POST_CHECK_MUL(x99);
        x112=x15*x53; POST_CHECK_MUL(x112);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x123=x122*x15; POST_CHECK_MUL(x123);
        x175=-x94 - x96; POST_CHECK_ADD(x175);
        x176=x123*x175 + x53*x93; POST_CHECK_ADD(x176);

        /* Assignment result[1, 4]=-x112*x98 + x112*x99 + x176*x25 - x176*x27 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x92=x8*rn; POST_CHECK_MUL(x92);
        x93=mu*rt1; POST_CHECK_MUL(x93);
        x94=x15*x93; POST_CHECK_MUL(x94);
        x95=mu*rt2; POST_CHECK_MUL(x95);
        x96=x18*x95; POST_CHECK_MUL(x96);
        x97=x53*(x94 + x96); POST_CHECK_MUL(x97);
        x98=x50*(x92 + x97); POST_CHECK_MUL(x98);
        x99=x56*(x92 - x97); POST_CHECK_MUL(x99);
        x112=x15*x53; POST_CHECK_MUL(x112);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x123=x122*x15; POST_CHECK_MUL(x123);
        x175=-x94 - x96; POST_CHECK_ADD(x175);
        x176=x123*x175 + x53*x93; POST_CHECK_ADD(x176);
        result[13] = -x112*x98 + x112*x99 + x176*x25 - x176*x27;

    }
    else if (x177)
    {
        DEBUG_PRINT("Case (x177) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x92=x8*rn; POST_CHECK_MUL(x92);
        x93=mu*rt1; POST_CHECK_MUL(x93);
        x94=x15*x93; POST_CHECK_MUL(x94);
        x95=mu*rt2; POST_CHECK_MUL(x95);
        x96=x18*x95; POST_CHECK_MUL(x96);
        x97=x53*(x94 + x96); POST_CHECK_MUL(x97);
        x98=x50*(x92 + x97); POST_CHECK_MUL(x98);
        x99=x56*(x92 - x97); POST_CHECK_MUL(x99);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);

        /* Assignment result[1, 4]=-x114*x98 + x114*x99 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x92=x8*rn; POST_CHECK_MUL(x92);
        x93=mu*rt1; POST_CHECK_MUL(x93);
        x94=x15*x93; POST_CHECK_MUL(x94);
        x95=mu*rt2; POST_CHECK_MUL(x95);
        x96=x18*x95; POST_CHECK_MUL(x96);
        x97=x53*(x94 + x96); POST_CHECK_MUL(x97);
        x98=x50*(x92 + x97); POST_CHECK_MUL(x98);
        x99=x56*(x92 - x97); POST_CHECK_MUL(x99);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);
        result[13] = -x114*x98 + x114*x99;

    }


    /* Assignment result[2, 4]=Piecewise((x174, x39), (-x188*x98 + x188*x99 + x201*x25 - x201*x27, x100), (-x189*x98 + x189*x99, x177)) */
    double x201;
    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x32=0.5*x28; POST_CHECK_MUL(x32);
        x33=0.5*x29; POST_CHECK_MUL(x33);
        x101=0.353553390593273786368655464684707112610340118408203125*x28; POST_CHECK_MUL(x101);
        x102=0.353553390593273786368655464684707112610340118408203125*x29; POST_CHECK_MUL(x102);
        x174=x31*(-x32*mu - x33*mu - x101*x8 + x102*x8); POST_CHECK_MUL(x174);

        /* Assignment result[2, 4]=x174 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x32=0.5*x28; POST_CHECK_MUL(x32);
        x33=0.5*x29; POST_CHECK_MUL(x33);
        x101=0.353553390593273786368655464684707112610340118408203125*x28; POST_CHECK_MUL(x101);
        x102=0.353553390593273786368655464684707112610340118408203125*x29; POST_CHECK_MUL(x102);
        x174=x31*(-x32*mu - x33*mu - x101*x8 + x102*x8); POST_CHECK_MUL(x174);
        result[14] = x174;

    }
    else if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x92=x8*rn; POST_CHECK_MUL(x92);
        x93=mu*rt1; POST_CHECK_MUL(x93);
        x94=x15*x93; POST_CHECK_MUL(x94);
        x95=mu*rt2; POST_CHECK_MUL(x95);
        x96=x18*x95; POST_CHECK_MUL(x96);
        x97=x53*(x94 + x96); POST_CHECK_MUL(x97);
        x98=x50*(x92 + x97); POST_CHECK_MUL(x98);
        x99=x56*(x92 - x97); POST_CHECK_MUL(x99);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x175=-x94 - x96; POST_CHECK_ADD(x175);
        x188=x18*x53; POST_CHECK_MUL(x188);
        x190=x122*x18; POST_CHECK_MUL(x190);
        x201=x175*x190 + x53*x95; POST_CHECK_ADD(x201);

        /* Assignment result[2, 4]=-x188*x98 + x188*x99 + x201*x25 - x201*x27 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x92=x8*rn; POST_CHECK_MUL(x92);
        x93=mu*rt1; POST_CHECK_MUL(x93);
        x94=x15*x93; POST_CHECK_MUL(x94);
        x95=mu*rt2; POST_CHECK_MUL(x95);
        x96=x18*x95; POST_CHECK_MUL(x96);
        x97=x53*(x94 + x96); POST_CHECK_MUL(x97);
        x98=x50*(x92 + x97); POST_CHECK_MUL(x98);
        x99=x56*(x92 - x97); POST_CHECK_MUL(x99);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x175=-x94 - x96; POST_CHECK_ADD(x175);
        x188=x18*x53; POST_CHECK_MUL(x188);
        x190=x122*x18; POST_CHECK_MUL(x190);
        x201=x175*x190 + x53*x95; POST_CHECK_ADD(x201);
        result[14] = -x188*x98 + x188*x99 + x201*x25 - x201*x27;

    }
    else if (x177)
    {
        DEBUG_PRINT("Case (x177) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x92=x8*rn; POST_CHECK_MUL(x92);
        x93=mu*rt1; POST_CHECK_MUL(x93);
        x94=x15*x93; POST_CHECK_MUL(x94);
        x95=mu*rt2; POST_CHECK_MUL(x95);
        x96=x18*x95; POST_CHECK_MUL(x96);
        x97=x53*(x94 + x96); POST_CHECK_MUL(x97);
        x98=x50*(x92 + x97); POST_CHECK_MUL(x98);
        x99=x56*(x92 - x97); POST_CHECK_MUL(x99);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x189=random2*x113; POST_CHECK_MUL(x189);

        /* Assignment result[2, 4]=-x189*x98 + x189*x99 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x92=x8*rn; POST_CHECK_MUL(x92);
        x93=mu*rt1; POST_CHECK_MUL(x93);
        x94=x15*x93; POST_CHECK_MUL(x94);
        x95=mu*rt2; POST_CHECK_MUL(x95);
        x96=x18*x95; POST_CHECK_MUL(x96);
        x97=x53*(x94 + x96); POST_CHECK_MUL(x97);
        x98=x50*(x92 + x97); POST_CHECK_MUL(x98);
        x99=x56*(x92 - x97); POST_CHECK_MUL(x99);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x189=random2*x113; POST_CHECK_MUL(x189);
        result[14] = -x189*x98 + x189*x99;

    }


    /* Assignment result[0, 5]=Piecewise((x103, x39), (-x106 - x107, x100)) */
    double x103;
    double x104;
    double x105;
    double x106;
    double x107;
    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x32=0.5*x28; POST_CHECK_MUL(x32);
        x33=0.5*x29; POST_CHECK_MUL(x33);
        x34=-x32 - x33; POST_CHECK_ADD(x34);
        x101=0.353553390593273786368655464684707112610340118408203125*x28; POST_CHECK_MUL(x101);
        x102=0.353553390593273786368655464684707112610340118408203125*x29; POST_CHECK_MUL(x102);
        x103=x31*(-x101*mu + x102*mu + x34); POST_CHECK_MUL(x103);

        /* Assignment result[0, 5]=x103 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x32=0.5*x28; POST_CHECK_MUL(x32);
        x33=0.5*x29; POST_CHECK_MUL(x33);
        x34=-x32 - x33; POST_CHECK_ADD(x34);
        x101=0.353553390593273786368655464684707112610340118408203125*x28; POST_CHECK_MUL(x101);
        x102=0.353553390593273786368655464684707112610340118408203125*x29; POST_CHECK_MUL(x102);
        x103=x31*(-x101*mu + x102*mu + x34); POST_CHECK_MUL(x103);
        result[15] = x103;

    }
    else if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x105=x104*x15; POST_CHECK_MUL(x105);
        x106=x50*(rt1 + x105); POST_CHECK_MUL(x106);
        x107=x56*(rt1 - x105); POST_CHECK_MUL(x107);

        /* Assignment result[0, 5]=-x106 - x107 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x105=x104*x15; POST_CHECK_MUL(x105);
        x106=x50*(rt1 + x105); POST_CHECK_MUL(x106);
        x107=x56*(rt1 - x105); POST_CHECK_MUL(x107);
        result[15] = -x106 - x107;

    }


    /* Assignment result[1, 5]=Piecewise((x181, x39), (1 - x106*x112 + x107*x112 + x183*x25 - x183*x27, x100), (1 - x106*x114 + x107*x114, x177)) */
    double x178;
    double x179;
    double x180;
    double x181;
    double x182;
    double x183;
    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x132=x28*x8; POST_CHECK_MUL(x132);
        x178=9.8994949366116653521885382360778748989105224609375*x28; POST_CHECK_MUL(x178);
        x179=4.0*x130; POST_CHECK_MUL(x179);
        x180=1.4142135623730951454746218587388284504413604736328125*x29; POST_CHECK_MUL(x180);
        x181=x127*(x126 - 15.5563491861040450459086059709079563617706298828125*x132 - x178 - x179 - x180*x8 + 9.8994949366116653521885382360778748989105224609375*x29 - 20.0*x90); POST_CHECK_MUL(x181);

        /* Assignment result[1, 5]=x181 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x132=x28*x8; POST_CHECK_MUL(x132);
        x178=9.8994949366116653521885382360778748989105224609375*x28; POST_CHECK_MUL(x178);
        x179=4.0*x130; POST_CHECK_MUL(x179);
        x180=1.4142135623730951454746218587388284504413604736328125*x29; POST_CHECK_MUL(x180);
        x181=x127*(x126 - 15.5563491861040450459086059709079563617706298828125*x132 - x178 - x179 - x180*x8 + 9.8994949366116653521885382360778748989105224609375*x29 - 20.0*x90); POST_CHECK_MUL(x181);
        result[16] = x181;

    }
    else if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x105=x104*x15; POST_CHECK_MUL(x105);
        x106=x50*(rt1 + x105); POST_CHECK_MUL(x106);
        x107=x56*(rt1 - x105); POST_CHECK_MUL(x107);
        x112=x15*x53; POST_CHECK_MUL(x112);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x182=x122*mu*rn; POST_CHECK_MUL(x182);
        x183=x104 - x16*x182; POST_CHECK_ADD(x183);

        /* Assignment result[1, 5]=1 - x106*x112 + x107*x112 + x183*x25 - x183*x27 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x105=x104*x15; POST_CHECK_MUL(x105);
        x106=x50*(rt1 + x105); POST_CHECK_MUL(x106);
        x107=x56*(rt1 - x105); POST_CHECK_MUL(x107);
        x112=x15*x53; POST_CHECK_MUL(x112);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x182=x122*mu*rn; POST_CHECK_MUL(x182);
        x183=x104 - x16*x182; POST_CHECK_ADD(x183);
        result[16] = 1 - x106*x112 + x107*x112 + x183*x25 - x183*x27;

    }
    else if (x177)
    {
        DEBUG_PRINT("Case (x177) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x105=x104*x15; POST_CHECK_MUL(x105);
        x106=x50*(rt1 + x105); POST_CHECK_MUL(x106);
        x107=x56*(rt1 - x105); POST_CHECK_MUL(x107);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);

        /* Assignment result[1, 5]=1 - x106*x114 + x107*x114 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x105=x104*x15; POST_CHECK_MUL(x105);
        x106=x50*(rt1 + x105); POST_CHECK_MUL(x106);
        x107=x56*(rt1 - x105); POST_CHECK_MUL(x107);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);
        result[16] = 1 - x106*x114 + x107*x114;

    }


    /* Assignment result[2, 5]=Piecewise((x184, x39), (-x106*x188 + x107*x188 + x186, x100), (-x106*x189 + x107*x189, x177)) */
    double x184;
    double x185;
    double x186;
    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x35=1.4142135623730951454746218587388284504413604736328125*x28; POST_CHECK_MUL(x35);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x129=4.24264068711928477029005080112256109714508056640625*x29; POST_CHECK_MUL(x129);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x134=x129*x8; POST_CHECK_MUL(x134);
        x135=-x134; POST_CHECK_MUL(x135);
        x178=9.8994949366116653521885382360778748989105224609375*x28; POST_CHECK_MUL(x178);
        x179=4.0*x130; POST_CHECK_MUL(x179);
        x180=1.4142135623730951454746218587388284504413604736328125*x29; POST_CHECK_MUL(x180);
        x184=x127*(x135 + x178*x8 + x179 + x180 - x35 + 4.0*x90); POST_CHECK_MUL(x184);

        /* Assignment result[2, 5]=x184 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x35=1.4142135623730951454746218587388284504413604736328125*x28; POST_CHECK_MUL(x35);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x129=4.24264068711928477029005080112256109714508056640625*x29; POST_CHECK_MUL(x129);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x134=x129*x8; POST_CHECK_MUL(x134);
        x135=-x134; POST_CHECK_MUL(x135);
        x178=9.8994949366116653521885382360778748989105224609375*x28; POST_CHECK_MUL(x178);
        x179=4.0*x130; POST_CHECK_MUL(x179);
        x180=1.4142135623730951454746218587388284504413604736328125*x29; POST_CHECK_MUL(x180);
        x184=x127*(x135 + x178*x8 + x179 + x180 - x35 + 4.0*x90); POST_CHECK_MUL(x184);
        result[17] = x184;

    }
    else if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x105=x104*x15; POST_CHECK_MUL(x105);
        x106=x50*(rt1 + x105); POST_CHECK_MUL(x106);
        x107=x56*(rt1 - x105); POST_CHECK_MUL(x107);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x185=x122*x15*x18*mu*rn; POST_CHECK_MUL(x185);
        x186=-x185*x25 + x185*x27; POST_CHECK_ADD(x186);
        x188=x18*x53; POST_CHECK_MUL(x188);

        /* Assignment result[2, 5]=-x106*x188 + x107*x188 + x186 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x105=x104*x15; POST_CHECK_MUL(x105);
        x106=x50*(rt1 + x105); POST_CHECK_MUL(x106);
        x107=x56*(rt1 - x105); POST_CHECK_MUL(x107);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x185=x122*x15*x18*mu*rn; POST_CHECK_MUL(x185);
        x186=-x185*x25 + x185*x27; POST_CHECK_ADD(x186);
        x188=x18*x53; POST_CHECK_MUL(x188);
        result[17] = -x106*x188 + x107*x188 + x186;

    }
    else if (x177)
    {
        DEBUG_PRINT("Case (x177) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x105=x104*x15; POST_CHECK_MUL(x105);
        x106=x50*(rt1 + x105); POST_CHECK_MUL(x106);
        x107=x56*(rt1 - x105); POST_CHECK_MUL(x107);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x189=random2*x113; POST_CHECK_MUL(x189);

        /* Assignment result[2, 5]=-x106*x189 + x107*x189 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x105=x104*x15; POST_CHECK_MUL(x105);
        x106=x50*(rt1 + x105); POST_CHECK_MUL(x106);
        x107=x56*(rt1 - x105); POST_CHECK_MUL(x107);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x189=random2*x113; POST_CHECK_MUL(x189);
        result[17] = -x106*x189 + x107*x189;

    }


    /* Assignment result[0, 6]=Piecewise((x103, x39), (-x109 - x110, x100)) */
    double x108;
    double x109;
    double x110;
    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x32=0.5*x28; POST_CHECK_MUL(x32);
        x33=0.5*x29; POST_CHECK_MUL(x33);
        x34=-x32 - x33; POST_CHECK_ADD(x34);
        x101=0.353553390593273786368655464684707112610340118408203125*x28; POST_CHECK_MUL(x101);
        x102=0.353553390593273786368655464684707112610340118408203125*x29; POST_CHECK_MUL(x102);
        x103=x31*(-x101*mu + x102*mu + x34); POST_CHECK_MUL(x103);

        /* Assignment result[0, 6]=x103 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        /*@ assert (x29) != 0.;*/
        x30=1.0/x29; POST_CHECK_POW(x30);
        /*@ assert (x28) != 0.;*/
        x31=x30/x28; POST_CHECK_MUL(x31);
        x32=0.5*x28; POST_CHECK_MUL(x32);
        x33=0.5*x29; POST_CHECK_MUL(x33);
        x34=-x32 - x33; POST_CHECK_ADD(x34);
        x101=0.353553390593273786368655464684707112610340118408203125*x28; POST_CHECK_MUL(x101);
        x102=0.353553390593273786368655464684707112610340118408203125*x29; POST_CHECK_MUL(x102);
        x103=x31*(-x101*mu + x102*mu + x34); POST_CHECK_MUL(x103);
        result[18] = x103;

    }
    else if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x108=x104*x18; POST_CHECK_MUL(x108);
        x109=x50*(rt2 + x108); POST_CHECK_MUL(x109);
        x110=x56*(rt2 - x108); POST_CHECK_MUL(x110);

        /* Assignment result[0, 6]=-x109 - x110 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x108=x104*x18; POST_CHECK_MUL(x108);
        x109=x50*(rt2 + x108); POST_CHECK_MUL(x109);
        x110=x56*(rt2 - x108); POST_CHECK_MUL(x110);
        result[18] = -x109 - x110;

    }


    /* Assignment result[1, 6]=Piecewise((x184, x39), (-x109*x112 + x110*x112 + x186, x100), (-x109*x114 + x110*x114, x177)) */

    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x35=1.4142135623730951454746218587388284504413604736328125*x28; POST_CHECK_MUL(x35);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x129=4.24264068711928477029005080112256109714508056640625*x29; POST_CHECK_MUL(x129);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x134=x129*x8; POST_CHECK_MUL(x134);
        x135=-x134; POST_CHECK_MUL(x135);
        x178=9.8994949366116653521885382360778748989105224609375*x28; POST_CHECK_MUL(x178);
        x179=4.0*x130; POST_CHECK_MUL(x179);
        x180=1.4142135623730951454746218587388284504413604736328125*x29; POST_CHECK_MUL(x180);
        x184=x127*(x135 + x178*x8 + x179 + x180 - x35 + 4.0*x90); POST_CHECK_MUL(x184);

        /* Assignment result[1, 6]=x184 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x35=1.4142135623730951454746218587388284504413604736328125*x28; POST_CHECK_MUL(x35);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x129=4.24264068711928477029005080112256109714508056640625*x29; POST_CHECK_MUL(x129);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x134=x129*x8; POST_CHECK_MUL(x134);
        x135=-x134; POST_CHECK_MUL(x135);
        x178=9.8994949366116653521885382360778748989105224609375*x28; POST_CHECK_MUL(x178);
        x179=4.0*x130; POST_CHECK_MUL(x179);
        x180=1.4142135623730951454746218587388284504413604736328125*x29; POST_CHECK_MUL(x180);
        x184=x127*(x135 + x178*x8 + x179 + x180 - x35 + 4.0*x90); POST_CHECK_MUL(x184);
        result[19] = x184;

    }
    else if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x108=x104*x18; POST_CHECK_MUL(x108);
        x109=x50*(rt2 + x108); POST_CHECK_MUL(x109);
        x110=x56*(rt2 - x108); POST_CHECK_MUL(x110);
        x112=x15*x53; POST_CHECK_MUL(x112);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x185=x122*x15*x18*mu*rn; POST_CHECK_MUL(x185);
        x186=-x185*x25 + x185*x27; POST_CHECK_ADD(x186);

        /* Assignment result[1, 6]=-x109*x112 + x110*x112 + x186 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x108=x104*x18; POST_CHECK_MUL(x108);
        x109=x50*(rt2 + x108); POST_CHECK_MUL(x109);
        x110=x56*(rt2 - x108); POST_CHECK_MUL(x110);
        x112=x15*x53; POST_CHECK_MUL(x112);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x185=x122*x15*x18*mu*rn; POST_CHECK_MUL(x185);
        x186=-x185*x25 + x185*x27; POST_CHECK_ADD(x186);
        result[19] = -x109*x112 + x110*x112 + x186;

    }
    else if (x177)
    {
        DEBUG_PRINT("Case (x177) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x108=x104*x18; POST_CHECK_MUL(x108);
        x109=x50*(rt2 + x108); POST_CHECK_MUL(x109);
        x110=x56*(rt2 - x108); POST_CHECK_MUL(x110);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);

        /* Assignment result[1, 6]=-x109*x114 + x110*x114 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x108=x104*x18; POST_CHECK_MUL(x108);
        x109=x50*(rt2 + x108); POST_CHECK_MUL(x109);
        x110=x56*(rt2 - x108); POST_CHECK_MUL(x110);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x114=random1*x113; POST_CHECK_MUL(x114);
        result[19] = -x109*x114 + x110*x114;

    }


    /* Assignment result[2, 6]=Piecewise((x181, x39), (1 - x109*x188 + x110*x188 + x202*x25 - x202*x27, x100), (1 - x109*x189 + x110*x189, x177)) */
    double x202;
    if (x39)
    {
        DEBUG_PRINT("Case (x39) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x132=x28*x8; POST_CHECK_MUL(x132);
        x178=9.8994949366116653521885382360778748989105224609375*x28; POST_CHECK_MUL(x178);
        x179=4.0*x130; POST_CHECK_MUL(x179);
        x180=1.4142135623730951454746218587388284504413604736328125*x29; POST_CHECK_MUL(x180);
        x181=x127*(x126 - 15.5563491861040450459086059709079563617706298828125*x132 - x178 - x179 - x180*x8 + 9.8994949366116653521885382360778748989105224609375*x29 - 20.0*x90); POST_CHECK_MUL(x181);

        /* Assignment result[2, 6]=x181 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8) >= 0.;*/
        x28=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x8); POST_CHECK_POW(x28);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8) >= 0.;*/
        x29=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x8); POST_CHECK_POW(x29);
        x36=x28*x29; POST_CHECK_MUL(x36);
        x90=x28*mu; POST_CHECK_MUL(x90);
        x91=x29*x90; POST_CHECK_MUL(x91);
        x126=16.0*x36 + 11.3137084989847611637969748699106276035308837890625*x91; POST_CHECK_ADD(x126);
        /*@ assert (x126) != 0.;*/
        x127=1.0/x126; POST_CHECK_POW(x127);
        x130=x29*mu; POST_CHECK_MUL(x130);
        x132=x28*x8; POST_CHECK_MUL(x132);
        x178=9.8994949366116653521885382360778748989105224609375*x28; POST_CHECK_MUL(x178);
        x179=4.0*x130; POST_CHECK_MUL(x179);
        x180=1.4142135623730951454746218587388284504413604736328125*x29; POST_CHECK_MUL(x180);
        x181=x127*(x126 - 15.5563491861040450459086059709079563617706298828125*x132 - x178 - x179 - x180*x8 + 9.8994949366116653521885382360778748989105224609375*x29 - 20.0*x90); POST_CHECK_MUL(x181);
        result[20] = x181;

    }
    else if (x100)
    {
        DEBUG_PRINT("Case (x100) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x108=x104*x18; POST_CHECK_MUL(x108);
        x109=x50*(rt2 + x108); POST_CHECK_MUL(x109);
        x110=x56*(rt2 - x108); POST_CHECK_MUL(x110);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x182=x122*mu*rn; POST_CHECK_MUL(x182);
        x188=x18*x53; POST_CHECK_MUL(x188);
        x202=x104 - x182*x19; POST_CHECK_ADD(x202);

        /* Assignment result[2, 6]=1 - x109*x188 + x110*x188 + x202*x25 - x202*x27 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x108=x104*x18; POST_CHECK_MUL(x108);
        x109=x50*(rt2 + x108); POST_CHECK_MUL(x109);
        x110=x56*(rt2 - x108); POST_CHECK_MUL(x110);
        /*@ assert (x20) >= 0.;*/
        /*@ assert (x20) != 0.;*/
        x122=pow(x20, -3.0/2.0); POST_CHECK_POW(x122);
        x182=x122*mu*rn; POST_CHECK_MUL(x182);
        x188=x18*x53; POST_CHECK_MUL(x188);
        x202=x104 - x182*x19; POST_CHECK_ADD(x202);
        result[20] = 1 - x109*x188 + x110*x188 + x202*x25 - x202*x27;

    }
    else if (x177)
    {
        DEBUG_PRINT("Case (x177) is True.\n");
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x108=x104*x18; POST_CHECK_MUL(x108);
        x109=x50*(rt2 + x108); POST_CHECK_MUL(x109);
        x110=x56*(rt2 - x108); POST_CHECK_MUL(x110);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x189=random2*x113; POST_CHECK_MUL(x189);

        /* Assignment result[2, 6]=1 - x109*x189 + x110*x189 */
        /*@ assert (x26) != 0.;*/
        x50=0.5/x26; POST_CHECK_MUL(x50);
        /*@ assert (x21) != 0.;*/
        x53=1.0/x21; POST_CHECK_POW(x53);
        /*@ assert (x24) != 0.;*/
        x56=0.5/x24; POST_CHECK_MUL(x56);
        x104=x5*x53; POST_CHECK_MUL(x104);
        x108=x104*x18; POST_CHECK_MUL(x108);
        x109=x50*(rt2 + x108); POST_CHECK_MUL(x109);
        x110=x56*(rt2 - x108); POST_CHECK_MUL(x110);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x113=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x113);
        x189=random2*x113; POST_CHECK_MUL(x189);
        result[20] = 1 - x109*x189 + x110*x189;

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
    double x4;
    double x7;
    double x8;
    double x9;
    double x10;
    double x11;
    int x18;
    double x5;
    double x6;
    double x12;
    double x13;
    double x14;
    double x15;
    double x16;
    double x17;
    x1=ut1*ut1; POST_CHECK_POW(x1);
    /*@ assert (x1) >= 0.;*/
    x2=ut2*ut2; POST_CHECK_POW(x2);
    /*@ assert (x2) >= 0.;*/
    /*@ assert (x1 + x2) >= 0.;*/
    x3=mu*sqrt(x1 + x2) + un; POST_CHECK_ADD(x3);
    x4=mu*rn; POST_CHECK_MUL(x4);
    /*@ assert (x4) >= 0.;*/
    /*@ assert (x4) != 0.;*/
    x7=mu*ut1; POST_CHECK_MUL(x7);
    x8=rt1*x4 + x7*x3; POST_CHECK_ADD(x8);
    x9=mu*ut2; POST_CHECK_MUL(x9);
    x10=rt2*x4 + x9*x3; POST_CHECK_ADD(x10);
    /*@ assert (x10*x10 + x8*x8) >= 0.;*/
    x11=sqrt(x10*x10 + x8*x8); POST_CHECK_POW(x11);
    x18=x11 > 0; POST_CHECK(x18);
    int x21;
    double x19;
    double x20;
    x21=x11 <= ZERO; POST_CHECK(x21);
    if (x18)
    {
        x5=mu*mu; POST_CHECK_POW(x5);
        x6=rn*rn*x5 + rt1*rt1 + rt2*rt2 + x1*x5 + x2*x5 + x3*x3; POST_CHECK_ADD(x6);
        x12=2*x11; POST_CHECK_MUL(x12);
        /*@ assert (-x12 + x6) >= 0.;*/
        x13=0.5*sqrt(-x12 + x6); POST_CHECK_MUL(x13);
        /*@ assert (x12 + x6) >= 0.;*/
        x14=0.5*sqrt(x12 + x6); POST_CHECK_MUL(x14);
        x15=rt1 + x7; POST_CHECK_ADD(x15);
        /*@ assert (x11) != 0.;*/
        x16=1.0/x11; POST_CHECK_POW(x16);
        x17=x8*x16; POST_CHECK_MUL(x17);
    }
    else if (x21)
    {
        x5=mu*mu; POST_CHECK_POW(x5);
        x6=rn*rn*x5 + rt1*rt1 + rt2*rt2 + x1*x5 + x2*x5 + x3*x3; POST_CHECK_ADD(x6);
        x12=2*x11; POST_CHECK_MUL(x12);
        /*@ assert (-x12 + x6) >= 0.;*/
        x13=0.5*sqrt(-x12 + x6); POST_CHECK_MUL(x13);
        /*@ assert (x12 + x6) >= 0.;*/
        x14=0.5*sqrt(x12 + x6); POST_CHECK_MUL(x14);
        x15=rt1 + x7; POST_CHECK_ADD(x15);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x19=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x19);
        x20=random1*x19; POST_CHECK_MUL(x20);
    }
    /* Assignment result[0, 0]=x4 - x13 - x14 + x3 */
    x5=mu*mu; POST_CHECK_POW(x5);
    x6=rn*rn*x5 + rt1*rt1 + rt2*rt2 + x1*x5 + x2*x5 + x3*x3; POST_CHECK_ADD(x6);
    x12=2*x11; POST_CHECK_MUL(x12);
    /*@ assert (-x12 + x6) >= 0.;*/
    x13=0.5*sqrt(-x12 + x6); POST_CHECK_MUL(x13);
    /*@ assert (x12 + x6) >= 0.;*/
    x14=0.5*sqrt(x12 + x6); POST_CHECK_MUL(x14);
    result[0] = x4 - x13 - x14 + x3;


    /* Assignment result[1, 0]=Piecewise((x15 + x17*x13 - x17*x14, x18), (x15 + x20*x13 - x20*x14, x21)) */

    if (x18)
    {
        DEBUG_PRINT("Case (x18) is True.\n");
        x15=rt1 + x7; POST_CHECK_ADD(x15);
        /*@ assert (x11) != 0.;*/
        x16=1.0/x11; POST_CHECK_POW(x16);
        x17=x8*x16; POST_CHECK_MUL(x17);

        /* Assignment result[1, 0]=x15 + x17*x13 - x17*x14 */
        x15=rt1 + x7; POST_CHECK_ADD(x15);
        /*@ assert (x11) != 0.;*/
        x16=1.0/x11; POST_CHECK_POW(x16);
        x17=x8*x16; POST_CHECK_MUL(x17);
        result[1] = x15 + x17*x13 - x17*x14;

    }
    else if (x21)
    {
        DEBUG_PRINT("Case (x21) is True.\n");
        x15=rt1 + x7; POST_CHECK_ADD(x15);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x19=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x19);
        x20=random1*x19; POST_CHECK_MUL(x20);

        /* Assignment result[1, 0]=x15 + x20*x13 - x20*x14 */
        x15=rt1 + x7; POST_CHECK_ADD(x15);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x19=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x19);
        x20=random1*x19; POST_CHECK_MUL(x20);
        result[1] = x15 + x20*x13 - x20*x14;

    }
    /*@ assert (result[1]) >= 0.;*/

    /* Assignment result[2, 0]=Piecewise((x22 + x23*x13 - x23*x14, x18), (x22 + x24*x13 - x24*x14, x21)) */
    double x22;
    double x23;
    double x24;
    if (x18)
    {
        DEBUG_PRINT("Case (x18) is True.\n");
        /*@ assert (x11) != 0.;*/
        x16=1.0/x11; POST_CHECK_POW(x16);
        x22=rt2 + x9; POST_CHECK_ADD(x22);
        x23=x10*x16; POST_CHECK_MUL(x23);

        /* Assignment result[2, 0]=x22 + x23*x13 - x23*x14 */
        /*@ assert (x11) != 0.;*/
        x16=1.0/x11; POST_CHECK_POW(x16);
        x22=rt2 + x9; POST_CHECK_ADD(x22);
        x23=x10*x16; POST_CHECK_MUL(x23);
        result[2] = x22 + x23*x13 - x23*x14;

    }
    else if (x21)
    {
        DEBUG_PRINT("Case (x21) is True.\n");
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x19=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x19);
        x22=rt2 + x9; POST_CHECK_ADD(x22);
        x24=random2*x19; POST_CHECK_MUL(x24);

        /* Assignment result[2, 0]=x22 + x24*x13 - x24*x14 */
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x19=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x19);
        x22=rt2 + x9; POST_CHECK_ADD(x22);
        x24=random2*x19; POST_CHECK_MUL(x24);
        result[2] = x22 + x24*x13 - x24*x14;

    }
    /*@ assert (result[2]) >= 0.;*/
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
    double x11;
    double x12;
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
    int x23;
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
    int x34;
    double x2;
    double x3;
    double x4;
    double x5;
    double x6;
    double x7;
    double x8;
    double x9;
    double x10;
    x1=mu*mu; POST_CHECK_POW(x1);
    x11=mu*rn; POST_CHECK_MUL(x11);
    x12=ut1*ut1; POST_CHECK_POW(x12);
    /*@ assert (x12) >= 0.;*/
    x13=ut2*ut2; POST_CHECK_POW(x13);
    /*@ assert (x13) >= 0.;*/
    /*@ assert (x12 + x13) >= 0.;*/
    x14=sqrt(x12 + x13); POST_CHECK_POW(x14);
    x15=mu*x14 + un; POST_CHECK_ADD(x15);
    x16=mu*x15; POST_CHECK_MUL(x16);
    x17=rt1*x11 + ut1*x16; POST_CHECK_ADD(x17);
    x18=x17*x17; POST_CHECK_POW(x18);
    x19=rt2*x11 + ut2*x16; POST_CHECK_ADD(x19);
    x20=x19*x19; POST_CHECK_POW(x20);
    x21=x18 + x20; POST_CHECK_ADD(x21);
    /*@ assert (x21) >= 0.;*/
    x22=sqrt(x21); POST_CHECK_POW(x22);
    x23=x22 <= ZERO; POST_CHECK(x23);
    x24=rt1*rt1; POST_CHECK_POW(x24);
    /*@ assert (x24) >= 0.;*/
    x25=rt2*rt2; POST_CHECK_POW(x25);
    /*@ assert (x25) >= 0.;*/
    x26=rn*rn; POST_CHECK_POW(x26);
    /*@ assert (x26) != 0.;*/
    x27=x1*x26; POST_CHECK_MUL(x27);
    x28=x1*x12; POST_CHECK_MUL(x28);
    x29=x1*x13; POST_CHECK_MUL(x29);
    x30=x15*x15 + x24 + x25 + x27 + x28 + x29; POST_CHECK_ADD(x30);
    x31=2*x22; POST_CHECK_MUL(x31);
    x32=x30 - x31; POST_CHECK_ADD(x32);
    x33=fabs(x32); POST_CHECK(x33);
    /*@ assert (x33) >= 0.;*/
    x34=x23 || x33 <= ZERO; POST_CHECK(x34);
    int x44;
    double x35;
    double x36;
    double x37;
    double x38;
    double x39;
    double x40;
    double x41;
    double x42;
    double x43;
    x44=x14 <= ZERO; POST_CHECK(x44);
    int x57;
    int x58;
    int x59;
    int x60;
    double x45;
    double x46;
    double x47;
    double x48;
    double x49;
    double x50;
    double x51;
    double x52;
    double x53;
    double x54;
    double x55;
    double x56;
    x57=x14 > 0; POST_CHECK(x57);
    x58=x22 > 0; POST_CHECK(x58);
    x59=x33 > 0; POST_CHECK(x59);
    x60=x57 && x58 && x59; POST_CHECK(x60);
    int x97;
    double x91;
    double x92;
    double x93;
    double x94;
    double x95;
    double x96;
    x97=x58 && x59; POST_CHECK(x97);
    int x120;
    double x118;
    double x119;
    x120=x23 && x57 && x58 && x59; POST_CHECK(x120);
    int x172;
    x172=x23 && x58 && x59; POST_CHECK(x172);
    if (x34)
    {
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x6=0.5*x2; POST_CHECK_MUL(x6);
        x7=0.5*x3; POST_CHECK_MUL(x7);
        x8=-x6 - x7; POST_CHECK_ADD(x8);
        x9=1.4142135623730951454746218587388284504413604736328125*x2; POST_CHECK_MUL(x9);
        x10=x2*x3; POST_CHECK_MUL(x10);
    }
    else if (x44)
    {
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (x37) != 0.;*/
        x38=1.0/x37; POST_CHECK_POW(x38);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        /*@ assert (x39) != 0.;*/
        x40=1.0/x39; POST_CHECK_POW(x40);
        x41=0.5*un; POST_CHECK_MUL(x41);
        x42=x41*x37; POST_CHECK_MUL(x42);
        x43=x39*x37; POST_CHECK_MUL(x43);
    }
    else if (x60)
    {
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x48=ut1*x47; POST_CHECK_MUL(x48);
        x49=mu*x19; POST_CHECK_MUL(x49);
        x50=ut2*x49; POST_CHECK_MUL(x50);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        x52=(x48 + x50)*x51; POST_CHECK_MUL(x52);
        x53=x46*(x15 + x52); POST_CHECK_MUL(x53);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x56=x55*(x15 - x52); POST_CHECK_MUL(x56);
    }
    else if (x97)
    {
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x49=mu*x19; POST_CHECK_MUL(x49);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x91=x1*rn; POST_CHECK_MUL(x91);
        x92=rt1*x47; POST_CHECK_MUL(x92);
        x93=rt2*x49; POST_CHECK_MUL(x93);
        x94=(x92 + x93)*x51; POST_CHECK_MUL(x94);
        x95=x46*(x91 + x94); POST_CHECK_MUL(x95);
        x96=x55*(x91 - x94); POST_CHECK_MUL(x96);
    }
    else if (x120)
    {
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x48=ut1*x47; POST_CHECK_MUL(x48);
        x49=mu*x19; POST_CHECK_MUL(x49);
        x50=ut2*x49; POST_CHECK_MUL(x50);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        x52=(x48 + x50)*x51; POST_CHECK_MUL(x52);
        x53=x46*(x15 + x52); POST_CHECK_MUL(x53);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x56=x55*(x15 - x52); POST_CHECK_MUL(x56);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x119=random1*x118; POST_CHECK_MUL(x119);
    }
    else if (x172)
    {
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x49=mu*x19; POST_CHECK_MUL(x49);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x91=x1*rn; POST_CHECK_MUL(x91);
        x92=rt1*x47; POST_CHECK_MUL(x92);
        x93=rt2*x49; POST_CHECK_MUL(x93);
        x94=(x92 + x93)*x51; POST_CHECK_MUL(x94);
        x95=x46*(x91 + x94); POST_CHECK_MUL(x95);
        x96=x55*(x91 - x94); POST_CHECK_MUL(x96);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x119=random1*x118; POST_CHECK_MUL(x119);
    }
    /* Assignment result[0, 0]=Piecewise((x5*(-mu*x9 + x10 + x8), x34), (x38*x40*(-x39*x41 - x42 + x43), x44), (1 - x53 - x56, x60)) */

    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x6=0.5*x2; POST_CHECK_MUL(x6);
        x7=0.5*x3; POST_CHECK_MUL(x7);
        x8=-x6 - x7; POST_CHECK_ADD(x8);
        x9=1.4142135623730951454746218587388284504413604736328125*x2; POST_CHECK_MUL(x9);
        x10=x2*x3; POST_CHECK_MUL(x10);

        /* Assignment result[0, 0]=x5*(-mu*x9 + x10 + x8) */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x6=0.5*x2; POST_CHECK_MUL(x6);
        x7=0.5*x3; POST_CHECK_MUL(x7);
        x8=-x6 - x7; POST_CHECK_ADD(x8);
        x9=1.4142135623730951454746218587388284504413604736328125*x2; POST_CHECK_MUL(x9);
        x10=x2*x3; POST_CHECK_MUL(x10);
        result[0] = x5*(-mu*x9 + x10 + x8);

    }
    else if (x44)
    {
        DEBUG_PRINT("Case (x44) is True.\n");
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (x37) != 0.;*/
        x38=1.0/x37; POST_CHECK_POW(x38);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        /*@ assert (x39) != 0.;*/
        x40=1.0/x39; POST_CHECK_POW(x40);
        x41=0.5*un; POST_CHECK_MUL(x41);
        x42=x41*x37; POST_CHECK_MUL(x42);
        x43=x39*x37; POST_CHECK_MUL(x43);

        /* Assignment result[0, 0]=x38*x40*(-x39*x41 - x42 + x43) */
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (x37) != 0.;*/
        x38=1.0/x37; POST_CHECK_POW(x38);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        /*@ assert (x39) != 0.;*/
        x40=1.0/x39; POST_CHECK_POW(x40);
        x41=0.5*un; POST_CHECK_MUL(x41);
        x42=x41*x37; POST_CHECK_MUL(x42);
        x43=x39*x37; POST_CHECK_MUL(x43);
        result[0] = x38*x40*(-x39*x41 - x42 + x43);

    }
    else if (x60)
    {
        DEBUG_PRINT("Case (x60) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x48=ut1*x47; POST_CHECK_MUL(x48);
        x49=mu*x19; POST_CHECK_MUL(x49);
        x50=ut2*x49; POST_CHECK_MUL(x50);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        x52=(x48 + x50)*x51; POST_CHECK_MUL(x52);
        x53=x46*(x15 + x52); POST_CHECK_MUL(x53);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x56=x55*(x15 - x52); POST_CHECK_MUL(x56);

        /* Assignment result[0, 0]=1 - x53 - x56 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x48=ut1*x47; POST_CHECK_MUL(x48);
        x49=mu*x19; POST_CHECK_MUL(x49);
        x50=ut2*x49; POST_CHECK_MUL(x50);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        x52=(x48 + x50)*x51; POST_CHECK_MUL(x52);
        x53=x46*(x15 + x52); POST_CHECK_MUL(x53);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x56=x55*(x15 - x52); POST_CHECK_MUL(x56);
        result[0] = 1 - x53 - x56;

    }
    /*@ assert (result[0]) >= 0.;*/

    /* Assignment result[1, 0]=Piecewise((x108, x34), (x110*(0.5*x39*rt1*un*x111*x112*x25 - rt1*x111*x112*x25*x42 + 0.5*x39*un*x111*x112*x113 - x111*x112*x113*x42), x44), (-x117*x53 + x117*x56 - x45*x116 + x54*x116, x60), (-x119*x53 + x119*x56, x120)) */
    double x64;
    double x89;
    double x98;
    double x99;
    double x108;
    double x109;
    double x110;
    double x111;
    double x112;
    double x113;
    double x114;
    double x115;
    double x116;
    double x117;
    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x98=0.353553390593273786368655464684707112610340118408203125*x2; POST_CHECK_MUL(x98);
        x99=0.353553390593273786368655464684707112610340118408203125*x3; POST_CHECK_MUL(x99);
        x108=x5*(-1.0*x89 - x98 + x99); POST_CHECK_MUL(x108);

        /* Assignment result[1, 0]=x108 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x98=0.353553390593273786368655464684707112610340118408203125*x2; POST_CHECK_MUL(x98);
        x99=0.353553390593273786368655464684707112610340118408203125*x3; POST_CHECK_MUL(x99);
        x108=x5*(-1.0*x89 - x98 + x99); POST_CHECK_MUL(x108);
        result[1] = x108;

    }
    else if (x44)
    {
        DEBUG_PRINT("Case (x44) is True.\n");
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (x37) != 0.;*/
        x38=1.0/x37; POST_CHECK_POW(x38);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        /*@ assert (x39) != 0.;*/
        x40=1.0/x39; POST_CHECK_POW(x40);
        x41=0.5*un; POST_CHECK_MUL(x41);
        x42=x41*x37; POST_CHECK_MUL(x42);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x109=pow(x64, 3.0/2.0); POST_CHECK_POW(x109);
        /*@ assert (x109) != 0.;*/
        x110=x38*1.0/x109*x40; POST_CHECK_MUL(x110);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x112=pow(rn, 3); POST_CHECK_POW(x112);
        /*@ assert (x112) >= 0.;*/
        /*@ assert (x112) != 0.;*/
        x113=pow(rt1, 3); POST_CHECK_POW(x113);

        /* Assignment result[1, 0]=x110*(0.5*x39*rt1*un*x111*x112*x25 - rt1*x111*x112*x25*x42 + 0.5*x39*un*x111*x112*x113 - x111*x112*x113*x42) */
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (x37) != 0.;*/
        x38=1.0/x37; POST_CHECK_POW(x38);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        /*@ assert (x39) != 0.;*/
        x40=1.0/x39; POST_CHECK_POW(x40);
        x41=0.5*un; POST_CHECK_MUL(x41);
        x42=x41*x37; POST_CHECK_MUL(x42);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x109=pow(x64, 3.0/2.0); POST_CHECK_POW(x109);
        /*@ assert (x109) != 0.;*/
        x110=x38*1.0/x109*x40; POST_CHECK_MUL(x110);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x112=pow(rn, 3); POST_CHECK_POW(x112);
        /*@ assert (x112) >= 0.;*/
        /*@ assert (x112) != 0.;*/
        x113=pow(rt1, 3); POST_CHECK_POW(x113);
        result[1] = x110*(0.5*x39*rt1*un*x111*x112*x25 - rt1*x111*x112*x25*x42 + 0.5*x39*un*x111*x112*x113 - x111*x112*x113*x42);

    }
    else if (x60)
    {
        DEBUG_PRINT("Case (x60) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x48=ut1*x47; POST_CHECK_MUL(x48);
        x49=mu*x19; POST_CHECK_MUL(x49);
        x50=ut2*x49; POST_CHECK_MUL(x50);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        x52=(x48 + x50)*x51; POST_CHECK_MUL(x52);
        x53=x46*(x15 + x52); POST_CHECK_MUL(x53);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x56=x55*(x15 - x52); POST_CHECK_MUL(x56);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x115=-x48 - x50; POST_CHECK_ADD(x115);
        x116=0.5*mu*ut1*x51 + 0.5*x17*x114*x115; POST_CHECK_ADD(x116);
        x117=x17*x51; POST_CHECK_MUL(x117);

        /* Assignment result[1, 0]=-x117*x53 + x117*x56 - x45*x116 + x54*x116 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x48=ut1*x47; POST_CHECK_MUL(x48);
        x49=mu*x19; POST_CHECK_MUL(x49);
        x50=ut2*x49; POST_CHECK_MUL(x50);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        x52=(x48 + x50)*x51; POST_CHECK_MUL(x52);
        x53=x46*(x15 + x52); POST_CHECK_MUL(x53);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x56=x55*(x15 - x52); POST_CHECK_MUL(x56);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x115=-x48 - x50; POST_CHECK_ADD(x115);
        x116=0.5*mu*ut1*x51 + 0.5*x17*x114*x115; POST_CHECK_ADD(x116);
        x117=x17*x51; POST_CHECK_MUL(x117);
        result[1] = -x117*x53 + x117*x56 - x45*x116 + x54*x116;

    }
    else if (x120)
    {
        DEBUG_PRINT("Case (x120) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x48=ut1*x47; POST_CHECK_MUL(x48);
        x49=mu*x19; POST_CHECK_MUL(x49);
        x50=ut2*x49; POST_CHECK_MUL(x50);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        x52=(x48 + x50)*x51; POST_CHECK_MUL(x52);
        x53=x46*(x15 + x52); POST_CHECK_MUL(x53);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x56=x55*(x15 - x52); POST_CHECK_MUL(x56);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x119=random1*x118; POST_CHECK_MUL(x119);

        /* Assignment result[1, 0]=-x119*x53 + x119*x56 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x48=ut1*x47; POST_CHECK_MUL(x48);
        x49=mu*x19; POST_CHECK_MUL(x49);
        x50=ut2*x49; POST_CHECK_MUL(x50);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        x52=(x48 + x50)*x51; POST_CHECK_MUL(x52);
        x53=x46*(x15 + x52); POST_CHECK_MUL(x53);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x56=x55*(x15 - x52); POST_CHECK_MUL(x56);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x119=random1*x118; POST_CHECK_MUL(x119);
        result[1] = -x119*x53 + x119*x56;

    }
    /*@ assert (result[1]) >= 0.;*/

    /* Assignment result[2, 0]=Piecewise((x108, x34), (x110*(0.5*x39*rt2*un*x111*x112*x24 - rt2*x111*x112*x24*x42 + 0.5*x39*un*x111*x112*x164 - x111*x112*x164*x42), x44), (-x183*x53 + x183*x56 - x45*x182 + x54*x182, x60), (-x184*x53 + x184*x56, x120)) */
    double x164;
    double x182;
    double x183;
    double x184;
    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x98=0.353553390593273786368655464684707112610340118408203125*x2; POST_CHECK_MUL(x98);
        x99=0.353553390593273786368655464684707112610340118408203125*x3; POST_CHECK_MUL(x99);
        x108=x5*(-1.0*x89 - x98 + x99); POST_CHECK_MUL(x108);

        /* Assignment result[2, 0]=x108 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x98=0.353553390593273786368655464684707112610340118408203125*x2; POST_CHECK_MUL(x98);
        x99=0.353553390593273786368655464684707112610340118408203125*x3; POST_CHECK_MUL(x99);
        x108=x5*(-1.0*x89 - x98 + x99); POST_CHECK_MUL(x108);
        result[2] = x108;

    }
    else if (x44)
    {
        DEBUG_PRINT("Case (x44) is True.\n");
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (x37) != 0.;*/
        x38=1.0/x37; POST_CHECK_POW(x38);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        /*@ assert (x39) != 0.;*/
        x40=1.0/x39; POST_CHECK_POW(x40);
        x41=0.5*un; POST_CHECK_MUL(x41);
        x42=x41*x37; POST_CHECK_MUL(x42);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x109=pow(x64, 3.0/2.0); POST_CHECK_POW(x109);
        /*@ assert (x109) != 0.;*/
        x110=x38*1.0/x109*x40; POST_CHECK_MUL(x110);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x112=pow(rn, 3); POST_CHECK_POW(x112);
        /*@ assert (x112) >= 0.;*/
        /*@ assert (x112) != 0.;*/
        x164=pow(rt2, 3); POST_CHECK_POW(x164);

        /* Assignment result[2, 0]=x110*(0.5*x39*rt2*un*x111*x112*x24 - rt2*x111*x112*x24*x42 + 0.5*x39*un*x111*x112*x164 - x111*x112*x164*x42) */
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (x37) != 0.;*/
        x38=1.0/x37; POST_CHECK_POW(x38);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        /*@ assert (x39) != 0.;*/
        x40=1.0/x39; POST_CHECK_POW(x40);
        x41=0.5*un; POST_CHECK_MUL(x41);
        x42=x41*x37; POST_CHECK_MUL(x42);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x109=pow(x64, 3.0/2.0); POST_CHECK_POW(x109);
        /*@ assert (x109) != 0.;*/
        x110=x38*1.0/x109*x40; POST_CHECK_MUL(x110);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x112=pow(rn, 3); POST_CHECK_POW(x112);
        /*@ assert (x112) >= 0.;*/
        /*@ assert (x112) != 0.;*/
        x164=pow(rt2, 3); POST_CHECK_POW(x164);
        result[2] = x110*(0.5*x39*rt2*un*x111*x112*x24 - rt2*x111*x112*x24*x42 + 0.5*x39*un*x111*x112*x164 - x111*x112*x164*x42);

    }
    else if (x60)
    {
        DEBUG_PRINT("Case (x60) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x48=ut1*x47; POST_CHECK_MUL(x48);
        x49=mu*x19; POST_CHECK_MUL(x49);
        x50=ut2*x49; POST_CHECK_MUL(x50);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        x52=(x48 + x50)*x51; POST_CHECK_MUL(x52);
        x53=x46*(x15 + x52); POST_CHECK_MUL(x53);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x56=x55*(x15 - x52); POST_CHECK_MUL(x56);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x115=-x48 - x50; POST_CHECK_ADD(x115);
        x182=0.5*mu*ut2*x51 + 0.5*x19*x114*x115; POST_CHECK_ADD(x182);
        x183=x19*x51; POST_CHECK_MUL(x183);

        /* Assignment result[2, 0]=-x183*x53 + x183*x56 - x45*x182 + x54*x182 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x48=ut1*x47; POST_CHECK_MUL(x48);
        x49=mu*x19; POST_CHECK_MUL(x49);
        x50=ut2*x49; POST_CHECK_MUL(x50);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        x52=(x48 + x50)*x51; POST_CHECK_MUL(x52);
        x53=x46*(x15 + x52); POST_CHECK_MUL(x53);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x56=x55*(x15 - x52); POST_CHECK_MUL(x56);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x115=-x48 - x50; POST_CHECK_ADD(x115);
        x182=0.5*mu*ut2*x51 + 0.5*x19*x114*x115; POST_CHECK_ADD(x182);
        x183=x19*x51; POST_CHECK_MUL(x183);
        result[2] = -x183*x53 + x183*x56 - x45*x182 + x54*x182;

    }
    else if (x120)
    {
        DEBUG_PRINT("Case (x120) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x48=ut1*x47; POST_CHECK_MUL(x48);
        x49=mu*x19; POST_CHECK_MUL(x49);
        x50=ut2*x49; POST_CHECK_MUL(x50);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        x52=(x48 + x50)*x51; POST_CHECK_MUL(x52);
        x53=x46*(x15 + x52); POST_CHECK_MUL(x53);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x56=x55*(x15 - x52); POST_CHECK_MUL(x56);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x184=random2*x118; POST_CHECK_MUL(x184);

        /* Assignment result[2, 0]=-x184*x53 + x184*x56 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x48=ut1*x47; POST_CHECK_MUL(x48);
        x49=mu*x19; POST_CHECK_MUL(x49);
        x50=ut2*x49; POST_CHECK_MUL(x50);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        x52=(x48 + x50)*x51; POST_CHECK_MUL(x52);
        x53=x46*(x15 + x52); POST_CHECK_MUL(x53);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x56=x55*(x15 - x52); POST_CHECK_MUL(x56);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x184=random2*x118; POST_CHECK_MUL(x184);
        result[2] = -x184*x53 + x184*x56;

    }
    /*@ assert (result[2]) >= 0.;*/

    /* Assignment result[0, 1]=Piecewise((x63, x34), (x66*(0.5*x39*rt1*x1*rn*un - rt1*x1*rn*x42 + x68), x44), (x71 - x79 - x80, x60)) */
    double x61;
    double x62;
    double x63;
    double x65;
    double x66;
    double x67;
    double x68;
    double x69;
    double x70;
    double x71;
    double x72;
    double x73;
    double x74;
    double x75;
    double x76;
    double x77;
    double x78;
    double x79;
    double x80;
    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        x61=0.70710678118654757273731092936941422522068023681640625*mu; POST_CHECK_MUL(x61);
        x62=x61*x3; POST_CHECK_MUL(x62);
        x63=x4*(-2.0*x1 - x61 + x62); POST_CHECK_MUL(x63);

        /* Assignment result[0, 1]=x63 */
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        x61=0.70710678118654757273731092936941422522068023681640625*mu; POST_CHECK_MUL(x61);
        x62=x61*x3; POST_CHECK_MUL(x62);
        x63=x4*(-2.0*x1 - x61 + x62); POST_CHECK_MUL(x63);
        result[3] = x63;

    }
    else if (x44)
    {
        DEBUG_PRINT("Case (x44) is True.\n");
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (x37) != 0.;*/
        x38=1.0/x37; POST_CHECK_POW(x38);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        /*@ assert (x39) != 0.;*/
        x40=1.0/x39; POST_CHECK_POW(x40);
        x41=0.5*un; POST_CHECK_MUL(x41);
        x42=x41*x37; POST_CHECK_MUL(x42);
        x43=x39*x37; POST_CHECK_MUL(x43);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x65=sqrt(x64); POST_CHECK_POW(x65);
        /*@ assert (x65) != 0.;*/
        x66=x38*1.0/x65*x40; POST_CHECK_MUL(x66);
        x67=0.353553390593273786368655464684707112610340118408203125*mu*un*x65; POST_CHECK_MUL(x67);
        x68=0.70710678118654757273731092936941422522068023681640625*mu*x65*x43 - x39*x67 - x67*x37; POST_CHECK_ADD(x68);

        /* Assignment result[0, 1]=x66*(0.5*x39*rt1*x1*rn*un - rt1*x1*rn*x42 + x68) */
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (x37) != 0.;*/
        x38=1.0/x37; POST_CHECK_POW(x38);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        /*@ assert (x39) != 0.;*/
        x40=1.0/x39; POST_CHECK_POW(x40);
        x41=0.5*un; POST_CHECK_MUL(x41);
        x42=x41*x37; POST_CHECK_MUL(x42);
        x43=x39*x37; POST_CHECK_MUL(x43);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x65=sqrt(x64); POST_CHECK_POW(x65);
        /*@ assert (x65) != 0.;*/
        x66=x38*1.0/x65*x40; POST_CHECK_MUL(x66);
        x67=0.353553390593273786368655464684707112610340118408203125*mu*un*x65; POST_CHECK_MUL(x67);
        x68=0.70710678118654757273731092936941422522068023681640625*mu*x65*x43 - x39*x67 - x67*x37; POST_CHECK_ADD(x68);
        result[3] = x66*(0.5*x39*rt1*x1*rn*un - rt1*x1*rn*x42 + x68);

    }
    else if (x60)
    {
        DEBUG_PRINT("Case (x60) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x71=ut1*x70; POST_CHECK_MUL(x71);
        x72=ut1*x1 + x71*x15; POST_CHECK_ADD(x72);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x74=x73*x19; POST_CHECK_MUL(x74);
        x75=2*x16; POST_CHECK_MUL(x75);
        x76=x28*x69; POST_CHECK_MUL(x76);
        x77=(1.0/2.0)*x17*(x75 + 2*x76); POST_CHECK_MUL(x77);
        x78=x51*(x74 + x77); POST_CHECK_MUL(x78);
        x79=x46*(x72 + x78); POST_CHECK_MUL(x79);
        x80=x55*(x72 - x78); POST_CHECK_MUL(x80);

        /* Assignment result[0, 1]=x71 - x79 - x80 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x71=ut1*x70; POST_CHECK_MUL(x71);
        x72=ut1*x1 + x71*x15; POST_CHECK_ADD(x72);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x74=x73*x19; POST_CHECK_MUL(x74);
        x75=2*x16; POST_CHECK_MUL(x75);
        x76=x28*x69; POST_CHECK_MUL(x76);
        x77=(1.0/2.0)*x17*(x75 + 2*x76); POST_CHECK_MUL(x77);
        x78=x51*(x74 + x77); POST_CHECK_MUL(x78);
        x79=x46*(x72 + x78); POST_CHECK_MUL(x79);
        x80=x55*(x72 - x78); POST_CHECK_MUL(x80);
        result[3] = x71 - x79 - x80;

    }
    /*@ assert (result[3]) >= 0.;*/

    /* Assignment result[1, 1]=Piecewise((x133, x34), (x137*(-2.0*un*x111*x142*x26*x37 + 2.0*x39*un*x111*x26*x142 + 4.0*x39*x147*un*x111*x26*x24 + 2.0*x39*x152*un*x111*x26*x25 - 2.0*un*x138*x139*x147*x37 + 2.0*x39*x147*un*x138*x139 - 2.0*x111*x140*x147*x26*x37 + 2.0*x39*x147*x111*x26*x140 - x131*x153 + x141 + x145 + x146 + x148 + x149 - x150*x151 + x154 + x155 - 2.0*x156 - 2.0*x157 - 4.0*x158 - 4.0*x159), x44), (mu - x117*x79 + x117*x80 - x45*x161 + x54*x161, x60), (mu - x119*x79 + x119*x80, x120)) */
    double x90;
    double x121;
    double x122;
    double x123;
    double x124;
    double x125;
    double x126;
    double x127;
    double x128;
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
    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x123=4.24264068711928477029005080112256109714508056640625*x2; POST_CHECK_MUL(x123);
        x124=4.24264068711928477029005080112256109714508056640625*x3; POST_CHECK_MUL(x124);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x126=2.0*x125; POST_CHECK_MUL(x126);
        x127=x1*x2; POST_CHECK_MUL(x127);
        x128=x2*x111; POST_CHECK_MUL(x128);
        x129=x1*x124; POST_CHECK_MUL(x129);
        x130=-x129; POST_CHECK_MUL(x130);
        x131=2.0*x111; POST_CHECK_MUL(x131);
        x132=x131*x3; POST_CHECK_MUL(x132);
        x133=x122*(-x123 + x124 + x126 - 57.9827560572968963015227927826344966888427734375*x127 + 11.3137084989847611637969748699106276035308837890625*x127*x3 - 34.0*x128 + x130 + x132 - 26.0*x89 + 16.0*x90); POST_CHECK_MUL(x133);

        /* Assignment result[1, 1]=x133 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x123=4.24264068711928477029005080112256109714508056640625*x2; POST_CHECK_MUL(x123);
        x124=4.24264068711928477029005080112256109714508056640625*x3; POST_CHECK_MUL(x124);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x126=2.0*x125; POST_CHECK_MUL(x126);
        x127=x1*x2; POST_CHECK_MUL(x127);
        x128=x2*x111; POST_CHECK_MUL(x128);
        x129=x1*x124; POST_CHECK_MUL(x129);
        x130=-x129; POST_CHECK_MUL(x130);
        x131=2.0*x111; POST_CHECK_MUL(x131);
        x132=x131*x3; POST_CHECK_MUL(x132);
        x133=x122*(-x123 + x124 + x126 - 57.9827560572968963015227927826344966888427734375*x127 + 11.3137084989847611637969748699106276035308837890625*x127*x3 - 34.0*x128 + x130 + x132 - 26.0*x89 + 16.0*x90); POST_CHECK_MUL(x133);
        result[4] = x133;

    }
    else if (x44)
    {
        DEBUG_PRINT("Case (x44) is True.\n");
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x109=pow(x64, 3.0/2.0); POST_CHECK_POW(x109);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x112=pow(rn, 3); POST_CHECK_POW(x112);
        /*@ assert (x112) >= 0.;*/
        /*@ assert (x112) != 0.;*/
        x113=pow(rt1, 3); POST_CHECK_POW(x113);
        x131=2.0*x111; POST_CHECK_MUL(x131);
        x134=4.0*x39*x109*x37; POST_CHECK_MUL(x134);
        x135=x24*x134; POST_CHECK_MUL(x135);
        x136=x25*x134; POST_CHECK_MUL(x136);
        /*@ assert (x135 + x136) != 0.;*/
        x137=1.0/(x135 + x136); POST_CHECK_POW(x137);
        x138=pow(mu, 5); POST_CHECK_POW(x138);
        /*@ assert (x138) >= 0.;*/
        /*@ assert (x138) != 0.;*/
        x139=pow(rn, 4); POST_CHECK_POW(x139);
        /*@ assert (x139) >= 0.;*/
        /*@ assert (x139) != 0.;*/
        x140=pow(un, 3); POST_CHECK_POW(x140);
        /*@ assert (x140) >= 0.;*/
        /*@ assert (x140) != 0.;*/
        x141=mu*x135 + mu*x136 - 2.0*un*x138*x139*x24*x25*x37 + 2.0*x39*un*x138*x139*x24*x25 - 2.0*x111*x140*x24*x25*x26*x37 + 2.0*x39*x111*x26*x24*x25*x140; POST_CHECK_ADD(x141);
        x142=pow(rt2, 6); POST_CHECK_POW(x142);
        /*@ assert (x142) >= 0.;*/
        x143=pow(mu, 4); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x144=pow(rt1, 5); POST_CHECK_POW(x144);
        x145=-1.4142135623730951454746218587388284504413604736328125*un*x112*x143*x144*x37; POST_CHECK_MUL(x145);
        x146=1.4142135623730951454746218587388284504413604736328125*x39*un*x143*x144*x112; POST_CHECK_MUL(x146);
        x147=pow(rt2, 4); POST_CHECK_POW(x147);
        /*@ assert (x147) >= 0.;*/
        x148=-1.4142135623730951454746218587388284504413604736328125*rt1*un*x112*x143*x147*x37; POST_CHECK_MUL(x148);
        x149=1.4142135623730951454746218587388284504413604736328125*x39*rt1*x147*un*x143*x112; POST_CHECK_MUL(x149);
        x150=4.0*x111; POST_CHECK_MUL(x150);
        x151=un*x147*x24*x26*x37; POST_CHECK_MUL(x151);
        x152=pow(rt1, 4); POST_CHECK_POW(x152);
        /*@ assert (x152) >= 0.;*/
        x153=un*x152*x25*x26*x37; POST_CHECK_MUL(x153);
        x154=-2.828427124746190290949243717477656900882720947265625*un*x112*x113*x143*x25*x37; POST_CHECK_MUL(x154);
        x155=2.828427124746190290949243717477656900882720947265625*x39*un*x143*x112*x113*x25; POST_CHECK_MUL(x155);
        x156=mu*un*x24*x109*x37; POST_CHECK_MUL(x156);
        x157=x39*mu*un*x24*x109; POST_CHECK_MUL(x157);
        x158=mu*un*x25*x109*x37; POST_CHECK_MUL(x158);
        x159=x39*mu*un*x25*x109; POST_CHECK_MUL(x159);

        /* Assignment result[1, 1]=x137*(-2.0*un*x111*x142*x26*x37 + 2.0*x39*un*x111*x26*x142 + 4.0*x39*x147*un*x111*x26*x24 + 2.0*x39*x152*un*x111*x26*x25 - 2.0*un*x138*x139*x147*x37 + 2.0*x39*x147*un*x138*x139 - 2.0*x111*x140*x147*x26*x37 + 2.0*x39*x147*x111*x26*x140 - x131*x153 + x141 + x145 + x146 + x148 + x149 - x150*x151 + x154 + x155 - 2.0*x156 - 2.0*x157 - 4.0*x158 - 4.0*x159) */
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x109=pow(x64, 3.0/2.0); POST_CHECK_POW(x109);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x112=pow(rn, 3); POST_CHECK_POW(x112);
        /*@ assert (x112) >= 0.;*/
        /*@ assert (x112) != 0.;*/
        x113=pow(rt1, 3); POST_CHECK_POW(x113);
        x131=2.0*x111; POST_CHECK_MUL(x131);
        x134=4.0*x39*x109*x37; POST_CHECK_MUL(x134);
        x135=x24*x134; POST_CHECK_MUL(x135);
        x136=x25*x134; POST_CHECK_MUL(x136);
        /*@ assert (x135 + x136) != 0.;*/
        x137=1.0/(x135 + x136); POST_CHECK_POW(x137);
        x138=pow(mu, 5); POST_CHECK_POW(x138);
        /*@ assert (x138) >= 0.;*/
        /*@ assert (x138) != 0.;*/
        x139=pow(rn, 4); POST_CHECK_POW(x139);
        /*@ assert (x139) >= 0.;*/
        /*@ assert (x139) != 0.;*/
        x140=pow(un, 3); POST_CHECK_POW(x140);
        /*@ assert (x140) >= 0.;*/
        /*@ assert (x140) != 0.;*/
        x141=mu*x135 + mu*x136 - 2.0*un*x138*x139*x24*x25*x37 + 2.0*x39*un*x138*x139*x24*x25 - 2.0*x111*x140*x24*x25*x26*x37 + 2.0*x39*x111*x26*x24*x25*x140; POST_CHECK_ADD(x141);
        x142=pow(rt2, 6); POST_CHECK_POW(x142);
        /*@ assert (x142) >= 0.;*/
        x143=pow(mu, 4); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x144=pow(rt1, 5); POST_CHECK_POW(x144);
        x145=-1.4142135623730951454746218587388284504413604736328125*un*x112*x143*x144*x37; POST_CHECK_MUL(x145);
        x146=1.4142135623730951454746218587388284504413604736328125*x39*un*x143*x144*x112; POST_CHECK_MUL(x146);
        x147=pow(rt2, 4); POST_CHECK_POW(x147);
        /*@ assert (x147) >= 0.;*/
        x148=-1.4142135623730951454746218587388284504413604736328125*rt1*un*x112*x143*x147*x37; POST_CHECK_MUL(x148);
        x149=1.4142135623730951454746218587388284504413604736328125*x39*rt1*x147*un*x143*x112; POST_CHECK_MUL(x149);
        x150=4.0*x111; POST_CHECK_MUL(x150);
        x151=un*x147*x24*x26*x37; POST_CHECK_MUL(x151);
        x152=pow(rt1, 4); POST_CHECK_POW(x152);
        /*@ assert (x152) >= 0.;*/
        x153=un*x152*x25*x26*x37; POST_CHECK_MUL(x153);
        x154=-2.828427124746190290949243717477656900882720947265625*un*x112*x113*x143*x25*x37; POST_CHECK_MUL(x154);
        x155=2.828427124746190290949243717477656900882720947265625*x39*un*x143*x112*x113*x25; POST_CHECK_MUL(x155);
        x156=mu*un*x24*x109*x37; POST_CHECK_MUL(x156);
        x157=x39*mu*un*x24*x109; POST_CHECK_MUL(x157);
        x158=mu*un*x25*x109*x37; POST_CHECK_MUL(x158);
        x159=x39*mu*un*x25*x109; POST_CHECK_MUL(x159);
        result[4] = x137*(-2.0*un*x111*x142*x26*x37 + 2.0*x39*un*x111*x26*x142 + 4.0*x39*x147*un*x111*x26*x24 + 2.0*x39*x152*un*x111*x26*x25 - 2.0*un*x138*x139*x147*x37 + 2.0*x39*x147*un*x138*x139 - 2.0*x111*x140*x147*x26*x37 + 2.0*x39*x147*x111*x26*x140 - x131*x153 + x141 + x145 + x146 + x148 + x149 - x150*x151 + x154 + x155 - 2.0*x156 - 2.0*x157 - 4.0*x158 - 4.0*x159);

    }
    else if (x60)
    {
        DEBUG_PRINT("Case (x60) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x71=ut1*x70; POST_CHECK_MUL(x71);
        x72=ut1*x1 + x71*x15; POST_CHECK_ADD(x72);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x74=x73*x19; POST_CHECK_MUL(x74);
        x75=2*x16; POST_CHECK_MUL(x75);
        x76=x28*x69; POST_CHECK_MUL(x76);
        x77=(1.0/2.0)*x17*(x75 + 2*x76); POST_CHECK_MUL(x77);
        x78=x51*(x74 + x77); POST_CHECK_MUL(x78);
        x79=x46*(x72 + x78); POST_CHECK_MUL(x79);
        x80=x55*(x72 - x78); POST_CHECK_MUL(x80);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x117=x17*x51; POST_CHECK_MUL(x117);
        x160=-x74 - x77; POST_CHECK_ADD(x160);
        x161=0.5*x17*x114*x160 + 0.5*(x16 + x76)*x51; POST_CHECK_ADD(x161);

        /* Assignment result[1, 1]=mu - x117*x79 + x117*x80 - x45*x161 + x54*x161 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x71=ut1*x70; POST_CHECK_MUL(x71);
        x72=ut1*x1 + x71*x15; POST_CHECK_ADD(x72);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x74=x73*x19; POST_CHECK_MUL(x74);
        x75=2*x16; POST_CHECK_MUL(x75);
        x76=x28*x69; POST_CHECK_MUL(x76);
        x77=(1.0/2.0)*x17*(x75 + 2*x76); POST_CHECK_MUL(x77);
        x78=x51*(x74 + x77); POST_CHECK_MUL(x78);
        x79=x46*(x72 + x78); POST_CHECK_MUL(x79);
        x80=x55*(x72 - x78); POST_CHECK_MUL(x80);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x117=x17*x51; POST_CHECK_MUL(x117);
        x160=-x74 - x77; POST_CHECK_ADD(x160);
        x161=0.5*x17*x114*x160 + 0.5*(x16 + x76)*x51; POST_CHECK_ADD(x161);
        result[4] = mu - x117*x79 + x117*x80 - x45*x161 + x54*x161;

    }
    else if (x120)
    {
        DEBUG_PRINT("Case (x120) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x71=ut1*x70; POST_CHECK_MUL(x71);
        x72=ut1*x1 + x71*x15; POST_CHECK_ADD(x72);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x74=x73*x19; POST_CHECK_MUL(x74);
        x75=2*x16; POST_CHECK_MUL(x75);
        x76=x28*x69; POST_CHECK_MUL(x76);
        x77=(1.0/2.0)*x17*(x75 + 2*x76); POST_CHECK_MUL(x77);
        x78=x51*(x74 + x77); POST_CHECK_MUL(x78);
        x79=x46*(x72 + x78); POST_CHECK_MUL(x79);
        x80=x55*(x72 - x78); POST_CHECK_MUL(x80);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x119=random1*x118; POST_CHECK_MUL(x119);

        /* Assignment result[1, 1]=mu - x119*x79 + x119*x80 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x71=ut1*x70; POST_CHECK_MUL(x71);
        x72=ut1*x1 + x71*x15; POST_CHECK_ADD(x72);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x74=x73*x19; POST_CHECK_MUL(x74);
        x75=2*x16; POST_CHECK_MUL(x75);
        x76=x28*x69; POST_CHECK_MUL(x76);
        x77=(1.0/2.0)*x17*(x75 + 2*x76); POST_CHECK_MUL(x77);
        x78=x51*(x74 + x77); POST_CHECK_MUL(x78);
        x79=x46*(x72 + x78); POST_CHECK_MUL(x79);
        x80=x55*(x72 - x78); POST_CHECK_MUL(x80);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x119=random1*x118; POST_CHECK_MUL(x119);
        result[4] = mu - x119*x79 + x119*x80;

    }
    /*@ assert (result[4]) >= 0.;*/

    /* Assignment result[2, 1]=Piecewise((x162, x34), (x137*(x165 + x185 + x186 + x187 + x188 + x189 + x190), x44), (-x183*x79 + x183*x80 - x45*x191 + x54*x191, x60), (-x184*x79 + x184*x80, x120)) */
    double x162;
    double x163;
    double x165;
    double x166;
    double x185;
    double x186;
    double x187;
    double x188;
    double x189;
    double x190;
    double x191;
    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x9=1.4142135623730951454746218587388284504413604736328125*x2; POST_CHECK_MUL(x9);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x123=4.24264068711928477029005080112256109714508056640625*x2; POST_CHECK_MUL(x123);
        x124=4.24264068711928477029005080112256109714508056640625*x3; POST_CHECK_MUL(x124);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x126=2.0*x125; POST_CHECK_MUL(x126);
        x128=x2*x111; POST_CHECK_MUL(x128);
        x129=x1*x124; POST_CHECK_MUL(x129);
        x131=2.0*x111; POST_CHECK_MUL(x131);
        x132=x131*x3; POST_CHECK_MUL(x132);
        x162=x122*(x1*x9 + x123 - x124 - x126 + 2.0*x128 + x129 - x132 + 10.0*x89); POST_CHECK_MUL(x162);

        /* Assignment result[2, 1]=x162 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x9=1.4142135623730951454746218587388284504413604736328125*x2; POST_CHECK_MUL(x9);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x123=4.24264068711928477029005080112256109714508056640625*x2; POST_CHECK_MUL(x123);
        x124=4.24264068711928477029005080112256109714508056640625*x3; POST_CHECK_MUL(x124);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x126=2.0*x125; POST_CHECK_MUL(x126);
        x128=x2*x111; POST_CHECK_MUL(x128);
        x129=x1*x124; POST_CHECK_MUL(x129);
        x131=2.0*x111; POST_CHECK_MUL(x131);
        x132=x131*x3; POST_CHECK_MUL(x132);
        x162=x122*(x1*x9 + x123 - x124 - x126 + 2.0*x128 + x129 - x132 + 10.0*x89); POST_CHECK_MUL(x162);
        result[5] = x162;

    }
    else if (x44)
    {
        DEBUG_PRINT("Case (x44) is True.\n");
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x109=pow(x64, 3.0/2.0); POST_CHECK_POW(x109);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x112=pow(rn, 3); POST_CHECK_POW(x112);
        /*@ assert (x112) >= 0.;*/
        /*@ assert (x112) != 0.;*/
        x113=pow(rt1, 3); POST_CHECK_POW(x113);
        x134=4.0*x39*x109*x37; POST_CHECK_MUL(x134);
        x135=x24*x134; POST_CHECK_MUL(x135);
        x136=x25*x134; POST_CHECK_MUL(x136);
        /*@ assert (x135 + x136) != 0.;*/
        x137=1.0/(x135 + x136); POST_CHECK_POW(x137);
        x138=pow(mu, 5); POST_CHECK_POW(x138);
        /*@ assert (x138) >= 0.;*/
        /*@ assert (x138) != 0.;*/
        x139=pow(rn, 4); POST_CHECK_POW(x139);
        /*@ assert (x139) >= 0.;*/
        /*@ assert (x139) != 0.;*/
        x140=pow(un, 3); POST_CHECK_POW(x140);
        /*@ assert (x140) >= 0.;*/
        /*@ assert (x140) != 0.;*/
        x143=pow(mu, 4); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x144=pow(rt1, 5); POST_CHECK_POW(x144);
        x152=pow(rt1, 4); POST_CHECK_POW(x152);
        /*@ assert (x152) >= 0.;*/
        x163=pow(rt2, 5); POST_CHECK_POW(x163);
        x164=pow(rt2, 3); POST_CHECK_POW(x164);
        x165=2.0*mu*rt1*rt2*un*x109*x37 + 2.0*x39*rt1*rt2*mu*un*x109 + 2.0*rt1*un*x111*x26*x163*x37 - 2.0*x39*rt1*un*x111*x26*x163 + 2.0*rt1*un*x138*x139*x164*x37 - 2.0*x39*rt1*un*x138*x139*x164 + 2.0*rt1*x111*x26*x164*x140*x37 - 2.0*x39*rt1*x111*x26*x164*x140 + 2.0*rt2*un*x111*x144*x26*x37 - 2.0*x39*rt2*un*x144*x111*x26 + 2.0*rt2*un*x113*x138*x139*x37 - 2.0*x39*rt2*un*x138*x139*x113 + 2.0*rt2*x111*x113*x140*x26*x37 - 2.0*x39*rt2*x111*x26*x113*x140 + 4.0*un*x111*x113*x26*x164*x37 - 4.0*x39*un*x111*x26*x113*x164; POST_CHECK_ADD(x165);
        x185=-1.4142135623730951454746218587388284504413604736328125*un*x112*x143*x163*x37; POST_CHECK_MUL(x185);
        x186=1.4142135623730951454746218587388284504413604736328125*x39*un*x143*x112*x163; POST_CHECK_MUL(x186);
        x187=-1.4142135623730951454746218587388284504413604736328125*rt2*un*x112*x143*x152*x37; POST_CHECK_MUL(x187);
        x188=1.4142135623730951454746218587388284504413604736328125*x39*rt2*x152*un*x143*x112; POST_CHECK_MUL(x188);
        x189=-2.828427124746190290949243717477656900882720947265625*un*x112*x143*x24*x164*x37; POST_CHECK_MUL(x189);
        x190=2.828427124746190290949243717477656900882720947265625*x39*un*x143*x112*x24*x164; POST_CHECK_MUL(x190);

        /* Assignment result[2, 1]=x137*(x165 + x185 + x186 + x187 + x188 + x189 + x190) */
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x109=pow(x64, 3.0/2.0); POST_CHECK_POW(x109);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x112=pow(rn, 3); POST_CHECK_POW(x112);
        /*@ assert (x112) >= 0.;*/
        /*@ assert (x112) != 0.;*/
        x113=pow(rt1, 3); POST_CHECK_POW(x113);
        x134=4.0*x39*x109*x37; POST_CHECK_MUL(x134);
        x135=x24*x134; POST_CHECK_MUL(x135);
        x136=x25*x134; POST_CHECK_MUL(x136);
        /*@ assert (x135 + x136) != 0.;*/
        x137=1.0/(x135 + x136); POST_CHECK_POW(x137);
        x138=pow(mu, 5); POST_CHECK_POW(x138);
        /*@ assert (x138) >= 0.;*/
        /*@ assert (x138) != 0.;*/
        x139=pow(rn, 4); POST_CHECK_POW(x139);
        /*@ assert (x139) >= 0.;*/
        /*@ assert (x139) != 0.;*/
        x140=pow(un, 3); POST_CHECK_POW(x140);
        /*@ assert (x140) >= 0.;*/
        /*@ assert (x140) != 0.;*/
        x143=pow(mu, 4); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x144=pow(rt1, 5); POST_CHECK_POW(x144);
        x152=pow(rt1, 4); POST_CHECK_POW(x152);
        /*@ assert (x152) >= 0.;*/
        x163=pow(rt2, 5); POST_CHECK_POW(x163);
        x164=pow(rt2, 3); POST_CHECK_POW(x164);
        x165=2.0*mu*rt1*rt2*un*x109*x37 + 2.0*x39*rt1*rt2*mu*un*x109 + 2.0*rt1*un*x111*x26*x163*x37 - 2.0*x39*rt1*un*x111*x26*x163 + 2.0*rt1*un*x138*x139*x164*x37 - 2.0*x39*rt1*un*x138*x139*x164 + 2.0*rt1*x111*x26*x164*x140*x37 - 2.0*x39*rt1*x111*x26*x164*x140 + 2.0*rt2*un*x111*x144*x26*x37 - 2.0*x39*rt2*un*x144*x111*x26 + 2.0*rt2*un*x113*x138*x139*x37 - 2.0*x39*rt2*un*x138*x139*x113 + 2.0*rt2*x111*x113*x140*x26*x37 - 2.0*x39*rt2*x111*x26*x113*x140 + 4.0*un*x111*x113*x26*x164*x37 - 4.0*x39*un*x111*x26*x113*x164; POST_CHECK_ADD(x165);
        x185=-1.4142135623730951454746218587388284504413604736328125*un*x112*x143*x163*x37; POST_CHECK_MUL(x185);
        x186=1.4142135623730951454746218587388284504413604736328125*x39*un*x143*x112*x163; POST_CHECK_MUL(x186);
        x187=-1.4142135623730951454746218587388284504413604736328125*rt2*un*x112*x143*x152*x37; POST_CHECK_MUL(x187);
        x188=1.4142135623730951454746218587388284504413604736328125*x39*rt2*x152*un*x143*x112; POST_CHECK_MUL(x188);
        x189=-2.828427124746190290949243717477656900882720947265625*un*x112*x143*x24*x164*x37; POST_CHECK_MUL(x189);
        x190=2.828427124746190290949243717477656900882720947265625*x39*un*x143*x112*x24*x164; POST_CHECK_MUL(x190);
        result[5] = x137*(x165 + x185 + x186 + x187 + x188 + x189 + x190);

    }
    else if (x60)
    {
        DEBUG_PRINT("Case (x60) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x71=ut1*x70; POST_CHECK_MUL(x71);
        x72=ut1*x1 + x71*x15; POST_CHECK_ADD(x72);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x74=x73*x19; POST_CHECK_MUL(x74);
        x75=2*x16; POST_CHECK_MUL(x75);
        x76=x28*x69; POST_CHECK_MUL(x76);
        x77=(1.0/2.0)*x17*(x75 + 2*x76); POST_CHECK_MUL(x77);
        x78=x51*(x74 + x77); POST_CHECK_MUL(x78);
        x79=x46*(x72 + x78); POST_CHECK_MUL(x79);
        x80=x55*(x72 - x78); POST_CHECK_MUL(x80);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x160=-x74 - x77; POST_CHECK_ADD(x160);
        x166=0.5*ut1*ut2*x1*x69*x51; POST_CHECK_MUL(x166);
        x183=x19*x51; POST_CHECK_MUL(x183);
        x191=x166 + 0.5*x19*x114*x160; POST_CHECK_ADD(x191);

        /* Assignment result[2, 1]=-x183*x79 + x183*x80 - x45*x191 + x54*x191 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x71=ut1*x70; POST_CHECK_MUL(x71);
        x72=ut1*x1 + x71*x15; POST_CHECK_ADD(x72);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x74=x73*x19; POST_CHECK_MUL(x74);
        x75=2*x16; POST_CHECK_MUL(x75);
        x76=x28*x69; POST_CHECK_MUL(x76);
        x77=(1.0/2.0)*x17*(x75 + 2*x76); POST_CHECK_MUL(x77);
        x78=x51*(x74 + x77); POST_CHECK_MUL(x78);
        x79=x46*(x72 + x78); POST_CHECK_MUL(x79);
        x80=x55*(x72 - x78); POST_CHECK_MUL(x80);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x160=-x74 - x77; POST_CHECK_ADD(x160);
        x166=0.5*ut1*ut2*x1*x69*x51; POST_CHECK_MUL(x166);
        x183=x19*x51; POST_CHECK_MUL(x183);
        x191=x166 + 0.5*x19*x114*x160; POST_CHECK_ADD(x191);
        result[5] = -x183*x79 + x183*x80 - x45*x191 + x54*x191;

    }
    else if (x120)
    {
        DEBUG_PRINT("Case (x120) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x71=ut1*x70; POST_CHECK_MUL(x71);
        x72=ut1*x1 + x71*x15; POST_CHECK_ADD(x72);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x74=x73*x19; POST_CHECK_MUL(x74);
        x75=2*x16; POST_CHECK_MUL(x75);
        x76=x28*x69; POST_CHECK_MUL(x76);
        x77=(1.0/2.0)*x17*(x75 + 2*x76); POST_CHECK_MUL(x77);
        x78=x51*(x74 + x77); POST_CHECK_MUL(x78);
        x79=x46*(x72 + x78); POST_CHECK_MUL(x79);
        x80=x55*(x72 - x78); POST_CHECK_MUL(x80);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x184=random2*x118; POST_CHECK_MUL(x184);

        /* Assignment result[2, 1]=-x184*x79 + x184*x80 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x71=ut1*x70; POST_CHECK_MUL(x71);
        x72=ut1*x1 + x71*x15; POST_CHECK_ADD(x72);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x74=x73*x19; POST_CHECK_MUL(x74);
        x75=2*x16; POST_CHECK_MUL(x75);
        x76=x28*x69; POST_CHECK_MUL(x76);
        x77=(1.0/2.0)*x17*(x75 + 2*x76); POST_CHECK_MUL(x77);
        x78=x51*(x74 + x77); POST_CHECK_MUL(x78);
        x79=x46*(x72 + x78); POST_CHECK_MUL(x79);
        x80=x55*(x72 - x78); POST_CHECK_MUL(x80);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x184=random2*x118; POST_CHECK_MUL(x184);
        result[5] = -x184*x79 + x184*x80;

    }
    /*@ assert (result[5]) >= 0.;*/

    /* Assignment result[0, 2]=Piecewise((x63, x34), (x66*(0.5*x39*rt2*x1*rn*un - rt2*x1*rn*x42 + x68), x44), (x81 - x87 - x88, x60)) */
    double x81;
    double x82;
    double x83;
    double x84;
    double x85;
    double x86;
    double x87;
    double x88;
    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        x61=0.70710678118654757273731092936941422522068023681640625*mu; POST_CHECK_MUL(x61);
        x62=x61*x3; POST_CHECK_MUL(x62);
        x63=x4*(-2.0*x1 - x61 + x62); POST_CHECK_MUL(x63);

        /* Assignment result[0, 2]=x63 */
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        x61=0.70710678118654757273731092936941422522068023681640625*mu; POST_CHECK_MUL(x61);
        x62=x61*x3; POST_CHECK_MUL(x62);
        x63=x4*(-2.0*x1 - x61 + x62); POST_CHECK_MUL(x63);
        result[6] = x63;

    }
    else if (x44)
    {
        DEBUG_PRINT("Case (x44) is True.\n");
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (x37) != 0.;*/
        x38=1.0/x37; POST_CHECK_POW(x38);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        /*@ assert (x39) != 0.;*/
        x40=1.0/x39; POST_CHECK_POW(x40);
        x41=0.5*un; POST_CHECK_MUL(x41);
        x42=x41*x37; POST_CHECK_MUL(x42);
        x43=x39*x37; POST_CHECK_MUL(x43);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x65=sqrt(x64); POST_CHECK_POW(x65);
        /*@ assert (x65) != 0.;*/
        x66=x38*1.0/x65*x40; POST_CHECK_MUL(x66);
        x67=0.353553390593273786368655464684707112610340118408203125*mu*un*x65; POST_CHECK_MUL(x67);
        x68=0.70710678118654757273731092936941422522068023681640625*mu*x65*x43 - x39*x67 - x67*x37; POST_CHECK_ADD(x68);

        /* Assignment result[0, 2]=x66*(0.5*x39*rt2*x1*rn*un - rt2*x1*rn*x42 + x68) */
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (x37) != 0.;*/
        x38=1.0/x37; POST_CHECK_POW(x38);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        /*@ assert (x39) != 0.;*/
        x40=1.0/x39; POST_CHECK_POW(x40);
        x41=0.5*un; POST_CHECK_MUL(x41);
        x42=x41*x37; POST_CHECK_MUL(x42);
        x43=x39*x37; POST_CHECK_MUL(x43);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x65=sqrt(x64); POST_CHECK_POW(x65);
        /*@ assert (x65) != 0.;*/
        x66=x38*1.0/x65*x40; POST_CHECK_MUL(x66);
        x67=0.353553390593273786368655464684707112610340118408203125*mu*un*x65; POST_CHECK_MUL(x67);
        x68=0.70710678118654757273731092936941422522068023681640625*mu*x65*x43 - x39*x67 - x67*x37; POST_CHECK_ADD(x68);
        result[6] = x66*(0.5*x39*rt2*x1*rn*un - rt2*x1*rn*x42 + x68);

    }
    else if (x60)
    {
        DEBUG_PRINT("Case (x60) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x75=2*x16; POST_CHECK_MUL(x75);
        x81=ut2*x70; POST_CHECK_MUL(x81);
        x82=ut2*x1 + x81*x15; POST_CHECK_ADD(x82);
        x83=x73*x17; POST_CHECK_MUL(x83);
        x84=x29*x69; POST_CHECK_MUL(x84);
        x85=(1.0/2.0)*x19*(x75 + 2*x84); POST_CHECK_MUL(x85);
        x86=x51*(x83 + x85); POST_CHECK_MUL(x86);
        x87=x46*(x82 + x86); POST_CHECK_MUL(x87);
        x88=x55*(x82 - x86); POST_CHECK_MUL(x88);

        /* Assignment result[0, 2]=x81 - x87 - x88 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x75=2*x16; POST_CHECK_MUL(x75);
        x81=ut2*x70; POST_CHECK_MUL(x81);
        x82=ut2*x1 + x81*x15; POST_CHECK_ADD(x82);
        x83=x73*x17; POST_CHECK_MUL(x83);
        x84=x29*x69; POST_CHECK_MUL(x84);
        x85=(1.0/2.0)*x19*(x75 + 2*x84); POST_CHECK_MUL(x85);
        x86=x51*(x83 + x85); POST_CHECK_MUL(x86);
        x87=x46*(x82 + x86); POST_CHECK_MUL(x87);
        x88=x55*(x82 - x86); POST_CHECK_MUL(x88);
        result[6] = x81 - x87 - x88;

    }
    /*@ assert (result[6]) >= 0.;*/

    /* Assignment result[1, 2]=Piecewise((x162, x34), (x137*(x145 + x146 + x148 + x149 + x154 + x155 + x165), x44), (-x117*x87 + x117*x88 - x45*x168 + x54*x168, x60), (-x119*x87 + x119*x88, x120)) */
    double x167;
    double x168;
    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x9=1.4142135623730951454746218587388284504413604736328125*x2; POST_CHECK_MUL(x9);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x123=4.24264068711928477029005080112256109714508056640625*x2; POST_CHECK_MUL(x123);
        x124=4.24264068711928477029005080112256109714508056640625*x3; POST_CHECK_MUL(x124);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x126=2.0*x125; POST_CHECK_MUL(x126);
        x128=x2*x111; POST_CHECK_MUL(x128);
        x129=x1*x124; POST_CHECK_MUL(x129);
        x131=2.0*x111; POST_CHECK_MUL(x131);
        x132=x131*x3; POST_CHECK_MUL(x132);
        x162=x122*(x1*x9 + x123 - x124 - x126 + 2.0*x128 + x129 - x132 + 10.0*x89); POST_CHECK_MUL(x162);

        /* Assignment result[1, 2]=x162 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x9=1.4142135623730951454746218587388284504413604736328125*x2; POST_CHECK_MUL(x9);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x123=4.24264068711928477029005080112256109714508056640625*x2; POST_CHECK_MUL(x123);
        x124=4.24264068711928477029005080112256109714508056640625*x3; POST_CHECK_MUL(x124);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x126=2.0*x125; POST_CHECK_MUL(x126);
        x128=x2*x111; POST_CHECK_MUL(x128);
        x129=x1*x124; POST_CHECK_MUL(x129);
        x131=2.0*x111; POST_CHECK_MUL(x131);
        x132=x131*x3; POST_CHECK_MUL(x132);
        x162=x122*(x1*x9 + x123 - x124 - x126 + 2.0*x128 + x129 - x132 + 10.0*x89); POST_CHECK_MUL(x162);
        result[7] = x162;

    }
    else if (x44)
    {
        DEBUG_PRINT("Case (x44) is True.\n");
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x109=pow(x64, 3.0/2.0); POST_CHECK_POW(x109);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x112=pow(rn, 3); POST_CHECK_POW(x112);
        /*@ assert (x112) >= 0.;*/
        /*@ assert (x112) != 0.;*/
        x113=pow(rt1, 3); POST_CHECK_POW(x113);
        x134=4.0*x39*x109*x37; POST_CHECK_MUL(x134);
        x135=x24*x134; POST_CHECK_MUL(x135);
        x136=x25*x134; POST_CHECK_MUL(x136);
        /*@ assert (x135 + x136) != 0.;*/
        x137=1.0/(x135 + x136); POST_CHECK_POW(x137);
        x138=pow(mu, 5); POST_CHECK_POW(x138);
        /*@ assert (x138) >= 0.;*/
        /*@ assert (x138) != 0.;*/
        x139=pow(rn, 4); POST_CHECK_POW(x139);
        /*@ assert (x139) >= 0.;*/
        /*@ assert (x139) != 0.;*/
        x140=pow(un, 3); POST_CHECK_POW(x140);
        /*@ assert (x140) >= 0.;*/
        /*@ assert (x140) != 0.;*/
        x143=pow(mu, 4); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x144=pow(rt1, 5); POST_CHECK_POW(x144);
        x145=-1.4142135623730951454746218587388284504413604736328125*un*x112*x143*x144*x37; POST_CHECK_MUL(x145);
        x146=1.4142135623730951454746218587388284504413604736328125*x39*un*x143*x144*x112; POST_CHECK_MUL(x146);
        x147=pow(rt2, 4); POST_CHECK_POW(x147);
        /*@ assert (x147) >= 0.;*/
        x148=-1.4142135623730951454746218587388284504413604736328125*rt1*un*x112*x143*x147*x37; POST_CHECK_MUL(x148);
        x149=1.4142135623730951454746218587388284504413604736328125*x39*rt1*x147*un*x143*x112; POST_CHECK_MUL(x149);
        x154=-2.828427124746190290949243717477656900882720947265625*un*x112*x113*x143*x25*x37; POST_CHECK_MUL(x154);
        x155=2.828427124746190290949243717477656900882720947265625*x39*un*x143*x112*x113*x25; POST_CHECK_MUL(x155);
        x163=pow(rt2, 5); POST_CHECK_POW(x163);
        x164=pow(rt2, 3); POST_CHECK_POW(x164);
        x165=2.0*mu*rt1*rt2*un*x109*x37 + 2.0*x39*rt1*rt2*mu*un*x109 + 2.0*rt1*un*x111*x26*x163*x37 - 2.0*x39*rt1*un*x111*x26*x163 + 2.0*rt1*un*x138*x139*x164*x37 - 2.0*x39*rt1*un*x138*x139*x164 + 2.0*rt1*x111*x26*x164*x140*x37 - 2.0*x39*rt1*x111*x26*x164*x140 + 2.0*rt2*un*x111*x144*x26*x37 - 2.0*x39*rt2*un*x144*x111*x26 + 2.0*rt2*un*x113*x138*x139*x37 - 2.0*x39*rt2*un*x138*x139*x113 + 2.0*rt2*x111*x113*x140*x26*x37 - 2.0*x39*rt2*x111*x26*x113*x140 + 4.0*un*x111*x113*x26*x164*x37 - 4.0*x39*un*x111*x26*x113*x164; POST_CHECK_ADD(x165);

        /* Assignment result[1, 2]=x137*(x145 + x146 + x148 + x149 + x154 + x155 + x165) */
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x109=pow(x64, 3.0/2.0); POST_CHECK_POW(x109);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x112=pow(rn, 3); POST_CHECK_POW(x112);
        /*@ assert (x112) >= 0.;*/
        /*@ assert (x112) != 0.;*/
        x113=pow(rt1, 3); POST_CHECK_POW(x113);
        x134=4.0*x39*x109*x37; POST_CHECK_MUL(x134);
        x135=x24*x134; POST_CHECK_MUL(x135);
        x136=x25*x134; POST_CHECK_MUL(x136);
        /*@ assert (x135 + x136) != 0.;*/
        x137=1.0/(x135 + x136); POST_CHECK_POW(x137);
        x138=pow(mu, 5); POST_CHECK_POW(x138);
        /*@ assert (x138) >= 0.;*/
        /*@ assert (x138) != 0.;*/
        x139=pow(rn, 4); POST_CHECK_POW(x139);
        /*@ assert (x139) >= 0.;*/
        /*@ assert (x139) != 0.;*/
        x140=pow(un, 3); POST_CHECK_POW(x140);
        /*@ assert (x140) >= 0.;*/
        /*@ assert (x140) != 0.;*/
        x143=pow(mu, 4); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x144=pow(rt1, 5); POST_CHECK_POW(x144);
        x145=-1.4142135623730951454746218587388284504413604736328125*un*x112*x143*x144*x37; POST_CHECK_MUL(x145);
        x146=1.4142135623730951454746218587388284504413604736328125*x39*un*x143*x144*x112; POST_CHECK_MUL(x146);
        x147=pow(rt2, 4); POST_CHECK_POW(x147);
        /*@ assert (x147) >= 0.;*/
        x148=-1.4142135623730951454746218587388284504413604736328125*rt1*un*x112*x143*x147*x37; POST_CHECK_MUL(x148);
        x149=1.4142135623730951454746218587388284504413604736328125*x39*rt1*x147*un*x143*x112; POST_CHECK_MUL(x149);
        x154=-2.828427124746190290949243717477656900882720947265625*un*x112*x113*x143*x25*x37; POST_CHECK_MUL(x154);
        x155=2.828427124746190290949243717477656900882720947265625*x39*un*x143*x112*x113*x25; POST_CHECK_MUL(x155);
        x163=pow(rt2, 5); POST_CHECK_POW(x163);
        x164=pow(rt2, 3); POST_CHECK_POW(x164);
        x165=2.0*mu*rt1*rt2*un*x109*x37 + 2.0*x39*rt1*rt2*mu*un*x109 + 2.0*rt1*un*x111*x26*x163*x37 - 2.0*x39*rt1*un*x111*x26*x163 + 2.0*rt1*un*x138*x139*x164*x37 - 2.0*x39*rt1*un*x138*x139*x164 + 2.0*rt1*x111*x26*x164*x140*x37 - 2.0*x39*rt1*x111*x26*x164*x140 + 2.0*rt2*un*x111*x144*x26*x37 - 2.0*x39*rt2*un*x144*x111*x26 + 2.0*rt2*un*x113*x138*x139*x37 - 2.0*x39*rt2*un*x138*x139*x113 + 2.0*rt2*x111*x113*x140*x26*x37 - 2.0*x39*rt2*x111*x26*x113*x140 + 4.0*un*x111*x113*x26*x164*x37 - 4.0*x39*un*x111*x26*x113*x164; POST_CHECK_ADD(x165);
        result[7] = x137*(x145 + x146 + x148 + x149 + x154 + x155 + x165);

    }
    else if (x60)
    {
        DEBUG_PRINT("Case (x60) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x75=2*x16; POST_CHECK_MUL(x75);
        x81=ut2*x70; POST_CHECK_MUL(x81);
        x82=ut2*x1 + x81*x15; POST_CHECK_ADD(x82);
        x83=x73*x17; POST_CHECK_MUL(x83);
        x84=x29*x69; POST_CHECK_MUL(x84);
        x85=(1.0/2.0)*x19*(x75 + 2*x84); POST_CHECK_MUL(x85);
        x86=x51*(x83 + x85); POST_CHECK_MUL(x86);
        x87=x46*(x82 + x86); POST_CHECK_MUL(x87);
        x88=x55*(x82 - x86); POST_CHECK_MUL(x88);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x117=x17*x51; POST_CHECK_MUL(x117);
        x166=0.5*ut1*ut2*x1*x69*x51; POST_CHECK_MUL(x166);
        x167=-x83 - x85; POST_CHECK_ADD(x167);
        x168=x166 + 0.5*x17*x114*x167; POST_CHECK_ADD(x168);

        /* Assignment result[1, 2]=-x117*x87 + x117*x88 - x45*x168 + x54*x168 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x75=2*x16; POST_CHECK_MUL(x75);
        x81=ut2*x70; POST_CHECK_MUL(x81);
        x82=ut2*x1 + x81*x15; POST_CHECK_ADD(x82);
        x83=x73*x17; POST_CHECK_MUL(x83);
        x84=x29*x69; POST_CHECK_MUL(x84);
        x85=(1.0/2.0)*x19*(x75 + 2*x84); POST_CHECK_MUL(x85);
        x86=x51*(x83 + x85); POST_CHECK_MUL(x86);
        x87=x46*(x82 + x86); POST_CHECK_MUL(x87);
        x88=x55*(x82 - x86); POST_CHECK_MUL(x88);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x117=x17*x51; POST_CHECK_MUL(x117);
        x166=0.5*ut1*ut2*x1*x69*x51; POST_CHECK_MUL(x166);
        x167=-x83 - x85; POST_CHECK_ADD(x167);
        x168=x166 + 0.5*x17*x114*x167; POST_CHECK_ADD(x168);
        result[7] = -x117*x87 + x117*x88 - x45*x168 + x54*x168;

    }
    else if (x120)
    {
        DEBUG_PRINT("Case (x120) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x75=2*x16; POST_CHECK_MUL(x75);
        x81=ut2*x70; POST_CHECK_MUL(x81);
        x82=ut2*x1 + x81*x15; POST_CHECK_ADD(x82);
        x83=x73*x17; POST_CHECK_MUL(x83);
        x84=x29*x69; POST_CHECK_MUL(x84);
        x85=(1.0/2.0)*x19*(x75 + 2*x84); POST_CHECK_MUL(x85);
        x86=x51*(x83 + x85); POST_CHECK_MUL(x86);
        x87=x46*(x82 + x86); POST_CHECK_MUL(x87);
        x88=x55*(x82 - x86); POST_CHECK_MUL(x88);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x119=random1*x118; POST_CHECK_MUL(x119);

        /* Assignment result[1, 2]=-x119*x87 + x119*x88 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x75=2*x16; POST_CHECK_MUL(x75);
        x81=ut2*x70; POST_CHECK_MUL(x81);
        x82=ut2*x1 + x81*x15; POST_CHECK_ADD(x82);
        x83=x73*x17; POST_CHECK_MUL(x83);
        x84=x29*x69; POST_CHECK_MUL(x84);
        x85=(1.0/2.0)*x19*(x75 + 2*x84); POST_CHECK_MUL(x85);
        x86=x51*(x83 + x85); POST_CHECK_MUL(x86);
        x87=x46*(x82 + x86); POST_CHECK_MUL(x87);
        x88=x55*(x82 - x86); POST_CHECK_MUL(x88);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x119=random1*x118; POST_CHECK_MUL(x119);
        result[7] = -x119*x87 + x119*x88;

    }
    /*@ assert (result[7]) >= 0.;*/

    /* Assignment result[2, 2]=Piecewise((x133, x34), (x137*(2.0*x39*x147*un*x111*x26*x24 + 4.0*x39*x152*un*x111*x26*x25 - 2.0*un*x111*x26*x192*x37 + 2.0*x39*un*x111*x26*x192 - 2.0*un*x138*x139*x152*x37 + 2.0*x39*x152*un*x138*x139 - 2.0*x111*x140*x152*x26*x37 + 2.0*x39*x152*x111*x26*x140 - x131*x151 + x141 - x150*x153 - 4.0*x156 - 4.0*x157 - 2.0*x158 - 2.0*x159 + x185 + x186 + x187 + x188 + x189 + x190), x44), (mu - x183*x87 + x183*x88 - x45*x193 + x54*x193, x60), (mu - x184*x87 + x184*x88, x120)) */
    double x192;
    double x193;
    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x123=4.24264068711928477029005080112256109714508056640625*x2; POST_CHECK_MUL(x123);
        x124=4.24264068711928477029005080112256109714508056640625*x3; POST_CHECK_MUL(x124);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x126=2.0*x125; POST_CHECK_MUL(x126);
        x127=x1*x2; POST_CHECK_MUL(x127);
        x128=x2*x111; POST_CHECK_MUL(x128);
        x129=x1*x124; POST_CHECK_MUL(x129);
        x130=-x129; POST_CHECK_MUL(x130);
        x131=2.0*x111; POST_CHECK_MUL(x131);
        x132=x131*x3; POST_CHECK_MUL(x132);
        x133=x122*(-x123 + x124 + x126 - 57.9827560572968963015227927826344966888427734375*x127 + 11.3137084989847611637969748699106276035308837890625*x127*x3 - 34.0*x128 + x130 + x132 - 26.0*x89 + 16.0*x90); POST_CHECK_MUL(x133);

        /* Assignment result[2, 2]=x133 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x123=4.24264068711928477029005080112256109714508056640625*x2; POST_CHECK_MUL(x123);
        x124=4.24264068711928477029005080112256109714508056640625*x3; POST_CHECK_MUL(x124);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x126=2.0*x125; POST_CHECK_MUL(x126);
        x127=x1*x2; POST_CHECK_MUL(x127);
        x128=x2*x111; POST_CHECK_MUL(x128);
        x129=x1*x124; POST_CHECK_MUL(x129);
        x130=-x129; POST_CHECK_MUL(x130);
        x131=2.0*x111; POST_CHECK_MUL(x131);
        x132=x131*x3; POST_CHECK_MUL(x132);
        x133=x122*(-x123 + x124 + x126 - 57.9827560572968963015227927826344966888427734375*x127 + 11.3137084989847611637969748699106276035308837890625*x127*x3 - 34.0*x128 + x130 + x132 - 26.0*x89 + 16.0*x90); POST_CHECK_MUL(x133);
        result[8] = x133;

    }
    else if (x44)
    {
        DEBUG_PRINT("Case (x44) is True.\n");
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x109=pow(x64, 3.0/2.0); POST_CHECK_POW(x109);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x112=pow(rn, 3); POST_CHECK_POW(x112);
        /*@ assert (x112) >= 0.;*/
        /*@ assert (x112) != 0.;*/
        x131=2.0*x111; POST_CHECK_MUL(x131);
        x134=4.0*x39*x109*x37; POST_CHECK_MUL(x134);
        x135=x24*x134; POST_CHECK_MUL(x135);
        x136=x25*x134; POST_CHECK_MUL(x136);
        /*@ assert (x135 + x136) != 0.;*/
        x137=1.0/(x135 + x136); POST_CHECK_POW(x137);
        x138=pow(mu, 5); POST_CHECK_POW(x138);
        /*@ assert (x138) >= 0.;*/
        /*@ assert (x138) != 0.;*/
        x139=pow(rn, 4); POST_CHECK_POW(x139);
        /*@ assert (x139) >= 0.;*/
        /*@ assert (x139) != 0.;*/
        x140=pow(un, 3); POST_CHECK_POW(x140);
        /*@ assert (x140) >= 0.;*/
        /*@ assert (x140) != 0.;*/
        x141=mu*x135 + mu*x136 - 2.0*un*x138*x139*x24*x25*x37 + 2.0*x39*un*x138*x139*x24*x25 - 2.0*x111*x140*x24*x25*x26*x37 + 2.0*x39*x111*x26*x24*x25*x140; POST_CHECK_ADD(x141);
        x143=pow(mu, 4); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x147=pow(rt2, 4); POST_CHECK_POW(x147);
        /*@ assert (x147) >= 0.;*/
        x150=4.0*x111; POST_CHECK_MUL(x150);
        x151=un*x147*x24*x26*x37; POST_CHECK_MUL(x151);
        x152=pow(rt1, 4); POST_CHECK_POW(x152);
        /*@ assert (x152) >= 0.;*/
        x153=un*x152*x25*x26*x37; POST_CHECK_MUL(x153);
        x156=mu*un*x24*x109*x37; POST_CHECK_MUL(x156);
        x157=x39*mu*un*x24*x109; POST_CHECK_MUL(x157);
        x158=mu*un*x25*x109*x37; POST_CHECK_MUL(x158);
        x159=x39*mu*un*x25*x109; POST_CHECK_MUL(x159);
        x163=pow(rt2, 5); POST_CHECK_POW(x163);
        x164=pow(rt2, 3); POST_CHECK_POW(x164);
        x185=-1.4142135623730951454746218587388284504413604736328125*un*x112*x143*x163*x37; POST_CHECK_MUL(x185);
        x186=1.4142135623730951454746218587388284504413604736328125*x39*un*x143*x112*x163; POST_CHECK_MUL(x186);
        x187=-1.4142135623730951454746218587388284504413604736328125*rt2*un*x112*x143*x152*x37; POST_CHECK_MUL(x187);
        x188=1.4142135623730951454746218587388284504413604736328125*x39*rt2*x152*un*x143*x112; POST_CHECK_MUL(x188);
        x189=-2.828427124746190290949243717477656900882720947265625*un*x112*x143*x24*x164*x37; POST_CHECK_MUL(x189);
        x190=2.828427124746190290949243717477656900882720947265625*x39*un*x143*x112*x24*x164; POST_CHECK_MUL(x190);
        x192=pow(rt1, 6); POST_CHECK_POW(x192);

        /* Assignment result[2, 2]=x137*(2.0*x39*x147*un*x111*x26*x24 + 4.0*x39*x152*un*x111*x26*x25 - 2.0*un*x111*x26*x192*x37 + 2.0*x39*un*x111*x26*x192 - 2.0*un*x138*x139*x152*x37 + 2.0*x39*x152*un*x138*x139 - 2.0*x111*x140*x152*x26*x37 + 2.0*x39*x152*x111*x26*x140 - x131*x151 + x141 - x150*x153 - 4.0*x156 - 4.0*x157 - 2.0*x158 - 2.0*x159 + x185 + x186 + x187 + x188 + x189 + x190) */
        /*@ assert (x24 + x25) >= 0.;*/
        x35=sqrt(x24 + x25); POST_CHECK_POW(x35);
        x36=un*un + x24 + x25 + x27; POST_CHECK_ADD(x36);
        /*@ assert (-2.0*x11*x35 + x36) >= 0.;*/
        x37=sqrt(-2.0*x11*x35 + x36); POST_CHECK_POW(x37);
        /*@ assert (2.0*mu*rn*x35 + x36) >= 0.;*/
        x39=sqrt(2.0*mu*rn*x35 + x36); POST_CHECK_POW(x39);
        x64=x1*x26*x24 + x1*x26*x25; POST_CHECK_ADD(x64);
        /*@ assert (x64) >= 0.;*/
        x109=pow(x64, 3.0/2.0); POST_CHECK_POW(x109);
        x111=pow(mu, 3); POST_CHECK_POW(x111);
        /*@ assert (x111) >= 0.;*/
        /*@ assert (x111) != 0.;*/
        x112=pow(rn, 3); POST_CHECK_POW(x112);
        /*@ assert (x112) >= 0.;*/
        /*@ assert (x112) != 0.;*/
        x131=2.0*x111; POST_CHECK_MUL(x131);
        x134=4.0*x39*x109*x37; POST_CHECK_MUL(x134);
        x135=x24*x134; POST_CHECK_MUL(x135);
        x136=x25*x134; POST_CHECK_MUL(x136);
        /*@ assert (x135 + x136) != 0.;*/
        x137=1.0/(x135 + x136); POST_CHECK_POW(x137);
        x138=pow(mu, 5); POST_CHECK_POW(x138);
        /*@ assert (x138) >= 0.;*/
        /*@ assert (x138) != 0.;*/
        x139=pow(rn, 4); POST_CHECK_POW(x139);
        /*@ assert (x139) >= 0.;*/
        /*@ assert (x139) != 0.;*/
        x140=pow(un, 3); POST_CHECK_POW(x140);
        /*@ assert (x140) >= 0.;*/
        /*@ assert (x140) != 0.;*/
        x141=mu*x135 + mu*x136 - 2.0*un*x138*x139*x24*x25*x37 + 2.0*x39*un*x138*x139*x24*x25 - 2.0*x111*x140*x24*x25*x26*x37 + 2.0*x39*x111*x26*x24*x25*x140; POST_CHECK_ADD(x141);
        x143=pow(mu, 4); POST_CHECK_POW(x143);
        /*@ assert (x143) >= 0.;*/
        /*@ assert (x143) != 0.;*/
        x147=pow(rt2, 4); POST_CHECK_POW(x147);
        /*@ assert (x147) >= 0.;*/
        x150=4.0*x111; POST_CHECK_MUL(x150);
        x151=un*x147*x24*x26*x37; POST_CHECK_MUL(x151);
        x152=pow(rt1, 4); POST_CHECK_POW(x152);
        /*@ assert (x152) >= 0.;*/
        x153=un*x152*x25*x26*x37; POST_CHECK_MUL(x153);
        x156=mu*un*x24*x109*x37; POST_CHECK_MUL(x156);
        x157=x39*mu*un*x24*x109; POST_CHECK_MUL(x157);
        x158=mu*un*x25*x109*x37; POST_CHECK_MUL(x158);
        x159=x39*mu*un*x25*x109; POST_CHECK_MUL(x159);
        x163=pow(rt2, 5); POST_CHECK_POW(x163);
        x164=pow(rt2, 3); POST_CHECK_POW(x164);
        x185=-1.4142135623730951454746218587388284504413604736328125*un*x112*x143*x163*x37; POST_CHECK_MUL(x185);
        x186=1.4142135623730951454746218587388284504413604736328125*x39*un*x143*x112*x163; POST_CHECK_MUL(x186);
        x187=-1.4142135623730951454746218587388284504413604736328125*rt2*un*x112*x143*x152*x37; POST_CHECK_MUL(x187);
        x188=1.4142135623730951454746218587388284504413604736328125*x39*rt2*x152*un*x143*x112; POST_CHECK_MUL(x188);
        x189=-2.828427124746190290949243717477656900882720947265625*un*x112*x143*x24*x164*x37; POST_CHECK_MUL(x189);
        x190=2.828427124746190290949243717477656900882720947265625*x39*un*x143*x112*x24*x164; POST_CHECK_MUL(x190);
        x192=pow(rt1, 6); POST_CHECK_POW(x192);
        result[8] = x137*(2.0*x39*x147*un*x111*x26*x24 + 4.0*x39*x152*un*x111*x26*x25 - 2.0*un*x111*x26*x192*x37 + 2.0*x39*un*x111*x26*x192 - 2.0*un*x138*x139*x152*x37 + 2.0*x39*x152*un*x138*x139 - 2.0*x111*x140*x152*x26*x37 + 2.0*x39*x152*x111*x26*x140 - x131*x151 + x141 - x150*x153 - 4.0*x156 - 4.0*x157 - 2.0*x158 - 2.0*x159 + x185 + x186 + x187 + x188 + x189 + x190);

    }
    else if (x60)
    {
        DEBUG_PRINT("Case (x60) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x75=2*x16; POST_CHECK_MUL(x75);
        x81=ut2*x70; POST_CHECK_MUL(x81);
        x82=ut2*x1 + x81*x15; POST_CHECK_ADD(x82);
        x83=x73*x17; POST_CHECK_MUL(x83);
        x84=x29*x69; POST_CHECK_MUL(x84);
        x85=(1.0/2.0)*x19*(x75 + 2*x84); POST_CHECK_MUL(x85);
        x86=x51*(x83 + x85); POST_CHECK_MUL(x86);
        x87=x46*(x82 + x86); POST_CHECK_MUL(x87);
        x88=x55*(x82 - x86); POST_CHECK_MUL(x88);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x167=-x83 - x85; POST_CHECK_ADD(x167);
        x183=x19*x51; POST_CHECK_MUL(x183);
        x193=0.5*x19*x114*x167 + 0.5*(x16 + x84)*x51; POST_CHECK_ADD(x193);

        /* Assignment result[2, 2]=mu - x183*x87 + x183*x88 - x45*x193 + x54*x193 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x75=2*x16; POST_CHECK_MUL(x75);
        x81=ut2*x70; POST_CHECK_MUL(x81);
        x82=ut2*x1 + x81*x15; POST_CHECK_ADD(x82);
        x83=x73*x17; POST_CHECK_MUL(x83);
        x84=x29*x69; POST_CHECK_MUL(x84);
        x85=(1.0/2.0)*x19*(x75 + 2*x84); POST_CHECK_MUL(x85);
        x86=x51*(x83 + x85); POST_CHECK_MUL(x86);
        x87=x46*(x82 + x86); POST_CHECK_MUL(x87);
        x88=x55*(x82 - x86); POST_CHECK_MUL(x88);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x167=-x83 - x85; POST_CHECK_ADD(x167);
        x183=x19*x51; POST_CHECK_MUL(x183);
        x193=0.5*x19*x114*x167 + 0.5*(x16 + x84)*x51; POST_CHECK_ADD(x193);
        result[8] = mu - x183*x87 + x183*x88 - x45*x193 + x54*x193;

    }
    else if (x120)
    {
        DEBUG_PRINT("Case (x120) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x75=2*x16; POST_CHECK_MUL(x75);
        x81=ut2*x70; POST_CHECK_MUL(x81);
        x82=ut2*x1 + x81*x15; POST_CHECK_ADD(x82);
        x83=x73*x17; POST_CHECK_MUL(x83);
        x84=x29*x69; POST_CHECK_MUL(x84);
        x85=(1.0/2.0)*x19*(x75 + 2*x84); POST_CHECK_MUL(x85);
        x86=x51*(x83 + x85); POST_CHECK_MUL(x86);
        x87=x46*(x82 + x86); POST_CHECK_MUL(x87);
        x88=x55*(x82 - x86); POST_CHECK_MUL(x88);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x184=random2*x118; POST_CHECK_MUL(x184);

        /* Assignment result[2, 2]=mu - x184*x87 + x184*x88 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        /*@ assert (x14) != 0.;*/
        x69=1.0/x14; POST_CHECK_POW(x69);
        x70=mu*x69; POST_CHECK_MUL(x70);
        x73=ut1*ut2*x1*x69; POST_CHECK_MUL(x73);
        x75=2*x16; POST_CHECK_MUL(x75);
        x81=ut2*x70; POST_CHECK_MUL(x81);
        x82=ut2*x1 + x81*x15; POST_CHECK_ADD(x82);
        x83=x73*x17; POST_CHECK_MUL(x83);
        x84=x29*x69; POST_CHECK_MUL(x84);
        x85=(1.0/2.0)*x19*(x75 + 2*x84); POST_CHECK_MUL(x85);
        x86=x51*(x83 + x85); POST_CHECK_MUL(x86);
        x87=x46*(x82 + x86); POST_CHECK_MUL(x87);
        x88=x55*(x82 - x86); POST_CHECK_MUL(x88);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x184=random2*x118; POST_CHECK_MUL(x184);
        result[8] = mu - x184*x87 + x184*x88;

    }
    /*@ assert (result[8]) >= 0.;*/

    /* Assignment result[0, 3]=Piecewise((x5*(-x1*x6 - x1*x7 - x61*x2 + x62 + x90), x34), (mu - x95 - x96, x97)) */

    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x6=0.5*x2; POST_CHECK_MUL(x6);
        x7=0.5*x3; POST_CHECK_MUL(x7);
        x61=0.70710678118654757273731092936941422522068023681640625*mu; POST_CHECK_MUL(x61);
        x62=x61*x3; POST_CHECK_MUL(x62);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);

        /* Assignment result[0, 3]=x5*(-x1*x6 - x1*x7 - x61*x2 + x62 + x90) */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x6=0.5*x2; POST_CHECK_MUL(x6);
        x7=0.5*x3; POST_CHECK_MUL(x7);
        x61=0.70710678118654757273731092936941422522068023681640625*mu; POST_CHECK_MUL(x61);
        x62=x61*x3; POST_CHECK_MUL(x62);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        result[9] = x5*(-x1*x6 - x1*x7 - x61*x2 + x62 + x90);

    }
    else if (x97)
    {
        DEBUG_PRINT("Case (x97) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x49=mu*x19; POST_CHECK_MUL(x49);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x91=x1*rn; POST_CHECK_MUL(x91);
        x92=rt1*x47; POST_CHECK_MUL(x92);
        x93=rt2*x49; POST_CHECK_MUL(x93);
        x94=(x92 + x93)*x51; POST_CHECK_MUL(x94);
        x95=x46*(x91 + x94); POST_CHECK_MUL(x95);
        x96=x55*(x91 - x94); POST_CHECK_MUL(x96);

        /* Assignment result[0, 3]=mu - x95 - x96 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x49=mu*x19; POST_CHECK_MUL(x49);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x91=x1*rn; POST_CHECK_MUL(x91);
        x92=rt1*x47; POST_CHECK_MUL(x92);
        x93=rt2*x49; POST_CHECK_MUL(x93);
        x94=(x92 + x93)*x51; POST_CHECK_MUL(x94);
        x95=x46*(x91 + x94); POST_CHECK_MUL(x95);
        x96=x55*(x91 - x94); POST_CHECK_MUL(x96);
        result[9] = mu - x95 - x96;

    }
    /*@ assert (result[9]) >= 0.;*/

    /* Assignment result[1, 3]=Piecewise((x169, x34), (-x117*x95 + x117*x96 - x45*x171 + x54*x171, x97), (-x119*x95 + x119*x96, x172)) */
    double x169;
    double x170;
    double x171;
    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x6=0.5*x2; POST_CHECK_MUL(x6);
        x7=0.5*x3; POST_CHECK_MUL(x7);
        x98=0.353553390593273786368655464684707112610340118408203125*x2; POST_CHECK_MUL(x98);
        x99=0.353553390593273786368655464684707112610340118408203125*x3; POST_CHECK_MUL(x99);
        x169=x5*(-mu*x6 - mu*x7 - x1*x98 + x1*x99); POST_CHECK_MUL(x169);

        /* Assignment result[1, 3]=x169 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x6=0.5*x2; POST_CHECK_MUL(x6);
        x7=0.5*x3; POST_CHECK_MUL(x7);
        x98=0.353553390593273786368655464684707112610340118408203125*x2; POST_CHECK_MUL(x98);
        x99=0.353553390593273786368655464684707112610340118408203125*x3; POST_CHECK_MUL(x99);
        x169=x5*(-mu*x6 - mu*x7 - x1*x98 + x1*x99); POST_CHECK_MUL(x169);
        result[10] = x169;

    }
    else if (x97)
    {
        DEBUG_PRINT("Case (x97) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x49=mu*x19; POST_CHECK_MUL(x49);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x91=x1*rn; POST_CHECK_MUL(x91);
        x92=rt1*x47; POST_CHECK_MUL(x92);
        x93=rt2*x49; POST_CHECK_MUL(x93);
        x94=(x92 + x93)*x51; POST_CHECK_MUL(x94);
        x95=x46*(x91 + x94); POST_CHECK_MUL(x95);
        x96=x55*(x91 - x94); POST_CHECK_MUL(x96);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x117=x17*x51; POST_CHECK_MUL(x117);
        x170=-x92 - x93; POST_CHECK_ADD(x170);
        x171=0.5*mu*rt1*x51 + 0.5*x17*x114*x170; POST_CHECK_ADD(x171);

        /* Assignment result[1, 3]=-x117*x95 + x117*x96 - x45*x171 + x54*x171 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x49=mu*x19; POST_CHECK_MUL(x49);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x91=x1*rn; POST_CHECK_MUL(x91);
        x92=rt1*x47; POST_CHECK_MUL(x92);
        x93=rt2*x49; POST_CHECK_MUL(x93);
        x94=(x92 + x93)*x51; POST_CHECK_MUL(x94);
        x95=x46*(x91 + x94); POST_CHECK_MUL(x95);
        x96=x55*(x91 - x94); POST_CHECK_MUL(x96);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x117=x17*x51; POST_CHECK_MUL(x117);
        x170=-x92 - x93; POST_CHECK_ADD(x170);
        x171=0.5*mu*rt1*x51 + 0.5*x17*x114*x170; POST_CHECK_ADD(x171);
        result[10] = -x117*x95 + x117*x96 - x45*x171 + x54*x171;

    }
    else if (x172)
    {
        DEBUG_PRINT("Case (x172) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x49=mu*x19; POST_CHECK_MUL(x49);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x91=x1*rn; POST_CHECK_MUL(x91);
        x92=rt1*x47; POST_CHECK_MUL(x92);
        x93=rt2*x49; POST_CHECK_MUL(x93);
        x94=(x92 + x93)*x51; POST_CHECK_MUL(x94);
        x95=x46*(x91 + x94); POST_CHECK_MUL(x95);
        x96=x55*(x91 - x94); POST_CHECK_MUL(x96);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x119=random1*x118; POST_CHECK_MUL(x119);

        /* Assignment result[1, 3]=-x119*x95 + x119*x96 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x49=mu*x19; POST_CHECK_MUL(x49);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x91=x1*rn; POST_CHECK_MUL(x91);
        x92=rt1*x47; POST_CHECK_MUL(x92);
        x93=rt2*x49; POST_CHECK_MUL(x93);
        x94=(x92 + x93)*x51; POST_CHECK_MUL(x94);
        x95=x46*(x91 + x94); POST_CHECK_MUL(x95);
        x96=x55*(x91 - x94); POST_CHECK_MUL(x96);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x119=random1*x118; POST_CHECK_MUL(x119);
        result[10] = -x119*x95 + x119*x96;

    }
    /*@ assert (result[10]) >= 0.;*/

    /* Assignment result[2, 3]=Piecewise((x169, x34), (-x183*x95 + x183*x96 - x45*x194 + x54*x194, x97), (-x184*x95 + x184*x96, x172)) */
    double x194;
    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x6=0.5*x2; POST_CHECK_MUL(x6);
        x7=0.5*x3; POST_CHECK_MUL(x7);
        x98=0.353553390593273786368655464684707112610340118408203125*x2; POST_CHECK_MUL(x98);
        x99=0.353553390593273786368655464684707112610340118408203125*x3; POST_CHECK_MUL(x99);
        x169=x5*(-mu*x6 - mu*x7 - x1*x98 + x1*x99); POST_CHECK_MUL(x169);

        /* Assignment result[2, 3]=x169 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x6=0.5*x2; POST_CHECK_MUL(x6);
        x7=0.5*x3; POST_CHECK_MUL(x7);
        x98=0.353553390593273786368655464684707112610340118408203125*x2; POST_CHECK_MUL(x98);
        x99=0.353553390593273786368655464684707112610340118408203125*x3; POST_CHECK_MUL(x99);
        x169=x5*(-mu*x6 - mu*x7 - x1*x98 + x1*x99); POST_CHECK_MUL(x169);
        result[11] = x169;

    }
    else if (x97)
    {
        DEBUG_PRINT("Case (x97) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x49=mu*x19; POST_CHECK_MUL(x49);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x91=x1*rn; POST_CHECK_MUL(x91);
        x92=rt1*x47; POST_CHECK_MUL(x92);
        x93=rt2*x49; POST_CHECK_MUL(x93);
        x94=(x92 + x93)*x51; POST_CHECK_MUL(x94);
        x95=x46*(x91 + x94); POST_CHECK_MUL(x95);
        x96=x55*(x91 - x94); POST_CHECK_MUL(x96);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x170=-x92 - x93; POST_CHECK_ADD(x170);
        x183=x19*x51; POST_CHECK_MUL(x183);
        x194=0.5*mu*rt2*x51 + 0.5*x19*x114*x170; POST_CHECK_ADD(x194);

        /* Assignment result[2, 3]=-x183*x95 + x183*x96 - x45*x194 + x54*x194 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x49=mu*x19; POST_CHECK_MUL(x49);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x91=x1*rn; POST_CHECK_MUL(x91);
        x92=rt1*x47; POST_CHECK_MUL(x92);
        x93=rt2*x49; POST_CHECK_MUL(x93);
        x94=(x92 + x93)*x51; POST_CHECK_MUL(x94);
        x95=x46*(x91 + x94); POST_CHECK_MUL(x95);
        x96=x55*(x91 - x94); POST_CHECK_MUL(x96);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x170=-x92 - x93; POST_CHECK_ADD(x170);
        x183=x19*x51; POST_CHECK_MUL(x183);
        x194=0.5*mu*rt2*x51 + 0.5*x19*x114*x170; POST_CHECK_ADD(x194);
        result[11] = -x183*x95 + x183*x96 - x45*x194 + x54*x194;

    }
    else if (x172)
    {
        DEBUG_PRINT("Case (x172) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x49=mu*x19; POST_CHECK_MUL(x49);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x91=x1*rn; POST_CHECK_MUL(x91);
        x92=rt1*x47; POST_CHECK_MUL(x92);
        x93=rt2*x49; POST_CHECK_MUL(x93);
        x94=(x92 + x93)*x51; POST_CHECK_MUL(x94);
        x95=x46*(x91 + x94); POST_CHECK_MUL(x95);
        x96=x55*(x91 - x94); POST_CHECK_MUL(x96);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x184=random2*x118; POST_CHECK_MUL(x184);

        /* Assignment result[2, 3]=-x184*x95 + x184*x96 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        x47=mu*x17; POST_CHECK_MUL(x47);
        x49=mu*x19; POST_CHECK_MUL(x49);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x91=x1*rn; POST_CHECK_MUL(x91);
        x92=rt1*x47; POST_CHECK_MUL(x92);
        x93=rt2*x49; POST_CHECK_MUL(x93);
        x94=(x92 + x93)*x51; POST_CHECK_MUL(x94);
        x95=x46*(x91 + x94); POST_CHECK_MUL(x95);
        x96=x55*(x91 - x94); POST_CHECK_MUL(x96);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x184=random2*x118; POST_CHECK_MUL(x184);
        result[11] = -x184*x95 + x184*x96;

    }
    /*@ assert (result[11]) >= 0.;*/

    /* Assignment result[0, 4]=Piecewise((x100, x34), (-x103 - x104, x97)) */
    double x100;
    double x101;
    double x102;
    double x103;
    double x104;
    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x6=0.5*x2; POST_CHECK_MUL(x6);
        x7=0.5*x3; POST_CHECK_MUL(x7);
        x8=-x6 - x7; POST_CHECK_ADD(x8);
        x98=0.353553390593273786368655464684707112610340118408203125*x2; POST_CHECK_MUL(x98);
        x99=0.353553390593273786368655464684707112610340118408203125*x3; POST_CHECK_MUL(x99);
        x100=x5*(-mu*x98 + mu*x99 + x8); POST_CHECK_MUL(x100);

        /* Assignment result[0, 4]=x100 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x6=0.5*x2; POST_CHECK_MUL(x6);
        x7=0.5*x3; POST_CHECK_MUL(x7);
        x8=-x6 - x7; POST_CHECK_ADD(x8);
        x98=0.353553390593273786368655464684707112610340118408203125*x2; POST_CHECK_MUL(x98);
        x99=0.353553390593273786368655464684707112610340118408203125*x3; POST_CHECK_MUL(x99);
        x100=x5*(-mu*x98 + mu*x99 + x8); POST_CHECK_MUL(x100);
        result[12] = x100;

    }
    else if (x97)
    {
        DEBUG_PRINT("Case (x97) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x102=x17*x101; POST_CHECK_MUL(x102);
        x103=x46*(rt1 + x102); POST_CHECK_MUL(x103);
        x104=(rt1 - x102)*x55; POST_CHECK_MUL(x104);

        /* Assignment result[0, 4]=-x103 - x104 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x102=x17*x101; POST_CHECK_MUL(x102);
        x103=x46*(rt1 + x102); POST_CHECK_MUL(x103);
        x104=(rt1 - x102)*x55; POST_CHECK_MUL(x104);
        result[12] = -x103 - x104;

    }
    /*@ assert (result[12]) >= 0.;*/

    /* Assignment result[1, 4]=Piecewise((x176, x34), (1 - x117*x103 + x117*x104 - x45*x178 + x54*x178, x97), (1 - x119*x103 + x119*x104, x172)) */
    double x173;
    double x174;
    double x175;
    double x176;
    double x177;
    double x178;
    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x127=x1*x2; POST_CHECK_MUL(x127);
        x173=9.8994949366116653521885382360778748989105224609375*x2; POST_CHECK_MUL(x173);
        x174=4.0*x125; POST_CHECK_MUL(x174);
        x175=1.4142135623730951454746218587388284504413604736328125*x3; POST_CHECK_MUL(x175);
        x176=x122*(-x1*x175 + x121 - 15.5563491861040450459086059709079563617706298828125*x127 - x173 - x174 + 9.8994949366116653521885382360778748989105224609375*x3 - 20.0*x89); POST_CHECK_MUL(x176);

        /* Assignment result[1, 4]=x176 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x127=x1*x2; POST_CHECK_MUL(x127);
        x173=9.8994949366116653521885382360778748989105224609375*x2; POST_CHECK_MUL(x173);
        x174=4.0*x125; POST_CHECK_MUL(x174);
        x175=1.4142135623730951454746218587388284504413604736328125*x3; POST_CHECK_MUL(x175);
        x176=x122*(-x1*x175 + x121 - 15.5563491861040450459086059709079563617706298828125*x127 - x173 - x174 + 9.8994949366116653521885382360778748989105224609375*x3 - 20.0*x89); POST_CHECK_MUL(x176);
        result[13] = x176;

    }
    else if (x97)
    {
        DEBUG_PRINT("Case (x97) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x102=x17*x101; POST_CHECK_MUL(x102);
        x103=x46*(rt1 + x102); POST_CHECK_MUL(x103);
        x104=(rt1 - x102)*x55; POST_CHECK_MUL(x104);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x117=x17*x51; POST_CHECK_MUL(x117);
        x177=0.5*mu*rn*x51; POST_CHECK_MUL(x177);
        x178=-0.5*mu*rn*x18*x114 + x177; POST_CHECK_ADD(x178);

        /* Assignment result[1, 4]=1 - x117*x103 + x117*x104 - x45*x178 + x54*x178 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x102=x17*x101; POST_CHECK_MUL(x102);
        x103=x46*(rt1 + x102); POST_CHECK_MUL(x103);
        x104=(rt1 - x102)*x55; POST_CHECK_MUL(x104);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x117=x17*x51; POST_CHECK_MUL(x117);
        x177=0.5*mu*rn*x51; POST_CHECK_MUL(x177);
        x178=-0.5*mu*rn*x18*x114 + x177; POST_CHECK_ADD(x178);
        result[13] = 1 - x117*x103 + x117*x104 - x45*x178 + x54*x178;

    }
    else if (x172)
    {
        DEBUG_PRINT("Case (x172) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x102=x17*x101; POST_CHECK_MUL(x102);
        x103=x46*(rt1 + x102); POST_CHECK_MUL(x103);
        x104=(rt1 - x102)*x55; POST_CHECK_MUL(x104);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x119=random1*x118; POST_CHECK_MUL(x119);

        /* Assignment result[1, 4]=1 - x119*x103 + x119*x104 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x102=x17*x101; POST_CHECK_MUL(x102);
        x103=x46*(rt1 + x102); POST_CHECK_MUL(x103);
        x104=(rt1 - x102)*x55; POST_CHECK_MUL(x104);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x119=random1*x118; POST_CHECK_MUL(x119);
        result[13] = 1 - x119*x103 + x119*x104;

    }
    /*@ assert (result[13]) >= 0.;*/

    /* Assignment result[2, 4]=Piecewise((x179, x34), (x181 - x183*x103 + x183*x104, x97), (-x184*x103 + x184*x104, x172)) */
    double x179;
    double x180;
    double x181;
    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x9=1.4142135623730951454746218587388284504413604736328125*x2; POST_CHECK_MUL(x9);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x124=4.24264068711928477029005080112256109714508056640625*x3; POST_CHECK_MUL(x124);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x129=x1*x124; POST_CHECK_MUL(x129);
        x130=-x129; POST_CHECK_MUL(x130);
        x173=9.8994949366116653521885382360778748989105224609375*x2; POST_CHECK_MUL(x173);
        x174=4.0*x125; POST_CHECK_MUL(x174);
        x175=1.4142135623730951454746218587388284504413604736328125*x3; POST_CHECK_MUL(x175);
        x179=x122*(x1*x173 + x130 + x174 + x175 + 4.0*x89 - x9); POST_CHECK_MUL(x179);

        /* Assignment result[2, 4]=x179 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x9=1.4142135623730951454746218587388284504413604736328125*x2; POST_CHECK_MUL(x9);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x124=4.24264068711928477029005080112256109714508056640625*x3; POST_CHECK_MUL(x124);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x129=x1*x124; POST_CHECK_MUL(x129);
        x130=-x129; POST_CHECK_MUL(x130);
        x173=9.8994949366116653521885382360778748989105224609375*x2; POST_CHECK_MUL(x173);
        x174=4.0*x125; POST_CHECK_MUL(x174);
        x175=1.4142135623730951454746218587388284504413604736328125*x3; POST_CHECK_MUL(x175);
        x179=x122*(x1*x173 + x130 + x174 + x175 + 4.0*x89 - x9); POST_CHECK_MUL(x179);
        result[14] = x179;

    }
    else if (x97)
    {
        DEBUG_PRINT("Case (x97) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x102=x17*x101; POST_CHECK_MUL(x102);
        x103=x46*(rt1 + x102); POST_CHECK_MUL(x103);
        x104=(rt1 - x102)*x55; POST_CHECK_MUL(x104);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x180=0.5*mu*rn*x17*x19*x114; POST_CHECK_MUL(x180);
        x181=x180*x45 - x180*x54; POST_CHECK_ADD(x181);
        x183=x19*x51; POST_CHECK_MUL(x183);

        /* Assignment result[2, 4]=x181 - x183*x103 + x183*x104 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x102=x17*x101; POST_CHECK_MUL(x102);
        x103=x46*(rt1 + x102); POST_CHECK_MUL(x103);
        x104=(rt1 - x102)*x55; POST_CHECK_MUL(x104);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x180=0.5*mu*rn*x17*x19*x114; POST_CHECK_MUL(x180);
        x181=x180*x45 - x180*x54; POST_CHECK_ADD(x181);
        x183=x19*x51; POST_CHECK_MUL(x183);
        result[14] = x181 - x183*x103 + x183*x104;

    }
    else if (x172)
    {
        DEBUG_PRINT("Case (x172) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x102=x17*x101; POST_CHECK_MUL(x102);
        x103=x46*(rt1 + x102); POST_CHECK_MUL(x103);
        x104=(rt1 - x102)*x55; POST_CHECK_MUL(x104);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x184=random2*x118; POST_CHECK_MUL(x184);

        /* Assignment result[2, 4]=-x184*x103 + x184*x104 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x102=x17*x101; POST_CHECK_MUL(x102);
        x103=x46*(rt1 + x102); POST_CHECK_MUL(x103);
        x104=(rt1 - x102)*x55; POST_CHECK_MUL(x104);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x184=random2*x118; POST_CHECK_MUL(x184);
        result[14] = -x184*x103 + x184*x104;

    }
    /*@ assert (result[14]) >= 0.;*/

    /* Assignment result[0, 5]=Piecewise((x100, x34), (-x106 - x107, x97)) */
    double x105;
    double x106;
    double x107;
    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x6=0.5*x2; POST_CHECK_MUL(x6);
        x7=0.5*x3; POST_CHECK_MUL(x7);
        x8=-x6 - x7; POST_CHECK_ADD(x8);
        x98=0.353553390593273786368655464684707112610340118408203125*x2; POST_CHECK_MUL(x98);
        x99=0.353553390593273786368655464684707112610340118408203125*x3; POST_CHECK_MUL(x99);
        x100=x5*(-mu*x98 + mu*x99 + x8); POST_CHECK_MUL(x100);

        /* Assignment result[0, 5]=x100 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        /*@ assert (x3) != 0.;*/
        x4=1.0/x3; POST_CHECK_POW(x4);
        /*@ assert (x2) != 0.;*/
        x5=x4/x2; POST_CHECK_MUL(x5);
        x6=0.5*x2; POST_CHECK_MUL(x6);
        x7=0.5*x3; POST_CHECK_MUL(x7);
        x8=-x6 - x7; POST_CHECK_ADD(x8);
        x98=0.353553390593273786368655464684707112610340118408203125*x2; POST_CHECK_MUL(x98);
        x99=0.353553390593273786368655464684707112610340118408203125*x3; POST_CHECK_MUL(x99);
        x100=x5*(-mu*x98 + mu*x99 + x8); POST_CHECK_MUL(x100);
        result[15] = x100;

    }
    else if (x97)
    {
        DEBUG_PRINT("Case (x97) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x105=x19*x101; POST_CHECK_MUL(x105);
        x106=x46*(rt2 + x105); POST_CHECK_MUL(x106);
        x107=(rt2 - x105)*x55; POST_CHECK_MUL(x107);

        /* Assignment result[0, 5]=-x106 - x107 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x105=x19*x101; POST_CHECK_MUL(x105);
        x106=x46*(rt2 + x105); POST_CHECK_MUL(x106);
        x107=(rt2 - x105)*x55; POST_CHECK_MUL(x107);
        result[15] = -x106 - x107;

    }
    /*@ assert (result[15]) >= 0.;*/

    /* Assignment result[1, 5]=Piecewise((x179, x34), (-x117*x106 + x117*x107 + x181, x97), (-x119*x106 + x119*x107, x172)) */

    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x9=1.4142135623730951454746218587388284504413604736328125*x2; POST_CHECK_MUL(x9);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x124=4.24264068711928477029005080112256109714508056640625*x3; POST_CHECK_MUL(x124);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x129=x1*x124; POST_CHECK_MUL(x129);
        x130=-x129; POST_CHECK_MUL(x130);
        x173=9.8994949366116653521885382360778748989105224609375*x2; POST_CHECK_MUL(x173);
        x174=4.0*x125; POST_CHECK_MUL(x174);
        x175=1.4142135623730951454746218587388284504413604736328125*x3; POST_CHECK_MUL(x175);
        x179=x122*(x1*x173 + x130 + x174 + x175 + 4.0*x89 - x9); POST_CHECK_MUL(x179);

        /* Assignment result[1, 5]=x179 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x9=1.4142135623730951454746218587388284504413604736328125*x2; POST_CHECK_MUL(x9);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x124=4.24264068711928477029005080112256109714508056640625*x3; POST_CHECK_MUL(x124);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x129=x1*x124; POST_CHECK_MUL(x129);
        x130=-x129; POST_CHECK_MUL(x130);
        x173=9.8994949366116653521885382360778748989105224609375*x2; POST_CHECK_MUL(x173);
        x174=4.0*x125; POST_CHECK_MUL(x174);
        x175=1.4142135623730951454746218587388284504413604736328125*x3; POST_CHECK_MUL(x175);
        x179=x122*(x1*x173 + x130 + x174 + x175 + 4.0*x89 - x9); POST_CHECK_MUL(x179);
        result[16] = x179;

    }
    else if (x97)
    {
        DEBUG_PRINT("Case (x97) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x105=x19*x101; POST_CHECK_MUL(x105);
        x106=x46*(rt2 + x105); POST_CHECK_MUL(x106);
        x107=(rt2 - x105)*x55; POST_CHECK_MUL(x107);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x117=x17*x51; POST_CHECK_MUL(x117);
        x180=0.5*mu*rn*x17*x19*x114; POST_CHECK_MUL(x180);
        x181=x180*x45 - x180*x54; POST_CHECK_ADD(x181);

        /* Assignment result[1, 5]=-x117*x106 + x117*x107 + x181 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x105=x19*x101; POST_CHECK_MUL(x105);
        x106=x46*(rt2 + x105); POST_CHECK_MUL(x106);
        x107=(rt2 - x105)*x55; POST_CHECK_MUL(x107);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x117=x17*x51; POST_CHECK_MUL(x117);
        x180=0.5*mu*rn*x17*x19*x114; POST_CHECK_MUL(x180);
        x181=x180*x45 - x180*x54; POST_CHECK_ADD(x181);
        result[16] = -x117*x106 + x117*x107 + x181;

    }
    else if (x172)
    {
        DEBUG_PRINT("Case (x172) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x105=x19*x101; POST_CHECK_MUL(x105);
        x106=x46*(rt2 + x105); POST_CHECK_MUL(x106);
        x107=(rt2 - x105)*x55; POST_CHECK_MUL(x107);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x119=random1*x118; POST_CHECK_MUL(x119);

        /* Assignment result[1, 5]=-x119*x106 + x119*x107 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x105=x19*x101; POST_CHECK_MUL(x105);
        x106=x46*(rt2 + x105); POST_CHECK_MUL(x106);
        x107=(rt2 - x105)*x55; POST_CHECK_MUL(x107);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x119=random1*x118; POST_CHECK_MUL(x119);
        result[16] = -x119*x106 + x119*x107;

    }
    /*@ assert (result[16]) >= 0.;*/

    /* Assignment result[2, 5]=Piecewise((x176, x34), (1 - x183*x106 + x183*x107 - x45*x195 + x54*x195, x97), (1 - x184*x106 + x184*x107, x172)) */
    double x195;
    if (x34)
    {
        DEBUG_PRINT("Case (x34) is True.\n");
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x127=x1*x2; POST_CHECK_MUL(x127);
        x173=9.8994949366116653521885382360778748989105224609375*x2; POST_CHECK_MUL(x173);
        x174=4.0*x125; POST_CHECK_MUL(x174);
        x175=1.4142135623730951454746218587388284504413604736328125*x3; POST_CHECK_MUL(x175);
        x176=x122*(-x1*x175 + x121 - 15.5563491861040450459086059709079563617706298828125*x127 - x173 - x174 + 9.8994949366116653521885382360778748989105224609375*x3 - 20.0*x89); POST_CHECK_MUL(x176);

        /* Assignment result[2, 5]=x176 */
        /*@ assert (-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1) >= 0.;*/
        x2=sqrt(-2.828427124746190290949243717477656900882720947265625*mu + 3.0 + x1); POST_CHECK_POW(x2);
        /*@ assert (8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1) >= 0.;*/
        x3=sqrt(8.4852813742385695405801016022451221942901611328125*mu + 3.0 + 9.0*x1); POST_CHECK_POW(x3);
        x10=x2*x3; POST_CHECK_MUL(x10);
        x89=x2*mu; POST_CHECK_MUL(x89);
        x90=x89*x3; POST_CHECK_MUL(x90);
        x121=16.0*x10 + 11.3137084989847611637969748699106276035308837890625*x90; POST_CHECK_ADD(x121);
        /*@ assert (x121) != 0.;*/
        x122=1.0/x121; POST_CHECK_POW(x122);
        x125=mu*x3; POST_CHECK_MUL(x125);
        x127=x1*x2; POST_CHECK_MUL(x127);
        x173=9.8994949366116653521885382360778748989105224609375*x2; POST_CHECK_MUL(x173);
        x174=4.0*x125; POST_CHECK_MUL(x174);
        x175=1.4142135623730951454746218587388284504413604736328125*x3; POST_CHECK_MUL(x175);
        x176=x122*(-x1*x175 + x121 - 15.5563491861040450459086059709079563617706298828125*x127 - x173 - x174 + 9.8994949366116653521885382360778748989105224609375*x3 - 20.0*x89); POST_CHECK_MUL(x176);
        result[17] = x176;

    }
    else if (x97)
    {
        DEBUG_PRINT("Case (x97) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x105=x19*x101; POST_CHECK_MUL(x105);
        x106=x46*(rt2 + x105); POST_CHECK_MUL(x106);
        x107=(rt2 - x105)*x55; POST_CHECK_MUL(x107);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x177=0.5*mu*rn*x51; POST_CHECK_MUL(x177);
        x183=x19*x51; POST_CHECK_MUL(x183);
        x195=-0.5*mu*rn*x20*x114 + x177; POST_CHECK_ADD(x195);

        /* Assignment result[2, 5]=1 - x183*x106 + x183*x107 - x45*x195 + x54*x195 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x105=x19*x101; POST_CHECK_MUL(x105);
        x106=x46*(rt2 + x105); POST_CHECK_MUL(x106);
        x107=(rt2 - x105)*x55; POST_CHECK_MUL(x107);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x114=pow(x21, -3.0/2.0); POST_CHECK_POW(x114);
        x177=0.5*mu*rn*x51; POST_CHECK_MUL(x177);
        x183=x19*x51; POST_CHECK_MUL(x183);
        x195=-0.5*mu*rn*x20*x114 + x177; POST_CHECK_ADD(x195);
        result[17] = 1 - x183*x106 + x183*x107 - x45*x195 + x54*x195;

    }
    else if (x172)
    {
        DEBUG_PRINT("Case (x172) is True.\n");
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x105=x19*x101; POST_CHECK_MUL(x105);
        x106=x46*(rt2 + x105); POST_CHECK_MUL(x106);
        x107=(rt2 - x105)*x55; POST_CHECK_MUL(x107);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x184=random2*x118; POST_CHECK_MUL(x184);

        /* Assignment result[2, 5]=1 - x184*x106 + x184*x107 */
        /*@ assert (x30 + x31) >= 0.;*/
        x45=sqrt(x30 + x31); POST_CHECK_POW(x45);
        /*@ assert (x45) != 0.;*/
        x46=0.5*1.0/x45; POST_CHECK_MUL(x46);
        /*@ assert (x22) != 0.;*/
        x51=1.0/x22; POST_CHECK_POW(x51);
        /*@ assert (x32) >= 0.;*/
        x54=sqrt(x32); POST_CHECK_POW(x54);
        /*@ assert (x54) != 0.;*/
        x55=0.5*1.0/x54; POST_CHECK_MUL(x55);
        x101=x11*x51; POST_CHECK_MUL(x101);
        x105=x19*x101; POST_CHECK_MUL(x105);
        x106=x46*(rt2 + x105); POST_CHECK_MUL(x106);
        x107=(rt2 - x105)*x55; POST_CHECK_MUL(x107);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x118=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x118);
        x184=random2*x118; POST_CHECK_MUL(x184);
        result[17] = 1 - x184*x106 + x184*x107;

    }
    /*@ assert (result[17]) >= 0.;*/
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

    /* Assignment result=0.5*(rt1 + x13*x19 - x14*x19 + x7)**2 + 0.5*(rt2 + x13*x20 - x14*x20 + x9)**2 + 0.5*(x4 - x13 - x14 + x3)**2 */
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
    double x13;
    double x14;
    double x15;
    int x16;
    double x17;
    int x18;
    double x19;
    double x20;x1=ut1*ut1; POST_CHECK_POW(x1);
    /*@ assert (x1) >= 0.;*/
    x2=ut2*ut2; POST_CHECK_POW(x2);
    /*@ assert (x2) >= 0.;*/
    /*@ assert (x1 + x2) >= 0.;*/
    x3=mu*sqrt(x1 + x2) + un; POST_CHECK_ADD(x3);
    x4=mu*rn; POST_CHECK_MUL(x4);
    /*@ assert (x4) >= 0.;*/
    /*@ assert (x4) != 0.;*/
    x5=mu*mu; POST_CHECK_POW(x5);
    x6=rn*rn*x5 + rt1*rt1 + rt2*rt2 + x5*x1 + x5*x2 + x3*x3; POST_CHECK_ADD(x6);
    x7=mu*ut1; POST_CHECK_MUL(x7);
    x8=rt1*x4 + x7*x3; POST_CHECK_ADD(x8);
    x9=mu*ut2; POST_CHECK_MUL(x9);
    x10=rt2*x4 + x9*x3; POST_CHECK_ADD(x10);
    /*@ assert (x10*x10 + x8*x8) >= 0.;*/
    x11=sqrt(x10*x10 + x8*x8); POST_CHECK_POW(x11);
    x12=2*x11; POST_CHECK_MUL(x12);
    /*@ assert (-x12 + x6) >= 0.;*/
    x13=0.5*sqrt(-x12 + x6); POST_CHECK_MUL(x13);
    /*@ assert (x12 + x6) >= 0.;*/
    x14=0.5*sqrt(x12 + x6); POST_CHECK_MUL(x14);
    /*@ assert (x11) != 0.;*/
    x15=1.0/x11; POST_CHECK_POW(x15);
    x16=x11 > 0; POST_CHECK(x16);
    /*@ assert (random1*random1 + random2*random2) >= 0.;*/
    /*@ assert (random1*random1 + random2*random2) != 0.;*/
    x17=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x17);
    x18=x11 <= ZERO; POST_CHECK(x18);
    x19=((x16) ? (x8*x15): (random1*x17)); POST_CHECK(x19);
    x20=((x16) ? (x10*x15): (random2*x17)); POST_CHECK(x20);
    *result = 0.5*(rt1 + x13*x19 - x14*x19 + x7)*(rt1 + x13*x19 - x14*x19 + x7) + 0.5*(rt2 + x13*x20 - x14*x20 + x9)*(rt2 + x13*x20 - x14*x20 + x9) + 0.5*(x4 - x13 - x14 + x3)*(x4 - x13 - x14 + x3);
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
    double x13;
    int x46;
    double x4;
    double x6;
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
    double x45;
    x1=ut1*ut1; POST_CHECK_POW(x1);
    /*@ assert (x1) >= 0.;*/
    x2=ut2*ut2; POST_CHECK_POW(x2);
    /*@ assert (x2) >= 0.;*/
    /*@ assert (x1 + x2) >= 0.;*/
    x3=sqrt(x1 + x2); POST_CHECK_POW(x3);
    /*@ assert (x3) >= 0.;*/
    x5=mu*x3 + un; POST_CHECK_ADD(x5);
    x7=mu*rn; POST_CHECK_MUL(x7);
    x8=mu*ut1*x5 + rt1*x7; POST_CHECK_ADD(x8);
    x9=x8*x8; POST_CHECK_POW(x9);
    x10=mu*ut2*x5 + rt2*x7; POST_CHECK_ADD(x10);
    x11=x10*x10; POST_CHECK_POW(x11);
    x12=x11 + x9; POST_CHECK_ADD(x12);
    /*@ assert (x12) >= 0.;*/
    x13=sqrt(x12); POST_CHECK_POW(x13);
    x46=x13 > 0; POST_CHECK(x46);
    int x54;
    double x47;
    double x48;
    double x49;
    double x50;
    double x51;
    double x52;
    double x53;
    x54=x13 <= ZERO; POST_CHECK(x54);
    if (x46)
    {
        x4=mu*mu; POST_CHECK_POW(x4);
        x6=x4*rn*rn + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5; POST_CHECK_ADD(x6);
        x14=2*x13; POST_CHECK_MUL(x14);
        /*@ assert (-x14 + x6) >= 0.;*/
        x15=sqrt(-x14 + x6); POST_CHECK_POW(x15);
        /*@ assert (x14 + x6) >= 0.;*/
        x16=sqrt(x14 + x6); POST_CHECK_POW(x16);
        x17=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x15 - 0.25*x16; POST_CHECK_ADD(x17);
        /*@ assert (x16) != 0.;*/
        x18=1.0*1.0/x16; POST_CHECK_MUL(x18);
        x19=rn*x4; POST_CHECK_MUL(x19);
        x20=mu*rt1; POST_CHECK_MUL(x20);
        x21=x20*x8; POST_CHECK_MUL(x21);
        x22=mu*rt2; POST_CHECK_MUL(x22);
        x23=x22*x10; POST_CHECK_MUL(x23);
        /*@ assert (x13) != 0.;*/
        x24=1.0/x13; POST_CHECK_POW(x24);
        x25=x24*(x21 + x23); POST_CHECK_MUL(x25);
        x26=x18*(x19 + x25); POST_CHECK_MUL(x26);
        /*@ assert (x15) != 0.;*/
        x27=1.0/x15; POST_CHECK_POW(x27);
        x28=1.0*x27; POST_CHECK_MUL(x28);
        x29=x19 - x25; POST_CHECK_ADD(x29);
        x30=x17*(2*mu - x26 - x28*x29); POST_CHECK_MUL(x30);
        x31=0.5*rt1; POST_CHECK_MUL(x31);
        x32=0.5*mu*ut1; POST_CHECK_MUL(x32);
        x33=0.25*x24*x8*x15 - 0.25*x24*x8*x16 + x31 + x32; POST_CHECK_ADD(x33);
        /*@ assert (x12) >= 0.;*/
        /*@ assert (x12) != 0.;*/
        x34=pow(x12, -3.0/2.0); POST_CHECK_POW(x34);
        x35=x34*(-x21 - x23); POST_CHECK_MUL(x35);
        x36=x24*x20 + x8*x35; POST_CHECK_ADD(x36);
        x37=1.0*x16; POST_CHECK_MUL(x37);
        x38=x24*x8*x27; POST_CHECK_MUL(x38);
        x39=x24*x8; POST_CHECK_MUL(x39);
        x40=0.5*rt2; POST_CHECK_MUL(x40);
        x41=0.5*mu*ut2; POST_CHECK_MUL(x41);
        x42=0.25*x24*x10*x15 - 0.25*x24*x10*x16 + x40 + x41; POST_CHECK_ADD(x42);
        x43=x10*x35 + x24*x22; POST_CHECK_ADD(x43);
        x44=x24*x10*x27; POST_CHECK_MUL(x44);
        x45=x24*x10; POST_CHECK_MUL(x45);
    }
    else if (x54)
    {
        x4=mu*mu; POST_CHECK_POW(x4);
        x6=x4*rn*rn + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5; POST_CHECK_ADD(x6);
        x14=2*x13; POST_CHECK_MUL(x14);
        /*@ assert (-x14 + x6) >= 0.;*/
        x15=sqrt(-x14 + x6); POST_CHECK_POW(x15);
        /*@ assert (x14 + x6) >= 0.;*/
        x16=sqrt(x14 + x6); POST_CHECK_POW(x16);
        x17=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x15 - 0.25*x16; POST_CHECK_ADD(x17);
        /*@ assert (x16) != 0.;*/
        x18=1.0*1.0/x16; POST_CHECK_MUL(x18);
        x19=rn*x4; POST_CHECK_MUL(x19);
        x20=mu*rt1; POST_CHECK_MUL(x20);
        x21=x20*x8; POST_CHECK_MUL(x21);
        x22=mu*rt2; POST_CHECK_MUL(x22);
        x23=x22*x10; POST_CHECK_MUL(x23);
        /*@ assert (x13) != 0.;*/
        x24=1.0/x13; POST_CHECK_POW(x24);
        x25=x24*(x21 + x23); POST_CHECK_MUL(x25);
        x26=x18*(x19 + x25); POST_CHECK_MUL(x26);
        /*@ assert (x15) != 0.;*/
        x27=1.0/x15; POST_CHECK_POW(x27);
        x28=1.0*x27; POST_CHECK_MUL(x28);
        x29=x19 - x25; POST_CHECK_ADD(x29);
        x30=x17*(2*mu - x26 - x28*x29); POST_CHECK_MUL(x30);
        x31=0.5*rt1; POST_CHECK_MUL(x31);
        x32=0.5*mu*ut1; POST_CHECK_MUL(x32);
        x40=0.5*rt2; POST_CHECK_MUL(x40);
        x41=0.5*mu*ut2; POST_CHECK_MUL(x41);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x47=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x47);
        x48=0.25*random1*x47*x15 - 0.25*random1*x47*x16 + x31 + x32; POST_CHECK_ADD(x48);
        x49=random1*x47*x27; POST_CHECK_MUL(x49);
        x50=random1*x47; POST_CHECK_MUL(x50);
        x51=0.25*random2*x47*x15 - 0.25*random2*x47*x16 + x40 + x41; POST_CHECK_ADD(x51);
        x52=random2*x47*x27; POST_CHECK_MUL(x52);
        x53=random2*x47; POST_CHECK_MUL(x53);
    }
    /* Assignment result[0, 0]=Piecewise((x30 + x33*(x15*x36 - x26*x39 + x29*x38 - x37*x36) + x42*(x15*x43 - x26*x45 + x29*x44 - x37*x43), x46), (x30 + x48*(-x26*x50 + x29*x49) + x51*(-x26*x53 + x29*x52), x54)) */

    if (x46)
    {
        DEBUG_PRINT("Case (x46) is True.\n");
        x4=mu*mu; POST_CHECK_POW(x4);
        x6=x4*rn*rn + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5; POST_CHECK_ADD(x6);
        x14=2*x13; POST_CHECK_MUL(x14);
        /*@ assert (-x14 + x6) >= 0.;*/
        x15=sqrt(-x14 + x6); POST_CHECK_POW(x15);
        /*@ assert (x14 + x6) >= 0.;*/
        x16=sqrt(x14 + x6); POST_CHECK_POW(x16);
        x17=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x15 - 0.25*x16; POST_CHECK_ADD(x17);
        /*@ assert (x16) != 0.;*/
        x18=1.0*1.0/x16; POST_CHECK_MUL(x18);
        x19=rn*x4; POST_CHECK_MUL(x19);
        x20=mu*rt1; POST_CHECK_MUL(x20);
        x21=x20*x8; POST_CHECK_MUL(x21);
        x22=mu*rt2; POST_CHECK_MUL(x22);
        x23=x22*x10; POST_CHECK_MUL(x23);
        /*@ assert (x13) != 0.;*/
        x24=1.0/x13; POST_CHECK_POW(x24);
        x25=x24*(x21 + x23); POST_CHECK_MUL(x25);
        x26=x18*(x19 + x25); POST_CHECK_MUL(x26);
        /*@ assert (x15) != 0.;*/
        x27=1.0/x15; POST_CHECK_POW(x27);
        x28=1.0*x27; POST_CHECK_MUL(x28);
        x29=x19 - x25; POST_CHECK_ADD(x29);
        x30=x17*(2*mu - x26 - x28*x29); POST_CHECK_MUL(x30);
        x31=0.5*rt1; POST_CHECK_MUL(x31);
        x32=0.5*mu*ut1; POST_CHECK_MUL(x32);
        x33=0.25*x24*x8*x15 - 0.25*x24*x8*x16 + x31 + x32; POST_CHECK_ADD(x33);
        /*@ assert (x12) >= 0.;*/
        /*@ assert (x12) != 0.;*/
        x34=pow(x12, -3.0/2.0); POST_CHECK_POW(x34);
        x35=x34*(-x21 - x23); POST_CHECK_MUL(x35);
        x36=x24*x20 + x8*x35; POST_CHECK_ADD(x36);
        x37=1.0*x16; POST_CHECK_MUL(x37);
        x38=x24*x8*x27; POST_CHECK_MUL(x38);
        x39=x24*x8; POST_CHECK_MUL(x39);
        x40=0.5*rt2; POST_CHECK_MUL(x40);
        x41=0.5*mu*ut2; POST_CHECK_MUL(x41);
        x42=0.25*x24*x10*x15 - 0.25*x24*x10*x16 + x40 + x41; POST_CHECK_ADD(x42);
        x43=x10*x35 + x24*x22; POST_CHECK_ADD(x43);
        x44=x24*x10*x27; POST_CHECK_MUL(x44);
        x45=x24*x10; POST_CHECK_MUL(x45);

        /* Assignment result[0, 0]=x30 + x33*(x15*x36 - x26*x39 + x29*x38 - x37*x36) + x42*(x15*x43 - x26*x45 + x29*x44 - x37*x43) */
        x4=mu*mu; POST_CHECK_POW(x4);
        x6=x4*rn*rn + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5; POST_CHECK_ADD(x6);
        x14=2*x13; POST_CHECK_MUL(x14);
        /*@ assert (-x14 + x6) >= 0.;*/
        x15=sqrt(-x14 + x6); POST_CHECK_POW(x15);
        /*@ assert (x14 + x6) >= 0.;*/
        x16=sqrt(x14 + x6); POST_CHECK_POW(x16);
        x17=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x15 - 0.25*x16; POST_CHECK_ADD(x17);
        /*@ assert (x16) != 0.;*/
        x18=1.0*1.0/x16; POST_CHECK_MUL(x18);
        x19=rn*x4; POST_CHECK_MUL(x19);
        x20=mu*rt1; POST_CHECK_MUL(x20);
        x21=x20*x8; POST_CHECK_MUL(x21);
        x22=mu*rt2; POST_CHECK_MUL(x22);
        x23=x22*x10; POST_CHECK_MUL(x23);
        /*@ assert (x13) != 0.;*/
        x24=1.0/x13; POST_CHECK_POW(x24);
        x25=x24*(x21 + x23); POST_CHECK_MUL(x25);
        x26=x18*(x19 + x25); POST_CHECK_MUL(x26);
        /*@ assert (x15) != 0.;*/
        x27=1.0/x15; POST_CHECK_POW(x27);
        x28=1.0*x27; POST_CHECK_MUL(x28);
        x29=x19 - x25; POST_CHECK_ADD(x29);
        x30=x17*(2*mu - x26 - x28*x29); POST_CHECK_MUL(x30);
        x31=0.5*rt1; POST_CHECK_MUL(x31);
        x32=0.5*mu*ut1; POST_CHECK_MUL(x32);
        x33=0.25*x24*x8*x15 - 0.25*x24*x8*x16 + x31 + x32; POST_CHECK_ADD(x33);
        /*@ assert (x12) >= 0.;*/
        /*@ assert (x12) != 0.;*/
        x34=pow(x12, -3.0/2.0); POST_CHECK_POW(x34);
        x35=x34*(-x21 - x23); POST_CHECK_MUL(x35);
        x36=x24*x20 + x8*x35; POST_CHECK_ADD(x36);
        x37=1.0*x16; POST_CHECK_MUL(x37);
        x38=x24*x8*x27; POST_CHECK_MUL(x38);
        x39=x24*x8; POST_CHECK_MUL(x39);
        x40=0.5*rt2; POST_CHECK_MUL(x40);
        x41=0.5*mu*ut2; POST_CHECK_MUL(x41);
        x42=0.25*x24*x10*x15 - 0.25*x24*x10*x16 + x40 + x41; POST_CHECK_ADD(x42);
        x43=x10*x35 + x24*x22; POST_CHECK_ADD(x43);
        x44=x24*x10*x27; POST_CHECK_MUL(x44);
        x45=x24*x10; POST_CHECK_MUL(x45);
        result[0] = x30 + x33*(x15*x36 - x26*x39 + x29*x38 - x37*x36) + x42*(x15*x43 - x26*x45 + x29*x44 - x37*x43);

    }
    else if (x54)
    {
        DEBUG_PRINT("Case (x54) is True.\n");
        x4=mu*mu; POST_CHECK_POW(x4);
        x6=x4*rn*rn + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5; POST_CHECK_ADD(x6);
        x14=2*x13; POST_CHECK_MUL(x14);
        /*@ assert (-x14 + x6) >= 0.;*/
        x15=sqrt(-x14 + x6); POST_CHECK_POW(x15);
        /*@ assert (x14 + x6) >= 0.;*/
        x16=sqrt(x14 + x6); POST_CHECK_POW(x16);
        x17=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x15 - 0.25*x16; POST_CHECK_ADD(x17);
        /*@ assert (x16) != 0.;*/
        x18=1.0*1.0/x16; POST_CHECK_MUL(x18);
        x19=rn*x4; POST_CHECK_MUL(x19);
        x20=mu*rt1; POST_CHECK_MUL(x20);
        x21=x20*x8; POST_CHECK_MUL(x21);
        x22=mu*rt2; POST_CHECK_MUL(x22);
        x23=x22*x10; POST_CHECK_MUL(x23);
        /*@ assert (x13) != 0.;*/
        x24=1.0/x13; POST_CHECK_POW(x24);
        x25=x24*(x21 + x23); POST_CHECK_MUL(x25);
        x26=x18*(x19 + x25); POST_CHECK_MUL(x26);
        /*@ assert (x15) != 0.;*/
        x27=1.0/x15; POST_CHECK_POW(x27);
        x28=1.0*x27; POST_CHECK_MUL(x28);
        x29=x19 - x25; POST_CHECK_ADD(x29);
        x30=x17*(2*mu - x26 - x28*x29); POST_CHECK_MUL(x30);
        x31=0.5*rt1; POST_CHECK_MUL(x31);
        x32=0.5*mu*ut1; POST_CHECK_MUL(x32);
        x40=0.5*rt2; POST_CHECK_MUL(x40);
        x41=0.5*mu*ut2; POST_CHECK_MUL(x41);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x47=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x47);
        x48=0.25*random1*x47*x15 - 0.25*random1*x47*x16 + x31 + x32; POST_CHECK_ADD(x48);
        x49=random1*x47*x27; POST_CHECK_MUL(x49);
        x50=random1*x47; POST_CHECK_MUL(x50);
        x51=0.25*random2*x47*x15 - 0.25*random2*x47*x16 + x40 + x41; POST_CHECK_ADD(x51);
        x52=random2*x47*x27; POST_CHECK_MUL(x52);
        x53=random2*x47; POST_CHECK_MUL(x53);

        /* Assignment result[0, 0]=x30 + x48*(-x26*x50 + x29*x49) + x51*(-x26*x53 + x29*x52) */
        x4=mu*mu; POST_CHECK_POW(x4);
        x6=x4*rn*rn + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5; POST_CHECK_ADD(x6);
        x14=2*x13; POST_CHECK_MUL(x14);
        /*@ assert (-x14 + x6) >= 0.;*/
        x15=sqrt(-x14 + x6); POST_CHECK_POW(x15);
        /*@ assert (x14 + x6) >= 0.;*/
        x16=sqrt(x14 + x6); POST_CHECK_POW(x16);
        x17=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x15 - 0.25*x16; POST_CHECK_ADD(x17);
        /*@ assert (x16) != 0.;*/
        x18=1.0*1.0/x16; POST_CHECK_MUL(x18);
        x19=rn*x4; POST_CHECK_MUL(x19);
        x20=mu*rt1; POST_CHECK_MUL(x20);
        x21=x20*x8; POST_CHECK_MUL(x21);
        x22=mu*rt2; POST_CHECK_MUL(x22);
        x23=x22*x10; POST_CHECK_MUL(x23);
        /*@ assert (x13) != 0.;*/
        x24=1.0/x13; POST_CHECK_POW(x24);
        x25=x24*(x21 + x23); POST_CHECK_MUL(x25);
        x26=x18*(x19 + x25); POST_CHECK_MUL(x26);
        /*@ assert (x15) != 0.;*/
        x27=1.0/x15; POST_CHECK_POW(x27);
        x28=1.0*x27; POST_CHECK_MUL(x28);
        x29=x19 - x25; POST_CHECK_ADD(x29);
        x30=x17*(2*mu - x26 - x28*x29); POST_CHECK_MUL(x30);
        x31=0.5*rt1; POST_CHECK_MUL(x31);
        x32=0.5*mu*ut1; POST_CHECK_MUL(x32);
        x40=0.5*rt2; POST_CHECK_MUL(x40);
        x41=0.5*mu*ut2; POST_CHECK_MUL(x41);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x47=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x47);
        x48=0.25*random1*x47*x15 - 0.25*random1*x47*x16 + x31 + x32; POST_CHECK_ADD(x48);
        x49=random1*x47*x27; POST_CHECK_MUL(x49);
        x50=random1*x47; POST_CHECK_MUL(x50);
        x51=0.25*random2*x47*x15 - 0.25*random2*x47*x16 + x40 + x41; POST_CHECK_ADD(x51);
        x52=random2*x47*x27; POST_CHECK_MUL(x52);
        x53=random2*x47; POST_CHECK_MUL(x53);
        result[0] = x30 + x48*(-x26*x50 + x29*x49) + x51*(-x26*x53 + x29*x52);

    }


    /* Assignment result[0, 1]=Piecewise((x33*(2 + x15*x63 - x37*x63 - x39*x57 + x58*x38) + x42*(-x45*x57 + x58*x44 + x61) + x59, x46), (x48*(2 - x50*x57 + x58*x49) + x51*(-x53*x57 + x58*x52) + x59, x54)) */
    double x55;
    double x56;
    double x57;
    double x58;
    double x59;
    double x60;
    double x61;
    double x62;
    double x63;
    if (x46)
    {
        DEBUG_PRINT("Case (x46) is True.\n");
        x4=mu*mu; POST_CHECK_POW(x4);
        x6=x4*rn*rn + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5; POST_CHECK_ADD(x6);
        x14=2*x13; POST_CHECK_MUL(x14);
        /*@ assert (-x14 + x6) >= 0.;*/
        x15=sqrt(-x14 + x6); POST_CHECK_POW(x15);
        /*@ assert (x14 + x6) >= 0.;*/
        x16=sqrt(x14 + x6); POST_CHECK_POW(x16);
        x17=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x15 - 0.25*x16; POST_CHECK_ADD(x17);
        /*@ assert (x16) != 0.;*/
        x18=1.0*1.0/x16; POST_CHECK_MUL(x18);
        /*@ assert (x13) != 0.;*/
        x24=1.0/x13; POST_CHECK_POW(x24);
        /*@ assert (x15) != 0.;*/
        x27=1.0/x15; POST_CHECK_POW(x27);
        x28=1.0*x27; POST_CHECK_MUL(x28);
        x31=0.5*rt1; POST_CHECK_MUL(x31);
        x32=0.5*mu*ut1; POST_CHECK_MUL(x32);
        x33=0.25*x24*x8*x15 - 0.25*x24*x8*x16 + x31 + x32; POST_CHECK_ADD(x33);
        /*@ assert (x12) >= 0.;*/
        /*@ assert (x12) != 0.;*/
        x34=pow(x12, -3.0/2.0); POST_CHECK_POW(x34);
        x37=1.0*x16; POST_CHECK_MUL(x37);
        x38=x24*x8*x27; POST_CHECK_MUL(x38);
        x39=x24*x8; POST_CHECK_MUL(x39);
        x40=0.5*rt2; POST_CHECK_MUL(x40);
        x41=0.5*mu*ut2; POST_CHECK_MUL(x41);
        x42=0.25*x24*x10*x15 - 0.25*x24*x10*x16 + x40 + x41; POST_CHECK_ADD(x42);
        x44=x24*x10*x27; POST_CHECK_MUL(x44);
        x45=x24*x10; POST_CHECK_MUL(x45);
        x55=x24*x7; POST_CHECK_MUL(x55);
        x56=x8*x55; POST_CHECK_MUL(x56);
        x57=(rt1 + x56)*x18; POST_CHECK_MUL(x57);
        x58=rt1 - x56; POST_CHECK_ADD(x58);
        x59=x17*(-x57 - x58*x28); POST_CHECK_MUL(x59);
        x60=x34*mu*rn*x8*x10; POST_CHECK_MUL(x60);
        x61=-x60*x15 + x60*x37; POST_CHECK_ADD(x61);
        x62=x34*mu*rn; POST_CHECK_MUL(x62);
        x63=x55 - x9*x62; POST_CHECK_ADD(x63);

        /* Assignment result[0, 1]=x33*(2 + x15*x63 - x37*x63 - x39*x57 + x58*x38) + x42*(-x45*x57 + x58*x44 + x61) + x59 */
        x4=mu*mu; POST_CHECK_POW(x4);
        x6=x4*rn*rn + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5; POST_CHECK_ADD(x6);
        x14=2*x13; POST_CHECK_MUL(x14);
        /*@ assert (-x14 + x6) >= 0.;*/
        x15=sqrt(-x14 + x6); POST_CHECK_POW(x15);
        /*@ assert (x14 + x6) >= 0.;*/
        x16=sqrt(x14 + x6); POST_CHECK_POW(x16);
        x17=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x15 - 0.25*x16; POST_CHECK_ADD(x17);
        /*@ assert (x16) != 0.;*/
        x18=1.0*1.0/x16; POST_CHECK_MUL(x18);
        /*@ assert (x13) != 0.;*/
        x24=1.0/x13; POST_CHECK_POW(x24);
        /*@ assert (x15) != 0.;*/
        x27=1.0/x15; POST_CHECK_POW(x27);
        x28=1.0*x27; POST_CHECK_MUL(x28);
        x31=0.5*rt1; POST_CHECK_MUL(x31);
        x32=0.5*mu*ut1; POST_CHECK_MUL(x32);
        x33=0.25*x24*x8*x15 - 0.25*x24*x8*x16 + x31 + x32; POST_CHECK_ADD(x33);
        /*@ assert (x12) >= 0.;*/
        /*@ assert (x12) != 0.;*/
        x34=pow(x12, -3.0/2.0); POST_CHECK_POW(x34);
        x37=1.0*x16; POST_CHECK_MUL(x37);
        x38=x24*x8*x27; POST_CHECK_MUL(x38);
        x39=x24*x8; POST_CHECK_MUL(x39);
        x40=0.5*rt2; POST_CHECK_MUL(x40);
        x41=0.5*mu*ut2; POST_CHECK_MUL(x41);
        x42=0.25*x24*x10*x15 - 0.25*x24*x10*x16 + x40 + x41; POST_CHECK_ADD(x42);
        x44=x24*x10*x27; POST_CHECK_MUL(x44);
        x45=x24*x10; POST_CHECK_MUL(x45);
        x55=x24*x7; POST_CHECK_MUL(x55);
        x56=x8*x55; POST_CHECK_MUL(x56);
        x57=(rt1 + x56)*x18; POST_CHECK_MUL(x57);
        x58=rt1 - x56; POST_CHECK_ADD(x58);
        x59=x17*(-x57 - x58*x28); POST_CHECK_MUL(x59);
        x60=x34*mu*rn*x8*x10; POST_CHECK_MUL(x60);
        x61=-x60*x15 + x60*x37; POST_CHECK_ADD(x61);
        x62=x34*mu*rn; POST_CHECK_MUL(x62);
        x63=x55 - x9*x62; POST_CHECK_ADD(x63);
        result[1] = x33*(2 + x15*x63 - x37*x63 - x39*x57 + x58*x38) + x42*(-x45*x57 + x58*x44 + x61) + x59;

    }
    else if (x54)
    {
        DEBUG_PRINT("Case (x54) is True.\n");
        x4=mu*mu; POST_CHECK_POW(x4);
        x6=x4*rn*rn + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5; POST_CHECK_ADD(x6);
        x14=2*x13; POST_CHECK_MUL(x14);
        /*@ assert (-x14 + x6) >= 0.;*/
        x15=sqrt(-x14 + x6); POST_CHECK_POW(x15);
        /*@ assert (x14 + x6) >= 0.;*/
        x16=sqrt(x14 + x6); POST_CHECK_POW(x16);
        x17=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x15 - 0.25*x16; POST_CHECK_ADD(x17);
        /*@ assert (x16) != 0.;*/
        x18=1.0*1.0/x16; POST_CHECK_MUL(x18);
        /*@ assert (x13) != 0.;*/
        x24=1.0/x13; POST_CHECK_POW(x24);
        /*@ assert (x15) != 0.;*/
        x27=1.0/x15; POST_CHECK_POW(x27);
        x28=1.0*x27; POST_CHECK_MUL(x28);
        x31=0.5*rt1; POST_CHECK_MUL(x31);
        x32=0.5*mu*ut1; POST_CHECK_MUL(x32);
        x40=0.5*rt2; POST_CHECK_MUL(x40);
        x41=0.5*mu*ut2; POST_CHECK_MUL(x41);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x47=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x47);
        x48=0.25*random1*x47*x15 - 0.25*random1*x47*x16 + x31 + x32; POST_CHECK_ADD(x48);
        x49=random1*x47*x27; POST_CHECK_MUL(x49);
        x50=random1*x47; POST_CHECK_MUL(x50);
        x51=0.25*random2*x47*x15 - 0.25*random2*x47*x16 + x40 + x41; POST_CHECK_ADD(x51);
        x52=random2*x47*x27; POST_CHECK_MUL(x52);
        x53=random2*x47; POST_CHECK_MUL(x53);
        x55=x24*x7; POST_CHECK_MUL(x55);
        x56=x8*x55; POST_CHECK_MUL(x56);
        x57=(rt1 + x56)*x18; POST_CHECK_MUL(x57);
        x58=rt1 - x56; POST_CHECK_ADD(x58);
        x59=x17*(-x57 - x58*x28); POST_CHECK_MUL(x59);

        /* Assignment result[0, 1]=x48*(2 - x50*x57 + x58*x49) + x51*(-x53*x57 + x58*x52) + x59 */
        x4=mu*mu; POST_CHECK_POW(x4);
        x6=x4*rn*rn + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5; POST_CHECK_ADD(x6);
        x14=2*x13; POST_CHECK_MUL(x14);
        /*@ assert (-x14 + x6) >= 0.;*/
        x15=sqrt(-x14 + x6); POST_CHECK_POW(x15);
        /*@ assert (x14 + x6) >= 0.;*/
        x16=sqrt(x14 + x6); POST_CHECK_POW(x16);
        x17=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x15 - 0.25*x16; POST_CHECK_ADD(x17);
        /*@ assert (x16) != 0.;*/
        x18=1.0*1.0/x16; POST_CHECK_MUL(x18);
        /*@ assert (x13) != 0.;*/
        x24=1.0/x13; POST_CHECK_POW(x24);
        /*@ assert (x15) != 0.;*/
        x27=1.0/x15; POST_CHECK_POW(x27);
        x28=1.0*x27; POST_CHECK_MUL(x28);
        x31=0.5*rt1; POST_CHECK_MUL(x31);
        x32=0.5*mu*ut1; POST_CHECK_MUL(x32);
        x40=0.5*rt2; POST_CHECK_MUL(x40);
        x41=0.5*mu*ut2; POST_CHECK_MUL(x41);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x47=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x47);
        x48=0.25*random1*x47*x15 - 0.25*random1*x47*x16 + x31 + x32; POST_CHECK_ADD(x48);
        x49=random1*x47*x27; POST_CHECK_MUL(x49);
        x50=random1*x47; POST_CHECK_MUL(x50);
        x51=0.25*random2*x47*x15 - 0.25*random2*x47*x16 + x40 + x41; POST_CHECK_ADD(x51);
        x52=random2*x47*x27; POST_CHECK_MUL(x52);
        x53=random2*x47; POST_CHECK_MUL(x53);
        x55=x24*x7; POST_CHECK_MUL(x55);
        x56=x8*x55; POST_CHECK_MUL(x56);
        x57=(rt1 + x56)*x18; POST_CHECK_MUL(x57);
        x58=rt1 - x56; POST_CHECK_ADD(x58);
        x59=x17*(-x57 - x58*x28); POST_CHECK_MUL(x59);
        result[1] = x48*(2 - x50*x57 + x58*x49) + x51*(-x53*x57 + x58*x52) + x59;

    }


    /* Assignment result[0, 2]=Piecewise((x33*(-x39*x65 + x61 + x66*x38) + x42*(2 + x15*x68 - x37*x68 - x45*x65 + x66*x44) + x67, x46), (x48*(-x50*x65 + x66*x49) + x51*(2 - x53*x65 + x66*x52) + x67, x54)) */
    double x64;
    double x65;
    double x66;
    double x67;
    double x68;
    if (x46)
    {
        DEBUG_PRINT("Case (x46) is True.\n");
        x4=mu*mu; POST_CHECK_POW(x4);
        x6=x4*rn*rn + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5; POST_CHECK_ADD(x6);
        x14=2*x13; POST_CHECK_MUL(x14);
        /*@ assert (-x14 + x6) >= 0.;*/
        x15=sqrt(-x14 + x6); POST_CHECK_POW(x15);
        /*@ assert (x14 + x6) >= 0.;*/
        x16=sqrt(x14 + x6); POST_CHECK_POW(x16);
        x17=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x15 - 0.25*x16; POST_CHECK_ADD(x17);
        /*@ assert (x16) != 0.;*/
        x18=1.0*1.0/x16; POST_CHECK_MUL(x18);
        /*@ assert (x13) != 0.;*/
        x24=1.0/x13; POST_CHECK_POW(x24);
        /*@ assert (x15) != 0.;*/
        x27=1.0/x15; POST_CHECK_POW(x27);
        x28=1.0*x27; POST_CHECK_MUL(x28);
        x31=0.5*rt1; POST_CHECK_MUL(x31);
        x32=0.5*mu*ut1; POST_CHECK_MUL(x32);
        x33=0.25*x24*x8*x15 - 0.25*x24*x8*x16 + x31 + x32; POST_CHECK_ADD(x33);
        /*@ assert (x12) >= 0.;*/
        /*@ assert (x12) != 0.;*/
        x34=pow(x12, -3.0/2.0); POST_CHECK_POW(x34);
        x37=1.0*x16; POST_CHECK_MUL(x37);
        x38=x24*x8*x27; POST_CHECK_MUL(x38);
        x39=x24*x8; POST_CHECK_MUL(x39);
        x40=0.5*rt2; POST_CHECK_MUL(x40);
        x41=0.5*mu*ut2; POST_CHECK_MUL(x41);
        x42=0.25*x24*x10*x15 - 0.25*x24*x10*x16 + x40 + x41; POST_CHECK_ADD(x42);
        x44=x24*x10*x27; POST_CHECK_MUL(x44);
        x45=x24*x10; POST_CHECK_MUL(x45);
        x55=x24*x7; POST_CHECK_MUL(x55);
        x60=x34*mu*rn*x8*x10; POST_CHECK_MUL(x60);
        x61=-x60*x15 + x60*x37; POST_CHECK_ADD(x61);
        x62=x34*mu*rn; POST_CHECK_MUL(x62);
        x64=x10*x55; POST_CHECK_MUL(x64);
        x65=(rt2 + x64)*x18; POST_CHECK_MUL(x65);
        x66=rt2 - x64; POST_CHECK_ADD(x66);
        x67=x17*(-x65 - x66*x28); POST_CHECK_MUL(x67);
        x68=-x11*x62 + x55; POST_CHECK_ADD(x68);

        /* Assignment result[0, 2]=x33*(-x39*x65 + x61 + x66*x38) + x42*(2 + x15*x68 - x37*x68 - x45*x65 + x66*x44) + x67 */
        x4=mu*mu; POST_CHECK_POW(x4);
        x6=x4*rn*rn + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5; POST_CHECK_ADD(x6);
        x14=2*x13; POST_CHECK_MUL(x14);
        /*@ assert (-x14 + x6) >= 0.;*/
        x15=sqrt(-x14 + x6); POST_CHECK_POW(x15);
        /*@ assert (x14 + x6) >= 0.;*/
        x16=sqrt(x14 + x6); POST_CHECK_POW(x16);
        x17=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x15 - 0.25*x16; POST_CHECK_ADD(x17);
        /*@ assert (x16) != 0.;*/
        x18=1.0*1.0/x16; POST_CHECK_MUL(x18);
        /*@ assert (x13) != 0.;*/
        x24=1.0/x13; POST_CHECK_POW(x24);
        /*@ assert (x15) != 0.;*/
        x27=1.0/x15; POST_CHECK_POW(x27);
        x28=1.0*x27; POST_CHECK_MUL(x28);
        x31=0.5*rt1; POST_CHECK_MUL(x31);
        x32=0.5*mu*ut1; POST_CHECK_MUL(x32);
        x33=0.25*x24*x8*x15 - 0.25*x24*x8*x16 + x31 + x32; POST_CHECK_ADD(x33);
        /*@ assert (x12) >= 0.;*/
        /*@ assert (x12) != 0.;*/
        x34=pow(x12, -3.0/2.0); POST_CHECK_POW(x34);
        x37=1.0*x16; POST_CHECK_MUL(x37);
        x38=x24*x8*x27; POST_CHECK_MUL(x38);
        x39=x24*x8; POST_CHECK_MUL(x39);
        x40=0.5*rt2; POST_CHECK_MUL(x40);
        x41=0.5*mu*ut2; POST_CHECK_MUL(x41);
        x42=0.25*x24*x10*x15 - 0.25*x24*x10*x16 + x40 + x41; POST_CHECK_ADD(x42);
        x44=x24*x10*x27; POST_CHECK_MUL(x44);
        x45=x24*x10; POST_CHECK_MUL(x45);
        x55=x24*x7; POST_CHECK_MUL(x55);
        x60=x34*mu*rn*x8*x10; POST_CHECK_MUL(x60);
        x61=-x60*x15 + x60*x37; POST_CHECK_ADD(x61);
        x62=x34*mu*rn; POST_CHECK_MUL(x62);
        x64=x10*x55; POST_CHECK_MUL(x64);
        x65=(rt2 + x64)*x18; POST_CHECK_MUL(x65);
        x66=rt2 - x64; POST_CHECK_ADD(x66);
        x67=x17*(-x65 - x66*x28); POST_CHECK_MUL(x67);
        x68=-x11*x62 + x55; POST_CHECK_ADD(x68);
        result[2] = x33*(-x39*x65 + x61 + x66*x38) + x42*(2 + x15*x68 - x37*x68 - x45*x65 + x66*x44) + x67;

    }
    else if (x54)
    {
        DEBUG_PRINT("Case (x54) is True.\n");
        x4=mu*mu; POST_CHECK_POW(x4);
        x6=x4*rn*rn + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5; POST_CHECK_ADD(x6);
        x14=2*x13; POST_CHECK_MUL(x14);
        /*@ assert (-x14 + x6) >= 0.;*/
        x15=sqrt(-x14 + x6); POST_CHECK_POW(x15);
        /*@ assert (x14 + x6) >= 0.;*/
        x16=sqrt(x14 + x6); POST_CHECK_POW(x16);
        x17=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x15 - 0.25*x16; POST_CHECK_ADD(x17);
        /*@ assert (x16) != 0.;*/
        x18=1.0*1.0/x16; POST_CHECK_MUL(x18);
        /*@ assert (x13) != 0.;*/
        x24=1.0/x13; POST_CHECK_POW(x24);
        /*@ assert (x15) != 0.;*/
        x27=1.0/x15; POST_CHECK_POW(x27);
        x28=1.0*x27; POST_CHECK_MUL(x28);
        x31=0.5*rt1; POST_CHECK_MUL(x31);
        x32=0.5*mu*ut1; POST_CHECK_MUL(x32);
        x40=0.5*rt2; POST_CHECK_MUL(x40);
        x41=0.5*mu*ut2; POST_CHECK_MUL(x41);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x47=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x47);
        x48=0.25*random1*x47*x15 - 0.25*random1*x47*x16 + x31 + x32; POST_CHECK_ADD(x48);
        x49=random1*x47*x27; POST_CHECK_MUL(x49);
        x50=random1*x47; POST_CHECK_MUL(x50);
        x51=0.25*random2*x47*x15 - 0.25*random2*x47*x16 + x40 + x41; POST_CHECK_ADD(x51);
        x52=random2*x47*x27; POST_CHECK_MUL(x52);
        x53=random2*x47; POST_CHECK_MUL(x53);
        x55=x24*x7; POST_CHECK_MUL(x55);
        x64=x10*x55; POST_CHECK_MUL(x64);
        x65=(rt2 + x64)*x18; POST_CHECK_MUL(x65);
        x66=rt2 - x64; POST_CHECK_ADD(x66);
        x67=x17*(-x65 - x66*x28); POST_CHECK_MUL(x67);

        /* Assignment result[0, 2]=x48*(-x50*x65 + x66*x49) + x51*(2 - x53*x65 + x66*x52) + x67 */
        x4=mu*mu; POST_CHECK_POW(x4);
        x6=x4*rn*rn + rt1*rt1 + rt2*rt2 + x1*x4 + x2*x4 + x5*x5; POST_CHECK_ADD(x6);
        x14=2*x13; POST_CHECK_MUL(x14);
        /*@ assert (-x14 + x6) >= 0.;*/
        x15=sqrt(-x14 + x6); POST_CHECK_POW(x15);
        /*@ assert (x14 + x6) >= 0.;*/
        x16=sqrt(x14 + x6); POST_CHECK_POW(x16);
        x17=0.5*mu*rn + 0.5*mu*x3 + 0.5*un - 0.25*x15 - 0.25*x16; POST_CHECK_ADD(x17);
        /*@ assert (x16) != 0.;*/
        x18=1.0*1.0/x16; POST_CHECK_MUL(x18);
        /*@ assert (x13) != 0.;*/
        x24=1.0/x13; POST_CHECK_POW(x24);
        /*@ assert (x15) != 0.;*/
        x27=1.0/x15; POST_CHECK_POW(x27);
        x28=1.0*x27; POST_CHECK_MUL(x28);
        x31=0.5*rt1; POST_CHECK_MUL(x31);
        x32=0.5*mu*ut1; POST_CHECK_MUL(x32);
        x40=0.5*rt2; POST_CHECK_MUL(x40);
        x41=0.5*mu*ut2; POST_CHECK_MUL(x41);
        /*@ assert (random1*random1 + random2*random2) >= 0.;*/
        /*@ assert (random1*random1 + random2*random2) != 0.;*/
        x47=pow(random1*random1 + random2*random2, -1.0/2.0); POST_CHECK_POW(x47);
        x48=0.25*random1*x47*x15 - 0.25*random1*x47*x16 + x31 + x32; POST_CHECK_ADD(x48);
        x49=random1*x47*x27; POST_CHECK_MUL(x49);
        x50=random1*x47; POST_CHECK_MUL(x50);
        x51=0.25*random2*x47*x15 - 0.25*random2*x47*x16 + x40 + x41; POST_CHECK_ADD(x51);
        x52=random2*x47*x27; POST_CHECK_MUL(x52);
        x53=random2*x47; POST_CHECK_MUL(x53);
        x55=x24*x7; POST_CHECK_MUL(x55);
        x64=x10*x55; POST_CHECK_MUL(x64);
        x65=(rt2 + x64)*x18; POST_CHECK_MUL(x65);
        x66=rt2 - x64; POST_CHECK_ADD(x66);
        x67=x17*(-x65 - x66*x28); POST_CHECK_MUL(x67);
        result[2] = x48*(-x50*x65 + x66*x49) + x51*(2 - x53*x65 + x66*x52) + x67;

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
