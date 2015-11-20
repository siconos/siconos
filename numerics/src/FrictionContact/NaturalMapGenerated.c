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
#define ZERO 0
#define FIX(x) do { if (!isfinite(x)) { DEBUG_PRINTF("%s is not finite\n", #x); x=0; }} while(0)
#define NOT_ZERO(x) fabs((double) x) > ZERO
#define IS_NOT_ZERO(x) fabs((double) x) > ZERO
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
#define pow(x, y) pow((double)x, (double)y)
#endif

#pragma GCC diagnostic ignored "-Wconversion"

// hack, should be prevented in sage/sympy/maple or in code generation
#define sqrt(x) ((x < 0) && ( x > - sqrt(DBL_EPSILON)) ? 0 : (assert(x>=0),sqrt(x)))

// ./nm2.py --ccode --ccodefac --ccodeAB --wrapper
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
    x3=sqrt(ut1*ut1 + ut2*ut2); FIX(x3);
    /*@ assert (x3) >= 0.;*/
    x31=x3 <= 0; FIX(x31);
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
    x36=x3 > 0; FIX(x36);
    int x72;
    double x14;
    double x18;
    double x45;
    double x60;
    double x70;
    double x71;
    x5=mu*ut1; FIX(x5);
    x6=-rt1 + x5; FIX(x6);
    x7=x6*x6; FIX(x7);
    x8=mu*ut2; FIX(x8);
    x9=-rt2 + x8; FIX(x9);
    x10=x9*x9; FIX(x10);
    x11=x10 + x7; FIX(x11);
    /*@ assert (x11) >= 0.;*/
    x12=sqrt(x11); FIX(x12);
    x72=x12 > 0; FIX(x72);
    int x73;
    x73=x12 <= 0; FIX(x73);
    int x80;
    double x61;
    double x62;
    double x64;
    double x79;
    x80=x36 && x72; FIX(x80);
    int x81;
    x81=x36 && x73; FIX(x81);
    if (x31)
    {
        x1=mu*rn; FIX(x1);
        x20=-1.0*un; FIX(x20);
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x26=Heaviside(x20 + x25); FIX(x26);
        x27=0.5*x26; FIX(x27);
        x28=-1.0*x24; FIX(x28);
        x29=Heaviside(x1 + x20 + x28); FIX(x29);
        x30=0.5*x29; FIX(x30);
    }
    else if (x36)
    {
        x1=mu*rn; FIX(x1);
        x2=-un; FIX(x2);
        /*@ assert (-un) >= 0.;*/
        x4=-mu*x3 + x1 + x2; FIX(x4);
        x13=x12 + x4; FIX(x13);
        x17=-x12 + x4; FIX(x17);
        x32=Heaviside(x13); FIX(x32);
        x33=0.5*x32; FIX(x33);
        x34=Heaviside(x17); FIX(x34);
        x35=0.5*x34; FIX(x35);
    }
    else if (x72)
    {
        x1=mu*rn; FIX(x1);
        x2=-un; FIX(x2);
        /*@ assert (-un) >= 0.;*/
        x4=-mu*x3 + x1 + x2; FIX(x4);
        x13=x12 + x4; FIX(x13);
        x14=Max(0, x13); FIX(x14);
        x17=-x12 + x4; FIX(x17);
        x18=Max(0, x17); FIX(x18);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x60=-mu*ut1 + rt1; FIX(x60);
        x70=0.5*x14*x45; FIX(x70);
        x71=0.5*x18*x45; FIX(x71);
    }
    else if (x80)
    {
        x1=mu*rn; FIX(x1);
        x2=-un; FIX(x2);
        /*@ assert (-un) >= 0.;*/
        x4=-mu*x3 + x1 + x2; FIX(x4);
        x13=x12 + x4; FIX(x13);
        x17=-x12 + x4; FIX(x17);
        x32=Heaviside(x13); FIX(x32);
        x34=Heaviside(x17); FIX(x34);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x60=-mu*ut1 + rt1; FIX(x60);
        x61=0.5*x32*x45; FIX(x61);
        x62=x60*x61; FIX(x62);
        x64=0.5*x34*x45; FIX(x64);
        x79=x6*x64; FIX(x79);
    }
    /* Assignment result[0, 0]=x1 + x16 - x19 */
    double x15;
    double x16;
    double x19;x1=mu*rn; FIX(x1);
    x2=-un; FIX(x2);
    /*@ assert (-un) >= 0.;*/
    x4=-mu*x3 + x1 + x2; FIX(x4);
    x13=x12 + x4; FIX(x13);
    x14=Max(0, x13); FIX(x14);
    x15=0.5*x14; FIX(x15);
    x16=-x15; FIX(x16);
    x17=-x12 + x4; FIX(x17);
    x18=Max(0, x17); FIX(x18);
    x19=0.5*x18; FIX(x19);
    result[0] = x1 + x16 - x19;


    /* Assignment result[1, 0]=Piecewise((rt1 - x6*x71 - x60*x70, x72), (rt1, x73)) */

    if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x60=-mu*ut1 + rt1; FIX(x60);
        x70=0.5*x14*x45; FIX(x70);
        x71=0.5*x18*x45; FIX(x71);

        /* Assignment result[1, 0]=rt1 - x6*x71 - x60*x70 */
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x60=-mu*ut1 + rt1; FIX(x60);
        x70=0.5*x14*x45; FIX(x70);
        x71=0.5*x18*x45; FIX(x71);
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
        x45=1.0/x12; FIX(x45);
        x66=-mu*ut2 + rt2; FIX(x66);
        x70=0.5*x14*x45; FIX(x70);
        x71=0.5*x18*x45; FIX(x71);

        /* Assignment result[2, 0]=rt2 - x66*x70 - x71*x9 */
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x66=-mu*ut2 + rt2; FIX(x66);
        x70=0.5*x14*x45; FIX(x70);
        x71=0.5*x18*x45; FIX(x71);
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
        x20=-1.0*un; FIX(x20);
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x26=Heaviside(x20 + x25); FIX(x26);
        x27=0.5*x26; FIX(x27);
        x28=-1.0*x24; FIX(x28);
        x29=Heaviside(x1 + x20 + x28); FIX(x29);
        x30=0.5*x29; FIX(x30);

        /* Assignment result[0, 1]=x27 + x30 */
        x20=-1.0*un; FIX(x20);
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x26=Heaviside(x20 + x25); FIX(x26);
        x27=0.5*x26; FIX(x27);
        x28=-1.0*x24; FIX(x28);
        x29=Heaviside(x1 + x20 + x28); FIX(x29);
        x30=0.5*x29; FIX(x30);
        result[3] = x27 + x30;

    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x32=Heaviside(x13); FIX(x32);
        x33=0.5*x32; FIX(x33);
        x34=Heaviside(x17); FIX(x34);
        x35=0.5*x34; FIX(x35);

        /* Assignment result[0, 1]=x33 + x35 */
        x32=Heaviside(x13); FIX(x32);
        x33=0.5*x32; FIX(x33);
        x34=Heaviside(x17); FIX(x34);
        x35=0.5*x34; FIX(x35);
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
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x28=-1.0*x24; FIX(x28);
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24; FIX(x37);
        x74=Heaviside(x24); FIX(x74);
        x75=rt1*x74; FIX(x75);
        x76=un - 1.0*x1; FIX(x76);
        x77=Heaviside(x28 + x76); FIX(x77);
        x78=Heaviside(x24 + x76); FIX(x78);

        /* Assignment result[1, 1]=-0.5*x37*x75*(x77 - 1.0*x78) */
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x28=-1.0*x24; FIX(x28);
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24; FIX(x37);
        x74=Heaviside(x24); FIX(x74);
        x75=rt1*x74; FIX(x75);
        x76=un - 1.0*x1; FIX(x76);
        x77=Heaviside(x28 + x76); FIX(x77);
        x78=Heaviside(x24 + x76); FIX(x78);
        result[4] = -0.5*x37*x75*(x77 - 1.0*x78);

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");
        x32=Heaviside(x13); FIX(x32);
        x34=Heaviside(x17); FIX(x34);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x60=-mu*ut1 + rt1; FIX(x60);
        x61=0.5*x32*x45; FIX(x61);
        x62=x60*x61; FIX(x62);
        x64=0.5*x34*x45; FIX(x64);
        x79=x6*x64; FIX(x79);

        /* Assignment result[1, 1]=x62 + x79 */
        x32=Heaviside(x13); FIX(x32);
        x34=Heaviside(x17); FIX(x34);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x60=-mu*ut1 + rt1; FIX(x60);
        x61=0.5*x32*x45; FIX(x61);
        x62=x60*x61; FIX(x62);
        x64=0.5*x34*x45; FIX(x64);
        x79=x6*x64; FIX(x79);
        result[4] = x62 + x79;

    }
    else if (x81)
    {
        DEBUG_PRINT("Case (x81) is True.\n");

        /* Assignment result[1, 1]=0 */

        result[4] = 0;
        /*@ assert (result[4]) >= 0.;*/
    }


    /* Assignment result[2, 1]=Piecewise((x37*(x131*x77 - x131*x78 + x114*x130 - x116*x130 + x24*x27 - x24*x30), x31), (x132 + x67, x80), (x33 - x35, x81)) */
    double x67;
    double x114;
    double x116;
    double x130;
    double x131;
    double x132;
    if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x20=-1.0*un; FIX(x20);
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x26=Heaviside(x20 + x25); FIX(x26);
        x27=0.5*x26; FIX(x27);
        x28=-1.0*x24; FIX(x28);
        x29=Heaviside(x1 + x20 + x28); FIX(x29);
        x30=0.5*x29; FIX(x30);
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24; FIX(x37);
        x74=Heaviside(x24); FIX(x74);
        x76=un - 1.0*x1; FIX(x76);
        x77=Heaviside(x28 + x76); FIX(x77);
        x78=Heaviside(x24 + x76); FIX(x78);
        x114=rt2*x78; FIX(x114);
        x116=rt2*x77; FIX(x116);
        x130=0.5*x74; FIX(x130);
        x131=0.5*x24*x74; FIX(x131);

        /* Assignment result[2, 1]=x37*(x131*x77 - x131*x78 + x114*x130 - x116*x130 + x24*x27 - x24*x30) */
        x20=-1.0*un; FIX(x20);
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x26=Heaviside(x20 + x25); FIX(x26);
        x27=0.5*x26; FIX(x27);
        x28=-1.0*x24; FIX(x28);
        x29=Heaviside(x1 + x20 + x28); FIX(x29);
        x30=0.5*x29; FIX(x30);
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24; FIX(x37);
        x74=Heaviside(x24); FIX(x74);
        x76=un - 1.0*x1; FIX(x76);
        x77=Heaviside(x28 + x76); FIX(x77);
        x78=Heaviside(x24 + x76); FIX(x78);
        x114=rt2*x78; FIX(x114);
        x116=rt2*x77; FIX(x116);
        x130=0.5*x74; FIX(x130);
        x131=0.5*x24*x74; FIX(x131);
        result[5] = x37*(x131*x77 - x131*x78 + x114*x130 - x116*x130 + x24*x27 - x24*x30);

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");
        x32=Heaviside(x13); FIX(x32);
        x34=Heaviside(x17); FIX(x34);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x61=0.5*x32*x45; FIX(x61);
        x64=0.5*x34*x45; FIX(x64);
        x66=-mu*ut2 + rt2; FIX(x66);
        x67=x61*x66; FIX(x67);
        x132=x64*x9; FIX(x132);

        /* Assignment result[2, 1]=x132 + x67 */
        x32=Heaviside(x13); FIX(x32);
        x34=Heaviside(x17); FIX(x34);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x61=0.5*x32*x45; FIX(x61);
        x64=0.5*x34*x45; FIX(x64);
        x66=-mu*ut2 + rt2; FIX(x66);
        x67=x61*x66; FIX(x67);
        x132=x64*x9; FIX(x132);
        result[5] = x132 + x67;

    }
    else if (x81)
    {
        DEBUG_PRINT("Case (x81) is True.\n");
        x32=Heaviside(x13); FIX(x32);
        x33=0.5*x32; FIX(x33);
        x34=Heaviside(x17); FIX(x34);
        x35=0.5*x34; FIX(x35);

        /* Assignment result[2, 1]=x33 - x35 */
        x32=Heaviside(x13); FIX(x32);
        x33=0.5*x32; FIX(x33);
        x34=Heaviside(x17); FIX(x34);
        x35=0.5*x34; FIX(x35);
        result[5] = x33 - x35;

    }


    /* Assignment result[0, 2]=Piecewise((x38*(rt1*x39 - rt1*x40 + x42), x31), (x49 - x51, x36)) */
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
        x20=-1.0*un; FIX(x20);
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x26=Heaviside(x20 + x25); FIX(x26);
        x28=-1.0*x24; FIX(x28);
        x29=Heaviside(x1 + x20 + x28); FIX(x29);
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24; FIX(x37);
        x38=0.25*mu*x37; FIX(x38);
        x39=2.0*x26; FIX(x39);
        x40=2.0*x29; FIX(x40);
        x41=1.4142135623730951454746218587388284504413604736328125*x24; FIX(x41);
        x42=x26*x41 + x29*x41; FIX(x42);

        /* Assignment result[0, 2]=x38*(rt1*x39 - rt1*x40 + x42) */
        x20=-1.0*un; FIX(x20);
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x26=Heaviside(x20 + x25); FIX(x26);
        x28=-1.0*x24; FIX(x28);
        x29=Heaviside(x1 + x20 + x28); FIX(x29);
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24; FIX(x37);
        x38=0.25*mu*x37; FIX(x38);
        x39=2.0*x26; FIX(x39);
        x40=2.0*x29; FIX(x40);
        x41=1.4142135623730951454746218587388284504413604736328125*x24; FIX(x41);
        x42=x26*x41 + x29*x41; FIX(x42);
        result[6] = x38*(rt1*x39 - rt1*x40 + x42);

    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x32=Heaviside(x13); FIX(x32);
        x33=0.5*x32; FIX(x33);
        x34=Heaviside(x17); FIX(x34);
        x35=0.5*x34; FIX(x35);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        x44=-x43*x5; FIX(x44);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x47=x46*x6; FIX(x47);
        x48=x44 + x47; FIX(x48);
        x49=-x33*x48; FIX(x49);
        x50=x44 - x47; FIX(x50);
        x51=x35*x50; FIX(x51);

        /* Assignment result[0, 2]=x49 - x51 */
        x32=Heaviside(x13); FIX(x32);
        x33=0.5*x32; FIX(x33);
        x34=Heaviside(x17); FIX(x34);
        x35=0.5*x34; FIX(x35);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        x44=-x43*x5; FIX(x44);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x47=x46*x6; FIX(x47);
        x48=x44 + x47; FIX(x48);
        x49=-x33*x48; FIX(x49);
        x50=x44 - x47; FIX(x50);
        x51=x35*x50; FIX(x51);
        result[6] = x49 - x51;

    }


    /* Assignment result[1, 2]=Piecewise((-x74*x83*(x102*x77 + x102*x78 + x100*x94 + x101 - x21*x85 - x86*x87 + x86*x90 - x87*x89 - x88*x92 + x88*x97 + x89*x90 - x94*x95 - x99), x31), (-x19*(-mu*x104 + x46) - x15*(-mu*x106 + x105) - x48*x62 - x50*x79, x80), (0, x81)) */
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
    double x106;
    if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x28=-1.0*x24; FIX(x28);
        x74=Heaviside(x24); FIX(x74);
        x76=un - 1.0*x1; FIX(x76);
        x77=Heaviside(x28 + x76); FIX(x77);
        x78=Heaviside(x24 + x76); FIX(x78);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=pow(x23, -5.0/2.0); FIX(x82);
        x83=0.25*mu*x82; FIX(x83);
        /*@ assert (x23) >= 0.;*/
        x84=pow(x23, 3.0/2.0); FIX(x84);
        x85=4.0*x84; FIX(x85);
        x86=pow(rt1, 5); FIX(x86);
        x87=1.4142135623730951454746218587388284504413604736328125*x78; FIX(x87);
        x88=pow(rt2, 4); FIX(x88);
        x89=rt1*x88; FIX(x89);
        x90=1.4142135623730951454746218587388284504413604736328125*x77; FIX(x90);
        x91=Max(0, x2 + x25); FIX(x91);
        x92=2.0*x91; FIX(x92);
        x93=pow(rt1, 3); FIX(x93);
        x94=x22*x93; FIX(x94);
        x95=2.828427124746190290949243717477656900882720947265625*x78; FIX(x95);
        x96=Max(0, x1 + x2 - x24); FIX(x96);
        x97=2.0*x96; FIX(x97);
        x98=x21*x22; FIX(x98);
        x99=x92*x98; FIX(x99);
        x100=2.828427124746190290949243717477656900882720947265625*x77; FIX(x100);
        x101=x97*x98; FIX(x101);
        x102=2.0*x21*x84; FIX(x102);

        /* Assignment result[1, 2]=-x74*x83*(x102*x77 + x102*x78 + x100*x94 + x101 - x21*x85 - x86*x87 + x86*x90 - x87*x89 - x88*x92 + x88*x97 + x89*x90 - x94*x95 - x99) */
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x28=-1.0*x24; FIX(x28);
        x74=Heaviside(x24); FIX(x74);
        x76=un - 1.0*x1; FIX(x76);
        x77=Heaviside(x28 + x76); FIX(x77);
        x78=Heaviside(x24 + x76); FIX(x78);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=pow(x23, -5.0/2.0); FIX(x82);
        x83=0.25*mu*x82; FIX(x83);
        /*@ assert (x23) >= 0.;*/
        x84=pow(x23, 3.0/2.0); FIX(x84);
        x85=4.0*x84; FIX(x85);
        x86=pow(rt1, 5); FIX(x86);
        x87=1.4142135623730951454746218587388284504413604736328125*x78; FIX(x87);
        x88=pow(rt2, 4); FIX(x88);
        x89=rt1*x88; FIX(x89);
        x90=1.4142135623730951454746218587388284504413604736328125*x77; FIX(x90);
        x91=Max(0, x2 + x25); FIX(x91);
        x92=2.0*x91; FIX(x92);
        x93=pow(rt1, 3); FIX(x93);
        x94=x22*x93; FIX(x94);
        x95=2.828427124746190290949243717477656900882720947265625*x78; FIX(x95);
        x96=Max(0, x1 + x2 - x24); FIX(x96);
        x97=2.0*x96; FIX(x97);
        x98=x21*x22; FIX(x98);
        x99=x92*x98; FIX(x99);
        x100=2.828427124746190290949243717477656900882720947265625*x77; FIX(x100);
        x101=x97*x98; FIX(x101);
        x102=2.0*x21*x84; FIX(x102);
        result[7] = -x74*x83*(x102*x77 + x102*x78 + x100*x94 + x101 - x21*x85 - x86*x87 + x86*x90 - x87*x89 - x88*x92 + x88*x97 + x89*x90 - x94*x95 - x99);

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");
        x32=Heaviside(x13); FIX(x32);
        x34=Heaviside(x17); FIX(x34);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        x44=-x43*x5; FIX(x44);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x47=x46*x6; FIX(x47);
        x48=x44 + x47; FIX(x48);
        x50=x44 - x47; FIX(x50);
        x60=-mu*ut1 + rt1; FIX(x60);
        x61=0.5*x32*x45; FIX(x61);
        x62=x60*x61; FIX(x62);
        x64=0.5*x34*x45; FIX(x64);
        x79=x6*x64; FIX(x79);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        x104=x103*x7; FIX(x104);
        x105=-x46; FIX(x105);
        x106=x103*x6*x60; FIX(x106);

        /* Assignment result[1, 2]=-x19*(-mu*x104 + x46) - x15*(-mu*x106 + x105) - x48*x62 - x50*x79 */
        x32=Heaviside(x13); FIX(x32);
        x34=Heaviside(x17); FIX(x34);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        x44=-x43*x5; FIX(x44);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x47=x46*x6; FIX(x47);
        x48=x44 + x47; FIX(x48);
        x50=x44 - x47; FIX(x50);
        x60=-mu*ut1 + rt1; FIX(x60);
        x61=0.5*x32*x45; FIX(x61);
        x62=x60*x61; FIX(x62);
        x64=0.5*x34*x45; FIX(x64);
        x79=x6*x64; FIX(x79);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        x104=x103*x7; FIX(x104);
        x105=-x46; FIX(x105);
        x106=x103*x6*x60; FIX(x106);
        result[7] = -x19*(-mu*x104 + x46) - x15*(-mu*x106 + x105) - x48*x62 - x50*x79;

    }
    else if (x81)
    {
        DEBUG_PRINT("Case (x81) is True.\n");

        /* Assignment result[1, 2]=0 */

        result[7] = 0;
        /*@ assert (result[7]) >= 0.;*/
    }


    /* Assignment result[2, 2]=Piecewise((-x83*(-x150*x77 - x150*x78 + x100*x147 - x107*x75 + x110*x75 - x112*x75 + x115*x75 + x117*x75 - x136*x144 - x136*x145 + x137*x144 + x137*x145 + x138 + x139*x86 + x139*x89 + x140*x94 - x141*x94 - x142*x94 - x143*x86 - x143*x89 - x146*x86 - x146*x89 - x147*x95 + x148*x149 - x148*x151 + x152*x87 - x152*x90 - x39*x86 - x39*x89 - x40*x86 - x40*x89), x31), (mu*x153 + x122 - x132*x50 - x48*x67, x80), (x49 + x51, x81)) */
    double x107;
    double x108;
    double x109;
    double x110;
    double x112;
    double x113;
    double x115;
    double x117;
    double x119;
    double x121;
    double x122;
    double x129;
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
    if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x20=-1.0*un; FIX(x20);
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x26=Heaviside(x20 + x25); FIX(x26);
        x28=-1.0*x24; FIX(x28);
        x29=Heaviside(x1 + x20 + x28); FIX(x29);
        x39=2.0*x26; FIX(x39);
        x40=2.0*x29; FIX(x40);
        x74=Heaviside(x24); FIX(x74);
        x75=rt1*x74; FIX(x75);
        x76=un - 1.0*x1; FIX(x76);
        x77=Heaviside(x28 + x76); FIX(x77);
        x78=Heaviside(x24 + x76); FIX(x78);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=pow(x23, -5.0/2.0); FIX(x82);
        x83=0.25*mu*x82; FIX(x83);
        /*@ assert (x23) >= 0.;*/
        x84=pow(x23, 3.0/2.0); FIX(x84);
        x85=4.0*x84; FIX(x85);
        x86=pow(rt1, 5); FIX(x86);
        x87=1.4142135623730951454746218587388284504413604736328125*x78; FIX(x87);
        x88=pow(rt2, 4); FIX(x88);
        x89=rt1*x88; FIX(x89);
        x90=1.4142135623730951454746218587388284504413604736328125*x77; FIX(x90);
        x91=Max(0, x2 + x25); FIX(x91);
        x92=2.0*x91; FIX(x92);
        x93=pow(rt1, 3); FIX(x93);
        x94=x22*x93; FIX(x94);
        x95=2.828427124746190290949243717477656900882720947265625*x78; FIX(x95);
        x96=Max(0, x1 + x2 - x24); FIX(x96);
        x97=2.0*x96; FIX(x97);
        x100=2.828427124746190290949243717477656900882720947265625*x77; FIX(x100);
        x107=rt2*x85; FIX(x107);
        x108=pow(rt1, 4); FIX(x108);
        x109=pow(rt2, 3); FIX(x109);
        x110=x109*x92; FIX(x110);
        x112=x109*x97; FIX(x112);
        x113=2.0*x84; FIX(x113);
        x114=rt2*x78; FIX(x114);
        x115=x113*x114; FIX(x115);
        x116=rt2*x77; FIX(x116);
        x117=x113*x116; FIX(x117);
        x133=1.4142135623730951454746218587388284504413604736328125*x26*x84; FIX(x133);
        x134=1.4142135623730951454746218587388284504413604736328125*x29*x84; FIX(x134);
        x135=x21*x84; FIX(x135);
        x136=1.4142135623730951454746218587388284504413604736328125*x74*x78; FIX(x136);
        x137=1.4142135623730951454746218587388284504413604736328125*x74*x77; FIX(x137);
        x138=-x133*x21 - x133*x22 + x134*x21 + x134*x22 + x135*x136 - x135*x137; FIX(x138);
        x139=4.0*x74; FIX(x139);
        x140=8.0*x74; FIX(x140);
        x141=4.0*x26; FIX(x141);
        x142=4.0*x29; FIX(x142);
        x143=2.0*x74*x78; FIX(x143);
        x144=pow(rt2, 5); FIX(x144);
        x145=rt2*x108; FIX(x145);
        x146=2.0*x74*x77; FIX(x146);
        x147=x109*x21*x74; FIX(x147);
        x148=rt2*x93; FIX(x148);
        x149=2.0*x74*x91; FIX(x149);
        x150=4.0*x22*x74*x93; FIX(x150);
        x151=2.0*x74*x96; FIX(x151);
        x152=x22*x74*x84; FIX(x152);

        /* Assignment result[2, 2]=-x83*(-x150*x77 - x150*x78 + x100*x147 - x107*x75 + x110*x75 - x112*x75 + x115*x75 + x117*x75 - x136*x144 - x136*x145 + x137*x144 + x137*x145 + x138 + x139*x86 + x139*x89 + x140*x94 - x141*x94 - x142*x94 - x143*x86 - x143*x89 - x146*x86 - x146*x89 - x147*x95 + x148*x149 - x148*x151 + x152*x87 - x152*x90 - x39*x86 - x39*x89 - x40*x86 - x40*x89) */
        x20=-1.0*un; FIX(x20);
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x26=Heaviside(x20 + x25); FIX(x26);
        x28=-1.0*x24; FIX(x28);
        x29=Heaviside(x1 + x20 + x28); FIX(x29);
        x39=2.0*x26; FIX(x39);
        x40=2.0*x29; FIX(x40);
        x74=Heaviside(x24); FIX(x74);
        x75=rt1*x74; FIX(x75);
        x76=un - 1.0*x1; FIX(x76);
        x77=Heaviside(x28 + x76); FIX(x77);
        x78=Heaviside(x24 + x76); FIX(x78);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=pow(x23, -5.0/2.0); FIX(x82);
        x83=0.25*mu*x82; FIX(x83);
        /*@ assert (x23) >= 0.;*/
        x84=pow(x23, 3.0/2.0); FIX(x84);
        x85=4.0*x84; FIX(x85);
        x86=pow(rt1, 5); FIX(x86);
        x87=1.4142135623730951454746218587388284504413604736328125*x78; FIX(x87);
        x88=pow(rt2, 4); FIX(x88);
        x89=rt1*x88; FIX(x89);
        x90=1.4142135623730951454746218587388284504413604736328125*x77; FIX(x90);
        x91=Max(0, x2 + x25); FIX(x91);
        x92=2.0*x91; FIX(x92);
        x93=pow(rt1, 3); FIX(x93);
        x94=x22*x93; FIX(x94);
        x95=2.828427124746190290949243717477656900882720947265625*x78; FIX(x95);
        x96=Max(0, x1 + x2 - x24); FIX(x96);
        x97=2.0*x96; FIX(x97);
        x100=2.828427124746190290949243717477656900882720947265625*x77; FIX(x100);
        x107=rt2*x85; FIX(x107);
        x108=pow(rt1, 4); FIX(x108);
        x109=pow(rt2, 3); FIX(x109);
        x110=x109*x92; FIX(x110);
        x112=x109*x97; FIX(x112);
        x113=2.0*x84; FIX(x113);
        x114=rt2*x78; FIX(x114);
        x115=x113*x114; FIX(x115);
        x116=rt2*x77; FIX(x116);
        x117=x113*x116; FIX(x117);
        x133=1.4142135623730951454746218587388284504413604736328125*x26*x84; FIX(x133);
        x134=1.4142135623730951454746218587388284504413604736328125*x29*x84; FIX(x134);
        x135=x21*x84; FIX(x135);
        x136=1.4142135623730951454746218587388284504413604736328125*x74*x78; FIX(x136);
        x137=1.4142135623730951454746218587388284504413604736328125*x74*x77; FIX(x137);
        x138=-x133*x21 - x133*x22 + x134*x21 + x134*x22 + x135*x136 - x135*x137; FIX(x138);
        x139=4.0*x74; FIX(x139);
        x140=8.0*x74; FIX(x140);
        x141=4.0*x26; FIX(x141);
        x142=4.0*x29; FIX(x142);
        x143=2.0*x74*x78; FIX(x143);
        x144=pow(rt2, 5); FIX(x144);
        x145=rt2*x108; FIX(x145);
        x146=2.0*x74*x77; FIX(x146);
        x147=x109*x21*x74; FIX(x147);
        x148=rt2*x93; FIX(x148);
        x149=2.0*x74*x91; FIX(x149);
        x150=4.0*x22*x74*x93; FIX(x150);
        x151=2.0*x74*x96; FIX(x151);
        x152=x22*x74*x84; FIX(x152);
        result[8] = -x83*(-x150*x77 - x150*x78 + x100*x147 - x107*x75 + x110*x75 - x112*x75 + x115*x75 + x117*x75 - x136*x144 - x136*x145 + x137*x144 + x137*x145 + x138 + x139*x86 + x139*x89 + x140*x94 - x141*x94 - x142*x94 - x143*x86 - x143*x89 - x146*x86 - x146*x89 - x147*x95 + x148*x149 - x148*x151 + x152*x87 - x152*x90 - x39*x86 - x39*x89 - x40*x86 - x40*x89);

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");
        x32=Heaviside(x13); FIX(x32);
        x34=Heaviside(x17); FIX(x34);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        x44=-x43*x5; FIX(x44);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x47=x46*x6; FIX(x47);
        x48=x44 + x47; FIX(x48);
        x50=x44 - x47; FIX(x50);
        x61=0.5*x32*x45; FIX(x61);
        x64=0.5*x34*x45; FIX(x64);
        x66=-mu*ut2 + rt2; FIX(x66);
        x67=x61*x66; FIX(x67);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        x119=0.5*x103*x14; FIX(x119);
        x121=x103*x19*x6*x9; FIX(x121);
        x122=mu*x121; FIX(x122);
        x129=x6*x66; FIX(x129);
        x132=x64*x9; FIX(x132);
        x153=x119*x129; FIX(x153);

        /* Assignment result[2, 2]=mu*x153 + x122 - x132*x50 - x48*x67 */
        x32=Heaviside(x13); FIX(x32);
        x34=Heaviside(x17); FIX(x34);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        x44=-x43*x5; FIX(x44);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x47=x46*x6; FIX(x47);
        x48=x44 + x47; FIX(x48);
        x50=x44 - x47; FIX(x50);
        x61=0.5*x32*x45; FIX(x61);
        x64=0.5*x34*x45; FIX(x64);
        x66=-mu*ut2 + rt2; FIX(x66);
        x67=x61*x66; FIX(x67);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        x119=0.5*x103*x14; FIX(x119);
        x121=x103*x19*x6*x9; FIX(x121);
        x122=mu*x121; FIX(x122);
        x129=x6*x66; FIX(x129);
        x132=x64*x9; FIX(x132);
        x153=x119*x129; FIX(x153);
        result[8] = mu*x153 + x122 - x132*x50 - x48*x67;

    }
    else if (x81)
    {
        DEBUG_PRINT("Case (x81) is True.\n");
        x32=Heaviside(x13); FIX(x32);
        x33=0.5*x32; FIX(x33);
        x34=Heaviside(x17); FIX(x34);
        x35=0.5*x34; FIX(x35);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        x44=-x43*x5; FIX(x44);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x47=x46*x6; FIX(x47);
        x48=x44 + x47; FIX(x48);
        x49=-x33*x48; FIX(x49);
        x50=x44 - x47; FIX(x50);
        x51=x35*x50; FIX(x51);

        /* Assignment result[2, 2]=x49 + x51 */
        x32=Heaviside(x13); FIX(x32);
        x33=0.5*x32; FIX(x33);
        x34=Heaviside(x17); FIX(x34);
        x35=0.5*x34; FIX(x35);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        x44=-x43*x5; FIX(x44);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x47=x46*x6; FIX(x47);
        x48=x44 + x47; FIX(x48);
        x49=-x33*x48; FIX(x49);
        x50=x44 - x47; FIX(x50);
        x51=x35*x50; FIX(x51);
        result[8] = x49 + x51;

    }


    /* Assignment result[0, 3]=Piecewise((x38*(rt2*x39 - rt2*x40 + x42), x31), (x55 - x57, x36)) */
    double x52;
    double x53;
    double x54;
    double x55;
    double x56;
    double x57;
    if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x20=-1.0*un; FIX(x20);
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x26=Heaviside(x20 + x25); FIX(x26);
        x28=-1.0*x24; FIX(x28);
        x29=Heaviside(x1 + x20 + x28); FIX(x29);
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24; FIX(x37);
        x38=0.25*mu*x37; FIX(x38);
        x39=2.0*x26; FIX(x39);
        x40=2.0*x29; FIX(x40);
        x41=1.4142135623730951454746218587388284504413604736328125*x24; FIX(x41);
        x42=x26*x41 + x29*x41; FIX(x42);

        /* Assignment result[0, 3]=x38*(rt2*x39 - rt2*x40 + x42) */
        x20=-1.0*un; FIX(x20);
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x26=Heaviside(x20 + x25); FIX(x26);
        x28=-1.0*x24; FIX(x28);
        x29=Heaviside(x1 + x20 + x28); FIX(x29);
        /*@ assert (x24) != 0.;*/
        x37=1.0/x24; FIX(x37);
        x38=0.25*mu*x37; FIX(x38);
        x39=2.0*x26; FIX(x39);
        x40=2.0*x29; FIX(x40);
        x41=1.4142135623730951454746218587388284504413604736328125*x24; FIX(x41);
        x42=x26*x41 + x29*x41; FIX(x42);
        result[9] = x38*(rt2*x39 - rt2*x40 + x42);

    }
    else if (x36)
    {
        DEBUG_PRINT("Case (x36) is True.\n");
        x32=Heaviside(x13); FIX(x32);
        x33=0.5*x32; FIX(x33);
        x34=Heaviside(x17); FIX(x34);
        x35=0.5*x34; FIX(x35);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x52=-x43*x8; FIX(x52);
        x53=x46*x9; FIX(x53);
        x54=x52 + x53; FIX(x54);
        x55=-x33*x54; FIX(x55);
        x56=x52 - x53; FIX(x56);
        x57=x35*x56; FIX(x57);

        /* Assignment result[0, 3]=x55 - x57 */
        x32=Heaviside(x13); FIX(x32);
        x33=0.5*x32; FIX(x33);
        x34=Heaviside(x17); FIX(x34);
        x35=0.5*x34; FIX(x35);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x52=-x43*x8; FIX(x52);
        x53=x46*x9; FIX(x53);
        x54=x52 + x53; FIX(x54);
        x55=-x33*x54; FIX(x55);
        x56=x52 - x53; FIX(x56);
        x57=x35*x56; FIX(x57);
        result[9] = x55 - x57;

    }


    /* Assignment result[1, 3]=Piecewise((0.25*mu*rt1*x74*x82*(-x100*x98 + x107 + x108*x87 - x108*x90 - x110 - x111*x92 + x111*x97 + x112 - x115 - x117 + x87*x88 - x88*x90 + x95*x98), x31), (mu*x120 + x122 - x54*x62 - x56*x79, x80), (0, x81)) */
    double x111;
    double x118;
    double x120;
    if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x28=-1.0*x24; FIX(x28);
        x74=Heaviside(x24); FIX(x74);
        x76=un - 1.0*x1; FIX(x76);
        x77=Heaviside(x28 + x76); FIX(x77);
        x78=Heaviside(x24 + x76); FIX(x78);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=pow(x23, -5.0/2.0); FIX(x82);
        /*@ assert (x23) >= 0.;*/
        x84=pow(x23, 3.0/2.0); FIX(x84);
        x85=4.0*x84; FIX(x85);
        x87=1.4142135623730951454746218587388284504413604736328125*x78; FIX(x87);
        x88=pow(rt2, 4); FIX(x88);
        x90=1.4142135623730951454746218587388284504413604736328125*x77; FIX(x90);
        x91=Max(0, x2 + x25); FIX(x91);
        x92=2.0*x91; FIX(x92);
        x95=2.828427124746190290949243717477656900882720947265625*x78; FIX(x95);
        x96=Max(0, x1 + x2 - x24); FIX(x96);
        x97=2.0*x96; FIX(x97);
        x98=x21*x22; FIX(x98);
        x100=2.828427124746190290949243717477656900882720947265625*x77; FIX(x100);
        x107=rt2*x85; FIX(x107);
        x108=pow(rt1, 4); FIX(x108);
        x109=pow(rt2, 3); FIX(x109);
        x110=x109*x92; FIX(x110);
        x111=rt2*x21; FIX(x111);
        x112=x109*x97; FIX(x112);
        x113=2.0*x84; FIX(x113);
        x114=rt2*x78; FIX(x114);
        x115=x113*x114; FIX(x115);
        x116=rt2*x77; FIX(x116);
        x117=x113*x116; FIX(x117);

        /* Assignment result[1, 3]=0.25*mu*rt1*x74*x82*(-x100*x98 + x107 + x108*x87 - x108*x90 - x110 - x111*x92 + x111*x97 + x112 - x115 - x117 + x87*x88 - x88*x90 + x95*x98) */
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x28=-1.0*x24; FIX(x28);
        x74=Heaviside(x24); FIX(x74);
        x76=un - 1.0*x1; FIX(x76);
        x77=Heaviside(x28 + x76); FIX(x77);
        x78=Heaviside(x24 + x76); FIX(x78);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=pow(x23, -5.0/2.0); FIX(x82);
        /*@ assert (x23) >= 0.;*/
        x84=pow(x23, 3.0/2.0); FIX(x84);
        x85=4.0*x84; FIX(x85);
        x87=1.4142135623730951454746218587388284504413604736328125*x78; FIX(x87);
        x88=pow(rt2, 4); FIX(x88);
        x90=1.4142135623730951454746218587388284504413604736328125*x77; FIX(x90);
        x91=Max(0, x2 + x25); FIX(x91);
        x92=2.0*x91; FIX(x92);
        x95=2.828427124746190290949243717477656900882720947265625*x78; FIX(x95);
        x96=Max(0, x1 + x2 - x24); FIX(x96);
        x97=2.0*x96; FIX(x97);
        x98=x21*x22; FIX(x98);
        x100=2.828427124746190290949243717477656900882720947265625*x77; FIX(x100);
        x107=rt2*x85; FIX(x107);
        x108=pow(rt1, 4); FIX(x108);
        x109=pow(rt2, 3); FIX(x109);
        x110=x109*x92; FIX(x110);
        x111=rt2*x21; FIX(x111);
        x112=x109*x97; FIX(x112);
        x113=2.0*x84; FIX(x113);
        x114=rt2*x78; FIX(x114);
        x115=x113*x114; FIX(x115);
        x116=rt2*x77; FIX(x116);
        x117=x113*x116; FIX(x117);
        result[10] = 0.25*mu*rt1*x74*x82*(-x100*x98 + x107 + x108*x87 - x108*x90 - x110 - x111*x92 + x111*x97 + x112 - x115 - x117 + x87*x88 - x88*x90 + x95*x98);

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");
        x32=Heaviside(x13); FIX(x32);
        x34=Heaviside(x17); FIX(x34);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x52=-x43*x8; FIX(x52);
        x53=x46*x9; FIX(x53);
        x54=x52 + x53; FIX(x54);
        x56=x52 - x53; FIX(x56);
        x60=-mu*ut1 + rt1; FIX(x60);
        x61=0.5*x32*x45; FIX(x61);
        x62=x60*x61; FIX(x62);
        x64=0.5*x34*x45; FIX(x64);
        x79=x6*x64; FIX(x79);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        x118=x60*x9; FIX(x118);
        x119=0.5*x103*x14; FIX(x119);
        x120=x118*x119; FIX(x120);
        x121=x103*x19*x6*x9; FIX(x121);
        x122=mu*x121; FIX(x122);

        /* Assignment result[1, 3]=mu*x120 + x122 - x54*x62 - x56*x79 */
        x32=Heaviside(x13); FIX(x32);
        x34=Heaviside(x17); FIX(x34);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x52=-x43*x8; FIX(x52);
        x53=x46*x9; FIX(x53);
        x54=x52 + x53; FIX(x54);
        x56=x52 - x53; FIX(x56);
        x60=-mu*ut1 + rt1; FIX(x60);
        x61=0.5*x32*x45; FIX(x61);
        x62=x60*x61; FIX(x62);
        x64=0.5*x34*x45; FIX(x64);
        x79=x6*x64; FIX(x79);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        x118=x60*x9; FIX(x118);
        x119=0.5*x103*x14; FIX(x119);
        x120=x118*x119; FIX(x120);
        x121=x103*x19*x6*x9; FIX(x121);
        x122=mu*x121; FIX(x122);
        result[10] = mu*x120 + x122 - x54*x62 - x56*x79;

    }
    else if (x81)
    {
        DEBUG_PRINT("Case (x81) is True.\n");

        /* Assignment result[1, 3]=0 */

        result[10] = 0;
        /*@ assert (result[10]) >= 0.;*/
    }


    /* Assignment result[2, 3]=Piecewise((-x83*(-1.17157287525381*x109*x21*x74*x77 - 6.82842712474619*x109*x21*x74*x78 + x101*x74 - x22*x74*x85 - x74*x99 + 0.585786437626905*x152*x77 + 3.41421356237309*x152*x78 - x108*x149 + x108*x151 + x138 + x139*x144 + x139*x145 + x140*x154 - x141*x154 - x142*x154 - x144*x155 - x144*x156 - x144*x39 - x144*x40 - x145*x155 - x145*x156 - x145*x39 - x145*x40), x31), (-x19*(-mu*x157 + x46) - x15*(-mu*x158 + x105) - x132*x56 - x54*x67, x80), (x55 + x57, x81)) */
    double x154;
    double x155;
    double x156;
    double x157;
    double x158;
    if (x31)
    {
        DEBUG_PRINT("Case (x31) is True.\n");
        x20=-1.0*un; FIX(x20);
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x26=Heaviside(x20 + x25); FIX(x26);
        x28=-1.0*x24; FIX(x28);
        x29=Heaviside(x1 + x20 + x28); FIX(x29);
        x39=2.0*x26; FIX(x39);
        x40=2.0*x29; FIX(x40);
        x74=Heaviside(x24); FIX(x74);
        x76=un - 1.0*x1; FIX(x76);
        x77=Heaviside(x28 + x76); FIX(x77);
        x78=Heaviside(x24 + x76); FIX(x78);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=pow(x23, -5.0/2.0); FIX(x82);
        x83=0.25*mu*x82; FIX(x83);
        /*@ assert (x23) >= 0.;*/
        x84=pow(x23, 3.0/2.0); FIX(x84);
        x85=4.0*x84; FIX(x85);
        x91=Max(0, x2 + x25); FIX(x91);
        x92=2.0*x91; FIX(x92);
        x96=Max(0, x1 + x2 - x24); FIX(x96);
        x97=2.0*x96; FIX(x97);
        x98=x21*x22; FIX(x98);
        x99=x92*x98; FIX(x99);
        x101=x97*x98; FIX(x101);
        x108=pow(rt1, 4); FIX(x108);
        x109=pow(rt2, 3); FIX(x109);
        x133=1.4142135623730951454746218587388284504413604736328125*x26*x84; FIX(x133);
        x134=1.4142135623730951454746218587388284504413604736328125*x29*x84; FIX(x134);
        x135=x21*x84; FIX(x135);
        x136=1.4142135623730951454746218587388284504413604736328125*x74*x78; FIX(x136);
        x137=1.4142135623730951454746218587388284504413604736328125*x74*x77; FIX(x137);
        x138=-x133*x21 - x133*x22 + x134*x21 + x134*x22 + x135*x136 - x135*x137; FIX(x138);
        x139=4.0*x74; FIX(x139);
        x140=8.0*x74; FIX(x140);
        x141=4.0*x26; FIX(x141);
        x142=4.0*x29; FIX(x142);
        x144=pow(rt2, 5); FIX(x144);
        x145=rt2*x108; FIX(x145);
        x149=2.0*x74*x91; FIX(x149);
        x151=2.0*x74*x96; FIX(x151);
        x152=x22*x74*x84; FIX(x152);
        x154=x109*x21; FIX(x154);
        x155=3.41421356237309492343001693370752036571502685546875*x74*x78; FIX(x155);
        x156=0.58578643762690496554768060377682559192180633544921875*x74*x77; FIX(x156);

        /* Assignment result[2, 3]=-x83*(-1.17157287525381*x109*x21*x74*x77 - 6.82842712474619*x109*x21*x74*x78 + x101*x74 - x22*x74*x85 - x74*x99 + 0.585786437626905*x152*x77 + 3.41421356237309*x152*x78 - x108*x149 + x108*x151 + x138 + x139*x144 + x139*x145 + x140*x154 - x141*x154 - x142*x154 - x144*x155 - x144*x156 - x144*x39 - x144*x40 - x145*x155 - x145*x156 - x145*x39 - x145*x40) */
        x20=-1.0*un; FIX(x20);
        x21=rt1*rt1; FIX(x21);
        x22=rt2*rt2; FIX(x22);
        x23=x21 + x22; FIX(x23);
        /*@ assert (x23) >= 0.;*/
        x24=sqrt(x23); FIX(x24);
        x25=x1 + x24; FIX(x25);
        x26=Heaviside(x20 + x25); FIX(x26);
        x28=-1.0*x24; FIX(x28);
        x29=Heaviside(x1 + x20 + x28); FIX(x29);
        x39=2.0*x26; FIX(x39);
        x40=2.0*x29; FIX(x40);
        x74=Heaviside(x24); FIX(x74);
        x76=un - 1.0*x1; FIX(x76);
        x77=Heaviside(x28 + x76); FIX(x77);
        x78=Heaviside(x24 + x76); FIX(x78);
        /*@ assert (x23) >= 0.;*/
        /*@ assert (x23) != 0.;*/
        x82=pow(x23, -5.0/2.0); FIX(x82);
        x83=0.25*mu*x82; FIX(x83);
        /*@ assert (x23) >= 0.;*/
        x84=pow(x23, 3.0/2.0); FIX(x84);
        x85=4.0*x84; FIX(x85);
        x91=Max(0, x2 + x25); FIX(x91);
        x92=2.0*x91; FIX(x92);
        x96=Max(0, x1 + x2 - x24); FIX(x96);
        x97=2.0*x96; FIX(x97);
        x98=x21*x22; FIX(x98);
        x99=x92*x98; FIX(x99);
        x101=x97*x98; FIX(x101);
        x108=pow(rt1, 4); FIX(x108);
        x109=pow(rt2, 3); FIX(x109);
        x133=1.4142135623730951454746218587388284504413604736328125*x26*x84; FIX(x133);
        x134=1.4142135623730951454746218587388284504413604736328125*x29*x84; FIX(x134);
        x135=x21*x84; FIX(x135);
        x136=1.4142135623730951454746218587388284504413604736328125*x74*x78; FIX(x136);
        x137=1.4142135623730951454746218587388284504413604736328125*x74*x77; FIX(x137);
        x138=-x133*x21 - x133*x22 + x134*x21 + x134*x22 + x135*x136 - x135*x137; FIX(x138);
        x139=4.0*x74; FIX(x139);
        x140=8.0*x74; FIX(x140);
        x141=4.0*x26; FIX(x141);
        x142=4.0*x29; FIX(x142);
        x144=pow(rt2, 5); FIX(x144);
        x145=rt2*x108; FIX(x145);
        x149=2.0*x74*x91; FIX(x149);
        x151=2.0*x74*x96; FIX(x151);
        x152=x22*x74*x84; FIX(x152);
        x154=x109*x21; FIX(x154);
        x155=3.41421356237309492343001693370752036571502685546875*x74*x78; FIX(x155);
        x156=0.58578643762690496554768060377682559192180633544921875*x74*x77; FIX(x156);
        result[11] = -x83*(-1.1715728752538099310953612075536511838436126708984375*x109*x21*x74*x77 - 6.8284271247461898468600338674150407314300537109375*x109*x21*x74*x78 + x101*x74 - x22*x74*x85 - x74*x99 + 0.58578643762690496554768060377682559192180633544921875*x152*x77 + 3.41421356237309492343001693370752036571502685546875*x152*x78 - x108*x149 + x108*x151 + x138 + x139*x144 + x139*x145 + x140*x154 - x141*x154 - x142*x154 - x144*x155 - x144*x156 - x144*x39 - x144*x40 - x145*x155 - x145*x156 - x145*x39 - x145*x40);

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");
        x32=Heaviside(x13); FIX(x32);
        x34=Heaviside(x17); FIX(x34);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x52=-x43*x8; FIX(x52);
        x53=x46*x9; FIX(x53);
        x54=x52 + x53; FIX(x54);
        x56=x52 - x53; FIX(x56);
        x61=0.5*x32*x45; FIX(x61);
        x64=0.5*x34*x45; FIX(x64);
        x66=-mu*ut2 + rt2; FIX(x66);
        x67=x61*x66; FIX(x67);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        x105=-x46; FIX(x105);
        x132=x64*x9; FIX(x132);
        x157=x10*x103; FIX(x157);
        x158=x103*x66*x9; FIX(x158);

        /* Assignment result[2, 3]=-x19*(-mu*x157 + x46) - x15*(-mu*x158 + x105) - x132*x56 - x54*x67 */
        x32=Heaviside(x13); FIX(x32);
        x34=Heaviside(x17); FIX(x34);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x52=-x43*x8; FIX(x52);
        x53=x46*x9; FIX(x53);
        x54=x52 + x53; FIX(x54);
        x56=x52 - x53; FIX(x56);
        x61=0.5*x32*x45; FIX(x61);
        x64=0.5*x34*x45; FIX(x64);
        x66=-mu*ut2 + rt2; FIX(x66);
        x67=x61*x66; FIX(x67);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        x105=-x46; FIX(x105);
        x132=x64*x9; FIX(x132);
        x157=x10*x103; FIX(x157);
        x158=x103*x66*x9; FIX(x158);
        result[11] = -x19*(-mu*x157 + x46) - x15*(-mu*x158 + x105) - x132*x56 - x54*x67;

    }
    else if (x81)
    {
        DEBUG_PRINT("Case (x81) is True.\n");
        x32=Heaviside(x13); FIX(x32);
        x33=0.5*x32; FIX(x33);
        x34=Heaviside(x17); FIX(x34);
        x35=0.5*x34; FIX(x35);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x52=-x43*x8; FIX(x52);
        x53=x46*x9; FIX(x53);
        x54=x52 + x53; FIX(x54);
        x55=-x33*x54; FIX(x55);
        x56=x52 - x53; FIX(x56);
        x57=x35*x56; FIX(x57);

        /* Assignment result[2, 3]=x55 + x57 */
        x32=Heaviside(x13); FIX(x32);
        x33=0.5*x32; FIX(x33);
        x34=Heaviside(x17); FIX(x34);
        x35=0.5*x34; FIX(x35);
        /*@ assert (x3) != 0.;*/
        x43=1.0/x3; FIX(x43);
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x52=-x43*x8; FIX(x52);
        x53=x46*x9; FIX(x53);
        x54=x52 + x53; FIX(x54);
        x55=-x33*x54; FIX(x55);
        x56=x52 - x53; FIX(x56);
        x57=x35*x56; FIX(x57);
        result[11] = x55 + x57;

    }


    /* Assignment result[0, 4]=mu + x58 - x59 */
    double x58;
    double x59;x32=Heaviside(x13); FIX(x32);
    x33=0.5*x32; FIX(x33);
    x34=Heaviside(x17); FIX(x34);
    x35=0.5*x34; FIX(x35);
    x58=-mu*x33; FIX(x58);
    x59=mu*x35; FIX(x59);
    result[12] = mu + x58 - x59;


    /* Assignment result[1, 4]=Piecewise((-x123*x60 - x35*x47, x72), (0, x73)) */
    double x123;
    if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x47=x46*x6; FIX(x47);
        x60=-mu*ut1 + rt1; FIX(x60);
        x123=0.5*mu*x32*x45; FIX(x123);

        /* Assignment result[1, 4]=-x123*x60 - x35*x47 */
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x47=x46*x6; FIX(x47);
        x60=-mu*ut1 + rt1; FIX(x60);
        x123=0.5*mu*x32*x45; FIX(x123);
        result[13] = -x123*x60 - x35*x47;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[1, 4]=0 */

        result[13] = 0;
        /*@ assert (result[13]) >= 0.;*/
    }


    /* Assignment result[2, 4]=Piecewise((-x123*x66 - x35*x53, x72), (x58 + x59, x73)) */

    if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x53=x46*x9; FIX(x53);
        x66=-mu*ut2 + rt2; FIX(x66);
        x123=0.5*mu*x32*x45; FIX(x123);

        /* Assignment result[2, 4]=-x123*x66 - x35*x53 */
        /*@ assert (x12) != 0.;*/
        x45=1.0/x12; FIX(x45);
        x46=mu*x45; FIX(x46);
        x53=x46*x9; FIX(x53);
        x66=-mu*ut2 + rt2; FIX(x66);
        x123=0.5*mu*x32*x45; FIX(x123);
        result[14] = -x123*x66 - x35*x53;

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
    x45=1.0/x12; FIX(x45);
    x60=-mu*ut1 + rt1; FIX(x60);
    x61=0.5*x32*x45; FIX(x61);
    x62=x60*x61; FIX(x62);
    x63=-x62; FIX(x63);
    x64=0.5*x34*x45; FIX(x64);
    x65=x60*x64; FIX(x65);
    result[15] = x63 + x65;


    /* Assignment result[1, 5]=Piecewise((-x125*x60**2 - x19*(x104 + x127) - x15*(x106 + x45) + 1 + x126*x6*x60, x72), (1, x73)) */
    double x124;
    double x125;
    double x126;
    double x127;
    if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        x104=x103*x7; FIX(x104);
        x106=x103*x6*x60; FIX(x106);
        /*@ assert (x11) != 0.;*/
        x124=1.0/x11; FIX(x124);
        x125=0.5*x124*x32; FIX(x125);
        x126=0.5*x124*x34; FIX(x126);
        x127=-x45; FIX(x127);

        /* Assignment result[1, 5]=-x125*x60**2 - x19*(x104 + x127) - x15*(x106 + x45) + 1 + x126*x6*x60 */
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        x104=x103*x7; FIX(x104);
        x106=x103*x6*x60; FIX(x106);
        /*@ assert (x11) != 0.;*/
        x124=1.0/x11; FIX(x124);
        x125=0.5*x124*x32; FIX(x125);
        x126=0.5*x124*x34; FIX(x126);
        x127=-x45; FIX(x127);
        result[16] = -x125*x60*x60 - x19*(x104 + x127) - x15*(x106 + x45) + 1 + x126*x6*x60;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[1, 5]=1 */

        result[16] = 1;
        /*@ assert (result[16]) >= 0.;*/
        /*@ assert (result[16]) != 0.;*/
    }


    /* Assignment result[2, 5]=Piecewise((x118*x126 + x128 - x153, x72), (x63 - x65, x73)) */
    double x128;
    if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        x66=-mu*ut2 + rt2; FIX(x66);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        x118=x60*x9; FIX(x118);
        x119=0.5*x103*x14; FIX(x119);
        x121=x103*x19*x6*x9; FIX(x121);
        /*@ assert (x11) != 0.;*/
        x124=1.0/x11; FIX(x124);
        x125=0.5*x124*x32; FIX(x125);
        x126=0.5*x124*x34; FIX(x126);
        x128=-x121 - x125*x60*x66; FIX(x128);
        x129=x6*x66; FIX(x129);
        x153=x119*x129; FIX(x153);

        /* Assignment result[2, 5]=x118*x126 + x128 - x153 */
        x66=-mu*ut2 + rt2; FIX(x66);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        x118=x60*x9; FIX(x118);
        x119=0.5*x103*x14; FIX(x119);
        x121=x103*x19*x6*x9; FIX(x121);
        /*@ assert (x11) != 0.;*/
        x124=1.0/x11; FIX(x124);
        x125=0.5*x124*x32; FIX(x125);
        x126=0.5*x124*x34; FIX(x126);
        x128=-x121 - x125*x60*x66; FIX(x128);
        x129=x6*x66; FIX(x129);
        x153=x119*x129; FIX(x153);
        result[17] = x118*x126 + x128 - x153;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[2, 5]=x63 - x65 */

        result[17] = x63 - x65;

    }


    /* Assignment result[0, 6]=x68 + x69 */
    double x68;
    double x69;x66=-mu*ut2 + rt2; FIX(x66);
    x67=x61*x66; FIX(x67);
    x68=-x67; FIX(x68);
    x69=x64*x66; FIX(x69);
    result[18] = x68 + x69;


    /* Assignment result[1, 6]=Piecewise((-x120 + x126*x129 + x128, x72), (0, x73)) */

    if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        x118=x60*x9; FIX(x118);
        x119=0.5*x103*x14; FIX(x119);
        x120=x118*x119; FIX(x120);
        x121=x103*x19*x6*x9; FIX(x121);
        /*@ assert (x11) != 0.;*/
        x124=1.0/x11; FIX(x124);
        x125=0.5*x124*x32; FIX(x125);
        x126=0.5*x124*x34; FIX(x126);
        x128=-x121 - x125*x60*x66; FIX(x128);
        x129=x6*x66; FIX(x129);

        /* Assignment result[1, 6]=-x120 + x126*x129 + x128 */
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        x118=x60*x9; FIX(x118);
        x119=0.5*x103*x14; FIX(x119);
        x120=x118*x119; FIX(x120);
        x121=x103*x19*x6*x9; FIX(x121);
        /*@ assert (x11) != 0.;*/
        x124=1.0/x11; FIX(x124);
        x125=0.5*x124*x32; FIX(x125);
        x126=0.5*x124*x34; FIX(x126);
        x128=-x121 - x125*x60*x66; FIX(x128);
        x129=x6*x66; FIX(x129);
        result[19] = -x120 + x126*x129 + x128;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[1, 6]=0 */

        result[19] = 0;
        /*@ assert (result[19]) >= 0.;*/
    }


    /* Assignment result[2, 6]=Piecewise((-x125*x66**2 - x19*(x127 + x157) - x15*(x158 + x45) + 1 + x126*x66*x9, x72), (1 + x68 - x69, x73)) */

    if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        /*@ assert (x11) != 0.;*/
        x124=1.0/x11; FIX(x124);
        x125=0.5*x124*x32; FIX(x125);
        x126=0.5*x124*x34; FIX(x126);
        x127=-x45; FIX(x127);
        x157=x10*x103; FIX(x157);
        x158=x103*x66*x9; FIX(x158);

        /* Assignment result[2, 6]=-x125*x66**2 - x19*(x127 + x157) - x15*(x158 + x45) + 1 + x126*x66*x9 */
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x103=pow(x11, -3.0/2.0); FIX(x103);
        /*@ assert (x11) != 0.;*/
        x124=1.0/x11; FIX(x124);
        x125=0.5*x124*x32; FIX(x125);
        x126=0.5*x124*x34; FIX(x126);
        x127=-x45; FIX(x127);
        x157=x10*x103; FIX(x157);
        x158=x103*x66*x9; FIX(x158);
        result[20] = -x125*x66*x66 - x19*(x127 + x157) - x15*(x158 + x45) + 1 + x126*x66*x9;

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
    int x13;
    double x1;
    double x2;
    double x6;
    double x8;
    double x10;
    double x11;
    double x12;
    x3=mu*ut1 - rt1; FIX(x3);
    x4=mu*ut2 - rt2; FIX(x4);
    /*@ assert (x3*x3 + x4*x4) >= 0.;*/
    x5=sqrt(x3*x3 + x4*x4); FIX(x5);
    x13=x5 > 0; FIX(x13);
    int x14;
    x14=x5 <= 0; FIX(x14);
    if (x13)
    {
        x1=mu*rn; FIX(x1);
        /*@ assert (ut1*ut1 + ut2*ut2) >= 0.;*/
        x2=-mu*sqrt(ut1*ut1 + ut2*ut2) - un + x1; FIX(x2);
        x6=Max(0, x2 + x5); FIX(x6);
        x8=Max(0, x2 - x5); FIX(x8);
        /*@ assert (x5) != 0.;*/
        x10=1.0/x5; FIX(x10);
        x11=0.5*x10*x6; FIX(x11);
        x12=0.5*x8*x10; FIX(x12);
    }
    /* Assignment result[0, 0]=x1 + x7 - x9 */
    double x7;
    double x9;x1=mu*rn; FIX(x1);
    /*@ assert (ut1*ut1 + ut2*ut2) >= 0.;*/
    x2=-mu*sqrt(ut1*ut1 + ut2*ut2) - un + x1; FIX(x2);
    x6=Max(0, x2 + x5); FIX(x6);
    x7=-0.5*x6; FIX(x7);
    x8=Max(0, x2 - x5); FIX(x8);
    x9=0.5*x8; FIX(x9);
    result[0] = x1 + x7 - x9;


    /* Assignment result[1, 0]=Piecewise((rt1 - x3*x12 - (-mu*ut1 + rt1)*x11, x13), (rt1, x14)) */

    if (x13)
    {
        DEBUG_PRINT("Case (x13) is True.\n");
        /*@ assert (x5) != 0.;*/
        x10=1.0/x5; FIX(x10);
        x11=0.5*x10*x6; FIX(x11);
        x12=0.5*x8*x10; FIX(x12);

        /* Assignment result[1, 0]=rt1 - x3*x12 - (-mu*ut1 + rt1)*x11 */
        /*@ assert (x5) != 0.;*/
        x10=1.0/x5; FIX(x10);
        x11=0.5*x10*x6; FIX(x11);
        x12=0.5*x8*x10; FIX(x12);
        result[1] = rt1 - x3*x12 - (-mu*ut1 + rt1)*x11;

    }
    else if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");

        /* Assignment result[1, 0]=rt1 */

        result[1] = rt1;

    }
    /*@ assert (result[1]) >= 0.;*/

    /* Assignment result[2, 0]=Piecewise((rt2 - x4*x12 - (-mu*ut2 + rt2)*x11, x13), (rt2 + x7 + x9, x14)) */

    if (x13)
    {
        DEBUG_PRINT("Case (x13) is True.\n");
        /*@ assert (x5) != 0.;*/
        x10=1.0/x5; FIX(x10);
        x11=0.5*x10*x6; FIX(x11);
        x12=0.5*x8*x10; FIX(x12);

        /* Assignment result[2, 0]=rt2 - x4*x12 - (-mu*ut2 + rt2)*x11 */
        /*@ assert (x5) != 0.;*/
        x10=1.0/x5; FIX(x10);
        x11=0.5*x10*x6; FIX(x11);
        x12=0.5*x8*x10; FIX(x12);
        result[2] = rt2 - x4*x12 - (-mu*ut2 + rt2)*x11;

    }
    else if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");

        /* Assignment result[2, 0]=rt2 + x7 + x9 */

        result[2] = rt2 + x7 + x9;

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
    x13=sqrt(ut1*ut1 + ut2*ut2); FIX(x13);
    /*@ assert (x13) >= 0.;*/
    x14=x13 <= 0; FIX(x14);
    int x29;
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
    x29=x13 > 0; FIX(x29);
    int x71;
    int x72;
    double x36;
    double x39;
    double x55;
    double x56;
    double x57;
    double x59;
    double x70;
    x17=mu*ut1 - rt1; FIX(x17);
    x18=x17*x17; FIX(x18);
    x19=mu*ut2 - rt2; FIX(x19);
    x20=x19*x19; FIX(x20);
    x21=x18 + x20; FIX(x21);
    /*@ assert (x21) >= 0.;*/
    x22=sqrt(x21); FIX(x22);
    x71=x22 > 0; FIX(x71);
    x72=x29 && x71; FIX(x72);
    int x73;
    int x74;
    x73=x22 <= 0; FIX(x73);
    x74=x29 && x73; FIX(x74);
    double x40;
    double x41;
    double x120;
    if (x14)
    {
        x1=-1.0*un; FIX(x1);
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x8=Heaviside(x1 + x7); FIX(x8);
        x9=0.5*x8; FIX(x9);
        x10=-1.0*x6; FIX(x10);
        x11=Heaviside(x2 + x1 + x10); FIX(x11);
        x12=0.5*x11; FIX(x12);
    }
    else if (x29)
    {
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x25=0.5*x24; FIX(x25);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x28=0.5*x27; FIX(x28);
    }
    else if (x72)
    {
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x36=mu*ut1; FIX(x36);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x55=rt1 - x36; FIX(x55);
        x56=0.5*x39*x24; FIX(x56);
        x57=x55*x56; FIX(x57);
        x59=0.5*x39*x27; FIX(x59);
        x70=x17*x59; FIX(x70);
    }
    else if (x71)
    {
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x28=0.5*x27; FIX(x28);
        x36=mu*ut1; FIX(x36);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x41=x17*x40; FIX(x41);
        x55=rt1 - x36; FIX(x55);
        x120=0.5*mu*x39*x24; FIX(x120);
    }
    /* Assignment result[0, 0]=Piecewise((x12 + x9, x14), (x25 + x28, x29)) */

    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x1=-1.0*un; FIX(x1);
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x8=Heaviside(x1 + x7); FIX(x8);
        x9=0.5*x8; FIX(x9);
        x10=-1.0*x6; FIX(x10);
        x11=Heaviside(x2 + x1 + x10); FIX(x11);
        x12=0.5*x11; FIX(x12);

        /* Assignment result[0, 0]=x12 + x9 */
        x1=-1.0*un; FIX(x1);
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x8=Heaviside(x1 + x7); FIX(x8);
        x9=0.5*x8; FIX(x9);
        x10=-1.0*x6; FIX(x10);
        x11=Heaviside(x2 + x1 + x10); FIX(x11);
        x12=0.5*x11; FIX(x12);
        result[0] = x12 + x9;

    }
    else if (x29)
    {
        DEBUG_PRINT("Case (x29) is True.\n");
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x25=0.5*x24; FIX(x25);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x28=0.5*x27; FIX(x28);

        /* Assignment result[0, 0]=x25 + x28 */
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x25=0.5*x24; FIX(x25);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x28=0.5*x27; FIX(x28);
        result[0] = x25 + x28;

    }
    /*@ assert (result[0]) >= 0.;*/

    /* Assignment result[1, 0]=Piecewise((-0.5*x66*x30*(x68 - 1.0*x69), x14), (x57 + x70, x72), (0, x74)) */
    double x30;
    double x65;
    double x66;
    double x67;
    double x68;
    double x69;
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x10=-1.0*x6; FIX(x10);
        /*@ assert (x6) != 0.;*/
        x30=1.0/x6; FIX(x30);
        x65=Heaviside(x6); FIX(x65);
        x66=rt1*x65; FIX(x66);
        x67=un - 1.0*x2; FIX(x67);
        x68=Heaviside(x10 + x67); FIX(x68);
        x69=Heaviside(x6 + x67); FIX(x69);

        /* Assignment result[1, 0]=-0.5*x66*x30*(x68 - 1.0*x69) */
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x10=-1.0*x6; FIX(x10);
        /*@ assert (x6) != 0.;*/
        x30=1.0/x6; FIX(x30);
        x65=Heaviside(x6); FIX(x65);
        x66=rt1*x65; FIX(x66);
        x67=un - 1.0*x2; FIX(x67);
        x68=Heaviside(x10 + x67); FIX(x68);
        x69=Heaviside(x6 + x67); FIX(x69);
        result[1] = -0.5*x66*x30*(x68 - 1.0*x69);

    }
    else if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x36=mu*ut1; FIX(x36);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x55=rt1 - x36; FIX(x55);
        x56=0.5*x39*x24; FIX(x56);
        x57=x55*x56; FIX(x57);
        x59=0.5*x39*x27; FIX(x59);
        x70=x17*x59; FIX(x70);

        /* Assignment result[1, 0]=x57 + x70 */
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x36=mu*ut1; FIX(x36);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x55=rt1 - x36; FIX(x55);
        x56=0.5*x39*x24; FIX(x56);
        x57=x55*x56; FIX(x57);
        x59=0.5*x39*x27; FIX(x59);
        x70=x17*x59; FIX(x70);
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

    /* Assignment result[2, 0]=Piecewise((x30*(x68*x128 - x69*x128 + x127*x111 - x127*x113 - x6*x12 + x6*x9), x14), (x129 + x62, x72), (x25 - x28, x74)) */
    double x46;
    double x61;
    double x62;
    double x111;
    double x113;
    double x127;
    double x128;
    double x129;
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x1=-1.0*un; FIX(x1);
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x8=Heaviside(x1 + x7); FIX(x8);
        x9=0.5*x8; FIX(x9);
        x10=-1.0*x6; FIX(x10);
        x11=Heaviside(x2 + x1 + x10); FIX(x11);
        x12=0.5*x11; FIX(x12);
        /*@ assert (x6) != 0.;*/
        x30=1.0/x6; FIX(x30);
        x65=Heaviside(x6); FIX(x65);
        x67=un - 1.0*x2; FIX(x67);
        x68=Heaviside(x10 + x67); FIX(x68);
        x69=Heaviside(x6 + x67); FIX(x69);
        x111=rt2*x69; FIX(x111);
        x113=rt2*x68; FIX(x113);
        x127=0.5*x65; FIX(x127);
        x128=0.5*x6*x65; FIX(x128);

        /* Assignment result[2, 0]=x30*(x68*x128 - x69*x128 + x127*x111 - x127*x113 - x6*x12 + x6*x9) */
        x1=-1.0*un; FIX(x1);
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x8=Heaviside(x1 + x7); FIX(x8);
        x9=0.5*x8; FIX(x9);
        x10=-1.0*x6; FIX(x10);
        x11=Heaviside(x2 + x1 + x10); FIX(x11);
        x12=0.5*x11; FIX(x12);
        /*@ assert (x6) != 0.;*/
        x30=1.0/x6; FIX(x30);
        x65=Heaviside(x6); FIX(x65);
        x67=un - 1.0*x2; FIX(x67);
        x68=Heaviside(x10 + x67); FIX(x68);
        x69=Heaviside(x6 + x67); FIX(x69);
        x111=rt2*x69; FIX(x111);
        x113=rt2*x68; FIX(x113);
        x127=0.5*x65; FIX(x127);
        x128=0.5*x6*x65; FIX(x128);
        result[2] = x30*(x68*x128 - x69*x128 + x127*x111 - x127*x113 - x6*x12 + x6*x9);

    }
    else if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x46=mu*ut2; FIX(x46);
        x56=0.5*x39*x24; FIX(x56);
        x59=0.5*x39*x27; FIX(x59);
        x61=rt2 - x46; FIX(x61);
        x62=x61*x56; FIX(x62);
        x129=x19*x59; FIX(x129);

        /* Assignment result[2, 0]=x129 + x62 */
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x46=mu*ut2; FIX(x46);
        x56=0.5*x39*x24; FIX(x56);
        x59=0.5*x39*x27; FIX(x59);
        x61=rt2 - x46; FIX(x61);
        x62=x61*x56; FIX(x62);
        x129=x19*x59; FIX(x129);
        result[2] = x129 + x62;

    }
    else if (x74)
    {
        DEBUG_PRINT("Case (x74) is True.\n");
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x25=0.5*x24; FIX(x25);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x28=0.5*x27; FIX(x28);

        /* Assignment result[2, 0]=x25 - x28 */
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x25=0.5*x24; FIX(x25);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x28=0.5*x27; FIX(x28);
        result[2] = x25 - x28;

    }
    /*@ assert (result[2]) >= 0.;*/

    /* Assignment result[0, 1]=Piecewise((x31*(rt1*x32 - rt1*x33 + x35), x14), (x43 - x45, x29)) */
    double x31;
    double x32;
    double x33;
    double x34;
    double x35;
    double x37;
    double x38;
    double x42;
    double x43;
    double x44;
    double x45;
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x1=-1.0*un; FIX(x1);
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x8=Heaviside(x1 + x7); FIX(x8);
        x10=-1.0*x6; FIX(x10);
        x11=Heaviside(x2 + x1 + x10); FIX(x11);
        /*@ assert (x6) != 0.;*/
        x30=1.0/x6; FIX(x30);
        x31=0.25*mu*x30; FIX(x31);
        x32=2.0*x8; FIX(x32);
        x33=2.0*x11; FIX(x33);
        x34=1.4142135623730951454746218587388284504413604736328125*x6; FIX(x34);
        x35=x34*x11 + x8*x34; FIX(x35);

        /* Assignment result[0, 1]=x31*(rt1*x32 - rt1*x33 + x35) */
        x1=-1.0*un; FIX(x1);
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x8=Heaviside(x1 + x7); FIX(x8);
        x10=-1.0*x6; FIX(x10);
        x11=Heaviside(x2 + x1 + x10); FIX(x11);
        /*@ assert (x6) != 0.;*/
        x30=1.0/x6; FIX(x30);
        x31=0.25*mu*x30; FIX(x31);
        x32=2.0*x8; FIX(x32);
        x33=2.0*x11; FIX(x33);
        x34=1.4142135623730951454746218587388284504413604736328125*x6; FIX(x34);
        x35=x34*x11 + x8*x34; FIX(x35);
        result[3] = x31*(rt1*x32 - rt1*x33 + x35);

    }
    else if (x29)
    {
        DEBUG_PRINT("Case (x29) is True.\n");
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x25=0.5*x24; FIX(x25);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x28=0.5*x27; FIX(x28);
        x36=mu*ut1; FIX(x36);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        x38=-x36*x37; FIX(x38);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x41=x17*x40; FIX(x41);
        x42=x38 + x41; FIX(x42);
        x43=-x25*x42; FIX(x43);
        x44=x38 - x41; FIX(x44);
        x45=x44*x28; FIX(x45);

        /* Assignment result[0, 1]=x43 - x45 */
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x25=0.5*x24; FIX(x25);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x28=0.5*x27; FIX(x28);
        x36=mu*ut1; FIX(x36);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        x38=-x36*x37; FIX(x38);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x41=x17*x40; FIX(x41);
        x42=x38 + x41; FIX(x42);
        x43=-x25*x42; FIX(x43);
        x44=x38 - x41; FIX(x44);
        x45=x44*x28; FIX(x45);
        result[3] = x43 - x45;

    }
    /*@ assert (result[3]) >= 0.;*/

    /* Assignment result[1, 1]=Piecewise((-x76*x65*(x95*x68 + x95*x69 - x3*x78 + x79*x83 - x80*x79 - x80*x82 - x81*x85 + x81*x90 + x82*x83 - x87*x88 + x87*x93 - x92 + x94), x14), (-x101*(-mu*x103 + x102) - x42*x57 - x44*x70 - x97*(-mu*x99 + x40), x72), (0, x74)) */
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
    double x103;
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x10=-1.0*x6; FIX(x10);
        x15=-un; FIX(x15);
        x65=Heaviside(x6); FIX(x65);
        x67=un - 1.0*x2; FIX(x67);
        x68=Heaviside(x10 + x67); FIX(x68);
        x69=Heaviside(x6 + x67); FIX(x69);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=pow(x5, -5.0/2.0); FIX(x75);
        x76=0.25*mu*x75; FIX(x76);
        /*@ assert (x5) >= 0.;*/
        x77=pow(x5, 3.0/2.0); FIX(x77);
        x78=4.0*x77; FIX(x78);
        x79=pow(rt1, 5); FIX(x79);
        x80=1.4142135623730951454746218587388284504413604736328125*x69; FIX(x80);
        x81=pow(rt2, 4); FIX(x81);
        x82=rt1*x81; FIX(x82);
        x83=1.4142135623730951454746218587388284504413604736328125*x68; FIX(x83);
        x84=Max(0, x15 + x7); FIX(x84);
        x85=2.0*x84; FIX(x85);
        x86=pow(rt1, 3); FIX(x86);
        x87=x4*x86; FIX(x87);
        x88=2.828427124746190290949243717477656900882720947265625*x69; FIX(x88);
        x89=Max(0, x2 + x15 - x6); FIX(x89);
        x90=2.0*x89; FIX(x90);
        x91=x3*x4; FIX(x91);
        x92=x91*x85; FIX(x92);
        x93=2.828427124746190290949243717477656900882720947265625*x68; FIX(x93);
        x94=x91*x90; FIX(x94);
        x95=2.0*x3*x77; FIX(x95);

        /* Assignment result[1, 1]=-x76*x65*(x95*x68 + x95*x69 - x3*x78 + x79*x83 - x80*x79 - x80*x82 - x81*x85 + x81*x90 + x82*x83 - x87*x88 + x87*x93 - x92 + x94) */
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x10=-1.0*x6; FIX(x10);
        x15=-un; FIX(x15);
        x65=Heaviside(x6); FIX(x65);
        x67=un - 1.0*x2; FIX(x67);
        x68=Heaviside(x10 + x67); FIX(x68);
        x69=Heaviside(x6 + x67); FIX(x69);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=pow(x5, -5.0/2.0); FIX(x75);
        x76=0.25*mu*x75; FIX(x76);
        /*@ assert (x5) >= 0.;*/
        x77=pow(x5, 3.0/2.0); FIX(x77);
        x78=4.0*x77; FIX(x78);
        x79=pow(rt1, 5); FIX(x79);
        x80=1.4142135623730951454746218587388284504413604736328125*x69; FIX(x80);
        x81=pow(rt2, 4); FIX(x81);
        x82=rt1*x81; FIX(x82);
        x83=1.4142135623730951454746218587388284504413604736328125*x68; FIX(x83);
        x84=Max(0, x15 + x7); FIX(x84);
        x85=2.0*x84; FIX(x85);
        x86=pow(rt1, 3); FIX(x86);
        x87=x4*x86; FIX(x87);
        x88=2.828427124746190290949243717477656900882720947265625*x69; FIX(x88);
        x89=Max(0, x2 + x15 - x6); FIX(x89);
        x90=2.0*x89; FIX(x90);
        x91=x3*x4; FIX(x91);
        x92=x91*x85; FIX(x92);
        x93=2.828427124746190290949243717477656900882720947265625*x68; FIX(x93);
        x94=x91*x90; FIX(x94);
        x95=2.0*x3*x77; FIX(x95);
        result[4] = -x76*x65*(x95*x68 + x95*x69 - x3*x78 + x79*x83 - x80*x79 - x80*x82 - x81*x85 + x81*x90 + x82*x83 - x87*x88 + x87*x93 - x92 + x94);

    }
    else if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x36=mu*ut1; FIX(x36);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        x38=-x36*x37; FIX(x38);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x41=x17*x40; FIX(x41);
        x42=x38 + x41; FIX(x42);
        x44=x38 - x41; FIX(x44);
        x55=rt1 - x36; FIX(x55);
        x56=0.5*x39*x24; FIX(x56);
        x57=x55*x56; FIX(x57);
        x59=0.5*x39*x27; FIX(x59);
        x70=x17*x59; FIX(x70);
        x96=Max(0, x26); FIX(x96);
        x97=0.5*x96; FIX(x97);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x99=x18*x98; FIX(x99);
        x100=Max(0, x23); FIX(x100);
        x101=0.5*x100; FIX(x101);
        x102=-x40; FIX(x102);
        x103=x55*x17*x98; FIX(x103);

        /* Assignment result[1, 1]=-x101*(-mu*x103 + x102) - x42*x57 - x44*x70 - x97*(-mu*x99 + x40) */
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x36=mu*ut1; FIX(x36);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        x38=-x36*x37; FIX(x38);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x41=x17*x40; FIX(x41);
        x42=x38 + x41; FIX(x42);
        x44=x38 - x41; FIX(x44);
        x55=rt1 - x36; FIX(x55);
        x56=0.5*x39*x24; FIX(x56);
        x57=x55*x56; FIX(x57);
        x59=0.5*x39*x27; FIX(x59);
        x70=x17*x59; FIX(x70);
        x96=Max(0, x26); FIX(x96);
        x97=0.5*x96; FIX(x97);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x99=x18*x98; FIX(x99);
        x100=Max(0, x23); FIX(x100);
        x101=0.5*x100; FIX(x101);
        x102=-x40; FIX(x102);
        x103=x55*x17*x98; FIX(x103);
        result[4] = -x101*(-mu*x103 + x102) - x42*x57 - x44*x70 - x97*(-mu*x99 + x40);

    }
    else if (x74)
    {
        DEBUG_PRINT("Case (x74) is True.\n");

        /* Assignment result[1, 1]=0 */

        result[4] = 0;
        /*@ assert (result[4]) >= 0.;*/
    }
    /*@ assert (result[4]) >= 0.;*/

    /* Assignment result[2, 1]=Piecewise((-x76*(-x147*x68 - x69*x147 + x135 - x141*x133 + x141*x134 - x142*x133 + x142*x134 - x144*x88 + x144*x93 + x145*x146 - x145*x148 - x32*x79 - x32*x82 - x66*x104 + x66*x107 - x66*x109 + x66*x112 + x66*x114 + x79*x136 - x79*x140 - x79*x143 - x79*x33 + x80*x149 + x82*x136 - x82*x140 - x82*x143 - x82*x33 - x83*x149 + x87*x137 - x87*x138 - x87*x139), x14), (mu*x150 + x119 - x42*x62 - x44*x129, x72), (x43 + x45, x74)) */
    double x104;
    double x105;
    double x106;
    double x107;
    double x109;
    double x110;
    double x112;
    double x114;
    double x116;
    double x118;
    double x119;
    double x126;
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
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x1=-1.0*un; FIX(x1);
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x8=Heaviside(x1 + x7); FIX(x8);
        x10=-1.0*x6; FIX(x10);
        x11=Heaviside(x2 + x1 + x10); FIX(x11);
        x15=-un; FIX(x15);
        x32=2.0*x8; FIX(x32);
        x33=2.0*x11; FIX(x33);
        x65=Heaviside(x6); FIX(x65);
        x66=rt1*x65; FIX(x66);
        x67=un - 1.0*x2; FIX(x67);
        x68=Heaviside(x10 + x67); FIX(x68);
        x69=Heaviside(x6 + x67); FIX(x69);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=pow(x5, -5.0/2.0); FIX(x75);
        x76=0.25*mu*x75; FIX(x76);
        /*@ assert (x5) >= 0.;*/
        x77=pow(x5, 3.0/2.0); FIX(x77);
        x78=4.0*x77; FIX(x78);
        x79=pow(rt1, 5); FIX(x79);
        x80=1.4142135623730951454746218587388284504413604736328125*x69; FIX(x80);
        x81=pow(rt2, 4); FIX(x81);
        x82=rt1*x81; FIX(x82);
        x83=1.4142135623730951454746218587388284504413604736328125*x68; FIX(x83);
        x84=Max(0, x15 + x7); FIX(x84);
        x85=2.0*x84; FIX(x85);
        x86=pow(rt1, 3); FIX(x86);
        x87=x4*x86; FIX(x87);
        x88=2.828427124746190290949243717477656900882720947265625*x69; FIX(x88);
        x89=Max(0, x2 + x15 - x6); FIX(x89);
        x90=2.0*x89; FIX(x90);
        x93=2.828427124746190290949243717477656900882720947265625*x68; FIX(x93);
        x104=rt2*x78; FIX(x104);
        x105=pow(rt1, 4); FIX(x105);
        x106=pow(rt2, 3); FIX(x106);
        x107=x106*x85; FIX(x107);
        x109=x106*x90; FIX(x109);
        x110=2.0*x77; FIX(x110);
        x111=rt2*x69; FIX(x111);
        x112=x110*x111; FIX(x112);
        x113=rt2*x68; FIX(x113);
        x114=x110*x113; FIX(x114);
        x130=1.4142135623730951454746218587388284504413604736328125*x8*x77; FIX(x130);
        x131=1.4142135623730951454746218587388284504413604736328125*x77*x11; FIX(x131);
        x132=x3*x77; FIX(x132);
        x133=1.4142135623730951454746218587388284504413604736328125*x65*x69; FIX(x133);
        x134=1.4142135623730951454746218587388284504413604736328125*x65*x68; FIX(x134);
        x135=x132*x133 - x132*x134 - x3*x130 + x3*x131 - x4*x130 + x4*x131; FIX(x135);
        x136=4.0*x65; FIX(x136);
        x137=8.0*x65; FIX(x137);
        x138=4.0*x8; FIX(x138);
        x139=4.0*x11; FIX(x139);
        x140=2.0*x65*x69; FIX(x140);
        x141=pow(rt2, 5); FIX(x141);
        x142=rt2*x105; FIX(x142);
        x143=2.0*x65*x68; FIX(x143);
        x144=x3*x106*x65; FIX(x144);
        x145=rt2*x86; FIX(x145);
        x146=2.0*x65*x84; FIX(x146);
        x147=4.0*x4*x86*x65; FIX(x147);
        x148=2.0*x65*x89; FIX(x148);
        x149=x4*x77*x65; FIX(x149);

        /* Assignment result[2, 1]=-x76*(-x147*x68 - x69*x147 + x135 - x141*x133 + x141*x134 - x142*x133 + x142*x134 - x144*x88 + x144*x93 + x145*x146 - x145*x148 - x32*x79 - x32*x82 - x66*x104 + x66*x107 - x66*x109 + x66*x112 + x66*x114 + x79*x136 - x79*x140 - x79*x143 - x79*x33 + x80*x149 + x82*x136 - x82*x140 - x82*x143 - x82*x33 - x83*x149 + x87*x137 - x87*x138 - x87*x139) */
        x1=-1.0*un; FIX(x1);
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x8=Heaviside(x1 + x7); FIX(x8);
        x10=-1.0*x6; FIX(x10);
        x11=Heaviside(x2 + x1 + x10); FIX(x11);
        x15=-un; FIX(x15);
        x32=2.0*x8; FIX(x32);
        x33=2.0*x11; FIX(x33);
        x65=Heaviside(x6); FIX(x65);
        x66=rt1*x65; FIX(x66);
        x67=un - 1.0*x2; FIX(x67);
        x68=Heaviside(x10 + x67); FIX(x68);
        x69=Heaviside(x6 + x67); FIX(x69);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=pow(x5, -5.0/2.0); FIX(x75);
        x76=0.25*mu*x75; FIX(x76);
        /*@ assert (x5) >= 0.;*/
        x77=pow(x5, 3.0/2.0); FIX(x77);
        x78=4.0*x77; FIX(x78);
        x79=pow(rt1, 5); FIX(x79);
        x80=1.4142135623730951454746218587388284504413604736328125*x69; FIX(x80);
        x81=pow(rt2, 4); FIX(x81);
        x82=rt1*x81; FIX(x82);
        x83=1.4142135623730951454746218587388284504413604736328125*x68; FIX(x83);
        x84=Max(0, x15 + x7); FIX(x84);
        x85=2.0*x84; FIX(x85);
        x86=pow(rt1, 3); FIX(x86);
        x87=x4*x86; FIX(x87);
        x88=2.828427124746190290949243717477656900882720947265625*x69; FIX(x88);
        x89=Max(0, x2 + x15 - x6); FIX(x89);
        x90=2.0*x89; FIX(x90);
        x93=2.828427124746190290949243717477656900882720947265625*x68; FIX(x93);
        x104=rt2*x78; FIX(x104);
        x105=pow(rt1, 4); FIX(x105);
        x106=pow(rt2, 3); FIX(x106);
        x107=x106*x85; FIX(x107);
        x109=x106*x90; FIX(x109);
        x110=2.0*x77; FIX(x110);
        x111=rt2*x69; FIX(x111);
        x112=x110*x111; FIX(x112);
        x113=rt2*x68; FIX(x113);
        x114=x110*x113; FIX(x114);
        x130=1.4142135623730951454746218587388284504413604736328125*x8*x77; FIX(x130);
        x131=1.4142135623730951454746218587388284504413604736328125*x77*x11; FIX(x131);
        x132=x3*x77; FIX(x132);
        x133=1.4142135623730951454746218587388284504413604736328125*x65*x69; FIX(x133);
        x134=1.4142135623730951454746218587388284504413604736328125*x65*x68; FIX(x134);
        x135=x132*x133 - x132*x134 - x3*x130 + x3*x131 - x4*x130 + x4*x131; FIX(x135);
        x136=4.0*x65; FIX(x136);
        x137=8.0*x65; FIX(x137);
        x138=4.0*x8; FIX(x138);
        x139=4.0*x11; FIX(x139);
        x140=2.0*x65*x69; FIX(x140);
        x141=pow(rt2, 5); FIX(x141);
        x142=rt2*x105; FIX(x142);
        x143=2.0*x65*x68; FIX(x143);
        x144=x3*x106*x65; FIX(x144);
        x145=rt2*x86; FIX(x145);
        x146=2.0*x65*x84; FIX(x146);
        x147=4.0*x4*x86*x65; FIX(x147);
        x148=2.0*x65*x89; FIX(x148);
        x149=x4*x77*x65; FIX(x149);
        result[5] = -x76*(-x147*x68 - x69*x147 + x135 - x141*x133 + x141*x134 - x142*x133 + x142*x134 - x144*x88 + x144*x93 + x145*x146 - x145*x148 - x32*x79 - x32*x82 - x66*x104 + x66*x107 - x66*x109 + x66*x112 + x66*x114 + x79*x136 - x79*x140 - x79*x143 - x79*x33 + x80*x149 + x82*x136 - x82*x140 - x82*x143 - x82*x33 - x83*x149 + x87*x137 - x87*x138 - x87*x139);

    }
    else if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x36=mu*ut1; FIX(x36);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        x38=-x36*x37; FIX(x38);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x41=x17*x40; FIX(x41);
        x42=x38 + x41; FIX(x42);
        x44=x38 - x41; FIX(x44);
        x46=mu*ut2; FIX(x46);
        x56=0.5*x39*x24; FIX(x56);
        x59=0.5*x39*x27; FIX(x59);
        x61=rt2 - x46; FIX(x61);
        x62=x61*x56; FIX(x62);
        x96=Max(0, x26); FIX(x96);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x100=Max(0, x23); FIX(x100);
        x116=0.5*x98*x100; FIX(x116);
        x118=0.5*x17*x19*x96*x98; FIX(x118);
        x119=mu*x118; FIX(x119);
        x126=x61*x17; FIX(x126);
        x129=x19*x59; FIX(x129);
        x150=x126*x116; FIX(x150);

        /* Assignment result[2, 1]=mu*x150 + x119 - x42*x62 - x44*x129 */
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x36=mu*ut1; FIX(x36);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        x38=-x36*x37; FIX(x38);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x41=x17*x40; FIX(x41);
        x42=x38 + x41; FIX(x42);
        x44=x38 - x41; FIX(x44);
        x46=mu*ut2; FIX(x46);
        x56=0.5*x39*x24; FIX(x56);
        x59=0.5*x39*x27; FIX(x59);
        x61=rt2 - x46; FIX(x61);
        x62=x61*x56; FIX(x62);
        x96=Max(0, x26); FIX(x96);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x100=Max(0, x23); FIX(x100);
        x116=0.5*x98*x100; FIX(x116);
        x118=0.5*x17*x19*x96*x98; FIX(x118);
        x119=mu*x118; FIX(x119);
        x126=x61*x17; FIX(x126);
        x129=x19*x59; FIX(x129);
        x150=x126*x116; FIX(x150);
        result[5] = mu*x150 + x119 - x42*x62 - x44*x129;

    }
    else if (x74)
    {
        DEBUG_PRINT("Case (x74) is True.\n");
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x25=0.5*x24; FIX(x25);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x28=0.5*x27; FIX(x28);
        x36=mu*ut1; FIX(x36);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        x38=-x36*x37; FIX(x38);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x41=x17*x40; FIX(x41);
        x42=x38 + x41; FIX(x42);
        x43=-x25*x42; FIX(x43);
        x44=x38 - x41; FIX(x44);
        x45=x44*x28; FIX(x45);

        /* Assignment result[2, 1]=x43 + x45 */
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x25=0.5*x24; FIX(x25);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x28=0.5*x27; FIX(x28);
        x36=mu*ut1; FIX(x36);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        x38=-x36*x37; FIX(x38);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x41=x17*x40; FIX(x41);
        x42=x38 + x41; FIX(x42);
        x43=-x25*x42; FIX(x43);
        x44=x38 - x41; FIX(x44);
        x45=x44*x28; FIX(x45);
        result[5] = x43 + x45;

    }
    /*@ assert (result[5]) >= 0.;*/

    /* Assignment result[0, 2]=Piecewise((x31*(rt2*x32 - rt2*x33 + x35), x14), (x50 - x52, x29)) */
    double x47;
    double x48;
    double x49;
    double x50;
    double x51;
    double x52;
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x1=-1.0*un; FIX(x1);
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x8=Heaviside(x1 + x7); FIX(x8);
        x10=-1.0*x6; FIX(x10);
        x11=Heaviside(x2 + x1 + x10); FIX(x11);
        /*@ assert (x6) != 0.;*/
        x30=1.0/x6; FIX(x30);
        x31=0.25*mu*x30; FIX(x31);
        x32=2.0*x8; FIX(x32);
        x33=2.0*x11; FIX(x33);
        x34=1.4142135623730951454746218587388284504413604736328125*x6; FIX(x34);
        x35=x34*x11 + x8*x34; FIX(x35);

        /* Assignment result[0, 2]=x31*(rt2*x32 - rt2*x33 + x35) */
        x1=-1.0*un; FIX(x1);
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x8=Heaviside(x1 + x7); FIX(x8);
        x10=-1.0*x6; FIX(x10);
        x11=Heaviside(x2 + x1 + x10); FIX(x11);
        /*@ assert (x6) != 0.;*/
        x30=1.0/x6; FIX(x30);
        x31=0.25*mu*x30; FIX(x31);
        x32=2.0*x8; FIX(x32);
        x33=2.0*x11; FIX(x33);
        x34=1.4142135623730951454746218587388284504413604736328125*x6; FIX(x34);
        x35=x34*x11 + x8*x34; FIX(x35);
        result[6] = x31*(rt2*x32 - rt2*x33 + x35);

    }
    else if (x29)
    {
        DEBUG_PRINT("Case (x29) is True.\n");
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x25=0.5*x24; FIX(x25);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x28=0.5*x27; FIX(x28);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x46=mu*ut2; FIX(x46);
        x47=-x46*x37; FIX(x47);
        x48=x19*x40; FIX(x48);
        x49=x47 + x48; FIX(x49);
        x50=-x25*x49; FIX(x50);
        x51=x47 - x48; FIX(x51);
        x52=x51*x28; FIX(x52);

        /* Assignment result[0, 2]=x50 - x52 */
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x25=0.5*x24; FIX(x25);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x28=0.5*x27; FIX(x28);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x46=mu*ut2; FIX(x46);
        x47=-x46*x37; FIX(x47);
        x48=x19*x40; FIX(x48);
        x49=x47 + x48; FIX(x49);
        x50=-x25*x49; FIX(x50);
        x51=x47 - x48; FIX(x51);
        x52=x51*x28; FIX(x52);
        result[6] = x50 - x52;

    }
    /*@ assert (result[6]) >= 0.;*/

    /* Assignment result[1, 2]=Piecewise((0.25*mu*rt1*x75*x65*(x104 - x105*x83 - x107 - x108*x85 + x108*x90 + x109 - x112 - x114 + x80*x105 + x80*x81 - x81*x83 + x91*x88 - x91*x93), x14), (mu*x117 + x119 - x49*x57 - x51*x70, x72), (0, x74)) */
    double x108;
    double x115;
    double x117;
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x10=-1.0*x6; FIX(x10);
        x15=-un; FIX(x15);
        x65=Heaviside(x6); FIX(x65);
        x67=un - 1.0*x2; FIX(x67);
        x68=Heaviside(x10 + x67); FIX(x68);
        x69=Heaviside(x6 + x67); FIX(x69);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=pow(x5, -5.0/2.0); FIX(x75);
        /*@ assert (x5) >= 0.;*/
        x77=pow(x5, 3.0/2.0); FIX(x77);
        x78=4.0*x77; FIX(x78);
        x80=1.4142135623730951454746218587388284504413604736328125*x69; FIX(x80);
        x81=pow(rt2, 4); FIX(x81);
        x83=1.4142135623730951454746218587388284504413604736328125*x68; FIX(x83);
        x84=Max(0, x15 + x7); FIX(x84);
        x85=2.0*x84; FIX(x85);
        x88=2.828427124746190290949243717477656900882720947265625*x69; FIX(x88);
        x89=Max(0, x2 + x15 - x6); FIX(x89);
        x90=2.0*x89; FIX(x90);
        x91=x3*x4; FIX(x91);
        x93=2.828427124746190290949243717477656900882720947265625*x68; FIX(x93);
        x104=rt2*x78; FIX(x104);
        x105=pow(rt1, 4); FIX(x105);
        x106=pow(rt2, 3); FIX(x106);
        x107=x106*x85; FIX(x107);
        x108=rt2*x3; FIX(x108);
        x109=x106*x90; FIX(x109);
        x110=2.0*x77; FIX(x110);
        x111=rt2*x69; FIX(x111);
        x112=x110*x111; FIX(x112);
        x113=rt2*x68; FIX(x113);
        x114=x110*x113; FIX(x114);

        /* Assignment result[1, 2]=0.25*mu*rt1*x75*x65*(x104 - x105*x83 - x107 - x108*x85 + x108*x90 + x109 - x112 - x114 + x80*x105 + x80*x81 - x81*x83 + x91*x88 - x91*x93) */
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x10=-1.0*x6; FIX(x10);
        x15=-un; FIX(x15);
        x65=Heaviside(x6); FIX(x65);
        x67=un - 1.0*x2; FIX(x67);
        x68=Heaviside(x10 + x67); FIX(x68);
        x69=Heaviside(x6 + x67); FIX(x69);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=pow(x5, -5.0/2.0); FIX(x75);
        /*@ assert (x5) >= 0.;*/
        x77=pow(x5, 3.0/2.0); FIX(x77);
        x78=4.0*x77; FIX(x78);
        x80=1.4142135623730951454746218587388284504413604736328125*x69; FIX(x80);
        x81=pow(rt2, 4); FIX(x81);
        x83=1.4142135623730951454746218587388284504413604736328125*x68; FIX(x83);
        x84=Max(0, x15 + x7); FIX(x84);
        x85=2.0*x84; FIX(x85);
        x88=2.828427124746190290949243717477656900882720947265625*x69; FIX(x88);
        x89=Max(0, x2 + x15 - x6); FIX(x89);
        x90=2.0*x89; FIX(x90);
        x91=x3*x4; FIX(x91);
        x93=2.828427124746190290949243717477656900882720947265625*x68; FIX(x93);
        x104=rt2*x78; FIX(x104);
        x105=pow(rt1, 4); FIX(x105);
        x106=pow(rt2, 3); FIX(x106);
        x107=x106*x85; FIX(x107);
        x108=rt2*x3; FIX(x108);
        x109=x106*x90; FIX(x109);
        x110=2.0*x77; FIX(x110);
        x111=rt2*x69; FIX(x111);
        x112=x110*x111; FIX(x112);
        x113=rt2*x68; FIX(x113);
        x114=x110*x113; FIX(x114);
        result[7] = 0.25*mu*rt1*x75*x65*(x104 - x105*x83 - x107 - x108*x85 + x108*x90 + x109 - x112 - x114 + x80*x105 + x80*x81 - x81*x83 + x91*x88 - x91*x93);

    }
    else if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x36=mu*ut1; FIX(x36);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x46=mu*ut2; FIX(x46);
        x47=-x46*x37; FIX(x47);
        x48=x19*x40; FIX(x48);
        x49=x47 + x48; FIX(x49);
        x51=x47 - x48; FIX(x51);
        x55=rt1 - x36; FIX(x55);
        x56=0.5*x39*x24; FIX(x56);
        x57=x55*x56; FIX(x57);
        x59=0.5*x39*x27; FIX(x59);
        x70=x17*x59; FIX(x70);
        x96=Max(0, x26); FIX(x96);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x100=Max(0, x23); FIX(x100);
        x115=x55*x19; FIX(x115);
        x116=0.5*x98*x100; FIX(x116);
        x117=x115*x116; FIX(x117);
        x118=0.5*x17*x19*x96*x98; FIX(x118);
        x119=mu*x118; FIX(x119);

        /* Assignment result[1, 2]=mu*x117 + x119 - x49*x57 - x51*x70 */
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x36=mu*ut1; FIX(x36);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x46=mu*ut2; FIX(x46);
        x47=-x46*x37; FIX(x47);
        x48=x19*x40; FIX(x48);
        x49=x47 + x48; FIX(x49);
        x51=x47 - x48; FIX(x51);
        x55=rt1 - x36; FIX(x55);
        x56=0.5*x39*x24; FIX(x56);
        x57=x55*x56; FIX(x57);
        x59=0.5*x39*x27; FIX(x59);
        x70=x17*x59; FIX(x70);
        x96=Max(0, x26); FIX(x96);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x100=Max(0, x23); FIX(x100);
        x115=x55*x19; FIX(x115);
        x116=0.5*x98*x100; FIX(x116);
        x117=x115*x116; FIX(x117);
        x118=0.5*x17*x19*x96*x98; FIX(x118);
        x119=mu*x118; FIX(x119);
        result[7] = mu*x117 + x119 - x49*x57 - x51*x70;

    }
    else if (x74)
    {
        DEBUG_PRINT("Case (x74) is True.\n");

        /* Assignment result[1, 2]=0 */

        result[7] = 0;
        /*@ assert (result[7]) >= 0.;*/
    }
    /*@ assert (result[7]) >= 0.;*/

    /* Assignment result[2, 2]=Piecewise((-x76*(-1.17157287525381*x3*x106*x65*x68 - 6.82842712474619*x3*x106*x65*x69 - x4*x78*x65 - x65*x92 + x65*x94 + 0.585786437626905*x68*x149 + 3.41421356237309*x69*x149 - x105*x146 + x105*x148 + x135 + x141*x136 - x141*x152 - x141*x153 - x141*x33 + x142*x136 - x142*x152 - x142*x153 - x142*x33 + x151*x137 - x151*x138 - x151*x139 - x32*x141 - x32*x142), x14), (-x101*(-mu*x155 + x102) - x49*x62 - x51*x129 - x97*(-mu*x154 + x40), x72), (x50 + x52, x74)) */
    double x151;
    double x152;
    double x153;
    double x154;
    double x155;
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");
        x1=-1.0*un; FIX(x1);
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x8=Heaviside(x1 + x7); FIX(x8);
        x10=-1.0*x6; FIX(x10);
        x11=Heaviside(x2 + x1 + x10); FIX(x11);
        x15=-un; FIX(x15);
        x32=2.0*x8; FIX(x32);
        x33=2.0*x11; FIX(x33);
        x65=Heaviside(x6); FIX(x65);
        x67=un - 1.0*x2; FIX(x67);
        x68=Heaviside(x10 + x67); FIX(x68);
        x69=Heaviside(x6 + x67); FIX(x69);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=pow(x5, -5.0/2.0); FIX(x75);
        x76=0.25*mu*x75; FIX(x76);
        /*@ assert (x5) >= 0.;*/
        x77=pow(x5, 3.0/2.0); FIX(x77);
        x78=4.0*x77; FIX(x78);
        x84=Max(0, x15 + x7); FIX(x84);
        x85=2.0*x84; FIX(x85);
        x89=Max(0, x2 + x15 - x6); FIX(x89);
        x90=2.0*x89; FIX(x90);
        x91=x3*x4; FIX(x91);
        x92=x91*x85; FIX(x92);
        x94=x91*x90; FIX(x94);
        x105=pow(rt1, 4); FIX(x105);
        x106=pow(rt2, 3); FIX(x106);
        x130=1.4142135623730951454746218587388284504413604736328125*x8*x77; FIX(x130);
        x131=1.4142135623730951454746218587388284504413604736328125*x77*x11; FIX(x131);
        x132=x3*x77; FIX(x132);
        x133=1.4142135623730951454746218587388284504413604736328125*x65*x69; FIX(x133);
        x134=1.4142135623730951454746218587388284504413604736328125*x65*x68; FIX(x134);
        x135=x132*x133 - x132*x134 - x3*x130 + x3*x131 - x4*x130 + x4*x131; FIX(x135);
        x136=4.0*x65; FIX(x136);
        x137=8.0*x65; FIX(x137);
        x138=4.0*x8; FIX(x138);
        x139=4.0*x11; FIX(x139);
        x141=pow(rt2, 5); FIX(x141);
        x142=rt2*x105; FIX(x142);
        x146=2.0*x65*x84; FIX(x146);
        x148=2.0*x65*x89; FIX(x148);
        x149=x4*x77*x65; FIX(x149);
        x151=x3*x106; FIX(x151);
        x152=3.41421356237309492343001693370752036571502685546875*x65*x69; FIX(x152);
        x153=0.58578643762690496554768060377682559192180633544921875*x65*x68; FIX(x153);

        /* Assignment result[2, 2]=-x76*(-1.17157287525381*x3*x106*x65*x68 - 6.82842712474619*x3*x106*x65*x69 - x4*x78*x65 - x65*x92 + x65*x94 + 0.585786437626905*x68*x149 + 3.41421356237309*x69*x149 - x105*x146 + x105*x148 + x135 + x141*x136 - x141*x152 - x141*x153 - x141*x33 + x142*x136 - x142*x152 - x142*x153 - x142*x33 + x151*x137 - x151*x138 - x151*x139 - x32*x141 - x32*x142) */
        x1=-1.0*un; FIX(x1);
        x2=mu*rn; FIX(x2);
        x3=rt1*rt1; FIX(x3);
        x4=rt2*rt2; FIX(x4);
        x5=x3 + x4; FIX(x5);
        /*@ assert (x5) >= 0.;*/
        x6=sqrt(x5); FIX(x6);
        x7=x2 + x6; FIX(x7);
        x8=Heaviside(x1 + x7); FIX(x8);
        x10=-1.0*x6; FIX(x10);
        x11=Heaviside(x2 + x1 + x10); FIX(x11);
        x15=-un; FIX(x15);
        x32=2.0*x8; FIX(x32);
        x33=2.0*x11; FIX(x33);
        x65=Heaviside(x6); FIX(x65);
        x67=un - 1.0*x2; FIX(x67);
        x68=Heaviside(x10 + x67); FIX(x68);
        x69=Heaviside(x6 + x67); FIX(x69);
        /*@ assert (x5) >= 0.;*/
        /*@ assert (x5) != 0.;*/
        x75=pow(x5, -5.0/2.0); FIX(x75);
        x76=0.25*mu*x75; FIX(x76);
        /*@ assert (x5) >= 0.;*/
        x77=pow(x5, 3.0/2.0); FIX(x77);
        x78=4.0*x77; FIX(x78);
        x84=Max(0, x15 + x7); FIX(x84);
        x85=2.0*x84; FIX(x85);
        x89=Max(0, x2 + x15 - x6); FIX(x89);
        x90=2.0*x89; FIX(x90);
        x91=x3*x4; FIX(x91);
        x92=x91*x85; FIX(x92);
        x94=x91*x90; FIX(x94);
        x105=pow(rt1, 4); FIX(x105);
        x106=pow(rt2, 3); FIX(x106);
        x130=1.4142135623730951454746218587388284504413604736328125*x8*x77; FIX(x130);
        x131=1.4142135623730951454746218587388284504413604736328125*x77*x11; FIX(x131);
        x132=x3*x77; FIX(x132);
        x133=1.4142135623730951454746218587388284504413604736328125*x65*x69; FIX(x133);
        x134=1.4142135623730951454746218587388284504413604736328125*x65*x68; FIX(x134);
        x135=x132*x133 - x132*x134 - x3*x130 + x3*x131 - x4*x130 + x4*x131; FIX(x135);
        x136=4.0*x65; FIX(x136);
        x137=8.0*x65; FIX(x137);
        x138=4.0*x8; FIX(x138);
        x139=4.0*x11; FIX(x139);
        x141=pow(rt2, 5); FIX(x141);
        x142=rt2*x105; FIX(x142);
        x146=2.0*x65*x84; FIX(x146);
        x148=2.0*x65*x89; FIX(x148);
        x149=x4*x77*x65; FIX(x149);
        x151=x3*x106; FIX(x151);
        x152=3.41421356237309492343001693370752036571502685546875*x65*x69; FIX(x152);
        x153=0.58578643762690496554768060377682559192180633544921875*x65*x68; FIX(x153);
        result[8] = -x76*(-1.1715728752538099310953612075536511838436126708984375*x3*x106*x65*x68 - 6.8284271247461898468600338674150407314300537109375*x3*x106*x65*x69 - x4*x78*x65 - x65*x92 + x65*x94 + 0.58578643762690496554768060377682559192180633544921875*x68*x149 + 3.41421356237309492343001693370752036571502685546875*x69*x149 - x105*x146 + x105*x148 + x135 + x141*x136 - x141*x152 - x141*x153 - x141*x33 + x142*x136 - x142*x152 - x142*x153 - x142*x33 + x151*x137 - x151*x138 - x151*x139 - x32*x141 - x32*x142);

    }
    else if (x72)
    {
        DEBUG_PRINT("Case (x72) is True.\n");
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x46=mu*ut2; FIX(x46);
        x47=-x46*x37; FIX(x47);
        x48=x19*x40; FIX(x48);
        x49=x47 + x48; FIX(x49);
        x51=x47 - x48; FIX(x51);
        x56=0.5*x39*x24; FIX(x56);
        x59=0.5*x39*x27; FIX(x59);
        x61=rt2 - x46; FIX(x61);
        x62=x61*x56; FIX(x62);
        x96=Max(0, x26); FIX(x96);
        x97=0.5*x96; FIX(x97);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x100=Max(0, x23); FIX(x100);
        x101=0.5*x100; FIX(x101);
        x102=-x40; FIX(x102);
        x129=x19*x59; FIX(x129);
        x154=x20*x98; FIX(x154);
        x155=x61*x19*x98; FIX(x155);

        /* Assignment result[2, 2]=-x101*(-mu*x155 + x102) - x49*x62 - x51*x129 - x97*(-mu*x154 + x40) */
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x46=mu*ut2; FIX(x46);
        x47=-x46*x37; FIX(x47);
        x48=x19*x40; FIX(x48);
        x49=x47 + x48; FIX(x49);
        x51=x47 - x48; FIX(x51);
        x56=0.5*x39*x24; FIX(x56);
        x59=0.5*x39*x27; FIX(x59);
        x61=rt2 - x46; FIX(x61);
        x62=x61*x56; FIX(x62);
        x96=Max(0, x26); FIX(x96);
        x97=0.5*x96; FIX(x97);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x100=Max(0, x23); FIX(x100);
        x101=0.5*x100; FIX(x101);
        x102=-x40; FIX(x102);
        x129=x19*x59; FIX(x129);
        x154=x20*x98; FIX(x154);
        x155=x61*x19*x98; FIX(x155);
        result[8] = -x101*(-mu*x155 + x102) - x49*x62 - x51*x129 - x97*(-mu*x154 + x40);

    }
    else if (x74)
    {
        DEBUG_PRINT("Case (x74) is True.\n");
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x25=0.5*x24; FIX(x25);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x28=0.5*x27; FIX(x28);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x46=mu*ut2; FIX(x46);
        x47=-x46*x37; FIX(x47);
        x48=x19*x40; FIX(x48);
        x49=x47 + x48; FIX(x49);
        x50=-x25*x49; FIX(x50);
        x51=x47 - x48; FIX(x51);
        x52=x51*x28; FIX(x52);

        /* Assignment result[2, 2]=x50 + x52 */
        x2=mu*rn; FIX(x2);
        x15=-un; FIX(x15);
        x16=-mu*x13 + x2 + x15; FIX(x16);
        x23=x16 + x22; FIX(x23);
        x24=Heaviside(x23); FIX(x24);
        x25=0.5*x24; FIX(x25);
        x26=x16 - x22; FIX(x26);
        x27=Heaviside(x26); FIX(x27);
        x28=0.5*x27; FIX(x28);
        /*@ assert (x13) != 0.;*/
        x37=1.0/x13; FIX(x37);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x46=mu*ut2; FIX(x46);
        x47=-x46*x37; FIX(x47);
        x48=x19*x40; FIX(x48);
        x49=x47 + x48; FIX(x49);
        x50=-x25*x49; FIX(x50);
        x51=x47 - x48; FIX(x51);
        x52=x51*x28; FIX(x52);
        result[8] = x50 + x52;

    }
    /*@ assert (result[8]) >= 0.;*/

    /* Assignment result[0, 3]=mu + x53 - x54 */
    double x53;
    double x54;x2=mu*rn; FIX(x2);
    x15=-un; FIX(x15);
    x16=-mu*x13 + x2 + x15; FIX(x16);
    x23=x16 + x22; FIX(x23);
    x24=Heaviside(x23); FIX(x24);
    x25=0.5*x24; FIX(x25);
    x26=x16 - x22; FIX(x26);
    x27=Heaviside(x26); FIX(x27);
    x28=0.5*x27; FIX(x28);
    x53=-mu*x25; FIX(x53);
    x54=mu*x28; FIX(x54);
    result[9] = mu + x53 - x54;


    /* Assignment result[1, 3]=Piecewise((-x41*x28 - x55*x120, x71), (0, x73)) */

    if (x71)
    {
        DEBUG_PRINT("Case (x71) is True.\n");
        x36=mu*ut1; FIX(x36);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x41=x17*x40; FIX(x41);
        x55=rt1 - x36; FIX(x55);
        x120=0.5*mu*x39*x24; FIX(x120);

        /* Assignment result[1, 3]=-x41*x28 - x55*x120 */
        x36=mu*ut1; FIX(x36);
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x41=x17*x40; FIX(x41);
        x55=rt1 - x36; FIX(x55);
        x120=0.5*mu*x39*x24; FIX(x120);
        result[10] = -x41*x28 - x55*x120;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[1, 3]=0 */

        result[10] = 0;
        /*@ assert (result[10]) >= 0.;*/
    }
    /*@ assert (result[10]) >= 0.;*/

    /* Assignment result[2, 3]=Piecewise((-x48*x28 - x61*x120, x71), (x53 + x54, x73)) */

    if (x71)
    {
        DEBUG_PRINT("Case (x71) is True.\n");
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x46=mu*ut2; FIX(x46);
        x48=x19*x40; FIX(x48);
        x61=rt2 - x46; FIX(x61);
        x120=0.5*mu*x39*x24; FIX(x120);

        /* Assignment result[2, 3]=-x48*x28 - x61*x120 */
        /*@ assert (x22) != 0.;*/
        x39=1.0/x22; FIX(x39);
        x40=mu*x39; FIX(x40);
        x46=mu*ut2; FIX(x46);
        x48=x19*x40; FIX(x48);
        x61=rt2 - x46; FIX(x61);
        x120=0.5*mu*x39*x24; FIX(x120);
        result[11] = -x48*x28 - x61*x120;

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
    double x60;x36=mu*ut1; FIX(x36);
    /*@ assert (x22) != 0.;*/
    x39=1.0/x22; FIX(x39);
    x55=rt1 - x36; FIX(x55);
    x56=0.5*x39*x24; FIX(x56);
    x57=x55*x56; FIX(x57);
    x58=-x57; FIX(x58);
    x59=0.5*x39*x27; FIX(x59);
    x60=x55*x59; FIX(x60);
    result[12] = x58 + x60;


    /* Assignment result[1, 4]=Piecewise((x55*x17*x123 + 1 - x101*(x103 + x39) - x55**2*x122 - x97*(x124 + x99), x71), (1, x73)) */
    double x121;
    double x122;
    double x123;
    double x124;
    if (x71)
    {
        DEBUG_PRINT("Case (x71) is True.\n");
        x96=Max(0, x26); FIX(x96);
        x97=0.5*x96; FIX(x97);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x99=x18*x98; FIX(x99);
        x100=Max(0, x23); FIX(x100);
        x101=0.5*x100; FIX(x101);
        x103=x55*x17*x98; FIX(x103);
        /*@ assert (x21) != 0.;*/
        x121=1.0/x21; FIX(x121);
        x122=0.5*x121*x24; FIX(x122);
        x123=0.5*x121*x27; FIX(x123);
        x124=-x39; FIX(x124);

        /* Assignment result[1, 4]=x55*x17*x123 + 1 - x101*(x103 + x39) - x55**2*x122 - x97*(x124 + x99) */
        x96=Max(0, x26); FIX(x96);
        x97=0.5*x96; FIX(x97);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x99=x18*x98; FIX(x99);
        x100=Max(0, x23); FIX(x100);
        x101=0.5*x100; FIX(x101);
        x103=x55*x17*x98; FIX(x103);
        /*@ assert (x21) != 0.;*/
        x121=1.0/x21; FIX(x121);
        x122=0.5*x121*x24; FIX(x122);
        x123=0.5*x121*x27; FIX(x123);
        x124=-x39; FIX(x124);
        result[13] = x55*x17*x123 + 1 - x101*(x103 + x39) - x55*x55*x122 - x97*(x124 + x99);

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

    /* Assignment result[2, 4]=Piecewise((x115*x123 + x125 - x150, x71), (x58 - x60, x73)) */
    double x125;
    if (x71)
    {
        DEBUG_PRINT("Case (x71) is True.\n");
        x46=mu*ut2; FIX(x46);
        x61=rt2 - x46; FIX(x61);
        x96=Max(0, x26); FIX(x96);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x100=Max(0, x23); FIX(x100);
        x115=x55*x19; FIX(x115);
        x116=0.5*x98*x100; FIX(x116);
        x118=0.5*x17*x19*x96*x98; FIX(x118);
        /*@ assert (x21) != 0.;*/
        x121=1.0/x21; FIX(x121);
        x122=0.5*x121*x24; FIX(x122);
        x123=0.5*x121*x27; FIX(x123);
        x125=-x118 - x55*x61*x122; FIX(x125);
        x126=x61*x17; FIX(x126);
        x150=x126*x116; FIX(x150);

        /* Assignment result[2, 4]=x115*x123 + x125 - x150 */
        x46=mu*ut2; FIX(x46);
        x61=rt2 - x46; FIX(x61);
        x96=Max(0, x26); FIX(x96);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x100=Max(0, x23); FIX(x100);
        x115=x55*x19; FIX(x115);
        x116=0.5*x98*x100; FIX(x116);
        x118=0.5*x17*x19*x96*x98; FIX(x118);
        /*@ assert (x21) != 0.;*/
        x121=1.0/x21; FIX(x121);
        x122=0.5*x121*x24; FIX(x122);
        x123=0.5*x121*x27; FIX(x123);
        x125=-x118 - x55*x61*x122; FIX(x125);
        x126=x61*x17; FIX(x126);
        x150=x126*x116; FIX(x150);
        result[14] = x115*x123 + x125 - x150;

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
    double x64;x46=mu*ut2; FIX(x46);
    x61=rt2 - x46; FIX(x61);
    x62=x61*x56; FIX(x62);
    x63=-x62; FIX(x63);
    x64=x61*x59; FIX(x64);
    result[15] = x63 + x64;


    /* Assignment result[1, 5]=Piecewise((-x117 + x125 + x126*x123, x71), (0, x73)) */

    if (x71)
    {
        DEBUG_PRINT("Case (x71) is True.\n");
        x96=Max(0, x26); FIX(x96);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x100=Max(0, x23); FIX(x100);
        x115=x55*x19; FIX(x115);
        x116=0.5*x98*x100; FIX(x116);
        x117=x115*x116; FIX(x117);
        x118=0.5*x17*x19*x96*x98; FIX(x118);
        /*@ assert (x21) != 0.;*/
        x121=1.0/x21; FIX(x121);
        x122=0.5*x121*x24; FIX(x122);
        x123=0.5*x121*x27; FIX(x123);
        x125=-x118 - x55*x61*x122; FIX(x125);
        x126=x61*x17; FIX(x126);

        /* Assignment result[1, 5]=-x117 + x125 + x126*x123 */
        x96=Max(0, x26); FIX(x96);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x100=Max(0, x23); FIX(x100);
        x115=x55*x19; FIX(x115);
        x116=0.5*x98*x100; FIX(x116);
        x117=x115*x116; FIX(x117);
        x118=0.5*x17*x19*x96*x98; FIX(x118);
        /*@ assert (x21) != 0.;*/
        x121=1.0/x21; FIX(x121);
        x122=0.5*x121*x24; FIX(x122);
        x123=0.5*x121*x27; FIX(x123);
        x125=-x118 - x55*x61*x122; FIX(x125);
        x126=x61*x17; FIX(x126);
        result[16] = -x117 + x125 + x126*x123;

    }
    else if (x73)
    {
        DEBUG_PRINT("Case (x73) is True.\n");

        /* Assignment result[1, 5]=0 */

        result[16] = 0;
        /*@ assert (result[16]) >= 0.;*/
    }
    /*@ assert (result[16]) >= 0.;*/

    /* Assignment result[2, 5]=Piecewise((x61*x19*x123 + 1 - x101*(x155 + x39) - x61**2*x122 - x97*(x124 + x154), x71), (1 + x63 - x64, x73)) */

    if (x71)
    {
        DEBUG_PRINT("Case (x71) is True.\n");
        x96=Max(0, x26); FIX(x96);
        x97=0.5*x96; FIX(x97);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x100=Max(0, x23); FIX(x100);
        x101=0.5*x100; FIX(x101);
        /*@ assert (x21) != 0.;*/
        x121=1.0/x21; FIX(x121);
        x122=0.5*x121*x24; FIX(x122);
        x123=0.5*x121*x27; FIX(x123);
        x124=-x39; FIX(x124);
        x154=x20*x98; FIX(x154);
        x155=x61*x19*x98; FIX(x155);

        /* Assignment result[2, 5]=x61*x19*x123 + 1 - x101*(x155 + x39) - x61**2*x122 - x97*(x124 + x154) */
        x96=Max(0, x26); FIX(x96);
        x97=0.5*x96; FIX(x97);
        /*@ assert (x21) >= 0.;*/
        /*@ assert (x21) != 0.;*/
        x98=pow(x21, -3.0/2.0); FIX(x98);
        x100=Max(0, x23); FIX(x100);
        x101=0.5*x100; FIX(x101);
        /*@ assert (x21) != 0.;*/
        x121=1.0/x21; FIX(x121);
        x122=0.5*x121*x24; FIX(x122);
        x123=0.5*x121*x27; FIX(x123);
        x124=-x39; FIX(x124);
        x154=x20*x98; FIX(x154);
        x155=x61*x19*x98; FIX(x155);
        result[17] = x61*x19*x123 + 1 - x101*(x155 + x39) - x61*x61*x122 - x97*(x124 + x154);

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
