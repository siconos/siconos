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
#define FIX(x) do { assert(isfinite(x)); } while(0)
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
    // Heaviside
    double x3;
    double x23;
    double x24;
    double x25;
    double x26;
    int x27;
    double x20;
    double x21;
    double x22;
    /*@ assert (ut1*ut1 + ut2*ut2) >= 0.;*/
    x3=sqrt(ut1*ut1 + ut2*ut2); FIX(x3);
    /*@ assert (x3) >= 0.;*/
    x23=rt1*rt1; FIX(x23);
    /*@ assert (x23) >= 0.;*/
    x24=rt2*rt2; FIX(x24);
    /*@ assert (x24) >= 0.;*/
    x25=x23 + x24; FIX(x25);
    /*@ assert (x25) >= 0.;*/
    /*@ assert (x25) >= 0.;*/
    x26=sqrt(x25); FIX(x26);
    x27=x3 + x26 <= 0; FIX(x27);
    int x33;
    double x28;
    double x29;
    double x30;
    double x31;
    double x32;
    x33=x3 <= 0; FIX(x33);
    int x38;
    double x1;
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
    double x34;
    double x35;
    double x36;
    double x37;
    x38=x3 > 0; FIX(x38);
    int x64;
    x5=mu*ut1; FIX(x5);
    x6=-rt1 + x5; FIX(x6);
    x7=x6*x6; FIX(x7);
    x8=mu*ut2; FIX(x8);
    x9=-rt2 + x8; FIX(x9);
    x10=x9*x9; FIX(x10);
    x11=x10 + x7; FIX(x11);
    x64=x11 <= 0; FIX(x64);
    int x67;
    double x65;
    double x66;
    x67=x11 > 0; FIX(x67);
    int x82;
    double x14;
    double x18;
    double x46;
    double x49;
    double x69;
    double x79;
    double x80;
    double x81;
    /*@ assert (x11) >= 0.;*/
    x12=sqrt(x11); FIX(x12);
    x82=x12 > 0; FIX(x82);
    int x83;
    x83=x12 <= 0; FIX(x83);
    int x90;
    double x70;
    double x71;
    double x73;
    double x89;
    x90=x38 && x82; FIX(x90);
    int x91;
    x91=x38 && x83; FIX(x91);
    int x144;
    double x142;
    double x143;
    x144=x67 && x82; FIX(x144);
    int x145;
    x145=x67 && x83; FIX(x145);
    if (x27)
    {
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
    }
    else if (x33)
    {
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x32=0.5*x31; FIX(x32);
    }
    else if (x38)
    {
        x1=mu*rn; FIX(x1);
        x2=-un; FIX(x2);
        x4=-mu*x3 + x1 + x2; FIX(x4);
        x13=x12 + x4; FIX(x13);
        x17=-x12 + x4; FIX(x17);
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
    }
    else if (x64)
    {
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x32=0.5*x31; FIX(x32);
    }
    else if (x67)
    {
        x1=mu*rn; FIX(x1);
        x2=-un; FIX(x2);
        x4=-mu*x3 + x1 + x2; FIX(x4);
        x13=x12 + x4; FIX(x13);
        x17=-x12 + x4; FIX(x17);
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
        x65=-mu*x35; FIX(x65);
        x66=mu*x37; FIX(x66);
    }
    else if (x82)
    {
        x1=mu*rn; FIX(x1);
        x2=-un; FIX(x2);
        x4=-mu*x3 + x1 + x2; FIX(x4);
        x13=x12 + x4; FIX(x13);
        x14=Max(0, x13); FIX(x14);
        x17=-x12 + x4; FIX(x17);
        x18=Max(0, x17); FIX(x18);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x69=rt1 - x46; FIX(x69);
        x79=0.5*x14*x49; FIX(x79);
        x80=-rt1 + x5; FIX(x80);
        x81=0.5*x18*x49; FIX(x81);
    }
    else if (x90)
    {
        x1=mu*rn; FIX(x1);
        x2=-un; FIX(x2);
        x4=-mu*x3 + x1 + x2; FIX(x4);
        x13=x12 + x4; FIX(x13);
        x17=-x12 + x4; FIX(x17);
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x69=rt1 - x46; FIX(x69);
        x70=0.5*x34*x49; FIX(x70);
        x71=x69*x70; FIX(x71);
        x73=0.5*x36*x49; FIX(x73);
        x80=-rt1 + x5; FIX(x80);
        x89=x73*x80; FIX(x89);
    }
    else if (x144)
    {
        x1=mu*rn; FIX(x1);
        x2=-un; FIX(x2);
        x4=-mu*x3 + x1 + x2; FIX(x4);
        x13=x12 + x4; FIX(x13);
        x17=-x12 + x4; FIX(x17);
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x69=rt1 - x46; FIX(x69);
        x80=-rt1 + x5; FIX(x80);
        x142=0.5*mu*x34*x49; FIX(x142);
        x143=0.5*mu*x36*x49; FIX(x143);
    }
    /* Assignment result[0, 0]=x1 + x16 - x19 */
    double x15;
    double x16;
    double x19;x1=mu*rn; FIX(x1);
    x2=-un; FIX(x2);
    x4=-mu*x3 + x1 + x2; FIX(x4);
    x13=x12 + x4; FIX(x13);
    x14=Max(0, x13); FIX(x14);
    x15=0.5*x14; FIX(x15);
    x16=-x15; FIX(x16);
    x17=-x12 + x4; FIX(x17);
    x18=Max(0, x17); FIX(x18);
    x19=0.5*x18; FIX(x19);
    result[0] = x1 + x16 - x19;


    /* Assignment result[1, 0]=Piecewise((rt1 - x69*x79 - x80*x81, x82), (rt1, x83)) */

    if (x82)
    {
        DEBUG_PRINT("Case (x82) is True.\n");
        x46=mu*ut1; FIX(x46);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x69=rt1 - x46; FIX(x69);
        x79=0.5*x14*x49; FIX(x79);
        x80=-rt1 + x5; FIX(x80);
        x81=0.5*x18*x49; FIX(x81);

        /* Assignment result[1, 0]=rt1 - x69*x79 - x80*x81 */
        x46=mu*ut1; FIX(x46);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x69=rt1 - x46; FIX(x69);
        x79=0.5*x14*x49; FIX(x79);
        x80=-rt1 + x5; FIX(x80);
        x81=0.5*x18*x49; FIX(x81);
        result[1] = rt1 - x69*x79 - x80*x81;

    }
    else if (x83)
    {
        DEBUG_PRINT("Case (x83) is True.\n");

        /* Assignment result[1, 0]=rt1 */

        result[1] = rt1;

    }


    /* Assignment result[2, 0]=Piecewise((rt2 - x160*x81 - x75*x79, x82), (rt2 + x16 + x19, x83)) */
    double x56;
    double x75;
    double x160;
    if (x82)
    {
        DEBUG_PRINT("Case (x82) is True.\n");
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x56=mu*ut2; FIX(x56);
        x75=rt2 - x56; FIX(x75);
        x79=0.5*x14*x49; FIX(x79);
        x81=0.5*x18*x49; FIX(x81);
        x160=-rt2 + x8; FIX(x160);

        /* Assignment result[2, 0]=rt2 - x160*x81 - x75*x79 */
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x56=mu*ut2; FIX(x56);
        x75=rt2 - x56; FIX(x75);
        x79=0.5*x14*x49; FIX(x79);
        x81=0.5*x18*x49; FIX(x81);
        x160=-rt2 + x8; FIX(x160);
        result[2] = rt2 - x160*x81 - x75*x79;

    }
    else if (x83)
    {
        DEBUG_PRINT("Case (x83) is True.\n");

        /* Assignment result[2, 0]=rt2 + x16 + x19 */

        result[2] = rt2 + x16 + x19;

    }


    /* Assignment result[0, 1]=Piecewise((x22, x27), (x29 + x32, x33), (x35 + x37, x38)) */

    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);

        /* Assignment result[0, 1]=x22 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        result[3] = x22;

    }
    else if (x33)
    {
        DEBUG_PRINT("Case (x33) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x32=0.5*x31; FIX(x32);

        /* Assignment result[0, 1]=x29 + x32 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x32=0.5*x31; FIX(x32);
        result[3] = x29 + x32;

    }
    else if (x38)
    {
        DEBUG_PRINT("Case (x38) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);

        /* Assignment result[0, 1]=x35 + x37 */
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
        result[3] = x35 + x37;

    }


    /* Assignment result[1, 1]=Piecewise((0.0, x27), (-0.5*x40*x85*(x87 - 1.0*x88), x33), (x71 + x89, x90), (0, x91)) */
    double x40;
    double x63;
    double x84;
    double x85;
    double x86;
    double x87;
    double x88;
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");

        /* Assignment result[1, 1]=0.0 */

        result[4] = 0.0;
        /*@ assert (result[4]) >= 0.;*/
    }
    else if (x33)
    {
        DEBUG_PRINT("Case (x33) is True.\n");
        x30=-1.0*x26; FIX(x30);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x63=1.0*mu; FIX(x63);
        x84=Heaviside(x26); FIX(x84);
        x85=rt1*x84; FIX(x85);
        x86=-rn*x63 + un; FIX(x86);
        x87=Heaviside(x30 + x86); FIX(x87);
        x88=Heaviside(x26 + x86); FIX(x88);

        /* Assignment result[1, 1]=-0.5*x40*x85*(x87 - 1.0*x88) */
        x30=-1.0*x26; FIX(x30);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x63=1.0*mu; FIX(x63);
        x84=Heaviside(x26); FIX(x84);
        x85=rt1*x84; FIX(x85);
        x86=-rn*x63 + un; FIX(x86);
        x87=Heaviside(x30 + x86); FIX(x87);
        x88=Heaviside(x26 + x86); FIX(x88);
        result[4] = -0.5*x40*x85*(x87 - 1.0*x88);

    }
    else if (x90)
    {
        DEBUG_PRINT("Case (x90) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x69=rt1 - x46; FIX(x69);
        x70=0.5*x34*x49; FIX(x70);
        x71=x69*x70; FIX(x71);
        x73=0.5*x36*x49; FIX(x73);
        x80=-rt1 + x5; FIX(x80);
        x89=x73*x80; FIX(x89);

        /* Assignment result[1, 1]=x71 + x89 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x69=rt1 - x46; FIX(x69);
        x70=0.5*x34*x49; FIX(x70);
        x71=x69*x70; FIX(x71);
        x73=0.5*x36*x49; FIX(x73);
        x80=-rt1 + x5; FIX(x80);
        x89=x73*x80; FIX(x89);
        result[4] = x71 + x89;

    }
    else if (x91)
    {
        DEBUG_PRINT("Case (x91) is True.\n");

        /* Assignment result[1, 1]=0 */

        result[4] = 0;
        /*@ assert (result[4]) >= 0.;*/
    }


    /* Assignment result[2, 1]=Piecewise((0.0, x27), (x40*(x162*x87 - x162*x88 + x133*x161 - x135*x161 + x26*x29 - x26*x32), x33), (x163 + x76, x90), (x35 - x37, x91)) */
    double x76;
    double x133;
    double x135;
    double x161;
    double x162;
    double x163;
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");

        /* Assignment result[2, 1]=0.0 */

        result[5] = 0.0;
        /*@ assert (result[5]) >= 0.;*/
    }
    else if (x33)
    {
        DEBUG_PRINT("Case (x33) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x32=0.5*x31; FIX(x32);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x63=1.0*mu; FIX(x63);
        x84=Heaviside(x26); FIX(x84);
        x86=-rn*x63 + un; FIX(x86);
        x87=Heaviside(x30 + x86); FIX(x87);
        x88=Heaviside(x26 + x86); FIX(x88);
        x133=rt2*x88; FIX(x133);
        x135=rt2*x87; FIX(x135);
        x161=0.5*x84; FIX(x161);
        x162=0.5*x26*x84; FIX(x162);

        /* Assignment result[2, 1]=x40*(x162*x87 - x162*x88 + x133*x161 - x135*x161 + x26*x29 - x26*x32) */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x32=0.5*x31; FIX(x32);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x63=1.0*mu; FIX(x63);
        x84=Heaviside(x26); FIX(x84);
        x86=-rn*x63 + un; FIX(x86);
        x87=Heaviside(x30 + x86); FIX(x87);
        x88=Heaviside(x26 + x86); FIX(x88);
        x133=rt2*x88; FIX(x133);
        x135=rt2*x87; FIX(x135);
        x161=0.5*x84; FIX(x161);
        x162=0.5*x26*x84; FIX(x162);
        result[5] = x40*(x162*x87 - x162*x88 + x133*x161 - x135*x161 + x26*x29 - x26*x32);

    }
    else if (x90)
    {
        DEBUG_PRINT("Case (x90) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x56=mu*ut2; FIX(x56);
        x70=0.5*x34*x49; FIX(x70);
        x73=0.5*x36*x49; FIX(x73);
        x75=rt2 - x56; FIX(x75);
        x76=x70*x75; FIX(x76);
        x160=-rt2 + x8; FIX(x160);
        x163=x160*x73; FIX(x163);

        /* Assignment result[2, 1]=x163 + x76 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x56=mu*ut2; FIX(x56);
        x70=0.5*x34*x49; FIX(x70);
        x73=0.5*x36*x49; FIX(x73);
        x75=rt2 - x56; FIX(x75);
        x76=x70*x75; FIX(x76);
        x160=-rt2 + x8; FIX(x160);
        x163=x160*x73; FIX(x163);
        result[5] = x163 + x76;

    }
    else if (x91)
    {
        DEBUG_PRINT("Case (x91) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);

        /* Assignment result[2, 1]=x35 - x37 */
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
        result[5] = x35 - x37;

    }


    /* Assignment result[0, 2]=Piecewise((x39, x27), (x41*(rt1*x42 - rt1*x43 + x45), x33), (x53 - x55, x38)) */
    double x39;
    double x41;
    double x42;
    double x43;
    double x44;
    double x45;
    double x47;
    double x48;
    double x50;
    double x51;
    double x52;
    double x53;
    double x54;
    double x55;
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x39=mu*x22; FIX(x39);

        /* Assignment result[0, 2]=x39 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x39=mu*x22; FIX(x39);
        result[6] = x39;

    }
    else if (x33)
    {
        DEBUG_PRINT("Case (x33) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x41=0.25*mu*x40; FIX(x41);
        x42=2.0*x28; FIX(x42);
        x43=2.0*x31; FIX(x43);
        x44=1.4142135623730951454746218587388284504413604736328125*x26; FIX(x44);
        x45=x28*x44 + x31*x44; FIX(x45);

        /* Assignment result[0, 2]=x41*(rt1*x42 - rt1*x43 + x45) */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x41=0.25*mu*x40; FIX(x41);
        x42=2.0*x28; FIX(x42);
        x43=2.0*x31; FIX(x43);
        x44=1.4142135623730951454746218587388284504413604736328125*x26; FIX(x44);
        x45=x28*x44 + x31*x44; FIX(x45);
        result[6] = x41*(rt1*x42 - rt1*x43 + x45);

    }
    else if (x38)
    {
        DEBUG_PRINT("Case (x38) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        x48=-x46*x47; FIX(x48);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x51=x50*x6; FIX(x51);
        x52=x48 + x51; FIX(x52);
        x53=-x35*x52; FIX(x53);
        x54=x48 - x51; FIX(x54);
        x55=x37*x54; FIX(x55);

        /* Assignment result[0, 2]=x53 - x55 */
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        x48=-x46*x47; FIX(x48);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x51=x50*x6; FIX(x51);
        x52=x48 + x51; FIX(x52);
        x53=-x35*x52; FIX(x53);
        x54=x48 - x51; FIX(x54);
        x55=x37*x54; FIX(x55);
        result[6] = x53 - x55;

    }


    /* Assignment result[1, 2]=Piecewise((x95*x97 - x97*x98, x27), (-x100*x84*(-x102*x23 - x116*x24 + x118*x24 + 2.0*x119*x87 + 2.0*x119*x88 - x103*x104 + x103*x107 - x104*x106 - x105*x110 + x105*x115 + x106*x107 - x112*x113 + x112*x117), x33), (-x19*(-mu*x121 + x50) - x15*(-mu*x123 + x122) - x52*x71 - x54*x89, x90), (0, x91)) */
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
    double x107;
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
    double x118;
    double x119;
    double x120;
    double x121;
    double x122;
    double x123;
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x63=1.0*mu; FIX(x63);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) >= 0.;*/
        x94=sqrt(x93); FIX(x94);
        x95=mu + x94; FIX(x95);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x97=0.5*x22*x92*x96; FIX(x97);
        x98=-x63 + x94; FIX(x98);

        /* Assignment result[1, 2]=x95*x97 - x97*x98 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x63=1.0*mu; FIX(x63);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) >= 0.;*/
        x94=sqrt(x93); FIX(x94);
        x95=mu + x94; FIX(x95);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x97=0.5*x22*x92*x96; FIX(x97);
        x98=-x63 + x94; FIX(x98);
        result[7] = x95*x97 - x97*x98;

    }
    else if (x33)
    {
        DEBUG_PRINT("Case (x33) is True.\n");
        x20=mu*rn; FIX(x20);
        x30=-1.0*x26; FIX(x30);
        x63=1.0*mu; FIX(x63);
        x84=Heaviside(x26); FIX(x84);
        x86=-rn*x63 + un; FIX(x86);
        x87=Heaviside(x30 + x86); FIX(x87);
        x88=Heaviside(x26 + x86); FIX(x88);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        x100=0.25*mu*x99; FIX(x100);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x102=4.0*x101; FIX(x102);
        x103=pow(rt1, 5); FIX(x103);
        x104=1.4142135623730951454746218587388284504413604736328125*x88; FIX(x104);
        x105=pow(rt2, 4); FIX(x105);
        x106=rt1*x105; FIX(x106);
        x107=1.4142135623730951454746218587388284504413604736328125*x87; FIX(x107);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x110=2.0*x109; FIX(x110);
        x111=pow(rt1, 3); FIX(x111);
        x112=x111*x24; FIX(x112);
        x113=2.828427124746190290949243717477656900882720947265625*x88; FIX(x113);
        x114=Max(0, x108 - x26); FIX(x114);
        x115=2.0*x114; FIX(x115);
        x116=2.0*x109*x23; FIX(x116);
        x117=2.828427124746190290949243717477656900882720947265625*x87; FIX(x117);
        x118=2.0*x114*x23; FIX(x118);
        x119=x101*x23; FIX(x119);

        /* Assignment result[1, 2]=-x100*x84*(-x102*x23 - x116*x24 + x118*x24 + 2.0*x119*x87 + 2.0*x119*x88 - x103*x104 + x103*x107 - x104*x106 - x105*x110 + x105*x115 + x106*x107 - x112*x113 + x112*x117) */
        x20=mu*rn; FIX(x20);
        x30=-1.0*x26; FIX(x30);
        x63=1.0*mu; FIX(x63);
        x84=Heaviside(x26); FIX(x84);
        x86=-rn*x63 + un; FIX(x86);
        x87=Heaviside(x30 + x86); FIX(x87);
        x88=Heaviside(x26 + x86); FIX(x88);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        x100=0.25*mu*x99; FIX(x100);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x102=4.0*x101; FIX(x102);
        x103=pow(rt1, 5); FIX(x103);
        x104=1.4142135623730951454746218587388284504413604736328125*x88; FIX(x104);
        x105=pow(rt2, 4); FIX(x105);
        x106=rt1*x105; FIX(x106);
        x107=1.4142135623730951454746218587388284504413604736328125*x87; FIX(x107);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x110=2.0*x109; FIX(x110);
        x111=pow(rt1, 3); FIX(x111);
        x112=x111*x24; FIX(x112);
        x113=2.828427124746190290949243717477656900882720947265625*x88; FIX(x113);
        x114=Max(0, x108 - x26); FIX(x114);
        x115=2.0*x114; FIX(x115);
        x116=2.0*x109*x23; FIX(x116);
        x117=2.828427124746190290949243717477656900882720947265625*x87; FIX(x117);
        x118=2.0*x114*x23; FIX(x118);
        x119=x101*x23; FIX(x119);
        result[7] = -x100*x84*(-x102*x23 - x116*x24 + x118*x24 + 2.0*x119*x87 + 2.0*x119*x88 - x103*x104 + x103*x107 - x104*x106 - x105*x110 + x105*x115 + x106*x107 - x112*x113 + x112*x117);

    }
    else if (x90)
    {
        DEBUG_PRINT("Case (x90) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        x48=-x46*x47; FIX(x48);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x51=x50*x6; FIX(x51);
        x52=x48 + x51; FIX(x52);
        x54=x48 - x51; FIX(x54);
        x69=rt1 - x46; FIX(x69);
        x70=0.5*x34*x49; FIX(x70);
        x71=x69*x70; FIX(x71);
        x73=0.5*x36*x49; FIX(x73);
        x80=-rt1 + x5; FIX(x80);
        x89=x73*x80; FIX(x89);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        x121=x120*x7; FIX(x121);
        x122=-x50; FIX(x122);
        x123=x120*x6*x69; FIX(x123);

        /* Assignment result[1, 2]=-x19*(-mu*x121 + x50) - x15*(-mu*x123 + x122) - x52*x71 - x54*x89 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        x48=-x46*x47; FIX(x48);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x51=x50*x6; FIX(x51);
        x52=x48 + x51; FIX(x52);
        x54=x48 - x51; FIX(x54);
        x69=rt1 - x46; FIX(x69);
        x70=0.5*x34*x49; FIX(x70);
        x71=x69*x70; FIX(x71);
        x73=0.5*x36*x49; FIX(x73);
        x80=-rt1 + x5; FIX(x80);
        x89=x73*x80; FIX(x89);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        x121=x120*x7; FIX(x121);
        x122=-x50; FIX(x122);
        x123=x120*x6*x69; FIX(x123);
        result[7] = -x19*(-mu*x121 + x50) - x15*(-mu*x123 + x122) - x52*x71 - x54*x89;

    }
    else if (x91)
    {
        DEBUG_PRINT("Case (x91) is True.\n");

        /* Assignment result[1, 2]=0 */

        result[7] = 0;
        /*@ assert (result[7]) >= 0.;*/
    }


    /* Assignment result[2, 2]=Piecewise((-x164*x95 + x164*x98, x27), (-x100*(-x181*x87 - x181*x88 + x103*x169 - x103*x173 - x103*x177 - x103*x42 - x103*x43 - x104*x175 - x104*x176 + x106*x169 - x106*x173 - x106*x177 - x106*x42 - x106*x43 + x107*x175 + x107*x176 - x107*x183 + x112*x170 - x112*x171 - x112*x172 - x113*x178 + x117*x178 - x126*x85 + x129*x85 - x131*x85 + x134*x85 + x136*x85 + x146*x167 + x168 + x179*x180 - x179*x182), x33), (mu*x184 + x140 - x163*x54 - x52*x76, x90), (x53 + x55, x91)) */
    double x126;
    double x127;
    double x128;
    double x129;
    double x131;
    double x132;
    double x134;
    double x136;
    double x137;
    double x139;
    double x140;
    double x146;
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
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x63=1.0*mu; FIX(x63);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) >= 0.;*/
        x94=sqrt(x93); FIX(x94);
        x95=mu + x94; FIX(x95);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x98=-x63 + x94; FIX(x98);
        x164=0.5*mu*x22*x96; FIX(x164);

        /* Assignment result[2, 2]=-x164*x95 + x164*x98 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x63=1.0*mu; FIX(x63);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) >= 0.;*/
        x94=sqrt(x93); FIX(x94);
        x95=mu + x94; FIX(x95);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x98=-x63 + x94; FIX(x98);
        x164=0.5*mu*x22*x96; FIX(x164);
        result[8] = -x164*x95 + x164*x98;

    }
    else if (x33)
    {
        DEBUG_PRINT("Case (x33) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x42=2.0*x28; FIX(x42);
        x43=2.0*x31; FIX(x43);
        x63=1.0*mu; FIX(x63);
        x84=Heaviside(x26); FIX(x84);
        x85=rt1*x84; FIX(x85);
        x86=-rn*x63 + un; FIX(x86);
        x87=Heaviside(x30 + x86); FIX(x87);
        x88=Heaviside(x26 + x86); FIX(x88);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        x100=0.25*mu*x99; FIX(x100);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x102=4.0*x101; FIX(x102);
        x103=pow(rt1, 5); FIX(x103);
        x104=1.4142135623730951454746218587388284504413604736328125*x88; FIX(x104);
        x105=pow(rt2, 4); FIX(x105);
        x106=rt1*x105; FIX(x106);
        x107=1.4142135623730951454746218587388284504413604736328125*x87; FIX(x107);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x110=2.0*x109; FIX(x110);
        x111=pow(rt1, 3); FIX(x111);
        x112=x111*x24; FIX(x112);
        x113=2.828427124746190290949243717477656900882720947265625*x88; FIX(x113);
        x114=Max(0, x108 - x26); FIX(x114);
        x115=2.0*x114; FIX(x115);
        x117=2.828427124746190290949243717477656900882720947265625*x87; FIX(x117);
        x119=x101*x23; FIX(x119);
        x126=rt2*x102; FIX(x126);
        x127=pow(rt1, 4); FIX(x127);
        x128=pow(rt2, 3); FIX(x128);
        x129=x110*x128; FIX(x129);
        x131=x115*x128; FIX(x131);
        x132=2.0*x101; FIX(x132);
        x133=rt2*x88; FIX(x133);
        x134=x132*x133; FIX(x134);
        x135=rt2*x87; FIX(x135);
        x136=x132*x135; FIX(x136);
        x146=x101*x24; FIX(x146);
        x165=1.4142135623730951454746218587388284504413604736328125*x28; FIX(x165);
        x166=1.4142135623730951454746218587388284504413604736328125*x31; FIX(x166);
        x167=1.4142135623730951454746218587388284504413604736328125*x84*x88; FIX(x167);
        x168=-x101*x107*x23*x84 - x119*x165 + x119*x166 + x119*x167 - x146*x165 + x146*x166; FIX(x168);
        x169=4.0*x84; FIX(x169);
        x170=8.0*x84; FIX(x170);
        x171=4.0*x28; FIX(x171);
        x172=4.0*x31; FIX(x172);
        x173=2.0*x84*x88; FIX(x173);
        x174=pow(rt2, 5); FIX(x174);
        x175=x174*x84; FIX(x175);
        x176=rt2*x127*x84; FIX(x176);
        x177=2.0*x84*x87; FIX(x177);
        x178=x128*x23*x84; FIX(x178);
        x179=rt2*x111; FIX(x179);
        x180=2.0*x109*x84; FIX(x180);
        x181=4.0*x111*x24*x84; FIX(x181);
        x182=2.0*x114*x84; FIX(x182);
        x183=x101*x24*x84; FIX(x183);

        /* Assignment result[2, 2]=-x100*(-x181*x87 - x181*x88 + x103*x169 - x103*x173 - x103*x177 - x103*x42 - x103*x43 - x104*x175 - x104*x176 + x106*x169 - x106*x173 - x106*x177 - x106*x42 - x106*x43 + x107*x175 + x107*x176 - x107*x183 + x112*x170 - x112*x171 - x112*x172 - x113*x178 + x117*x178 - x126*x85 + x129*x85 - x131*x85 + x134*x85 + x136*x85 + x146*x167 + x168 + x179*x180 - x179*x182) */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x42=2.0*x28; FIX(x42);
        x43=2.0*x31; FIX(x43);
        x63=1.0*mu; FIX(x63);
        x84=Heaviside(x26); FIX(x84);
        x85=rt1*x84; FIX(x85);
        x86=-rn*x63 + un; FIX(x86);
        x87=Heaviside(x30 + x86); FIX(x87);
        x88=Heaviside(x26 + x86); FIX(x88);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        x100=0.25*mu*x99; FIX(x100);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x102=4.0*x101; FIX(x102);
        x103=pow(rt1, 5); FIX(x103);
        x104=1.4142135623730951454746218587388284504413604736328125*x88; FIX(x104);
        x105=pow(rt2, 4); FIX(x105);
        x106=rt1*x105; FIX(x106);
        x107=1.4142135623730951454746218587388284504413604736328125*x87; FIX(x107);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x110=2.0*x109; FIX(x110);
        x111=pow(rt1, 3); FIX(x111);
        x112=x111*x24; FIX(x112);
        x113=2.828427124746190290949243717477656900882720947265625*x88; FIX(x113);
        x114=Max(0, x108 - x26); FIX(x114);
        x115=2.0*x114; FIX(x115);
        x117=2.828427124746190290949243717477656900882720947265625*x87; FIX(x117);
        x119=x101*x23; FIX(x119);
        x126=rt2*x102; FIX(x126);
        x127=pow(rt1, 4); FIX(x127);
        x128=pow(rt2, 3); FIX(x128);
        x129=x110*x128; FIX(x129);
        x131=x115*x128; FIX(x131);
        x132=2.0*x101; FIX(x132);
        x133=rt2*x88; FIX(x133);
        x134=x132*x133; FIX(x134);
        x135=rt2*x87; FIX(x135);
        x136=x132*x135; FIX(x136);
        x146=x101*x24; FIX(x146);
        x165=1.4142135623730951454746218587388284504413604736328125*x28; FIX(x165);
        x166=1.4142135623730951454746218587388284504413604736328125*x31; FIX(x166);
        x167=1.4142135623730951454746218587388284504413604736328125*x84*x88; FIX(x167);
        x168=-x101*x107*x23*x84 - x119*x165 + x119*x166 + x119*x167 - x146*x165 + x146*x166; FIX(x168);
        x169=4.0*x84; FIX(x169);
        x170=8.0*x84; FIX(x170);
        x171=4.0*x28; FIX(x171);
        x172=4.0*x31; FIX(x172);
        x173=2.0*x84*x88; FIX(x173);
        x174=pow(rt2, 5); FIX(x174);
        x175=x174*x84; FIX(x175);
        x176=rt2*x127*x84; FIX(x176);
        x177=2.0*x84*x87; FIX(x177);
        x178=x128*x23*x84; FIX(x178);
        x179=rt2*x111; FIX(x179);
        x180=2.0*x109*x84; FIX(x180);
        x181=4.0*x111*x24*x84; FIX(x181);
        x182=2.0*x114*x84; FIX(x182);
        x183=x101*x24*x84; FIX(x183);
        result[8] = -x100*(-x181*x87 - x181*x88 + x103*x169 - x103*x173 - x103*x177 - x103*x42 - x103*x43 - x104*x175 - x104*x176 + x106*x169 - x106*x173 - x106*x177 - x106*x42 - x106*x43 + x107*x175 + x107*x176 - x107*x183 + x112*x170 - x112*x171 - x112*x172 - x113*x178 + x117*x178 - x126*x85 + x129*x85 - x131*x85 + x134*x85 + x136*x85 + x146*x167 + x168 + x179*x180 - x179*x182);

    }
    else if (x90)
    {
        DEBUG_PRINT("Case (x90) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        x48=-x46*x47; FIX(x48);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x51=x50*x6; FIX(x51);
        x52=x48 + x51; FIX(x52);
        x54=x48 - x51; FIX(x54);
        x56=mu*ut2; FIX(x56);
        x70=0.5*x34*x49; FIX(x70);
        x73=0.5*x36*x49; FIX(x73);
        x75=rt2 - x56; FIX(x75);
        x76=x70*x75; FIX(x76);
        x80=-rt1 + x5; FIX(x80);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        x137=0.5*x120*x14; FIX(x137);
        x139=x120*x19*x80*x9; FIX(x139);
        x140=mu*x139; FIX(x140);
        x160=-rt2 + x8; FIX(x160);
        x163=x160*x73; FIX(x163);
        x184=x137*x6*x75; FIX(x184);

        /* Assignment result[2, 2]=mu*x184 + x140 - x163*x54 - x52*x76 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        x48=-x46*x47; FIX(x48);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x51=x50*x6; FIX(x51);
        x52=x48 + x51; FIX(x52);
        x54=x48 - x51; FIX(x54);
        x56=mu*ut2; FIX(x56);
        x70=0.5*x34*x49; FIX(x70);
        x73=0.5*x36*x49; FIX(x73);
        x75=rt2 - x56; FIX(x75);
        x76=x70*x75; FIX(x76);
        x80=-rt1 + x5; FIX(x80);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        x137=0.5*x120*x14; FIX(x137);
        x139=x120*x19*x80*x9; FIX(x139);
        x140=mu*x139; FIX(x140);
        x160=-rt2 + x8; FIX(x160);
        x163=x160*x73; FIX(x163);
        x184=x137*x6*x75; FIX(x184);
        result[8] = mu*x184 + x140 - x163*x54 - x52*x76;

    }
    else if (x91)
    {
        DEBUG_PRINT("Case (x91) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        x48=-x46*x47; FIX(x48);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x51=x50*x6; FIX(x51);
        x52=x48 + x51; FIX(x52);
        x53=-x35*x52; FIX(x53);
        x54=x48 - x51; FIX(x54);
        x55=x37*x54; FIX(x55);

        /* Assignment result[2, 2]=x53 + x55 */
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        x48=-x46*x47; FIX(x48);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x51=x50*x6; FIX(x51);
        x52=x48 + x51; FIX(x52);
        x53=-x35*x52; FIX(x53);
        x54=x48 - x51; FIX(x54);
        x55=x37*x54; FIX(x55);
        result[8] = x53 + x55;

    }


    /* Assignment result[0, 3]=Piecewise((0.0, x27), (x41*(rt2*x42 - rt2*x43 + x45), x33), (x60 - x62, x38)) */
    double x57;
    double x58;
    double x59;
    double x60;
    double x61;
    double x62;
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");

        /* Assignment result[0, 3]=0.0 */

        result[9] = 0.0;
        /*@ assert (result[9]) >= 0.;*/
    }
    else if (x33)
    {
        DEBUG_PRINT("Case (x33) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x41=0.25*mu*x40; FIX(x41);
        x42=2.0*x28; FIX(x42);
        x43=2.0*x31; FIX(x43);
        x44=1.4142135623730951454746218587388284504413604736328125*x26; FIX(x44);
        x45=x28*x44 + x31*x44; FIX(x45);

        /* Assignment result[0, 3]=x41*(rt2*x42 - rt2*x43 + x45) */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x41=0.25*mu*x40; FIX(x41);
        x42=2.0*x28; FIX(x42);
        x43=2.0*x31; FIX(x43);
        x44=1.4142135623730951454746218587388284504413604736328125*x26; FIX(x44);
        x45=x28*x44 + x31*x44; FIX(x45);
        result[9] = x41*(rt2*x42 - rt2*x43 + x45);

    }
    else if (x38)
    {
        DEBUG_PRINT("Case (x38) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x56=mu*ut2; FIX(x56);
        x57=-x47*x56; FIX(x57);
        x58=x50*x9; FIX(x58);
        x59=x57 + x58; FIX(x59);
        x60=-x35*x59; FIX(x60);
        x61=x57 - x58; FIX(x61);
        x62=x37*x61; FIX(x62);

        /* Assignment result[0, 3]=x60 - x62 */
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x56=mu*ut2; FIX(x56);
        x57=-x47*x56; FIX(x57);
        x58=x50*x9; FIX(x58);
        x59=x57 + x58; FIX(x59);
        x60=-x35*x59; FIX(x60);
        x61=x57 - x58; FIX(x61);
        x62=x37*x61; FIX(x62);
        result[9] = x60 - x62;

    }


    /* Assignment result[1, 3]=Piecewise((x125, x27), (-rt1*x100*x84*(rt2*x116 - rt2*x118 - x104*x105 - x104*x127 + x105*x107 + x107*x127 - x113*x130 + x117*x130 - x126 + x129 - x131 + x134 + x136), x33), (mu*x138 + x140 - x59*x71 - x61*x89, x90), (0, x91)) */
    double x124;
    double x125;
    double x130;
    double x138;
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x124=1.0*x22*x96; FIX(x124);
        x125=-x124*x92; FIX(x125);

        /* Assignment result[1, 3]=x125 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x124=1.0*x22*x96; FIX(x124);
        x125=-x124*x92; FIX(x125);
        result[10] = x125;

    }
    else if (x33)
    {
        DEBUG_PRINT("Case (x33) is True.\n");
        x20=mu*rn; FIX(x20);
        x30=-1.0*x26; FIX(x30);
        x63=1.0*mu; FIX(x63);
        x84=Heaviside(x26); FIX(x84);
        x86=-rn*x63 + un; FIX(x86);
        x87=Heaviside(x30 + x86); FIX(x87);
        x88=Heaviside(x26 + x86); FIX(x88);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        x100=0.25*mu*x99; FIX(x100);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x102=4.0*x101; FIX(x102);
        x104=1.4142135623730951454746218587388284504413604736328125*x88; FIX(x104);
        x105=pow(rt2, 4); FIX(x105);
        x107=1.4142135623730951454746218587388284504413604736328125*x87; FIX(x107);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x110=2.0*x109; FIX(x110);
        x113=2.828427124746190290949243717477656900882720947265625*x88; FIX(x113);
        x114=Max(0, x108 - x26); FIX(x114);
        x115=2.0*x114; FIX(x115);
        x116=2.0*x109*x23; FIX(x116);
        x117=2.828427124746190290949243717477656900882720947265625*x87; FIX(x117);
        x118=2.0*x114*x23; FIX(x118);
        x126=rt2*x102; FIX(x126);
        x127=pow(rt1, 4); FIX(x127);
        x128=pow(rt2, 3); FIX(x128);
        x129=x110*x128; FIX(x129);
        x130=x23*x24; FIX(x130);
        x131=x115*x128; FIX(x131);
        x132=2.0*x101; FIX(x132);
        x133=rt2*x88; FIX(x133);
        x134=x132*x133; FIX(x134);
        x135=rt2*x87; FIX(x135);
        x136=x132*x135; FIX(x136);

        /* Assignment result[1, 3]=-rt1*x100*x84*(rt2*x116 - rt2*x118 - x104*x105 - x104*x127 + x105*x107 + x107*x127 - x113*x130 + x117*x130 - x126 + x129 - x131 + x134 + x136) */
        x20=mu*rn; FIX(x20);
        x30=-1.0*x26; FIX(x30);
        x63=1.0*mu; FIX(x63);
        x84=Heaviside(x26); FIX(x84);
        x86=-rn*x63 + un; FIX(x86);
        x87=Heaviside(x30 + x86); FIX(x87);
        x88=Heaviside(x26 + x86); FIX(x88);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        x100=0.25*mu*x99; FIX(x100);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x102=4.0*x101; FIX(x102);
        x104=1.4142135623730951454746218587388284504413604736328125*x88; FIX(x104);
        x105=pow(rt2, 4); FIX(x105);
        x107=1.4142135623730951454746218587388284504413604736328125*x87; FIX(x107);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x110=2.0*x109; FIX(x110);
        x113=2.828427124746190290949243717477656900882720947265625*x88; FIX(x113);
        x114=Max(0, x108 - x26); FIX(x114);
        x115=2.0*x114; FIX(x115);
        x116=2.0*x109*x23; FIX(x116);
        x117=2.828427124746190290949243717477656900882720947265625*x87; FIX(x117);
        x118=2.0*x114*x23; FIX(x118);
        x126=rt2*x102; FIX(x126);
        x127=pow(rt1, 4); FIX(x127);
        x128=pow(rt2, 3); FIX(x128);
        x129=x110*x128; FIX(x129);
        x130=x23*x24; FIX(x130);
        x131=x115*x128; FIX(x131);
        x132=2.0*x101; FIX(x132);
        x133=rt2*x88; FIX(x133);
        x134=x132*x133; FIX(x134);
        x135=rt2*x87; FIX(x135);
        x136=x132*x135; FIX(x136);
        result[10] = -rt1*x100*x84*(rt2*x116 - rt2*x118 - x104*x105 - x104*x127 + x105*x107 + x107*x127 - x113*x130 + x117*x130 - x126 + x129 - x131 + x134 + x136);

    }
    else if (x90)
    {
        DEBUG_PRINT("Case (x90) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x56=mu*ut2; FIX(x56);
        x57=-x47*x56; FIX(x57);
        x58=x50*x9; FIX(x58);
        x59=x57 + x58; FIX(x59);
        x61=x57 - x58; FIX(x61);
        x69=rt1 - x46; FIX(x69);
        x70=0.5*x34*x49; FIX(x70);
        x71=x69*x70; FIX(x71);
        x73=0.5*x36*x49; FIX(x73);
        x80=-rt1 + x5; FIX(x80);
        x89=x73*x80; FIX(x89);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        x137=0.5*x120*x14; FIX(x137);
        x138=x137*x69*x9; FIX(x138);
        x139=x120*x19*x80*x9; FIX(x139);
        x140=mu*x139; FIX(x140);

        /* Assignment result[1, 3]=mu*x138 + x140 - x59*x71 - x61*x89 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x56=mu*ut2; FIX(x56);
        x57=-x47*x56; FIX(x57);
        x58=x50*x9; FIX(x58);
        x59=x57 + x58; FIX(x59);
        x61=x57 - x58; FIX(x61);
        x69=rt1 - x46; FIX(x69);
        x70=0.5*x34*x49; FIX(x70);
        x71=x69*x70; FIX(x71);
        x73=0.5*x36*x49; FIX(x73);
        x80=-rt1 + x5; FIX(x80);
        x89=x73*x80; FIX(x89);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        x137=0.5*x120*x14; FIX(x137);
        x138=x137*x69*x9; FIX(x138);
        x139=x120*x19*x80*x9; FIX(x139);
        x140=mu*x139; FIX(x140);
        result[10] = mu*x138 + x140 - x59*x71 - x61*x89;

    }
    else if (x91)
    {
        DEBUG_PRINT("Case (x91) is True.\n");

        /* Assignment result[1, 3]=0 */

        result[10] = 0;
        /*@ assert (result[10]) >= 0.;*/
    }


    /* Assignment result[2, 3]=Piecewise((x155, x27), (-x100*(-1.17157287525381*x178*x87 - 6.82842712474619*x178*x88 - x102*x187 - x116*x187 + x118*x187 - x127*x180 + x127*x182 + x168 + x169*x174 + x169*x185 + x170*x186 - x171*x186 - x172*x186 - x174*x42 - x174*x43 - x175*x188 - x175*x189 - x176*x188 - x176*x189 + x183*x188 + x183*x189 - x185*x42 - x185*x43), x33), (-x19*(-mu*x190 + x50) - x15*(-mu*x191 + x122) - x163*x61 - x59*x76, x90), (x60 + x62, x91)) */
    double x155;
    double x185;
    double x186;
    double x187;
    double x188;
    double x189;
    double x190;
    double x191;
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x39=mu*x22; FIX(x39);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x155=x39*x96; FIX(x155);

        /* Assignment result[2, 3]=x155 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x39=mu*x22; FIX(x39);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x155=x39*x96; FIX(x155);
        result[11] = x155;

    }
    else if (x33)
    {
        DEBUG_PRINT("Case (x33) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x42=2.0*x28; FIX(x42);
        x43=2.0*x31; FIX(x43);
        x63=1.0*mu; FIX(x63);
        x84=Heaviside(x26); FIX(x84);
        x86=-rn*x63 + un; FIX(x86);
        x87=Heaviside(x30 + x86); FIX(x87);
        x88=Heaviside(x26 + x86); FIX(x88);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        x100=0.25*mu*x99; FIX(x100);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x102=4.0*x101; FIX(x102);
        x107=1.4142135623730951454746218587388284504413604736328125*x87; FIX(x107);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x114=Max(0, x108 - x26); FIX(x114);
        x116=2.0*x109*x23; FIX(x116);
        x118=2.0*x114*x23; FIX(x118);
        x119=x101*x23; FIX(x119);
        x127=pow(rt1, 4); FIX(x127);
        x128=pow(rt2, 3); FIX(x128);
        x146=x101*x24; FIX(x146);
        x165=1.4142135623730951454746218587388284504413604736328125*x28; FIX(x165);
        x166=1.4142135623730951454746218587388284504413604736328125*x31; FIX(x166);
        x167=1.4142135623730951454746218587388284504413604736328125*x84*x88; FIX(x167);
        x168=-x101*x107*x23*x84 - x119*x165 + x119*x166 + x119*x167 - x146*x165 + x146*x166; FIX(x168);
        x169=4.0*x84; FIX(x169);
        x170=8.0*x84; FIX(x170);
        x171=4.0*x28; FIX(x171);
        x172=4.0*x31; FIX(x172);
        x174=pow(rt2, 5); FIX(x174);
        x175=x174*x84; FIX(x175);
        x176=rt2*x127*x84; FIX(x176);
        x178=x128*x23*x84; FIX(x178);
        x180=2.0*x109*x84; FIX(x180);
        x182=2.0*x114*x84; FIX(x182);
        x183=x101*x24*x84; FIX(x183);
        x185=rt2*x127; FIX(x185);
        x186=x128*x23; FIX(x186);
        x187=x24*x84; FIX(x187);
        x188=3.41421356237309492343001693370752036571502685546875*x88; FIX(x188);
        x189=0.58578643762690496554768060377682559192180633544921875*x87; FIX(x189);

        /* Assignment result[2, 3]=-x100*(-1.17157287525381*x178*x87 - 6.82842712474619*x178*x88 - x102*x187 - x116*x187 + x118*x187 - x127*x180 + x127*x182 + x168 + x169*x174 + x169*x185 + x170*x186 - x171*x186 - x172*x186 - x174*x42 - x174*x43 - x175*x188 - x175*x189 - x176*x188 - x176*x189 + x183*x188 + x183*x189 - x185*x42 - x185*x43) */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x42=2.0*x28; FIX(x42);
        x43=2.0*x31; FIX(x43);
        x63=1.0*mu; FIX(x63);
        x84=Heaviside(x26); FIX(x84);
        x86=-rn*x63 + un; FIX(x86);
        x87=Heaviside(x30 + x86); FIX(x87);
        x88=Heaviside(x26 + x86); FIX(x88);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        x100=0.25*mu*x99; FIX(x100);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x102=4.0*x101; FIX(x102);
        x107=1.4142135623730951454746218587388284504413604736328125*x87; FIX(x107);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x114=Max(0, x108 - x26); FIX(x114);
        x116=2.0*x109*x23; FIX(x116);
        x118=2.0*x114*x23; FIX(x118);
        x119=x101*x23; FIX(x119);
        x127=pow(rt1, 4); FIX(x127);
        x128=pow(rt2, 3); FIX(x128);
        x146=x101*x24; FIX(x146);
        x165=1.4142135623730951454746218587388284504413604736328125*x28; FIX(x165);
        x166=1.4142135623730951454746218587388284504413604736328125*x31; FIX(x166);
        x167=1.4142135623730951454746218587388284504413604736328125*x84*x88; FIX(x167);
        x168=-x101*x107*x23*x84 - x119*x165 + x119*x166 + x119*x167 - x146*x165 + x146*x166; FIX(x168);
        x169=4.0*x84; FIX(x169);
        x170=8.0*x84; FIX(x170);
        x171=4.0*x28; FIX(x171);
        x172=4.0*x31; FIX(x172);
        x174=pow(rt2, 5); FIX(x174);
        x175=x174*x84; FIX(x175);
        x176=rt2*x127*x84; FIX(x176);
        x178=x128*x23*x84; FIX(x178);
        x180=2.0*x109*x84; FIX(x180);
        x182=2.0*x114*x84; FIX(x182);
        x183=x101*x24*x84; FIX(x183);
        x185=rt2*x127; FIX(x185);
        x186=x128*x23; FIX(x186);
        x187=x24*x84; FIX(x187);
        x188=3.41421356237309492343001693370752036571502685546875*x88; FIX(x188);
        x189=0.58578643762690496554768060377682559192180633544921875*x87; FIX(x189);
        result[11] = -x100*(-1.1715728752538099310953612075536511838436126708984375*x178*x87 - 6.8284271247461898468600338674150407314300537109375*x178*x88 - x102*x187 - x116*x187 + x118*x187 - x127*x180 + x127*x182 + x168 + x169*x174 + x169*x185 + x170*x186 - x171*x186 - x172*x186 - x174*x42 - x174*x43 - x175*x188 - x175*x189 - x176*x188 - x176*x189 + x183*x188 + x183*x189 - x185*x42 - x185*x43);

    }
    else if (x90)
    {
        DEBUG_PRINT("Case (x90) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x56=mu*ut2; FIX(x56);
        x57=-x47*x56; FIX(x57);
        x58=x50*x9; FIX(x58);
        x59=x57 + x58; FIX(x59);
        x61=x57 - x58; FIX(x61);
        x70=0.5*x34*x49; FIX(x70);
        x73=0.5*x36*x49; FIX(x73);
        x75=rt2 - x56; FIX(x75);
        x76=x70*x75; FIX(x76);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        x122=-x50; FIX(x122);
        x160=-rt2 + x8; FIX(x160);
        x163=x160*x73; FIX(x163);
        x190=x10*x120; FIX(x190);
        x191=x120*x75*x9; FIX(x191);

        /* Assignment result[2, 3]=-x19*(-mu*x190 + x50) - x15*(-mu*x191 + x122) - x163*x61 - x59*x76 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x56=mu*ut2; FIX(x56);
        x57=-x47*x56; FIX(x57);
        x58=x50*x9; FIX(x58);
        x59=x57 + x58; FIX(x59);
        x61=x57 - x58; FIX(x61);
        x70=0.5*x34*x49; FIX(x70);
        x73=0.5*x36*x49; FIX(x73);
        x75=rt2 - x56; FIX(x75);
        x76=x70*x75; FIX(x76);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        x122=-x50; FIX(x122);
        x160=-rt2 + x8; FIX(x160);
        x163=x160*x73; FIX(x163);
        x190=x10*x120; FIX(x190);
        x191=x120*x75*x9; FIX(x191);
        result[11] = -x19*(-mu*x190 + x50) - x15*(-mu*x191 + x122) - x163*x61 - x59*x76;

    }
    else if (x91)
    {
        DEBUG_PRINT("Case (x91) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x56=mu*ut2; FIX(x56);
        x57=-x47*x56; FIX(x57);
        x58=x50*x9; FIX(x58);
        x59=x57 + x58; FIX(x59);
        x60=-x35*x59; FIX(x60);
        x61=x57 - x58; FIX(x61);
        x62=x37*x61; FIX(x62);

        /* Assignment result[2, 3]=x60 + x62 */
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
        /*@ assert (x3) != 0.;*/
        x47=1.0/x3; FIX(x47);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x50=mu*x49; FIX(x50);
        x56=mu*ut2; FIX(x56);
        x57=-x47*x56; FIX(x57);
        x58=x50*x9; FIX(x58);
        x59=x57 + x58; FIX(x59);
        x60=-x35*x59; FIX(x60);
        x61=x57 - x58; FIX(x61);
        x62=x37*x61; FIX(x62);
        result[11] = x60 + x62;

    }


    /* Assignment result[0, 4]=Piecewise((mu - x22*x63, x27), (mu - mu*x29 - mu*x32, x64), (mu + x65 - x66, x67)) */

    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x63=1.0*mu; FIX(x63);

        /* Assignment result[0, 4]=mu - x22*x63 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x63=1.0*mu; FIX(x63);
        result[12] = mu - x22*x63;

    }
    else if (x64)
    {
        DEBUG_PRINT("Case (x64) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x32=0.5*x31; FIX(x32);

        /* Assignment result[0, 4]=mu - mu*x29 - mu*x32 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x32=0.5*x31; FIX(x32);
        result[12] = mu - mu*x29 - mu*x32;

    }
    else if (x67)
    {
        DEBUG_PRINT("Case (x67) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
        x65=-mu*x35; FIX(x65);
        x66=mu*x37; FIX(x66);

        /* Assignment result[0, 4]=mu + x65 - x66 */
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
        x65=-mu*x35; FIX(x65);
        x66=mu*x37; FIX(x66);
        result[12] = mu + x65 - x66;

    }


    /* Assignment result[1, 4]=Piecewise((0.0, x27), (-rt1*x141, x64), (-x142*x69 - x143*x80, x144), (0, x145)) */
    double x141;
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");

        /* Assignment result[1, 4]=0.0 */

        result[13] = 0.0;
        /*@ assert (result[13]) >= 0.;*/
    }
    else if (x64)
    {
        DEBUG_PRINT("Case (x64) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x141=0.5*mu*x40*(x28 - 1.0*x31); FIX(x141);

        /* Assignment result[1, 4]=-rt1*x141 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x141=0.5*mu*x40*(x28 - 1.0*x31); FIX(x141);
        result[13] = -rt1*x141;

    }
    else if (x144)
    {
        DEBUG_PRINT("Case (x144) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x69=rt1 - x46; FIX(x69);
        x80=-rt1 + x5; FIX(x80);
        x142=0.5*mu*x34*x49; FIX(x142);
        x143=0.5*mu*x36*x49; FIX(x143);

        /* Assignment result[1, 4]=-x142*x69 - x143*x80 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x69=rt1 - x46; FIX(x69);
        x80=-rt1 + x5; FIX(x80);
        x142=0.5*mu*x34*x49; FIX(x142);
        x143=0.5*mu*x36*x49; FIX(x143);
        result[13] = -x142*x69 - x143*x80;

    }
    else if (x145)
    {
        DEBUG_PRINT("Case (x145) is True.\n");

        /* Assignment result[1, 4]=0 */

        result[13] = 0;
        /*@ assert (result[13]) >= 0.;*/
    }


    /* Assignment result[2, 4]=Piecewise((0.0, x27), (-rt2*x141, x64), (-x142*x75 - x143*x160, x144), (x65 + x66, x145)) */

    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");

        /* Assignment result[2, 4]=0.0 */

        result[14] = 0.0;
        /*@ assert (result[14]) >= 0.;*/
    }
    else if (x64)
    {
        DEBUG_PRINT("Case (x64) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x141=0.5*mu*x40*(x28 - 1.0*x31); FIX(x141);

        /* Assignment result[2, 4]=-rt2*x141 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x141=0.5*mu*x40*(x28 - 1.0*x31); FIX(x141);
        result[14] = -rt2*x141;

    }
    else if (x144)
    {
        DEBUG_PRINT("Case (x144) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x56=mu*ut2; FIX(x56);
        x75=rt2 - x56; FIX(x75);
        x142=0.5*mu*x34*x49; FIX(x142);
        x143=0.5*mu*x36*x49; FIX(x143);
        x160=-rt2 + x8; FIX(x160);

        /* Assignment result[2, 4]=-x142*x75 - x143*x160 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x56=mu*ut2; FIX(x56);
        x75=rt2 - x56; FIX(x75);
        x142=0.5*mu*x34*x49; FIX(x142);
        x143=0.5*mu*x36*x49; FIX(x143);
        x160=-rt2 + x8; FIX(x160);
        result[14] = -x142*x75 - x143*x160;

    }
    else if (x145)
    {
        DEBUG_PRINT("Case (x145) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
        x65=-mu*x35; FIX(x65);
        x66=mu*x37; FIX(x66);

        /* Assignment result[2, 4]=x65 + x66 */
        x34=Heaviside(x13); FIX(x34);
        x35=0.5*x34; FIX(x35);
        x36=Heaviside(x17); FIX(x36);
        x37=0.5*x36; FIX(x37);
        x65=-mu*x35; FIX(x65);
        x66=mu*x37; FIX(x66);
        result[14] = x65 + x66;

    }


    /* Assignment result[0, 5]=Piecewise((0.0, x27), (rt1*x68, x64), (x72 + x74, x67)) */
    double x68;
    double x72;
    double x74;
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");

        /* Assignment result[0, 5]=0.0 */

        result[15] = 0.0;
        /*@ assert (result[15]) >= 0.;*/
    }
    else if (x64)
    {
        DEBUG_PRINT("Case (x64) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x68=0.5*x40*(-1.0*x28 + x31); FIX(x68);

        /* Assignment result[0, 5]=rt1*x68 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x68=0.5*x40*(-1.0*x28 + x31); FIX(x68);
        result[15] = rt1*x68;

    }
    else if (x67)
    {
        DEBUG_PRINT("Case (x67) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x69=rt1 - x46; FIX(x69);
        x70=0.5*x34*x49; FIX(x70);
        x71=x69*x70; FIX(x71);
        x72=-x71; FIX(x72);
        x73=0.5*x36*x49; FIX(x73);
        x74=x69*x73; FIX(x74);

        /* Assignment result[0, 5]=x72 + x74 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x69=rt1 - x46; FIX(x69);
        x70=0.5*x34*x49; FIX(x70);
        x71=x69*x70; FIX(x71);
        x72=-x71; FIX(x72);
        x73=0.5*x36*x49; FIX(x73);
        x74=x69*x73; FIX(x74);
        result[15] = x72 + x74;

    }


    /* Assignment result[1, 5]=Piecewise((1.0 + x125, x27), (x99*(-x105*x147 + x105*x150 - x119*x29 - x119*x32 + x149), x64), (-x152*x69**2 - x19*(x121 + x154) - x15*(x123 + x49) + 1 + x153*x80, x144), (1, x145)) */
    double x147;
    double x148;
    double x149;
    double x150;
    double x151;
    double x152;
    double x153;
    double x154;
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x124=1.0*x22*x96; FIX(x124);
        x125=-x124*x92; FIX(x125);

        /* Assignment result[1, 5]=1.0 + x125 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x124=1.0*x22*x96; FIX(x124);
        x125=-x124*x92; FIX(x125);
        result[16] = 1.0 + x125;

    }
    else if (x64)
    {
        DEBUG_PRINT("Case (x64) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x32=0.5*x31; FIX(x32);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x105=pow(rt2, 4); FIX(x105);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x114=Max(0, x108 - x26); FIX(x114);
        x119=x101*x23; FIX(x119);
        x130=x23*x24; FIX(x130);
        x146=x101*x24; FIX(x146);
        x147=0.5*x109; FIX(x147);
        x148=x114*x23; FIX(x148);
        x149=0.5*x148*x24 + x119 - x130*x147 + x146; FIX(x149);
        x150=0.5*x114; FIX(x150);

        /* Assignment result[1, 5]=x99*(-x105*x147 + x105*x150 - x119*x29 - x119*x32 + x149) */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x32=0.5*x31; FIX(x32);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x105=pow(rt2, 4); FIX(x105);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x114=Max(0, x108 - x26); FIX(x114);
        x119=x101*x23; FIX(x119);
        x130=x23*x24; FIX(x130);
        x146=x101*x24; FIX(x146);
        x147=0.5*x109; FIX(x147);
        x148=x114*x23; FIX(x148);
        x149=0.5*x148*x24 + x119 - x130*x147 + x146; FIX(x149);
        x150=0.5*x114; FIX(x150);
        result[16] = x99*(-x105*x147 + x105*x150 - x119*x29 - x119*x32 + x149);

    }
    else if (x144)
    {
        DEBUG_PRINT("Case (x144) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x69=rt1 - x46; FIX(x69);
        x80=-rt1 + x5; FIX(x80);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        x121=x120*x7; FIX(x121);
        x123=x120*x6*x69; FIX(x123);
        /*@ assert (x11) != 0.;*/
        x151=1.0/x11; FIX(x151);
        x152=0.5*x151*x34; FIX(x152);
        x153=0.5*x151*x36*x69; FIX(x153);
        x154=-x49; FIX(x154);

        /* Assignment result[1, 5]=-x152*x69**2 - x19*(x121 + x154) - x15*(x123 + x49) + 1 + x153*x80 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x69=rt1 - x46; FIX(x69);
        x80=-rt1 + x5; FIX(x80);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        x121=x120*x7; FIX(x121);
        x123=x120*x6*x69; FIX(x123);
        /*@ assert (x11) != 0.;*/
        x151=1.0/x11; FIX(x151);
        x152=0.5*x151*x34; FIX(x152);
        x153=0.5*x151*x36*x69; FIX(x153);
        x154=-x49; FIX(x154);
        result[16] = -x152*x69*x69 - x19*(x121 + x154) - x15*(x123 + x49) + 1 + x153*x80;

    }
    else if (x145)
    {
        DEBUG_PRINT("Case (x145) is True.\n");

        /* Assignment result[1, 5]=1 */

        result[16] = 1;
        /*@ assert (result[16]) >= 0.;*/
        /*@ assert (result[16]) != 0.;*/
    }


    /* Assignment result[2, 5]=Piecewise((x155, x27), (x157, x64), (x153*x160 + x158 - x184, x144), (x72 - x74, x145)) */
    double x156;
    double x157;
    double x158;
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x39=mu*x22; FIX(x39);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x155=x39*x96; FIX(x155);

        /* Assignment result[2, 5]=x155 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x39=mu*x22; FIX(x39);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x155=x39*x96; FIX(x155);
        result[17] = x155;

    }
    else if (x64)
    {
        DEBUG_PRINT("Case (x64) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x114=Max(0, x108 - x26); FIX(x114);
        x148=x114*x23; FIX(x148);
        x156=1.0*x109; FIX(x156);
        x157=-0.5*rt1*rt2*x99*(x114*x24 - x156*x23 - x156*x24 + x101*x28 + x101*x31 + x148); FIX(x157);

        /* Assignment result[2, 5]=x157 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x114=Max(0, x108 - x26); FIX(x114);
        x148=x114*x23; FIX(x148);
        x156=1.0*x109; FIX(x156);
        x157=-0.5*rt1*rt2*x99*(x114*x24 - x156*x23 - x156*x24 + x101*x28 + x101*x31 + x148); FIX(x157);
        result[17] = x157;

    }
    else if (x144)
    {
        DEBUG_PRINT("Case (x144) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        x56=mu*ut2; FIX(x56);
        x69=rt1 - x46; FIX(x69);
        x75=rt2 - x56; FIX(x75);
        x80=-rt1 + x5; FIX(x80);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        x137=0.5*x120*x14; FIX(x137);
        x139=x120*x19*x80*x9; FIX(x139);
        /*@ assert (x11) != 0.;*/
        x151=1.0/x11; FIX(x151);
        x152=0.5*x151*x34; FIX(x152);
        x153=0.5*x151*x36*x69; FIX(x153);
        x158=-x139 - x152*x69*x75; FIX(x158);
        x160=-rt2 + x8; FIX(x160);
        x184=x137*x6*x75; FIX(x184);

        /* Assignment result[2, 5]=x153*x160 + x158 - x184 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        x56=mu*ut2; FIX(x56);
        x69=rt1 - x46; FIX(x69);
        x75=rt2 - x56; FIX(x75);
        x80=-rt1 + x5; FIX(x80);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        x137=0.5*x120*x14; FIX(x137);
        x139=x120*x19*x80*x9; FIX(x139);
        /*@ assert (x11) != 0.;*/
        x151=1.0/x11; FIX(x151);
        x152=0.5*x151*x34; FIX(x152);
        x153=0.5*x151*x36*x69; FIX(x153);
        x158=-x139 - x152*x69*x75; FIX(x158);
        x160=-rt2 + x8; FIX(x160);
        x184=x137*x6*x75; FIX(x184);
        result[17] = x153*x160 + x158 - x184;

    }
    else if (x145)
    {
        DEBUG_PRINT("Case (x145) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x69=rt1 - x46; FIX(x69);
        x70=0.5*x34*x49; FIX(x70);
        x71=x69*x70; FIX(x71);
        x72=-x71; FIX(x72);
        x73=0.5*x36*x49; FIX(x73);
        x74=x69*x73; FIX(x74);

        /* Assignment result[2, 5]=x72 - x74 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x69=rt1 - x46; FIX(x69);
        x70=0.5*x34*x49; FIX(x70);
        x71=x69*x70; FIX(x71);
        x72=-x71; FIX(x72);
        x73=0.5*x36*x49; FIX(x73);
        x74=x69*x73; FIX(x74);
        result[17] = x72 - x74;

    }


    /* Assignment result[0, 6]=Piecewise((0.0, x27), (rt2*x68, x64), (x77 + x78, x67)) */
    double x77;
    double x78;
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");

        /* Assignment result[0, 6]=0.0 */

        result[18] = 0.0;
        /*@ assert (result[18]) >= 0.;*/
    }
    else if (x64)
    {
        DEBUG_PRINT("Case (x64) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x68=0.5*x40*(-1.0*x28 + x31); FIX(x68);

        /* Assignment result[0, 6]=rt2*x68 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x26) != 0.;*/
        x40=1.0/x26; FIX(x40);
        x68=0.5*x40*(-1.0*x28 + x31); FIX(x68);
        result[18] = rt2*x68;

    }
    else if (x67)
    {
        DEBUG_PRINT("Case (x67) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x56=mu*ut2; FIX(x56);
        x70=0.5*x34*x49; FIX(x70);
        x73=0.5*x36*x49; FIX(x73);
        x75=rt2 - x56; FIX(x75);
        x76=x70*x75; FIX(x76);
        x77=-x76; FIX(x77);
        x78=x73*x75; FIX(x78);

        /* Assignment result[0, 6]=x77 + x78 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x56=mu*ut2; FIX(x56);
        x70=0.5*x34*x49; FIX(x70);
        x73=0.5*x36*x49; FIX(x73);
        x75=rt2 - x56; FIX(x75);
        x76=x70*x75; FIX(x76);
        x77=-x76; FIX(x77);
        x78=x73*x75; FIX(x78);
        result[18] = x77 + x78;

    }


    /* Assignment result[1, 6]=Piecewise((x155, x27), (x157, x64), (-x138 + x158 + x159*x80, x144), (0, x145)) */
    double x159;
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x39=mu*x22; FIX(x39);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x155=x39*x96; FIX(x155);

        /* Assignment result[1, 6]=x155 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x39=mu*x22; FIX(x39);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x155=x39*x96; FIX(x155);
        result[19] = x155;

    }
    else if (x64)
    {
        DEBUG_PRINT("Case (x64) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x114=Max(0, x108 - x26); FIX(x114);
        x148=x114*x23; FIX(x148);
        x156=1.0*x109; FIX(x156);
        x157=-0.5*rt1*rt2*x99*(x114*x24 - x156*x23 - x156*x24 + x101*x28 + x101*x31 + x148); FIX(x157);

        /* Assignment result[1, 6]=x157 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x114=Max(0, x108 - x26); FIX(x114);
        x148=x114*x23; FIX(x148);
        x156=1.0*x109; FIX(x156);
        x157=-0.5*rt1*rt2*x99*(x114*x24 - x156*x23 - x156*x24 + x101*x28 + x101*x31 + x148); FIX(x157);
        result[19] = x157;

    }
    else if (x144)
    {
        DEBUG_PRINT("Case (x144) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        x56=mu*ut2; FIX(x56);
        x69=rt1 - x46; FIX(x69);
        x75=rt2 - x56; FIX(x75);
        x80=-rt1 + x5; FIX(x80);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        x137=0.5*x120*x14; FIX(x137);
        x138=x137*x69*x9; FIX(x138);
        x139=x120*x19*x80*x9; FIX(x139);
        /*@ assert (x11) != 0.;*/
        x151=1.0/x11; FIX(x151);
        x152=0.5*x151*x34; FIX(x152);
        x158=-x139 - x152*x69*x75; FIX(x158);
        x159=0.5*x151*x36*x75; FIX(x159);

        /* Assignment result[1, 6]=-x138 + x158 + x159*x80 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        x46=mu*ut1; FIX(x46);
        x56=mu*ut2; FIX(x56);
        x69=rt1 - x46; FIX(x69);
        x75=rt2 - x56; FIX(x75);
        x80=-rt1 + x5; FIX(x80);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        x137=0.5*x120*x14; FIX(x137);
        x138=x137*x69*x9; FIX(x138);
        x139=x120*x19*x80*x9; FIX(x139);
        /*@ assert (x11) != 0.;*/
        x151=1.0/x11; FIX(x151);
        x152=0.5*x151*x34; FIX(x152);
        x158=-x139 - x152*x69*x75; FIX(x158);
        x159=0.5*x151*x36*x75; FIX(x159);
        result[19] = -x138 + x158 + x159*x80;

    }
    else if (x145)
    {
        DEBUG_PRINT("Case (x145) is True.\n");

        /* Assignment result[1, 6]=0 */

        result[19] = 0;
        /*@ assert (result[19]) >= 0.;*/
    }


    /* Assignment result[2, 6]=Piecewise((1.0 - x124, x27), (x99*(-x127*x147 + x127*x150 - x146*x29 - x146*x32 + x149), x64), (-x152*x75**2 - x19*(x154 + x190) - x15*(x191 + x49) + 1 + x159*x160, x144), (1 + x77 - x78, x145)) */

    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x124=1.0*x22*x96; FIX(x124);

        /* Assignment result[2, 6]=1.0 - x124 */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x22=Heaviside(x21); FIX(x22);
        x92=mu*mu; FIX(x92);
        x93=1.0 + x92; FIX(x93);
        /*@ assert (x93) != 0.;*/
        x96=1.0/x93; FIX(x96);
        x124=1.0*x22*x96; FIX(x124);
        result[20] = 1.0 - x124;

    }
    else if (x64)
    {
        DEBUG_PRINT("Case (x64) is True.\n");
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x32=0.5*x31; FIX(x32);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x114=Max(0, x108 - x26); FIX(x114);
        x119=x101*x23; FIX(x119);
        x127=pow(rt1, 4); FIX(x127);
        x130=x23*x24; FIX(x130);
        x146=x101*x24; FIX(x146);
        x147=0.5*x109; FIX(x147);
        x148=x114*x23; FIX(x148);
        x149=0.5*x148*x24 + x119 - x130*x147 + x146; FIX(x149);
        x150=0.5*x114; FIX(x150);

        /* Assignment result[2, 6]=x99*(-x127*x147 + x127*x150 - x146*x29 - x146*x32 + x149) */
        x20=mu*rn; FIX(x20);
        x21=-1.0*un + x20; FIX(x21);
        x28=Heaviside(x21 + x26); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x30=-1.0*x26; FIX(x30);
        x31=Heaviside(x21 + x30); FIX(x31);
        x32=0.5*x31; FIX(x32);
        /*@ assert (x25) >= 0.;*/
        /*@ assert (x25) != 0.;*/
        x99=pow(x25, -5.0/2.0); FIX(x99);
        /*@ assert (x25) >= 0.;*/
        x101=pow(x25, 3.0/2.0); FIX(x101);
        x108=x20 + x2; FIX(x108);
        x109=Max(0, x108 + x26); FIX(x109);
        x114=Max(0, x108 - x26); FIX(x114);
        x119=x101*x23; FIX(x119);
        x127=pow(rt1, 4); FIX(x127);
        x130=x23*x24; FIX(x130);
        x146=x101*x24; FIX(x146);
        x147=0.5*x109; FIX(x147);
        x148=x114*x23; FIX(x148);
        x149=0.5*x148*x24 + x119 - x130*x147 + x146; FIX(x149);
        x150=0.5*x114; FIX(x150);
        result[20] = x99*(-x127*x147 + x127*x150 - x146*x29 - x146*x32 + x149);

    }
    else if (x144)
    {
        DEBUG_PRINT("Case (x144) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x56=mu*ut2; FIX(x56);
        x75=rt2 - x56; FIX(x75);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        /*@ assert (x11) != 0.;*/
        x151=1.0/x11; FIX(x151);
        x152=0.5*x151*x34; FIX(x152);
        x154=-x49; FIX(x154);
        x159=0.5*x151*x36*x75; FIX(x159);
        x160=-rt2 + x8; FIX(x160);
        x190=x10*x120; FIX(x190);
        x191=x120*x75*x9; FIX(x191);

        /* Assignment result[2, 6]=-x152*x75**2 - x19*(x154 + x190) - x15*(x191 + x49) + 1 + x159*x160 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x56=mu*ut2; FIX(x56);
        x75=rt2 - x56; FIX(x75);
        /*@ assert (x11) >= 0.;*/
        /*@ assert (x11) != 0.;*/
        x120=pow(x11, -3.0/2.0); FIX(x120);
        /*@ assert (x11) != 0.;*/
        x151=1.0/x11; FIX(x151);
        x152=0.5*x151*x34; FIX(x152);
        x154=-x49; FIX(x154);
        x159=0.5*x151*x36*x75; FIX(x159);
        x160=-rt2 + x8; FIX(x160);
        x190=x10*x120; FIX(x190);
        x191=x120*x75*x9; FIX(x191);
        result[20] = -x152*x75*x75 - x19*(x154 + x190) - x15*(x191 + x49) + 1 + x159*x160;

    }
    else if (x145)
    {
        DEBUG_PRINT("Case (x145) is True.\n");
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x56=mu*ut2; FIX(x56);
        x70=0.5*x34*x49; FIX(x70);
        x73=0.5*x36*x49; FIX(x73);
        x75=rt2 - x56; FIX(x75);
        x76=x70*x75; FIX(x76);
        x77=-x76; FIX(x77);
        x78=x73*x75; FIX(x78);

        /* Assignment result[2, 6]=1 + x77 - x78 */
        x34=Heaviside(x13); FIX(x34);
        x36=Heaviside(x17); FIX(x36);
        /*@ assert (x12) != 0.;*/
        x49=1.0/x12; FIX(x49);
        x56=mu*ut2; FIX(x56);
        x70=0.5*x34*x49; FIX(x70);
        x73=0.5*x36*x49; FIX(x73);
        x75=rt2 - x56; FIX(x75);
        x76=x70*x75; FIX(x76);
        x77=-x76; FIX(x77);
        x78=x73*x75; FIX(x78);
        result[20] = 1 + x77 - x78;

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
    // Heaviside
    double x4;
    double x5;
    double x6;
    double x7;
    double x8;
    int x9;
    double x1;
    double x2;
    double x3;
    x4=rt1*rt1; FIX(x4);
    /*@ assert (x4) >= 0.;*/
    x5=rt2*rt2; FIX(x5);
    /*@ assert (x5) >= 0.;*/
    x6=x4 + x5; FIX(x6);
    /*@ assert (x6) >= 0.;*/
    x7=sqrt(x6); FIX(x7);
    /*@ assert (ut1*ut1 + ut2*ut2) >= 0.;*/
    x8=sqrt(ut1*ut1 + ut2*ut2); FIX(x8);
    /*@ assert (x8) >= 0.;*/
    x9=x8 + x7 <= 0; FIX(x9);
    int x15;
    double x10;
    double x11;
    double x12;
    double x13;
    double x14;
    x15=x8 <= 0; FIX(x15);
    int x30;
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
    x30=x8 > 0; FIX(x30);
    int x56;
    x18=mu*ut1 - rt1; FIX(x18);
    x19=x18*x18; FIX(x19);
    x20=mu*ut2 - rt2; FIX(x20);
    x21=x20*x20; FIX(x21);
    x22=x19 + x21; FIX(x22);
    x56=x22 <= 0; FIX(x56);
    int x59;
    double x57;
    double x58;
    x59=x22 > 0; FIX(x59);
    int x77;
    int x78;
    double x38;
    double x41;
    double x61;
    double x62;
    double x63;
    double x65;
    double x76;
    /*@ assert (x22) >= 0.;*/
    x23=sqrt(x22); FIX(x23);
    x77=x23 > 0; FIX(x77);
    x78=x30 && x77; FIX(x78);
    int x79;
    int x80;
    x79=x23 <= 0; FIX(x79);
    x80=x30 && x79; FIX(x80);
    int x137;
    double x42;
    double x43;
    double x136;
    x137=x59 && x77; FIX(x137);
    int x138;
    x138=x59 && x79; FIX(x138);
    if (x9)
    {
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
    }
    else if (x15)
    {
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x11=0.5*x10; FIX(x11);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x14=0.5*x13; FIX(x14);
    }
    else if (x30)
    {
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
    }
    else if (x56)
    {
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x11=0.5*x10; FIX(x11);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x14=0.5*x13; FIX(x14);
    }
    else if (x59)
    {
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x57=-mu*x26; FIX(x57);
        x58=mu*x29; FIX(x58);
    }
    else if (x78)
    {
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x61=rt1 - x38; FIX(x61);
        x62=0.5*x41*x25; FIX(x62);
        x63=x61*x62; FIX(x63);
        x65=0.5*x28*x41; FIX(x65);
        x76=x18*x65; FIX(x76);
    }
    else if (x137)
    {
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x43=x18*x42; FIX(x43);
        x61=rt1 - x38; FIX(x61);
        x136=0.5*mu*x41*x25; FIX(x136);
    }
    /* Assignment result[0, 0]=Piecewise((x3, x9), (x11 + x14, x15), (x26 + x29, x30)) */

    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);

        /* Assignment result[0, 0]=x3 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        result[0] = x3;

    }
    else if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x11=0.5*x10; FIX(x11);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x14=0.5*x13; FIX(x14);

        /* Assignment result[0, 0]=x11 + x14 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x11=0.5*x10; FIX(x11);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x14=0.5*x13; FIX(x14);
        result[0] = x11 + x14;

    }
    else if (x30)
    {
        DEBUG_PRINT("Case (x30) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);

        /* Assignment result[0, 0]=x26 + x29 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        result[0] = x26 + x29;

    }
    /*@ assert (result[0]) >= 0.;*/

    /* Assignment result[1, 0]=Piecewise((0.0, x9), (-0.5*x72*x32*(x74 - 1.0*x75), x15), (x63 + x76, x78), (0, x80)) */
    double x32;
    double x55;
    double x71;
    double x72;
    double x73;
    double x74;
    double x75;
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");

        /* Assignment result[1, 0]=0.0 */

        result[1] = 0.0;
        /*@ assert (result[1]) >= 0.;*/
    }
    else if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");
        x12=-1.0*x7; FIX(x12);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x55=1.0*mu; FIX(x55);
        x71=Heaviside(x7); FIX(x71);
        x72=rt1*x71; FIX(x72);
        x73=-rn*x55 + un; FIX(x73);
        x74=Heaviside(x12 + x73); FIX(x74);
        x75=Heaviside(x7 + x73); FIX(x75);

        /* Assignment result[1, 0]=-0.5*x72*x32*(x74 - 1.0*x75) */
        x12=-1.0*x7; FIX(x12);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x55=1.0*mu; FIX(x55);
        x71=Heaviside(x7); FIX(x71);
        x72=rt1*x71; FIX(x72);
        x73=-rn*x55 + un; FIX(x73);
        x74=Heaviside(x12 + x73); FIX(x74);
        x75=Heaviside(x7 + x73); FIX(x75);
        result[1] = -0.5*x72*x32*(x74 - 1.0*x75);

    }
    else if (x78)
    {
        DEBUG_PRINT("Case (x78) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x61=rt1 - x38; FIX(x61);
        x62=0.5*x41*x25; FIX(x62);
        x63=x61*x62; FIX(x63);
        x65=0.5*x28*x41; FIX(x65);
        x76=x18*x65; FIX(x76);

        /* Assignment result[1, 0]=x63 + x76 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x61=rt1 - x38; FIX(x61);
        x62=0.5*x41*x25; FIX(x62);
        x63=x61*x62; FIX(x63);
        x65=0.5*x28*x41; FIX(x65);
        x76=x18*x65; FIX(x76);
        result[1] = x63 + x76;

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");

        /* Assignment result[1, 0]=0 */

        result[1] = 0;
        /*@ assert (result[1]) >= 0.;*/
    }
    /*@ assert (result[1]) >= 0.;*/

    /* Assignment result[2, 0]=Piecewise((0.0, x9), (x32*(x74*x154 - x75*x154 + x153*x126 - x153*x128 + x7*x11 - x7*x14), x15), (x155 + x68, x78), (x26 - x29, x80)) */
    double x48;
    double x67;
    double x68;
    double x126;
    double x128;
    double x153;
    double x154;
    double x155;
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");

        /* Assignment result[2, 0]=0.0 */

        result[2] = 0.0;
        /*@ assert (result[2]) >= 0.;*/
    }
    else if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x11=0.5*x10; FIX(x11);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x14=0.5*x13; FIX(x14);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x55=1.0*mu; FIX(x55);
        x71=Heaviside(x7); FIX(x71);
        x73=-rn*x55 + un; FIX(x73);
        x74=Heaviside(x12 + x73); FIX(x74);
        x75=Heaviside(x7 + x73); FIX(x75);
        x126=rt2*x75; FIX(x126);
        x128=rt2*x74; FIX(x128);
        x153=0.5*x71; FIX(x153);
        x154=0.5*x7*x71; FIX(x154);

        /* Assignment result[2, 0]=x32*(x74*x154 - x75*x154 + x153*x126 - x153*x128 + x7*x11 - x7*x14) */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x11=0.5*x10; FIX(x11);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x14=0.5*x13; FIX(x14);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x55=1.0*mu; FIX(x55);
        x71=Heaviside(x7); FIX(x71);
        x73=-rn*x55 + un; FIX(x73);
        x74=Heaviside(x12 + x73); FIX(x74);
        x75=Heaviside(x7 + x73); FIX(x75);
        x126=rt2*x75; FIX(x126);
        x128=rt2*x74; FIX(x128);
        x153=0.5*x71; FIX(x153);
        x154=0.5*x7*x71; FIX(x154);
        result[2] = x32*(x74*x154 - x75*x154 + x153*x126 - x153*x128 + x7*x11 - x7*x14);

    }
    else if (x78)
    {
        DEBUG_PRINT("Case (x78) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x48=mu*ut2; FIX(x48);
        x62=0.5*x41*x25; FIX(x62);
        x65=0.5*x28*x41; FIX(x65);
        x67=rt2 - x48; FIX(x67);
        x68=x67*x62; FIX(x68);
        x155=x20*x65; FIX(x155);

        /* Assignment result[2, 0]=x155 + x68 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x48=mu*ut2; FIX(x48);
        x62=0.5*x41*x25; FIX(x62);
        x65=0.5*x28*x41; FIX(x65);
        x67=rt2 - x48; FIX(x67);
        x68=x67*x62; FIX(x68);
        x155=x20*x65; FIX(x155);
        result[2] = x155 + x68;

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);

        /* Assignment result[2, 0]=x26 - x29 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        result[2] = x26 - x29;

    }
    /*@ assert (result[2]) >= 0.;*/

    /* Assignment result[0, 1]=Piecewise((x31, x9), (x33*(rt1*x34 - rt1*x35 + x37), x15), (x45 - x47, x30)) */
    double x31;
    double x33;
    double x34;
    double x35;
    double x36;
    double x37;
    double x39;
    double x40;
    double x44;
    double x45;
    double x46;
    double x47;
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x31=mu*x3; FIX(x31);

        /* Assignment result[0, 1]=x31 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x31=mu*x3; FIX(x31);
        result[3] = x31;

    }
    else if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x33=0.25*mu*x32; FIX(x33);
        x34=2.0*x10; FIX(x34);
        x35=2.0*x13; FIX(x35);
        x36=1.4142135623730951454746218587388284504413604736328125*x7; FIX(x36);
        x37=x36*x10 + x36*x13; FIX(x37);

        /* Assignment result[0, 1]=x33*(rt1*x34 - rt1*x35 + x37) */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x33=0.25*mu*x32; FIX(x33);
        x34=2.0*x10; FIX(x34);
        x35=2.0*x13; FIX(x35);
        x36=1.4142135623730951454746218587388284504413604736328125*x7; FIX(x36);
        x37=x36*x10 + x36*x13; FIX(x37);
        result[3] = x33*(rt1*x34 - rt1*x35 + x37);

    }
    else if (x30)
    {
        DEBUG_PRINT("Case (x30) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        x40=-x38*x39; FIX(x40);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x43=x18*x42; FIX(x43);
        x44=x40 + x43; FIX(x44);
        x45=-x26*x44; FIX(x45);
        x46=x40 - x43; FIX(x46);
        x47=x46*x29; FIX(x47);

        /* Assignment result[0, 1]=x45 - x47 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        x40=-x38*x39; FIX(x40);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x43=x18*x42; FIX(x43);
        x44=x40 + x43; FIX(x44);
        x45=-x26*x44; FIX(x45);
        x46=x40 - x43; FIX(x46);
        x47=x46*x29; FIX(x47);
        result[3] = x45 - x47;

    }
    /*@ assert (result[3]) >= 0.;*/

    /* Assignment result[1, 1]=Piecewise((x84*x86 - x87*x86, x9), (-x89*x71*(-x4*x91 - x5*x105 + x5*x107 + 2.0*x108*x74 + 2.0*x108*x75 - x101*x102 + x101*x106 - x92*x93 + x92*x96 + x94*x104 - x94*x99 - x95*x93 + x95*x96), x15), (-x110*(-mu*x112 + x42) - x114*(-mu*x116 + x115) - x44*x63 - x46*x76, x78), (0, x80)) */
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
    double x104;
    double x105;
    double x106;
    double x107;
    double x108;
    double x109;
    double x110;
    double x111;
    double x112;
    double x113;
    double x114;
    double x115;
    double x116;
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x55=1.0*mu; FIX(x55);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) >= 0.;*/
        x83=sqrt(x82); FIX(x83);
        x84=mu + x83; FIX(x84);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x86=0.5*x3*x81*x85; FIX(x86);
        x87=-x55 + x83; FIX(x87);

        /* Assignment result[1, 1]=x84*x86 - x87*x86 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x55=1.0*mu; FIX(x55);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) >= 0.;*/
        x83=sqrt(x82); FIX(x83);
        x84=mu + x83; FIX(x84);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x86=0.5*x3*x81*x85; FIX(x86);
        x87=-x55 + x83; FIX(x87);
        result[4] = x84*x86 - x87*x86;

    }
    else if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");
        x1=mu*rn; FIX(x1);
        x12=-1.0*x7; FIX(x12);
        x16=-un; FIX(x16);
        x55=1.0*mu; FIX(x55);
        x71=Heaviside(x7); FIX(x71);
        x73=-rn*x55 + un; FIX(x73);
        x74=Heaviside(x12 + x73); FIX(x74);
        x75=Heaviside(x7 + x73); FIX(x75);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        x89=0.25*mu*x88; FIX(x89);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x91=4.0*x90; FIX(x91);
        x92=pow(rt1, 5); FIX(x92);
        x93=1.4142135623730951454746218587388284504413604736328125*x75; FIX(x93);
        x94=pow(rt2, 4); FIX(x94);
        x95=rt1*x94; FIX(x95);
        x96=1.4142135623730951454746218587388284504413604736328125*x74; FIX(x96);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x99=2.0*x98; FIX(x99);
        x100=pow(rt1, 3); FIX(x100);
        x101=x100*x5; FIX(x101);
        x102=2.828427124746190290949243717477656900882720947265625*x75; FIX(x102);
        x103=Max(0, -x7 + x97); FIX(x103);
        x104=2.0*x103; FIX(x104);
        x105=2.0*x4*x98; FIX(x105);
        x106=2.828427124746190290949243717477656900882720947265625*x74; FIX(x106);
        x107=2.0*x103*x4; FIX(x107);
        x108=x4*x90; FIX(x108);

        /* Assignment result[1, 1]=-x89*x71*(-x4*x91 - x5*x105 + x5*x107 + 2.0*x108*x74 + 2.0*x108*x75 - x101*x102 + x101*x106 - x92*x93 + x92*x96 + x94*x104 - x94*x99 - x95*x93 + x95*x96) */
        x1=mu*rn; FIX(x1);
        x12=-1.0*x7; FIX(x12);
        x16=-un; FIX(x16);
        x55=1.0*mu; FIX(x55);
        x71=Heaviside(x7); FIX(x71);
        x73=-rn*x55 + un; FIX(x73);
        x74=Heaviside(x12 + x73); FIX(x74);
        x75=Heaviside(x7 + x73); FIX(x75);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        x89=0.25*mu*x88; FIX(x89);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x91=4.0*x90; FIX(x91);
        x92=pow(rt1, 5); FIX(x92);
        x93=1.4142135623730951454746218587388284504413604736328125*x75; FIX(x93);
        x94=pow(rt2, 4); FIX(x94);
        x95=rt1*x94; FIX(x95);
        x96=1.4142135623730951454746218587388284504413604736328125*x74; FIX(x96);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x99=2.0*x98; FIX(x99);
        x100=pow(rt1, 3); FIX(x100);
        x101=x100*x5; FIX(x101);
        x102=2.828427124746190290949243717477656900882720947265625*x75; FIX(x102);
        x103=Max(0, -x7 + x97); FIX(x103);
        x104=2.0*x103; FIX(x104);
        x105=2.0*x4*x98; FIX(x105);
        x106=2.828427124746190290949243717477656900882720947265625*x74; FIX(x106);
        x107=2.0*x103*x4; FIX(x107);
        x108=x4*x90; FIX(x108);
        result[4] = -x89*x71*(-x4*x91 - x5*x105 + x5*x107 + 2.0*x108*x74 + 2.0*x108*x75 - x101*x102 + x101*x106 - x92*x93 + x92*x96 + x94*x104 - x94*x99 - x95*x93 + x95*x96);

    }
    else if (x78)
    {
        DEBUG_PRINT("Case (x78) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        x40=-x38*x39; FIX(x40);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x43=x18*x42; FIX(x43);
        x44=x40 + x43; FIX(x44);
        x46=x40 - x43; FIX(x46);
        x61=rt1 - x38; FIX(x61);
        x62=0.5*x41*x25; FIX(x62);
        x63=x61*x62; FIX(x63);
        x65=0.5*x28*x41; FIX(x65);
        x76=x18*x65; FIX(x76);
        x109=Max(0, x27); FIX(x109);
        x110=0.5*x109; FIX(x110);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x112=x19*x111; FIX(x112);
        x113=x24; FIX(x113);
        x114=0.5*x113; FIX(x114);
        x115=-x42; FIX(x115);
        x116=x18*x61*x111; FIX(x116);

        /* Assignment result[1, 1]=-x110*(-mu*x112 + x42) - x114*(-mu*x116 + x115) - x44*x63 - x46*x76 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        x40=-x38*x39; FIX(x40);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x43=x18*x42; FIX(x43);
        x44=x40 + x43; FIX(x44);
        x46=x40 - x43; FIX(x46);
        x61=rt1 - x38; FIX(x61);
        x62=0.5*x41*x25; FIX(x62);
        x63=x61*x62; FIX(x63);
        x65=0.5*x28*x41; FIX(x65);
        x76=x18*x65; FIX(x76);
        x109=Max(0, x27); FIX(x109);
        x110=0.5*x109; FIX(x110);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x112=x19*x111; FIX(x112);
        x113=x24; FIX(x113);
        x114=0.5*x113; FIX(x114);
        x115=-x42; FIX(x115);
        x116=x18*x61*x111; FIX(x116);
        result[4] = -x110*(-mu*x112 + x42) - x114*(-mu*x116 + x115) - x44*x63 - x46*x76;

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");

        /* Assignment result[1, 1]=0 */

        result[4] = 0;
        /*@ assert (result[4]) >= 0.;*/
    }
    /*@ assert (result[4]) >= 0.;*/

    /* Assignment result[2, 1]=Piecewise((-x84*x156 + x87*x156, x9), (-x89*(-x173*x74 - x75*x173 + x101*x162 - x101*x163 - x101*x164 + x139*x159 + x160 - x167*x93 + x167*x96 - x168*x93 + x168*x96 - x170*x102 + x170*x106 + x171*x172 - x171*x174 - x34*x92 - x34*x95 - x72*x119 + x72*x122 - x72*x124 + x72*x127 + x72*x129 + x92*x161 - x92*x165 - x92*x169 - x92*x35 + x95*x161 - x95*x165 - x95*x169 - x95*x35 - x96*x175), x15), (mu*x176 + x134 - x44*x68 - x46*x155, x78), (x45 + x47, x80)) */
    double x119;
    double x120;
    double x121;
    double x122;
    double x124;
    double x125;
    double x127;
    double x129;
    double x131;
    double x133;
    double x134;
    double x139;
    double x152;
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
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x55=1.0*mu; FIX(x55);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) >= 0.;*/
        x83=sqrt(x82); FIX(x83);
        x84=mu + x83; FIX(x84);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x87=-x55 + x83; FIX(x87);
        x156=0.5*mu*x3*x85; FIX(x156);

        /* Assignment result[2, 1]=-x84*x156 + x87*x156 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x55=1.0*mu; FIX(x55);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) >= 0.;*/
        x83=sqrt(x82); FIX(x83);
        x84=mu + x83; FIX(x84);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x87=-x55 + x83; FIX(x87);
        x156=0.5*mu*x3*x85; FIX(x156);
        result[5] = -x84*x156 + x87*x156;

    }
    else if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x16=-un; FIX(x16);
        x34=2.0*x10; FIX(x34);
        x35=2.0*x13; FIX(x35);
        x55=1.0*mu; FIX(x55);
        x71=Heaviside(x7); FIX(x71);
        x72=rt1*x71; FIX(x72);
        x73=-rn*x55 + un; FIX(x73);
        x74=Heaviside(x12 + x73); FIX(x74);
        x75=Heaviside(x7 + x73); FIX(x75);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        x89=0.25*mu*x88; FIX(x89);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x91=4.0*x90; FIX(x91);
        x92=pow(rt1, 5); FIX(x92);
        x93=1.4142135623730951454746218587388284504413604736328125*x75; FIX(x93);
        x94=pow(rt2, 4); FIX(x94);
        x95=rt1*x94; FIX(x95);
        x96=1.4142135623730951454746218587388284504413604736328125*x74; FIX(x96);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x99=2.0*x98; FIX(x99);
        x100=pow(rt1, 3); FIX(x100);
        x101=x100*x5; FIX(x101);
        x102=2.828427124746190290949243717477656900882720947265625*x75; FIX(x102);
        x103=Max(0, -x7 + x97); FIX(x103);
        x104=2.0*x103; FIX(x104);
        x106=2.828427124746190290949243717477656900882720947265625*x74; FIX(x106);
        x108=x4*x90; FIX(x108);
        x119=rt2*x91; FIX(x119);
        x120=pow(rt1, 4); FIX(x120);
        x121=pow(rt2, 3); FIX(x121);
        x122=x121*x99; FIX(x122);
        x124=x121*x104; FIX(x124);
        x125=2.0*x90; FIX(x125);
        x126=rt2*x75; FIX(x126);
        x127=x125*x126; FIX(x127);
        x128=rt2*x74; FIX(x128);
        x129=x125*x128; FIX(x129);
        x139=x90*x5; FIX(x139);
        x157=1.4142135623730951454746218587388284504413604736328125*x10; FIX(x157);
        x158=1.4142135623730951454746218587388284504413604736328125*x13; FIX(x158);
        x159=1.4142135623730951454746218587388284504413604736328125*x71*x75; FIX(x159);
        x160=-x4*x90*x96*x71 - x108*x157 + x108*x158 + x108*x159 - x139*x157 + x139*x158; FIX(x160);
        x161=4.0*x71; FIX(x161);
        x162=8.0*x71; FIX(x162);
        x163=4.0*x10; FIX(x163);
        x164=4.0*x13; FIX(x164);
        x165=2.0*x71*x75; FIX(x165);
        x166=pow(rt2, 5); FIX(x166);
        x167=x166*x71; FIX(x167);
        x168=rt2*x120*x71; FIX(x168);
        x169=2.0*x71*x74; FIX(x169);
        x170=x4*x121*x71; FIX(x170);
        x171=rt2*x100; FIX(x171);
        x172=2.0*x71*x98; FIX(x172);
        x173=4.0*x100*x5*x71; FIX(x173);
        x174=2.0*x103*x71; FIX(x174);
        x175=x90*x5*x71; FIX(x175);

        /* Assignment result[2, 1]=-x89*(-x173*x74 - x75*x173 + x101*x162 - x101*x163 - x101*x164 + x139*x159 + x160 - x167*x93 + x167*x96 - x168*x93 + x168*x96 - x170*x102 + x170*x106 + x171*x172 - x171*x174 - x34*x92 - x34*x95 - x72*x119 + x72*x122 - x72*x124 + x72*x127 + x72*x129 + x92*x161 - x92*x165 - x92*x169 - x92*x35 + x95*x161 - x95*x165 - x95*x169 - x95*x35 - x96*x175) */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x16=-un; FIX(x16);
        x34=2.0*x10; FIX(x34);
        x35=2.0*x13; FIX(x35);
        x55=1.0*mu; FIX(x55);
        x71=Heaviside(x7); FIX(x71);
        x72=rt1*x71; FIX(x72);
        x73=-rn*x55 + un; FIX(x73);
        x74=Heaviside(x12 + x73); FIX(x74);
        x75=Heaviside(x7 + x73); FIX(x75);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        x89=0.25*mu*x88; FIX(x89);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x91=4.0*x90; FIX(x91);
        x92=pow(rt1, 5); FIX(x92);
        x93=1.4142135623730951454746218587388284504413604736328125*x75; FIX(x93);
        x94=pow(rt2, 4); FIX(x94);
        x95=rt1*x94; FIX(x95);
        x96=1.4142135623730951454746218587388284504413604736328125*x74; FIX(x96);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x99=2.0*x98; FIX(x99);
        x100=pow(rt1, 3); FIX(x100);
        x101=x100*x5; FIX(x101);
        x102=2.828427124746190290949243717477656900882720947265625*x75; FIX(x102);
        x103=Max(0, -x7 + x97); FIX(x103);
        x104=2.0*x103; FIX(x104);
        x106=2.828427124746190290949243717477656900882720947265625*x74; FIX(x106);
        x108=x4*x90; FIX(x108);
        x119=rt2*x91; FIX(x119);
        x120=pow(rt1, 4); FIX(x120);
        x121=pow(rt2, 3); FIX(x121);
        x122=x121*x99; FIX(x122);
        x124=x121*x104; FIX(x124);
        x125=2.0*x90; FIX(x125);
        x126=rt2*x75; FIX(x126);
        x127=x125*x126; FIX(x127);
        x128=rt2*x74; FIX(x128);
        x129=x125*x128; FIX(x129);
        x139=x90*x5; FIX(x139);
        x157=1.4142135623730951454746218587388284504413604736328125*x10; FIX(x157);
        x158=1.4142135623730951454746218587388284504413604736328125*x13; FIX(x158);
        x159=1.4142135623730951454746218587388284504413604736328125*x71*x75; FIX(x159);
        x160=-x4*x90*x96*x71 - x108*x157 + x108*x158 + x108*x159 - x139*x157 + x139*x158; FIX(x160);
        x161=4.0*x71; FIX(x161);
        x162=8.0*x71; FIX(x162);
        x163=4.0*x10; FIX(x163);
        x164=4.0*x13; FIX(x164);
        x165=2.0*x71*x75; FIX(x165);
        x166=pow(rt2, 5); FIX(x166);
        x167=x166*x71; FIX(x167);
        x168=rt2*x120*x71; FIX(x168);
        x169=2.0*x71*x74; FIX(x169);
        x170=x4*x121*x71; FIX(x170);
        x171=rt2*x100; FIX(x171);
        x172=2.0*x71*x98; FIX(x172);
        x173=4.0*x100*x5*x71; FIX(x173);
        x174=2.0*x103*x71; FIX(x174);
        x175=x90*x5*x71; FIX(x175);
        result[5] = -x89*(-x173*x74 - x75*x173 + x101*x162 - x101*x163 - x101*x164 + x139*x159 + x160 - x167*x93 + x167*x96 - x168*x93 + x168*x96 - x170*x102 + x170*x106 + x171*x172 - x171*x174 - x34*x92 - x34*x95 - x72*x119 + x72*x122 - x72*x124 + x72*x127 + x72*x129 + x92*x161 - x92*x165 - x92*x169 - x92*x35 + x95*x161 - x95*x165 - x95*x169 - x95*x35 - x96*x175);

    }
    else if (x78)
    {
        DEBUG_PRINT("Case (x78) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        x40=-x38*x39; FIX(x40);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x43=x18*x42; FIX(x43);
        x44=x40 + x43; FIX(x44);
        x46=x40 - x43; FIX(x46);
        x48=mu*ut2; FIX(x48);
        x62=0.5*x41*x25; FIX(x62);
        x65=0.5*x28*x41; FIX(x65);
        x67=rt2 - x48; FIX(x67);
        x68=x67*x62; FIX(x68);
        x109=Max(0, x27); FIX(x109);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x113=x24; FIX(x113);
        x131=0.5*x111*x113; FIX(x131);
        x133=0.5*x109*x18*x20*x111; FIX(x133);
        x134=mu*x133; FIX(x134);
        x152=x18*x67; FIX(x152);
        x155=x20*x65; FIX(x155);
        x176=x152*x131; FIX(x176);

        /* Assignment result[2, 1]=mu*x176 + x134 - x44*x68 - x46*x155 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        x40=-x38*x39; FIX(x40);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x43=x18*x42; FIX(x43);
        x44=x40 + x43; FIX(x44);
        x46=x40 - x43; FIX(x46);
        x48=mu*ut2; FIX(x48);
        x62=0.5*x41*x25; FIX(x62);
        x65=0.5*x28*x41; FIX(x65);
        x67=rt2 - x48; FIX(x67);
        x68=x67*x62; FIX(x68);
        x109=Max(0, x27); FIX(x109);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x113=x24; FIX(x113);
        x131=0.5*x111*x113; FIX(x131);
        x133=0.5*x109*x18*x20*x111; FIX(x133);
        x134=mu*x133; FIX(x134);
        x152=x18*x67; FIX(x152);
        x155=x20*x65; FIX(x155);
        x176=x152*x131; FIX(x176);
        result[5] = mu*x176 + x134 - x44*x68 - x46*x155;

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        x40=-x38*x39; FIX(x40);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x43=x18*x42; FIX(x43);
        x44=x40 + x43; FIX(x44);
        x45=-x26*x44; FIX(x45);
        x46=x40 - x43; FIX(x46);
        x47=x46*x29; FIX(x47);

        /* Assignment result[2, 1]=x45 + x47 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        x40=-x38*x39; FIX(x40);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x43=x18*x42; FIX(x43);
        x44=x40 + x43; FIX(x44);
        x45=-x26*x44; FIX(x45);
        x46=x40 - x43; FIX(x46);
        x47=x46*x29; FIX(x47);
        result[5] = x45 + x47;

    }
    /*@ assert (result[5]) >= 0.;*/

    /* Assignment result[0, 2]=Piecewise((0.0, x9), (x33*(rt2*x34 - rt2*x35 + x37), x15), (x52 - x54, x30)) */
    double x49;
    double x50;
    double x51;
    double x52;
    double x53;
    double x54;
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");

        /* Assignment result[0, 2]=0.0 */

        result[6] = 0.0;
        /*@ assert (result[6]) >= 0.;*/
    }
    else if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x33=0.25*mu*x32; FIX(x33);
        x34=2.0*x10; FIX(x34);
        x35=2.0*x13; FIX(x35);
        x36=1.4142135623730951454746218587388284504413604736328125*x7; FIX(x36);
        x37=x36*x10 + x36*x13; FIX(x37);

        /* Assignment result[0, 2]=x33*(rt2*x34 - rt2*x35 + x37) */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x33=0.25*mu*x32; FIX(x33);
        x34=2.0*x10; FIX(x34);
        x35=2.0*x13; FIX(x35);
        x36=1.4142135623730951454746218587388284504413604736328125*x7; FIX(x36);
        x37=x36*x10 + x36*x13; FIX(x37);
        result[6] = x33*(rt2*x34 - rt2*x35 + x37);

    }
    else if (x30)
    {
        DEBUG_PRINT("Case (x30) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x48=mu*ut2; FIX(x48);
        x49=-x48*x39; FIX(x49);
        x50=x20*x42; FIX(x50);
        x51=x49 + x50; FIX(x51);
        x52=-x26*x51; FIX(x52);
        x53=x49 - x50; FIX(x53);
        x54=x53*x29; FIX(x54);

        /* Assignment result[0, 2]=x52 - x54 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x48=mu*ut2; FIX(x48);
        x49=-x48*x39; FIX(x49);
        x50=x20*x42; FIX(x50);
        x51=x49 + x50; FIX(x51);
        x52=-x26*x51; FIX(x52);
        x53=x49 - x50; FIX(x53);
        x54=x53*x29; FIX(x54);
        result[6] = x52 - x54;

    }
    /*@ assert (result[6]) >= 0.;*/

    /* Assignment result[1, 2]=Piecewise((x118, x9), (-rt1*x89*x71*(rt2*x105 - rt2*x107 - x119 - x120*x93 + x120*x96 + x122 - x123*x102 + x123*x106 - x124 + x127 + x129 - x94*x93 + x94*x96), x15), (mu*x132 + x134 - x51*x63 - x53*x76, x78), (0, x80)) */
    double x117;
    double x118;
    double x123;
    double x130;
    double x132;
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x117=1.0*x3*x85; FIX(x117);
        x118=-x81*x117; FIX(x118);

        /* Assignment result[1, 2]=x118 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x117=1.0*x3*x85; FIX(x117);
        x118=-x81*x117; FIX(x118);
        result[7] = x118;

    }
    else if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");
        x1=mu*rn; FIX(x1);
        x12=-1.0*x7; FIX(x12);
        x16=-un; FIX(x16);
        x55=1.0*mu; FIX(x55);
        x71=Heaviside(x7); FIX(x71);
        x73=-rn*x55 + un; FIX(x73);
        x74=Heaviside(x12 + x73); FIX(x74);
        x75=Heaviside(x7 + x73); FIX(x75);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        x89=0.25*mu*x88; FIX(x89);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x91=4.0*x90; FIX(x91);
        x93=1.4142135623730951454746218587388284504413604736328125*x75; FIX(x93);
        x94=pow(rt2, 4); FIX(x94);
        x96=1.4142135623730951454746218587388284504413604736328125*x74; FIX(x96);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x99=2.0*x98; FIX(x99);
        x102=2.828427124746190290949243717477656900882720947265625*x75; FIX(x102);
        x103=Max(0, -x7 + x97); FIX(x103);
        x104=2.0*x103; FIX(x104);
        x105=2.0*x4*x98; FIX(x105);
        x106=2.828427124746190290949243717477656900882720947265625*x74; FIX(x106);
        x107=2.0*x103*x4; FIX(x107);
        x119=rt2*x91; FIX(x119);
        x120=pow(rt1, 4); FIX(x120);
        x121=pow(rt2, 3); FIX(x121);
        x122=x121*x99; FIX(x122);
        x123=x4*x5; FIX(x123);
        x124=x121*x104; FIX(x124);
        x125=2.0*x90; FIX(x125);
        x126=rt2*x75; FIX(x126);
        x127=x125*x126; FIX(x127);
        x128=rt2*x74; FIX(x128);
        x129=x125*x128; FIX(x129);

        /* Assignment result[1, 2]=-rt1*x89*x71*(rt2*x105 - rt2*x107 - x119 - x120*x93 + x120*x96 + x122 - x123*x102 + x123*x106 - x124 + x127 + x129 - x94*x93 + x94*x96) */
        x1=mu*rn; FIX(x1);
        x12=-1.0*x7; FIX(x12);
        x16=-un; FIX(x16);
        x55=1.0*mu; FIX(x55);
        x71=Heaviside(x7); FIX(x71);
        x73=-rn*x55 + un; FIX(x73);
        x74=Heaviside(x12 + x73); FIX(x74);
        x75=Heaviside(x7 + x73); FIX(x75);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        x89=0.25*mu*x88; FIX(x89);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x91=4.0*x90; FIX(x91);
        x93=1.4142135623730951454746218587388284504413604736328125*x75; FIX(x93);
        x94=pow(rt2, 4); FIX(x94);
        x96=1.4142135623730951454746218587388284504413604736328125*x74; FIX(x96);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x99=2.0*x98; FIX(x99);
        x102=2.828427124746190290949243717477656900882720947265625*x75; FIX(x102);
        x103=Max(0, -x7 + x97); FIX(x103);
        x104=2.0*x103; FIX(x104);
        x105=2.0*x4*x98; FIX(x105);
        x106=2.828427124746190290949243717477656900882720947265625*x74; FIX(x106);
        x107=2.0*x103*x4; FIX(x107);
        x119=rt2*x91; FIX(x119);
        x120=pow(rt1, 4); FIX(x120);
        x121=pow(rt2, 3); FIX(x121);
        x122=x121*x99; FIX(x122);
        x123=x4*x5; FIX(x123);
        x124=x121*x104; FIX(x124);
        x125=2.0*x90; FIX(x125);
        x126=rt2*x75; FIX(x126);
        x127=x125*x126; FIX(x127);
        x128=rt2*x74; FIX(x128);
        x129=x125*x128; FIX(x129);
        result[7] = -rt1*x89*x71*(rt2*x105 - rt2*x107 - x119 - x120*x93 + x120*x96 + x122 - x123*x102 + x123*x106 - x124 + x127 + x129 - x94*x93 + x94*x96);

    }
    else if (x78)
    {
        DEBUG_PRINT("Case (x78) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x48=mu*ut2; FIX(x48);
        x49=-x48*x39; FIX(x49);
        x50=x20*x42; FIX(x50);
        x51=x49 + x50; FIX(x51);
        x53=x49 - x50; FIX(x53);
        x61=rt1 - x38; FIX(x61);
        x62=0.5*x41*x25; FIX(x62);
        x63=x61*x62; FIX(x63);
        x65=0.5*x28*x41; FIX(x65);
        x76=x18*x65; FIX(x76);
        x109=Max(0, x27); FIX(x109);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x113=x24; FIX(x113);
        x130=x20*x61; FIX(x130);
        x131=0.5*x111*x113; FIX(x131);
        x132=x130*x131; FIX(x132);
        x133=0.5*x109*x18*x20*x111; FIX(x133);
        x134=mu*x133; FIX(x134);

        /* Assignment result[1, 2]=mu*x132 + x134 - x51*x63 - x53*x76 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x48=mu*ut2; FIX(x48);
        x49=-x48*x39; FIX(x49);
        x50=x20*x42; FIX(x50);
        x51=x49 + x50; FIX(x51);
        x53=x49 - x50; FIX(x53);
        x61=rt1 - x38; FIX(x61);
        x62=0.5*x41*x25; FIX(x62);
        x63=x61*x62; FIX(x63);
        x65=0.5*x28*x41; FIX(x65);
        x76=x18*x65; FIX(x76);
        x109=Max(0, x27); FIX(x109);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x113=x24; FIX(x113);
        x130=x20*x61; FIX(x130);
        x131=0.5*x111*x113; FIX(x131);
        x132=x130*x131; FIX(x132);
        x133=0.5*x109*x18*x20*x111; FIX(x133);
        x134=mu*x133; FIX(x134);
        result[7] = mu*x132 + x134 - x51*x63 - x53*x76;

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");

        /* Assignment result[1, 2]=0 */

        result[7] = 0;
        /*@ assert (result[7]) >= 0.;*/
    }
    /*@ assert (result[7]) >= 0.;*/

    /* Assignment result[2, 2]=Piecewise((x148, x9), (-x89*(-1.17157287525381*x170*x74 - 6.82842712474619*x170*x75 - x120*x172 + x120*x174 + x160 + x166*x161 - x166*x35 - x167*x180 - x167*x181 - x168*x180 - x168*x181 + x177*x161 - x177*x35 + x178*x162 - x178*x163 - x178*x164 - x179*x105 + x179*x107 + x180*x175 + x181*x175 - x34*x166 - x34*x177 - x91*x179), x15), (-x110*(-mu*x182 + x42) - x114*(-mu*x183 + x115) - x51*x68 - x53*x155, x78), (x52 + x54, x80)) */
    double x148;
    double x177;
    double x178;
    double x179;
    double x180;
    double x181;
    double x182;
    double x183;
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x31=mu*x3; FIX(x31);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x148=x31*x85; FIX(x148);

        /* Assignment result[2, 2]=x148 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x31=mu*x3; FIX(x31);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x148=x31*x85; FIX(x148);
        result[8] = x148;

    }
    else if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x16=-un; FIX(x16);
        x34=2.0*x10; FIX(x34);
        x35=2.0*x13; FIX(x35);
        x55=1.0*mu; FIX(x55);
        x71=Heaviside(x7); FIX(x71);
        x73=-rn*x55 + un; FIX(x73);
        x74=Heaviside(x12 + x73); FIX(x74);
        x75=Heaviside(x7 + x73); FIX(x75);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        x89=0.25*mu*x88; FIX(x89);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x91=4.0*x90; FIX(x91);
        x96=1.4142135623730951454746218587388284504413604736328125*x74; FIX(x96);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x103=Max(0, -x7 + x97); FIX(x103);
        x105=2.0*x4*x98; FIX(x105);
        x107=2.0*x103*x4; FIX(x107);
        x108=x4*x90; FIX(x108);
        x120=pow(rt1, 4); FIX(x120);
        x121=pow(rt2, 3); FIX(x121);
        x139=x90*x5; FIX(x139);
        x157=1.4142135623730951454746218587388284504413604736328125*x10; FIX(x157);
        x158=1.4142135623730951454746218587388284504413604736328125*x13; FIX(x158);
        x159=1.4142135623730951454746218587388284504413604736328125*x71*x75; FIX(x159);
        x160=-x4*x90*x96*x71 - x108*x157 + x108*x158 + x108*x159 - x139*x157 + x139*x158; FIX(x160);
        x161=4.0*x71; FIX(x161);
        x162=8.0*x71; FIX(x162);
        x163=4.0*x10; FIX(x163);
        x164=4.0*x13; FIX(x164);
        x166=pow(rt2, 5); FIX(x166);
        x167=x166*x71; FIX(x167);
        x168=rt2*x120*x71; FIX(x168);
        x170=x4*x121*x71; FIX(x170);
        x172=2.0*x71*x98; FIX(x172);
        x174=2.0*x103*x71; FIX(x174);
        x175=x90*x5*x71; FIX(x175);
        x177=rt2*x120; FIX(x177);
        x178=x4*x121; FIX(x178);
        x179=x5*x71; FIX(x179);
        x180=3.41421356237309492343001693370752036571502685546875*x75; FIX(x180);
        x181=0.58578643762690496554768060377682559192180633544921875*x74; FIX(x181);

        /* Assignment result[2, 2]=-x89*(-1.17157287525381*x170*x74 - 6.82842712474619*x170*x75 - x120*x172 + x120*x174 + x160 + x166*x161 - x166*x35 - x167*x180 - x167*x181 - x168*x180 - x168*x181 + x177*x161 - x177*x35 + x178*x162 - x178*x163 - x178*x164 - x179*x105 + x179*x107 + x180*x175 + x181*x175 - x34*x166 - x34*x177 - x91*x179) */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x16=-un; FIX(x16);
        x34=2.0*x10; FIX(x34);
        x35=2.0*x13; FIX(x35);
        x55=1.0*mu; FIX(x55);
        x71=Heaviside(x7); FIX(x71);
        x73=-rn*x55 + un; FIX(x73);
        x74=Heaviside(x12 + x73); FIX(x74);
        x75=Heaviside(x7 + x73); FIX(x75);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        x89=0.25*mu*x88; FIX(x89);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x91=4.0*x90; FIX(x91);
        x96=1.4142135623730951454746218587388284504413604736328125*x74; FIX(x96);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x103=Max(0, -x7 + x97); FIX(x103);
        x105=2.0*x4*x98; FIX(x105);
        x107=2.0*x103*x4; FIX(x107);
        x108=x4*x90; FIX(x108);
        x120=pow(rt1, 4); FIX(x120);
        x121=pow(rt2, 3); FIX(x121);
        x139=x90*x5; FIX(x139);
        x157=1.4142135623730951454746218587388284504413604736328125*x10; FIX(x157);
        x158=1.4142135623730951454746218587388284504413604736328125*x13; FIX(x158);
        x159=1.4142135623730951454746218587388284504413604736328125*x71*x75; FIX(x159);
        x160=-x4*x90*x96*x71 - x108*x157 + x108*x158 + x108*x159 - x139*x157 + x139*x158; FIX(x160);
        x161=4.0*x71; FIX(x161);
        x162=8.0*x71; FIX(x162);
        x163=4.0*x10; FIX(x163);
        x164=4.0*x13; FIX(x164);
        x166=pow(rt2, 5); FIX(x166);
        x167=x166*x71; FIX(x167);
        x168=rt2*x120*x71; FIX(x168);
        x170=x4*x121*x71; FIX(x170);
        x172=2.0*x71*x98; FIX(x172);
        x174=2.0*x103*x71; FIX(x174);
        x175=x90*x5*x71; FIX(x175);
        x177=rt2*x120; FIX(x177);
        x178=x4*x121; FIX(x178);
        x179=x5*x71; FIX(x179);
        x180=3.41421356237309492343001693370752036571502685546875*x75; FIX(x180);
        x181=0.58578643762690496554768060377682559192180633544921875*x74; FIX(x181);
        result[8] = -x89*(-1.1715728752538099310953612075536511838436126708984375*x170*x74 - 6.8284271247461898468600338674150407314300537109375*x170*x75 - x120*x172 + x120*x174 + x160 + x166*x161 - x166*x35 - x167*x180 - x167*x181 - x168*x180 - x168*x181 + x177*x161 - x177*x35 + x178*x162 - x178*x163 - x178*x164 - x179*x105 + x179*x107 + x180*x175 + x181*x175 - x34*x166 - x34*x177 - x91*x179);

    }
    else if (x78)
    {
        DEBUG_PRINT("Case (x78) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x48=mu*ut2; FIX(x48);
        x49=-x48*x39; FIX(x49);
        x50=x20*x42; FIX(x50);
        x51=x49 + x50; FIX(x51);
        x53=x49 - x50; FIX(x53);
        x62=0.5*x41*x25; FIX(x62);
        x65=0.5*x28*x41; FIX(x65);
        x67=rt2 - x48; FIX(x67);
        x68=x67*x62; FIX(x68);
        x109=Max(0, x27); FIX(x109);
        x110=0.5*x109; FIX(x110);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x113=x24; FIX(x113);
        x114=0.5*x113; FIX(x114);
        x115=-x42; FIX(x115);
        x155=x20*x65; FIX(x155);
        x182=x21*x111; FIX(x182);
        x183=x20*x67*x111; FIX(x183);

        /* Assignment result[2, 2]=-x110*(-mu*x182 + x42) - x114*(-mu*x183 + x115) - x51*x68 - x53*x155 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x48=mu*ut2; FIX(x48);
        x49=-x48*x39; FIX(x49);
        x50=x20*x42; FIX(x50);
        x51=x49 + x50; FIX(x51);
        x53=x49 - x50; FIX(x53);
        x62=0.5*x41*x25; FIX(x62);
        x65=0.5*x28*x41; FIX(x65);
        x67=rt2 - x48; FIX(x67);
        x68=x67*x62; FIX(x68);
        x109=Max(0, x27); FIX(x109);
        x110=0.5*x109; FIX(x110);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x113=x24; FIX(x113);
        x114=0.5*x113; FIX(x114);
        x115=-x42; FIX(x115);
        x155=x20*x65; FIX(x155);
        x182=x21*x111; FIX(x182);
        x183=x20*x67*x111; FIX(x183);
        result[8] = -x110*(-mu*x182 + x42) - x114*(-mu*x183 + x115) - x51*x68 - x53*x155;

    }
    else if (x80)
    {
        DEBUG_PRINT("Case (x80) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x48=mu*ut2; FIX(x48);
        x49=-x48*x39; FIX(x49);
        x50=x20*x42; FIX(x50);
        x51=x49 + x50; FIX(x51);
        x52=-x26*x51; FIX(x52);
        x53=x49 - x50; FIX(x53);
        x54=x53*x29; FIX(x54);

        /* Assignment result[2, 2]=x52 + x54 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        /*@ assert (x8) != 0.;*/
        x39=1.0/x8; FIX(x39);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x48=mu*ut2; FIX(x48);
        x49=-x48*x39; FIX(x49);
        x50=x20*x42; FIX(x50);
        x51=x49 + x50; FIX(x51);
        x52=-x26*x51; FIX(x52);
        x53=x49 - x50; FIX(x53);
        x54=x53*x29; FIX(x54);
        result[8] = x52 + x54;

    }
    /*@ assert (result[8]) >= 0.;*/

    /* Assignment result[0, 3]=Piecewise((mu - x3*x55, x9), (mu - mu*x11 - mu*x14, x56), (mu + x57 - x58, x59)) */

    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x55=1.0*mu; FIX(x55);

        /* Assignment result[0, 3]=mu - x3*x55 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x55=1.0*mu; FIX(x55);
        result[9] = mu - x3*x55;

    }
    else if (x56)
    {
        DEBUG_PRINT("Case (x56) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x11=0.5*x10; FIX(x11);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x14=0.5*x13; FIX(x14);

        /* Assignment result[0, 3]=mu - mu*x11 - mu*x14 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x11=0.5*x10; FIX(x11);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x14=0.5*x13; FIX(x14);
        result[9] = mu - mu*x11 - mu*x14;

    }
    else if (x59)
    {
        DEBUG_PRINT("Case (x59) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x57=-mu*x26; FIX(x57);
        x58=mu*x29; FIX(x58);

        /* Assignment result[0, 3]=mu + x57 - x58 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x57=-mu*x26; FIX(x57);
        x58=mu*x29; FIX(x58);
        result[9] = mu + x57 - x58;

    }
    /*@ assert (result[9]) >= 0.;*/

    /* Assignment result[1, 3]=Piecewise((0.0, x9), (-rt1*x135, x56), (-x43*x29 - x61*x136, x137), (0, x138)) */
    double x135;
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");

        /* Assignment result[1, 3]=0.0 */

        result[10] = 0.0;
        /*@ assert (result[10]) >= 0.;*/
    }
    else if (x56)
    {
        DEBUG_PRINT("Case (x56) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x135=0.5*mu*x32*(x10 - 1.0*x13); FIX(x135);

        /* Assignment result[1, 3]=-rt1*x135 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x135=0.5*mu*x32*(x10 - 1.0*x13); FIX(x135);
        result[10] = -rt1*x135;

    }
    else if (x137)
    {
        DEBUG_PRINT("Case (x137) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x43=x18*x42; FIX(x43);
        x61=rt1 - x38; FIX(x61);
        x136=0.5*mu*x41*x25; FIX(x136);

        /* Assignment result[1, 3]=-x43*x29 - x61*x136 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x43=x18*x42; FIX(x43);
        x61=rt1 - x38; FIX(x61);
        x136=0.5*mu*x41*x25; FIX(x136);
        result[10] = -x43*x29 - x61*x136;

    }
    else if (x138)
    {
        DEBUG_PRINT("Case (x138) is True.\n");

        /* Assignment result[1, 3]=0 */

        result[10] = 0;
        /*@ assert (result[10]) >= 0.;*/
    }
    /*@ assert (result[10]) >= 0.;*/

    /* Assignment result[2, 3]=Piecewise((0.0, x9), (-rt2*x135, x56), (-x50*x29 - x67*x136, x137), (x57 + x58, x138)) */

    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");

        /* Assignment result[2, 3]=0.0 */

        result[11] = 0.0;
        /*@ assert (result[11]) >= 0.;*/
    }
    else if (x56)
    {
        DEBUG_PRINT("Case (x56) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x135=0.5*mu*x32*(x10 - 1.0*x13); FIX(x135);

        /* Assignment result[2, 3]=-rt2*x135 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x135=0.5*mu*x32*(x10 - 1.0*x13); FIX(x135);
        result[11] = -rt2*x135;

    }
    else if (x137)
    {
        DEBUG_PRINT("Case (x137) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x48=mu*ut2; FIX(x48);
        x50=x20*x42; FIX(x50);
        x67=rt2 - x48; FIX(x67);
        x136=0.5*mu*x41*x25; FIX(x136);

        /* Assignment result[2, 3]=-x50*x29 - x67*x136 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x42=mu*x41; FIX(x42);
        x48=mu*ut2; FIX(x48);
        x50=x20*x42; FIX(x50);
        x67=rt2 - x48; FIX(x67);
        x136=0.5*mu*x41*x25; FIX(x136);
        result[11] = -x50*x29 - x67*x136;

    }
    else if (x138)
    {
        DEBUG_PRINT("Case (x138) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x57=-mu*x26; FIX(x57);
        x58=mu*x29; FIX(x58);

        /* Assignment result[2, 3]=x57 + x58 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x26=0.5*x25; FIX(x26);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x29=0.5*x28; FIX(x29);
        x57=-mu*x26; FIX(x57);
        x58=mu*x29; FIX(x58);
        result[11] = x57 + x58;

    }
    /*@ assert (result[11]) >= 0.;*/

    /* Assignment result[0, 4]=Piecewise((0.0, x9), (rt1*x60, x56), (x64 + x66, x59)) */
    double x60;
    double x64;
    double x66;
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");

        /* Assignment result[0, 4]=0.0 */

        result[12] = 0.0;
        /*@ assert (result[12]) >= 0.;*/
    }
    else if (x56)
    {
        DEBUG_PRINT("Case (x56) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x60=0.5*x32*(-1.0*x10 + x13); FIX(x60);

        /* Assignment result[0, 4]=rt1*x60 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x60=0.5*x32*(-1.0*x10 + x13); FIX(x60);
        result[12] = rt1*x60;

    }
    else if (x59)
    {
        DEBUG_PRINT("Case (x59) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x61=rt1 - x38; FIX(x61);
        x62=0.5*x41*x25; FIX(x62);
        x63=x61*x62; FIX(x63);
        x64=-x63; FIX(x64);
        x65=0.5*x28*x41; FIX(x65);
        x66=x61*x65; FIX(x66);

        /* Assignment result[0, 4]=x64 + x66 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x61=rt1 - x38; FIX(x61);
        x62=0.5*x41*x25; FIX(x62);
        x63=x61*x62; FIX(x63);
        x64=-x63; FIX(x64);
        x65=0.5*x28*x41; FIX(x65);
        x66=x61*x65; FIX(x66);
        result[12] = x64 + x66;

    }
    /*@ assert (result[12]) >= 0.;*/

    /* Assignment result[1, 4]=Piecewise((1.0 + x118, x9), (x88*(-x108*x11 - x108*x14 + x142 - x94*x140 + x94*x143), x56), (x18*x61*x146 + 1 - x110*(x112 + x147) - x114*(x116 + x41) - x145*x61**2, x137), (1, x138)) */
    double x140;
    double x141;
    double x142;
    double x143;
    double x144;
    double x145;
    double x146;
    double x147;
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x117=1.0*x3*x85; FIX(x117);
        x118=-x81*x117; FIX(x118);

        /* Assignment result[1, 4]=1.0 + x118 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x117=1.0*x3*x85; FIX(x117);
        x118=-x81*x117; FIX(x118);
        result[13] = 1.0 + x118;

    }
    else if (x56)
    {
        DEBUG_PRINT("Case (x56) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x11=0.5*x10; FIX(x11);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x14=0.5*x13; FIX(x14);
        x16=-un; FIX(x16);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x94=pow(rt2, 4); FIX(x94);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x103=Max(0, -x7 + x97); FIX(x103);
        x108=x4*x90; FIX(x108);
        x123=x4*x5; FIX(x123);
        x139=x90*x5; FIX(x139);
        x140=0.5*x98; FIX(x140);
        x141=x103*x4; FIX(x141);
        x142=0.5*x5*x141 + x108 - x123*x140 + x139; FIX(x142);
        x143=0.5*x103; FIX(x143);

        /* Assignment result[1, 4]=x88*(-x108*x11 - x108*x14 + x142 - x94*x140 + x94*x143) */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x11=0.5*x10; FIX(x11);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x14=0.5*x13; FIX(x14);
        x16=-un; FIX(x16);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x94=pow(rt2, 4); FIX(x94);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x103=Max(0, -x7 + x97); FIX(x103);
        x108=x4*x90; FIX(x108);
        x123=x4*x5; FIX(x123);
        x139=x90*x5; FIX(x139);
        x140=0.5*x98; FIX(x140);
        x141=x103*x4; FIX(x141);
        x142=0.5*x5*x141 + x108 - x123*x140 + x139; FIX(x142);
        x143=0.5*x103; FIX(x143);
        result[13] = x88*(-x108*x11 - x108*x14 + x142 - x94*x140 + x94*x143);

    }
    else if (x137)
    {
        DEBUG_PRINT("Case (x137) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x61=rt1 - x38; FIX(x61);
        x109=Max(0, x27); FIX(x109);
        x110=0.5*x109; FIX(x110);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x112=x19*x111; FIX(x112);
        x113=x24; FIX(x113);
        x114=0.5*x113; FIX(x114);
        x116=x18*x61*x111; FIX(x116);
        /*@ assert (x22) != 0.;*/
        x144=1.0/x22; FIX(x144);
        x145=0.5*x144*x25; FIX(x145);
        x146=0.5*x144*x28; FIX(x146);
        x147=-x41; FIX(x147);

        /* Assignment result[1, 4]=x18*x61*x146 + 1 - x110*(x112 + x147) - x114*(x116 + x41) - x145*x61**2 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x61=rt1 - x38; FIX(x61);
        x109=Max(0, x27); FIX(x109);
        x110=0.5*x109; FIX(x110);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x112=x19*x111; FIX(x112);
        x113=x24; FIX(x113);
        x114=0.5*x113; FIX(x114);
        x116=x18*x61*x111; FIX(x116);
        /*@ assert (x22) != 0.;*/
        x144=1.0/x22; FIX(x144);
        x145=0.5*x144*x25; FIX(x145);
        x146=0.5*x144*x28; FIX(x146);
        x147=-x41; FIX(x147);
        result[13] = x18*x61*x146 + 1 - x110*(x112 + x147) - x114*(x116 + x41) - x145*x61*x61;

    }
    else if (x138)
    {
        DEBUG_PRINT("Case (x138) is True.\n");

        /* Assignment result[1, 4]=1 */

        result[13] = 1;
        /*@ assert (result[13]) >= 0.;*/
        /*@ assert (result[13]) != 0.;*/
    }
    /*@ assert (result[13]) >= 0.;*/

    /* Assignment result[2, 4]=Piecewise((x148, x9), (x150, x56), (x130*x146 + x151 - x176, x137), (x64 - x66, x138)) */
    double x149;
    double x150;
    double x151;
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x31=mu*x3; FIX(x31);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x148=x31*x85; FIX(x148);

        /* Assignment result[2, 4]=x148 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x31=mu*x3; FIX(x31);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x148=x31*x85; FIX(x148);
        result[14] = x148;

    }
    else if (x56)
    {
        DEBUG_PRINT("Case (x56) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x16=-un; FIX(x16);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x103=Max(0, -x7 + x97); FIX(x103);
        x141=x103*x4; FIX(x141);
        x149=1.0*x98; FIX(x149);
        x150=-0.5*rt1*rt2*x88*(x90*x10 + x103*x5 + x90*x13 - x4*x149 - x5*x149 + x141); FIX(x150);

        /* Assignment result[2, 4]=x150 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x16=-un; FIX(x16);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x103=Max(0, -x7 + x97); FIX(x103);
        x141=x103*x4; FIX(x141);
        x149=1.0*x98; FIX(x149);
        x150=-0.5*rt1*rt2*x88*(x90*x10 + x103*x5 + x90*x13 - x4*x149 - x5*x149 + x141); FIX(x150);
        result[14] = x150;

    }
    else if (x137)
    {
        DEBUG_PRINT("Case (x137) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        x48=mu*ut2; FIX(x48);
        x61=rt1 - x38; FIX(x61);
        x67=rt2 - x48; FIX(x67);
        x109=Max(0, x27); FIX(x109);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x113=x24; FIX(x113);
        x130=x20*x61; FIX(x130);
        x131=0.5*x111*x113; FIX(x131);
        x133=0.5*x109*x18*x20*x111; FIX(x133);
        /*@ assert (x22) != 0.;*/
        x144=1.0/x22; FIX(x144);
        x145=0.5*x144*x25; FIX(x145);
        x146=0.5*x144*x28; FIX(x146);
        x151=-x133 - x145*x61*x67; FIX(x151);
        x152=x18*x67; FIX(x152);
        x176=x152*x131; FIX(x176);

        /* Assignment result[2, 4]=x130*x146 + x151 - x176 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        x48=mu*ut2; FIX(x48);
        x61=rt1 - x38; FIX(x61);
        x67=rt2 - x48; FIX(x67);
        x109=Max(0, x27); FIX(x109);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x113=x24; FIX(x113);
        x130=x20*x61; FIX(x130);
        x131=0.5*x111*x113; FIX(x131);
        x133=0.5*x109*x18*x20*x111; FIX(x133);
        /*@ assert (x22) != 0.;*/
        x144=1.0/x22; FIX(x144);
        x145=0.5*x144*x25; FIX(x145);
        x146=0.5*x144*x28; FIX(x146);
        x151=-x133 - x145*x61*x67; FIX(x151);
        x152=x18*x67; FIX(x152);
        x176=x152*x131; FIX(x176);
        result[14] = x130*x146 + x151 - x176;

    }
    else if (x138)
    {
        DEBUG_PRINT("Case (x138) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x61=rt1 - x38; FIX(x61);
        x62=0.5*x41*x25; FIX(x62);
        x63=x61*x62; FIX(x63);
        x64=-x63; FIX(x64);
        x65=0.5*x28*x41; FIX(x65);
        x66=x61*x65; FIX(x66);

        /* Assignment result[2, 4]=x64 - x66 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x61=rt1 - x38; FIX(x61);
        x62=0.5*x41*x25; FIX(x62);
        x63=x61*x62; FIX(x63);
        x64=-x63; FIX(x64);
        x65=0.5*x28*x41; FIX(x65);
        x66=x61*x65; FIX(x66);
        result[14] = x64 - x66;

    }
    /*@ assert (result[14]) >= 0.;*/

    /* Assignment result[0, 5]=Piecewise((0.0, x9), (rt2*x60, x56), (x69 + x70, x59)) */
    double x69;
    double x70;
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");

        /* Assignment result[0, 5]=0.0 */

        result[15] = 0.0;
        /*@ assert (result[15]) >= 0.;*/
    }
    else if (x56)
    {
        DEBUG_PRINT("Case (x56) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x60=0.5*x32*(-1.0*x10 + x13); FIX(x60);

        /* Assignment result[0, 5]=rt2*x60 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        /*@ assert (x7) != 0.;*/
        x32=1.0/x7; FIX(x32);
        x60=0.5*x32*(-1.0*x10 + x13); FIX(x60);
        result[15] = rt2*x60;

    }
    else if (x59)
    {
        DEBUG_PRINT("Case (x59) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x48=mu*ut2; FIX(x48);
        x62=0.5*x41*x25; FIX(x62);
        x65=0.5*x28*x41; FIX(x65);
        x67=rt2 - x48; FIX(x67);
        x68=x67*x62; FIX(x68);
        x69=-x68; FIX(x69);
        x70=x67*x65; FIX(x70);

        /* Assignment result[0, 5]=x69 + x70 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x48=mu*ut2; FIX(x48);
        x62=0.5*x41*x25; FIX(x62);
        x65=0.5*x28*x41; FIX(x65);
        x67=rt2 - x48; FIX(x67);
        x68=x67*x62; FIX(x68);
        x69=-x68; FIX(x69);
        x70=x67*x65; FIX(x70);
        result[15] = x69 + x70;

    }
    /*@ assert (result[15]) >= 0.;*/

    /* Assignment result[1, 5]=Piecewise((x148, x9), (x150, x56), (-x132 + x151 + x152*x146, x137), (0, x138)) */

    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x31=mu*x3; FIX(x31);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x148=x31*x85; FIX(x148);

        /* Assignment result[1, 5]=x148 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x31=mu*x3; FIX(x31);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x148=x31*x85; FIX(x148);
        result[16] = x148;

    }
    else if (x56)
    {
        DEBUG_PRINT("Case (x56) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x16=-un; FIX(x16);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x103=Max(0, -x7 + x97); FIX(x103);
        x141=x103*x4; FIX(x141);
        x149=1.0*x98; FIX(x149);
        x150=-0.5*rt1*rt2*x88*(x90*x10 + x103*x5 + x90*x13 - x4*x149 - x5*x149 + x141); FIX(x150);

        /* Assignment result[1, 5]=x150 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x16=-un; FIX(x16);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x103=Max(0, -x7 + x97); FIX(x103);
        x141=x103*x4; FIX(x141);
        x149=1.0*x98; FIX(x149);
        x150=-0.5*rt1*rt2*x88*(x90*x10 + x103*x5 + x90*x13 - x4*x149 - x5*x149 + x141); FIX(x150);
        result[16] = x150;

    }
    else if (x137)
    {
        DEBUG_PRINT("Case (x137) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        x48=mu*ut2; FIX(x48);
        x61=rt1 - x38; FIX(x61);
        x67=rt2 - x48; FIX(x67);
        x109=Max(0, x27); FIX(x109);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x113=x24; FIX(x113);
        x130=x20*x61; FIX(x130);
        x131=0.5*x111*x113; FIX(x131);
        x132=x130*x131; FIX(x132);
        x133=0.5*x109*x18*x20*x111; FIX(x133);
        /*@ assert (x22) != 0.;*/
        x144=1.0/x22; FIX(x144);
        x145=0.5*x144*x25; FIX(x145);
        x146=0.5*x144*x28; FIX(x146);
        x151=-x133 - x145*x61*x67; FIX(x151);
        x152=x18*x67; FIX(x152);

        /* Assignment result[1, 5]=-x132 + x151 + x152*x146 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        x38=mu*ut1; FIX(x38);
        x48=mu*ut2; FIX(x48);
        x61=rt1 - x38; FIX(x61);
        x67=rt2 - x48; FIX(x67);
        x109=Max(0, x27); FIX(x109);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x113=x24; FIX(x113);
        x130=x20*x61; FIX(x130);
        x131=0.5*x111*x113; FIX(x131);
        x132=x130*x131; FIX(x132);
        x133=0.5*x109*x18*x20*x111; FIX(x133);
        /*@ assert (x22) != 0.;*/
        x144=1.0/x22; FIX(x144);
        x145=0.5*x144*x25; FIX(x145);
        x146=0.5*x144*x28; FIX(x146);
        x151=-x133 - x145*x61*x67; FIX(x151);
        x152=x18*x67; FIX(x152);
        result[16] = -x132 + x151 + x152*x146;

    }
    else if (x138)
    {
        DEBUG_PRINT("Case (x138) is True.\n");

        /* Assignment result[1, 5]=0 */

        result[16] = 0;
        /*@ assert (result[16]) >= 0.;*/
    }
    /*@ assert (result[16]) >= 0.;*/

    /* Assignment result[2, 5]=Piecewise((1.0 - x117, x9), (x88*(-x120*x140 + x120*x143 - x139*x11 - x139*x14 + x142), x56), (x20*x67*x146 + 1 - x110*(x147 + x182) - x114*(x183 + x41) - x145*x67**2, x137), (1 + x69 - x70, x138)) */

    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x117=1.0*x3*x85; FIX(x117);

        /* Assignment result[2, 5]=1.0 - x117 */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x3=Heaviside(x2); FIX(x3);
        x81=mu*mu; FIX(x81);
        x82=1.0 + x81; FIX(x82);
        /*@ assert (x82) != 0.;*/
        x85=1.0/x82; FIX(x85);
        x117=1.0*x3*x85; FIX(x117);
        result[17] = 1.0 - x117;

    }
    else if (x56)
    {
        DEBUG_PRINT("Case (x56) is True.\n");
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x11=0.5*x10; FIX(x11);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x14=0.5*x13; FIX(x14);
        x16=-un; FIX(x16);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x103=Max(0, -x7 + x97); FIX(x103);
        x108=x4*x90; FIX(x108);
        x120=pow(rt1, 4); FIX(x120);
        x123=x4*x5; FIX(x123);
        x139=x90*x5; FIX(x139);
        x140=0.5*x98; FIX(x140);
        x141=x103*x4; FIX(x141);
        x142=0.5*x5*x141 + x108 - x123*x140 + x139; FIX(x142);
        x143=0.5*x103; FIX(x143);

        /* Assignment result[2, 5]=x88*(-x120*x140 + x120*x143 - x139*x11 - x139*x14 + x142) */
        x1=mu*rn; FIX(x1);
        x2=-1.0*un + x1; FIX(x2);
        x10=Heaviside(x2 + x7); FIX(x10);
        x11=0.5*x10; FIX(x11);
        x12=-1.0*x7; FIX(x12);
        x13=Heaviside(x12 + x2); FIX(x13);
        x14=0.5*x13; FIX(x14);
        x16=-un; FIX(x16);
        /*@ assert (x6) >= 0.;*/
        /*@ assert (x6) != 0.;*/
        x88=pow(x6, -5.0/2.0); FIX(x88);
        /*@ assert (x6) >= 0.;*/
        x90=pow(x6, 3.0/2.0); FIX(x90);
        x97=x1 + x16; FIX(x97);
        x98=Max(0, x7 + x97); FIX(x98);
        x103=Max(0, -x7 + x97); FIX(x103);
        x108=x4*x90; FIX(x108);
        x120=pow(rt1, 4); FIX(x120);
        x123=x4*x5; FIX(x123);
        x139=x90*x5; FIX(x139);
        x140=0.5*x98; FIX(x140);
        x141=x103*x4; FIX(x141);
        x142=0.5*x5*x141 + x108 - x123*x140 + x139; FIX(x142);
        x143=0.5*x103; FIX(x143);
        result[17] = x88*(-x120*x140 + x120*x143 - x139*x11 - x139*x14 + x142);

    }
    else if (x137)
    {
        DEBUG_PRINT("Case (x137) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x48=mu*ut2; FIX(x48);
        x67=rt2 - x48; FIX(x67);
        x109=Max(0, x27); FIX(x109);
        x110=0.5*x109; FIX(x110);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x113=x24; FIX(x113);
        x114=0.5*x113; FIX(x114);
        /*@ assert (x22) != 0.;*/
        x144=1.0/x22; FIX(x144);
        x145=0.5*x144*x25; FIX(x145);
        x146=0.5*x144*x28; FIX(x146);
        x147=-x41; FIX(x147);
        x182=x21*x111; FIX(x182);
        x183=x20*x67*x111; FIX(x183);

        /* Assignment result[2, 5]=x20*x67*x146 + 1 - x110*(x147 + x182) - x114*(x183 + x41) - x145*x67**2 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x48=mu*ut2; FIX(x48);
        x67=rt2 - x48; FIX(x67);
        x109=Max(0, x27); FIX(x109);
        x110=0.5*x109; FIX(x110);
        /*@ assert (x22) >= 0.;*/
        /*@ assert (x22) != 0.;*/
        x111=pow(x22, -3.0/2.0); FIX(x111);
        x113=x24; FIX(x113);
        x114=0.5*x113; FIX(x114);
        /*@ assert (x22) != 0.;*/
        x144=1.0/x22; FIX(x144);
        x145=0.5*x144*x25; FIX(x145);
        x146=0.5*x144*x28; FIX(x146);
        x147=-x41; FIX(x147);
        x182=x21*x111; FIX(x182);
        x183=x20*x67*x111; FIX(x183);
        result[17] = x20*x67*x146 + 1 - x110*(x147 + x182) - x114*(x183 + x41) - x145*x67*x67;

    }
    else if (x138)
    {
        DEBUG_PRINT("Case (x138) is True.\n");
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x48=mu*ut2; FIX(x48);
        x62=0.5*x41*x25; FIX(x62);
        x65=0.5*x28*x41; FIX(x65);
        x67=rt2 - x48; FIX(x67);
        x68=x67*x62; FIX(x68);
        x69=-x68; FIX(x69);
        x70=x67*x65; FIX(x70);

        /* Assignment result[2, 5]=1 + x69 - x70 */
        x16=-un; FIX(x16);
        x17=mu*rn - mu*x8 + x16; FIX(x17);
        x24=x17 + x23; FIX(x24);
        x25=Heaviside(x24); FIX(x25);
        x27=x17 - x23; FIX(x27);
        x28=Heaviside(x27); FIX(x28);
        /*@ assert (x23) != 0.;*/
        x41=1.0/x23; FIX(x41);
        x48=mu*ut2; FIX(x48);
        x62=0.5*x41*x25; FIX(x62);
        x65=0.5*x28*x41; FIX(x65);
        x67=rt2 - x48; FIX(x67);
        x68=x67*x62; FIX(x68);
        x69=-x68; FIX(x69);
        x70=x67*x65; FIX(x70);
        result[17] = 1 + x69 - x70;

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
