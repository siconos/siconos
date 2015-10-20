#include <math.h>
#include <assert.h>
#include <op3x3.h>
#include <stdlib.h>

//#define DEBUG_MESSAGES 1
//#include <stdio.h>
#include <debug.h>
#include "FischerBurmeisterGenerated.h"

#pragma GCC diagnostic ignored "-Wmissing-prototypes"
#ifdef __clang__
#pragma clang diagnostic ignored "-Wlogical-op-parentheses"
#endif

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

// this file consists of generated code
//#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

void fc3d_AlartCurnierFABGenerated(
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
    int x3;
    x1=rhon*un;
    x2=rn - x1;
    x3=x2 <= 0;
    int x4;
    x4=x2 > 0;
    double x5;
    double x6;
    double x7;
    double x8;
    double x9;
    double x10;
    int x11;
    double x12;
    int x13;
    int x14;
    x5=rhot1*ut1;
    x6=-rt1 + x5;
    x7=rhot2*ut2;
    x8=-rt2 + x7;
    x9=x6*x6 + x8*x8;
    x10=(assert(IS_POSITIVE(x9)), sqrt(x9));
    x11=x10 <= 0;
    x12=mu*x2;
    x13=x10 <= x12;
    x14=x11 && x3 || x13 && x4;
    int x15;
    int x16;
    int x17;
    int x18;
    x15=x10 > 0;
    x16=x10 > x12;
    x17=x15 && x3 || x16 && x4;
    x18=x17 && x3;
    int x22;
    double x19;
    double x20;
    double x21;
    x22=x17 && x4;
    int x23;
    x23=(x13 || x3) && (x11 || x13 || x3) && (x11 || x3 || x4) && (x13 || x16 || x3) && (x13 || x3 || x4) && (x11 || x13 || x16 || x3) && (x11 || x13 || x3 || x4) && (x13 || x15 || x16 || x3) && (x13 || x15 || x3 || x4);
    /* Assignment result[0, 0]=Piecewise((rn, x3), (x1, x4)) */
    if (x3)
    {
        DEBUG_PRINT("Case (x3) is True.\n");

        /* Assignment result[0, 0]=rn */
        result[0] = rn;
    }
    else if (x4)
    {
        DEBUG_PRINT("Case (x4) is True.\n");

        /* Assignment result[0, 0]=x1 */
        result[0] = x1;
    }

    /* Assignment result[1, 0]=Piecewise((x5, x14), (rt1, x18), (rt1 - x19*x21, x22)) */
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");

        /* Assignment result[1, 0]=x5 */
        result[1] = x5;
    }
    else if (x18)
    {
        DEBUG_PRINT("Case (x18) is True.\n");

        /* Assignment result[1, 0]=rt1 */
        result[1] = rt1;
    }
    else if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");
        x19=rt1 - x5;
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x21=x12*x20;

        /* Assignment result[1, 0]=rt1 - x19*x21 */
        result[1] = rt1 - x19*x21;
    }

    /* Assignment result[2, 0]=Piecewise((x7, x14), (rt2, x18), (rt2 - x21*x31, x22)) */
    double x31;if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");

        /* Assignment result[2, 0]=x7 */
        result[2] = x7;
    }
    else if (x18)
    {
        DEBUG_PRINT("Case (x18) is True.\n");

        /* Assignment result[2, 0]=rt2 */
        result[2] = rt2;
    }
    else if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x21=x12*x20;
        x31=rt2 - x7;

        /* Assignment result[2, 0]=rt2 - x21*x31 */
        result[2] = rt2 - x21*x31;
    }

    /* Assignment result[0, 1]=Piecewise((0, x3), (rhon, x4)) */
    if (x3)
    {
        DEBUG_PRINT("Case (x3) is True.\n");

        /* Assignment result[0, 1]=0 */
        result[3] = 0;
    }
    else if (x4)
    {
        DEBUG_PRINT("Case (x4) is True.\n");

        /* Assignment result[0, 1]=rhon */
        result[3] = rhon;
    }

    /* Assignment result[1, 1]=Piecewise((0, x23), (rhon*x25, x22)) */
    double x24;
    double x25;if (x23)
    {
        DEBUG_PRINT("Case (x23) is True.\n");

        /* Assignment result[1, 1]=0 */
        result[4] = 0;
    }
    else if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");
        x19=rt1 - x5;
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x24=mu*x20;
        x25=x19*x24;

        /* Assignment result[1, 1]=rhon*x25 */
        result[4] = rhon*x25;
    }

    /* Assignment result[2, 1]=Piecewise((0, x23), (rhon*x32, x22)) */
    double x32;if (x23)
    {
        DEBUG_PRINT("Case (x23) is True.\n");

        /* Assignment result[2, 1]=0 */
        result[5] = 0;
    }
    else if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x24=mu*x20;
        x31=rt2 - x7;
        x32=x24*x31;

        /* Assignment result[2, 1]=rhon*x32 */
        result[5] = rhon*x32;
    }

    /* Assignment result[0, 2]=0 */
    result[6] = 0;

    /* Assignment result[1, 2]=Piecewise((rhot1, x14), (0, x18), (rhot1*x21 + rhot1*x28, x22)) */
    double x26;
    double x27;
    double x28;if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");

        /* Assignment result[1, 2]=rhot1 */
        result[7] = rhot1;
    }
    else if (x18)
    {
        DEBUG_PRINT("Case (x18) is True.\n");

        /* Assignment result[1, 2]=0 */
        result[7] = 0;
    }
    else if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");
        x19=rt1 - x5;
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x21=x12*x20;
        x26=pow(x9, -3.0/2.0);
        x27=mu*x19*x2*x26;
        x28=x27*x6;

        /* Assignment result[1, 2]=rhot1*x21 + rhot1*x28 */
        result[7] = rhot1*x21 + rhot1*x28;
    }

    /* Assignment result[2, 2]=Piecewise((0, x23), (rhot1*x34, x22)) */
    double x33;
    double x34;if (x23)
    {
        DEBUG_PRINT("Case (x23) is True.\n");

        /* Assignment result[2, 2]=0 */
        result[8] = 0;
    }
    else if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");
        x26=pow(x9, -3.0/2.0);
        x31=rt2 - x7;
        x33=mu*x2*x26*x31;
        x34=x33*x6;

        /* Assignment result[2, 2]=rhot1*x34 */
        result[8] = rhot1*x34;
    }

    /* Assignment result[0, 3]=0 */
    result[9] = 0;

    /* Assignment result[1, 3]=Piecewise((0, x23), (rhot2*x29, x22)) */
    double x29;if (x23)
    {
        DEBUG_PRINT("Case (x23) is True.\n");

        /* Assignment result[1, 3]=0 */
        result[10] = 0;
    }
    else if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");
        x19=rt1 - x5;
        x26=pow(x9, -3.0/2.0);
        x27=mu*x19*x2*x26;
        x29=x27*x8;

        /* Assignment result[1, 3]=rhot2*x29 */
        result[10] = rhot2*x29;
    }

    /* Assignment result[2, 3]=Piecewise((rhot2, x14), (0, x18), (rhot2*x21 + rhot2*x35, x22)) */
    double x35;if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");

        /* Assignment result[2, 3]=rhot2 */
        result[11] = rhot2;
    }
    else if (x18)
    {
        DEBUG_PRINT("Case (x18) is True.\n");

        /* Assignment result[2, 3]=0 */
        result[11] = 0;
    }
    else if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x21=x12*x20;
        x26=pow(x9, -3.0/2.0);
        x31=rt2 - x7;
        x33=mu*x2*x26*x31;
        x35=x33*x8;

        /* Assignment result[2, 3]=rhot2*x21 + rhot2*x35 */
        result[11] = rhot2*x21 + rhot2*x35;
    }

    /* Assignment result[0, 4]=Piecewise((1, x3), (0, x4)) */
    if (x3)
    {
        DEBUG_PRINT("Case (x3) is True.\n");

        /* Assignment result[0, 4]=1 */
        result[12] = 1;
    }
    else if (x4)
    {
        DEBUG_PRINT("Case (x4) is True.\n");

        /* Assignment result[0, 4]=0 */
        result[12] = 0;
    }

    /* Assignment result[1, 4]=Piecewise((0, x23), (-x25, x22)) */
    if (x23)
    {
        DEBUG_PRINT("Case (x23) is True.\n");

        /* Assignment result[1, 4]=0 */
        result[13] = 0;
    }
    else if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");
        x19=rt1 - x5;
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x24=mu*x20;
        x25=x19*x24;

        /* Assignment result[1, 4]=-x25 */
        result[13] = -x25;
    }

    /* Assignment result[2, 4]=Piecewise((0, x23), (-x32, x22)) */
    if (x23)
    {
        DEBUG_PRINT("Case (x23) is True.\n");

        /* Assignment result[2, 4]=0 */
        result[14] = 0;
    }
    else if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x24=mu*x20;
        x31=rt2 - x7;
        x32=x24*x31;

        /* Assignment result[2, 4]=-x32 */
        result[14] = -x32;
    }

    /* Assignment result[0, 5]=0 */
    result[15] = 0;

    /* Assignment result[1, 5]=Piecewise((0, x14), (1, x18), (-x28 + x30, x22)) */
    double x30;if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");

        /* Assignment result[1, 5]=0 */
        result[16] = 0;
    }
    else if (x18)
    {
        DEBUG_PRINT("Case (x18) is True.\n");

        /* Assignment result[1, 5]=1 */
        result[16] = 1;
    }
    else if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");
        x19=rt1 - x5;
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x21=x12*x20;
        x26=pow(x9, -3.0/2.0);
        x27=mu*x19*x2*x26;
        x28=x27*x6;
        x30=-x21 + 1;

        /* Assignment result[1, 5]=-x28 + x30 */
        result[16] = -x28 + x30;
    }

    /* Assignment result[2, 5]=Piecewise((0, x23), (-x34, x22)) */
    if (x23)
    {
        DEBUG_PRINT("Case (x23) is True.\n");

        /* Assignment result[2, 5]=0 */
        result[17] = 0;
    }
    else if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");
        x26=pow(x9, -3.0/2.0);
        x31=rt2 - x7;
        x33=mu*x2*x26*x31;
        x34=x33*x6;

        /* Assignment result[2, 5]=-x34 */
        result[17] = -x34;
    }

    /* Assignment result[0, 6]=0 */
    result[18] = 0;

    /* Assignment result[1, 6]=Piecewise((0, x23), (-x29, x22)) */
    if (x23)
    {
        DEBUG_PRINT("Case (x23) is True.\n");

        /* Assignment result[1, 6]=0 */
        result[19] = 0;
    }
    else if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");
        x19=rt1 - x5;
        x26=pow(x9, -3.0/2.0);
        x27=mu*x19*x2*x26;
        x29=x27*x8;

        /* Assignment result[1, 6]=-x29 */
        result[19] = -x29;
    }

    /* Assignment result[2, 6]=Piecewise((0, x14), (1, x18), (x30 - x35, x22)) */
    if (x14)
    {
        DEBUG_PRINT("Case (x14) is True.\n");

        /* Assignment result[2, 6]=0 */
        result[20] = 0;
    }
    else if (x18)
    {
        DEBUG_PRINT("Case (x18) is True.\n");

        /* Assignment result[2, 6]=1 */
        result[20] = 1;
    }
    else if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x21=x12*x20;
        x26=pow(x9, -3.0/2.0);
        x30=-x21 + 1;
        x31=rt2 - x7;
        x33=mu*x2*x26*x31;
        x35=x33*x8;

        /* Assignment result[2, 6]=x30 - x35 */
        result[20] = x30 - x35;
    }
}
void fc3d_AlartCurnierFGenerated(
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
    int x3;
    x1=rhon*un;
    x2=rn - x1;
    x3=x2 <= 0;
    int x4;
    x4=x2 > 0;
    double x5;
    double x6;
    double x7;
    double x8;
    int x9;
    x5=rhot1*ut1;
    x6=rhot2*ut2;
    x7=(assert(IS_POSITIVE((-rt1 + x5)*(-rt1 + x5) + (-rt2 + x6)*(-rt2 + x6))), sqrt((-rt1 + x5)*(-rt1 + x5) + (-rt2 + x6)*(-rt2 + x6)));
    x8=mu*x2;
    x9=x3 && x7 <= 0 || x4 && x7 <= x8;
    int x10;
    int x11;
    x10=x3 && x7 > 0 || x4 && x7 > x8;
    x11=x10 && x3;
    int x13;
    double x12;
    x13=x10 && x4;
    /* Assignment result[0, 0]=Piecewise((rn, x3), (x1, x4)) */
    if (x3)
    {
        DEBUG_PRINT("Case (x3) is True.\n");

        /* Assignment result[0, 0]=rn */
        result[0] = rn;
    }
    else if (x4)
    {
        DEBUG_PRINT("Case (x4) is True.\n");

        /* Assignment result[0, 0]=x1 */
        result[0] = x1;
    }

    /* Assignment result[1, 0]=Piecewise((x5, x9), (rt1, x11), (rt1 - x12*(rt1 - x5), x13)) */
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");

        /* Assignment result[1, 0]=x5 */
        result[1] = x5;
    }
    else if (x11)
    {
        DEBUG_PRINT("Case (x11) is True.\n");

        /* Assignment result[1, 0]=rt1 */
        result[1] = rt1;
    }
    else if (x13)
    {
        DEBUG_PRINT("Case (x13) is True.\n");
        x12=(assert(IS_NOT_ZERO(x7)), mu*x2/x7);

        /* Assignment result[1, 0]=rt1 - x12*(rt1 - x5) */
        result[1] = rt1 - x12*(rt1 - x5);
    }

    /* Assignment result[2, 0]=Piecewise((x6, x9), (rt2, x11), (rt2 - x12*(rt2 - x6), x13)) */
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");

        /* Assignment result[2, 0]=x6 */
        result[2] = x6;
    }
    else if (x11)
    {
        DEBUG_PRINT("Case (x11) is True.\n");

        /* Assignment result[2, 0]=rt2 */
        result[2] = rt2;
    }
    else if (x13)
    {
        DEBUG_PRINT("Case (x13) is True.\n");
        x12=(assert(IS_NOT_ZERO(x7)), mu*x2/x7);

        /* Assignment result[2, 0]=rt2 - x12*(rt2 - x6) */
        result[2] = rt2 - x12*(rt2 - x6);
    }
}
void fc3d_AlartCurnierABGenerated(
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
    int x2;
    x1=-rhon*un + rn;
    x2=x1 <= 0;
    int x3;
    x3=x1 > 0;
    double x4;
    double x5;
    double x6;
    double x7;
    double x8;
    double x9;
    double x10;
    int x11;
    int x12;
    int x13;
    int x14;
    int x15;
    x4=rhot1*ut1;
    x5=-rt1 + x4;
    x6=rhot2*ut2;
    x7=-rt2 + x6;
    x8=x5*x5 + x7*x7;
    x9=(assert(IS_POSITIVE(x8)), sqrt(x8));
    x10=mu*x1;
    x11=x9 <= x10;
    x12=x9 <= 0;
    x13=x9 > x10;
    x14=x9 > 0;
    x15=(x11 || x2) && (x11 || x12 || x2) && (x11 || x13 || x2) && (x11 || x2 || x3) && (x12 || x2 || x3) && (x11 || x12 || x13 || x2) && (x11 || x12 || x2 || x3) && (x11 || x13 || x14 || x2) && (x11 || x14 || x2 || x3);
    int x20;
    int x21;
    double x16;
    double x17;
    double x18;
    double x19;
    x20=x13 && x3 || x14 && x2;
    x21=x20 && x3;
    int x22;
    x22=x11 && x3 || x12 && x2;
    int x23;
    x23=x2 && x20;
    /* Assignment result[0, 0]=Piecewise((0, x2), (rhon, x3)) */
    if (x2)
    {
        DEBUG_PRINT("Case (x2) is True.\n");

        /* Assignment result[0, 0]=0 */
        result[0] = 0;
    }
    else if (x3)
    {
        DEBUG_PRINT("Case (x3) is True.\n");

        /* Assignment result[0, 0]=rhon */
        result[0] = rhon;
    }

    /* Assignment result[1, 0]=Piecewise((0, x15), (rhon*x19, x21)) */
    if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");

        /* Assignment result[1, 0]=0 */
        result[1] = 0;
    }
    else if (x21)
    {
        DEBUG_PRINT("Case (x21) is True.\n");
        x16=rt1 - x4;
        x17=1.0/(assert(IS_NOT_ZERO(x9)), x9);
        x18=mu*x17;
        x19=x16*x18;

        /* Assignment result[1, 0]=rhon*x19 */
        result[1] = rhon*x19;
    }

    /* Assignment result[2, 0]=Piecewise((0, x15), (rhon*x31, x21)) */
    double x30;
    double x31;if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");

        /* Assignment result[2, 0]=0 */
        result[2] = 0;
    }
    else if (x21)
    {
        DEBUG_PRINT("Case (x21) is True.\n");
        x17=1.0/(assert(IS_NOT_ZERO(x9)), x9);
        x18=mu*x17;
        x30=rt2 - x6;
        x31=x18*x30;

        /* Assignment result[2, 0]=rhon*x31 */
        result[2] = rhon*x31;
    }

    /* Assignment result[0, 1]=0 */
    result[3] = 0;

    /* Assignment result[1, 1]=Piecewise((rhot1, x22), (0, x23), (rhot1*x24 + rhot1*x27, x21)) */
    double x24;
    double x25;
    double x26;
    double x27;if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");

        /* Assignment result[1, 1]=rhot1 */
        result[4] = rhot1;
    }
    else if (x23)
    {
        DEBUG_PRINT("Case (x23) is True.\n");

        /* Assignment result[1, 1]=0 */
        result[4] = 0;
    }
    else if (x21)
    {
        DEBUG_PRINT("Case (x21) is True.\n");
        x16=rt1 - x4;
        x17=1.0/(assert(IS_NOT_ZERO(x9)), x9);
        x24=x10*x17;
        x25=pow(x8, -3.0/2.0);
        x26=mu*x1*x16*x25;
        x27=x26*x5;

        /* Assignment result[1, 1]=rhot1*x24 + rhot1*x27 */
        result[4] = rhot1*x24 + rhot1*x27;
    }

    /* Assignment result[2, 1]=Piecewise((0, x15), (rhot1*x33, x21)) */
    double x32;
    double x33;if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");

        /* Assignment result[2, 1]=0 */
        result[5] = 0;
    }
    else if (x21)
    {
        DEBUG_PRINT("Case (x21) is True.\n");
        x25=pow(x8, -3.0/2.0);
        x30=rt2 - x6;
        x32=mu*x1*x25*x30;
        x33=x32*x5;

        /* Assignment result[2, 1]=rhot1*x33 */
        result[5] = rhot1*x33;
    }

    /* Assignment result[0, 2]=0 */
    result[6] = 0;

    /* Assignment result[1, 2]=Piecewise((0, x15), (rhot2*x28, x21)) */
    double x28;if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");

        /* Assignment result[1, 2]=0 */
        result[7] = 0;
    }
    else if (x21)
    {
        DEBUG_PRINT("Case (x21) is True.\n");
        x16=rt1 - x4;
        x25=pow(x8, -3.0/2.0);
        x26=mu*x1*x16*x25;
        x28=x26*x7;

        /* Assignment result[1, 2]=rhot2*x28 */
        result[7] = rhot2*x28;
    }

    /* Assignment result[2, 2]=Piecewise((rhot2, x22), (0, x23), (rhot2*x24 + rhot2*x34, x21)) */
    double x34;if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");

        /* Assignment result[2, 2]=rhot2 */
        result[8] = rhot2;
    }
    else if (x23)
    {
        DEBUG_PRINT("Case (x23) is True.\n");

        /* Assignment result[2, 2]=0 */
        result[8] = 0;
    }
    else if (x21)
    {
        DEBUG_PRINT("Case (x21) is True.\n");
        x17=1.0/(assert(IS_NOT_ZERO(x9)), x9);
        x24=x10*x17;
        x25=pow(x8, -3.0/2.0);
        x30=rt2 - x6;
        x32=mu*x1*x25*x30;
        x34=x32*x7;

        /* Assignment result[2, 2]=rhot2*x24 + rhot2*x34 */
        result[8] = rhot2*x24 + rhot2*x34;
    }

    /* Assignment result[0, 3]=Piecewise((1, x2), (0, x3)) */
    if (x2)
    {
        DEBUG_PRINT("Case (x2) is True.\n");

        /* Assignment result[0, 3]=1 */
        result[9] = 1;
    }
    else if (x3)
    {
        DEBUG_PRINT("Case (x3) is True.\n");

        /* Assignment result[0, 3]=0 */
        result[9] = 0;
    }

    /* Assignment result[1, 3]=Piecewise((0, x15), (-x19, x21)) */
    if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");

        /* Assignment result[1, 3]=0 */
        result[10] = 0;
    }
    else if (x21)
    {
        DEBUG_PRINT("Case (x21) is True.\n");
        x16=rt1 - x4;
        x17=1.0/(assert(IS_NOT_ZERO(x9)), x9);
        x18=mu*x17;
        x19=x16*x18;

        /* Assignment result[1, 3]=-x19 */
        result[10] = -x19;
    }

    /* Assignment result[2, 3]=Piecewise((0, x15), (-x31, x21)) */
    if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");

        /* Assignment result[2, 3]=0 */
        result[11] = 0;
    }
    else if (x21)
    {
        DEBUG_PRINT("Case (x21) is True.\n");
        x17=1.0/(assert(IS_NOT_ZERO(x9)), x9);
        x18=mu*x17;
        x30=rt2 - x6;
        x31=x18*x30;

        /* Assignment result[2, 3]=-x31 */
        result[11] = -x31;
    }

    /* Assignment result[0, 4]=0 */
    result[12] = 0;

    /* Assignment result[1, 4]=Piecewise((0, x22), (1, x23), (-x27 + x29, x21)) */
    double x29;if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");

        /* Assignment result[1, 4]=0 */
        result[13] = 0;
    }
    else if (x23)
    {
        DEBUG_PRINT("Case (x23) is True.\n");

        /* Assignment result[1, 4]=1 */
        result[13] = 1;
    }
    else if (x21)
    {
        DEBUG_PRINT("Case (x21) is True.\n");
        x16=rt1 - x4;
        x17=1.0/(assert(IS_NOT_ZERO(x9)), x9);
        x24=x10*x17;
        x25=pow(x8, -3.0/2.0);
        x26=mu*x1*x16*x25;
        x27=x26*x5;
        x29=-x24 + 1;

        /* Assignment result[1, 4]=-x27 + x29 */
        result[13] = -x27 + x29;
    }

    /* Assignment result[2, 4]=Piecewise((0, x15), (-x33, x21)) */
    if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");

        /* Assignment result[2, 4]=0 */
        result[14] = 0;
    }
    else if (x21)
    {
        DEBUG_PRINT("Case (x21) is True.\n");
        x25=pow(x8, -3.0/2.0);
        x30=rt2 - x6;
        x32=mu*x1*x25*x30;
        x33=x32*x5;

        /* Assignment result[2, 4]=-x33 */
        result[14] = -x33;
    }

    /* Assignment result[0, 5]=0 */
    result[15] = 0;

    /* Assignment result[1, 5]=Piecewise((0, x15), (-x28, x21)) */
    if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");

        /* Assignment result[1, 5]=0 */
        result[16] = 0;
    }
    else if (x21)
    {
        DEBUG_PRINT("Case (x21) is True.\n");
        x16=rt1 - x4;
        x25=pow(x8, -3.0/2.0);
        x26=mu*x1*x16*x25;
        x28=x26*x7;

        /* Assignment result[1, 5]=-x28 */
        result[16] = -x28;
    }

    /* Assignment result[2, 5]=Piecewise((0, x22), (1, x23), (x29 - x34, x21)) */
    if (x22)
    {
        DEBUG_PRINT("Case (x22) is True.\n");

        /* Assignment result[2, 5]=0 */
        result[17] = 0;
    }
    else if (x23)
    {
        DEBUG_PRINT("Case (x23) is True.\n");

        /* Assignment result[2, 5]=1 */
        result[17] = 1;
    }
    else if (x21)
    {
        DEBUG_PRINT("Case (x21) is True.\n");
        x17=1.0/(assert(IS_NOT_ZERO(x9)), x9);
        x24=x10*x17;
        x25=pow(x8, -3.0/2.0);
        x29=-x24 + 1;
        x30=rt2 - x6;
        x32=mu*x1*x25*x30;
        x34=x32*x7;

        /* Assignment result[2, 5]=x29 - x34 */
        result[17] = x29 - x34;
    }
}

void fc3d_AlartCurnierFunctionGenerated(
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

    fc3d_AlartCurnierFABGenerated(
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
      fc3d_AlartCurnierFGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3(result, f);
    }

    if (A && B)
    {
      fc3d_AlartCurnierABGenerated(
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

void fc3d_AlartCurnierJeanMoreauFABGenerated(
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
    int x3;
    x1=rhon*un;
    x2=rn - x1;
    x3=x2 <= 0;
    int x4;
    x4=x2 > 0;
    double x5;
    int x6;
    double x7;
    double x8;
    double x9;
    double x10;
    double x11;
    int x12;
    int x13;
    double x14;
    int x15;
    int x16;
    x5=rhot1*ut1;
    x6=rn <= 0;
    x7=-rt1 + x5;
    x8=rhot2*ut2;
    x9=-rt2 + x8;
    x10=x7*x7 + x9*x9;
    x11=(assert(IS_POSITIVE(x10)), sqrt(x10));
    x12=x11 <= 0;
    x13=rn > 0;
    x14=mu*rn;
    x15=x11 <= x14;
    x16=x12 && x6 || x13 && x15;
    int x17;
    int x18;
    int x19;
    int x20;
    x17=x11 > 0;
    x18=x11 > x14;
    x19=x13 && x18 || x17 && x6;
    x20=x19 && x6;
    int x24;
    double x21;
    double x22;
    double x23;
    x24=x13 && x19;
    int x28;
    x28=(x15 || x6) && (x12 || x13 || x6) && (x12 || x15 || x6) && (x13 || x15 || x6) && (x15 || x18 || x6) && (x12 || x13 || x15 || x6) && (x12 || x15 || x18 || x6) && (x13 || x15 || x17 || x6) && (x15 || x17 || x18 || x6);
    /* Assignment result[0, 0]=Piecewise((rn, x3), (x1, x4)) */
    if (x3)
    {
        DEBUG_PRINT("Case (x3) is True.\n");

        /* Assignment result[0, 0]=rn */
        result[0] = rn;
    }
    else if (x4)
    {
        DEBUG_PRINT("Case (x4) is True.\n");

        /* Assignment result[0, 0]=x1 */
        result[0] = x1;
    }

    /* Assignment result[1, 0]=Piecewise((x5, x16), (rt1, x20), (rt1 - x21*x23, x24)) */
    if (x16)
    {
        DEBUG_PRINT("Case (x16) is True.\n");

        /* Assignment result[1, 0]=x5 */
        result[1] = x5;
    }
    else if (x20)
    {
        DEBUG_PRINT("Case (x20) is True.\n");

        /* Assignment result[1, 0]=rt1 */
        result[1] = rt1;
    }
    else if (x24)
    {
        DEBUG_PRINT("Case (x24) is True.\n");
        x21=rt1 - x5;
        x22=1.0/(assert(IS_NOT_ZERO(x11)), x11);
        x23=x14*x22;

        /* Assignment result[1, 0]=rt1 - x21*x23 */
        result[1] = rt1 - x21*x23;
    }

    /* Assignment result[2, 0]=Piecewise((x8, x16), (rt2, x20), (rt2 - x23*x32, x24)) */
    double x32;if (x16)
    {
        DEBUG_PRINT("Case (x16) is True.\n");

        /* Assignment result[2, 0]=x8 */
        result[2] = x8;
    }
    else if (x20)
    {
        DEBUG_PRINT("Case (x20) is True.\n");

        /* Assignment result[2, 0]=rt2 */
        result[2] = rt2;
    }
    else if (x24)
    {
        DEBUG_PRINT("Case (x24) is True.\n");
        x22=1.0/(assert(IS_NOT_ZERO(x11)), x11);
        x23=x14*x22;
        x32=rt2 - x8;

        /* Assignment result[2, 0]=rt2 - x23*x32 */
        result[2] = rt2 - x23*x32;
    }

    /* Assignment result[0, 1]=Piecewise((0, x3), (rhon, x4)) */
    if (x3)
    {
        DEBUG_PRINT("Case (x3) is True.\n");

        /* Assignment result[0, 1]=0 */
        result[3] = 0;
    }
    else if (x4)
    {
        DEBUG_PRINT("Case (x4) is True.\n");

        /* Assignment result[0, 1]=rhon */
        result[3] = rhon;
    }

    /* Assignment result[1, 1]=0 */
    result[4] = 0;

    /* Assignment result[2, 1]=0 */
    result[5] = 0;

    /* Assignment result[0, 2]=0 */
    result[6] = 0;

    /* Assignment result[1, 2]=Piecewise((rhot1, x16), (0, x20), (rhot1*x23 + rhot1*x27, x24)) */
    double x25;
    double x26;
    double x27;if (x16)
    {
        DEBUG_PRINT("Case (x16) is True.\n");

        /* Assignment result[1, 2]=rhot1 */
        result[7] = rhot1;
    }
    else if (x20)
    {
        DEBUG_PRINT("Case (x20) is True.\n");

        /* Assignment result[1, 2]=0 */
        result[7] = 0;
    }
    else if (x24)
    {
        DEBUG_PRINT("Case (x24) is True.\n");
        x21=rt1 - x5;
        x22=1.0/(assert(IS_NOT_ZERO(x11)), x11);
        x23=x14*x22;
        x25=pow(x10, -3.0/2.0);
        x26=mu*rn*x21*x25;
        x27=x26*x7;

        /* Assignment result[1, 2]=rhot1*x23 + rhot1*x27 */
        result[7] = rhot1*x23 + rhot1*x27;
    }

    /* Assignment result[2, 2]=Piecewise((0, x28), (rhot1*x34, x24)) */
    double x33;
    double x34;if (x28)
    {
        DEBUG_PRINT("Case (x28) is True.\n");

        /* Assignment result[2, 2]=0 */
        result[8] = 0;
    }
    else if (x24)
    {
        DEBUG_PRINT("Case (x24) is True.\n");
        x25=pow(x10, -3.0/2.0);
        x32=rt2 - x8;
        x33=mu*rn*x25*x32;
        x34=x33*x7;

        /* Assignment result[2, 2]=rhot1*x34 */
        result[8] = rhot1*x34;
    }

    /* Assignment result[0, 3]=0 */
    result[9] = 0;

    /* Assignment result[1, 3]=Piecewise((0, x28), (rhot2*x29, x24)) */
    double x29;if (x28)
    {
        DEBUG_PRINT("Case (x28) is True.\n");

        /* Assignment result[1, 3]=0 */
        result[10] = 0;
    }
    else if (x24)
    {
        DEBUG_PRINT("Case (x24) is True.\n");
        x21=rt1 - x5;
        x25=pow(x10, -3.0/2.0);
        x26=mu*rn*x21*x25;
        x29=x26*x9;

        /* Assignment result[1, 3]=rhot2*x29 */
        result[10] = rhot2*x29;
    }

    /* Assignment result[2, 3]=Piecewise((rhot2, x16), (0, x20), (rhot2*x23 + rhot2*x35, x24)) */
    double x35;if (x16)
    {
        DEBUG_PRINT("Case (x16) is True.\n");

        /* Assignment result[2, 3]=rhot2 */
        result[11] = rhot2;
    }
    else if (x20)
    {
        DEBUG_PRINT("Case (x20) is True.\n");

        /* Assignment result[2, 3]=0 */
        result[11] = 0;
    }
    else if (x24)
    {
        DEBUG_PRINT("Case (x24) is True.\n");
        x22=1.0/(assert(IS_NOT_ZERO(x11)), x11);
        x23=x14*x22;
        x25=pow(x10, -3.0/2.0);
        x32=rt2 - x8;
        x33=mu*rn*x25*x32;
        x35=x33*x9;

        /* Assignment result[2, 3]=rhot2*x23 + rhot2*x35 */
        result[11] = rhot2*x23 + rhot2*x35;
    }

    /* Assignment result[0, 4]=Piecewise((1, x3), (0, x4)) */
    if (x3)
    {
        DEBUG_PRINT("Case (x3) is True.\n");

        /* Assignment result[0, 4]=1 */
        result[12] = 1;
    }
    else if (x4)
    {
        DEBUG_PRINT("Case (x4) is True.\n");

        /* Assignment result[0, 4]=0 */
        result[12] = 0;
    }

    /* Assignment result[1, 4]=Piecewise((0, x28), (-x21*x30, x24)) */
    double x30;if (x28)
    {
        DEBUG_PRINT("Case (x28) is True.\n");

        /* Assignment result[1, 4]=0 */
        result[13] = 0;
    }
    else if (x24)
    {
        DEBUG_PRINT("Case (x24) is True.\n");
        x21=rt1 - x5;
        x22=1.0/(assert(IS_NOT_ZERO(x11)), x11);
        x30=mu*x22;

        /* Assignment result[1, 4]=-x21*x30 */
        result[13] = -x21*x30;
    }

    /* Assignment result[2, 4]=Piecewise((0, x28), (-x30*x32, x24)) */
    if (x28)
    {
        DEBUG_PRINT("Case (x28) is True.\n");

        /* Assignment result[2, 4]=0 */
        result[14] = 0;
    }
    else if (x24)
    {
        DEBUG_PRINT("Case (x24) is True.\n");
        x22=1.0/(assert(IS_NOT_ZERO(x11)), x11);
        x30=mu*x22;
        x32=rt2 - x8;

        /* Assignment result[2, 4]=-x30*x32 */
        result[14] = -x30*x32;
    }

    /* Assignment result[0, 5]=0 */
    result[15] = 0;

    /* Assignment result[1, 5]=Piecewise((0, x16), (1, x20), (-x27 + x31, x24)) */
    double x31;if (x16)
    {
        DEBUG_PRINT("Case (x16) is True.\n");

        /* Assignment result[1, 5]=0 */
        result[16] = 0;
    }
    else if (x20)
    {
        DEBUG_PRINT("Case (x20) is True.\n");

        /* Assignment result[1, 5]=1 */
        result[16] = 1;
    }
    else if (x24)
    {
        DEBUG_PRINT("Case (x24) is True.\n");
        x21=rt1 - x5;
        x22=1.0/(assert(IS_NOT_ZERO(x11)), x11);
        x23=x14*x22;
        x25=pow(x10, -3.0/2.0);
        x26=mu*rn*x21*x25;
        x27=x26*x7;
        x31=-x23 + 1;

        /* Assignment result[1, 5]=-x27 + x31 */
        result[16] = -x27 + x31;
    }

    /* Assignment result[2, 5]=Piecewise((0, x28), (-x34, x24)) */
    if (x28)
    {
        DEBUG_PRINT("Case (x28) is True.\n");

        /* Assignment result[2, 5]=0 */
        result[17] = 0;
    }
    else if (x24)
    {
        DEBUG_PRINT("Case (x24) is True.\n");
        x25=pow(x10, -3.0/2.0);
        x32=rt2 - x8;
        x33=mu*rn*x25*x32;
        x34=x33*x7;

        /* Assignment result[2, 5]=-x34 */
        result[17] = -x34;
    }

    /* Assignment result[0, 6]=0 */
    result[18] = 0;

    /* Assignment result[1, 6]=Piecewise((0, x28), (-x29, x24)) */
    if (x28)
    {
        DEBUG_PRINT("Case (x28) is True.\n");

        /* Assignment result[1, 6]=0 */
        result[19] = 0;
    }
    else if (x24)
    {
        DEBUG_PRINT("Case (x24) is True.\n");
        x21=rt1 - x5;
        x25=pow(x10, -3.0/2.0);
        x26=mu*rn*x21*x25;
        x29=x26*x9;

        /* Assignment result[1, 6]=-x29 */
        result[19] = -x29;
    }

    /* Assignment result[2, 6]=Piecewise((0, x16), (1, x20), (x31 - x35, x24)) */
    if (x16)
    {
        DEBUG_PRINT("Case (x16) is True.\n");

        /* Assignment result[2, 6]=0 */
        result[20] = 0;
    }
    else if (x20)
    {
        DEBUG_PRINT("Case (x20) is True.\n");

        /* Assignment result[2, 6]=1 */
        result[20] = 1;
    }
    else if (x24)
    {
        DEBUG_PRINT("Case (x24) is True.\n");
        x22=1.0/(assert(IS_NOT_ZERO(x11)), x11);
        x23=x14*x22;
        x25=pow(x10, -3.0/2.0);
        x31=-x23 + 1;
        x32=rt2 - x8;
        x33=mu*rn*x25*x32;
        x35=x33*x9;

        /* Assignment result[2, 6]=x31 - x35 */
        result[20] = x31 - x35;
    }
}
void fc3d_AlartCurnierJeanMoreauFGenerated(
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
    x1=rhon*un;
    x2=rn - x1;
    double x3;
    int x4;
    double x5;
    double x6;
    int x7;
    double x8;
    int x9;
    x3=rhot1*ut1;
    x4=rn <= 0;
    x5=rhot2*ut2;
    x6=(assert(IS_POSITIVE((-rt1 + x3)*(-rt1 + x3) + (-rt2 + x5)*(-rt2 + x5))), sqrt((-rt1 + x3)*(-rt1 + x3) + (-rt2 + x5)*(-rt2 + x5)));
    x7=rn > 0;
    x8=mu*rn;
    x9=x4 && x6 <= 0 || x7 && x6 <= x8;
    int x10;
    int x11;
    x10=x4 && x6 > 0 || x7 && x6 > x8;
    x11=x10 && x4;
    int x13;
    double x12;
    x13=x10 && x7;
    /* Assignment result[0, 0]=Piecewise((rn, x2 <= 0), (x1, x2 > 0)) */
    if (x2 <= 0)
    {
        DEBUG_PRINT("Case (x2 <= 0) is True.\n");

        /* Assignment result[0, 0]=rn */
        result[0] = rn;
    }
    else if (x2 > 0)
    {
        DEBUG_PRINT("Case (x2 > 0) is True.\n");

        /* Assignment result[0, 0]=x1 */
        result[0] = x1;
    }

    /* Assignment result[1, 0]=Piecewise((x3, x9), (rt1, x11), (rt1 - x12*(rt1 - x3), x13)) */
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");

        /* Assignment result[1, 0]=x3 */
        result[1] = x3;
    }
    else if (x11)
    {
        DEBUG_PRINT("Case (x11) is True.\n");

        /* Assignment result[1, 0]=rt1 */
        result[1] = rt1;
    }
    else if (x13)
    {
        DEBUG_PRINT("Case (x13) is True.\n");
        x12=(assert(IS_NOT_ZERO(x6)), mu*rn/x6);

        /* Assignment result[1, 0]=rt1 - x12*(rt1 - x3) */
        result[1] = rt1 - x12*(rt1 - x3);
    }

    /* Assignment result[2, 0]=Piecewise((x5, x9), (rt2, x11), (rt2 - x12*(rt2 - x5), x13)) */
    if (x9)
    {
        DEBUG_PRINT("Case (x9) is True.\n");

        /* Assignment result[2, 0]=x5 */
        result[2] = x5;
    }
    else if (x11)
    {
        DEBUG_PRINT("Case (x11) is True.\n");

        /* Assignment result[2, 0]=rt2 */
        result[2] = rt2;
    }
    else if (x13)
    {
        DEBUG_PRINT("Case (x13) is True.\n");
        x12=(assert(IS_NOT_ZERO(x6)), mu*rn/x6);

        /* Assignment result[2, 0]=rt2 - x12*(rt2 - x5) */
        result[2] = rt2 - x12*(rt2 - x5);
    }
}
void fc3d_AlartCurnierJeanMoreauABGenerated(
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
    int x2;
    x1=-rhon*un + rn;
    x2=x1 <= 0;
    int x3;
    x3=x1 > 0;
    int x4;
    double x5;
    double x6;
    double x7;
    double x8;
    double x9;
    double x10;
    int x11;
    int x12;
    double x13;
    int x14;
    int x15;
    x4=rn <= 0;
    x5=rhot1*ut1;
    x6=-rt1 + x5;
    x7=rhot2*ut2;
    x8=-rt2 + x7;
    x9=x6*x6 + x8*x8;
    x10=(assert(IS_POSITIVE(x9)), sqrt(x9));
    x11=x10 <= 0;
    x12=rn > 0;
    x13=mu*rn;
    x14=x10 <= x13;
    x15=x11 && x4 || x12 && x14;
    int x16;
    int x17;
    int x18;
    int x19;
    x16=x10 > 0;
    x17=x10 > x13;
    x18=x12 && x17 || x16 && x4;
    x19=x18 && x4;
    int x26;
    double x20;
    double x21;
    double x22;
    double x23;
    double x24;
    double x25;
    x26=x12 && x18;
    int x27;
    x27=(x14 || x4) && (x11 || x12 || x4) && (x11 || x14 || x4) && (x12 || x14 || x4) && (x14 || x17 || x4) && (x11 || x12 || x14 || x4) && (x11 || x14 || x17 || x4) && (x12 || x14 || x16 || x4) && (x14 || x16 || x17 || x4);
    /* Assignment result[0, 0]=Piecewise((0, x2), (rhon, x3)) */
    if (x2)
    {
        DEBUG_PRINT("Case (x2) is True.\n");

        /* Assignment result[0, 0]=0 */
        result[0] = 0;
    }
    else if (x3)
    {
        DEBUG_PRINT("Case (x3) is True.\n");

        /* Assignment result[0, 0]=rhon */
        result[0] = rhon;
    }

    /* Assignment result[1, 0]=0 */
    result[1] = 0;

    /* Assignment result[2, 0]=0 */
    result[2] = 0;

    /* Assignment result[0, 1]=0 */
    result[3] = 0;

    /* Assignment result[1, 1]=Piecewise((rhot1, x15), (0, x19), (rhot1*x21 + rhot1*x25, x26)) */
    if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");

        /* Assignment result[1, 1]=rhot1 */
        result[4] = rhot1;
    }
    else if (x19)
    {
        DEBUG_PRINT("Case (x19) is True.\n");

        /* Assignment result[1, 1]=0 */
        result[4] = 0;
    }
    else if (x26)
    {
        DEBUG_PRINT("Case (x26) is True.\n");
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x21=x13*x20;
        x22=rt1 - x5;
        x23=pow(x9, -3.0/2.0);
        x24=mu*rn*x22*x23;
        x25=x24*x6;

        /* Assignment result[1, 1]=rhot1*x21 + rhot1*x25 */
        result[4] = rhot1*x21 + rhot1*x25;
    }

    /* Assignment result[2, 1]=Piecewise((0, x27), (rhot1*x33, x26)) */
    double x31;
    double x32;
    double x33;if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");

        /* Assignment result[2, 1]=0 */
        result[5] = 0;
    }
    else if (x26)
    {
        DEBUG_PRINT("Case (x26) is True.\n");
        x23=pow(x9, -3.0/2.0);
        x31=rt2 - x7;
        x32=mu*rn*x23*x31;
        x33=x32*x6;

        /* Assignment result[2, 1]=rhot1*x33 */
        result[5] = rhot1*x33;
    }

    /* Assignment result[0, 2]=0 */
    result[6] = 0;

    /* Assignment result[1, 2]=Piecewise((0, x27), (rhot2*x28, x26)) */
    double x28;if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");

        /* Assignment result[1, 2]=0 */
        result[7] = 0;
    }
    else if (x26)
    {
        DEBUG_PRINT("Case (x26) is True.\n");
        x22=rt1 - x5;
        x23=pow(x9, -3.0/2.0);
        x24=mu*rn*x22*x23;
        x28=x24*x8;

        /* Assignment result[1, 2]=rhot2*x28 */
        result[7] = rhot2*x28;
    }

    /* Assignment result[2, 2]=Piecewise((rhot2, x15), (0, x19), (rhot2*x21 + rhot2*x34, x26)) */
    double x34;if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");

        /* Assignment result[2, 2]=rhot2 */
        result[8] = rhot2;
    }
    else if (x19)
    {
        DEBUG_PRINT("Case (x19) is True.\n");

        /* Assignment result[2, 2]=0 */
        result[8] = 0;
    }
    else if (x26)
    {
        DEBUG_PRINT("Case (x26) is True.\n");
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x21=x13*x20;
        x23=pow(x9, -3.0/2.0);
        x31=rt2 - x7;
        x32=mu*rn*x23*x31;
        x34=x32*x8;

        /* Assignment result[2, 2]=rhot2*x21 + rhot2*x34 */
        result[8] = rhot2*x21 + rhot2*x34;
    }

    /* Assignment result[0, 3]=Piecewise((1, x2), (0, x3)) */
    if (x2)
    {
        DEBUG_PRINT("Case (x2) is True.\n");

        /* Assignment result[0, 3]=1 */
        result[9] = 1;
    }
    else if (x3)
    {
        DEBUG_PRINT("Case (x3) is True.\n");

        /* Assignment result[0, 3]=0 */
        result[9] = 0;
    }

    /* Assignment result[1, 3]=Piecewise((0, x27), (-x22*x29, x26)) */
    double x29;if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");

        /* Assignment result[1, 3]=0 */
        result[10] = 0;
    }
    else if (x26)
    {
        DEBUG_PRINT("Case (x26) is True.\n");
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x22=rt1 - x5;
        x29=mu*x20;

        /* Assignment result[1, 3]=-x22*x29 */
        result[10] = -x22*x29;
    }

    /* Assignment result[2, 3]=Piecewise((0, x27), (-x29*x31, x26)) */
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");

        /* Assignment result[2, 3]=0 */
        result[11] = 0;
    }
    else if (x26)
    {
        DEBUG_PRINT("Case (x26) is True.\n");
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x29=mu*x20;
        x31=rt2 - x7;

        /* Assignment result[2, 3]=-x29*x31 */
        result[11] = -x29*x31;
    }

    /* Assignment result[0, 4]=0 */
    result[12] = 0;

    /* Assignment result[1, 4]=Piecewise((0, x15), (1, x19), (-x25 + x30, x26)) */
    double x30;if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");

        /* Assignment result[1, 4]=0 */
        result[13] = 0;
    }
    else if (x19)
    {
        DEBUG_PRINT("Case (x19) is True.\n");

        /* Assignment result[1, 4]=1 */
        result[13] = 1;
    }
    else if (x26)
    {
        DEBUG_PRINT("Case (x26) is True.\n");
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x21=x13*x20;
        x22=rt1 - x5;
        x23=pow(x9, -3.0/2.0);
        x24=mu*rn*x22*x23;
        x25=x24*x6;
        x30=-x21 + 1;

        /* Assignment result[1, 4]=-x25 + x30 */
        result[13] = -x25 + x30;
    }

    /* Assignment result[2, 4]=Piecewise((0, x27), (-x33, x26)) */
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");

        /* Assignment result[2, 4]=0 */
        result[14] = 0;
    }
    else if (x26)
    {
        DEBUG_PRINT("Case (x26) is True.\n");
        x23=pow(x9, -3.0/2.0);
        x31=rt2 - x7;
        x32=mu*rn*x23*x31;
        x33=x32*x6;

        /* Assignment result[2, 4]=-x33 */
        result[14] = -x33;
    }

    /* Assignment result[0, 5]=0 */
    result[15] = 0;

    /* Assignment result[1, 5]=Piecewise((0, x27), (-x28, x26)) */
    if (x27)
    {
        DEBUG_PRINT("Case (x27) is True.\n");

        /* Assignment result[1, 5]=0 */
        result[16] = 0;
    }
    else if (x26)
    {
        DEBUG_PRINT("Case (x26) is True.\n");
        x22=rt1 - x5;
        x23=pow(x9, -3.0/2.0);
        x24=mu*rn*x22*x23;
        x28=x24*x8;

        /* Assignment result[1, 5]=-x28 */
        result[16] = -x28;
    }

    /* Assignment result[2, 5]=Piecewise((0, x15), (1, x19), (x30 - x34, x26)) */
    if (x15)
    {
        DEBUG_PRINT("Case (x15) is True.\n");

        /* Assignment result[2, 5]=0 */
        result[17] = 0;
    }
    else if (x19)
    {
        DEBUG_PRINT("Case (x19) is True.\n");

        /* Assignment result[2, 5]=1 */
        result[17] = 1;
    }
    else if (x26)
    {
        DEBUG_PRINT("Case (x26) is True.\n");
        x20=1.0/(assert(IS_NOT_ZERO(x10)), x10);
        x21=x13*x20;
        x23=pow(x9, -3.0/2.0);
        x30=-x21 + 1;
        x31=rt2 - x7;
        x32=mu*rn*x23*x31;
        x34=x32*x8;

        /* Assignment result[2, 5]=x30 - x34 */
        result[17] = x30 - x34;
    }
}

void fc3d_AlartCurnierJeanMoreauFunctionGenerated(
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

    fc3d_AlartCurnierJeanMoreauFABGenerated(
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
      fc3d_AlartCurnierJeanMoreauFGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3(result, f);
    }

    if (A && B)
    {
      fc3d_AlartCurnierJeanMoreauABGenerated(
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
