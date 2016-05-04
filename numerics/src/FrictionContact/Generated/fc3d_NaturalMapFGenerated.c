#include "fc3d_NaturalMapFGenerated.h"
#include "funcodegen.h"
/*@
requires (0.0 <= rn <= 1.0e+6);
requires (-1.0e+6 <= rt1 <= 1.0e+6);
requires (-1.0e+6 <= rt2 <= 1.0e+6);
requires (-1.0e+6 <= un <= 1.0e+6);
requires (-1.0e+6 <= ut1 <= 1.0e+6);
requires (-1.0e+6 <= ut2 <= 1.0e+6);
requires (0.0 <= mu <= 1.0);
requires (-1.0e+6 <= rhon <= 1.0e+6);
requires (-1.0e+6 <= rhot1 <= 1.0e+6);
requires (-1.0e+6 <= rhot2 <= 1.0e+6);
assigns result[0..2];
ensures \is_finite((double) result[0]);
ensures \is_finite((double) result[1]);
ensures \is_finite((double) result[2]);*/
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
    /*@ assert \is_finite((double) ut1); */
    /*@ assert \is_finite((double) epsilon); */
    /*@ assert \is_finite((double) un); */
    /*@ assert \is_finite((double) rt1); */
    /*@ assert \is_finite((double) mu); */
    /*@ assert \is_finite((double) ut2); */
    /*@ assert \is_finite((double) rn); */
    /*@ assert \is_finite((double) rt2); */
    double x9 = 0.;
    double x10 = 0.;
    double x11 = 0.;
    double x12 = 0.;
    double x13 = 0.;
    double x14 = 0.;
    double x15 = 0.;
    double x24 = 0.;
    int x25 = 0;

    double x2 = 0.;
    double x6 = 0.;
    double x7 = 0.;
    double x8 = 0.;
    double x16 = 0.;
    double x17 = 0.;
    double x19 = 0.;
    double x21 = 0.;
    double x22 = 0.;
    double x23 = 0.;

    x9 = mu*ut1;
    x10 = rt1 - x9;
    x11 = x10*x10;
    /*@ assert (x11) >= 0; */
    x12 = mu*ut2;
    x13 = rt2 - x12;
    x14 = x13*x13;
    /*@ assert (x14) >= 0; */
    x15 = x11 + x14;
    x24 = epsilon*(mu + 1);
    /*@ assert (x24) >= 0; */
    /*@ assert (x24) != 0; */
    x25 = x15 > x24;
    /*@ assert x25 <==> (x15 > x24); */

    int x26 = 0;

    double x18 = 0.;
    double x20 = 0.;

    x26 = x15 <= x24;
    /*@ assert x26 <==> (x15 <= x24); */

    if (x25)
    {
        x2 = ut1*ut1 + ut2*ut2;
        /*@ assert (x2) >= 0; */
        x6 = mu*rn;
        /*@ assert (x6) >= 0; */
        /*@ assert (x2 >= 0); */
        x7 = sqrt(x2);
        x8 = -mu*x7 - un + x6;
        /*@ assert (x15 >= 0); */
        x16 = sqrt(x15);
        x17 = Max(0, x16 + x8);
        /*@ assert (x17) >= 0; */
        x19 = Max(0, -x16 + x8);
        /*@ assert (x19) >= 0; */
        /*@ assert (x16 < -epsilon*epsilon*epsilon || x16 > epsilon*epsilon*epsilon); */
        x21 = 1.0/x16;
        x22 = 0.5*x17*x21;
        x23 = 0.5*x19*x21;

    }
    if (x26)
    {
        x2 = ut1*ut1 + ut2*ut2;
        /*@ assert (x2) >= 0; */
        x6 = mu*rn;
        /*@ assert (x6) >= 0; */
        /*@ assert (x2 >= 0); */
        x7 = sqrt(x2);
        x8 = -mu*x7 - un + x6;
        /*@ assert (x15 >= 0); */
        x16 = sqrt(x15);
        x17 = Max(0, x16 + x8);
        /*@ assert (x17) >= 0; */
        x18 = -0.5*x17;
        x19 = Max(0, -x16 + x8);
        /*@ assert (x19) >= 0; */
        x20 = 0.5*x19;

    }x2 = ut1*ut1 + ut2*ut2;
    /*@ assert (x2) >= 0; */
    x6 = mu*rn;
    /*@ assert (x6) >= 0; */
    /*@ assert (x2 >= 0); */
    x7 = sqrt(x2);
    x8 = -mu*x7 - un + x6;
    /*@ assert (x15 >= 0); */
    x16 = sqrt(x15);
    x17 = Max(0, x16 + x8);
    /*@ assert (x17) >= 0; */
    x18 = -0.5*x17;
    x19 = Max(0, -x16 + x8);
    /*@ assert (x19) >= 0; */
    x20 = 0.5*x19;
    /*@ assigns result[0]; */
    result[0] = x18 - x20 + x6;


    /*@ assert x25 || x26; */
    if (x25)
    {
        /*@ assigns result[1]; */
        result[1] = rt1 - x10*x22 + x10*x23;

    }
    else if (x26)
    {
        /*@ assigns result[1]; */
        result[1] = rt1;

    }


    /*@ assert x25 || x26; */
    if (x25)
    {
        /*@ assigns result[2]; */
        result[2] = rt2 - x13*x22 + x13*x23;

    }
    else if (x26)
    {
        /*@ assigns result[2]; */
        result[2] = rt2 + x18 + x20;

    }
}
#ifdef __FRAMAC__
int main()
{
    double rn =  Frama_C_double_interval(0.0, 1.0e+6);
    double rt1 =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double rt2 =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double un =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double ut1 =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double ut2 =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double mu =  Frama_C_double_interval(0.0, 1.0);
    double rhon =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double rhot1 =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double rhot2 =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double result[3];
    fc3d_NaturalMapFGenerated(rn, rt1, rt2, un, ut1, ut2, mu, rhon, rhot1, rhot2, result);
    return(0);
}
#endif