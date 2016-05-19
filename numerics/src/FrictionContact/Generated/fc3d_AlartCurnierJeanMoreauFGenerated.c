#include "fc3d_AlartCurnierJeanMoreauFGenerated.h"
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
    /*@ assert \is_finite((double) un); */
    /*@ assert \is_finite((double) rn); */
    /*@ assert \is_finite((double) rhon); */
    /*@ assert \is_finite((double) rhot2); */
    /*@ assert \is_finite((double) rhot1); */
    /*@ assert \is_finite((double) ut1); */
    /*@ assert \is_finite((double) rt1); */
    /*@ assert \is_finite((double) mu); */
    /*@ assert \is_finite((double) ut2); */
    /*@ assert \is_finite((double) rt2); */
    double x3 = 0.;
    double x4 = 0.;
    double x5 = 0.;
    double x6 = 0.;
    double x7 = 0.;
    double x8 = 0.;
    double x9 = 0.;
    double x10 = 0.;
    double x11 = 0.;
    double x12 = 0.;
    int x13 = 0;

    x3 = rhot1*ut1;
    x4 = -rt1 + x3;
    x5 = x4*x4;
    /*@ assert x5 >= 0; */
    x6 = rhot2*ut2;
    x7 = -rt2 + x6;
    x8 = x7*x7;
    /*@ assert x8 >= 0; */
    x9 = x5 + x8;
    /*@ assert x9 >= 0; */
    x10 = sqrt(x9);
    x11 = Max(0, rn);
    /*@ assert x11 >= 0; */
    x12 = Max(0.0000000000000002220446049250313080847263336181640625, mu*x11);
    /*@ assert x12 >= 0; */
    /*@ assert x12 != 0; */
    x13 = x10 <= x12;
    /*@ assert x13 <==> (x10 <= x12); */

    int x16 = 0;

    double x14 = 0.;
    double x15 = 0.;

    x16 = x10 > x12;
    /*@ assert x16 <==> (x10 > x12); */

    if (x13)
    {
    }
    if (x16)
    {
        /*@ assert x10 < -1.09476442525e-47 || x10 > 1.09476442525e-47; */
        x14 = 1.0/x10;
        x15 = x12*x14;

    }double x2 = 0.;
    x2 = Max(0, -rhon*un + rn);
    /*@ assert x2 >= 0; */
    /*@ assigns result[0]; */
    result[0] = rn - x2;


    /*@ assert x13 || x16; */
    if (x13)
    {
        /*@ assigns result[1]; */
        result[1] = x3;

    }
    else if (x16)
    {
        /*@ assigns result[1]; */
        result[1] = rt1 + x15*x4;

    }


    if (x13)
    {
        /*@ assigns result[2]; */
        result[2] = x6;

    }
    else if (x16)
    {
        /*@ assigns result[2]; */
        result[2] = rt2 + x15*x7;

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
    fc3d_AlartCurnierJeanMoreauFGenerated(rn, rt1, rt2, un, ut1, ut2, mu, rhon, rhot1, rhot2, result);
    return(0);
}
#endif