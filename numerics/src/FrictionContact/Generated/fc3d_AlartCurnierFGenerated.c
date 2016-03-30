#include "fc3d_AlartCurnierFGenerated.h"
#include "funcodegen.h"
/*@
requires (0.0 <= rn <= 1.0e+6) && (-1.0e+6 <= rt1 <= 1.0e+6) && (-1.0e+6 <= rt2 <= 1.0e+6) && (-1.0e+6 <= un <= 1.0e+6) && (-1.0e+6 <= ut1 <= 1.0e+6) && (-1.0e+6 <= ut2 <= 1.0e+6) && (0.0 <= mu <= 1.0) && (-1.0e+6 <= rhon <= 1.0e+6) && (-1.0e+6 <= rhot1 <= 1.0e+6) && (-1.0e+6 <= rhot2 <= 1.0e+6);
assigns result[0..2];
ensures \is_finite((double) result[0]) && \is_finite((double) result[1]) && \is_finite((double) result[2]);
*/
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
    /*@
    assert \is_finite((double) un) && \is_finite((double) rn) && \is_finite((double) rhon) && \is_finite((double) rhot2) && \is_finite((double) rhot1) && \is_finite((double) ut1) && \is_finite((double) rt1) && \is_finite((double) mu) && \is_finite((double) ut2) && \is_finite((double) rt2);
    */
    double x1 = 0.;
    double x2 = 0.;
    double x3 = 0.;
    int x4 = 0;
    double x5 = 0.;
    double x6 = 0.;
    double x7 = 0.;
    double x8 = 0.;
    int x9 = 0;

    x1 = Max(0, -rhon*un + rn);
    /*@ assert \is_finite((double) (x1)) && (x1) >= 0; */
    x2 = rhot1*ut1;
    /*@ assert \is_finite((double) (x2)); */
    x3 = mu*x1;
    /*@ assert \is_finite((double) (x3)); */
    x4 = x3 > 0.0000000000000002220446049250313080847263336181640625;
    /*@ assert x4 <==> (x3 > 0.0000000000000002220446049250313080847263336181640625); */

    /*@ assert \is_finite((double) (-rt1 + x2)); */
    x5 = (-rt1 + x2)*(-rt1 + x2);
    /*@ assert \is_finite((double) (x5)) && (x5) >= 0; */
    x6 = rhot2*ut2;
    /*@ assert \is_finite((double) (x6)); */

    /*@ assert \is_finite((double) (-rt2 + x6)); */
    x7 = (-rt2 + x6)*(-rt2 + x6);
    /*@ assert \is_finite((double) (x7)) && (x7) >= 0; */

    /*@ assert \is_finite((double) (x5 + x7)); */
    /*@ assert x5 + x7 >= 0; */
    x8 = sqrt(x5 + x7);
    /*@ assert \is_finite((double) (x8)) && (x8) >= 0; */
    x9 = x4 && x8 <= x3;
    /*@ assert x9 <==> (x4 && x8 <= x3); */

    int x12 = 0;

    double x10 = 0.;
    double x11 = 0.;

    x12 = x4 && x8 > x3;
    /*@ assert x12 <==> (x4 && x8 > x3); */

    int x13 = 0;
    int x14 = 0;

    x13 = x3 <= 0.0000000000000002220446049250313080847263336181640625;
    /*@ assert x13 <==> (x3 <= 0.0000000000000002220446049250313080847263336181640625); */
    x14 = x13 && x8 <= 0;
    /*@ assert x14 <==> (x13 && x8 <= 0); */

    int x15 = 0;

    x15 = x13 && x8 > 0;
    /*@ assert x15 <==> (x13 && x8 > 0); */

    /*@ assert x9 || x12 || x14 || x15; */
    if (x9)
    {
    }
    else if (x12)
    {

        /*@ assert \is_finite((double) (x8)); */
        /*@ assert x8 < -1.09476442525e-47 || x8 > 1.09476442525e-47; */
        x10 = 1.0/x8;
        /*@ assert \is_finite((double) (x10)); */
        x11 = 1.0*mu*x1*x10;
        /*@ assert \is_finite((double) (x11)); */

    }
    else if (x14)
    {
    }/*@ assigns result[0]; */
    result[0] = rn - x1;
    /*@ assert \is_finite((double) (result[0])); */

    /*@ assert x9 || x12 || x14 || x15; */
    if (x9)
    {
        /*@ assigns result[1]; */
        result[1] = x2;
        /*@ assert \is_finite((double) (result[1])); */
    }
    else if (x12)
    {
        /*@ assigns result[1]; */
        result[1] = rt1 - x11*(rt1 - x2);
        /*@ assert \is_finite((double) (result[1])); */
    }
    else if (x14)
    {
        /*@ assigns result[1]; */
        result[1] = x2;
        /*@ assert \is_finite((double) (result[1])); */
    }
    else if (x15)
    {
        /*@ assigns result[1]; */
        result[1] = rt1;
        /*@ assert \is_finite((double) (result[1])); */
    }
    /*@ assert \is_finite((double) (result[1])); */

    /*@ assert x9 || x12 || x14 || x15; */
    if (x9)
    {
        /*@ assigns result[2]; */
        result[2] = x6;
        /*@ assert \is_finite((double) (result[2])); */
    }
    else if (x12)
    {
        /*@ assigns result[2]; */
        result[2] = rt2 - x11*(rt2 - x6);
        /*@ assert \is_finite((double) (result[2])); */
    }
    else if (x14)
    {
        /*@ assigns result[2]; */
        result[2] = x6;
        /*@ assert \is_finite((double) (result[2])); */
    }
    else if (x15)
    {
        /*@ assigns result[2]; */
        result[2] = rt2;
        /*@ assert \is_finite((double) (result[2])); */
    }
    /*@ assert \is_finite((double) (result[2])); */
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
    fc3d_AlartCurnierFGenerated(rn, rt1, rt2, un, ut1, ut2, mu, rhon, rhot1, rhot2, result);
}
#endif