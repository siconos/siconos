#include "fc3d_AlartCurnierJeanMoreauABGenerated.h"
#include "funcodegen.h"
/*@
requires (0.0 <= rn <= 1.0e+6) && (-1.0e+6 <= rt1 <= 1.0e+6) && (-1.0e+6 <= rt2 <= 1.0e+6) && (-1.0e+6 <= un <= 1.0e+6) && (-1.0e+6 <= ut1 <= 1.0e+6) && (-1.0e+6 <= ut2 <= 1.0e+6) && (0.0 <= mu <= 1.0) && (-1.0e+6 <= rhon <= 1.0e+6) && (-1.0e+6 <= rhot1 <= 1.0e+6) && (-1.0e+6 <= rhot2 <= 1.0e+6);
assigns result[0..17];
ensures \is_finite((double) result[0]) && \is_finite((double) result[1]) && \is_finite((double) result[2]) && \is_finite((double) result[3]) && \is_finite((double) result[4]) && \is_finite((double) result[5]) && \is_finite((double) result[6]) && \is_finite((double) result[7]) && \is_finite((double) result[8]) && \is_finite((double) result[9]) && \is_finite((double) result[10]) && \is_finite((double) result[11]) && \is_finite((double) result[12]) && \is_finite((double) result[13]) && \is_finite((double) result[14]) && \is_finite((double) result[15]) && \is_finite((double) result[16]) && \is_finite((double) result[17]);
*/
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
    /*@
    assert \is_finite((double) un) && \is_finite((double) rn) && \is_finite((double) rhon) && \is_finite((double) rhot2) && \is_finite((double) rhot1) && \is_finite((double) ut1) && \is_finite((double) rt1) && \is_finite((double) mu) && \is_finite((double) ut2) && \is_finite((double) rt2);
    */
    double x2 = 0.;
    double x3 = 0.;
    int x4 = 0;
    double x5 = 0.;
    double x6 = 0.;
    double x7 = 0.;
    double x8 = 0.;
    double x9 = 0.;
    double x10 = 0.;
    double x11 = 0.;
    double x12 = 0.;
    int x13 = 0;

    x2 = Max(0, rn);
    /*@ assert \is_finite((double) (x2)) && (x2) >= 0; */
    x3 = mu*x2;
    /*@ assert \is_finite((double) (x3)); */
    x4 = x3 > 0.0000000000000002220446049250313080847263336181640625;
    /*@ assert x4 <==> (x3 > 0.0000000000000002220446049250313080847263336181640625); */
    x5 = rhot1*ut1;
    /*@ assert \is_finite((double) (x5)); */
    x6 = -rt1 + x5;
    /*@ assert \is_finite((double) (x6)); */

    /*@ assert \is_finite((double) (x6)); */
    x7 = x6*x6;
    /*@ assert \is_finite((double) (x7)) && (x7) >= 0; */
    x8 = rhot2*ut2;
    /*@ assert \is_finite((double) (x8)); */
    x9 = -rt2 + x8;
    /*@ assert \is_finite((double) (x9)); */

    /*@ assert \is_finite((double) (x9)); */
    x10 = x9*x9;
    /*@ assert \is_finite((double) (x10)) && (x10) >= 0; */
    x11 = x10 + x7;
    /*@ assert \is_finite((double) (x11)); */

    /*@ assert \is_finite((double) (x11)); */
    /*@ assert x11 >= 0; */
    x12 = sqrt(x11);
    /*@ assert \is_finite((double) (x12)) && (x12) >= 0; */
    x13 = x4 && x12 <= x3;
    /*@ assert x13 <==> (x4 && x12 <= x3); */

    int x19 = 0;

    double x14 = 0.;
    double x15 = 0.;
    double x16 = 0.;
    double x17 = 0.;
    double x18 = 0.;

    x19 = x4 && x12 > x3;
    /*@ assert x19 <==> (x4 && x12 > x3); */

    int x20 = 0;
    int x21 = 0;

    x20 = x3 <= 0.0000000000000002220446049250313080847263336181640625;
    /*@ assert x20 <==> (x3 <= 0.0000000000000002220446049250313080847263336181640625); */
    x21 = x20 && x12 <= 0;
    /*@ assert x21 <==> (x20 && x12 <= 0); */

    int x22 = 0;

    x22 = x20 && x12 > 0;
    /*@ assert x22 <==> (x20 && x12 > 0); */

    /*@ assert x13 || x19 || x21 || x22 || x20; */
    if (x19)
    {

        /*@ assert \is_finite((double) (x12)); */
        /*@ assert x12 < -1.09476442525e-47 || x12 > 1.09476442525e-47; */
        x14 = 1.0/x12;
        /*@ assert \is_finite((double) (x14)); */
        x15 = 1.0*x14*x3;
        /*@ assert \is_finite((double) (x15)); */
        x16 = rt1 - x5;
        /*@ assert \is_finite((double) (x16)); */

        /*@ assert \is_finite((double) (x11)); */
        /*@ assert x11 < -1.09476442525e-47 || x11 > 1.09476442525e-47; */
        x17 = 1.0/x11;
        /*@ assert \is_finite((double) (x17)); */
        x18 = 1.0*mu*rhot1*x14*x17*x2*x6;
        /*@ assert \is_finite((double) (x18)); */

    }double x1 = 0.;
    x1 = Heaviside(-rhon*un + rn);
    /*@ assert \is_finite((double) (x1)); */
    /*@ assigns result[0]; */
    result[0] = rhon*x1;
    /*@ assert \is_finite((double) (result[0])); */
    /*@ assigns result[1]; */
    result[1] = 0;
    /*@ assert \is_finite((double) (result[1])) && (result[1]) >= 0; */
    /*@ assigns result[2]; */
    result[2] = 0;
    /*@ assert \is_finite((double) (result[2])) && (result[2]) >= 0; */
    /*@ assigns result[3]; */
    result[3] = 0;
    /*@ assert \is_finite((double) (result[3])) && (result[3]) >= 0; */

    /*@ assert x13 || x19 || x21 || x22; */
    if (x13)
    {
        /*@ assigns result[4]; */
        result[4] = rhot1;
        /*@ assert \is_finite((double) (result[4])); */
    }
    else if (x19)
    {
        /*@ assigns result[4]; */
        result[4] = rhot1*x15 + x16*x18;
        /*@ assert \is_finite((double) (result[4])); */
    }
    else if (x21)
    {
        /*@ assigns result[4]; */
        result[4] = rhot1;
        /*@ assert \is_finite((double) (result[4])); */
    }
    else if (x22)
    {
        /*@ assigns result[4]; */
        result[4] = 0;
        /*@ assert \is_finite((double) (result[4])) && (result[4]) >= 0; */
    }
    /*@ assert \is_finite((double) (result[4])); */
    double x28 = 0.;

    /*@ assert x13 || x19 || x20; */
    if (x13)
    {
        /*@ assigns result[5]; */
        result[5] = 0;
        /*@ assert \is_finite((double) (result[5])) && (result[5]) >= 0; */
    }
    else if (x19)
    {
        x28 = rt2 - x8;
        /*@ assert \is_finite((double) (x28)); */

        /*@ assigns result[5]; */
        result[5] = x18*x28;
        /*@ assert \is_finite((double) (result[5])); */
    }
    else if (x20)
    {
        /*@ assigns result[5]; */
        result[5] = 0;
        /*@ assert \is_finite((double) (result[5])) && (result[5]) >= 0; */
    }
    /*@ assert \is_finite((double) (result[5])); */
    /*@ assigns result[6]; */
    result[6] = 0;
    /*@ assert \is_finite((double) (result[6])) && (result[6]) >= 0; */
    double x23 = 0.;

    /*@ assert x13 || x19 || x20; */
    if (x13)
    {
        /*@ assigns result[7]; */
        result[7] = 0;
        /*@ assert \is_finite((double) (result[7])) && (result[7]) >= 0; */
    }
    else if (x19)
    {
        x23 = 1.0*mu*rhot2*x14*x17*x2*x9;
        /*@ assert \is_finite((double) (x23)); */

        /*@ assigns result[7]; */
        result[7] = x16*x23;
        /*@ assert \is_finite((double) (result[7])); */
    }
    else if (x20)
    {
        /*@ assigns result[7]; */
        result[7] = 0;
        /*@ assert \is_finite((double) (result[7])) && (result[7]) >= 0; */
    }
    /*@ assert \is_finite((double) (result[7])); */

    /*@ assert x13 || x19 || x21 || x22; */
    if (x13)
    {
        /*@ assigns result[8]; */
        result[8] = rhot2;
        /*@ assert \is_finite((double) (result[8])); */
    }
    else if (x19)
    {
        /*@ assigns result[8]; */
        result[8] = rhot2*x15 + x23*x28;
        /*@ assert \is_finite((double) (result[8])); */
    }
    else if (x21)
    {
        /*@ assigns result[8]; */
        result[8] = rhot2;
        /*@ assert \is_finite((double) (result[8])); */
    }
    else if (x22)
    {
        /*@ assigns result[8]; */
        result[8] = 0;
        /*@ assert \is_finite((double) (result[8])) && (result[8]) >= 0; */
    }
    /*@ assert \is_finite((double) (result[8])); */
    /*@ assigns result[9]; */
    result[9] = -x1 + 1;
    /*@ assert \is_finite((double) (result[9])); */
    double x24 = 0.;
    double x25 = 0.;

    /*@ assert x13 || x19 || x20; */
    if (x13)
    {
        /*@ assigns result[10]; */
        result[10] = 0;
        /*@ assert \is_finite((double) (result[10])) && (result[10]) >= 0; */
    }
    else if (x19)
    {
        x24 = Heaviside(rn);
        /*@ assert \is_finite((double) (x24)); */
        x25 = 1.0*mu*x14*x24;
        /*@ assert \is_finite((double) (x25)); */

        /*@ assigns result[10]; */
        result[10] = -x16*x25;
        /*@ assert \is_finite((double) (result[10])); */
    }
    else if (x20)
    {
        /*@ assigns result[10]; */
        result[10] = 0;
        /*@ assert \is_finite((double) (result[10])) && (result[10]) >= 0; */
    }
    /*@ assert \is_finite((double) (result[10])); */

    /*@ assert x13 || x19 || x20; */
    if (x13)
    {
        /*@ assigns result[11]; */
        result[11] = 0;
        /*@ assert \is_finite((double) (result[11])) && (result[11]) >= 0; */
    }
    else if (x19)
    {
        /*@ assigns result[11]; */
        result[11] = -x25*x28;
        /*@ assert \is_finite((double) (result[11])); */
    }
    else if (x20)
    {
        /*@ assigns result[11]; */
        result[11] = 0;
        /*@ assert \is_finite((double) (result[11])) && (result[11]) >= 0; */
    }
    /*@ assert \is_finite((double) (result[11])); */
    /*@ assigns result[12]; */
    result[12] = 0;
    /*@ assert \is_finite((double) (result[12])) && (result[12]) >= 0; */
    double x26 = 0.;
    double x27 = 0.;

    /*@ assert x13 || x19 || x21 || x22; */
    if (x13)
    {
        /*@ assigns result[13]; */
        result[13] = 0;
        /*@ assert \is_finite((double) (result[13])) && (result[13]) >= 0; */
    }
    else if (x19)
    {
        x26 = -x15 + 1;
        /*@ assert \is_finite((double) (x26)); */
        x27 = 1.0*mu*x14*x16*x17*x2;
        /*@ assert \is_finite((double) (x27)); */

        /*@ assigns result[13]; */
        result[13] = x26 - x27*x6;
        /*@ assert \is_finite((double) (result[13])); */
    }
    else if (x21)
    {
        /*@ assigns result[13]; */
        result[13] = 0;
        /*@ assert \is_finite((double) (result[13])) && (result[13]) >= 0; */
    }
    else if (x22)
    {
        /*@ assigns result[13]; */
        result[13] = 1;
        /*@ assert \is_finite((double) (result[13])) && (result[13]) >= 0 && (result[13]) != 0; */
    }
    /*@ assert \is_finite((double) (result[13])); */
    double x29 = 0.;

    /*@ assert x13 || x19 || x20; */
    if (x13)
    {
        /*@ assigns result[14]; */
        result[14] = 0;
        /*@ assert \is_finite((double) (result[14])) && (result[14]) >= 0; */
    }
    else if (x19)
    {
        x29 = 1.0*mu*x14*x17*x2*x28;
        /*@ assert \is_finite((double) (x29)); */

        /*@ assigns result[14]; */
        result[14] = -x29*x6;
        /*@ assert \is_finite((double) (result[14])); */
    }
    else if (x20)
    {
        /*@ assigns result[14]; */
        result[14] = 0;
        /*@ assert \is_finite((double) (result[14])) && (result[14]) >= 0; */
    }
    /*@ assert \is_finite((double) (result[14])); */
    /*@ assigns result[15]; */
    result[15] = 0;
    /*@ assert \is_finite((double) (result[15])) && (result[15]) >= 0; */

    /*@ assert x13 || x19 || x20; */
    if (x13)
    {
        /*@ assigns result[16]; */
        result[16] = 0;
        /*@ assert \is_finite((double) (result[16])) && (result[16]) >= 0; */
    }
    else if (x19)
    {
        /*@ assigns result[16]; */
        result[16] = -x27*x9;
        /*@ assert \is_finite((double) (result[16])); */
    }
    else if (x20)
    {
        /*@ assigns result[16]; */
        result[16] = 0;
        /*@ assert \is_finite((double) (result[16])) && (result[16]) >= 0; */
    }
    /*@ assert \is_finite((double) (result[16])); */

    /*@ assert x13 || x19 || x21 || x22; */
    if (x13)
    {
        /*@ assigns result[17]; */
        result[17] = 0;
        /*@ assert \is_finite((double) (result[17])) && (result[17]) >= 0; */
    }
    else if (x19)
    {
        /*@ assigns result[17]; */
        result[17] = x26 - x29*x9;
        /*@ assert \is_finite((double) (result[17])); */
    }
    else if (x21)
    {
        /*@ assigns result[17]; */
        result[17] = 0;
        /*@ assert \is_finite((double) (result[17])) && (result[17]) >= 0; */
    }
    else if (x22)
    {
        /*@ assigns result[17]; */
        result[17] = 1;
        /*@ assert \is_finite((double) (result[17])) && (result[17]) >= 0 && (result[17]) != 0; */
    }
    /*@ assert \is_finite((double) (result[17])); */
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
    double result[18];
    fc3d_AlartCurnierJeanMoreauABGenerated(rn, rt1, rt2, un, ut1, ut2, mu, rhon, rhot1, rhot2, result);
}
#endif