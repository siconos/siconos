#include "fc3d_AlartCurnierJeanMoreauFABGenerated.h"
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
assigns result[0..20];
ensures \is_finite((double) result[0]);
ensures \is_finite((double) result[1]);
ensures \is_finite((double) result[2]);
ensures \is_finite((double) result[3]);
ensures \is_finite((double) result[4]);
ensures \is_finite((double) result[5]);
ensures \is_finite((double) result[6]);
ensures \is_finite((double) result[7]);
ensures \is_finite((double) result[8]);
ensures \is_finite((double) result[9]);
ensures \is_finite((double) result[10]);
ensures \is_finite((double) result[11]);
ensures \is_finite((double) result[12]);
ensures \is_finite((double) result[13]);
ensures \is_finite((double) result[14]);
ensures \is_finite((double) result[15]);
ensures \is_finite((double) result[16]);
ensures \is_finite((double) result[17]);
ensures \is_finite((double) result[18]);
ensures \is_finite((double) result[19]);
ensures \is_finite((double) result[20]);*/
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
    double x5 = 0.;
    double x6 = 0.;
    double x7 = 0.;
    double x8 = 0.;
    double x9 = 0.;
    double x10 = 0.;
    double x11 = 0.;
    double x12 = 0.;
    double x13 = 0.;
    double x14 = 0.;
    double x15 = 0.;
    int x16 = 0;

    x5 = rhot1*ut1;
    x6 = -rt1 + x5;
    x7 = x6*x6;
    /*@ assert x7 >= 0; */
    x8 = rhot2*ut2;
    x9 = -rt2 + x8;
    x10 = x9*x9;
    /*@ assert x10 >= 0; */
    x11 = x10 + x7;
    /*@ assert x11 >= 0; */
    x12 = sqrt(x11);
    x13 = Max(0, rn);
    /*@ assert x13 >= 0; */
    x14 = mu*x13;
    x15 = Max(0.0000000000000002220446049250313080847263336181640625, x14);
    /*@ assert x15 >= 0; */
    /*@ assert x15 != 0; */
    x16 = x12 <= x15;
    /*@ assert x16 <==> (x12 <= x15); */

    int x20 = 0;

    double x17 = 0.;
    double x18 = 0.;
    double x19 = 0.;
    double x21 = 0.;
    double x22 = 0.;
    double x23 = 0.;
    double x24 = 0.;
    double x25 = 0.;
    double x26 = 0.;
    double x27 = 0.;
    double x28 = 0.;
    double x29 = 0.;
    double x30 = 0.;
    double x31 = 0.;
    double x32 = 0.;

    x20 = x12 > x15;
    /*@ assert x20 <==> (x12 > x15); */

    if (x16)
    {
    }
    if (x20)
    {
        x17 = rt1 - x5;
        /*@ assert x12 < -1.09476442525e-47 || x12 > 1.09476442525e-47; */
        x18 = 1.0/x12;
        x19 = x15*x18;
        /*@ assert x11 < -1.09476442525e-47 || x11 > 1.09476442525e-47; */
        x21 = 1.0/((sqrt(x11))*(sqrt(x11))*(sqrt(x11)));
        x22 = x15*x17*x21;
        x23 = x22*x6;
        x24 = x22*x9;
        x25 = Heaviside(rn);
        x26 = Heaviside(x14 - 0.0000000000000002220446049250313080847263336181640625);
        x27 = mu*x18*x25*x26;
        x28 = -x19 + 1;
        x29 = rt2 - x8;
        x30 = x15*x21*x29;
        x31 = x30*x6;
        x32 = x30*x9;

    }double x2 = 0.;
    double x3 = 0.;
    x2 = -rhon*un + rn;
    x3 = Max(0, x2);
    /*@ assert x3 >= 0; */
    /*@ assigns result[0]; */
    result[0] = rn - x3;


    /*@ assert x16 || x20; */
    if (x16)
    {
        /*@ assigns result[1]; */
        result[1] = x5;

    }
    else if (x20)
    {
        /*@ assigns result[1]; */
        result[1] = rt1 - x17*x19;

    }


    if (x16)
    {
        /*@ assigns result[2]; */
        result[2] = x8;

    }
    else if (x20)
    {
        /*@ assigns result[2]; */
        result[2] = rt2 - x19*x29;

    }

    double x4 = 0.;
    x4 = Heaviside(x2);
    /*@ assigns result[3]; */
    result[3] = rhon*x4;

    /*@ assigns result[4]; */
    result[4] = 0;
    /*@ assert result[4] >= 0; */
    /*@ assigns result[5]; */
    result[5] = 0;
    /*@ assert result[5] >= 0; */
    /*@ assigns result[6]; */
    result[6] = 0;
    /*@ assert result[6] >= 0; */

    if (x16)
    {
        /*@ assigns result[7]; */
        result[7] = rhot1;

    }
    else if (x20)
    {
        /*@ assigns result[7]; */
        result[7] = rhot1*x19 + rhot1*x23;

    }


    if (x16)
    {
        /*@ assigns result[8]; */
        result[8] = 0;
        /*@ assert result[8] >= 0; */
    }
    else if (x20)
    {
        /*@ assigns result[8]; */
        result[8] = rhot1*x31;

    }

    /*@ assigns result[9]; */
    result[9] = 0;
    /*@ assert result[9] >= 0; */

    if (x16)
    {
        /*@ assigns result[10]; */
        result[10] = 0;
        /*@ assert result[10] >= 0; */
    }
    else if (x20)
    {
        /*@ assigns result[10]; */
        result[10] = rhot2*x24;

    }


    if (x16)
    {
        /*@ assigns result[11]; */
        result[11] = rhot2;

    }
    else if (x20)
    {
        /*@ assigns result[11]; */
        result[11] = rhot2*x19 + rhot2*x32;

    }

    /*@ assigns result[12]; */
    result[12] = -x4 + 1;


    if (x16)
    {
        /*@ assigns result[13]; */
        result[13] = 0;
        /*@ assert result[13] >= 0; */
    }
    else if (x20)
    {
        /*@ assigns result[13]; */
        result[13] = -x17*x27;

    }


    if (x16)
    {
        /*@ assigns result[14]; */
        result[14] = 0;
        /*@ assert result[14] >= 0; */
    }
    else if (x20)
    {
        /*@ assigns result[14]; */
        result[14] = -x27*x29;

    }

    /*@ assigns result[15]; */
    result[15] = 0;
    /*@ assert result[15] >= 0; */

    if (x16)
    {
        /*@ assigns result[16]; */
        result[16] = 0;
        /*@ assert result[16] >= 0; */
    }
    else if (x20)
    {
        /*@ assigns result[16]; */
        result[16] = -x23 + x28;

    }


    if (x16)
    {
        /*@ assigns result[17]; */
        result[17] = 0;
        /*@ assert result[17] >= 0; */
    }
    else if (x20)
    {
        /*@ assigns result[17]; */
        result[17] = -x31;

    }

    /*@ assigns result[18]; */
    result[18] = 0;
    /*@ assert result[18] >= 0; */

    if (x16)
    {
        /*@ assigns result[19]; */
        result[19] = 0;
        /*@ assert result[19] >= 0; */
    }
    else if (x20)
    {
        /*@ assigns result[19]; */
        result[19] = -x24;

    }


    if (x16)
    {
        /*@ assigns result[20]; */
        result[20] = 0;
        /*@ assert result[20] >= 0; */
    }
    else if (x20)
    {
        /*@ assigns result[20]; */
        result[20] = x28 - x32;

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
    double result[21];
    fc3d_AlartCurnierJeanMoreauFABGenerated(rn, rt1, rt2, un, ut1, ut2, mu, rhon, rhot1, rhot2, result);
    return(0);
}
#endif