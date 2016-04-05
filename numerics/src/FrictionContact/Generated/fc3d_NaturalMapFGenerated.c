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
    /*@ assert \is_finite((double) rn); */
    /*@ assert \is_finite((double) ut1); */
    /*@ assert \is_finite((double) mu); */
    /*@ assert \is_finite((double) un); */
    /*@ assert \is_finite((double) rt1); */
    /*@ assert \is_finite((double) ut2); */
    /*@ assert \is_finite((double) rt2); */
    double x8 = 0.;
    double x9 = 0.;
    double x10 = 0.;
    double x11 = 0.;
    double x12 = 0.;
    double x13 = 0.;
    double x14 = 0.;
    double x15 = 0.;
    double x16 = 0.;
    int x30 = 0;

    double x1 = 0.;
    double x2 = 0.;
    double x3 = 0.;
    double x4 = 0.;
    double x5 = 0.;
    double x6 = 0.;
    double x7 = 0.;
    double x17 = 0.;
    double x19 = 0.;
    double x20 = 0.;
    double x23 = 0.;
    double x24 = 0.;
    double x25 = 0.;
    double x26 = 0.;
    double x27 = 0.;
    double x28 = 0.;
    double x29 = 0.;
    double x32 = 0.;
    double x33 = 0.;
    double x34 = 0.;
    double x35 = 0.;

    x8 = -rt1;
    /*@ assert \is_finite((double) (x8)); */
    x9 = mu*ut1;
    /*@ assert \is_finite((double) (x9)); */
    x10 = x8 + x9;
    /*@ assert \is_finite((double) (x10)); */

    /*@ assert (\is_finite((double) (x10))); */
    x11 = x10*x10;
    /*@ assert \is_finite((double) (x11)); */
    /*@ assert (x11) >= 0; */
    x12 = -rt2;
    /*@ assert \is_finite((double) (x12)); */
    x13 = mu*ut2;
    /*@ assert \is_finite((double) (x13)); */
    x14 = x12 + x13;
    /*@ assert \is_finite((double) (x14)); */

    /*@ assert (\is_finite((double) (x14))); */
    x15 = x14*x14;
    /*@ assert \is_finite((double) (x15)); */
    /*@ assert (x15) >= 0; */

    /*@ assert (\is_finite((double) (x11 + x15))); */
    /*@ assert (x11 + x15 >= 0); */
    x16 = sqrt(x11 + x15);
    /*@ assert \is_finite((double) (x16)); */
    /*@ assert (x16) >= 0; */
    /*@ assert (x16) > 2.22044604925e-16 ==> x11 + x15 > 4.930380657631323783823303533017413935457540219431394e-32; */
    x30 = x16 > 0.0000000000000002220446049250313080847263336181640625;
    /*@ assert x30 <==> (x16 > 0.0000000000000002220446049250313080847263336181640625); */

    int x31 = 0;

    double x18 = 0.;
    double x21 = 0.;

    x31 = x16 <= 0.0000000000000002220446049250313080847263336181640625;
    /*@ assert x31 <==> (x16 <= 0.0000000000000002220446049250313080847263336181640625); */

    if (x30)
    {
        x1 = mu*rn;
        /*@ assert \is_finite((double) (x1)); */
        x2 = -un;
        /*@ assert \is_finite((double) (x2)); */

        /*@ assert (\is_finite((double) (ut1))); */
        x3 = ut1*ut1;
        /*@ assert \is_finite((double) (x3)); */
        /*@ assert (x3) >= 0; */

        /*@ assert (\is_finite((double) (ut2))); */
        x4 = ut2*ut2;
        /*@ assert \is_finite((double) (x4)); */
        /*@ assert (x4) >= 0; */

        /*@ assert (\is_finite((double) (x3 + x4))); */
        /*@ assert (x3 + x4 >= 0); */
        x5 = sqrt(x3 + x4);
        /*@ assert \is_finite((double) (x5)); */
        /*@ assert (x5) >= 0; */
        /*@ assert (x5) > 2.22044604925e-16 ==> x3 + x4 > 4.930380657631323783823303533017413935457540219431394e-32; */
        x6 = -mu*x5;
        /*@ assert \is_finite((double) (x6)); */
        x7 = x1 + x2 + x6;
        /*@ assert \is_finite((double) (x7)); */
        x17 = Max(0, x16 + x7);
        /*@ assert \is_finite((double) (x17)); */
        /*@ assert (x17) >= 0; */
        x19 = -x16;
        /*@ assert \is_finite((double) (x19)); */
        x20 = Max(0, x19 + x7);
        /*@ assert \is_finite((double) (x20)); */
        /*@ assert (x20) >= 0; */
        x23 = -mu*ut1;
        /*@ assert \is_finite((double) (x23)); */
        x24 = rt1 + x23;
        /*@ assert \is_finite((double) (x24)); */

        /*@ assert (\is_finite((double) (x16))); */
        /*@ assert (x16 < -1.09476442525e-47 || x16 > 1.09476442525e-47); */
        x25 = 1.0/x16;
        /*@ assert \is_finite((double) (x25)); */
        x26 = 0.5*x17*x25;
        /*@ assert \is_finite((double) (x26)); */
        x27 = -x24*x26;
        /*@ assert \is_finite((double) (x27)); */
        x28 = 0.5*x20*x25;
        /*@ assert \is_finite((double) (x28)); */
        x29 = -x10*x28;
        /*@ assert \is_finite((double) (x29)); */
        x32 = -mu*ut2;
        /*@ assert \is_finite((double) (x32)); */
        x33 = rt2 + x32;
        /*@ assert \is_finite((double) (x33)); */
        x34 = -x26*x33;
        /*@ assert \is_finite((double) (x34)); */
        x35 = -x14*x28;
        /*@ assert \is_finite((double) (x35)); */

    }
    if (x31)
    {
        x1 = mu*rn;
        /*@ assert \is_finite((double) (x1)); */
        x2 = -un;
        /*@ assert \is_finite((double) (x2)); */

        /*@ assert (\is_finite((double) (ut1))); */
        x3 = ut1*ut1;
        /*@ assert \is_finite((double) (x3)); */
        /*@ assert (x3) >= 0; */

        /*@ assert (\is_finite((double) (ut2))); */
        x4 = ut2*ut2;
        /*@ assert \is_finite((double) (x4)); */
        /*@ assert (x4) >= 0; */

        /*@ assert (\is_finite((double) (x3 + x4))); */
        /*@ assert (x3 + x4 >= 0); */
        x5 = sqrt(x3 + x4);
        /*@ assert \is_finite((double) (x5)); */
        /*@ assert (x5) >= 0; */
        /*@ assert (x5) > 2.22044604925e-16 ==> x3 + x4 > 4.930380657631323783823303533017413935457540219431394e-32; */
        x6 = -mu*x5;
        /*@ assert \is_finite((double) (x6)); */
        x7 = x1 + x2 + x6;
        /*@ assert \is_finite((double) (x7)); */
        x17 = Max(0, x16 + x7);
        /*@ assert \is_finite((double) (x17)); */
        /*@ assert (x17) >= 0; */
        x18 = -0.5*x17;
        /*@ assert \is_finite((double) (x18)); */
        x19 = -x16;
        /*@ assert \is_finite((double) (x19)); */
        x20 = Max(0, x19 + x7);
        /*@ assert \is_finite((double) (x20)); */
        /*@ assert (x20) >= 0; */
        x21 = 0.5*x20;
        /*@ assert \is_finite((double) (x21)); */

    }double x22 = 0.;
    x1 = mu*rn;
    /*@ assert \is_finite((double) (x1)); */
    x2 = -un;
    /*@ assert \is_finite((double) (x2)); */

    /*@ assert (\is_finite((double) (ut1))); */
    x3 = ut1*ut1;
    /*@ assert \is_finite((double) (x3)); */
    /*@ assert (x3) >= 0; */

    /*@ assert (\is_finite((double) (ut2))); */
    x4 = ut2*ut2;
    /*@ assert \is_finite((double) (x4)); */
    /*@ assert (x4) >= 0; */

    /*@ assert (\is_finite((double) (x3 + x4))); */
    /*@ assert (x3 + x4 >= 0); */
    x5 = sqrt(x3 + x4);
    /*@ assert \is_finite((double) (x5)); */
    /*@ assert (x5) >= 0; */
    /*@ assert (x5) > 2.22044604925e-16 ==> x3 + x4 > 4.930380657631323783823303533017413935457540219431394e-32; */
    x6 = -mu*x5;
    /*@ assert \is_finite((double) (x6)); */
    x7 = x1 + x2 + x6;
    /*@ assert \is_finite((double) (x7)); */
    x17 = Max(0, x16 + x7);
    /*@ assert \is_finite((double) (x17)); */
    /*@ assert (x17) >= 0; */
    x18 = -0.5*x17;
    /*@ assert \is_finite((double) (x18)); */
    x19 = -x16;
    /*@ assert \is_finite((double) (x19)); */
    x20 = Max(0, x19 + x7);
    /*@ assert \is_finite((double) (x20)); */
    /*@ assert (x20) >= 0; */
    x21 = 0.5*x20;
    /*@ assert \is_finite((double) (x21)); */
    x22 = -x21;
    /*@ assert \is_finite((double) (x22)); */
    /*@ assigns result[0]; */
    result[0] = x1 + x18 + x22;
    /*@ assert \is_finite((double) (result[0])); */

    /*@ assert x30 || x31; */
    if (x30)
    {
        /*@ assigns result[1]; */
        result[1] = rt1 + x27 + x29;
        /*@ assert \is_finite((double) (result[1])); */
    }
    else if (x31)
    {
        /*@ assigns result[1]; */
        result[1] = rt1;
        /*@ assert \is_finite((double) (result[1])); */
    }
    /*@ assert \is_finite((double) (result[1])); */

    /*@ assert x30 || x31; */
    if (x30)
    {
        /*@ assigns result[2]; */
        result[2] = rt2 + x34 + x35;
        /*@ assert \is_finite((double) (result[2])); */
    }
    else if (x31)
    {
        /*@ assigns result[2]; */
        result[2] = rt2 + x18 + x21;
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
    fc3d_NaturalMapFGenerated(rn, rt1, rt2, un, ut1, ut2, mu, rhon, rhot1, rhot2, result);
}
#endif