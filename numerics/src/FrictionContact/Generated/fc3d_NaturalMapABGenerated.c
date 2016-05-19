#include "fc3d_NaturalMapABGenerated.h"
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
assigns result[0..17];
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
ensures \is_finite((double) result[17]);*/
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
    /*@ assert \is_finite((double) ut1); */
    /*@ assert \is_finite((double) epsilon); */
    /*@ assert \is_finite((double) un); */
    /*@ assert \is_finite((double) rt1); */
    /*@ assert \is_finite((double) mu); */
    /*@ assert \is_finite((double) ut2); */
    /*@ assert \is_finite((double) rn); */
    /*@ assert \is_finite((double) rt2); */
    double x1 = 0.;
    double x2 = 0.;
    int x29 = 0;
    int x30 = 0;
    int x31 = 0;

    double x4 = 0.;
    double x7 = 0.;
    double x24 = 0.;
    double x25 = 0.;
    double x26 = 0.;
    double x27 = 0.;
    double x28 = 0.;
    double x74 = 0.;
    double x75 = 0.;
    double x76 = 0.;
    double x100 = 0.;
    double x101 = 0.;
    double x121 = 0.;
    double x122 = 0.;
    double x123 = 0.;

    x1 = rt1*rt1 + rt2*rt2;
    /*@ assert x1 >= 0; */
    x2 = ut1*ut1 + ut2*ut2;
    /*@ assert x2 >= 0; */
    x29 = x1 <= epsilon;
    /*@ assert x29 <==> (x1 <= epsilon); */
    x30 = x2 <= epsilon;
    /*@ assert x30 <==> (x2 <= epsilon); */
    x31 = x29 && x30;
    /*@ assert x31 <==> (x29 && x30); */

    int x42 = 0;
    int x43 = 0;

    double x6 = 0.;
    double x32 = 0.;
    double x33 = 0.;
    double x34 = 0.;
    double x35 = 0.;
    double x36 = 0.;
    double x37 = 0.;
    double x38 = 0.;
    double x39 = 0.;
    double x40 = 0.;
    double x41 = 0.;
    double x60 = 0.;
    double x77 = 0.;
    double x78 = 0.;
    double x79 = 0.;
    double x80 = 0.;
    double x81 = 0.;
    double x82 = 0.;
    double x83 = 0.;
    double x84 = 0.;
    double x85 = 0.;
    double x86 = 0.;
    double x87 = 0.;
    double x88 = 0.;
    double x89 = 0.;
    double x90 = 0.;
    double x91 = 0.;
    double x102 = 0.;
    double x103 = 0.;
    double x104 = 0.;
    double x105 = 0.;
    double x106 = 0.;
    double x124 = 0.;

    x42 = x1 > epsilon;
    /*@ assert x42 <==> (x1 > epsilon); */
    x43 = x30 && x42;
    /*@ assert x43 <==> (x30 && x42); */

    double x10 = 0.;
    double x11 = 0.;
    double x12 = 0.;
    double x13 = 0.;
    double x14 = 0.;
    double x15 = 0.;
    double x16 = 0.;
    int x44 = 0;
    double x45 = 0.;
    int x46 = 0;
    int x47 = 0;

    x10 = mu*ut1;
    x11 = rt1 - x10;
    x12 = x11*x11;
    /*@ assert x12 >= 0; */
    x13 = mu*ut2;
    x14 = rt2 - x13;
    x15 = x14*x14;
    /*@ assert x15 >= 0; */
    x16 = x12 + x15;
    x44 = x2 >= epsilon;
    /*@ assert x44 <==> (x2 >= epsilon); */
    x45 = epsilon*(mu + 1);
    /*@ assert x45 >= 0; */
    /*@ assert x45 != 0; */
    x46 = x16 <= x45;
    /*@ assert x46 <==> (x16 <= x45); */
    x47 = x29 && x44 && x46;
    /*@ assert x47 <==> (x29 && x44 && x46); */

    int x49 = 0;

    double x48 = 0.;

    x49 = x42 && x44 && x46;
    /*@ assert x49 <==> (x42 && x44 && x46); */

    int x57 = 0;
    int x58 = 0;
    int x59 = 0;

    double x8 = 0.;
    double x9 = 0.;
    double x17 = 0.;
    double x18 = 0.;
    double x19 = 0.;
    double x20 = 0.;
    double x21 = 0.;
    double x22 = 0.;
    double x23 = 0.;
    double x50 = 0.;
    double x51 = 0.;
    double x52 = 0.;
    double x53 = 0.;
    double x54 = 0.;
    double x55 = 0.;
    double x56 = 0.;
    double x61 = 0.;
    double x62 = 0.;
    double x63 = 0.;
    double x64 = 0.;
    double x67 = 0.;
    double x68 = 0.;
    double x69 = 0.;
    double x71 = 0.;
    double x72 = 0.;
    double x73 = 0.;
    double x92 = 0.;
    double x93 = 0.;
    double x94 = 0.;
    double x95 = 0.;
    double x96 = 0.;
    double x97 = 0.;
    double x98 = 0.;
    double x99 = 0.;
    double x107 = 0.;
    double x108 = 0.;
    double x109 = 0.;
    double x116 = 0.;
    double x118 = 0.;
    double x120 = 0.;
    double x125 = 0.;

    x57 = x44 && x46;
    /*@ assert x57 <==> (x44 && x46); */
    x58 = x30 || x57;
    /*@ assert x58 <==> (x30 || x57); */
    x59 = !x58;
    /*@ assert x59 <==> (!x58); */

    double x65 = 0.;
    double x66 = 0.;

    int x70 = 0;

    double x110 = 0.;
    double x111 = 0.;
    double x112 = 0.;
    double x113 = 0.;
    double x114 = 0.;
    double x115 = 0.;
    double x117 = 0.;
    double x119 = 0.;
    double x126 = 0.;

    x70 = x16 > x45;
    /*@ assert x70 <==> (x16 > x45); */

    if (x31)
    {
        x4 = 2*mu*mu - 2*mu + 1;
        x7 = mu*rn;
        /*@ assert x7 >= 0; */
        /*@ assert 2 >= 0; */
        x24 = sqrt(2);
        /*@ assert x24 >= 0; */
        /*@ assert x24 != 0; */
        x25 = -1.0*un + x7;
        x26 = Heaviside(x25);
        x27 = mu*x26;
        x28 = (1.0/2.0)*x24*x27;
        /*@ assert x4 < -epsilon*epsilon*epsilon || x4 > epsilon*epsilon*epsilon; */
        x74 = 1.0/x4;
        x75 = mu - 1;
        x76 = x75*x75;
        /*@ assert x76 >= 0; */
        x100 = mu*mu;
        /*@ assert x100 >= 0; */
        x101 = x26*x74;
        x121 = -2.0*mu + 2*x100 + 1.0;
        /*@ assert x121 < -epsilon*epsilon*epsilon || x121 > epsilon*epsilon*epsilon; */
        x122 = 1.0/x121;
        x123 = mu*mu*mu;
        /*@ assert x123 >= 0; */

    }
    if (x43)
    {
        x6 = -un;
        x7 = mu*rn;
        /*@ assert 2 >= 0; */
        x24 = sqrt(2);
        x25 = -1.0*un + x7;
        x32 = sqrt(x1);
        /*@ assert x32 < -epsilon*epsilon*epsilon || x32 > epsilon*epsilon*epsilon; */
        x33 = 1.0/x32;
        x34 = 0.25*mu*x33;
        x35 = Heaviside(x25 + x32);
        x36 = 2*x35;
        x37 = rt1*x36;
        x38 = Heaviside(x25 - 1.0*x32);
        x39 = 2.0*x38;
        x40 = x24*x32;
        x41 = x35*x40 + x38*x40;
        x60 = rt2*x36;
        /*@ assert x1 < -epsilon*epsilon*epsilon || x1 > epsilon*epsilon*epsilon; */
        x77 = 1.0/((sqrt(x1))*(sqrt(x1))*(sqrt(x1)));
        x78 = 0.25*mu*x77;
        x79 = rt2*rt2;
        /*@ assert x79 >= 0; */
        x80 = x6 + x7;
        x81 = Max(0, x32 + x80);
        /*@ assert x81 >= 0; */
        x82 = 2*x81;
        x83 = Max(0, -x32 + x80);
        /*@ assert x83 >= 0; */
        x84 = 2.0*x83;
        x85 = rt1*rt1;
        /*@ assert x85 >= 0; */
        x86 = 2*x32*x35;
        x87 = 2*x32;
        x88 = x38*x85;
        x89 = rt1*rt1*rt1;
        x90 = x35*x79;
        x91 = x38*x79;
        x102 = 2.0*x81;
        x103 = 2*x83;
        x104 = 2*x32*x38;
        x105 = x35*x85;
        x106 = x24*(x105 - x88 + x90 - x91);
        x124 = rt2*rt2*rt2;

    }
    if (x47)
    {
        x7 = mu*rn;
        /*@ assert 2 >= 0; */
        x24 = sqrt(2);
        x25 = -1.0*un + x7;
        x26 = Heaviside(x25);
        x27 = mu*x26;
        x28 = (1.0/2.0)*x24*x27;

    }
    if (x49)
    {
        x7 = mu*rn;
        x25 = -1.0*un + x7;
        x32 = sqrt(x1);
        /*@ assert x32 < -epsilon*epsilon*epsilon || x32 > epsilon*epsilon*epsilon; */
        x33 = 1.0/x32;
        x38 = Heaviside(x25 - 1.0*x32);
        x48 = mu*x33*x38;

    }
    if (x59)
    {
        x6 = -un;
        x7 = mu*rn;
        x8 = sqrt(x2);
        x9 = -mu*x8 + x6 + x7;
        /*@ assert x16 >= 0; */
        x17 = sqrt(x16);
        x18 = x17 + x9;
        x19 = Heaviside(x18);
        x20 = 0.5*x19;
        x21 = -x17 + x9;
        x22 = Heaviside(x21);
        x23 = 0.5*x22;
        /*@ assert x8 < -epsilon*epsilon*epsilon || x8 > epsilon*epsilon*epsilon; */
        x50 = 1.0/x8;
        x51 = -x10*x50;
        /*@ assert x17 < -epsilon*epsilon*epsilon || x17 > epsilon*epsilon*epsilon; */
        x52 = 1.0/x17;
        x53 = mu*x52;
        x54 = x11*x53;
        x55 = x51 - x54;
        x56 = x51 + x54;
        x61 = -x13*x50;
        x62 = x14*x53;
        x63 = x61 - x62;
        x64 = x61 + x62;
        x67 = 0.5*x19*x52;
        x68 = x11*x67;
        x69 = 0.5*x22*x52;
        x71 = x14*x67;
        x72 = -rt1 + x10;
        x73 = x69*x72;
        x92 = Max(0, x18);
        /*@ assert x92 >= 0; */
        x93 = 0.5*x92;
        x94 = -x53;
        /*@ assert x16 < -epsilon*epsilon*epsilon || x16 > epsilon*epsilon*epsilon; */
        x95 = 1.0/((sqrt(x16))*(sqrt(x16))*(sqrt(x16)));
        x96 = mu*x95;
        x97 = Max(0, x21);
        /*@ assert x97 >= 0; */
        x98 = 0.5*x97;
        x99 = x11*x72;
        x107 = -0.5*mu*x11*x14*x92*x95;
        x108 = x14*x72;
        x109 = 0.5*mu*x95*x97;
        x116 = -rt2 + x13;
        x118 = x11*x116;
        x120 = x116*x69;
        x125 = x116*x14;

    }
    if (x46)
    {
        x6 = -un;
        x7 = mu*rn;
        x8 = sqrt(x2);
        x9 = -mu*x8 + x6 + x7;
        /*@ assert x16 >= 0; */
        x17 = sqrt(x16);
        x18 = x17 + x9;
        x19 = Heaviside(x18);
        x20 = 0.5*x19;
        x21 = -x17 + x9;
        x22 = Heaviside(x21);
        x23 = 0.5*x22;
        x25 = -1.0*un + x7;
        x32 = sqrt(x1);
        x38 = Heaviside(x25 - 1.0*x32);
        x65 = -mu*x20;
        x66 = mu*x23;

    }
    if (x70)
    {
        x6 = -un;
        x7 = mu*rn;
        x8 = sqrt(x2);
        x9 = -mu*x8 + x6 + x7;
        /*@ assert x16 >= 0; */
        x17 = sqrt(x16);
        x18 = x17 + x9;
        x19 = Heaviside(x18);
        x20 = 0.5*x19;
        x21 = -x17 + x9;
        x22 = Heaviside(x21);
        /*@ assert x17 < -epsilon*epsilon*epsilon || x17 > epsilon*epsilon*epsilon; */
        x52 = 1.0/x17;
        x53 = mu*x52;
        x54 = x11*x53;
        x62 = x14*x53;
        x67 = 0.5*x19*x52;
        x68 = x11*x67;
        x69 = 0.5*x22*x52;
        x71 = x14*x67;
        x72 = -rt1 + x10;
        x73 = x69*x72;
        x92 = Max(0, x18);
        x93 = 0.5*x92;
        /*@ assert x16 < -epsilon*epsilon*epsilon || x16 > epsilon*epsilon*epsilon; */
        x95 = 1.0/((sqrt(x16))*(sqrt(x16))*(sqrt(x16)));
        x97 = Max(0, x21);
        x98 = 0.5*x97;
        x99 = x11*x72;
        x108 = x14*x72;
        x110 = 0.5*mu*x22*x52;
        x111 = 1.0/x16;
        x112 = 0.5*x111*x19;
        x113 = 0.5*x111*x22;
        x114 = -x52;
        x115 = x72*x72;
        /*@ assert x115 >= 0; */
        x116 = -rt2 + x13;
        x117 = -x11*x112*x14 - x116*x72*x95*x98;
        x118 = x11*x116;
        x119 = 0.5*x92*x95;
        x120 = x116*x69;
        x125 = x116*x14;
        x126 = x116*x116;
        /*@ assert x126 >= 0; */

    }
    if (x57)
    {
        x7 = mu*rn;
        x25 = -1.0*un + x7;
        x32 = sqrt(x1);
        x38 = Heaviside(x25 - 1.0*x32);

    }x6 = -un;
    x7 = mu*rn;
    x8 = sqrt(x2);
    x9 = -mu*x8 + x6 + x7;
    /*@ assert x16 >= 0; */
    x17 = sqrt(x16);
    x18 = x17 + x9;
    x19 = Heaviside(x18);
    x20 = 0.5*x19;
    x21 = -x17 + x9;
    x22 = Heaviside(x21);
    x23 = 0.5*x22;
    /*@ assigns result[0]; */
    result[0] = x20 + x23;


    /*@ assert x70 || x46; */
    if (x70)
    {
        /*@ assigns result[1]; */
        result[1] = x68 + x73;

    }
    else if (x46)
    {
        /*@ assigns result[1]; */
        result[1] = 0;
        /*@ assert result[1] >= 0; */
    }


    if (x70)
    {
        /*@ assigns result[2]; */
        result[2] = x120 + x71;

    }
    else if (x46)
    {
        /*@ assigns result[2]; */
        result[2] = x20 - x23;

    }


    /*@ assert x31 || x43 || x47 || x49 || x59; */
    if (x31)
    {
        /*@ assigns result[3]; */
        result[3] = x28;

    }
    else if (x43)
    {
        /*@ assigns result[3]; */
        result[3] = x34*(-rt1*x39 + x37 + x41);

    }
    else if (x47)
    {
        /*@ assigns result[3]; */
        result[3] = x28;

    }
    else if (x49)
    {
        /*@ assigns result[3]; */
        result[3] = rt1*x48;

    }
    else if (x59)
    {
        /*@ assigns result[3]; */
        result[3] = -x20*x55 - x23*x56;

    }


    /*@ assert x31 || x43 || x57 || x59; */
    if (x31)
    {
        /*@ assigns result[4]; */
        result[4] = 1.0*x27*x74*x76;

    }
    else if (x43)
    {
        /*@ assigns result[4]; */
        result[4] = x78*(x24*(rt1*x90 - rt1*x91 + x35*x89 - x38*x89) + x79*x82 - x79*x84 + x85*x86 + x87*x88);

    }
    else if (x57)
    {
        /*@ assigns result[4]; */
        result[4] = 0.0;
        /*@ assert result[4] >= 0; */
    }
    else if (x59)
    {
        /*@ assigns result[4]; */
        result[4] = -x55*x68 - x56*x73 - x93*(x12*x96 + x94) - x98*(-x94 + x96*x99);

    }


    if (x31)
    {
        /*@ assigns result[5]; */
        result[5] = 1.0*x100*x122*x26*x75;

    }
    else if (x43)
    {
        /*@ assigns result[5]; */
        result[5] = rt2*x78*(-rt1*x102 + rt1*x103 + rt1*x104 + x106 + x32*x37);

    }
    else if (x57)
    {
        /*@ assigns result[5]; */
        result[5] = 0.0;
        /*@ assert result[5] >= 0; */
    }
    else if (x59)
    {
        /*@ assigns result[5]; */
        result[5] = x107 - x109*x118 - x120*x56 - x55*x71;

    }


    if (x31)
    {
        /*@ assigns result[6]; */
        result[6] = x28;

    }
    else if (x43)
    {
        /*@ assigns result[6]; */
        result[6] = x34*(-rt2*x39 + x41 + x60);

    }
    else if (x47)
    {
        /*@ assigns result[6]; */
        result[6] = x28;

    }
    else if (x49)
    {
        /*@ assigns result[6]; */
        result[6] = rt2*x48;

    }
    else if (x59)
    {
        /*@ assigns result[6]; */
        result[6] = -x20*x63 - x23*x64;

    }


    if (x31)
    {
        /*@ assigns result[7]; */
        result[7] = x100*x101*(mu - 1.0);

    }
    else if (x43)
    {
        /*@ assigns result[7]; */
        result[7] = rt1*x78*(-rt2*x102 + rt2*x103 + rt2*x104 + x106 + x32*x60);

    }
    else if (x57)
    {
        /*@ assigns result[7]; */
        result[7] = 0.0;
        /*@ assert result[7] >= 0; */
    }
    else if (x59)
    {
        /*@ assigns result[7]; */
        result[7] = x107 - x108*x109 - x63*x68 - x64*x73;

    }


    if (x31)
    {
        /*@ assigns result[8]; */
        result[8] = x101*x123;

    }
    else if (x43)
    {
        /*@ assigns result[8]; */
        result[8] = x78*(x24*(rt2*x105 - rt2*x88 + x124*x35 - x124*x38) + x79*x86 + x82*x85 - x84*x85 + x87*x91);

    }
    else if (x57)
    {
        /*@ assigns result[8]; */
        result[8] = mu*x38;

    }
    else if (x59)
    {
        /*@ assigns result[8]; */
        result[8] = -x120*x64 - x63*x71 - x93*(x15*x96 + x94) - x98*(x125*x96 - x94);

    }

    x65 = -mu*x20;
    x66 = mu*x23;
    x75 = mu - 1;
    /*@ assigns result[9]; */
    result[9] = x65 - x66 + x75 + 1;


    if (x70)
    {
        /*@ assigns result[10]; */
        result[10] = -x110*x72 - x20*x54;

    }
    else if (x46)
    {
        /*@ assigns result[10]; */
        result[10] = 0;
        /*@ assert result[10] >= 0; */
    }


    if (x70)
    {
        /*@ assigns result[11]; */
        result[11] = -x110*x116 - x20*x62;

    }
    else if (x46)
    {
        /*@ assigns result[11]; */
        result[11] = x65 + x66;

    }


    /*@ assert x46 || x70; */
    if (x46)
    {
        /*@ assigns result[12]; */
        result[12] = 0.0;
        /*@ assert result[12] >= 0; */
    }
    else if (x70)
    {
        /*@ assigns result[12]; */
        result[12] = x11*x69 - x68;

    }


    if (x70)
    {
        /*@ assigns result[13]; */
        result[13] = -x112*x12 + x113*x99 - x93*(-x114 + x95*x99) - x98*(x114 + x115*x95) + 1;

    }
    else if (x46)
    {
        /*@ assigns result[13]; */
        result[13] = 1;
        /*@ assert result[13] >= 0; */
        /*@ assert result[13] != 0; */
    }


    if (x46)
    {
        /*@ assigns result[14]; */
        result[14] = 0.0;
        /*@ assert result[14] >= 0; */
    }
    else if (x70)
    {
        /*@ assigns result[14]; */
        result[14] = -x108*x119 + x113*x118 + x117;

    }


    if (x46)
    {
        /*@ assigns result[15]; */
        result[15] = 0.0;
        /*@ assert result[15] >= 0; */
    }
    else if (x70)
    {
        /*@ assigns result[15]; */
        result[15] = x14*x69 - x71;

    }


    if (x70)
    {
        /*@ assigns result[16]; */
        result[16] = x108*x113 + x117 - x118*x119;

    }
    else if (x46)
    {
        /*@ assigns result[16]; */
        result[16] = 0;
        /*@ assert result[16] >= 0; */
    }


    if (x46)
    {
        /*@ assigns result[17]; */
        result[17] = -1.0*x38 + 1.0;

    }
    else if (x70)
    {
        /*@ assigns result[17]; */
        result[17] = -x112*x15 + x113*x125 - x93*(-x114 + x125*x95) - x98*(x114 + x126*x95) + 1;

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
    double result[18];
    fc3d_NaturalMapABGenerated(rn, rt1, rt2, un, ut1, ut2, mu, rhon, rhot1, rhot2, result);
    return(0);
}
#endif