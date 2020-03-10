#include "fc3d_NaturalMapFABGenerated.h"
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
  int x34 = 0;
  int x35 = 0;
  int x36 = 0;

  double x4 = 0.;
  double x6 = 0.;
  double x29 = 0.;
  double x30 = 0.;
  double x31 = 0.;
  double x32 = 0.;
  double x33 = 0.;
  double x81 = 0.;
  double x82 = 0.;
  double x83 = 0.;
  double x103 = 0.;
  double x104 = 0.;
  double x125 = 0.;
  double x126 = 0.;
  double x127 = 0.;

  x1 = rt1*rt1 + rt2*rt2;
  /*@ assert x1 >= 0; */
  x2 = ut1*ut1 + ut2*ut2;
  /*@ assert x2 >= 0; */
  x34 = x1 <= epsilon;
  /*@ assert x34 <==> (x1 <= epsilon); */
  x35 = x2 <= epsilon;
  /*@ assert x35 <==> (x2 <= epsilon); */
  x36 = x34 && x35;
  /*@ assert x36 <==> (x34 && x35); */

  int x47 = 0;
  int x48 = 0;

  double x7 = 0.;
  double x37 = 0.;
  double x38 = 0.;
  double x39 = 0.;
  double x40 = 0.;
  double x41 = 0.;
  double x42 = 0.;
  double x43 = 0.;
  double x44 = 0.;
  double x45 = 0.;
  double x46 = 0.;
  double x65 = 0.;
  double x84 = 0.;
  double x85 = 0.;
  double x86 = 0.;
  double x87 = 0.;
  double x88 = 0.;
  double x89 = 0.;
  double x90 = 0.;
  double x91 = 0.;
  double x92 = 0.;
  double x93 = 0.;
  double x94 = 0.;
  double x95 = 0.;
  double x96 = 0.;
  double x97 = 0.;
  double x98 = 0.;
  double x105 = 0.;
  double x106 = 0.;
  double x107 = 0.;
  double x108 = 0.;
  double x109 = 0.;
  double x128 = 0.;

  x47 = x1 > epsilon;
  /*@ assert x47 <==> (x1 > epsilon); */
  x48 = x35 && x47;
  /*@ assert x48 <==> (x35 && x47); */

  double x10 = 0.;
  double x11 = 0.;
  double x12 = 0.;
  double x13 = 0.;
  double x14 = 0.;
  double x15 = 0.;
  double x16 = 0.;
  int x49 = 0;
  double x50 = 0.;
  int x51 = 0;
  int x52 = 0;

  x10 = mu*ut1;
  x11 = rt1 - x10;
  x12 = x11*x11;
  /*@ assert x12 >= 0; */
  x13 = mu*ut2;
  x14 = rt2 - x13;
  x15 = x14*x14;
  /*@ assert x15 >= 0; */
  x16 = x12 + x15;
  x49 = x2 >= epsilon;
  /*@ assert x49 <==> (x2 >= epsilon); */
  x50 = epsilon*(mu + 1);
  /*@ assert x50 >= 0; */
  /*@ assert x50 != 0; */
  x51 = x16 <= x50;
  /*@ assert x51 <==> (x16 <= x50); */
  x52 = x34 && x49 && x51;
  /*@ assert x52 <==> (x34 && x49 && x51); */

  int x54 = 0;

  double x53 = 0.;

  x54 = x47 && x49 && x51;
  /*@ assert x54 <==> (x47 && x49 && x51); */

  int x62 = 0;
  int x63 = 0;
  int x64 = 0;

  double x8 = 0.;
  double x9 = 0.;
  double x17 = 0.;
  double x18 = 0.;
  double x19 = 0.;
  double x20 = 0.;
  double x21 = 0.;
  double x22 = 0.;
  double x23 = 0.;
  double x24 = 0.;
  double x25 = 0.;
  double x26 = 0.;
  double x27 = 0.;
  double x28 = 0.;
  double x55 = 0.;
  double x56 = 0.;
  double x57 = 0.;
  double x58 = 0.;
  double x59 = 0.;
  double x60 = 0.;
  double x61 = 0.;
  double x66 = 0.;
  double x67 = 0.;
  double x68 = 0.;
  double x69 = 0.;
  double x72 = 0.;
  double x73 = 0.;
  double x74 = 0.;
  double x76 = 0.;
  double x78 = 0.;
  double x80 = 0.;
  double x99 = 0.;
  double x100 = 0.;
  double x101 = 0.;
  double x102 = 0.;
  double x110 = 0.;
  double x111 = 0.;
  double x112 = 0.;
  double x113 = 0.;
  double x120 = 0.;
  double x122 = 0.;
  double x124 = 0.;
  double x129 = 0.;

  x62 = x49 && x51;
  /*@ assert x62 <==> (x49 && x51); */
  x63 = x35 || x62;
  /*@ assert x63 <==> (x35 || x62); */
  x64 = !x63;
  /*@ assert x64 <==> (!x63); */

  double x70 = 0.;
  double x71 = 0.;

  int x75 = 0;

  double x77 = 0.;
  double x79 = 0.;
  double x114 = 0.;
  double x115 = 0.;
  double x116 = 0.;
  double x117 = 0.;
  double x118 = 0.;
  double x119 = 0.;
  double x121 = 0.;
  double x123 = 0.;
  double x130 = 0.;

  x75 = x16 > x50;
  /*@ assert x75 <==> (x16 > x50); */

  if(x36)
  {
    x4 = 2*mu*mu - 2*mu + 1;
    x6 = mu*rn;
    /*@ assert x6 >= 0; */
    /*@ assert 2 >= 0; */
    x29 = sqrt(2);
    /*@ assert x29 >= 0; */
    /*@ assert x29 != 0; */
    x30 = -1.0*un + x6;
    x31 = Heaviside(x30);
    x32 = mu*x31;
    x33 = (1.0/2.0)*x29*x32;
    /*@ assert x4 < -epsilon || x4 > epsilon; */
    x81 = 1.0/x4;
    x82 = mu - 1;
    x83 = x82*x82;
    /*@ assert x83 >= 0; */
    x103 = mu*mu;
    /*@ assert x103 >= 0; */
    x104 = x31*x81;
    x125 = -2.0*mu + 2*x103 + 1.0;
    /*@ assert x125 < -epsilon || x125 > epsilon; */
    x126 = 1.0/x125;
    x127 = mu*mu*mu;
    /*@ assert x127 >= 0; */

  }
  if(x48)
  {
    x6 = mu*rn;
    x7 = -un;
    /*@ assert 2 >= 0; */
    x29 = sqrt(2);
    x30 = -1.0*un + x6;
    x37 = sqrt(x1);
    /*@ assert x37 < -epsilon || x37 > epsilon; */
    x38 = 1.0/x37;
    x39 = 0.25*mu*x38;
    x40 = Heaviside(x30 + x37);
    x41 = 2*x40;
    x42 = rt1*x41;
    x43 = Heaviside(x30 - 1.0*x37);
    x44 = 2.0*x43;
    x45 = x29*x37;
    x46 = x40*x45 + x43*x45;
    x65 = rt2*x41;
    /*@ assert x1 < -epsilon || x1 > epsilon; */
    x84 = 1.0/((sqrt(x1))*(sqrt(x1))*(sqrt(x1)));
    x85 = 0.25*mu*x84;
    x86 = rt2*rt2;
    /*@ assert x86 >= 0; */
    x87 = x6 + x7;
    x88 = Max(0, x37 + x87);
    /*@ assert x88 >= 0; */
    x89 = 2*x88;
    x90 = Max(0, -x37 + x87);
    /*@ assert x90 >= 0; */
    x91 = 2.0*x90;
    x92 = rt1*rt1;
    /*@ assert x92 >= 0; */
    x93 = 2*x37*x40;
    x94 = 2*x37;
    x95 = x43*x92;
    x96 = rt1*rt1*rt1;
    x97 = x40*x86;
    x98 = x43*x86;
    x105 = 2.0*x88;
    x106 = 2*x90;
    x107 = 2*x37*x43;
    x108 = x40*x92;
    x109 = x29*(x108 - x95 + x97 - x98);
    x128 = rt2*rt2*rt2;

  }
  if(x52)
  {
    x6 = mu*rn;
    /*@ assert 2 >= 0; */
    x29 = sqrt(2);
    x30 = -1.0*un + x6;
    x31 = Heaviside(x30);
    x32 = mu*x31;
    x33 = (1.0/2.0)*x29*x32;

  }
  if(x54)
  {
    x6 = mu*rn;
    x30 = -1.0*un + x6;
    x37 = sqrt(x1);
    /*@ assert x37 < -epsilon || x37 > epsilon; */
    x38 = 1.0/x37;
    x43 = Heaviside(x30 - 1.0*x37);
    x53 = mu*x38*x43;

  }
  if(x64)
  {
    x6 = mu*rn;
    x7 = -un;
    x8 = sqrt(x2);
    x9 = -mu*x8 + x6 + x7;
    /*@ assert x16 >= 0; */
    x17 = sqrt(x16);
    x18 = x17 + x9;
    x19 = Max(0, x18);
    /*@ assert x19 >= 0; */
    x20 = 0.5*x19;
    x21 = -x20;
    x22 = -x17 + x9;
    x23 = Max(0, x22);
    /*@ assert x23 >= 0; */
    x24 = 0.5*x23;
    x25 = Heaviside(x18);
    x26 = 0.5*x25;
    x27 = Heaviside(x22);
    x28 = 0.5*x27;
    /*@ assert x8 < -epsilon || x8 > epsilon; */
    x55 = 1.0/x8;
    x56 = -x10*x55;
    /*@ assert x17 < -epsilon || x17 > epsilon; */
    x57 = 1.0/x17;
    x58 = mu*x57;
    x59 = x11*x58;
    x60 = x56 - x59;
    x61 = x56 + x59;
    x66 = -x13*x55;
    x67 = x14*x58;
    x68 = x66 - x67;
    x69 = x66 + x67;
    x72 = 0.5*x25*x57;
    x73 = x11*x72;
    x74 = 0.5*x27*x57;
    x76 = x14*x72;
    x78 = -rt1 + x10;
    x80 = x74*x78;
    x99 = -x58;
    /*@ assert x16 < -epsilon || x16 > epsilon; */
    x100 = 1.0/((sqrt(x16))*(sqrt(x16))*(sqrt(x16)));
    x101 = mu*x100;
    x102 = x11*x78;
    x110 = x11*x14;
    x111 = -0.5*mu*x100*x110*x19;
    x112 = x14*x78;
    x113 = 0.5*mu*x100*x23;
    x120 = -rt2 + x13;
    x122 = x11*x120;
    x124 = x120*x74;
    x129 = x120*x14;

  }
  if(x51)
  {
    x6 = mu*rn;
    x7 = -un;
    x8 = sqrt(x2);
    x9 = -mu*x8 + x6 + x7;
    /*@ assert x16 >= 0; */
    x17 = sqrt(x16);
    x18 = x17 + x9;
    x19 = Max(0, x18);
    x20 = 0.5*x19;
    x21 = -x20;
    x22 = -x17 + x9;
    x23 = Max(0, x22);
    x24 = 0.5*x23;
    x25 = Heaviside(x18);
    x26 = 0.5*x25;
    x27 = Heaviside(x22);
    x28 = 0.5*x27;
    x30 = -1.0*un + x6;
    x37 = sqrt(x1);
    x43 = Heaviside(x30 - 1.0*x37);
    x70 = -mu*x26;
    x71 = mu*x28;

  }
  if(x75)
  {
    x6 = mu*rn;
    x7 = -un;
    x8 = sqrt(x2);
    x9 = -mu*x8 + x6 + x7;
    /*@ assert x16 >= 0; */
    x17 = sqrt(x16);
    x18 = x17 + x9;
    x19 = Max(0, x18);
    x20 = 0.5*x19;
    x21 = -x20;
    x22 = -x17 + x9;
    x23 = Max(0, x22);
    x24 = 0.5*x23;
    x25 = Heaviside(x18);
    x26 = 0.5*x25;
    x27 = Heaviside(x22);
    /*@ assert x17 < -epsilon || x17 > epsilon; */
    x57 = 1.0/x17;
    x58 = mu*x57;
    x59 = x11*x58;
    x67 = x14*x58;
    x72 = 0.5*x25*x57;
    x73 = x11*x72;
    x74 = 0.5*x27*x57;
    x76 = x14*x72;
    x77 = 0.5*x19*x57;
    x78 = -rt1 + x10;
    x79 = 0.5*x23*x57;
    x80 = x74*x78;
    /*@ assert x16 < -epsilon || x16 > epsilon; */
    x100 = 1.0/((sqrt(x16))*(sqrt(x16))*(sqrt(x16)));
    x102 = x11*x78;
    x110 = x11*x14;
    x112 = x14*x78;
    x114 = 0.5*mu*x27*x57;
    x115 = 1.0/x16;
    x116 = 0.5*x115*x25;
    x117 = 0.5*x115*x27;
    x118 = -x57;
    x119 = x78*x78;
    /*@ assert x119 >= 0; */
    x120 = -rt2 + x13;
    x121 = -x100*x120*x24*x78 - x110*x116;
    x122 = x11*x120;
    x123 = 0.5*x100*x19;
    x124 = x120*x74;
    x129 = x120*x14;
    x130 = x120*x120;
    /*@ assert x130 >= 0; */

  }
  if(x62)
  {
    x6 = mu*rn;
    x30 = -1.0*un + x6;
    x37 = sqrt(x1);
    x43 = Heaviside(x30 - 1.0*x37);

  }
  x6 = mu*rn;
  x7 = -un;
  x8 = sqrt(x2);
  x9 = -mu*x8 + x6 + x7;
  /*@ assert x16 >= 0; */
  x17 = sqrt(x16);
  x18 = x17 + x9;
  x19 = Max(0, x18);
  x20 = 0.5*x19;
  x21 = -x20;
  x22 = -x17 + x9;
  x23 = Max(0, x22);
  x24 = 0.5*x23;
  /*@ assigns result[0]; */
  result[0] = x21 - x24 + x6;


  /*@ assert x75 || x51; */
  if(x75)
  {
    /*@ assigns result[1]; */
    result[1] = rt1 - x11*x77 - x78*x79;

  }
  else if(x51)
  {
    /*@ assigns result[1]; */
    result[1] = rt1;

  }


  if(x75)
  {
    /*@ assigns result[2]; */
    result[2] = rt2 - x120*x79 - x14*x77;

  }
  else if(x51)
  {
    /*@ assigns result[2]; */
    result[2] = rt2 + x21 + x24;

  }

  x25 = Heaviside(x18);
  x26 = 0.5*x25;
  x27 = Heaviside(x22);
  x28 = 0.5*x27;
  /*@ assigns result[3]; */
  result[3] = x26 + x28;


  if(x75)
  {
    /*@ assigns result[4]; */
    result[4] = x73 + x80;

  }
  else if(x51)
  {
    /*@ assigns result[4]; */
    result[4] = 0;
    /*@ assert result[4] >= 0; */
  }


  if(x75)
  {
    /*@ assigns result[5]; */
    result[5] = x124 + x76;

  }
  else if(x51)
  {
    /*@ assigns result[5]; */
    result[5] = x26 - x28;

  }


  /*@ assert x36 || x48 || x52 || x54 || x64; */
  if(x36)
  {
    /*@ assigns result[6]; */
    result[6] = x33;

  }
  else if(x48)
  {
    /*@ assigns result[6]; */
    result[6] = x39*(-rt1*x44 + x42 + x46);

  }
  else if(x52)
  {
    /*@ assigns result[6]; */
    result[6] = x33;

  }
  else if(x54)
  {
    /*@ assigns result[6]; */
    result[6] = rt1*x53;

  }
  else if(x64)
  {
    /*@ assigns result[6]; */
    result[6] = -x26*x60 - x28*x61;

  }


  /*@ assert x36 || x48 || x62 || x64; */
  if(x36)
  {
    /*@ assigns result[7]; */
    result[7] = 1.0*x32*x81*x83;

  }
  else if(x48)
  {
    /*@ assigns result[7]; */
    result[7] = x85*(x29*(rt1*x97 - rt1*x98 + x40*x96 - x43*x96) + x86*x89 - x86*x91 + x92*x93 + x94*x95);

  }
  else if(x62)
  {
    /*@ assigns result[7]; */
    result[7] = 0.0;
    /*@ assert result[7] >= 0; */
  }
  else if(x64)
  {
    /*@ assigns result[7]; */
    result[7] = x21*(x101*x12 + x99) - x24*(x101*x102 - x99) - x60*x73 - x61*x80;

  }


  if(x36)
  {
    /*@ assigns result[8]; */
    result[8] = 1.0*x103*x126*x31*x82;

  }
  else if(x48)
  {
    /*@ assigns result[8]; */
    result[8] = rt2*x85*(-rt1*x105 + rt1*x106 + rt1*x107 + x109 + x37*x42);

  }
  else if(x62)
  {
    /*@ assigns result[8]; */
    result[8] = 0.0;
    /*@ assert result[8] >= 0; */
  }
  else if(x64)
  {
    /*@ assigns result[8]; */
    result[8] = x111 - x113*x122 - x124*x61 - x60*x76;

  }


  if(x36)
  {
    /*@ assigns result[9]; */
    result[9] = x33;

  }
  else if(x48)
  {
    /*@ assigns result[9]; */
    result[9] = x39*(-rt2*x44 + x46 + x65);

  }
  else if(x52)
  {
    /*@ assigns result[9]; */
    result[9] = x33;

  }
  else if(x54)
  {
    /*@ assigns result[9]; */
    result[9] = rt2*x53;

  }
  else if(x64)
  {
    /*@ assigns result[9]; */
    result[9] = -x26*x68 - x28*x69;

  }


  if(x36)
  {
    /*@ assigns result[10]; */
    result[10] = x103*x104*(mu - 1.0);

  }
  else if(x48)
  {
    /*@ assigns result[10]; */
    result[10] = rt1*x85*(-rt2*x105 + rt2*x106 + rt2*x107 + x109 + x37*x65);

  }
  else if(x62)
  {
    /*@ assigns result[10]; */
    result[10] = 0.0;
    /*@ assert result[10] >= 0; */
  }
  else if(x64)
  {
    /*@ assigns result[10]; */
    result[10] = x111 - x112*x113 - x68*x73 - x69*x80;

  }


  if(x36)
  {
    /*@ assigns result[11]; */
    result[11] = x104*x127;

  }
  else if(x48)
  {
    /*@ assigns result[11]; */
    result[11] = x85*(x29*(rt2*x108 - rt2*x95 + x128*x40 - x128*x43) + x86*x93 + x89*x92 - x91*x92 + x94*x98);

  }
  else if(x62)
  {
    /*@ assigns result[11]; */
    result[11] = mu*x43;

  }
  else if(x64)
  {
    /*@ assigns result[11]; */
    result[11] = -x124*x69 + x21*(x101*x15 + x99) - x24*(x101*x129 - x99) - x68*x76;

  }

  x70 = -mu*x26;
  x71 = mu*x28;
  x82 = mu - 1;
  /*@ assigns result[12]; */
  result[12] = x70 - x71 + x82 + 1;


  if(x75)
  {
    /*@ assigns result[13]; */
    result[13] = -x114*x78 - x26*x59;

  }
  else if(x51)
  {
    /*@ assigns result[13]; */
    result[13] = 0;
    /*@ assert result[13] >= 0; */
  }


  if(x75)
  {
    /*@ assigns result[14]; */
    result[14] = -x114*x120 - x26*x67;

  }
  else if(x51)
  {
    /*@ assigns result[14]; */
    result[14] = x70 + x71;

  }


  /*@ assert x51 || x75; */
  if(x51)
  {
    /*@ assigns result[15]; */
    result[15] = 0.0;
    /*@ assert result[15] >= 0; */
  }
  else if(x75)
  {
    /*@ assigns result[15]; */
    result[15] = x11*x74 - x73;

  }


  if(x75)
  {
    /*@ assigns result[16]; */
    result[16] = x102*x117 - x116*x12 + x21*(x100*x102 + x57) - x24*(x100*x119 + x118) + 1;

  }
  else if(x51)
  {
    /*@ assigns result[16]; */
    result[16] = 1;
    /*@ assert result[16] >= 0; */
    /*@ assert result[16] != 0; */
  }


  if(x51)
  {
    /*@ assigns result[17]; */
    result[17] = 0.0;
    /*@ assert result[17] >= 0; */
  }
  else if(x75)
  {
    /*@ assigns result[17]; */
    result[17] = -x112*x123 + x117*x122 + x121;

  }


  if(x51)
  {
    /*@ assigns result[18]; */
    result[18] = 0.0;
    /*@ assert result[18] >= 0; */
  }
  else if(x75)
  {
    /*@ assigns result[18]; */
    result[18] = x14*x74 - x76;

  }


  if(x75)
  {
    /*@ assigns result[19]; */
    result[19] = x112*x117 + x121 - x122*x123;

  }
  else if(x51)
  {
    /*@ assigns result[19]; */
    result[19] = 0;
    /*@ assert result[19] >= 0; */
  }


  if(x51)
  {
    /*@ assigns result[20]; */
    result[20] = -1.0*x43 + 1.0;

  }
  else if(x75)
  {
    /*@ assigns result[20]; */
    result[20] = -x116*x15 + x117*x129 + x21*(x100*x129 + x57) - x24*(x100*x130 + x118) + 1;

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
  fc3d_NaturalMapFABGenerated(rn, rt1, rt2, un, ut1, ut2, mu, rhon, rhot1, rhot2, result);
  return(0);
}
#endif