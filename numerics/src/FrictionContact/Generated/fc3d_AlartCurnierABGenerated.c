#include "fc3d_AlartCurnierABGenerated.h"
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
  double x2 = 0.;
  double x4 = 0.;
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
  int x15 = 0;

  x2 = -rhon*un + rn;
  x4 = rhot1*ut1;
  x5 = -rt1 + x4;
  x6 = x5*x5;
  /*@ assert x6 >= 0; */
  x7 = rhot2*ut2;
  x8 = -rt2 + x7;
  x9 = x8*x8;
  /*@ assert x9 >= 0; */
  x10 = x6 + x9;
  /*@ assert x10 >= 0; */
  x11 = sqrt(x10);
  x12 = Max(0, x2);
  /*@ assert x12 >= 0; */
  x13 = mu*x12;
  x14 = Max(0.0000000000000002220446049250313080847263336181640625, x13);
  /*@ assert x14 >= 0; */
  /*@ assert x14 != 0; */
  x15 = x11 <= x14;
  /*@ assert x15 <==> (x11 <= x14); */

  int x20 = 0;

  double x3 = 0.;
  double x16 = 0.;
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

  x20 = x11 > x14;
  /*@ assert x20 <==> (x11 > x14); */

  if(x20)
  {
    x3 = Heaviside(x2);
    x16 = rt1 - x4;
    x17 = Heaviside(x13 - 0.0000000000000002220446049250313080847263336181640625);
    /*@ assert x11 < -1.09476442525e-47 || x11 > 1.09476442525e-47; */
    x18 = 1.0/x11;
    x19 = mu*rhon*x17*x18*x3;
    x21 = x14*x18;
    /*@ assert x10 < -1.09476442525e-47 || x10 > 1.09476442525e-47; */
    x22 = 1.0/((sqrt(x10))*(sqrt(x10))*(sqrt(x10)));
    x23 = x14*x16*x22;
    x24 = x23*x5;
    x25 = x23*x8;
    x26 = mu*x17*x18*x3;
    x27 = -x21 + 1;
    x28 = rt2 - x7;
    x29 = x14*x22*x28;
    x30 = x29*x5;
    x31 = x29*x8;

  }
  x3 = Heaviside(x2);
  /*@ assigns result[0]; */
  result[0] = rhon*x3;


  /*@ assert x15 || x20; */
  if(x15)
  {
    /*@ assigns result[1]; */
    result[1] = 0;
    /*@ assert result[1] >= 0; */
  }
  else if(x20)
  {
    /*@ assigns result[1]; */
    result[1] = x16*x19;

  }


  if(x15)
  {
    /*@ assigns result[2]; */
    result[2] = 0;
    /*@ assert result[2] >= 0; */
  }
  else if(x20)
  {
    /*@ assigns result[2]; */
    result[2] = x19*x28;

  }

  /*@ assigns result[3]; */
  result[3] = 0;
  /*@ assert result[3] >= 0; */

  if(x15)
  {
    /*@ assigns result[4]; */
    result[4] = rhot1;

  }
  else if(x20)
  {
    /*@ assigns result[4]; */
    result[4] = rhot1*x21 + rhot1*x24;

  }


  if(x15)
  {
    /*@ assigns result[5]; */
    result[5] = 0;
    /*@ assert result[5] >= 0; */
  }
  else if(x20)
  {
    /*@ assigns result[5]; */
    result[5] = rhot1*x30;

  }

  /*@ assigns result[6]; */
  result[6] = 0;
  /*@ assert result[6] >= 0; */

  if(x15)
  {
    /*@ assigns result[7]; */
    result[7] = 0;
    /*@ assert result[7] >= 0; */
  }
  else if(x20)
  {
    /*@ assigns result[7]; */
    result[7] = rhot2*x25;

  }


  if(x15)
  {
    /*@ assigns result[8]; */
    result[8] = rhot2;

  }
  else if(x20)
  {
    /*@ assigns result[8]; */
    result[8] = rhot2*x21 + rhot2*x31;

  }

  /*@ assigns result[9]; */
  result[9] = -x3 + 1;


  if(x15)
  {
    /*@ assigns result[10]; */
    result[10] = 0;
    /*@ assert result[10] >= 0; */
  }
  else if(x20)
  {
    /*@ assigns result[10]; */
    result[10] = -x16*x26;

  }


  if(x15)
  {
    /*@ assigns result[11]; */
    result[11] = 0;
    /*@ assert result[11] >= 0; */
  }
  else if(x20)
  {
    /*@ assigns result[11]; */
    result[11] = -x26*x28;

  }

  /*@ assigns result[12]; */
  result[12] = 0;
  /*@ assert result[12] >= 0; */

  if(x15)
  {
    /*@ assigns result[13]; */
    result[13] = 0;
    /*@ assert result[13] >= 0; */
  }
  else if(x20)
  {
    /*@ assigns result[13]; */
    result[13] = -x24 + x27;

  }


  if(x15)
  {
    /*@ assigns result[14]; */
    result[14] = 0;
    /*@ assert result[14] >= 0; */
  }
  else if(x20)
  {
    /*@ assigns result[14]; */
    result[14] = -x30;

  }

  /*@ assigns result[15]; */
  result[15] = 0;
  /*@ assert result[15] >= 0; */

  if(x15)
  {
    /*@ assigns result[16]; */
    result[16] = 0;
    /*@ assert result[16] >= 0; */
  }
  else if(x20)
  {
    /*@ assigns result[16]; */
    result[16] = -x25;

  }


  if(x15)
  {
    /*@ assigns result[17]; */
    result[17] = 0;
    /*@ assert result[17] >= 0; */
  }
  else if(x20)
  {
    /*@ assigns result[17]; */
    result[17] = x27 - x31;

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
  fc3d_AlartCurnierABGenerated(rn, rt1, rt2, un, ut1, ut2, mu, rhon, rhot1, rhot2, result);
  return(0);
}
#endif
