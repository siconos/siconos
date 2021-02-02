#include "fc3d_AlartCurnierJeanMoreauABGenerated.h"
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
  double x13 = 0.;
  int x14 = 0;

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
  x12 = mu*x11;
  x13 = Max(0.0000000000000002220446049250313080847263336181640625, x12);
  /*@ assert x13 >= 0; */
  /*@ assert x13 != 0; */
  x14 = x10 <= x13;
  /*@ assert x14 <==> (x10 <= x13); */

  int x21 = 0;

  double x15 = 0.;
  double x16 = 0.;
  double x17 = 0.;
  double x18 = 0.;
  double x19 = 0.;
  double x20 = 0.;
  double x22 = 0.;
  double x23 = 0.;
  double x24 = 0.;
  double x25 = 0.;
  double x26 = 0.;
  double x27 = 0.;
  double x28 = 0.;
  double x29 = 0.;
  double x30 = 0.;

  x21 = x10 > x13;
  /*@ assert x21 <==> (x10 > x13); */

  if(x21)
  {
    /*@ assert x10 < -1.09476442525e-47 || x10 > 1.09476442525e-47; */
    x15 = 1.0/x10;
    x16 = x13*x15;
    x17 = rt1 - x3;
    /*@ assert x9 < -1.09476442525e-47 || x9 > 1.09476442525e-47; */
    x18 = 1.0/((sqrt(x9))*(sqrt(x9))*(sqrt(x9)));
    x19 = x13*x17*x18;
    x20 = x19*x4;
    x22 = x19*x7;
    x23 = Heaviside(rn);
    x24 = Heaviside(x12 - 0.0000000000000002220446049250313080847263336181640625);
    x25 = mu*x15*x23*x24;
    x26 = -x16 + 1;
    x27 = rt2 - x6;
    x28 = x13*x18*x27;
    x29 = x28*x4;
    x30 = x28*x7;

  }
  double x2 = 0.;
  x2 = Heaviside(-rhon*un + rn);
  /*@ assigns result[0]; */
  result[0] = rhon*x2;

  /*@ assigns result[1]; */
  result[1] = 0;
  /*@ assert result[1] >= 0; */
  /*@ assigns result[2]; */
  result[2] = 0;
  /*@ assert result[2] >= 0; */
  /*@ assigns result[3]; */
  result[3] = 0;
  /*@ assert result[3] >= 0; */

  /*@ assert x14 || x21; */
  if(x14)
  {
    /*@ assigns result[4]; */
    result[4] = rhot1;

  }
  else if(x21)
  {
    /*@ assigns result[4]; */
    result[4] = rhot1*x16 + rhot1*x20;

  }


  if(x14)
  {
    /*@ assigns result[5]; */
    result[5] = 0;
    /*@ assert result[5] >= 0; */
  }
  else if(x21)
  {
    /*@ assigns result[5]; */
    result[5] = rhot1*x29;

  }

  /*@ assigns result[6]; */
  result[6] = 0;
  /*@ assert result[6] >= 0; */

  if(x14)
  {
    /*@ assigns result[7]; */
    result[7] = 0;
    /*@ assert result[7] >= 0; */
  }
  else if(x21)
  {
    /*@ assigns result[7]; */
    result[7] = rhot2*x22;

  }


  if(x14)
  {
    /*@ assigns result[8]; */
    result[8] = rhot2;

  }
  else if(x21)
  {
    /*@ assigns result[8]; */
    result[8] = rhot2*x16 + rhot2*x30;

  }

  /*@ assigns result[9]; */
  result[9] = -x2 + 1;


  if(x14)
  {
    /*@ assigns result[10]; */
    result[10] = 0;
    /*@ assert result[10] >= 0; */
  }
  else if(x21)
  {
    /*@ assigns result[10]; */
    result[10] = -x17*x25;

  }


  if(x14)
  {
    /*@ assigns result[11]; */
    result[11] = 0;
    /*@ assert result[11] >= 0; */
  }
  else if(x21)
  {
    /*@ assigns result[11]; */
    result[11] = -x25*x27;

  }

  /*@ assigns result[12]; */
  result[12] = 0;
  /*@ assert result[12] >= 0; */

  if(x14)
  {
    /*@ assigns result[13]; */
    result[13] = 0;
    /*@ assert result[13] >= 0; */
  }
  else if(x21)
  {
    /*@ assigns result[13]; */
    result[13] = -x20 + x26;

  }


  if(x14)
  {
    /*@ assigns result[14]; */
    result[14] = 0;
    /*@ assert result[14] >= 0; */
  }
  else if(x21)
  {
    /*@ assigns result[14]; */
    result[14] = -x29;

  }

  /*@ assigns result[15]; */
  result[15] = 0;
  /*@ assert result[15] >= 0; */

  if(x14)
  {
    /*@ assigns result[16]; */
    result[16] = 0;
    /*@ assert result[16] >= 0; */
  }
  else if(x21)
  {
    /*@ assigns result[16]; */
    result[16] = -x22;

  }


  if(x14)
  {
    /*@ assigns result[17]; */
    result[17] = 0;
    /*@ assert result[17] >= 0; */
  }
  else if(x21)
  {
    /*@ assigns result[17]; */
    result[17] = x26 - x30;

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
  fc3d_AlartCurnierJeanMoreauABGenerated(rn, rt1, rt2, un, ut1, ut2, mu, rhon, rhot1, rhot2, result);
  return(0);
}
#endif