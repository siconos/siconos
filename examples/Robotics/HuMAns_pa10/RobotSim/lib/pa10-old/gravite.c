#include "robotmodel.h"



double gravity[N_DOF];

void modele_gravite(const double Q[N_DOF], double gravity[N_DOF])
{

  int i;

  /* Intermediate cos and sin values */
  double C1, C2, C3, C4, C5, C6, C7;
  double C11, C22, C33, C44, C55, C66, C77;
  double S1, S2, S3, S4, S5, S6, S7;

  C1 = cos(Q[0]);
  C2 = cos(Q[1]);
  C3 = cos(Q[2]);
  C4 = cos(Q[3]);
  C5 = cos(Q[4]);
  C6 = cos(Q[5]);
  C7 = cos(Q[6]);

  C11 = cos(2 * Q[0]);
  C22 = cos(2 * Q[1]);
  C33 = cos(2 * Q[2]);
  C44 = cos(2 * Q[3]);
  C55 = cos(2 * Q[4]);
  C66 = cos(2 * Q[5]);
  C77 = cos(2 * Q[6]);

  S1 = sin(Q[0]);
  S2 = sin(Q[1]);
  S3 = sin(Q[2]);
  S4 = sin(Q[3]);
  S5 = sin(Q[4]);
  S6 = sin(Q[5]);
  S7 = sin(Q[6]);

  gravity[0] = 0;

  gravity[1] = -(\
                 g * (S2 * (d3 * (m3 + m4 + m5 + m6 + m7) - m3 * ry3 + m2 * rz2 + \
                            C4 * (d5 * (m5 + m6 + m7) - m5 * ry5 + m4 * rz4 + \
                                  C6 * (m6 * rz6 + m7 * (d7 + rz7))) - \
                            C5 * (m6 * rz6 + m7 * (d7 + rz7)) * S4 * S6) + \
                      C2 * (C3 * (d5 * (m5 + m6 + m7) - m5 * ry5 + m4 * rz4 + \
                                  C6 * (m6 * rz6 + m7 * (d7 + rz7)))*\
                            S4 + (m6 * rz6 + m7 * (d7 + rz7)) * (C3 * C4 * C5 - S3 * S5) * S6)));

  gravity[2] = g*\
               S2 * ((d5 * (m5 + m6 + m7) - m5 * ry5 + m4 * rz4 + C6 * (m6 * rz6 + m7 * (d7 + rz7)))*\
                     S3 * S4 + (m6 * rz6 + m7 * (d7 + rz7)) * (C4 * C5 * S3 + C3 * S5) * S6);

  gravity[3] = -(\
                 g * (d5 * (m5 + m6 + m7) - m5 * ry5 + m4 * rz4 + \
                      C6 * (m6 * rz6 + m7 * (d7 + rz7))) * (C3 * C4 * S2 + C2 * S4)) + \
               C5 * g * (m6 * rz6 + m7 * (d7 + rz7)) * (-(C2 * C4) + C3 * S2 * S4) * S6;

  gravity[4] = g * (m6 * rz6 + m7 * (d7 + rz7)) * (C5 * S2 * S3 + (C3 * C4 * S2 + C2 * S4) * S5) * S6;

  gravity[5] = -(\
                 g * (m6 * rz6 + m7 * (d7 + rz7)) * (\
                     C6 * (C5 * (C3 * C4 * S2 + C2 * S4) - S2 * S3 * S5) + (C2 * C4 - C3 * S2 * S4) * S6));
  gravity[6] = 0;

}


