// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "NSSTools.h"  // for diffns
#include "numerics_errors.h"

void diffns(int* na, int* a, int* nb, int* b, int* nc, int* c) {
  int pta, ptb, ptc;
  int aa, i;

  pta = 0;
  ptb = 0;
  ptc = 0;

  if (*nb == 0) {
    for (i = 0; i < *na; i++) c[i] = a[i];
    *nc = *na;

  }

  else {
    for (i = 0; i < *na; i++) c[i] = -1;

    while ((pta < *na) && (ptb < *nb)) {
      aa = a[pta];

      if (b[ptb] > aa) {
        c[ptc] = aa;
        ptc = ptc + 1;
        pta = pta + 1;
      } else if (b[ptb] == aa) {
        pta = pta + 1;

      } else {
        while ((b[ptb] < aa) && (ptb < *nb)) {
          ptb = ptb + 1;

          if (ptb >= *nb) {
            c[ptc] = aa;
            ptc = ptc + 1;

            break;
          }
        }
      }
    }

    for (i = pta + 1; i < *na; i++) {
      c[ptc] = a[i];
      ptc = ptc + 1;
    }

    *nc = ptc;
  }
}
