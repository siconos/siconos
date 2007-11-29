#include "SiconosKernel.h"
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

#define PI 3.1415926535


/* Sum of beads by floor */
unsigned int sum(unsigned int et)
{
  return et * (et + 1) / 2;
}

/* Sum of beads of all floors */
unsigned int SUM(unsigned int et)
{
  unsigned int i, SUM;
  SUM = 0;
  if (et == 0) SUM = 0;
  else
  {
    for (i = 1 ; i <= et ; ++i)
      SUM += sum(i);
  }
  return SUM;
}


unsigned int qlq(unsigned int et, unsigned int j)
{
  unsigned int k;
  k = 0;
  if (j == 0) return 0;
  else
  {
    if (j == SUM(et - 1) - 1) return 0;
    else
    {
      for (k = 1 ; k < et ; k++)
      {
        if (SUM(et - 1) - 1 - j < sum(k + 1) && SUM(et - 1) - 1 - j >= sum(k)) return k;
      }
    }
  }
}

unsigned int L(unsigned int et, unsigned j)
{
  unsigned int k;
  if (j == 0) return et - 1;
  else
  {
    for (k = 1 ; k < et ; k++)
    {
      for (j = SUM(k) ; j < SUM(k + 1) ; ++j)
        return k;
    }
  }
}

/* Function to compute the 3 beads under each bead */
void Q(unsigned int et, unsigned int j, unsigned int i, double * x, double a, double b, double c, double d)
{
  double ti;
  ti = PI * (9 - 4 * i) / 6;
  x[0] = 0.12 * cos(ti) + a;
  x[1] = 0.12 * sin(ti) + b;
  x[2] = c - L(et, j) * d;
}
