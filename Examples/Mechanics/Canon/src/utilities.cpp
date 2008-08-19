#include "math.h"
#include "utilities.h"
unsigned int sum(unsigned int et)
{
  if (et == 0) return 0;
  else
    return et * (et + 1) / 2;
}

unsigned int SUM(unsigned int et)
{
  if (et == 0) return 0;
  else
  {
    unsigned int s = 0;
    for (int i = 1 ; i <= et ; ++i)
      s += sum(i);
    return s;
  }
}

unsigned int qlq(unsigned int et, unsigned int j)
{
  int val = SUM(et - 1) - 1;
  if (j == 0 || j == val)
    return 0;
  else
  {
    for (unsigned int k = 1 ; k < et ; k++)
    {
      if (val - j < sum(k + 1) && val - j >= sum(k)) return k;
    }
  }
}

void Q(unsigned int et, unsigned int j, unsigned int i, double * x, double qx, double qy, double qz, double RR)
{
  double coef = 2.0 / sqrt(3.0) * RR;
  double ti = PI * (9 - 4 * i) / 6.0;
  x[0] = coef * cos(ti) + qx;
  x[1] = coef * sin(ti) + qy;
  x[2] = qz - 2 * sqrt(2.0 / 3.0) * RR;
}
