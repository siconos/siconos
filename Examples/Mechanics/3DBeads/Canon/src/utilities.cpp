#include "math.h"
#include "utilities.h"
unsigned int sum(unsigned int et)
{
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
  if (j == 0 || j == SUM(et - 1) - 1)
    return 0;
  else
  {
    for (unsigned int k = 1 ; k < et ; k++)
    {
      if (SUM(et - 1) - 1 - j < sum(k + 1) && SUM(et - 1) - 1 - j >= sum(k)) return k;
    }
  }
}

unsigned int L(unsigned int et, unsigned j)
{
  if (j == 0) return et - 1;
  else
  {
    for (unsigned int k = 1 ; k < et ; k++)
    {
      for (j = SUM(k) ; j < SUM(k + 1) ; ++j)
        return k;
    }
  }
}

void Q(unsigned int et, unsigned int j, unsigned int i, double * x, double a, double b, double c, double d)
{
  double ti;
  ti = PI * (9 - 4 * i) / 6;
  x[0] = 0.12 * cos(ti) + a;
  x[1] = 0.12 * sin(ti) + b;
  x[2] = c - L(et, j) * d;
}
