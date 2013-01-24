#include <stdio.h>
#include <math.h>


extern "C" double computeDefault(double time)
{
  double u;
  u = sin(50 * time);
  return u;
}

extern "C"   void computeB(double time, unsigned int sizeOfB, double* b, unsigned int sizeOfZ, double* z)
{
  b[0] = z[0];
  b[1] = z[1];
}
