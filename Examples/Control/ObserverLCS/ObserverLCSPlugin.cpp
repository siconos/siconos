#include <stdio.h>
#include <math.h>


extern "C" double computeControl(double time)
{
  double u;
  double alpha = 1.0;
  double T = 0.1;
  int oddoreven = int(time / T);
  printf("oodoreven = %i \n ", int(time));
  printf("oodoreven = %i \n ", int(T));
  printf("oodoreven = %i \n ", oddoreven);

  if ((oddoreven / 2) == 0) u = alpha;
  else u = 0;
  u = 30 * sin(50 * time);
  return u;
}

extern "C" void uProcess(double time, unsigned int sizeOfB, double* b, unsigned int sizeOfZ, double* z)
{
  double u = computeControl(time);
  b[0] = u ;
  b[1] = 2.0 * u;
  printf("End u(t) \n ");
}

extern "C"   void uObserver(double time, unsigned int sizeOfB, double* b, unsigned int sizeOfZ, double* z)
{
  printf("Compute u for observer \n ");
  double u = computeControl(time);
  double L[2];
  L[0] = 1.0;
  L[1] = 1.0;

  b[0] = u + L[0] * z[0];
  b[1] = 2.0 * u + L[1] * z[0];
  printf("End u observer \n ");
}

extern "C"   void computeE(double time, unsigned int sizeOfB, double* e, unsigned int sizeOfZ, double* z)
{
  printf("ComputeE\n ");
  e[0] = computeControl(time);
  e[1] = 0.0;
}
