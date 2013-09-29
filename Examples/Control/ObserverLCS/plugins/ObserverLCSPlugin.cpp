#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <stdio.h>
#include <math.h>


extern "C" double computeControl(double time)
{
  double u;
  double alpha = 1.0;
  double T = 0.1;
  int oddoreven = int(time / T);

  if ((oddoreven / 2) == 0) u = alpha;
  else u = 0;
  u = 30 * sin(50 * time);
  return u;
}

SICONOS_EXPORT void uProcess(double time, unsigned int sizeOfB, double* b, unsigned int sizeOfZ, double* z)
{
  double u = computeControl(time);
  b[0] = u ;
  b[1] = 2.0 * u;
}

SICONOS_EXPORT void uObserver(double time, unsigned int sizeOfB, double* b, unsigned int sizeOfZ, double* z)
{
  double u = computeControl(time);
//  double L[2];
//  L[0] = 1.0;
//  L[1] = 1.0;

  //   b[0] = u + L[0]*z[0];
  //   b[1] = 2.0* u + L[1]*z[0];
  double coeff = 2 * (z[0] + z[1]);
  b[0] = u + coeff;
  b[1] = 2.0 * u + coeff;
}

SICONOS_EXPORT void computeE(double time, unsigned int sizeOfB, double* e, unsigned int sizeOfZ, double* z)
{
  e[0] = computeControl(time);
}
