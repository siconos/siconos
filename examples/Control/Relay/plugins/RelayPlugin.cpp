#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <stdio.h>
#include <math.h>


extern "C" double computeDefault(double time)
{
  double u;
  u = sin(50 * time);
  return u;
}


SICONOS_EXPORT void computeE(double time, unsigned int sizeOfB, double* e, unsigned int sizeOfZ, double* z)
{
  e[0] = 0.0;//computeDefault(time);
}

SICONOS_EXPORT void computeB(double time, unsigned int sizeOfB, double* b, unsigned int sizeOfZ, double* z)
{
  b[0] = computeDefault(time) + z[0];
  b[1] = -1.0 * computeDefault(time) + z[1];
}
