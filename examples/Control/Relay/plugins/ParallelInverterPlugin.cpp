#ifdef _WIN32
#define SICONOS_EXPORT extern "C" __declspec(dllexport)
#else
#define SICONOS_EXPORT extern "C"
#endif

// for M_PI
#define _USE_MATH_DEFINES
#include <math.h>

SICONOS_EXPORT void eLDS(double t, unsigned int N, double* e, unsigned int z, double*zz)
{
  double Z[3][2];

  Z[0][0] = -1;
  Z[0][1] = -6e-5;
  Z[1][0] = -1;
  Z[1][1] = -6e-5;
  Z[2][0] = -1;
  Z[2][1] = -6e-5;
  zz[0] = 55 * sin(100 * M_PI * t);

  e[0] = Z[0][0]*zz[0];
  e[1] = Z[1][0]*zz[0];
  e[2] = Z[2][0]*zz[0];
}


