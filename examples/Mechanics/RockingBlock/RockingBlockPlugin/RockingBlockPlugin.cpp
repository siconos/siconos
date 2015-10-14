// Define plugin functions
//=============================================================================================================
#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <cmath>
#include <stdio.h>
using namespace std;
//===========================================================================================================
extern "C" double LengthBlock;
extern "C" double HeightBlock;
//1. Plugin function to calculate the gap function h1 at contact point 1 and h2 at contact point 2
SICONOS_EXPORT void h1(unsigned int sizeOfq, const double* q, unsigned int sizeOfy, double* y, unsigned int sizeOfZ, double* z)
{
  y[0] = q[1] - 0.5 * LengthBlock * sin(q[2]) - 0.5 * HeightBlock * cos(q[2]);
}
//
SICONOS_EXPORT void h2(unsigned int sizeOfq, const double* q, unsigned int sizeOfy, double* y, unsigned int sizeOfZ, double* z)
{
  y[0] = q[1] + 0.5 * LengthBlock * sin(q[2]) - 0.5 * HeightBlock * cos(q[2]);
}
//2. Plugin function to calculate the gradient of h1 (G1) and the one of h2 (G2)
SICONOS_EXPORT void G1(unsigned int sizeOfq, const double* q, unsigned int sizeOfy, double* Wn, unsigned int sizeOfZ, double* z)
{
  Wn[0] = 0.0;
  Wn[1] = 1.0;
  Wn[2] = -0.5 * LengthBlock * cos(q[2]) + 0.5 * HeightBlock * sin(q[2]);
}
//
SICONOS_EXPORT void G2(unsigned int sizeOfq, const double* q, unsigned int sizeOfy, double* Wn, unsigned int sizeOfZ, double* z)
{
  Wn[0] = 0.0;
  Wn[1] = 1.0;
  Wn[2] = 0.5 * LengthBlock * cos(q[2]) + 0.5 * HeightBlock * sin(q[2]);
}
// Plugin function to calculate the derivative of the Jacobian H with respect to the time (G1dot for contact 1 and G2dot for contact 2)
SICONOS_EXPORT void G1dot(unsigned int sizeOfq, const double* q, unsigned int sizeOfqdot, const double* qdot, double* S1, unsigned int sizeOfZ, double* z)
{
  S1[0] = 0.0;
  S1[1] = 0.0;
  S1[2] = (0.5 * LengthBlock * sin(q[2]) + 0.5 * HeightBlock * cos(q[2])) * qdot[2];
}
//
SICONOS_EXPORT void G2dot(unsigned int sizeOfq, const double* q, unsigned int sizeOfqdot, const double* qdot, double* S2, unsigned int sizeOfZ, double* z)
{
  S2[0] = 0.0;
  S2[1] = 0.0;
  S2[2] = (-0.5 * LengthBlock * sin(q[2]) + 0.5 * HeightBlock * cos(q[2])) * qdot[2];
}
