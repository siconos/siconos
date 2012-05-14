// Define plugin functions
//=============================================================================================================
//#include "SiconosKernel.hpp"
#include <math.h>
#include <stdio.h>
using namespace std;
//===========================================================================================================
extern double LengthBlock;
extern double HeightBlock;
//1. Plugin function to calculate the gap function h1 at contact point 1 and h2 at contact point 2
extern "C" void h1(unsigned int sizeOfq, const double* q, unsigned int sizeOfy, double* y, unsigned int sizeOfZ, double* z)
{
  y[0] = q[1] - 0.5 * LengthBlock * sin(q[2]) - 0.5 * HeightBlock * cos(q[2]);
}
//
extern "C" void h2(unsigned int sizeOfq, const double* q, unsigned int sizeOfy, double* y, unsigned int sizeOfZ, double* z)
{
  y[0] = q[1] + 0.5 * LengthBlock * sin(q[2]) - 0.5 * HeightBlock * cos(q[2]);
}
//2. Plugin function to calculate the gradient of h1 (G1) and the one of h2 (G2)
extern "C" void G1(unsigned int sizeOfq, const double* q, unsigned int sizeOfy, double* Wn, unsigned int sizeOfZ, double* z)
{
  Wn[0] = 0.0;
  Wn[1] = 1.0;
  Wn[2] = -0.5 * LengthBlock * cos(q[2]) + 0.5 * HeightBlock * sin(q[2]);
}
//
extern "C" void G2(unsigned int sizeOfq, const double* q, unsigned int sizeOfy, double* Wn, unsigned int sizeOfZ, double* z)
{
  Wn[0] = 0.0;
  Wn[1] = 1.0;
  Wn[2] = 0.5 * LengthBlock * cos(q[2]) + 0.5 * HeightBlock * sin(q[2]);
}