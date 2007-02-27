#include <stdio.h>
#include <math.h>

// function to compute u
extern "C"   void SinPo(unsigned int sizeOfU, double time, double* UPtr, double* z)
{
  /* input parameter :
   *  sizeOfU (size of the vector u)
   *  sizeOfX (size of the vector x)
   *  time
   *  pointer to x
   *  pointer to xDot
   *  pointer to u (in-out parameter)
   */

  double omega = 1e4;
  double Voffset = 0.0;
  double amplitude = 10.0;
  double phase = 0.0;

  unsigned int size = sizeOfU;
  for (unsigned int i = 0; i < size; i++) UPtr[i] = Voffset + (amplitude * cos((omega * time) + phase));
}
