#include <stdio.h>

extern "C" void vectorField(int *sizeOfX, double *time, double *x, double *xdot)
{
  /* input parameter : sizeOfX (size of the vector X); time ; x (pointer to X vector);
   * output parameter : xdot (pointer to Xdot vector)
   */

  printf("Call of the function 'vectorField' of the oscillator plugin.\n");

  if (sizeOfX[0] == 2)
  {
    xdot[0] = x[1];
    xdot[1] = -x[0];
  }
  else
  {
    printf("OscillatorPlugin:vectorField --- Bad size of x. %i -- %i\n", *sizeOfX, sizeOfX);
  }

}


extern "C"void computeJacobianX(int *sizeOfX, double *time, double *x, double *jacob)
{
  /* input parameter : sizeOfX (size of the vector X); time; x (pointer to x vector);
   * output parameter : jacob (pointer to JacobianX matrix)
   */

  printf("Call of the function 'ballJacobianX' of the oscillator plugin.\n");

  if (*sizeOfX == 2)
  {
    jacob[0] = 0;
    jacob[1] = -1;
    jacob[2] = 1;
    jacob[3] = 0;
  }
  else
  {
    printf("OscillatorPlugin:computeJacobianX --- Bad size of x. \n");
  }
}
