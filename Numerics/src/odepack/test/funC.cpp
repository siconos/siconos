#include "g2c.h"
#include <stdio.h>

#define F77NAME(x) x##_


extern "C" void F77NAME(f1)(integer *sizeOfX, doublereal *time, doublereal *x, doublereal *xdot)
{
  /* input parameter : sizeOfX (size of the vector X); time ; x (pointer to X vector);
   * output parameter : xdot (pointer to Xdot vector)
   */

  //        printf("Call of the function 'f1' of the file funC.\n");

  if (sizeOfX[0] == 2)
  {
    xdot[0] = x[1];
    xdot[1] = 3.0 * (1.0 - x[0] * x[0]) * x[1] - x[0];
  }
  else
  {
    printf("OscillatorPlugin:vectorField --- Bad size of x. %i -- %i\n", *sizeOfX, sizeOfX);
  }

}

extern "C"void F77NAME(jac1)(integer *sizeOfX, doublereal *time, doublereal *x, integer* ml, integer *mu,  doublereal *jacob, integer *nrowpd)
{
  /* input parameter : sizeOfX (size of the vector X); time; x (pointer to x vector);
   * output parameter : jacob (pointer to JacobianX matrix)
   */

  //      printf("Call of the function 'jac1' of the  the file funC.c .\n");

  if (*sizeOfX == 2)
  {
    jacob[0] = 0.0;
    jacob[1] = -6.0 * x[0] * x[1] - 1.0;
    jacob[2] = 1.0;
    jacob[3] =  3.0 * (1.0 - x[0] * x[0]);
  }
  else
  {
    printf("OscillatorPlugin:computeJacobianX --- Bad size of x. \n");
  }
}
