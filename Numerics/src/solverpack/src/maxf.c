#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/** \fn void maxf(double *a,double *b,double *c)
 *  \param double* :
 *  \param double* :
 *  \param double* :
 */
void maxf(double *a, double *b, double *c)
{
  double aa = *a, bb = *b;

  if (aa < bb)
  {
    *c = bb;
  }
  else
  {
    *c = aa;
  }
}
