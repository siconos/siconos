#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void abs_part(double x[], double sol[], int *nn)

{

  int      i, n = *nn;
  double   eps;


  eps = 1.e-12;


  for (i = 0 ; i < n ; i++)
  {


    if (x[i] >= eps)
    {
      sol[i] = x[i];
    }
    else
    {
      sol[i] =  -x[i];
    }

  }

}
