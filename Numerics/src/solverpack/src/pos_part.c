#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void pos_part(double x[], double sol[], int *nn)

{

  int     i, n = *nn;





  for (i = 0 ; i < n ; i ++)
  {


    if (x[i] > 0.0)
    {
      sol[i] = x[i];
    }
    else
    {
      sol[i] =  0.0;
    }

  }
}

