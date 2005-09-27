#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


max_part(double x[], double *sol, int *nn)

{

  int     i, n = *nn;
  double  max;


  max = x[0];

  for (i = 1 ; i < n ; i++)
  {

    if (max < x[i]) max = x[i] ;

  }


  *sol = max;

}
