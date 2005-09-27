#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


min_part(double x[], double *sol, int *nn)

{

  int     i, n = *nn;
  double  minx;



  minx = x[0];

  for (i = 1 ; i < n ; i++)
  {

    if (minx > x[i]) minx = x[i] ;

  }


  *sol = minx;

}
