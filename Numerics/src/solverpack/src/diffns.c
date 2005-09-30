#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
Input na, a, nb, b
Output nc, c
a and b: interger vectors in increasing order
c : vector of integers of a that are not in b.

 */

void diffns(int *na, int *a, int *nb, int * b, int *nc, int *c)
{

  int pta, ptb, ptc;
  int aa, i;





  pta = 0;
  ptb = 0;
  ptc = 0;


  if (*nb == 0)
  {

    for (i = 0 ; i < *na ; i++)
      c[i] = a[i];
    *nc  = *na;


  }
  else
  {

    for (i = 0 ; i < *na ; i++)
      c[i] = -1;


    //printf("bb %d \n", b[ptb]);




    while ((pta < *na) && (ptb < *nb))
    {

      aa  = a[pta];

      // printf("aa %d \n", aa );

      if (b[ptb] > aa)
      {

        c[ptc] = aa ;
        //printf("cc %d \n", c[ptc]);
        ptc    = ptc + 1 ;
        pta = pta + 1;
      }
      else if (b[ptb] == aa)
      {

        pta = pta + 1;

      }
      else
      {

        while ((b[ptb] < aa) && (ptb < *nb))
        {

          ptb = ptb + 1;

          if (ptb < *nb)
            //printf("bb %d \n", b[ptb]);

            if (ptb >= *nb)
            {

              c[ptc] = aa;
              //printf("cc %d \n", c[ptc]);
              ptc    = ptc + 1;

              break;

            }
        }

      }



    }



    for (i = pta + 1; i < *na ; i++)
    {


      c[ptc] = a[i];
      ptc = ptc + 1;
    }
    /**/





    *nc = ptc;

  }

}
