#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


projc(double xi[], int *nn, int statusi[], double pi[], double fric[], double *projc1, int *projc2)
{

  int i, nc, n = *nn, stat, k;
  double mina, maxa, bb, mu1, eps;



  //  n=sizeof(xi)/sizeof(xi[0]);
  nc = n / 2;
  bb = 0.0;
  eps = 1.e-10;

  //  printf("dans projc.c n %d \n",n);

  /*  for(k=0;k<nc;k++){
  printf("dans proj k %d statusi %d  xi %g\n",k,statusi[k],xi[k]);

  }*/


  for (i = 0; i < nc; i++)
  {
    mu1 = fric[i];
    stat = statusi[i];
    // printf("dans for \n");
    if (xi[2 * i] <= 0.0) // !status de non contact
    {
      projc1[2 * i] = 0.0;
      projc1[2 * i + 1] = 0.0;
      projc2[i] = 0;
    }
    //printf("if\n");}
    else
    {
      projc1[2 * i] = xi[2 * i]; //else ppl
      //       printf("else\n");
      if (xi[2 * i + 1] <= -mu1 * xi[2 * i]) //!slide backward
      {
        projc1[2 * i + 1] = -mu1 * xi[2 * i] ;
        projc2[i] = 1;
      }
      else if (xi[2 * i + 1] >= mu1 * xi[2 * i]) //!slide forward
      {
        projc1[2 * i + 1] = mu1 * xi[2 * i];
        projc2[i] = 3;
      }
      else
      {
        if (pi[2 * i + 1] == 0.0)
        {
          if (stat == 1) // !slide backward
          {
            projc1[2 * i + 1] = -mu1 * xi[2 * i];
            projc2[i] = 1;
          }
          else if (stat == 3) //!slide forward
          {
            projc1[2 * i + 1] = mu1 * xi[2 * i];
            projc2[i] = 3;
          }
          else
            //     !stick contact
          {
            projc1[2 * i + 1] = xi[2 * i + 1];
            projc2[i] = 2;
          }
        } // stat
        else
          //      !stick contact
        {
          projc1[2 * i + 1] = xi[2 * i + 1];
          projc2[i] = 2;
        }
      }//pi

    } // else ppl


  }
  /*
    for(k=0;k<nc;k++){
  printf("projc1_out %g projc2_out %d xi %g\n",projc1[k],projc2[k],xi[k]);

  }*/

  /*   for(k=0;k<nc;k++){
  printf("dans projc k %d statusi %d \n",k,statusi[k]);

  }*/


}

