#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



projf(int etat[], int *nn, double y[], double fric[], double projf1[])
{

  int i, nc, n = *nn;
  double mina, maxa, bb;



  //  n=sizeof(y)/sizeof(y[0]);
  nc = n / 2;
  bb = 0.0;

  for (i = 0; i < nc; i++)
  {
    if (etat[i] == 0) // no contact etat
    {
      if (y[2 * i] <= 0.0)
      {
        projf1[2 * i] = 0.0;
        projf1[2 * i + 1] = 0.0;
      }
      else
      {
        projf1[2 * i] = y[2 * i];
        projf1[2 * i + 1] = y[2 * i + 1];
      }
    }
    else if (etat[i] == 3) // !etat de contact glissant+
    {
      projf1[2 * i] = y[2 * i];
      minf(&y[2 * i + 1], &bb, &mina);
      projf1[2 * i + 1] = mina;
    }
    else if (etat[i] == 1) // !etat de contact glissant-
    {
      projf1[2 * i] = y[2 * i];
      maxf(&y[2 * i + 1], &bb, &maxa);
      projf1[2 * i + 1] = maxa;
    }
    else
      //     !etat de contact adhérent
    {
      projf1[2 * i] = y[2 * i];
      projf1[2 * i + 1] = y[2 * i + 1];
    }
  }


  /*  for(i=0;i<nc;i++)
      {  printf("i %d y %g projf1 %g etat %d\n",i,y[i],projf1[i],etat[i]);}*/

}




















