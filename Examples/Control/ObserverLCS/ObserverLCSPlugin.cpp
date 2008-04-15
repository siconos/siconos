#include <stdio.h>
#include <math.h>


extern "C"   void computeU(double time, unsigned int sizeOfB, double* b, unsigned int sizeOfZ, double* z)
{

  // Si b est de taille 2

  printf("start ComputeU\n ");
  double alpha = 1.0;
  double u;
  double T = 0.1;



  int oddoreven = 1 ;
  printf("oodoreven = %i \n ", oddoreven);
  printf("oodoreven = %i \n ", int(time));
  printf("oodoreven = %i \n ", int(T));

  oddoreven = int(time / T);

  printf("oodoreven = %i \n ", oddoreven);

  if ((oddoreven / 2) == 0) u = alpha;
  else u = 0;
  u = 1.0;

  b[0] = u ;
  b[1] = 2.0 * u;
  printf("End ComputeU\n ");

}
extern "C"   void computebobserver(double time, unsigned int sizeOfB, double* b, unsigned int sizeOfZ, double* z)
{

  // Si b est de taille 2

  printf("start ComputeU\n ");
  double alpha = 1.0;
  double u;
  double T = 0.1;
  double L[2];


  L[0] = 1.0;
  L[1] = 1.0;


  int oddoreven = 1 ;
  printf("oodoreven = %i \n ", oddoreven);
  printf("oodoreven = %i \n ", int(time));
  printf("oodoreven = %i \n ", int(T));

  oddoreven = int(time / T);

  printf("oodoreven = %i \n ", oddoreven);

  if ((oddoreven / 2) == 0) u = alpha;
  else u = 0;
  u = 1.0;

  b[0] = u + L[0] * z[0];
  b[1] = 2.0 * u + L[1] * z[0];
  printf("End ComputeU\n ");

}

extern "C"   void computeE(double time, unsigned int sizeOfB, double* e, unsigned int sizeOfZ, double* z)
{
  printf("ComputeE\n ");
  // Si b est de taille 2
  double alpha = 1.0;
  double u;
  double T = 0.1;
  int oddoreven = 1 ;
  oddoreven = int(time / T);

  printf("oodoreven = %i \n ", oddoreven);

  if ((oddoreven / 2) == 0) u = alpha;
  else u = 0;

  e[0] = u;
}
