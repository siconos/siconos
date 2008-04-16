#include <stdio.h>
#include <math.h>


extern "C"   double computeControl(double time)
{


  double u;
  double alpha = 10000.0;
  double T = 0.0001;
  int oddoreven = 1 ;

  oddoreven = int(time / T);


  if ((oddoreven / 2) == 0) u = alpha;
  else u = -alpha;
  u = 30 * sin(50 * time);

  return u;
}


extern "C"   void computeU(double time, unsigned int sizeOfB, double* b, unsigned int sizeOfZ, double* z)
{

  double u;
  u = computeControl(time);
  b[0] = u ;
  b[1] = 2.0 * u;
  b[2] = u ;
  b[3] = 2.0 * u;

}


extern "C"   void computeE(double time, unsigned int sizeOfB, double* e, unsigned int sizeOfZ, double* z)
{

  // Si b est de taille 2
  e[0] = computeControl(time) ;
  e[1] = computeControl(time) ;
}
