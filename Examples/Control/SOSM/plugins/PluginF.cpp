#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <stdio.h>
#include <iostream>
#include <math.h>


extern "C" double computeDefault(double time)
{
  double u;
  u = time;
  return u;
}


extern "C" void computef1(double time, unsigned int sizeofx, double* x, double* f, unsigned int sizeofz , double* z)
{ 
  //printf("we are in computef1 \n");
  int sgn1;
  int sgn2;
  if(x[0]>0)
    sgn1 = 1;
  else if (x[0]<0)
    sgn1 = -1;
  else
    sgn1 = 0;

  if(x[1]>0)
    sgn2 = 1;
  else if (x[1]<0)
    sgn2 = -1;
  else
    sgn2 = 0;
  f[0] = x[1];
  f[1] = 0;//pow(fabs(x[0]),1.5) + pow(fabs(x[1]),1.5);	
  //printf("we are out computef1 \n");
}

extern "C" void computeJacf1(double time, unsigned int sizeofx, double* x, double* Jacf, unsigned int sizeofz , double* z)
{ 
  //printf("we are in computeJacf1 \n");
  int sgn1;
  int sgn2;
  
  Jacf[0] = 0.0;
  Jacf[2] = 1.0;
  if(x[0]>0)
    sgn1 = 1;
  else if (x[0]<0)
    sgn1 = -1;
  else
    sgn1 = 0;

  if(x[1]>0)
    sgn2 = 1;
  else if (x[1]<0)
    sgn2 = -1;
  else
    sgn2 = 0;

  Jacf[1] = 0;//0.5*pow(fabs(x[0]),0.5)*sgn1;
  Jacf[3] = 0;//0.5*pow(fabs(x[1]),0.5)*sgn2;
  //printf("we are out computeJacf1 \n");
}


