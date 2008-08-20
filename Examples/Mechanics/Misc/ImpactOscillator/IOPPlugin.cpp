#include <stdio.h>
#include <math.h>

const double omega = 2.6;

extern "C" double FextFunction(const double& time, const double& p)
{
  return cos(p * time);
}


extern "C" void IOFExt(const unsigned int& sizeOfq, const double *time, double *fExt, double *param)
{
  double t = *time;
  //  double p = q[1];
  double p = 3;

  for (unsigned int i = 0; i < sizeOfq; i++)
  {
    fExt[i] = 0.0;
    //    printf("This is loop: %f --- %f \n",sizeOfq,q[1]);
  }

  fExt[0] = FextFunction(t, p);
  fExt[1] = 0;

  //  printf("This is fEXT: %f --- %f --- %f \n", t, fExt[0], p);

}


extern "C" void groundFExt(const unsigned int& sizeOfq, const double *time, double *fExt, double *param)
{
  for (unsigned int i = 0; i < sizeOfq; i++)
    fExt[i] = 0.0;
}

